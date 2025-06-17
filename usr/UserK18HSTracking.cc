// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iostream>
#include <sstream>

#include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "ConfMan.hh"
#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCRawHit.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoCluster.hh"
#include "HodoAnalyzer.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "K18TrackD2U.hh"
#include "HSTrack.hh"
#include "BH1Match.hh"
#include "BH2Filter.hh"
#include "KuramaLib.hh"
#include "MathTools.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UnpackerManager.hh"

#define HodoCut 0 // with BH1/BH2
#define TimeCut 1 // in cluster analysis
#define Chi2Cut 1 // for BcOut tracking
#define SaveBft 1
#define BH1MatchCut 1
#define BH2MatchCut 1

#define KmBeam 0 // BeamA event selection

namespace
{
using namespace root;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
auto& gFilter = BH2Filter::GetInstance();
auto& gBH1Mth = BH1Match::GetInstance();
}

//_____________________________________________________________________________
class UserK18HSTracking : public VEvent
{
private:
  RawData*      rawData;
  HodoAnalyzer* hodoAna;
  DCAnalyzer*   DCAna;

public:
  UserK18HSTracking();
  ~UserK18HSTracking();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserK18HSTracking::ClassName()
{
  static TString s_name("UserK18HSTracking");
  return s_name;
}

//_____________________________________________________________________________
UserK18HSTracking::UserK18HSTracking()
  : VEvent(),
    rawData(new RawData),
    hodoAna(new HodoAnalyzer),
    DCAna(new DCAnalyzer)
{
}

//_____________________________________________________________________________
UserK18HSTracking::~UserK18HSTracking()
{
  if(rawData) delete rawData;
  if(hodoAna) delete hodoAna;
  if(DCAna) delete DCAna;
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;

  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  // Time0
  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;

  // Btof0 BH1
  Double_t Btof0Seg;
  Double_t deBtof0;
  Double_t Btof0;
  Double_t CBtof0;

  //BH1
  Int_t    nhBh1;

  //BH2
  Int_t    nhBh2;

  // BFT
  Int_t    bft_ncl;
  Int_t    bft_ncl_bh1mth;
  Int_t    bft_clsize[NumOfSegBFT];
  Double_t bft_ctime[NumOfSegBFT];
  Double_t bft_clpos[NumOfSegBFT];
  Int_t    bft_bh1mth[NumOfSegBFT];

  // BcOut
  Int_t nlBcOut;
  Int_t ntBcOut;
  Int_t nhBcOut[MaxHits];
  Double_t chisqrBcOut[MaxHits];
  Double_t x0BcOut[MaxHits];
  Double_t y0BcOut[MaxHits];
  Double_t u0BcOut[MaxHits];
  Double_t v0BcOut[MaxHits];

  // K18
  Int_t ntK18;
  Int_t nhK18[MaxHits];
  Double_t p_2nd[MaxHits];
  Double_t p_3rd[MaxHits];
  Double_t delta_2nd[MaxHits];
  Double_t delta_3rd[MaxHits];

  Double_t xin[MaxHits];
  Double_t yin[MaxHits];
  Double_t uin[MaxHits];
  Double_t vin[MaxHits];

  Double_t xout[MaxHits];
  Double_t yout[MaxHits];
  Double_t uout[MaxHits];
  Double_t vout[MaxHits];

  Double_t chisqrK18[MaxHits];
  Double_t xtgtK18[MaxHits];
  Double_t ytgtK18[MaxHits];
  Double_t utgtK18[MaxHits];
  Double_t vtgtK18[MaxHits];

  Double_t theta[MaxHits];
  Double_t phi[MaxHits];
  Double_t layerK18[MaxHits][NumOfLayersBcOut];
  Double_t wireK18[MaxHits][NumOfLayersBcOut];
  Double_t localhitposK18[MaxHits][NumOfLayersBcOut];
  Double_t wposK18[MaxHits][NumOfLayersBcOut];

  //Double_t qHS[MaxHits];
  Double_t initmomHS[MaxHits];

  Double_t xbcHS[MaxHits][NumOfLayersBcOut];
  Double_t ybcHS[MaxHits][NumOfLayersBcOut];
  Double_t zbcHS[MaxHits][NumOfLayersBcOut];
  Double_t ubcHS[MaxHits][NumOfLayersBcOut];
  Double_t vbcHS[MaxHits][NumOfLayersBcOut];

  Double_t xbacHS[MaxHits];
  Double_t ybacHS[MaxHits];
  Double_t zbacHS[MaxHits];
  Double_t ubacHS[MaxHits];
  Double_t vbacHS[MaxHits];
  Double_t pbacHS[MaxHits];

  Double_t xbh2HS[MaxHits];
  Double_t ybh2HS[MaxHits];
  Double_t zbh2HS[MaxHits];
  Double_t ubh2HS[MaxHits];
  Double_t vbh2HS[MaxHits];
  Double_t pbh2HS[MaxHits];

  Double_t xgasvesselHS[MaxHits];
  Double_t ygasvesselHS[MaxHits];
  Double_t zgasvesselHS[MaxHits];
  Double_t ugasvesselHS[MaxHits];
  Double_t vgasvesselHS[MaxHits];
  Double_t pgasvesselHS[MaxHits];

  Double_t xvp1HS[MaxHits];
  Double_t yvp1HS[MaxHits];
  Double_t zvp1HS[MaxHits];
  Double_t uvp1HS[MaxHits];
  Double_t vvp1HS[MaxHits];

  Double_t xvp2HS[MaxHits];
  Double_t yvp2HS[MaxHits];
  Double_t zvp2HS[MaxHits];
  Double_t uvp2HS[MaxHits];
  Double_t vvp2HS[MaxHits];

  Double_t xvp3HS[MaxHits];
  Double_t yvp3HS[MaxHits];
  Double_t zvp3HS[MaxHits];
  Double_t uvp3HS[MaxHits];
  Double_t vvp3HS[MaxHits];

  Double_t xvp4HS[MaxHits];
  Double_t yvp4HS[MaxHits];
  Double_t zvp4HS[MaxHits];
  Double_t uvp4HS[MaxHits];
  Double_t vvp4HS[MaxHits];

  Double_t xhtofHS[MaxHits];
  Double_t yhtofHS[MaxHits];
  Double_t zhtofHS[MaxHits];
  Double_t uhtofHS[MaxHits];
  Double_t vhtofHS[MaxHits];

  Double_t xtgtHS[MaxHits];
  Double_t ytgtHS[MaxHits];
  Double_t ztgtHS[MaxHits];
  Double_t utgtHS[MaxHits];
  Double_t vtgtHS[MaxHits];
  Double_t pHS[MaxHits];
  Double_t thetaHS[MaxHits];
  Double_t phiHS[MaxHits];
  Double_t pathHS[MaxHits];
  Double_t m2[MaxHits];

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum          =  0;
  bft_ncl        =  0;
  bft_ncl_bh1mth =  0;
  nlBcOut        =  0;
  ntBcOut        =  0;
  ntK18          =  0;
  Time0Seg       = -1;
  deTime0        = qnan;
  Time0          = qnan;
  CTime0         = qnan;
  nhBh1          = 0;
  nhBh2          = 0;
  Btof0Seg       = -1;
  deBtof0        = qnan;
  Btof0          = qnan;
  CBtof0         = qnan;

  for(Int_t it=0; it<NumOfSegTrig; it++){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<NumOfSegBFT; it++){
    bft_clsize[it] = qnan;
    bft_ctime[it]  = qnan;
    bft_clpos[it]  = qnan;
    bft_bh1mth[it] = -1;
  }

  for(Int_t i = 0; i<MaxHits; ++i){
    nhBcOut[i] = 0;
    chisqrBcOut[i] = qnan;
    x0BcOut[i] = qnan;
    y0BcOut[i] = qnan;
    u0BcOut[i] = qnan;
    v0BcOut[i] = qnan;
    p_2nd[i] = qnan;
    p_3rd[i] = qnan;
    delta_2nd[i] = qnan;
    delta_3rd[i] = qnan;
    xin[i] = qnan;
    yin[i] = qnan;
    uin[i] = qnan;
    vin[i] = qnan;
    xout[i] = qnan;
    yout[i] = qnan;
    uout[i] = qnan;
    vout[i] = qnan;
    nhK18[i] = 0;
    chisqrK18[i] = qnan;
    xtgtK18[i] = qnan;
    ytgtK18[i] = qnan;
    utgtK18[i] = qnan;
    vtgtK18[i] = qnan;
    theta[i] = qnan;
    phi[i] = qnan;

    //qHS[i] = qnan;
    initmomHS[i] = qnan;

    xbacHS[i] = qnan;
    ybacHS[i] = qnan;
    zbacHS[i] = qnan;
    ubacHS[i] = qnan;
    vbacHS[i] = qnan;
    pbacHS[i] = qnan;

    xbh2HS[i] = qnan;
    ybh2HS[i] = qnan;
    zbh2HS[i] = qnan;
    ubh2HS[i] = qnan;
    vbh2HS[i] = qnan;
    pbh2HS[i] = qnan;

    xgasvesselHS[i] = qnan;
    ygasvesselHS[i] = qnan;
    zgasvesselHS[i] = qnan;
    ugasvesselHS[i] = qnan;
    vgasvesselHS[i] = qnan;
    pgasvesselHS[i] = qnan;

    xvp1HS[i] = qnan;
    yvp1HS[i] = qnan;
    zvp1HS[i] = qnan;
    uvp1HS[i] = qnan;
    vvp1HS[i] = qnan;

    xvp2HS[i] = qnan;
    yvp2HS[i] = qnan;
    zvp2HS[i] = qnan;
    ubh2HS[i] = qnan;
    vbh2HS[i] = qnan;

    xvp3HS[i] = qnan;
    yvp3HS[i] = qnan;
    zvp3HS[i] = qnan;
    uvp3HS[i] = qnan;
    vvp3HS[i] = qnan;

    xvp4HS[i] = qnan;
    yvp4HS[i] = qnan;
    zvp4HS[i] = qnan;
    uvp4HS[i] = qnan;
    vvp4HS[i] = qnan;

    xhtofHS[i] = qnan;
    yhtofHS[i] = qnan;
    zhtofHS[i] = qnan;
    uhtofHS[i] = qnan;
    vhtofHS[i] = qnan;

    xtgtHS[i] = qnan;
    ytgtHS[i] = qnan;
    ztgtHS[i] = qnan;
    utgtHS[i] = qnan;
    vtgtHS[i] = qnan;
    pHS[i] = qnan;
    thetaHS[i] = qnan;
    phiHS[i] = qnan;
    pathHS[i] = qnan;
    m2[i] = qnan;

    for(Int_t j=0; j<NumOfLayersBcOut; j++){
      layerK18[i][j] = qnan;
      wireK18[i][j] = qnan;
      localhitposK18[i][j] = qnan;
      wposK18[i][j] = qnan;

      xbcHS[i][j] = qnan;
      ybcHS[i][j] = qnan;
      zbcHS[i][j] = qnan;
      ubcHS[i][j] = qnan;
      vbcHS[i][j] = qnan;
    }
  }
}

//_____________________________________________________________________________
namespace root
{
Event event;
TH1   *h[MaxHist];
TTree *tree;
const Int_t BFTHid = 10000;
enum eParticle
{
  kKaon, kPion, nParticle
};
}

//_____________________________________________________________________________
Bool_t
UserK18HSTracking::ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserK18HSTracking::ProcessingNormal()
{
#if HodoCut
  static const auto MinDeBH2   = gUser.GetParameter("DeBH2", 0);
  static const auto MaxDeBH2   = gUser.GetParameter("DeBH2", 1);
  static const auto MinDeBH1   = gUser.GetParameter("DeBH1", 0);
  static const auto MaxDeBH1   = gUser.GetParameter("DeBH1", 1);
  static const auto MinBeamToF = gUser.GetParameter("BTOF",  0);
  static const auto MaxBeamToF = gUser.GetParameter("BTOF",  1);
#endif
#if TimeCut
  static const auto MinTimeBFT = gUser.GetParameter("TimeBFT", 0);
  static const auto MaxTimeBFT = gUser.GetParameter("TimeBFT", 1);
#endif
  static const auto MaxChisqrBcOut = gUser.GetParameter("MaxChisqrBcOut");
  static const auto MinDriftTimeBC34 = gUser.GetParameter("DriftTimeBC34", 0);
  static const auto MaxDriftTimeBC34 = gUser.GetParameter("DriftTimeBC34", 1);
  static const auto MinTotBcOut = gUser.GetParameter("MinTotBcOut");
  static const auto MaxMultiHitBcOut = gUser.GetParameter("MaxMultiHitBcOut");

  const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);

  rawData->DecodeHits();

  event.evnum = gUnpacker.get_event_number();

  HF1(1, 0.);

  ///// Trigger Flag
  std::bitset<NumOfSegTrig> trigger_flag;
  {
    for(const auto& hit: rawData->GetTrigRawHC()){
      Int_t seg = hit->SegmentId();
      Int_t tdc = hit->GetTdc1();
      if(tdc > 0){
	event.trigpat[trigger_flag.count()] = seg;
	event.trigflag[seg]  = tdc;
        trigger_flag.set(seg);
	HF1(10, seg);
	HF1(10+seg, tdc);
      }
    }
  }

  if(trigger_flag[trigger::kSpillEnd]) return true;

  HF1(1, 1);

#if KmBeam
  if(event.trigflag[14] < 0) return false; //Select BeamA events only
#endif

  ////////// BH2 time 0
  hodoAna->DecodeBH2Hits(rawData);
#if HodoCut
  Int_t nhBh2 = hodoAna->GetNHitsBH2();
  if(nhBh2==0) return true;
#endif
  HF1(1, 2);

  //////////////BH2 Analysis
  Double_t t0_seg = -1;
  BH2Cluster *cl_time0 = hodoAna->GetTime0BH2Cluster();
  if(cl_time0){
    t0_seg = cl_time0->MeanSeg();
    event.Time0Seg = cl_time0->MeanSeg()+1;
    event.deTime0  = cl_time0->DeltaE();
    event.Time0    = cl_time0->Time0();
    event.CTime0   = cl_time0->CTime0();
  }else{
    return true;
  }
  event.nhBh2 = hodoAna->GetNClustersBH2();

  HF1(1, 3);

  ////////// BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData);
#if HodoCut
  Int_t nhBh1 = hodoAna->GetNHitsBH1();
  if(nhBh1==0) return true;
#endif
  HF1(1, 4);

  Double_t btof0_seg = -1;
  HodoCluster* cl_btof0 = event.Time0Seg > 0? hodoAna->GetBtof0BH1Cluster(event.CTime0) : NULL;
  if(cl_btof0){
    btof0_seg = cl_btof0->MeanSeg();
    event.Btof0Seg = cl_btof0->MeanSeg()+1;
    event.deBtof0 = cl_btof0->DeltaE();
    event.Btof0 = cl_btof0->MeanTime() - event.Time0;
    event.CBtof0 = cl_btof0->CMeanTime() - event.CTime0;
  }
  event.nhBh1 = hodoAna->GetNClustersBH1();

  HF1(1, 5);

  HF1(1, 6);

  std::vector<Double_t> xCand;
  ////////// BFT
  {
    hodoAna->DecodeBFTHits(rawData);
    // Fiber Cluster
    Int_t ncl_raw = hodoAna->GetNClustersBFT();
    for(Int_t i=0; i<ncl_raw; ++i){
      FiberCluster *cl = hodoAna->GetClusterBFT(i);
      if(!cl) continue;
      Double_t ctime  = cl->CMeanTime();
      HF1(BFTHid +103, ctime);
    }
#if TimeCut
    hodoAna->TimeCutBFT(MinTimeBFT, MaxTimeBFT);
#endif
    Int_t ncl = hodoAna->GetNClustersBFT();
    event.bft_ncl = ncl;
    HF1(BFTHid +100, ncl_raw);
    HF1(BFTHid +101, ncl);
    for(Int_t i=0; i<ncl; ++i){
      FiberCluster *cl = hodoAna->GetClusterBFT(i);
      if(!cl) continue;
      Double_t clsize = cl->ClusterSize();
      Double_t ctime  = cl->CMeanTime();
      Double_t pos    = cl->MeanPosition();
      // Double_t width  = cl->Width();

      event.bft_clsize[i] = clsize;
      event.bft_ctime[i]  = ctime;
      event.bft_clpos[i]  = pos;
#if BH1MatchCut
      //Veto events which BH1-BFT are not matched
      if(btof0_seg > 0 && ncl != 0){
	if(gBH1Mth.Judge(pos, btof0_seg)){
	  event.bft_bh1mth[i] = 1;
	  xCand.push_back(pos);
	}
      }
#else
      if(btof0_seg > 0 && ncl != 1){
	if(gBH1Mth.Judge(pos, btof0_seg)){
	  event.bft_bh1mth[i] = 1;
	  xCand.push_back(pos);
	}
      }else{
	xCand.push_back(pos);
      }
#endif
      HF1(BFTHid +102, clsize);
      HF1(BFTHid +104, ctime);
      HF1(BFTHid +105, pos);
    }

    event.bft_ncl_bh1mth = xCand.size();
    HF1(BFTHid + 106, event.bft_ncl_bh1mth);
  }

  HF1(1, 7.);
  if(xCand.size()!=1) return true;
  HF1(1, 10.);

  DCAna->DecodeRawHits(rawData);
  DCAna->TotCutBCOut(MinTotBcOut);
  DCAna->DriftTimeCutBC34(MinDriftTimeBC34, MaxDriftTimeBC34);

  ////////////// BC3&4 number of hit in one layer not 0
  Double_t multi_BcOut=0.;
  {
    Int_t nlBcOut = 0;
    for(Int_t layer=1; layer<=NumOfLayersBcOut; ++layer){
      const DCHitContainer &contBcOut = DCAna->GetBcOutHC(layer);
      Int_t nhBcOut = contBcOut.size();
      multi_BcOut += Double_t(nhBcOut);
      if(nhBcOut>0) nlBcOut++;
    }
    event.nlBcOut = nlBcOut;
  }

  if(multi_BcOut/Double_t(NumOfLayersBcOut) > MaxMultiHitBcOut) return true;

  HF1(1, 11.);

  //////////////BCOut tracking
  // BH2Filter::FilterList cands;
  // gFilter.Apply((Int_t)event.Time0Seg-1, *DCAna, cands);
  //DCAna->TrackSearchBcOut(cands, event.Time0Seg-1);
  //  DCAna->TrackSearchBcOut(-1);
#if BH2MatchCut
  if(t0_seg<0) return true;
  DCAna->TrackSearchBcOut(t0_seg);
#else
  DCAna->TrackSearchBcOut();
#endif

 #if Chi2Cut
  DCAna->ChiSqrCutBcOut(MaxChisqrBcOut);
 #endif

  Int_t ntBcOut = DCAna->GetNtracksBcOut();
  event.ntBcOut = ntBcOut;
  if(ntBcOut > MaxHits){
    std::cout << "#W too many BcOut tracks : ntBcOut = "
	      << ntBcOut << std::endl;
    ntBcOut = MaxHits;
  }
  HF1(50, Double_t(ntBcOut));
  for(Int_t it=0; it<ntBcOut; ++it){
    DCLocalTrack *tp = DCAna->GetTrackBcOut(it);
    Int_t nh = tp->GetNHit();
    Double_t chisqr = tp->GetChiSquare();
    Double_t u0 = tp->GetU0(),  v0 = tp->GetV0();
    Double_t x0 = tp->GetX(0.), y0 = tp->GetY(0.);

    HF1(51, Double_t(nh));
    HF1(52, chisqr);
    HF1(54, x0); HF1(55, y0);
    HF1(56, u0); HF1(57, v0);
    HF2(58, x0, u0); HF2(39, y0, v0);
    HF2(60, x0, y0);

    event.nhBcOut[it] = nh;
    event.chisqrBcOut[it] = chisqr;
    event.x0BcOut[it] = x0;
    event.y0BcOut[it] = y0;
    event.u0BcOut[it] = u0;
    event.v0BcOut[it] = v0;

    for(Int_t ih=0; ih<nh; ++ih){
      DCLTrackHit *hit=tp->GetHit(ih);
      Int_t layerId=hit->GetLayer()-100;
      HF1(53, layerId);
    }
  }
  if(ntBcOut==0) return true;

  HF1(1, 12.);

  HF1(1, 20.);
  ////////// K18HSTracking D2U
  DCAna->TrackSearchK18D2U(xCand);
  Int_t ntK18 = DCAna->GetNTracksK18D2U();
  if(ntK18 > MaxHits){
    std::cout << "#W too many ntK18 "
	      << ntK18 << "/" << MaxHits << std::endl;
    ntK18 = MaxHits;
  }
  event.ntK18 = ntK18;
  HF1(70, Double_t(ntK18));
  for(Int_t i=0; i<ntK18; ++i){
    K18TrackD2U *tp=DCAna->GetK18TrackD2U(i);
    if(!tp) continue;
    DCLocalTrack *track = tp->TrackOut();
    std::size_t nh = track->GetNHit();
    Double_t chisqr = track->GetChiSquare();

    Double_t xin=tp->Xin(), yin=tp->Yin();
    Double_t uin=tp->Uin(), vin=tp->Vin();

    Double_t xt=tp->Xtgt(), yt=tp->Ytgt();
    Double_t ut=tp->Utgt(), vt=tp->Vtgt();

    Double_t xout=tp->Xout(), yout=tp->Yout();
    Double_t uout=tp->Uout(), vout=tp->Vout();

    Double_t p_2nd=tp->P();
    Double_t p_3rd=tp->P3rd();
    Double_t delta_2nd=tp->Delta();
    Double_t delta_3rd=tp->Delta3rd();
    Double_t theta = track->GetTheta();
    Double_t phi   = track->GetPhi();

    HF1(81, p_3rd); HF1(82, delta_3rd);

    event.p_2nd[i] = p_2nd;
    event.p_3rd[i] = p_3rd;
    event.delta_2nd[i] = delta_2nd;
    event.delta_3rd[i] = delta_3rd;

    event.xin[i] = xin;
    event.yin[i] = yin;
    event.uin[i] = uin;
    event.vin[i] = vin;

    event.xout[i] = xout;
    event.yout[i] = yout;
    event.uout[i] = uout;
    event.vout[i] = vout;

    event.nhK18[i]     = nh;
    event.chisqrK18[i] = chisqr;
    event.xtgtK18[i]   = xt;
    event.ytgtK18[i]   = yt;
    event.utgtK18[i]   = ut;
    event.vtgtK18[i]   = vt;
    event.theta[i] = theta;
    event.phi[i]   = phi;
    for(Int_t j=0; j<nh; ++j){
      DCLTrackHit *hit=track->GetHit(j);
      event.layerK18[i][j] = hit->GetLayer();
      event.wireK18[i][j] = hit->GetWire();
      event.localhitposK18[i][j] = hit->GetLocalHitPos();
      event.wposK18[i][j] = hit->GetWirePosition();
    }

    ////////// K18HSTracking Propagation to Tgt
    static const auto StofOffset = abs(gUser.GetParameter("StofOffset"));

    const Double_t& pK18 = ConfMan::Get<Double_t>("PK18");
    auto trHS = new HSTrack(xout, yout, uout, vout, p_3rd);
    if(!trHS) continue;
    if(BeamThroughTPC){
      Int_t pikp = gUser.GetParameter("BeamThroughPID");
      trHS -> SetPID(pikp);
    }
    trHS->Propagate();
    const auto& PosTgt = trHS->TgtPosition();
    const auto& MomTgt = trHS->TgtMomentum();
    Double_t xtgt = PosTgt.x(), ytgt = PosTgt.y(), ztgt = PosTgt.z();
    Double_t pHS = MomTgt.Mag();
    //Double_t qHS = trHS->Polarity();

    Double_t utgt = MomTgt.x()/MomTgt.z(), vtgt = MomTgt.y()/MomTgt.z();
    Double_t costHS = 1./TMath::Sqrt(1.+utgt*utgt+vtgt*vtgt);
    Double_t thetaHS = TMath::ACos(costHS)*TMath::RadToDeg();
    Double_t phiHS = TMath::ATan2(utgt, vtgt);
    Double_t initial_momentum = trHS->GetInitialMomentum();

    HF1(74, xtgt); HF1(75, ytgt); HF1(76, utgt); HF1(77, vtgt);
    HF2(78, xtgt, utgt); HF2(79, ytgt, vtgt); HF2(80, xtgt, ytgt);

    const auto& PosBAC = trHS->BACPosition();
    const auto& MomBAC = trHS->BACMomentum();
    Double_t xBAC = PosBAC.x(), yBAC = PosBAC.y(), zBAC = PosBAC.z();
    Double_t uBAC = MomBAC.x()/MomBAC.z(), vBAC = MomBAC.y()/MomBAC.z();
    Double_t pBACHS = MomBAC.Mag();

    const auto& PosBH2 = trHS->BH2Position();
    const auto& MomBH2 = trHS->BH2Momentum();
    Double_t xBH2 = PosBH2.x(), yBH2 = PosBH2.y(), zBH2 = PosBH2.z();
    Double_t uBH2 = MomBH2.x()/MomBH2.z(), vBH2 = MomBH2.y()/MomBH2.z();
    Double_t pBH2HS = MomBH2.Mag();

    const auto& PosVP1 = trHS->VP1Position();
    const auto& MomVP1 = trHS->VP1Momentum();
    Double_t xVP1 = PosVP1.x(), yVP1 = PosVP1.y(), zVP1 = PosVP1.z();
    Double_t uVP1 = MomVP1.x()/MomVP1.z(), vVP1 = MomVP1.y()/MomVP1.z();

    const auto& PosVP2 = trHS->VP2Position();
    const auto& MomVP2 = trHS->VP2Momentum();
    Double_t xVP2 = PosVP2.x(), yVP2 = PosVP2.y(), zVP2 = PosVP2.z();
    Double_t uVP2 = MomVP2.x()/MomVP1.z(), vVP2 = MomVP2.y()/MomVP2.z();

    const auto& PosVP3 = trHS->VP3Position();
    const auto& MomVP3 = trHS->VP3Momentum();
    Double_t xVP3 = PosVP3.x(), yVP3 = PosVP3.y(), zVP3 = PosVP3.z();
    Double_t uVP3 = MomVP3.x()/MomVP3.z(), vVP3 = MomVP3.y()/MomVP3.z();

    const auto& PosVP4 = trHS->VP4Position();
    const auto& MomVP4 = trHS->VP4Momentum();
    Double_t xVP4 = PosVP4.x(), yVP4 = PosVP4.y(), zVP4 = PosVP4.z();
    Double_t uVP4 = MomVP4.x()/MomVP4.z(), vVP4 = MomVP4.y()/MomVP4.z();

    const auto& PosHtof = trHS->HtofPosition();
    const auto& MomHtof = trHS->HtofMomentum();
    Double_t xHtof = PosHtof.x(), yHtof = PosHtof.y(), zHtof = PosHtof.z();
    Double_t uHtof = MomHtof.x()/MomHtof.z(), vHtof = MomHtof.y()/MomHtof.z();

    const auto& PosGasVesselU = trHS->GasVesselUPosition();
    const auto& MomGasVesselU = trHS->GasVesselUMomentum();
    Double_t xGasVesselU = PosGasVesselU.x(), yGasVesselU = PosGasVesselU.y(), zGasVesselU = PosGasVesselU.z();
    Double_t uGasVesselU = MomGasVesselU.x()/MomGasVesselU.z(), vGasVesselU = MomGasVesselU.y()/MomGasVesselU.z();
    Double_t pGasVesselUHS = MomGasVesselU.Mag();

    const auto& PosGasVesselD = trHS->GasVesselDPosition();
    const auto& MomGasVesselD = trHS->GasVesselDMomentum();
    Double_t xGasVesselD = PosGasVesselD.x(), yGasVesselD = PosGasVesselD.y(), zGasVesselD = PosGasVesselD.z();
    Double_t uGasVesselD = MomGasVesselD.x()/MomGasVesselD.z(), vGasVesselD = MomGasVesselD.y()/MomGasVesselD.z();

    HF2(90, xBH2, yBH2);
    HF2(91, xGasVesselU, yGasVesselU);
    HF2(92, xGasVesselD, yGasVesselD);

    for(Int_t j=0; j<NumOfLayersBcOut; ++j){
      const auto& PosBC = trHS->BCPosition(j);
      const auto& MomBC = trHS->BCMomentum(j);
      Double_t xBC = PosBC.x(), yBC = PosBC.y(), zBC = PosBC.z();
      Double_t uBC = MomBC.x()/MomBC.z(), vBC = MomBC.y()/MomBC.z();

      event.xbcHS[i][j] = xBC;
      event.ybcHS[i][j] = yBC;
      event.zbcHS[i][j] = zBC;
      event.ubcHS[i][j] = uBC;
      event.vbcHS[i][j] = vBC;
    }

    //BH2 - Tgt
    Double_t path = trHS->PathLength();
    Double_t m2 = Kinematics::MassSquare(pHS, path, StofOffset);

    event.xbacHS[i] = xBAC;
    event.ybacHS[i] = yBAC;
    event.zbacHS[i] = zBAC;
    event.ubacHS[i] = uBAC;
    event.vbacHS[i] = vBAC;
    event.pbacHS[i] = pBACHS;

    event.xbh2HS[i] = xBH2;
    event.ybh2HS[i] = yBH2;
    event.zbh2HS[i] = zBH2;
    event.ubh2HS[i] = uBH2;
    event.vbh2HS[i] = vBH2;
    event.pbh2HS[i] = pBH2HS;

    event.xgasvesselHS[i] = xGasVesselU;
    event.ygasvesselHS[i] = yGasVesselU;
    event.zgasvesselHS[i] = zGasVesselU;
    event.ugasvesselHS[i] = uGasVesselU;
    event.vgasvesselHS[i] = vGasVesselU;
    event.pgasvesselHS[i] = pGasVesselUHS;

    event.xvp1HS[i] = xVP1;
    event.yvp1HS[i] = yVP1;
    event.zvp1HS[i] = zVP1;
    event.uvp1HS[i] = uVP1;
    event.vvp1HS[i] = vVP1;

    event.xvp2HS[i] = xVP2;
    event.yvp2HS[i] = yVP2;
    event.zvp2HS[i] = zVP2;
    event.uvp2HS[i] = uVP2;
    event.vvp2HS[i] = vVP2;

    event.xvp3HS[i] = xVP3;
    event.yvp3HS[i] = yVP3;
    event.zvp3HS[i] = zVP3;
    event.uvp3HS[i] = uVP3;
    event.vvp3HS[i] = vVP3;

    event.xvp4HS[i] = xVP4;
    event.yvp4HS[i] = yVP4;
    event.zvp4HS[i] = zVP4;
    event.uvp4HS[i] = uVP4;
    event.vvp4HS[i] = vVP4;

    event.xhtofHS[i] = xHtof;
    event.yhtofHS[i] = yHtof;
    event.zhtofHS[i] = zHtof;
    event.uhtofHS[i] = uHtof;
    event.vhtofHS[i] = vHtof;

    event.xtgtHS[i] = xtgt;
    event.ytgtHS[i] = ytgt;
    event.ztgtHS[i] = ztgt;
    event.utgtHS[i] = utgt;
    event.vtgtHS[i] = vtgt;
    event.pHS[i] = std::abs(pHS);
    event.thetaHS[i] = thetaHS;
    event.phiHS[i] = phiHS;
    event.pathHS[i] = path;
    event.m2[i] = m2;
    //event.qHS[i] = qHS;
    event.initmomHS[i] = initial_momentum;
  }
  HF1(1, 22.);

  return true;

}

//_____________________________________________________________________________
Bool_t
UserK18HSTracking::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserK18HSTracking;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
#if KmBeam == 0
  const Int_t    NbinTot =  136;
  const Double_t MinTot  =   -8.;
  const Double_t MaxTot  =  128.;
  const Int_t    NbinTime = 200;
  const Double_t MinTime  = -10.;
  const Double_t MaxTime  =  10.;

  HB1(1, "Status",  30,   0., 30.);
  HB1(10, "Trigger HitPat", NumOfSegTrig, 0, NumOfSegTrig);
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1(10+i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000);
  }

#if SaveBft
  //BFT
  HB1(BFTHid +21, "BFT CTime U",     NbinTime, MinTime, MaxTime);
  HB1(BFTHid +22, "BFT CTime D",     NbinTime, MinTime, MaxTime);
  HB2(BFTHid +23, "BFT CTime/Tot U", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB2(BFTHid +24, "BFT CTime/Tot D", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);

  HB1(BFTHid +100, "BFT NCluster [Raw]", 100, 0, 100);
  HB1(BFTHid +101, "BFT NCluster [TimeCut]", 10, 0, 10);
  HB1(BFTHid +102, "BFT Cluster Size", 5, 0, 5);
  HB1(BFTHid +103, "BFT CTime (Cluster) [Raw]", NbinTime, MinTime, MaxTime);
  HB1(BFTHid +104, "BFT CTime (Cluster) [TimeCut]", NbinTime, MinTime, MaxTime);
  HB1(BFTHid +105, "BFT Cluster Position",
      NumOfSegBFT, -0.5*(Double_t)NumOfSegBFT, 0.5*(Double_t)NumOfSegBFT);
  HB1(BFTHid +106, "BFT NCluster [TimeCut && BH1Matching]", 10, 0, 10);
#endif

  // BcOut
  HB1(50, "#Tracks BcOut", 10, 0., 10.);
  HB1(51, "#Hits of Track BcOut", 20, 0., 20.);
  HB1(52, "Chisqr BcOut", 500, 0., 50.);
  HB1(53, "LayerId BcOut", 15, 12., 27.);
  HB1(54, "X0 BcOut", 400, -100., 100.);
  HB1(55, "Y0 BcOut", 400, -100., 100.);
  HB1(56, "U0 BcOut",  200, -0.20, 0.20);
  HB1(57, "V0 BcOut",  200, -0.20, 0.20);
  HB2(58, "U0%X0 BcOut", 100, -100., 100., 100, -0.20, 0.20);
  HB2(59, "V0%Y0 BcOut", 100, -100., 100., 100, -0.20, 0.20);
  HB2(60, "X0%Y0 BcOut", 100, -100., 100., 100, -100., 100.);

  // K18
  HB1(70, "#Tracks K18", 20, 0., 20.);
  HB1(71, "#Hits of K18Track", 30, 0., 30.);
  HB1(72, "Chisqr K18Track", 500, 0., 100.);
  HB1(73, "LayerId K18Track", 50, 0., 50.);
  HB1(74, "Xtgt K18Track", 200, -100., 100.);
  HB1(75, "Ytgt K18Track", 200, -100., 100.);
  HB1(76, "Utgt K18Track", 300, -0.30, 0.30);
  HB1(77, "Vtgt K18Track", 300, -0.20, 0.20);
  HB2(78, "U%Xtgt K18Track", 100, -100., 100., 100, -0.25, 0.25);
  HB2(79, "V%Ytgt K18Track", 100, -100., 100., 100, -0.10, 0.10);
  HB2(80, "Y%Xtgt K18Track", 100, -100., 100., 100, -100., 100.);
  HB1(81, "P K18Track", 500, 0.50, 2.0);
  HB1(82, "dP K18Track", 200, -0.1, 0.1);

  HB2(90, "Y%X@BH2", 150, -150., 150., 150, -150., 150.);
  HB2(91, "Y%X@GasVessel-U", 500, -500., 500., 400, -400., 400.);
  HB2(92, "Y%X@GasVessel-D", 500, -500., 500., 400, -400., 400.);
#endif
  //tree
  HBTree("k18track","Data Summary Table of K18HSTracking");
  // Trigger Flag
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  tree->Branch("Time0Seg", &event.Time0Seg,  "Time0Seg/D");
  tree->Branch("deTime0",  &event.deTime0,   "deTime0/D");
  tree->Branch("Time0",    &event.Time0,     "Time0/D");
  tree->Branch("CTime0",   &event.CTime0,    "CTime0/D");

  tree->Branch("Btof0Seg", &event.Btof0Seg,  "Btof0Seg/D");
  tree->Branch("deBtof0",  &event.deBtof0,   "deBtof0/D");
  tree->Branch("Btof0",    &event.Btof0,     "Btof0/D");
  tree->Branch("CBtof0",   &event.CBtof0,    "CBtof0/D");

  tree->Branch("nhBh1",    &event.nhBh1,     "nhBh1/I");
  tree->Branch("nhBh2",    &event.nhBh2,     "nhBh2/I");

#if SaveBft
  //BFT
  tree->Branch("bft_ncl",        &event.bft_ncl,    "bft_ncl/I");
  tree->Branch("bft_ncl_bh1mth", &event.bft_ncl_bh1mth, "bft_ncl_bh1mth/I");
  tree->Branch("bft_clsize",      event.bft_clsize, "bft_clsize[bft_ncl]/I");
  tree->Branch("bft_ctime",       event.bft_ctime,  "bft_ctime[bft_ncl]/D");
  tree->Branch("bft_clpos",       event.bft_clpos,  "bft_clpos[bft_ncl]/D");
  tree->Branch("bft_bh1mth",      event.bft_bh1mth, "bft_bh1mth[bft_ncl]/I");
#endif

  // BcOut
  tree->Branch("nlBcOut",   &event.nlBcOut,     "nlBcOut/I");
  tree->Branch("ntBcOut",   &event.ntBcOut,     "ntBcOut/I");
  tree->Branch("nhBcOut",    event.nhBcOut,     "nhBcOut[ntBcOut]/I");
  tree->Branch("chisqrBcOut",event.chisqrBcOut, "chisqrIn[ntBcOut]/D");
  tree->Branch("x0BcOut",    event.x0BcOut,     "x0BcOut[ntBcOut]/D");
  tree->Branch("y0BcOut",    event.y0BcOut,     "y0BcOut[ntBcOut]/D");
  tree->Branch("u0BcOut",    event.u0BcOut,     "u0BcOut[ntBcOut]/D");
  tree->Branch("v0BcOut",    event.v0BcOut,     "v0BcOut[ntBcOut]/D");

  // K18
  tree->Branch("ntK18",      &event.ntK18,     "ntK18/I");
  tree->Branch("nhK18",       event.nhK18,     "nhK18[ntK18]/I");
  tree->Branch("chisqrK18",   event.chisqrK18, "chisqrK18[ntK18]/D");
  tree->Branch("p_2nd",       event.p_2nd,     "p_2nd[ntK18]/D");
  tree->Branch("p_3rd",       event.p_3rd,     "p_3rd[ntK18]/D");
  tree->Branch("delta_2nd",   event.delta_2nd, "delta_2nd[ntK18]/D");
  tree->Branch("delta_3rd",   event.delta_3rd, "delta_3rd[ntK18]/D");
  tree->Branch("xin",    event.xin,   "xin[ntK18]/D");
  tree->Branch("yin",    event.yin,   "yin[ntK18]/D");
  tree->Branch("uin",    event.uin,   "uin[ntK18]/D");
  tree->Branch("vin",    event.vin,   "vin[ntK18]/D");
  tree->Branch("xout",    event.xout,   "xout[ntK18]/D");
  tree->Branch("yout",    event.yout,   "yout[ntK18]/D");
  tree->Branch("uout",    event.uout,   "uout[ntK18]/D");
  tree->Branch("vout",    event.vout,   "vout[ntK18]/D");
  tree->Branch("xtgtK18",    event.xtgtK18,   "xtgtK18[ntK18]/D");
  tree->Branch("ytgtK18",    event.ytgtK18,   "ytgtK18[ntK18]/D");
  tree->Branch("utgtK18",    event.utgtK18,   "utgtK18[ntK18]/D");
  tree->Branch("vtgtK18",    event.vtgtK18,   "vtgtK18[ntK18]/D");
  tree->Branch("theta",   event.theta,  "theta[ntK18]/D");
  tree->Branch("phi",     event.phi,    "phi[ntK18]/D");

  tree->Branch("layerK18", event.layerK18, "layerK18[ntK18][12]/D");
  tree->Branch("wireK18", event.wireK18, "wireK18[ntK18][12]/D");
  tree->Branch("localhitposK18", event.localhitposK18, "localhitposK18[ntK18][12]/D");
  tree->Branch("wposK18", event.wposK18, "wposK18[ntK18][12]/D");

  // HS Propagation
  //tree->Branch("qHS",   event.qHS, "qHS[ntK18]/D");
  tree->Branch("xbcHS",   event.xbcHS, "xbcHS[ntK18][12]/D");
  tree->Branch("ybcHS",   event.ybcHS, "ybcHS[ntK18][12]/D");
  tree->Branch("zbcHS",   event.zbcHS, "zbcHS[ntK18][12]/D");
  tree->Branch("ubcHS",   event.ubcHS, "ubcHS[ntK18][12]/D");
  tree->Branch("vbcHS",   event.vbcHS, "vbcHS[ntK18][12]/D");

  tree->Branch("xbacHS",   event.xbacHS, "xbacHS[ntK18]/D");
  tree->Branch("ybacHS",   event.ybacHS, "ybacHS[ntK18]/D");
  tree->Branch("zbacHS",   event.zbacHS, "zbacHS[ntK18]/D");
  tree->Branch("ubacHS",   event.ubacHS, "ubacHS[ntK18]/D");
  tree->Branch("vbacHS",   event.vbacHS, "vbacHS[ntK18]/D");
  tree->Branch("pbacHS",   event.pbacHS, "pbacHS[ntK18]/D");

  tree->Branch("xbh2HS",   event.xbh2HS, "xbh2HS[ntK18]/D");
  tree->Branch("ybh2HS",   event.ybh2HS, "ybh2HS[ntK18]/D");
  tree->Branch("zbh2HS",   event.zbh2HS, "zbh2HS[ntK18]/D");
  tree->Branch("ubh2HS",   event.ubh2HS, "ubh2HS[ntK18]/D");
  tree->Branch("vbh2HS",   event.vbh2HS, "vbh2HS[ntK18]/D");
  tree->Branch("pbh2HS",   event.pbh2HS, "pbh2HS[ntK18]/D");

  tree->Branch("xgasvesselHS",   event.xgasvesselHS, "xgasvesselHS[ntK18]/D");
  tree->Branch("ygasvesselHS",   event.ygasvesselHS, "ygasvesselHS[ntK18]/D");
  tree->Branch("zgasvesselHS",   event.zgasvesselHS, "zgasvesselHS[ntK18]/D");
  tree->Branch("ugasvesselHS",   event.ugasvesselHS, "ugasvesselHS[ntK18]/D");
  tree->Branch("vgasvesselHS",   event.vgasvesselHS, "vgasvesselHS[ntK18]/D");
  tree->Branch("pgasvesselHS",   event.pgasvesselHS, "pgasvesselHS[ntK18]/D");

  tree->Branch("xvp1HS",   event.xvp1HS, "xvp1HS[ntK18]/D");
  tree->Branch("yvp1HS",   event.yvp1HS, "yvp1HS[ntK18]/D");
  tree->Branch("zvp1HS",   event.zvp1HS, "zvp1HS[ntK18]/D");
  tree->Branch("uvp1HS",   event.uvp1HS, "uvp1HS[ntK18]/D");
  tree->Branch("vvp1HS",   event.vvp1HS, "vvp1HS[ntK18]/D");

  tree->Branch("xvp2HS",   event.xvp2HS, "xvp2HS[ntK18]/D");
  tree->Branch("yvp2HS",   event.yvp2HS, "yvp2HS[ntK18]/D");
  tree->Branch("zvp2HS",   event.zvp2HS, "zvp2HS[ntK18]/D");
  tree->Branch("uvp2HS",   event.uvp2HS, "uvp2HS[ntK18]/D");
  tree->Branch("vvp2HS",   event.vvp2HS, "vvp2HS[ntK18]/D");

  tree->Branch("xvp3HS",   event.xvp3HS, "xvp3HS[ntK18]/D");
  tree->Branch("yvp3HS",   event.yvp3HS, "yvp3HS[ntK18]/D");
  tree->Branch("zvp3HS",   event.zvp3HS, "zvp3HS[ntK18]/D");
  tree->Branch("uvp3HS",   event.uvp3HS, "uvp3HS[ntK18]/D");
  tree->Branch("vvp3HS",   event.vvp3HS, "vvp3HS[ntK18]/D");

  tree->Branch("xvp4HS",   event.xvp4HS, "xvp4HS[ntK18]/D");
  tree->Branch("yvp4HS",   event.yvp4HS, "yvp4HS[ntK18]/D");
  tree->Branch("zvp4HS",   event.zvp4HS, "zvp4HS[ntK18]/D");
  tree->Branch("uvp4HS",   event.uvp4HS, "uvp4HS[ntK18]/D");
  tree->Branch("vvp4HS",   event.vvp4HS, "vvp4HS[ntK18]/D");

  tree->Branch("xhtofHS",  event.xhtofHS, "xhtofHS[ntK18]/D");
  tree->Branch("yhtofHS",  event.yhtofHS, "yhtofHS[ntK18]/D");
  tree->Branch("zhtofHS",  event.zhtofHS, "zhtofHS[ntK18]/D");
  tree->Branch("uhtofHS",  event.uhtofHS, "uhtofHS[ntK18]/D");
  tree->Branch("vhtofHS",  event.vhtofHS, "vhtofHS[ntK18]/D");

  tree->Branch("xtgtHS",   event.xtgtHS, "xtgtHS[ntK18]/D");
  tree->Branch("ytgtHS",   event.ytgtHS, "ytgtHS[ntK18]/D");
  tree->Branch("ztgtHS",   event.ztgtHS, "ztgtHS[ntK18]/D");
  tree->Branch("utgtHS",   event.utgtHS, "utgtHS[ntK18]/D");
  tree->Branch("vtgtHS",   event.vtgtHS, "vtgtHS[ntK18]/D");
  tree->Branch("pHS"   ,   event.pHS,    "pHS[ntK18]/D");
  tree->Branch("thetaHS",  event.thetaHS,"thetaHS[ntK18]/D");
  tree->Branch("phiHS",    event.phiHS,  "phiHS[ntK18]/D");
  tree->Branch("pathHS",   event.pathHS, "pathHS[ntK18]/D");
  tree->Branch("m2",       event.m2, "m2[ntK18]/D");
  tree->Branch("initmomHS", event.initmomHS, "initmomHS[ntK18]/D");
  HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")        &&
     InitializeParameter<DCDriftParamMan>("DCDRFT") &&
     InitializeParameter<DCTdcCalibMan>("DCTDC")    &&
     InitializeParameter<HodoParamMan>("HDPRM")     &&
     InitializeParameter<HodoPHCMan>("HDPHC")       &&
     InitializeParameter<BH2Filter>("BH2FLT")       &&
     InitializeParameter<BH1Match>("BH1MTH")        &&
     InitializeParameter<K18TransMatrix>("K18TM")   &&
     InitializeParameter<FieldMan>("FLDMAP", "HSFLDMAP") &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

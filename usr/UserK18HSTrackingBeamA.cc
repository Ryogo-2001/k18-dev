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
class UserK18HSTrackingBeamA : public VEvent
{
private:
  RawData*      rawData;
  HodoAnalyzer* hodoAna;
  DCAnalyzer*   DCAna;

public:
  UserK18HSTrackingBeamA();
  ~UserK18HSTrackingBeamA();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserK18HSTrackingBeamA::ClassName()
{
  static TString s_name("UserK18HSTrackingBeamABeamA");
  return s_name;
}

//_____________________________________________________________________________
UserK18HSTrackingBeamA::UserK18HSTrackingBeamA()
  : VEvent(),
    rawData(new RawData),
    hodoAna(new HodoAnalyzer),
    DCAna(new DCAnalyzer)
{
}

//_____________________________________________________________________________
UserK18HSTrackingBeamA::~UserK18HSTrackingBeamA()
{
  if(rawData) delete rawData;
  if(hodoAna) delete hodoAna;
  if(DCAna) delete DCAna;
}

//_____________________________________________________________________________
struct Event
{
  Int_t runnum;
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
  runnum         =  0;
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
UserK18HSTrackingBeamA::ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserK18HSTrackingBeamA::ProcessingNormal()
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

  event.runnum = gUnpacker.get_run_number();
  event.evnum = gUnpacker.get_event_number();

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
      }
    }
  }

  if(trigger_flag[trigger::kSpillEnd]) return true;

  if(event.trigflag[14] < 0) return false; //Select BeamA events only

  ////////// BH2 time 0
  hodoAna->DecodeBH2Hits(rawData);
#if HodoCut
  Int_t nhBh2 = hodoAna->GetNHitsBH2();
  if(nhBh2==0) return true;
#endif

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

  ////////// BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData);
#if HodoCut
  Int_t nhBh1 = hodoAna->GetNHitsBH1();
  if(nhBh1==0) return true;
#endif

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
    }
#if TimeCut
    hodoAna->TimeCutBFT(MinTimeBFT, MaxTimeBFT);
#endif
    Int_t ncl = hodoAna->GetNClustersBFT();
    event.bft_ncl = ncl;
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
    }

    event.bft_ncl_bh1mth = xCand.size();
  }

  if(xCand.size()!=1) return true;

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

  for(Int_t it=0; it<ntBcOut; ++it){
    DCLocalTrack *tp = DCAna->GetTrackBcOut(it);
    Int_t nh = tp->GetNHit();
    Double_t chisqr = tp->GetChiSquare();
    Double_t u0 = tp->GetU0(),  v0 = tp->GetV0();
    Double_t x0 = tp->GetX(0.), y0 = tp->GetY(0.);

    event.nhBcOut[it] = nh;
    event.chisqrBcOut[it] = chisqr;
    event.x0BcOut[it] = x0;
    event.y0BcOut[it] = y0;
    event.u0BcOut[it] = u0;
    event.v0BcOut[it] = v0;

    for(Int_t ih=0; ih<nh; ++ih){
      DCLTrackHit *hit=tp->GetHit(ih);
      Int_t layerId=hit->GetLayer()-100;
    }
  }
  if(ntBcOut==0) return true;

  ////////// K18HSTracking D2U
  DCAna->TrackSearchK18D2U(xCand);
  Int_t ntK18 = DCAna->GetNTracksK18D2U();
  if(ntK18 > MaxHits){
    std::cout << "#W too many ntK18 "
	      << ntK18 << "/" << MaxHits << std::endl;
    ntK18 = MaxHits;
  }
  event.ntK18 = ntK18;
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

  if(event.ntK18==1){
    HF1(10, 1);
    HF1(11, event.xtgtHS[0]);
    HF1(12, event.xvp4HS[0]);
    HF1(13, event.ytgtHS[0]);
    HF1(14, event.yvp4HS[0]);
    HF2(15, event.xtgtHS[0], event.ytgtHS[0]);
    HF2(16, event.xvp4HS[0], event.yvp4HS[0]);
    if(event.trigflag[20] > 0){
      HF1(20, 1);
      HF1(21, event.xtgtHS[0]);
      HF1(22, event.xvp4HS[0]);
      HF1(23, event.ytgtHS[0]);
      HF1(24, event.yvp4HS[0]);
      HF2(25, event.xtgtHS[0], event.ytgtHS[0]);
      HF2(26, event.xvp4HS[0], event.yvp4HS[0]);
    }
    if(event.trigflag[21] > 0){
      HF1(30, 1);
      HF1(31, event.xtgtHS[0]);
      HF1(32, event.xvp4HS[0]);
      HF1(33, event.ytgtHS[0]);
      HF1(34, event.yvp4HS[0]);
      HF2(35, event.xtgtHS[0], event.ytgtHS[0]);
      HF2(36, event.xvp4HS[0], event.yvp4HS[0]);
    }
    if(event.trigflag[22] > 0){
      HF1(40, 1);
      HF1(41, event.xtgtHS[0]);
      HF1(42, event.xvp4HS[0]);
      HF1(43, event.ytgtHS[0]);
      HF1(44, event.yvp4HS[0]);
      HF2(45, event.xtgtHS[0], event.ytgtHS[0]);
      HF2(46, event.xvp4HS[0], event.yvp4HS[0]);
    }
    if(event.trigflag[23] > 0){
      HF1(50, 1);
      HF1(51, event.xtgtHS[0]);
      HF1(52, event.xvp4HS[0]);
      HF1(53, event.ytgtHS[0]);
      HF1(54, event.yvp4HS[0]);
      HF2(55, event.xtgtHS[0], event.ytgtHS[0]);
      HF2(56, event.xvp4HS[0], event.yvp4HS[0]);
    }
  }

  return true;

}

//_____________________________________________________________________________
Bool_t
UserK18HSTrackingBeamA::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserK18HSTrackingBeamA;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{

  HB1(10, "Trig All", 2, 0, 2);
  HB1(11, "Xtgt K18Track", 200, -100., 100.);
  HB1(12, "Ytgt K18Track", 200, -100., 100.);
  HB1(13, "Upstream Xtgt K18Track", 200, -100., 100.);
  HB1(14, "Upstream Ytgt K18Track", 200, -100., 100.);
  HB2(15, "Y%X@Tgt", 400, -100., 100., 400, -100., 100.);
  HB2(16, "Y%X@TgtUpstream", 400, -100., 100., 400, -100., 100.);

  HB1(20, "Trig A", 2, 0, 2);
  HB1(21, "Xtgt K18Track Trig.A", 200, -100., 100.);
  HB1(22, "Ytgt K18Track Trig.A", 200, -100., 100.);
  HB1(23, "Upstream Xtgt K18Track Trig.A", 200, -100., 100.);
  HB1(24, "Upstream Ytgt K18Track Trig.A", 200, -100., 100.);
  HB2(25, "Y%X@Tgt Trig.A", 400, -100., 100., 400, -100., 100.);
  HB2(26, "Y%X@TgtUpstream Trig.A", 400, -100., 100., 400, -100., 100.);

  HB1(30, "Trig B", 2, 0, 2);
  HB1(31, "Xtgt K18Track Trig.B", 200, -100., 100.);
  HB1(32, "Ytgt K18Track Trig.B", 200, -100., 100.);
  HB1(33, "Upstream Xtgt K18Track Trig.B", 200, -100., 100.);
  HB1(34, "Upstream Ytgt K18Track Trig.B", 200, -100., 100.);
  HB2(35, "Y%X@Tgt Trig.B", 400, -100., 100., 400, -100., 100.);
  HB2(36, "Y%X@TgtUpstream Trig.B", 400, -100., 100., 400, -100., 100.);

  HB1(40, "Trig C", 2, 0, 2);
  HB1(41, "Xtgt K18Track Trig.C", 200, -100., 100.);
  HB1(42, "Ytgt K18Track Trig.C", 200, -100., 100.);
  HB1(43, "Upstream Xtgt K18Track Trig.C", 200, -100., 100.);
  HB1(44, "Upstream Ytgt K18Track Trig.C", 200, -100., 100.);
  HB2(45, "Y%X@Tgt Trig.C", 400, -100., 100., 400, -100., 100.);
  HB2(46, "Y%X@TgtUpstream Trig.C", 400, -100., 100., 400, -100., 100.);

  HB1(50, "Trig D", 2, 0, 2);
  HB1(51, "Xtgt K18Track Trig.D", 200, -100., 100.);
  HB1(52, "Ytgt K18Track Trig.D", 200, -100., 100.);
  HB1(53, "Upstream Xtgt K18Track Trig.D", 200, -100., 100.);
  HB1(54, "Upstream Ytgt K18Track Trig.D", 200, -100., 100.);
  HB2(55, "Y%X@Tgt Trig.D", 400, -100., 100., 400, -100., 100.);
  HB2(56, "Y%X@TgtUpstream Trig.D", 400, -100., 100., 400, -100., 100.);

  //tree
  HBTree("k18track","Data Summary Table of K18HSTracking(BeamA)");
  tree->Branch("runnum", &event.runnum, "runnum/I");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("ntK18",      &event.ntK18,     "ntK18/I");
  tree->Branch("chisqrK18",   event.chisqrK18, "chisqrK18[ntK18]/D");
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));
  tree->Branch("xvpHS",   event.xvp4HS, "xvpHS[ntK18]/D");
  tree->Branch("yvpHS",   event.yvp4HS, "yvpHS[ntK18]/D");
  tree->Branch("zvpHS",   event.zvp4HS, "zvpHS[ntK18]/D");
  tree->Branch("xtgtHS",   event.xtgtHS, "xtgtHS[ntK18]/D");
  tree->Branch("ytgtHS",   event.ytgtHS, "ytgtHS[ntK18]/D");
  tree->Branch("ztgtHS",   event.ztgtHS, "ztgtHS[ntK18]/D");

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

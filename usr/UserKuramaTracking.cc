// -*- C++ -*-

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>

#include <TMath.h>

#include "ConfMan.hh"
#include "DatabasePDG.hh"
#include "DCRawHit.hh"
#include "DebugTimer.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "RMAnalyzer.hh"
#include "KuramaLib.hh"
#include "MathTools.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UnpackerManager.hh"

#define HodoCut 0
#define UseTOF  1
#define Chi2Cut 1 // for Kurama tracking

namespace
{
using namespace root;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& zTOF = gGeom.LocalZ("TOF");

const Int_t &IdHTOF = gGeom.DetectorId("VPHTOF");
const Double_t &zHTOF = gGeom.LocalZ("VPHTOF");
const Int_t &IdVPTPC1 = gGeom.DetectorId("VPTPC1");
const Double_t& zVPTPC1 = gGeom.LocalZ("VPTPC1");
const Double_t& zVPTPC2 = gGeom.LocalZ("VPTPC2");
const Double_t& zVPTPC3 = gGeom.LocalZ("VPTPC3");
const Double_t& zVPTPC4 = gGeom.LocalZ("VPTPC4");
const Double_t& zVPTPC5 = gGeom.LocalZ("VPTPC5");
const Int_t &IdVP1 = gGeom.DetectorId("VP1");
const Int_t &IdWC = gGeom.DetectorId("WC");
}

//_____________________________________________________________________________
class UserKuramaTracking : public VEvent
{
private:
  RawData*      rawData;
  HodoAnalyzer* hodoAna;
  DCAnalyzer*   DCAna;

public:
  UserKuramaTracking();
  ~UserKuramaTracking();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserKuramaTracking::ClassName()
{
  static TString s_name("UserKuramaTracking");
  return s_name;
}

//_____________________________________________________________________________
UserKuramaTracking::UserKuramaTracking()
  : VEvent(),
    rawData(new RawData),
    hodoAna(new HodoAnalyzer),
    DCAna(new DCAnalyzer)
{
}

//_____________________________________________________________________________
UserKuramaTracking::~UserKuramaTracking()
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

  Double_t btof;
  Double_t time0;

  Int_t nhBh2;
  Double_t Bh2Seg[MaxHits];
  Double_t tBh2[MaxHits];
  Double_t t0Bh2[MaxHits];
  Double_t deBh2[MaxHits];

  Int_t nhBh1;
  Double_t Bh1Seg[MaxHits];
  Double_t tBh1[MaxHits];
  Double_t deBh1[MaxHits];

  Int_t nhTof;
  Double_t TofSeg[MaxHits];
  Double_t tTof[MaxHits];
  Double_t dtTof[MaxHits];
  Double_t deTof[MaxHits];

  Int_t ntSdcIn;
  Int_t much; // debug
  Int_t nlSdcIn;
  Int_t nhSdcIn[MaxHits];
  Double_t wposSdcIn[NumOfLayersSdcIn];
  Double_t chisqrSdcIn[MaxHits];
  Double_t u0SdcIn[MaxHits];
  Double_t v0SdcIn[MaxHits];
  Double_t x0SdcIn[MaxHits];
  Double_t y0SdcIn[MaxHits];

  Int_t ntSdcOut;
  Int_t nlSdcOut;
  Int_t nhSdcOut[MaxHits];
  Double_t wposSdcOut[NumOfLayersSdcOut];
  Double_t chisqrSdcOut[MaxHits];
  Double_t u0SdcOut[MaxHits];
  Double_t v0SdcOut[MaxHits];
  Double_t x0SdcOut[MaxHits];
  Double_t y0SdcOut[MaxHits];

  Int_t ntKurama;
  Int_t nlKurama;
  Int_t nhKurama[MaxHits];
  Double_t chisqrKurama[MaxHits];
  Double_t path[MaxHits];
  Double_t stof[MaxHits];
  Double_t pKurama[MaxHits];
  Double_t qKurama[MaxHits];
  Double_t m2[MaxHits];
  Double_t resP[MaxHits];
  Double_t vpx[NumOfLayersVP];
  Double_t vpy[NumOfLayersVP];
  Double_t vpxtpc[MaxHits][NumOfLayersVP];
  Double_t vpytpc[MaxHits][NumOfLayersVP];
  Double_t vpztpc[MaxHits][NumOfLayersVP];
  Double_t vpxhtof[MaxHits];
  Double_t vpyhtof[MaxHits];
  Double_t vpzhtof[MaxHits];

  Double_t xtgtKurama[MaxHits];
  Double_t ytgtKurama[MaxHits];
  Double_t utgtKurama[MaxHits];
  Double_t vtgtKurama[MaxHits];
  Double_t thetaKurama[MaxHits];
  Double_t phiKurama[MaxHits];

  Double_t xtofKurama[MaxHits];
  Double_t ytofKurama[MaxHits];
  Double_t utofKurama[MaxHits];
  Double_t vtofKurama[MaxHits];
  Double_t tofsegKurama[MaxHits];

  Double_t xwcKurama[MaxHits];
  Double_t ywcKurama[MaxHits];
  Double_t pathwcKurama[MaxHits];

  std::vector< std::vector<Double_t> > resL;
  std::vector< std::vector<Double_t> > resG;

  // Calib
  enum eParticle { Pion, Kaon, Proton, nParticle };
  Double_t tTofCalc[nParticle];
  Double_t utTofSeg[NumOfSegTOF];
  Double_t dtTofSeg[NumOfSegTOF];
  Double_t udeTofSeg[NumOfSegTOF];
  Double_t ddeTofSeg[NumOfSegTOF];
  Double_t tofua[NumOfSegTOF];
  Double_t tofda[NumOfSegTOF];

  // For TPC + Kurama RK
  Int_t nh[MaxHits];
  Double_t xin[MaxHits];
  Double_t yin[MaxHits];
  Double_t zin[MaxHits];
  Double_t pxin[MaxHits];
  Double_t pyin[MaxHits];
  Double_t pzin[MaxHits];
  Double_t xout[MaxHits];
  Double_t yout[MaxHits];
  Double_t zout[MaxHits];
  Double_t pxout[MaxHits];
  Double_t pyout[MaxHits];
  Double_t pzout[MaxHits];
  Double_t xtof[MaxHits];
  Double_t ytof[MaxHits];
  Double_t ztof[MaxHits];
  Double_t pxtof[MaxHits];
  Double_t pytof[MaxHits];
  Double_t pztof[MaxHits];
  Int_t layer[MaxHits][22];
  Double_t wire[MaxHits][22];
  Double_t localhitpos[MaxHits][22];
  Double_t wpos[MaxHits][22];

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  ntSdcIn  = 0;
  nlSdcIn  = 0;
  ntSdcOut = 0;
  nlSdcOut = 0;
  ntKurama = 0;
  nlKurama = 0;
  nhBh2    = 0;
  nhBh1    = 0;
  nhTof    = 0;
  much     = -1;

  time0 = qnan;
  btof  = qnan;

  for(Int_t i = 0; i<NumOfLayersVP; ++i){
    vpx[i] = qnan;
    vpy[i] = qnan;
  }

  for(Int_t it=0; it<MaxHits; ++it){
    vpxhtof[it] = qnan;
    vpyhtof[it] = qnan;
    vpzhtof[it] = qnan;
    for(Int_t il = 0; il<NumOfLayersVP; ++il){
      vpxtpc[it][il] = qnan;
      vpytpc[it][il] = qnan;
      vpztpc[it][il] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegTrig; it++){
    trigpat[it] = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<MaxHits; it++){
    Bh2Seg[it] = -1;
    tBh2[it] = qnan;
    deBh2[it] = qnan;
    Bh1Seg[it] = -1;
    tBh1[it] = qnan;
    deBh1[it] = qnan;
    TofSeg[it] = -1;
    tTof[it] = qnan;
    deTof[it] = qnan;
  }

  for(Int_t it=0; it<NumOfLayersSdcIn; ++it){
    wposSdcIn[it] = qnan;
  }
  for(Int_t it=0; it<NumOfLayersSdcOut; ++it){
    wposSdcOut[it] = qnan;
  }

  for(Int_t it=0; it<MaxHits; it++){
    nhSdcIn[it] = 0;
    chisqrSdcIn[it] = qnan;
    x0SdcIn[it] = qnan;
    y0SdcIn[it] = qnan;
    u0SdcIn[it] = qnan;
    v0SdcIn[it] = qnan;
  }

  for(Int_t it=0; it<MaxHits; it++){
    nhSdcOut[it] = 0;
    chisqrSdcOut[it] = qnan;
    x0SdcOut[it] = qnan;
    y0SdcOut[it] = qnan;
    u0SdcOut[it] = qnan;
    v0SdcOut[it] = qnan;
  }

  for(Int_t it=0; it<MaxHits; it++){
    nhKurama[it]     = 0;
    chisqrKurama[it] = qnan;
    stof[it]         = qnan;
    path[it]         = qnan;
    pKurama[it]      = qnan;
    qKurama[it]      = qnan;
    m2[it]           = qnan;
    xtgtKurama[it]  = qnan;
    ytgtKurama[it]  = qnan;
    utgtKurama[it]  = qnan;
    vtgtKurama[it]  = qnan;
    thetaKurama[it] = qnan;
    phiKurama[it]   = qnan;
    resP[it]        = qnan;
    xtofKurama[it]  = qnan;
    ytofKurama[it]  = qnan;
    utofKurama[it]  = qnan;
    vtofKurama[it]  = qnan;
    tofsegKurama[it] = qnan;
    xwcKurama[it]   = qnan;
    ywcKurama[it]   = qnan;
    pathwcKurama[it] = qnan;
  }

  for(Int_t it=0; it<MaxHits; it++){
    nh[it] = qnan;
    xin[it] = qnan;
    yin[it] = qnan;
    zin[it] = qnan;
    pxin[it] = qnan;
    pyin[it] = qnan;
    pzin[it] = qnan;
    xout[it] = qnan;
    yout[it] = qnan;
    zout[it] = qnan;
    pxout[it] = qnan;
    pyout[it] = qnan;
    pzout[it] = qnan;
    xtof[it] = qnan;
    ytof[it] = qnan;
    ztof[it] = qnan;
    pxtof[it] = qnan;
    pytof[it] = qnan;
    pztof[it] = qnan;

    for(Int_t i=0; i<NumOfLayersSdcIn+NumOfLayersSdcOut+4; ++i){
      layer[it][i] = -1;
      wire[it][i] = qnan;
      localhitpos[it][i] = qnan;
      wpos[it][i] = qnan;
    }
  }
  for(Int_t i=0; i<NumOfLayersSdcIn+NumOfLayersSdcOut+4; ++i){
    resL[i].clear();
    resG[i].clear();
  }

  for(Int_t i=0; i<Event::nParticle; ++i){
    tTofCalc[i] = qnan;
  }

  for(Int_t i=0; i<NumOfSegTOF; ++i){
    // tofmt[i] = qnan;
    utTofSeg[i]  = qnan;
    dtTofSeg[i]  = qnan;
    // ctuTofSeg[i] = qnan;
    // ctdTofSeg[i] = qnan;
    // ctTofSeg[i]  = qnan;
    udeTofSeg[i] = qnan;
    ddeTofSeg[i] = qnan;
    // deTofSeg[i]  = qnan;
    tofua[i]     = qnan;
    tofda[i]     = qnan;
  }

}

//_____________________________________________________________________________
namespace root
{
Event  event;
TH1   *h[MaxHist];
TTree *tree;
}


//_____________________________________________________________________________
Bool_t
UserKuramaTracking::ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserKuramaTracking::ProcessingNormal()
{
#if HodoCut
  static const auto MinDeBH2   = gUser.GetParameter("DeBH2", 0);
  static const auto MaxDeBH2   = gUser.GetParameter("DeBH2", 1);
  static const auto MinDeBH1   = gUser.GetParameter("DeBH1", 0);
  static const auto MaxDeBH1   = gUser.GetParameter("DeBH1", 1);
  static const auto MinBeamToF = gUser.GetParameter("BTOF", 0);
  static const auto MaxBeamToF = gUser.GetParameter("BTOF", 1);
  static const auto MinDeTOF   = gUser.GetParameter("DeTOF", 0);
  static const auto MaxDeTOF   = gUser.GetParameter("DeTOF", 1);
  static const auto MinTimeTOF = gUser.GetParameter("TimeTOF", 0);
  static const auto MaxTimeTOF = gUser.GetParameter("TimeTOF", 1);
#endif

  static const auto MinStopTimingSdcOut = gUser.GetParameter("StopTimingSdcOut", 0);
  static const auto MaxStopTimingSdcOut = gUser.GetParameter("StopTimingSdcOut", 1);
  // static const auto StopTimeDiffSdcOut = gUser.GetParameter("StopTimeDiffSdcOut");

  static const auto MaxMultiHitSdcIn  = gUser.GetParameter("MaxMultiHitSdcIn");
  static const auto MaxMultiHitSdcOut = gUser.GetParameter("MaxMultiHitSdcOut");

  static const auto MinTotSDC1 = gUser.GetParameter("MinTotSDC1");
  static const auto MinTotSDC2 = gUser.GetParameter("MinTotSDC2");
  static const auto MinTotSDC3 = gUser.GetParameter("MinTotSDC3");
  static const auto MinTotSDC4 = gUser.GetParameter("MinTotSDC4");

  static const auto MinDriftTimeSDC1 = gUser.GetParameter("DriftTimeSDC1", 0);
  static const auto MaxDriftTimeSDC1 = gUser.GetParameter("DriftTimeSDC1", 1);
  static const auto MinDriftTimeSDC2 = gUser.GetParameter("DriftTimeSDC2", 0);
  static const auto MaxDriftTimeSDC2 = gUser.GetParameter("DriftTimeSDC2", 1);
  static const auto MinDriftTimeSDC3 = gUser.GetParameter("DriftTimeSDC3", 0);
  static const auto MaxDriftTimeSDC3 = gUser.GetParameter("DriftTimeSDC3", 1);
  static const auto MinDriftTimeSDC4 = gUser.GetParameter("DriftTimeSDC4", 0);
  static const auto MaxDriftTimeSDC4 = gUser.GetParameter("DriftTimeSDC4", 1);

  static const auto MaxChisqrSdcIn = gUser.GetParameter("MaxChisqrSdcIn");
  static const auto MaxChisqrSdcOut = gUser.GetParameter("MaxChisqrSdcOut");
  static const auto MaxChisqrKurama = gUser.GetParameter("MaxChisqrKurama");
  rawData->DecodeHits();

  event.evnum = gUnpacker.get_event_number();

  Double_t common_stop_tdc = qnan;

  // Trigger Flag
  std::bitset<NumOfSegTrig> trigger_flag;
  for(const auto& hit: rawData->GetTrigRawHC()){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc1();
    if(tdc > 0){
      event.trigpat[trigger_flag.count()] = seg;
      event.trigflag[seg] = tdc;
      trigger_flag.set(seg);
      if(seg == trigger::kCommonStopSdcOut){
        common_stop_tdc = tdc;
      }
    }
  }

  HF1(1, 0.);

  if(trigger_flag[trigger::kSpillEnd]) return true;

  HF1(1, 1.);

  //////////////BH2 Analysis
  hodoAna->DecodeBH2Hits(rawData);
  Int_t nhBh2 = hodoAna->GetNHitsBH2();
  event.nhBh2 = nhBh2;
#if HodoCut
  if(nhBh2==0) return true;
#endif
  Double_t time0 = -999.;
  //////////////BH2 Analysis
  Double_t min_time = -999.;
  for(Int_t i=0; i<nhBh2; ++i){
    BH2Hit *hit = hodoAna->GetHitBH2(i);
    if(!hit) continue;
    Double_t seg = hit->SegmentId()+1;
    Double_t mt  = hit->MeanTime();
    Double_t cmt = hit->CMeanTime();
    Double_t ct0 = hit->CTime0();
    Double_t de  = hit->DeltaE();
#if HodoCut
    if(de<MinDeBH2 || MaxDeBH2<de) continue;
#endif
    event.tBh2[i]   = cmt;
    event.t0Bh2[i]  = ct0;
    event.deBh2[i]  = de;
    event.Bh2Seg[i] = seg;
    if(std::abs(mt)<std::abs(min_time)){
      min_time = mt;
      time0    = ct0;
    }
  }
  BH2Cluster* cl_time0 = hodoAna->GetTime0BH2Cluster();
  if(cl_time0) time0 = cl_time0->CTime0();
  event.time0 = time0;

  //////////////BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData);
  Int_t nhBh1 = hodoAna->GetNHitsBH1();
  event.nhBh1 = nhBh1;
#if HodoCut
  if(nhBh1==0) return true;
#endif

  HF1(1, 2.);

  Double_t btof0 = -999.;
  for(Int_t i=0; i<nhBh1; ++i){
    Hodo2Hit *hit = hodoAna->GetHitBH1(i);
    if(!hit) continue;
    Int_t    seg  = hit->SegmentId()+1;
    Double_t cmt  = hit->CMeanTime();
    Double_t dE   = hit->DeltaE();
    Double_t btof = cmt - time0;
#if HodoCut
    if(dE<MinDeBH1 || MaxDeBH1<dE) continue;
    if(btof<MinBeamToF && MaxBeamToF<btof) continue;
#endif
    event.Bh1Seg[i] = seg;
    event.tBh1[i]   = cmt;
    event.deBh1[i]  = dE;
    if(std::abs(btof)<std::abs(btof0)){
      btof0 = btof;
    }
  }

  event.btof = btof0;

  HF1(1, 3.);

  HodoClusterContainer TOFCont;
  //////////////Tof Analysis
  hodoAna->DecodeTOFHits(rawData);
  //hodoAna->TimeCutTOF(7, 25);
  Int_t nhTof = hodoAna->GetNClustersTOF();
  event.nhTof = nhTof;
  {
#if HodoCut
    Int_t nhOk = 0;
#endif
    for(Int_t i=0; i<nhTof; ++i){
      HodoCluster *hit = hodoAna->GetClusterTOF(i);
      Double_t seg = hit->MeanSeg()+1;
      Double_t cmt = hit->CMeanTime();
      Double_t dt  = hit->TimeDif();
      Double_t de   = hit->DeltaE();
      event.TofSeg[i] = seg;
      event.tTof[i]   = cmt;//stof;
      event.dtTof[i]  = dt;
      event.deTof[i]  = de;
      // for PHC
      // HF2(100*seg+30000+81, ua, stof);
      // HF2(100*seg+30000+82, da, stof);
      // HF2(100*seg+30000+83, ua, (time0-OffsetTof)-ut);
      // HF2(100*seg+30000+84, da, (time0-OffsetTof)-dt);
      TOFCont.push_back(hit);
#if HodoCut
      if(MinDeTOF<de  && de<MaxDeTOF  &&
         MinTimeTOF<stof && stof<MaxTimeTOF){
	++nhOk;
      }
#endif
    }

#if HodoCut
    if(nhOk<1) return true;
#endif
  }

  HF1(1, 4.);

  // Common stop timing
  Bool_t common_stop_is_tof = (common_stop_tdc < MinStopTimingSdcOut
                               || MaxStopTimingSdcOut < common_stop_tdc);
  if(!common_stop_is_tof) return true;

  HF1(1, 6.);

  HF1(1, 10.);

  DCAna->DecodeSdcInHits(rawData);
  DCAna->TotCutSDC1(MinTotSDC1);
  DCAna->TotCutSDC2(MinTotSDC2);
  DCAna->DriftTimeCutSDC1(MinDriftTimeSDC1, MaxDriftTimeSDC1);
  DCAna->DriftTimeCutSDC2(MinDriftTimeSDC2, MaxDriftTimeSDC2);

  // Double_t offset = common_stop_is_tof ? 0 : StopTimeDiffSdcOut;
  DCAna->DecodeSdcOutHits(rawData);
  DCAna->TotCutSDC3(MinTotSDC3);
  DCAna->TotCutSDC4(MinTotSDC4);
  DCAna->DriftTimeCutSDC3(MinDriftTimeSDC3, MaxDriftTimeSDC3);
  DCAna->DriftTimeCutSDC4(MinDriftTimeSDC4, MaxDriftTimeSDC4);

  Double_t multi_SdcIn  = 0.;
  ////////////// SdcIn number of hit layer
  {
    Int_t nlSdcIn = 0;
    for(Int_t layer=1; layer<=NumOfLayersSdcIn; ++layer){
      const DCHitContainer &contSdcIn =DCAna->GetSdcInHC(layer);
      Int_t nhSdcIn = contSdcIn.size();
      if(nhSdcIn==1){
	DCHit *hit = contSdcIn[0];
	Double_t wpos = hit->GetWirePosition();
	event.wposSdcIn[layer-1] = wpos;
      }
      multi_SdcIn += Double_t(nhSdcIn);
      if(nhSdcIn>0) nlSdcIn++;
    }
    event.nlSdcIn   = nlSdcIn;
    event.nlKurama += nlSdcIn;
  }
  if(multi_SdcIn/Double_t(NumOfLayersSdcIn) > MaxMultiHitSdcIn){
    // return true;
  }

  Double_t multi_SdcOut = 0.;
  ////////////// SdcOut number of hit layer
  {
    Int_t nlSdcOut = 0;
    for(Int_t layer=1; layer<=NumOfLayersSdcOut; ++layer){
      const DCHitContainer &contSdcOut =DCAna->GetSdcOutHC(layer);
      Int_t nhSdcOut=contSdcOut.size();
      if(nhSdcOut==1){
	DCHit *hit = contSdcOut[0];
	Double_t wpos = hit->GetWirePosition();
	event.wposSdcOut[layer-1] = wpos;
      }
      multi_SdcOut += Double_t(nhSdcOut);
      if(nhSdcOut>0) nlSdcOut++;
    }
    event.nlSdcOut = nlSdcOut;
    event.nlKurama += nlSdcOut;
  }
  if(multi_SdcOut/Double_t(NumOfLayersSdcOut) > MaxMultiHitSdcOut){
    // return true;
  }

  HF1(1, 11.);

  // std::cout << "==========TrackSearch SdcIn============" << std::endl;
  DCAna->TrackSearchSdcIn();
  DCAna->ChiSqrCutSdcIn(MaxChisqrSdcIn);

  Int_t ntSdcIn = DCAna->GetNtracksSdcIn();
  if(MaxHits<ntSdcIn){
    std::cout << "#W " << FUNC_NAME << " "
	      << "too many ntSdcIn " << ntSdcIn << "/" << MaxHits << std::endl;
    ntSdcIn = MaxHits;
  }
  event.ntSdcIn = ntSdcIn;
  HF1(10, Double_t(ntSdcIn));
  Int_t much_combi = DCAna->MuchCombinationSdcIn();
  event.much = much_combi;
  for(Int_t it=0; it<ntSdcIn; ++it){
    DCLocalTrack *tp=DCAna->GetTrackSdcIn(it);
    Int_t nh=tp->GetNHit();
    Double_t chisqr=tp->GetChiSquare();
    Double_t x0=tp->GetX0(), y0=tp->GetY0();
    Double_t u0=tp->GetU0(), v0=tp->GetV0();
    HF1(11, Double_t(nh));
    HF1(12, chisqr);
    HF1(14, x0); HF1(15, y0);
    HF1(16, u0); HF1(17, v0);
    HF2(18, x0, u0); HF2(19, y0, v0);
    HF2(20, x0, y0);
    event.nhSdcIn[it] = nh;
    event.chisqrSdcIn[it] = chisqr;
    event.x0SdcIn[it] = x0;
    event.y0SdcIn[it] = y0;
    event.u0SdcIn[it] = u0;
    event.v0SdcIn[it] = v0;
    for(Int_t ih=0; ih<nh; ++ih){
      DCLTrackHit *hit=tp->GetHit(ih);

      Int_t layerId = hit->GetLayer();
      HF1(13, hit->GetLayer());
      Double_t wire = hit->GetWire();
      Double_t dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1(100*layerId+1, wire-0.5);
      HF1(100*layerId+2, dt);
      HF1(100*layerId+3, dl);
      Double_t xcal=hit->GetXcal(), ycal=hit->GetYcal();
      Double_t pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1(100*layerId+4, pos);
      HF1(100*layerId+5, res);
      HF2(100*layerId+6, pos, res);
      HF2(100*layerId+7, xcal, ycal);
      for(Int_t i=0; i<NumOfLayersSdcIn; ++i){
	if(i==layerId-1) event.resL[i].push_back(res);
      }
    }
  }
  if(ntSdcIn<1) return true;
  //  if(!(ntSdcIn==1)) return true;

  HF1(1, 12.);

#if 0
  //////////////SdcIn vs Tof cut for Proton event
  {
    Int_t ntOk=0;
    for(Int_t it=0; it<ntSdcIn; ++it){
      DCLocalTrack *tp=DCAna->GetTrackSdcIn(it);
      if(!tp) continue;

      Int_t nh=tp->GetNHit();
      Double_t chisqr=tp->GetChiSquare();
      Double_t u0=tp->GetU0(), v0=tp->GetV0();
      Double_t x0=tp->GetX0(), y0=tp->GetY0();

      Bool_t condTof=false;
      for(Int_t j=0; j<ncTof; ++j){
	HodoCluster *clTof=hodoAna->GetClusterTOF(j);
	if(!clTof || !clTof->GoodForAnalysis()) continue;
	Double_t ttof=clTof->CMeanTime()-time0;
	//------------------------Cut
	if(MinModTimeTof< ttof+14.0*u0 && ttof+14.0*u0 <MaxModTimeTof)
	  condTof=true;
      }
      if(condTof){
	++ntOk;
	for(Int_t j=0; j<ncTof; ++j){
	  HodoCluster *clTof=hodoAna->GetClusterTOF(j);
	  if(!clTof || !clTof->GoodForAnalysis()) continue;
	  Double_t ttof=clTof->CMeanTime()-time0;
	}
	// if(ntOk>0) tp->GoodForTracking(false);
      }
      else {
	// tp->GoodForTracking(false);
      }
    }
    // if(ntOk<1) return true;
  }
#endif

  HF1(1, 13.);

  //////////////SdcOut tracking
  //std::cout << "==========TrackSearch SdcOut============" << std::endl;

#if UseTOF
  DCAna->TrackSearchSdcOut(TOFCont);
#else
  DCAna->TrackSearchSdcOut();
#endif
  DCAna->ChiSqrCutSdcOut(MaxChisqrSdcOut);
  Int_t ntSdcOut = DCAna->GetNtracksSdcOut();
  if(MaxHits<ntSdcOut){
    std::cout << "#W " << FUNC_NAME << " "
	      << "too many ntSdcOut " << ntSdcOut << "/" << MaxHits << std::endl;
    ntSdcOut = MaxHits;
  }
  event.ntSdcOut=ntSdcOut;
  HF1(30, Double_t(ntSdcOut));
  for(Int_t it=0; it<ntSdcOut; ++it){
    DCLocalTrack *tp=DCAna->GetTrackSdcOut(it);
    Int_t nh=tp->GetNHit();
    Double_t chisqr=tp->GetChiSquare();
    Double_t u0=tp->GetU0(), v0=tp->GetV0();
    Double_t x0=tp->GetX0(), y0=tp->GetY0();
    HF1(31, Double_t(nh));
    HF1(32, chisqr);
    HF1(34, x0); HF1(35, y0);
    HF1(36, u0); HF1(37, v0);
    HF2(38, x0, u0); HF2(39, y0, v0);
    HF2(40, x0, y0);
    event.nhSdcOut[it] = nh;
    event.chisqrSdcOut[it] = chisqr;
    event.x0SdcOut[it] = x0;
    event.y0SdcOut[it] = y0;
    event.u0SdcOut[it] = u0;
    event.v0SdcOut[it] = v0;
    for(Int_t ih=0; ih<nh; ++ih){
      DCLTrackHit *hit=tp->GetHit(ih);
      Int_t layerId = hit->GetLayer();
      if(hit->GetLayer()>79) layerId -= 62;
      else if(hit->GetLayer()>40) layerId -= 22;
      else if(hit->GetLayer()>30) layerId -= 20;
      HF1(33, hit->GetLayer());
      Double_t wire=hit->GetWire();
      Double_t dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1(100*layerId+1, wire-0.5);
      HF1(100*layerId+2, dt);
      HF1(100*layerId+3, dl);
      Double_t xcal=hit->GetXcal(), ycal=hit->GetYcal();
      Double_t pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1(100*layerId+4, pos);
      HF1(100*layerId+5, res);
      HF2(100*layerId+6, pos, res);
      HF2(100*layerId+7, xcal, ycal);
      for(Int_t i=0; i<NumOfLayersSdcOut+4; ++i){
	if(i==layerId-NumOfLayersSdcIn-1){
	  event.resL[i+NumOfLayersSdcIn].push_back(res);
	}
      }
    }
  }

  if(ntSdcOut<1) return true;

  HF1(1, 14.);

  for(Int_t i1=0; i1<ntSdcIn; ++i1){
    DCLocalTrack *trSdcIn=DCAna->GetTrackSdcIn(i1);
    Double_t xin=trSdcIn->GetX0(), yin=trSdcIn->GetY0();
    Double_t uin=trSdcIn->GetU0(), vin=trSdcIn->GetV0();
    for(Int_t i2=0; i2<ntSdcOut; ++i2){
      DCLocalTrack *trSdcOut=DCAna->GetTrackSdcOut(i2);
      Double_t xout=trSdcOut->GetX0(), yout=trSdcOut->GetY0();
      Double_t uout=trSdcOut->GetU0(), vout=trSdcOut->GetV0();
      HF2(20001, xin, xout); HF2(20002, yin, yout);
      HF2(20003, uin, uout); HF2(20004, vin, vout);
    }
  }

  HF1(1, 20.);

  // if(ntSdcIn*ntSdcOut > 4) return true;

  HF1(1, 21.);

  ///// BTOF BH2-Target
  static const auto StofOffset = gUser.GetParameter("StofOffset");

  //////////////KURAMA Tracking
  DCAna->TrackSearchKurama();
#if Chi2Cut
  DCAna->ChiSqrCutKurama(MaxChisqrKurama);
#endif
  Int_t ntKurama = DCAna->GetNTracksKurama();
  if(MaxHits < ntKurama){
    std::cout << "#W " << FUNC_NAME << " "
	      << "too many ntKurama " << ntKurama << "/" << MaxHits << std::endl;
    ntKurama = MaxHits;
  }
  event.ntKurama = ntKurama;
  HF1(50, ntKurama);

  for(Int_t i=0; i<ntKurama; ++i){
    auto track = DCAna->GetKuramaTrack(i);
    if(!track) continue;
    // track->Print();
    Int_t nh = track->GetNHits();
    Double_t chisqr = track->ChiSquare();
    const auto& Pos = track->PrimaryPosition();
    const auto& Mom = track->PrimaryMomentum();
    // hddaq::cout << std::fixed
    // 		<< "Pos = " << Pos << std::endl
    // 		<< "Mom = " << Mom << std::endl;
    Double_t path = track->PathLengthToTOF();
    Double_t xt = Pos.x(), yt = Pos.y();
    Double_t p = Mom.Mag();
    Double_t q = track->Polarity();
    Double_t ut = Mom.x()/Mom.z(), vt = Mom.y()/Mom.z();
    Double_t cost = 1./TMath::Sqrt(1.+ut*ut+vt*vt);
    Double_t theta = TMath::ACos(cost)*TMath::RadToDeg();
    Double_t phi = TMath::ATan2(ut, vt);
    Double_t initial_momentum = track->GetInitialMomentum();
    HF1(51, nh);
    HF1(52, chisqr);
    HF1(54, xt); HF1(55, yt); HF1(56, ut); HF1(57,vt);
    HF2(58, xt, ut); HF2(59, yt, vt); HF2(60, xt, yt);
    HF1(61, p); HF1(62, path);
    event.nhKurama[i] = nh;
    event.chisqrKurama[i] = chisqr;
    event.path[i] = path;
    event.pKurama[i] = p;
    event.qKurama[i] = q;
    event.xtgtKurama[i] = xt;
    event.ytgtKurama[i] = yt;
    event.utgtKurama[i] = ut;
    event.vtgtKurama[i] = vt;
    event.thetaKurama[i] = theta;
    event.phiKurama[i] = phi;
    event.resP[i] = p - initial_momentum;
    if(ntKurama == 1){
      Double_t x, y;
      for(Int_t l = 0; l<NumOfLayersVP; ++l){
	track->GetTrajectoryLocalPosition(IdVP1 + l, x, y);
	event.vpx[l] = x;
	event.vpy[l] = y;
      }// for(l)
    }
    Double_t x, y;
    for(Int_t l = 0; l<NumOfLayersVPTPC; ++l){
      track->GetTrajectoryLocalPosition(IdVPTPC1 + l, x, y);
      event.vpxtpc[i][l] = x;
      event.vpytpc[i][l] = y;
    }// for(l)
    event.vpztpc[i][0] = zVPTPC1;
    event.vpztpc[i][1] = zVPTPC2;
    event.vpztpc[i][2] = zVPTPC3;
    event.vpztpc[i][3] = zVPTPC4;
    event.vpztpc[i][4] = zVPTPC5;

    track->GetTrajectoryLocalPosition(IdHTOF, x, y);
    event.vpxhtof[i] = x;
    event.vpyhtof[i] = y;
    event.vpzhtof[i] = zHTOF;

    Double_t pathwc;
    track->GetTrajectoryLocalPosition(IdWC, pathwc, x, y);
    event.xwcKurama[i] = x;
    event.ywcKurama[i] = y;
    event.pathwcKurama[i] = pathwc;

    const auto& posTof = track->TofPos();
    const auto& momTof = track->TofMom();
    event.xtofKurama[i] = posTof.x();
    event.ytofKurama[i] = posTof.y();
    event.utofKurama[i] = momTof.x()/momTof.z();
    event.vtofKurama[i] = momTof.y()/momTof.z();
#if UseTOF
    Double_t tof_seg = track->TofSeg();
    event.tofsegKurama[i] = tof_seg;
#else
    Double_t tof_x = track->GetLocalTrackOut()->GetX(zTOF);
    Double_t tof_seg = MathTools::Round(tof_x/75. + (NumOfSegTOF + 1)*0.5);
    event.tofsegKurama[i] = tof_seg;
# if 0
    std::cout << "posTof " << posTof << std::endl;
    std::cout << "momTof " << momTof << std::endl;
    std::cout << std::setw(10) << vecTof.X()
	      << std::setw(10) << vecTof.Y()
	      << std::setw(10) << sign*vecTof.Mod()
	      << std::setw(10) << TofSegKurama << std::endl;
# endif
    Double_t minres = 1.0e10;
#endif
    Double_t time = qnan;
    for(const auto& hit: hodoAna->GetHitsTOF()){
      if(!hit) continue;
      Int_t seg = hit->SegmentId() + 1;
#if UseTOF
      if((Int_t)tof_seg == seg){
	time = hit->CMeanTime() - time0 + StofOffset;
      }
#else
      Double_t res = TMath::Abs(tof_seg - seg);
      if(res < minres){
	minres = res;
	time = hit->CMeanTime() - time0 + StofOffset;
      }
#endif
    }
    event.stof[i] = time;
    if(time > 0.){
      Double_t m2 = Kinematics::MassSquare(p, path, time);
      HF1(63, m2);
      event.m2[i] = m2;
# if 0
      std::ios::fmtflags pre_flags     = std::cout.flags();
      std::size_t        pre_precision = std::cout.precision();
      std::cout.setf(std::ios::fixed);
      std::cout.precision(5);
      std::cout << FUNC_NAME << std::endl
		<< "   Mom  = " << p     << std::endl
		<< "   Path = " << path << std::endl
		<< "   Time = " << time  << std::endl
		<< "   m2   = " << m2    << std::endl;
      std::cout.flags(pre_flags);
      std::cout.precision(pre_precision);
# endif
    }

    static const auto gposTOFDX = gGeom.GetGlobalPosition("TOF-DX");
    const Double_t xOut = track->GetLocalTrackOut()->GetX(gposTOFDX.z());
    const Double_t yOut = track->GetLocalTrackOut()->GetY(gposTOFDX.z());
    const Double_t uOut = track->GetLocalTrackOut()->GetU0();
    const Double_t vOut = track->GetLocalTrackOut()->GetV0();
    const Double_t pzOut = track->GetInitialMomentum()/std::sqrt(1.+uOut*uOut+vOut*vOut);
    const ThreeVector posOut(xOut, yOut, gposTOFDX.z());
    const ThreeVector momOut(pzOut*uOut, pzOut*vOut, pzOut);
# if 0
    std::cout<<" evnum : "<<event.evnum<<std::endl;
    std::cout<<" kurama "<<i<<"/"<<ntKurama
	     <<" q : "<<track->Polarity()
	     <<" mom : "<<Mom.x()<<" "<<Mom.y()<<" "<<Mom.z()<<" mag : "<<Mom.Mag()
	     <<" init mom : "<<track->GetInitialMomentum()<<" "<<momOut.Mag()<<std::endl;
    std::cout<<" posOut : "<<posOut.x()<<" "<<posOut.y()<<" "<<posOut.z()<<std::endl;
    std::cout<<" momOut : "<<momOut.x()<<" "<<momOut.y()<<" "<<momOut.z()<<std::endl;
    std::cout<<" posTof : "<<posTof.x()<<" "<<posTof.y()<<" "<<posTof.z()<<std::endl;
    std::cout<<" momTof : "<<momTof.x()<<" "<<momTof.y()<<" "<<momTof.z()<<std::endl;
# endif

    static const auto gposTgt = gGeom.GetGlobalPosition("Target");
    const Double_t xIn = track->GetLocalTrackIn()->GetX(gposTgt.z());
    const Double_t yIn = track->GetLocalTrackIn()->GetY(gposTgt.z());
    const Double_t uIn = track->GetLocalTrackIn()->GetU0();
    const Double_t vIn = track->GetLocalTrackIn()->GetV0();
    const Double_t pzIn = track->GetInitialMomentum()/std::sqrt(1.+uIn*uIn+vIn*vIn);
    const ThreeVector posIn(xIn, yIn, gposTgt.z());
    const ThreeVector momIn(pzIn*uIn, pzIn*vIn, pzIn);

    event.nh[i]=nh;
    event.xout[i]=posOut.x();
    event.yout[i]=posOut.y();
    event.zout[i]=posOut.z();
    event.pxout[i]=momOut.x();
    event.pyout[i]=momOut.y();
    event.pzout[i]=momOut.z();
    event.xtof[i]=posTof.x();
    event.ytof[i]=posTof.y();
    event.ztof[i]=posTof.z();
    event.pxtof[i]=momTof.x();
    event.pytof[i]=momTof.y();
    event.pztof[i]=momTof.z();
    event.xin[i]=posIn.x();
    event.yin[i]=posIn.y();
    event.zin[i]=posIn.z();
    event.pxin[i]=momIn.x();
    event.pyin[i]=momIn.y();
    event.pzin[i]=momIn.z();

    for(Int_t j=0; j<nh; ++j){
      TrackHit *hit=track->GetHit(j);
      if(!hit) continue;
      Int_t layerId = hit->GetLayer();
      if(hit->GetLayer()>79) layerId -= 62;
      else if(hit->GetLayer()>40) layerId -= 22;
      else if(hit->GetLayer()>30) layerId -= 20;
      HF1(53, hit->GetLayer());
      Double_t wire = hit->GetHit()->GetWire();
      Double_t dt   = hit->GetHit()->GetDriftTime();
      Double_t dl   = hit->GetHit()->GetDriftLength();
      Double_t pos  = hit->GetCalLPos();
      Double_t res  = hit->GetResidual();
      DCLTrackHit *lhit = hit->GetHit();
      Double_t xcal = lhit->GetXcal();
      Double_t ycal = lhit->GetYcal();
      HF1(100*layerId+11, Double_t(wire)-0.5);
      Double_t wp   = lhit->GetWirePosition();
      HF1(100*layerId+12, dt);
      HF1(100*layerId+13, dl);
      HF1(100*layerId+14, pos);

      event.layer[i][j] = hit->GetLayer();
      event.wire[i][j] = wire;
      event.localhitpos[i][j] = lhit -> GetLocalHitPos();
      event.wpos[i][j] = wp;
      if(nh>17 && q<=0. && chisqr<200.){
      	HF1(100*layerId+15, res);
	HF2(100*layerId+16, pos, res);
      }
      HF2(100*layerId+17, xcal, ycal);
      if (std::abs(dl-std::abs(xcal-wp))<2.0){
	HF2(100*layerId+22, dt, std::abs(xcal-wp));
      }
      for(Int_t l=0; l<NumOfLayersSdcIn; ++l){
	if(l==layerId-1)
	  event.resG[l].push_back(res);
      }
      for(Int_t l=0; l<NumOfLayersSdcOut+4; ++l){
	if(l==layerId-1-NumOfLayersSdcIn)
	  event.resG[l+NumOfLayersSdcIn].push_back(res);
      }
    }

    DCLocalTrack *trSdcIn  = track->GetLocalTrackIn();
    DCLocalTrack *trSdcOut = track->GetLocalTrackOut();

    if(trSdcIn){
      Int_t nhSdcIn=trSdcIn->GetNHit();
      Double_t x0in=trSdcIn->GetX0(), y0in=trSdcIn->GetY0();
      Double_t u0in=trSdcIn->GetU0(), v0in=trSdcIn->GetV0();
      Double_t chiin=trSdcIn->GetChiSquare();
      HF1(21, Double_t(nhSdcIn)); HF1(22, chiin);
      HF1(24, x0in); HF1(25, y0in); HF1(26, u0in); HF1(27, v0in);
      HF2(28, x0in, u0in); HF2(29, y0in, v0in);
      for(Int_t jin=0; jin<nhSdcIn; ++jin){
	Int_t layer=trSdcIn->GetHit(jin)->GetLayer();
	HF1(23, layer);
      }
    }
    if(trSdcOut){
      Int_t nhSdcOut=trSdcOut->GetNHit();
      Double_t x0out=trSdcOut->GetX(zTOF), y0out=trSdcOut->GetY(zTOF);
      Double_t u0out=trSdcOut->GetU0(), v0out=trSdcOut->GetV0();
      Double_t chiout=trSdcOut->GetChiSquare();
      HF1(41, Double_t(nhSdcOut)); HF1(42, chiout);
      HF1(44, x0out); HF1(45, y0out); HF1(46, u0out); HF1(47, v0out);
      HF2(48, x0out, u0out); HF2(49, y0out, v0out);
      for(Int_t jout=0; jout<nhSdcOut; ++jout){
	Int_t layer=trSdcOut->GetHit(jout)->GetLayer();
	HF1(43, layer);
      }
    }
  }

  for(Int_t i=0; i<ntKurama; ++i){
    KuramaTrack *track=DCAna->GetKuramaTrack(i);
    if(!track) continue;
    DCLocalTrack *trSdcIn =track->GetLocalTrackIn();
    DCLocalTrack *trSdcOut=track->GetLocalTrackOut();
    if(!trSdcIn || !trSdcOut) continue;
    Double_t yin=trSdcIn->GetY(500.), vin=trSdcIn->GetV0();
    Double_t yout=trSdcOut->GetY(3800.), vout=trSdcOut->GetV0();
    HF2(20021, yin, yout); HF2(20022, vin, vout);
    HF2(20023, vin, yout); HF2(20024, vout, yin);
  }

  if(ntKurama==0) return true;
  KuramaTrack *track = DCAna->GetKuramaTrack(0);
  Double_t path = track->PathLengthToTOF();
  Double_t p    = track->PrimaryMomentum().Mag();
  //Double_t tTof[Event::nParticle];
  Double_t calt[Event::nParticle] = {
    Kinematics::CalcTimeOfFlight(p, path, pdg::PionMass()),
    Kinematics::CalcTimeOfFlight(p, path, pdg::KaonMass()),
    Kinematics::CalcTimeOfFlight(p, path, pdg::ProtonMass())
  };
  for(Int_t i=0; i<Event::nParticle; ++i){
    event.tTofCalc[i] = calt[i];
  }

  // TOF
  {
    Int_t nh = hodoAna->GetNHitsTOF();
    for(Int_t i=0; i<nh; ++i){
      Hodo2Hit *hit=hodoAna->GetHitTOF(i);
      if(!hit) continue;
      Int_t seg=hit->SegmentId()+1;
      Double_t tu = hit->GetTUp(), td=hit->GetTDown();
      // Double_t ctu=hit->GetCTUp(), ctd=hit->GetCTDown();
      Double_t cmt=hit->CMeanTime();//, t= cmt-time0+StofOffset;//cmt-time0;
      Double_t ude=hit->GetAUp(), dde=hit->GetADown();
      // Double_t de=hit->DeltaE();
      // Double_t m2 = Kinematics::MassSquare(p, path, t);
      // event.tofmt[seg-1] = hit->MeanTime();
      event.utTofSeg[seg-1] = tu - time0 + StofOffset;
      event.dtTofSeg[seg-1] = td - time0 + StofOffset;
      // event.uctTofSeg[seg-1] = ctu - time0 + offset;
      // event.dctTofSeg[seg-1] = ctd - time0 + offset;
      event.udeTofSeg[seg-1] = ude;
      event.ddeTofSeg[seg-1] = dde;
      // event.ctTofSeg[seg-1]  = t;
      // event.deTofSeg[seg-1]  = de;
      HF2(30000+100*seg+83, ude, calt[Event::Pion]+time0-StofOffset-cmt);
      HF2(30000+100*seg+84, dde, calt[Event::Pion]+time0-StofOffset-cmt);
      HF2(30000+100*seg+83, ude, calt[Event::Pion]+time0-StofOffset-tu);
      HF2(30000+100*seg+84, dde, calt[Event::Pion]+time0-StofOffset-td);
    }

    const HodoRHitContainer &cont=rawData->GetTOFRawHC();
    Int_t NofHit = cont.size();
    for(Int_t i = 0; i<NofHit; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId();
      event.tofua[seg] = hit->GetAdcUp();
      event.tofda[seg] = hit->GetAdcDown();
    }
  }

  HF1(1, 22.);

  return true;
}

//_____________________________________________________________________________
Bool_t
UserKuramaTracking::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserKuramaTracking;
}

//_____________________________________________________________________________
Bool_t
ConfMan:: InitializeHistograms()
{
  const Int_t    NBinDTSDC1 = 180;
  const Double_t MinDTSDC1  = -30.;
  const Double_t MaxDTSDC1  = 120.;
  const Int_t    NBinDLSDC1 =  80;
  const Double_t MinDLSDC1  = -0.5;
  const Double_t MaxDLSDC1  =  3.5;

  const Int_t    NBinDTSDC2 = 360;
  const Double_t MinDTSDC2  = -50.;
  const Double_t MaxDTSDC2  = 250.;
  const Int_t    NBinDLSDC2 = 130;
  const Double_t MinDLSDC2  = -0.5;
  const Double_t MaxDLSDC2  =  6.;

  const Int_t    NBinDTSDC3 = 360;
  const Double_t MinDTSDC3  = -50.;
  const Double_t MaxDTSDC3  = 250.;
  const Int_t    NBinDLSDC3 = 150;
  const Double_t MinDLSDC3  = -2.0;
  const Double_t MaxDLSDC3  =  5.5;

  const Int_t    NBinDTSDC4 = 720;
  const Double_t MinDTSDC4  = -50.;
  const Double_t MaxDTSDC4  = 550.;
  const Int_t    NBinDLSDC4 = 300;
  const Double_t MinDLSDC4  = -3.0;
  const Double_t MaxDLSDC4  = 12.0;

  HB1(1, "Status", 30, 0., 30.);

  HB1(10, "#Tracks SdcIn", 10, 0., 10.);
  HB1(11, "#Hits of Track SdcIn", 15, 0., 15.);
  HB1(12, "Chisqr SdcIn", 500, 0., 50.);
  HB1(13, "LayerId SdcIn", 15, 0., 15.);
  HB1(14, "X0 SdcIn", 400, -100., 100.);
  HB1(15, "Y0 SdcIn", 400, -100., 100.);
  HB1(16, "U0 SdcIn", 200, -0.20, 0.20);
  HB1(17, "V0 SdcIn", 200, -0.20, 0.20);
  HB2(18, "U0%X0 SdcIn", 100, -100., 100., 100, -0.20, 0.20);
  HB2(19, "V0%Y0 SdcIn", 100, -100., 100., 100, -0.20, 0.20);
  HB2(20, "X0%Y0 SdcIn", 100, -100., 100., 100, -100., 100.);

  HB1(21, "#Hits of Track SdcIn [KuramaTrack]", 15, 0., 15.);
  HB1(22, "Chisqr SdcIn [KuramaTrack]", 500, 0., 50.);
  HB1(23, "LayerId SdcIn [KuramaTrack]", 15, 0., 15.);
  HB1(24, "X0 SdcIn [KuramaTrack]", 400, -100., 100.);
  HB1(25, "Y0 SdcIn [KuramaTrack]", 400, -100., 100.);
  HB1(26, "U0 SdcIn [KuramaTrack]", 200, -0.20, 0.20);
  HB1(27, "V0 SdcIn [KuramaTrack]", 200, -0.20, 0.20);
  HB2(28, "U0%X0 SdcIn [KuramaTrack]", 100, -100., 100., 100, -0.20, 0.20);
  HB2(29, "V0%Y0 SdcIn [KuramaTrack]", 100, -100., 100., 100, -0.20, 0.20);
  //HB2(30, "X0%Y0 SdcIn [KuramaTrack]", 100, -100., 100., 100, -100., 100.);

  HB1(30, "#Tracks SdcOut", 10, 0., 10.);
  HB1(31, "#Hits of Track SdcOut", 20, 0., 20.);
  HB1(32, "Chisqr SdcOut", 500, 0., 50.);
  HB1(33, "LayerId SdcOut", 20, 30., 50.);
  HB1(34, "X0 SdcOut", 1400, -1200., 1200.);
  HB1(35, "Y0 SdcOut", 1000, -500., 500.);
  HB1(36, "U0 SdcOut",  700, -0.35, 0.35);
  HB1(37, "V0 SdcOut",  200, -0.20, 0.20);
  HB2(38, "U0%X0 SdcOut", 120, -1200., 1200., 100, -0.40, 0.40);
  HB2(39, "V0%Y0 SdcOut", 100, -500., 500., 100, -0.20, 0.20);
  HB2(40, "X0%Y0 SdcOut", 100, -1200., 1200., 100, -500., 500.);

  HB1(41, "#Hits of Track SdcOut [KuramaTrack]", 20, 0., 20.);
  HB1(42, "Chisqr SdcOut [KuramaTrack]", 500, 0., 50.);
  HB1(43, "LayerId SdcOut [KuramaTrack]", 20, 30., 50.);
  HB1(44, "X0 SdcOut [KuramaTrack]", 1400, -1200., 1200.);
  HB1(45, "Y0 SdcOut [KuramaTrack]", 1000, -500., 500.);
  HB1(46, "U0 SdcOut [KuramaTrack]",  700, -0.35, 0.35);
  HB1(47, "V0 SdcOut [KuramaTrack]",  200, -0.10, 0.10);
  HB2(48, "U0%X0 SdcOut [KuramaTrack]", 120, -600., 600., 100, -0.40, 0.40);
  HB2(49, "V0%Y0 SdcOut [KuramaTrack]", 100, -500., 500., 100, -0.10, 0.10);
  //HB2(50, "X0%Y0 SdcOut [KuramaTrack]", 100, -700., 700., 100, -500., 500.);

  HB1(50, "#Tracks KURAMA", 10, 0., 10.);
  HB1(51, "#Hits of KuramaTrack", 50, 0., 50.);
  HB1(52, "Chisqr KuramaTrack", 500, 0., 50.);
  HB1(53, "LayerId KuramaTrack", 90, 0., 90.);
  HB1(54, "Xtgt KuramaTrack", 200, -100., 100.);
  HB1(55, "Ytgt KuramaTrack", 200, -100., 100.);
  HB1(56, "Utgt KuramaTrack", 300, -0.30, 0.30);
  HB1(57, "Vtgt KuramaTrack", 300, -0.20, 0.20);
  HB2(58, "U%Xtgt KuramaTrack", 100, -100., 100., 100, -0.25, 0.25);
  HB2(59, "V%Ytgt KuramaTrack", 100, -100., 100., 100, -0.10, 0.10);
  HB2(60, "Y%Xtgt KuramaTrack", 100, -100., 100., 100, -100., 100.);
  HB1(61, "P KuramaTrack", 500, 0.00, 2.50);
  HB1(62, "PathLength KuramaTrack", 600, 3000., 4000.);
  HB1(63, "MassSqr", 600, -0.4, 1.4);

  // SDC1
  for(Int_t i=1; i<=NumOfLayersSDC1; ++i){
    TString title1 = Form("HitPat Sdc1_%d", i);
    TString title2 = Form("DriftTime Sdc1_%d", i);
    TString title3 = Form("DriftLength Sdc1_%d", i);
    TString title4 = Form("Position Sdc1_%d", i);
    TString title5 = Form("Residual Sdc1_%d", i);
    TString title6 = Form("Resid%%Pos Sdc1_%d", i);
    TString title7 = Form("Y%%Xcal Sdc1_%d", i);
    HB1(100*i+1, title1, 96, 0., 96.);
    HB1(100*i+2, title2, NBinDTSDC1, MinDTSDC1, MaxDTSDC1);
    HB1(100*i+3, title3, NBinDLSDC1, MinDLSDC1, MaxDLSDC1);
    HB1(100*i+4, title4, 500, -100., 100.);
    HB1(100*i+5, title5, 200, -2.0, 2.0);
    HB2(100*i+6, title6, 100, -100., 100., 100, -2.0, 2.0);
    HB2(100*i+7, title7, 100, -100., 100., 100, -50., 50.);
    title1 += " [KuramaTrack]";
    title2 += " [KuramaTrack]";
    title3 += " [KuramaTrack]";
    title4 += " [KuramaTrack]";
    title5 += " [KuramaTrack]";
    title6 += " [KuramaTrack]";
    title7 += " [KuramaTrack]";
    TString title22 = Form("DriftLength%%DriftTime Sdc%d [KuramaTrack]", i);
    HB1(100*i+11, title1, 96, 0., 96.);
    HB1(100*i+12, title2, NBinDTSDC1, MinDTSDC1, MaxDTSDC1);
    HB1(100*i+13, title3, NBinDLSDC1, MinDLSDC1, MaxDLSDC1);
    HB1(100*i+14, title4, 500, -100., 100.);
    HB1(100*i+15, title5, 200, -2.0, 2.0);
    HB2(100*i+16, title6, 100, -100., 100., 100, -2.0, 2.0);
    HB2(100*i+17, title7, 100, -100., 100., 100, -50., 50.);
    HB2(100*i+22, title22, NBinDTSDC1, MinDTSDC1, MaxDTSDC1, NBinDLSDC1, MinDLSDC1, MaxDLSDC1);
  }

  //SDC2
  for(Int_t i=NumOfLayersSDC1+1; i<=NumOfLayersSdcIn; ++i){
    TString title1 = Form("HitPat Sdc2%d", i-NumOfLayersSDC1);
    TString title2 = Form("DriftTime Sdc2_%d", i-NumOfLayersSDC1);
    TString title3 = Form("DriftLength Sdc2_%d", i-NumOfLayersSDC1);
    TString title4 = Form("Position Sdc2%d", i-NumOfLayersSDC1);
    TString title5 = Form("Residual Sdc2%d", i-NumOfLayersSDC1);
    TString title6 = Form("Resid%%Pos Sdc2%d", i-NumOfLayersSDC1);
    TString title7 = Form("Y%%Xcal Sdc2%d", i-NumOfLayersSDC1);
    HB1(100*i+1, title1, 70, 0., 70.);
    HB1(100*i+2, title2, NBinDTSDC2, MinDTSDC2, MaxDTSDC2);
    HB1(100*i+3, title3, NBinDLSDC2, MinDLSDC2, MaxDLSDC2);
    HB1(100*i+4, title4, 800, -400., 400.);
    HB1(100*i+5, title5, 500, -5.0, 5.0);
    HB2(100*i+6, title6, 100, -600., 600., 100, -5.0, 5.0);
    HB2(100*i+7, title7, 100, -600., 600., 100, -600., 600.);
    title1 += " [KuramaTrack]";
    title4 += " [KuramaTrack]";
    title5 += " [KuramaTrack]";
    title6 += " [KuramaTrack]";
    title7 += " [KuramaTrack]";
    HB1(100*i+11, title1, 70, 0., 70.);
    HB1(100*i+14, title4, 800, -400., 400.);
    HB1(100*i+15, title5, 500, -5.0, 5.0);
    HB2(100*i+16, title6, 100, -600., 600., 100, -5.0, 5.0);
    HB2(100*i+17, title7, 100, -600., 600., 100, -600., 600.);
  }

  // SDC3
  for(Int_t i=NumOfLayersSdcIn+1; i<=(NumOfLayersSdcIn+NumOfLayersSDC3); ++i){
    TString title1 = Form("HitPat Sdc3_%d", i-NumOfLayersSdcIn);
    TString title2 = Form("DriftTime Sdc3_%d", i-NumOfLayersSdcIn);
    TString title3 = Form("DriftLength Sdc3_%d", i-NumOfLayersSdcIn);
    TString title4 = Form("Position Sdc3_%d", i-NumOfLayersSdcIn);
    TString title5 = Form("Residual Sdc3_%d", i-NumOfLayersSdcIn);
    TString title6 = Form("Resid%%Pos Sdc3_%d", i-NumOfLayersSdcIn);
    TString title7 = Form("Y%%Xcal Sdc3_%d", i-NumOfLayersSdcIn);
    HB1(100*i+1, title1, 112, 0., 112.);
    HB1(100*i+2, title2, NBinDTSDC3, MinDTSDC3, MaxDTSDC3);
    HB1(100*i+3, title3, NBinDLSDC3, MinDLSDC3, MaxDLSDC3);
    HB1(100*i+4, title4, 1000, -600., 600.);
    HB1(100*i+5, title5, 200, -2.0, 2.0);
    HB2(100*i+6, title6, 100, -600., 600., 100, -2.0, 2.0);
    HB2(100*i+7, title7, 100, -600., 600., 100, -600., 600.);
    title1 += " [KuramaTrack]";
    title2 += " [KuramaTrack]";
    title3 += " [KuramaTrack]";
    title4 += " [KuramaTrack]";
    title5 += " [KuramaTrack]";
    title6 += " [KuramaTrack]";
    title7 += " [KuramaTrack]";
    TString title22 = Form("DriftLength%%DriftTime Sdc3_%d [KuramaTrack]", i-NumOfLayersSdcIn);
    HB1(100*i+11, title1, 112, 0., 112.);
    HB1(100*i+12, title2, NBinDTSDC3, MinDTSDC3, MaxDTSDC3);
    HB1(100*i+13, title3, NBinDLSDC3, MinDLSDC3, MaxDLSDC3);
    HB1(100*i+14, title4, 1000, -600., 600.);
    HB1(100*i+15, title5, 200, -2.0, 2.0);
    HB2(100*i+16, title6, 100, -600., 600., 100, -2.0, 2.0);
    HB2(100*i+17, title7, 100, -600., 600., 100, -600., 600.);
    HB2(100*i+22, title22, NBinDTSDC3, MinDTSDC3, MaxDTSDC3, NBinDLSDC3, MinDLSDC3, MaxDLSDC3);
  }

  // SDC4
  for(Int_t i=NumOfLayersSdcIn+NumOfLayersSDC3+1;
      i<=(NumOfLayersSdcIn+NumOfLayersSDC3+NumOfLayersSDC4) ; ++i){
    TString title1 = Form("HitPat Sdc4_%d", i-(NumOfLayersSdcIn+NumOfLayersSDC3));
    TString title2 = Form("DriftTime Sdc4_%d", i-(NumOfLayersSdcIn+NumOfLayersSDC3));
    TString title3 = Form("DriftLength Sdc4_%d", i-(NumOfLayersSdcIn+NumOfLayersSDC3));
    TString title4 = Form("Position Sdc4_%d", i-(NumOfLayersSdcIn+NumOfLayersSDC3));
    TString title5 = Form("Residual Sdc4_%d", i-(NumOfLayersSdcIn+NumOfLayersSDC3));
    TString title6 = Form("Resid%%Pos Sdc4_%d", i-(NumOfLayersSdcIn+NumOfLayersSDC3));
    TString title7 = Form("Y%%Xcal Sdc4_%d", i-(NumOfLayersSdcIn+NumOfLayersSDC3));
    HB1(100*i+1, title1, 120, 0., 120.);
    HB1(100*i+2, title2, NBinDTSDC4, MinDTSDC4, MaxDTSDC4);
    HB1(100*i+3, title3, NBinDLSDC4, MinDLSDC4, MaxDLSDC4);
    HB1(100*i+4, title4, 1000, -600., 600.);
    HB1(100*i+5, title5, 200, -2.0, 2.0);
    HB2(100*i+6, title6, 100, -600., 600., 100, -1.0, 1.0);
    HB2(100*i+7, title7, 100, -600., 600., 100, -600., 600.);
    title1 += " [KuramaTrack]";
    title2 += " [KuramaTrack]";
    title3 += " [KuramaTrack]";
    title4 += " [KuramaTrack]";
    title5 += " [KuramaTrack]";
    title6 += " [KuramaTrack]";
    title7 += " [KuramaTrack]";
    TString title22 = Form("DriftLength%%DriftTime Sdc4_%d [KuramaTrack]", i-(NumOfLayersSdcIn+NumOfLayersSDC3));
    HB1(100*i+11, title1, 120, 0., 120.);
    HB1(100*i+12, title2, NBinDTSDC4, MinDTSDC4, MaxDTSDC4);
    HB1(100*i+13, title3, NBinDLSDC4, MinDLSDC4, MaxDLSDC4);
    HB1(100*i+14, title4, 1000, -600., 600.);
    HB1(100*i+15, title5, 200, -2.0, 2.0);
    HB2(100*i+16, title6, 100, -600., 600., 100, -2.0, 2.0);
    HB2(100*i+17, title7, 100, -600., 600., 100, -600., 600.);
    HB2(100*i+22, title22, NBinDTSDC4, MinDTSDC4, MaxDTSDC4, NBinDLSDC4, MinDLSDC4, MaxDLSDC4);
  }

  // TOF in SdcOut/KuramaTracking
  for(Int_t i=NumOfLayersSdcIn+NumOfLayersSdcOut+1;
      i<=NumOfLayersSdcIn+NumOfLayersSdcOut+4; ++i){
    TString title1 = Form("HitPat Tof%d", i-(NumOfLayersSdcIn+NumOfLayersSdcOut));
    TString title4 = Form("Position Tof%d", i-(NumOfLayersSdcIn+NumOfLayersSdcOut));
    TString title5 = Form("Residual Tof%d", i-(NumOfLayersSdcIn+NumOfLayersSdcOut));
    TString title6 = Form("Resid%%Pos Tof%d", i-(NumOfLayersSdcIn+NumOfLayersSdcOut));
    TString title7 = Form("Y%%Xcal Tof%d", i-(NumOfLayersSdcIn+NumOfLayersSdcOut));
    HB1(100*i+1, title1, 200, 0., 200.);
    HB1(100*i+4, title4, 1000, -1000., 1000.);
    HB1(100*i+5, title5, 200, -20., 20.);
    HB2(100*i+6, title6, 100, -1000., 1000., 100, -200., 200.);
    HB2(100*i+7, title7, 100, -1000., 1000., 100, -1000., 1000.);
    title1 += " [KuramaTrack]";
    title4 += " [KuramaTrack]";
    title5 += " [KuramaTrack]";
    title6 += " [KuramaTrack]";
    title7 += " [KuramaTrack]";
    HB1(100*i+11, title1, 200, 0., 200.);
    HB1(100*i+14, title4, 1000, -1000., 1000.);
    HB1(100*i+15, title5, 200, -20., 20.);
    HB2(100*i+16, title6, 100, -1000., 1000., 100, -200., 200.);
    HB2(100*i+17, title7, 100, -1000., 1000., 100, -1000., 1000.);
  }

  HB2(20001, "Xout%Xin", 100, -200., 200., 100, -200., 200.);
  HB2(20002, "Yout%Yin", 100, -200., 200., 100, -200., 200.);
  HB2(20003, "Uout%Uin", 100, -0.5,  0.5,  100, -0.5,  0.5);
  HB2(20004, "Vin%Vout", 100, -0.1,  0.1,  100, -0.1,  0.1);

  HB2(20021, "Yout%Yin [KuramaTrack]", 100, -150., 150., 120, -300., 300.);
  HB2(20022, "Vout%Vin [KuramaTrack]", 100, -0.05, 0.05, 100, -0.1, 0.1);
  HB2(20023, "Yout%Vin [KuramaTrack]", 100, -0.05, 0.05, 100, -300., 300.);
  HB2(20024, "Yin%Vout [KuramaTrack]", 100, -0.10, 0.10, 100, -150., 150.);

  ////////////////////////////////////////////
  //Tree
  HBTree("kurama","tree of KuramaTracking");
  tree->Branch("evnum",     &event.evnum,    "evnum/I");
  tree->Branch("trigpat",    event.trigpat,  Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag, Form("trigflag[%d]/I", NumOfSegTrig));

  //Hodoscope
  tree->Branch("nhBh2",   &event.nhBh2,   "nhBh2/I");
  tree->Branch("Bh2Seg",   event.Bh2Seg,  "Bh2Seg[nhBh2]/D");
  tree->Branch("tBh2",     event.tBh2,    "tBh2[nhBh2]/D");
  tree->Branch("deBh2",    event.deBh2,   "deBh2[nhBh2]/D");
  tree->Branch("time0",   &event.time0,   "time0/D");

  tree->Branch("nhBh1",   &event.nhBh1,   "nhBh1/I");
  tree->Branch("Bh1Seg",   event.Bh1Seg,  "Bh1Seg[nhBh1]/D");
  tree->Branch("tBh1",     event.tBh1,    "tBh1[nhBh1]/D");
  tree->Branch("deBh1",    event.deBh1,   "deBh1[nhBh1]/D");
  tree->Branch("btof",    &event.btof,    "btof/D");

  tree->Branch("nhTof",   &event.nhTof,   "nhTof/I");
  tree->Branch("TofSeg",   event.TofSeg,  "TofSeg[nhTof]/D");
  tree->Branch("tTof",     event.tTof,    "tTof[nhTof]/D");
  tree->Branch("dtTof",    event.dtTof,   "dtTof[nhTof]/D");
  tree->Branch("deTof",    event.deTof,   "deTof[nhTof]/D");

  tree->Branch("wposSdcIn",  event.wposSdcIn,  Form("wposSdcIn[%d]/D", NumOfLayersSdcIn));
  tree->Branch("wposSdcOut", event.wposSdcOut, Form("wposSdcOut[%d]/D", NumOfLayersSdcOut));

  //Tracking
  tree->Branch("ntSdcIn",    &event.ntSdcIn,     "ntSdcIn/I");
  tree->Branch("much",       &event.much,        "much/I");
  tree->Branch("nlSdcIn",    &event.nlSdcIn,     "nlSdcIn/I");
  tree->Branch("nhSdcIn",     event.nhSdcIn,     "nhSdcIn[ntSdcIn]/I");
  tree->Branch("chisqrSdcIn", event.chisqrSdcIn, "chisqrSdcIn[ntSdcIn]/D");
  tree->Branch("x0SdcIn",     event.x0SdcIn,     "x0SdcIn[ntSdcIn]/D");
  tree->Branch("y0SdcIn",     event.y0SdcIn,     "y0SdcIn[ntSdcIn]/D");
  tree->Branch("u0SdcIn",     event.u0SdcIn,     "u0SdcIn[ntSdcIn]/D");
  tree->Branch("v0SdcIn",     event.v0SdcIn,     "v0SdcIn[ntSdcIn]/D");

  tree->Branch("ntSdcOut",   &event.ntSdcOut,     "ntSdcOut/I");
  tree->Branch("nlSdcOut",   &event.nlSdcOut,     "nlSdcOut/I");
  tree->Branch("nhSdcOut",    event.nhSdcOut,     "nhSdcOut[ntSdcOut]/I");
  tree->Branch("chisqrSdcOut",event.chisqrSdcOut, "chisqrSdcOut[ntSdcOut]/D");
  tree->Branch("x0SdcOut",    event.x0SdcOut,     "x0SdcOut[ntSdcOut]/D");
  tree->Branch("y0SdcOut",    event.y0SdcOut,     "y0SdcOut[ntSdcOut]/D");
  tree->Branch("u0SdcOut",    event.u0SdcOut,     "u0SdcOut[ntSdcOut]/D");
  tree->Branch("v0SdcOut",    event.v0SdcOut,     "v0SdcOut[ntSdcOut]/D");

  // KURAMA Tracking
  tree->Branch("ntKurama",    &event.ntKurama,     "ntKurama/I");
  tree->Branch("nlKurama",    &event.nlKurama,     "nlKurama/I");
  tree->Branch("nhKurama",     event.nhKurama,     "nhKurama[ntKurama]/I");
  tree->Branch("chisqrKurama", event.chisqrKurama, "chisqrKurama[ntKurama]/D");
  tree->Branch("stof",         event.stof,         "stof[ntKurama]/D");
  tree->Branch("path",         event.path,         "path[ntKurama]/D");
  tree->Branch("pKurama",      event.pKurama,      "pKurama[ntKurama]/D");
  tree->Branch("qKurama",      event.qKurama,      "qKurama[ntKurama]/D");
  tree->Branch("m2",           event.m2,           "m2[ntKurama]/D");

  tree->Branch("xtgtKurama",   event.xtgtKurama,   "xtgtKurama[ntKurama]/D");
  tree->Branch("ytgtKurama",   event.ytgtKurama,   "ytgtKurama[ntKurama]/D");
  tree->Branch("utgtKurama",   event.utgtKurama,   "utgtKurama[ntKurama]/D");
  tree->Branch("vtgtKurama",   event.vtgtKurama,   "vtgtKurama[ntKurama]/D");

  tree->Branch("thetaKurama",  event.thetaKurama,  "thetaKurama[ntKurama]/D");
  tree->Branch("phiKurama",    event.phiKurama,    "phiKurama[ntKurama]/D");
  tree->Branch("resP",         event.resP,         "resP[ntKurama]/D");

  tree->Branch("xtofKurama",   event.xtofKurama,   "xtofKurama[ntKurama]/D");
  tree->Branch("ytofKurama",   event.ytofKurama,   "ytofKurama[ntKurama]/D");
  tree->Branch("utofKurama",   event.utofKurama,   "utofKurama[ntKurama]/D");
  tree->Branch("vtofKurama",   event.vtofKurama,   "vtofKurama[ntKurama]/D");
  tree->Branch("tofsegKurama", event.tofsegKurama, "tofsegKurama[ntKurama]/D");

  tree->Branch("xwcKurama",    event.xwcKurama,    "xwcKurama[ntKurama]/D");
  tree->Branch("ywcKurama",    event.ywcKurama,    "ywcKurama[ntKurama]/D");
  tree->Branch("pathwcKurama", event.pathwcKurama, "pathwcKurama[ntKurama]/D");

  tree->Branch("vpx",          event.vpx,          Form("vpx[%d]/D", NumOfLayersVP));
  tree->Branch("vpy",          event.vpy,          Form("vpy[%d]/D", NumOfLayersVP));
  tree->Branch("vpxtpc",       event.vpxtpc,       Form("vpxtpc[ntKurama][%d]/D", NumOfLayersVP));
  tree->Branch("vpytpc",       event.vpytpc,       Form("vpytpc[ntKurama][%d]/D", NumOfLayersVP));
  tree->Branch("vpztpc",       event.vpztpc,       Form("vpztpc[ntKurama][%d]/D", NumOfLayersVP));
  tree->Branch("vpxhtof",      event.vpxhtof,      "vpxhtof[ntKurama]/D");
  tree->Branch("vpyhtof",      event.vpyhtof,      "vpyhtof[ntKurama]/D");
  tree->Branch("vpzhtof",      event.vpzhtof,      "vpzhtof[ntKurama]/D");

  event.resL.resize(NumOfLayersSdcIn+NumOfLayersSdcOut+4);
  event.resG.resize(NumOfLayersSdcIn+NumOfLayersSdcOut+4);
  for(Int_t i=0; i<NumOfLayersSdcIn; ++i){
    tree->Branch(Form("ResL%d",i+1), &event.resL[i]);
  }

  for(Int_t i=0; i<NumOfLayersSdcOut; ++i){
    tree->Branch(Form("ResL%d",i+31), &event.resL[i+NumOfLayersSdcIn]);
  }

  tree->Branch("ResL41", &event.resL[NumOfLayersSdcIn+NumOfLayersSdcOut]);
  tree->Branch("ResL42", &event.resL[NumOfLayersSdcIn+NumOfLayersSdcOut+1]);
  tree->Branch("ResL43", &event.resL[NumOfLayersSdcIn+NumOfLayersSdcOut+2]);
  tree->Branch("ResL44", &event.resL[NumOfLayersSdcIn+NumOfLayersSdcOut+3]);

  for(Int_t i=0; i<NumOfLayersSdcIn; ++i){
    tree->Branch(Form("ResG%d",i+ 1), &event.resG[i]);
  }

  for(Int_t i=0; i<NumOfLayersSdcOut; ++i){
    tree->Branch(Form("ResG%d",i+31), &event.resG[i+NumOfLayersSdcIn]);
  }

  tree->Branch("ResG41", &event.resG[NumOfLayersSdcIn+NumOfLayersSdcOut]);
  tree->Branch("ResG42", &event.resG[NumOfLayersSdcIn+NumOfLayersSdcOut+1]);
  tree->Branch("ResG43", &event.resG[NumOfLayersSdcIn+NumOfLayersSdcOut+2]);
  tree->Branch("ResG44", &event.resG[NumOfLayersSdcIn+NumOfLayersSdcOut+3]);

  tree->Branch("tTofCalc",  event.tTofCalc,  "tTofCalc[3]/D");
  tree->Branch("utTofSeg",  event.utTofSeg,  Form("utTofSeg[%d]/D", NumOfSegTOF));
  tree->Branch("dtTofSeg",  event.dtTofSeg,  Form("dtTofSeg[%d]/D", NumOfSegTOF));
  tree->Branch("udeTofSeg", event.udeTofSeg, Form("udeTofSeg[%d]/D", NumOfSegTOF));
  tree->Branch("ddeTofSeg", event.ddeTofSeg, Form("ddeTofSeg[%d]/D", NumOfSegTOF));
  tree->Branch("tofua",     event.tofua,     Form("tofua[%d]/D", NumOfSegTOF));
  tree->Branch("tofda",     event.tofda,     Form("tofda[%d]/D", NumOfSegTOF));

  tree->Branch("nh", event.nh, "[ntKurama]/I");
  tree->Branch("xout", event.xout, "[ntKurama]/D");
  tree->Branch("yout", event.yout, "[ntKurama]/D");
  tree->Branch("zout", event.zout, "[ntKurama]/D");
  tree->Branch("pxout", event.pxout, "[ntKurama]/D");
  tree->Branch("pyout", event.pyout, "[ntKurama]/D");
  tree->Branch("pzout", event.pzout, "[ntKurama]/D");
  tree->Branch("xin", event.xin, "[ntKurama]/D");
  tree->Branch("yin", event.yin, "[ntKurama]/D");
  tree->Branch("zin", event.zin, "[ntKurama]/D");
  tree->Branch("pxin", event.pxin, "[ntKurama]/D");
  tree->Branch("pyin", event.pyin, "[ntKurama]/D");
  tree->Branch("pzin", event.pzin, "[ntKurama]/D");
  tree->Branch("xtof", event.xtof, "[ntKurama]/D");
  tree->Branch("ytof", event.ytof, "[ntKurama]/D");
  tree->Branch("ztof", event.ztof, "[ntKurama]/D");
  tree->Branch("pxtof", event.pxtof, "[ntKurama]/D");
  tree->Branch("pytof", event.pytof, "[ntKurama]/D");
  tree->Branch("pztof", event.pztof, "[ntKurama]/D");
  tree->Branch("layer", event.layer, "[ntKurama][22]/I");
  tree->Branch("wire", event.wire, "[ntKurama][22]/D");
  tree->Branch("localhitpos", event.localhitpos, "[ntKurama][22]/D");
  tree->Branch("wpos", event.wpos, "[ntKurama][22]/D");

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
     InitializeParameter<FieldMan>("FLDMAP", "HSFLDMAP") &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

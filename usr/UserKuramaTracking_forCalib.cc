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

namespace
{
using namespace root;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& zTOF = gGeom.LocalZ("TOF");
const Bool_t exclusive_tracking = true;
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
  static TString s_name("UserKuramaTracking_forCalib");
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

  // Calib
  enum eParticle { Pion, Kaon, Proton, nParticle };
  Double_t tTofCalc[nParticle];
  Double_t utTofSeg[NumOfSegTOF];
  Double_t dtTofSeg[NumOfSegTOF];
  Double_t udeTofSeg[NumOfSegTOF];
  Double_t ddeTofSeg[NumOfSegTOF];
  Double_t tofua[NumOfSegTOF];
  Double_t tofda[NumOfSegTOF];

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

  static const auto MaxMultiHitSdcIn  = gUser.GetParameter("MaxMultiHitSdcIn");
  static const auto MaxMultiHitSdcOut = gUser.GetParameter("MaxMultiHitSdcOut");

  static const auto MinStopTimingSdcOut = gUser.GetParameter("StopTimingSdcOut", 0);
  static const auto MaxStopTimingSdcOut = gUser.GetParameter("StopTimingSdcOut", 1);
  // static const auto StopTimeDiffSdcOut = gUser.GetParameter("StopTimeDiffSdcOut");

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
      TOFCont.push_back(hit);
#if HodoCut
      Double_t stof = cmt-time0;
      if(MinDeTOF<de  && de<MaxDeTOF  &&
         MinTimeTOF<stof && stof<MaxTimeTOF){
	++nhOk;
      }
#endif
    }

    HF1(1, 4.);

#if HodoCut
    if(nhOk<1) return true;
#endif
  }

  // Common stop timing
  Bool_t common_stop_is_tof = (common_stop_tdc < MinStopTimingSdcOut
                               || MaxStopTimingSdcOut < common_stop_tdc);
  if(!common_stop_is_tof) return true;

  HF1(1, 5.);

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
	Int_t layerId = hit->GetLayer();
	Double_t wire = hit->GetWire();
	Double_t wpos = hit->GetWirePosition();
	event.wposSdcIn[layer-1] = wpos;

	Int_t nhtdc = hit->GetTdcSize();
	Int_t tdc1st = -1;
	for(Int_t k=0; k<nhtdc; k++){
	  Int_t tdc = hit->GetTdcVal(k);
	  if(tdc > tdc1st){
	    tdc1st = tdc;
	  }
	}
	HF1(10000*layerId+Int_t(wire), -0.833*tdc1st);
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
	Int_t layerId = hit->GetLayer();
	Double_t wire = hit->GetWire();
	Double_t wpos = hit->GetWirePosition();
	event.wposSdcOut[layer-1] = wpos;

	Int_t nhtdc = hit->GetTdcSize();
	Int_t tdc1st = -1;
	for(Int_t k=0; k<nhtdc; k++){
	  Int_t tdc = hit->GetTdcVal(k);
	  if(tdc > tdc1st){
	    tdc1st = tdc;
	  }
	}

	layerId -= 20;
	HF1(10000*layerId+Int_t(wire), -0.833*tdc1st);
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

  HF1(1, 6.);

  //std::cout << "==========TrackSearch SdcIn============" << std::endl;
  //DCAna->TrackSearchSdcIn();
  DCAna->TrackSearchSdcIn(exclusive_tracking);
  // DCAna->ChiSqrCutSdcIn(50.);
  //std::cout << "==========End of TrackSearch SdcIn============" << std::endl;
  Int_t ntSdcIn = DCAna->GetNtracksSdcIn();
  if(MaxHits<ntSdcIn){
    std::cout << "#W " << FUNC_NAME << " "
	      << "too many ntSdcIn " << ntSdcIn << "/" << MaxHits << std::endl;
    ntSdcIn = MaxHits;
  }
  if(ntSdcIn!=1) return true;

  event.ntSdcIn = ntSdcIn;
  Int_t much_combi = DCAna->MuchCombinationSdcIn();
  event.much = much_combi;
  for(Int_t it=0; it<ntSdcIn; ++it){
    DCLocalTrack *tp=DCAna->GetTrackSdcIn(it);
    Int_t nh=tp->GetNHit();
    Double_t chisqr=tp->GetChiSquare();
    Double_t x0=tp->GetX0(), y0=tp->GetY0();
    Double_t u0=tp->GetU0(), v0=tp->GetV0();
    event.nhSdcIn[it] = nh;
    event.chisqrSdcIn[it] = chisqr;
    event.x0SdcIn[it] = x0;
    event.y0SdcIn[it] = y0;
    event.u0SdcIn[it] = u0;
    event.v0SdcIn[it] = v0;
  }

  HF1(1, 7.);

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

  HF1(1, 8.);

  //////////////SdcOut tracking
  //std::cout << "==========TrackSearch SdcOut============" << std::endl;
#if UseTOF
  DCAna->TrackSearchSdcOut(TOFCont, exclusive_tracking);
#else
  DCAna->TrackSearchSdcOut(exclusive_tracking);
#endif
  //std::cout << "==========End of TrackSearch SdcOut============" << std::endl;

  // DCAna->ChiSqrCutSdcOut(50.);
  Int_t ntSdcOut = DCAna->GetNtracksSdcOut();
  if(ntSdcOut!=1) return true;
  if(MaxHits<ntSdcOut){
    std::cout << "#W " << FUNC_NAME << " "
	      << "too many ntSdcOut " << ntSdcOut << "/" << MaxHits << std::endl;
    ntSdcOut = MaxHits;
  }
  event.ntSdcOut=ntSdcOut;
  for(Int_t it=0; it<ntSdcOut; ++it){
    DCLocalTrack *tp=DCAna->GetTrackSdcOut(it);
    Int_t nh=tp->GetNHit();
    Double_t chisqr=tp->GetChiSquare();
    Double_t u0=tp->GetU0(), v0=tp->GetV0();
    Double_t x0=tp->GetX0(), y0=tp->GetY0();
    event.nhSdcOut[it] = nh;
    event.chisqrSdcOut[it] = chisqr;
    event.x0SdcOut[it] = x0;
    event.y0SdcOut[it] = y0;
    event.u0SdcOut[it] = u0;
    event.v0SdcOut[it] = v0;
    Double_t z_tof = qnan;
    for(Int_t itof=0; itof<event.nhTof; ++itof){
      Int_t lnum = 0;
      TVector3 gposTof;
      if((Int_t)event.TofSeg[itof]%2 == 0){
        lnum = gGeom.GetDetectorId("TOF-UX");
        gposTof = gGeom.GetGlobalPosition("TOF-UX");
        z_tof = gGeom.LocalZ("TOF-UX");
      }
      if((Int_t)event.TofSeg[itof]%2 == 1){
        lnum = gGeom.GetDetectorId("TOF-DX");
        gposTof = gGeom.GetGlobalPosition("TOF-DX");
	z_tof = gGeom.LocalZ("TOF-DX");
      }
      Double_t wpos = gGeom.CalcWirePosition(lnum, event.TofSeg[itof]);
      TVector3 w(wpos, 0, 0);
      Double_t ytTof = event.dtTof[itof]*77.7481;
      TVector3 posTof = gposTof + w + TVector3(0, ytTof, 0);
      Double_t xtof = tp->GetX(z_tof), ytof = tp->GetY(z_tof);
      if(event.nhTof == 1){
	HF2(40, event.TofSeg[itof], xtof);
	HF2(41, event.TofSeg[itof], posTof.X() - xtof);
	HF2(42, posTof.X(), xtof);
	HF1(43, posTof.X() - xtof);
	HF2(44, posTof.X(), posTof.X() - xtof);
        HF2(45, event.dtTof[itof], ytof);
        HF1(46, ytTof - ytof);
        HF2(47, event.dtTof[itof], ytTof - ytof);
	HF2(48, event.TofSeg[itof], ytTof - ytof);

	Int_t id = (Int_t)event.TofSeg[itof];
	HF1(100*id+50, posTof.X() - xtof);
	HF2(100*id+51, event.dtTof[itof], ytof);
      }
    }
  }

  HF1(1, 9.);

  ///// BTOF BH2-Target
  static const auto StofOffset = gUser.GetParameter("StofOffset");

  //////////////KURAMA Tracking
  DCAna->TrackSearchKurama();
  Int_t ntKurama = DCAna->GetNTracksKurama();
  if(MaxHits < ntKurama){
    std::cout << "#W " << FUNC_NAME << " "
	      << "too many ntKurama " << ntKurama << "/" << MaxHits << std::endl;
    ntKurama = MaxHits;
  }
  event.ntKurama = ntKurama;
  HF1(10, ntKurama);

  //std::cout << "==========Kurama============" << std::endl;
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
    HF1(11, nh);
    HF1(12, chisqr);
    HF1(14, p);
    HF1(15, path);
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
      for(Int_t l = 0; l<NumOfLayersVP; ++l){
	Double_t x, y;
	track->GetTrajectoryLocalPosition(21 + l, x, y);
	event.vpx[l] = x;
	event.vpy[l] = y;
      }// for(l)
    }
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
      HF1(16, m2);
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

    if(chisqr>200.||nh!=20) continue;
    for(Int_t j=0; j<nh; ++j){
      TrackHit *hit=track->GetHit(j);
      if(!hit) continue;
      Int_t layerId = hit->GetLayer();
      if(hit->GetLayer()>79) layerId -= 62;
      else if(hit->GetLayer()>40) layerId -= 22;
      else if(hit->GetLayer()>30) layerId -= 20;

      HF1(13, layerId);
      //Double_t pos  = hit->GetCalLPos();
      Double_t res  = hit->GetResidual();
      Double_t wire = hit->GetHit()->GetWire();
      Double_t dt   = hit->GetHit()->GetDriftTime();
      Double_t dl   = hit->GetHit()->GetDriftLength();
      DCLTrackHit *lhit = hit->GetHit();
      Double_t wp   = lhit->GetWirePosition();
      Double_t xlcal = lhit->GetLocalCalPos();
      Double_t xlres = lhit->GetResidual();
      Double_t xlexres = lhit->GetResidualExclusive();

      HF1(100*layerId+11, Double_t(wire)-0.5);

      if(layerId > NumOfLayersSdcIn+NumOfLayersSdcOut) continue;
      Int_t tdc = lhit->GetTdcVal();
      Int_t tot = lhit->GetTot();
      HF1(100*layerId+1, tdc);
      HF1(100*layerId+2, tot);
      if(theta>=0 && theta<15){
	HF1(100*layerId+8, res);
	HF1(100*layerId+9, xlres);
	HF1(100*layerId+10, xlexres);
      }
      HF1(100*layerId+12, dt);
      HF1(100*layerId+13, dl);
      HF1(100*layerId+14, xlres);
      HF1(100*layerId+15, res);

      HFProf(100*layerId+16, dt, std::abs(xlcal-wp));
      HF2(100*layerId+18, dt, std::abs(xlcal-wp));
      HF2(100*layerId+19, dt, xlcal-wp);
      if(std::abs(dl-std::abs(xlcal-wp))<3.0) {
	HFProf(100*layerId+20, dt, std::abs(xlcal-wp));
	HF2(100*layerId+22, dt, std::abs(xlcal-wp));
	HF2(100*layerId+23, dt , xlcal-wp);
      }
    }

    DCLocalTrack *tpin = track->GetLocalTrackIn();
    Double_t chisqrsdcin = tpin->GetChiSquare();
    HF1(2, chisqrsdcin);
    DCLocalTrack *tp = track->GetLocalTrackOut();
    Double_t chisqrsdcout = tp->GetChiSquare();
    HF1(3, chisqrsdcout);
    HF1(31,tp->GetX0() - tpin->GetX0());
    HF1(32,tp->GetY0() - tpin->GetY0());
    HF1(33,tp->GetU0() - tpin->GetU0());
    HF1(34,tp->GetV0() - tpin->GetV0());

    Double_t z_tof = qnan;
    for(Int_t itof=0; itof<event.nhTof; ++itof){
      Int_t lnum = 0;
      TVector3 gposTof;
      if((Int_t)event.TofSeg[itof]%2 == 0){
        lnum = gGeom.GetDetectorId("TOF-UX");
        gposTof = gGeom.GetGlobalPosition("TOF-UX");
        z_tof = gGeom.LocalZ("TOF-UX");
      }
      if((Int_t)event.TofSeg[itof]%2 == 1){
        lnum = gGeom.GetDetectorId("TOF-DX");
        gposTof = gGeom.GetGlobalPosition("TOF-DX");
	z_tof = gGeom.LocalZ("TOF-DX");
      }
      Double_t wpos = gGeom.CalcWirePosition(lnum, event.TofSeg[itof]);
      TVector3 w(wpos, 0, 0);
      Double_t ytTof = event.dtTof[itof]*77.7481;
      TVector3 posTof = gposTof + w + TVector3(0, ytTof, 0);
      Double_t xtof=tp->GetX(z_tof), ytof=tp->GetY(z_tof);
      if(event.nhTof == 1){
	HF2(140, event.TofSeg[itof], xtof);
	HF2(141, event.TofSeg[itof], posTof.X() - xtof);
	HF2(142, posTof.X(), xtof);
	HF1(143, posTof.X() - xtof);
	HF2(144, posTof.X(), posTof.X() - xtof);
        HF2(145, event.dtTof[itof], ytof);
        HF1(146, ytTof - ytof);
        HF2(147, event.dtTof[itof], ytTof - ytof);
	HF2(148, event.TofSeg[itof], ytTof - ytof);

	Int_t id = (Int_t)event.TofSeg[itof];
	HF1(100*id+52, posTof.X() - xtof);
	HF2(100*id+53, event.dtTof[itof], ytof);
      }
    }

    {
      for(Int_t j=0; j<nh; ++j){
	TrackHit *hit=track->GetHit(j);
	if(!hit) continue;
	Int_t layerId = hit->GetLayer();
	if(hit->GetLayer()>79) continue;
	else if(hit->GetLayer()>40) layerId -= 22;
	else if(hit->GetLayer()>30) layerId -= 20;

	Double_t wire = hit->GetHit()->GetWire();
	Double_t dt   = hit->GetHit()->GetDriftTime();
	Double_t dl   = hit->GetHit()->GetDriftLength();
	DCLTrackHit *lhit = hit->GetHit();
	Double_t z = lhit->GetZ();
	Double_t wp   = lhit->GetWirePosition();
	Double_t tilt = lhit->GetTiltAngle()*TMath::DegToRad();
	Double_t u0_out=tp->GetU0(), v0_out=tp->GetV0();
	Double_t localcalpos_out = tp -> GetX(z)*TMath::Cos(tilt) + tp -> GetY(z)*TMath::Sin(tilt);
	Double_t u0_in=tpin->GetU0(), v0_in=tpin->GetV0();
	Double_t localcalpos_in = tpin -> GetX(z)*TMath::Cos(tilt) + tpin -> GetY(z)*TMath::Sin(tilt);
	if(layerId<11){
	  Double_t dsdz = u0_out*TMath::Cos(tilt)+v0_out*TMath::Sin(tilt);
	  Double_t coss = TMath::Cos(TMath::ATan(dsdz));
	  Double_t scal_out = localcalpos_out;
	  Double_t ss   = wp+dl/coss;
	  if(scal_out<wp) ss = wp-dl/coss;
	  //Double_t resi_out = (ss-scal_out)*coss;
	  Double_t resi_out = ss-scal_out;
	  HF1(100*layerId+54, resi_out);

	  HFProf(100*layerId+56, dt, std::abs(localcalpos_out-wp));
	  HF2(100*layerId+58, dt, std::abs(localcalpos_out-wp));
	  HF2(100*layerId+59, dt, localcalpos_out-wp);
	  if(std::abs(dl-std::abs(localcalpos_out-wp))<3.0) {
	    HFProf(100*layerId+60, dt, std::abs(localcalpos_out-wp));
	    HF2(100*layerId+62, dt, std::abs(localcalpos_out-wp));
	    HF2(100*layerId+63, dt , localcalpos_out-wp);
	  }
	  HF1(100*layerId+64, wp - localcalpos_in);
	  HF1(100*layerId+65, wp - localcalpos_out);
	  HF2(100*layerId+66, wp, localcalpos_in);
	  HF2(100*layerId+67, wp, localcalpos_out);
	  HF2(100*layerId+68, wp, wp - localcalpos_in);
	  HF2(100*layerId+69, wp, wp - localcalpos_out);
	}
	else{
	  Double_t dsdz = u0_in*TMath::Cos(tilt)+v0_in*TMath::Sin(tilt);
	  Double_t coss = TMath::Cos(TMath::ATan(dsdz));
	  Double_t scal_in = localcalpos_in;
	  Double_t ss   = wp+dl/coss;
	  if(scal_in<wp) ss = wp-dl/coss;
	  //Double_t resi_in = (ss-scal_in)*coss;
	  Double_t resi_in = ss-scal_in;

	  HF1(100*layerId+54, resi_in);
	  HFProf(100*layerId+56, dt, std::abs(localcalpos_in-wp));
	  HF2(100*layerId+58, dt, std::abs(localcalpos_in-wp));
	  HF2(100*layerId+59, dt, localcalpos_in-wp);
	  if(std::abs(dl-std::abs(localcalpos_in-wp))<3.0) {
	    HFProf(100*layerId+60, dt, std::abs(localcalpos_in-wp));
	    HF2(100*layerId+62, dt, std::abs(localcalpos_in-wp));
	    HF2(100*layerId+63, dt , localcalpos_in-wp);
	  }
	  HF1(100*layerId+64, wp - localcalpos_in);
	  HF1(100*layerId+65, wp - localcalpos_out);
	  HF2(100*layerId+66, wp, localcalpos_in);
	  HF2(100*layerId+67, wp, localcalpos_out);
	  HF2(100*layerId+68, wp, wp - localcalpos_in);
	  HF2(100*layerId+69, wp, wp - localcalpos_out);
	}
      }
    }
  }

  //std::cout << "==========TOF============" << std::endl;
  // TOF
  {
    Int_t nh = hodoAna->GetNHitsTOF();
    for(Int_t i=0; i<nh; ++i){
      Hodo2Hit *hit=hodoAna->GetHitTOF(i);
      if(!hit) continue;
      Int_t seg=hit->SegmentId()+1;
      Double_t tu = hit->GetTUp(), td=hit->GetTDown();
      // Double_t ctu=hit->GetCTUp(), ctd=hit->GetCTDown();
      //Double_t cmt=hit->CMeanTime();//, t= cmt-time0+StofOffset;//cmt-time0;
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

  HF1(1, 10.);

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
ConfMan::InitializeHistograms()
{

  const Int_t    NbinSdcTdc = 2000;
  const Double_t MinSdcTdc  = 0.;
  const Double_t MaxSdcTdc  = 2000.;

  const Int_t    NbinSdcTime = 2000*0.833;
  const Double_t MinSdcTime  = -0.833*2000.;
  const Double_t MaxSdcTime  = 0.;

  const Int_t    NBinDTSDC1 = 180;
  const Double_t MinDTSDC1  = -30.;
  const Double_t MaxDTSDC1  = 120.;
  const Int_t    NBinDLSDC1 =  80;
  const Double_t MinDLSDC1  = -0.5;
  const Double_t MaxDLSDC1  =  3.5;

  const Int_t    NBinDTSDC2 = 560;
  const Double_t MinDTSDC2  = -50.;
  const Double_t MaxDTSDC2  = 250.;
  const Int_t    NBinDLSDC2 = 130;
  const Double_t MinDLSDC2  = -0.5;
  const Double_t MaxDLSDC2  =  6.;

  const Int_t    NBinDTSDC3 = 560;
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
  HB1(2, "Chisqr SdcInTrack", 2500, 0., 500.);
  HB1(3, "Chisqr SdcOutTrack", 5000, 0., 1000.);
  HB1(10, "#Tracks KURAMA", 10, 0., 10.);
  HB1(11, "#Hits of KuramaTrack", 30, 0., 30.);
  HB1(12, "Chisqr KuramaTrack", 2500, 0., 500.);
  HB1(13, "LayerId KuramaTrack", 30, 0., 30.);
  HB1(14, "P KuramaTrack", 500, 0.00, 2.50);
  HB1(15, "PathLength KuramaTrack", 600, 3000., 4000.);
  HB1(16, "MassSqr", 600, -0.4, 1.4);

  HB2(20, "xSdcOut % TofSeg", NumOfSegTOF, 1, NumOfSegTOF + 1, 400, -1000, 1000.);
  HB1(31, "xSdcOut - xSdcIn", 300, -15, 15);
  HB1(32, "ySdcOut - ySdcIn", 300, -15, 15);
  HB1(33, "uSdcOut - uSdcIn", 300, -0.3, 0.3);
  HB1(34, "vSdcOut - vSdcIn", 300, -0.3, 0.3);

  // ySdcOut vs dtTof
  HB2(40, "xSdcOut % TofSeg", NumOfSegTOF, 1, NumOfSegTOF + 1, 400, -1000, 1000.);
  HB2(41, "TofPos - xSdcOut % TofSeg", NumOfSegTOF, 1, NumOfSegTOF + 1, 400, -200, 200);
  HB2(42, "xSdcOut % TofPos", 40, -1500, 1500, 400, -1000, 1000.);
  HB1(43, "TofPos - xSdcOut", 400, -200, 200);
  HB2(44, "TofPos - xSdcOut % TofPos", 40, -1500, 1500., 400, -200, 200);
  HB2(45, "ySdcOut % dtTof", 300, -15., 15., 400, -1000, 1000.);
  HB1(46, "ytTof - ySdcOut", 300, -300., 300.);
  HB2(47, "ytTof - ySdcOut % dtTof", 300, -15., 15., 300, -300, 300.);
  HB2(48, "ytTof - ySdcOut % TofSeg", NumOfSegTOF, 1, NumOfSegTOF + 1, 300, -300, 300.);

  HB2(140, "xSdcOut % TofSeg", NumOfSegTOF, 1, NumOfSegTOF+1, 400, -1000, 1000.);
  HB2(141, "TofPos - xSdcOut % TofSeg", NumOfSegTOF, 1, NumOfSegTOF + 1, 400, -200, 200);
  HB2(142, "xSdcOut % TofPos", 40, -1500, 1500, 400, -1000, 1000.);
  HB1(143, "TofPos - xSdcOut", 400, -200, 200);
  HB2(144, "TofPos - xSdcOut % TofPos", 40, -1500, 1500., 400, -200, 200);
  HB2(145, "ySdcOut % dtTof", 300, -15., 15., 400, -1000, 1000.);
  HB1(146, "ytTof - ySdcOut", 300, -300., 300.);
  HB2(147, "ytTof - ySdcOut % dtTof", 300, -15., 15., 300, -300, 300.);
  HB2(148, "ytTof - ySdcOut % TofSeg", NumOfSegTOF, 1, NumOfSegTOF + 1, 300, -300, 300.);

  for(Int_t i=1; i<=NumOfSegTOF; ++i){
    HB1(100*i+50, Form("TofPos - xSdcOut, tof_%d",i), 400, -200, 200);
    HB2(100*i+51, Form("ySdcOut %% dtTof, tof_%d",i), 300, -15., 15., 400, -1000, 1000.);
    HB1(100*i+52, Form("TofPos - xSdcOut, tof_%d",i), 400, -200, 200);
    HB2(100*i+53, Form("ySdcOut %% dtTof, tof_%d",i), 300, -15., 15., 400, -1000, 1000.);
  }

  // SDC1
  for(Int_t i=1; i<=NumOfLayersSDC1; ++i){
    Int_t nwire  = MaxWireSDC1;
    Int_t nbindt = NBinDTSDC1;
    Double_t mindt = MinDTSDC1;
    Double_t maxdt = MaxDTSDC1;
    Int_t nbindl = NBinDLSDC1;
    Double_t mindl = MinDLSDC1;
    Double_t maxdl = MaxDLSDC1;

    TString title1 = Form("TDC Sdc1_%d", i);
    TString title2 = Form("TOT Sdc1_%d", i);
    TString title8 = Form("Kurama Residual Sdc1_%d (0<theta<15)", i);
    TString title9 = Form("SdcIn Residual Sdc1_%d (0<theta<15)", i);
    TString title10 = Form("SdcIn Exclusive Residual Sdc1_%d (0<theta<15)", i);
    TString title11 = Form("HitPat Sdc1_%d", i);
    TString title12 = Form("DriftTime Sdc1_%d", i);
    TString title13 = Form("DriftLength Sdc1_%d", i);
    TString title14 = Form("SdcIn Residual Sdc1_%d", i);
    TString title15 = Form("Residual Sdc1_%d", i);
    TString title16 = Form("DriftLength%%DriftTime Sdc1_%d", i);
    TString title54 = Form("Residual Sdc1_%d (SdcOut)", i);
    TString title56 = Form("DriftLength%%DriftTime Sdc1_%d (SdcOut)", i);
    TString title64 = Form("Wire Position - SdcIn LocalCalPosition Sdc1_%d", i);
    TString title65 = Form("Wire Position - SdcOut LocalCalPosition Sdc1_%d", i);
    TString title66 = Form("SdcIn LocalCalPosition%%Wire Position Sdc1_%d", i);
    TString title67 = Form("SdcOut LocalCalPosition%%Wire Position Sdc1_%d", i);
    TString title68 = Form("Wire Position - SdcIn LocalCalPosition%%Wire Position Sdc1_%d", i);
    TString title69 = Form("Wire Position - SdcOut LocalCalPosition%%Wire Position Sdc1_%d", i);

    HB1(100*i+1, title1, NbinSdcTdc, MinSdcTdc, MaxSdcTdc);
    HB1(100*i+2, title2, 500, 0., 500.);
    HB1(100*i+8, title8, 500, -5.0, 5.0);
    HB1(100*i+9, title9, 500, -5.0, 5.0);
    HB1(100*i+10, title10, 500, -5.0, 5.0);
    HB1(100*i+11, title11, 130, 0., 130.);
    HB1(100*i+12, title12, NBinDTSDC1, MinDTSDC1, MaxDTSDC1);
    HB1(100*i+13, title13, NBinDLSDC1, MinDLSDC1, MaxDLSDC1);
    HB1(100*i+14, title14, 500, -5.0, 5.0);
    HB1(100*i+15, title15, 500, -5.0, 5.0);
    HBProf(100*i+16, title16, nbindt, mindt, maxdt, mindl, maxdl);
    HB2(100*i+18, title16, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB2(100*i+19, title16, nbindt, mindt, maxdt, 2*nbindl,-maxdl, maxdl);
    HBProf(100*i+20, title16, nbindt, mindt, maxdt, mindl, maxdl);
    HB2(100*i+22, title16, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB2(100*i+23, title16, nbindt, mindt, maxdt, 2*nbindl, -maxdl, maxdl);
    HB1(100*i+54, title54, 500, -5.0, 5.0);
    HBProf(100*i+56, title56, nbindt, mindt, maxdt, mindl, maxdl);
    HB2(100*i+58, title56, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB2(100*i+59, title56, nbindt, mindt, maxdt, 2*nbindl,-maxdl, maxdl);
    HBProf(100*i+60, title56, nbindt, mindt, maxdt, mindl, maxdl);
    HB2(100*i+62, title56, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB2(100*i+63, title56, nbindt, mindt, maxdt, 2*nbindl, -maxdl, maxdl);
    //HB1(100*i+64, title64, 1500, -15.0, 15.0);
    //HB1(100*i+65, title65, 1500, -15.0, 15.0);
    HB1(100*i+64, title64, 1500, -300.0, 300.0);
    HB1(100*i+65, title65, 1500, -300.0, 300.0);
    HB2(100*i+66, title66, 1000, -1000.0, 1000.0, 1000, -1000.0, 1000.0);
    HB2(100*i+67, title67, 1000, -1000.0, 1000.0, 1000, -1000.0, 1000.0);
    HB2(100*i+68, title68, 1000, -1000.0, 1000.0, 150, -15.0, 15.0);
    HB2(100*i+69, title69, 1000, -1000.0, 1000.0, 150, -15.0, 15.0);

    for(Int_t j=1; j<=nwire; j++){
      TString title_t0 = Form("t0 Layer %d Wire #%4d", i, j);
      HB1(10000*i+j, title_t0, NbinSdcTime, MinSdcTime, MaxSdcTime);
    }
  }

  // SDC2
  for(Int_t i=NumOfLayersSDC1+1; i<=NumOfLayersSdcIn; ++i){
    Int_t nwire = (i==NumOfLayersSDC1+1 || i==NumOfLayersSDC1+2) ? MaxWireSDC2X : MaxWireSDC2Y;
    Int_t nbindt = NBinDTSDC2;
    Double_t mindt = MinDTSDC2;
    Double_t maxdt = MaxDTSDC2;
    Int_t nbindl = NBinDLSDC2;
    Double_t mindl = MinDLSDC2;
    Double_t maxdl = MaxDLSDC2;

    TString title1 = Form("TDC Sdc2_%d", i);
    TString title2 = Form("TOT Sdc2_%d", i);
    TString title8 = Form("Kurama Residual Sdc2_%d (0<theta<15)", i);
    TString title9 = Form("SdcIn Residual Sdc2_%d (0<theta<15)", i);
    TString title10 = Form("SdcIn Exclusive Residual Sdc2_%d (0<theta<15)", i);
    TString title11 = Form("HitPat Sdc2_%d", i);
    TString title12 = Form("DriftTime Sdc2_%d", i);
    TString title13 = Form("DriftLength Sdc2_%d", i);
    TString title14 = Form("SdcIn Residual Sdc2_%d", i);
    TString title15 = Form("Residual Sdc2_%d", i);
    TString title16 = Form("DriftLength%%DriftTime Sdc2_%d", i);
    TString title54 = Form("Residual Sdc2_%d (SdcOut)", i);
    TString title56 = Form("DriftLength%%DriftTime Sdc2_%d (SdcOut)", i);
    TString title64 = Form("Wire Position - SdcIn LocalCalPosition Sdc2_%d", i);
    TString title65 = Form("Wire Position - SdcOut LocalCalPosition Sdc2_%d", i);
    TString title66 = Form("SdcIn LocalCalPosition%%Wire Position Sdc2_%d", i);
    TString title67 = Form("SdcOut LocalCalPosition%%Wire Position Sdc2_%d", i);
    TString title68 = Form("Wire Position - SdcIn LocalCalPosition%%Wire Position Sdc2_%d", i);
    TString title69 = Form("Wire Position - SdcOut LocalCalPosition%%Wire Position Sdc2_%d", i);

    HB1(100*i+1, title1, NbinSdcTdc, MinSdcTdc, MaxSdcTdc);
    HB1(100*i+2, title2, 500, 0., 500.);
    HB1(100*i+8, title8, 500, -5.0, 5.0);
    HB1(100*i+9, title9, 500, -5.0, 5.0);
    HB1(100*i+10, title10, 500, -5.0, 5.0);
    HB1(100*i+11, title11, 130, 0., 130.);
    HB1(100*i+12, title12, NBinDTSDC2, MinDTSDC2, MaxDTSDC2);
    HB1(100*i+13, title13, NBinDLSDC2, MinDLSDC2, MaxDLSDC2);
    HB1(100*i+14, title14, 500, -5.0, 5.0);
    HB1(100*i+15, title15, 500, -5.0, 5.0);
    HBProf(100*i+16, title16, nbindt, mindt, maxdt, mindl, maxdl);
    HB2(100*i+18, title16, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB2(100*i+19, title16, nbindt, mindt, maxdt, 2*nbindl,-maxdl, maxdl);
    HBProf(100*i+20, title16, nbindt, mindt, maxdt, mindl, maxdl);
    HB2(100*i+22, title16, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB2(100*i+23, title16, nbindt, mindt, maxdt, 2*nbindl, -maxdl, maxdl);
    HB1(100*i+54, title54, 500, -5.0, 5.0);
    HBProf(100*i+56, title56, nbindt, mindt, maxdt, mindl, maxdl);
    HB2(100*i+58, title56, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB2(100*i+59, title56, nbindt, mindt, maxdt, 2*nbindl,-maxdl, maxdl);
    HBProf(100*i+60, title56, nbindt, mindt, maxdt, mindl, maxdl);
    HB2(100*i+62, title56, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB2(100*i+63, title56, nbindt, mindt, maxdt, 2*nbindl, -maxdl, maxdl);
    //HB1(100*i+64, title64, 1500, -15.0, 15.0);
    //HB1(100*i+65, title65, 1500, -15.0, 15.0);
    HB1(100*i+64, title64, 1500, -300.0, 300.0);
    HB1(100*i+65, title65, 1500, -300.0, 300.0);
    HB2(100*i+66, title66, 1000, -1000.0, 1000.0, 1000, -1000.0, 1000.0);
    HB2(100*i+67, title67, 1000, -1000.0, 1000.0, 1000, -1000.0, 1000.0);
    HB2(100*i+68, title68, 1000, -1000.0, 1000.0, 150, -15.0, 15.0);
    HB2(100*i+69, title69, 1000, -1000.0, 1000.0, 150, -15.0, 15.0);
    for(Int_t j=1; j<=nwire; j++){
      TString title_t0 = Form("t0 Layer %d Wire #%4d", i, j);
      HB1(10000*i+j, title_t0, NbinSdcTime, MinSdcTime, MaxSdcTime);
    }
  }

  // SDC3
  for(Int_t i=NumOfLayersSdcIn+1; i<=(NumOfLayersSdcIn+NumOfLayersSDC3); ++i){
    Int_t nwire  = MaxWireSDC3;
    Int_t nbindt = NBinDTSDC3;
    Double_t mindt = MinDTSDC3;
    Double_t maxdt = MaxDTSDC3;
    Int_t nbindl = NBinDLSDC3;
    Double_t mindl = MinDLSDC3;
    Double_t maxdl = MaxDLSDC3;

    TString title1 = Form("TDC Sdc3_%d", i);
    TString title2 = Form("TOT 1st Sdc3_%d", i);
    TString title8 = Form("Kurama Residual Sdc3_%d (0<theta<15)", i);
    TString title9 = Form("SdcOut Residual Sdc3_%d (0<theta<15)", i);
    TString title10 = Form("SdcOut Exclusive Residual Sdc3_%d (0<theta<15)", i);
    TString title11 = Form("HitPat Sdc3_%d", i);
    TString title12 = Form("DriftTime Sdc3_%d", i);
    TString title13 = Form("DriftLength Sdc3_%d", i);
    TString title14 = Form("SdcOut Residual Sdc3_%d", i);
    TString title15 = Form("Residual Sdc3_%d", i);
    TString title16 = Form("DriftLength%%DriftTime Sdc3_%d", i);
    TString title54 = Form("Residual Sdc3_%d (SdcIn)", i);
    TString title56 = Form("DriftLength%%DriftTime Sdc3_%d (SdcIn)", i);
    TString title64 = Form("Wire Position - SdcIn LocalCalPosition Sdc3_%d", i);
    TString title65 = Form("Wire Position - SdcOut LocalCalPosition Sdc3_%d", i);
    TString title66 = Form("SdcIn LocalCalPosition%%Wire Position Sdc3_%d", i);
    TString title67 = Form("SdcOut LocalCalPosition%%Wire Position Sdc3_%d", i);
    TString title68 = Form("Wire Position - SdcIn LocalCalPosition%%Wire Position Sdc3_%d", i);
    TString title69 = Form("Wire Position - SdcOut LocalCalPosition%%Wire Position Sdc3_%d", i);

    HB1(100*i+1, title1, NbinSdcTdc, MinSdcTdc, MaxSdcTdc);
    HB1(100*i+2, title2, 500, 0., 500.);
    HB1(100*i+8, title8, 500, -5.0, 5.0);
    HB1(100*i+9, title9, 500, -5.0, 5.0);
    HB1(100*i+10, title10, 500, -5.0, 5.0);
    HB1(100*i+11, title11, 130, 0., 130.);
    HB1(100*i+12, title12, NBinDTSDC3, MinDTSDC3, MaxDTSDC3);
    HB1(100*i+13, title13, NBinDLSDC3, MinDLSDC3, MaxDLSDC3);
    HB1(100*i+14, title14, 500, -5.0, 5.0);
    HB1(100*i+15, title15, 500, -5.0, 5.0);
    HBProf(100*i+16, title16, nbindt, mindt, maxdt, mindl, maxdl);
    HB2(100*i+18, title16, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB2(100*i+19, title16, nbindt, mindt, maxdt, 2*nbindl,-maxdl, maxdl);
    HBProf(100*i+20, title16, nbindt, mindt, maxdt, mindl, maxdl);
    HB2(100*i+22, title16, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB2(100*i+23, title16, nbindt, mindt, maxdt, 2*nbindl, -maxdl, maxdl);
    HB1(100*i+54, title54, 500, -5.0, 5.0);
    HBProf(100*i+56, title56, nbindt, mindt, maxdt, mindl, maxdl);
    HB2(100*i+58, title56, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB2(100*i+59, title56, nbindt, mindt, maxdt, 2*nbindl,-maxdl, maxdl);
    HBProf(100*i+60, title56, nbindt, mindt, maxdt, mindl, maxdl);
    HB2(100*i+62, title56, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB2(100*i+63, title56, nbindt, mindt, maxdt, 2*nbindl, -maxdl, maxdl);
    //HB1(100*i+64, title64, 1500, -15.0, 15.0);
    //HB1(100*i+65, title65, 1500, -15.0, 15.0);
    HB1(100*i+64, title64, 1500, -300.0, 300.0);
    HB1(100*i+65, title65, 1500, -300.0, 300.0);
    HB2(100*i+66, title66, 1000, -1000.0, 1000.0, 1000, -1000.0, 1000.0);
    HB2(100*i+67, title67, 1000, -1000.0, 1000.0, 1000, -1000.0, 1000.0);
    HB2(100*i+68, title68, 1000, -1000.0, 1000.0, 150, -15.0, 15.0);
    HB2(100*i+69, title69, 1000, -1000.0, 1000.0, 150, -15.0, 15.0);
    for(Int_t j=1; j<=nwire; j++){
      TString title_t0 = Form("t0 Layer %d Wire #%4d", i, j);
      HB1(10000*i+j, title_t0, NbinSdcTime, MinSdcTime, MaxSdcTime);
    }
  }

  // SDC4
  for(Int_t i=NumOfLayersSdcIn+NumOfLayersSDC3+1;i<=NumOfLayersSdcIn+NumOfLayersSdcOut ; ++i){
    Int_t nwire = (i==NumOfLayersSdcIn+NumOfLayersSDC3+1 || i==NumOfLayersSdcIn+NumOfLayersSDC3+2) ? MaxWireSDC4Y : MaxWireSDC4X;
    Int_t nbindt = NBinDTSDC4;
    Double_t mindt = MinDTSDC4;
    Double_t maxdt = MaxDTSDC4;
    Int_t nbindl = NBinDLSDC4;
    Double_t mindl = MinDLSDC4;
    Double_t maxdl = MaxDLSDC4;

    TString title1 = Form("TDC Sdc4_%d", i);
    TString title2 = Form("TOT Sdc4_%d", i);
    TString title8 = Form("Kurama Residual Sdc4_%d (0<theta<15)", i);
    TString title9 = Form("SdcOut Residual Sdc4_%d (0<theta<15)", i);
    TString title10 = Form("SdcOut Exclusive Residual Sdc4_%d (0<theta<15)", i);
    TString title11 = Form("HitPat Sdc4_%d", i);
    TString title12 = Form("DriftTime Sdc4_%d", i);
    TString title13 = Form("DriftLength Sdc4_%d", i);
    TString title14 = Form("SdcOut Residual Sdc4_%d", i);
    TString title15 = Form("Residual Sdc4_%d", i);
    TString title16 = Form("DriftLength%%DriftTime Sdc4_%d", i);
    TString title54 = Form("Residual Sdc4_%d (SdcIn)", i);
    TString title56 = Form("DriftLength%%DriftTime Sdc4_%d (SdcIn)", i);
    TString title64 = Form("Wire Position - SdcIn LocalCalPosition Sdc4_%d", i);
    TString title65 = Form("Wire Position - SdcOut LocalCalPosition Sdc4_%d", i);
    TString title66 = Form("SdcIn LocalCalPosition%%Wire Position Sdc4_%d", i);
    TString title67 = Form("SdcOut LocalCalPosition%%Wire Position Sdc4_%d", i);
    TString title68 = Form("Wire Position - SdcIn LocalCalPosition%%Wire Position Sdc4_%d", i);
    TString title69 = Form("Wire Position - SdcOut LocalCalPosition%%Wire Position Sdc4_%d", i);

    HB1(100*i+1, title1, NbinSdcTdc, MinSdcTdc, MaxSdcTdc);
    HB1(100*i+2, title2, 500, 0., 500.);
    HB1(100*i+8, title8, 500, -5.0, 5.0);
    HB1(100*i+9, title9, 500, -5.0, 5.0);
    HB1(100*i+10, title10, 500, -5.0, 5.0);
    HB1(100*i+11, title11, 130, 0., 130.);
    HB1(100*i+12, title12, NBinDTSDC4, MinDTSDC4, MaxDTSDC4);
    HB1(100*i+13, title13, NBinDLSDC4, MinDLSDC4, MaxDLSDC4);
    HB1(100*i+14, title14, 500, -5.0, 5.0);
    HB1(100*i+15, title15, 500, -5.0, 5.0);
    HBProf(100*i+16, title16, nbindt, mindt, maxdt, mindl, maxdl);
    HB2(100*i+18, title16, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB2(100*i+19, title16, nbindt, mindt, maxdt, 2*nbindl,-maxdl, maxdl);
    HBProf(100*i+20, title16, nbindt, mindt, maxdt, mindl, maxdl);
    HB2(100*i+22, title16, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB2(100*i+23, title16, nbindt, mindt, maxdt, 2*nbindl, -maxdl, maxdl);
    HB1(100*i+54, title54, 500, -5.0, 5.0);
    HBProf(100*i+56, title56, nbindt, mindt, maxdt, mindl, maxdl);
    HB2(100*i+58, title56, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB2(100*i+59, title56, nbindt, mindt, maxdt, 2*nbindl,-maxdl, maxdl);
    HBProf(100*i+60, title56, nbindt, mindt, maxdt, mindl, maxdl);
    HB2(100*i+62, title56, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB2(100*i+63, title56, nbindt, mindt, maxdt, 2*nbindl, -maxdl, maxdl);
    //HB1(100*i+64, title64, 1500, -15.0, 15.0);
    //HB1(100*i+65, title65, 1500, -15.0, 15.0);
    HB1(100*i+64, title64, 1500, -300.0, 300.0);
    HB1(100*i+65, title65, 1500, -300.0, 300.0);
    HB2(100*i+66, title66, 1000, -1000.0, 1000.0, 1000, -1000.0, 1000.0);
    HB2(100*i+67, title67, 1000, -1000.0, 1000.0, 1000, -1000.0, 1000.0);
    HB2(100*i+68, title68, 1000, -1000.0, 1000.0, 150, -15.0, 15.0);
    HB2(100*i+69, title69, 1000, -1000.0, 1000.0, 150, -15.0, 15.0);
    for(Int_t j=1; j<=nwire; j++){
      TString title_t0 = Form("t0 Layer %d Wire #%4d", i, j);
      HB1(10000*i+j, title_t0, NbinSdcTime, MinSdcTime, MaxSdcTime);
    }
  }

  // TOF in SdcOut/KuramaTracking
  for(Int_t i=NumOfLayersSdcIn+NumOfLayersSdcOut+1;
      i<=NumOfLayersSdcIn+NumOfLayersSdcOut+4; ++i){
    TString title11 = Form("HitPat Tof%d", i-(NumOfLayersSdcIn+NumOfLayersSdcOut));
    HB1(100*i+11, title11, 30, 0., 30.);
  }

  //Tree
  HBTree("kurama","tree of KuramaTracking");
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

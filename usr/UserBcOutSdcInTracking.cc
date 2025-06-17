// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

#include <TMath.h>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DCRawHit.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "KuramaLib.hh"
#include "MathTools.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UnpackerManager.hh"

#define HodoCut 0

namespace
{
using namespace root;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& zK18Target = gGeom.LocalZ("K18Target");
}

//_____________________________________________________________________________
class UserBcOutSdcInTracking : public VEvent
{
private:
  RawData*      rawData;
  DCAnalyzer*   DCAna;
  HodoAnalyzer* hodoAna;

public:
  UserBcOutSdcInTracking();
  ~UserBcOutSdcInTracking();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserBcOutSdcInTracking::ClassName()
{
  static TString s_name("UserBcOutSdcInTracking");
  return s_name;
}

//_____________________________________________________________________________
UserBcOutSdcInTracking::UserBcOutSdcInTracking()
  : VEvent(),
    rawData(new RawData),
    DCAna(new DCAnalyzer),
    hodoAna(new HodoAnalyzer)
{
}

//_____________________________________________________________________________
UserBcOutSdcInTracking::~UserBcOutSdcInTracking()
{
  if(rawData) delete rawData;
  if(DCAna) delete DCAna;
  if(hodoAna) delete hodoAna;
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;
  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  Double_t btof;
  Double_t time0;

  Int_t nhBh1;
  Double_t tBh1[MaxHits];
  Double_t deBh1[MaxHits];

  Int_t nhBh2;
  Double_t tBh2[MaxHits];
  Double_t deBh2[MaxHits];
  Double_t Bh2Seg[MaxHits];

  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;

  Int_t ntrack;
  Int_t nh[MaxHits];
  Double_t chisqr[MaxHits];
  Double_t x_Bh2[MaxHits];
  Double_t x0[MaxHits];
  Double_t y0[MaxHits];
  Double_t u0[MaxHits];
  Double_t v0[MaxHits];

  Int_t ntBcOut;
  Int_t nhBcOut[MaxHits];
  Double_t chisqrBcOut[MaxHits];
  Double_t x0BcOut[MaxHits];
  Double_t y0BcOut[MaxHits];
  Double_t u0BcOut[MaxHits];
  Double_t v0BcOut[MaxHits];

  Int_t ntSdcIn;
  Int_t nhSdcIn[MaxHits];
  Double_t chisqrSdcIn[MaxHits];
  Double_t x0SdcIn[MaxHits];
  Double_t y0SdcIn[MaxHits];
  Double_t u0SdcIn[MaxHits];
  Double_t v0SdcIn[MaxHits];

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum     =  0;
  ntrack    =  0;
  nhBh2     =  0;
  nhBh1     =  0;
  Time0Seg  = -1;
  deTime0   = -1;
  Time0     = qnan;
  CTime0    = qnan;

  time0 = qnan;
  btof  = qnan;

  ntBcOut = 0;
  ntSdcIn = 0;
  for(Int_t it=0; it<MaxHits; it++){
    nh[it]     = 0;
    tBh1[it]   = qnan;
    deBh1[it]  = qnan;
    Bh2Seg[it] = -1;
    tBh2[it]   = qnan;
    deBh2[it]  = qnan;
    chisqr[it] = qnan;
    x0[it]     = qnan;
    y0[it]     = qnan;
    u0[it]     = qnan;
    v0[it]     = qnan;

    nhBcOut[it] = 0;
    chisqrBcOut[it] = qnan;
    x0BcOut[it] = qnan;
    y0BcOut[it] = qnan;
    u0BcOut[it] = qnan;
    v0BcOut[it] = qnan;

    nhSdcIn[it] = 0;
    chisqrSdcIn[it] = qnan;
    x0SdcIn[it] = qnan;
    y0SdcIn[it] = qnan;
    u0SdcIn[it] = qnan;
    v0SdcIn[it] = qnan;
  }

  for(Int_t it=0; it<NumOfSegTrig; it++){
    trigpat[it] = -1;
    trigflag[it] = -1;
  }
}

//_____________________________________________________________________________
namespace root
{
Event event;
TH1   *h[MaxHist];
TTree *tree;
}

//_____________________________________________________________________________
Bool_t
UserBcOutSdcInTracking::ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserBcOutSdcInTracking::ProcessingNormal()
{
#if HodoCut
  static const auto MinDeBH2 = gUser.GetParameter("DeBH2", 0);
  static const auto MaxDeBH2 = gUser.GetParameter("DeBH2", 1);
  static const auto MinDeBH1 = gUser.GetParameter("DeBH1", 0);
  static const auto MaxDeBH1 = gUser.GetParameter("DeBH1", 1);
  static const auto MinBeamToF = gUser.GetParameter("BTOF", 0);
  static const auto MaxBeamToF = gUser.GetParameter("BTOF", 1);
#endif
  static const auto MinTotBcOut = gUser.GetParameter("MinTotBcOut");
  static const auto MinTotSDC1 = gUser.GetParameter("MinTotSDC1");
  static const auto MinTotSDC2 = gUser.GetParameter("MinTotSDC2");

  static const auto MaxChisqrBcOut = gUser.GetParameter("MaxChisqrBcOut");
  static const auto MaxChisqrSdcIn = gUser.GetParameter("MaxChisqrSdcIn");

  static const auto MaxMultiHitBcOut = gUser.GetParameter("MaxMultiHitBcOut");
  static const auto MaxMultiHitSdcIn = gUser.GetParameter("MaxMultiHitSdcIn");

  static const auto MinDriftTimeBC34 = gUser.GetParameter("DriftTimeBC34", 0);
  static const auto MaxDriftTimeBC34 = gUser.GetParameter("DriftTimeBC34", 1);
  static const auto MinDriftTimeSDC1 = gUser.GetParameter("DriftTimeSDC1", 0);
  static const auto MaxDriftTimeSDC1 = gUser.GetParameter("DriftTimeSDC1", 1);
  static const auto MinDriftTimeSDC2 = gUser.GetParameter("DriftTimeSDC2", 0);
  static const auto MaxDriftTimeSDC2 = gUser.GetParameter("DriftTimeSDC2", 1);

  rawData->DecodeHits();

  event.evnum = gUnpacker.get_event_number();

  // Trigger Flag
  std::bitset<NumOfSegTrig> trigger_flag;
  for(const auto& hit: rawData->GetTrigRawHC()){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc1();
    if(tdc > 0){
      event.trigpat[trigger_flag.count()] = seg;
      event.trigflag[seg] = tdc;
      trigger_flag.set(seg);
    }
  }

  HF1(1, 0.);

  if(trigger_flag[trigger::kSpillEnd]) return true;

  HF1(1, 1.);

  //////////////BH2 time 0
  hodoAna->DecodeBH2Hits(rawData);
  Int_t nhBh2 = hodoAna->GetNHitsBH2();
  event.nhBh2 = nhBh2;
#if HodoCut
  if(nhBh2==0) return true;
#endif
  HF1(1, 2);

  Double_t time0 = qnan;
  //////////////BH2 Analysis
  for(Int_t i=0; i<nhBh2; ++i){
    BH2Hit* hit = hodoAna->GetHitBH2(i);
    if(!hit) continue;
    Double_t seg = hit->SegmentId()+1;
    Double_t cmt = hit->CMeanTime();
    Double_t dE  = hit->DeltaE();

#if HodoCut
    if(dE<MinDeBH2 || MaxDeBH2<dE) continue;
#endif
    event.tBh2[i]   = cmt;
    event.deBh2[i]  = dE;
    event.Bh2Seg[i] = seg;
  }

  BH2Cluster *cl_time0 = hodoAna->GetTime0BH2Cluster();
  if(cl_time0){
    event.Time0Seg = cl_time0->MeanSeg()+1;
    event.deTime0  = cl_time0->DeltaE();
    event.Time0    = cl_time0->Time0();
    event.CTime0   = cl_time0->CTime0();
    time0          = cl_time0->CTime0();
  } else {
#if HodoCut
    return true;
#endif
  }

  HF1(1, 3.);

  //////////////BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData);
  Int_t nhBh1 = hodoAna->GetNHitsBH1();
  event.nhBh1 = nhBh1;
#if HodoCut
  if(nhBh1==0) return true;
#endif

  HF1(1, 4);

  for(Int_t i=0; i<nhBh1; ++i){
    Hodo2Hit *hit = hodoAna->GetHitBH1(i);
    if(!hit) continue;
    Double_t cmt = hit->CMeanTime();
    Double_t dE  = hit->DeltaE();
#if HodoCut
    Double_t btof = cmt - time0;
    if(dE<MinDeBH1 || MaxDeBH1<dE) continue;
    if(btof<MinBeamToF || MaxBeamToF<btof) continue;
#endif
    event.tBh1[i]  = cmt;
    event.deBh1[i] = dE;
  }

  Double_t btof0 = qnan;
  HodoCluster* cl_btof0 = event.Time0Seg > 0 ?
	  hodoAna->GetBtof0BH1Cluster(event.CTime0) : nullptr;
  if(cl_btof0) btof0 = cl_btof0->CMeanTime() - time0;
  event.btof = btof0;

  HF1(1, 5.);

  DCAna->DecodeRawHits(rawData);

  ///// BcOut number of hit layer
  DCAna->TotCutBCOut(MinTotBcOut);
  DCAna->TotCutBCOut(MinTotBcOut);
  DCAna->DriftTimeCutBC34(MinDriftTimeBC34, MaxDriftTimeBC34);
  // DCAna->MakeBH2DCHit(event.Time0Seg-1);

  Double_t multi_BcOut = 0.;
  for(Int_t layer=1; layer<=NumOfLayersBcOut; ++layer){
    multi_BcOut += DCAna->GetBcOutHC(layer).size();
  }
  if(multi_BcOut/Double_t(NumOfLayersBcOut) > MaxMultiHitBcOut){
    std::cout << "#W " << __FILE__ << " L" << __LINE__ << std::endl
	      << "multi_BcOut is too many: " << multi_BcOut << std::endl;
    return true;
  }

  HF1(1, 11.);

  ///// SdcIn number of hit layer
  DCAna->TotCutSDC1(MinTotSDC1);
  DCAna->TotCutSDC2(MinTotSDC2);
  DCAna->DriftTimeCutSDC1(MinDriftTimeSDC1, MaxDriftTimeSDC1);
  DCAna->DriftTimeCutSDC2(MinDriftTimeSDC2, MaxDriftTimeSDC2);

  Double_t multi_SdcIn = 0.;
  for(Int_t layer=1; layer<=NumOfLayersSdcIn; ++layer){
    multi_SdcIn += DCAna->GetSdcInHC(layer).size();
  }
  if(multi_SdcIn/Double_t(NumOfLayersSdcIn) > MaxMultiHitSdcIn){
    std::cout << "#W " << __FILE__ << " L" << __LINE__ << std::endl
	      << "multi_SdcIn is too many: " << multi_SdcIn << std::endl;
    return true;
  }
  HF1(1, 12.);

  //  std::cout << "==========TrackSearch BcOutSdcIn============" << std::endl;
  DCAna->TrackSearchBcOutSdcIn();
  auto track_cont = DCAna->GetTrackContainerBcOutSdcIn();
  Int_t nt = track_cont.size();
  event.ntrack = nt;
  HF1(10, Double_t(nt));
  for(Int_t it=0; it<nt; ++it){
    const auto& track = track_cont.at(it);
    // track->Print();
    Int_t nh = track->GetNHit();
    Double_t chisqr=track->GetChiSquare();
    Double_t x0=track->GetX0(), y0=track->GetY0();
    Double_t u0=track->GetU0(), v0=track->GetV0();
    Double_t theta = track->GetTheta();
    event.nh[it] = nh;
    event.chisqr[it] = chisqr;
    event.x0[it] = x0;
    event.y0[it] = y0;
    event.u0[it] = u0;
    event.v0[it] = v0;

    HF1(11, Double_t(nh));
    HF1(12, chisqr);
    HF1(14, x0); HF1(15, y0);
    HF1(16, u0); HF1(17, v0);
    HF2(18, x0, u0); HF2(19, y0, v0);
    HF2(20, x0, y0);

    Double_t xtgt=track->GetX(zK18Target), ytgt=track->GetY(zK18Target);
    Double_t utgt=u0, vtgt=v0;
    HF1(21, xtgt); HF1(22, ytgt);
    HF1(23, utgt); HF1(24, vtgt);
    HF2(25, xtgt, utgt); HF2(26, ytgt, vtgt);
    HF2(27, xtgt, ytgt);
    /*
    Double_t Xangle = -1000*TMath::ATan(u0);
    HF2(51, -track->GetX(245.), Xangle);
    HF2(52, -track->GetX(600.), Xangle);
    HF2(53, -track->GetX(1200.), Xangle);
    HF2(54, -track->GetX(1600.), Xangle);

    HF2(61, track->GetX(600.), track->GetX(245.));
    HF2(62, track->GetX(1200.), track->GetX(245.));
    HF2(63, track->GetX(1600.), track->GetX(245.));
    */
    for(Int_t ih=0; ih<nh; ++ih){
      const auto hit = track->GetHit(ih);
      Int_t layerId = hit->GetLayer();//-112;
      if(layerId > 112) layerId -= 112;
      else layerId += 12;
      HF1(13, layerId);
      Double_t wire=hit->GetWire();
      Double_t dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1(100*layerId+11, wire-0.5);
      HF1(100*layerId+12, dt);
      HF1(100*layerId+13, dl);
      Double_t xcal=hit->GetXcal(), ycal=hit->GetYcal();
      Double_t pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      if(layerId > 12){
        res += 50.*TMath::Cos(hit->GetTiltAngle()*TMath::DegToRad());
      }
      HF1(100*layerId+14, pos);
      HF1(100*layerId+15, res);
      HF2(100*layerId+16, pos, res);
      HF2(100*layerId+17, xcal, ycal);
      Double_t wp=hit->GetWirePosition();
      Double_t sign=1.;
      if(pos-wp<0.) sign=-1;
      HF2(100*layerId+18, sign*dl, res);
      Double_t xlcal=hit->GetLocalCalPos();
      HF2(100*layerId+19, dt, xlcal-wp);
      HFProf(100*layerId+20, dt, dl);
      HF2(100*layerId+22, dt, dl);
      Double_t tot = hit->GetTot();
      HF1(100*layerId+40, tot);
      if (theta>=0 && theta<15)
        HF1(100*layerId+71, res);
      else if (theta>=15 && theta<30)
        HF1(100*layerId+72, res);
      else if (theta>=30 && theta<45)
        HF1(100*layerId+73, res);
      else if (theta>=45)
        HF1(100*layerId+74, res);
    }
  }

  HF1(1, 15.);

#if 1
  // Bc Out
  //  std::cout << "==========TrackSearch BcOut============" << std::endl;
  DCAna->TrackSearchBcOut();
  DCAna->ChiSqrCutBcOut(MaxChisqrBcOut);

  Int_t ntBcOut = DCAna->GetNtracksBcOut();
  event.ntBcOut = ntBcOut;
  for(Int_t it=0; it<ntBcOut; ++it){
    DCLocalTrack *tp=DCAna->GetTrackBcOut(it);
    Int_t nh=tp->GetNHit();
    Double_t chisqr=tp->GetChiSquare();
    Double_t x0=tp->GetX0(), y0=tp->GetY0();
    Double_t u0=tp->GetU0(), v0=tp->GetV0();
    event.nhBcOut[it] = nh;
    event.chisqrBcOut[it] = chisqr;
    event.x0BcOut[it] = x0;
    event.y0BcOut[it] = y0;
    event.u0BcOut[it] = u0;
    event.v0BcOut[it] = v0;
  }
#endif

#if 1
  // std::cout << "==========TrackSearch SdcIn============" << std::endl;
  DCAna->TrackSearchSdcIn();
  DCAna->ChiSqrCutSdcIn(MaxChisqrSdcIn);

  Int_t ntSdcIn = DCAna->GetNtracksSdcIn();
  event.ntSdcIn = ntSdcIn;
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
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
UserBcOutSdcInTracking::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserBcOutSdcInTracking;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  const Int_t nl = NumOfLayersBcOut+NumOfLayersSdcIn;

  const Int_t    NbinBcOutDT = 96;
  const Double_t MinBcOutDT  = -10.;
  const Double_t MaxBcOutDT  = 70.;
  const Int_t    NbinBcOutDL = 80;
  const Double_t MinBcOutDL  = -0.5;
  const Double_t MaxBcOutDL  = 3.5;

  const Int_t    NbinSDC1DT = 180;
  const Double_t MinSDC1DT  = -30.;
  const Double_t MaxSDC1DT  = 120.;
  const Int_t    NbinSDC1DL = 80;
  const Double_t MinSDC1DL  = -0.5;
  const Double_t MaxSDC1DL  = 3.5;

  const Int_t    NbinSDC2DT = 360;
  const Double_t MinSDC2DT  = -50.;
  const Double_t MaxSDC2DT  = 250.;
  const Int_t    NbinSDC2DL = 130;
  const Double_t MinSDC2DL  = -0.5;
  const Double_t MaxSDC2DL  = 6.;

  HB1(1, "Status", 20, 0., 20.);

  //***********************Chamber
  for(Int_t i=1; i<=nl; ++i){
    TString tag;
    Int_t nwire;
    Int_t nbindt;
    Double_t mindt;
    Double_t maxdt;
    Int_t nbindl;
    Double_t mindl;
    Double_t maxdl;
    if(i <= NumOfLayersBc){
      tag = "BC3";
      nwire = MaxWireBC3;
      nbindt = NbinBcOutDT;
      mindt = MinBcOutDT;
      maxdt = MaxBcOutDT;
      nbindl = NbinBcOutDL;
      mindl = MinBcOutDL;
      maxdl = MaxBcOutDL;
    }else if(i <= NumOfLayersBcOut){
      tag = "BC4";
      nwire = MaxWireBC4;
      nbindt = NbinBcOutDT;
      mindt = MinBcOutDT;
      maxdt = MaxBcOutDT;
      nbindl = NbinBcOutDL;
      mindl = MinBcOutDL;
      maxdl = MaxBcOutDL;
    }else if(i <= NumOfLayersBcOut+NumOfLayersSDC1){
      tag = "SDC1";
      nwire = MaxWireSDC1;
      nbindt = NbinSDC1DT;
      mindt = MinSDC1DT;
      maxdt = MaxSDC1DT;
      nbindl = NbinSDC1DL;
      mindl = MinSDC1DL;
      maxdl = MaxSDC1DL;
    }else if(i <= NumOfLayersBcOut+NumOfLayersSdcIn){
      tag = "SDC2";
      if(i <= NumOfLayersBcOut+NumOfLayersSDC1+2){
        nwire = MaxWireSDC2X;
      }else{
        nwire = MaxWireSDC2Y;
      }
      nbindt = NbinSDC2DT;
      mindt = MinSDC2DT;
      maxdt = MaxSDC2DT;
      nbindl = NbinSDC2DL;
      mindl = MinSDC2DL;
      maxdl = MaxSDC2DL;
    }

    // Tracking Histgrams
    TString title11 = Form("HitPat BcOutSdcIn L%2d [Track]", i);
    TString title12 = Form("DriftTime BcOutSdcIn L%2d [Track]", i);
    TString title13 = Form("DriftLength BcOutSdcIn L%2d [Track]", i);
    TString title14 = Form("Position BcOutSdcIn L%2d", i);
    TString title15 = Form("Residual BcOutSdcIn L%2d", i);
    TString title16 = Form("Resid%%Pos BcOutSdcIn L%2d", i);
    TString title17 = Form("Y%%Xcal BcOutSdcIn L%2d", i);
    TString title18 = Form("Res%%dl BcOutSdcIn L%2d", i);
    TString title19 = Form("HitPos%%DriftTime BcOutSdcIn L%2d", i);
    TString title20 = Form("DriftLength%%DriftTime BcOutSdcIn L%2d", i);
    TString title21 = title15 + " [w/o Self]";
    TString title22 = title20;
    TString title40 = Form("TOT BcOutSdcIn L%2d [Track]", i);
    TString title71 = Form("Residual BcOutSdcIn L%2d (0<theta<15)", i);
    TString title72 = Form("Residual BcOutSdcIn L%2d (15<theta<30)", i);
    TString title73 = Form("Residual BcOutSdcIn L%2d (30<theta<45)", i);
    TString title74 = Form("Residual BcOutSdcIn L%2d (45<theta)", i);
    HB1(100*i+11, title11, nwire, 0., nwire);
    HB1(100*i+12, title12, nbindt, mindt, maxdt);
    HB1(100*i+13, title13, nbindl, mindl, maxdl);
    HB1(100*i+14, title14, 100, -250., 250.);
    HB1(100*i+15, title15, 400, -2.0, 2.0);
    HB2(100*i+16, title16, 250, -250., 250., 100, -1.0, 1.0);
    HB2(100*i+17, title17, 100, -250., 250., 100, -250., 250.);
    HB2(100*i+18, title18, 50, 3., 3., 100, -1.0, 1.0);
    HB2(100*i+19, title19, nbindt, mindt, maxdt,
        Int_t(maxdl*20), -maxdl, maxdl);
    HBProf(100*i+20, title20, 100, -5., 50., 0., 12.);
    HB2(100*i+22, title22, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB1(100*i+21, title21, 200, -5.0, 5.0);
    HB1(100*i+40, title40, 300,    0, 300);
    HB1(100*i+71, title71, 200, -5.0, 5.0);
    HB1(100*i+72, title72, 200, -5.0, 5.0);
    HB1(100*i+73, title73, 200, -5.0, 5.0);
    HB1(100*i+74, title74, 200, -5.0, 5.0);
  }

  // Tracking Histgrams
  HB1(10, "#Tracks BcOutSdcIn", 10, 0., 10.);
  HB1(11, "#Hits of Track BcOutSdcIn", 24, 0., 24.);
  HB1(12, "Chisqr BcOutSdcIn", 500, 0., 50.);
  HB1(13, "LayerId BcOutSdcIn", 24, 0., 24);
  HB1(14, "X0 BcOutSdcIn", 400, -100., 100.);
  HB1(15, "Y0 BcOutSdcIn", 400, -100., 100.);
  HB1(16, "U0 BcOutSdcIn", 200, -0.20, 0.20);
  HB1(17, "V0 BcOutSdcIn", 200, -0.20, 0.20);
  HB2(18, "U0%X0 BcOutSdcIn", 100, -100., 100., 100, -0.20, 0.20);
  HB2(19, "V0%Y0 BcOutSdcIn", 100, -100., 100., 100, -0.20, 0.20);
  HB2(20, "X0%Y0 BcOutSdcIn", 100, -100., 100., 100, -100, 100);

  HB1(21, "Xtgt BcOutSdcIn", 400, -100., 100.);
  HB1(22, "Ytgt BcOutSdcIn", 400, -100., 100.);
  HB1(23, "Utgt BcOutSdcIn", 200, -0.20, 0.20);
  HB1(24, "Vtgt BcOutSdcIn", 200, -0.20, 0.20);
  HB2(25, "Utgt%Xtgt BcOutSdcIn", 100, -100., 100., 100, -0.20, 0.20);
  HB2(26, "Vtgt%Ytgt BcOutSdcIn", 100, -100., 100., 100, -0.20, 0.20);
  HB2(27, "Xtgt%Ytgt BcOutSdcIn", 100, -100., 100., 100, -100, 100);

  //HB1(38, "Plane Eff", 36, 0, 36);
  /*
  HB2(51, "X-X' 245 BcOutSdcIn", 400, -100., 100., 120, -60, 60);
  HB2(52, "X-X' 600 BcOutSdcIn", 400, -100., 100., 120, -60, 60);
  HB2(53, "X-X' 1200 BcOutSdcIn", 400, -100., 100., 120, -60, 60);
  HB2(54, "X-X' 1600 BcOutSdcIn", 400, -100., 100., 120, -60, 60);

  HB2(61, "X-X 600 BcOutSdcIn", 400, -100., 100., 400, -100, 100);
  HB2(62, "X-X 1200 BcOutSdcIn", 400, -100., 100., 400, -100, 100);
  HB2(63, "X-X 1600 BcOutSdcIn", 400, -100., 100., 400, -100, 100);
  */
  // Analysis status
  HB1(40, "Tacking status", 11, -1., 10.);
  HB2(41, "BC3X0/BC4X1", 20, 0, 20, 20, 0, 20);

  ////////////////////////////////////////////
  //Tree
  HBTree("bcsdc", "tree of BcOutSdcInSdcInTracking");
  tree->Branch("evnum", &event.evnum, "evnum/I");
  tree->Branch("trigpat", event.trigpat, Form("trigpat[%d]/I", MaxHits));
  tree->Branch("trigflag", event.trigflag, Form("trigflag[%d]/I", NumOfSegTrig));

  //Hodoscope
  tree->Branch("nhBh1",    &event.nhBh1,   "nhBh1/I");
  tree->Branch("tBh1",      event.tBh1,    Form("tBh1[%d]/D",   MaxHits));
  tree->Branch("deBh1",     event.deBh1,   Form("deBh1[%d]/D",  MaxHits));

  tree->Branch("nhBh2",    &event.nhBh2,   "nhBh2/I");
  tree->Branch("tBh2",      event.tBh2,    Form("tBh2[%d]/D",   MaxHits));
  tree->Branch("deBh2",     event.deBh2,   Form("deBh2[%d]/D",  MaxHits));
  tree->Branch("Bh2Seg",    event.Bh2Seg,  Form("Bh2Seg[%d]/D", MaxHits));

  tree->Branch("btof",     &event.btof,      "btof/D");
  tree->Branch("Time0Seg", &event.Time0Seg,  "Time0Seg/D");
  tree->Branch("deTime0",  &event.deTime0,   "deTime0/D");
  tree->Branch("Time0",    &event.Time0,     "Time0/D");
  tree->Branch("CTime0",   &event.CTime0,    "CTime0/D");

  tree->Branch("ntrack", &event.ntrack, "ntrack/I");
  tree->Branch("nh",     &event.nh, "nh[ntrack]/I");
  tree->Branch("chisqr",  event.chisqr, "chisqr[ntrack]/D");
  tree->Branch("x0",      event.x0,     "x0[ntrack]/D");
  tree->Branch("y0",      event.y0,     "y0[ntrack]/D");
  tree->Branch("u0",      event.u0,     "u0[ntrack]/D");
  tree->Branch("v0",      event.v0,     "v0[ntrack]/D");

  tree->Branch("ntBcOut",    &event.ntBcOut,     "ntBcOut/I");
  tree->Branch("nhBcOut",     event.nhBcOut,     "nhBcOut[ntBcOut]/I");
  tree->Branch("chisqrBcOut", event.chisqrBcOut, "chisqrBcOut[ntBcOut]/D");
  tree->Branch("x0BcOut",     event.x0BcOut,     "x0BcOut[ntBcOut]/D");
  tree->Branch("y0BcOut",     event.y0BcOut,     "y0BcOut[ntBcOut]/D");
  tree->Branch("u0BcOut",     event.u0BcOut,     "u0BcOut[ntBcOut]/D");
  tree->Branch("v0BcOut",     event.v0BcOut,     "v0BcOut[ntBcOut]/D");

  tree->Branch("ntSdcIn",    &event.ntSdcIn,     "ntSdcIn/I");
  tree->Branch("nhSdcIn",     event.nhSdcIn,     "nhSdcIn[ntSdcIn]/I");
  tree->Branch("chisqrSdcIn", event.chisqrSdcIn, "chisqrSdcIn[ntSdcIn]/D");
  tree->Branch("x0SdcIn",     event.x0SdcIn,     "x0SdcIn[ntSdcIn]/D");
  tree->Branch("y0SdcIn",     event.y0SdcIn,     "y0SdcIn[ntSdcIn]/D");
  tree->Branch("u0SdcIn",     event.u0SdcIn,     "u0SdcIn[ntSdcIn]/D");
  tree->Branch("v0SdcIn",     event.v0SdcIn,     "v0SdcIn[ntSdcIn]/D");

  // HPrint();

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
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

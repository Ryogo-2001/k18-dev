// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <UnpackerManager.hh>

#include "ConfMan.hh"
#include "DCRawHit.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "RMAnalyzer.hh"
#include "RootHelper.hh"
#include "KuramaLib.hh"
#include "RawData.hh"

#define HodoCut     0
#define TotCut      1
#define Chi2Cut     0
#define MaxMultiCut 0
#define UseTOF      0 // use or not TOF for tracking

namespace
{
using namespace root;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& zTOF = gGeom.LocalZ("TOF");
const auto& zTOFU = gGeom.LocalZ("TOF-UX");
const auto& zTOFD = gGeom.LocalZ("TOF-DX");
}

//_____________________________________________________________________________
class UserSdcOutTracking : public VEvent
{
private:
  RawData*      rawData;
  HodoAnalyzer* hodoAna;
  DCAnalyzer*   DCAna;

public:
  UserSdcOutTracking();
  ~UserSdcOutTracking();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserSdcOutTracking::ClassName()
{
  static TString s_name("UserSdcOutTracking");
  return s_name;
}

//_____________________________________________________________________________
UserSdcOutTracking::UserSdcOutTracking()
  : VEvent(),
    rawData(new RawData),
    hodoAna(new HodoAnalyzer),
    DCAna(new DCAnalyzer)
{
}

//_____________________________________________________________________________
UserSdcOutTracking::~UserSdcOutTracking()
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

  Double_t btof;
  Double_t stof[MaxHits];

  Int_t nhTof;
  Double_t TofSeg[MaxHits];
  Double_t tTof[MaxHits];
  Double_t dtTof[MaxHits];
  Double_t deTof[MaxHits];

  Int_t nhit[NumOfLayersSdcOut+2];
  Int_t nlayer;
  Double_t pos[NumOfLayersSdcOut+2][MaxHits];
  Double_t wire[NumOfLayersSdcOut+2][MaxHits];

  Int_t ntrack;
  Double_t chisqr[MaxHits];
  Double_t x0[MaxHits];
  Double_t y0[MaxHits];
  Double_t u0[MaxHits];
  Double_t v0[MaxHits];
  Double_t xTof[MaxHits];
  Double_t yTof[MaxHits];
  Double_t uTof[MaxHits];
  Double_t vTof[MaxHits];

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum     =  0;
  nlayer    =  0;
  ntrack    =  0;
  nhBh2     =  0;
  nhBh1     =  0;
  nhTof     =  0;

  btof      = qnan;
  Time0Seg  = qnan;
  deTime0   = qnan;
  Time0     = qnan;
  CTime0    = qnan;

  for(Int_t it=0; it<NumOfSegTrig; it++){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<MaxHits; it++){
    tBh1[it]   = qnan;
    deBh1[it]  = qnan;
    Bh2Seg[it] = -1;
    tBh2[it]   = qnan;
    deBh2[it]  = qnan;
    stof[it]   = qnan;
    TofSeg[it] = -1;
    tTof[it]   = qnan;
    dtTof[it]  = qnan;
    deTof[it]  = qnan;
    chisqr[it] = qnan;
    x0[it] = qnan;
    y0[it] = qnan;
    u0[it] = qnan;
    v0[it] = qnan;
    xTof[it] = qnan;
    yTof[it] = qnan;
    uTof[it] = qnan;
    vTof[it] = qnan;
  }

  for(Int_t it=0; it<NumOfLayersSdcOut; ++it){
    nhit[it] = 0;
    for(Int_t that=0; that<MaxHits; ++that){
      pos[it][that] = qnan;
      wire[it][that] = qnan;
    }
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
UserSdcOutTracking::ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserSdcOutTracking::ProcessingNormal()
{

#if HodoCut
  static const auto MinDeBH2 = gUser.GetParameter("DeBH2", 0);
  static const auto MaxDeBH2 = gUser.GetParameter("DeBH2", 1);
  static const auto MinDeBH1 = gUser.GetParameter("DeBH1", 0);
  static const auto MaxDeBH1 = gUser.GetParameter("DeBH1", 1);
  static const auto MinBeamToF = gUser.GetParameter("BTOF", 0);
  static const auto MaxBeamToF = gUser.GetParameter("BTOF", 1);
#endif
  static const auto MinDeTOF = gUser.GetParameter("DeTOF", 0);
  static const auto MaxDeTOF = gUser.GetParameter("DeTOF", 1);
  static const auto MinTimeTOF = gUser.GetParameter("TimeTOF", 0);
  static const auto MaxTimeTOF = gUser.GetParameter("TimeTOF", 1);
  static const auto MinStopTimingSdcOut = gUser.GetParameter("StopTimingSdcOut", 0);
  static const auto MaxStopTimingSdcOut = gUser.GetParameter("StopTimingSdcOut", 1);
  // static const auto StopTimeDiffSdcOut = gUser.GetParameter("StopTimeDiffSdcOut");
#if TotCut
  static const auto MinTotSDC3 = gUser.GetParameter("MinTotSDC3");
  static const auto MinTotSDC4 = gUser.GetParameter("MinTotSDC4");
  static const auto MinTotSDC5 = gUser.GetParameter("MinTotSDC5");
#endif
#if MaxMultiCut
  static const auto MaxMultiHitSdcOut = gUser.GetParameter("MaxMultiHitSdcOut");
#endif

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

  //////////////BH2 time 0
  hodoAna->DecodeBH2Hits(rawData);
  Int_t nhBh2 = hodoAna->GetNHitsBH2();
  event.nhBh2 = nhBh2;
#if HodoCut
  if(nhBh2==0) return true;
#endif
  HF1(1, 2);

  Double_t time0 = -9999.;
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
    if(dE<MinDeBH1 || MaxDeBH1<dE) continue;
    if(btof<MinBeamToF || MaxBeamToF<btof) continue;
#endif
    event.tBh1[i]  = cmt;
    event.deBh1[i] = dE;
  }

  Double_t btof0 = -999.;
  HodoCluster* cl_btof0 = event.Time0Seg > 0 ?
    hodoAna->GetBtof0BH1Cluster(event.CTime0) : nullptr;
  if(cl_btof0) btof0 = cl_btof0->CMeanTime() - time0;
  event.btof = btof0;

  HF1(1, 5.);


  //////////////Tof Analysis
  HodoClusterContainer TOFCont;
  hodoAna->DecodeTOFHits(rawData);
  // hodoAna->TimeCutTOF(7, 25);
  Int_t nhTof = hodoAna->GetNClustersTOF();
#if HodoCut
  if(nhTof!=0) return true;
#endif
  event.nhTof = nhTof;
  {
    Int_t nhOk = 0;
    for(Int_t i=0; i<nhTof; ++i){
      HodoCluster *hit = hodoAna->GetClusterTOF(i);
      if(!hit) continue;
      Double_t cmt  = hit->CMeanTime();
      Double_t dt   = hit->TimeDif();
      Double_t de   = hit->DeltaE();
      Double_t stof = cmt-time0;
      event.TofSeg[i] = hit->MeanSeg()+1;
      event.tTof[i]   = cmt;
      event.dtTof[i]  = dt;
      event.deTof[i]  = de;
      event.stof[i]   = stof;
      TOFCont.push_back(hit);
      if(MinDeTOF<de  && de<MaxDeTOF  &&
         MinTimeTOF<stof && stof<MaxTimeTOF){
        ++nhOk;
      }
    }
#if HodoCut
    if(nhOk==0) return true;
#endif
  }

  // if(event.nhTof==0) return true;

  HF1(1, 6.);

  // Common stop timing
  HF1(41, common_stop_tdc);
  if(nhTof > 0) HF1(42, common_stop_tdc);
  if(nhTof == 0) HF1(43, common_stop_tdc);
  Bool_t common_stop_is_tof = (common_stop_tdc < MinStopTimingSdcOut
                               || MaxStopTimingSdcOut < common_stop_tdc);
 // if(!common_stop_is_tof) return true;

  HF1(1, 7);

  //////SdcIn
#if 0
  static const Double_t MaxMultiHitSdcIn = gUser.GetParameter("MaxMultiHitSdcIn");
  //////////////SdcIn number of hit layer
  DCAna->DecodeSdcInHits(rawData);
  Double_t multi_SdcIn=0.;
  for(Int_t layer=1; layer<=NumOfLayersSdcIn; ++layer){
    Int_t nhIn = DCAna->GetSdcInHC(layer).size();
    multi_SdcIn += Double_t(nhIn);
  }
# if MaxMultiCut
  if(multi_SdcIn/Double_t(NumOfLayersSdcIn) > MaxMultiHitSdcIn)
    return true;
# endif
#endif

  HF1(1, 10.);
  // Double_t offset = common_stop_is_tof ? 0 : StopTimeDiffSdcOut;
  DCAna->DecodeSdcOutHits(rawData);

#if TotCut
  //DCAna->TotCutSDC3(MinTotSDC3);
  //DCAna->TotCutSDC3(150);
  //DCAna->TotCutSDC4(MinTotSDC4);
  //DCAna->TotCutSDC5(MinTotSDC5);
#endif
#if 0
  const auto& contOut_3y = DCAna->GetSdcOutHC(3);
  const auto& contOut_5y = DCAna->GetSdcOutHC(9);
  if(contOut_3y.size()==0 || contOut_5y.size()==0){
    return true;
  }
  const auto& contOut_3x = DCAna->GetSdcOutHC(1);
  const auto& contOut_5x = DCAna->GetSdcOutHC(11);
  if(contOut_3x.size()==0 || contOut_5x.size()==0){
    return true;
  }
  DCHit *hit_3y=contOut_3y[0];
  Double_t wire_3y=hit_3y->GetWire();
  DCHit *hit_5y=contOut_5y[0];
  Double_t wire_5y=hit_5y->GetWire();
  if(wire_3y < wire_5y-27 || wire_3y > wire_5y+59){
    return true;
  }
  DCHit *hit_3x=contOut_3x[0];
  Double_t wire_3x=hit_3x->GetWire();
  DCHit *hit_5x=contOut_5x[0];
  Double_t wire_5x=hit_5x->GetWire();
  if(wire_3x < wire_5x-43 || wire_3x > wire_5x+43){
    return true;
  }
#endif
  Double_t multi_SdcOut = 0.;
  for(Int_t layer=1; layer<=NumOfLayersSdcOut; ++layer) {
    const auto& contOut = DCAna->GetSdcOutHC(layer);
    Int_t nhOut = contOut.size();
    event.nhit[layer-1] = nhOut;
    if(nhOut>0) event.nlayer++;
    multi_SdcOut += Double_t(nhOut);
    HF1(100*layer, nhOut);
    Int_t plane_eff = (layer-1)*3;
    Bool_t fl_valid_sig = false;
    Int_t tdc1st_2 = -1;
    for(Int_t i=0; i<nhOut; ++i){
      DCHit *hit=contOut[i];
      Double_t wire=hit->GetWire();
      HF1(100*layer+1, wire-0.5);
      event.wire[layer-1][i] = wire-0.5;
      Int_t nhtdc = hit->GetTdcSize();
      Int_t tdc1st = -1;
      for(Int_t k=0; k<nhtdc; k++){
        Int_t tdc = hit->GetTdcVal(k);
        Int_t tot_1 = hit->GetTot(k);
        HF1(100*layer+2, tdc);
        HF1(10000*layer+Int_t(wire), tdc);
        HF2(100*layer+51, tdc, wire-0.5);
        HF2(100*layer+53, tdc, tot_1);
        if(tdc > tdc1st){
          tdc1st = tdc;
          fl_valid_sig = true;
        }
      }
      HF1(100*layer+6, tdc1st);
      HF2(100*layer+52, tdc1st, wire-0.5);
      if(tdc1st > tdc1st_2){
        tdc1st_2 = tdc1st;
      }
      for(Int_t k=0, n=hit->GetTdcTrailingSize(); k<n; ++k){
        Int_t trailing = hit->GetTdcTrailing(k);
        HF1(100*layer+10, trailing);
      }

      if(i<MaxHits){
        event.pos[layer-1][i] = hit->GetWirePosition();
      }

      Int_t nhdt = hit->GetDriftTimeSize();
      Int_t tot1st = -1;
      for(Int_t k=0; k<nhdt; k++){
        Double_t dt = hit->GetDriftTime(k);
        if(common_stop_is_tof) HF1(100*layer+3, dt);
        else                   HF1(100*layer+8, dt);
        HF1(10000*layer+1000+Int_t(wire), dt);

        Double_t tot = hit->GetTot(k);
        HF1(100*layer+5, tot);
        if(tot > tot1st){
          tot1st = tot;
        }
      }
      HF1(100*layer+7, tot1st);
      Int_t nhdl = hit->GetDriftTimeSize();
      for(Int_t k=0; k<nhdl; k++){
        Double_t dl = hit->GetDriftLength(k);
        Double_t tot_2 = hit->GetTot(k);
        HF1(100*layer+4, dl);
        HF2(100*layer+54, dl, tot_2);
      }
    }
    HF1(100*layer+9, tdc1st_2);
    if(fl_valid_sig) ++plane_eff;
    HF1(38, plane_eff);
  }

#if MaxMultiCut
  if(multi_SdcOut/Double_t(NumOfLayersSdcOut) > MaxMultiHitSdcOut)
    return true;
#endif

  HF1(1, 11.);

  // std::cout << "==========TrackSearch SdcOut============" << std::endl;
#if UseTOF
  DCAna->TrackSearchSdcOut(TOFCont);
#else
  DCAna->TrackSearchSdcOut();
#endif

#if 1
#if Chi2Cut
  DCAna->ChiSqrCutSdcOut(30.);
#endif

  Int_t nt=DCAna->GetNtracksSdcOut();
  // std::cout << nt << std::endl;
  if(MaxHits<nt){
    std::cout << "#W " << FUNC_NAME << " "
              << "too many ntSdcOut " << nt << "/" << MaxHits << std::endl;
    nt = MaxHits;
  }
  event.ntrack=nt;
  HF1(10, Double_t(nt));
  for(Int_t it=0; it<nt; ++it){
    const auto track = DCAna->GetTrackSdcOut(it);
    if(!track) continue;
    Int_t nh=track->GetNHit();
    Double_t chisqr    = track->GetChiSquare();
    Double_t chisqr1st = track->GetChiSquare1st();
    Double_t x0=track->GetX0(), y0=track->GetY0();
    Double_t u0=track->GetU0(), v0=track->GetV0();
    Double_t theta = track->GetTheta();
    Int_t    nitr  = track->GetNIteration();
    event.chisqr[it]=chisqr;
    event.x0[it]=x0;
    event.y0[it]=y0;
    event.u0[it]=u0;
    event.v0[it]=v0;

    HF1(11, Double_t(nh));
    HF1(12, chisqr);
    HF1(14, x0); HF1(15, y0);
    HF1(16, u0); HF1(17, v0);
    HF2(18, x0, u0); HF2(19, y0, v0);
    HF2(20, x0, y0);

    HF1(28, nitr);
    HF1(29, chisqr1st);
    HF1(30, chisqr1st-chisqr);
    if(theta<=10.)
      HF1(31, chisqr1st-chisqr);
    if(10.<theta && theta<=20.)
      HF1(32, chisqr1st-chisqr);
    if(20.<theta && theta<=30.)
      HF1(33, chisqr1st-chisqr);
    if(30.<theta && theta<=40.)
      HF1(34, chisqr1st-chisqr);
    if(40.<theta)
      HF1(35, chisqr1st-chisqr);

    ///// TOF
#if UseTOF
#else
    Double_t z_tof = qnan;
    for(Int_t itof=0; itof<event.nhTof; ++itof){
      Int_t lnum = 0;
      TVector3 gposTof;
      if((Int_t)event.TofSeg[itof]%2 == 0){
        lnum = gGeom.GetDetectorId("TOF-UX");
        gposTof = gGeom.GetGlobalPosition("TOF-UX");
        z_tof = zTOFU;
      }
      if((Int_t)event.TofSeg[itof]%2 == 1){
        lnum = gGeom.GetDetectorId("TOF-DX");
        gposTof = gGeom.GetGlobalPosition("TOF-DX");
        z_tof = zTOFD;
      }
      Double_t wpos = gGeom.CalcWirePosition(lnum, event.TofSeg[itof]);
      TVector3 w(wpos, 0, 0);
      Double_t ytTof = event.dtTof[itof]*77.3511;
      TVector3 posTof = gposTof + w + TVector3(0, ytTof, 0);
      Double_t xtof=track->GetX(z_tof), ytof=track->GetY(z_tof);
      Double_t utof=u0, vtof=v0;
      if(nhTof == 1){
        event.xTof[it] = xtof;
        event.yTof[it] = ytof;
        event.uTof[it] = utof;
        event.vTof[it] = vtof;
        HF2(51, event.TofSeg[itof], xtof);
        HF2(52, posTof.X(), xtof);
        HF1(53, posTof.X() - xtof);
        HF2(54, posTof.X(), posTof.X() - xtof);
        HF2(55, event.dtTof[itof], ytof);
        HF1(56, ytTof - ytof);
        HF2(57, event.dtTof[itof], ytTof - ytof);
      }
      HF1(21, xtof); HF1(22, ytof);
      HF1(23, utof); HF1(24, vtof);
      HF2(25, xtof, utof); HF2(26, ytof, vtof);
      HF2(27, xtof, ytof);
    }
#endif

    for(const auto& hit: track->GetHitArray()){
      if(!hit) continue;
      Int_t layerId = 0;
      layerId = hit->GetLayer();
      if(layerId <= PlMaxSdcOut)
        layerId -= PlOffsSdcOut;
      else
        layerId -= PlOffsTOF - NumOfLayersSdcOut;

      HF1(13, hit->GetLayer());
      HF1(36, Double_t(nh));
      HF1(37, chisqr);

      Double_t wire=hit->GetWire();
      Double_t dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1(100*layerId+11, wire-0.5);
      HF1(100*layerId+12, dt);
      HF1(100*layerId+13, dl);
      Double_t xcal=hit->GetXcal(), ycal=hit->GetYcal();
      Double_t pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1(100*layerId+14, pos);
      HF1(100*layerId+15, res);
      HF2(100*layerId+16, pos, res);
      HF2(100*layerId+17, xcal, ycal);
      //      HF1(100000*layerId+50000+wire, res);
      Double_t wp=hit->GetWirePosition();
      Double_t sign=1.;
      if(pos-wp<0.) sign=-1;
      HF2(100*layerId+18, sign*dl, res);
      Double_t xlcal=hit->GetLocalCalPos();
      HF2(100*layerId+19, dt, xlcal-wp);
      HF1(100*layerId+23, std::abs(dl)-std::abs(pos-wp));
      HF1(100*layerId+24, std::abs(dl)-std::abs(xlcal-wp));

      HF2(100*layerId+31, xcal, res);
      HF2(100*layerId+32, ycal, res);
      HF2(100*layerId+33, u0, res);
      HF2(100*layerId+34, v0, res);

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

      if (std::abs(dl-std::abs(xlcal-wp))<2.0){
        HFProf(100*layerId+20, dt, std::abs(xlcal-wp));
        HF2(100*layerId+22, dt, std::abs(xlcal-wp));
        HFProf(100000*layerId+3000+Int_t(wire), xlcal-wp,dt);
        HF2(100000*layerId+4000+Int_t(wire), xlcal-wp,dt);
      }
    }
  }
#endif

  HF1(1, 12.);

  return true;
}

//_____________________________________________________________________________
Bool_t
UserSdcOutTracking::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserSdcOutTracking;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  const Int_t    NbinSdcOutTdc = 2000;
  const Double_t MinSdcOutTdc  =    0.;
  const Double_t MaxSdcOutTdc  = 2000.;

  const Int_t    NbinSDC3DT = 360;
  const Double_t MinSDC3DT  = -50.;
  const Double_t MaxSDC3DT  = 250.;
  const Int_t    NbinSDC3DL =  90;
  const Double_t MinSDC3DL  =  -2.;
  const Double_t MaxSDC3DL  =   7.;

  const Int_t    NbinSDC4DT = 360;
  const Double_t MinSDC4DT  = -50.;
  const Double_t MaxSDC4DT  = 250.;
  const Int_t    NbinSDC4DL =  90;
  const Double_t MinSDC4DL  =  -2.;
  const Double_t MaxSDC4DL  =   7.;

  const Int_t    NbinSDC5DT = 360;
  const Double_t MinSDC5DT  = -50.;
  const Double_t MaxSDC5DT  = 250.;
  const Int_t    NbinSDC5DL =  90;
  const Double_t MinSDC5DL  =  -2.;
  const Double_t MaxSDC5DL  =   7.;

  HB1(1, "Status", 20, 0., 20.);

  for(Int_t i=1; i<=NumOfLayersSdcOut+NumOfLayersTOF; ++i){
    TString tag;
    Int_t nwire = 0, nbindt = 1, nbindl = 1;
    Double_t mindt = 0., maxdt = 1., mindl = 0., maxdl = 1.;
    if(i<=NumOfLayersSDC3){
      tag    = "SDC3";
      nwire  = MaxWireSDC3;
      nbindt = NbinSDC3DT;
      mindt  = MinSDC3DT;
      maxdt  = MaxSDC3DT;
      nbindl = NbinSDC3DL;
      mindl  = MinSDC3DL;
      maxdl  = MaxSDC3DL;
    }else if(i<=NumOfLayersSDC3+NumOfLayersSDC4){
      tag = "SDC4";
      nwire  = MaxWireSDC4;
      nbindt = NbinSDC4DT;
      mindt  = MinSDC4DT;
      maxdt  = MaxSDC4DT;
      nbindl = NbinSDC4DL;
      mindl  = MinSDC4DL;
      maxdl  = MaxSDC4DL;
    }else if(i<=NumOfLayersSdcOut){
      tag = "SDC5";
      nwire   = (i==9 || i==10) ? MaxWireSDC5Y : MaxWireSDC5X;
      nbindt = NbinSDC5DT;
      mindt  = MinSDC5DT;
      maxdt  = MaxSDC5DT;
      nbindl = NbinSDC5DL;
      mindl  = MinSDC5DL;
      maxdl  = MaxSDC5DL;
    }else if(i<=NumOfLayersSdcOut+NumOfLayersTOF){
      tag = "TOF";
      nwire = NumOfSegTOF;
      nbindt = NbinSDC4DT;
      mindt  = MinSDC4DT;
      maxdt  = MaxSDC4DT;
      nbindl = NbinSDC4DL;
      mindl  = MinSDC4DL;
      maxdl  = MaxSDC4DL;
    }

    if(i<=NumOfLayersSdcOut){
      TString title0 = Form("#Hits %s#%2d", tag.Data(), i);
      TString title1 = Form("Hitpat %s#%2d", tag.Data(), i);
      TString title2 = Form("Tdc %s#%2d", tag.Data(), i);
      TString title3 = Form("Drift Time %s#%2d", tag.Data(), i);
      TString title4 = Form("Drift Length %s#%2d", tag.Data(), i);
      TString title5 = Form("TOT %s#%2d", tag.Data(), i);
      TString title6 = Form("Tdc 1st %s#%2d", tag.Data(), i);
      TString title7 = Form("TOT 1st %s#%2d", tag.Data(), i);
      TString title8 = Form("Drift Time %s#%2d (BH2 timing)", tag.Data(), i);
      TString title9 = Form("Tdc 1st-2 %s#%2d", tag.Data(), i);
      TString title10 = Form("Trailing %s#%2d", tag.Data(), i);
      TString title51 = Form("Tdc vs wire %s#%2d", tag.Data(), i);
      TString title52 = Form("Tdc 1st vs wire %s#%2d", tag.Data(), i);
      TString title53 = Form("Tdc vs tot %s#%2d", tag.Data(), i);
      TString title54 = Form("Drift Length vs tot %s#%2d", tag.Data(), i);
      HB1(100*i+0, title0, nwire+1, 0., nwire+1);
      HB1(100*i+1, title1, nwire, 0., nwire);
      HB1(100*i+2, title2, NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc);
      HB1(100*i+3, title3, nbindt, mindt, maxdt);
      HB1(100*i+4, title4, nbindl, mindl, maxdl);
      HB1(100*i+5, title5, 500,    0, 500);
      HB1(100*i+6, title6, NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc);
      HB1(100*i+7, title7, 500,    0, 500);
      HB1(100*i+8, title8, nbindt, mindt, maxdt);
      HB1(100*i+9, title9, NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc);
      HB1(100*i+10, title10, NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc);
      HB2(100*i+51, title51, NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc, nwire, 0., nwire);
      HB2(100*i+52, title52, NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc, nwire, 0., nwire);
      HB2(100*i+53, title53, NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc, 500, 0., 500);
      HB2(100*i+54, title54, nbindl, mindl, maxdl, 500, 0., 500);
      for (Int_t wire=1; wire<=nwire; wire++){
        TString title11 = Form("Tdc %s#%2d  Wire#%4d", tag.Data(), i, wire);
        TString title12 = Form("DriftTime %s#%2d Wire#%4d", tag.Data(), i, wire);
        TString title13 = Form("DriftLength %s#%2d Wire#%d", tag.Data(), i, wire);
        TString title14 = Form("DriftTime %s#%2d Wire#%4d [Track]", tag.Data(), i, wire);
        HB1(10000*i+wire, title11, NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc);
        HB1(10000*i+1000+wire, title12, nbindt, mindt, maxdt);
        HB1(10000*i+2000+wire, title13, nbindl, mindl, maxdl);
        HB1(10000*i+5000+wire, title14, nbindt, mindt, maxdt);
      }
    }

    // Tracking Histgrams
    TString title11 = Form("HitPat SdcOut%2d [Track]", i);
    TString title12 = Form("DriftTime SdcOut%2d [Track]", i);
    TString title13 = Form("DriftLength SdcOut%2d [Track]", i);
    TString title14 = Form("Position SdcOut%2d", i);
    TString title15 = Form("Residual SdcOut%2d", i);
    TString title16 = Form("Resid%%Pos SdcOut%2d", i);
    TString title17 = Form("Y%%Xcal SdcOut%2d", i);
    TString title18 = Form("Res%%DL SdcOut%2d", i);
    TString title19 = Form("HitPos%%DriftTime SdcOut%2d", i);
    TString title20 = Form("DriftLength%%DriftTime SdcOut%2d", i);
    TString title21 = title15 + " [w/o Self]";
    TString title22 = title20;
    TString title23 = Form("DriftLength-Hit Position SdcOut%2d [Track]", i);
    TString title24 = Form("DriftLength-Cal Position SdcOut%2d [Track]", i);
    TString title40 = Form("TOT SdcOut%2d [Track]", i);
    TString title71 = Form("Residual SdcOut%2d (0<theta<15)", i);
    TString title72 = Form("Residual SdcOut%2d (15<theta<30)", i);
    TString title73 = Form("Residual SdcOut%2d (30<theta<45)", i);
    TString title74 = Form("Residual SdcOut%2d (45<theta)", i);
    HB1(100*i+11, title11, nwire, 0., nwire);
    HB1(100*i+12, title12, nbindt, mindt, maxdt);
    HB1(100*i+13, title13, 100, -5, maxdl);
    HB1(100*i+14, title14, 100, -1000., 1000.);
    if(i<=NumOfLayersSdcOut)
      HB1(100*i+15, title15, 1000, -5.0, 5.0);
    else
      HB1(100*i+15, title15, 1000, -200.0, 200.0);
    if(i<=NumOfLayersSdcOut)
      HB2(100*i+16, title16, 400, -1000., 1000., 100, -1.0, 1.0);
    else
      HB2(100*i+16, title16, 100, -1000., 1000., 100, -200.0, 200.0);
    HB2(100*i+17, title17, 100, -1000., 1000., 100, -1000., 1000.);
    //if(i<=NumOfLayersSDC3)
    if(i<=NumOfLayersSdcOut)
      HB2(100*i+18, title18, 110, -5.5, 5.5, 100, -1.0, 1.0);
    else
      HB2(100*i+18, title18, 110, -11., 11., 100, -1.0, 1.0);
    HB2(100*i+19, title19, nbindt, mindt, maxdt, Int_t(maxdl*20), -maxdl, maxdl);
    HBProf(100*i+20, title20, 100, -50, 300, 0, maxdl);
    HB2(100*i+22, title22, nbindt, mindt, maxdt, 100, -0.5, maxdl);
    HB1(100*i+21, title21, 200, -5.0, 5.0);
    HB1(100*i+23, title23, 1000, -5.0, 5.0);
    HB1(100*i+24, title24, 1000, -5.0, 5.0);
    HB2(100*i+31, Form("Resid%%X SdcOut %d", i), 100, -1000., 1000., 100, -2., 2.);
    HB2(100*i+32, Form("Resid%%Y SdcOut %d", i), 100, -1000., 1000., 100, -2., 2.);
    HB2(100*i+33, Form("Resid%%U SdcOut %d", i), 100, -0.5, 0.5, 100, -2., 2.);
    HB2(100*i+34, Form("Resid%%V SdcOut %d", i), 100, -0.5, 0.5, 100, -2., 2.);
    HB1(100*i+40, title40, 500,    0, 500);
    HB1(100*i+71, title71, 200, -5.0, 5.0);
    HB1(100*i+72, title72, 200, -5.0, 5.0);
    HB1(100*i+73, title73, 200, -5.0, 5.0);
    HB1(100*i+74, title74, 200, -5.0, 5.0);

    for (Int_t j=1; j<=nwire; j++) {
      TString title = Form("XT of Layer %2d Wire #%4d", i, j);
      HBProf(100000*i+3000+j, title, 100, -12., 12., -30, 300);
      HB2(100000*i+4000+j, title, 100, -12., 12., 100, -30., 300.);
    }

  }

  // Tracking Histgrams
  HB1(10, "#Tracks SdcOut", 10, 0., 10.);
  HB1(11, "#Hits of Track SdcOut", 20, 0., 20.);
  HB1(12, "Chisqr SdcOut", 500, 0., 50.);
  HB1(13, "LayerId SdcOut", 60, 30., 90.);
  HB1(14, "X0 SdcOut", 1400, -1200., 1200.);
  HB1(15, "Y0 SdcOut", 1000, -500., 500.);
  HB1(16, "U0 SdcOut", 200, -0.35, 0.35);
  HB1(17, "V0 SdcOut", 200, -0.20, 0.20);
  HB2(18, "U0%X0 SdcOut", 120, -1200., 1200., 100, -0.40, 0.40);
  HB2(19, "V0%Y0 SdcOut", 100, -500., 500., 100, -0.20, 0.20);
  HB2(20, "X0%Y0 SdcOut", 100, -1200., 1200., 100, -500, 500);
  HB1(21, "Xtof SdcOut", 1400, -1200., 1200.);
  HB1(22, "Ytof SdcOut", 1000, -500., 500.);
  HB1(23, "Utof SdcOut", 200, -0.35, 0.35);
  HB1(24, "Vtof SdcOut", 200, -0.20, 0.20);
  HB2(25, "Utof%Xtof SdcOut", 100, -1200., 1200., 100, -0.40, 0.40);
  HB2(26, "Vtof%Ytof SdcOut", 100, -500., 500., 100, -0.20, 0.20);
  HB2(27, "Xtof%Ytof SdcOut", 100, -1200., 1200., 100, -500, 500);

  HB1(28, "NIteration SdcOut", 100, 0., 100.);
  HB1(29, "Chisqr1st SdcOut", 500, 0., 50.);
  HB1(30, "Chisqr1st-Chisqr SdcOut", 500, 0., 10.);
  HB1(31, "Chisqr1st-Chisqr SdcOut (0<theta<10)", 500, 0., 10.);
  HB1(32, "Chisqr1st-Chisqr SdcOut (10<theta<20)", 500, 0., 10.);
  HB1(33, "Chisqr1st-Chisqr SdcOut (20<theta<30)", 500, 0., 10.);
  HB1(34, "Chisqr1st-Chisqr SdcOut (30<theta<40)", 500, 0., 10.);
  HB1(35, "Chisqr1st-Chisqr SdcOut (40<theta)", 500, 0., 10.);
  HB1(36, "#Hits of Track SdcOut(SDC)", 20, 0., 20.);
  HB1(37, "Chisqr SdcOut(SDC)", 500, 0., 50.);

  // Plane eff
  HB1(38, "Plane Eff", 30, 0, 30);

  // Common Stop Timing
  HB1(41, "Common Stop TDC", 0x1000, 0, 0x1000);
  HB1(42, "Common Stop TDC (w/TOF)", 0x1000, 0, 0x1000);
  HB1(43, "Common Stop TDC (w/oTOF)", 0x1000, 0, 0x1000);

  // ySdcOut vs dtTof
  HB2(51, "xSdcOut % TofSeg", NumOfSegTOF, 1, NumOfSegTOF+1, 300, -1500, 1500.);
  HB2(52, "xSdcOut % TofPos", 40, -1500, 1500, 300, -1500, 1500.);
  HB1(53, "TofPos - xSdcOut", 400, -400, 400);
  HB2(54, "TofPos - xSdcOut % TofPos", 40, -1500, 1500., 400, -400, 400);
  HB2(55, "ySdcOut % dtTof", 300, -15., 15., 300, -1500, 1500.);
  HB1(56, "ytTof - ySdcOut", 300, -300., 300.);
  HB2(57, "ytTof - ySdcOut % dtTof", 300, -15., 15., 300, -300, 300.);

  ////////////////////////////////////////////
  //Tree
  HBTree("sdcout", "tree of SdcOutTracking");
  tree->Branch("evnum", &event.evnum, "evnum/I");
  tree->Branch("trigpat", event.trigpat, Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag", event.trigflag, Form("trigflag[%d]/I", NumOfSegTrig));

  tree->Branch("nhBh1", &event.nhBh1, "nhBh1/I");
  tree->Branch("tBh1", event.tBh1, Form("tBh1[%d]/D", MaxHits));
  tree->Branch("deBh1", event.deBh1, Form("deBh1[%d]/D", MaxHits));

  tree->Branch("nhBh2",    &event.nhBh2,   "nhBh2/I");
  tree->Branch("tBh2",      event.tBh2,    Form("tBh2[%d]/D",   MaxHits));
  tree->Branch("deBh2",     event.deBh2,   Form("deBh2[%d]/D",  MaxHits));
  tree->Branch("Bh2Seg",    event.Bh2Seg,  Form("Bh2Seg[%d]/D", MaxHits));

  tree->Branch("Time0Seg", &event.Time0Seg,  "Time0Seg/D");
  tree->Branch("deTime0",  &event.deTime0,   "deTime0/D");
  tree->Branch("Time0",    &event.Time0,     "Time0/D");
  tree->Branch("CTime0",   &event.CTime0,    "CTime0/D");

  tree->Branch("btof",     &event.btof,     "btof/D");
  tree->Branch("stof",     event.stof,    Form("stof[%d]/D", MaxHits));

  tree->Branch("nhTof",   &event.nhTof,   "nhTof/I");
  tree->Branch("TofSeg",   event.TofSeg,  "TofSeg[nhTof]/D");
  tree->Branch("tTof",     event.tTof,    "tTof[nhTof]/D");
  tree->Branch("dtTof",    event.dtTof,   "dtTof[nhTof]/D");
  tree->Branch("deTof",    event.deTof,   "deTof[nhTof]/D");

  tree->Branch("nhit",     &event.nhit,     Form("nhit[%d]/I", NumOfLayersSdcOut));
  tree->Branch("nlayer",   &event.nlayer,   "nlayer/I");
  tree->Branch("pos",      &event.pos,     Form("pos[%d][%d]/D", NumOfLayersSdcOut, MaxHits));
  tree->Branch("wire",      &event.wire,     Form("wire[%d][%d]/D", NumOfLayersSdcOut, MaxHits));
  tree->Branch("ntrack",   &event.ntrack,   "ntrack/I");
  tree->Branch("chisqr",    event.chisqr,   "chisqr[ntrack]/D");
  tree->Branch("x0",        event.x0,       "x0[ntrack]/D");
  tree->Branch("y0",        event.y0,       "y0[ntrack]/D");
  tree->Branch("u0",        event.u0,       "u0[ntrack]/D");
  tree->Branch("v0",        event.v0,       "v0[ntrack]/D");
  tree->Branch("xTof",      event.xTof,     "xTof[ntrack]/D");
  tree->Branch("yTof",      event.yTof,     "yTof[ntrack]/D");
  tree->Branch("uTof",      event.uTof,     "uTof[ntrack]/D");
  tree->Branch("vTof",      event.vTof,     "vTof[ntrack]/D");
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

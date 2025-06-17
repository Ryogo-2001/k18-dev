// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <TSystem.h>

#include "ConfMan.hh"
#include "DCRawHit.hh"
#include "DatabasePDG.hh"
#include "DetectorID.hh"
#include "EventDisplay.hh"
#include "RMAnalyzer.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "FLHit.hh"
#include "FuncName.hh"
#include "KuramaLib.hh"
#include "RawData.hh"
//#include "RootHelper.hh"
#include "UnpackerManager.hh"
#include "BH2Filter.hh"
#include "TPCHit.hh"
#include "TPCPadHelper.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "TPCRawHit.hh"

#define SAVEPDF 0
#define DRAWTEXT_DETECTOR 0
#define DRAWTEXT_TRACKING 1

namespace
{
const auto& gGeom   = DCGeomMan::GetInstance();
auto&       gEvDisp = EventDisplay::GetInstance();
const auto& gUser   = UserParamMan::GetInstance();
auto&       gFilter = BH2Filter::GetInstance();
auto&       gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
const Double_t PionMass   = pdg::PionMass();
const Double_t KaonMass   = pdg::KaonMass();
const Double_t ProtonMass = pdg::ProtonMass();
}

//_____________________________________________________________________________
class UserEventDisplay : public VEvent
{
private:
  RawData*      rawData;
  DCAnalyzer*   DCAna;
  HodoAnalyzer* hodoAna;
public:
  UserEventDisplay();
  ~UserEventDisplay();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserEventDisplay::ClassName()
{
  static TString s_name("UserEventDisplay");
  return s_name;
}

//_____________________________________________________________________________
UserEventDisplay::UserEventDisplay()
  : VEvent(),
    rawData(new RawData),
    DCAna(new DCAnalyzer),
    hodoAna(new HodoAnalyzer)
{
}

//_____________________________________________________________________________
UserEventDisplay::~UserEventDisplay()
{
  if(DCAna) delete DCAna;
  if(hodoAna) delete hodoAna;
  if(rawData) delete rawData;
}

//_____________________________________________________________________________
Bool_t
UserEventDisplay::ProcessingBegin()
{
  return true;
}

//_____________________________________________________________________________
Bool_t
UserEventDisplay::ProcessingNormal()
{
  static const auto MaxMultiHitBcOut = gUser.GetParameter("MaxMultiHitBcOut");
  static const auto MaxMultiHitSdcIn = gUser.GetParameter("MaxMultiHitSdcIn");
  static const auto MaxMultiHitSdcOut = gUser.GetParameter("MaxMultiHitSdcOut");

  static const auto MinTdcBH2 = gUser.GetParameter("TdcBH2", 0);
  static const auto MaxTdcBH2 = gUser.GetParameter("TdcBH2", 1);
  static const auto MinTdcBH1 = gUser.GetParameter("TdcBH1", 0);
  static const auto MaxTdcBH1 = gUser.GetParameter("TdcBH1", 1);
  static const auto MinTimeBFT = gUser.GetParameter("TimeBFT", 0);
  static const auto MaxTimeBFT = gUser.GetParameter("TimeBFT", 1);
  static const auto MinTdcSCH  = gUser.GetParameter("TdcSCH", 0);
  static const auto MaxTdcSCH  = gUser.GetParameter("TdcSCH", 1);
  static const auto MinTdcTOF  = gUser.GetParameter("TdcTOF", 0);
  static const auto MaxTdcTOF  = gUser.GetParameter("TdcTOF", 1);
  static const auto MinTdcHTOF = gUser.GetParameter("TdcHTOF", 0);
  static const auto MaxTdcHTOF = gUser.GetParameter("TdcHTOF", 1);
  static const auto MinTdcWC = gUser.GetParameter("TdcWC", 0);
  static const auto MaxTdcWC = gUser.GetParameter("TdcWC", 1);

  // static const auto StopTimeDiffSdcOut = gUser.GetParameter("StopTimeDiffSdcOut");
  // static const auto MinStopTimingSdcOut = gUser.GetParameter("StopTimingSdcOut", 0);
  // static const auto MaxStopTimingSdcOut = gUser.GetParameter("StopTimingSdcOut", 1);
  static const auto MinTotBcOut = gUser.GetParameter("MinTotBcOut");
  static const auto MinTotSDC3 = gUser.GetParameter("MinTotSDC3");
  static const auto MinTotSDC4 = gUser.GetParameter("MinTotSDC4");

  static const Int_t NumOfTimeBucket = gUser.GetParameter("NumOfTimeBucket");

  // static const Int_t IdBH2 = gGeom.GetDetectorId("BH2");
  static const Int_t IdSCH = gGeom.GetDetectorId("SCH");
  static const Int_t IdTOF = gGeom.GetDetectorId("TOF");
  static const Int_t IdWC = gGeom.GetDetectorId("WC");

  static TString evinfo;
  evinfo = Form("Run# %5d%4sEvent# %6d",
                gUnpacker.get_run_number(), "",
                gUnpacker.get_event_number());
  hddaq::cout << "\033c" << TString('=', 80) << std::endl
              << "[Info] " << evinfo << std::endl;
  gEvDisp.DrawRunEvent(0.04, 0.5, evinfo);
  rawData->DecodeHits();

  //________________________________________________________
  //___ TrigRawHit
  std::bitset<NumOfSegTrig> trigger_flag;
  for(auto& hit: rawData->GetTrigRawHC()){
    if(hit->GetTdc1() > 0) trigger_flag.set(hit->SegmentId());
  }
  if(trigger_flag[trigger::kSpillEnd]) return true;
  // if(!trigger_flag[trigger::kTrigBPS]) return true;
  hddaq::cout << "[Info] TrigPat = " << trigger_flag << std::endl;

  //________________________________________________________
  //___ BH2RawHit
  std::vector<Int_t> BH2SegCont;
  for(const auto& hit: rawData->GetBH2RawHC()){
    Int_t seg = hit->SegmentId();
    Bool_t is_hit_u = false;
    Bool_t is_hit_d = false;
    for(Int_t j=0, m=hit->GetSizeTdcUp(); j<m; ++j){
      Int_t Tu = hit->GetTdcUp(j);
      if(MinTdcBH2 < Tu && Tu < MaxTdcBH2){
        is_hit_u = true;
      }
      // gEvDisp.FillBH2(seg, Tu);
    }
    for(Int_t j=0, m=hit->GetSizeTdcDown(); j<m; ++j){
      Int_t Td = hit->GetTdcDown(j);
      if(MinTdcBH2 < Td && Td < MaxTdcBH2){
        is_hit_d = true;
      }
      // gEvDisp.FillBH2(seg, Td);
    }
    if(is_hit_u && is_hit_d){
      hddaq::cout << "[Info] Bh2Seg = " << seg << std::endl;
      BH2SegCont.push_back(seg);
    }
  }

  //________________________________________________________
  //___ BH1RawHit
  std::vector<Int_t> BH1SegCont;
  for(const auto& hit: rawData->GetBH1RawHC()){
    Int_t seg = hit->SegmentId();
    Bool_t is_hit_u = false;
    Bool_t is_hit_d = false;
    for(Int_t j=0, m=hit->GetSizeTdcUp(); j<m; ++j){
      Int_t Tu = hit->GetTdcUp(j);
      if(MinTdcBH1 < Tu && Tu < MaxTdcBH1){
        is_hit_u = true;
      }
      // gEvDisp.FillBH1(seg, Tu);
    }
    for(Int_t j=0, m=hit->GetSizeTdcDown(); j<m; ++j){
      Int_t Td = hit->GetTdcDown(j);
      if(MinTdcBH1 < Td && Td < MaxTdcBH1){
        is_hit_d = true;
      }
      // gEvDisp.FillBH1(seg, Td);
    }
    if(is_hit_u && is_hit_d){
      hddaq::cout << "[Info] Bh1Seg = " << seg << std::endl;
      BH1SegCont.push_back(seg);
    }
  }

  //________________________________________________________
  //___ BH2HodoCluster
  hodoAna->DecodeBH2Hits(rawData);
  const auto Time0Cl = hodoAna->GetTime0BH2Cluster();
  if(!Time0Cl){
    hddaq::cout << "[Warning] Time0Cl is null!" << std::endl;
    // gEvDisp.GetCommand();
    return true;
  }
  Double_t ctime0 = Time0Cl->CTime0();
  // Double_t t0seg = Time0Cl->MeanSeg();

  //________________________________________________________
  //___ BH1HodoCluster
  hodoAna->DecodeBH1Hits(rawData);
  const auto Btof0Cl = hodoAna->GetBtof0BH1Cluster(ctime0);
  Double_t btof = (Btof0Cl)
    ? Btof0Cl->CMeanTime() - ctime0
    : TMath::QuietNaN();

  //________________________________________________________
  //___ TOFRawHit
  std::vector<Int_t> TOFSegCont;
  for(const auto& hit: rawData->GetTOFRawHC()){
    Int_t seg = hit->SegmentId();
    Bool_t is_hit_u = false;
    Bool_t is_hit_d = false;
    for(Int_t j=0, m=hit->GetSizeTdcUp(); j<m; ++j){
      Int_t Tu = hit->GetTdcUp(j);
      if(MinTdcTOF < Tu && Tu < MaxTdcTOF){
        is_hit_u = true;
      }
      // gEvDisp.FillTOF(seg, Tu);
    }
    for(Int_t j=0, m=hit->GetSizeTdcDown(); j<m; ++j){
      Int_t Td = hit->GetTdcDown(j);
      if(MinTdcTOF < Td && Td < MaxTdcTOF){
        is_hit_d = true;
      }
      // gEvDisp.FillTOF(seg, Td);
    }
    if(is_hit_u || is_hit_d){
      gEvDisp.DrawHitHodoscope(IdTOF, seg, is_hit_u, is_hit_d);
    }
    if(is_hit_u && is_hit_d){
      hddaq::cout << "[Info] TofSeg = " << seg << std::endl;
      TOFSegCont.push_back(seg);
    }
  }

  //________________________________________________________
  //___ TOFHodoHit
  hodoAna->DecodeTOFHits(rawData);
  const auto& TOFCont = hodoAna->GetHitsTOF();
  if(TOFCont.empty()){
    hddaq::cout << "[Warning] TOFCont is empty!" << std::endl;
    //gEvDisp.GetCommand();
    return true;
  }

  //________________________________________________________
  //___ SCHRawHit
  for(const auto& hit: rawData->GetSCHRawHC()){
    Int_t seg = hit->SegmentId();
    Bool_t is_hit = false;
    for(Int_t j=0, m=hit->GetSizeTdcUp(); j<m; ++j){
      Int_t Tu = hit->GetTdcUp(j);
      if(Tu > 0) gEvDisp.FillSCH(seg, Tu);
      if(MinTdcSCH < Tu && Tu < MaxTdcSCH){
        is_hit = true;
      }
    }
    if(is_hit) gEvDisp.DrawHitHodoscope(IdSCH, seg);
  }

  //________________________________________________________
  //___ WCRawHit
  std::vector<Int_t> WCSegCont;
  for(const auto& hit: rawData->GetWCSUMRawHC()){
    Int_t seg = hit->SegmentId();
    Bool_t is_hit = false;
    for(Int_t j=0, m=hit->GetSizeTdcUp(); j<m; ++j){
      Int_t Tu = hit->GetTdcUp(j);
      if(MinTdcWC < Tu && Tu < MaxTdcWC){
        is_hit = true;
      }
    }
    if(is_hit){
      gEvDisp.DrawHitHodoscope(IdWC, seg, is_hit, is_hit);
      hddaq::cout << "[Info] WcSeg = " << seg << std::endl;
      WCSegCont.push_back(seg);
    }
  }

#if 0
  static const Int_t IdSDC1 = gGeom.DetectorId("SDC1-X1");
  // static const Int_t IdSDC2 = gGeom.DetectorId("SDC2-X1");
  static const Int_t IdSDC3 = gGeom.DetectorId("SDC3-X1");
  static const Int_t IdSDC4 = gGeom.DetectorId("SDC4-X1");
  //________________________________________________________
  //___ BcOutRawHit
  for(Int_t layer=1; layer<=NumOfLayersBcOut; ++layer){
    for(const auto& hit: rawData->GetBcOutRawHC(layer)){
      Int_t wire = hit->WireId();
      for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
        Int_t tdc = hit->GetTdc(j);
        if(tdc > 0) gEvDisp.FillBcOutHit(layer, wire, tdc);
      }
    }
  }
  //________________________________________________________
  //___ SdcInRawHit
  for(const auto& hit: rawData->GetSdcInRawHC(IdSDC1)){
    Int_t wire = hit->WireId();
    for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
      Int_t tdc = hit->GetTdc(j);
      if(tdc>0) gEvDisp.FillSDC1(wire, tdc);
    }
  }
  for(const auto& hit: rawData->GetSdcInRawHC(IdSDC1+1)){
    Int_t wire = hit->WireId();
    for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
      Int_t tdc = hit->GetTdc(j);
      if(tdc>0) gEvDisp.FillSDC1p(wire, tdc);
    }
  }
  //________________________________________________________
  //___ SdcOutRawHit
  for(const auto& hit: rawData->GetSdcOutRawHC(IdSDC3-30)){
    Int_t wire = hit->WireId();
    for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
      Int_t tdc = hit->GetTdc(j);
      if(tdc>0) gEvDisp.FillSDC3_Leading(wire, tdc);
    }
    for(Int_t j=0, m=hit->GetTrailingSize(); j<m; ++j){
      Int_t tdc = hit->GetTrailing(j);
      if(tdc>0) gEvDisp.FillSDC3_Trailing(wire, tdc);
    }
  }
  for(const auto& hit: rawData->GetSdcOutRawHC(IdSDC3+1-30)){
    Int_t wire = hit->WireId();
    for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
      Int_t tdc = hit->GetTdc(j);
      if(tdc>0) gEvDisp.FillSDC3p_Leading(wire, tdc);
    }
    for(Int_t j=0, m=hit->GetTrailingSize(); j<m; ++j){
      Int_t tdc = hit->GetTrailing(j);
      if(tdc>0) gEvDisp.FillSDC3p_Trailing(wire, tdc);
    }
  }
  for(const auto& hit: rawData->GetSdcOutRawHC(IdSDC4-30)){
    Int_t wire = hit->WireId();
    for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
      Int_t tdc = hit->GetTdc(j);
      if(tdc>0) gEvDisp.FillSDC4_Leading(wire, tdc);
    }
    for(Int_t j=0, m=hit->GetTrailingSize(); j<m; ++j){
      Int_t tdc = hit->GetTrailing(j);
      if(tdc>0) gEvDisp.FillSDC4_Trailing(wire, tdc);
    }
  }
  for(const auto& hit: rawData->GetSdcOutRawHC(IdSDC4+1-30)){
    Int_t wire = hit->WireId();
    for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
      Int_t tdc = hit->GetTdc(j);
      if(tdc>0) gEvDisp.FillSDC4p_Leading(wire, tdc);
    }
    for(Int_t j=0, m=hit->GetTrailingSize(); j<m; ++j){
      Int_t tdc = hit->GetTrailing(j);
      if(tdc>0) gEvDisp.FillSDC4p_Trailing(wire, tdc);
    }
  }
  //________________________________________________________
  //___ BFTRawHit
  for(Int_t layer=0; layer<NumOfPlaneBFT; layer++){
    for(const auto& hit: rawData->GetBFTRawHC(layer)){
      Int_t seg = hit->SegmentId();
      for(Int_t j=0, m=hit->GetSizeTdcUp(); j<m; ++j){
        Int_t Tu = hit->GetTdcUp(j);
        // hddaq::cout << "[Info] BFTxSeg = " << seg << std::endl;
        if(Tu>0) gEvDisp.FillBFT(layer, seg, Tu);
      }
    }
  }
#endif

  //________________________________________________________
  //___ BcOutDCHit
  DCAna->DecodeBcOutHits(rawData);
  DCAna->TotCutBCOut(MinTotBcOut);
  Double_t multi_BcOut = 0.;
  for(Int_t layer=1; layer<=NumOfLayersBcOut; ++layer){
    const auto& cont = DCAna->GetBcOutHC(layer);
    multi_BcOut += cont.size();
    // for(const auto& hit: cont){
    //   Double_t wire = hit->GetWire();
    //   Bool_t goodFlag = false;
    //   for(Int_t j=0, mh=hit->GetTdcSize(); j<mh; j++){
    //     if(hit->IsWithinRange(j)){
    //       goodFlag = true;
    //       break;
    //     }
    //   }
    //   if(goodFlag){
    //     gEvDisp.DrawHitWire(layer+112, Int_t(wire));
    //   }else{
    //     gEvDisp.DrawHitWire(layer+112, Int_t(wire), false, false);
    //   }
    // }
  }
  multi_BcOut /= (Double_t)NumOfLayersBcOut;
  if(multi_BcOut > MaxMultiHitBcOut){
    hddaq::cout << "[Warning] BcOutHits exceed MaxMultiHit "
                << multi_BcOut << "/" << MaxMultiHitBcOut << std::endl;
    // gEvDisp.GetCommand();
    return true;
  }

  //________________________________________________________
  //___ BcOutTracking
  DCAna->TrackSearchBcOut();
  Int_t ntBcOut = DCAna->GetNtracksBcOut();
  hddaq::cout << "[Info] ntBcOut = " << ntBcOut << std::endl;

  for(Int_t it=0; it<ntBcOut; ++it){
    auto track = DCAna->GetTrackBcOut(it);
    auto chisqr = track->GetChiSquare();
    hddaq::cout << "       " << it << "-th track, chi2 = "
                << chisqr << std::endl;
    // auto nh = track->GetNHit();
    // for(Int_t ih=0; ih<nh; ++ih){
    //   auto hit = track->GetHit(ih);
    //   Int_t layerId = hit->GetLayer();
    //   Double_t wire = hit->GetWire();
    //   Double_t res = hit->GetResidual();
    //   hddaq::cout << "       layer = " << layerId << ", wire = "
    //             << wire << ", res = " << res << std::endl;
    // }
    gEvDisp.DrawBcOutLocalTrack(track);
  }
  if(ntBcOut==0) {
    hddaq::cout << "[Warning] BcOutTrack is empty!" << std::endl;
    return true;
  }

  //________________________________________________________
  //___ BFTCluster
  hodoAna->DecodeBFTHits(rawData);
  hodoAna->TimeCutBFT(MinTimeBFT, MaxTimeBFT);
  std::vector<Double_t> BftXCont;
  for(const auto& cl: hodoAna->GetClustersBFT()){
    BftXCont.push_back(cl->MeanPosition());
  }
  if(BftXCont.empty()){
    hddaq::cout << "[Warning] BftXCont is empty!" << std::endl;
    return true;
  }

  std::vector<ThreeVector> KmPCont, KmXCont;

  //________________________________________________________
  //___ K18Tracking
  DCAna->TrackSearchK18D2U(BftXCont);
  Int_t ntK18 = DCAna->GetNTracksK18D2U();
  for(Int_t i=0; i<ntK18; ++i){
    auto track = DCAna->GetK18TrackD2U(i);
    if(!track) continue;
    Double_t x = track->Xtgt(), y = track->Ytgt();
    Double_t u = track->Utgt(), v = track->Vtgt();
    Double_t p = track->P3rd();
    Double_t pt = p/TMath::Sqrt(1.+u*u+v*v);
    ThreeVector Pos(x, y, 0.);
    ThreeVector Mom(pt*u, pt*v, pt);
    KmPCont.push_back(Mom);
    KmXCont.push_back(Pos);
  }
  if(ntK18 == 0){
    hddaq::cout << "[Warning] Km is empty!" << std::endl;
    // gEvDisp.GetCommand();
    return true;
  }

  //________________________________________________________
  //___ SdcInDCHit
  DCAna->DecodeSdcInHits(rawData);
  Double_t multi_SdcIn = 0.;
  for(Int_t layer=1; layer<=NumOfLayersSdcIn; ++layer){
    const auto& cont = DCAna->GetSdcInHC(layer);
    Int_t n = cont.size();
    multi_SdcIn += n;
    if(n > MaxMultiHitSdcIn) continue;
    for(const auto& hit: cont){
      Int_t wire = hit->GetWire();
      Int_t mhit = hit->GetTdcSize();
      Bool_t goodFlag = false;
      for(Int_t j=0; j<mhit; j++){
        if(hit->IsWithinRange(j)){
          goodFlag = true;
          break;
        }
      }
      if(goodFlag) gEvDisp.DrawHitWire(layer, wire);
      else gEvDisp.DrawHitWire(layer, wire, false, false);
    }
  }
  multi_SdcIn /= (Double_t)NumOfLayersSdcIn;
  if(multi_SdcIn > MaxMultiHitSdcIn){
    hddaq::cout << "[Warning] SdcInHits exceed MaxMultiHit "
                << multi_SdcIn << "/" << MaxMultiHitSdcIn << std::endl;
    // gEvDisp.GetCommand();
    return true;
  }

  //________________________________________________________
  //___ SdcInTracking
  DCAna->TrackSearchSdcIn();
  Int_t ntSdcIn = DCAna->GetNtracksSdcIn();
  hddaq::cout << "[Info] ntSdcIn = " << ntSdcIn << std::endl;

  for(Int_t it=0; it<ntSdcIn; ++it){
    auto track = DCAna->GetTrackSdcIn(it);
    auto chisqr = track->GetChiSquare();
    hddaq::cout << "       " << it << "-th track, chi2 = "
                << chisqr << std::endl;
    // auto nh = track->GetNHit();
    // for(Int_t ih=0; ih<nh; ++ih){
    //   auto hit = track->GetHit(ih);
    //   Int_t layerId = hit->GetLayer();
    //   Double_t wire = hit->GetWire();
    //   Double_t res = hit->GetResidual();
    //   hddaq::cout << "       layer = " << layerId << ", wire = "
    //             << wire << ", res = " << res << std::endl;
    // }
    gEvDisp.DrawSdcInLocalTrack(track);
  }
  if(ntSdcIn != 1){
    hddaq::cout << "[Warning] SdcInTrack is empty!" << std::endl;
    return true;
  }

  //________________________________________________________
  //___ SdcOutDCHit
  DCAna->DecodeSdcOutHits(rawData);
  DCAna->TotCutSDC3(MinTotSDC3);
  DCAna->TotCutSDC4(MinTotSDC4);
  Double_t multi_SdcOut = 0.;
  for(Int_t layer=1; layer<=NumOfLayersSdcOut; ++layer){
    const auto& cont = DCAna->GetSdcOutHC(layer);
    Int_t n = cont.size();
    multi_SdcOut += n;
    if(n > MaxMultiHitSdcOut) continue;
    for(const auto& hit: cont){
      Int_t  wire = hit->GetWire();
      Int_t  mhit = hit->GetDriftTimeSize();
      Bool_t is_good = false;
      for(Int_t j=0; j<mhit && !is_good; ++j){
        is_good = hit->IsWithinRange(j);
      }
      if(is_good) gEvDisp.DrawHitWire(layer+30, wire);
      else gEvDisp.DrawHitWire(layer+30, wire, false, false);
    }
  }
  multi_SdcOut /= (Double_t)NumOfLayersSdcOut;
  if(multi_SdcOut > MaxMultiHitSdcOut){
    hddaq::cout << "[Warning] SdcOutHits exceed MaxMultiHit "
                << multi_SdcOut << "/" << MaxMultiHitSdcOut << std::endl;
    // gEvDisp.GetCommand();
    return true;
  }

  //________________________________________________________
  //___ SdcOutTracking
  DCAna->TrackSearchSdcOut(TOFCont);
  Int_t ntSdcOut = DCAna->GetNtracksSdcOut();
  hddaq::cout << "[Info] ntSdcOut = " << ntSdcOut << std::endl;
  for(Int_t it=0; it<ntSdcOut; ++it){
    auto track = DCAna->GetTrackSdcOut(it);
    auto chisqr = track->GetChiSquare();
    hddaq::cout << "       " << it << "-th track, chi2 = "
                << chisqr << std::endl;
    // auto nh = track->GetNHit();
    // for(Int_t ih=0; ih<nh; ++ih){
    //   auto hit = track->GetHit(ih);
    //   Int_t layerId = hit->GetLayer();
    //   Double_t wire = hit->GetWire();
    //   Double_t res = hit->GetResidual();
    //   hddaq::cout << "       layer = " << layerId << ", wire = "
    //             << wire << ", res = " << res << std::endl;
    // }
    gEvDisp.DrawSdcOutLocalTrack(track);
  }
  if(ntSdcOut != 1){
    hddaq::cout << "[Warning] SdcInTrack is empty!" << std::endl;
    return true;
  }
  //________________________________________________________
  //___ HTOFRawHit
  std::vector<Int_t> HTOFSegCont;
  Int_t nhHtof = 0;
  for(const auto& hit: rawData->GetHTOFRawHC()){
    if(!hit) continue;
    Int_t seg = hit->SegmentId();
    Bool_t is_hit_u = false;
    for(const auto& tdc: hit->GetArrayTdc1()){
      if(MinTdcHTOF < tdc && tdc < MaxTdcHTOF) is_hit_u = true;
    }
    Bool_t is_hit_d = false;
    for(const auto& tdc: hit->GetArrayTdc2()){
      if(MinTdcHTOF < tdc && tdc < MaxTdcHTOF) is_hit_d = true;
    }
    if(is_hit_u && is_hit_d){
      ++nhHtof;
      Int_t binid = 0;
      if(seg == 0) binid=1;
      else if(seg == 1 || seg == 2) binid = 2;
      else if(seg == 3 || seg == 4) binid = 3;
      else binid = seg-1;
      gEvDisp.FillHTOF(binid);
      HTOFSegCont.push_back(seg);
    }
  }
  if(nhHtof==0) {
    hddaq::cout << "[Warning] HtofCont is empty!" << std::endl;
    return true;
  }
  //________________________________________________________
  //___ TPCRawHit
  rawData->DecodeTPCHits();
  Int_t npadTpc = 0;
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    const auto hc = rawData->GetTPCRawHC(layer);
    const auto nhit = hc.size();
    npadTpc += nhit;
    for (const auto& rhit : hc) {
      Int_t layer = rhit->LayerId();
      Int_t row = rhit->RowId();
      auto mean = rhit->Mean(0, NumOfTimeBucket);
      auto max_adc = rhit->MaxAdc(0, NumOfTimeBucket);
      // auto rms = rhit->RMS(0, NumOfTimeBucket);
      auto loc_max = rhit->LocMax(0, NumOfTimeBucket);
      if(loc_max < 25 || 155 <loc_max)
        continue;
      gEvDisp.FillTPCADC(layer, row, max_adc - mean);
      gEvDisp.FillTPCTDC(layer, row, loc_max);
      // TVector3 pos = tpc::getPosition(layer, row);
      // pos.SetY((loc_max - 76.75)*80.0*0.05);
      // gEvDisp.SetTPCMarker(pos);
    }
  }
  gEvDisp.Update();

  // if(nhHtof >= 4){
  //   gEvDisp.GetCommand();
  // }
  std::vector<ThreeVector> KpPCont, KpXCont;
  std::vector<Double_t> M2Cont;
  std::vector<Double_t> Chi2KuramaCont;
  //________________________________________________________
  //___ KuramaTracking
  static const auto StofOffset = gUser.GetParameter("StofOffset");
  // DCAna->SetMaxV0Diff(10.);
  std::cout<<"kurama check1"<<std::endl;
  DCAna->TrackSearchKurama();
  std::cout<<"kurama check2"<<std::endl;
  Bool_t through_target = false;

  Int_t ntKurama = DCAna->GetNTracksKurama();
  hddaq::cout << "[Info] ntKurama = " << ntKurama << std::endl;
  for(Int_t it=0; it<ntKurama; ++it){
    auto track = DCAna->GetKuramaTrack(it);
    // track->Print();
    auto chisqr = track->GetChiSquare();
    hddaq::cout << "       " << it << "-th track, chi2 = "
                << chisqr << std::endl;
    const auto& postgt = track->PrimaryPosition();
    const auto& momtgt = track->PrimaryMomentum();
    Double_t path = track->PathLengthToTOF();
    Double_t p = momtgt.Mag();
    gEvDisp.FillMomentum(p);
    if(chisqr > 20.) continue;
    if(TMath::Abs(postgt.x()) < 30.
       && TMath::Abs(postgt.y()) < 20.){
      through_target = true;
    }
    // MassSquare
    Double_t tofseg = track->TofSeg();
    for(const auto& hit: TOFCont){
      Double_t seg = hit->SegmentId()+1;
      if(tofseg != seg) continue;
      Double_t stof = hit->CMeanTime()-ctime0+StofOffset;
      if(stof <= 0) continue;
      Double_t m2 = Kinematics::MassSquare(p, path, stof);
      gEvDisp.FillMassSquare(m2);
      KpPCont.push_back(momtgt);
      KpXCont.push_back(postgt);
      M2Cont.push_back(m2);
      Chi2KuramaCont.push_back(chisqr);
    }
  }
  if(KpPCont.size() == 0){
    hddaq::cout << "[Warning] Kp is empty!" << std::endl;
    // gEvDisp.GetCommand();
    return true;
  }
  if(through_target) gEvDisp.DrawTarget();

#if(DRAWTEXT_DETECTOR)
  TString buf;
  std::stringstream ss; ss << trigger_flag;
  buf = ss.str();
  buf.ReplaceAll("0", ".").ReplaceAll("1", "!");
  gEvDisp.DrawText(0.040, 0.960, Form("TrigFlag   %s", buf.Data()));
  buf = "BH1Seg  ";
  for(const auto& seg: BH1SegCont){
    buf += Form(" %d", seg);
  }
  gEvDisp.DrawText(0.040, 0.920, buf);
  buf = "BH2Seg  ";
  for(const auto& seg: BH2SegCont){
    buf += Form(" %d", seg);
  }
  gEvDisp.DrawText(0.040, 0.880, buf);
  buf = "HTOFSeg  ";
  for(const auto& seg: HTOFSegCont){
    buf += Form(" %d", seg);
  }
  gEvDisp.DrawText(0.040, 0.840, buf);
  buf = "TOFSeg  ";
  for(const auto& seg: TOFSegCont){
    buf += Form(" %d", seg);
  }
  gEvDisp.DrawText(0.040, 0.800, buf);
  buf = "WCSeg  ";
  for(const auto& seg: WCSegCont){
    buf += Form(" %d", seg);
  }
  gEvDisp.DrawText(0.040, 0.760, buf);
  buf = "BcOut"; gEvDisp.DrawText(0.040, 0.280, buf);
  buf= "#chi^{2} = ";
  for(Int_t i=0; i<ntBcOut; ++i){
    buf += Form(" %.3f", DCAna->GetTrackBcOut(i)->GetChiSquare());
  }
  gEvDisp.DrawText(0.130, 0.280, buf);
  buf = "SdcIn"; gEvDisp.DrawText(0.040, 0.240, buf);
  buf = "#chi^{2} = ";
  for(Int_t i=0; i<ntSdcIn; ++i){
    buf += Form(" %.3f", DCAna->GetTrackSdcIn(i)->GetChiSquare());
  }
  gEvDisp.DrawText(0.130, 0.240, buf);
  buf = "SdcOut"; gEvDisp.DrawText(0.040, 0.20, buf);
  buf = "#chi^{2} = ";
  for(Int_t i=0; i<ntSdcOut; ++i){
    buf += Form(" %.3f", DCAna->GetTrackSdcOut(i)->GetChiSquare());
  }
  gEvDisp.DrawText(0.130, 0.200, buf);
  buf = "Kurama"; gEvDisp.DrawText(0.040, 0.16, buf);
  buf = "#chi^{2} = ";
  for(Int_t i=0; i<ntKurama; ++i){
    buf += Form(" %.3f", DCAna->GetKuramaTrack(i)->GetChiSquare());
  }
  gEvDisp.DrawText(0.130, 0.160, buf);
  gEvDisp.DrawText(0.680, 0.960, "BTOF");
  gEvDisp.DrawText(0.860, 0.960, Form("%.3f", btof));
#endif
  //________________________________________________________
  //___ Reaction
  Bool_t is_good = false;
  if(KmPCont.size()==1 && KpPCont.size()==1){
    ThreeVector pkp = KpPCont[0];
    ThreeVector pkm = KmPCont[0];
    ThreeVector xkp = KpXCont[0];
    ThreeVector xkm = KmXCont[0];
    Double_t m2 = M2Cont[0];
    Double_t mass = TMath::QuietNaN();
    if(TMath::Abs(m2) < 0.15 && pkp.Mag() < 1.5) mass = PionMass;
    if(m2 > 0.15 && m2 < 0.35 && pkp.Mag() < 1.4) mass = KaonMass;
    if(m2 > 0.55) mass = ProtonMass;
    ThreeVector vertex = Kinematics::VertexPoint(xkm, xkp, pkm, pkp);
    Double_t closedist = Kinematics::CloseDist(xkm, xkp, pkm, pkp);
    LorentzVector LvKm(KmPCont[0], TMath::Sqrt(mass*mass+pkm.Mag2()));
    LorentzVector LvKp(KpPCont[0], TMath::Sqrt(mass*mass+pkp.Mag2()));
    LorentzVector LvP(0., 0., 0., ProtonMass);
    LorentzVector LvRp = LvKm+LvP-LvKp;
    ThreeVector MissMom = LvRp.Vect();
    Double_t MissMass = LvRp.Mag();

    hddaq::cout << "[Info] Vertex = " << vertex << std::endl;
    hddaq::cout << "[Info] MissingMomentum = " << MissMom << std::endl;

    gEvDisp.DrawVertex(vertex);
    gEvDisp.DrawMissingMomentum(MissMom, vertex);
    if(true
       && TMath::Abs(vertex.x()-8.9) < 25.
       && TMath::Abs(vertex.y()) < 20.
       && TMath::Abs(vertex.z()+70) < 100
       && closedist < 20.
       // && through_target
    ){
#if(DRAWTEXT_TRACKING)
      gEvDisp.DrawText(0.680, 0.920, "pK18");
      gEvDisp.DrawText(0.860, 0.920, Form("%.3f", pkm.Mag()));
      gEvDisp.DrawText(0.680, 0.880, "pKurama");
      gEvDisp.DrawText(0.860, 0.880, Form("%.3f", pkp.Mag()));
      gEvDisp.DrawText(0.680, 0.840, "MassSquared");
      gEvDisp.DrawText(0.860, 0.840, Form("%.3f", m2));
      gEvDisp.DrawText(0.660, 0.280, "CloseDist");
      gEvDisp.DrawText(0.770, 0.280, Form("%.2f", closedist));
      gEvDisp.DrawText(0.660, 0.240, "Vertex");
      gEvDisp.DrawText(0.770, 0.240, Form("(%.2f, %.2f, %.2f)",
                                          vertex.X(), vertex.Y(), vertex.Z()));
      gEvDisp.DrawText(0.660, 0.200, "MissMom");
      gEvDisp.DrawText(0.770, 0.200, Form("(%.3f, %.3f, %.3f)",
                                          MissMom.X(), MissMom.Y(), MissMom.Z()));
      gEvDisp.DrawText(0.660, 0.160, "MissMass");
      gEvDisp.DrawText(0.770, 0.160, Form("%.4f", MissMass));
#endif
      if(true
         && mass == KaonMass
         && pkp.z() > 0
      ){
        // gEvDisp.GetCommand();
        is_good = true;
      }
    }
  }

  gEvDisp.Update();
  // gEvDisp.GetCommand();
  hddaq::cout << "[Info] IsGood = " << is_good << std::endl;

  if(is_good){
#if SAVEPDF
    gEvDisp.Print(gUnpacker.get_run_number(),
                  gUnpacker.get_event_number());
#else
    gSystem->Sleep(10000);
#endif
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
UserEventDisplay::ProcessingEnd()
{
  // gEvDisp.GetCommand();
  gEvDisp.EndOfEvent();
  // if(utility::UserStop()) gEvDisp.Run();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserEventDisplay;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  gUnpacker.disable_istream_bookmark();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
#if SAVEPDF
  gEvDisp.SetSaveMode();
#endif
  return
    (InitializeParameter<DCGeomMan>("DCGEO")
     && InitializeParameter<DCDriftParamMan>("DCDRFT")
     && InitializeParameter<DCTdcCalibMan>("DCTDC")
     && InitializeParameter<HodoParamMan>("HDPRM")
     && InitializeParameter<HodoPHCMan>("HDPHC")
     && InitializeParameter<FieldMan>("FLDMAP", "HSFLDMAP")
     && InitializeParameter<K18TransMatrix>("K18TM")
     && InitializeParameter<BH2Filter>("BH2FLT")
     && InitializeParameter<UserParamMan>("USER")
     && InitializeParameter<EventDisplay>());
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

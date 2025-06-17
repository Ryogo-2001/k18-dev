// -*- C++ -*-

#include "VEvent.hh"

#include <iostream>

#include <TMath.h>

#include <Unpacker.hh>
#include <UnpackerManager.hh>

#include "ConfMan.hh"
#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "TPCHit.hh"
#include "TPCPadHelper.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "TPCRawHit.hh"
#include "UserParamMan.hh"

//#define GateCalib 1
#define GateCalib 0
#define GainCalib 0
#define Srdata 0
//#define GainCalib 0

namespace
{
using namespace root;
const auto qnan = TMath::QuietNaN();
using hddaq::unpacker::GUnpacker;
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser     = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
class UserTPCHit : public VEvent
{
private:
  RawData*      rawData;
  HodoAnalyzer* hodoAna;
  DCAnalyzer*   DCAna;

public:
  UserTPCHit();
  ~UserTPCHit();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserTPCHit::ClassName()
{
  static TString s_name("UserTPCHit");
  return s_name;
}

//_____________________________________________________________________________
UserTPCHit::UserTPCHit()
  : VEvent(),
    rawData(new RawData),
    hodoAna(new HodoAnalyzer),
    DCAna(new DCAnalyzer)
{
}

//_____________________________________________________________________________
UserTPCHit::~UserTPCHit()
{
  if(rawData) delete rawData;
  if(hodoAna) delete hodoAna;
  if(DCAna) delete DCAna;
}

//_____________________________________________________________________________
struct Event
{
  Int_t                 runnum;
  Int_t                 evnum;
  std::vector<Int_t>    trigpat;
  std::vector<Int_t>    trigflag;
  Int_t                 npadTpc;   // number of pads
  Int_t                 nhTpc;     // number of hits

  Int_t browTpc;
  Int_t blayerTpc;
  Double_t brmsTpc;//Baseline RMS;
  // vector (size=nhTpc)
  std::vector<Int_t>    layerTpc;  // layer id
  std::vector<Int_t>    rowTpc;    // row id
  std::vector<Int_t>    padTpc;    // pad id
  std::vector<Double_t> pedTpc;    // pedestal
  std::vector<Double_t> rmsTpc;    // rms
  std::vector<Double_t> rawrmsTpc; // rawrms
  std::vector<Double_t> deTpc;     // dE
  std::vector<Double_t> sigmaTpc;  // sigma
  std::vector<Double_t> tTpc;      // time
  std::vector<Double_t> chisqrTpc; // chi^2 of signal fitting
  std::vector<Double_t> cdeTpc;    // dE
  std::vector<Double_t> ctTpc;     // time
  std::vector<Double_t> dlTpc;     // time
  std::vector<Double_t> clkTpc;    // clock timing

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  runnum  = 0;
  evnum   = 0;
  npadTpc = 0;
  nhTpc   = 0;
  browTpc   = -1;
  blayerTpc = -1;
  brmsTpc   = 0;
  trigpat.clear();
  trigflag.clear();
  layerTpc.clear();
  rowTpc.clear();
  padTpc.clear();
  pedTpc.clear();
  rmsTpc.clear();
  rawrmsTpc.clear();
  deTpc.clear();
  sigmaTpc.clear();
  tTpc.clear();
  chisqrTpc.clear();
  cdeTpc.clear();
  ctTpc.clear();
  dlTpc.clear();
  clkTpc.clear();
}

//_____________________________________________________________________________
namespace root
{
Event event;
TH1   *h[MaxHist];
TTree *tree;
enum eDetHid {
  PadHid    = 100000,
};
}

//_____________________________________________________________________________
Bool_t
UserTPCHit::ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserTPCHit::ProcessingNormal()
{
  static const Int_t MaxMultiHitTPC = gUser.GetParameter("MaxMultiHitTPC");
  static const Int_t NumOfTimeBucket = gUser.GetParameter("NumOfTimeBucket");
  const Int_t run_number   = gUnpacker.get_root()->get_run_number();
  const Int_t event_number = gUnpacker.get_event_number();

  event.runnum = run_number;
  event.evnum  = event_number;

  rawData->DecodeHits();
  event.trigpat.resize(NumOfSegTrig);
  event.trigflag.resize(NumOfSegTrig);
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

  if(trigger_flag[trigger::kSpillEnd]) return true;

  rawData->DecodeTPCHits();

  //___ TPC Clock
  for(const auto& hit: rawData->GetTPCClockRawHC()){
    for(const auto& tdc: hit->GetArrayTdc1()){
      HF1(501, tdc);
    }
  }

  Double_t clock_timing = 0.;
  if(hodoAna->DecodeTPCClock(rawData)){
    clock_timing = hodoAna->GetHitTPCClock()->CTime();
  }

  event.clkTpc.push_back(clock_timing);
  HF1(502, clock_timing);

  // {
  //   bool maxadccut = true; // max adccut
  //   //bool maxadctbcut = true; // max adccut time bucket cut
  //   // bool maxadccut = false; // max adccut
  //   bool maxadctbcut = false; // max adccut time bucket cut
  //   rawData->SelectTPCHits(maxadccut, maxadctbcut);
  // }

  HF1(1, 0);

  //________________________________________________________
  //___ TPCRawHit
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    auto hc = rawData->GetTPCRawHC(layer);
    for(const auto& rhit : hc){
      auto mean    = rhit->Mean(0, NumOfTimeBucket);
      auto max_adc = rhit->MaxAdc(0, NumOfTimeBucket);
      auto min_adc = rhit->MinAdc(0, NumOfTimeBucket);
      auto rms     = rhit->RMS(0, NumOfTimeBucket);
      auto loc_max = rhit->LocMax(0, NumOfTimeBucket);
      HF1(11, mean);
      HF1(12, max_adc);
      HF1(13, rms);
      HF1(14, loc_max);
      HF1(15, min_adc);
      auto fadc = rhit->Fadc();
      for(Int_t tb=0, ntb=fadc.size(); tb<ntb; ++tb){
        HF2(121, tb, fadc.at(tb));
      }
    }
  }

  //________________________________________________________
  //___ TPCRawHit after baseline correction
  auto baseline = rawData->GetBaselineTPC();
  if(baseline){
    auto fadc = baseline->Fadc();
    for(Int_t tb=0, ntb=fadc.size(); tb<ntb; ++tb){
      HF2(39, tb, fadc.at(tb));
    }
    auto browTpc = baseline->RowId();
    auto blayerTpc = baseline->LayerId();
    auto brmsTpc = baseline->RMS(0, NumOfTimeBucket);
    event.browTpc = browTpc;
    event.blayerTpc = blayerTpc;
    event.brmsTpc = brmsTpc;
  }
  Int_t npadTpc = 0;
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    auto hc = rawData->GetTPCCorHC(layer);
    const auto nhit = hc.size();
    npadTpc += nhit;
    for(const auto& rhit : hc){
      auto mean    = rhit->Mean(0, NumOfTimeBucket);
      auto max_adc = rhit->MaxAdc(0, NumOfTimeBucket);
      auto min_adc = rhit->MinAdc(0, NumOfTimeBucket);
      auto rms     = rhit->RMS(0, NumOfTimeBucket);
      auto loc_max = rhit->LocMax(0, NumOfTimeBucket);
      auto pars    = rhit->GetParameters();
      HF1(31, mean);
      HF1(32, max_adc);
      HF1(33, rms);
      HF1(34, loc_max);
      HF1(35, min_adc);
      HF1(36, pars.at(0));
      HF1(37, pars.at(1));
      HF1(38, pars.at(2));
      auto fadc = rhit->Fadc();
      for(Int_t tb=0, ntb=fadc.size(); tb<ntb; ++tb){
        HF2(122, tb, fadc.at(tb));
      }
    }
  }
  event.npadTpc = npadTpc;
  HF1(10, npadTpc);

  HF1(1, 19);

  //________________________________________________________
  //___ TPCHit
  if(MaxMultiHitTPC>0 && npadTpc>MaxMultiHitTPC){
    std::cout << "#W Too many hits found, npadTpc = " << npadTpc << std::endl;
    return true;
  }
  DCAna->DecodeTPCHits(rawData, clock_timing);

  Int_t nhTpc = 0;
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    auto hc = DCAna->GetTPCHC(layer);
    for(const auto& hit : hc){
      if(!hit || !hit->IsGood())
        continue;
      // hit->Print();
      //Int_t layer = hit->GetLayer();
      Int_t row = hit->GetWire();
      Int_t pad = tpc::GetPadId(layer, row);
      Double_t ped = hit->GetPedestal();
      Double_t rms = hit->GetRMS();
      Double_t rawrms = hit->GetRawRMS();
      HF1(101, ped);
      HF1(103, rms);
      const auto& vec = tpc::getPosition(pad);
      HF2Poly(1001, vec.Z(), vec.X());
      Int_t nhit = hit->GetNHits();
      Bool_t good_for_analysis = false;
      for(Int_t i=0; i<nhit; ++i){
	Double_t cde = hit->GetCDe(i);
	Double_t de = hit->GetDe(i);
        Double_t time = hit->GetTime(i);
        Double_t chisqr = hit->GetChisqr(i);
	Double_t ctime = hit->GetCTime(i);
        Double_t dl = hit->GetDriftLength(i);
        Double_t sigma = hit->GetSigma(i);
        event.layerTpc.push_back(layer);
        event.rowTpc.push_back(row);
        event.padTpc.push_back(pad);
        event.pedTpc.push_back(ped);
        event.rmsTpc.push_back(rms);
	event.rawrmsTpc.push_back(rawrms);
        event.deTpc.push_back(de);
        event.tTpc.push_back(time);
        event.chisqrTpc.push_back(chisqr);
        event.cdeTpc.push_back(cde);
        event.ctTpc.push_back(ctime);
        event.dlTpc.push_back(dl);
	event.sigmaTpc.push_back(sigma);
	// if(69.<time&&time<85.&&nhit==1)
	HF1(102, de);
        HF1(104, time);
        HF1(105, chisqr);
        HF1(106, cde);
        HF1(107, ctime);
        HF1(108, dl);
        HF1(109, sigma);
        HF2(110, de, sigma);
        HF2(111, de, time);

#if GateCalib
	HF1(PadHid + layer*1000 + row, time);
#endif

#if GainCalib
	//	if(69.<time&&time<85.&&nhit==1)
	HF1(2*PadHid + layer*1000 + row, de);
#endif

#if Srdata
	if(69.<time&&time<85.&&nhit==1&&layer==31)
	  HF1(3*PadHid + layer*1000 + row, de);
#endif

        good_for_analysis = true;
        ++nhTpc;
      }
      if(good_for_analysis){
        auto fadc = hit->GetRawHit()->Fadc();
        for(Int_t tb=0, ntb=fadc.size(); tb<ntb; ++tb){
          HF2(123, tb, fadc.at(tb));
        }
      }
    }
  }

  event.nhTpc = nhTpc;
  HF1(100, nhTpc);

  return true;
}

//_____________________________________________________________________________
Bool_t
UserTPCHit::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserTPCHit;
}

//_____________________________________________________________________________
Bool_t
ConfMan:: InitializeHistograms()
{
  const Int_t    NbinAdc = 4096;
  const Double_t MinAdc  =    0.;
  const Double_t MaxAdc  = 4096.;
  const Int_t    NbinRms = 1000;
  const Double_t MinRms  =    0.;
  const Double_t MaxRms  = 1000.;
  const Int_t    NbinDe = 1000;
  const Double_t MinDe  =    0.;
  const Double_t MaxDe  = 1000.;
  const Int_t    NbinChisqr = 1000;
  const Double_t MinChisqr  =    0.;
  const Double_t MaxChisqr  = 1000.;
  const Int_t    NbinTime = 1000;
  const Double_t MinTime  = -8000.;
  const Double_t MaxTime  =  8000.;
  const Int_t    NbinDL = 800;
  const Double_t MinDL  = -400.;
  const Double_t MaxDL  =  400.;
  const Int_t    NbinSigma = 500;
  const Double_t MinSigma  = 0.;
  const Double_t MaxSigma  =  50.;
  const Int_t    NTimeBucket = 170;
  HB1(1, "Status", 20, 0., 20.);
  HB1(10, "TPC Multiplicity (Raw)", NumOfPadTPC+1, 0, NumOfPadTPC+1);
  HB1(11, "TPC FADC Mean", NbinAdc, MinAdc, MaxAdc);
  HB1(12, "TPC FADC Max", NbinAdc, MinAdc, MaxAdc);
  HB1(13, "TPC FADC RMS", NbinRms, MinRms, MaxRms);
  HB1(14, "TPC FADC LocMax", NTimeBucket+1, 0, NTimeBucket+1);
  HB1(15, "TPC FADC Min", NbinAdc, MinAdc, MaxAdc);
  HB1(31, "TPC FADC Mean Cor", NbinAdc, MinAdc, MaxAdc);
  HB1(32, "TPC FADC Max Cor", NbinAdc, MinAdc, MaxAdc);
  HB1(33, "TPC FADC RMS Cor", NbinRms, MinRms, MaxRms);
  HB1(34, "TPC FADC LocMax Cor", NTimeBucket+1, 0, NTimeBucket+1);
  HB1(35, "TPC FADC Min Cor", NbinAdc, MinAdc, MaxAdc);
  HB1(36, "TPC FADC Baseline p0", NbinAdc, MinAdc, MaxAdc);
  HB1(37, "TPC FADC Baseline p1", 120, -6, 6);
  HB1(38, "TPC FADC Baseline p2", 120, -12, 12);
  HB2(39, "TPC FADC Baseline",
      NTimeBucket+1, 0, NTimeBucket+1, NbinAdc, MinAdc, MaxAdc);

  HB1(100, "TPC Multiplicity (TPCHit)", NumOfPadTPC+1, 0, NumOfPadTPC+1);
  HB1(101, "TPC Pedestal", NbinAdc, MinAdc, MaxAdc);
  HB1(102, "TPC DeltaE", NbinDe, MinDe, MaxDe);
  HB1(103, "TPC RMS", NbinRms, MinRms, MaxRms);
  HB1(104, "TPC Time", (NTimeBucket+1)*30, 0, NTimeBucket+1);
  HB1(105, "TPC Chisqr", NbinChisqr, MinChisqr, MaxChisqr);
  HB1(106, "TPC CDeltaE", NbinDe, MinDe, MaxDe);
  HB1(107, "TPC CTime", NbinTime, MinTime, MaxTime);
  HB1(108, "TPC DriftLength", NbinDL, MinDL, MaxDL);
  HB1(109, "TPC sigma", NbinSigma, MinSigma, MaxSigma);
  HB2(110, "TPC sigma%de", NbinDe, MinDe, MaxDe, NbinSigma, MinSigma, MaxSigma);
  HB2(111, "TPC time%de", NbinDe, MinDe, MaxDe, NbinTime, MinTime, MaxTime);
  HB2(121, "TPC FADC (Before)",
      NTimeBucket+1, 0, NTimeBucket+1, NbinAdc, MinAdc, MaxAdc);
  HB2(122, "TPC FADC (After)",
      NTimeBucket+1, 0, NTimeBucket+1, NbinAdc, MinAdc-500, MaxAdc-500);
  HB2(123, "TPC FADC (Good)",
      NTimeBucket+1, 0, NTimeBucket+1, NbinAdc, MinAdc, MaxAdc);

  HB1(501, "TPC Clock TDC", 100000, 0., 1000000.);
  HB1(502, "TPC Clock Time", 20000, -100., 100.);

  HB2Poly(1001, "TPC HitPat");

  tpc::InitializeHistograms();

#if GateCalib
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
    for(Int_t r=0; r<NumOfRow; ++r){
      HB1(PadHid + layer*1000 + r , "TPC Time", (NTimeBucket+1)*30, 0, NTimeBucket+1);
    }
  }
#endif

#if GainCalib
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
    for(Int_t r=0; r<NumOfRow; ++r){
      HB1(2*PadHid + layer*1000 + r , "TPC DeltaE", NbinDe, MinDe, MaxDe);
    }
  }
#endif

#if Srdata
  //  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
  Int_t layer= 31;
  const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
  for(Int_t r=0; r<NumOfRow; ++r){
    HB1(3*PadHid + layer*1000 + r , "TPC DeltaE", NbinDe, MinDe, MaxDe);
  }
  //}
#endif

  // Tree
  HBTree("tpc", "tree of TPCHit");
  tree->Branch("runnum", &event.runnum);
  tree->Branch("evnum", &event.evnum);
  tree->Branch("trigpat", &event.trigpat);
  tree->Branch("trigflag", &event.trigflag);
  tree->Branch("npadTpc", &event.npadTpc);
  tree->Branch("nhTpc", &event.nhTpc);
  tree->Branch("browTpc", &event.browTpc);
  tree->Branch("blayerTpc", &event.blayerTpc);
  tree->Branch("brmsTpc", &event.brmsTpc);
  tree->Branch("layerTpc", &event.layerTpc);
  tree->Branch("rowTpc", &event.rowTpc);
  tree->Branch("padTpc", &event.padTpc);
  tree->Branch("pedTpc", &event.pedTpc);
  tree->Branch("rmsTpc", &event.rmsTpc);
  tree->Branch("rawrmsTpc", &event.rawrmsTpc);
  tree->Branch("deTpc", &event.deTpc);
  tree->Branch("tTpc", &event.tTpc);
  tree->Branch("chisqrTpc", &event.chisqrTpc);
  tree->Branch("cdeTpc", &event.cdeTpc);
  tree->Branch("ctTpc", &event.ctTpc);
  tree->Branch("dlTpc", &event.dlTpc);
  tree->Branch("sigmaTpc", &event.sigmaTpc);
  tree->Branch("clkTpc", &event.clkTpc);

  // HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO") &&
     InitializeParameter<HodoParamMan>("HDPRM") &&
     InitializeParameter<HodoPHCMan>("HDPHC") &&
     InitializeParameter<TPCParamMan>("TPCPRM") &&
     InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

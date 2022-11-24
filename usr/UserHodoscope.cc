// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iostream>
#include <sstream>

#include <UnpackerManager.hh>

#include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "KuramaLib.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "DCGeomMan.hh"

namespace
{
using namespace root;
const auto qnan = TMath::QuietNaN();
auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
auto& gRM = RMAnalyzer::GetInstance();
auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
class UserHodoscope : public VEvent
{
private:
  RawData*      rawData;
  HodoAnalyzer* hodoAna;

public:
  UserHodoscope();
  ~UserHodoscope();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserHodoscope::ClassName()
{
  static TString s_name("UserHodoscope");
  return s_name;
}

//_____________________________________________________________________________
UserHodoscope::UserHodoscope()
  : VEvent(),
    rawData(new RawData),
    hodoAna(new HodoAnalyzer)
{
}

//_____________________________________________________________________________
UserHodoscope::~UserHodoscope()
{
  if(hodoAna) delete hodoAna;
  if(rawData) delete rawData;
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;
  Int_t spill;

  Double_t adc[NumOfSegCaenV792];

  void clear()
  {
    evnum      = 0;
    spill      = 0;
    for(Int_t i=0; i<NumOfSegCaenV792; ++i){
      adc[i] = qnan;
    }
  }
};

//_____________________________________________________________________________
namespace root
{
Event  event;
TH1*   h[MaxHist];
TTree* tree;
enum EHid {
  kCaenV792Hid = 10000
};
}

//_____________________________________________________________________________
Bool_t
UserHodoscope::ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserHodoscope::ProcessingNormal()
{
  rawData->DecodeHits();

  event.evnum = gUnpacker.get_event_number();

  HF1(1, 0);

  HF1(1, 1);

  {
    const HodoRHitContainer& cont = rawData->GetCaenV792RawHC();
    for(Int_t i=0, n=cont.size(); i<n; ++i){
      const auto& hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Int_t Au = hit->GetAdcUp();
      HF1(kCaenV792Hid +100*seg +1, Au);
      event.adc[seg-1] = Au;
      // Bool_t is_hit_u = false;
      // for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
      //   Int_t T = hit->GetTdcUp(m);
      //   HF1(BH1Hid +100*seg +3, T);
      //   event.bh1ut[seg-1][m] = T;
      //   if(MinTdcBH1 < T && T < MaxTdcBH1) is_hit_u = true;
      // }
      // if(is_hit_u) HF1(BH1Hid +100*seg +5, Au);
      // else         HF1(BH1Hid +100*seg +7, Au);
      // // Down
      // Int_t Ad = hit->GetAdcDown();
      // HF1(BH1Hid +100*seg +2, Ad);
      // event.bh1da[seg-1] = Ad;
      // Bool_t is_hit_d = false;
      // for(Int_t m=0, n_mhit=hit->GetSizeTdcDown(); m<n_mhit; ++m){
      //   Int_t T = hit->GetTdcDown(m);
      //   HF1(BH1Hid +100*seg +4, T);
      //   event.bh1dt[seg-1][m] = T;
      //   if(MinTdcBH1 < T && T < MaxTdcBH1) is_hit_d = true;
      // }
      // if(is_hit_d) HF1(BH1Hid +100*seg +6, Ad);
      // else         HF1(BH1Hid +100*seg +8, Ad);
      // // HitPat
      // if(is_hit_u || is_hit_d){
      //   ++nh1; HF1(BH1Hid +3, seg-0.5);
      // }
      // if(is_hit_u && is_hit_d){
      //   event.bh1hitpat[bh1_nhits++] = seg;
      //   ++nh2; HF1(BH1Hid +5, seg-0.5);
      // }
    }
    // HF1(BH1Hid +2, nh1); HF1(BH1Hid +4, nh2);
    // event.bh1nhits = bh1_nhits;
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
UserHodoscope::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserHodoscope;
}

//_____________________________________________________________________________
namespace
{
const Int_t    NbinAdc = 4096;
const Double_t MinAdc  =    0.;
const Double_t MaxAdc  = 4096.;

const Int_t    NbinTdc = 4096;
const Double_t MinTdc  =    0.;
const Double_t MaxTdc  = 4096.;

const Int_t    NbinTdcHr = 1e6/10;
const Double_t MinTdcHr  =  0.;
const Double_t MaxTdcHr  = 1e6;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  HB1(1, "Status", 20, 0., 20.);

  // Rawdata
  HB1(kCaenV792Hid +0, "#Hits CaenV792",
      NumOfSegCaenV792+1, 0., Double_t(NumOfSegCaenV792+1));
  HB1(kCaenV792Hid +1, "Hitpat CaenV792",
      NumOfSegCaenV792,   0., Double_t(NumOfSegCaenV792));
  HB1(kCaenV792Hid +2, "#Hits CaenV792(Tor)",
      NumOfSegCaenV792+1, 0., Double_t(NumOfSegCaenV792+1));
  HB1(kCaenV792Hid +3, "Hitpat CaenV792(Tor)",
      NumOfSegCaenV792,   0., Double_t(NumOfSegCaenV792));
  HB1(kCaenV792Hid +4, "#Hits CaenV792(Tand)",
      NumOfSegCaenV792+1, 0., Double_t(NumOfSegCaenV792+1));
  HB1(kCaenV792Hid +5, "Hitpat CaenV792(Tand)",
      NumOfSegCaenV792,   0., Double_t(NumOfSegCaenV792));

  for(Int_t i=1; i<=NumOfSegCaenV792; ++i){
    TString title1 = Form("CaenV792-%d UpAdc", i);
    TString title2 = Form("CaenV792-%d DownAdc", i);
    TString title3 = Form("CaenV792-%d UpTdc", i);
    TString title4 = Form("CaenV792-%d DownTdc", i);
    TString title5 = Form("CaenV792-%d UpAdc(w Tdc)", i);
    TString title6 = Form("CaenV792-%d DownAdc(w Tdc)", i);
    TString title7 = Form("CaenV792-%d UpAdc(w/o Tdc)", i);
    TString title8 = Form("CaenV792-%d DownAdc(w/o Tdc)", i);
    TString title9 = Form("CaenV792-%d UpTdc (Time0Seg==4)", i);
    TString title10= Form("CaenV792-%d DownTdc (Time0Seg==4)", i);
    HB1(kCaenV792Hid +100*i +1, title1, NbinAdc,   MinAdc,   MaxAdc);
    HB1(kCaenV792Hid +100*i +2, title2, NbinAdc,   MinAdc,   MaxAdc);
    HB1(kCaenV792Hid +100*i +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(kCaenV792Hid +100*i +4, title4, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(kCaenV792Hid +100*i +5, title5, NbinAdc,   MinAdc,   MaxAdc);
    HB1(kCaenV792Hid +100*i +6, title6, NbinAdc,   MinAdc,   MaxAdc);
    HB1(kCaenV792Hid +100*i +7, title7, NbinAdc,   MinAdc,   MaxAdc);
    HB1(kCaenV792Hid +100*i +8, title8, NbinAdc,   MinAdc,   MaxAdc);
    HB1(kCaenV792Hid +100*i +9, title9, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(kCaenV792Hid +100*i +10,title10,NbinTdcHr, MinTdcHr, MaxTdcHr);
  }

  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("spill",     &event.spill,     "spill/I");

  //BH1
  tree->Branch("adc", event.adc, Form("adc[%d]/D", NumOfSegCaenV792));

  // HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return true;
    // (InitializeParameter<DCGeomMan>("DCGEO") &&
    //  InitializeParameter<HodoParamMan>("HDPRM") &&
    //  InitializeParameter<HodoPHCMan>("HDPHC")   &&
    //  InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

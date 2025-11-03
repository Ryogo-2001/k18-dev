
// -*- C++ -*-
//#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>   // ★ 追加：max_element 用

#include "BH2Hit.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "RootHelper.hh"
#include "HodoHit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "S2sLib.hh"
#include "RawData.hh"
#include "UnpackerManager.hh"

#include "RayrawWaveformHit.hh"
#include "TemplateFitMan.hh"

namespace
{
  using namespace root;
  using hddaq::unpacker::GUnpacker;
  const auto qnan       = TMath::QuietNaN();
  const auto& gUnpacker = GUnpacker::get_instance();
  const auto& gUser     = UserParamMan::GetInstance();
  const auto& gHodo     = HodoParamMan::GetInstance();
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;

  // [seg][pulse]
  std::vector<std::vector<Double_t>> fit_time;
  std::vector<std::vector<Double_t>> fit_height;
  std::vector<std::vector<Double_t>> fit_dE;

  // [seg][sample]
  std::vector<std::vector<Double_t>> waveform_time;
  std::vector<std::vector<Double_t>> waveform_volt;

  void clear();
};

void Event::clear()
{
  evnum = 0;
  fit_time.clear();
  fit_height.clear();
  fit_dE.clear();
  waveform_time.clear();
  waveform_volt.clear();
}

//_____________________________________________________________________________
namespace root
{
  Event  event;
  TH1   *h[MaxHist];
  TTree *tree;
  enum eDetHid
    {
      RAYRAWHid  = 100000,
    };
}

//_____________________________________________________________________________
Bool_t
ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingNormal()
{
  RawData rawData;
  HodoAnalyzer hodoAna(rawData);

  
  event.evnum = gUnpacker.get_event_number();
  HF1(1, 0);

  
  rawData.DecodeHits("RAYRAW");
  hodoAna.DecodeHits<RayrawWaveformHit>("RAYRAW");

  const Int_t nh = hodoAna.GetNHits("RAYRAW");
  const Int_t nSeg = NumOfSegRayraw;

  
  event.fit_time.assign(nSeg, {});
  event.fit_height.assign(nSeg, {});
  event.fit_dE.assign(nSeg, {});
  event.waveform_time.assign(nSeg, {});
  event.waveform_volt.assign(nSeg, {});

  for (Int_t i = 0; i < nh; ++i) {
    const auto& hit = hodoAna.GetHit<RayrawWaveformHit>("RAYRAW", i);
    if (!hit) continue;

    Int_t seg = hit->SegmentId();
    if (seg < 0 || seg >= nSeg) continue;

    
    auto &wf_t = event.waveform_time[seg];
    auto &wf_v = event.waveform_volt[seg];

    const int nWave = hit->GetWaveformEntries(HodoRawHit::kUp);
    wf_t.reserve(nWave);
    wf_v.reserve(nWave);
    for (int k = 0; k < nWave; ++k) {
      auto wf = hit->GetWaveform(HodoRawHit::kUp, k);
      wf_t.push_back(wf.first);
      wf_v.push_back(wf.second);
    }

    
    int npulse = hit->GetNPulse(HodoRawHit::kUp);
    if (npulse > 0) {
      for (int k = 0; k < npulse; ++k) {
        const Double_t time   = hit->GetPulseTime  (HodoRawHit::kUp, k);
        const Double_t height = hit->GetPulseHeight(HodoRawHit::kUp, k);
        const Double_t dE     = hit->DeltaE(k);

        event.fit_time[seg].push_back(time);
        event.fit_height[seg].push_back(height);
        event.fit_dE[seg].push_back(dE);
      }
    }
    else {
      
      if (!wf_v.empty()) {
        auto it  = std::max_element(wf_v.begin(), wf_v.end());
        double maxval = *it;
        if (maxval > 1.0) {               
          size_t idx = std::distance(wf_v.begin(), it);
          double t0  = wf_t[idx];

          event.fit_time[seg].push_back(t0);
          event.fit_height[seg].push_back(maxval);
          event.fit_dE[seg].push_back(0.0); 
        }
      }
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  const Int_t    NbinFADC_X     = 500;
  const Int_t    NbinFADC_Y     = 1024;
  const Double_t MinIntegral    = 0;
  const Double_t MaxIntegral    = 10000;
  const Int_t    NbinIntegral   = (Int_t)(MaxIntegral - MinIntegral);
  const Double_t MinTDC         = 0.;
  const Double_t MaxTDC         = 3000000.;
  const Int_t    NbinTDC        = (Int_t)(MaxTDC - MinTDC);
  const Double_t MindE         = -0.5;
  const Double_t MaxdE         = 29.5;
  const Int_t    NbindE        = (Int_t)(MaxdE - MindE) * 100;

  HB1( 1, "Status",  20,   0., 20.);

  for (Int_t seg=0; seg<NumOfSegRayraw; ++seg) {
    TString title0  = Form("RAYRAW seg%d - Raw Waveform",      seg);
    TString title1  = Form("RAYRAW seg%d - Max ADC",           seg);
    TString title2  = Form("RAYRAW seg%d - Integral",          seg);
    TString title3  = Form("RAYRAW seg%d - Pedestal Integral", seg);
    TString title4  = Form("RAYRAW seg%d - TDC Leading",       seg);
    TString title5  = Form("RAYRAW seg%d - TDC Trailing",      seg);
    TString title6  = Form("RAYRAW seg%d - TDC First",         seg);
    TString title7  = Form("RAYRAW seg%d - TOT",               seg);
    TString title10 = Form("RAYRAW seg%d - dE(ADC)",           seg);
    TString title11 = Form("RAYRAW seg%d - dE(Integral)",      seg);

    HB2( RAYRAWHid + (seg+1)*1000 + 0,  title0,  NbinFADC_X, 0., (Double_t)NbinFADC_X, NbinFADC_Y, 0., (Double_t)NbinFADC_Y );
    HB1( RAYRAWHid + (seg+1)*1000 + 1,  title1,  NbinFADC_Y,   0.,          (Double_t)NbinFADC_Y );
    HB1( RAYRAWHid + (seg+1)*1000 + 2,  title2,  NbinIntegral, MinIntegral, MaxIntegral );
    HB1( RAYRAWHid + (seg+1)*1000 + 3,  title3,  NbinIntegral, MinIntegral, MaxIntegral );
    HB1( RAYRAWHid + (seg+1)*1000 + 4,  title4,  NbinTDC,      MinTDC,      MaxTDC );
    HB1( RAYRAWHid + (seg+1)*1000 + 5,  title5,  NbinTDC,      MinTDC,      MaxTDC );
    HB1( RAYRAWHid + (seg+1)*1000 + 6,  title6,  NbinTDC,      MinTDC,      MaxTDC );
    HB1( RAYRAWHid + (seg+1)*1000 + 7,  title7,  NbinTDC,      MinTDC,      MaxTDC );
    HB1( RAYRAWHid + (seg+1)*1000 + 10, title10, NbindE,       MindE,       MaxdE );
    HB1( RAYRAWHid + (seg+1)*1000 + 11, title11, NbindE,       MindE,       MaxdE );
  }

  // TTree
  HBTree( "tree", "tree" );
  tree->Branch("evnum",          &event.evnum,     "evnum/I");
  tree->Branch("fit_time",       &event.fit_time);
  tree->Branch("fit_height",     &event.fit_height);
  tree->Branch("fit_dE",         &event.fit_dE);
  tree->Branch("waveform_time",  &event.waveform_time);
  tree->Branch("waveform_volt",  &event.waveform_volt);

  HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")    &&
     InitializeParameter<HodoParamMan>("HDPRM") &&
     InitializeParameter<HodoPHCMan>("HDPHC")   &&
     InitializeParameter<UserParamMan>("USER")  &&
     InitializeParameter<TemplateFitMan>("RAYRAWTEMP")
     );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

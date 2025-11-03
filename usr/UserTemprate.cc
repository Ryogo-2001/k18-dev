// -*- C++ -*-

#include <iostream>
#include <sstream>
#include <cmath>
#include <numeric>
#include <algorithm>

#include "ConfMan.hh"
#include "RawData.hh"
#include "HodoAnalyzer.hh"
#include "HodoRawHit.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "UserParamMan.hh"
#include "DCGeomMan.hh"
#include "UnpackerManager.hh"   // GUnpacker

#include "RootHelper.hh"

#include "TROOT.h"
#include "TDirectory.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"

namespace
{
  using namespace root;
  using hddaq::unpacker::GUnpacker;

 
  const int    kSegToInspect      = 0;          
  const int    kMaxEventsToFill   = 100000000;  
  const double kRelativeThreshold = 60.0;       
  const int    kBaselineSamples   = 8;          

  
  const int    kXbins = 200;
  const double kXmin  = -50.0;
  const double kXmax  = 150.0;
  const int    kYbins = 140;
  const double kYmin  = -0.2;
  const double kYmax  =  1.2;

  
  const double kFitMin         = -20.0;
  const double kFitMax         =  80.0;
  const double kFixedAmplitude =  1.05;

  TH2F *hAligned     = nullptr;
  long  gFilledEvents = 0;

  
  double estimate_baseline_mean(const std::vector<Double_t> &wf, int n_samples)
  {
    if (wf.empty() || n_samples <= 0) return 0.0;
    int n = std::min<int>((int)wf.size(), n_samples);
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += wf[i];
    return s / n;
  }

} 


Bool_t
ProcessingBegin()
{
  return true;
}


Bool_t
ProcessingNormal()
{
  static const auto &gUnpacker  = GUnpacker::get_instance();

  if (gFilledEvents >= kMaxEventsToFill)
    return true;

  RawData rawData;
  rawData.DecodeHits("RAYRAW");

  const auto &cont = rawData.GetHodoRawHC("RAYRAW");
  if (cont.empty())
    return true;


  int evnum = gUnpacker.get_event_number();
  (void)evnum;

  HF1(1, 0);  // Status

  for (auto *hit : cont) {
    if (!hit) continue;
    int seg = hit->SegmentId();
    if (seg != kSegToInspect) continue;

    const auto &fadc = hit->GetArrayAdc();   // vector<Double_t>
    if ((int)fadc.size() <= kBaselineSamples) continue;

    // baseline
    double baseline = estimate_baseline_mean(fadc, kBaselineSamples);

    
    double max_val = -1e9;
    int    max_idx = -1;
    for (int i = kBaselineSamples; i < (int)fadc.size(); ++i) {
      double v = fadc[i];
      if (v > max_val) {
        max_val = v;
        max_idx = i;
      }
    }
    if (max_idx < 0) continue;

    double pulse_height = max_val - baseline;
    if (pulse_height < kRelativeThreshold)
      continue; 

    
    for (int i = 0; i < (int)fadc.size(); ++i) {
      double v_sub  = fadc[i] - baseline;
      double v_norm = (pulse_height > 1e-9) ? (v_sub / pulse_height) : 0.0;
      double rel_t  = (double)(i - max_idx);
      hAligned->Fill(rel_t, v_norm);
    }

    ++gFilledEvents;
  }

  return true;
}


Bool_t
ProcessingEnd()
{
  return true;
}


Bool_t
ConfMan::InitializeHistograms()
{
  HB1(1, "Status", 20, 0., 20.);

  
  hAligned = new TH2F("h_aligned_norm_2d",
                      Form("Aligned & Normalized Waveforms (seg%d, thr=%.1f);Relative Time from Peak (sample);Normalized Amplitude",
                           kSegToInspect, kRelativeThreshold),
                      kXbins, kXmin, kXmax,
                      kYbins, kYmin, kYmax);
  hAligned->SetStats(0);

  HPrint();
  return true;
}


Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")    &&
     InitializeParameter<HodoParamMan>("HDPRM") &&
     InitializeParameter<HodoPHCMan>("HDPHC")   &&
     InitializeParameter<UserParamMan>("USER"));
}


Bool_t
ConfMan::FinalizeProcess()
{
  if (!hAligned) {
    std::cerr << "[UserTemprate] no hAligned, skip." << std::endl;
    return true;
  }

  std::cout << "[UserTemprate] finalize: FitSlicesY ..." << std::endl;

  
  hAligned->FitSlicesY(0, 0, -1, 0, "QNR");

 

  TH1D *h_mean  = (TH1D*)gROOT->FindObject(Form("%s_1", hAligned->GetName()));
  TH1D *h_sigma = (TH1D*)gROOT->FindObject(Form("%s_2", hAligned->GetName()));
  TH1D *h_chi2  = (TH1D*)gROOT->FindObject(Form("%s_0", hAligned->GetName()));

  if (!h_mean) {
    std::cerr << "[UserTemprate] mean hist not found. (FitSlicesY failed?)" << std::endl;
    return true;
  }


  TF1 *func_conv = new TF1("rayraw_template_func_seg0",
    "[0]/2.0 * TMath::Exp( ([2]*[2]) / (2.0*[3]*[3]) - (x-[1])/[3] ) * (1.0 - TMath::Erf( ([2])/(sqrt(2.0)*[3]) - (x-[1])/(sqrt(2.0)*[2]) ))",
    kFitMin, kFitMax);

  func_conv->FixParameter(0, kFixedAmplitude); 
  func_conv->SetParameter(1, 0.0);  // t0
  func_conv->SetParameter(2, 1.5);  // sigma
  func_conv->SetParameter(3, 10.0); // tau
  func_conv->SetParName(0, "Amplitude(Fixed)");
  func_conv->SetParName(1, "t0");
  func_conv->SetParName(2, "sigma");
  func_conv->SetParName(3, "tau");
  func_conv->SetParLimits(2, 0.01, 10.0);
  func_conv->SetParLimits(3, 0.1, 100.0);

  h_mean->Fit(func_conv, "RQN");

 
  int nb = 400;  
  double xmin = kXmin;
  double xmax = kXmax;
  TH1D *h_template = new TH1D(Form("rayraw_template_seg%d", kSegToInspect),
                              Form("Template from conv-fit (seg%d);Relative Time (sample);Normalized Amp", kSegToInspect),
                              nb, xmin, xmax);
  for (int i = 1; i <= nb; ++i) {
    double x = h_template->GetBinCenter(i);
    double y = func_conv->Eval(x);
    h_template->SetBinContent(i, y);
  }

  
  std::cout << "[UserTemprate] objects to be written:" << std::endl;
  std::cout << "  - " << hAligned->GetName() << std::endl;
  std::cout << "  - " << h_mean->GetName() << std::endl;
  if (h_sigma) std::cout << "  - " << h_sigma->GetName() << std::endl;
  if (h_chi2)  std::cout << "  - " << h_chi2->GetName() << std::endl;
  std::cout << "  - " << h_template->GetName() << std::endl;
  std::cout << "  - " << func_conv->GetName() << " (TF1)" << std::endl;

  return true;
}

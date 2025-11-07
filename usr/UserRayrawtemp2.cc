// -*- C++ -*-
// #include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

#include "TMath.h"  // QuietNaN()

//#include "BH2Cluster.hh"
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
// #include "RAYRAWCalibMan.hh" // Not implemented

// ---------------------------------------------------------------------------------------
// makeHWF: 1 → Normalization for template "creation", 0 → Template "usage visualization"
// ---------------------------------------------------------------------------------------
#define makeHWF 1

namespace
{
using namespace root;
using hddaq::unpacker::GUnpacker;

const auto  qnan       = TMath::QuietNaN();
const auto& gUnpacker  = GUnpacker::get_instance();
const auto& gUser      = UserParamMan::GetInstance();
const auto& gHodo      = HodoParamMan::GetInstance();
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;

  std::vector<std::vector<Double_t>> waveform;       // [seg][sample] ADC
  std::vector<Double_t>              max_adc;        // [seg]
  std::vector<Double_t>              integral;       // [seg]
  std::vector<Double_t>              integral_ped;   // [seg]
  std::vector<std::vector<Double_t>> leading;        // [seg][k]
  std::vector<std::vector<Double_t>> trailing;       // [seg][k]
  std::vector<std::vector<Double_t>> leading_corr;   // [seg][k]
  std::vector<Double_t>              tdc_first;      // [seg]
  std::vector<Double_t>              tdc_first_corr; // [seg]
  std::vector<Double_t>              dE;             // [seg]
  std::vector<Double_t>              dE_q;           // [seg]

  Int_t Npulse[NumOfSegRayraw];                    

  void clear();
};

void Event::clear()
{
  evnum = 0;
  waveform.clear();
  max_adc.clear();
  integral.clear();
  integral_ped.clear();
  leading.clear();
  trailing.clear();
  leading_corr.clear();
  tdc_first.clear();
  tdc_first_corr.clear();
  dE.clear();
  dE_q.clear();

  for (Int_t i=0; i<NumOfSegRayraw; ++i) Npulse[i] = 0;
}

//_____________________________________________________________________________
namespace root
{
  Event  event;
  TH1   *h[MaxHist];
  TTree *tree;

  enum eDetHid { RAYRAWHid = 100000 };
}

//_____________________________________________________________________________
Bool_t ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t ProcessingNormal()
{
  
  static const auto MinRange         = gUser.GetParameter("RangeRAYRAW", 0);
  static const auto MaxRange         = gUser.GetParameter("RangeRAYRAW", 1);
  static const auto CorrectionTdcMin = gUser.GetParameter("TdcRAYRAW",   0);
  static const auto CorrectionTdcMax = gUser.GetParameter("TdcRAYRAW",   1);
  static const Int_t PedInitial      = 500;

  RawData rawData;
  rawData.DecodeHits("RAYRAW");         
  HodoAnalyzer hodoAna(rawData);

  event.evnum = GUnpacker::get_instance().get_event_number();
  HF1(1, 0);

  const auto& phcMan = HodoPHCMan::GetInstance();


  {
    const auto& cont = rawData.GetHodoRawHC("RAYRAW");
    const Int_t nh   = cont.size();

    
    event.waveform.resize(NumOfSegRayraw);
    event.leading.resize(NumOfSegRayraw);
    event.trailing.resize(NumOfSegRayraw);
    event.leading_corr.resize(NumOfSegRayraw);
    event.max_adc.resize(NumOfSegRayraw, qnan);
    event.integral.resize(NumOfSegRayraw, qnan);
    event.integral_ped.resize(NumOfSegRayraw, qnan);
    event.tdc_first.resize(NumOfSegRayraw, qnan);
    event.tdc_first_corr.resize(NumOfSegRayraw, qnan);
    event.dE.resize(NumOfSegRayraw, qnan);
    event.dE_q.resize(NumOfSegRayraw, qnan);

    for (Int_t i=0; i<nh; ++i) {
      HodoRawHit* hit = cont[i];
      const Int_t cid  = hit->DetectorId();
      const Int_t plid = hit->PlaneId();
      const Int_t seg  = hit->SegmentId();
      if (seg < 0 || seg >= NumOfSegRayraw) continue;

      const Double_t ped    = gHodo.GetP0(cid, plid, seg, 0);
      const Double_t gain   = gHodo.GetP1(cid, plid, seg, 0);
      const Double_t ped_q  = gHodo.GetP0(cid, plid, seg, 2);
      const Double_t gain_q = gHodo.GetP1(cid, plid, seg, 2);

      Int_t max_adc      = -10;
      Int_t integral     = -10;
      Int_t integral_ped = -10;
      Int_t tdc_first    = -10;
      Int_t nsample      = 0;

      const Int_t hid_wf           = RAYRAWHid + (seg+1)*1000 + 0;
      const Int_t hid_adc          = RAYRAWHid + (seg+1)*1000 + 1;
      const Int_t hid_integral     = RAYRAWHid + (seg+1)*1000 + 2;
      const Int_t hid_integral_ped = RAYRAWHid + (seg+1)*1000 + 3;
      const Int_t hid_tdc_l        = RAYRAWHid + (seg+1)*1000 + 4;
      const Int_t hid_tdc_t        = RAYRAWHid + (seg+1)*1000 + 5;
      const Int_t hid_tdc_first    = RAYRAWHid + (seg+1)*1000 + 6;
      const Int_t hid_dE           = RAYRAWHid + (seg+1)*1000 + 10;
      const Int_t hid_dE_q         = RAYRAWHid + (seg+1)*1000 + 11;

      for (const auto& fadc : hit->GetArrayAdc()) {
        HF2(hid_wf, nsample, fadc);
        event.waveform[seg].push_back(fadc);
        if (MinRange <= nsample && nsample <= MaxRange) {
          integral += fadc - PedInitial;
          if (fadc > max_adc) max_adc = fadc;
        }
        if (nsample < MinRange) integral_ped += fadc - PedInitial;
        ++nsample;
      }

      HF1(hid_adc, max_adc);
      event.max_adc[seg] = max_adc;

      HF1(hid_integral, integral);
      event.integral[seg] = integral;

      HF1(hid_integral_ped, integral_ped);
      event.integral_ped[seg] = integral_ped;

      const Double_t dE   = ( (Double_t)max_adc - ped   ) / (gain   - ped);
      const Double_t dE_q = ( (Double_t)integral - ped_q) / (gain_q - ped_q);
      HF1(hid_dE,   dE);
      HF1(hid_dE_q, dE_q);
      event.dE[seg]   = dE;
      event.dE_q[seg] = dE_q;

      // Leading/Trailing TDC
      for (const auto& tdc_l : hit->GetArrayTdcLeading()) {
        HF1(hid_tdc_l, tdc_l);
        event.leading[seg].push_back(tdc_l);
        if (tdc_l > tdc_first) tdc_first = tdc_l;

        Double_t corrected_tdc = tdc_l;
        if (tdc_l >= CorrectionTdcMin && tdc_l <= CorrectionTdcMax) {
          Double_t out;
          if (HodoPHCMan::GetInstance().DoCorrection(cid, plid, seg, 0, tdc_l, integral, out)) {
            corrected_tdc = out;
          }
        }
        event.leading_corr[seg].push_back(corrected_tdc);
      }
      HF1(hid_tdc_first, tdc_first);
      event.tdc_first[seg] = tdc_first;

      Double_t tdc_first_corr = qnan;
      if (tdc_first != -10) {
        tdc_first_corr = tdc_first;
        if (tdc_first >= CorrectionTdcMin && tdc_first <= CorrectionTdcMax) {
          Double_t out;
          if (HodoPHCMan::GetInstance().DoCorrection(cid, plid, seg, 0, tdc_first, integral, out)) {
            tdc_first_corr = out;
          }
        }
      }
      event.tdc_first_corr[seg] = tdc_first_corr;

      for (const auto& tdc_t : hit->GetArrayTdcTrailing()) {
        HF1(hid_tdc_t, tdc_t);
        event.trailing[seg].push_back(tdc_t);
      }
    }
  }

  // ---------------------------------------------------------------------------------------------
  hodoAna.DecodeHits<RayrawWaveformHit>("RAYRAW");
  {
    const auto& U  = HodoRawHit::kUp;
    const Int_t nh = hodoAna.GetNHits("RAYRAW");

    for (Int_t i=0; i<nh; ++i) {
      const auto& hit = hodoAna.GetHit<RayrawWaveformHit>("RAYRAW", i);
      if (!hit) continue;

      const Int_t seg   = hit->SegmentId();
      const Int_t plane = hit->PlaneId(); (void)plane;
      if (seg < 0 || seg >= NumOfSegRayraw) continue;

      const Int_t hid_wf_raw       = RAYRAWHid + (seg+1)*1000 + 201; // HWF
      const Int_t hid_pulse_height = RAYRAWHid + (seg+1)*1000 + 202; // Fit height
      const Int_t hid_pulse_time   = RAYRAWHid + (seg+1)*1000 + 203; // Fit time
      const Int_t hid_wf_fail      = RAYRAWHid + (seg+1)*1000 + 205; // Fit failure waveform
      const Int_t hid_wf_temp      = RAYRAWHid + (seg+1)*1000 + 206; // normalized template

      const Int_t NhitWF = hit->GetWaveformEntries(U);

      // Packing the HWF raw waveform into “time x ADC”（X:0..~数百ns, Y:0..ADC）
      for (Int_t m=0; m<NhitWF; ++m) {
        const auto wf = hit->GetWaveform(U, m); // (time[ns], adc)
        HF2(hid_wf_raw, wf.first, wf.second);
      }

      const Int_t Npulse = hit->GetNPulse(U);
      event.Npulse[seg]  = Npulse;

      Double_t pulse_height = -1.0;
      Double_t pulse_time   = -1.0;
      for (Int_t m=0; m<Npulse; ++m) {
        pulse_height = hit->GetPulseHeight(U, m);
        pulse_time   = hit->GetPulseTime(U, m);
        HF1(hid_pulse_height, pulse_height);
        HF1(hid_pulse_time,   pulse_time);
      }

      
      bool useThis = false;
    #if makeHWF
      // Template creation mode: Set the conditions as you like (here we use Npulse==1)
      useThis = (Npulse == 1);
    #else
      useThis = (Npulse == 1);
    #endif

      if (useThis && pulse_height > 0) {
      //pedestal subtraction
        double pede = 0.0; int NPede = 0;
        for (Int_t m=0; m<NhitWF; ++m) {
          const auto wf = hit->GetWaveform(U, m);
          if (wf.first < pulse_time - 5.0) { 
            pede += wf.second; ++NPede;
          }
        }
        if (NPede > 0) pede /= NPede;

      #if makeHWF
        for (Int_t m=0; m<NhitWF; ++m) {
          const auto wf = hit->GetWaveform(U, m);
          HF2(hid_wf_temp, wf.first - pulse_time, (wf.second - pede)/pulse_height);
        }
      #else
        // Fit height normalization 
        for (Int_t m=0; m<NhitWF; ++m) {
          const auto wf = hit->GetWaveform(U, m);
          HF2(hid_wf_temp, wf.first - pulse_time, (wf.second - pede)/pulse_height);
        }
      #endif
      }

      if (Npulse == 0) {
        for (Int_t m=0; m<NhitWF; ++m) {
          const auto wf = hit->GetWaveform(U, m);
          HF2(hid_wf_fail, wf.first, wf.second);
        }
      }
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
Bool_t ConfMan::InitializeHistograms()
{
  HB1(1, "Status", 20, 0., 20.);

  // Raw(HodoRawHit)
  const Int_t    NbinFADC_X   = 500;
  const Int_t    NbinFADC_Y   = 4096;
  const Double_t Xmin_wf_raw  = 0.0;
  const Double_t Xmax_wf_raw  = (Double_t)NbinFADC_X;
  const Double_t Ymin_wf_raw  = 0.0;              
  const Double_t Ymax_wf_raw  = (Double_t)NbinFADC_Y;

  // HWF range（time[ns] × ADC）
  const Int_t    NbinTime     = 500;
  const Double_t Tmin         = 0.0;
  const Double_t Tmax         = 250.0;             //  dt=0.5ns × 500 sample
  const Int_t    NbinADC      = 1024;
  const Double_t ADCmin       = 0.0;
  const Double_t ADCmax       = 1024.0;

  // Fit monitor 
  const Double_t MinPH = 0.;
  const Double_t MaxPH = 20000.;
  const Int_t    NbinPH = 4000;
  const Double_t MinPulseTime = -2.0;
  const Double_t MaxPulseTime =  2.0;
  const Int_t    NbinPulseTime = 400;

  for (Int_t seg=0; seg<NumOfSegRayraw; ++seg) {
    // ---- Raw(HodoRawHit) ----
    TString t0  = Form("RAYRAW seg%d - Raw Waveform (HodoRawHit)",      seg);
    TString t1  = Form("RAYRAW seg%d - Max ADC (HodoRawHit)",           seg);
    TString t2  = Form("RAYRAW seg%d - Integral (HodoRawHit)",          seg);
    TString t3  = Form("RAYRAW seg%d - Pedestal Integral (HodoRawHit)", seg);
    TString t4  = Form("RAYRAW seg%d - TDC Leading (HodoRawHit)",       seg);
    TString t5  = Form("RAYRAW seg%d - TDC Trailing (HodoRawHit)",      seg);
    TString t6  = Form("RAYRAW seg%d - TDC First (HodoRawHit)",         seg);
    TString t7  = Form("RAYRAW seg%d - TOT (HodoRawHit)",               seg);
    TString t10 = Form("RAYRAW seg%d - dE(ADC) (HodoRawHit)",           seg);
    TString t11 = Form("RAYRAW seg%d - dE(Integral) (HodoRawHit)",      seg);

    HB2(RAYRAWHid + (seg+1)*1000 + 0, t0,
        NbinFADC_X, Xmin_wf_raw, Xmax_wf_raw,
        NbinFADC_Y, Ymin_wf_raw, Ymax_wf_raw);
    HB1(RAYRAWHid + (seg+1)*1000 + 1,  t1,  NbinFADC_Y,   0., (Double_t)NbinFADC_Y);
    HB1(RAYRAWHid + (seg+1)*1000 + 2,  t2,  10000, 0., 10000.);
    HB1(RAYRAWHid + (seg+1)*1000 + 3,  t3,  10000, 0., 10000.);
    HB1(RAYRAWHid + (seg+1)*1000 + 4,  t4,  3000000, 0., 3000000.);
    HB1(RAYRAWHid + (seg+1)*1000 + 5,  t5,  3000000, 0., 3000000.);
    HB1(RAYRAWHid + (seg+1)*1000 + 6,  t6,  3000000, 0., 3000000.);
    HB1(RAYRAWHid + (seg+1)*1000 + 7,  t7,  3000000, 0., 3000000.);
    HB1(RAYRAWHid + (seg+1)*1000 + 10, t10, 3000, -0.5, 29.5);
    HB1(RAYRAWHid + (seg+1)*1000 + 11, t11, 3000, -0.5, 29.5);

    // ---- HWF(RayrawWaveformHit) ----
    TString w201 = Form("RAYRAW seg %d : Waveform (from HWF)", seg);
    TString w202 = Form("RAYRAW seg %d : Pulse Height (from Fit)", seg);
    TString w203 = Form("RAYRAW seg %d : Pulse Time (from Fit)", seg);
    TString w205 = Form("RAYRAW seg %d : Waveform (Failure at Pulse Search)", seg);
    TString w206 = Form("RAYRAW seg %d : Template Waveform (normalized)", seg);

    HB2(RAYRAWHid + (seg+1)*1000 + 201, w201,200, 0.0, 250.0,500, 0.0, 2000.0); 
    HB2(RAYRAWHid + (seg+1)*1000 + 205, w205, 200, 0.0, 250.0,500, 0.0, 2000.0); 
    HB1(RAYRAWHid + (seg+1)*1000 + 202, w202,NbinPulseTime,0,200);
    HB1(RAYRAWHid + (seg+1)*1000 + 203, w203, NbinPulseTime, MinPulseTime, MaxPulseTime);
    HB2(RAYRAWHid + (seg+1)*1000 + 206, w206,200, -50, 60, 500, -0.5, 3);
  }

  // tree
  HBTree("tree", "tree");
  tree->Branch("evnum",          &event.evnum,     "evnum/I");
  tree->Branch("waveform",       &event.waveform);
  tree->Branch("max_adc",        &event.max_adc);
  tree->Branch("integral",       &event.integral);
  tree->Branch("integral_ped",   &event.integral_ped);
  tree->Branch("leading",        &event.leading);
  tree->Branch("trailing",       &event.trailing);
  tree->Branch("leading_corr",   &event.leading_corr);
  tree->Branch("tdc_first",      &event.tdc_first);
  tree->Branch("tdc_first_corr", &event.tdc_first_corr);
  tree->Branch("dE(MaxADC)",     &event.dE);
  tree->Branch("dE(Integral)",   &event.dE_q);
  tree->Branch("Npulse",         event.Npulse,  Form("Npulse[%d]/I", NumOfSegRayraw));

  HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")    &&
     InitializeParameter<HodoParamMan>("HDPRM") &&
     InitializeParameter<HodoPHCMan>("HDPHC")   &&
     InitializeParameter<UserParamMan>("USER")  &&
     InitializeParameter<TemplateFitMan>("RAYRAWTEMP")
     // && InitializeParameter<RAYRAWCalibMan>("RAYRAWCALIB")
    );
}

//_____________________________________________________________________________
Bool_t ConfMan::FinalizeProcess()
{
  return true;
}

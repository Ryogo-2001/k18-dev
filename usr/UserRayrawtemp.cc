//#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>

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

// Added from UserGBO.cc logic
#include "RayrawWaveformHit.hh" 
#include "TemplateFitMan.hh"
// #include "RAYRAWCalibMan.hh" // <-- don't need 

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

  std::vector<std::vector<Double_t>> waveform;
  std::vector<Double_t>              max_adc;
  std::vector<Double_t>              integral;
  std::vector<Double_t>              integral_ped;
  std::vector<std::vector<Double_t>> leading;
  std::vector<std::vector<Double_t>> trailing;
  std::vector<std::vector<Double_t>> leading_corr;
  std::vector<Double_t>              tdc_first;
  std::vector<Double_t>              tdc_first_corr;
  std::vector<Double_t>              dE;
  std::vector<Double_t>              dE_q;

  //
  Int_t Npulse[NumOfSegRayraw]; // Use NumOfSegRayraw

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
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

  // template fit pulses
  for(Int_t it=0; it<NumOfSegRayraw; ++it){
    Npulse[it] = 0;
  }
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
  // Parameters from UserRayraw.cc
  static const auto MinRange    = gUser.GetParameter("RangeRAYRAW", 0);
  static const auto MaxRange    = gUser.GetParameter("RangeRAYRAW", 1);
  static const auto CorrectionTdcMin = gUser.GetParameter("TdcRAYRAW", 0);
  static const auto CorrectionTdcMax = gUser.GetParameter("TdcRAYRAW", 1);
  static const Int_t PedInitial = 500;

  RawData rawData;
  HodoAnalyzer hodoAna(rawData);

  event.evnum = gUnpacker.get_event_number();

  HF1(1, 0);

  const auto& phcMan = HodoPHCMan::GetInstance();

  rawData.DecodeHits("RAYRAW");

  //_____________________________________________________________________________
  {
    const auto& cont = rawData.GetHodoRawHC("RAYRAW");
    Int_t nh = cont.size();

    
    event.waveform.resize(nh);
    event.leading.resize(nh);
    event.trailing.resize(nh);
    event.leading_corr.resize(nh);

    for( int i=0; i<nh; ++i ){

      HodoRawHit *hit = cont[i];

      Int_t cid  = hit->DetectorId();
      Int_t plid = hit->PlaneId();
      Int_t seg  = hit->SegmentId();
      
      
      if(seg < 0 || seg >= NumOfSegRayraw) continue;

     
      if(event.waveform.size() <= seg) event.waveform.resize(NumOfSegRayraw);
      if(event.leading.size() <= seg) event.leading.resize(NumOfSegRayraw);
      if(event.trailing.size() <= seg) event.trailing.resize(NumOfSegRayraw);
      if(event.leading_corr.size() <= seg) event.leading_corr.resize(NumOfSegRayraw);
      if(event.max_adc.size() <= seg) event.max_adc.resize(NumOfSegRayraw, qnan);
      if(event.integral.size() <= seg) event.integral.resize(NumOfSegRayraw, qnan);
      if(event.integral_ped.size() <= seg) event.integral_ped.resize(NumOfSegRayraw, qnan);
      if(event.tdc_first.size() <= seg) event.tdc_first.resize(NumOfSegRayraw, qnan);
      if(event.tdc_first_corr.size() <= seg) event.tdc_first_corr.resize(NumOfSegRayraw, qnan);
      if(event.dE.size() <= seg) event.dE.resize(NumOfSegRayraw, qnan);
      if(event.dE_q.size() <= seg) event.dE_q.resize(NumOfSegRayraw, qnan);


      Double_t ped    = gHodo.GetP0(cid, plid, seg, 0);
      Double_t gain   = gHodo.GetP1(cid, plid, seg, 0);
      Double_t ped_q  = gHodo.GetP0(cid, plid, seg, 2);
      Double_t gain_q = gHodo.GetP1(cid, plid, seg, 2);

      Int_t max_adc       = -10;
      Int_t integral      = -10;
      Int_t integral_ped  = -10;
      Int_t tdc_first     = -10;
      Int_t nsample        = 0;

      Int_t hid_wf           = RAYRAWHid + (seg+1)*1000 + 0;
      Int_t hid_adc          = RAYRAWHid + (seg+1)*1000 + 1;
      Int_t hid_integral     = RAYRAWHid + (seg+1)*1000 + 2;
      Int_t hid_integral_ped = RAYRAWHid + (seg+1)*1000 + 3;
      Int_t hid_tdc_l        = RAYRAWHid + (seg+1)*1000 + 4;
      Int_t hid_tdc_t        = RAYRAWHid + (seg+1)*1000 + 5;
      Int_t hid_tdc_first    = RAYRAWHid + (seg+1)*1000 + 6;
      Int_t hid_dE           = RAYRAWHid + (seg+1)*1000 + 10;
      Int_t hid_dE_q         = RAYRAWHid + (seg+1)*1000 + 11;
      
      for(const auto& fadc : hit->GetArrayAdc()){
        HF2(hid_wf, nsample, fadc);
        event.waveform[seg].push_back(fadc);
        if(MinRange <= nsample && nsample <= MaxRange){
          integral += fadc - PedInitial;
          if (fadc > max_adc) max_adc = fadc;
        }
        if(nsample < MinRange){
          integral_ped += fadc - PedInitial;
        }
        ++nsample;
      }

      HF1(hid_adc, max_adc);
      event.max_adc[seg] = max_adc; //  [seg]
      HF1(hid_integral, integral);
      event.integral[seg] = integral; //  [seg]
      HF1(hid_integral_ped, integral_ped);
      event.integral_ped[seg] = integral_ped; //  [seg]

      Double_t dE = ((Double_t)max_adc - ped) / (gain - ped);
      HF1(hid_dE, dE);
      event.dE[seg] = dE; //  [seg]
      Double_t dE_q = ((Double_t)integral - ped_q) / (gain_q - ped_q);
      HF1(hid_dE_q, dE_q);
      event.dE_q[seg] = dE_q; //  [seg]

      // TDC Leading
      for (const auto& tdc_l : hit->GetArrayTdcLeading() ) {
        HF1(hid_tdc_l, tdc_l);
        event.leading[seg].push_back(tdc_l);
        if(tdc_l > tdc_first)
          tdc_first = tdc_l;

        Double_t corrected_tdc = tdc_l;
        if (tdc_l >= CorrectionTdcMin && tdc_l <= CorrectionTdcMax) {
          Double_t temp_corrected_tdc;
          bool success = phcMan.DoCorrection(cid, plid, seg, 0,
                                              tdc_l,
                                              integral,
                                              temp_corrected_tdc);
          if (success) {
            corrected_tdc = temp_corrected_tdc; 
          }
        }
        event.leading_corr[seg].push_back(corrected_tdc);
      }
      HF1(hid_tdc_first, tdc_first);
      event.tdc_first[seg] = tdc_first; //  [seg]

      Double_t tdc_first_corr = qnan;
      if (tdc_first != -10) { 
        tdc_first_corr = tdc_first; 
        if (tdc_first >= CorrectionTdcMin && tdc_first <= CorrectionTdcMax) {
          Double_t corrected_val;
          bool success = phcMan.DoCorrection(cid, plid, seg, 0,
                                              tdc_first,
                                              integral,
                                              corrected_val);
          if (success) {
            tdc_first_corr = corrected_val; 
          }
        }
      }
      event.tdc_first_corr[seg] = tdc_first_corr; //  [seg]

      // TDC Trailing
      for (const auto& tdc_t : hit->GetArrayTdcTrailing() ) {
        HF1(hid_tdc_t, tdc_t);
        event.trailing[seg].push_back(tdc_t);
      }
    } 
  }
  
  // ------------------------------------------------------------
  // RAYRAW Waveform Hit Analysis
  
  hodoAna.DecodeHits<RayrawWaveformHit>("RAYRAW"); 
  {
    const auto& U = HodoRawHit::kUp; // Use this as Channel (ch)
    Int_t nh=hodoAna.GetNHits("RAYRAW"); 
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit<RayrawWaveformHit>("RAYRAW", i); 
      if(!hit) continue;
      Int_t seg   = hit->SegmentId();
      Int_t plane = hit->PlaneId(); // plane

      
      Int_t hid_wf_raw       = RAYRAWHid + (seg+1)*1000 + 201; 
      Int_t hid_pulse_height = RAYRAWHid + (seg+1)*1000 + 202; 
      Int_t hid_pulse_time   = RAYRAWHid + (seg+1)*1000 + 203; 
      // Int_t hid_pulse_de     = RAYRAWHid + (seg+1)*1000 + 204; // do not use
      Int_t hid_wf_fail      = RAYRAWHid + (seg+1)*1000 + 205; 
      Int_t hid_wf_temp      = RAYRAWHid + (seg+1)*1000 + 206; 
      // Int_t hid_chi2_res     = RAYRAWHid + (seg+1)*1000 + 302; // do not use
      // Int_t hid_chi2_ph      = RAYRAWHid + (seg+1)*1000 + 303; // do not use

      Int_t NhitWF = hit->GetWaveformEntries(U);
      // Int_t NDiscri = hit->GetNDiscriPulse(); 
      // Int_t NDiffDiscri = hit->GetNDiscriDiffPulse(); 
      Double_t peak = 999;
      Double_t time = 999;
      Double_t pede = 0;
      for(Int_t m = 0; m<NhitWF; ++m){
        std::pair<Double_t, Double_t> waveform = hit->GetWaveform(U, m);
        HF2 (hid_wf_raw, waveform.first, waveform.second);
      }

      Int_t Npulse = hit->GetNPulse(U);
      Double_t pulse_height = -1;
      Double_t pulse_time = -1;
      // Double_t de = -1; 

      // Double_t chi2 = -999; 
      // Double_t max_res = -999; 

      event.Npulse[seg] = Npulse;
      for(Int_t m = 0; m<Npulse; ++m){
        pulse_height = hit->GetPulseHeight(U, m);
        pulse_time   = hit->GetPulseTime(U, m);
        // de           = hit->DeltaE(m); 
        
        // chi2         = hit->GetChi2();
        // max_res      = hit->GetMaxRes(); 

        HF1 (hid_pulse_height, pulse_height);
        HF1 (hid_pulse_time, pulse_time);
        // HF1 (hid_pulse_de, de); 
        // HF2 (hid_chi2_res, max_res, chi2); 
        // HF2 (hid_chi2_ph, pulse_height, chi2); 
      }
        /*
       if(chi2>=0 && max_res>=0){ 
         HF2 (RAYRAWHid + 1, seg, chi2);
         HF2 (RAYRAWHid + 2, seg, max_res);
       }
        */
      Int_t NPede = 0;
      Bool_t flag = false;
#if makeHWF
      // if(NDiscri==1 && NDiffDiscri==1) // <-- This logic is now broken
      //   flag = true;
#else
      if(Npulse==1)
        flag = true;
#endif
      if(flag){
        Int_t m0 = -1;
        for(Int_t m=0; m<NhitWF; m++){
          std::pair<Double_t, Double_t> waveform = hit->GetWaveform(U, m);
          if(waveform.second<peak){
            time = waveform.first;
            peak = waveform.second;
            m0 = m;
          }
        }

        Int_t m1 = 10;
        Double_t range = 50;
        
        
        if(seg==14 || seg==15) 
          m1 = 30;
        if(seg==0){
          m1 = 40;
          range = 150;
        }
        if(seg==1){
          m1 = 40;
          range = 100;
        }

        for(Int_t m=0; m<m0-2; m++){
          if(m>m0-m1){
            std::pair<Double_t, Double_t> height = hit->GetWaveform(U, m);
            if(std::abs(height.second) < range){
              pede += height.second;
              NPede++;
            }
          }
        }
        if(NPede>=1){
          pede /= NPede;

          if(pulse_height>1000){
            for(Int_t m=0; m<NhitWF; m++){
              std::pair<Double_t, Double_t> waveform = hit->GetWaveform(U, m);
#if makeTWF
              HF2 (hid_wf_temp, waveform.first - time, (waveform.second - pede)/(peak - pede));
#else
              HF2 (hid_wf_temp, waveform.first - time, (waveform.second - pede)/pulse_height);
#endif
            }
          }
        }
      }

      if (Npulse == 0) {
        for(Int_t m = 0; m<NhitWF; ++m){
          std::pair<Double_t, Double_t> waveform = hit->GetWaveform(U, m);
          HF2 (hid_wf_fail, waveform.first, waveform.second);
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
  //same UserRayraw.cc
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

  // waveform histograms
  const Double_t MinTime = -100.;
  const Double_t MaxTime = 100.;
  const Int_t    NbinTime = 4000;
  const Double_t MinPH = 0.;
  const Double_t MaxPH = 20000.;
  const Int_t    NbinPH = 4000;
  const Double_t MinPulseTime = -2.0;
  const Double_t MaxPulseTime = 2.0;
  const Int_t    NbinPulseTime = 400;
  const Double_t MinDE_fit = 0.0;
  const Double_t MaxDE_fit = 400.0;
  const Int_t    NbinDE_fit = 400;
  const Double_t MinChi2 = 0.;
  const Double_t MaxChi2 = 30.;
  const Int_t    NbinChi2 = 300;
  const Double_t MinRes = 0.;
  const Double_t MaxRes = 50.;
  const Int_t    NbinRes = 500;


  for (Int_t seg=0; seg<NumOfSegRayraw; ++seg) {
    // --- same UserRayraw.cc ---
    TString title0  = Form("RAYRAW seg%d - Raw Waveform (HodoRawHit)",      seg);
    TString title1  = Form("RAYRAW seg%d - Max ADC (HodoRawHit)",           seg);
    TString title2  = Form("RAYRAW seg%d - Integral (HodoRawHit)",          seg);
    TString title3  = Form("RAYRAW seg%d - Pedestal Integral (HodoRawHit)", seg);
    TString title4  = Form("RAYRAW seg%d - TDC Leading (HodoRawHit)",       seg);
    TString title5  = Form("RAYRAW seg%d - TDC Trailing (HodoRawHit)",      seg);
    TString title6  = Form("RAYRAW seg%d - TDC First (HodoRawHit)",         seg);
    TString title7  = Form("RAYRAW seg%d - TOT (HodoRawHit)",               seg);
    TString title10 = Form("RAYRAW seg%d - dE(ADC) (HodoRawHit)",           seg);
    TString title11 = Form("RAYRAW seg%d - dE(Integral) (HodoRawHit)",      seg);

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
    
    // --- WaveformFit ---
    TString title201 = Form("RAYRAW seg %d : Waveform (from HWF)", seg);
    TString title202 = Form("RAYRAW seg %d : Pulse Height (from Fit)", seg);
    TString title203 = Form("RAYRAW seg %d : Pulse Time (us) (from Fit)", seg);
    // TString title204 = Form("RAYRAW seg %d : dE (MeV) (from Fit)", seg); 
    TString title205 = Form("RAYRAW seg %d : Waveform (Failure at Pulse Search)", seg);
    TString title206 = Form("RAYRAW seg %d : Template Waveform", seg);
    // TString title302 = Form("RAYRAW seg %d : chi2 vs max_res", seg); 
    // TString title303 = Form("RAYRAW seg %d : chi2 vs pulse height", seg);

    HB2( RAYRAWHid + (seg+1)*1000 + 201, title201, 200, -5, 5, 500, -20000, 10000);
    HB1( RAYRAWHid + (seg+1)*1000 + 202, title202, NbinPH, MinPH, MaxPH);
    HB1( RAYRAWHid + (seg+1)*1000 + 203, title203, NbinPulseTime, MinPulseTime, MaxPulseTime);
    // HB1( RAYRAWHid + (seg+1)*1000 + 204, title204, NbinDE_fit, MinDE_fit, MaxDE_fit); 
    HB2( RAYRAWHid + (seg+1)*1000 + 205, title205, 200, -5, 5, 500, -20000, 10000);
    HB2( RAYRAWHid + (seg+1)*1000 + 206, title206, 200, -5, 5, 500, -20000, 10000);
    // HB2( RAYRAWHid + (seg+1)*1000 + 302, title302, NbinRes, MinRes, MaxRes, NbinChi2, MinChi2, MaxChi2); 
    // HB2( RAYRAWHid + (seg+1)*1000 + 303, title303, NbinPH, MinPH, MaxPH, NbinChi2, MinChi2, MaxChi2); 
  }

  // Summary Histograms from UserGBO.cc logic
  // TString title300 = Form("RAYRAW chi2 vs seg"); 
  // TString title301 = Form("RAYRAW max_res vs seg"); 
  // HB2( RAYRAWHid + 1, title300, NumOfSegRayraw, 0, NumOfSegRayraw, NbinChi2, MinChi2, MaxChi2); 
  // HB2( RAYRAWHid + 2, title301, NumOfSegRayraw, 0, NumOfSegRayraw, NbinRes, MinRes, MaxRes); 


  //Tree
  HBTree( "tree", "tree" );
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

  // add Npulse
  tree->Branch("Npulse",   event.Npulse,  Form("Npulse[%d]/I", NumOfSegRayraw));

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
     // Added from UserGBO.cc
     InitializeParameter<TemplateFitMan>("RAYRAWTEMP") 
     // && InitializeParameter<RAYRAWCalibMan>("RAYRAWCALIB") // <- don't need
     );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
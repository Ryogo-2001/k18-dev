// -*- C++ -*-
// Author: Rintaro Kurata

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#include "BH2Cluster.hh"
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

  std::vector<std::vector<Double_t>> waveform; // [seg][nsample]
  std::vector<Double_t>              max_adc;
  std::vector<Double_t>              integral;
  std::vector<Double_t>              integral_ped;
  std::vector<std::vector<Double_t>> leading;  // [seg][depth]
  std::vector<std::vector<Double_t>> trailing; // [seg][depth]
  std::vector<Double_t>              tdc_first;
  std::vector<Double_t>              dE;
  std::vector<Double_t>              dE_q;
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
  tdc_first.clear();
  dE.clear();
  dE_q.clear();
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
  static const auto MinRange    = gUser.GetParameter("RangeRAYRAW", 0);
  static const auto MaxRange    = gUser.GetParameter("RangeRAYRAW", 1);
  static const Int_t PedInitial = 500;
  // static const Double_t V_per_ch = 3.3/1024.; // [V]
  // static const Double_t T_per_ch = 13.33; // [ns]

  RawData rawData;
  rawData.DecodeHits();
  HodoAnalyzer hodoAna(rawData);

  event.evnum = gUnpacker.get_event_number();

  HF1(1, 0);

    {
    const auto& cont = rawData.GetHodoRawHC("RAYRAW");
    Int_t nh = cont.size();
    // Double_t true_hit = 0;

    event.waveform.resize(nh);
    event.leading.resize(nh);
    event.trailing.resize(nh);

    for( int i=0; i<nh; ++i ){

      HodoRawHit *hit = cont[i];


      Int_t cid  = hit->DetectorId();
      Int_t plid = hit->PlaneId();
      Int_t seg  = hit->SegmentId();

      Double_t ped    = gHodo.GetP0(cid, plid, seg, 0); // for MaxADC (HodoParam p0 HighGain)
      Double_t gain   = gHodo.GetP1(cid, plid, seg, 0); // for MaxADC (HodoParam p1 HighGain)
      Double_t ped_q  = gHodo.GetP0(cid, plid, seg, 2); // for Integral (HodoParam p0 LowhGain)
      Double_t gain_q = gHodo.GetP1(cid, plid, seg, 2); // for Integral (HodoParam p1 LowGain)

      // if(seg == 0 && i%1000 == 0)
      // 	std::cout << "Pedestal: " << ped << ", Gain: " << gain << std::endl;


      Int_t max_adc       = -10;
      Int_t integral      = -10;
      Int_t integral_ped  = -10;
      Int_t tdc_first     = -10;
      // Int_t leading_hit_in = 0;
      // Int_t trailing_hit_in = 0;
      // Int_t leading_hit_out = 0;
      // Int_t trailing_hit_out = 0;

      Int_t nsample        = 0;

      Int_t hid_wf           = RAYRAWHid + (seg+1)*1000 + 0; // Raw Waveform
      Int_t hid_adc          = RAYRAWHid + (seg+1)*1000 + 1; // Max ADC
      Int_t hid_integral     = RAYRAWHid + (seg+1)*1000 + 2; // Integral
      Int_t hid_integral_ped = RAYRAWHid + (seg+1)*1000 + 3; // Integral of pedestal (before signal region)
      Int_t hid_tdc_l        = RAYRAWHid + (seg+1)*1000 + 4; // TDC Leading
      Int_t hid_tdc_t        = RAYRAWHid + (seg+1)*1000 + 5; // TDC Trailing
      Int_t hid_tdc_first    = RAYRAWHid + (seg+1)*1000 + 6; // TDC first Leading
      // Int_t hid_tot       = RAYRAWHid + (seg+1)*1000 + 7; // TOT

      Int_t hid_dE           = RAYRAWHid + (seg+1)*1000 + 10; // maxADC
      Int_t hid_dE_q         = RAYRAWHid + (seg+1)*1000 + 11; // Integral

      // TDC Leading
      for (const auto& tdc_l : hit->GetArrayTdcLeading() ) {
	HF1(hid_tdc_l, tdc_l);
	event.leading[seg].push_back(tdc_l);
	if(tdc_l > tdc_first)
	  tdc_first = tdc_l;

	// if(tdc_l != 0 ){
	//   if(min_l < tdc_l && tdc_l < max_l){
	//     leading_hit_in += 1;
	//   }else{
	//     leading_hit_out += 1;
	//   }
	// }
      }
      HF1(hid_tdc_first, tdc_first);
      event.tdc_first.push_back(tdc_first);

      // TDC Trailing
      for (const auto& tdc_t : hit->GetArrayTdcTrailing() ) {
	HF1(hid_tdc_t, tdc_t);
	event.trailing[seg].push_back(tdc_t);

	// if(tdc_t != 0 ){
	//   if(min_t < tdc_t && tdc_t < max_t){
	//     trailing_hit_in += 1;
	//   }else{
	//     trailing_hit_out += 1;
	//   }
	// }
      }

      // ADC block
      for(const auto& fadc : hit->GetArrayAdc()){
	HF2(hid_wf, nsample, fadc);
	event.waveform[seg].push_back(fadc);
	// Max ADC
	if(MinRange <= nsample && nsample <= MaxRange){
	  integral += fadc - PedInitial; // Reduce QDC value
	  if (fadc > max_adc)
	    max_adc = fadc;
	}
	if(nsample < MinRange){
	  integral_ped += fadc - PedInitial; // Reduce QDC value
	}
	++nsample;
      }

      HF1(hid_adc, max_adc);
      event.max_adc.push_back(max_adc);

      HF1(hid_integral, integral);
      event.integral.push_back(integral);

      HF1(hid_integral_ped, integral_ped);
      event.integral_ped.push_back(integral_ped);


      // Calculate dE using HodoParam
      Double_t dE = ((Double_t)max_adc - ped) / (gain - ped);
      HF1(hid_dE, dE);
      event.dE.push_back(dE);
      Double_t dE_q = ((Double_t)integral - ped_q) / (gain_q - ped_q);
      HF1(hid_dE_q, dE_q);
      event.dE_q.push_back(dE_q);

    } // for nh

    return true;
  }
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
  //HRTDC
  const Double_t MaxTDC         = 3000000.;
  //RAYRAW Multi Hit TDC
  //const Double_t MaxTDC         = 4096.;
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

  //Tree
  HBTree( "tree", "tree" );
  tree->Branch("evnum",          &event.evnum,     "evnum/I");
  tree->Branch("waveform",       &event.waveform);
  tree->Branch("max_adc",        &event.max_adc);
  tree->Branch("integral",       &event.integral);
  tree->Branch("integral_ped",   &event.integral_ped);
  tree->Branch("leading",        &event.leading);
  tree->Branch("trailing",       &event.trailing);
  tree->Branch("tdc_first",      &event.tdc_first);
  tree->Branch("dE(MaxADC)",     &event.dE);
  tree->Branch("dE(Integral)",   &event.dE_q);

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
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

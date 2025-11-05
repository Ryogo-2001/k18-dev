// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iostream>
#include <sstream>

#include <UnpackerManager.hh>

// #include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "HodoHit.hh"
#include "HodoWaveformHit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "HypsLib.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "DCGeomMan.hh"
#include "TemplateFitMan.hh"
#include "BGOCalibMan.hh"

// #define TimeCut    1 // in cluster analysis
#define FHitBranch 0 // make FiberHit branches (becomes heavy)
#define HodoHitPos 0
#define makeHWF    0 // 1:make TWF, 0:use HWF

namespace
{
using namespace root;
const auto qnan = TMath::QuietNaN();
auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
auto& gRM       = RMAnalyzer::GetInstance();
auto& gUser     = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;
  Int_t spill;

  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  Int_t Npulse[NumOfSegBGO];

  ////////// Normalized

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum     = 0;
  spill     = 0;

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<NumOfSegBGO; ++it){
    Npulse[it] = 0;
  }
}

//_____________________________________________________________________________
struct Dst
{
  Int_t evnum;
  Int_t spill;

  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  Int_t    nhBh2;
  Int_t    csBh2[NumOfSegBH2*MaxDepth];
  Double_t Bh2Seg[NumOfSegBH2*MaxDepth];
  Double_t tBh2[NumOfSegBH2*MaxDepth];
  Double_t t0Bh2[NumOfSegBH2*MaxDepth];
  Double_t dtBh2[NumOfSegBH2*MaxDepth];
  Double_t deBh2[NumOfSegBH2*MaxDepth];
  Double_t posBh2[NumOfSegBH2*MaxDepth];

  Int_t    nhT0;
  Int_t    csT0[NumOfSegT0*MaxDepth];
  Double_t T0Seg[NumOfSegT0*MaxDepth];
  Double_t tT0[NumOfSegT0*MaxDepth];

  // for HodoParam
  Double_t t0la[NumOfSegT0];
  Double_t t0lt[NumOfSegT0][MaxDepth];
  Double_t t0ra[NumOfSegT0];
  Double_t t0rt[NumOfSegT0][MaxDepth];
  Double_t ltT0Seg[NumOfSegT0][MaxDepth];
  Double_t rtT0Seg[NumOfSegT0][MaxDepth];

  Int_t    nhSac;
  Int_t    csSac[NumOfSegSAC*MaxDepth];
  Double_t SacSeg[NumOfSegSAC*MaxDepth];
  Double_t tSac[NumOfSegSAC*MaxDepth];
  Double_t deSac[NumOfSegSAC*MaxDepth];

  // Time0
  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;

  Int_t    nhTof;
  Int_t    csTof[NumOfSegTOF*MaxDepth];
  Double_t TofSeg[NumOfSegTOF*MaxDepth];
  Double_t tTof[NumOfSegTOF*MaxDepth];
  Double_t dtTof[NumOfSegTOF*MaxDepth];
  Double_t deTof[NumOfSegTOF*MaxDepth];

  // for HodoParam
  Double_t tofua[NumOfSegTOF];
  Double_t tofut[NumOfSegTOF][MaxDepth];
  Double_t tofda[NumOfSegTOF];
  Double_t tofdt[NumOfSegTOF][MaxDepth];
  Double_t utTofSeg[NumOfSegTOF][MaxDepth];
  Double_t dtTofSeg[NumOfSegTOF][MaxDepth];
  Double_t udeTofSeg[NumOfSegTOF];
  Double_t ddeTofSeg[NumOfSegTOF];

  void clear();
};

//_____________________________________________________________________________
void
Dst::clear()
{
  nhT0     = 0;
  nhSac    = 0;
  nhBh2    = 0;
  nhTof    = 0;
  evnum    = 0;
  spill    = 0;
  Time0Seg = qnan;
  deTime0  = qnan;
  Time0    = qnan;
  CTime0   = qnan;

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<NumOfSegBH2; ++it){
    for(Int_t m=0; m<MaxDepth; ++m){
      csBh2[MaxDepth*it + m]  = 0;
      Bh2Seg[MaxDepth*it + m] = qnan;
      tBh2[MaxDepth*it + m]   = qnan;
      t0Bh2[MaxDepth*it + m]  = qnan;
      dtBh2[MaxDepth*it + m]  = qnan;
      deBh2[MaxDepth*it + m]  = qnan;
      posBh2[MaxDepth*it + m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegT0; it++){
    t0la[it]  = qnan;
    t0ra[it]  = qnan;
    T0Seg[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      t0lt[it][m]        = qnan;
      t0rt[it][m]        = qnan;
      tT0[MaxDepth*it+m] = qnan;
      ltT0Seg[it][m]     = qnan;
      rtT0Seg[it][m]     = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegSAC; it++){
    SacSeg[it] = qnan;
    deSac[it]  = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      tSac[MaxDepth*it+m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegTOF; it++){
    tofua[it]     = qnan;
    tofda[it]     = qnan;
    udeTofSeg[it] = qnan;
    ddeTofSeg[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      tofut[it][m]            = qnan;
      tofdt[it][m]            = qnan;
      utTofSeg[it][m]         = qnan;
      dtTofSeg[it][m]         = qnan;
      csTof[MaxDepth*it + m]  = 0;
      TofSeg[MaxDepth*it + m] = qnan;
      tTof[MaxDepth*it + m]   = qnan;
      dtTof[MaxDepth*it + m]  = qnan;
      deTof[MaxDepth*it + m]  = qnan;
    }
  }

}


//_____________________________________________________________________________
namespace root
{
Event  event;
Dst    dst;
TH1*   h[MaxHist];
TTree* tree;
// TTree* hodo;
enum eDetHid {
  BH1Hid  = 10000,
  BH2Hid  = 20000,
  T0Hid   = 40000,
  SACHid  = 30000,
  TOFHid  = 60000,
};
}

//_____________________________________________________________________________
Bool_t
ProcessingBegin()
{
  event.clear();
  dst.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingNormal()
{
  static const auto MinTdcCFT = gUser.GetParameter("TdcCFT", 0);
  static const auto MaxTdcCFT = gUser.GetParameter("TdcCFT", 1);
  //static const auto MinTdcT0  = gUser.GetParameter("TdcT0", 0);
  //static const auto MaxTdcT0  = gUser.GetParameter("TdcT0", 1);
  //static const auto MinTdcSAC = gUser.GetParameter("TdcSAC", 0);
  //static const auto MaxTdcSAC = gUser.GetParameter("TdcSAC", 1);
  //static const auto MinTdcTOF = gUser.GetParameter("TdcTOF", 0);
  //static const auto MaxTdcTOF = gUser.GetParameter("TdcTOF", 1);
#if HodoHitPos
  static const auto PropVelBH2 = gUser.GetParameter("PropagationBH2");
#endif

  RawData rawData;
  HodoAnalyzer hodoAna(rawData);


  gRM.Decode();

  // event.evnum = gRM.EventNumber();
  event.evnum = gUnpacker.get_event_number();
  event.spill = gRM.SpillNumber();
  dst.evnum   = gRM.EventNumber();
  dst.spill   = gRM.SpillNumber();

  HF1(1, 0);


  //*****************  RawData  *****************

  // Trigger Flag
  rawData.DecodeHits("TFlag");
  std::bitset<NumOfSegTrig> trigger_flag;
  for(const auto& hit: rawData.GetHodoRawHC("TFlag")){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc();
    if(tdc > 0){
      event.trigpat[trigger_flag.count()] = seg;
      event.trigflag[seg] = tdc;
      dst.trigpat[trigger_flag.count()] = seg;
      dst.trigflag[seg] = tdc;
      trigger_flag.set(seg);
      HF1(10, seg);
      HF1(10+seg, tdc);
    }
  }

  if(trigger_flag[trigger::kSpillOnEnd] || trigger_flag[trigger::kSpillOffEnd])
    return true;

  HF1(1, 1);

  rawData.DecodeHits("BGO");
  // BGO
  {
    const auto& cont = rawData.GetHodoRawHC("BGO");
    const auto& U = HodoRawHit::kUp;
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t plane = hit->PlaneId();
      Int_t seg   = hit->SegmentId();

      Int_t NhitA = hit->GetSizeAdcHigh();
      Int_t NhitT = hit->GetSizeTdcUp();

      bool flag_tdc = false;
      for(Int_t m = 0; m<NhitT; ++m){
	Int_t tdc_l = hit->GetTdcLeading(U, m);
	HF1 (1000*(seg+1)+100, tdc_l);  //TDC Nhits

	//if (tdc_l>MinTdcCFT && tdc_l < MaxTdcCFT) {
	//flag_tdc = true;
	//event.cfttdc[plane][seg] = tdc_l;
	//}
      }
      //if (flag_tdc) {
      //HF1 (1000*(plane+1)+102, seg);
      //}


      //FADC
      for(Int_t m = 0; m<NhitA; ++m){
	Int_t adc = hit->GetAdcUp(m);
	HF2 (1000*(seg+1)+101, m, adc);

	//if (flag_tdc) {
	//HF2 (1000*(plane+1)+106, seg, adcH);
	//}
      }
    }
  }

  hodoAna.DecodeHits<HodoWaveformHit>("BGO");
  {
    const auto& U = HodoRawHit::kUp;
    Int_t nh=hodoAna.GetNHits("BGO");
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit<HodoWaveformHit>("BGO", i);
      if(!hit) continue;
      Int_t seg   = hit->SegmentId();
      Int_t plane = hit->PlaneId();

      Int_t NhitT = hit->GetEntries(U);
      bool flagTime = false;
      for(Int_t m = 0; m<NhitT; ++m){
	Double_t time = hit->GetTUp(m);
	HF1 (1000*(seg+1)+200, time);
	//if (std::abs(time)<100)
	//flagTime = true;
      }

      Int_t NhitWF = hit->GetWaveformEntries(U);
      Int_t NDiscri = hit->GetNDiscriPulse();
      Int_t NDiffDiscri = hit->GetNDiscriDiffPulse();
      Double_t peak = 999;
      Double_t time = 999;
      Double_t pede = 0;
      for(Int_t m = 0; m<NhitWF; ++m){
	std::pair<Double_t, Double_t> waveform = hit->GetWaveform(m);
	HF2 (1000*(seg+1)+201, waveform.first, waveform.second);
      }


      Int_t Npulse = hit->GetNPulse(U);
      Double_t pulse_height = -1;
      Double_t pulse_time = -1;
      Double_t de = -1;

      Double_t chi2 = -999;
      Double_t max_res = -999;

      event.Npulse[seg] = Npulse;
      for(Int_t m = 0; m<Npulse; ++m){
	pulse_height = hit->GetPulseHeight(m);
	pulse_time   = hit->GetPulseTime(m);
	de           = hit->DeltaE(m);
	chi2         = hit->GetChi2();
	max_res      = hit->GetMaxRes();
	HF1 (1000*(seg+1)+202, pulse_height);
	HF1 (1000*(seg+1)+203, pulse_time);
	HF1 (1000*(seg+1)+204, de);
	HF2 (1000*(seg+1)+302, max_res, chi2);
	HF2 (1000*(seg+1)+303, pulse_height, chi2);
      }

      if(chi2>=0 && max_res>=0){
	HF2 (300, seg, chi2);
	HF2 (301, seg, max_res);
      }

      Int_t NPede = 0;
      Bool_t flag = false;
#if makeHWF
      if(NDiscri==1 && NDiffDiscri==1)
	flag = true;
#else
      if(Npulse==1)
	flag = true;
#endif
      if(flag){
	Int_t m0 = -1;
	for(Int_t m=0; m<NhitWF; m++){
	  std::pair<Double_t, Double_t> waveform = hit->GetWaveform(m);
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
	    std::pair<Double_t, Double_t> height = hit->GetWaveform(m);
	    if(std::abs(height.second) < range){
	      pede += height.second;
	      NPede++;
	    }
	  }
	}
	if(NPede>=1){
	  pede /= NPede;
	  if(seg==0 || seg==1)
	    std::cout << seg << " " << pede << " "  << NPede << std::endl;

	  // if(peak<-200){
	  if(pulse_height>1000){
	    for(Int_t m=0; m<NhitWF; m++){
	      std::pair<Double_t, Double_t> waveform = hit->GetWaveform(m);

#if makeTWF
	      HF2 (1000*(seg+1)+206, waveform.first - time, (waveform.second - pede)/(peak - pede));
#else
	      HF2 (1000*(seg+1)+206, waveform.first - time, (waveform.second - pede)/pulse_height);
#endif
	      if (seg==0||seg==1) {
		// std::cout << seg << " : " << waveform.first
		// 	      << " - " << pulse_time << ", "
		// 	      << waveform.second << "-" << pede
		// 	      << " / " << pulse_height << std::endl;
	      }
	    }
	  }
	}
      }

      if (Npulse == 0) {
	for(Int_t m = 0; m<NhitWF; ++m){
	  std::pair<Double_t, Double_t> waveform = hit->GetWaveform(m);
	  HF2 (1000*(seg+1)+205, waveform.first, waveform.second);
	}
      }
      //if (flagTime)
      //HF2 (1000*(plane+1)+206, seg, adccorHi);
    }
  }



  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingEnd()
{
  tree->Fill();
  // hodo->Fill();
  return true;
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
  HB1(10, "Trigger HitPat", NumOfSegTrig, 0., Double_t(NumOfSegTrig));
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1(10+i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000);
  }

  for(Int_t i=0; i<NumOfSegBGO; i++){
    TString title100 = Form("BGO seg %d : Tdc(Leading)", i);
    TString title101 = Form("BGO seg %d : FADC", i);
    TString title200 = Form("BGO seg %d : Time(Leading)", i);
    TString title201 = Form("BGO seg %d : Waveform", i);
    TString title202 = Form("BGO seg %d : Pulse Height", i);
    TString title203 = Form("BGO seg %d : Pulse Time (us)", i);
    TString title204 = Form("BGO seg %d : dE (MeV)", i);
    TString title205 = Form("BGO seg %d : Waveform (Failure at Pulse Search)", i);
    TString title206 = Form("BGO seg %d : Template Waveform", i);
    TString title302 = Form("BGO seg %d : chi2 vs max_res", i);
    TString title303 = Form("BGO seg %d : chi2 vs pulse height", i);
    HB1( 1000*(i+1)+100, title100, 4000, 0, 4000);
    HB2( 1000*(i+1)+101, title101, 200, 0, 400, 400, 0, 20000);
    HB1( 1000*(i+1)+200, title200, 4000, -100, 100);
    HB2( 1000*(i+1)+201, title201, 200, -5, 5, 500, -20000, 10000);
    HB1( 1000*(i+1)+202, title202, 4000, 0, 20000);
    HB1( 1000*(i+1)+203, title203, 400, -2.0, 2.0);
    HB1( 1000*(i+1)+204, title204, 400, 0.0, 400);
    HB2( 1000*(i+1)+205, title205, 200, -5, 5, 500, -20000, 10000);
    HB2( 1000*(i+1)+206, title206, 200, -5, 5, 500, -20000, 10000);
    HB2( 1000*(i+1)+302, title302, 500, 0, 50, 300, 0, 30);
    HB2( 1000*(i+1)+303, title303, 4000, 0, 20000, 300, 0, 30);
  }
  TString title300 = Form("BGO chi2 vs seg");
  TString title301 = Form("BGO max_res vs seg");
  HB2( 300, title300, 24, 0, 24, 300, 0, 30);
  HB2( 301, title301, 24, 0, 24, 500, 0, 50);


  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("spill",     &event.spill,     "spill/I");
  //Trig
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  tree->Branch("Npulse",   event.Npulse,  Form("Npulse[%d]/I", NumOfSegBGO));


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
     InitializeParameter<HodoPHCMan>("HDPHC")   &&
     InitializeParameter<UserParamMan>("USER")   &&
     InitializeParameter<TemplateFitMan>("BGOTEMP") &&
     InitializeParameter<BGOCalibMan>("BGOCALIB")
     );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

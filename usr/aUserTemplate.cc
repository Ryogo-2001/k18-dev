// -*- C++ -*-
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

// ROOT includes for template saving
#include "TFile.h"
#include "TH2.h"
#include "TProfile.h"

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
      ALIGNEDHid = 200000  // テンプレート作成用ID
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
  static const auto CorrectionTdcMin = gUser.GetParameter("TdcRAYRAW", 0);
  static const auto CorrectionTdcMax = gUser.GetParameter("TdcRAYRAW", 1);
  static const Int_t PedInitial = 500;

  RawData rawData;
  rawData.DecodeHits();
  HodoAnalyzer hodoAna(rawData);

  event.evnum = gUnpacker.get_event_number();

  HF1(1, 0);

  const auto& phcMan = HodoPHCMan::GetInstance();

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

      // =======================================================
      // ===== テンプレート作成ロジック (最終修正版) =====
      // =======================================================
      {
        // 0. 対応するヒストグラムを取得
        Int_t hid_aligned = ALIGNEDHid + (seg+1)*1000 + 0;
        if (!root::h[hid_aligned] || !root::h[hid_aligned]->InheritsFrom(TH2::Class())) {
            continue;
        }
        TH2D* h_aligned = (TH2D*)root::h[hid_aligned];
        
        // 1. このイベントの生波形を取得
        const std::vector<Double_t>& raw_wf = event.waveform[seg];
        if(raw_wf.empty()) continue; 

        // 2a. ベースライン計算 (立ち上がり(x=15)より手前の 10 サンプルで計算)
        const Int_t ped_samples = 5;
        if(raw_wf.size() <= ped_samples) continue; 

        double baseline = 0.;
        for(int k=0; k<ped_samples; ++k) {
            baseline += raw_wf[k];
        }
        baseline /= ped_samples; 

        // 2b. ベースライン補正済みの波形(wf)を作成
        std::vector<Double_t> wf(raw_wf.size());
        double max_val = -1.e9; // ベースライン補正後の最大値
        int    t_peak  = -1;    // ピークのサンプル時間
        
        for(size_t k=0; k<wf.size(); ++k) {
            wf[k] = raw_wf[k] - baseline;
            
            // ピーク探索 (ベースライン領域 *直後* から探索を開始する)
            if(k >= ped_samples && wf[k] > max_val) {
                max_val = wf[k];
                t_peak  = k; 
            }
        }

        // 2c. "クリーンな波形" の選別
        if(t_peak == -1) continue; 
        if(max_val + baseline >= 1023.0) continue; 
        
        const double noise_threshold = 50.0; // (要調整)
        if(max_val < noise_threshold) continue;
        
        // 3. アライメントして TH2D に充填
        for(size_t k=0; k<wf.size(); ++k) {
            double rel_time = k - t_peak; 
            h_aligned->Fill(rel_time, wf[k]);
        }
      }
      // =======================================================
      // ===== 波形テンプレートロジックここまで =====
      // =======================================================


      HF1(hid_adc, max_adc);
      event.max_adc.push_back(max_adc);
      HF1(hid_integral, integral);
      event.integral.push_back(integral);
      HF1(hid_integral_ped, integral_ped);
      event.integral_ped.push_back(integral_ped);

      Double_t dE = ((Double_t)max_adc - ped) / (gain - ped);
      HF1(hid_dE, dE);
      event.dE.push_back(dE);
      Double_t dE_q = ((Double_t)integral - ped_q) / (gain_q - ped_q);
      HF1(hid_dE_q, dE_q);
      event.dE_q.push_back(dE_q);

      // =======================================================
      // ★★★ バグ修正 ★★★
      // TDC Leading (元のロジックに戻す)
      // =======================================================
      for (const auto& tdc_l : hit->GetArrayTdcLeading() ) {
        HF1(hid_tdc_l, tdc_l);
        event.leading[seg].push_back(tdc_l);
        if(tdc_l > tdc_first)
          tdc_first = tdc_l;

        Double_t corrected_tdc = tdc_l; // <-- ★ここを 't_peak' から 'tdc_l' に戻しました
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
      // =======================================================
      // ★★★ 修正ここまで ★★★
      // =======================================================

      HF1(hid_tdc_first, tdc_first);
      event.tdc_first.push_back(tdc_first);

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
      event.tdc_first_corr.push_back(tdc_first_corr);

      // TDC Trailing
      for (const auto& tdc_t : hit->GetArrayTdcTrailing() ) {
        HF1(hid_tdc_t, tdc_t);
        event.trailing[seg].push_back(tdc_t);
      }
    } 

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

    // --- テンプレート作成用のヒストグラム (★修正点3) ---
    Int_t n_pre = 50;  // (x=20がピークなので 50 あれば x=15 の立ち上がりも十分入る)
    Int_t n_post = 100; // (x=50でテールが終わるので 100 あれば十分)
    Int_t n_total = n_pre + n_post;
    
    TString title_aligned = Form("RAYRAW seg%d - Aligned Waveforms", seg);
    Int_t hid_aligned = ALIGNEDHid + (seg+1)*1000 + 0;
    
    // (Y軸の範囲も、生ベースライン500を引くことを考慮し -100 ~ 600 程度に)
    HB2( hid_aligned, title_aligned, n_total, -n_pre, n_post, // X軸: ピークからの相対時間
         700, -100., 600.);                              // Y軸: ベースライン補正後のADC
    // --- 追加ここまで ---
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
  tree->Branch("leading_corr",   &event.leading_corr);
  tree->Branch("tdc_first",      &event.tdc_first);
  tree->Branch("tdc_first_corr", &event.tdc_first_corr);
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
  // =======================================================
  // ===== テンプレートをファイルに保存する処理 (SegV修正版) =====
  // =======================================================
  
  TFile* template_file = new TFile("waveform_template.root", "RECREATE");
  if(!template_file || template_file->IsZombie()) {
      std::cerr << "Error: Could not create waveform_template.root" << std::endl;
  } else {
      std::cout << "Saving waveform templates to waveform_template.root ..." << std::endl;

      for (Int_t seg=0; seg<NumOfSegRayraw; ++seg) {
          Int_t hid_aligned = ALIGNEDHid + (seg+1)*1000 + 0;
          TH1* h1 = root::h[hid_aligned];
          if(!h1) continue; 
          
          if (!h1->InheritsFrom(TH2::Class())) continue; 
          TH2D* h_aligned = (TH2D*)h1;

          if(h_aligned->GetEntries() == 0) {
              std::cout << "Seg " << seg << ": Aligned histogram is empty. (GOOD EVENTS WERE 0)" << std::endl;
              continue; 
          }

          TProfile* prof_raw = h_aligned->ProfileX();
          if(!prof_raw) continue;
          
          prof_raw->SetDirectory(nullptr); 

          TString prof_name = Form("F_template_seg%d", seg);
          prof_raw->SetName(prof_name);

          // 4. 最大値で正規化 (Normalize)
          double max_val = prof_raw->GetMaximum();
          if(max_val > 0) {
              prof_raw->Scale(1.0 / max_val);
          }

          template_file->cd();
          prof_raw->Write();
          
          // (デバッグ用にアライメント済みTH2Dも保存)
          h_aligned->SetDirectory(nullptr); 
          h_aligned->Write(); 

          delete prof_raw;
      }

      template_file->Close();
      delete template_file; 
      
      std::cout << "Waveform template generation complete." << std::endl;
  }
  // =======================================================
  // ===== 追加処理ここまで =====
  // =======================================================

  return true;
}
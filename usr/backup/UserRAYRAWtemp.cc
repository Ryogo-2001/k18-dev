// -*- C++ -*-
#include "VEvent.hh"

#include <cmath>
#include <iostream>
#include <sstream>
#include <bitset>
#include <algorithm> // std::max

#include "BH2Hit.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "RootHelper.hh"
#include "HodoHit.hh"
#include "HodoWaveformHit.hh"
#include "HodoAnalyzer.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "RawData.hh"
#include "UserParamMan.hh"
#include "UnpackerManager.hh"
// #include "TemplateFitMan.hh"  // ← テンプレートを“使う”解析(makeHWF=0)に進めるときに有効化

#define makeHWF 0  // 1: make Template Waveform (TWF), 0: use existing HWF

namespace
{
using namespace root;
using hddaq::unpacker::GUnpacker;

const auto  qnan       = TMath::QuietNaN();
auto&       gUnpacker  = GUnpacker::get_instance();
auto&       gRM        = RMAnalyzer::GetInstance();
auto&       gUser      = UserParamMan::GetInstance();
auto&       gHodo      = HodoParamMan::GetInstance();
auto&       gPHC       = HodoPHCMan::GetInstance();
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
    };
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
  // いまは未使用（必要になったら復活してOK）
  // static const auto MinRange          = gUser.GetParameter("RangeRAYRAW", 0);
  // static const auto MaxRange          = gUser.GetParameter("RangeRAYRAW", 1);
  // static const auto CorrectionTdcMin  = gUser.GetParameter("TdcRAYRAW", 0);
  // static const auto CorrectionTdcMax  = gUser.GetParameter("TdcRAYRAW", 1);

  RawData rawData;
  rawData.DecodeHits();                            // ★ 全検出器まとめて1回だけデコード
  HodoAnalyzer hodoAna(rawData);
  hodoAna.DecodeHits<HodoWaveformHit>("RAYRAW");   // RAYRAWの波形ヒットを生成

  event.evnum = gUnpacker.get_event_number();

  HF1(1, 0);

  // ------------------------------
  // RAYRAW (FADC Waveform)
  // ------------------------------
  const auto& U = HodoRawHit::kUp;
  Int_t nh = hodoAna.GetNHits("RAYRAW");

  for (Int_t i = 0; i < nh; ++i) {
    const auto& hit = hodoAna.GetHit<HodoWaveformHit>("RAYRAW", i);
    if (!hit) continue;

    const Int_t seg   = hit->SegmentId();
    const Int_t plane = hit->PlaneId(); // 使うなら活用
    (void)plane;

    // セグメント範囲ガード（安全策）
    if (seg < 0 || seg >= NumOfSegRayraw) continue;

    // ヒストID
    const Int_t hid_wf  = RAYRAWHid + (seg+1)*1000 + 200; // Waveform (t, V)
    const Int_t hid_ph  = RAYRAWHid + (seg+1)*1000 + 202; // Pulse Height
    const Int_t hid_pt  = RAYRAWHid + (seg+1)*1000 + 203; // Pulse Time
    const Int_t hid_de  = RAYRAWHid + (seg+1)*1000 + 204; // dE (from fit)
    const Int_t hid_twf = RAYRAWHid + (seg+1)*1000 + 206; // Template Waveform (normalized)

    // まずは波形そのものを保存
    const Int_t NhitWF = hit->GetWaveformEntries(U);
    for (Int_t m=0; m<NhitWF; ++m) {
      const auto wf = hit->GetWaveform(m); // (time, voltage)
      HF2(hid_wf, wf.first, wf.second);
    }

    // パルス情報（テンプレートフィット結果: makeHWF=0 で有効になる）
    const Int_t Npulse = hit->GetNPulse(U);
    event.Npulse[seg] = Npulse;

    for (Int_t m = 0; m < Npulse; ++m) {
      const Double_t pulse_height = hit->GetPulseHeight(m);
      const Double_t pulse_time   = hit->GetPulseTime(m);
      const Double_t de           = hit->DeltaE(m);
      HF1(hid_ph, pulse_height);
      HF1(hid_pt, pulse_time);
      HF1(hid_de, de);
    }

    // ---- テンプレート波形（正規化）の生成：BGOと同じ分岐 ----
#if makeHWF
    // テンプレート“作成”モード：ディスクリ1本のクリーン波形のみ採用
    const int NDiscri     = hit->GetNDiscriPulse();
    const int NDiffDiscri = hit->GetNDiscriDiffPulse();
    const bool useThis = (NDiscri==1 && NDiffDiscri==1);
#else
    // テンプレート“利用”モード：1パルスイベントの可視化（任意だがBGO準拠）
    const bool useThis = (Npulse==1);
#endif

    if (useThis) {
      // ピーク（最小値）と時刻の推定
      double peak = -1e9;
      double time = -1e9;
      int    m0   = -1;
    for (Int_t m=0; m<NhitWF; ++m) {
      const auto wf = hit->GetWaveform(m);
  if (wf.second > peak) {  // ★正パルス → 最大値をピークとみなす
    peak = wf.second;
    time = wf.first;
    m0 = m;
  }
}     

      // ペデスタル推定（ピークより前の一定窓）
      double pede = 0.0; int NPede = 0;
      const int m1   = 10; // 窓幅はBGO実装に倣い適度に
      const int mend = std::max(0, m0 - m1);
      for (Int_t m=0; m<mend; ++m) {
        const auto wf = hit->GetWaveform(m);
        pede += wf.second; ++NPede;
      }
      if (NPede>0) pede /= NPede;

#if makeHWF
      // テンプレート作成時： (peak - pede) で正規化（高さ=1）
      const double amp = (peak - pede);
      if (amp == 0) continue;
      for (Int_t m=0; m<NhitWF; ++m) {
        const auto wf = hit->GetWaveform(m);
        HF2(hid_twf, wf.first - time, (wf.second - pede)/amp);
      }
#else
      // 通常解析時：フィット済みの pulse_height で正規化
      if (Npulse >= 1) {
        const double pulse_height = hit->GetPulseHeight(0);
        if (pulse_height <= 0) continue;
        for (Int_t m=0; m<NhitWF; ++m) {
          const auto wf = hit->GetWaveform(m);
          HF2(hid_twf, wf.first - time, (wf.second - pede)/pulse_height);
        }
      }
#endif
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
  HB1(1, "Status", 20, 0., 20.);

  for (Int_t seg=0; seg<NumOfSegRayraw; ++seg) {
    TString title200 = Form("RAYRAW seg %d : Waveform", seg);
    TString title202 = Form("RAYRAW seg %d : Pulse Height", seg);
    TString title203 = Form("RAYRAW seg %d : Pulse Time", seg);
    TString title204 = Form("RAYRAW seg %d : dE (MeV)", seg);
    TString title206 = Form("RAYRAW seg %d : Template Waveform", seg);

    // ビンや範囲はBGO版と整合（必要なら調整）
    HB2(RAYRAWHid + (seg+1)*1000 + 200, title200, 200, -5, 5, 500, -20000, 10000);
    HB1(RAYRAWHid + (seg+1)*1000 + 202, title202, 4000, 0, 20000);
    HB1(RAYRAWHid + (seg+1)*1000 + 203, title203, 400, -2.0, 2.0);
    HB1(RAYRAWHid + (seg+1)*1000 + 204, title204, 400, 0.0, 400);
    HB2(RAYRAWHid + (seg+1)*1000 + 206, title206, 200, -5, 5, 500, -20000, 10000);
  }

  HBTree("tree", "tree of RAYRAW");
  tree->Branch("evnum",  &event.evnum,  "evnum/I");
  //tree->Branch("spill",  &event.spill,  "spill/I");
  tree->Branch("Npulse", event.Npulse,  Form("Npulse[%d]/I", NumOfSegRayraw));

  return true;
}

//_____________________________________________________________________________
Bool_t ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")   &&
     InitializeParameter<HodoParamMan>("HDPRM") &&
     InitializeParameter<HodoPHCMan>("HDPHC")   &&
     InitializeParameter<UserParamMan>("USER")
     // && InitializeParameter<TemplateFitMan>("RAYRAWTEMP") // ← テンプレート“利用”段階で有効化
     );
}

//_____________________________________________________________________________
Bool_t ConfMan::FinalizeProcess()
{
  return true;
}

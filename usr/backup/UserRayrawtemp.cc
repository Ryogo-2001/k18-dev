

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "HodoAnalyzer.hh"
#include "HodoRawHit.hh"
#include "HodoHit.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "UserParamMan.hh"
#include "TemplateFitMan.hh"
#include "RayrawWaveformHit.hh"

#include "DCGeomMan.hh"
#include "RMAnalyzer.hh"
#include "S2sLib.hh"
#include "UnpackerManager.hh"

#include "RootHelper.hh"

#include <TMath.h>
#include <TGraph.h>

namespace
{
  using namespace root;
  using hddaq::unpacker::GUnpacker;

  const auto  qnan       = TMath::QuietNaN();
  const auto& gUnpacker  = GUnpacker::get_instance();
  const auto& gUser      = UserParamMan::GetInstance();
  const auto& gHodo      = HodoParamMan::GetInstance();

  // 再構成に使う時間範囲と刻み
  const double kRecoTmin = 0.0;    // [ns]
  const double kRecoTmax = 70.0;   // [ns]
  const int    kOverSample = 10;   // (0.5ns を 0.05ns にしたいので 10倍)

  // テンプレートファイルの実体をキャッシュしておく配列
  // segごとに  param/RAYRAWTEMP/RAYRAWTEMP.<seg> を読む
  std::vector<TGraph*> gTemplateGraphs;
}

/*================================================================
 *  Event 構造体
 *===============================================================*/
struct Event
{
  Int_t evnum;

  // [seg][pulse] テンプレートフィットで出たパルス列
  std::vector<std::vector<Double_t>> fit_time;
  std::vector<std::vector<Double_t>> fit_height;

  // [seg][sample] 生波形 (ADCそのまま)
  std::vector<std::vector<Double_t>> waveform_time;
  std::vector<std::vector<Double_t>> waveform_adc;

  // [seg][sample] テンプレートを足し合わせて作った細かい波形
  // 時間軸は segment ごとに共通でなくてもいいので [seg][i]
  std::vector<std::vector<Double_t>> fitted_time;
  std::vector<std::vector<Double_t>> fitted_adc;

  void clear();
};

void Event::clear()
{
  evnum = 0;
  fit_time.clear();
  fit_height.clear();
  waveform_time.clear();
  waveform_adc.clear();
  fitted_time.clear();
  fitted_adc.clear();
}

/*================================================================
 *  ROOT 名前空間にグローバル実体を置く (既存Userと同じやり方)
 *===============================================================*/
namespace root
{
  Event  event;
  TH1   *h[MaxHist];
  TTree *tree;

  enum eDetHid {
    RAYRAWHid = 100000
  };
}

/*================================================================
 *  テンプレートファイルを読むヘルパー
 *===============================================================*/
namespace
{
  // テキストで置いてある param/RAYRAWTEMP/RAYRAWTEMP.<seg> を読む
  // フォーマットは "time  amp" の2列を想定 (あなたがmacroで使ってたやつと同じ)
  TGraph* LoadTemplateGraphForSeg(int seg)
  {
    // すでに読み込み済みならそれを返す
    if (seg < (int)gTemplateGraphs.size() && gTemplateGraphs[seg])
      return gTemplateGraphs[seg];

    // パスはいくつか候補を試す
    // 実際に使ってるのが ../RAYRAWTEMP/RAYRAWTEMP.<seg> だったのでまずそれ
    std::vector<std::string> candidates = {
      Form("param/RAYRAWTEMP/RAYRAWTEMP.%d", seg),
      Form("../RAYRAWTEMP/RAYRAWTEMP.%d", seg),
      Form("RAYRAWTEMP/RAYRAWTEMP.%d", seg)
    };

    std::ifstream fin;
    std::string   used_path;
    for (auto &p : candidates) {
      fin.open(p.c_str());
      if (fin.is_open()) {
        used_path = p;
        break;
      }
      fin.clear();
    }

    if (!fin.is_open()) {
      std::cerr << "[UserRayrawtemp] ERROR: template file for seg " << seg
                << " not found." << std::endl;
      return nullptr;
    }

    TGraph *gr = new TGraph();
    double t, v;
    int    n = 0;
    while (fin >> t >> v) {
      gr->SetPoint(n++, t, v);
    }
    fin.close();

    if ((int)gTemplateGraphs.size() <= seg)
      gTemplateGraphs.resize(seg+1, nullptr);
    gTemplateGraphs[seg] = gr;

    // std::cout << "[UserRayrawtemp] loaded template: " << used_path
    //           << "  N=" << n << std::endl;

    return gr;
  }

  // TGraph から単純線形補間で値を取る
  double EvalTemplate(const TGraph *gr, double t)
  {
    if (!gr) return 0.0;

    int n = gr->GetN();
    if (n == 0) return 0.0;

    const double *x = gr->GetX();
    const double *y = gr->GetY();

    if (t <= x[0])          return y[0];
    if (t >= x[n-1])        return y[n-1];

    // 2分探索でもいいけどn小さいので線形で
    for (int i = 0; i < n-1; ++i) {
      if (x[i] <= t && t < x[i+1]) {
        double u = (t - x[i]) / (x[i+1] - x[i]);
        return y[i] + u * (y[i+1] - y[i]);
      }
    }
    return 0.0;
  }

} // anonymous namespace

/*================================================================
 *  ProcessingBegin
 *===============================================================*/
Bool_t
ProcessingBegin()
{
  root::event.clear();
  return true;
}

/*================================================================
 *  ProcessingNormal : イベントごとに呼ばれる
 *===============================================================*/
Bool_t
ProcessingNormal()
{
  RawData      rawData;
  HodoAnalyzer hodoAna(rawData);

  // イベント番号
  root::event.evnum = gUnpacker.get_event_number();
  HF1(1, 0);  // status histo  (元のUserTemplateと合わせておく)

  // RAYRAW をデコードして RayrawWaveformHit として作ってもらう
  rawData.DecodeHits("RAYRAW");
  hodoAna.DecodeHits<RayrawWaveformHit>("RAYRAW");

  const Int_t nSeg = NumOfSegRayraw;
  const Int_t nHit = hodoAna.GetNHits("RAYRAW");

  // セグメントごとの入れ物を初期化
  root::event.fit_time.assign(nSeg, {});
  root::event.fit_height.assign(nSeg, {});
  root::event.waveform_time.assign(nSeg, {});
  root::event.waveform_adc.assign(nSeg, {});
  root::event.fitted_time.assign(nSeg, {});
  root::event.fitted_adc.assign(nSeg, {});

  // 再構成の刻み幅 = (元のサンプリング間隔)/kOverSample
  // UserParamMan から SamplingInterval をもらう (あなたの環境だと 0.5 ns)
  const double sampling = gUser.GetParameter("SamplingInterval", 0.5);
  const double dt_reco  = sampling / (double)kOverSample;

  // ============================================================
  //   ヒットを1本ずつ処理
  // ============================================================
  for (Int_t ihit = 0; ihit < nHit; ++ihit) {

    const auto *hit = hodoAna.GetHit<RayrawWaveformHit>("RAYRAW", ihit);
    if (!hit) continue;

    Int_t seg = hit->SegmentId();
    if (seg < 0 || seg >= nSeg) continue;

    // ------------------------------------------------------------
    // 生波形 (ADC) をコピー
    // ------------------------------------------------------------
    auto &wf_t = root::event.waveform_time[seg];
    auto &wf_a = root::event.waveform_adc[seg];

    const Int_t nWave = hit->GetWaveformEntries(HodoRawHit::kUp);
    wf_t.reserve(nWave);
    wf_a.reserve(nWave);

    for (Int_t i = 0; i < nWave; ++i) {
      auto p = hit->GetWaveform(HodoRawHit::kUp, i);
      wf_t.push_back(p.first);   // time [ns]
      wf_a.push_back(p.second);  // ADC (pedestal-subtracted)
    }

    // ------------------------------------------------------------
    // テンプレートフィットで見つかったパルスをコピー
    // ------------------------------------------------------------
    auto &ft = root::event.fit_time[seg];
    auto &fh = root::event.fit_height[seg];

    const Int_t nPulse = hit->GetNPulse(HodoRawHit::kUp);
    for (Int_t ip = 0; ip < nPulse; ++ip) {
      ft.push_back( hit->GetPulseTime  (HodoRawHit::kUp, ip) );
      fh.push_back( hit->GetPulseHeight(HodoRawHit::kUp, ip) );
    }

    // ------------------------------------------------------------
    // フィットがなかった場合のフォールバック
    //  (一番高いADCのところを1パルスとして入れておく)
    // ------------------------------------------------------------
    if (nPulse == 0 && !wf_a.empty()) {
      auto it = std::max_element(wf_a.begin(), wf_a.end());
      double maxadc = *it;
      if (maxadc > 1.0) { // 閾値は適当
        int idx = std::distance(wf_a.begin(), it);
        double t0 = wf_t[idx];
        ft.push_back(t0);
        fh.push_back(maxadc);
      }
    }

    // ------------------------------------------------------------
    // ここから“細かいテンプレート波形”を作る
    // ------------------------------------------------------------
    // まずこのseg用のテンプレートをロード
    TGraph *gr_temp = LoadTemplateGraphForSeg(seg);

    auto &fw_t = root::event.fitted_time[seg];
    auto &fw_a = root::event.fitted_adc[seg];

    // パルスが1つもなければ空のままでOK
    if (!ft.empty() && gr_temp) {

      const int nstep = (int)((kRecoTmax - kRecoTmin) / dt_reco) + 1;
      fw_t.reserve(nstep);
      fw_a.reserve(nstep);

      for (int istep = 0; istep < nstep; ++istep) {
        double t = kRecoTmin + istep * dt_reco;
        double val = 0.0;

        // パルスを全部足し合わせる
        for (size_t ip = 0; ip < ft.size(); ++ip) {
          double t0 = ft[ip];
          double A  = fh[ip];
          double templ = EvalTemplate(gr_temp, t - t0);  // t-t0 で評価
          val += A * templ;
        }

        fw_t.push_back(t);
        fw_a.push_back(val);
      }
    }

  } // for ihit

  return true;
}

/*================================================================
 *  ProcessingEnd : 1イベントぶん終わったらTTreeに書く
 *===============================================================*/
Bool_t
ProcessingEnd()
{
  root::tree->Fill();
  return true;
}

/*================================================================
 *  ConfMan::InitializeHistograms
 *===============================================================*/
Bool_t
ConfMan::InitializeHistograms()
{
  // ---- ヒストはとりあえず UserTemplate.cc とほぼ同じ ----
  const Int_t    NbinFADC_X     = 500;
  const Int_t    NbinFADC_Y     = 1024;
  const Double_t MinIntegral    = 0;
  const Double_t MaxIntegral    = 10000;
  const Int_t    NbinIntegral   = (Int_t)(MaxIntegral - MinIntegral);
  const Double_t MinTDC         = 0.;
  const Double_t MaxTDC         = 3000000.;
  const Int_t    NbinTDC        = (Int_t)(MaxTDC - MinTDC);
  const Double_t MindE          = -0.5;
  const Double_t MaxdE          = 29.5;
  const Int_t    NbindE         = (Int_t)(MaxdE - MindE) * 100;

  HB1(1, "Status", 20, 0., 20.);

  for (Int_t seg = 0; seg < NumOfSegRayraw; ++seg) {
    TString title0  = Form("RAYRAW seg%d - Raw Waveform (ADC)", seg);
    TString title1  = Form("RAYRAW seg%d - Max ADC",            seg);
    TString title2  = Form("RAYRAW seg%d - Integral",           seg);
    TString title3  = Form("RAYRAW seg%d - Pedestal Integral",  seg);
    TString title4  = Form("RAYRAW seg%d - TDC Leading",        seg);
    TString title5  = Form("RAYRAW seg%d - TDC Trailing",       seg);
    TString title6  = Form("RAYRAW seg%d - TDC First",          seg);
    TString title7  = Form("RAYRAW seg%d - TOT",                seg);
    TString title10 = Form("RAYRAW seg%d - dE(ADC)",            seg);
    TString title11 = Form("RAYRAW seg%d - dE(Integral)",       seg);

    HB2( RAYRAWHid + (seg+1)*1000 + 0,  title0,
         NbinFADC_X, 0., (Double_t)NbinFADC_X,
         NbinFADC_Y, 0., (Double_t)NbinFADC_Y );
    HB1( RAYRAWHid + (seg+1)*1000 + 1,  title1,
         NbinFADC_Y,   0.,          (Double_t)NbinFADC_Y );
    HB1( RAYRAWHid + (seg+1)*1000 + 2,  title2,
         NbinIntegral, MinIntegral, MaxIntegral );
    HB1( RAYRAWHid + (seg+1)*1000 + 3,  title3,
         NbinIntegral, MinIntegral, MaxIntegral );
    HB1( RAYRAWHid + (seg+1)*1000 + 4,  title4,
         NbinTDC,      MinTDC,      MaxTDC );
    HB1( RAYRAWHid + (seg+1)*1000 + 5,  title5,
         NbinTDC,      MinTDC,      MaxTDC );
    HB1( RAYRAWHid + (seg+1)*1000 + 6,  title6,
         NbinTDC,      MinTDC,      MaxTDC );
    HB1( RAYRAWHid + (seg+1)*1000 + 7,  title7,
         NbinTDC,      MinTDC,      MaxTDC );
    HB1( RAYRAWHid + (seg+1)*1000 + 10, title10,
         NbindE,       MindE,       MaxdE );
    HB1( RAYRAWHid + (seg+1)*1000 + 11, title11,
         NbindE,       MindE,       MaxdE );
  }

  // ---- TTree の準備 ----
  HBTree("tree", "tree");
  root::tree->Branch("evnum",          &root::event.evnum,        "evnum/I");
  root::tree->Branch("fit_time",       &root::event.fit_time);
  root::tree->Branch("fit_height",     &root::event.fit_height);
  root::tree->Branch("waveform_time",  &root::event.waveform_time);
  root::tree->Branch("waveform_adc",   &root::event.waveform_adc);
  // 新しく追加した「細かく作ったテンプレート波形」
  root::tree->Branch("fitted_time",    &root::event.fitted_time);
  root::tree->Branch("fitted_adc",     &root::event.fitted_adc);

  HPrint();
  return true;
}

/*================================================================
 *  ConfMan::InitializeParameterFiles
 *===============================================================*/
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    ( InitializeParameter<DCGeomMan>   ("DCGEO")     &&
      InitializeParameter<HodoParamMan>("HDPRM")     &&
      InitializeParameter<HodoPHCMan>  ("HDPHC")     &&
      InitializeParameter<UserParamMan>("USER")      &&
      InitializeParameter<TemplateFitMan>("RAYRAWTEMP")
    );
}

/*================================================================
 *  ConfMan::FinalizeProcess
 *===============================================================*/
Bool_t
ConfMan::FinalizeProcess()
{
  // 特にやることがなければ true でOK
  return true;
}

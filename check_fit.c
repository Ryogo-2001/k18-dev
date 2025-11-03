#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm> // std::min, std::max

// ===== 設定 =====
const char* ROOT_FILE_NAME   = "tempfit1358.root";
const int   SEGMENT_TO_CHECK = 0;
// -1 なら「そのsegでfit_timeが空じゃない最初のイベント」を自動で探す
Long64_t    EVENT_TO_CHECK   = -1;
// =================

// テンプレートを1つ読む
TGraph* LoadTemplateGraph(const char* filename)
{
  std::ifstream infile(filename);
  if (!infile) {
    std::cerr << "!!! エラー: テンプレートファイル " << filename << " が見つかりません。" << std::endl;
    return nullptr;
  }

  auto g = new TGraph();
  double t, v;
  int n = 0;
  while (infile >> t >> v) {
    g->SetPoint(n++, t, v);
  }
  return g;
}

void check_fit()
{
  // 1. テンプレート読み込み
  TString template_path = Form("param/RAYRAWTEMP/RAYRAWTEMP.%d", SEGMENT_TO_CHECK);
  TGraph* gr_template = LoadTemplateGraph(template_path.Data());
  if (!gr_template || gr_template->GetN() == 0) {
    std::cerr << "テンプレートが読めていません。path=" << template_path << std::endl;
    return;
  }

  // 2. ROOTファイルを開く
  TFile* f_in = TFile::Open(ROOT_FILE_NAME);
  if (!f_in || f_in->IsZombie()) {
    std::cerr << "エラー: " << ROOT_FILE_NAME << " が開けません。" << std::endl;
    return;
  }

  TTree* tree = (TTree*)f_in->Get("tree");
  if (!tree) {
    std::cerr << "エラー: TTree 'tree' が見つかりません。" << std::endl;
    f_in->Close();
    return;
  }

  // 3. ブランチセット
  std::vector<std::vector<double>>* waveform_time = nullptr;
  std::vector<std::vector<double>>* waveform_volt = nullptr;
  std::vector<std::vector<double>>* fit_time      = nullptr;
  std::vector<std::vector<double>>* fit_height    = nullptr;

  tree->SetBranchAddress("waveform_time", &waveform_time);
  tree->SetBranchAddress("waveform_volt", &waveform_volt);
  tree->SetBranchAddress("fit_time",      &fit_time);
  tree->SetBranchAddress("fit_height",    &fit_height);

  // 4. イベントを探す
  if (EVENT_TO_CHECK == -1) {
    std::cout << "セグメント " << SEGMENT_TO_CHECK
              << " でfit済みパルスを含む最初のイベントを探します..." << std::endl;
    Long64_t nent = tree->GetEntries();
    for (Long64_t iEvt = 0; iEvt < nent; ++iEvt) {
      tree->GetEntry(iEvt);
      if (fit_time && fit_time->size() > (size_t)SEGMENT_TO_CHECK) {
        if (!fit_time->at(SEGMENT_TO_CHECK).empty()) {
          EVENT_TO_CHECK = iEvt;
          std::cout << "→ イベント " << iEvt << " にパルスあり" << std::endl;
          break;
        }
      }
    }
    if (EVENT_TO_CHECK == -1) {
      std::cerr << "該当するイベントが見つかりませんでした。" << std::endl;
      f_in->Close();
      return;
    }
  }

  // 念のためもう一度読んでおく
  tree->GetEntry(EVENT_TO_CHECK);

  // 5. 波形チェック
  if (!waveform_time || !waveform_volt ||
      waveform_time->size() <= (size_t)SEGMENT_TO_CHECK ||
      waveform_volt->size()  <= (size_t)SEGMENT_TO_CHECK ||
      waveform_time->at(SEGMENT_TO_CHECK).empty() ||
      waveform_volt->at(SEGMENT_TO_CHECK).empty()) {

    std::cerr << "イベント " << EVENT_TO_CHECK
              << " のセグメント " << SEGMENT_TO_CHECK
              << " の生波形が空です。" << std::endl;
    f_in->Close();
    return;
  }

  const auto& wf_t = waveform_time->at(SEGMENT_TO_CHECK);
  const auto& wf_v = waveform_volt->at(SEGMENT_TO_CHECK);

  auto* gr_raw = new TGraph(wf_t.size(), wf_t.data(), wf_v.data());
  gr_raw->SetLineColor(kGray+2);
  gr_raw->SetLineWidth(1);
  gr_raw->SetTitle(Form("Fit Check (Event %lld, Seg %d);Time [ns];ADC", EVENT_TO_CHECK, SEGMENT_TO_CHECK));

  // 6. フィット結果を重ねる
  const auto& pulses_t = fit_time->at(SEGMENT_TO_CHECK);
  const auto& pulses_h = fit_height->at(SEGMENT_TO_CHECK);

  std::cout << pulses_t.size() << " 個のパルスが見つかりました。" << std::endl;
  for (size_t i = 0; i < pulses_t.size(); ++i) {
    std::cout << "  Pulse " << i
              << ": t0 = " << pulses_t[i]
              << " ns, A = " << pulses_h[i] << std::endl;
  }

  // 7. 描画
  gStyle->SetOptStat(0);
  auto* c1 = new TCanvas("c1", "Fit Check Canvas", 900, 600);
  gr_raw->Draw("AL");

  // Yレンジ用のmin/maxを生波形から取る
  double ymin = +1e9;
  double ymax = -1e9;
  for (size_t i = 0; i < wf_v.size(); ++i) {
    ymin = std::min(ymin, wf_v[i]);
    ymax = std::max(ymax, wf_v[i]);
  }

  // パルスを重ねる
  for (size_t i = 0; i < pulses_t.size(); ++i) {
    double t0 = pulses_t[i];
    double A  = pulses_h[i];

    auto* gr_pulse = (TGraph*)gr_template->Clone(Form("pulse_%zu", i));
    for (int j = 0; j < gr_pulse->GetN(); ++j) {
      double xt = gr_pulse->GetX()[j];
      double yt = gr_pulse->GetY()[j];
      gr_pulse->SetPoint(j, xt + t0, yt * A);
      ymin = std::min(ymin, yt * A);
      ymax = std::max(ymax, yt * A);
    }
    gr_pulse->SetLineColor(kRed + (int)i % 9);
    gr_pulse->SetLineWidth(2);
    gr_pulse->Draw("L SAME");
  }

  // ちょっとだけマージン
  gr_raw->GetYaxis()->SetRangeUser(ymin - 0.1*(std::abs(ymax)+1),
                                   ymax + 0.1*(std::abs(ymax)+1));

  c1->SaveAs("fit_check.pdf");
  std::cout << "→ fit_check.pdf を保存しました" << std::endl;

  f_in->Close();
}

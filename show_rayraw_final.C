// show_rayraw_final.C
// tempfitXXXXX.root に入っている
//   waveform_time (raw)
//   waveform_adc  (raw)
//   fitted_time   (resampled)
//   fitted_adc    (resampled)
//   fit_time      (pulse time)
//   fit_height    (pulse height)
// を1枚にまとめて表示・保存する

#include <vector>
#include <iostream>

void show_rayraw_final(int evt = 0, int seg = 0,
                       const char *outname = "")
{
  // tree を取る
  TTree *t = (TTree*)gDirectory->Get("tree");
  if (!t) {
    std::cerr << "[show_rayraw_final] TTree \"tree\" が見つかりません。先に root -l tempfit1358.root してください。\n";
    return;
  }

  // ブランチ用ポインタ
  std::vector<std::vector<double>> *waveform_time = nullptr;
  std::vector<std::vector<double>> *waveform_adc  = nullptr;
  std::vector<std::vector<double>> *fit_time      = nullptr;
  std::vector<std::vector<double>> *fit_height    = nullptr;
  std::vector<std::vector<double>> *fitted_time   = nullptr;
  std::vector<std::vector<double>> *fitted_adc    = nullptr;
  Int_t evnum = -1;

  t->SetBranchAddress("evnum",         &evnum);
  t->SetBranchAddress("waveform_time", &waveform_time);
  t->SetBranchAddress("waveform_adc",  &waveform_adc);
  t->SetBranchAddress("fit_time",      &fit_time);
  t->SetBranchAddress("fit_height",    &fit_height);
  t->SetBranchAddress("fitted_time",   &fitted_time);
  t->SetBranchAddress("fitted_adc",    &fitted_adc);

  Long64_t nent = t->GetEntries();
  if (evt < 0 || evt >= nent) {
    std::cerr << "[show_rayraw_final] evt=" << evt << " は範囲外です (0.." << (nent-1) << ")\n";
    return;
  }

  // イベントを読む
  t->GetEntry(evt);

  // seg の存在チェック
  auto badseg = [&](const std::vector<std::vector<double>> *v){
    return (!v || seg < 0 || seg >= (int)v->size());
  };

  if (badseg(waveform_time) || badseg(waveform_adc)) {
    std::cerr << "[show_rayraw_final] このイベントには seg=" << seg << " の生波形がありません\n";
    return;
  }

  // それぞれ取り出す
  auto &x_raw = waveform_time->at(seg);
  auto &y_raw = waveform_adc->at(seg);

  // 描画用キャンバス
  TCanvas *c = new TCanvas(Form("c_evt%d_seg%d", evt, seg),
                           Form("evt %d (file evnum=%d), seg %d", evt, evnum, seg),
                           900, 600);
  c->SetGrid();

  // まず生波形
  TGraph *g_raw = nullptr;
  double ymin =  1e9;
  double ymax = -1e9;
  if (!x_raw.empty()) {
    g_raw = new TGraph(x_raw.size(), x_raw.data(), y_raw.data());
    g_raw->SetLineColor(kBlack);
    g_raw->SetLineWidth(1);
    g_raw->SetTitle(Form("evt %d (tree entry) / evnum %d / seg %d;time [ns];ADC", evt, evnum, seg));
    g_raw->Draw("AL");

    for (size_t i=0; i<y_raw.size(); ++i) {
      if (y_raw[i] < ymin) ymin = y_raw[i];
      if (y_raw[i] > ymax) ymax = y_raw[i];
    }
  }

  // 次にフィットで再サンプルした波形（細かい赤）
  TGraph *g_fitw = nullptr;
  bool has_fitw = false;
  if (!badseg(fitted_time) && !badseg(fitted_adc)
      && !fitted_time->at(seg).empty()) {

    auto &x_fitw = fitted_time->at(seg);
    auto &y_fitw = fitted_adc->at(seg);

    g_fitw = new TGraph(x_fitw.size(), x_fitw.data(), y_fitw.data());
    g_fitw->SetLineColor(kRed);
    g_fitw->SetLineWidth(2);

    if (g_raw) g_fitw->Draw("L SAME");
    else       g_fitw->Draw("AL");

    has_fitw = true;

    for (size_t i=0; i<y_fitw.size(); ++i) {
      if (y_fitw[i] < ymin) ymin = y_fitw[i];
      if (y_fitw[i] > ymax) ymax = y_fitw[i];
    }
  }

  // パルス位置 (青丸)
  bool has_pulse = false;
  if (!badseg(fit_time) && !badseg(fit_height)) {
    auto &pt = fit_time->at(seg);
    auto &ph = fit_height->at(seg);
    if (!pt.empty() && pt.size() == ph.size()) {
      has_pulse = true;
      for (size_t i=0; i<pt.size(); ++i) {
        TMarker *m = new TMarker(pt[i], ph[i], 20);
        m->SetMarkerColor(kBlue+1);
        m->SetMarkerSize(1.1);
        m->Draw();

        // y のレンジ更新
        if (ph[i] < ymin) ymin = ph[i];
        if (ph[i] > ymax) ymax = ph[i];

        // ラベルもつけておく（pulse0, pulse1,…）
        TLatex *lat = new TLatex(pt[i], ph[i] + 0.03*(ymax - ymin), Form("p%d", (int)i));
        lat->SetTextSize(0.03);
        lat->SetTextColor(kBlue+2);
        lat->Draw();
      }
    }
  }

  // y 範囲をちょっと広げる
  if (ymax > ymin) {
    double dy = ymax - ymin;
    gPad->Update();
    gPad->GetUymin();
    gPad->GetUymax();
    gPad->Modified();
    gPad->Update();
    // axis を書き換えるには
    gPad->GetFrame()->SetY1(ymin - 0.05*dy);
    gPad->GetFrame()->SetY2(ymax + 0.10*dy);
  }

  // 凡例
  TLegend *leg = new TLegend(0.65, 0.70, 0.88, 0.90);
  if (g_raw)  leg->AddEntry(g_raw,  "raw ADC", "l");
  if (has_fitw) leg->AddEntry(g_fitw, "fitted waveform", "l");
  if (has_pulse) leg->AddEntry((TObject*)0, "blue = pulse", "");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  // 保存: outname が指定されてたらそれを使う
  if (outname && std::string(outname).size()) {
    c->Print(outname);
    std::cout << "[show_rayraw_final] saved to " << outname << "\n";
  } else {
    // 自動でいい感じの名前にする
    TString fname = Form("rayraw_evt%05d_seg%02d.pdf", evt, seg);
    c->Print(fname);
    std::cout << "[show_rayraw_final] saved to " << fname << "\n";
  }
}

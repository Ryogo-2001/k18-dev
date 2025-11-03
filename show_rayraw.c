void show_rayraw(int evt = 0, int seg = 0)
{
  TTree *t = (TTree*)gDirectory->Get("tree");
  if (!t) { std::cerr << "no tree\n"; return; }

  std::vector<std::vector<double>> *wf_t = nullptr, *wf_a = nullptr;
  std::vector<std::vector<double>> *ft_t = nullptr, *ft_a = nullptr;
  std::vector<std::vector<double>> *rt_t = nullptr, *rt_a = nullptr;

  t->SetBranchAddress("waveform_time", &wf_t);
  t->SetBranchAddress("waveform_adc",  &wf_a);
  t->SetBranchAddress("fit_time",      &ft_t);
  t->SetBranchAddress("fit_height",    &ft_a);
  t->SetBranchAddress("fitted_time",   &rt_t);
  t->SetBranchAddress("fitted_adc",    &rt_a);

  Long64_t nent = t->GetEntries();
  if (evt < 0 || evt >= nent) {
    std::cerr << "evt out of range\n";
    return;
  }

  t->GetEntry(evt);

  if (seg < 0 || seg >= (int)wf_t->size()) {
    std::cerr << "seg out of range\n";
    return;
  }

  auto &x_raw = wf_t->at(seg);
  auto &y_raw = wf_a->at(seg);

  TCanvas *c = new TCanvas("c","rayraw view",900,600);
  c->SetGrid();

  // 生波形（ADCそのまま）
  TGraph *g_raw = nullptr;
  if (!x_raw.empty()) {
    g_raw = new TGraph(x_raw.size(), x_raw.data(), y_raw.data());
    g_raw->SetLineColor(kBlack);
    g_raw->SetLineWidth(1);
    g_raw->SetTitle(Form("evt %d, seg %d;time [ns];ADC", evt, seg));
    g_raw->Draw("AL");
  } else {
    std::cout << "no raw waveform for this seg\n";
  }

  // テンプレで再サンプリングしたやつ
  if (rt_t && rt_a && seg < (int)rt_t->size() && !rt_t->at(seg).empty()) {
    auto &x_fit = rt_t->at(seg);
    auto &y_fit = rt_a->at(seg);
    TGraph *g_fit = new TGraph(x_fit.size(), x_fit.data(), y_fit.data());
    g_fit->SetLineColor(kRed);
    g_fit->SetLineWidth(2);
    if (g_raw) g_fit->Draw("L SAME");
    else       g_fit->Draw("AL");
  } else {
    std::cout << "no fitted waveform for this seg\n";
  }

  // パルスの位置だけマーカーで見たいとき
  if (ft_t && ft_a && seg < (int)ft_t->size()) {
    auto &pt = ft_t->at(seg);
    auto &ph = ft_a->at(seg);
    for (size_t i=0; i<pt.size(); ++i) {
      TMarker *m = new TMarker(pt[i], ph[i], 20);
      m->SetMarkerColor(kBlue+1);
      m->SetMarkerSize(1.0);
      m->Draw();
    }
  }
}

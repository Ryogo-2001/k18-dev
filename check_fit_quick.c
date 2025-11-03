void check_fit_quick()
{
  const char* ROOT_FILE_NAME   = "tempfit1358.root";
  const int   SEGMENT_TO_CHECK = 0;

  TFile* f = TFile::Open(ROOT_FILE_NAME);
  if (!f || f->IsZombie()) {
    std::cerr << "ファイルが開けません: " << ROOT_FILE_NAME << std::endl;
    return;
  }
  TTree* tree = (TTree*)f->Get("tree");
  if (!tree) {
    std::cerr << "tree がありません" << std::endl;
    return;
  }

  std::vector<std::vector<double>>* waveform_time = nullptr;
  std::vector<std::vector<double>>* waveform_volt = nullptr;
  std::vector<std::vector<double>>* fit_time      = nullptr;

  tree->SetBranchAddress("waveform_time", &waveform_time);
  tree->SetBranchAddress("waveform_volt", &waveform_volt);
  tree->SetBranchAddress("fit_time",      &fit_time);

  Long64_t nent = tree->GetEntries();
  Long64_t n_with_wave = 0;
  Long64_t n_with_fit  = 0;

  for (Long64_t i=0; i<nent; ++i) {
    tree->GetEntry(i);

    bool has_wave = waveform_volt &&
                    waveform_volt->size() > (size_t)SEGMENT_TO_CHECK &&
                    !waveform_volt->at(SEGMENT_TO_CHECK).empty();

    bool has_fit  = fit_time &&
                    fit_time->size() > (size_t)SEGMENT_TO_CHECK &&
                    !fit_time->at(SEGMENT_TO_CHECK).empty();

    if (has_wave) n_with_wave++;
    if (has_fit)  n_with_fit++;
  }

  std::cout << "総イベント数      : " << nent << std::endl;
  std::cout << "波形があるイベント: " << n_with_wave << std::endl;
  std::cout << "fitがあるイベント : " << n_with_fit  << std::endl;

  // とりあえず最初の波形だけ描く
  int first_evt = -1;
  for (Long64_t i=0; i<nent; ++i) {
    tree->GetEntry(i);
    if (waveform_time->size() > (size_t)SEGMENT_TO_CHECK &&
        !waveform_time->at(SEGMENT_TO_CHECK).empty()) {
      first_evt = i;
      break;
    }
  }

  if (first_evt == -1) {
    std::cerr << "描画できる波形がありませんでした" << std::endl;
    return;
  }

  tree->GetEntry(first_evt);
  const auto& wf_t = waveform_time->at(SEGMENT_TO_CHECK);
  const auto& wf_v = waveform_volt->at(SEGMENT_TO_CHECK);

  auto* gr = new TGraph(wf_t.size(), wf_t.data(), wf_v.data());
  gr->SetTitle(Form("Event %d, Seg %d;time [ns];ADC", first_evt, SEGMENT_TO_CHECK));
  gr->SetLineColor(kBlue+1);

  gStyle->SetOptStat(0);
  auto* c1 = new TCanvas("c1","waveform",900,600);
  gr->Draw("AL");
  c1->SaveAs("waveform_first.pdf");

  std::cout << "waveform_first.pdf を保存しました" << std::endl;

  f->Close();
}

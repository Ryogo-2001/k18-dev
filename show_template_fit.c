// show_template_fit.c
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TLine.h>
#include <iostream>
#include <fstream>
#include <vector>


const char* ROOT_FILE_NAME   = "tempfit1358.root";
const int   SEG              = 0; // 見たいセグメント
const char* TEMPLATE_BASENAME = "param/RAYRAWTEMP/RAYRAWTEMP"; 



bool LoadTemplate(int seg, std::vector<double>& xt, std::vector<double>& yt)
{
  TString fname = Form("%s.%d", TEMPLATE_BASENAME, seg);
  std::ifstream fin(fname.Data());
  if (!fin) {
    std::cerr << "テンプレートファイルが開けません: " << fname << std::endl;
    return false;
  }
  double x,y;
  while (fin >> x >> y) {
    xt.push_back(x);
    yt.push_back(y);
  }
  return !xt.empty();
}

void show_template_fit()
{
  
  TFile* f = TFile::Open(ROOT_FILE_NAME);
  if (!f || f->IsZombie()) {
    std::cerr << "cannot open " << ROOT_FILE_NAME << std::endl;
    return;
  }
  TTree* tree = (TTree*)f->Get("tree");
  if (!tree) {
    std::cerr << "tree がありません" << std::endl;
    f->Close();
    return;
  }

  
  std::vector<std::vector<double>> *waveform_time = nullptr;
  std::vector<std::vector<double>> *waveform_volt = nullptr;
  std::vector<std::vector<double>> *fit_time      = nullptr;
  std::vector<std::vector<double>> *fit_height    = nullptr;

  tree->SetBranchAddress("waveform_time", &waveform_time);
  tree->SetBranchAddress("waveform_volt", &waveform_volt);
  tree->SetBranchAddress("fit_time",      &fit_time);
  tree->SetBranchAddress("fit_height",    &fit_height);

  Long64_t nent = tree->GetEntries();

 
  Long64_t draw_evt = -1;
  Long64_t n_fit_evt = 0;
  for (Long64_t i = 0; i < nent; ++i) {
    tree->GetEntry(i);
    bool has_fit = (fit_time && fit_time->size() > (size_t)SEG &&
                    !fit_time->at(SEG).empty());
    if (has_fit) {
      if (draw_evt < 0) draw_evt = i;
      n_fit_evt++;
    }
  }

  std::cout << "総イベント数: " << nent << std::endl;
  std::cout << "SEG " << SEG << " でfitがあるイベント: " << n_fit_evt << std::endl;

  if (draw_evt < 0) {
    std::cerr << "このセグメントでfit済みイベントが見つかりませんでした。" << std::endl;
    f->Close();
    return;
  }

  
  tree->GetEntry(draw_evt);

  
  if (!waveform_time || waveform_time->size() <= (size_t)SEG ||
      waveform_time->at(SEG).empty()) {
    std::cerr << "選んだイベントに波形がありません。" << std::endl;
    f->Close();
    return;
  }

  const auto& wf_t = waveform_time->at(SEG);
  const auto& wf_v = waveform_volt->at(SEG);

  
  std::vector<double> tmpl_t, tmpl_y;
  if (!LoadTemplate(SEG, tmpl_t, tmpl_y)) {
    f->Close();
    return;
  }

  
  std::vector<double> fit_sum_y(wf_t.size(), 0.0);

  
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1", "Template fit view", 900, 600);

  
  TGraph *gr_raw = new TGraph(wf_t.size(), wf_t.data(), wf_v.data());
  gr_raw->SetLineColor(kBlack);
  gr_raw->SetLineWidth(2);
  gr_raw->SetTitle(Form("Event %lld, seg %d;time [ns];ADC", draw_evt, SEG));
  gr_raw->Draw("AL");

  
  double ymin = gr_raw->GetYaxis()->GetXmin();
  double ymax = gr_raw->GetYaxis()->GetXmax();

  
  const auto &ft = fit_time->at(SEG);
  const auto &fh = fit_height->at(SEG);
  std::cout << "このイベントのテンプレートパルス数: " << ft.size() << std::endl;

  for (size_t ip = 0; ip < ft.size(); ++ip) {
    double t0 = ft[ip];   // ns
    double A  = fh[ip];   // scale

    
    TGraph *gpulse = new TGraph(tmpl_t.size());
    for (size_t j = 0; j < tmpl_t.size(); ++j) {
      double x = tmpl_t[j] + t0;
      double y = tmpl_y[j] * A;
      gpulse->SetPoint(j, x, y);

      
      double dt = wf_t[1] - wf_t[0]; 
      int idx = (int)std::round((x - wf_t[0]) / dt);
      if (0 <= idx && idx < (int)fit_sum_y.size()) {
        fit_sum_y[idx] += y;
      }
    }
    gpulse->SetLineColor(kRed + (ip % 5));
    gpulse->SetLineWidth(2);
    gpulse->Draw("L SAME");

    
    for (int j = 0; j < gpulse->GetN(); ++j) {
      double xx, yy;
      gpulse->GetPoint(j, xx, yy);
      if (yy < ymin) ymin = yy;
      if (yy > ymax) ymax = yy;
    }

    std::cout << "  pulse " << ip
              << " : t0 = " << t0
              << " ns, A = " << A << std::endl;
  }

 
  TGraph *gr_sum = new TGraph(wf_t.size());
  for (size_t i = 0; i < wf_t.size(); ++i) {
    gr_sum->SetPoint(i, wf_t[i], fit_sum_y[i]);
    if (fit_sum_y[i] < ymin) ymin = fit_sum_y[i];
    if (fit_sum_y[i] > ymax) ymax = fit_sum_y[i];
  }
  gr_sum->SetLineColor(kGreen+2);
  gr_sum->SetLineWidth(2);
  gr_sum->Draw("L SAME");

  
  gr_raw->GetYaxis()->SetRangeUser(ymin - 0.1*(std::abs(ymax-ymin)+1),
                                   ymax + 0.1*(std::abs(ymax-ymin)+1));

  
  TLegend *leg = new TLegend(0.12, 0.70, 0.40, 0.88);
  leg->AddEntry(gr_raw, "raw waveform", "l");
  leg->AddEntry(gr_sum, "sum of fitted pulses", "l");
  leg->AddEntry((TObject*)0, Form("event %lld", draw_evt), "");
  leg->AddEntry((TObject*)0, Form("seg %d", SEG), "");
  leg->Draw();

  c1->SaveAs("template_fit_view.pdf");
  std::cout << "template_fit_view.pdf を保存しました" << std::endl;

  f->Close();
}

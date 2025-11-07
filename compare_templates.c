#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <vector>

TGraph* loadGraph(const char* filename)
{
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return nullptr;
    }

    TGraph *g = new TGraph();
    double x, y;
    int n = 0;
    
    while (infile >> x >> y) {
        g->SetPoint(n, x, y);
        n++;
    }
    
    infile.close();
    
    if (n == 0) {
        std::cerr << "Warning: No data loaded from " << filename << std::endl;
        delete g;
        return nullptr;
    }
    
    std::cout << "Loaded " << n << " points from " << filename << std::endl;
    return g;
}

void compare_templates()
{
    const char* file_old = "param/RAYRAWTEMP/RAYRAWTEMP.0";
    const char* file_new = "param/NEW_RAYRAWTEMP/NEW_RAYRAWTEMP.0";

    TGraph *g_old = loadGraph(file_old);
    TGraph *g_new = loadGraph(file_new);

    if (!g_old || !g_new) {
        std::cerr << "Error loading graphs. Aborting." << std::endl;
        return;
    }

    TCanvas *c1 = new TCanvas("c_compare", "Template Comparison (Seg 0)", 800, 600);
    c1->SetGrid();

    // ★★★ 青の線と点を太く修正 ★★★
    g_old->SetLineColor(kBlue);
    g_old->SetLineWidth(3);     // 2 -> 3
    g_old->SetMarkerColor(kBlue);
    g_old->SetMarkerStyle(20);
    g_old->SetMarkerSize(1.0);  // 0.5 -> 1.0
    g_old->SetName("g_old");
    g_old->SetTitle("Original Template (RAYRAWTEMP.0)");

    g_new->SetLineColor(kRed);
    g_new->SetLineWidth(2);
    g_new->SetMarkerColor(kRed);
    g_new->SetMarkerStyle(24);
    g_new->SetMarkerSize(0.5);
    g_new->SetName("g_new");
    g_new->SetTitle("New Template (NEW_RAYRAWTEMP.0)");

    g_old->Draw("ALP");

    g_old->GetXaxis()->SetTitle("Time (relative to T_fit) [ns]");
    g_old->GetYaxis()->SetTitle("Normalized Amplitude");

   
    g_old->GetXaxis()->SetRangeUser(-50.0, 60.0); 
    
    // 縦軸（Y軸）の描画レンジを指定
    g_old->GetYaxis()->SetRangeUser(-0.2, 1.1);
    // ---------------------------------

    g_new->Draw("LP SAME");

    TLegend *leg = new TLegend(0.55, 0.75, 0.9, 0.9);
    leg->AddEntry(g_old, "Original (RAYRAWTEMP.0)", "lp");
    leg->AddEntry(g_new, "New (NEW_RAYRAWTEMP.0)", "lp");
    leg->Draw();

    c1->Update();
    
    c1->SaveAs("template_comparison.pdf");
    std::cout << "Saved canvas to template_comparison.pdf" << std::endl;
}
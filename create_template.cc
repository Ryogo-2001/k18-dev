#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TString.h>
#include <TF1.h>
#include <TROOT.h>

#include <iostream>
#include <fstream>
#include <string>

const int NumOfSegRayraw = 32;
const int RAYRAWHid = 100000;

void create_template(const char* in_filename = "run1358WaveformTemp.root")
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);

    TFile *file = new TFile(in_filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open input file: " << in_filename << std::endl;
        return;
    }

    gSystem->mkdir("param/NEW_RAYRAWTEMP", kTRUE);
    std::cout << "Output directory 'param/NEW_RAYRAWTEMP/' is ready." << std::endl;

    std::cout << "Generating PDF (rayraw_histograms.pdf)..." << std::endl;
    TCanvas *c_pdf = new TCanvas("c_pdf", "RAYRAW HWF Histograms", 1200, 1200);
    TString pdf_filename = "rayraw_histograms.pdf";
    c_pdf->Print(pdf_filename + "[");

    for (int seg = 0; seg < NumOfSegRayraw; ++seg) {
        c_pdf->Clear();
        c_pdf->Divide(2, 2);

        const int hid_201 = RAYRAWHid + (seg + 1) * 1000 + 201;
        const int hid_202 = RAYRAWHid + (seg + 1) * 1000 + 202;
        const int hid_203 = RAYRAWHid + (seg + 1) * 1000 + 203;
        const int hid_205 = RAYRAWHid + (seg + 1) * 1000 + 205;

        TH2* h201 = (TH2*)file->Get(Form("h%d", hid_201));
        TH1* h202 = (TH1*)file->Get(Form("h%d", hid_202));
        TH1* h203 = (TH1*)file->Get(Form("h%d", hid_203));
        TH2* h205 = (TH2*)file->Get(Form("h%d", hid_205));

        c_pdf->cd(1);
        if (h201) {
            h201->Draw("COLZ");
            h201->SetTitle(Form("HWF Waveform (Seg %d)", seg));
        }

        c_pdf->cd(2);
        if (h202) {
            h202->Draw();
            h202->SetTitle(Form("Pulse Height (Seg %d)", seg));
        }

        c_pdf->cd(3);
        if (h203) {
            h203->Draw();
            h203->SetTitle(Form("Pulse Time (Seg %d)", seg));
        }

        c_pdf->cd(4);
        if (h205) {
            h205->Draw("COLZ");
            h205->SetTitle(Form("Fit Fail Waveform (Seg %d)", seg));
        }

        c_pdf->Print(pdf_filename);
    }
    c_pdf->Print(pdf_filename + "]");
    delete c_pdf;
    std::cout << "PDF generation complete." << std::endl;

    std::cout << "Generating new template files in param/NEW_RAYRAWTEMP/ ..." << std::endl;
    TCanvas *c_template = new TCanvas("c_template", "Template Creation", 1200, 600);
    gStyle->SetOptStat(0);

    for (int seg = 0; seg < NumOfSegRayraw; ++seg) {
        const int hid_206 = RAYRAWHid + (seg + 1) * 1000 + 206;
        TH2* h206 = (TH2*)file->Get(Form("h%d", hid_206));

        if (!h206 || h206->GetEntries() == 0) {
            std::cout << "Skipping Seg " << seg << " (Hist 206 not found or empty)" << std::endl;
            continue;
        }

        h206->FitSlicesY(0, 0, -1, 1, "QNRG");

        TH1D *h_template = (TH1D*)gROOT->FindObject(Form("h%d_1", hid_206));

        if (!h_template) {
            std::cerr << "Error: FitSlicesY failed to create histogram for Seg " << seg << std::endl;
            continue;
        }
        
        h_template->SetTitle(Form("New Template by FitSlicesY (Seg %d)", seg));
        h_template->GetXaxis()->SetTitle("Time (relative to T_fit) [ns]");
        h_template->GetYaxis()->SetTitle("Gaussian Mean of Norm. Amplitude");

        c_template->Clear();
        c_template->Divide(2, 1);
        c_template->cd(1);
        h206->Draw("COLZ");
        c_template->cd(2);
        h_template->Draw();
        
        std::ofstream outfile(Form("param/NEW_RAYRAWTEMP/NEW_RAYRAWTEMP.%d", seg));
        if (!outfile.is_open()) {
            std::cerr << "Error: Cannot write to param/NEW_RAYRAWTEMP/NEW_RAYRAWTEMP." << seg << std::endl;
            continue;
        }

        for (int bin = 1; bin <= h_template->GetNbinsX(); ++bin) {
            double x_time = h_template->GetBinCenter(bin);
            double y_adc = h_template->GetBinContent(bin);

            if (y_adc != 0) {
                outfile << x_time << " " << y_adc << std::endl;
            }
        }
        outfile.close();
        
        delete h_template;

    }

    delete c_template;
    file->Close();
    delete file;
    std::cout << "All tasks complete." << std::endl;
}
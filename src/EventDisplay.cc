// -*- C++ -*-

#include "EventDisplay.hh"

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include <TApplication.h>
#include <TBRIK.h>
#include <TArc.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TEnv.h>
#include <TFile.h>
#include <TGeometry.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH2Poly.h>
#include <TF1.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TMarker3DBox.h>
#include <TMixture.h>
#include <TNode.h>
#include <TPad.h>
#include <TSystem.h>
#include <TPave.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TPolyLine.h>
#include <TPolyLine3D.h>
#include <TPolyMarker.h>
#include <TPolyMarker3D.h>
#include <TROOT.h>
#include <TRint.h>
#include <TRotMatrix.h>
#include <TStyle.h>
#include <TTRD1.h>
#include <TTRD2.h>
#include <TTUBS.h>
#include <TTUBS.h>
#include <TView.h>
#include <TLine.h>

#include <std_ostream.hh>

#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "Exception.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "DeleteUtility.hh"
#include "HodoParamMan.hh"
#include "DCTdcCalibMan.hh"

#define BH2        0
#define BcOut      1
#define KURAMA     1
#define SdcIn      1
#define SdcOut     1
#define SCH        1
#define TOF        1
#define WC         1
#define Vertex     0
#define Hist       0
#define Hist_Timing 0
#define Hist_SdcOut 0
#define Hist_BcIn   0

namespace
{
const DCGeomMan& gGeom = DCGeomMan::GetInstance();
const Int_t& IdBH1 = gGeom.DetectorId("BH1");
const Int_t& IdBH2 = gGeom.DetectorId("BH2");
const Int_t& IdSCH = gGeom.DetectorId("SCH");
const Int_t& IdTOF = gGeom.DetectorId("TOF");
const Int_t& IdWC = gGeom.DetectorId("WC");
const Int_t& IdTarget = gGeom.DetectorId("Target");
const Double_t& zTarget = gGeom.LocalZ("Target");
const Double_t& zHS = gGeom.LocalZ("HS");
const Double_t& zK18Target = gGeom.LocalZ("K18Target");
const Double_t& gzK18Target = gGeom.GlobalZ("K18Target");
// const Double_t& gxK18Target = gGeom.GetGlobalPosition("K18Target").x();
const Double_t& gxK18Target = -240.;
const Double_t& zBFT = gGeom.LocalZ("BFT");

//const Double_t BeamAxis = -150.; //E07
// const Double_t BeamAxis = -240.; //E40
const Double_t BeamAxis = -50.; //E42
#if Vertex
const Double_t MinX = -50.;
const Double_t MaxX =  50.;
const Double_t MinY = -50.;
const Double_t MaxY =  50.;
const Double_t MinZ = -25.;
#endif
const Double_t MaxZ =  50.;

const HodoParamMan& gHodo = HodoParamMan::GetInstance();
const DCTdcCalibMan& gTdc = DCTdcCalibMan::GetInstance();
}

//_____________________________________________________________________________
EventDisplay::EventDisplay()
  : m_is_ready(false),
    m_is_save_mode(),
    m_theApp(),
    m_geometry(),
    m_node(),
    m_canvas(),
    m_canvas_vertex(),
    m_canvas_hist(),
    m_canvas_hist2(),
    m_canvas_hist3(),
    m_canvas_hist4(),
    m_canvas_hist5(),
    m_canvas_hist6(),
    m_canvas_hist7(),
    m_canvas_hist8(),
    m_canvas_hist9(),
    m_hist_vertex_x(),
    m_hist_vertex_y(),
    m_hist_p(),
    m_hist_m2(),
    m_hist_missmass(),
    m_hist_bh1(),
    m_hist_bft(),
    m_hist_bft_p(),
    m_hist_bcIn(),
    m_hist_bcOut(),
    m_BH1box_cont(),
    m_BH2box_cont(),
    m_hist_bh2(),
    m_hist_bcOut_sdcIn(),
    m_hist_sdcIn_predict(),
    m_hist_sdcIn_predict2(),
    m_TargetXZ_box2(),
    m_TargetYZ_box2(),
    m_hist_sch(),
    m_hist_tof(),
    m_hist_sdc1(),
    m_hist_sdc1p(),
    m_hist_sdc3_l(),
    m_hist_sdc3_t(),
    m_hist_sdc3p_l(),
    m_hist_sdc3p_t(),
    m_hist_sdc3y_l(),
    m_hist_sdc3y_t(),
    m_hist_sdc3yp_l(),
    m_hist_sdc3yp_t(),
    m_hist_sdc4_l(),
    m_hist_sdc4_t(),
    m_hist_sdc4p_l(),
    m_hist_sdc4p_t(),
    m_hist_sdc4y_l(),
    m_hist_sdc4y_t(),
    m_hist_sdc4yp_l(),
    m_hist_sdc4yp_t(),
    m_hist_bc3(),
    m_hist_bc3p(),
    m_hist_bc3u(),
    m_hist_bc3up(),
    m_hist_bc3v(),
    m_hist_bc3vp(),
    m_hist_bc4(),
    m_hist_bc4p(),
    m_hist_bc4u(),
    m_hist_bc4up(),
    m_hist_bc4v(),
    m_hist_bc4vp(),
    m_hist_bc3_time(),
    m_hist_bc3p_time(),
    m_hist_bc4_time(),
    m_hist_bc4p_time(),
    m_target_node(),
    m_kurama_inner_node(),
    m_kurama_outer_node(),
    m_BH2wall_node(),
    m_SCHwall_node(),
    m_TOFwall_node(),
    m_WCwall_node(),
    m_BcOutTrack(),
    m_BcOutTrackShs(),
    m_SdcInTrack(),
    m_init_step_mark(),
    m_hs_step_mark(),
    m_kurama_step_mark(),
    m_TargetXZ_box(),
    m_TargetYZ_box(),
    m_VertexPointXZ(),
    m_VertexPointYZ(),
    m_HSMarkVertexXShs(),
    m_KuramaMarkVertexXShs(),
    m_KuramaMarkVertexX(),
    m_KuramaMarkVertexY(),
    m_MissMomXZ_line(),
    m_MissMomYZ_line()
{
}

//_____________________________________________________________________________
EventDisplay::~EventDisplay()
{
}

//_____________________________________________________________________________
// local function
void
ConstructionDone(const TString& name, std::ostream& ost=hddaq::cout)
{
  const Int_t n = 20;
  const Int_t s = name.Length();
  ost << " " << name << " ";
  for(Int_t i=0; i<n-s; ++i) ost << ".";
  ost << " done" << std::endl;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::Initialize()
{
  if(m_is_ready){
    hddaq::cerr << "#W " << FUNC_NAME
                << " already initialied" << std::endl;
    return false;
  }

  gStyle->SetOptStat(0);
  gStyle->SetStatH(0.040);
  gStyle->SetStatX(0.900);
  gStyle->SetStatY(0.900);
  Int_t myfont = 42;
  gStyle->SetTextFont(myfont);
  gStyle->SetLabelFont(myfont, "xyz");
  gStyle->SetLabelSize(0.025, "xyz");
  gStyle->SetTitleFont(myfont, "xyz");
  // gStyle->SetTitleFont(myfont, "p");
  gStyle->SetTitleSize(0.025, "xy");
  gStyle->SetTitleSize(0.016, "z");
  // gStyle->SetTitleSize(0.036, "p");
  gStyle->SetTitleOffset(1.5, "y");
  gStyle->SetTitleOffset(-1.36, "z");
  // gStyle->SetTitleOffset(-0.02, "p");

  gSystem->MakeDirectory("fig/evdisp");

  m_theApp = new TApplication("App", 0, 0);

#if ROOT_VERSION_CODE > ROOT_VERSION(6,4,0)
  // gStyle->SetPalette(kCool);
#endif
  gStyle->SetNumberContours(255);

  m_geometry = new TGeometry("evdisp", "K1.8 Event Display");

  ThreeVector worldSize(1000., 1000., 1000.); /*mm*/
  new TBRIK("world", "world", "void",
            worldSize.x(), worldSize.y(), worldSize.z());

  m_node = new TNode("node", "node", "world", 0., 0., 0.);
  m_geometry->GetNode("node")->SetVisibility(0);

#if BH2
  ConstructBH2();
#endif

  ConstructTarget();

#if KURAMA
  ConstructKURAMA();
#endif

#if BcOut
  ConstructBcOut();
#endif

#if SdcIn
  ConstructSdcIn();
#endif

#if SdcOut
  ConstructSdcOut();
#endif

#if SCH
  ConstructSCH();
#endif

#if TOF
  ConstructTOF();
#endif

#if WC
  ConstructWC();
#endif

  m_canvas = new TCanvas("canvas", "K1.8 Event Display",
                         1800, 900);
  m_canvas->Divide(2, 1);
  m_canvas->cd(1)->Divide(1, 2);
  m_canvas->cd(1)->cd(1)->SetPad(0.00, 0.92, 1.00, 1.00);
  m_canvas->cd(1)->cd(2)->SetPad(0.00, 0.00, 1.00, 0.92);
  // m_canvas->cd(1)->SetPad(0.00, 0.72, 1.00, 1.00);
  // m_canvas->cd(2)->SetPad(0.00, 0.00, 1.00, 0.72);
  m_canvas->cd(1)->cd(2);

  m_geometry->Draw();

  gPad->SetPhi(180.);
  gPad->SetTheta(-20.);
  gPad->GetView()->ZoomIn();

  m_canvas->cd(2);
  gPad->Update();

  m_canvas->Modified();
  m_canvas->Update();

#if Vertex

  m_canvas_vertex = new TCanvas("canvas_vertex", "K1.8 Event Display (Vertex)",
                                1000, 800);
  m_canvas_vertex->Divide(1,2);
  m_canvas_vertex->cd(1);
  gPad->DrawFrame(MinZ, MinX, MaxZ, MaxX, "Vertex XZ projection");

  m_TargetXZ_box->Draw("L");

  gPad->Update();

  m_canvas_vertex->cd(2);
  gPad->DrawFrame(MinZ, MinY, MaxZ, MaxY, "Vertex YZ projection");
  m_TargetYZ_box->Draw("L");
  gPad->Update();
#endif

#if Hist
  m_canvas_hist = new TCanvas("canvas_hist", "EventDisplay Hist",
                              400, 800);
  m_canvas_hist->Divide(1,3);
  m_hist_p  = new TH1F("hist1", "Momentum", 100, 0., 3.);
  m_hist_m2 = new TH1F("hist2", "Mass Square", 200, -0.5, 1.5);
  m_hist_missmass = new TH1F("hist3", "Missing Mass", 100, 1, 1.4);
  m_canvas_hist->cd(1)->SetGrid();
  m_hist_p->Draw();
  m_canvas_hist->cd(2)->SetGrid();
  m_hist_m2->Draw();
  m_canvas_hist->cd(3)->SetGrid();
  m_hist_missmass->Draw();
  // gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1111110);
#endif

#if Hist_Timing
  m_canvas_hist2 = new TCanvas("canvas_hist2", "EventDisplay Detector Timing",
                               800, 1000);
  m_canvas_hist2->Divide(3,3);
  m_hist_bh2  = new TH2F("hist_bh2", "BH2", NumOfSegBH2, 0., NumOfSegBH2, 500, -500, 500);
  m_hist_bh2->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_sch = new TH2F("hist_sch", "SCH", NumOfSegSCH, 0, NumOfSegSCH, 500, -500, 500);
  m_hist_sch->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_tof = new TH2F("hist_tof", "TOF", NumOfSegTOF, 0, NumOfSegTOF, 500, -500, 500);
  m_hist_tof->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_sdc1 = new TH2F("hist_sdc1", "SDC1", MaxWireSDC1, 0, MaxWireSDC1, 500, -500, 500);
  m_hist_sdc1->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_sdc1p = new TH2F("hist_sdc1p", "SDC1 Xp", MaxWireSDC1, 0, MaxWireSDC1, 500, -500, 500);

  m_hist_sdc1p->SetFillColor(kBlack);
  m_canvas_hist2->cd(1)->SetGrid();
  m_hist_bh2->Draw("box");
  m_canvas_hist2->cd(2)->SetGrid();
  m_hist_sch->Draw("box");
  m_canvas_hist2->cd(8)->SetGrid();
  m_hist_tof->Draw("box");
  m_canvas_hist2->cd(3)->SetGrid();
  m_hist_sdc1->Draw("box");
  m_hist_sdc1p->Draw("samebox");

  m_canvas_hist3 = new TCanvas("canvas_hist3", "EventDisplay Detector Timing",
                               800, 1000);
  m_canvas_hist3->Divide(3,2);

  m_hist_bc3 = new TH2F("hist_bc3", "BC3 X", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500);
  m_hist_bc3p = new TH2F("hist_bc3p", "BC3 Xp", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500);
  m_hist_bc3p->SetFillColor(kBlack);
  m_hist_bc3->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc3p->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_bc3_time = new TH1F("hist_bc3_time", "BC3 X", 500, -500, 500);
  m_hist_bc3p_time = new TH1F("hist_bc3p_time", "BC3 X", 500, -500, 500);
  m_hist_bc3p_time->SetFillStyle(3001);
  m_hist_bc3p_time->SetFillColor(kBlack);

  m_hist_bc3u = new TH2F("hist_bc3u", "BC3 U", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500);
  m_hist_bc3up = new TH2F("hist_bc3up", "BC3 Up", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500);
  m_hist_bc3up->SetFillColor(kBlack);
  m_hist_bc3u->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc3up->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_bc3v = new TH2F("hist_bc3v", "BC3 V", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500);
  m_hist_bc3vp = new TH2F("hist_bc3vp", "BC3 Vp", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500);
  m_hist_bc3vp->SetFillColor(kBlack);
  m_hist_bc3v->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc3vp->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_bc4 = new TH2F("hist_bc4", "BC4 X", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500);
  m_hist_bc4p = new TH2F("hist_bc4p", "BC4 Xp", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500);
  m_hist_bc4p->SetFillColor(kBlack);
  m_hist_bc4->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc4p->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_bc4_time = new TH1F("hist_bc4_time", "BC4 X", 500, -500, 500);
  m_hist_bc4p_time = new TH1F("hist_bc4p_time", "BC4 X", 500, -500, 500);
  m_hist_bc4p_time->SetFillStyle(3001);
  m_hist_bc4p_time->SetFillColor(kBlack);

  m_hist_bc4u = new TH2F("hist_bc4u", "BC4 U", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500);
  m_hist_bc4up = new TH2F("hist_bc4up", "BC4 Up", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500);
  m_hist_bc4up->SetFillColor(kBlack);
  m_hist_bc4u->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc4up->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_bc4v = new TH2F("hist_bc4v", "BC4 V", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500);
  m_hist_bc4vp = new TH2F("hist_bc4vp", "BC4 Vp", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500);
  m_hist_bc4vp->SetFillColor(kBlack);
  m_hist_bc4v->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc4vp->GetYaxis()->SetRangeUser(-100, 100);

  m_canvas_hist3->cd(1)->SetGrid();
  m_hist_bc3->Draw("box");
  m_hist_bc3p->Draw("samebox");
  m_canvas_hist3->cd(2)->SetGrid();
  m_hist_bc3v->Draw("box");
  m_hist_bc3vp->Draw("samebox");
  m_canvas_hist3->cd(3)->SetGrid();
  m_hist_bc3u->Draw("box");
  m_hist_bc3up->Draw("samebox");
  m_canvas_hist3->cd(4)->SetGrid();
  m_hist_bc4->Draw("box");
  m_hist_bc4p->Draw("samebox");
  m_canvas_hist3->cd(5)->SetGrid();
  m_hist_bc4v->Draw("box");
  m_hist_bc4vp->Draw("samebox");
  m_canvas_hist3->cd(6)->SetGrid();
  m_hist_bc4u->Draw("box");
  m_hist_bc4up->Draw("samebox");
  /*
    m_canvas_hist3->cd(3)->SetGrid();
    m_hist_bc3_time->Draw("");
    m_hist_bc3p_time->Draw("same");
    m_canvas_hist3->cd(4)->SetGrid();
    m_hist_bc4_time->Draw();
    m_hist_bc4p_time->Draw("same");
  */
  // gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

#endif

#if Hist_SdcOut
  m_canvas_hist4 = new TCanvas("canvas_hist4", "EventDisplay Detector Timing (SdcOut)", 800, 800);
  m_canvas_hist4->Divide(2,4);

  m_hist_sdc3_l = new TH2F("hist_sdc3_l", "SDC3 (leading)", MaxWireSDC3, 0, MaxWireSDC3, 500, -500, 500);
  m_hist_sdc3_t = new TH2F("hist_sdc3_t", "SDC3 (trailing)", MaxWireSDC3, 0, MaxWireSDC3, 500, -500, 500);
  m_hist_sdc3_t->SetFillColor(kRed);

  m_hist_sdc3p_l = new TH2F("hist_sdc3p_l", "SDC3 Xp (leading)", MaxWireSDC3, 0, MaxWireSDC3, 500, -500, 500);
  m_hist_sdc3p_t = new TH2F("hist_sdc3p_t", "SDC3 Xp(trailing)", MaxWireSDC3, 0, MaxWireSDC3, 500, -500, 500);
  m_hist_sdc3p_l->SetFillColor(kBlack);
  m_hist_sdc3p_t->SetFillColor(kGreen);

  m_hist_sdc3y_l = new TH2F("hist_sdc3y_l", "SDC3 Y (leading)", MaxWireSDC3, 0, MaxWireSDC3, 500, -500, 500);
  m_hist_sdc3y_t = new TH2F("hist_sdc3y_t", "SDC3 Y (trailing)", MaxWireSDC3, 0, MaxWireSDC3, 500, -500, 500);
  m_hist_sdc3y_t->SetFillColor(kRed);

  m_hist_sdc3yp_l = new TH2F("hist_sdc3yp_l", "SDC3 Yp (leading)", MaxWireSDC3, 0, MaxWireSDC3, 500, -500, 500);
  m_hist_sdc3yp_t = new TH2F("hist_sdc3yp_t", "SDC3 Yp(trailing)", MaxWireSDC3, 0, MaxWireSDC3, 500, -500, 500);
  m_hist_sdc3yp_l->SetFillColor(kBlack);
  m_hist_sdc3yp_t->SetFillColor(kGreen);


  m_hist_sdc4_l = new TH2F("hist_sdc4_l", "SDC4 (leading)", MaxWireSDC4X, 0, MaxWireSDC4X, 500, -500, 500);
  m_hist_sdc4_t = new TH2F("hist_sdc4_t", "SDC4 (trailing)", MaxWireSDC4X, 0, MaxWireSDC4X, 500, -500, 500);
  m_hist_sdc4_t->SetFillColor(kRed);

  m_hist_sdc4p_l = new TH2F("hist_sdc4p_l", "SDC4 Xp(leading)", MaxWireSDC4X, 0, MaxWireSDC4X, 500, -500, 500);
  m_hist_sdc4p_t = new TH2F("hist_sdc4p_t", "SDC4 Xp(trailing)", MaxWireSDC4X, 0, MaxWireSDC4X, 500, -500, 500);
  m_hist_sdc4p_l->SetFillColor(kBlack);
  m_hist_sdc4p_t->SetFillColor(kGreen);

  m_hist_sdc4y_l = new TH2F("hist_sdc4y_l", "SDC4 Y (leading)", MaxWireSDC4Y, 0, MaxWireSDC4Y, 500, -500, 500);
  m_hist_sdc4y_t = new TH2F("hist_sdc4y_t", "SDC4 Y (trailing)", MaxWireSDC4Y, 0, MaxWireSDC4Y, 500, -500, 500);
  m_hist_sdc4y_t->SetFillColor(kRed);

  m_hist_sdc4yp_l = new TH2F("hist_sdc4yp_l", "SDC4 Yp(leading)", MaxWireSDC4Y, 0, MaxWireSDC4Y, 500, -500, 500);
  m_hist_sdc4yp_t = new TH2F("hist_sdc4yp_t", "SDC4 Yp(trailing)", MaxWireSDC4Y, 0, MaxWireSDC4Y, 500, -500, 500);
  m_hist_sdc4yp_l->SetFillColor(kBlack);
  m_hist_sdc4yp_t->SetFillColor(kGreen);

  m_canvas_hist4->cd(5)->SetGrid();
  m_hist_sdc3_l->Draw("box");
  m_hist_sdc3_t->Draw("samebox");
  m_hist_sdc3p_l->Draw("samebox");
  m_hist_sdc3p_t->Draw("samebox");

  m_canvas_hist4->cd(6)->SetGrid();
  m_hist_sdc3y_l->Draw("box");
  m_hist_sdc3y_t->Draw("samebox");
  m_hist_sdc3yp_l->Draw("samebox");
  m_hist_sdc3yp_t->Draw("samebox");

  m_canvas_hist4->cd(7)->SetGrid();
  m_hist_sdc4_l->Draw("box");
  m_hist_sdc4_t->Draw("samebox");
  m_hist_sdc4p_l->Draw("samebox");
  m_hist_sdc4p_t->Draw("samebox");

  m_canvas_hist4->cd(8)->SetGrid();
  m_hist_sdc4y_l->Draw("box");
  m_hist_sdc4y_t->Draw("samebox");
  m_hist_sdc4yp_l->Draw("samebox");
  m_hist_sdc4yp_t->Draw("samebox");


#endif

#if Hist_BcIn
  m_canvas_hist5 = new TCanvas("canvas_hist5", "EventDisplay Detector Timing (BcIn)",
                               800, 1000);
  m_canvas_hist5->Divide(2,2);
  m_hist_bh1  = new TH2F("hist_bh1", "BH1", NumOfSegBH1, 0., NumOfSegBH1, 500, -500, 500);
  m_hist_bh1->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bft  = new TH2F("hist_bft", "BFT", NumOfSegBFT, 0., NumOfSegBFT, 500, -500, 500);
  m_hist_bft_p  = new TH2F("hist_bft_p", "BFT prime", NumOfSegBFT, 0., NumOfSegBFT, 500, -500, 500);
  m_hist_bft->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bft_p->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bft_p->SetFillColor(kBlack);

  m_hist_bcIn  = new TH2F("hist_bcIn", "BcIn Tracking", 200, -100, 100, 200, -150, 50);

  m_hist_bcOut  = new TH2F("hist_bcOut", "BcOut Tracking", 200, -100, 100, 200, 0, 600);


  m_canvas_hist5->cd(1)->SetGrid();
  m_hist_bh1->Draw("box");
  m_canvas_hist5->cd(2)->SetGrid();
  m_hist_bft->Draw("box");
  m_hist_bft_p->Draw("samebox");
  m_canvas_hist5->cd(3)->SetGrid();
  m_hist_bcIn->Draw("box");

  Double_t Bh1SegX[NumOfSegBH1] = {30./2., 20./2., 16./2., 12./2., 8./2., 8./2., 8./2., 12./2., 16./2., 20./2., 30./2.};
  Double_t Bh1SegY[NumOfSegBH1] = {5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2.};

  //Double_t localPosBh1Z = -515.;
  Double_t localPosBh1Z = gGeom.GetLocalZ(IdBH1);
  Double_t localPosBh1X_dX = 0.;
  Double_t localPosBh1X[NumOfSegBH1] = {-70. + localPosBh1X_dX,
                                        -46. + localPosBh1X_dX,
                                        -29. + localPosBh1X_dX,
                                        -16. + localPosBh1X_dX,
                                        -7. + localPosBh1X_dX,
                                        0. + localPosBh1X_dX,
                                        7. + localPosBh1X_dX,
                                        16. + localPosBh1X_dX,
                                        29. + localPosBh1X_dX,
                                        46. + localPosBh1X_dX,
                                        70. + localPosBh1X_dX};
  Double_t localPosBh1_dZ[NumOfSegBH1] = {4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5};

  for (Int_t i=0; i<NumOfSegBH1; i++) {
    m_BH1box_cont.push_back(new TBox(localPosBh1X[i]-Bh1SegX[i],
                                     localPosBh1Z+localPosBh1_dZ[i]-Bh1SegY[i],
                                     localPosBh1X[i]+Bh1SegX[i],
                                     localPosBh1Z+localPosBh1_dZ[i]+Bh1SegY[i]));
  }
  for (Int_t i=0; i<NumOfSegBH1; i++) {
    m_BH1box_cont[i]->SetFillColor(kWhite);
    m_BH1box_cont[i]->SetLineColor(kBlack);
    m_BH1box_cont[i]->Draw("L");
  }

  m_canvas_hist5->cd(4)->SetGrid();
  m_hist_bcOut->Draw("box");

  Double_t Bh2SegX[NumOfSegBH2] = {35./2., 10./2., 7./2., 7./2., 7./2., 7./2., 10./2., 35./2.};
  Double_t Bh2SegY[NumOfSegBH2] = {5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2.};

  Double_t localPosBh2X[NumOfSegBH2] = {-41.5, -19.0, -10.5, -3.5, 3.5, 10.5, 19.0, 41.5};

  Double_t localPosBh2Z = gGeom.GetLocalZ(IdBH2);

  for (Int_t i=0; i<NumOfSegBH2; i++) {
    m_BH2box_cont.push_back(new TBox(localPosBh2X[i]-Bh2SegX[i],
                                     localPosBh2Z-Bh2SegY[i],
                                     localPosBh2X[i]+Bh2SegX[i],
                                     localPosBh2Z+Bh2SegY[i]));
  }
  for (Int_t i=0; i<NumOfSegBH2; i++) {
    m_BH2box_cont[i]->SetFillColor(kWhite);
    m_BH2box_cont[i]->SetLineColor(kBlack);
    m_BH2box_cont[i]->Draw("L");
  }

  gStyle->SetOptStat(0);


  m_canvas_hist6 = new TCanvas("canvas_hist6", "EventDisplay Detector Timing (BcOut SdcIn)",
                               800, 1000);

  m_hist_sdcIn_predict  = new TH2F("hist_sdcIn_predict", "BcOut-SdcIn Tracking", 200, -400, 400, 750, -3000, 0);
  m_hist_sdcIn_predict->SetFillColor(kRed);
  m_hist_sdcIn_predict->SetLineColor(kRed);
  m_hist_sdcIn_predict->Draw("box");

  m_hist_sdcIn_predict2  = new TH2F("hist_sdcIn_predict2", "BcOut-SdcIn Tracking", 200, -400, 400, 750, -3000, 0);
  m_hist_sdcIn_predict2->SetFillColor(kMagenta);
  m_hist_sdcIn_predict2->SetLineColor(kMagenta);
  m_hist_sdcIn_predict2->Draw("samebox");

  m_hist_bcOut_sdcIn  = new TH2F("hist_bcOut_sdcIn", "BcOut-SdcIn Tracking", 400, -400, 400, 1500, -3000, 0);
  m_canvas_hist6->SetGrid();
  m_hist_bcOut_sdcIn->Draw("samebox");

  Double_t globalPosBh2Z = gGeom.GetGlobalPosition(IdBH2).z();
  Double_t globalPosBh2X = gGeom.GetGlobalPosition(IdBH2).x();

  for (Int_t i=0; i<NumOfSegBH2; i++) {
    m_BH2box_cont2.push_back(new TBox(globalPosBh2X+localPosBh2X[i]-Bh2SegX[i],
                                      globalPosBh2Z-Bh2SegY[i]-20,
                                      globalPosBh2X+localPosBh2X[i]+Bh2SegX[i],
                                      globalPosBh2Z+Bh2SegY[i]+20));
  }
  for (Int_t i=0; i<NumOfSegBH2; i++) {
    m_BH2box_cont2[i]->SetFillColor(kWhite);
    m_BH2box_cont2[i]->SetLineColor(kBlack);
    m_BH2box_cont2[i]->Draw("L");
  }

  Double_t globalPosSchZ = gGeom.GetGlobalPosition(IdSCH).z();

  Double_t SchSegX = 11.5/2.;
  Double_t SchSegY = 2.0/2.;
  for (Int_t i=1; i<=NumOfSegSCH; i++) {
    Double_t globalPosSchX = gGeom.CalcWirePosition(IdSCH, (Double_t)i);
    m_SCHbox_cont.push_back(new TBox(globalPosSchX-SchSegX,
                                     globalPosSchZ-SchSegY-20,
                                     globalPosSchX+SchSegX,
                                     globalPosSchZ+SchSegY+20));
  }
  for (Int_t i=0; i<NumOfSegSCH; i++) {
    m_SCHbox_cont[i]->SetFillColor(kWhite);
    m_SCHbox_cont[i]->SetLineColor(kBlack);
    m_SCHbox_cont[i]->Draw("L");
  }


  Double_t globalPosTarget_x = gGeom.GetGlobalPosition(IdTarget).x();
  Double_t globalPosTarget_y = gGeom.GetGlobalPosition(IdTarget).y();
  Double_t globalPosTarget_z = gGeom.GetGlobalPosition(IdTarget).z();
  Double_t target_r = 20.;
  Double_t target_z = 300./2;
  m_TargetXZ_box2 =  new TBox(globalPosTarget_x-target_r,
                              globalPosTarget_z-target_z,
                              globalPosTarget_x+target_r,
                              globalPosTarget_z+target_z);

  m_TargetXZ_box2->SetFillColor(kWhite);
  m_TargetXZ_box2->SetLineColor(kBlack);
  m_TargetXZ_box2->Draw("L");

  m_TargetYZ_box2 =  new TBox(globalPosTarget_y-target_r,
                              globalPosTarget_z-target_z,
                              globalPosTarget_y+target_r,
                              globalPosTarget_z+target_z);

  m_TargetYZ_box2->SetFillColor(kWhite);
  m_TargetYZ_box2->SetLineColor(kBlack);
  m_TargetYZ_box2->Draw("L");

#endif

  ResetVisibility();

  m_is_ready = true;
  return m_is_ready;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructBH2()
{
  const Int_t lid = gGeom.GetDetectorId("BH2");

  Double_t rotMatBH2[9] = {};
  Double_t BH2wallX = 120.0/2.; // X
  Double_t BH2wallY =   6.0/2.; // Z
  Double_t BH2wallZ =  40.0/2.; // Y
  Double_t BH2SizeX[NumOfSegBH2] = { 120./2. }; // X
  Double_t BH2SizeY[NumOfSegBH2] = {   6./2. }; // Z
  Double_t BH2SizeZ[NumOfSegBH2] = {  40./2. }; // Y
  Double_t BH2PosX[NumOfSegBH2]  = { 0./2. };
  Double_t BH2PosY[NumOfSegBH2]  = { 0./2. };
  Double_t BH2PosZ[NumOfSegBH2]  = { 0./2. };

  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatBH2);

  new TRotMatrix("rotBH2", "rotBH2", rotMatBH2);
  const ThreeVector& BH2wallPos = gGeom.GetGlobalPosition(lid);
  new TBRIK("BH2wall_brik", "BH2wall_brik", "void",
            BH2wallX, BH2wallY, BH2wallZ);
  m_BH2wall_node = new TNode("BH2wall_node", "BH2wall_node", "BH2wall_brik",
                             BH2wallPos.x(),
                             BH2wallPos.y(),
                             BH2wallPos.z(), "rotBH2", "void");
  m_BH2wall_node->SetVisibility(0);
  m_BH2wall_node->cd();

  for(Int_t i=0; i<NumOfSegBH2; ++i){
    new TBRIK(Form("BH2seg_brik_%d", i),
              Form("BH2seg_brik_%d", i),
              "void", BH2SizeX[i], BH2SizeY[i], BH2SizeZ[i]);
    m_BH2seg_node.push_back(new TNode(Form("BH2seg_node_%d", i),
                                      Form("BH2seg_node_%d", i),
                                      Form("BH2seg_brik_%d", i),
                                      BH2PosX[i], BH2PosY[i], BH2PosZ[i]));
  }
  m_node->cd();
  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructKURAMA()
{
  Double_t Matrix[9] = {};

  Double_t inner_x = 1400.0/2.; // X
  Double_t inner_y =  800.0/2.; // Z
  Double_t inner_z =  800.0/2.; // Y

  Double_t outer_x = 2200.0/2.; // X
  Double_t outer_y =  800.0/2.; // Z
  Double_t outer_z = 1540.0/2.; // Y

  Double_t uguard_inner_x = 1600.0/2.; // X
  Double_t uguard_inner_y =  100.0/2.; // Z
  Double_t uguard_inner_z = 1940.0/2.; // Y
  Double_t uguard_outer_x =  600.0/2.; // X
  Double_t uguard_outer_y =  100.0/2.; // Z
  Double_t uguard_outer_z =  300.0/2.; // Y

  Double_t dguard_inner_x = 1600.0/2.; // X
  Double_t dguard_inner_y =  100.0/2.; // Z
  Double_t dguard_inner_z = 1940.0/2.; // Y
  Double_t dguard_outer_x = 1100.0/2.; // X
  Double_t dguard_outer_y =  100.0/2.; // Z
  Double_t dguard_outer_z = 1100.0/2.; // Y

  CalcRotMatrix(0., 0., 0., Matrix);

  new TRotMatrix("rotKURAMA", "rotKURAMA", Matrix);

  new TBRIK("kurama_inner_brik", "kurama_inner_brik",
            "void", inner_x, inner_y, inner_z);

  new TBRIK("kurama_outer_brik", "kurama_outer_brik",
            "void", outer_x, outer_y, outer_z);

  new TBRIK("uguard_inner_brik", "uguard_inner_brik",
            "void", uguard_inner_x, uguard_inner_y, uguard_inner_z);

  new TBRIK("uguard_outer_brik", "uguard_outer_brik",
            "void", uguard_outer_x, uguard_outer_y, uguard_outer_z);

  new TBRIK("dguard_inner_brik", "dguard_inner_brik",
            "void", dguard_inner_x, dguard_inner_y, dguard_inner_z);

  new TBRIK("dguard_outer_brik", "dguard_outer_brik",
            "void", dguard_outer_x, dguard_outer_y, dguard_outer_z);


  m_kurama_inner_node = new TNode("kurama_inner_node",
                                  "kurama_inner_node",
                                  "kurama_inner_brik",
                                  0., 0., 0., "rotKURAMA", "void");
  m_kurama_outer_node = new TNode("kurama_outer_node",
                                  "kurama_outer_node",
                                  "kurama_outer_brik",
                                  0., 0., 0., "rotKURAMA", "void");

  TNode *uguard_inner = new TNode("uguard_inner_node",
                                  "uguard_inner_node",
                                  "uguard_inner_brik",
                                  0., 0., -820.,
                                  "rotKURAMA", "void");
  TNode *uguard_outer = new TNode("uguard_outer_node",
                                  "uguard_outer_node",
                                  "uguard_outer_brik",
                                  0., 0., -820.,
                                  "rotKURAMA", "void");

  TNode *dguard_inner = new TNode("dguard_inner_node",
                                  "dguard_inner_node",
                                  "dguard_inner_brik",
                                  0., 0., 820.,
                                  "rotKURAMA", "void");
  TNode *dguard_outer = new TNode("dguard_outer_node",
                                  "dguard_outer_node",
                                  "dguard_outer_brik",
                                  0., 0., 820.,
                                  "rotKURAMA", "void");

  const Color_t color = kBlack;

  m_kurama_inner_node->SetLineColor(color);
  m_kurama_outer_node->SetLineColor(color);
  uguard_inner->SetLineColor(color);
  uguard_outer->SetLineColor(color);
  dguard_inner->SetLineColor(color);
  dguard_outer->SetLineColor(color);

  m_node->cd();
  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructBcOut()
{
  const Double_t wireL = 200.0;

  const Double_t offsetZ = -3068.4;

  // BC3 X1
  {
    const Int_t lid = gGeom.GetDetectorId("BC3-X1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC3X1", "rotBC3X1", Matrix);
    new TTUBE("BC3X1Tube", "BC3X1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC3x1_node.push_back(new TNode(Form("BC3x1_Node_%d", wire),
                                       Form("BC3x1_Node_%d", wire),
                                       "BC3X1Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC3X1", "void"));
    }
  }

  // BC3 X2
  {
    const Int_t lid = gGeom.GetDetectorId("BC3-X2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC3X2", "rotBC3X2", Matrix);
    new TTUBE("BC3X2Tube", "BC3X2Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC3x2_node.push_back(new TNode(Form("BC3x2_Node_%d", wire),
                                       Form("BC3x2_Node_%d", wire),
                                       "BC3X2Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC3X2", "void"));
    }
  }

  // BC3 V1
  {
    const Int_t lid = gGeom.GetDetectorId("BC3-V1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC3V1", "rotBC3V1", Matrix);
    new TTUBE("BC3V1Tube", "BC3V1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC3v1_node.push_back(new TNode(Form("BC3v1_Node_%d", wire),
                                       Form("BC3v1_Node_%d", wire),
                                       "BC3V1Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC3V1", "void"));
    }
  }

  // BC3 V2
  {
    const Int_t lid = gGeom.GetDetectorId("BC3-V2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC3V2", "rotBC3V2", Matrix);
    new TTUBE("BC3V2Tube", "BC3V2Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC3v2_node.push_back(new TNode(Form("BC3v2_Node_%d", wire),
                                       Form("BC3v2_Node_%d", wire),
                                       "BC3V2Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC3V2", "void"));
    }
  }

  // BC3 U1
  {
    const Int_t lid = gGeom.GetDetectorId("BC3-U1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC3U1", "rotBC3U1", Matrix);
    new TTUBE("BC3U1Tube", "BC3U1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC3u1_node.push_back(new TNode(Form("BC3u1_Node_%d", wire),
                                       Form("BC3u1_Node_%d", wire),
                                       "BC3U1Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC3U1", "void"));
    }
  }

  // BC3 U2
  {
    const Int_t lid = gGeom.GetDetectorId("BC3-U2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC3U2", "rotBC3U2", Matrix);
    new TTUBE("BC3U2Tube", "BC3U2Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC3u2_node.push_back(new TNode(Form("BC3u2_Node_%d", wire),
                                       Form("BC3u2_Node_%d", wire),
                                       "BC3U2Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC3U2", "void"));
    }
  }

  // BC4 U1
  {
    const Int_t lid = gGeom.GetDetectorId("BC4-U1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC4U1", "rotBC4U1", Matrix);
    new TTUBE("BC4U1Tube", "BC4U1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC4; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC4u1_node.push_back(new TNode(Form("BC4u1_Node_%d", wire),
                                       Form("BC4u1_Node_%d", wire),
                                       "BC4U1Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC4U1", "void"));
    }
  }

  // BC4 U2
  {
    const Int_t lid = gGeom.GetDetectorId("BC4-U2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC4U2", "rotBC4U2", Matrix);
    new TTUBE("BC4U2Tube", "BC4U2Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC4; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC4u2_node.push_back(new TNode(Form("BC4u2_Node_%d", wire),
                                       Form("BC4u2_Node_%d", wire),
                                       "BC4U2Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC4U2", "void"));
    }
  }

  // BC4 V1
  {
    const Int_t lid = gGeom.GetDetectorId("BC4-V1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC4V1", "rotBC4V1", Matrix);
    new TTUBE("BC4V1Tube", "BC4V1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC4; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC4v1_node.push_back(new TNode(Form("BC4v1_Node_%d", wire),
                                       Form("BC4v1_Node_%d", wire),
                                       "BC4V1Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC4V1", "void"));
    }
  }

  // BC4 V2
  {
    const Int_t lid = gGeom.GetDetectorId("BC4-V2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC4V2", "rotBC4V2", Matrix);
    new TTUBE("BC4V2Tube", "BC4V2Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC4; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC4v2_node.push_back(new TNode(Form("BC4v2_Node_%d", wire),
                                       Form("BC4v2_Node_%d", wire),
                                       "BC4V2Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC4V2", "void"));
    }
  }

  // BC4 X1
  {
    const Int_t lid = gGeom.GetDetectorId("BC4-X1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC4X1", "rotBC4X1", Matrix);
    new TTUBE("BC4X1Tube", "BC4X1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC4; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC4x1_node.push_back(new TNode(Form("BC4x1_Node_%d", wire),
                                       Form("BC4x1_Node_%d", wire),
                                       "BC4X1Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC4X1", "void"));
    }
  }

  // BC4 X2
  {
    const Int_t lid = gGeom.GetDetectorId("BC4-X2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC4X2", "rotBC4X2", Matrix);
    new TTUBE("BC4X2Tube", "BC4X2Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC4; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC4x2_node.push_back(new TNode(Form("BC4x2_Node_%d", wire),
                                       Form("BC4x2_Node_%d", wire),
                                       "BC4X2Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC4X2", "void"));
    }
  }

  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructSdcIn()
{
  const Double_t wireLSDC1 = 200.0;
  const Double_t wireLSDC2X = 400.0;
  const Double_t wireLSDC2Y = 700.0;

  // SDC1 V1
  {
    const Int_t lid = gGeom.GetDetectorId("SDC1-V1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireLSDC1/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotSDC1V1", "rotSDC1V1", Matrix);
    new TTUBE("SDC1V1Tube", "SDC1V1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireSDC1; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_SDC1v1_node.push_back(new TNode(Form("SDC1v1_Node_%d", wire),
                                        Form("SDC1v1_Node_%d", wire),
                                        "SDC1V1Tube",
                                        wireGlobalPos.x()+localPos,
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotSDC1V1", "void"));
    }
  }

  // SDC1 V2
  {
    const Int_t lid = gGeom.GetDetectorId("SDC1-V2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t Z    = wireLSDC1/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotV2", "rotV2", Matrix);
    new TTUBE("SDC1V2Tube", "SDC1V2Tube", "void", Rmin, Rmax, Z);
    for(Int_t wire=1; wire<=MaxWireSDC1; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_SDC1v2_node.push_back(new TNode(Form("SDC1v2_Node_%d", wire),
                                        Form("SDC1v2_Node_%d", wire),
                                        "SDC1V2Tube",
                                        wireGlobalPos.x()+localPos,
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotV2", "void"));
    }
  }

  // SDC1 X1
  {
    const Int_t lid = gGeom.GetDetectorId("SDC1-X1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t Z    = wireLSDC1/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotX1", "rotX1", Matrix);
    new TTUBE("SDC1X1Tube", "SDC1X1Tube", "void", Rmin, Rmax, Z);
    for(Int_t wire=1; wire<=MaxWireSDC1; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_SDC1x1_node.push_back(new TNode(Form("SDC1x1_Node_%d", wire),
                                        Form("SDC1x1_Node_%d", wire),
                                        "SDC1X1Tube",
                                        wireGlobalPos.x()+localPos,
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotX1", "void"));
    }
  }

  // SDC1 X2
  {
    const Int_t lid = gGeom.GetDetectorId("SDC1-X2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t Z    = wireLSDC1/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotX2", "rotX2", Matrix);
    new TTUBE("SDC1X2Tube", "SDC1X2Tube", "void", Rmin, Rmax, Z);
    for(Int_t wire=1; wire<=MaxWireSDC1; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_SDC1x2_node.push_back(new TNode(Form("SDC1x2_Node_%d", wire),
                                        Form("SDC1x2_Node_%d", wire),
                                        "SDC1X2Tube",
                                        wireGlobalPos.x()+localPos,
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotX2", "void"));
    }
  }

  // SDC1 U1
  {
    const Int_t lid = gGeom.GetDetectorId("SDC1-U1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t Z    = wireLSDC1/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotU1", "rotU1", Matrix);
    new TTUBE("SDC1U1Tube", "SDC1U1Tube", "void", Rmin, Rmax, Z);
    for(Int_t wire=1; wire<=MaxWireSDC1; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_SDC1u1_node.push_back(new TNode(Form("SDC1u1_Node_%d", wire),
                                        Form("SDC1u1_Node_%d", wire),
                                        "SDC1U1Tube",
                                        wireGlobalPos.x()+localPos,
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotU1", "void"));
    }
  }

  // SDC1 U2
  {
    const Int_t lid = gGeom.GetDetectorId("SDC1-U2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t Z    = wireLSDC1/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotU2", "rotU2", Matrix);
    new TTUBE("SDC1U2Tube", "SDC1U2Tube", "void", Rmin, Rmax, Z);
    for(Int_t wire=1; wire<=MaxWireSDC1; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_SDC1u2_node.push_back(new TNode(Form("SDC1u2_Node_%d", wire),
                                        Form("SDC1u2_Node_%d", wire),
                                        "SDC1U2Tube",
                                        wireGlobalPos.x()+localPos,
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotU2", "void"));
    }
  }

  // SDC2 X1
  {
    const Int_t lid = gGeom.GetDetectorId("SDC2-X1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireLSDC2X/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotSDC2X1", "rotSDC2X1", Matrix);
    new TTUBE("SDC2X1Tube", "SDC2X1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireSDC2X; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0, 0);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC2x1_node.push_back(new TNode(Form("SDC2x1_Node_%d", wire),
                                        Form("SDC2x1_Node_%d", wire),
                                        "SDC2X1Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotSDC2X1", "void"));
    }
  }

  // SDC2 X2
  {
    const Int_t lid = gGeom.GetDetectorId("SDC2-X2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t Z    = wireLSDC2X/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotSDC2X2", "rotSDC2X2", Matrix);
    new TTUBE("SDC2X2Tube", "SDC2X2Tube", "void", Rmin, Rmax, Z);
    for(Int_t wire=1; wire<=MaxWireSDC2X; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0, 0);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC2x2_node.push_back(new TNode(Form("SDC2x2_Node_%d", wire),
                                        Form("SDC2x2_Node_%d", wire),
                                        "SDC2X2Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotSDC2X2", "void"));
    }
  }

  // SDC2 Y1
  {
    const Int_t lid = gGeom.GetDetectorId("SDC2-Y1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t Z    = wireLSDC2Y/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotSDC2Y1", "rotSDC2Y1", Matrix);
    new TTUBE("SDC2Y1Tube", "SDC2Y1Tube", "void", Rmin, Rmax, Z);
    for(Int_t wire=1; wire<=MaxWireSDC2Y; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0, 0);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC2y1_node.push_back(new TNode(Form("SDC2y1_Node_%d", wire),
                                        Form("SDC2y1_Node_%d", wire),
                                        "SDC2Y1Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotSDC2Y1", "void"));
    }
  }

  // SDC2 Y2
  {
    const Int_t lid = gGeom.GetDetectorId("SDC2-Y2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t Z    = wireLSDC2Y/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotSDC2Y2", "rotSDC2Y2", Matrix);
    new TTUBE("SDC2Y2Tube", "SDC2Y2Tube", "void", Rmin, Rmax, Z);
    for(Int_t wire=1; wire<=MaxWireSDC2Y; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0, 0);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC2y2_node.push_back(new TNode(Form("SDC2y2_Node_%d", wire),
                                        Form("SDC2y2_Node_%d", wire),
                                        "SDC2Y2Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotSDC2Y2", "void"));
    }
  }

  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructSdcOut()
{
  const Double_t wireLSDC3  = 1152.0;
  const Double_t wireLSDC4Y = 1920.0;
  const Double_t wireLSDC4X = 1280.0;

  // SDC3 X1
  {
    const Int_t lid = gGeom.GetDetectorId("SDC3-X1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireLSDC3/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotX1", "rotX1", Matrix);
    new TTUBE("SDC3X1Tube", "SDC3X1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<= MaxWireSDC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC3x1_node.push_back(new TNode(Form("SDC3x1_Node_%d", wire),
                                        Form("SDC3x1_Node_%d", wire),
                                        "SDC3X1Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotX1", "void"));
    }
  }

  // SDC3 X2
  {
    const Int_t lid = gGeom.GetDetectorId("SDC3-X2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireLSDC3/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotX2", "rotX2", Matrix);
    new TTUBE("SDC3X2Tube", "SDC3X2Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<= MaxWireSDC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC3x2_node.push_back(new TNode(Form("SDC3x2_Node_%d", wire),
                                        Form("SDC3x2_Node_%d", wire),
                                        "SDC3X2Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotX2", "void"));
    }
  }

  // SDC3 Y1
  {
    const Int_t lid = gGeom.GetDetectorId("SDC3-Y1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireLSDC3/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotY1", "rotY1", Matrix);
    new TTUBE("SDC3Y1Tube", "SDC3Y1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<= MaxWireSDC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC3y1_node.push_back(new TNode(Form("SDC3y1_Node_%d", wire),
                                        Form("SDC3y1_Node_%d", wire),
                                        "SDC3Y1Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotY1", "void"));
    }
  }

  // SDC3 Y2
  {
    const Int_t lid = gGeom.GetDetectorId("SDC3-Y2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireLSDC3/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotY2", "rotY2", Matrix);
    new TTUBE("SDC3Y2Tube", "SDC3Y2Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<= MaxWireSDC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC3y2_node.push_back(new TNode(Form("SDC3y2_Node_%d", wire),
                                        Form("SDC3y2_Node_%d", wire),
                                        "SDC3Y2Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotY2", "void"));
    }
  }

  // SDC4 Y1
  {
    const Int_t lid = gGeom.GetDetectorId("SDC4-Y1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireLSDC4Y/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid), Matrix);
    new TRotMatrix("rotY1", "rotY1", Matrix);
    new TTUBE("SDC4Y1Tube", "SDC4Y1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<= MaxWireSDC4Y; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC4y1_node.push_back(new TNode(Form("SDC4y1_Node_%d", wire),
                                        Form("SDC4y1_Node_%d", wire),
                                        "SDC4Y1Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotY1", "void"));
    }
  }

  // SDC4 Y2
  {
    const Int_t lid = gGeom.GetDetectorId("SDC4-Y2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireLSDC4Y/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotY2", "rotY2", Matrix);
    new TTUBE("SDC4Y2Tube", "SDC4Y2Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<= MaxWireSDC4Y; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC4y2_node.push_back(new TNode(Form("SDC4y2_Node_%d", wire),
                                        Form("SDC4y2_Node_%d", wire),
                                        "SDC4Y2Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotY2", "void"));
    }
  }

  // SDC4 X1
  {
    const Int_t lid = gGeom.GetDetectorId("SDC4-X1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireLSDC4X/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotX1", "rotX1", Matrix);
    new TTUBE("SDC4X1Tube", "SDC4X1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<= MaxWireSDC4X; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC4x1_node.push_back(new TNode(Form("SDC4x1_Node_%d", wire),
                                        Form("SDC4x1_Node_%d", wire),
                                        "SDC4X1Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotX1", "void"));
    }
  }

  // SDC4 X2
  {
    const Int_t lid = gGeom.GetDetectorId("SDC4-X2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireLSDC4X/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotX2", "rotX2", Matrix);
    new TTUBE("SDC4X2Tube", "SDC4X2Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<= MaxWireSDC4X; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC4x2_node.push_back(new TNode(Form("SDC4x2_Node_%d", wire),
                                        Form("SDC4x2_Node_%d", wire),
                                        "SDC4X2Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotX2", "void"));
    }
  }

  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructTarget()
{
  Double_t TargetX = 50.0/2.0;
  Double_t TargetY = 30.0/2.0; // Z
  Double_t TargetZ = 30.0/2.0; // Y
  Double_t rotMatTarget[9];

  new TBRIK("target_brik", "target_brik", "void",
            TargetX, TargetY, TargetZ);

  CalcRotMatrix(0., 0., 0., rotMatTarget);
  new TRotMatrix("rotTarget", "rotTarget", rotMatTarget);

  const Int_t lid = IdTarget;
  ThreeVector GlobalPos = gGeom.GetGlobalPosition(lid);
  m_target_node = new TNode("target_node", "target_node", "target_brik",
                            GlobalPos.x(), GlobalPos.y(), GlobalPos.z(),
                            "rotTarget", "void");

#if Vertex
  m_TargetXZ_box = new TBox(-TargetY, -TargetX, TargetY, TargetX);
  m_TargetXZ_box->SetFillColor(kWhite);
  m_TargetXZ_box->SetLineColor(kBlack);
  m_TargetYZ_box = new TBox(-TargetY, -TargetZ, TargetY, TargetZ);
  m_TargetYZ_box->SetFillColor(kWhite);
  m_TargetYZ_box->SetLineColor(kBlack);
#endif

  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructSCH()
{
  const Int_t lid = gGeom.GetDetectorId("SCH");

  Double_t rotMatSCH[9] = {};
  Double_t SCHwallX =   11.5/2.0*NumOfSegSCH; // X
  Double_t SCHwallY =    2.0/2.0*2.0; // Z
  Double_t SCHwallZ =  450.0/2.0; // Y

  Double_t SCHSegX =   11.5/2.0; // X
  Double_t SCHSegY =    2.0/2.0; // Z
  Double_t SCHSegZ =  450.0/2.0; // Y

  Double_t overlap = 1.0;

  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatSCH);

  new TRotMatrix("rotSCH", "rotSCH", rotMatSCH);
  ThreeVector  SCHwallPos = gGeom.GetGlobalPosition(lid);
  Double_t offset =  gGeom.CalcWirePosition(lid, (Double_t)NumOfSegSCH/2.-0.5);
  new TBRIK("SCHwall_brik", "SCHwall_brik", "void",
            SCHwallX, SCHwallY, SCHwallZ);
  m_SCHwall_node = new TNode("SCHwall_node", "SCHwall_node", "SCHwall_brik",
                             SCHwallPos.x() + offset,
                             SCHwallPos.y(),
                             SCHwallPos.z(),
                             "rotSCH", "void");

  m_SCHwall_node->SetVisibility(0);
  m_SCHwall_node->cd();


  new TBRIK("SCHseg_brik", "SCHseg_brik", "void",
            SCHSegX, SCHSegY, SCHSegZ);
  for(Int_t i=0; i<NumOfSegSCH; i++){
    ThreeVector schSegLocalPos((-NumOfSegSCH/2.+i)*(SCHSegX*2.-overlap)+10.5/2.,
                               (-(i%2)*2+1)*SCHSegY-(i%2)*2+1,
                               0.);
    m_SCHseg_node.push_back(new TNode(Form("SCHseg_node_%d", i),
                                      Form("SCHseg_node_%d", i),
                                      "SCHseg_brik",
                                      schSegLocalPos.x(),
                                      schSegLocalPos.y(),
                                      schSegLocalPos.z()));
  }

  m_node->cd();
  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructTOF()
{
  const Int_t lid = gGeom.GetDetectorId("TOF");

  Double_t rotMatTOF[9] = {};
  Double_t TOFwallX =   80.0/2.0*NumOfSegTOF; // X
  Double_t TOFwallY =   30.0/2.0*2.0; // Z
  Double_t TOFwallZ = 1800.0/2.0; // Y

  Double_t TOFSegX =   80.0/2.0; // X
  Double_t TOFSegY =   30.0/2.0; // Z
  Double_t TOFSegZ = 1800.0/2.0; // Y

  Double_t overlap = 5.0;

  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatTOF);

  new TRotMatrix("rotTOF", "rotTOF", rotMatTOF);
  ThreeVector  TOFwallPos = gGeom.GetGlobalPosition(lid);
  // Double_t offset = gGeom.CalcWirePosition(lid, (Double_t)NumOfSegTOF/2.-0.5);
  new TBRIK("TOFwall_brik", "TOFwall_brik", "void",
            TOFwallX, TOFwallY, TOFwallZ);
  m_TOFwall_node = new TNode("TOFwall_node", "TOFwall_node", "TOFwall_brik",
                             TOFwallPos.x(),// + offset,
                             TOFwallPos.y(),
                             TOFwallPos.z(),
                             "rotTOF", "void");

  m_TOFwall_node->SetVisibility(0);
  m_TOFwall_node->cd();

  new TBRIK("TOFseg_brik", "TOFseg_brik", "void",
            TOFSegX, TOFSegY, TOFSegZ);
  for(Int_t i=0; i<NumOfSegTOF; i++){
    ThreeVector tofSegLocalPos((-NumOfSegTOF/2.+i)*(TOFSegX*2.-overlap)+75./2.,
                               (-(i%2)*2+1)*TOFSegY-(i%2)*2+1,
                               0.);
    std::cout << i << " " << tofSegLocalPos << std::endl;
    m_TOFseg_node.push_back(new TNode(Form("TOFseg_node_%d", i),
                                      Form("TOFseg_node_%d", i),
                                      "TOFseg_brik",
                                      tofSegLocalPos.x(),
                                      -tofSegLocalPos.y(),
                                      tofSegLocalPos.z()));
  }

  m_node->cd();
  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructWC()
{
  const Int_t lid = gGeom.GetDetectorId("WC");

  Double_t rotMatWC[9] = {};
  Double_t WCwallX =  255.0/2.0*NumOfSegWC; // X
  Double_t WCwallY =  205.0/2.0*2.0; // Z
  Double_t WCwallZ = 1840.0/2.0; // Y

  Double_t WCSegX =  255.0/2.0; // X
  Double_t WCSegY =  205.0/2.0; // Z
  Double_t WCSegZ = 1840.0/2.0; // Y

  Double_t overlap = 255.0/2.0;

  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatWC);

  new TRotMatrix("rotWC", "rotWC", rotMatWC);
  ThreeVector  WCwallPos = gGeom.GetGlobalPosition(lid);
  // Double_t offset = gGeom.CalcWirePosition(lid, (Double_t)NumOfSegWC/2.-0.5);
  new TBRIK("WCwall_brik", "WCwall_brik", "void",
            WCwallX, WCwallY, WCwallZ);
  m_WCwall_node = new TNode("WCwall_node", "WCwall_node", "WCwall_brik",
                             WCwallPos.x(),// + offset,
                             WCwallPos.y(),
                             WCwallPos.z(),
                             "rotWC", "void");

  m_WCwall_node->SetVisibility(0);
  m_WCwall_node->cd();

  new TBRIK("WCseg_brik", "WCseg_brik", "void",
            WCSegX, WCSegY, WCSegZ);
  for(Int_t i=0; i<NumOfSegWC; i++){
    ThreeVector wcSegLocalPos((-NumOfSegWC/2.+i)*(WCSegX*2.-overlap)+255./4.,
                               (-(i%2)*2+1)*WCSegY-(i%2)*2+1,
                               0.);
    std::cout << i << " " << wcSegLocalPos << std::endl;
    m_WCseg_node.push_back(new TNode(Form("WCseg_node_%d", i),
                                      Form("WCseg_node_%d", i),
                                      "WCseg_brik",
                                      wcSegLocalPos.x(),
                                      -wcSegLocalPos.y(),
                                      wcSegLocalPos.z()));
  }

  m_node->cd();
  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
void
EventDisplay::DrawInitTrack()
{
  m_canvas->cd(1)->cd(2);
  if(m_init_step_mark) m_init_step_mark->Draw();

  m_canvas->Update();

}

//_____________________________________________________________________________
void
EventDisplay::DrawHitWire(Int_t lid, Int_t hit_wire,
                          Bool_t range_check, Bool_t tdc_check)
{
  if(hit_wire<=0) return;

  TString node_name;

  const TString bcout_node_name[NumOfLayersBcOut]
    = { Form("BC3x1_Node_%d", hit_wire),
        Form("BC3x2_Node_%d", hit_wire),
        Form("BC3u1_Node_%d", hit_wire),
        Form("BC3u2_Node_%d", hit_wire),
        Form("BC3v1_Node_%d", hit_wire),
        Form("BC3v2_Node_%d", hit_wire),
        Form("BC4x1_Node_%d", hit_wire),
        Form("BC4x2_Node_%d", hit_wire),
        Form("BC4u1_Node_%d", hit_wire),
        Form("BC4u2_Node_%d", hit_wire),
        Form("BC4v1_Node_%d", hit_wire),
        Form("BC4v2_Node_%d", hit_wire) };

  const TString sdcin_node_name[NumOfLayersSdcIn]
    = { Form("SDC1v1_Node_%d", hit_wire),
        Form("SDC1v2_Node_%d", hit_wire),
        Form("SDC1x1_Node_%d", hit_wire),
        Form("SDC1x2_Node_%d", hit_wire),
        Form("SDC1u1_Node_%d", hit_wire),
        Form("SDC1u2_Node_%d", hit_wire),
        Form("SDC2x1_Node_%d", hit_wire),
        Form("SDC2x2_Node_%d", hit_wire),
        Form("SDC2y1_Node_%d", hit_wire),
        Form("SDC2y2_Node_%d", hit_wire) };

  const TString sdcout_node_name[NumOfLayersSdcOut]
    = { Form("SDC3x1_Node_%d", hit_wire),
        Form("SDC3x2_Node_%d", hit_wire),
        Form("SDC3y1_Node_%d", hit_wire),
        Form("SDC3y2_Node_%d", hit_wire),
        Form("SDC4y1_Node_%d", hit_wire),
        Form("SDC4y2_Node_%d", hit_wire),
        Form("SDC4x1_Node_%d", hit_wire),
        Form("SDC4x2_Node_%d", hit_wire) };

  switch (lid) {

    // SDC1
  case 1: case 2: case 3: case 4: case 5: case 6:
    if(hit_wire>MaxWireSDC1) return;
    node_name = sdcin_node_name[lid-1];
    break;

    // SDC2X
  case 7: case 8:
    if(hit_wire>MaxWireSDC2X) return;
    node_name = sdcin_node_name[lid-1];
    break;

    // SDC2Y
  case 9: case 10:
    if(hit_wire>MaxWireSDC2Y) return;
    node_name = sdcin_node_name[lid-1];
    break;


    // SDC3
  case 31: case 32: case 33: case 34:
    if(hit_wire>MaxWireSDC3) return;
    node_name = sdcout_node_name[lid-31];
    break;

    // SDC4Y
  case 35: case 36:
    if(hit_wire>MaxWireSDC4Y) return;
    node_name = sdcout_node_name[lid-31];
    break;

    // SDC4X
  case 37: case 38:
    if(hit_wire>MaxWireSDC4X) return;
    node_name = sdcout_node_name[lid-31];
    break;

    // BC3
  case 113: case 114: case 115: case 116: case 117: case 118:
    if(hit_wire>MaxWireBC4) return;
    node_name = bcout_node_name[lid-113];
    break;

    // BC4
  case 119: case 120: case 121: case 122: case 123: case 124:
    if(hit_wire>MaxWireBC3) return;
    node_name = bcout_node_name[lid-113];
    break;

  default:
    throw Exception(FUNC_NAME + Form(" no such plane : %d", lid));
  }

  auto node = m_geometry->GetNode(node_name);
  if(!node){
    throw Exception(FUNC_NAME + Form(" no such node : %s", node_name.Data()));
  }

  node->SetVisibility(1);
  if(range_check && tdc_check)
    node->SetLineColor(kBlue);
  else if(range_check && !tdc_check)
    node->SetLineColor(28);
  else
    node->SetLineColor(kBlack);

  // m_canvas->cd(2);
  // m_geometry->Draw();
  m_canvas->Update();
}

//_____________________________________________________________________________
void
EventDisplay::DrawHitHodoscope(Int_t lid, Int_t seg, Int_t Tu, Int_t Td)
{
  TString node_name;
  if(seg<0) return;

  if(lid == IdBH2){
    node_name = Form("BH2seg_node_%d", seg);
  }else if(lid == IdSCH){
    node_name = Form("SCHseg_node_%d", seg);
  }else if(lid == IdTOF){
    node_name = Form("TOFseg_node_%d", seg);
  }else if(lid == IdWC){
    node_name = Form("WCseg_node_%d", seg);
  }else{
    throw Exception(FUNC_NAME + Form(" no such plane : %d", lid));
  }

  auto node = m_geometry->GetNode(node_name);
  if(!node){
    throw Exception(FUNC_NAME + Form(" no such node : %s",
                                     node_name.Data()));
  }

  node->SetVisibility(1);

  if(Tu>0 && Td>0){
    node->SetLineWidth(2);
    node->SetLineColor(kBlue);
  }
  else {
    node->SetLineWidth(2);
    node->SetLineColor(kGreen);
  }

  m_canvas->cd(1)->cd(2);
  m_geometry->Draw();
  m_canvas->Update();
}

//_____________________________________________________________________________
void
EventDisplay::DrawBcOutLocalTrack(DCLocalTrack *tp)
{
  const Double_t offsetZ = -3108.4;

  Double_t z0 = gGeom.GetLocalZ("BC3-X1") - 100.;
  Double_t x0 = tp->GetX(z0) + BeamAxis;
  Double_t y0 = tp->GetY(z0);
  z0 += offsetZ;

  Double_t z1 = zK18Target + 100.;
  Double_t x1 = tp->GetX(z1) + BeamAxis;
  Double_t y1 = tp->GetY(z1);
  z1 += offsetZ;

  ThreeVector gPos0(x0, y0, z0);
  ThreeVector gPos1(x1, y1, z1);

  auto p = new TPolyLine3D(2);
  p->SetLineColor(kRed);
  p->SetLineWidth(1);
  p->SetPoint(0, gPos0.x(), gPos0.y(), gPos0.z());
  p->SetPoint(1, gPos1.x(), gPos1.y(), gPos1.z());
  m_BcOutTrack.push_back(p);
  m_canvas->cd(1)->cd(2);
  p->Draw();

#if Vertex
  z0 = zK18Target + MinZ;
  z1 = zK18Target;
  x0 = tp->GetX(z0); y0 = tp->GetY(z0);
  x1 = tp->GetX(z1); y1 = tp->GetY(z1);
  z0 -= zK18Target;
  z1 -= zK18Target;
  {
    m_canvas_vertex->cd(1);
    TPolyLine *line = new TPolyLine(2);
    line->SetPoint(0, z0, x0);
    line->SetPoint(1, z1, x1);
    line->SetLineColor(kRed);
    line->SetLineWidth(1);
    line->Draw();
    m_BcOutXZ_line.push_back(line);
  }
  {
    m_canvas_vertex->cd(2);
    TPolyLine *line = new TPolyLine(2);
    line->SetPoint(0, z0, y0);
    line->SetPoint(1, z1, y1);
    line->SetLineColor(kRed);
    line->SetLineWidth(1);
    line->Draw();
    m_BcOutYZ_line.push_back(line);
  }
#endif
}

//_____________________________________________________________________________
// void
// EventDisplay::DrawBcOutLocalTrack(Double_t x0, Double_t y0, Double_t u0, Double_t v0)
// {
//   const Double_t offsetZ = -3108.4;

//   Double_t z0 = offsetZ;

//   Double_t z1 = zK18Target + 100.;
//   Double_t x1 = x0+u0*z1;
//   Double_t y1 = y0+v0*z1;
//   z1 += offsetZ;

//   ThreeVector gPos0(x0+BeamAxis, y0, z0);
//   ThreeVector gPos1(x1+BeamAxis, y1, z1);

//   TPolyLine3D *p = new TPolyLine3D(2);
//   p->SetLineColor(kRed);
//   p->SetLineWidth(1);
//   p->SetPoint(0, gPos0.x(), gPos0.y(), gPos0.z());
//   p->SetPoint(1, gPos1.x(), gPos1.y(), gPos1.z());
//   m_BcOutTrack.push_back(p);
//   m_canvas->cd(1)->cd(2);
//   p->Draw();
//   gPad->Update();

// #if Vertex
//   Double_t z2 = zK18Target + MinZ;
//   z1 = zK18Target;
//   Double_t x2 = x0+u0*z0, y2 = y0+v0*z0;
//   x1 = x0+u0*z1; y1 = y0+v0*z1;
//   z2 -= zK18Target;
//   z1 -= zK18Target;
//   {
//     TPolyLine *line = new TPolyLine(2);
//     line->SetPoint(0, z2, x2);
//     line->SetPoint(1, z1, x1);
//     line->SetLineColor(kRed);
//     line->SetLineWidth(1);
//     m_BcOutXZ_line.push_back(line);
//     m_canvas_vertex->cd(1);
//     line->Draw();
//   }
//   {
//     m_canvas_vertex->cd(2);
//     TPolyLine *line = new TPolyLine(2);
//     line->SetPoint(0, z2, y2);
//     line->SetPoint(1, z1, y1);
//     line->SetLineColor(kRed);
//     line->SetLineWidth(1);
//     line->Draw();
//     m_BcOutYZ_line.push_back(line);
//   }
//   m_canvas_vertex->Update();
// #endif

// }

//_____________________________________________________________________________
void
EventDisplay::DrawSdcInLocalTrack(DCLocalTrack *tp)
{
#if SdcIn
  static const Double_t zSdc1x1 = gGeom.GetLocalZ("SDC1-X1") - 100.;
  Double_t x0 = tp->GetX0(), y0 = tp->GetY0();
  Double_t x1 = tp->GetX(zSdc1x1), y1 = tp->GetY(zSdc1x1);
  // ThreeVector gPos0(x0+BeamAxis, y0, 0.);
  // ThreeVector gPos1(x1+BeamAxis, y1, zSdc1x1);
  ThreeVector gPos0(x0, y0, 0.);
  ThreeVector gPos1(x1, y1, zSdc1x1);
  auto p = new TPolyLine3D(2);
  p->SetLineColor(kRed);
  p->SetLineWidth(1);
  p->SetPoint(0, gPos0.x(), gPos0.y(), gPos0.z());
  p->SetPoint(1, gPos1.x(), gPos1.y(), gPos1.z());
  m_SdcInTrack.push_back(p);
  m_canvas->cd(1)->cd(2);
  p->Draw();
  gPad->Update();
#endif

#if Vertex
  Double_t z0 = zTarget;
  Double_t z1 = zTarget + MaxZ;
  x0 = tp->GetX(z0); y0 = tp->GetY(z0);
  x1 = tp->GetX(z1); y1 = tp->GetY(z1);
  z0 -= zTarget;
  z1 -= zTarget;
  {
    m_canvas_vertex->cd(1);
    TPolyLine *line = new TPolyLine(2);
    line->SetPoint(0, z0, x0);
    line->SetPoint(1, z1, x1);
    line->SetLineColor(kRed);
    line->SetLineWidth(1);
    line->Draw();
    m_SdcInXZ_line.push_back(line);
  }
  {
    m_canvas_vertex->cd(2);
    TPolyLine *line = new TPolyLine(2);
    line->SetPoint(0, z0, y0);
    line->SetPoint(1, z1, y1);
    line->SetLineColor(kRed);
    line->SetLineWidth(1);
    line->Draw();
    m_SdcInYZ_line.push_back(line);
  }
  m_canvas_vertex->Update();
#endif

}

//_____________________________________________________________________________
void
EventDisplay::DrawSdcOutLocalTrack(DCLocalTrack *tp)
{
#if SdcOut
  Double_t x0 = tp->GetX0(), y0 = tp->GetY0();
  // Double_t zSdcOut = gGeom.GetLocalZ("TOF-DX");
  Double_t zSdcOut = gGeom.GetLocalZ("RKINIT");
  Double_t x1 = tp->GetX(zSdcOut), y1 = tp->GetY(zSdcOut);

  ThreeVector gPos0(x0, y0, 0.);
  ThreeVector gPos1(x1, y1, zSdcOut);

  TPolyLine3D *p = new TPolyLine3D(2);
  p->SetLineColor(kRed);
  p->SetLineWidth(1);
  p->SetPoint(0, gPos0.x(), gPos0.y(), gPos0.z());
  p->SetPoint(1, gPos1.x(), gPos1.y(), gPos1.z());
  m_SdcOutTrack.push_back(p);
  m_canvas->cd(1)->cd(2);
  p->Draw();
  gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawVertex(const ThreeVector& vertex)
{
#if Vertex
  del::DeleteObject(m_VertexPointXZ);
  del::DeleteObject(m_VertexPointYZ);

  Double_t x = vertex.x();
  Double_t y = vertex.y();
  Double_t z = vertex.z();

  m_VertexPointXZ = new TMarker(z, x, 34);
  m_VertexPointXZ->SetMarkerSize(1);
  m_VertexPointXZ->SetMarkerColor(kBlue);
  m_VertexPointYZ = new TMarker(z, y, 34);
  m_VertexPointYZ->SetMarkerSize(1);
  m_VertexPointYZ->SetMarkerColor(kBlue);

  m_canvas_vertex->cd(1);
  m_VertexPointXZ->Draw();
  m_canvas_vertex->cd(2);
  m_VertexPointYZ->Draw();
  m_canvas_vertex->Update();
  hddaq::cout <<"Draw Vertex! " << x << " " << y << " " << z << std::endl;
#endif

}

//_____________________________________________________________________________
void
EventDisplay::DrawMissingMomentum(const ThreeVector& mom, const ThreeVector& pos)
{
#if Vertex
  if(std::abs(pos.x())>40.) return;
  if(std::abs(pos.y())>20.) return;
  if(std::abs(pos.z())>20.) return;

  Double_t z0 = pos.z();
  Double_t x0 = pos.x();
  Double_t y0 = pos.y();
  Double_t u0 = mom.x()/mom.z();
  Double_t v0 = mom.y()/mom.z();
  Double_t z1 = MaxZ;
  Double_t x1 = x0 + u0*z1;
  Double_t y1 = y0 + v0*z1;
  m_canvas_vertex->cd(1);
  m_MissMomXZ_line = new TPolyLine(2);
  m_MissMomXZ_line->SetPoint(0, z0, x0);
  m_MissMomXZ_line->SetPoint(1, z1, x1);
  m_MissMomXZ_line->SetLineColor(kBlue);
  m_MissMomXZ_line->SetLineWidth(1);
  m_MissMomXZ_line->Draw();
  m_canvas_vertex->cd(2);
  m_MissMomYZ_line = new TPolyLine(2);
  m_MissMomYZ_line->SetPoint(0, z0, y0);
  m_MissMomYZ_line->SetPoint(1, z1, y1);
  m_MissMomYZ_line->SetLineColor(kBlue);
  m_MissMomYZ_line->SetLineWidth(1);
  m_MissMomYZ_line->Draw();
  m_canvas_vertex->Update();
  ::sleep(3);
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawKuramaTrack(Int_t nStep, const std::vector<TVector3>& StepPoint,
                              Double_t q)
{
  del::DeleteObject(m_kurama_step_mark);

  m_kurama_step_mark = new TPolyMarker3D(nStep);
  for(Int_t i=0; i<nStep; ++i){
    m_kurama_step_mark->SetPoint(i,
                                 StepPoint[i].x(),
                                 StepPoint[i].y(),
                                 StepPoint[i].z());
  }

  Color_t color = (q > 0) ? kRed : kBlue;

  if(!m_is_save_mode){
    m_kurama_step_mark->SetMarkerSize(1);
    m_kurama_step_mark->SetMarkerStyle(6);
  }
  m_kurama_step_mark->SetMarkerColor(color);

  m_canvas->cd(1)->cd(2);
  m_kurama_step_mark->Draw();
  m_canvas->Update();

#if Vertex
  del::DeleteObject(m_KuramaMarkVertexX);
  m_KuramaMarkVertexX = new TPolyMarker(nStep);
  for(Int_t i=0; i<nStep; ++i){
    Double_t x = StepPoint[i].x()-BeamAxis;
    Double_t z = StepPoint[i].z()-zTarget;
    m_KuramaMarkVertexX->SetPoint(i, z, x);
  }
  m_KuramaMarkVertexX->SetMarkerSize(0.4);
  m_KuramaMarkVertexX->SetMarkerColor(color);
  m_KuramaMarkVertexX->SetMarkerStyle(6);
  m_canvas_vertex->cd(1);
  m_KuramaMarkVertexX->Draw();
  del::DeleteObject(m_KuramaMarkVertexY);
  m_KuramaMarkVertexY = new TPolyMarker(nStep);
  for(Int_t i=0; i<nStep; ++i){
    Double_t y = StepPoint[i].y();
    Double_t z = StepPoint[i].z()-zTarget;
    m_KuramaMarkVertexY->SetPoint(i, z, y);
  }
  m_KuramaMarkVertexY->SetMarkerSize(0.4);
  m_KuramaMarkVertexY->SetMarkerColor(color);
  m_KuramaMarkVertexY->SetMarkerStyle(6);
  m_canvas_vertex->cd(2);
  m_KuramaMarkVertexY->Draw();
  m_canvas_vertex->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawTarget()
{
  if(!m_target_node){
    hddaq::cout << FUNC_NAME << " " << "node is null" << std::endl;
    return;
  }
  m_target_node->SetLineColor(kMagenta);
  m_canvas->cd(1)->cd(2);
  m_canvas->Update();
}

//_____________________________________________________________________________
void
EventDisplay::DrawRunEvent(Double_t xpos, Double_t ypos, const TString& arg)
{
  if(arg.Contains("Run")){
    // std::cout << arg << " find " << std::endl;
    m_canvas->cd(1)->cd(1)->Clear();
  }
  m_canvas->cd(1)->cd(1);
  TLatex tex;
  tex.SetTextAlign(12);
  tex.SetTextSize(0.5);
  // tex.SetTextSize(0.15);
  tex.SetNDC();
  tex.DrawLatex(xpos, ypos, arg);
  m_canvas->Update();
}

//_____________________________________________________________________________
void
EventDisplay::DrawText(Double_t xpos, Double_t ypos, const TString& arg)
{
  m_canvas->cd(1)->cd(2);
  TLatex tex;
  tex.SetTextAlign(12);
  tex.SetTextFont(42);
  tex.SetTextSize(0.025);
  tex.SetNDC();
  tex.DrawLatex(xpos, ypos, arg);
  m_canvas->Update();
}

//_____________________________________________________________________________
void
EventDisplay::EndOfEvent()
{
  del::DeleteObject(m_init_step_mark);
  del::DeleteObject(m_BcOutXZ_line);
  del::DeleteObject(m_BcOutYZ_line);
  del::DeleteObject(m_SdcInXZ_line);
  del::DeleteObject(m_SdcInYZ_line);
  del::DeleteObject(m_BcInTrack);
  del::DeleteObject(m_BcOutTrack);
  del::DeleteObject(m_BcOutTrackShs);
  del::DeleteObject(m_BcOutTrack2);
  del::DeleteObject(m_BcOutTrack3);
  del::DeleteObject(m_SdcInTrack);
  del::DeleteObject(m_SdcInTrack2);
  del::DeleteObject(m_SdcOutTrack);
  del::DeleteObject(m_hs_step_mark);
  del::DeleteObject(m_kurama_step_mark);
  del::DeleteObject(m_VertexPointXZ);
  del::DeleteObject(m_VertexPointYZ);
  del::DeleteObject(m_MissMomXZ_line);
  del::DeleteObject(m_MissMomYZ_line);
  del::DeleteObject(m_HSMarkVertexXShs);
  del::DeleteObject(m_KuramaMarkVertexXShs);
  del::DeleteObject(m_KuramaMarkVertexX);
  del::DeleteObject(m_KuramaMarkVertexY);
  ResetVisibility();
  ResetHist();
}

//_____________________________________________________________________________
void
EventDisplay::ResetVisibility(TNode *& node, Color_t c)
{
  if(!node) return;
  node->SetLineWidth(1);
  if(c==kWhite)
    node->SetVisibility(kFALSE);
  else
    node->SetLineColor(c);
}

//_____________________________________________________________________________
void
EventDisplay::ResetVisibility(std::vector<TNode*>& node, Color_t c)
{
  const std::size_t n = node.size();
  for(std::size_t i=0; i<n; ++i){
    ResetVisibility(node[i], c);
  }
}

//_____________________________________________________________________________
void
EventDisplay::ResetVisibility()
{
  ResetVisibility(m_BC3x1_node);
  ResetVisibility(m_BC3x2_node);
  ResetVisibility(m_BC3v1_node);
  ResetVisibility(m_BC3v2_node);
  ResetVisibility(m_BC3u1_node);
  ResetVisibility(m_BC3u2_node);
  ResetVisibility(m_BC4u1_node);
  ResetVisibility(m_BC4u2_node);
  ResetVisibility(m_BC4v1_node);
  ResetVisibility(m_BC4v2_node);
  ResetVisibility(m_BC4x1_node);
  ResetVisibility(m_BC4x2_node);
  ResetVisibility(m_SDC1v1_node);
  ResetVisibility(m_SDC1v2_node);
  ResetVisibility(m_SDC1x1_node);
  ResetVisibility(m_SDC1x2_node);
  ResetVisibility(m_SDC1u1_node);
  ResetVisibility(m_SDC1u2_node);
  ResetVisibility(m_SDC2x1_node);
  ResetVisibility(m_SDC2x2_node);
  ResetVisibility(m_SDC2y1_node);
  ResetVisibility(m_SDC2y2_node);
  ResetVisibility(m_SDC3x1_node);
  ResetVisibility(m_SDC3x2_node);
  ResetVisibility(m_SDC3y1_node);
  ResetVisibility(m_SDC3y2_node);
  ResetVisibility(m_SDC4y1_node);
  ResetVisibility(m_SDC4y2_node);
  ResetVisibility(m_SDC4x1_node);
  ResetVisibility(m_SDC4x2_node);
  ResetVisibility(m_BH2seg_node, kBlack);
  ResetVisibility(m_SCHseg_node, kBlack);
  ResetVisibility(m_TOFseg_node, kBlack);
  ResetVisibility(m_WCseg_node, kBlack);
  ResetVisibility(m_target_node, kBlack);
}

//_____________________________________________________________________________
void
EventDisplay::ResetHist()
{
#if Hist_Timing
  m_hist_bh2->Reset();
  m_hist_sch->Reset();
  m_hist_tof->Reset();
  m_hist_sdc1->Reset();
  m_hist_sdc1p->Reset();

  m_hist_bc3->Reset();
  m_hist_bc3p->Reset();
  m_hist_bc3u->Reset();
  m_hist_bc3up->Reset();
  m_hist_bc3v->Reset();
  m_hist_bc3vp->Reset();

  m_hist_bc4->Reset();
  m_hist_bc4p->Reset();
  m_hist_bc4u->Reset();
  m_hist_bc4up->Reset();
  m_hist_bc4v->Reset();
  m_hist_bc4vp->Reset();

  m_hist_bc3_time->Reset();
  m_hist_bc3_time->SetMaximum(-1111);
  m_hist_bc3p_time->Reset();

  m_hist_bc4_time->Reset();
  m_hist_bc4_time->SetMaximum(-1111);
  m_hist_bc4p_time->Reset();
#endif

#if Hist_SdcOut
  m_hist_sdc3_l->Reset();
  m_hist_sdc3_t->Reset();
  m_hist_sdc3p_l->Reset();
  m_hist_sdc3p_t->Reset();

  m_hist_sdc3y_l->Reset();
  m_hist_sdc3y_t->Reset();
  m_hist_sdc3yp_l->Reset();
  m_hist_sdc3yp_t->Reset();

  m_hist_sdc4_l->Reset();
  m_hist_sdc4_t->Reset();
  m_hist_sdc4p_l->Reset();
  m_hist_sdc4p_t->Reset();

  m_hist_sdc4y_l->Reset();
  m_hist_sdc4y_t->Reset();
  m_hist_sdc4yp_l->Reset();
  m_hist_sdc4yp_t->Reset();

#endif

#if Hist_Timing
  m_hist_bh1->Reset();
  m_hist_bft->Reset();
  m_hist_bft_p->Reset();
  m_hist_bcIn->Reset();
  m_hist_bcOut->Reset();

  Int_t nc=m_BH1box_cont.size();
  for (Int_t i=0; i<nc; i++)
    m_BH1box_cont[i]->SetFillColor(kWhite);

  Int_t ncBh2=m_BH2box_cont.size();
  for (Int_t i=0; i<ncBh2; i++) {
    m_BH2box_cont[i]->SetFillColor(kWhite);
    m_BH2box_cont2[i]->SetFillColor(kWhite);
  }

  Int_t ncSch=m_SCHbox_cont.size();
  for (Int_t i=0; i<ncSch; i++)
    m_SCHbox_cont[i]->SetFillColor(kWhite);

  m_hist_bcOut_sdcIn->Reset();
  m_hist_sdcIn_predict->Reset();
  m_hist_sdcIn_predict2->Reset();
#endif

}

//_____________________________________________________________________________
void
EventDisplay::FillMomentum(Double_t momentum)
{
#if Hist
  m_hist_p->Fill(momentum);
  m_canvas_hist->cd(1);
  gPad->Modified();
  gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillMassSquare(Double_t mass_square)
{
#if Hist
  m_hist_m2->Fill(mass_square);
  m_canvas_hist->cd(2);
  gPad->Modified();
  gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillMissMass(Double_t missmass)
{
#if Hist
  m_hist_missmass->Fill(missmass);
  m_canvas_hist->cd(3);
  gPad->Modified();
  gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillBH1(Int_t seg, Int_t tdc)
{
#if Hist_Timing
  Double_t p0 = gHodo.GetOffset(DetIdBH1, 0, seg, 0);
  Double_t p1 = gHodo.GetGain(DetIdBH1, 0, seg, 0);
  m_hist_bh1->Fill(seg, p1*((Double_t)tdc-p0));
#endif
}

//_____________________________________________________________________________
void
EventDisplay::SetCorrectTimeBH1(Int_t seg, Double_t de)
{
#if Hist_BcIn
  Int_t color;
  if (de<0.5)
    color = kBlue;
  else if (de >= 0.5 && de <= 1.5)
    color = kOrange;
  else
    color = kRed;
  m_BH1box_cont[seg]->SetFillColor(color);
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillBFT(Int_t layer, Int_t seg, Int_t tdc)
{
#if Hist_Timing
  Double_t p0 = gHodo.GetOffset(DetIdBFT, layer, seg, 0);
  Double_t p1 = gHodo.GetGain(DetIdBFT, layer, seg, 0);

  TH2 *hp=0;

  if (layer == 0)
    hp = m_hist_bft;
  else
    hp = m_hist_bft_p;

  hp->Fill(seg, p1*((Double_t)tdc-p0));
  //m_canvas_hist2->cd(1);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawBcInTrack(Double_t x0, Double_t u0)
{
#if Hist_Timing

  Double_t z1 = -150, z2 = 50;
  TLine *l = new TLine(x0+u0*(z1-zBFT), z1, x0+u0*(z2-zBFT), z2);
  m_BcInTrack.push_back(l);
#endif
}

//_____________________________________________________________________________
void
EventDisplay::SetCorrectTimeBFT(Double_t pos)
{
#if Hist_BcIn
  m_hist_bcIn->Fill(pos, zBFT);
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillBH2(Int_t seg, Int_t tdc)
{
#if Hist_Timing
  Double_t p0 = gHodo.GetOffset(DetIdBH2, 0, seg, 0);
  Double_t p1 = gHodo.GetGain(DetIdBH2, 0, seg, 0);
  m_hist_bh2->Fill(seg, p1*((Double_t)tdc-p0));
  //m_canvas_hist2->cd(1);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::SetCorrectTimeBH2(Int_t seg, Double_t de)
{
#if Hist_BcIn
  Int_t color;
  if (de<0.5)
    color = kBlue;
  else if (de >= 0.5 && de <= 1.5)
    color = kOrange;
  else
    color = kRed;
  m_BH2box_cont[seg]->SetFillColor(color);
  m_BH2box_cont2[seg]->SetFillColor(color);
#endif
}

//_____________________________________________________________________________
void
EventDisplay::SetCorrectTimeBcOut(Int_t layer, Double_t pos)
{
#if Hist_BcIn
  Double_t z = gGeom.GetLocalZ(layer);
  m_hist_bcOut->Fill(pos, z);

  m_hist_bcOut_sdcIn->Fill(pos + gGeom.GetGlobalPosition(layer).x(),
                           gGeom.GetGlobalPosition(layer).z());
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawBcOutTrack(Double_t x0, Double_t u0, Double_t y0, Double_t v0, Bool_t flagGoodForTracking )
{
#if Hist_BcIn

  Double_t z1 = 0, z2 = 600;
  TLine *l = new TLine(x0+u0*z1, z1, x0+u0*z2, z2);
  if (! flagGoodForTracking)
    l->SetLineColor(kOrange);

  m_BcOutTrack2.push_back(l);

  Double_t z3 = -50, z4 = 2950;
  TLine *l2 = new TLine(x0+u0*z3 + gxK18Target, z3 - (zK18Target-gzK18Target),
                        x0+u0*z4 + gxK18Target, z4 - (zK18Target-gzK18Target));
  if (! flagGoodForTracking)
    l2->SetLineColor(kOrange);

  m_BcOutTrack3.push_back(l2);

  TLine *l3 = new TLine(y0+v0*z3, z3 - (zK18Target-gzK18Target),
                        y0+v0*z4, z4 - (zK18Target-gzK18Target));
  if (! flagGoodForTracking)
    l3->SetLineColor(kOrange);

  m_BcOutTrack3.push_back(l3);

  for (Int_t layer=1; layer<=9; layer++) {
    Double_t z_bcOut = gGeom.GetLocalZ(layer) - gzK18Target + zK18Target;
    Double_t x = x0 + u0*z_bcOut + gxK18Target;
    Double_t y = y0 + v0*z_bcOut;
    Double_t z = gGeom.GetGlobalPosition(layer).z();

    ThreeVector gloPos(x, y, z);
    ThreeVector localPos = gGeom.Global2LocalPos(layer, gloPos);
    m_hist_sdcIn_predict->Fill(localPos.x(), z);

  }



#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawSdcInTrack(Double_t x0, Double_t u0, Double_t y0, Double_t v0, Bool_t flagKurama, Bool_t flagBeam)
{
#if Hist_BcIn

  Double_t z1 = -1500, z2 = 0;
  TLine *l = new TLine(x0+u0*z1, z1, x0+u0*z2, z2);
  if (flagKurama)
    l->SetLineColor(kRed);
  if (flagBeam)
    l->SetLineColor(kOrange);
  m_SdcInTrack2.push_back(l);

  TLine *l2 = new TLine(y0+v0*z1, z1, y0+v0*z2, z2);
  if (flagKurama)
    l2->SetLineColor(kRed);
  if (flagBeam)
    l2->SetLineColor(kOrange);

  m_SdcInTrack2.push_back(l2);


  for (Int_t layer=1; layer<=9; layer++) {
    Double_t z = gGeom.GetLocalZ(layer);
    Double_t x = x0 + u0*z;
    Double_t y = y0 + v0*z;
    Double_t gz = gGeom.GetGlobalPosition(layer).z();

    ThreeVector gloPos(x, y, gz);
    ThreeVector localPos = gGeom.Global2LocalPos(layer, gloPos);
    m_hist_sdcIn_predict2->Fill(localPos.x(), z);

  }

#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillSCH(Int_t seg, Int_t tdc)
{
#if Hist_Timing
  Double_t p0 = gHodo.GetOffset(DetIdSCH, 0, seg, 0);
  Double_t p1 = gHodo.GetGain(DetIdSCH, 0, seg, 0);

  m_hist_sch->Fill(seg, p1*((Double_t)tdc-p0));
  //m_canvas_hist2->cd(2);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::SetCorrectTimeSCH(Int_t seg, Double_t de)
{
#if Hist_BcIn
  Int_t color;
  if (de<30)
    color = kBlue;
  else if (de >= 30 && de <= 60)
    color = kOrange;
  else
    color = kRed;
  m_SCHbox_cont[seg]->SetFillColor(color);
#endif
}

//_____________________________________________________________________________
void
EventDisplay::SetCorrectTimeSdcIn(Int_t layer, Double_t pos)
{
#if Hist_BcIn
  Double_t z = gGeom.GetLocalZ(layer);
  m_hist_bcOut_sdcIn->Fill(pos, z);

#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillTOF(Int_t seg, Int_t tdc)
{
#if Hist_Timing
  Double_t p0 = gHodo.GetOffset(DetIdTOF, 0, seg, 0);
  Double_t p1 = gHodo.GetGain(DetIdTOF, 0, seg, 0);

  m_hist_tof->Fill(seg, p1*((Double_t)tdc-p0));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillBcOutHit(Int_t layer,  Int_t wire, Int_t tdc)
{
#if Hist_Timing
  TH2 *hp=0;

  if (layer == 1)
    hp = m_hist_bc3;
  else if (layer == 2)
    hp = m_hist_bc3p;
  else if (layer == 3)
    hp = m_hist_bc3v;
  else if (layer == 4)
    hp = m_hist_bc3vp;
  else if (layer == 5)
    hp = m_hist_bc3u;
  else if (layer == 6)
    hp = m_hist_bc3up;
  else if (layer == 7)
    hp = m_hist_bc4u;
  else if (layer == 8)
    hp = m_hist_bc4up;
  else if (layer == 9)
    hp = m_hist_bc4v;
  else if (layer == 10)
    hp = m_hist_bc4vp;
  else if (layer == 11)
    hp = m_hist_bc4;
  else if (layer == 12)
    hp = m_hist_bc4p;

  Double_t p0=0.0, p1=-0.;
  gTdc.GetParameter(layer + 112, wire, p0, p1);

  hp->Fill(wire, -p1*(tdc+p0));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawBC3(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_bc3->Fill(wire, -1.*(tdc-351));
  m_hist_bc3_time->Fill(-1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawBC3p(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_bc3p->Fill(wire, -1.*(tdc-351));
  m_hist_bc3p_time->Fill(-1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//_____________________________________________________________________________
void
EventDisplay::DrawBC4(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_bc4->Fill(wire, -1.*(tdc-351));
  m_hist_bc4_time->Fill(-1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawBC4p(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_bc4p->Fill(wire, -1.*(tdc-351));
  m_hist_bc4p_time->Fill(-1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillSDC1(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc1->Fill(wire, -1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//_____________________________________________________________________________
void
EventDisplay::FillSDC1p(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc1p->Fill(wire, -1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//_____________________________________________________________________________
void
EventDisplay::FillSdcOutHit(Int_t layer,  Int_t wire, Int_t LorT, Int_t tdc)
{
#if Hist_SdcOut
  TH2 *hp=0;

  if (layer == 1 && LorT == 0)
    hp = m_hist_sdc3_l;
  else if (layer == 1 && LorT == 1)
    hp = m_hist_sdc3_t;
  else if (layer == 2 && LorT == 0)
    hp = m_hist_sdc3p_l;
  else if (layer == 2 && LorT == 1)
    hp = m_hist_sdc3p_t;
  else if (layer == 3 && LorT == 0)
    hp = m_hist_sdc3y_l;
  else if (layer == 3 && LorT == 1)
    hp = m_hist_sdc3y_t;
  else if (layer == 4 && LorT == 0)
    hp = m_hist_sdc3yp_l;
  else if (layer == 4 && LorT == 1)
    hp = m_hist_sdc3yp_t;
  else if (layer == 5 && LorT == 0)
    hp = m_hist_sdc4y_l;
  else if (layer == 5 && LorT == 1)
    hp = m_hist_sdc4y_t;
  else if (layer == 6 && LorT == 0)
    hp = m_hist_sdc4yp_l;
  else if (layer == 6 && LorT == 1)
    hp = m_hist_sdc4yp_t;
  else if (layer == 7 && LorT == 0)
    hp = m_hist_sdc4_l;
  else if (layer == 7 && LorT == 1)
    hp = m_hist_sdc4_t;
  else if (layer == 8 && LorT == 0)
    hp = m_hist_sdc4p_l;
  else if (layer == 8 && LorT == 1)
    hp = m_hist_sdc4p_t;

  Double_t p0=0.0, p1=-0.;
  gTdc.GetParameter(layer + 30, wire, p0, p1);

  hp->Fill(wire, -p1*(tdc+p0));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//_____________________________________________________________________________
void
EventDisplay::FillSDC3_Leading(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc3_l->Fill(wire, -0.833*(tdc-890));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//_____________________________________________________________________________
void
EventDisplay::FillSDC3_Trailing(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc3_t->Fill(wire, -0.833*(tdc-890));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillSDC3p_Leading(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc3p_l->Fill(wire, -0.833*(tdc-890));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillSDC3p_Trailing(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc3p_t->Fill(wire, -0.833*(tdc-890));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//_____________________________________________________________________________
void
EventDisplay::FillSDC4_Leading(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc4_l->Fill(wire, -0.833*(tdc-885));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillSDC4_Trailing(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc4_t->Fill(wire, -0.833*(tdc-885));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillSDC4p_Leading(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc4p_l->Fill(wire, -0.833*(tdc-885));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillSDC4p_Trailing(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc4p_t->Fill(wire, -0.833*(tdc-885));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::Update()
{
  TIter canvas_iterator(gROOT->GetListOfCanvases());
  while (true) {
    auto canvas = dynamic_cast<TCanvas*>(canvas_iterator.Next());
    if (!canvas) break;
    canvas->UseCurrentStyle();
    canvas->cd(1)->SetLogz();
    canvas->Modified();
    canvas->Update();
  }

// #if Hist_Timing
//   for (Int_t i=0; i<9; i++) {
//     m_canvas_hist2->cd(i+1);
//     gPad->Modified();
//     gPad->Update();
//   }

//   Int_t max1 = m_hist_bc3_time->GetMaximum();
//   Int_t max2 = m_hist_bc3p_time->GetMaximum();
//   if (max2 > max1)
//     m_hist_bc3_time->SetMaximum(max2*1.1);

//   max1 = m_hist_bc4_time->GetMaximum();
//   max2 = m_hist_bc4p_time->GetMaximum();
//   if (max2 > max1)
//     m_hist_bc4_time->SetMaximum(max2*1.1);

//   for (Int_t i=0; i<6; i++) {
//     m_canvas_hist3->cd(i+1);
//     gPad->Modified();
//     gPad->Update();
//   }
// #endif

// #if Hist_SdcOut
//   for (Int_t i=0; i<8; i++) {
//     m_canvas_hist4->cd(i+1);
//     gPad->Modified();
//     gPad->Update();
//   }
// #endif

// #if Hist_BcIn
//   for (Int_t i=0; i<4; i++) {
//     m_canvas_hist5->cd(i+1);


//     if (i==2) {
//       Int_t nc=m_BcInTrack.size();
//       for (Int_t i=0; i<nc; i++)
//  m_BcInTrack[i]->Draw("same");
//     } else if (i==3) {
//       Int_t nc=m_BcOutTrack2.size();
//       for (Int_t i=0; i<nc; i++)
//  m_BcOutTrack2[i]->Draw("same");
//     }


//     gPad->Modified();
//     gPad->Update();
//   }

//   {
//     m_canvas_hist6->cd();

//     Int_t nc=m_BcOutTrack3.size();
//     for (Int_t i=0; i<nc; i++)
//       m_BcOutTrack3[i]->Draw("same");

//     Int_t ncSdcIn=m_SdcInTrack2.size();
//     for (Int_t i=0; i<ncSdcIn; i++)
//  m_SdcInTrack2[i]->Draw("same");

//     gPad->Modified();
//     gPad->Update();
//   }

// #endif

}

//_____________________________________________________________________________
void
EventDisplay::CalcRotMatrix(Double_t TA, Double_t RA1, Double_t RA2, Double_t *rotMat)
{
  Double_t ct0 = TMath::Cos(TA*TMath::DegToRad() );
  Double_t st0 = TMath::Sin(TA*TMath::DegToRad() );
  Double_t ct1 = TMath::Cos(RA1*TMath::DegToRad());
  Double_t st1 = TMath::Sin(RA1*TMath::DegToRad());
  Double_t ct2 = TMath::Cos(RA2*TMath::DegToRad());
  Double_t st2 = TMath::Sin(RA2*TMath::DegToRad());

  Double_t rotMat1[3][3], rotMat2[3][3];

  /* rotation matrix which is same as DCGeomRecord.cc */

#if 1
  /* new definition of RA2 */
  rotMat1[0][0] =  ct0*ct2+st0*st1*st2;
  rotMat1[0][1] = -st0*ct2+ct0*st1*st2;
  rotMat1[0][2] =  ct1*st2;

  rotMat1[1][0] =  st0*ct1;
  rotMat1[1][1] =  ct0*ct1;
  rotMat1[1][2] = -st1;

  rotMat1[2][0] = -ct0*st2+st0*st1*ct2;
  rotMat1[2][1] =  st0*st2+ct0*st1*ct2;
  rotMat1[2][2] =  ct1*ct2;

#else
  /* old definition of RA2 */
  rotMat1[0][0] =  ct0*ct2-st0*ct1*st2;
  rotMat1[0][1] = -st0*ct2-ct0*ct1*st2;
  rotMat1[0][2] =  st1*st2;

  rotMat1[1][0] =  ct0*st2+st0*ct1*ct2;
  rotMat1[1][1] = -st0*st2+ct0*ct1*ct2;
  rotMat1[1][2] = -st1*ct2;

  rotMat1[2][0] =  st0*st1;
  rotMat1[2][1] =  ct0*st1;
  rotMat1[2][2] =  ct1;
#endif

  /* rotation matrix which rotate -90 deg at x axis*/
  rotMat2[0][0] =  1.0;
  rotMat2[0][1] =  0.0;
  rotMat2[0][2] =  0.0;

  rotMat2[1][0] =  0.0;
  rotMat2[1][1] =  0.0;
  rotMat2[1][2] =  1.0;

  rotMat2[2][0] =  0.0;
  rotMat2[2][1] = -1.0;
  rotMat2[2][2] =  0.0;

  for (Int_t i=0; i<9; i++)
    rotMat[i]=0.0;

  for (Int_t i=0; i<3; i++) {
    for (Int_t j=0; j<3; j++) {
      for (Int_t k=0; k<3; k++) {
        //rotMat[3*i+j] += rotMat1[i][k]*rotMat2[k][j];
        rotMat[i+3*j] += rotMat1[i][k]*rotMat2[k][j];
      }
    }
  }

}

//_____________________________________________________________________________
Int_t
EventDisplay::GetCommand()
{
  Update();
  char ch;
  char data[100];
  static Int_t stat   = 0;
  static Int_t Nevent = 0;
  static Int_t ev     = 0;

  const Int_t kSkip = 1;
  const Int_t kNormal = 0;

  if (stat == 1 && Nevent > 0 && ev<Nevent) {
    ev++;
    return kSkip;
  }
  if (ev==Nevent) {
    stat=0;
    ev=0;
  }

  if (stat == 0) {
    hddaq::cout << "q|n|p>" << std::endl;

    scanf("%c",&ch);
    if (ch!='\n')
      while(getchar() != '\n');

    switch (ch) {
    case 'q':
      std::exit(EXIT_SUCCESS);
    case 'n':
      stat = 1;
      do {
        printf("event#>");
        scanf("%s",data);
      } while ((Nevent=atoi(data))<=0);
      //hddaq::cout << "Continue " << Nevent << "event" << std::endl;
      hddaq::cout << "Skip " << Nevent << "event" << std::endl;
      break;
    case 'p':
      m_theApp->Run(kTRUE);
      break;
    }
  }

  if (stat == 1)
    return kSkip;

  return kNormal;
}

//_____________________________________________________________________________
void
EventDisplay::Print(Int_t run_number, Int_t event_number)
{
  static Int_t prev_run_number = 0;
  static TString fig_dir;
  if(run_number != prev_run_number){
    fig_dir = Form("fig/evdisp/run%05d", run_number);
    gSystem->MakeDirectory(fig_dir);
  }
  m_canvas->Print(Form("%s/evdisp_run%05d_ev%d.png",
                       fig_dir.Data(), run_number, event_number));
  prev_run_number = run_number;
}

//_____________________________________________________________________________
void
EventDisplay::Run(Bool_t flag)
{
  hddaq::cout << FUNC_NAME << " TApplication is running" << std::endl;

  m_theApp->Run(flag);
}


//_____________________________________________________________________________
void
EventDisplay::DrawHSTrack(Int_t nStep, const std::vector<TVector3>& StepPoint,
                              Double_t q)
{
  del::DeleteObject(m_hs_step_mark);

  m_hs_step_mark = new TPolyMarker3D(nStep);
  for(Int_t i=0; i<nStep; ++i){
    m_hs_step_mark->SetPoint(i,
			     StepPoint[i].x(),
			     StepPoint[i].y(),
			     StepPoint[i].z());
  }

  Color_t color = (q > 0) ? kRed : kBlue;

  if(!m_is_save_mode){
    m_hs_step_mark->SetMarkerSize(1);
    m_hs_step_mark->SetMarkerStyle(6);
  }
  m_hs_step_mark->SetMarkerColor(color);

  m_canvas->cd(1)->cd(2);
  m_hs_step_mark->Draw();
  m_canvas->Update();

}

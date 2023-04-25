// -*- C++ -*-

#include "EventDisplayShs.hh"

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
#include <TF2.h>
#include <TGeometry.h>
#include <TGraph.h>
#include <TGraph2D.h>
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
#include <TPave.h>
#include <TRandom.h>
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
#include <TSystem.h>
#include <TPGON.h>
#include <TTRD1.h>
#include <TTRD2.h>
#include <TTUBS.h>
#include <TTUBS.h>
#include <TView.h>
#include <TLine.h>

#include <std_ostream.hh>

#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "DeleteUtility.hh"
#include "HodoParamMan.hh"
#include "DCTdcCalibMan.hh"
#include "TPCPadHelper.hh"

namespace
{
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gTdc = DCTdcCalibMan::GetInstance();
const auto& IdTarget = gGeom.DetectorId("Target");
const auto& zTarget = gGeom.LocalZ("Target");
const auto& zK18Target = gGeom.LocalZ("K18Target");
const auto& gzK18Target = gGeom.GlobalZ("K18Target");
//const Double_t BeamAxis = -150.; //E07
const Double_t BeamAxis = -240.; //E40

//_____________________________________________________________________________
// local function
[[maybe_unused]] void
ConstructionDone(const TString& name, std::ostream& ost=hddaq::cout)
{
  const Int_t n = 20;
  Int_t s = name.Length();
  ost << " " << name << " ";
  for (Int_t i=0; i<n-s; ++i) ost << ".";
  ost << " done" << std::endl;
}
}

//_____________________________________________________________________________
EventDisplayShs::EventDisplayShs( void )
  : m_is_ready(false),
    m_theApp(nullptr),
    m_geometry(nullptr),
    m_node(nullptr),
    m_canvas(nullptr),
    m_tpc_adc(nullptr),
    m_tpc_tdc(nullptr),
    m_tpc_adc2d(nullptr),
    m_tpc_tdc2d(nullptr),
    m_htof_2d(nullptr),
    m_tpc_pgon(nullptr)
{
}

//_____________________________________________________________________________
EventDisplayShs::~EventDisplayShs( void )
{
}

//_____________________________________________________________________________
void
EventDisplayShs::EndOfEvent( void )
{
  // m_tpc_mark3d->SetPolyMarker(0, (Double_t*)nullptr, 0);
  // if( m_tpc_mark3d ) delete m_tpc_mark3d;
}

//_____________________________________________________________________________
void
EventDisplayShs::FillTPCADC(Int_t layer, Int_t row, Double_t adc)
{
  Int_t pad = tpc::GetPadId(layer, row);
  m_tpc_adc->Fill(adc);
  m_tpc_adc2d->SetBinContent(pad + 1, adc);
}

//_____________________________________________________________________________
void
EventDisplayShs::FillTPCTDC(Int_t layer, Int_t row, Double_t tdc)
{
  Int_t pad = tpc::GetPadId(layer, row);
  m_tpc_tdc->Fill(tdc);
  m_tpc_tdc2d->SetBinContent(pad + 1, tdc);
}

//_____________________________________________________________________________
void
EventDisplayShs::FillHTOF(Int_t segment)
{
  m_htof_2d->SetBinContent(segment, 1);
}

//_____________________________________________________________________________
Bool_t
EventDisplayShs::Initialize( void )
{
  if (m_is_ready) {
    hddaq::cerr << FUNC_NAME << " already initialied" << std::endl;
    return false;
  }

  m_theApp = new TApplication("App", 0, 0);

#if ROOT_VERSION_CODE > ROOT_VERSION(6,4,0)
  gStyle->SetPalette(kBird);
#endif
  gStyle->SetNumberContours(255);

  m_canvas = new TCanvas("Event Display", "Event Display", 1800, 900);

  m_tpc_adc = new TH1D("h_tpc_adc", "TPC ADC", 4096, 0, 4096);
  m_tpc_tdc = new TH1D("h_tpc_adc", "TPC TDC",
                       NumOfTimeBucket, 0, NumOfTimeBucket);
  const Double_t MinX = -400.;
  const Double_t MaxX =  400.;
  const Double_t MinZ = -400.;
  const Double_t MaxZ =  400.;
  m_tpc_adc2d = new TH2Poly("h_tpc_adc2d", "TPC ADC;Z;X", MinZ, MaxZ, MinX, MaxX);
  m_tpc_tdc2d = new TH2Poly("h_tpc_tdc2d", "TPC TDC;Z;X", MinZ, MaxZ, MinX, MaxX);
  m_htof_2d = new TH2Poly("h_htof_2d", "HTOF;Z;X", MinZ, MaxZ, MinX, MaxX);

  Double_t X[5];
  Double_t Y[5];
  for (Int_t l=0; l<NumOfLayersTPC; ++l) {
    Double_t pLength = tpc::padParameter[l][5];
    Double_t st      = (180.-(360./tpc::padParameter[l][3]) *
                        tpc::padParameter[l][1]/2.);
    Double_t sTheta  = (-1+st/180.)*TMath::Pi();
    Double_t dTheta  = (360./tpc::padParameter[l][3])/180.*TMath::Pi();
    Double_t cRad    = tpc::padParameter[l][2];
    Int_t    nPad    = tpc::padParameter[l][1];
    for (Int_t j=0; j<nPad; ++j) {
      X[1] = (cRad+(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
      X[2] = (cRad+(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
      X[3] = (cRad-(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
      X[4] = (cRad-(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
      X[0] = X[4];
      Y[1] = (cRad+(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
      Y[2] = (cRad+(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
      Y[3] = (cRad-(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
      Y[4] = (cRad-(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
      Y[0] = Y[4];
      for (Int_t k=0; k<5; ++k) X[k] += tpc::ZTarget;
      m_tpc_adc2d->AddBin(5, X, Y);
      m_tpc_tdc2d->AddBin(5, X, Y);
    }
  }
  {
    const Double_t L = 337;
    const Double_t t = 10;
    const Double_t w = 68;
    Double_t theta[8];
    Double_t X[5];
    Double_t Y[5];
    Double_t seg_X[5];
    Double_t seg_Y[5];
    for( Int_t i=0; i<8; i++ ){
      theta[i] = (-180+45*i)*acos(-1)/180.;
      for( Int_t j=0; j<4; j++ ){
	seg_X[1] = L-t/2.;
	seg_X[2] = L+t/2.;
	seg_X[3] = L+t/2.;
	seg_X[4] = L-t/2.;
	seg_X[0] = seg_X[4];
	seg_Y[1] = w*j-2*w;
	seg_Y[2] = w*j-2*w;
	seg_Y[3] = w*j-1*w;
	seg_Y[4] = w*j-1*w;
	seg_Y[0] = seg_Y[4];
	for( Int_t k=0; k<5; k++ ){
	  X[k] = cos(theta[i])*seg_X[k]-sin(theta[i])*seg_Y[k];
	  Y[k] = sin(theta[i])*seg_X[k]+cos(theta[i])*seg_Y[k];
	}
	m_htof_2d->AddBin(5, X, Y);
      }
    }
  }
  m_tpc_adc2d->SetStats(0);
  m_tpc_tdc2d->SetStats(0);
  m_htof_2d->SetStats(0);
  m_tpc_adc2d->SetMaximum(4000);
  m_tpc_tdc2d->SetMaximum(NumOfTimeBucket);

  m_canvas->Divide(2, 1);
  m_canvas->cd(1)->SetLogz();
  m_tpc_adc2d->Draw("colz");
  m_htof_2d->Draw("same col");
  m_canvas->cd(2);
  m_tpc_tdc2d->Draw("colz");
  m_htof_2d->Draw("same col");

  // m_canvas->cd(3);
  // m_tpc_adc->Draw("colz");
  // m_canvas->cd(4);
  // m_tpc_tdc->Draw("colz");

  // m_tpc_pgon = new TPGON("TPC_PGON", "TPC_PGON", "void",
  //                        22.5, 360, 8, 2);
  // const Double_t r = 250./TMath::Cos(TMath::DegToRad()*360/8/2);
  // m_tpc_pgon->DefineSection(0, -300, r, r);
  // m_tpc_pgon->DefineSection(1,  300, r, r);
  // m_tpc_pgon->SetFillColorAlpha(kWhite, 0);
  // m_tpc_pgon->SetLineColor(kGray);
  // m_tpc_node = new TNode("NODE1", "NODE1", "TPC_PGON", 0, 500, 500);
  // m_tpc_node->cd();
  // m_tpc_node->Draw("gl");

  // m_tpc_node->cd();
  // m_tpc_mark3d = new TPolyMarker3D;
  // m_tpc_mark3d->SetMarkerSize(10);
  // m_tpc_mark3d->SetMarkerColor(kRed+1);
  // m_tpc_mark3d->SetMarkerStyle(6);
  // m_tpc_mark3d->Draw();
  // gPad->GetView()->SetAutoRange();
  // gPad->GetView()->SetRange(-250,-250,-350,250,250,350);
  // Int_t irep = 0;
  // gPad->GetView()->SetView(0, 0, 0, irep);
  // gPad->GetView()->ShowAxis();

  m_canvas->Update();

  m_is_ready = true;
  return m_is_ready;
}

//_____________________________________________________________________________
void
EventDisplayShs::Reset( void )
{
  m_tpc_adc2d->Reset("");
  m_tpc_tdc2d->Reset("");
  m_htof_2d->Reset("");
}

//_____________________________________________________________________________
Int_t
EventDisplayShs::GetCommand( void ) const
{
  char ch;
  char data[100];
  static int stat   = 0;
  static int Nevent = 0;
  static int ev     = 0;

  const int kSkip = 1;
  const int kNormal = 0;

  if (stat == 1 && Nevent > 0 && ev<Nevent) {
    ev++;
    return kSkip;
  }
  if (ev==Nevent) {
    stat=0;
    ev=0;
  }

  if (stat == 0) {
    hddaq::cout << std::endl
                << "q: Quit." << std::endl
                << "n: Process Nevent." << std::endl
                << "p: Run TApplication." << std::endl
                << std::endl
                << "q|n|p> " << std::flush;
    std::scanf("%c",&ch);
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
      } while ( (Nevent=atoi(data))<=0 );
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
EventDisplayShs::Run(Bool_t flag)
{
  hddaq::cout << FUNC_NAME << " TApplication is running" << std::endl;
  m_theApp->Run(flag);
}

//_____________________________________________________________________________
void
EventDisplayShs::SetTPCMarker(Double_t x, Double_t y, Double_t z)
{
  m_tpc_mark3d->SetNextPoint(x, z, y);
}


//_____________________________________________________________________________
void
EventDisplayShs::SetTPCMarker(const TVector3& pos)
{
  SetTPCMarker(pos.X(), pos.Y(), pos.Z());
}

//_____________________________________________________________________________
void
EventDisplayShs::Update(void)
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
}

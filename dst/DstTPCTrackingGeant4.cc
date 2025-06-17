// -*- C++ -*-

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "Kinematics.hh"
#include "DCGeomMan.hh"
#include "DatabasePDG.hh"
#include "DetectorID.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "HodoPHCMan.hh"
#include "DCAnalyzer.hh"
#include "DCHit.hh"
#include "TPCCluster.hh"
#include "TPCLocalTrack.hh"
#include "TPCHit.hh"

#include "DstHelper.hh"
#include <TRandom.h>
#include "DebugCounter.hh"

namespace
{
using namespace root;
using namespace dst;
const std::string& class_name("DstTPCTrackingGeant4");
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();
debug::ObjectCounter& gCounter  = debug::ObjectCounter::GetInstance();

const Int_t MaxTPCHits = 10000;
const Int_t MaxTPCTracks = 100;
const Int_t MaxTPCnHits = 50;
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kTPCGeant, kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[TPCGeant]",
  "[OutFile]" };
std::vector<TString> TreeName =
{ "", "", "TPC_g", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________
struct Event
{

  Int_t evnum;
  Int_t status;
  Int_t nhittpc;                 // Number of Hits
  Int_t ntTpc;                   // Number of Tracks
  Int_t nhit_track[MaxTPCTracks]; // Number of Hits (in 1 tracks)
  Double_t chisqr[MaxTPCTracks];
  Double_t x0[MaxTPCTracks];
  Double_t y0[MaxTPCTracks];
  Double_t u0[MaxTPCTracks];
  Double_t v0[MaxTPCTracks];
  Double_t theta[MaxTPCTracks];
  Int_t hitlayer[MaxTPCTracks][MaxTPCnHits];
  Double_t hitpos_x[MaxTPCTracks][MaxTPCnHits];
  Double_t hitpos_y[MaxTPCTracks][MaxTPCnHits];
  Double_t hitpos_z[MaxTPCTracks][MaxTPCnHits];
  Double_t calpos_x[MaxTPCTracks][MaxTPCnHits];
  Double_t calpos_y[MaxTPCTracks][MaxTPCnHits];
  Double_t calpos_z[MaxTPCTracks][MaxTPCnHits];
  Double_t residual[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_x[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_y[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_z[MaxTPCTracks][MaxTPCnHits];
};

//_____________________________________________________________________
struct Src
{
  Int_t evnum;
  Int_t nhittpc;                 // Number of Hits

  Int_t nhPrm;
  Double_t xPrm[MaxTPCTracks];
  Double_t yPrm[MaxTPCTracks];
  Double_t zPrm[MaxTPCTracks];
  Double_t pxPrm[MaxTPCTracks];
  Double_t pyPrm[MaxTPCTracks];
  Double_t pzPrm[MaxTPCTracks];

	Int_t NumberOfTracks;
	Int_t PIDOfTrack[MaxTPCTracks];
	Double_t VertexOfTrack_x[MaxTPCTracks];
	Double_t VertexOfTrack_y[MaxTPCTracks];
	Double_t VertexOfTrack_z[MaxTPCTracks];
	Double_t MomentumOfTrack_x[MaxTPCTracks];
	Double_t MomentumOfTrack_y[MaxTPCTracks];
	Double_t MomentumOfTrack_z[MaxTPCTracks];


  Int_t ititpc[MaxTPCHits];
  Int_t idtpc[MaxTPCHits];
  Double_t xtpc[MaxTPCHits];//with resolution
  Double_t ytpc[MaxTPCHits];//with resolution
  Double_t ztpc[MaxTPCHits];//with resolution
  Double_t x0tpc[MaxTPCHits];//w/o resolution
  Double_t y0tpc[MaxTPCHits];//w/o resolution
  Double_t z0tpc[MaxTPCHits];//w/o resolution
  //  Double_t resoX[MaxTPCHits];
  // Double_t pxtpc[MaxTPCHits];
  // Double_t pytpc[MaxTPCHits];
  // Double_t pztpc[MaxTPCHits];
  // Double_t pptpc[MaxTPCHits];   // total mometum
  // Double_t masstpc[MaxTPCHits];   // mass TPC
  // Double_t betatpc[MaxTPCHits];
  Double_t edeptpc[MaxTPCHits];
  // Double_t dedxtpc[MaxTPCHits];
  // Double_t slengthtpc[MaxTPCHits];
  // Int_t laytpc[MaxTPCHits];
  // Int_t rowtpc[MaxTPCHits];
  	 Int_t parentID[MaxTPCHits];
  // Int_t iPadtpc[MaxTPCHits];//Pad number (0 origin)

  // Double_t xtpc_pad[MaxTPCHits];//pad center
  // Double_t ytpc_pad[MaxTPCHits];//pad center(dummy)
  // Double_t ztpc_pad[MaxTPCHits];//pad center

  // Double_t dxtpc_pad[MaxTPCHits];//x0tpc - xtpc
  // Double_t dytpc_pad[MaxTPCHits];//y0tpc - ytpc = 0 (dummy)
  // Double_t dztpc_pad[MaxTPCHits];//z0tpc - ztpc


};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
}

//_____________________________________________________________________
int
main(int argc, char **argv)
{
  std::vector<std::string> arg(argv, argv+argc);

  if(!CheckArg(arg))
    return EXIT_FAILURE;
  if(!DstOpen(arg))
    return EXIT_FAILURE;
  if(!gConf.Initialize(arg[kConfFile]))
    return EXIT_FAILURE;
  if(!gConf.InitializeUnpacker())
    return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries(TTreeCont);
  if (max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();

  Int_t ievent = skip;
  for(; ievent<nevent && !CatchSignal::Stop(); ++ievent){
    gCounter.check();
    InitializeEvent();
    if(DstRead(ievent)) tree->Fill();
  }
  std::cout << "#D Event Number: " << std::setw(6)
	    << ievent << std::endl;


  DstClose();

  return EXIT_SUCCESS;
}

//_____________________________________________________________________
Bool_t
dst::InitializeEvent()
{
  event.status   = 0;
  event.evnum = 0;
  event.nhittpc = 0;
  event.ntTpc = 0;

  for(Int_t i=0; i<MaxTPCTracks; ++i){
    event.nhit_track[i] =0;
    event.chisqr[i] =-9999.;
    event.x0[i] =-9999.;
    event.y0[i] =-9999.;
    event.u0[i] =-9999.;
    event.v0[i] =-9999.;
    event.theta[i] =-9999.;
    for(Int_t j=0; j<MaxTPCnHits; ++j){
      event.hitlayer[i][j] =-999;
      event.hitpos_x[i][j] =-9999.;
      event.hitpos_y[i][j] =-9999.;
      event.hitpos_z[i][j] =-9999.;
      event.calpos_x[i][j] =-9999.;
      event.calpos_y[i][j] =-9999.;
      event.calpos_z[i][j] =-9999.;
      event.residual[i][j] =-9999.;
      event.residual_x[i][j] =-9999.;
      event.residual_y[i][j] =-9999.;
      event.residual_z[i][j] =-9999.;
    }
  }


  return true;
}

//_____________________________________________________________________
Bool_t
dst::DstOpen(std::vector<std::string> arg)
{
  Int_t open_file = 0;
  Int_t open_tree = 0;
  for(std::size_t i=0; i<nArgc; ++i){
    if(i==kProcess || i==kConfFile || i==kOutFile) continue;
    open_file += OpenFile(TFileCont[i], arg[i]);
    open_tree += OpenTree(TFileCont[i], TTreeCont[i], TreeName[i]);
  }

  if(open_file!=open_tree || open_file!=nArgc-3)
    return false;
  if(!CheckEntries(TTreeCont))
    return false;

  TFileCont[kOutFile] = new TFile(arg[kOutFile].c_str(), "recreate");

  return true;
}

//_____________________________________________________________________
Bool_t
dst::DstRead(Int_t ievent)
{
  static const std::string func_name("["+class_name+"::"+__func__+"]");


  if(ievent%10000==0){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  HF1(1, event.status++);

  event.evnum = src.evnum;
  event.nhittpc = src.nhittpc;

  // Double_t u = src.pxPrm[0]/src.pzPrm[0];
  // Double_t v = src.pyPrm[0]/src.pzPrm[0];
  // Double_t cost = 1./std::sqrt(1.+u*u+v*v);
  // Double_t theta=std::acos(cost)*math::Rad2Deg();
  // if(theta>20.)
  //   return true;
  DCAnalyzer DCAna;

  //with stable resolution
  // for(Int_t it=0; it<src.nhittpc; ++it){
  //   src.x0tpc[it] += gRandom->Gaus(0,0.2);
  //   src.y0tpc[it] += gRandom->Gaus(0,0.5);
  //   src.z0tpc[it] += gRandom->Gaus(0,0.2);
  // }
  //for test
  // for(Int_t it=0; it<src.nhittpc; ++it){
  //   src.xtpc[it] += gRandom->Gaus(0,0.1);
  //   src.ztpc[it] += gRandom->Gaus(0,0.1);
  // }

  DCAna.DecodeTPCHitsGeant4(src.nhittpc,
                            src.x0tpc, src.y0tpc, src.z0tpc, src.edeptpc,src.idtpc);
                            // src.xtpc, src.ytpc, src.ztpc, src.edeptpc);
  DCAna.TrackSearchTPC();

  Int_t ntTpc = DCAna.GetNTracksTPC();
  HF1(10, ntTpc);
  if(MaxHits<ntTpc){
    std::cout << "#W " << func_name << " "
      	      << "too many ntTpc " << ntTpc << "/" << MaxHits << std::endl;
    ntTpc = MaxHits;
  }
  event.ntTpc = ntTpc;
  for(Int_t it=0; it<ntTpc; ++it){
    auto track= DCAna.GetTrackTPC(it);
    if(!track) continue;
    Int_t nh=track->GetNHit();
    Double_t chisqr    = track->GetChiSquare();
    Double_t x0=track->GetX0(), y0=track->GetY0();
    Double_t u0=track->GetU0(), v0=track->GetV0();
    Double_t theta = track->GetTheta();
    HF1(11, nh);
    HF1(12, chisqr);
    HF1(14, x0);
    HF1(15, y0);
    HF1(16, u0);
    HF1(17, v0);
    HF2(18, x0, u0);
    HF2(19, y0, v0);
    HF2(20, x0, y0);
    event.nhit_track[it]=nh;
    event.chisqr[it]=chisqr;
    event.x0[it]=x0;
    event.y0[it]=y0;
    event.u0[it]=u0;
    event.v0[it]=v0;
    event.theta[it]=theta;

    for(Int_t ih=0; ih<nh; ++ih){
      auto hit = track->GetHit(ih);
      if(!hit) continue;
      Int_t layer = hit->GetLayer();
      TVector3 hitpos = hit->GetLocalHitPos();
      TVector3 calpos = hit->GetLocalCalPos();
      Double_t residual = hit->GetResidual();
      TVector3 res_vect = hit->GetResidualVect();
      HF1(13, layer);
      HF1(100*(layer+1)+15, residual);
      HF1(100*(layer+1)+31, res_vect.X());
      HF1(100*(layer+1)+32, res_vect.Y());
      HF1(100*(layer+1)+33, res_vect.Z());
      event.hitlayer[it][ih] = layer;
      event.hitpos_x[it][ih] = hitpos.x();
      event.hitpos_y[it][ih] = hitpos.y();
      event.hitpos_z[it][ih] = hitpos.z();
      event.calpos_x[it][ih] = calpos.x();
      event.calpos_y[it][ih] = calpos.y();
      event.calpos_z[it][ih] = calpos.z();
      event.residual[it][ih] = residual;
      event.residual_x[it][ih] = res_vect.x();
      event.residual_y[it][ih] = res_vect.y();
      event.residual_z[it][ih] = res_vect.z();
    }
  }


#if 0
  std::cout<<"[event]: "<<std::setw(6)<<ievent<<" ";
  std::cout<<"[nhittpc]: "<<std::setw(2)<<src.nhittpc<<" "<<std::endl;
#endif

  // if(event.nhBh1<=0) return true;
  // HF1(1, event.status++);

  // if(event.nhBh2<=0) return true;
  // HF1(1, event.status++);

  // if(event.nhTof<=0) return true;
  // HF1(1, event.status++);

  // if(event.ntKurama<=0) return true;
  // if(event.ntKurama>MaxTPCHits)
  //   event.ntKurama = MaxTPCHits;

  //HF1(1, event.status++);

  return true;
}

//_____________________________________________________________________
Bool_t
dst::DstClose()
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  const std::size_t n = TFileCont.size();
  for(std::size_t i=0; i<n; ++i){
    if(TTreeCont[i]) delete TTreeCont[i];
    if(TFileCont[i]) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  HB1(1, "Status", 21, 0., 21.);
  HB1(10, "#Tracks TPC", 40, 0., 40.);
  HB1(11, "#Hits of Track TPC", 50, 0., 50.);
  HB1(12, "Chisqr TPC", 500, 0., 500.);
  HB1(13, "LayerId TPC", 35, 0., 35.);
  HB1(14, "X0 TPC", 400, -100., 100.);
  HB1(15, "Y0 TPC", 400, -100., 100.);
  HB1(16, "U0 TPC", 200, -0.20, 0.20);
  HB1(17, "V0 TPC", 200, -0.20, 0.20);
  HB2(18, "U0%X0 TPC", 100, -100., 100., 100, -0.20, 0.20);
  HB2(19, "V0%Y0 TPC", 100, -100., 100., 100, -0.20, 0.20);
  HB2(20, "X0%Y0 TPC", 100, -100., 100., 100, -100, 100);

  for(Int_t i=1; i<=NumOfLayersTPC; ++i){
    // Tracking Histgrams
    TString title11 = Form("HitPat TPC%2d [Track]", i);
    TString title14 = Form("Position TPC%2d", i);
    TString title15 = Form("Residual TPC%2d", i);
    TString title16 = Form("Resid%%Pos TPC%2d", i);
    TString title17 = Form("Y%%Xcal TPC%2d", i);
    TString title31 = Form("ResidualX TPC%2d", i);
    TString title32 = Form("ResidualY TPC%2d", i);
    TString title33 = Form("ResidualZ TPC%2d", i);
    HB1(100*i+11, title11, 400, 0., 400.);
    HB1(100*i+14, title14, 200, -250., 250.);
    HB1(100*i+15, title15, 200, 0.0, 10.0);
    HB2(100*i+16, title16, 250, -250., 250., 100, -1.0, 1.0);
    HB2(100*i+17, title17, 100, -250., 250., 100, -250., 250.);
    HB1(100*i+31, title31, 200, -2.0, 2.0);
    HB1(100*i+32, title32, 200, -2.0, 2.0);
    HB1(100*i+33, title33, 200, -2.0, 2.0);
  }

  HBTree("g4tpc", "tree of DstTPCTracking using Geant4 input");

  tree->Branch("status", &event.status, "status/I");
  tree->Branch("evnum", &event.evnum, "evnum/I");
  tree->Branch("nhittpc",&event.nhittpc,"nhittpc/I");
  tree->Branch("ntTpc",&event.ntTpc,"ntTpc/I");

  tree->Branch("nhit_track",event.nhit_track,"nhit_track[ntTpc]/I");
  tree->Branch("chisqr",event.chisqr,"chisqr[ntTpc]/D");
  tree->Branch("x0",event.x0,"x0[ntTpc]/D");
  tree->Branch("y0",event.y0,"y0[ntTpc]/D");
  tree->Branch("u0",event.u0,"u0[ntTpc]/D");
  tree->Branch("v0",event.v0,"v0[ntTpc]/D");
  tree->Branch("theta",event.theta,"theta[ntTpc]/D");

  tree->Branch("hitlayer",event.hitlayer,"hitlayer[ntTpc][32]/I");
  tree->Branch("hitpos_x",event.hitpos_x,"hitpos_x[ntTpc][32]/D");
  tree->Branch("hitpos_y",event.hitpos_y,"hitpos_y[ntTpc][32]/D");
  tree->Branch("hitpos_z",event.hitpos_z,"hitpos_z[ntTpc][32]/D");
  tree->Branch("calpos_x",event.calpos_x,"calpos_x[ntTpc][32]/D");
  tree->Branch("calpos_y",event.calpos_y,"calpos_y[ntTpc][32]/D");
  tree->Branch("calpos_z",event.calpos_z,"calpos_z[ntTpc][32]/D");
  tree->Branch("residual",event.residual,"residual[ntTpc][32]/D");
  tree->Branch("residual_x",event.residual_x,"residual_x[ntTpc][32]/D");
  tree->Branch("residual_y",event.residual_y,"residual_y[ntTpc][32]/D");
  tree->Branch("residual_z",event.residual_z,"residual_z[ntTpc][32]/D");

  tree->Branch("nPrm",&src.nhittpc,"nPrm/I");
  tree->Branch("xPrm",src.xPrm,"xPrm[nPrm]/D");
  tree->Branch("yPrm",src.yPrm,"yPrm[nPrm]/D");
  tree->Branch("zPrm",src.zPrm,"zPrm[nPrm]/D");
  tree->Branch("pxPrm",src.pxPrm,"pxPrm[nPrm]/D");
  tree->Branch("pyPrm",src.pyPrm,"pyPrm[nPrm]/D");
  tree->Branch("pzPrm",src.pzPrm,"pzPrm[nPrm]/D");



  ////////// Bring Address From Dst
  TTreeCont[kTPCGeant]->SetBranchStatus("*", 0);
  TTreeCont[kTPCGeant]->SetBranchStatus("evnum",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("nhittpc",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("nhPrm",  1);

  TTreeCont[kTPCGeant]->SetBranchStatus("xPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("yPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("zPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pxPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pyPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pzPrm",  1);

  TTreeCont[kTPCGeant]->SetBranchStatus("ititpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("idtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("xtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ytpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ztpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("x0tpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("y0tpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("z0tpc", 1);
  //  TTreeCont[kTPCGeant]->SetBranchStatus("resoX", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("pxtpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("pytpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("pztpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("pptpc", 1);   // total mometum
  // TTreeCont[kTPCGeant]->SetBranchStatus("masstpc", 1);   // mass TPC
  // TTreeCont[kTPCGeant]->SetBranchStatus("betatpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("edeptpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("dedxtpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("slengthtpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("laytpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("rowtpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("parentID", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("iPadtpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("xtpc_pad", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("ytpc_pad", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("ztpc_pad", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("dxtpc_pad", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("dytpc_pad", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("dztpc_pad", 1);


  TTreeCont[kTPCGeant]->SetBranchAddress("evnum", &src.evnum);
  TTreeCont[kTPCGeant]->SetBranchAddress("nhittpc", &src.nhittpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("nhPrm", &src.nhPrm);

  TTreeCont[kTPCGeant]->SetBranchAddress("xPrm", src.xPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("yPrm", src.yPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("zPrm", src.zPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxPrm", src.pxPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("pyPrm", src.pyPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzPrm", src.pzPrm);

  TTreeCont[kTPCGeant]->SetBranchAddress("ititpc", src.ititpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("idtpc", src.idtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("xtpc", src.xtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("ytpc", src.ytpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("ztpc", src.ztpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("x0tpc", src.x0tpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("y0tpc", src.y0tpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("z0tpc", src.z0tpc);
  //  TTreeCont[kTPCGeant]->SetBranchAddress("resoX", src.resoX);
  // TTreeCont[kTPCGeant]->SetBranchAddress("pxtpc", src.pxtpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("pytpc", src.pytpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("pztpc", src.pztpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("pptpc", src.pptpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("masstpc", src.masstpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("betatpc", src.betatpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("edeptpc", src.edeptpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("dedxtpc", src.dedxtpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("slengthtpc", src.slengthtpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("laytpc", src.laytpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("rowtpc", src.rowtpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("parentID", src.parentID);
  // TTreeCont[kTPCGeant]->SetBranchAddress("iPadtpc", src.iPadtpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("xtpc_pad", src.xtpc_pad);
  // TTreeCont[kTPCGeant]->SetBranchAddress("ytpc_pad", src.ytpc_pad);
  // TTreeCont[kTPCGeant]->SetBranchAddress("ztpc_pad", src.ztpc_pad);
  // TTreeCont[kTPCGeant]->SetBranchAddress("dxtpc_pad", src.dxtpc_pad);
  // TTreeCont[kTPCGeant]->SetBranchAddress("dytpc_pad", src.dytpc_pad);
  // TTreeCont[kTPCGeant]->SetBranchAddress("dztpc_pad", src.dztpc_pad);


  return true;
}

//_____________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")   &&
     InitializeParameter<UserParamMan>("USER") &&
     InitializeParameter<HodoPHCMan>("HDPHC"));
}

//_____________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

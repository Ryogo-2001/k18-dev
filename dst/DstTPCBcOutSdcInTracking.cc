// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DatabasePDG.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DstHelper.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCCluster.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrack.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "UserParamMan.hh"

#define Exclusive 1

namespace
{
  using namespace root;
  using namespace dst;
  using hddaq::unpacker::GUnpacker;
  const auto qnan = TMath::QuietNaN();
  const auto& gUnpacker = GUnpacker::get_instance();
  auto&       gConf = ConfMan::GetInstance();
  const auto& gGeom = DCGeomMan::GetInstance();
  const auto& gUser = UserParamMan::GetInstance();
  const auto& gCounter = debug::ObjectCounter::GetInstance();
  const Double_t& zK18Target = gGeom.LocalZ("K18Target");
  const Double_t& zTarget = gGeom.LocalZ("Target");
  const Int_t& IdHS = gGeom.DetectorId("HS");
}

namespace dst
{
  enum kArgc
    {
      kProcess, kConfFile,
      kTpc, kBcSdc, kOutFile, nArgc
    };
  std::vector<TString> ArgName =
    { "[Process]", "[ConfFile]", "[TPCTracking]", "[BcOutSdcInTracking]", "[OutFile]" };
  std::vector<TString> TreeName = { "", "", "tpc", "bcsdc", "" };
  std::vector<TFile*> TFileCont;
  std::vector<TTree*> TTreeCont;
  std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________________
struct Event
{
  Int_t status;
  Int_t runnum;
  Int_t evnum;
  std::vector<Int_t> trigpat;
  std::vector<Int_t> trigflag;

  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Double_t> chisqrTpc;
  std::vector<Double_t> x0Tpc;
  std::vector<Double_t> y0Tpc;
  std::vector<Double_t> u0Tpc;
  std::vector<Double_t> v0Tpc;
  std::vector<Double_t> theta;
  std::vector<Double_t> dE;
  std::vector<Double_t> dEdx; //reference dedx
  std::vector<std::vector<Double_t>> hitlayer;
  std::vector<std::vector<Double_t>> hitpos_x;
  std::vector<std::vector<Double_t>> hitpos_y;
  std::vector<std::vector<Double_t>> hitpos_z;
  std::vector<std::vector<Double_t>> calpos_x;
  std::vector<std::vector<Double_t>> calpos_y;
  std::vector<std::vector<Double_t>> calpos_z;
  std::vector<std::vector<Double_t>> residual;
  std::vector<std::vector<Double_t>> residual_x;
  std::vector<std::vector<Double_t>> residual_y;
  std::vector<std::vector<Double_t>> residual_z;
  std::vector<std::vector<Double_t>> residual_horizontal;
  std::vector<std::vector<Double_t>> residual_vertical;
  std::vector<std::vector<Double_t>> resolution_x;
  std::vector<std::vector<Double_t>> resolution_y;
  std::vector<std::vector<Double_t>> resolution_z;
  std::vector<std::vector<Double_t>> resolution_horizontal;
  std::vector<std::vector<Double_t>> resolution_vertical;
  std::vector<std::vector<Double_t>> pathhit;
  std::vector<std::vector<Double_t>> theta_diff;
  std::vector<std::vector<Double_t>> track_cluster_de;
  std::vector<std::vector<Double_t>> track_cluster_size;
  std::vector<std::vector<Double_t>> track_cluster_mrow;
  std::vector<std::vector<Double_t>> track_cluster_de_center;
  std::vector<std::vector<Double_t>> track_cluster_x_center;
  std::vector<std::vector<Double_t>> track_cluster_y_center;
  std::vector<std::vector<Double_t>> track_cluster_z_center;
  std::vector<std::vector<Double_t>> track_cluster_row_center;

  //exclusive
  std::vector<std::vector<Double_t>> exresidual;
  std::vector<std::vector<Double_t>> exresidual_x;
  std::vector<std::vector<Double_t>> exresidual_y;
  std::vector<std::vector<Double_t>> exresidual_z;
  std::vector<std::vector<Double_t>> exresidual_horizontal;
  std::vector<std::vector<Double_t>> exresidual_vertical;

  Double_t btof;

  Int_t ntBcSdc;
  std::vector<Int_t> nhBcSdc;
  std::vector<Double_t> chisqrBcSdc;
  std::vector<Double_t> x0BcSdc;
  std::vector<Double_t> y0BcSdc;
  std::vector<Double_t> u0BcSdc;
  std::vector<Double_t> v0BcSdc;
  std::vector<Double_t> xtgtBcSdc;
  std::vector<Double_t> ytgtBcSdc;
  std::vector<Double_t> utgtBcSdc;
  std::vector<Double_t> vtgtBcSdc;

  Int_t ntBcOut;
  std::vector<Int_t> nhBcOut;
  std::vector<Double_t> chisqrBcOut;
  std::vector<Double_t> x0BcOut;
  std::vector<Double_t> y0BcOut;
  std::vector<Double_t> u0BcOut;
  std::vector<Double_t> v0BcOut;
  std::vector<Double_t> xtgtBcOut;
  std::vector<Double_t> ytgtBcOut;
  std::vector<Double_t> utgtBcOut;
  std::vector<Double_t> vtgtBcOut;

  Int_t ntSdcIn;
  std::vector<Int_t> nhSdcIn;
  std::vector<Double_t> chisqrSdcIn;
  std::vector<Double_t> x0SdcIn;
  std::vector<Double_t> y0SdcIn;
  std::vector<Double_t> u0SdcIn;
  std::vector<Double_t> v0SdcIn;
  std::vector<Double_t> xtgtSdcIn;
  std::vector<Double_t> ytgtSdcIn;
  std::vector<Double_t> utgtSdcIn;
  std::vector<Double_t> vtgtSdcIn;

  void clear()
  {
    runnum = 0;
    evnum = 0;
    status = 0;
    trigpat.clear();
    trigflag.clear();

    ntTpc = 0;
    nhtrack.clear();
    chisqrTpc.clear();
    x0Tpc.clear();
    y0Tpc.clear();
    u0Tpc.clear();
    v0Tpc.clear();
    theta.clear();
    hitlayer.clear();
    hitpos_x.clear();
    hitpos_y.clear();
    hitpos_z.clear();
    calpos_x.clear();
    calpos_y.clear();
    calpos_z.clear();
    residual.clear();
    residual_x.clear();
    residual_y.clear();
    residual_z.clear();
    residual_horizontal.clear();
    residual_vertical.clear();
    resolution_x.clear();
    resolution_y.clear();
    resolution_z.clear();
    resolution_horizontal.clear();
    resolution_vertical.clear();

    dE.clear();
    dEdx.clear();
    pathhit.clear();
    theta_diff.clear();
    track_cluster_de.clear();
    track_cluster_size.clear();
    track_cluster_mrow.clear();
    track_cluster_de_center.clear();
    track_cluster_x_center.clear();
    track_cluster_y_center.clear();
    track_cluster_z_center.clear();
    track_cluster_row_center.clear();

    exresidual.clear();
    exresidual_x.clear();
    exresidual_y.clear();
    exresidual_z.clear();
    exresidual_horizontal.clear();
    exresidual_vertical.clear();

    btof = qnan;

    ntBcSdc = 0;
    nhBcSdc.clear();
    chisqrBcSdc.clear();
    x0BcSdc.clear();
    y0BcSdc.clear();
    u0BcSdc.clear();
    v0BcSdc.clear();
    xtgtBcSdc.clear();
    ytgtBcSdc.clear();
    utgtBcSdc.clear();
    vtgtBcSdc.clear();

    ntBcOut = 0;
    nhBcOut.clear();
    chisqrBcOut.clear();
    x0BcOut.clear();
    y0BcOut.clear();
    u0BcOut.clear();
    v0BcOut.clear();
    xtgtBcOut.clear();
    ytgtBcOut.clear();
    utgtBcOut.clear();
    vtgtBcOut.clear();

    ntSdcIn = 0;
    nhSdcIn.clear();
    chisqrSdcIn.clear();
    x0SdcIn.clear();
    y0SdcIn.clear();
    u0SdcIn.clear();
    v0SdcIn.clear();
    xtgtSdcIn.clear();
    ytgtSdcIn.clear();
    utgtSdcIn.clear();
    vtgtSdcIn.clear();
  }
};

//_____________________________________________________________________________
struct Src
{

  TTreeReaderValue<Int_t>* status;
  TTreeReaderValue<Int_t>* runnum;
  TTreeReaderValue<Int_t>* evnum;
  TTreeReaderValue<std::vector<Int_t>>* trigpat;
  TTreeReaderValue<std::vector<Int_t>>* trigflag;

  TTreeReaderValue<Int_t>* ntTpc;
  TTreeReaderValue<std::vector<Int_t>>* nhtrack;
  TTreeReaderValue<std::vector<Double_t>>* chisqrTpc;
  TTreeReaderValue<std::vector<Double_t>>* x0Tpc;
  TTreeReaderValue<std::vector<Double_t>>* y0Tpc;
  TTreeReaderValue<std::vector<Double_t>>* u0Tpc;
  TTreeReaderValue<std::vector<Double_t>>* v0Tpc;
  TTreeReaderValue<std::vector<Double_t>>* theta;
  TTreeReaderValue<std::vector<Double_t>>* dE;
  TTreeReaderValue<std::vector<Double_t>>* dEdx;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* theta_diff;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* pathhit;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_horizontal;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_vertical;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_horizontal;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_vertical;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitlayer;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitpos_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitpos_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitpos_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* calpos_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* calpos_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* calpos_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_size;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_mrow;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_x_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_y_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_z_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_row_center;

  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual_horizontal;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual_vertical;

  Double_t btof;

  Int_t ntBcSdc;
  Int_t nhBcSdc[MaxHits];
  Double_t chisqrBcSdc[MaxHits];
  Double_t x0BcSdc[MaxHits];
  Double_t y0BcSdc[MaxHits];
  Double_t u0BcSdc[MaxHits];
  Double_t v0BcSdc[MaxHits];

  Int_t ntBcOut;
  Int_t nhBcOut[MaxHits];
  Double_t chisqrBcOut[MaxHits];
  Double_t x0BcOut[MaxHits];
  Double_t y0BcOut[MaxHits];
  Double_t u0BcOut[MaxHits];
  Double_t v0BcOut[MaxHits];

  Int_t ntSdcIn;
  Int_t nhSdcIn[MaxHits];
  Double_t chisqrSdcIn[MaxHits];
  Double_t x0SdcIn[MaxHits];
  Double_t y0SdcIn[MaxHits];
  Double_t u0SdcIn[MaxHits];
  Double_t v0SdcIn[MaxHits];

};

namespace root
{
  Event  event;
  Src    src;
  TH1   *h[MaxHist];
  TTree *tree;
  enum eDetHid {
    PadHid    = 100000,
    TPCCalibHid =1000,
  };

  TVector3 Global2LocalPos(TVector3 GlobalPos, TVector3 DetectorGlobalPos,
			   Double_t tilt, Double_t RA1, Double_t RA2){

    Double_t ct0 = TMath::Cos(tilt*TMath::DegToRad());
    Double_t st0 = TMath::Sin(tilt*TMath::DegToRad());
    Double_t ct1 = TMath::Cos(RA1*TMath::DegToRad());
    Double_t st1 = TMath::Sin(RA1*TMath::DegToRad());
    Double_t ct2 = TMath::Cos(RA2*TMath::DegToRad());
    Double_t st2 = TMath::Sin(RA2*TMath::DegToRad());

    Double_t dsdx =  ct0*ct2+st0*st1*st2;
    Double_t dsdy =  st0*ct1;
    Double_t dsdz = -ct0*st2+st0*st1*ct2;

    Double_t dtdx = -st0*ct2+ct0*st1*st2;
    Double_t dtdy =  ct0*ct1;
    Double_t dtdz =  st0*st2+ct0*st1*ct2;

    Double_t dudx =  ct1*st2;
    Double_t dudy = -st1;
    Double_t dudz =  ct1*ct2;

    Double_t x
      = dsdx*(GlobalPos.x()-DetectorGlobalPos.x())
      + dsdy*(GlobalPos.y()-DetectorGlobalPos.y())
      + dsdz*(GlobalPos.z()-DetectorGlobalPos.z());
    Double_t y
      = dtdx*(GlobalPos.x()-DetectorGlobalPos.x())
      + dtdy*(GlobalPos.y()-DetectorGlobalPos.y())
      + dtdz*(GlobalPos.z()-DetectorGlobalPos.z());
    Double_t z
      = dudx*(GlobalPos.x()-DetectorGlobalPos.x())
      + dudy*(GlobalPos.y()-DetectorGlobalPos.y())
      + dudz*(GlobalPos.z()-DetectorGlobalPos.z());

    return TVector3(x, y, z);
  }

  TVector3 Local2GlobalPos(TVector3 LocalPos, TVector3 DetectorGlobalPos,
			    Double_t tilt, Double_t RA1, Double_t RA2){

    Double_t ct0 = TMath::Cos(tilt*TMath::DegToRad());
    Double_t st0 = TMath::Sin(tilt*TMath::DegToRad());
    Double_t ct1 = TMath::Cos(RA1*TMath::DegToRad());
    Double_t st1 = TMath::Sin(RA1*TMath::DegToRad());
    Double_t ct2 = TMath::Cos(RA2*TMath::DegToRad());
    Double_t st2 = TMath::Sin(RA2*TMath::DegToRad());

    Double_t dxds =  ct0*ct2+st0*st1*st2;
    Double_t dxdt = -st0*ct2+ct0*st1*st2;
    Double_t dxdu =  ct1*st2;

    Double_t dyds =  st0*ct1;
    Double_t dydt =  ct0*ct1;
    Double_t dydu = -st1;

    Double_t dzds = -ct0*st2+st0*st1*ct2;
    Double_t dzdt =  st0*st2+ct0*st1*ct2;
    Double_t dzdu =  ct1*ct2;

    Double_t x = dxds*LocalPos.x() + dxdt*LocalPos.y()
      + dxdu*LocalPos.z() + DetectorGlobalPos.x();
    Double_t y = dyds*LocalPos.x() + dydt*LocalPos.y()
      + dydu*LocalPos.z() + DetectorGlobalPos.y();
    Double_t z = dzds*LocalPos.x() + dzdt*LocalPos.y()
      + dzdu*LocalPos.z() + DetectorGlobalPos.z();

    return TVector3(x, y, z);
  }

}

//_____________________________________________________________________________
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
  if(skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries(TTreeCont);
  if(max_loop > 0) nevent = skip + max_loop;

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

//_____________________________________________________________________________
Bool_t
dst::InitializeEvent()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstOpen(std::vector<std::string> arg)
{
  int open_file = 0;
  int open_tree = 0;
  for(Int_t i=0; i<nArgc; ++i){
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

//_____________________________________________________________________________
Bool_t
dst::DstRead(int ievent)
{

  const Double_t& zBC3 = gGeom.LocalZ("BC3-X1");
  const ThreeVector& posBcOut = gGeom.GetGlobalPosition("BC3-X1");
  TVector3 posV0 = posBcOut + TVector3(0., 0., -zBC3);

  const Double_t& zSdcIn = gGeom.LocalZ("SDC1-U2");
  const ThreeVector& posSdcIn = gGeom.GetGlobalPosition("SDC1-U2");
  TVector3 posKurama = posSdcIn + TVector3(0., 0., -zSdcIn);

  const ThreeVector& posHS = gGeom.GetGlobalPosition("K18HS");
  const Double_t tilt = gGeom.GetTiltAngle("K18HS");
  const Double_t RA1 = gGeom.GetRotAngle1("K18HS");
  const Double_t RA2 = gGeom.GetRotAngle2("K18HS");

  if(ievent%100==0){
  //if(ievent%1==0){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;

  HF1(1, event.status++);

  event.ntTpc = **src.ntTpc;
  event.nhtrack = **src.nhtrack;
  event.chisqrTpc = **src.chisqrTpc;
  event.x0Tpc = **src.x0Tpc;
  event.y0Tpc = **src.y0Tpc;
  event.u0Tpc = **src.u0Tpc;
  event.v0Tpc = **src.v0Tpc;
  event.theta = **src.theta;
  event.dE = **src.dE;
  event.dEdx = **src.dEdx;
  event.theta_diff = **src.theta_diff;
  event.pathhit = **src.pathhit;
  event.residual = **src.residual;
  event.residual_x = **src.residual_x;
  event.residual_y = **src.residual_y;
  event.residual_z = **src.residual_z;
  event.residual_horizontal = **src.residual_horizontal;
  event.residual_vertical = **src.residual_vertical;
  event.resolution_x = **src.resolution_x;
  event.resolution_y = **src.resolution_y;
  event.resolution_z = **src.resolution_z;
  event.resolution_horizontal = **src.resolution_horizontal;
  event.resolution_vertical = **src.resolution_vertical;

  event.hitlayer = **src.hitlayer;
  event.hitpos_x = **src.hitpos_x;
  event.hitpos_y = **src.hitpos_y;
  event.hitpos_z = **src.hitpos_z;
  event.calpos_x = **src.calpos_x;
  event.calpos_y = **src.calpos_y;
  event.calpos_z = **src.calpos_z;
  event.track_cluster_de = **src.track_cluster_de;
  event.track_cluster_size = **src.track_cluster_size;
  event.track_cluster_mrow = **src.track_cluster_mrow;
  event.track_cluster_de_center = **src.track_cluster_de_center;
  event.track_cluster_x_center = **src.track_cluster_x_center;
  event.track_cluster_y_center = **src.track_cluster_y_center;
  event.track_cluster_z_center = **src.track_cluster_z_center;
  event.track_cluster_row_center = **src.track_cluster_row_center;
  event.exresidual = **src.exresidual;
  event.exresidual_x = **src.exresidual_x;
  event.exresidual_y = **src.exresidual_y;
  event.exresidual_z = **src.exresidual_z;
  event.exresidual_horizontal = **src.exresidual_horizontal;
  event.exresidual_vertical = **src.exresidual_vertical;

  HF1(1, event.status++);

  event.btof = src.btof;
  event.ntBcOut = src.ntBcOut;
  for(Int_t it=0; it<src.ntBcOut; ++it){
    event.nhBcOut.push_back(src.nhBcOut[it]);
    event.chisqrBcOut.push_back(src.chisqrBcOut[it]);
    event.x0BcOut.push_back(src.x0BcOut[it]);
    event.y0BcOut.push_back(src.y0BcOut[it]);
    event.u0BcOut.push_back(src.u0BcOut[it]);
    event.v0BcOut.push_back(src.v0BcOut[it]);

    Double_t utgt = src.u0BcOut[it];
    Double_t vtgt = src.v0BcOut[it];
    Double_t xtgt = utgt*zK18Target + src.x0BcOut[it];
    Double_t ytgt = vtgt*zK18Target + src.y0BcOut[it];
    event.xtgtBcOut.push_back(xtgt);
    event.ytgtBcOut.push_back(ytgt);
    event.utgtBcOut.push_back(utgt);
    event.vtgtBcOut.push_back(vtgt);
  }

  HF1(1, event.status++);

  event.ntBcSdc = src.ntBcSdc;
  for(Int_t it=0; it<src.ntBcSdc; ++it){
    event.nhBcSdc.push_back(src.nhBcSdc[it]);
    event.chisqrBcSdc.push_back(src.chisqrBcSdc[it]);
    event.x0BcSdc.push_back(src.x0BcSdc[it]);
    event.y0BcSdc.push_back(src.y0BcSdc[it]);
    event.u0BcSdc.push_back(src.u0BcSdc[it]);
    event.v0BcSdc.push_back(src.v0BcSdc[it]);

    Double_t utgt = src.u0BcSdc[it];
    Double_t vtgt = src.v0BcSdc[it];
    Double_t xtgt = utgt*zK18Target + src.x0BcSdc[it];
    Double_t ytgt = vtgt*zK18Target + src.y0BcSdc[it];
    event.xtgtBcSdc.push_back(xtgt);
    event.ytgtBcSdc.push_back(ytgt);
    event.utgtBcSdc.push_back(utgt);
    event.vtgtBcSdc.push_back(vtgt);
  }

  HF1(1, event.status++);

  event.ntSdcIn = src.ntSdcIn;
  for(Int_t it=0; it<src.ntSdcIn; ++it){
    event.nhSdcIn.push_back(src.nhSdcIn[it]);
    event.chisqrSdcIn.push_back(src.chisqrSdcIn[it]);
    event.x0SdcIn.push_back(src.x0SdcIn[it]);
    event.y0SdcIn.push_back(src.y0SdcIn[it]);
    event.u0SdcIn.push_back(src.u0SdcIn[it]);
    event.v0SdcIn.push_back(src.v0SdcIn[it]);

    Double_t utgt = src.u0SdcIn[it];
    Double_t vtgt = src.v0SdcIn[it];
    Double_t xtgt = utgt*zTarget + src.x0SdcIn[it];
    Double_t ytgt = vtgt*zTarget + src.y0SdcIn[it];
    event.xtgtSdcIn.push_back(xtgt);
    event.ytgtSdcIn.push_back(ytgt);
    event.utgtSdcIn.push_back(utgt);
    event.vtgtSdcIn.push_back(vtgt);
  }

  HF1(1, event.status++);

  if(event.ntTpc!=1 || event.ntBcOut!=1 || event.ntSdcIn!=1) return true;
  if(event.chisqrTpc[0]>10. || event.chisqrBcOut[0]>5. || event.chisqrSdcIn[0]>5.) return true;

  TVector3 xytgt(event.x0Tpc[0], event.y0Tpc[0], tpc::ZTarget); //Target local position in the HS
  TVector3 xytgtGlobal = Local2GlobalPos(xytgt, posHS, tilt, RA1, RA2);
  std::vector<Double_t> xyuvBcOut_HScoor(4);
  {
    //On the HS coordinate
    Double_t zpos = Global2LocalPos(xytgtGlobal, posV0, 0., 0., 0.).z();
    TVector3 localpos(src.u0BcOut[0]*zpos + src.x0BcOut[0], src.v0BcOut[0]*zpos + src.y0BcOut[0], zpos);
    TVector3 globalpos = Local2GlobalPos(localpos, posV0, 0., 0., 0.);
    xyuvBcOut_HScoor[0] = Global2LocalPos(globalpos, posHS, tilt, RA1, RA2).x();
    xyuvBcOut_HScoor[1] = Global2LocalPos(globalpos, posHS, tilt, RA1, RA2).y();
    xyuvBcOut_HScoor[2] = src.u0BcOut[0] - TMath::Tan(RA2*TMath::DegToRad());
    xyuvBcOut_HScoor[3] = src.v0BcOut[0] + TMath::Tan(RA1*TMath::DegToRad());
  }

  std::vector<Double_t> xyuvSdcIn_HScoor(4);
  {
    //On the HS coordinate
    Double_t zpos = Global2LocalPos(xytgtGlobal, posKurama, 0., 0., 0.).z();
    TVector3 localpos(src.u0SdcIn[0]*zpos + src.x0SdcIn[0], src.v0SdcIn[0]*zpos + src.y0SdcIn[0], zpos);
    TVector3 globalpos = Local2GlobalPos(localpos, posKurama, 0., 0., 0.);
    xyuvSdcIn_HScoor[0] = Global2LocalPos(globalpos, posHS, tilt, RA1, RA2).x();
    xyuvSdcIn_HScoor[1] = Global2LocalPos(globalpos, posHS, tilt, RA1, RA2).y();
    xyuvSdcIn_HScoor[2] = src.u0SdcIn[0] - TMath::Tan(RA2*TMath::DegToRad());
    xyuvSdcIn_HScoor[3] = src.v0SdcIn[0] + TMath::Tan(RA1*TMath::DegToRad());
  }

  HF1(10, event.ntTpc);
  for(Int_t it=0; it<event.ntTpc; ++it){
    HF1(11, event.nhtrack[it]);
    HF1(12, event.chisqrTpc[it]);
    HF1(14, event.x0Tpc[it]);
    HF1(15, event.y0Tpc[it]);
    HF1(16, event.u0Tpc[it]);
    HF1(17, event.v0Tpc[it]);
    HF2(18, event.x0Tpc[it], event.u0Tpc[it]);
    HF2(19, event.y0Tpc[it], event.v0Tpc[it]);
    HF2(20, event.x0Tpc[it], event.y0Tpc[it]);
    HF2(21, event.x0Tpc[it], atan(event.u0Tpc[it])*180./acos(-1.));
    HF2(22, event.y0Tpc[it], atan(event.v0Tpc[it])*180./acos(-1.));
  }

  for(Int_t it=0; it<event.ntTpc; ++it){
    for(Int_t ih=0; ih<event.nhtrack[it]; ++ih){
      HF2(TPCCalibHid, event.hitlayer[it][ih], event.residual_y[it][ih]);
      if(event.resolution_x[it][ih] > 0.9e5 && event.resolution_y[it][ih] > 0.9e5 && event.resolution_x[it][ih] > 0.9e5) continue; // exclude dummy hits
      HF2(1000*(event.hitlayer[it][ih]+1) + 1, event.xtgtSdcIn[0], event.residual_x[it][ih]);
      HF2(1000*(event.hitlayer[it][ih]+1) + 2, event.ytgtSdcIn[0], event.residual_y[it][ih]);
      HF2(1000*(event.hitlayer[it][ih]+1) + 3, event.xtgtSdcIn[0], event.residual_z[it][ih]);
      HF2(1000*(event.hitlayer[it][ih]+1) + 4, event.xtgtSdcIn[0], event.residual_horizontal[it][ih]);
      HF2(1000*(event.hitlayer[it][ih]+1) + 5, event.xtgtSdcIn[0], event.residual_vertical[it][ih]);
      HF2(1000*(event.hitlayer[it][ih]+1) + 6, event.xtgtSdcIn[0], TMath::Sqrt(event.residual_x[it][ih]*event.residual_x[it][ih] + event.residual_z[it][ih]*event.residual_z[it][ih]));

      HF1(1000*(event.hitlayer[it][ih]+1) + 11, event.residual_x[it][ih]/event.resolution_x[it][ih]);
      HF1(1000*(event.hitlayer[it][ih]+1) + 12, event.residual_y[it][ih]/event.resolution_y[it][ih]);
      HF1(1000*(event.hitlayer[it][ih]+1) + 13, event.residual_z[it][ih]/event.resolution_z[it][ih]);
      HF1(1000*(event.hitlayer[it][ih]+1) + 14, event.residual_horizontal[it][ih]/event.resolution_horizontal[it][ih]);
      HF1(1000*(event.hitlayer[it][ih]+1) + 15, event.residual_vertical[it][ih]/event.resolution_vertical[it][ih]);

      HF1(1000*(event.hitlayer[it][ih]+1) + 101, event.residual_x[it][ih]);
      HF1(1000*(event.hitlayer[it][ih]+1) + 102, event.residual_y[it][ih]);
      HF1(1000*(event.hitlayer[it][ih]+1) + 103, event.residual_z[it][ih]);
      HF1(1000*(event.hitlayer[it][ih]+1) + 104, event.residual_horizontal[it][ih]);
      HF1(1000*(event.hitlayer[it][ih]+1) + 105, event.residual_vertical[it][ih]);
      HF1(1000*(event.hitlayer[it][ih]+1) + 106, TMath::Sqrt(event.residual_x[it][ih]*event.residual_x[it][ih] + event.residual_z[it][ih]*event.residual_z[it][ih]));

      ThreeVector lpos_tpc(event.hitpos_x[it][ih],  event.hitpos_y[it][ih],  event.hitpos_z[it][ih]);
      ThreeVector gpos = gGeom.Local2GlobalPos(IdHS, lpos_tpc);
      ThreeVector lpos_bc = gpos - posV0;
      ThreeVector lpos_sdcin = gpos - posKurama;

      HF1(1000*(event.hitlayer[it][ih]+1) + 201, -(src.u0BcOut[it]*lpos_bc.z() + src.x0BcOut[it]) + lpos_bc.x());
      HF1(1000*(event.hitlayer[it][ih]+1) + 202, -(src.v0BcOut[it]*lpos_bc.z() + src.y0BcOut[it]) + lpos_bc.y());
      HF1(1000*(event.hitlayer[it][ih]+1) + 203, -(src.u0SdcIn[it]*lpos_sdcin.z() + src.x0SdcIn[it]) + lpos_sdcin.x());
      HF1(1000*(event.hitlayer[it][ih]+1) + 204, -(src.v0SdcIn[it]*lpos_sdcin.z() + src.y0SdcIn[it]) + lpos_sdcin.y());
    }
  }

  HF1(110, event.ntBcOut);
  for(Int_t it=0; it<src.ntBcOut; ++it){
    HF1(111, event.nhBcOut[it]);
    HF1(112, event.chisqrBcOut[it]);
    HF1(114, event.x0BcOut[it]);
    HF1(115, event.y0BcOut[it]);
    HF1(116, event.u0BcOut[it]);
    HF1(117, event.v0BcOut[it]);
    HF2(118, event.utgtBcOut[it], event.xtgtBcOut[it]);
    HF2(119, event.vtgtBcOut[it], event.ytgtBcOut[it]);
    HF2(120, event.xtgtBcOut[it], event.ytgtBcOut[it]);
    HF1(121, event.xtgtBcOut[it]);
    HF1(122, event.ytgtBcOut[it]);
    HF1(123, event.utgtBcOut[it]);
    HF1(124, event.vtgtBcOut[it]);
  }

  HF1(210, event.ntSdcIn);
  for(Int_t it=0; it<src.ntBcSdc; ++it){
    HF1(211, event.nhSdcIn[it]);
    HF1(212, event.chisqrSdcIn[it]);
    HF1(214, event.x0SdcIn[it]);
    HF1(215, event.y0SdcIn[it]);
    HF1(216, event.u0SdcIn[it]);
    HF1(217, event.v0SdcIn[it]);
    HF2(218, event.utgtSdcIn[it], event.xtgtSdcIn[it]);
    HF2(219, event.vtgtSdcIn[it], event.ytgtSdcIn[it]);
    HF2(220, event.xtgtSdcIn[it], event.ytgtSdcIn[it]);
    HF1(221, event.xtgtSdcIn[it]);
    HF1(222, event.ytgtSdcIn[it]);
    HF1(223, event.utgtSdcIn[it]);
    HF1(224, event.vtgtSdcIn[it]);
  }

  HF1(310, event.ntBcSdc);
  for(Int_t it=0; it<src.ntBcSdc; ++it){
    HF1(311, event.nhBcSdc[it]);
    HF1(312, event.chisqrBcSdc[it]);
    HF1(314, event.x0BcSdc[it]);
    HF1(315, event.y0BcSdc[it]);
    HF1(316, event.u0BcSdc[it]);
    HF1(317, event.v0BcSdc[it]);
    HF2(318, event.utgtBcSdc[it], event.xtgtBcSdc[it]);
    HF2(319, event.vtgtBcSdc[it], event.ytgtBcSdc[it]);
    HF2(320, event.xtgtBcSdc[it], event.ytgtBcSdc[it]);
    HF1(321, event.xtgtBcSdc[it]);
    HF1(322, event.ytgtBcSdc[it]);
    HF1(323, event.utgtBcSdc[it]);
    HF1(324, event.vtgtBcSdc[it]);
  }

  HF2(401, event.x0Tpc[0], xyuvBcOut_HScoor[0]);
  HF2(402, event.y0Tpc[0], xyuvBcOut_HScoor[1]);
  HF2(403, event.u0Tpc[0], xyuvBcOut_HScoor[2]);
  HF2(404, event.v0Tpc[0], xyuvBcOut_HScoor[3]);
  HF1(405, xyuvBcOut_HScoor[0]-event.x0Tpc[0]);
  HF1(406, xyuvBcOut_HScoor[1]-event.y0Tpc[0]);
  HF1(407, xyuvBcOut_HScoor[2]-event.u0Tpc[0]);
  HF1(408, xyuvBcOut_HScoor[3]-event.v0Tpc[0]);
  HF2(409, xyuvBcOut_HScoor[0], xyuvBcOut_HScoor[0]-event.x0Tpc[0]);
  HF2(410, xyuvBcOut_HScoor[1], xyuvBcOut_HScoor[1]-event.y0Tpc[0]);
  HF2(411, xyuvBcOut_HScoor[0], xyuvBcOut_HScoor[2]-event.u0Tpc[0]);
  HF2(412, xyuvBcOut_HScoor[1], xyuvBcOut_HScoor[3]-event.v0Tpc[0]);

  HF2(501, event.x0Tpc[0], xyuvSdcIn_HScoor[0]);
  HF2(502, event.y0Tpc[0], xyuvSdcIn_HScoor[1]);
  HF2(503, event.u0Tpc[0], xyuvSdcIn_HScoor[2]);
  HF2(504, event.v0Tpc[0], xyuvSdcIn_HScoor[3]);
  HF1(505, xyuvSdcIn_HScoor[0]-event.x0Tpc[0]);
  HF1(506, xyuvSdcIn_HScoor[1]-event.y0Tpc[0]);
  HF1(507, xyuvSdcIn_HScoor[2]-event.u0Tpc[0]);
  HF1(508, xyuvSdcIn_HScoor[3]-event.v0Tpc[0]);
  HF2(509, xyuvSdcIn_HScoor[0], xyuvSdcIn_HScoor[0]-event.x0Tpc[0]);
  HF2(510, xyuvSdcIn_HScoor[1], xyuvSdcIn_HScoor[1]-event.y0Tpc[0]);
  HF2(511, xyuvSdcIn_HScoor[0], xyuvSdcIn_HScoor[2]-event.u0Tpc[0]);
  HF2(512, xyuvSdcIn_HScoor[1], xyuvSdcIn_HScoor[3]-event.v0Tpc[0]);

  HF2(601, xyuvBcOut_HScoor[0], xyuvSdcIn_HScoor[0]);
  HF2(602, xyuvBcOut_HScoor[1], xyuvSdcIn_HScoor[1]);
  HF2(603, xyuvBcOut_HScoor[2], xyuvSdcIn_HScoor[2]);
  HF2(604, xyuvBcOut_HScoor[3], xyuvSdcIn_HScoor[3]);
  HF1(605, xyuvSdcIn_HScoor[0]-xyuvBcOut_HScoor[0]);
  HF1(606, xyuvSdcIn_HScoor[1]-xyuvBcOut_HScoor[1]);
  HF1(607, xyuvSdcIn_HScoor[2]-xyuvBcOut_HScoor[2]);
  HF1(608, xyuvSdcIn_HScoor[3]-xyuvBcOut_HScoor[3]);
  HF2(609, xyuvSdcIn_HScoor[0], xyuvSdcIn_HScoor[0]-xyuvBcOut_HScoor[0]);
  HF2(610, xyuvSdcIn_HScoor[1], xyuvSdcIn_HScoor[1]-xyuvBcOut_HScoor[1]);
  HF2(611, xyuvSdcIn_HScoor[0], xyuvSdcIn_HScoor[2]-xyuvBcOut_HScoor[2]);
  HF2(612, xyuvSdcIn_HScoor[1], xyuvSdcIn_HScoor[3]-xyuvBcOut_HScoor[3]);

  if(event.ntBcSdc==1){

    std::vector<Double_t> xyuvBcSdc_HScoor(4);
    {
      //On the HS coordinate
      Double_t zpos = Global2LocalPos(xytgtGlobal, posV0, 0., 0., 0.).z();
      TVector3 localpos(src.u0BcSdc[0]*zpos + src.x0BcSdc[0], src.v0BcSdc[0]*zpos + src.y0BcSdc[0], zpos);
      TVector3 globalpos = Local2GlobalPos(localpos, posV0, 0., 0., 0.);
      xyuvBcSdc_HScoor[0] = Global2LocalPos(globalpos, posHS, tilt, RA1, RA2).x();
      xyuvBcSdc_HScoor[1] = Global2LocalPos(globalpos, posHS, tilt, RA1, RA2).y();
      xyuvBcSdc_HScoor[2] = src.u0BcSdc[0] - TMath::Tan(RA2*TMath::DegToRad());
      xyuvBcSdc_HScoor[3] = src.v0BcSdc[0] + TMath::Tan(RA1*TMath::DegToRad());
    }

    HF2(701, event.xtgtBcOut[0], event.xtgtBcSdc[0]);
    HF2(702, event.ytgtBcOut[0], event.ytgtBcSdc[0]);
    HF2(703, event.utgtBcOut[0], event.utgtBcSdc[0]);
    HF2(704, event.vtgtBcOut[0], event.vtgtBcSdc[0]);
    HF1(705, event.xtgtBcSdc[0]-event.xtgtBcOut[0]);
    HF1(706, event.ytgtBcSdc[0]-event.ytgtBcOut[0]);
    HF1(707, event.utgtBcSdc[0]-event.utgtBcOut[0]);
    HF1(708, event.vtgtBcSdc[0]-event.vtgtBcOut[0]);
    HF2(709, event.xtgtBcSdc[0], event.xtgtBcSdc[0]-event.xtgtBcOut[0]);
    HF2(710, event.ytgtBcSdc[0], event.ytgtBcSdc[0]-event.ytgtBcOut[0]);
    HF2(711, event.xtgtBcSdc[0], event.utgtBcSdc[0]-event.utgtBcOut[0]);
    HF2(712, event.ytgtBcSdc[0], event.vtgtBcSdc[0]-event.vtgtBcOut[0]);

    HF2(721, event.x0BcOut[0], event.x0BcSdc[0]);
    HF2(722, event.y0BcOut[0], event.y0BcSdc[0]);
    HF2(723, event.u0BcOut[0], event.u0BcSdc[0]);
    HF2(724, event.v0BcOut[0], event.v0BcSdc[0]);
    HF1(725, event.x0BcSdc[0]-event.x0BcOut[0]);
    HF1(726, event.y0BcSdc[0]-event.y0BcOut[0]);
    HF1(727, event.u0BcSdc[0]-event.u0BcOut[0]);
    HF1(728, event.v0BcSdc[0]-event.v0BcOut[0]);
    HF2(729, event.x0BcSdc[0], event.x0BcSdc[0]-event.x0BcOut[0]);
    HF2(730, event.y0BcSdc[0], event.y0BcSdc[0]-event.y0BcOut[0]);
    HF2(731, event.x0BcSdc[0], event.u0BcSdc[0]-event.u0BcOut[0]);
    HF2(732, event.y0BcSdc[0], event.v0BcSdc[0]-event.v0BcOut[0]);

    HF2(801, xyuvSdcIn_HScoor[0], xyuvBcSdc_HScoor[0]);
    HF2(802, xyuvSdcIn_HScoor[1], xyuvBcSdc_HScoor[1]);
    HF2(803, xyuvSdcIn_HScoor[2], xyuvBcSdc_HScoor[2]);
    HF2(804, xyuvSdcIn_HScoor[3], xyuvBcSdc_HScoor[3]);
    HF1(805, xyuvBcSdc_HScoor[0]-xyuvSdcIn_HScoor[0]);
    HF1(806, xyuvBcSdc_HScoor[1]-xyuvSdcIn_HScoor[1]);
    HF1(807, xyuvBcSdc_HScoor[2]-xyuvSdcIn_HScoor[2]);
    HF1(808, xyuvBcSdc_HScoor[3]-xyuvSdcIn_HScoor[3]);
    HF2(809, xyuvBcSdc_HScoor[0], xyuvBcSdc_HScoor[0]-xyuvSdcIn_HScoor[0]);
    HF2(810, xyuvBcSdc_HScoor[1], xyuvBcSdc_HScoor[1]-xyuvSdcIn_HScoor[1]);
    HF2(811, xyuvBcSdc_HScoor[0], xyuvBcSdc_HScoor[2]-xyuvSdcIn_HScoor[2]);
    HF2(812, xyuvBcSdc_HScoor[1], xyuvBcSdc_HScoor[3]-xyuvSdcIn_HScoor[3]);

    HF2(901, event.x0Tpc[0], xyuvBcSdc_HScoor[0]);
    HF2(902, event.y0Tpc[0], xyuvBcSdc_HScoor[1]);
    HF2(903, event.u0Tpc[0], xyuvBcSdc_HScoor[2]);
    HF2(904, event.v0Tpc[0], xyuvBcSdc_HScoor[3]);
    HF1(905, xyuvBcSdc_HScoor[0]-event.x0Tpc[0]);
    HF1(906, xyuvBcSdc_HScoor[1]-event.y0Tpc[0]);
    HF1(907, xyuvBcSdc_HScoor[2]-event.u0Tpc[0]);
    HF1(908, xyuvBcSdc_HScoor[3]-event.v0Tpc[0]);
    HF2(909, xyuvBcSdc_HScoor[0], xyuvBcSdc_HScoor[0]-event.x0Tpc[0]);
    HF2(910, xyuvBcSdc_HScoor[1], xyuvBcSdc_HScoor[1]-event.y0Tpc[0]);
    HF2(911, xyuvBcSdc_HScoor[0], xyuvBcSdc_HScoor[2]-event.u0Tpc[0]);
    HF2(912, xyuvBcSdc_HScoor[1], xyuvBcSdc_HScoor[3]-event.v0Tpc[0]);

    //Double_t BCSdcZoffset = - zK18Target + zTarget;
    HF2(821, xyuvSdcIn_HScoor[0], xyuvBcSdc_HScoor[0]);
    HF2(822, xyuvSdcIn_HScoor[1], xyuvBcSdc_HScoor[1]);
    HF2(823, xyuvSdcIn_HScoor[2], xyuvBcSdc_HScoor[2]);
    HF2(824, xyuvSdcIn_HScoor[3], xyuvBcSdc_HScoor[3]);
    HF1(825, xyuvBcSdc_HScoor[0] - xyuvSdcIn_HScoor[0]);
    HF1(826, xyuvBcSdc_HScoor[1] - xyuvSdcIn_HScoor[1]);
    HF1(827, xyuvBcSdc_HScoor[2] - xyuvSdcIn_HScoor[2]);
    HF1(828, xyuvBcSdc_HScoor[3] - xyuvSdcIn_HScoor[3]);
    HF2(829, xyuvBcSdc_HScoor[0], xyuvBcSdc_HScoor[0] - xyuvSdcIn_HScoor[0]);
    HF2(830, xyuvBcSdc_HScoor[1], xyuvBcSdc_HScoor[1] - xyuvSdcIn_HScoor[1]);
    HF2(831, xyuvBcSdc_HScoor[0], xyuvBcSdc_HScoor[2] - xyuvSdcIn_HScoor[2]);
    HF2(832, xyuvBcSdc_HScoor[1], xyuvBcSdc_HScoor[3] - xyuvSdcIn_HScoor[3]);

    for(Int_t it=0; it<event.ntTpc; ++it){
      for(Int_t ih=0; ih<event.nhtrack[it]; ++ih){
	if(event.resolution_x[it][ih] > 0.9e5 && event.resolution_y[it][ih] > 0.9e5 && event.resolution_x[it][ih] > 0.9e5) continue; // exclude dummy hits
	ThreeVector lpos_tpc(event.hitpos_x[it][ih],  event.hitpos_y[it][ih],  event.hitpos_z[it][ih]);
	ThreeVector gpos = gGeom.Local2GlobalPos(IdHS, lpos_tpc);
	ThreeVector lpos_bc = gpos - posV0;
	HF1(1000*(event.hitlayer[it][ih]+1) + 205, -(src.u0BcSdc[it]*lpos_bc.z() + src.x0BcSdc[it]) + lpos_bc.x());
	HF1(1000*(event.hitlayer[it][ih]+1) + 206, -(src.v0BcSdc[it]*lpos_bc.z() + src.y0BcSdc[it]) + lpos_bc.y());
      }
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstClose()
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  const Int_t n = TFileCont.size();
  for(Int_t i=0; i<n; ++i){
    if(TTreeReaderCont[i]) delete TTreeReaderCont[i];
    if(TTreeCont[i]) delete TTreeCont[i];
    if(TFileCont[i]) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  HB1(1, "Status", 21, 0., 21.);
  HB1(10, "#Tracks TPC", 40, 0., 40.);
  HB1(11, "#Hits of Track TPC", 50, 0., 50.);
  HB1(12, "Chisqr TPC", 500, 0., 100.);
  HB1(14, "X0 TPC", 600, -300., 300.);
  HB1(15, "Y0 TPC", 200, -100., 100.);
  HB1(16, "U0 TPC", 200, -0.20, 0.20);
  HB1(17, "V0 TPC", 200, -0.20, 0.20);
  HB2(18, "U0%X0 TPC", 300, -300., 300., 100, -0.20, 0.20);
  HB2(19, "V0%Y0 TPC", 100, -100., 100., 100, -0.20, 0.20);
  HB2(20, "Y0%X0 TPC", 300, -300., 300., 100, -100, 100);
  HB2(21, "atan(U0)%X0 TPC", 300, -300., 300., 100, -20, 20);
  HB2(22, "atan(V0)%Y0 TPC", 300, -300., 300., 100, -20, 20);

  HB1(110, "#Tracks BcOut", 40, 0., 40.);
  HB1(111, "#Hits of Track BcOut", 50, 0., 50.);
  HB1(112, "Chisqr BcOut", 500, 0., 500.);
  HB1(114, "X0 BcOut", 400, -100., 100.);
  HB1(115, "Y0 BcOut", 400, -100., 100.);
  HB1(116, "U0 BcOut", 200, -0.20, 0.20);
  HB1(117, "V0 BcOut", 200, -0.20, 0.20);
  HB2(118, "U%Xtgt BcOut", 100, -0.20, 0.20, 100, -100., 100.);
  HB2(119, "V%Ytgt BcOut", 100, -0.20, 0.20, 100, -100., 100.);
  HB2(120, "X%Ytgt BcOut", 100, -100., 100., 100, -100, 100);
  HB1(121, "Xtgt BcOut", 400, -100., 100.);
  HB1(122, "Ytgt BcOut", 400, -100., 100.);
  HB1(123, "Utgt BcOut", 200, -0.20, 0.20);
  HB1(124, "Vtgt BcOut", 200, -0.20, 0.20);

  HB1(210, "#Tracks SdcIn", 40, 0., 40.);
  HB1(211, "#Hits of Track SdcIn", 50, 0., 50.);
  HB1(212, "Chisqr SdcIn", 500, 0., 500.);
  HB1(214, "X0 SdcIn", 400, -100., 100.);
  HB1(215, "Y0 SdcIn", 400, -100., 100.);
  HB1(216, "U0 SdcIn", 200, -0.20, 0.20);
  HB1(217, "V0 SdcIn", 200, -0.20, 0.20);
  HB2(218, "U%Xtgt SdcIn", 100, -0.20, 0.20, 100, -100., 100.);
  HB2(219, "V%Ytgt SdcIn", 100, -0.20, 0.20, 100, -100., 100.);
  HB2(220, "X%Ytgt SdcIn", 100, -100., 100., 100, -100, 100);
  HB1(221, "Xtgt SdcIn", 400, -100., 100.);
  HB1(222, "Ytgt SdcIn", 400, -100., 100.);
  HB1(223, "Utgt SdcIn", 200, -0.20, 0.20);
  HB1(224, "Vtgt SdcIn", 200, -0.20, 0.20);

  HB1(310, "#Tracks BcSdc", 40, 0., 40.);
  HB1(311, "#Hits of Track BcSdc", 50, 0., 50.);
  HB1(312, "Chisqr BcSdc", 500, 0., 500.);
  HB1(314, "X0 BcSdc", 400, -100., 100.);
  HB1(315, "Y0 BcSdc", 400, -100., 100.);
  HB1(316, "U0 BcSdc", 200, -0.20, 0.20);
  HB1(317, "V0 BcSdc", 200, -0.20, 0.20);
  HB2(318, "U%Xtgt BcSdc", 100, -0.20, 0.20, 100, -100., 100.);
  HB2(319, "V%Ytgt BcSdc", 100, -0.20, 0.20, 100, -100., 100.);
  HB2(320, "X%Ytgt BcSdc", 100, -100., 100., 100, -100, 100);
  HB1(321, "Xtgt BcSdc", 400, -100., 100.);
  HB1(322, "Ytgt BcSdc", 400, -100., 100.);
  HB1(323, "Utgt BcSdc", 200, -0.20, 0.20);
  HB1(324, "Vtgt BcSdc", 200, -0.20, 0.20);

  HB2(401, "Xtgt BcOut%%Tpc", 400, -200., 200., 400, -200., 200.);
  HB2(402, "Ytgt BcOut%%Tpc", 400, -100., 100., 400, -100., 100.);
  HB2(403, "Utgt BcOut%%Tpc", 400, -0.15, 0.15, 400, -0.15, 0.15);
  HB2(404, "Vtgt BcOut%%Tpc", 400, -0.05, 0.05, 400, -0.05, 0.05);
  HB1(405, "Xtgt, BcOut-Tpc", 640, -16., 16.);
  HB1(406, "Ytgt, BcOut-Tpc", 640, -16., 16.);
  HB1(407, "Utgt, BcOut-Tpc", 200, -0.05, 0.05);
  HB1(408, "Vtgt, BcOut-Tpc", 200, -0.05, 0.05);
  HB2(409, "Xtgt BcOut-Tpc%%Xtgt BcOut", 400, -150., 150., 400, -10., 10.);
  HB2(410, "Ytgt BcOut-Tpc%%Ytgt BcOut", 400, -100., 100., 400, -15., 15.);
  HB2(411, "Utgt BcOut-Tpc%%Xtgt BcOut", 400, -150., 150., 400, -0.03, 0.03);
  HB2(412, "Vtgt BcOut-Tpc%%Ytgt BcOut", 400, -100., 100., 400, -0.03, 0.03);

  HB2(501, "Xtgt SdcIn%%Tpc", 400, -200., 200., 400, -200., 200.);
  HB2(502, "Ytgt SdcIn%%Tpc", 400, -100., 100., 400, -100., 100.);
  HB2(503, "Utgt SdcIn%%Tpc", 400, -0.15, 0.15, 400, -0.15, 0.15);
  HB2(504, "Vtgt SdcIn%%Tpc", 400, -0.05, 0.05, 400, -0.05, 0.05);
  HB1(505, "Xtgt, SdcIn-Tpc", 640, -16., 16.);
  HB1(506, "Ytgt, SdcIn-Tpc", 640, -16., 16.);
  HB1(507, "Utgt, SdcIn-Tpc", 200, -0.05, 0.05);
  HB1(508, "Vtgt, SdcIn-Tpc", 200, -0.05, 0.05);
  HB2(509, "Xtgt SdcIn-Tpc%%Xtgt SdcIn", 400, -150., 150., 400, -10., 10.);
  HB2(510, "Ytgt SdcIn-Tpc%%Ytgt SdcIn", 400, -100., 100., 400, -15., 15.);
  HB2(511, "Utgt SdcIn-Tpc%%Xtgt SdcIn", 400, -150., 150., 400, -0.03, 0.03);
  HB2(512, "Vtgt SdcIn-Tpc%%Ytgt SdcIn", 400, -100., 100., 400, -0.03, 0.03);

  HB2(601, "Xtgt SdcIn%%BcOut", 400, -200., 200., 400, -200., 200.);
  HB2(602, "Ytgt SdcIn%%BcOut", 400, -100., 100., 400, -100., 100.);
  HB2(603, "Utgt SdcIn%%BcOut", 400, -0.15, 0.15, 400, -0.15, 0.15);
  HB2(604, "Vtgt SdcIn%%BcOut", 400, -0.05, 0.05, 400, -0.05, 0.05);
  HB1(605, "Xtgt, SdcIn-BcOut", 640, -16., 16.);
  HB1(606, "Ytgt, SdcIn-BcOut", 640, -16., 16.);
  HB1(607, "Utgt, SdcIn-BcOut", 200, -0.05, 0.05);
  HB1(608, "Vtgt, SdcIn-BcOut", 200, -0.05, 0.05);
  HB2(609, "Xtgt SdcIn-BcOut%%Xtgt SdcIn", 400, -150., 150., 400, -10., 10.);
  HB2(610, "Ytgt SdcIn-BcOut%%Ytgt SdcIn", 400, -100., 100., 400, -15., 15.);
  HB2(611, "Utgt SdcIn-BcOut%%Xtgt SdcIn", 400, -150., 150., 400, -0.03, 0.03);
  HB2(612, "Vtgt SdcIn-BcOut%%Ytgt SdcIn", 400, -100., 100., 400, -0.03, 0.03);

  HB2(701, "Xtgt BcSdc%%BcOut", 400, -200., 200., 400, -200., 200.);
  HB2(702, "Ytgt BcSdc%%BcOut", 400, -100., 100., 400, -100., 100.);
  HB2(703, "Utgt BcSdc%%BcOut", 400, -0.15, 0.15, 400, -0.15, 0.15);
  HB2(704, "Vtgt BcSdc%%BcOut", 400, -0.05, 0.05, 400, -0.05, 0.05);
  HB1(705, "Xtgt, BcSdc-BcOut", 640, -16., 16.);
  HB1(706, "Ytgt, BcSdc-BcOut", 640, -16., 16.);
  HB1(707, "Utgt, BcSdc-BcOut", 200, -0.05, 0.05);
  HB1(708, "Vtgt, BcSdc-BcOut", 200, -0.05, 0.05);
  HB2(709, "Xtgt BcSdc-BcOut%%Xtgt BcSdc",  400, -150., 150., 400, -10., 10.);
  HB2(710, "Ytgt BcSdc-BcOut%%Ytgt BcSdc", 400, -100., 100., 400, -15., 15.);
  HB2(711, "Utgt BcSdc-BcOut%%Xtgt BcSdc", 400, -150., 150., 400, -0.03, 0.03);
  HB2(712, "Vtgt BcSdc-BcOut%%Ytgt BcSdc", 400, -100., 100., 400, -0.03, 0.03);
  HB2(721, "X0 BcSdc%%BcOut", 400, -200., 200., 400, -200., 200.);
  HB2(722, "Y0 BcSdc%%BcOut", 400, -100., 100., 400, -100., 100.);
  HB2(723, "U0 BcSdc%%BcOut", 400, -0.15, 0.15, 400, -0.15, 0.15);
  HB2(724, "V0 BcSdc%%BcOut", 400, -0.05, 0.05, 400, -0.05, 0.05);
  HB1(725, "X0, BcSdc-BcOut", 640, -16., 16.);
  HB1(726, "Y0, BcSdc-BcOut", 640, -16., 16.);
  HB1(727, "U0, BcSdc-BcOut", 200, -0.05, 0.05);
  HB1(728, "V0, BcSdc-BcOut", 200, -0.05, 0.05);
  HB2(729, "X0 BcSdc-BcOut%%X0 BcSdc",  400, -150., 150., 400, -10., 10.);
  HB2(730, "Y0 BcSdc-BcOut%%Y0 BcSdc", 400, -100., 100., 400, -15., 15.);
  HB2(731, "U0 BcSdc-BcOut%%X0 BcSdc", 400, -150., 150., 400, -0.03, 0.03);
  HB2(732, "V0 BcSdc-BcOut%%Y0 BcSdc", 400, -100., 100., 400, -0.03, 0.03);

  HB2(801, "Xtgt BcSdc%%SdcIn", 400, -200., 200., 400, -200., 200.);
  HB2(802, "Ytgt BcSdc%%SdcIn", 400, -100., 100., 400, -100., 100.);
  HB2(803, "Utgt BcSdc%%SdcIn", 400, -0.15, 0.15, 400, -0.15, 0.15);
  HB2(804, "Vtgt BcSdc%%SdcIn", 400, -0.05, 0.05, 400, -0.05, 0.05);
  HB1(805, "Xtgt, BcSdc-SdcIn", 640, -16., 16.);
  HB1(806, "Ytgt, BcSdc-SdcIn", 640, -16., 16.);
  HB1(807, "Utgt, BcSdc-SdcIn", 200, -0.05, 0.05);
  HB1(808, "Vtgt, BcSdc-SdcIn", 200, -0.05, 0.05);
  HB2(809, "Xtgt BcSdc-SdcIn%%Xtgt BcSdc", 400, -150., 150., 400, -10., 10.);
  HB2(810, "Ytgt BcSdc-SdcIn%%Ytgt BcSdc", 400, -100., 100., 400, -15., 15.);
  HB2(811, "Utgt BcSdc-SdcIn%%Xtgt BcSdc", 400, -150., 150., 400, -0.03, 0.03);
  HB2(812, "Vtgt BcSdc-SdcIn%%Ytgt BcSdc", 400, -100., 100., 400, -0.03, 0.03);
  HB2(821, "X0 BcSdc%%SdcIn", 400, -200., 200., 400, -200., 200.);
  HB2(822, "Y0 BcSdc%%SdcIn", 400, -100., 100., 400, -100., 100.);
  HB2(823, "U0 BcSdc%%SdcIn", 400, -0.15, 0.15, 400, -0.15, 0.15);
  HB2(824, "V0 BcSdc%%SdcIn", 400, -0.05, 0.05, 400, -0.05, 0.05);
  HB1(825, "X0, BcSdc-SdcIn", 640, -16., 16.);
  HB1(826, "Y0, BcSdc-SdcIn", 640, -16., 16.);
  HB1(827, "U0, BcSdc-SdcIn", 200, -0.05, 0.05);
  HB1(828, "V0, BcSdc-SdcIn", 200, -0.05, 0.05);
  HB2(829, "X0 BcSdc-SdcIn%%X0 BcSdc", 400, -150., 150., 400, -10., 10.);
  HB2(830, "Y0 BcSdc-SdcIn%%Y0 BcSdc", 400, -100., 100., 400, -15., 15.);
  HB2(831, "U0 BcSdc-SdcIn%%X0 BcSdc", 400, -150., 150., 400, -0.03, 0.03);
  HB2(832, "V0 BcSdc-SdcIn%%Y0 BcSdc", 400, -100., 100., 400, -0.03, 0.03);

  HB2(901, "Xtgt BcSdc%%Tpc", 400, -200., 200., 400, -200., 200.);
  HB2(902, "Ytgt BcSdc%%Tpc", 400, -100., 100., 400, -100., 100.);
  HB2(903, "Utgt BcSdc%%Tpc", 400, -0.15, 0.15, 400, -0.15, 0.15);
  HB2(904, "Vtgt BcSdc%%Tpc", 400, -0.05, 0.05, 400, -0.05, 0.05);
  HB1(905, "Xtgt, BcSdc-Tpc", 640, -16., 16.);
  HB1(906, "Ytgt, BcSdc-Tpc", 640, -16., 16.);
  HB1(907, "Utgt, BcSdc-Tpc", 200, -0.05, 0.05);
  HB1(908, "Vtgt, BcSdc-Tpc", 200, -0.05, 0.05);
  HB2(909, "Xtgt BcSdc-Tpc%%Xtgt BcSdc", 400, -150., 150., 400, -10., 10.);
  HB2(910, "Ytgt BcSdc-Tpc%%Ytgt BcSdc", 400, -100., 100., 400, -15., 15.);
  HB2(911, "Utgt BcSdc-Tpc%%Xtgt BcSdc", 400, -150., 150., 400, -0.03, 0.03);
  HB2(912, "Vtgt BcSdc-Tpc%%Ytgt BcSdc", 400, -100., 100., 400, -0.03, 0.03);

  for(Int_t i=0; i<NumOfLayersTPC; ++i){
    HB2(1000*(i+1) + 1, Form("TPC L%d X Residual%%SdcIn X", i+1), 400, -200., 200., 400, -20., 20.);
    HB2(1000*(i+1) + 2, Form("TPC L%d Y Residual%%SdcIn Y", i+1), 400, -200., 200., 400, -20., 20.);
    HB2(1000*(i+1) + 3, Form("TPC L%d Z Residual%%SdcIn X", i+1), 400, -200., 200., 400, -20., 20.);
    HB2(1000*(i+1) + 4, Form("TPC L%d Local X Residual%%SdcIn X", i+1), 800, -200., 200., 400, -20., 20.);
    HB2(1000*(i+1) + 5, Form("TPC L%d Local Y Residual%%SdcIn X", i+1), 800, -200., 200., 400, -20., 20.);
    HB2(1000*(i+1) + 6, Form("TPC L%d XZ Residual%%SdcIn X", i+1), 400, -200., 200., 200, 0., 20.);

    HB1(1000*(i+1) + 11, Form("TPC L%d X Pull", i+1), 200, -5., 5.);
    HB1(1000*(i+1) + 12, Form("TPC L%d Y Pull", i+1), 200, -5., 5.);
    HB1(1000*(i+1) + 13, Form("TPC L%d Z Pull", i+1), 200, -5., 5.);
    HB1(1000*(i+1) + 14, Form("TPC L%d Local X Pull", i+1), 200, -5., 5.);
    HB1(1000*(i+1) + 15, Form("TPC L%d Local Y Pull", i+1), 200, -5., 5.);

    HB1(1000*(i+1) + 101, Form("TPC L%d X Residual", i+1), 200, -8., 8.);
    HB1(1000*(i+1) + 102, Form("TPC L%d Y Residual", i+1), 200, -8., 8.);
    HB1(1000*(i+1) + 103, Form("TPC L%d Z Residual", i+1), 200, -8., 8.);
    HB1(1000*(i+1) + 104, Form("TPC L%d Local X Residual", i+1), 400, -16., 16.);
    HB1(1000*(i+1) + 105, Form("TPC L%d Local Y Residual", i+1), 200, -8., 8.);
    HB1(1000*(i+1) + 106, Form("TPC L%d XZ Residual", i+1), 100, 0, 8.);

    HB1(1000*(i+1) + 201, Form("TPC L%d, BcOut X Residual", i+1), 200, -8., 8.);
    HB1(1000*(i+1) + 202, Form("TPC L%d, BcOut Y Residual", i+1), 200, -8., 8.);
    HB1(1000*(i+1) + 203, Form("TPC L%d, SdcIn X Residual", i+1), 200, -8., 8.);
    HB1(1000*(i+1) + 204, Form("TPC L%d, SdcIn Y Residual", i+1), 200, -8., 8.);
    HB1(1000*(i+1) + 205, Form("TPC L%d, BcSdc X Residual", i+1), 200, -8., 8.);
    HB1(1000*(i+1) + 206, Form("TPC L%d, BcSdc Y Residual", i+1), 200, -8., 8.);
  }
  HB2(TPCCalibHid, "Layer%ResY", 32, 0, 32, 400, -5, 5);

  HBTree("tpc", "tree of DstTPCBcOutSdcInTracking");
  tree->Branch( "status", &event.status );
  tree->Branch( "runnum", &event.runnum );
  tree->Branch( "evnum", &event.evnum );
  tree->Branch( "trigpat", &event.trigpat );
  tree->Branch( "trigflag", &event.trigflag );

  tree->Branch("ntTpc", &event.ntTpc);
  tree->Branch("nhtrack", &event.nhtrack);
  tree->Branch("chisqrTpc", &event.chisqrTpc);
  tree->Branch("x0Tpc", &event.x0Tpc);
  tree->Branch("y0Tpc", &event.y0Tpc);
  tree->Branch("u0Tpc", &event.u0Tpc);
  tree->Branch("v0Tpc", &event.v0Tpc);
  tree->Branch("theta", &event.theta);
  tree->Branch("hitlayer", &event.hitlayer);
  tree->Branch("hitpos_x", &event.hitpos_x);
  tree->Branch("hitpos_y", &event.hitpos_y);
  tree->Branch("hitpos_z", &event.hitpos_z);
  tree->Branch("calpos_x", &event.calpos_x);
  tree->Branch("calpos_y", &event.calpos_y);
  tree->Branch("calpos_z", &event.calpos_z);
  tree->Branch("residual", &event.residual);
  tree->Branch("residual_x", &event.residual_x);
  tree->Branch("residual_y", &event.residual_y);
  tree->Branch("residual_z", &event.residual_z);
  tree->Branch("residual_horizontal", &event.residual_horizontal);
  tree->Branch("residual_vertical", &event.residual_vertical);
  tree->Branch("resolution_x", &event.resolution_x);
  tree->Branch("resolution_y", &event.resolution_y);
  tree->Branch("resolution_z", &event.resolution_z);
  tree->Branch("resolution_horizontal", &event.resolution_horizontal);
  tree->Branch("resolution_vertical", &event.resolution_vertical);
  tree->Branch("dE", &event.dE );
  tree->Branch("dEdx", &event.dEdx );
  tree->Branch("pathhit", &event.pathhit);
  tree->Branch("theta_diff", &event.theta_diff);
  tree->Branch("track_cluster_de", &event.track_cluster_de);
  tree->Branch("track_cluster_size", &event.track_cluster_size);
  tree->Branch("track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch("track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch("track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch("track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch("track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch("track_cluster_row_center", &event.track_cluster_row_center);

#if Exclusive
  tree->Branch("exresidual", &event.exresidual);
  tree->Branch("exresidual_x", &event.exresidual_x);
  tree->Branch("exresidual_y", &event.exresidual_y);
  tree->Branch("exresidual_z", &event.exresidual_z);
  tree->Branch("exresidual_horizontal", &event.exresidual_horizontal);
  tree->Branch("exresidual_vertical", &event.exresidual_vertical);
#endif

  tree->Branch("btof", &event.btof);
  tree->Branch("ntBcSdc", &event.ntBcSdc);
  tree->Branch("nhBcSdc", &event.nhBcSdc);
  tree->Branch("chisqrBcSdc", &event.chisqrBcSdc);
  tree->Branch("x0BcSdc", &event.x0BcSdc);
  tree->Branch("y0BcSdc", &event.y0BcSdc);
  tree->Branch("u0BcSdc", &event.u0BcSdc);
  tree->Branch("v0BcSdc", &event.v0BcSdc);
  tree->Branch("xtgtBcSdc", &event.xtgtBcSdc);
  tree->Branch("ytgtBcSdc", &event.ytgtBcSdc);
  tree->Branch("utgtBcSdc", &event.utgtBcSdc);
  tree->Branch("vtgtBcSdc", &event.vtgtBcSdc);

  tree->Branch("ntBcOut", &event.ntBcOut);
  tree->Branch("nhBcOut", &event.nhBcOut);
  tree->Branch("chisqrBcOut", &event.chisqrBcOut);
  tree->Branch("x0BcOut", &event.x0BcOut);
  tree->Branch("y0BcOut", &event.y0BcOut);
  tree->Branch("u0BcOut", &event.u0BcOut);
  tree->Branch("v0BcOut", &event.v0BcOut);
  tree->Branch("xtgtBcOut", &event.xtgtBcOut);
  tree->Branch("ytgtBcOut", &event.ytgtBcOut);
  tree->Branch("utgtBcOut", &event.utgtBcOut);
  tree->Branch("vtgtBcOut", &event.vtgtBcOut);

  tree->Branch("ntSdcIn", &event.ntSdcIn);
  tree->Branch("nhSdcIn", &event.nhSdcIn);
  tree->Branch("chisqrSdcIn", &event.chisqrSdcIn);
  tree->Branch("x0SdcIn", &event.x0SdcIn);
  tree->Branch("y0SdcIn", &event.y0SdcIn);
  tree->Branch("u0SdcIn", &event.u0SdcIn);
  tree->Branch("v0SdcIn", &event.v0SdcIn);
  tree->Branch("xtgtSdcIn", &event.xtgtSdcIn);
  tree->Branch("ytgtSdcIn", &event.ytgtSdcIn);
  tree->Branch("utgtSdcIn", &event.utgtSdcIn);
  tree->Branch("vtgtSdcIn", &event.vtgtSdcIn);

  TTreeReaderCont[kTpc] = new TTreeReader("tpc", TFileCont[kTpc]);
  const auto& reader = TTreeReaderCont[kTpc];
  src.runnum = new TTreeReaderValue<Int_t>(*reader, "runnum");
  src.evnum = new TTreeReaderValue<Int_t>(*reader, "evnum");
  src.trigpat = new TTreeReaderValue<std::vector<Int_t>>(*reader, "trigpat");
  src.trigflag = new TTreeReaderValue<std::vector<Int_t>>(*reader, "trigflag");

  src.ntTpc = new TTreeReaderValue<Int_t>(*reader, "ntTpc");
  src.nhtrack = new TTreeReaderValue<std::vector<Int_t>>(*reader, "nhtrack");
  src.chisqrTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "chisqrTpc");
  src.x0Tpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "x0Tpc");
  src.y0Tpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "y0Tpc");
  src.u0Tpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "u0Tpc");
  src.v0Tpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "v0Tpc");
  src.theta = new TTreeReaderValue<std::vector<Double_t>>(*reader, "theta");
  src.dE = new TTreeReaderValue<std::vector<Double_t>>(*reader, "dE");
  src.dEdx = new TTreeReaderValue<std::vector<Double_t>>(*reader, "dEdx");
  src.hitlayer = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "hitlayer");
  src.hitpos_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "hitpos_x");
  src.hitpos_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "hitpos_y");
  src.hitpos_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "hitpos_z");
  src.calpos_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "calpos_x");
  src.calpos_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "calpos_y");
  src.calpos_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "calpos_z");
  src.residual = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "residual");
  src.residual_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "residual_x");
  src.residual_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "residual_y");
  src.residual_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "residual_z");
  src.residual_horizontal = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "residual_horizontal");
  src.residual_vertical = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "residual_vertical");
  src.resolution_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "resolution_x");
  src.resolution_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "resolution_y");
  src.resolution_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "resolution_z");
  src.resolution_horizontal = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "resolution_horizontal");
  src.resolution_vertical = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "resolution_vertical");
  src.pathhit = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "pathhit");
  src.theta_diff = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "theta_diff");
  src.track_cluster_de = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "track_cluster_de");
  src.track_cluster_size = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "track_cluster_size");
  src.track_cluster_mrow = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "track_cluster_mrow");
  src.track_cluster_de_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "track_cluster_de_center");
  src.track_cluster_x_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "track_cluster_x_center");
  src.track_cluster_y_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "track_cluster_y_center");
  src.track_cluster_z_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "track_cluster_z_center");
  src.track_cluster_row_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "track_cluster_row_center");

  src.exresidual = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "exresidual");
  src.exresidual_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "exresidual_x");
  src.exresidual_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "exresidual_y");
  src.exresidual_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "exresidual_z");
  src.exresidual_horizontal = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "exresidual_horizontal");
  src.exresidual_vertical = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "exresidual_vertical");

  TTreeCont[kBcSdc]->SetBranchStatus("btof",               1);
  TTreeCont[kBcSdc]->SetBranchStatus("ntrack",             1);
  TTreeCont[kBcSdc]->SetBranchStatus("nh",                 1);
  TTreeCont[kBcSdc]->SetBranchStatus("chisqr",             1);
  TTreeCont[kBcSdc]->SetBranchStatus("x0",                 1);
  TTreeCont[kBcSdc]->SetBranchStatus("y0",                 1);
  TTreeCont[kBcSdc]->SetBranchStatus("u0",                 1);
  TTreeCont[kBcSdc]->SetBranchStatus("v0",                 1);
  TTreeCont[kBcSdc]->SetBranchStatus("ntBcOut",            1);
  TTreeCont[kBcSdc]->SetBranchStatus("nhBcOut",            1);
  TTreeCont[kBcSdc]->SetBranchStatus("chisqrBcOut",        1);
  TTreeCont[kBcSdc]->SetBranchStatus("x0BcOut",            1);
  TTreeCont[kBcSdc]->SetBranchStatus("y0BcOut",            1);
  TTreeCont[kBcSdc]->SetBranchStatus("u0BcOut",            1);
  TTreeCont[kBcSdc]->SetBranchStatus("v0BcOut",            1);
  TTreeCont[kBcSdc]->SetBranchStatus("ntSdcIn",            1);
  TTreeCont[kBcSdc]->SetBranchStatus("ntSdcIn",            1);
  TTreeCont[kBcSdc]->SetBranchStatus("chisqrSdcIn",        1);
  TTreeCont[kBcSdc]->SetBranchStatus("x0SdcIn",            1);
  TTreeCont[kBcSdc]->SetBranchStatus("y0SdcIn",            1);
  TTreeCont[kBcSdc]->SetBranchStatus("u0SdcIn",            1);
  TTreeCont[kBcSdc]->SetBranchStatus("v0SdcIn",            1);

  TTreeCont[kBcSdc]->SetBranchAddress("btof", &src.btof);
  TTreeCont[kBcSdc]->SetBranchAddress("ntrack", &src.ntBcSdc);
  TTreeCont[kBcSdc]->SetBranchAddress("nh", &src.nhBcSdc);
  TTreeCont[kBcSdc]->SetBranchAddress("chisqr", src.chisqrBcSdc);
  TTreeCont[kBcSdc]->SetBranchAddress("x0", src.x0BcSdc);
  TTreeCont[kBcSdc]->SetBranchAddress("y0", src.y0BcSdc);
  TTreeCont[kBcSdc]->SetBranchAddress("u0", src.u0BcSdc);
  TTreeCont[kBcSdc]->SetBranchAddress("v0", src.v0BcSdc);
  TTreeCont[kBcSdc]->SetBranchAddress("ntBcOut", &src.ntBcOut);
  TTreeCont[kBcSdc]->SetBranchAddress("nhBcOut", &src.nhBcOut);
  TTreeCont[kBcSdc]->SetBranchAddress("chisqrBcOut", src.chisqrBcOut);
  TTreeCont[kBcSdc]->SetBranchAddress("x0BcOut", src.x0BcOut);
  TTreeCont[kBcSdc]->SetBranchAddress("y0BcOut", src.y0BcOut);
  TTreeCont[kBcSdc]->SetBranchAddress("u0BcOut", src.u0BcOut);
  TTreeCont[kBcSdc]->SetBranchAddress("v0BcOut", src.v0BcOut);
  TTreeCont[kBcSdc]->SetBranchAddress("ntSdcIn", &src.ntSdcIn);
  TTreeCont[kBcSdc]->SetBranchAddress("nhSdcIn", &src.nhSdcIn);
  TTreeCont[kBcSdc]->SetBranchAddress("chisqrSdcIn", src.chisqrSdcIn);
  TTreeCont[kBcSdc]->SetBranchAddress("x0SdcIn", src.x0SdcIn);
  TTreeCont[kBcSdc]->SetBranchAddress("y0SdcIn", src.y0SdcIn);
  TTreeCont[kBcSdc]->SetBranchAddress("u0SdcIn", src.u0SdcIn);
  TTreeCont[kBcSdc]->SetBranchAddress("v0SdcIn", src.v0SdcIn);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")   &&
     InitializeParameter<TPCParamMan>("TPCPRM") &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

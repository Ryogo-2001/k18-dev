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
#include "TPCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DstHelper.hh"
#include "HodoPHCMan.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCCluster.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrack.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#define Exclusive 1

#define RawHit 0
#define RawCluster 1
#define TrackSearchFailed 1

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto& gUnpacker = GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
const double truncatedMean = 0.8; //80%
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kTpcHit,  kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[TPCHit]",  "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc", "" };
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
  std::vector<Double_t> clkTpc;

  Int_t nhTpc;
  std::vector<Double_t> raw_hitpos_x;
  std::vector<Double_t> raw_hitpos_y;
  std::vector<Double_t> raw_hitpos_z;
  std::vector<Double_t> raw_de;
  std::vector<Int_t> raw_padid;
  std::vector<Int_t> raw_layer;
  std::vector<Int_t> raw_row;

  Int_t nclTpc;
  std::vector<Double_t> cluster_x;
  std::vector<Double_t> cluster_y;
  std::vector<Double_t> cluster_z;
  std::vector<Double_t> cluster_de;
  std::vector<Int_t> cluster_size;
  std::vector<Int_t> cluster_layer;
  std::vector<Double_t> cluster_mrow;
  std::vector<Double_t> cluster_de_center;
  std::vector<Double_t> cluster_x_center;
  std::vector<Double_t> cluster_y_center;
  std::vector<Double_t> cluster_z_center;
  std::vector<Int_t> cluster_row_center;
  std::vector<Int_t> cluster_houghflag;

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

  Int_t ntTpc_inside;
  Double_t prodvtx_x;
  Double_t prodvtx_y;
  Double_t prodvtx_z;

  //failed track container
  Int_t failed_ntTpc;
  std::vector<Int_t> failed_nhtrack;
  std::vector<Double_t> failed_x0Tpc;
  std::vector<Double_t> failed_y0Tpc;
  std::vector<Double_t> failed_u0Tpc;
  std::vector<Double_t> failed_v0Tpc;

  std::vector<std::vector<Double_t>> failed_hitlayer;
  std::vector<std::vector<Double_t>> failed_hitpos_x;
  std::vector<std::vector<Double_t>> failed_hitpos_y;
  std::vector<std::vector<Double_t>> failed_hitpos_z;
  std::vector<std::vector<Double_t>> failed_calpos_x;
  std::vector<std::vector<Double_t>> failed_calpos_y;
  std::vector<std::vector<Double_t>> failed_calpos_z;

  void clear()
  {
    runnum = 0;
    evnum = 0;
    status = 0;
    trigpat.clear();
    trigflag.clear();
    clkTpc.clear();

    nhTpc = 0;
    raw_hitpos_x.clear();
    raw_hitpos_y.clear();
    raw_hitpos_z.clear();
    raw_de.clear();
    raw_padid.clear();
    raw_layer.clear();
    raw_row.clear();

    nclTpc = 0;
    cluster_x.clear();
    cluster_y.clear();
    cluster_z.clear();
    cluster_de.clear();
    cluster_size.clear();
    cluster_layer.clear();
    cluster_mrow.clear();
    cluster_houghflag.clear();
    cluster_de_center.clear();
    cluster_x_center.clear();
    cluster_y_center.clear();
    cluster_z_center.clear();
    cluster_row_center.clear();

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

    ntTpc_inside = 0;
    prodvtx_x = TMath::QuietNaN();
    prodvtx_y = TMath::QuietNaN();
    prodvtx_z = TMath::QuietNaN();

    exresidual.clear();
    exresidual_x.clear();
    exresidual_y.clear();
    exresidual_z.clear();
    exresidual_horizontal.clear();
    exresidual_vertical.clear();

    failed_ntTpc = 0;
    failed_nhtrack.clear();
    failed_x0Tpc.clear();
    failed_y0Tpc.clear();
    failed_u0Tpc.clear();
    failed_v0Tpc.clear();
    failed_hitlayer.clear();
    failed_hitpos_x.clear();
    failed_hitpos_y.clear();
    failed_hitpos_z.clear();
    failed_calpos_x.clear();
    failed_calpos_y.clear();
    failed_calpos_z.clear();
  }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* runnum;
  TTreeReaderValue<Int_t>* evnum;
  TTreeReaderValue<std::vector<Int_t>>* trigpat;
  TTreeReaderValue<std::vector<Int_t>>* trigflag;
  TTreeReaderValue<Int_t>* npadTpc;   // number of pads
  TTreeReaderValue<Int_t>* nhTpc;     // number of hits
  // vector (size=nhTpc)
  TTreeReaderValue<std::vector<Int_t>>* layerTpc;     // layer id
  TTreeReaderValue<std::vector<Int_t>>* rowTpc;       // row id
  TTreeReaderValue<std::vector<Int_t>>* padTpc;       // pad id
  TTreeReaderValue<std::vector<Double_t>>* pedTpc;    // pedestal
  TTreeReaderValue<std::vector<Double_t>>* rmsTpc;    // rms
  TTreeReaderValue<std::vector<Double_t>>* deTpc;     // dE
  TTreeReaderValue<std::vector<Double_t>>* tTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* chisqrTpc; // chi^2 of signal fitting
  TTreeReaderValue<std::vector<Double_t>>* clkTpc;    // clock time
};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
enum eDetHid {
  PadHid    = 100000,
  TPCClHid  = 200000,
};

  Double_t
  TranseverseDistance(Double_t x_center, Double_t z_center, Double_t x, Double_t z)
  {
    Double_t dummy = TMath::Sqrt((x-x_center)*(x-x_center) + (z-z_center)*(z-z_center));
    Double_t dist;
    if(x_center-x<0) dist=-1.*dummy;
    else dist=dummy;
    return dist;
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
  //skip = 309;
  //skip = 209;
  //skip = 215;
  //max_loop = 100;
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
  //if(ievent%1000==0){
  if(ievent%1==0){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;
  event.clkTpc = **src.clkTpc;

  HF1(1, event.status++);

  if(**src.nhTpc == 0)
    return true;

  HF1(1, event.status++);

  if(event.clkTpc.size() != 1){
    std::cerr << "something is wrong" << std::endl;
    return true;
  }

  Double_t clock = event.clkTpc.at(0);
  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, clock);
#if Exclusive
  TPCAna.TrackSearchTPC(true);
#else
  TPCAna.TrackSearchTPC();
#endif

#if RawHit
  Int_t nhTpc = 0;
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    auto hc = TPCAna.GetTPCHC(layer);
    for(const auto& hit : hc){
      if(!hit || !hit->IsGood())
	continue;
      const auto& pos = hit->GetPosition();
      Double_t x = pos.X();
      Double_t y = pos.Y();
      Double_t z = pos.Z();
      Double_t de = hit->GetCDe();
      Int_t pad = hit->GetPad();
      Int_t row = hit->GetRow();
      event.raw_hitpos_x.push_back(x);
      event.raw_hitpos_y.push_back(y);
      event.raw_hitpos_z.push_back(z);
      event.raw_de.push_back(de);
      event.raw_padid.push_back(pad);
      event.raw_layer.push_back(layer);
      event.raw_row.push_back(row);
      ++nhTpc;
    }
  }
  event.nhTpc = nhTpc;
  HF1(1, event.status++);
#endif

#if RawCluster
  Int_t nclTpc = 0;
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    auto hc = TPCAna.GetTPCClCont( layer );
    for( const auto& cl : hc ){
      if( !cl || !cl->IsGood() )
        continue;
      Double_t x = cl->GetX();
      Double_t y = cl->GetY();
      Double_t z = cl->GetZ();
      Double_t de = cl->GetDe();
      Int_t cl_size = cl->GetClusterSize();
      Double_t mrow = cl->MeanRow();
      Int_t houghFlag = cl->GetHoughFlag();
      TPCHit* centerHit = cl->GetCenterHit();
      const TVector3& centerPos = centerHit->GetPosition();
      Double_t centerDe = centerHit->GetCDe();
      Int_t centerRow = centerHit->GetRow();

      event.cluster_x.push_back(x);
      event.cluster_y.push_back(y);
      event.cluster_z.push_back(z);
      event.cluster_de.push_back(de);
      event.cluster_size.push_back(cl_size);
      event.cluster_layer.push_back(layer);
      event.cluster_mrow.push_back(mrow);
      event.cluster_houghflag.push_back(houghFlag);
      event.cluster_de_center.push_back(centerDe);
      event.cluster_x_center.push_back(centerPos.X());
      event.cluster_y_center.push_back(centerPos.Y());
      event.cluster_z_center.push_back(centerPos.Z());
      event.cluster_row_center.push_back(centerRow);

      ++nclTpc;
    }
  }
  event.nclTpc = nclTpc;
  HF1(1, event.status++);
#endif

  Int_t ntTpc = TPCAna.GetNTracksTPC();
  event.ntTpc = ntTpc;

  HF1(10, ntTpc);
  if(event.ntTpc == 0)
     return true;

  HF1(1, event.status++);
  event.nhtrack.resize(ntTpc);
  event.chisqrTpc.resize(ntTpc);
  event.x0Tpc.resize(ntTpc);
  event.y0Tpc.resize(ntTpc);
  event.u0Tpc.resize(ntTpc);
  event.v0Tpc.resize(ntTpc);
  event.theta.resize(ntTpc);
  event.hitlayer.resize(ntTpc);
  event.hitpos_x.resize(ntTpc);
  event.hitpos_y.resize(ntTpc);
  event.hitpos_z.resize(ntTpc);
  event.calpos_x.resize(ntTpc);
  event.calpos_y.resize(ntTpc);
  event.calpos_z.resize(ntTpc);
  event.residual.resize(ntTpc);
  event.residual_x.resize(ntTpc);
  event.residual_y.resize(ntTpc);
  event.residual_z.resize(ntTpc);
  event.residual_horizontal.resize(ntTpc);
  event.residual_vertical.resize(ntTpc);
  event.resolution_x.resize(ntTpc);
  event.resolution_y.resize(ntTpc);
  event.resolution_z.resize(ntTpc);
  event.resolution_horizontal.resize(ntTpc);
  event.resolution_vertical.resize(ntTpc);
  event.exresidual.resize(ntTpc);
  event.exresidual_x.resize(ntTpc);
  event.exresidual_y.resize(ntTpc);
  event.exresidual_z.resize(ntTpc);
  event.exresidual_horizontal.resize(ntTpc);
  event.exresidual_vertical.resize(ntTpc);
  event.dE.resize( ntTpc );
  event.dEdx.resize( ntTpc );
  event.pathhit.resize(ntTpc);
  event.theta_diff.resize(ntTpc);
  event.track_cluster_de.resize(ntTpc);
  event.track_cluster_size.resize(ntTpc);
  event.track_cluster_mrow.resize(ntTpc);
  event.track_cluster_de_center.resize(ntTpc);
  event.track_cluster_x_center.resize(ntTpc);
  event.track_cluster_y_center.resize(ntTpc);
  event.track_cluster_z_center.resize(ntTpc);
  event.track_cluster_row_center.resize(ntTpc);

  Int_t ntrack_intarget = 0;
  Double_t x0_vtx[100] = {0};
  Double_t y0_vtx[100] = {0};
  Double_t u0_vtx[100] = {0};
  Double_t v0_vtx[100] = {0};

  for(Int_t it=0; it<ntTpc; ++it){
    auto track = TPCAna.GetTrackTPC(it);
    if(!track) continue;
    Int_t nh = track->GetNHit();
    Double_t chisqr = track->GetChiSquare();
    Double_t x0=track->GetX0(), y0=track->GetY0();
    Double_t u0=track->GetU0(), v0=track->GetV0();
    Double_t theta = track->GetTheta();
    if(TMath::Abs(x0)<50. && TMath::Abs(y0)<50.){
      x0_vtx[ntrack_intarget] = x0;
      y0_vtx[ntrack_intarget] = y0;
      u0_vtx[ntrack_intarget] = u0;
      v0_vtx[ntrack_intarget] = v0;
      ntrack_intarget++;
    }

    //Tracking information
    Int_t niter = track->GetNIteration();
    Int_t searchtime = track->GetSearchTime();
    Int_t fittime = track->GetFitTime();
    Int_t fitflag = track->GetFitFlag();
    Int_t minuit_status = track->GetMinuitStatus();

    HF1(4, niter);
    HF1(5, fitflag);
    HF1(6, searchtime);
    HF1(7, fittime);
    HF1(8, minuit_status);

    HF1(11, nh);
    HF1(12, chisqr);
    HF1(14, x0);
    HF1(15, y0);
    HF1(16, u0);
    HF1(17, v0);
    HF2(18, x0, u0);
    HF2(19, y0, v0);
    HF2(20, x0, y0);
    event.nhtrack[it] = nh;
    event.chisqrTpc[it] = chisqr;
    event.x0Tpc[it] = x0;
    event.y0Tpc[it] = y0;
    event.u0Tpc[it] = u0;
    event.v0Tpc[it] = v0;
    event.theta[it] = theta;
    event.hitlayer[it].resize(nh);
    event.hitpos_x[it].resize(nh);
    event.hitpos_y[it].resize(nh);
    event.hitpos_z[it].resize(nh);
    event.calpos_x[it].resize(nh);
    event.calpos_y[it].resize(nh);
    event.calpos_z[it].resize(nh);
    event.residual[it].resize(nh);
    event.residual_x[it].resize(nh);
    event.residual_y[it].resize(nh);
    event.residual_z[it].resize(nh);
    event.residual_horizontal[it].resize(nh);
    event.residual_vertical[it].resize(nh);
    event.resolution_x[it].resize(nh);
    event.resolution_y[it].resize(nh);
    event.resolution_z[it].resize(nh);
    event.resolution_horizontal[it].resize(nh);
    event.resolution_vertical[it].resize(nh);
    event.exresidual[it].resize(nh);
    event.exresidual_x[it].resize(nh);
    event.exresidual_y[it].resize(nh);
    event.exresidual_z[it].resize(nh);
    event.exresidual_horizontal[it].resize(nh);
    event.exresidual_vertical[it].resize(nh);
    event.pathhit[it].resize(nh);
    event.theta_diff[it].resize(nh);
    event.track_cluster_de[it].resize(nh);
    event.track_cluster_size[it].resize(nh);
    event.track_cluster_mrow[it].resize(nh);
    event.track_cluster_de_center[it].resize(nh);
    event.track_cluster_x_center[it].resize(nh);
    event.track_cluster_y_center[it].resize(nh);
    event.track_cluster_z_center[it].resize(nh);
    event.track_cluster_row_center[it].resize(nh);

    Double_t de=0.;
    std::vector<Double_t> dEdx_vect;
    for(int ih=0; ih<nh; ++ih){
      auto hit = track->GetHit(ih);
      if(!hit) continue;
      HF1(2, hit->GetHoughDist());
      HF1(3, hit->GetHoughDistY());

      Int_t layer = hit->GetLayer();
      const TVector3& hitpos = hit->GetLocalHitPos();
      const TVector3& calpos = hit->GetLocalCalPos();
      const TVector3& resi_vect = hit->GetResidualVect();
      const TVector3& res_vect = hit->GetResolutionVect();

      TPCHit *clhit = hit -> GetHit();
      TPCCluster *cl = clhit -> GetParentCluster();
      Int_t clsize = cl->GetClusterSize();
      Double_t clde = cl->GetDe();
      Double_t mrow = cl->MeanRow(); // same
      TPCHit* centerHit = cl->GetCenterHit();
      const TVector3& centerPos = centerHit->GetPosition();
      Double_t centerDe = centerHit->GetCDe();
      Int_t centerRow = centerHit->GetRow();
      Double_t residual = hit->GetResidual();

      HF1(TPCClHid, clsize);
      HF1(TPCClHid+(layer+1)*1000, clsize);
      HF1(TPCClHid+1, clde);
      HF1(TPCClHid+(layer+1)*1000+1, clde);
      const TPCHitContainer& hc = cl -> GetHitContainer();
      for(const auto& hits : hc){
	if(!hits || !hits->IsGood()) continue;
	const TVector3& pos = hits->GetPosition();
	Double_t de = hits->GetCDe();
	Double_t transDist = TranseverseDistance(hitpos.x(), hitpos.z(), pos.x(), pos.z());
	Double_t ratio = de/clde;
	HF2(TPCClHid+2, transDist, ratio);
	HF2(TPCClHid+(layer+1)*1000+2, transDist, ratio);
      }
      event.track_cluster_de[it][ih] = clde;
      event.track_cluster_size[it][ih] = clsize;
      event.track_cluster_mrow[it][ih] = mrow;
      event.track_cluster_de_center[it][ih] = centerDe;
      event.track_cluster_x_center[it][ih] = centerPos.X();
      event.track_cluster_y_center[it][ih] = centerPos.Y();
      event.track_cluster_z_center[it][ih] = centerPos.Z();
      event.track_cluster_row_center[it][ih] = centerRow;

      HF1(13, layer);
      HF1(100*(layer+1)+11, centerRow);
      HF1(100*(layer+1)+14, hitpos.x());
      HF1(100*(layer+1)+15, residual);
      HF2(100*(layer+1)+16, hitpos.x(), residual);
      HF2(100*(layer+1)+17, calpos.x(), hitpos.y());
      HF1(100*(layer+1)+31, resi_vect.X());
      HF1(100*(layer+1)+32, resi_vect.Y());
      HF1(100*(layer+1)+33, resi_vect.Z());
      event.hitlayer[it][ih] = layer;
      event.hitpos_x[it][ih] = hitpos.x();
      event.hitpos_y[it][ih] = hitpos.y();
      event.hitpos_z[it][ih] = hitpos.z();
      event.calpos_x[it][ih] = calpos.x();
      event.calpos_y[it][ih] = calpos.y();
      event.calpos_z[it][ih] = calpos.z();
      event.residual[it][ih] = residual;
      event.residual_x[it][ih] = resi_vect.x();
      event.residual_y[it][ih] = resi_vect.y();
      event.residual_z[it][ih] = resi_vect.z();
      event.residual_horizontal[it][ih] = track->GetHorizontalResidual(ih);
      event.residual_vertical[it][ih] = track->GetVerticalResidual(ih);
      event.resolution_x[it][ih] = res_vect.x();
      event.resolution_y[it][ih] = res_vect.y();
      event.resolution_z[it][ih] = res_vect.z();
      event.resolution_horizontal[it][ih] = track->GetHorizontalResolution(ih);
      event.resolution_vertical[it][ih] = track->GetVerticalResolution(ih);
      event.theta_diff[it][ih] = track->GetAlpha(ih);
      event.pathhit[it][ih] = track->GetHitLength(ih);

      de += clde;
      Double_t dEdx = clde/event.pathhit[it][ih];
      dEdx_vect.push_back(dEdx);

#if Exclusive
      const TVector3& exres_vect = hit->GetResidualVectExclusive();
      Double_t exresidual = hit->GetResidualExclusive();
      event.exresidual[it][ih] = exresidual;
      event.exresidual_x[it][ih] = exres_vect.x();
      event.exresidual_y[it][ih] = exres_vect.y();
      event.exresidual_z[it][ih] = exres_vect.z();
      event.exresidual_horizontal[it][ih] = track->GetHorizontalResidualExclusive(ih);
      event.exresidual_vertical[it][ih] = track->GetVerticalResidualExclusive(ih);
#endif
    }

    event.dE[it] = de;
    std::sort(dEdx_vect.begin(), dEdx_vect.end());
    int n_truncated = (int)(dEdx_vect.size()*truncatedMean);
    for( int ih=0; ih<dEdx_vect.size(); ++ih ){
      if(ih<n_truncated) event.dEdx[it]+=dEdx_vect[ih];
    }
    event.dEdx[it]/=(double)n_truncated;
  }

  TVector3 vertex = Kinematics::MultitrackVertex(ntrack_intarget,
						 x0_vtx,
						 y0_vtx,
						 u0_vtx,
						 v0_vtx);

  event.ntTpc_inside = ntrack_intarget;
  event.prodvtx_x = vertex.x();
  event.prodvtx_y = vertex.y();
  event.prodvtx_z = vertex.z();

#if TrackSearchFailed
  Int_t failed_ntTpc = TPCAna.GetNTracksTPCFailed();
  event.failed_ntTpc = failed_ntTpc;
  event.failed_nhtrack.resize( failed_ntTpc );
  event.failed_x0Tpc.resize( failed_ntTpc );
  event.failed_y0Tpc.resize( failed_ntTpc );
  event.failed_u0Tpc.resize( failed_ntTpc );
  event.failed_v0Tpc.resize( failed_ntTpc );
  event.failed_hitlayer.resize( failed_ntTpc );
  event.failed_hitpos_x.resize( failed_ntTpc );
  event.failed_hitpos_y.resize( failed_ntTpc );
  event.failed_hitpos_z.resize( failed_ntTpc );
  event.failed_calpos_x.resize( failed_ntTpc );
  event.failed_calpos_y.resize( failed_ntTpc );
  event.failed_calpos_z.resize( failed_ntTpc );

  for(Int_t it=0; it<failed_ntTpc; ++it){
    auto track = TPCAna.GetTrackTPCFailed(it);
    if(!track) continue;
    Int_t nh = track->GetNHit();
    Double_t x0 = track->GetX0(), y0 = track->GetY0();
    Double_t u0 = track->GetU0(), v0 = track->GetV0();
    //Int_t fitflag = track->GetFitFlag();

    event.failed_nhtrack[it] = nh;
    event.failed_x0Tpc[it] = x0;
    event.failed_y0Tpc[it] = y0;
    event.failed_u0Tpc[it] = u0;
    event.failed_v0Tpc[it] = v0;
    event.failed_hitlayer[it].resize(nh);
    event.failed_hitpos_x[it].resize(nh);
    event.failed_hitpos_y[it].resize(nh);
    event.failed_hitpos_z[it].resize(nh);
    event.failed_calpos_x[it].resize(nh);
    event.failed_calpos_y[it].resize(nh);
    event.failed_calpos_z[it].resize(nh);

    for(int ih=0; ih<nh; ++ih){
      auto hit = track->GetHit(ih);
      if(!hit) continue;
      Int_t layer = hit->GetLayer();
      const TVector3& hitpos = hit->GetLocalHitPos();
      const TVector3& calpos = hit->GetLocalCalPos();

      event.failed_hitlayer[it][ih] = layer;
      event.failed_hitpos_x[it][ih] = hitpos.x();
      event.failed_hitpos_y[it][ih] = hitpos.y();
      event.failed_hitpos_z[it][ih] = hitpos.z();
      event.failed_calpos_x[it][ih] = calpos.x();
      event.failed_calpos_y[it][ih] = calpos.y();
      event.failed_calpos_z[it][ih] = calpos.z();
    }
  }
#endif

  HF1(1, event.status++);

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
  HB1(2, "HoughDist", 500, 0., 50.);
  HB1(3, "HoughDistY", 500, 0., 50.);
  HB1(4, "#Iteration of Tracking", 100, 0., 100.);
  HB1(5, "Fitting Flag", 10, 0., 10.);
  HB1(6, "Track Searching Time [ms]", 100, 0., 100.);
  HB1(7, "Track Fitting Time [ms]", 100, 0., 100.);
  HB1(8, "Minuit Output Status", 5, 0., 5.);

  HB1(10, "#Tracks TPC", 40, 0., 40.);
  HB1(11, "#Hits of Track TPC", 50, 0., 50.);
  HB1(12, "Chisqr TPC", 500, 0., 100.);
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

  const Int_t    NbinDe = 1000;
  const Double_t MinDe  =    0.;
  const Double_t MaxDe  = 2000.;

  const Int_t NbinClSize = 25;
  const Double_t MinClSize = 0;
  const Double_t MaxClSize = 25;
  const Int_t NbinDist = 60;
  const Double_t MinDist = -15.;
  const Double_t MaxDist = 15.;
  const Int_t NbinRatio = 100;
  const Double_t MinRatio = 0.;
  const Double_t MaxRatio = 1.;

  HB1(TPCClHid, "Cluster size;Cluster size;Counts", NbinClSize, MinClSize, MaxClSize);
  HB1(TPCClHid+1, "Cluster dE;Cluster dE;Counts", NbinDe, MinDe, MaxDe);
  HB2(TPCClHid+2, "Transverse diffusion;X_{cluster_center}-X_{pad};A/A_{sum}", NbinDist, MinDist, MaxDist, NbinRatio, MinRatio, MaxRatio);
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    HB1(TPCClHid+(layer+1)*1000, Form("Cluster size layer%d;Cluster size;Counts",layer), NbinClSize, MinClSize, MaxClSize);
    HB1(TPCClHid+(layer+1)*1000+1, Form("Cluster dE layer%d;Cluster dE;Counts",layer), NbinDe, MinDe, MaxDe);
    HB2(TPCClHid+(layer+1)*1000+2, Form("Transverse diffusion Layer%d;X_{cluster_center}-X_{pad};A/A_{sum}",layer), NbinDist, MinDist, MaxDist, NbinRatio, MinRatio, MaxRatio);
  }

  HBTree("tpc", "tree of DstTPCTracking");

  tree->Branch( "status", &event.status );
  tree->Branch( "runnum", &event.runnum );
  tree->Branch( "evnum", &event.evnum );
  tree->Branch( "trigpat", &event.trigpat );
  tree->Branch( "trigflag", &event.trigflag );
  tree->Branch( "clkTpc", &event.clkTpc);

#if RawHit
  tree->Branch( "nhTpc", &event.nhTpc );
  tree->Branch( "raw_hitpos_x", &event.raw_hitpos_x );
  tree->Branch( "raw_hitpos_y", &event.raw_hitpos_y );
  tree->Branch( "raw_hitpos_z", &event.raw_hitpos_z );
  tree->Branch( "raw_de", &event.raw_de );
  tree->Branch( "raw_padid", &event.raw_padid );
  tree->Branch( "raw_layer", &event.raw_layer );
  tree->Branch( "raw_row", &event.raw_row );
#endif

#if RawCluster
  tree->Branch( "nclTpc", &event.nclTpc );
  tree->Branch( "cluster_x", &event.cluster_x );
  tree->Branch( "cluster_y", &event.cluster_y );
  tree->Branch( "cluster_z", &event.cluster_z );
  tree->Branch( "cluster_de", &event.cluster_de );
  tree->Branch( "cluster_size", &event.cluster_size );
  tree->Branch( "cluster_layer", &event.cluster_layer );
  tree->Branch( "cluster_row_center", &event.cluster_row_center );
  tree->Branch( "cluster_mrow", &event.cluster_mrow );
  tree->Branch( "cluster_houghflag", &event.cluster_houghflag );
  tree->Branch( "cluster_de_center", &event.cluster_de_center );
  tree->Branch( "cluster_x_center", &event.cluster_x_center );
  tree->Branch( "cluster_y_center", &event.cluster_y_center );
  tree->Branch( "cluster_z_center", &event.cluster_z_center );
#endif

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

  tree->Branch("ntTpc_target", &event.ntTpc_inside);
  tree->Branch("prodvtx_x", &event.prodvtx_x);
  tree->Branch("prodvtx_y", &event.prodvtx_y);
  tree->Branch("prodvtx_z", &event.prodvtx_z);

#if Exclusive
  tree->Branch("exresidual", &event.exresidual);
  tree->Branch("exresidual_x", &event.exresidual_x);
  tree->Branch("exresidual_y", &event.exresidual_y);
  tree->Branch("exresidual_z", &event.exresidual_z);
  tree->Branch("exresidual_horizontal", &event.exresidual_horizontal);
  tree->Branch("exresidual_vertical", &event.exresidual_vertical);
#endif

#if TrackSearchFailed
  tree->Branch("failed_ntTpc", &event.failed_ntTpc);
  tree->Branch("failed_nhtrack", &event.failed_nhtrack);
  tree->Branch("failed_x0Tpc", &event.failed_x0Tpc);
  tree->Branch("failed_y0Tpc", &event.failed_y0Tpc);
  tree->Branch("failed_u0Tpc", &event.failed_u0Tpc);
  tree->Branch("failed_v0Tpc", &event.failed_v0Tpc);
  tree->Branch("failed_hitlayer", &event.failed_hitlayer);
  tree->Branch("failed_hitpos_x", &event.failed_hitpos_x);
  tree->Branch("failed_hitpos_y", &event.failed_hitpos_y);
  tree->Branch("failed_hitpos_z", &event.failed_hitpos_z);
  tree->Branch("failed_calpos_x", &event.failed_calpos_x);
  tree->Branch("failed_calpos_y", &event.failed_calpos_y);
  tree->Branch("failed_calpos_z", &event.failed_calpos_z);
#endif

  TTreeReaderCont[kTpcHit] = new TTreeReader("tpc", TFileCont[kTpcHit]);
  const auto& reader = TTreeReaderCont[kTpcHit];
  src.runnum = new TTreeReaderValue<Int_t>(*reader, "runnum");
  src.evnum = new TTreeReaderValue<Int_t>(*reader, "evnum");
  src.trigpat = new TTreeReaderValue<std::vector<Int_t>>(*reader, "trigpat");
  src.trigflag = new TTreeReaderValue<std::vector<Int_t>>(*reader, "trigflag");
  src.npadTpc = new TTreeReaderValue<Int_t>(*reader, "npadTpc");
  src.nhTpc = new TTreeReaderValue<Int_t>(*reader, "nhTpc");
  src.layerTpc = new TTreeReaderValue<std::vector<Int_t>>(*reader, "layerTpc");
  src.rowTpc = new TTreeReaderValue<std::vector<Int_t>>(*reader, "rowTpc");
  src.padTpc = new TTreeReaderValue<std::vector<Int_t>>(*reader, "padTpc");
  src.pedTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "pedTpc");
  src.rmsTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "rmsTpc");
  src.deTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "deTpc");
  src.tTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "tTpc");
  src.chisqrTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "chisqrTpc");
  src.clkTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "clkTpc");

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")   &&
     InitializeParameter<TPCParamMan>("TPCPRM") &&
     InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
     InitializeParameter<UserParamMan>("USER") &&
     InitializeParameter<HodoPHCMan>("HDPHC"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

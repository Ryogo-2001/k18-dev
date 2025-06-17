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
#include "FieldMan.hh"
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
#include "TPCPadHelper.hh"

#include "HypTPCFitter.hh"
#include "HypTPCTask.hh"

#include "DstHelper.hh"
#include <TRandom.h>
#include "DebugCounter.hh"

namespace
{
using namespace root;
using namespace dst;
const auto qnan = TMath::QuietNaN();
const std::string& class_name("GenfitLinearGeant4");
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();
debug::ObjectCounter& gCounter  = debug::ObjectCounter::GetInstance();

const Int_t MaxTPCHits = 500;
const Int_t MaxTPCTracks = 100;

//const bool IsWithRes = false;
const bool IsWithRes = true;

//For GenFit Setting
const bool Const_field = false; //Must be false for linear tracking
const Int_t verbosity = 0;//3;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");
const Double_t mom_smear = 0.034; //dp/p ~ 3.4%

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

  //Geant4 Generated Hits & Tracks
  Int_t g4nhTpc;               // Number of Hits
  Int_t g4ntTpc;               // Number of Tracks
  Int_t g4tidTpc[MaxTPCTracks];  // TrackId
  Int_t g4pidTpc[MaxTPCTracks];  // TrackPid
  Int_t g4nhtrack[MaxTPCHits]; // Number of Hits (in 1 tracks)
  Int_t g4nhHtof;
  Int_t g4tidHtof[MaxTPCHits];
  Int_t g4pidHtof[MaxTPCHits];
  Int_t g4didHtof[MaxTPCHits];
  Double_t g4xHtof[MaxTPCHits];
  Double_t g4yHtof[MaxTPCHits];
  Double_t g4zHtof[MaxTPCHits];
  Double_t g4lengthHtof[MaxTPCHits];
  Double_t g4tHtof[MaxTPCHits];

  //TrackSearchTPC()
  Int_t ntTpc;                   // Number of Tracks
  Int_t nhtrack[MaxTPCTracks]; // Number of Hits (in 1 tracks)
  Double_t chisqr[MaxTPCTracks];
  Double_t x0[MaxTPCTracks];
  Double_t y0[MaxTPCTracks];
  Double_t u0[MaxTPCTracks];
  Double_t v0[MaxTPCTracks];
  Double_t theta[MaxTPCTracks];

  Double_t hitresolution_x[MaxTPCTracks][MaxTPCHits];
  Double_t hitresolution_y[MaxTPCTracks][MaxTPCHits];
  Double_t hitresolution_z[MaxTPCTracks][MaxTPCHits];

  Int_t hitlayer[MaxTPCTracks][MaxTPCHits];
  Double_t hitpos_x[MaxTPCTracks][MaxTPCHits];
  Double_t hitpos_y[MaxTPCTracks][MaxTPCHits];
  Double_t hitpos_z[MaxTPCTracks][MaxTPCHits];
  Double_t calpos_x[MaxTPCTracks][MaxTPCHits];
  Double_t calpos_y[MaxTPCTracks][MaxTPCHits];
  Double_t calpos_z[MaxTPCTracks][MaxTPCHits];
  Double_t mom_x[MaxTPCTracks][MaxTPCHits];
  Double_t mom_y[MaxTPCTracks][MaxTPCHits];
  Double_t mom_z[MaxTPCTracks][MaxTPCHits];
  Double_t mom[MaxTPCTracks][MaxTPCHits];
  Double_t residual[MaxTPCTracks][MaxTPCHits];
  Double_t residual_x[MaxTPCTracks][MaxTPCHits];
  Double_t residual_y[MaxTPCTracks][MaxTPCHits];
  Double_t residual_z[MaxTPCTracks][MaxTPCHits];
  Double_t residual_px[MaxTPCTracks][MaxTPCHits];
  Double_t residual_py[MaxTPCTracks][MaxTPCHits];
  Double_t residual_pz[MaxTPCTracks][MaxTPCHits];
  Double_t residual_p[MaxTPCTracks][MaxTPCHits];

  Double_t g4pos_x[MaxTPCTracks][MaxTPCHits];
  Double_t g4pos_y[MaxTPCTracks][MaxTPCHits];
  Double_t g4pos_z[MaxTPCTracks][MaxTPCHits];
  Double_t g4mom_x[MaxTPCTracks][MaxTPCHits];
  Double_t g4mom_y[MaxTPCTracks][MaxTPCHits];
  Double_t g4mom_z[MaxTPCTracks][MaxTPCHits];
  Double_t g4mom[MaxTPCTracks][MaxTPCHits];
  Int_t g4tid[MaxTPCTracks][MaxTPCHits];
  Int_t g4pid[MaxTPCTracks][MaxTPCHits];
  Double_t tracklen[MaxTPCTracks][MaxTPCHits];
  Double_t tof[MaxTPCTracks][MaxTPCHits];

  //GenFit outputs
  Int_t GFstatus;
  Int_t GFntTpc;
  Int_t GFfitstatus[MaxTPCTracks];
  Int_t GFinside[MaxTPCTracks];
  Int_t GFndf[MaxTPCTracks];
  Int_t GFnhits[MaxTPCTracks];
  Double_t GFtracklen[MaxTPCTracks];
  Double_t GFchisqr[MaxTPCTracks];
  Double_t GFpval[MaxTPCTracks];
  Double_t GFtof[MaxTPCTracks];
  Double_t GFcharge[MaxTPCTracks];
  Double_t GFresidual_p[MaxTPCTracks];

  Int_t GFnhHtof[MaxTPCTracks];
  Int_t GFsegHtof[MaxTPCTracks][8];
  Double_t GFxHtof[MaxTPCTracks][8];
  Double_t GFyHtof[MaxTPCTracks][8];
  Double_t GFzHtof[MaxTPCTracks][8];
  Double_t GFextraHtof[MaxTPCTracks][8];
  Double_t GFtofHtof[MaxTPCTracks][8];

  Double_t GFpos_x[MaxTPCTracks][MaxTPCHits];
  Double_t GFpos_y[MaxTPCTracks][MaxTPCHits];
  Double_t GFpos_z[MaxTPCTracks][MaxTPCHits];
  Double_t GFmom_x[MaxTPCTracks][MaxTPCHits];
  Double_t GFmom_y[MaxTPCTracks][MaxTPCHits];
  Double_t GFmom_z[MaxTPCTracks][MaxTPCHits];
  Double_t GFmom[MaxTPCTracks][MaxTPCHits];
  Double_t GFpulls[MaxTPCTracks][5]; //extrapolation to the reference point
  Double_t GFresiduals[MaxTPCTracks][5]; //extrapolation to the reference point
  Double_t GFresidual_x[MaxTPCTracks][MaxTPCHits];
  Double_t GFresidual_y[MaxTPCTracks][MaxTPCHits];
  Double_t GFresidual_z[MaxTPCTracks][MaxTPCHits];

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

  Int_t ititpc[MaxTPCHits];
  Int_t idtpc[MaxTPCHits];
  Double_t xtpc[MaxTPCHits];//with resolution
  Double_t ytpc[MaxTPCHits];//with resolution
  Double_t ztpc[MaxTPCHits];//with resolution
  Double_t x0tpc[MaxTPCHits];//w/o resolution
  Double_t y0tpc[MaxTPCHits];//w/o resolution
  Double_t z0tpc[MaxTPCHits];//w/o resolution
  //  Double_t resoX[MaxTPCHits];
  Double_t pxtpc[MaxTPCHits];
  Double_t pytpc[MaxTPCHits];
  Double_t pztpc[MaxTPCHits];
  Double_t pptpc[MaxTPCHits];   // total mometum
  // Double_t masstpc[MaxTPCHits];   // mass TPC
  Double_t timetpc[MaxTPCHits];
  Double_t betatpc[MaxTPCHits];
  Double_t edeptpc[MaxTPCHits];
  // Double_t dedxtpc[MaxTPCHits];
  Double_t slengthtpc[MaxTPCHits];
  Double_t tlengthtpc[MaxTPCHits];
  // Int_t laytpc[MaxTPCHits];
  // Int_t rowtpc[MaxTPCHits];
  // Int_t parentID[MaxTPCHits];
  // Int_t iPadtpc[MaxTPCHits];//Pad number (0 origin)

  // Double_t xtpc_pad[MaxTPCHits];//pad center
  // Double_t ytpc_pad[MaxTPCHits];//pad center(dummy)
  // Double_t ztpc_pad[MaxTPCHits];//pad center

  // Double_t dxtpc_pad[MaxTPCHits];//x0tpc - xtpc
  // Double_t dytpc_pad[MaxTPCHits];//y0tpc - ytpc = 0 (dummy)
  // Double_t dztpc_pad[MaxTPCHits];//z0tpc - ztpc

  Int_t nhHtof;
  Int_t tidHtof[MaxTPCHits];
  Int_t pidHtof[MaxTPCHits];
  Int_t didHtof[MaxTPCHits];
  Double_t xHtof[MaxTPCHits];
  Double_t yHtof[MaxTPCHits];
  Double_t zHtof[MaxTPCHits];
  Double_t lengthHtof[MaxTPCHits];
  Double_t tHtof[MaxTPCHits];

};

namespace root
{
  Event  event;
  Src    src;
  TH1   *h[MaxHist];
  TTree *tree;
  enum eDetHid {
    g4Hid = 100,
    k18Hid = 200,
    genfitHid = 300
  };

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

  //Initiallize Geometry, Field, Fitter
  HypTPCFitter* fitter = new HypTPCFitter(tpcGeo.Data(),Const_field);
  //Initiallize the genfit track container
  HypTPCTask& GFtracks = HypTPCTask::GetInstance();
  GFtracks.SetVerbosity(verbosity);
  std::cout<<"GenFit verbosity = "<<"-1: Silent, 0: Minimum, 1: Errors only, 2: Errors and Warnings, 3: Verbose mode, long term debugging(default)"<<std::endl;
  std::cout<<"Current verbosity = "<<GFtracks.GetVerbosity()<<std::endl;

#if 0
  GFtracks.DebugMode();
#endif

  Int_t ievent = skip;
  for(; ievent<nevent && !CatchSignal::Stop(); ++ievent){
    gCounter.check();
    InitializeEvent();
    if(DstRead(ievent)) tree->Fill();
  }
  std::cout << "#D Event Number: " << std::setw(6)
	    << ievent << std::endl;
  DstClose();

  delete fitter;
  return EXIT_SUCCESS;
}

//_____________________________________________________________________
Bool_t
dst::InitializeEvent()
{
  event.status   = 0;
  event.evnum = 0;
  event.g4nhTpc = 0;
  event.g4ntTpc = 0;
  event.ntTpc = 0;

  event.GFstatus = 0;
  event.GFntTpc = 0;

  for(Int_t i=0; i<MaxTPCTracks; ++i){
    event.g4tidTpc[i] =0;
    event.g4nhtrack[i] =0;
    event.nhtrack[i] =0;
    event.GFinside[i] = -9999;
    event.chisqr[i] =qnan;
    event.x0[i] =qnan;
    event.y0[i] =qnan;
    event.u0[i] =qnan;
    event.v0[i] =qnan;
    event.theta[i] =qnan;
    event.GFfitstatus[i] =0;
    event.GFchisqr[i] =qnan;
    event.GFpval[i] =qnan;
    event.GFtracklen[i] =qnan;
    event.GFtof[i] =qnan;
    event.GFcharge[i] =qnan;
    event.GFresidual_p[i] =qnan;
    event.GFndf[i] =0;
    event.GFnhits[i] =0;
    for(Int_t j=0; j<5; ++j){
      event.GFpulls[i][j] =qnan;
      event.GFresiduals[i][j] =qnan;
    }
    event.GFnhHtof[i] = 0;
    for(Int_t j=0; j<8; ++j){
      event.GFsegHtof[i][j] =-9999;
      event.GFxHtof[i][j] =qnan;
      event.GFyHtof[i][j] =qnan;
      event.GFzHtof[i][j] =qnan;
      event.GFextraHtof[i][j] =qnan;
      event.GFtofHtof[i][j] =qnan;
    }

    for(Int_t j=0; j<MaxTPCHits; ++j){
      event.hitlayer[i][j] =-999;
      event.mom_x[i][j] =qnan;
      event.mom_y[i][j] =qnan;
      event.mom_z[i][j] =qnan;
      event.mom[i][j] =qnan;
      event.hitresolution_x[i][j] =qnan;
      event.hitresolution_y[i][j] =qnan;
      event.hitresolution_z[i][j] =qnan;
      event.g4mom_x[i][j] =qnan;
      event.g4mom_y[i][j] =qnan;
      event.g4mom_z[i][j] =qnan;
      event.g4pos_x[i][j] =qnan;
      event.g4pos_y[i][j] =qnan;
      event.g4pos_z[i][j] =qnan;
      event.g4mom[i][j] =qnan;
      event.hitpos_x[i][j] =qnan;
      event.hitpos_y[i][j] =qnan;
      event.hitpos_z[i][j] =qnan;
      event.calpos_x[i][j] =qnan;
      event.calpos_y[i][j] =qnan;
      event.calpos_z[i][j] =qnan;
      event.tracklen[i][j] =qnan;
      event.tof[i][j] =qnan;
      event.residual[i][j] =qnan;
      event.residual_x[i][j] =qnan;
      event.residual_y[i][j] =qnan;
      event.residual_z[i][j] =qnan;
      event.residual_px[i][j] =qnan;
      event.residual_py[i][j] =qnan;
      event.residual_pz[i][j] =qnan;
      event.residual_p[i][j] =qnan;
      event.g4tid[i][j] =0;
      event.g4pid[i][j] =0;

      event.GFmom_x[i][j] =qnan;
      event.GFmom_y[i][j] =qnan;
      event.GFmom_z[i][j] =qnan;
      event.GFmom[i][j] =qnan;
      event.GFpos_x[i][j] =qnan;
      event.GFpos_y[i][j] =qnan;
      event.GFpos_z[i][j] =qnan;
      event.GFresidual_x[i][j] =qnan;
      event.GFresidual_y[i][j] =qnan;
      event.GFresidual_z[i][j] =qnan;
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
  event.g4nhTpc = src.nhittpc;
  for(Int_t ihit=0; ihit<event.g4nhTpc; ++ihit){
    event.g4tidTpc[ihit] = src.ititpc[ihit];
    event.g4pidTpc[ihit] = src.idtpc[ihit];
    if(event.g4ntTpc<src.ititpc[ihit])
      event.g4ntTpc = src.ititpc[ihit];
    ++event.g4nhtrack[src.ititpc[ihit]-1];
  }

  event.g4nhHtof = src.nhHtof;
  for(Int_t ihit=0; ihit<event.g4nhHtof; ++ihit){
    event.g4tidHtof[ihit] = src.tidHtof[ihit];
    event.g4xHtof[ihit] = src.xHtof[ihit];
    event.g4yHtof[ihit] = src.yHtof[ihit];
    event.g4zHtof[ihit] = src.zHtof[ihit];
    event.g4lengthHtof[ihit] = src.lengthHtof[ihit];
    event.g4tHtof[ihit] = src.tHtof[ihit];
  }

  // Double_t u = src.pxPrm[0]/src.pzPrm[0];
  // Double_t v = src.pyPrm[0]/src.pzPrm[0];
  // Double_t cost = 1./std::sqrt(1.+u*u+v*v);
  // Double_t theta=std::acos(cost)*math::Rad2Deg();
  // if(theta>20.)
  //   return true;

  HF1( 1, event.status++ );

  if(src.nhittpc<5)
    return true;

  HF1( 1, event.status++ );

  HypTPCTask& GFtracks = HypTPCTask::GetInstance();

  HF1( 2, event.GFstatus++ );

  TPCAnalyzer TPCAna;
  if(IsWithRes) TPCAna.DecodeTPCHitsGeant4(src.nhittpc, src.xtpc, src.ytpc, src.ztpc, src.edeptpc, src.idtpc);
  else TPCAna.DecodeTPCHitsGeant4(src.nhittpc, src.x0tpc, src.y0tpc, src.z0tpc, src.edeptpc, src.idtpc);
  TPCAna.TrackSearchTPC();
  HF1( 1, event.status++ );
  HF1( 2, event.GFstatus++ );

  Int_t ntTpc = TPCAna.GetNTracksTPC();
  if(MaxTPCHits<ntTpc){
    std::cout << "#W " << func_name << " "
      	      << "too many ntTpc " << ntTpc << "/" << MaxTPCHits << std::endl;
    ntTpc = MaxTPCHits;
  }

  event.ntTpc = ntTpc;
  HF1( k18Hid, event.ntTpc );
  for(Int_t it=0; it<ntTpc; ++it){
    TPCLocalTrack *track= TPCAna.GetTrackTPC(it);
    if(!track) continue;
    Int_t nh = track->GetNHit();
    Double_t chisqr = track->GetChiSquare();
    Double_t x0 = track->GetX0(), y0 = track->GetY0();
    Double_t u0 = track->GetU0(), v0 = track->GetV0();
    Double_t theta = track->GetTheta();
    HF1(k18Hid+14, x0);
    HF1(k18Hid+15, y0);
    HF1(k18Hid+16, u0);
    HF1(k18Hid+17, v0);
    HF2(k18Hid+18, x0, u0);
    HF2(k18Hid+19, y0, v0);
    HF2(k18Hid+20, y0, x0);
    event.nhtrack[it] = nh;
    event.chisqr[it] = chisqr;
    event.x0[it] = x0;
    event.y0[it] = y0;
    event.u0[it] = u0;
    event.v0[it] = v0;
    event.theta[it] = theta;
    Double_t offset_tracklen = 0; Double_t offset_tof = 0;
    TVector3 mom;
    for(Int_t ih=0; ih<nh; ++ih){
      TPCLTrackHit *hit = track->GetHit(ih);
      if(!hit) continue;
      Int_t layerId = 0;
      layerId = hit->GetLayer();
      TVector3 hitpos = hit->GetLocalHitPos();
      TVector3 calpos = hit->GetLocalCalPos();
      //Double_t residual = hit->GetResidual();
      //TVector3 res_vect = hit->GetResidualVect();
      HF1(k18Hid+13, layerId);
      Double_t checker=9999;
      for( Int_t ih2=0; ih2<src.nhittpc; ++ih2 ){
	TVector3 setpos;
	if(IsWithRes) setpos = TVector3(src.xtpc[ih2], src.ytpc[ih2], src.ztpc[ih2]);
        else setpos = TVector3(src.x0tpc[ih2], src.y0tpc[ih2], src.z0tpc[ih2]);
	TVector3 d = setpos - hitpos;
	if(d.Mag()<checker){
	  checker = d.Mag();

	  event.g4pos_x[it][ih] = src.x0tpc[ih2];
	  event.g4pos_y[it][ih] = src.y0tpc[ih2];
	  event.g4pos_z[it][ih] = src.z0tpc[ih2];
	  event.g4mom_x[it][ih] = src.pxtpc[ih2];
	  event.g4mom_y[it][ih] = src.pytpc[ih2];
	  event.g4mom_z[it][ih] = src.pztpc[ih2];
	  event.g4mom[it][ih] = src.pptpc[ih2];
	  event.g4tid[it][ih] = src.ititpc[ih2];
	  event.g4pid[it][ih] = src.idtpc[ih2];

	  event.residual_x[it][ih] = calpos.x() - event.g4pos_x[it][ih];
	  event.residual_y[it][ih] = calpos.y() - event.g4pos_y[it][ih];
	  event.residual_z[it][ih] = calpos.z() - event.g4pos_z[it][ih];
	  event.residual[it][ih] = TMath::Sqrt(event.residual_x[it][ih]*event.residual_x[it][ih]+
					       event.residual_y[it][ih]*event.residual_y[it][ih]+
					       event.residual_z[it][ih]*event.residual_z[it][ih]);

	  event.tracklen[it][ih] = src.tlengthtpc[ih2];
	  event.tof[it][ih] = src.timetpc[ih2];
	}
      } //ih2
      if(ih==0){
	offset_tracklen = event.tracklen[it][ih];
	offset_tof = event.tof[it][ih];
	double smear_factor = gRandom->Gaus(1.0, mom_smear);
	mom = smear_factor*TVector3(event.g4mom_x[it][ih],event.g4mom_y[it][ih],event.g4mom_z[it][ih]);
      }
      HF2(1000*(layerId+1)+1, calpos.X(), calpos.Y());
      HF1(1000*(layerId+1)+2, calpos.X());
      HF1(1000*(layerId+1)+3, calpos.Y());
      HF1(1000*(layerId+1)+4, calpos.Z());
      HF1(1000*(layerId+1)+5, hitpos.X());
      HF1(1000*(layerId+1)+6, hitpos.Y());
      HF1(1000*(layerId+1)+7, hitpos.Z());
      HF1(1000*(layerId+1)+8, event.residual[it][ih]);
      HF1(1000*(layerId+1)+9, event.residual_x[it][ih]);
      HF1(1000*(layerId+1)+10, event.residual_y[it][ih]);
      HF1(1000*(layerId+1)+11, event.residual_z[it][ih]);

      event.residual_px[it][ih] = mom.x() - event.g4mom_x[it][it];
      event.residual_py[it][ih] = mom.y() - event.g4mom_y[it][it];
      event.residual_pz[it][ih] = mom.z() - event.g4mom_z[it][it];
      event.residual_p[it][ih] = mom.Mag() - event.g4mom[it][it];

      event.hitresolution_x[it][ih] = hit->GetResolutionVect().x();
      event.hitresolution_y[it][ih] = hit->GetResolutionVect().y();
      event.hitresolution_z[it][ih] = hit->GetResolutionVect().z();
      event.tracklen[it][ih] -= offset_tracklen;
      event.tof[it][ih] -= offset_tof;
      event.hitlayer[it][ih] = layerId;
      event.hitpos_x[it][ih] = hitpos.x();
      event.hitpos_y[it][ih] = hitpos.y();
      event.hitpos_z[it][ih] = hitpos.z();
      event.calpos_x[it][ih] = calpos.x();
      event.calpos_y[it][ih] = calpos.y();
      event.calpos_z[it][ih] = calpos.z();
      event.mom_x[it][ih] = mom.x();
      event.mom_y[it][ih] = mom.y();
      event.mom_z[it][ih] = mom.z();
      event.mom[it][ih] = mom.Mag();
      //event.residual[it][ih] = residual;
      //event.residual_x[it][ih] = res_vect.x();
      //event.residual_y[it][ih] = res_vect.y();
      //event.residual_z[it][ih] = res_vect.z();
    } //ih
    //Add tracks into the GenFit TrackCand
    if(event.g4pid[it][0]!=0) GFtracks.AddLinearTrack(event.g4pid[it][0], track, event.mom[it][0]);
  }//it

  HF1( g4Hid, event.g4ntTpc );
  for(Int_t i=0; i<event.g4ntTpc; ++i){
    HF1( g4Hid+1, event.g4nhtrack[i] );
  }

  HF1( 2, event.GFstatus++ );

  GFtracks.FitTracks();
  event.GFntTpc=GFtracks.GetNTrack();
  for( Int_t igf=0; igf<GFtracks.GetNTrack(); ++igf ){
    event.GFfitstatus[igf] = (int)GFtracks.TrackCheck(igf);
    HF1( 3, event.GFfitstatus[igf]);
    if(!GFtracks.TrackCheck(igf)) continue;
    for( Int_t ihit=0; ihit<GFtracks.GetNHits(igf); ++ihit ){
      TVector3 hit = GFtracks.GetPos(igf, ihit);
      TVector3 mom = GFtracks.GetMom(igf, ihit);
      event.GFmom_x[igf][ihit] = mom.x();
      event.GFmom_y[igf][ihit] = mom.y();
      event.GFmom_z[igf][ihit] = mom.z();
      event.GFmom[igf][ihit] = mom.Mag();
      event.GFpos_x[igf][ihit] = hit.x();
      event.GFpos_y[igf][ihit] = hit.y();
      event.GFpos_z[igf][ihit] = hit.z();
      event.GFresidual_x[igf][ihit] = hit.x() - event.g4pos_x[igf][ihit];
      event.GFresidual_y[igf][ihit] = hit.y() - event.g4pos_y[igf][ihit];
      event.GFresidual_z[igf][ihit] = hit.z() - event.g4pos_z[igf][ihit];
    }
    event.GFchisqr[igf]=GFtracks.GetChi2NDF(igf);
    event.GFcharge[igf]=GFtracks.GetCharge(igf);
    event.GFtof[igf]=GFtracks.GetTrackTOF(igf, 0, -1);
    event.GFtracklen[igf]=GFtracks.GetTrackLength(igf, 0, -1);
    event.GFpval[igf]=GFtracks.GetPvalue(igf);
    event.GFndf[igf]=GFtracks.GetNDF(igf);
    event.GFnhits[igf]=GFtracks.GetNHits(igf);

    int candidates = 0; int htofid[8]; TVector3 htofpos[8]; TVector3 htofmom[8]; double htoflen[8]; double htoftof[8];
    if(GFtracks.ExtrapolateToHTOF(igf, candidates, htofid, htofpos, htofmom, htoflen, htoftof)){
      event.GFnhHtof[igf]=candidates;
      for( Int_t ihit=0; ihit<candidates; ++ihit ){
	event.GFsegHtof[igf][ihit]=htofid[ihit];
	event.GFxHtof[igf][ihit]=htofpos[ihit].x();
	event.GFyHtof[igf][ihit]=htofpos[ihit].y();
	event.GFzHtof[igf][ihit]=htofpos[ihit].z();
	event.GFextraHtof[igf][ihit]=htoflen[ihit];
	event.GFtofHtof[igf][ihit]=htoftof[ihit];
      }
    }
    if(!GFtracks.IsInsideTarget(igf)) event.GFinside[igf]=1;
    else event.GFinside[igf]=0;

    TVector3 posv; TVector3 momv; double len; double tof;
    GFtracks.ExtrapolateToTarget(igf,posv,momv,len,tof);
    HF1( genfitHid+7, posv.x());
    HF1( genfitHid+8, posv.y());
    HF1( genfitHid+9, posv.z());
    HF1( genfitHid+15, event.GFpval[igf]);
  }
  HF1( genfitHid, event.GFntTpc );

  for(int it=0;it<event.g4ntTpc;it++){
    if(event.ntTpc==event.g4ntTpc && event.GFntTpc==event.g4ntTpc){
      TVector3 resolution(event.hitresolution_x[it][0],event.hitresolution_y[it][0],event.hitresolution_z[it][0]);
      TVector3 pos(event.g4pos_x[it][0],event.g4pos_y[it][0],event.g4pos_z[it][0]);
      TVector3 mom(event.g4mom_x[it][0],event.g4mom_y[it][0],event.g4mom_z[it][0]);
      double residual[5]; double pull[5]; double residual6D[6]; double pull6D[6];
      if(GFtracks.GetTrackPull(it,GFtracks.GetPDGcode(it),resolution,event.tracklen[it][event.nhtrack[it]-1],mom,pos,residual,pull,residual6D,pull6D)){
	for(Int_t i=0; i<5; ++i){
	  HF1( genfitHid+10+i, pull[i]);
	  event.GFpulls[it][i]=pull[i];
	  event.GFresiduals[it][i]=residual[i];
	  HF1( genfitHid+1, event.GFnhits[it]);
	  if(i!=0) HF1( genfitHid+15+i, residual[i]);
	}
	event.GFresidual_p[it]=residual[0];
	for(Int_t i=0; i<6; ++i){
	  HF1( genfitHid+20+i, pull6D[i]);
	  HF1( genfitHid+26+i, residual6D[i]);
	}
      }
      HF1( k18Hid+1, event.nhtrack[it]);
      HF1( k18Hid+2, event.mom[it][0]);
      HF1( k18Hid+3, event.residual_p[it][0]);
      HF1( k18Hid+4, event.chisqr[it]);
      HF1( k18Hid+5, event.tracklen[it][event.nhtrack[it]-1]);
      HF1( k18Hid+6, event.tof[it][event.nhtrack[it]-1]);
      HF1( k18Hid+7, event.residual_x[it][0]);
      HF1( k18Hid+8, event.residual_y[it][0]);
      HF1( k18Hid+9, event.residual_z[it][0]);
      HF1( k18Hid+10, event.residual_px[it][0]);
      HF1( k18Hid+11, event.residual_py[it][0]);
      HF1( k18Hid+12, event.residual_pz[it][0]);
      HF1( genfitHid+1, event.GFnhits[it]);
      HF1( genfitHid+2, event.GFmom[it][0]);
      HF1( genfitHid+3, event.GFresidual_p[it]);
      HF1( genfitHid+4, event.GFchisqr[it]);
      HF1( genfitHid+5, event.GFtracklen[it]);
      HF1( genfitHid+6, event.GFtof[it]);
    }
  }

  HF1( 2, event.GFstatus++ );

#if 0
  std::cout<<"[event]: "<<std::setw(6)<<ievent<<" ";
  std::cout<<"[g4nhTpc]: "<<std::setw(2)<<src.nhittpc<<" "<<std::endl;
#endif

  HF1(1, event.status++);

  GFtracks.Clear();
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

  for(Int_t i=1; i<=NumOfLayersTPC; ++i){
    // Tracking Histgrams
    TString title1 = Form("[K1.8] LayerId%2d Y%%Xcal;X cal. [mm]; Y cal. [mm]", i);
    TString title2 = Form("[K1.8] LayerId%2d Calc. Position x;Calpos x [mm]; Number of Tracks", i);
    TString title3 = Form("[K1.8] LayerId%2d Calc. Position y;Calpos y [mm]; Number of Tracks", i);
    TString title4 = Form("[K1.8] LayerId%2d Calc. Position z;Calpos z [mm]; Number of Tracks", i);
    TString title5 = Form("[K1.8] LayerId%2d Hit Position x;Hitpos x [mm]; Number of Tracks", i);
    TString title6 = Form("[K1.8] LayerId%2d Hit Position y;Hitpos y [mm]; Number of Tracks", i);
    TString title7 = Form("[K1.8] LayerId%2d Hit Position z;Hitpos z [mm]; Number of Tracks", i);
    TString title8 = Form("[K1.8] LayerId%2d Residual;Residual [mm]; Number of Tracks", i);
    TString title9 = Form("[K1.8] LayerId%2d Residual x;Residual x [mm]; Number of Tracks", i);
    TString title10 = Form("[K1.8] LayerId%2d Residual y;Residual y [mm]; Number of Tracks", i);
    TString title11 = Form("[K1.8] LayerId%2d Residual z;Residual z [mm]; Number of Tracks", i);
    HB2(1000*i+1, title1, 100, -250., 250., 100, -250., 250.);
    HB1(1000*i+2, title2, 200, -250., 250.);
    HB1(1000*i+3, title3, 200, -250., 250.);
    HB1(1000*i+4, title4, 200, -250., 250.);
    HB1(1000*i+5, title5, 200, -250., 250.);
    HB1(1000*i+6, title6, 200, -250., 250.);
    HB1(1000*i+7, title7, 200, -250., 250.);
    HB1(1000*i+8, title8, 200, -2.0, 2.0);
    HB1(1000*i+9, title9, 200, -2.0, 2.0);
    HB1(1000*i+10, title10, 200, -2.0, 2.0);
    HB1(1000*i+11, title11, 200, -2.0, 2.0);
  }

  HB1(1, "Status", 21, 0., 21. );
  HB1(2, "[GenFit] Status", 21, 0., 21. );
  HB1(3, "[Genfit] Fit Status", 2, 0., 2. );
  HB1(g4Hid, "[Geant4] #Track TPC", 10, 0., 10. );
  HB1(g4Hid+1, "[Geant4] #Hits of Track TPC", 33, 0., 33. );
  HB1(k18Hid, "[K1.8] #Track TPC", 10, 0., 10. );
  HB1(k18Hid+1, "[K1.8] #Hits of Track TPC", 33, 0., 33. );
  HB1(k18Hid+2, "[K1.8] Reconstructed P; P [GeV/c]; Counts [/0.001 GeV/c]", 1500, 0., 1.5 );
  HB1(k18Hid+3, "[K1.8] Reconstructed P Residual; P Residual [GeV/c]; Counts [/0.001 GeV/c]", 400, -0.2, 0.2 );
  HB1(k18Hid+4, "[K1.8] Chisqr/ndf;", 1000, 0, 50 );
  HB1(k18Hid+5, "[K1.8] Track Length; Length [mm]; Counts [/1 mm]", 500, 0, 500 );
  HB1(k18Hid+6, "[K1.8] Tof; Tof [ns]; Counts [/0.01 ns]", 500, 0, 5 );
  HB1(k18Hid+7, "[K1.8] Residual x;Residual x [mm]; Number of Tracks", 500, -1, 1);
  HB1(k18Hid+8, "[K1.8] Residual y;Residual y [mm]; Number of Tracks", 500, -1, 1);
  HB1(k18Hid+9, "[K1.8] Residual z;Residual z [mm]; Number of Tracks", 500, -1, 1);
  HB1(k18Hid+10, "[K1.8] p_x Residual;Residual [GeV/c]; Number of Tracks", 500, -0.1, 0.1);
  HB1(k18Hid+11, "[K1.8] p_y Residual;Residual [GeV/c]; Number of Tracks", 500, -0.1, 0.1);
  HB1(k18Hid+12, "[K1.8] p_z Residual;Residual [GeV/c]; Number of Tracks", 500, -0.1, 0.1);
  HB1(k18Hid+13, "[K1.8] LayerId TPC", 35, 0., 35.);
  HB1(k18Hid+14, "[K1.8] X0 TPC", 400, -100., 100.);
  HB1(k18Hid+15, "[K1.8] Y0 TPC", 400, -100., 100.);
  HB1(k18Hid+16, "[K1.8] U0 TPC", 200, -10, 10);
  HB1(k18Hid+17, "[K1.8] V0 TPC", 200, -10, 10);
  HB2(k18Hid+18, "[K1.8] U0%X0 TPC", 100, -100., 100., 100, -10, 10);
  HB2(k18Hid+19, "[K1.8] V0%Y0 TPC", 100, -100., 100., 100, -10, 10);
  HB2(k18Hid+20, "[K1.8] X0%Y0 TPC", 100, -100., 100., 100, -100, 100);

  HB1(genfitHid, "[GenFit] #Track TPC", 10, 0., 10. );
  HB1(genfitHid+1, "[GenFit] #Hits of Track TPC", 33, 0., 33. );
  HB1(genfitHid+2, "[GenFit] Reconstructed P; P [GeV/c]; Counts [/0.001 GeV/c]", 1500, 0., 1.5 );
  HB1(genfitHid+3, "[GenFit] Reconstructed P Residual; P Residual [GeV/c]; Counts [/0.001 GeV/c]", 400, -0.2, 0.2 );
  HB1(genfitHid+4, "[GenFit] Chisqr/ndf;", 1000, 0, 50 );
  HB1(genfitHid+5, "[GenFit] Track Length; Length [mm]; Counts [/1 mm]", 500, 0, 500 );
  HB1(genfitHid+6, "[GenFit] Tof; Tof [ns]; Counts [/0.01 ns]", 500, 0, 5 );
  HB1(genfitHid+7, "[GenFit] Vertex X; Vertex X [mm]; Counts [/0.1 mm]", 500, -25, 25 );
  HB1(genfitHid+8, "[GenFit] Vertex Y; Vertex Y [mm]; Counts [/0.1 mm]", 500, -25, 25 );
  HB1(genfitHid+9, "[GenFit] Vertex Z; Vertex Z [mm]; Counts [/0.1 mm]", 500, -143-25, -143+25 );
  HB1(genfitHid+10, "[GenFit] Track pull q/p;pull; Number of Tracks", 1000, -6, 6);
  HB1(genfitHid+11, "[GenFit] Track pull u';pull; Number of Tracks", 1000, -6, 6);
  HB1(genfitHid+12, "[GenFit] Track pull v';pull; Number of Tracks", 1000, -6, 6);
  HB1(genfitHid+13, "[GenFit] Track pull u;pull; Number of Tracks", 1000, -6, 6);
  HB1(genfitHid+14, "[GenFit] Track pull v;pull; Number of Tracks", 1000, -6, 6);
  HB1(genfitHid+15, "[GenFit] p-value;p-value; Number of Tracks", 100, -0.05, 1.05);
  HB1(genfitHid+16, "[GenFit] Residual u';Residual; Number of Tracks", 100, -0.03, 0.03);
  HB1(genfitHid+17, "[GenFit] Residual v';Residual; Number of Tracks", 100, -0.03, 0.03);
  HB1(genfitHid+18, "[GenFit] Residual u;Residual [cm]; Number of Tracks", 100, -0.1, 0.1);
  HB1(genfitHid+19, "[GenFit] Residual v;Residual [cm]; Number of Tracks", 100, -0.1, 0.1);
  HB1(genfitHid+20, "[GenFit] pull x;pull; Number of Tracks", 1000, -6, 6);
  HB1(genfitHid+21, "[GenFit] pull y;pull; Number of Tracks", 1000, -6, 6);
  HB1(genfitHid+22, "[GenFit] pull z;pull; Number of Tracks", 1000, -6, 6);
  HB1(genfitHid+23, "[GenFit] pull P_x;pull; Number of Tracks", 1000, -6, 6);
  HB1(genfitHid+24, "[GenFit] pull P_y;pull; Number of Tracks", 1000, -6, 6);
  HB1(genfitHid+25, "[GenFit] pull P_z;pull; Number of Tracks", 1000, -6, 6);
  HB1(genfitHid+26, "[GenFit] Residual x;Residual x [cm]; Number of Tracks", 500, -0.1, 0.1);
  HB1(genfitHid+27, "[GenFit] Residual y;Residual y [cm]; Number of Tracks", 500, -0.1, 0.1);
  HB1(genfitHid+28, "[GenFit] Residual z;Residual z [cm]; Number of Tracks", 500, -0.1, 0.1);
  HB1(genfitHid+29, "[GenFit] Residual P_x;Residual P_x [GeV/c]; Number of Tracks", 500, -0.1, 0.1);
  HB1(genfitHid+30, "[GenFit] Residual P_y;Residual P_y [GeV/c]; Number of Tracks", 500, -0.1, 0.1);
  HB1(genfitHid+31, "[GenFit] Residual P_z;Residual P_z [GeV/c]; Number of Tracks", 500, -0.1, 0.1);

  HBTree("tpc", "tree of TPCTracking using Geant4 input");

  tree->Branch("status", &event.status, "status/I");
  tree->Branch("evnum", &event.evnum, "evnum/I");

  //Geant4 Generated Hits & Tracks
  tree->Branch("g4ntTpc",&event.g4ntTpc,"g4ntTpc/I");
  tree->Branch("g4nhTpc",&event.g4nhTpc,"g4nhTpc/I");
  tree->Branch("g4tidTpc",event.g4tidTpc,"g4tidTpc[g4nhTpc]/I");
  tree->Branch("g4pidTpc",event.g4pidTpc,"g4pidTpc[g4nhTpc]/I");
  tree->Branch("g4nhtrack",event.g4nhtrack,"g4nhtrack[g4ntTpc]/I");
  tree->Branch("g4nhHtof",&event.g4nhHtof,"g4nhHtof/I");
  tree->Branch("g4tidHtof",event.g4tidHtof,"g4tidHtof[g4nhHtof]/I");
  tree->Branch("g4pidHtof",event.g4pidHtof,"g4pidHtof[g4nhHtof]/I");
  tree->Branch("g4xHtof",event.g4xHtof,"g4xHtof[g4nhHtof]/D");
  tree->Branch("g4yHtof",event.g4yHtof,"g4yHtof[g4nhHtof]/D");
  tree->Branch("g4zHtof",event.g4zHtof,"g4zHtof[g4nhHtof]/D");
  tree->Branch("g4lengthHtof",event.g4lengthHtof,"g4lengthHtof[g4nhHtof]/D");
  tree->Branch("g4tHtof",event.g4tHtof,"g4tHtof[g4nhHtof]/D");

  //TrackSearchTPCHelix()
  tree->Branch("ntTpc",&event.ntTpc,"ntTpc/I");
  tree->Branch("nhtrack",event.nhtrack,"nhtrack[ntTpc]/I");
  tree->Branch("chisqr",event.chisqr,"chisqr[ntTpc]/D");
  tree->Branch("x0",event.x0,"x0[ntTpc]/D");
  tree->Branch("y0",event.y0,"y0[ntTpc]/D");
  tree->Branch("u0",event.u0,"u0[ntTpc]/D");
  tree->Branch("v0",event.v0,"v0[ntTpc]/D");
  tree->Branch("theta",event.theta,"theta[ntTpc]/D");
  tree->Branch("g4pos_x",event.g4pos_x,Form("g4pos_x[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("g4pos_y",event.g4pos_y,Form("g4pos_y[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("g4pos_z",event.g4pos_z,Form("g4pos_z[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("g4mom_x",event.g4mom_x,Form("g4mom_x[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("g4mom_y",event.g4mom_y,Form("g4mom_y[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("g4mom_z",event.g4mom_z,Form("g4mom_z[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("g4mom",event.g4mom,Form("g4mom[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("g4tid",event.g4tid,Form("g4tid[ntTpc][%d]/I",MaxTPCHits));
  tree->Branch("g4pid",event.g4pid,Form("g4pid[ntTpc][%d]/I",MaxTPCHits));
  tree->Branch("hitlayer",event.hitlayer,Form("hitlayer[ntTpc][%d]/I",MaxTPCHits));
  tree->Branch("hitpos_x",event.hitpos_x,Form("hitpos_x[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("hitpos_y",event.hitpos_y,Form("hitpos_y[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("hitpos_z",event.hitpos_z,Form("hitpos_z[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("calpos_x",event.calpos_x,Form("calpos_x[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("calpos_y",event.calpos_y,Form("calpos_y[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("calpos_z",event.calpos_z,Form("calpos_z[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("mom_x",event.mom_x,Form("mom_x[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("mom_y",event.mom_y,Form("mom_y[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("mom_z",event.mom_z,Form("mom_z[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("hitresolution_x",event.hitresolution_x,Form("hitresolution_x[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("hitresolution_y",event.hitresolution_y,Form("hitresolution_y[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("hitresolution_z",event.hitresolution_z,Form("hitresolution_z[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("mom",event.mom,Form("mom[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("tracklen",event.tracklen,Form("tracklen[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("tof",event.tof,Form("tof[ntTpc][%d]/D",MaxTPCHits));

  //Residuals
  tree->Branch("residual",event.residual,Form("residual[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("residual_x",event.residual_x,Form("residual_x[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("residual_y",event.residual_y,Form("residual_y[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("residual_z",event.residual_z,Form("residual_z[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("residual_px",event.residual_px,Form("residual_px[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("residual_py",event.residual_py,Form("residual_py[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("residual_pz",event.residual_pz,Form("residual_pz[ntTpc][%d]/D",MaxTPCHits));
  tree->Branch("residual_p",event.residual_p,Form("residual_p[ntTpc][%d]/D",MaxTPCHits));

  //GenFit fit results
  tree->Branch("GFstatus",&event.GFstatus,"GFstatus/I");
  tree->Branch("GFntTpc",&event.GFntTpc,"GFntTpc/I");
  tree->Branch("GFinside",event.GFinside,"GFinside[GFntTpc]/I");
  tree->Branch("GFresidual_p",event.GFresidual_p,"GFresidual_p[GFntTpc]/D");
  tree->Branch("GFchisqr",event.GFchisqr,"GFchisqr[GFntTpc]/D");
  tree->Branch("GFpval",event.GFpval,"GFpval[GFntTpc]/D");
  tree->Branch("GFfitstatus", event.GFfitstatus,"GFfitstatus[GFntTpc]/I");
  tree->Branch("GFcharge",event.GFcharge,"GFcharge[GFntTpc]/D");
  tree->Branch("GFtracklen",event.GFtracklen,"GFtracklen[GFntTpc]/D");
  tree->Branch("GFtof",event.GFtof,"GFtof[GFntTpc]/D");
  tree->Branch("GFndf",event.GFndf,"GFndf[GFntTpc]/I");
  tree->Branch("GFnhits",event.GFnhits,"GFnhits[GFntTpc]/I");
  tree->Branch("GFpulls",event.GFpulls,"GFpulls[GFntTpc][5]/D");
  tree->Branch("GFresiduals",event.GFresiduals,"GFresiduals[GFntTpc][5]/D");
  tree->Branch("GFresidual_x",event.GFresidual_x,Form("GFresidual_x[GFntTpc][%d]/D",MaxTPCHits));
  tree->Branch("GFresidual_y",event.GFresidual_y,Form("GFresidual_y[GFntTpc][%d]/D",MaxTPCHits));
  tree->Branch("GFresidual_z",event.GFresidual_z,Form("GFresidual_z[GFntTpc][%d]/D",MaxTPCHits));
  tree->Branch("GFpos_x",event.GFpos_x,Form("GFpos_x[GFntTpc][%d]/D",MaxTPCHits));
  tree->Branch("GFpos_y",event.GFpos_y,Form("GFpos_y[GFntTpc][%d]/D",MaxTPCHits));
  tree->Branch("GFpos_z",event.GFpos_z,Form("GFpos_z[GFntTpc][%d]/D",MaxTPCHits));
  tree->Branch("GFmom_x",event.GFmom_x,Form("GFmom_x[GFntTpc][%d]/D",MaxTPCHits));
  tree->Branch("GFmom_y",event.GFmom_y,Form("GFmom_y[GFntTpc][%d]/D",MaxTPCHits));
  tree->Branch("GFmom_z",event.GFmom_z,Form("GFmom_z[GFntTpc][%d]/D",MaxTPCHits));
  tree->Branch("GFmom",event.GFmom,Form("GFmom[GFntTpc][%d]/D",MaxTPCHits));

  tree->Branch("GFnhHtof",event.GFnhHtof,"GFnhHtof[GFntTpc]/I");
  tree->Branch("GFsegHtof",event.GFsegHtof,Form("GFsegHtof[GFntTpc][%d]/I",8));
  tree->Branch("GFxHtof",event.GFxHtof,Form("GFxHtof[GFntTpc][%d]/D",8));
  tree->Branch("GFyHtof",event.GFyHtof,Form("GFyHtof[GFntTpc][%d]/D",8));
  tree->Branch("GFzHtof",event.GFzHtof,Form("GFzHtof[GFntTpc][%d]/D",8));
  tree->Branch("GFextraHtof",event.GFextraHtof,Form("GFextraHtof[GFntTpc][%d]/D",8));
  tree->Branch("GFtofHtof",event.GFtofHtof,Form("GFtofHtof[GFntTpc][%d]/D",8));

  /*
    tree->Branch("nPrm",&src.nhittpc,"nPrm/I");
    tree->Branch("xPrm",src.xPrm,"xPrm[nPrm]/D");
    tree->Branch("yPrm",src.yPrm,"yPrm[nPrm]/D");
    tree->Branch("zPrm",src.zPrm,"zPrm[nPrm]/D");
    tree->Branch("pxPrm",src.pxPrm,"pxPrm[nPrm]/D");
    tree->Branch("pyPrm",src.pyPrm,"pyPrm[nPrm]/D");
    tree->Branch("pzPrm",src.pzPrm,"pzPrm[nPrm]/D");
  */

  ////////// Bring Address From Dst
  TTreeCont[kTPCGeant]->SetBranchStatus("*", 0);
  TTreeCont[kTPCGeant]->SetBranchStatus("evnum",  1);
  /*
    TTreeCont[kTPCGeant]->SetBranchStatus("nhPrm",  1);
    TTreeCont[kTPCGeant]->SetBranchStatus("pidPrm",  1);
    TTreeCont[kTPCGeant]->SetBranchStatus("xPrm",  1);
    TTreeCont[kTPCGeant]->SetBranchStatus("yPrm",  1);
    TTreeCont[kTPCGeant]->SetBranchStatus("zPrm",  1);
    TTreeCont[kTPCGeant]->SetBranchStatus("pxPrm",  1);
    TTreeCont[kTPCGeant]->SetBranchStatus("pyPrm",  1);
    TTreeCont[kTPCGeant]->SetBranchStatus("pzPrm",  1);
  */
  TTreeCont[kTPCGeant]->SetBranchStatus("nhittpc",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ititpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("idtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("xtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ytpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ztpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("x0tpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("y0tpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("z0tpc", 1);
  //  TTreeCont[kTPCGeant]->SetBranchStatus("resoX", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pxtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pytpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pztpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pptpc", 1);   // total mometum
  // TTreeCont[kTPCGeant]->SetBranchStatus("masstpc", 1);   // mass TPC
  // TTreeCont[kTPCGeant]->SetBranchStatus("betatpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("timetpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("edeptpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("dedxtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("slengthtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("tlengthtpc", 1);
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

  TTreeCont[kTPCGeant]->SetBranchStatus("nhHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("tidHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pidHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("didHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("xHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("yHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("zHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("lengthHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("tHtof", 1);
  //TTreeCont[kTPCGeant]->SetBranchStatus("pxHtof", 1);
  //TTreeCont[kTPCGeant]->SetBranchStatus("pyHtof", 1);
  //TTreeCont[kTPCGeant]->SetBranchStatus("pzHtof", 1);
  //TTreeCont[kTPCGeant]->SetBranchStatus("ppHtof", 1);

  TTreeCont[kTPCGeant]->SetBranchAddress("evnum", &src.evnum);
  /*
    TTreeCont[kTPCGeant]->SetBranchAddress("nhPrm", &src.nhPrm);
    TTreeCont[kTPCGeant]->SetBranchAddress("pidPrm", &src.pidPrm);
    TTreeCont[kTPCGeant]->SetBranchAddress("xPrm", src.xPrm);
    TTreeCont[kTPCGeant]->SetBranchAddress("yPrm", src.yPrm);
    TTreeCont[kTPCGeant]->SetBranchAddress("zPrm", src.zPrm);
    TTreeCont[kTPCGeant]->SetBranchAddress("pxPrm", src.pxPrm);
    TTreeCont[kTPCGeant]->SetBranchAddress("pyPrm", src.pyPrm);
    TTreeCont[kTPCGeant]->SetBranchAddress("pzPrm", src.pzPrm);
  */
  TTreeCont[kTPCGeant]->SetBranchAddress("nhittpc", &src.nhittpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("ititpc", src.ititpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("idtpc", src.idtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("xtpc", src.xtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("ytpc", src.ytpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("ztpc", src.ztpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("x0tpc", src.x0tpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("y0tpc", src.y0tpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("z0tpc", src.z0tpc);
  //  TTreeCont[kTPCGeant]->SetBranchAddress("resoX", src.resoX);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxtpc", src.pxtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pytpc", src.pytpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pztpc", src.pztpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pptpc", src.pptpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("masstpc", src.masstpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("betatpc", src.betatpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("timetpc", src.timetpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("edeptpc", src.edeptpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("dedxtpc", src.dedxtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("slengthtpc", src.slengthtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("tlengthtpc", src.tlengthtpc);
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

  TTreeCont[kTPCGeant]->SetBranchAddress("nhHtof", &src.nhHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("tidHtof", src.tidHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pidHtof", src.pidHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("didHtof", src.didHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("xHtof", src.xHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("yHtof", src.yHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("zHtof", src.zHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("lengthHtof", src.lengthHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("tHtof", src.tHtof);

  return true;
}

//_____________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")   &&
     InitializeParameter<UserParamMan>("USER") &&
     InitializeParameter<FieldMan>("FLDMAP", "HSFLDMAP") &&
     InitializeParameter<HodoPHCMan>("HDPHC"));
}

//_____________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

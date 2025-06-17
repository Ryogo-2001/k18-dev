// -*- C++ -*-

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "TPCParamMan.hh"
#include "TPCHit.hh"
#include "HodoPHCMan.hh"
#include "DCAnalyzer.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCPositionCorrector.hh"

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
const std::string& class_name("GenfitSkeleton");
using hddaq::unpacker::GUnpacker;
const auto& gUnpacker = GUnpacker::get_instance();
ConfMan& gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();
debug::ObjectCounter& gCounter  = debug::ObjectCounter::GetInstance();

const Int_t MaxTPCHits = 500;
const Int_t MaxTPCTracks = 100;

//const bool IsWithRes = false;
const bool IsWithRes = true;

//For GenFit Setting
const bool Const_field = true;
const Int_t verbosity = 3;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kTpcHit, kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[TPCHit]",
  "[OutFile]" };
std::vector<TString> TreeName =
{ "", "", "tpc", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________
struct Event
{

  Int_t status;
  Int_t runnum;
  Int_t evnum;
  std::vector<Int_t> trigpat;
  std::vector<Int_t> trigflag;
  std::vector<Double_t> clkTpc;
  Int_t nhTpc;

  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)

  //GenFit outputs
  Int_t GFstatus;
  Int_t GFntTpc;
  Int_t GFfitstatus[MaxTPCTracks];
  Double_t GFtracklen[MaxTPCTracks];
  Double_t GFchisqr[MaxTPCTracks];
  Double_t GFtof[MaxTPCTracks];
  Double_t GFpos_x[MaxTPCTracks][MaxTPCHits];
  Double_t GFpos_y[MaxTPCTracks][MaxTPCHits];
  Double_t GFpos_z[MaxTPCTracks][MaxTPCHits];
  Double_t GFmom_x[MaxTPCTracks][MaxTPCHits];
  Double_t GFmom_y[MaxTPCTracks][MaxTPCHits];
  Double_t GFmom_z[MaxTPCTracks][MaxTPCHits];
  Double_t GFmom[MaxTPCTracks][MaxTPCHits];

  void clear()
  {
    runnum = 0;
    evnum = 0;
    status = 0;
    clkTpc.clear();
    nhTpc = 0;
    ntTpc = 0;
    trigpat.clear();
    trigflag.clear();
    nhtrack.clear();
  }

};

//_____________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* runnum;
  TTreeReaderValue<Int_t>* evnum;
  TTreeReaderValue<std::vector<Int_t>>* trigpat;
  TTreeReaderValue<std::vector<Int_t>>* trigflag;
  TTreeReaderValue<Int_t>* npadTpc;   // number of pads
  TTreeReaderValue<Int_t>* nhTpc;     // number of hits
  TTreeReaderValue<std::vector<Int_t>>* layerTpc;     // layer id
  TTreeReaderValue<std::vector<Int_t>>* rowTpc;       // row id
  TTreeReaderValue<std::vector<Int_t>>* padTpc;       // pad id
  TTreeReaderValue<std::vector<Double_t>>* pedTpc;    // pedestal
  TTreeReaderValue<std::vector<Double_t>>* rmsTpc;    // rms
  TTreeReaderValue<std::vector<Double_t>>* deTpc;     // dE
  TTreeReaderValue<std::vector<Double_t>>* tTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* ctTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* chisqrTpc; // chi^2 of signal fitting
  TTreeReaderValue<std::vector<Double_t>>* clkTpc;

};

namespace root
{
  Event  event;
  Src    src;
  TH1   *h[MaxHist];
  TTree *tree;
  enum eDetHid {
    TPCHid = 100000,
    genfitHid = 200000,
  };
}

//_____________________________________________________________________
Int_t
main( Int_t argc, char **argv )
{
  std::vector<std::string> arg( argv, argv+argc );

  if( !CheckArg( arg ) )
    return EXIT_FAILURE;
  if( !DstOpen( arg ) )
    return EXIT_FAILURE;
  if( !gConf.Initialize( arg[kConfFile] ) )
    return EXIT_FAILURE;
  if( !gConf.InitializeUnpacker() )
    return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries( TTreeCont );
  if (max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();

  //Initiallize Geometry, Field, Fitter
  HypTPCFitter* fitter = new HypTPCFitter(tpcGeo.Data(),Const_field);
  //Initiallize the genfit track container
  HypTPCTask& GFtracks = HypTPCTask::GetInstance();
  GFtracks.SetVerbosity(verbosity);
  std::cout<<"GenFit verbosity = "<<"-1: Silent, 0: Minimum, 1: Errors only, 2: Errors and Warnings, 3: Verbose mode, long term debugging(default)"<<std::endl;
  std::cout<<"Current verbosity = "<<GFtracks.GetVerbosity()<<std::endl;

  Int_t ievent = skip;
  for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
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
bool
dst::InitializeEvent( void )
{
  event.clear();

  event.GFstatus = 0;
  event.GFntTpc = 0;
  for(Int_t i=0; i<MaxTPCTracks; ++i){
    event.GFfitstatus[i] = 0;
    event.GFchisqr[i] =qnan;
    event.GFtracklen[i] =qnan;
    event.GFtof[i] =qnan;
    for(Int_t j=0; j<MaxTPCHits; ++j){
      event.GFmom_x[i][j] =qnan;
      event.GFmom_y[i][j] =qnan;
      event.GFmom_z[i][j] =qnan;
      event.GFmom[i][j] =qnan;
      event.GFpos_x[i][j] =qnan;
      event.GFpos_y[i][j] =qnan;
      event.GFpos_z[i][j] =qnan;
    }
  }

  return true;
}

//_____________________________________________________________________
bool
dst::DstOpen( std::vector<std::string> arg )
{
  Int_t open_file = 0;
  Int_t open_tree = 0;
  for( std::size_t i=0; i<nArgc; ++i ){
    if( i==kProcess || i==kConfFile || i==kOutFile ) continue;
    open_file += OpenFile( TFileCont[i], arg[i] );
    open_tree += OpenTree( TFileCont[i], TTreeCont[i], TreeName[i] );
  }

  if( open_file!=open_tree || open_file!=nArgc-3 )
    return false;
  if( !CheckEntries( TTreeCont ) )
    return false;

  TFileCont[kOutFile] = new TFile( arg[kOutFile].c_str(), "recreate" );

  return true;
}

//_____________________________________________________________________
bool
dst::DstRead( Int_t ievent )
{
  static const std::string func_name("["+class_name+"::"+__func__+"]");


  if( ievent%10000==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }


  GetEntry(ievent);

  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;

  HF1( 1, event.status++ );

  if( **src.nhTpc == 0 )
    return true;


  HF1( 2, event.GFstatus++ );

  event.clkTpc = **src.clkTpc;
  Double_t clock = event.clkTpc.at(0);

  DCAnalyzer DCAna;
  DCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, clock);
#if 0 //HoughYcut
  DCAna.HoughYCut(min_ycut, max_ycut);
#endif

  Int_t ntTpc = DCAna.GetNTracksTPCHelix();
  if( MaxHits<ntTpc ){
    std::cout << "#W " << func_name << " "
      	      << "too many ntTpc " << ntTpc << "/" << MaxHits << std::endl;
    ntTpc = MaxHits;
  }

  HypTPCTask& GFtracks = HypTPCTask::GetInstance();

  HF1( 2, event.GFstatus++ );

  event.ntTpc = ntTpc;
  event.nhtrack.resize( ntTpc );

  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp= DCAna.GetTrackTPCHelix(it);
    if(!tp) continue;
    int nh = tp -> GetNHit();
    event.nhtrack[it] = nh;

    //Add tracks into the GenFit TrackCand
#if 0
    event.pid[it]=tpc -> GetPid();
    std::vector<Int_t> pdgcode;
    Kinematics::HypTPCPID_PDGCode(event.charge[it], particleID, pdgcode);
    GFtracks.AddHelixTrack(pdgcode, tp);
#endif
  } //it

  HF1( 2, event.GFstatus++ );

  GFtracks.FitTracks();
  HF1( genfitHid, event.GFntTpc );
  for( Int_t igf=0; igf<GFtracks.GetNTrack(); ++igf ){
    event.GFfitstatus[igf] = (int)GFtracks.TrackCheck(igf);
    HF1( 3, event.GFfitstatus[igf]);
    if(!GFtracks.TrackCheck(igf)) continue;
    event.GFntTpc++;
    event.GFchisqr[igf]=GFtracks.GetChi2NDF(igf);
    event.GFtracklen[igf]=GFtracks.GetTrackLength(igf, 0, -1);
    event.GFtof[igf]=GFtracks.GetTrackTOF(igf, 0, -1);
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
    }
  }
  HF1( 2, event.GFstatus++ );

#if 0
  std::cout<<"[event]: "<<std::setw(6)<<ievent<<" ";
  std::cout<<"[g4nhTpc]: "<<std::setw(2)<<src.nhittpc<<" "<<std::endl;
#endif

  HF1( 1, event.status++ );

  GFtracks.Clear();
  return true;
}

//_____________________________________________________________________
bool
dst::DstClose( void )
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  const std::size_t n = TFileCont.size();
  for( std::size_t i=0; i<n; ++i ){
    if( TTreeCont[i] ) delete TTreeCont[i];
    if( TFileCont[i] ) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________
bool
ConfMan::InitializeHistograms( void )
{
  HB1(1, "Status", 20, 0., 20. );
  HB1(2, "[GenFit] Status", 20, 0., 20. );
  HB1(3, "[Genfit] Fit Status", 2, 0., 2. );
  HBTree( "dst", "tree of GenFitSkeleton" );

  //TrackSearchTPCHelix()
  tree->Branch("ntTpc",&event.ntTpc,"ntTpc/I");
  tree->Branch("nhtrack",&event.nhtrack);

  //GenFit fit results
  tree->Branch("GFstatus",&event.GFstatus,"GFstatus/I");
  tree->Branch("GFntTpc",&event.GFntTpc,"GFntTpc/I");
  tree->Branch("GFfitstatus",event.GFfitstatus,"GFfitstatus[GFntTpc]/I");
  tree->Branch("GFchisqr",event.GFchisqr,"GFchisqr[GFntTpc]/D");
  tree->Branch("GFmom",event.GFmom,"GFmom[GFntTpc]/D");
  tree->Branch("GFtracklen",event.GFtracklen,"GFtracklen[GFntTpc]/D");
  tree->Branch("GFtof",event.GFtof,"GFtof[GFntTpc]/D");
  tree->Branch("GFpos_x",event.GFpos_x,Form("GFpos_x[GFntTpc][%d]/D",MaxTPCHits));
  tree->Branch("GFpos_y",event.GFpos_y,Form("GFpos_y[GFntTpc][%d]/D",MaxTPCHits));
  tree->Branch("GFpos_z",event.GFpos_z,Form("GFpos_z[GFntTpc][%d]/D",MaxTPCHits));
  tree->Branch("GFmom_x",event.GFmom_x,Form("GFmom_x[GFntTpc][%d]/D",MaxTPCHits));
  tree->Branch("GFmom_y",event.GFmom_y,Form("GFmom_y[GFntTpc][%d]/D",MaxTPCHits));
  tree->Branch("GFmom_z",event.GFmom_z,Form("GFmom_z[GFntTpc][%d]/D",MaxTPCHits));
  tree->Branch("GFmom",event.GFmom,Form("GFmom[GFntTpc][%d]/D",MaxTPCHits));

  ////////// Bring Address From Dst
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
  src.ctTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ctTpc");
  src.chisqrTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "chisqrTpc");
  src.clkTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "clkTpc");

  return true;
}

//_____________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")   &&
     InitializeParameter<TPCParamMan>("TPCPRM") &&
     InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
     InitializeParameter<UserParamMan>("USER") &&
     InitializeParameter<HodoPHCMan>("HDPHC"));
}

//_____________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}

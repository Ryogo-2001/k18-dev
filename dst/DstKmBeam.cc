// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <TLorentzVector.h>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "FieldMan.hh"
#include "DatabasePDG.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DstHelper.hh"
#include "HodoPHCMan.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCAnalyzer.hh"
#include "TPCCluster.hh"
#include "TPCVertex.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kTpc, kK18HSTracking, kOutFile, nArgc
};
std::vector<TString> ArgName =
  { "[Process]", "[ConfFile]", "[GenfitCarbon]", "[K18HSTracking]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc", "k18track","" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________________
struct Event
{

  Int_t runnum;
  Int_t evnum;
  std::vector<Int_t> trigpat;
  std::vector<Int_t> trigflag;
  Bool_t xiflag;
  Double_t ximass;

  Int_t ntK18;
  std::vector<Double_t> chisqrK18;
  std::vector<Double_t> xvpK18;
  std::vector<Double_t> yvpK18;
  std::vector<Double_t> zvpK18;
  std::vector<Double_t> xtgtK18;
  std::vector<Double_t> ytgtK18;
  std::vector<Double_t> ztgtK18;

  void clear( void )
  {
    runnum = 0;
    evnum = 0;
    trigpat.clear();
    trigflag.clear();
    xiflag = false;
    ximass = qnan;

    ntK18 = 0;
    chisqrK18.clear();
    xvpK18.clear();
    yvpK18.clear();
    zvpK18.clear();
    xtgtK18.clear();
    ytgtK18.clear();
    ztgtK18.clear();
  }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* runnum;
  TTreeReaderValue<Int_t>* evnum;
  TTreeReaderValue<std::vector<Int_t>>* trigpat;
  TTreeReaderValue<std::vector<Int_t>>* trigflag;

  TTreeReaderValue<Bool_t>* xiflag;
  TTreeReaderValue<Double_t>* ximass;

  Int_t    ntK18;
  Double_t chisqrK18[MaxHits];
  Double_t xtgtK18[MaxHits];
  Double_t ytgtK18[MaxHits];
  Double_t ztgtK18[MaxHits];
  Double_t xvpK18[MaxHits];
  Double_t yvpK18[MaxHits];
  Double_t zvpK18[MaxHits];

};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
  enum eDetHid {
    TPCHid    = 100000,
  };
}

//_____________________________________________________________________________
int
main( int argc, char **argv )
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
  if(skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries( TTreeCont );
  if(max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();

  Int_t ievent = skip;
  for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
    gCounter.check();
    InitializeEvent();
    if( DstRead( ievent ) ) tree->Fill();
  }

  std::cout << "#D Event Number: " << std::setw(6)
            << ievent << std::endl;

  DstClose();

  return EXIT_SUCCESS;
}

//_____________________________________________________________________________
Bool_t
dst::InitializeEvent( void )
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstOpen( std::vector<std::string> arg )
{
  int open_file = 0;
  int open_tree = 0;
  for( Int_t i=0; i<nArgc; ++i ){
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

//_____________________________________________________________________________
Bool_t
dst::DstRead( int ievent )
{
  static const auto MaxChisqrBcOut = gUser.GetParameter("MaxChisqrBcOut");
  static const auto MaxChisqrKurama = gUser.GetParameter("MaxChisqrKurama");

  //if( ievent%1000==0 ){
  if( ievent%1==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;

  event.xiflag = **src.xiflag;
  event.ximass = **src.ximass;

  event.ntK18 = src.ntK18;
  for(Int_t itK18=0; itK18<src.ntK18; ++itK18){
    event.chisqrK18.push_back(src.chisqrK18[itK18]);
    event.xtgtK18.push_back(src.xtgtK18[itK18]);
    event.ytgtK18.push_back(src.ytgtK18[itK18]);
    event.ztgtK18.push_back(src.ztgtK18[itK18]);
    event.xvpK18.push_back(src.xvpK18[itK18]);
    event.yvpK18.push_back(src.yvpK18[itK18]);
    event.zvpK18.push_back(src.zvpK18[itK18]);
  }

  if(event.ntK18==1){
    HF1(10, 1);
    HF1(11, event.xtgtK18[0]);
    HF1(12, event.xvpK18[0]);
    HF1(13, event.ytgtK18[0]);
    HF1(14, event.yvpK18[0]);
    HF2(15, event.xtgtK18[0], event.ytgtK18[0]);
    HF2(16, event.xvpK18[0], event.yvpK18[0]);
    if(event.trigflag[20] > 0){
      HF1(20, 1);
      HF1(21, event.xtgtK18[0]);
      HF1(22, event.xvpK18[0]);
      HF1(23, event.ytgtK18[0]);
      HF1(24, event.yvpK18[0]);
      HF2(25, event.xtgtK18[0], event.ytgtK18[0]);
      HF2(26, event.xvpK18[0], event.yvpK18[0]);
    }
    if(event.trigflag[21] > 0){
      HF1(30, 1);
      HF1(31, event.xtgtK18[0]);
      HF1(32, event.xvpK18[0]);
      HF1(33, event.ytgtK18[0]);
      HF1(34, event.yvpK18[0]);
      HF2(35, event.xtgtK18[0], event.ytgtK18[0]);
      HF2(36, event.xvpK18[0], event.yvpK18[0]);
    }
    if(event.trigflag[22] > 0){
      HF1(40, 1);
      HF1(41, event.xtgtK18[0]);
      HF1(42, event.xvpK18[0]);
      HF1(43, event.ytgtK18[0]);
      HF1(44, event.yvpK18[0]);
      HF2(45, event.xtgtK18[0], event.ytgtK18[0]);
      HF2(46, event.xvpK18[0], event.yvpK18[0]);
    }
    if(event.trigflag[23] > 0){
      HF1(50, 1);
      HF1(51, event.xtgtK18[0]);
      HF1(52, event.xvpK18[0]);
      HF1(53, event.ytgtK18[0]);
      HF1(54, event.yvpK18[0]);
      HF2(55, event.xtgtK18[0], event.ytgtK18[0]);
      HF2(56, event.xvpK18[0], event.yvpK18[0]);
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstClose( void )
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  const Int_t n = TFileCont.size();
  for( Int_t i=0; i<n; ++i ){
    if( TTreeReaderCont[i] ) delete TTreeReaderCont[i];
    if( TTreeCont[i] ) delete TTreeCont[i];
    if( TFileCont[i] ) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms( void )
{

  HBTree( "tpc", "tree of DstKmBeam" );
  tree->Branch( "runnum", &event.runnum );
  tree->Branch( "evnum", &event.evnum );
  tree->Branch( "trigpat", &event.trigpat );
  tree->Branch( "trigflag", &event.trigflag );

  tree->Branch( "Xiflag", &event.xiflag );
  tree->Branch( "XiMass", &event.ximass );

  tree->Branch( "ntK18", &event.ntK18);
  tree->Branch( "chisqrK18", &event.chisqrK18);
  tree->Branch( "xvpK18", &event.xvpK18);
  tree->Branch( "yvpK18", &event.yvpK18);
  tree->Branch( "zvpK18", &event.zvpK18);
  tree->Branch( "xtgtK18", &event.xtgtK18);
  tree->Branch( "ytgtK18", &event.ytgtK18);
  tree->Branch( "ztgtK18", &event.ztgtK18);

  TTreeReaderCont[kTpc] = new TTreeReader( "tpc", TFileCont[kTpc] );
  const auto& reader = TTreeReaderCont[kTpc];

  src.runnum = new TTreeReaderValue<Int_t>( *reader, "runnum" );
  src.evnum = new TTreeReaderValue<Int_t>( *reader, "evnum" );
  src.trigpat = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trigpat" );
  src.trigflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trigflag" );
  src.xiflag = new TTreeReaderValue<Bool_t>( *reader, "Xiflag" );
  src.ximass = new TTreeReaderValue<Double_t>( *reader, "XiMass" );

  TTreeCont[kK18HSTracking]->SetBranchStatus("*", 0);
  TTreeCont[kK18HSTracking]->SetBranchStatus("ntK18",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("chisqrK18",  1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("xtgtHS",     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("ytgtHS",     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("ztgtHS",     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("xvpHS",     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("yvpHS",     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("zvpHS",     1);

  TTreeCont[kK18HSTracking]->SetBranchAddress("ntK18", &src.ntK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("chisqrK18", src.chisqrK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("xtgtHS", src.xtgtK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("ytgtHS", src.ytgtK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("ztgtHS", src.ztgtK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("xvpHS", src.xvpK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("yvpHS", src.yvpK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("zvpHS", src.zvpK18);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<UserParamMan>("USER") );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess( void )
{
  return true;
}

/**
 *  file: DstTPCTrackingHelixGeant4.cc
 *  date: 2021.01.22
 *
 */

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include <filesystem_util.hh>

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
#include "TPCLocalTrack_Helix.hh"
#include "TPCHit.hh"

#include "DstHelper.hh"
#include <TRandom.h>
#include "DebugCounter.hh"

namespace
{
  using namespace root;
  using namespace dst;
  const std::string& class_name("DstTPCTrackingHelixGeant4");
  ConfMan&            gConf = ConfMan::GetInstance();
  const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const HodoPHCMan&   gPHC  = HodoPHCMan::GetInstance(); 
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
}

//_____________________________________________________________________
struct Event
{
  
  Int_t evnum;
  Int_t status;
  Int_t nhittpc;                 // Number of Hits
  Int_t nttpc;                   // Number of Tracks
  Int_t nhit_track[MaxTPCTracks]; // Number of Hits (in 1 tracks)
  Double_t chisqr[MaxTPCTracks];
  Double_t helix_cx[MaxTPCTracks];
  Double_t helix_cy[MaxTPCTracks];
  Double_t helix_z0[MaxTPCTracks];
  Double_t helix_r[MaxTPCTracks];
  Double_t helix_dz[MaxTPCTracks];
  Double_t mom0_x[MaxTPCTracks];//Helix momentum at Y = 0
  Double_t mom0_y[MaxTPCTracks];
  Double_t mom0_z[MaxTPCTracks];
  Int_t hitlayer[MaxTPCTracks][MaxTPCnHits];
  Double_t hitpos_x[MaxTPCTracks][MaxTPCnHits];
  Double_t hitpos_y[MaxTPCTracks][MaxTPCnHits];
  Double_t hitpos_z[MaxTPCTracks][MaxTPCnHits];
  Double_t calpos_x[MaxTPCTracks][MaxTPCnHits];
  Double_t calpos_y[MaxTPCTracks][MaxTPCnHits];
  Double_t calpos_z[MaxTPCTracks][MaxTPCnHits];
  Double_t mom_x[MaxTPCTracks][MaxTPCnHits];
  Double_t mom_y[MaxTPCTracks][MaxTPCnHits];
  Double_t mom_z[MaxTPCTracks][MaxTPCnHits];
  Double_t momg_x[MaxTPCTracks][MaxTPCnHits];
  Double_t momg_y[MaxTPCTracks][MaxTPCnHits];
  Double_t momg_z[MaxTPCTracks][MaxTPCnHits];

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
  // Double_t pptpc[MaxTPCHits];   // total mometum
  // Double_t masstpc[MaxTPCHits];   // mass TPC
  // Double_t betatpc[MaxTPCHits];
  Double_t edeptpc[MaxTPCHits];
  // Double_t dedxtpc[MaxTPCHits];
  // Double_t slengthtpc[MaxTPCHits];
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
main( int argc, char **argv )
{
  std::vector<std::string> arg( argv, argv+argc );

  if( !CheckArg( arg ) )
    return EXIT_FAILURE;
  if( !DstOpen( arg ) )
    return EXIT_FAILURE;
  if( !gConf.Initialize( arg[kConfFile] ) )
    return EXIT_FAILURE;

  int nevent = GetEntries( TTreeCont );

  CatchSignal::Set();

  int ievent = 0;

  // for(int ii=0; ii<100; ++ii){
  //   std::cout<<"ii="<<ii<<std::endl;
  //   ievent = 0;
  for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
    gCounter.check();
    InitializeEvent();
    if( DstRead( ievent ) ) tree->Fill();
  }
  //  }

  std::cout << "#D Event Number: " << std::setw(6)
	    << ievent << std::endl;
  DstClose();

  return EXIT_SUCCESS;
}

//_____________________________________________________________________
bool
dst::InitializeEvent( void )
{
  event.status   = 0;
  event.evnum = 0;
  event.nhittpc = 0; 
  event.nttpc = 0; 

  for(int i=0; i<MaxTPCTracks; ++i){
    event.nhit_track[i] =0;
    event.chisqr[i] =-9999.;
    event.helix_cx[i] =-9999.;
    event.helix_cy[i] =-9999.;
    event.helix_z0[i] =-9999.;
    event.helix_r[i] =-9999.;
    event.helix_dz[i] =-9999.;
    event.mom0_x[i] =-9999.;
    event.mom0_y[i] =-9999.;
    event.mom0_z[i] =-9999.;
    for(int j=0; j<MaxTPCnHits; ++j){
      event.hitlayer[i][j] =-999;    
      event.mom_x[i][j] =-9999.;    
      event.mom_y[i][j] =-9999.;    
      event.mom_z[i][j] =-9999.;   
      event.momg_x[i][j] =-9999.;    
      event.momg_y[i][j] =-9999.;    
      event.momg_z[i][j] =-9999.;   
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
bool
dst::DstOpen( std::vector<std::string> arg )
{
  int open_file = 0;
  int open_tree = 0;
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
dst::DstRead( int ievent )
{
  static const std::string func_name("["+class_name+"::"+__func__+"]");


  if( ievent%10000==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  HF1( 1, event.status++ );
  
  event.evnum = src.evnum;
  event.nhittpc = src.nhittpc;
 
  double u = src.pxPrm[0]/src.pzPrm[0];
  double v = src.pyPrm[0]/src.pzPrm[0];
  double cost = 1./std::sqrt(1.+u*u+v*v);
  double theta=std::acos(cost)*math::Rad2Deg();
  // if(theta>20.)
  //   return true;
  DCAnalyzer *DCAna = new DCAnalyzer();
  
  //with stable resolution
  // for(int it=0; it<src.nhittpc; ++it){
  //   src.x0tpc[it] += gRandom->Gaus(0,0.2);
  //   src.y0tpc[it] += gRandom->Gaus(0,0.5);
  //   src.z0tpc[it] += gRandom->Gaus(0,0.2);
  // }
  //for test
   // for(int it=0; it<src.nhittpc; ++it){
   //   src.xtpc[it] += gRandom->Gaus(0,0.1);
   //   src.ztpc[it] += gRandom->Gaus(0,0.1);
   // }

  //  std::cout<<"nhittpc"<<src.nhittpc<<std::endl;
  DCAna->DecodeTPCHitsGeant4(src.nhittpc, 
			     src.x0tpc, src.y0tpc, src.z0tpc, src.edeptpc);
   			     //src.xtpc, src.ytpc, src.ztpc, src.edeptpc);
  DCAna->TrackSearchTPC_Helix();
  
  int nttpc = DCAna->GetNTracksTPC_Helix();
  if( MaxHits<nttpc ){
    std::cout << "#W " << func_name << " " 
      	      << "too many nttpc " << nttpc << "/" << MaxHits << std::endl;
    nttpc = MaxHits;
  }
  event.nttpc = nttpc;
  for( int it=0; it<nttpc; ++it ){
    TPCLocalTrack_Helix *tp= DCAna->GetTrackTPC_Helix(it); 
    if(!tp) continue;
    int nh=tp->GetNHit();
    double chisqr    = tp->GetChiSquare();

    double cx=tp->Getcx(), cy=tp->Getcy();
    double z0=tp->Getz0(), r=tp->Getr();
    double dz = tp->Getdz();
    TVector3 mom0 = tp->GetMom0();

    event.nhit_track[it] = nh;
    event.chisqr[it] = chisqr;
    event.helix_cx[it] = cx;
    event.helix_cy[it] = cy;
    event.helix_z0[it] = z0;
    event.helix_r[it]  = r;
    event.helix_dz[it] = dz;
    event.mom0_x[it] = mom0.X();
    event.mom0_y[it] = mom0.Y();
    event.mom0_z[it] = mom0.Z();
    
    for( int ih=0; ih<nh; ++ih ){
      TPCLTrackHit *hit=tp->GetHit(ih);
      
      if(!hit) continue;
      int layerId = 0;
      layerId = hit->GetLayer();
      TVector3 hitpos = hit->GetLocalHitPos();
      TVector3 calpos = hit->GetLocalCalPos_Helix();
      double residual = hit->GetResidual();
      TVector3 res_vect = hit->GetResidualVect();
      TVector3 mom = hit->GetMomentum_Helix();      

      for( int ih2=0; ih2<src.nhittpc; ++ih2 ){
	TVector3 setpos(src.xtpc[ih2], src.ytpc[ih2], src.ztpc[ih2]);
	//TVector3 setpos(src.xtpc0[ih], src.ytpc0[ih], src.ztpc0[ih]);
	TVector3 d = setpos - hitpos;
	if(fabs(d.Mag()<0.1)){
	  event.momg_x[it][ih] = src.pxtpc[ih2];
	  event.momg_y[it][ih] = src.pytpc[ih2];
	  event.momg_z[it][ih] = src.pztpc[ih2];
	  break;
	}
      }

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

  // if( event.nhBh1<=0 ) return true;
  // HF1( 1, event.status++ );

  // if( event.nhBh2<=0 ) return true;
  // HF1( 1, event.status++ );

  // if( event.nhTof<=0 ) return true;
  // HF1( 1, event.status++ );

  // if( event.ntKurama<=0 ) return true;
  // if( event.ntKurama>MaxTPCHits )
  //   event.ntKurama = MaxTPCHits;

  //HF1( 1, event.status++ );

  delete DCAna;
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
  HB1( 1, "Status", 21, 0., 21. );


  HBTree( "ktpc_g", "tree of DstTPC_g" );
  
  tree->Branch("status", &event.status, "status/I" );  
  tree->Branch("evnum", &event.evnum, "evnum/I" );  
  tree->Branch("nhittpc",&event.nhittpc,"nhittpc/I");
  tree->Branch("nttpc",&event.nttpc,"nttpc/I");
  
  tree->Branch("nhit_track",event.nhit_track,"nhit_track[nttpc]/I");
  tree->Branch("chisqr",event.chisqr,"chisqr[nttpc]/D");
  tree->Branch("helix_cx",event.helix_cx,"helix_cx[nttpc]/D");
  tree->Branch("helix_cy",event.helix_cy,"helix_cy[nttpc]/D");
  tree->Branch("helix_z0",event.helix_z0,"helix_z0[nttpc]/D");
  tree->Branch("helix_r" ,event.helix_r, "helix_r[nttpc]/D");
  tree->Branch("helix_dz",event.helix_dz,"helix_dz[nttpc]/D");
  // Momentum at Y = 0
  tree->Branch("mom0_x",event.mom0_x,"mom0_x[nttpc]/D");
  tree->Branch("mom0_y",event.mom0_y,"mom0_y[nttpc]/D");
  tree->Branch("mom0_z",event.mom0_z,"mom0_z[nttpc]/D");
  
  tree->Branch("hitlayer",event.hitlayer,"hitlayer[nttpc][32]/I");
  tree->Branch("hitpos_x",event.hitpos_x,"hitpos_x[nttpc][32]/D");
  tree->Branch("hitpos_y",event.hitpos_y,"hitpos_y[nttpc][32]/D");
  tree->Branch("hitpos_z",event.hitpos_z,"hitpos_z[nttpc][32]/D");
  tree->Branch("calpos_x",event.calpos_x,"calpos_x[nttpc][32]/D");
  tree->Branch("calpos_y",event.calpos_y,"calpos_y[nttpc][32]/D");
  tree->Branch("calpos_z",event.calpos_z,"calpos_z[nttpc][32]/D");
  tree->Branch("mom_x",event.mom_x,"mom_x[nttpc][32]/D");
  tree->Branch("mom_y",event.mom_y,"mom_y[nttpc][32]/D");
  tree->Branch("mom_z",event.mom_z,"mom_z[nttpc][32]/D");
  tree->Branch("momg_x",event.momg_x,"momg_x[nttpc][32]/D");
  tree->Branch("momg_y",event.momg_y,"momg_y[nttpc][32]/D");
  tree->Branch("momg_z",event.momg_z,"momg_z[nttpc][32]/D");

  tree->Branch("residual",event.residual,"residual[nttpc][32]/D");
  tree->Branch("residual_x",event.residual_x,"residual_x[nttpc][32]/D");
  tree->Branch("residual_y",event.residual_y,"residual_y[nttpc][32]/D");
  tree->Branch("residual_z",event.residual_z,"residual_z[nttpc][32]/D");

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
  TTreeCont[kTPCGeant]->SetBranchStatus("pxtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pytpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pztpc", 1);
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
  TTreeCont[kTPCGeant]->SetBranchAddress("pxtpc", src.pxtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pytpc", src.pytpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pztpc", src.pztpc);
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
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")   &&
      InitializeParameter<UserParamMan>("USER") &&
	  InitializeParameter<HodoPHCMan>("HDPHC") ); 
}

//_____________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}

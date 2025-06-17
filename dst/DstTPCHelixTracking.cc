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
#include "TPCCluster.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#define TrackSearch 1
#define TrackCluster 1
#define TruncatedMean 0
#define TrigA 0 //if 1, TrigA is required
#define TrigB 0
#define TrigC 0
#define TrigD 0

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

  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Int_t> isBeam; // isBeam: 1 = Beam, 0 = Scat
  std::vector<Double_t> chisqr;
  std::vector<Double_t> helix_cx;
  std::vector<Double_t> helix_cy;
  std::vector<Double_t> helix_z0;
  std::vector<Double_t> helix_r;
  std::vector<Double_t> helix_dz;
  std::vector<Double_t> dE;
  std::vector<Double_t> dEdx; //reference dedx

  std::vector<Double_t> dEdx_0;
  std::vector<Double_t> dEdx_10;
  std::vector<Double_t> dEdx_20;
  std::vector<Double_t> dEdx_30;
  std::vector<Double_t> dEdx_40;
  std::vector<Double_t> dEdx_50;
  std::vector<Double_t> dEdx_60;
  std::vector<Double_t> dEdx_cor_0;
  std::vector<Double_t> dEdx_cor_10;
  std::vector<Double_t> dEdx_cor_20;
  std::vector<Double_t> dEdx_cor_30;
  std::vector<Double_t> dEdx_cor_40;
  std::vector<Double_t> dEdx_cor_50;
  std::vector<Double_t> dEdx_cor_60;

  std::vector<Double_t> dz_factor;
  std::vector<Double_t> mom0_x;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_y;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_z;//Helix momentum at Y = 0
  std::vector<Double_t> mom0;//Helix momentum at Y = 0

  std::vector<Int_t> charge;//Helix charge
  std::vector<Double_t> path;//Helix path

  // under dev
  std::vector<Int_t> pid;

  std::vector<std::vector<Double_t>> combi_id; //track number of combi
  std::vector<std::vector<Double_t>> closeDistTpc;
  std::vector<std::vector<Double_t>> vtxTpc;
  std::vector<std::vector<Double_t>> vtyTpc;
  std::vector<std::vector<Double_t>> vtzTpc;
  std::vector<std::vector<Double_t>> mom_vtx;//Helix momentum at vtx
  std::vector<std::vector<Double_t>> mom_vty;//Helix momentum at vtx
  std::vector<std::vector<Double_t>> mom_vtz;//Helix momentum at vtx

  std::vector<Double_t> M_Lambda;
  std::vector<Double_t> Lvtx;
  std::vector<Double_t> Lvty;
  std::vector<Double_t> Lvtz;
  std::vector<Double_t> LcloseDist;
  std::vector<Double_t> Mom_Lambda;
  std::vector<Double_t> Mom_Lambdax;
  std::vector<Double_t> Mom_Lambday;
  std::vector<Double_t> Mom_Lambdaz;

  std::vector<Double_t> M_Ks;
  std::vector<Double_t> Ksvtx;
  std::vector<Double_t> Ksvty;
  std::vector<Double_t> Ksvtz;
  std::vector<Double_t> KscloseDist;
  std::vector<Double_t> Mom_Ks;
  std::vector<Double_t> Mom_Ksx;
  std::vector<Double_t> Mom_Ksy;
  std::vector<Double_t> Mom_Ksz;

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
  std::vector<std::vector<Double_t>> helix_t;
  std::vector<std::vector<Double_t>> pathhit;
  std::vector<std::vector<Double_t>> pathhit_cor;
  std::vector<std::vector<Double_t>> theta_diff;
  std::vector<std::vector<Double_t>> track_cluster_de;
  std::vector<std::vector<Double_t>> track_cluster_size;
  std::vector<std::vector<Double_t>> track_cluster_mrow;
  std::vector<std::vector<Double_t>> track_cluster_de_center;
  std::vector<std::vector<Double_t>> track_cluster_x_center;
  std::vector<std::vector<Double_t>> track_cluster_y_center;
  std::vector<std::vector<Double_t>> track_cluster_z_center;
  std::vector<std::vector<Double_t>> track_cluster_row_center;

  void clear( void )
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
    cluster_de_center.clear();
    cluster_x_center.clear();
    cluster_y_center.clear();
    cluster_z_center.clear();
    cluster_row_center.clear();

    ntTpc = 0;
    nhtrack.clear();
    isBeam.clear();
    chisqr.clear();
    helix_cx.clear();
    helix_cy.clear();
    helix_z0.clear();
    helix_r.clear();
    helix_dz.clear();
    dE.clear();
    dEdx.clear();
#if TruncatedMean
    dEdx_0.clear();
    dEdx_10.clear();
    dEdx_20.clear();
    dEdx_30.clear();
    dEdx_40.clear();
    dEdx_50.clear();
    dEdx_60.clear();
    dEdx_cor_0.clear();
    dEdx_cor_10.clear();
    dEdx_cor_20.clear();
    dEdx_cor_30.clear();
    dEdx_cor_40.clear();
    dEdx_cor_50.clear();
    dEdx_cor_60.clear();
#endif
    dz_factor.clear();

    mom0_x.clear();
    mom0_y.clear();
    mom0_z.clear();
    mom0.clear();

    charge.clear();
    path.clear();

    combi_id.clear();
    mom_vtx.clear();
    mom_vty.clear();
    mom_vtz.clear();

    pid.clear();
    vtxTpc.clear();
    vtyTpc.clear();
    vtzTpc.clear();
    closeDistTpc.clear();

    M_Lambda.clear();
    Lvtx.clear();
    Lvty.clear();
    Lvtz.clear();
    LcloseDist.clear();
    Mom_Lambda.clear();
    Mom_Lambdax.clear();
    Mom_Lambday.clear();
    Mom_Lambdaz.clear();

    M_Ks.clear();
    Ksvtx.clear();
    Ksvty.clear();
    Ksvtz.clear();
    KscloseDist.clear();
    Mom_Ks.clear();
    Mom_Ksx.clear();
    Mom_Ksy.clear();
    Mom_Ksz.clear();

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
    helix_t.clear();

    pathhit.clear();
    pathhit_cor.clear();
    theta_diff.clear();
    track_cluster_de.clear();
    track_cluster_size.clear();
    track_cluster_mrow.clear();
    track_cluster_de_center.clear();
    track_cluster_x_center.clear();
    track_cluster_y_center.clear();
    track_cluster_z_center.clear();
    track_cluster_row_center.clear();
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
  TTreeReaderValue<std::vector<Double_t>>* cdeTpc;    // cdE
  TTreeReaderValue<std::vector<Double_t>>* tTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* ctTpc;     // time
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
    TPCGainHid = 100000,
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
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries( TTreeCont );
  if (max_loop > 0) nevent = skip + max_loop;

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
  event.clkTpc = **src.clkTpc;

#if TrigA
  if(event.trigflag[20]<0) return true;
#endif
#if TrigB
  if(event.trigflag[21]<0) return true;
#endif
#if TrigC
  if(event.trigflag[22]<0) return true;
#endif
#if TrigD
  if(event.trigflag[23]<0) return true;
#endif

  HF1( 1, event.status++ );

  if( **src.nhTpc == 0 )
    return true;

  HF1( 1, event.status++ );

  if(event.clkTpc.size() != 1){
    std::cerr << "something is wrong" << std::endl;
    return true;
  }
  Double_t clock = event.clkTpc.at(0);
  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, clock);

  const Double_t VertexScanRange = gUser.GetParameter("VertexScanRange"); //mm

  HF1(1, event.status++);
  Int_t nh_Tpc = 0;
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    auto hc = TPCAna.GetTPCHC( layer );
    for( const auto& hit : hc ){
      if( !hit || !hit->IsGood() )
        continue;
      Double_t x = hit->GetX();
      Double_t y = hit->GetY();
      Double_t z = hit->GetZ();
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
      ++nh_Tpc;
    }
  }
  event.nhTpc = nh_Tpc;
  HF1(1, event.status++);

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
      event.cluster_de_center.push_back(centerDe);
      event.cluster_x_center.push_back(centerPos.X());
      event.cluster_y_center.push_back(centerPos.Y());
      event.cluster_z_center.push_back(centerPos.Z());
      event.cluster_row_center.push_back(centerRow);

      ++nclTpc;
    }
  }
  event.nclTpc = nclTpc;
  HF1( 1, event.status++ );

#if TrackSearch
  TPCAna.TrackSearchTPCHelix();
#endif

  Int_t ntTpc = TPCAna.GetNTracksTPCHelix();
  event.ntTpc = ntTpc;
  HF1( 10, ntTpc );
  if( event.ntTpc == 0 )
    return true;

  HF1( 1, event.status++ );

  event.nhtrack.resize( ntTpc );
  event.isBeam.resize( ntTpc );
  event.chisqr.resize( ntTpc );
  event.helix_cx.resize( ntTpc );
  event.helix_cy.resize( ntTpc );
  event.helix_z0.resize( ntTpc );
  event.helix_r.resize( ntTpc );
  event.helix_dz.resize( ntTpc );
  event.mom0_x.resize( ntTpc );
  event.mom0_y.resize( ntTpc );
  event.mom0_z.resize( ntTpc );
  event.mom0.resize( ntTpc );

  event.dE.resize( ntTpc );
  event.dEdx.resize( ntTpc );
#if TruncatedMean
  event.dEdx_0.resize( ntTpc );
  event.dEdx_10.resize( ntTpc );
  event.dEdx_20.resize( ntTpc );
  event.dEdx_30.resize( ntTpc );
  event.dEdx_40.resize( ntTpc );
  event.dEdx_50.resize( ntTpc );
  event.dEdx_60.resize( ntTpc );
  event.dEdx_cor_0.resize( ntTpc );
  event.dEdx_cor_10.resize( ntTpc );
  event.dEdx_cor_20.resize( ntTpc );
  event.dEdx_cor_30.resize( ntTpc );
  event.dEdx_cor_40.resize( ntTpc );
  event.dEdx_cor_50.resize( ntTpc );
  event.dEdx_cor_60.resize( ntTpc );
#endif
  event.dz_factor.resize( ntTpc );
  event.charge.resize( ntTpc );
  event.path.resize( ntTpc );

  event.pid.resize( ntTpc );
  event.combi_id.resize( ntTpc );
  event.closeDistTpc.resize( ntTpc );
  event.vtxTpc.resize( ntTpc );
  event.vtyTpc.resize( ntTpc );
  event.vtzTpc.resize( ntTpc );
  event.mom_vtx.resize( ntTpc );
  event.mom_vty.resize( ntTpc );
  event.mom_vtz.resize( ntTpc );

  event.hitlayer.resize( ntTpc );
  event.hitpos_x.resize( ntTpc );
  event.hitpos_y.resize( ntTpc );
  event.hitpos_z.resize( ntTpc );
  event.calpos_x.resize( ntTpc );
  event.calpos_y.resize( ntTpc );
  event.calpos_z.resize( ntTpc );
  event.residual.resize( ntTpc );
  event.residual_x.resize( ntTpc );
  event.residual_y.resize( ntTpc );
  event.residual_z.resize( ntTpc );
  event.helix_t.resize( ntTpc );
  event.pathhit.resize(ntTpc);
  event.pathhit_cor.resize(ntTpc);
  event.theta_diff.resize(ntTpc);

  event.track_cluster_de.resize(ntTpc);
  event.track_cluster_size.resize(ntTpc);
  event.track_cluster_mrow.resize(ntTpc);
  event.track_cluster_de_center.resize(ntTpc);
  event.track_cluster_x_center.resize(ntTpc);
  event.track_cluster_y_center.resize(ntTpc);
  event.track_cluster_z_center.resize(ntTpc);
  event.track_cluster_row_center.resize(ntTpc);

  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( it );
    if( !tp ) continue;
    Int_t nh = tp->GetNHit();
    Double_t chisqr = tp->GetChiSquare();
    Double_t helix_cx = tp->Getcx(), helix_cy = tp->Getcy();
    Double_t helix_z0 = tp->Getz0(), helix_r = tp->Getr();
    Double_t helix_dz = tp->Getdz();
    TVector3 Mom0 = tp->GetMom0();
    Int_t isbeam = tp->GetIsBeam();

    HF1(11, nh);
    HF1(12, chisqr);

    event.nhtrack[it] = nh;
    event.isBeam[it] = isbeam;
    event.chisqr[it] = chisqr;
    event.helix_cx[it] = helix_cx;
    event.helix_cy[it] = helix_cy;
    event.helix_z0[it] = helix_z0;
    event.helix_r[it] = helix_r ;
    event.helix_dz[it] = helix_dz;
    event.mom0_x[it] = Mom0.x();
    event.mom0_y[it] = Mom0.y();
    event.mom0_z[it] = Mom0.z();
    event.mom0[it] = Mom0.Mag();

    event.combi_id[it].resize( ntTpc );
    event.closeDistTpc[it].resize( ntTpc );
    event.vtxTpc[it].resize( ntTpc );
    event.vtyTpc[it].resize( ntTpc );
    event.vtzTpc[it].resize( ntTpc );
    event.mom_vtx[it].resize( ntTpc );
    event.mom_vty[it].resize( ntTpc );
    event.mom_vtz[it].resize( ntTpc );

    event.hitlayer[it].resize( nh );
    event.hitpos_x[it].resize( nh );
    event.hitpos_y[it].resize( nh );
    event.hitpos_z[it].resize( nh );
    event.calpos_x[it].resize( nh );
    event.calpos_y[it].resize( nh );
    event.calpos_z[it].resize( nh );
    event.residual[it].resize( nh );
    event.residual_x[it].resize( nh );
    event.residual_y[it].resize( nh );
    event.residual_z[it].resize( nh );
    event.helix_t[it].resize( nh );
    event.pathhit[it].resize(nh);
    event.pathhit_cor[it].resize(nh);
    event.theta_diff[it].resize(nh);

    event.track_cluster_de[it].resize(nh);
    event.track_cluster_size[it].resize(nh);
    event.track_cluster_mrow[it].resize(nh);
    event.track_cluster_de_center[it].resize(nh);
    event.track_cluster_x_center[it].resize(nh);
    event.track_cluster_y_center[it].resize(nh);
    event.track_cluster_z_center[it].resize(nh);
    event.track_cluster_row_center[it].resize(nh);

    Double_t par1[5]={helix_cx, helix_cy, helix_z0,
		    helix_r, helix_dz};

    Double_t scantheta1 = VertexScanRange/par1[3]; //mm -> rad.
    Double_t range_theta1[2] = {tp -> GetMint() - scantheta1,
				tp -> GetMaxt() + scantheta1};

    for( Int_t it2=0; it2<ntTpc; ++it2 ){
      if(it2==it) continue;
      TPCLocalTrackHelix *tp2 = TPCAna.GetTrackTPCHelix( it2 );
      if( !tp2 ) continue;
      Double_t helix_cx2 = tp2->Getcx(), helix_cy2 = tp2->Getcy();
      Double_t helix_z02 = tp2->Getz0(), helix_r2 = tp2->Getr();
      Double_t helix_dz2 = tp2->Getdz();

      Double_t par2[5]={helix_cx2, helix_cy2, helix_z02,
			helix_r2, helix_dz2};

      Double_t scantheta2 = VertexScanRange/par2[3]; //mm -> rad.
      Double_t range_theta2[2] = {tp -> GetMint() - scantheta2,
				  tp -> GetMaxt() + scantheta2};

      double closeDistTpc, t1, t2;
      TVector3 vert
	= Kinematics::VertexPointHelix(par1, par2,
				       range_theta1[0], range_theta1[1],
				       range_theta2[0], range_theta2[1],
				       t1, t2, closeDistTpc);

      event.combi_id[it][it2] = (double)it2;
      event.closeDistTpc[it][it2] = closeDistTpc;
      event.vtxTpc[it][it2] = vert.x();
      event.vtyTpc[it][it2] = vert.y();
      event.vtzTpc[it][it2] = vert.z();

      TVector3 mom_vtx = tp->CalcHelixMom(par1, t1);
      event.mom_vtx[it][it2] = mom_vtx.x();
      event.mom_vty[it][it2] = mom_vtx.y();
      event.mom_vtz[it][it2] = mom_vtx.z();
    }

    Double_t min_t = 10000.; Double_t max_t = -10000.;
    Double_t min_layer_t = 0., max_layer_t = 0.;
    Int_t min_layer = 33, max_layer = -1;
    Double_t de=0.;
    std::vector<Double_t> dEdx_vect; std::vector<Double_t> dEdx_cor_vect;
    for( int ih=0; ih<nh; ++ih ){
      TPCLTrackHit *hit = tp->GetHit( ih );
      if( !hit ) continue;
      HF1( 2, hit->GetHoughDist());
      HF1( 3, hit->GetHoughDistY());

      Int_t layer = hit->GetLayer();
      const TVector3& hitpos = hit->GetLocalHitPos();
      const TVector3& calpos = hit->GetLocalCalPosHelix();
      const TVector3& resi_vect = hit->GetResidualVect();

      HF1(13, layer);

      TPCHit *clhit = hit -> GetHit();
      TPCCluster *cl = clhit -> GetParentCluster();
      Int_t clsize = cl->GetClusterSize();
      Double_t clde = cl->GetDe();
      Double_t mrow = cl->MeanRow(); // same
      TPCHit* centerHit = cl->GetCenterHit();
      const TVector3& centerPos = centerHit->GetPosition();
      Double_t centerDe = centerHit->GetCDe();
      Int_t centerRow = centerHit->GetRow();

#if TrackCluster
      if(ntTpc==1&&nh>=20){
	HF1(TPCClHid, clsize);
	HF1(TPCClHid+(layer+1)*1000, clsize);
	HF1(TPCClHid+1, clde);
	HF1(TPCClHid+(layer+1)*1000+1, clde);
	const TPCHitContainer& hc = cl -> GetHitContainer();
	for(const auto& hits : hc){
	  if(!hits || !hits->IsGood() || !hit->IsGoodForTracking()) continue;
	  const TVector3& pos = hits->GetPosition();
	  Double_t de = hits->GetCDe();
	  Double_t transDist = TranseverseDistance(hitpos.x(), hitpos.z(), pos.x(), pos.z());
	  Double_t ratio = de/clde;
	  HF2(TPCClHid+2, transDist, ratio);
	  HF2(TPCClHid+(layer+1)*1000+2, transDist, ratio);
	}
      }
#endif
      event.track_cluster_de[it][ih] = clde;
      event.track_cluster_size[it][ih] = clsize;
      event.track_cluster_mrow[it][ih] = mrow;
      event.track_cluster_de_center[it][ih] = centerDe;
      event.track_cluster_x_center[it][ih] = centerPos.X();
      event.track_cluster_y_center[it][ih] = centerPos.Y();
      event.track_cluster_z_center[it][ih] = centerPos.Z();
      event.track_cluster_row_center[it][ih] = centerRow;

      Double_t padTheta = tpc::getTheta(layer, mrow)*acos(-1)/180.;
      Double_t t_cal = hit->GetTheta();
      Double_t thetaDiff = t_cal - padTheta;
      event.theta_diff[it][ih] = thetaDiff;
      Double_t pathHit = tpc::padParameter[layer][5];
      event.pathhit[it][ih] = pathHit;

      //Approximation of pathHit correction
      Double_t pathHit_cor = (pathHit/fabs(cos(thetaDiff)))*sqrt(1.+(pow(helix_dz,2)));
      event.pathhit_cor[it][ih] = pathHit_cor;

      //Charge
      if(min_t>t_cal) min_t = t_cal;
      if(max_t<t_cal) max_t = t_cal;
      if(layer<min_layer){
	min_layer = layer;
	min_layer_t = t_cal;
      }
      if(layer>max_layer){
	max_layer = layer;
	max_layer_t = t_cal;
      }

      de += clde;
      Double_t dEdx = clde/pathHit;
      Double_t dEdx_cor = clde/pathHit_cor;
      dEdx_vect.push_back(dEdx);
      dEdx_cor_vect.push_back(dEdx_cor);

      Double_t residual = hit->GetResidual();
      event.hitlayer[it][ih] = (double)layer;
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
      event.helix_t[it][ih] = t_cal;
    }
    if(min_layer_t<max_layer_t) event.charge[it] = 1;
    else event.charge[it] = -1;

    Double_t pathlen = (max_t - min_t)*sqrt(helix_r*helix_r*(1.+helix_dz*helix_dz));
    event.path[it] = pathlen;
    event.dE[it] = de;

    std::sort(dEdx_vect.begin(), dEdx_vect.end());
    std::sort(dEdx_cor_vect.begin(), dEdx_cor_vect.end());
    int n_truncated = (int)(dEdx_vect.size()*truncatedMean);
    for( int ih=0; ih<dEdx_vect.size(); ++ih ){
      if(ih<n_truncated) event.dEdx[it]+=dEdx_cor_vect[ih];
    }
    event.dEdx[it]/=(double)n_truncated;

#if TruncatedMean
    event.dEdx_0[it] = TMath::Mean(dEdx_vect.size(),dEdx_vect.data());
    event.dEdx_cor_0[it] = TMath::Mean(dEdx_cor_vect.size(),dEdx_cor_vect.data());
    event.dEdx_10[it]=0.;
    event.dEdx_20[it]=0.;
    event.dEdx_30[it]=0.;
    event.dEdx_40[it]=0.;
    event.dEdx_50[it]=0.;
    event.dEdx_60[it]=0.;
    event.dEdx_cor_10[it]=0.;
    event.dEdx_cor_20[it]=0.;
    event.dEdx_cor_30[it]=0.;
    event.dEdx_cor_40[it]=0.;
    event.dEdx_cor_50[it]=0.;
    event.dEdx_cor_60[it]=0.;

    int n_10 = (int)(dEdx_vect.size()*0.9);
    int n_20 = (int)(dEdx_vect.size()*0.8);
    int n_30 = (int)(dEdx_vect.size()*0.7);
    int n_40 = (int)(dEdx_vect.size()*0.6);
    int n_50 = (int)(dEdx_vect.size()*0.5);
    int n_60 = (int)(dEdx_vect.size()*0.4);
    for( int ih=0; ih<dEdx_vect.size(); ++ih ){
      if(ih<n_10){
	event.dEdx_10[it]+=dEdx_vect[ih];
	event.dEdx_cor_10[it]+=dEdx_cor_vect[ih];
      }
      if(ih<n_20){
	event.dEdx_20[it]+=dEdx_vect[ih];
	event.dEdx_cor_20[it] +=dEdx_cor_vect[ih];
      }
      if(ih<n_30){
	event.dEdx_30[it]+=dEdx_vect[ih];
	event.dEdx_cor_30[it]+=dEdx_cor_vect[ih];
      }
      if(ih<n_40){
	event.dEdx_40[it]+=dEdx_vect[ih];
	event.dEdx_cor_40[it]+=dEdx_cor_vect[ih];
      }
      if(ih<n_50){
	event.dEdx_50[it]+=dEdx_vect[ih];
	event.dEdx_cor_50[it]+=dEdx_cor_vect[ih];
      }
      if(ih<n_60){
	event.dEdx_60[it]+=dEdx_vect[ih];
	event.dEdx_cor_60[it]+=dEdx_cor_vect[ih];
      }
    }
    event.dEdx_10[it]/=(double)n_10;
    event.dEdx_20[it]/=(double)n_20;
    event.dEdx_30[it]/=(double)n_30;
    event.dEdx_40[it]/=(double)n_40;
    event.dEdx_50[it]/=(double)n_50;
    event.dEdx_60[it]/=(double)n_60;
    event.dEdx_cor_10[it]/=(double)n_10;
    event.dEdx_cor_20[it]/=(double)n_20;
    event.dEdx_cor_30[it]/=(double)n_30;
    event.dEdx_cor_40[it]/=(double)n_40;
    event.dEdx_cor_50[it]/=(double)n_50;
    event.dEdx_cor_60[it]/=(double)n_60;
#endif
    event.dz_factor[it] = sqrt(1.+(pow(helix_dz,2)));
    HF1(14, event.mom0[it]);

    int particleID= Kinematics::HypTPCdEdxPID(event.dEdx[it], event.mom0[it]*event.charge[it]);
    event.pid[it]=particleID;
  }

#if 0
  // static const auto KaonMass    = pdg::KaonMass();
  static const auto PionMass    = pdg::PionMass();
  static const auto ProtonMass  = pdg::ProtonMass();

  // rough estimation for Lambda and Ks event
  if(ntTpc>1){
    for( Int_t it=0; it<ntTpc-1; ++it ){
      for( Int_t it2=it+1; it2<ntTpc; ++it2 ){
	if(event.isBeam[it]==0&&event.isBeam[it2]==0
	   &&event.combi_id[it]==it2&&event.combi_id[it2]==it
	   &&event.charge[it]==-1*event.charge[it2]){
	  if(event.charge[it]==1){
	    TVector3 mom_pos(event.mom_cor_vtx[it], event.mom_cor_vty[it], event.mom_cor_vtz[it]);
	    TLorentzVector Lp(mom_pos, std::sqrt(ProtonMass*ProtonMass+mom_pos.Mag2()));
	    TLorentzVector Lpip(mom_pos, std::sqrt(PionMass*PionMass+mom_pos.Mag2()));
	    TVector3 mom_neg(event.mom_cor_vtx[it2], event.mom_cor_vty[it2], event.mom_cor_vtz[it2]);
	    TLorentzVector Lpi(mom_neg, std::sqrt(PionMass*PionMass+mom_neg.Mag2()));
	    TLorentzVector LLambda = Lp + Lpi;
	    TLorentzVector LKs = Lpip + Lpi;
	    //	    if(event.dEdx_20_cor[it]>event.dEdx_20_cor[it2]*1.5){
	    if(event.pid[it]==2&&event.pid[it2]==-1){
	      event.M_Lambda.push_back(LLambda.M());
	      event.Lvtx.push_back(event.vtx[it]);
	      event.Lvty.push_back(event.vty[it]);
	      event.Lvtz.push_back(event.vtz[it]);
	      event.LcloseDist.push_back(event.closeDist[it]);
	      event.Mom_Lambda.push_back(LLambda.P());
	      event.Mom_Lambdax.push_back(LLambda.Px());
	      event.Mom_Lambday.push_back(LLambda.Py());
	      event.Mom_Lambdaz.push_back(LLambda.Pz());
	    }
	    if(event.pid[it]==1&&event.pid[it2]==-1){
	      event.M_Ks.push_back(LKs.M());
	      event.Ksvtx.push_back(event.vtx[it]);
	      event.Ksvty.push_back(event.vty[it]);
	      event.Ksvtz.push_back(event.vtz[it]);
	      event.KscloseDist.push_back(event.closeDist[it]);
	      event.Mom_Ks.push_back(LKs.P());
	      event.Mom_Ksx.push_back(LKs.Px());
	      event.Mom_Ksx.push_back(LKs.Py());
	      event.Mom_Ksx.push_back(LKs.Pz());
	    }
	  }
	  else{
	    TVector3 mom_pos(event.mom_cor_vtx[it2], event.mom_cor_vty[it2], event.mom_cor_vtz[it2]);
	    TLorentzVector Lp(mom_pos, std::sqrt(ProtonMass*ProtonMass+mom_pos.Mag2()));
	    TLorentzVector Lpip(mom_pos, std::sqrt(PionMass*PionMass+mom_pos.Mag2()));
	    TVector3 mom_neg(event.mom_cor_vtx[it], event.mom_cor_vty[it], event.mom_cor_vtz[it]);
	    TLorentzVector Lpi(mom_neg, std::sqrt(PionMass*PionMass+mom_neg.Mag2()));
	    TLorentzVector LLambda = Lp + Lpi;
	    TLorentzVector LKs = Lpip + Lpi;
	    //	    if(event.dEdx_20_cor[it2]>event.dEdx_20_cor[it]*1.5){
	    if(event.pid[it2]==2&&event.pid[it]==-1){
	      event.M_Lambda.push_back(LLambda.M());
	      event.Lvtx.push_back(event.mom_vtx[it]);
	      event.Lvty.push_back(event.mom_vty[it]);
	      event.Lvtz.push_back(event.mom_vtz[it]);
	      event.LcloseDist.push_back(event.closeDist[it]);
	      event.Mom_Lambda.push_back(LLambda.P());
	      event.Mom_Lambdax.push_back(LLambda.Px());
	      event.Mom_Lambday.push_back(LLambda.Py());
	      event.Mom_Lambdaz.push_back(LLambda.Pz());
	    }
	    if(event.pid[it2]==1&&event.pid[it]==-1){
	      event.M_Ks.push_back(LKs.M());
	      event.Ksvtx.push_back(event.mom_vtx[it]);
	      event.Ksvty.push_back(event.mom_vty[it]);
	      event.Ksvtz.push_back(event.mom_vtz[it]);
	      event.KscloseDist.push_back(event.closeDist[it]);
	      event.Mom_Ks.push_back(LKs.P());
	      event.Mom_Ksx.push_back(LKs.Px());
	      event.Mom_Ksy.push_back(LKs.Py());
	      event.Mom_Ksz.push_back(LKs.Pz());
	    }
	  }
	}
      }
    }
  }
#endif
  HF1( 1, event.status++ );

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

  HB1(1, "Status", 21, 0., 21. );
  HB1(2, "HoughDist", 500, 0., 50.);
  HB1(3, "HoughDistY", 500, 0., 50.);
  HB1(10, "#Tracks TPC", 40, 0., 40. );
  HB1(11, "#Hits of Track TPC", 50, 0., 50.);
  HB1(12, "Chisqr TPC", 500, 0., 500.);
  HB1(13, "LayerId TPC", 35, 0., 35.);
  HB1(14, "mom0", 1000, 0., 2.5);

#if TrackCluster
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
#endif

  HBTree( "tpc", "tree of DstTPCTracking" );

  tree->Branch( "status", &event.status );
  tree->Branch( "runnum", &event.runnum );
  tree->Branch( "evnum", &event.evnum );
  tree->Branch( "trigpat", &event.trigpat );
  tree->Branch( "trigflag", &event.trigflag );
  tree->Branch( "clkTpc", &event.clkTpc);

  tree->Branch( "nhTpc", &event.nhTpc );
  tree->Branch( "raw_hitpos_x", &event.raw_hitpos_x );
  tree->Branch( "raw_hitpos_y", &event.raw_hitpos_y );
  tree->Branch( "raw_hitpos_z", &event.raw_hitpos_z );
  tree->Branch( "raw_de", &event.raw_de );
  tree->Branch( "raw_padid", &event.raw_padid );
  tree->Branch( "raw_layer", &event.raw_layer );
  tree->Branch( "raw_row", &event.raw_row );

  tree->Branch( "nclTpc", &event.nclTpc );
  tree->Branch( "cluster_x", &event.cluster_x );
  tree->Branch( "cluster_y", &event.cluster_y );
  tree->Branch( "cluster_z", &event.cluster_z );
  tree->Branch( "cluster_de", &event.cluster_de );
  tree->Branch( "cluster_size", &event.cluster_size );
  tree->Branch( "cluster_layer", &event.cluster_layer );
  tree->Branch( "cluster_row_center", &event.cluster_row_center );
  tree->Branch( "cluster_mrow", &event.cluster_mrow );
  tree->Branch( "cluster_de_center", &event.cluster_de_center );
  tree->Branch( "cluster_x_center", &event.cluster_x_center );
  tree->Branch( "cluster_y_center", &event.cluster_y_center );
  tree->Branch( "cluster_z_center", &event.cluster_z_center );

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "isBeam", &event.isBeam );
  tree->Branch( "chisqr", &event.chisqr );
  tree->Branch( "helix_cx", &event.helix_cx );
  tree->Branch( "helix_cy", &event.helix_cy );
  tree->Branch( "helix_z0", &event.helix_z0 );
  tree->Branch( "helix_r", &event.helix_r );
  tree->Branch( "helix_dz", &event.helix_dz );
  tree->Branch( "mom0_x", &event.mom0_x );
  tree->Branch( "mom0_y", &event.mom0_y );
  tree->Branch( "mom0_z", &event.mom0_z );
  tree->Branch( "mom0", &event.mom0 );
  tree->Branch( "dE", &event.dE );
  tree->Branch( "dEdx", &event.dEdx );

#if TruncatedMean
  tree->Branch( "dEdx_0", &event.dEdx_0 );
  tree->Branch( "dEdx_10", &event.dEdx_10 );
  tree->Branch( "dEdx_20", &event.dEdx_20 );
  tree->Branch( "dEdx_30", &event.dEdx_30 );
  tree->Branch( "dEdx_40", &event.dEdx_40 );
  tree->Branch( "dEdx_50", &event.dEdx_50 );
  tree->Branch( "dEdx_60", &event.dEdx_60 );
  tree->Branch( "dEdx_cor_0", &event.dEdx_cor_0 );
  tree->Branch( "dEdx_cor_10", &event.dEdx_cor_10 );
  tree->Branch( "dEdx_cor_20", &event.dEdx_cor_20 );
  tree->Branch( "dEdx_cor_30", &event.dEdx_cor_30 );
  tree->Branch( "dEdx_cor_40", &event.dEdx_cor_40 );
  tree->Branch( "dEdx_cor_50", &event.dEdx_cor_50 );
  tree->Branch( "dEdx_cor_60", &event.dEdx_cor_60 );
#endif
  tree->Branch( "dz_factor", &event.dz_factor );
  tree->Branch( "charge", &event.charge );
  tree->Branch( "path", &event.path );

  tree->Branch( "pid", &event.pid );
  tree->Branch( "combi_id", &event.combi_id );
  tree->Branch( "closeDistTpc", &event.closeDistTpc );
  tree->Branch( "vtxTpc", &event.vtxTpc );
  tree->Branch( "vtyTpc", &event.vtyTpc );
  tree->Branch( "vtzTpc", &event.vtzTpc );
  tree->Branch( "mom_vtx", &event.mom_vtx );
  tree->Branch( "mom_vty", &event.mom_vty );
  tree->Branch( "mom_vtz", &event.mom_vtz );
#if 0
  tree->Branch( "M_Lambda", &event.M_Lambda );
  tree->Branch( "Lvtx", &event.Lvtx );
  tree->Branch( "Lvty", &event.Lvty );
  tree->Branch( "Lvtz", &event.Lvtz );
  tree->Branch( "LcloseDist", &event.LcloseDist );
  tree->Branch( "Mom_Lambda", &event.Mom_Lambda );
  tree->Branch( "Mom_Lambdax", &event.Mom_Lambdax );
  tree->Branch( "Mom_Lambday", &event.Mom_Lambday );
  tree->Branch( "Mom_Lambdaz", &event.Mom_Lambdaz );
  tree->Branch( "M_Ks", &event.M_Ks );
  tree->Branch( "Ksvtx", &event.Ksvtx );
  tree->Branch( "Ksvty", &event.Ksvty );
  tree->Branch( "Ksvtz", &event.Ksvtz );
  tree->Branch( "KscloseDist", &event.KscloseDist );
  tree->Branch( "Mom_Ks", &event.Mom_Ks );
  tree->Branch( "Mom_Ksx", &event.Mom_Ksx );
  tree->Branch( "Mom_Ksy", &event.Mom_Ksy );
  tree->Branch( "Mom_Ksz", &event.Mom_Ksz );
#endif
  tree->Branch( "hitlayer", &event.hitlayer );
  tree->Branch( "hitpos_x", &event.hitpos_x );
  tree->Branch( "hitpos_y", &event.hitpos_y );
  tree->Branch( "hitpos_z", &event.hitpos_z );
  tree->Branch( "calpos_x", &event.calpos_x );
  tree->Branch( "calpos_y", &event.calpos_y );
  tree->Branch( "calpos_z", &event.calpos_z );
  tree->Branch( "residual", &event.residual );
  tree->Branch( "residual_x", &event.residual_x );
  tree->Branch( "residual_y", &event.residual_y );
  tree->Branch( "residual_z", &event.residual_z );
  tree->Branch( "helix_t", &event.helix_t );
  tree->Branch( "theta_diff", &event.theta_diff);
  tree->Branch( "pathhit", &event.pathhit);
  tree->Branch( "pathhit_cor", &event.pathhit_cor);
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_size", &event.track_cluster_size);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch( "track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch( "track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch( "track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch( "track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch( "track_cluster_row_center", &event.track_cluster_row_center);

  TTreeReaderCont[kTpcHit] = new TTreeReader( "tpc", TFileCont[kTpcHit] );
  const auto& reader = TTreeReaderCont[kTpcHit];
  src.runnum = new TTreeReaderValue<Int_t>( *reader, "runnum" );
  src.evnum = new TTreeReaderValue<Int_t>( *reader, "evnum" );
  src.trigpat = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trigpat" );
  src.trigflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trigflag" );
  src.npadTpc = new TTreeReaderValue<Int_t>( *reader, "npadTpc" );
  src.nhTpc = new TTreeReaderValue<Int_t>( *reader, "nhTpc" );
  src.layerTpc = new TTreeReaderValue<std::vector<Int_t>>( *reader, "layerTpc" );
  src.rowTpc = new TTreeReaderValue<std::vector<Int_t>>( *reader, "rowTpc" );
  src.padTpc = new TTreeReaderValue<std::vector<Int_t>>( *reader, "padTpc" );
  src.pedTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pedTpc" );
  src.rmsTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "rmsTpc" );
  src.deTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "deTpc" );
  src.cdeTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cdeTpc" );
  src.tTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "tTpc" );
  src.ctTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ctTpc" );
  src.chisqrTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrTpc" );
  src.clkTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "clkTpc");

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")   &&
      InitializeParameter<TPCParamMan>("TPCPRM") &&
      InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
      InitializeParameter<UserParamMan>("USER") &&
      InitializeParameter<HodoPHCMan>("HDPHC") );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess( void )
{
  return true;
}

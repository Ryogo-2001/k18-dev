// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <Math/ProbFunc.h>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "FieldMan.hh"
#include "DatabasePDG.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "DCAnalyzer.hh"
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
#include "TPCLocalTrackHelix.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"
#include "FourVectorFitter.hh"

#include "HypTPCFitter.hh"
#include "HypTPCTask.hh"

#define SaveTPCK18 0
#define MMcut_Xi 1
#define SkipGenfit 0
#define WithKurama 0
#define DoKinematicFitLdXi 1
namespace
{
using namespace root;
using namespace std;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
const auto& gTPC  = TPCParamMan::GetInstance();
const Int_t MaxTPCHits = 10000;

//For GenFit Setting
const Bool_t Const_field = false; //Must be false for linear tracking
const Int_t verbosity = -1;//0~3;
//const Int_t verbosity = 1;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");

//const Double_t l_masscut = 0.15; //ref
//const Double_t l_masscut = 0.05;
const Double_t xi_masscut = 0.15;
//const Double_t lambda_masscut = 0.1; //ref
//const Double_t lambda_masscut = 0.05;
const Double_t lambda_masscut = 0.1;
const Double_t p_vtx_distcut = 300;
const Double_t pi_vtx_distcut = 300;
const Double_t pi2_vtx_distcut = 300;
//const Double_t l_vtx_distcut = 100;
//const Double_t ppi_distcut = 10.;
//const Double_t lpi_distcut = 10.;
const Double_t ppi_distcut = 10.;//ref
const Double_t lpi_distcut = 10.;//ref
const Double_t vtx_scan_range = 150.;
//const Double_t vtx_scan_range = 50.;
//const Double_t vtx_scan_rangeInsideL = 50.;
//const Double_t vtx_scan_rangeInsidePi = 50.;
//const Double_t vtx_scan_rangeInsideL = 200.;
//const Double_t vtx_scan_rangeInsideL = 300.;
const Double_t vtx_scan_rangeInsideL = 50.;
const Double_t vtx_scan_rangeInsidePi = 50.;
const Double_t GFppi_distcut = 15.;
const Double_t GFlpi_distcut = 15.;

const Double_t& HS_field_0 = ConfMan::Get<Double_t>("HSFLDCALIB");
const Double_t& HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
const Double_t& HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");
int TPCToG4TrackID(std::vector<TVector3>TPCHit, int nhG4,int* tidG4, double* xG4,double* yG4,double* zG4 ,int& nhits){
  std::vector<TVector3> G4Hits;
  for(int ih=0;ih<nhG4;++ih){
    TVector3 G4Hit(xG4[ih],yG4[ih],zG4[ih]);
    G4Hits.push_back(G4Hit);
  }
  int MaxTracks = 1000;
  TH1I Counter("counter","counter",MaxTracks,0,MaxTracks);
  for(auto hit:TPCHit){
    double dl = 5000;
    int G4ID = -1;
    for(int ih=0;ih<nhG4;++ih){
      auto G4Hit = G4Hits.at(ih);
      double dist = (G4Hit - hit).Mag();
      if(dist < dl){
        dl = dist;
        G4ID = tidG4[ih];
      }
    }
    Counter.Fill(G4ID);
  }
  nhits = Counter.GetMaximum();
  int G4id = Counter.GetMaximumBin()-1;
  return G4id;
}
TVector3 GetG4Mom(TVector3 TPCHit, vector<TVector3> G4Hits,vector<TVector3>G4Moms){
  int nh = G4Hits.size();
  double dl = 5000;
  TVector3 mom;
  for(int ih=0;ih<nh;++ih){
    auto G4Hit = G4Hits.at(ih);
    double dist = (G4Hit - TPCHit).Mag();
    if(dist < dl){
      dl = dist;
      mom = G4Moms.at(ih);
    }
  }
  return mom;
}
int CountHits(int id, double* posx,double* posz,int* tidG4, int nh){
  int count = 0;
  for(int ih=0;ih<nh;++ih){
    double x = posx[ih],z=posz[ih];
    int pad = tpc::findPadID(z,x);
    int layer = tpc::getLayerID(pad);
    int row = tpc::getRowID(pad);
    double val = 0;
    gTPC.GetCDe(layer,row,1,val);
    if(val==0) continue;
    if(tidG4[ih] == id) count++;
  }
  return count;
}
TLorentzVector ToHelix(TLorentzVector GlobalLV){
  double E = GlobalLV.E();
  double X = -GlobalLV.X();
  double Y = GlobalLV.Z();
  double Z = GlobalLV.Y();
  return TLorentzVector(X,Y,Z,E);
}
TLorentzVector ToGlobal(TLorentzVector HelixLV){
  return ToHelix(HelixLV);
}
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kE42, kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[DstTPCTrackingHelixGeant4]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________________
struct Event
{
  Int_t status;
  Int_t evnum;

  Int_t nhHtof;
  std::vector<Double_t> HtofSeg;
  std::vector<Double_t> tHtof;
  std::vector<Double_t> deHtof;
  std::vector<Double_t> posHtof;

  Int_t nhFtof;
  std::vector<Double_t> FtofSeg;
  std::vector<Double_t> tFtof;
  std::vector<Double_t> deFtof;
  std::vector<Double_t> posFtof;



  int G4kmid;
  int G4kmtid;
  double G4kmvtx_x;// Production vertex, identical to xi vert.
  double G4kmvtx_y;
  double G4kmvtx_z;
  double G4kmmom;
  double G4kmmom_x;
  double G4kmmom_y;
  double G4kmmom_z;

  int G4kpid;
  int G4kptid;
  double G4kpvtx_x;// Production vertex, identical to xi vert.
  double G4kpvtx_y;
  double G4kpvtx_z;
  double G4kpmom;
  double G4kpmom_x;
  double G4kpmom_y;
  double G4kpmom_z;

  int G4xiid;
  double G4xivtx_x;//Production vertex
  double G4xivtx_y;
  double G4xivtx_z;
  double G4ximom;
  double G4ximom_x;//Momentum at production vtx
  double G4ximom_y;
  double G4ximom_z;

  double xivtx_x;
  double xivtx_y;
  double xivtx_z;
  double ximom;

  int G4lid;
  double G4lvtx_x;//Production vertex, identical to xi decay vert
  double G4lvtx_y;
  double G4lvtx_z;
  double G4lmom;
  double G4lmom_x;
  double G4lmom_y;
  double G4lmom_z;

  double lvtx_x;
  double lvtx_y;
  double lvtx_z;
  double lmom;

  int G4pid;//Id of Geant4 Track
  int G4ptid;//G4 Id of TPC Track
  int G4pnh;//Number of G4 Hits
  int G4ptnh;//Number of correct G4 Hits
  double G4pvtx_x;//Production vertex, identical to l decay vert
  double G4pvtx_y;
  double G4pvtx_z;
  double G4pmom;
  double G4pmom_x;
  double G4pmom_y;
  double G4pmom_z;

  int ptid;
  int pnh;
  double pvtx_x;
  double pvtx_y;
  double pvtx_z;
  double pmom;
  double pmom_x;
  double pmom_y;
  double pmom_z;
  double pres_x;
  double pres_y;
  double pres_z;
  double pcov_xy;
  double pcov_yz;
  double pcov_zx;
  double GFpmom;
  double GFpmom_x;
  double GFpmom_y;
  double GFpmom_z;


  int G4pi1id;
  int G4pi1tid;
  int G4pi1nh;
  int G4pi1tnh;
  double G4pi1vtx_x;
  double G4pi1vtx_y;
  double G4pi1vtx_z;
  double G4pi1mom;
  double G4pi1mom_x;
  double G4pi1mom_y;
  double G4pi1mom_z;

  int pi1tid;
  int pi1nh;
  double pi1vtx_x;
  double pi1vtx_y;
  double pi1vtx_z;
  double pi1mom;
  double pi1mom_x;
  double pi1mom_y;
  double pi1mom_z;
  double pi1res_x;
  double pi1res_y;
  double pi1res_z;
  double pi1cov_xy;
  double pi1cov_yz;
  double pi1cov_zx;
  double GFpi1mom;
  double GFpi1mom_x;
  double GFpi1mom_y;
  double GFpi1mom_z;


  int G4pi2id;
  int G4pi2tid;
  int G4pi2nh;
  int G4pi2tnh;
  double G4pi2vtx_x;
  double G4pi2vtx_y;
  double G4pi2vtx_z;
  double G4pi2mom;
  double G4pi2mom_x;
  double G4pi2mom_y;
  double G4pi2mom_z;

  int pi2tid;
  int pi2nh;
  double pi2vtx_x;
  double pi2vtx_y;
  double pi2vtx_z;
  double pi2mom;
  double pi2mom_x;
  double pi2mom_y;
  double pi2mom_z;
  double pi2res_x;
  double pi2res_y;
  double pi2res_z;
  double pi2cov_xy;
  double pi2cov_yz;
  double pi2cov_zx;
  double GFpi2mom;
  double GFpi2mom_x;
  double GFpi2mom_y;
  double GFpi2mom_z;







  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Int_t> trackid; //for Kurama K1.8 tracks
  std::vector<Int_t> isBeam;
  std::vector<Int_t> isKurama;
  std::vector<Int_t> isK18;
  std::vector<Int_t> isAccidental;
  std::vector<Int_t> charge; //Helix charge
  std::vector<Int_t> pid;
  std::vector<Double_t> chisqr;
  std::vector<Double_t> pval;
  std::vector<Double_t> purity;
  std::vector<Double_t> efficiency;
  std::vector<Int_t> G4tid;
  std::vector<Double_t> helix_cx;
  std::vector<Double_t> helix_cy;
  std::vector<Double_t> helix_z0;
  std::vector<Double_t> helix_r;
  std::vector<Double_t> helix_dz;
  std::vector<Double_t> dE;
  std::vector<Double_t> dEdx; //reference dedx
  std::vector<Double_t> mom0; //Helix momentum at Y = 0
  std::vector<Double_t> path; //Helix path

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
  std::vector<std::vector<Double_t>> resolution_x;
  std::vector<std::vector<Double_t>> resolution_y;
  std::vector<std::vector<Double_t>> resolution_z;
  std::vector<std::vector<Double_t>> helix_t;
  std::vector<std::vector<Double_t>> alpha;
  std::vector<std::vector<Double_t>> track_cluster_de;
  std::vector<std::vector<Double_t>> track_cluster_mrow;


  int best;
  Bool_t xiflag;
  Double_t ximass;
  Double_t xidecayvtx_x;
  Double_t xidecayvtx_y;
  Double_t xidecayvtx_z;
  Double_t ximom_x;
  Double_t ximom_y;
  Double_t ximom_z;
  Double_t lpi_dist;

  Bool_t lflag;
  Double_t lmass;
  Double_t ldecayvtx_x;
  Double_t ldecayvtx_y;
  Double_t ldecayvtx_z;
  Double_t lmom_x;
  Double_t lmom_y;
  Double_t lmom_z;
  Double_t ppi_dist;
  std::vector<Int_t> decays_id;
  std::vector<Double_t> decays_purity;
  std::vector<Double_t> decays_efficiency;
  std::vector<Int_t> decays_G4tid;
  std::vector<Double_t> decays_mom;
  std::vector<Double_t> decays_mom_x;
  std::vector<Double_t> decays_mom_y;
  std::vector<Double_t> decays_mom_z;
  std::vector<Double_t> decays_res_mom;
  std::vector<Double_t> decays_res_mom_x;
  std::vector<Double_t> decays_res_mom_y;
  std::vector<Double_t> decays_res_mom_z;
  std::vector<Double_t> decays_res_mom_t;
  std::vector<Double_t> decays_res_th;
  std::vector<Double_t> decays_res_ph;
  std::vector<Double_t> decays_cov_mom_th;
  std::vector<Double_t> decays_cov_mom_ph;
  std::vector<Double_t> decays_cov_mom_xy;
  std::vector<Double_t> decays_cov_mom_yz;
  std::vector<Double_t> decays_cov_mom_zx;

  Bool_t GFxiflag;
  Double_t GFximass;
  Double_t GFxidecayvtx_x;
  Double_t GFxidecayvtx_y;
  Double_t GFxidecayvtx_z;
  Double_t GFximom;
  Double_t GFximom_x;
  Double_t GFximom_y;
  Double_t GFximom_z;
  Double_t GFlpi_dist;
  Double_t KFlpi_dist;

  Bool_t GFlflag;
  Double_t GFlmass;
  Double_t GFldecayvtx_x;
  Double_t GFldecayvtx_y;
  Double_t GFldecayvtx_z;
  Double_t GFlmom;
  Double_t GFlmom_x;
  Double_t GFlmom_y;
  Double_t GFlmom_z;
  Double_t GFppi_dist;
  std::vector<Double_t> GFdecays_mass;
  std::vector<Double_t> GFdecays_mom;
  std::vector<Double_t> GFdecays_mom_x;
  std::vector<Double_t> GFdecays_mom_y;
  std::vector<Double_t> GFdecays_mom_z;
  std::vector<Double_t> GFmomloss;
  std::vector<Double_t> GFeloss;

  Int_t GFstatus;
  Int_t GFntTpc;
  std::vector<Int_t> GFfitstatus;
  std::vector<Int_t> GFpdgcode;
  std::vector<Int_t> GFnhtrack;
  std::vector<Double_t> GFcharge;
  std::vector<Double_t> GFchisqr;
  std::vector<Double_t> GFtof;
  std::vector<Double_t> GFtracklen;
  std::vector<Double_t> GFpval;
  std::vector<Double_t> GFchisqrPos;
  std::vector<Double_t> GFpvalPos;
  std::vector<std::vector<Double_t>> GFlayer;
  std::vector<std::vector<Double_t>> GFrow;
  std::vector<std::vector<Double_t>> GFpos_x;
  std::vector<std::vector<Double_t>> GFpos_y;
  std::vector<std::vector<Double_t>> GFpos_z;
  std::vector<std::vector<Double_t>> GFmom;
  std::vector<std::vector<Double_t>> GFmom_x;
  std::vector<std::vector<Double_t>> GFmom_y;
  std::vector<std::vector<Double_t>> GFmom_z;
  std::vector<std::vector<Double_t>> GFresidual_x;
  std::vector<std::vector<Double_t>> GFresidual_y;
  std::vector<std::vector<Double_t>> GFresidual_z;
  std::vector<std::vector<Double_t>> GFresidual_p;
  std::vector<std::vector<Double_t>> GFresidual_px;
  std::vector<std::vector<Double_t>> GFresidual_py;
  std::vector<std::vector<Double_t>> GFresidual_pz;
  std::vector<std::vector<Double_t>> GFresolution_x;
  std::vector<std::vector<Double_t>> GFresolution_y;
  std::vector<std::vector<Double_t>> GFresolution_z;
  std::vector<std::vector<Double_t>> GFresolution_p;
  std::vector<std::vector<Double_t>> GFresolution_px;
  std::vector<std::vector<Double_t>> GFresolution_py;
  std::vector<std::vector<Double_t>> GFresolution_pz;
  std::vector<std::vector<Double_t>> GFpull_x;
  std::vector<std::vector<Double_t>> GFpull_y;
  std::vector<std::vector<Double_t>> GFpull_z;
  std::vector<std::vector<Double_t>> GFpull_p;
  std::vector<std::vector<Double_t>> GFpull_px;
  std::vector<std::vector<Double_t>> GFpull_py;
  std::vector<std::vector<Double_t>> GFpull_pz;

  Double_t KFchisqrxi;
  Double_t KFpvalxi;
  vector<double> KFxidecays_pull;
  Double_t KFximass;//MassBefore;
  Double_t KFximom;
  Double_t KFximom_x;
  Double_t KFximom_y;
  Double_t KFximom_z;

  Double_t KFchisqrl;
  Double_t KFpvall;
  vector<double> KFldecays_pull;
  Double_t KFlmass;//MassBefore;
  Double_t KFlmom;
  Double_t KFlmom_x;
  Double_t KFlmom_y;
  Double_t KFlmom_z;

  vector<Double_t> KFdecays_mom;
  vector<Double_t> KFdecays_mom_x;
  vector<Double_t> KFdecays_mom_y;
  vector<Double_t> KFdecays_mom_z;

  bool XiFlight;
  bool XiProd;

  double xtgtXi;
  double ytgtXi;
  double utgtXi;
  double vtgtXi;

  double vtxKKXi;
  double vtyKKXi;
  double vtzKKXi;

  double xiprodvtx_x;
  double xiprodvtx_y;
  double xiprodvtx_z;

  double xiprodmom_x;
  double xiprodmom_y;
  double xiprodmom_z;

  vector<double>xtgtHS;
  vector<double>ytgtHS;
  vector<double>xtgtKurama;
  vector<double>ytgtKurama;



  void clear( void )
  {
    evnum = 0;
    status = 0;

    nhHtof = 0;
    HtofSeg.clear();
    tHtof.clear();
    deHtof.clear();
    posHtof.clear();

    nhFtof = 0;
    FtofSeg.clear();
    tFtof.clear();
    deFtof.clear();
    posFtof.clear();


    G4kmid = -1;
    G4kmtid = -1;
    G4kmvtx_x = qnan;
    G4kmvtx_y = qnan;
    G4kmvtx_z = qnan;
    G4kmmom = qnan;
    G4kmmom_x = qnan;
    G4kmmom_y = qnan;
    G4kmmom_z = qnan;

    G4kpid = -1;
    G4kptid = -1;
    G4kpvtx_x = qnan;
    G4kpvtx_y = qnan;
    G4kpvtx_z = qnan;
    G4kpmom = qnan;
    G4kpmom_x = qnan;
    G4kpmom_y = qnan;
    G4kpmom_z = qnan;

    G4xiid = -1;
    G4xivtx_x = qnan;
    G4xivtx_y = qnan;
    G4xivtx_z = qnan;
    G4ximom = qnan;
    G4ximom_x = qnan;
    G4ximom_y = qnan;
    G4ximom_z = qnan;
    xivtx_x = qnan;
    xivtx_y = qnan;
    xivtx_z = qnan;
    ximom = qnan;

    G4lid = -1;
    G4lvtx_x = qnan;
    G4lvtx_y = qnan;
    G4lvtx_z = qnan;
    G4lmom = qnan;
    G4lmom_x = qnan;
    G4lmom_y = qnan;
    G4lmom_z = qnan;
    lvtx_x = qnan;
    lvtx_y = qnan;
    lvtx_z = qnan;
    lmom = qnan;

    G4pid = -1;
    G4ptid = -1;
    G4pnh = -1;
    G4ptnh = -1;
    G4pvtx_x = qnan;
    G4pvtx_y = qnan;
    G4pvtx_z = qnan;
    G4pmom = qnan;
    G4pmom_x = qnan;
    G4pmom_y = qnan;
    G4pmom_z = qnan;
    ptid = -1;
    pnh = -1;
    pvtx_x = qnan;
    pvtx_y = qnan;
    pvtx_z = qnan;
    pmom = qnan;
    pmom_x = qnan;
    pmom_y = qnan;
    pmom_z = qnan;
    pres_x = qnan;
    pres_y = qnan;
    pres_z = qnan;
    pcov_xy = qnan;
    pcov_yz = qnan;
    pcov_zx = qnan;
    GFpmom = qnan;
    GFpmom_x = qnan;
    GFpmom_y = qnan;
    GFpmom_z = qnan;

    G4pi1id = -1;
    G4pi1tid = -1;
    G4pi1nh = -1;
    G4pi1tnh = -1;
    G4pi1vtx_x = qnan;
    G4pi1vtx_y = qnan;
    G4pi1vtx_z = qnan;
    G4pi1mom = qnan;
    G4pi1mom_x = qnan;
    G4pi1mom_y = qnan;
    G4pi1mom_z = qnan;
    pi1tid = -1;
    pi1nh = -1;
    pi1vtx_x = qnan;
    pi1vtx_y = qnan;
    pi1vtx_z = qnan;
    pi1mom = qnan;
    pi1mom_x = qnan;
    pi1mom_y = qnan;
    pi1mom_z = qnan;
    pi1res_x = qnan;
    pi1res_y = qnan;
    pi1res_z = qnan;
    pi1cov_xy = qnan;
    pi1cov_yz = qnan;
    pi1cov_zx = qnan;
    GFpi1mom = qnan;
    GFpi1mom_x = qnan;
    GFpi1mom_y = qnan;
    GFpi1mom_z = qnan;

    G4pi2id = -1;
    G4pi2tid = -1;
    G4pi2nh = -1;
    G4pi2tnh = -1;
    G4pi2vtx_x = qnan;
    G4pi2vtx_y = qnan;
    G4pi2vtx_z = qnan;
    G4pi2mom = qnan;
    G4pi2mom_x = qnan;
    G4pi2mom_y = qnan;
    G4pi2mom_z = qnan;
    pi2tid = -1;
    pi2nh = -1;
    pi2vtx_x = qnan;
    pi2vtx_y = qnan;
    pi2vtx_z = qnan;
    pi2mom = qnan;
    pi2mom_x = qnan;
    pi2mom_y = qnan;
    pi2mom_z = qnan;
    pi2res_x = qnan;
    pi2res_y = qnan;
    pi2res_z = qnan;
    pi2cov_xy = qnan;
    pi2cov_yz = qnan;
    pi2cov_zx = qnan;
    GFpi2mom = qnan;
    GFpi2mom_x = qnan;
    GFpi2mom_y = qnan;
    GFpi2mom_z = qnan;

    ntTpc = 0;
    nhtrack.clear();
    trackid.clear();
    isBeam.clear();
    isKurama.clear();
    isK18.clear();
    isAccidental.clear();
    charge.clear();
    pid.clear();

    purity.clear();
    efficiency.clear();
    G4tid.clear();
    chisqr.clear();
    pval.clear();
    helix_cx.clear();
    helix_cy.clear();
    helix_z0.clear();
    helix_r.clear();
    helix_dz.clear();
    dE.clear();
    dEdx.clear();
    mom0.clear();
    path.clear();

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
    resolution_x.clear();
    resolution_y.clear();
    resolution_z.clear();
    helix_t.clear();
    alpha.clear();
    track_cluster_de.clear();
    track_cluster_mrow.clear();

    best = -1;
    xiflag = false;
    ximass = qnan;
    xidecayvtx_x = qnan;
    xidecayvtx_y = qnan;
    xidecayvtx_z = qnan;
    ximom_x = qnan;
    ximom_y = qnan;
    ximom_z = qnan;
    lpi_dist = qnan;

    lflag = false;
    lmass = qnan;
    ldecayvtx_x = qnan;
    ldecayvtx_y = qnan;
    ldecayvtx_z = qnan;
    lmom_x = qnan;
    lmom_y = qnan;
    lmom_z = qnan;
    ppi_dist = qnan;
    decays_id.clear();
    decays_purity.clear();
    decays_efficiency.clear();
    decays_G4tid.clear();
    decays_mom.clear();
    decays_mom_x.clear();
    decays_mom_y.clear();
    decays_mom_z.clear();
    decays_res_mom.clear();
    decays_res_mom_x.clear();
    decays_res_mom_y.clear();
    decays_res_mom_z.clear();
    decays_res_mom_t.clear();
    decays_res_th.clear();
    decays_res_ph.clear();
    decays_cov_mom_th.clear();
    decays_cov_mom_ph.clear();
    decays_cov_mom_xy.clear();
    decays_cov_mom_yz.clear();
    decays_cov_mom_zx.clear();

    GFxiflag = false;
    GFximass = qnan;
    GFxidecayvtx_x = qnan;
    GFxidecayvtx_y = qnan;
    GFxidecayvtx_z = qnan;
    GFximom = qnan;
    GFximom_x = qnan;
    GFximom_y = qnan;
    GFximom_z = qnan;
    GFlpi_dist = qnan;
    KFlpi_dist = qnan;

    GFlflag = false;
    GFlmass = qnan;
    GFldecayvtx_x = qnan;
    GFldecayvtx_y = qnan;
    GFldecayvtx_z = qnan;
    GFlmom = qnan;
    GFlmom_x = qnan;
    GFlmom_y = qnan;
    GFlmom_z = qnan;
    GFppi_dist = qnan;

    GFdecays_mass.clear();
    GFdecays_mom.clear();
    GFdecays_mom_x.clear();
    GFdecays_mom_y.clear();
    GFdecays_mom_z.clear();
    GFmomloss.clear();
    GFeloss.clear();

    GFstatus = 0;
    GFntTpc = 0;
    GFcharge.clear();
    GFchisqr.clear();
    GFpval.clear();
    GFchisqrPos.clear();
    GFpvalPos.clear();
    GFtof.clear();
    GFtracklen.clear();
    GFpval.clear();
    GFfitstatus.clear();
    GFpdgcode.clear();
    GFnhtrack.clear();
    GFlayer.clear();
    GFrow.clear();
    GFpos_x.clear();
    GFpos_y.clear();
    GFpos_z.clear();
    GFmom.clear();
    GFmom_x.clear();
    GFmom_y.clear();
    GFmom_z.clear();
    GFresidual_x.clear();
    GFresidual_y.clear();
    GFresidual_z.clear();
    GFresidual_p.clear();
    GFresidual_px.clear();
    GFresidual_py.clear();
    GFresidual_pz.clear();
    GFresolution_x.clear();
    GFresolution_y.clear();
    GFresolution_z.clear();
    GFresolution_p.clear();
    GFresolution_px.clear();
    GFresolution_py.clear();
    GFresolution_pz.clear();
    GFpull_x.clear();
    GFpull_y.clear();
    GFpull_z.clear();
    GFpull_p.clear();
    GFpull_px.clear();
    GFpull_py.clear();
    GFpull_pz.clear();

    KFchisqrxi = qnan;
    KFpvalxi = qnan;
    KFxidecays_pull.clear();
    KFximass = qnan;
    KFximom = qnan;
    KFximom_x = qnan;
    KFximom_y = qnan;
    KFximom_z = qnan;

    KFchisqrl = qnan;
    KFpvall = qnan;
    KFldecays_pull.clear();
    KFlmass = qnan;
    KFlmom = qnan;
    KFlmom_x = qnan;
    KFlmom_y = qnan;
    KFlmom_z = qnan;

    KFdecays_mom.clear();
    KFdecays_mom_x.clear();
    KFdecays_mom_y.clear();
    KFdecays_mom_z.clear();

    XiFlight = false;
    XiProd = false;
    xtgtXi = qnan;
    ytgtXi = qnan;
    utgtXi = qnan;
    vtgtXi = qnan;
    vtxKKXi = qnan;
    vtyKKXi = qnan;
    vtzKKXi = qnan;
    xiprodvtx_x = qnan;
    xiprodvtx_y = qnan;
    xiprodvtx_z = qnan;
    xiprodmom_x = qnan;
    xiprodmom_y = qnan;
    xiprodmom_z = qnan;

    xtgtHS.clear();
    ytgtHS.clear();
    xtgtKurama.clear();
    ytgtKurama.clear();

  }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* evnum;

  /*
  TTreeReaderValue<Int_t>* nhHtof;
  TTreeReaderValue<std::vector<Double_t>>* HtofSeg;
  TTreeReaderValue<std::vector<Double_t>>* tHtof;
  TTreeReaderValue<std::vector<Double_t>>* dtHtof;
  TTreeReaderValue<std::vector<Double_t>>* deHtof;
  TTreeReaderValue<std::vector<Double_t>>* posHtof;
*/
  int nhHtof;
  int didHtof[500];
  double tHtof[500];
  double deHtof[500];
  double yHtof[500];//->posHtof

  int nhFtof;
  int didFtof[500];
  double tFtof[500];
  double deFtof[500];
  double yFtof[500];//->posHtof

  int nhittpc;
  Int_t ititpc[MaxTPCHits];
  Double_t xtpc[MaxTPCHits];//with resolution
  Double_t ytpc[MaxTPCHits];//with resolution
  Double_t ztpc[MaxTPCHits];//with resolution
  Double_t pxtpc[MaxTPCHits];//with resolution
  Double_t pytpc[MaxTPCHits];//with resolution
  Double_t pztpc[MaxTPCHits];//with resolution

  Int_t NumberOfTracks;
  Int_t PIDOfTrack[1000];
  Int_t ParentIDOfTrack[1000];
  Double_t VertexOfTrack_x[1000];
  Double_t VertexOfTrack_y[1000];
  Double_t VertexOfTrack_z[1000];
  Double_t MomentumOfTrack[1000];
  Double_t MomentumOfTrack_x[1000];
  Double_t MomentumOfTrack_y[1000];
  Double_t MomentumOfTrack_z[1000];



  TTreeReaderValue<Int_t>* ntTpc; // Number of Tracks
  TTreeReaderValue<std::vector<Int_t>>* nhtrack; // Number of Hits (in 1 tracks)
  TTreeReaderValue<std::vector<Int_t>>* trackid; //for Kurama K1.8 tracks
  TTreeReaderValue<std::vector<Int_t>>* isBeam;
  TTreeReaderValue<std::vector<Int_t>>* isKurama;
  TTreeReaderValue<std::vector<Int_t>>* isK18;
  TTreeReaderValue<std::vector<Int_t>>* isAccidental;

  TTreeReaderValue<std::vector<Int_t>>* charge;//Helix charge
  TTreeReaderValue<std::vector<Int_t>>* pid;
  TTreeReaderValue<std::vector<Double_t>>* purity;
  TTreeReaderValue<std::vector<Double_t>>* efficiency;
  TTreeReaderValue<std::vector<Int_t>>* G4tid;
  TTreeReaderValue<std::vector<Double_t>>* chisqr;
  TTreeReaderValue<std::vector<Double_t>>* pval;
  TTreeReaderValue<std::vector<Double_t>>* helix_cx;
  TTreeReaderValue<std::vector<Double_t>>* helix_cy;
  TTreeReaderValue<std::vector<Double_t>>* helix_z0;
  TTreeReaderValue<std::vector<Double_t>>* helix_r;
  TTreeReaderValue<std::vector<Double_t>>* helix_dz;
  TTreeReaderValue<std::vector<Double_t>>* dE;
  TTreeReaderValue<std::vector<Double_t>>* dEdx; //reference dedx
  TTreeReaderValue<std::vector<Double_t>>* mom0;//Helix momentum at Y = 0
  TTreeReaderValue<std::vector<Double_t>>* path;//Helix path
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitlayer;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitpos_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitpos_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitpos_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* calpos_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* calpos_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* calpos_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* helix_t;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* alpha;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_mrow;
	TTreeReaderValue<vector<double>>* xtgtHS=nullptr;
	TTreeReaderValue<vector<double>>* ytgtHS=nullptr;
	TTreeReaderValue<vector<double>>* xtgtKurama=nullptr;
	TTreeReaderValue<vector<double>>* ytgtKurama=nullptr;



};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
  enum eDetHid {
    TPCHid    = 100000,
    genfitHid = 200000
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
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries( TTreeCont );
  if (max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();

  //Initiallize Geometry, Field, Fitter
  HypTPCFitter* fitter = new HypTPCFitter(tpcGeo.Data(), Const_field);
  //Initiallize the genfit track container
  HypTPCTask& GFTrackCont = HypTPCTask::GetInstance();
  GFTrackCont.SetVerbosity(verbosity);
  std::cout<<"GenFit verbosity = "<<"-1: Silent, 0: Minimum, 1: Errors only, 2: Errors and Warnings, 3: Verbose mode, long term debugging(default)"<<std::endl;
  std::cout<<"Current verbosity = "<<GFTrackCont.GetVerbosity()<<std::endl;

#if 0
  GFTrackCont.DebugMode();
#endif

  Int_t ievent = skip;
  for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
    gCounter.check();
    InitializeEvent();
    if( DstRead( ievent ) ) tree->Fill();
  }

  std::cout << "#D Event Number: " << std::setw(6)
            << ievent << std::endl;

  DstClose();

  delete fitter;
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
  Int_t open_file = 0;
  Int_t open_tree = 0;
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
dst::DstRead( Int_t ievent )
{
  static const auto PionMass    = pdg::PionMass();
  static const auto ProtonMass  = pdg::ProtonMass();
  static const auto LambdaMass  = pdg::LambdaMass();
  static const auto XiMinusMass = pdg::XiMinusMass();
  static const int XiMinusPdgCode = 3312;
  int cout_scale = 1000;
  //if( ievent%1000==0 ){
  if( ievent%cout_scale==0 ){
  //if( ievent%100000==0 ){
    std::cout << "#D Event Number: "
        << std::setw(6) << ievent << std::endl;
  }
  TVector3 tgtpos(0, 0, tpc::ZTarget);

  GetEntry(ievent);
  event.evnum = **src.evnum;

  HF1( 1, event.status++ );

  event.nhHtof = src.nhHtof;
  for(int ih=0;ih<event.nhHtof;++ih){
    event.HtofSeg.push_back(src.didHtof[ih]);
    event.tHtof.push_back(src.tHtof[ih]);
    event.deHtof.push_back(src.deHtof[ih]);
    event.posHtof.push_back(src.yHtof[ih]);
  }

  event.nhFtof = src.nhFtof;
  for(int ih=0;ih<event.nhFtof;++ih){
    event.FtofSeg.push_back(src.didFtof[ih]);
    event.tFtof.push_back(src.tFtof[ih]);
    event.deFtof.push_back(src.deFtof[ih]);
    event.posFtof.push_back(src.yFtof[ih]);
  }

  Int_t ntTpc = **src.ntTpc;
  if( ntTpc == 0 ) return true;
    std::cout << "#D Event Number: "
        << std::setw(6) << ievent << std::endl;
  HF1( 1, event.status++ );

  event.ntTpc = ntTpc;
  event.nhtrack = **src.nhtrack;
  event.trackid = **src.trackid;
  event.isBeam = **src.isBeam;
  event.isKurama = **src.isKurama;
  event.isK18 = **src.isK18;
  event.isAccidental = **src.isAccidental;
  event.charge = **src.charge;
  event.pid = **src.pid;
  event.purity = **src.purity;
  event.G4tid = **src.G4tid;
  event.efficiency = **src.efficiency;
  event.chisqr = **src.chisqr;
  event.pval = **src.pval;
  event.helix_cx = **src.helix_cx;
  event.helix_cy = **src.helix_cy;
  event.helix_z0 = **src.helix_z0;
  event.helix_r = **src.helix_r;
  event.helix_dz = **src.helix_dz;
  event.dE = **src.dE;
  event.dEdx = **src.dEdx;
  event.mom0 = **src.mom0;
  event.path = **src.path;
  event.hitlayer = **src.hitlayer;
  event.hitpos_x = **src.hitpos_x;
  event.hitpos_y = **src.hitpos_y;
  event.hitpos_z = **src.hitpos_z;
  event.calpos_x = **src.calpos_x;
  event.calpos_y = **src.calpos_y;
  event.calpos_z = **src.calpos_z;
  event.residual = **src.residual;
  event.residual_x = **src.residual_x;
  event.residual_y = **src.residual_y;
  event.residual_z = **src.residual_z;
  event.resolution_x = **src.resolution_x;
  event.resolution_y = **src.resolution_y;
  event.resolution_z = **src.resolution_z;
  event.helix_t = **src.helix_t;
  event.alpha = **src.alpha;
  event.track_cluster_de = **src.track_cluster_de;
  event.track_cluster_mrow = **src.track_cluster_mrow;
	event.xtgtHS = **src.xtgtHS;
	event.ytgtHS = **src.ytgtHS;
	event.xtgtKurama = **src.xtgtKurama;
	event.ytgtKurama = **src.ytgtKurama;


  vector<TVector3> G4Hits;
  vector<TVector3> G4Moms;

  for(int ih=0;ih<src.nhittpc;++ih){
    TVector3 G4Hit(src.xtpc[ih],src.ytpc[ih],src.ztpc[ih]);
    TVector3 G4Mom(src.pxtpc[ih],src.pytpc[ih],src.pztpc[ih]);
    G4Hits.push_back(G4Hit);
    G4Moms.push_back(G4Mom);
  }
  Int_t l_candidates = 0;
  //Xi candidates searching
  Int_t xi_candidates = 0;

  Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  std::vector<Int_t> xi_p_container, xi_pi_container, xi_pi2_container;
  std::vector<TVector3> xi_mom_container, xi_vert_container;
  std::vector<TVector3> l_mom_container, l_vert_container;
  std::vector<Double_t> xi_mass_container, lambda_mass_container;
  std::vector<TVector3> p_mom_container, pi_mom_container, pi2_mom_container;
  std::vector<TVector3> p_res_container, pi_res_container, pi2_res_container;
  std::vector<TVector3> p_cov_container, pi_cov_container, pi2_cov_container;
  std::vector<Double_t> ppi_closedist; std::vector<Double_t> lpi_closedist;
  std::vector<Double_t> ppi_closedist_l;
  std::vector<Double_t> l_mass_container_l;
  std::vector<Int_t> l_p_container, l_pi_container;
  std::vector<TVector3> l_mom_container_l, l_vert_container_l;
  std::vector<TVector3> p_mom_container_l, pi_mom_container_l;


  TPCAnalyzer TPCAna;
#if WithKurama
  TPCAna.ReCalcTPCTracks(**src.ntTpc, **src.isK18, **src.isKurama,
        **src.charge, **src.nhtrack, **src.helix_cx,
        **src.helix_cy, **src.helix_z0, **src.helix_r,
        **src.helix_dz, **src.hitlayer, **src.track_cluster_mrow,
        **src.helix_t, **src.track_cluster_de, **src.resolution_x,
        **src.resolution_y, **src.resolution_z, **src.hitpos_x,
        **src.hitpos_y, **src.hitpos_z);
#else
    TPCAna.ReCalcTPCTracksGeant4(**src.ntTpc,
       **src.charge, **src.nhtrack, **src.helix_cx,
       **src.helix_cy, **src.helix_z0, **src.helix_r,
       **src.helix_dz, **src.hitlayer, **src.track_cluster_mrow,
       **src.helix_t, **src.track_cluster_de, **src.resolution_x,
       **src.resolution_y, **src.resolution_z, **src.hitpos_x,
       **src.hitpos_y, **src.hitpos_z);
#endif
  for(int it=1;it<=src.NumberOfTracks;++it){
    int parent = src.ParentIDOfTrack[it];
    if(parent<0) continue;
    int pid = src.PIDOfTrack[it];
    double mom = src.MomentumOfTrack[it]/1000;//MeV to GeV
    double mom_x = src.MomentumOfTrack_x[it]/1000;
    double mom_y = src.MomentumOfTrack_y[it]/1000;
    double mom_z = src.MomentumOfTrack_z[it]/1000;
    double vert_x = src.VertexOfTrack_x[it];
    double vert_y = src.VertexOfTrack_y[it];
    double vert_z = src.VertexOfTrack_z[it];
    if(abs(pid) == 321 and parent ==0){
      if(mom_z < 0){
        event.G4kmid = it;
        event.G4kmvtx_x = vert_x;
        event.G4kmvtx_y = vert_y;
        event.G4kmvtx_z = vert_z;
        event.G4kmmom = mom;
        event.G4kmmom_x = -mom_x;
        event.G4kmmom_y = -mom_y;
        event.G4kmmom_z = -mom_z;
      }
      else {
        event.G4kpid = it;
        event.G4kpvtx_x = vert_x;
        event.G4kpvtx_y = vert_y;
        event.G4kpvtx_z = vert_z;
        event.G4kpmom = mom;
        event.G4kpmom_x = mom_x;
        event.G4kpmom_y = mom_y;
        event.G4kpmom_z = mom_z;
      }
    }
    if(pid==3312){//Xi
        event.G4xiid = it;
        event.G4xivtx_x = vert_x;
        event.G4xivtx_y = vert_y;
        event.G4xivtx_z = vert_z;
        event.G4ximom = mom;
        event.G4ximom_x = mom_x;
        event.G4ximom_y = mom_y;
        event.G4ximom_z = mom_z;
    }
    if(pid==3122 and src.PIDOfTrack[parent] == 3312){//Lambda
        event.G4lid = it;
        event.G4lvtx_x = vert_x;
        event.G4lvtx_y = vert_y;
        event.G4lvtx_z = vert_z;
        event.G4lmom = mom;
        event.G4lmom_x = mom_x;
        event.G4lmom_y = mom_y;
        event.G4lmom_z = mom_z;
    }
    if(pid == 2212 and src.PIDOfTrack[parent] == 3122){//Proton
        event.G4pid = it;
        event.G4pvtx_x = vert_x;
        event.G4pvtx_y = vert_y;
        event.G4pvtx_z = vert_z;
        event.G4pmom = mom;
        event.G4pmom_x = mom_x;
        event.G4pmom_y = mom_y;
        event.G4pmom_z = mom_z;
    }
    if(abs(pid) == 211){
      if(src.PIDOfTrack[parent] == 3312){//Pi form Xi
        event.G4pi2id = it;
        event.G4pi2vtx_x = vert_x;
        event.G4pi2vtx_y = vert_y;
        event.G4pi2vtx_z = vert_z;
        event.G4pi2mom = mom;
        event.G4pi2mom_x = mom_x;
        event.G4pi2mom_y = mom_y;
        event.G4pi2mom_z = mom_z;
      }
      else if(src.PIDOfTrack[parent] == 3122){//Pi from L
        event.G4pi1id = it;
        event.G4pi1vtx_x = vert_x;
        event.G4pi1vtx_y = vert_y;
        event.G4pi1vtx_z = vert_z;
        event.G4pi1mom = mom;
        event.G4pi1mom_x = mom_x;
        event.G4pi1mom_y = mom_y;
        event.G4pi1mom_z = mom_z;
      }
    }
  }
  vector<int> G4TrackID;
  vector<int> PureHits;

  event.G4pnh = CountHits(event.G4pid,src.xtpc,src.ztpc,src.ititpc,src.nhittpc);
  event.G4pi1nh = CountHits(event.G4pi1id,src.xtpc,src.ztpc,src.ititpc,src.nhittpc);
  event.G4pi2nh = CountHits(event.G4pi2id,src.xtpc,src.ztpc,src.ititpc,src.nhittpc);
  for(int ih=0;ih<src.nhittpc;++ih){
    double x = src.xtpc[ih],z=src.ztpc[ih];
    int pad = tpc::findPadID(z,x);
    int layer = tpc::getLayerID(pad);
    int row = tpc::getRowID(pad);
    double val = 0;
    gTPC.GetCDe(layer,row,1,val);
    if(val==0) continue;
    if(src.ititpc[ih]==event.G4pid){
      HF2(2001,z,x);
    }
    if(src.ititpc[ih]==event.G4pi1id){
      HF2(2002,z,x);
    }
    if(src.ititpc[ih]==event.G4pi2id){
      HF2(2003,z,x);
    }
  }



  for(Int_t it1=0;it1<ntTpc;it1++){
    TPCLocalTrackHelix *tp1 = TPCAna.GetTrackTPCHelix(it1);

    vector<TVector3> TPCHit;
    for(int ih=0;ih<event.hitpos_x.at(it1).size();++ih){
      TPCHit.push_back(TVector3(
            event.hitpos_x.at(it1).at(ih),
            event.hitpos_y.at(it1).at(ih),
            event.hitpos_z.at(it1).at(ih)));
    }
    int nPureHit;
    int G4tid = TPCToG4TrackID(TPCHit,src.nhittpc,src.ititpc,src.xtpc,src.ytpc,src.ztpc,nPureHit);
    G4TrackID.push_back(G4tid);
    PureHits.push_back(nPureHit);
    if(event.isK18[it1]==1 || event.isKurama[it1]==1) continue;
    if((event.pid[it1]&4)!=4) continue;
    Bool_t p_like = false;
    if(event.charge[it1]==1) p_like = true;

    Double_t p_par[5];
    p_par[0] = event.helix_cx[it1];
    p_par[1] = event.helix_cy[it1];
    p_par[2] = event.helix_z0[it1];
    p_par[3] = event.helix_r[it1];
    p_par[4] = event.helix_dz[it1];
    Int_t p_nh = event.helix_t[it1].size();
    Double_t p_theta_min = event.helix_t[it1][0] - vtx_scan_range/p_par[3];
    Double_t p_theta_max = TMath::Min(event.helix_t[it1][0] + vtx_scan_rangeInsideL/p_par[3], event.helix_t[it1][p_nh-1]);

    TVector3 p_start = TVector3(event.calpos_x[it1][0], event.calpos_y[it1][0], event.calpos_z[it1][0]);
    TVector3 p_end = TVector3(event.calpos_x[it1][p_nh-1], event.calpos_y[it1][p_nh-1], event.calpos_z[it1][p_nh-1]);
    if(!p_like){
      p_theta_min = event.helix_t[it1][p_nh-1] - vtx_scan_range/p_par[3];
      p_theta_max = TMath::Min(event.helix_t[it1][p_nh-1] + vtx_scan_rangeInsideL/p_par[3], event.helix_t[it1][0]);
      p_start = TVector3(event.calpos_x[it1][p_nh-1], event.calpos_y[it1][p_nh-1], event.calpos_z[it1][p_nh-1]);
      p_end = TVector3(event.calpos_x[it1][0], event.calpos_y[it1][0], event.calpos_z[it1][0]);
    }

    for(Int_t it2=0;it2<ntTpc;it2++){
      if(it1==it2) continue;
      if(event.isK18[it1]==1 || event.isKurama[it1]==1) continue;
      TPCLocalTrackHelix *tp2 = TPCAna.GetTrackTPCHelix(it2);
      if((event.pid[it2]&1)!=1) continue; //select pi like
      Bool_t pim_like = false;
      if(event.charge[it2]==-1) pim_like = true;

      Double_t pi_par[5];
      pi_par[0] = event.helix_cx[it2];
      pi_par[1] = event.helix_cy[it2];
      pi_par[2] = event.helix_z0[it2];
      pi_par[3] = event.helix_r[it2];
      pi_par[4] = event.helix_dz[it2];

      Int_t pi_nh = event.helix_t[it2].size();
      Double_t pi_theta_min = TMath::Max(event.helix_t[it2][0] - vtx_scan_rangeInsideL/pi_par[3], event.helix_t[it2][pi_nh-1]);
      Double_t pi_theta_max = event.helix_t[it2][0] + vtx_scan_range/pi_par[3];

      TVector3 pi_start = TVector3(event.calpos_x[it2][0], event.calpos_y[it2][0], event.calpos_z[it2][0]);
      TVector3 pi_end = TVector3(event.calpos_x[it2][pi_nh-1], event.calpos_y[it2][pi_nh-1], event.calpos_z[it2][pi_nh-1]);
      if(!pim_like){
  pi_theta_min = TMath::Max(event.helix_t[it2][pi_nh-1] - vtx_scan_rangeInsideL/pi_par[3], event.helix_t[it2][0]);
  pi_theta_max = event.helix_t[it2][pi_nh-1] + vtx_scan_range/pi_par[3];
  pi_start = TVector3(event.calpos_x[it2][pi_nh-1], event.calpos_y[it2][pi_nh-1], event.calpos_z[it2][pi_nh-1]);
  pi_end = TVector3(event.calpos_x[it2][0], event.calpos_y[it2][0], event.calpos_z[it2][0]);
      }
      Double_t ppi_dist = 10000.;
      TVector3 p_mom; TVector3 pi_mom; TVector3 lambda_mom;
      TVector3 lambda_vert = Kinematics::LambdaVertex(dMagneticField, p_par, pi_par, p_theta_min, p_theta_max, pi_theta_min, pi_theta_max, p_mom, pi_mom, lambda_mom, ppi_dist);
      if(TMath::IsNaN(ppi_dist)) continue;
      if(!pim_like) pi_mom = -1.*pi_mom;
      if(!p_like) pi_mom = -1.*p_mom;
      lambda_mom = pi_mom + p_mom;

      TLorentzVector Lpi(pi_mom, TMath::Sqrt(pi_mom.Mag()*pi_mom.Mag() + PionMass*PionMass));
      TLorentzVector Lp(p_mom, TMath::Sqrt(p_mom.Mag()*p_mom.Mag() + ProtonMass*ProtonMass));
      TLorentzVector Llambda = Lp + Lpi;
      if(TMath::Abs(lambda_vert.x()) > 250. || TMath::Abs(lambda_vert.z()) > 250. || TMath::Abs(lambda_vert.y()) > 250.) continue; //Vertex cut
      if(ppi_dist > ppi_distcut) continue;
      Double_t pi_vertex_dist; Double_t p_vertex_dist;
      if(!Kinematics::HelixDirection(lambda_vert, p_start, p_end, p_vertex_dist) ||
   !Kinematics::HelixDirection(lambda_vert, pi_start, pi_end, pi_vertex_dist)) continue;

      if(pi_vertex_dist > pi_vtx_distcut) continue;
      if(p_vertex_dist > p_vtx_distcut) continue;

      if(TMath::Abs(Llambda.M() - LambdaMass) > lambda_masscut) continue;
      event.lflag = true;
      l_candidates++;
      ppi_closedist_l.push_back(ppi_dist);
      l_mom_container_l.push_back(lambda_mom);
      l_vert_container_l.push_back(lambda_vert);
      l_p_container.push_back(it1);
      l_pi_container.push_back(it2);
      l_mass_container_l.push_back(Llambda.M());
      p_mom_container_l.push_back(p_mom);
      pi_mom_container_l.push_back(pi_mom);
      for(int it3=0;it3<ntTpc;it3++){
   TPCLocalTrackHelix *tp3 = TPCAna.GetTrackTPCHelix(it3);

  if(it3==it2 || it3==it1) continue;
	if(event.isK18[it3]==1 || event.isKurama[it3]==1) continue;
	if((event.pid[it3]&1)!=1) continue; //select pi like
  Bool_t pim_like2 = false;
  if(event.charge[it3]==-1) pim_like2 = true;

  Double_t pi2_par[5];
  pi2_par[0] = event.helix_cx[it3];
  pi2_par[1] = event.helix_cy[it3];
  pi2_par[2] = event.helix_z0[it3];
  pi2_par[3] = event.helix_r[it3];
  pi2_par[4] = event.helix_dz[it3];
  ////Tag

  Int_t pi2_nh = event.helix_t[it3].size();
  Double_t pi2_theta_min = TMath::Max(event.helix_t[it3][0] - vtx_scan_rangeInsidePi/pi2_par[3], event.helix_t[it3][pi2_nh-1]);
  Double_t pi2_theta_max = event.helix_t[it3][0] + vtx_scan_range/pi2_par[3];

  TVector3 pi2_start = TVector3(event.calpos_x[it3][0], event.calpos_y[it3][0], event.calpos_z[it3][0]);
  TVector3 pi2_end = TVector3(event.calpos_x[it3][pi2_nh-1], event.calpos_y[it3][pi2_nh-1], event.calpos_z[it3][pi2_nh-1]);
  if(!pim_like2){
    pi2_theta_min = TMath::Max(event.helix_t[it3][pi2_nh-1] - vtx_scan_rangeInsidePi/pi2_par[3], event.helix_t[it3][0]);
    pi2_theta_max = event.helix_t[it3][pi2_nh-1] + vtx_scan_range/pi2_par[3];
    pi2_start = TVector3(event.calpos_x[it3][pi2_nh-1], event.calpos_y[it3][pi2_nh-1], event.calpos_z[it3][pi2_nh-1]);
    pi2_end = TVector3(event.calpos_x[it3][0], event.calpos_y[it3][0], event.calpos_z[it3][0]);
  }
  TVector3 pi2_mom; Double_t lpi_dist;
  TVector3 xi_vert = Kinematics::XiVertex(dMagneticField, pi2_par, pi2_theta_min, pi2_theta_max, lambda_vert, lambda_mom, pi2_mom, lpi_dist);
  //std::cout<<"<<xi vertex "<<xi_vert<<std::endl;
  if(xi_vert.z() < -200.) continue; //Vertex cut
  if(TMath::IsNaN(lpi_dist)) continue;

  if(!pim_like) pi2_mom = -1.*pi2_mom;
  TLorentzVector Lpi2(pi2_mom, TMath::Sqrt(pi2_mom.Mag()*pi2_mom.Mag() + PionMass*PionMass));
  TLorentzVector Llambda_fixedmass(lambda_mom, TMath::Sqrt(lambda_mom.Mag()*lambda_mom.Mag() + LambdaMass*LambdaMass));
  TLorentzVector Lxi = Llambda_fixedmass + Lpi2;
  TVector3 xi_mom = TVector3(Lxi.Px(), Lxi.Py(), Lxi.Pz());
  Double_t pi2_vertex_dist;
  //std::cout<<"it3 "<<it3<<" mom "<<pi2_mom.Mag()<<" "<<pi2_mom<<" "<<" lpi_dist "<<lpi_dist<<std::endl;
  if(lpi_dist > lpi_distcut) continue;
  if(!Kinematics::HelixDirection(xi_vert, pi2_start, pi2_end, pi2_vertex_dist)) continue;
  if(pi2_vertex_dist > pi2_vtx_distcut) continue;
  if(TMath::Abs(Lxi.M() - XiMinusMass) > xi_masscut) continue;

  ppi_closedist.push_back(ppi_dist);
  lpi_closedist.push_back(lpi_dist);

  xi_p_container.push_back(it1);
  xi_pi_container.push_back(it2);
  xi_pi2_container.push_back(it3);


  xi_mass_container.push_back(Lxi.M());
  lambda_mass_container.push_back(Llambda.M());
  p_mom_container.push_back(p_mom);
  pi_mom_container.push_back(pi_mom);
  pi2_mom_container.push_back(pi2_mom);
  xi_mom_container.push_back(xi_mom);
  xi_vert_container.push_back(xi_vert);
  l_mom_container.push_back(lambda_mom);
  l_vert_container.push_back(lambda_vert);

  xi_candidates++;
  event.xiflag = true;
      } //it3
    } //it2
  } //it1
if(!event.lflag){
  cout<<"NoLambda"<<endl;
  return true;
} 
if( !event.xiflag ){
  cout<<"NoXi"<<endl;
  Int_t best = -1; Double_t prev_massdiff = 9999.;
  for(Int_t candi=0;candi<l_candidates;candi++){
    Double_t diff = TMath::Abs(l_mass_container_l[candi] - LambdaMass);
    if(prev_massdiff > diff){
      prev_massdiff = diff;
      best = candi;
    }
  }
  event.best = best;
  event.lmass = l_mass_container_l[best];
  event.ldecayvtx_x = l_vert_container_l[best].x();
  event.ldecayvtx_y = l_vert_container_l[best].y();
  event.ldecayvtx_z = l_vert_container_l[best].z();
  event.lmom_x = l_mom_container_l[best].x();
  event.lmom_y = l_mom_container_l[best].y();
  event.lmom_z = l_mom_container_l[best].z();
  event.ppi_dist = ppi_closedist_l[best];
  event.decays_id.push_back(l_p_container[best]);
  event.decays_purity.push_back(event.purity[l_p_container[best]]);
  event.decays_efficiency.push_back(event.efficiency[l_p_container[best]]);
  event.decays_G4tid.push_back(event.G4tid[l_p_container[best]]);
  event.decays_mom.push_back(p_mom_container_l[best].Mag());
  event.decays_mom_x.push_back(p_mom_container_l[best].x());
  event.decays_mom_y.push_back(p_mom_container_l[best].y());
  event.decays_mom_z.push_back(p_mom_container_l[best].z());

  event.decays_id.push_back(l_pi_container[best]);
  event.decays_purity.push_back(event.purity[l_pi_container[best]]);
  event.decays_efficiency.push_back(event.efficiency[l_pi_container[best]]);
  event.decays_G4tid.push_back(event.G4tid[l_pi_container[best]]);
  event.decays_mom.push_back(pi_mom_container_l[best].Mag());
  event.decays_mom_x.push_back(pi_mom_container_l[best].x());
  event.decays_mom_y.push_back(pi_mom_container_l[best].y());
  event.decays_mom_z.push_back(pi_mom_container_l[best].z());

  int id_p = l_p_container[best];
  int id_pi = l_pi_container[best];
  auto Track_p = TPCAna.GetTrackTPCHelix(id_p);
  auto Track_pi = TPCAna.GetTrackTPCHelix(id_pi);
  for(int ih=0;ih<Track_p->GetNHit();++ih){
    auto pos = Track_p->GetHitInOrder(ih)->GetLocalHitPos();
    HF2(1011,pos.z(),pos.x());
  }
  for(int ih=0;ih<Track_pi->GetNHit();++ih){
    auto pos = Track_pi->GetHitInOrder(ih)->GetLocalHitPos();
    HF2(1012,pos.z(),pos.x());
  }
  event.ptid = id_p;
  event.pi1tid = id_pi;
  int G4ptid = G4TrackID.at(id_p);
  int G4pi1tid = G4TrackID.at(id_pi);
  event.G4ptid = G4ptid;
  event.G4pi1tid = G4pi1tid;
  int G4ptnh = PureHits.at(id_p);
  int G4pi1tnh = PureHits.at(id_pi);
  event.G4ptnh = G4ptnh;
  event.G4pi1tnh = G4pi1tnh;

  event.pnh =  event.nhtrack.at(id_p);
  event.pi1nh =  event.nhtrack.at(id_pi);

  event.pmom = p_mom_container_l[best].Mag();
  event.pmom_x = p_mom_container_l[best].x();
  event.pmom_y = p_mom_container_l[best].y();
  event.pmom_z = p_mom_container_l[best].z();

  event.pi1mom = pi_mom_container_l[best].Mag();
  event.pi1mom_x = pi_mom_container_l[best].x();
  event.pi1mom_y = pi_mom_container_l[best].y();
  event.pi1mom_z = pi_mom_container_l[best].z();
#if SkipGenfit
  return true;
#endif
  HF1( 1, event.status++ );
  HypTPCTask& GFTrackCont = HypTPCTask::GetInstance();
  cout<<"Tracks for GF"<<endl;
  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix(it);
    if( !tp ) continue;

//    event.pid[it] = tp -> GetPid();
    std::vector<Int_t> pdgcode;
    Kinematics::HypTPCPID_PDGCode(event.charge[it], event.pid[it], pdgcode);

    GFTrackCont.AddHelixTrack(pdgcode, tp);
  }
  HF1( 2, event.GFstatus++ );
  HF1( 1, event.status++ );
  cout<<"FittingGF"<<endl;
  GFTrackCont.FitTracks();
  cout<<"GF Fitted"<<endl;
  Int_t GFntTpc = GFTrackCont.GetNTrack();
  if(GFntTpc!=event.ntTpc){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    return true;
  }
  HF1( 2, event.GFstatus++ );
  HF1( 1, event.status++ );
  std::vector<Int_t> GFl_p_id_container(l_candidates, -1);
  std::vector<Int_t> GFl_pi_id_container(l_candidates, -1);
  std::vector<Int_t> GFl_p_rep_container(l_candidates, -1);
  std::vector<Int_t> GFl_pi_rep_container(l_candidates, -1);
  std::vector<TVector3> GFl_mom_container(l_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFp_mom_container(l_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFpi_mom_container(l_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFl_vert_container(l_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFl_mass_container(l_candidates, qnan);
  std::vector<Double_t> GFp_mass_container(l_candidates, qnan);
  std::vector<Double_t> GFpi_mass_container(l_candidates, qnan);
  std::vector<Double_t> GFppi_closedist_container(l_candidates, qnan);
  std::vector<Double_t> KFchisqrl_container(l_candidates, qnan); 
  std::vector<Double_t> KFpvall_container(l_candidates, qnan);
  std::vector<std::vector<Double_t>> KFlpull_container(l_candidates, std::vector<Double_t>(6, qnan));  
  std::vector<TVector3> KFp_mom_container(l_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> KFpi_mom_container(l_candidates, TVector3(qnan, qnan, qnan));
  
  Int_t gfbest = -1; prev_massdiff = 9999.;
  for(Int_t candi=0;candi<l_candidates;candi++){
    Int_t trackid_p = l_p_container[candi];
    Int_t trackid_pi = l_pi_container[candi];
    Int_t flag = 1; Int_t repid_p = 0;
    for(Int_t i=0;i<2;i++){
      Int_t temp = flag&event.pid[trackid_p];
      if(temp==flag) repid_p += 1;
      flag*=2;
    }
    Int_t repid_pi = 0;
    Double_t GFextrapolation_decays[3];
    Double_t GFmass2_decays[3] = {qnan, qnan, qnan};
    TVector3 GFmom_decays[3]; TVector3 GFlambda_vert; double GFppi_dist=qnan;
    cout<<"FindingVertex"<<endl;
    if(!GFTrackCont.FindVertex(trackid_p, trackid_pi,
             repid_p, repid_pi,
             GFextrapolation_decays[0], GFextrapolation_decays[1],
             GFmom_decays[0], GFmom_decays[1],
             GFppi_dist, GFlambda_vert,
             vtx_scan_range)
       || GFppi_dist > GFppi_distcut){
          std::cout<<"DistCut! "<<GFppi_dist<<std::endl;
         continue;
       }
    cout<<"VertexFound"<<endl;
    TLorentzVector GFLp(GFmom_decays[0], TMath::Sqrt(GFmom_decays[0].Mag()*GFmom_decays[0].Mag() + ProtonMass*ProtonMass));
    TLorentzVector GFLpi(GFmom_decays[1], TMath::Sqrt(GFmom_decays[1].Mag()*GFmom_decays[1].Mag() + PionMass*PionMass));
    TLorentzVector GFLlambda = GFLp + GFLpi;
    TVector3 GFlambda_mom = GFmom_decays[0] + GFmom_decays[1];
#if DoKinematicFitLdXi
    Double_t KFchisqrl;
    Double_t KFpvall;
    auto Vp = Track_p->GetCovarianceMatrix();
    auto Vpi1 = Track_pi->GetCovarianceMatrix();
    auto HLVP = TLorentzVector(GFLp.X(),GFLp.Z(),GFLp.Y(),GFLp.E());  
    auto HLVPi1 = TLorentzVector(GFLpi.X(),GFLpi.Z(),GFLpi.Y(),GFLpi.E());
    auto HLVLd = TLorentzVector(GFLlambda.X(),GFLlambda.Z(),GFLlambda.Y(),GFLlambda.E());  
    FourVectorFitter KFLd(HLVP,HLVPi1,HLVLd);
    KFLd.SetInvMass(LambdaMass);
    KFLd.SetMaximumStep(5);
    double VarianceLd[6] =
    {Vp(0,0),Vp(1,1),Vp(2,2),
    Vpi1(0,0),Vpi1(1,1),Vpi1(2,2)};
    double OffdiagElemLd[36]={0};
    auto OffdiagLd = MathTools::MergeOffdiagonals(Vp,Vpi1);
    KFLd.SetVariance(VarianceLd);
    KFLd.AddOffdiagonals(OffdiagLd);
    KFchisqrl = KFLd.DoKinematicFit();
    cout<<Form("KFLambda done:: chi2 = %g",KFchisqrl)<<endl;
    KFpvall = KFLd.GetPValue();
    auto HcontLd = KFLd.GetFittedLV();
    auto PullLd = KFLd.GetPull();
    auto KFHLVP = HcontLd.at(0);
    auto KFHLVPi1 = HcontLd.at(1);
    auto KFHLVLd = HcontLd.at(2);
    auto VLd = KFLd.GetUnmeasuredCovariance();
    KFchisqrl_container[candi] = KFchisqrl;
    KFpvall_container[candi] = KFpvall;
    for(int i=0;i<6;i++){
      KFlpull_container[candi][i] = PullLd.at(i);
    }
    KFp_mom_container[candi] = TVector3(KFHLVP.X(),KFHLVP.Z(),KFHLVP.Y());  
    KFpi_mom_container[candi] = TVector3(KFHLVPi1.X(),KFHLVPi1.Z(),KFHLVPi1.Y()); 
    TVector3 KFlambda_mom = TVector3(KFHLVLd.X(),KFHLVLd.Z(),KFHLVLd.Y());  
    GFlambda_mom = KFlambda_mom;   
#endif
    event.GFlflag = true;
    GFl_p_id_container[candi] = trackid_p;
    GFl_pi_id_container[candi] = trackid_pi;
    GFl_p_rep_container[candi] = repid_p;
    GFl_pi_rep_container[candi] = repid_pi;
    GFl_mom_container[candi] = GFlambda_mom;
    GFp_mom_container[candi] = GFmom_decays[0];
    GFpi_mom_container[candi] = GFmom_decays[1];
    GFl_vert_container[candi] = GFlambda_vert;
    GFl_mass_container[candi] = GFLlambda.M();
    GFp_mass_container[candi] = GFmass2_decays[0];
    GFpi_mass_container[candi] = GFmass2_decays[1];
    GFppi_closedist_container[candi] = GFppi_dist;
    Double_t diff = TMath::Abs(GFLlambda.M() - LambdaMass);
    if(prev_massdiff > diff){
      prev_massdiff = diff;
      gfbest = candi;
    }
  }//candi
  std::cout<<"Candidate Loop end"<<std::endl;
  if(!event.GFlflag) return true;
  event.GFntTpc = 2;
  event.GFlmass = GFl_mass_container[gfbest];
  event.GFldecayvtx_x = GFl_vert_container[gfbest].x();
  event.GFldecayvtx_y = GFl_vert_container[gfbest].y();
  event.GFldecayvtx_z = GFl_vert_container[gfbest].z();
  event.GFlmom = GFl_mom_container[gfbest].Mag();
  event.GFlmom_x = GFl_mom_container[gfbest].x();
  event.GFlmom_y = GFl_mom_container[gfbest].y();
  event.GFlmom_z = GFl_mom_container[gfbest].z();
  event.GFppi_dist = GFppi_closedist_container[gfbest];
  event.GFdecays_mass.resize(2);
  event.GFdecays_mom.resize(2);
  event.GFdecays_mom_x.resize(2);
  event.GFdecays_mom_y.resize(2);
  event.GFdecays_mom_z.resize(2);
  event.GFmomloss.resize(2);
  event.GFeloss.resize(2);

  event.GFcharge.resize(2);
  event.GFchisqr.resize(2);
  event.GFtof.resize(2);
  event.GFtracklen.resize(2);
  event.GFpval.resize(2);
  event.GFpdgcode.resize(2);
  event.GFchisqrPos.resize(2);
  event.GFpvalPos.resize(2);
  event.GFfitstatus.resize(2);
  event.GFnhtrack.resize(2);
  event.GFlayer.resize(2);
  event.GFrow.resize(2);
  event.GFpos_x.resize(2);
  event.GFpos_y.resize(2);
  event.GFpos_z.resize(2);
  event.GFmom.resize(2);
  event.GFmom_x.resize(2);
  event.GFmom_y.resize(2);
  event.GFmom_z.resize(2);
  event.GFresidual_x.resize(2);
  event.GFresidual_y.resize(2);
  event.GFresidual_z.resize(2);
  event.GFresidual_p.resize(2);
  event.GFresidual_px.resize(2);
  event.GFresidual_py.resize(2);
  event.GFresidual_pz.resize(2);
  event.GFresolution_x.resize(2);
  event.GFresolution_y.resize(2);
  event.GFresolution_z.resize(2);
  event.GFresolution_p.resize(2);
  event.GFresolution_px.resize(2);
  event.GFresolution_py.resize(2);
  event.GFresolution_pz.resize(2);
  event.GFpull_x.resize(2);
  event.GFpull_y.resize(2);
  event.GFpull_z.resize(2);
  event.GFpull_p.resize(2);
  event.GFpull_px.resize(2);
  event.GFpull_py.resize(2);
  event.GFpull_pz.resize(2);
  for( Int_t j=0; j<2; ++j ){
    Int_t igf = GFl_p_id_container[gfbest];
    if(j==1) igf = GFl_pi_id_container[gfbest];

    event.GFfitstatus[j] = (Int_t) GFTrackCont.TrackCheck(igf);
    HF1( 3, event.GFfitstatus[j]);
    if(!GFTrackCont.TrackCheck(igf)) continue;
    Int_t nh = GFTrackCont.GetNHits(igf);
    Int_t repid = GFl_p_rep_container[gfbest];
    if(j==1) repid = GFl_pi_rep_container[gfbest];

    event.GFnhtrack[j] = nh;
    event.GFchisqr[j] = GFTrackCont.GetChi2NDF(igf, repid);
    event.GFcharge[j] = GFTrackCont.GetCharge(igf, repid);
    event.GFtof[j] = GFTrackCont.GetTrackTOF(igf, 0, -1, repid);
    event.GFtracklen[j] = GFTrackCont.GetTrackLength(igf, 0, -1, repid);
    event.GFpval[j] = GFTrackCont.GetPvalue(igf, repid);
    event.GFpdgcode[j] = GFTrackCont.GetPDGcode(igf, repid);
    HF1(600,event.GFpval[j]);
    TVector3 GFmom_decays = GFp_mom_container[gfbest];
    if(j==1) GFmom_decays = GFpi_mom_container[gfbest];

    Double_t GFmass_decays = GFp_mass_container[gfbest];
    if(j==1) GFmass_decays = GFpi_mass_container[gfbest];

    event.GFdecays_mass[j] = GFmass_decays;
    event.GFdecays_mom[j] = GFmom_decays.Mag();
    event.GFdecays_mom_x[j] = GFmom_decays.x();
    event.GFdecays_mom_y[j] = GFmom_decays.y();
    event.GFdecays_mom_z[j] = GFmom_decays.z();
    event.GFmomloss[j] = GFmom_decays.Mag() - GFTrackCont.GetMom(igf, 0, repid).Mag();
    Double_t pdgmass[3] = {ProtonMass, PionMass, PionMass};
    event.GFeloss[j] = TMath::Sqrt(GFmom_decays.Mag()*GFmom_decays.Mag() + pdgmass[j]*pdgmass[j]) - TMath::Sqrt(GFTrackCont.GetMom(igf, 0, repid).Mag()*GFTrackCont.GetMom(igf, 0, repid).Mag() + pdgmass[j]*pdgmass[j]);

    event.GFlayer[j].resize(nh);
    event.GFrow[j].resize(nh);
    event.GFpos_x[j].resize(nh);
    event.GFpos_y[j].resize(nh);
    event.GFpos_z[j].resize(nh);
    event.GFmom[j].resize(nh);
    event.GFmom_x[j].resize(nh);
    event.GFmom_y[j].resize(nh);
    event.GFmom_z[j].resize(nh);
    event.GFresidual_x[j].resize(nh);
    event.GFresidual_y[j].resize(nh);
    event.GFresidual_z[j].resize(nh);
    event.GFresidual_p[j].resize(nh);
    event.GFresidual_px[j].resize(nh);
    event.GFresidual_py[j].resize(nh);
    event.GFresidual_pz[j].resize(nh);
    event.GFresolution_x[j].resize(nh);
    event.GFresolution_y[j].resize(nh);
    event.GFresolution_z[j].resize(nh);
    event.GFresolution_p[j].resize(nh);
    event.GFresolution_px[j].resize(nh);
    event.GFresolution_py[j].resize(nh);
    event.GFresolution_pz[j].resize(nh);
    event.GFpull_x[j].resize(nh);
    event.GFpull_y[j].resize(nh);
    event.GFpull_z[j].resize(nh);
    event.GFpull_p[j].resize(nh);
    event.GFpull_px[j].resize(nh);
    event.GFpull_py[j].resize(nh);
    event.GFpull_pz[j].resize(nh);

    Int_t id = GFl_p_id_container[gfbest];
    if(j==1) id = GFl_pi_id_container[gfbest];
    event.GFpmom = GFp_mom_container[gfbest].Mag();
    event.GFpmom_x = GFp_mom_container[gfbest].x();
    event.GFpmom_y = GFp_mom_container[gfbest].y();
    event.GFpmom_z = GFp_mom_container[gfbest].z();
    event.GFpi1mom = GFpi_mom_container[gfbest].Mag();
    event.GFpi1mom_x = GFpi_mom_container[gfbest].x();
    event.GFpi1mom_y = GFpi_mom_container[gfbest].y();
    event.GFpi1mom_z = GFpi_mom_container[gfbest].z();



    int ihit = 0;
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( id );

    double GFchisqrPos=0;
    for( Int_t ih=0; ih<tp->GetNHit(); ++ih ){
      Int_t layer = (Int_t) event.hitlayer[id][ih];
      TPCLTrackHit *helix_point = tp -> GetHitInOrder(ih);
      if(!helix_point -> IsGoodForTracking()) continue;
      const TVector3 &hit0 = helix_point -> GetLocalHitPos();
      TVector3 mom0 = helix_point -> GetMomentumHelix(event.GFcharge[j]);
      double row = helix_point->GetMRow()+0.5;
      TVector3 hit = GFTrackCont.GetPos(igf, ihit, repid);
      TVector3 mom = GFTrackCont.GetMom(igf, ihit, repid);
      double residual_[5];
      double pull_[5];
      double GFresidual6D[6];
      double GFresolution6D[6];
      double GFpull6D[6];
      TVector3 dumV;
      double dumd;
      GFTrackCont.GetTrackPull(igf,event.GFpdgcode[j],dumV,
             dumd,mom0,hit0,residual_,
             pull_,GFresidual6D,GFpull6D);

      event.GFlayer[j][ihit] = layer;
      event.GFrow[j][ihit] = row;
      event.GFmom_x[j][ihit] = mom.x();
      event.GFmom_y[j][ihit] = mom.y();
      event.GFmom_z[j][ihit] = mom.z();
      event.GFmom[j][ihit] = mom.Mag();
      event.GFpos_x[j][ihit] = hit.x();
      event.GFpos_y[j][ihit] = hit.y();
      event.GFpos_z[j][ihit] = hit.z();
      event.GFresidual_x[j][ihit] = hit.x() - hit0.x();
      event.GFresidual_y[j][ihit] = hit.y() - hit0.y();
      event.GFresidual_z[j][ihit] = hit.z() - hit0.z();
      event.GFresidual_p[j][ihit] = mom.Mag() - mom0.Mag();
      event.GFresidual_px[j][ihit] = mom.x() - mom0.x();
      event.GFresidual_py[j][ihit] = mom.y() - mom0.y();
      event.GFresidual_pz[j][ihit] = mom.z() - mom0.z();

      event.GFresolution_x[j][ihit] = GFresidual6D[0];
       event.GFresolution_y[j][ihit] = GFresidual6D[1];
       event.GFresolution_z[j][ihit] = GFresidual6D[2];
       event.GFresolution_px[j][ihit] = GFresidual6D[3];
       event.GFresolution_py[j][ihit] = GFresidual6D[4];
       event.GFresolution_pz[j][ihit] = GFresidual6D[5];
      event.GFpull_x[j][ihit] = GFpull6D[0];
       event.GFpull_y[j][ihit] = GFpull6D[1];
       event.GFpull_z[j][ihit] = GFpull6D[2];
       event.GFpull_px[j][ihit] = GFpull6D[3];
       event.GFpull_py[j][ihit] = GFpull6D[4];
       event.GFpull_pz[j][ihit] = GFpull6D[5];
      HF1(601,event.GFpull_x[j][ihit]);
      HF1(602,event.GFpull_y[j][ihit]);
      HF1(603,event.GFpull_z[j][ihit]);
      HF1(604,event.GFpull_px[j][ihit]);
      HF1(605,event.GFpull_py[j][ihit]);
      HF1(606,event.GFpull_pz[j][ihit]);
      if(event.GFpval[j]>0.01){
        HF1(611,event.GFpull_x[j][ihit]);
        HF1(612,event.GFpull_y[j][ihit]);
        HF1(613,event.GFpull_z[j][ihit]);
        HF1(614,event.GFpull_px[j][ihit]);
        HF1(615,event.GFpull_py[j][ihit]);
        HF1(616,event.GFpull_pz[j][ihit]);
      }
      ihit++;
    } //ih

    double GFndf = 2*ihit - 5; //Effective number of clusters
    double GFpvalPos = -1;
    if(GFndf > 0){
      GFchisqrPos /= GFndf;
      GFpvalPos = 1-ROOT::Math::chisquared_cdf(GFchisqrPos*GFndf, GFndf);
    }
    else GFchisqrPos = -1;

    event.GFchisqrPos[j]=GFchisqrPos;
    event.GFpvalPos[j]=GFpvalPos;
  } //igf
#if DoKinematicFitLdXi
  Double_t KFchisqrl = KFchisqrl_container[gfbest];
  Double_t KFpvall = KFpvall_container[gfbest];
  std::vector<Double_t> PullLd = KFlpull_container[gfbest];  
  TVector3 KFTVP = KFp_mom_container[gfbest];
  TVector3 KFTVPi1 = KFpi_mom_container[gfbest];
  TVector3 KFTVLd = KFTVP+KFTVPi1;  
  
  event.KFchisqrl = KFchisqrl;
  event.KFpvall = KFpvall;
  event.KFlmom = KFTVLd.Mag();
  event.KFlmom_x = KFTVLd.x();
  event.KFlmom_y = KFTVLd.y();
  event.KFlmom_z = KFTVLd.z();
  event.KFdecays_mom.push_back(KFTVP.Mag());
  event.KFdecays_mom_x.push_back(KFTVP.x());
  event.KFdecays_mom_y.push_back(KFTVP.y());
  event.KFdecays_mom_z.push_back(KFTVP.z());
  event.KFdecays_mom.push_back(KFTVPi1.Mag());
  event.KFdecays_mom_x.push_back(KFTVPi1.x());
  event.KFdecays_mom_y.push_back(KFTVPi1.y());
  event.KFdecays_mom_z.push_back(KFTVPi1.z());
#endif

  return true;
}
#if 1 //Select the best L mass combination
  Int_t best = -1; Double_t prev_massdiff = 9999.;
  for(Int_t candi=0;candi<xi_candidates;candi++){
    Double_t diff = TMath::Abs(lambda_mass_container[candi] - LambdaMass);
    if(prev_massdiff > diff){
      prev_massdiff = diff;
      best = candi;
    }
  }
#else //Select the best Xi mass combination
  Int_t best = -1; Double_t prev_massdiff = 9999.;
  for(Int_t candi=0;candi<xi_candidates;candi++){
    Double_t diff = TMath::Abs(xi_mass_container[candi] - XiMinusMass);
    if(prev_massdiff > diff){
      prev_massdiff = diff;
      best = candi;
    }
  }
#endif
  


  event.best = best;
  event.ximass = xi_mass_container[best];
  event.xidecayvtx_x = xi_vert_container[best].x();
  event.xidecayvtx_y = xi_vert_container[best].y();
  event.xidecayvtx_z = xi_vert_container[best].z();
  event.ximom_x = xi_mom_container[best].x();
  event.ximom_y = xi_mom_container[best].y();
  event.ximom_z = xi_mom_container[best].z();
  event.lpi_dist = lpi_closedist[best];

  event.lmass = lambda_mass_container[best];
  event.ldecayvtx_x = l_vert_container[best].x();
  event.ldecayvtx_y = l_vert_container[best].y();
  event.ldecayvtx_z = l_vert_container[best].z();
  event.lmom_x = l_mom_container[best].x();
  event.lmom_y = l_mom_container[best].y();
  event.lmom_z = l_mom_container[best].z();
  event.ppi_dist = ppi_closedist[best];

  event.decays_id.push_back(xi_p_container[best]);
  event.decays_purity.push_back(event.purity[xi_p_container[best]]);
  event.decays_efficiency.push_back(event.efficiency[xi_p_container[best]]);
  event.decays_G4tid.push_back(event.G4tid[xi_p_container[best]]);
  event.decays_mom.push_back(p_mom_container[best].Mag());
  event.decays_mom_x.push_back(p_mom_container[best].x());
  event.decays_mom_y.push_back(p_mom_container[best].y());
  event.decays_mom_z.push_back(p_mom_container[best].z());


  event.decays_id.push_back(xi_pi_container[best]);
  event.decays_purity.push_back(event.purity[xi_pi_container[best]]);
  event.decays_efficiency.push_back(event.efficiency[xi_pi_container[best]]);
  event.decays_G4tid.push_back(event.G4tid[xi_pi_container[best]]);
  event.decays_mom.push_back(pi_mom_container[best].Mag());
  event.decays_mom_x.push_back(pi_mom_container[best].x());
  event.decays_mom_y.push_back(pi_mom_container[best].y());
  event.decays_mom_z.push_back(pi_mom_container[best].z());


  event.decays_id.push_back(xi_pi2_container[best]);
  event.decays_purity.push_back(event.purity[xi_pi2_container[best]]);
  event.decays_efficiency.push_back(event.efficiency[xi_pi2_container[best]]);
  event.decays_G4tid.push_back(event.G4tid[xi_pi2_container[best]]);
  event.decays_mom.push_back(pi2_mom_container[best].Mag());
  event.decays_mom_x.push_back(pi2_mom_container[best].x());
  event.decays_mom_y.push_back(pi2_mom_container[best].y());
  event.decays_mom_z.push_back(pi2_mom_container[best].z());

  int id_p = xi_p_container[best];
  int id_pi = xi_pi_container[best];
  int id_pi2 = xi_pi2_container[best];
  auto Track_p = TPCAna.GetTrackTPCHelix(id_p);
  auto Track_pi = TPCAna.GetTrackTPCHelix(id_pi);
  auto Track_pi2 = TPCAna.GetTrackTPCHelix(id_pi2);
  for(int ih=0;ih<Track_p->GetNHit();++ih){
    auto pos = Track_p->GetHitInOrder(ih)->GetLocalHitPos();
    HF2(1001,pos.z(),pos.x());
  }
  for(int ih=0;ih<Track_pi->GetNHit();++ih){
    auto pos = Track_pi->GetHitInOrder(ih)->GetLocalHitPos();
    HF2(1002,pos.z(),pos.x());
  }
  for(int ih=0;ih<Track_pi2->GetNHit();++ih){
    auto pos = Track_pi2->GetHitInOrder(ih)->GetLocalHitPos();
    HF2(1003,pos.z(),pos.x());
  }

  event.decays_res_mom.push_back(Track_p->GetMomentumResolution());
  event.decays_res_mom_x.push_back(Track_p->GetMomentumResolutionVect().X());
  event.decays_res_mom_y.push_back(Track_p->GetMomentumResolutionVect().Y());
  event.decays_res_mom_z.push_back(Track_p->GetMomentumResolutionVect().Z());
  event.decays_res_mom_t.push_back(Track_p->GetTransverseMomentumResolution());
  event.decays_res_th.push_back(Track_p->GetThetaResolution());
  event.decays_res_ph.push_back(Track_p->GetTransverseAngularResolution());

  event.decays_cov_mom_ph.push_back(Track_p->GetTransverseMomentumAngularCovariance());
  event.decays_cov_mom_th.push_back(Track_p->GetMomentumPitchAngleCovariance());
  event.decays_cov_mom_xy.push_back(Track_p->GetMomentumCovarianceVect().x());
  event.decays_cov_mom_yz.push_back(Track_p->GetMomentumCovarianceVect().y());
  event.decays_cov_mom_zx.push_back(Track_p->GetMomentumCovarianceVect().z());

  event.decays_res_mom.push_back(Track_pi->GetMomentumResolution());
  event.decays_res_mom_x.push_back(Track_pi->GetMomentumResolutionVect().X());
  event.decays_res_mom_y.push_back(Track_pi->GetMomentumResolutionVect().Y());
  event.decays_res_mom_z.push_back(Track_pi->GetMomentumResolutionVect().Z());
  event.decays_res_mom_t.push_back(Track_pi->GetTransverseMomentumResolution());
  event.decays_res_th.push_back(Track_pi->GetThetaResolution());
  event.decays_res_ph.push_back(Track_pi->GetTransverseAngularResolution());

  event.decays_cov_mom_ph.push_back(Track_pi->GetTransverseMomentumAngularCovariance());
  event.decays_cov_mom_th.push_back(Track_pi->GetMomentumPitchAngleCovariance());
  event.decays_cov_mom_xy.push_back(Track_pi->GetMomentumCovarianceVect().x());
  event.decays_cov_mom_yz.push_back(Track_pi->GetMomentumCovarianceVect().y());
  event.decays_cov_mom_zx.push_back(Track_pi->GetMomentumCovarianceVect().z());

  event.decays_res_mom.push_back(Track_pi2->GetMomentumResolution());
  event.decays_res_mom_x.push_back(Track_pi2->GetMomentumResolutionVect().X());
  event.decays_res_mom_y.push_back(Track_pi2->GetMomentumResolutionVect().Y());
  event.decays_res_mom_z.push_back(Track_pi2->GetMomentumResolutionVect().Z());
  event.decays_res_mom_t.push_back(Track_pi2->GetTransverseMomentumResolution());
  event.decays_res_th.push_back(Track_pi2->GetThetaResolution());
  event.decays_res_ph.push_back(Track_pi2->GetTransverseAngularResolution());

  event.decays_cov_mom_ph.push_back(Track_pi2->GetTransverseMomentumAngularCovariance());
  event.decays_cov_mom_th.push_back(Track_pi2->GetMomentumPitchAngleCovariance());
  event.decays_cov_mom_xy.push_back(Track_pi2->GetMomentumCovarianceVect().x());
  event.decays_cov_mom_yz.push_back(Track_pi2->GetMomentumCovarianceVect().y());
  event.decays_cov_mom_zx.push_back(Track_pi2->GetMomentumCovarianceVect().z());


  event.ptid = id_p;
  event.pi1tid = id_pi;
  event.pi2tid = id_pi2;


  int G4ptid = G4TrackID.at(id_p);
  int G4pi1tid = G4TrackID.at(id_pi);
  int G4pi2tid = G4TrackID.at(id_pi2);
  event.G4ptid = G4ptid;
  event.G4pi1tid = G4pi1tid;
  event.G4pi2tid = G4pi2tid;

  int G4ptnh = PureHits.at(id_p);
  int G4pi1tnh = PureHits.at(id_pi);
  int G4pi2tnh = PureHits.at(id_pi2);
  event.G4ptnh = G4ptnh;
  event.G4pi1tnh = G4pi1tnh;
  event.G4pi2tnh = G4pi2tnh;

  event.pnh =  event.nhtrack.at(id_p);
  event.pi1nh =  event.nhtrack.at(id_pi);
  event.pi2nh =  event.nhtrack.at(id_pi2);

  event.pmom = p_mom_container[best].Mag();
  event.pmom_x = p_mom_container[best].x();
  event.pmom_y = p_mom_container[best].y();
  event.pmom_z = p_mom_container[best].z();

  event.pi1mom = pi_mom_container[best].Mag();
  event.pi1mom_x = pi_mom_container[best].x();
  event.pi1mom_y = pi_mom_container[best].y();
  event.pi1mom_z = pi_mom_container[best].z();

  event.pi2mom = pi2_mom_container[best].Mag();
  event.pi2mom_x = pi2_mom_container[best].x();
  event.pi2mom_y = pi2_mom_container[best].y();
  event.pi2mom_z = pi2_mom_container[best].z();


  HF1( 11, event.lmass );
  HF1( 12, event.ximass );
  HF1( 13, event.decays_mom[0] );
  HF1( 14, event.decays_mom[1] );
  HF1( 15, event.decays_mom[2] );
  HF1( 16, TMath::Sqrt( event.lmom_x*event.lmom_x +
      event.lmom_y*event.lmom_y +
      event.lmom_z*event.lmom_z) );
  HF1( 17, TMath::Sqrt( event.ximom_x*event.ximom_x +
      event.ximom_y*event.ximom_y +
      event.ximom_z*event.ximom_z) );
  HF1( 18, event.ppi_dist);
  HF1( 19, event.lpi_dist);

  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix(it);
    if( !tp ) continue;
    for(int ih=0;ih<event.hitpos_x.at(it).size();++ih){
      auto hit = tp->GetHitInOrder(ih);
      auto hitpos = hit->GetLocalHitPos();
      double t = hit-> GetTheta();
      TVector3 mom = hit -> GetMomentumHelix(tp->GetCharge());
      TVector3 mom0 = GetG4Mom(hitpos,G4Hits,G4Moms);
      auto MomResVect = tp->GetMomentumResolutionVect(ih);
      auto MomResiVect = mom - mom0;
      HF1(101,MomResiVect.x());
      HF1(102,MomResVect.x());
      HF1(103,MomResiVect.x()/MomResVect.x());
      HF1(201,MomResiVect.y());
      HF1(202,MomResVect.y());
      HF1(203,MomResiVect.y()/MomResVect.y());
      HF1(301,MomResiVect.z());
      HF1(302,MomResVect.z());
      HF1(303,MomResiVect.z()/MomResVect.z());
      HF1(400,tp->GetTransverseMomentumResolution());
      HF1(401,tp->GetTransverseAngularResolution(t));
      HF1(402,tp->GetdZResolution());
      auto ResVect = tp->GetResolutionVect(ih,false);
      HF1(500,ResVect.Mag());
      HF1(501,ResVect.X());
      HF1(502,ResVect.Y());
      HF1(503,ResVect.Z());
    }
  }
#if SkipGenfit
  return true;
#endif
  HF1( 1, event.status++ );

  HypTPCTask& GFTrackCont = HypTPCTask::GetInstance();
  cout<<"Tracks for GF"<<endl;
  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix(it);
    if( !tp ) continue;

//    event.pid[it] = tp -> GetPid();
    std::vector<Int_t> pdgcode;
    Kinematics::HypTPCPID_PDGCode(event.charge[it], event.pid[it], pdgcode);

    GFTrackCont.AddHelixTrack(pdgcode, tp);
  }
  HF1( 2, event.GFstatus++ );
  HF1( 1, event.status++ );
  GFTrackCont.FitTracks();

  Int_t GFntTpc = GFTrackCont.GetNTrack();
  if(GFntTpc!=event.ntTpc){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    return true;
  }

  HF1( 2, event.GFstatus++ );
  HF1( 1, event.status++ );

  std::vector<Int_t> GFxi_p_id_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_pi_id_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_pi2_id_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_p_rep_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_pi_rep_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_pi2_rep_container(xi_candidates, -1);
  std::vector<TVector3> GFxi_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFl_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFp_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFpi_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFpi2_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFxi_vert_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFl_vert_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFxi_mass_container(xi_candidates, qnan);
  std::vector<Double_t> GFl_mass_container(xi_candidates, qnan);
  std::vector<Double_t> GFp_mass_container(xi_candidates, qnan);
  std::vector<Double_t> GFpi_mass_container(xi_candidates, qnan);
  std::vector<Double_t> GFpi2_mass_container(xi_candidates, qnan);
  std::vector<Double_t> GFppi_closedist_container(xi_candidates, qnan);
  std::vector<Double_t> GFlpi_closedist_container(xi_candidates, qnan);
  std::vector<Double_t> KFlpi_closedist_container(xi_candidates, qnan);
  std::vector<Double_t> KFchisqrl_container(xi_candidates, qnan); 
  std::vector<Double_t> KFpvall_container(xi_candidates, qnan);
  std::vector<std::vector<Double_t>> KFlpull_container(xi_candidates, std::vector<Double_t>(6, qnan));  
  std::vector<TMatrixD> KFVarianceLd_container(xi_candidates, TMatrixD(3,3)); 
  std::vector<TVector3> KFp_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> KFpi_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));

  Int_t gfbest = -1; prev_massdiff = 9999.;
  for(Int_t candi=0;candi<xi_candidates;candi++){
    Int_t trackid_p = xi_p_container[candi];
    Int_t trackid_pi = xi_pi_container[candi];
    Int_t trackid_pi2 = xi_pi2_container[candi];

    Int_t flag = 1; Int_t repid_p = 0;

    for(Int_t i=0;i<2;i++){
      Int_t temp = flag&event.pid[trackid_p];
      if(temp==flag) repid_p += 1;
      flag*=2;
    }

    Int_t repid_pi = 0; Int_t repid_pi2 = 0;

    Double_t GFextrapolation_decays[3];
    Double_t GFmass2_decays[3] = {qnan, qnan, qnan};
    TVector3 GFmom_decays[3]; TVector3 GFlambda_vert; double GFppi_dist=qnan;
    if(!GFTrackCont.FindVertex(trackid_p, trackid_pi,
             repid_p, repid_pi,
             GFextrapolation_decays[0], GFextrapolation_decays[1],
             GFmom_decays[0], GFmom_decays[1],
             GFppi_dist, GFlambda_vert,
             vtx_scan_range)
       || GFppi_dist > GFppi_distcut) continue;

    TLorentzVector GFLp(GFmom_decays[0], TMath::Sqrt(GFmom_decays[0].Mag()*GFmom_decays[0].Mag() + ProtonMass*ProtonMass));
    TLorentzVector GFLpi(GFmom_decays[1], TMath::Sqrt(GFmom_decays[1].Mag()*GFmom_decays[1].Mag() + PionMass*PionMass));
    TLorentzVector GFLlambda = GFLp + GFLpi;
    TVector3 GFlambda_mom = GFmom_decays[0] + GFmom_decays[1];
#if DoKinematicFitLdXi
  Double_t KFchisqrl;
  Double_t KFpvall;
  auto Vp = Track_p->GetCovarianceMatrix();
  auto Vpi1 = Track_pi->GetCovarianceMatrix();
  auto HLVP = TLorentzVector(GFLp.X(),GFLp.Z(),GFLp.Y(),GFLp.E());  
  auto HLVPi1 = TLorentzVector(GFLpi.X(),GFLpi.Z(),GFLpi.Y(),GFLpi.E());
  auto HLVLd = TLorentzVector(GFLlambda.X(),GFLlambda.Z(),GFLlambda.Y(),GFLlambda.E());  
  FourVectorFitter KFLd(HLVP,HLVPi1,HLVLd);
  KFLd.SetInvMass(LambdaMass);
  KFLd.SetMaximumStep(5);
  double VarianceLd[6] =
  {Vp(0,0),Vp(1,1),Vp(2,2),
  Vpi1(0,0),Vpi1(1,1),Vpi1(2,2)};
  double OffdiagElemLd[36]={0};
  auto OffdiagLd = MathTools::MergeOffdiagonals(Vp,Vpi1);
  KFLd.SetVariance(VarianceLd);
  KFLd.AddOffdiagonals(OffdiagLd);
  KFchisqrl = KFLd.DoKinematicFit();
  cout<<Form("KFLambda done:: chi2 = %g",KFchisqrl)<<endl;
  KFpvall = KFLd.GetPValue();
  auto HcontLd = KFLd.GetFittedLV();
  auto PullLd = KFLd.GetPull();
  auto KFHLVP = HcontLd.at(0);
  auto KFHLVPi1 = HcontLd.at(1);
  auto KFHLVLd = HcontLd.at(2);
  auto VLd = KFLd.GetUnmeasuredCovariance();
  KFchisqrl_container[candi] = KFchisqrl;
  KFpvall_container[candi] = KFpvall;
  for(int i=0;i<6;i++){
    KFlpull_container[candi][i] = PullLd.at(i);
  }
  KFVarianceLd_container[candi] = VLd;
  KFp_mom_container[candi] = TVector3(KFHLVP.X(),KFHLVP.Z(),KFHLVP.Y());  
  KFpi_mom_container[candi] = TVector3(KFHLVPi1.X(),KFHLVPi1.Z(),KFHLVPi1.Y()); 
  TVector3 GFxi_vert; Double_t GFlpi_dist = qnan; Double_t GFlambda_tracklen;
  if(!GFTrackCont.FindVertexXi(trackid_pi2, repid_pi2, GFlambda_vert, GFlambda_mom,
  GFlambda_tracklen, GFextrapolation_decays[2], GFmom_decays[2], GFlpi_dist, GFxi_vert, vtx_scan_range) 
    || GFlpi_dist > GFlpi_distcut) continue;
  TVector3 KFlambda_mom = TVector3(KFHLVLd.X(),KFHLVLd.Z(),KFHLVLd.Y());  
  GFlambda_mom = KFlambda_mom;   
#endif
    TLorentzVector GFLlambda_fixed(GFlambda_mom, TMath::Sqrt(GFlambda_mom.Mag()*GFlambda_mom.Mag() + LambdaMass*LambdaMass));
    TVector3 KFxi_vert; Double_t KFlpi_dist = qnan; Double_t KFlambda_tracklen;
    #if DoKinematicFitLdXi
    double l_res_x,l_res_y,l_phi;
    MathTools::DecomposeResolution(VLd,KFlambda_mom,l_res_x,l_res_y,l_phi); 
    if(!GFTrackCont.FindVertexXi(trackid_pi2, repid_pi2, GFlambda_vert, GFlambda_mom,
    KFlambda_tracklen, GFextrapolation_decays[2], GFmom_decays[2], GFlpi_dist, KFxi_vert, vtx_scan_range,l_res_x,l_res_y,l_phi) 
      || KFlpi_dist > GFlpi_distcut) continue;
    #else
    if(!GFTrackCont.FindVertexXi(trackid_pi2, repid_pi2,
         GFlambda_vert, GFlambda_mom, GFlambda_tracklen,
         GFextrapolation_decays[0], GFmom_decays[2],
         GFlpi_dist, GFxi_vert,
         vtx_scan_range)
       || GFlpi_dist > GFlpi_distcut) continue;
    #endif
    Double_t GFlambda_tof = Kinematics::CalcTimeOfFlight(GFlambda_mom.Mag(), KFlambda_tracklen, pdg::LambdaMass());

    Int_t htofhitid_p; Double_t tracklen_p; Double_t tof; TVector3 pos; Double_t track2tgt_dist;
    Bool_t htofextrapolation_p = GFTrackCont.TPCHTOFTrackMatching(trackid_p, repid_p, GFlambda_vert, event.HtofSeg, event.posHtof, htofhitid_p, tof, tracklen_p, pos, track2tgt_dist);
    if(htofextrapolation_p){
      GFmass2_decays[0] = Kinematics::MassSquare(GFmom_decays[0].Mag(),
						 tracklen_p  - GFextrapolation_decays[0],
						 event.tHtof[htofhitid_p] - GFlambda_tof);
      //if(GFmass2_decays[0] < 0.25) continue;
    }

    Int_t htofhitid_pi; Double_t tracklen_pi;
    Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(trackid_pi, repid_pi, GFlambda_vert, event.HtofSeg, event.posHtof,htofhitid_pi, tof, tracklen_pi, pos, track2tgt_dist);
    if(htofextrapolation_pi){
      GFmass2_decays[1] = Kinematics::MassSquare(GFmom_decays[1].Mag(),
						 tracklen_pi - GFextrapolation_decays[1],
						 event.tHtof[htofhitid_pi] - GFlambda_tof);
      //if(GFmass2_decays[1] > 0.25) continue;
    }

    TLorentzVector GFLpi2(GFmom_decays[2], TMath::Sqrt(GFmom_decays[2].Mag()*GFmom_decays[2].Mag() + PionMass*PionMass));
    TLorentzVector GFLxi = GFLlambda_fixed + GFLpi2;
    TVector3 GFxi_mom = GFlambda_mom + GFmom_decays[2];

    Int_t htofhitid_pi2; Double_t tracklen_pi2;
    Bool_t htofextrapolation_pi2 = GFTrackCont.TPCHTOFTrackMatching(trackid_pi2, repid_pi2, tgtpos, event.HtofSeg, event.posHtof, htofhitid_pi2, tof, tracklen_pi2, pos, track2tgt_dist);
    if(htofextrapolation_pi2){
      GFmass2_decays[2] = Kinematics::MassSquare(GFmom_decays[2].Mag(),
						 tracklen_pi2 - GFextrapolation_decays[2],
						 event.tHtof[htofhitid_pi2]);
      //if(GFmass2_decays[2] > 0.25) continue;
    }

    event.GFxiflag = true;
    GFxi_p_id_container[candi] = trackid_p;
    GFxi_pi_id_container[candi] = trackid_pi;
    GFxi_pi2_id_container[candi] = trackid_pi2;
    GFxi_p_rep_container[candi] = repid_p;
    GFxi_pi_rep_container[candi] = repid_pi;
    GFxi_pi2_rep_container[candi] = repid_pi2;
    GFxi_mom_container[candi] = GFxi_mom;
    GFl_mom_container[candi] = GFlambda_mom;
    GFp_mom_container[candi] = GFmom_decays[0];
    GFpi_mom_container[candi] = GFmom_decays[1];
    GFpi2_mom_container[candi] = GFmom_decays[2];
    GFxi_vert_container[candi] = KFxi_vert;
    GFl_vert_container[candi] = GFlambda_vert;
    GFxi_mass_container[candi] = GFLxi.M();
    GFl_mass_container[candi] = GFLlambda.M();
    GFp_mass_container[candi] = GFmass2_decays[0];
    GFpi_mass_container[candi] = GFmass2_decays[1];
    GFpi2_mass_container[candi] = GFmass2_decays[2];
    GFppi_closedist_container[candi] = GFppi_dist;
    GFlpi_closedist_container[candi] = GFlpi_dist;
    KFlpi_closedist_container[candi] = KFlpi_dist;

    Double_t diff = TMath::Abs(GFLxi.M() - XiMinusMass);
    if(prev_massdiff > diff){
      prev_massdiff = diff;
      gfbest = candi;
    }
  }//candi

  if(!event.GFxiflag) return true;
  HF1( genfitHid, GFntTpc);

  HF1( 21, event.lmass );
  HF1( 22, event.ximass );
  HF1( 23, event.decays_mom[0] );
  HF1( 24, event.decays_mom[1] );
  HF1( 25, event.decays_mom[2] );
  HF1( 26, TMath::Sqrt( event.lmom_x*event.lmom_x +
      event.lmom_y*event.lmom_y +
      event.lmom_z*event.lmom_z) );
  HF1( 27, TMath::Sqrt( event.ximom_x*event.ximom_x +
      event.ximom_y*event.ximom_y +
      event.ximom_z*event.ximom_z) );
  HF1( 28, event.ppi_dist);
  HF1( 29, event.lpi_dist);

  event.GFntTpc = 3;
  event.GFximass = GFxi_mass_container[gfbest];
  event.GFxidecayvtx_x = GFxi_vert_container[gfbest].x();
  event.GFxidecayvtx_y = GFxi_vert_container[gfbest].y();
  event.GFxidecayvtx_z = GFxi_vert_container[gfbest].z();
  event.GFximom = GFxi_mom_container[gfbest].Mag();
  event.GFximom_x = GFxi_mom_container[gfbest].x();
  event.GFximom_y = GFxi_mom_container[gfbest].y();
  event.GFximom_z = GFxi_mom_container[gfbest].z();
  event.GFlpi_dist = GFlpi_closedist_container[gfbest];
  event.KFlpi_dist = KFlpi_closedist_container[gfbest];

  event.GFlmass = GFl_mass_container[gfbest];
  event.GFldecayvtx_x = GFl_vert_container[gfbest].x();
  event.GFldecayvtx_y = GFl_vert_container[gfbest].y();
  event.GFldecayvtx_z = GFl_vert_container[gfbest].z();
  event.GFlmom = GFl_mom_container[gfbest].Mag();
  event.GFlmom_x = GFl_mom_container[gfbest].x();
  event.GFlmom_y = GFl_mom_container[gfbest].y();
  event.GFlmom_z = GFl_mom_container[gfbest].z();
  event.GFppi_dist = GFppi_closedist_container[gfbest];

  event.GFdecays_mass.resize(3);
  event.GFdecays_mom.resize(3);
  event.GFdecays_mom_x.resize(3);
  event.GFdecays_mom_y.resize(3);
  event.GFdecays_mom_z.resize(3);
  event.GFmomloss.resize(3);
  event.GFeloss.resize(3);

  event.GFcharge.resize(3);
  event.GFchisqr.resize(3);
  event.GFtof.resize(3);
  event.GFtracklen.resize(3);
  event.GFpval.resize(3);
  event.GFpdgcode.resize(3);
  event.GFchisqrPos.resize(3);
  event.GFpvalPos.resize(3);

  event.GFfitstatus.resize(3);
  event.GFnhtrack.resize(3);
  event.GFlayer.resize(3);
  event.GFrow.resize(3);
  event.GFpos_x.resize(3);
  event.GFpos_y.resize(3);
  event.GFpos_z.resize(3);
  event.GFmom.resize(3);
  event.GFmom_x.resize(3);
  event.GFmom_y.resize(3);
  event.GFmom_z.resize(3);
  event.GFresidual_x.resize(3);
  event.GFresidual_y.resize(3);
  event.GFresidual_z.resize(3);
  event.GFresidual_p.resize(3);
  event.GFresidual_px.resize(3);
  event.GFresidual_py.resize(3);
  event.GFresidual_pz.resize(3);
  event.GFresolution_x.resize(3);
  event.GFresolution_y.resize(3);
  event.GFresolution_z.resize(3);
  event.GFresolution_p.resize(3);
  event.GFresolution_px.resize(3);
  event.GFresolution_py.resize(3);
  event.GFresolution_pz.resize(3);
  event.GFpull_x.resize(3);
  event.GFpull_y.resize(3);
  event.GFpull_z.resize(3);
  event.GFpull_p.resize(3);
  event.GFpull_px.resize(3);
  event.GFpull_py.resize(3);
  event.GFpull_pz.resize(3);
  for( Int_t j=0; j<3; ++j ){
    Int_t igf = GFxi_p_id_container[gfbest];
    if(j==1) igf = GFxi_pi_id_container[gfbest];
    if(j==2) igf = GFxi_pi2_id_container[gfbest];

    event.GFfitstatus[j] = (Int_t) GFTrackCont.TrackCheck(igf);
    HF1( 3, event.GFfitstatus[j]);
    if(!GFTrackCont.TrackCheck(igf)) continue;
    Int_t nh = GFTrackCont.GetNHits(igf);
    Int_t repid = GFxi_p_rep_container[gfbest];
    if(j==1) repid = GFxi_pi_rep_container[gfbest];
    if(j==2) repid = GFxi_pi2_rep_container[gfbest];

    event.GFnhtrack[j] = nh;
    event.GFchisqr[j] = GFTrackCont.GetChi2NDF(igf, repid);
    event.GFcharge[j] = GFTrackCont.GetCharge(igf, repid);
    event.GFtof[j] = GFTrackCont.GetTrackTOF(igf, 0, -1, repid);
    event.GFtracklen[j] = GFTrackCont.GetTrackLength(igf, 0, -1, repid);
    event.GFpval[j] = GFTrackCont.GetPvalue(igf, repid);
    event.GFpdgcode[j] = GFTrackCont.GetPDGcode(igf, repid);
    HF1(600,event.GFpval[j]);
    TVector3 GFmom_decays = GFp_mom_container[gfbest];
    if(j==1) GFmom_decays = GFpi_mom_container[gfbest];
    if(j==2) GFmom_decays = GFpi2_mom_container[gfbest];

    Double_t GFmass_decays = GFp_mass_container[gfbest];
    if(j==1) GFmass_decays = GFpi_mass_container[gfbest];
    if(j==2) GFmass_decays = GFpi2_mass_container[gfbest];

    event.GFdecays_mass[j] = GFmass_decays;
    event.GFdecays_mom[j] = GFmom_decays.Mag();
    event.GFdecays_mom_x[j] = GFmom_decays.x();
    event.GFdecays_mom_y[j] = GFmom_decays.y();
    event.GFdecays_mom_z[j] = GFmom_decays.z();
    event.GFmomloss[j] = GFmom_decays.Mag() - GFTrackCont.GetMom(igf, 0, repid).Mag();
    Double_t pdgmass[3] = {ProtonMass, PionMass, PionMass};
    event.GFeloss[j] = TMath::Sqrt(GFmom_decays.Mag()*GFmom_decays.Mag() + pdgmass[j]*pdgmass[j]) - TMath::Sqrt(GFTrackCont.GetMom(igf, 0, repid).Mag()*GFTrackCont.GetMom(igf, 0, repid).Mag() + pdgmass[j]*pdgmass[j]);

    event.GFlayer[j].resize(nh);
    event.GFrow[j].resize(nh);
    event.GFpos_x[j].resize(nh);
    event.GFpos_y[j].resize(nh);
    event.GFpos_z[j].resize(nh);
    event.GFmom[j].resize(nh);
    event.GFmom_x[j].resize(nh);
    event.GFmom_y[j].resize(nh);
    event.GFmom_z[j].resize(nh);
    event.GFresidual_x[j].resize(nh);
    event.GFresidual_y[j].resize(nh);
    event.GFresidual_z[j].resize(nh);
    event.GFresidual_p[j].resize(nh);
    event.GFresidual_px[j].resize(nh);
    event.GFresidual_py[j].resize(nh);
    event.GFresidual_pz[j].resize(nh);
    event.GFresolution_x[j].resize(nh);
    event.GFresolution_y[j].resize(nh);
    event.GFresolution_z[j].resize(nh);
    event.GFresolution_p[j].resize(nh);
    event.GFresolution_px[j].resize(nh);
    event.GFresolution_py[j].resize(nh);
    event.GFresolution_pz[j].resize(nh);
    event.GFpull_x[j].resize(nh);
    event.GFpull_y[j].resize(nh);
    event.GFpull_z[j].resize(nh);
    event.GFpull_p[j].resize(nh);
    event.GFpull_px[j].resize(nh);
    event.GFpull_py[j].resize(nh);
    event.GFpull_pz[j].resize(nh);

    Int_t id = GFxi_p_id_container[gfbest];
    if(j==1) id = GFxi_pi_id_container[gfbest];
    if(j==2) id = GFxi_pi2_id_container[gfbest];
    event.GFpmom = GFp_mom_container[gfbest].Mag();
    event.GFpmom_x = GFp_mom_container[gfbest].x();
    event.GFpmom_y = GFp_mom_container[gfbest].y();
    event.GFpmom_z = GFp_mom_container[gfbest].z();
    event.GFpi1mom = GFpi_mom_container[gfbest].Mag();
    event.GFpi1mom_x = GFpi_mom_container[gfbest].x();
    event.GFpi1mom_y = GFpi_mom_container[gfbest].y();
    event.GFpi1mom_z = GFpi_mom_container[gfbest].z();
    event.GFpi2mom = GFpi2_mom_container[gfbest].Mag();
    event.GFpi2mom_x = GFpi2_mom_container[gfbest].x();
    event.GFpi2mom_y = GFpi2_mom_container[gfbest].y();
    event.GFpi2mom_z = GFpi2_mom_container[gfbest].z();



    int ihit = 0;
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( id );

    double GFchisqrPos=0;
    for( Int_t ih=0; ih<tp->GetNHit(); ++ih ){
      Int_t layer = (Int_t) event.hitlayer[id][ih];
      TPCLTrackHit *helix_point = tp -> GetHitInOrder(ih);
      if(!helix_point -> IsGoodForTracking()) continue;
      const TVector3 &hit0 = helix_point -> GetLocalHitPos();
      TVector3 mom0 = helix_point -> GetMomentumHelix(event.GFcharge[j]);
      double row = helix_point->GetMRow()+0.5;
      TVector3 hit = GFTrackCont.GetPos(igf, ihit, repid);
      TVector3 mom = GFTrackCont.GetMom(igf, ihit, repid);
      double residual_[5];
      double pull_[5];
      double GFresidual6D[6];
      double GFresolution6D[6];
      double GFpull6D[6];
      TVector3 dumV;
      double dumd;
      GFTrackCont.GetTrackPull(igf,event.GFpdgcode[j],dumV,
             dumd,mom0,hit0,residual_,
             pull_,GFresidual6D,GFpull6D);

      event.GFlayer[j][ihit] = layer;
      event.GFrow[j][ihit] = row;
      event.GFmom_x[j][ihit] = mom.x();
      event.GFmom_y[j][ihit] = mom.y();
      event.GFmom_z[j][ihit] = mom.z();
      event.GFmom[j][ihit] = mom.Mag();
      event.GFpos_x[j][ihit] = hit.x();
      event.GFpos_y[j][ihit] = hit.y();
      event.GFpos_z[j][ihit] = hit.z();
      event.GFresidual_x[j][ihit] = hit.x() - hit0.x();
      event.GFresidual_y[j][ihit] = hit.y() - hit0.y();
      event.GFresidual_z[j][ihit] = hit.z() - hit0.z();
      event.GFresidual_p[j][ihit] = mom.Mag() - mom0.Mag();
      event.GFresidual_px[j][ihit] = mom.x() - mom0.x();
      event.GFresidual_py[j][ihit] = mom.y() - mom0.y();
      event.GFresidual_pz[j][ihit] = mom.z() - mom0.z();

      event.GFresolution_x[j][ihit] = GFresidual6D[0];
       event.GFresolution_y[j][ihit] = GFresidual6D[1];
       event.GFresolution_z[j][ihit] = GFresidual6D[2];
       event.GFresolution_px[j][ihit] = GFresidual6D[3];
       event.GFresolution_py[j][ihit] = GFresidual6D[4];
       event.GFresolution_pz[j][ihit] = GFresidual6D[5];
      event.GFpull_x[j][ihit] = GFpull6D[0];
       event.GFpull_y[j][ihit] = GFpull6D[1];
       event.GFpull_z[j][ihit] = GFpull6D[2];
       event.GFpull_px[j][ihit] = GFpull6D[3];
       event.GFpull_py[j][ihit] = GFpull6D[4];
       event.GFpull_pz[j][ihit] = GFpull6D[5];
      HF1(601,event.GFpull_x[j][ihit]);
      HF1(602,event.GFpull_y[j][ihit]);
      HF1(603,event.GFpull_z[j][ihit]);
      HF1(604,event.GFpull_px[j][ihit]);
      HF1(605,event.GFpull_py[j][ihit]);
      HF1(606,event.GFpull_pz[j][ihit]);
      if(event.GFpval[j]>0.01){
        HF1(611,event.GFpull_x[j][ihit]);
        HF1(612,event.GFpull_y[j][ihit]);
        HF1(613,event.GFpull_z[j][ihit]);
        HF1(614,event.GFpull_px[j][ihit]);
        HF1(615,event.GFpull_py[j][ihit]);
        HF1(616,event.GFpull_pz[j][ihit]);
      }
      ihit++;
    } //ih

    double GFndf = 2*ihit - 5; //Effective number of clusters
    double GFpvalPos = -1;
    if(GFndf > 0){
      GFchisqrPos /= GFndf;
      GFpvalPos = 1-ROOT::Math::chisquared_cdf(GFchisqrPos*GFndf, GFndf);
    }
    else GFchisqrPos = -1;

    event.GFchisqrPos[j]=GFchisqrPos;
    event.GFpvalPos[j]=GFpvalPos;
  } //igf
TVector3 XiMom = GFxi_mom_container[gfbest];
TVector3 XiVert = GFxi_vert_container[gfbest];
#if DoKinematicFitLdXi
  Double_t KFchisqrl = KFchisqrl_container[gfbest];
  Double_t KFpvall = KFpvall_container[gfbest];
  std::vector<Double_t> PullLd = KFlpull_container[gfbest];  
  Double_t KFchisqrxi;
  Double_t KFpvalxi;


  TVector3 KFTVP = KFp_mom_container[gfbest];
  TVector3 KFTVPi1 = KFpi_mom_container[gfbest];
  
  TVector3 KFHTVP = TVector3(KFTVP.X(),KFTVP.Z(),KFTVP.Y());
  TVector3 KFHTVPi1 = TVector3(KFTVPi1.X(),KFTVPi1.Z(),KFTVPi1.Y());  
  TVector3 HTVPi2(event.GFdecays_mom_x[2],event.GFdecays_mom_z[2],event.GFdecays_mom_y[2]); 
  TVector3 KFHTVLd = KFHTVP+KFHTVPi1;  
  TVector3 HTVXi = KFHTVLd+HTVPi2;


  TLorentzVector KFHLVLd(KFHTVLd,hypot(LambdaMass,KFHTVLd.Mag()));
  TLorentzVector HLVPi2(HTVPi2,hypot(PionMass,HTVPi2.Mag()));
  TLorentzVector HLVXi(HTVXi,hypot(XiMinusMass,HTVXi.Mag()));

  TMatrixD VLd = KFVarianceLd_container[gfbest]; 
  auto Vpi2 = Track_pi2->GetCovarianceMatrix();
  
  double VarianceXi[6] =
  {VLd(0,0),VLd(1,1),VLd(2,2),
  Vpi2(0,0),Vpi2(1,1),Vpi2(2,2)};

  auto OffdiagXi = MathTools::MergeOffdiagonals(VLd,Vpi2);  
  HLVXi = KFHLVLd+HLVPi2;
  FourVectorFitter KFXi(KFHLVLd,HLVPi2,HLVXi);
  KFXi.SetInvMass(XiMinusMass);
  KFXi.SetMaximumStep(5);
  KFXi.SetVariance(VarianceXi);
  KFXi.AddOffdiagonals(OffdiagXi);
  KFchisqrxi = KFXi.DoKinematicFit();
  cout<<Form("KFXi done:: chi2 = %g",KFchisqrxi)<<endl;
  auto VXi = KFXi.GetUnmeasuredCovariance();
  double rU,rV;
  MathTools::DecomposeResolutionUV(VXi,HTVXi,rU,rV);
  cout<<Form("Resolution U = %g, V = %g",rU,rV)<<endl;
  KFpvalxi = KFXi.GetPValue();  
  auto HcontXi = KFXi.GetFittedLV();
  auto PullXi = KFXi.GetPull();
  auto KFKFHLVLd = HcontXi.at(0);
  auto KFHLVPi2 = HcontXi.at(1);
  auto KFHLVXi = HcontXi.at(2);

  auto KFHLVP = TLorentzVector(KFHTVP.X(),KFHTVP.Z(),KFHTVP.Y(),hypot(ProtonMass,KFHTVP.Mag()));  
  auto KFHLVPi1 = TLorentzVector(KFHTVPi1.X(),KFHTVPi1.Z(),KFHTVPi1.Y(),hypot(PionMass,KFHTVPi1.Mag()));
  auto KFLVPi2 = TLorentzVector(KFHLVPi2.X(),KFHLVPi2.Z(),KFHLVPi2.Y(),KFHLVPi2.E());
  auto KFLVLd = TLorentzVector(KFKFHLVLd.X(),KFKFHLVLd.Z(),KFKFHLVLd.Y(),KFKFHLVLd.E());
  auto KFLVXi = TLorentzVector(KFHLVXi.X(),KFHLVXi.Z(),KFHLVXi.Y(),KFHLVXi.E());

  auto KFTVPi2 = KFLVPi2.Vect();
  auto KFTVLd = KFLVLd.Vect();
  auto KFTVXi = KFLVXi.Vect();

  event.KFchisqrxi = KFchisqrxi;
  event.KFpvalxi = KFpvalxi;
  event.KFximom = KFTVXi.Mag();
  event.KFximom_x = KFTVXi.x();
  event.KFximom_y = KFTVXi.y();
  event.KFximom_z = KFTVXi.z();

  event.KFchisqrl = KFchisqrl;
  event.KFpvall = KFpvall;
  event.KFlmom = KFTVLd.Mag();
  event.KFlmom_x = KFTVLd.x();
  event.KFlmom_y = KFTVLd.y();
  event.KFlmom_z = KFTVLd.z();

  event.KFdecays_mom.push_back(KFTVP.Mag());
  event.KFdecays_mom_x.push_back(KFTVP.x());
  event.KFdecays_mom_y.push_back(KFTVP.y());
  event.KFdecays_mom_z.push_back(KFTVP.z());
  event.KFdecays_mom.push_back(KFTVPi1.Mag());
  event.KFdecays_mom_x.push_back(KFTVPi1.x());
  event.KFdecays_mom_y.push_back(KFTVPi1.y());
  event.KFdecays_mom_z.push_back(KFTVPi1.z());
  event.KFdecays_mom.push_back(KFTVPi2.Mag());
  event.KFdecays_mom_x.push_back(KFTVPi2.x());
  event.KFdecays_mom_y.push_back(KFTVPi2.y());
  event.KFdecays_mom_z.push_back(KFTVPi2.z());


  HF1(10000,KFpvall);
  HF1(10001,KFchisqrl);
  HF1(10002,KFLVLd.Mag());

  if(KFchisqrxi>-1) XiMom = KFTVXi;
  for(int i=0;i<PullLd.size();++i){
    int num = 10010+i;
    HF1(num,PullLd.at(i));
  }
  TVector3 HTVP = TVector3(event.GFdecays_mom_x[0],event.GFdecays_mom_z[0],event.GFdecays_mom_y[0]);  
  TVector3 HTVPi1 = TVector3(event.GFdecays_mom_x[1],event.GFdecays_mom_z[1],event.GFdecays_mom_y[1]);
  TLorentzVector HLVP(HTVP,hypot(ProtonMass,HTVP.Mag())); 
  TLorentzVector HLVPi1(HTVPi1,hypot(PionMass,HTVPi1.Mag()));
  TLorentzVector HLVLd = HLVP+HLVPi1;
  HF1(10020,HLVP.P()-KFHLVP.P());
  HF1(10021,HLVP.Theta()-KFHLVP.Theta());
  HF1(10022,HLVP.Phi()-KFHLVP.Phi());
  HF1(10023,HLVPi1.P()-KFHLVPi1.P());
  HF1(10024,HLVPi1.Theta()-KFHLVPi1.Theta());
  HF1(10025,HLVPi1.Phi()-KFHLVPi1.Phi());
  HF1(10026,HLVLd.Phi()-KFHLVLd.Phi());
  HF1(10027,HLVLd.Phi()-KFHLVLd.Phi());
  HF1(10028,HLVLd.Phi()-KFHLVLd.Phi());

  HF1(20000,KFpvalxi);
  HF1(20001,KFchisqrxi);
  HF1(20002,KFLVXi.Mag());
  HF1(20020,HLVLd.P()-KFKFHLVLd.P());
  HF1(20021,HLVLd.Theta()-KFKFHLVLd.Theta());
  HF1(20022,HLVLd.Phi()-KFKFHLVLd.Phi());
  HF1(20023,HLVPi2.P()-KFHLVPi2.P());
  HF1(20024,HLVPi2.Theta()-KFHLVPi2.Theta());
  HF1(20025,HLVPi2.Phi()-KFHLVPi2.Phi());
  HF1(20026,HLVXi.Phi()-KFHLVXi.Phi());
  HF1(20027,HLVXi.Phi()-KFHLVXi.Phi());
  HF1(20028,HLVXi.Phi()-KFHLVXi.Phi());
  for(int i=0;i<PullXi.size();++i){
    int num = 20010+i;
    HF1(num,PullXi.at(i));
  }
  TVector3 G4HTVP(event.G4pmom_x,event.G4pmom_z,event.G4pmom_y);
  TVector3 G4HTVPi1(event.G4pi1mom_x,event.G4pi1mom_z,event.G4pi1mom_y);
  TVector3 G4HTVPi2(event.G4pi2mom_x,event.G4pi2mom_z,event.G4pi2mom_y);
  TVector3 G4HTVLd = G4HTVP + G4HTVPi1;
  TVector3 G4HTVXi = G4HTVLd + G4HTVPi2;

  auto KFHTVXi = KFHLVXi.Vect();
  auto KFKFHTVLd = KFKFHLVLd.Vect();
  auto KFHTVPi2 = KFHLVPi2.Vect();

  HF1(11010,(KFHTVP.Mag()-G4HTVP.Mag())/event.decays_res_mom.at(0));
  HF1(11011,(KFHTVP.Theta()-G4HTVP.Theta())/event.decays_res_th.at(0));
  HF1(11012,(KFHTVP.Phi()-G4HTVP.Phi())/event.decays_res_ph.at(0));
  HF1(11013,(KFHTVPi1.Mag()-G4HTVPi1.Mag())/event.decays_res_mom.at(1));
  HF1(11014,(KFHTVPi1.Theta()-G4HTVPi1.Theta())/event.decays_res_th.at(1));
  HF1(11015,(KFHTVPi1.Phi()-G4HTVPi1.Phi())/event.decays_res_ph.at(1));
  HF1(11016,(KFHTVLd.Mag()-G4HTVLd.Mag())/sqrt(VLd(0,0)));
  HF1(11017,(KFHTVLd.Theta()-G4HTVLd.Theta())/sqrt(VLd(1,1)));
  HF1(11018,(KFHTVLd.Phi()-G4HTVLd.Phi())/sqrt(VLd(2,2)));

  HF1(11020,KFHTVP.Mag()-G4HTVP.Mag());
  HF1(11021,KFHTVP.Theta()-G4HTVP.Theta());
  HF1(11022,KFHTVP.Phi()-G4HTVP.Phi());
  HF1(11023,KFHTVPi1.Mag()-G4HTVPi1.Mag());
  HF1(11024,KFHTVPi1.Theta()-G4HTVPi1.Theta());
  HF1(11025,KFHTVPi1.Phi()-G4HTVPi1.Phi());
  HF1(11026,KFHTVLd.Mag()-G4HTVLd.Mag());
  HF1(11027,KFHTVLd.Theta()-G4HTVLd.Theta());
  HF1(11028,KFHTVLd.Phi()-G4HTVLd.Phi());

  HF1(21010,(KFKFHTVLd.Mag()-G4HTVLd.Mag())/sqrt(VLd(0,0)));
  HF1(21011,(KFKFHTVLd.Theta()-G4HTVLd.Theta())/sqrt(VLd(1,1)));
  HF1(21012,(KFKFHTVLd.Phi()-G4HTVLd.Phi())/sqrt(VLd(2,2)));
  HF1(21013,(KFHTVPi2.Mag()-G4HTVPi2.Mag())/event.decays_res_mom.at(2));
  HF1(21014,(KFHTVPi2.Theta()-G4HTVPi2.Theta())/event.decays_res_th.at(2));
  HF1(21015,(KFHTVPi2.Phi()-G4HTVPi2.Phi())/event.decays_res_ph.at(2));

  HF1(21020,KFKFHTVLd.Mag()-G4HTVLd.Mag());
  HF1(21021,KFKFHTVLd.Theta()-G4HTVLd.Theta());
  HF1(21022,KFKFHTVLd.Phi()-G4HTVLd.Phi());
  HF1(21023,KFHTVPi2.Mag()-G4HTVPi2.Mag());
  HF1(21024,KFHTVPi2.Theta()-G4HTVPi2.Theta());
  HF1(21025,KFHTVPi2.Phi()-G4HTVPi2.Phi());
  HF1(21026,KFHTVXi.Mag()-G4HTVXi.Mag());
  HF1(21027,KFHTVXi.Theta()-G4HTVXi.Theta());
  HF1(21028,KFHTVXi.Phi()-G4HTVXi.Phi());

#endif
  GFTrackCont.AddReconstructedTrack(XiMinusPdgCode,XiVert,XiMom);
  GFTrackCont.FitTrack(GFTrackCont.GetNTrack()-1);
  TVector3 XiTgtVert, XiTgtMom;
  double XiTgtLen,XiTgtTof;
  bool XiFlight = GFTrackCont.ExtrapolateToTargetCenter(GFTrackCont.GetNTrack()-1
  ,XiTgtVert,XiTgtMom,XiTgtLen,XiTgtTof);
  if(XiFlight){
    event.xtgtXi = XiTgtVert.x();
    event.ytgtXi = XiTgtVert.y();
    event.utgtXi = XiTgtMom.x()/XiTgtMom.Z();
    event.vtgtXi = XiTgtMom.y()/XiTgtMom.Z();
  }
  const int ntrack = 3;
  TVector3 G4Km(event.G4kmmom_x,event.G4kmmom_z,event.G4kmmom_y);
  TVector3 G4Kp(event.G4kpmom_x,event.G4kpmom_z,event.G4kpmom_y);
  double uKm = G4Km.x()/G4Km.z();
  double vKm = G4Km.y()/G4Km.z();
  double uKp = G4Kp.x()/G4Kp.z();
  double vKp = G4Kp.y()/G4Kp.z();
  double x0track[ntrack]={event.xtgtHS.at(0),event.xtgtKurama.at(0),event.xtgtXi};
  double y0track[ntrack]={event.ytgtHS.at(0),event.ytgtKurama.at(0),event.ytgtXi};
  double u0track[ntrack]={uKm,uKp,event.utgtXi};
  double v0track[ntrack]={vKm,vKp,event.vtgtXi};
  TVector3 KKXiVert = Kinematics::MultitrackVertex(ntrack,x0track,y0track,u0track,v0track);
  event.XiFlight = XiFlight;
  if(XiFlight){
    event.vtxKKXi= KKXiVert.x();
    event.vtyKKXi= KKXiVert.y();
    event.vtzKKXi= KKXiVert.z();
  }
  TVector3 XiProdVert,XiProdMom;
  double XiProdLen,XiProdTof;
  bool XiProd = false;
  if(XiFlight){
    XiProd = GFTrackCont.XiDecayToProdVertex(GFTrackCont.GetNTrack()-1
  ,KKXiVert,XiProdVert,XiProdMom,XiProdLen,XiProdTof);
  }
  event.XiProd = XiProd;
  TVector3 G4XiProdMom(event.G4ximom_x,event.G4ximom_z,event.G4ximom_y);
  if(XiProd){
    event.xiprodvtx_x = XiProdVert.x();
    event.xiprodvtx_y = XiProdVert.y();
    event.xiprodvtx_z = XiProdVert.z();
    event.xiprodmom_x = XiProdMom.x();
    event.xiprodmom_y = XiProdMom.y();
    event.xiprodmom_z = XiProdMom.z();
    HF1(3000,XiMom.Mag()-XiProdMom.Mag());
    HF1(3001,G4XiProdMom.Mag()-XiProdMom.Mag());
  }


  auto HTVLd = HLVP+HLVPi1;
  HF1(12010,(HTVP.Mag()-G4HTVP.Mag())/event.decays_res_mom.at(0));
  HF1(12011,(HTVP.Theta()-G4HTVP.Theta())/event.decays_res_th.at(0));
  HF1(12012,(HTVP.Phi()-G4HTVP.Phi())/event.decays_res_ph.at(0));
  HF1(12013,(HTVPi1.Mag()-G4HTVPi1.Mag())/event.decays_res_mom.at(1));
  HF1(12014,(HTVPi1.Theta()-G4HTVPi1.Theta())/event.decays_res_th.at(1));
  HF1(12015,(HTVPi1.Phi()-G4HTVPi1.Phi())/event.decays_res_ph.at(1));
  HF1(12020,HTVP.Mag()-G4HTVP.Mag());
  HF1(12021,HTVP.Theta()-G4HTVP.Theta());
  HF1(12022,HTVP.Phi()-G4HTVP.Phi());
  HF1(12023,HTVPi1.Mag()-G4HTVPi1.Mag());
  HF1(12024,HTVPi1.Theta()-G4HTVPi1.Theta());
  HF1(12025,HTVPi1.Phi()-G4HTVPi1.Phi());
  HF1(12026,HTVLd.Mag()-G4HTVLd.Mag());
  HF1(12027,HTVLd.Theta()-G4HTVLd.Theta());
  HF1(12028,HTVLd.Phi()-G4HTVLd.Phi());

  HF1(22013,(HTVPi2.Mag()-G4HTVPi2.Mag())/event.decays_res_mom.at(2));
  HF1(22014,(HTVPi2.Theta()-G4HTVPi2.Theta())/event.decays_res_th.at(2));
  HF1(22015,(HTVPi2.Phi()-G4HTVPi2.Phi())/event.decays_res_ph.at(2));
  HF1(22020,HTVLd.Mag()-G4HTVLd.Mag());
  HF1(22021,HTVLd.Theta()-G4HTVLd.Theta());
  HF1(22022,HTVLd.Phi()-G4HTVLd.Phi());
  HF1(22023,HTVPi2.Mag()-G4HTVPi2.Mag());
  HF1(22024,HTVPi2.Theta()-G4HTVPi2.Theta());
  HF1(22025,HTVPi2.Phi()-G4HTVPi2.Phi());
  HF1(22026,HTVXi.Mag()-G4HTVXi.Mag());
  HF1(22027,HTVXi.Theta()-G4HTVXi.Theta());
  HF1(22028,HTVXi.Phi()-G4HTVXi.Phi());

  HF1( 31, event.GFlmass);
  HF1( 32, event.GFximass);
  HF1( 33, event.GFdecays_mom[0]);
  HF1( 34, event.GFdecays_mom[1]);
  HF1( 35, event.GFdecays_mom[2]);
  HF1( 36, event.GFlmom);
  HF1( 37, event.GFximom);
  HF1( 38, event.GFppi_dist);
  HF1( 39, event.GFlpi_dist);
  HF1( 40, event.GFdecays_mass[0]);
  HF1( 41, event.GFdecays_mass[1]);
  HF1( 42, event.GFdecays_mass[2]);
  HF1( 2, event.GFstatus++);
  HF1( 1, event.status++ );

  GFTrackCont.Clear();
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

  HB1( 1, "Status", 21, 0., 21. );
  HB1( 2, "Genfit Status", 20, 0., 20. );
  HB1( 3, "Genfit Fit Status", 2, 0., 2. );
  HB1( 11," #Lambda Invariant Mass; M_{p#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 80, 1.04, 1.2);
  HB1( 12," #Xi Invariant Mass; M_{#Lambda#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 125, 1.2 ,1.45);
  HB1( 13, "p Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 14, "#pi Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  HB1( 15, "#pi Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  HB1( 16, "#Lambda Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 17, "#Xi Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 18, "Closest Dist. p#pi;Distance [mm]; Counts [/0.1 mm]", 1000, 0., 100.);
  HB1( 19, "Closest Dist. #Lambda#pi;Distance [mm]; Counts [/0.1 mm]", 1000, 0., 100.);

  HB1( 21," #Lambda Invariant Mass; M_{p#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 80, 1.04, 1.2);
  HB1( 22," #Xi Invariant Mass; M_{#Lambda#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 125, 1.2 ,1.45);
  HB1( 23, "p Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 24, "#pi Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  HB1( 25, "#pi Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  HB1( 26, "#Lambda Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 27, "#Xi Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 28, "Closest Dist. p#pi;Distance [mm]; Counts [/0.1 mm]", 1000, 0., 100.);
  HB1( 29, "Closest Dist. #Lambda#pi;Distance [mm]; Counts [/0.1 mm]", 1000, 0., 100.);








  HB1( 30, "[GenFit] Missing Mass [KURAMA]; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 160, 1., 1.8);
  HB1( 31, "[GenFit] #Lambda Invariant Mass; M_{p#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 80, 1.04, 1.2);
  HB1( 32, "[GenFit] #Xi Invariant Mass; M_{#Lambda#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 125, 1.2 ,1.45);
  HB1( 33, "[GenFit] p_{#Lambda} Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 34, "[GenFit] #pi_{#Lambda} Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  HB1( 35, "[GenFit] #pi_{#Xi} Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  HB1( 36, "[GenFit] #Lambda Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 37, "[GenFit] #Xi Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 38, "[GenFit] Closest Dist. p#pi;Distance [mm]; Counts [/0.1 mm]", 1000, 0., 100.);
  HB1( 39, "[GenFit] Closest Dist. #Lambda#pi;Distance [mm]; Counts [/0.1 mm]", 1000, 0., 100.);
  HB1( 40, "[GenFit] p_{#Lambda} Mass; M [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 750, 0., 1.5);
  HB1( 41, "[GenFit] pi_{#Lambda} Mass; M [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 750, 0., 1.5);
  HB1( 42, "[GenFit] pi_{#Xi} Mass; M [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 750, 0., 1.5);

  HB1( 101, "p_{x} Residuals",1000,-0.5,0.5);
  HB1( 102, "p_{x} Resolutions",1000,0.,0.1);
  HB1( 103, "p_{x} Pull",1000,-5,5);
  HB1( 201, "p_{y} Residuals",1000,-0.5,0.5);
  HB1( 202, "p_{y} Resolutions",1000,0.,0.1);
  HB1( 203, "p_{y} Pull",1000,-5,5);
  HB1( 301, "p_{z} Residuals",1000,-0.5,0.5);
  HB1( 302, "p_{z} Resolutions",1000,0.,0.1);
  HB1( 303, "p_{z} Pull",1000,-5,5);
  HB1( 400, "p_{T} Resolutions",10000,0,0.1);
  HB1( 401, "#phi Resolutions",10000,0,0.1);
  HB1( 402, "dz Resolutions",10000,0,0.01);

  HB1( 500, "pos Resolutions",10000,0,0.5);
  HB1( 501, "pos_{x} Resolutions",10000,0,0.5);
  HB1( 502, "pos_{y} Resolutions",10000,0,0.5);
  HB1( 503, "pos_{z} Resolutions",10000,0,0.5);
  HB1(600,"GF p-value",1000,0,1);
  HB1(601,"GFx Pull",1000,-5,5);
  HB1(602,"GFy Pull",1000,-5,5);
  HB1(603,"GFz Pull",1000,-5,5);
  HB1(604,"GFp_{x} Pull",1000,-5,5);
  HB1(605,"GFp_{y} Pull",1000,-5,5);
  HB1(606,"GFp_{z} Pull",1000,-5,5);
  HB1(611,"GFx Pull pval>0.01",1000,-5,5);
  HB1(612,"GFy Pull pval>0.01",1000,-5,5);
  HB1(613,"GFz Pull pval>0.01",1000,-5,5);
  HB1(614,"GFp_{x} Pull pval>0.01",1000,-5,5);
  HB1(615,"GFp_{y} Pull pval>0.01",1000,-5,5);
  HB1(616,"GFp_{z} Pull pval>0.01",1000,-5,5);


  HB2(1001, "p hit pattern",100,-250,250,100,-250,250);
  HB2(1002, "#pi_{#Lambda} hit pattern",100,-250,250,100,-250,250);
  HB2(1011, "p hit pattern_{WO#Xi}",100,-250,250,100,-250,250);
  HB2(1012, "#pi_{#Lambda} hit pattern_{WO#Xi}",100,-250,250,100,-250,250);
  HB2(1003, "#pi_{#Xi} hit pattern",100,-250,250,100,-250,250);
  HB2(2001, "p hit patternG4",100,-250,250,100,-250,250);
  HB2(2002, "#pi_{#Lambda} hit patternG4",100,-250,250,100,-250,250);
  HB2(2003, "#pi_{#Xi} hit patternG4",100,-250,250,100,-250,250);

  HB1(3000, "#Xi Decay mom - Prod mom; #Delta p [GeV/#font[12]{c}]; Counts [/ 2MeV/#font[12]{c}]", 300, -0.3, 0.3);
  HB1(3001, "#Xi G4 Prod mom - Prod mom; #Delta p [GeV/#font[12]{c}]; Counts [/ 2MeV/#font[12]{c}]", 300, -0.3, 0.3);
#if DoKinematicFitLdXi
  HB1(10000,"KF#{Lambda} pvalue",100,0,1);
  HB1(10001,"KF#{Lambda} chisqr",1000,0,15);
  HB1(10002,"KF#{Lambda} mass",1000,pdg::LambdaMass()-0.1,pdg::LambdaMass()+0.1);
  HB1(10010,"KF#{Lambda} pull p_{p}",100,-5,5);
  HB1(10011,"KF#{Lambda} pull p_{#theta}",100,-5,5);
  HB1(10012,"KF#{Lambda} pull p_{#phi}",100,-5,5);
  HB1(10013,"KF#{Lambda} pull #pi_{p}",100,-5,5);
  HB1(10014,"KF#{Lambda} pull #pi_{#theta}",100,-5,5);
  HB1(10015,"KF#{Lambda} pull #pi_{#phi}",100,-5,5);
  HB1(10016,"KF#{Lambda} pull #Lambda_{p}",100,-5,5);
  HB1(10017,"KF#{Lambda} pull #Lambda_{#theta}",100,-5,5);
  HB1(10018,"KF#{Lambda} pull #Lambda_{#phi}",100,-5,5);

  HB1(10020,"KF#{Lambda} residual p_{p}",1000,-1,1);
  HB1(10021,"KF#{Lambda} residual p_{#theta}",1000,-1,1);
  HB1(10022,"KF#{Lambda} residual p_{#phi}",1000,-1,1);
  HB1(10023,"KF#{Lambda} residual #pi_{p}",1000,-1,1);
  HB1(10024,"KF#{Lambda} residual #pi_{#theta}",1000,-1,1);
  HB1(10025,"KF#{Lambda} residual #pi_{#phi}",1000,-1,1);


  HB1(20000,"KF#{Xi}^{-} pvalue",100,0,1);
  HB1(20001,"KF#{Xi}^{-} chisqr",1000,0,15);
  HB1(20002,"KF#{Xi} mass",1000,pdg::XiMinusMass()-0.1,pdg::XiMinusMass()+0.1);
  HB1(20010,"KF#{Xi}^{-} pull #Lambda_{p}",100,-5,5);
  HB1(20011,"KF#{Xi}^{-} pull #Lambda_{#theta}",100,-5,5);
  HB1(20012,"KF#{Xi}^{-} pull #Lambda_{#phi}",100,-5,5);
  HB1(20013,"KF#{Xi}^{-} pull #pi_{p}",100,-5,5);
  HB1(20014,"KF#{Xi}^{-} pull #pi_{#theta}",100,-5,5);
  HB1(20015,"KF#{Xi}^{-} pull #pi_{#phi}",100,-5,5);

  HB1(20020,"KF#{Xi}^{-} residual #Lambda_{p}",1000,-1,1);
  HB1(20021,"KF#{Xi}^{-} residual #Lambda_{#theta}",1000,-0.5,0.5);
  HB1(20022,"KF#{Xi}^{-} residual #Lambda_{#phi}",1000,-0.5,0.5);
  HB1(20023,"KF#{Xi}^{-} residual #pi_{p}",1000,-0.3,0.3);
  HB1(20024,"KF#{Xi}^{-} residual #pi_{#theta}",1000,-1,1);
  HB1(20025,"KF#{Xi}^{-} residual #pi_{#phi}",1000,-1,1);
  HB1(20026,"KF#{Xi}^{-} residual #Xi_{p}",1000,-0.3,0.3);
  HB1(20027,"KF#{Xi}^{-} residual #Xi_{#theta}",1000,-0.1,0.1);
  HB1(20028,"KF#{Xi}^{-} residual #Xi_{#phi}",1000,-0.3,0.3);

  HB1(11010,"KF#{Lambda} G4pull p_{p} ",100,-5,5);
  HB1(11011,"KF#{Lambda} G4pull p_{#theta} ",100,-5,5);
  HB1(11012,"KF#{Lambda} G4pull p_{#phi} ",100,-5,5);
  HB1(11013,"KF#{Lambda} G4pull #pi_{p} ",100,-5,5);
  HB1(11014,"KF#{Lambda} G4pull #pi_{#theta} ",100,-5,5);
  HB1(11015,"KF#{Lambda} G4pull #pi_{#phi} ",100,-5,5);
  HB1(11016,"KF#{Lambda} G4pull #Lambda_{p} ",100,-5,5);
  HB1(11017,"KF#{Lambda} G4pull #Lambda_{#theta} ",100,-5,5);
  HB1(11018,"KF#{Lambda} G4pull #Lambda_{#phi} ",100,-5,5);
  HB1(11020,"KF#{Lambda} G4residual p_{p}",1000,-1,1);
  HB1(11021,"KF#{Lambda} G4residual p_{#theta}",1000,-1,1);
  HB1(11022,"KF#{Lambda} G4residual p_{#phi}",1000,-1,1);
  HB1(11023,"KF#{Lambda} G4residual #pi_{p}",1000,-1,1);
  HB1(11024,"KF#{Lambda} G4residual #pi_{#theta}",1000,-1,1);
  HB1(11025,"KF#{Lambda} G4residual #pi_{#phi}",1000,-1,1);
  HB1(11026,"KF#{Lambda} G4residual #Lambda_{p}",1000,-1,1);
  HB1(11027,"KF#{Lambda} G4residual #Lambda_{#theta}",1000,-1,1);
  HB1(11028,"KF#{Lambda} G4residual #Lambda_{#phi}",1000,-1,1);

  HB1(21010,"KF#{Xi}^{-} G4pull #Lambda_{p} ",100,-5,5);
  HB1(21011,"KF#{Xi}^{-} G4pull #Lambda_{#theta} ",100,-5,5);
  HB1(21012,"KF#{Xi}^{-} G4pull #Lambda_{#phi} ",100,-5,5);
  HB1(21013,"KF#{Xi}^{-} G4pull #pi_{p} ",100,-5,5);
  HB1(21014,"KF#{Xi}^{-} G4pull #pi_{#theta} ",100,-5,5);
  HB1(21015,"KF#{Xi}^{-} G4pull #pi_{#phi} ",100,-5,5);
  HB1(21020,"KF#{Xi}^{-} G4residual #Lambda_{p}",1000,-1,1);
  HB1(21021,"KF#{Xi}^{-} G4residual #Lambda_{#theta}",1000,-0.5,0.5);
  HB1(21022,"KF#{Xi}^{-} G4residual #Lambda_{#phi}",1000,-0.5,0.5);
  HB1(21023,"KF#{Xi}^{-} G4residual #pi_{p}",1000,-0.3,0.3);
  HB1(21024,"KF#{Xi}^{-} G4residual #pi_{#theta}",1000,-1,1);
  HB1(21025,"KF#{Xi}^{-} G4residual #pi_{#phi}",1000,-1,1);
  HB1(21026,"KF#{Xi}^{-} G4residual #Xi_{p}",1000,-0.3,0.3);
  HB1(21027,"KF#{Xi}^{-} G4residual #Xi_{#theta}",1000,-0.1,0.1);
  HB1(21028,"KF#{Xi}^{-} G4residual #Xi_{#phi}",1000,-0.3,0.3);


#endif
  HB1(12010,"#{Lambda} G4pull p_{p} ",100,-5,5);
  HB1(12011,"#{Lambda} G4pull p_{#theta} ",100,-5,5);
  HB1(12012,"#{Lambda} G4pull p_{#phi} ",100,-5,5);
  HB1(12013,"#{Lambda} G4pull #pi_{p} ",100,-5,5);
  HB1(12014,"#{Lambda} G4pull #pi_{#theta} ",100,-5,5);
  HB1(12015,"#{Lambda} G4pull #pi_{#phi} ",100,-5,5);
  HB1(12016,"#{Lambda} G4pull #Lambda_{p} ",100,-5,5);
  HB1(12017,"#{Lambda} G4pull #Lambda_{#theta} ",100,-5,5);
  HB1(12018,"#{Lambda} G4pull #Lambda_{#phi} ",100,-5,5);
  HB1(12020,"#{Lambda} G4residual p_{p}",1000,-1,1);
  HB1(12021,"#{Lambda} G4residual p_{#theta}",1000,-1,1);
  HB1(12022,"#{Lambda} G4residual p_{#phi}",1000,-1,1);
  HB1(12023,"#{Lambda} G4residual #pi_{p}",1000,-1,1);
  HB1(12024,"#{Lambda} G4residual #pi_{#theta}",1000,-1,1);
  HB1(12025,"#{Lambda} G4residual #pi_{#phi}",1000,-1,1);
  HB1(12026,"#{Lambda} G4residual #Lambda_{p}",1000,-1,1);
  HB1(12027,"#{Lambda} G4residual #Lambda_{#theta}",1000,-1,1);
  HB1(12028,"#{Lambda} G4residual #Lambda_{#phi}",1000,-1,1);

  HB1(22010,"#{Xi}^{-} G4pull #Lambda_{p} ",100,-5,5);
  HB1(22011,"#{Xi}^{-} G4pull #Lambda_{#theta} ",100,-5,5);
  HB1(22012,"#{Xi}^{-} G4pull #Lambda_{#phi} ",100,-5,5);
  HB1(22013,"#{Xi}^{-} G4pull #pi_{p} ",100,-5,5);
  HB1(22014,"#{Xi}^{-} G4pull #pi_{#theta} ",100,-5,5);
  HB1(22015,"#{Xi}^{-} G4pull #pi_{#phi} ",100,-5,5);
  HB1(22020,"#{Xi}^{-} G4residual #Lambda_{p}",1000,-1,1);
  HB1(22021,"#{Xi}^{-} G4residual #Lambda_{#theta}",1000,-0.5,0.5);
  HB1(22022,"#{Xi}^{-} G4residual #Lambda_{#phi}",1000,-0.5,0.5);
  HB1(22023,"#{Xi}^{-} G4residual #pi_{p}",1000,-0.3,0.3);
  HB1(22024,"#{Xi}^{-} G4residual #pi_{#theta}",1000,-1,1);
  HB1(22025,"#{Xi}^{-} G4residual #pi_{#phi}",1000,-1,1);
  HB1(22026,"#{Xi}^{-} G4residual #Xi_{p}",1000,-0.3,0.3);
  HB1(22027,"#{Xi}^{-} G4residual #Xi_{#theta}",1000,-0.1,0.1);
  HB1(22028,"#{Xi}^{-} G4residual #Xi_{#phi}",1000,-0.3,0.3);



  HBTree( "tpc", "tree of DstTPCTracking" );
  tree->Branch( "status", &event.status );
  tree->Branch( "evnum", &event.evnum );

  tree->Branch( "nhHtof", &event.nhHtof );
  tree->Branch( "HtofSeg", &event.HtofSeg );
  tree->Branch( "tHtof", &event.tHtof );
  tree->Branch( "deHtof", &event.deHtof );
  tree->Branch( "posHtof", &event.posHtof );

  tree->Branch( "nhFtof", &event.nhFtof );
  tree->Branch( "FtofSeg", &event.FtofSeg );
  tree->Branch( "tFtof", &event.tFtof );
  tree->Branch( "deFtof", &event.deFtof );
  tree->Branch( "posFtof", &event.posFtof );

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "trackid", &event.trackid );
  tree->Branch( "isBeam", &event.isBeam );
  tree->Branch( "isKurama", &event.isKurama );
  tree->Branch( "isK18", &event.isK18 );
  tree->Branch( "isAccidental", &event.isAccidental );
  tree->Branch( "charge", &event.charge );
  tree->Branch( "pid", &event.pid );
  tree->Branch( "purity", &event.purity );
  tree->Branch( "efficiency", &event.efficiency );
  tree->Branch(  "G4tid", &event.G4tid );
  tree->Branch( "chisqr", &event.chisqr );
  tree->Branch( "pval", &event.pval );
  tree->Branch( "helix_cx", &event.helix_cx );
  tree->Branch( "helix_cy", &event.helix_cy );
  tree->Branch( "helix_z0", &event.helix_z0 );
  tree->Branch( "helix_r", &event.helix_r );
  tree->Branch( "helix_dz", &event.helix_dz );
  tree->Branch( "dE", &event.dE );
  tree->Branch( "dEdx", &event.dEdx );
  tree->Branch( "mom0", &event.mom0 );
  tree->Branch( "path", &event.path );
  tree->Branch( "helix_t", &event.helix_t );
  tree->Branch( "resolution_x", &event.resolution_x);
  tree->Branch( "resolution_y", &event.resolution_y);
  tree->Branch( "resolution_z", &event.resolution_z);

  tree->Branch("G4kmid",&event.G4kmid);
  tree->Branch("G4kmtid",&event.G4kmtid);
  tree->Branch("G4kmvtx_x",&event.G4kmvtx_x);
  tree->Branch("G4kmvtx_y",&event.G4kmvtx_y);
  tree->Branch("G4kmvtx_z",&event.G4kmvtx_z);
  tree->Branch("G4kmmom",&event.G4kmmom);
  tree->Branch("G4kmmom_x",&event.G4kmmom_x);
  tree->Branch("G4kmmom_y",&event.G4kmmom_y);
  tree->Branch("G4kmmom_z",&event.G4kmmom_z);

  tree->Branch("G4kpid",&event.G4kpid);
  tree->Branch("G4kptid",&event.G4kptid);
  tree->Branch("G4kpvtx_x",&event.G4kpvtx_x);
  tree->Branch("G4kpvtx_y",&event.G4kpvtx_y);
  tree->Branch("G4kpvtx_z",&event.G4kpvtx_z);
  tree->Branch("G4kpmom",&event.G4kpmom);
  tree->Branch("G4kpmom_x",&event.G4kpmom_x);
  tree->Branch("G4kpmom_y",&event.G4kpmom_y);
  tree->Branch("G4kpmom_z",&event.G4kpmom_z);


  tree->Branch("G4xiid",&event.G4xiid);
  tree->Branch("G4xivtx_x",&event.G4xivtx_x);
  tree->Branch("G4xivtx_y",&event.G4xivtx_y);
  tree->Branch("G4xivtx_z",&event.G4xivtx_z);
  tree->Branch("G4ximom",&event.G4ximom);
  tree->Branch("G4ximom_x",&event.G4ximom_x);
  tree->Branch("G4ximom_y",&event.G4ximom_y);
  tree->Branch("G4ximom_z",&event.G4ximom_z);

  tree->Branch("G4lid",&event.G4lid);
  tree->Branch("G4lvtx_x",&event.G4lvtx_x);
  tree->Branch("G4lvtx_y",&event.G4lvtx_y);
  tree->Branch("G4lvtx_z",&event.G4lvtx_z);
  tree->Branch("G4lmom",&event.G4lmom);
  tree->Branch("G4lmom_x",&event.G4lmom_x);
  tree->Branch("G4lmom_y",&event.G4lmom_y);
  tree->Branch("G4lmom_z",&event.G4lmom_z);


  tree->Branch("G4pid",&event.G4pid);
  tree->Branch("G4ptid",&event.G4ptid);
  tree->Branch("G4pnh",&event.G4pnh);
  tree->Branch("G4ptnh",&event.G4ptnh);
  tree->Branch("G4pvtx_x",&event.G4pvtx_x);
  tree->Branch("G4pvtx_y",&event.G4pvtx_y);
  tree->Branch("G4pvtx_z",&event.G4pvtx_z);
  tree->Branch("G4pmom",&event.G4pmom);
  tree->Branch("G4pmom_x",&event.G4pmom_x);
  tree->Branch("G4pmom_y",&event.G4pmom_y);
  tree->Branch("G4pmom_z",&event.G4pmom_z);
  tree->Branch("ptid",&event.ptid);
  tree->Branch("pnh",&event.pnh);
  tree->Branch("pvtx_x",&event.pvtx_x);
  tree->Branch("pvtx_y",&event.pvtx_y);
  tree->Branch("pvtx_z",&event.pvtx_z);
  tree->Branch("pmom",&event.pmom);
  tree->Branch("pmom_x",&event.pmom_x);
  tree->Branch("pmom_y",&event.pmom_y);
  tree->Branch("pmom_z",&event.pmom_z);
  tree->Branch("pres_x",&event.pres_x);
  tree->Branch("pres_y",&event.pres_y);
  tree->Branch("pres_z",&event.pres_z);
  tree->Branch("pcov_xy",&event.pcov_xy);
  tree->Branch("pcov_yz",&event.pcov_yz);
  tree->Branch("pcov_zx",&event.pcov_zx);
  tree->Branch("GFpmom",&event.GFpmom);
  tree->Branch("GFpmom_x",&event.GFpmom_x);
  tree->Branch("GFpmom_y",&event.GFpmom_y);
  tree->Branch("GFpmom_z",&event.GFpmom_z);

  tree->Branch("G4pi1id",&event.G4pi1id);
  tree->Branch("G4pi1tid",&event.G4pi1tid);
  tree->Branch("G4pi1nh",&event.G4pi1nh);
  tree->Branch("G4pi1tnh",&event.G4pi1tnh);
  tree->Branch("G4pi1vtx_x",&event.G4pi1vtx_x);
  tree->Branch("G4pi1vtx_y",&event.G4pi1vtx_y);
  tree->Branch("G4pi1vtx_z",&event.G4pi1vtx_z);
  tree->Branch("G4pi1mom",&event.G4pi1mom);
  tree->Branch("G4pi1mom_x",&event.G4pi1mom_x);
  tree->Branch("G4pi1mom_y",&event.G4pi1mom_y);
  tree->Branch("G4pi1mom_z",&event.G4pi1mom_z);
  tree->Branch("pi1tid",&event.pi1tid);
  tree->Branch("pi1nh",&event.pi1nh);
  tree->Branch("pi1vtx_x",&event.pi1vtx_x);
  tree->Branch("pi1vtx_y",&event.pi1vtx_y);
  tree->Branch("pi1vtx_z",&event.pi1vtx_z);
  tree->Branch("pi1mom",&event.pi1mom);
  tree->Branch("pi1mom_x",&event.pi1mom_x);
  tree->Branch("pi1mom_y",&event.pi1mom_y);
  tree->Branch("pi1mom_z",&event.pi1mom_z);
  tree->Branch("pi1res_x",&event.pi1res_x);
  tree->Branch("pi1res_y",&event.pi1res_y);
  tree->Branch("pi1res_z",&event.pi1res_z);
  tree->Branch("pi1cov_xy",&event.pi1cov_xy);
  tree->Branch("pi1cov_yz",&event.pi1cov_yz);
  tree->Branch("pi1cov_zx",&event.pi1cov_zx);
  tree->Branch("GFpi1mom",&event.GFpi1mom);
  tree->Branch("GFpi1mom_x",&event.GFpi1mom_x);
  tree->Branch("GFpi1mom_y",&event.GFpi1mom_y);
  tree->Branch("GFpi1mom_z",&event.GFpi1mom_z);


  tree->Branch("G4pi2id",&event.G4pi2id);
  tree->Branch("G4pi2tid",&event.G4pi2tid);
  tree->Branch("G4pi2nh",&event.G4pi2nh);
  tree->Branch("G4pi2tnh",&event.G4pi2tnh);
  tree->Branch("G4pi2vtx_x",&event.G4pi2vtx_x);
  tree->Branch("G4pi2vtx_y",&event.G4pi2vtx_y);
  tree->Branch("G4pi2vtx_z",&event.G4pi2vtx_z);
  tree->Branch("G4pi2mom",&event.G4pi2mom);
  tree->Branch("G4pi2mom_x",&event.G4pi2mom_x);
  tree->Branch("G4pi2mom_y",&event.G4pi2mom_y);
  tree->Branch("G4pi2mom_z",&event.G4pi2mom_z);
  tree->Branch("pi2tid",&event.pi2tid);
  tree->Branch("pi2nh",&event.pi2nh);
  tree->Branch("pi2vtx_x",&event.pi2vtx_x);
  tree->Branch("pi2vtx_y",&event.pi2vtx_y);
  tree->Branch("pi2vtx_z",&event.pi2vtx_z);
  tree->Branch("pi2mom",&event.pi2mom);
  tree->Branch("pi2mom_x",&event.pi2mom_x);
  tree->Branch("pi2mom_y",&event.pi2mom_y);
  tree->Branch("pi2mom_z",&event.pi2mom_z);
  tree->Branch("pi2res_x",&event.pi2res_x);
  tree->Branch("pi2res_y",&event.pi2res_y);
  tree->Branch("pi2res_z",&event.pi2res_z);
  tree->Branch("pi2cov_xy",&event.pi2cov_xy);
  tree->Branch("pi2cov_yz",&event.pi2cov_yz);
  tree->Branch("pi2cov_zx",&event.pi2cov_zx);
  tree->Branch("GFpi2mom",&event.GFpi2mom);
  tree->Branch("GFpi2mom_x",&event.GFpi2mom_x);
  tree->Branch("GFpi2mom_y",&event.GFpi2mom_y);
  tree->Branch("GFpi2mom_z",&event.GFpi2mom_z);




  tree->Branch("Lflag", &event.lflag);
  tree->Branch("Xiflag", &event.xiflag);
  tree->Branch("best", &event.best);
  tree->Branch("XiMass", &event.ximass);
  tree->Branch("XiDecayVtx_x", &event.xidecayvtx_x);
  tree->Branch("XiDecayVtx_y", &event.xidecayvtx_y);
  tree->Branch("XiDecayVtx_z", &event.xidecayvtx_z);
  tree->Branch("XiMom_x", &event.ximom_x);
  tree->Branch("XiMom_y", &event.ximom_y);
  tree->Branch("XiMom_z", &event.ximom_z);
  tree->Branch("XiVtxCloseDist", &event.lpi_dist);
  tree->Branch("LambdaMass", &event.lmass);
  tree->Branch("LambdaDecayVtx_x", &event.ldecayvtx_x);
  tree->Branch("LambdaDecayVtx_y", &event.ldecayvtx_y);
  tree->Branch("LambdaDecayVtx_z", &event.ldecayvtx_z);
  tree->Branch("LambdaMom_x", &event.lmom_x);
  tree->Branch("LambdaMom_y", &event.lmom_y);
  tree->Branch("LambdaMom_z", &event.lmom_z);
  tree->Branch("LambdaVtxCloseDist", &event.ppi_dist);
  tree->Branch("DecaysTrackId", &event.decays_id);
  tree->Branch("DecaysPurity", &event.decays_purity);
  tree->Branch("DecaysEfficiency", &event.decays_efficiency);
  tree->Branch("DecaysG4tid", &event.decays_G4tid);
  tree->Branch("DecaysMom", &event.decays_mom);
  tree->Branch("DecaysMom_x", &event.decays_mom_x);
  tree->Branch("DecaysMom_y", &event.decays_mom_y);
  tree->Branch("DecaysMom_z", &event.decays_mom_z);
  tree->Branch("DecaysMomRes", &event.decays_res_mom);
  tree->Branch("DecaysMomRes_x", &event.decays_res_mom_x);
  tree->Branch("DecaysMomRes_y", &event.decays_res_mom_y);
  tree->Branch("DecaysMomRes_z", &event.decays_res_mom_z);
  tree->Branch("DecaysMomRes_t", &event.decays_res_mom_t);
  tree->Branch("DecaysThRes", &event.decays_res_th);
  tree->Branch("DecaysPhRes", &event.decays_res_ph);
  tree->Branch("DecaysThCov", &event.decays_cov_mom_th);
  tree->Branch("DecaysPhCov", &event.decays_cov_mom_ph);
  tree->Branch("DecaysMomCov_xy", &event.decays_cov_mom_xy);
  tree->Branch("DecaysMomCov_yz", &event.decays_cov_mom_yz);
  tree->Branch("DecaysMomCov_zx", &event.decays_cov_mom_zx);


  tree->Branch("GFXiflag", &event.GFxiflag);
  tree->Branch("GFXiMass", &event.GFximass);
  tree->Branch("GFXiDecayVtx_x", &event.GFxidecayvtx_x);
  tree->Branch("GFXiDecayVtx_y", &event.GFxidecayvtx_y);
  tree->Branch("GFXiDecayVtx_z", &event.GFxidecayvtx_z);
  tree->Branch("GFXiMom_x", &event.GFximom_x);
  tree->Branch("GFXiMom_y", &event.GFximom_y);
  tree->Branch("GFXiMom_z", &event.GFximom_z);
  tree->Branch("GFXiVtxCloseDist", &event.GFlpi_dist);
  tree->Branch("GFLambdaflag", &event.GFlflag);
  tree->Branch("GFLambdaMass", &event.GFlmass);
  tree->Branch("GFLambdaDecayVtx_x", &event.GFldecayvtx_x);
  tree->Branch("GFLambdaDecayVtx_y", &event.GFldecayvtx_y);
  tree->Branch("GFLambdaDecayVtx_z", &event.GFldecayvtx_z);
  tree->Branch("GFLambdaMom_x", &event.GFlmom_x);
  tree->Branch("GFLambdaMom_y", &event.GFlmom_y);
  tree->Branch("GFLambdaMom_z", &event.GFlmom_z);
  tree->Branch("GFLambdaVtxCloseDist", &event.GFppi_dist);
  tree->Branch("GFDecaysMass", &event.GFdecays_mass);
  tree->Branch("GFDecaysMom", &event.GFdecays_mom);
  tree->Branch("GFDecaysMom_x", &event.GFdecays_mom_x);
  tree->Branch("GFDecaysMom_y", &event.GFdecays_mom_y);
  tree->Branch("GFDecaysMom_z", &event.GFdecays_mom_z);
  tree->Branch("GFMomLoss", &event.GFmomloss);
  tree->Branch("GFELoss", &event.GFeloss);

  //track fitting results
  tree->Branch("GFstatus", &event.GFstatus);
  tree->Branch("GFntTpc", &event.GFntTpc);
  tree->Branch("GFcharge", &event.GFcharge);
  tree->Branch("GFchisqr", &event.GFchisqr);
  tree->Branch("GFchisqrPos", &event.GFchisqrPos);
  tree->Branch("GFtof", &event.GFtof);
  tree->Branch("GFtracklen", &event.GFtracklen);
  tree->Branch("GFpval", &event.GFpval);
  tree->Branch("GFpvalPos", &event.GFpvalPos);
  tree->Branch("GFfitstatus", &event.GFfitstatus);
  tree->Branch("GFpdgcode", &event.GFpdgcode);
  tree->Branch("GFnhtrack", &event.GFnhtrack);
  tree->Branch("GFlayer", &event.GFlayer);
  tree->Branch("GFrow", &event.GFrow);
  tree->Branch("GFpos_x", &event.GFpos_x);
  tree->Branch("GFpos_y", &event.GFpos_y);
  tree->Branch("GFpos_z", &event.GFpos_z);
  tree->Branch("GFmom", &event.GFmom);
  tree->Branch("GFmom_x", &event.GFmom_x);
  tree->Branch("GFmom_y", &event.GFmom_y);
  tree->Branch("GFmom_z", &event.GFmom_z);
  tree->Branch("GFresidual_x", &event.GFresidual_x);
  tree->Branch("GFresidual_y", &event.GFresidual_y);
  tree->Branch("GFresidual_z", &event.GFresidual_z);
  tree->Branch("GFresidual_p", &event.GFresidual_p);
  tree->Branch("GFresidual_px", &event.GFresidual_px);
  tree->Branch("GFresidual_py", &event.GFresidual_py);
  tree->Branch("GFresidual_pz", &event.GFresidual_pz);
  tree->Branch("GFresolution_x", &event.GFresolution_x);
  tree->Branch("GFresolution_y", &event.GFresolution_y);
  tree->Branch("GFresolution_z", &event.GFresolution_z);
  tree->Branch("GFresolution_p", &event.GFresolution_p);
  tree->Branch("GFresolution_px", &event.GFresolution_px);
  tree->Branch("GFresolution_py", &event.GFresolution_py);
  tree->Branch("GFresolution_pz", &event.GFresolution_pz);
  tree->Branch("GFpull_x", &event.GFpull_x);
  tree->Branch("GFpull_y", &event.GFpull_y);
  tree->Branch("GFpull_z", &event.GFpull_z);
  tree->Branch("GFpull_p", &event.GFpull_p);
  tree->Branch("GFpull_px", &event.GFpull_px);
  tree->Branch("GFpull_py", &event.GFpull_py);
  tree->Branch("GFpull_pz", &event.GFpull_pz);
#if DoKinematicFitLdXi

  tree->Branch("KFXiVtxCloseDist", &event.KFlpi_dist);
  tree->Branch("KFchisqrXi",&event.KFchisqrxi);
  tree->Branch("KFpvalXi",&event.KFpvalxi);
  tree->Branch("KFXimom",&event.KFximom);
  tree->Branch("KFXimom_x",&event.KFximom_x);
  tree->Branch("KFXimom_y",&event.KFximom_y);
  tree->Branch("KFXimom_z",&event.KFximom_z);

  tree->Branch("KFchisqrLambda",&event.KFchisqrl);
  tree->Branch("KFpvalLambda",&event.KFpvall);
  tree->Branch("KFLambdamom",&event.KFlmom);
  tree->Branch("KFLambdamom_x",&event.KFlmom_x);
  tree->Branch("KFLambdamom_y",&event.KFlmom_y);
  tree->Branch("KFLambdamom_z",&event.KFlmom_z);
  tree->Branch("KFDecaysMom", &event.KFdecays_mom);
  tree->Branch("KFDecaysMom_x", &event.KFdecays_mom_x);
  tree->Branch("KFDecaysMom_y", &event.KFdecays_mom_y);
  tree->Branch("KFDecaysMom_z", &event.KFdecays_mom_z);


  tree->Branch("XiFlight", &event.XiFlight);
  tree->Branch("XiProd", &event.XiProd);
  tree->Branch("xtgtXi", &event.xtgtXi);
  tree->Branch("ytgtXi", &event.ytgtXi);
  tree->Branch("utgtXi", &event.utgtXi);
  tree->Branch("vtgtXi", &event.vtgtXi);
  tree->Branch("vtxKKXi", &event.vtxKKXi);
  tree->Branch("vtyKKXi", &event.vtyKKXi);
  tree->Branch("vtzKKXi", &event.vtzKKXi);
  tree->Branch("xiprodvtx_x", &event.xiprodvtx_x);
  tree->Branch("xiprodvtx_y", &event.xiprodvtx_y);
  tree->Branch("xiprodvtx_z", &event.xiprodvtx_z);
  tree->Branch("xiprodmom_x", &event.xiprodmom_x);
  tree->Branch("xiprodmom_y", &event.xiprodmom_y);
  tree->Branch("xiprodmom_z", &event.xiprodmom_z);

#endif

  TTreeReaderCont[kE42] = new TTreeReader( "tpc", TFileCont[kE42] );
  const auto& reader = TTreeReaderCont[kE42];
  src.evnum = new TTreeReaderValue<Int_t>( *reader, "evnum" );

  TTreeCont[kE42]->SetBranchAddress("nhHtof",&src.nhHtof);
  TTreeCont[kE42]->SetBranchAddress("didHtof",src.didHtof);
  TTreeCont[kE42]->SetBranchAddress("tHtof",src.tHtof);
  TTreeCont[kE42]->SetBranchAddress("deHtof",src.deHtof);
  TTreeCont[kE42]->SetBranchAddress("yHtof",src.yHtof);
  TTreeCont[kE42]->SetBranchAddress("nhFtof",&src.nhFtof);
  TTreeCont[kE42]->SetBranchAddress("didFtof",src.didFtof);
  TTreeCont[kE42]->SetBranchAddress("tFtof",src.tFtof);
  TTreeCont[kE42]->SetBranchAddress("deFtof",src.deFtof);
  TTreeCont[kE42]->SetBranchAddress("yFtof",src.yFtof);
  TTreeCont[kE42]->SetBranchAddress("nhittpc",&src.nhittpc);
  TTreeCont[kE42]->SetBranchAddress("ititpc",src.ititpc);
  TTreeCont[kE42]->SetBranchAddress("xtpc",src.xtpc);
  TTreeCont[kE42]->SetBranchAddress("ytpc",src.ytpc);
  TTreeCont[kE42]->SetBranchAddress("ztpc",src.ztpc);
  TTreeCont[kE42]->SetBranchAddress("pxtpc",src.pxtpc);
  TTreeCont[kE42]->SetBranchAddress("pytpc",src.pytpc);
  TTreeCont[kE42]->SetBranchAddress("pztpc",src.pztpc);

  TTreeCont[kE42]->SetBranchAddress("NumberOfTracks",&src.NumberOfTracks);
  TTreeCont[kE42]->SetBranchAddress("PIDOfTrack",src.PIDOfTrack);
  TTreeCont[kE42]->SetBranchAddress("ParentIDOfTrack",src.ParentIDOfTrack);
  TTreeCont[kE42]->SetBranchAddress("VertexOfTrack_x",src.VertexOfTrack_x);
  TTreeCont[kE42]->SetBranchAddress("VertexOfTrack_y",src.VertexOfTrack_y);
  TTreeCont[kE42]->SetBranchAddress("VertexOfTrack_z",src.VertexOfTrack_z);
  TTreeCont[kE42]->SetBranchAddress("MomentumOfTrack",src.MomentumOfTrack);
  TTreeCont[kE42]->SetBranchAddress("MomentumOfTrack_x",src.MomentumOfTrack_x);
  TTreeCont[kE42]->SetBranchAddress("MomentumOfTrack_y",src.MomentumOfTrack_y);
  TTreeCont[kE42]->SetBranchAddress("MomentumOfTrack_z",src.MomentumOfTrack_z);

  src.ntTpc = new TTreeReaderValue<Int_t>( *reader, "ntTpc" );
  src.nhtrack = new TTreeReaderValue<std::vector<Int_t>>( *reader, "nhtrack" );
  src.trackid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trackid" );
  src.isBeam = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isBeam" );
  src.isKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isKurama" );
  src.isK18 = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isK18" );
  src.isAccidental = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isAccidental" );
  src.charge = new TTreeReaderValue<std::vector<Int_t>>( *reader, "charge" );
  src.pid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "pid" );
  src.purity = new TTreeReaderValue<std::vector<Double_t>>( *reader, "purity" );
  src.efficiency = new TTreeReaderValue<std::vector<Double_t>>( *reader, "efficiency" );
  src.G4tid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "G4tid" );
  src.chisqr = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqr" );
  src.pval = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pval" );
  src.helix_cx = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_cx" );
  src.helix_cy = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_cy" );
  src.helix_z0 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_z0" );
  src.helix_r = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_r" );
  src.helix_dz = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_dz" );
  src.dE = new TTreeReaderValue<std::vector<Double_t>>( *reader, "dE" );
  src.dEdx = new TTreeReaderValue<std::vector<Double_t>>( *reader, "dEdx" );
  src.mom0 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "mom0" );
  src.path = new TTreeReaderValue<std::vector<Double_t>>( *reader, "path" );
  src.hitlayer = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "hitlayer" );
  src.hitpos_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "hitpos_x" );
  src.hitpos_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "hitpos_y" );
  src.hitpos_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "hitpos_z" );
  src.calpos_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "calpos_x" );
  src.calpos_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "calpos_y" );
  src.calpos_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "calpos_z" );
  src.residual = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "residual" );
  src.residual_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "residual_x" );
  src.residual_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "residual_y" );
  src.residual_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "residual_z" );
  src.resolution_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "resolution_x" );
  src.resolution_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "resolution_y" );
  src.resolution_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "resolution_z" );
  src.helix_t = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "helix_t" );
  src.alpha = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "alpha" );
  src.track_cluster_de = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_de" );
  src.track_cluster_mrow = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_mrow" );
  src.xtgtHS = new TTreeReaderValue<vector<Double_t>>( *reader, "xtgtHS" );
  src.ytgtHS = new TTreeReaderValue<vector<Double_t>>( *reader, "ytgtHS" );
  src.xtgtKurama = new TTreeReaderValue<vector<Double_t>>( *reader, "xtgtKurama" );
  src.ytgtKurama = new TTreeReaderValue<vector<Double_t>>( *reader, "ytgtKurama" );



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
      InitializeParameter<FieldMan>("FLDMAP", "HSFLDMAP") &&
      InitializeParameter<HodoPHCMan>("HDPHC") );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess( void )
{
  return true;
}

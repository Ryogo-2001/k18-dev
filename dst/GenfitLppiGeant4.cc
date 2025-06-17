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
#include "TPCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DstHelper.hh"
#include "HodoPHCMan.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCCluster.hh"
#include "TPCVertex.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"
#include "FourVectorFitter.hh"

#include "HypTPCFitter.hh"
#include "HypTPCTask.hh"

#define SaveRawData 1
#define DebugDisp 0
#define SaveHistograms 0
#define KinematicFit 0

namespace
{
using namespace root;
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
const bool Const_field = false; //Must be false for linear tracking
const Int_t verbosity = 0;//0~3;
//const Int_t verbosity = 1;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");

const Double_t vtx_scan_range = 150.; //ref
const Double_t vtx_scan_rangeInsideL = 50.;
const Double_t vtx_scan_rangeInsidePi = 50.;

const Double_t xi_masscut = 0.1; const Double_t lambda_masscut = 0.2; //final
//const Double_t xi_masscut = 0.15; const Double_t lambda_masscut = 0.1; //ref
const Double_t p_vtx_distcut = 300;
const Double_t pi_vtx_distcut = 300;
const Double_t pi2_vtx_distcut = 300;
const Double_t pipi_distcut = 10.; //ref
const Double_t ppi_distcut = 10.; //ref
//const Double_t lpi_distcut = 10.; //ref
const Double_t lpi_distcut = 15.; //ref
const Double_t xitarget_distcut = 50.; //ref
const Double_t ltarget_distcut = 50.;

const Double_t GFppi_distcut = 10.;
const Double_t GFlpi_distcut = 10.;
//const Double_t GFlpi_distcut = 15.;
const Double_t GFxitarget_distcut = 50.;
const Double_t GFltarget_distcut = 50.;
const Double_t GFxitarget_ycut = 20.;
const Double_t GFltarget_ycut = 20.;

//Measured resolutions for multi-track vertexing
const Double_t duCh2 = 0.001381;
const Double_t dvCh2 = duCh2;
const Double_t duDiamond = 0.2796;
const Double_t dvDiamond = duDiamond;

Double_t res_uK18 = 0.00288/sqrt(2);
Double_t res_xK18 = 1.432/sqrt(2);
Double_t res_vK18 = 0.00334/sqrt(2);
Double_t res_yK18 = 2.836/sqrt(2);
Double_t res_uKurama = hypot(res_uK18,duCh2);
Double_t res_xKurama = res_xK18;
Double_t res_vKurama = hypot(res_vK18,dvCh2);
Double_t res_yKurama = res_yK18;
Double_t res_xXiVtx = 0.6;
Double_t res_yXiVtx = 0.5;
Double_t res_xLdVtx = 0.5;
Double_t res_yLdVtx = 0.5;

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
{ "[Process]", "[ConfFile]", "[DstTPCTrackingHelixgeant4]", "[OutFile]" };
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
  std::vector<Int_t> trigpat;
  std::vector<Int_t> trigflag;

  Int_t nhHtof;
  std::vector<Double_t> HtofSeg;
  std::vector<Double_t> tHtof;
  std::vector<Double_t> dtHtof;
  std::vector<Double_t> deHtof;
  std::vector<Double_t> posHtof;
  std::vector<Int_t> G4tidHtof;

  Int_t nhFtof;
  std::vector<Double_t> FtofSeg;
  std::vector<Double_t> tFtof;
  std::vector<Double_t> deFtof;
  std::vector<Double_t> posFtof;

  Int_t ntK18;
  std::vector<Double_t> pK18;
  std::vector<Double_t> chisqrK18;
  std::vector<Double_t> xtgtK18;
  std::vector<Double_t> ytgtK18;
  std::vector<Double_t> utgtK18;
  std::vector<Double_t> vtgtK18;

  Int_t ntKurama;
  std::vector<Double_t> chisqrKurama;
  std::vector<Double_t> pKurama;
  std::vector<Double_t> qKurama;
  std::vector<Double_t> m2;
  std::vector<Double_t> thetaKurama;
  std::vector<Double_t> xtgtKurama;
  std::vector<Double_t> ytgtKurama;
  std::vector<Double_t> utgtKurama;
  std::vector<Double_t> vtgtKurama;
  std::vector<Double_t> pathwcKurama;
  std::vector<Double_t> xin;
  std::vector<Double_t> yin;
  std::vector<Double_t> zin;
  std::vector<Double_t> pxin;
  std::vector<Double_t> pyin;
  std::vector<Double_t> pzin;
  std::vector<Double_t> xout;
  std::vector<Double_t> yout;
  std::vector<Double_t> zout;
  std::vector<Double_t> pxout;
  std::vector<Double_t> pyout;
  std::vector<Double_t> pzout;

  Int_t nKm;
  Int_t nKp;
  Int_t nKK;
  std::vector<Int_t> inside;
  std::vector<Double_t> vtx;
  std::vector<Double_t> vty;
  std::vector<Double_t> vtz;
  std::vector<Double_t> closeDist;
  std::vector<Double_t> MissMass;
  std::vector<Double_t> MissMassCorr;
  std::vector<Double_t> MissMassCorrDE;
  std::vector<Double_t> pOrg;
  std::vector<Double_t> pCalc;
  std::vector<Double_t> pCorr;
  std::vector<Double_t> pCorrDE;
  std::vector<Double_t> xb;
  std::vector<Double_t> yb;
  std::vector<Double_t> ub;
  std::vector<Double_t> vb;
  std::vector<Double_t> xs;
  std::vector<Double_t> ys;
  std::vector<Double_t> us;
  std::vector<Double_t> vs;
  std::vector<Int_t> Kflag;

  //TPC RK
  std::vector<Int_t> isgoodTPCK18;
  std::vector<Double_t> chisqrTPCK18;
  std::vector<Double_t> pTPCK18;
  std::vector<Double_t> qTPCK18;
  std::vector<Double_t> xtgtTPCK18;
  std::vector<Double_t> ytgtTPCK18;
  std::vector<Double_t> utgtTPCK18;
  std::vector<Double_t> vtgtTPCK18;
  std::vector<Double_t> thetaTPCK18;
  std::vector<Double_t> lhtofTPCK18;
  std::vector<Double_t> xhtofTPCK18;
  std::vector<Double_t> yhtofTPCK18;
  std::vector<std::vector<Double_t>> lvpTPCK18;
  std::vector<std::vector<Double_t>> xvpTPCK18;
  std::vector<std::vector<Double_t>> yvpTPCK18;

  std::vector<Int_t> isgoodTPCKurama;
  std::vector<Int_t> kflagTPCKurama;
  std::vector<Double_t> chisqrTPCKurama;
  std::vector<Double_t> pTPCKurama;
  std::vector<Double_t> qTPCKurama;
  std::vector<Double_t> m2TPCKurama;
  std::vector<Double_t> xtgtTPCKurama;
  std::vector<Double_t> ytgtTPCKurama;
  std::vector<Double_t> utgtTPCKurama;
  std::vector<Double_t> vtgtTPCKurama;
  std::vector<Double_t> thetaTPCKurama;
  std::vector<Double_t> pathTPCKurama;
  std::vector<Double_t> lhtofTPCKurama;
  std::vector<Double_t> xhtofTPCKurama;
  std::vector<Double_t> yhtofTPCKurama;
  std::vector<std::vector<Double_t>> lvpTPCKurama;
  std::vector<std::vector<Double_t>> xvpTPCKurama;
  std::vector<std::vector<Double_t>> yvpTPCKurama;

  std::vector<Int_t> isgoodTPC;
  std::vector<Int_t> insideTPC;
  std::vector<Double_t> vtxTPC;
  std::vector<Double_t> vtyTPC;
  std::vector<Double_t> vtzTPC;
  std::vector<Double_t> closeDistTPC;
  std::vector<Double_t> MissMassTPC;
  std::vector<Double_t> MissMassCorrTPC;
  std::vector<Double_t> MissMassCorrDETPC;
  std::vector<Double_t> pOrgTPC;
  std::vector<Double_t> pCalcTPC;
  std::vector<Double_t> pCorrTPC;
  std::vector<Double_t> pCorrDETPC;
  std::vector<Double_t> thetaTPC;
  std::vector<Double_t> thetaCMTPC;
  std::vector<Double_t> costCMTPC;
  std::vector<Double_t> xbTPC;
  std::vector<Double_t> ybTPC;
  std::vector<Double_t> ubTPC;
  std::vector<Double_t> vbTPC;
  std::vector<Double_t> xsTPC;
  std::vector<Double_t> ysTPC;
  std::vector<Double_t> usTPC;
  std::vector<Double_t> vsTPC;

  std::vector<Double_t> BE;
  std::vector<Double_t> BETPC;
  std::vector<Double_t> BE_LL;
  std::vector<Double_t> BETPC_LL;
std::vector<Double_t> km_mom_x;
  std::vector<Double_t> km_mom_y;
  std::vector<Double_t> km_mom_z;
  std::vector<Double_t> kp_mom_x;
  std::vector<Double_t> kp_mom_y;
  std::vector<Double_t> kp_mom_z;
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

  Int_t nclTpc; // Number of clusters
  std::vector<Double_t> cluster_x;
  std::vector<Double_t> cluster_y;
  std::vector<Double_t> cluster_z;
  std::vector<Double_t> cluster_de;
  std::vector<Int_t> cluster_size;
  std::vector<Int_t> cluster_layer;
  std::vector<Double_t> cluster_mrow;
  std::vector<Double_t> cluster_de_center;
  std::vector<Int_t> cluster_houghflag;
  std::vector<Int_t> cluster_G4tid;
  std::vector<Int_t> cluster_G4pid;

  Int_t remain_nclTpc; // Number of remain clusters not occupied in the tracks
  std::vector<Double_t> remain_cluster_x;
  std::vector<Double_t> remain_cluster_y;
  std::vector<Double_t> remain_cluster_z;
  std::vector<Double_t> remain_cluster_de;
  std::vector<Int_t> remain_cluster_size;
  std::vector<Int_t> remain_cluster_layer;
  std::vector<Int_t> remain_cluster_houghflag;
  std::vector<Int_t> remain_cluster_G4tid;
  std::vector<Int_t> remain_cluster_G4pid;

  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Int_t> trackid; //for Kurama K1.8 tracks
  std::vector<Int_t> isXi;
  std::vector<Int_t> isBeam;
  std::vector<Int_t> isKurama;
  std::vector<Int_t> isK18;
  std::vector<Int_t> isAccidental;
  std::vector<Int_t> isMultiloop;
  std::vector<Int_t> isInTarget;
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
  std::vector<std::vector<Double_t>> mom_x;
  std::vector<std::vector<Double_t>> mom_y;
  std::vector<std::vector<Double_t>> mom_z;
  std::vector<std::vector<Double_t>> residual;
  std::vector<std::vector<Double_t>> residual_x;
  std::vector<std::vector<Double_t>> residual_y;
  std::vector<std::vector<Double_t>> residual_z;
  std::vector<std::vector<Double_t>> resolution_x;
  std::vector<std::vector<Double_t>> resolution_y;
  std::vector<std::vector<Double_t>> resolution_z;
  std::vector<std::vector<Double_t>> helix_t;
  std::vector<std::vector<Double_t>> pathhit;
  std::vector<std::vector<Double_t>> alpha;
  std::vector<std::vector<Double_t>> track_cluster_de;
  std::vector<std::vector<Double_t>> track_cluster_size;
  std::vector<std::vector<Double_t>> track_cluster_mrow;

  std::vector<Int_t> chargeIndistinguishable;
  std::vector<Double_t> chisqr_inverted;
  std::vector<Double_t> pval_inverted;
  std::vector<Double_t> helix_cx_inverted;
  std::vector<Double_t> helix_cy_inverted;
  std::vector<Double_t> helix_z0_inverted;
  std::vector<Double_t> helix_r_inverted;
  std::vector<Double_t> helix_dz_inverted;
  std::vector<Double_t> mom0_inverted;//Helix momentum at Y = 0
  std::vector<Int_t> pid_inverted;

  Int_t nvtxTpc;
  std::vector<Double_t> vtx_x;
  std::vector<Double_t> vtx_y;
  std::vector<Double_t> vtx_z;
  std::vector<Double_t> vtx_dist;
  std::vector<Double_t> vtx_angle;
  std::vector<std::vector<Double_t>> vtxid;
  std::vector<std::vector<Double_t>> vtxmom_theta;
  std::vector<std::vector<Double_t>> vtxpos_x;
  std::vector<std::vector<Double_t>> vtxpos_y;
  std::vector<std::vector<Double_t>> vtxpos_z;
  std::vector<std::vector<Double_t>> vtxmom_x;
  std::vector<std::vector<Double_t>> vtxmom_y;
  std::vector<std::vector<Double_t>> vtxmom_z;

  std::vector<Int_t> isLambda;
  std::vector<Int_t> ncombiLambda;
  std::vector<Double_t> distLambda;
  std::vector<Double_t> angleLambda;
  std::vector<Double_t> bestmassLambda;
  std::vector<std::vector<Double_t>> massLambda;
  std::vector<std::vector<Double_t>> vtxLambda_x;
  std::vector<std::vector<Double_t>> vtxLambda_y;
  std::vector<std::vector<Double_t>> vtxLambda_z;
  std::vector<std::vector<Double_t>> momLambda;
  std::vector<std::vector<Double_t>> momLambda_x;
  std::vector<std::vector<Double_t>> momLambda_y;
  std::vector<std::vector<Double_t>> momLambda_z;
  std::vector<std::vector<Double_t>> decaysidLambda;
  std::vector<std::vector<Double_t>> decaysmomLambda;
  std::vector<std::vector<Double_t>> decaysmomLambda_x;
  std::vector<std::vector<Double_t>> decaysmomLambda_y;
  std::vector<std::vector<Double_t>> decaysmomLambda_z;

  Int_t nvtxTpcClustered;
  std::vector<Double_t> clusteredVtx_x;
  std::vector<Double_t> clusteredVtx_y;
  std::vector<Double_t> clusteredVtx_z;
  std::vector<std::vector<Double_t>> clusteredVtxid;

  Int_t ncombiLreconfailed;
  std::vector<Int_t> pidLreconfailed;
  std::vector<Int_t> piidLreconfailed;
  std::vector<Double_t> LdecayvtxLreconfailed_x;
  std::vector<Double_t> LdecayvtxLreconfailed_y;
  std::vector<Double_t> LdecayvtxLreconfailed_z;
  std::vector<Double_t> LmassLreconfailed;
  std::vector<Double_t> LmomLreconfailed;
  std::vector<Double_t> LmomLreconfailed_x;
  std::vector<Double_t> LmomLreconfailed_y;
  std::vector<Double_t> LmomLreconfailed_z;
  std::vector<Double_t> pmomLreconfailed;
  std::vector<Double_t> pmomLreconfailed_x;
  std::vector<Double_t> pmomLreconfailed_y;
  std::vector<Double_t> pmomLreconfailed_z;
  std::vector<Double_t> pimomLreconfailed;
  std::vector<Double_t> pimomLreconfailed_x;
  std::vector<Double_t> pimomLreconfailed_y;
  std::vector<Double_t> pimomLreconfailed_z;
  std::vector<Double_t> ppidistLreconfailed;

  Int_t ncombiPipair;
  std::vector<Int_t> pipidPipair;
  std::vector<Int_t> pimidPipair;
  std::vector<Double_t> pipmomPipair;
  std::vector<Double_t> pipmomPipair_x;
  std::vector<Double_t> pipmomPipair_y;
  std::vector<Double_t> pipmomPipair_z;
  std::vector<Double_t> pimmomPipair;
  std::vector<Double_t> pimmomPipair_x;
  std::vector<Double_t> pimmomPipair_y;
  std::vector<Double_t> pimmomPipair_z;
  std::vector<Double_t> momPipair;
  std::vector<Double_t> momPipair_x;
  std::vector<Double_t> momPipair_y;
  std::vector<Double_t> momPipair_z;
  std::vector<Double_t> reconLmassPipair;
  std::vector<Double_t> pipidistPipair;

  //Multi-track production vertex
  Double_t GFprodvtx_x_ll;
  Double_t GFprodvtx_y_ll;
  Double_t GFprodvtx_z_ll;
  Double_t GFprodvtx_x_l1;
  Double_t GFprodvtx_y_l1;
  Double_t GFprodvtx_z_l1;
  Double_t GFprodvtx_x_l2;
  Double_t GFprodvtx_y_l2;
  Double_t GFprodvtx_z_l2;
  Double_t GFprodvtx_x_l;
  Double_t GFprodvtx_y_l;
  Double_t GFprodvtx_z_l;

  Bool_t emptyflag;
  Bool_t pimflag;
  Bool_t lpiflag;
  Bool_t lflag;

  Bool_t llflag;
  Double_t ltarget_dist1;
  Double_t ltargetvtx_x1;
  Double_t ltargetvtx_y1;
  Double_t ltargetvtx_z1;
  Double_t lmass1;
  Double_t ldecayvtx_x1;
  Double_t ldecayvtx_y1;
  Double_t ldecayvtx_z1;
  Double_t lmom1;
  Double_t lmom_x1;
  Double_t lmom_y1;
  Double_t lmom_z1;
  Double_t ppi_dist1;
  Double_t ltarget_dist2;
  Double_t ltargetvtx_x2;
  Double_t ltargetvtx_y2;
  Double_t ltargetvtx_z2;
  Double_t lmass2;
  Double_t ldecayvtx_x2;
  Double_t ldecayvtx_y2;
  Double_t ldecayvtx_z2;
  Double_t lmom2;
  Double_t lmom_x2;
  Double_t lmom_y2;
  Double_t lmom_z2;
  Double_t ppi_dist2;

  Double_t GFllexcitation;
  Double_t GFlmass1;
  Double_t GFldecayvtx_x1;
  Double_t GFldecayvtx_y1;
  Double_t GFldecayvtx_z1;
  Double_t GFlmom1;
  Double_t GFlmom_x1;
  Double_t GFlmom_y1;
  Double_t GFlmom_z1;
  Double_t GFppi_dist1;
  Double_t GFltarget_dist1;
  Double_t GFltargetvtx_x1;
  Double_t GFltargetvtx_y1;
  Double_t GFltargetvtx_z1;
  Double_t GFltargetcenter_x1;
  Double_t GFltargetcenter_y1;
  Double_t GFltargetcenter_z1;
  Double_t GFltargetcenter_dist1;
  Double_t GFlprodvtx_x1;
  Double_t GFlprodvtx_y1;
  Double_t GFlprodvtx_z1;
  Double_t GFlprodvtx_dist1;
  Double_t GFltracklen1;
  Double_t GFltof1;

  Double_t GFlmass2;
  Double_t GFldecayvtx_x2;
  Double_t GFldecayvtx_y2;
  Double_t GFldecayvtx_z2;
  Double_t GFlmom2;
  Double_t GFlmom_x2;
  Double_t GFlmom_y2;
  Double_t GFlmom_z2;
  Double_t GFppi_dist2;
  Double_t GFltarget_dist2;
  Double_t GFltargetvtx_x2;
  Double_t GFltargetvtx_y2;
  Double_t GFltargetvtx_z2;
  Double_t GFltargetcenter_x2;
  Double_t GFltargetcenter_y2;
  Double_t GFltargetcenter_z2;
  Double_t GFltargetcenter_dist2;
  Double_t GFlprodvtx_x2;
  Double_t GFlprodvtx_y2;
  Double_t GFlprodvtx_z2;
  Double_t GFlprodvtx_dist2;
  Double_t GFltracklen2;
  Double_t GFltof2;

  Double_t llvtx_x;
  Double_t llvtx_y;
  Double_t llvtx_z;
  Double_t lldist;
  Double_t GFllvtx_x;
  Double_t GFllvtx_y;
  Double_t GFllvtx_z;
  Double_t GFlldist;

  Double_t KFllexcitation;
  Double_t KFlmom1;
  Double_t KFlmom_x1;
  Double_t KFlmom_y1;
  Double_t KFlmom_z1;
  Double_t KFlmom2;
  Double_t KFlmom_x2;
  Double_t KFlmom_y2;
  Double_t KFlmom_z2;
  Double_t KFlchisqr1;
  Double_t KFlchisqr2;

  Double_t KFprodvtx_chisqr_ll;
  Double_t KFprodvtx_x_ll;
  Double_t KFprodvtx_y_ll;
  Double_t KFprodvtx_z_ll;
  Double_t KFprodvtx_x_l1;
  Double_t KFprodvtx_y_l1;
  Double_t KFprodvtx_z_l1;
  Double_t KFprodvtx_x_l2;
  Double_t KFprodvtx_y_l2;
  Double_t KFprodvtx_z_l2;
  Double_t KFprodvtx_x_l;
  Double_t KFprodvtx_y_l;
  Double_t KFprodvtx_z_l;

  Double_t KFllvtx_x;
  Double_t KFllvtx_y;
  Double_t KFllvtx_z;
  Double_t KFlldist;

  Double_t KFlprodvtx_x1;
  Double_t KFlprodvtx_y1;
  Double_t KFlprodvtx_z1;
  Double_t KFlprodvtx_dist1;
  Double_t KFltracklen1;
  Double_t KFltof1;

  Double_t KFlprodvtx_x2;
  Double_t KFlprodvtx_y2;
  Double_t KFlprodvtx_z2;
  Double_t KFlprodvtx_dist2;
  Double_t KFltracklen2;
  Double_t KFltof2;

  Double_t KFlprodvtx_x;
  Double_t KFlprodvtx_y;
  Double_t KFlprodvtx_z;
  Double_t KFlprodvtx_dist;
  Double_t KFltracklen;
  Double_t KFltof;

  Bool_t xiflag;
  Double_t ximass;
  Double_t xidecayvtx_x;
  Double_t xidecayvtx_y;
  Double_t xidecayvtx_z;
  Double_t ximom;
  Double_t ximom_x;
  Double_t ximom_y;
  Double_t ximom_z;
  Double_t lpi_dist;
  Double_t xitargetvtx_x;
  Double_t xitargetvtx_y;
  Double_t xitargetvtx_z;
  Double_t xitargetmom;
  Double_t xitargetmom_x;
  Double_t xitargetmom_y;
  Double_t xitargetmom_z;
  Double_t xitarget_dist;

  Double_t lmass;
  Int_t l_intarget;
  Double_t ldecayvtx_x;
  Double_t ldecayvtx_y;
  Double_t ldecayvtx_z;
  Double_t lmom;
  Double_t lmom_x;
  Double_t lmom_y;
  Double_t lmom_z;
  Double_t ppi_dist;
  Double_t ltarget_dist;
  Double_t ltargetvtx_x;
  Double_t ltargetvtx_y;
  Double_t ltargetvtx_z;

  Double_t lmass_vtx;
  Double_t ldecayvtx_x_vtx;
  Double_t ldecayvtx_y_vtx;
  Double_t ldecayvtx_z_vtx;
  Double_t lmom_vtx;
  Double_t lmom_x_vtx;
  Double_t lmom_y_vtx;
  Double_t lmom_z_vtx;
  Double_t ppi_dist_vtx;

  Double_t GFximass;
  Double_t GFxidecayvtx_x;
  Double_t GFxidecayvtx_y;
  Double_t GFxidecayvtx_z;
  Double_t GFximom;
  Double_t GFximom_x;
  Double_t GFximom_y;
  Double_t GFximom_z;

  Double_t GFxikkvtx_x;
  Double_t GFxikkvtx_y;
  Double_t GFxikkvtx_z;
  Double_t GFxikkmom;
  Double_t GFxikkmom_x;
  Double_t GFxikkmom_y;
  Double_t GFxikkmom_z;
  Double_t GFxikkvtx_dist;

  Double_t GFxiprodvtx_x;
  Double_t GFxiprodvtx_y;
  Double_t GFxiprodvtx_z;
  Double_t GFxiprodmom;
  Double_t GFxiprodmom_x;
  Double_t GFxiprodmom_y;
  Double_t GFxiprodmom_z;
  Double_t GFxiprodvtx_dist;
  Double_t GFxitracklen;
  Double_t GFxitof;
  Double_t GFlpi_dist;
  Double_t GFximomloss;
  Double_t GFxiexcitation;

  Double_t GFxitargetvtx_x;
  Double_t GFxitargetvtx_y;
  Double_t GFxitargetvtx_z;
  Double_t GFxitargetmom;
  Double_t GFxitargetmom_x;
  Double_t GFxitargetmom_y;
  Double_t GFxitargetmom_z;
  Double_t GFxitarget_dist;

  Double_t GFxitargetcenter_x;
  Double_t GFxitargetcenter_y;
  Double_t GFxitargetcenter_z;
  Double_t GFxitargetcentermom;
  Double_t GFxitargetcentermom_x;
  Double_t GFxitargetcentermom_y;
  Double_t GFxitargetcentermom_z;
  Double_t GFxitargetcenter_dist;

  Double_t GFlmass;
  Double_t GFldecayvtx_x;
  Double_t GFldecayvtx_y;
  Double_t GFldecayvtx_z;
  Double_t GFlmom;
  Double_t GFlmom_x;
  Double_t GFlmom_y;
  Double_t GFlmom_z;
  Double_t GFppi_dist;
  Double_t GFltarget_dist;
  Double_t GFltargetvtx_x;
  Double_t GFltargetvtx_y;
  Double_t GFltargetvtx_z;
  Double_t GFltargetcenter_x;
  Double_t GFltargetcenter_y;
  Double_t GFltargetcenter_z;
  Double_t GFltargetcenter_dist;
  Double_t GFlprodvtx_x;
  Double_t GFlprodvtx_y;
  Double_t GFlprodvtx_z;
  Double_t GFlprodvtx_dist;
  Double_t GFltracklen;
  Double_t GFltof;

  //Multi-track production vertex
  Double_t GFprodvtx_x_kkxi;
  Double_t GFprodvtx_y_kkxi;
  Double_t GFprodvtx_z_kkxi;

  Bool_t lphiflag;
  Double_t phimass;
  Double_t phidecayvtx_x;
  Double_t phidecayvtx_y;
  Double_t phidecayvtx_z;
  Double_t phimom;
  Double_t phimom_x;
  Double_t phimom_y;
  Double_t phimom_z;
  Double_t kk_dist;
  Double_t phi_km_mass2;

  Double_t GFphimass;
  Double_t GFphidecayvtx_x;
  Double_t GFphidecayvtx_y;
  Double_t GFphidecayvtx_z;
  Double_t GFphimom;
  Double_t GFphimom_x;
  Double_t GFphimom_y;
  Double_t GFphimom_z;
  Double_t GFkk_dist;
  Double_t GFphiprodvtx_dist;
  std::vector<Int_t> phidecays_id;
  std::vector<Double_t> phidecays_mom;
  std::vector<Double_t> phidecays_mom_x;
  std::vector<Double_t> phidecays_mom_y;
  std::vector<Double_t> phidecays_mom_z;
  std::vector<Double_t> GFphidecays_mom;
  std::vector<Double_t> GFphidecays_mom_x;
  std::vector<Double_t> GFphidecays_mom_y;
  std::vector<Double_t> GFphidecays_mom_z;

  Int_t GFntdecays;
  std::vector<Int_t> GFdecays_htofid;
  std::vector<Double_t> GFdecays_tracklen;
  std::vector<Double_t> GFdecays_tof;
  std::vector<Double_t> GFdecays_mass2;
  std::vector<Double_t> GFdecays_mom;
  std::vector<Double_t> GFdecays_mom_x;
  std::vector<Double_t> GFdecays_mom_y;
  std::vector<Double_t> GFdecays_mom_z;
  std::vector<Double_t> GFdecays_CMmom;
  std::vector<Double_t> GFdecays_CMmom_x;
  std::vector<Double_t> GFdecays_CMmom_y;
  std::vector<Double_t> GFdecays_CMmom_z;
  std::vector<Double_t> GFmomloss;
  std::vector<Double_t> GFeloss;

  Int_t GFntTpc;
  std::vector<Int_t> GFpdgcode;
  std::vector<Int_t> GFnhtrack;
  std::vector<Double_t> GFcharge;
  std::vector<Double_t> GFchisqr;
  std::vector<Double_t> GFtof;
  std::vector<Double_t> GFtracklen;
  std::vector<Double_t> GFpval;

  std::vector<Int_t> decays_id;
  std::vector<Double_t> decays_purity;
  std::vector<Double_t> decays_efficiency;
  std::vector<Int_t> decays_G4tid;
  std::vector<Double_t> decays_mom;
  std::vector<Double_t> decays_mom_x;
  std::vector<Double_t> decays_mom_y;
  std::vector<Double_t> decays_mom_z;
  std::vector<Double_t> decays_CMmom;
  std::vector<Double_t> decays_CMmom_x;
  std::vector<Double_t> decays_CMmom_y;
  std::vector<Double_t> decays_CMmom_z;

  Bool_t pipiflag;

  Int_t residual_multi;
  Int_t pim_multi;
  Int_t pip_multi;
  Int_t p_multi;
  Int_t ppip_multi;
  Int_t accident_multi;

  std::vector<Int_t> residual_id;
  std::vector<Double_t> residual_dist2tgt;
  std::vector<Double_t> residual_mass2;
  std::vector<Double_t> residual_mom;
  std::vector<Double_t> residual_mom_x;
  std::vector<Double_t> residual_mom_y;
  std::vector<Double_t> residual_mom_z;
  std::vector<Double_t> residual_charge;

  //Kinematic fitting
  Double_t KFlmom0;
  Double_t KFlmom_x0;
  Double_t KFlmom_y0;
  Double_t KFlmom_z0;
  Double_t KFlmom;
  Double_t KFlmom_x;
  Double_t KFlmom_y;
  Double_t KFlmom_z;
  Double_t KFlchisqr;
  Double_t KFlpval;
  Double_t KFlpi_dist;
  Double_t KFximom;
  Double_t KFximom_x;
  Double_t KFximom_y;
  Double_t KFximom_z;
  Double_t KFxichisqr;
  Double_t KFxipval;
  Double_t KFximass;
  Double_t KFxidecayvtx_x;
  Double_t KFxidecayvtx_y;
  Double_t KFxidecayvtx_z;
  std::vector<Double_t> KFlpull;
  std::vector<Double_t> KFxipull;

  //Multi-track vertex
  Double_t KFprodvtx_chisqr_kkxi;
  Double_t KFprodvtx_x_kkxi;
  Double_t KFprodvtx_y_kkxi;
  Double_t KFprodvtx_z_kkxi;
  Double_t KFprodvtx_x_kpxi;
  Double_t KFprodvtx_y_kpxi;
  Double_t KFprodvtx_z_kpxi;

  Double_t KFxiprodvtx_x;
  Double_t KFxiprodvtx_y;
  Double_t KFxiprodvtx_z;
  Double_t KFxiprodmom;
  Double_t KFxiprodmom_x;
  Double_t KFxiprodmom_y;
  Double_t KFxiprodmom_z;
  Double_t KFxiprodvtx_dist;
  Double_t KFxitracklen;
  Double_t KFxitof;
  Double_t KFximomloss;
  Double_t KFxiexcitation;

  Double_t KFxi_kkvtx_x;
  Double_t KFxi_kkvtx_y;
  Double_t KFxi_kkvtx_z;
  Double_t KFxi_kkvtx_mom;
  Double_t KFxi_kkvtx_mom_x;
  Double_t KFxi_kkvtx_mom_y;
  Double_t KFxi_kkvtx_mom_z;
  Double_t KFxi_kkvtx_dist;

  Double_t KFxi_kpxiprodvtx_x;
  Double_t KFxi_kpxiprodvtx_y;
  Double_t KFxi_kpxiprodvtx_z;
  Double_t KFxi_kpxiprodmom;
  Double_t KFxi_kpxiprodmom_x;
  Double_t KFxi_kpxiprodmom_y;
  Double_t KFxi_kpxiprodmom_z;
  Double_t KFxi_kpxiprodvtx_dist;

  Double_t KFxitargetvtx_x;
  Double_t KFxitargetvtx_y;
  Double_t KFxitargetvtx_z;
  Double_t KFxitargetmom;
  Double_t KFxitargetmom_x;
  Double_t KFxitargetmom_y;
  Double_t KFxitargetmom_z;
  Double_t KFxitarget_dist;

  Double_t KFxitargetcenter_x;
  Double_t KFxitargetcenter_y;
  Double_t KFxitargetcenter_z;
  Double_t KFxitargetcentermom;
  Double_t KFxitargetcentermom_x;
  Double_t KFxitargetcentermom_y;
  Double_t KFxitargetcentermom_z;
  Double_t KFxitargetcenter_dist;

  std::vector<Double_t> KFdecays_mom;
  std::vector<Double_t> KFdecays_mom_x;
  std::vector<Double_t> KFdecays_mom_y;
  std::vector<Double_t> KFdecays_mom_z;
  std::vector<Double_t> KFdecays_CMmom;
  std::vector<Double_t> KFdecays_CMmom_x;
  std::vector<Double_t> KFdecays_CMmom_y;
  std::vector<Double_t> KFdecays_CMmom_z;

  //for p mom correction test
  Double_t GFximass_pmomcorr;
  Double_t GFlmass_pmomcorr;
  Double_t GFximom_pmomcorr;
  Double_t GFlmom_pmomcorr;
  std::vector<Double_t> GFdecays_mom_pmomcorr;
  std::vector<Double_t> GFdecays_mom_x_pmomcorr;
  std::vector<Double_t> GFdecays_mom_y_pmomcorr;
  std::vector<Double_t> GFdecays_mom_z_pmomcorr;
  std::vector<Double_t> GFpval_pmomcorr;
  std::vector<Double_t> GFchisqr_pmomcorr;

  int G4l1id;
  int G4l1vtx_x;
  int G4l1vtx_y;
  int G4l1vtx_z;
  double G4l1mom;
  double G4l1mom_x;
  double G4l1mom_y;
  double G4l1mom_z;

  double l1vtx_x;
  double l1vtx_y;
  double l1vtx_z;

  int G4l2id;
  int G4l2vtx_x;
  int G4l2vtx_y;
  int G4l2vtx_z;
  double G4l2mom;
  double G4l2mom_x;
  double G4l2mom_y;
  double G4l2mom_z;

  double l2vtx_x;
  double l2vtx_y;
  double l2vtx_z;

  double G4llmass;

  int G4p1id;
  int G4p1tid;
  int G4p1nh;
  int G4p1tnh;
  double G4p1vtx_x;
  double G4p1vtx_y;
  double G4p1vtx_z;
  double G4p1mom;
  double G4p1mom_x;
  double G4p1mom_y;
  double G4p1mom_z;

  int p1tid;
  int p1nh;
  double p1vtx_x;
  double p1vtx_y;
  double p1vtx_z;
  double p1mom;
  double p1mom_x;
  double p1mom_y;
  double p1mom_z;
  double GFp1mom;
  double GFp1mom_x;
  double GFp1mom_y;
  double GFp1mom_z;

  int G4p2id;
  int G4p2tid;
  int G4p2nh;
  int G4p2tnh;
  double G4p2vtx_x;
  double G4p2vtx_y;
  double G4p2vtx_z;
  double G4p2mom;
  double G4p2mom_x;
  double G4p2mom_y;
  double G4p2mom_z;

  int p2tid;
  int p2nh;
  double p2vtx_x;
  double p2vtx_y;
  double p2vtx_z;
  double p2mom;
  double p2mom_x;
  double p2mom_y;
  double p2mom_z;
  double GFp2mom;
  double GFp2mom_x;
  double GFp2mom_y;
  double GFp2mom_z;


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
  double GFpi2mom;
  double GFpi2mom_x;
  double GFpi2mom_y;
  double GFpi2mom_z;

  bool l1good,l2good,llswap,piswap;
  bool p1_tracked,p2_tracked,pi1_tracked,pi2_tracked;
  double p1t_mom0,p2t_mom0,pi1t_mom0,pi2t_mom0;

  void clear( void )
  {
    evnum = 0;
    status = 0;
    trigpat.clear();
    trigflag.clear();

    nhHtof = 0;
    HtofSeg.clear();
    tHtof.clear();
    dtHtof.clear();
    deHtof.clear();
    posHtof.clear();
    G4tidHtof.clear();

    ntK18 = 0;
    pK18.clear();
    chisqrK18.clear();
    xtgtK18.clear();
    ytgtK18.clear();
    utgtK18.clear();
    vtgtK18.clear();

    ntKurama = 0;
    chisqrKurama.clear();
    pKurama.clear();
    qKurama.clear();
    m2.clear();
    thetaKurama.clear();
    xtgtKurama.clear();
    ytgtKurama.clear();
    utgtKurama.clear();
    vtgtKurama.clear();
    pathwcKurama.clear();
    xin.clear();
    yin.clear();
    zin.clear();
    pxin.clear();
    pyin.clear();
    pzin.clear();
    xout.clear();
    yout.clear();
    zout.clear();
    pxout.clear();
    pyout.clear();
    pzout.clear();

    nKm = 0;
    nKp = 0;
    nKK = 0;
    inside.clear();
    vtx.clear();
    vty.clear();
    vtz.clear();
    closeDist.clear();
    MissMass.clear();
    MissMassCorr.clear();
    MissMassCorrDE.clear();
    pOrg.clear();
    pCalc.clear();
    pCorr.clear();
    pCorrDE.clear();
    xb.clear();
    yb.clear();
    ub.clear();
    vb.clear();
    xs.clear();
    ys.clear();
    us.clear();
    vs.clear();
    Kflag.clear();

    isgoodTPCK18.clear();
    chisqrTPCK18.clear();
    qTPCK18.clear();
    pTPCK18.clear();
    xtgtTPCK18.clear();
    ytgtTPCK18.clear();
    utgtTPCK18.clear();
    vtgtTPCK18.clear();
    thetaTPCK18.clear();
    lhtofTPCK18.clear();
    xhtofTPCK18.clear();
    yhtofTPCK18.clear();
    lvpTPCK18.clear();
    xvpTPCK18.clear();
    yvpTPCK18.clear();

    isgoodTPCKurama.clear();
    kflagTPCKurama.clear();
    chisqrTPCKurama.clear();
    pTPCKurama.clear();
    qTPCKurama.clear();
    m2TPCKurama.clear();
    xtgtTPCKurama.clear();
    ytgtTPCKurama.clear();
    utgtTPCKurama.clear();
    vtgtTPCKurama.clear();
    thetaTPCKurama.clear();
    pathTPCKurama.clear();
    lhtofTPCKurama.clear();
    xhtofTPCKurama.clear();
    yhtofTPCKurama.clear();
    lvpTPCKurama.clear();
    xvpTPCKurama.clear();
    yvpTPCKurama.clear();

    isgoodTPC.clear();
    insideTPC.clear();
    vtxTPC.clear();
    vtyTPC.clear();
    vtzTPC.clear();
    closeDistTPC.clear();
    MissMassTPC.clear();
    MissMassCorrTPC.clear();
    MissMassCorrDETPC.clear();
    pOrgTPC.clear();
    pCalcTPC.clear();
    pCorrTPC.clear();
    pCorrDETPC.clear();
    thetaTPC.clear();
    thetaCMTPC.clear();
    costCMTPC.clear();
    xbTPC.clear();
    ybTPC.clear();
    ubTPC.clear();
    vbTPC.clear();
    xsTPC.clear();
    ysTPC.clear();
    usTPC.clear();
    vsTPC.clear();

    BE.clear();
    BETPC.clear();
    BE_LL.clear();
    BETPC_LL.clear();
    km_mom_x.clear();
    km_mom_y.clear();
    km_mom_z.clear();
    kp_mom_x.clear();
    kp_mom_y.clear();
    kp_mom_z.clear();

    nclTpc = 0;
    cluster_x.clear();
    cluster_y.clear();
    cluster_z.clear();
    cluster_de.clear();
    cluster_size.clear();
    cluster_layer.clear();
    cluster_mrow.clear();
    cluster_de_center.clear();
    cluster_houghflag.clear();
    cluster_G4tid.clear();
    cluster_G4pid.clear();

    remain_nclTpc = 0;
    remain_cluster_x.clear();
    remain_cluster_y.clear();
    remain_cluster_z.clear();
    remain_cluster_de.clear();
    remain_cluster_size.clear();
    remain_cluster_layer.clear();
    remain_cluster_houghflag.clear();
    remain_cluster_G4tid.clear();
    remain_cluster_G4pid.clear();

    ntTpc = 0;
    nhtrack.clear();
    trackid.clear();
    isXi.clear();
    isBeam.clear();
    isKurama.clear();
    isK18.clear();
    isAccidental.clear();
    isInTarget.clear();
    isMultiloop.clear();
    charge.clear();
    pid.clear();

    chisqr.clear();
    pval.clear();
    purity.clear();
    efficiency.clear();
    G4tid.clear();
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
    mom_x.clear();
    mom_y.clear();
    mom_z.clear();
    residual.clear();
    residual_x.clear();
    residual_y.clear();
    residual_z.clear();
    resolution_x.clear();
    resolution_y.clear();
    resolution_z.clear();
    helix_t.clear();
    pathhit.clear();
    alpha.clear();
    track_cluster_de.clear();
    track_cluster_size.clear();
    track_cluster_mrow.clear();

    chargeIndistinguishable.clear();
    chisqr_inverted.clear();
    pval_inverted.clear();
    helix_cx_inverted.clear();
    helix_cy_inverted.clear();
    helix_z0_inverted.clear();
    helix_r_inverted.clear();
    helix_dz_inverted.clear();
    mom0_inverted.clear();
    pid_inverted.clear();

    nvtxTpc = 0;
    vtx_x.clear();
    vtx_y.clear();
    vtx_z.clear();
    vtx_dist.clear();
    vtx_angle.clear();
    vtxid.clear();
    vtxmom_theta.clear();
    vtxpos_x.clear();
    vtxpos_y.clear();
    vtxpos_z.clear();
    vtxmom_x.clear();
    vtxmom_y.clear();
    vtxmom_z.clear();

    nvtxTpcClustered = 0;
    clusteredVtx_x.clear();
    clusteredVtx_y.clear();
    clusteredVtx_z.clear();
    clusteredVtxid.clear();

    ncombiLreconfailed = 0;
    pidLreconfailed.clear();
    piidLreconfailed.clear();
    LdecayvtxLreconfailed_x.clear();
    LdecayvtxLreconfailed_y.clear();
    LdecayvtxLreconfailed_z.clear();
    LmassLreconfailed.clear();
    LmomLreconfailed.clear();
    LmomLreconfailed_x.clear();
    LmomLreconfailed_y.clear();
    LmomLreconfailed_z.clear();
    pmomLreconfailed.clear();
    pmomLreconfailed_x.clear();
    pmomLreconfailed_y.clear();
    pmomLreconfailed_z.clear();
    pimomLreconfailed.clear();
    pimomLreconfailed_x.clear();
    pimomLreconfailed_y.clear();
    pimomLreconfailed_z.clear();
    ppidistLreconfailed.clear();

    ncombiPipair = 0;
    pipidPipair.clear();
    pimidPipair.clear();
    pipmomPipair.clear();
    pipmomPipair_x.clear();
    pipmomPipair_y.clear();
    pipmomPipair_z.clear();
    pimmomPipair.clear();
    pimmomPipair_x.clear();
    pimmomPipair_y.clear();
    pimmomPipair_z.clear();
    momPipair.clear();
    momPipair_x.clear();
    momPipair_y.clear();
    momPipair_z.clear();
    reconLmassPipair.clear();
    pipidistPipair.clear();

    isLambda.clear();
    ncombiLambda.clear();
    distLambda.clear();
    angleLambda.clear();
    bestmassLambda.clear();
    massLambda.clear();
    vtxLambda_x.clear();
    vtxLambda_y.clear();
    vtxLambda_z.clear();
    momLambda.clear();
    momLambda_x.clear();
    momLambda_y.clear();
    momLambda_z.clear();
    decaysidLambda.clear();
    decaysmomLambda.clear();
    decaysmomLambda_x.clear();
    decaysmomLambda_y.clear();
    decaysmomLambda_z.clear();

    GFprodvtx_x_ll = qnan;
    GFprodvtx_y_ll = qnan;
    GFprodvtx_z_ll = qnan;
    GFprodvtx_x_l1 = qnan;
    GFprodvtx_y_l1 = qnan;
    GFprodvtx_z_l1 = qnan;
    GFprodvtx_x_l2 = qnan;
    GFprodvtx_y_l2 = qnan;
    GFprodvtx_z_l2 = qnan;
    GFprodvtx_x_l = qnan;
    GFprodvtx_y_l = qnan;
    GFprodvtx_z_l = qnan;

    GFntTpc = 0;
    GFcharge.clear();
    GFchisqr.clear();
    GFtof.clear();
    GFtracklen.clear();
    GFpval.clear();
    GFpdgcode.clear();
    GFnhtrack.clear();

    GFntdecays = 0;
    GFdecays_htofid.clear();
    GFdecays_tracklen.clear();
    GFdecays_tof.clear();
    GFdecays_mass2.clear();
    GFdecays_mom.clear();
    GFdecays_mom_x.clear();
    GFdecays_mom_y.clear();
    GFdecays_mom_z.clear();
    GFdecays_CMmom.clear();
    GFdecays_CMmom_x.clear();
    GFdecays_CMmom_y.clear();
    GFdecays_CMmom_z.clear();
    GFmomloss.clear();
    GFeloss.clear();

    xiflag = false;
    ximass = qnan;
    xidecayvtx_x = qnan;
    xidecayvtx_y = qnan;
    xidecayvtx_z = qnan;
    ximom = qnan;
    ximom_x = qnan;
    ximom_y = qnan;
    ximom_z = qnan;
    lpi_dist = qnan;

    lmass = qnan;
    l_intarget = 0;
    ldecayvtx_x = qnan;
    ldecayvtx_y = qnan;
    ldecayvtx_z = qnan;
    lmom = qnan;
    lmom_x = qnan;
    lmom_y = qnan;
    lmom_z = qnan;
    ppi_dist = qnan;
    ltarget_dist = qnan;
    ltargetvtx_x = qnan;
    ltargetvtx_y = qnan;
    ltargetvtx_z = qnan;

    lmass_vtx = qnan;
    ldecayvtx_x_vtx = qnan;
    ldecayvtx_y_vtx = qnan;
    ldecayvtx_z_vtx = qnan;
    lmom_vtx = qnan;
    lmom_x_vtx = qnan;
    lmom_y_vtx = qnan;
    lmom_z_vtx = qnan;
    ppi_dist_vtx = qnan;

    xitargetvtx_x = qnan;
    xitargetvtx_y = qnan;
    xitargetvtx_z = qnan;
    xitargetmom = qnan;
    xitargetmom_x = qnan;
    xitargetmom_y = qnan;
    xitargetmom_z = qnan;
    xitarget_dist = qnan;

    GFximass = qnan;
    GFxidecayvtx_x = qnan;
    GFxidecayvtx_y = qnan;
    GFxidecayvtx_z = qnan;
    GFximom = qnan;
    GFximom_x = qnan;
    GFximom_y = qnan;
    GFximom_z = qnan;
    GFlpi_dist = qnan;

    GFlmass = qnan;
    GFldecayvtx_x = qnan;
    GFldecayvtx_y = qnan;
    GFldecayvtx_z = qnan;
    GFlmom = qnan;
    GFlmom_x = qnan;
    GFlmom_y = qnan;
    GFlmom_z = qnan;
    GFltracklen = qnan;
    GFltof = qnan;
    GFppi_dist = qnan;
    GFltarget_dist = qnan;
    GFltargetvtx_x = qnan;
    GFltargetvtx_y = qnan;
    GFltargetvtx_z = qnan;
    GFltargetcenter_dist = qnan;
    GFltargetcenter_x = qnan;
    GFltargetcenter_y = qnan;
    GFltargetcenter_z = qnan;
    GFlprodvtx_x = qnan;
    GFlprodvtx_y = qnan;
    GFlprodvtx_z = qnan;
    GFlprodvtx_dist = qnan;
    GFltracklen = qnan;
    GFltof = qnan;

    GFxitargetvtx_x = qnan;
    GFxitargetvtx_y = qnan;
    GFxitargetvtx_z = qnan;
    GFxitargetmom = qnan;
    GFxitargetmom_x = qnan;
    GFxitargetmom_y = qnan;
    GFxitargetmom_z = qnan;
    GFxitarget_dist = qnan;

    GFxitargetcenter_x = qnan;
    GFxitargetcenter_y = qnan;
    GFxitargetcenter_z = qnan;
    GFxitargetcentermom = qnan;
    GFxitargetcentermom_x = qnan;
    GFxitargetcentermom_y = qnan;
    GFxitargetcentermom_z = qnan;
    GFxitargetcenter_dist = qnan;

    GFxikkvtx_x = qnan;
    GFxikkvtx_y = qnan;
    GFxikkvtx_z = qnan;
    GFxikkmom = qnan;
    GFxikkmom_x = qnan;
    GFxikkmom_y = qnan;
    GFxikkmom_z = qnan;
    GFxikkvtx_dist = qnan;

    GFxiprodvtx_x = qnan;
    GFxiprodvtx_y = qnan;
    GFxiprodvtx_z = qnan;
    GFxiprodmom = qnan;
    GFxiprodmom_x = qnan;
    GFxiprodmom_y = qnan;
    GFxiprodmom_z = qnan;
    GFxiprodvtx_dist = qnan;
    GFxitracklen = qnan;
    GFxitof = qnan;
    GFximomloss = qnan;
    GFxiexcitation = qnan;

    GFprodvtx_x_kkxi = qnan;
    GFprodvtx_y_kkxi = qnan;
    GFprodvtx_z_kkxi = qnan;

    GFximass_pmomcorr = qnan;
    GFlmass_pmomcorr = qnan;
    GFximom_pmomcorr = qnan;
    GFlmom_pmomcorr = qnan;
    GFdecays_mom_pmomcorr.clear();
    GFdecays_mom_x_pmomcorr.clear();
    GFdecays_mom_y_pmomcorr.clear();
    GFdecays_mom_z_pmomcorr.clear();
    GFpval_pmomcorr.clear();
    GFchisqr_pmomcorr.clear();

    lphiflag = false;
    phimass  = qnan;
    phidecayvtx_x  = qnan;
    phidecayvtx_y = qnan;
    phidecayvtx_z = qnan;
    phimom = qnan;
    phimom_x = qnan;
    phimom_y = qnan;
    phimom_z = qnan;
    kk_dist = qnan;
    phi_km_mass2 = qnan;
    GFphimass = qnan;
    GFphidecayvtx_x = qnan;
    GFphidecayvtx_y = qnan;
    GFphidecayvtx_z = qnan;
    GFphimom = qnan;
    GFphimom_x = qnan;
    GFphimom_y = qnan;
    GFphimom_z = qnan;
    GFkk_dist = qnan;
    GFphiprodvtx_dist = qnan;

    phidecays_id.clear();
    phidecays_mom.clear();
    phidecays_mom_x.clear();
    phidecays_mom_y.clear();
    phidecays_mom_z.clear();
    GFphidecays_mom.clear();
    GFphidecays_mom_x.clear();
    GFphidecays_mom_y.clear();
    GFphidecays_mom_z.clear();

    emptyflag = false;
    pimflag = false;
    lpiflag = false;
    lflag = false;
    llflag = false;
    pipiflag = false;

    lmass1 = qnan;
    ltarget_dist1 = qnan;
    ltargetvtx_x1 = qnan;
    ltargetvtx_y1 = qnan;
    ltargetvtx_z1 = qnan;
    ldecayvtx_x1 = qnan;
    ldecayvtx_y1 = qnan;
    ldecayvtx_z1 = qnan;
    lmom1 = qnan;
    lmom_x1 = qnan;
    lmom_y1 = qnan;
    lmom_z1 = qnan;
    ppi_dist1 = qnan;
    lmass2 = qnan;
    ltarget_dist2 = qnan;
    ltargetvtx_x2 = qnan;
    ltargetvtx_y2 = qnan;
    ltargetvtx_z2 = qnan;

    ldecayvtx_x2 = qnan;
    ldecayvtx_y2 = qnan;
    ldecayvtx_z2 = qnan;
    lmom2 = qnan;
    lmom_x2 = qnan;
    lmom_y2 = qnan;
    lmom_z2 = qnan;
    ppi_dist2 = qnan;

    GFllexcitation = qnan;
    GFlmass1 = qnan;
    GFldecayvtx_x1 = qnan;
    GFldecayvtx_y1 = qnan;
    GFldecayvtx_z1 = qnan;
    GFlmom1 = qnan;
    GFlmom_x1 = qnan;
    GFlmom_y1 = qnan;
    GFlmom_z1 = qnan;
    GFppi_dist1 = qnan;
    GFltarget_dist1 = qnan;
    GFltargetvtx_x1 = qnan;
    GFltargetvtx_y1 = qnan;
    GFltargetvtx_z1 = qnan;
    GFltargetcenter_dist1 = qnan;
    GFltargetcenter_x1 = qnan;
    GFltargetcenter_y1 = qnan;
    GFltargetcenter_z1 = qnan;
    GFlprodvtx_x1 = qnan;
    GFlprodvtx_y1 = qnan;
    GFlprodvtx_z1 = qnan;
    GFlprodvtx_dist1 = qnan;
    GFltracklen1 = qnan;
    GFltof1 = qnan;

    GFlmass2 = qnan;
    GFldecayvtx_x2 = qnan;
    GFldecayvtx_y2 = qnan;
    GFldecayvtx_z2 = qnan;
    GFlmom2 = qnan;
    GFlmom_x2 = qnan;
    GFlmom_y2 = qnan;
    GFlmom_z2 = qnan;
    GFppi_dist2 = qnan;
    GFltarget_dist2 = qnan;
    GFltargetvtx_x2 = qnan;
    GFltargetvtx_y2 = qnan;
    GFltargetvtx_z2 = qnan;
    GFltargetcenter_dist2 = qnan;
    GFltargetcenter_x2 = qnan;
    GFltargetcenter_y2 = qnan;
    GFltargetcenter_z2 = qnan;
    GFlprodvtx_x2 = qnan;
    GFlprodvtx_y2 = qnan;
    GFlprodvtx_z2 = qnan;
    GFlprodvtx_dist2 = qnan;
    GFltracklen2 = qnan;
    GFltof2 = qnan;

    llvtx_x = qnan;
    llvtx_y = qnan;
    llvtx_z = qnan;
    lldist = qnan;
    GFllvtx_x = qnan;
    GFllvtx_y = qnan;
    GFllvtx_z = qnan;
    GFlldist = qnan;

    KFllexcitation = qnan;
    KFlmom1 = qnan;
    KFlmom_x1 = qnan;
    KFlmom_y1 = qnan;
    KFlmom_z1 = qnan;
    KFlmom2 = qnan;
    KFlmom_x2 = qnan;
    KFlmom_y2 = qnan;
    KFlmom_z2 = qnan;
    KFlchisqr1 = qnan;
    KFlchisqr2 = qnan;

    KFprodvtx_chisqr_ll = qnan;
    KFprodvtx_x_ll = qnan;
    KFprodvtx_y_ll = qnan;
    KFprodvtx_z_ll = qnan;
    KFprodvtx_x_l1 = qnan;
    KFprodvtx_y_l1 = qnan;
    KFprodvtx_z_l1 = qnan;
    KFprodvtx_x_l2 = qnan;
    KFprodvtx_y_l2 = qnan;
    KFprodvtx_z_l2 = qnan;
    KFprodvtx_x_l = qnan;
    KFprodvtx_y_l = qnan;
    KFprodvtx_z_l = qnan;

    KFllvtx_x = qnan;
    KFllvtx_y = qnan;
    KFllvtx_z = qnan;
    KFlldist = qnan;

    KFlprodvtx_x1 = qnan;
    KFlprodvtx_y1 = qnan;
    KFlprodvtx_z1 = qnan;
    KFlprodvtx_dist1 = qnan;
    KFltracklen1 = qnan;
    KFltof1 = qnan;

    KFlprodvtx_x2 = qnan;
    KFlprodvtx_y2 = qnan;
    KFlprodvtx_z2 = qnan;
    KFlprodvtx_dist2 = qnan;
    KFltracklen2 = qnan;
    KFltof2 = qnan;

    KFlprodvtx_x = qnan;
    KFlprodvtx_y = qnan;
    KFlprodvtx_z = qnan;
    KFlprodvtx_dist = qnan;
    KFltracklen = qnan;
    KFltof = qnan;

    decays_id.clear();
    decays_purity.clear();
    decays_efficiency.clear();
    decays_G4tid.clear();
    decays_mom.clear();
    decays_mom_x.clear();
    decays_mom_y.clear();
    decays_mom_z.clear();
    decays_CMmom.clear();
    decays_CMmom_x.clear();
    decays_CMmom_y.clear();
    decays_CMmom_z.clear();

    residual_multi = 0;
    pim_multi = 0;
    pip_multi = 0;
    p_multi = 0;
    ppip_multi = 0;
    accident_multi = 0;

    residual_id.clear();
    residual_dist2tgt.clear();
    residual_mass2.clear();
    residual_mom.clear();
    residual_mom_x.clear();
    residual_mom_y.clear();
    residual_mom_z.clear();
    residual_charge.clear();

    KFlmom0 = qnan;
    KFlmom_x0 = qnan;
    KFlmom_y0 = qnan;
    KFlmom_z0 = qnan;
    KFlmom = qnan;
    KFlmom_x = qnan;
    KFlmom_y = qnan;
    KFlmom_z = qnan;
    KFlchisqr = qnan;
    KFlpval = qnan;
    KFlpi_dist = qnan;
    KFximom = qnan;
    KFximom_x = qnan;
    KFximom_y = qnan;
    KFximom_z = qnan;
    KFxichisqr = qnan;
    KFxipval = qnan;
    KFximass = qnan;
    KFxidecayvtx_x = qnan;
    KFxidecayvtx_y = qnan;
    KFxidecayvtx_z = qnan;
    KFdecays_mom.clear();
    KFdecays_mom_x.clear();
    KFdecays_mom_y.clear();
    KFdecays_mom_z.clear();
    KFdecays_CMmom.clear();
    KFdecays_CMmom_x.clear();
    KFdecays_CMmom_y.clear();
    KFdecays_CMmom_z.clear();
    KFlpull.clear();
    KFxipull.clear();

    KFprodvtx_chisqr_kkxi = qnan;
    KFprodvtx_x_kkxi = qnan;
    KFprodvtx_y_kkxi = qnan;
    KFprodvtx_z_kkxi = qnan;
    KFprodvtx_x_kpxi = qnan;
    KFprodvtx_y_kpxi = qnan;
    KFprodvtx_z_kpxi = qnan;

    KFxi_kkvtx_x = qnan;
    KFxi_kkvtx_y = qnan;
    KFxi_kkvtx_z = qnan;
    KFxi_kkvtx_mom = qnan;
    KFxi_kkvtx_mom_x = qnan;
    KFxi_kkvtx_mom_y = qnan;
    KFxi_kkvtx_mom_z = qnan;
    KFxi_kkvtx_dist = qnan;

    KFxiprodvtx_x = qnan;
    KFxiprodvtx_y = qnan;
    KFxiprodvtx_z = qnan;
    KFxiprodmom = qnan;
    KFxiprodmom_x = qnan;
    KFxiprodmom_y = qnan;
    KFxiprodmom_z = qnan;
    KFxiprodvtx_dist = qnan;
    KFxitracklen = qnan;
    KFxitof = qnan;
    KFximomloss = qnan;
    KFxiexcitation = qnan;

    KFxi_kpxiprodvtx_x = qnan;
    KFxi_kpxiprodvtx_y = qnan;
    KFxi_kpxiprodvtx_z = qnan;
    KFxi_kpxiprodmom = qnan;
    KFxi_kpxiprodmom_x = qnan;
    KFxi_kpxiprodmom_y = qnan;
    KFxi_kpxiprodmom_z = qnan;
    KFxi_kpxiprodvtx_dist = qnan;

    KFxitargetvtx_x = qnan;
    KFxitargetvtx_y = qnan;
    KFxitargetvtx_z = qnan;
    KFxitargetmom = qnan;
    KFxitargetmom_x = qnan;
    KFxitargetmom_y = qnan;
    KFxitargetmom_z = qnan;
    KFxitarget_dist = qnan;

    KFxitargetcenter_x = qnan;
    KFxitargetcenter_y = qnan;
    KFxitargetcenter_z = qnan;
    KFxitargetcentermom = qnan;
    KFxitargetcentermom_x = qnan;
    KFxitargetcentermom_y = qnan;
    KFxitargetcentermom_z = qnan;
    KFxitargetcenter_dist = qnan;

//G4 Initialization
    G4kmid = qnan,G4kmtid = qnan;
    G4kmvtx_x = qnan,G4kmvtx_y = qnan,G4kmvtx_z = qnan;
    G4kmmom = qnan,G4kmmom_x = qnan,G4kmmom_y = qnan,G4kmmom_z = qnan;

    G4kpid = qnan,G4kptid = qnan;
    G4kpvtx_x = qnan,G4kpvtx_y = qnan,G4kpvtx_z = qnan;
    G4kpmom = qnan,G4kpmom_x = qnan,G4kpmom_y = qnan,G4kpmom_z = qnan;

    G4l1id =-1;
    G4l1vtx_x = qnan,G4l1vtx_y = qnan,G4l1vtx_z = qnan;
    G4l1mom = qnan,G4l1mom_x = qnan,G4l1mom_y = qnan,G4l1mom_z = qnan;
    l1vtx_x = qnan,l1vtx_y = qnan,l1vtx_z = qnan;

    G4l2id =-1;
    G4l2vtx_x = qnan,G4l2vtx_y = qnan,G4l2vtx_z = qnan;
    G4l2mom = qnan,G4l2mom_x = qnan,G4l2mom_y = qnan,G4l2mom_z = qnan;
    l2vtx_x = qnan,l2vtx_y = qnan,l2vtx_z = qnan;

    G4llmass = qnan;

    G4p1id = -1,G4p1tid = -1,G4p1nh = -1,G4p1tnh = -1;
    G4p1vtx_x = qnan,G4p1vtx_y = qnan,G4p1vtx_z = qnan;
    G4p1mom = qnan,G4p1mom_x = qnan,G4p1mom_y = qnan,G4p1mom_z = qnan;

    G4p2id = -1,G4p2tid = -1,G4p2nh = -1,G4p2tnh = -1;
    G4p2vtx_x = qnan,G4p2vtx_y = qnan,G4p2vtx_z = qnan;
    G4p2mom = qnan,G4p2mom_x = qnan,G4p2mom_y = qnan,G4p2mom_z = qnan;

    G4pi1id = -1,G4pi1tid = -1,G4pi1nh = -1,G4pi1tnh = -1;
    G4pi1vtx_x = qnan,G4pi1vtx_y = qnan,G4pi1vtx_z = qnan;
    G4pi1mom = qnan,G4pi1mom_x = qnan,G4pi1mom_y = qnan,G4pi1mom_z = qnan;

    G4pi2id = -1,G4pi2tid = -1,G4pi2nh = -1,G4pi2tnh = -1;
    G4pi2vtx_x = qnan,G4pi2vtx_y = qnan,G4pi2vtx_z = qnan;
    G4pi2mom = qnan,G4pi2mom_x = qnan,G4pi2mom_y = qnan,G4pi2mom_z = qnan;

    l1good = false;l2good = false; llswap = false; piswap = false;
    p1_tracked = false;p2_tracked = false;
    pi1_tracked = false;pi2_tracked = false;
    p1t_mom0 = qnan;p2t_mom0 = qnan;
    pi1t_mom0 = qnan;pi2t_mom0 = qnan;
  }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* evnum;

  TTreeReaderValue<Int_t>* nhHtof;
  TTreeReaderValue<std::vector<Double_t>>* HtofSeg;
  TTreeReaderValue<std::vector<Double_t>>* tHtof;
  TTreeReaderValue<std::vector<Double_t>>* dtHtof;
  TTreeReaderValue<std::vector<Double_t>>* deHtof;
  TTreeReaderValue<std::vector<Double_t>>* posHtof;
  TTreeReaderValue<std::vector<Int_t>>* G4tidHtof;

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

  TTreeReaderValue<Int_t>* nclTpc; // Number of clusters
  TTreeReaderValue<Int_t>* remain_nclTpc; // Number of clusters without tracks
  TTreeReaderValue<std::vector<Double_t>>* cluster_x;
  TTreeReaderValue<std::vector<Double_t>>* cluster_y;
  TTreeReaderValue<std::vector<Double_t>>* cluster_z;
  TTreeReaderValue<std::vector<Double_t>>* cluster_de;
  TTreeReaderValue<std::vector<Int_t>>* cluster_layer;
  TTreeReaderValue<std::vector<Double_t>>* cluster_mrow;
  TTreeReaderValue<std::vector<Int_t>>* cluster_houghflag;
//  TTreeReaderValue<std::vector<Int_t>>* cluster_G4tid;

  TTreeReaderValue<Int_t>* ntTpc; // Number of Tracks
  TTreeReaderValue<std::vector<Int_t>>* nhtrack; // Number of Hits (in 1 tracks)
  TTreeReaderValue<std::vector<Int_t>>* trackid; //for Kurama K1.8 tracks
  TTreeReaderValue<std::vector<Int_t>>* isXi;
  TTreeReaderValue<std::vector<Int_t>>* isBeam;
  TTreeReaderValue<std::vector<Int_t>>* isKurama;
  TTreeReaderValue<std::vector<Int_t>>* isK18;
  TTreeReaderValue<std::vector<Int_t>>* isAccidental;
  TTreeReaderValue<std::vector<Int_t>>* isMultiloop;
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
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* pathhit;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* alpha;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de;
//  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_size;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_mrow;

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
    GFTrackCont.Clear();
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
  static const auto KaonMass    = pdg::KaonMass();
  static const auto ProtonMass  = pdg::ProtonMass();
  static const auto LambdaMass  = pdg::LambdaMass();
  static const auto XiMinusMass = pdg::XiMinusMass();
  static const auto m12C = 11.174864;
  static const auto m11B = 10.252548;
  static const auto m10Be = 9.325504;
  static const auto me = 0.001*0.5109989461;
  static const int XiMinusPdgCode = 3312;
  Double_t pdgmass[3] = {ProtonMass, KaonMass, PionMass};

  TVector3 tgtpos(0, 0, tpc::ZTarget);

  //if( ievent%100000==0 ){
  //if( ievent%1000==0 ){
  //if( ievent%100==0 ){
  //if( ievent%10==0 ){
  if( ievent%1==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);
  event.evnum = **src.evnum;


  vector<TVector3> G4Hits;
  vector<TVector3> G4Moms;
  for(int ih=0;ih<src.nhittpc;++ih){
    TVector3 G4Hit(src.xtpc[ih],src.ytpc[ih],src.ztpc[ih]);
    TVector3 G4Mom(src.pxtpc[ih],src.pytpc[ih],src.pztpc[ih]);
    G4Hits.push_back(G4Hit);
    G4Moms.push_back(G4Mom);
  }
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
    if(pid == 3122 && it == 1){
      event.G4l1id = it;
      event.G4l1vtx_x = vert_x;
      event.G4l1vtx_y = vert_y;
      event.G4l1vtx_z = vert_z;
      event.G4l1mom = mom;
      event.G4l1mom_x = mom_x;
      event.G4l1mom_y = mom_y;
      event.G4l1mom_z = mom_z;
    }
    if(pid == 2212 && parent == 1){
      event.G4p1id = it;
      event.G4p1vtx_x = vert_x;
      event.G4p1vtx_y = vert_y;
      event.G4p1vtx_z = vert_z;
      event.G4p1mom = mom;
      event.G4p1mom_x = mom_x;
      event.G4p1mom_y = mom_y;
      event.G4p1mom_z = mom_z;
    }
    if(pid == -211 && parent == 1){
      event.G4pi1id = it;
      event.G4pi1vtx_x = vert_x;
      event.G4pi1vtx_y = vert_y;
      event.G4pi1vtx_z = vert_z;
      event.G4pi1mom = mom;
      event.G4pi1mom_x = mom_x;
      event.G4pi1mom_y = mom_y;
      event.G4pi1mom_z = mom_z;
    }
    if(pid == 2212 && it == 2){
      event.G4p2id = it;
      event.G4p2vtx_x = vert_x;
      event.G4p2vtx_y = vert_y;
      event.G4p2vtx_z = vert_z;
      event.G4p2mom = mom;
      event.G4p2mom_x = mom_x;
      event.G4p2mom_y = mom_y;
      event.G4p2mom_z = mom_z;
    }
    if(pid == -211 && it == 3){
      event.G4pi2id = it;
      event.G4pi2vtx_x = vert_x;
      event.G4pi2vtx_y = vert_y;
      event.G4pi2vtx_z = vert_z;
      event.G4pi2mom = mom;
      event.G4pi2mom_x = mom_x;
      event.G4pi2mom_y = mom_y;
      event.G4pi2mom_z = mom_z;
    }
  }
  vector<int> G4TrackID;
  vector<int> PureHits;
  event.G4p1nh = CountHits(event.G4p1id,src.xtpc,src.ztpc,src.ititpc,src.nhittpc);
  event.G4p2nh = CountHits(event.G4p2id,src.xtpc,src.ztpc,src.ititpc,src.nhittpc);
  event.G4pi1nh = CountHits(event.G4pi1id,src.xtpc,src.ztpc,src.ititpc,src.nhittpc);
  event.G4pi2nh = CountHits(event.G4pi2id,src.xtpc,src.ztpc,src.ititpc,src.nhittpc);

  TLorentzVector G4LvL(event.G4l1mom_x,event.G4l1mom_y,event.G4l1mom_z,hypot(event.G4l1mom,LambdaMass));
  TLorentzVector G4Lvp2(event.G4p2mom_x,event.G4p2mom_y,event.G4p2mom_z,hypot(event.G4p2mom,ProtonMass));
  TLorentzVector G4Lvpi2(event.G4pi2mom_x,event.G4pi2mom_y,event.G4pi2mom_z,hypot(event.G4pi2mom,PionMass));
  TLorentzVector G4LvLL = G4LvL + G4Lvp2 + G4Lvpi2;
  event.G4llmass = G4LvLL.M();
  TVector3 kkvtxTPC(event.G4kpvtx_x, event.G4kpvtx_y, event.G4kpvtx_z);

  Double_t pKp = event.G4kpmom;
  Double_t uKp = event.G4kpmom_x/event.G4kpmom_z;
  Double_t vKp = event.G4kpmom_y/event.G4kpmom_z;
  Double_t ptKp = pKp/std::sqrt(1.+uKp*uKp+vKp*vKp);
  TVector3 kp_mom(ptKp*uKp, ptKp*vKp, ptKp);

  Double_t pKm = event.G4kmmom;
  Double_t uKm = event.G4kmmom_x/event.G4kmmom_z;
  Double_t vKm = event.G4kmmom_y/event.G4kmmom_z;
  Double_t ptKm = pKm/std::sqrt(1.+uKm*uKm+vKm*vKm);
  TVector3 km_mom(ptKm*uKm, ptKm*vKm, ptKm);
  //Double_t cost = km_mom*kp_mom/(pKm*pKp);
  //Double_t theta = TMath::ACos(cost)*TMath::RadToDeg();
  TLorentzVector LvKm(km_mom, TMath::Hypot(km_mom.Mag(), KaonMass));
  TLorentzVector LvKp(kp_mom, TMath::Hypot(kp_mom.Mag(), KaonMass));

  TLorentzVector LvC(0., 0., 0., m12C);
  TLorentzVector LvP(0., 0., 0., ProtonMass);
  TLorentzVector LvRproton = LvKm + LvP - LvKp;
  TLorentzVector LvRc = LvKm + LvC - LvKp;
  Double_t mm_12C = LvRc.Mag();
  Double_t binding_energy = m11B + XiMinusMass - mm_12C; //GeV/c2
  TLorentzVector LvMM = LvKm + LvP - LvKp;
  auto veolcityMM = LvMM.BoostVector();
  TLorentzVector LvKmCM = LvKm;
  TLorentzVector LvKpCM = LvKp;
  LvKmCM.Boost(-veolcityMM);
  LvKpCM.Boost(-veolcityMM);

  event.nKK = 1;
  event.Kflag.push_back(1);
  event.MissMass.push_back((LvKm+LvP-LvKp).M());
  event.MissMassCorr.push_back((LvKm+LvP-LvKp).M());
  event.MissMassCorrDE.push_back((LvKm+LvP-LvKp).M());
  event.pOrg.push_back(pKm);
  event.pCalc.push_back(pKm);
  event.pCorr.push_back(pKm);
  event.pCorrDE.push_back(pKm);
  event.ub.push_back(uKm);
  event.vb.push_back(vKm);
  event.us.push_back(uKp);
  event.vs.push_back(vKp);
  event.vtx.push_back(event.G4kpvtx_x);
  event.vty.push_back(event.G4kpvtx_y);
  event.vtz.push_back(event.G4kpvtx_z);
  event.closeDist.push_back(0);
  event.inside.push_back(1);

  event.pKurama.push_back(pKp);
  event.qKurama.push_back(1);
  event.chisqrKurama.push_back(1);
  event.thetaKurama.push_back(kp_mom.Theta()*TMath::RadToDeg());
  event.xtgtKurama.push_back(event.G4kpvtx_x);
  event.ytgtKurama.push_back(event.G4kpvtx_y);
  event.utgtKurama.push_back(uKp);
  event.vtgtKurama.push_back(vKp);
  event.xin.push_back(event.G4kpvtx_x);
  event.yin.push_back(event.G4kpvtx_y);
  event.zin.push_back(event.G4kpvtx_z);
  event.pxin.push_back(event.G4kmmom_x);
  event.pyin.push_back(event.G4kmmom_y);
  event.pzin.push_back(event.G4kmmom_z);
  event.xout.push_back(event.G4kpvtx_x);
  event.yout.push_back(event.G4kpvtx_y);
  event.zout.push_back(event.G4kpvtx_z);
  event.pxout.push_back(event.G4kpmom_x);
  event.pyout.push_back(event.G4kpmom_y);
  event.pzout.push_back(event.G4kpmom_z);

  event.ntK18 = 1;
  event.pK18.push_back(pKp);
  event.chisqrK18.push_back(1);
  event.xtgtK18.push_back(event.G4kpvtx_x);
  event.ytgtK18.push_back(event.G4kpvtx_y);
  event.utgtK18.push_back(uKp);
  event.vtgtK18.push_back(vKp);

  event.ntKurama = 1;
  event.isgoodTPCKurama.push_back(1);
  event.kflagTPCKurama.push_back(1);
  event.chisqrTPCKurama.push_back(1);
  event.pTPCKurama.push_back(pKp);
  event.qTPCKurama.push_back(1);
  event.m2TPCKurama.push_back((LvKp.M2()-KaonMass*KaonMass));
  event.xtgtTPCKurama.push_back(event.G4kpvtx_x);
  event.ytgtTPCKurama.push_back(event.G4kpvtx_y);
  event.utgtTPCKurama.push_back(uKp);
  event.vtgtTPCKurama.push_back(vKp);
  event.thetaTPCKurama.push_back(kp_mom.Theta()*TMath::RadToDeg());
  event.isgoodTPC.push_back(1);
  event.insideTPC.push_back(1);
  event.vtxTPC.push_back(event.G4kpvtx_x);
  event.vtyTPC.push_back(event.G4kpvtx_y);
  event.vtzTPC.push_back(event.G4kpvtx_z-tpc::ZTarget);
  event.closeDistTPC.push_back(0);
  event.MissMassTPC.push_back((LvKm+LvP-LvKp).M());
  event.MissMassCorrTPC.push_back((LvKm+LvP-LvKp).M());
  event.MissMassCorrDETPC.push_back((LvKm+LvP-LvKp).M());
  event.pOrgTPC.push_back(pKm);
  event.pCalcTPC.push_back(pKm);
  event.pCorrTPC.push_back(pKm);
  event.pCorrDETPC.push_back(pKm);
  event.thetaTPC.push_back(km_mom.Theta()*TMath::RadToDeg());
  event.thetaCMTPC.push_back(LvKmCM.Theta()*TMath::RadToDeg());
  event.costCMTPC.push_back(cos((LvKmCM.Vect()).Angle(LvKpCM.Vect())));
  event.xbTPC.push_back(event.G4kpvtx_x);
  event.ybTPC.push_back(event.G4kpvtx_y);
  event.ubTPC.push_back(uKm);
  event.vbTPC.push_back(vKm);
  event.xsTPC.push_back(event.G4kpvtx_x);
  event.ysTPC.push_back(event.G4kpvtx_y);
  event.usTPC.push_back(uKp);
  event.vsTPC.push_back(vKp);

  event.BE.resize(event.nKK);
  event.BETPC.resize(event.nKK);
  event.BE_LL.resize(event.nKK);
  event.BETPC_LL.resize(event.nKK);

  event.km_mom_x.resize(event.nKK);
  event.km_mom_y.resize(event.nKK);
  event.km_mom_z.resize(event.nKK);
  event.kp_mom_x.resize(event.nKK);
  event.kp_mom_y.resize(event.nKK);
  event.kp_mom_z.resize(event.nKK);

  event.BE[0] = 1000.*binding_energy; //MeV/c2
  Double_t binding_energy_LL = m10Be + 2.*LambdaMass - mm_12C; //GeV/c2
  event.BE_LL[0] = 1000.*binding_energy_LL; //MeV/c2

  event.km_mom_x[0] = km_mom.x();
  event.km_mom_y[0] = km_mom.y();
  event.km_mom_z[0] = km_mom.z();
  event.kp_mom_x[0] = kp_mom.x();
  event.kp_mom_y[0] = kp_mom.y();
  event.kp_mom_z[0] = kp_mom.z();

  HF1( 100, -event.BE[0]);

  TLorentzVector LvRcTPC;
  TVector3 km_unit = TVector3(event.ubTPC[0], event.vbTPC[0], 1.).Unit();
  TVector3 km_momTPC = km_unit*event.pK18[0];

  TVector3 kp_unit = TVector3(event.usTPC[0], event.vsTPC[0], 1.).Unit();
  TVector3 kp_momTPC = kp_unit*event.pCorrDETPC[0];
  //Double_t thetaTPC = event.thetaTPC[0];

  TLorentzVector LvKmTPC(km_momTPC, TMath::Hypot(km_momTPC.Mag(), KaonMass));
  TLorentzVector LvKpTPC(kp_momTPC, TMath::Hypot(kp_momTPC.Mag(), KaonMass));
  TLorentzVector LvCTPC(0., 0., 0., m12C);
  TLorentzVector LvPTPC(0., 0., 0., ProtonMass);
  TLorentzVector LvRprotonTPC = LvKmTPC + LvPTPC - LvKpTPC;
  LvRcTPC = LvKmTPC + LvCTPC - LvKpTPC;

  double mm_12CTPC = LvRcTPC.M();
  double binding_energyTPC = m11B + XiMinusMass - mm_12CTPC; //GeV/c2
  event.BETPC[0] = 1000.*binding_energyTPC; //MeV/c2
  Double_t binding_energyTPC_LL = m10Be + 2.*LambdaMass - mm_12CTPC; //GeV/c2
  event.BETPC_LL[0] = 1000.*binding_energyTPC_LL; //MeV/c2
  HF1( 102, -event.BETPC[0]);

  HF1( 1, event.status++ );

  event.nhHtof = **src.nhHtof;
  event.HtofSeg = **src.HtofSeg;
  event.tHtof = **src.tHtof;
  event.dtHtof = **src.dtHtof;
  event.deHtof = **src.deHtof;
  event.posHtof = **src.posHtof;
  event.G4tidHtof = **src.G4tidHtof;
  Int_t ntTpc = **src.ntTpc;
  if( ntTpc == 0 ) return true;
  HF1( 1, event.status++ );

  event.nclTpc = **src.nclTpc;
  event.remain_nclTpc = **src.remain_nclTpc;
#if SaveRawData
  event.cluster_x = **src.cluster_x;
  event.cluster_y = **src.cluster_y;
  event.cluster_z = **src.cluster_z;
  event.cluster_de = **src.cluster_de;
  event.cluster_layer = **src.cluster_layer;
  event.cluster_mrow = **src.cluster_mrow;
  event.cluster_houghflag = **src.cluster_houghflag;
//  event.cluster_G4tid = **src.cluster_G4tid;

  event.remain_cluster_x.resize(event.remain_nclTpc);
  event.remain_cluster_y.resize(event.remain_nclTpc);
  event.remain_cluster_z.resize(event.remain_nclTpc);
  event.remain_cluster_de.resize(event.remain_nclTpc);
//  event.remain_cluster_size.resize(event.remain_nclTpc);
  event.remain_cluster_layer.resize(event.remain_nclTpc);
//  event.remain_cluster_mrow.resize(event.remain_nclTpc);
  event.remain_cluster_houghflag.resize(event.remain_nclTpc);
//  event.remain_cluster_G4tid.resize(event.remain_nclTpc);
  Int_t icl_remain = 0;
  for( Int_t icl=0; icl<ntTpc; ++icl ){
    if(event.cluster_houghflag[icl]!=0) continue;
    event.remain_cluster_x[icl_remain] = event.cluster_x[icl];
    event.remain_cluster_y[icl_remain] = event.cluster_y[icl];
    event.remain_cluster_z[icl_remain] = event.cluster_z[icl];
    event.remain_cluster_de[icl_remain] = event.cluster_de[icl];
    event.remain_cluster_layer[icl_remain] = event.cluster_layer[icl];
//    event.remain_cluster_mrow[icl_remain] = event.cluster_mrow[icl];
    event.remain_cluster_houghflag[icl_remain] = event.cluster_houghflag[icl];
//    event.remain_cluster_G4tid.push_back(event.cluster_G4tid[icl]);
    icl_remain++;
  }
#endif

  event.ntTpc = ntTpc;
  event.nhtrack = **src.nhtrack;
  event.trackid = **src.trackid;
  event.isXi = **src.isXi;
  event.isBeam = **src.isBeam;
  event.isKurama = **src.isKurama;
  event.isK18 = **src.isK18;
  event.isAccidental = **src.isAccidental;
  event.isMultiloop = **src.isMultiloop;
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
  event.pathhit = **src.pathhit;
  event.alpha = **src.alpha;
  event.track_cluster_de = **src.track_cluster_de;
//  event.track_cluster_size = **src.track_cluster_size;
  event.track_cluster_mrow = **src.track_cluster_mrow;

  HF1( 1, event.status++ );
  if( ntTpc == 0 ) return true;
  HF1( 1, event.status++ );

  Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCTracksGeant4(**src.ntTpc,
      **src.charge, **src.nhtrack, **src.helix_cx,
      **src.helix_cy, **src.helix_z0, **src.helix_r,
      **src.helix_dz, **src.hitlayer, **src.track_cluster_mrow,
      **src.helix_t, **src.track_cluster_de, **src.resolution_x,
      **src.resolution_y, **src.resolution_z, **src.hitpos_x,
      **src.hitpos_y, **src.hitpos_z);
  HF1( 1, event.status++ );
  HypTPCTask& GFTrackCont = HypTPCTask::GetInstance();
  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix(it);
    if( !tp ) continue;
    std::vector<Int_t> pdgcode;
    Kinematics::HypTPCPID_PDGCode(event.charge[it], event.pid[it], pdgcode);
    GFTrackCont.AddHelixTrack(pdgcode, tp);
    vector<TVector3> TPCHit;
    for(int ih=0;ih<event.hitpos_x.at(it).size();++ih){
      TPCHit.push_back(TVector3(
            event.hitpos_x.at(it).at(ih),
            event.hitpos_y.at(it).at(ih),
            event.hitpos_z.at(it).at(ih)));
    }
    int nPureHit;
    int G4tid = TPCToG4TrackID(TPCHit,src.nhittpc,src.ititpc,src.xtpc,src.ytpc,src.ztpc,nPureHit);
    G4TrackID.push_back(G4tid);
    PureHits.push_back(nPureHit);
    if(G4tid<0) continue;
    if(G4tid == event.G4p1id){
      event.p1_tracked = true;
      event.p1t_mom0 = event.mom0[it];
    }
    if(G4tid == event.G4p2id){
      event.p2_tracked = true;
      event.p2t_mom0 = event.mom0[it];
    }
    if(G4tid == event.G4pi1id){
      event.pi1_tracked = true;
      event.pi1t_mom0 = event.mom0[it];
    }
    if(G4tid == event.G4pi2id){
      event.pi2_tracked = true;
      event.pi2t_mom0 = event.mom0[it];
    }

  }
  GFTrackCont.FitTracks();
  HF1( 1, event.status++ );

  Int_t GFntTpc = GFTrackCont.GetNTrack();
  if(GFntTpc!=event.ntTpc){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    return true;
  }

#if DebugDisp
  std::cout<<"0. Before particle reconstruction, checking each track is originated from the target or not."<<std::endl;
#endif

  std::vector<Int_t> target_accidental_id_container;
  std::vector<Int_t> target_k_id_container;
  std::vector<Double_t> target_k_mass2_container;
  std::vector<TVector3> target_k_mom_container;
  std::vector<Int_t> target_p_id_container, target_pip_id_container, target_pim_id_container, target_ppip_id_container;
  std::vector<Double_t> target_p_mass2_container, target_pip_mass2_container, target_pim_mass2_container, target_ppip_mass2_container;
  std::vector<TVector3> target_p_mom_container, target_pip_mom_container, target_pim_mom_container, target_ppip_mom_container;
  std::vector<Double_t> target_p_dist2tgt_container, target_pip_dist2tgt_container, target_pim_dist2tgt_container, target_ppip_dist2tgt_container;

  std::vector<Int_t> intarget_id_container;
  std::vector<Int_t> notintarget_id_container;
  event.isInTarget.resize(ntTpc,0);
  for(Int_t it=0;it<ntTpc;it++){ //1st searching
    if(event.isK18[it]==1) continue;
    if(event.isXi[it]==1) continue;
    if(event.isKurama[it]==1) continue;
    if(event.isAccidental[it]==1) continue;
    if(event.isBeam[it]==1) continue;

    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( it );
    Int_t nhit = event.nhtrack[it];
    TVector3 start(event.hitpos_x[it][0], event.hitpos_y[it][0], event.hitpos_z[it][0]);
    TVector3 end(event.hitpos_x[it][nhit-1], event.hitpos_y[it][nhit-1], event.hitpos_z[it][nhit-1]);
    if(nhit<=6 && TMath::Abs(start.x()) < 25. && TMath::Abs(end.x()) < 25. &&
       TMath::Abs(start.y()) < 25. && TMath::Abs(end.y()) < 25. &&
       start.z() < tpc::ZTarget && end.z() < tpc::ZTarget){
      event.isAccidental[it] = 1; //Accidental K-
      target_accidental_id_container.push_back(it); //Accidental beam on the target
    }
  }

  for(Int_t it=0;it<ntTpc;it++){ //1st searching
    if(event.isK18[it]==1) continue;
    if(event.isXi[it]==1) continue;
    if(event.isKurama[it]==1){
      if(event.isBeam[it]==1) target_accidental_id_container.push_back(it); //Accidental beam on the target
      continue;
    }
    if(event.isAccidental[it]==1) continue;
    if(!GFTrackCont.IsInsideTarget(it)){
      notintarget_id_container.push_back(it);
      continue;
    }
    else{
      event.isInTarget[it] = 1;
    }
    if(event.isBeam[it]==1){
      target_accidental_id_container.push_back(it); //Accidental beam on the target
      continue;
    }

    intarget_id_container.push_back(it);
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( it );
    Int_t nhit = event.nhtrack[it];
    TVector3 start(event.hitpos_x[it][0], event.hitpos_y[it][0], event.hitpos_z[it][0]);
    TVector3 end(event.hitpos_x[it][nhit-1], event.hitpos_y[it][nhit-1], event.hitpos_z[it][nhit-1]);
    if((event.pid[it]&2)==2 && event.charge[it]==-1){ //k-
      Bool_t km_pid = Kinematics::HypTPCdEdxPID_IsKaonTemp(event.dEdx[it], event.mom0[it]);
      if(km_pid){
	Int_t repid = 1;
	if(GFTrackCont.TrackCheck(it, repid)){
	  Int_t htofhitid_k; Double_t tracklen_k;
	  Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	  Bool_t htofextrapolation_k = GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
									event.HtofSeg, event.posHtof,
									htofhitid_k, tof, tracklen_k,
									pos, track2tgt_dist);
	  Double_t mass2 = qnan;
	  TVector3 mom = GFTrackCont.GetMom(it, 0, repid);
	  if(htofextrapolation_k) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_k, event.tHtof[htofhitid_k]);
	  target_k_id_container.push_back(it);
	  target_k_mass2_container.push_back(mass2);
	  target_k_mom_container.push_back(mom);
	}
      }
    }

    if((event.pid[it]&4)==4 && (event.pid[it]&1)!=1 && event.charge[it]==1){ //proton

      Double_t mass2 = qnan;
      TVector3 mom(qnan, qnan, qnan);
      target_p_id_container.push_back(it);
      target_p_dist2tgt_container.push_back(tp -> GetClosestDist());

      Int_t repid = 0;
      Int_t flag = 1;
      for(Int_t i=0;i<2;i++){
	Int_t temp = flag&event.pid[it];
	if(temp==flag) repid += 1;
	flag*=2;
      }
      if(!GFTrackCont.TrackCheck(it, repid)){
	target_p_mass2_container.push_back(mass2);
	target_p_mom_container.push_back(mom);
	continue;
      }

      Int_t htofhitid_p; Double_t tracklen_p;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation_p = GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
								    event.HtofSeg, event.posHtof,
								    htofhitid_p, tof, tracklen_p,
								    pos, track2tgt_dist);
      mom = GFTrackCont.GetMom(it, 0, repid);
      if(htofextrapolation_p) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_p, event.tHtof[htofhitid_p]);
      target_p_mass2_container.push_back(mass2);
      target_p_mom_container.push_back(mom);
    }
    else if((event.pid[it]&1)==1 && event.charge[it]==-1){ //pi-

      Double_t slope = event.helix_dz[it];
      Double_t helixmom = event.mom0[it];
      //if(TMath::Abs(slope)<0.05 && TMath::Abs(helixmom)>0.5){
      if(TMath::Abs(slope)<0.1 && TMath::Abs(helixmom)>0.5 &&
	 (start.x()-end.x())>-10 && (start.x()-end.x())<50.){
	event.isAccidental[it] = 1;
	continue; //Accidental K-
      }

      Double_t mass2 = qnan;
      TVector3 mom(qnan, qnan, qnan);
      target_pim_id_container.push_back(it);
      target_pim_dist2tgt_container.push_back(tp -> GetClosestDist());

      Int_t repid = 0;
      if(!GFTrackCont.TrackCheck(it, repid)){
	target_pim_mass2_container.push_back(mass2);
	target_pim_mom_container.push_back(mom);
	continue;
      }

      Int_t htofhitid_pi; Double_t tracklen_pi;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
								     event.HtofSeg, event.posHtof,
								     htofhitid_pi, tof, tracklen_pi,
								     pos, track2tgt_dist);
      mom = GFTrackCont.GetMom(it, 0, repid);
      if(htofextrapolation_pi) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_pi, event.tHtof[htofhitid_pi]);
      //if(htofextrapolation_pi && mass2 > 0.25) continue;
      target_pim_mass2_container.push_back(mass2);
      target_pim_mom_container.push_back(mom);
    }
    else if((event.pid[it]&4)!=4 && (event.pid[it]&1)==1 && event.charge[it]==1){ //pi+

      Double_t mass2 = qnan;
      TVector3 mom(qnan, qnan, qnan);
      target_pip_id_container.push_back(it);
      target_pip_dist2tgt_container.push_back(tp -> GetClosestDist());

      Int_t repid = 0;
      if(!GFTrackCont.TrackCheck(it, repid)){
	target_pip_mass2_container.push_back(mass2);
	target_pip_mom_container.push_back(mom);
	continue;
      }

      Int_t htofhitid_pi; Double_t tracklen_pi;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
								     event.HtofSeg, event.posHtof,
								     htofhitid_pi, tof, tracklen_pi,
								     pos, track2tgt_dist);
      mom = GFTrackCont.GetMom(it, 0, repid);
      if(htofextrapolation_pi) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_pi, event.tHtof[htofhitid_pi]);
      target_pip_mass2_container.push_back(mass2);
      target_pip_mom_container.push_back(mom);
    }
    else{ //p or pi

      Double_t mass2 = qnan;
      TVector3 mom(qnan, qnan, qnan);
      target_ppip_id_container.push_back(it);
      target_ppip_dist2tgt_container.push_back(tp -> GetClosestDist());

      Int_t repid = -1;
      if(!GFTrackCont.TrackCheck(it, repid)){
	target_ppip_mass2_container.push_back(mass2);
	target_ppip_mom_container.push_back(mom);
	continue;
      }

      Int_t htofhitid_ppi; Double_t tracklen_ppi;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation = GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
								  event.HtofSeg, event.posHtof,
								  htofhitid_ppi, tof, tracklen_ppi,
								  pos, track2tgt_dist);
      mom = GFTrackCont.GetMom(it, 0, repid);
      if(htofextrapolation) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_ppi, event.tHtof[htofhitid_ppi]);
      target_ppip_mass2_container.push_back(mass2);
      target_ppip_mom_container.push_back(mom);
    }
  }

  //Veto accidental tracks near the clustered vertex of accidental dummy tracks
  for(Int_t ivtx=0;ivtx<event.nvtxTpcClustered;ivtx++){
    TVector3 vtx = TVector3(event.clusteredVtx_x[ivtx], event.clusteredVtx_y[ivtx], event.clusteredVtx_z[ivtx]+tpc::ZTarget);
    if(TMath::Abs(vtx.y()) < 15.) continue;
    for(Int_t ith=0;ith<notintarget_id_container.size();ith++){
      Int_t it = notintarget_id_container[ith];
      if(event.isK18[it]==1) continue;
      if(event.isXi[it]==1) continue;
      if(event.isKurama[it]==1) continue;
      if(event.isBeam[it]==1) continue;
      if(event.isAccidental[it]==1) continue;

      Double_t par[5];
      par[0] = event.helix_cx[it];
      par[1] = event.helix_cy[it];
      par[2] = event.helix_z0[it];
      par[3] = event.helix_r[it];
      par[4] = event.helix_dz[it];

      Double_t theta_min = event.helix_t[it][0] - 100./par[3];
      Double_t theta_max = event.helix_t[it][0] + 100./par[3];
      Double_t dist = Kinematics::CalcHelixCloseDist(vtx, par, theta_min, theta_max);
      Int_t nhtrack = event.nhtrack[it];
      TVector3 start = TVector3(event.calpos_x[it][0], event.calpos_y[it][0], event.calpos_z[it][0]);
      TVector3 end = TVector3(event.calpos_x[it][nhtrack-1], event.calpos_y[it][nhtrack-1], event.calpos_z[it][nhtrack-1]);

      TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( it );
      if(dist<20.) event.isAccidental[it] = 1;
    }
  }

  std::vector<Int_t> decays_id_container;
  for(Int_t ivtx=0;ivtx<event.nvtxTpc;ivtx++){
    if(event.vtx_dist[ivtx]>10) continue;
    TVector3 vtx = TVector3(event.vtx_x[ivtx], event.vtx_y[ivtx], event.vtx_z[ivtx]);
    Int_t trackid1 = event.vtxid[ivtx][0];
    Int_t trackid2 = event.vtxid[ivtx][1];
    if(event.isK18[trackid1]==1) continue;
    if(event.isKurama[trackid1]==1) continue;
    if(event.isBeam[trackid1]==1) continue;
    if(event.isXi[trackid1]==1) continue;
    if(event.isAccidental[trackid1]==1) continue;
    if(event.isK18[trackid2]==1) continue;
    if(event.isKurama[trackid2]==1) continue;
    if(event.isBeam[trackid2]==1) continue;
    if(event.isXi[trackid2]==1) continue;
    if(event.isAccidental[trackid2]==1) continue;

    if(std::find(intarget_id_container.begin(), intarget_id_container.end(), trackid1) != intarget_id_container.end() &&
       std::find(intarget_id_container.begin(), intarget_id_container.end(), trackid2) != intarget_id_container.end()) continue;
    if(std::find(notintarget_id_container.begin(), notintarget_id_container.end(), trackid1) == notintarget_id_container.end() &&
       std::find(notintarget_id_container.begin(), notintarget_id_container.end(), trackid2) == notintarget_id_container.end()) continue;

    Int_t nhtrack1 = event.nhtrack[trackid1];
    Int_t nhtrack2 = event.nhtrack[trackid2];
    TVector3 start1 = TVector3(event.calpos_x[trackid1][0], event.calpos_y[trackid1][0], event.calpos_z[trackid1][0]);
    TVector3 start2 = TVector3(event.calpos_x[trackid2][0], event.calpos_y[trackid2][0], event.calpos_z[trackid2][0]);
    TVector3 end1 = TVector3(event.calpos_x[trackid1][nhtrack1-1], event.calpos_y[trackid1][nhtrack1-1], event.calpos_z[trackid1][nhtrack1-1]);
    TVector3 end2 = TVector3(event.calpos_x[trackid2][nhtrack2-1], event.calpos_y[trackid2][nhtrack2-1], event.calpos_z[trackid2][nhtrack2-1]);

    Double_t vertex_dist1; Double_t vertex_dist2;
    if(!Kinematics::HelixDirection(vtx, start1, end1, vertex_dist1) ||
       !Kinematics::HelixDirection(vtx, start2, end2, vertex_dist2)) continue;

    Bool_t parallel = false;
    if(event.vtx_angle[ivtx] < 0.1*TMath::Pi() ||
       event.vtx_angle[ivtx] > (1. - 0.1)*TMath::Pi()) parallel = true;
    if(vertex_dist1<20 && vertex_dist2<20 && parallel) continue; //parallel tracks(merging process has been failed)

    decays_id_container.push_back(trackid1);
    decays_id_container.push_back(trackid2);
  }

  std::sort(decays_id_container.begin(), decays_id_container.end());
  decays_id_container.erase(std::unique(decays_id_container.begin(), decays_id_container.end()), decays_id_container.end());
  for(Int_t id=0;id<decays_id_container.size();id++){ //1st searching
    Int_t it = decays_id_container[id];
    if(std::find(intarget_id_container.begin(), intarget_id_container.end(), it) !=
       intarget_id_container.end()) continue;

    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( it );
    Int_t nhit = event.nhtrack[it];
    TVector3 start(event.hitpos_x[it][0], event.hitpos_y[it][0], event.hitpos_z[it][0]);
    TVector3 end(event.hitpos_x[it][nhit-1], event.hitpos_y[it][nhit-1], event.hitpos_z[it][nhit-1]);
    if((event.pid[it]&2)==2 && event.charge[it]==-1){ //k-
      Bool_t km_pid = Kinematics::HypTPCdEdxPID_IsKaonTemp(event.dEdx[it], event.mom0[it]);
      if(km_pid){
	Int_t repid = 1;
	if(GFTrackCont.TrackCheck(it, repid)){
	  Int_t htofhitid_k; Double_t tracklen_k;
	  Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	  Bool_t htofextrapolation_k = GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
									event.HtofSeg, event.posHtof,
									htofhitid_k, tof, tracklen_k,
									pos, track2tgt_dist);
	  Double_t mass2 = qnan;
	  TVector3 mom = GFTrackCont.GetMom(it, 0, repid);
	  if(htofextrapolation_k) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_k, event.tHtof[htofhitid_k]);
	  target_k_id_container.push_back(it);
	  target_k_mass2_container.push_back(mass2);
	  target_k_mom_container.push_back(mom);
	}
      }
    }

    if((event.pid[it]&4)==4 && (event.pid[it]&1)!=1 && event.charge[it]==1){ //proton

      Double_t mass2 = qnan;
      TVector3 mom(qnan, qnan, qnan);
      target_p_id_container.push_back(it);
      target_p_dist2tgt_container.push_back(tp -> GetClosestDist());

      Int_t repid = 0;
      Int_t flag = 1;
      for(Int_t i=0;i<2;i++){
	Int_t temp = flag&event.pid[it];
	if(temp==flag) repid += 1;
	flag*=2;
      }
      if(!GFTrackCont.TrackCheck(it, repid)){
	target_p_mass2_container.push_back(mass2);
	target_p_mom_container.push_back(mom);
	continue;
      }

      Int_t htofhitid_p; Double_t tracklen_p;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation_p = GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
								    event.HtofSeg, event.posHtof,
								    htofhitid_p, tof, tracklen_p,
								    pos, track2tgt_dist);
      mom = GFTrackCont.GetMom(it, 0, repid);
      if(htofextrapolation_p) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_p, event.tHtof[htofhitid_p]);
      target_p_mass2_container.push_back(mass2);
      target_p_mom_container.push_back(mom);
    }
    else if((event.pid[it]&1)==1 && event.charge[it]==-1){ //pi-

      Double_t slope = event.helix_dz[it];
      Double_t helixmom = event.mom0[it];
      if(TMath::Abs(slope)<0.1 && TMath::Abs(helixmom)>0.5 &&
	 (start.x()-end.x())>-10 && (start.x()-end.x())<50. &&
	 TMath::Abs(start.x())<50. && TMath::Abs(end.x())<50.){
	event.isAccidental[it] = 1;
	continue; //Accidental K-
      }

      Double_t mass2 = qnan;
      TVector3 mom(qnan, qnan, qnan);
      target_pim_id_container.push_back(it);
      target_pim_dist2tgt_container.push_back(tp -> GetClosestDist());

      Int_t repid = 0;
      if(!GFTrackCont.TrackCheck(it, repid)){
	target_pim_mass2_container.push_back(mass2);
	target_pim_mom_container.push_back(mom);
	continue;
      }

      Int_t htofhitid_pi; Double_t tracklen_pi;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
								     event.HtofSeg, event.posHtof,
								     htofhitid_pi, tof, tracklen_pi,
								     pos, track2tgt_dist);
      mom = GFTrackCont.GetMom(it, 0, repid);
      if(htofextrapolation_pi) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_pi, event.tHtof[htofhitid_pi]);
      //if(htofextrapolation_pi && mass2 > 0.25) continue;
      target_pim_mass2_container.push_back(mass2);
      target_pim_mom_container.push_back(mom);
    }
    else if((event.pid[it]&4)!=4 && (event.pid[it]&1)==1 && event.charge[it]==1){ //pi+

      Double_t mass2 = qnan;
      TVector3 mom(qnan, qnan, qnan);
      target_pip_id_container.push_back(it);
      target_pip_dist2tgt_container.push_back(tp -> GetClosestDist());

      Int_t repid = 0;
      if(!GFTrackCont.TrackCheck(it, repid)){
	target_pip_mass2_container.push_back(mass2);
	target_pip_mom_container.push_back(mom);
	continue;
      }

      Int_t htofhitid_pi; Double_t tracklen_pi;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
								     event.HtofSeg, event.posHtof,
								     htofhitid_pi, tof, tracklen_pi,
								     pos, track2tgt_dist);
      mom = GFTrackCont.GetMom(it, 0, repid);
      if(htofextrapolation_pi) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_pi, event.tHtof[htofhitid_pi]);
      target_pip_mass2_container.push_back(mass2);
      target_pip_mom_container.push_back(mom);
    }
    else{ //p or pi

      Double_t mass2 = qnan;
      TVector3 mom(qnan, qnan, qnan);
      target_ppip_id_container.push_back(it);
      target_ppip_dist2tgt_container.push_back(tp -> GetClosestDist());

      Int_t repid = -1;
      if(!GFTrackCont.TrackCheck(it, repid)){
	target_ppip_mass2_container.push_back(mass2);
	target_ppip_mom_container.push_back(mom);
	continue;
      }

      Int_t htofhitid_ppi; Double_t tracklen_ppi;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation = GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
								  event.HtofSeg, event.posHtof,
								  htofhitid_ppi, tof, tracklen_ppi,
								  pos, track2tgt_dist);
      mom = GFTrackCont.GetMom(it, 0, repid);
      if(htofextrapolation) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_ppi, event.tHtof[htofhitid_ppi]);
      target_ppip_mass2_container.push_back(mass2);
      target_ppip_mom_container.push_back(mom);
    }
  }


  //for L, LL
  std::vector<Int_t> L_p_id_container, L_pi_id_container;
  std::vector<Int_t> L_p_repid_container, L_pi_repid_container;
  std::vector<Int_t> L_p_G4tid_container, L_pi_G4tid_container;
  std::vector<Double_t> L_p_purity_container, L_pi_purity_container;
  std::vector<Double_t> L_p_efficiency_container, L_pi_efficiency_container;
  std::vector<TVector3> L_p_mom_container, L_pi_mom_container;
  std::vector<Double_t> L_mass_container;
  std::vector<Double_t> L_ppidist_container;
  std::vector<Double_t> L_targetdist_container;
  std::vector<TVector3> L_mom_container, L_vtx_container, L_targetvtx_container;
  std::vector<Int_t> L_intarget_container;

  std::vector<Int_t> reconfailed_L_p_id_container, reconfailed_L_pi_id_container;
  std::vector<Double_t> reconfailed_L_mass_container, reconfailed_L_ppidist_container;
  std::vector<TVector3> reconfailed_L_mom_container, reconfailed_L_p_mom_container,
    reconfailed_L_pi_mom_container, reconfailed_L_vtx_container;

  //Xi/L candidates searching
#if DebugDisp
  std::cout<<"1. Xi/L candidate searching"<<std::endl;
#endif
  Int_t l_candidates = 0;
  Int_t xi_candidates = 0;

  std::vector<Int_t> xi_l_container;
  std::vector<Int_t> xi_l_intarget_container;
  std::vector<Int_t> xi_p_container, xi_pi_container, xi_pi2_container;
  std::vector<Int_t> p_repid_container, pi_repid_container, pi2_repid_container;
  std::vector<Int_t> xi_p_G4tid_container, xi_pi_G4tid_container, xi_pi2_G4tid_container;
  std::vector<TVector3> xi_mom_container, xi_decayvertex_container;
  std::vector<Double_t> xi_p_purity_container, xi_pi_purity_container, xi_pi2_purity_container;
  std::vector<Double_t> xi_p_efficiency_container, xi_pi_efficiency_container, xi_pi2_efficiency_container;
  std::vector<TVector3> l_mom_container, l_vert_container;
  std::vector<Double_t> xi_mass_container, lambda_mass_container;
  std::vector<TVector3> xi_p_mom_container, xi_pi_mom_container, xi_pi2_mom_container;
  std::vector<Double_t> ppi_closedist; std::vector<Double_t> lpi_closedist;
  std::vector<Double_t> xi_targetdist_container;
  std::vector<TVector3> xi_targetvtx_container, xi_targetmom_container;
  {
    for(Int_t it1=0;it1<ntTpc;it1++){
      if(event.isK18[it1]==1) continue;
      if(event.isKurama[it1]==1) continue;
      if(event.isBeam[it1]==1) continue;
      if(event.isXi[it1]==1) continue;
      if(event.isAccidental[it1]==1) continue;
      if((event.pid[it1]&4)!=4) continue;
      if(event.charge[it1]!=1) continue;

      Int_t repid_p = 0;
      Int_t flag = 1;
      for(Int_t i=0;i<2;i++){
	Int_t temp = flag&event.pid[it1];
	if(temp==flag) repid_p += 1;
	flag*=2;
      }
      if(!GFTrackCont.TrackCheck(it1, repid_p)) continue;

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
      for(Int_t it2=0;it2<ntTpc;it2++){
	if(it1==it2) continue;
	if(event.isK18[it2]==1) continue;
	if(event.isKurama[it2]==1) continue;
	if(event.isBeam[it2]==1) continue;
	if(event.isXi[it2]==1) continue;
	if(event.isAccidental[it2]==1) continue;
	if((event.pid[it2]&1)!=1) continue; //select pi-like
	//if((event.pid[it2]&4)==4) continue; //veto p-like
	if(event.charge[it2]!=-1) continue;

	Int_t repid_pi = 0;
	if(!GFTrackCont.TrackCheck(it2, repid_pi)) continue;

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

	Double_t ppi_dist = 10000.;
	TVector3 p_mom; TVector3 pi_mom; TVector3 lambda_mom;
	TVector3 lambda_vert = Kinematics::LambdaVertex(dMagneticField, p_par, pi_par, p_theta_min, p_theta_max, pi_theta_min, pi_theta_max, p_mom, pi_mom, lambda_mom, ppi_dist);
	if(TMath::IsNaN(ppi_dist)) continue;
	lambda_mom = pi_mom + p_mom;
	TLorentzVector Lp(p_mom, TMath::Hypot(p_mom.Mag(), ProtonMass));
	TLorentzVector Lpi(pi_mom, TMath::Hypot(pi_mom.Mag(), PionMass));
	TLorentzVector Llambda = Lp + Lpi;

	if(TMath::Abs(lambda_vert.x()) > 250. ||
	   TMath::Abs(lambda_vert.z()) > 250. ||
	   TMath::Abs(lambda_vert.y()) > 250.) continue; //Vertex cut

	Double_t pi_vertex_dist; Double_t p_vertex_dist;
	if(!Kinematics::HelixDirection(lambda_vert, p_start, p_end, p_vertex_dist) ||
	   !Kinematics::HelixDirection(lambda_vert, pi_start, pi_end, pi_vertex_dist)) continue;

	if(pi_vertex_dist > pi_vtx_distcut) continue;
 	if(p_vertex_dist > p_vtx_distcut) continue;
	if(ppi_dist > ppi_distcut || TMath::Abs(Llambda.M() - LambdaMass) > lambda_masscut) continue;
/*
	if(ppi_dist > ppi_distcut ||
	   TMath::Abs(Llambda.M() - LambdaMass) > lambda_masscut){
	  if(ppi_dist<ppi_distcut*2 && TMath::Abs(Llambda.M() - LambdaMass) < lambda_masscut*2.){
	    reconfailed_L_p_id_container.push_back(it1);
	    reconfailed_L_pi_id_container.push_back(it2);
	    reconfailed_L_mass_container.push_back(Llambda.M());
	    reconfailed_L_vtx_container.push_back(lambda_vert);
	    reconfailed_L_mom_container.push_back(lambda_mom);
	    reconfailed_L_p_mom_container.push_back(p_mom);
	    reconfailed_L_pi_mom_container.push_back(pi_mom);
	    reconfailed_L_ppidist_container.push_back(ppi_dist);
	  }
	  continue;
	}
*/
	Double_t ltarget_dist;
	TVector3 ltarget_vtx = Kinematics::CalcCloseDistLambda(tgtpos,
							       lambda_vert,
							       lambda_mom,
							       ltarget_dist);
	Double_t l_targetcenter_dist;
	TVector3 l_pos_tgtcenter = Kinematics::LambdaTargetCenter(lambda_vert, lambda_mom, l_targetcenter_dist);

	if(ltarget_dist > GFltarget_distcut or TMath::Abs(l_pos_tgtcenter.y()) > GFltarget_ycut){
	  L_intarget_container.push_back(0);
	}
	else{
	  L_intarget_container.push_back(1);
	}
	L_p_id_container.push_back(it1);
	L_pi_id_container.push_back(it2);
	L_p_purity_container.push_back(event.purity[it1]);
	L_pi_purity_container.push_back(event.purity[it2]);
	L_p_efficiency_container.push_back(event.efficiency[it1]);
	L_pi_efficiency_container.push_back(event.efficiency[it2]);
	L_p_G4tid_container.push_back(G4TrackID[it1]);
	L_pi_G4tid_container.push_back(G4TrackID[it2]);
	L_mass_container.push_back(Llambda.M());
	L_mom_container.push_back(lambda_mom);
	L_p_mom_container.push_back(p_mom);
	L_pi_mom_container.push_back(pi_mom);
	L_ppidist_container.push_back(ppi_dist);
	L_vtx_container.push_back(lambda_vert);
	L_p_repid_container.push_back(repid_p);
	L_pi_repid_container.push_back(repid_pi);
	L_targetdist_container.push_back(ltarget_dist);
	L_targetvtx_container.push_back(ltarget_vtx);
	l_candidates++;

	for(int it3=0;it3<ntTpc;it3++){
	  if(it3==it2 || it3==it1) continue;
	  if(event.isK18[it3]==1) continue;
	  if(event.isKurama[it3]==1) continue;
	  if(event.isBeam[it3]==1) continue;
	  if(event.isXi[it3]==1) continue;
	  if(event.isAccidental[it3]==1) continue;
	  if((event.pid[it3]&1)!=1) continue; //select pi like
	  //if((event.pid[it3]&4)==4) continue; //veto p-like
	  if(event.charge[it3]!=-1) continue;

	  Int_t repid_pi2 = 0;
	  if(!GFTrackCont.TrackCheck(it3, repid_pi2)) continue;

	  Double_t pi2_par[5];
	  pi2_par[0] = event.helix_cx[it3];
	  pi2_par[1] = event.helix_cy[it3];
	  pi2_par[2] = event.helix_z0[it3];
	  pi2_par[3] = event.helix_r[it3];
	  pi2_par[4] = event.helix_dz[it3];

	  Int_t pi2_nh = event.helix_t[it3].size();
	  Double_t pi2_theta_min = TMath::Max(event.helix_t[it3][0] - vtx_scan_rangeInsidePi/pi2_par[3], event.helix_t[it3][pi2_nh-1]);
	  Double_t pi2_theta_max = event.helix_t[it3][0] + vtx_scan_range/pi2_par[3];

	  TVector3 pi2_start = TVector3(event.calpos_x[it3][0], event.calpos_y[it3][0], event.calpos_z[it3][0]);
	  TVector3 pi2_end = TVector3(event.calpos_x[it3][pi2_nh-1], event.calpos_y[it3][pi2_nh-1], event.calpos_z[it3][pi2_nh-1]);

	  TVector3 pi2_mom; Double_t lpi_dist;

	  TVector3 xi_vert = Kinematics::XiVertex(dMagneticField,
						  pi2_par,
						  pi2_theta_min,
						  pi2_theta_max,
						  lambda_vert,
						  lambda_mom,
						  pi2_mom,
						  lpi_dist);
	  if(TMath::IsNaN(lpi_dist)) continue;

	  TLorentzVector Lpi2(pi2_mom, TMath::Sqrt(pi2_mom.Mag()*pi2_mom.Mag() + PionMass*PionMass));
	  TLorentzVector Llambda_fixedmass(lambda_mom, TMath::Sqrt(lambda_mom.Mag()*lambda_mom.Mag() + LambdaMass*LambdaMass));
	  TLorentzVector Lxi = Llambda_fixedmass + Lpi2;
	  TVector3 xi_mom = TVector3(Lxi.Px(), Lxi.Py(), Lxi.Pz());
	  Double_t pi2_vertex_dist;
	  if(lpi_dist > lpi_distcut) continue;
	  if(!Kinematics::HelixDirection(xi_vert, pi2_start, pi2_end, pi2_vertex_dist)) continue;
	  if(pi2_vertex_dist > pi2_vtx_distcut) continue;

	  if(TMath::Abs(Lxi.M() - XiMinusMass) > xi_masscut) continue; //Check reconstructed mass cut
	  Double_t xitarget_dist; TVector3 xi_mom_target;
	  TVector3 xi_vert_target = Kinematics::CalcCloseDistXi(tgtpos,
								dMagneticField,
								xi_vert,
								xi_mom,
								xi_mom_target,
								xitarget_dist);
	  if(xitarget_dist > xitarget_distcut) continue; //Closest distance between the xi and the target cut
	  xi_targetdist_container.push_back(xitarget_dist);
	  xi_targetvtx_container.push_back(xi_vert_target);
	  xi_targetmom_container.push_back(xi_mom_target);

	  ppi_closedist.push_back(ppi_dist);
	  lpi_closedist.push_back(lpi_dist);

	  xi_l_container.push_back(l_candidates - 1);
    xi_l_intarget_container.push_back(L_intarget_container[l_candidates - 1]);
	  xi_p_container.push_back(it1);
	  xi_pi_container.push_back(it2);
	  xi_pi2_container.push_back(it3);

    xi_p_purity_container.push_back(event.purity[it1]);
    xi_pi_purity_container.push_back(event.purity[it2]);
    xi_pi2_purity_container.push_back(event.purity[it3]);
    xi_p_efficiency_container.push_back(event.efficiency[it1]);
    xi_pi_efficiency_container.push_back(event.efficiency[it2]);
    xi_pi2_efficiency_container.push_back(event.efficiency[it3]);
    xi_p_G4tid_container.push_back(G4TrackID[it1]);
    xi_pi_G4tid_container.push_back(G4TrackID[it2]);
    xi_pi2_G4tid_container.push_back(G4TrackID[it3]);

	  p_repid_container.push_back(repid_p);
	  pi_repid_container.push_back(repid_pi);
	  pi2_repid_container.push_back(repid_pi2);

	  xi_mass_container.push_back(Lxi.M());
	  lambda_mass_container.push_back(Llambda.M());

	  xi_mom_container.push_back(xi_mom);
	  l_mom_container.push_back(lambda_mom);
	  xi_p_mom_container.push_back(p_mom);
	  xi_pi_mom_container.push_back(pi_mom);
	  xi_pi2_mom_container.push_back(pi2_mom);

	  xi_decayvertex_container.push_back(xi_vert);
	  l_vert_container.push_back(lambda_vert);

	  xi_candidates++;
	} //it3
      } //it2
    } //it1
  }

#if DebugDisp
  std::cout<<"1-1. pi- pi+ pair candidate searching"<<std::endl;
#endif

  //pi+&pi- pair
  std::vector<Int_t> pipair_pip_id_container, pipair_pim_id_container;
  std::vector<Double_t> pipair_reconL_mass_container, pipair_pipidist_container;
  std::vector<TVector3> pipair_mom_container, pipair_pip_mom_container, pipair_pim_mom_container;
  {
    for(Int_t it1=0;it1<ntTpc;it1++){
      if(event.isK18[it1]==1) continue;
      if(event.isKurama[it1]==1) continue;
      if(event.isBeam[it1]==1) continue;
      if(event.isXi[it1]==1) continue;
      if(event.isAccidental[it1]==1) continue;
      if((event.pid[it1]&4)==4) continue; //not proton, pi+ like
      if(event.charge[it1]!=1) continue;

      Int_t repid_pip = 0;
      if(!GFTrackCont.TrackCheck(it1, repid_pip)) continue;

      Double_t pip_par[5];
      pip_par[0] = event.helix_cx[it1];
      pip_par[1] = event.helix_cy[it1];
      pip_par[2] = event.helix_z0[it1];
      pip_par[3] = event.helix_r[it1];
      pip_par[4] = event.helix_dz[it1];
      Int_t pip_nh = event.helix_t[it1].size();
      Double_t pip_theta_min = event.helix_t[it1][0] - vtx_scan_range/pip_par[3];
      Double_t pip_theta_max = TMath::Min(event.helix_t[it1][0] + vtx_scan_rangeInsideL/pip_par[3], event.helix_t[it1][pip_nh-1]);
      TVector3 pip_start = TVector3(event.calpos_x[it1][0], event.calpos_y[it1][0], event.calpos_z[it1][0]);
      TVector3 pip_end = TVector3(event.calpos_x[it1][pip_nh-1], event.calpos_y[it1][pip_nh-1], event.calpos_z[it1][pip_nh-1]);
      for(Int_t it2=0;it2<ntTpc;it2++){
	if(it1==it2) continue;
	if(event.isK18[it2]==1) continue;
	if(event.isKurama[it2]==1) continue;
	if(event.isBeam[it2]==1) continue;
	if(event.isXi[it2]==1) continue;
	if(event.isAccidental[it2]==1) continue;
	if((event.pid[it2]&1)!=1) continue; //select pi-like
	//if((event.pid[it2]&4)==4) continue; //veto p-like
	if(event.charge[it2]!=-1) continue;

	Int_t repid_pim = 0;
	if(!GFTrackCont.TrackCheck(it2, repid_pim)) continue;

	Double_t pim_par[5];
	pim_par[0] = event.helix_cx[it2];
	pim_par[1] = event.helix_cy[it2];
	pim_par[2] = event.helix_z0[it2];
	pim_par[3] = event.helix_r[it2];
	pim_par[4] = event.helix_dz[it2];

	Int_t pim_nh = event.helix_t[it2].size();
	Double_t pim_theta_min = TMath::Max(event.helix_t[it2][0] - vtx_scan_rangeInsideL/pim_par[3], event.helix_t[it2][pim_nh-1]);
	Double_t pim_theta_max = event.helix_t[it2][0] + vtx_scan_range/pim_par[3];
	TVector3 pim_start = TVector3(event.calpos_x[it2][0], event.calpos_y[it2][0], event.calpos_z[it2][0]);
	TVector3 pim_end = TVector3(event.calpos_x[it2][pim_nh-1], event.calpos_y[it2][pim_nh-1], event.calpos_z[it2][pim_nh-1]);

	Double_t pipi_dist = 10000.;
	TVector3 pip_mom; TVector3 pim_mom; TVector3 lambda_mom;
	TVector3 lambda_vert = Kinematics::LambdaVertex(dMagneticField, pip_par, pim_par, pip_theta_min, pip_theta_max, pim_theta_min, pim_theta_max, pip_mom, pim_mom, lambda_mom, pipi_dist);
	if(TMath::IsNaN(pipi_dist)) continue;
	lambda_mom = pim_mom + pip_mom;
	TLorentzVector Lpip(pip_mom, TMath::Hypot(pip_mom.Mag(), ProtonMass));
	TLorentzVector Lpim(pim_mom, TMath::Hypot(pim_mom.Mag(), PionMass));
	TLorentzVector Llambda = Lpip + Lpim;

	if(TMath::Abs(lambda_vert.x()) > 250. ||
	   TMath::Abs(lambda_vert.z()) > 250. ||
	   TMath::Abs(lambda_vert.y()) > 250.) continue; //Vertex cut

	Double_t pip_vertex_dist; Double_t pim_vertex_dist;
	if(!Kinematics::HelixDirection(lambda_vert, pip_start, pip_end, pip_vertex_dist) ||
	   !Kinematics::HelixDirection(lambda_vert, pim_start, pim_end, pim_vertex_dist)) continue;

	if(pip_vertex_dist > pi_vtx_distcut) continue;
 	if(pim_vertex_dist > pi_vtx_distcut) continue;

	if(pipi_dist < pipi_distcut){
	  pipair_pip_id_container.push_back(it1);
	  pipair_pim_id_container.push_back(it2);
	  pipair_reconL_mass_container.push_back(Llambda.M());
	  pipair_pipidist_container.push_back(pipi_dist);
	  pipair_mom_container.push_back(lambda_mom);
	  pipair_pip_mom_container.push_back(pip_mom);
	  pipair_pim_mom_container.push_back(pim_mom);
	}
      } //it2
    } //it1
  }

#if DebugDisp
  std::cout<<"1-2. p-bar like pi-  pair candidate searching"<<std::endl;
#endif
    {
    for(Int_t it1=0;it1<ntTpc;it1++){
      if(event.isK18[it1]==1) continue;
      if(event.isKurama[it1]==1) continue;
      if(event.isBeam[it1]==1) continue;
      if(event.isXi[it1]==1) continue;
      if(event.isAccidental[it1]==1) continue;
      if(event.chargeIndistinguishable[it1]!=1) continue;
      if((event.pid_inverted[it1]&4)!=4) continue;
      continue;
      if(event.charge[it1]!=-1) continue;

      Double_t p_par[5];
      p_par[0] = event.helix_cx_inverted[it1];
      p_par[1] = event.helix_cy_inverted[it1];
      p_par[2] = event.helix_z0_inverted[it1];
      p_par[3] = event.helix_r_inverted[it1];
      p_par[4] = event.helix_dz_inverted[it1];

      Int_t p_nh = event.helix_t[it1].size();
      Double_t p_theta_min = TMath::Min((event.calpos_y[it1][0]-p_par[2])/(p_par[3]*p_par[4]), (event.calpos_y[it1][p_nh-1]-p_par[2])/(p_par[3]*p_par[4])) - vtx_scan_rangeInsideL/p_par[3];
      Double_t p_theta_max = TMath::Max((event.calpos_y[it1][0]-p_par[2])/(p_par[3]*p_par[4]), (event.calpos_y[it1][p_nh-1]-p_par[2])/(p_par[3]*p_par[4])) + vtx_scan_rangeInsideL/p_par[3];
      TVector3 p_start = TVector3(event.calpos_x[it1][0], event.calpos_y[it1][0], event.calpos_z[it1][0]);
      TVector3 p_end = TVector3(event.calpos_x[it1][p_nh-1], event.calpos_y[it1][p_nh-1], event.calpos_z[it1][p_nh-1]);
      for(Int_t it2=0;it2<ntTpc;it2++){
	if(it1==it2) continue;
	if(event.isK18[it2]==1) continue;
	if(event.isKurama[it2]==1) continue;
	if(event.isBeam[it2]==1) continue;
	if(event.isXi[it2]==1) continue;
	if(event.isAccidental[it2]==1) continue;
	if((event.pid[it2]&1)!=1) continue; //select pi-like
	//if((event.pid[it2]&4)==4) continue; //veto p-like
	if(event.charge[it2]!=-1) continue;

	Int_t repid_pi = 0;
	if(!GFTrackCont.TrackCheck(it2, repid_pi)) continue;

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

	Double_t ppi_dist = 10000.;
	TVector3 p_mom; TVector3 pi_mom; TVector3 lambda_mom;
	TVector3 lambda_vert = Kinematics::LambdaVertex(dMagneticField, p_par, pi_par, p_theta_min, p_theta_max, pi_theta_min, pi_theta_max, p_mom, pi_mom, lambda_mom, ppi_dist);
	if(TMath::IsNaN(ppi_dist)) continue;
	lambda_mom = pi_mom + p_mom;
	TLorentzVector Lp(p_mom, TMath::Hypot(p_mom.Mag(), ProtonMass));
	TLorentzVector Lpi(pi_mom, TMath::Hypot(pi_mom.Mag(), PionMass));
	TLorentzVector Llambda = Lp + Lpi;

	if(TMath::Abs(lambda_vert.x()) > 250. ||
	   TMath::Abs(lambda_vert.z()) > 250. ||
	   TMath::Abs(lambda_vert.y()) > 250.) continue; //Vertex cut

	Double_t pi_vertex_dist; Double_t p_vertex_dist;
	if(!Kinematics::HelixDirection(lambda_vert, p_start, p_end, p_vertex_dist) ||
	   !Kinematics::HelixDirection(lambda_vert, pi_start, pi_end, pi_vertex_dist)) continue;

	if(pi_vertex_dist > pi_vtx_distcut) continue;
 	if(p_vertex_dist > p_vtx_distcut) continue;

	Double_t ltarget_dist;
	TVector3 ltarget_vtx = Kinematics::CalcCloseDistLambda(tgtpos,
							       lambda_vert,
							       lambda_mom,
							       ltarget_dist);
	if(ltarget_dist > ltarget_distcut) continue;
	if(ppi_dist<ppi_distcut && TMath::Abs(Llambda.M() - LambdaMass) < lambda_masscut){
	  reconfailed_L_p_id_container.push_back(it1);
	  reconfailed_L_pi_id_container.push_back(it2);
	  reconfailed_L_mass_container.push_back(Llambda.M());
	  reconfailed_L_vtx_container.push_back(lambda_vert);
	  reconfailed_L_mom_container.push_back(lambda_mom);
	  reconfailed_L_p_mom_container.push_back(p_mom);
	  reconfailed_L_pi_mom_container.push_back(pi_mom);
	  reconfailed_L_ppidist_container.push_back(ppi_dist);
	}
      } //it2
    } //it1
  }

#if DebugDisp
  std::cout<<"1-3. p  p like pair candidate searching"<<std::endl;
#endif
    {
    for(Int_t it1=0;it1<ntTpc;it1++){
      if(event.isK18[it1]==1) continue;
      if(event.isKurama[it1]==1) continue;
      if(event.isBeam[it1]==1) continue;
      if(event.isXi[it1]==1) continue;
      if(event.isAccidental[it1]==1) continue;
      if((event.pid[it1]&4)!=4) continue;
      if(event.charge[it1]!=1) continue;

      Int_t repid_p = 0;
      Int_t flag = 1;
      for(Int_t i=0;i<2;i++){
	Int_t temp = flag&event.pid[it1];
	if(temp==flag) repid_p += 1;
	flag*=2;
      }
      if(!GFTrackCont.TrackCheck(it1, repid_p)) continue;

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
      for(Int_t it2=0;it2<ntTpc;it2++){
	if(it1==it2) continue;
	if(event.isK18[it2]==1) continue;
	if(event.isKurama[it2]==1) continue;
	if(event.isBeam[it2]==1) continue;
	if(event.isXi[it2]==1) continue;
	if(event.isAccidental[it2]==1) continue;
	if(event.chargeIndistinguishable[it1]!=1) continue;
	if((event.pid_inverted[it1]&1)!=1) continue;
  continue;
	if(event.charge[it1]!=1) continue;

	Double_t pi_par[5];
	pi_par[0] = event.helix_cx_inverted[it2];
	pi_par[1] = event.helix_cy_inverted[it2];
	pi_par[2] = event.helix_z0_inverted[it2];
	pi_par[3] = event.helix_r_inverted[it2];
	pi_par[4] = event.helix_dz_inverted[it2];

	Int_t pi_nh = event.helix_t[it2].size();
	Double_t pi_theta_min = TMath::Min((event.calpos_y[it2][0]-pi_par[2])/(pi_par[3]*pi_par[4]), (event.calpos_y[it2][pi_nh-1]-pi_par[2])/(pi_par[3]*pi_par[4])) - vtx_scan_rangeInsideL/TMath::Abs(pi_par[3]);
	Double_t pi_theta_max = TMath::Max((event.calpos_y[it2][0]-pi_par[2])/(pi_par[3]*pi_par[4]), (event.calpos_y[it2][pi_nh-1]-pi_par[2])/(pi_par[3]*pi_par[4])) + vtx_scan_rangeInsideL/TMath::Abs(pi_par[3]);
	TVector3 pi_start = TVector3(event.calpos_x[it2][0], event.calpos_y[it2][0], event.calpos_z[it2][0]);
	TVector3 pi_end = TVector3(event.calpos_x[it2][pi_nh-1], event.calpos_y[it2][pi_nh-1], event.calpos_z[it2][pi_nh-1]);

	Double_t ppi_dist = 10000.;
	TVector3 p_mom; TVector3 pi_mom; TVector3 lambda_mom;
	TVector3 lambda_vert = Kinematics::LambdaVertex(dMagneticField, p_par, pi_par, p_theta_min, p_theta_max, pi_theta_min, pi_theta_max, p_mom, pi_mom, lambda_mom, ppi_dist);
	if(TMath::IsNaN(ppi_dist)) continue;
	lambda_mom = pi_mom + p_mom;
	TLorentzVector Lp(p_mom, TMath::Hypot(p_mom.Mag(), ProtonMass));
	TLorentzVector Lpi(pi_mom, TMath::Hypot(pi_mom.Mag(), PionMass));
	TLorentzVector Llambda = Lp + Lpi;

	if(TMath::Abs(lambda_vert.x()) > 250. ||
	   TMath::Abs(lambda_vert.z()) > 250. ||
	   TMath::Abs(lambda_vert.y()) > 250.) continue; //Vertex cut

	Double_t pi_vertex_dist; Double_t p_vertex_dist;
	if(!Kinematics::HelixDirection(lambda_vert, p_start, p_end, p_vertex_dist) ||
	   !Kinematics::HelixDirection(lambda_vert, pi_start, pi_end, pi_vertex_dist)) continue;

	if(pi_vertex_dist > pi_vtx_distcut) continue;
 	if(p_vertex_dist > p_vtx_distcut) continue;
	Double_t ltarget_dist;
	TVector3 ltarget_vtx = Kinematics::CalcCloseDistLambda(tgtpos,
							       lambda_vert,
							       lambda_mom,
							       ltarget_dist);
	if(ltarget_dist > ltarget_distcut) continue;
	if(ppi_dist<ppi_distcut && TMath::Abs(Llambda.M() - LambdaMass) < lambda_masscut){
	  reconfailed_L_p_id_container.push_back(it1);
	  reconfailed_L_pi_id_container.push_back(it2);
	  reconfailed_L_mass_container.push_back(Llambda.M());
	  reconfailed_L_vtx_container.push_back(lambda_vert);
	  reconfailed_L_mom_container.push_back(lambda_mom);
	  reconfailed_L_p_mom_container.push_back(p_mom);
	  reconfailed_L_pi_mom_container.push_back(pi_mom);
	  reconfailed_L_ppidist_container.push_back(ppi_dist);
	}
      } //it2
    } //it1
  }

#if DebugDisp
  std::cout<<"Xi/L candidates searching ends"<<std::endl;
#endif

  std::vector<Double_t> GFL_mass_container(l_candidates, qnan);
  std::vector<Double_t> GFL_ppidist_container(l_candidates, qnan);
  std::vector<Double_t> GFL_targetdist_container(l_candidates, qnan);
  std::vector<TVector3> GFL_targetvtx_container(l_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFL_targetcenterdist_container(l_candidates, qnan);
  std::vector<TVector3> GFL_targetcentervtx_container(l_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFL_mom_container(l_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFL_vtx_container(l_candidates, TVector3(qnan, qnan, qnan));

  std::vector<Int_t> GFL_p_id_container(l_candidates, qnan);
  std::vector<Int_t> GFL_pi_id_container(l_candidates, qnan);
  std::vector<Int_t> GFL_p_repid_container(l_candidates, qnan);
  std::vector<Int_t> GFL_pi_repid_container(l_candidates, qnan);
  std::vector<TVector3> GFL_p_mom_container(l_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFL_pi_mom_container(l_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFL_p_extrapolation_container(l_candidates, qnan);
  std::vector<Double_t> GFL_pi_extrapolation_container(l_candidates, qnan);
  std::vector<Int_t> GFL_p_htofid_container(l_candidates, qnan);
  std::vector<Int_t> GFL_pi_htofid_container(l_candidates, qnan);
  std::vector<Double_t> GFL_p_tracklen_container(l_candidates, qnan);
  std::vector<Double_t> GFL_pi_tracklen_container(l_candidates, qnan);
  std::vector<Double_t> GFL_p_tof_container(l_candidates, qnan);
  std::vector<Double_t> GFL_pi_tof_container(l_candidates, qnan);
  std::vector<Double_t> GFL_p_mass2_container(l_candidates, qnan);
  std::vector<Double_t> GFL_pi_mass2_container(l_candidates, qnan);
  if(l_candidates>0){
#if DebugDisp
    std::cout<<"2. Single L, LL candidates searching starts"<<std::endl;
#endif
    //Reconstructed real Lambdas
    for(int idp=0;idp<l_candidates;idp++){
      if(L_targetdist_container[idp] > ltarget_distcut) continue;
      Int_t p_id = L_p_id_container[idp];
      Int_t pi_id = L_pi_id_container[idp];
      Int_t p_repid = L_p_repid_container[idp];
      Int_t pi_repid = L_pi_repid_container[idp];
      Double_t p_extrapolation; Double_t pi_extrapolation;
      TVector3 p_mom; TVector3 pi_mom;
      Double_t ppi_dist; TVector3 l_vertex;

      Bool_t vtxcut = (GFTrackCont.FindVertex(p_id, pi_id,
					      p_repid, pi_repid,
					      p_extrapolation, pi_extrapolation,
					      p_mom, pi_mom,
					      ppi_dist, l_vertex,
					      vtx_scan_range)
		       && ppi_dist < GFppi_distcut);
      if(!vtxcut) continue;

      TVector3 l_mom = p_mom + pi_mom;
      Double_t l_target_dist;
      TVector3 l_pos_tgt = Kinematics::CalcCloseDistLambda(tgtpos, l_vertex, l_mom, l_target_dist);
      TVector3 l_fight = l_vertex - l_pos_tgt;
      Double_t l_tof = Kinematics::CalcTimeOfFlight(l_mom.Mag(), l_fight.Mag(), pdg::LambdaMass());
      if(l_target_dist > GFltarget_distcut) continue;

      Double_t l_targetcenter_dist;
      TVector3 l_pos_tgtcenter = Kinematics::LambdaTargetCenter(l_vertex, l_mom, l_targetcenter_dist);
      if(TMath::Abs(l_pos_tgtcenter.y()) > GFltarget_ycut) continue;

      Int_t hitid_htof; Double_t tof_htof; Double_t tracklen_htof; TVector3 pos_htof; Double_t track2tgt_dist;
      Bool_t htofextrapolation = GFTrackCont.TPCHTOFTrackMatching(p_id, p_repid, l_vertex,
								  event.HtofSeg, event.posHtof,
								  hitid_htof, tof_htof,
								  tracklen_htof, pos_htof, track2tgt_dist);
      if(htofextrapolation){
	GFL_p_htofid_container[idp] = hitid_htof;
	GFL_p_tracklen_container[idp] = tracklen_htof;
	GFL_p_tof_container[idp] = event.tHtof[hitid_htof] - l_tof;
	GFL_p_mass2_container[idp] = Kinematics::MassSquare(p_mom.Mag(), tracklen_htof, event.tHtof[hitid_htof] - l_tof);
	//if(GFL_p_mass2_container[idp] < 0.25) continue;
      }

      htofextrapolation = GFTrackCont.TPCHTOFTrackMatching(pi_id, pi_repid, l_vertex,
							   event.HtofSeg, event.posHtof,
							   hitid_htof, tof_htof,
							   tracklen_htof, pos_htof, track2tgt_dist);
      if(htofextrapolation){
	GFL_pi_htofid_container[idp] = hitid_htof;
	GFL_pi_tracklen_container[idp] = tracklen_htof;
	GFL_pi_tof_container[idp] = event.tHtof[hitid_htof] - l_tof;
	GFL_pi_mass2_container[idp] = Kinematics::MassSquare(pi_mom.Mag(), tracklen_htof, event.tHtof[hitid_htof] - l_tof);
	//if(GFL_pi_mass2_container[idp] > 0.25) continue;
      }

      TLorentzVector GFLp(p_mom, TMath::Hypot(p_mom.Mag(), ProtonMass));
      TLorentzVector GFLpi(pi_mom, TMath::Hypot(pi_mom.Mag(), PionMass));
      TLorentzVector GFLlambda = GFLp + GFLpi;

      GFL_p_id_container[idp] = p_id;
      GFL_pi_id_container[idp] = pi_id;
      GFL_p_repid_container[idp] = p_repid;
      GFL_pi_repid_container[idp] = pi_repid;
      GFL_p_mom_container[idp] = p_mom;
      GFL_pi_mom_container[idp] = pi_mom;
      GFL_p_extrapolation_container[idp] = p_extrapolation;
      GFL_pi_extrapolation_container[idp] = pi_extrapolation;
      GFL_mass_container[idp] = GFLlambda.M();
      GFL_mom_container[idp] = l_mom;
      GFL_ppidist_container[idp] = ppi_dist;
      GFL_vtx_container[idp] = l_vertex;
      GFL_targetdist_container[idp] = l_target_dist;
      GFL_targetvtx_container[idp] = l_pos_tgt;
      GFL_targetcenterdist_container[idp] = l_targetcenter_dist;
      GFL_targetcentervtx_container[idp] = l_pos_tgtcenter;
    }

    //order : L1, L2
    std::vector<std::vector<Int_t>> GFLLid_container;
    std::vector<std::vector<TVector3>> GFLL_Ldecayvtx_container;
    std::vector<std::vector<TVector3>> GFLLmom_container;
    std::vector<std::vector<Double_t>> GFLLmass_container;
    std::vector<std::vector<Double_t>> GFLLppidist_container;
    std::vector<std::vector<Double_t>> GFLL_Ltargetdist_container;
    std::vector<std::vector<TVector3>> GFLL_Ltargetvtx_container;
    std::vector<std::vector<Double_t>> GFLL_Ltargetcenterdist_container;
    std::vector<std::vector<TVector3>> GFLL_Ltargetcentervtx_container;
    //order : p, pi(L1), p, pi(L2)
    std::vector<std::vector<Int_t>> GFLLdecays_trackid_container;
    std::vector<std::vector<Int_t>> GFLLdecays_repid_container;
    std::vector<std::vector<Int_t>> GFLLdecays_htofid_container;
    std::vector<std::vector<Double_t>> GFLLdecays_tof_container;
    std::vector<std::vector<TVector3>> GFLLdecays_mom_container;
    std::vector<std::vector<Double_t>> GFLLdecays_tracklen_container;
    std::vector<std::vector<Double_t>> GFLLdecays_mass2_container;
    Int_t LL_best = -1; Int_t LLcount = 0;
    Double_t prev_LLfit_chisqr = 9999.;
    for(Int_t idp1=0;idp1<l_candidates;idp1++){
      if(TMath::IsNaN(GFL_mass_container[idp1])) continue;
      //For LL searching
      for(Int_t idp2=idp1+1;idp2<l_candidates;idp2++){
	if(TMath::IsNaN(GFL_mass_container[idp2])) continue;
	if(GFL_p_id_container[idp1]==GFL_p_id_container[idp2] ||
	   GFL_p_id_container[idp1]==GFL_pi_id_container[idp2]) continue;
	if(GFL_pi_id_container[idp1]==GFL_p_id_container[idp2] ||
	   GFL_pi_id_container[idp1]==GFL_pi_id_container[idp2]) continue;
	std::cout<<"idp "<<idp1<<" "<<idp2<<std::endl;
	Int_t trackcount = 0;
	if(std::find(intarget_id_container.begin(), intarget_id_container.end(), GFL_p_id_container[idp1]) != intarget_id_container.end()) trackcount++;
	if(std::find(intarget_id_container.begin(), intarget_id_container.end(), GFL_p_id_container[idp2]) != intarget_id_container.end()) trackcount++;
	if(std::find(intarget_id_container.begin(), intarget_id_container.end(), GFL_pi_id_container[idp1]) != intarget_id_container.end()) trackcount++;
	if(std::find(intarget_id_container.begin(), intarget_id_container.end(), GFL_pi_id_container[idp2]) != intarget_id_container.end()) trackcount++;

	if(std::find(decays_id_container.begin(), decays_id_container.end(), GFL_p_id_container[idp1]) != decays_id_container.end()) trackcount++;
	if(std::find(decays_id_container.begin(), decays_id_container.end(), GFL_p_id_container[idp2]) != decays_id_container.end()) trackcount++;
	if(std::find(decays_id_container.begin(), decays_id_container.end(), GFL_pi_id_container[idp1]) != decays_id_container.end()) trackcount++;
	if(std::find(decays_id_container.begin(), decays_id_container.end(), GFL_pi_id_container[idp2]) != decays_id_container.end()) trackcount++;
	if(intarget_id_container.size() + decays_id_container.size()-trackcount!=0) continue;
	std::cout<<"check1 "<<"idp "<<idp1<<" "<<idp2<<std::endl;
	//if extra pi-, pi+ pair exist, it's not LL events
	if(std::find(pipair_pip_id_container.begin(), pipair_pip_id_container.end(), GFL_p_id_container[idp1]) != pipair_pip_id_container.end()) continue;
	if(std::find(pipair_pip_id_container.begin(), pipair_pip_id_container.end(), GFL_p_id_container[idp2]) != pipair_pip_id_container.end()) continue;
	if(std::find(pipair_pim_id_container.begin(), pipair_pim_id_container.end(), GFL_pi_id_container[idp1]) != pipair_pim_id_container.end()) continue;
	if(std::find(pipair_pim_id_container.begin(), pipair_pim_id_container.end(), GFL_pi_id_container[idp2]) != pipair_pim_id_container.end()) continue;
	std::cout<<"check2 "<<"idp "<<idp1<<" "<<idp2<<std::endl;
  //if extra p, pi- pair exist, it's not LL events
	if(std::find(reconfailed_L_p_id_container.begin(), reconfailed_L_p_id_container.end(), GFL_p_id_container[idp1]) != pipair_pip_id_container.end()) continue;
	if(std::find(reconfailed_L_p_id_container.begin(), reconfailed_L_p_id_container.end(), GFL_p_id_container[idp2]) != pipair_pip_id_container.end()) continue;
	if(std::find(reconfailed_L_pi_id_container.begin(), reconfailed_L_pi_id_container.end(), GFL_pi_id_container[idp1]) != pipair_pim_id_container.end()) continue;
	if(std::find(reconfailed_L_pi_id_container.begin(), reconfailed_L_pi_id_container.end(), GFL_pi_id_container[idp2]) != pipair_pim_id_container.end()) continue;
	std::cout<<"check3 "<<"idp "<<idp1<<" "<<idp2<<std::endl;
	std::vector<Int_t> GFlambda_containerid(2);
	std::vector<TVector3> GFlambda_vert(2);
	std::vector<TVector3> GFlambda_mom(2);
	std::vector<Double_t> GFlambda_mass(2);
	std::vector<Double_t> GFppi_dist(2);
	std::vector<Double_t> GFlambdatgt_dist(2);
	std::vector<TVector3> GFlambdatgt_vtx(2);
	std::vector<Double_t> GFlambdatgtcenter_dist(2);
	std::vector<TVector3> GFlambdatgtcenter_vtx(2);

	std::vector<Int_t> GFtrackid_decays(4);
	std::vector<Int_t> GFrepid_decays(4);
	std::vector<Int_t> GFhtofid_decays(4);
	std::vector<Double_t> GFtof_decays(4);
	std::vector<TVector3> GFmom_decays(4);
	std::vector<Double_t> GFtracklen_decays(4);
	std::vector<Double_t> GFmass2_decays(4);

	GFlambda_containerid[0] = idp1;
	GFlambda_containerid[1] = idp2;
	GFlambda_vert[0] = GFL_vtx_container[idp1];
	GFlambda_vert[1] = GFL_vtx_container[idp2];
	GFlambda_mom[0] = GFL_mom_container[idp1];
	GFlambda_mom[1] = GFL_mom_container[idp2];
	GFlambda_mass[0] = GFL_mass_container[idp1];
	GFlambda_mass[1] = GFL_mass_container[idp2];
	GFppi_dist[0] = GFL_ppidist_container[idp1];
	GFppi_dist[1] = GFL_ppidist_container[idp2];
	GFlambdatgt_dist[0] = GFL_targetdist_container[idp1];
	GFlambdatgt_dist[1] = GFL_targetdist_container[idp2];
	GFlambdatgt_vtx[0] = GFL_targetvtx_container[idp1];
	GFlambdatgt_vtx[1] = GFL_targetvtx_container[idp2];
	GFlambdatgtcenter_dist[0] = GFL_targetcenterdist_container[idp1];
	GFlambdatgtcenter_dist[1] = GFL_targetcenterdist_container[idp2];
	GFlambdatgtcenter_vtx[0] = GFL_targetcentervtx_container[idp1];
	GFlambdatgtcenter_vtx[1] = GFL_targetcentervtx_container[idp2];

	GFtrackid_decays[0] = GFL_p_id_container[idp1];
	GFtrackid_decays[1] = GFL_pi_id_container[idp1];
	GFtrackid_decays[2] = GFL_p_id_container[idp2];
	GFtrackid_decays[3] = GFL_pi_id_container[idp2];
	GFrepid_decays[0] = GFL_p_repid_container[idp1];
	GFrepid_decays[1] = GFL_pi_repid_container[idp1];
	GFrepid_decays[2] = GFL_p_repid_container[idp2];
	GFrepid_decays[3] = GFL_pi_repid_container[idp2];
	GFhtofid_decays[0] = GFL_p_htofid_container[idp1];
	GFhtofid_decays[1] = GFL_pi_htofid_container[idp1];
	GFhtofid_decays[2] = GFL_p_htofid_container[idp2];
	GFhtofid_decays[3] = GFL_pi_htofid_container[idp2];
	GFtof_decays[0] = GFL_p_tof_container[idp1];
	GFtof_decays[1] = GFL_pi_tof_container[idp1];
	GFtof_decays[2] = GFL_p_tof_container[idp2];
	GFtof_decays[3] = GFL_pi_tof_container[idp2];
	GFmom_decays[0] = GFL_p_mom_container[idp1];
	GFmom_decays[1] = GFL_pi_mom_container[idp1];
	GFmom_decays[2] = GFL_p_mom_container[idp2];
	GFmom_decays[3] = GFL_pi_mom_container[idp2];
	GFtracklen_decays[0] = GFL_p_tracklen_container[idp1];
	GFtracklen_decays[1] = GFL_pi_tracklen_container[idp1];
	GFtracklen_decays[2] = GFL_p_tracklen_container[idp2];
	GFtracklen_decays[3] = GFL_pi_tracklen_container[idp2];
	GFmass2_decays[0] = GFL_p_mass2_container[idp1];
	GFmass2_decays[1] = GFL_pi_mass2_container[idp1];
	GFmass2_decays[2] = GFL_p_mass2_container[idp2];
	GFmass2_decays[3] = GFL_pi_mass2_container[idp2];

	//order : L1, L2
	GFLLid_container.push_back(GFlambda_containerid);
	GFLL_Ldecayvtx_container.push_back(GFlambda_vert);
	GFLLmom_container.push_back(GFlambda_mom);
	GFLLmass_container.push_back(GFlambda_mass);
	GFLLppidist_container.push_back(GFppi_dist);
	GFLL_Ltargetdist_container.push_back(GFlambdatgt_dist);
	GFLL_Ltargetvtx_container.push_back(GFlambdatgt_vtx);
	GFLL_Ltargetcenterdist_container.push_back(GFlambdatgtcenter_dist);
	GFLL_Ltargetcentervtx_container.push_back(GFlambdatgtcenter_vtx);

	//order : p, pi(L1), p, pi(L2)
	GFLLdecays_trackid_container.push_back(GFtrackid_decays);
	GFLLdecays_repid_container.push_back(GFrepid_decays);
	GFLLdecays_htofid_container.push_back(GFhtofid_decays);
	GFLLdecays_tof_container.push_back(GFtof_decays);
	GFLLdecays_mom_container.push_back(GFmom_decays);
	GFLLdecays_tracklen_container.push_back(GFtracklen_decays);
	GFLLdecays_mass2_container.push_back(GFmass2_decays);

	const Int_t ntrack_ll = 4;
#if KinematicFit
	event.KFdecays_mom.resize(ntrack_ll);
	event.KFdecays_mom_x.resize(ntrack_ll);
	event.KFdecays_mom_y.resize(ntrack_ll);
	event.KFdecays_mom_z.resize(ntrack_ll);
	event.KFdecays_CMmom.resize(ntrack_ll);
	event.KFdecays_CMmom_x.resize(ntrack_ll);
	event.KFdecays_CMmom_y.resize(ntrack_ll);
	event.KFdecays_CMmom_z.resize(ntrack_ll);

	int trackid_p1 = GFL_p_id_container[idp1];
	int trackid_pi1 = GFL_pi_id_container[idp1];
	auto track_p1 = TPCAna.GetTrackTPCHelix(trackid_p1);
	auto Vp1 = track_p1->GetCovarianceMatrix();
	auto track_pi1 = TPCAna.GetTrackTPCHelix(trackid_pi1);
	auto Vpi1 = track_pi1->GetCovarianceMatrix();
	double Diag_ppi1[6]={
	  Vp1(0,0),Vp1(1,1),Vp1(2,2),Vpi1(0,0),Vpi1(1,1),Vpi1(2,2)
	};
	auto Offdiag_ppi1 = MathTools::MergeOffdiagonals(Vp1, Vpi1);
	TVector3 HTVP1(GFL_p_mom_container[idp1].x(),
		       GFL_p_mom_container[idp1].z(),
		       GFL_p_mom_container[idp1].y());
	TVector3 HTVPi1(GFL_pi_mom_container[idp1].x(),
			GFL_pi_mom_container[idp1].z(),
			GFL_pi_mom_container[idp1].y());
	TVector3 HTVLd1 = HTVP1+HTVPi1;
	TLorentzVector HLVP1(HTVP1, TMath::Hypot(HTVP1.Mag(), pdg::ProtonMass()));
	TLorentzVector HLVPi1(HTVPi1, TMath::Hypot(HTVPi1.Mag(), pdg::PionMass()));
	TLorentzVector HLVLd1(HTVLd1, TMath::Hypot(HTVLd1.Mag(), pdg::LambdaMass()));
	Double_t KFchisqrl1=-1;
	Double_t KFpvall1=-1;
	FourVectorFitter KFLd1(HLVP1, HLVPi1, HLVLd1);
	KFLd1.SetInvMass(LambdaMass);
	KFLd1.SetMaximumStep(5);
	KFLd1.SetVariance(Diag_ppi1);
	KFLd1.AddOffdiagonals(Offdiag_ppi1);
	KFchisqrl1 = KFLd1.DoKinematicFit();
	KFpvall1 = KFLd1.GetPValue();
	auto HcontLd1 = KFLd1.GetFittedLV();
	auto PullLd1 = KFLd1.GetPull();
	auto KFHLVP1 = HcontLd1.at(0);
	auto KFHLVPi1 = HcontLd1.at(1);
	auto KFHLVLd1 = HcontLd1.at(2);
	auto KFTVP1 = TVector3(KFHLVP1.X(),KFHLVP1.Z(),KFHLVP1.Y());
	auto KFTVPi1 = TVector3(KFHLVPi1.X(),KFHLVPi1.Z(),KFHLVPi1.Y());
	auto KFTVLd1 = TVector3(KFHLVLd1.X(),KFHLVLd1.Z(),KFHLVLd1.Y());

	int trackid_p2 = GFL_p_id_container[idp2];
	int trackid_pi2 = GFL_pi_id_container[idp2];
	auto track_p2 = TPCAna.GetTrackTPCHelix(trackid_p2);
	auto Vp2 = track_p2->GetCovarianceMatrix();
	auto track_pi2 = TPCAna.GetTrackTPCHelix(trackid_pi2);
	auto Vpi2 = track_pi2->GetCovarianceMatrix();
	double Diag_ppi2[6]={
	  Vp2(0,0),Vp2(1,1),Vp2(2,2),Vpi2(0,0),Vpi2(1,1),Vpi2(2,2)
	};
	auto Offdiag_ppi2 = MathTools::MergeOffdiagonals(Vp2, Vpi2);
	TVector3 HTVP2(GFL_p_mom_container[idp2].x(),
		       GFL_p_mom_container[idp2].z(),
		       GFL_p_mom_container[idp2].y());
	TVector3 HTVPi2(GFL_pi_mom_container[idp2].x(),
			GFL_pi_mom_container[idp2].z(),
			GFL_pi_mom_container[idp2].y());
	TVector3 HTVLd2 = HTVP2+HTVPi2;
	TLorentzVector HLVP2(HTVP2, TMath::Hypot(HTVP2.Mag(), pdg::ProtonMass()));
	TLorentzVector HLVPi2(HTVPi2, TMath::Hypot(HTVPi2.Mag(), pdg::PionMass()));
	TLorentzVector HLVLd2(HTVLd2, TMath::Hypot(HTVLd2.Mag(), pdg::LambdaMass()));
	Double_t KFchisqrl2=-1;
	Double_t KFpvall2=-1;
	FourVectorFitter KFLd2(HLVP2, HLVPi2, HLVLd2);
	KFLd2.SetInvMass(LambdaMass);
	KFLd2.SetMaximumStep(5);
	KFLd2.SetVariance(Diag_ppi2);
	KFLd2.AddOffdiagonals(Offdiag_ppi2);
	KFchisqrl2 = KFLd2.DoKinematicFit();
	KFpvall2 = KFLd2.GetPValue();
	auto HcontLd2 = KFLd2.GetFittedLV();
	auto PullLd2 = KFLd2.GetPull();
	auto KFHLVP2 = HcontLd2.at(0);
	auto KFHLVPi2 = HcontLd2.at(1);
	auto KFHLVLd2 = HcontLd2.at(2);
	auto KFTVP2 = TVector3(KFHLVP2.X(),KFHLVP2.Z(),KFHLVP2.Y());
	auto KFTVPi2 = TVector3(KFHLVPi2.X(),KFHLVPi2.Z(),KFHLVPi2.Y());
	auto KFTVLd2 = TVector3(KFHLVLd2.X(),KFHLVLd2.Z(),KFHLVLd2.Y());

	Double_t KFll_dist;
	TVector3 KFll_vtx1, KFll_vtx2;
	TVector3 KFll_vtx
	  = Kinematics::LambdaLambdaVertex(GFL_vtx_container[idp1],
					   KFTVLd1,
					   GFL_vtx_container[idp2],
					   KFTVLd2,
					   KFll_vtx1, KFll_vtx2, KFll_dist);

	Double_t KFl_targetcenter_dist1;
	TVector3 KFl_pos_tgtcenter1 =
	  Kinematics::LambdaTargetCenter(GFL_vtx_container[idp1],
					 KFTVLd1,
					 KFl_targetcenter_dist1);
	Double_t KFl_targetcenter_dist2;
	TVector3 KFl_pos_tgtcenter2 =
	  Kinematics::LambdaTargetCenter(GFL_vtx_container[idp2],
					 KFTVLd2,
					 KFl_targetcenter_dist2);

	auto VLd1 = KFLd1.GetUnmeasuredCovariance();
	auto VLd2 = KFLd2.GetUnmeasuredCovariance();
	Double_t KFx0[ntrack_ll] =
	  {event.xtgtK18[0], event.xtgtTPCKurama[0],
	   KFl_pos_tgtcenter1.x(), KFl_pos_tgtcenter2.x()};
	Double_t KFy0[ntrack_ll] =
	  {event.ytgtK18[0], event.ytgtTPCKurama[0],
	   KFl_pos_tgtcenter1.y(), KFl_pos_tgtcenter2.y()};
	Double_t KFu0[ntrack_ll] =
	  {event.utgtK18[0], event.utgtTPCKurama[0],
	   KFTVLd1.x()/KFTVLd1.z(), KFTVLd2.x()/KFTVLd2.z()};
	Double_t KFv0[ntrack_ll] =
	  {event.vtgtK18[0], event.vtgtTPCKurama[0],
	   KFTVLd1.y()/KFTVLd1.z(), KFTVLd2.y()/KFTVLd2.z()};
	double resl1_u,resl1_v,resl2_u,resl2_v;
	MathTools::DecomposeResolutionUV(VLd1, KFTVLd1, resl1_u, resl1_v);
	MathTools::DecomposeResolutionUV(VLd2, KFTVLd2, resl2_u, resl2_v);
	std::vector<double> res_x =
	  {res_xKurama, res_xK18, res_xLdVtx, res_xLdVtx};
	std::vector<double> res_y =
	  {res_yKurama, res_yK18, res_yLdVtx, res_yLdVtx};
	std::vector<double> res_u =
	  {res_uKurama, res_uK18, resl1_u, resl2_u};
	std::vector<double> res_v =
	  {res_vKurama, res_vK18, resl1_v, resl2_v};

	Double_t chisqr_KFkk_ll_vertex = qnan;
	TVector3 KFkk_ll_vertex =
	  Kinematics::MultitrackVertex(ntrack_ll,
				       KFx0, KFy0,
				       KFu0, KFv0,
				       res_x, res_y,
				       res_u, res_v,
				       chisqr_KFkk_ll_vertex);
	Double_t KFprodvtx_closedist1 = qnan;

	TVector3 KFprodvtx_closest1 =
	  Kinematics::CalcCloseDistLambda(KFkk_ll_vertex,
					  GFL_vtx_container[idp1],
					  KFTVLd1,
					  KFprodvtx_closedist1);
	Double_t KFprodvtx_closedist2 = qnan;
	TVector3 KFprodvtx_closest2 =
	  Kinematics::CalcCloseDistLambda(KFkk_ll_vertex,
					  GFL_vtx_container[idp2],
					  KFTVLd2,
					  KFprodvtx_closedist2);

	Double_t KFx0_l1[3] = {event.xtgtK18[0], event.xtgtTPCKurama[0],
	  KFl_pos_tgtcenter1.x()};
	Double_t KFy0_l1[3] = {event.ytgtK18[0], event.ytgtTPCKurama[0],
	  KFl_pos_tgtcenter1.y()};
	Double_t KFu0_l1[3] = {event.utgtK18[0], event.utgtTPCKurama[0],
	  KFTVP1.x()/KFTVP1.z()};
	Double_t KFv0_l1[3] = {event.vtgtK18[0], event.vtgtTPCKurama[0],
	  KFTVP1.y()/KFTVP1.z()};
	std::vector<double> res_x_l1 = {res_xKurama, res_xK18, res_xLdVtx};
	std::vector<double> res_y_l1 = {res_yKurama, res_yK18, res_yLdVtx};
	std::vector<double> res_u_l1 = {res_uKurama, res_uK18, resl1_u};
	std::vector<double> res_v_l1 = {res_vKurama, res_vK18, resl1_v};
	TVector3 KFkk_ll_vertex_l1 = Kinematics::MultitrackVertex(3,
								  KFx0_l1, KFy0_l1,
								  KFu0_l1, KFv0_l1,
								  res_x_l1, res_y_l1,
								  res_u_l1, res_v_l1);

	Double_t KFx0_l2[3] = {event.xtgtK18[0], event.xtgtTPCKurama[0],
	  KFl_pos_tgtcenter2.x()};
	Double_t KFy0_l2[3] = {event.ytgtK18[0], event.ytgtTPCKurama[0],
	  KFl_pos_tgtcenter2.y()};
	Double_t KFu0_l2[3] = {event.utgtK18[0], event.utgtTPCKurama[0],
	  KFTVP2.x()/KFTVP2.z()};
	Double_t KFv0_l2[3] = {event.vtgtK18[0], event.vtgtTPCKurama[0],
	  KFTVP2.y()/KFTVP2.z()};
	std::vector<double> res_x_l2 = {res_xKurama, res_xK18, res_xLdVtx};
	std::vector<double> res_y_l2 = {res_yKurama, res_yK18, res_yLdVtx};
	std::vector<double> res_u_l2 = {res_uKurama, res_uK18, resl2_u};
	std::vector<double> res_v_l2 = {res_vKurama, res_vK18, resl2_v};
	TVector3 KFkk_ll_vertex_l2 = Kinematics::MultitrackVertex(3,
								  KFx0_l2, KFy0_l2,
								  KFu0_l2, KFv0_l2,
								  res_x_l2, res_y_l2,
								  res_u_l2, res_v_l2);

	/*
	 */
	//Double_t diff = TMath::Hypot(GFlambda_mass[0] - LambdaMass, GFlambda_mass[1] - LambdaMass);
	//Double_t diff = TMath::Hypot(L_mass_container[idp1] - LambdaMass, L_mass_container[idp2] - LambdaMass);
	if(prev_LLfit_chisqr > chisqr_KFkk_ll_vertex){
	  event.llflag = true;
	  prev_LLfit_chisqr = chisqr_KFkk_ll_vertex;
	  LL_best = LLcount;
	}
#else
	Double_t GFll_dist;
	TVector3 GFll_vtx1, GFll_vtx2;
	TVector3 GFll_vtx
	  = Kinematics::LambdaLambdaVertex(GFL_vtx_container.at(idp1),
					   GFL_mom_container.at(idp1),
					   GFL_vtx_container.at(idp2),
					   GFL_mom_container.at(idp2),
					   GFll_vtx1, GFll_vtx2, GFll_dist);
	Double_t diff = GFll_dist;
	if(prev_LLfit_chisqr > diff){
	  event.llflag = true;
	  prev_LLfit_chisqr = diff;
	  LL_best = LLcount;
	}
#endif
	LLcount++;
      } //L combi1
    } //L combi2

    //for the LL event
    if(event.llflag){
#if DebugDisp
    std::cout<<"3. Calculate and save the best LL pair"<<std::endl;
#endif
      Int_t L1_id = GFLLid_container[LL_best][0];
      Int_t L2_id = GFLLid_container[LL_best][1];

      event.lmass1 = L_mass_container.at(L1_id);
      event.ldecayvtx_x1 = L_vtx_container.at(L1_id).x();
      event.ldecayvtx_y1 = L_vtx_container.at(L1_id).y();
      event.ldecayvtx_z1 = L_vtx_container.at(L1_id).z();
      event.lmom1 = L_mom_container.at(L1_id).Mag();
      event.lmom_x1 = L_mom_container.at(L1_id).x();
      event.lmom_y1 = L_mom_container.at(L1_id).y();
      event.lmom_z1 = L_mom_container.at(L1_id).z();
      event.ppi_dist1 = L_ppidist_container.at(L1_id);
      event.ltarget_dist1 = L_targetdist_container.at(L1_id);
      event.ltargetvtx_x1 = L_targetvtx_container.at(L1_id).x();
      event.ltargetvtx_y1 = L_targetvtx_container.at(L1_id).y();
      event.ltargetvtx_z1 = L_targetvtx_container.at(L1_id).z();

      event.lmass2 = L_mass_container.at(L2_id);
      event.ldecayvtx_x2 = L_vtx_container.at(L2_id).x();
      event.ldecayvtx_y2 = L_vtx_container.at(L2_id).y();
      event.ldecayvtx_z2 = L_vtx_container.at(L2_id).z();
      event.lmom2 = L_mom_container.at(L2_id).Mag();
      event.lmom_x2 = L_mom_container.at(L2_id).x();
      event.lmom_y2 = L_mom_container.at(L2_id).y();
      event.lmom_z2 = L_mom_container.at(L2_id).z();
      event.ppi_dist2 = L_ppidist_container.at(L2_id);
      event.ltarget_dist2 = L_targetdist_container.at(L2_id);
      event.ltargetvtx_x2 = L_targetvtx_container.at(L2_id).x();
      event.ltargetvtx_y2 = L_targetvtx_container.at(L2_id).y();
      event.ltargetvtx_z2 = L_targetvtx_container.at(L2_id).z();

      event.GFlmass1 = GFLLmass_container[LL_best][0];
      event.GFldecayvtx_x1 = GFLL_Ldecayvtx_container[LL_best][0].x();
      event.GFldecayvtx_y1 = GFLL_Ldecayvtx_container[LL_best][0].y();
      event.GFldecayvtx_z1 = GFLL_Ldecayvtx_container[LL_best][0].z();
      event.GFlmom1 = GFLLmom_container[LL_best][0].Mag();
      event.GFlmom_x1 = GFLLmom_container[LL_best][0].x();
      event.GFlmom_y1 = GFLLmom_container[LL_best][0].y();
      event.GFlmom_z1 = GFLLmom_container[LL_best][0].z();
      event.GFppi_dist1 = GFLLppidist_container[LL_best][0];
      event.GFltarget_dist1 = GFLL_Ltargetdist_container[LL_best][0];
      event.GFltargetvtx_x1 = GFLL_Ltargetvtx_container[LL_best][0].x();
      event.GFltargetvtx_y1 = GFLL_Ltargetvtx_container[LL_best][0].y();
      event.GFltargetvtx_z1 = GFLL_Ltargetvtx_container[LL_best][0].z();
      event.GFltargetcenter_dist1 = GFLL_Ltargetcenterdist_container[LL_best][0];
      event.GFltargetcenter_x1 = GFLL_Ltargetcentervtx_container[LL_best][0].x();
      event.GFltargetcenter_y1 = GFLL_Ltargetcentervtx_container[LL_best][0].y();
      event.GFltargetcenter_z1 = GFLL_Ltargetcentervtx_container[LL_best][0].z();

      event.GFlmass2 = GFLLmass_container[LL_best][1];
      event.GFldecayvtx_x2 = GFLL_Ldecayvtx_container[LL_best][1].x();
      event.GFldecayvtx_y2 = GFLL_Ldecayvtx_container[LL_best][1].y();
      event.GFldecayvtx_z2 = GFLL_Ldecayvtx_container[LL_best][1].z();
      event.GFlmom2 = GFLLmom_container[LL_best][1].Mag();
      event.GFlmom_x2 = GFLLmom_container[LL_best][1].x();
      event.GFlmom_y2 = GFLLmom_container[LL_best][1].y();
      event.GFlmom_z2 = GFLLmom_container[LL_best][1].z();
      event.GFppi_dist2 = GFLLppidist_container[LL_best][1];
      event.GFltarget_dist2 = GFLL_Ltargetdist_container[LL_best][1];
      event.GFltargetvtx_x2 = GFLL_Ltargetvtx_container[LL_best][1].x();
      event.GFltargetvtx_y2 = GFLL_Ltargetvtx_container[LL_best][1].y();
      event.GFltargetvtx_z2 = GFLL_Ltargetvtx_container[LL_best][1].z();
      event.GFltargetcenter_dist2 = GFLL_Ltargetcenterdist_container[LL_best][1];
      event.GFltargetcenter_x2 = GFLL_Ltargetcentervtx_container[LL_best][1].x();
      event.GFltargetcenter_y2 = GFLL_Ltargetcentervtx_container[LL_best][1].y();
      event.GFltargetcenter_z2 = GFLL_Ltargetcentervtx_container[LL_best][1].z();

      Double_t ll_dist;
      TVector3 ll_vtx1, ll_vtx2;
      TVector3 ll_vtx
	= Kinematics::LambdaLambdaVertex(L_vtx_container.at(L1_id),
					 L_mom_container.at(L1_id),
					 L_vtx_container.at(L2_id),
					 L_mom_container.at(L2_id),
					 ll_vtx1, ll_vtx2, ll_dist);
      event.llvtx_x = ll_vtx.x();
      event.llvtx_y = ll_vtx.y();
      event.llvtx_z = ll_vtx.z();
      event.lldist = ll_dist;

      Double_t GFll_dist;
      TVector3 GFll_vtx1, GFll_vtx2;
      TVector3 GFll_vtx
	= Kinematics::LambdaLambdaVertex(GFLL_Ldecayvtx_container[LL_best][0],
					 GFLLmom_container[LL_best][0],
					 GFLL_Ldecayvtx_container[LL_best][1],
					 GFLLmom_container[LL_best][1],
					 GFll_vtx1, GFll_vtx2, GFll_dist);
      event.GFllvtx_x = GFll_vtx.x();
      event.GFllvtx_y = GFll_vtx.y();
      event.GFllvtx_z = GFll_vtx.z();
      event.GFlldist = GFll_dist;

      const Int_t ntrack_ll = 4;
      Double_t x0[ntrack_ll] = {event.xtgtK18[0], event.xtgtTPCKurama[0],
				event.GFltargetcenter_x1, event.GFltargetcenter_x2};
      Double_t y0[ntrack_ll] = {event.ytgtK18[0], event.ytgtTPCKurama[0],
				event.GFltargetcenter_y1, event.GFltargetcenter_y2};
      Double_t u0[ntrack_ll] = {event.utgtK18[0], event.utgtTPCKurama[0],
				event.GFlmom_x1/event.GFlmom_z1,
				event.GFlmom_x2/event.GFlmom_z2};
      Double_t v0[ntrack_ll] = {event.vtgtK18[0], event.vtgtTPCKurama[0],
				event.GFlmom_y1/event.GFlmom_z1,
				event.GFlmom_y2/event.GFlmom_z2};
      TVector3 GFkk_ll_vertex = Kinematics::MultitrackVertex(ntrack_ll, x0, y0, u0, v0);

      Double_t GFprodvtx_closedist1 = qnan;
      TVector3 GFprodvtx_closest1 = Kinematics::CalcCloseDistLambda(GFkk_ll_vertex,
								    GFLL_Ldecayvtx_container[LL_best][0],
								    GFLLmom_container[LL_best][0],
								    GFprodvtx_closedist1);

      Double_t GFprodvtx_closedist2 = qnan;
      TVector3 GFprodvtx_closest2 = Kinematics::CalcCloseDistLambda(GFkk_ll_vertex,
								    GFLL_Ldecayvtx_container[LL_best][1],
								    GFLLmom_container[LL_best][1],
								    GFprodvtx_closedist2);

    int id_p1 = L_p_id_container.at(L1_id);
    int id_pi1 = L_pi_id_container.at(L1_id);
    int id_p2 = L_p_id_container.at(L2_id);
    int id_pi2 = L_pi_id_container.at(L2_id);
    int G4p1tid = G4TrackID.at(id_p1);
    int G4pi1tid = G4TrackID.at(id_pi1);
    int G4p2tid = G4TrackID.at(id_p2);
    int G4pi2tid = G4TrackID.at(id_pi2);

    if(G4p1tid == event.G4p1id and G4p2tid == event.G4p2id ){
      if(G4pi1tid == event.G4pi1id){
        event.l1good = true;
      }
      if(G4pi2tid == event.G4pi2id){
        event.l2good = true;
      }
      if(G4pi1tid == event.G4pi2id and G4pi2tid == event.G4pi1id ){
        event.piswap = true;
      }

    }
    else if(G4p1tid == event.G4p2id and G4p2tid == event.G4p1id ){
      event.llswap = true;
      if(G4pi1tid == event.G4pi2id){
	event.l1good = true;
      }
      if(G4pi2tid == event.G4pi1id){
	event.l2good = true;
      }
      if(G4pi1tid == event.G4p1id and G4pi2tid == event.G4p2id ){
	event.piswap = true;
      }
    }

    event.p1tid = id_p1;
    event.p2tid = id_p2;
    event.pi1tid = id_pi1;
    event.pi2tid = id_pi2;
    event.G4p1tid = G4p1tid;
    event.G4p2tid = G4p2tid;
    event.G4pi1tid = G4pi1tid;
    event.G4pi2tid = G4pi2tid;

    event.p1nh = event.nhtrack.at(id_p1);
    event.p2nh = event.nhtrack.at(id_p2);
    event.pi1nh = event.nhtrack.at(id_pi1);
    event.pi2nh = event.nhtrack.at(id_pi2);

    event.p1mom = L_p_mom_container.at(L1_id).Mag();
    event.p1mom_x = L_p_mom_container.at(L1_id).x();
    event.p1mom_y = L_p_mom_container.at(L1_id).y();
    event.p1mom_z = L_p_mom_container.at(L1_id).z();
    event.pi1mom = L_pi_mom_container.at(L1_id).Mag();
    event.pi1mom_x = L_pi_mom_container.at(L1_id).x();
    event.pi1mom_y = L_pi_mom_container.at(L1_id).y();
    event.pi1mom_z = L_pi_mom_container.at(L1_id).z();
    event.p2mom = L_p_mom_container.at(L2_id).Mag();
    event.p2mom_x = L_p_mom_container.at(L2_id).x();
    event.p2mom_y = L_p_mom_container.at(L2_id).y();
    event.p2mom_z = L_p_mom_container.at(L2_id).z();
    event.pi2mom = L_pi_mom_container.at(L2_id).Mag();
    event.pi2mom_x = L_pi_mom_container.at(L2_id).x();
    event.pi2mom_y = L_pi_mom_container.at(L2_id).y();
    event.pi2mom_z = L_pi_mom_container.at(L2_id).z();

#if DebugDisp
      std::cout<<"K-K+LL vertex "<<GFkk_ll_vertex
	       <<" Lambda1's closest point to the vertex"<<GFprodvtx_closest1
	       <<" Lambda2's closest point to the vertex"<<GFprodvtx_closest2
	       <<std::endl;
#endif

      event.GFprodvtx_x_ll = GFkk_ll_vertex.x();
      event.GFprodvtx_y_ll = GFkk_ll_vertex.y();
      event.GFprodvtx_z_ll = GFkk_ll_vertex.z();

      Double_t x0_l1[3] = {event.xtgtK18[0], event.xtgtTPCKurama[0],
			   event.GFltargetcenter_x1};
      Double_t y0_l1[3] = {event.ytgtK18[0], event.ytgtTPCKurama[0],
			   event.GFltargetcenter_y1};
      Double_t u0_l1[3] = {event.utgtK18[0], event.utgtTPCKurama[0],
			   event.GFlmom_x1/event.GFlmom_z1};
      Double_t v0_l1[3] = {event.vtgtK18[0], event.vtgtTPCKurama[0],
			   event.GFlmom_y1/event.GFlmom_z1};
      TVector3 GFprodvtxkk_ll_vertex_l1 = Kinematics::MultitrackVertex(3, x0_l1, y0_l1, u0_l1, v0_l1);

      Double_t x0_l2[3] = {event.xtgtK18[0], event.xtgtTPCKurama[0],
			   event.GFltargetcenter_x2};
      Double_t y0_l2[3] = {event.ytgtK18[0], event.ytgtTPCKurama[0],
			   event.GFltargetcenter_y2};
      Double_t u0_l2[3] = {event.utgtK18[0], event.utgtTPCKurama[0],
			   event.GFlmom_x2/event.GFlmom_z2};
      Double_t v0_l2[3] = {event.vtgtK18[0], event.vtgtTPCKurama[0],
			   event.GFlmom_y2/event.GFlmom_z2};
      TVector3 GFprodvtxkk_ll_vertex_l2 = Kinematics::MultitrackVertex(3, x0_l2, y0_l2, u0_l2, v0_l2);

      event.GFprodvtx_x_l1 = GFprodvtxkk_ll_vertex_l1.x();
      event.GFprodvtx_y_l1 = GFprodvtxkk_ll_vertex_l1.y();
      event.GFprodvtx_z_l1 = GFprodvtxkk_ll_vertex_l1.z();
      event.GFprodvtx_x_l2 = GFprodvtxkk_ll_vertex_l2.x();
      event.GFprodvtx_y_l2 = GFprodvtxkk_ll_vertex_l2.y();
      event.GFprodvtx_z_l2 = GFprodvtxkk_ll_vertex_l2.z();

      TLorentzVector LvL1_fixedmass(GFLLmom_container[LL_best][0],
				    TMath::Hypot(GFLLmom_container[LL_best][0].Mag(), LambdaMass));
      TLorentzVector LvL2_fixedmass(GFLLmom_container[LL_best][1],
				    TMath::Hypot(GFLLmom_container[LL_best][1].Mag(), LambdaMass));
      TLorentzVector LvRcLL_fixedmass = LvRcTPC - LvL1_fixedmass - LvL2_fixedmass;
      event.GFllexcitation = LvRcLL_fixedmass.M() - m10Be;

      TVector3 lambda_tracklen = GFkk_ll_vertex - GFLL_Ldecayvtx_container[LL_best][0];
      event.GFltracklen1 = lambda_tracklen.Mag();
      event.GFltof1 = Kinematics::CalcTimeOfFlight(event.GFlmom1,
						   lambda_tracklen.Mag(),
						   pdg::LambdaMass());

      lambda_tracklen = GFkk_ll_vertex - GFLL_Ldecayvtx_container[LL_best][1];
      event.GFltracklen2 = lambda_tracklen.Mag();
      event.GFltof2 = Kinematics::CalcTimeOfFlight(event.GFlmom2,
						   lambda_tracklen.Mag(),
						   pdg::LambdaMass());

      //p, pi(L1), p, pi(L2) tracks
      event.GFntdecays = 4;
      event.GFnhtrack.resize(event.GFntdecays);
      event.GFchisqr.resize(event.GFntdecays);
      event.GFcharge.resize(event.GFntdecays);
      event.GFtof.resize(event.GFntdecays);
      event.GFtracklen.resize(event.GFntdecays);
      event.GFpval.resize(event.GFntdecays);
      event.GFpdgcode.resize(event.GFntdecays);

      event.GFdecays_htofid.resize(event.GFntdecays);
      event.GFdecays_tracklen.resize(event.GFntdecays);
      event.GFdecays_tof.resize(event.GFntdecays);
      event.GFdecays_mass2.resize(event.GFntdecays);
      event.GFdecays_mom.resize(event.GFntdecays);
      event.GFdecays_mom_x.resize(event.GFntdecays);
      event.GFdecays_mom_y.resize(event.GFntdecays);
      event.GFdecays_mom_z.resize(event.GFntdecays);
      event.GFdecays_CMmom.resize(event.GFntdecays);
      event.GFdecays_CMmom_x.resize(event.GFntdecays);
      event.GFdecays_CMmom_y.resize(event.GFntdecays);
      event.GFdecays_CMmom_z.resize(event.GFntdecays);
      event.GFmomloss.resize(event.GFntdecays);
      event.GFeloss.resize(event.GFntdecays);

      event.decays_id.push_back(GFLLdecays_trackid_container[LL_best][0]);
      event.decays_id.push_back(GFLLdecays_trackid_container[LL_best][1]);
      event.decays_id.push_back(GFLLdecays_trackid_container[LL_best][2]);
      event.decays_id.push_back(GFLLdecays_trackid_container[LL_best][3]);

      event.decays_mom.push_back(L_p_mom_container[L1_id].Mag());
      event.decays_mom.push_back(L_pi_mom_container[L1_id].Mag());
      event.decays_mom.push_back(L_p_mom_container[L2_id].Mag());
      event.decays_mom.push_back(L_pi_mom_container[L2_id].Mag());
      event.decays_mom_x.push_back(L_p_mom_container[L1_id].x());
      event.decays_mom_x.push_back(L_pi_mom_container[L1_id].x());
      event.decays_mom_x.push_back(L_p_mom_container[L2_id].x());
      event.decays_mom_x.push_back(L_pi_mom_container[L2_id].x());
      event.decays_mom_y.push_back(L_p_mom_container[L1_id].y());
      event.decays_mom_y.push_back(L_pi_mom_container[L1_id].y());
      event.decays_mom_y.push_back(L_p_mom_container[L2_id].y());
      event.decays_mom_y.push_back(L_pi_mom_container[L2_id].y());
      event.decays_mom_z.push_back(L_p_mom_container[L1_id].z());
      event.decays_mom_z.push_back(L_pi_mom_container[L1_id].z());
      event.decays_mom_z.push_back(L_p_mom_container[L2_id].z());
      event.decays_mom_z.push_back(L_pi_mom_container[L2_id].z());

      for(int j=0;j<event.GFntdecays;j++){
	Int_t igf = GFLLdecays_trackid_container[LL_best][j];
	Int_t repid = GFLLdecays_repid_container[LL_best][j];
	event.GFnhtrack[j] = GFTrackCont.GetNHits(igf);
	event.GFchisqr[j] = GFTrackCont.GetChi2NDF(igf, repid);
	event.GFcharge[j] = GFTrackCont.GetCharge(igf, repid);
	event.GFtof[j] = GFTrackCont.GetTrackTOF(igf, 0, -1, repid);
	event.GFtracklen[j] = GFTrackCont.GetTrackLength(igf, 0, -1, repid);
	event.GFpval[j] = GFTrackCont.GetPvalue(igf, repid);
	event.GFpdgcode[j] = GFTrackCont.GetPDGcode(igf, repid);

	event.GFdecays_htofid[j] = GFLLdecays_htofid_container[LL_best][j];
	event.GFdecays_tracklen[j] = GFLLdecays_tracklen_container[LL_best][j];
	event.GFdecays_tof[j] = GFLLdecays_tof_container[LL_best][j];
	event.GFdecays_mass2[j] = GFLLdecays_mass2_container[LL_best][j];
	event.GFdecays_mom[j] = GFLLdecays_mom_container[LL_best][j].Mag();
	event.GFdecays_mom_x[j] = GFLLdecays_mom_container[LL_best][j].x();
	event.GFdecays_mom_y[j] = GFLLdecays_mom_container[LL_best][j].y();
	event.GFdecays_mom_z[j] = GFLLdecays_mom_container[LL_best][j].z();
	event.GFmomloss[j] = qnan;
	event.GFeloss[j] = qnan;

	Double_t mass = PionMass;
	if(j==0 || j==2) mass = ProtonMass;

	TLorentzVector GFLv_decays(GFLLdecays_mom_container[LL_best][j],
				   TMath::Hypot(event.GFdecays_mom[j], mass));
	if(j<2) GFLv_decays.Boost(-GFLLmom_container[LL_best][0]);
	else GFLv_decays.Boost(-GFLLmom_container[LL_best][1]);

	event.GFdecays_CMmom[j] = GFLv_decays.P();
	event.GFdecays_CMmom_x[j] = GFLv_decays.Px();
	event.GFdecays_CMmom_y[j] = GFLv_decays.Py();
	event.GFdecays_CMmom_z[j] = GFLv_decays.Pz();

      } //GFntdecays

      HF2( 502, event.lmass1, event.lmass2);
      HF2( 503, event.GFlmass1, event.GFlmass2);
      HF1( 504, event.lmass1);
      HF1( 504, event.lmass2);
      HF1( 505, event.GFlmass1);
      HF1( 505, event.GFlmass2);
#if KinematicFit
      TVector3 KFTVLd1(event.KFlmom_x1, event.KFlmom_y1, event.KFlmom_z1);
      TVector3 KFTVLd2(event.KFlmom_x2, event.KFlmom_y2, event.KFlmom_z2);
      TVector3 KFTVP1(event.KFdecays_mom_x[0], event.KFdecays_mom_y[0], event.KFdecays_mom_z[0]);
      TVector3 KFTVPi1(event.KFdecays_mom_x[1], event.KFdecays_mom_y[1], event.KFdecays_mom_z[1]);
      TVector3 KFTVP2(event.KFdecays_mom_x[2], event.KFdecays_mom_y[2], event.KFdecays_mom_z[2]);
      TVector3 KFTVPi2(event.KFdecays_mom_x[3], event.KFdecays_mom_y[3], event.KFdecays_mom_z[3]);
      TVector3 KFkk_ll_vertex(event.KFprodvtx_x_ll, event.KFprodvtx_y_ll, event.KFprodvtx_z_ll);

      TLorentzVector KFLvL1_fixedmass(KFTVLd1,
				      TMath::Hypot(KFTVLd1.Mag(), LambdaMass));
      TLorentzVector KFLvL2_fixedmass(KFTVLd2,
				      TMath::Hypot(KFTVLd2.Mag(), LambdaMass));
      TLorentzVector KFLvRcLL_fixedmass = LvRcTPC - KFLvL1_fixedmass - KFLvL2_fixedmass;
      event.KFllexcitation = KFLvRcLL_fixedmass.M() - m10Be;

      TVector3 KFlambda_tracklen = KFkk_ll_vertex - GFLL_Ldecayvtx_container[LL_best][0];
      event.KFltracklen1 = KFlambda_tracklen.Mag();
      event.KFltof1 = Kinematics::CalcTimeOfFlight(KFTVLd1.Mag(),
						   KFlambda_tracklen.Mag(),
						   pdg::LambdaMass());

      KFlambda_tracklen = KFkk_ll_vertex - GFLL_Ldecayvtx_container[LL_best][1];
      event.KFltracklen2 = KFlambda_tracklen.Mag();
      event.KFltof2 = Kinematics::CalcTimeOfFlight(KFTVLd2.Mag(),
						   KFlambda_tracklen.Mag(),
						   pdg::LambdaMass());

      TLorentzVector KFLv_p1(KFTVP1, TMath::Hypot(KFTVP1.Mag(), ProtonMass));
      TLorentzVector KFLv_pi1(KFTVPi1, TMath::Hypot(KFTVPi1.Mag(), PionMass));
      TLorentzVector KFLv_p2(KFTVP2, TMath::Hypot(KFTVP2.Mag(), ProtonMass));
      TLorentzVector KFLv_pi2(KFTVPi2, TMath::Hypot(KFTVPi2.Mag(), PionMass));
      KFLv_p1.Boost(-KFTVLd1);
      KFLv_pi1.Boost(-KFTVLd1);
      KFLv_p2.Boost(-KFTVLd2);
      KFLv_pi2.Boost(-KFTVLd2);

      event.KFdecays_CMmom[0] = KFLv_p1.P();
      event.KFdecays_CMmom_x[0] = KFLv_p1.Px();
      event.KFdecays_CMmom_y[0] = KFLv_p1.Py();
      event.KFdecays_CMmom_z[0] = KFLv_p1.Pz();
      event.KFdecays_CMmom[1] = KFLv_pi1.P();
      event.KFdecays_CMmom_x[1] = KFLv_pi1.Px();
      event.KFdecays_CMmom_y[1] = KFLv_pi1.Py();
      event.KFdecays_CMmom_z[1] = KFLv_pi1.Pz();
      event.KFdecays_CMmom[2] = KFLv_p2.P();
      event.KFdecays_CMmom_x[2] = KFLv_p2.Px();
      event.KFdecays_CMmom_y[2] = KFLv_p2.Py();
      event.KFdecays_CMmom_z[2] = KFLv_p2.Pz();
      event.KFdecays_CMmom[3] = KFLv_pi2.P();
      event.KFdecays_CMmom_x[3] = KFLv_pi2.Px();
      event.KFdecays_CMmom_y[3] = KFLv_pi2.Py();
      event.KFdecays_CMmom_z[3] = KFLv_pi2.Pz();
#endif
      //Remaining tracks
#if DebugDisp
      std::cout<<"LL event, Save remaining tracks"<<std::endl;
#endif

      std::vector<Int_t> more_ppip_residuals;
      std::vector<Int_t> more_pim_residuals;
      std::vector<Int_t> more_pip_residuals;

      for(Int_t icombi=0; icombi<reconfailed_L_p_id_container.size(); ++icombi){
	Int_t p_id = reconfailed_L_p_id_container[icombi];
	Int_t pi_id = reconfailed_L_pi_id_container[icombi];

	if(p_id!=event.decays_id[0] && p_id!=event.decays_id[1] &&
	   p_id!=event.decays_id[2] && p_id!=event.decays_id[3]){
	  if(event.charge[p_id]==1){
	    if(std::find(target_p_id_container.begin(), target_p_id_container.end(), p_id)
	       == target_p_id_container.end() &&
	       std::find(target_ppip_id_container.begin(), target_ppip_id_container.end(), p_id)
	       == target_ppip_id_container.end()) more_ppip_residuals.push_back(p_id);
	  }
	  else{
	    if(std::find(target_pim_id_container.begin(), target_pim_id_container.end(), p_id)
	       == target_pim_id_container.end()) more_pim_residuals.push_back(p_id);
	  }
	}

	if(pi_id!=event.decays_id[0] && pi_id!=event.decays_id[1] &&
	   pi_id!=event.decays_id[2] && pi_id!=event.decays_id[3]){
	  if(event.charge[pi_id]==-1){
	    if(std::find(target_pim_id_container.begin(), target_pim_id_container.end(), pi_id)
	       == target_pim_id_container.end()) more_pim_residuals.push_back(pi_id);
	  }
	  else{
	    if(std::find(target_p_id_container.begin(), target_p_id_container.end(), pi_id)
	       == target_p_id_container.end() &&
	       std::find(target_ppip_id_container.begin(), target_ppip_id_container.end(), pi_id)
	       == target_ppip_id_container.end()) more_ppip_residuals.push_back(pi_id);
	  }
	}
      }

      std::vector<Int_t> delete_combi;
      for(Int_t icombi=0; icombi<reconfailed_L_p_id_container.size(); ++icombi){
	Int_t p_id = reconfailed_L_p_id_container[icombi];
	Int_t pi_id = reconfailed_L_pi_id_container[icombi];
	if(p_id==event.decays_id[0] || p_id==event.decays_id[1] ||
	   p_id==event.decays_id[2] || p_id==event.decays_id[3] ||
	   pi_id==event.decays_id[0] || pi_id==event.decays_id[1] ||
	   pi_id==event.decays_id[2] || pi_id==event.decays_id[3]) delete_combi.push_back(icombi);
      }

      for(Int_t icombi=0; icombi<delete_combi.size(); ++icombi){
	reconfailed_L_p_id_container.erase(reconfailed_L_p_id_container.begin()+delete_combi[icombi]-icombi);
	reconfailed_L_pi_id_container.erase(reconfailed_L_pi_id_container.begin()+delete_combi[icombi]-icombi);
	reconfailed_L_mass_container.erase(reconfailed_L_mass_container.begin()+delete_combi[icombi]-icombi);
	reconfailed_L_mom_container.erase(reconfailed_L_mom_container.begin()+delete_combi[icombi]-icombi);
	reconfailed_L_vtx_container.erase(reconfailed_L_vtx_container.begin()+delete_combi[icombi]-icombi);
	reconfailed_L_p_mom_container.erase(reconfailed_L_p_mom_container.begin()+delete_combi[icombi]-icombi);
	reconfailed_L_pi_mom_container.erase(reconfailed_L_pi_mom_container.begin()+delete_combi[icombi]-icombi);
	reconfailed_L_ppidist_container.erase(reconfailed_L_ppidist_container.begin()+delete_combi[icombi]-icombi);
      }

      event.ncombiLreconfailed = reconfailed_L_p_id_container.size();
      for(Int_t icombi=0; icombi<event.ncombiLreconfailed; ++icombi){
	event.pidLreconfailed.push_back(reconfailed_L_p_id_container[icombi]);
	event.piidLreconfailed.push_back(reconfailed_L_pi_id_container[icombi]);
	event.LmassLreconfailed.push_back(reconfailed_L_mass_container[icombi]);
	event.LdecayvtxLreconfailed_x.push_back(reconfailed_L_vtx_container[icombi].x());
	event.LdecayvtxLreconfailed_y.push_back(reconfailed_L_vtx_container[icombi].y());
	event.LdecayvtxLreconfailed_z.push_back(reconfailed_L_vtx_container[icombi].z());
	event.LmomLreconfailed.push_back(reconfailed_L_pi_mom_container[icombi].Mag());
	event.LmomLreconfailed.push_back(reconfailed_L_pi_mom_container[icombi].x());
	event.LmomLreconfailed.push_back(reconfailed_L_pi_mom_container[icombi].y());
	event.LmomLreconfailed.push_back(reconfailed_L_pi_mom_container[icombi].z());
	event.pmomLreconfailed.push_back(reconfailed_L_p_mom_container[icombi].Mag());
	event.pmomLreconfailed.push_back(reconfailed_L_p_mom_container[icombi].x());
	event.pmomLreconfailed.push_back(reconfailed_L_p_mom_container[icombi].y());
	event.pmomLreconfailed.push_back(reconfailed_L_p_mom_container[icombi].z());
	event.pimomLreconfailed.push_back(reconfailed_L_mom_container[icombi].Mag());
	event.pimomLreconfailed.push_back(reconfailed_L_mom_container[icombi].x());
	event.pimomLreconfailed.push_back(reconfailed_L_mom_container[icombi].y());
	event.pimomLreconfailed.push_back(reconfailed_L_mom_container[icombi].z());
	event.ppidistLreconfailed.push_back(reconfailed_L_ppidist_container[icombi]);
      }

      for(Int_t icombi=0; icombi<pipair_pip_id_container.size(); ++icombi){
	Int_t pip_id = pipair_pip_id_container[icombi];
	Int_t pim_id = pipair_pim_id_container[icombi];
	if(pip_id!=event.decays_id[0] && pip_id!=event.decays_id[2]){
	  if(std::find(target_pip_id_container.begin(), target_pip_id_container.end(), pip_id) ==
	     target_pip_id_container.end()) more_pip_residuals.push_back(pip_id);
	}
	if(pim_id!=event.decays_id[1] && pim_id!=event.decays_id[3]){
	  if(std::find(target_pim_id_container.begin(), target_pim_id_container.end(), pim_id) ==
	     target_pim_id_container.end()) more_pim_residuals.push_back(pim_id);
	}
      }

      delete_combi.clear();
      for(Int_t icombi=0; icombi<pipair_pip_id_container.size(); ++icombi){
	Int_t pip_id = pipair_pip_id_container[icombi];
	Int_t pim_id = pipair_pim_id_container[icombi];
	if(pip_id==event.decays_id[0] || pip_id==event.decays_id[2] ||
	   pim_id==event.decays_id[1] || pim_id==event.decays_id[3]) delete_combi.push_back(icombi);
      }

      for(Int_t icombi=0; icombi<delete_combi.size(); ++icombi){
	pipair_pip_id_container.erase(pipair_pip_id_container.begin()+delete_combi[icombi]-icombi);
	pipair_pim_id_container.erase(pipair_pim_id_container.begin()+delete_combi[icombi]-icombi);
	pipair_reconL_mass_container.erase(pipair_reconL_mass_container.begin()+delete_combi[icombi]-icombi);
	pipair_mom_container.erase(pipair_mom_container.begin()+delete_combi[icombi]-icombi);
	pipair_pip_mom_container.erase(pipair_pip_mom_container.begin()+delete_combi[icombi]-icombi);
	pipair_pim_mom_container.erase(pipair_pim_mom_container.begin()+delete_combi[icombi]-icombi);
	pipair_pipidist_container.erase(pipair_pipidist_container.begin()+delete_combi[icombi]-icombi);
      }

      event.ncombiPipair = pipair_pip_id_container.size();
      for(Int_t icombi=0; icombi<event.ncombiPipair; ++icombi){
	event.pipidPipair.push_back(pipair_pip_id_container[icombi]);
	event.pimidPipair.push_back(pipair_pim_id_container[icombi]);
	event.pipmomPipair.push_back(pipair_pip_mom_container[icombi].Mag());
	event.pipmomPipair_x.push_back(pipair_pip_mom_container[icombi].x());
	event.pipmomPipair_y.push_back(pipair_pip_mom_container[icombi].y());
	event.pipmomPipair_z.push_back(pipair_pip_mom_container[icombi].z());
	event.pimmomPipair.push_back(pipair_pim_mom_container[icombi].Mag());
	event.pimmomPipair_x.push_back(pipair_pim_mom_container[icombi].x());
	event.pimmomPipair_y.push_back(pipair_pim_mom_container[icombi].y());
	event.pimmomPipair_z.push_back(pipair_pim_mom_container[icombi].z());
	event.momPipair.push_back(pipair_mom_container[icombi].Mag());
	event.momPipair_x.push_back(pipair_mom_container[icombi].x());
	event.momPipair_y.push_back(pipair_mom_container[icombi].y());
	event.momPipair_z.push_back(pipair_mom_container[icombi].z());
	event.reconLmassPipair.push_back(pipair_reconL_mass_container[icombi]);
	event.pipidistPipair.push_back(pipair_pipidist_container[icombi]);
      }

      for(Int_t icombi=0; icombi<more_ppip_residuals.size(); ++icombi){
	Int_t id = more_ppip_residuals[icombi];
	if(std::find(target_p_id_container.begin(), target_p_id_container.end(), id) == target_p_id_container.end() &&
	   std::find(target_ppip_id_container.begin(), target_ppip_id_container.end(), id) == target_ppip_id_container.end()){
	  TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( id );
	  if((event.pid[id]&4)==4 && (event.pid[id]&1)!=1 && event.charge[id]==1){ //proton

	    Double_t mass2 = qnan;
	    TVector3 mom(qnan, qnan, qnan);
	    target_p_id_container.push_back(id);
	    target_p_dist2tgt_container.push_back(tp -> GetClosestDist());

	    Int_t repid = 0;
	    Int_t flag = 1;
	    for(Int_t i=0;i<2;i++){
	      Int_t temp = flag&event.pid[id];
	      if(temp==flag) repid += 1;
	      flag*=2;
	    }
	    if(!GFTrackCont.TrackCheck(id, repid)){
	      target_p_mass2_container.push_back(mass2);
	      target_p_mom_container.push_back(mom);
	      continue;
	    }

	    Int_t htofhitid_p; Double_t tracklen_p;
	    Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	    Bool_t htofextrapolation_p = GFTrackCont.TPCHTOFTrackMatching(id, repid, tgtpos,
									  event.HtofSeg, event.posHtof,
									  htofhitid_p, tof, tracklen_p,
									  pos, track2tgt_dist);
	    mom = GFTrackCont.GetMom(id, 0, repid);
	    if(htofextrapolation_p) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_p, event.tHtof[htofhitid_p]);
	    target_p_mass2_container.push_back(mass2);
	    target_p_mom_container.push_back(mom);
	  }
	  else{
	    Double_t mass2 = qnan;
	    TVector3 mom(qnan, qnan, qnan);
	    target_ppip_id_container.push_back(id);
	    target_ppip_dist2tgt_container.push_back(tp -> GetClosestDist());

	    Int_t repid = -1;
	    if(!GFTrackCont.TrackCheck(id, repid)){
	      target_ppip_mass2_container.push_back(mass2);
	      target_ppip_mom_container.push_back(mom);
	      continue;
	    }

	    Int_t htofhitid_ppi; Double_t tracklen_ppi;
	    Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	    Bool_t htofextrapolation = GFTrackCont.TPCHTOFTrackMatching(id, repid, tgtpos,
									event.HtofSeg, event.posHtof,
									htofhitid_ppi, tof, tracklen_ppi,
									pos, track2tgt_dist);
	    mom = GFTrackCont.GetMom(id, 0, repid);
	    if(htofextrapolation) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_ppi, event.tHtof[htofhitid_ppi]);
	    target_ppip_mass2_container.push_back(mass2);
	    target_ppip_mom_container.push_back(mom);
	  }
	}
      }

      for(Int_t icombi=0; icombi<more_pim_residuals.size(); ++icombi){
	Int_t id = more_pim_residuals[icombi];
	if(std::find(target_pim_id_container.begin(), target_pim_id_container.end(), id) ==
	   target_pim_id_container.end()){
	  TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( id );

	  Double_t mass2 = qnan;
	  TVector3 mom(qnan, qnan, qnan);
	  target_pim_id_container.push_back(id);
	  target_pim_dist2tgt_container.push_back(tp -> GetClosestDist());

	  Int_t repid = 0;
	  if(!GFTrackCont.TrackCheck(id, repid)){
	    target_pim_mass2_container.push_back(mass2);
	    target_pim_mom_container.push_back(mom);
	    continue;
	  }

	  Int_t htofhitid_pi; Double_t tracklen_pi;
	  Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	  Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(id, repid, tgtpos,
									 event.HtofSeg, event.posHtof,
									 htofhitid_pi, tof, tracklen_pi,
									 pos, track2tgt_dist);
	  mom = GFTrackCont.GetMom(id, 0, repid);
	  if(htofextrapolation_pi) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_pi, event.tHtof[htofhitid_pi]);
	  //if(htofextrapolation_pi && mass2 > 0.25) continue;
	  target_pim_mass2_container.push_back(mass2);
	  target_pim_mom_container.push_back(mom);
	}
      }

      for(Int_t icombi=0; icombi<more_pip_residuals.size(); ++icombi){
	Int_t id = more_pip_residuals[icombi];
	if(std::find(target_pip_id_container.begin(), target_pip_id_container.end(), id) ==
	   target_pip_id_container.end()){
	  TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( id );

	  Double_t mass2 = qnan;
	  TVector3 mom(qnan, qnan, qnan);
	  target_pip_id_container.push_back(id);
	  target_pip_dist2tgt_container.push_back(tp -> GetClosestDist());

	  Int_t repid = 0;
	  if(!GFTrackCont.TrackCheck(id, repid)){
	    target_pip_mass2_container.push_back(mass2);
	    target_pip_mom_container.push_back(mom);
	    continue;
	  }

	  Int_t htofhitid_pi; Double_t tracklen_pi;
	  Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	  Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(id, repid, tgtpos,
									 event.HtofSeg, event.posHtof,
									 htofhitid_pi, tof, tracklen_pi,
									 pos, track2tgt_dist);
	  mom = GFTrackCont.GetMom(id, 0, repid);
	  if(htofextrapolation_pi) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_pi, event.tHtof[htofhitid_pi]);
	  target_pip_mass2_container.push_back(mass2);
	  target_pip_mom_container.push_back(mom);
	}
      }

      event.pip_multi = target_pip_id_container.size();
      event.p_multi = target_p_id_container.size();
      event.pim_multi = target_pim_id_container.size();
      event.ppip_multi = target_ppip_id_container.size();
      event.accident_multi = target_accidental_id_container.size();
      for(Int_t itp=0; itp<target_p_id_container.size(); ++itp){
	Int_t id_p = target_p_id_container[itp];
	if(id_p==event.decays_id[0] || id_p==event.decays_id[2]){
	  event.p_multi -= 1;
	  continue;
	}
	event.residual_id.push_back(id_p);
	event.residual_dist2tgt.push_back(target_p_dist2tgt_container[itp]);
	event.residual_mass2.push_back(target_p_mass2_container[itp]);
	event.residual_mom.push_back(target_p_mom_container[itp].Mag());
	event.residual_mom_x.push_back(target_p_mom_container[itp].x());
	event.residual_mom_y.push_back(target_p_mom_container[itp].y());
	event.residual_mom_z.push_back(target_p_mom_container[itp].z());
	event.residual_charge.push_back(1);
	HF1( 1110, target_p_mom_container[itp].Mag());
      }
      for(Int_t itpip=0; itpip<target_pip_id_container.size(); ++itpip){
	Int_t id_pip = target_pip_id_container[itpip];
	if(id_pip==event.decays_id[0] || id_pip==event.decays_id[2]){
	  event.pip_multi -= 1;
	  continue;
	}
	event.residual_id.push_back(id_pip);
	event.residual_dist2tgt.push_back(target_pip_dist2tgt_container[itpip]);
	event.residual_mass2.push_back(target_pip_mass2_container[itpip]);
	event.residual_mom.push_back(target_pip_mom_container[itpip].Mag());
	event.residual_mom_x.push_back(target_pip_mom_container[itpip].x());
	event.residual_mom_y.push_back(target_pip_mom_container[itpip].y());
	event.residual_mom_z.push_back(target_pip_mom_container[itpip].z());
	event.residual_charge.push_back(1);
	HF1( 1111, target_pip_mom_container[itpip].Mag());
      }
      for(Int_t itpim=0; itpim<target_pim_id_container.size(); ++itpim){
	Int_t id_pim = target_pim_id_container[itpim];
	if(id_pim==event.decays_id[1] || id_pim==event.decays_id[3]){
	  event.pim_multi -= 1;
	  continue;
	}
	event.residual_id.push_back(id_pim);
	event.residual_dist2tgt.push_back(target_pim_dist2tgt_container[itpim]);
	event.residual_mass2.push_back(target_pim_mass2_container[itpim]);
	event.residual_mom.push_back(target_pim_mom_container[itpim].Mag());
	event.residual_mom_x.push_back(target_pim_mom_container[itpim].x());
	event.residual_mom_y.push_back(target_pim_mom_container[itpim].y());
	event.residual_mom_z.push_back(target_pim_mom_container[itpim].z());
	event.residual_charge.push_back(-1);
	HF1( 1112, target_pim_mom_container[itpim].Mag());
      }
      for(Int_t itppip=0; itppip<target_ppip_id_container.size(); ++itppip){
	Int_t id_ppip = target_ppip_id_container[itppip];
	if(id_ppip==event.decays_id[0] || id_ppip==event.decays_id[2]){
	  event.ppip_multi -= 1;
	  continue;
	}
	event.residual_id.push_back(id_ppip);
	event.residual_dist2tgt.push_back(target_ppip_dist2tgt_container[itppip]);
	event.residual_mass2.push_back(target_ppip_mass2_container[itppip]);
	event.residual_mom.push_back(target_ppip_mom_container[itppip].Mag());
	event.residual_mom_x.push_back(target_ppip_mom_container[itppip].x());
	event.residual_mom_y.push_back(target_ppip_mom_container[itppip].y());
	event.residual_mom_z.push_back(target_ppip_mom_container[itppip].z());
	event.residual_charge.push_back(event.charge[id_ppip]);
      }
      event.residual_multi = event.pip_multi + event.p_multi + event.pim_multi + event.ppip_multi + event.accident_multi;

      HF1( 1010, event.pim_multi);
      HF1( 1011, event.p_multi);
      HF1( 1012, event.pip_multi);

      HF1( 110, -event.BE[0]);
      HF1( 112, -event.BETPC[0]);
      return true;
    } //LL flag
  } //l_candidates

  std::vector<Int_t> GFxi_p_id_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_pi_id_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_pi2_id_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_p_rep_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_pi_rep_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_pi2_rep_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_id_container(xi_candidates, -1);
  std::vector<TVector3> GFxi_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFxi_decayvertex_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFxi_mass_container(xi_candidates, qnan);
  std::vector<TVector3> GFxi_pos_targetcenter_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFxi_mom_targetcenter_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFxi_dist_targetcenter_container(xi_candidates, qnan);
  std::vector<TVector3> GFl_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFl_vert_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFl_mass_container(xi_candidates, qnan);
  std::vector<Double_t> GFl_tracklen_container(xi_candidates, qnan);
  std::vector<Double_t> GFl_tof_container(xi_candidates, qnan);
  std::vector<TVector3> GFp_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFpi_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFpi2_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFp_mass2_container(xi_candidates, qnan);
  std::vector<Double_t> GFpi_mass2_container(xi_candidates, qnan);
  std::vector<Double_t> GFpi2_mass2_container(xi_candidates, qnan);
  std::vector<Double_t> GFp_tof_container(xi_candidates, qnan);
  std::vector<Double_t> GFpi_tof_container(xi_candidates, qnan);
  std::vector<Double_t> GFpi2_tof_container(xi_candidates, qnan);
  std::vector<Double_t> GFp_tracklen_container(xi_candidates, qnan);
  std::vector<Double_t> GFpi_tracklen_container(xi_candidates, qnan);
  std::vector<Double_t> GFpi2_tracklen_container(xi_candidates, qnan);
  std::vector<Int_t> GFp_htofid_container(xi_candidates, qnan);
  std::vector<Int_t> GFpi_htofid_container(xi_candidates, qnan);
  std::vector<Int_t> GFpi2_htofid_container(xi_candidates, qnan);
  std::vector<Double_t> GFppi_closedist_container(xi_candidates, qnan);
  std::vector<Double_t> GFlpi_closedist_container(xi_candidates, qnan);

  std::vector<Int_t> KFxi_id_container(xi_candidates, -1);
  std::vector<TVector3> KFl_mom_container0(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> KFxi_decayvertex_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> KFxi_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> KFl_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> KFxi_pos_targetcenter_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> KFxi_mom_targetcenter_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> KFxi_dist_targetcenter_container(xi_candidates, qnan);

  std::vector<Double_t> KFlchisqr_container(xi_candidates, qnan);
  std::vector<Double_t> KFlpval_container(xi_candidates, qnan);
  std::vector<Double_t> KFxi_chisqr_container(xi_candidates, qnan);
  std::vector<Double_t> KFxi_pval_container(xi_candidates, qnan);
  std::vector<Double_t> KFxi_mass_container(xi_candidates, qnan);
  std::vector<Double_t> KFlpi_closedist_container(xi_candidates, qnan);
  std::vector<std::vector<Double_t>> KFlpull_container(xi_candidates, std::vector<Double_t>(6, qnan));
  std::vector<std::vector<Double_t>> KFxipull_container(xi_candidates, std::vector<Double_t>(6, qnan));
  std::vector<TMatrixD> VXiContainer(xi_candidates, TMatrixD(3, 3));
  std::vector<TVector3> KFp_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> KFpi_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> KFpi2_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));

#if DebugDisp
  if(xi_candidates>0) std::cout<<"4. Xi candidates searching starts"<<std::endl;
#endif
  Int_t best_xi = -1; Double_t prev_Lmassdiff = 9999.;
  for(Int_t candi=0;candi<xi_candidates;candi++){
    Int_t trackid_p = xi_p_container[candi];
    Int_t trackid_pi = xi_pi_container[candi];
    Int_t trackid_pi2 = xi_pi2_container[candi];
    Int_t repid_p = p_repid_container[candi];
    Int_t repid_pi = pi_repid_container[candi];
    Int_t repid_pi2 = pi2_repid_container[candi];
    Int_t l_id = xi_l_container[candi];
    Int_t l_in_target = xi_l_intarget_container[candi];
    if(TMath::IsNaN(GFL_ppidist_container[l_id])) continue; //Genfit's fitting was succeeded.

    double GFppi_dist = GFL_ppidist_container[l_id];
    TVector3 GFlambda_vert = GFL_vtx_container[l_id];
    Double_t GFextrapolation_decays[3];
    GFextrapolation_decays[0] = GFL_p_extrapolation_container[l_id];
    GFextrapolation_decays[1] = GFL_pi_extrapolation_container[l_id];
    TVector3 GFmom_decays[3];
    GFmom_decays[0] = GFL_p_mom_container[l_id];
    GFmom_decays[1] = GFL_pi_mom_container[l_id];
    TLorentzVector GFLp(GFmom_decays[0], TMath::Hypot(GFmom_decays[0].Mag(), ProtonMass));
    TLorentzVector GFLpi(GFmom_decays[1], TMath::Hypot(GFmom_decays[1].Mag(), PionMass));
    TLorentzVector GFLlambda = GFLp + GFLpi;
    TVector3 GFlambda_mom = GFmom_decays[0] + GFmom_decays[1];

    TLorentzVector GFLlambda_fixed(GFlambda_mom, TMath::Sqrt(GFlambda_mom.Mag()*GFlambda_mom.Mag() + LambdaMass*LambdaMass));

    TVector3 GFxi_vert; Double_t GFlpi_dist = qnan; Double_t GFlambda_tracklen;
    if(!GFTrackCont.FindVertexXi(trackid_pi2, repid_pi2, GFlambda_vert, GFlambda_mom, GFlambda_tracklen,
				 GFextrapolation_decays[2], GFmom_decays[2], GFlpi_dist, GFxi_vert, vtx_scan_range)
       || GFlpi_dist > GFlpi_distcut) continue;

    TLorentzVector GFLpi2(GFmom_decays[2], TMath::Sqrt(GFmom_decays[2].Mag()*GFmom_decays[2].Mag() + PionMass*PionMass));
    TLorentzVector GFLxi = GFLlambda_fixed + GFLpi2;
    TVector3 GFxi_mom = GFlambda_mom + GFmom_decays[2];

    GFTrackCont.AddReconstructedTrack(XiMinusPdgCode, GFxi_vert, GFxi_mom);
    Int_t xi_trackid = GFTrackCont.GetNTrack() - 1;
    GFTrackCont.FitTrack(xi_trackid);

    TVector3 GFxipos_tgtcenter; TVector3 GFximom_tgtcenter;
    Double_t GFxitracklen_tgtcenter; Double_t GFxitof_tgtcenter;
    Bool_t xi_extrapolation = GFTrackCont.ExtrapolateToTargetCenter(xi_trackid,
								    GFxipos_tgtcenter,
								    GFximom_tgtcenter,
								    GFxitracklen_tgtcenter,
								    GFxitof_tgtcenter);
    if(!xi_extrapolation) continue;
    if(TMath::Abs(GFxipos_tgtcenter.y()) > GFxitarget_ycut) continue;
    event.xiflag = true; //Xi event

    TVector3 dist = tgtpos - GFxipos_tgtcenter;
    GFxi_dist_targetcenter_container[candi] = dist.Mag();
    GFxi_pos_targetcenter_container[candi] = GFxipos_tgtcenter;
    GFxi_mom_targetcenter_container[candi] = GFximom_tgtcenter;
    GFxi_id_container[candi] = xi_trackid;
#if KinematicFit
#if DebugDisp
    std::cout<<"KFLd"<<std::endl;
#endif
    Double_t KFchisqrxi = -1;
    Double_t KFpvalxi = -1;
    Double_t KFchisqrl = -1;
    Double_t KFpvall = -1;
    TPCLocalTrackHelix *track_p = TPCAna.GetTrackTPCHelix(trackid_p);
    auto Vp = track_p->GetCovarianceMatrix();
    TPCLocalTrackHelix *track_pi = TPCAna.GetTrackTPCHelix(trackid_pi);
    auto Vpi1 = track_pi->GetCovarianceMatrix();
    double Diag_ppi1[6]={
      Vp(0,0),Vp(1,1),Vp(2,2),Vpi1(0,0),Vpi1(1,1),Vpi1(2,2)
    };
    auto Offdiag_ppi1 = MathTools::MergeOffdiagonals(Vp,Vpi1);
    //Y and Z coordinates should be swapped in KF.
    TVector3 HTVP(GFLp.X(),GFLp.Z(),GFLp.Y());
    TVector3 HTVPi1(GFLpi.X(),GFLpi.Z(),GFLpi.Y());
    TVector3 HTVLd = HTVP+HTVPi1;
    TLorentzVector HLVP(HTVP,hypot(ProtonMass,HTVP.Mag()));
    TLorentzVector HLVPi1(HTVPi1,hypot(PionMass,HTVPi1.Mag()));
    TLorentzVector HLVLd(HTVLd,hypot(LambdaMass,HTVLd.Mag()));

    FourVectorFitter KFLd(HLVP,HLVPi1,HLVLd);
    KFLd.SetInvMass(LambdaMass);
    KFLd.SetMaximumStep(5);
    KFLd.SetVariance(Diag_ppi1);
    KFLd.AddOffdiagonals(Offdiag_ppi1);
    KFchisqrl = KFLd.DoKinematicFit();
    KFpvall = KFLd.GetPValue();
    auto HcontLd = KFLd.GetFittedLV();
    auto PullLd = KFLd.GetPull();
    auto KFHLVP = HcontLd.at(0);
    auto KFHLVPi1 = HcontLd.at(1);
    auto KFHLVLd = HcontLd.at(2);
    auto KFlambda_mom = TVector3(KFHLVLd.X(),KFHLVLd.Z(),KFHLVLd.Y());
    TLorentzVector KFLlambda_fixed(KFlambda_mom,hypot(KFlambda_mom.Mag(),LambdaMass));
    auto VLd = KFLd.GetUnmeasuredCovariance();

    TVector3 KFxi_vert; Double_t KFlpi_dist = qnan; Double_t KFlambda_tracklen;
    Double_t l_res_x, l_res_y, l_phi;
    MathTools::DecomposeResolution(VLd, KFlambda_mom, l_res_x, l_res_y, l_phi);
    GFTrackCont.FindVertexXi(trackid_pi2, repid_pi2, GFlambda_vert, KFlambda_mom,
			     KFlambda_tracklen, GFextrapolation_decays[2], GFmom_decays[2],
			     KFlpi_dist, KFxi_vert, vtx_scan_range, l_res_x, l_res_y, l_phi);
    /*
    if(!GFTrackCont.FindVertexXi(trackid_pi2, repid_pi2, GFlambda_vert, KFlambda_mom,
				 KFlambda_tracklen, GFextrapolation_decays[2], GFmom_decays[2],
				 KFlpi_dist, KFxi_vert, vtx_scan_range, l_res_x, l_res_y, l_phi)
      || KFlpi_dist > GFlpi_distcut) continue;
    */

#if DebugDisp
    std::cout<<Form("Resolution: %f %f %f",l_res_x,l_res_y,l_phi)<<std::endl;
#endif

    Double_t KFlambda_tof = Kinematics::CalcTimeOfFlight(KFlambda_mom.Mag(), KFlambda_tracklen, pdg::LambdaMass());
    TLorentzVector KFLpi2(GFmom_decays[2], TMath::Sqrt(GFmom_decays[2].Mag()*GFmom_decays[2].Mag() + PionMass*PionMass));
    TLorentzVector KFLxi = KFLlambda_fixed + KFLpi2;
    //TVector3 KFxi_mom = KFlambda_mom + GFmom_decays[2];

    auto* track_pi2 = TPCAna.GetTrackTPCHelix(trackid_pi2);
    auto VPi2 = track_pi2->GetCovarianceMatrix();
    double Diag_lpi2[6] =
      {VLd(0,0),VLd(1,1),VLd(2,2),VPi2(0,0),VPi2(1,1),VPi2(2,2)};
    auto Offdiag_lpi2 = MathTools::MergeOffdiagonals(VLd,VPi2);

    TVector3 HTVPi2(GFLpi2.X(),GFLpi2.Z(),GFLpi2.Y());
    TLorentzVector HLVPi2(HTVPi2,hypot(PionMass,HTVPi2.Mag()));
    auto HLVXi = KFHLVLd + HLVPi2;
    FourVectorFitter KFXi(KFHLVLd,HLVPi2,HLVXi);
    KFXi.SetInvMass(XiMinusMass);

#if DebugDisp
    std::cout<<"KFXi"<<std::endl;
#endif
    KFXi.SetMaximumStep(5);
    KFXi.SetVariance(Diag_lpi2);
    KFXi.AddOffdiagonals(Offdiag_lpi2);
    KFchisqrxi = KFXi.DoKinematicFit();
    KFpvalxi = KFXi.GetPValue();
    auto HcontXi = KFXi.GetFittedLV();
    auto PullXi = KFXi.GetPull();
    auto KFKFHLVLd = HcontXi.at(0);
    auto KFHLVPi2 = HcontXi.at(1);
    auto KFHLVXi = HcontXi.at(2);

    auto KFLVP = TLorentzVector(KFHLVP.X(),KFHLVP.Z(),KFHLVP.Y(),KFHLVP.E());
    auto KFLVPi1 = TLorentzVector(KFHLVPi1.X(),KFHLVPi1.Z(),KFHLVPi1.Y(),KFHLVPi1.E());
    auto KFLVPi2 = TLorentzVector(KFHLVPi2.X(),KFHLVPi2.Z(),KFHLVPi2.Y(),KFHLVPi2.E());
    auto KFLVLd = TLorentzVector(KFKFHLVLd.X(),KFKFHLVLd.Z(),KFKFHLVLd.Y(),KFKFHLVLd.E());
    auto KFLVXi = TLorentzVector(KFHLVXi.X(),KFHLVXi.Z(),KFHLVXi.Y(),KFHLVXi.E());

    auto KFTVP = KFLVP.Vect();
    auto KFTVPi1 = KFLVPi1.Vect();
    auto KFTVPi2 = KFLVPi2.Vect();
    auto KFTVLd = KFLVLd.Vect();
    auto KFTVXi = KFLVXi.Vect();
    auto VXi = KFXi.GetUnmeasuredCovariance();

    TVector3 KFxi_mom = KFTVXi;
    TVector3 KFlambda_mom_KFXi = KFTVLd;
    GFTrackCont.AddReconstructedTrack(XiMinusPdgCode, KFxi_vert, KFxi_mom);
    Int_t KFxi_trackid = GFTrackCont.GetNTrack() - 1;
    GFTrackCont.FitTrack(KFxi_trackid);
    TVector3 KFxipos_tgtcenter; TVector3 KFximom_tgtcenter;
    Double_t KFxitracklen_tgtcenter; Double_t KFxitof_tgtcenter;
    Bool_t KFxi_extrapolation = GFTrackCont.ExtrapolateToTargetCenter(KFxi_trackid,
								      KFxipos_tgtcenter,
								      KFximom_tgtcenter,
								      KFxitracklen_tgtcenter,
								      KFxitof_tgtcenter);
    KFxi_id_container[candi] = KFxi_trackid;
    KFxi_decayvertex_container[candi] = KFxi_vert;

    TVector3 KFdist = tgtpos - KFxipos_tgtcenter;
    KFxi_dist_targetcenter_container[candi] = KFdist.Mag();
    KFxi_pos_targetcenter_container[candi] = KFxipos_tgtcenter;
    KFxi_mom_targetcenter_container[candi] = KFximom_tgtcenter;
    KFxi_id_container[candi] = KFxi_trackid;

    KFxipull_container[candi] = PullXi;
    KFlpull_container[candi] = PullLd;
    KFl_mom_container0[candi] = KFlambda_mom;
    KFl_mom_container[candi] = KFlambda_mom_KFXi;
    KFlchisqr_container[candi] = KFchisqrl;
    KFlpval_container[candi] = KFpvall;
    KFxi_mom_container[candi] = KFxi_mom;
    KFxi_chisqr_container[candi] = KFchisqrxi;
    KFxi_pval_container[candi] = KFpvalxi;
    KFxi_mass_container[candi] = KFLxi.M();
    KFlpi_closedist_container[candi] = KFlpi_dist;
    KFp_mom_container[candi] = KFTVP;
    KFpi_mom_container[candi] = KFTVPi1;
    KFpi2_mom_container[candi] = KFTVPi2;
#endif
    Double_t GFlambda_tof = Kinematics::CalcTimeOfFlight(GFlambda_mom.Mag(), GFlambda_tracklen, pdg::LambdaMass());

    Int_t hitid_htof; Double_t tof_htof; Double_t tracklen_htof; TVector3 pos_htof; Double_t track2tgt_dist;
    Bool_t htofextrapolation_p = GFTrackCont.TPCHTOFTrackMatching(trackid_p, repid_p, GFlambda_vert,
								  event.HtofSeg, event.posHtof,
								  hitid_htof, tof_htof,
								  tracklen_htof, pos_htof, track2tgt_dist);
    if(htofextrapolation_p){
      GFp_htofid_container[candi] = hitid_htof;
      //GFp_tof_container[candi] = event.tHtof[hitid_htof] - GFlambda_tof - GFxi_tof;
      GFp_tof_container[candi] = event.tHtof[hitid_htof] - GFlambda_tof;
      GFp_tracklen_container[candi] = tracklen_htof - GFextrapolation_decays[0];
      GFp_mass2_container[candi] = Kinematics::MassSquare(GFmom_decays[0].Mag(),
							  GFp_tracklen_container[candi],
							  GFp_tof_container[candi]);
      //if(GFL_p_mass2_container[candi] < 0.25) continue;
    }

    Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(trackid_pi, repid_pi, GFlambda_vert,
								   event.HtofSeg, event.posHtof,
								   hitid_htof, tof_htof,
								   tracklen_htof, pos_htof, track2tgt_dist);
    if(htofextrapolation_pi){
      GFpi_htofid_container[candi] = hitid_htof;
      //GFpi_tof_container[candi] = event.tHtof[hitid_htof] - GFlambda_tof - GFxi_tof;
      GFpi_tof_container[candi] = event.tHtof[hitid_htof] - GFlambda_tof;
      GFpi_tracklen_container[candi] = tracklen_htof - GFextrapolation_decays[1];
      GFpi_mass2_container[candi] = Kinematics::MassSquare(GFmom_decays[1].Mag(),
							   GFpi_tracklen_container[candi],
							   GFpi_tof_container[candi]);
      //if(GFL_pi_mass2_container[candi] > 0.25) continue;
    }

    Bool_t htofextrapolation_pi2 = GFTrackCont.TPCHTOFTrackMatching(trackid_pi2, repid_pi2, tgtpos,
								    event.HtofSeg, event.posHtof,
								    hitid_htof, tof_htof,
								    tracklen_htof, pos_htof, track2tgt_dist);
    if(htofextrapolation_pi2){
      GFpi2_htofid_container[candi] = hitid_htof;
      //GFpi2_tof_container[candi] = event.tHtof[hitid_htof] - GFxi_tof;
      GFpi2_tof_container[candi] = event.tHtof[hitid_htof];
      GFpi2_tracklen_container[candi] = tracklen_htof - GFextrapolation_decays[2];
      GFpi2_mass2_container[candi] = Kinematics::MassSquare(GFmom_decays[2].Mag(),
							    GFpi2_tracklen_container[candi],
							    GFpi2_tof_container[candi]);
      //if(GFL_pi2_mass2_container[candi] > 0.25) continue;
    }

    GFxi_p_id_container[candi] = trackid_p;
    GFxi_pi_id_container[candi] = trackid_pi;
    GFxi_pi2_id_container[candi] = trackid_pi2;
    GFxi_p_rep_container[candi] = repid_p;
    GFxi_pi_rep_container[candi] = repid_pi;
    GFxi_pi2_rep_container[candi] = repid_pi2;
    GFxi_mom_container[candi] = GFxi_mom;
    GFxi_decayvertex_container[candi] = GFxi_vert;
    GFxi_mass_container[candi] = GFLxi.M();

    GFl_mom_container[candi] = GFlambda_mom;
    GFl_vert_container[candi] = GFlambda_vert;
    GFl_mass_container[candi] = GFLlambda.M();
    GFl_tracklen_container[candi] = GFlambda_tracklen;
    GFl_tof_container[candi] = GFlambda_tof;
    GFp_mom_container[candi] = GFmom_decays[0];
    GFpi_mom_container[candi] = GFmom_decays[1];
    GFpi2_mom_container[candi] = GFmom_decays[2];
    GFppi_closedist_container[candi] = GFppi_dist;
    GFlpi_closedist_container[candi] = GFlpi_dist;
    Double_t diff = TMath::Abs(GFLlambda.M() - LambdaMass);
    //Double_t diff = TMath::Abs(lambda_mass_container[candi] - LambdaMass);
    if(prev_Lmassdiff > diff){
      prev_Lmassdiff = diff;
      best_xi = candi;
    }
  } //candi

  if(event.xiflag){
#if DebugDisp
    std::cout<<"5. Calculate and save the best Xi combination"<<std::endl;
#endif
    event.ximass = xi_mass_container[best_xi];
    event.xidecayvtx_x = xi_decayvertex_container[best_xi].x();
    event.xidecayvtx_y = xi_decayvertex_container[best_xi].y();
    event.xidecayvtx_z = xi_decayvertex_container[best_xi].z();
    event.ximom = xi_mom_container[best_xi].Mag();
    event.ximom_x = xi_mom_container[best_xi].x();
    event.ximom_y = xi_mom_container[best_xi].y();
    event.ximom_z = xi_mom_container[best_xi].z();
    event.lpi_dist = lpi_closedist[best_xi];
    event.xitarget_dist = xi_targetdist_container[best_xi];
    event.xitargetvtx_x = xi_targetvtx_container[best_xi].x();
    event.xitargetvtx_y = xi_targetvtx_container[best_xi].y();
    event.xitargetvtx_z = xi_targetvtx_container[best_xi].z();
    event.xitargetmom = xi_targetmom_container[best_xi].Mag();
    event.xitargetmom_x = xi_targetmom_container[best_xi].x();
    event.xitargetmom_y = xi_targetmom_container[best_xi].y();
    event.xitargetmom_z = xi_targetmom_container[best_xi].z();
    event.lmass = lambda_mass_container[best_xi];
    event.ldecayvtx_x = l_vert_container[best_xi].x();
    event.ldecayvtx_y = l_vert_container[best_xi].y();
    event.ldecayvtx_z = l_vert_container[best_xi].z();
    event.l_intarget = xi_l_intarget_container[best_xi];
    event.lmom = l_mom_container[best_xi].Mag();
    event.lmom_x = l_mom_container[best_xi].x();
    event.lmom_y = l_mom_container[best_xi].y();
    event.lmom_z = l_mom_container[best_xi].z();
    event.ppi_dist = ppi_closedist[best_xi];

    event.decays_mom.push_back(xi_p_mom_container[best_xi].Mag());
    event.decays_mom.push_back(xi_pi_mom_container[best_xi].Mag());
    event.decays_mom.push_back(xi_pi2_mom_container[best_xi].Mag());
    event.decays_mom_x.push_back(xi_p_mom_container[best_xi].x());
    event.decays_mom_x.push_back(xi_pi_mom_container[best_xi].x());
    event.decays_mom_x.push_back(xi_pi2_mom_container[best_xi].x());
    event.decays_mom_y.push_back(xi_p_mom_container[best_xi].y());
    event.decays_mom_y.push_back(xi_pi_mom_container[best_xi].y());
    event.decays_mom_y.push_back(xi_pi2_mom_container[best_xi].y());
    event.decays_mom_z.push_back(xi_p_mom_container[best_xi].z());
    event.decays_mom_z.push_back(xi_pi_mom_container[best_xi].z());
    event.decays_mom_z.push_back(xi_pi2_mom_container[best_xi].z());
    event.decays_id.push_back(xi_p_container[best_xi]);
    event.decays_id.push_back(xi_pi_container[best_xi]);
    event.decays_id.push_back(xi_pi2_container[best_xi]);
    event.decays_G4tid.push_back(G4TrackID.at(xi_p_container[best_xi]));
    event.decays_G4tid.push_back(G4TrackID.at(xi_pi_container[best_xi]));
    event.decays_G4tid.push_back(G4TrackID.at(xi_pi2_container[best_xi]));
    event.decays_purity.push_back(event.purity[xi_p_container[best_xi]]);
    event.decays_purity.push_back(event.purity[xi_pi_container[best_xi]]);
    event.decays_purity.push_back(event.purity[xi_pi2_container[best_xi]]);
    event.decays_efficiency.push_back(event.efficiency[xi_p_container[best_xi]]);
    event.decays_efficiency.push_back(event.efficiency[xi_pi_container[best_xi]]);
    event.decays_efficiency.push_back(event.efficiency[xi_pi2_container[best_xi]]);

    TLorentzVector Lv_p(xi_p_mom_container[best_xi],
			TMath::Hypot(xi_p_mom_container[best_xi].Mag(), ProtonMass));
    TLorentzVector Lv_pi(xi_pi_mom_container[best_xi],
			 TMath::Hypot(xi_pi_mom_container[best_xi].Mag(), PionMass));
    TLorentzVector Lv_pi2(xi_pi2_mom_container[best_xi],
			  TMath::Hypot(xi_pi2_mom_container[best_xi].Mag(), PionMass));
    Lv_p.Boost(-xi_mom_container[best_xi]);
    Lv_pi.Boost(-xi_mom_container[best_xi]);
    Lv_pi2.Boost(-xi_mom_container[best_xi]);
    event.decays_CMmom.push_back(Lv_p.P());
    event.decays_CMmom_x.push_back(Lv_p.Px());
    event.decays_CMmom_y.push_back(Lv_p.Py());
    event.decays_CMmom_z.push_back(Lv_p.Pz());
    event.decays_CMmom.push_back(Lv_pi.P());
    event.decays_CMmom_x.push_back(Lv_pi.Px());
    event.decays_CMmom_y.push_back(Lv_pi.Py());
    event.decays_CMmom_z.push_back(Lv_pi.Pz());
    event.decays_CMmom.push_back(Lv_pi2.P());
    event.decays_CMmom_x.push_back(Lv_pi2.Px());
    event.decays_CMmom_y.push_back(Lv_pi2.Py());
    event.decays_CMmom_z.push_back(Lv_pi2.Pz());

    event.GFximass = GFxi_mass_container[best_xi];
    event.GFxidecayvtx_x = GFxi_decayvertex_container[best_xi].x();
    event.GFxidecayvtx_y = GFxi_decayvertex_container[best_xi].y();
    event.GFxidecayvtx_z = GFxi_decayvertex_container[best_xi].z();
    event.GFximom = GFxi_mom_container[best_xi].Mag();
    event.GFximom_x = GFxi_mom_container[best_xi].x();
    event.GFximom_y = GFxi_mom_container[best_xi].y();
    event.GFximom_z = GFxi_mom_container[best_xi].z();
    event.GFlpi_dist = GFlpi_closedist_container[best_xi];
    event.GFlmass = GFl_mass_container[best_xi];
    event.GFldecayvtx_x = GFl_vert_container[best_xi].x();
    event.GFldecayvtx_y = GFl_vert_container[best_xi].y();
    event.GFldecayvtx_z = GFl_vert_container[best_xi].z();
    event.GFlmom = GFl_mom_container[best_xi].Mag();
    event.GFlmom_x = GFl_mom_container[best_xi].x();
    event.GFlmom_y = GFl_mom_container[best_xi].y();
    event.GFlmom_z = GFl_mom_container[best_xi].z();
    event.GFltracklen = GFl_tracklen_container[best_xi];
    event.GFltof = GFl_tof_container[best_xi];
    event.GFppi_dist = GFppi_closedist_container[best_xi];

    event.GFntdecays = 3;
    event.GFnhtrack.resize(event.GFntdecays);
    event.GFchisqr.resize(event.GFntdecays);
    event.GFcharge.resize(event.GFntdecays);
    event.GFtof.resize(event.GFntdecays);
    event.GFtracklen.resize(event.GFntdecays);
    event.GFpval.resize(event.GFntdecays);
    event.GFpdgcode.resize(event.GFntdecays);

    event.GFdecays_htofid.resize(event.GFntdecays);
    event.GFdecays_tracklen.resize(event.GFntdecays);
    event.GFdecays_tof.resize(event.GFntdecays);
    event.GFdecays_mass2.resize(event.GFntdecays);
    event.GFdecays_mom.resize(event.GFntdecays);
    event.GFdecays_mom_x.resize(event.GFntdecays);
    event.GFdecays_mom_y.resize(event.GFntdecays);
    event.GFdecays_mom_z.resize(event.GFntdecays);
    event.GFdecays_CMmom.resize(event.GFntdecays);
    event.GFdecays_CMmom_x.resize(event.GFntdecays);
    event.GFdecays_CMmom_y.resize(event.GFntdecays);
    event.GFdecays_CMmom_z.resize(event.GFntdecays);
    event.GFmomloss.resize(event.GFntdecays);
    event.GFeloss.resize(event.GFntdecays);

    for(Int_t j=0; j<event.GFntdecays; ++j){
      Int_t igf = GFxi_p_id_container[best_xi];
      if(j==1) igf = GFxi_pi_id_container[best_xi];
      if(j==2) igf = GFxi_pi2_id_container[best_xi];

      Int_t repid = GFxi_p_rep_container[best_xi];
      if(j==1) repid = GFxi_pi_rep_container[best_xi];
      if(j==2) repid = GFxi_pi2_rep_container[best_xi];

      Int_t GFhtofid_decays = GFp_htofid_container[best_xi];
      if(j==1) GFhtofid_decays = GFpi_htofid_container[best_xi];
      if(j==2) GFhtofid_decays = GFpi2_htofid_container[best_xi];

      Double_t GFtof_decays = GFp_tof_container[best_xi];
      if(j==1) GFtof_decays = GFpi_tof_container[best_xi];
      if(j==2) GFtof_decays = GFpi2_tof_container[best_xi];

      Double_t GFtracklen_decays = GFp_tracklen_container[best_xi];
      if(j==1) GFtracklen_decays = GFpi_tracklen_container[best_xi];
      if(j==2) GFtracklen_decays = GFpi2_tracklen_container[best_xi];

      TVector3 GFmom_decays = GFp_mom_container[best_xi];
      if(j==1) GFmom_decays = GFpi_mom_container[best_xi];
      if(j==2) GFmom_decays = GFpi2_mom_container[best_xi];

      Double_t GFmass2_decays = GFp_mass2_container[best_xi];
      if(j==1) GFmass2_decays = GFpi_mass2_container[best_xi];
      if(j==2) GFmass2_decays = GFpi2_mass2_container[best_xi];

      event.GFdecays_htofid[j] = GFhtofid_decays;
      event.GFdecays_tracklen[j] = GFtracklen_decays;
      event.GFdecays_tof[j] = GFtof_decays;
      event.GFdecays_mass2[j] = GFmass2_decays;
      event.GFdecays_mom[j] = GFmom_decays.Mag();
      event.GFdecays_mom_x[j] = GFmom_decays.x();
      event.GFdecays_mom_y[j] = GFmom_decays.y();
      event.GFdecays_mom_z[j] = GFmom_decays.z();
      event.GFmomloss[j] = GFmom_decays.Mag() - GFTrackCont.GetMom(igf, 0, repid).Mag();
      event.GFeloss[j] = TMath::Hypot(GFmom_decays.Mag(), pdgmass[j]) - TMath::Hypot(GFTrackCont.GetMom(igf, 0, repid).Mag(), pdgmass[j]);

      Double_t mass = PionMass;
      if(j==0) mass = ProtonMass;

      TLorentzVector GFLv_decays(GFmom_decays, TMath::Hypot(GFmom_decays.Mag(), mass));
      GFLv_decays.Boost(-GFxi_mom_container[best_xi]);
      event.GFdecays_CMmom[j] = GFLv_decays.P();
      event.GFdecays_CMmom_x[j] = GFLv_decays.Px();
      event.GFdecays_CMmom_y[j] = GFLv_decays.Py();
      event.GFdecays_CMmom_z[j] = GFLv_decays.Pz();

      event.GFnhtrack[j] = GFTrackCont.GetNHits(igf);
      event.GFchisqr[j] = GFTrackCont.GetChi2NDF(igf, repid);
      event.GFcharge[j] = GFTrackCont.GetCharge(igf, repid);
      event.GFtof[j] = GFTrackCont.GetTrackTOF(igf, 0, -1, repid);
      event.GFtracklen[j] = GFTrackCont.GetTrackLength(igf, 0, -1, repid);
      event.GFpval[j] = GFTrackCont.GetPvalue(igf, repid);
      event.GFpdgcode[j] = GFTrackCont.GetPDGcode(igf, repid);
    } //j : decays

    TVector3 xipos_tgtcenter = GFxi_pos_targetcenter_container[best_xi];
    TVector3 ximom_tgtcenter = GFxi_mom_targetcenter_container[best_xi];
    event.GFxitargetcenter_x = xipos_tgtcenter.x();
    event.GFxitargetcenter_y = xipos_tgtcenter.y();
    event.GFxitargetcenter_z = xipos_tgtcenter.z();
    event.GFxitargetcentermom = ximom_tgtcenter.Mag();
    event.GFxitargetcentermom_x = ximom_tgtcenter.x();
    event.GFxitargetcentermom_y = ximom_tgtcenter.y();
    event.GFxitargetcentermom_z = ximom_tgtcenter.z();
    event.GFxitargetcenter_dist = GFxi_dist_targetcenter_container[best_xi];

    const Int_t ntrack_xi = 3;
    Double_t x0[ntrack_xi] = {event.xtgtTPCKurama[0], event.xtgtK18[0],
			      xipos_tgtcenter.x()};
    Double_t y0[ntrack_xi] = {event.ytgtTPCKurama[0], event.ytgtK18[0],
			      xipos_tgtcenter.y()};
    Double_t u0[ntrack_xi] = {event.utgtTPCKurama[0], event.utgtK18[0],
			      ximom_tgtcenter.x()/ximom_tgtcenter.z()};
    Double_t v0[ntrack_xi] = {event.vtgtTPCKurama[0], event.vtgtK18[0],
			      ximom_tgtcenter.y()/ximom_tgtcenter.z()};

    TVector3 kk_xi_vertex = Kinematics::MultitrackVertex(ntrack_xi, x0, y0, u0, v0);
    event.GFprodvtx_x_kkxi = kk_xi_vertex.x();
    event.GFprodvtx_y_kkxi = kk_xi_vertex.y();
    event.GFprodvtx_z_kkxi = kk_xi_vertex.z();

    TVector3 xipos_prodvtx; TVector3 ximom_prodvtx;
    Double_t xitracklen_prodvtx; Double_t xitof_prodvtx;
    Bool_t xi_extrapolation_prodvtx = GFTrackCont.XiDecayToProdVertex(GFxi_id_container[best_xi],
								      kk_xi_vertex, xipos_prodvtx,
								      ximom_prodvtx,
								      xitracklen_prodvtx,
								      xitof_prodvtx);

    if(xi_extrapolation_prodvtx){
      TVector3 dist = xipos_prodvtx - kk_xi_vertex;
      event.GFxiprodvtx_dist = dist.Mag();
      event.GFxiprodvtx_x = xipos_prodvtx.x();
      event.GFxiprodvtx_y = xipos_prodvtx.y();
      event.GFxiprodvtx_z = xipos_prodvtx.z();
      event.GFxiprodmom = ximom_prodvtx.Mag();
      event.GFxiprodmom_x = ximom_prodvtx.x();
      event.GFxiprodmom_y = ximom_prodvtx.y();
      event.GFxiprodmom_z = ximom_prodvtx.z();
      event.GFxitracklen = xitracklen_prodvtx;
      event.GFxitof = xitof_prodvtx;
      event.GFximomloss = ximom_prodvtx.Mag() - event.GFximom;

      TLorentzVector LvXi_fixedmass(ximom_prodvtx, TMath::Hypot(ximom_prodvtx.Mag(), XiMinusMass));
      TLorentzVector LvRcXi_fixedmass = LvRcTPC - LvXi_fixedmass;
      event.GFxiexcitation = LvRcXi_fixedmass.M() - m11B - me;

      HF1( 200, LvRcXi_fixedmass.M() - m11B - me);

#if DebugDisp
      std::cout<<"K-K+Xi vertex "<<kk_xi_vertex
	       <<" Xi extrapolate to the vertex "<<xipos_prodvtx<<std::endl;
#endif
    }

    TVector3 xipos_kk; TVector3 ximom_kk;
    Double_t xitracklen_kk; Double_t xitof_kk;
    Bool_t xi_extrapolation = GFTrackCont.XiDecayToProdVertex(GFxi_id_container[best_xi],
							      kkvtxTPC,
							      xipos_kk,
							      ximom_kk,
							      xitracklen_kk,
							      xitof_kk);
    if(xi_extrapolation){
      TVector3 dist = xipos_kk - kkvtxTPC;
      event.GFxikkvtx_x = xipos_kk.x();
      event.GFxikkvtx_y = xipos_kk.y();
      event.GFxikkvtx_z = xipos_kk.z();
      event.GFxikkmom = ximom_kk.Mag();
      event.GFxikkmom_x = ximom_kk.x();
      event.GFxikkmom_y = ximom_kk.y();
      event.GFxikkmom_z = ximom_kk.z();
      event.GFxikkvtx_dist = dist.Mag();
    }

    TVector3 xipos_tgt; TVector3 ximom_tgt;
    Double_t xitracklen_tgt; Double_t xitof_tgt;
    xi_extrapolation = GFTrackCont.XiDecayToProdVertex(GFxi_id_container[best_xi],
						       tgtpos,
						       xipos_tgt,
						       ximom_tgt,
						       xitracklen_tgt,
						       xitof_tgt);
    if(xi_extrapolation){
      TVector3 dist = tgtpos - xipos_tgt;
      event.GFxitargetvtx_x = xipos_tgt.x();
      event.GFxitargetvtx_y = xipos_tgt.y();
      event.GFxitargetvtx_z = xipos_tgt.z();
      event.GFxitargetmom = ximom_tgt.Mag();
      event.GFxitargetmom_x = ximom_tgt.x();
      event.GFxitargetmom_y = ximom_tgt.y();
      event.GFxitargetmom_z = ximom_tgt.z();
      event.GFxitarget_dist = dist.Mag();
    }
#if KinematicFit
    event.KFlmom0 = KFl_mom_container0[best_xi].Mag();
    event.KFlmom_x0 = KFl_mom_container0[best_xi].x();
    event.KFlmom_y0 = KFl_mom_container0[best_xi].y();
    event.KFlmom_z0 = KFl_mom_container0[best_xi].z();
    event.KFlmom = KFl_mom_container[best_xi].Mag();
    event.KFlmom_x = KFl_mom_container[best_xi].x();
    event.KFlmom_y = KFl_mom_container[best_xi].y();
    event.KFlmom_z = KFl_mom_container[best_xi].z();
    event.KFlchisqr = KFlchisqr_container[best_xi];
    event.KFlpval = KFlpval_container[best_xi];
    event.KFlpull = KFlpull_container[best_xi];
    event.KFlpi_dist = KFlpi_closedist_container[best_xi];
    event.KFximom = KFxi_mom_container[best_xi].Mag();
    event.KFximom_x = KFxi_mom_container[best_xi].x();
    event.KFximom_y = KFxi_mom_container[best_xi].y();
    event.KFximom_z = KFxi_mom_container[best_xi].z();
    event.KFxichisqr = KFxi_chisqr_container[best_xi];
    event.KFxipval = KFxi_pval_container[best_xi];
    event.KFximass = KFxi_mass_container[best_xi];
    event.KFxidecayvtx_x = KFxi_decayvertex_container[best_xi].x();
    event.KFxidecayvtx_y = KFxi_decayvertex_container[best_xi].y();
    event.KFxidecayvtx_z = KFxi_decayvertex_container[best_xi].z();
    event.KFxipull = KFxipull_container[best_xi];
    event.KFdecays_mom.push_back(KFp_mom_container[best_xi].Mag());
    event.KFdecays_mom_x.push_back(KFp_mom_container[best_xi].x());
    event.KFdecays_mom_y.push_back(KFp_mom_container[best_xi].y());
    event.KFdecays_mom_z.push_back(KFp_mom_container[best_xi].z());
    event.KFdecays_mom.push_back(KFpi_mom_container[best_xi].Mag());
    event.KFdecays_mom_x.push_back(KFpi_mom_container[best_xi].x());
    event.KFdecays_mom_y.push_back(KFpi_mom_container[best_xi].y());
    event.KFdecays_mom_z.push_back(KFpi_mom_container[best_xi].z());
    event.KFdecays_mom.push_back(KFpi2_mom_container[best_xi].Mag());
    event.KFdecays_mom_x.push_back(KFpi2_mom_container[best_xi].x());
    event.KFdecays_mom_y.push_back(KFpi2_mom_container[best_xi].y());
    event.KFdecays_mom_z.push_back(KFpi2_mom_container[best_xi].z());

    TLorentzVector KFLv_p(KFp_mom_container[best_xi],
			TMath::Hypot(KFp_mom_container[best_xi].Mag(), ProtonMass));
    TLorentzVector KFLv_pi(KFpi_mom_container[best_xi],
			 TMath::Hypot(KFpi_mom_container[best_xi].Mag(), PionMass));
    TLorentzVector KFLv_pi2(KFpi2_mom_container[best_xi],
			  TMath::Hypot(KFpi2_mom_container[best_xi].Mag(), PionMass));
    KFLv_p.Boost(-KFxi_mom_container[best_xi]);
    KFLv_pi.Boost(-KFxi_mom_container[best_xi]);
    KFLv_pi2.Boost(-KFxi_mom_container[best_xi]);
    event.KFdecays_CMmom.push_back(KFLv_p.P());
    event.KFdecays_CMmom_x.push_back(KFLv_p.Px());
    event.KFdecays_CMmom_y.push_back(KFLv_p.Py());
    event.KFdecays_CMmom_z.push_back(KFLv_p.Pz());
    event.KFdecays_CMmom.push_back(KFLv_pi.P());
    event.KFdecays_CMmom_x.push_back(KFLv_pi.Px());
    event.KFdecays_CMmom_y.push_back(KFLv_pi.Py());
    event.KFdecays_CMmom_z.push_back(KFLv_pi.Pz());
    event.KFdecays_CMmom.push_back(KFLv_pi2.P());
    event.KFdecays_CMmom_x.push_back(KFLv_pi2.Px());
    event.KFdecays_CMmom_y.push_back(KFLv_pi2.Py());
    event.KFdecays_CMmom_z.push_back(KFLv_pi2.Pz());

    for(int i=0;i<event.KFlpull.size();++i){
      int hn = 10000 + 10 + i;
      auto pull = event.KFlpull.at(i);
      HF1(hn,pull);
    }
    for(int i=0;i<event.KFxipull.size();++i){
      int hn = 20000 + 10 + i;
      auto pull = event.KFxipull.at(i);
      HF1(hn,pull);
    }

    HF1(10000,event.KFlpval);
    HF1(10001,event.KFlchisqr);
    HF1(20002,event.KFximass);
    HF1(10003,event.KFlpi_dist);

    TVector3 KFxipos_kk; TVector3 KFximom_kk;
    Double_t KFxitracklen_kk; Double_t KFxitof_kk;
    Bool_t KFxi_extrapolation_kkvtx = GFTrackCont.XiDecayToProdVertex(KFxi_id_container[best_xi],
								      kkvtxTPC,
								      KFxipos_kk,
								      KFximom_kk,
								      KFxitracklen_kk,
								      KFxitof_kk);
    if(KFxi_extrapolation_kkvtx){
      TVector3 dist = KFxipos_kk - kkvtxTPC;
      event.KFxi_kkvtx_x = KFxipos_kk.x();
      event.KFxi_kkvtx_y = KFxipos_kk.y();
      event.KFxi_kkvtx_z = KFxipos_kk.z();
      event.KFxi_kkvtx_mom = KFximom_kk.Mag();
      event.KFxi_kkvtx_mom_x = KFximom_kk.x();
      event.KFxi_kkvtx_mom_y = KFximom_kk.y();
      event.KFxi_kkvtx_mom_z = KFximom_kk.z();
      event.KFxi_kkvtx_dist = dist.Mag();
    }

    TVector3 KFxipos_tgt; TVector3 KFximom_tgt;
    Double_t KFxitracklen_tgt; Double_t KFxitof_tgt;
    Bool_t KFxi_extrapolation_tgt = GFTrackCont.XiDecayToProdVertex(KFxi_id_container[best_xi],
								    tgtpos,
								    KFxipos_tgt,
								    KFximom_tgt,
								    KFxitracklen_tgt,
								    KFxitof_tgt);
    if(KFxi_extrapolation_tgt){
      TVector3 dist = tgtpos - KFxipos_tgt;
      event.KFxitargetvtx_x = KFxipos_tgt.x();
      event.KFxitargetvtx_y = KFxipos_tgt.y();
      event.KFxitargetvtx_z = KFxipos_tgt.z();
      event.KFxitargetmom = KFximom_tgt.Mag();
      event.KFxitargetmom_x = KFximom_tgt.x();
      event.KFxitargetmom_y = KFximom_tgt.y();
      event.KFxitargetmom_z = KFximom_tgt.z();
      event.KFxitarget_dist = dist.Mag();
    }

    TVector3 KFxipos_tgtcenter = KFxi_pos_targetcenter_container[best_xi];
    TVector3 KFximom_tgtcenter = KFxi_mom_targetcenter_container[best_xi];
    event.KFxitargetcenter_x = KFxipos_tgtcenter.x();
    event.KFxitargetcenter_y = KFxipos_tgtcenter.y();
    event.KFxitargetcenter_z = KFxipos_tgtcenter.z();
    event.KFxitargetcentermom = KFximom_tgtcenter.Mag();
    event.KFxitargetcentermom_x = KFximom_tgtcenter.x();
    event.KFxitargetcentermom_y = KFximom_tgtcenter.y();
    event.KFxitargetcentermom_z = KFximom_tgtcenter.z();
    event.KFxitargetcenter_dist = KFxi_dist_targetcenter_container[best_xi];
    //const Int_t ntrack_xi = 3;
    Double_t KFxtgt[ntrack_xi] = {event.xtgtTPCKurama[0], event.xtgtK18[0],
				  KFxipos_tgtcenter.x()};
    Double_t KFytgt[ntrack_xi] = {event.ytgtTPCKurama[0], event.ytgtK18[0],
				  KFxipos_tgtcenter.y()};
    Double_t KFutgt[ntrack_xi] = {event.utgtTPCKurama[0], event.utgtK18[0],
				  KFximom_tgtcenter.x()/KFximom_tgtcenter.z()};
    Double_t KFvtgt[ntrack_xi] = {event.vtgtTPCKurama[0], event.vtgtK18[0],
				  KFximom_tgtcenter.y()/KFximom_tgtcenter.z()};
    Double_t resxi_u, resxi_v;
    TVector3 KFXimom(event.KFximom_x, event.KFximom_y, event.KFximom_z);
    MathTools::DecomposeResolutionUV(VXiContainer[best_xi], KFXimom, resxi_u, resxi_v);
    std::vector<double> res_x = {res_xKurama, res_xK18, res_xXiVtx};
    std::vector<double> res_y = {res_yKurama, res_yK18, res_yXiVtx};
    std::vector<double> res_u = {res_uKurama, res_uK18, resxi_u};
    std::vector<double> res_v = {res_vKurama, res_vK18, resxi_v};

    Double_t chisqr_xikk;
    TVector3 KFkkxi_prodvertex = Kinematics::MultitrackVertex(ntrack_xi,
							      KFxtgt, KFytgt,
							      KFutgt, KFvtgt,
							      res_x, res_y,
							      res_u, res_v,
							      chisqr_xikk);
    event.KFprodvtx_chisqr_kkxi = chisqr_xikk;
    event.KFprodvtx_x_kkxi = KFkkxi_prodvertex.x();
    event.KFprodvtx_y_kkxi = KFkkxi_prodvertex.y();
    event.KFprodvtx_z_kkxi = KFkkxi_prodvertex.z();

    TVector3 KFxipos_prodvtx; TVector3 KFximom_prodvtx;
    Double_t KFxitracklen_prodvtx; Double_t KFxitof_prodvtx;
    Bool_t KFxi_extrapolation_prodvtx = GFTrackCont.XiDecayToProdVertex(KFxi_id_container[best_xi],
									KFkkxi_prodvertex,
									KFxipos_prodvtx,
									KFximom_prodvtx,
									KFxitracklen_prodvtx,
									KFxitof_prodvtx);
    if(KFxi_extrapolation_prodvtx){
      TVector3 dist = KFxipos_prodvtx - KFkkxi_prodvertex;
      event.KFxiprodvtx_dist = dist.Mag();
      event.KFxiprodvtx_x = KFxipos_prodvtx.x();
      event.KFxiprodvtx_y = KFxipos_prodvtx.y();
      event.KFxiprodvtx_z = KFxipos_prodvtx.z();
      event.KFxiprodmom = KFximom_prodvtx.Mag();
      event.KFxiprodmom_x = KFximom_prodvtx.x();
      event.KFxiprodmom_y = KFximom_prodvtx.y();
      event.KFxiprodmom_z = KFximom_prodvtx.z();
      event.KFxitracklen = KFxitracklen_prodvtx;
      event.KFxitof = KFxitof_prodvtx;
      event.KFximomloss = KFximom_prodvtx.Mag() - event.KFximom;
      TLorentzVector LvXi_fixedmass(KFximom_prodvtx, TMath::Hypot(KFximom_prodvtx.Mag(), XiMinusMass));
      TLorentzVector LvRcXi_fixedmass = LvRcTPC - LvXi_fixedmass;
      event.KFxiexcitation = LvRcXi_fixedmass.M() - m11B - me;
    }

    //Kp, Xi vertex
    Double_t KFxtgtkpxi[2] = {KFxtgt[0], KFxtgt[2]};
    Double_t KFytgtkpxi[2] = {KFytgt[0], KFytgt[2]};
    Double_t KFutgtkpxi[2] = {KFutgt[0], KFutgt[2]};
    Double_t KFvtgtkpxi[2] = {KFvtgt[0], KFvtgt[2]};
    std::vector<double> KFkpxi_res_x = {res_x[0], res_x[2]};
    std::vector<double> KFkpxi_res_y = {res_y[0], res_y[2]};
    std::vector<double> KFkpxi_res_u = {res_u[0], res_u[2]};
    std::vector<double> KFkpxi_res_v = {res_v[0], res_v[2]};

    TVector3 KFkpxi_prodvtx = Kinematics::MultitrackVertex(2, KFxtgtkpxi, KFytgtkpxi,
							   KFutgtkpxi, KFvtgtkpxi,
							   KFkpxi_res_x, KFkpxi_res_y,
							   KFkpxi_res_u, KFkpxi_res_v);
    event.KFprodvtx_x_kpxi = KFkpxi_prodvtx.x();
    event.KFprodvtx_y_kpxi = KFkpxi_prodvtx.y();
    event.KFprodvtx_z_kpxi = KFkpxi_prodvtx.z();

    TVector3 KFkpxipos_prodvtx; TVector3 KFkpximom_prodvtx;
    Double_t KFkpxitracklen_prodvtx; Double_t KFkpxitof_prodvtx;
    Bool_t KFkpxi_extrapolation_prodvtx = GFTrackCont.XiDecayToProdVertex(KFxi_id_container[best_xi],
									  KFkpxi_prodvtx,
									  KFkpxipos_prodvtx,
									  KFkpximom_prodvtx,
									  KFkpxitracklen_prodvtx,
									  KFkpxitof_prodvtx);
    if(KFkpxi_extrapolation_prodvtx){
      TVector3 dist = KFkpxipos_prodvtx - KFkpxi_prodvtx;
      event.KFxi_kpxiprodvtx_dist = dist.Mag();
      event.KFxi_kpxiprodvtx_x = KFkpxipos_prodvtx.x();
      event.KFxi_kpxiprodvtx_y = KFkpxipos_prodvtx.y();
      event.KFxi_kpxiprodvtx_z = KFkpxipos_prodvtx.z();
      event.KFxi_kpxiprodmom = KFkpximom_prodvtx.Mag();
      event.KFxi_kpxiprodmom_x = KFkpximom_prodvtx.x();
      event.KFxi_kpxiprodmom_y = KFkpximom_prodvtx.y();
      event.KFxi_kpxiprodmom_z = KFkpximom_prodvtx.z();
    }
#endif
    HF1( 10, event.MissMass[0] );
    HF1( 11, event.lmass );
    HF1( 12, event.ximass );
    HF1( 13, event.decays_mom[0] );
    HF1( 14, event.decays_mom[1] );
    HF1( 15, event.decays_mom[2] );
    HF1( 16, event.lmom );
    HF1( 17, event.ximom );
    HF1( 18, event.ppi_dist);
    HF1( 19, event.lpi_dist);

    HF1( 21, event.GFlmass);
    HF1( 22, event.GFximass);
    HF1( 23, event.GFdecays_mom[0]);
    HF1( 24, event.GFdecays_mom[1]);
    HF1( 25, event.GFdecays_mom[2]);
    HF1( 26, event.GFlmom);
    HF1( 27, event.GFximom);
    HF1( 28, event.GFppi_dist);
    HF1( 29, event.GFlpi_dist);
    HF1( 30, event.GFdecays_mass2[0]);
    HF1( 31, event.GFdecays_mass2[1]);
    HF1( 32, event.GFdecays_mass2[2]);
    HF1( 33, event.GFdecays_tof[0] - Kinematics::CalcTimeOfFlight(event.GFdecays_mom[0],
								  event.GFdecays_tracklen[0],
								  pdg::ProtonMass()));
    HF1( 34, event.GFdecays_tof[1] - Kinematics::CalcTimeOfFlight(event.GFdecays_mom[1],
								  event.GFdecays_tracklen[1],
								  pdg::PionMass()));
    HF1( 35, event.GFdecays_tof[2] - Kinematics::CalcTimeOfFlight(event.GFdecays_mom[2],
								  event.GFdecays_tracklen[2],
								  pdg::PionMass()));
    HF2( 36, event.GFdecays_mass2[0], event.GFdecays_mom[0]);
    HF2( 37, -event.GFdecays_mass2[1], event.GFdecays_mom[1]);
    HF2( 38, -event.GFdecays_mass2[2], event.GFdecays_mom[2]);

    //Remaining tracks
#if DebugDisp
    std::cout<<"Xi event, Save remaining tracks"<<std::endl;
#endif

    std::vector<Int_t> more_ppip_residuals;
    std::vector<Int_t> more_pim_residuals;
    std::vector<Int_t> more_pip_residuals;

    for(Int_t icombi=0; icombi<reconfailed_L_p_id_container.size(); ++icombi){
      Int_t p_id = reconfailed_L_p_id_container[icombi];
      Int_t pi_id = reconfailed_L_pi_id_container[icombi];

      if(p_id!=event.decays_id[0] && p_id!=event.decays_id[1] &&
	 p_id!=event.decays_id[2]){
	if(event.charge[p_id]==1){
	  if(std::find(target_p_id_container.begin(), target_p_id_container.end(), p_id)
	     == target_p_id_container.end() &&
	     std::find(target_ppip_id_container.begin(), target_ppip_id_container.end(), p_id)
	     == target_ppip_id_container.end()) more_ppip_residuals.push_back(p_id);
	}
	else{
	  if(std::find(target_pim_id_container.begin(), target_pim_id_container.end(), p_id)
	     == target_pim_id_container.end()) more_pim_residuals.push_back(p_id);
	}
      }

      if(pi_id!=event.decays_id[0] && pi_id!=event.decays_id[1] &&
	 pi_id!=event.decays_id[2]){
	if(event.charge[pi_id]==-1){
	  if(std::find(target_pim_id_container.begin(), target_pim_id_container.end(), pi_id)
	     == target_pim_id_container.end()) more_pim_residuals.push_back(pi_id);
	}
	else{
	  if(std::find(target_p_id_container.begin(), target_p_id_container.end(), pi_id)
	     == target_p_id_container.end() &&
	     std::find(target_ppip_id_container.begin(), target_ppip_id_container.end(), pi_id)
	     == target_ppip_id_container.end()) more_ppip_residuals.push_back(pi_id);
	}
      }
    }

    std::vector<Int_t> delete_combi;
    for(Int_t icombi=0; icombi<reconfailed_L_p_id_container.size(); ++icombi){
      Int_t p_id = reconfailed_L_p_id_container[icombi];
      Int_t pi_id = reconfailed_L_pi_id_container[icombi];
      if(p_id==event.decays_id[0] || p_id==event.decays_id[1] ||
	 p_id==event.decays_id[2] || pi_id==event.decays_id[0] ||
	 pi_id==event.decays_id[1] || pi_id==event.decays_id[2]) delete_combi.push_back(icombi);
    }

    for(Int_t icombi=0; icombi<delete_combi.size(); ++icombi){
      reconfailed_L_p_id_container.erase(reconfailed_L_p_id_container.begin()+delete_combi[icombi]-icombi);
      reconfailed_L_pi_id_container.erase(reconfailed_L_pi_id_container.begin()+delete_combi[icombi]-icombi);
      reconfailed_L_mass_container.erase(reconfailed_L_mass_container.begin()+delete_combi[icombi]-icombi);
      reconfailed_L_mom_container.erase(reconfailed_L_mom_container.begin()+delete_combi[icombi]-icombi);
      reconfailed_L_vtx_container.erase(reconfailed_L_vtx_container.begin()+delete_combi[icombi]-icombi);
      reconfailed_L_p_mom_container.erase(reconfailed_L_p_mom_container.begin()+delete_combi[icombi]-icombi);
      reconfailed_L_pi_mom_container.erase(reconfailed_L_pi_mom_container.begin()+delete_combi[icombi]-icombi);
      reconfailed_L_ppidist_container.erase(reconfailed_L_ppidist_container.begin()+delete_combi[icombi]-icombi);
    }

    event.ncombiLreconfailed = reconfailed_L_p_id_container.size();
    for(Int_t icombi=0; icombi<event.ncombiLreconfailed; ++icombi){
      event.pidLreconfailed.push_back(reconfailed_L_p_id_container[icombi]);
      event.piidLreconfailed.push_back(reconfailed_L_pi_id_container[icombi]);
      event.LmassLreconfailed.push_back(reconfailed_L_mass_container[icombi]);
      event.LdecayvtxLreconfailed_x.push_back(reconfailed_L_vtx_container[icombi].x());
      event.LdecayvtxLreconfailed_y.push_back(reconfailed_L_vtx_container[icombi].y());
      event.LdecayvtxLreconfailed_z.push_back(reconfailed_L_vtx_container[icombi].z());
      event.LmomLreconfailed.push_back(reconfailed_L_pi_mom_container[icombi].Mag());
      event.LmomLreconfailed.push_back(reconfailed_L_pi_mom_container[icombi].x());
      event.LmomLreconfailed.push_back(reconfailed_L_pi_mom_container[icombi].y());
      event.LmomLreconfailed.push_back(reconfailed_L_pi_mom_container[icombi].z());
      event.pmomLreconfailed.push_back(reconfailed_L_p_mom_container[icombi].Mag());
      event.pmomLreconfailed.push_back(reconfailed_L_p_mom_container[icombi].x());
      event.pmomLreconfailed.push_back(reconfailed_L_p_mom_container[icombi].y());
      event.pmomLreconfailed.push_back(reconfailed_L_p_mom_container[icombi].z());
      event.pimomLreconfailed.push_back(reconfailed_L_mom_container[icombi].Mag());
      event.pimomLreconfailed.push_back(reconfailed_L_mom_container[icombi].x());
      event.pimomLreconfailed.push_back(reconfailed_L_mom_container[icombi].y());
      event.pimomLreconfailed.push_back(reconfailed_L_mom_container[icombi].z());
      event.ppidistLreconfailed.push_back(reconfailed_L_ppidist_container[icombi]);
    }

    for(Int_t icombi=0; icombi<pipair_pip_id_container.size(); ++icombi){
      Int_t pip_id = pipair_pip_id_container[icombi];
      Int_t pim_id = pipair_pim_id_container[icombi];
      if(pip_id!=event.decays_id[0]){
	if(std::find(target_pip_id_container.begin(), target_pip_id_container.end(), pip_id) ==
	   target_pip_id_container.end()) more_pip_residuals.push_back(pip_id);
      }
      if(pim_id!=event.decays_id[1] && pim_id!=event.decays_id[2]){
	if(std::find(target_pim_id_container.begin(), target_pim_id_container.end(), pim_id) ==
	   target_pim_id_container.end()) more_pim_residuals.push_back(pim_id);
      }
    }

    delete_combi.clear();
    for(Int_t icombi=0; icombi<pipair_pip_id_container.size(); ++icombi){
      Int_t pip_id = pipair_pip_id_container[icombi];
      Int_t pim_id = pipair_pim_id_container[icombi];
      if(pip_id==event.decays_id[0] || pim_id==event.decays_id[1] ||
	 pim_id==event.decays_id[2]) delete_combi.push_back(icombi);
    }

    for(Int_t icombi=0; icombi<delete_combi.size(); ++icombi){
      pipair_pip_id_container.erase(pipair_pip_id_container.begin()+delete_combi[icombi]-icombi);
      pipair_pim_id_container.erase(pipair_pim_id_container.begin()+delete_combi[icombi]-icombi);
      pipair_reconL_mass_container.erase(pipair_reconL_mass_container.begin()+delete_combi[icombi]-icombi);
      pipair_mom_container.erase(pipair_mom_container.begin()+delete_combi[icombi]-icombi);
      pipair_pip_mom_container.erase(pipair_pip_mom_container.begin()+delete_combi[icombi]-icombi);
      pipair_pim_mom_container.erase(pipair_pim_mom_container.begin()+delete_combi[icombi]-icombi);
      pipair_pipidist_container.erase(pipair_pipidist_container.begin()+delete_combi[icombi]-icombi);
    }

    event.ncombiPipair = pipair_pip_id_container.size();
    for(Int_t icombi=0; icombi<event.ncombiPipair; ++icombi){
      event.pipidPipair.push_back(pipair_pip_id_container[icombi]);
      event.pimidPipair.push_back(pipair_pim_id_container[icombi]);
      event.pipmomPipair.push_back(pipair_pip_mom_container[icombi].Mag());
      event.pipmomPipair_x.push_back(pipair_pip_mom_container[icombi].x());
      event.pipmomPipair_y.push_back(pipair_pip_mom_container[icombi].y());
      event.pipmomPipair_z.push_back(pipair_pip_mom_container[icombi].z());
      event.pimmomPipair.push_back(pipair_pim_mom_container[icombi].Mag());
      event.pimmomPipair_x.push_back(pipair_pim_mom_container[icombi].x());
      event.pimmomPipair_y.push_back(pipair_pim_mom_container[icombi].y());
      event.pimmomPipair_z.push_back(pipair_pim_mom_container[icombi].z());
      event.momPipair.push_back(pipair_mom_container[icombi].Mag());
      event.momPipair_x.push_back(pipair_mom_container[icombi].x());
      event.momPipair_y.push_back(pipair_mom_container[icombi].y());
      event.momPipair_z.push_back(pipair_mom_container[icombi].z());
      event.reconLmassPipair.push_back(pipair_reconL_mass_container[icombi]);
      event.pipidistPipair.push_back(pipair_pipidist_container[icombi]);
    }

    for(Int_t icombi=0; icombi<more_ppip_residuals.size(); ++icombi){
      Int_t id = more_ppip_residuals[icombi];
      if(std::find(target_p_id_container.begin(), target_p_id_container.end(), id) == target_p_id_container.end() &&
	 std::find(target_ppip_id_container.begin(), target_ppip_id_container.end(), id) == target_ppip_id_container.end()){
	TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( id );
	if((event.pid[id]&4)==4 && (event.pid[id]&1)!=1 && event.charge[id]==1){ //proton
	  Double_t mass2 = qnan;
	  TVector3 mom(qnan, qnan, qnan);
	  target_p_id_container.push_back(id);
	  target_p_dist2tgt_container.push_back(tp -> GetClosestDist());

	  Int_t repid = 0;
	  Int_t flag = 1;
	  for(Int_t i=0;i<2;i++){
	    Int_t temp = flag&event.pid[id];
	    if(temp==flag) repid += 1;
	    flag*=2;
	  }
	  if(!GFTrackCont.TrackCheck(id, repid)){
	    target_p_mass2_container.push_back(mass2);
	    target_p_mom_container.push_back(mom);
	    continue;
	  }

	  Int_t htofhitid_p; Double_t tracklen_p;
	  Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	  Bool_t htofextrapolation_p = GFTrackCont.TPCHTOFTrackMatching(id, repid, tgtpos,
									event.HtofSeg, event.posHtof,
									htofhitid_p, tof, tracklen_p,
									pos, track2tgt_dist);
	  mom = GFTrackCont.GetMom(id, 0, repid);
	  if(htofextrapolation_p) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_p, event.tHtof[htofhitid_p]);
	  target_p_mass2_container.push_back(mass2);
	  target_p_mom_container.push_back(mom);
	}
	else{
	    Double_t mass2 = qnan;
	    TVector3 mom(qnan, qnan, qnan);
	    target_ppip_id_container.push_back(id);
	    target_ppip_dist2tgt_container.push_back(tp -> GetClosestDist());

	    Int_t repid = -1;
	    if(!GFTrackCont.TrackCheck(id, repid)){
	      target_ppip_mass2_container.push_back(mass2);
	      target_ppip_mom_container.push_back(mom);
	      continue;
	    }

	    Int_t htofhitid_ppi; Double_t tracklen_ppi;
	    Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	    Bool_t htofextrapolation = GFTrackCont.TPCHTOFTrackMatching(id, repid, tgtpos,
									event.HtofSeg, event.posHtof,
									htofhitid_ppi, tof, tracklen_ppi,
									pos, track2tgt_dist);
	    mom = GFTrackCont.GetMom(id, 0, repid);
	    if(htofextrapolation) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_ppi, event.tHtof[htofhitid_ppi]);
	    target_ppip_mass2_container.push_back(mass2);
	    target_ppip_mom_container.push_back(mom);
	}
      }
    }

    for(Int_t icombi=0; icombi<more_pim_residuals.size(); ++icombi){
      Int_t id = more_pim_residuals[icombi];
      if(std::find(target_pim_id_container.begin(), target_pim_id_container.end(), id) ==
	 target_pim_id_container.end()){
	TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( id );

	Double_t mass2 = qnan;
	TVector3 mom(qnan, qnan, qnan);
	target_pim_id_container.push_back(id);
	target_pim_dist2tgt_container.push_back(tp -> GetClosestDist());

	Int_t repid = 0;
	if(!GFTrackCont.TrackCheck(id, repid)){
	  target_pim_mass2_container.push_back(mass2);
	  target_pim_mom_container.push_back(mom);
	  continue;
	}

	Int_t htofhitid_pi; Double_t tracklen_pi;
	Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(id, repid, tgtpos,
								       event.HtofSeg, event.posHtof,
								       htofhitid_pi, tof, tracklen_pi,
								       pos, track2tgt_dist);
	mom = GFTrackCont.GetMom(id, 0, repid);
	if(htofextrapolation_pi) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_pi, event.tHtof[htofhitid_pi]);
	//if(htofextrapolation_pi && mass2 > 0.25) continue;
	target_pim_mass2_container.push_back(mass2);
	target_pim_mom_container.push_back(mom);
      }
    }

    for(Int_t icombi=0; icombi<more_pip_residuals.size(); ++icombi){
      Int_t id = more_pip_residuals[icombi];
      if(std::find(target_pip_id_container.begin(), target_pip_id_container.end(), id) ==
	 target_pip_id_container.end()){
	TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( id );

	Double_t mass2 = qnan;
	TVector3 mom(qnan, qnan, qnan);
	target_pip_id_container.push_back(id);
	target_pip_dist2tgt_container.push_back(tp -> GetClosestDist());

	Int_t repid = 0;
	if(!GFTrackCont.TrackCheck(id, repid)){
	  target_pip_mass2_container.push_back(mass2);
	  target_pip_mom_container.push_back(mom);
	  continue;
	}

	Int_t htofhitid_pi; Double_t tracklen_pi;
	Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(id, repid, tgtpos,
								       event.HtofSeg, event.posHtof,
								       htofhitid_pi, tof, tracklen_pi,
								       pos, track2tgt_dist);
	mom = GFTrackCont.GetMom(id, 0, repid);
	if(htofextrapolation_pi) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_pi, event.tHtof[htofhitid_pi]);
	target_pip_mass2_container.push_back(mass2);
	target_pip_mom_container.push_back(mom);
      }
    }

    event.pip_multi = target_pip_id_container.size();
    event.p_multi = target_p_id_container.size();
    event.pim_multi = target_pim_id_container.size();
    event.ppip_multi = target_ppip_id_container.size();
    event.accident_multi = target_accidental_id_container.size();
    for(Int_t itp=0; itp<target_p_id_container.size(); ++itp){
      Int_t id_p = target_p_id_container[itp];
      if(id_p==event.decays_id[0]){
	event.p_multi -= 1;
	continue;
      }
      event.residual_id.push_back(id_p);
      event.residual_dist2tgt.push_back(target_p_dist2tgt_container[itp]);
      event.residual_mass2.push_back(target_p_mass2_container[itp]);
      event.residual_mom.push_back(target_p_mom_container[itp].Mag());
      event.residual_mom_x.push_back(target_p_mom_container[itp].x());
      event.residual_mom_y.push_back(target_p_mom_container[itp].y());
      event.residual_mom_z.push_back(target_p_mom_container[itp].z());
      event.residual_charge.push_back(1);
      HF1( 1120, target_p_mom_container[itp].Mag());
    }
    for(Int_t itpip=0; itpip<target_pip_id_container.size(); ++itpip){
      Int_t id_pip = target_pip_id_container[itpip];
      if(id_pip==event.decays_id[0]){
	event.pip_multi -= 1;
	continue;
      }
      event.residual_id.push_back(id_pip);
      event.residual_dist2tgt.push_back(target_pip_dist2tgt_container[itpip]);
      event.residual_mass2.push_back(target_pip_mass2_container[itpip]);
      event.residual_mom.push_back(target_pip_mom_container[itpip].Mag());
      event.residual_mom_x.push_back(target_pip_mom_container[itpip].x());
      event.residual_mom_y.push_back(target_pip_mom_container[itpip].y());
      event.residual_mom_z.push_back(target_pip_mom_container[itpip].z());
      event.residual_charge.push_back(1);
      HF1( 1121, target_pip_mom_container[itpip].Mag());
    }
    for(Int_t itpim=0; itpim<target_pim_id_container.size(); ++itpim){
      Int_t id_pim = target_pim_id_container[itpim];
      if(id_pim==event.decays_id[1] || id_pim==event.decays_id[2]){
	event.pim_multi -= 1;
	continue;
      }
      event.residual_id.push_back(id_pim);
      event.residual_dist2tgt.push_back(target_pim_dist2tgt_container[itpim]);
      event.residual_mass2.push_back(target_pim_mass2_container[itpim]);
      event.residual_mom.push_back(target_pim_mom_container[itpim].Mag());
      event.residual_mom_x.push_back(target_pim_mom_container[itpim].x());
      event.residual_mom_y.push_back(target_pim_mom_container[itpim].y());
      event.residual_mom_z.push_back(target_pim_mom_container[itpim].z());
      event.residual_charge.push_back(-1);
      HF1( 1122, target_pim_mom_container[itpim].Mag());
    }
    for(Int_t itppip=0; itppip<target_ppip_id_container.size(); ++itppip){
      Int_t id_ppip = target_ppip_id_container[itppip];
      if(id_ppip==event.decays_id[0]){
	event.ppip_multi -= 1;
	continue;
      }
      event.residual_id.push_back(id_ppip);
      event.residual_dist2tgt.push_back(target_ppip_dist2tgt_container[itppip]);
      event.residual_mass2.push_back(target_ppip_mass2_container[itppip]);
      event.residual_mom.push_back(target_ppip_mom_container[itppip].Mag());
      event.residual_mom_x.push_back(target_ppip_mom_container[itppip].x());
      event.residual_mom_y.push_back(target_ppip_mom_container[itppip].y());
      event.residual_mom_z.push_back(target_ppip_mom_container[itppip].z());
      event.residual_charge.push_back(event.charge[id_ppip]);
    }
    event.residual_multi = event.pip_multi + event.p_multi + event.pim_multi + event.ppip_multi + event.accident_multi;

    HF1( 1020, event.pim_multi);
    HF1( 1021, event.p_multi);
    HF1( 1022, event.pip_multi);

    //find Xi-p event
    if(event.residual_multi==1 && event.p_multi==1){
      HF1( 500, -event.BE[0]);
    }

    HF1( 120, -event.BE[0]);
    HF1( 122, -event.BETPC[0]);
    int id_p = GFxi_p_id_container.at(best_xi);
    int id_pi1 = GFxi_pi_id_container.at(best_xi);
    int id_pi2 = GFxi_pi2_id_container.at(best_xi);
    event.pi1tid = id_pi1;
    event.pi2tid = id_pi2;
    int G4ptid = G4TrackID.at(id_p);
    int G4pi1tid = G4TrackID.at(id_pi1);
    int G4pi2tid = G4TrackID.at(id_pi2);
    int G4ptnh = PureHits.at(id_p);
    int G4pi1tnh = PureHits.at(id_pi1);
    int G4pi2tnh = PureHits.at(id_pi2);
    event.p1tid = id_p;
    if(event.G4p1id == G4ptid and event.G4pi1id == G4pi1tid){
      event.l1good = true;
    }
    if(event.G4p2id == G4ptid and event.G4pi2id == G4pi1tid){
      event.l2good = true;
      event.llswap = true;
    }
    if((event.G4p1id == G4ptid and event.G4pi2id == G4pi1tid ) or
       (event.G4p2id == G4ptid and event.G4pi1id == G4pi1tid)){
      event.piswap = true;
    }
    return true;
  } //if(event.xiflag)

#if DebugDisp
  std::cout<<"6. Single L searching starts (No Xi-, LL)"<<std::endl;
#endif
  Int_t GFtrackid_decays[2] = {-1, -1}; Int_t GFrepid_decays[2];

  Int_t best_l = -1; Double_t prev_massdiff = 9999.;
  for(Int_t id=0; id<l_candidates; ++id){
    if(TMath::IsNaN(GFL_mass_container[id])) continue;
    if(TMath::IsNaN(GFL_ppidist_container[id])) continue; //Genfit's fitting was succeeded.
    if(L_targetdist_container[id] > ltarget_distcut) continue; //Select Lambda from the traget
    if(GFL_targetdist_container[id] > GFltarget_distcut) continue;

    event.lflag = true;
    Double_t diff = TMath::Abs(GFL_mass_container[id] - LambdaMass);
    if(prev_massdiff > diff){
      prev_massdiff = diff;
      best_l = id;
    }
  }

  if(event.lflag){
    Int_t id = best_l;
    event.lmass = L_mass_container[id];
    event.ldecayvtx_x = L_vtx_container[id].x();
    event.ldecayvtx_y = L_vtx_container[id].y();
    event.ldecayvtx_z = L_vtx_container[id].z();
    event.lmom = L_mom_container[id].Mag();
    event.lmom_x = L_mom_container[id].x();
    event.lmom_y = L_mom_container[id].y();
    event.lmom_z = L_mom_container[id].z();
    event.ppi_dist = L_ppidist_container[id];
    event.ltarget_dist = L_targetdist_container[id];
    event.ltargetvtx_x = L_targetvtx_container[id].x();
    event.ltargetvtx_y = L_targetvtx_container[id].y();
    event.ltargetvtx_z = L_targetvtx_container[id].z();

    event.decays_mom.push_back(L_p_mom_container[id].Mag());
    event.decays_mom.push_back(L_pi_mom_container[id].Mag());
    event.decays_mom_x.push_back(L_p_mom_container[id].x());
    event.decays_mom_x.push_back(L_pi_mom_container[id].x());
    event.decays_mom_y.push_back(L_p_mom_container[id].y());
    event.decays_mom_y.push_back(L_pi_mom_container[id].y());
    event.decays_mom_z.push_back(L_p_mom_container[id].z());
    event.decays_mom_z.push_back(L_pi_mom_container[id].z());
    event.decays_id.push_back(L_p_id_container[id]);
    event.decays_id.push_back(L_pi_id_container[id]);

    TLorentzVector Lv_p(L_p_mom_container[id],
			TMath::Hypot(L_p_mom_container[id].Mag(), ProtonMass));
    Lv_p.Boost(-L_mom_container[id]);
    event.decays_CMmom.push_back(Lv_p.P());
    event.decays_CMmom_x.push_back(Lv_p.Px());
    event.decays_CMmom_y.push_back(Lv_p.Py());
    event.decays_CMmom_z.push_back(Lv_p.Pz());

    TLorentzVector Lv_pi(L_pi_mom_container[id],
			 TMath::Hypot(L_pi_mom_container[id].Mag(), PionMass));
    Lv_pi.Boost(-L_mom_container[id]);
    event.decays_CMmom.push_back(Lv_pi.P());
    event.decays_CMmom_x.push_back(Lv_pi.Px());
    event.decays_CMmom_y.push_back(Lv_pi.Py());
    event.decays_CMmom_z.push_back(Lv_pi.Pz());

    event.GFlmass = GFL_mass_container[id];
    event.GFldecayvtx_x = GFL_vtx_container[id].x();
    event.GFldecayvtx_y = GFL_vtx_container[id].y();
    event.GFldecayvtx_z = GFL_vtx_container[id].z();
    event.GFlmom = GFL_mom_container[id].Mag();
    event.GFlmom_x = GFL_mom_container[id].x();
    event.GFlmom_y = GFL_mom_container[id].y();
    event.GFlmom_z = GFL_mom_container[id].z();
    event.GFppi_dist = GFL_ppidist_container[id];
    event.GFltarget_dist = GFL_targetdist_container[id];
    event.GFltargetvtx_x = GFL_targetvtx_container[id].x();
    event.GFltargetvtx_y = GFL_targetvtx_container[id].y();
    event.GFltargetvtx_z = GFL_targetvtx_container[id].z();
    event.GFltargetcenter_dist = GFL_targetcenterdist_container[id];
    event.GFltargetcenter_x = GFL_targetcentervtx_container[id].x();
    event.GFltargetcenter_y = GFL_targetcentervtx_container[id].y();
    event.GFltargetcenter_z = GFL_targetcentervtx_container[id].z();

    GFtrackid_decays[0] = GFL_p_id_container[id];
    GFtrackid_decays[1] = GFL_pi_id_container[id];
    GFrepid_decays[0] = GFL_p_repid_container[id];
    GFrepid_decays[1] = GFL_pi_repid_container[id];

    const Int_t ntrack_l = 3;
    Double_t x0[ntrack_l] = {event.xtgtK18[0], event.xtgtTPCKurama[0], event.GFltargetcenter_x};
    Double_t y0[ntrack_l] = {event.ytgtK18[0], event.ytgtTPCKurama[0], event.GFltargetcenter_y};
    Double_t u0[ntrack_l] = {event.utgtK18[0], event.utgtTPCKurama[0], event.GFlmom_x/event.GFlmom_z};
    Double_t v0[ntrack_l] = {event.vtgtK18[0], event.vtgtTPCKurama[0], event.GFlmom_y/event.GFlmom_z};
    TVector3 kk_l_vertex = Kinematics::MultitrackVertex(ntrack_l, x0, y0, u0, v0);

    event.GFprodvtx_x_l = kk_l_vertex.x();
    event.GFprodvtx_y_l = kk_l_vertex.y();
    event.GFprodvtx_z_l = kk_l_vertex.z();
    Double_t prodvtx_closedist = qnan;
    TVector3 prodvtx_closest = Kinematics::CalcCloseDistLambda(kk_l_vertex,
							       GFL_vtx_container[id],
							       GFL_mom_container[id],
							       prodvtx_closedist);

    event.GFlprodvtx_x = prodvtx_closest.x();
    event.GFlprodvtx_y = prodvtx_closest.y();
    event.GFlprodvtx_z = prodvtx_closest.z();
    event.GFlprodvtx_dist = prodvtx_closedist;

    TVector3 lambda_tracklen = kk_l_vertex - GFL_vtx_container[id];
    event.GFltracklen = lambda_tracklen.Mag();
    event.GFltof = Kinematics::CalcTimeOfFlight(event.GFlmom, lambda_tracklen.Mag(), pdg::LambdaMass());

    event.GFdecays_htofid.push_back(GFL_p_htofid_container[id]);
    event.GFdecays_htofid.push_back(GFL_pi_htofid_container[id]);
    event.GFdecays_tracklen.push_back(GFL_p_tracklen_container[id]);
    event.GFdecays_tracklen.push_back(GFL_pi_tracklen_container[id]);
    event.GFdecays_tof.push_back(GFL_p_tof_container[id]);
    event.GFdecays_tof.push_back(GFL_pi_tof_container[id]);
    event.GFdecays_mass2.push_back(GFL_p_mass2_container[id]);
    event.GFdecays_mass2.push_back(GFL_pi_mass2_container[id]);

    event.GFdecays_mom.push_back(GFL_p_mom_container[id].Mag());
    event.GFdecays_mom.push_back(GFL_pi_mom_container[id].Mag());
    event.GFdecays_mom_x.push_back(GFL_p_mom_container[id].x());
    event.GFdecays_mom_x.push_back(GFL_pi_mom_container[id].x());
    event.GFdecays_mom_y.push_back(GFL_p_mom_container[id].y());
    event.GFdecays_mom_y.push_back(GFL_pi_mom_container[id].y());
    event.GFdecays_mom_z.push_back(GFL_p_mom_container[id].z());
    event.GFdecays_mom_z.push_back(GFL_pi_mom_container[id].z());
    event.GFmomloss.push_back(qnan);
    event.GFmomloss.push_back(qnan);
    event.GFeloss.push_back(qnan);
    event.GFeloss.push_back(qnan);

    TLorentzVector GFLv_p(GFL_p_mom_container[id],
			  TMath::Hypot(GFL_p_mom_container[id].Mag(), ProtonMass));
    GFLv_p.Boost(-GFL_mom_container[id]);
    event.GFdecays_CMmom.push_back(GFLv_p.P());
    event.GFdecays_CMmom_x.push_back(GFLv_p.Px());
    event.GFdecays_CMmom_y.push_back(GFLv_p.Py());
    event.GFdecays_CMmom_z.push_back(GFLv_p.Pz());

    TLorentzVector GFLv_pi(GFL_pi_mom_container[id],
			   TMath::Hypot(GFL_pi_mom_container[id].Mag(), PionMass));
    GFLv_pi.Boost(-GFL_mom_container[id]);
    event.GFdecays_CMmom.push_back(GFLv_pi.P());
    event.GFdecays_CMmom_x.push_back(GFLv_pi.Px());
    event.GFdecays_CMmom_y.push_back(GFLv_pi.Py());
    event.GFdecays_CMmom_z.push_back(GFLv_pi.Pz());

    int trackid_p1 = GFL_p_id_container[id];
    int trackid_pi1 = GFL_pi_id_container[id];
#if KinematicFit
    auto track_p1 = TPCAna.GetTrackTPCHelix(trackid_p1);
    auto Vp1 = track_p1->GetCovarianceMatrix();
    auto track_pi1 = TPCAna.GetTrackTPCHelix(trackid_pi1);
    auto Vpi1 = track_pi1->GetCovarianceMatrix();
    double Diag_ppi1[6]={
      Vp1(0,0),Vp1(1,1),Vp1(2,2),Vpi1(0,0),Vpi1(1,1),Vpi1(2,2)
    };
    auto Offdiag_ppi1 = MathTools::MergeOffdiagonals(Vp1, Vpi1);
    TVector3 HTVP1(GFL_p_mom_container[id].x(),
		   GFL_p_mom_container[id].z(),
		   GFL_p_mom_container[id].y());
    TVector3 HTVPi1(GFL_pi_mom_container[id].x(),
		    GFL_pi_mom_container[id].z(),
		    GFL_pi_mom_container[id].y());
    TVector3 HTVLd1 = HTVP1+HTVPi1;
    TLorentzVector HLVP1(HTVP1, TMath::Hypot(HTVP1.Mag(), pdg::ProtonMass()));
    TLorentzVector HLVPi1(HTVPi1, TMath::Hypot(HTVPi1.Mag(), pdg::PionMass()));
    TLorentzVector HLVLd1(HTVLd1, TMath::Hypot(HTVLd1.Mag(), pdg::LambdaMass()));
    Double_t KFchisqrl1=-1;
    Double_t KFpvall1=-1;
    FourVectorFitter KFLd1(HLVP1, HLVPi1, HLVLd1);
    KFLd1.SetInvMass(LambdaMass);
    KFLd1.SetMaximumStep(5);
    KFLd1.SetVariance(Diag_ppi1);
    KFLd1.AddOffdiagonals(Offdiag_ppi1);
    KFchisqrl1 = KFLd1.DoKinematicFit();
    KFpvall1 = KFLd1.GetPValue();
    auto HcontLd1 = KFLd1.GetFittedLV();
    auto PullLd1 = KFLd1.GetPull();
    auto KFHLVP1 = HcontLd1.at(0);
    auto KFHLVPi1 = HcontLd1.at(1);
    auto KFHLVLd1 = HcontLd1.at(2);
    auto KFTVP1 = TVector3(KFHLVP1.X(),KFHLVP1.Z(),KFHLVP1.Y());
    auto KFTVPi1 = TVector3(KFHLVPi1.X(),KFHLVPi1.Z(),KFHLVPi1.Y());
    auto KFTVLd1 = TVector3(KFHLVLd1.X(),KFHLVLd1.Z(),KFHLVLd1.Y());
    auto VLd1 = KFLd1.GetUnmeasuredCovariance();

    Double_t KFl_targetcenter_dist;
    TVector3 KFl_pos_tgtcenter = Kinematics::LambdaTargetCenter(GFL_vtx_container[id],
								KFTVLd1,
								KFl_targetcenter_dist);

    Double_t KFx0[ntrack_l] = {event.xtgtK18[0], event.xtgtTPCKurama[0],
			       KFl_pos_tgtcenter.x()};
    Double_t KFy0[ntrack_l] = {event.ytgtK18[0], event.ytgtTPCKurama[0],
			       KFl_pos_tgtcenter.y()};
    Double_t KFu0[ntrack_l] = {event.utgtK18[0], event.utgtTPCKurama[0],
			       KFTVP1.x()/KFTVP1.z()};
    Double_t KFv0[ntrack_l] = {event.vtgtK18[0], event.vtgtTPCKurama[0],
			       KFTVP1.y()/KFTVP1.z()};
    double resl1_u,resl1_v;
    MathTools::DecomposeResolutionUV(VLd1, KFTVLd1, resl1_u, resl1_v);
    std::vector<double> res_x = {res_xKurama, res_xK18, res_xLdVtx};
    std::vector<double> res_y = {res_yKurama, res_yK18, res_yLdVtx};
    std::vector<double> res_u = {res_uKurama, res_uK18, resl1_u};
    std::vector<double> res_v = {res_vKurama, res_vK18, resl1_v};

    Double_t chisqr_kkl;
    TVector3 KFkk_l_vertex = Kinematics::MultitrackVertex(ntrack_l,
							  KFx0, KFy0,
							  KFu0, KFv0,
							  res_x, res_y,
							  res_u, res_v,
							  chisqr_kkl);
    Double_t KFprodvtx_closedist = qnan;
    TVector3 KFprodvtx_closest = Kinematics::CalcCloseDistLambda(KFkk_l_vertex,
								 GFL_vtx_container[id],
								 KFTVLd1,
								 KFprodvtx_closedist);

    event.KFlchisqr = chisqr_kkl;
    event.KFlmom = KFTVLd1.Mag();
    event.KFlmom_x = KFTVLd1.x();
    event.KFlmom_y = KFTVLd1.y();
    event.KFlmom_z = KFTVLd1.z();

    event.KFprodvtx_x_l = KFkk_l_vertex.x();
    event.KFprodvtx_y_l = KFkk_l_vertex.y();
    event.KFprodvtx_z_l = KFkk_l_vertex.z();

    event.KFlprodvtx_x = KFprodvtx_closest.x();
    event.KFlprodvtx_y = KFprodvtx_closest.y();
    event.KFlprodvtx_z = KFprodvtx_closest.z();
    event.KFlprodvtx_dist = KFprodvtx_closedist;

    TVector3 KFlambda_tracklen = KFkk_l_vertex - GFL_vtx_container[id];
    event.KFltracklen = KFlambda_tracklen.Mag();
    event.KFltof = Kinematics::CalcTimeOfFlight(KFTVLd1.Mag(),
						KFlambda_tracklen.Mag(),
						pdg::LambdaMass());

    TLorentzVector KFLv_p1(KFTVP1, TMath::Hypot(KFTVP1.Mag(), ProtonMass));
    TLorentzVector KFLv_pi1(KFTVPi1, TMath::Hypot(KFTVPi1.Mag(), PionMass));
    event.KFdecays_mom.push_back(KFLv_p1.P());
    event.KFdecays_mom_x.push_back(KFLv_p1.Px());
    event.KFdecays_mom_y.push_back(KFLv_p1.Py());
    event.KFdecays_mom_z.push_back(KFLv_p1.Pz());
    event.KFdecays_mom.push_back(KFLv_pi1.P());
    event.KFdecays_mom_x.push_back(KFLv_pi1.Px());
    event.KFdecays_mom_y.push_back(KFLv_pi1.Py());
    event.KFdecays_mom_z.push_back(KFLv_pi1.Pz());

    KFLv_p1.Boost(-KFTVLd1);
    KFLv_pi1.Boost(-KFTVLd1);
    event.KFdecays_CMmom.push_back(KFLv_p1.P());
    event.KFdecays_CMmom_x.push_back(KFLv_p1.Px());
    event.KFdecays_CMmom_y.push_back(KFLv_p1.Py());
    event.KFdecays_CMmom_z.push_back(KFLv_p1.Pz());
    event.KFdecays_CMmom.push_back(KFLv_pi1.P());
    event.KFdecays_CMmom_x.push_back(KFLv_pi1.Px());
    event.KFdecays_CMmom_y.push_back(KFLv_pi1.Py());
    event.KFdecays_CMmom_z.push_back(KFLv_pi1.Pz());
#endif
    //order : p, pi
    event.GFntdecays = 2;
    event.GFnhtrack.resize(event.GFntdecays);
    event.GFchisqr.resize(event.GFntdecays);
    event.GFcharge.resize(event.GFntdecays);
    event.GFtof.resize(event.GFntdecays);
    event.GFtracklen.resize(event.GFntdecays);
    event.GFpval.resize(event.GFntdecays);
    event.GFpdgcode.resize(event.GFntdecays);

    for(int j=0;j<event.GFntdecays;j++){
      Int_t igf = GFtrackid_decays[j];
      Int_t repid = GFrepid_decays[j];
      event.GFnhtrack[j] = GFTrackCont.GetNHits(igf);
      event.GFchisqr[j] = GFTrackCont.GetChi2NDF(igf, repid);
      event.GFcharge[j] = GFTrackCont.GetCharge(igf, repid);
      event.GFtof[j] = GFTrackCont.GetTrackTOF(igf, 0, -1, repid);
      event.GFtracklen[j] = GFTrackCont.GetTrackLength(igf, 0, -1, repid);
      event.GFpval[j] = GFTrackCont.GetPvalue(igf, repid);
      event.GFpdgcode[j] = GFTrackCont.GetPDGcode(igf, repid);
    } //GFntdecays
    HF1( 130, -event.BE[0]);
    HF1( 132, -event.BETPC[0]);

    {
      auto iter = find(target_p_id_container.begin(), target_p_id_container.end(), GFtrackid_decays[0]);
      if(iter != target_p_id_container.end()){
	target_p_id_container.erase(iter);
      }
    }
    {
      auto iter = find(target_ppip_id_container.begin(), target_ppip_id_container.end(), GFtrackid_decays[0]);
      if(iter != target_ppip_id_container.end()){
	target_ppip_id_container.erase(iter);
      }
    }
    {
      auto iter = find(target_pim_id_container.begin(), target_pim_id_container.end(), GFtrackid_decays[1]);
      if(iter != target_pim_id_container.end()){
	target_pim_id_container.erase(iter);
      }
    }
    {
      auto iter = find(target_k_id_container.begin(), target_k_id_container.end(), GFtrackid_decays[1]);
      if(iter != target_k_id_container.end()){
	target_k_id_container.erase(iter);
      }
    }

    //Lphi searching
    if(target_k_id_container.size()==1 &&
       target_p_id_container.size()==0 && target_pip_id_container.size()==0 &&
       target_ppip_id_container.size()==0 && target_accidental_id_container.size()==0){

#if DebugDisp
      std::cout<<"optional: L phi searching "<<std::endl;
#endif
      Int_t id_kp = -1;
      for(Int_t k=0;k<ntTpc;k++){
	if(event.isKurama[k]==1) id_kp = k;
      }

      Int_t id_km = target_k_id_container[0];
      Double_t km_mass2 = target_k_mass2_container[0];
      //Double_t km_mom_helix = target_k_mom_container[0];

      Bool_t pim_km_container_overlap = (find(target_pim_id_container.begin(), target_pim_id_container.end(), id_km) != target_pim_id_container.end());
      if((pim_km_container_overlap && target_pim_id_container.size()==1) || target_pim_id_container.size()==0){

	Double_t km_par[5];
	km_par[0] = event.helix_cx[id_km];
	km_par[1] = event.helix_cy[id_km];
	km_par[2] = event.helix_z0[id_km];
	km_par[3] = event.helix_r[id_km];
	km_par[4] = event.helix_dz[id_km];

	Int_t km_nh = event.helix_t[id_km].size();
	Double_t km_theta_min = TMath::Max(event.helix_t[id_km][0] - vtx_scan_rangeInsideL/km_par[3], event.helix_t[id_km][km_nh-1]);
	Double_t km_theta_max = event.helix_t[id_km][0] + vtx_scan_range/km_par[3];
	TVector3 km_start = TVector3(event.calpos_x[id_km][0], event.calpos_y[id_km][0], event.calpos_z[id_km][0]);
	TVector3 km_end = TVector3(event.calpos_x[id_km][km_nh-1], event.calpos_y[id_km][km_nh-1], event.calpos_z[id_km][km_nh-1]);

	if(id_kp!=-1 && (event.pid[id_kp]&2)==2 && event.charge[id_kp]==1){
	  Int_t kp_repid = 0;
	  Int_t flag = 1;
	  for(Int_t i=0;i<1;i++){
	    Int_t temp = flag&event.pid[id_kp];
	    if(temp==flag) kp_repid += 1;
	    flag*=2;
	  }

	  Double_t kp_par[5];
	  kp_par[0] = event.helix_cx[id_kp];
	  kp_par[1] = event.helix_cy[id_kp];
	  kp_par[2] = event.helix_z0[id_kp];
	  kp_par[3] = event.helix_r[id_kp];
	  kp_par[4] = event.helix_dz[id_kp];
	  Int_t kp_nh = event.helix_t[id_kp].size();
	  Double_t kp_theta_min = event.helix_t[id_kp][0] - vtx_scan_range/kp_par[3];
	  Double_t kp_theta_max = TMath::Min(event.helix_t[id_kp][0] + vtx_scan_rangeInsideL/kp_par[3], event.helix_t[id_kp][kp_nh-1]);
	  TVector3 kp_start = TVector3(event.calpos_x[id_kp][0], event.calpos_y[id_kp][0], event.calpos_z[id_kp][0]);
	  TVector3 kp_end = TVector3(event.calpos_x[id_kp][kp_nh-1], event.calpos_y[id_kp][kp_nh-1], event.calpos_z[id_kp][kp_nh-1]);

	  Double_t kk_dist = 10000.;
	  TVector3 kp_mom; TVector3 km_mom;
	  TVector3 phi_mom;
	  TVector3 phi_vert = Kinematics::LambdaVertex(dMagneticField, kp_par, km_par,
						       kp_theta_min, kp_theta_max,
						       km_theta_min, km_theta_max,
						       kp_mom, km_mom,
						       phi_mom, kk_dist);
	  phi_mom = km_mom + kp_mom;

	  TLorentzVector Lkp(kp_mom, TMath::Hypot(kp_mom.Mag(), KaonMass));
	  TLorentzVector Lkm(km_mom, TMath::Hypot(km_mom.Mag(), KaonMass));
	  TLorentzVector Lphi = Lkp + Lkm;

	  if(!TMath::IsNaN(kk_dist) && TMath::Abs(phi_vert.x()) < 30. && TMath::Abs(phi_vert.z() - tpc::ZTarget) < 30. && TMath::Abs(phi_vert.y()) < 30.){
	    if(GFTrackCont.TrackCheck(id_km, 1) && GFTrackCont.TrackCheck(id_kp, kp_repid)){
	      Double_t GFkk_dist = 10000.;
	      Double_t GFextrapolation_kk[2];
	      TVector3 GFmom_kk[2];
	      TVector3 GFphi_vert;
	      GFTrackCont.FindVertex(id_kp, id_km,
				     1, 1,
				     GFextrapolation_kk[0], GFextrapolation_kk[1],
				     GFmom_kk[0], GFmom_kk[1],
				     GFkk_dist, GFphi_vert,
				     vtx_scan_range);

	      TLorentzVector GFkm(GFmom_kk[0],
				  TMath::Hypot(GFmom_kk[0].Mag(), KaonMass));
	      TLorentzVector GFkp(GFmom_kk[1],
				  TMath::Hypot(GFmom_kk[1].Mag(), KaonMass));
	      TLorentzVector GFphi = GFkm + GFkp;
	      TVector3 GFmom_phi = GFmom_kk[0] + GFmom_kk[1];
	      TVector3 phivtx_dist(event.GFprodvtx_x_l - GFphi_vert.x(),
				   event.GFprodvtx_y_l - GFphi_vert.y(),
				   event.GFprodvtx_z_l - GFphi_vert.z());

	      event.lphiflag = true;
	      event.GFphimass = GFphi.M();
	      event.GFphidecayvtx_x = GFphi_vert.x();
	      event.GFphidecayvtx_y = GFphi_vert.y();
	      event.GFphidecayvtx_z = GFphi_vert.z();
	      event.GFphimom = GFmom_phi.Mag();
	      event.GFphimom_x = GFmom_phi.x();
	      event.GFphimom_y = GFmom_phi.y();
	      event.GFphimom_z = GFmom_phi.z();
	      event.GFkk_dist = GFkk_dist;
	      event.GFphiprodvtx_dist = phivtx_dist.Mag();

	      event.phimass = Lphi.M();
	      event.phidecayvtx_x = phi_vert.x();
	      event.phidecayvtx_y = phi_vert.y();
	      event.phidecayvtx_z = phi_vert.z();
	      event.phimom = phi_mom.Mag();
	      event.phimom_x = phi_mom.x();
	      event.phimom_y = phi_mom.y();
	      event.phimom_z = phi_mom.z();
	      event.kk_dist = kk_dist;
	      event.phi_km_mass2 = km_mass2;

	      event.phidecays_id.push_back(id_kp);
	      event.phidecays_id.push_back(id_km);
	      event.phidecays_mom.push_back(kp_mom.Mag());
	      event.phidecays_mom.push_back(km_mom.Mag());
	      event.phidecays_mom_x.push_back(kp_mom.x());
	      event.phidecays_mom_x.push_back(km_mom.x());
	      event.phidecays_mom_y.push_back(kp_mom.y());
	      event.phidecays_mom_y.push_back(km_mom.y());
	      event.phidecays_mom_z.push_back(kp_mom.z());
	      event.phidecays_mom_z.push_back(km_mom.z());

	      event.GFphidecays_mom.push_back(GFmom_kk[1].Mag());
	      event.GFphidecays_mom.push_back(GFmom_kk[0].Mag());
	      event.GFphidecays_mom_x.push_back(GFmom_kk[1].x());
	      event.GFphidecays_mom_x.push_back(GFmom_kk[0].x());
	      event.GFphidecays_mom_y.push_back(GFmom_kk[1].y());
	      event.GFphidecays_mom_y.push_back(GFmom_kk[0].y());
	      event.GFphidecays_mom_z.push_back(GFmom_kk[1].z());
	      event.GFphidecays_mom_z.push_back(GFmom_kk[0].z());
	    }
	  }
	}
      }
    } //lphiflag
    int id_p = L_p_id_container[id];
    int id_pi = L_pi_id_container[id];
    event.p1tid = id_p;
    event.pi1tid = id_pi;
    int G4ptid = G4TrackID.at(id_p);
    int G4pitid = G4TrackID.at(id_pi);
    if(event.G4p1id == G4ptid and event.G4pi1id == G4pitid){
      event.l1good = true;
    }
    if(event.G4p2id == G4ptid and event.G4pi2id == G4pitid){
      event.l2good = true;
    }
    if((event.G4p1id == G4ptid and event.G4pi2id == G4pitid)or
       (event.G4p2id == G4ptid and event.G4pi1id == G4pitid)){
      event.piswap = true;
    }
  } //lflag

  //Remaining tracks
#if DebugDisp
  std::cout<<"One L event or others, Save remaining tracks"<<std::endl;
#endif

  std::vector<Int_t> more_ppip_residuals;
  std::vector<Int_t> more_pim_residuals;
  std::vector<Int_t> more_pip_residuals;

  for(Int_t icombi=0; icombi<reconfailed_L_p_id_container.size(); ++icombi){
    Int_t p_id = reconfailed_L_p_id_container[icombi];
    Int_t pi_id = reconfailed_L_pi_id_container[icombi];

    if(p_id!=GFtrackid_decays[0] && p_id!=GFtrackid_decays[1]){
      if(event.charge[p_id]==1){
	if(std::find(target_p_id_container.begin(), target_p_id_container.end(), p_id)
	   == target_p_id_container.end() &&
	   std::find(target_ppip_id_container.begin(), target_ppip_id_container.end(), p_id)
	   == target_ppip_id_container.end()) more_ppip_residuals.push_back(p_id);
      }
      else{
	if(std::find(target_pim_id_container.begin(), target_pim_id_container.end(), p_id)
	   == target_pim_id_container.end()) more_pim_residuals.push_back(p_id);
      }
    }

    if(pi_id!=GFtrackid_decays[0] && pi_id!=GFtrackid_decays[1]){
      if(event.charge[pi_id]==-1){
	if(std::find(target_pim_id_container.begin(), target_pim_id_container.end(), pi_id)
	   == target_pim_id_container.end()) more_pim_residuals.push_back(pi_id);
      }
      else{
	if(std::find(target_p_id_container.begin(), target_p_id_container.end(), pi_id)
	   == target_p_id_container.end() &&
	   std::find(target_ppip_id_container.begin(), target_ppip_id_container.end(), pi_id)
	   == target_ppip_id_container.end()) more_ppip_residuals.push_back(pi_id);
      }
    }
  }

  std::vector<Int_t> delete_combi;
  for(Int_t icombi=0; icombi<reconfailed_L_p_id_container.size(); ++icombi){
    Int_t p_id = reconfailed_L_p_id_container[icombi];
    Int_t pi_id = reconfailed_L_pi_id_container[icombi];
    if(p_id==GFtrackid_decays[0] || p_id==GFtrackid_decays[1] ||
       pi_id==GFtrackid_decays[0] || pi_id==GFtrackid_decays[1]) delete_combi.push_back(icombi);
  }

  for(Int_t icombi=0; icombi<delete_combi.size(); ++icombi){
    reconfailed_L_p_id_container.erase(reconfailed_L_p_id_container.begin()+delete_combi[icombi]-icombi);
    reconfailed_L_pi_id_container.erase(reconfailed_L_pi_id_container.begin()+delete_combi[icombi]-icombi);
    reconfailed_L_mass_container.erase(reconfailed_L_mass_container.begin()+delete_combi[icombi]-icombi);
    reconfailed_L_mom_container.erase(reconfailed_L_mom_container.begin()+delete_combi[icombi]-icombi);
    reconfailed_L_vtx_container.erase(reconfailed_L_vtx_container.begin()+delete_combi[icombi]-icombi);
    reconfailed_L_p_mom_container.erase(reconfailed_L_p_mom_container.begin()+delete_combi[icombi]-icombi);
    reconfailed_L_pi_mom_container.erase(reconfailed_L_pi_mom_container.begin()+delete_combi[icombi]-icombi);
    reconfailed_L_ppidist_container.erase(reconfailed_L_ppidist_container.begin()+delete_combi[icombi]-icombi);
  }

  event.ncombiLreconfailed = reconfailed_L_p_id_container.size();
  for(Int_t icombi=0; icombi<event.ncombiLreconfailed; ++icombi){
    event.pidLreconfailed.push_back(reconfailed_L_p_id_container[icombi]);
    event.piidLreconfailed.push_back(reconfailed_L_pi_id_container[icombi]);
    event.LmassLreconfailed.push_back(reconfailed_L_mass_container[icombi]);
    event.LdecayvtxLreconfailed_x.push_back(reconfailed_L_vtx_container[icombi].x());
    event.LdecayvtxLreconfailed_y.push_back(reconfailed_L_vtx_container[icombi].y());
    event.LdecayvtxLreconfailed_z.push_back(reconfailed_L_vtx_container[icombi].z());
    event.LmomLreconfailed.push_back(reconfailed_L_pi_mom_container[icombi].Mag());
    event.LmomLreconfailed.push_back(reconfailed_L_pi_mom_container[icombi].x());
    event.LmomLreconfailed.push_back(reconfailed_L_pi_mom_container[icombi].y());
    event.LmomLreconfailed.push_back(reconfailed_L_pi_mom_container[icombi].z());
    event.pmomLreconfailed.push_back(reconfailed_L_p_mom_container[icombi].Mag());
    event.pmomLreconfailed.push_back(reconfailed_L_p_mom_container[icombi].x());
    event.pmomLreconfailed.push_back(reconfailed_L_p_mom_container[icombi].y());
    event.pmomLreconfailed.push_back(reconfailed_L_p_mom_container[icombi].z());
    event.pimomLreconfailed.push_back(reconfailed_L_mom_container[icombi].Mag());
    event.pimomLreconfailed.push_back(reconfailed_L_mom_container[icombi].x());
    event.pimomLreconfailed.push_back(reconfailed_L_mom_container[icombi].y());
    event.pimomLreconfailed.push_back(reconfailed_L_mom_container[icombi].z());
    event.ppidistLreconfailed.push_back(reconfailed_L_ppidist_container[icombi]);
  }

  for(Int_t icombi=0; icombi<pipair_pip_id_container.size(); ++icombi){
    Int_t pip_id = pipair_pip_id_container[icombi];
    Int_t pim_id = pipair_pim_id_container[icombi];
    if(pip_id!=GFtrackid_decays[0]){
      if(std::find(target_pip_id_container.begin(), target_pip_id_container.end(), pip_id) ==
	 target_pip_id_container.end()) more_pip_residuals.push_back(pip_id);
    }
    if(pim_id!=GFtrackid_decays[1]){
      if(std::find(target_pim_id_container.begin(), target_pim_id_container.end(), pim_id) ==
	 target_pim_id_container.end()) more_pim_residuals.push_back(pim_id);
    }
  }

  delete_combi.clear();
  for(Int_t icombi=0; icombi<pipair_pip_id_container.size(); ++icombi){
    Int_t pip_id = pipair_pip_id_container[icombi];
    Int_t pim_id = pipair_pim_id_container[icombi];
    if(pip_id==GFtrackid_decays[0] || pim_id==GFtrackid_decays[1]) delete_combi.push_back(icombi);
  }

  for(Int_t icombi=0; icombi<delete_combi.size(); ++icombi){
    pipair_pip_id_container.erase(pipair_pip_id_container.begin()+delete_combi[icombi]-icombi);
    pipair_pim_id_container.erase(pipair_pim_id_container.begin()+delete_combi[icombi]-icombi);
    pipair_reconL_mass_container.erase(pipair_reconL_mass_container.begin()+delete_combi[icombi]-icombi);
    pipair_mom_container.erase(pipair_mom_container.begin()+delete_combi[icombi]-icombi);
    pipair_pip_mom_container.erase(pipair_pip_mom_container.begin()+delete_combi[icombi]-icombi);
    pipair_pim_mom_container.erase(pipair_pim_mom_container.begin()+delete_combi[icombi]-icombi);
    pipair_pipidist_container.erase(pipair_pipidist_container.begin()+delete_combi[icombi]-icombi);
  }

  event.ncombiPipair = pipair_pip_id_container.size();
  for(Int_t icombi=0; icombi<event.ncombiPipair; ++icombi){
    event.pipidPipair.push_back(pipair_pip_id_container[icombi]);
    event.pimidPipair.push_back(pipair_pim_id_container[icombi]);
    event.pipmomPipair.push_back(pipair_pip_mom_container[icombi].Mag());
    event.pipmomPipair_x.push_back(pipair_pip_mom_container[icombi].x());
    event.pipmomPipair_y.push_back(pipair_pip_mom_container[icombi].y());
    event.pipmomPipair_z.push_back(pipair_pip_mom_container[icombi].z());
    event.pimmomPipair.push_back(pipair_pim_mom_container[icombi].Mag());
    event.pimmomPipair_x.push_back(pipair_pim_mom_container[icombi].x());
    event.pimmomPipair_y.push_back(pipair_pim_mom_container[icombi].y());
    event.pimmomPipair_z.push_back(pipair_pim_mom_container[icombi].z());
    event.momPipair.push_back(pipair_mom_container[icombi].Mag());
    event.momPipair_x.push_back(pipair_mom_container[icombi].x());
    event.momPipair_y.push_back(pipair_mom_container[icombi].y());
    event.momPipair_z.push_back(pipair_mom_container[icombi].z());
    event.reconLmassPipair.push_back(pipair_reconL_mass_container[icombi]);
    event.pipidistPipair.push_back(pipair_pipidist_container[icombi]);
  }

  for(Int_t icombi=0; icombi<more_ppip_residuals.size(); ++icombi){
    Int_t id = more_ppip_residuals[icombi];
    if(std::find(target_p_id_container.begin(), target_p_id_container.end(), id) == target_p_id_container.end() &&
       std::find(target_ppip_id_container.begin(), target_ppip_id_container.end(), id) == target_ppip_id_container.end()){
      TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( id );
      if((event.pid[id]&4)==4 && (event.pid[id]&1)!=1 && event.charge[id]==1){ //proton
	Double_t mass2 = qnan;
	TVector3 mom(qnan, qnan, qnan);
	target_p_id_container.push_back(id);
	target_p_dist2tgt_container.push_back(tp -> GetClosestDist());

	Int_t repid = 0;
	Int_t flag = 1;
	for(Int_t i=0;i<2;i++){
	  Int_t temp = flag&event.pid[id];
	  if(temp==flag) repid += 1;
	  flag*=2;
	}
	if(!GFTrackCont.TrackCheck(id, repid)){
	  target_p_mass2_container.push_back(mass2);
	  target_p_mom_container.push_back(mom);
	  continue;
	}

	Int_t htofhitid_p; Double_t tracklen_p;
	Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	Bool_t htofextrapolation_p = GFTrackCont.TPCHTOFTrackMatching(id, repid, tgtpos,
								      event.HtofSeg, event.posHtof,
								      htofhitid_p, tof, tracklen_p,
								      pos, track2tgt_dist);
	mom = GFTrackCont.GetMom(id, 0, repid);
	if(htofextrapolation_p) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_p, event.tHtof[htofhitid_p]);
	target_p_mass2_container.push_back(mass2);
	target_p_mom_container.push_back(mom);
      }
      else{
	Int_t repid = 0;
	if(!GFTrackCont.TrackCheck(id, repid)) continue;
	Int_t htofhitid_ppi; Double_t tracklen_ppi;
	Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	Bool_t htofextrapolation = GFTrackCont.TPCHTOFTrackMatching(id, repid, tgtpos,
								    event.HtofSeg, event.posHtof,
								    htofhitid_ppi, tof, tracklen_ppi,
								    pos, track2tgt_dist);
	Double_t mass2 = qnan;
	TVector3 mom = GFTrackCont.GetMom(id, 0, repid);
	if(htofextrapolation) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_ppi, event.tHtof[htofhitid_ppi]);
	target_ppip_id_container.push_back(id);
	target_ppip_mass2_container.push_back(mass2);
	target_ppip_mom_container.push_back(mom);
	target_ppip_dist2tgt_container.push_back(tp -> GetClosestDist());
      }
    }
  }

  for(Int_t icombi=0; icombi<more_pim_residuals.size(); ++icombi){
    Int_t id = more_pim_residuals[icombi];
    if(std::find(target_pim_id_container.begin(), target_pim_id_container.end(), id) ==
       target_pim_id_container.end()){
      TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( id );

      Double_t mass2 = qnan;
      TVector3 mom(qnan, qnan, qnan);
      target_pim_id_container.push_back(id);
      target_pim_dist2tgt_container.push_back(tp -> GetClosestDist());

      Int_t repid = 0;
      if(!GFTrackCont.TrackCheck(id, repid)){
	target_pim_mass2_container.push_back(mass2);
	target_pim_mom_container.push_back(mom);
	continue;
      }

      Int_t htofhitid_pi; Double_t tracklen_pi;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(id, repid, tgtpos,
								     event.HtofSeg, event.posHtof,
								     htofhitid_pi, tof, tracklen_pi,
								     pos, track2tgt_dist);
      mom = GFTrackCont.GetMom(id, 0, repid);
      if(htofextrapolation_pi) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_pi, event.tHtof[htofhitid_pi]);
      target_pim_mass2_container.push_back(mass2);
      target_pim_mom_container.push_back(mom);
    }
  }

  for(Int_t icombi=0; icombi<more_pip_residuals.size(); ++icombi){
    Int_t id = more_pip_residuals[icombi];
    if(std::find(target_pip_id_container.begin(), target_pip_id_container.end(), id) ==
       target_pip_id_container.end()){
      TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( id );

      Double_t mass2 = qnan;
      TVector3 mom(qnan, qnan, qnan);
      target_pip_id_container.push_back(id);
      target_pip_dist2tgt_container.push_back(tp -> GetClosestDist());

      Int_t repid = 0;
      if(!GFTrackCont.TrackCheck(id, repid)){
	target_pip_mass2_container.push_back(mass2);
	target_pip_mom_container.push_back(mom);
	continue;
      }

      Int_t htofhitid_pi; Double_t tracklen_pi;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(id, repid, tgtpos,
								     event.HtofSeg, event.posHtof,
								     htofhitid_pi, tof, tracklen_pi,
								     pos, track2tgt_dist);
      mom = GFTrackCont.GetMom(id, 0, repid);
      if(htofextrapolation_pi) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_pi, event.tHtof[htofhitid_pi]);
      target_pip_mass2_container.push_back(mass2);
      target_pip_mom_container.push_back(mom);
    }
  }

  std::vector<Double_t> twopi;
  event.pip_multi = target_pip_id_container.size();
  event.p_multi = target_p_id_container.size();
  event.pim_multi = target_pim_id_container.size();
  event.ppip_multi = target_ppip_id_container.size();
  event.accident_multi = target_accidental_id_container.size();
  for(Int_t itp=0; itp<target_p_id_container.size(); ++itp){
    Int_t id_p = target_p_id_container[itp];
    if(id_p==GFtrackid_decays[0]){
      event.p_multi -= 1;
      continue;
    }
    event.residual_id.push_back(id_p);
    event.residual_dist2tgt.push_back(target_p_dist2tgt_container[itp]);
    event.residual_mass2.push_back(target_p_mass2_container[itp]);
    event.residual_mom.push_back(target_p_mom_container[itp].Mag());
    event.residual_mom_x.push_back(target_p_mom_container[itp].x());
    event.residual_mom_y.push_back(target_p_mom_container[itp].y());
    event.residual_mom_z.push_back(target_p_mom_container[itp].z());
    event.residual_charge.push_back(1);
    if(event.lflag) HF1( 1130, target_p_mom_container[itp].Mag());
    else HF1( 1100, target_p_mom_container[itp].Mag());
  }

  for(Int_t itpip=0; itpip<target_pip_id_container.size(); ++itpip){
    Int_t id_pip = target_pip_id_container[itpip];
    if(id_pip==GFtrackid_decays[0]){
      event.pip_multi -= 1;
      continue;
    }
    event.residual_id.push_back(id_pip);
    event.residual_dist2tgt.push_back(target_pip_dist2tgt_container[itpip]);
    event.residual_mass2.push_back(target_pip_mass2_container[itpip]);
    event.residual_mom.push_back(target_pip_mom_container[itpip].Mag());
    event.residual_mom_x.push_back(target_pip_mom_container[itpip].x());
    event.residual_mom_y.push_back(target_pip_mom_container[itpip].y());
    event.residual_mom_z.push_back(target_pip_mom_container[itpip].z());
    event.residual_charge.push_back(1);
    if(event.lflag) HF1( 1131, target_pip_mom_container[itpip].Mag());
    else HF1( 1101, target_pip_mom_container[itpip].Mag());
  }

  for(Int_t itpim=0; itpim<target_pim_id_container.size(); ++itpim){
    Int_t id_pim = target_pim_id_container[itpim];
    if(id_pim==GFtrackid_decays[1]){
      event.pim_multi -= 1;
      continue;
    }
    double KE = 1000.*(TMath::Hypot(target_pim_mom_container[itpim].Mag(), PionMass) - PionMass); //MeV/c2
    if(KE < 60.) twopi.push_back(KE);
    event.residual_id.push_back(id_pim);
    event.residual_dist2tgt.push_back(target_pim_dist2tgt_container[itpim]);
    event.residual_mass2.push_back(target_pim_mass2_container[itpim]);
    event.residual_mom.push_back(target_pim_mom_container[itpim].Mag());
    event.residual_mom_x.push_back(target_pim_mom_container[itpim].x());
    event.residual_mom_y.push_back(target_pim_mom_container[itpim].y());
    event.residual_mom_z.push_back(target_pim_mom_container[itpim].z());
    event.residual_charge.push_back(-1);
    if(event.lflag) HF1( 1132, target_pim_mom_container[itpim].Mag());
    else HF1( 1102, target_pim_mom_container[itpim].Mag());
  }

  for(Int_t itppip=0; itppip<target_ppip_id_container.size(); ++itppip){
    Int_t id_ppip = target_ppip_id_container[itppip];
    if(id_ppip==GFtrackid_decays[0]){
      event.ppip_multi -= 1;
      continue;
    }
    event.residual_id.push_back(id_ppip);
    event.residual_dist2tgt.push_back(target_ppip_dist2tgt_container[itppip]);
    event.residual_mass2.push_back(target_ppip_mass2_container[itppip]);
    event.residual_mom.push_back(target_ppip_mom_container[itppip].Mag());
    event.residual_mom_x.push_back(target_ppip_mom_container[itppip].x());
    event.residual_mom_y.push_back(target_ppip_mom_container[itppip].y());
    event.residual_mom_z.push_back(target_ppip_mom_container[itppip].z());
    event.residual_charge.push_back(event.charge[id_ppip]);
  }

  event.residual_multi = event.pip_multi + event.p_multi + event.pim_multi + event.ppip_multi + event.accident_multi;

  if(event.lflag){
    HF1( 1030, event.pim_multi);
    HF1( 1031, event.p_multi);
    HF1( 1032, event.pip_multi);

    HF1( 130, -event.BE[0]);
    HF1( 132, -event.BETPC[0]);

    if(event.pim_multi>0){
      HF1( 134, -event.BE[0]);
      HF1( 136, -event.BETPC[0]);
    }
    if(event.pim_multi==1) event.lpiflag = true;
  }
  else if(twopi.size()>=2){
    for(int pi1 = 0;pi1<twopi.size();pi1++){
      Double_t prevKE = twopi[pi1];
      for(int pi2 = pi1+1;pi1<twopi.size();pi1++){
	Double_t currentKE = twopi[pi2];
	HF2( 501, TMath::Min(prevKE, currentKE), TMath::Max(prevKE, currentKE) );
	event.pipiflag = true;
      }
    }
    HF1( 140, -event.BE[0]);
    HF1( 142, -event.BETPC[0]);
  }
  else if(event.residual_multi==1 && event.pim_multi==1) event.pimflag = true;
  else if(event.residual_multi==0) event.emptyflag = true;
  else{
    HF1( 1000, event.pim_multi);
    HF1( 1001, event.p_multi);
    HF1( 1002, event.pip_multi);
  }

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
#if SaveHistograms
  /*
  HB1( 1, "Status", 21, 0., 21. );
  HB1( 3, "Genfit Fit Status", 2, 0., 2. );

  HB1( 10, "Missing Mass [KURAMA]; Missing mass (GeV/c^{2}); Counts (/10 MeV/#font[12]{c}^{2})", 320, 1., 1.8);
  HB1( 11," #Lambda Invariant Mass; M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); Counts (/0.002 GeV/#font[12]{c}^{2})", 80, 1.04, 1.2);
  HB1( 12," #Xi Invariant Mass; M_{#Lambda#pi^{-}} (GeV/#font[12]{c}^{2}); Counts (/0.002 GeV/#font[12]{c}^{2})", 125, 1.2 ,1.45);
  HB1( 13, "p Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 14, "#pi Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  HB1( 15, "#pi Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  HB1( 16, "#Lambda Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 17, "#Xi Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 18, "Closest Dist. p#pi;Distance [mm]; Counts (/0.1 mm)", 1000, 0., 100.);
  HB1( 19, "Closest Dist. #Lambda#pi;Distance [mm]; Counts (/0.1 mm)", 1000, 0., 100.);

  HB1( 21, "[GenFit] #Lambda Invariant Mass; M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 80, 1.04, 1.2);
  HB1( 22, "[GenFit] #Xi Invariant Mass; M_{#Lambda#pi^{-}} (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 125, 1.2 ,1.45);
  HB1( 23, "[GenFit] p_{#Lambda} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 24, "[GenFit] #pi_{#Lambda} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  HB1( 25, "[GenFit] #pi_{#Xi} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  HB1( 26, "[GenFit] #Lambda Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 27, "[GenFit] #Xi Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 28, "[GenFit] Closest Dist. p#pi;Distance [mm]; Counts (/0.1 mm)", 1000, 0., 100.);
  HB1( 29, "[GenFit] Closest Dist. #Lambda#pi;Distance [mm]; Counts (/0.1 mm)", 1000, 0., 100.);
  HB1( 30, "[GenFit] p_{#Lambda} Mass; M (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 350, 0.5, 1.2);
  HB1( 31, "[GenFit] pi_{#Lambda} Mass; M (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 250, 0.0, 0.5);
  HB1( 32, "[GenFit] pi_{#Xi} Mass; M (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 250, 0., 0.5);
  HB1( 33, "[GenFit] p_{#Lambda} TOF_{HTOF mt.} - TOF_{pi calc.}; TOF_{HTOF mt.} - TOF_{pi calc.} (ns);Counts [0.01 ns]", 1000, -5, 5);
  HB1( 34, "[GenFit] pi_{#Lambda} TOF_{HTOF mt.} - TOF_{pi calc.}; TOF_{HTOF mt.} - TOF_{pi calc.} (ns);Counts [0.01 ns]", 1000, -5, 5);
  HB1( 35, "[GenFit] pi_{#Xi} TOF_{HTOF mt.} - TOF_{pi calc.}; TOF_{HTOF mt.} - TOF_{pi calc.} (ns);Counts [0.01 ns]", 1000, -5, 5);
  HB2( 36, "[GenFit] p_{#Lambda} Mass/Charge;Mass/Charge (GeV/#font[12]{c}^{2}/q);Momentum (GeV/#font[12]{c})",1000, -1.0, 1.5, 1000, 0, 2.0);
  HB2( 37, "[GenFit] pi_{#Lambda} Mass/Charge;Mass/Charge (GeV/#font[12]{c}^{2}/q);Momentum (GeV/#font[12]{c})",1000, -1.0, 1.5, 1000, 0, 2.0);
  HB2( 38, "[GenFit] pi_{#Xi} Mass/Charge;Mass/Charge (GeV/#font[12]{c}^{2}/q);Momentum (GeV/#font[12]{c})",1000, -1.0, 1.5, 1000, 0, 2.0);

  HB1( 100, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 101, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 102, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 103, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);

  HB1( 110, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 111, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 112, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 113, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);

  HB1( 120, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 121, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 122, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 123, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);

  HB1( 130, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 131, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 132, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 133, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);

  HB1( 134, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 135, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 136, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 137, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);

  HB1( 140, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 141, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 142, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 143, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);

  HB1( 200, ";^{12}C(K^{-}, K^{+}#Xi^{-})X MM - M(^{11}B) (GeV/#font[12]{c}^{2});Counts (/5 MeV/#font[12]{c}^{2})", 80, -0.1, 0.3);

  HB1( 500, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/5 MeV/#font[12]{c}^{2})", 280, -400., 1000.);
  HB2( 501, ";T_{#piL} (MeV/#font[12]{c}); T_{#piH} (MeV/#font[12]{c})", 60, 0, 60, 60, 0, 60);
  HB2( 502, ";M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); ;M_{p#pi^{-}} (GeV/#font[12]{c}^{2})", 80, 1.04, 1.2, 80, 1.04, 1.2);
  HB2( 503, "[GenFit] L Vs L Invariant Mass;M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); ;M_{p#pi^{-}} (GeV/#font[12]{c}^{2})", 80, 1.04, 1.2, 80, 1.04, 1.2);
  HB1( 504, "#Lambda Invariant Mass; M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 80, 1.04, 1.2);
  HB1( 505, "[GenFit] #Lambda Invariant Mass; M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 80, 1.04, 1.2);

  HB1( 506, "p_{#Lambda} MassSquare; M^{2} (GeV/#font[12]{c}^{2})^{2}; Counts (/2 MeV/#font[12]{c}^{2})", 350, 0.5, 1.2);
  HB1( 507, "pi_{#Lambda} MassSquare; M^{2} (GeV/#font[12]{c}^{2})^{2}; Counts (/2 MeV/#font[12]{c}^{2})", 250, 0.0, 0.5);
  HB1( 508, "[GenFit] p_{#Lambda} MassSquare; M^{2} (GeV/#font[12]{c}^{2})^{2}; Counts (/2 MeV/#font[12]{c}^{2})", 350, 0.5, 1.2);
  HB1( 509, "[GenFit] pi_{#Lambda} MassSquare; M^{2} (GeV/#font[12]{c}^{2})^{2}; Counts (/2 MeV/#font[12]{c}^{2})", 250, 0.0, 0.5);

  HB1( 1000, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1001, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1002, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1100, "p Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 1101, "#pi^{+} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  HB1( 1102, "#pi^{-} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);

  HB1( 1010, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1011, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1012, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1110, "p Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 1111, "#pi^{+} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  HB1( 1112, "#pi^{-} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);

  HB1( 1020, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1021, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1022, "; Multiplicity;Counts ", 5,0,5);

  HB1( 1030, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1031, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1032, "; Multiplicity;Counts ", 5,0,5);
  */
  HB1(10000,"KF#{Lambda} pvalue",100,0,1);
  HB1(10001,"KF#{Lambda} chisqr",1000,0,15);
  HB1(10002,"KF#{Lambda} mass",1000,pdg::LambdaMass()-0.2,pdg::LambdaMass()+0.2);
  HB1(10003,"KF#{Lambda} dist",1000,0,20);
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
  HB1(10021,"KF#{Lambda} residual p_{#theta}",1000,-0.1,0.1);
  HB1(10022,"KF#{Lambda} residual p_{#phi}",1000,-0.1,0.1);
  HB1(10023,"KF#{Lambda} residual #pi_{p}",1000,-0.3,0.3);
  HB1(10024,"KF#{Lambda} residual #pi_{#theta}",1000,-0.1,0.1);
  HB1(10025,"KF#{Lambda} residual #pi_{#phi}",1000,-0.1,0.1);

  HB1(20000,"KF#{Xi}^{-} pvalue",100,0,1);
  HB1(20001,"KF#{Xi}^{-} chisqr",1000,0,15);
  HB1(20002,"KF#{Xi} mass",1000,pdg::XiMinusMass()-0.2,pdg::XiMinusMass()+0.2);
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
#endif
  HBTree( "tpc", "tree of GenfitCarbon" );
  tree->Branch( "status", &event.status );
  tree->Branch( "evnum", &event.evnum );
  tree->Branch( "trigpat", &event.trigpat );
  tree->Branch( "trigflag", &event.trigflag );

  tree->Branch( "nhHtof", &event.nhHtof );
  tree->Branch( "HtofSeg", &event.HtofSeg );
  tree->Branch( "tHtof", &event.tHtof );
  tree->Branch( "dtHtof", &event.dtHtof );
  tree->Branch( "deHtof", &event.deHtof );
  tree->Branch( "posHtof", &event.posHtof );
  tree->Branch( "G4tidHtof", &event.G4tidHtof );

  tree->Branch( "nclTpc", &event.nclTpc );
  tree->Branch( "remain_nclTpc", &event.nclTpc );
#if SaveRawData
  tree->Branch( "remain_cluster_x", &event.remain_cluster_x );
  tree->Branch( "remain_cluster_y", &event.remain_cluster_y );
  tree->Branch( "remain_cluster_z", &event.remain_cluster_z );
  tree->Branch( "remain_cluster_de", &event.remain_cluster_de );
//  tree->Branch( "remain_cluster_size", &event.remain_cluster_size );
  tree->Branch( "remain_cluster_layer", &event.remain_cluster_layer );
#if 0
  tree->Branch( "remain_cluster_row_center", &event.remain_cluster_row_center );
  tree->Branch( "remain_cluster_mrow", &event.remain_cluster_mrow );
  tree->Branch( "remain_cluster_de_center", &event.remain_cluster_de_center );
  tree->Branch( "remain_cluster_x_center", &event.remain_cluster_x_center );
  tree->Branch( "remain_cluster_y_center", &event.remain_cluster_y_center );
  tree->Branch( "remain_cluster_z_center", &event.remain_cluster_z_center );
#endif
  tree->Branch( "remain_cluster_houghflag", &event.remain_cluster_houghflag );
#endif

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "isInTarget", &event.isInTarget );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "trackid", &event.trackid );
  tree->Branch( "isXi", &event.isXi );
  tree->Branch( "isBeam", &event.isBeam );
  tree->Branch( "isKurama", &event.isKurama );
  tree->Branch( "isK18", &event.isK18 );
  tree->Branch( "isAccidental", &event.isAccidental );
  tree->Branch( "isMultiloop", &event.isMultiloop );
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
#if SaveRawData
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
  tree->Branch( "resolution_x", &event.resolution_x);
  tree->Branch( "resolution_y", &event.resolution_y);
  tree->Branch( "resolution_z", &event.resolution_z);
  tree->Branch( "pathhit", &event.pathhit);
  tree->Branch( "alpha", &event.alpha);
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
//  tree->Branch( "track_cluster_size", &event.track_cluster_size);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
#if 0
  tree->Branch( "track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch( "track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch( "track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch( "track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch( "track_cluster_row_center", &event.track_cluster_row_center);
#endif
#endif
  tree->Branch( "helix_t", &event.helix_t );
  tree->Branch( "chargeIndistinguishable", &event.chargeIndistinguishable );
  tree->Branch( "chisqr_inverted", &event.chisqr_inverted );
  tree->Branch( "pval_inverted", &event.pval_inverted );
  tree->Branch( "helix_cx_inverted", &event.helix_cx_inverted );
  tree->Branch( "helix_cy_inverted", &event.helix_cy_inverted );
  tree->Branch( "helix_z0_inverted", &event.helix_z0_inverted );
  tree->Branch( "helix_r_inverted", &event.helix_r_inverted );
  tree->Branch( "helix_dz_inverted", &event.helix_dz_inverted );
  tree->Branch( "mom0_inverted", &event.mom0_inverted );
  tree->Branch( "pid_inverted", &event.pid_inverted );

  tree->Branch( "nvtxTpc", &event.nvtxTpc );
  tree->Branch( "vtx_x", &event.vtx_x );
  tree->Branch( "vtx_y", &event.vtx_y );
  tree->Branch( "vtx_z", &event.vtx_z );
  tree->Branch( "vtx_dist", &event.vtx_dist );
  tree->Branch( "vtx_angle", &event.vtx_angle );
  tree->Branch( "vtxid", &event.vtxid );
  tree->Branch( "vtxmom_theta", &event.vtxmom_theta );
  tree->Branch( "vtxpos_x", &event.vtxpos_x );
  tree->Branch( "vtxpos_y", &event.vtxpos_y );
  tree->Branch( "vtxpos_z", &event.vtxpos_z );
  tree->Branch( "vtxmom_x", &event.vtxmom_x );
  tree->Branch( "vtxmom_y", &event.vtxmom_y );
  tree->Branch( "vtxmom_z", &event.vtxmom_z );

  tree->Branch( "isLambda", &event.isLambda );
  tree->Branch( "ncombiLambda", &event.ncombiLambda );
  tree->Branch( "distLambda", &event.distLambda );
  tree->Branch( "angleLambda", &event.angleLambda );
  tree->Branch( "bestmassLambda", &event.bestmassLambda );
  tree->Branch( "massLambda", &event.massLambda );
  tree->Branch( "vtxLambda_x", &event.vtxLambda_x );
  tree->Branch( "vtxLambda_y", &event.vtxLambda_y );
  tree->Branch( "vtxLambda_z", &event.vtxLambda_z );
  tree->Branch( "momLambda", &event.momLambda );
  tree->Branch( "momLambda_x", &event.momLambda_x );
  tree->Branch( "momLambda_y", &event.momLambda_y );
  tree->Branch( "momLambda_z", &event.momLambda_z );
  tree->Branch( "decaysidLambda", &event.decaysidLambda );
  tree->Branch( "decaysmomLambda", &event.decaysmomLambda );
  tree->Branch( "decaysmomLambda_x", &event.decaysmomLambda_x );
  tree->Branch( "decaysmomLambda_y", &event.decaysmomLambda_y );
  tree->Branch( "decaysmomLambda_z", &event.decaysmomLambda_z );

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

  tree->Branch("G4l1id",&event.G4l1id);
  tree->Branch("G4l1vtx_x",&event.G4l1vtx_x);
  tree->Branch("G4l1vtx_y",&event.G4l1vtx_y);
  tree->Branch("G4l1vtx_z",&event.G4l1vtx_z);
  tree->Branch("G4l1mom",&event.G4l1mom);
  tree->Branch("G4l1mom_x",&event.G4l1mom_x);
  tree->Branch("G4l1mom_y",&event.G4l1mom_y);
  tree->Branch("G4l1mom_z",&event.G4l1mom_z);

  tree->Branch("G4l2id",&event.G4l2id);
  tree->Branch("G4l2vtx_x",&event.G4l2vtx_x);
  tree->Branch("G4l2vtx_y",&event.G4l2vtx_y);
  tree->Branch("G4l2vtx_z",&event.G4l2vtx_z);
  tree->Branch("G4l2mom",&event.G4l2mom);
  tree->Branch("G4l2mom_x",&event.G4l2mom_x);
  tree->Branch("G4l2mom_y",&event.G4l2mom_y);
  tree->Branch("G4l2mom_z",&event.G4l2mom_z);

  tree->Branch("G4llmass",&event.G4llmass);

  tree->Branch("G4p1id",&event.G4p1id);
  tree->Branch("G4p1tid",&event.G4p1tid);
  tree->Branch("G4p1nh",&event.G4p1nh);
  tree->Branch("G4p1tnh",&event.G4p1tnh);
  tree->Branch("p1tid",&event.p1tid);
  tree->Branch("G4p1vtx_x",&event.G4p1vtx_x);
  tree->Branch("G4p1vtx_y",&event.G4p1vtx_y);
  tree->Branch("G4p1vtx_z",&event.G4p1vtx_z);
  tree->Branch("G4p1mom",&event.G4p1mom);
  tree->Branch("G4p1mom_x",&event.G4p1mom_x);
  tree->Branch("G4p1mom_y",&event.G4p1mom_y);
  tree->Branch("G4p1mom_z",&event.G4p1mom_z);

  tree->Branch("G4p2id",&event.G4p2id);
  tree->Branch("G4p2tid",&event.G4p2tid);
  tree->Branch("G4p2nh",&event.G4p2nh);
  tree->Branch("G4p2tnh",&event.G4p2tnh);
  tree->Branch("p2tid",&event.p2tid);
  tree->Branch("G4p2vtx_x",&event.G4p2vtx_x);
  tree->Branch("G4p2vtx_y",&event.G4p2vtx_y);
  tree->Branch("G4p2vtx_z",&event.G4p2vtx_z);
  tree->Branch("G4p2mom",&event.G4p2mom);
  tree->Branch("G4p2mom_x",&event.G4p2mom_x);
  tree->Branch("G4p2mom_y",&event.G4p2mom_y);
  tree->Branch("G4p2mom_z",&event.G4p2mom_z);

  tree->Branch("l1good", &event.l1good);
  tree->Branch("l2good", &event.l2good);
  tree->Branch("llswapped", &event.llswap);

  tree->Branch("p1_tracked", &event.p1_tracked);
  tree->Branch("p2_tracked", &event.p2_tracked);
  tree->Branch("piswapped", &event.piswap);
  tree->Branch("p1t_mom0", &event.p1t_mom0);
  tree->Branch("p2t_mom0", &event.p2t_mom0);

  tree->Branch("G4pi1id",&event.G4pi1id);
  tree->Branch("G4pi1tid",&event.G4pi1tid);
  tree->Branch("G4pi1nh",&event.G4pi1nh);
  tree->Branch("G4pi1tnh",&event.G4pi1tnh);
  tree->Branch("pi1tid",&event.pi1tid);
  tree->Branch("G4pi1vtx_x",&event.G4pi1vtx_x);
  tree->Branch("G4pi1vtx_y",&event.G4pi1vtx_y);
  tree->Branch("G4pi1vtx_z",&event.G4pi1vtx_z);
  tree->Branch("G4pi1mom",&event.G4pi1mom);
  tree->Branch("G4pi1mom_x",&event.G4pi1mom_x);
  tree->Branch("G4pi1mom_y",&event.G4pi1mom_y);
  tree->Branch("G4pi1mom_z",&event.G4pi1mom_z);

  tree->Branch("G4pi2id",&event.G4pi2id);
  tree->Branch("G4pi2tid",&event.G4pi2tid);
  tree->Branch("G4pi2nh",&event.G4pi2nh);
  tree->Branch("G4pi2tnh",&event.G4pi2tnh);
  tree->Branch("pi2tid",&event.pi2tid);
  tree->Branch("G4pi2vtx_x",&event.G4pi2vtx_x);
  tree->Branch("G4pi2vtx_y",&event.G4pi2vtx_y);
  tree->Branch("G4pi2vtx_z",&event.G4pi2vtx_z);
  tree->Branch("G4pi2mom",&event.G4pi2mom);
  tree->Branch("G4pi2mom_x",&event.G4pi2mom_x);
  tree->Branch("G4pi2mom_y",&event.G4pi2mom_y);
  tree->Branch("G4pi2mom_z",&event.G4pi2mom_z);

  tree->Branch("pi1_tracked", &event.pi1_tracked);
  tree->Branch("pi2_tracked", &event.pi2_tracked);
  tree->Branch("pi1t_mom0", &event.pi1t_mom0);
  tree->Branch("pi2t_mom0", &event.pi2t_mom0);

  tree->Branch( "nvtxTpcClustered", &event.nvtxTpcClustered );
  tree->Branch( "clusteredVtx_x", &event.clusteredVtx_x );
  tree->Branch( "clusteredVtx_y", &event.clusteredVtx_y );
  tree->Branch( "clusteredVtx_z", &event.clusteredVtx_z );
  tree->Branch( "clusteredVtxid", &event.clusteredVtxid );

  tree->Branch( "ncombiReconFaildLambda", &event.ncombiLreconfailed  );
  tree->Branch( "ReconFailedLambdaPId", &event.pidLreconfailed);
  tree->Branch( "ReconFailedLambdaPiId", &event.piidLreconfailed);
  tree->Branch( "ReconFailedLmassdMass", &event.LmassLreconfailed);
  tree->Branch( "ReconFailedLambdaDecayVtx_x", &event.LdecayvtxLreconfailed_x);
  tree->Branch( "ReconFailedLambdaDecayVtx_y", &event.LdecayvtxLreconfailed_y);
  tree->Branch( "ReconFailedLambdaDecayVtx_z", &event.LdecayvtxLreconfailed_z);
  tree->Branch( "ReconFailedLambdaMom", &event.LmomLreconfailed );
  tree->Branch( "ReconFailedLambdaMom_x", &event.LmomLreconfailed_x );
  tree->Branch( "ReconFailedLambdaMom_y", &event.LmomLreconfailed_y );
  tree->Branch( "ReconFailedLambdaMom_z", &event.LmomLreconfailed_z );
  tree->Branch( "ReconFailedLambdaPMom", &event.pmomLreconfailed );
  tree->Branch( "ReconFailedLambdaPMom_x", &event.pmomLreconfailed_x );
  tree->Branch( "ReconFailedLambdaPMom_y", &event.pmomLreconfailed_y );
  tree->Branch( "ReconFailedLambdaPMom_z", &event.pmomLreconfailed_z );
  tree->Branch( "ReconFailedLambdaPiMom", &event.pimomLreconfailed );
  tree->Branch( "ReconFailedLambdaPiMom_x", &event.pimomLreconfailed_x );
  tree->Branch( "ReconFailedLambdaPiMom_y", &event.pimomLreconfailed_y );
  tree->Branch( "ReconFailedLambdaPiMom_z", &event.pimomLreconfailed_z );
  tree->Branch( "ReconFailedLambdavtxCloseDist", & event.ppidistLreconfailed);

  tree->Branch( "ncombiPiPair", &event.ncombiPipair);
  tree->Branch( "PiPairPipId", &event.pipidPipair);
  tree->Branch( "PiPairPimId", &event.pimidPipair);
  tree->Branch( "PiPairPipMom", &event.pipmomPipair);
  tree->Branch( "PiPairPipMom_x", &event.pipmomPipair_x);
  tree->Branch( "PiPairPipMom_y", &event.pipmomPipair_y);
  tree->Branch( "PiPairPipMom_z", &event.pipmomPipair_z);
  tree->Branch( "PiPairPimMom", &event.pimmomPipair);
  tree->Branch( "PiPairPimMom_x", &event.pimmomPipair_x);
  tree->Branch( "PiPairPimMom_y", &event.pimmomPipair_y);
  tree->Branch( "PiPairPimMom_z", &event.pimmomPipair_z);
  tree->Branch( "PiPairMom", & event.momPipair);
  tree->Branch( "PiPairMom_x", & event.momPipair_x);
  tree->Branch( "PiPairMom_y", & event.momPipair_y);
  tree->Branch( "PiPairMom_z", & event.momPipair_z);
  tree->Branch( "PiPairReconLambdaMass", & event.reconLmassPipair);
  tree->Branch( "PiPairCloseDist", & event.pipidistPipair);

  tree->Branch( "ntK18", &event.ntK18);
  tree->Branch( "pK18", &event.pK18);
  tree->Branch( "chisqrK18", &event.chisqrK18);
  tree->Branch( "xtgtK18", &event.xtgtK18);
  tree->Branch( "ytgtK18", &event.ytgtK18);
  tree->Branch( "utgtK18", &event.utgtK18);
  tree->Branch( "vtgtK18", &event.vtgtK18);

  tree->Branch( "isgoodTPCK18", &event.isgoodTPCK18);
  tree->Branch( "chisqrTPCK18", &event.chisqrTPCK18);
  tree->Branch( "qTPCK18", &event.qTPCK18);
  tree->Branch( "pTPCK18", &event.pTPCK18);
  tree->Branch( "xtgtTPCK18", &event.xtgtTPCK18);
  tree->Branch( "ytgtTPCK18", &event.ytgtTPCK18);
  tree->Branch( "utgtTPCK18", &event.utgtTPCK18);
  tree->Branch( "vtgtTPCK18", &event.vtgtTPCK18);
  tree->Branch( "thetaTPCK18", &event.thetaTPCK18);
  tree->Branch( "lhtofTPCK18", &event.lhtofTPCK18);
  tree->Branch( "xhtofTPCK18", &event.xhtofTPCK18);
  tree->Branch( "yhtofTPCK18", &event.yhtofTPCK18);
  tree->Branch( "lvpTPCK18", &event.lvpTPCK18);
  tree->Branch( "xvpTPCK18", &event.xvpTPCK18);
  tree->Branch( "yvpTPCK18", &event.yvpTPCK18);

  tree->Branch( "ntKurama", &event.ntKurama);
  tree->Branch( "pKurama", &event.pKurama);
  tree->Branch( "qKurama", &event.qKurama);
  tree->Branch( "chisqrKurama", &event.chisqrKurama);
  tree->Branch( "m2Kurama", &event.m2);
  tree->Branch( "xtgtKurama", &event.xtgtKurama);
  tree->Branch( "ytgtKurama", &event.ytgtKurama);
  tree->Branch( "utgtKurama", &event.utgtKurama);
  tree->Branch( "vtgtKurama", &event.vtgtKurama);
  tree->Branch( "thetaKurama", &event.thetaKurama);
  tree->Branch( "pathwcKurama", &event.pathwcKurama);
  tree->Branch( "xin", &event.xin);
  tree->Branch( "yin", &event.yin);
  tree->Branch( "zin", &event.zin);
  tree->Branch( "pxin", &event.pxin);
  tree->Branch( "pyin", &event.pyin);
  tree->Branch( "pzin", &event.pzin);
  tree->Branch( "isgoodTPCKurama", &event.isgoodTPCKurama);
  tree->Branch( "kflagTPCKurama", &event.kflagTPCKurama);
  tree->Branch( "chisqrTPCKurama", &event.chisqrTPCKurama);
  tree->Branch( "pTPCKurama", &event.pTPCKurama);
  tree->Branch( "qTPCKurama", &event.qTPCKurama);
  tree->Branch( "m2TPCKurama", &event.m2TPCKurama);
  tree->Branch( "xtgtTPCKurama", &event.xtgtTPCKurama);
  tree->Branch( "ytgtTPCKurama", &event.ytgtTPCKurama);
  tree->Branch( "utgtTPCKurama", &event.utgtTPCKurama);
  tree->Branch( "vtgtTPCKurama", &event.vtgtTPCKurama);
  tree->Branch( "thetaTPCKurama", &event.thetaTPCKurama);
  tree->Branch( "pathTPCKurama", &event.pathTPCKurama);
  tree->Branch( "lhtofTPCKurama", &event.lhtofTPCKurama);
  tree->Branch( "xhtofTPCKurama", &event.xhtofTPCKurama);
  tree->Branch( "yhtofTPCKurama", &event.yhtofTPCKurama);
  tree->Branch( "lvpTPCKurama", &event.lvpTPCKurama);
  tree->Branch( "xvpTPCKurama", &event.xvpTPCKurama);
  tree->Branch( "yvpTPCKurama", &event.yvpTPCKurama);

  tree->Branch( "isgoodTPC", &event.isgoodTPC);
  tree->Branch( "insideTPC", &event.insideTPC);
  tree->Branch( "vtxTPC", &event.vtxTPC);
  tree->Branch( "vtyTPC", &event.vtyTPC);
  tree->Branch( "vtzTPC", &event.vtzTPC);
  tree->Branch( "closeDistTPC", &event.closeDistTPC);
  tree->Branch( "MissMassTPC", &event.MissMassTPC);
  tree->Branch( "MissMassCorrTPC", &event.MissMassCorrTPC);
  tree->Branch( "MissMassCorrDETPC", &event.MissMassCorrDETPC);
  tree->Branch( "pOrgTPC", &event.pOrgTPC);
  tree->Branch( "pCalcTPC", &event.pCalcTPC);
  tree->Branch( "pCorrTPC", &event.pCorrTPC);
  tree->Branch( "pCorrDETPC", &event.pCorrDETPC);
  tree->Branch("xbTPC", &event.xbTPC);
  tree->Branch("ybTPC", &event.ybTPC);
  tree->Branch("ubTPC", &event.ubTPC);
  tree->Branch("vbTPC", &event.vbTPC);
  tree->Branch("xsTPC", &event.xsTPC);
  tree->Branch("ysTPC", &event.ysTPC);
  tree->Branch("usTPC", &event.usTPC);
  tree->Branch("vsTPC", &event.vsTPC);

  tree->Branch("nKK", &event.nKK);
  tree->Branch("Kflag", &event.Kflag);
  tree->Branch("MissMass", &event.MissMass);
  tree->Branch("MissMassCorr", &event.MissMassCorr);
  tree->Branch("MissMassCorrDE", &event.MissMassCorrDE);
  tree->Branch("vtx", &event.vtx);
  tree->Branch("vty", &event.vty);
  tree->Branch("vtz", &event.vtz);
  tree->Branch("pOrg", &event.pOrg);
  tree->Branch("pCalc", &event.pCalc);
  tree->Branch("pCorr", &event.pCorr);
  tree->Branch("pCorrDE", &event.pCorrDE);
  tree->Branch("xb", &event.xb);
  tree->Branch("yb", &event.yb);
  tree->Branch("ub", &event.ub);
  tree->Branch("vb", &event.vb);
  tree->Branch("xs", &event.xs);
  tree->Branch("ys", &event.ys);
  tree->Branch("us", &event.us);
  tree->Branch("vs", &event.vs);

  tree->Branch("KmMom_x", &event.km_mom_x);
  tree->Branch("KmMom_y", &event.km_mom_y);
  tree->Branch("KmMom_z", &event.km_mom_z);
  tree->Branch("KpMom_x", &event.kp_mom_x);
  tree->Branch("KpMom_y", &event.kp_mom_y);
  tree->Branch("KpMom_z", &event.kp_mom_z);

  tree->Branch("BE", &event.BE);
  tree->Branch("BETPC", &event.BETPC);
  tree->Branch("BE_LL", &event.BE_LL);
  tree->Branch("BETPC_LL", &event.BETPC_LL);

  //Multi-track vertex
  tree->Branch("GFKKXiProductionVtx_x", &event.GFprodvtx_x_kkxi);
  tree->Branch("GFKKXiProductionVtx_y", &event.GFprodvtx_y_kkxi);
  tree->Branch("GFKKXiProductionVtx_z", &event.GFprodvtx_z_kkxi);
  tree->Branch("GFKKLLProductionVtx_x", &event.GFprodvtx_x_ll);
  tree->Branch("GFKKLLProductionVtx_y", &event.GFprodvtx_y_ll);
  tree->Branch("GFKKLLProductionVtx_z", &event.GFprodvtx_z_ll);
  tree->Branch("GFKKLProductionVtx_x1", &event.GFprodvtx_x_l1);
  tree->Branch("GFKKLProductionVtx_y1", &event.GFprodvtx_y_l1);
  tree->Branch("GFKKLProductionVtx_z1", &event.GFprodvtx_z_l1);
  tree->Branch("GFKKLProductionVtx_x2", &event.GFprodvtx_x_l2);
  tree->Branch("GFKKLProductionVtx_y2", &event.GFprodvtx_y_l2);
  tree->Branch("GFKKLProductionVtx_z2", &event.GFprodvtx_z_l2);
  tree->Branch("GFKKLProductionVtx_x", &event.GFprodvtx_x_l);
  tree->Branch("GFKKLProductionVtx_y", &event.GFprodvtx_y_l);
  tree->Branch("GFKKLProductionVtx_z", &event.GFprodvtx_z_l);

  //track fitting results
  //for all tracks
  //tree->Branch("GFntTpc", &event.GFntTpc);
  //for decay particles
  tree->Branch("GFntDecays", &event.GFntdecays);
  tree->Branch("GFnhtrack", &event.GFnhtrack);
  tree->Branch("GFchisqr", &event.GFchisqr);
  tree->Branch("GFcharge", &event.GFcharge);
  tree->Branch("GFtof", &event.GFtof);
  tree->Branch("GFtracklen", &event.GFtracklen);
  tree->Branch("GFpval", &event.GFpval);
  tree->Branch("GFpdgcode", &event.GFpdgcode);

  tree->Branch("GFDecaysHtofId", &event.GFdecays_htofid);
  tree->Branch("GFDecaysTrackLen", &event.GFdecays_tracklen);
  tree->Branch("GFDecaysTof", &event.GFdecays_tof);
  tree->Branch("GFDecaysMassSquare", &event.GFdecays_mass2);
  tree->Branch("GFDecaysMom", &event.GFdecays_mom);
  tree->Branch("GFDecaysMom_x", &event.GFdecays_mom_x);
  tree->Branch("GFDecaysMom_y", &event.GFdecays_mom_y);
  tree->Branch("GFDecaysMom_z", &event.GFdecays_mom_z);
  tree->Branch("GFDecaysMomCM", &event.GFdecays_CMmom);
  tree->Branch("GFDecaysMomCM_x", &event.GFdecays_CMmom_x);
  tree->Branch("GFDecaysMomCM_y", &event.GFdecays_CMmom_y);
  tree->Branch("GFDecaysMomCM_z", &event.GFdecays_CMmom_z);
  tree->Branch("GFMomLoss", &event.GFmomloss);
  tree->Branch("GFELoss", &event.GFeloss);

  tree->Branch("Xiflag", &event.xiflag);
  tree->Branch("XiMass", &event.ximass);
  tree->Branch("XiDecayVtx_x", &event.xidecayvtx_x);
  tree->Branch("XiDecayVtx_y", &event.xidecayvtx_y);
  tree->Branch("XiDecayVtx_z", &event.xidecayvtx_z);
  tree->Branch("XiMom", &event.ximom);
  tree->Branch("XiMom_x", &event.ximom_x);
  tree->Branch("XiMom_y", &event.ximom_y);
  tree->Branch("XiMom_z", &event.ximom_z);
  tree->Branch("XiVtxCloseDist", &event.lpi_dist);
  tree->Branch("LambdaMass", &event.lmass);
  tree->Branch("LambdaInTarget", &event.l_intarget);
  tree->Branch("LambdaDecayVtx_x", &event.ldecayvtx_x);
  tree->Branch("LambdaDecayVtx_y", &event.ldecayvtx_y);
  tree->Branch("LambdaDecayVtx_z", &event.ldecayvtx_z);
  tree->Branch("LambdaMom", &event.lmom);
  tree->Branch("LambdaMom_x", &event.lmom_x);
  tree->Branch("LambdaMom_y", &event.lmom_y);
  tree->Branch("LambdaMom_z", &event.lmom_z);
  tree->Branch("LambdaVtxCloseDist", &event.ppi_dist);
  tree->Branch("XiTarget_x", &event.xitargetvtx_x);
  tree->Branch("XiTarget_y", &event.xitargetvtx_y);
  tree->Branch("XiTarget_z", &event.xitargetvtx_z);
  tree->Branch("XiTargetMom", &event.xitargetmom);
  tree->Branch("XiTargetMom_x", &event.xitargetmom_x);
  tree->Branch("XiTargetMom_y", &event.xitargetmom_y);
  tree->Branch("XiTargetMom_z", &event.xitargetmom_z);
  tree->Branch("XiTargetCloseDist", &event.xitarget_dist);

  tree->Branch("GFXiMass", &event.GFximass);
  tree->Branch("GFXiDecayVtx_x", &event.GFxidecayvtx_x);
  tree->Branch("GFXiDecayVtx_y", &event.GFxidecayvtx_y);
  tree->Branch("GFXiDecayVtx_z", &event.GFxidecayvtx_z);
  tree->Branch("GFXiMom", &event.GFximom);
  tree->Branch("GFXiMom_x", &event.GFximom_x);
  tree->Branch("GFXiMom_y", &event.GFximom_y);
  tree->Branch("GFXiMom_z", &event.GFximom_z);
  tree->Branch("GFXiVtxCloseDist", &event.GFlpi_dist);

  tree->Branch("GFXiMass_RefitProton", &event.GFximass_pmomcorr);
  tree->Branch("GFXiMom_RefitProton", &event.GFximom_pmomcorr);
  tree->Branch("GFLambdaMass_RefitProton", &event.GFlmass_pmomcorr);
  tree->Branch("GFLambdaMom_RefitProton", &event.GFlmom_pmomcorr);
  tree->Branch("GFDecaysMom_RefitProton", &event.GFdecays_mom_pmomcorr);
  tree->Branch("GFDecaysMom_x_RefitProton", &event.GFdecays_mom_x_pmomcorr);
  tree->Branch("GFDecaysMom_y_RefitProton", &event.GFdecays_mom_y_pmomcorr);
  tree->Branch("GFDecaysMom_z_RefitProton", &event.GFdecays_mom_z_pmomcorr);
  tree->Branch("GFpval_RefitProton", &event.GFpval_pmomcorr);
  tree->Branch("GFchisqr_RefitProton", &event.GFchisqr_pmomcorr);

  //extrapolation
  //to the K, K vertex
  tree->Branch("GFXiKKVtx_x", &event.GFxikkvtx_x);
  tree->Branch("GFXiKKVtx_y", &event.GFxikkvtx_y);
  tree->Branch("GFXiKKVtx_z", &event.GFxikkvtx_z);
  tree->Branch("GFXiKKVtxMom", &event.GFxikkmom);
  tree->Branch("GFXiKKVtxMom_x", &event.GFxikkmom_x);
  tree->Branch("GFXiKKVtxMom_y", &event.GFxikkmom_y);
  tree->Branch("GFXiKKVtxMom_z", &event.GFxikkmom_z);
  tree->Branch("GFXiKKVtxCloseDist", &event.GFxikkvtx_dist);
  //to the production vertex
  tree->Branch("GFXiProductionVtx_x", &event.GFxiprodvtx_x);
  tree->Branch("GFXiProductionVtx_y", &event.GFxiprodvtx_y);
  tree->Branch("GFXiProductionVtx_z", &event.GFxiprodvtx_z);
  tree->Branch("GFXiProductionVtxMom", &event.GFxiprodmom);
  tree->Branch("GFXiProductionVtxMom_x", &event.GFxiprodmom_x);
  tree->Branch("GFXiProductionVtxMom_y", &event.GFxiprodmom_y);
  tree->Branch("GFXiProductionVtxMom_z", &event.GFxiprodmom_z);
  tree->Branch("GFXiProductionVtxCloseDist", &event.GFxiprodvtx_dist);
  tree->Branch("GFXiTrackLen", &event.GFxitracklen);
  tree->Branch("GFXiTof", &event.GFxitof);
  tree->Branch("GFXiMomLoss", &event.GFximomloss);
  tree->Branch("GFXiExcitation", &event.GFxiexcitation);
  //closest point to the target
  tree->Branch("GFXiTarget_x", &event.GFxitargetvtx_x);
  tree->Branch("GFXiTarget_y", &event.GFxitargetvtx_y);
  tree->Branch("GFXiTarget_z", &event.GFxitargetvtx_z);
  tree->Branch("GFXiTargetMom", &event.GFxitargetmom);
  tree->Branch("GFXiTargetMom_x", &event.GFxitargetmom_x);
  tree->Branch("GFXiTargetMom_y", &event.GFxitargetmom_y);
  tree->Branch("GFXiTargetMom_z", &event.GFxitargetmom_z);
  tree->Branch("GFXiTargetCloseDist", &event.GFxitarget_dist);
  //at z=z_target
  tree->Branch("GFXiTargetCenter_x", &event.GFxitargetcenter_x);
  tree->Branch("GFXiTargetCenter_y", &event.GFxitargetcenter_y);
  tree->Branch("GFXiTargetCenter_z", &event.GFxitargetcenter_z);
  tree->Branch("GFXiTargetCenterMom", &event.GFxitargetcentermom);
  tree->Branch("GFXiTargetCenterMom_x", &event.GFxitargetcentermom_x);
  tree->Branch("GFXiTargetCenterMom_y", &event.GFxitargetcentermom_y);
  tree->Branch("GFXiTargetCenterMom_z", &event.GFxitargetcentermom_z);
  tree->Branch("GFXiTargetCenterCloseDist", &event.GFxitargetcenter_dist);

  tree->Branch("LLflag", &event.llflag);
  tree->Branch("LambdaLambdaVtx_x", &event.llvtx_x);
  tree->Branch("LambdaLambdaVtx_y", &event.llvtx_y);
  tree->Branch("LambdaLambdaVtx_z", &event.llvtx_z);
  tree->Branch("LambdaLambdaCloseDist", &event.lldist);

  tree->Branch("LambdaMass1", &event.lmass1);
  tree->Branch("LambdaDecayVtx_x1", &event.ldecayvtx_x1);
  tree->Branch("LambdaDecayVtx_y1", &event.ldecayvtx_y1);
  tree->Branch("LambdaDecayVtx_z1", &event.ldecayvtx_z1);
  tree->Branch("LambdaMom1", &event.lmom1);
  tree->Branch("LambdaMom_x1", &event.lmom_x1);
  tree->Branch("LambdaMom_y1", &event.lmom_y1);
  tree->Branch("LambdaMom_z1", &event.lmom_z1);
  tree->Branch("LambdaVtxCloseDist1", &event.ppi_dist1);
  tree->Branch("LambdaTargetCloseDist1", &event.ltarget_dist1);
  tree->Branch("LambdaTargetCloseVtx_x1", &event.ltargetvtx_x1);
  tree->Branch("LambdaTargetCloseVtx_y1", &event.ltargetvtx_y1);
  tree->Branch("LambdaTargetCloseVtx_z1", &event.ltargetvtx_z1);
  tree->Branch("KFLambdaChisqr1", &event.KFlchisqr1);

  tree->Branch("LambdaMass2", &event.lmass2);
  tree->Branch("LambdaDecayVtx_x2", &event.ldecayvtx_x2);
  tree->Branch("LambdaDecayVtx_y2", &event.ldecayvtx_y2);
  tree->Branch("LambdaDecayVtx_z2", &event.ldecayvtx_z2);
  tree->Branch("LambdaMom2", &event.lmom2);
  tree->Branch("LambdaMom_x2", &event.lmom_x2);
  tree->Branch("LambdaMom_y2", &event.lmom_y2);
  tree->Branch("LambdaMom_z2", &event.lmom_z2);
  tree->Branch("LambdaVtxCloseDist2", &event.ppi_dist2);
  tree->Branch("LambdaTargetCloseDist2", &event.ltarget_dist2);
  tree->Branch("LambdaTargetCloseVtx_x2", &event.ltargetvtx_x2);
  tree->Branch("LambdaTargetCloseVtx_y2", &event.ltargetvtx_y2);
  tree->Branch("LambdaTargetCloseVtx_z2", &event.ltargetvtx_z2);
  tree->Branch("KFLambdaChisqr2", &event.KFlchisqr2);

  tree->Branch("GFLLExcitation", &event.GFllexcitation);
  tree->Branch("GFLambdaMass1", &event.GFlmass1);
  tree->Branch("GFLambdaDecayVtx_x1", &event.GFldecayvtx_x1);
  tree->Branch("GFLambdaDecayVtx_y1", &event.GFldecayvtx_y1);
  tree->Branch("GFLambdaDecayVtx_z1", &event.GFldecayvtx_z1);
  tree->Branch("GFLambdaMom1", &event.GFlmom1);
  tree->Branch("GFLambdaMom_x1", &event.GFlmom_x1);
  tree->Branch("GFLambdaMom_y1", &event.GFlmom_y1);
  tree->Branch("GFLambdaMom_z1", &event.GFlmom_z1);
  tree->Branch("GFLambdaVtxCloseDist1", &event.GFppi_dist1);

  tree->Branch("GFLambdaMass2", &event.GFlmass2);
  tree->Branch("GFLambdaDecayVtx_x2", &event.GFldecayvtx_x2);
  tree->Branch("GFLambdaDecayVtx_y2", &event.GFldecayvtx_y2);
  tree->Branch("GFLambdaDecayVtx_z2", &event.GFldecayvtx_z2);
  tree->Branch("GFLambdaMom2", &event.GFlmom2);
  tree->Branch("GFLambdaMom_x2", &event.GFlmom_x2);
  tree->Branch("GFLambdaMom_y2", &event.GFlmom_y2);
  tree->Branch("GFLambdaMom_z2", &event.GFlmom_z2);
  tree->Branch("GFLambdaVtxCloseDist2", &event.GFppi_dist2);

  //L, L vertex
  tree->Branch("GFLambdaLambdaVtx_x", &event.GFllvtx_x);
  tree->Branch("GFLambdaLambdaVtx_y", &event.GFllvtx_y);
  tree->Branch("GFLambdaLambdaVtx_z", &event.GFllvtx_z);
  tree->Branch("GFLambdaLambdaCloseDist", &event.GFlldist);

  //extrapolation
  //to closest point to the target
  tree->Branch("GFLambdaTarget_x1", &event.GFltargetvtx_x1);
  tree->Branch("GFLambdaTarget_y1", &event.GFltargetvtx_y1);
  tree->Branch("GFLambdaTarget_z1", &event.GFltargetvtx_z1);
  tree->Branch("GFLambdaTargetCloseDist1", &event.GFltarget_dist1);
  tree->Branch("GFLambdaTarget_x2", &event.GFltargetvtx_x2);
  tree->Branch("GFLambdaTarget_y2", &event.GFltargetvtx_y2);
  tree->Branch("GFLambdaTarget_z2", &event.GFltargetvtx_z2);
  tree->Branch("GFLambdaTargetCloseDist2", &event.GFltarget_dist2);
  //at z=z_target
  tree->Branch("GFLambdaTargetCenter_x1", &event.GFltargetcenter_x1);
  tree->Branch("GFLambdaTargetCenter_y1", &event.GFltargetcenter_y1);
  tree->Branch("GFLambdaTargetCenter_z1", &event.GFltargetcenter_z1);
  tree->Branch("GFLambdaTargetCenterCloseDist1", &event.GFltargetcenter_dist1);
  tree->Branch("GFLambdaTargetCenter_x2", &event.GFltargetcenter_x2);
  tree->Branch("GFLambdaTargetCenter_y2", &event.GFltargetcenter_y2);
  tree->Branch("GFLambdaTargetCenter_z2", &event.GFltargetcenter_z2);
  tree->Branch("GFLambdaTargetCenterCloseDist2", &event.GFltargetcenter_dist2);
  //to the production vertex
  tree->Branch("GFLambdaProductionVtx_x1", &event.GFlprodvtx_x1);
  tree->Branch("GFLambdaProductionVtx_y1", &event.GFlprodvtx_y1);
  tree->Branch("GFLambdaProductionVtx_z1", &event.GFlprodvtx_z1);
  tree->Branch("GFLambdaProductionVtxCloseDist1", &event.GFlprodvtx_dist1);
  tree->Branch("GFLambdaTrackLen1", &event.GFltracklen1);
  tree->Branch("GFLambdaTof1", &event.GFltof1);
  tree->Branch("GFLambdaProductionVtx_x2", &event.GFlprodvtx_x2);
  tree->Branch("GFLambdaProductionVtx_y2", &event.GFlprodvtx_y2);
  tree->Branch("GFLambdaProductionVtx_z2", &event.GFlprodvtx_z2);
  tree->Branch("GFLambdaProductionVtxCloseDist2", &event.GFlprodvtx_dist2);
  tree->Branch("GFLambdaTrackLen2", &event.GFltracklen2);
  tree->Branch("GFLambdaTof2", &event.GFltof2);

  //Km, Kp, L1, L2 vertex
  tree->Branch("KFLLExcitation", &event.KFllexcitation);
  tree->Branch("KFKKLLProductionVtxChisqr", &event.KFprodvtx_chisqr_ll);
  tree->Branch("KFKKLLProductionVtx_x", &event.KFprodvtx_x_ll);
  tree->Branch("KFKKLLProductionVtx_y", &event.KFprodvtx_y_ll);
  tree->Branch("KFKKLLProductionVtx_z", &event.KFprodvtx_z_ll);
  //Km, Kp, L1 vertex
  tree->Branch("KFKKLProductionVtx_x1", &event.KFprodvtx_x_l1);
  tree->Branch("KFKKLProductionVtx_y1", &event.KFprodvtx_y_l1);
  tree->Branch("KFKKLProductionVtx_z1", &event.KFprodvtx_z_l1);
  //Km, Kp, L2 vertex
  tree->Branch("KFKKLProductionVtx_x2", &event.KFprodvtx_x_l2);
  tree->Branch("KFKKLProductionVtx_y2", &event.KFprodvtx_y_l2);
  tree->Branch("KFKKLProductionVtx_z2", &event.KFprodvtx_z_l2);
  //L, L vertex
  tree->Branch("KFLambdaLambdaVtx_x", &event.KFllvtx_x);
  tree->Branch("KFLambdaLambdaVtx_y", &event.KFllvtx_y);
  tree->Branch("KFLambdaLambdaVtx_z", &event.KFllvtx_z);
  tree->Branch("KFLambdaLambdaCloseDist", &event.KFlldist);
  //to the production vertex
  tree->Branch("KFLambdaProductionVtx_x1", &event.KFlprodvtx_x1);
  tree->Branch("KFLambdaProductionVtx_y1", &event.KFlprodvtx_y1);
  tree->Branch("KFLambdaProductionVtx_z1", &event.KFlprodvtx_z1);
  tree->Branch("KFLambdaProductionVtxCloseDist1", &event.KFlprodvtx_dist1);
  tree->Branch("KFLambdaTrackLen1", &event.KFltracklen1);
  tree->Branch("KFLambdaTof1", &event.KFltof1);
  tree->Branch("KFLambdaProductionVtx_x2", &event.KFlprodvtx_x2);
  tree->Branch("KFLambdaProductionVtx_y2", &event.KFlprodvtx_y2);
  tree->Branch("KFLambdaProductionVtx_z2", &event.KFlprodvtx_z2);
  tree->Branch("KFLambdaProductionVtxCloseDist2", &event.KFlprodvtx_dist2);
  tree->Branch("KFLambdaTrackLen2", &event.KFltracklen2);
  tree->Branch("KFLambdaTof2", &event.KFltof2);
  //L, L mom
  tree->Branch("KFLambdaMom1", &event.KFlmom1);
  tree->Branch("KFLambdaMom_x1", &event.KFlmom_x1);
  tree->Branch("KFLambdaMom_y1", &event.KFlmom_y1);
  tree->Branch("KFLambdaMom_z1", &event.KFlmom_z1);
  tree->Branch("KFLambdaMom2", &event.KFlmom2);
  tree->Branch("KFLambdaMom_x2", &event.KFlmom_x2);
  tree->Branch("KFLambdaMom_y2", &event.KFlmom_y2);
  tree->Branch("KFLambdaMom_z2", &event.KFlmom_z2);

  tree->Branch("Lflag", &event.lflag);
  tree->Branch("LambdaTargetCloseDist", &event.ltarget_dist);
  tree->Branch("LambdaTargetCloseVtx_x", &event.ltargetvtx_x);
  tree->Branch("LambdaTargetCloseVtx_y", &event.ltargetvtx_y);
  tree->Branch("LambdaTargetCloseVtx_z", &event.ltargetvtx_z);

  tree->Branch("GFLambdaMass", &event.GFlmass);
  tree->Branch("GFLambdaDecayVtx_x", &event.GFldecayvtx_x);
  tree->Branch("GFLambdaDecayVtx_y", &event.GFldecayvtx_y);
  tree->Branch("GFLambdaDecayVtx_z", &event.GFldecayvtx_z);
  tree->Branch("GFLambdaMom", &event.GFlmom);
  tree->Branch("GFLambdaMom_x", &event.GFlmom_x);
  tree->Branch("GFLambdaMom_y", &event.GFlmom_y);
  tree->Branch("GFLambdaMom_z", &event.GFlmom_z);
  tree->Branch("GFLambdaVtxCloseDist", &event.GFppi_dist);
  tree->Branch("GFLambdaTargetCloseDist", &event.GFltarget_dist);
  tree->Branch("GFLambdaTarget_x", &event.GFltargetvtx_x);
  tree->Branch("GFLambdaTarget_y", &event.GFltargetvtx_y);
  tree->Branch("GFLambdaTarget_z", &event.GFltargetvtx_z);
  tree->Branch("GFLambdaTargetCenter_x", &event.GFltargetcenter_x);
  tree->Branch("GFLambdaTargetCenter_y", &event.GFltargetcenter_y);
  tree->Branch("GFLambdaTargetCenter_z", &event.GFltargetcenter_z);
  tree->Branch("GFLambdaTargetCenterCloseDist", &event.GFltargetcenter_dist);
  tree->Branch("GFLambdaProductionVtx_x", &event.GFlprodvtx_x);
  tree->Branch("GFLambdaProductionVtx_y", &event.GFlprodvtx_y);
  tree->Branch("GFLambdaProductionVtx_z", &event.GFlprodvtx_z);
  tree->Branch("GFLambdaProductionVtxCloseDist", &event.GFlprodvtx_dist);
  tree->Branch("GFLambdaTrackLen", &event.GFltracklen);
  tree->Branch("GFLambdaTof", &event.GFltof);

  tree->Branch("KFKKLProductionVtx_x", &event.KFprodvtx_x_l);
  tree->Branch("KFKKLProductionVtx_y", &event.KFprodvtx_y_l);
  tree->Branch("KFKKLProductionVtx_z", &event.KFprodvtx_z_l);
  tree->Branch("KFLambdaProductionVtx_x", &event.KFlprodvtx_x);
  tree->Branch("KFLambdaProductionVtx_y", &event.KFlprodvtx_y);
  tree->Branch("KFLambdaProductionVtx_z", &event.KFlprodvtx_z);
  tree->Branch("KFLambdaProductionVtxCloseDist", &event.KFlprodvtx_dist);
  tree->Branch("KFLambdaTrackLen", &event.KFltracklen);
  tree->Branch("KFLambdaTof", &event.KFltof);

  tree->Branch("LPhiflag", &event.lphiflag);
  tree->Branch("PhiMass", &event.phimass);
  tree->Branch("PhiDecayVtx_x", &event.phidecayvtx_x);
  tree->Branch("PhiDecayVtx_y", &event.phidecayvtx_y);
  tree->Branch("PhiDecayVtx_z", &event.phidecayvtx_z);
  tree->Branch("PhiMom", &event.phimom);
  tree->Branch("PhiMom_x", &event.phimom_x);
  tree->Branch("PhiMom_y", &event.phimom_y);
  tree->Branch("PhiMom_z", &event.phimom_z);
  tree->Branch("PhiVtxCloseDist", &event.kk_dist);
  tree->Branch("PhiDecaysTrackId", &event.phidecays_id);
  tree->Branch("PhiDecaysMom", &event.phidecays_mom);
  tree->Branch("PhiDecaysMom_x", &event.phidecays_mom_x);
  tree->Branch("PhiDecaysMom_y", &event.phidecays_mom_y);
  tree->Branch("PhiDecaysMom_z", &event.phidecays_mom_z);

  tree->Branch("GFPhiMass", &event.GFphimass);
  tree->Branch("GFPhiDecayVtx_x", &event.GFphidecayvtx_x);
  tree->Branch("GFPhiDecayVtx_y", &event.GFphidecayvtx_y);
  tree->Branch("GFPhiDecayVtx_z", &event.GFphidecayvtx_z);
  tree->Branch("GFPhiMom", &event.GFphimom);
  tree->Branch("GFPhiMom_x", &event.GFphimom_x);
  tree->Branch("GFPhiMom_y", &event.GFphimom_y);
  tree->Branch("GFPhiMom_z", &event.GFphimom_z);
  tree->Branch("GFPhiVtxCloseDist", &event.GFkk_dist);
  tree->Branch("GFPhiProductionVtxCloseDist", &event.GFphiprodvtx_dist);
  tree->Branch("KmMass2", &event.phi_km_mass2);
  tree->Branch("GFPhiDecaysMom", &event.GFphidecays_mom);
  tree->Branch("GFPhiDecaysMom_x", &event.GFphidecays_mom_x);
  tree->Branch("GFPhiDecaysMom_y", &event.GFphidecays_mom_y);
  tree->Branch("GFPhiDecaysMom_z", &event.GFphidecays_mom_z);

  tree->Branch("LPiflag", &event.lpiflag);

  tree->Branch("DecaysTrackId", &event.decays_id);
  tree->Branch("DecaysMom", &event.decays_mom);
  tree->Branch("DecaysMom_x", &event.decays_mom_x);
  tree->Branch("DecaysMom_y", &event.decays_mom_y);
  tree->Branch("DecaysMom_z", &event.decays_mom_z);
  tree->Branch("DecaysMomCM", &event.decays_CMmom);
  tree->Branch("DecaysMomCM_x", &event.decays_CMmom_x);
  tree->Branch("DecaysMomCM_y", &event.decays_CMmom_y);
  tree->Branch("DecaysMomCM_z", &event.decays_CMmom_z);

  tree->Branch("PiPiflag", &event.pipiflag);
  tree->Branch("Pimflag", &event.pimflag);
  tree->Branch("Emptyflag", &event.emptyflag);

  //Remaining p, pi after Xi, L searching
  //Multiplicity means tracks comes from the target
  tree->Branch("ResidualsMultiplicity", &event.residual_multi);
  tree->Branch("pipMultiplicity", &event.pip_multi);
  tree->Branch("pMultiplicity", &event.p_multi);
  tree->Branch("pimMultiplicity", &event.pim_multi);
  tree->Branch("ppipMultiplicity", &event.ppip_multi);
  tree->Branch("AccidentalMultiplicity", &event.accident_multi);
  tree->Branch("ResidualsTrackId", &event.residual_id);
  tree->Branch("ResidualsMassSquare", &event.residual_mass2);
  tree->Branch("ResidualsCloseDistTgt", &event.residual_dist2tgt);
  tree->Branch("ResidualsMom", &event.residual_mom);
  tree->Branch("ResidualsMom_x", &event.residual_mom_x);
  tree->Branch("ResidualsMom_y", &event.residual_mom_y);
  tree->Branch("ResidualsMom_z", &event.residual_mom_z);
  tree->Branch("ResidualsCharge", &event.residual_charge);

  tree->Branch("KFLambdaMomPpi", &event.KFlmom0);
  tree->Branch("KFLambdaMomPpi_x", &event.KFlmom_x0);
  tree->Branch("KFLambdaMomPpi_y0", &event.KFlmom_y0);
  tree->Branch("KFLambdaMomPpi_z0", &event.KFlmom_z0);

  tree->Branch("KFLambdaMom", &event.KFlmom);
  tree->Branch("KFLambdaMom_x", &event.KFlmom_x);
  tree->Branch("KFLambdaMom_y", &event.KFlmom_y);
  tree->Branch("KFLambdaMom_z", &event.KFlmom_z);
  tree->Branch("KFLambdaChisqr", &event.KFlchisqr);
  tree->Branch("KFLambdaPval", &event.KFlpval);
  tree->Branch("KFLambdaPull",&event.KFlpull);

  tree->Branch("KFXiVtxCloseDist", &event.KFlpi_dist);
  tree->Branch("KFXiMom", &event.KFximom);
  tree->Branch("KFXiMom_x", &event.KFximom_x);
  tree->Branch("KFXiMom_y", &event.KFximom_y);
  tree->Branch("KFXiMom_z", &event.KFximom_z);
  tree->Branch("KFXiChisqr", &event.KFxichisqr);
  tree->Branch("KFXiPval", &event.KFxipval);
  tree->Branch("KFXiMass",&event.KFximass);
  tree->Branch("KFXiDecayVtx_x", &event.KFxidecayvtx_x);
  tree->Branch("KFXiDecayVtx_y", &event.KFxidecayvtx_y);
  tree->Branch("KFXiDecayVtx_z", &event.KFxidecayvtx_z);
  tree->Branch("KFXiPull",&event.KFxipull);

  //K, K, Xi vertex
  tree->Branch("KFKKXiProductionVtxChisqr", &event.KFprodvtx_chisqr_kkxi);
  tree->Branch("KFKKXiProductionVtx_x", &event.KFprodvtx_x_kkxi);
  tree->Branch("KFKKXiProductionVtx_y", &event.KFprodvtx_y_kkxi);
  tree->Branch("KFKKXiProductionVtx_z", &event.KFprodvtx_z_kkxi);
  //extrapolation to the production vertex
  tree->Branch("KFXiProductionVtx_x", &event.KFxiprodvtx_x);
  tree->Branch("KFXiProductionVtx_y", &event.KFxiprodvtx_y);
  tree->Branch("KFXiProductionVtx_z", &event.KFxiprodvtx_z);
  tree->Branch("KFXiProductionVtxMom", &event.KFxiprodmom);
  tree->Branch("KFXiProductionVtxMom_x", &event.KFxiprodmom_x);
  tree->Branch("KFXiProductionVtxMom_y", &event.KFxiprodmom_y);
  tree->Branch("KFXiProductionVtxMom_z", &event.KFxiprodmom_z);
  tree->Branch("KFXiProductionVtxCloseDist", &event.KFxiprodvtx_dist);
  tree->Branch("KFXiTrackLen", &event.KFxitracklen);
  tree->Branch("KFXiTof", &event.KFxitof);
  tree->Branch("KFXiMomLoss", &event.KFximomloss);
  tree->Branch("KFXiExcitation", &event.KFxiexcitation);

  //Kp, Xi vertex
  tree->Branch("KFKpXiProductionVtx_x", &event.KFprodvtx_x_kpxi);
  tree->Branch("KFKpXiProductionVtx_y", &event.KFprodvtx_y_kpxi);
  tree->Branch("KFKpXiProductionVtx_z", &event.KFprodvtx_z_kpxi);
  //extrapolation to the production vertex
  tree->Branch("KFXiProductionVtx_x_KpXi", &event.KFxi_kpxiprodvtx_x);
  tree->Branch("KFXiProductionVtx_y_KpXi", &event.KFxi_kpxiprodvtx_y);
  tree->Branch("KFXiProductionVtx_z_KpXi", &event.KFxi_kpxiprodvtx_z);
  tree->Branch("KFXiProductionVtxMom_KpXi", &event.KFxi_kpxiprodmom);
  tree->Branch("KFXiProductionVtxMom_x_KpXi", &event.KFxi_kpxiprodmom_x);
  tree->Branch("KFXiProductionVtxMom_y_KpXi", &event.KFxi_kpxiprodmom_y);
  tree->Branch("KFXiProductionVtxMom_z_KpXi", &event.KFxi_kpxiprodmom_z);
  tree->Branch("KFXiProductionVtxCloseDist_KpXi", &event.KFxi_kpxiprodvtx_dist);

  //Km, Kp vertex
  tree->Branch("KFXiProductionVtx_x_KK", &event.KFxi_kkvtx_x);
  tree->Branch("KFXiProductionVtx_y_KK", &event.KFxi_kkvtx_y);
  tree->Branch("KFXiProductionVtx_z_KK", &event.KFxi_kkvtx_z);
  tree->Branch("KFXiProductionVtxMom_KK", &event.KFxi_kkvtx_mom);
  tree->Branch("KFXiProductionVtxMom_x_KK", &event.KFxi_kkvtx_mom_x);
  tree->Branch("KFXiProductionVtxMom_y_KK", &event.KFxi_kkvtx_mom_y);
  tree->Branch("KFXiProductionVtxMom_z_KK", &event.KFxi_kkvtx_mom_z);
  tree->Branch("KFXiProductionVtxCloseDist_KK", &event.KFxi_kkvtx_dist);

  tree->Branch("KFDecaysMom", &event.KFdecays_mom);
  tree->Branch("KFDecaysMom_x", &event.KFdecays_mom_x);
  tree->Branch("KFDecaysMom_y", &event.KFdecays_mom_y);
  tree->Branch("KFDecaysMom_z", &event.KFdecays_mom_z);
  tree->Branch("KFDecaysMomCM", &event.KFdecays_CMmom);
  tree->Branch("KFDecaysMomCM_x", &event.KFdecays_CMmom_x);
  tree->Branch("KFDecaysMomCM_y", &event.KFdecays_CMmom_y);
  tree->Branch("KFDecaysMomCM_z", &event.KFdecays_CMmom_z);

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

  TTreeReaderCont[kE42] = new TTreeReader( "tpc", TFileCont[kE42] );
  const auto& reader = TTreeReaderCont[kE42];
  src.evnum = new TTreeReaderValue<Int_t>( *reader, "evnum" );

  src.nhHtof = new TTreeReaderValue<Int_t>( *reader, "nhHtof" );
  src.HtofSeg = new TTreeReaderValue<std::vector<Double_t>>( *reader, "HtofSeg" );
  src.tHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "tHtof" );
  src.dtHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "dtHtof" );
  src.deHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "deHtof" );
  src.posHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "posHtof" );
  src.G4tidHtof = new TTreeReaderValue<std::vector<Int_t>>( *reader, "G4tidHtof" );

  src.nclTpc = new TTreeReaderValue<Int_t>( *reader, "nclTpc" );
  src.remain_nclTpc = new TTreeReaderValue<Int_t>( *reader, "remain_nclTpc" );
#if SaveRawData
  src.cluster_x = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_x" );
  src.cluster_y = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_y" );
  src.cluster_z = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_z" );
  src.cluster_de = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_de" );
//  src.cluster_size = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_size" );
  src.cluster_layer = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_layer" );
  src.cluster_mrow = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_mrow" );
#if 0
	src.cluster_de_center = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_de_center" );
  src.cluster_x_center = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_x_center" );
  src.cluster_y_center = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_y_center" );
  src.cluster_z_center = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_z_center" );
  src.cluster_row_center = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_row_center" );
#endif
	src.cluster_houghflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_houghflag" );
#endif

  src.ntTpc = new TTreeReaderValue<Int_t>( *reader, "ntTpc" );
  src.nhtrack = new TTreeReaderValue<std::vector<Int_t>>( *reader, "nhtrack" );
  src.trackid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trackid" );
  src.isXi = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isXi" );
  src.isBeam = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isBeam" );
  src.isKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isKurama" );
  src.isK18 = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isK18" );
  src.isAccidental = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isAccidental" );
  src.isMultiloop = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isMultiloop" );
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
  src.pathhit = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "pathhit" );
  src.alpha = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "alpha" );
  src.track_cluster_de = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_de" );
//  src.track_cluster_size = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_size" );
  src.track_cluster_mrow = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_mrow" );

#if SaveTPCK18
  src.isgoodTPCK18 = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isgoodTPCK18" );
  src.chisqrTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrTPCK18" );
  src.qTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "qTPCK18");
  src.pTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pTPCK18");
  src.xtgtTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtTPCK18" );
  src.ytgtTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtTPCK18" );
  src.utgtTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtTPCK18" );
  src.vtgtTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtTPCK18" );
  src.thetaTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaTPCK18" );
  src.lhtofTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "lhtofTPCK18" );
  src.xhtofTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xhtofTPCK18" );
  src.yhtofTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yhtofTPCK18" );
  src.lvpTPCK18 = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "lvpTPCK18" );
  src.xvpTPCK18 = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "xvpTPCK18" );
  src.yvpTPCK18 = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "yvpTPCK18" );
  src.xhtofK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xhtofHS" );
  src.yhtofK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yhtofHS" );
#endif

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

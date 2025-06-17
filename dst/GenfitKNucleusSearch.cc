// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <string>
#include <fstream>
#include <filesystem>
#include <Math/ProbFunc.h>

#include <TGeoPhysicalConstants.h>

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
#include "FourVectorCartesianFitter.hh"
#include "FourVectorFitter.hh"

#include "HypTPCFitter.hh"
#include "HypTPCTask.hh"


#define SaveTPCK18 1
#define MMcut_Xi 0
#define DoKinematicFitLd 0
#define Debug 1
#define WriteOutput 0
#define SelectEvents 0
#define MomCutTemp 1

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

//For GenFit Setting
const Bool_t Const_field = false; //Must be false for linear tracking
const Int_t verbosity = 0;//0~3;
//const Int_t verbosity = 1;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");

const Double_t lambda_masscut = 0.02;
const Double_t p_vtx_distcut = 300;
const Double_t pi_vtx_distcut = 300;
const Double_t p2_vtx_distcut = 300;
//const Double_t xi_vtx_distcut = 100;
const Double_t ppi_distcut = 10.;
const Double_t pipi_distcut = 10.;
const Double_t kmdecay_distcut = 10.;
const Double_t k0target_distcut = 50.;
const Double_t lp2_distcut = 50.;
const Double_t vtx_scan_range = 150.;
//const Double_t vtx_scan_range = 50.;
//const Double_t vtx_scan_rangeInsideL = 50.;
//const Double_t vtx_scan_rangeInsidePi = 50.;
//const Double_t vtx_scan_rangeInsideL = 200.;
//const Double_t vtx_scan_rangeInsideL = 300.;
const Double_t vtx_scan_rangeInsideL = 50.;
const Double_t vtx_scan_rangeInsidePi = 50.;
const Double_t GFppi_distcut = 100.;
const Double_t GFlpi_distcut = 100.;
const Double_t GFlp2_distcut = 100.;
const Double_t GFkmdecay_distcut = 10.;
const Double_t GFkmdecay_hitcut = 30.;
const Double_t GFxitarget_distcut = 50.;
const Double_t GFltarget_distcut = 50.;
const Double_t GFxitarget_ycut = 20.;
const Double_t GFltarget_ycut = 20.;

const Double_t& HS_field_0 = ConfMan::Get<Double_t>("HSFLDCALIB");
const Double_t& HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
const Double_t& HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kE42, kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[E42]", "[OutFile]" };
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

  Int_t nhHtof;
  std::vector<Double_t> HtofSeg;
  std::vector<Double_t> tHtof;
  std::vector<Double_t> dtHtof;
  std::vector<Double_t> deHtof;
  std::vector<Double_t> posHtof;

  Int_t ntK18;
  std::vector<Double_t> pK18;
  std::vector<Double_t> thetaK18;  
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
  std::vector<Double_t> xtgtKurama;
  std::vector<Double_t> ytgtKurama;
  std::vector<Double_t> utgtKurama;
  std::vector<Double_t> vtgtKurama;

  std::vector<Double_t> momTransfer;
  std::vector<Double_t> thetaK18Kurama;  
  Double_t momDiff;

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
  std::vector<Int_t> Kflag;
  std::vector<Int_t> Pflag;

  //TPC RK
  Int_t ntTPCK18;
  std::vector<Int_t> isgoodTPCK18;
  std::vector<Int_t> tpcidTPCK18;
  std::vector<Double_t> chisqrTPCK18;
  std::vector<Double_t> pTPCK18;
  std::vector<Double_t> qTPCK18;
  std::vector<Double_t> thetaTPCK18;  
  std::vector<Double_t> xtgtTPCK18;
  std::vector<Double_t> ytgtTPCK18;
  std::vector<Double_t> utgtTPCK18;
  std::vector<Double_t> vtgtTPCK18;
  std::vector<Double_t> lhtofTPCK18;
  std::vector<Double_t> xhtofTPCK18;
  std::vector<Double_t> yhtofTPCK18;
  std::vector<std::vector<Double_t>> lvpTPCK18;
  std::vector<std::vector<Double_t>> xvpTPCK18;
  std::vector<std::vector<Double_t>> yvpTPCK18;

  std::vector<Int_t> isgoodTPCKurama;
  std::vector<Int_t> tpcidTPCKurama;
  std::vector<Int_t> kflagTPCKurama;
  std::vector<Int_t> pflagTPCKurama;
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
  std::vector<Double_t> MissMassNuclTPC;
  std::vector<Double_t> MissMassNuclCorrTPC;
  std::vector<Double_t> MissMassNuclCorrDETPC;
  std::vector<Double_t> BEkaonTPC;
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
  std::vector<Double_t> helix_cx;
  std::vector<Double_t> helix_cy;
  std::vector<Double_t> helix_z0;
  std::vector<Double_t> helix_r;
  std::vector<Double_t> helix_dz;
  std::vector<Double_t> dE;
  std::vector<Double_t> dEdx; //reference dedx
  std::vector<Double_t> mom0; //Helix momentum at Y = 0
  std::vector<Double_t> path; //Helix path
  std::vector<Int_t> isElectron;
  std::vector<Double_t> nsigma_triton;
  std::vector<Double_t> nsigma_deutron;
  std::vector<Double_t> nsigma_proton;
  std::vector<Double_t> nsigma_kaon;
  std::vector<Double_t> nsigma_pion;
  std::vector<Double_t> nsigma_electron;
  
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
  std::vector<std::vector<Double_t>> pathhit;
  std::vector<std::vector<Double_t>> alpha;
  std::vector<std::vector<Double_t>> track_cluster_de;
  std::vector<std::vector<Double_t>> track_cluster_mrow;

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

  std::vector<Double_t> angleLambda;

  // Bool_t xiflag;
  // Double_t ximass;
  // Double_t xidecayvtx_x;
  // Double_t xidecayvtx_y;
  // Double_t xidecayvtx_z;
  // Double_t ximom_x;
  // Double_t ximom_y;
  // Double_t ximom_z;
  // Double_t lpi_dist;

  Bool_t lflag;
  Bool_t forwardlflag;  
  Double_t lmass;
  Double_t ldecayvtx_x;
  Double_t ldecayvtx_y;
  Double_t ldecayvtx_z;
  Double_t lmom_x;
  Double_t lmom_y;
  Double_t lmom_z;
  Double_t ppi_dist;
  Double_t ppiangle;    
  std::vector<Int_t> ldecays_id;
  std::vector<Double_t> ldecays_mom;
  std::vector<Double_t> ldecays_mom_x;
  std::vector<Double_t> ldecays_mom_y;
  std::vector<Double_t> ldecays_mom_z;
  std::vector<Double_t> ldecays_res_mom;
  std::vector<Double_t> ldecays_res_mom_x;
  std::vector<Double_t> ldecays_res_mom_y;
  std::vector<Double_t> ldecays_res_mom_z;
  std::vector<Double_t> ldecays_res_mom_t;
  std::vector<Double_t> ldecays_res_th;
  std::vector<Double_t> ldecays_res_ph;
  std::vector<Double_t> ldecays_cov_mom_th;
  std::vector<Double_t> ldecays_cov_mom_ph;
  std::vector<Double_t> ldecays_cov_mom_xy;
  std::vector<Double_t> ldecays_cov_mom_yz;
  std::vector<Double_t> ldecays_cov_mom_zx;

  Bool_t lp2flag;
  Double_t lmass_lp2;
  Double_t lp2mass;
  Double_t lp2decayvtx_x;
  Double_t lp2decayvtx_y;
  Double_t lp2decayvtx_z;
  Double_t lp2mom_x;
  Double_t lp2mom_y;
  Double_t lp2mom_z;
  Double_t lmom_lp2_x;
  Double_t lmom_lp2_y;
  Double_t lmom_lp2_z;  
  Double_t lp2_dist;
  Double_t lp2angle;  
  std::vector<Int_t> lp2decays_id;
  std::vector<Double_t> lp2decays_mom;
  std::vector<Double_t> lp2decays_mom_x;
  std::vector<Double_t> lp2decays_mom_y;
  std::vector<Double_t> lp2decays_mom_z;
  std::vector<Double_t> lp2decays_res_mom;
  std::vector<Double_t> lp2decays_res_mom_x;
  std::vector<Double_t> lp2decays_res_mom_y;
  std::vector<Double_t> lp2decays_res_mom_z;
  std::vector<Double_t> lp2decays_res_mom_t;
  std::vector<Double_t> lp2decays_res_th;
  std::vector<Double_t> lp2decays_res_ph;
  std::vector<Double_t> lp2decays_cov_mom_th;
  std::vector<Double_t> lp2decays_cov_mom_ph;
  std::vector<Double_t> lp2decays_cov_mom_xy;
  std::vector<Double_t> lp2decays_cov_mom_yz;
  std::vector<Double_t> lp2decays_cov_mom_zx;

			 

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

  std::vector<Double_t> GFdecays_m2;
  std::vector<Double_t> GFdecays_mom;
  std::vector<Double_t> GFdecays_mom_x;
  std::vector<Double_t> GFdecays_mom_y;
  std::vector<Double_t> GFdecays_mom_z;
  std::vector<Double_t> GFmomloss;
  std::vector<Double_t> GFeloss;

  Int_t GFstatus;
  Int_t GFntTpc;
  std::vector<Int_t> GFinside;
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
  std::vector<std::vector<Double_t>> GFresidual_px;
  std::vector<std::vector<Double_t>> GFresidual_py;
  std::vector<std::vector<Double_t>> GFresidual_pz;
  std::vector<std::vector<Double_t>> GFresidual_p;
  std::vector<std::vector<Double_t>> GFresidual6D_x;
  std::vector<std::vector<Double_t>> GFresidual6D_y;
  std::vector<std::vector<Double_t>> GFresidual6D_z;
  std::vector<std::vector<Double_t>> GFresidual6D_px;
  std::vector<std::vector<Double_t>> GFresidual6D_py;
  std::vector<std::vector<Double_t>> GFresidual6D_pz;
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

  Bool_t GFlp2flag;
  Double_t GFlp2mass;
  Double_t GFlmass_lp2;  
  Double_t GFlp2decayvtx_x;
  Double_t GFlp2decayvtx_y;
  Double_t GFlp2decayvtx_z;
  Double_t GFlp2mom;
  Double_t GFlp2mom_x;
  Double_t GFlp2mom_y;
  Double_t GFlp2mom_z;
  Double_t GFlmom_lp2;  
  Double_t GFlp2_dist;
  Double_t GFlp2_angle;

  std::vector<Double_t> GFdecays_m2_lp2;
  std::vector<Double_t> GFdecays_mom_lp2;
  std::vector<Double_t> GFdecays_mom_x_lp2;
  std::vector<Double_t> GFdecays_mom_y_lp2;
  std::vector<Double_t> GFdecays_mom_z_lp2;
  std::vector<Double_t> GFmomloss_lp2;
  std::vector<Double_t> GFeloss_lp2;

  Int_t GFstatus_lp2;
  Int_t GFntTpc_lp2;
  std::vector<Int_t> GFfitstatus_lp2;
  std::vector<Int_t> GFpdgcode_lp2;
  std::vector<Int_t> GFnhtrack_lp2;
  std::vector<Double_t> GFcharge_lp2;
  std::vector<Double_t> GFchisqr_lp2;
  std::vector<Double_t> GFtof_lp2;
  std::vector<Double_t> GFtracklen_lp2;
  std::vector<Double_t> GFpval_lp2;
  std::vector<Double_t> GFchisqrPos_lp2;
  std::vector<Double_t> GFpvalPos_lp2;
  std::vector<std::vector<Double_t>> GFlayer_lp2;
  std::vector<std::vector<Double_t>> GFrow_lp2;
  std::vector<std::vector<Double_t>> GFpos_x_lp2;
  std::vector<std::vector<Double_t>> GFpos_y_lp2;
  std::vector<std::vector<Double_t>> GFpos_z_lp2;
  std::vector<std::vector<Double_t>> GFmom_lp2;
  std::vector<std::vector<Double_t>> GFmom_x_lp2;
  std::vector<std::vector<Double_t>> GFmom_y_lp2;
  std::vector<std::vector<Double_t>> GFmom_z_lp2;
  std::vector<std::vector<Double_t>> GFresidual_x_lp2;
  std::vector<std::vector<Double_t>> GFresidual_y_lp2;
  std::vector<std::vector<Double_t>> GFresidual_z_lp2;
  std::vector<std::vector<Double_t>> GFresidual_px_lp2;
  std::vector<std::vector<Double_t>> GFresidual_py_lp2;
  std::vector<std::vector<Double_t>> GFresidual_pz_lp2;
  std::vector<std::vector<Double_t>> GFresidual_p_lp2;
  std::vector<std::vector<Double_t>> GFresidual6D_x_lp2;
  std::vector<std::vector<Double_t>> GFresidual6D_y_lp2;
  std::vector<std::vector<Double_t>> GFresidual6D_z_lp2;
  std::vector<std::vector<Double_t>> GFresidual6D_px_lp2;
  std::vector<std::vector<Double_t>> GFresidual6D_py_lp2;
  std::vector<std::vector<Double_t>> GFresidual6D_pz_lp2;
  std::vector<std::vector<Double_t>> GFresolution_x_lp2;
  std::vector<std::vector<Double_t>> GFresolution_y_lp2;
  std::vector<std::vector<Double_t>> GFresolution_z_lp2;
  std::vector<std::vector<Double_t>> GFresolution_p_lp2;
  std::vector<std::vector<Double_t>> GFresolution_px_lp2;
  std::vector<std::vector<Double_t>> GFresolution_py_lp2;
  std::vector<std::vector<Double_t>> GFresolution_pz_lp2;
  std::vector<std::vector<Double_t>> GFpull_x_lp2;
  std::vector<std::vector<Double_t>> GFpull_y_lp2;
  std::vector<std::vector<Double_t>> GFpull_z_lp2;
  std::vector<std::vector<Double_t>> GFpull_p_lp2;
  std::vector<std::vector<Double_t>> GFpull_px_lp2;
  std::vector<std::vector<Double_t>> GFpull_py_lp2;
  std::vector<std::vector<Double_t>> GFpull_pz_lp2;

  std::vector<Double_t> GFtrack2vtxdist;
  std::vector<Double_t> GFcalctof;
  std::vector<Double_t> GFsegHtof;
  std::vector<Double_t> GFtofHtof;
  std::vector<Double_t> GFtdiffHtof;
  std::vector<Double_t> GFposHtof;
  std::vector<Double_t> GFposx;
  std::vector<Double_t> GFposy;
  std::vector<Double_t> GFposz;
  std::vector<Double_t> GFinvbeta;
  std::vector<Double_t> GFm2;  
  std::vector<Double_t> nsigma_tritonHtof;
  std::vector<Double_t> nsigma_deutronHtof;
  std::vector<Double_t> nsigma_protonHtof;
  std::vector<Double_t> nsigma_kaonHtof;
  std::vector<Double_t> nsigma_pionHtof;
  std::vector<Double_t> nsigma_electronHtof;

  Bool_t kmflag;
  Bool_t kmdecayflag;
  Bool_t kmmommatch;
  Int_t kmntrack;  
  Double_t kmmass;
  Double_t kmmass2;
  Double_t kmdecayvtx_x;
  Double_t kmdecayvtx_y;
  Double_t kmdecayvtx_z;
  std::vector<Double_t> kmmomdiff;  
  std::vector<Double_t> kmscatpmom;  
  std::vector<Double_t> kmscatpmom_x;
  std::vector<Double_t> kmscatpmom_y;
  std::vector<Double_t> kmscatpmom_z;      
  Double_t kmmom;  
  Double_t kmmom_x;
  Double_t kmmom_y;
  Double_t kmmom_z;
  std::vector<Int_t> km_id;
  std::vector<Int_t> kmdecays_id;    
  
  Double_t KFchisqrxi;
  Double_t KFpvalxi;  
  vector<double> KFxidecays_pull;
  // Double_t KFximass;//MassBefore;
  // Double_t KFximom;
  // Double_t KFximom_x;
  // Double_t KFximom_y;
  // Double_t KFximom_z;

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

  // bool XiAccidental;
  // bool XiFlight;
  // bool XiProd;

  // double xtgtXi;
  // double ytgtXi;
  // double utgtXi;
  // double vtgtXi;

  // double vtxKKXi;
  // double vtyKKXi;
  // double vtzKKXi;

  // double xiprodvtx_x;
  // double xiprodvtx_y;
  // double xiprodvtx_z;

  // double xiprodmom_x;
  // double xiprodmom_y;
  // double xiprodmom_z;


  void clear( void )
  {
    runnum = 0;
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

    ntK18 = 0;
    pK18.clear();
    thetaK18.clear();    
    chisqrK18.clear();
    xtgtK18.clear();
    ytgtK18.clear();
    utgtK18.clear();
    vtgtK18.clear();
    // xtgtHS.clear();
    // ytgtHS.clear();
    // ztgtHS.clear();
    // xvpHS.clear();
    // yvpHS.clear();
    // zvpHS.clear();

    ntKurama = 0;
    chisqrKurama.clear();
    pKurama.clear();
    qKurama.clear();
    m2.clear();
    xtgtKurama.clear();
    ytgtKurama.clear();
    utgtKurama.clear();
    vtgtKurama.clear();

    momTransfer.clear();
    thetaK18Kurama.clear();    
    momDiff = 0.;

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
    Kflag.clear();
    Pflag.clear();

    isgoodTPCK18.clear();
    tpcidTPCK18.clear();
    chisqrTPCK18.clear();
    qTPCK18.clear();
    pTPCK18.clear();
    thetaTPCK18.clear();    
    xtgtTPCK18.clear();
    ytgtTPCK18.clear();
    utgtTPCK18.clear();
    vtgtTPCK18.clear();
    lhtofTPCK18.clear();
    xhtofTPCK18.clear();
    yhtofTPCK18.clear();
    lvpTPCK18.clear();
    xvpTPCK18.clear();
    yvpTPCK18.clear();

    isgoodTPCKurama.clear();
    tpcidTPCKurama.clear();
    kflagTPCKurama.clear();
    pflagTPCKurama.clear();
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
    MissMassNuclTPC.clear();
    MissMassNuclCorrTPC.clear();
    MissMassNuclCorrDETPC.clear();
    BEkaonTPC.clear();
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

    ntTpc = 0;
    nhtrack.clear();
    trackid.clear();
    isBeam.clear();
    isKurama.clear();
    isK18.clear();
    isAccidental.clear();
    charge.clear();
    pid.clear();

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
    isElectron.clear();
    nsigma_triton.clear();
    nsigma_deutron.clear();
    nsigma_proton.clear();
    nsigma_kaon.clear();
    nsigma_pion.clear();
    nsigma_electron.clear();

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
    pathhit.clear();
    alpha.clear();
    track_cluster_de.clear();
    track_cluster_mrow.clear();

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

    angleLambda.clear();

    // xiflag = false;
    // ximass = qnan;
    // xidecayvtx_x = qnan;
    // xidecayvtx_y = qnan;
    // xidecayvtx_z = qnan;
    // ximom_x = qnan;
    // ximom_y = qnan;
    // ximom_z = qnan;
    // lpi_dist = qnan;

    lflag = false;
    forwardlflag = false;    
    lmass = qnan;
    ldecayvtx_x = qnan;
    ldecayvtx_y = qnan;
    ldecayvtx_z = qnan;
    lmom_x = qnan;
    lmom_y = qnan;
    lmom_z = qnan;
    ppi_dist = qnan;

    lp2flag = false;
    lmass_lp2 = qnan;
    lp2mass = qnan;
    lp2decayvtx_x = qnan;
    lp2decayvtx_y = qnan;
    lp2decayvtx_z = qnan;
    lp2mom_x = qnan;
    lp2mom_y = qnan;
    lp2mom_z = qnan;
    lmom_lp2_x = qnan;
    lmom_lp2_y = qnan;
    lmom_lp2_z = qnan;    
    lp2_dist = qnan;
    lp2angle = qnan;
    
    ldecays_id.clear();
    ldecays_mom.clear();
    ldecays_mom_x.clear();
    ldecays_mom_y.clear();
    ldecays_mom_z.clear();
    ldecays_res_mom.clear();
    ldecays_res_mom_x.clear();
    ldecays_res_mom_y.clear();
    ldecays_res_mom_z.clear();
    ldecays_res_mom_t.clear();
    ldecays_res_th.clear();
    ldecays_res_ph.clear();
    ldecays_cov_mom_th.clear();
    ldecays_cov_mom_ph.clear();
    ldecays_cov_mom_xy.clear();
    ldecays_cov_mom_yz.clear();
    ldecays_cov_mom_zx.clear();

    lp2decays_id.clear();
    lp2decays_mom.clear();
    lp2decays_mom_x.clear();
    lp2decays_mom_y.clear();
    lp2decays_mom_z.clear();
    lp2decays_res_mom.clear();
    lp2decays_res_mom_x.clear();
    lp2decays_res_mom_y.clear();
    lp2decays_res_mom_z.clear();
    lp2decays_res_mom_t.clear();
    lp2decays_res_th.clear();
    lp2decays_res_ph.clear();
    lp2decays_cov_mom_th.clear();
    lp2decays_cov_mom_ph.clear();
    lp2decays_cov_mom_xy.clear();
    lp2decays_cov_mom_yz.clear();
    lp2decays_cov_mom_zx.clear();

     
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

    GFdecays_m2.clear();
    GFdecays_mom.clear();
    GFdecays_mom_x.clear();
    GFdecays_mom_y.clear();
    GFdecays_mom_z.clear();
    GFmomloss.clear();
    GFeloss.clear();

    GFstatus = 0;
    GFntTpc = 0;
    GFinside.clear();
    GFcharge.clear();
    GFchisqr.clear();
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
    GFresidual_px.clear();
    GFresidual_py.clear();
    GFresidual_pz.clear();
    GFresidual_p.clear();
    GFresidual6D_x.clear();
    GFresidual6D_y.clear();
    GFresidual6D_z.clear();
    GFresidual6D_px.clear();
    GFresidual6D_py.clear();
    GFresidual6D_pz.clear();
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

    GFlp2flag = false;
    GFlp2mass = qnan;
    GFlmass_lp2 = qnan;    
    GFlp2decayvtx_x = qnan;
    GFlp2decayvtx_y = qnan;
    GFlp2decayvtx_z = qnan;
    GFlp2mom = qnan;
    GFlp2mom_x = qnan;
    GFlp2mom_y = qnan;
    GFlp2mom_z = qnan;
    GFlmom_lp2 = qnan;    
    GFlp2_dist = qnan;
    GFlp2_angle = qnan;

    GFdecays_m2_lp2.clear();
    GFdecays_mom_lp2.clear();
    GFdecays_mom_x_lp2.clear();
    GFdecays_mom_y_lp2.clear();
    GFdecays_mom_z_lp2.clear();
    GFmomloss_lp2.clear();
    GFeloss_lp2.clear();

    GFstatus_lp2 = 0;
    GFntTpc_lp2 = 0;
    GFcharge_lp2.clear();
    GFchisqr_lp2.clear();
    GFtof_lp2.clear();
    GFtracklen_lp2.clear();
    GFpval_lp2.clear();
    GFfitstatus_lp2.clear();
    GFpdgcode_lp2.clear();
    GFnhtrack_lp2.clear();
    GFlayer_lp2.clear();
    GFrow_lp2.clear();
    GFpos_x_lp2.clear();
    GFpos_y_lp2.clear();
    GFpos_z_lp2.clear();
    GFmom_lp2.clear();
    GFmom_x_lp2.clear();
    GFmom_y_lp2.clear();
    GFmom_z_lp2.clear();
    GFresidual_x_lp2.clear();
    GFresidual_y_lp2.clear();
    GFresidual_z_lp2.clear();
    GFresidual_px_lp2.clear();
    GFresidual_py_lp2.clear();
    GFresidual_pz_lp2.clear();
    GFresidual_p_lp2.clear();
    GFresidual6D_x_lp2.clear();
    GFresidual6D_y_lp2.clear();
    GFresidual6D_z_lp2.clear();
    GFresidual6D_px_lp2.clear();
    GFresidual6D_py_lp2.clear();
    GFresidual6D_pz_lp2.clear();
    GFresolution_x_lp2.clear();
    GFresolution_y_lp2.clear();
    GFresolution_z_lp2.clear();
    GFresolution_p_lp2.clear();
    GFresolution_px_lp2.clear();
    GFresolution_py_lp2.clear();
    GFresolution_pz_lp2.clear();
    GFpull_x_lp2.clear();
    GFpull_y_lp2.clear();
    GFpull_z_lp2.clear();
    GFpull_p_lp2.clear();
    GFpull_px_lp2.clear();
    GFpull_py_lp2.clear();
    GFpull_pz_lp2.clear();

    GFtrack2vtxdist.clear();
    GFcalctof.clear();
    GFsegHtof.clear();
    GFtofHtof.clear();
    GFtdiffHtof.clear();
    GFposHtof.clear();
    GFposx.clear();
    GFposy.clear();
    GFposz.clear();
    GFinvbeta.clear();
    GFm2.clear();
    nsigma_tritonHtof.clear();
    nsigma_deutronHtof.clear();
    nsigma_protonHtof.clear();
    nsigma_kaonHtof.clear();
    nsigma_pionHtof.clear();
    nsigma_electronHtof.clear();    

    kmflag = false;
    kmdecayflag = false;
    kmmommatch = false;
    kmntrack = qnan;
    kmmass = qnan;
    kmmass2 = qnan;
    kmdecayvtx_x = qnan;
    kmdecayvtx_y = qnan;
    kmdecayvtx_z = qnan;
    kmmom = qnan;    
    kmmom_x = qnan;
    kmmom_y = qnan;
    kmmom_z = qnan;

    kmmomdiff.clear();        
    kmscatpmom.clear();    
    kmscatpmom_x.clear();
    kmscatpmom_y.clear();
    kmscatpmom_z.clear();

    km_id.clear();
    kmdecays_id.clear();    

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

    // XiAccidental = false;
    // XiFlight = false;
    // XiProd = false;
    // xtgtXi = qnan;
    // ytgtXi = qnan;
    // utgtXi = qnan;
    // vtgtXi = qnan;
    // vtxKKXi = qnan;
    // vtyKKXi = qnan;
    // vtzKKXi = qnan;
    // xiprodvtx_x = qnan;
    // xiprodvtx_y = qnan;
    // xiprodvtx_z = qnan;
    // xiprodmom_x = qnan;
    // xiprodmom_y = qnan;
    // xiprodmom_z = qnan;
  }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* runnum;
  TTreeReaderValue<Int_t>* evnum;
  TTreeReaderValue<std::vector<Int_t>>* trigpat;
  TTreeReaderValue<std::vector<Int_t>>* trigflag;

  TTreeReaderValue<Int_t>* nhHtof;
  TTreeReaderValue<std::vector<Double_t>>* HtofSeg;
  TTreeReaderValue<std::vector<Double_t>>* tHtof;
  TTreeReaderValue<std::vector<Double_t>>* dtHtof;
  TTreeReaderValue<std::vector<Double_t>>* deHtof;
  TTreeReaderValue<std::vector<Double_t>>* posHtof;

  TTreeReaderValue<Int_t>* ntK18;
  TTreeReaderValue<std::vector<Double_t>>* pK18;
  TTreeReaderValue<std::vector<Double_t>>* thetaK18;  
  TTreeReaderValue<std::vector<Double_t>>* chisqrK18;
  TTreeReaderValue<std::vector<Double_t>>* xtgtK18;
  TTreeReaderValue<std::vector<Double_t>>* ytgtK18;
  TTreeReaderValue<std::vector<Double_t>>* utgtK18;
  TTreeReaderValue<std::vector<Double_t>>* vtgtK18;

  // TTreeReaderValue<std::vector<std::vector<Double_t>>>* xvpHS;
  // TTreeReaderValue<std::vector<std::vector<Double_t>>>* yvpHS;
  // TTreeReaderValue<std::vector<std::vector<Double_t>>>* zvpHS;
  // TTreeReaderValue<std::vector<Double_t>>* xtgtHS;
  // TTreeReaderValue<std::vector<Double_t>>* ytgtHS;
  // TTreeReaderValue<std::vector<Double_t>>* ztgtHS;

  TTreeReaderValue<Int_t>* ntKurama;
  TTreeReaderValue<std::vector<Double_t>>* chisqrKurama;
  TTreeReaderValue<std::vector<Double_t>>* pKurama;
  TTreeReaderValue<std::vector<Double_t>>* qKurama;
  TTreeReaderValue<std::vector<Double_t>>* m2;
  TTreeReaderValue<std::vector<Double_t>>* xtgtKurama;
  TTreeReaderValue<std::vector<Double_t>>* ytgtKurama;
  TTreeReaderValue<std::vector<Double_t>>* utgtKurama;
  TTreeReaderValue<std::vector<Double_t>>* vtgtKurama;

  TTreeReaderValue<Int_t>* nKm;
  TTreeReaderValue<Int_t>* nKp;
  TTreeReaderValue<Int_t>* nKK;
  TTreeReaderValue<std::vector<Int_t>>* inside;
  TTreeReaderValue<std::vector<Double_t>>* vtx;
  TTreeReaderValue<std::vector<Double_t>>* vty;
  TTreeReaderValue<std::vector<Double_t>>* vtz;
  TTreeReaderValue<std::vector<Double_t>>* closeDist;
  TTreeReaderValue<std::vector<Double_t>>* MissMass;
  TTreeReaderValue<std::vector<Double_t>>* MissMassCorr;
  TTreeReaderValue<std::vector<Double_t>>* MissMassCorrDE;
  TTreeReaderValue<std::vector<Int_t>>* Kflag;
  TTreeReaderValue<std::vector<Int_t>>* Pflag;

  TTreeReaderValue<Int_t>* ntTPCK18; // Number of Tracks
  TTreeReaderValue<std::vector<Int_t>>* tpcidTPCK18;
  TTreeReaderValue<std::vector<Int_t>>* isgoodTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* chisqrTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* qTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* pTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* thetaTPCK18;  
  TTreeReaderValue<std::vector<Double_t>>* xtgtTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* ytgtTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* utgtTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* vtgtTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* lhtofTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* xhtofTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* yhtofTPCK18;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* lvpTPCK18;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* xvpTPCK18;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* yvpTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* xhtofK18;
  TTreeReaderValue<std::vector<Double_t>>* yhtofK18;
  // std::vector<std::vector<Double_t>> xvpHS;
  // std::vector<std::vector<Double_t>> yvpHS;
  // std::vector<std::vector<Double_t>> zvpHS;
  // std::vector<Double_t> xtgtHS;
  // std::vector<Double_t> ytgtHS;
  // std::vector<Double_t> ztgtHS;

  TTreeReaderValue<Int_t>* ntTPCKurama; // Number of Tracks
  TTreeReaderValue<std::vector<Int_t>>* tpcidTPCKurama;
  TTreeReaderValue<std::vector<Int_t>>* isgoodTPCKurama;
  TTreeReaderValue<std::vector<Int_t>>* kflagTPCKurama;
  TTreeReaderValue<std::vector<Int_t>>* pflagTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* chisqrTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* pTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* qTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* m2TPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* xtgtTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* ytgtTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* utgtTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* vtgtTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* thetaTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* pathTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* lhtofTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* xhtofTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* yhtofTPCKurama;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* lvpTPCKurama;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* xvpTPCKurama;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* yvpTPCKurama;

  TTreeReaderValue<std::vector<Int_t>>* isgoodTPC;
  TTreeReaderValue<std::vector<Int_t>>* insideTPC;
  TTreeReaderValue<std::vector<Double_t>>* vtxTPC;
  TTreeReaderValue<std::vector<Double_t>>* vtyTPC;
  TTreeReaderValue<std::vector<Double_t>>* vtzTPC;
  TTreeReaderValue<std::vector<Double_t>>* closeDistTPC;
  TTreeReaderValue<std::vector<Double_t>>* MissMassTPC;
  TTreeReaderValue<std::vector<Double_t>>* MissMassCorrTPC;
  TTreeReaderValue<std::vector<Double_t>>* MissMassCorrDETPC;
  TTreeReaderValue<std::vector<Double_t>>* MissMassNuclTPC;
  TTreeReaderValue<std::vector<Double_t>>* MissMassNuclCorrTPC;
  TTreeReaderValue<std::vector<Double_t>>* MissMassNuclCorrDETPC;
  TTreeReaderValue<std::vector<Double_t>>* pOrgTPC;
  TTreeReaderValue<std::vector<Double_t>>* pCalcTPC;
  TTreeReaderValue<std::vector<Double_t>>* pCorrTPC;
  TTreeReaderValue<std::vector<Double_t>>* pCorrDETPC;
  TTreeReaderValue<std::vector<Double_t>>* thetaTPC;
  TTreeReaderValue<std::vector<Double_t>>* thetaCMTPC;
  TTreeReaderValue<std::vector<Double_t>>* costCMTPC;
  TTreeReaderValue<std::vector<Double_t>>* ubTPC;
  TTreeReaderValue<std::vector<Double_t>>* vbTPC;
  TTreeReaderValue<std::vector<Double_t>>* usTPC;
  TTreeReaderValue<std::vector<Double_t>>* vsTPC;
  
  TTreeReaderValue<Int_t>* ntTpc; // Number of Tracks
  TTreeReaderValue<std::vector<Int_t>>* nhtrack; // Number of Hits (in 1 tracks)
  TTreeReaderValue<std::vector<Int_t>>* trackid; //for Kurama K1.8 tracks
  TTreeReaderValue<std::vector<Int_t>>* isBeam;
  TTreeReaderValue<std::vector<Int_t>>* isKurama;
  TTreeReaderValue<std::vector<Int_t>>* isK18;
  TTreeReaderValue<std::vector<Int_t>>* isAccidental;
  TTreeReaderValue<std::vector<Int_t>>* charge;//Helix charge
  TTreeReaderValue<std::vector<Int_t>>* pid;
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
  TTreeReaderValue<std::vector<Int_t>>* isElectron;  
  TTreeReaderValue<std::vector<Double_t>>* nsigma_triton;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_deutron;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_proton;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_kaon;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_pion;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_electron;  
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
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_mrow;

  TTreeReaderValue<Int_t>* nvtxTpc;
  TTreeReaderValue<std::vector<Double_t>>* vtx_x;
  TTreeReaderValue<std::vector<Double_t>>* vtx_y;
  TTreeReaderValue<std::vector<Double_t>>* vtx_z;
  TTreeReaderValue<std::vector<Double_t>>* vtx_dist;
  TTreeReaderValue<std::vector<Double_t>>* vtx_angle;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxid;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxmom_theta;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxpos_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxpos_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxpos_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxmom_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxmom_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxmom_z;
  
  TTreeReaderValue<std::vector<Double_t>>* angleLambda;

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
Bool_t SelectEvent( Int_t event )
{
  // input file
  std::stringstream file;
  file.str("");
  file << "output/tpc/proHws/run05588_GenfitKNucleusSearch_output.txt";
  std::string inputfile = file.str().c_str();
  std::vector<std::vector<int>> values;
  std::ifstream ifs(inputfile);
  if(!ifs){
    std::cerr << "Failed to open file: " << file.str().c_str() << std::endl;
    return true;
  }
  std::string line;
  while ( std::getline(ifs, line) ) {
    std::stringstream ss(line);
    std::vector<int> rowValues;
    int value;
    while ( ss >> value ) {
      rowValues.push_back(value);
    }
    values.push_back(rowValues);
  }
  
  for(const auto& row: values){
    if( row[0]==event ) return true;
  }
  
  return false;
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
  
#if WriteOutput
  std::cout << "Write output file mode" << std::endl;
#endif
#if SelectEvents
  std::cout << "select event mode" << std::endl;
#endif
  
  Int_t ievent = skip;
  for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
    gCounter.check();
    InitializeEvent();
#if SelectEvents
    if( !SelectEvent( ievent ) ) continue;
    std::cout << "Selected Event : event " << std::setw(6)
	      << ievent << std::endl;
#endif   
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
  static const auto ElectronMass = pdg::ElectronMass();
  static const Double_t Carbon12Mass = 12.*TGeoUnit::amu_c2 - 6.*ElectronMass;
  static const Double_t Boron11Mass  = 11.009305167*TGeoUnit::amu_c2 - 5.*ElectronMass;
  
  //if( ievent%1000==0 ){
  if( ievent%10==0 ){
    std::cout << "#D Event Number: "
        << std::setw(6) << ievent << std::endl;
  }  

  TVector3 tgtpos(0, 0, tpc::ZTarget);

  GetEntry(ievent);
  event.nKm = **src.nKm;
  
  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;

  event.nKK = **src.nKK;
  event.Pflag = **src.Pflag;
  event.MissMass = **src.MissMass;
  event.MissMassCorr = **src.MissMassCorr;
  event.MissMassCorrDE = **src.MissMassCorrDE;

  event.vtx = **src.vtx;
  event.vty = **src.vty;
  event.vtz = **src.vtz;
  event.closeDist = **src.closeDist;
  
  event.ntK18 = **src.ntK18;
  event.ntTPCK18 = **src.ntTPCK18;
  event.chisqrK18 = **src.chisqrK18;
  event.pK18 = **src.pK18;
  event.thetaK18 = **src.thetaK18;  
  event.xtgtK18 = **src.xtgtK18;
  event.ytgtK18 = **src.ytgtK18;
  event.utgtK18 = **src.utgtK18;
  event.vtgtK18 = **src.vtgtK18;
  // event.xtgtHS = **src.xtgtHS;
  // event.ytgtHS = **src.ytgtHS;
  // event.ztgtHS = **src.ztgtHS;
  // event.xvpHS = **src.xvpHS;
  // event.yvpHS = **src.yvpHS;
  // event.zvpHS = **src.zvpHS;
#if SaveTPCK18
  event.isgoodTPCK18 = **src.isgoodTPCK18;
  event.chisqrTPCK18 = **src.chisqrTPCK18;
  event.pTPCK18 = **src.pTPCK18;
  event.qTPCK18 = **src.qTPCK18;
  event.thetaTPCK18 = **src.thetaTPCK18;  
  event.xtgtTPCK18 = **src.xtgtTPCK18;
  event.ytgtTPCK18 = **src.ytgtTPCK18;
  event.utgtTPCK18 = **src.utgtTPCK18;
  event.vtgtTPCK18 = **src.vtgtTPCK18;
  event.lhtofTPCK18 = **src.lhtofTPCK18;
  event.xhtofTPCK18 = **src.xhtofTPCK18;
  event.yhtofTPCK18 = **src.yhtofTPCK18;
  event.lvpTPCK18 = **src.lvpTPCK18;
  event.xvpTPCK18 = **src.xvpTPCK18;
  event.yvpTPCK18 = **src.yvpTPCK18;
#endif
  event.ntKurama = **src.ntKurama;
  event.chisqrKurama = **src.chisqrKurama;
  event.pKurama = **src.pKurama;
  event.qKurama = **src.qKurama;
  event.xtgtKurama = **src.xtgtKurama;
  event.ytgtKurama = **src.ytgtKurama;
  event.utgtKurama = **src.utgtKurama;
  event.vtgtKurama = **src.vtgtKurama;
  event.tpcidTPCKurama = **src.tpcidTPCKurama;
  event.isgoodTPCKurama = **src.isgoodTPCKurama;
  event.kflagTPCKurama = **src.kflagTPCKurama;
  event.pflagTPCKurama = **src.pflagTPCKurama;
  event.chisqrTPCKurama = **src.chisqrTPCKurama;
  event.pTPCKurama = **src.pTPCKurama;
  event.qTPCKurama = **src.qTPCKurama;
  event.m2TPCKurama = **src.m2TPCKurama;
  event.xtgtTPCKurama = **src.xtgtTPCKurama;
  event.ytgtTPCKurama = **src.ytgtTPCKurama;
  event.utgtTPCKurama = **src.utgtTPCKurama;
  event.vtgtTPCKurama = **src.vtgtTPCKurama;
  event.thetaTPCKurama = **src.thetaTPCKurama;

  event.pathTPCKurama = **src.pathTPCKurama;
  event.lhtofTPCKurama = **src.lhtofTPCKurama;
  event.xhtofTPCKurama = **src.xhtofTPCKurama;
  event.yhtofTPCKurama = **src.yhtofTPCKurama;
  event.lvpTPCKurama = **src.lvpTPCKurama;
  event.xvpTPCKurama = **src.xvpTPCKurama;
  event.yvpTPCKurama = **src.yvpTPCKurama;

  event.isgoodTPC = **src.isgoodTPC;
  event.insideTPC = **src.insideTPC;
  event.vtxTPC = **src.vtxTPC;
  event.vtyTPC = **src.vtyTPC;
  event.vtzTPC = **src.vtzTPC;
  event.closeDistTPC = **src.closeDistTPC;
  event.MissMassTPC = **src.MissMassTPC;
  event.MissMassCorrTPC = **src.MissMassCorrTPC;
  event.MissMassCorrDETPC = **src.MissMassCorrDETPC;
  event.MissMassNuclTPC = **src.MissMassNuclTPC;
  event.MissMassNuclCorrTPC = **src.MissMassNuclCorrTPC;
  event.MissMassNuclCorrDETPC = **src.MissMassNuclCorrDETPC;
  event.pOrgTPC = **src.pOrgTPC;
  event.pCalcTPC = **src.pCalcTPC;
  event.pCorrTPC = **src.pCorrTPC;
  event.pCorrDETPC = **src.pCorrDETPC;
  event.thetaTPC = **src.thetaTPC;
  event.ubTPC = **src.ubTPC;
  event.vbTPC = **src.vbTPC;
  event.usTPC = **src.usTPC;
  event.vsTPC = **src.vsTPC;


  Bool_t FlagOutput = false;
  std::stringstream file;
  file.str("");
  std::filesystem::path full_path = __FILE__;
  std::string file_name = full_path.stem().string();
  
  file << "output/tpc/proHws/run0" << event.runnum << "_" << file_name << "_output.txt"; 
  std::string outputfile = file.str().c_str();
#if WriteOutput
  FlagOutput = true;
#endif
  if(FlagOutput){
    if(ievent==0){
      std::ofstream ofs(outputfile);
      ofs.open(outputfile);
      if(!ofs.is_open()){
	std::cerr << "Failed to open file: " << file.str().c_str() << std::endl;
	//return 1;
      }
      ofs << __FILE__ << " run:" << event.runnum <<" output file" << std::endl;
      ofs.close();
    }
  }
  Bool_t momcut = false;
#if MomCutTemp  
  momcut = true;
#endif

  Bool_t FlagDebug = false;
#if Debug
  std::cout << "Debug: Start searching p1,pi,p2" << std::endl;
  FlagDebug = true;
#endif
  
  
  if(FlagDebug) std::cout << "debug " << __LINE__  << std::endl;
  
  if( event.nKK != 1 || event.Pflag[0]!=1 ) return true;
  HF1( 4, event.MissMassTPC[0] );
  HF1( 5, event.MissMass[0] );
  HF1( 6, event.MissMassCorrTPC[0] );
  HF1( 7, event.MissMassCorrDETPC[0] );
  event.BEkaonTPC.resize(event.nKK);
  for(Int_t id=0; id<event.nKK; id++){
    event.BEkaonTPC[id] = event.MissMassNuclCorrDETPC[id] - KaonMass - Boron11Mass - 0.075;
  }
  
  if(FlagDebug) std::cout << "debug " << __LINE__  << std::endl;
  
  Int_t KuramaTrackId = -1; Int_t K18TrackId = -1;
  for(Int_t it=0;it<event.ntKurama;it++){
    if(event.isgoodTPCKurama[it]==1){
      KuramaTrackId = event.tpcidTPCKurama[it];
    }
  }
  event.momTransfer.resize(event.nKK);
  event.thetaK18Kurama.resize(event.nKK);  
  TVector3 momKurama[event.nKK];
  TVector3 momK18[event.nKK];
  TVector3 momTrans[event.nKK];
  //  Double_t momTransMag[event.nKK];
  std::vector<Double_t> momTransMag;
  //  for(int itK18=0; itK18<event.ntK18; itK18++){
  for(int itK18=0; itK18<event.ntK18; itK18++){
    Double_t pmagK18 = event.pK18[itK18];
    if(pmagK18>0.025) HF1(50,pmagK18);
    //    HF2( 90, pmagK18, event.thetaK18[itK18] );
    Double_t thetaK18 = event.thetaK18[itK18]; //to be used 
    //Double_t thetaK18 = event.thetaTPCK18[itK18];
    Double_t thetaK18Rad = thetaK18*TMath::DegToRad();
    HF2( 90, pmagK18, thetaK18 );   
    Double_t uK18 = event.utgtK18[itK18];
    Double_t vK18 = event.vtgtK18[itK18];
    for(int itKurama=0; itKurama<event.ntKurama; itKurama++){
      Int_t j=itK18*event.ntKurama+itKurama;
      momK18[j].SetXYZ( pmagK18*TMath::Cos(thetaK18Rad)*uK18,
			pmagK18*TMath::Cos(thetaK18Rad)*vK18,
			pmagK18*TMath::Cos(thetaK18Rad)*1.0  );      
      Double_t pmagKrm = event.pTPCKurama[itKurama];
      if(itK18==0 and pmagKrm>0.025) HF1(51,pmagKrm);
      if(itK18==0) HF2( 91, pmagKrm, event.thetaTPCKurama[itKurama] );
      Double_t thetaKrm = event.thetaTPCKurama[itKurama];
      Double_t thetaKrmRad = thetaKrm*TMath::DegToRad();
      Double_t uKrm = event.utgtTPCKurama[itKurama];
      Double_t vKrm = event.vtgtTPCKurama[itKurama];
      momKurama[j].SetXYZ( pmagKrm*TMath::Cos(thetaKrmRad)*uKrm,
			   pmagKrm*TMath::Cos(thetaKrmRad)*vKrm,
			   pmagKrm*TMath::Cos(thetaKrmRad)*1.0  );
      momTrans[j] = momK18[j] - momKurama[j];
      double mt = momTrans[j].Mag();
      momTransMag.push_back(mt);
      if(pmagK18>0.025 and pmagKrm>0.025) HF1( 52, mt );
      double thetaK18Kurama = momKurama[j].Angle(momK18[j])*TMath::RadToDeg();
      if(pmagK18>0.025 and pmagKrm>0.025){
	HF2( 92, mt, thetaK18Kurama );
	event.momTransfer.push_back(mt);
	event.thetaK18Kurama.push_back(thetaK18Kurama);
      }
    }
  }
  
  if(FlagDebug) std::cout << "debug " << __LINE__  << std::endl;
  
  // Fill histograms without lflag or lpflag
  for(int i=0; i<event.nKK; ++i){
    if(momK18[i].Mag()<0.025) continue;      
    if(momKurama[i].Mag()<0.025) continue;
    Double_t bek = event.BEkaonTPC[i];
    Double_t theta = event.thetaTPC[i];
    HF1( 1003, theta ); //Scattering Angle of Kp 
    HF1( 1200, theta ); //Scattering Angle of Kp
    HF1( 1201, bek ); // Inclusive Spectrum
    HF1( 11201, bek ); // Inclusive Spectrum
    if( theta>3.5&&theta<4.5){
      HF1( 1190, bek );
    }
    if( theta<5. ){
      HF1( 1202, bek );
      HF1( 11202, bek );
    }
    if( theta<10. ){
      HF1( 1203, bek );
      HF1( 11203, bek );
    }
  }
  
  if(FlagDebug) std::cout << "debug " << __LINE__  << std::endl;
  
  event.nhHtof = **src.nhHtof;
  event.HtofSeg = **src.HtofSeg;
  event.tHtof = **src.tHtof;
  event.dtHtof = **src.dtHtof;
  event.deHtof = **src.deHtof;
  event.posHtof = **src.posHtof;
 
  Int_t ntTpc = **src.ntTpc;
  if( ntTpc == 0 ) return true;
  //  HF1( 1, event.status++ ); // status 1
  HF1( 1, 1 ); // status 1

  event.ntTpc = ntTpc;
  event.nhtrack = **src.nhtrack;
  event.trackid = **src.trackid;
  event.isBeam = **src.isBeam;
  event.isKurama = **src.isKurama;
  event.isK18 = **src.isK18;
  event.isAccidental = **src.isAccidental;
  event.charge = **src.charge;
  event.pid = **src.pid;
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
  event.isElectron = **src.isElectron;
  event.nsigma_triton = **src.nsigma_triton;
  event.nsigma_deutron = **src.nsigma_deutron;
  event.nsigma_proton = **src.nsigma_proton;
  event.nsigma_kaon = **src.nsigma_kaon;
  event.nsigma_pion = **src.nsigma_pion;
  event.nsigma_electron = **src.nsigma_electron;
  
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
  event.track_cluster_mrow = **src.track_cluster_mrow;

  event.nvtxTpc = **src.nvtxTpc;
  event.vtx_x = **src.vtx_x;
  event.vtx_y = **src.vtx_y;
  event.vtx_z = **src.vtx_z;
  event.vtx_dist = **src.vtx_dist;
  event.vtx_angle = **src.vtx_angle;
  event.vtxid = **src.vtxid;
  event.vtxmom_theta = **src.vtxmom_theta;
  event.vtxpos_x = **src.vtxpos_x;
  event.vtxpos_y = **src.vtxpos_y;
  event.vtxpos_z = **src.vtxpos_z;
  event.vtxmom_x = **src.vtxmom_x;
  event.vtxmom_y = **src.vtxmom_y;
  event.vtxmom_z = **src.vtxmom_z;

  event.angleLambda = **src.angleLambda; 

  const Double_t ub_off = 0.0;
  const Double_t vb_off = 0.0;
  const Double_t xb_off = 0.0;
  const Double_t yb_off = 0.0;

  //Kaonic nucleus candidates searching
  Int_t l_candidates = 0;
  Int_t lp2_candidates = 0;
  Int_t forwardl_candidates = 0;
  
  Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  std::vector<Int_t> lp2_p_container, lp2_pi_container, lp2_p2_container, lp2_l_container;
  std::vector<Int_t> l_p_container, l_pi_container;
  std::vector<TVector3> lp2_mom_container, lp2_vert_container;
  std::vector<TVector3> l_mom_container, l_vert_container;
  //  std::vector<Double_t> xi_mass_container, lambda_mass_container;
  std::vector<Double_t> lambda_mass_container, lp2mass_container, lambda_lp2_mass_container;
  std::vector<Double_t> ppiangle_container,lp2angle_container;
  std::vector<TVector3> l_p_mom_container, l_pi_mom_container, l_l_mom_container,
    lp2_p_mom_container, lp2_pi_mom_container, lp2_p2_mom_container, lp2_l_mom_container;
  std::vector<TVector3> kurama_p_mom_container[event.ntKurama];
  std::vector<Double_t> ppi_closedist; std::vector<Double_t> lp2_closedist;

  Int_t counter[6]; // K+/-, pi+/-, p+/-
  for(int i=0; i<6; i++){
    counter[i]=0;
  } 
  for(int it1=0; it1<ntTpc; it1++){
    if(event.isElectron[it1]==1) continue;
    if(event.isK18[it1]==1 || event.isKurama[it1]==1 || event.isAccidental[it1]==1 ) continue;
    if((event.pid[it1]&4)==4){ // proton like
      if(event.charge[it1]<0) counter[0]++; // anti-proton
      else counter[1]++; // proton
    } 
    if((event.pid[it1]&2)==2){ // kaon like
      if(event.charge[it1]<0) counter[2]++; // k-
      else counter[3]++; // k+
    }
    if((event.pid[it1]&1)==1){ // pion like
      if(event.charge[it1]<0) counter[4]++; // pi-
      else counter[5]++; // pi+
    }
    // for dEdx vs mom
    for(int i=0; i<event.nKK; ++i){
      Double_t bek = event.BEkaonTPC[it1];
      Double_t theta = event.thetaTPC[it1];
      if(momK18[i].Mag()<0.025) continue;      
      if(momKurama[i].Mag()<0.025) continue;
      if(bek<-0.1){
	HF2(170, event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
      } else if(bek<0.){
	HF2(171, event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
      } else if(bek<0.1){
	HF2(172, event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
      } else if(bek<0.2){
	HF2(173, event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
      } else if(bek<0.3){
	HF2(174, event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
      }
      if(bek<0.3){
	HF2(180, event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
	if(event.charge[it1]<0){
	  HF2(181, -event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
	} else {
	  HF2(182, event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
	}
	if(theta>3.5&&theta<4.5){
	  HF2(183, event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
	  if(event.charge[it1]<0){
	    HF2(184, -event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
	  } else {
	    HF2(185, event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
	  }
	}
      }
    }    
    if(event.nsigma_proton[it1]>3 || event.nsigma_proton[it1]<-2) continue;
    for(int it2=0; it2<ntTpc; it2++){
      if(event.isElectron[it2]==1) continue;      
      if(event.isK18[it2]==1 || event.isKurama[it2]==1 || event.isAccidental[it2]==1 ) continue;
      if(event.nsigma_pion[it2]>3  || event.nsigma_pion[it2]<-3) continue;
      for(int it3=0; it3<ntTpc; it3++){
	if(event.isElectron[it3]==1) continue;
	if(event.isK18[it3]==1 || event.isKurama[it3]==1 || event.isAccidental[it3]==1 ) continue;
	if(event.nsigma_proton[it3]>3  || event.nsigma_proton[it3]<-2) continue;
	HF2(120, event.mom0[it1]*event.charge[it1], event.dEdx[it1]);      
	HF2(123, event.mom0[it2]*event.charge[it2], event.dEdx[it2]);
	HF2(126, event.mom0[it3]*event.charge[it3], event.dEdx[it3]);
	if ( event.charge[it1]<0 ){
	  HF2(121, -event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
	} else {
	  HF2(122, event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
	}
	if ( event.charge[it2]<0 ){
	  HF2(124, -event.mom0[it2]*event.charge[it2], event.dEdx[it2]);
	} else {
	  HF2(125, event.mom0[it2]*event.charge[it2], event.dEdx[it2]);
	}
	if ( event.charge[it3]<0 ){
	  HF2(127, -event.mom0[it3]*event.charge[it3], event.dEdx[it3]);
	} else {
	  HF2(128, event.mom0[it3]*event.charge[it3], event.dEdx[it3]);
	}      	
      }    
    }
  }  

  Bool_t FlagPidNsigma = true;
  //  Bool_t FlagPidNsigma = false;  
  Bool_t FlagKuramaProton = false;
  // proton1 loop
  for(Int_t it1=0;it1<ntTpc;it1++){
    Int_t itp1 = it1;
    Int_t itpi = -1;
    Int_t itp2 = -1;
    HF1( 1051, event.nhtrack[it1] );
    Bool_t p2fromLambda = false;
    if(event.isK18[it1]==1 || event.isKurama[it1]==1 || event.isAccidental[it1]==1 ) continue;
    Bool_t p_like = false;
    if(!FlagPidNsigma){
      if((event.pid[it1]&4)!=4) continue; // proton like
      if(event.charge[it1]==1) p_like = true;
    } else {
      if(event.nsigma_proton[it1]>5 || event.nsigma_proton[it1]<-2) continue;
      if(event.charge[it1]==1) p_like = true;      
    }
    if(FlagDebug) std::cout << "debug: p1 candidate found" << std::endl;
    if(event.charge[it1]==1) p_like = true; // select proton like
    // dEdx with one proton meaured
    HF2(140, event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
    if(event.charge[it1]<0){
      HF2(141, -event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
    } else {
      HF2(142, event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
    }    
    HF1( 1052, event.nhtrack[it1] );
    
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
    // pion1 loop
    for(Int_t it2=0;it2<ntTpc;it2++){
      itpi = it2;
      if(it1==it2) continue;
      if(event.isK18[it2]==1 || event.isKurama[it2]==1 || event.isAccidental[it2]==1 ) continue;
      HF1( 1053, event.nhtrack[it2] );
      Bool_t pim_like = false;
      if(!FlagPidNsigma){
	if((event.pid[it2]&1)!=1) continue; //select pi like	
	if(event.charge[it2]==-1) pim_like = true;
      } else {
	if(event.nsigma_pion[it2]>3 || event.nsigma_pion[it2]<-3) continue;
	if(event.charge[it2]==-1) pim_like = true;
      }
      if(FlagDebug) std::cout << "debug: pi candidate found" << std::endl;      
      // dEdx with one proton meaured and one pion 
      HF2(150, event.mom0[it2]*event.charge[it2], event.dEdx[it2]);
      if(event.charge[it2]<0){
	HF2(151, -event.mom0[it2]*event.charge[it2], event.dEdx[it2]);
      } else {
	HF2(152, event.mom0[it2]*event.charge[it2], event.dEdx[it2]);
      }
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
      if(TMath::Abs(lambda_vert.x()) > 250. || TMath::Abs(lambda_vert.z()) > 250. || TMath::Abs(lambda_vert.y()) > 250.) continue; //Vertex Cut
      if(ppi_dist > ppi_distcut) continue;
      Double_t pi_vertex_dist; Double_t p_vertex_dist;
      if(!Kinematics::HelixDirection(lambda_vert, p_start, p_end, p_vertex_dist) ||
	 !Kinematics::HelixDirection(lambda_vert, pi_start, pi_end, pi_vertex_dist)) continue;
      if(pi_vertex_dist > pi_vtx_distcut) continue;
      if(p_vertex_dist > p_vtx_distcut) continue;
      if(TMath::Abs(Llambda.M() - LambdaMass) > lambda_masscut) continue;
      event.lflag = true;
      if(FlagDebug) std::cout << "debug: Lambda reconstructed " << std::endl;
      Double_t thetap1 = p_mom.Theta()*TMath::RadToDeg();
      Double_t thetapi = pi_mom.Theta()*TMath::RadToDeg();
      Double_t thetaL = lambda_mom.Theta()*TMath::RadToDeg();
      // HF2( 2601, thetap1, p_mom.Mag() );
      // HF2( 2602, thetapi, pi_mom.Mag() );      
      // HF2( 2603, thetaL, lambda_mom.Mag() );
      // // theta corr
      // HF2( 2411, thetap1, thetap1 );
      // HF2( 2412, thetap1, thetapi );
      // HF2( 2413, thetap1, thetaL );
      // HF2( 2421, thetapi, thetapi );
      // HF2( 2422, thetapi, thetaL );
      // HF2( 2431, thetaL, thetaL );
      // //mom corr
      // HF2( 2211, p_mom.Mag(), p_mom.Mag() );
      // HF2( 2212, p_mom.Mag(), pi_mom.Mag() );
      // HF2( 2213, p_mom.Mag(), lambda_mom.Mag() );
      // HF2( 2221, pi_mom.Mag(), pi_mom.Mag() );
      // HF2( 2222, pi_mom.Mag(), lambda_mom.Mag() );
      // HF2( 2231, lambda_mom.Mag(), lambda_mom.Mag());
      if(FlagKuramaProton){
	// for(int itKurama=0; itKurama<event.ntKurama; itKurama++){
	//   Double_t thetap0 = momPKurama[itKurama].Theta()*TMath::RadToDeg();
	//   HF2( 2604, thetap0, momPKurama[itKurama].Mag() );
	//   Double_t anglep0p1 = p_mom.Angle(momPKurama[itKurama])*TMath::RadToDeg();
	//   HF1( 1160, anglep0p1 );
	//   Double_t anglep0pi = pi_mom.Angle(momPKurama[itKurama])*TMath::RadToDeg();
	//   HF1( 1161, anglep0pi );
	//   Double_t anglep0L = lambda_mom.Angle(momPKurama[itKurama])*TMath::RadToDeg();
	//   HF1( 1162, anglep0L );
	//   // theta corr
	//   HF2( 2401, thetap0, thetap0);
	//   HF2( 2402, thetap0, thetap1);
	//   HF2( 2403, thetap0, thetapi);
	//   HF2( 2404, thetap0, thetaL);
	//   // mom corr
	//   HF2( 2201, momPKurama[itKurama].Mag(), momPKurama[itKurama].Mag());
	//   HF2( 2202, momPKurama[itKurama].Mag(), p_mom.Mag());
	//   HF2( 2203, momPKurama[itKurama].Mag(), pi_mom.Mag());
	//   HF2( 2204, momPKurama[itKurama].Mag(),  lambda_mom.Mag());	  
	//  }
      }
      {
	Double_t anglep1pi = p_mom.Angle(pi_mom)*TMath::RadToDeg();
	Double_t anglep1L = p_mom.Angle(lambda_mom)*TMath::RadToDeg();
	Double_t anglepiL = pi_mom.Angle(lambda_mom)*TMath::RadToDeg();
	// HF1( 1163, anglep1pi );
	// HF1( 1164, anglep1L );	
	// HF1( 1165, anglepiL );
	// 
      }
      // HF2( 2501, thetap1, event.nhtrack[it1] );
      // HF2( 2502, thetapi, event.nhtrack[it2] );
      // dEdx one proton and pion 
      // proton2 loop
      for(Int_t it3=0;it3<ntTpc;it3++){
	p2fromLambda = false;
	itp2=it3;
	HF1( 1, 2 );  // status 2
	if(it3==it2 || it3==it1) continue;	
	HF1( 1, 3 );  // for debug // status 3
	if(event.isK18[it3]==1 || event.isKurama[it3]==1 || event.isAccidental[it3]==1 ) continue;
	HF1( 1, 4 );  // status 4
	Bool_t p_like2 = false;
	if(!FlagPidNsigma){
	  if((event.pid[it3]&4)!=4) continue; //select proton like	  
	  if(event.charge[it3]==1) p_like2 = true;	  
	} else {
	  if(event.nsigma_proton[it3]>5 || event.nsigma_proton[it3]<-2) continue;
	  //	  if(event.nsigma_pion[it3]<3) continue;	  
	  if(event.charge[it3]==1) p_like2 = true;
	}
	if(FlagDebug) std::cout << "debug: selected p2 " << std::endl;
	HF1( 1054, event.nhtrack[it3] );
	HF1( 1, 5 );  // status 5	
	//	if(event.charge[it3]==1) p_like2 = true;
	// dEdx with two protons and one pion measured
	HF2(160, event.mom0[it3]*event.charge[it3], event.dEdx[it3]);
	if(event.charge[it3]<0){
	  HF2(161, -event.mom0[it3]*event.charge[it3], event.dEdx[it3]);
	} else {
	  HF2(162, event.mom0[it3]*event.charge[it3], event.dEdx[it3]);
	}
	// ------- proton2 -------	
	Double_t p2_par[5];
	p2_par[0] = event.helix_cx[it3];
	p2_par[1] = event.helix_cy[it3];
	p2_par[2] = event.helix_z0[it3];
	p2_par[3] = event.helix_r[it3];
	p2_par[4] = event.helix_dz[it3];
	Int_t p2_nh = event.helix_t[it3].size();
	Double_t p2_theta_min = event.helix_t[it3][0] - vtx_scan_range/p2_par[3];
	Double_t p2_theta_max = TMath::Min(event.helix_t[it3][0] + vtx_scan_rangeInsideL/p2_par[3], event.helix_t[it3][p2_nh-1]);
	TVector3 p2_start = TVector3(event.calpos_x[it3][0], event.calpos_y[it3][0], event.calpos_z[it3][0]);
	TVector3 p2_end = TVector3(event.calpos_x[it3][p2_nh-1], event.calpos_y[it3][p2_nh-1], event.calpos_z[it3][p2_nh-1]);
	if(!p_like2){
	  p2_theta_min = event.helix_t[it3][p2_nh-1] - vtx_scan_range/p2_par[3];
	  p2_theta_max = TMath::Min(event.helix_t[it3][p2_nh-1] + vtx_scan_rangeInsideL/p2_par[3], event.helix_t[it3][0]);
	  p2_start = TVector3(event.calpos_x[it3][p2_nh-1], event.calpos_y[it3][p2_nh-1], event.calpos_z[it3][p2_nh-1]);
	  p2_end = TVector3(event.calpos_x[it3][0], event.calpos_y[it3][0], event.calpos_z[it3][0]);
	}
	TVector3 p2_mom; Double_t lp2_dist;
	TVector3 lp2_vert
	  = Kinematics::LambdaPVertex(dMagneticField, p2_par, p2_theta_min, p2_theta_max,
				      lambda_vert, lambda_mom, p2_mom, lp2_dist);	
	// compare the invariant masses of Lp1 and Lp2
	Double_t p2pi_dist = 10000.;
	TVector3 pi_mom2; TVector3 lambda2_mom; TVector3 p2_mom_lambda;
	TVector3 lambda2_vert
	  = Kinematics::LambdaVertex(dMagneticField, p2_par, pi_par, p2_theta_min, p2_theta_max, pi_theta_min, pi_theta_max,
				     p2_mom_lambda, pi_mom2, lambda2_mom, p2pi_dist);
	Double_t angp1p2 = p2_mom.Angle(p_mom)*TMath::RadToDeg();
	Double_t angp1p2l = p2_mom_lambda.Angle(p_mom)*TMath::RadToDeg();	    
	if( FlagDebug ) std::cout << "Before checking p2 froma Lambda or not" << std::endl;
	if( FlagDebug ) std::cout << "angle b/w p1p2 " << angp1p2 << ", angle b/w p1p2L " << angp1p2l << std::endl;
	if( FlagDebug ) std::cout << "Mag of momenta,  p1:" << p_mom.Mag() << "  p2:" << p2_mom.Mag() << std::endl;
	if( FlagDebug ) std::cout << "track id p1:" << it1 << " p2:" << it3 << std::endl;
	
	if(!TMath::IsNaN(p2pi_dist)){
	  if(FlagDebug) std::cout << "debug: NOT IsNaN(p2pi_dist) " << std::endl;
	  HF1( 1, 7 );// for debug // status 7
	  if(!pim_like) pi_mom2 = -1.*pi_mom2;
	  if(!p_like2) p2_mom_lambda = -1.*p2_mom_lambda;
	  lambda2_mom = pi_mom2 + p2_mom_lambda;
	  TLorentzVector Lpi_2(pi_mom2, TMath::Sqrt(pi_mom2.Mag()*pi_mom2.Mag() + PionMass*PionMass)); // the same pion
	  TLorentzVector Lp2l(p2_mom_lambda, TMath::Sqrt(p2_mom_lambda.Mag()*p2_mom_lambda.Mag() + ProtonMass*ProtonMass));
	  TLorentzVector Llambda2 = Lp2l + Lpi_2;
	  Double_t pi_2_vertex_dist; Double_t p2_vertex_dist;
	  if( TMath::Abs(lambda2_vert.x())<250. && TMath::Abs(lambda2_vert.z())<250. && TMath::Abs(lambda2_vert.y())<250.
	      && p2pi_dist < ppi_distcut // p2 Pi closest distance cut
	      && Kinematics::HelixDirection(lambda2_vert, p2_start, p2_end, p2_vertex_dist)
	      && Kinematics::HelixDirection(lambda2_vert, pi_start, pi_end, pi_2_vertex_dist)
	      && pi_2_vertex_dist < pi_vtx_distcut
	      && p2_vertex_dist < p_vtx_distcut
	      && TMath::Abs(Llambda2.M() - LambdaMass) < lambda_masscut
	      && TMath::Abs(Llambda2.M() - LambdaMass) < TMath::Abs(Llambda.M() - LambdaMass) ){
	    if(FlagDebug) std::cout << "debug: p2 from Lambda " << std::endl;
	    p2fromLambda = true;
	    Bool_t smallanglep1p2 = false;
	    angp1p2 = p2_mom.Angle(p_mom)*TMath::RadToDeg();
	    angp1p2l = p2_mom_lambda.Angle(p_mom)*TMath::RadToDeg();	    	    
	    if( FlagDebug ) std::cout << "anlge b/w p1p2 " << angp1p2 << ", angle b/w p1p2L " << angp1p2l << std::endl;
	    if( FlagDebug ) std::cout << "Mag of momenta,  p1:" << p_mom.Mag() << "  p2:" << p2_mom.Mag() << std::endl;
	    if( FlagDebug ) std::cout << "track id p1:" << it1 << " p2:" << it3 << std::endl;
	    if( FlagDebug ) std::cout << "p2_mom_lambda.Angle(p_mom)*TMath::RadToDeg()<0.01 : "
				      << smallanglep1p2 << std::endl;
	    if(p2_mom_lambda.Angle(p_mom)*TMath::RadToDeg()<0.01) smallanglep1p2 = true;
	    if(FlagDebug&&p2fromLambda&&smallanglep1p2){
	      std::cout << " p2 Misidentifcation. p1 track can be segemented into two. skipped this proton. " << std::endl;
	      std::cout << " Please push [Enter] " << std::endl;
	      // getchar();
	      //	      continue; // temporary
	    }
	    // replace p1 and p2
	    {
	      Int_t itp1temp=itp2;
	      Int_t itp2temp=itp1;
	      itp1=itp1temp;
	      itp2=itp2temp;
	      // Mom
	      if( FlagDebug ) std::cout << "Mom mag p_mom:" << p_mom.Mag()
					<< " p2_mom:" << p2_mom.Mag() << std::endl;	      	      
	      TVector3 p2_mom_temp = p_mom;
	      TVector3 p_mom_temp = p2_mom_lambda;
	      if( FlagDebug ) std::cout << "Mom mag p_mom_temp:" << p_mom_temp.Mag()
					<< " p2_mom_temp:" << p2_mom_temp.Mag() << std::endl;
	      p_mom = p_mom_temp;
	      p2_mom = p2_mom_temp;	      
	      std::cout << "Mom mag p_mom:" << p_mom.Mag()
			<< " p2_mom:" << p2_mom.Mag() << std::endl;	      
	      pi_mom = pi_mom2;
	      lambda_mom = lambda2_mom;
	      // LorentzVector
	      TLorentzVector Lp_temp = Lp2l;
	      Lp = Lp_temp;
	      Llambda=Llambda2;	    
	      // vertex
	      Double_t p_vertex_dist_temp = p2_vertex_dist;
	      Double_t p2_vertex_dist_temp = p_vertex_dist;
	      p_vertex_dist = p_vertex_dist_temp;
	      p2_vertex_dist = p2_vertex_dist_temp;	      
	      pi_vertex_dist = pi_2_vertex_dist;
	      lambda_vert = lambda2_vert;
	      // distance 
	      ppi_dist = p2pi_dist;
	      // helix parameter
	      Double_t p_par_temp[5];
	      Double_t p2_par_temp[5];	      
	      std::copy(std::begin(p2_par),std::end(p2_par),std::begin(p_par_temp));
	      std::copy(std::begin(p_par),std::end(p_par),std::begin(p2_par_temp));
	      
	      std::copy(std::begin(p_par_temp),std::end(p_par_temp),std::begin(p_par));	      
	      std::copy(std::begin(p2_par_temp),std::end(p2_par_temp),std::begin(p2_par));
	    	    
	      TVector3 p_start_temp = p2_start;
	      TVector3 p_end_temp = p2_end;
	      TVector3 p2_start_temp = p_start;
	      TVector3 p2_end_temp = p_end;	      
	      p_start = p_start_temp;
	      p_end = p_end_temp;	      
	      p2_start = p2_start_temp;
	      p2_end = p2_end_temp;

	      Double_t p2_theta_min_temp = p_theta_min;
	      Double_t p2_theta_max_temp = p_theta_max;
	      p2_theta_min = p2_theta_min_temp;
	      p2_theta_max = p2_theta_max_temp;
	    }
	    {
	      Double_t lp_dist;
	      if(FlagDebug) std::cout << "debug: " << __FILE__ << " line " << __LINE__ << std::endl;	
	      TVector3 lp_vert
		= Kinematics::LambdaPVertex(dMagneticField, p2_par, p2_theta_min, p2_theta_max,
					    lambda_vert, lambda_mom, p2_mom, lp_dist);
	      if(FlagDebug) std::cout << "debug: " << __FILE__ << " line " << __LINE__ << std::endl;
	      lp2_dist = lp_dist;
	      lp2_vert = lp_vert;
	    }
	  }
	}
	if(FlagDebug&&!p2fromLambda) std::cout << "debug: p2 not from Lambda" << std::endl;
	if(TMath::IsNaN(lp2_dist)) continue;
	if(FlagDebug) std::cout << "debug: NOT IsNaN(lp2_dist)" << std::endl;
	if(lp2_dist>lp2_distcut) continue;
	//	if(FlagDebug) std::cout << "debug: lp2_dist<lp2_distcut" << std::endl;
	//	if(FlagDebug) std::cout << "debug: lp2_dist=" << lp2_dist << std::endl;	
	if( TMath::Abs(lp2_vert.x()) > 250.
	    || TMath::Abs(lp2_vert.z()) > 250.
	    || TMath::Abs(lp2_vert.y()) > 250. ) continue;
	if(FlagDebug) std::cout << "debug: vertex cut applied" << std::endl;
	if(lp2_vert.z() < -200.) continue;
	if(lp2_vert.z() > 200.) continue;
	if(FlagDebug) std::cout << "debug: vertex Z cut applied" << std::endl;
	TLorentzVector Lp2(p2_mom, TMath::Sqrt(p2_mom.Mag()*p2_mom.Mag() + ProtonMass*ProtonMass));
	TLorentzVector Llambda_fixedmass(lambda_mom, TMath::Sqrt(lambda_mom.Mag()*lambda_mom.Mag() + LambdaMass*LambdaMass));
	TLorentzVector Llp2 = Llambda_fixedmass + Lp2;
	TVector3 lp2_mom = TVector3(Llp2.Px(), Llp2.Py(), Llp2.Pz());
	Double_t p2_vertex_dist = 1000;
	Double_t lp2_opening_angle;	
	if(!Kinematics::HelixDirection(lp2_vert, p2_start, p2_end, p2_vertex_dist)) continue;
	if(FlagDebug) std::cout << "debug: vertex Z cut applied" << std::endl;
	if(p2_vertex_dist > p2_vtx_distcut) continue;
	if(FlagDebug) std::cout << "debug p2_vertex_dist : " << p2_vertex_dist << std::endl;
	std::cout << "debug ProtonMass:" << ProtonMass << " LambdaMass:" << LambdaMass << std::endl;
	std::cout << "debug ProtonMom:" << p2_mom.Mag() << " LambdaMom:" << lambda_mom.Mag() << std::endl;
	std::cout << "debug Lp2Mass:" << Llp2.M() << " Lp2Mom:" << Llp2.P() << std::endl;
	HF1( 1, 15 );// for debug // status 15
	Double_t thetap2 = p2_mom.Theta()*TMath::RadToDeg();
	
	HF1( 1302, lp2_mom.Mag() );
	
	lp2_opening_angle = lambda_mom.Angle(p2_mom)*TMath::RadToDeg();
	lp2_l_container.push_back(l_candidates);
	
	lp2_p_container.push_back(itp1);
	lp2_pi_container.push_back(itpi);
	lp2_p2_container.push_back(itp2);
	
	HF2(130, event.mom0[itp1]*event.charge[itp1], event.dEdx[itp1]);
	HF2(130, event.mom0[itpi]*event.charge[itpi], event.dEdx[itpi]);
	HF2(130, event.mom0[itp2]*event.charge[itp2], event.dEdx[itp2]);		
	if(event.charge[itp1]<0){
	  HF2(131, -event.mom0[itp1]*event.charge[itp1], event.dEdx[itp1]);
	} else {
	  HF2(132, event.mom0[itp1]*event.charge[itp1], event.dEdx[itp1]);
	}
	if(event.charge[itpi]<0){
	  HF2(131, -event.mom0[itpi]*event.charge[itpi], event.dEdx[itpi]);
	} else {
	  HF2(132, event.mom0[itpi]*event.charge[itpi], event.dEdx[itpi]);
	}
	if(event.charge[itp2]<0){
	  HF2(131, -event.mom0[itp2]*event.charge[itp2], event.dEdx[itp2]);
	} else {
	  HF2(132, event.mom0[itp2]*event.charge[itp2], event.dEdx[itp2]);
	}
	
	lp2mass_container.push_back(Llp2.M());
	lambda_lp2_mass_container.push_back(Llambda.M());
	lp2_p_mom_container.push_back(p_mom);
	lp2_pi_mom_container.push_back(pi_mom);
	lp2_p2_mom_container.push_back(p2_mom);
	lp2_l_mom_container.push_back(lambda_mom);
	lp2_mom_container.push_back(lp2_mom);
	lp2_closedist.push_back(lp2_dist);
	lp2_vert_container.push_back(lp2_vert);
	lp2angle_container.push_back(lp2_opening_angle);
	Double_t thetapi = pi_mom.Theta()*TMath::RadToDeg();
	Double_t thetaL = lambda_mom.Theta()*TMath::RadToDeg();	
	if(FlagKuramaProton){	    
	  // for(int itKurama=0; itKurama<event.ntKurama; itKurama++){
	  //   Double_t thetap0 = momPKurama[itKurama].Theta()*TMath::RadToDeg();
	  //   HF2( 2615, thetap0, momPKurama[itKurama].Mag() );	      
	  //   Double_t anglep0p1 = p_mom.Angle(momPKurama[itKurama])*TMath::RadToDeg();
	  //   HF1( 1170, anglep0p1 );
	  //   Double_t anglep0pi = pi_mom.Angle(momPKurama[itKurama])*TMath::RadToDeg();
	  //   HF1( 1171, anglep0pi );
	  //   Double_t anglep0p2 = p2_mom.Angle(momPKurama[itKurama])*TMath::RadToDeg();
	  //   HF1( 1172, anglep0p2 );
	  //   Double_t anglep0L = lambda_mom.Angle(momPKurama[itKurama])*TMath::RadToDeg();
	  //   HF1( 1173, anglep0p2 );
	  //   HF2( 2451, thetap0, thetap0);
	  //   HF2( 2452, thetap0, thetap1);
	  //   HF2( 2453, thetap0, thetapi);
	  //   HF2( 2454, thetap0, thetaL);
	  //   HF2( 2455, thetap0, thetap2);
	  //   // mom corr
	  //   HF2( 2251, momPKurama[itKurama].Mag(), momPKurama[itKurama].Mag());
	  //   HF2( 2252, momPKurama[itKurama].Mag(), p_mom.Mag());
	  //   HF2( 2253, momPKurama[itKurama].Mag(), pi_mom.Mag());
	  //   HF2( 2254, momPKurama[itKurama].Mag(), lambda_mom.Mag());
	  //   HF2( 2255, momPKurama[itKurama].Mag(), p2_mom.Mag());
	  // }
	}
	{
	  Double_t anglep1pi = p_mom.Angle(pi_mom)*TMath::RadToDeg();
	  Double_t anglep1p2 = p_mom.Angle(p2_mom)*TMath::RadToDeg();
	  Double_t anglep1L = p_mom.Angle(lambda_mom)*TMath::RadToDeg();
	  Double_t anglepip2 = pi_mom.Angle(p2_mom)*TMath::RadToDeg();
	  Double_t anglepiL = pi_mom.Angle(lambda_mom)*TMath::RadToDeg();
	  Double_t anglep2L = p2_mom.Angle(lambda_mom)*TMath::RadToDeg();
	  Double_t diffP1P2 = TMath::Abs(p_mom.Mag()-p2_mom.Mag());
	  if( FlagDebug ) std::cout << "Mag of momenta,  p1:" << p_mom.Mag() << "  p2:" << p2_mom.Mag() << std::endl;
	  if( FlagDebug ) std::cout << "track id p1:" << itp1 << " p2:" << itp2 << std::endl;	    
	  if( FlagDebug ) std::cout << " |p1-p2|=" << diffP1P2 << std::endl;
	  if( FlagDebug ) std::cout << " anglep1p2=" << anglep1p2 << std::endl;	    
	  if( FlagDebug ) std::cout << "p2fromLambda:" << p2fromLambda << std::endl;
	  if( FlagDebug ) std::cout << "p1track:" << itp1 << " p2track:" << itp2 << std::endl;
	  if( FlagDebug ) std::cout << "lp2OpAngle:" << lp2_opening_angle << std::endl;
	  if( FlagDebug ) std::cout << "anglep1L:" << anglep1L << std::endl;
	}
	event.lp2flag = true;
	lp2_candidates++;
	if(FlagDebug) std::cout << "debug: LambdaP hist filled " << std::endl;
      } //it 3
      Double_t ppi_opening_angle = p_mom.Angle(pi_mom)*TMath::RadToDeg();;	  
      ppi_closedist.push_back(ppi_dist);
      ppiangle_container.push_back(ppi_opening_angle);
      l_p_container.push_back(itp1);
      l_pi_container.push_back(itpi);
      lambda_mass_container.push_back(Llambda.M());
      l_p_mom_container.push_back(p_mom);
      l_pi_mom_container.push_back(pi_mom);
      l_mom_container.push_back(lambda_mom);	  
      l_vert_container.push_back(lambda_vert);
      l_candidates++;
    } //it2
  } //it1

// Kurama proton loop 
  // for(Int_t it0=0; it0<ntTpc; it0++){
  //   Bool_t FlagKuramaProton = false;
  //   Bool_t p0_like = false;
  //   TVector3 momPKurama[event.ntKurama]; 
  //   if(event.isKurama[it0]==1 && (event.pid[it0]&4)==4 && event.charge[it0]==1 ){ // proton like
  //     FlagKuramaProton=true;
  //     p0_like = true;
  //     for(int itKurama=0; itKurama<event.ntKurama; itKurama++){
  // 	Double_t pmag = event.pTPCKurama[itKurama];
  // 	Double_t theta = event.thetaTPCKurama[itKurama];
  // 	Double_t u = event.utgtTPCKurama[itKurama];
  // 	Double_t v = event.vtgtTPCKurama[itKurama];
  // 	momPKurama[itKurama].SetXYZ(pmag*TMath::Cos(theta)*u,
  // 				    pmag*TMath::Cos(theta)*v,
  // 				    pmag*TMath::Cos(theta)*1.0);
  // 	kurama_p_mom_container[itKurama].push_back(momPKurama[itKurama]);
  //     }
  //     Double_t p0_par[5];
  //     p0_par[0] = event.helix_cx[it0];
  //     p0_par[1] = event.helix_cy[it0];
  //     p0_par[2] = event.helix_z0[it0];
  //     p0_par[3] = event.helix_r[it0];
  //     p0_par[4] = event.helix_dz[it0];
  //     Int_t p0_nh = event.helix_t[it0].size();
  //     Double_t p0_theta_min = event.helix_t[it0][0] - vtx_scan_range/p0_par[3];
  //     Double_t p0_theta_max = TMath::Min(event.helix_t[it0][0] + vtx_scan_rangeInsideL/p0_par[3], event.helix_t[it0][p0_nh-1]);
  //     TVector3 p0_start = TVector3(event.calpos_x[it0][0],
  // 				   event.calpos_y[it0][0],
  // 				   event.calpos_z[it0][0]);
  //     TVector3 p0_end = TVector3(event.calpos_x[it0][p0_nh-1],
  // 				 event.calpos_y[it0][p0_nh-1],
  // 				 event.calpos_z[it0][p0_nh-1]);
  //     if(!p0_like){
  // 	p0_theta_min = event.helix_t[it0][p0_nh-1] - vtx_scan_range/p0_par[3];
  // 	p0_theta_max = TMath::Min(event.helix_t[it0][p0_nh-1] + vtx_scan_rangeInsideL/p0_par[3], event.helix_t[it0][0]);
  // 	p0_start = TVector3(event.calpos_x[it0][p0_nh-1],
  // 			    event.calpos_y[it0][p0_nh-1],
  // 			    event.calpos_z[it0][p0_nh-1]);
  // 	p0_end = TVector3(event.calpos_x[it0][0],
  // 			  event.calpos_y[it0][0],
  // 			  event.calpos_z[it0][0]);
  //     }    
  //     // pion from Forward Lambda
  //     for(Int_t itPi=0;itPi<ntTpc;itPi++){
  // 	if(it0==itPi) continue;
  // 	if(event.isK18[itPi]==1 || event.isAccidental[itPi]==1 || event.isKurama[itPi]==1 ) continue;
  // 	if((event.pid[itPi]&1)!=1) continue; //select pi like
  // 	Bool_t pimForwardL_like = false;
  // 	if(event.charge[itPi]==-1) pimForwardL_like = true;
  // 	Double_t piForward_par[5];
  // 	piForward_par[0] = event.helix_cx[itPi];
  // 	piForward_par[1] = event.helix_cy[itPi];
  // 	piForward_par[2] = event.helix_z0[itPi];
  // 	piForward_par[3] = event.helix_r[itPi];
  // 	piForward_par[4] = event.helix_dz[itPi];
  // 	Int_t piForward_nh = event.helix_t[itPi].size();
  // 	Double_t piForward_theta_min = TMath::Max(event.helix_t[itPi][0] - vtx_scan_rangeInsideL/piForward_par[3],
  // 						  event.helix_t[itPi][piForward_nh-1]);
  // 	Double_t piForward_theta_max = event.helix_t[itPi][0] + vtx_scan_range/piForward_par[3];
  // 	TVector3 piForward_start = TVector3(event.calpos_x[itPi][0],
  // 					    event.calpos_y[itPi][0],
  // 					    event.calpos_z[itPi][0]);
  // 	TVector3 piForward_end = TVector3(event.calpos_x[itPi][piForward_nh-1],
  // 					  event.calpos_y[itPi][piForward_nh-1],
  // 					  event.calpos_z[itPi][piForward_nh-1]);
  // 	if(!pimForwardL_like){
  // 	  piForward_theta_min = TMath::Max(event.helix_t[itPi][piForward_nh-1] - vtx_scan_rangeInsideL/piForward_par[3], event.helix_t[itPi][0]);
  // 	  piForward_theta_max = event.helix_t[itPi][piForward_nh-1] + vtx_scan_range/piForward_par[3];
  // 	  piForward_start = TVector3(event.calpos_x[itPi][piForward_nh-1],
  // 				     event.calpos_y[itPi][piForward_nh-1],
  // 				     event.calpos_z[itPi][piForward_nh-1]);
  // 	  piForward_end = TVector3(event.calpos_x[itPi][0],
  // 				   event.calpos_y[itPi][0],
  // 				   event.calpos_z[itPi][0]);
  // 	}
  // 	Double_t ppiForward_dist = 10000.;
  // 	TVector3 p0_mom; TVector3 piForward_mom; TVector3 lambdaForward_mom;
  // 	TVector3 lambdaForward_vert = Kinematics::LambdaVertex(dMagneticField,
  // 							       p0_par,
  // 							       piForward_par,
  // 							       p0_theta_min,
  // 							       p0_theta_max,
  // 							       piForward_theta_min,
  // 							       piForward_theta_max,
  // 							       p0_mom, piForward_mom,
  // 							       lambdaForward_mom,
  // 							       ppiForward_dist);
  // 	if(TMath::IsNaN(ppiForward_dist)) continue;
  // 	if(!pimForwardL_like) piForward_mom = -1.*piForward_mom;
  // 	if(!p0_like) piForward_mom = -1.*p0_mom;
  // 	lambdaForward_mom = piForward_mom + p0_mom;

  // 	TLorentzVector LpiForward(piForward_mom, TMath::Sqrt(piForward_mom.Mag()*piForward_mom.Mag() + PionMass*PionMass));
  // 	TLorentzVector Lp0(p0_mom, TMath::Sqrt(p0_mom.Mag()*p0_mom.Mag() + ProtonMass*ProtonMass));
  // 	TLorentzVector LlambdaForward = Lp0 + LpiForward;
  // 	if(TMath::Abs(lambdaForward_vert.x()) > 250.
  // 	   || TMath::Abs(lambdaForward_vert.z()) > 250.
  // 	   || TMath::Abs(lambdaForward_vert.y()) > 250.) continue; //Vertex Cut
  // 	if(ppiForward_dist > ppi_distcut) continue;
  // 	Double_t piForward_vertex_dist; Double_t p0_vertex_dist;
  // 	if(!Kinematics::HelixDirection(lambdaForward_vert, p0_start, p0_end, p0_vertex_dist)
  // 	   || !Kinematics::HelixDirection(lambdaForward_vert, piForward_start, piForward_end, piForward_vertex_dist)) continue;
  // 	if(piForward_vertex_dist > pi_vtx_distcut) continue;
  // 	if(p0_vertex_dist > p_vtx_distcut) continue;
  // 	if(TMath::Abs(LlambdaForward.M() - LambdaMass) > lambda_masscut) continue;
  // 	event.forwardlflag = true;
  //     }
  //   }  
    
  if(FlagDebug) std::cout << "Lambda candidates: " << l_candidates << " and LambdaP candidates: " << lp2_candidates << std::endl;
  HF1(4,l_candidates);
  HF1(5,lp2_candidates);

#if Debug
  std::cout << "Debug: select best Lambda and LambdaP" << std::endl; 
#endif
  
#if 1 //Select the best L mass combination
  Int_t best = -1; Double_t prev_massdiff = 9999.;
  for(Int_t candi=0;candi<l_candidates;candi++){
    Double_t diff = TMath::Abs(lambda_mass_container[candi] - LambdaMass);
    if(prev_massdiff > diff){
      prev_massdiff = diff;
      best = candi;
      std::cout << "best lambda: " << best << std::endl;
    }
  }
#endif
  // #if 1
  //   Int_t best_lp2 = -1; Double_t prev_distdiff_lp2 = 9999.;
  //   for(Int_t candi=0;candi<lp2_candidates;candi++){    
  //     if(lp2_l_container.at(candi) != best) continue;
  //     Double_t diff = TMath::Abs(lp2_closedist.at(candi));
  //     //Double_t diff = TMath::Abs(lambda_lp2_mass_container[candi] - LambdaMass);
  //     //    std::cout << "invariant mass Lp2 : "
  //     //	      << lp2mass_container[candi]
  //     //	      << std::endl;
  //     if(prev_distdiff_lp2 > diff){
  //       prev_distdiff_lp2 = diff;
  //       best_lp2 = candi;
  //       std::cout << "best lambda proton: " << best_lp2 << std::endl;
  //     }
  //   }
  // #endif
  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;  
  //  if(best==-1) return true;  

  Int_t pcount = std::pow(10,6);
  for(int i=0; i<6; i++){
    if(counter[i]>9) std::cout << "over limit counter" << i << ": " << counter[i] << std::endl;    
    std::cout << i << "loop, counter:" << counter[i] << std::endl;
    pcount = pcount+counter[i]*std::pow(10,i);
    std::cout << "total count: " << pcount << endl;    
  }
  std::cout << "total count: " << pcount << std::endl;

  if(best!=-1){
  event.lmass = lambda_mass_container[best];
  event.ldecayvtx_x = l_vert_container[best].x();
  event.ldecayvtx_y = l_vert_container[best].y();
  event.ldecayvtx_z = l_vert_container[best].z();  
  event.lmom_x = l_mom_container[best].x();
  event.lmom_y = l_mom_container[best].y();
  event.lmom_z = l_mom_container[best].z();
  event.ppi_dist = ppi_closedist[best];
  event.ppiangle = ppiangle_container[best];  
  
  event.ldecays_id.push_back(l_p_container[best]);
  event.ldecays_mom.push_back(l_p_mom_container[best].Mag());
  event.ldecays_mom_x.push_back(l_p_mom_container[best].x());
  event.ldecays_mom_y.push_back(l_p_mom_container[best].y());
  event.ldecays_mom_z.push_back(l_p_mom_container[best].z());
  //  event.ldecays_theta.push_back(l_p_mom_container[best].Theta()*TMath::RadToDeg());
  
  event.ldecays_id.push_back(l_pi_container[best]); 
  event.ldecays_mom.push_back(l_pi_mom_container[best].Mag()); 
  event.ldecays_mom_x.push_back(l_pi_mom_container[best].x()); 
  event.ldecays_mom_y.push_back(l_pi_mom_container[best].y()); 
  event.ldecays_mom_z.push_back(l_pi_mom_container[best].z()); 
  //  event.ldecays_theta.push_back(l_pi_mom_container[best].Theta()*TMath::RadToDeg()); 
  
  Double_t thetaP1 = l_p_mom_container[best].Theta()*TMath::RadToDeg();
  Double_t thetaPi = l_pi_mom_container[best].Theta()*TMath::RadToDeg();
  Double_t thetaL = l_mom_container[best].Theta()*TMath::RadToDeg();
  Double_t angleP1Pi = l_p_mom_container[best].Angle(l_pi_mom_container[best])*TMath::RadToDeg();
  Double_t angleP1L = l_p_mom_container[best].Angle(l_mom_container[best])*TMath::RadToDeg();
  Double_t anglePiL = l_pi_mom_container[best].Angle(l_mom_container[best])*TMath::RadToDeg();    
  Double_t momP1 = l_p_mom_container[best].Mag();
  Double_t momPi = l_pi_mom_container[best].Mag();
  Double_t momL = l_mom_container[best].Mag();
  
  //HF1( 1106, thetaP );
  HF1( 1107, thetaP1 );
  //  HF1( 1111, thetaPi );
  HF1( 1113, thetaPi );  
  //  HF1( 1116, thetaL );
  HF1( 1119, thetaL );    
  // Fill histograms with lflag
  HF1( 211, event.ldecayvtx_x);
  HF1( 212, event.ldecayvtx_y);
  HF1( 213, event.ldecayvtx_z);

  HF2( 2601, thetaP1, momP1 );    
  HF2( 2602, thetaPi, momPi );
  HF2( 2603, thetaL, momL );        

  HF1( 1163, angleP1Pi );
  HF1( 1164, angleP1L );
  HF1( 1165, anglePiL );
  
  HF2( 2411, thetaP1, thetaP1 );
  HF2( 2412, thetaP1, thetaPi );
  HF2( 2413, thetaP1, thetaL );  
  HF2( 2421, thetaPi, thetaPi );
  HF2( 2422, thetaPi, thetaL );
  HF2( 2431, thetaL, thetaL );
  
  HF2( 2211, momP1, momP1 );
  HF2( 2212, momP1, momPi );
  HF2( 2213, momP1, momL );  
  HF2( 2221, momPi, momPi );
  HF2( 2222, momPi, momL );
  HF2( 2231, momL, momL );

  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
	
  for(int i=0; i<event.nKK; ++i){
    //    if(momTransMag.at(i)<0.025) continue;
    if(momK18[i].Mag()<0.025) continue;    
    if(momKurama[i].Mag()<0.025) continue;
    Double_t bek = event.BEkaonTPC[i];
    Double_t theta = event.thetaTPC[i];
    HF1( 1204, bek );
    HF1( 11204, bek );
    HF1( 1401, event.ldecays_mom[0] );
    HF1( 1411, event.ldecays_mom[1] );
    HF1( 1431, l_mom_container[best].Mag() );
    HF1( 62, momTransMag.at(i) );
    HF2( 2001, event.ppiangle, l_p_mom_container[best].Mag() );
    HF2( 2002, event.ppiangle, l_pi_mom_container[best].Mag() );
    HF2( 2003, event.ppiangle, l_mom_container[best].Mag() );
    HF2( 2011, event.ppiangle, bek );
    HF2( 2012, event.ppiangle, event.lmass );
    
    if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
    
    if( theta<5. ){
      HF1( 1205, bek );
    }
    if( theta<10. ){
      HF1( 1206, bek );
    }
    HF1( 1211, event.lmass );
    if ( bek<-0.1 ){
      //  HF1( 1107, thetaP );
      HF1( 1108, thetaP1 );      
      //      HF1( 1112, thetaPi );
      HF1( 1114, thetaPi );      
      //      HF1( 1117, thetaL );
      HF1( 1120, thetaL );                  
      HF1( 1212, event.lmass );
      HF1( 1402, event.ldecays_mom[0] );
      HF1( 1412, event.ldecays_mom[1] );
      HF1( 1432, l_mom_container[best].Mag() );
    } else if ( bek<0. ){
      HF1( 1109, thetaP1 );      
      HF1( 1115, thetaPi );
      HF1( 1121, thetaL );                
      HF1( 1213, event.lmass );
      HF1( 1403, event.ldecays_mom[0] );
      HF1( 1413, event.ldecays_mom[1] );
      HF1( 1433, l_mom_container[best].Mag() );
    } else if ( bek<0.1 ){
      HF1( 1110, thetaP1 );
      HF1( 1116, thetaPi );
      HF1( 1122, thetaL );      
      HF1( 1214, event.lmass );
      HF1( 1404, event.ldecays_mom[0] );
      HF1( 1414, event.ldecays_mom[1] );
      HF1( 1434, l_mom_container[best].Mag() );                
    } else if ( bek<0.2 ){
      HF1( 1111, thetaP1 );
      HF1( 1117, thetaPi );
      HF1( 1123, thetaL );          
      HF1( 1215, event.lmass );
      HF1( 1405, event.ldecays_mom[0] );
      HF1( 1415, event.ldecays_mom[1] );
      HF1( 1435, l_mom_container[best].Mag() );                
    } else if ( bek<0.3 ){
      HF1( 1112, thetaP1 );
      HF1( 1118, thetaPi );
      HF1( 1124, thetaL );          
      HF1( 1216, event.lmass );
      HF1( 1406, event.ldecays_mom[0] );
      HF1( 1416, event.ldecays_mom[1] );
      HF1( 1436, l_mom_container[best].Mag() );                
    }
  }
  }
  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
  //  if(best_lp2 != -1){
  if(best!=-1 && lp2_candidates == 1){
    Int_t best_lp2 = 0;
    //    for(int best_lp2=0; best_lp2<lp2_candidates; best_lp2++){
    event.lmass_lp2 = lambda_lp2_mass_container[best_lp2];
    event.lp2mass = lp2mass_container[best_lp2];
    event.lp2decayvtx_x = lp2_vert_container[best_lp2].x();
    event.lp2decayvtx_y = lp2_vert_container[best_lp2].y();
    event.lp2decayvtx_z = lp2_vert_container[best_lp2].z();
    event.lp2mom_x = lp2_mom_container[best_lp2].x();
    event.lp2mom_y = lp2_mom_container[best_lp2].y();
    event.lp2mom_z = lp2_mom_container[best_lp2].z();

    //    event.lmom_lp2 = lp2_l_mom_container[best_lp2].Mag();
    event.lmom_lp2_x = lp2_l_mom_container[best_lp2].x();
    event.lmom_lp2_y = lp2_l_mom_container[best_lp2].y();
    event.lmom_lp2_z = lp2_l_mom_container[best_lp2].z();
    event.lp2_dist = lp2_closedist[best_lp2];
    event.lp2angle = lp2angle_container[best_lp2];
    event.lp2decays_id.push_back(lp2_p_container[best_lp2]);
    event.lp2decays_mom.push_back(lp2_p_mom_container[best_lp2].Mag());
    event.lp2decays_mom_x.push_back(lp2_p_mom_container[best_lp2].x());
    event.lp2decays_mom_y.push_back(lp2_p_mom_container[best_lp2].y());
    event.lp2decays_mom_z.push_back(lp2_p_mom_container[best_lp2].z());
    event.lp2decays_id.push_back(lp2_pi_container[best_lp2]);
    event.lp2decays_mom.push_back(lp2_pi_mom_container[best_lp2].Mag());
    event.lp2decays_mom_x.push_back(lp2_pi_mom_container[best_lp2].x());
    event.lp2decays_mom_y.push_back(lp2_pi_mom_container[best_lp2].y());
    event.lp2decays_mom_z.push_back(lp2_pi_mom_container[best_lp2].z());
    event.lp2decays_id.push_back(lp2_p2_container[best_lp2]);
    event.lp2decays_mom.push_back(lp2_p2_mom_container[best_lp2].Mag());
    event.lp2decays_mom_x.push_back(lp2_p2_mom_container[best_lp2].x());
    event.lp2decays_mom_y.push_back(lp2_p2_mom_container[best_lp2].y());
    event.lp2decays_mom_z.push_back(lp2_p2_mom_container[best_lp2].z());
    // Fill histograms
    HF1( 200, event.lp2_dist);
    HF1( 201, event.lp2decayvtx_x);
    HF1( 202, event.lp2decayvtx_y);
    HF1( 203, event.lp2decayvtx_z);
    HF1( 113, event.lp2decays_mom[0] );
    HF1( 114, event.lp2decays_mom[1] );
    HF1( 115, event.lp2decays_mom[2] );
    HF1( 1441, event.lp2decays_mom[0] );
    HF1( 1451, event.lp2decays_mom[1] );
    HF1( 1461, event.lp2decays_mom[2] );
    HF1( 1471, lp2_l_mom_container[best_lp2].Mag() );
    HF1( 1101, event.lp2angle );
    HF1( 1217, event.lp2mass );
    HF1( 1231, event.lmass );
    Double_t thetaP1 = lp2_p_mom_container[best_lp2].Theta()*TMath::RadToDeg();
    Double_t thetaPi = lp2_pi_mom_container[best_lp2].Theta()*TMath::RadToDeg();
    Double_t thetaP2 = lp2_p2_mom_container[best_lp2].Theta()*TMath::RadToDeg();
    Double_t thetaL = lp2_l_mom_container[best_lp2].Theta()*TMath::RadToDeg();

    Double_t angleP1Pi = lp2_p_mom_container[best].Angle(lp2_pi_mom_container[best])*TMath::RadToDeg();
    Double_t angleP1P2 = lp2_p_mom_container[best].Angle(lp2_p2_mom_container[best])*TMath::RadToDeg();    
    Double_t angleP1L = lp2_p_mom_container[best].Angle(lp2_l_mom_container[best])*TMath::RadToDeg();
    Double_t anglePiP2 = lp2_pi_mom_container[best].Angle(lp2_p2_mom_container[best])*TMath::RadToDeg();    
    Double_t anglePiL = lp2_pi_mom_container[best].Angle(lp2_l_mom_container[best])*TMath::RadToDeg();    
    Double_t angleP2L = lp2_p2_mom_container[best].Angle(lp2_l_mom_container[best])*TMath::RadToDeg();    
    Double_t momP1 = lp2_p_mom_container[best].Mag();
    Double_t momPi = lp2_pi_mom_container[best].Mag();
    Double_t momP2 = lp2_p2_mom_container[best].Mag();        
    Double_t momL = lp2_l_mom_container[best].Mag();

    if(FlagOutput&&angleP1P2<0.01){
      std::ofstream appendfile(outputfile, std::ios::app);
      if (!appendfile.is_open()) {
	std::cerr << "Failed to open file for appending!" << std::endl;
	return 1;
      }
      std::cout << "start writing to " << outputfile << std::endl;      
      std::cout  << ievent
		 << " " << lp2_p_container[best_lp2]
		 << " " << lp2_pi_container[best_lp2]
		 << " " << lp2_p2_container[best_lp2]
		 << " " << event.lp2angle
		 << std::endl;
      appendfile << ievent
		 << " " << lp2_p_container[best_lp2]
		 << " " << lp2_pi_container[best_lp2]
		 << " " << lp2_p2_container[best_lp2]
		 << std::endl;
      appendfile.close();
    }    

  
    HF1( 1174, angleP1Pi );
    HF1( 1175, angleP1P2 );  
    HF1( 1176, angleP1L );
    HF1( 1177, anglePiP2 );  
    HF1( 1178, anglePiL );
    HF1( 1179, angleP2L );


    HF2( 2611, thetaP1, momP1 );    
    HF2( 2612, thetaPi, momPi );
    HF2( 2613, thetaP2, momP2 );
    HF2( 2614, thetaL, momL );        
    
    HF2( 2461, thetaP1, thetaP1 );
    HF2( 2462, thetaP1, thetaPi );
    HF2( 2463, thetaP1, thetaL );
    HF2( 2464, thetaP1, thetaP2 );
    HF2( 2471, thetaPi, thetaPi );
    HF2( 2472, thetaPi, thetaL );
    HF2( 2473, thetaPi, thetaP2 );
    HF2( 2481, thetaL, thetaL );
    HF2( 2482, thetaL, thetaP2 );
    HF2( 1180, thetaL, thetaP2 );    
    HF2( 2491, thetaP2, thetaP2 );
  
    HF2( 2261, momP1, momP1 );
    HF2( 2262, momP1, momPi );
    HF2( 2263, momP1, momL );
    HF2( 2264, momP1, momP2 );      
    HF2( 2271, momPi, momPi );
    HF2( 2272, momPi, momL );
    HF2( 2273, momPi, momP2 );    
    HF2( 2281, momL, momL );
    HF2( 2282, momL, momP2 );
    HF2( 2291, momP2, momP2 );        

    HF1( 1125, thetaP1 );
    HF1( 1131, thetaPi );
    HF1( 1137, thetaP2 );
    HF1( 1143, thetaL );
    for(int i=0; i<event.nKK; ++i){
      if(momK18[i].Mag()<0.025) continue;      
      if(momKurama[i].Mag()<0.025) continue;
      Double_t bek = event.BEkaonTPC[i];
      Double_t theta = event.thetaTPC[i];
      Double_t missMom = (lp2_mom_container[best_lp2] - momTrans[i])*momTrans[i]/momTransMag.at(i);
      HF2( 4001, bek, angleP1P2 );
      HF1( 72, momTransMag.at(i) );
      HF1( 1481, missMom );
      HF2( 3000, event.lp2angle, lp2_mom_container[best_lp2].Mag() );
      HF2( 3010, event.lp2angle, event.lp2mass );      
      HF2( 3020, lp2_mom_container[best_lp2].Mag(), momTransMag.at(i) );
      if( bek<-0.1 ){	
	HF1( 1102, event.lp2angle );
	HF1( 1218, event.lp2mass );
	HF1( 1232, event.lmass );	
	HF1( 1126, thetaP1 );
	HF1( 1132, thetaPi );
	HF1( 1138, thetaP2 );
	HF1( 1144, thetaL );
	HF2( 1181, thetaL, thetaP2 );    	
	HF1( 1442, event.lp2decays_mom[0] );
	HF1( 1452, event.lp2decays_mom[1] );
	HF1( 1462, event.lp2decays_mom[2] );
	HF1( 1472, lp2_l_mom_container[best_lp2].Mag() );
	HF1( 1482, missMom );	
	HF2( 3001, event.lp2angle, lp2_mom_container[best_lp2].Mag() );
	HF2( 3011, event.lp2angle, event.lp2mass );
	HF2( 3021, lp2_mom_container[best_lp2].Mag(), momTransMag.at(i) );	
      } else if ( bek<0. ){
	HF1( 1103, event.lp2angle );
	HF1( 1219, event.lp2mass );
	HF1( 1233, event.lmass );		
	HF1( 1127, thetaP1 );
	HF1( 1133, thetaPi );
	HF1( 1139, thetaP2 );
	HF1( 1145, thetaL );
	HF2( 1182, thetaL, thetaP2 );    		
	HF1( 1443, event.lp2decays_mom[0] );
	HF1( 1453, event.lp2decays_mom[1] );
	HF1( 1463, event.lp2decays_mom[2] );
	HF1( 1473, lp2_l_mom_container[best_lp2].Mag() );
	HF1( 1483, missMom );		
	HF2( 3002, event.lp2angle, lp2_mom_container[best_lp2].Mag() );
	HF2( 3012, event.lp2angle, event.lp2mass );
	HF2( 3022, lp2_mom_container[best_lp2].Mag(), momTransMag.at(i) );		
      } else if ( bek<0.1 ){
	HF1( 1104, event.lp2angle );
	HF1( 1220, event.lp2mass );
	HF1( 1234, event.lmass );			
	HF1( 1128, thetaP1 );
	HF1( 1134, thetaPi );
	HF1( 1140, thetaP2 );
	HF1( 1146, thetaL );
	HF2( 1183, thetaL, thetaP2 );    			
	HF1( 1444, event.lp2decays_mom[0] );
	HF1( 1454, event.lp2decays_mom[1] );
	HF1( 1464, event.lp2decays_mom[2] );
	HF1( 1474, lp2_l_mom_container[best_lp2].Mag() );
	HF1( 1484, missMom );			
	HF2( 3003, event.lp2angle, lp2_mom_container[best_lp2].Mag() );
	HF2( 3013, event.lp2angle, event.lp2mass );
	HF2( 3023, lp2_mom_container[best_lp2].Mag(), momTransMag.at(i) );			
      } else if ( bek<0.2 ){
	HF1( 1105, event.lp2angle );
	HF1( 1221, event.lp2mass );
	HF1( 1235, event.lmass );				
	HF1( 1129, thetaP1 );
	HF1( 1135, thetaPi );
	HF1( 1141, thetaP2 );
	HF1( 1147, thetaL );
	HF2( 1184, thetaL, thetaP2 );    				
	HF1( 1445, event.lp2decays_mom[0] );
	HF1( 1455, event.lp2decays_mom[1] );
	HF1( 1465, event.lp2decays_mom[2] );
	HF1( 1475, lp2_l_mom_container[best_lp2].Mag() );
	HF1( 1485, missMom );				
	HF2( 3004, event.lp2angle, lp2_mom_container[best_lp2].Mag() );
	HF2( 3014, event.lp2angle, event.lp2mass );
	HF2( 3024, lp2_mom_container[best_lp2].Mag(), momTransMag.at(i) );			
      } else if ( bek<0.3 ){
	HF1( 1106, event.lp2angle );
	HF1( 1222, event.lp2mass );
	HF1( 1236, event.lmass );
	HF1( 1130, thetaP1 );
	HF1( 1136, thetaPi );
	HF1( 1142, thetaP2 );
	HF1( 1148, thetaL );
	HF2( 1185, thetaL, thetaP2 );    					
	HF1( 1446, event.lp2decays_mom[0] );
	HF1( 1456, event.lp2decays_mom[1] );
	HF1( 1466, event.lp2decays_mom[2] );
	HF1( 1476, lp2_l_mom_container[best_lp2].Mag() );
	HF1( 1486, missMom );					
	HF2( 3005, event.lp2angle, lp2_mom_container[best_lp2].Mag() );
	HF2( 3015, event.lp2angle, event.lp2mass );
	HF2( 3025, lp2_mom_container[best_lp2].Mag(), momTransMag.at(i) );				
      }

      HF1( 1207, bek );
      HF1( 1250, bek );
      HF1( 1260, momTransMag.at(i) );
      HF1( 1270, lp2_mom_container[best_lp2].Mag() );
      HF1( 1280, missMom );
      
      HF2( 3040, bek, momTransMag.at(i) );
      HF2( 3050, bek, event.lp2mass );
      HF1( 1240, event.lp2mass );      
      if( theta<5. ){
	HF1( 1310, bek );
	HF1( 1208, bek );
	if( thetaL<20 ){
	  HF1( 1311, bek );
	} else if ( thetaL< 100 ){
	  HF1( 1312, bek );
	} else if ( thetaL< 180 ){
	  HF1( 1313, bek );
	}
	if ( thetaP2<20 ){
	  HF1( 1314, bek );
	} else if ( thetaP2<90 ){
	  HF1( 1315, bek );
	} else if ( thetaP2< 160 ){
	  HF1( 1316, bek );
	} else if ( thetaP2<180 ){
	  HF1( 1317, bek );
	}
      }
      if( theta<10. ) HF1( 1209, bek );
      if( event.lp2angle<20.0 ){
	HF1( 1241, event.lp2mass );      	
	HF1( 1251, bek );
	HF1( 1261, momTransMag.at(i) );
	HF1( 1271, lp2_mom_container[best_lp2].Mag() );
	HF1( 1281, missMom );	
	HF2( 3031, lp2_mom_container[best_lp2].Mag(), momTransMag.at(i) );
	HF2( 3041, bek, momTransMag.at(i) );
	HF2( 3051, bek, event.lp2mass );      
      }	else if( event.lp2angle<40.0 ){
	HF1( 1242, event.lp2mass );      		
	HF1( 1252, bek );
	HF1( 1262, momTransMag.at(i) );	
	HF1( 1272, lp2_mom_container[best_lp2].Mag() );
	HF1( 1282, missMom );
	HF2( 3032, lp2_mom_container[best_lp2].Mag(), momTransMag.at(i) );
	HF2( 3042, bek, momTransMag.at(i) );
	HF2( 3052, bek, event.lp2mass );      	
      }	else if( event.lp2angle<60.0 ){
	HF1( 1243, event.lp2mass );      			
	HF1( 1253, bek );
	HF1( 1263, momTransMag.at(i) );
	HF1( 1273, lp2_mom_container[best_lp2].Mag() );
	HF1( 1283, missMom );
	HF2( 3033, lp2_mom_container[best_lp2].Mag(), momTransMag.at(i) );
	HF2( 3043, bek, momTransMag.at(i) );
	HF2( 3053, bek, event.lp2mass );      	
      }	else if( event.lp2angle<80.0 ){
	HF1( 1244, event.lp2mass );      			
	HF1( 1254, bek );
	HF1( 1264, momTransMag.at(i) );
	HF1( 1274, lp2_mom_container[best_lp2].Mag() );
	HF1( 1284, missMom );			
	HF2( 3034, lp2_mom_container[best_lp2].Mag(), momTransMag.at(i) );
	HF2( 3044, bek, momTransMag.at(i) );
	HF2( 3054, bek, event.lp2mass );      	
      }	else if( event.lp2angle<100.0 ){
	HF1( 1245, event.lp2mass );      			
	HF1( 1255, bek );
	HF1( 1265, momTransMag.at(i) );
	HF1( 1275, lp2_mom_container[best_lp2].Mag() );
	HF1( 1285, missMom );			
	HF2( 3035, lp2_mom_container[best_lp2].Mag(), momTransMag.at(i) );
	HF2( 3045, bek, momTransMag.at(i) );
	HF2( 3055, bek, event.lp2mass );      	
      }	else if( event.lp2angle<120.0 ){
	HF1( 1246, event.lp2mass );      			
	HF1( 1256, bek );
	HF1( 1266, momTransMag.at(i) );
	HF1( 1276, lp2_mom_container[best_lp2].Mag() );
	HF1( 1286, missMom );				
	HF2( 3036, lp2_mom_container[best_lp2].Mag(), momTransMag.at(i) );
	HF2( 3046, bek, momTransMag.at(i) );
	HF2( 3056, bek, event.lp2mass );      	
      }	else if( event.lp2angle<140.0 ){
	HF1( 1247, event.lp2mass );      			
	HF1( 1257, bek );
	HF1( 1267, momTransMag.at(i) );
	HF1( 1277, lp2_mom_container[best_lp2].Mag() );
	HF1( 1287, missMom );					
	HF2( 3037, lp2_mom_container[best_lp2].Mag(), momTransMag.at(i) );
	HF2( 3047, bek, momTransMag.at(i) );
	HF2( 3057, bek, event.lp2mass );      	
      }	else if( event.lp2angle<160.0 ){
	HF1( 1248, event.lp2mass );      			
	HF1( 1258, bek );
	HF1( 1268, momTransMag.at(i) );
	HF1( 1278, lp2_mom_container[best_lp2].Mag() );
	HF1( 1288, missMom );						
	HF2( 3038, lp2_mom_container[best_lp2].Mag(), momTransMag.at(i) );
	HF2( 3048, bek, momTransMag.at(i) );
	HF2( 3058, bek, event.lp2mass );      
      }	else if( event.lp2angle<180.0 ){
	HF1( 1249, event.lp2mass );      			
	HF1( 1259, bek );
	HF1( 1269, momTransMag.at(i) );
	HF1( 1279, lp2_mom_container[best_lp2].Mag() );
	HF1( 1289, missMom );
	HF2( 3039, lp2_mom_container[best_lp2].Mag(), momTransMag.at(i) );
	HF2( 3049, bek, momTransMag.at(i) );
	HF2( 3059, bek, event.lp2mass );      	
      }
      HF2( 2101, event.lp2angle, lp2_p_mom_container[best_lp2].Mag() );
      HF2( 2102, event.lp2angle, lp2_pi_mom_container[best_lp2].Mag() );
      HF2( 2103, event.lp2angle, lp2_p2_mom_container[best_lp2].Mag() );
      HF2( 2104, event.lp2angle, l_mom_container[best_lp2].Mag() );
      HF2( 2105, event.lp2angle, lp2_mom_container[best_lp2].Mag() );
      HF2( 2111, event.lp2angle, bek );
      HF2( 2112, event.lp2angle, event.lmass );
      HF2( 2113, event.lp2angle, event.lp2mass );

      Int_t itP1 = lp2_p_container[best_lp2];
      Int_t itPi = lp2_pi_container[best_lp2];
      Int_t itP2 = lp2_p2_container[best_lp2];
      HF2(230, event.mom0[itP1]*event.charge[itP1], event.dEdx[itP1]);      
      HF2(230, event.mom0[itPi]*event.charge[itPi], event.dEdx[itPi]);
      HF2(230, event.mom0[itP2]*event.charge[itP2], event.dEdx[itP2]);
      HF2(240, event.mom0[itP1]*event.charge[itP1], event.dEdx[itP1]);
      HF2(250, event.mom0[itPi]*event.charge[itPi], event.dEdx[itPi]);
      HF2(260, event.mom0[itP2]*event.charge[itP2], event.dEdx[itP2]);
      if ( event.charge[itP1]<0 ){
	HF2(231, -event.mom0[itP1]*event.charge[itP1], event.dEdx[itP1]);
	HF2(241, -event.mom0[itP1]*event.charge[itP1], event.dEdx[itP1]);	
      } else {
	HF2(232, event.mom0[itP1]*event.charge[itP1], event.dEdx[itP1]);
	HF2(242, event.mom0[itP1]*event.charge[itP1], event.dEdx[itP1]);	
      }
      if ( event.charge[itPi]<0 ){
	HF2(231, -event.mom0[itPi]*event.charge[itPi], event.dEdx[itPi]);
	HF2(251, -event.mom0[itPi]*event.charge[itPi], event.dEdx[itPi]);	
      } else {
	HF2(232, event.mom0[itPi]*event.charge[itPi], event.dEdx[itPi]);
	HF2(252, event.mom0[itPi]*event.charge[itPi], event.dEdx[itPi]);	
      }
      if ( event.charge[itP2]<0 ){
	HF2(231, -event.mom0[itP2]*event.charge[itP2], event.dEdx[itP2]);
	HF2(261, -event.mom0[itP2]*event.charge[itP2], event.dEdx[itP2]);	
      } else {
	HF2(232, event.mom0[itP2]*event.charge[itP2], event.dEdx[itP2]);
	HF2(262, event.mom0[itP2]*event.charge[itP2], event.dEdx[itP2]);	
      }      
    } // nKK
  } // bestLp2
  // }
  
  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;  
  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCTracks(**src.ntTpc, **src.isK18, **src.isKurama,
			 **src.charge, **src.nhtrack, **src.helix_cx,
			 **src.helix_cy, **src.helix_z0, **src.helix_r,
			 **src.helix_dz, **src.hitlayer, **src.track_cluster_mrow,
			 **src.helix_t, **src.track_cluster_de, **src.resolution_x,
			 **src.resolution_y, **src.resolution_z, **src.hitpos_x,
			 **src.hitpos_y, **src.hitpos_z);

  //  HF1( 1, event.status++ ); // status 12
  
  HF1( 1, 16 );// for debug // status  12
  
  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
  
  if(best!=-1){
    Int_t id_p = l_p_container[best];
    Int_t id_pi = l_pi_container[best];
    auto Track_p = TPCAna.GetTrackTPCHelix(id_p);
    auto Track_pi = TPCAna.GetTrackTPCHelix(id_pi);
    for(Int_t ih=0;ih<Track_p->GetNHit();++ih){
      auto pos = Track_p->GetHitInOrder(ih)->GetLocalHitPos();
      HF2(1001,pos.z(),pos.x());
    }
    for(Int_t ih=0;ih<Track_pi->GetNHit();++ih){
      auto pos = Track_pi->GetHitInOrder(ih)->GetLocalHitPos();
      HF2(1002,pos.z(),pos.x());
    }
    if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;  
    event.ldecays_res_mom.push_back(Track_p->GetMomentumResolution());
    event.ldecays_res_mom_x.push_back(Track_p->GetMomentumResolutionVect().X());
    event.ldecays_res_mom_y.push_back(Track_p->GetMomentumResolutionVect().Y());
    event.ldecays_res_mom_z.push_back(Track_p->GetMomentumResolutionVect().Z());
    event.ldecays_res_mom_t.push_back(Track_p->GetTransverseMomentumResolution());
    event.ldecays_res_th.push_back(Track_p->GetThetaResolution());
    event.ldecays_res_ph.push_back(Track_p->GetTransverseAngularResolution());

    event.ldecays_cov_mom_ph.push_back(Track_p->GetTransverseMomentumAngularCovariance());
    event.ldecays_cov_mom_th.push_back(Track_p->GetMomentumPitchAngleCovariance());
    event.ldecays_cov_mom_xy.push_back(Track_p->GetMomentumCovarianceVect().x());
    event.ldecays_cov_mom_yz.push_back(Track_p->GetMomentumCovarianceVect().y());
    event.ldecays_cov_mom_zx.push_back(Track_p->GetMomentumCovarianceVect().z());

    event.ldecays_res_mom.push_back(Track_pi->GetMomentumResolution());
    event.ldecays_res_mom_x.push_back(Track_pi->GetMomentumResolutionVect().X());
    event.ldecays_res_mom_y.push_back(Track_pi->GetMomentumResolutionVect().Y());
    event.ldecays_res_mom_z.push_back(Track_pi->GetMomentumResolutionVect().Z());
    event.ldecays_res_mom_t.push_back(Track_pi->GetTransverseMomentumResolution());
    event.ldecays_res_th.push_back(Track_pi->GetThetaResolution());
    event.ldecays_res_ph.push_back(Track_pi->GetTransverseAngularResolution());

    event.ldecays_cov_mom_ph.push_back(Track_pi->GetTransverseMomentumAngularCovariance());
    event.ldecays_cov_mom_th.push_back(Track_pi->GetMomentumPitchAngleCovariance());
    event.ldecays_cov_mom_xy.push_back(Track_pi->GetMomentumCovarianceVect().x());
    event.ldecays_cov_mom_yz.push_back(Track_pi->GetMomentumCovarianceVect().y());
    event.ldecays_cov_mom_zx.push_back(Track_pi->GetMomentumCovarianceVect().z());
  } // bestL
  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
  
  //  if(best_lp2 != -1){
  if(best!=-1 && lp2_candidates==1){
    Int_t best_lp2=0;
    //for(int best_lp2=0; best_lp2<lp2_candidates; best_lp2++){
    Int_t id_p_lp2 = lp2_p_container[best_lp2];
    Int_t id_pi_lp2 = lp2_pi_container[best_lp2];
    Int_t id_p2_lp2 = lp2_p2_container[best_lp2];

    auto Track_p_lp2 = TPCAna.GetTrackTPCHelix(id_p_lp2);
    auto Track_pi_lp2 = TPCAna.GetTrackTPCHelix(id_pi_lp2);
    auto Track_p2_lp2 = TPCAna.GetTrackTPCHelix(id_p2_lp2);
    for(Int_t ih=0;ih<Track_p_lp2->GetNHit();++ih){
      auto pos = Track_p_lp2->GetHitInOrder(ih)->GetLocalHitPos();
      //    HF2(1001,pos.z(),pos.x());
    }
    for(Int_t ih=0;ih<Track_pi_lp2->GetNHit();++ih){
      auto pos = Track_pi_lp2->GetHitInOrder(ih)->GetLocalHitPos();
      //    HF2(1002,pos.z(),pos.x());
    }
    for(Int_t ih=0;ih<Track_p2_lp2->GetNHit();++ih){
      auto pos = Track_p2_lp2->GetHitInOrder(ih)->GetLocalHitPos();
      //    HF2(1001,pos.z(),pos.x());
    }
    event.lp2decays_res_mom.push_back(Track_p_lp2->GetMomentumResolution());
    event.lp2decays_res_mom_x.push_back(Track_p_lp2->GetMomentumResolutionVect().X());
    event.lp2decays_res_mom_y.push_back(Track_p_lp2->GetMomentumResolutionVect().Y());
    event.lp2decays_res_mom_z.push_back(Track_p_lp2->GetMomentumResolutionVect().Z());
    event.lp2decays_res_mom_t.push_back(Track_p_lp2->GetTransverseMomentumResolution());
    event.lp2decays_res_th.push_back(Track_p_lp2->GetThetaResolution());
    event.lp2decays_res_ph.push_back(Track_p_lp2->GetTransverseAngularResolution());

    event.lp2decays_cov_mom_ph.push_back(Track_p_lp2->GetTransverseMomentumAngularCovariance());
    event.lp2decays_cov_mom_th.push_back(Track_p_lp2->GetMomentumPitchAngleCovariance());
    event.lp2decays_cov_mom_xy.push_back(Track_p_lp2->GetMomentumCovarianceVect().x());
    event.lp2decays_cov_mom_yz.push_back(Track_p_lp2->GetMomentumCovarianceVect().y());
    event.lp2decays_cov_mom_zx.push_back(Track_p_lp2->GetMomentumCovarianceVect().z());

    event.lp2decays_res_mom.push_back(Track_pi_lp2->GetMomentumResolution());
    event.lp2decays_res_mom_x.push_back(Track_pi_lp2->GetMomentumResolutionVect().X());
    event.lp2decays_res_mom_y.push_back(Track_pi_lp2->GetMomentumResolutionVect().Y());
    event.lp2decays_res_mom_z.push_back(Track_pi_lp2->GetMomentumResolutionVect().Z());
    event.lp2decays_res_mom_t.push_back(Track_pi_lp2->GetTransverseMomentumResolution());
    event.lp2decays_res_th.push_back(Track_pi_lp2->GetThetaResolution());
    event.lp2decays_res_ph.push_back(Track_pi_lp2->GetTransverseAngularResolution());

    event.lp2decays_cov_mom_ph.push_back(Track_pi_lp2->GetTransverseMomentumAngularCovariance());
    event.lp2decays_cov_mom_th.push_back(Track_pi_lp2->GetMomentumPitchAngleCovariance());
    event.lp2decays_cov_mom_xy.push_back(Track_pi_lp2->GetMomentumCovarianceVect().x());
    event.lp2decays_cov_mom_yz.push_back(Track_pi_lp2->GetMomentumCovarianceVect().y());
    event.lp2decays_cov_mom_zx.push_back(Track_pi_lp2->GetMomentumCovarianceVect().z());

    event.lp2decays_res_mom.push_back(Track_p2_lp2->GetMomentumResolution());
    event.lp2decays_res_mom_x.push_back(Track_p2_lp2->GetMomentumResolutionVect().X());
    event.lp2decays_res_mom_y.push_back(Track_p2_lp2->GetMomentumResolutionVect().Y());
    event.lp2decays_res_mom_z.push_back(Track_p2_lp2->GetMomentumResolutionVect().Z());
    event.lp2decays_res_mom_t.push_back(Track_p2_lp2->GetTransverseMomentumResolution());
    event.lp2decays_res_th.push_back(Track_p2_lp2->GetThetaResolution());
    event.lp2decays_res_ph.push_back(Track_p2_lp2->GetTransverseAngularResolution());

    event.lp2decays_cov_mom_ph.push_back(Track_p2_lp2->GetTransverseMomentumAngularCovariance());
    event.lp2decays_cov_mom_th.push_back(Track_p2_lp2->GetMomentumPitchAngleCovariance());
    event.lp2decays_cov_mom_xy.push_back(Track_p2_lp2->GetMomentumCovarianceVect().x());
    event.lp2decays_cov_mom_yz.push_back(Track_p2_lp2->GetMomentumCovarianceVect().y());
    event.lp2decays_cov_mom_zx.push_back(Track_p2_lp2->GetMomentumCovarianceVect().z());
  } // bestLp2
  
  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
  
  HypTPCTask& GFTrackCont = HypTPCTask::GetInstance();
  
  if(best!=-1){
    HF1( 10, event.MissMass[0] );
    HF1( 11, event.lmass );
    HF1( 12, event.lp2mass );
    HF1( 13, event.ldecays_mom[0] );
    HF1( 14, event.ldecays_mom[1] );
    HF1( 15, event.ldecays_mom[2] );
    HF1( 18, event.ppi_dist );  
    HF1( 103, event.ldecays_mom[0] );
    HF1( 104, event.ldecays_mom[1] );
    HF1( 105, event.ldecays_mom[2] );  
    HF1( 16, TMath::Sqrt( event.lmom_x*event.lmom_x +
			  event.lmom_y*event.lmom_y +
			  event.lmom_z*event.lmom_z) );
    // HF1( 17, TMath::Sqrt( event.ximom_x*event.ximom_x +
    // 			event.ximom_y*event.ximom_y +
    // 			event.ximom_z*event.ximom_z) );
    
  } // bestL
  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix(it);
    if( !tp ) continue;
    std::vector<Int_t> pdgcode;
    Kinematics::HypTPCPID_PDGCode(event.charge[it], event.pid[it], pdgcode);
    GFTrackCont.AddHelixTrack(pdgcode, tp);
  }
  HF1( 2, event.GFstatus++ );
  HF1( 1, event.status++ );  
  GFTrackCont.FitTracks();
  
  Int_t GFntTpc = GFTrackCont.GetNTrack();
  if(GFntTpc!=event.ntTpc && best!=-1){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    return true;
  }

  event.GFinside.resize(GFntTpc);  

  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
    
  if(best!=-1){
    HF1( 2, event.GFstatus++ );
    HF1( 1, event.status++ );  
    // lambda
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
    // lambda proton
    std::vector<Int_t> GFlp2_p_id_container(lp2_candidates, -1);
    std::vector<Int_t> GFlp2_pi_id_container(lp2_candidates, -1);
    std::vector<Int_t> GFlp2_p2_id_container(lp2_candidates, -1);
    std::vector<Int_t> GFlp2_p_rep_container(lp2_candidates, -1);
    std::vector<Int_t> GFlp2_pi_rep_container(lp2_candidates, -1);
    std::vector<Int_t> GFlp2_p2_rep_container(lp2_candidates, -1);
    std::vector<TVector3> GFlp2_mom_container(lp2_candidates, TVector3(qnan, qnan, qnan));
    std::vector<TVector3> GFlp2_p_mom_container(lp2_candidates, TVector3(qnan, qnan, qnan));
    std::vector<TVector3> GFlp2_pi_mom_container(lp2_candidates, TVector3(qnan, qnan, qnan));
    std::vector<TVector3> GFlp2_p2_mom_container(lp2_candidates, TVector3(qnan, qnan, qnan));
    std::vector<TVector3> GFlp2_l_mom_container(lp2_candidates, TVector3(qnan, qnan, qnan));  
    std::vector<TVector3> GFlp2_vert_container(lp2_candidates, TVector3(qnan, qnan, qnan));
    std::vector<Double_t> GFl_lp2_mass_container(lp2_candidates, qnan);
    std::vector<Double_t> GFlp2_mass_container(lp2_candidates, qnan);
    std::vector<Double_t> GFlp2_p_mass_container(lp2_candidates, qnan);
    std::vector<Double_t> GFlp2_pi_mass_container(lp2_candidates, qnan);
    std::vector<Double_t> GFlp2_p2_mass_container(lp2_candidates, qnan);
    std::vector<Double_t> GFlp2_closedist_container(lp2_candidates, qnan);
    std::vector<Double_t> GFlp2_opangle_container(lp2_candidates, qnan);
    // kinematical fitting
    std::vector<Double_t> KFchisqrl_container(l_candidates, qnan);
    std::vector<Double_t> KFpvall_container(l_candidates, qnan);
    std::vector<std::vector<Double_t>> KFlpull_container(l_candidates, std::vector<Double_t>(6, qnan));
    std::vector<TMatrixD> KFVarianceLd_container(l_candidates, TMatrixD(3,3));
    std::vector<TVector3> KFp_mom_container(l_candidates, TVector3(qnan, qnan, qnan));
    std::vector<TVector3> KFpi_mom_container(l_candidates, TVector3(qnan, qnan, qnan));
    Int_t gfbest = -1; prev_massdiff = 9999.;

    if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
  
    for(Int_t candi=0;candi<l_candidates;candi++){
      Int_t trackid_p = l_p_container[candi];
      Int_t trackid_pi = l_pi_container[candi];
      Int_t flag = 1; Int_t repid_p = 0;
      for(Int_t i=0;i<2;i++){
	Int_t temp = flag&event.pid[trackid_p];
	if(temp==flag) repid_p += 1;
	flag*=2;
      }
      //    Int_t repid_pi = 0; Int_t repid_pi2 = 0;
      Int_t repid_pi = 0; Int_t repid_p2 = 0;    

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
      TLorentzVector GFLlambda_fixed(GFlambda_mom, TMath::Sqrt(GFlambda_mom.Mag()*GFlambda_mom.Mag() + LambdaMass*LambdaMass));
      //    TVector3 GFxi_vert; Double_t GFlpi_dist = qnan; Double_t GFlambda_tracklen;
      TVector3 GFlp_vert; Double_t GFlp_dist = qnan; Double_t GFlambda_tracklen;
      Double_t GFlambda_tof = Kinematics::CalcTimeOfFlight(GFlambda_mom.Mag(), GFlambda_tracklen, pdg::LambdaMass());
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
      Double_t diff = TMath::Abs(lambda_mass_container[candi] - LambdaMass);
      if(best!=-1 && prev_massdiff > diff){
	prev_massdiff = diff;
	gfbest = candi;
      }
    }
    //  } // gfbest!=-1
    if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
  
    // lambda proton 
    Int_t gfbest_lp2 = -1;
    //  prev_massdiff = 9999.;
    Double_t gfprev_distdiff_lp2 = 9999.;  
    //  for(Int_t candi=0;candi<lp2_candidates;candi++){
    //  if(gfbest!=-1 && lp2_candidates==1){
    Bool_t findVtxGFppi = false;
    Bool_t findVtxGFp2L = false;
    if(lp2_candidates==1){  
      Int_t candi=0;
      if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;    
      Int_t trackid_p = lp2_p_container[candi];
      Int_t trackid_pi = lp2_pi_container[candi];
      Int_t trackid_p2 = lp2_p2_container[candi];
      if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;          
      Int_t flag = 1; Int_t repid_p = 0;
      for(Int_t i=0;i<2;i++){
	Int_t temp = flag&event.pid[trackid_p];
	if(temp==flag) repid_p += 1;
	flag*=2;
      }
      if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;                
      Int_t repid_pi = 0; Int_t repid_p2 = 0;
      Double_t GFextrapolation_decays[3];
      Double_t GFmass2_decays[3] = {qnan, qnan, qnan};
      TVector3 GFmom_decays[3]; TVector3 GFlambda_vert; double GFppi_dist=qnan;
      if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
      findVtxGFppi = (GFTrackCont.FindVertex(trackid_p, trackid_pi,
					     repid_p, repid_pi,
					     GFextrapolation_decays[0], GFextrapolation_decays[1],
					     GFmom_decays[0], GFmom_decays[1],
					     GFppi_dist, GFlambda_vert,
					     vtx_scan_range) ) && (GFppi_dist < GFppi_distcut);
      // if(GFTrackCont.FindVertex(trackid_p, trackid_pi,
      // 				repid_p, repid_pi,
      // 				GFextrapolation_decays[0], GFextrapolation_decays[1],
      // 				GFmom_decays[0], GFmom_decays[1],
      // 				GFppi_dist, GFlambda_vert,
      // 				vtx_scan_range)
      // 	 || GFppi_dist > GFppi_distcut){
      if(findVtxGFppi){
	if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
	TLorentzVector GFLp(GFmom_decays[0], TMath::Sqrt(GFmom_decays[0].Mag()*GFmom_decays[0].Mag() + ProtonMass*ProtonMass));
	TLorentzVector GFLpi(GFmom_decays[1], TMath::Sqrt(GFmom_decays[1].Mag()*GFmom_decays[1].Mag() + PionMass*PionMass));
	//    std::cout << "debug GFmom_decays[2].Mag() " << GFmom_decays[2].Mag() << std::endl;
	TLorentzVector GFLlambda = GFLp + GFLpi;
	TVector3 GFlambda_mom = GFmom_decays[0] + GFmom_decays[1];
	TLorentzVector GFLlambda_fixed(GFlambda_mom, TMath::Sqrt(GFlambda_mom.Mag()*GFlambda_mom.Mag() + LambdaMass*LambdaMass));
	TVector3 GFlp2_vert; Double_t GFlp2_dist = qnan; Double_t GFlambda_tracklen;
	findVtxGFp2L = (GFTrackCont.FindVertexXi(trackid_p2, repid_p2,
						 GFlambda_vert, GFlambda_mom, GFlambda_tracklen,
						 GFextrapolation_decays[2], GFmom_decays[2],
						 GFlp2_dist, GFlp2_vert, vtx_scan_range) ) && (GFlp2_dist < GFlp2_distcut);
	// if(GFTrackCont.FindVertexXi(trackid_p2, repid_p2,
	// 			     GFlambda_vert, GFlambda_mom, GFlambda_tracklen,
	// 			     GFextrapolation_decays[2], GFmom_decays[2],
	// 			     GFlp2_dist, GFlp2_vert, vtx_scan_range)
	//    || GFlp2_dist > GFlp2_distcut){
	if(findVtxGFp2L){
	  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
	  Double_t GFlambda_tof = Kinematics::CalcTimeOfFlight(GFlambda_mom.Mag(), GFlambda_tracklen, pdg::LambdaMass());    
	  Int_t htofhitid_p; Double_t tracklen_p; Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	  Bool_t htofextrapolation_p = GFTrackCont.TPCHTOFTrackMatching(trackid_p, repid_p, GFlambda_vert, event.HtofSeg, event.posHtof, htofhitid_p, tof, tracklen_p, pos, track2tgt_dist);
	  if(htofextrapolation_p){
	    if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
	    GFmass2_decays[0] = Kinematics::MassSquare(GFmom_decays[0].Mag(),
						       tracklen_p  - GFextrapolation_decays[0],
						       event.tHtof[htofhitid_p] - GFlambda_tof);
	    //if(GFmass2_decays[0] < 0.25) continue;
	  }
	  Int_t htofhitid_pi; Double_t tracklen_pi;
	  Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(trackid_pi, repid_pi, GFlambda_vert, event.HtofSeg, event.posHtof,htofhitid_pi, tof, tracklen_pi, pos, track2tgt_dist);
	  if(htofextrapolation_pi){
	    if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;                	    
	    GFmass2_decays[1] = Kinematics::MassSquare(GFmom_decays[1].Mag(),
						       tracklen_pi - GFextrapolation_decays[1],
						       event.tHtof[htofhitid_pi] - GFlambda_tof);
	    // if(GFmass2_decays[1] > 0.25) continue;
	  }    
	  TLorentzVector GFLp2(GFmom_decays[2], TMath::Sqrt(GFmom_decays[2].Mag()*GFmom_decays[2].Mag() + ProtonMass*ProtonMass));
	  TLorentzVector GFLlp2 = GFLlambda_fixed + GFLp2;
	  TVector3 GFlp2_mom = GFmom_decays[0] + GFmom_decays[1] + GFmom_decays[2];
    
	  Int_t htofhitid_p2; Double_t tracklen_p2; 
	  Bool_t htofextrapolation_p2 = GFTrackCont.TPCHTOFTrackMatching(trackid_p2, repid_p2, GFlambda_vert, event.HtofSeg, event.posHtof, htofhitid_p2, tof, tracklen_p2, pos, track2tgt_dist);
	  if(htofextrapolation_p2){
	    if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;                	    
	    GFmass2_decays[2] = Kinematics::MassSquare(GFmom_decays[2].Mag(),
						       tracklen_p2  - GFextrapolation_decays[2],
						       event.tHtof[htofhitid_p2]);
	    //if(GFmass2_decays[0] < 0.25) continue;
	  }
	  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;                	  
	  Double_t lp2_opening_angle = GFlambda_mom.Angle(GFmom_decays[2])*TMath::RadToDeg();    
	  event.GFlp2flag = true;
	  GFlp2_p_id_container[candi] = trackid_p;	  
	  GFlp2_pi_id_container[candi] = trackid_pi;
	  GFlp2_p2_id_container[candi] = trackid_p2;
	  if(FlagDebug) std::cout << "debug " << __LINE__ <<   " trackid_p:" << trackid_p << ", trackid_pi:" << trackid_pi << ", trackid_p2:" << trackid_p2 << std::endl;                	  	 
	  GFlp2_p_rep_container[candi] = repid_p;
	  GFlp2_pi_rep_container[candi] = repid_pi;
	  GFlp2_p2_rep_container[candi] = repid_p2;
	  GFlp2_mom_container[candi] = GFlp2_mom; // 
	  GFlp2_p_mom_container[candi] = GFmom_decays[0]; 
	  GFlp2_pi_mom_container[candi] = GFmom_decays[1]; 
	  GFlp2_p2_mom_container[candi] = GFmom_decays[2];
	  GFlp2_l_mom_container[candi] = GFlambda_mom;
	  GFlp2_vert_container[candi] = GFlp2_vert;
	  GFl_lp2_mass_container[candi] = GFLlambda.M();
	  GFlp2_mass_container[candi] = GFLlp2.M();
	  GFlp2_p_mass_container[candi] = GFmass2_decays[0];
	  GFlp2_pi_mass_container[candi] = GFmass2_decays[1];
	  GFlp2_p2_mass_container[candi] = GFmass2_decays[2];
	  GFlp2_closedist_container[candi] = GFlp2_dist;
	  GFlp2_opangle_container[candi] = lp2_opening_angle;
	  //    Double_t diff = TMath::Abs(lambda_mass_container[candi] - LambdaMass); // to be changed
	  Double_t diff = TMath::Abs(GFlp2_closedist_container.at(candi));
	  // if(gfprev_distdiff_lp2 > diff){
	  //   gfprev_distdiff_lp2 = diff;
	  //   gfbest_lp2 = candi;
	  //   std::cout << "GFbest lambda proton: " << best_lp2 << std::endl;      
	  // }
	}//
      }
    }
    if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;                
  
    if(event.GFlflag){
      if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;                      
      HF1( genfitHid, GFntTpc);
      HF1( 20, event.MissMass[0] );
      HF1( 21, event.lmass );
      HF1( 22, event.lp2mass );
      HF1( 23, event.ldecays_mom[0] );
      HF1( 24, event.ldecays_mom[1] );
      HF1( 26, TMath::Sqrt( event.lmom_x*event.lmom_x +
			    event.lmom_y*event.lmom_y +
			    event.lmom_z*event.lmom_z) );
      HF1( 27, TMath::Sqrt( event.lp2mom_x*event.lp2mom_x +
			    event.lp2mom_y*event.lp2mom_y +
			    event.lp2mom_z*event.lp2mom_z) );
      HF1( 28, event.ppi_dist);
      event.GFntTpc = 3; // is this correct? 
      if(event.GFlp2flag&&lp2_candidates==1){
	if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;                	
	Int_t gfbest_lp2=0;
	event.GFlp2mass = GFlp2_mass_container[gfbest_lp2];
	event.GFlmass_lp2 = GFl_lp2_mass_container[gfbest_lp2];
	event.GFlp2decayvtx_x = GFlp2_vert_container[gfbest_lp2].x();
	event.GFlp2decayvtx_y = GFlp2_vert_container[gfbest_lp2].y();
	event.GFlp2decayvtx_z = GFlp2_vert_container[gfbest_lp2].z();
	event.GFlp2mom = GFlp2_mom_container[gfbest_lp2].Mag();
	event.GFlp2mom_x = GFlp2_mom_container[gfbest_lp2].x();
	event.GFlp2mom_y = GFlp2_mom_container[gfbest_lp2].y();
	event.GFlp2mom_z = GFlp2_mom_container[gfbest_lp2].z();
	event.GFlmom_lp2 = GFlp2_l_mom_container[gfbest_lp2].Mag();
	event.GFlp2_dist = GFlp2_closedist_container[gfbest_lp2];
	event.GFlp2_angle = GFlp2_opangle_container[gfbest_lp2];
    
	event.GFdecays_m2_lp2.resize(3);
	event.GFdecays_mom_lp2.resize(3);
	event.GFdecays_mom_x_lp2.resize(3);
	event.GFdecays_mom_y_lp2.resize(3);
	event.GFdecays_mom_z_lp2.resize(3);
	event.GFmomloss_lp2.resize(3);
	event.GFeloss_lp2.resize(3);

	event.GFcharge_lp2.resize(3);
	event.GFchisqr_lp2.resize(3);
	event.GFtof_lp2.resize(3);
	event.GFtracklen_lp2.resize(3);
	event.GFpval_lp2.resize(3);
	event.GFchisqrPos_lp2.resize(3);
	event.GFpvalPos_lp2.resize(3);
	event.GFpdgcode_lp2.resize(3);

	event.GFfitstatus_lp2.resize(3);
	event.GFnhtrack_lp2.resize(3);
	event.GFlayer_lp2.resize(3);
	event.GFrow_lp2.resize(3);
	event.GFpos_x_lp2.resize(3);
	event.GFpos_y_lp2.resize(3);
	event.GFpos_z_lp2.resize(3);
	event.GFmom_lp2.resize(3);
	event.GFmom_x_lp2.resize(3);
	event.GFmom_y_lp2.resize(3);
	event.GFmom_z_lp2.resize(3);
	event.GFresidual_x_lp2.resize(3);
	event.GFresidual_y_lp2.resize(3);
	event.GFresidual_z_lp2.resize(3);
	event.GFresidual_px_lp2.resize(3);
	event.GFresidual_py_lp2.resize(3);
	event.GFresidual_pz_lp2.resize(3);
	event.GFresidual_p_lp2.resize(3);
	event.GFresidual6D_x_lp2.resize(3);
	event.GFresidual6D_y_lp2.resize(3);
	event.GFresidual6D_z_lp2.resize(3);
	event.GFresidual6D_px_lp2.resize(3);
	event.GFresidual6D_py_lp2.resize(3);
	event.GFresidual6D_pz_lp2.resize(3);
	event.GFresolution_x_lp2.resize(3);
	event.GFresolution_y_lp2.resize(3);
	event.GFresolution_z_lp2.resize(3);
	event.GFresolution_p_lp2.resize(3);
	event.GFresolution_px_lp2.resize(3);
	event.GFresolution_py_lp2.resize(3);
	event.GFresolution_pz_lp2.resize(3);
	event.GFpull_x_lp2.resize(3);
	event.GFpull_y_lp2.resize(3);
	event.GFpull_z_lp2.resize(3);
	event.GFpull_p_lp2.resize(3);
	event.GFpull_px_lp2.resize(3);
	event.GFpull_py_lp2.resize(3);
	event.GFpull_pz_lp2.resize(3);
      } // GFlp2flag && Lp2candidate==1
      if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
      
      event.GFlmass = GFl_mass_container[gfbest];
      event.GFldecayvtx_x = GFl_vert_container[gfbest].x();
      event.GFldecayvtx_y = GFl_vert_container[gfbest].y();
      event.GFldecayvtx_z = GFl_vert_container[gfbest].z();
      event.GFlmom = GFl_mom_container[gfbest].Mag();
      event.GFlmom_x = GFl_mom_container[gfbest].x();
      event.GFlmom_y = GFl_mom_container[gfbest].y();
      event.GFlmom_z = GFl_mom_container[gfbest].z();
      event.GFppi_dist = GFppi_closedist_container[gfbest];
  
      event.GFdecays_m2.resize(3);
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
      event.GFchisqrPos.resize(3);
      event.GFpvalPos.resize(3);
      event.GFpdgcode.resize(3);

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
      event.GFresidual_px.resize(3);
      event.GFresidual_py.resize(3);
      event.GFresidual_pz.resize(3);
      event.GFresidual_p.resize(3);
      event.GFresidual6D_x.resize(3);
      event.GFresidual6D_y.resize(3);
      event.GFresidual6D_z.resize(3);
      event.GFresidual6D_px.resize(3);
      event.GFresidual6D_py.resize(3);
      event.GFresidual6D_pz.resize(3);
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
      
      if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
      
      // lambda
      for( Int_t j=0; j<2; ++j ){
	Int_t igf = GFl_p_id_container[gfbest];    
	if(j==1) igf = GFl_pi_id_container[gfbest];

	Int_t repid = GFl_p_rep_container[gfbest];
	if(j==1) repid = GFl_pi_rep_container[gfbest];

	event.GFfitstatus[j] = (Int_t) GFTrackCont.TrackCheck(igf, repid);
	HF1( 3, event.GFfitstatus[j]);
	if(!event.GFfitstatus[j]) continue;

	Int_t nh = GFTrackCont.GetNHits(igf);
	event.GFnhtrack[j] = nh;
	event.GFchisqr[j] = GFTrackCont.GetChi2NDF(igf, repid);
	event.GFcharge[j] = GFTrackCont.GetCharge(igf, repid);
	event.GFtof[j] = GFTrackCont.GetTrackTOF(igf, 0, -1, repid);
	event.GFtracklen[j] = GFTrackCont.GetTrackLength(igf, 0, -1, repid);
	event.GFpval[j] = GFTrackCont.GetPvalue(igf, repid);
	event.GFpdgcode[j] = GFTrackCont.GetPDGcode(igf, repid);

	TVector3 GFmom_decays = GFp_mom_container[gfbest];
	if(j==1) GFmom_decays = GFpi_mom_container[gfbest];

	Double_t GFmass_decays = GFp_mass_container[gfbest];
	if(j==1) GFmass_decays = GFpi_mass_container[gfbest];

	event.GFdecays_m2[j] = GFmass_decays;
	event.GFdecays_mom[j] = GFmom_decays.Mag();

	//    if(j==1) HF2( )
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
	event.GFresidual_px[j].resize(nh);
	event.GFresidual_py[j].resize(nh);
	event.GFresidual_pz[j].resize(nh);
	event.GFresidual_p[j].resize(nh);
	event.GFresidual6D_x[j].resize(nh);
	event.GFresidual6D_y[j].resize(nh);
	event.GFresidual6D_z[j].resize(nh);
	event.GFresidual6D_px[j].resize(nh);
	event.GFresidual6D_py[j].resize(nh);
	event.GFresidual6D_pz[j].resize(nh);
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

	HF1( 600, event.GFpval[j]);

	Int_t id = GFl_p_id_container[gfbest];
	if(j==1) id = GFl_pi_id_container[gfbest];

	Int_t ihit = 0;
	TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( id );
	double GFchisqrPos=0;
	for( Int_t ih=0; ih<tp -> GetNHit(); ++ih ){
	  Int_t layer = (Int_t) event.hitlayer[id][ih];
	  TPCLTrackHit *helix_point = tp -> GetHitInOrder(ih);
	  if(!helix_point -> IsGoodForTracking()) continue;

	  const TVector3 &hit0 = helix_point -> GetLocalHitPos();
	  double row = helix_point->GetMRow()+0.5;
	  TVector3 mom0 = helix_point -> GetMomentumHelix(event.GFcharge[j]);
	  TVector3 hit = GFTrackCont.GetPos(igf, ihit, repid);
	  TVector3 mom = GFTrackCont.GetMom(igf, ihit, repid);

	  Double_t residual_[5];
	  Double_t pull_[5];
	  Double_t GFresidual6D[6];
	  Double_t GFpull6D[6];
	  TVector3 dumV;
	  Double_t dumd;
	  GFTrackCont.GetTrackPull(igf, event.GFpdgcode[j], dumV,
				   dumd, mom0, hit0, residual_,
				   pull_, GFresidual6D, GFpull6D);

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
	  event.GFresidual_px[j][ihit] = mom.x() - mom0.x();
	  event.GFresidual_py[j][ihit] = mom.y() - mom0.y();
	  event.GFresidual_pz[j][ihit] = mom.z() - mom0.z();
	  event.GFresidual_p[j][ihit] = mom.Mag() - mom0.Mag();

	  event.GFresidual6D_x[j][ihit] =GFresidual6D[0];
	  event.GFresidual6D_y[j][ihit] =GFresidual6D[1];
	  event.GFresidual6D_z[j][ihit] =GFresidual6D[2];
	  event.GFresidual6D_px[j][ihit] =GFresidual6D[3];
	  event.GFresidual6D_py[j][ihit] =GFresidual6D[4];
	  event.GFresidual6D_pz[j][ihit] =GFresidual6D[5];

	  event.GFresolution_x[j][ihit] = GFresidual6D[0]/GFpull6D[0];
	  event.GFresolution_y[j][ihit] = GFresidual6D[1]/GFpull6D[1];
	  event.GFresolution_z[j][ihit] = GFresidual6D[2]/GFpull6D[2];
	  event.GFresolution_px[j][ihit] = GFresidual6D[3]/GFpull6D[3];
	  event.GFresolution_py[j][ihit] = GFresidual6D[4]/GFpull6D[4];
	  event.GFresolution_pz[j][ihit] = GFresidual6D[5]/GFpull6D[5];
	  event.GFpull_x[j][ihit] = GFpull6D[0];
	  event.GFpull_y[j][ihit] = GFpull6D[1];
	  event.GFpull_z[j][ihit] = GFpull6D[2];
	  event.GFpull_px[j][ihit] = GFpull6D[3];
	  event.GFpull_py[j][ihit] = GFpull6D[4];
	  event.GFpull_pz[j][ihit] = GFpull6D[5];
	  double GFresolution_t = hypot(event.GFresolution_x[j][ihit],event.GFresolution_z[j][ihit]);
	  double GFresidual_t = hypot(event.GFresidual_x[j][ihit],event.GFresidual_z[j][ihit]);
	  GFchisqrPos+= hypot(GFresidual_t/GFresolution_t,event.GFresidual_y[j][ihit]/(event.GFresolution_y[j][ihit]));

	  int hn = 2001 + j;
	  HF2(hn,hit.z(),hit.x());
	  HF1( 601, event.GFpull_x[j][ihit]);
	  HF1( 602, event.GFpull_y[j][ihit]);
	  HF1( 603, event.GFpull_z[j][ihit]);
	  HF1( 604, event.GFpull_px[j][ihit]);
	  HF1( 605, event.GFpull_py[j][ihit]);
	  HF1( 606, event.GFpull_pz[j][ihit]);
	  if(event.GFpval[j]>0.01){
	    HF1( 611, event.GFpull_x[j][ihit]);
	    HF1( 612, event.GFpull_y[j][ihit]);
	    HF1( 613, event.GFpull_z[j][ihit]);
	    HF1( 614, event.GFpull_px[j][ihit]);
	    HF1( 615, event.GFpull_py[j][ihit]);
	    HF1( 616, event.GFpull_pz[j][ihit]);
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
  
	HF1( 10103+j, GFmom_decays.Mag() );
	HF1( 10121+j, event.GFdecays_m2[j] );

	if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;                  
	for(int i=0; i<event.nKK; ++i){
	  if(momK18[i].Mag()<0.025) continue;      
	  if(momKurama[i].Mag()<0.025) continue;
	  Double_t bek = event.BEkaonTPC[i];
	  Double_t theta = event.thetaTPC[i];
	  HF1( 11401+j*10, event.GFdecays_mom[j] );
	  for( Int_t i_bek=0; i_bek<5; i_bek++){
	    if( bek<-0.1+0.1*i_bek ){
	      HF1( 11402+j*10+i_bek, event.GFdecays_mom[j] );
	      break;
	    }
	  }
	} // nKp
      } //igf

      for(int i=0; i<event.nKK; ++i){
	if(momK18[i].Mag()<0.025) continue;      
	if(momKurama[i].Mag()<0.025) continue;
	Double_t bek = event.BEkaonTPC[i];
	Double_t theta = event.thetaTPC[i];
	if(theta<5.0){
	  HF1( 11205, bek);
	}
	if(theta<10.0){
	  HF1( 11206, bek);
	}
	for(Int_t i_bek=0; i_bek<5; i_bek++){
	  if( bek<-0.1+0.1*i_bek ){
	    HF1( 11212+i_bek, event.GFlmass );
	    break;
	  }
	}    
      }

      if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
      // lambda proton
      if(event.GFlp2flag&&lp2_candidates==1){
	Int_t gfbest_lp2=0;
	for( Int_t j=0; j<3; ++j ){
	  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;                	  
	  Int_t igf = GFlp2_p_id_container[gfbest_lp2];
	  if(j==1) igf = GFlp2_pi_id_container[gfbest_lp2];
	  if(j==2) igf = GFlp2_p2_id_container[gfbest_lp2];
	  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;                	  
	  Int_t repid = GFlp2_p_rep_container[gfbest_lp2];
	  if(j==1) repid = GFlp2_pi_rep_container[gfbest_lp2];
	  if(j==2) repid = GFlp2_p2_rep_container[gfbest_lp2];
	  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
	  if(FlagDebug) std::cout << "debug igf:" << igf << " repid:" << repid << std::endl;
	  event.GFfitstatus_lp2[j] = (Int_t) GFTrackCont.TrackCheck(igf, repid);
	  //    HF1( 3, event.GFfitstatus[j]);
	  if(!event.GFfitstatus_lp2[j]) continue;
	  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;                	  
	  Int_t nh = GFTrackCont.GetNHits(igf);
	  event.GFnhtrack_lp2[j] = nh;
	  event.GFchisqr_lp2[j] = GFTrackCont.GetChi2NDF(igf, repid);
	  event.GFcharge_lp2[j] = GFTrackCont.GetCharge(igf, repid);
	  event.GFtof_lp2[j] = GFTrackCont.GetTrackTOF(igf, 0, -1, repid); // to be checked  3rd argument is true or not
	  event.GFtracklen_lp2[j] = GFTrackCont.GetTrackLength(igf, 0, -1, repid); // to be checked  3rd argument is true or not
	  event.GFpval_lp2[j] = GFTrackCont.GetPvalue(igf, repid);
	  event.GFpdgcode_lp2[j] = GFTrackCont.GetPDGcode(igf, repid);
	  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;                	  
	  TVector3 GFmom_decays = GFlp2_p_mom_container[gfbest_lp2];
	  if(j==1) GFmom_decays = GFlp2_pi_mom_container[gfbest_lp2];
	  if(j==2) GFmom_decays = GFlp2_p2_mom_container[gfbest_lp2];
	  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;                	  
	  Double_t GFmass_decays = GFlp2_p_mass_container[gfbest_lp2];
	  if(j==1) GFmass_decays = GFlp2_pi_mass_container[gfbest_lp2];
	  if(j==2) GFmass_decays = GFlp2_p2_mass_container[gfbest_lp2];
	  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
	  
	  event.GFdecays_m2_lp2[j] = GFmass_decays;
	  event.GFdecays_mom_lp2[j] = GFmom_decays.Mag();
	  event.GFdecays_mom_x_lp2[j] = GFmom_decays.x();
	  event.GFdecays_mom_y_lp2[j] = GFmom_decays.y();
	  event.GFdecays_mom_z_lp2[j] = GFmom_decays.z();
	  event.GFmomloss_lp2[j] = GFmom_decays.Mag() - GFTrackCont.GetMom(igf, 0, repid).Mag();
	  Double_t pdgmass[3] = {ProtonMass, PionMass, ProtonMass};
	  event.GFeloss_lp2[j] = TMath::Sqrt(GFmom_decays.Mag()*GFmom_decays.Mag() + pdgmass[j]*pdgmass[j]) - TMath::Sqrt(GFTrackCont.GetMom(igf, 0, repid).Mag()*GFTrackCont.GetMom(igf, 0, repid).Mag() + pdgmass[j]*pdgmass[j]);

	  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
	  
	  event.GFlayer_lp2[j].resize(nh);
	  event.GFrow_lp2[j].resize(nh);
	  event.GFpos_x_lp2[j].resize(nh);
	  event.GFpos_y_lp2[j].resize(nh);
	  event.GFpos_z_lp2[j].resize(nh);
	  event.GFmom_lp2[j].resize(nh);
	  event.GFmom_x_lp2[j].resize(nh);
	  event.GFmom_y_lp2[j].resize(nh);
	  event.GFmom_z_lp2[j].resize(nh);
	  event.GFresidual_x_lp2[j].resize(nh);
	  event.GFresidual_y_lp2[j].resize(nh);
	  event.GFresidual_z_lp2[j].resize(nh);
	  event.GFresidual_px_lp2[j].resize(nh);
	  event.GFresidual_py_lp2[j].resize(nh);
	  event.GFresidual_pz_lp2[j].resize(nh);
	  event.GFresidual_p_lp2[j].resize(nh);
	  event.GFresidual6D_x_lp2[j].resize(nh);
	  event.GFresidual6D_y_lp2[j].resize(nh);
	  event.GFresidual6D_z_lp2[j].resize(nh);
	  event.GFresidual6D_px_lp2[j].resize(nh);
	  event.GFresidual6D_py_lp2[j].resize(nh);
	  event.GFresidual6D_pz_lp2[j].resize(nh);
	  event.GFresolution_x_lp2[j].resize(nh);
	  event.GFresolution_y_lp2[j].resize(nh);
	  event.GFresolution_z_lp2[j].resize(nh);
	  event.GFresolution_p_lp2[j].resize(nh);
	  event.GFresolution_px_lp2[j].resize(nh);
	  event.GFresolution_py_lp2[j].resize(nh);
	  event.GFresolution_pz_lp2[j].resize(nh);
	  event.GFpull_x_lp2[j].resize(nh);
	  event.GFpull_y_lp2[j].resize(nh);
	  event.GFpull_z_lp2[j].resize(nh);
	  event.GFpull_p_lp2[j].resize(nh);
	  event.GFpull_px_lp2[j].resize(nh);
	  event.GFpull_py_lp2[j].resize(nh);
	  event.GFpull_pz_lp2[j].resize(nh);
	  //    HF1( 600, event.GFpval_lp2[j]);

	  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;                	  

	  Int_t id = GFlp2_p_id_container[gfbest_lp2];
	  if(j==1) id = GFlp2_pi_id_container[gfbest_lp2];
	  if(j==2) id = GFlp2_p2_id_container[gfbest_lp2];
	  
	  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
	  
	  Int_t ihit = 0;
	  TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( id );
	  double GFchisqrPos=0;
	  for( Int_t ih=0; ih<tp -> GetNHit(); ++ih ){
	    
	    if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
	    
	    Int_t layer = (Int_t) event.hitlayer[id][ih];
	    TPCLTrackHit *helix_point = tp -> GetHitInOrder(ih);
	    if(!helix_point -> IsGoodForTracking()) continue;

	    const TVector3 &hit0 = helix_point -> GetLocalHitPos();
	    double row = helix_point->GetMRow()+0.5;
	    TVector3 mom0 = helix_point -> GetMomentumHelix(event.GFcharge_lp2[j]);
	    TVector3 hit = GFTrackCont.GetPos(igf, ihit, repid);
	    TVector3 mom = GFTrackCont.GetMom(igf, ihit, repid);

	    Double_t residual_[5];
	    Double_t pull_[5];
	    Double_t GFresidual6D[6];
	    Double_t GFpull6D[6];
	    TVector3 dumV;
	    Double_t dumd;
	    GFTrackCont.GetTrackPull(igf, event.GFpdgcode_lp2[j], dumV,
				     dumd, mom0, hit0, residual_,
				     pull_, GFresidual6D, GFpull6D);

	    event.GFlayer_lp2[j][ihit] = layer;
	    event.GFrow_lp2[j][ihit] = row;
	    event.GFmom_x_lp2[j][ihit] = mom.x();
	    event.GFmom_y_lp2[j][ihit] = mom.y();
	    event.GFmom_z_lp2[j][ihit] = mom.z();
	    event.GFmom_lp2[j][ihit] = mom.Mag();
	    event.GFpos_x_lp2[j][ihit] = hit.x();
	    event.GFpos_y_lp2[j][ihit] = hit.y();
	    event.GFpos_z_lp2[j][ihit] = hit.z();
	    event.GFresidual_x_lp2[j][ihit] = hit.x() - hit0.x();
	    event.GFresidual_y_lp2[j][ihit] = hit.y() - hit0.y();
	    event.GFresidual_z_lp2[j][ihit] = hit.z() - hit0.z();
	    event.GFresidual_px_lp2[j][ihit] = mom.x() - mom0.x();
	    event.GFresidual_py_lp2[j][ihit] = mom.y() - mom0.y();
	    event.GFresidual_pz_lp2[j][ihit] = mom.z() - mom0.z();
	    event.GFresidual_p_lp2[j][ihit] = mom.Mag() - mom0.Mag();

	    event.GFresidual6D_x_lp2[j][ihit] =GFresidual6D[0];
	    event.GFresidual6D_y_lp2[j][ihit] =GFresidual6D[1];
	    event.GFresidual6D_z_lp2[j][ihit] =GFresidual6D[2];
	    event.GFresidual6D_px_lp2[j][ihit] =GFresidual6D[3];
	    event.GFresidual6D_py_lp2[j][ihit] =GFresidual6D[4];
	    event.GFresidual6D_pz_lp2[j][ihit] =GFresidual6D[5];

	    event.GFresolution_x_lp2[j][ihit] = GFresidual6D[0]/GFpull6D[0];
	    event.GFresolution_y_lp2[j][ihit] = GFresidual6D[1]/GFpull6D[1];
	    event.GFresolution_z_lp2[j][ihit] = GFresidual6D[2]/GFpull6D[2];
	    event.GFresolution_px_lp2[j][ihit] = GFresidual6D[3]/GFpull6D[3];
	    event.GFresolution_py_lp2[j][ihit] = GFresidual6D[4]/GFpull6D[4];
	    event.GFresolution_pz_lp2[j][ihit] = GFresidual6D[5]/GFpull6D[5];
	    event.GFpull_x_lp2[j][ihit] = GFpull6D[0];
	    event.GFpull_y_lp2[j][ihit] = GFpull6D[1];
	    event.GFpull_z_lp2[j][ihit] = GFpull6D[2];
	    event.GFpull_px_lp2[j][ihit] = GFpull6D[3];
	    event.GFpull_py_lp2[j][ihit] = GFpull6D[4];
	    event.GFpull_pz_lp2[j][ihit] = GFpull6D[5];
	    double GFresolution_t = hypot(event.GFresolution_x_lp2[j][ihit],event.GFresolution_z_lp2[j][ihit]);
	    double GFresidual_t = hypot(event.GFresidual_x_lp2[j][ihit],event.GFresidual_z_lp2[j][ihit]);
	    GFchisqrPos+= hypot(GFresidual_t/GFresolution_t,event.GFresidual_y_lp2[j][ihit]/(event.GFresolution_y_lp2[j][ihit]));

	    int hn = 2001 + j;
	    HF2(hn,hit.z(),hit.x());
	    HF1( 601, event.GFpull_x_lp2[j][ihit]);  
	    HF1( 602, event.GFpull_y_lp2[j][ihit]);
	    HF1( 603, event.GFpull_z_lp2[j][ihit]);
	    HF1( 604, event.GFpull_px_lp2[j][ihit]);
	    HF1( 605, event.GFpull_py_lp2[j][ihit]);
	    HF1( 606, event.GFpull_pz_lp2[j][ihit]);
	    if(event.GFpval_lp2[j]>0.01){
	      HF1( 611, event.GFpull_x_lp2[j][ihit]);
	      HF1( 612, event.GFpull_y_lp2[j][ihit]);
	      HF1( 613, event.GFpull_z_lp2[j][ihit]);
	      HF1( 614, event.GFpull_px_lp2[j][ihit]);
	      HF1( 615, event.GFpull_py_lp2[j][ihit]);
	      HF1( 616, event.GFpull_pz_lp2[j][ihit]);
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

	  event.GFchisqrPos_lp2[j]=GFchisqrPos;
	  event.GFpvalPos_lp2[j]=GFpvalPos;

	  if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
		
	  HF1( 10131+j, event.GFdecays_m2_lp2[j] );
	  for(int i=0; i<event.nKK; ++i){
	    Double_t bek = event.BEkaonTPC[i];
	    Double_t theta = event.thetaTPC[i];
	    if(momK18[i].Mag()<0.025) continue;      
	    if(momKurama[i].Mag()<0.025) continue;
	    Double_t missMom = (GFlp2_mom_container[gfbest_lp2] - momTrans[i])*momTrans[i]/momTransMag.at(i);	
	    HF1( 10113+j, event.GFdecays_mom_lp2[j] );
	    HF1( 11441+j*10, event.GFdecays_mom_lp2[j] );
	    for( Int_t i_bek=0; i_bek<5; i_bek++){
	      if( bek<-0.1+0.1*i_bek ){
		HF1( 11442+j*10+i_bek, event.GFdecays_mom_lp2[j] );
		break;
	      }
	    }
	  } // nKp
	} //igf Lp2
	TVector3 KPVert(event.vtx[0],event.vty[0],event.vtz[0]);
	TVector3 LP2Mom = GFlp2_mom_container[gfbest_lp2];
	TVector3 LP2Vert = GFlp2_vert_container[gfbest_lp2];
	if(event.isgoodTPC[0] == 1){
	  KPVert = TVector3(event.vtxTPC[0],event.vtyTPC[0],event.vtzTPC[0]);
	}
	for(int i=0; i<event.nKK; ++i){
	  if(momK18[i].Mag()<0.025) continue;    
	  if(momKurama[i].Mag()<0.025) continue;
	  Double_t bek = event.BEkaonTPC[i];
	  Double_t theta = event.thetaTPC[i];
	  if(momK18[i].Mag()<0.025) continue;      
	  if(momKurama[i].Mag()<0.025) continue;
	  Double_t missMom = (GFlp2_mom_container[gfbest_lp2] - momTrans[i])*momTrans[i]/momTransMag.at(i);      
	  HF1( 11101, event.GFlp2_angle );
	  HF1( 11211, event.GFlmass_lp2 );
	  HF1( 11240, event.GFlp2mass );
	  HF1( 11231, event.GFlmass_lp2 );
	  HF1( 11471, event.GFlmom_lp2 );
	  HF1( 11481, missMom );
	  HF1( 11217, event.GFlp2mass );
	  HF2( 13000, event.GFlp2_angle, event.GFlp2mom );
	  HF2( 13010, event.GFlp2_angle, event.GFlp2mass );
	  HF2( 13020, event.GFlmom_lp2, momTransMag.at(i) );
	  for(Int_t i_bek=0; i_bek<5; i_bek++){
	    if( bek<-0.1+0.1*i_bek ){
	      HF1( 11102+i_bek, event.GFlp2_angle );
	      HF1( 11218+i_bek, event.GFlp2mass );
	      HF1( 11232+i_bek, event.GFlmass_lp2 );	
	      HF1( 11472+i_bek, event.GFlmom_lp2 );
	      HF1( 11482+i_bek, missMom );	
	      HF2( 13001+i_bek, event.GFlp2_angle, event.GFlp2mom );
	      HF2( 13011+i_bek, event.GFlp2_angle, event.GFlp2mass );
	      HF2( 13021+i_bek, event.GFlp2mom, momTransMag.at(i) );
	      break;
	    }
	  }
	  HF1( 11207, bek );
	  HF1( 11250, bek );
	  HF1( 11260, momTransMag.at(i) );
	  HF1( 11270, event.GFlp2mom );
	  HF1( 11280, missMom );
	  HF2( 13040, bek, momTransMag.at(i) );
	  HF2( 13050, bek, event.GFlp2mass );
	  HF1( 11240, event.GFlp2mass );
	  if( theta<5. )  HF1( 11208, bek );
	  if( theta<10. ) HF1( 11209, bek );
	  for(Int_t i_angle=0; i_angle<9; i_angle++){
	    if( event.GFlp2_angle<20.0*(i_angle+1) ){
	      HF1( 11241+i_angle, event.GFlp2mass );
	      HF1( 11251+i_angle, bek );
	      HF1( 11261+i_angle, momTransMag.at(i) );
	      HF1( 11271+i_angle, event.GFlp2mom );
	      HF1( 11281+i_angle, missMom);	
	      HF2( 13031+i_angle, event.GFlp2mom, momTransMag.at(i) );
	      HF2( 13041+i_angle, bek, momTransMag.at(i) );
	      HF2( 13051+i_angle, bek, event.GFlp2mass );
	      break;
	    }
	  }      
	} // nKp
      } // GFlp2flag&&lp2candidate==1
      
      if(FlagDebug) std::cout << "debug " << __LINE__ << std::endl;
      
      HF1( 30, event.MissMass[0]);
      HF1( 31, event.GFlmass);
      //  HF1( 32, event.GFlpmass);
      HF1( 33, event.GFdecays_mom[0]);
      HF1( 34, event.GFdecays_mom[1]);
      //  HF1( 35, event.GFdecays_mom[2]);
      HF1( 36, event.GFlmom);
      //  HF1( 37, event.GFximom);
      HF1( 38, event.GFppi_dist);
      //  HF1( 39, event.GFlpi_dist);
      HF1( 40, event.GFdecays_m2[0]);
      HF1( 41, event.GFdecays_m2[1]);
      //  HF1( 42, event.GFdecays_m2[2]);
      HF1( 2, event.GFstatus++);
      HF1( 1, event.status++ );    
    }
  } // best!=-1

  {
    for(int it=0; it<ntTpc; it++){
      if(event.isK18[it]==1 || event.isKurama[it]==1 || event.isAccidental[it]==1 ) continue; 
      if(!GFTrackCont.IsInsideTarget(it)) continue;     
      TVector3 mom = GFTrackCont.GetMom(it, 0);
      Int_t repid=-1;
      Int_t hitid_htof; Double_t tof; Double_t len;
      TVector3 pos_htof; Double_t track2tgt_dist;
      Bool_t htofextrapolation =
        GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
                                         event.HtofSeg, event.posHtof,
                                         hitid_htof, tof,
                                         len, pos_htof, track2tgt_dist);
      if(htofextrapolation){
	Double_t beta = len / event.tHtof[hitid_htof]/MathTools::C();
	Double_t mass2 = Kinematics::MassSquare(mom.Mag(),len,event.tHtof[hitid_htof]);
	for(int ihhtof=0; ihhtof<event.nhHtof; ihhtof++){
	  HF2( 520, event.HtofSeg[ihhtof], mass2 );
	}

	for(int ikk=0; ikk<event.nKK; ++ikk){
	  Double_t bek = event.BEkaonTPC[ikk];
	  Double_t theta = event.thetaTPC[ikk];
	  if(momK18[ikk].Mag()<0.025) continue;      
	  if(momKurama[ikk].Mag()<0.025) continue;
	  if(bek>0.3) continue;
	  if(event.charge[it]<0){
	    HF2( 23000, mass2, event.dEdx[it] );
	    HF1( 500, mass2 );
	    HF1( 510, beta );	
	    for(int i=0; i<5; i++){
	      double minmom=0.1*i;
	      double maxmom=minmom+0.1;
	      if(mom.Mag()>minmom&&mom.Mag()<maxmom) HF2( 23001+i, mass2, event.dEdx[it] );
	    }
	  } else {
	    HF2( 23010, mass2, event.dEdx[it] );
	    HF1( 501, mass2 );
	    HF1( 511, beta );
	    for(int i=0; i<5; i++){
	      double minmom=0.1*i;
	      double maxmom=minmom+0.1;
	      if(mom.Mag()>minmom&&mom.Mag()<maxmom) HF2( 23011+i, mass2, event.dEdx[it] );
	    }
	  }	  
	}
      }
    }
  }
  
  // Kp->Kp reaction start
  if(FlagDebug) std::cout << "debug: start searching scat K- " << std::endl;
  //  HypTPCTask& GFTrackCont = HypTPCTask::GetInstance();  
  std::vector<Int_t> GFkm_id_container;
  std::vector<Int_t> GFkm_rep_container;
  std::vector<TVector3> GFkm_mom_container;
  std::vector<Double_t> GFkm_mass2_container;
  Int_t km_candidates=0;
  Int_t trackkp=0;
  if(FlagDebug) std::cout << "debug: line "  << __LINE__ << std::endl;  
  for(Int_t it1=0;it1<ntTpc;it1++){  
    if(event.isK18[it1]==1 || event.isKurama[it1]==1 || event.isAccidental[it1]==1 ) continue;
    if(!GFTrackCont.IsInsideTarget(it1)) continue;    
    trackkp++;
    if(event.nsigma_kaon[it1]>3 || event.nsigma_kaon[it1]<-3) continue;
    if(event.nsigma_pion[it1]<2.5) continue;
    if(event.mom0[it1]>0.5) continue;
    if(event.charge[it1]>0) continue;
    HF2(20110, event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
    if(event.charge[it1]<0){
      HF2(20111, -event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
    } else {
      HF2(20112, event.mom0[it1]*event.charge[it1], event.dEdx[it1]);
    }
    Double_t mass2 = qnan;
    Double_t inverse_beta = qnan;
    TVector3 mom(qnan, qnan, qnan);
    //    Int_t repid = -1;
    //    Int_t repid = +1;
    Int_t repid = 0;
    if(!GFTrackCont.TrackCheck(it1, repid)) continue;
    if(!GFTrackCont.IsInsideTarget(it1)) continue;
    event.GFinside[it1] = 1;
    Int_t htofhitid_km; Double_t tracklen_km;
    Double_t tof; TVector3 pos; Double_t track2tgt_dist;
    Bool_t htofextrapolation =
      GFTrackCont.TPCHTOFTrackMatching(it1, repid, tgtpos,
				       event.HtofSeg, event.posHtof,
				       htofhitid_km, tof, tracklen_km,
				       pos, track2tgt_dist);
    mom = GFTrackCont.GetMom(it1, -1, repid);
    if(htofextrapolation){
      mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_km, event.tHtof[htofhitid_km]);
      inverse_beta = MathTools::C()*event.tHtof[htofhitid_km]/tracklen_km;
    }
    if(FlagDebug) std::cout << "debug: line "  << __LINE__ << std::endl;      
    GFkm_id_container.push_back(it1);
    GFkm_rep_container.push_back(repid);
    GFkm_mom_container.push_back(mom);
    GFkm_mass2_container.push_back(mass2);    
    km_candidates++;
    if(FlagDebug) std::cout << "debug: line "  << __LINE__ << std::endl;      
  } // it1
  HF1(20011, trackkp);
  if(FlagDebug) std::cout << "debug: line "  << __LINE__ << std::endl;  
  Int_t best_km = -1; Double_t prev_massdiff_km = 9999.;
  for(Int_t candi=0;candi<km_candidates;candi++){
    Double_t diff = TMath::Abs(GFkm_mass2_container[candi] - KaonMass);
    if(prev_massdiff_km > diff){
      prev_massdiff_km = diff;
      best_km = candi;
      std::cout << "best K-: " << best_km << std::endl;
    }
  }
  if(FlagDebug) std::cout << "debug: line "  << __LINE__ << std::endl;  
  if(best_km != -1){
    Int_t it= GFkm_id_container[best_km];
    if(event.isK18[it]!=1 && event.isKurama[it]!=1 && event.isAccidental[it]!=1 ){      
      for(Int_t it2=0; it2<ntTpc; it2++){
	if(event.isK18[it2]==1 || event.isKurama[it2]==1 || event.isAccidental[it2]==1 ) continue;
	if(it2==it) continue;
	if(GFTrackCont.IsInsideTarget(it2)) continue;
	if(event.charge[it2]>0)	continue;
	//Km endpoint
	TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix(it);
	Int_t nhit = tp->GetNHit();
	TPCLTrackHit *helix_point = tp -> GetHitInOrder(nhit-1);
	const TVector3 &endhit = helix_point -> GetLocalHitPos();
	
	Int_t trackid_km = it;
	Int_t trackid_decay = it2;
	Int_t repid_km = 0; Int_t repid_decay = 0;
	Double_t GFextrapolation_decays[3];
	Double_t GFmass2_decays[3] = {qnan, qnan, qnan};
	TVector3 GFmom_decays[3]; TVector3 GFkmdecay_vert; double GFkmdecay_dist=qnan;
	if(!GFTrackCont.FindVertex(trackid_km, trackid_decay,
				   repid_km, repid_decay,
				   GFextrapolation_decays[0], GFextrapolation_decays[1],
				   GFmom_decays[0], GFmom_decays[1],
				   GFkmdecay_dist, GFkmdecay_vert,
				   vtx_scan_range)
	   || GFkmdecay_dist > GFkmdecay_distcut) continue;
	if(!((endhit-GFkmdecay_vert).Mag()<GFkmdecay_hitcut)) continue;
	event.kmdecays_id.push_back(it2);
	event.kmdecayflag=true;	
      }
      
      Double_t mom_km = GFkm_mom_container[best_km].Mag();
      Double_t mass2_km = GFkm_mass2_container[best_km];      
      HF1(21101, mom_km);
      HF1(21201, mass2_km);

      event.kmntrack = trackkp;
      event.kmflag=true;
      event.km_id.push_back(GFkm_id_container[best_km]);
      event.kmmass2 = mass2_km;
      event.kmmom = GFkm_mom_container[best_km].Mag();      
      event.kmmom_x = GFkm_mom_container[best_km].X();
      event.kmmom_y = GFkm_mom_container[best_km].Y();
      event.kmmom_z = GFkm_mom_container[best_km].Z();
      TVector3 momKmScatp[event.nKK];
      for(int i=0; i<event.nKK; ++i){
	if(momK18[i].Mag()<0.025) continue;
	if(momKurama[i].Mag()<0.025) continue;	
	momKmScatp[i] = GFkm_mom_container[best_km] + momKurama[i];
	event.kmmomdiff.push_back((momKmScatp[i]-momK18[i]).Mag());
	event.kmscatpmom.push_back(momKmScatp[i].Mag());
	event.kmscatpmom_x.push_back(momKmScatp[i].X());
	event.kmscatpmom_y.push_back(momKmScatp[i].Y());
	event.kmscatpmom_z.push_back(momKmScatp[i].Z());
	if((momKmScatp[i]-momK18[i]).Mag()<0.30) event.kmmommatch=true;	
	//	if(mass2_km<0.2 || mass2_km>0.4) continue;
	if(mass2_km<0.15 || mass2_km>0.35) continue;
	Double_t bek = event.BEkaonTPC[i];
	Double_t theta = event.thetaTPC[i];
	HF1(21001,bek);
	if(theta<5){
	  HF1(21002,bek);
	  for(int i=0; i<5; i++){
	    if(mass2_km>0.05*(i+1)) HF1( 21040+i, bek ); 
	  }	  
	}
	if(theta<10){
	  HF1(21003,bek);
	  for(int i=0; i<5; i++){
	    if(mass2_km>0.05*(i+1)) HF1( 21050+i, bek ); 
	  }
	}
	if(theta>3.5&&theta<4.5){
	  HF1(21004,bek);
	  for(int i=0; i<4; i++){
	    if(trackkp==i) HF1(21020+i,bek);
	  }
	  for(int i=0; i<5; i++){
	    if(mass2_km>0.05*(i+1)) HF1( 21030+i, bek ); 
	  }
	}
	for(int i=0; i<4; i++){
	  if(trackkp==i) HF1(21010+i,bek);
	}
      }
    }
  }
  if(FlagDebug) std::cout << "debug: line "  << __LINE__ << std::endl; 
  // Kp->Kp reaction end

  // Kbar0 search 
  if(FlagDebug) std::cout << "debug: start searching scat K0 " << std::endl;
  Int_t trackK0=0;
  Int_t k0_candidates=0;
  std::vector<Int_t> GFk0_pip_id_container, GFk0_pim_id_container;
  std::vector<Double_t> GFk0_reconL_mass_container, GFk0_recon_mass_container, GFk0_miss_mass_container;
  std::vector<Double_t> GFk0_pipidist_container;
  std::vector<Double_t> GFk0_pip_mass2_container, GFk0_pim_mass2_container;
  std::vector<TVector3> GFk0_mom_container, GFk0_pip_mom_container, GFk0_pim_mom_container, GFk0_miss_mom_container;  
  {
    for(Int_t it1=0;it1<ntTpc;it1++){
      if(event.isElectron[it1]==1) continue;
      if(event.isK18[it1]==1) continue;
      if(event.isKurama[it1]==1) continue; //correct or not?
      if(event.isBeam[it1]==1) continue;
      if(event.isAccidental[it1]==1) continue;
      //      if((event.pid[it1]&1)!=1) continue; //select pi+like      
      if(event.charge[it1]!=1) continue;
      if(event.nsigma_pion[it1]<-3 || event.nsigma_pion[it1]>3) continue;      

      Double_t mass2pip = qnan;
      Double_t inverse_beta_pip = qnan;
      TVector3 mom_pip(qnan, qnan, qnan);
      
      Int_t repid_pip = 0;
      if(!GFTrackCont.TrackCheck(it1, repid_pip)) continue;

      Int_t htofhitid_pip; Double_t tracklen_pip;
      Double_t tof_pip; TVector3 pos_pip; Double_t track2tgt_dist_pip;
      Bool_t htofextrapolation_pip =
	GFTrackCont.TPCHTOFTrackMatching(it1, repid_pip, tgtpos,
					 event.HtofSeg, event.posHtof,
					 htofhitid_pip, tof_pip, tracklen_pip,
					 pos_pip, track2tgt_dist_pip);
      mom_pip = GFTrackCont.GetMom(it1, -1, repid_pip);
      if(htofextrapolation_pip){
	mass2pip = Kinematics::MassSquare(mom_pip.Mag(), tracklen_pip, event.tHtof[htofhitid_pip]);
	inverse_beta_pip = MathTools::C()*event.tHtof[htofhitid_pip]/tracklen_pip;
      }      
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
	if(event.isElectron[it2]==1) continue;
	if(event.isK18[it2]==1) continue;
	if(event.isKurama[it2]==1) continue; //correct or not?
	if(event.isBeam[it2]==1) continue;
	if(event.isAccidental[it2]==1) continue;
	//	if((event.pid[it2]&1)!=1) continue; //select pi-like
	if(event.charge[it2]!=-1) continue;
	if(event.nsigma_pion[it2]<-3 || event.nsigma_pion[it2]>3) continue;

	Double_t mass2pim = qnan;
	Double_t inverse_beta_pim = qnan;
	TVector3 mom_pim(qnan, qnan, qnan);

	Int_t repid_pim = 0;
	if(!GFTrackCont.TrackCheck(it2, repid_pim)) continue;

	Int_t htofhitid_pim; Double_t tracklen_pim;
	Double_t tof_pim; TVector3 pos_pim; Double_t track2tgt_dist_pim;
	Bool_t htofextrapolation_pim =
	  GFTrackCont.TPCHTOFTrackMatching(it2, repid_pim, tgtpos,
					   event.HtofSeg, event.posHtof,
					   htofhitid_pim, tof_pim, tracklen_pim,
					   pos_pim, track2tgt_dist_pim);
	mom_pim = GFTrackCont.GetMom(it2, -1, repid_pim);
	if(htofextrapolation_pim){
	  mass2pim = Kinematics::MassSquare(mom_pim.Mag(), tracklen_pim, event.tHtof[htofhitid_pim]);
	  inverse_beta_pim = MathTools::C()*event.tHtof[htofhitid_pim]/tracklen_pim;
	}      	
	
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
	TVector3 pip_mom; TVector3 pim_mom; TVector3 k0_mom;
	TVector3 k0_vert = Kinematics::LambdaVertex(dMagneticField, pip_par, pim_par, pip_theta_min, pip_theta_max, pim_theta_min, pim_theta_max, pip_mom, pim_mom, k0_mom, pipi_dist);
	if(TMath::IsNaN(pipi_dist)) continue;
	k0_mom = pim_mom + pip_mom;

	TLorentzVector Lpim(pim_mom, TMath::Hypot(pim_mom.Mag(), PionMass));
	TLorentzVector Lpip(pip_mom, TMath::Hypot(pip_mom.Mag(), PionMass));
	TLorentzVector Lpipi = Lpip + Lpim;

	TLorentzVector Lpip_plike(pip_mom, TMath::Hypot(pip_mom.Mag(), ProtonMass));
	TLorentzVector Llambda = Lpip_plike + Lpim;

	if(TMath::Abs(k0_vert.x()) > 250. ||
	   TMath::Abs(k0_vert.z()) > 250. ||
	   TMath::Abs(k0_vert.y()) > 250.) continue; //Vertex cut

	Double_t pip_vertex_dist; Double_t pim_vertex_dist;
	if(!Kinematics::HelixDirection(k0_vert, pip_start, pip_end, pip_vertex_dist) ||
	   !Kinematics::HelixDirection(k0_vert, pim_start, pim_end, pim_vertex_dist)) continue;

	if(pip_vertex_dist > pi_vtx_distcut) continue;
	if(pim_vertex_dist > pi_vtx_distcut) continue;

	Double_t k0target_dist;
	TVector3 k0target_vtx = Kinematics::CalcCloseDistLambda(tgtpos,
								k0_vert,
								k0_mom,
								k0target_dist);
	if(k0target_dist > k0target_distcut) continue;
	if(pipi_dist < pipi_distcut){
	  trackK0++;	  
	  GFk0_pip_id_container.push_back(it1);
	  GFk0_pim_id_container.push_back(it2);
	  GFk0_reconL_mass_container.push_back(Llambda.M());
	  GFk0_recon_mass_container.push_back(Lpipi.M());
	  GFk0_pipidist_container.push_back(pipi_dist);
	  GFk0_mom_container.push_back(k0_mom);
	  GFk0_pip_mom_container.push_back(pip_mom);
	  GFk0_pim_mom_container.push_back(pim_mom);
	  GFk0_pip_mass2_container.push_back(mass2pip);
	  GFk0_pim_mass2_container.push_back(mass2pim);
	  k0_candidates++;
	}
      } //it2
    } //it1
  }

  HF1(31001,trackK0);
  //Select the best L mass combination
  Int_t bestk0 = -1; Double_t prev_massdiff_k0 = 9999.;
  for(Int_t candi=0;candi<k0_candidates;candi++){
    Double_t diff = TMath::Abs(GFk0_recon_mass_container[candi] - KaonMass);
    if(prev_massdiff_k0 > diff){
      prev_massdiff_k0 = diff;
      bestk0 = candi;
      std::cout << "best K0: " << bestk0 << std::endl;
    }
  }
  if(bestk0!=-1){
    Int_t itpip = GFk0_pip_id_container[bestk0];
    Int_t itpim = GFk0_pim_id_container[bestk0];    
    HF2(30100, event.mom0[itpip]*event.charge[itpip], event.dEdx[itpip]);
    HF2(30100, event.mom0[itpim]*event.charge[itpim], event.dEdx[itpim]);
    if(event.charge[itpip]<0){
      HF2(30101, -event.mom0[itpip]*event.charge[itpip], event.dEdx[itpip]);
    } else {
      HF2(30102, event.mom0[itpip]*event.charge[itpip], event.dEdx[itpip]);
    }
    if(event.charge[itpim]<0){
      HF2(30101, -event.mom0[itpim]*event.charge[itpim], event.dEdx[itpim]);
    } else {
      HF2(30102, event.mom0[itpim]*event.charge[itpim], event.dEdx[itpim]);
    }
    HF1( 33201, GFk0_pim_mass2_container[bestk0] );
    HF1( 34201, GFk0_pip_mass2_container[bestk0] );
    Double_t mass2pim=GFk0_pim_mass2_container[bestk0];
    Double_t mass2pip=GFk0_pip_mass2_container[bestk0];
    if(mass2pim>0.1 || mass2pip>0.1){
      HF1( 31101, GFk0_mom_container[bestk0].Mag() );
      HF1( 31301, GFk0_recon_mass_container[bestk0] );
      HF1( 31401, GFk0_reconL_mass_container[bestk0] );    
      for(int ikk=0; ikk<event.nKK; ikk++){
	Double_t bek = event.BEkaonTPC[ikk];
	Double_t theta = event.thetaTPC[ikk];
	if(momK18[ikk].Mag()<0.025) continue;      
	if(momKurama[ikk].Mag()<0.025) continue;
	HF1(31001,bek);
	if(theta<5){
	  HF1(31002,bek);
	}
	if(theta<10){
	  HF1(31003,bek);	    
	}
	if(theta>3.5&&theta<4.5){
	  HF1(31004,bek);
	}
	Double_t momRecoKmMag = momTrans[ikk].Mag();
	TLorentzVector LrecoKm(momTrans[ikk], TMath::Hypot(momRecoKmMag, KaonMass));
	TLorentzVector Lpipi(GFk0_mom_container[bestk0], TMath::Hypot(GFk0_mom_container[bestk0].Mag(), KaonMass));
	TLorentzVector LmissMom = LrecoKm - Lpipi;
	HF1(32101, LmissMom.P());
	HF1(32301, LmissMom.M());
      }
    }
  }
 

  // event.ncombiPipair = pipair_pip_id_container.size();
  // for(Int_t icombi=0; icombi<event.ncombiPipair; ++icombi){
  //   event.pipidPipair.push_back(pipair_pip_id_container[icombi]);
  //   event.pimidPipair.push_back(pipair_pim_id_container[icombi]);
  //   event.pipmomPipair.push_back(pipair_pip_mom_container[icombi].Mag());
  //   event.pipmomPipair_x.push_back(pipair_pip_mom_container[icombi].x());
  //   event.pipmomPipair_y.push_back(pipair_pip_mom_container[icombi].y());
  //   event.pipmomPipair_z.push_back(pipair_pip_mom_container[icombi].z());
  //   event.pimmomPipair.push_back(pipair_pim_mom_container[icombi].Mag());
  //   event.pimmomPipair_x.push_back(pipair_pim_mom_container[icombi].x());
  //   event.pimmomPipair_y.push_back(pipair_pim_mom_container[icombi].y());
  //   event.pimmomPipair_z.push_back(pipair_pim_mom_container[icombi].z());
  //   event.momPipair.push_back(pipair_mom_container[icombi].Mag());
  //   event.momPipair_x.push_back(pipair_mom_container[icombi].x());
  //   event.momPipair_y.push_back(pipair_mom_container[icombi].y());
  //   event.momPipair_z.push_back(pipair_mom_container[icombi].z());
  //   event.reconLmassPipair.push_back(pipair_reconL_mass_container[icombi]);
  //   event.reconmassPipair.push_back(pipair_recon_mass_container[icombi]);
  //   event.pipidistPipair.push_back(pipair_pipidist_container[icombi]);
  // }
  // Kbar0 search
  
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
  HB1( 1, "Status", 41, 0., 41. );
  HB1( 2, "Genfit Status", 20, 0., 20. );
  HB1( 3, "Genfit Fit Status", 2, 0., 2. );
  HB1( 4, "#Lambda candidates", 10, 0., 10. );
  HB1( 5, "#Lambdap candidates", 10, 0., 10. );    
  // HB1( 4, "Missing Mass [TPC+KURAMA]; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 200, 0.3, 1.3);
  // HB1( 5, "Missing Mass [KURAMA]; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 200, 0.3, 1.3);
  // HB1( 6, "Missing Mass [TPC+KURAMA] with MomCorr; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 200, 0.3, 1.3);
  // HB1( 7, "Missing Mass, #Lamda tagged [KURAMA]; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 200, 0.3, 1.3);
  // HB1( 8, "Missing Mass [TPC+KURAMA]; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 200, 0.3, 1.3);
  // HB1( 9, "Missing Mass [KURAMA]; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 200, 0.3, 1.3);
  // //  HB1( 10, "Missing Mass, #Xi^{-} tagged [KURAMA]; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 160, 1., 1.8);
  // HB1( 10, "Missing Mass, #Lambda{p} tagged [KURAMA]; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 200, 0.3, 1.3);
  // HB1( 11," #Lambda Invariant Mass; M_{p#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 80, 1.04, 1.2);
  // //  HB1( 12," #Xi Invariant Mass; M_{#Lambda#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 125, 1.2 ,1.45);
  // HB1( 12," #Lambda{p} Invariant Mass; M_{#Lambda{p}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1000, 1.0, 3.0);
  // HB1( 13, "p Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  // HB1( 14, "#pi Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  // HB1( 15, "p2 Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  // HB1( 16, "#Lambda Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  // //  HB1( 17, "#Xi Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  // HB1( 17, "#Lambda{p} Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  // HB1( 18, "Closest Dist. p#pi;Distance [mm]; Counts [/0.1 mm]", 1000, 0., 100.);
  // HB1( 19, "Closest Dist. #Lambda#p2;Distance [mm]; Counts [/0.1 mm]", 1000, 0., 100.);

  // HB1( 20, "Missing Mass [KURAMA]; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 200, 0.3, 1.3);
  // HB1( 21," #Lambda Invariant Mass; M_{p#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 80, 1.04, 1.2);
  // //  HB1( 22," #Xi Invariant Mass; M_{#Lambda#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 125, 1.2 ,1.45);
  // HB1( 22," #Lambda{p} Invariant Mass; M_{#Lambda#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1000, 1.0, 3.0);
  // HB1( 23, "[GenFit] p Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  // HB1( 24, "[GenFit] #pi Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  // HB1( 25, "[GenFit] p2 Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  // HB1( 26, "#Lambda Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  // //  HB1( 27, "#Xi Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  // HB1( 27, "#Lambda{p} Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  // HB1( 28, "Closest Dist. p#pi;Distance [mm]; Counts [/0.1 mm]", 1000, 0., 100.);
  // HB1( 29, "Closest Dist. p#Lambda;Distance [mm]; Counts [/0.1 mm]", 1000, 0., 100.);
  
  HB1( 30, "[GenFit] Missing Mass [KURAMA]; Missing mass [GeV/#font[12]{c}^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 200, 0.3, 1.3);
  HB1( 31, "[GenFit] #Lambda Invariant Mass; M_{p#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 80, 1.04, 1.2);
  //  HB1( 32, "[GenFit] #Xi Invariant Mass; M_{#Lambda#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 125, 1.2 ,1.45);
  HB1( 32, "[GenFit] {#Lambda}p Invariant Mass; M_{#Lambda#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1000, 1.0, 3.0);
  HB1( 33, "[GenFit] p_{#Lambda} Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 34, "[GenFit] #pi_{#Lambda} Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  //  HB1( 35, "[GenFit] #pi_{#Xi} Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  HB1( 35, "[GenFit] #pi_{{#Lambda}p} Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  HB1( 36, "[GenFit] #Lambda Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  //  HB1( 37, "[GenFit] #Xi Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 37, "[GenFit] {#Lambda}p Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 38, "[GenFit] Closest Dist. p#pi;Distance [mm]; Counts [/0.1 mm]", 1000, 0., 100.);
  //  HB1( 39, "[GenFit] Closest Dist. #Lambda#pi;Distance [mm]; Counts [/0.1 mm]", 1000, 0., 100.);
  HB1( 39, "[GenFit] Closest Dist. {#Lambda}p;Distance [mm]; Counts [/0.1 mm]", 1000, 0., 100.);
  HB1( 40, "[GenFit] p_{#Lambda} Mass; M^{2} [(GeV/#font[12]{c}^{2})^{2}]; Counts [/0.002 (GeV/#font[12]{c}^{2})^{2}]", 750, 0., 1.5);
  HB1( 41, "[GenFit] pi_{#Lambda} Mass; M^{2} [(GeV/#font[12]{c}^{2})^{2}]; Counts [/0.002 (GeV/#font[12]{c}^{2})^{2}]", 200, -0.1, 0.3);
  HB1( 42, "[GenFit] p_{#Lambda{p}} Mass; M^{2} [(GeV/#font[12]{c}^{2})^{2}]; Counts [/0.002 (GeV/#font[12]{c}^{2})^{2}]", 200, 0.1, 0.3);  

  HB1( 50, "beam momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 150, 0., 3.0);
  HB1( 51, "scattered momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 150, 0., 3.0);
  HB1( 52, "momentum transfer; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 200, 0., 4.0);
  HB2( 53, "MomTransfer:pK18; Momentum Transfer [GeV/#font[12]{c}]; K18 Momentum [GeV/#font[12]{c}]", 200, 0., 4.0, 200, 0., 4.0);
  HB2( 54, "MomTransfer:pKurama; Momentum Transfer [GeV/#font[12]{c}]; Kurama Momentum [GeV/#font[12]{c}]", 200, 0., 4.0, 200, 0., 4.0);    
  HB1( 60, "beam momentum[#Lambda]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 61, "scattered momentum[#Lambda]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 150, 0., 3.0);
  HB1( 62, "momentum transfer[#Lambda]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 150, 0., 3.0);
  HB1( 70, "beam momentum [#Lambdap]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 71, "scattered momentum [#Lambdap]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 150, 0., 3.0);
  HB1( 72, "momentum transfer [#Lambdap]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 150, 0., 3.0);

  HB2( 90, "pK18:thetaK18; Momentum[GeV/#font[12]{c}]; theta[deg.]", 200, 0., 4.0, 180, 0, 180);
  HB2( 91, "pTPCKurama:thetaTPCKurama; Momentum[GeV/#font[12]{c}]; theta[deg.]", 200, 0., 4.0, 180, 0, 180);
  HB2( 92, "MomTransfer:angleK18Kurama; Momentum[GeV/#font[12]{c}]; theta[deg.]", 200, 0., 4.0, 180, 0, 180);
  
  
  const Int_t nbinpoq = 1000;
  const Double_t minpoq = -2.0;
  const Double_t maxpoq = 2.0;
  const Int_t nbindedx = 1000;
  const Double_t mindedx = 0.;
  const Double_t maxdedx = 350.;

  
  HB1( 103, "p Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 104, "#pi Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  HB1( 105, "p2 Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  HB1( 113, "p Momentum [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 114, "#pi Momentum [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  HB1( 115, "p2 Momentum [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  //  HB1( 116, "#Lambda p Momentum [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  HB2( 120, "<dE/dx> [p1] nsigma; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 121, "<dE/dx> [p1] nsigma; -p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 122, "<dE/dx> [p1] nsigma; +p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 123, "<dE/dx> [pi] nsigma; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 124, "<dE/dx> [pi] nsigma; -p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 125, "<dE/dx> [pi] nsigma; +p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 126, "<dE/dx> [p2] nsigma; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 127, "<dE/dx> [p2] nsigma; -p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 128, "<dE/dx> [p2] nsigma; +p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);  
  
  HB2( 130, "<dE/dx> [p1,pi,p2];p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 131, "<dE/dx> [p1,pi,p2];-p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 132, "<dE/dx> [p1,pi,p2];+p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 140, "<dE/dx> [p1]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 141, "<dE/dx> [p1];-p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 142, "<dE/dx> [p1];+p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 150, "<dE/dx> [pi];p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 151, "<dE/dx> [pi];-p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 152, "<dE/dx> [pi];+p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 160, "<dE/dx> [p2];p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 161, "<dE/dx> [p2];-p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 162, "<dE/dx> [p2];+p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 170, "<dE/dx> [-BEk<-0.1(GeV)];p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 171, "<dE/dx> [-0.1<-BEk<0.(GeV)];p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 172, "<dE/dx> [0.<-BEk<0.1(GeV)];p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 173, "<dE/dx> [0.1<-BEk<0.2(GeV)];p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 174, "<dE/dx> [0.2<-BEk<0.3(GeV)];p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);

  HB2( 180, "<dE/dx> [-BEk<0.3(GeV)]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 181, "<dE/dx> [-BEk<0.3(GeV)]; -p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 182, "<dE/dx> [-BEk<0.3(GeV)]; +p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 183, "<dE/dx> [-BEk<0.3(GeV),3.5#circ<thetaTPC<4.5#circ]; p/q [GeV/#font[12]{c}];<dE/dx>[arb.]",nbinpoq,minpoq,maxpoq,nbindedx,mindedx,maxdedx);
  HB2( 184, "<dE/dx> [-BEk<0.3(GeV),3.5#circ<thetaTPC<4.5#circ]; -p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 185, "<dE/dx> [-BEk<0.3(GeV),3.5#circ<thetaTPC<4.5#circ]; +p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);  

  HB1( 200, "Closest Dist. #Lambda#p2;Distance [mm]; Counts [/0.1 mm]", 1000, 0., 1000.);
  HB1( 201, "VertexX #Lambdap2;Distance [mm]; Counts [/0.1 mm]", 1000, -500, 500);
  HB1( 202, "VertexY #Lambdap2;Distance [mm]; Counts [/0.1 mm]", 1000, -500, 500);
  HB1( 203, "VertexZ #Lambdap2;Distance [mm]; Counts [/0.1 mm]", 1000, -500, 500);
  HB1( 211, "VertexX p#pi;Distance [mm]; Counts [/0.1 mm]", 1000, -500, 500);
  HB1( 212, "VertexY p#pi;Distance [mm]; Counts [/0.1 mm]", 1000, -500, 500);
  HB1( 213, "VertexZ p#pi;Distance [mm]; Counts [/0.1 mm]", 1000, -500, 500);
  HB1( 221, "VertexX Kp;Distance [mm]; Counts [/0.1 mm]", 1000, -500, 500);
  HB1( 222, "VertexY Kp;Distance [mm]; Counts [/0.1 mm]", 1000, -500, 500);
  HB1( 223, "VertexZ Kp;Distance [mm]; Counts [/0.1 mm]", 1000, -500, 500);   

  HB2( 230, "<dE/dx> [p1,pi,p2] [#Lambdap]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 231, "<dE/dx> [p1,pi,p2] [#Lambdap];-p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 232, "<dE/dx> [p1,pi,p2] [#Lambdap];+p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 240, "<dE/dx> [p1] [#Lambdap]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 241, "<dE/dx> [p1] [#Lambdap];-p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 242, "<dE/dx> [p1] [#Lambdap];+p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 250, "<dE/dx> [pi] [#Lambdap]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 251, "<dE/dx> [pi] [#Lambdap];-p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 252, "<dE/dx> [pi] [#Lambdap];+p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 260, "<dE/dx> [p2] [#Lambdap]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 261, "<dE/dx> [p2] [#Lambdap];-p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 262, "<dE/dx> [p2] [#Lambdap];+p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);

  HB1( 500, "Mass2 [nega]; mass2[GeV]; counts", 1000, -5, 5);
  HB1( 501, "Mass2 [posi]; mass2[GeV]; counts", 1000, -5, 5);  
  HB1( 510, "beta [nega]; beta; counts", 100, 0, 1);
  HB1( 511, "beta [posi]; beta; counts", 100, 0, 1);  
  HB2( 550, "HTOF seg:m2; segment; m2[GeV]", 35, 0, 35, 1000, -5, 5); // segHtof=34

  HB1( 1001, "Binding Energy of Kaon with TPC; B_{K} [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 200, -0.5, 0.5);
  HB1( 1002, "Binding Energy of Kaon with TPC; B_{K} [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 200, -0.5, 0.5);
  HB1( 1003, "Kp scattering angle; Kp scattering angle [degree]; Counts ", 360., 0., 180.);
  HB1( 1051, "nhtrack of TPC; nhtrack ; Counts", 50, 0, 50);
  HB1( 1052, "nhtrack of TPC [p1]; nhtrack ; Counts", 50, 0, 50);
  HB1( 1053, "nhtrack of TPC [#pi]; nhtrack ; Counts", 50, 0, 50);
  HB1( 1054, "nhtrack of TPC [p2]; nhtrack ; Counts", 50, 0, 50);
  
  HB1( 1101, "#Lambda p opening angle; #Lambda p opening angle [degree]; Counts ", 360., 0., 180.);
  HB1( 1102, "#Lambda p opening angle [-BEk<-0.1(GeV)]; #Lambda p opening angle [degree]; Counts ", 360., 0., 180.);
  HB1( 1103, "#Lambda p opening angle [-0.1<-BEk<0(GeV)]; #Lambda p opening angle [degree]; Counts ", 360., 0., 180.);
  HB1( 1104, "#Lambda p opening angle [0<-BEk<0.1(GeV)]; #Lambda p opening angle [degree]; Counts ", 360., 0., 180.);
  HB1( 1105, "#Lambda p opening angle [0.1<-BEk<0.2(GeV)]; #Lambda p opening angle [degree]; Counts ", 360., 0., 180.);
  HB1( 1106, "#Lambda p opening angle [0.2<-BEk<0.3(GeV)]; #Lambda p opening angle [degree]; Counts ", 360., 0., 180.);
  HB1( 1107, "p1 theta; p1 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1108, "p1 theta [-BEk<-0.1(GeV)]; p1 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1109, "p1 theta [-0.1<-BEk<0(GeV)]; p1 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1110, "p1 theta [0<-BEk<0.1(GeV)]; p1 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1111, "p1 theta [0.1<-BEk<0.2(GeV)]; p1 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1112, "p1 theta [0.2<-BEk<0.3(GeV)]; p1 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1113, "#pi theta; #pi theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1114, "#pi theta [-BEk<-0.1(GeV)]; #pi theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1115, "#pi theta [-0.1<-BEk<0(GeV)]; #pi theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1116, "#pi theta [0<-BEk<0.1(GeV)]; #pi theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1117, "#pi theta [0.1<-BEk<0.2(GeV)]; #pi theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1118, "#pi theta [0.2<-BEk<0.3(GeV)]; #pi theta [degree]; Counts ", 360., 0., 180.);  
  HB1( 1119, "#Lambda theta; #Lambda theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1120, "#Lambda theta [-BEk<-0.1(GeV)]   ; #Lambda theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1121, "#Lambda theta [-0.1<-BEk<0(GeV)] ; #Lambda theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1122, "#Lambda theta [0<-BEk<0.1(GeV)]  ; #Lambda theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1123, "#Lambda theta [0.1<-BEk<0.2(GeV)]; #Lambda theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1124, "#Lambda theta [0.2<-BEk<0.3(GeV)]; #Lambda theta [degree]; Counts ", 360., 0., 180.);  
  HB1( 1125, "p1 theta (#Lambdap); p1 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1126, "p1 theta (#Lambdap)[-BEk<-0.1(GeV)]; p1 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1127, "p1 theta (#Lambdap)[-0.1<-BEk<0(GeV)]; p1 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1128, "p1 theta (#Lambdap)[0<-BEk<0.1(GeV)]; p1 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1129, "p1 theta (#Lambdap)[0.1<-BEk<0.2(GeV)]; p1 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1130, "p1 theta (#Lambdap)[0.2<-BEk<0.3(GeV)]; p1 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1131, "#pi theta (#Lambdap); #pi theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1132, "#pi theta (#Lambdap)[-BEk<-0.1(GeV)]; #pi theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1133, "#pi theta (#Lambdap)[-0.1<-BEk<0(GeV)]; #pi theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1134, "#pi theta (#Lambdap)[0<-BEk<0.1(GeV)]; #pi theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1135, "#pi theta (#Lambdap)[0.1<-BEk<0.2(GeV)]; #pi theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1136, "#pi theta (#Lambdap)[0.2<-BEk<0.3(GeV)]; #pi theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1137, "p2 theta (#Lambdap); p2 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1138, "p2 theta (#Lambdap)[-BEk<-0.1(GeV)]   ; p2 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1139, "p2 theta (#Lambdap)[-0.1<-BEk<0(GeV)] ; p2 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1140, "p2 theta (#Lambdap)[0<-BEk<0.1(GeV)]  ; p2 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1141, "p2 theta (#Lambdap)[0.1<-BEk<0.2(GeV)]; p2 theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1142, "p2 theta (#Lambdap)[0.2<-BEk<0.3(GeV)]; p2 theta [degree]; Counts ", 360., 0., 180.);  
  HB1( 1143, "#Lambda theta (#Lambdap); #Lambda theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1144, "#Lambda theta (#Lambdap)[-BEk<-0.1(GeV)]   ; #Lambda theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1145, "#Lambda theta (#Lambdap)[-0.1<-BEk<0(GeV)] ; #Lambda theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1146, "#Lambda theta (#Lambdap)[0<-BEk<0.1(GeV)]  ; #Lambda theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1147, "#Lambda theta (#Lambdap)[0.1<-BEk<0.2(GeV)]; #Lambda theta [degree]; Counts ", 360., 0., 180.);
  HB1( 1148, "#Lambda theta (#Lambdap)[0.2<-BEk<0.3(GeV)]; #Lambda theta [degree]; Counts ", 360., 0., 180.);
  
  // angle b/w decay particles
  HB1( 1160, "angle p0 and p1 [#Lambda]; angle [degree]; Counts", 360., 0., 180.);
  HB1( 1161, "angle p0 and #pi [#Lambda]; angle [degree]; Counts", 360., 0., 180.);
  HB1( 1162, "angle p0 and #Lambda [#Lambda]; angle [degree]; Counts", 360., 0., 180.);
  HB1( 1163, "angle p1 and #pi [#Lambda]; angle [degree]; Counts", 360., 0., 180.);
  HB1( 1164, "angle p1 and #Lambda [#Lambda]; angle [degree]; Counts", 360., 0., 180.);
  HB1( 1165, "angle #pi and #Lambda [#Lambda]; angle [degree]; Counts", 360., 0., 180.);
  HB1( 1170, "angle p0 and p1 [#Lambdap]; angle [degree]; Counts", 360., 0., 180.);
  HB1( 1171, "angle p0 and #pi [#Lambdap]; angle [degree]; Counts", 360., 0., 180.);
  HB1( 1172, "angle p0 and p2 [#Lambdap]; angle [degree]; Counts", 360., 0., 180.);
  HB1( 1173, "angle p0 and #Lambda [#Lambdap]; angle [degree]; Counts", 360., 0., 180.);
  HB1( 1174, "angle p1 and #pi [#Lambdap]; angle [degree]; Counts", 360., 0., 180.);
  HB1( 1175, "angle p1 and p2 [#Lambdap]; angle [degree]; Counts", 360., 0., 180.);
  HB1( 1176, "angle p1 and #Lambda [#Lambdap]; angle [degree]; Counts", 360., 0., 180.);
  HB1( 1177, "angle #pi and p2 [#Lambdap]; angle [degree]; Counts", 360., 0., 180.);
  HB1( 1178, "angle #pi and #Lambda [#Lambdap]; angle [degree]; Counts", 360., 0., 180.);
  HB1( 1179, "angle p2 and #Lambda [#Lambdap]; angle [degree]; Counts", 360., 0., 180.);
  HB2( 1180, "theta correlation b/w p2 and #Lambda [#Lambdap][-BEk<-0.1(GeV)]   ; angle [degree]; Counts", 360., 0., 180., 360., 0., 180.);
  HB2( 1181, "theta correlation b/w p2 and #Lambda [#Lambdap][-0.1<-BEk<0(GeV)] ; angle [degree]; Counts", 360., 0., 180., 360., 0., 180.);
  HB2( 1182, "theta correlation b/w p2 and #Lambda [#Lambdap][0<-BEk<0.1(GeV)]  ; angle [degree]; Counts", 360., 0., 180., 360., 0., 180.);
  HB2( 1183, "theta correlation b/w p2 and #Lambda [#Lambdap][0.1<-BEk<0.2(GeV)]; angle [degree]; Counts", 360., 0., 180., 360., 0., 180.);
  HB2( 1184, "theta correlation b/w p2 and #Lambda [#Lambdap][0.2<-BEk<0.3(GeV)]; angle [degree]; Counts", 360., 0., 180., 360., 0., 180.);
  HB2( 1185, "theta correlation b/w p2 and #Lambda [#Lambdap][0.2<-BEk<0.3(GeV)]; angle [degree]; Counts", 360., 0., 180., 360., 0., 180.);   

  HB1( 1190, "Binding Energy of Kaon [3.5#circ<thetaTPC<4.5#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1191, "Binding Energy of Kaon [3.5#circ<thetaTPC<4.5#circ && w/ #Lambda]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1192, "Binding Energy of Kaon [3.5#circ<thetaTPC<4.5#circ && w/ #Lambdap]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);    
  HB1( 1200, "K p scattering angle; Kp scattering angle [degree]; Counts ", 360., 0., 180.);
  HB1( 1201, "Binding Energy of Kaon ; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1202, "Binding Energy of Kaon [thetaTPC<5#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1203, "Binding Energy of Kaon [thetaTPC<10#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1204, "Binding Energy of Kaon [w/ #Lambda]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1205, "Binding Energy of Kaon [thetaTPC<5#circ && w/ #Lambda]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1206, "Binding Energy of Kaon [thetaTPC<10#circ && w/ #Lambda]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1207, "Binding Energy of Kaon [w/ #Lambdap]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);  
  HB1( 1208, "Binding Energy of Kaon [thetaTPC<5#circ && w/ #Lambdap]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1209, "Binding Energy of Kaon [thetaTPC<10#circ && w/ #Lambdap]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1211, "#Lambda Invariant Mass ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 1212, "#Lambda Invariant Mass [-BEk<-0.1(GeV)]   ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 1213, "#Lambda Invariant Mass [-0.1<-BEk<0(GeV)] ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 1214, "#Lambda Invariant Mass [0<-BEk<0.1(GeV)]  ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 1215, "#Lambda Invariant Mass [0.1<-BEk<0.2(GeV)]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 1216, "#Lambda Invariant Mass [0.2<-BEk<0.3(GeV)]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 1217, "#LambdaP Invariant Mass ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 1218, "#LambdaP Invariant Mass [-BEk<-0.1(GeV)]   ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 1219, "#LambdaP Invariant Mass [-0.1<-BEk<0(GeV)] ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 1220, "#LambdaP Invariant Mass [0<-BEk<0.1(GeV)]  ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 1221, "#LambdaP Invariant Mass [0.1<-BEk<0.2(GeV)]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 1222, "#LambdaP Invariant Mass [0.2<-BEk<0.3(GeV)]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 1231, "#Lambda Invariant Mass [one more p]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 1232, "#Lambda Invariant Mass [one more p][-BEk<-0.1(GeV)]   ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 1233, "#Lambda Invariant Mass [one more p][-0.1<-BEk<0(GeV)] ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 1234, "#Lambda Invariant Mass [one more p][0<-BEk<0.1(GeV)]  ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 1235, "#Lambda Invariant Mass [one more p][0.1<-BEk<0.2(GeV)]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 1236, "#Lambda Invariant Mass [one more p][0.2<-BEk<0.3(GeV)]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  
  HB1( 1240, "#LambdaP Invariant Mass ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 1241, "#LambdaP Invariant Mass [0#circ<LPOpAngle<20#circ]   ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 1242, "#LambdaP Invariant Mass [20#circ<LPOpAngle<40#circ]  ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 1243, "#LambdaP Invariant Mass [40#circ<LPOpAngle<60#circ]  ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 1244, "#LambdaP Invariant Mass [60#circ<LPOpAngle<80#circ]  ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 1245, "#LambdaP Invariant Mass [80#circ<LPOpAngle<100#circ] ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 1246, "#LambdaP Invariant Mass [100#circ<LPOpAngle<120#circ]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 1247, "#LambdaP Invariant Mass [120#circ<LPOpAngle<140#circ]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 1248, "#LambdaP Invariant Mass [140#circ<LPOpAngle<160#circ]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 1249, "#LambdaP Invariant Mass [160#circ<LPOpAngle<180#circ]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);    
  
  HB1( 1250, "Binding Energy of Kaon ; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1251, "Binding Energy of Kaon [0#circ<LPOpAngle<20#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1252, "Binding Energy of Kaon [20#circ<LPOpAngle<40#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1253, "Binding Energy of Kaon [40#circ<LPOpAngle<60#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1254, "Binding Energy of Kaon [60#circ<LPOpAngle<80#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1255, "Binding Energy of Kaon [80#circ<LPOpAngle<100#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1256, "Binding Energy of Kaon [100#circ<LPOpAngle<120#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1257, "Binding Energy of Kaon [120#circ<LPOpAngle<140#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1258, "Binding Energy of Kaon [140#circ<LPOpAngle<160#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1259, "Binding Energy of Kaon [160#circ<LPOpAngle<180#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1260, "mom transfer ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 200, -2.0, 2.0);
  HB1( 1261, "mom transfer [0#circ<LPOpAngle<20#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 200, -2.0, 2.0);
  HB1( 1262, "mom transfer [20#circ<LPOpAngle<40#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 200, -2.0, 2.0);
  HB1( 1263, "mom transfer [40#circ<LPOpAngle<60#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 200, -2.0, 2.0);
  HB1( 1264, "mom transfer [60#circ<LPOpAngle<80#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 200, -2.0, 2.0);
  HB1( 1265, "mom transfer [80#circ<LPOpAngle<100#circ];  Mom [GeV/#font[12]{c}];  Counts [/0.02 GeV/#font[12]{c}]", 200, -2.0, 2.0);
  HB1( 1266, "mom transfer [100#circ<LPOpAngle<120#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 200, -2.0, 2.0);
  HB1( 1267, "mom transfer [120#circ<LPOpAngle<140#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 200, -2.0, 2.0);
  HB1( 1268, "mom transfer [140#circ<LPOpAngle<160#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 200, -2.0, 2.0);
  HB1( 1269, "mom transfer [160#circ<LPOpAngle<180#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 200, -2.0, 2.0);  
  HB1( 1270, "Lp Mom; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 1271, "Lp Mom [0#circ<LPOpAngle<20#circ]   ;  Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 1272, "Lp Mom [20#circ<LPOpAngle<40#circ]  ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 1273, "Lp Mom [40#circ<LPOpAngle<60#circ]  ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 1274, "Lp Mom [60#circ<LPOpAngle<80#circ]  ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 1275, "Lp Mom [80#circ<LPOpAngle<100#circ] ;  Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 1276, "Lp Mom [100#circ<LPOpAngle<120#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 1277, "Lp Mom [120#circ<LPOpAngle<140#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 1278, "Lp Mom [140#circ<LPOpAngle<160#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 1279, "Lp Mom [160#circ<LPOpAngle<180#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 1280, "MissMom; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 300, -2.0, 4.0); // missMom = (momTrans-momLp)*unitvec_momTrans
  HB1( 1281, "MissMom [0#circ<LPOpAngle<20#circ]   ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 1282, "MissMom [20#circ<LPOpAngle<40#circ]  ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 1283, "MissMom [40#circ<LPOpAngle<60#circ]  ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);    
  HB1( 1284, "MissMom [60#circ<LPOpAngle<80#circ]  ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 1285, "MissMom [80#circ<LPOpAngle<100#circ] ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 1286, "MissMom [100#circ<LPOpAngle<120#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 1287, "MissMom [120#circ<LPOpAngle<140#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 1288, "MissMom [140#circ<LPOpAngle<160#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 1289, "MissMom [160#circ<LPOpAngle<180#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);

  HB1( 1301, "Momentum Transfer; Momentum Transfer [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1302, "#Lambda p Momentum [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1310, "Binding Energy of Kaon ; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1311, "Binding Energy of Kaon [thetaTPC<5#circ][theta_{#Lambda}<20#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1312, "Binding Energy of Kaon [thetaTPC<5#circ][20#circ<theta_{#Lambda}<100#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1313, "Binding Energy of Kaon [thetaTPC<5#circ][100#circ<theta_{#Lambda}<180#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1314, "Binding Energy of Kaon [thetaTPC<5#circ][theta_p<20#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1315, "Binding Energy of Kaon [thetaTPC<5#circ][20#circ<theta_{p}<90#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1316, "Binding Energy of Kaon [thetaTPC<5#circ][90#circ<theta_{p}<160#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 1317, "Binding Energy of Kaon [thetaTPC<5#circ][160#circ<theta_{p}<180#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  
  HB1( 1401, "p Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1402, "p Momentum [#Lambda decay] [-BEk<-0.1(GeV)]   ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1403, "p Momentum [#Lambda decay] [-0.1<-BEk<0(GeV)] ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1404, "p Momentum [#Lambda decay] [0<-BEk<0.1(GeV)]  ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1405, "p Momentum [#Lambda decay] [0.1<-BEk<0.2(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1406, "p Momentum [#Lambda decay] [0.2<-BEk<0.3(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);      
  HB1( 1411, "#pi Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1412, "#pi Momentum [#Lambda decay] [-BEk<-0.1(GeV)]   ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1413, "#pi Momentum [#Lambda decay] [-0.1<-BEk<0(GeV)] ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1414, "#pi Momentum [#Lambda decay] [0<-BEk<0.1(GeV)]  ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1415, "#pi Momentum [#Lambda decay] [0.1<-BEk<0.2(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1416, "#pi Momentum [#Lambda decay] [0.2<-BEk<0.3(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);      
  HB1( 1431, "#Lambda Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1432, "#Lambda Momentum [#Lambda decay] [-BEk<-0.1(GeV)]   ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1433, "#Lambda Momentum [#Lambda decay] [-0.1<-BEk<0(GeV)] ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1434, "#Lambda Momentum [#Lambda decay] [0<-BEk<0.1(GeV)]  ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1435, "#Lambda Momentum [#Lambda decay] [0.1<-BEk<0.2(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1436, "#Lambda Momentum [#Lambda decay] [0.2<-BEk<0.3(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);      
  HB1( 1441, "p Momentum [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1442, "p Momentum [#Lambdap decay] [-BEk<-0.1(GeV)]   ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1443, "p Momentum [#Lambdap decay] [-0.1<-BEk<0(GeV)] ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1444, "p Momentum [#Lambdap decay] [0<-BEk<0.1(GeV)]  ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1445, "p Momentum [#Lambdap decay] [0.1<-BEk<0.2(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1446, "p Momentum [#Lambdap decay] [0.2<-BEk<0.3(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);      
  HB1( 1451, "#pi Momentum [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1452, "#pi Momentum [#Lambdap decay] [-BEk<-0.1(GeV)]   ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1453, "#pi Momentum [#Lambdap decay] [-0.1<-BEk<0(GeV)] ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1454, "#pi Momentum [#Lambdap decay] [0<-BEk<0.1(GeV)]  ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1455, "#pi Momentum [#Lambdap decay] [0.1<-BEk<0.2(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1456, "#pi Momentum [#Lambdap decay] [0.2<-BEk<0.3(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);  
  HB1( 1461, "p2 Momentum [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1462, "p2 Momentum [#Lambdap decay] [-BEk<-0.1(GeV)]   ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1463, "p2 Momentum [#Lambdap decay] [-0.1<-BEk<0(GeV)] ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1464, "p2 Momentum [#Lambdap decay] [0<-BEk<0.1(GeV)]  ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1465, "p2 Momentum [#Lambdap decay] [0.1<-BEk<0.2(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1466, "p2 Momentum [#Lambdap decay] [0.2<-BEk<0.3(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1471, "#Lambda Momentum [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1472, "#Lambda Momentum [#Lambdap decay] [-BEk<-0.1(GeV)]   ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1473, "#Lambda Momentum [#Lambdap decay] [-0.1<-BEk<0(GeV)] ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1474, "#Lambda Momentum [#Lambdap decay] [0<-BEk<0.1(GeV)]  ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1475, "#Lambda Momentum [#Lambdap decay] [0.1<-BEk<0.2(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1476, "#Lambda Momentum [#Lambdap decay] [0.2<-BEk<0.3(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 1481, "MissMom [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 200, -2.0, 2.0); // missMom = (momTrans-momLp)*unitvec_momTrans
  HB1( 1482, "MissMom [#Lambdap decay] [-BEk<-0.1(GeV)]   ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 1483, "MissMom [#Lambdap decay] [-0.1<-BEk<0(GeV)] ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);  
  HB1( 1484, "MissMom [#Lambdap decay] [0<-BEk<0.1(GeV)]  ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 1485, "MissMom [#Lambdap decay] [0.1<-BEk<0.2(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 1486, "MissMom [#Lambdap decay] [0.2<-BEk<0.3(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  
  // vtx
  //  HB1(1301, " ; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  // - thetaTPC
  // - Vtx xyz
  // - Vtx xyz w/ thetaTPC<5
  // - Vtx xyz w/ thetaTPC<10
  // - Vtx xyz                w/ Lambda
  // - Vtx xyz w/ thetaTPC<5  w/ Lambda
  // - Vtx xyz w/ thetaTPC<10 w/ Lambda
  // 2D

  // PPi2OpeningAngle:mom of p1,pi,L
  HB2( 2001, "PPiOpAngle:MomP1; PPiOpAngle[degree]; Momentum p1  [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);
  HB2( 2002, "PPiOpAngle:MomPi; PPiOpAngle[degree]; Momentum #pi  [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);
  HB2( 2003, "PPiOpAngle:MomL; PPiOpAngle[degree]; Momentum #Lambda [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);
  HB2( 2011, "PPiOpAngle:BEkaonTPC; LPOpAngle[degree]; Binding Energy [GeV/#font[12]{c}^2]", 360, 0.0, 180., 200, -0.5, 0.5);
  HB2( 2012, "PPiOpAngle:LambdaInvariantMass; PPiOpAngle[degree]; Invariant Mass #Lambda [/0.002 GeV/#font[12]{c}^{2}]", 360, 0.0, 180., 180, 1.04, 1.4);
  // LP2OpeningAngle:mom of p1,pi,p2,L  
  HB2( 2101, "LPOpAngle:MomP1; LPOpAngle[degree]; Momentum p1  [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);  
  HB2( 2102, "LPOpAngle:MomPi; LPOpAngle[degree]; Momentum #pi [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);  
  HB2( 2103, "LPOpAngle:MomP2; LPOpAngle[degree]; Momentum p2  [GeV/#font[12]{c}]", 360, 0.0, 180.,  100, 0., 2.0);  
  HB2( 2104, "LPOpAngle:MomL ; LPOpAngle[degree]; Momentum #Lambda [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);
  HB2( 2105, "LPOpAngle:MomLp ; LPOpAngle[degree]; Momentum #Lambdap [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);  
  HB2( 2111, "LPOpAngle:BEkaonTPC; LPOpAngle[degree]; Binding Energy [GeV/#font[12]{c}^2]", 360, 0.0, 180., 200, -0.5, 0.5);
  HB2( 2112, "LPOpAngle:LambdaInvariantMass; PPiOpAngle[degree]; Invariant Mass #Lambda [/0.002 GeV/#font[12]{c}^{2}]", 360, 0.0, 180., 180, 1.04, 1.4);
  HB2( 2113, "LPOpAngle:LPInvariantMass; LPOpAngle[degree]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 360, 0.0, 180., 1750, 2.0, 5.5);
  HB2( 2114, "LPOpAngle:LPMom; LPOpAngle[degree]; Mom #Lambdap [GeV/#font[12]{c}^2]", 360, 0.0, 180., 150, 0., 3.0);  
    
  // Mom p1,pi,p2,L : mom of p1,pi,p2,L
  HB2( 2201, "MomCorr p0 and p0 [#Lambda tagged]; Mom p0[GeV/#font[12]{c}]; Mom p0[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2202, "MomCorr p0 and p1 [#Lambda tagged]; Mom p0[GeV/#font[12]{c}]; Mom p1[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0); 
  HB2( 2203, "MomCorr p0 and #pi [#Lambda tagged]; Mom p0[GeV/#font[12]{c}]; Mom #pi[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2204, "MomCorr p0 and #Lambda [#Lambda tagged]; Mom p0[MeV/#font[12]{c}]; Mom #Lambda[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2205, "MomCorr p0 and p2 [#Lambda tagged]; Mom p0[GeV/#font[12]{c}]; Mom p2[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2211, "MomCorr p1 and p1 [#Lambda tagged]; Mom p1[GeV/#font[12]{c}]; Mom p1[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2212, "MomCorr p1 and #pi [#Lambda tagged]; Mom p1[GeV/#font[12]{c}]; Mom #pi[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);  
  HB2( 2213, "MomCorr p1 and #Lambda [#Lambda tagged]; Mom p1[GeV/#font[12]{c}]; Mom #Lambda[GeV/#font[12]{c}]",100, 0., 2.0, 100, 0., 2.0);
  HB2( 2214, "MomCorr p1 and p2 [#Lambda tagged]; Mom p1[MeV/#font[12]{c}]; Mom p2[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2221, "MomCorr #pi and pi [#Lambda tagged]; Mom #pi[GeV/#font[12]{c}]; Mom #pi[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2222, "MomCorr #pi and #Lambda [#Lambda tagged]; Mom #pi[GeV/#font[12]{c}]; Mom #Lambda[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2223, "MomCorr #pi and p2 [#Lambda tagged]; Mom #pi[MeV/#font[12]{c}]; Mom p2[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2231, "MomCorr #Lambda and #Lambda [#Lambda tagged]; Mom #Lambda[GeV/#font[12]{c}]; Mom #Lambda[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2232, "MomCorr #Lambda and p2 [#Lambda tagged]; Mom #Lambda[MeV/#font[12]{c}]; Mom p2[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2241, "MomCorr p2 and p2 [#Lambda tagged]; Mom p2[MeV/#font[12]{c}]; Mom p2[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2251, "MomCorr p0 and p0 [#Lambdap tagged]; Mom p0[GeV/#font[12]{c}]; Mom p0[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2252, "MomCorr p0 and p1 [#Lambdap tagged]; Mom p0[GeV/#font[12]{c}]; Mom p1[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2253, "MomCorr p0 and #pi [#Lambdap tagged]; Mom p0[GeV/#font[12]{c}]; Mom #pi[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2254, "MomCorr p0 and #Lambda [#Lambdap tagged]; Mom p0[MeV/#font[12]{c}]; Mom #Lambda[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2255, "MomCorr p0 and p2 [#Lambdap tagged]; Mom p0[GeV/#font[12]{c}]; Mom p2[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2261, "MomCorr p1 and p1 [#Lambdap tagged]; Mom p1[GeV/#font[12]{c}]; Mom p1[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2262, "MomCorr p1 and #pi [#Lambdap tagged]; Mom p1[GeV/#font[12]{c}]; Mom #pi[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);  
  HB2( 2263, "MomCorr p1 and #Lambda [#Lambdap tagged]; Mom p1[GeV/#font[12]{c}]; Mom #Lambda[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2264, "MomCorr p1 and p2 [#Lambdap tagged]; Mom p1[MeV/#font[12]{c}]; Mom p2[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2271, "MomCorr #pi and pi [#Lambdap tagged]; Mom #pi[GeV/#font[12]{c}]; Mom #pi[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2272, "MomCorr #pi and #Lambda [#Lambdap tagged]; Mom #pi[GeV/#font[12]{c}]; Mom #Lambda[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2273, "MomCorr #pi and p2 [#Lambdap tagged]; Mom #pi[MeV/#font[12]{c}]; Mom p2[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2281, "MomCorr #Lambda and #Lambda [#Lambdap tagged]; Mom #Lambda[GeV/#font[12]{c}]; Mom #Lambda[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2282, "MomCorr #Lambda and p2 [#Lambdap tagged]; Mom #Lambda[MeV/#font[12]{c}]; Mom p2[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 2291, "MomCorr p2 and p2 [#Lambdap tagged]; Mom p2[MeV/#font[12]{c}]; Mom p2[GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  
  HB2( 2301, "MomCorr; Mom #Lambdap[GeV/#font[12]{c}]; Mom Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);

  HB2( 2401, "TheataCorr p0 and p0 [#Lambda tagged]; theta p0 [degree]; theta p0 [degree]", 360, 0., 180, 360, 0., 180);
  HB2( 2402, "TheataCorr p0 and p1 [#Lambda tagged]; theta p0 [degree]; theta p1 [degree]", 360, 0., 180, 360, 0., 180);  
  HB2( 2403, "TheataCorr p0 and #pi [#Lambda tagged]; theta p0 [degree]; theta #pi [degree]", 360, 0., 180, 360, 0., 180);  
  HB2( 2404, "TheataCorr p0 and #Lambda [#Lambda tagged]; theta p0 [degree]; theta #Lambda [degree]", 360, 0., 180, 360, 0., 180);
  HB2( 2405, "TheataCorr p0 and p2 [#Lambda tagged]; theta p0 [degree]; theta p2 [degree]", 360, 0., 180, 360, 0., 180);
  HB2( 2411, "TheataCorr p1 and p1 [#Lambda tagged]; theta p1 [degree]; theta p1 [degree]", 360, 0., 180, 360, 0., 180);  
  HB2( 2412, "TheataCorr p1 and #pi [#Lambda tagged]; theta p1 [degree]; theta #pi [degree]", 360, 0., 180, 360, 0., 180);  
  HB2( 2413, "TheataCorr p1 and #Lambda [#Lambda tagged]; theta p1 [degree]; theta #Lambda [degree]", 360, 0., 180, 360, 0., 180);
  HB2( 2414, "TheataCorr p1 and p2 [#Lambda tagged]; theta p1 [degree]; theta p2 [degree]", 360, 0., 180, 360, 0., 180);    
  HB2( 2421, "TheataCorr #pi and #pi [#Lambda tagged]; theta #pi [degree]; theta #pi [degree]", 360, 0., 180, 360, 0., 180);  
  HB2( 2422, "TheataCorr #pi and #Lambda [#Lambda tagged]; theta #pi [degree]; theta #Lambda [degree]", 360, 0., 180, 360, 0., 180);
  HB2( 2423, "TheataCorr #pi and p2 [#Lambda tagged]; theta #pi [degree]; theta p2 [degree]", 360, 0., 180, 360, 0., 180);    
  HB2( 2431, "TheataCorr #Lambda and #Lambda [#Lambda tagged]; theta #Lambda [degree]; theta #Lambda [degree]", 360, 0., 180, 360, 0., 180);
  HB2( 2432, "TheataCorr #Lambda and p2 [#Lambda tagged]; theta #Lambda [degree]; theta p2 [degree]", 360, 0., 180, 360, 0., 180);    
  HB2( 2441, "TheataCorr p2 and p2 [#Lambda tagged]; theta p2 [degree]; theta p2 [degree]", 360, 0., 180, 360, 0., 180);
  HB2( 2451, "TheataCorr p0 and p0 [#Lambdap tagged]; theta p0 [degree]; theta p0 [degree]", 360, 0., 180, 360, 0., 180);
  HB2( 2452, "TheataCorr p0 and p1 [#Lambdap tagged]; theta p0 [degree]; theta p1 [degree]", 360, 0., 180, 360, 0., 180);  
  HB2( 2453, "TheataCorr p0 and #pi [#Lambdap tagged]; theta p0 [degree]; theta #pi [degree]", 360, 0., 180, 360, 0., 180);  
  HB2( 2454, "TheataCorr p0 and #Lambda [#Lambdap tagged]; theta p0 [degree]; theta #Lambda [degree]", 360, 0., 180, 360, 0., 180);
  HB2( 2455, "TheataCorr p0 and p2 [#Lambdap tagged]; theta p0 [degree]; theta p2 [degree]", 360, 0., 180, 360, 0., 180);
  HB2( 2461, "TheataCorr p1 and p1 [#Lambdap tagged]; theta p1 [degree]; theta p1 [degree]", 360, 0., 180, 360, 0., 180);  
  HB2( 2462, "TheataCorr p1 and #pi [#Lambdap tagged]; theta p1 [degree]; theta #pi [degree]", 360, 0., 180, 360, 0., 180);  
  HB2( 2463, "TheataCorr p1 and #Lambda [#Lambdap tagged]; theta p1 [degree]; theta #Lambda [degree]", 360, 0., 180, 360, 0., 180);
  HB2( 2464, "TheataCorr p1 and p2 [#Lambdap tagged]; theta p1 [degree]; theta p2 [degree]", 360, 0., 180, 360, 0., 180);    
  HB2( 2471, "TheataCorr #pi and #pi [#Lambdap tagged]; theta #pi [degree]; theta #pi [degree]", 360, 0., 180, 360, 0., 180);  
  HB2( 2472, "TheataCorr #pi and #Lambda [#Lambdap tagged]; theta #pi [degree]; theta #Lambda [degree]", 360, 0., 180, 360, 0., 180);
  HB2( 2473, "TheataCorr #pi and p2 [#Lambdap tagged]; theta #pi [degree]; theta p2 [degree]", 360, 0., 180, 360, 0., 180);    
  HB2( 2481, "TheataCorr #Lambda and #Lambda [#Lambdap tagged]; theta #Lambda [degree]; theta #Lambda [degree]", 360, 0., 180, 360, 0., 180);
  HB2( 2482, "TheataCorr #Lambda and p2 [#Lambdap tagged]; theta #Lambda [degree]; theta p2 [degree]", 360, 0., 180, 360, 0., 180);    
  HB2( 2491, "TheataCorr p2 and p2 [#Lambdap tagged]; theta p2 [degree]; theta p2 [degree]", 360, 0., 180, 360, 0., 180);
  
  HB2( 2501, "Theata:nhtrack [#Lambda tagged]; theta p1 [degree]; nhtrack ", 360, 0., 180, 50, 0, 50);
  HB2( 2502, "Theata:nhtrack [#Lambda tagged]; theta #pi [degree]; nhtrack ", 360, 0., 180, 50, 0, 50);
  HB2( 2511, "Theata:nhtrack [#Lambdap tagged]; theta p1 [degree]; nhtrack ", 360, 0., 180, 50, 0, 50);
  HB2( 2512, "Theata:nhtrack [#Lambdap tagged]; theta #pi [degree]; nhtrack ", 360, 0., 180, 50, 0, 50);
  HB2( 2513, "Theata:nhtrack [#Lambdap tagged]; theta p2 [degree]; nhtrack ", 360, 0., 180, 50, 0, 50);
  
  HB2( 2601, "Theata:Mom p1 [#Lambda tagged]; #theta [degree]; Mom [GeV/#font[12]{c}]", 360, 0., 180, 100, 0, 2.0);
  HB2( 2602, "Theata:Mom #pi [#Lambda tagged]; #theta [degree]; Mom [GeV/#font[12]{c}]", 360, 0., 180, 100, 0, 2.0);
  HB2( 2603, "Theata:Mom #Lambda [#Lambda tagged]; #theta [degree]; Mom [GeV/#font[12]{c}]", 360, 0., 180, 100, 0, 2.0);
  HB2( 2604, "Theata:Mom p0 [#Lambda tagged]; #theta [degree]; Mom [GeV/#font[12]{c}]", 360, 0., 180, 100, 0, 2.0);
  HB2( 2611, "Theata:Mom p1 [#Lambdap tagged]; #theta [degree]; Mom [GeV/#font[12]{c}]", 360, 0., 180, 100, 0, 2.0);
  HB2( 2612, "Theata:Mom #pi [#Lambdap tagged]; #theta [degree]; Mom [GeV/#font[12]{c}]", 360, 0., 180, 100, 0, 2.0);
  HB2( 2613, "Theata:Mom p2 [#Lambdap tagged]; #theta [degree]; Mom [GeV/#font[12]{c}]", 360, 0., 180, 100, 0, 2.0);
  HB2( 2614, "Theata:Mom #Lambda [#Lambdap tagged]; #theta [degree]; Mom [GeV/#font[12]{c}]", 360, 0., 180, 100, 0, 2.0);
  HB2( 2615, "Theata:Mom p0 [#Lambdap tagged]; #theta [degree]; Mom [GeV/#font[12]{c}]", 360, 0., 180, 100, 0, 2.0);

  HB2( 3000, "LPOpAngle:MomLp ; LPOpAngle[degree]; Momentum #Lambdap [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);    
  HB2( 3001, "LPOpAngle:MomLp [-BEk<-0.1(GeV)]   ; LPOpAngle[degree]; Momentum #Lambdap [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);
  HB2( 3002, "LPOpAngle:MomLp [-0.1<-BEk<0(GeV)] ; LPOpAngle[degree]; Momentum #Lambdap [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);
  HB2( 3003, "LPOpAngle:MomLp [0<-BEk<0.1(GeV)]  ; LPOpAngle[degree]; Momentum #Lambdap [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);
  HB2( 3004, "LPOpAngle:MomLp [0.1<-BEk<0.2(GeV)]; LPOpAngle[degree]; Momentum #Lambdap [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);
  HB2( 3005, "LPOpAngle:MomLp [0.2<-BEk<0.3(GeV)]; LPOpAngle[degree]; Momentum #Lambdap [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);
  
  HB2( 3010, "LPOpAngle:LPInvMass ; LPOpAngle[degree]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]",360, 0.0, 180., 1750, 2.0, 5.5);
  HB2( 3011, "LPOpAngle:LPInvMass [-BEk<-0.1(GeV)]   ; LPOpAngle[degree]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]",360, 0.0, 180., 1750, 2.0, 5.5);
  HB2( 3012, "LPOpAngle:LPInvMass [-0.1<-BEk<0(GeV)] ; LPOpAngle[degree]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]",360, 0.0, 180., 1750, 2.0, 5.5);
  HB2( 3013, "LPOpAngle:LPInvMass [0<-BEk<0.1(GeV)]  ; LPOpAngle[degree]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]",360, 0.0, 180., 1750, 2.0, 5.5);
  HB2( 3014, "LPOpAngle:LPInvMass [0.1<-BEk<0.2(GeV)]; LPOpAngle[degree]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]",360, 0.0, 180., 1750, 2.0, 5.5);
  HB2( 3015, "LPOpAngle:LPInvMass [0.2<-BEk<0.3(GeV)]; LPOpAngle[degree]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]",360, 0.0, 180., 1750, 2.0, 5.5);
  
  HB2( 3020, "MomTransfer:MomLp ; Momentum #Lambdap [GeV/#font[12]{c}] ;Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 3021, "MomTransfer:MomLp [-BEk<-0.1(GeV)]  ; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 3022, "MomTransfer:MomLp [-0.1<-BEk<0(GeV)]; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);  
  HB2( 3023, "MomTransfer:MomLp [0<-BEk<0.1(GeV)] ; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 3024, "MomTransfer:MomLp [0.1<-BEk<0.2(GeV)]; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 3025, "MomTransfer:MomLp [0.2<-BEk<0.3(GeV)]; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 3031, "MomTransfer:MomLp [0#circ<LPOpAngle<20#circ]   ; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 3032, "MomTransfer:MomLp [20#circ<LPOpAngle<40#circ]  ; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 3033, "MomTransfer:MomLp [40#circ<LPOpAngle<60#circ]  ; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 3034, "MomTransfer:MomLp [60#circ<LPOpAngle<80#circ]  ; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 3035, "MomTransfer:MomLp [80#circ<LPOpAngle<100#circ] ; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 3036, "MomTransfer:MomLp [100#circ<LPOpAngle<120#circ]; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 3037, "MomTransfer:MomLp [120#circ<LPOpAngle<140#circ]; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 3038, "MomTransfer:MomLp [140#circ<LPOpAngle<160#circ]; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);  
  HB2( 3039, "MomTransfer:MomLp [160#circ<LPOpAngle<180#circ]; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);

  HB2( 3040, "BE_{K}:MomTransfer; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 3041, "BE_{K}:MomTransfer [0#circ<LPOpAngle<20#circ]   ; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 3042, "BE_{K}:MomTransfer [20#circ<LPOpAngle<40#circ]  ; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 3043, "BE_{K}:MomTransfer [40#circ<LPOpAngle<60#circ]  ; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 3044, "BE_{K}:MomTransfer [60#circ<LPOpAngle<80#circ]  ; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 3045, "BE_{K}:MomTransfer [80#circ<LPOpAngle<100#circ] ; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 3046, "BE_{K}:MomTransfer [100#circ<LPOpAngle<120#circ]; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 3047, "BE_{K}:MomTransfer [120#circ<LPOpAngle<140#circ]; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 3048, "BE_{K}:MomTransfer [140#circ<LPOpAngle<160#circ]; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 3049, "BE_{K}:MomTransfer [160#circ<LPOpAngle<180#circ]; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);

  HB2( 3050, "BE_{K}:LPInvMass; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 3051, "BE_{K}:LPInvMass [0#circ<LPOpAngle<20#circ]   ; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 3052, "BE_{K}:LPInvMass [20#circ<LPOpAngle<40#circ]  ; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 3053, "BE_{K}:LPInvMass [40#circ<LPOpAngle<60#circ]  ; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 3054, "BE_{K}:LPInvMass [60#circ<LPOpAngle<80#circ]  ; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 3055, "BE_{K}:LPInvMass [80#circ<LPOpAngle<100#circ] ; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 3056, "BE_{K}:LPInvMass [100#circ<LPOpAngle<120#circ]; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 3057, "BE_{K}:LPInvMass [120#circ<LPOpAngle<140#circ]; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 3058, "BE_{K}:LPInvMass [140#circ<LPOpAngle<160#circ]; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 3059, "BE_{K}:LPInvMass [160#circ<LPOpAngle<180#circ]; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);

  HB2( 4001, "BE_{K}:angleP1P2; B_{K} [GeV]; angle p1-p2", 200,-0.5,0.5,360,0.,180. );
    
  HB1( 600, "[Genfit] p-value",1000,0,1);
  HB1( 601, "[Genfit] x Pull",1000,-5,5);
  HB1( 602, "[Genfit] y Pull",1000,-5,5);
  HB1( 603, "[Genfit] z Pull",1000,-5,5);
  HB1( 604, "[Genfit] p_{x} Pull",1000,-5,5);
  HB1( 605, "[Genfit] p_{y} Pull",1000,-5,5);
  HB1( 606, "[Genfit] p_{z} Pull",1000,-5,5);
  HB1( 611, "[Genfit] x Pull pval>0.01",1000,-5,5);
  HB1( 612, "[Genfit] y Pull pval>0.01",1000,-5,5);
  HB1( 613, "[Genfit] z Pull pval>0.01",1000,-5,5);
  HB1( 614, "[Genfit] p_{x} Pull pval>0.01",1000,-5,5);
  HB1( 615, "[Genfit] p_{y} Pull pval>0.01",1000,-5,5);
  HB1( 616, "[Genfit] p_{z} Pull pval>0.01",1000,-5,5);

  HB1( 800, "[Genfit] p-value (LP)",1000,0,1);
  HB1( 801, "[Genfit] x Pull (LP)",1000,-5,5);
  HB1( 802, "[Genfit] y Pull (LP)",1000,-5,5);
  HB1( 803, "[Genfit] z Pull (LP)",1000,-5,5);
  HB1( 804, "[Genfit] p_{x} Pull (LP)",1000,-5,5);
  HB1( 805, "[Genfit] p_{y} Pull (LP)",1000,-5,5);
  HB1( 806, "[Genfit] p_{z} Pull (LP)",1000,-5,5);
  HB1( 811, "[Genfit] x Pull pval>0.01 (LP)",1000,-5,5);
  HB1( 812, "[Genfit] y Pull pval>0.01 (LP)",1000,-5,5);
  HB1( 813, "[Genfit] z Pull pval>0.01 (LP)",1000,-5,5);
  HB1( 814, "[Genfit] p_{x} Pull pval>0.01 (LP)",1000,-5,5);
  HB1( 815, "[Genfit] p_{y} Pull pval>0.01 (LP)",1000,-5,5);
  HB1( 816, "[Genfit] p_{z} Pull pval>0.01 (LP)",1000,-5,5);
  
  HB1( 10103, "[Genfit] p Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 10104, "[Genfit] #pi Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  HB1( 10105, "[Genfit] p2 Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0); 
  HB1( 10113, "[Genfit] p Momentum [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 10114, "[Genfit] #pi Momentum [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  HB1( 10115, "[Genfit] p2 Momentum [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  
  HB1( 10121, "[Genfit] p m2   [#Lambda decay];  M2 [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 630, 0., 1.26);
  HB1( 10122, "[Genfit] #pi m2 [#Lambda decay];  M2 [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 630, 0., 1.26);
  HB1( 10123, "[Genfit] p2 m2  [#Lambda decay];  M2 [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 630, 0., 1.26);
  HB1( 10131, "[Genfit] p m2   [#Lambdap decay]; M2 [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 630, 0., 1.26);
  HB1( 10132, "[Genfit] #pi m2 [#Lambdap decay]; M2 [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 630, 0., 1.26);
  HB1( 10133, "[Genfit] p2 m2  [#Lambdap decay]; M2 [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 630, 0., 1.26);
  
  HB1( 11101, "[Genfit] #Lambda p opening angle; #Lambda p opening angle [degree]; Counts ", 360., 0., 180.);
  HB1( 11102, "[Genfit] #Lambda p opening angle [-BEk<-0.1(GeV)]; #Lambda p opening angle [degree]; Counts ", 360., 0., 180.);
  HB1( 11103, "[Genfit] #Lambda p opening angle [-0.1<-BEk<0(GeV)]; #Lambda p opening angle [degree]; Counts ", 360., 0., 180.);
  HB1( 11104, "[Genfit] #Lambda p opening angle [0<-BEk<0.1(GeV)]; #Lambda p opening angle [degree]; Counts ", 360., 0., 180.);
  HB1( 11105, "[Genfit] #Lambda p opening angle [0.1<-BEk<0.2(GeV)]; #Lambda p opening angle [degree]; Counts ", 360., 0., 180.);
  HB1( 11106, "[Genfit] #Lambda p opening angle [0.2<-BEk<0.3(GeV)]; #Lambda p opening angle [degree]; Counts ", 360., 0., 180.);

  HB1( 11201, "[Genfit] Binding Energy of Kaon ; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11202, "[Genfit] Binding Energy of Kaon [thetaTPC<5#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11203, "[Genfit] Binding Energy of Kaon [thetaTPC<10#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11204, "[Genfit] Binding Energy of Kaon [w/ #Lambda]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11205, "[Genfit] Binding Energy of Kaon [thetaTPC<5#circ && w/ #Lambda]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11206, "[Genfit] Binding Energy of Kaon [thetaTPC<10#circ && w/ #Lambda]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11207, "[Genfit] Binding Energy of Kaon [w/ #Lambdap]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);  
  HB1( 11208, "[Genfit] Binding Energy of Kaon [thetaTPC<5#circ && w/ #Lambdap]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11209, "[Genfit] Binding Energy of Kaon [thetaTPC<10#circ && w/ #Lambdap]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11211, "[Genfit] #Lambda Invariant Mass ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 11212, "[Genfit] #Lambda Invariant Mass [-BEk<-0.1(GeV)]   ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 11213, "[Genfit] #Lambda Invariant Mass [-0.1<-BEk<0(GeV)] ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 11214, "[Genfit] #Lambda Invariant Mass [0<-BEk<0.1(GeV)]  ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 11215, "[Genfit] #Lambda Invariant Mass [0.1<-BEk<0.2(GeV)]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 11216, "[Genfit] #Lambda Invariant Mass [0.2<-BEk<0.3(GeV)]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 11217, "[Genfit] #LambdaP Invariant Mass ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 11218, "[Genfit] #LambdaP Invariant Mass [-BEk<-0.1(GeV)]   ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 11219, "[Genfit] #LambdaP Invariant Mass [-0.1<-BEk<0(GeV)] ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 11220, "[Genfit] #LambdaP Invariant Mass [0<-BEk<0.1(GeV)]  ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 11221, "[Genfit] #LambdaP Invariant Mass [0.1<-BEk<0.2(GeV)]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 11222, "[Genfit] #LambdaP Invariant Mass [0.2<-BEk<0.3(GeV)]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 11231, "[Genfit] #Lambda Invariant Mass [one more p]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 11232, "[Genfit] #Lambda Invariant Mass [one more p][-BEk<-0.1(GeV)]   ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 11233, "[Genfit] #Lambda Invariant Mass [one more p][-0.1<-BEk<0(GeV)] ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 11234, "[Genfit] #Lambda Invariant Mass [one more p][0<-BEk<0.1(GeV)]  ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 11235, "[Genfit] #Lambda Invariant Mass [one more p][0.1<-BEk<0.2(GeV)]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 11236, "[Genfit] #Lambda Invariant Mass [one more p][0.2<-BEk<0.3(GeV)]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);

  HB1( 11250, "[Genfit] Binding Energy of Kaon ; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11251, "[Genfit] Binding Energy of Kaon [0#circ<LPOpAngle<20#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11252, "[Genfit] Binding Energy of Kaon [20#circ<LPOpAngle<40#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11253, "[Genfit] Binding Energy of Kaon [40#circ<LPOpAngle<60#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11254, "[Genfit] Binding Energy of Kaon [60#circ<LPOpAngle<80#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11255, "[Genfit] Binding Energy of Kaon [80#circ<LPOpAngle<100#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11256, "[Genfit] Binding Energy of Kaon [100#circ<LPOpAngle<120#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11257, "[Genfit] Binding Energy of Kaon [120#circ<LPOpAngle<140#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11258, "[Genfit] Binding Energy of Kaon [140#circ<LPOpAngle<160#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11259, "[Genfit] Binding Energy of Kaon [160#circ<LPOpAngle<180#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11260, "[Genfit] mom transfer ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11261, "[Genfit] mom transfer [0#circ<LPOpAngle<20#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11262, "[Genfit] mom transfer [20#circ<LPOpAngle<40#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11263, "[Genfit] mom transfer [40#circ<LPOpAngle<60#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11264, "[Genfit] mom transfer [60#circ<LPOpAngle<80#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11265, "[Genfit] mom transfer [80#circ<LPOpAngle<100#circ];  Mom [GeV/#font[12]{c}];  Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11266, "[Genfit] mom transfer [100#circ<LPOpAngle<120#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11267, "[Genfit] mom transfer [120#circ<LPOpAngle<140#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11268, "[Genfit] mom transfer [140#circ<LPOpAngle<160#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11269, "[Genfit] mom transfer [160#circ<LPOpAngle<180#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);  
  HB1( 11270, "[Genfit] Lp Mom; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11271, "[Genfit] Lp Mom [0#circ<LPOpAngle<20#circ]   ;  Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11272, "[Genfit] Lp Mom [20#circ<LPOpAngle<40#circ]  ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11273, "[Genfit] Lp Mom [40#circ<LPOpAngle<60#circ]  ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11274, "[Genfit] Lp Mom [60#circ<LPOpAngle<80#circ]  ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11275, "[Genfit] Lp Mom [80#circ<LPOpAngle<100#circ] ;  Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11276, "[Genfit] Lp Mom [100#circ<LPOpAngle<120#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11277, "[Genfit] Lp Mom [120#circ<LPOpAngle<140#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11278, "[Genfit] Lp Mom [140#circ<LPOpAngle<160#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11279, "[Genfit] Lp Mom [160#circ<LPOpAngle<180#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0, 2.0);
  HB1( 11280, "[Genfit] MissMom; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0); // missMom = (momTrans-momLp)*unitvec_momTrans
  HB1( 11281, "[Genfit] MissMom [0#circ<LPOpAngle<20#circ]   ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 11282, "[Genfit] MissMom [20#circ<LPOpAngle<40#circ]  ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 11283, "[Genfit] MissMom [40#circ<LPOpAngle<60#circ]  ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);    
  HB1( 11284, "[Genfit] MissMom [60#circ<LPOpAngle<80#circ]  ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 11285, "[Genfit] MissMom [80#circ<LPOpAngle<100#circ] ; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 11286, "[Genfit] MissMom [100#circ<LPOpAngle<120#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 11287, "[Genfit] MissMom [120#circ<LPOpAngle<140#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 11288, "[Genfit] MissMom [140#circ<LPOpAngle<160#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 11289, "[Genfit] MissMom [160#circ<LPOpAngle<180#circ]; Mom [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  
  HB1( 11240, "[Genfit] #LambdaP Invariant Mass ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 11241, "[Genfit] #LambdaP Invariant Mass [0#circ<LPOpAngle<20#circ]   ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 11242, "[Genfit] #LambdaP Invariant Mass [20#circ<LPOpAngle<40#circ]  ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 11243, "[Genfit] #LambdaP Invariant Mass [40#circ<LPOpAngle<60#circ]  ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 11244, "[Genfit] #LambdaP Invariant Mass [60#circ<LPOpAngle<80#circ]  ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 11245, "[Genfit] #LambdaP Invariant Mass [80#circ<LPOpAngle<100#circ] ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 11246, "[Genfit] #LambdaP Invariant Mass [100#circ<LPOpAngle<120#circ]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 11247, "[Genfit] #LambdaP Invariant Mass [120#circ<LPOpAngle<140#circ]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 11248, "[Genfit] #LambdaP Invariant Mass [140#circ<LPOpAngle<160#circ]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 11249, "[Genfit] #LambdaP Invariant Mass [160#circ<LPOpAngle<180#circ]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);    
  
  HB1( 11401, "[Genfit] p Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11402, "[Genfit] p Momentum [#Lambda decay] [-BEk<-0.1(GeV)]   ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11403, "[Genfit] p Momentum [#Lambda decay] [-0.1<-BEk<0(GeV)] ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11404, "[Genfit] p Momentum [#Lambda decay] [0<-BEk<0.1(GeV)]  ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11405, "[Genfit] p Momentum [#Lambda decay] [0.1<-BEk<0.2(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11406, "[Genfit] p Momentum [#Lambda decay] [0.2<-BEk<0.3(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);      
  HB1( 11411, "[Genfit] #pi Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11412, "[Genfit] #pi Momentum [#Lambda decay] [-BEk<-0.1(GeV)]   ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11413, "[Genfit] #pi Momentum [#Lambda decay] [-0.1<-BEk<0(GeV)] ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11414, "[Genfit] #pi Momentum [#Lambda decay] [0<-BEk<0.1(GeV)]  ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11415, "[Genfit] #pi Momentum [#Lambda decay] [0.1<-BEk<0.2(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11416, "[Genfit] #pi Momentum [#Lambda decay] [0.2<-BEk<0.3(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);      
  HB1( 11431, "[Genfit] #Lambda Momentum [#Lambda decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11432, "[Genfit] #Lambda Momentum [#Lambda decay] [-BEk<-0.1(GeV)]   ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11433, "[Genfit] #Lambda Momentum [#Lambda decay] [-0.1<-BEk<0(GeV)] ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11434, "[Genfit] #Lambda Momentum [#Lambda decay] [0<-BEk<0.1(GeV)]  ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11435, "[Genfit] #Lambda Momentum [#Lambda decay] [0.1<-BEk<0.2(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11436, "[Genfit] #Lambda Momentum [#Lambda decay] [0.2<-BEk<0.3(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);  
  HB1( 11441, "[Genfit] p Momentum [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11442, "[Genfit] p Momentum [#Lambdap decay] [-BEk<-0.1(GeV)]   ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11443, "[Genfit] p Momentum [#Lambdap decay] [-0.1<-BEk<0(GeV)] ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11444, "[Genfit] p Momentum [#Lambdap decay] [0<-BEk<0.1(GeV)]  ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11445, "[Genfit] p Momentum [#Lambdap decay] [0.1<-BEk<0.2(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11446, "[Genfit] p Momentum [#Lambdap decay] [0.2<-BEk<0.3(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);      
  HB1( 11451, "[Genfit] #pi Momentum [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11452, "[Genfit] #pi Momentum [#Lambdap decay] [-BEk<-0.1(GeV)]   ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11453, "[Genfit] #pi Momentum [#Lambdap decay] [-0.1<-BEk<0(GeV)] ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11454, "[Genfit] #pi Momentum [#Lambdap decay] [0<-BEk<0.1(GeV)]  ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11455, "[Genfit] #pi Momentum [#Lambdap decay] [0.1<-BEk<0.2(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11456, "[Genfit] #pi Momentum [#Lambdap decay] [0.2<-BEk<0.3(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);  
  HB1( 11461, "[Genfit] p2 Momentum [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11462, "[Genfit] p2 Momentum [#Lambdap decay] [-BEk<-0.1(GeV)]   ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11463, "[Genfit] p2 Momentum [#Lambdap decay] [-0.1<-BEk<0(GeV)] ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11464, "[Genfit] p2 Momentum [#Lambdap decay] [0<-BEk<0.1(GeV)]  ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11465, "[Genfit] p2 Momentum [#Lambdap decay] [0.1<-BEk<0.2(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11466, "[Genfit] p2 Momentum [#Lambdap decay] [0.2<-BEk<0.3(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);  
  HB1( 11471, "[Genfit] #Lambda Momentum [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11472, "[Genfit] #Lambda Momentum [#Lambdap decay] [-BEk<-0.1(GeV)]   ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11473, "[Genfit] #Lambda Momentum [#Lambdap decay] [-0.1<-BEk<0(GeV)] ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11474, "[Genfit] #Lambda Momentum [#Lambdap decay] [0<-BEk<0.1(GeV)]  ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11475, "[Genfit] #Lambda Momentum [#Lambdap decay] [0.1<-BEk<0.2(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11476, "[Genfit] #Lambda Momentum [#Lambdap decay] [0.2<-BEk<0.3(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 11481, "[Genfit] MissMom [#Lambdap decay]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);// missMom = (momTrans-momLp)*unitvec_momTrans
  HB1( 11482, "[Genfit] MissMom [#Lambdap decay] [-BEk<-0.1(GeV)]   ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 11483, "[Genfit] MissMom [#Lambdap decay] [-0.1<-BEk<0(GeV)] ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);  
  HB1( 11484, "[Genfit] MissMom [#Lambdap decay] [0<-BEk<0.1(GeV)]  ; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 11485, "[Genfit] MissMom [#Lambdap decay] [0.1<-BEk<0.2(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  HB1( 11486, "[Genfit] MissMom [#Lambdap decay] [0.2<-BEk<0.3(GeV)]; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 400, -4.0, 4.0);
  
  
  // 2D 
  // LPOpeningAngle:mom of p1,pi,p2,L
  HB2( 12101, "[Genfit] LPOpAngle:MomP1; LPOpAngle[degree]; Momentum p1  [GeV/#font[12]{c}]", 360, 0.0, 180., 50, 0., 1.0);
  HB2( 12102, "[Genfit] LPOpAngle:MomPi; LPOpAngle[degree]; Momentum #pi [GeV/#font[12]{c}]", 360, 0.0, 180., 50, 0., 1.0);
  HB2( 12103, "[Genfit] LPOpAngle:MomP2; LPOpAngle[degree]; Momentum p2  [GeV/#font[12]{c}]", 360, 0.0, 180., 50, 0., 1.0);
  HB2( 12104, "[Genfit] LPOpAngle:MomL ; LPOpAngle[degree]; Momentum #Lambda [GeV/#font[12]{c}]", 360, 0.0, 180., 50, 0., 1.0);
  // HB2( 12111, "[Genfit] LPOpAngle:BEkaonTPC; LPOpAngle[degree]; Binding Energy [GeV/#font[12]{c}^2]", 360, 0.0, 180., 200, -0.5, 0.5);
  // Mom p1,pi,p2,L : mom of p1,pi,p2,L
  HB2( 12201, "[Genfit] MomCorr [#Lambda tagged]; Mom p1[GeV/#font[12]{c}]; Mom p1[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12202, "[Genfit] MomCorr [#Lambda tagged]; Mom p1[GeV/#font[12]{c}]; Mom pi[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12203, "[Genfit] MomCorr [#Lambda tagged]; Mom p1[GeV/#font[12]{c}]; Mom p2[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12204, "[Genfit] MomCorr [#Lambda tagged]; Mom p1[GeV/#font[12]{c}]; Mom #Lambda[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12205, "[Genfit] MomCorr [#Lambda tagged]; Mom pi[GeV/#font[12]{c}]; Mom pi[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12206, "[Genfit] MomCorr [#Lambda tagged]; Mom pi[GeV/#font[12]{c}]; Mom p2[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12207, "[Genfit] MomCorr [#Lambda tagged]; Mom pi[GeV/#font[12]{c}]; Mom #Lambda[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12208, "[Genfit] MomCorr [#Lambda tagged]; Mom p2[GeV/#font[12]{c}]; Mom p2[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12209, "[Genfit] MomCorr [#Lambda tagged]; Mom p2[GeV/#font[12]{c}]; Mom #Lambda[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12210, "[Genfit] MomCorr [#Lambda tagged]; Mom #Lambda[MeV/#font[12]{c}]; Mom #Lambda[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12211, "[Genfit] MomCorr [#Lambdap tagged]; Mom p1[GeV/#font[12]{c}]; Mom p1[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12212, "[Genfit] MomCorr [#Lambdap tagged]; Mom p1[GeV/#font[12]{c}]; Mom pi[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12213, "[Genfit] MomCorr [#Lambdap tagged]; Mom p1[GeV/#font[12]{c}]; Mom p2[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12214, "[Genfit] MomCorr [#Lambdap tagged]; Mom p1[GeV/#font[12]{c}]; Mom #Lambda[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12215, "[Genfit] MomCorr [#Lambdap tagged]; Mom pi[GeV/#font[12]{c}]; Mom pi[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12216, "[Genfit] MomCorr [#Lambdap tagged]; Mom pi[GeV/#font[12]{c}]; Mom p2[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12217, "[Genfit] MomCorr [#Lambdap tagged]; Mom pi[GeV/#font[12]{c}]; Mom #Lambda[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12218, "[Genfit] MomCorr [#Lambdap tagged]; Mom p2[GeV/#font[12]{c}]; Mom p2[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12219, "[Genfit] MomCorr [#Lambdap tagged]; Mom p2[GeV/#font[12]{c}]; Mom #Lambda[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);
  HB2( 12220, "[Genfit] MomCorr [#Lambdap tagged]; Mom #Lambda[MeV/#font[12]{c}]; Mom #Lambda[GeV/#font[12]{c}]", 50, 0., 1.0, 50, 0., 1.0);  
  HB2( 12301, "[Genfit] MomCorr; Mom #Lambdap[GeV/#font[12]{c}]; Mom Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0); 

  HB2( 13000, "[Genfit] LPOpAngle:MomLp ; LPOpAngle[degree]; Momentum #Lambdap [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);
  HB2( 13001, "[Genfit] LPOpAngle:MomLp [-BEk<-0.1(GeV)]   ; LPOpAngle[degree]; Momentum #Lambdap [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);
  HB2( 13002, "[Genfit] LPOpAngle:MomLp [-0.1<-BEk<0(GeV)] ; LPOpAngle[degree]; Momentum #Lambdap [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);
  HB2( 13003, "[Genfit] LPOpAngle:MomLp [0<-BEk<0.1(GeV)]  ; LPOpAngle[degree]; Momentum #Lambdap [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);
  HB2( 13004, "[Genfit] LPOpAngle:MomLp [0.1<-BEk<0.2(GeV)]; LPOpAngle[degree]; Momentum #Lambdap [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);
  HB2( 13005, "[Genfit] LPOpAngle:MomLp [0.2<-BEk<0.3(GeV)]; LPOpAngle[degree]; Momentum #Lambdap [GeV/#font[12]{c}]", 360, 0.0, 180., 100, 0., 2.0);  
  HB2( 13010, "[Genfit] LPOpAngle:LPInvMass ; LPOpAngle[degree]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]",360, 0.0, 180., 1750, 2.0, 5.5);
  HB2( 13011, "[Genfit] LPOpAngle:LPInvMass [-BEk<-0.1(GeV)]   ; LPOpAngle[degree]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]",360, 0.0, 180., 1750, 2.0, 5.5);
  HB2( 13012, "[Genfit] LPOpAngle:LPInvMass [-0.1<-BEk<0(GeV)] ; LPOpAngle[degree]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]",360, 0.0, 180., 1750, 2.0, 5.5);
  HB2( 13013, "[Genfit] LPOpAngle:LPInvMass [0<-BEk<0.1(GeV)]  ; LPOpAngle[degree]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]",360, 0.0, 180., 1750, 2.0, 5.5);
  HB2( 13014, "[Genfit] LPOpAngle:LPInvMass [0.1<-BEk<0.2(GeV)]; LPOpAngle[degree]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]",360, 0.0, 180., 1750, 2.0, 5.5);
  HB2( 13015, "[Genfit] LPOpAngle:LPInvMass [0.2<-BEk<0.3(GeV)]; LPOpAngle[degree]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]",360, 0.0, 180., 1750, 2.0, 5.5); 
  HB2( 13020, "[Genfit] MomTransfer:MomLp ; Momentum #Lambdap [GeV/#font[12]{c}] ;Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 100, 0., 2.0);
  HB2( 13021, "[Genfit] MomTransfer:MomLp [-BEk<-0.1(GeV)]  ; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 13022, "[Genfit] MomTransfer:MomLp [-0.1<-BEk<0(GeV)]; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);  
  HB2( 13023, "[Genfit] MomTransfer:MomLp [0<-BEk<0.1(GeV)] ; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 13024, "[Genfit] MomTransfer:MomLp [0.1<-BEk<0.2(GeV)]; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 13025, "[Genfit] MomTransfer:MomLp [0.2<-BEk<0.3(GeV)]; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 13031, "[Genfit] MomTransfer:MomLp [0#circ<LPOpAngle<20#circ]   ; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 13032, "[Genfit] MomTransfer:MomLp [20#circ<LPOpAngle<40#circ]  ; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 13033, "[Genfit] MomTransfer:MomLp [40#circ<LPOpAngle<60#circ]  ; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 13034, "[Genfit] MomTransfer:MomLp [60#circ<LPOpAngle<80#circ]  ; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 13035, "[Genfit] MomTransfer:MomLp [80#circ<LPOpAngle<100#circ] ; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 13036, "[Genfit] MomTransfer:MomLp [100#circ<LPOpAngle<120#circ]; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 13037, "[Genfit] MomTransfer:MomLp [120#circ<LPOpAngle<140#circ]; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 13038, "[Genfit] MomTransfer:MomLp [140#circ<LPOpAngle<160#circ]; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);  
  HB2( 13039, "[Genfit] MomTransfer:MomLp [160#circ<LPOpAngle<180#circ]; Momentum #Lambdap [GeV/#font[12]{c}]; Momentum Transfer [GeV/#font[12]{c}]", 100, 0., 2.0, 150, 0., 3.0);
  HB2( 13040, "[Genfit] BE_{K}:MomTransfer; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 13041, "[Genfit] BE_{K}:MomTransfer [0#circ<LPOpAngle<20#circ]   ; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 13042, "[Genfit] BE_{K}:MomTransfer [20#circ<LPOpAngle<40#circ]  ; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 13043, "[Genfit] BE_{K}:MomTransfer [40#circ<LPOpAngle<60#circ]  ; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 13044, "[Genfit] BE_{K}:MomTransfer [60#circ<LPOpAngle<80#circ]  ; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 13045, "[Genfit] BE_{K}:MomTransfer [80#circ<LPOpAngle<100#circ] ; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 13046, "[Genfit] BE_{K}:MomTransfer [100#circ<LPOpAngle<120#circ]; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 13047, "[Genfit] BE_{K}:MomTransfer [120#circ<LPOpAngle<140#circ]; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 13048, "[Genfit] BE_{K}:MomTransfer [140#circ<LPOpAngle<160#circ]; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 13049, "[Genfit] BE_{K}:MomTransfer [160#circ<LPOpAngle<180#circ]; B_{K} [GeV]; Momentum Transfer [GeV/#font[12]{c}]", 200, -0.5, 0.5, 100, 0., 2.0);
  HB2( 13050, "[Genfit] BE_{K}:LPInvMass; B_{K} [GeV]; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 13051, "[Genfit] BE_{K}:LPInvMass [0#circ<LPOpAngle<20#circ]   ; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 13052, "[Genfit] BE_{K}:LPInvMass [20#circ<LPOpAngle<40#circ]  ; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 13053, "[Genfit] BE_{K}:LPInvMass [40#circ<LPOpAngle<60#circ]  ; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 13054, "[Genfit] BE_{K}:LPInvMass [60#circ<LPOpAngle<80#circ]  ; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 13055, "[Genfit] BE_{K}:LPInvMass [80#circ<LPOpAngle<100#circ] ; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 13056, "[Genfit] BE_{K}:LPInvMass [100#circ<LPOpAngle<120#circ]; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 13057, "[Genfit] BE_{K}:LPInvMass [120#circ<LPOpAngle<140#circ]; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 13058, "[Genfit] BE_{K}:LPInvMass [140#circ<LPOpAngle<160#circ]; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);
  HB2( 13059, "[Genfit] BE_{K}:LPInvMass [160#circ<LPOpAngle<180#circ]; B_{K} [GeV]; Invariant Mass #Lambdap [/0.002 GeV/#font[12]{c}^{2}]", 200, -0.5, 0.5,  1750, 2.0, 5.5);

  // scat K
  HB1( 20011, "TPC tracks; track number; Counts", 10, 0, 10);
  HB2( 20021, "PID vs TPC tracks; track number; Counts", 20, -10, 10, 10, 0, 10);      
  HB2( 20100, "<dE/dx> [scat pK] nsigma; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 20101, "<dE/dx> [scat pK] nsigma; -p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 20102, "<dE/dx> [scat pK] nsigma; +p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 20110, "<dE/dx> [scat K] nsigma; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 20111, "<dE/dx> [scat K] nsigma; -p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 20112, "<dE/dx> [scat K] nsigma; +p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB1( 21001, "[Genfit] Binding Energy of Kaon [scat K]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);  
  HB1( 21002, "[Genfit] Binding Energy of Kaon [scat K][thetaTPC<5#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21003, "[Genfit] Binding Energy of Kaon [scat K][thetaTPC<10#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21004, "[Genfit] Binding Energy of Kaon [scat K][3.5#circ<thetaTPC<4.5#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21010, "[Genfit] Binding Energy of Kaon [scat K][t=0]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21011, "[Genfit] Binding Energy of Kaon [scat K][t=1]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);  
  HB1( 21012, "[Genfit] Binding Energy of Kaon [scat K][t=2]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21013, "[Genfit] Binding Energy of Kaon [scat K][t=3]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21020, "[Genfit] Binding Energy of Kaon [scat K][3.5#circ<thetaTPC<4.5#circ][t=0]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21021, "[Genfit] Binding Energy of Kaon [scat K][3.5#circ<thetaTPC<4.5#circ][t=1]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);  
  HB1( 21022, "[Genfit] Binding Energy of Kaon [scat K][3.5#circ<thetaTPC<4.5#circ][t=2]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21023, "[Genfit] Binding Energy of Kaon [scat K][3.5#circ<thetaTPC<4.5#circ][t=3]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21030, "[Genfit] Binding Energy of Kaon [scat K][3.5#circ<thetaTPC<4.5#circ][m2>0.05]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21031, "[Genfit] Binding Energy of Kaon [scat K][3.5#circ<thetaTPC<4.5#circ][m2>0.10]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5); 
  HB1( 21032, "[Genfit] Binding Energy of Kaon [scat K][3.5#circ<thetaTPC<4.5#circ][m2>0.15]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21033, "[Genfit] Binding Energy of Kaon [scat K][3.5#circ<thetaTPC<4.5#circ][m2>0.20]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21034, "[Genfit] Binding Energy of Kaon [scat K][3.5#circ<thetaTPC<4.5#circ][m2>0.25]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5); 
  HB1( 21040, "[Genfit] Binding Energy of Kaon [scat K][thetaTPC<5.0#circ][m2>0.05]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21041, "[Genfit] Binding Energy of Kaon [scat K][thetaTPC<5.0#circ][m2>0.10]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);  
  HB1( 21042, "[Genfit] Binding Energy of Kaon [scat K][thetaTPC<5.0#circ][m2>0.15]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21043, "[Genfit] Binding Energy of Kaon [scat K][thetaTPC<5.0#circ][m2>0.20]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21044, "[Genfit] Binding Energy of Kaon [scat K][thetaTPC<5.0#circ][m2>0.25]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21050, "[Genfit] Binding Energy of Kaon [scat K][thetaTPC<10.0#circ][m2>0.05]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21051, "[Genfit] Binding Energy of Kaon [scat K][thetaTPC<10.0#circ][m2>0.10]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);  
  HB1( 21052, "[Genfit] Binding Energy of Kaon [scat K][thetaTPC<10.0#circ][m2>0.15]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21053, "[Genfit] Binding Energy of Kaon [scat K][thetaTPC<10.0#circ][m2>0.20]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 21054, "[Genfit] Binding Energy of Kaon [scat K][thetaTPC<10.0#circ][m2>0.25]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);  
 
  HB1( 21101, "[Genfit] K Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 21201, "[Genfit] K Mass2; mass2 [GeV]; Counts [/0.02 GeV]", 800, -2.0, 2.0);
  HB2( 23000, "[Genfit] Mass2:dEdx [nega][-B_{K}<0.3(GeV), 3.5#circ<#theta_{Kp}<4.5#circ]; mass2 [GeV]; <dEdx> [arb.]", 1000,-5.0, 5.0, nbindedx, mindedx, maxdedx);
  HB2( 23001, "[Genfit] Mass2:dEdx [nega][0.0<mom<0.1(GeV/#font[12]{c}), -B_{K}<0.3(GeV), 3.5#circ<#theta_{Kp}<4.5#circ]; mass2 [GeV]; <dEdx> [arb.]", 1000,-5.0, 5.0, nbindedx, mindedx, maxdedx);
  HB2( 23002, "[Genfit] Mass2:dEdx [nega][0.1<mom<0.2(GeV/#font[12]{c}), -B_{K}<0.3(GeV), 3.5#circ<#theta_{Kp}<4.5#circ]; mass2 [GeV]; <dEdx> [arb.]", 1000,-5.0, 5.0, nbindedx, mindedx, maxdedx);
  HB2( 23003, "[Genfit] Mass2:dEdx [nega][0.2<mom<0.3(GeV/#font[12]{c}), -B_{K}<0.3(GeV), 3.5#circ<#theta_{Kp}<4.5#circ]; mass2 [GeV]; <dEdx> [arb.]", 1000,-5.0, 5.0, nbindedx, mindedx, maxdedx);
  HB2( 23004, "[Genfit] Mass2:dEdx [nega][0.3<mom<0.4(GeV/#font[12]{c}), -B_{K}<0.3(GeV), 3.5#circ<#theta_{Kp}<4.5#circ]; mass2 [GeV]; <dEdx> [arb.]", 1000,-5.0, 5.0, nbindedx, mindedx, maxdedx);
  HB2( 23005, "[Genfit] Mass2:dEdx [nega][0.4<mom<0.5(GeV/#font[12]{c}), -B_{K}<0.3(GeV), 3.5#circ<#theta_{Kp}<4.5#circ]; mass2 [GeV]; <dEdx> [arb.]", 1000,-5.0, 5.0, nbindedx, mindedx, maxdedx);
  HB2( 23010, "[Genfit] Mass2:dEdx [posi][-B_{K}<0.3(GeV), 3.5#circ<#theta_{Kp}<4.5#circ]; mass2 [GeV]; <dEdx> [arb.]", 1000,-5.0, 5.0, nbindedx, mindedx, maxdedx);
  HB2( 23011, "[Genfit] Mass2:dEdx [posi][0.0<mom<0.1(GeV/#font[12]{c}), -B_{K}<0.3(GeV), 3.5#circ<#theta_{Kp}<4.5#circ]; mass2 [GeV]; <dEdx> [arb.]", 1000,-5.0, 5.0, nbindedx, mindedx, maxdedx);
  HB2( 23012, "[Genfit] Mass2:dEdx [posi][0.1<mom<0.2(GeV/#font[12]{c}), -B_{K}<0.3(GeV), 3.5#circ<#theta_{Kp}<4.5#circ]; mass2 [GeV]; <dEdx> [arb.]", 1000,-5.0, 5.0, nbindedx, mindedx, maxdedx);
  HB2( 23013, "[Genfit] Mass2:dEdx [posi][0.2<mom<0.3(GeV/#font[12]{c}), -B_{K}<0.3(GeV), 3.5#circ<#theta_{Kp}<4.5#circ]; mass2 [GeV]; <dEdx> [arb.]", 1000,-5.0, 5.0, nbindedx, mindedx, maxdedx);
  HB2( 23014, "[Genfit] Mass2:dEdx [posi][0.3<mom<0.4(GeV/#font[12]{c}), -B_{K}<0.3(GeV), 3.5#circ<#theta_{Kp}<4.5#circ]; mass2 [GeV]; <dEdx> [arb.]", 1000,-5.0, 5.0, nbindedx, mindedx, maxdedx);
  HB2( 23015, "[Genfit] Mass2:dEdx [posi][0.4<mom<0.5(GeV/#font[12]{c}), -B_{K}<0.3(GeV), 3.5#circ<#theta_{Kp}<4.5#circ]; mass2 [GeV]; <dEdx> [arb.]", 1000,-5.0, 5.0, nbindedx, mindedx, maxdedx);

  // Kbar0 search
  HB1( 30011, "TPC tracks [Kbar0]; track number; Counts", 10, 0, 10); 
  HB2( 30100, "<dE/dx> [Kbar0] nsigma; p/q  [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 30101, "<dE/dx> [Kbar0] nsigma; -p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 30102, "<dE/dx> [Kbar0] nsigma; +p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);  
  HB1( 31001, "[Genfit] Binding Energy of Kaon [Kbar0]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);  
  HB1( 31002, "[Genfit] Binding Energy of Kaon [Kbar0][thetaTPC<5#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 31003, "[Genfit] Binding Energy of Kaon [Kbar0][thetaTPC<10#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 31004, "[Genfit] Binding Energy of Kaon [Kbar0][3.5#circ<thetaTPC<4.5#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);    
  HB1( 31101, "[Genfit] Kbar0 Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 31301, "[Genfit] Kbar0 InvariantMass; invariant mass [GeV]; Counts [/0.002 GeV]", 500, 0., 1.0);
  HB1( 31401, "[Genfit] Kbar0 InvariantMass LambdaLike; invariant mass [GeV]; Counts [/0.002 GeV]", 1000, 0., 2.0);
  HB1( 32101, "[Genfit] Missing Momentum of K-p->Kb0X; momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV]", 400, 0., 2.0);
  HB1( 32301, "[Genfit] Missing Mass of K-p->Kb0X; mass [GeV]; Counts [/0.002 GeV]", 750, 0., 1.5);
  
  HB1( 33101, "[Genfit] pi- Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 33201, "[Genfit] pi- Mass2; mass2 [GeV]; Counts [/0.02 GeV]", 400, 0., 2.0);  
  HB1( 34101, "[Genfit] pi+ Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);  
  HB1( 34201, "[Genfit] pi+ Mass2; mass2 [GeV]; Counts [/0.02 GeV]", 400, 0., 2.0);  

  HBTree( "tpc", "tree of DstTPCTracking" );
  tree->Branch( "status", &event.status );
  tree->Branch( "runnum", &event.runnum );
  tree->Branch( "evnum", &event.evnum );
  tree->Branch( "trigpat", &event.trigpat );
  tree->Branch( "trigflag", &event.trigflag );

  tree->Branch( "nhHtof", &event.nhHtof );
  tree->Branch( "HtofSeg", &event.HtofSeg );
  tree->Branch( "tHtof", &event.tHtof );
  tree->Branch( "dtHtof", &event.dtHtof );
  tree->Branch( "deHtof", &event.deHtof );
  tree->Branch( "posHtof", &event.posHtof );

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "trackid", &event.trackid );
  tree->Branch( "isBeam", &event.isBeam );
  tree->Branch( "isKurama", &event.isKurama );
  tree->Branch( "isK18", &event.isK18 );
  tree->Branch( "isAccidental", &event.isAccidental );
  tree->Branch( "charge", &event.charge );
  tree->Branch( "pid", &event.pid );
  /*** temp
  tree->Branch( "chisqr", &event.chisqr );
  tree->Branch( "pval", &event.pval );
  ***/
  tree->Branch( "helix_cx", &event.helix_cx );
  tree->Branch( "helix_cy", &event.helix_cy );
  tree->Branch( "helix_z0", &event.helix_z0 );
  tree->Branch( "helix_r", &event.helix_r );
  tree->Branch( "helix_dz", &event.helix_dz );

  tree->Branch( "dE", &event.dE );
  tree->Branch( "dEdx", &event.dEdx );
  tree->Branch( "mom0", &event.mom0 );
  tree->Branch( "path", &event.path );

  tree->Branch( "hitlayer", &event.hitlayer );
  tree->Branch( "hitpos_x", &event.hitpos_x );
  tree->Branch( "hitpos_y", &event.hitpos_y );
  tree->Branch( "hitpos_z", &event.hitpos_z );
  tree->Branch( "calpos_x", &event.calpos_x );
  tree->Branch( "calpos_y", &event.calpos_y );
  tree->Branch( "calpos_z", &event.calpos_z );
  /***  
  tree->Branch( "residual", &event.residual );
  tree->Branch( "residual_x", &event.residual_x );
  tree->Branch( "residual_y", &event.residual_y );
  tree->Branch( "residual_z", &event.residual_z );
  ***/
  tree->Branch( "resolution_x", &event.resolution_x);
  tree->Branch( "resolution_y", &event.resolution_y);
  tree->Branch( "resolution_z", &event.resolution_z);

  tree->Branch( "helix_t", &event.helix_t );
  /***  
  tree->Branch( "alpha", &event.alpha);
  tree->Branch( "pathhit", &event.pathhit);
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
  ***/

  tree->Branch( "ntK18", &event.ntK18);
  tree->Branch( "chisqrK18", &event.chisqrK18);
  tree->Branch( "pK18", &event.pK18);
  tree->Branch( "thetaK18", &event.thetaK18);  
  tree->Branch( "xtgtK18", &event.xtgtK18);
  tree->Branch( "ytgtK18", &event.ytgtK18);
  tree->Branch( "utgtK18", &event.utgtK18);
  tree->Branch( "vtgtK18", &event.vtgtK18);
  // tree->Branch( "xvpHS", &event.xvpHS);
  // tree->Branch( "yvpHS", &event.yvpHS);
  // tree->Branch( "zvpHS", &event.zvpHS);
  // tree->Branch( "xtgtHS", &event.xtgtHS);
  // tree->Branch( "ytgtHS", &event.ytgtHS);
  // tree->Branch( "ztgtHS", &event.ztgtHS);

#if SaveTPCK18
  tree->Branch( "isgoodTPCK18", &event.isgoodTPCK18);
  tree->Branch( "chisqrTPCK18", &event.chisqrTPCK18);
  tree->Branch( "qTPCK18", &event.qTPCK18);
  tree->Branch( "pTPCK18", &event.pTPCK18);
  tree->Branch( "thetaTPCK18", &event.thetaTPCK18);  
  tree->Branch( "xtgtTPCK18", &event.xtgtTPCK18);
  tree->Branch( "ytgtTPCK18", &event.ytgtTPCK18);
  tree->Branch( "utgtTPCK18", &event.utgtTPCK18);
  tree->Branch( "vtgtTPCK18", &event.vtgtTPCK18);
  tree->Branch( "lhtofTPCK18", &event.lhtofTPCK18);
  tree->Branch( "xhtofTPCK18", &event.xhtofTPCK18);
  tree->Branch( "yhtofTPCK18", &event.yhtofTPCK18);
  tree->Branch( "lvpTPCK18", &event.lvpTPCK18);
  tree->Branch( "xvpTPCK18", &event.xvpTPCK18);
  tree->Branch( "yvpTPCK18", &event.yvpTPCK18);
#endif  
  tree->Branch( "ntKurama",     &event.ntKurama);
  tree->Branch( "chisqrKurama", &event.chisqrKurama);
  tree->Branch( "pKurama", &event.pKurama);
  tree->Branch( "qKurama", &event.qKurama);
  tree->Branch( "m2Kurama", &event.m2);
  tree->Branch( "xtgtKurama", &event.xtgtKurama);
  tree->Branch( "ytgtKurama", &event.ytgtKurama);
  tree->Branch( "utgtKurama", &event.utgtKurama);
  tree->Branch( "vtgtKurama", &event.vtgtKurama);
  tree->Branch( "momTransfer", &event.momTransfer);
  tree->Branch( "thetaK18Kurama", &event.thetaK18Kurama);    
  tree->Branch( "isgoodTPCKurama", &event.isgoodTPCKurama);
  tree->Branch( "kflagTPCKurama", &event.kflagTPCKurama);
  tree->Branch( "pflagTPCKurama", &event.pflagTPCKurama);
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

  tree->Branch("nKK", &event.nKK);
  tree->Branch("Kflag", &event.Kflag);
  tree->Branch("Pflag", &event.Pflag);
  tree->Branch("MissMass", &event.MissMass);
  tree->Branch("MissMassCorr", &event.MissMassCorr);
  tree->Branch("MissMassCorrDE", &event.MissMassCorrDE);
  tree->Branch("vtx", &event.vtx);
  tree->Branch("vty", &event.vty);
  tree->Branch("vtz", &event.vtz);
  tree->Branch("closeDist", &event.closeDist);

  tree->Branch( "isgoodTPC", &event.isgoodTPC);
  tree->Branch( "insideTPC", &event.insideTPC);
  tree->Branch( "vtxTPC", &event.vtxTPC);
  tree->Branch( "vtyTPC", &event.vtyTPC);
  tree->Branch( "vtzTPC", &event.vtzTPC);
  tree->Branch( "closeDistTPC", &event.closeDistTPC);
  tree->Branch( "MissMassTPC", &event.MissMassTPC);
  tree->Branch( "MissMassCorrTPC", &event.MissMassCorrTPC);
  tree->Branch( "MissMassCorrDETPC", &event.MissMassCorrDETPC);
  tree->Branch( "MissMassNuclTPC", &event.MissMassNuclTPC);
  tree->Branch( "MissMassNuclCorrTPC", &event.MissMassNuclCorrTPC);
  tree->Branch( "MissMassNuclCorrDETPC", &event.MissMassNuclCorrDETPC);
  tree->Branch( "BEkaonTPC", &event.BEkaonTPC);
  tree->Branch( "pOrgTPC", &event.pOrgTPC);
  tree->Branch( "pCalcTPC", &event.pCalcTPC);
  tree->Branch( "pCorrTPC", &event.pCorrTPC);
  tree->Branch( "pCorrDETPC", &event.pCorrDETPC);
  tree->Branch( "pCalcTPC", &event.pCalcTPC);
  tree->Branch( "thetaCMTPC", &event.thetaCMTPC);
  tree->Branch( "costCMTPC", &event.costCMTPC);
  //  tree->Branch( "pCalcDETPC", &event.pCalcDETPC);
  //  tree->Branch( "thetaCMDETPC", &event.thetaCMDETPC);
  //  tree->Branch( "costCMDETPC", &event.costCMDETPC);
  //  tree->Branch( "xistarpCalcDETPC", &event.xistarpCalcDETPC);
  // tree->Branch( "xistarthetaCMDETPC", &event.xistarthetaCMDETPC);
  // tree->Branch( "xistarcostCMDETPC", &event.xistarcostCMDETPC);
  // tree->Branch( "kpscatpCalcTPC", &event.kpscatpCalcTPC);
  // tree->Branch( "kpscatthetaCMTPC", &event.kpscatthetaCMTPC);
  // tree->Branch( "kpscatcostCMTPC", &event.kpscatcostCMTPC);
  // tree->Branch( "kpscatpCalcDETPC", &event.kpscatpCalcDETPC);
  // tree->Branch( "kpscatthetaCMDETPC", &event.kpscatthetaCMDETPC);
  // tree->Branch( "kpscatcostCMDETPC", &event.kpscatcostCMDETPC);
  tree->Branch( "thetaTPC", &event.thetaTPC);
  tree->Branch( "xbTPC", &event.xbTPC);
  tree->Branch( "ybTPC", &event.ybTPC);
  tree->Branch( "ubTPC", &event.ubTPC);
  tree->Branch( "vbTPC", &event.vbTPC);
  tree->Branch( "xsTPC", &event.xsTPC);
  tree->Branch( "ysTPC", &event.ysTPC);
  tree->Branch( "usTPC", &event.usTPC);
  tree->Branch( "vsTPC", &event.vsTPC);

  tree->Branch("Lflag", &event.lflag);
  tree->Branch("ForwardLflag", &event.forwardlflag);
  // tree->Branch("Xiflag", &event.xiflag);
  // tree->Branch("XiMass", &event.ximass);
  // tree->Branch("XiDecayVtx_x", &event.xidecayvtx_x);
  // tree->Branch("XiDecayVtx_y", &event.xidecayvtx_y);
  // tree->Branch("XiDecayVtx_z", &event.xidecayvtx_z);
  // tree->Branch("XiMom_x", &event.ximom_x);
  // tree->Branch("XiMom_y", &event.ximom_y);
  // tree->Branch("XiMom_z", &event.ximom_z);
  // tree->Branch("XiVtxCloseDist", &event.lpi_dist);
  tree->Branch("LambdaMass", &event.lmass);
  tree->Branch("LambdaDecayVtx_x", &event.ldecayvtx_x);
  tree->Branch("LambdaDecayVtx_y", &event.ldecayvtx_y);
  tree->Branch("LambdaDecayVtx_z", &event.ldecayvtx_z);
  tree->Branch("LambdaMom_x", &event.lmom_x);
  tree->Branch("LambdaMom_y", &event.lmom_y);
  tree->Branch("LambdaMom_z", &event.lmom_z);
  tree->Branch("LambdaVtxCloseDist", &event.ppi_dist);
  tree->Branch("LDecaysTrackId", &event.ldecays_id);
  tree->Branch("LDecaysMom", &event.ldecays_mom);
  tree->Branch("LDecaysMom_x", &event.ldecays_mom_x);
  tree->Branch("LDecaysMom_y", &event.ldecays_mom_y);
  tree->Branch("LDecaysMom_z", &event.ldecays_mom_z);
  tree->Branch("LDecaysMomRes", &event.ldecays_res_mom);
  tree->Branch("LDecaysMomRes_x", &event.ldecays_res_mom_x);
  tree->Branch("LDecaysMomRes_y", &event.ldecays_res_mom_y);
  tree->Branch("LDecaysMomRes_z", &event.ldecays_res_mom_z);
  tree->Branch("LDecaysMomRes_t", &event.ldecays_res_mom_t);
  tree->Branch("LDecaysThRes", &event.ldecays_res_th);
  tree->Branch("LDecaysPhRes", &event.ldecays_res_ph);
  tree->Branch("LDecaysThCov", &event.ldecays_cov_mom_th);
  tree->Branch("LDecaysPhCov", &event.ldecays_cov_mom_ph);
  tree->Branch("LDecaysMomCov_xy", &event.ldecays_cov_mom_xy);
  tree->Branch("LDecaysMomCov_yz", &event.ldecays_cov_mom_yz);
  tree->Branch("LDecaysMomCov_zx", &event.ldecays_cov_mom_zx);
  tree->Branch("angleLambda", &event.angleLambda);

  tree->Branch("LPflag", &event.lp2flag);
  tree->Branch("LambdaMassLP", &event.lmass_lp2);
  tree->Branch("LPDecayVtx_x", &event.lp2decayvtx_x);
  tree->Branch("LPDecayVtx_y", &event.lp2decayvtx_y);
  tree->Branch("LPDecayVtx_z", &event.lp2decayvtx_z);
  tree->Branch("LPMom_x", &event.lp2mom_x);
  tree->Branch("LPMom_y", &event.lp2mom_y);
  tree->Branch("LPMom_z", &event.lp2mom_z);
  tree->Branch("LPVtxCloseDist", &event.lp2_dist);
  tree->Branch("LPDecaysTrackId", &event.lp2decays_id);
  tree->Branch("LPDecaysMom", &event.lp2decays_mom);
  tree->Branch("LPDecaysMom_x", &event.lp2decays_mom_x);
  tree->Branch("LPDecaysMom_y", &event.lp2decays_mom_y);
  tree->Branch("LPDecaysMom_z", &event.lp2decays_mom_z);
  tree->Branch("LPDecaysMomRes", &event.lp2decays_res_mom);
  tree->Branch("LPDecaysMomRes_x", &event.lp2decays_res_mom_x);
  tree->Branch("LPDecaysMomRes_y", &event.lp2decays_res_mom_y);
  tree->Branch("LPDecaysMomRes_z", &event.lp2decays_res_mom_z);
  tree->Branch("LPDecaysMomRes_t", &event.lp2decays_res_mom_t);
  tree->Branch("LPDecaysThRes", &event.lp2decays_res_th);
  tree->Branch("LPDecaysPhRes", &event.lp2decays_res_ph);
  tree->Branch("LPDecaysThCov", &event.lp2decays_cov_mom_th);
  tree->Branch("LPDecaysPhCov", &event.lp2decays_cov_mom_ph);
  tree->Branch("LPDecaysMomCov_xy", &event.lp2decays_cov_mom_xy);
  tree->Branch("LPDecaysMomCov_yz", &event.lp2decays_cov_mom_yz);
  tree->Branch("LPDecaysMomCov_zx", &event.lp2decays_cov_mom_zx);
  tree->Branch("LPOpeningAngle", &event.lp2angle);
  tree->Branch("LPMass", &event.lp2mass);

  // tree->Branch("GFXiflag", &event.GFxiflag);
  // tree->Branch("GFXiMass", &event.GFximass);
  // tree->Branch("GFXiDecayVtx_x", &event.GFxidecayvtx_x);
  // tree->Branch("GFXiDecayVtx_y", &event.GFxidecayvtx_y);
  // tree->Branch("GFXiDecayVtx_z", &event.GFxidecayvtx_z);
  // tree->Branch("GFXiMom_x", &event.GFximom_x);
  // tree->Branch("GFXiMom_y", &event.GFximom_y);
  // tree->Branch("GFXiMom_z", &event.GFximom_z);
  // tree->Branch("GFXiVtxCloseDist", &event.GFlpi_dist);
  tree->Branch("GFLambdaMass", &event.GFlmass);
  tree->Branch("GFLambdaDecayVtx_x", &event.GFldecayvtx_x);
  tree->Branch("GFLambdaDecayVtx_y", &event.GFldecayvtx_y);
  tree->Branch("GFLambdaDecayVtx_z", &event.GFldecayvtx_z);
  tree->Branch("GFLambdaMom_x", &event.GFlmom_x);
  tree->Branch("GFLambdaMom_y", &event.GFlmom_y);
  tree->Branch("GFLambdaMom_z", &event.GFlmom_z);
  tree->Branch("GFLambdaVtxCloseDist", &event.GFppi_dist);
  tree->Branch("GFDecaysMassSquare", &event.GFdecays_m2);
  tree->Branch("GFDecaysMom", &event.GFdecays_mom);
  tree->Branch("GFDecaysMom_x", &event.GFdecays_mom_x);
  tree->Branch("GFDecaysMom_y", &event.GFdecays_mom_y);
  tree->Branch("GFDecaysMom_z", &event.GFdecays_mom_z);
  tree->Branch("GFMomLoss", &event.GFmomloss);
  tree->Branch("GFELoss", &event.GFeloss);

  tree->Branch("GFLambdaPMass", &event.GFlp2mass);
  tree->Branch("GFLambdaPDecayVtx_x", &event.GFlp2decayvtx_x);
  tree->Branch("GFLambdaPDecayVtx_y", &event.GFlp2decayvtx_y);
  tree->Branch("GFLambdaPDecayVtx_z", &event.GFlp2decayvtx_z);
  tree->Branch("GFLambdaPMom_x", &event.GFlp2mom_x);
  tree->Branch("GFLambdaPMom_y", &event.GFlp2mom_y);
  tree->Branch("GFLambdaPMom_z", &event.GFlp2mom_z);
  tree->Branch("GFLambdaPVtxCloseDist", &event.GFlp2_dist);
  tree->Branch("GFDecaysMassSquareLP", &event.GFdecays_m2_lp2);
  tree->Branch("GFDecaysMomLP", &event.GFdecays_mom_lp2);
  tree->Branch("GFDecaysMom_x_LP", &event.GFdecays_mom_x_lp2);
  tree->Branch("GFDecaysMom_y_LP", &event.GFdecays_mom_y_lp2);
  tree->Branch("GFDecaysMom_z_LP", &event.GFdecays_mom_z_lp2);
  tree->Branch("GFMomLossLP", &event.GFmomloss_lp2);
  tree->Branch("GFELossLP", &event.GFeloss_lp2);

  //track fitting results
  tree->Branch("GFstatus", &event.GFstatus);
  tree->Branch("GFntTpc", &event.GFntTpc);
  tree->Branch("GFinside", &event.GFinside);  
  tree->Branch("GFcharge", &event.GFcharge);
  tree->Branch("GFchisqr", &event.GFchisqr);
  tree->Branch("GFtof", &event.GFtof);
  tree->Branch("GFtracklen", &event.GFtracklen);
  tree->Branch("GFpval", &event.GFpval);
  tree->Branch("GFchisqrPos", &event.GFchisqrPos);
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
  
  /*** temp
  tree->Branch("GFresidual_x", &event.GFresidual_x);
  tree->Branch("GFresidual_y", &event.GFresidual_y);
  tree->Branch("GFresidual_z", &event.GFresidual_z);
  tree->Branch("GFresidual_px", &event.GFresidual_px);
  tree->Branch("GFresidual_py", &event.GFresidual_py);
  tree->Branch("GFresidual_pz", &event.GFresidual_pz);
  tree->Branch("GFresidual_p", &event.GFresidual_p);
  tree->Branch("GFresidual6D_x", &event.GFresidual6D_x);
  tree->Branch("GFresidual6D_y", &event.GFresidual6D_y);
  tree->Branch("GFresidual6D_z", &event.GFresidual6D_z);
  tree->Branch("GFresidual6D_px", &event.GFresidual6D_px);
  tree->Branch("GFresidual6D_py", &event.GFresidual6D_py);
  tree->Branch("GFresidual6D_pz", &event.GFresidual6D_pz);
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
  ***/
  
  //track fitting results
  tree->Branch("GFstatusLP", &event.GFstatus_lp2);
  tree->Branch("GFntTpcLP", &event.GFntTpc_lp2);
  tree->Branch("GFchargeLP", &event.GFcharge_lp2);
  tree->Branch("GFchisqrLP", &event.GFchisqr_lp2);
  tree->Branch("GFtofLP", &event.GFtof_lp2);
  tree->Branch("GFtracklenLP", &event.GFtracklen_lp2);
  tree->Branch("GFpvalLP", &event.GFpval_lp2);
  tree->Branch("GFchisqrPosLP", &event.GFchisqrPos_lp2);
  tree->Branch("GFpvalPosLP", &event.GFpvalPos_lp2);
  tree->Branch("GFfitstatusLP", &event.GFfitstatus_lp2);
  tree->Branch("GFpdgcodeLP", &event.GFpdgcode_lp2);
  tree->Branch("GFnhtrackLP", &event.GFnhtrack_lp2);
  tree->Branch("GFlayerLP", &event.GFlayer_lp2);
  tree->Branch("GFrowLP", &event.GFrow_lp2);
  tree->Branch("GFposLP_x", &event.GFpos_x_lp2);
  tree->Branch("GFposLP_y", &event.GFpos_y_lp2);
  tree->Branch("GFposLP_z", &event.GFpos_z_lp2);
  tree->Branch("GFmomLP", &event.GFmom_lp2);
  tree->Branch("GFmomLP_x", &event.GFmom_x_lp2);
  tree->Branch("GFmomLP_y", &event.GFmom_y_lp2);
  tree->Branch("GFmomLP_z", &event.GFmom_z_lp2);
  /***
  tree->Branch("GFresidualLP_x", &event.GFresidual_x_lp2);
  tree->Branch("GFresidualLP_y", &event.GFresidual_y_lp2);
  tree->Branch("GFresidualLP_z", &event.GFresidual_z_lp2);
  tree->Branch("GFresidualLP_px", &event.GFresidual_px_lp2);
  tree->Branch("GFresidualLP_py", &event.GFresidual_py_lp2);
  tree->Branch("GFresidualLP_pz", &event.GFresidual_pz_lp2);
  tree->Branch("GFresidualLP_p", &event.GFresidual_p_lp2);
  tree->Branch("GFresidual6DLP_x", &event.GFresidual6D_x_lp2);
  tree->Branch("GFresidual6DLP_y", &event.GFresidual6D_y_lp2);
  tree->Branch("GFresidual6DLP_z", &event.GFresidual6D_z_lp2);
  tree->Branch("GFresidual6DLP_px", &event.GFresidual6D_px_lp2);
  tree->Branch("GFresidual6DLP_py", &event.GFresidual6D_py_lp2);
  tree->Branch("GFresidual6DLP_pz", &event.GFresidual6D_pz_lp2);
  tree->Branch("GFresolutionLP_x", &event.GFresolution_x_lp2);
  tree->Branch("GFresolutionLP_y", &event.GFresolution_y_lp2);
  tree->Branch("GFresolutionLP_z", &event.GFresolution_z_lp2);
  tree->Branch("GFresolutionLP_p", &event.GFresolution_p_lp2);
  tree->Branch("GFresolutionLP_px", &event.GFresolution_px_lp2);
  tree->Branch("GFresolutionLP_py", &event.GFresolution_py_lp2);
  tree->Branch("GFresolutionLP_pz", &event.GFresolution_pz_lp2);
  tree->Branch("GFpullLP_x", &event.GFpull_x_lp2);
  tree->Branch("GFpullLP_y", &event.GFpull_y_lp2);
  tree->Branch("GFpullLP_z", &event.GFpull_z_lp2);
  tree->Branch("GFpullLP_p", &event.GFpull_p_lp2);
  tree->Branch("GFpullLP_px", &event.GFpull_px_lp2);
  tree->Branch("GFpullLP_py", &event.GFpull_py_lp2);
  tree->Branch("GFpullLP_pz", &event.GFpull_pz_lp2);
  ***/
  
  // K-
  tree->Branch("Kmflag", &event.kmflag);
  tree->Branch("KmDecayflag", &event.kmdecayflag);
  tree->Branch("KmMomMatchflag", &event.kmmommatch);
  tree->Branch("KmMomDiff", &event.kmmomdiff);      
  tree->Branch("KmNTrack", &event.kmntrack);  
  tree->Branch("KmMass", &event.kmmass);
  tree->Branch("KmMassSquare", &event.kmmass2);
  tree->Branch("KmDecayVtx_x", &event.kmdecayvtx_x);
  tree->Branch("KmDecayVtx_y", &event.kmdecayvtx_y);
  tree->Branch("KmDecayVtx_z", &event.kmdecayvtx_z);
  tree->Branch("KmMom", &event.kmmom);  
  tree->Branch("KmMom_x", &event.kmmom_x);
  tree->Branch("KmMom_y", &event.kmmom_y);
  tree->Branch("KmMom_z", &event.kmmom_z);
  tree->Branch("KmScatPMom", &event.kmscatpmom);  
  tree->Branch("KmScatPMom_x", &event.kmscatpmom_x);
  tree->Branch("KmScatPMom_y", &event.kmscatpmom_y);
  tree->Branch("KmScatPMom_z", &event.kmscatpmom_z);  
  tree->Branch("KmTrackId", &event.km_id);  
  tree->Branch("KmDecaysTrackId", &event.kmdecays_id);  
  
#if DoKinematicFitLd
  // tree->Branch("KFchisqrXi",&event.KFchisqrxi);
  // tree->Branch("KFpvalXi",&event.KFpvalxi);
  // tree->Branch("KFXimom",&event.KFximom);
  // tree->Branch("KFXimom_x",&event.KFximom_x);
  // tree->Branch("KFXimom_y",&event.KFximom_y);
  // tree->Branch("KFXimom_z",&event.KFximom_z);
  /***
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
  ***/
  // tree->Branch("XiAccidentals", &event.XiAccidental); 
  // tree->Branch("XiFlight", &event.XiFlight);
  // tree->Branch("XiProd", &event.XiProd);
  // tree->Branch("xtgtXi", &event.xtgtXi);
  // tree->Branch("ytgtXi", &event.ytgtXi);
  // tree->Branch("utgtXi", &event.utgtXi);
  // tree->Branch("vtgtXi", &event.vtgtXi);
  // tree->Branch("vtxKKXi", &event.vtxKKXi);
  // tree->Branch("vtyKKXi", &event.vtyKKXi);
  // tree->Branch("vtzKKXi", &event.vtzKKXi);
  // tree->Branch("xiprodvtx_x", &event.xiprodvtx_x);
  // tree->Branch("xiprodvtx_y", &event.xiprodvtx_y);
  // tree->Branch("xiprodvtx_z", &event.xiprodvtx_z);
  // tree->Branch("xiprodmom_x", &event.xiprodmom_x);
  // tree->Branch("xiprodmom_y", &event.xiprodmom_y);
  // tree->Branch("xiprodmom_z", &event.xiprodmom_z);


#endif
  TTreeReaderCont[kE42] = new TTreeReader( "tpc", TFileCont[kE42] );
  const auto& reader = TTreeReaderCont[kE42];
  src.runnum = new TTreeReaderValue<Int_t>( *reader, "runnum" );
  src.evnum = new TTreeReaderValue<Int_t>( *reader, "evnum" );
  src.trigpat = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trigpat" );
  src.trigflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trigflag" );

  src.nhHtof = new TTreeReaderValue<Int_t>( *reader, "nhHtof" );
  src.HtofSeg = new TTreeReaderValue<std::vector<Double_t>>( *reader, "HtofSeg" );
  src.tHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "tHtof" );
  src.dtHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "dtHtof" );
  src.deHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "deHtof" );
  src.posHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "posHtof" );

  src.ntK18 = new TTreeReaderValue<Int_t>( *reader, "ntK18" );
  src.pK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pK18" );
  src.thetaK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaK18" );  
  src.chisqrK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrK18" );
  src.xtgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtK18" );
  src.ytgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtK18" );
  src.utgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtK18" );
  src.vtgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtK18" );
  // src.xvpHS  =  new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "xvpHS" );
  // src.yvpHS  =  new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "yvpHS" );
  // src.zvpHS  =  new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "zvpHS" );
  // src.xtgtHS = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtHS" );
  // src.ytgtHS = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtHS" );
  // src.ztgtHS = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ztgtHS" );

  src.ntKurama = new TTreeReaderValue<Int_t>( *reader, "ntKurama" );
  src.chisqrKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrKurama" );
  src.pKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pKurama" );
  src.qKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "qKurama" );
  //  src.m2 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "m2" );
  src.xtgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtKurama" );
  src.ytgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtKurama" );
  src.utgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtKurama" );
  src.vtgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtKurama" );

  src.nKm = new TTreeReaderValue<Int_t>( *reader, "nKm" );
  src.nKp = new TTreeReaderValue<Int_t>( *reader, "nKp" );
  src.nKK = new TTreeReaderValue<Int_t>( *reader, "nKK" );
  src.inside = new TTreeReaderValue<std::vector<Int_t>>( *reader, "inside" );
  src.vtx = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx" );
  src.vty = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vty" );
  src.vtz = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtz" );
  src.closeDist = new TTreeReaderValue<std::vector<Double_t>>( *reader, "closeDist" );
  src.MissMass = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMass" );
  src.MissMassCorr = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassCorr" );
  src.MissMassCorrDE = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassCorrDE" );
  src.Kflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "Kflag" );
  src.Pflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "Pflag" );

  src.ntTpc = new TTreeReaderValue<Int_t>( *reader, "ntTpc" );
  src.nhtrack = new TTreeReaderValue<std::vector<Int_t>>( *reader, "nhtrack" );
  src.trackid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trackid" );
  src.isBeam = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isBeam" );
  src.isKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isKurama" );
  src.isK18 = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isK18" );
  src.isAccidental = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isAccidental" );
  src.charge = new TTreeReaderValue<std::vector<Int_t>>( *reader, "charge" );
  src.pid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "pid" );
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
  src.isElectron = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isElectron" );
  src.nsigma_triton = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_triton" );
  src.nsigma_deutron = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_deutron" );
  src.nsigma_proton = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_proton" );
  src.nsigma_kaon = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_kaon" );
  src.nsigma_pion = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_pion" );
  src.nsigma_electron = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_electron" );  
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
  src.track_cluster_mrow = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_mrow" );

  src.ntTPCK18 = new TTreeReaderValue<Int_t>( *reader, "ntK18" );
  src.chisqrK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrK18" );
  //  src.pK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pK18");
  src.xtgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtK18" );
  src.ytgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtK18" );
  src.utgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtK18" );
  src.vtgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtK18" );
#if SaveTPCK18
  src.isgoodTPCK18 = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isgoodTPCK18" );
  src.chisqrTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrTPCK18" );
  src.qTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "qTPCK18");
  src.pTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pTPCK18");
  src.thetaTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaTPCK18");  
  src.xtgtTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtTPCK18" );
  src.ytgtTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtTPCK18" );
  src.utgtTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtTPCK18" );
  src.vtgtTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtTPCK18" );
  src.lhtofTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "lhtofTPCK18" );
  src.xhtofTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xhtofTPCK18" );
  src.yhtofTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yhtofTPCK18" );
  src.lvpTPCK18 = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "lvpTPCK18" );
  src.xvpTPCK18 = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "xvpTPCK18" );
  src.yvpTPCK18 = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "yvpTPCK18" );
  // src.xhtofK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xhtofHS" );
  // src.yhtofK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yhtofHS" );
#endif

  src.ntTPCKurama = new TTreeReaderValue<Int_t>( *reader, "ntKurama" );
  src.chisqrKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrKurama" );
  src.pKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pKurama" );
  src.qKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "qKurama" );
  //  src.m2  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "m2" );
  src.xtgtKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtKurama" );
  src.ytgtKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtKurama" );
  src.utgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtKurama" );
  src.vtgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtKurama" );
  src.tpcidTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "tpcidTPCKurama" );
  src.isgoodTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isgoodTPCKurama" );
  src.kflagTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "kflagTPCKurama" );
  src.pflagTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "pflagTPCKurama" );
  src.chisqrTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrTPCKurama" );
  src.pTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pTPCKurama" );
  src.qTPCKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "qTPCKurama" );
  src.m2TPCKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "m2TPCKurama" );
  src.xtgtTPCKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtTPCKurama" );
  src.ytgtTPCKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtTPCKurama" );
  src.utgtTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtTPCKurama" );
  src.vtgtTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtTPCKurama" );
  src.thetaTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaTPCKurama" );
  src.pathTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pathTPCKurama" );
  src.lhtofTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "lhtofTPCKurama" );
  src.xhtofTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xhtofTPCKurama" );
  src.yhtofTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yhtofTPCKurama" );
  src.lvpTPCKurama = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "lvpTPCKurama" );
  src.xvpTPCKurama = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "xvpTPCKurama" );
  src.yvpTPCKurama = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "yvpTPCKurama" );

  src.isgoodTPC = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isgoodTPC" );
  src.insideTPC = new TTreeReaderValue<std::vector<Int_t>>( *reader, "insideTPC" );
  src.vtxTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtxTPC" );
  src.vtyTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtyTPC" );
  src.vtzTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtzTPC" );
  src.closeDistTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "closeDistTPC" );
  src.MissMassTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassTPC" );
  src.MissMassCorrTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassCorrTPC" );
  src.MissMassCorrDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassCorrDETPC" );
  src.MissMassNuclTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassNuclTPC" );
  src.MissMassNuclCorrTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassNuclCorrTPC" );
  src.MissMassNuclCorrDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassNuclCorrDETPC" );
  src.pOrgTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pOrgTPC" );
  src.pCalcTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCalcTPC" );
  src.pCorrTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCorrTPC" );
  src.pCorrDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCorrDETPC" );
  src.thetaTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaTPC" );
  //src.thetaCMTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaCMTPC" );
  //src.costCMTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "costCMTPC" );
  src.ubTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ubTPC" );
  src.vbTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vbTPC" );
  src.usTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "usTPC" );
  src.vsTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vsTPC" );
  src.nvtxTpc = new TTreeReaderValue<Int_t>(*reader,"nvtxTpc");
  src.vtx_x = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx_x" );
  src.vtx_y = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx_y" );
  src.vtx_z = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx_z" );
  src.vtx_dist = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx_dist" );
  src.vtx_angle = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx_angle" );
  src.vtxid = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxid" );
  src.vtxmom_theta = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxmom_theta" );
  src.vtxpos_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxpos_x" );
  src.vtxpos_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxpos_y" );
  src.vtxpos_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxpos_z" );
  src.vtxmom_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxmom_x" );
  src.vtxmom_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxmom_y" );
  src.vtxmom_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxmom_z" );

  src.angleLambda = new TTreeReaderValue<std::vector<Double_t>>( *reader, "angleLambda" );

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

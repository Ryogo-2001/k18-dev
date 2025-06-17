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
#include "FourVectorCartesianFitter.hh"
#include "FourVectorFitter.hh"

#include "HypTPCFitter.hh"
#include "HypTPCTask.hh"

#define SaveTPCK18 0
#define MMcut_Xi 1
#define DoKinematicFitLdXi 1

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

const Double_t xi_masscut = 0.15;
const Double_t lambda_masscut = 0.1;
//const Double_t xi_masscut = 0.03; const Double_t lambda_masscut = 0.03;
const Double_t p_vtx_distcut = 300;
const Double_t pi_vtx_distcut = 300;
const Double_t pi2_vtx_distcut = 300;
//const Double_t xi_vtx_distcut = 100;
const Double_t ppi_distcut = 10.;
const Double_t lpi_distcut = 10.;
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
const Double_t GFxitarget_distcut = 50.;
const Double_t GFxitarget_ycut = 20.;

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

  //TPC RK
  std::vector<Int_t> isgoodTPCK18;
  std::vector<Int_t> tpcidTPCK18;
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
  std::vector<Int_t> tpcidTPCKurama;
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

  bool XiAccidental;
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
    xtgtKurama.clear();
    ytgtKurama.clear();
    utgtKurama.clear();
    vtgtKurama.clear();

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

    isgoodTPCK18.clear();
    tpcidTPCK18.clear();
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
    tpcidTPCKurama.clear();
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

    XiAccidental = false;
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
  TTreeReaderValue<std::vector<Double_t>>* chisqrK18;
  TTreeReaderValue<std::vector<Double_t>>* xtgtK18;
  TTreeReaderValue<std::vector<Double_t>>* ytgtK18;
  TTreeReaderValue<std::vector<Double_t>>* utgtK18;
  TTreeReaderValue<std::vector<Double_t>>* vtgtK18;

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

  TTreeReaderValue<Int_t>* ntTPCK18; // Number of Tracks
  TTreeReaderValue<std::vector<Int_t>>* tpcidTPCK18;
  TTreeReaderValue<std::vector<Int_t>>* isgoodTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* chisqrTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* qTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* pTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* xtgtTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* ytgtTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* utgtTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* vtgtTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* thetaTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* lhtofTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* xhtofTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* yhtofTPCK18;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* lvpTPCK18;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* xvpTPCK18;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* yvpTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* xhtofK18;
  TTreeReaderValue<std::vector<Double_t>>* yhtofK18;

  TTreeReaderValue<Int_t>* ntTPCKurama; // Number of Tracks
  TTreeReaderValue<std::vector<Int_t>>* tpcidTPCKurama;
  TTreeReaderValue<std::vector<Int_t>>* isgoodTPCKurama;
  TTreeReaderValue<std::vector<Int_t>>* kflagTPCKurama;
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
  //if( ievent%1000==0 ){
  if( ievent%100==0 ){
    std::cout << "#D Event Number: "
        << std::setw(6) << ievent << std::endl;
  }

  TVector3 tgtpos(0, 0, tpc::ZTarget);

  GetEntry(ievent);

  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;

  event.nKK = **src.nKK;
  event.Kflag = **src.Kflag;
  event.MissMass = **src.MissMass;
  event.MissMassCorr = **src.MissMassCorr;
  event.MissMassCorrDE = **src.MissMassCorrDE;

  event.vtx = **src.vtx;
  event.vty = **src.vty;
  event.vtz = **src.vtz;
  event.closeDist = **src.closeDist;

  event.ntK18 = **src.ntK18;
  event.chisqrK18 = **src.chisqrK18;
  event.pK18 = **src.pK18;
  event.xtgtK18 = **src.xtgtK18;
  event.ytgtK18 = **src.ytgtK18;
  event.utgtK18 = **src.utgtK18;
  event.vtgtK18 = **src.vtgtK18;
#if SaveTPCK18
  event.isgoodTPCK18 = **src.isgoodTPCK18;
  event.chisqrTPCK18 = **src.chisqrTPCK18;
  event.pTPCK18 = **src.pTPCK18;
  event.qTPCK18 = **src.qTPCK18;
  event.xtgtTPCK18 = **src.xtgtTPCK18;
  event.ytgtTPCK18 = **src.ytgtTPCK18;
  event.utgtTPCK18 = **src.utgtTPCK18;
  event.vtgtTPCK18 = **src.vtgtTPCK18;
  event.thetaTPCK18 = **src.thetaTPCK18;
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
  event.pOrgTPC = **src.pOrgTPC;
  event.pCalcTPC = **src.pCalcTPC;
  event.pCorrTPC = **src.pCorrTPC;
  event.pCorrDETPC = **src.pCorrDETPC;
  event.thetaTPC = **src.thetaTPC;
  event.ubTPC = **src.ubTPC;
  event.vbTPC = **src.vbTPC;
  event.usTPC = **src.usTPC;
  event.vsTPC = **src.vsTPC;

  if( event.nKK != 1 || event.Kflag[0]!=1 ) return true;
  HF1( 4, event.MissMassTPC[0] );
  HF1( 5, event.MissMass[0] );

#if MMcut_Xi
  if( event.MissMass[0] > XiMinusMass + 0.1 || event.MissMass[0] < XiMinusMass - 0.1 ){
    std::cout<<"cut by MM"<<std::endl;
    return true;
  }
#endif
  HF1( 8, event.MissMassTPC[0] );
  HF1( 9, event.MissMass[0] );

  Int_t KuramaTrackId = -1; Int_t K18TrackId = -1;
  for(Int_t it=0;it<event.ntKurama;it++){
    if(event.isgoodTPCKurama[it]==1){
      KuramaTrackId = event.tpcidTPCKurama[it];
    }
  }
#if SaveTPCK18
  for(Int_t it=0;it<event.ntK18;it++){
    if(event.isgoodTPCKurama[it]==1){
      K18TrackId = event.tpcidTPCK18[it];
    }
  }
#endif

  HF1( 1, event.status++ );

  event.nhHtof = **src.nhHtof;
  event.HtofSeg = **src.HtofSeg;
  event.tHtof = **src.tHtof;
  event.dtHtof = **src.dtHtof;
  event.deHtof = **src.deHtof;
  event.posHtof = **src.posHtof;

  Int_t ntTpc = **src.ntTpc;
  if( ntTpc == 0 ) return true;
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
  event.track_cluster_mrow = **src.track_cluster_mrow;

  event.ntTpc = **src.ntTpc;
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

  //Xi candidates searching
  Int_t xi_candidates = 0;

  Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  std::vector<Int_t> xi_p_container, xi_pi_container, xi_pi2_container;
  std::vector<TVector3> xi_mom_container, xi_vert_container;
  std::vector<TVector3> l_mom_container, l_vert_container;
  std::vector<Double_t> xi_mass_container, lambda_mass_container;
  std::vector<TVector3> p_mom_container, pi_mom_container, pi2_mom_container;
  std::vector<Double_t> ppi_closedist; std::vector<Double_t> lpi_closedist;

  for(Int_t it1=0;it1<ntTpc;it1++){
    //std::cout<<"it1 "<<it1<<" pid "<<event.pid[it1]<<" charge "<<event.charge[it1]<<std::endl;
    if(event.isK18[it1]==1 || event.isKurama[it1]==1) continue;
    //if(event.isK18[it1]==1) continue;
    //if(event.isK18[it1]==1 || event.isgoodTPCKurama[it1]==1) continue;
    //if((event.pid[it1]&4)!=4 && event.charge[it1]==1) continue; //select proton like
    //if((event.pid[it1]&4)!=4 || event.charge[it1]!=1) continue; //select proton like
    //if(event.charge[it1]==1 && ((event.pid[it1]&4)==4)) p_like = true;
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
      //std::cout<<"!p_like"<<std::endl;
    }

    for(Int_t it2=0;it2<ntTpc;it2++){
      if(it1==it2) continue;
      if(event.isK18[it2]==1 || event.isKurama[it2]==1) continue;
      //if(event.isK18[it2]==1 || event.isgoodTPCKurama[it2]==1) continue;
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
      //std::cout<<"recon L vtx "<<lambda_vert<<" id "<<it1<<" "<<it2<<std::endl;
      if(TMath::Abs(lambda_vert.x()) > 250. || TMath::Abs(lambda_vert.z()) > 250. || TMath::Abs(lambda_vert.y()) > 250.) continue; //Vertex cut
      //std::cout<<"ppi_dist "<<ppi_dist<<std::endl;
      if(ppi_dist > ppi_distcut) continue;
      Double_t pi_vertex_dist; Double_t p_vertex_dist;
      /*
	Kinematics::HelixDirection(lambda_vert, p_start, p_end, p_vertex_dist);
	Kinematics::HelixDirection(lambda_vert, pi_start, pi_end, pi_vertex_dist);
      */
      //std::cout<<"before helix direction"<<std::endl;
      if(!Kinematics::HelixDirection(lambda_vert, p_start, p_end, p_vertex_dist) ||
	 !Kinematics::HelixDirection(lambda_vert, pi_start, pi_end, pi_vertex_dist)) continue;

      //std::cout<<"p, pi Lvertex dist "<<p_vertex_dist<<" pi "<<pi_vertex_dist<<std::endl;
      if(pi_vertex_dist > pi_vtx_distcut) continue;
      if(p_vertex_dist > p_vtx_distcut) continue;

      //std::cout<<"before l mass cut "<<TMath::Abs(Llambda.M() - LambdaMass)<<std::endl;
      if(TMath::Abs(Llambda.M() - LambdaMass) > lambda_masscut) continue;
      event.lflag = true;
      //std::cout<<"!!!recon L mass "<<Llambda.M()<<" diff "<<TMath::Abs(Llambda.M() - LambdaMass)<<" id "<<it1<<" "<<it2<<" "<<std::endl;
      //std::cout<<"!!! L reconstructed!!! "<<it1<<" "<<it2<<" "<<std::endl;
      for(Int_t it3=0;it3<ntTpc;it3++){
	if(it3==it2 || it3==it1) continue;
	if(event.isK18[it3]==1 || event.isKurama[it3]==1) continue;
	//if(event.isK18[it3]==1 || event.isgoodTPCKurama[it3]==1) continue;
	if((event.pid[it3]&1)!=1) continue; //select pi like
	Bool_t pim_like2 = false;
	if(event.charge[it3]==-1) pim_like2 = true;

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
	/*
	  Double_t xi_vertex_dist;
	  Kinematics::HelixDirection(xi_vert, p_start, p_end, xi_vertex_dist);
	  Kinematics::HelixDirection(xi_vert, pi2_start, pi2_end, pi2_vertex_dist);
	*/
	/*
	  std::cout<<"before helixdirection "<<std::endl;
	  if(!Kinematics::HelixDirection(xi_vert, p_start, p_end, xi_vertex_dist)) continue;
	  if(xi_vertex_dist > xi_vtx_distcut) continue;
	  std::cout<<"<<xi_vertex_dist<<" vertex "<<xi_vert<<std::endl;
	*/
	if(!Kinematics::HelixDirection(xi_vert, pi2_start, pi2_end, pi2_vertex_dist)) continue;
	//std::cout<<"pi2 start "<<pi2_start<<" end "<<pi2_end<<std::endl;
	//std::cout<<"pi2_vertex_dist "<<pi2_vertex_dist<<" vertex "<<xi_vert<<std::endl;
	if(pi2_vertex_dist > pi2_vtx_distcut) continue;
	//std::cout<<"before xi masscut "<<TMath::Abs(Lxi.M() - XiMinusMass)<<std::endl;
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
	//std::cout<<"!!! recon Xi mass "<<Lxi.M()<<" id "<<it1<<" "<<it2<<" "<<it3<<std::endl;
      } //it3
    } //it2
  } //it1

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

  if( event.lflag) HF1( 7, event.MissMass[0] );

  if( !event.xiflag ) return true;

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
  event.decays_mom.push_back(p_mom_container[best].Mag());
  event.decays_mom_x.push_back(p_mom_container[best].x());
  event.decays_mom_y.push_back(p_mom_container[best].y());
  event.decays_mom_z.push_back(p_mom_container[best].z());

  event.decays_id.push_back(xi_pi_container[best]);
  event.decays_mom.push_back(pi_mom_container[best].Mag());
  event.decays_mom_x.push_back(pi_mom_container[best].x());
  event.decays_mom_y.push_back(pi_mom_container[best].y());
  event.decays_mom_z.push_back(pi_mom_container[best].z());

  event.decays_id.push_back(xi_pi2_container[best]);
  event.decays_mom.push_back(pi2_mom_container[best].Mag());
  event.decays_mom_x.push_back(pi2_mom_container[best].x());
  event.decays_mom_y.push_back(pi2_mom_container[best].y());
  event.decays_mom_z.push_back(pi2_mom_container[best].z());

  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCTracks(**src.ntTpc, **src.isK18, **src.isKurama,
			 **src.charge, **src.nhtrack, **src.helix_cx,
			 **src.helix_cy, **src.helix_z0, **src.helix_r,
			 **src.helix_dz, **src.hitlayer, **src.track_cluster_mrow,
			 **src.helix_t, **src.track_cluster_de, **src.resolution_x,
			 **src.resolution_y, **src.resolution_z, **src.hitpos_x,
			 **src.hitpos_y, **src.hitpos_z);
  HF1( 1, event.status++ );

  Int_t id_p = xi_p_container[best];
  Int_t id_pi = xi_pi_container[best];
  Int_t id_pi2 = xi_pi2_container[best];
  auto Track_p = TPCAna.GetTrackTPCHelix(id_p);
  auto Track_pi = TPCAna.GetTrackTPCHelix(id_pi);
  auto Track_pi2 = TPCAna.GetTrackTPCHelix(id_pi2);
  for(Int_t ih=0;ih<Track_p->GetNHit();++ih){
    auto pos = Track_p->GetHitInOrder(ih)->GetLocalHitPos();
    HF2(1001,pos.z(),pos.x());
  }
  for(Int_t ih=0;ih<Track_pi->GetNHit();++ih){
    auto pos = Track_pi->GetHitInOrder(ih)->GetLocalHitPos();
    HF2(1002,pos.z(),pos.x());
  }
  for(Int_t ih=0;ih<Track_pi2->GetNHit();++ih){
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



  HypTPCTask& GFTrackCont = HypTPCTask::GetInstance();

  HF1( 10, event.MissMass[0] );
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
  TVector3 KFlambda_mom = TVector3(KFHLVLd.X(),KFHLVLd.Z(),KFHLVLd.Y());  
  GFlambda_mom = KFlambda_mom;   
#endif
    TLorentzVector GFLlambda_fixed(GFlambda_mom, TMath::Sqrt(GFlambda_mom.Mag()*GFlambda_mom.Mag() + LambdaMass*LambdaMass));

    TVector3 GFxi_vert; Double_t GFlpi_dist = qnan; Double_t GFlambda_tracklen;
    #if DoKinematicFitLdXi
    double l_res_x,l_res_y,l_phi;
    MathTools::DecomposeResolution(VLd,KFlambda_mom,l_res_x,l_res_y,l_phi); 
    if(!GFTrackCont.FindVertexXi(trackid_pi2, repid_pi2, GFlambda_vert, GFlambda_mom,
    GFlambda_tracklen, GFextrapolation_decays[2], GFmom_decays[2], GFlpi_dist, GFxi_vert, vtx_scan_range,l_res_x,l_res_y,l_phi) 
      || GFlpi_dist > GFlpi_distcut) continue;
    #else
    if(!GFTrackCont.FindVertexXi(trackid_pi2, repid_pi2,
         GFlambda_vert, GFlambda_mom, GFlambda_tracklen,
         GFextrapolation_decays[0], GFmom_decays[2],
         GFlpi_dist, GFxi_vert,
         vtx_scan_range)
       || GFlpi_dist > GFlpi_distcut) continue;
    #endif

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
    GFxi_vert_container[candi] = GFxi_vert;
    GFl_vert_container[candi] = GFlambda_vert;
    GFxi_mass_container[candi] = GFLxi.M();
    GFl_mass_container[candi] = GFLlambda.M();
    GFp_mass_container[candi] = GFmass2_decays[0];
    GFpi_mass_container[candi] = GFmass2_decays[1];
    GFpi2_mass_container[candi] = GFmass2_decays[2];
    GFppi_closedist_container[candi] = GFppi_dist;
    GFlpi_closedist_container[candi] = GFlpi_dist;

    Double_t diff = TMath::Abs(GFLxi.M() - XiMinusMass);
    if(prev_massdiff > diff){
      prev_massdiff = diff;
      gfbest = candi;
    }
  }//

  if(!event.GFxiflag) return true;
  HF1( genfitHid, GFntTpc);

  HF1( 20, event.MissMass[0] );
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
  for( Int_t j=0; j<3; ++j ){
    Int_t igf = GFxi_p_id_container[gfbest];
    if(j==1) igf = GFxi_pi_id_container[gfbest];
    if(j==2) igf = GFxi_pi2_id_container[gfbest];

    Int_t repid = GFxi_p_rep_container[gfbest];
    if(j==1) repid = GFxi_pi_rep_container[gfbest];
    if(j==2) repid = GFxi_pi2_rep_container[gfbest];

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
    if(j==2) GFmom_decays = GFpi2_mom_container[gfbest];

    Double_t GFmass_decays = GFp_mass_container[gfbest];
    if(j==1) GFmass_decays = GFpi_mass_container[gfbest];
    if(j==2) GFmass_decays = GFpi2_mass_container[gfbest];

    event.GFdecays_m2[j] = GFmass_decays;
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

    Int_t id = GFxi_p_id_container[gfbest];
    if(j==1) id = GFxi_pi_id_container[gfbest];
    if(j==2) id = GFxi_pi2_id_container[gfbest];

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
  } //igf
TVector3 KKVert(event.vtx[0],event.vty[0],event.vtz[0]);
TVector3 XiMom = GFxi_mom_container[gfbest];
TVector3 XiVert = GFxi_vert_container[gfbest];
if(event.isgoodTPC[0] == 1){
  KKVert = TVector3(event.vtxTPC[0],event.vtyTPC[0],event.vtzTPC[0]);
}
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

  TVector3 HTVP(event.GFdecays_mom_x[0],event.GFdecays_mom_z[0],event.GFdecays_mom_y[0]);
  TVector3 HTVPi1(event.GFdecays_mom_x[1],event.GFdecays_mom_z[1],event.GFdecays_mom_y[1]);
  TVector3 HTVLd = HTVP+HTVPi1;
  TLorentzVector HLVP(HTVP,hypot(ProtonMass,HTVP.Mag()));
  TLorentzVector HLVPi1(HTVPi1,hypot(PionMass,HTVPi1.Mag()));
  TLorentzVector HLVLd(HTVLd,hypot(LambdaMass,HTVLd.Mag()));

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
  if(TMath::Abs(XiTgtVert.y()) > GFxitarget_ycut) event.XiAccidental = true;
  const int ntrack = 3;
  double x0track[ntrack]={event.xtgtTPCKurama[0],event.xtgtK18[0],event.xtgtXi};
  double y0track[ntrack]={event.ytgtTPCKurama[0],event.ytgtK18[0],event.ytgtXi};
  double u0track[ntrack]={event.utgtTPCKurama[0],event.utgtK18[0],event.utgtXi};
  double v0track[ntrack]={event.vtgtTPCKurama[0],event.vtgtK18[0],event.vtgtXi};
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
  if(XiProd){
    event.xiprodvtx_x = XiProdVert.x();
    event.xiprodvtx_y = XiProdVert.y();
    event.xiprodvtx_z = XiProdVert.z();
    event.xiprodmom_x = XiProdMom.x();
    event.xiprodmom_y = XiProdMom.y();
    event.xiprodmom_z = XiProdMom.z();
    HF1(100,XiMom.Mag()-XiProdMom.Mag());
  }
  HF1( 30, event.MissMass[0]);
  HF1( 31, event.GFlmass);
  HF1( 32, event.GFximass);
  HF1( 33, event.GFdecays_mom[0]);
  HF1( 34, event.GFdecays_mom[1]);
  HF1( 35, event.GFdecays_mom[2]);
  HF1( 36, event.GFlmom);
  HF1( 37, event.GFximom);
  HF1( 38, event.GFppi_dist);
  HF1( 39, event.GFlpi_dist);
  HF1( 40, event.GFdecays_m2[0]);
  HF1( 41, event.GFdecays_m2[1]);
  HF1( 42, event.GFdecays_m2[2]);
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
  HB1( 4, "Missing Mass [TPC+KURAMA]; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 160, 1., 1.8);
  HB1( 5, "Missing Mass [KURAMA]; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 160, 1., 1.8);
  HB1( 7, "Missing Mass, #L tagged [KURAMA]; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 160, 1., 1.8);
  HB1( 8, "Missing Mass [TPC+KURAMA]; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 160, 1., 1.8);
  HB1( 9, "Missing Mass [KURAMA]; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 160, 1., 1.8);
  HB1( 10, "Missing Mass, #Xi^{-} tagged [KURAMA]; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 160, 1., 1.8);
  HB1( 11," #Lambda Invariant Mass; M_{p#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 80, 1.04, 1.2);
  HB1( 12," #Xi Invariant Mass; M_{#Lambda#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 125, 1.2 ,1.45);
  HB1( 13, "p Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 14, "#pi Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  HB1( 15, "#pi Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 50, 0., 1.0);
  HB1( 16, "#Lambda Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 17, "#Xi Momentum; Momentum [GeV/#font[12]{c}]; Counts [/0.02 GeV/#font[12]{c}]", 100, 0., 2.0);
  HB1( 18, "Closest Dist. p#pi;Distance [mm]; Counts [/0.1 mm]", 1000, 0., 100.);
  HB1( 19, "Closest Dist. #Lambda#pi;Distance [mm]; Counts [/0.1 mm]", 1000, 0., 100.);

  HB1( 20, "Missing Mass [KURAMA]; Missing mass [GeV/c^{2}]; Counts [/5 MeV/#font[12]{c}^{2}]", 160, 1., 1.8);
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
  HB1( 40, "[GenFit] p_{#Lambda} Mass; M^{2} [(GeV/#font[12]{c}^{2})^{2}]; Counts [/0.002 (GeV/#font[12]{c}^{2})^{2}]", 750, 0., 1.5);
  HB1( 41, "[GenFit] pi_{#Lambda} Mass; M^{2} [(GeV/#font[12]{c}^{2})^{2}]; Counts [/0.002 (GeV/#font[12]{c}^{2})^{2}]", 200, -0.1, 0.3);
  HB1( 42, "[GenFit] pi_{#Xi} Mass; M^{2} [(GeV/#font[12]{c}^{2})^{2}]; Counts [/0.002 (GeV/#font[12]{c}^{2})^{2}]", 200, 0.1, 0.3);


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

  HB2(1001, "#Xi^{-} decay, p hit pattern",100,-250,250,100,-250,250);
  HB2(1002, "#Xi^{-} decay, #pi_{#Lambda} hit pattern",100,-250,250,100,-250,250);
  HB2(1003, "#Xi^{-} decay, #pi_{#Xi} hit pattern",100,-250,250,100,-250,250);
  HB2(2001, "#Xi^{-} decay, p hit patternGF",100,-250,250,100,-250,250);
  HB2(2002, "#Xi^{-} decay, #pi_{#Lambda} hit patternGF",100,-250,250,100,-250,250);
  HB2(2003, "#Xi^{-} decay, #pi_{#Xi} hit patternGF",100,-250,250,100,-250,250);

  HB1(3000, "#Xi Decay mom - Prod mom; #Delta p [GeV/#font[12]{c}]; Counts [/ 2MeV/#font[12]{c}]", 300, -0.3, 0.3);

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
  HB1(10021,"KF#{Lambda} residual p_{#theta}",1000,-0.1,0.1);
  HB1(10022,"KF#{Lambda} residual p_{#phi}",1000,-0.1,0.1);
  HB1(10023,"KF#{Lambda} residual #pi_{p}",1000,-0.3,0.3);
  HB1(10024,"KF#{Lambda} residual #pi_{#theta}",1000,-0.1,0.1);
  HB1(10025,"KF#{Lambda} residual #pi_{#phi}",1000,-0.1,0.1);


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
#endif


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
  tree->Branch( "helix_t", &event.helix_t );
  tree->Branch( "alpha", &event.alpha);
  tree->Branch( "pathhit", &event.pathhit);
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);

  tree->Branch( "ntK18", &event.ntK18);
  tree->Branch( "chisqrK18", &event.chisqrK18);
  tree->Branch( "pK18", &event.pK18);
  tree->Branch( "xtgtK18", &event.xtgtK18);
  tree->Branch( "ytgtK18", &event.ytgtK18);
  tree->Branch( "utgtK18", &event.utgtK18);
  tree->Branch( "vtgtK18", &event.vtgtK18);

#if SaveTPCK18
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

  tree->Branch("nKK", &event.nKK);
  tree->Branch("Kflag", &event.Kflag);
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
  tree->Branch( "pOrgTPC", &event.pOrgTPC);
  tree->Branch( "pCalcTPC", &event.pCalcTPC);
  tree->Branch( "pCorrTPC", &event.pCorrTPC);
  tree->Branch( "pCorrDETPC", &event.pCorrDETPC);

  tree->Branch("Lflag", &event.lflag);
  tree->Branch("Xiflag", &event.xiflag);
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

  //track fitting results
  tree->Branch("GFstatus", &event.GFstatus);
  tree->Branch("GFntTpc", &event.GFntTpc);
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
#if DoKinematicFitLdXi

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

  tree->Branch("XiAccidentals", &event.XiAccidental); 
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
  src.chisqrK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrK18" );
  src.xtgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtK18" );
  src.ytgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtK18" );
  src.utgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtK18" );
  src.vtgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtK18" );

  src.ntKurama = new TTreeReaderValue<Int_t>( *reader, "ntKurama" );
  src.chisqrKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrKurama" );
  src.pKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pKurama" );
  src.qKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "qKurama" );
  src.m2 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "m2" );
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
  src.pK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pK18");
  src.xtgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtK18" );
  src.ytgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtK18" );
  src.utgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtK18" );
  src.vtgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtK18" );
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

  src.ntTPCKurama = new TTreeReaderValue<Int_t>( *reader, "ntKurama" );
  src.chisqrKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrKurama" );
  src.pKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pKurama" );
  src.qKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "qKurama" );
  src.m2  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "m2" );
  src.xtgtKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtKurama" );
  src.ytgtKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtKurama" );
  src.utgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtKurama" );
  src.vtgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtKurama" );
  src.tpcidTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "tpcidTPCKurama" );
  src.isgoodTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isgoodTPCKurama" );
  src.kflagTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "kflagTPCKurama" );
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

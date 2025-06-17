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

#define SaveRawData 0
#define DebugDisp 0
#define SaveHistograms 0
#define kkevent 0

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
  const bool Const_field = false; //Must be false for linear tracking
  const Int_t verbosity = 0;//0~3;
  //const Int_t verbosity = 1;
  const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");

  const Double_t vtx_scan_range = 150.; //ref
  const Double_t vtx_scan_rangeInsideL = 50.;
  const Double_t vtx_scan_rangeInsidePi = 50.;

  const Double_t xi_masscut = 0.1; const Double_t lambda_masscut = 0.1; //final
  //const Double_t xi_masscut = 0.15; const Double_t lambda_masscut = 0.1; //ref
  const Double_t p_vtx_distcut = 300;
  const Double_t pi_vtx_distcut = 300;
  const Double_t pi2_vtx_distcut = 300;
  const Double_t e_vtx_distcut = 300;
  const Double_t pipi_distcut = 10.; //ref
  const Double_t ppi_distcut = 10.; //ref
  //const Double_t lpi_distcut = 10.;
  const Double_t lpi_distcut = 15.; //ref
  const Double_t phi_kk_distcut = 15.; //ref
  const Double_t xitarget_distcut = 50.; //ref
  const Double_t ltarget_distcut = 25.;

  const Double_t GFppi_distcut = 10.;
  const Double_t GFlpi_distcut = 10.;
  //const Double_t GFlpi_distcut = 15.;
  const Double_t GFphi_kk_distcut = 10.; //ref
  const Double_t GFltarget_distcut = 25.;
  const Double_t GFxitarget_ycut = 20.;
  const Double_t GFltarget_ycut = 20.;

  const Double_t residual_track_distcut = 25.;

  //For gamma reconstruction
  const Double_t gammatarget_distcut = 25.; //temp
  //const Double_t gamma_ecut = 0.5; //for Xi* decay pi0 has momentum < 0.35 GeV/c
  const Double_t gamma_ecut = 2.5; //for Xi* decay pi0 has momentum < 0.35 GeV/c
  const Double_t epem_distcut = 10.; //temp

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
  Double_t res_xXiVtx = 0.6*2;
  Double_t res_yXiVtx = 0.5*2;
  Double_t res_xLdVtx = 0.5*2;
  Double_t res_yLdVtx = 0.5*2;

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
  std::vector<Double_t> m2Kurama;
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

  Int_t nclTpc; // Number of clusters
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

  Int_t remain_nclTpc; // Number of remain clusters not occupied in the tracks
  std::vector<Double_t> remain_cluster_x;
  std::vector<Double_t> remain_cluster_y;
  std::vector<Double_t> remain_cluster_z;
  std::vector<Double_t> remain_cluster_de;
  std::vector<Int_t> remain_cluster_size;
  std::vector<Int_t> remain_cluster_layer;
  std::vector<Double_t> remain_cluster_mrow;
  std::vector<Double_t> remain_cluster_de_center;
  std::vector<Double_t> remain_cluster_x_center;
  std::vector<Double_t> remain_cluster_y_center;
  std::vector<Double_t> remain_cluster_z_center;
  std::vector<Int_t> remain_cluster_row_center;
  std::vector<Int_t> remain_cluster_houghflag;

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
  std::vector<std::vector<Double_t>> track_cluster_de_center;
  std::vector<std::vector<Double_t>> track_cluster_x_center;
  std::vector<std::vector<Double_t>> track_cluster_y_center;
  std::vector<std::vector<Double_t>> track_cluster_z_center;
  std::vector<std::vector<Double_t>> track_cluster_row_center;

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
  std::vector<Double_t> reconmassPipair;
  std::vector<Double_t> pipidistPipair;

  //Genfit Track's information
  std::vector<Int_t> GFfitstatus;
  std::vector<Int_t> GFpdgcode;
  std::vector<Int_t> GFnhtrack;
  std::vector<Double_t> GFcharge;
  std::vector<Double_t> GFchisqr;
  std::vector<Double_t> GFtof;
  //std::vector<Double_t> GFtracklen;
  std::vector<Double_t> GFpval;
  std::vector<std::vector<Double_t>> GFlayer;
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

  std::vector<Int_t> GFinside;

  Int_t GFntTpc_inside;
  Double_t GFprodvtx_x;
  Double_t GFprodvtx_y;
  Double_t GFprodvtx_z;

  std::vector<Double_t> GFtracklen;
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
  Double_t GFltargetcenter_dist1;
  Double_t GFltargetcenter_x1;
  Double_t GFltargetcenter_y1;
  Double_t GFltargetcenter_z1;

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
  Double_t GFltargetcenter_dist2;
  Double_t GFltargetcenter_x2;
  Double_t GFltargetcenter_y2;
  Double_t GFltargetcenter_z2;

  Double_t GFlprodvtx_x2;
  Double_t GFlprodvtx_y2;
  Double_t GFlprodvtx_z2;
  Double_t GFlprodvtx_dist2;
  Double_t GFltracklen2;
  Double_t GFltof2;

  Double_t GFlmass_alter1;
  Double_t GFldecayvtx_x_alter1;
  Double_t GFldecayvtx_y_alter1;
  Double_t GFldecayvtx_z_alter1;
  Double_t GFlmom_alter1;
  Double_t GFlmom_x_alter1;
  Double_t GFlmom_y_alter1;
  Double_t GFlmom_z_alter1;
  Double_t GFppi_dist_alter1;
  Double_t GFltarget_dist_alter1;
  Double_t GFltargetvtx_x_alter1;
  Double_t GFltargetvtx_y_alter1;
  Double_t GFltargetvtx_z_alter1;
  Double_t GFltargetcenter_dist_alter1;
  Double_t GFltargetcenter_x_alter1;
  Double_t GFltargetcenter_y_alter1;
  Double_t GFltargetcenter_z_alter1;

  Double_t GFlmass_alter2;
  Double_t GFldecayvtx_x_alter2;
  Double_t GFldecayvtx_y_alter2;
  Double_t GFldecayvtx_z_alter2;
  Double_t GFlmom_alter2;
  Double_t GFlmom_x_alter2;
  Double_t GFlmom_y_alter2;
  Double_t GFlmom_z_alter2;
  Double_t GFppi_dist_alter2;
  Double_t GFltarget_dist_alter2;
  Double_t GFltargetvtx_x_alter2;
  Double_t GFltargetvtx_y_alter2;
  Double_t GFltargetvtx_z_alter2;
  Double_t GFltargetcenter_dist_alter2;
  Double_t GFltargetcenter_x_alter2;
  Double_t GFltargetcenter_y_alter2;
  Double_t GFltargetcenter_z_alter2;

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
  Bool_t xipflag;
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
  Double_t phicosKK;
  Double_t phimom;
  Double_t phimom_x;
  Double_t phimom_y;
  Double_t phimom_z;
  Double_t kk_dist;
  Double_t GFphi_km_mass2;
  Double_t GFphi_km_invbeta;
  Double_t GFphi_kp_mass2;
  Double_t GFphi_kp_invbeta;

  Double_t GFphimass;
  Double_t GFphidecayvtx_x;
  Double_t GFphidecayvtx_y;
  Double_t GFphidecayvtx_z;
  Double_t GFphicosKK;
  Double_t GFphimom;
  Double_t GFphimom_x;
  Double_t GFphimom_y;
  Double_t GFphimom_z;
  Double_t GFkk_dist;
  Double_t GFphiprodvtx_dist;

  Double_t GFphimass_wKurama;
  Double_t GFphimom_wKurama;
  Double_t GFphimom_x_wKurama;
  Double_t GFphimom_y_wKurama;
  Double_t GFphimom_z_wKurama;
  std::vector<Int_t> phidecays_id;
  std::vector<Double_t> phidecays_mom;
  std::vector<Double_t> phidecays_mom_x;
  std::vector<Double_t> phidecays_mom_y;
  std::vector<Double_t> phidecays_mom_z;
  std::vector<Double_t> GFphidecays_mom;
  std::vector<Double_t> GFphidecays_mom_x;
  std::vector<Double_t> GFphidecays_mom_y;
  std::vector<Double_t> GFphidecays_mom_z;

  std::vector<Int_t> GFlldecays_pdgcode;
  std::vector<Int_t> GFlldecays_nhtrack;
  std::vector<Double_t> GFlldecays_charge;
  std::vector<Double_t> GFlldecays_chisqr;
  std::vector<Double_t> GFlldecays_pval;
  std::vector<Int_t> GFlldecays_htofid;
  std::vector<Double_t> GFlldecays_tracklen;
  std::vector<Double_t> GFlldecays_tof;
  std::vector<Double_t> GFlldecays_mass2;
  std::vector<Double_t> GFlldecays_invbeta;
  std::vector<Double_t> GFlldecays_mom;
  std::vector<Double_t> GFlldecays_mom_x;
  std::vector<Double_t> GFlldecays_mom_y;
  std::vector<Double_t> GFlldecays_mom_z;
  std::vector<Double_t> GFlldecays_CMmom;
  std::vector<Double_t> GFlldecays_CMmom_x;
  std::vector<Double_t> GFlldecays_CMmom_y;
  std::vector<Double_t> GFlldecays_CMmom_z;
  std::vector<Double_t> GFlldecays_momloss;
  std::vector<Double_t> GFlldecays_eloss;

  std::vector<Int_t> lldecays_id;
  std::vector<Double_t> lldecays_mom;
  std::vector<Double_t> lldecays_mom_x;
  std::vector<Double_t> lldecays_mom_y;
  std::vector<Double_t> lldecays_mom_z;
  std::vector<Double_t> lldecays_CMmom;
  std::vector<Double_t> lldecays_CMmom_x;
  std::vector<Double_t> lldecays_CMmom_y;
  std::vector<Double_t> lldecays_CMmom_z;

  std::vector<Int_t> GFxidecays_pdgcode;
  std::vector<Int_t> GFxidecays_nhtrack;
  std::vector<Double_t> GFxidecays_charge;
  std::vector<Double_t> GFxidecays_chisqr;
  std::vector<Double_t> GFxidecays_tracktof;
  std::vector<Double_t> GFxidecays_pval;
  std::vector<Int_t> GFxidecays_htofid;
  std::vector<Double_t> GFxidecays_tracklen;
  std::vector<Double_t> GFxidecays_tof;
  std::vector<Double_t> GFxidecays_mass2;
  std::vector<Double_t> GFxidecays_invbeta;
  std::vector<Double_t> GFxidecays_mom;
  std::vector<Double_t> GFxidecays_mom_x;
  std::vector<Double_t> GFxidecays_mom_y;
  std::vector<Double_t> GFxidecays_mom_z;
  std::vector<Double_t> GFxidecays_CMmom;
  std::vector<Double_t> GFxidecays_CMmom_x;
  std::vector<Double_t> GFxidecays_CMmom_y;
  std::vector<Double_t> GFxidecays_CMmom_z;
  std::vector<Double_t> GFxidecays_momloss;
  std::vector<Double_t> GFxidecays_eloss;

  std::vector<Int_t> xidecays_id;
  std::vector<Double_t> xidecays_mom;
  std::vector<Double_t> xidecays_mom_x;
  std::vector<Double_t> xidecays_mom_y;
  std::vector<Double_t> xidecays_mom_z;
  std::vector<Double_t> xidecays_CMmom;
  std::vector<Double_t> xidecays_CMmom_x;
  std::vector<Double_t> xidecays_CMmom_y;
  std::vector<Double_t> xidecays_CMmom_z;

  std::vector<Int_t> GFdecays_pdgcode;
  std::vector<Int_t> GFdecays_nhtrack;
  std::vector<Double_t> GFdecays_charge;
  std::vector<Double_t> GFdecays_chisqr;
  std::vector<Double_t> GFdecays_pval;
  std::vector<Int_t> GFdecays_htofid;
  std::vector<Double_t> GFdecays_tracklen;
  std::vector<Double_t> GFdecays_tof;
  std::vector<Double_t> GFdecays_mass2;
  std::vector<Double_t> GFdecays_invbeta;
  std::vector<Double_t> GFdecays_mom;
  std::vector<Double_t> GFdecays_mom_x;
  std::vector<Double_t> GFdecays_mom_y;
  std::vector<Double_t> GFdecays_mom_z;
  std::vector<Double_t> GFdecays_CMmom;
  std::vector<Double_t> GFdecays_CMmom_x;
  std::vector<Double_t> GFdecays_CMmom_y;
  std::vector<Double_t> GFdecays_CMmom_z;
  std::vector<Double_t> GFdecays_momloss;
  std::vector<Double_t> GFdecays_eloss;

  std::vector<Int_t> decays_id;
  std::vector<Double_t> decays_mom;
  std::vector<Double_t> decays_mom_x;
  std::vector<Double_t> decays_mom_y;
  std::vector<Double_t> decays_mom_z;
  std::vector<Double_t> decays_CMmom;
  std::vector<Double_t> decays_CMmom_x;
  std::vector<Double_t> decays_CMmom_y;
  std::vector<Double_t> decays_CMmom_z;

  Bool_t pipiflag;

  Int_t accident_multi;
  std::vector<Int_t> accident_id;

  Int_t xiresidual_multi;
  Int_t xipim_multi;
  Int_t xipip_multi;
  Int_t xiem_multi;
  Int_t xiep_multi;
  Int_t xip_multi;
  Int_t xippip_multi;

  std::vector<Int_t> xiresidual_id;
  std::vector<Double_t> xiresidual_dist2tgt;
  std::vector<Double_t> xiresidual_KFdist2prodvtx;
  std::vector<Double_t> xiresidual_GFdist2prodvtx;
  std::vector<Double_t> xiresidual_mass2;
  std::vector<Double_t> xiresidual_invbeta;
  std::vector<Double_t> xiresidual_mom;
  std::vector<Double_t> xiresidual_mom_x;
  std::vector<Double_t> xiresidual_mom_y;
  std::vector<Double_t> xiresidual_mom_z;
  std::vector<Double_t> xiresidual_charge;

  Int_t llresidual_multi;
  Int_t llpim_multi;
  Int_t llpip_multi;
  Int_t llem_multi;
  Int_t llep_multi;
  Int_t llp_multi;
  Int_t llppip_multi;

  std::vector<Int_t> llresidual_id;
  std::vector<Double_t> llresidual_dist2tgt;
  std::vector<Double_t> llresidual_GFdist2prodvtx;
  std::vector<Double_t> llresidual_KFdist2prodvtx;
  std::vector<Double_t> llresidual_mass2;
  std::vector<Double_t> llresidual_invbeta;
  std::vector<Double_t> llresidual_mom;
  std::vector<Double_t> llresidual_mom_x;
  std::vector<Double_t> llresidual_mom_y;
  std::vector<Double_t> llresidual_mom_z;
  std::vector<Double_t> llresidual_charge;

  Int_t residual_multi;
  Int_t pim_multi;
  Int_t pip_multi;
  Int_t em_multi;
  Int_t ep_multi;
  Int_t p_multi;
  Int_t ppip_multi;
  std::vector<Int_t> residual_id;
  std::vector<Double_t> residual_dist2tgt;
  std::vector<Double_t> residual_GFdist2prodvtx;
  std::vector<Double_t> residual_KFdist2prodvtx;
  std::vector<Double_t> residual_mass2;
  std::vector<Double_t> residual_invbeta;
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

  std::vector<Double_t> KFlldecays_mom;
  std::vector<Double_t> KFlldecays_mom_x;
  std::vector<Double_t> KFlldecays_mom_y;
  std::vector<Double_t> KFlldecays_mom_z;
  std::vector<Double_t> KFlldecays_CMmom;
  std::vector<Double_t> KFlldecays_CMmom_x;
  std::vector<Double_t> KFlldecays_CMmom_y;
  std::vector<Double_t> KFlldecays_CMmom_z;

  std::vector<Double_t> KFxidecays_mom;
  std::vector<Double_t> KFxidecays_mom_x;
  std::vector<Double_t> KFxidecays_mom_y;
  std::vector<Double_t> KFxidecays_mom_z;
  std::vector<Double_t> KFxidecays_CMmom;
  std::vector<Double_t> KFxidecays_CMmom_x;
  std::vector<Double_t> KFxidecays_CMmom_y;
  std::vector<Double_t> KFxidecays_CMmom_z;

  std::vector<Double_t> KFdecays_mom;
  std::vector<Double_t> KFdecays_mom_x;
  std::vector<Double_t> KFdecays_mom_y;
  std::vector<Double_t> KFdecays_mom_z;
  std::vector<Double_t> KFdecays_CMmom;
  std::vector<Double_t> KFdecays_CMmom_x;
  std::vector<Double_t> KFdecays_CMmom_y;
  std::vector<Double_t> KFdecays_CMmom_z;

  //For gamma searching
  Int_t llg_multi;
  std::vector<Int_t> llepidgamma;
  std::vector<Int_t> llemidgamma;
  std::vector<Double_t> llepmomgamma;
  std::vector<Double_t> llepmomgamma_x;
  std::vector<Double_t> llepmomgamma_y;
  std::vector<Double_t> llepmomgamma_z;
  std::vector<Double_t> llemmomgamma;
  std::vector<Double_t> llemmomgamma_x;
  std::vector<Double_t> llemmomgamma_y;
  std::vector<Double_t> llemmomgamma_z;
  std::vector<Double_t> llmomgamma;
  std::vector<Double_t> llmomgamma_x;
  std::vector<Double_t> llmomgamma_y;
  std::vector<Double_t> llmomgamma_z;
  std::vector<Double_t> llepidistgamma;
  std::vector<Double_t> llvtxgamma_x;
  std::vector<Double_t> llvtxgamma_y;
  std::vector<Double_t> llvtxgamma_z;

  Int_t xig_multi;
  std::vector<Int_t> xiepidgamma;
  std::vector<Int_t> xiemidgamma;
  std::vector<Double_t> xiepmomgamma;
  std::vector<Double_t> xiepmomgamma_x;
  std::vector<Double_t> xiepmomgamma_y;
  std::vector<Double_t> xiepmomgamma_z;
  std::vector<Double_t> xiemmomgamma;
  std::vector<Double_t> xiemmomgamma_x;
  std::vector<Double_t> xiemmomgamma_y;
  std::vector<Double_t> xiemmomgamma_z;
  std::vector<Double_t> ximomgamma;
  std::vector<Double_t> ximomgamma_x;
  std::vector<Double_t> ximomgamma_y;
  std::vector<Double_t> ximomgamma_z;
  std::vector<Double_t> xiepidistgamma;
  std::vector<Double_t> xivtxgamma_x;
  std::vector<Double_t> xivtxgamma_y;
  std::vector<Double_t> xivtxgamma_z;

  Int_t g_multi;
  std::vector<Int_t> epidgamma;
  std::vector<Int_t> emidgamma;
  std::vector<Double_t> epmomgamma;
  std::vector<Double_t> epmomgamma_x;
  std::vector<Double_t> epmomgamma_y;
  std::vector<Double_t> epmomgamma_z;
  std::vector<Double_t> emmomgamma;
  std::vector<Double_t> emmomgamma_x;
  std::vector<Double_t> emmomgamma_y;
  std::vector<Double_t> emmomgamma_z;
  std::vector<Double_t> momgamma;
  std::vector<Double_t> momgamma_x;
  std::vector<Double_t> momgamma_y;
  std::vector<Double_t> momgamma_z;
  std::vector<Double_t> epidistgamma;
  std::vector<Double_t> vtxgamma_x;
  std::vector<Double_t> vtxgamma_y;
  std::vector<Double_t> vtxgamma_z;

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
    m2Kurama.clear();
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
    cluster_x_center.clear();
    cluster_y_center.clear();
    cluster_z_center.clear();
    cluster_row_center.clear();
    cluster_houghflag.clear();

    remain_nclTpc = 0;
    remain_cluster_x.clear();
    remain_cluster_y.clear();
    remain_cluster_z.clear();
    remain_cluster_de.clear();
    remain_cluster_size.clear();
    remain_cluster_layer.clear();
    remain_cluster_mrow.clear();
    remain_cluster_de_center.clear();
    remain_cluster_x_center.clear();
    remain_cluster_y_center.clear();
    remain_cluster_z_center.clear();
    remain_cluster_row_center.clear();
    remain_cluster_houghflag.clear();

    ntTpc = 0;
    nhtrack.clear();
    trackid.clear();
    isXi.clear();
    isBeam.clear();
    isKurama.clear();
    isK18.clear();
    isAccidental.clear();
    isMultiloop.clear();
    isInTarget.clear();
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
    track_cluster_de_center.clear();
    track_cluster_x_center.clear();
    track_cluster_y_center.clear();
    track_cluster_z_center.clear();
    track_cluster_row_center.clear();

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
    reconmassPipair.clear();
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

    GFcharge.clear();
    GFchisqr.clear();
    GFtof.clear();
    GFtracklen.clear();
    GFpval.clear();
    GFfitstatus.clear();
    GFpdgcode.clear();
    GFnhtrack.clear();
    GFlayer.clear();
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

    GFinside.clear();

    GFntTpc_inside = 0;
    GFprodvtx_x = qnan;
    GFprodvtx_y = qnan;
    GFprodvtx_z = qnan;

    GFtracklen.clear();
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

    GFdecays_pdgcode.clear();
    GFdecays_nhtrack.clear();
    GFdecays_charge.clear();
    GFdecays_chisqr.clear();
    GFdecays_pval.clear();
    GFdecays_htofid.clear();
    GFdecays_tracklen.clear();
    GFdecays_tof.clear();
    GFdecays_mass2.clear();
    GFdecays_invbeta.clear();
    GFdecays_mom.clear();
    GFdecays_mom_x.clear();
    GFdecays_mom_y.clear();
    GFdecays_mom_z.clear();
    GFdecays_CMmom.clear();
    GFdecays_CMmom_x.clear();
    GFdecays_CMmom_y.clear();
    GFdecays_CMmom_z.clear();
    GFdecays_momloss.clear();
    GFdecays_eloss.clear();

    decays_id.clear();
    decays_mom.clear();
    decays_mom_x.clear();
    decays_mom_y.clear();
    decays_mom_z.clear();
    decays_CMmom.clear();
    decays_CMmom_x.clear();
    decays_CMmom_y.clear();
    decays_CMmom_z.clear();

    GFxidecays_pdgcode.clear();
    GFxidecays_nhtrack.clear();
    GFxidecays_charge.clear();
    GFxidecays_chisqr.clear();
    GFxidecays_tracktof.clear();
    GFxidecays_pval.clear();
    GFxidecays_htofid.clear();
    GFxidecays_tracklen.clear();
    GFxidecays_tof.clear();
    GFxidecays_mass2.clear();
    GFxidecays_invbeta.clear();
    GFxidecays_mom.clear();
    GFxidecays_mom_x.clear();
    GFxidecays_mom_y.clear();
    GFxidecays_mom_z.clear();
    GFxidecays_CMmom.clear();
    GFxidecays_CMmom_x.clear();
    GFxidecays_CMmom_y.clear();
    GFxidecays_CMmom_z.clear();
    GFxidecays_momloss.clear();
    GFxidecays_eloss.clear();

    xidecays_id.clear();
    xidecays_mom.clear();
    xidecays_mom_x.clear();
    xidecays_mom_y.clear();
    xidecays_mom_z.clear();
    xidecays_CMmom.clear();
    xidecays_CMmom_x.clear();
    xidecays_CMmom_y.clear();
    xidecays_CMmom_z.clear();

    GFlldecays_pdgcode.clear();
    GFlldecays_nhtrack.clear();
    GFlldecays_charge.clear();
    GFlldecays_chisqr.clear();
    GFlldecays_pval.clear();
    GFlldecays_htofid.clear();
    GFlldecays_tracklen.clear();
    GFlldecays_tof.clear();
    GFlldecays_mass2.clear();
    GFlldecays_invbeta.clear();
    GFlldecays_mom.clear();
    GFlldecays_mom_x.clear();
    GFlldecays_mom_y.clear();
    GFlldecays_mom_z.clear();
    GFlldecays_CMmom.clear();
    GFlldecays_CMmom_x.clear();
    GFlldecays_CMmom_y.clear();
    GFlldecays_CMmom_z.clear();
    GFlldecays_momloss.clear();
    GFlldecays_eloss.clear();

    lldecays_id.clear();
    lldecays_mom.clear();
    lldecays_mom_x.clear();
    lldecays_mom_y.clear();
    lldecays_mom_z.clear();
    lldecays_CMmom.clear();
    lldecays_CMmom_x.clear();
    lldecays_CMmom_y.clear();
    lldecays_CMmom_z.clear();

    xiflag = false;
    xipflag = false;
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

    lphiflag = false;
    phimass  = qnan;
    phidecayvtx_x  = qnan;
    phidecayvtx_y = qnan;
    phidecayvtx_z = qnan;
    phicosKK = qnan;
    phimom = qnan;
    phimom_x = qnan;
    phimom_y = qnan;
    phimom_z = qnan;
    kk_dist = qnan;
    GFphi_km_mass2 = qnan;
    GFphi_km_invbeta = qnan;
    GFphi_kp_mass2 = qnan;
    GFphi_kp_invbeta = qnan;
    GFphimass = qnan;
    GFphidecayvtx_x = qnan;
    GFphidecayvtx_y = qnan;
    GFphidecayvtx_z = qnan;
    GFphicosKK = qnan;
    GFphimom = qnan;
    GFphimom_x = qnan;
    GFphimom_y = qnan;
    GFphimom_z = qnan;
    GFkk_dist = qnan;
    GFphiprodvtx_dist = qnan;

    GFphimass_wKurama = qnan;
    GFphimom_wKurama = qnan;
    GFphimom_x_wKurama = qnan;
    GFphimom_y_wKurama = qnan;
    GFphimom_z_wKurama = qnan;

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

    GFlmass_alter1 = qnan;
    GFldecayvtx_x_alter1 = qnan;
    GFldecayvtx_y_alter1 = qnan;
    GFldecayvtx_z_alter1 = qnan;
    GFlmom_alter1 = qnan;
    GFlmom_x_alter1 = qnan;
    GFlmom_y_alter1 = qnan;
    GFlmom_z_alter1 = qnan;
    GFppi_dist_alter1 = qnan;
    GFltarget_dist_alter1 = qnan;
    GFltargetvtx_x_alter1 = qnan;
    GFltargetvtx_y_alter1 = qnan;
    GFltargetvtx_z_alter1 = qnan;
    GFltargetcenter_dist_alter1 = qnan;
    GFltargetcenter_x_alter1 = qnan;
    GFltargetcenter_y_alter1 = qnan;
    GFltargetcenter_z_alter1 = qnan;

    GFlmass_alter2 = qnan;
    GFldecayvtx_x_alter2 = qnan;
    GFldecayvtx_y_alter2 = qnan;
    GFldecayvtx_z_alter2 = qnan;
    GFlmom_alter2 = qnan;
    GFlmom_x_alter2 = qnan;
    GFlmom_y_alter2 = qnan;
    GFlmom_z_alter2 = qnan;
    GFppi_dist_alter2 = qnan;
    GFltarget_dist_alter2 = qnan;
    GFltargetvtx_x_alter2 = qnan;
    GFltargetvtx_y_alter2 = qnan;
    GFltargetvtx_z_alter2 = qnan;
    GFltargetcenter_dist_alter2 = qnan;
    GFltargetcenter_x_alter2 = qnan;
    GFltargetcenter_y_alter2 = qnan;
    GFltargetcenter_z_alter2 = qnan;

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

    KFlldecays_mom.clear();
    KFlldecays_mom_x.clear();
    KFlldecays_mom_y.clear();
    KFlldecays_mom_z.clear();
    KFlldecays_CMmom.clear();
    KFlldecays_CMmom_x.clear();
    KFlldecays_CMmom_y.clear();
    KFlldecays_CMmom_z.clear();

    KFxidecays_mom.clear();
    KFxidecays_mom_x.clear();
    KFxidecays_mom_y.clear();
    KFxidecays_mom_z.clear();
    KFxidecays_CMmom.clear();
    KFxidecays_CMmom_x.clear();
    KFxidecays_CMmom_y.clear();
    KFxidecays_CMmom_z.clear();

    KFdecays_mom.clear();
    KFdecays_mom_x.clear();
    KFdecays_mom_y.clear();
    KFdecays_mom_z.clear();
    KFdecays_CMmom.clear();
    KFdecays_CMmom_x.clear();
    KFdecays_CMmom_y.clear();
    KFdecays_CMmom_z.clear();

    accident_multi = 0;
    accident_id.clear();

    residual_multi = 0;
    pim_multi = 0;
    pip_multi = 0;
    em_multi = 0;
    ep_multi = 0;
    p_multi = 0;
    ppip_multi = 0;

    residual_id.clear();
    residual_dist2tgt.clear();
    residual_GFdist2prodvtx.clear();
    residual_KFdist2prodvtx.clear();
    residual_mass2.clear();
    residual_invbeta.clear();
    residual_mom.clear();
    residual_mom_x.clear();
    residual_mom_y.clear();
    residual_mom_z.clear();
    residual_charge.clear();

    xiresidual_multi = 0;
    xipim_multi = 0;
    xipip_multi = 0;
    xiem_multi = 0;
    xiep_multi = 0;
    xip_multi = 0;
    xippip_multi = 0;

    xiresidual_id.clear();
    xiresidual_dist2tgt.clear();
    xiresidual_KFdist2prodvtx.clear();
    xiresidual_GFdist2prodvtx.clear();
    xiresidual_mass2.clear();
    xiresidual_invbeta.clear();
    xiresidual_mom.clear();
    xiresidual_mom_x.clear();
    xiresidual_mom_y.clear();
    xiresidual_mom_z.clear();
    xiresidual_charge.clear();

    llresidual_multi = 0;
    llpim_multi = 0;
    llpip_multi = 0;
    llem_multi = 0;
    llep_multi = 0;
    llp_multi = 0;
    llppip_multi = 0;

    llresidual_id.clear();
    llresidual_dist2tgt.clear();
    llresidual_GFdist2prodvtx.clear();
    llresidual_KFdist2prodvtx.clear();
    llresidual_mass2.clear();
    llresidual_invbeta.clear();
    llresidual_mom.clear();
    llresidual_mom_x.clear();
    llresidual_mom_y.clear();
    llresidual_mom_z.clear();
    llresidual_charge.clear();

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

    llg_multi = 0;
    llepidgamma.clear();
    llemidgamma.clear();
    llepmomgamma.clear();
    llepmomgamma_x.clear();
    llepmomgamma_y.clear();
    llepmomgamma_z.clear();
    llemmomgamma.clear();
    llemmomgamma_x.clear();
    llemmomgamma_y.clear();
    llemmomgamma_z.clear();
    llmomgamma.clear();
    llmomgamma_x.clear();
    llmomgamma_y.clear();
    llmomgamma_z.clear();
    llepidistgamma.clear();
    llvtxgamma_x.clear();
    llvtxgamma_y.clear();
    llvtxgamma_z.clear();

    xig_multi = 0;
    xiepidgamma.clear();
    xiemidgamma.clear();
    xiepmomgamma.clear();
    xiepmomgamma_x.clear();
    xiepmomgamma_y.clear();
    xiepmomgamma_z.clear();
    xiemmomgamma.clear();
    xiemmomgamma_x.clear();
    xiemmomgamma_y.clear();
    xiemmomgamma_z.clear();
    ximomgamma.clear();
    ximomgamma_x.clear();
    ximomgamma_y.clear();
    ximomgamma_z.clear();
    xiepidistgamma.clear();
    xivtxgamma_x.clear();
    xivtxgamma_y.clear();
    xivtxgamma_z.clear();

    g_multi = 0;
    epidgamma.clear();
    emidgamma.clear();
    epmomgamma.clear();
    epmomgamma_x.clear();
    epmomgamma_y.clear();
    epmomgamma_z.clear();
    emmomgamma.clear();
    emmomgamma_x.clear();
    emmomgamma_y.clear();
    emmomgamma_z.clear();
    momgamma.clear();
    momgamma_x.clear();
    momgamma_y.clear();
    momgamma_z.clear();
    epidistgamma.clear();
    vtxgamma_x.clear();
    vtxgamma_y.clear();
    vtxgamma_z.clear();
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
  TTreeReaderValue<std::vector<Double_t>>* m2Org;
  TTreeReaderValue<std::vector<Double_t>>* thetaKurama;
  TTreeReaderValue<std::vector<Double_t>>* xtgtKurama;
  TTreeReaderValue<std::vector<Double_t>>* ytgtKurama;
  TTreeReaderValue<std::vector<Double_t>>* utgtKurama;
  TTreeReaderValue<std::vector<Double_t>>* vtgtKurama;
  TTreeReaderValue<std::vector<Double_t>>* pathwcKurama;
  TTreeReaderValue<std::vector<Double_t>>* xin;
  TTreeReaderValue<std::vector<Double_t>>* yin;
  TTreeReaderValue<std::vector<Double_t>>* zin;
  TTreeReaderValue<std::vector<Double_t>>* pxin;
  TTreeReaderValue<std::vector<Double_t>>* pyin;
  TTreeReaderValue<std::vector<Double_t>>* pzin;
  TTreeReaderValue<std::vector<Double_t>>* xout;
  TTreeReaderValue<std::vector<Double_t>>* yout;
  TTreeReaderValue<std::vector<Double_t>>* zout;
  TTreeReaderValue<std::vector<Double_t>>* pxout;
  TTreeReaderValue<std::vector<Double_t>>* pyout;
  TTreeReaderValue<std::vector<Double_t>>* pzout;

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
  TTreeReaderValue<std::vector<Double_t>>* pOrg;
  TTreeReaderValue<std::vector<Double_t>>* pCalc;
  TTreeReaderValue<std::vector<Double_t>>* pCorr;
  TTreeReaderValue<std::vector<Double_t>>* pCorrDE;
  TTreeReaderValue<std::vector<Double_t>>* xb;
  TTreeReaderValue<std::vector<Double_t>>* yb;
  TTreeReaderValue<std::vector<Double_t>>* ub;
  TTreeReaderValue<std::vector<Double_t>>* vb;
  TTreeReaderValue<std::vector<Double_t>>* xs;
  TTreeReaderValue<std::vector<Double_t>>* ys;
  TTreeReaderValue<std::vector<Double_t>>* us;
  TTreeReaderValue<std::vector<Double_t>>* vs;
  TTreeReaderValue<std::vector<Int_t>>* Kflag;

  TTreeReaderValue<Int_t>* ntTPCK18; // Number of Tracks
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

  //TTreeReaderValue<Int_t>* ntTPCKurama; // Number of Tracks
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
  TTreeReaderValue<std::vector<Double_t>>* xbTPC;
  TTreeReaderValue<std::vector<Double_t>>* ybTPC;
  TTreeReaderValue<std::vector<Double_t>>* ubTPC;
  TTreeReaderValue<std::vector<Double_t>>* vbTPC;
  TTreeReaderValue<std::vector<Double_t>>* xsTPC;
  TTreeReaderValue<std::vector<Double_t>>* ysTPC;
  TTreeReaderValue<std::vector<Double_t>>* usTPC;
  TTreeReaderValue<std::vector<Double_t>>* vsTPC;

  TTreeReaderValue<Int_t>* nclTpc; // Number of clusters
  TTreeReaderValue<Int_t>* remain_nclTpc; // Number of clusters without tracks
  TTreeReaderValue<std::vector<Double_t>>* cluster_x;
  TTreeReaderValue<std::vector<Double_t>>* cluster_y;
  TTreeReaderValue<std::vector<Double_t>>* cluster_z;
  TTreeReaderValue<std::vector<Double_t>>* cluster_de;
  TTreeReaderValue<std::vector<Int_t>>* cluster_size;
  TTreeReaderValue<std::vector<Int_t>>* cluster_layer;
  TTreeReaderValue<std::vector<Double_t>>* cluster_mrow;
  TTreeReaderValue<std::vector<Double_t>>* cluster_de_center;
  TTreeReaderValue<std::vector<Double_t>>* cluster_x_center;
  TTreeReaderValue<std::vector<Double_t>>* cluster_y_center;
  TTreeReaderValue<std::vector<Double_t>>* cluster_z_center;
  TTreeReaderValue<std::vector<Int_t>>* cluster_row_center;
  TTreeReaderValue<std::vector<Int_t>>* cluster_houghflag;

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
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* mom_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* mom_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* mom_z;
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
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_size;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_mrow;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_x_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_y_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_z_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_row_center;

  TTreeReaderValue<std::vector<Int_t>>* chargeIndistinguishable;
  TTreeReaderValue<std::vector<Int_t>>* pid_inverted;
  TTreeReaderValue<std::vector<Double_t>>* chisqr_inverted;
  TTreeReaderValue<std::vector<Double_t>>* pval_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_cx_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_cy_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_z0_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_r_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_dz_inverted;
  TTreeReaderValue<std::vector<Double_t>>* mom0_inverted;//Helix momentum at Y = 0

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

  TTreeReaderValue<std::vector<Int_t>>* isLambda;
  TTreeReaderValue<std::vector<Int_t>>* ncombiLambda;
  TTreeReaderValue<std::vector<Double_t>>* distLambda;
  TTreeReaderValue<std::vector<Double_t>>* angleLambda;
  TTreeReaderValue<std::vector<Double_t>>* bestmassLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* massLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxLambda_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxLambda_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxLambda_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* momLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* momLambda_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* momLambda_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* momLambda_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* decaysidLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* decaysmomLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* decaysmomLambda_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* decaysmomLambda_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* decaysmomLambda_z;

  TTreeReaderValue<Int_t>* nvtxTpcClustered;
  TTreeReaderValue<std::vector<Double_t>>* clusteredVtx_x;
  TTreeReaderValue<std::vector<Double_t>>* clusteredVtx_y;
  TTreeReaderValue<std::vector<Double_t>>* clusteredVtx_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* clusteredVtxid;

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

  static const auto ElectronMass = pdg::ElectronMass();
  static const auto PionMass = pdg::PionMass();
  static const auto KaonMass = pdg::KaonMass();
  static const auto PhiMass = pdg::Mass(333);
  static const auto ProtonMass = pdg::ProtonMass();
  static const auto LambdaMass = pdg::LambdaMass();
  static const auto XiMinusMass = pdg::XiMinusMass();
  static const auto m12C = 11.174864;
  static const auto m11B = 10.252548;
  static const auto m10Be = 9.325504;
  static const auto me = 0.001*0.5109989461;
  static const int XiMinusPdgCode = 3312;
  Double_t pdgmass[3] = {ProtonMass, KaonMass, PionMass};
  TVector3 tgtpos(0, 0, tpc::ZTarget);
  TVector3 qnan_vec = TVector3(qnan, qnan, qnan);

  //if( ievent%100000==0 ){
  //if( ievent%1000==0 ){
  //if( ievent%100==0 ){
  //if( ievent%10==0 ){
  if( ievent%1==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

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
  event.pOrg = **src.pOrg;
  event.pCalc = **src.pCalc;
  event.pCorr = **src.pCorr;
  event.pCorrDE = **src.pCorrDE;
  event.ub = **src.ub;
  event.vb = **src.vb;
  event.us = **src.us;
  event.vs = **src.vs;

  event.vtx = **src.vtx;
  event.vty = **src.vty;
  event.vtz = **src.vtz;
  event.closeDist = **src.closeDist;
  event.inside = **src.inside;

  event.pKurama = **src.pKurama;
  event.qKurama = **src.qKurama;
  event.chisqrKurama = **src.chisqrKurama;
  event.m2Kurama = **src.m2Org;
  event.thetaKurama = **src.thetaKurama;
  event.xtgtKurama = **src.xtgtKurama;
  event.ytgtKurama = **src.ytgtKurama;
  event.utgtKurama = **src.utgtKurama;
  event.vtgtKurama = **src.vtgtKurama;
  event.pathwcKurama = **src.pathwcKurama;
  event.xin = **src.xin;
  event.yin = **src.yin;
  event.zin = **src.zin;
  event.pxin = **src.pxin;
  event.pyin = **src.pyin;
  event.pzin = **src.pzin;
  event.xout = **src.xout;
  event.yout = **src.yout;
  event.zout = **src.zout;
  event.pxout = **src.pxout;
  event.pyout = **src.pyout;
  event.pzout = **src.pzout;

  event.ntK18 = **src.ntK18;
  event.pK18 = **src.pK18;
  event.chisqrK18 = **src.chisqrK18;
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
  event.thetaCMTPC = **src.thetaCMTPC;
  event.costCMTPC = **src.costCMTPC;
  event.xbTPC = **src.xbTPC;
  event.ybTPC = **src.ybTPC;
  event.ubTPC = **src.ubTPC;
  event.vbTPC = **src.vbTPC;
  event.xsTPC = **src.xsTPC;
  event.ysTPC = **src.ysTPC;
  event.usTPC = **src.usTPC;
  event.vsTPC = **src.vsTPC;

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

  if( event.nKK != 1 ) return false;
  if( event.isgoodTPCKurama.size()!=1 ) return false;
  if( event.isgoodTPCKurama[0]!=1 ) return false;
  if( event.insideTPC[0] != 1) return false;
#if kkevent
  if( event.kflagTPCKurama[0]!=1 ) return false;
#endif

  TVector3 kkvtxTPC(event.vtxTPC[0], event.vtyTPC[0], event.vtzTPC[0] + tpc::ZTarget);

  Double_t pKp = event.pCorr[0];
  Double_t uKp = event.us[0];
  Double_t vKp = event.vs[0];
  Double_t ptKp = pKp/std::sqrt(1.+uKp*uKp+vKp*vKp);
  TVector3 kp_mom(ptKp*uKp, ptKp*vKp, ptKp);

  Double_t pKm = event.pK18[0];
  Double_t uKm = event.ub[0];
  Double_t vKm = event.vb[0];
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
  event.BE[0] = 1000.*binding_energy; //MeV/c2
  Double_t binding_energy_LL = m10Be + 2.*LambdaMass - mm_12C; //GeV/c2
  event.BE_LL[0] = 1000.*binding_energy_LL; //MeV/c2

  event.km_mom_x[0] = km_mom.x();
  event.km_mom_y[0] = km_mom.y();
  event.km_mom_z[0] = km_mom.z();
  event.kp_mom_x[0] = kp_mom.x();
  event.kp_mom_y[0] = kp_mom.y();
  event.kp_mom_z[0] = kp_mom.z();

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

  event.nhHtof = **src.nhHtof;
  event.HtofSeg = **src.HtofSeg;
  event.tHtof = **src.tHtof;
  event.dtHtof = **src.dtHtof;
  event.deHtof = **src.deHtof;
  event.posHtof = **src.posHtof;
  Int_t ntTpc = **src.ntTpc;
  if( ntTpc == 0 ) return true;

  event.nclTpc = **src.nclTpc;
  event.remain_nclTpc = **src.remain_nclTpc;
#if SaveRawData
  event.cluster_x = **src.cluster_x;
  event.cluster_y = **src.cluster_y;
  event.cluster_z = **src.cluster_z;
  event.cluster_de = **src.cluster_de;
  event.cluster_size = **src.cluster_size;
  event.cluster_layer = **src.cluster_layer;
  event.cluster_mrow = **src.cluster_mrow;
  event.cluster_de_center = **src.cluster_de_center;
  event.cluster_x_center = **src.cluster_x_center;
  event.cluster_y_center = **src.cluster_y_center;
  event.cluster_z_center = **src.cluster_z_center;
  event.cluster_row_center = **src.cluster_row_center;
  event.cluster_houghflag = **src.cluster_houghflag;

  event.remain_cluster_x.resize(event.remain_nclTpc);
  event.remain_cluster_y.resize(event.remain_nclTpc);
  event.remain_cluster_z.resize(event.remain_nclTpc);
  event.remain_cluster_de.resize(event.remain_nclTpc);
  event.remain_cluster_size.resize(event.remain_nclTpc);
  event.remain_cluster_layer.resize(event.remain_nclTpc);
  event.remain_cluster_mrow.resize(event.remain_nclTpc);
  event.remain_cluster_de_center.resize(event.remain_nclTpc);
  event.remain_cluster_x_center.resize(event.remain_nclTpc);
  event.remain_cluster_y_center.resize(event.remain_nclTpc);
  event.remain_cluster_z_center.resize(event.remain_nclTpc);
  event.remain_cluster_row_center.resize(event.remain_nclTpc);
  event.remain_cluster_houghflag.resize(event.remain_nclTpc)
  Int_t icl_remain = 0;
  for( Int_t icl=0; icl<event.nclTpc; ++icl ){
    if(event.cluster_houghflag[icl]!=0) continue;
    event.remain_cluster_x[icl_remain] = event.cluster_x[icl];
    event.remain_cluster_y[icl_remain] = event.cluster_y[icl];
    event.remain_cluster_z[icl_remain] = event.cluster_z[icl];
    event.remain_cluster_de[icl_remain] = event.cluster_de[icl];
    event.remain_cluster_size[icl_remain] = event.cluster_size[icl];
    event.remain_cluster_layer[icl_remain] = event.cluster_layer[icl];
    event.remain_cluster_mrow[icl_remain] = event.cluster_mrow[icl];
    event.remain_cluster_de_center[icl_remain] = event.cluster_de_center[icl];
    event.remain_cluster_x_center[icl_remain] = event.cluster_x_center[icl];
    event.remain_cluster_y_center[icl_remain] = event.cluster_y_center[icl];
    event.remain_cluster_z_center[icl_remain] = event.cluster_z_center[icl];
    event.remain_cluster_row_center[icl_remain] = event.cluster_row_center[icl];
    event.remain_cluster_houghflag[icl_remain] = event.cluster_houghflag[icl];
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
  event.mom_x = **src.mom_x;
  event.mom_y = **src.mom_y;
  event.mom_z = **src.mom_z;
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
  event.track_cluster_size = **src.track_cluster_size;
  event.track_cluster_mrow = **src.track_cluster_mrow;
  event.track_cluster_de_center = **src.track_cluster_de_center;
  event.track_cluster_x_center = **src.track_cluster_x_center;
  event.track_cluster_y_center = **src.track_cluster_y_center;
  event.track_cluster_z_center = **src.track_cluster_z_center;
  event.track_cluster_row_center = **src.track_cluster_row_center;

  event.chargeIndistinguishable = **src.chargeIndistinguishable;
  event.pid_inverted = **src.pid_inverted;
  event.chisqr_inverted = **src.chisqr_inverted;
  event.pval_inverted = **src.pval_inverted;
  event.helix_cx_inverted = **src.helix_cx_inverted;
  event.helix_cy_inverted = **src.helix_cy_inverted;
  event.helix_z0_inverted = **src.helix_z0_inverted;
  event.helix_r_inverted = **src.helix_r_inverted;
  event.helix_dz_inverted = **src.helix_dz_inverted;
  event.mom0_inverted = **src.mom0_inverted;

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

  event.isLambda = **src.isLambda;
  event.ncombiLambda = **src.ncombiLambda;
  event.distLambda = **src.distLambda;
  event.angleLambda = **src.angleLambda;
  event.bestmassLambda = **src.bestmassLambda;
  event.massLambda = **src.massLambda;
  event.vtxLambda_x = **src.vtxLambda_x;
  event.vtxLambda_y = **src.vtxLambda_y;
  event.vtxLambda_z = **src.vtxLambda_z;
  event.momLambda = **src.momLambda;
  event.momLambda_x = **src.momLambda_x;
  event.momLambda_y = **src.momLambda_y;
  event.momLambda_z = **src.momLambda_z;
  event.decaysidLambda = **src.decaysidLambda;
  event.decaysmomLambda = **src.decaysmomLambda;
  event.decaysmomLambda_x = **src.decaysmomLambda_x;
  event.decaysmomLambda_y = **src.decaysmomLambda_y;
  event.decaysmomLambda_z = **src.decaysmomLambda_z;

  event.nvtxTpcClustered = **src.nvtxTpcClustered;
  event.clusteredVtx_x = **src.clusteredVtx_x;
  event.clusteredVtx_y = **src.clusteredVtx_y;
  event.clusteredVtx_z = **src.clusteredVtx_z;
  event.clusteredVtxid = **src.clusteredVtxid;

  Int_t id_kp = -1;
  for(Int_t k=0;k<ntTpc;k++){
    if(event.isKurama[k]==1) id_kp = k;
  }

  Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);

  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCTracks(**src.ntTpc, **src.isK18, **src.isKurama,
			 **src.charge, **src.nhtrack, **src.helix_cx,
			 **src.helix_cy, **src.helix_z0, **src.helix_r,
			 **src.helix_dz, **src.hitlayer, **src.track_cluster_mrow,
			 **src.helix_t, **src.track_cluster_de, **src.resolution_x,
			 **src.resolution_y, **src.resolution_z, **src.hitpos_x,
			 **src.hitpos_y, **src.hitpos_z);

  HypTPCTask& GFTrackCont = HypTPCTask::GetInstance();
  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix(it);
    if( !tp ) continue;
    event.pid[it] = tp -> GetPid();
    event.isElectron[it] = Kinematics::HypTPCdEdxElectron(event.dEdx[it], event.mom0[it]);
    event.nsigma_triton[it] = Kinematics::HypTPCdEdxNsigmaTriton(event.dEdx[it], event.mom0[it]);
    event.nsigma_deutron[it] = Kinematics::HypTPCdEdxNsigmaDeutron(event.dEdx[it], event.mom0[it]);
    event.nsigma_proton[it] = Kinematics::HypTPCdEdxNsigmaProton(event.dEdx[it], event.mom0[it]);
    event.nsigma_kaon[it] = Kinematics::HypTPCdEdxNsigmaKaon(event.dEdx[it], event.mom0[it]);
    event.nsigma_pion[it] = Kinematics::HypTPCdEdxNsigmaPion(event.dEdx[it], event.mom0[it]);
    event.nsigma_electron[it] = Kinematics::HypTPCdEdxNsigmaElectron(event.dEdx[it], event.mom0[it]);

    if(event.isKurama[it]==1) GFTrackCont.AddHelixTrack(321, tp);
    else if(event.isElectron[it]==1) GFTrackCont.AddHelixTrack(event.charge[it]*(-11), tp);
    else{
      std::vector<Int_t> pdgcode;
      Kinematics::HypTPCPID_PDGCode(event.charge[it], event.pid[it], pdgcode);
      if((event.pid[it]&4)!=4 && TMath::Abs(event.nsigma_electron[it]) < 3)
	pdgcode.push_back(event.charge[it]*(-11));
      GFTrackCont.AddHelixTrack(pdgcode, tp);
    }
  }
  GFTrackCont.FitTracks();
  HF1( 1, event.status++ );

  Int_t GFntTpc = GFTrackCont.GetNTrack();
  if(GFntTpc!=event.ntTpc){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    return true;
  }

  event.GFcharge.resize(GFntTpc);
  event.GFchisqr.resize(GFntTpc);
  event.GFtof.resize(GFntTpc);
  //event.GFtracklen.resize(GFntTpc);
  event.GFpval.resize(GFntTpc);
  event.GFfitstatus.resize(GFntTpc);
  event.GFpdgcode.resize(GFntTpc);
  event.GFnhtrack.resize(GFntTpc);
  event.GFlayer.resize(GFntTpc);
  event.GFpos_x.resize(GFntTpc);
  event.GFpos_y.resize(GFntTpc);
  event.GFpos_z.resize(GFntTpc);
  event.GFmom.resize(GFntTpc);
  event.GFmom_x.resize(GFntTpc);
  event.GFmom_y.resize(GFntTpc);
  event.GFmom_z.resize(GFntTpc);
  event.GFresidual_x.resize(GFntTpc);
  event.GFresidual_y.resize(GFntTpc);
  event.GFresidual_z.resize(GFntTpc);
  event.GFresidual_p.resize(GFntTpc);
  event.GFresidual_px.resize(GFntTpc);
  event.GFresidual_py.resize(GFntTpc);
  event.GFresidual_pz.resize(GFntTpc);
  event.GFinside.resize(GFntTpc);

  event.GFtracklen.resize(GFntTpc);
  event.GFtrack2vtxdist.resize(GFntTpc);
  event.GFcalctof.resize(GFntTpc);
  event.GFposx.resize(GFntTpc);
  event.GFposy.resize(GFntTpc);
  event.GFposz.resize(GFntTpc);
  event.GFsegHtof.resize(GFntTpc);
  event.GFtofHtof.resize(GFntTpc);
  event.GFtdiffHtof.resize(GFntTpc);
  event.GFposHtof.resize(GFntTpc);
  event.GFinvbeta.resize(GFntTpc);
  event.GFm2.resize(GFntTpc);
  event.nsigma_tritonHtof.resize(ntTpc);
  event.nsigma_deutronHtof.resize(ntTpc);
  event.nsigma_protonHtof.resize(ntTpc);
  event.nsigma_kaonHtof.resize(ntTpc);
  event.nsigma_pionHtof.resize(ntTpc);
  event.nsigma_electronHtof.resize(ntTpc);

  Int_t ntrack_intarget = 0;
  Double_t x0[100] = {0};
  Double_t y0[100] = {0};
  Double_t u0[100] = {0};
  Double_t v0[100] = {0};
  for( Int_t igf=0; igf<GFntTpc; ++igf ){
    event.GFfitstatus[igf] = (int)GFTrackCont.TrackCheck(igf);
    if(!GFTrackCont.TrackCheck(igf)) continue;
    int nh = GFTrackCont.GetNHits(igf);
    event.GFlayer[igf].resize(nh);
    event.GFpos_x[igf].resize(nh);
    event.GFpos_y[igf].resize(nh);
    event.GFpos_z[igf].resize(nh);
    event.GFmom[igf].resize(nh);
    event.GFmom_x[igf].resize(nh);
    event.GFmom_y[igf].resize(nh);
    event.GFmom_z[igf].resize(nh);
    event.GFresidual_x[igf].resize(nh);
    event.GFresidual_y[igf].resize(nh);
    event.GFresidual_z[igf].resize(nh);
    event.GFresidual_p[igf].resize(nh);
    event.GFresidual_px[igf].resize(nh);
    event.GFresidual_py[igf].resize(nh);
    event.GFresidual_pz[igf].resize(nh);

    event.GFchisqr[igf] = GFTrackCont.GetChi2NDF(igf);
    event.GFcharge[igf] = GFTrackCont.GetCharge(igf);
    event.GFtof[igf] = GFTrackCont.GetTrackTOF(igf, 0, -1);
    event.GFpval[igf] = GFTrackCont.GetPvalue(igf);
    event.GFnhtrack[igf] = GFTrackCont.GetNHits(igf);
    event.GFpdgcode[igf] = GFTrackCont.GetPDGcode(igf);
    for( Int_t ihit=0; ihit<nh; ++ihit ){
      TVector3 hit = GFTrackCont.GetPos(igf, ihit);
      TVector3 mom = GFTrackCont.GetMom(igf, ihit);
      Int_t layer = (int)event.hitlayer[igf][ihit];
      event.GFlayer[igf][ihit] = layer;
      event.GFmom_x[igf][ihit] = mom.x();
      event.GFmom_y[igf][ihit] = mom.y();
      event.GFmom_z[igf][ihit] = mom.z();
      event.GFmom[igf][ihit] = mom.Mag();
      event.GFpos_x[igf][ihit] = hit.x();
      event.GFpos_y[igf][ihit] = hit.y();
      event.GFpos_z[igf][ihit] = hit.z();

      event.GFresidual_x[igf][ihit] = hit.x() - event.hitpos_x[igf][ihit];
      event.GFresidual_y[igf][ihit] = hit.y() - event.hitpos_y[igf][ihit];
      event.GFresidual_z[igf][ihit] = hit.z() - event.hitpos_z[igf][ihit];

      double chargetest = event.GFcharge[igf]*event.charge[igf];
      event.GFresidual_p[igf][ihit] = mom.Mag() - event.mom0[igf];
      event.GFresidual_px[igf][ihit] = mom.x() - chargetest*event.mom_x[igf][ihit];
      event.GFresidual_py[igf][ihit] = mom.y() - chargetest*event.mom_y[igf][ihit];
      event.GFresidual_pz[igf][ihit] = mom.z() - chargetest*event.mom_z[igf][ihit];
    }//ih
    if(event.isBeam[igf]==1) continue;
    if(event.isK18[igf]==1) continue;
    if(event.isAccidental[igf]==1) continue;
    if(GFTrackCont.IsInsideTarget(igf)){
      event.GFinside[igf] = 1;
      TVector3 posv; TVector3 momv; double len; double tof;
      if(GFTrackCont.ExtrapolateToTargetCenter(igf, posv, momv, len, tof)){
	x0[ntrack_intarget] = posv.x();
	y0[ntrack_intarget] = posv.y();
	u0[ntrack_intarget] = momv.x()/momv.z();
	v0[ntrack_intarget] = momv.y()/momv.z();
	ntrack_intarget++;
      }
    }
    else event.GFinside[igf] = 0;

    TVector3 vertex = Kinematics::MultitrackVertex(ntrack_intarget, x0, y0, u0, v0);
    event.GFntTpc_inside = ntrack_intarget;
    event.GFprodvtx_x = vertex.x();
    event.GFprodvtx_y = vertex.y();
    event.GFprodvtx_z = vertex.z();
    {
      Int_t repid = -1;
      Int_t hitid_htof; Double_t tof; Double_t len;
      TVector3 pos_htof; Double_t track2tgt_dist;
      Bool_t htofextrapolation =
	GFTrackCont.TPCHTOFTrackMatching(igf, repid, vertex,
					 event.HtofSeg, event.posHtof,
					 hitid_htof, tof,
					 len, pos_htof, track2tgt_dist);
      if(htofextrapolation){
	event.GFtracklen[igf] = len;
	event.GFtrack2vtxdist[igf] = track2tgt_dist;
	event.GFcalctof[igf] = tof;
	event.GFposx[igf] = pos_htof.x();
	event.GFposy[igf] = pos_htof.y();
	event.GFposz[igf] = pos_htof.z();
	event.GFsegHtof[igf] = event.HtofSeg[hitid_htof];
	event.GFtofHtof[igf] = event.tHtof[hitid_htof];
	event.GFposHtof[igf] = event.posHtof[hitid_htof];

	Double_t beta = len/event.tHtof[hitid_htof]/MathTools::C();
	event.GFinvbeta[igf] = 1./beta;
	Double_t mass2 = Kinematics::MassSquare(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
	event.GFm2[igf] = mass2;

	event.nsigma_tritonHtof[igf] = Kinematics::HypTPCHTOFNsigmaTriton(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
	event.nsigma_deutronHtof[igf] = Kinematics::HypTPCHTOFNsigmaDeutron(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
	event.nsigma_protonHtof[igf] = Kinematics::HypTPCHTOFNsigmaProton(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
	event.nsigma_kaonHtof[igf] = Kinematics::HypTPCHTOFNsigmaKaon(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
	event.nsigma_pionHtof[igf] = Kinematics::HypTPCHTOFNsigmaPion(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
	event.nsigma_electronHtof[igf] = Kinematics::HypTPCHTOFNsigmaElectron(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
      }
    } //common
  }//igf

#if DebugDisp
  std::cout<<"0. Before particle reconstruction, checking each track is originated from the target or not."<<std::endl;
#endif

  std::vector<Double_t> target_kurama_mass2_container;
  std::vector<Double_t> target_kurama_invbeta_container;
  std::vector<TVector3> target_kurama_mom_container;

  std::vector<Int_t> target_accidental_id_container;
  std::vector<Int_t> target_km_id_container;
  std::vector<Double_t> target_km_mass2_container;
  std::vector<Double_t> target_km_invbeta_container;
  std::vector<TVector3> target_km_mom_container;
  std::vector<Int_t> target_km_repid_container;

  std::vector<Int_t> target_em_id_container, target_ep_id_container,
    target_empim_id_container, target_eppip_id_container;
  std::vector<Double_t> target_em_mass2_container, target_ep_mass2_container,
    target_empim_mass2_container, target_eppip_mass2_container;
  std::vector<Double_t> target_em_invbeta_container, target_ep_invbeta_container,
    target_empim_invbeta_container, target_eppip_invbeta_container;
  std::vector<TVector3> target_em_mom_container, target_ep_mom_container,
    target_empim_mom_container, target_eppip_mom_container;
  std::vector<Double_t> target_em_dist2tgt_container, target_ep_dist2tgt_container,
    target_empim_dist2tgt_container, target_eppip_dist2tgt_container;
  std::vector<Int_t> target_em_repid_container, target_ep_repid_container,
    target_empim_repid_container, target_eppip_repid_container;

  std::vector<Int_t> target_p_id_container, target_pip_id_container,
    target_pim_id_container, target_ppip_id_container;
  std::vector<Double_t> target_p_mass2_container, target_pip_mass2_container,
    target_pim_mass2_container, target_ppip_mass2_container;
  std::vector<Double_t> target_p_invbeta_container, target_pip_invbeta_container,
    target_pim_invbeta_container, target_ppip_invbeta_container;
  std::vector<TVector3> target_p_mom_container, target_pip_mom_container,
    target_pim_mom_container, target_ppip_mom_container;
  std::vector<Double_t> target_p_dist2tgt_container, target_pip_dist2tgt_container,
    target_pim_dist2tgt_container, target_ppip_dist2tgt_container;
  std::vector<Int_t> target_p_repid_container, target_pip_repid_container,
    target_pim_repid_container, target_ppip_repid_container;

  std::vector<Int_t> intarget_id_container;
  std::vector<Int_t> notintarget_id_container;
  event.isInTarget.resize(ntTpc, 0);
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
    if(event.isKurama[it]!=1) continue;
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( it );
    Int_t nhit = event.nhtrack[it];
    TVector3 start(event.hitpos_x[it][0], event.hitpos_y[it][0], event.hitpos_z[it][0]);
    TVector3 end(event.hitpos_x[it][nhit-1], event.hitpos_y[it][nhit-1], event.hitpos_z[it][nhit-1]);

    Double_t mass2 = qnan;
    Double_t inverse_beta = qnan;
    TVector3 mom(qnan, qnan, qnan);

    Int_t repid = -1;
    if(event.isBeam[it]==1 || !GFTrackCont.TrackCheck(it, repid)){
      target_kurama_mass2_container.push_back(mass2);
      target_kurama_invbeta_container.push_back(inverse_beta);
      target_kurama_mom_container.push_back(mom);
    }
    else{
      Int_t htofhitid_k; Double_t tracklen_k;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation_k =
	GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
					 event.HtofSeg, event.posHtof,
					 htofhitid_k, tof, tracklen_k,
					 pos, track2tgt_dist);
      mom = GFTrackCont.GetMom(it, -1, repid);
      if(htofextrapolation_k){
	mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_k, event.tHtof[htofhitid_k]);
	inverse_beta = MathTools::C()*event.tHtof[htofhitid_k]/tracklen_k;
      }
      target_kurama_mass2_container.push_back(mass2);
      target_kurama_invbeta_container.push_back(inverse_beta);
      target_kurama_mom_container.push_back(mom);
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
  }

  for(Int_t ith=0;ith<intarget_id_container.size();ith++){
    Int_t it = intarget_id_container[ith];
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( it );
    Int_t nhit = event.nhtrack[it];
    TVector3 start(event.hitpos_x[it][0], event.hitpos_y[it][0], event.hitpos_z[it][0]);
    TVector3 end(event.hitpos_x[it][nhit-1], event.hitpos_y[it][nhit-1], event.hitpos_z[it][nhit-1]);

    if((event.pid[it]&2)==2 && event.charge[it]==-1){ //k-
      Bool_t km_pid = Kinematics::HypTPCdEdxKaon(event.dEdx[it], event.mom0[it]);
      if(km_pid){
	Double_t mass2 = qnan;
	Double_t inverse_beta = qnan;
	TVector3 mom(qnan, qnan, qnan);

	Int_t repid = 0;
	if((1&event.pid[it])==1) repid += 1;
	target_km_id_container.push_back(it);
	target_km_repid_container.push_back(repid);
	if(!GFTrackCont.TrackCheck(it, repid)){
	  target_km_mass2_container.push_back(mass2);
	  target_km_invbeta_container.push_back(inverse_beta);
	  target_km_mom_container.push_back(mom);
	}
	else{
	  Int_t htofhitid_k; Double_t tracklen_k;
	  Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	  Bool_t htofextrapolation_k =
	    GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
					     event.HtofSeg, event.posHtof,
					     htofhitid_k, tof, tracklen_k,
					     pos, track2tgt_dist);
	  mom = GFTrackCont.GetMom(it, -1, repid);
	  if(htofextrapolation_k){
	    mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_k, event.tHtof[htofhitid_k]);
	    inverse_beta = MathTools::C()*event.tHtof[htofhitid_k]/tracklen_k;
	  }

	  target_km_mass2_container.push_back(mass2);
	  target_km_invbeta_container.push_back(inverse_beta);
	  target_km_mom_container.push_back(mom);
	}
      }
    }

    if(event.isElectron[it]==1){ //e+, e-
      Double_t mass2 = qnan;
      Double_t inverse_beta = qnan;
      TVector3 mom(qnan, qnan, qnan);
      Int_t repid = -1;
      if(event.charge[it]==1){
	target_ep_id_container.push_back(it);
	target_ep_dist2tgt_container.push_back(tp -> GetClosestDist());
	target_ep_repid_container.push_back(repid);
      }
      else{
	target_em_id_container.push_back(it);
	target_em_dist2tgt_container.push_back(tp -> GetClosestDist());
	target_em_repid_container.push_back(repid);
      }

      if(!GFTrackCont.TrackCheck(it, repid)){
	if(event.charge[it]==1){
	  target_ep_mass2_container.push_back(mass2);
	  target_ep_invbeta_container.push_back(inverse_beta);
	  target_ep_mom_container.push_back(mom);
	}
	else{
	  target_em_mass2_container.push_back(mass2);
	  target_em_invbeta_container.push_back(inverse_beta);
	  target_em_mom_container.push_back(mom);
	}
	continue;
      }

      Int_t htofhitid; Double_t tracklen;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation =
	GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
					 event.HtofSeg, event.posHtof,
					 htofhitid, tof, tracklen,
					 pos, track2tgt_dist);
      mom = GFTrackCont.GetMom(it, -1, repid);
      if(htofextrapolation){
	mass2 = Kinematics::MassSquare(mom.Mag(), tracklen, event.tHtof[htofhitid]);
	inverse_beta = MathTools::C()*event.tHtof[htofhitid]/tracklen;
      }
      if(event.charge[it]==1){
	target_ep_mass2_container.push_back(mass2);
	target_ep_invbeta_container.push_back(inverse_beta);
	target_ep_mom_container.push_back(mom);
      }
      else{
	target_em_mass2_container.push_back(mass2);
	target_em_invbeta_container.push_back(inverse_beta);
	target_em_mom_container.push_back(mom);
      }
    }
    else if((event.pid[it]&4)==4 && (event.pid[it]&1)!=1 && event.charge[it]==1){ //proton

      Double_t mass2 = qnan;
      Double_t inverse_beta = qnan;
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

      target_p_repid_container.push_back(repid);
      if(!GFTrackCont.TrackCheck(it, repid)){
	target_p_mass2_container.push_back(mass2);
	target_p_invbeta_container.push_back(inverse_beta);
	target_p_mom_container.push_back(mom);
	continue;
      }

      Int_t htofhitid_p; Double_t tracklen_p;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation_p =
	GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
					 event.HtofSeg, event.posHtof,
					 htofhitid_p, tof, tracklen_p,
					 pos, track2tgt_dist);
      mom = GFTrackCont.GetMom(it, -1, repid);
      if(htofextrapolation_p){
	mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_p, event.tHtof[htofhitid_p]);
	inverse_beta = MathTools::C()*event.tHtof[htofhitid_p]/tracklen_p;
      }
      target_p_mass2_container.push_back(mass2);
      target_p_invbeta_container.push_back(inverse_beta);
      target_p_mom_container.push_back(mom);
    }
    else if((event.pid[it]&1)==1 && event.charge[it]==-1){ //pi-

      Double_t slope = event.helix_dz[it];
      Double_t helixmom = event.mom0[it];
      //if(TMath::Abs(slope)<0.05 && TMath::Abs(helixmom)>0.5){
      if(TMath::Abs(slope)<0.1 && TMath::Abs(helixmom)>0.5 &&
	 (start.x()-end.x())>-10 && (start.x()-end.x())<50.){
	event.isAccidental[it] = 1;
	target_accidental_id_container.push_back(it); //Accidental beam on the target
	continue; //Accidental K-
      }

      Double_t mass2 = qnan;
      Double_t inverse_beta = qnan;
      TVector3 mom(qnan, qnan, qnan);
      target_pim_id_container.push_back(it);
      target_pim_dist2tgt_container.push_back(tp -> GetClosestDist());

      Int_t repid = 0;
      target_pim_repid_container.push_back(repid);
      if(!GFTrackCont.TrackCheck(it, repid)){
	target_pim_mass2_container.push_back(mass2);
	target_pim_invbeta_container.push_back(inverse_beta);
	target_pim_mom_container.push_back(mom);
	continue;
      }

      Int_t htofhitid_pi; Double_t tracklen_pi;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation_pi =
	GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
					 event.HtofSeg, event.posHtof,
					 htofhitid_pi, tof, tracklen_pi,
					 pos, track2tgt_dist);
      mom = GFTrackCont.GetMom(it, -1, repid);
      if(htofextrapolation_pi){
	mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_pi, event.tHtof[htofhitid_pi]);
	inverse_beta = MathTools::C()*event.tHtof[htofhitid_pi]/tracklen_pi;
      }
      target_pim_mass2_container.push_back(mass2);
      target_pim_invbeta_container.push_back(inverse_beta);
      target_pim_mom_container.push_back(mom);

      if(TMath::Abs(event.nsigma_electron[it]) < 3){ //e-
	Double_t mass2 = qnan;
	Double_t inverse_beta = qnan;
	TVector3 mom(qnan, qnan, qnan);
	target_empim_id_container.push_back(it);
	target_empim_dist2tgt_container.push_back(tp -> GetClosestDist());
	target_empim_repid_container.push_back(repid);

	if(!GFTrackCont.TrackCheck(it, repid)){
	  target_empim_mass2_container.push_back(mass2);
	  target_empim_invbeta_container.push_back(inverse_beta);
	  target_empim_mom_container.push_back(mom);
	  continue;
	}

	Int_t htofhitid; Double_t tracklen;
	Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	Bool_t htofextrapolation =
	  GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
					   event.HtofSeg, event.posHtof,
					   htofhitid, tof, tracklen,
					   pos, track2tgt_dist);
	mom = GFTrackCont.GetMom(it, -1, repid);
	if(htofextrapolation){
	  mass2 = Kinematics::MassSquare(mom.Mag(), tracklen, event.tHtof[htofhitid]);
	  inverse_beta = MathTools::C()*event.tHtof[htofhitid]/tracklen;
	}
	target_empim_mass2_container.push_back(mass2);
	target_empim_invbeta_container.push_back(inverse_beta);
	target_empim_mom_container.push_back(mom);
      }
    }
    else if((event.pid[it]&4)!=4 && (event.pid[it]&1)==1 && event.charge[it]==1){ //pi+

      Double_t mass2 = qnan;
      Double_t inverse_beta = qnan;
      TVector3 mom(qnan, qnan, qnan);
      target_pip_id_container.push_back(it);
      target_pip_dist2tgt_container.push_back(tp -> GetClosestDist());

      Int_t repid = 0;
      target_pip_repid_container.push_back(repid);
      if(!GFTrackCont.TrackCheck(it, repid)){
	target_pip_mass2_container.push_back(mass2);
	target_pip_invbeta_container.push_back(inverse_beta);
	target_pip_mom_container.push_back(mom);
	continue;
      }

      Int_t htofhitid_pi; Double_t tracklen_pi;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation_pi =
	GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
					 event.HtofSeg, event.posHtof,
					 htofhitid_pi, tof, tracklen_pi,
					 pos, track2tgt_dist);
      mom = GFTrackCont.GetMom(it, -1, repid);
      if(htofextrapolation_pi){
	mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_pi, event.tHtof[htofhitid_pi]);
	inverse_beta = MathTools::C()*event.tHtof[htofhitid_pi]/tracklen_pi;
      }
      target_pip_mass2_container.push_back(mass2);
      target_pip_invbeta_container.push_back(inverse_beta);
      target_pip_mom_container.push_back(mom);

      if(TMath::Abs(event.nsigma_electron[it]) < 3){ //e+
      	Double_t mass2 = qnan;
	Double_t inverse_beta = qnan;
	TVector3 mom(qnan, qnan, qnan);
	target_eppip_id_container.push_back(it);
	target_eppip_dist2tgt_container.push_back(tp -> GetClosestDist());
	target_eppip_repid_container.push_back(repid);

	if(!GFTrackCont.TrackCheck(it, repid)){
	  target_eppip_mass2_container.push_back(mass2);
	  target_eppip_invbeta_container.push_back(inverse_beta);
	  target_eppip_mom_container.push_back(mom);
	  continue;
	}

	Int_t htofhitid; Double_t tracklen;
	Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	Bool_t htofextrapolation =
	  GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
					   event.HtofSeg, event.posHtof,
					   htofhitid, tof, tracklen,
					   pos, track2tgt_dist);
	mom = GFTrackCont.GetMom(it, -1, repid);
	if(htofextrapolation){
	  mass2 = Kinematics::MassSquare(mom.Mag(), tracklen, event.tHtof[htofhitid]);
	  inverse_beta = MathTools::C()*event.tHtof[htofhitid]/tracklen;
	}
	target_eppip_mass2_container.push_back(mass2);
	target_eppip_invbeta_container.push_back(inverse_beta);
	target_eppip_mom_container.push_back(mom);
      }
    }
    else if(((event.pid[it]&4)==4 || (event.pid[it]&1)==1) && event.charge[it]==1){ //p or pi+ with high-mom
      Double_t mass2 = qnan;
      Double_t inverse_beta = qnan;
      TVector3 mom(qnan, qnan, qnan);
      target_ppip_id_container.push_back(it);
      target_ppip_dist2tgt_container.push_back(tp -> GetClosestDist());

      Int_t repid = -1;
      target_ppip_repid_container.push_back(repid);
      if(!GFTrackCont.TrackCheck(it, repid)){
	target_ppip_mass2_container.push_back(mass2);
	target_ppip_invbeta_container.push_back(inverse_beta);
	target_ppip_mom_container.push_back(mom);
	continue;
      }

      Int_t htofhitid_ppi; Double_t tracklen_ppi;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation =
	GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
					 event.HtofSeg, event.posHtof,
					 htofhitid_ppi, tof, tracklen_ppi,
					 pos, track2tgt_dist);
      mom = GFTrackCont.GetMom(it, -1, repid);
      if(htofextrapolation){
	mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_ppi, event.tHtof[htofhitid_ppi]);
	inverse_beta = MathTools::C()*event.tHtof[htofhitid_ppi]/tracklen_ppi;
      }
      target_ppip_mass2_container.push_back(mass2);
      target_ppip_invbeta_container.push_back(inverse_beta);
      target_ppip_mom_container.push_back(mom);
    }
  }

  //Veto accidental tracks near the clustered vertex of accidental dummy tracks
  for(Int_t ivtx=0;ivtx<event.nvtxTpcClustered;ivtx++){
    TVector3 vtx = TVector3(event.clusteredVtx_x[ivtx], event.clusteredVtx_y[ivtx], event.clusteredVtx_z[ivtx]+tpc::ZTarget);
    if(TMath::Abs(vtx.x()) < 15 && TMath::Abs(vtx.y()) < 35 && TMath::Abs(vtx.z()) < 15) continue; //LL, Xi- decays can make cluster within this region
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

  //for L, LL
  std::vector<Int_t> L_p_id_container, L_pi_id_container;
  std::vector<Int_t> L_p_repid_container, L_pi_repid_container;
  std::vector<TVector3> L_p_mom_container, L_pi_mom_container;
  std::vector<Double_t> L_mass_container;
  std::vector<Double_t> L_ppidist_container;
  std::vector<Double_t> L_targetdist_container;
  std::vector<TVector3> L_mom_container, L_vtx_container, L_targetvtx_container;

  //Xi/L candidates searching
#if DebugDisp
  std::cout<<"1. Xi/L candidate searching"<<std::endl;
#endif
  Int_t l_candidates = 0;
  Int_t xi_candidates = 0;

  std::vector<Int_t> xi_l_container;
  std::vector<Int_t> xi_p_container, xi_pi_container, xi_pi2_container;
  std::vector<Int_t> p_repid_container, pi_repid_container, pi2_repid_container;
  std::vector<TVector3> xi_mom_container, xi_decayvertex_container;
  std::vector<TVector3> l_mom_container, l_vert_container;
  std::vector<Double_t> xi_mass_container, lambda_mass_container;
  std::vector<TVector3> xi_p_mom_container, xi_pi_mom_container, xi_pi2_mom_container;
  std::vector<Double_t> ppi_closedist; std::vector<Double_t> lpi_closedist;
  std::vector<Double_t> xi_targetdist_container;
  std::vector<TVector3> xi_targetvtx_container, xi_targetmom_container;
  {
    for(Int_t it1=0;it1<ntTpc;it1++){
      if(event.isElectron[it1]==1) continue;
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
	if(event.isElectron[it2]==1) continue;
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

	Double_t ltarget_dist;
	TVector3 ltarget_vtx =
	  Kinematics::CalcCloseDistLambda(tgtpos,
					  lambda_vert,
					  lambda_mom,
					  ltarget_dist);

	L_p_id_container.push_back(it1);
	L_pi_id_container.push_back(it2);
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

	//Xi- searching
	for(int it3=0;it3<ntTpc;it3++){
	  if(it3==it2 || it3==it1) continue;
	  if(event.isElectron[it3]==1) continue;
	  if(event.isK18[it3]==1) continue;
	  if(event.isKurama[it3]==1) continue;
	  if(event.isBeam[it3]==1) continue;
	  if(event.isXi[it3]==1) continue;
	  if(event.isAccidental[it3]==1) continue;
	  if((event.pid[it3]&1)!=1) continue; //select pi like
	  //if((event.pid[it3]&4)==4) continue; //veto p-like
	  if(event.charge[it3]!=-1) continue;
	  if((event.pid[it3]&2)==2 && event.nsigma_pion[it3] > 5.) continue; //veto K- like

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
	  xi_p_container.push_back(it1);
	  xi_pi_container.push_back(it2);
	  xi_pi2_container.push_back(it3);

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
  std::cout<<"Xi/L candidates searching ends"<<std::endl;
#endif

#if DebugDisp
  std::cout<<"2. gamma (e+e- pair) candidate searching"<<std::endl;
#endif

  //e+&e- pair
  std::vector<Int_t> gamma_ep_id_container, gamma_em_id_container;
  std::vector<Double_t> gamma_epemdist_container;
  std::vector<TVector3> gamma_mom_container, gamma_ep_mom_container, gamma_em_mom_container;
  std::vector<TVector3> gamma_decayvtx_container;
  {
    for(Int_t it1=0;it1<ntTpc;it1++){
      if(event.isK18[it1]==1) continue;
      if(event.isKurama[it1]==1) continue;
      if(event.isBeam[it1]==1) continue;
      if(event.isXi[it1]==1) continue;
      if(event.isAccidental[it1]==1) continue;
      if((event.pid[it1]&4)==4) continue; //not proton
      if(event.isElectron[it1]!=1 && TMath::Abs(event.nsigma_electron[it1]) > 3) continue;
      if(event.charge[it1]!=1) continue;
      Int_t repid_ep = 0;
      if(event.isElectron[it1]!=1){
	Int_t flag = 1;
	for(Int_t i=0;i<3;i++){
	  Int_t temp = flag&event.pid[it1];
	  if(temp==flag) repid_ep += 1;
	  flag*=2;
	}
      }
      if(!GFTrackCont.TrackCheck(it1, repid_ep)) continue;

      Double_t ep_par[5];
      ep_par[0] = event.helix_cx[it1];
      ep_par[1] = event.helix_cy[it1];
      ep_par[2] = event.helix_z0[it1];
      ep_par[3] = event.helix_r[it1];
      ep_par[4] = event.helix_dz[it1];
      Int_t ep_nh = event.helix_t[it1].size();
      Double_t ep_theta_min = event.helix_t[it1][0] - vtx_scan_range/ep_par[3];
      Double_t ep_theta_max = TMath::Min(event.helix_t[it1][0] + vtx_scan_rangeInsideL/ep_par[3], event.helix_t[it1][ep_nh-1]);
      TVector3 ep_start = TVector3(event.calpos_x[it1][0], event.calpos_y[it1][0], event.calpos_z[it1][0]);
      TVector3 ep_end = TVector3(event.calpos_x[it1][ep_nh-1], event.calpos_y[it1][ep_nh-1], event.calpos_z[it1][ep_nh-1]);
      for(Int_t it2=0;it2<ntTpc;it2++){
	if(it1==it2) continue;
	if(event.isK18[it2]==1) continue;
	if(event.isKurama[it2]==1) continue;
	if(event.isBeam[it2]==1) continue;
	if(event.isXi[it2]==1) continue;
	if(event.isAccidental[it2]==1) continue;
	if((event.pid[it2]&4)==4) continue; //not proton
	if(event.isElectron[it2]!=1 && TMath::Abs(event.nsigma_electron[it2]) > 3) continue;
	if(event.charge[it2]!=-1) continue;
	Int_t repid_em = 0;
	if(event.isElectron[it2]!=1){
	  Int_t flag = 1;
	  for(Int_t i=0;i<3;i++){
	    Int_t temp = flag&event.pid[it2];
	    if(temp==flag) repid_em += 1;
	    flag*=2;
	  }
	}
	if(!GFTrackCont.TrackCheck(it2, repid_em)) continue;

	Double_t em_par[5];
	em_par[0] = event.helix_cx[it2];
	em_par[1] = event.helix_cy[it2];
	em_par[2] = event.helix_z0[it2];
	em_par[3] = event.helix_r[it2];
	em_par[4] = event.helix_dz[it2];

	Int_t em_nh = event.helix_t[it2].size();
	Double_t em_theta_min = TMath::Max(event.helix_t[it2][0] - vtx_scan_rangeInsideL/em_par[3], event.helix_t[it2][em_nh-1]);
	Double_t em_theta_max = event.helix_t[it2][0] + vtx_scan_range/em_par[3];
	TVector3 em_start = TVector3(event.calpos_x[it2][0], event.calpos_y[it2][0], event.calpos_z[it2][0]);
	TVector3 em_end = TVector3(event.calpos_x[it2][em_nh-1], event.calpos_y[it2][em_nh-1], event.calpos_z[it2][em_nh-1]);

	Double_t epem_dist = 10000.;
	TVector3 ep_mom; TVector3 em_mom; TVector3 epem_momsum;
	TVector3 gamma_vert =
	  Kinematics::LambdaVertex(dMagneticField, ep_par, em_par,
				   ep_theta_min, ep_theta_max,
				   em_theta_min, em_theta_max,
				   ep_mom, em_mom, epem_momsum, epem_dist);
	if(TMath::IsNaN(epem_dist)) continue;

	TLorentzVector Lem(em_mom, TMath::Hypot(em_mom.Mag(), ElectronMass));
	TLorentzVector Lep(ep_mom, TMath::Hypot(ep_mom.Mag(), ElectronMass));
	TLorentzVector Lepem = Lep + Lem;
	Double_t gamma_e = Lepem.E();
	TVector3 gamma_direction = epem_momsum.Unit();
	TVector3 gamma_mom = gamma_e*gamma_direction;

	if(gamma_e > gamma_ecut) continue;
	if(TMath::Abs(gamma_vert.x()) > 250. ||
	   TMath::Abs(gamma_vert.z()) > 250. ||
	   TMath::Abs(gamma_vert.y()) > 250.) continue; //Vertex cut

	Double_t ep_vertex_dist; Double_t em_vertex_dist;
	if(!Kinematics::HelixDirection(gamma_vert, ep_start, ep_end, ep_vertex_dist) ||
	   !Kinematics::HelixDirection(gamma_vert, em_start, em_end, em_vertex_dist)) continue;
	if(ep_vertex_dist > e_vtx_distcut) continue;
	if(em_vertex_dist > e_vtx_distcut) continue;

	Double_t gammatarget_dist;
	TVector3 gammatarget_vtx =
	  Kinematics::CalcCloseDistLambda(tgtpos,
					  gamma_vert,
					  gamma_mom,
					  gammatarget_dist);
	if(gammatarget_dist > gammatarget_distcut) continue;
	if(epem_dist < epem_distcut){
	  gamma_decayvtx_container.push_back(gamma_vert);
	  gamma_ep_id_container.push_back(it1);
	  gamma_em_id_container.push_back(it2);
	  gamma_epemdist_container.push_back(epem_dist);
	  gamma_mom_container.push_back(gamma_mom);
	  gamma_ep_mom_container.push_back(ep_mom);
	  gamma_em_mom_container.push_back(em_mom);
	}
      } //it2
    } //it1
  }

#if DebugDisp
  std::cout<<"gamma (e+e- pair) candidate searching ends"<<std::endl;
#endif

  std::vector<Double_t> GFL_mass_container(l_candidates, qnan);
  std::vector<Double_t> GFL_ppidist_container(l_candidates, qnan);
  std::vector<Double_t> GFL_targetdist_container(l_candidates, qnan);
  std::vector<TVector3> GFL_targetvtx_container(l_candidates, qnan_vec);
  std::vector<Double_t> GFL_targetcenterdist_container(l_candidates, qnan);
  std::vector<TVector3> GFL_targetcentervtx_container(l_candidates, qnan_vec);
  std::vector<TVector3> GFL_mom_container(l_candidates, qnan_vec);
  std::vector<TVector3> GFL_vtx_container(l_candidates, qnan_vec);

  std::vector<Int_t> GFL_p_id_container(l_candidates, qnan);
  std::vector<Int_t> GFL_pi_id_container(l_candidates, qnan);
  std::vector<Int_t> GFL_p_repid_container(l_candidates, qnan);
  std::vector<Int_t> GFL_pi_repid_container(l_candidates, qnan);
  std::vector<TVector3> GFL_p_mom_container(l_candidates, qnan_vec);
  std::vector<TVector3> GFL_pi_mom_container(l_candidates, qnan_vec);
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
  std::vector<Double_t> GFL_p_invbeta_container(l_candidates, qnan);
  std::vector<Double_t> GFL_pi_invbeta_container(l_candidates, qnan);
  if(l_candidates>0){
#if DebugDisp
    std::cout<<"3. Single L, LL candidates searching starts"<<std::endl;
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

      Bool_t vtxcut =
	(GFTrackCont.FindVertex(p_id, pi_id,
				p_repid, pi_repid,
				p_extrapolation,
				pi_extrapolation,
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
      TVector3 l_pos_tgtcenter =
	Kinematics::LambdaTargetCenter(l_vertex, l_mom, l_targetcenter_dist);
      if(TMath::Abs(l_pos_tgtcenter.y()) > GFltarget_ycut) continue;

      Int_t hitid_htof; Double_t tof_htof; Double_t tracklen_htof; TVector3 pos_htof; Double_t track2tgt_dist;
      Bool_t htofextrapolation =
	GFTrackCont.TPCHTOFTrackMatching(p_id, p_repid, l_vertex,
					 event.HtofSeg, event.posHtof,
					 hitid_htof, tof_htof,
					 tracklen_htof, pos_htof, track2tgt_dist);
      if(htofextrapolation){
	GFL_p_htofid_container[idp] = hitid_htof;
	GFL_p_tracklen_container[idp] = tracklen_htof;
	GFL_p_tof_container[idp] = event.tHtof[hitid_htof] - l_tof;
	GFL_p_mass2_container[idp] =
	  Kinematics::MassSquare(p_mom.Mag(), tracklen_htof, event.tHtof[hitid_htof] - l_tof);
	GFL_p_invbeta_container[idp] =
	  MathTools::C()*(event.tHtof[hitid_htof] - l_tof)/tracklen_htof;
      }

      htofextrapolation =
	GFTrackCont.TPCHTOFTrackMatching(pi_id, pi_repid, l_vertex,
					 event.HtofSeg, event.posHtof,
					 hitid_htof, tof_htof,
					 tracklen_htof, pos_htof, track2tgt_dist);
      if(htofextrapolation){
	GFL_pi_htofid_container[idp] = hitid_htof;
	GFL_pi_tracklen_container[idp] = tracklen_htof;
	GFL_pi_tof_container[idp] = event.tHtof[hitid_htof] - l_tof;
	GFL_pi_mass2_container[idp] =
	  Kinematics::MassSquare(pi_mom.Mag(), tracklen_htof, event.tHtof[hitid_htof] - l_tof);
	GFL_pi_invbeta_container[idp] =
	  MathTools::C()*(event.tHtof[hitid_htof] - l_tof)/tracklen_htof;

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

#if DebugDisp
    std::cout<<"Fitting with Genfit and Kinematic fitting for L/LL candidates"<<std::endl;
#endif

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
    std::vector<std::vector<Double_t>> GFLLdecays_invbeta_container;
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
	std::vector<Double_t> GFinvbeta_decays(4);

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
	GFinvbeta_decays[0] = GFL_p_invbeta_container[idp1];
	GFinvbeta_decays[1] = GFL_pi_invbeta_container[idp1];
	GFinvbeta_decays[2] = GFL_p_invbeta_container[idp2];
	GFinvbeta_decays[3] = GFL_pi_invbeta_container[idp2];

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
	GFLLdecays_invbeta_container.push_back(GFinvbeta_decays);

	const Int_t ntrack_ll = 4;
	event.KFlldecays_mom.resize(ntrack_ll);
	event.KFlldecays_mom_x.resize(ntrack_ll);
	event.KFlldecays_mom_y.resize(ntrack_ll);
	event.KFlldecays_mom_z.resize(ntrack_ll);
	event.KFlldecays_CMmom.resize(ntrack_ll);
	event.KFlldecays_CMmom_x.resize(ntrack_ll);
	event.KFlldecays_CMmom_y.resize(ntrack_ll);
	event.KFlldecays_CMmom_z.resize(ntrack_ll);

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
	TVector3 KFkk_ll_vertex_l1 =
	  Kinematics::MultitrackVertex(3,
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
	TVector3 KFkk_ll_vertex_l2 =
	  Kinematics::MultitrackVertex(3,
				       KFx0_l2, KFy0_l2,
				       KFu0_l2, KFv0_l2,
				       res_x_l2, res_y_l2,
				       res_u_l2, res_v_l2);

	/*
	  Double_t GFll_dist;
	  TVector3 GFll_vtx1, GFll_vtx2;
	  TVector3 GFll_vtx
	  = Kinematics::LambdaLambdaVertex(GFL_vtx_container.at(idp1),
	  GFL_mom_container.at(idp1),
	  GFL_vtx_container.at(idp2),
	  GFL_mom_container.at(idp2),
	  GFll_vtx1, GFll_vtx2, GFll_dist);
	  Double_t diff = GFll_dist;
	*/
	//Double_t diff = TMath::Hypot(GFlambda_mass[0] - LambdaMass, GFlambda_mass[1] - LambdaMass);
	//Double_t diff = TMath::Hypot(L_mass_container[idp1] - LambdaMass, L_mass_container[idp2] - LambdaMass);
	if(prev_LLfit_chisqr > chisqr_KFkk_ll_vertex){
	  event.llflag = true;

	  prev_LLfit_chisqr = chisqr_KFkk_ll_vertex;
	  LL_best = LLcount;

	  event.KFllvtx_x = KFll_vtx.x();
	  event.KFllvtx_y = KFll_vtx.y();
	  event.KFllvtx_z = KFll_vtx.z();
	  event.KFlldist = KFll_dist;

	  event.KFlmom1 = KFTVLd1.Mag();
	  event.KFlmom_x1 = KFTVLd1.x();
	  event.KFlmom_y1 = KFTVLd1.y();
	  event.KFlmom_z1 = KFTVLd1.z();
	  event.KFlmom2 = KFTVLd2.Mag();
	  event.KFlmom_x2 = KFTVLd2.x();
	  event.KFlmom_y2 = KFTVLd2.y();
	  event.KFlmom_z2 = KFTVLd2.z();
	  event.KFlchisqr1 = KFchisqrl1;
	  event.KFlchisqr2 = KFchisqrl2;

	  event.KFprodvtx_chisqr_ll = chisqr_KFkk_ll_vertex;
	  event.KFprodvtx_x_ll = KFkk_ll_vertex.x();
	  event.KFprodvtx_y_ll = KFkk_ll_vertex.y();
	  event.KFprodvtx_z_ll = KFkk_ll_vertex.z();

	  event.KFlprodvtx_x1 = KFprodvtx_closest1.x();
	  event.KFlprodvtx_y1 = KFprodvtx_closest1.y();
	  event.KFlprodvtx_z1 = KFprodvtx_closest1.z();
	  event.KFlprodvtx_dist1 = KFprodvtx_closedist1;

	  event.KFlprodvtx_x2 = KFprodvtx_closest2.x();
	  event.KFlprodvtx_y2 = KFprodvtx_closest2.y();
	  event.KFlprodvtx_z2 = KFprodvtx_closest2.z();
	  event.KFlprodvtx_dist2 = KFprodvtx_closedist2;

	  event.KFlldecays_mom[0] = KFTVP1.Mag();
	  event.KFlldecays_mom_x[0] = KFTVP1.x();
	  event.KFlldecays_mom_y[0] = KFTVP1.y();
	  event.KFlldecays_mom_z[0] = KFTVP1.z();
	  event.KFlldecays_mom[1] = KFTVPi1.Mag();
	  event.KFlldecays_mom_x[1] = KFTVPi1.x();
	  event.KFlldecays_mom_y[1] = KFTVPi1.y();
	  event.KFlldecays_mom_z[1] = KFTVPi1.z();
	  event.KFlldecays_mom[2] = KFTVP2.Mag();
	  event.KFlldecays_mom_x[2] = KFTVP2.x();
	  event.KFlldecays_mom_y[2] = KFTVP2.y();
	  event.KFlldecays_mom_z[2] = KFTVP2.z();
	  event.KFlldecays_mom[3] = KFTVPi2.Mag();
	  event.KFlldecays_mom_x[3] = KFTVPi2.x();
	  event.KFlldecays_mom_y[3] = KFTVPi2.y();
	  event.KFlldecays_mom_z[3] = KFTVPi2.z();

	  event.KFprodvtx_x_l1 = KFkk_ll_vertex_l1.x();
	  event.KFprodvtx_y_l1 = KFkk_ll_vertex_l1.y();
	  event.KFprodvtx_z_l1 = KFkk_ll_vertex_l1.z();
	  event.KFprodvtx_x_l2 = KFkk_ll_vertex_l2.x();
	  event.KFprodvtx_y_l2 = KFkk_ll_vertex_l2.y();
	  event.KFprodvtx_z_l2 = KFkk_ll_vertex_l2.z();
	}
	LLcount++;
      } //L combi1
    } //L combi2

    //for the LL event
    if(event.llflag){
#if DebugDisp
      std::cout<<"4. Detemine the best LL pair and save all"<<std::endl;
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

      //another p, pi pairing
      for(Int_t idll=0;idll<LLcount;idll++){
	if(idll==LL_best) continue;

	Int_t L1_id = GFLLid_container[LL_best][0];
	Int_t L2_id = GFLLid_container[LL_best][1];

	Int_t alter_L1_id = GFLLid_container[idll][0];
	Int_t alter_L2_id = GFLLid_container[idll][1];
	L_p_id_container.at(L1_id);
	L_p_id_container.at(L2_id);

	Bool_t another_ppi_pairing = false;
	if((L_p_id_container.at(L1_id) == L_p_id_container.at(alter_L2_id)) &&
	   (L_p_id_container.at(L2_id) == L_p_id_container.at(alter_L1_id)) &&
	   (L_pi_id_container.at(L1_id) == L_pi_id_container.at(alter_L1_id)) &&
	   (L_pi_id_container.at(L2_id) == L_pi_id_container.at(alter_L2_id))) another_ppi_pairing = true;

	if((L_p_id_container.at(L1_id) == L_p_id_container.at(alter_L1_id)) &&
	   (L_p_id_container.at(L2_id) == L_p_id_container.at(alter_L2_id)) &&
	   (L_pi_id_container.at(L1_id) == L_pi_id_container.at(alter_L2_id)) &&
	   (L_pi_id_container.at(L2_id) == L_pi_id_container.at(alter_L1_id))) another_ppi_pairing = true;

	if(!another_ppi_pairing) continue;

	Int_t LL_alter = idll;
	event.GFlmass_alter1 = GFLLmass_container[LL_alter][0];
	event.GFldecayvtx_x_alter1 = GFLL_Ldecayvtx_container[LL_alter][0].x();
	event.GFldecayvtx_y_alter1 = GFLL_Ldecayvtx_container[LL_alter][0].y();
	event.GFldecayvtx_z_alter1 = GFLL_Ldecayvtx_container[LL_alter][0].z();
	event.GFlmom_alter1 = GFLLmom_container[LL_alter][0].Mag();
	event.GFlmom_x_alter1 = GFLLmom_container[LL_alter][0].x();
	event.GFlmom_y_alter1 = GFLLmom_container[LL_alter][0].y();
	event.GFlmom_z_alter1 = GFLLmom_container[LL_alter][0].z();
	event.GFppi_dist_alter1 = GFLLppidist_container[LL_alter][0];
	event.GFltarget_dist_alter1 = GFLL_Ltargetdist_container[LL_alter][0];
	event.GFltargetvtx_x_alter1 = GFLL_Ltargetvtx_container[LL_alter][0].x();
	event.GFltargetvtx_y_alter1 = GFLL_Ltargetvtx_container[LL_alter][0].y();
	event.GFltargetvtx_z_alter1 = GFLL_Ltargetvtx_container[LL_alter][0].z();
	event.GFltargetcenter_dist_alter1 = GFLL_Ltargetcenterdist_container[LL_alter][0];
	event.GFltargetcenter_x_alter1 = GFLL_Ltargetcentervtx_container[LL_alter][0].x();
	event.GFltargetcenter_y_alter1 = GFLL_Ltargetcentervtx_container[LL_alter][0].y();
	event.GFltargetcenter_z_alter1 = GFLL_Ltargetcentervtx_container[LL_alter][0].z();

	event.GFlmass_alter2 = GFLLmass_container[LL_alter][1];
	event.GFldecayvtx_x_alter2 = GFLL_Ldecayvtx_container[LL_alter][1].x();
	event.GFldecayvtx_y_alter2 = GFLL_Ldecayvtx_container[LL_alter][1].y();
	event.GFldecayvtx_z_alter2 = GFLL_Ldecayvtx_container[LL_alter][1].z();
	event.GFlmom_alter2 = GFLLmom_container[LL_alter][1].Mag();
	event.GFlmom_x_alter2 = GFLLmom_container[LL_alter][1].x();
	event.GFlmom_y_alter2 = GFLLmom_container[LL_alter][1].y();
	event.GFlmom_z_alter2 = GFLLmom_container[LL_alter][1].z();
	event.GFppi_dist_alter2 = GFLLppidist_container[LL_alter][1];
	event.GFltarget_dist_alter2 = GFLL_Ltargetdist_container[LL_alter][1];
	event.GFltargetvtx_x_alter2 = GFLL_Ltargetvtx_container[LL_alter][1].x();
	event.GFltargetvtx_y_alter2 = GFLL_Ltargetvtx_container[LL_alter][1].y();
	event.GFltargetvtx_z_alter2 = GFLL_Ltargetvtx_container[LL_alter][1].z();
	event.GFltargetcenter_dist_alter2 = GFLL_Ltargetcenterdist_container[LL_alter][1];
	event.GFltargetcenter_x_alter2 = GFLL_Ltargetcentervtx_container[LL_alter][1].x();
	event.GFltargetcenter_y_alter2 = GFLL_Ltargetcentervtx_container[LL_alter][1].y();
	event.GFltargetcenter_z_alter2 = GFLL_Ltargetcentervtx_container[LL_alter][1].z();
      }

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
      TVector3 GFprodvtx_closest1 =
	Kinematics::CalcCloseDistLambda(GFkk_ll_vertex,
					GFLL_Ldecayvtx_container[LL_best][0],
					GFLLmom_container[LL_best][0],
					GFprodvtx_closedist1);

      Double_t GFprodvtx_closedist2 = qnan;
      TVector3 GFprodvtx_closest2 =
	Kinematics::CalcCloseDistLambda(GFkk_ll_vertex,
					GFLL_Ldecayvtx_container[LL_best][1],
					GFLLmom_container[LL_best][1],
					GFprodvtx_closedist2);

#if DebugDisp
      std::cout<<"K-K+LL vertex "<<GFkk_ll_vertex
	       <<" Lambda1's closest point to the vertex"<<GFprodvtx_closest1
	       <<" Lambda2's closest point to the vertex"<<GFprodvtx_closest2
	       <<std::endl;
#endif

      event.GFprodvtx_x_ll = GFkk_ll_vertex.x();
      event.GFprodvtx_y_ll = GFkk_ll_vertex.y();
      event.GFprodvtx_z_ll = GFkk_ll_vertex.z();

      event.GFlprodvtx_x1 = GFprodvtx_closest1.x();
      event.GFlprodvtx_y1 = GFprodvtx_closest1.y();
      event.GFlprodvtx_z1 = GFprodvtx_closest1.z();
      event.GFlprodvtx_dist1 = GFprodvtx_closedist1;

      event.GFlprodvtx_x2 = GFprodvtx_closest2.x();
      event.GFlprodvtx_y2 = GFprodvtx_closest2.y();
      event.GFlprodvtx_z2 = GFprodvtx_closest2.z();
      event.GFlprodvtx_dist2 = GFprodvtx_closedist2;

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
      Int_t GFntdecays = 4;
      event.GFlldecays_pdgcode.resize(GFntdecays);
      event.GFlldecays_nhtrack.resize(GFntdecays);
      event.GFlldecays_charge.resize(GFntdecays);
      event.GFlldecays_chisqr.resize(GFntdecays);
      event.GFlldecays_tof.resize(GFntdecays);
      event.GFlldecays_pval.resize(GFntdecays);
      event.GFlldecays_htofid.resize(GFntdecays);
      event.GFlldecays_tracklen.resize(GFntdecays);
      event.GFlldecays_tof.resize(GFntdecays);
      event.GFlldecays_mass2.resize(GFntdecays);
      event.GFlldecays_invbeta.resize(GFntdecays);
      event.GFlldecays_mom.resize(GFntdecays);
      event.GFlldecays_mom_x.resize(GFntdecays);
      event.GFlldecays_mom_y.resize(GFntdecays);
      event.GFlldecays_mom_z.resize(GFntdecays);
      event.GFlldecays_CMmom.resize(GFntdecays);
      event.GFlldecays_CMmom_x.resize(GFntdecays);
      event.GFlldecays_CMmom_y.resize(GFntdecays);
      event.GFlldecays_CMmom_z.resize(GFntdecays);
      event.GFlldecays_momloss.resize(GFntdecays);
      event.GFlldecays_eloss.resize(GFntdecays);

      event.lldecays_id.push_back(GFLLdecays_trackid_container[LL_best][0]);
      event.lldecays_id.push_back(GFLLdecays_trackid_container[LL_best][1]);
      event.lldecays_id.push_back(GFLLdecays_trackid_container[LL_best][2]);
      event.lldecays_id.push_back(GFLLdecays_trackid_container[LL_best][3]);
      event.lldecays_mom.push_back(L_p_mom_container[L1_id].Mag());
      event.lldecays_mom.push_back(L_pi_mom_container[L1_id].Mag());
      event.lldecays_mom.push_back(L_p_mom_container[L2_id].Mag());
      event.lldecays_mom.push_back(L_pi_mom_container[L2_id].Mag());
      event.lldecays_mom_x.push_back(L_p_mom_container[L1_id].x());
      event.lldecays_mom_x.push_back(L_pi_mom_container[L1_id].x());
      event.lldecays_mom_x.push_back(L_p_mom_container[L2_id].x());
      event.lldecays_mom_x.push_back(L_pi_mom_container[L2_id].x());
      event.lldecays_mom_y.push_back(L_p_mom_container[L1_id].y());
      event.lldecays_mom_y.push_back(L_pi_mom_container[L1_id].y());
      event.lldecays_mom_y.push_back(L_p_mom_container[L2_id].y());
      event.lldecays_mom_y.push_back(L_pi_mom_container[L2_id].y());
      event.lldecays_mom_z.push_back(L_p_mom_container[L1_id].z());
      event.lldecays_mom_z.push_back(L_pi_mom_container[L1_id].z());
      event.lldecays_mom_z.push_back(L_p_mom_container[L2_id].z());
      event.lldecays_mom_z.push_back(L_pi_mom_container[L2_id].z());

      {
	TLorentzVector Lv_decays1(L_p_mom_container[L1_id],
				  TMath::Hypot(L_p_mom_container[L1_id].Mag(), ProtonMass));
	TLorentzVector Lv_decays2(L_pi_mom_container[L1_id],
				  TMath::Hypot(L_pi_mom_container[L1_id].Mag(), PionMass));
	TLorentzVector Lv_decays3(L_p_mom_container[L2_id],
				  TMath::Hypot(L_p_mom_container[L2_id].Mag(), ProtonMass));
	TLorentzVector Lv_decays4(L_pi_mom_container[L2_id],
				  TMath::Hypot(L_pi_mom_container[L2_id].Mag(), PionMass));
	TLorentzVector LvL1(L_mom_container.at(L1_id),
			    TMath::Hypot(L_mom_container.at(L1_id).Mag(), LambdaMass));
	TLorentzVector LvL2(L_mom_container.at(L2_id),
			    TMath::Hypot(L_mom_container.at(L2_id).Mag(), LambdaMass));
	TVector3 BoostL1 = LvL1.BoostVector();
	TVector3 BoostL2 = LvL2.BoostVector();
	Lv_decays1.Boost(-BoostL1);
	Lv_decays2.Boost(-BoostL1);
	Lv_decays3.Boost(-BoostL2);
	Lv_decays4.Boost(-BoostL2);

	event.lldecays_CMmom.push_back(Lv_decays1.P());
	event.lldecays_CMmom_x.push_back(Lv_decays1.Px());
	event.lldecays_CMmom_y.push_back(Lv_decays1.Py());
	event.lldecays_CMmom_z.push_back(Lv_decays1.Pz());
	event.lldecays_CMmom.push_back(Lv_decays2.P());
	event.lldecays_CMmom_x.push_back(Lv_decays2.Px());
	event.lldecays_CMmom_y.push_back(Lv_decays2.Py());
	event.lldecays_CMmom_z.push_back(Lv_decays2.Pz());
	event.lldecays_CMmom.push_back(Lv_decays3.P());
	event.lldecays_CMmom_x.push_back(Lv_decays3.Px());
	event.lldecays_CMmom_y.push_back(Lv_decays3.Py());
	event.lldecays_CMmom_z.push_back(Lv_decays3.Pz());
	event.lldecays_CMmom.push_back(Lv_decays4.P());
	event.lldecays_CMmom_x.push_back(Lv_decays4.Px());
	event.lldecays_CMmom_y.push_back(Lv_decays4.Py());
	event.lldecays_CMmom_z.push_back(Lv_decays4.Pz());
      } //ntdecays

      for(int j=0;j<GFntdecays;j++){
	Int_t igf = GFLLdecays_trackid_container[LL_best][j];
	Int_t repid = GFLLdecays_repid_container[LL_best][j];
	event.GFlldecays_nhtrack[j] = GFTrackCont.GetNHits(igf);
	event.GFlldecays_chisqr[j] = GFTrackCont.GetChi2NDF(igf, repid);
	event.GFlldecays_charge[j] = GFTrackCont.GetCharge(igf, repid);
	event.GFlldecays_tof[j] = GFTrackCont.GetTrackTOF(igf, 0, -1, repid);
	event.GFlldecays_pval[j] = GFTrackCont.GetPvalue(igf, repid);
	event.GFlldecays_pdgcode[j] = GFTrackCont.GetPDGcode(igf, repid);

	event.GFlldecays_htofid[j] = GFLLdecays_htofid_container[LL_best][j];
	event.GFlldecays_tracklen[j] = GFLLdecays_tracklen_container[LL_best][j];
	event.GFlldecays_tof[j] = GFLLdecays_tof_container[LL_best][j];
	event.GFlldecays_mass2[j] = GFLLdecays_mass2_container[LL_best][j];
	event.GFlldecays_invbeta[j] = GFLLdecays_invbeta_container[LL_best][j];
	event.GFlldecays_mom[j] = GFLLdecays_mom_container[LL_best][j].Mag();
	event.GFlldecays_mom_x[j] = GFLLdecays_mom_container[LL_best][j].x();
	event.GFlldecays_mom_y[j] = GFLLdecays_mom_container[LL_best][j].y();
	event.GFlldecays_mom_z[j] = GFLLdecays_mom_container[LL_best][j].z();
	event.GFlldecays_momloss[j] = qnan;
	event.GFlldecays_eloss[j] = qnan;

	Double_t mass = PionMass;
	if(j==0 || j==2) mass = ProtonMass;

	TLorentzVector GFLv_decays(GFLLdecays_mom_container[LL_best][j],
				   TMath::Hypot(event.GFlldecays_mom[j], mass));

	TVector3 BoostL1 = LvL1_fixedmass.BoostVector();
	TVector3 BoostL2 = LvL2_fixedmass.BoostVector();
	if(j<2) GFLv_decays.Boost(-BoostL1);
	else GFLv_decays.Boost(-BoostL2);

	event.GFlldecays_CMmom[j] = GFLv_decays.P();
	event.GFlldecays_CMmom_x[j] = GFLv_decays.Px();
	event.GFlldecays_CMmom_y[j] = GFLv_decays.Py();
	event.GFlldecays_CMmom_z[j] = GFLv_decays.Pz();
      } //GFntdecays

      TVector3 KFTVLd1(event.KFlmom_x1, event.KFlmom_y1, event.KFlmom_z1);
      TVector3 KFTVLd2(event.KFlmom_x2, event.KFlmom_y2, event.KFlmom_z2);
      TVector3 KFTVP1(event.KFlldecays_mom_x[0], event.KFlldecays_mom_y[0], event.KFlldecays_mom_z[0]);
      TVector3 KFTVPi1(event.KFlldecays_mom_x[1], event.KFlldecays_mom_y[1], event.KFlldecays_mom_z[1]);
      TVector3 KFTVP2(event.KFlldecays_mom_x[2], event.KFlldecays_mom_y[2], event.KFlldecays_mom_z[2]);
      TVector3 KFTVPi2(event.KFlldecays_mom_x[3], event.KFlldecays_mom_y[3], event.KFlldecays_mom_z[3]);
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
      TVector3 BoostKFLd1 = KFLvL1_fixedmass.BoostVector();
      TVector3 BoostKFLd2 = KFLvL2_fixedmass.BoostVector();
      KFLv_p1.Boost(-BoostKFLd1);
      KFLv_pi1.Boost(-BoostKFLd1);
      KFLv_p2.Boost(-BoostKFLd2);
      KFLv_pi2.Boost(-BoostKFLd2);

      event.KFlldecays_CMmom[0] = KFLv_p1.P();
      event.KFlldecays_CMmom_x[0] = KFLv_p1.Px();
      event.KFlldecays_CMmom_y[0] = KFLv_p1.Py();
      event.KFlldecays_CMmom_z[0] = KFLv_p1.Pz();
      event.KFlldecays_CMmom[1] = KFLv_pi1.P();
      event.KFlldecays_CMmom_x[1] = KFLv_pi1.Px();
      event.KFlldecays_CMmom_y[1] = KFLv_pi1.Py();
      event.KFlldecays_CMmom_z[1] = KFLv_pi1.Pz();
      event.KFlldecays_CMmom[2] = KFLv_p2.P();
      event.KFlldecays_CMmom_x[2] = KFLv_p2.Px();
      event.KFlldecays_CMmom_y[2] = KFLv_p2.Py();
      event.KFlldecays_CMmom_z[2] = KFLv_p2.Pz();
      event.KFlldecays_CMmom[3] = KFLv_pi2.P();
      event.KFlldecays_CMmom_x[3] = KFLv_pi2.Px();
      event.KFlldecays_CMmom_y[3] = KFLv_pi2.Py();
      event.KFlldecays_CMmom_z[3] = KFLv_pi2.Pz();

      //Remaining tracks
#if DebugDisp
      std::cout<<"LL event, Save remaining tracks"<<std::endl;
#endif

      event.llp_multi = target_p_id_container.size();
      event.llpip_multi = target_pip_id_container.size();
      event.llpim_multi = target_pim_id_container.size();
      event.llep_multi = target_ep_id_container.size();
      event.llem_multi = target_em_id_container.size();
      event.llppip_multi = target_ppip_id_container.size();

      TVector3 closepoint, closepoint_mom;
      Double_t closepoint_tracklen, closepoint_tof;
      for(Int_t itp=0; itp<target_p_id_container.size(); ++itp){
	Int_t id_p = target_p_id_container[itp];
	if(id_p==event.lldecays_id[0] || id_p==event.lldecays_id[2]){
	  event.llp_multi -= 1;
	  continue;
	}
	if(target_p_dist2tgt_container[itp] > residual_track_distcut){
	  event.llp_multi -= 1;
	  continue;
	}
	event.llresidual_id.push_back(id_p);
	event.llresidual_dist2tgt.push_back(target_p_dist2tgt_container[itp]);
	event.llresidual_mass2.push_back(target_p_mass2_container[itp]);
	event.llresidual_invbeta.push_back(target_p_invbeta_container[itp]);
	event.llresidual_mom.push_back(target_p_mom_container[itp].Mag());
	event.llresidual_mom_x.push_back(target_p_mom_container[itp].x());
	event.llresidual_mom_y.push_back(target_p_mom_container[itp].y());
	event.llresidual_mom_z.push_back(target_p_mom_container[itp].z());
	event.llresidual_charge.push_back(1);

	if(GFTrackCont.ExtrapolateToPoint(id_p,
					  GFkk_ll_vertex,
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_p_repid_container[itp])){
	  closepoint -= GFkk_ll_vertex;
	  event.llresidual_GFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.llresidual_GFdist2prodvtx.push_back(qnan);
	if(GFTrackCont.ExtrapolateToPoint(id_p,
					  KFkk_ll_vertex,
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_p_repid_container[itp])){
	  closepoint -= KFkk_ll_vertex;
	  event.llresidual_KFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.llresidual_KFdist2prodvtx.push_back(qnan);
      }
      for(Int_t itpip=0; itpip<target_pip_id_container.size(); ++itpip){
	Int_t id_pip = target_pip_id_container[itpip];
	if(id_pip==event.lldecays_id[0] || id_pip==event.lldecays_id[2]){
	  event.llpip_multi -= 1;
	  continue;
	}
	if(target_pip_dist2tgt_container[itpip] > residual_track_distcut){
	  event.llpip_multi -= 1;
	  continue;
	}
	event.llresidual_id.push_back(id_pip);
	event.llresidual_dist2tgt.push_back(target_pip_dist2tgt_container[itpip]);
	event.llresidual_mass2.push_back(target_pip_mass2_container[itpip]);
	event.llresidual_invbeta.push_back(target_pip_invbeta_container[itpip]);
	event.llresidual_mom.push_back(target_pip_mom_container[itpip].Mag());
	event.llresidual_mom_x.push_back(target_pip_mom_container[itpip].x());
	event.llresidual_mom_y.push_back(target_pip_mom_container[itpip].y());
	event.llresidual_mom_z.push_back(target_pip_mom_container[itpip].z());
	event.llresidual_charge.push_back(1);
	if(GFTrackCont.ExtrapolateToPoint(id_pip,
					  GFkk_ll_vertex,
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_pip_repid_container[itpip])){
	  closepoint -= GFkk_ll_vertex;
	  event.llresidual_GFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.llresidual_GFdist2prodvtx.push_back(qnan);
	if(GFTrackCont.ExtrapolateToPoint(id_pip,
					  KFkk_ll_vertex,
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_pip_repid_container[itpip])){
	  closepoint -= KFkk_ll_vertex;
	  event.llresidual_KFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.llresidual_KFdist2prodvtx.push_back(qnan);
      }
      for(Int_t itpim=0; itpim<target_pim_id_container.size(); ++itpim){
	Int_t id_pim = target_pim_id_container[itpim];
	if(id_pim==event.lldecays_id[1] || id_pim==event.lldecays_id[3]){
	  event.llpim_multi -= 1;
	  continue;
	}
	if(target_pim_dist2tgt_container[itpim] > residual_track_distcut){
	  event.llpim_multi -= 1;
	  continue;
	}
	event.llresidual_id.push_back(id_pim);
	event.llresidual_dist2tgt.push_back(target_pim_dist2tgt_container[itpim]);
	event.llresidual_mass2.push_back(target_pim_mass2_container[itpim]);
	event.llresidual_invbeta.push_back(target_pim_invbeta_container[itpim]);
	event.llresidual_mom.push_back(target_pim_mom_container[itpim].Mag());
	event.llresidual_mom_x.push_back(target_pim_mom_container[itpim].x());
	event.llresidual_mom_y.push_back(target_pim_mom_container[itpim].y());
	event.llresidual_mom_z.push_back(target_pim_mom_container[itpim].z());
	event.llresidual_charge.push_back(-1);

	if(GFTrackCont.ExtrapolateToPoint(id_pim,
					  GFkk_ll_vertex,
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_pim_repid_container[itpim])){
	  closepoint -= GFkk_ll_vertex;
	  event.llresidual_GFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.llresidual_GFdist2prodvtx.push_back(qnan);
	if(GFTrackCont.ExtrapolateToPoint(id_pim,
					  KFkk_ll_vertex,
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_pim_repid_container[itpim])){
	  closepoint -= KFkk_ll_vertex;
	  event.llresidual_KFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.llresidual_KFdist2prodvtx.push_back(qnan);
      }
      for(Int_t itppip=0; itppip<target_ppip_id_container.size(); ++itppip){
	Int_t id_ppip = target_ppip_id_container[itppip];
	if(id_ppip==event.lldecays_id[0] || id_ppip==event.lldecays_id[2]){
	  event.llppip_multi -= 1;
	  continue;
	}
	if(target_ppip_dist2tgt_container[itppip] > residual_track_distcut){
	  event.llppip_multi -= 1;
	  continue;
	}
	event.llresidual_id.push_back(id_ppip);
	event.llresidual_dist2tgt.push_back(target_ppip_dist2tgt_container[itppip]);
	event.llresidual_mass2.push_back(target_ppip_mass2_container[itppip]);
	event.llresidual_invbeta.push_back(target_ppip_invbeta_container[itppip]);
	event.llresidual_mom.push_back(target_ppip_mom_container[itppip].Mag());
	event.llresidual_mom_x.push_back(target_ppip_mom_container[itppip].x());
	event.llresidual_mom_y.push_back(target_ppip_mom_container[itppip].y());
	event.llresidual_mom_z.push_back(target_ppip_mom_container[itppip].z());
	event.llresidual_charge.push_back(1);
	if(GFTrackCont.ExtrapolateToPoint(id_ppip,
					  GFkk_ll_vertex,
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_ppip_repid_container[itppip])){
	  closepoint -= GFkk_ll_vertex;
	  event.llresidual_GFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.llresidual_GFdist2prodvtx.push_back(qnan);
	if(GFTrackCont.ExtrapolateToPoint(id_ppip,
					  KFkk_ll_vertex,
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_ppip_repid_container[itppip])){
	  closepoint -= KFkk_ll_vertex;
	  event.llresidual_KFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.llresidual_KFdist2prodvtx.push_back(qnan);
      }
      for(Int_t itep=0; itep<target_ep_id_container.size(); ++itep){
	Int_t id_ep = target_ep_id_container[itep];
	if(id_ep==event.lldecays_id[0] || id_ep==event.lldecays_id[2]){
	  event.llep_multi -= 1;
	  continue;
	}
	if(target_ep_dist2tgt_container[itep] > residual_track_distcut){
	  event.llep_multi -= 1;
	  continue;
	}
	event.llresidual_id.push_back(id_ep);
	event.llresidual_dist2tgt.push_back(target_ep_dist2tgt_container[itep]);
	event.llresidual_mass2.push_back(target_ep_mass2_container[itep]);
	event.llresidual_invbeta.push_back(target_ep_invbeta_container[itep]);
	event.llresidual_mom.push_back(target_ep_mom_container[itep].Mag());
	event.llresidual_mom_x.push_back(target_ep_mom_container[itep].x());
	event.llresidual_mom_y.push_back(target_ep_mom_container[itep].y());
	event.llresidual_mom_z.push_back(target_ep_mom_container[itep].z());
	event.llresidual_charge.push_back(1);
	if(GFTrackCont.ExtrapolateToPoint(id_ep,
					  GFkk_ll_vertex,
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_ep_repid_container[itep])){
	  closepoint -= GFkk_ll_vertex;
	  event.llresidual_GFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.llresidual_GFdist2prodvtx.push_back(qnan);
	if(GFTrackCont.ExtrapolateToPoint(id_ep,
					  KFkk_ll_vertex,
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_ep_repid_container[itep])){
	  closepoint -= KFkk_ll_vertex;
	  event.llresidual_KFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.llresidual_KFdist2prodvtx.push_back(qnan);
      }
      for(Int_t item=0; item<target_em_id_container.size(); ++item){
	Int_t id_em = target_em_id_container[item];
	if(id_em==event.lldecays_id[1] || id_em==event.lldecays_id[3]){
	  event.llem_multi -= 1;
	  continue;
	}
	if(target_em_dist2tgt_container[item] > residual_track_distcut){
	  event.llem_multi -= 1;
	  continue;
	}
	event.llresidual_id.push_back(id_em);
	event.llresidual_dist2tgt.push_back(target_em_dist2tgt_container[item]);
	event.llresidual_mass2.push_back(target_em_mass2_container[item]);
	event.llresidual_invbeta.push_back(target_em_invbeta_container[item]);
	event.llresidual_mom.push_back(target_em_mom_container[item].Mag());
	event.llresidual_mom_x.push_back(target_em_mom_container[item].x());
	event.llresidual_mom_y.push_back(target_em_mom_container[item].y());
	event.llresidual_mom_z.push_back(target_em_mom_container[item].z());
	event.llresidual_charge.push_back(-1);
	if(GFTrackCont.ExtrapolateToPoint(id_em,
					  GFkk_ll_vertex,
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_em_repid_container[item])){
	  closepoint -= GFkk_ll_vertex;
	  event.llresidual_GFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.llresidual_GFdist2prodvtx.push_back(qnan);
	if(GFTrackCont.ExtrapolateToPoint(id_em,
					  KFkk_ll_vertex,
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_em_repid_container[item])){
	  closepoint -= KFkk_ll_vertex;
	  event.llresidual_KFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.llresidual_KFdist2prodvtx.push_back(qnan);
      }

      event.llresidual_multi =
	event.llep_multi + event.llem_multi +
	event.llpip_multi + event.llp_multi +
	event.llpim_multi + event.llppip_multi;

#if DebugDisp
      std::cout<<"LL event, Saving gamma candidates"<<std::endl;
#endif

      event.llg_multi = gamma_mom_container.size();
      for(Int_t itg=0; itg<gamma_mom_container.size(); ++itg){
	Int_t epid = gamma_ep_id_container[itg];
	Int_t emid = gamma_em_id_container[itg];
	if(epid==event.lldecays_id[0] || epid==event.lldecays_id[2] ||
	   emid==event.lldecays_id[1] || emid==event.lldecays_id[3]){
	  event.llg_multi -=1;
	  continue;
	}
	auto iter_ep = find(event.llepidgamma.begin(), event.llepidgamma.end(), epid);
	if(iter_ep != event.llepidgamma.end()){
	  event.llg_multi -=1;
	  continue;
	}
	auto iter_em = find(event.llemidgamma.begin(), event.llemidgamma.end(), emid);
	if(iter_em != event.llemidgamma.end()){
	  event.llg_multi -=1;
	  continue;
	}

	Double_t dist = gamma_epemdist_container[itg];
	TVector3 ep_mom = gamma_ep_mom_container[itg];
	TVector3 em_mom = gamma_em_mom_container[itg];
	TVector3 gamma_mom = gamma_mom_container[itg];
	TVector3 gamma_vtx = gamma_decayvtx_container[itg];

	event.llepidgamma.push_back(epid);
	event.llemidgamma.push_back(emid);
	event.llepmomgamma.push_back(ep_mom.Mag());
	event.llepmomgamma_x.push_back(ep_mom.x());
	event.llepmomgamma_y.push_back(ep_mom.y());
	event.llepmomgamma_z.push_back(ep_mom.z());
	event.llemmomgamma.push_back(em_mom.Mag());
	event.llemmomgamma_x.push_back(em_mom.x());
	event.llemmomgamma_y.push_back(em_mom.y());
	event.llemmomgamma_z.push_back(em_mom.z());
	event.llmomgamma.push_back(gamma_mom.Mag());
	event.llmomgamma_x.push_back(gamma_mom.x());
	event.llmomgamma_y.push_back(gamma_mom.y());
	event.llmomgamma_z.push_back(gamma_mom.z());
	event.llepidistgamma.push_back(dist);
	event.llvtxgamma_x.push_back(gamma_vtx.x());
	event.llvtxgamma_y.push_back(gamma_vtx.y());
	event.llvtxgamma_z.push_back(gamma_vtx.z());
      }
    } //LL flag
  } //l_candidates

#if DebugDisp
  std::cout<<"Genfit fitting & Kinematic fitting start"<<std::endl;
  std::cout<<"fitting with Genfit and Kinematic fitting for Xi- candidates"<<std::endl;
#endif

  std::vector<Int_t> GFxi_p_id_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_pi_id_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_pi2_id_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_p_rep_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_pi_rep_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_pi2_rep_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_gfid_container(xi_candidates, -1);
  std::vector<TVector3> GFxi_mom_container(xi_candidates, qnan_vec);
  std::vector<TVector3> GFxi_decayvertex_container(xi_candidates, qnan_vec);
  std::vector<Double_t> GFxi_mass_container(xi_candidates, qnan);
  std::vector<TVector3> GFxi_pos_targetcenter_container(xi_candidates, qnan_vec);
  std::vector<TVector3> GFxi_mom_targetcenter_container(xi_candidates, qnan_vec);
  std::vector<Double_t> GFxi_dist_targetcenter_container(xi_candidates, qnan);
  std::vector<TVector3> GFxi_l_mom_container(xi_candidates, qnan_vec);
  std::vector<TVector3> GFxi_l_vert_container(xi_candidates, qnan_vec);
  std::vector<Double_t> GFxi_l_mass_container(xi_candidates, qnan);
  std::vector<Double_t> GFxi_l_tracklen_container(xi_candidates, qnan);
  std::vector<Double_t> GFxi_l_tof_container(xi_candidates, qnan);
  std::vector<TVector3> GFxi_p_mom_container(xi_candidates, qnan_vec);
  std::vector<TVector3> GFxi_pi_mom_container(xi_candidates, qnan_vec);
  std::vector<TVector3> GFxi_pi2_mom_container(xi_candidates, qnan_vec);
  std::vector<Double_t> GFxi_p_mass2_container(xi_candidates, qnan);
  std::vector<Double_t> GFxi_pi_mass2_container(xi_candidates, qnan);
  std::vector<Double_t> GFxi_pi2_mass2_container(xi_candidates, qnan);
  std::vector<Double_t> GFxi_p_invbeta_container(xi_candidates, qnan);
  std::vector<Double_t> GFxi_pi_invbeta_container(xi_candidates, qnan);
  std::vector<Double_t> GFxi_pi2_invbeta_container(xi_candidates, qnan);

  std::vector<Double_t> GFxi_p_tof_container(xi_candidates, qnan);
  std::vector<Double_t> GFxi_pi_tof_container(xi_candidates, qnan);
  std::vector<Double_t> GFxi_pi2_tof_container(xi_candidates, qnan);
  std::vector<Double_t> GFxi_p_tracklen_container(xi_candidates, qnan);
  std::vector<Double_t> GFxi_pi_tracklen_container(xi_candidates, qnan);
  std::vector<Double_t> GFxi_pi2_tracklen_container(xi_candidates, qnan);
  std::vector<Int_t> GFxi_p_htofid_container(xi_candidates, qnan);
  std::vector<Int_t> GFxi_pi_htofid_container(xi_candidates, qnan);
  std::vector<Int_t> GFxi_pi2_htofid_container(xi_candidates, qnan);
  std::vector<Double_t> GFxi_ppi_closedist_container(xi_candidates, qnan);
  std::vector<Double_t> GFxi_lpi_closedist_container(xi_candidates, qnan);

  std::vector<Int_t> KFxi_gfid_container(xi_candidates, -1);
  std::vector<TVector3> KFxi_l_mom_container0(xi_candidates, qnan_vec);
  std::vector<TVector3> KFxi_decayvertex_container(xi_candidates, qnan_vec);
  std::vector<TVector3> KFxi_mom_container(xi_candidates, qnan_vec);
  std::vector<TVector3> KFxi_l_mom_container(xi_candidates, qnan_vec);
  std::vector<TVector3> KFxi_pos_targetcenter_container(xi_candidates, qnan_vec);
  std::vector<TVector3> KFxi_mom_targetcenter_container(xi_candidates, qnan_vec);
  std::vector<Double_t> KFxi_dist_targetcenter_container(xi_candidates, qnan);

  std::vector<Double_t> KFxi_lchisqr_container(xi_candidates, qnan);
  std::vector<Double_t> KFxi_lpval_container(xi_candidates, qnan);
  std::vector<Double_t> KFxi_chisqr_container(xi_candidates, qnan);
  std::vector<Double_t> KFxi_pval_container(xi_candidates, qnan);
  std::vector<Double_t> KFxi_mass_container(xi_candidates, qnan);
  std::vector<Double_t> KFxi_lpi_closedist_container(xi_candidates, qnan);
  std::vector<std::vector<Double_t>> KFxi_lpull_container(xi_candidates, std::vector<Double_t>(6, qnan));
  std::vector<std::vector<Double_t>> KFxi_pull_container(xi_candidates, std::vector<Double_t>(6, qnan));
  std::vector<TMatrixD> VXiContainer(xi_candidates, TMatrixD(3, 3));
  std::vector<TVector3> KFxi_p_mom_container(xi_candidates, qnan_vec);
  std::vector<TVector3> KFxi_pi_mom_container(xi_candidates, qnan_vec);
  std::vector<TVector3> KFxi_pi2_mom_container(xi_candidates, qnan_vec);

#if DebugDisp
  if(xi_candidates > 0) std::cout<<"5. Detemine the best Xi candidate and save all"<<std::endl;
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
    if(!GFTrackCont.FindVertexXi(trackid_pi2, repid_pi2,
				 GFlambda_vert, GFlambda_mom, GFlambda_tracklen,
				 GFextrapolation_decays[2], GFmom_decays[2],
				 GFlpi_dist, GFxi_vert, vtx_scan_range)
       || GFlpi_dist > GFlpi_distcut) continue;

    TLorentzVector GFLpi2(GFmom_decays[2], TMath::Sqrt(GFmom_decays[2].Mag()*GFmom_decays[2].Mag() + PionMass*PionMass));
    TLorentzVector GFLxi = GFLlambda_fixed + GFLpi2;
    TVector3 GFxi_mom = GFlambda_mom + GFmom_decays[2];

    GFTrackCont.AddReconstructedTrack(XiMinusPdgCode, GFxi_vert, GFxi_mom);
    Int_t xi_trackid = GFTrackCont.GetNTrack() - 1;
    GFTrackCont.FitTrack(xi_trackid);

    TVector3 GFxipos_tgtcenter; TVector3 GFximom_tgtcenter;
    Double_t GFxitracklen_tgtcenter; Double_t GFxitof_tgtcenter;
    Bool_t xi_extrapolation =
      GFTrackCont.ExtrapolateToTargetCenter(xi_trackid,
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
    GFxi_gfid_container[candi] = xi_trackid;

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
    Bool_t KFxi_extrapolation =
      GFTrackCont.ExtrapolateToTargetCenter(KFxi_trackid,
					    KFxipos_tgtcenter,
					    KFximom_tgtcenter,
					    KFxitracklen_tgtcenter,
					    KFxitof_tgtcenter);
    KFxi_gfid_container[candi] = KFxi_trackid;
    KFxi_decayvertex_container[candi] = KFxi_vert;

    TVector3 KFdist = tgtpos - KFxipos_tgtcenter;
    KFxi_dist_targetcenter_container[candi] = KFdist.Mag();
    KFxi_pos_targetcenter_container[candi] = KFxipos_tgtcenter;
    KFxi_mom_targetcenter_container[candi] = KFximom_tgtcenter;
    KFxi_gfid_container[candi] = KFxi_trackid;

    KFxi_pull_container[candi] = PullXi;
    KFxi_lpull_container[candi] = PullLd;
    KFxi_l_mom_container0[candi] = KFlambda_mom;
    KFxi_l_mom_container[candi] = KFlambda_mom_KFXi;
    KFxi_lchisqr_container[candi] = KFchisqrl;
    KFxi_lpval_container[candi] = KFpvall;
    KFxi_mom_container[candi] = KFxi_mom;
    KFxi_chisqr_container[candi] = KFchisqrxi;
    KFxi_pval_container[candi] = KFpvalxi;
    KFxi_mass_container[candi] = KFLxi.M();
    KFxi_lpi_closedist_container[candi] = KFlpi_dist;
    KFxi_p_mom_container[candi] = KFTVP;
    KFxi_pi_mom_container[candi] = KFTVPi1;
    KFxi_pi2_mom_container[candi] = KFTVPi2;

    Double_t GFlambda_tof =
      Kinematics::CalcTimeOfFlight(GFlambda_mom.Mag(), GFlambda_tracklen, pdg::LambdaMass());

    Int_t hitid_htof; Double_t tof_htof; Double_t tracklen_htof; TVector3 pos_htof; Double_t track2tgt_dist;
    Bool_t htofextrapolation_p =
      GFTrackCont.TPCHTOFTrackMatching(trackid_p, repid_p, GFlambda_vert,
				       event.HtofSeg, event.posHtof,
				       hitid_htof, tof_htof,
				       tracklen_htof, pos_htof, track2tgt_dist);
    if(htofextrapolation_p){
      GFxi_p_htofid_container[candi] = hitid_htof;
      //GFxi_p_tof_container[candi] = event.tHtof[hitid_htof] - GFlambda_tof - GFxi_tof;
      GFxi_p_tof_container[candi] = event.tHtof[hitid_htof] - GFlambda_tof;
      GFxi_p_tracklen_container[candi] = tracklen_htof;
      GFxi_p_mass2_container[candi] =
	Kinematics::MassSquare(GFmom_decays[0].Mag(),
			       GFxi_p_tracklen_container[candi],
			       GFxi_p_tof_container[candi]);
      GFxi_p_invbeta_container[candi] =
	MathTools::C()*(GFxi_p_tof_container[candi])/tracklen_htof;
    }

    Bool_t htofextrapolation_pi =
      GFTrackCont.TPCHTOFTrackMatching(trackid_pi, repid_pi, GFlambda_vert,
				       event.HtofSeg, event.posHtof,
				       hitid_htof, tof_htof,
				       tracklen_htof, pos_htof, track2tgt_dist);
    if(htofextrapolation_pi){
      GFxi_pi_htofid_container[candi] = hitid_htof;
      //GFxi_pi_tof_container[candi] = event.tHtof[hitid_htof] - GFlambda_tof - GFxi_tof;
      GFxi_pi_tof_container[candi] = event.tHtof[hitid_htof] - GFlambda_tof;
      GFxi_pi_tracklen_container[candi] = tracklen_htof;
      GFxi_pi_mass2_container[candi] =
	Kinematics::MassSquare(GFmom_decays[1].Mag(),
			       GFxi_pi_tracklen_container[candi],
			       GFxi_pi_tof_container[candi]);
      GFxi_pi_invbeta_container[candi] =
	MathTools::C()*(GFxi_pi_tof_container[candi])/tracklen_htof;
    }

    Bool_t htofextrapolation_pi2 =
      GFTrackCont.TPCHTOFTrackMatching(trackid_pi2, repid_pi2, tgtpos,
				       event.HtofSeg, event.posHtof,
				       hitid_htof, tof_htof,
				       tracklen_htof, pos_htof, track2tgt_dist);
    if(htofextrapolation_pi2){
      GFxi_pi2_htofid_container[candi] = hitid_htof;
      //GFxi_pi2_tof_container[candi] = event.tHtof[hitid_htof] - GFxi_tof;
      GFxi_pi2_tof_container[candi] = event.tHtof[hitid_htof];
      GFxi_pi2_tracklen_container[candi] = tracklen_htof;
      GFxi_pi2_mass2_container[candi] =
	Kinematics::MassSquare(GFmom_decays[2].Mag(),
			       GFxi_pi2_tracklen_container[candi],
			       GFxi_pi2_tof_container[candi]);
      GFxi_pi2_invbeta_container[candi] =
	MathTools::C()*(GFxi_pi2_tof_container[candi])/tracklen_htof;
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

    GFxi_l_mom_container[candi] = GFlambda_mom;
    GFxi_l_vert_container[candi] = GFlambda_vert;
    GFxi_l_mass_container[candi] = GFLlambda.M();
    GFxi_l_tracklen_container[candi] = GFlambda_tracklen;
    GFxi_l_tof_container[candi] = GFlambda_tof;
    GFxi_p_mom_container[candi] = GFmom_decays[0];
    GFxi_pi_mom_container[candi] = GFmom_decays[1];
    GFxi_pi2_mom_container[candi] = GFmom_decays[2];
    GFxi_ppi_closedist_container[candi] = GFppi_dist;
    GFxi_lpi_closedist_container[candi] = GFlpi_dist;
    Double_t diff = TMath::Abs(GFLlambda.M() - LambdaMass);
    if(prev_Lmassdiff > diff){
      prev_Lmassdiff = diff;
      best_xi = candi;
    }
  } //candi

#if DebugDisp
  if(xi_candidates > 0) std::cout<<"Fitting processes for Xi- candidates end"<<std::endl;
#endif

  if(event.xiflag){
#if DebugDisp
    std::cout<<"6. Calculate and save the best Xi combination"<<std::endl;
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
    event.lmom = l_mom_container[best_xi].Mag();
    event.lmom_x = l_mom_container[best_xi].x();
    event.lmom_y = l_mom_container[best_xi].y();
    event.lmom_z = l_mom_container[best_xi].z();
    event.ppi_dist = ppi_closedist[best_xi];

    event.xidecays_mom.push_back(xi_p_mom_container[best_xi].Mag());
    event.xidecays_mom.push_back(xi_pi_mom_container[best_xi].Mag());
    event.xidecays_mom.push_back(xi_pi2_mom_container[best_xi].Mag());
    event.xidecays_mom_x.push_back(xi_p_mom_container[best_xi].x());
    event.xidecays_mom_x.push_back(xi_pi_mom_container[best_xi].x());
    event.xidecays_mom_x.push_back(xi_pi2_mom_container[best_xi].x());
    event.xidecays_mom_y.push_back(xi_p_mom_container[best_xi].y());
    event.xidecays_mom_y.push_back(xi_pi_mom_container[best_xi].y());
    event.xidecays_mom_y.push_back(xi_pi2_mom_container[best_xi].y());
    event.xidecays_mom_z.push_back(xi_p_mom_container[best_xi].z());
    event.xidecays_mom_z.push_back(xi_pi_mom_container[best_xi].z());
    event.xidecays_mom_z.push_back(xi_pi2_mom_container[best_xi].z());
    event.xidecays_id.push_back(xi_p_container[best_xi]);
    event.xidecays_id.push_back(xi_pi_container[best_xi]);
    event.xidecays_id.push_back(xi_pi2_container[best_xi]);

    TLorentzVector Lv_p(xi_p_mom_container[best_xi],
			TMath::Hypot(xi_p_mom_container[best_xi].Mag(), ProtonMass));
    TLorentzVector Lv_pi(xi_pi_mom_container[best_xi],
			 TMath::Hypot(xi_pi_mom_container[best_xi].Mag(), PionMass));
    TLorentzVector Lv_pi2(xi_pi2_mom_container[best_xi],
			  TMath::Hypot(xi_pi2_mom_container[best_xi].Mag(), PionMass));
    TLorentzVector Lv_xi_fixedmass(xi_mom_container[best_xi],
				   TMath::Hypot(xi_mom_container[best_xi].Mag(), XiMinusMass));
    TVector3 BoostXi = Lv_xi_fixedmass.BoostVector();
    Lv_p.Boost(-BoostXi);
    Lv_pi.Boost(-BoostXi);
    Lv_pi2.Boost(-BoostXi);
    event.xidecays_CMmom.push_back(Lv_p.P());
    event.xidecays_CMmom_x.push_back(Lv_p.Px());
    event.xidecays_CMmom_y.push_back(Lv_p.Py());
    event.xidecays_CMmom_z.push_back(Lv_p.Pz());
    event.xidecays_CMmom.push_back(Lv_pi.P());
    event.xidecays_CMmom_x.push_back(Lv_pi.Px());
    event.xidecays_CMmom_y.push_back(Lv_pi.Py());
    event.xidecays_CMmom_z.push_back(Lv_pi.Pz());
    event.xidecays_CMmom.push_back(Lv_pi2.P());
    event.xidecays_CMmom_x.push_back(Lv_pi2.Px());
    event.xidecays_CMmom_y.push_back(Lv_pi2.Py());
    event.xidecays_CMmom_z.push_back(Lv_pi2.Pz());

    event.GFximass = GFxi_mass_container[best_xi];
    event.GFxidecayvtx_x = GFxi_decayvertex_container[best_xi].x();
    event.GFxidecayvtx_y = GFxi_decayvertex_container[best_xi].y();
    event.GFxidecayvtx_z = GFxi_decayvertex_container[best_xi].z();
    event.GFximom = GFxi_mom_container[best_xi].Mag();
    event.GFximom_x = GFxi_mom_container[best_xi].x();
    event.GFximom_y = GFxi_mom_container[best_xi].y();
    event.GFximom_z = GFxi_mom_container[best_xi].z();
    event.GFlpi_dist = GFxi_lpi_closedist_container[best_xi];
    event.GFlmass = GFxi_l_mass_container[best_xi];
    event.GFldecayvtx_x = GFxi_l_vert_container[best_xi].x();
    event.GFldecayvtx_y = GFxi_l_vert_container[best_xi].y();
    event.GFldecayvtx_z = GFxi_l_vert_container[best_xi].z();
    event.GFlmom = GFxi_l_mom_container[best_xi].Mag();
    event.GFlmom_x = GFxi_l_mom_container[best_xi].x();
    event.GFlmom_y = GFxi_l_mom_container[best_xi].y();
    event.GFlmom_z = GFxi_l_mom_container[best_xi].z();
    event.GFltracklen = GFxi_l_tracklen_container[best_xi];
    event.GFltof = GFxi_l_tof_container[best_xi];
    event.GFppi_dist = GFxi_ppi_closedist_container[best_xi];

    Int_t GFntdecays = 3;
    event.GFxidecays_nhtrack.resize(GFntdecays);
    event.GFxidecays_chisqr.resize(GFntdecays);
    event.GFxidecays_charge.resize(GFntdecays);
    event.GFxidecays_tracktof.resize(GFntdecays);
    event.GFxidecays_pval.resize(GFntdecays);
    event.GFxidecays_pdgcode.resize(GFntdecays);

    event.GFxidecays_htofid.resize(GFntdecays);
    event.GFxidecays_tracklen.resize(GFntdecays);
    event.GFxidecays_tof.resize(GFntdecays);
    event.GFxidecays_mass2.resize(GFntdecays);
    event.GFxidecays_invbeta.resize(GFntdecays);
    event.GFxidecays_mom.resize(GFntdecays);
    event.GFxidecays_mom_x.resize(GFntdecays);
    event.GFxidecays_mom_y.resize(GFntdecays);
    event.GFxidecays_mom_z.resize(GFntdecays);
    event.GFxidecays_CMmom.resize(GFntdecays);
    event.GFxidecays_CMmom_x.resize(GFntdecays);
    event.GFxidecays_CMmom_y.resize(GFntdecays);
    event.GFxidecays_CMmom_z.resize(GFntdecays);
    event.GFxidecays_momloss.resize(GFntdecays);
    event.GFxidecays_eloss.resize(GFntdecays);

    TLorentzVector GFLv_xi_fixedmass(GFxi_mom_container[best_xi],
				     TMath::Hypot(GFxi_mom_container[best_xi].Mag(), XiMinusMass));
    TVector3 BoostGFXi = GFLv_xi_fixedmass.BoostVector();
    for(Int_t j=0; j<GFntdecays; ++j){
      Int_t igf = GFxi_p_id_container[best_xi];
      if(j==1) igf = GFxi_pi_id_container[best_xi];
      if(j==2) igf = GFxi_pi2_id_container[best_xi];

      Int_t repid = GFxi_p_rep_container[best_xi];
      if(j==1) repid = GFxi_pi_rep_container[best_xi];
      if(j==2) repid = GFxi_pi2_rep_container[best_xi];

      Int_t GFhtofid_decays = GFxi_p_htofid_container[best_xi];
      if(j==1) GFhtofid_decays = GFxi_pi_htofid_container[best_xi];
      if(j==2) GFhtofid_decays = GFxi_pi2_htofid_container[best_xi];

      Double_t GFtof_decays = GFxi_p_tof_container[best_xi];
      if(j==1) GFtof_decays = GFxi_pi_tof_container[best_xi];
      if(j==2) GFtof_decays = GFxi_pi2_tof_container[best_xi];

      Double_t GFtracklen_decays = GFxi_p_tracklen_container[best_xi];
      if(j==1) GFtracklen_decays = GFxi_pi_tracklen_container[best_xi];
      if(j==2) GFtracklen_decays = GFxi_pi2_tracklen_container[best_xi];

      TVector3 GFmom_decays = GFxi_p_mom_container[best_xi];
      if(j==1) GFmom_decays = GFxi_pi_mom_container[best_xi];
      if(j==2) GFmom_decays = GFxi_pi2_mom_container[best_xi];

      Double_t GFmass2_decays = GFxi_p_mass2_container[best_xi];
      if(j==1) GFmass2_decays = GFxi_pi_mass2_container[best_xi];
      if(j==2) GFmass2_decays = GFxi_pi2_mass2_container[best_xi];

      Double_t GFinvbeta_decays = GFxi_p_invbeta_container[best_xi];
      if(j==1) GFinvbeta_decays = GFxi_pi_invbeta_container[best_xi];
      if(j==2) GFinvbeta_decays = GFxi_pi2_invbeta_container[best_xi];

      event.GFxidecays_htofid[j] = GFhtofid_decays;
      event.GFxidecays_tof[j] = GFtof_decays;
      event.GFxidecays_mass2[j] = GFmass2_decays;
      event.GFxidecays_invbeta[j] = GFinvbeta_decays;
      event.GFxidecays_mom[j] = GFmom_decays.Mag();
      event.GFxidecays_mom_x[j] = GFmom_decays.x();
      event.GFxidecays_mom_y[j] = GFmom_decays.y();
      event.GFxidecays_mom_z[j] = GFmom_decays.z();
      event.GFxidecays_momloss[j] = GFmom_decays.Mag() - GFTrackCont.GetMom(igf, 0, repid).Mag();
      event.GFxidecays_eloss[j] = TMath::Hypot(GFmom_decays.Mag(), pdgmass[j]) - TMath::Hypot(GFTrackCont.GetMom(igf, 0, repid).Mag(), pdgmass[j]);

      Double_t mass = PionMass;
      if(j==0) mass = ProtonMass;

      TLorentzVector GFLv_decays(GFmom_decays, TMath::Hypot(GFmom_decays.Mag(), mass));
      GFLv_decays.Boost(-BoostGFXi);
      event.GFxidecays_CMmom[j] = GFLv_decays.P();
      event.GFxidecays_CMmom_x[j] = GFLv_decays.Px();
      event.GFxidecays_CMmom_y[j] = GFLv_decays.Py();
      event.GFxidecays_CMmom_z[j] = GFLv_decays.Pz();

      event.GFxidecays_nhtrack[j] = GFTrackCont.GetNHits(igf);
      event.GFxidecays_chisqr[j] = GFTrackCont.GetChi2NDF(igf, repid);
      event.GFxidecays_charge[j] = GFTrackCont.GetCharge(igf, repid);
      event.GFxidecays_tracktof[j] = GFTrackCont.GetTrackTOF(igf, 0, -1, repid);
      event.GFxidecays_pval[j] = GFTrackCont.GetPvalue(igf, repid);
      event.GFxidecays_pdgcode[j] = GFTrackCont.GetPDGcode(igf, repid);
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
    Bool_t xi_extrapolation_prodvtx =
      GFTrackCont.XiDecayToProdVertex(GFxi_gfid_container[best_xi],
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
    Bool_t xi_extrapolation =
      GFTrackCont.XiDecayToProdVertex(GFxi_gfid_container[best_xi],
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
    xi_extrapolation =
      GFTrackCont.XiDecayToProdVertex(GFxi_gfid_container[best_xi],
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

    event.KFlmom0 = KFxi_l_mom_container0[best_xi].Mag();
    event.KFlmom_x0 = KFxi_l_mom_container0[best_xi].x();
    event.KFlmom_y0 = KFxi_l_mom_container0[best_xi].y();
    event.KFlmom_z0 = KFxi_l_mom_container0[best_xi].z();
    event.KFlmom = KFxi_l_mom_container[best_xi].Mag();
    event.KFlmom_x = KFxi_l_mom_container[best_xi].x();
    event.KFlmom_y = KFxi_l_mom_container[best_xi].y();
    event.KFlmom_z = KFxi_l_mom_container[best_xi].z();
    event.KFlchisqr = KFxi_lchisqr_container[best_xi];
    event.KFlpval = KFxi_lpval_container[best_xi];
    event.KFlpull = KFxi_lpull_container[best_xi];
    event.KFlpi_dist = KFxi_lpi_closedist_container[best_xi];
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
    event.KFxipull = KFxi_pull_container[best_xi];
    event.KFxidecays_mom.push_back(KFxi_p_mom_container[best_xi].Mag());
    event.KFxidecays_mom_x.push_back(KFxi_p_mom_container[best_xi].x());
    event.KFxidecays_mom_y.push_back(KFxi_p_mom_container[best_xi].y());
    event.KFxidecays_mom_z.push_back(KFxi_p_mom_container[best_xi].z());
    event.KFxidecays_mom.push_back(KFxi_pi_mom_container[best_xi].Mag());
    event.KFxidecays_mom_x.push_back(KFxi_pi_mom_container[best_xi].x());
    event.KFxidecays_mom_y.push_back(KFxi_pi_mom_container[best_xi].y());
    event.KFxidecays_mom_z.push_back(KFxi_pi_mom_container[best_xi].z());
    event.KFxidecays_mom.push_back(KFxi_pi2_mom_container[best_xi].Mag());
    event.KFxidecays_mom_x.push_back(KFxi_pi2_mom_container[best_xi].x());
    event.KFxidecays_mom_y.push_back(KFxi_pi2_mom_container[best_xi].y());
    event.KFxidecays_mom_z.push_back(KFxi_pi2_mom_container[best_xi].z());

    TLorentzVector KFLv_p(KFxi_p_mom_container[best_xi],
			  TMath::Hypot(KFxi_p_mom_container[best_xi].Mag(), ProtonMass));
    TLorentzVector KFLv_pi(KFxi_pi_mom_container[best_xi],
			   TMath::Hypot(KFxi_pi_mom_container[best_xi].Mag(), PionMass));
    TLorentzVector KFLv_pi2(KFxi_pi2_mom_container[best_xi],
			    TMath::Hypot(KFxi_pi2_mom_container[best_xi].Mag(), PionMass));
    TLorentzVector KFLv_Xi(KFxi_mom_container[best_xi],
			   TMath::Hypot(KFxi_mom_container[best_xi].Mag(), XiMinusMass));
    auto BoostKFXi = KFLv_Xi.BoostVector();
    KFLv_p.Boost(-BoostKFXi);
    KFLv_pi.Boost(-BoostKFXi);
    KFLv_pi2.Boost(-BoostKFXi);
    event.KFxidecays_CMmom.push_back(KFLv_p.P());
    event.KFxidecays_CMmom_x.push_back(KFLv_p.Px());
    event.KFxidecays_CMmom_y.push_back(KFLv_p.Py());
    event.KFxidecays_CMmom_z.push_back(KFLv_p.Pz());
    event.KFxidecays_CMmom.push_back(KFLv_pi.P());
    event.KFxidecays_CMmom_x.push_back(KFLv_pi.Px());
    event.KFxidecays_CMmom_y.push_back(KFLv_pi.Py());
    event.KFxidecays_CMmom_z.push_back(KFLv_pi.Pz());
    event.KFxidecays_CMmom.push_back(KFLv_pi2.P());
    event.KFxidecays_CMmom_x.push_back(KFLv_pi2.Px());
    event.KFxidecays_CMmom_y.push_back(KFLv_pi2.Py());
    event.KFxidecays_CMmom_z.push_back(KFLv_pi2.Pz());

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
    Bool_t KFxi_extrapolation_kkvtx =
      GFTrackCont.XiDecayToProdVertex(KFxi_gfid_container[best_xi],
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
    Bool_t KFxi_extrapolation_tgt =
      GFTrackCont.XiDecayToProdVertex(KFxi_gfid_container[best_xi],
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
    TVector3 KFkkxi_prodvertex =
      Kinematics::MultitrackVertex(ntrack_xi,
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
    Bool_t KFxi_extrapolation_prodvtx =
      GFTrackCont.XiDecayToProdVertex(KFxi_gfid_container[best_xi],
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

    TVector3 KFkpxi_prodvtx =
      Kinematics::MultitrackVertex(2, KFxtgtkpxi, KFytgtkpxi,
				   KFutgtkpxi, KFvtgtkpxi,
				   KFkpxi_res_x, KFkpxi_res_y,
				   KFkpxi_res_u, KFkpxi_res_v);
    event.KFprodvtx_x_kpxi = KFkpxi_prodvtx.x();
    event.KFprodvtx_y_kpxi = KFkpxi_prodvtx.y();
    event.KFprodvtx_z_kpxi = KFkpxi_prodvtx.z();

    TVector3 KFkpxipos_prodvtx; TVector3 KFkpximom_prodvtx;
    Double_t KFkpxitracklen_prodvtx; Double_t KFkpxitof_prodvtx;
    Bool_t KFkpxi_extrapolation_prodvtx =
      GFTrackCont.XiDecayToProdVertex(KFxi_gfid_container[best_xi],
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

    //Remaining tracks
#if DebugDisp
    std::cout<<"Xi event, Save remaining tracks"<<std::endl;
#endif

    event.xip_multi = target_p_id_container.size();
    event.xipip_multi = target_pip_id_container.size();
    event.xipim_multi = target_pim_id_container.size();
    event.xiep_multi = target_ep_id_container.size();
    event.xiem_multi = target_em_id_container.size();
    event.xippip_multi = target_ppip_id_container.size();
    TVector3 closepoint, closepoint_mom;
    Double_t closepoint_tracklen, closepoint_tof;
    for(Int_t itp=0; itp<target_p_id_container.size(); ++itp){
      Int_t id_p = target_p_id_container[itp];
      if(id_p==event.xidecays_id[0]){
	event.xip_multi -= 1;
	continue;
      }
      if(target_p_dist2tgt_container[itp] > residual_track_distcut){
	event.xip_multi -= 1;
	continue;
      }
      event.xiresidual_id.push_back(id_p);
      event.xiresidual_dist2tgt.push_back(target_p_dist2tgt_container[itp]);
      event.xiresidual_mass2.push_back(target_p_mass2_container[itp]);
      event.xiresidual_invbeta.push_back(target_p_invbeta_container[itp]);
      event.xiresidual_mom.push_back(target_p_mom_container[itp].Mag());
      event.xiresidual_mom_x.push_back(target_p_mom_container[itp].x());
      event.xiresidual_mom_y.push_back(target_p_mom_container[itp].y());
      event.xiresidual_mom_z.push_back(target_p_mom_container[itp].z());
      event.xiresidual_charge.push_back(1);

      if(GFTrackCont.ExtrapolateToPoint(id_p,
					kk_xi_vertex,
					closepoint,
					closepoint_mom,
					closepoint_tracklen,
					closepoint_tof,
					target_p_repid_container[itp])){
	closepoint -= kk_xi_vertex;
	event.xiresidual_GFdist2prodvtx.push_back(closepoint.Mag());
      }
      else event.xiresidual_GFdist2prodvtx.push_back(qnan);
      if(GFTrackCont.ExtrapolateToPoint(id_p,
					KFkkxi_prodvertex,
					closepoint,
					closepoint_mom,
					closepoint_tracklen,
					closepoint_tof,
					target_p_repid_container[itp])){
	closepoint -= KFkkxi_prodvertex;
	event.xiresidual_KFdist2prodvtx.push_back(closepoint.Mag());
      }
      else event.xiresidual_KFdist2prodvtx.push_back(qnan);
    }
    for(Int_t itpip=0; itpip<target_pip_id_container.size(); ++itpip){
      Int_t id_pip = target_pip_id_container[itpip];
      if(id_pip==event.xidecays_id[0]){
	event.xipip_multi -= 1;
	continue;
      }
      if(target_pip_dist2tgt_container[itpip] > residual_track_distcut){
	event.xipip_multi -= 1;
	continue;
      }
      event.xiresidual_id.push_back(id_pip);
      event.xiresidual_dist2tgt.push_back(target_pip_dist2tgt_container[itpip]);
      event.xiresidual_mass2.push_back(target_pip_mass2_container[itpip]);
      event.xiresidual_invbeta.push_back(target_pip_invbeta_container[itpip]);
      event.xiresidual_mom.push_back(target_pip_mom_container[itpip].Mag());
      event.xiresidual_mom_x.push_back(target_pip_mom_container[itpip].x());
      event.xiresidual_mom_y.push_back(target_pip_mom_container[itpip].y());
      event.xiresidual_mom_z.push_back(target_pip_mom_container[itpip].z());
      event.xiresidual_charge.push_back(1);

      if(GFTrackCont.ExtrapolateToPoint(id_pip,
					kk_xi_vertex,
					closepoint,
					closepoint_mom,
					closepoint_tracklen,
					closepoint_tof,
					target_pip_repid_container[itpip])){
	closepoint -= kk_xi_vertex;
	event.xiresidual_GFdist2prodvtx.push_back(closepoint.Mag());
      }
      else event.xiresidual_GFdist2prodvtx.push_back(qnan);
      if(GFTrackCont.ExtrapolateToPoint(id_pip,
					KFkkxi_prodvertex,
					closepoint,
					closepoint_mom,
					closepoint_tracklen,
					closepoint_tof,
					target_pip_repid_container[itpip])){
	closepoint -= KFkkxi_prodvertex;
	event.xiresidual_KFdist2prodvtx.push_back(closepoint.Mag());
      }
      else event.xiresidual_KFdist2prodvtx.push_back(qnan);
    }
    for(Int_t itpim=0; itpim<target_pim_id_container.size(); ++itpim){
      Int_t id_pim = target_pim_id_container[itpim];
      if(id_pim==event.xidecays_id[1] || id_pim==event.xidecays_id[2]){
	event.xipim_multi -= 1;
	continue;
      }
      if(target_pim_dist2tgt_container[itpim] > residual_track_distcut){
	event.xipim_multi -= 1;
	continue;
      }
      event.xiresidual_id.push_back(id_pim);
      event.xiresidual_dist2tgt.push_back(target_pim_dist2tgt_container[itpim]);
      event.xiresidual_mass2.push_back(target_pim_mass2_container[itpim]);
      event.xiresidual_invbeta.push_back(target_pim_invbeta_container[itpim]);
      event.xiresidual_mom.push_back(target_pim_mom_container[itpim].Mag());
      event.xiresidual_mom_x.push_back(target_pim_mom_container[itpim].x());
      event.xiresidual_mom_y.push_back(target_pim_mom_container[itpim].y());
      event.xiresidual_mom_z.push_back(target_pim_mom_container[itpim].z());
      event.xiresidual_charge.push_back(-1);

      if(GFTrackCont.ExtrapolateToPoint(id_pim,
					kk_xi_vertex,
					closepoint,
					closepoint_mom,
					closepoint_tracklen,
					closepoint_tof,
					target_pim_repid_container[itpim])){
	closepoint -= kk_xi_vertex;
	event.xiresidual_GFdist2prodvtx.push_back(closepoint.Mag());
      }
      else event.xiresidual_GFdist2prodvtx.push_back(qnan);
      if(GFTrackCont.ExtrapolateToPoint(id_pim,
					KFkkxi_prodvertex,
					closepoint,
					closepoint_mom,
					closepoint_tracklen,
					closepoint_tof,
					target_pim_repid_container[itpim])){
	closepoint -= KFkkxi_prodvertex;
	event.xiresidual_KFdist2prodvtx.push_back(closepoint.Mag());
      }
      else event.xiresidual_KFdist2prodvtx.push_back(qnan);
    }
    for(Int_t itppip=0; itppip<target_ppip_id_container.size(); ++itppip){
      Int_t id_ppip = target_ppip_id_container[itppip];
      if(id_ppip==event.xidecays_id[0]){
	event.xippip_multi -= 1;
	continue;
      }
      if(target_ppip_dist2tgt_container[itppip] > residual_track_distcut){
	event.xippip_multi -= 1;
	continue;
      }
      event.xiresidual_id.push_back(id_ppip);
      event.xiresidual_dist2tgt.push_back(target_ppip_dist2tgt_container[itppip]);
      event.xiresidual_mass2.push_back(target_ppip_mass2_container[itppip]);
      event.xiresidual_invbeta.push_back(target_ppip_invbeta_container[itppip]);
      event.xiresidual_mom.push_back(target_ppip_mom_container[itppip].Mag());
      event.xiresidual_mom_x.push_back(target_ppip_mom_container[itppip].x());
      event.xiresidual_mom_y.push_back(target_ppip_mom_container[itppip].y());
      event.xiresidual_mom_z.push_back(target_ppip_mom_container[itppip].z());
      event.xiresidual_charge.push_back(1);

      if(GFTrackCont.ExtrapolateToPoint(id_ppip,
					kk_xi_vertex,
					closepoint,
					closepoint_mom,
					closepoint_tracklen,
					closepoint_tof,
					target_ppip_repid_container[itppip])){
	closepoint -= kk_xi_vertex;
	event.xiresidual_GFdist2prodvtx.push_back(closepoint.Mag());
      }
      else event.xiresidual_GFdist2prodvtx.push_back(qnan);
      if(GFTrackCont.ExtrapolateToPoint(id_ppip,
					KFkkxi_prodvertex,
					closepoint,
					closepoint_mom,
					closepoint_tracklen,
					closepoint_tof,
					target_ppip_repid_container[itppip])){
	closepoint -= KFkkxi_prodvertex;
	event.xiresidual_KFdist2prodvtx.push_back(closepoint.Mag());
      }
      else event.xiresidual_KFdist2prodvtx.push_back(qnan);
    }
    for(Int_t itep=0; itep<target_ep_id_container.size(); ++itep){
      Int_t id_ep = target_ep_id_container[itep];
      if(id_ep==event.xidecays_id[0]){
	event.xiep_multi -= 1;
	continue;
      }
      if(target_ep_dist2tgt_container[itep] > residual_track_distcut){
	event.xiep_multi -= 1;
	continue;
      }
      event.xiresidual_id.push_back(id_ep);
      event.xiresidual_dist2tgt.push_back(target_ep_dist2tgt_container[itep]);
      event.xiresidual_mass2.push_back(target_ep_mass2_container[itep]);
      event.xiresidual_invbeta.push_back(target_ep_invbeta_container[itep]);
      event.xiresidual_mom.push_back(target_ep_mom_container[itep].Mag());
      event.xiresidual_mom_x.push_back(target_ep_mom_container[itep].x());
      event.xiresidual_mom_y.push_back(target_ep_mom_container[itep].y());
      event.xiresidual_mom_z.push_back(target_ep_mom_container[itep].z());
      event.xiresidual_charge.push_back(1);

      if(GFTrackCont.ExtrapolateToPoint(id_ep,
					kk_xi_vertex,
					closepoint,
					closepoint_mom,
					closepoint_tracklen,
					closepoint_tof,
					target_ep_repid_container[itep])){
	closepoint -= kk_xi_vertex;
	event.xiresidual_GFdist2prodvtx.push_back(closepoint.Mag());
      }
      else event.xiresidual_GFdist2prodvtx.push_back(qnan);
      if(GFTrackCont.ExtrapolateToPoint(id_ep,
					KFkkxi_prodvertex,
					closepoint,
					closepoint_mom,
					closepoint_tracklen,
					closepoint_tof,
					target_ep_repid_container[itep])){
	closepoint -= KFkkxi_prodvertex;
	event.xiresidual_KFdist2prodvtx.push_back(closepoint.Mag());
      }
      else event.xiresidual_KFdist2prodvtx.push_back(qnan);
    }
    for(Int_t item=0; item<target_em_id_container.size(); ++item){
      Int_t id_em = target_em_id_container[item];
      if(id_em==event.xidecays_id[1] || id_em==event.xidecays_id[2]){
	event.xiem_multi -= 1;
	continue;
      }
      if(target_em_dist2tgt_container[item] > residual_track_distcut){
	event.xiem_multi -= 1;
	continue;
      }
      event.xiresidual_id.push_back(id_em);
      event.xiresidual_dist2tgt.push_back(target_em_dist2tgt_container[item]);
      event.xiresidual_mass2.push_back(target_em_mass2_container[item]);
      event.xiresidual_invbeta.push_back(target_em_invbeta_container[item]);
      event.xiresidual_mom.push_back(target_em_mom_container[item].Mag());
      event.xiresidual_mom_x.push_back(target_em_mom_container[item].x());
      event.xiresidual_mom_y.push_back(target_em_mom_container[item].y());
      event.xiresidual_mom_z.push_back(target_em_mom_container[item].z());
      event.xiresidual_charge.push_back(-1);

      if(GFTrackCont.ExtrapolateToPoint(id_em,
					kk_xi_vertex,
					closepoint,
					closepoint_mom,
					closepoint_tracklen,
					closepoint_tof,
					target_em_repid_container[item])){
	closepoint -= kk_xi_vertex;
	event.xiresidual_GFdist2prodvtx.push_back(closepoint.Mag());
      }
      else event.xiresidual_GFdist2prodvtx.push_back(qnan);
      if(GFTrackCont.ExtrapolateToPoint(id_em,
					KFkkxi_prodvertex,
					closepoint,
					closepoint_mom,
					closepoint_tracklen,
					closepoint_tof,
					target_em_repid_container[item])){
	closepoint -= KFkkxi_prodvertex;
	event.xiresidual_KFdist2prodvtx.push_back(closepoint.Mag());
      }
      else event.xiresidual_KFdist2prodvtx.push_back(qnan);
    }

    event.xiresidual_multi =
      event.xiep_multi + event.xiep_multi +
      event.xipip_multi + event.xip_multi +
      event.xipim_multi + event.xippip_multi;

#if DebugDisp
    std::cout<<"Xi event, Saving gamma candidates"<<std::endl;
#endif

    event.xig_multi = gamma_mom_container.size();
    for(Int_t itg=0; itg<gamma_mom_container.size(); ++itg){
      Int_t epid = gamma_ep_id_container[itg];
      Int_t emid = gamma_em_id_container[itg];
      if(epid==event.xidecays_id[0] ||
	 emid==event.xidecays_id[1] ||
	 emid==event.xidecays_id[2]){
	event.xig_multi -= 1;
	continue;
      }
      auto iter_ep = find(event.xiepidgamma.begin(), event.xiepidgamma.end(), epid);
      if(iter_ep != event.xiepidgamma.end()){
	event.xig_multi -=1;
	continue;
      }
      auto iter_em = find(event.xiemidgamma.begin(), event.xiemidgamma.end(), emid);
      if(iter_em != event.xiemidgamma.end()){
	event.xig_multi -=1;
	continue;
      }

      Double_t dist = gamma_epemdist_container[itg];
      TVector3 ep_mom = gamma_ep_mom_container[itg];
      TVector3 em_mom = gamma_em_mom_container[itg];
      TVector3 gamma_mom = gamma_mom_container[itg];
      TVector3 gamma_vtx = gamma_decayvtx_container[itg];

      event.xiepidgamma.push_back(epid);
      event.xiemidgamma.push_back(emid);
      event.xiepmomgamma.push_back(ep_mom.Mag());
      event.xiepmomgamma_x.push_back(ep_mom.x());
      event.xiepmomgamma_y.push_back(ep_mom.y());
      event.xiepmomgamma_z.push_back(ep_mom.z());
      event.xiemmomgamma.push_back(em_mom.Mag());
      event.xiemmomgamma_x.push_back(em_mom.x());
      event.xiemmomgamma_y.push_back(em_mom.y());
      event.xiemmomgamma_z.push_back(em_mom.z());
      event.ximomgamma.push_back(gamma_mom.Mag());
      event.ximomgamma_x.push_back(gamma_mom.x());
      event.ximomgamma_y.push_back(gamma_mom.y());
      event.ximomgamma_z.push_back(gamma_mom.z());
      event.xiepidistgamma.push_back(dist);
      event.xivtxgamma_x.push_back(gamma_vtx.x());
      event.xivtxgamma_y.push_back(gamma_vtx.y());
      event.xivtxgamma_z.push_back(gamma_vtx.z());
    }

    //find Xi-p event
    if(event.xiresidual_multi==1 && event.xip_multi==1) event.xipflag = true; //Xi event
  } //if(event.xiflag)

#if DebugDisp
  std::cout<<"7. Single L searching starts (No Xi-, LL)"<<std::endl;
#endif
  Int_t GFtrackid_decays[2] = {-1, -1}; Int_t GFrepid_decays[2];

  if(!event.xiflag && !event.llflag){
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

      TLorentzVector Lv_L_fixedmass(L_mom_container[id],
				    TMath::Hypot(L_mom_container[id].Mag(), LambdaMass));
      auto BoostL = Lv_L_fixedmass.BoostVector();

      TLorentzVector Lv_p(L_p_mom_container[id],
			  TMath::Hypot(L_p_mom_container[id].Mag(), ProtonMass));
      Lv_p.Boost(-BoostL);
      event.decays_CMmom.push_back(Lv_p.P());
      event.decays_CMmom_x.push_back(Lv_p.Px());
      event.decays_CMmom_y.push_back(Lv_p.Py());
      event.decays_CMmom_z.push_back(Lv_p.Pz());

      TLorentzVector Lv_pi(L_pi_mom_container[id],
			   TMath::Hypot(L_pi_mom_container[id].Mag(), PionMass));
      Lv_pi.Boost(-BoostL);
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
      event.GFdecays_tof.push_back(GFL_p_tof_container[id]);
      event.GFdecays_tof.push_back(GFL_pi_tof_container[id]);
      event.GFdecays_mass2.push_back(GFL_p_mass2_container[id]);
      event.GFdecays_mass2.push_back(GFL_pi_mass2_container[id]);
      event.GFdecays_invbeta.push_back(GFL_p_invbeta_container[id]);
      event.GFdecays_invbeta.push_back(GFL_pi_invbeta_container[id]);

      event.GFdecays_mom.push_back(GFL_p_mom_container[id].Mag());
      event.GFdecays_mom.push_back(GFL_pi_mom_container[id].Mag());
      event.GFdecays_mom_x.push_back(GFL_p_mom_container[id].x());
      event.GFdecays_mom_x.push_back(GFL_pi_mom_container[id].x());
      event.GFdecays_mom_y.push_back(GFL_p_mom_container[id].y());
      event.GFdecays_mom_y.push_back(GFL_pi_mom_container[id].y());
      event.GFdecays_mom_z.push_back(GFL_p_mom_container[id].z());
      event.GFdecays_mom_z.push_back(GFL_pi_mom_container[id].z());
      event.GFdecays_momloss.push_back(qnan);
      event.GFdecays_momloss.push_back(qnan);
      event.GFdecays_eloss.push_back(qnan);
      event.GFdecays_eloss.push_back(qnan);

      TLorentzVector GFLv_L_fixedmass(GFL_mom_container[id],
				      TMath::Hypot(GFL_mom_container[id].Mag(), LambdaMass));
      auto BoostGFL = GFLv_L_fixedmass.BoostVector();
      TLorentzVector GFLv_p(GFL_p_mom_container[id],
			    TMath::Hypot(GFL_p_mom_container[id].Mag(), ProtonMass));
      GFLv_p.Boost(-BoostGFL);
      event.GFdecays_CMmom.push_back(GFLv_p.P());
      event.GFdecays_CMmom_x.push_back(GFLv_p.Px());
      event.GFdecays_CMmom_y.push_back(GFLv_p.Py());
      event.GFdecays_CMmom_z.push_back(GFLv_p.Pz());

      TLorentzVector GFLv_pi(GFL_pi_mom_container[id],
			     TMath::Hypot(GFL_pi_mom_container[id].Mag(), PionMass));
      GFLv_pi.Boost(-BoostGFL);
      event.GFdecays_CMmom.push_back(GFLv_pi.P());
      event.GFdecays_CMmom_x.push_back(GFLv_pi.Px());
      event.GFdecays_CMmom_y.push_back(GFLv_pi.Py());
      event.GFdecays_CMmom_z.push_back(GFLv_pi.Pz());

      int trackid_p1 = GFL_p_id_container[id];
      int trackid_pi1 = GFL_pi_id_container[id];
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
      auto KFLVLd1 = TLorentzVector(KFTVLd1,
				    TMath::Hypot(KFTVLd1.Mag(), pdg::LambdaMass()));
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

      TVector3 BoostKFLd = KFLVLd1.BoostVector();
      KFLv_p1.Boost(-BoostKFLd);
      KFLv_pi1.Boost(-BoostKFLd);
      event.KFdecays_CMmom.push_back(KFLv_p1.P());
      event.KFdecays_CMmom_x.push_back(KFLv_p1.Px());
      event.KFdecays_CMmom_y.push_back(KFLv_p1.Py());
      event.KFdecays_CMmom_z.push_back(KFLv_p1.Pz());
      event.KFdecays_CMmom.push_back(KFLv_pi1.P());
      event.KFdecays_CMmom_x.push_back(KFLv_pi1.Px());
      event.KFdecays_CMmom_y.push_back(KFLv_pi1.Py());
      event.KFdecays_CMmom_z.push_back(KFLv_pi1.Pz());

      //order : p, pi
      Int_t GFntdecays = 2;
      event.GFdecays_nhtrack.resize(GFntdecays);
      event.GFdecays_chisqr.resize(GFntdecays);
      event.GFdecays_charge.resize(GFntdecays);
      event.GFdecays_tof.resize(GFntdecays);
      event.GFdecays_tracklen.resize(GFntdecays);
      event.GFdecays_pval.resize(GFntdecays);
      event.GFdecays_pdgcode.resize(GFntdecays);

      for(int j=0;j<GFntdecays;j++){
	Int_t igf = GFtrackid_decays[j];
	Int_t repid = GFrepid_decays[j];
	event.GFdecays_nhtrack[j] = GFTrackCont.GetNHits(igf);
	event.GFdecays_chisqr[j] = GFTrackCont.GetChi2NDF(igf, repid);
	event.GFdecays_charge[j] = GFTrackCont.GetCharge(igf, repid);
	event.GFdecays_tof[j] = GFTrackCont.GetTrackTOF(igf, 0, -1, repid);
	event.GFdecays_tracklen[j] = GFTrackCont.GetTrackLength(igf, 0, -1, repid);
	event.GFdecays_pval[j] = GFTrackCont.GetPvalue(igf, repid);
	event.GFdecays_pdgcode[j] = GFTrackCont.GetPDGcode(igf, repid);
      } //GFntdecays

      std::vector<Int_t> p_id_container;
      p_id_container.assign(target_p_id_container.begin(), target_p_id_container.end());
      {
	auto iter = find(p_id_container.begin(), p_id_container.end(), GFtrackid_decays[0]);
	if(iter != p_id_container.end()){
	  p_id_container.erase(iter);
	}
      }
      std::vector<Int_t> ppip_id_container;
      ppip_id_container.assign(target_ppip_id_container.begin(), target_ppip_id_container.end());
      {
	auto iter = find(ppip_id_container.begin(), ppip_id_container.end(), GFtrackid_decays[0]);
	if(iter != ppip_id_container.end()){
	  ppip_id_container.erase(iter);
	}
      }
      std::vector<Int_t> pip_id_container;
      pip_id_container.assign(target_pip_id_container.begin(), target_pip_id_container.end());
      {
	auto iter = find(pip_id_container.begin(), pip_id_container.end(), GFtrackid_decays[0]);
	if(iter != pip_id_container.end()){
	  pip_id_container.erase(iter);
	}
      }

      std::vector<Int_t> km_id_container;
      std::vector<Int_t> km_repid_container;
      std::vector<TVector3> km_mom_container;
      std::vector<Double_t> km_mass2_container;
      std::vector<Double_t> km_invbeta_container;
      for(Int_t itkm=0; itkm<target_km_id_container.size(); ++itkm){
	Int_t id_km = target_km_id_container[itkm];
	if(id_km != GFtrackid_decays[1]){
	  km_id_container.push_back(target_km_id_container[itkm]);
	  km_repid_container.push_back(target_km_repid_container[itkm]);
	  km_mom_container.push_back(target_km_mom_container[itkm]);
	  km_mass2_container.push_back(target_km_mass2_container[itkm]);
	  km_invbeta_container.push_back(target_km_invbeta_container[itkm]);
	}
      }

      //Lphi searching
      if(event.MissMassCorrDETPC[0] > 1.5 &&
	 km_id_container.size()==1 &&
	 km_mom_container[0].Mag() < 1.0 &&
	 p_id_container.size()==0 &&
	 pip_id_container.size()==0 &&
	 ppip_id_container.size()==0){

#if DebugDisp
	std::cout<<"L phi searching "<<std::endl;
#endif

	Int_t id_km = km_id_container[0];
	Int_t km_repid = km_repid_container[0];
	Double_t km_mass2 = km_mass2_container[0];
	Double_t km_invbeta = km_invbeta_container[0];

	std::vector<Int_t> pim_id_container;
	pim_id_container.assign(target_pim_id_container.begin(), target_pim_id_container.end());
	{
	  auto iter = find(pim_id_container.begin(), pim_id_container.end(), GFtrackid_decays[1]);
	  if(iter != pim_id_container.end()){
	    pim_id_container.erase(iter);
	  }
	}
	{
	  auto iter = find(pim_id_container.begin(), pim_id_container.end(), id_km);
	  if(iter != pim_id_container.end()){
	    pim_id_container.erase(iter);
	  }
	}

	if(pim_id_container.size()==0){
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
	    Int_t kp_repid = -1; //K+'s pdgcode
	    Double_t kp_mass2 = target_kurama_mass2_container[0];
	    Double_t kp_invbeta = target_kurama_invbeta_container[0];

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
	    TVector3 phi_vert =
	      Kinematics::LambdaVertex(dMagneticField,
				       kp_par, km_par,
				       kp_theta_min, kp_theta_max,
				       km_theta_min, km_theta_max,
				       kp_mom, km_mom,
				       phi_mom, kk_dist);
	    phi_mom = km_mom + kp_mom;

	    TLorentzVector Lkp(kp_mom, TMath::Hypot(kp_mom.Mag(), KaonMass));
	    TLorentzVector Lkm(km_mom, TMath::Hypot(km_mom.Mag(), KaonMass));
	    TLorentzVector Lphi = Lkp + Lkm;
	    Double_t cosKK = TMath::Cos(km_mom.Angle(kp_mom));
	    if(!TMath::IsNaN(kk_dist) &&
	       kk_dist < phi_kk_distcut&&
	       TMath::Abs(phi_vert.x()) < 30. &&
	       TMath::Abs(phi_vert.z() - tpc::ZTarget) < 30. &&
	       TMath::Abs(phi_vert.y()) < 30. &&
	       GFTrackCont.TrackCheck(id_kp, kp_repid) &&
	       GFTrackCont.TrackCheck(id_km, km_repid)){

	      Double_t GFkk_dist = 10000.;
	      Double_t GFextrapolation_kk[2];
	      TVector3 GFmom_kk[2];
	      TVector3 GFphi_vert;
	      GFTrackCont.FindVertex(id_kp, id_km,
				     kp_repid, km_repid,
				     GFextrapolation_kk[0],
				     GFextrapolation_kk[1],
				     GFmom_kk[0], GFmom_kk[1],
				     GFkk_dist, GFphi_vert,
				     vtx_scan_range);

	      Double_t GFcosKK = TMath::Cos(GFmom_kk[1].Angle(GFmom_kk[0]));
	      TLorentzVector GFkm(GFmom_kk[0],
				  TMath::Hypot(GFmom_kk[0].Mag(), KaonMass));
	      TLorentzVector GFkp(GFmom_kk[1],
				  TMath::Hypot(GFmom_kk[1].Mag(), KaonMass));
	      TVector3 GFmom_phi = GFmom_kk[0] + GFmom_kk[1];
	      TLorentzVector LvPhi = GFkm + GFkp;
	      TLorentzVector LvPhi_fixedmass(GFmom_phi,
					     TMath::Hypot(GFmom_phi.Mag(), PhiMass));
	      TLorentzVector LvLambda(GFL_mom_container[id],
				      TMath::Hypot(GFL_mom_container[id].Mag(), LambdaMass));
	      TLorentzVector LvRc_phiLambda = LvKmTPC + LvP - LvPhi_fixedmass - LvLambda;

	      //Kp momentum from TPCKurama tracking
	      TVector3 KuramaKp = GFmom_kk[1].Unit();
	      KuramaKp *= event.pCorrDETPC[0];
	      TLorentzVector GFkp_wKurama(KuramaKp,
					  TMath::Hypot(KuramaKp.Mag(), KaonMass));
	      TVector3 GFmom_phi_wKurama = GFmom_kk[0] + KuramaKp;
	      TLorentzVector LvPhi_wKurama = GFkm + GFkp_wKurama;

	      if(GFkk_dist < GFphi_kk_distcut && GFcosKK > 0.8 && LvRc_phiLambda.Mag() < 0){
		TVector3 phivtx_dist(event.GFprodvtx_x_l - GFphi_vert.x(),
				     event.GFprodvtx_y_l - GFphi_vert.y(),
				     event.GFprodvtx_z_l - GFphi_vert.z());

		event.lphiflag = true;
		event.GFphimass = LvPhi.M();
		event.GFphidecayvtx_x = GFphi_vert.x();
		event.GFphidecayvtx_y = GFphi_vert.y();
		event.GFphidecayvtx_z = GFphi_vert.z();
		event.GFphicosKK = GFcosKK;
		event.GFphimom = GFmom_phi.Mag();
		event.GFphimom_x = GFmom_phi.x();
		event.GFphimom_y = GFmom_phi.y();
		event.GFphimom_z = GFmom_phi.z();
		event.GFkk_dist = GFkk_dist;
		event.GFphiprodvtx_dist = phivtx_dist.Mag();

		event.GFphimass_wKurama = LvPhi_wKurama.M();
		event.GFphimom_wKurama = GFmom_phi_wKurama.Mag();
		event.GFphimom_x_wKurama = GFmom_phi_wKurama.x();
		event.GFphimom_y_wKurama = GFmom_phi_wKurama.y();
		event.GFphimom_z_wKurama = GFmom_phi_wKurama.z();

		event.phimass = Lphi.M();
		event.phidecayvtx_x = phi_vert.x();
		event.phidecayvtx_y = phi_vert.y();
		event.phidecayvtx_z = phi_vert.z();
		event.phicosKK = cosKK;
		event.phimom = phi_mom.Mag();
		event.phimom_x = phi_mom.x();
		event.phimom_y = phi_mom.y();
		event.phimom_z = phi_mom.z();
		event.kk_dist = kk_dist;
		event.GFphi_km_mass2 = km_mass2;
		event.GFphi_km_invbeta = km_invbeta;
		event.GFphi_kp_mass2 = kp_mass2;
		event.GFphi_kp_invbeta = kp_invbeta;

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

		event.GFphidecays_mom.push_back(GFmom_kk[0].Mag());
		event.GFphidecays_mom.push_back(GFmom_kk[1].Mag());
		event.GFphidecays_mom_x.push_back(GFmom_kk[0].x());
		event.GFphidecays_mom_x.push_back(GFmom_kk[1].x());
		event.GFphidecays_mom_y.push_back(GFmom_kk[0].y());
		event.GFphidecays_mom_y.push_back(GFmom_kk[1].y());
		event.GFphidecays_mom_z.push_back(GFmom_kk[0].z());
		event.GFphidecays_mom_z.push_back(GFmom_kk[1].z());
	      }
	    }
	  }
	}
      } //lphiflag
    } //lflag

    //Remaining tracks
#if DebugDisp
    std::cout<<"One L event or others, Save remaining tracks"<<std::endl;
#endif

    std::vector<Double_t> twopi;
    event.p_multi = target_p_id_container.size();
    event.pip_multi = target_pip_id_container.size();
    event.pim_multi = target_pim_id_container.size();
    event.ep_multi = target_ep_id_container.size();
    event.em_multi = target_em_id_container.size();
    TVector3 closepoint, closepoint_mom;
    Double_t closepoint_tracklen, closepoint_tof;
    event.ppip_multi = target_ppip_id_container.size();
    for(Int_t itp=0; itp<target_p_id_container.size(); ++itp){
      Int_t id_p = target_p_id_container[itp];
      if(id_p==GFtrackid_decays[0]){
	event.p_multi -= 1;
	continue;
      }
      event.residual_id.push_back(id_p);
      event.residual_dist2tgt.push_back(target_p_dist2tgt_container[itp]);
      event.residual_mass2.push_back(target_p_mass2_container[itp]);
      event.residual_invbeta.push_back(target_p_invbeta_container[itp]);
      event.residual_mom.push_back(target_p_mom_container[itp].Mag());
      event.residual_mom_x.push_back(target_p_mom_container[itp].x());
      event.residual_mom_y.push_back(target_p_mom_container[itp].y());
      event.residual_mom_z.push_back(target_p_mom_container[itp].z());
      event.residual_charge.push_back(1);
      if(event.lflag){
	HF1( 1130, target_p_mom_container[itp].Mag());
	if(GFTrackCont.ExtrapolateToPoint(id_p,
					  TVector3(event.GFprodvtx_x_l,
						   event.GFprodvtx_y_l,
						   event.GFprodvtx_z_l),
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_p_repid_container[itp])){

	  closepoint -= TVector3(event.GFprodvtx_x_l,
				 event.GFprodvtx_y_l,
				 event.GFprodvtx_z_l);
	  event.residual_GFdist2prodvtx.push_back(closepoint.Mag());
	}
	else{
	  event.residual_GFdist2prodvtx.push_back(qnan);
	}
	if(GFTrackCont.ExtrapolateToPoint(id_p,
					  TVector3(event.KFprodvtx_x_l,
						   event.KFprodvtx_y_l,
						   event.KFprodvtx_z_l),
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_p_repid_container[itp])){
	  closepoint -= TVector3(event.KFprodvtx_x_l,
				 event.KFprodvtx_y_l,
				 event.KFprodvtx_z_l);
	  event.residual_KFdist2prodvtx.push_back(closepoint.Mag());
	}
	else{
	  event.residual_KFdist2prodvtx.push_back(qnan);
	}
      }
      else{
	HF1( 1100, target_p_mom_container[itp].Mag());
	event.residual_GFdist2prodvtx.push_back(qnan);
	event.residual_KFdist2prodvtx.push_back(qnan);
      }
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
      event.residual_invbeta.push_back(target_pip_invbeta_container[itpip]);
      event.residual_mom.push_back(target_pip_mom_container[itpip].Mag());
      event.residual_mom_x.push_back(target_pip_mom_container[itpip].x());
      event.residual_mom_y.push_back(target_pip_mom_container[itpip].y());
      event.residual_mom_z.push_back(target_pip_mom_container[itpip].z());
      event.residual_charge.push_back(1);
      if(event.lflag){
	HF1( 1131, target_pip_mom_container[itpip].Mag());

	if(GFTrackCont.ExtrapolateToPoint(id_pip,
					  TVector3(event.GFprodvtx_x_l,
						   event.GFprodvtx_y_l,
						   event.GFprodvtx_z_l),
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_pip_repid_container[itpip])){
	  closepoint -= TVector3(event.GFprodvtx_x_l,
				 event.GFprodvtx_y_l,
				 event.GFprodvtx_z_l);
	  event.residual_GFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.residual_GFdist2prodvtx.push_back(qnan);
	if(GFTrackCont.ExtrapolateToPoint(id_pip,
					  TVector3(event.KFprodvtx_x_l,
						   event.KFprodvtx_y_l,
						   event.KFprodvtx_z_l),
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_pip_repid_container[itpip])){
	  closepoint -= TVector3(event.KFprodvtx_x_l,
				 event.KFprodvtx_y_l,
				 event.KFprodvtx_z_l);
	  event.residual_KFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.residual_KFdist2prodvtx.push_back(qnan);
      }
      else{
	HF1( 1101, target_pip_mom_container[itpip].Mag());
	event.residual_GFdist2prodvtx.push_back(qnan);
	event.residual_KFdist2prodvtx.push_back(qnan);
      }
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
      event.residual_invbeta.push_back(target_pim_invbeta_container[itpim]);
      event.residual_mom.push_back(target_pim_mom_container[itpim].Mag());
      event.residual_mom_x.push_back(target_pim_mom_container[itpim].x());
      event.residual_mom_y.push_back(target_pim_mom_container[itpim].y());
      event.residual_mom_z.push_back(target_pim_mom_container[itpim].z());
      event.residual_charge.push_back(-1);
      if(event.lflag){
	HF1( 1132, target_pim_mom_container[itpim].Mag());

	if(GFTrackCont.ExtrapolateToPoint(id_pim,
					  TVector3(event.GFprodvtx_x_l,
						   event.GFprodvtx_y_l,
						   event.GFprodvtx_z_l),
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_pim_repid_container[itpim])){
	  closepoint -= TVector3(event.GFprodvtx_x_l,
				 event.GFprodvtx_y_l,
				 event.GFprodvtx_z_l);
	  event.residual_GFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.residual_GFdist2prodvtx.push_back(qnan);
	if(GFTrackCont.ExtrapolateToPoint(id_pim,
					  TVector3(event.KFprodvtx_x_l,
						   event.KFprodvtx_y_l,
						   event.KFprodvtx_z_l),
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_pim_repid_container[itpim])){
	  closepoint -= TVector3(event.KFprodvtx_x_l,
				 event.KFprodvtx_y_l,
				 event.KFprodvtx_z_l);
	  event.residual_KFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.residual_KFdist2prodvtx.push_back(qnan);
      }
      else{
	HF1( 1102, target_pim_mom_container[itpim].Mag());
	event.residual_GFdist2prodvtx.push_back(qnan);
	event.residual_KFdist2prodvtx.push_back(qnan);
      }
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
      event.residual_invbeta.push_back(target_ppip_invbeta_container[itppip]);
      event.residual_mom.push_back(target_ppip_mom_container[itppip].Mag());
      event.residual_mom_x.push_back(target_ppip_mom_container[itppip].x());
      event.residual_mom_y.push_back(target_ppip_mom_container[itppip].y());
      event.residual_mom_z.push_back(target_ppip_mom_container[itppip].z());
      event.residual_charge.push_back(1);

      if(event.lflag){
	if(GFTrackCont.ExtrapolateToPoint(id_ppip,
					  TVector3(event.GFprodvtx_x_l,
						   event.GFprodvtx_y_l,
						   event.GFprodvtx_z_l),
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_ppip_repid_container[itppip])){
	  closepoint -= TVector3(event.GFprodvtx_x_l,
				 event.GFprodvtx_y_l,
				 event.GFprodvtx_z_l);
	  event.residual_GFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.residual_GFdist2prodvtx.push_back(qnan);
	if(GFTrackCont.ExtrapolateToPoint(id_ppip,
					  TVector3(event.KFprodvtx_x_l,
						   event.KFprodvtx_y_l,
						   event.KFprodvtx_z_l),
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_ppip_repid_container[itppip])){
	  closepoint -= TVector3(event.KFprodvtx_x_l,
				 event.KFprodvtx_y_l,
				 event.KFprodvtx_z_l);
	  event.residual_KFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.residual_KFdist2prodvtx.push_back(qnan);
      }
      else{
	event.residual_GFdist2prodvtx.push_back(qnan);
	event.residual_KFdist2prodvtx.push_back(qnan);
      }
    }
    for(Int_t itep=0; itep<target_ep_id_container.size(); ++itep){
      Int_t id_ep = target_ep_id_container[itep];
      if(id_ep==GFtrackid_decays[0]){
	event.ep_multi -= 1;
	continue;
      }
      event.residual_id.push_back(id_ep);
      event.residual_dist2tgt.push_back(target_ep_dist2tgt_container[itep]);
      event.residual_mass2.push_back(target_ep_mass2_container[itep]);
      event.residual_invbeta.push_back(target_ep_invbeta_container[itep]);
      event.residual_mom.push_back(target_ep_mom_container[itep].Mag());
      event.residual_mom_x.push_back(target_ep_mom_container[itep].x());
      event.residual_mom_y.push_back(target_ep_mom_container[itep].y());
      event.residual_mom_z.push_back(target_ep_mom_container[itep].z());
      event.residual_charge.push_back(1);

      if(event.lflag){
	if(GFTrackCont.ExtrapolateToPoint(id_ep,
					  TVector3(event.GFprodvtx_x_l,
						   event.GFprodvtx_y_l,
						   event.GFprodvtx_z_l),
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_ep_repid_container[itep])){
	  closepoint -= TVector3(event.GFprodvtx_x_l,
				 event.GFprodvtx_y_l,
				 event.GFprodvtx_z_l);
	  event.residual_GFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.residual_GFdist2prodvtx.push_back(qnan);
	if(GFTrackCont.ExtrapolateToPoint(id_ep,
					  TVector3(event.KFprodvtx_x_l,
						   event.KFprodvtx_y_l,
						   event.KFprodvtx_z_l),
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_ep_repid_container[itep])){
	  closepoint -= TVector3(event.KFprodvtx_x_l,
				 event.KFprodvtx_y_l,
				 event.KFprodvtx_z_l);
	  event.residual_KFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.residual_KFdist2prodvtx.push_back(qnan);
      }
      else{
	event.residual_GFdist2prodvtx.push_back(qnan);
	event.residual_KFdist2prodvtx.push_back(qnan);
      }
    }
    for(Int_t item=0; item<target_em_id_container.size(); ++item){
      Int_t id_em = target_em_id_container[item];
      if(id_em==GFtrackid_decays[1]){
	event.em_multi -= 1;
	continue;
      }
      event.residual_id.push_back(id_em);
      event.residual_dist2tgt.push_back(target_em_dist2tgt_container[item]);
      event.residual_mass2.push_back(target_em_mass2_container[item]);
      event.residual_invbeta.push_back(target_em_invbeta_container[item]);
      event.residual_mom.push_back(target_em_mom_container[item].Mag());
      event.residual_mom_x.push_back(target_em_mom_container[item].x());
      event.residual_mom_y.push_back(target_em_mom_container[item].y());
      event.residual_mom_z.push_back(target_em_mom_container[item].z());
      event.residual_charge.push_back(-1);

      if(event.lflag){
	if(GFTrackCont.ExtrapolateToPoint(id_em,
					  TVector3(event.GFprodvtx_x_l,
						   event.GFprodvtx_y_l,
						   event.GFprodvtx_z_l),
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_em_repid_container[item])){
	  closepoint -= TVector3(event.GFprodvtx_x_l,
				 event.GFprodvtx_y_l,
				 event.GFprodvtx_z_l);
	  event.residual_GFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.residual_GFdist2prodvtx.push_back(qnan);
	if(GFTrackCont.ExtrapolateToPoint(id_em,
					  TVector3(event.KFprodvtx_x_l,
						   event.KFprodvtx_y_l,
						   event.KFprodvtx_z_l),
					  closepoint,
					  closepoint_mom,
					  closepoint_tracklen,
					  closepoint_tof,
					  target_em_repid_container[item])){
	  closepoint -= TVector3(event.KFprodvtx_x_l,
				 event.KFprodvtx_y_l,
				 event.KFprodvtx_z_l);
	  event.residual_KFdist2prodvtx.push_back(closepoint.Mag());
	}
	else event.residual_KFdist2prodvtx.push_back(qnan);
      }
      else{
	event.residual_GFdist2prodvtx.push_back(qnan);
	event.residual_KFdist2prodvtx.push_back(qnan);
      }
    }

    event.residual_multi = event.ep_multi + event.em_multi + event.pip_multi + event.p_multi + event.pim_multi + event.ppip_multi;

    if(event.lflag){
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
    }
    else if(event.residual_multi==1 && event.pim_multi==1) event.pimflag = true;
    else if(event.residual_multi==0) event.emptyflag = true;

#if DebugDisp
    std::cout<<"Saving gamma calndidates"<<std::endl;
#endif
    event.g_multi = gamma_mom_container.size();
    for(Int_t itg=0; itg<gamma_mom_container.size(); ++itg){
      Int_t epid = gamma_ep_id_container[itg];
      Int_t emid = gamma_em_id_container[itg];
      if(epid==GFtrackid_decays[0] || epid==GFtrackid_decays[1]){
	event.g_multi -=1;
	continue;
      }
      auto iter_ep = find(event.epidgamma.begin(), event.epidgamma.end(), epid);
      if(iter_ep != event.epidgamma.end()){
	event.g_multi -=1;
	continue;
      }
      auto iter_em = find(event.emidgamma.begin(), event.emidgamma.end(), emid);
      if(iter_em != event.emidgamma.end()){
	event.g_multi -=1;
	continue;
      }

      Double_t dist = gamma_epemdist_container[itg];
      TVector3 ep_mom = gamma_ep_mom_container[itg];
      TVector3 em_mom = gamma_em_mom_container[itg];
      TVector3 gamma_mom = gamma_mom_container[itg];
      TVector3 gamma_vtx = gamma_decayvtx_container[itg];

      event.epidgamma.push_back(epid);
      event.emidgamma.push_back(emid);
      event.epmomgamma.push_back(ep_mom.Mag());
      event.epmomgamma_x.push_back(ep_mom.x());
      event.epmomgamma_y.push_back(ep_mom.y());
      event.epmomgamma_z.push_back(ep_mom.z());
      event.emmomgamma.push_back(em_mom.Mag());
      event.emmomgamma_x.push_back(em_mom.x());
      event.emmomgamma_y.push_back(em_mom.y());
      event.emmomgamma_z.push_back(em_mom.z());
      event.momgamma.push_back(gamma_mom.Mag());
      event.momgamma_x.push_back(gamma_mom.x());
      event.momgamma_y.push_back(gamma_mom.y());
      event.momgamma_z.push_back(gamma_mom.z());
      event.epidistgamma.push_back(dist);
      event.vtxgamma_x.push_back(gamma_vtx.x());
      event.vtxgamma_y.push_back(gamma_vtx.y());
      event.vtxgamma_z.push_back(gamma_vtx.z());
    }
  } //if(!event.xiflag && !event.llflag)

  event.accident_multi = target_accidental_id_container.size();
  for(Int_t id=0; id<target_accidental_id_container.size(); ++id){
    event.accident_id.push_back(id);
  }

#if DebugDisp
  std::cout<<"Debugging 1. pi- pi+ pair candidate searching"<<std::endl;
#endif

  //pi+&pi- pair
  std::vector<Int_t> pipair_pip_id_container, pipair_pim_id_container;
  std::vector<Double_t> pipair_reconL_mass_container, pipair_recon_mass_container;
  std::vector<Double_t> pipair_pipidist_container;
  std::vector<TVector3> pipair_mom_container, pipair_pip_mom_container, pipair_pim_mom_container;
  {
    for(Int_t it1=0;it1<ntTpc;it1++){
      if(event.isElectron[it1]==1) continue;
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
	if(event.isElectron[it2]==1) continue;
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

	TLorentzVector Lpim(pim_mom, TMath::Hypot(pim_mom.Mag(), PionMass));
	TLorentzVector Lpip(pip_mom, TMath::Hypot(pip_mom.Mag(), PionMass));
	TLorentzVector Lpipi = Lpip + Lpim;

	TLorentzVector Lpip_plike(pip_mom, TMath::Hypot(pip_mom.Mag(), ProtonMass));
	TLorentzVector Llambda = Lpip_plike + Lpim;

	if(TMath::Abs(lambda_vert.x()) > 250. ||
	   TMath::Abs(lambda_vert.z()) > 250. ||
	   TMath::Abs(lambda_vert.y()) > 250.) continue; //Vertex cut

	Double_t pip_vertex_dist; Double_t pim_vertex_dist;
	if(!Kinematics::HelixDirection(lambda_vert, pip_start, pip_end, pip_vertex_dist) ||
	   !Kinematics::HelixDirection(lambda_vert, pim_start, pim_end, pim_vertex_dist)) continue;

	if(pip_vertex_dist > pi_vtx_distcut) continue;
	if(pim_vertex_dist > pi_vtx_distcut) continue;

	Double_t ltarget_dist;
	TVector3 ltarget_vtx = Kinematics::CalcCloseDistLambda(tgtpos,
							       lambda_vert,
							       lambda_mom,
							       ltarget_dist);
	if(ltarget_dist > ltarget_distcut) continue;
	if(pipi_dist < pipi_distcut){
	  pipair_pip_id_container.push_back(it1);
	  pipair_pim_id_container.push_back(it2);
	  pipair_reconL_mass_container.push_back(Llambda.M());
	  pipair_recon_mass_container.push_back(Lpipi.M());
	  pipair_pipidist_container.push_back(pipi_dist);
	  pipair_mom_container.push_back(lambda_mom);
	  pipair_pip_mom_container.push_back(pip_mom);
	  pipair_pim_mom_container.push_back(pim_mom);
	}
      } //it2
    } //it1
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
    event.reconmassPipair.push_back(pipair_recon_mass_container[icombi]);
    event.pipidistPipair.push_back(pipair_pipidist_container[icombi]);
  }

#if DebugDisp
  std::cout<<"Debugging 2. p-bar-like p and pi- pair candidate searching"<<std::endl;
#endif
  std::vector<Int_t> reconfailed_L_p_id_container, reconfailed_L_pi_id_container;
  std::vector<Double_t> reconfailed_L_mass_container, reconfailed_L_ppidist_container;
  std::vector<TVector3> reconfailed_L_mom_container, reconfailed_L_p_mom_container,
    reconfailed_L_pi_mom_container, reconfailed_L_vtx_container;

  {
    for(Int_t it1=0;it1<ntTpc;it1++){
      if(event.isElectron[it1]==1) continue;
      if(event.isK18[it1]==1) continue;
      if(event.isKurama[it1]==1) continue;
      if(event.isBeam[it1]==1) continue;
      if(event.isXi[it1]==1) continue;
      if(event.isAccidental[it1]==1) continue;
      if(event.chargeIndistinguishable[it1]!=1) continue;
      if((event.pid_inverted[it1]&4)!=4) continue;
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
	if(event.isElectron[it2]==1) continue;
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
  std::cout<<"Debugging 3. p and p-like pi- pair candidate searching"<<std::endl;
#endif
  {
    for(Int_t it1=0;it1<ntTpc;it1++){
      if(event.isElectron[it1]==1) continue;
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
	if(event.isElectron[it2]==1) continue;
	if(event.isK18[it2]==1) continue;
	if(event.isKurama[it2]==1) continue;
	if(event.isBeam[it2]==1) continue;
	if(event.isXi[it2]==1) continue;
	if(event.isAccidental[it2]==1) continue;
	if(event.chargeIndistinguishable[it1]!=1) continue;
	if((event.pid_inverted[it1]&1)!=1) continue;
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

  tree->Branch( "nclTpc", &event.nclTpc );
  tree->Branch( "remain_nclTpc", &event.nclTpc );
#if SaveRawData
  tree->Branch( "remain_cluster_x", &event.remain_cluster_x );
  tree->Branch( "remain_cluster_y", &event.remain_cluster_y );
  tree->Branch( "remain_cluster_z", &event.remain_cluster_z );
  tree->Branch( "remain_cluster_de", &event.remain_cluster_de );
  tree->Branch( "remain_cluster_size", &event.remain_cluster_size );
  tree->Branch( "remain_cluster_layer", &event.remain_cluster_layer );
  tree->Branch( "remain_cluster_row_center", &event.remain_cluster_row_center );
  tree->Branch( "remain_cluster_mrow", &event.remain_cluster_mrow );
  tree->Branch( "remain_cluster_de_center", &event.remain_cluster_de_center );
  tree->Branch( "remain_cluster_x_center", &event.remain_cluster_x_center );
  tree->Branch( "remain_cluster_y_center", &event.remain_cluster_y_center );
  tree->Branch( "remain_cluster_z_center", &event.remain_cluster_z_center );
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
  tree->Branch( "isElectron", &event.isElectron );
  tree->Branch( "nsigma_triton", &event.nsigma_triton );
  tree->Branch( "nsigma_deutron", &event.nsigma_deutron );
  tree->Branch( "nsigma_proton", &event.nsigma_proton );
  tree->Branch( "nsigma_kaon", &event.nsigma_kaon );
  tree->Branch( "nsigma_pion", &event.nsigma_pion );
  tree->Branch( "nsigma_electron", &event.nsigma_electron );
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
  tree->Branch( "track_cluster_size", &event.track_cluster_size);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch( "track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch( "track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch( "track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch( "track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch( "track_cluster_row_center", &event.track_cluster_row_center);
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

  tree->Branch("GFcharge", &event.GFcharge);
  tree->Branch("GFchisqr", &event.GFchisqr);
  tree->Branch("GFtof", &event.GFtof);
  //tree->Branch("GFtracklen", &event.GFtracklen);
  tree->Branch("GFpval", &event.GFpval);
  tree->Branch("GFfitstatus", &event.GFfitstatus);
  tree->Branch("GFpdgcode", &event.GFpdgcode);
  tree->Branch("GFnhtrack", &event.GFnhtrack);
  tree->Branch("GFlayer", &event.GFlayer);
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

  tree->Branch("GFntTpc_target", &event.GFntTpc_inside);
  tree->Branch("GFprodvtx_x", &event.GFprodvtx_x);
  tree->Branch("GFprodvtx_y", &event.GFprodvtx_y);
  tree->Branch("GFprodvtx_z", &event.GFprodvtx_z);

  //extrapolation
  tree->Branch("GFinside", &event.GFinside);
  tree->Branch("GFtracklen", &event.GFtracklen);
  tree->Branch("GFtrack2vtxdist", &event.GFtrack2vtxdist);
  tree->Branch("GFcalctof", &event.GFcalctof);
  tree->Branch("GFsegHtof", &event.GFsegHtof);
  tree->Branch("GFtofHtof", &event.GFtofHtof);
  tree->Branch("GFtdiffHtof", &event.GFtdiffHtof);
  tree->Branch("GFposHtof", &event.GFposHtof);
  tree->Branch("GFposx", &event.GFposx);
  tree->Branch("GFposy", &event.GFposy);
  tree->Branch("GFposz", &event.GFposz);
  tree->Branch("GFinvbeta", &event.GFinvbeta);
  tree->Branch("GFm2", &event.GFm2);
  tree->Branch("nsigma_tritonHtof", &event.nsigma_tritonHtof);
  tree->Branch("nsigma_deutronHtof", &event.nsigma_deutronHtof);
  tree->Branch("nsigma_protonHtof", &event.nsigma_protonHtof);
  tree->Branch("nsigma_kaonHtof", &event.nsigma_kaonHtof);
  tree->Branch("nsigma_pionHtof", &event.nsigma_pionHtof);
  tree->Branch("nsigma_electronHtof", &event.nsigma_electronHtof);

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

  tree->Branch( "nvtxTpcClustered", &event.nvtxTpcClustered );
  tree->Branch( "clusteredVtx_x", &event.clusteredVtx_x );
  tree->Branch( "clusteredVtx_y", &event.clusteredVtx_y );
  tree->Branch( "clusteredVtx_z", &event.clusteredVtx_z );
  tree->Branch( "clusteredVtxid", &event.clusteredVtxid );

  tree->Branch( "ncombiReconFailedLambda", &event.ncombiLreconfailed  );
  tree->Branch( "ReconFailedLambdaPId", &event.pidLreconfailed);
  tree->Branch( "ReconFailedLambdaPiId", &event.piidLreconfailed);
  tree->Branch( "ReconFailedLambdaMass", &event.LmassLreconfailed);
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
  tree->Branch( "PiPairReconMass", & event.reconmassPipair);
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
  tree->Branch( "m2Kurama", &event.m2Kurama);
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
  tree->Branch( "xbTPC", &event.xbTPC);
  tree->Branch( "ybTPC", &event.ybTPC);
  tree->Branch( "ubTPC", &event.ubTPC);
  tree->Branch( "vbTPC", &event.vbTPC);
  tree->Branch( "xsTPC", &event.xsTPC);
  tree->Branch( "ysTPC", &event.ysTPC);
  tree->Branch( "usTPC", &event.usTPC);
  tree->Branch( "vsTPC", &event.vsTPC);

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

  tree->Branch("Xiflag", &event.xiflag);
  tree->Branch("XiPflag", &event.xipflag);
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

  //Alternative p, pi pairing for LL events
  tree->Branch("GFLambdaMass1_Alt", &event.GFlmass_alter1);
  tree->Branch("GFLambdaDecayVtx_x1_Alt", &event.GFldecayvtx_x_alter1);
  tree->Branch("GFLambdaDecayVtx_y1_Alt", &event.GFldecayvtx_y_alter1);
  tree->Branch("GFLambdaDecayVtx_z1_Alt", &event.GFldecayvtx_z_alter1);
  tree->Branch("GFLambdaMom1_Alt", &event.GFlmom_alter1);
  tree->Branch("GFLambdaMom_x1_Alt", &event.GFlmom_x_alter1);
  tree->Branch("GFLambdaMom_y1_Alt", &event.GFlmom_y_alter1);
  tree->Branch("GFLambdaMom_z1_Alt", &event.GFlmom_z_alter1);
  tree->Branch("GFLambdaVtxCloseDist1_Alt", &event.GFppi_dist_alter1);
  tree->Branch("GFLambdaTargetCloseDist1_Alt", &event.GFltarget_dist_alter1);
  tree->Branch("GFLambdaTarget_x1_Alt", &event.GFltargetvtx_x_alter1);
  tree->Branch("GFLambdaTarget_y1_Alt", &event.GFltargetvtx_y_alter1);
  tree->Branch("GFLambdaTarget_z1_Alt", &event.GFltargetvtx_z_alter1);
  tree->Branch("GFLambdaTargetCenterCloseDist1_Alt", &event.GFltargetcenter_dist_alter1);
  tree->Branch("GFLambdaTargetCenter_x1_Alt", &event.GFltargetcenter_x_alter1);
  tree->Branch("GFLambdaTargetCenter_y1_Alt", &event.GFltargetcenter_y_alter1);
  tree->Branch("GFLambdaTargetCenter_z1_Alt", &event.GFltargetcenter_z_alter1);

  tree->Branch("GFLambdaMass2_Alt", &event.GFlmass_alter2);
  tree->Branch("GFLambdaDecayVtx_x2_Alt", &event.GFldecayvtx_x_alter2);
  tree->Branch("GFLambdaDecayVtx_y2_Alt", &event.GFldecayvtx_y_alter2);
  tree->Branch("GFLambdaDecayVtx_z2_Alt", &event.GFldecayvtx_z_alter2);
  tree->Branch("GFLambdaMom2_Alt", &event.GFlmom_alter2);
  tree->Branch("GFLambdaMom_x2_Alt", &event.GFlmom_x_alter2);
  tree->Branch("GFLambdaMom_y2_Alt", &event.GFlmom_y_alter2);
  tree->Branch("GFLambdaMom_z2_Alt", &event.GFlmom_z_alter2);
  tree->Branch("GFLambdaVtxCloseDist2_Alt", &event.GFppi_dist_alter2);
  tree->Branch("GFLambdaTargetCloseDist2_Alt", &event.GFltarget_dist_alter2);
  tree->Branch("GFLambdaTarget_x2_Alt", &event.GFltargetvtx_x_alter2);
  tree->Branch("GFLambdaTarget_y2_Alt", &event.GFltargetvtx_y_alter2);
  tree->Branch("GFLambdaTarget_z2_Alt", &event.GFltargetvtx_z_alter2);
  tree->Branch("GFLambdaTargetCenterCloseDist2_Alt", &event.GFltargetcenter_dist_alter2);
  tree->Branch("GFLambdaTargetCenter_x2_Alt", &event.GFltargetcenter_x_alter2);
  tree->Branch("GFLambdaTargetCenter_y2_Alt", &event.GFltargetcenter_y_alter2);
  tree->Branch("GFLambdaTargetCenter_z2_Alt", &event.GFltargetcenter_z_alter2);

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
  tree->Branch("PhiCosKK", &event.phicosKK);
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
  tree->Branch("GFPhiCosKK", &event.GFphicosKK);
  tree->Branch("GFPhiMom", &event.GFphimom);
  tree->Branch("GFPhiMom_x", &event.GFphimom_x);
  tree->Branch("GFPhiMom_y", &event.GFphimom_y);
  tree->Branch("GFPhiMom_z", &event.GFphimom_z);
  tree->Branch("GFPhiVtxCloseDist", &event.GFkk_dist);
  tree->Branch("GFPhiProductionVtxCloseDist", &event.GFphiprodvtx_dist);
  tree->Branch("GFPhiKmMassSquare", &event.GFphi_km_mass2);
  tree->Branch("GFPhiKmInvbeta", &event.GFphi_km_invbeta);
  tree->Branch("GFPhiKpMassSquare", &event.GFphi_kp_mass2);
  tree->Branch("GFPhiKpInvbeta", &event.GFphi_kp_invbeta);
  tree->Branch("GFPhiDecaysMom", &event.GFphidecays_mom);
  tree->Branch("GFPhiDecaysMom_x", &event.GFphidecays_mom_x);
  tree->Branch("GFPhiDecaysMom_y", &event.GFphidecays_mom_y);
  tree->Branch("GFPhiDecaysMom_z", &event.GFphidecays_mom_z);

  tree->Branch("GFPhiMassKurama", &event.GFphimass_wKurama);
  tree->Branch("GFPhiMomKurama", &event.GFphimom_wKurama);
  tree->Branch("GFPhiMomKurama_x", &event.GFphimom_x_wKurama);
  tree->Branch("GFPhiMomKurama_y", &event.GFphimom_y_wKurama);
  tree->Branch("GFPhiMomKurama_z", &event.GFphimom_z_wKurama);

  tree->Branch("LPiflag", &event.lpiflag);
  tree->Branch("PiPiflag", &event.pipiflag);
  tree->Branch("Pimflag", &event.pimflag);
  tree->Branch("Emptyflag", &event.emptyflag);

  //for decay particles
  tree->Branch("GFXiDecaysNhit", &event.GFxidecays_nhtrack);
  tree->Branch("GFXiDecaysChisqr", &event.GFxidecays_chisqr);
  tree->Branch("GFXiDecaysCharge", &event.GFxidecays_charge);
  tree->Branch("GFXiDecaysTof", &event.GFxidecays_tracktof);
  tree->Branch("GFXiDecaysPval", &event.GFxidecays_pval);
  tree->Branch("GFXiDecaysPdgcode", &event.GFxidecays_pdgcode);
  tree->Branch("GFXiDecaysHtofId", &event.GFxidecays_htofid);
  tree->Branch("GFXiDecaysTrackLen", &event.GFxidecays_tracklen);
  tree->Branch("GFXiDecaysTrackTof", &event.GFxidecays_tof);
  tree->Branch("GFXiDecaysMassSquare", &event.GFxidecays_mass2);
  tree->Branch("GFXiDecaysInvbeta", &event.GFxidecays_invbeta);
  tree->Branch("GFXiDecaysMom", &event.GFxidecays_mom);
  tree->Branch("GFXiDecaysMom_x", &event.GFxidecays_mom_x);
  tree->Branch("GFXiDecaysMom_y", &event.GFxidecays_mom_y);
  tree->Branch("GFXiDecaysMom_z", &event.GFxidecays_mom_z);
  tree->Branch("GFXiDecaysMomCM", &event.GFxidecays_CMmom);
  tree->Branch("GFXiDecaysMomCM_x", &event.GFxidecays_CMmom_x);
  tree->Branch("GFXiDecaysMomCM_y", &event.GFxidecays_CMmom_y);
  tree->Branch("GFXiDecaysMomCM_z", &event.GFxidecays_CMmom_z);
  tree->Branch("GFXiDecaysMomLoss", &event.GFxidecays_momloss);
  tree->Branch("GFXiDecaysELoss", &event.GFxidecays_eloss);

  tree->Branch("XiDecaysTrackId", &event.xidecays_id);
  tree->Branch("XiDecaysMom", &event.xidecays_mom);
  tree->Branch("XiDecaysMom_x", &event.xidecays_mom_x);
  tree->Branch("XiDecaysMom_y", &event.xidecays_mom_y);
  tree->Branch("XiDecaysMom_z", &event.xidecays_mom_z);
  tree->Branch("XiDecaysMomCM", &event.xidecays_CMmom);
  tree->Branch("XiDecaysMomCM_x", &event.xidecays_CMmom_x);
  tree->Branch("XiDecaysMomCM_y", &event.xidecays_CMmom_y);
  tree->Branch("XiDecaysMomCM_z", &event.xidecays_CMmom_z);

  tree->Branch("GFLLDecaysNhit", &event.GFlldecays_nhtrack);
  tree->Branch("GFLLDecaysChisqr", &event.GFlldecays_chisqr);
  tree->Branch("GFLLDecaysCharge", &event.GFlldecays_charge);
  tree->Branch("GFLLDecaysTof", &event.GFlldecays_tof);
  tree->Branch("GFLLDecaysPval", &event.GFlldecays_pval);
  tree->Branch("GFLLDecaysPdgcode", &event.GFlldecays_pdgcode);
  tree->Branch("GFLLDecaysHtofId", &event.GFlldecays_htofid);
  tree->Branch("GFLLDecaysTrackLen", &event.GFlldecays_tracklen);
  tree->Branch("GFLLDecaysTrackTof", &event.GFlldecays_tof);
  tree->Branch("GFLLDecaysMassSquare", &event.GFlldecays_mass2);
  tree->Branch("GFLLDecaysInvbeta", &event.GFlldecays_invbeta);
  tree->Branch("GFLLDecaysMom", &event.GFlldecays_mom);
  tree->Branch("GFLLDecaysMom_x", &event.GFlldecays_mom_x);
  tree->Branch("GFLLDecaysMom_y", &event.GFlldecays_mom_y);
  tree->Branch("GFLLDecaysMom_z", &event.GFlldecays_mom_z);
  tree->Branch("GFLLDecaysMomCM", &event.GFlldecays_CMmom);
  tree->Branch("GFLLDecaysMomCM_x", &event.GFlldecays_CMmom_x);
  tree->Branch("GFLLDecaysMomCM_y", &event.GFlldecays_CMmom_y);
  tree->Branch("GFLLDecaysMomCM_z", &event.GFlldecays_CMmom_z);
  tree->Branch("GFLLDecaysMomLoss", &event.GFlldecays_momloss);
  tree->Branch("GFLLDecaysELoss", &event.GFlldecays_eloss);

  tree->Branch("LLDecaysTrackId", &event.lldecays_id);
  tree->Branch("LLDecaysMom", &event.lldecays_mom);
  tree->Branch("LLDecaysMom_x", &event.lldecays_mom_x);
  tree->Branch("LLDecaysMom_y", &event.lldecays_mom_y);
  tree->Branch("LLDecaysMom_z", &event.lldecays_mom_z);
  tree->Branch("LLDecaysMomCM", &event.lldecays_CMmom);
  tree->Branch("LLDecaysMomCM_x", &event.lldecays_CMmom_x);
  tree->Branch("LLDecaysMomCM_y", &event.lldecays_CMmom_y);
  tree->Branch("LLDecaysMomCM_z", &event.lldecays_CMmom_z);

  tree->Branch("GFDecaysNhit", &event.GFdecays_nhtrack);
  tree->Branch("GFDecaysChisqr", &event.GFdecays_chisqr);
  tree->Branch("GFDecaysCharge", &event.GFdecays_charge);
  tree->Branch("GFDecaysTof", &event.GFdecays_tof);
  tree->Branch("GFDecaysPval", &event.GFdecays_pval);
  tree->Branch("GFDecaysPdgcode", &event.GFdecays_pdgcode);
  tree->Branch("GFDecaysHtofId", &event.GFdecays_htofid);
  tree->Branch("GFDecaysTrackLen", &event.GFdecays_tracklen);
  tree->Branch("GFDecaysTrackTof", &event.GFdecays_tof);
  tree->Branch("GFDecaysMassSquare", &event.GFdecays_mass2);
  tree->Branch("GFDecaysInvbeta", &event.GFdecays_invbeta);
  tree->Branch("GFDecaysMom", &event.GFdecays_mom);
  tree->Branch("GFDecaysMom_x", &event.GFdecays_mom_x);
  tree->Branch("GFDecaysMom_y", &event.GFdecays_mom_y);
  tree->Branch("GFDecaysMom_z", &event.GFdecays_mom_z);
  tree->Branch("GFDecaysMomCM", &event.GFdecays_CMmom);
  tree->Branch("GFDecaysMomCM_x", &event.GFdecays_CMmom_x);
  tree->Branch("GFDecaysMomCM_y", &event.GFdecays_CMmom_y);
  tree->Branch("GFDecaysMomCM_z", &event.GFdecays_CMmom_z);
  tree->Branch("GFDecaysMomLoss", &event.GFdecays_momloss);
  tree->Branch("GFDecaysELoss", &event.GFdecays_eloss);

  tree->Branch("DecaysTrackId", &event.decays_id);
  tree->Branch("DecaysMom", &event.decays_mom);
  tree->Branch("DecaysMom_x", &event.decays_mom_x);
  tree->Branch("DecaysMom_y", &event.decays_mom_y);
  tree->Branch("DecaysMom_z", &event.decays_mom_z);
  tree->Branch("DecaysMomCM", &event.decays_CMmom);
  tree->Branch("DecaysMomCM_x", &event.decays_CMmom_x);
  tree->Branch("DecaysMomCM_y", &event.decays_CMmom_y);
  tree->Branch("DecaysMomCM_z", &event.decays_CMmom_z);

  //Remaining p, pi after Xi, L searching
  //Multiplicity means tracks comes from the target
  tree->Branch("AccidentalMultiplicity", &event.accident_multi);
  tree->Branch("AccidentalTrackId", &event.accident_id);

  tree->Branch("LLResidualsMultiplicity", &event.llresidual_multi);
  tree->Branch("LLpMultiplicity", &event.llp_multi);
  tree->Branch("LLpipMultiplicity", &event.llpip_multi);
  tree->Branch("LLpimMultiplicity", &event.llpim_multi);
  tree->Branch("LLepMultiplicity", &event.llep_multi);
  tree->Branch("LLemMultiplicity", &event.llem_multi);
  tree->Branch("LLppipMultiplicity", &event.llppip_multi);
  tree->Branch("LLResidualsTrackId", &event.llresidual_id);
  tree->Branch("LLResidualsMassSquare", &event.llresidual_mass2);
  tree->Branch("LLResidualsInvbeta", &event.llresidual_invbeta);
  tree->Branch("LLResidualsCloseDistTgt", &event.llresidual_dist2tgt);
  tree->Branch("LLResidualsGFProductionVtxCloseDist", &event.llresidual_GFdist2prodvtx);
  tree->Branch("LLResidualsKFProductionVtxCloseDist", &event.llresidual_KFdist2prodvtx);
  tree->Branch("LLResidualsMom", &event.llresidual_mom);
  tree->Branch("LLResidualsMom_x", &event.llresidual_mom_x);
  tree->Branch("LLResidualsMom_y", &event.llresidual_mom_y);
  tree->Branch("LLResidualsMom_z", &event.llresidual_mom_z);
  tree->Branch("LLResidualsCharge", &event.llresidual_charge);

  tree->Branch("XiResidualsMultiplicity", &event.xiresidual_multi);
  tree->Branch("XipMultiplicity", &event.xip_multi);
  tree->Branch("XipipMultiplicity", &event.xipip_multi);
  tree->Branch("XipimMultiplicity", &event.xipim_multi);
  tree->Branch("XiepMultiplicity", &event.xiep_multi);
  tree->Branch("XiemMultiplicity", &event.xiem_multi);
  tree->Branch("XippipMultiplicity", &event.xippip_multi);
  tree->Branch("XiResidualsTrackId", &event.xiresidual_id);
  tree->Branch("XiResidualsMassSquare", &event.xiresidual_mass2);
  tree->Branch("XiResidualsInvbeta", &event.xiresidual_invbeta);
  tree->Branch("XiResidualsCloseDistTgt", &event.xiresidual_dist2tgt);
  tree->Branch("XiResidualsGFProductionVtxCloseDist", &event.xiresidual_GFdist2prodvtx);
  tree->Branch("XiResidualsKFProductionVtxCloseDist", &event.xiresidual_KFdist2prodvtx);
  tree->Branch("XiResidualsMom", &event.xiresidual_mom);
  tree->Branch("XiResidualsMom_x", &event.xiresidual_mom_x);
  tree->Branch("XiResidualsMom_y", &event.xiresidual_mom_y);
  tree->Branch("XiResidualsMom_z", &event.xiresidual_mom_z);
  tree->Branch("XiResidualsCharge", &event.xiresidual_charge);

  tree->Branch("ResidualsMultiplicity", &event.residual_multi);
  tree->Branch("pMultiplicity", &event.p_multi);
  tree->Branch("pipMultiplicity", &event.pip_multi);
  tree->Branch("pimMultiplicity", &event.pim_multi);
  tree->Branch("epMultiplicity", &event.ep_multi);
  tree->Branch("emMultiplicity", &event.em_multi);
  tree->Branch("ppipMultiplicity", &event.ppip_multi);
  tree->Branch("ResidualsTrackId", &event.residual_id);
  tree->Branch("ResidualsMassSquare", &event.residual_mass2);
  tree->Branch("ResidualsInvbeta", &event.residual_invbeta);
  tree->Branch("ResidualsCloseDistTgt", &event.residual_dist2tgt);
  tree->Branch("ResidualsGFProductionVtxCloseDist", &event.residual_GFdist2prodvtx);
  tree->Branch("ResidualsKFProductionVtxCloseDist", &event.residual_KFdist2prodvtx);
  tree->Branch("ResidualsMom", &event.residual_mom);
  tree->Branch("ResidualsMom_x", &event.residual_mom_x);
  tree->Branch("ResidualsMom_y", &event.residual_mom_y);
  tree->Branch("ResidualsMom_z", &event.residual_mom_z);
  tree->Branch("ResidualsCharge", &event.residual_charge);

  //Kinematic Fitting
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

  //for decay particles
  tree->Branch("KFLLDecaysMom", &event.KFlldecays_mom);
  tree->Branch("KFLLDecaysMom_x", &event.KFlldecays_mom_x);
  tree->Branch("KFLLDecaysMom_y", &event.KFlldecays_mom_y);
  tree->Branch("KFLLDecaysMom_z", &event.KFlldecays_mom_z);
  tree->Branch("KFLLDecaysMomCM", &event.KFlldecays_CMmom);
  tree->Branch("KFLLDecaysMomCM_x", &event.KFlldecays_CMmom_x);
  tree->Branch("KFLLDecaysMomCM_y", &event.KFlldecays_CMmom_y);
  tree->Branch("KFLLDecaysMomCM_z", &event.KFlldecays_CMmom_z);

  tree->Branch("KFXiDecaysMom", &event.KFxidecays_mom);
  tree->Branch("KFXiDecaysMom_x", &event.KFxidecays_mom_x);
  tree->Branch("KFXiDecaysMom_y", &event.KFxidecays_mom_y);
  tree->Branch("KFXiDecaysMom_z", &event.KFxidecays_mom_z);
  tree->Branch("KFXiDecaysMomCM", &event.KFxidecays_CMmom);
  tree->Branch("KFXiDecaysMomCM_x", &event.KFxidecays_CMmom_x);
  tree->Branch("KFXiDecaysMomCM_y", &event.KFxidecays_CMmom_y);
  tree->Branch("KFXiDecaysMomCM_z", &event.KFxidecays_CMmom_z);

  tree->Branch("KFDecaysMom", &event.KFdecays_mom);
  tree->Branch("KFDecaysMom_x", &event.KFdecays_mom_x);
  tree->Branch("KFDecaysMom_y", &event.KFdecays_mom_y);
  tree->Branch("KFDecaysMom_z", &event.KFdecays_mom_z);
  tree->Branch("KFDecaysMomCM", &event.KFdecays_CMmom);
  tree->Branch("KFDecaysMomCM_x", &event.KFdecays_CMmom_x);
  tree->Branch("KFDecaysMomCM_y", &event.KFdecays_CMmom_y);
  tree->Branch("KFDecaysMomCM_z", &event.KFdecays_CMmom_z);

  tree->Branch("LLnGamma", &event.llg_multi);
  tree->Branch("LLGammaEpTrackId", &event.llepidgamma);
  tree->Branch("LLGammaEmTrackId", &event.llemidgamma);
  tree->Branch("LLGammaMomId", &event.llepmomgamma);
  tree->Branch("LLGammaEpMom_x", &event.llepmomgamma_x);
  tree->Branch("LLGammaEpMom_y", &event.llepmomgamma_y);
  tree->Branch("LLGammaEpMom_z", &event.llepmomgamma_z);
  tree->Branch("LLGammaEmMom", &event.llemmomgamma);
  tree->Branch("LLGammaEmMom_x", &event.llemmomgamma_x);
  tree->Branch("LLGammaEmMom_y", &event.llemmomgamma_y);
  tree->Branch("LLGammaEmMom_z", &event.llemmomgamma_z);
  tree->Branch("LLGammaMom", &event.llmomgamma);
  tree->Branch("LLGammaMom_x", &event.llmomgamma_x);
  tree->Branch("LLGammaMom_y", &event.llmomgamma_y);
  tree->Branch("LLGammaMom_z", &event.llmomgamma_z);
  tree->Branch("LLGammaVtxCloseDist", &event.llepidistgamma);
  tree->Branch("LLGammaDecayVtx_x", &event.llvtxgamma_x);
  tree->Branch("LLGammaDecayVtx_y", &event.llvtxgamma_y);
  tree->Branch("LLGammaDecayVtx_z", &event.llvtxgamma_z);

  tree->Branch("XinGamma", &event.xig_multi);
  tree->Branch("XiGammaEpTrackId", &event.xiepidgamma);
  tree->Branch("XiGammaEmTrackId", &event.xiemidgamma);
  tree->Branch("XiGammaMomId", &event.xiepmomgamma);
  tree->Branch("XiGammaEpMom_x", &event.xiepmomgamma_x);
  tree->Branch("XiGammaEpMom_y", &event.xiepmomgamma_y);
  tree->Branch("XiGammaEpMom_z", &event.xiepmomgamma_z);
  tree->Branch("XiGammaEmMom", &event.xiemmomgamma);
  tree->Branch("XiGammaEmMom_x", &event.xiemmomgamma_x);
  tree->Branch("XiGammaEmMom_y", &event.xiemmomgamma_y);
  tree->Branch("XiGammaEmMom_z", &event.xiemmomgamma_z);
  tree->Branch("XiGammaMom", &event.ximomgamma);
  tree->Branch("XiGammaMom_x", &event.ximomgamma_x);
  tree->Branch("XiGammaMom_y", &event.ximomgamma_y);
  tree->Branch("XiGammaMom_z", &event.ximomgamma_z);
  tree->Branch("XiGammaVtxCloseDist", &event.xiepidistgamma);
  tree->Branch("XiGammaDecayVtx_x", &event.xivtxgamma_x);
  tree->Branch("XiGammaDecayVtx_y", &event.xivtxgamma_y);
  tree->Branch("XiGammaDecayVtx_z", &event.xivtxgamma_z);

  tree->Branch("nGamma", &event.g_multi);
  tree->Branch("GammaEpTrackId", &event.epidgamma);
  tree->Branch("GammaEmTrackId", &event.emidgamma);
  tree->Branch("GammaMomId", &event.epmomgamma);
  tree->Branch("GammaEpMom_x", &event.epmomgamma_x);
  tree->Branch("GammaEpMom_y", &event.epmomgamma_y);
  tree->Branch("GammaEpMom_z", &event.epmomgamma_z);
  tree->Branch("GammaEmMom", &event.emmomgamma);
  tree->Branch("GammaEmMom_x", &event.emmomgamma_x);
  tree->Branch("GammaEmMom_y", &event.emmomgamma_y);
  tree->Branch("GammaEmMom_z", &event.emmomgamma_z);
  tree->Branch("GammaMom", &event.momgamma);
  tree->Branch("GammaMom_x", &event.momgamma_x);
  tree->Branch("GammaMom_y", &event.momgamma_y);
  tree->Branch("GammaMom_z", &event.momgamma_z);
  tree->Branch("GammaVtxCloseDist", &event.epidistgamma);
  tree->Branch("GammaDecayVtx_x", &event.vtxgamma_x);
  tree->Branch("GammaDecayVtx_y", &event.vtxgamma_y);
  tree->Branch("GammaDecayVtx_z", &event.vtxgamma_z);

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
  src.m2Org = new TTreeReaderValue<std::vector<Double_t>>( *reader, "m2Org" );
  src.thetaKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaKurama" );
  src.xtgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtKurama" );
  src.ytgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtKurama" );
  src.utgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtKurama" );
  src.vtgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtKurama" );
  src.pathwcKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pathwcKurama" );
  src.xin = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xin" );
  src.yin = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yin" );
  src.zin = new TTreeReaderValue<std::vector<Double_t>>( *reader, "zin" );
  src.pxin = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pxin" );
  src.pyin = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pyin" );
  src.pzin = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pzin" );
  src.xout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xout" );
  src.yout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yout" );
  src.zout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "zout" );
  src.pxout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pxout" );
  src.pyout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pyout" );
  src.pzout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pzout" );

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

  src.pOrg = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pOrg" );
  src.pCalc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCalc" );
  src.pCorr = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCorr" );
  src.pCorrDE = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCorrDE" );
  src.xb = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xb" );
  src.yb = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yb" );
  src.ub = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ub" );
  src.vb = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vb" );
  src.xs = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xs" );
  src.ys = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ys" );
  src.us = new TTreeReaderValue<std::vector<Double_t>>( *reader, "us" );
  src.vs = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vs" );
  src.Kflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "Kflag" );

  src.nclTpc = new TTreeReaderValue<Int_t>( *reader, "nclTpc" );
  src.remain_nclTpc = new TTreeReaderValue<Int_t>( *reader, "remain_nclTpc" );
#if SaveRawData
  src.cluster_x = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_x" );
  src.cluster_y = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_y" );
  src.cluster_z = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_z" );
  src.cluster_de = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_de" );
  src.cluster_size = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_size" );
  src.cluster_layer = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_layer" );
  src.cluster_mrow = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_mrow" );
  src.cluster_de_center = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_de_center" );
  src.cluster_x_center = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_x_center" );
  src.cluster_y_center = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_y_center" );
  src.cluster_z_center = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_z_center" );
  src.cluster_row_center = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_row_center" );
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
  src.mom_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "mom_x" );
  src.mom_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "mom_y" );
  src.mom_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "mom_z" );
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
  src.track_cluster_size = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_size" );
  src.track_cluster_mrow = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_mrow" );
  src.track_cluster_de_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_de_center" );
  src.track_cluster_x_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_x_center" );
  src.track_cluster_y_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_y_center" );
  src.track_cluster_z_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_z_center" );
  src.track_cluster_row_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_row_center" );

  src.chargeIndistinguishable = new TTreeReaderValue<std::vector<Int_t>>( *reader, "chargeIndistinguishable" );
  src.pid_inverted = new TTreeReaderValue<std::vector<Int_t>>( *reader, "pid_inverted" );
  src.chisqr_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqr_inverted" );
  src.pval_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pval_inverted" );
  src.helix_cx_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_cx_inverted" );
  src.helix_cy_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_cy_inverted" );
  src.helix_z0_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_z0_inverted" );
  src.helix_r_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_r_inverted" );
  src.helix_dz_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_dz_inverted" );
  src.mom0_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "mom0_inverted" );

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

  src.isLambda = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isLambda" );
  src.ncombiLambda = new TTreeReaderValue<std::vector<Int_t>>( *reader, "ncombiLambda" );
  src.distLambda = new TTreeReaderValue<std::vector<Double_t>>( *reader, "distLambda" );
  src.angleLambda = new TTreeReaderValue<std::vector<Double_t>>( *reader, "angleLambda" );
  src.bestmassLambda = new TTreeReaderValue<std::vector<Double_t>>( *reader, "bestmassLambda" );
  src.massLambda = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "massLambda" );
  src.vtxLambda_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxLambda_x" );
  src.vtxLambda_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxLambda_y" );
  src.vtxLambda_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxLambda_z" );
  src.momLambda = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "momLambda" );
  src.momLambda_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "momLambda_x" );
  src.momLambda_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "momLambda_y" );
  src.momLambda_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "momLambda_z" );
  src.decaysidLambda = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "decaysidLambda" );
  src.decaysmomLambda = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "decaysmomLambda" );
  src.decaysmomLambda_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "decaysmomLambda_x" );
  src.decaysmomLambda_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "decaysmomLambda_y" );
  src.decaysmomLambda_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "decaysmomLambda_z" );

  src.nvtxTpcClustered = new TTreeReaderValue<Int_t>(*reader,"nvtxTpcClustered");
  src.clusteredVtx_x = new TTreeReaderValue<std::vector<Double_t>>( *reader, "clusteredVtx_x" );
  src.clusteredVtx_y = new TTreeReaderValue<std::vector<Double_t>>( *reader, "clusteredVtx_y" );
  src.clusteredVtx_z = new TTreeReaderValue<std::vector<Double_t>>( *reader, "clusteredVtx_z" );
  src.clusteredVtxid = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "clusteredVtxid" );

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
  src.thetaCMTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaCMTPC" );
  src.costCMTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "costCMTPC" );
  src.xbTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xbTPC" );
  src.ybTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ybTPC" );
  src.ubTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ubTPC" );
  src.vbTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vbTPC" );
  src.xsTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xsTPC" );
  src.ysTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ysTPC" );
  src.usTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "usTPC" );
  src.vsTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vsTPC" );

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

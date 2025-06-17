// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <TLorentzVector.h>
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
#include "DCGeomMan.hh"
#include "DCDriftParamMan.hh"
#include "DCParameters.hh"
#include "DCTdcCalibMan.hh"
#include "DCLTrackHit.hh"
#include "DCHit.hh"
#include "DstHelper.hh"
#include "HodoPHCMan.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "NuclearMass.hh"
#include "RootHelper.hh"
#include "LorentzVector.hh"
#include "TPCAnalyzer.hh"
#include "TPCCluster.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCVertex.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "TPCRKTrack.hh"
#include "UserParamMan.hh"

#define SaveHistograms 0
#define RawHit 0
#define RawCluster 1
#define TrackClusterHist 0
#define TruncatedMean 0
#define TrackSearchFailed 0

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
const auto& gGeom  = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& zHSCenter = gGeom.LocalZ("HS");
//const auto& gPHC  = HodoPHCMan::GetInstance();
const Double_t truncatedMean = 0.8; //80%
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kTpcHit, kKuramaTracking, kK18HSTracking, kOutFile, nArgc
};
std::vector<TString> ArgName =
  { "[Process]", "[ConfFile]", "[TPCHit]", "[KuramaTracking]", "[K18HSTracking]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc", "kurama", "k18track","" };
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

  // K18
  Int_t ntK18;
  std::vector<Int_t> nhK18;
  std::vector<Double_t> chisqrK18;
  std::vector<Double_t> pHS;
  std::vector<Double_t> p_3rd;
  std::vector<Double_t> delta_3rd;
  std::vector<Double_t> xoutK18;
  std::vector<Double_t> youtK18;
  std::vector<Double_t> uoutK18;
  std::vector<Double_t> voutK18;
  std::vector<std::vector<Double_t>> layerK18;
  std::vector<std::vector<Int_t>> layerK18Int;
  std::vector<std::vector<Double_t>> wireK18;
  std::vector<std::vector<Double_t>> localhitposK18;
  std::vector<std::vector<Double_t>> wposK18;
  std::vector<Double_t> xtgtHS;
  std::vector<Double_t> ytgtHS;
  std::vector<Double_t> ztgtHS;
  std::vector<Double_t> utgtHS;
  std::vector<Double_t> vtgtHS;

  std::vector<Double_t> xbh2HS;
  std::vector<Double_t> ybh2HS;
  std::vector<Double_t> zbh2HS;
  std::vector<Double_t> ubh2HS;
  std::vector<Double_t> vbh2HS;
  std::vector<Double_t> pbh2HS;

  std::vector<Double_t> xgasvesselHS;
  std::vector<Double_t> ygasvesselHS;
  std::vector<Double_t> zgasvesselHS;
  std::vector<Double_t> ugasvesselHS;
  std::vector<Double_t> vgasvesselHS;
  std::vector<Double_t> pgasvesselHS;

  std::vector<Double_t> xhtofHS;
  std::vector<Double_t> yhtofHS;
  std::vector<Double_t> zhtofHS;
  std::vector<Double_t> uhtofHS;
  std::vector<Double_t> vhtofHS;
  std::vector<Double_t> initmomHS;
  std::vector<std::vector<Double_t>> xvpHS;
  std::vector<std::vector<Double_t>> yvpHS;
  std::vector<std::vector<Double_t>> zvpHS;
  std::vector<std::vector<Double_t>> uvpHS;
  std::vector<std::vector<Double_t>> vvpHS;
  std::vector<std::vector<Double_t>> xbcHS;
  std::vector<std::vector<Double_t>> ybcHS;
  std::vector<std::vector<Double_t>> zbcHS;
  std::vector<std::vector<Double_t>> ubcHS;
  std::vector<std::vector<Double_t>> vbcHS;

  //Kurama
  Int_t ntKurama;
  std::vector<Double_t> chisqrKurama;
  std::vector<Double_t> pKurama;
  std::vector<Double_t> qKurama;
  std::vector<Double_t> xtgtKurama;
  std::vector<Double_t> ytgtKurama;
  std::vector<Double_t> utgtKurama;
  std::vector<Double_t> vtgtKurama;
  std::vector<Double_t> thetaKurama;
  std::vector<Double_t> phiKurama;
  std::vector<Double_t> m2;
  std::vector<std::vector<Double_t>> xvpKurama;
  std::vector<std::vector<Double_t>> yvpKurama;
  std::vector<std::vector<Double_t>> zvpKurama;
  std::vector<Double_t> xhtofKurama;
  std::vector<Double_t> yhtofKurama;
  std::vector<Double_t> zhtofKurama;
  // For TPC + Kurama RK
  std::vector<Int_t> nh;
  std::vector<Double_t> xout;
  std::vector<Double_t> yout;
  std::vector<Double_t> zout;
  std::vector<Double_t> pxout;
  std::vector<Double_t> pyout;
  std::vector<Double_t> pzout;
  std::vector<Double_t> xtof;
  std::vector<Double_t> ytof;
  std::vector<Double_t> ztof;
  std::vector<Double_t> pxtof;
  std::vector<Double_t> pytof;
  std::vector<Double_t> pztof;
  std::vector<std::vector<Double_t>> layer;
  std::vector<std::vector<Int_t>> layerInt;
  std::vector<std::vector<Double_t>> wire;
  std::vector<std::vector<Double_t>> localhitpos;
  std::vector<std::vector<Double_t>> wpos;

  Int_t nhTpc;
  std::vector<Double_t> raw_hitpos_x;
  std::vector<Double_t> raw_hitpos_y;
  std::vector<Double_t> raw_hitpos_z;
  std::vector<Double_t> raw_de;
  std::vector<Int_t> raw_padid;
  std::vector<Int_t> raw_layer;
  std::vector<Int_t> raw_row;

  Int_t nclTpc;
  Int_t remain_nclTpc;
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

  Int_t ntTpc; // Number of tracks
  Int_t ntKuramaCandidate; //Numer of tracks which are kurama track candidates(before TPCKurama tracking)
  std::vector<Int_t> isKuramaCandidate;
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Int_t> nhtrackEff; // Number of Hits actually used in tracking.
  std::vector<Int_t> trackid; //for Kurama & K1.8 tracks
  std::vector<Int_t> isBeam;
  std::vector<Int_t> isXi;
  std::vector<Int_t> isKurama;
  std::vector<Int_t> isK18;
  std::vector<Int_t> isAccidental;
  std::vector<Int_t> isMultiloop;
  std::vector<Int_t> flag;
  std::vector<Int_t> fittime;  //usec
  std::vector<Int_t> searchtime; //usec
  std::vector<Int_t> niteration; //usec
  std::vector<Double_t> chisqr;
  std::vector<Double_t> pval;
  std::vector<Double_t> distTgt;
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

  std::vector<Double_t> dz_factor;
  std::vector<Double_t> mom0_x;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_y;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_z;//Helix momentum at Y = 0
  std::vector<Double_t> mom0;//Helix momentum at Y = 0
  std::vector<Int_t> charge;//Helix charge
  std::vector<Double_t> path;//Helix path

  std::vector<Int_t> pid;
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
  std::vector<std::vector<Double_t>> residual_t;
  std::vector<std::vector<Double_t>> residual_x;
  std::vector<std::vector<Double_t>> residual_y;
  std::vector<std::vector<Double_t>> residual_z;
  std::vector<std::vector<Double_t>> resolution_x;
  std::vector<std::vector<Double_t>> resolution_y;
  std::vector<std::vector<Double_t>> resolution_z;
  std::vector<std::vector<Double_t>> pull;
  std::vector<std::vector<Double_t>> helix_t;
  std::vector<std::vector<Double_t>> pathhit;
  std::vector<std::vector<Double_t>> alpha;
  std::vector<std::vector<Double_t>> houghflag;
  std::vector<std::vector<Double_t>> track_cluster_de;
  std::vector<std::vector<Double_t>> track_cluster_size;
  std::vector<std::vector<Double_t>> track_cluster_mrow;
  std::vector<std::vector<Double_t>> track_cluster_de_center;
  std::vector<std::vector<Double_t>> track_cluster_x_center;
  std::vector<std::vector<Double_t>> track_cluster_y_center;
  std::vector<std::vector<Double_t>> track_cluster_z_center;
  std::vector<std::vector<Double_t>> track_cluster_row_center;
  //Exclusive tracking
  std::vector<std::vector<Double_t>> exresidual_t;
  std::vector<std::vector<Double_t>> exresidual_x;
  std::vector<std::vector<Double_t>> exresidual_y;
  std::vector<std::vector<Double_t>> exresidual_z;
  //Geometric mean of inclusive & exclusive residual
  std::vector<std::vector<Double_t>> intrinsic_residual_t;
  std::vector<std::vector<Double_t>> intrinsic_residual_x;
  std::vector<std::vector<Double_t>> intrinsic_residual_y;
  std::vector<std::vector<Double_t>> intrinsic_residual_z;

  //Inverted charge track
  std::vector<Int_t> chargeIndistinguishable;
  std::vector<Double_t> chisqr_inverted;
  std::vector<Double_t> pval_inverted;
  std::vector<Double_t> helix_cx_inverted;
  std::vector<Double_t> helix_cy_inverted;
  std::vector<Double_t> helix_z0_inverted;
  std::vector<Double_t> helix_r_inverted;
  std::vector<Double_t> helix_dz_inverted;
  std::vector<Double_t> mom0_x_inverted;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_y_inverted;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_z_inverted;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_inverted;//Helix momentum at Y = 0
  std::vector<Int_t> pid_inverted;

  Int_t vpntTpc; // Number of Tracks
  std::vector<Int_t> vpnhtrack;
  std::vector<Int_t> vptrackid; //for Kurama K1.8 tracks
  std::vector<Int_t> vpisKurama; // isKurama: 1 = Beam, 0 = Scat
  std::vector<Int_t> vpisK18;
  std::vector<Double_t> vphelix_cx;
  std::vector<Double_t> vphelix_cy;
  std::vector<Double_t> vphelix_z0;
  std::vector<Double_t> vphelix_r;
  std::vector<Double_t> vphelix_dz;
  std::vector<Double_t> vpmom0;//Helix momentum at Y = 0
  std::vector<Int_t> vpcharge;//Helix charge
  std::vector<std::vector<Double_t>> vphelix_t;
  std::vector<std::vector<Double_t>> vppos_x;
  std::vector<std::vector<Double_t>> vppos_y;
  std::vector<std::vector<Double_t>> vppos_z;
  std::vector<std::vector<Double_t>> residual_vppos_x;
  std::vector<std::vector<Double_t>> residual_vppos_y;
  std::vector<std::vector<Double_t>> residual_vppos_z;

  Int_t failed_ntTpc; // Number of Tracks
  std::vector<Int_t> failed_nhtrack;
  std::vector<Int_t> failed_trackid; //for Kurama K1.8 tracks
  std::vector<Int_t> failed_isBeam; // isBeam: 1 = Beam, 0 = Scat
  std::vector<Int_t> failed_isKurama; // isKurama: 1 = Beam, 0 = Scat
  std::vector<Int_t> failed_isK18;
  std::vector<Int_t> failed_nclbeforetgt;
  std::vector<Int_t> failed_isAccidental;
  std::vector<Int_t> failed_flag;
  std::vector<Int_t> failed_fittime; //usec
  std::vector<Int_t> failed_searchtime; //usec
  std::vector<Int_t> failed_niteration; //usec
  std::vector<Double_t> failed_helix_cx;
  std::vector<Double_t> failed_helix_cy;
  std::vector<Double_t> failed_helix_z0;
  std::vector<Double_t> failed_helix_r;
  std::vector<Double_t> failed_helix_dz;
  std::vector<Double_t> failed_mom0;//Helix momentum at Y = 0
  std::vector<Int_t> failed_charge;//Helix charge

  std::vector<std::vector<Double_t>> failed_hitlayer;
  std::vector<std::vector<Double_t>> failed_hitpos_x;
  std::vector<std::vector<Double_t>> failed_hitpos_y;
  std::vector<std::vector<Double_t>> failed_hitpos_z;
  std::vector<std::vector<Double_t>> failed_calpos_x;
  std::vector<std::vector<Double_t>> failed_calpos_y;
  std::vector<std::vector<Double_t>> failed_calpos_z;
  std::vector<std::vector<Double_t>> failed_helix_t;
  std::vector<std::vector<Double_t>> failed_residual;
  std::vector<std::vector<Double_t>> failed_residual_x;
  std::vector<std::vector<Double_t>> failed_residual_y;
  std::vector<std::vector<Double_t>> failed_residual_z;
  std::vector<std::vector<Double_t>> failed_track_cluster_de;
  std::vector<std::vector<Double_t>> failed_track_cluster_size;
  std::vector<std::vector<Double_t>> failed_track_cluster_mrow;

  std::vector<Int_t> isgoodTPCK18;
  std::vector<Int_t> tpcidTPCK18;
  std::vector<Int_t> niterationTPCK18;
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
  std::vector<Int_t> niterationTPCKurama;
  std::vector<Int_t> kflagTPCKurama;
  std::vector<Int_t> pflagTPCKurama;
  std::vector<Double_t> chisqrTPCKurama;
  std::vector<Double_t> pTPCKurama;
  std::vector<Double_t> qTPCKurama;
  std::vector<Double_t> m2TPCKurama;
  std::vector<Double_t> m2OrgTPCKurama;
  std::vector<Double_t> xtgtTPCKurama;
  std::vector<Double_t> ytgtTPCKurama;
  std::vector<Double_t> utgtTPCKurama;
  std::vector<Double_t> vtgtTPCKurama;
  std::vector<Double_t> thetaTPCKurama;
  std::vector<Double_t> pathTPCKurama;
  std::vector<Double_t> lhtofTPCKurama;
  std::vector<Double_t> xhtofTPCKurama;
  std::vector<Double_t> yhtofTPCKurama;
  std::vector<Double_t> phtofTPCKurama;
  std::vector<Double_t> lgasvesselTPCKurama;
  std::vector<Double_t> xgasvesselTPCKurama;
  std::vector<Double_t> ygasvesselTPCKurama;
  std::vector<Double_t> pgasvesselTPCKurama;
  std::vector<std::vector<Double_t>> lvpTPCKurama;
  std::vector<std::vector<Double_t>> xvpTPCKurama;
  std::vector<std::vector<Double_t>> yvpTPCKurama;

  std::vector<Int_t> isgoodTPC;
  std::vector<Int_t> insideTPC;
  std::vector<Double_t> vtxTPC;
  std::vector<Double_t> vtyTPC;
  std::vector<Double_t> vtzTPC;
  std::vector<Double_t> pxKmTPC;
  std::vector<Double_t> pyKmTPC;
  std::vector<Double_t> pzKmTPC;
  std::vector<Double_t> pxScatTPC;
  std::vector<Double_t> pyScatTPC;
  std::vector<Double_t> pzScatTPC;

  std::vector<Double_t> closeDistTPC;
  std::vector<Double_t> MissMassTPC;
  std::vector<Double_t> MissMassCorrTPC;
  std::vector<Double_t> MissMassCorrDETPC;
  std::vector<Double_t> MissMassNuclTPC;
  std::vector<Double_t> MissMassNuclCorrTPC;
  std::vector<Double_t> MissMassNuclCorrDETPC;
  std::vector<Double_t> pOrgTPC;
  std::vector<Double_t> pCorrTPC;
  std::vector<Double_t> pCorrDETPC;
  std::vector<Double_t> pCalcTPC;
  std::vector<Double_t> thetaCMTPC;
  std::vector<Double_t> costCMTPC;
  std::vector<Double_t> pCalcDETPC;
  std::vector<Double_t> thetaCMDETPC;
  std::vector<Double_t> costCMDETPC;
  std::vector<Double_t> xistarpCalcDETPC;
  std::vector<Double_t> xistarthetaCMDETPC;
  std::vector<Double_t> xistarcostCMDETPC;
  std::vector<Double_t> kpscatpCalcTPC;
  std::vector<Double_t> kpscatthetaCMTPC;
  std::vector<Double_t> kpscatcostCMTPC;
  std::vector<Double_t> kpscatpCalcDETPC;
  std::vector<Double_t> kpscatthetaCMDETPC;
  std::vector<Double_t> kpscatcostCMDETPC;
  std::vector<Double_t> thetaTPC;
  std::vector<Double_t> ubTPC;
  std::vector<Double_t> vbTPC;
  std::vector<Double_t> usTPC;
  std::vector<Double_t> vsTPC;

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
  std::vector<Double_t> Clusteredvtx_x;
  std::vector<Double_t> Clusteredvtx_y;
  std::vector<Double_t> Clusteredvtx_z;
  std::vector<std::vector<Double_t>> Clusteredvtxid;

  void clear( void )
  {
    trigpat.clear();
    trigflag.clear();
    clkTpc.clear();

    runnum = 0;
    evnum = 0;
    status = 0;

    ntK18 = 0;
    nhK18.clear();
    chisqrK18.clear();
    p_3rd.clear();
    pHS.clear();
    delta_3rd.clear();
    xoutK18.clear();
    youtK18.clear();
    uoutK18.clear();
    voutK18.clear();
    layerK18.clear();
    layerK18Int.clear();
    wireK18.clear();
    localhitposK18.clear();
    wposK18.clear();

    initmomHS.clear();
    xtgtHS.clear();
    ytgtHS.clear();
    ztgtHS.clear();
    utgtHS.clear();
    vtgtHS.clear();
    xbh2HS.clear();
    ybh2HS.clear();
    zbh2HS.clear();
    ubh2HS.clear();
    vbh2HS.clear();
    pbh2HS.clear();

    xgasvesselHS.clear();
    ygasvesselHS.clear();
    zgasvesselHS.clear();
    ugasvesselHS.clear();
    vgasvesselHS.clear();
    pgasvesselHS.clear();

    xhtofHS.clear();
    yhtofHS.clear();
    zhtofHS.clear();
    uhtofHS.clear();
    vhtofHS.clear();
    xvpHS.clear();
    yvpHS.clear();
    uvpHS.clear();
    vvpHS.clear();
    zvpHS.clear();
    xbcHS.clear();
    ybcHS.clear();
    zbcHS.clear();
    ubcHS.clear();
    vbcHS.clear();

    ntKurama = 0;
    chisqrKurama.clear();
    pKurama.clear();
    qKurama.clear();
    m2.clear();
    xtgtKurama.clear();
    ytgtKurama.clear();
    utgtKurama.clear();
    vtgtKurama.clear();
    thetaKurama.clear();
    phiKurama.clear();
    xvpKurama.clear();
    yvpKurama.clear();
    zvpKurama.clear();
    xhtofKurama.clear();
    yhtofKurama.clear();
    zhtofKurama.clear();

    nh.clear();
    xout.clear();
    yout.clear();
    zout.clear();
    pxout.clear();
    pyout.clear();
    pzout.clear();
    xtof.clear();
    ytof.clear();
    ztof.clear();
    pxtof.clear();
    pytof.clear();
    pztof.clear();
    layer.clear();
    layerInt.clear();
    wire.clear();
    localhitpos.clear();
    wpos.clear();

    nhTpc = 0;
    raw_hitpos_x.clear();
    raw_hitpos_y.clear();
    raw_hitpos_z.clear();
    raw_de.clear();
    raw_padid.clear();
    raw_layer.clear();
    raw_row.clear();

    nclTpc = 0;
    remain_nclTpc = 0;
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

    ntTpc = 0;
    ntKuramaCandidate = 0; //Numer of tracks which are kurama track candidates(before TPCKurama tracking)
    isKuramaCandidate.clear();
    nhtrack.clear();
    nhtrackEff.clear();
    trackid.clear();
    isBeam.clear();
    isXi.clear();
    isKurama.clear();
    isK18.clear();
    isAccidental.clear();
    isMultiloop.clear();
    flag.clear();
    fittime.clear();
    searchtime.clear();
    niteration.clear();
    chisqr.clear();
    pval.clear();
    distTgt.clear();
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
#endif
    dz_factor.clear();

    mom0_x.clear();
    mom0_y.clear();
    mom0_z.clear();
    mom0.clear();

    charge.clear();
    path.clear();

    pid.clear();
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
    residual_t.clear();
    residual_x.clear();
    residual_y.clear();
    residual_z.clear();
    resolution_x.clear();
    resolution_y.clear();
    resolution_z.clear();
    pull.clear();
    helix_t.clear();
    pathhit.clear();
    alpha.clear();
    houghflag.clear();
    track_cluster_de.clear();
    track_cluster_size.clear();
    track_cluster_mrow.clear();
    track_cluster_de_center.clear();
    track_cluster_x_center.clear();
    track_cluster_y_center.clear();
    track_cluster_z_center.clear();
    track_cluster_row_center.clear();

    exresidual_t.clear();
    exresidual_x.clear();
    exresidual_y.clear();
    exresidual_z.clear();
    intrinsic_residual_t.clear();
    intrinsic_residual_x.clear();
    intrinsic_residual_y.clear();
    intrinsic_residual_z.clear();

    chargeIndistinguishable.clear();
    chisqr_inverted.clear();
    pval_inverted.clear();
    helix_cx_inverted.clear();
    helix_cy_inverted.clear();
    helix_z0_inverted.clear();
    helix_r_inverted.clear();
    helix_dz_inverted.clear();
    mom0_x_inverted.clear();
    mom0_y_inverted.clear();
    mom0_z_inverted.clear();
    mom0_inverted.clear();
    pid_inverted.clear();

    vpntTpc = 0;
    vpnhtrack.clear();
    vptrackid.clear();
    vpisKurama.clear();
    vpisK18.clear();
    vphelix_cx.clear();
    vphelix_cy.clear();
    vphelix_z0.clear();
    vphelix_r.clear();
    vphelix_dz.clear();
    vpmom0.clear();
    vpcharge.clear();
    vphelix_t.clear();
    vppos_x.clear();
    vppos_y.clear();
    vppos_z.clear();
    residual_vppos_x.clear();
    residual_vppos_y.clear();
    residual_vppos_z.clear();

    failed_ntTpc = 0;
    failed_nhtrack.clear();
    failed_trackid.clear();
    failed_isBeam.clear();
    failed_isKurama.clear();
    failed_isK18.clear();
    failed_nclbeforetgt.clear();
    failed_isAccidental.clear();
    failed_flag.clear();
    failed_fittime.clear();
    failed_searchtime.clear();
    failed_niteration.clear();

    failed_helix_cx.clear();
    failed_helix_cy.clear();
    failed_helix_z0.clear();
    failed_helix_r.clear();
    failed_helix_dz.clear();
    failed_mom0.clear();
    failed_charge.clear();

    failed_hitlayer.clear();
    failed_hitpos_x.clear();
    failed_hitpos_y.clear();
    failed_hitpos_z.clear();
    failed_calpos_x.clear();
    failed_calpos_y.clear();
    failed_calpos_z.clear();
    failed_helix_t.clear();
    failed_residual.clear();
    failed_residual_x.clear();
    failed_residual_y.clear();
    failed_residual_z.clear();
    failed_track_cluster_de.clear();
    failed_track_cluster_size.clear();
    failed_track_cluster_mrow.clear();

    isgoodTPCK18.clear();
    tpcidTPCK18.clear();
    niterationTPCK18.clear();
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
    niterationTPCKurama.clear();
    kflagTPCKurama.clear();
    pflagTPCKurama.clear();
    chisqrTPCKurama.clear();
    pTPCKurama.clear();
    qTPCKurama.clear();
    m2TPCKurama.clear();
    m2OrgTPCKurama.clear();
    xtgtTPCKurama.clear();
    ytgtTPCKurama.clear();
    utgtTPCKurama.clear();
    vtgtTPCKurama.clear();
    thetaTPCKurama.clear();
    pathTPCKurama.clear();
    lhtofTPCKurama.clear();
    xhtofTPCKurama.clear();
    yhtofTPCKurama.clear();
    phtofTPCKurama.clear();
    lgasvesselTPCKurama.clear();
    xgasvesselTPCKurama.clear();
    ygasvesselTPCKurama.clear();
    pgasvesselTPCKurama.clear();
    lvpTPCKurama.clear();
    xvpTPCKurama.clear();
    yvpTPCKurama.clear();

    isgoodTPC.clear();
    insideTPC.clear();
    vtxTPC.clear();
    vtyTPC.clear();
    vtzTPC.clear();
    pxKmTPC.clear();
    pyKmTPC.clear();
    pzKmTPC.clear();
    pxScatTPC.clear();
    pyScatTPC.clear();
    pzScatTPC.clear();

    closeDistTPC.clear();
    MissMassTPC.clear();
    MissMassCorrTPC.clear();
    MissMassCorrDETPC.clear();
    MissMassNuclTPC.clear();
    MissMassNuclCorrTPC.clear();
    MissMassNuclCorrDETPC.clear();
    pOrgTPC.clear();
    pCorrTPC.clear();
    pCorrDETPC.clear();
    pCalcTPC.clear();
    thetaCMTPC.clear();
    costCMTPC.clear();
    pCalcDETPC.clear();
    thetaCMDETPC.clear();
    costCMDETPC.clear();
    xistarpCalcDETPC.clear();
    xistarthetaCMDETPC.clear();
    xistarcostCMDETPC.clear();
    kpscatpCalcTPC.clear();
    kpscatthetaCMTPC.clear();
    kpscatcostCMTPC.clear();
    kpscatpCalcDETPC.clear();
    kpscatthetaCMDETPC.clear();
    kpscatcostCMDETPC.clear();
    thetaTPC.clear();

    ubTPC.clear();
    vbTPC.clear();
    usTPC.clear();
    vsTPC.clear();

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
    Clusteredvtx_x.clear();
    Clusteredvtx_y.clear();
    Clusteredvtx_z.clear();
    Clusteredvtxid.clear();

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
  TTreeReaderValue<std::vector<Double_t>>* cdeTpc;     // cdE
  TTreeReaderValue<std::vector<Double_t>>* tTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* ctTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* chisqrTpc; // chi^2 of signal fitting
  TTreeReaderValue<std::vector<Double_t>>* clkTpc;      //clock time

  // K18HS
  Int_t ntK18;
  Int_t nhK18[MaxHits];
  Double_t chisqrK18[MaxHits];
  Double_t pHS[MaxHits];
  Double_t p_3rd[MaxHits];
  Double_t delta_3rd[MaxHits];
  Double_t xoutK18[MaxHits];
  Double_t youtK18[MaxHits];
  Double_t uoutK18[MaxHits];
  Double_t voutK18[MaxHits];
  Double_t layerK18[MaxHits][NumOfLayersBcOut];
  Double_t wireK18[MaxHits][NumOfLayersBcOut];
  Double_t localhitposK18[MaxHits][NumOfLayersBcOut];
  Double_t wposK18[MaxHits][NumOfLayersBcOut];
  Double_t xtgtHS[MaxHits];
  Double_t ytgtHS[MaxHits];
  Double_t ztgtHS[MaxHits];
  Double_t utgtHS[MaxHits];
  Double_t vtgtHS[MaxHits];
  Double_t initmomHS[MaxHits];

  Double_t xbh2HS[MaxHits];
  Double_t ybh2HS[MaxHits];
  Double_t zbh2HS[MaxHits];
  Double_t ubh2HS[MaxHits];
  Double_t vbh2HS[MaxHits];
  Double_t pbh2HS[MaxHits];

  Double_t xgasvesselHS[MaxHits];
  Double_t ygasvesselHS[MaxHits];
  Double_t zgasvesselHS[MaxHits];
  Double_t ugasvesselHS[MaxHits];
  Double_t vgasvesselHS[MaxHits];
  Double_t pgasvesselHS[MaxHits];

  Double_t xvp1HS[MaxHits];
  Double_t yvp1HS[MaxHits];
  Double_t zvp1HS[MaxHits];
  Double_t uvp1HS[MaxHits];
  Double_t vvp1HS[MaxHits];

  Double_t xvp2HS[MaxHits];
  Double_t yvp2HS[MaxHits];
  Double_t zvp2HS[MaxHits];
  Double_t uvp2HS[MaxHits];
  Double_t vvp2HS[MaxHits];

  Double_t xvp3HS[MaxHits];
  Double_t yvp3HS[MaxHits];
  Double_t zvp3HS[MaxHits];
  Double_t uvp3HS[MaxHits];
  Double_t vvp3HS[MaxHits];

  Double_t xvp4HS[MaxHits];
  Double_t yvp4HS[MaxHits];
  Double_t zvp4HS[MaxHits];
  Double_t uvp4HS[MaxHits];
  Double_t vvp4HS[MaxHits];

  Double_t xhtofHS[MaxHits];
  Double_t yhtofHS[MaxHits];
  Double_t zhtofHS[MaxHits];
  Double_t uhtofHS[MaxHits];
  Double_t vhtofHS[MaxHits];

  Double_t xbcHS[MaxHits][NumOfLayersBcOut];
  Double_t ybcHS[MaxHits][NumOfLayersBcOut];
  Double_t zbcHS[MaxHits][NumOfLayersBcOut];
  Double_t ubcHS[MaxHits][NumOfLayersBcOut];
  Double_t vbcHS[MaxHits][NumOfLayersBcOut];

  //Kurama
  Int_t ntKurama;
  Double_t chisqrKurama[MaxHits];
  Double_t m2[MaxHits];
  Double_t pKurama[MaxHits];
  Double_t qKurama[MaxHits];
  Double_t xtgtKurama[MaxHits];
  Double_t ytgtKurama[MaxHits];
  Double_t utgtKurama[MaxHits];
  Double_t vtgtKurama[MaxHits];
  Double_t xtofKurama[MaxHits];
  Double_t ytofKurama[MaxHits];
  Double_t utofKurama[MaxHits];
  Double_t vtofKurama[MaxHits];
  Double_t tofsegKurama[MaxHits];
  Double_t thetaKurama[MaxHits];
  Double_t phiKurama[MaxHits];
  Double_t path[MaxHits];
  Double_t stof[MaxHits];

  Double_t vpxtpc[MaxHits][NumOfLayersVP];
  Double_t vpytpc[MaxHits][NumOfLayersVP];
  Double_t vpztpc[MaxHits][NumOfLayersVP];
  Double_t vpxhtof[MaxHits];
  Double_t vpyhtof[MaxHits];
  Double_t vpzhtof[MaxHits];

  Int_t nh[MaxHits];
  Double_t xout[MaxHits];
  Double_t yout[MaxHits];
  Double_t zout[MaxHits];
  Double_t pxout[MaxHits];
  Double_t pyout[MaxHits];
  Double_t pzout[MaxHits];
  Double_t xtof[MaxHits];
  Double_t ytof[MaxHits];
  Double_t ztof[MaxHits];
  Double_t pxtof[MaxHits];
  Double_t pytof[MaxHits];
  Double_t pztof[MaxHits];
  Int_t layer[MaxHits][22];
  Double_t wire[MaxHits][22];
  Double_t localhitpos[MaxHits][22];
  Double_t wpos[MaxHits][22];

};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
  enum eDetHid {
    TPCK18VPHid = 100000,
    TPCKuramaVPHid = 200000,
    TPCK18RKHid = 300000,
    TPCKuramaRKHid = 400000,
    TPCClHid = 500000,
    TPCInclusiveHid = 600000,
    TPCExclusiveHid = 700000,
    TPCIntrinsicHid = 800000
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

const Double_t pK18_factor = 1.000;
const Double_t pK18_offset = 0.000;
const Double_t xb_off = 0.000;
const Double_t yb_off = 0.000;
const Double_t ub_off = 0.000;
const Double_t vb_off = 0.000;

Double_t pKuramaCorrection(Double_t u, Double_t v, Double_t pOrg)
{
  //Angular dependence correction
  Double_t pCorr = pOrg;
  //#if 1 // w/o Z + 6 mm correction
#if 0
  Double_t par_dxdz[4] = {0.00167653, 0.0329239, 0.39883, 1.37901};
  Double_t par_dydz[5] = {-0.00545227, -0.0424615, 1.01906, 1.25306, -11.627};

  if(u<-0.275) pCorr -=  -0.00589507;
  else if(u>0.025) pCorr -= 0.00277044;
  else{
    for(Int_t i=0; i<4; ++i) pCorr -= par_dxdz[i]*TMath::Power(u, i);
  }

  if(v<-0.135) pCorr -= 0.0119074;
  else if(v>0.125) pCorr -= 0.00477157;
  else{
    for(Int_t i=0; i<5; ++i) pCorr -= par_dydz[i]*TMath::Power(v, i);
  }
#elif 0 // w/ Z + 6 mm correction before TPC resolution debugging
  Double_t par_dxdz[4] = {0.00598859, 0.184888, 1.91559, 6.0048};
  Double_t par_dydz[5] = {-0.004968, -0.0187119, 0.956009, -0.620882, -10.637};

  if(u<-0.275) pCorr -= -0.02487;
  else if(u>0.025) pCorr -= 0.0119018;
  else{
    for(Int_t i=0; i<4; ++i) pCorr -= par_dxdz[i]*TMath::Power(u, i);
  }

  if(v<-0.135) pCorr -= 0.0129759;
  else if(v>0.125) pCorr -= 0.00382108;
  else{
    for(Int_t i=0; i<5; ++i) pCorr -= par_dydz[i]*TMath::Power(v, i);
  }
#else // w/ Z + 6 mm correction after TPC resolution debugging
  Double_t par_dxdz[4] = {0.00783981, 0.0858266, 0.511154, 2.10857};
  Double_t par_dydz[5] = {-0.00593035, -0.0209602, 1.1396, -0.269324, -23.126};

  if(u<-0.275) pCorr -= -0.0209581;
  else if(u>0.025) pCorr -= 0.0103379;
  else{
    for(Int_t i=0; i<4; ++i) pCorr -= par_dxdz[i]*TMath::Power(u, i);
  }

  if(v<-0.135) pCorr -= 0.0106498;
  else if(v>0.125) pCorr -= 0.00308382;
  else{
    for(Int_t i=0; i<5; ++i) pCorr -= par_dydz[i]*TMath::Power(v, i);
  }
#endif

  //momentum correction
  //P_measured = p2*P^2 + p1*P + p0
  //#if 1 // w/o Z + 6 mm correction
#if 0
  Double_t p0 = 0.00816054; Double_t p1 = 1.00456; Double_t p2 = -2.50295e-05;
  if(p1*p1+4.*(pCorr-p0)*p2<0) pCorr = 0.001;
  else pCorr = 0.5*(-p1+TMath::Sqrt(p1*p1+4.*(pCorr-p0)*p2))/p2;
#elif 0 // w/ Z + 6 mm correction
  Double_t p0 = 0.0412952; Double_t p1 = 0.955566; Double_t p2 = 0.017896;
  if(p1*p1+4.*(pCorr-p0)*p2<0) pCorr = 0.001;
  else pCorr = 0.5*(-p1+TMath::Sqrt(p1*p1+4.*(pCorr-p0)*p2))/p2;
#else // w/ Z + 6 mm correction after TPC resolution debugging
  Double_t p0 = 0.0463731; Double_t p1 = 0.955957; Double_t p2 = 0.0178651;
  if(p1*p1+4.*(pCorr-p0)*p2<0) pCorr = 0.001;
  else pCorr = 0.5*(-p1+TMath::Sqrt(p1*p1+4.*(pCorr-p0)*p2))/p2;
#endif
  return pCorr;
}

Double_t pKuramaCorrection_proton(Double_t u, Double_t v, Double_t pOrg)
{
  //Angular dependence correction
  Double_t pCorr = pOrg;
  //#if 1 // w/o Z + 6 mm correction
#if 0
  Double_t par_dxdz[4] = {0.00331135, -0.0838977, -0.508636, -0.516439};
  Double_t par_dydz[5] = {-0.0121156, -0.051397, 3.82854, 0.702847, -142.941};

  if(u<-0.245) pCorr -= 0.000930233;
  else if(u>0.105) pCorr -=  -0.0117035;
  else{
    for(Int_t i=0; i<4; ++i) pCorr -= par_dxdz[i]*TMath::Power(u, i);
  }

  if(v<-0.135) pCorr -= 0.0153907;
  else if(v>0.125) pCorr -= 0.00775557;
  else{
    for(Int_t i=0; i<5; ++i) pCorr -= par_dydz[i]*TMath::Power(v, i);
  }
#elif 0 // w/ Z + 6 mm correction before TPC resolution debugging
  Double_t par_dxdz[4] = {0.00889138, 0.0178901, -0.540299, -0.456336};
  Double_t par_dydz[5] = {-0.0137751, -0.0461785, 3.77421, -0.407101, -139.886};

  if(u<-0.245) pCorr -= -0.0212122;
  else if(u>0.105) pCorr -= 0.00428478;
  else{
    for(Int_t i=0; i<4; ++i) pCorr -= par_dxdz[i]*TMath::Power(u, i);
  }

  if(v<-0.135) pCorr -= 0.0157823;
  else if(v>0.125) pCorr -= 0.0044776;
  else{
    for(Int_t i=0; i<5; ++i) pCorr -= par_dydz[i]*TMath::Power(v, i);
  }
#else // w/ Z + 6 mm correction after TPC resolution debugging
  Double_t par_dxdz[4] = {0.00961677, 0.0628218, -0.251013, 0.262886};
  Double_t par_dydz[5] = {-0.0134842, -0.053486, 3.89929, -0.230734, -147.784};

  if(u<-0.245) pCorr -= -0.0247077;
  else if(u>0.105) pCorr -= 0.01375;
  else{
    for(Int_t i=0; i<4; ++i) pCorr -= par_dxdz[i]*TMath::Power(u, i);
  }

  if(v<-0.135) pCorr -= 0.0162822;
  else if(v>0.125) pCorr -= 0.00422581;
  else{
    for(Int_t i=0; i<5; ++i) pCorr -= par_dydz[i]*TMath::Power(v, i);
  }
#endif

  //momentum correction
  //P_measured = p2*P^2 + p1*P + p0
  //#if 1
#if 0 // w/o Z + 6 mm correction
  Double_t p0 = 0.00816054; Double_t p1 = 1.00456; Double_t p2 = -2.50295e-05;
  if(p1*p1+4.*(pCorr-p0)*p2<0) pCorr = 0.001;
  else pCorr = 0.5*(-p1+TMath::Sqrt(p1*p1+4.*(pCorr-p0)*p2))/p2;
#elif 0 // w/ Z + 6 mm correction before TPC resolution debugging
  Double_t p0 = 0.0412952; Double_t p1 = 0.955566; Double_t p2 = 0.017896;
  if(p1*p1+4.*(pCorr-p0)*p2<0) pCorr = 0.001;
  else pCorr = 0.5*(-p1+TMath::Sqrt(p1*p1+4.*(pCorr-p0)*p2))/p2;
#else // w/ Z + 6 mm correction after TPC resolution debugging
  Double_t p0 = 0.0463731; Double_t p1 = 0.955957; Double_t p2 = 0.0178651;
  if(p1*p1+4.*(pCorr-p0)*p2<0) pCorr = 0.001;
  else pCorr = 0.5*(-p1+TMath::Sqrt(p1*p1+4.*(pCorr-p0)*p2))/p2;
#endif
  return pCorr;
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
  auto start = std::chrono::high_resolution_clock::now();
  static const auto KaonMass    = pdg::KaonMass();
  static const auto PionMass    = pdg::PionMass();
  static const auto ProtonMass  = pdg::ProtonMass();
  static const auto LambdaMass  = pdg::LambdaMass();
  static const auto XiMass      = pdg::XiMinusMass();
  static const auto XiStarMass  = 1.5350;
  static const auto ElectronMass = pdg::ElectronMass();
  static const Double_t Carbon12Mass = 12.*TGeoUnit::amu_c2 - 6.*ElectronMass;
  static const Double_t Boron11Mass  = 11.009305167*TGeoUnit::amu_c2 - 5.*ElectronMass;
  static const auto MaxChisqrBcOut = gUser.GetParameter("MaxChisqrBcOut");
  static const auto MaxChisqrKurama = gUser.GetParameter("MaxChisqrKurama");
  static const Bool_t KKEvent = gUser.GetParameter("KKEvent");
  static const Bool_t KPEvent = gUser.GetParameter("KPEvent");
  static const Bool_t KHeavyEvent = gUser.GetParameter("KHeavyEvent");
  static const Bool_t ExclusiveTracking = gUser.GetParameter("ExclusiveTracking");
  static const auto xGlobalBcOut = gGeom.GetGlobalPosition("BC3-X1").X();
  static const auto yGlobalBcOut = gGeom.GetGlobalPosition("BC3-X1").Y();
  static const auto zGlobalBcOut = gGeom.GetGlobalPosition("BC3-X1").Z();
  static const auto zLocalBcOut = gGeom.GetLocalZ("BC3-X1");
  static const auto xGlobalSdcOut = gGeom.GetGlobalPosition("SDC4-X2").X();
  static const auto yGlobalSdcOut = gGeom.GetGlobalPosition("SDC4-X2").Y();
  const Int_t& IdTPCGasVessel_D = gGeom.DetectorId("VesselD");
  const Int_t& IdVPHTOF = gGeom.DetectorId("VPHTOF");

  GetEntry(ievent);

  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;
  event.clkTpc = **src.clkTpc;

  HF1( 1, event.status++ );

  if(src.ntKurama!=1 || src.ntK18!=1) return true;
  if(src.chisqrK18[0] > MaxChisqrBcOut || src.chisqrKurama[0] > MaxChisqrKurama) return true;

  if(KKEvent){
    if(src.m2[0] < 0.05 || src.m2[0] > 0.70) return true;
    if(src.qKurama[0] < 0 || src.pKurama[0] > 1.4) return true;
  }
  if(KPEvent){
    if(src.m2[0] < 0.50 || src.m2[0] > 1.40) return true;
    if(src.qKurama[0] < 0 || src.pKurama[0] < 0.0) return true;
  }
  if(KHeavyEvent){
    if(src.m2[0] < 2.50) return true;
    if(src.stof[0] > 45) return true;
    if(src.qKurama[0] < 0 || src.pKurama[0] < 0.0) return true;
  }

  if( ievent%1==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  HF1( 1, event.status++ );
  event.ntK18 = src.ntK18;
  std::vector<std::vector<TVector3>> vpK18;
  vpK18.resize( src.ntK18 );
  std::vector<TVector3> initPosK18;
  std::vector<TVector3> initMomK18;

  event.nhK18.resize(src.ntK18);
  event.chisqrK18.resize(src.ntK18);
  event.p_3rd.resize(src.ntK18);
  event.pHS.resize(src.ntK18);
  event.delta_3rd.resize(src.ntK18);
  event.xoutK18.resize(src.ntK18);
  event.youtK18.resize(src.ntK18);
  event.uoutK18.resize(src.ntK18);
  event.voutK18.resize(src.ntK18);
  event.layerK18.resize(src.ntK18);
  event.layerK18Int.resize(src.ntK18);
  event.wireK18.resize(src.ntK18);
  event.localhitposK18.resize(src.ntK18);
  event.wposK18.resize(src.ntK18);

  event.xtgtHS.resize(src.ntK18);
  event.ytgtHS.resize(src.ntK18);
  event.ztgtHS.resize(src.ntK18);
  event.utgtHS.resize(src.ntK18);
  event.vtgtHS.resize(src.ntK18);
  event.initmomHS.resize(src.ntK18);

  event.xbh2HS.resize(src.ntK18);
  event.ybh2HS.resize(src.ntK18);
  event.zbh2HS.resize(src.ntK18);
  event.ubh2HS.resize(src.ntK18);
  event.vbh2HS.resize(src.ntK18);
  event.pbh2HS.resize(src.ntK18);

  event.xgasvesselHS.resize(src.ntK18);
  event.ygasvesselHS.resize(src.ntK18);
  event.zgasvesselHS.resize(src.ntK18);
  event.ugasvesselHS.resize(src.ntK18);
  event.vgasvesselHS.resize(src.ntK18);
  event.pgasvesselHS.resize(src.ntK18);

  event.xhtofHS.resize(src.ntK18);
  event.yhtofHS.resize(src.ntK18);
  event.zhtofHS.resize(src.ntK18);
  event.uhtofHS.resize(src.ntK18);
  event.vhtofHS.resize(src.ntK18);

  event.xbcHS.resize(src.ntK18);
  event.ybcHS.resize(src.ntK18);
  event.zbcHS.resize(src.ntK18);
  event.ubcHS.resize(src.ntK18);
  event.vbcHS.resize(src.ntK18);
  event.xvpHS.resize(src.ntK18);
  event.yvpHS.resize(src.ntK18);
  event.zvpHS.resize(src.ntK18);
  event.uvpHS.resize(src.ntK18);
  event.vvpHS.resize(src.ntK18);

  for(int it=0; it<src.ntK18; ++it){
    event.nhK18[it] = src.nhK18[it];
    event.chisqrK18[it] = src.chisqrK18[it];
    event.p_3rd[it] = src.p_3rd[it];
    event.delta_3rd[it] = src.delta_3rd[it];
    event.pHS[it] = src.pHS[it];
    event.xoutK18[it] = src.xoutK18[it];
    event.youtK18[it] = src.youtK18[it];
    event.uoutK18[it] = src.uoutK18[it];
    event.voutK18[it] = src.voutK18[it];
    event.xtgtHS[it] = src.xtgtHS[it];
    event.ytgtHS[it] = src.ytgtHS[it];
    event.ztgtHS[it] = src.ztgtHS[it];
    event.utgtHS[it] = src.utgtHS[it];
    event.vtgtHS[it] = src.vtgtHS[it];
    event.initmomHS[it] = src.initmomHS[it];
    event.xbh2HS[it] = src.xbh2HS[it];
    event.ybh2HS[it] = src.ybh2HS[it];
    event.zbh2HS[it] = src.zbh2HS[it];
    event.ubh2HS[it] = src.ubh2HS[it];
    event.vbh2HS[it] = src.vbh2HS[it];
    event.pbh2HS[it] = src.pbh2HS[it];
    event.xgasvesselHS[it] = src.xgasvesselHS[it];
    event.ygasvesselHS[it] = src.ygasvesselHS[it];
    event.zgasvesselHS[it] = src.zgasvesselHS[it];
    event.ugasvesselHS[it] = src.ugasvesselHS[it];
    event.vgasvesselHS[it] = src.vgasvesselHS[it];
    event.pgasvesselHS[it] = src.pgasvesselHS[it];
    event.xhtofHS[it] = src.xhtofHS[it];
    event.yhtofHS[it] = src.yhtofHS[it];
    event.zhtofHS[it] = src.zhtofHS[it];
    event.uhtofHS[it] = src.uhtofHS[it];
    event.vhtofHS[it] = src.vhtofHS[it];
    event.xvpHS[it].push_back(src.xvp1HS[it]);
    event.xvpHS[it].push_back(src.xvp2HS[it]);
    event.xvpHS[it].push_back(src.xvp3HS[it]);
    event.xvpHS[it].push_back(src.xvp4HS[it]);
    event.yvpHS[it].push_back(src.yvp1HS[it]);
    event.yvpHS[it].push_back(src.yvp2HS[it]);
    event.yvpHS[it].push_back(src.yvp3HS[it]);
    event.yvpHS[it].push_back(src.yvp4HS[it]);
    event.zvpHS[it].push_back(src.zvp1HS[it]);
    event.zvpHS[it].push_back(src.zvp2HS[it]);
    event.zvpHS[it].push_back(src.zvp3HS[it]);
    event.zvpHS[it].push_back(src.zvp4HS[it]);
    event.uvpHS[it].push_back(src.uvp1HS[it]);
    event.uvpHS[it].push_back(src.uvp2HS[it]);
    event.uvpHS[it].push_back(src.uvp3HS[it]);
    event.uvpHS[it].push_back(src.uvp4HS[it]);
    event.vvpHS[it].push_back(src.vvp1HS[it]);
    event.vvpHS[it].push_back(src.vvp2HS[it]);
    event.vvpHS[it].push_back(src.vvp3HS[it]);
    event.vvpHS[it].push_back(src.vvp4HS[it]);
    for(int il=0; il<NumOfLayersBcOut; ++il){
      event.xbcHS[it].push_back(src.xbcHS[it][il]);
      event.ybcHS[it].push_back(src.ybcHS[it][il]);
      event.zbcHS[it].push_back(src.zbcHS[it][il]);
      event.ubcHS[it].push_back(src.ubcHS[it][il]);
      event.vbcHS[it].push_back(src.vbcHS[it][il]);
    }

    for(Int_t il=0; il<NumOfLayersVPHS; ++il){
      vpK18[it].push_back(TVector3(event.xvpHS[it][il], event.yvpHS[it][il], event.zvpHS[it][il]));
    }
    vpK18[it].push_back(TVector3(event.xtgtHS[it], event.ytgtHS[it], event.ztgtHS[it]));

    HF1(14, src.pHS[it]);
    TVector3 posOut(xGlobalBcOut + event.xoutK18[it],
		    yGlobalBcOut + event.youtK18[it],
		    zGlobalBcOut - zLocalBcOut);
    TVector3 momOut(event.uoutK18[it], event.voutK18[it], 1.);
    momOut *= 1./momOut.Mag();
    momOut *= src.p_3rd[it];

    initPosK18.push_back(posOut);
    initMomK18.push_back(momOut);
    for(int il=0; il<src.nhK18[it]; ++il){
      event.layerK18[it].push_back(src.layerK18[it][il]);
      event.layerK18Int[it].push_back(src.layerK18[it][il]);
      event.wireK18[it].push_back(src.wireK18[it][il]);
      event.localhitposK18[it].push_back(src.localhitposK18[it][il]);
      event.wposK18[it].push_back(src.wposK18[it][il]);
    }
  }

  std::vector<std::vector<TVector3>> vpKurama;
  vpKurama.resize( src.ntKurama );

  std::vector<Int_t> pidKurama;
  std::vector<TVector3> initPosKurama;
  std::vector<TVector3> initMomKurama;
  event.xvpKurama.resize( src.ntKurama );
  event.yvpKurama.resize( src.ntKurama );
  event.zvpKurama.resize( src.ntKurama );
  event.xhtofKurama.resize( src.ntKurama );
  event.yhtofKurama.resize( src.ntKurama );
  event.zhtofKurama.resize( src.ntKurama );
  for(Int_t it=0; it<src.ntKurama; ++it){
    event.chisqrKurama.push_back(src.chisqrKurama[it]);
    event.pKurama.push_back(src.pKurama[it]);
    event.qKurama.push_back(src.qKurama[it]);
    event.xtgtKurama.push_back(src.xtgtKurama[it]);
    event.ytgtKurama.push_back(src.ytgtKurama[it]);
    event.utgtKurama.push_back(src.utgtKurama[it]);
    event.vtgtKurama.push_back(src.vtgtKurama[it]);
    event.thetaKurama.push_back(src.thetaKurama[it]);
    event.phiKurama.push_back(src.phiKurama[it]);
    event.m2.push_back(src.m2[it]);
    event.xvpKurama[it].resize( NumOfLayersVPTPC );
    event.yvpKurama[it].resize( NumOfLayersVPTPC );
    event.zvpKurama[it].resize( NumOfLayersVPTPC );

    vpKurama[it].push_back(TVector3(event.xtgtKurama[it], event.ytgtKurama[it], tpc::ZTarget));
    for(Int_t il=0; il<NumOfLayersVPTPC; ++il){
      event.xvpKurama[it][il] = src.vpxtpc[it][il];
      event.yvpKurama[it][il] = src.vpytpc[it][il];
      event.zvpKurama[it][il] = src.vpztpc[it][il] - zHSCenter;
      vpKurama[it].push_back(TVector3(event.xvpKurama[it][il], event.yvpKurama[it][il], event.zvpKurama[it][il]));
    }
    event.xhtofKurama[it] = src.vpxhtof[it];
    event.yhtofKurama[it] = src.vpyhtof[it];
    event.zhtofKurama[it] = src.vpzhtof[it] - zHSCenter;

    Int_t pikp = -1;
    if(src.m2[it] > 0. && src.m2[it] < 0.10) pikp=0;
    else if(src.m2[it] > 0.10 && src.m2[it] < 0.4) pikp=1;
    else if(src.m2[it] > 0.4 && src.m2[it] < 1.5)  pikp=2;
    pidKurama.push_back(pikp);
  }

  event.ntKurama = src.ntKurama;
  event.layer.resize(event.ntKurama);
  event.layerInt.resize(event.ntKurama);
  event.wire.resize(event.ntKurama);
  event.localhitpos.resize(event.ntKurama);
  event.wpos.resize(event.ntKurama);
  for(Int_t i=0; i<src.ntKurama; ++i){
    event.nh.push_back(src.nh[i]);
    event.xout.push_back(src.xout[i]);
    event.yout.push_back(src.yout[i]);
    event.zout.push_back(src.zout[i]);
    event.pxout.push_back(src.pxout[i]);
    event.pyout.push_back(src.pyout[i]);
    event.pzout.push_back(src.pzout[i]);
    event.xtof.push_back(src.xtof[i]);
    event.ytof.push_back(src.ytof[i]);
    event.ztof.push_back(src.ztof[i]);
    event.pxtof.push_back(src.pxtof[i]);
    event.pytof.push_back(src.pytof[i]);
    event.pztof.push_back(src.pztof[i]);

    TVector3 posOut(xGlobalSdcOut + event.xout[i], yGlobalSdcOut + event.yout[i], event.zout[i]);
    TVector3 momOut(event.pxout[i], event.pyout[i], event.pzout[i]);
    initPosKurama.push_back(posOut);
    initMomKurama.push_back(momOut);

    for(Int_t j=0; j<event.nh[i]; ++j){
      event.layer[i].push_back(src.layer[i][j]);
      event.layerInt[i].push_back(src.layer[i][j]);
      event.wire[i].push_back(src.wire[i][j]);
      event.localhitpos[i].push_back(src.localhitpos[i][j]);
      event.wpos[i].push_back(src.wpos[i][j]);
    }
  }

  HF1( 1, event.status++ );

  if(event.clkTpc.size() != 1){
    std::cerr << "something is wrong" << std::endl;
    return true;
  }

  Double_t clock = event.clkTpc.at(0);
  TPCAnalyzer TPCAna;
  HF1( 1, event.status++ );
  TPCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, clock);

  HF1( 1, event.status++ );
  TPCAna.TrackSearchTPCHelix(vpK18, vpKurama, event.qKurama, ExclusiveTracking);

  Int_t ntTpc = TPCAna.GetNTracksTPCHelix();
  event.ntTpc = ntTpc;

  Int_t nkurama_candidates = 0;
  event.isKuramaCandidate.resize(ntTpc);
  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( it );
    if( !tp ) continue;
    Int_t iskurama = tp->GetIsKurama();
    if(iskurama==1) nkurama_candidates++;
    event.isKuramaCandidate[it] = iskurama;
  }
  event.ntKuramaCandidate = nkurama_candidates;

  HF1( 1, event.status++ );
  TPCAna.TrackSearchTPCKurama(pidKurama, initPosKurama, initMomKurama, event.layerInt, event.wire, event.localhitpos);

  HF1( 1, event.status++ );
  TPCAna.TrackSearchTPCK18(initPosK18, initMomK18, event.layerK18Int, event.wireK18, event.localhitposK18);

  std::vector<TPCLocalTrackHelix*> vptracks;
  Int_t vpntTpc = TPCAna.GetNTracksTPCVP();
  event.vpntTpc = vpntTpc;
  event.vpnhtrack.resize(vpntTpc);
  event.vptrackid.resize(vpntTpc);
  event.vpisKurama.resize(vpntTpc);
  event.vpisK18.resize(vpntTpc);
  event.vphelix_cx.resize(vpntTpc);
  event.vphelix_cy.resize(vpntTpc);
  event.vphelix_z0.resize(vpntTpc);
  event.vphelix_r.resize(vpntTpc);
  event.vphelix_dz.resize(vpntTpc);
  event.vpmom0.resize(vpntTpc);
  event.vpcharge.resize(vpntTpc);
  event.vphelix_t.resize(vpntTpc);
  event.vppos_x.resize(vpntTpc);
  event.vppos_y.resize(vpntTpc);
  event.vppos_z.resize(vpntTpc);
  event.residual_vppos_x.resize(vpntTpc);
  event.residual_vppos_y.resize(vpntTpc);
  event.residual_vppos_z.resize(vpntTpc);
  for( Int_t it=0; it<vpntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCVP( it );
    if( !tp ) continue;
    Int_t nh = tp->GetVPNHit();
    Double_t helix_cx=tp->Getcx(), helix_cy=tp->Getcy();
    Double_t helix_z0=tp->Getz0(), helix_r=tp->Getr();
    Double_t helix_dz=tp->Getdz();
    TVector3 mom0 = tp->GetMom0();
    Int_t trackid = tp->GetTrackID();
    Int_t iskurama = tp->GetIsKurama();
    Int_t isk18 = tp->GetIsK18();
    Int_t charge = tp->GetCharge();
    if(iskurama==1||isk18==1) vptracks.push_back(tp);
    event.vpnhtrack[it] = nh;
    event.vptrackid[it] = trackid;
    event.vpisKurama[it] = iskurama;
    event.vpisK18[it] = isk18;
    event.vphelix_cx[it] = helix_cx;
    event.vphelix_cy[it] = helix_cy;
    event.vphelix_z0[it] = helix_z0;
    event.vphelix_r[it] = helix_r ;
    event.vphelix_dz[it] = helix_dz;
    event.vpmom0[it] = mom0.Mag();
    event.vpcharge[it] = charge;

    event.vphelix_t[it].resize( nh );
    event.vppos_x[it].resize( nh );
    event.vppos_y[it].resize( nh );
    event.vppos_z[it].resize( nh );
    event.residual_vppos_x[it].resize( nh );
    event.residual_vppos_y[it].resize( nh );
    event.residual_vppos_z[it].resize( nh );
    for( int ih=0; ih<nh; ++ih ){
      Double_t par[5] = {helix_cx, helix_cy, helix_z0, helix_r, helix_dz};
      Double_t theta = tp->GetTheta( ih );
      const TVector3& fittmp = tp->GetPosition(par, theta);
      TVector3 calpos(-fittmp.X(), fittmp.Z(), fittmp.Y() + tpc::ZTarget);
      const TVector3& vppos = tp->GetVPPos(ih);
      event.vphelix_t[it][ih] = theta;
      event.vppos_x[it][ih] = vppos.x();
      event.vppos_y[it][ih] = vppos.y();
      event.vppos_z[it][ih] = vppos.z();
      event.residual_vppos_x[it][ih] = vppos.x() - calpos.x();
      event.residual_vppos_y[it][ih] = vppos.y() - calpos.y();
      event.residual_vppos_z[it][ih] = vppos.z() - calpos.z();
    }
  }

  event.isgoodTPCK18.resize( src.ntK18 );
  event.tpcidTPCK18.resize( src.ntK18 );
  event.niterationTPCK18.resize( src.ntK18 );
  event.chisqrTPCK18.resize( src.ntK18 );
  event.qTPCK18.resize( src.ntK18 );
  event.pTPCK18.resize( src.ntK18 );
  event.xtgtTPCK18.resize( src.ntK18 );
  event.ytgtTPCK18.resize( src.ntK18 );
  event.utgtTPCK18.resize( src.ntK18 );
  event.vtgtTPCK18.resize( src.ntK18 );
  event.thetaTPCK18.resize( src.ntK18 );
  event.lhtofTPCK18.resize( src.ntK18 );
  event.xhtofTPCK18.resize( src.ntK18 );
  event.yhtofTPCK18.resize( src.ntK18 );
  event.lvpTPCK18.resize( src.ntK18 );
  event.xvpTPCK18.resize( src.ntK18 );
  event.yvpTPCK18.resize( src.ntK18 );
  for(Int_t itk18=0; itk18<TPCAna.GetNTracksTPCK18(); ++itk18){
    TPCRKTrack* tr_km = TPCAna.GetTPCK18Track(itk18);
    Int_t niteration = tr_km -> Niteration();
    Double_t chisqr = tr_km -> GetChiSquare();
    const TVector3& tgtmom = tr_km -> TargetMomentum();
    const TVector3& tgtpos = tr_km -> TargetPosition();
    Double_t utgt = tgtmom.x()/tgtmom.z();
    Double_t vtgt = tgtmom.y()/tgtmom.z();
    Double_t pt = tgtmom.Mag()/TMath::Sqrt(1. + utgt*utgt + vtgt*vtgt);
    Double_t q = tr_km -> Polarity();
    ThreeVector posCorr(tgtpos.x()+xb_off, tgtpos.y()+yb_off, tgtpos.z());
    ThreeVector momCorr(pt*(utgt+ub_off), pt*(vtgt+vb_off), pt);
    Double_t uCorr = momCorr.x()/momCorr.z();
    Double_t vCorr = momCorr.y()/momCorr.z();
    Double_t cost = 1./TMath::Sqrt(1.+uCorr*uCorr+vCorr*vCorr);
    Double_t theta = TMath::ACos(cost)*TMath::RadToDeg();

    Int_t idtpc = tr_km -> GetTPCTrackID();
    Int_t idk18 = tr_km -> GetTrackID();
    event.isgoodTPCK18[idk18] = 1;
    event.tpcidTPCK18[idk18] = idtpc;
    event.niterationTPCK18[idk18] = niteration;
    event.chisqrTPCK18[idk18] = chisqr;
    event.qTPCK18[idk18] = q;
    event.pTPCK18[idk18] = tgtmom.Mag()*pK18_factor+pK18_offset;
    event.xtgtTPCK18[idk18] = posCorr.x();
    event.ytgtTPCK18[idk18] = posCorr.y();
    event.utgtTPCK18[idk18] = uCorr;
    event.vtgtTPCK18[idk18] = vCorr;
    event.thetaTPCK18[idk18] = theta;

    Double_t path, x, y;
    tr_km -> GetTrajectoryLocalPosition(20, path, x, y);
    event.lhtofTPCK18[idk18] = path;
    event.xhtofTPCK18[idk18] = x;
    event.yhtofTPCK18[idk18] = y;
    event.lvpTPCK18[idk18].resize( NumOfLayersVPHS );
    event.xvpTPCK18[idk18].resize( NumOfLayersVPHS );
    event.yvpTPCK18[idk18].resize( NumOfLayersVPHS );
    for(Int_t l = 0; l<NumOfLayersVPHS; ++l){
      tr_km -> GetTrajectoryLocalPosition(208 + l, path, x, y);
      event.lvpTPCK18[idk18][l] = path;
      event.xvpTPCK18[idk18][l] = x;
      event.yvpTPCK18[idk18][l] = y;
    }// for(l)
    for(Int_t i=0; i<NumOfLayersBcOut; ++i){
      Double_t resolution; Double_t residual;
      if(tr_km -> GetTrajectoryResidual(i+PlOffsBcOut+1, resolution, residual)){
	HF1(TPCK18RKHid+100+i, residual/resolution);
	HF1(TPCK18RKHid+200+i, residual);
      }
    }

    Int_t nhTpc = tr_km -> GetTPCTrack() -> GetNHit();
    for(Int_t i=0; i<nhTpc; ++i){
      TVector3 resolution; TVector3 residual;
      if(tr_km -> GetTrajectoryResidualTPC(i, resolution, residual)){
	if(resolution.x() > 0.9e+10 && resolution.y() > 0.9e+10 && resolution.z() > 0.9e+10) continue; // exclude bad hits
	Int_t layer = tr_km -> GetTPCTrack() -> GetHitInOrder(i) -> GetLayer();
	const TVector3& resi_vect = tr_km -> GetTPCTrack() -> GetHitInOrder(i) -> GetResidualVect();
	HF1(TPCK18RKHid+layer+300, residual.x()/resolution.x());
	HF1(TPCK18RKHid+layer+400, residual.y()/resolution.y());
	HF1(TPCK18RKHid+layer+500, residual.x());
	HF1(TPCK18RKHid+layer+600, residual.y());
	HF1(TPCK18RKHid+layer+700, resi_vect.x());
	HF1(TPCK18RKHid+layer+800, resi_vect.y());
	HF1(TPCK18RKHid+layer+900, resi_vect.z());
      }
    }

    for(Int_t it=0; it<vptracks.size(); ++it){
      TPCLocalTrackHelix *tr_vptrack = vptracks[it];
      if(tr_vptrack->GetTrackID()==idk18){
	TPCLocalTrackHelix *tr_tpc = tr_km -> GetTPCTrack();
	Int_t nh = tr_tpc -> GetNHit();
	for( int ih=0; ih<nh; ++ih ){
	  TPCLTrackHit *hit = tr_tpc -> GetHitInOrder(ih);
	  Int_t layer = hit -> GetLayer();
	  const TVector3& hitpos = hit -> GetLocalHitPos();
	  TVector3 resi = tr_vptrack -> CalcResidual(hitpos);
	  Double_t xz_resi = TMath::Sqrt(resi.x()*resi.x()+resi.z()*resi.z());
	  Double_t x_resi = resi.x();
	  Double_t y_resi = resi.y();
	  Double_t z_resi = resi.z();
	  HF1(TPCK18VPHid+layer, xz_resi);
	  HF1(TPCK18VPHid+100+layer, x_resi);
	  HF1(TPCK18VPHid+200+layer, y_resi);
	  HF1(TPCK18VPHid+300+layer, z_resi);
	}
      }
    }
  }

  event.isgoodTPCKurama.resize( src.ntKurama );
  event.tpcidTPCKurama.resize( src.ntKurama );
  event.niterationTPCKurama.resize( src.ntKurama );
  event.kflagTPCKurama.resize( src.ntKurama );
  event.pflagTPCKurama.resize( src.ntKurama );
  event.chisqrTPCKurama.resize( src.ntKurama );
  event.pTPCKurama.resize( src.ntKurama );
  event.qTPCKurama.resize( src.ntKurama );
  event.m2TPCKurama.resize( src.ntKurama );
  event.m2OrgTPCKurama.resize( src.ntKurama );
  event.xtgtTPCKurama.resize( src.ntKurama );
  event.ytgtTPCKurama.resize( src.ntKurama );
  event.utgtTPCKurama.resize( src.ntKurama );
  event.vtgtTPCKurama.resize( src.ntKurama );
  event.thetaTPCKurama.resize( src.ntKurama );
  event.pathTPCKurama.resize( src.ntKurama );
  event.lhtofTPCKurama.resize( src.ntKurama );
  event.xhtofTPCKurama.resize( src.ntKurama );
  event.yhtofTPCKurama.resize( src.ntKurama );
  event.phtofTPCKurama.resize( src.ntKurama );
  event.lgasvesselTPCKurama.resize( src.ntKurama );
  event.xgasvesselTPCKurama.resize( src.ntKurama );
  event.ygasvesselTPCKurama.resize( src.ntKurama );
  event.pgasvesselTPCKurama.resize( src.ntKurama );
  event.lvpTPCKurama.resize( src.ntKurama );
  event.xvpTPCKurama.resize( src.ntKurama );
  event.yvpTPCKurama.resize( src.ntKurama );
  for(Int_t ittpckurama=0; ittpckurama<TPCAna.GetNTracksTPCKurama(); ++ittpckurama){
    TPCRKTrack* tr_kp = TPCAna.GetTPCKuramaTrack(ittpckurama);
    Int_t idtpc = tr_kp -> GetTPCTrackID();
    Int_t idkurama = tr_kp -> GetTrackID();

    Int_t niteration = tr_kp -> Niteration();
    Double_t chisqr = tr_kp -> GetChiSquare();
    const TVector3& tgtpos = tr_kp -> TargetPosition();
    const TVector3& tgtmom = tr_kp -> TargetMomentum();
    Double_t q = tr_kp -> Polarity();
    Double_t utgt = tgtmom.x()/tgtmom.z();
    Double_t vtgt = tgtmom.y()/tgtmom.z();
    Double_t pOrg = tgtmom.Mag();
    Double_t pCorr = pKuramaCorrection(utgt, vtgt, pOrg);
    Double_t cost = 1./TMath::Sqrt(1.+utgt*utgt+vtgt*vtgt);
    Double_t theta = TMath::ACos(cost)*TMath::RadToDeg();
    Double_t pathtof = tr_kp -> PathLengthToTOF();

    Double_t cstof = src.stof[idkurama];
    if(cstof > 0.){
      event.m2TPCKurama[idkurama] = Kinematics::MassSquare(pCorr, pathtof, cstof);
      event.m2OrgTPCKurama[idkurama] = Kinematics::MassSquare(pOrg, pathtof, cstof);
      if(event.m2TPCKurama[idkurama] > 0.4 && event.m2TPCKurama[idkurama] < 1.5){
	pCorr = pKuramaCorrection_proton(utgt, vtgt, pOrg);
      }
    }
    else{
      event.m2TPCKurama[idkurama] = TMath::QuietNaN();
      event.m2OrgTPCKurama[idkurama] = TMath::QuietNaN();
    }

    event.isgoodTPCKurama[idkurama] = 1;
    event.tpcidTPCKurama[idkurama] = idtpc;
    event.niterationTPCKurama[idkurama] = niteration;
    event.chisqrTPCKurama[idkurama] = chisqr;
    event.pTPCKurama[idkurama] = pCorr;
    event.qTPCKurama[idkurama] = q;
    event.xtgtTPCKurama[idkurama] = tgtpos.x();
    event.ytgtTPCKurama[idkurama] = tgtpos.y();
    event.utgtTPCKurama[idkurama] = utgt;
    event.vtgtTPCKurama[idkurama] = vtgt;
    event.thetaTPCKurama[idkurama] = theta;
    event.pathTPCKurama[idkurama] = pathtof;
    if(q>0 && event.m2TPCKurama[idkurama] > 0.10 && event.m2TPCKurama[idkurama] < 0.40) event.kflagTPCKurama[idkurama] = 1;
    if(q>0 && event.m2TPCKurama[idkurama] > 0.5 && event.m2TPCKurama[idkurama] < 1.4) event.pflagTPCKurama[idkurama] = 1;

    Double_t path, x, y;
    TVector3 mom;
    tr_kp -> GetTrajectoryLocalPosition(IdVPHTOF, path, x, y);
    event.lhtofTPCKurama[idkurama] = path;
    event.xhtofTPCKurama[idkurama] = x;
    event.yhtofTPCKurama[idkurama] = y;
    tr_kp -> GetTrajectoryMomentum(IdVPHTOF, mom);
    event.phtofTPCKurama[idkurama] = mom.Mag();

    tr_kp -> GetTrajectoryLocalPosition(IdTPCGasVessel_D, path, x, y);
    event.lgasvesselTPCKurama[idkurama] = path;
    event.xgasvesselTPCKurama[idkurama] = x;
    event.ygasvesselTPCKurama[idkurama] = y;
    tr_kp -> GetTrajectoryMomentum(IdTPCGasVessel_D, mom);
    event.pgasvesselTPCKurama[idkurama] = mom.Mag();

    event.lvpTPCKurama[idkurama].resize( NumOfLayersVPTPC );
    event.xvpTPCKurama[idkurama].resize( NumOfLayersVPTPC );
    event.yvpTPCKurama[idkurama].resize( NumOfLayersVPTPC );
    for(Int_t l = 0; l<NumOfLayersVPTPC; ++l){
      tr_kp -> GetTrajectoryLocalPosition(21 + l, path, x, y);
      event.lvpTPCKurama[idkurama][l] = path;
      event.xvpTPCKurama[idkurama][l] = x;
      event.yvpTPCKurama[idkurama][l] = y;
    }// for(l)
  }

  int nkk = src.ntKurama*src.ntK18;
  event.isgoodTPC.resize(nkk);
  event.insideTPC.resize(nkk);
  event.vtxTPC.resize(nkk);
  event.vtyTPC.resize(nkk);
  event.vtzTPC.resize(nkk);
  event.pxKmTPC.resize(nkk);
  event.pyKmTPC.resize(nkk);
  event.pzKmTPC.resize(nkk);
  event.pxScatTPC.resize(nkk);
  event.pyScatTPC.resize(nkk);
  event.pzScatTPC.resize(nkk);

  event.closeDistTPC.resize(nkk);
  event.MissMassTPC.resize(nkk);
  event.MissMassCorrTPC.resize(nkk);
  event.MissMassCorrDETPC.resize(nkk);
  event.MissMassNuclTPC.resize(nkk);
  event.MissMassNuclCorrTPC.resize(nkk);
  event.MissMassNuclCorrDETPC.resize(nkk);
  event.pOrgTPC.resize(nkk);
  event.pCorrTPC.resize(nkk);
  event.pCorrDETPC.resize(nkk);
  event.pCalcTPC.resize(nkk);
  event.thetaCMTPC.resize(nkk);
  event.costCMTPC.resize(nkk);
  event.pCalcDETPC.resize(nkk);
  event.thetaCMDETPC.resize(nkk);
  event.costCMDETPC.resize(nkk);
  event.xistarpCalcDETPC.resize(nkk);
  event.xistarthetaCMDETPC.resize(nkk);
  event.xistarcostCMDETPC.resize(nkk);
  event.kpscatpCalcTPC.resize(nkk);
  event.kpscatthetaCMTPC.resize(nkk);
  event.kpscatcostCMTPC.resize(nkk);
  event.kpscatpCalcDETPC.resize(nkk);
  event.kpscatthetaCMDETPC.resize(nkk);
  event.kpscatcostCMDETPC.resize(nkk);
  event.thetaTPC.resize(nkk);
  event.ubTPC.resize(nkk);
  event.vbTPC.resize(nkk);
  event.usTPC.resize(nkk);
  event.vsTPC.resize(nkk);
  for(Int_t itkurama=0; itkurama<TPCAna.GetNTracksTPCKurama(); ++itkurama){
    TPCRKTrack* trScat = TPCAna.GetTPCKuramaTrack(itkurama);
    Int_t idScat = trScat -> GetTrackID();
    const TVector3& tgtmomScat = trScat -> TargetMomentum();
    const TVector3& tgtposScat = trScat -> TargetPosition();
    Double_t us = tgtmomScat.x()/tgtmomScat.z(), vs = tgtmomScat.y()/tgtmomScat.z();
    Double_t pOrg = tgtmomScat.Mag();
    Double_t pCorr = event.pTPCKurama[idScat];
    ThreeVector pScat = tgtmomScat;
    ThreeVector pScatCorr(pCorr/pOrg*tgtmomScat.x(), pCorr/pOrg*tgtmomScat.y(), pCorr/pOrg*tgtmomScat.z());
    ThreeVector xScat(tgtposScat.x(), tgtposScat.y(), tgtposScat.z());

    for(Int_t itk18=0; itk18<src.ntK18; ++itk18){
      Int_t idKm = itk18;
      Int_t id = idScat*src.ntK18 + idKm;
      Double_t ub = src.utgtHS[idKm], vb = src.vtgtHS[idKm];
      Double_t ptKm = src.pHS[idKm]/TMath::Sqrt(1. + ub*ub + vb*vb);
      ThreeVector pKm(ptKm*ub, ptKm*vb, ptKm);
      ThreeVector pKmCorr(ptKm*(ub+ub_off), ptKm*(vb+vb_off), ptKm);
      ThreeVector xKm(src.xtgtHS[idKm]+xb_off, src.ytgtHS[idKm]+xb_off, tgtposScat.z());

      //PID
      Double_t ScatMass = qnan;
      if(event.m2TPCKurama[idScat] > 0. &&
	 event.m2TPCKurama[idScat] < 0.10) ScatMass = PionMass;
      else if(event.m2TPCKurama[idScat] > 0.10 &&
	      event.m2TPCKurama[idScat] < 0.4) ScatMass = KaonMass;
      else if(event.m2TPCKurama[idScat] > 0.4 &&
	      event.m2TPCKurama[idScat] < 1.5) ScatMass = ProtonMass;

      //Reaction vertex
      TVector3 KKVertex; TVector3 KmMomVtx; TVector3 KpMomVtx;
      Double_t closeDist; Double_t KmPathInTgt; Double_t KpPathInTgt;
      //Bool_t IsKKVertexInTarget = TPCAna.GetProductionVertex(xKm, pKmCorr, xScat, pScatCorr, KKVertex, closeDist, KmPathInTgt, KpPathInTgt, KmMomVtx, KpMomVtx);
      TPCAna.GetProductionVertex(xKm, pKmCorr, xScat, pScatCorr, KKVertex, closeDist, KmPathInTgt, KpPathInTgt, KmMomVtx, KpMomVtx);

      Int_t inside = 0;
      if(TMath::Abs(KKVertex.x()) < 30.
         && TMath::Abs(KKVertex.y()) < 30.
         && TMath::Abs(KKVertex.z()) < 70.
         && closeDist < 30.) inside = 1;

      //Eloss correction
      Int_t target_material = 2; //Carbon
      if(event.runnum >= 5641 && event.runnum <= 5666) target_material = 1; //CH2
      ThreeVector pKmCorrDE = Kinematics::HypTPCCorrElossIn(target_material, pKmCorr, KmPathInTgt, KaonMass);
      ThreeVector pScatCorrDE = Kinematics::HypTPCCorrElossOut(target_material, pScatCorr, KpPathInTgt, ScatMass);

      //Missing-Mass reconstruction
      LorentzVector LvKm(pKm, std::sqrt(KaonMass*KaonMass+pKm.Mag2()));
      LorentzVector LvKmCorr(pKmCorr, std::sqrt(KaonMass*KaonMass+pKmCorr.Mag2()));
      LorentzVector LvKmCorrDE(pKmCorrDE, sqrt(KaonMass*KaonMass+pKmCorrDE.Mag2()));
      LorentzVector LvScat(pScat, std::sqrt(ScatMass*ScatMass+pScat.Mag2()));
      LorentzVector LvScatCorr(pScatCorr, std::sqrt(ScatMass*ScatMass+pScatCorr.Mag2()));
      LorentzVector LvScatCorrDE(pScatCorrDE, std::sqrt(ScatMass*ScatMass+pScatCorrDE.Mag2()));
      // proton
      LorentzVector LvC(0., 0., 0., ProtonMass);
      LorentzVector LvRc = LvKm + LvC - LvScat;
      LorentzVector LvRcCorr = LvKmCorr + LvC - LvScatCorr;
      LorentzVector LvRcCorrDE = LvKmCorrDE + LvC - LvScatCorrDE;

      Double_t MissMass = LvRc.Mag();
      Double_t MissMassCorr = LvRcCorr.Mag();
      Double_t MissMassCorrDE = LvRcCorrDE.Mag();//-LvC.Mag();

      // Carbon12 nucleus
      LorentzVector LvCNucl(0., 0., 0., Carbon12Mass);
      LorentzVector LvRcNucl = LvKm + LvCNucl - LvScat;
      LorentzVector LvRcNuclCorr = LvKmCorr + LvCNucl - LvScatCorr;
      LorentzVector LvRcNuclCorrDE = LvKmCorrDE + LvCNucl - LvScatCorrDE;
      Double_t MissMassNucl = LvRcNucl.Mag();
      Double_t MissMassNuclCorr = LvRcNuclCorr.Mag();
      Double_t MissMassNuclCorrDE = LvRcNuclCorrDE.Mag();//-LvC.Mag();

      Double_t cost = pKmCorr*pScatCorr/(pKmCorr.Mag()*pScatCorr.Mag());
      Double_t theta = TMath::ACos(cost)*TMath::RadToDeg();
      { //Xi-, dE in the target is not considered
	//CM

	//Primary frame
	LorentzVector PrimaryLv = LvKmCorr+LvC;
	Double_t TotalEnergyCM = PrimaryLv.Mag();
	ThreeVector beta(1/PrimaryLv.E()*PrimaryLv.Vect());

	Double_t TotalMomCM
	  = 0.5*std::sqrt((TotalEnergyCM*TotalEnergyCM
			   -(KaonMass+XiMass)*(KaonMass+XiMass))
			  *(TotalEnergyCM*TotalEnergyCM
			    -(KaonMass-XiMass)*(KaonMass-XiMass)))/TotalEnergyCM;

	Double_t costLab = cost;
	Double_t cottLab = costLab/std::sqrt(1.-costLab*costLab);
	Double_t bt = beta.Mag(), gamma = 1./std::sqrt(1.-bt*bt);
	Double_t gbep = gamma*bt*std::sqrt(TotalMomCM*TotalMomCM+KaonMass*KaonMass)/TotalMomCM;
	Double_t a  = gamma*gamma+cottLab*cottLab;
	Double_t bp = gamma*gbep;
	Double_t c  = gbep*gbep-cottLab*cottLab;
	Double_t dd = bp*bp-a*c;

	if(dd<0.){
	  std::cerr << "dd<0." << std::endl;
	  dd = 0.;
	}

	Double_t costCM = (std::sqrt(dd)-bp)/a;
	if(costCM>1. || costCM<-1.){
	  std::cerr << "costCM>1. || costCM<-1." << std::endl;
	  costCM=-1.;
	}
	Double_t sintCM  = std::sqrt(1.-costCM*costCM);
	Double_t KaonMom = TotalMomCM*sintCM/std::sqrt(1.-costLab*costLab);
	event.thetaCMTPC[id] = TMath::ACos(costCM)*TMath::RadToDeg();
	event.costCMTPC[id] = costCM;
	event.pCalcTPC[id] = KaonMom;
      }
      { //Xi-
	//CM

	//Primary frame
	LorentzVector PrimaryLv = LvKmCorrDE+LvC;
	Double_t TotalEnergyCM = PrimaryLv.Mag();
	ThreeVector beta(1/PrimaryLv.E()*PrimaryLv.Vect());

	Double_t TotalMomCM
	  = 0.5*std::sqrt((TotalEnergyCM*TotalEnergyCM
			   -(KaonMass+XiMass)*(KaonMass+XiMass))
			  *(TotalEnergyCM*TotalEnergyCM
			    -(KaonMass-XiMass)*(KaonMass-XiMass)))/TotalEnergyCM;

	Double_t costLab = cost;
	Double_t cottLab = costLab/std::sqrt(1.-costLab*costLab);
	Double_t bt = beta.Mag(), gamma = 1./std::sqrt(1.-bt*bt);
	Double_t gbep = gamma*bt*std::sqrt(TotalMomCM*TotalMomCM+KaonMass*KaonMass)/TotalMomCM;
	Double_t a  = gamma*gamma+cottLab*cottLab;
	Double_t bp = gamma*gbep;
	Double_t c  = gbep*gbep-cottLab*cottLab;
	Double_t dd = bp*bp-a*c;

	if(dd<0.){
	  std::cerr << "dd<0." << std::endl;
	  dd = 0.;
	}

	Double_t costCM = (std::sqrt(dd)-bp)/a;
	if(costCM>1. || costCM<-1.){
	  std::cerr << "costCM>1. || costCM<-1." << std::endl;
	  costCM=-1.;
	}
	Double_t sintCM  = std::sqrt(1.-costCM*costCM);
	Double_t KaonMom = TotalMomCM*sintCM/std::sqrt(1.-costLab*costLab);
	event.thetaCMDETPC[id] = TMath::ACos(costCM)*TMath::RadToDeg();
	event.costCMDETPC[id] = costCM;
	event.pCalcDETPC[id] = KaonMom;
      }
      {
	//Primary frame
	LorentzVector PrimaryLv = LvKmCorrDE+LvC;
	Double_t TotalEnergyCM = PrimaryLv.Mag();
	ThreeVector beta(1/PrimaryLv.E()*PrimaryLv.Vect());

	//CM Xi(1530)
	Double_t TotalMomCM
	  = 0.5*std::sqrt((TotalEnergyCM*TotalEnergyCM
			   -(KaonMass+XiStarMass)*(KaonMass+XiStarMass))
			  *(TotalEnergyCM*TotalEnergyCM
			    -(KaonMass-XiStarMass)*(KaonMass-XiStarMass)))/TotalEnergyCM;
	Double_t costLab = cost;
	Double_t cottLab = costLab/std::sqrt(1.-costLab*costLab);
	Double_t bt = beta.Mag(), gamma = 1./std::sqrt(1.-bt*bt);
	Double_t gbep = gamma*bt*std::sqrt(TotalMomCM*TotalMomCM+KaonMass*KaonMass)/TotalMomCM;
	Double_t a  = gamma*gamma+cottLab*cottLab;
	Double_t bp = gamma*gbep;
	Double_t c  = gbep*gbep-cottLab*cottLab;
	Double_t dd = bp*bp-a*c;

	if(dd<0.){
	  std::cerr << "dd<0." << std::endl;
	  dd = 0.;
	}

	Double_t costCM = (std::sqrt(dd)-bp)/a;
	if(costCM>1. || costCM<-1.){
	  std::cerr << "costCM>1. || costCM<-1." << std::endl;
	  costCM=-1.;
	}
	Double_t sintCM  = std::sqrt(1.-costCM*costCM);
	Double_t KaonMom = TotalMomCM*sintCM/std::sqrt(1.-costLab*costLab);
	event.xistarthetaCMDETPC[id] = TMath::ACos(costCM)*TMath::RadToDeg();
	event.xistarcostCMDETPC[id] = costCM;
	event.xistarpCalcDETPC[id] = KaonMom;
      }

      event.isgoodTPC[id] = 1;
      event.insideTPC[id] = inside;
      event.vtxTPC[id] = KKVertex.x();
      event.vtyTPC[id] = KKVertex.y();
      event.vtzTPC[id] = KKVertex.z();
      event.pxKmTPC[id] = KmMomVtx.x();
      event.pyKmTPC[id] = KmMomVtx.y();
      event.pzKmTPC[id] = KmMomVtx.z();
      event.pxScatTPC[id] = KpMomVtx.x();
      event.pyScatTPC[id] = KpMomVtx.y();
      event.pzScatTPC[id] = KpMomVtx.z();

      event.closeDistTPC[id] = closeDist;
      event.MissMassTPC[id] = MissMass;
      event.MissMassCorrTPC[id] = MissMassCorr;
      event.MissMassCorrDETPC[id] = MissMassCorrDE;
      event.MissMassNuclTPC[id] = MissMassNucl;
      event.MissMassNuclCorrTPC[id] = MissMassNuclCorr;
      event.MissMassNuclCorrDETPC[id] = MissMassNuclCorrDE;

      event.thetaTPC[id] = theta;
      event.pOrgTPC[id] = pOrg;
      event.pCorrTPC[id] = pCorr;
      event.pCorrDETPC[id] = pScatCorrDE.Mag();

      event.ubTPC[id] = ub;
      event.vbTPC[id] = vb;
      event.usTPC[id] = us;
      event.vsTPC[id] = vs;

      if(event.chisqrK18[idKm] < MaxChisqrBcOut && event.chisqrKurama[idScat] < MaxChisqrKurama){
	HF1(5001, event.pTPCK18[idKm]);
	HF1(5002, event.pTPCKurama[idScat]);
	HF1(5003, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]));
	HF2(5004, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]), event.pTPCKurama[idScat]);
	HF2(5010, KKVertex.z() + tpc::ZTarget, KKVertex.x());
	HF1(5011, closeDist);
	HF1(5012, KKVertex.x());
	HF1(5013, KKVertex.y());
	HF1(5014, KKVertex.z());

	if(inside==1){
	  {
	    TPCLocalTrackHelix *tr_kuramavp;
	    for(Int_t it=0; it<vptracks.size(); ++it){
	      if(vptracks[it]->GetTrackID()==idScat){
		tr_kuramavp = vptracks[it];
		break;
	      }
	    }
	    TPCLocalTrackHelix *tr_tpcKurama = trScat -> GetTPCTrack();
	    Int_t nhTpc = tr_tpcKurama -> GetNHit();
	    for(Int_t i=0; i<nhTpc; ++i){

	      TVector3 resolution; TVector3 resi_TPCKurama;
	      if(trScat -> GetTrajectoryResidualTPC(i, resolution, resi_TPCKurama)){
		if(resolution.x() > 0.9e+10 && resolution.y() > 0.9e+10 && resolution.z() > 0.9e+10) continue; // exclude bad hits

		Int_t layer = tr_tpcKurama -> GetHitInOrder(i) -> GetLayer();
		const TVector3& hitpos = tr_tpcKurama -> GetHitInOrder(i) -> GetLocalHitPos();

		TVector3 resi_Kurama = tr_kuramavp -> CalcResidual(hitpos);
		const TVector3& resi_TPC = tr_tpcKurama -> GetHitInOrder(i) -> GetResidualVect();
		if(event.qTPCKurama[idScat]>0.){
		  HF2(TPCKuramaRKHid+layer+6100, resi_Kurama.x(), resi_TPCKurama.x());
		  HF2(TPCKuramaRKHid+layer+6200, resi_Kurama.x(), resi_TPC.x());
		  HF2(TPCKuramaRKHid+layer+6300, resi_TPC.x(), resi_TPCKurama.x());
		  HF2(TPCKuramaRKHid+layer+6400, resi_Kurama.y(), resi_TPCKurama.y());
		  HF2(TPCKuramaRKHid+layer+6500, resi_Kurama.y(), resi_TPC.y());
		  HF2(TPCKuramaRKHid+layer+6600, resi_TPC.y(), resi_TPCKurama.y());
		}
		else{
		  HF2(TPCKuramaRKHid+layer+7100, resi_Kurama.x(), resi_TPCKurama.x());
		  HF2(TPCKuramaRKHid+layer+7200, resi_Kurama.x(), resi_TPC.x());
		  HF2(TPCKuramaRKHid+layer+7300, resi_TPC.x(), resi_TPCKurama.x());
		  HF2(TPCKuramaRKHid+layer+7400, resi_Kurama.y(), resi_TPCKurama.y());
		  HF2(TPCKuramaRKHid+layer+7500, resi_Kurama.y(), resi_TPC.y());
		  HF2(TPCKuramaRKHid+layer+7600, resi_TPC.y(), resi_TPCKurama.y());
		}
	      }
	    }
	  }

	  for(Int_t it=0; it<vptracks.size(); ++it){
	    TPCLocalTrackHelix *tr_vptrack = vptracks[it];
	    if(tr_vptrack->GetTrackID()==idScat){
	      TPCLocalTrackHelix *tr_tpc = trScat -> GetTPCTrack();
	      Int_t nh = tr_tpc -> GetNHit();
	      for( int ih=0; ih<nh; ++ih ){
		TPCLTrackHit *hit = tr_tpc -> GetHitInOrder(ih);
		Int_t layer = hit -> GetLayer();
		const TVector3& hitpos = hit -> GetLocalHitPos();
		TVector3 resi = tr_vptrack -> CalcResidual(hitpos);
		double xz_resi = TMath::Sqrt(resi.x()*resi.x()+resi.z()*resi.z());
		double x_resi = resi.x();
		double y_resi = resi.y();
		double z_resi = resi.z();
		if(event.qTPCKurama[idScat]>0.){
		  HF1(TPCKuramaVPHid+layer, xz_resi);
		  HF1(TPCKuramaVPHid+100+layer, x_resi);
		  HF1(TPCKuramaVPHid+200+layer, y_resi);
		  HF1(TPCKuramaVPHid+300+layer, z_resi);
		}
		else{
		  HF1(TPCKuramaVPHid+1000+layer, xz_resi);
		  HF1(TPCKuramaVPHid+1100+layer, x_resi);
		  HF1(TPCKuramaVPHid+1200+layer, y_resi);
		  HF1(TPCKuramaVPHid+1300+layer, z_resi);
		}

	      }
	    }
	  }

	  for(Int_t i=0; i<NumOfLayersSdcIn; ++i){
	    Double_t resolution; Double_t residual;
	    if(trScat -> GetTrajectoryResidual(i+PlOffsSdcIn+1, resolution, residual)){
	      if(event.qTPCKurama[idScat]>0.){
		HF1(TPCKuramaRKHid+100+i, residual/resolution);
		HF1(TPCKuramaRKHid+200+i, residual);
	      }
	      else{
		HF1(TPCKuramaRKHid+1100+i, residual/resolution);
		HF1(TPCKuramaRKHid+1200+i, residual);
	      }
	    }
	  }
	  for(Int_t i=0; i<NumOfLayersSdcOut; ++i){
	    Double_t resolution; Double_t residual;
	    if(trScat -> GetTrajectoryResidual(i+PlOffsSdcOut+1, resolution, residual)){
	      if(event.qTPCKurama[idScat]>0.){
		HF1(TPCKuramaRKHid+300+i, residual/resolution);
		HF1(TPCKuramaRKHid+400+i, residual);
	      }
	      else{
		HF1(TPCKuramaRKHid+1300+i, residual/resolution);
		HF1(TPCKuramaRKHid+1400+i, residual);
	      }
	    }
	  }

	  Int_t nhTpc = trScat -> GetTPCTrack() -> GetNHit();
	  for(Int_t i=0; i<nhTpc; ++i){
	    TVector3 resolution; TVector3 residual;
	    if(trScat -> GetTrajectoryResidualTPC(i, resolution, residual)){
	      if(resolution.x() > 0.9e+10 && resolution.y() > 0.9e+10 && resolution.z() > 0.9e+10) continue; // exclude bad hits
	      Int_t layer = trScat -> GetTPCTrack() -> GetHitInOrder(i) -> GetLayer();
	      const TVector3& resi_vect = trScat -> GetTPCTrack() -> GetHitInOrder(i) -> GetResidualVect();
	      if(event.qTPCKurama[idScat]>0.){
		HF1(TPCKuramaRKHid+layer+500, residual.x()/resolution.x());
		HF1(TPCKuramaRKHid+layer+600, residual.y()/resolution.y());
		HF1(TPCKuramaRKHid+layer+700, residual.x());
		HF1(TPCKuramaRKHid+layer+800, residual.y());

		HF1(TPCKuramaRKHid+layer+4100, resi_vect.x());
		HF1(TPCKuramaRKHid+layer+4200, resi_vect.y());
		HF1(TPCKuramaRKHid+layer+4300, resi_vect.z());
	      }
	      else{
		HF1(TPCKuramaRKHid+layer+1500, residual.x()/resolution.x());
		HF1(TPCKuramaRKHid+layer+1600, residual.y()/resolution.y());
		HF1(TPCKuramaRKHid+layer+1700, residual.x());
		HF1(TPCKuramaRKHid+layer+1800, residual.y());

		HF1(TPCKuramaRKHid+layer+5100, resi_vect.x());
		HF1(TPCKuramaRKHid+layer+5200, resi_vect.y());
		HF1(TPCKuramaRKHid+layer+5300, resi_vect.z());
	      }
	    }
	  }

	  HF1(6001, event.pTPCK18[idKm]);
	  HF1(6002, event.pTPCKurama[idScat]);
	  HF1(6003, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]));
	  HF2(6004, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]), event.pTPCKurama[idScat]);
	  HF2(6010, KKVertex.z() + tpc::ZTarget, KKVertex.x());
	  HF1(6011, closeDist);
	  HF1(6012, KKVertex.x());
	  HF1(6013, KKVertex.y());
	  HF1(6014, KKVertex.z());

	  Double_t cstof = src.stof[idScat];
	  Int_t tofseg = src.tofsegKurama[id];
	  if(event.pTPCKurama[idScat] < 1.4 && event.m2TPCKurama[idScat] > 0. && event.m2TPCKurama[idScat] < 0.12){
	    Double_t calctof = Kinematics::CalcTimeOfFlight(event.pTPCKurama[idScat], event.pathTPCKurama[idScat], pdg::PionMass());
	    HF1(10000+tofseg*100+1, cstof-calctof);
	  }
	  if(event.m2TPCKurama[idScat] > 0.15 && event.m2TPCKurama[idScat] < 0.35){
	    Double_t calctof = Kinematics::CalcTimeOfFlight(event.pTPCKurama[idScat], event.pathTPCKurama[idScat], pdg::KaonMass());
	    HF1(10000+tofseg*100+2, cstof-calctof);
	  }
	  if(event.m2TPCKurama[idScat] > 0.50 && event.m2TPCKurama[idScat] < 1.5){
	    Double_t calctof = Kinematics::CalcTimeOfFlight(event.pTPCKurama[idScat], event.pathTPCKurama[idScat], pdg::ProtonMass());
	    HF1(10000+tofseg*100+3, cstof-calctof);
	  }

	  //K+
	  if(event.qTPCKurama[idScat] > 0 && event.pTPCKurama[idScat] < 1.4 &&
	     event.m2TPCKurama[idScat] > 0.12 && event.m2TPCKurama[idScat] < 0.3){
	    HF1(7001, event.pTPCK18[idKm]);
	    HF1(7002, event.pTPCKurama[idScat]);
	    HF1(7003, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]));
	    HF2(7004, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]), event.pTPCKurama[idScat]);
	    HF1(7010, closeDist);
	    HF1(7011, KKVertex.x());
	    HF1(7012, KKVertex.y());
	    HF1(7013, KKVertex.z());
	    HF1(7014, MissMassCorr);
	    HF1(7015, MissMassCorrDE);
	    HF2(7016, us, pScat.Mag() - event.pCalcTPC[id]);
	    HF2(7017, us, pScatCorr.Mag() - event.pCalcTPC[id]);
	    HF2(7018, us, pScatCorrDE.Mag() - event.pCalcDETPC[id]);
	    HF2(7019, us, MissMass);
	    HF2(7020, us, MissMassCorr);
	    HF2(7021, us, MissMassCorrDE);
	    HF2(7022, vs, pScat.Mag() - event.pCalcTPC[id]);
	    HF2(7023, vs, pScatCorr.Mag() - event.pCalcTPC[id]);
	    HF2(7024, vs, pScatCorrDE.Mag() - event.pCalcDETPC[id]);
	    HF2(7025, vs, MissMass);
	    HF2(7026, vs, MissMassCorr);
	    HF2(7027, vs, MissMassCorrDE);
	    HF2(7028, event.pCalcTPC[id], pScat.Mag() - event.pCalcTPC[id]);
	    HF2(7029, event.pCalcTPC[id], pScatCorr.Mag() - event.pCalcTPC[id]);
	    HF2(7030, event.pCalcDETPC[id], pScatCorrDE.Mag() - event.pCalcDETPC[id]);
	    HF2(7031, event.pCalcTPC[id], MissMass);
	    HF2(7032, event.pCalcTPC[id], MissMassCorr);
	    HF2(7033, event.pCalcDETPC[id], MissMassCorrDE);
	    HF2(7036, event.xistarpCalcDETPC[id], pScatCorrDE.Mag() - event.xistarpCalcDETPC[id]);
	    HF2(7039, event.xistarpCalcDETPC[id], MissMassCorrDE);

	    HFProf(7116, us, pScat.Mag() - event.pCalcTPC[id]);
	    HFProf(7117, us, pScatCorr.Mag() - event.pCalcTPC[id]);
	    HFProf(7118, us, pScatCorrDE.Mag() - event.pCalcDETPC[id]);
	    HFProf(7119, us, MissMass);
	    HFProf(7120, us, MissMassCorr);
	    HFProf(7121, us, MissMassCorrDE);
	    HFProf(7122, vs, pScat.Mag() - event.pCalcTPC[id]);
	    HFProf(7123, vs, pScatCorr.Mag() - event.pCalcTPC[id]);
	    HFProf(7124, vs, pScatCorrDE.Mag() - event.pCalcDETPC[id]);
	    HFProf(7125, vs, MissMass);
	    HFProf(7126, vs, MissMassCorr);
	    HFProf(7127, vs, MissMassCorrDE);
	    HFProf(7128, event.pCalcTPC[id], pScat.Mag() - event.pCalcTPC[id]);
	    HFProf(7129, event.pCalcTPC[id], pScatCorr.Mag() - event.pCalcTPC[id]);
	    HFProf(7130, event.pCalcDETPC[id], pScatCorrDE.Mag() - event.pCalcDETPC[id]);
	    HFProf(7131, event.pCalcTPC[id], MissMass);
	    HFProf(7132, event.pCalcTPC[id], MissMassCorr);
	    HFProf(7133, event.pCalcDETPC[id], MissMassCorrDE);
	    HFProf(7136, event.xistarpCalcDETPC[id], pScatCorrDE.Mag() - event.xistarpCalcDETPC[id]);
	    HFProf(7139, event.xistarpCalcDETPC[id], MissMassCorrDE);

	    if(event.pTPCKurama[idScat] > 1.1){
	      HF1(8001, event.pTPCK18[idKm]);
	      HF1(8002, event.pTPCKurama[idScat]);
	      HF1(8003, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]));
	      HF2(8004, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]), event.pTPCKurama[idScat]);
	      HF1(8010, closeDist);
	      HF1(8011, KKVertex.x());
	      HF1(8012, KKVertex.y());
	      HF1(8013, KKVertex.z());
	      HF1(8014, MissMassCorr);
	      HF1(8015, MissMassCorrDE);
	      HF2(8016, us, pScat.Mag() - event.pCalcTPC[id]);
	      HF2(8017, us, pScatCorr.Mag() - event.pCalcTPC[id]);
	      HF2(8018, us, pScatCorrDE.Mag() - event.pCalcDETPC[id]);
	      HF2(8019, us, MissMass);
	      HF2(8020, us, MissMassCorr);
	      HF2(8021, us, MissMassCorrDE);
	      HF2(8022, vs, pScat.Mag() - event.pCalcTPC[id]);
	      HF2(8023, vs, pScatCorr.Mag() - event.pCalcTPC[id]);
	      HF2(8024, vs, pScatCorrDE.Mag() - event.pCalcDETPC[id]);
	      HF2(8025, vs, MissMass);
	      HF2(8026, vs, MissMassCorr);
	      HF2(8027, vs, MissMassCorrDE);
	      HF2(8028, event.pCalcTPC[id], pScat.Mag() - event.pCalcTPC[id]);
	      HF2(8029, event.pCalcTPC[id], pScatCorr.Mag() - event.pCalcTPC[id]);
	      HF2(8030, event.pCalcDETPC[id], pScatCorrDE.Mag() - event.pCalcDETPC[id]);
	      HF2(8031, event.pCalcTPC[id], MissMass);
	      HF2(8032, event.pCalcTPC[id], MissMassCorr);
	      HF2(8033, event.pCalcDETPC[id], MissMassCorrDE);
	      HF2(8036, event.xistarpCalcDETPC[id], pScatCorrDE.Mag() - event.xistarpCalcDETPC[id]);

	      HFProf(8116, us, pScat.Mag() - event.pCalcTPC[id]);
	      HFProf(8117, us, pScatCorr.Mag() - event.pCalcTPC[id]);
	      HFProf(8118, us, pScatCorrDE.Mag() - event.pCalcDETPC[id]);
	      HFProf(8119, us, MissMass);
	      HFProf(8120, us, MissMassCorr);
	      HFProf(8121, us, MissMassCorrDE);
	      HFProf(8122, vs, pScat.Mag() - event.pCalcTPC[id]);
	      HFProf(8123, vs, pScatCorr.Mag() - event.pCalcTPC[id]);
	      HFProf(8124, vs, pScatCorrDE.Mag() - event.pCalcDETPC[id]);
	      HFProf(8125, vs, MissMass);
	      HFProf(8126, vs, MissMassCorr);
	      HFProf(8127, vs, MissMassCorrDE);
	      HFProf(8128, event.pCalcTPC[id], pScat.Mag() - event.pCalcTPC[id]);
	      HFProf(8129, event.pCalcTPC[id], pScatCorr.Mag() - event.pCalcTPC[id]);
	      HFProf(8130, event.pCalcDETPC[id], pScatCorrDE.Mag() - event.pCalcDETPC[id]);
	      HFProf(8131, event.pCalcTPC[id], MissMass);
	      HFProf(8132, event.pCalcTPC[id], MissMassCorr);
	      HFProf(8133, event.pCalcDETPC[id], MissMassCorrDE);
	    }
	    else{
	      HF2(8036, event.xistarpCalcDETPC[id], pScatCorrDE.Mag() - event.xistarpCalcDETPC[id]);
	      HF2(8039, event.xistarpCalcDETPC[id], MissMassCorrDE);

	      HFProf(8136, event.xistarpCalcDETPC[id], pScatCorrDE.Mag() - event.xistarpCalcDETPC[id]);
	      HFProf(8139, event.xistarpCalcDETPC[id], MissMassCorrDE);
	    }
	  }
	  //Proton
	  if(event.qTPCKurama[idScat] > 0 &&
	     event.m2TPCKurama[idScat] > 0.5 && event.m2TPCKurama[idScat] < 1.5){

	    { //Kp Scattering, dE in the target is not considered
	      //Primary frame
	      LorentzVector PrimaryLv = LvKmCorr+LvC;
	      Double_t TotalEnergyCM = PrimaryLv.Mag();
	      ThreeVector beta(1/PrimaryLv.E()*PrimaryLv.Vect());

	      Double_t TotalMomCM
		= 0.5*std::sqrt((TotalEnergyCM*TotalEnergyCM
				 -(KaonMass+ProtonMass)*(KaonMass+ProtonMass))
				*(TotalEnergyCM*TotalEnergyCM
				  -(KaonMass-ProtonMass)*(KaonMass-ProtonMass)))/TotalEnergyCM;
	      Double_t costLab = cost;
	      Double_t cottLab = costLab/std::sqrt(1.-costLab*costLab);
	      Double_t bt = beta.Mag(), gamma = 1./std::sqrt(1.-bt*bt);
	      Double_t gbep = gamma*bt*std::sqrt(TotalMomCM*TotalMomCM+ProtonMass*ProtonMass)/TotalMomCM;
	      Double_t a  = gamma*gamma+cottLab*cottLab;
	      Double_t bp = gamma*gbep;
	      Double_t c  = gbep*gbep-cottLab*cottLab;
	      Double_t dd = bp*bp-a*c;
	      if(dd<0.){
		std::cerr << "dd<0." << std::endl;
		dd = 0.;
	      }

	      Double_t costCM = (std::sqrt(dd)-bp)/a;
	      if(costCM>1. || costCM<-1.){
		std::cerr << "costCM>1. || costCM<-1." << std::endl;
		costCM=-1.;
	      }
	      Double_t sintCM  = std::sqrt(1.-costCM*costCM);
	      Double_t ProtonMom = TotalMomCM*sintCM/std::sqrt(1.-costLab*costLab);
	      event.kpscatthetaCMTPC[id] = TMath::ACos(costCM)*TMath::RadToDeg();
	      event.kpscatcostCMTPC[id] = costCM;
	      event.kpscatpCalcTPC[id] = ProtonMom;
	    }
	    { //Kp Scattering
	      //Primary frame
	      LorentzVector PrimaryLv = LvKmCorrDE+LvC;
	      Double_t TotalEnergyCM = PrimaryLv.Mag();
	      ThreeVector beta(1/PrimaryLv.E()*PrimaryLv.Vect());

	      Double_t TotalMomCM
		= 0.5*std::sqrt((TotalEnergyCM*TotalEnergyCM
				 -(KaonMass+ProtonMass)*(KaonMass+ProtonMass))
				*(TotalEnergyCM*TotalEnergyCM
				  -(KaonMass-ProtonMass)*(KaonMass-ProtonMass)))/TotalEnergyCM;
	      Double_t costLab = cost;
	      Double_t cottLab = costLab/std::sqrt(1.-costLab*costLab);
	      Double_t bt = beta.Mag(), gamma = 1./std::sqrt(1.-bt*bt);
	      Double_t gbep = gamma*bt*std::sqrt(TotalMomCM*TotalMomCM+ProtonMass*ProtonMass)/TotalMomCM;
	      Double_t a  = gamma*gamma+cottLab*cottLab;
	      Double_t bp = gamma*gbep;
	      Double_t c  = gbep*gbep-cottLab*cottLab;
	      Double_t dd = bp*bp-a*c;
	      if(dd<0.){
		std::cerr << "dd<0." << std::endl;
		dd = 0.;
	      }

	      Double_t costCM = (std::sqrt(dd)-bp)/a;
	      if(costCM>1. || costCM<-1.){
		std::cerr << "costCM>1. || costCM<-1." << std::endl;
		costCM=-1.;
	      }
	      Double_t sintCM  = std::sqrt(1.-costCM*costCM);
	      Double_t ProtonMom = TotalMomCM*sintCM/std::sqrt(1.-costLab*costLab);
	      event.kpscatthetaCMDETPC[id] = TMath::ACos(costCM)*TMath::RadToDeg();
	      event.kpscatcostCMDETPC[id] = costCM;
	      event.kpscatpCalcDETPC[id] = ProtonMom;
	    }

	    HF1(8201, event.pTPCK18[idKm]);
	    HF1(8202, event.pTPCKurama[idScat]);
	    HF1(8203, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]));
	    HF2(8204, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]), event.pTPCKurama[idScat]);
	    HF1(8210, closeDist);
	    HF1(8211, KKVertex.x());
	    HF1(8212, KKVertex.y());
	    HF1(8213, KKVertex.z());
	    HF1(8214, MissMassCorr);
	    HF1(8215, MissMassCorrDE);
	    HF2(8216, us, pScat.Mag() - event.kpscatpCalcTPC[id]);
	    HF2(8217, us, pScatCorr.Mag() - event.kpscatpCalcTPC[id]);
	    HF2(8218, us, pScatCorrDE.Mag() - event.kpscatpCalcDETPC[id]);
	    HF2(8219, us, MissMass);
	    HF2(8220, us, MissMassCorr);
	    HF2(8221, us, MissMassCorrDE);
	    HF2(8222, us, pScat.Mag() - event.kpscatpCalcTPC[id]);
	    HF2(8223, us, pScatCorr.Mag() - event.kpscatpCalcTPC[id]);
	    HF2(8224, us, pScatCorrDE.Mag() - event.kpscatpCalcDETPC[id]);
	    HF2(8225, vs, MissMass);
	    HF2(8226, vs, MissMassCorr);
	    HF2(8227, vs, MissMassCorrDE);
	    HF2(8228, event.kpscatpCalcTPC[id], pScat.Mag() - event.kpscatpCalcTPC[id]);
	    HF2(8229, event.kpscatpCalcTPC[id], pScatCorr.Mag() - event.kpscatpCalcTPC[id]);
	    HF2(8230, event.kpscatpCalcDETPC[id], pScatCorrDE.Mag() - event.kpscatpCalcDETPC[id]);
	    HF2(8231, event.kpscatpCalcTPC[id], MissMass);
	    HF2(8232, event.kpscatpCalcTPC[id], MissMassCorr);
	    HF2(8233, event.kpscatpCalcDETPC[id], MissMassCorrDE);

	    HFProf(8316, us, pScat.Mag() - event.kpscatpCalcTPC[id]);
	    HFProf(8317, us, pScatCorr.Mag() - event.kpscatpCalcTPC[id]);
	    HFProf(8318, us, pScatCorrDE.Mag() - event.kpscatpCalcDETPC[id]);
	    HFProf(8319, us, MissMass);
	    HFProf(8320, us, MissMassCorr);
	    HFProf(8321, us, MissMassCorrDE);
	    HFProf(8322, us, pScat.Mag() - event.kpscatpCalcTPC[id]);
	    HFProf(8323, us, pScatCorr.Mag() - event.kpscatpCalcTPC[id]);
	    HFProf(8324, us, pScatCorrDE.Mag() - event.kpscatpCalcDETPC[id]);
	    HFProf(8325, vs, MissMass);
	    HFProf(8326, vs, MissMassCorr);
	    HFProf(8327, vs, MissMassCorrDE);
	    HFProf(8328, event.kpscatpCalcTPC[id], pScat.Mag() - event.kpscatpCalcTPC[id]);
	    HFProf(8329, event.kpscatpCalcTPC[id], pScatCorr.Mag() - event.kpscatpCalcTPC[id]);
	    HFProf(8330, event.kpscatpCalcDETPC[id], pScatCorrDE.Mag() - event.kpscatpCalcDETPC[id]);
	    HFProf(8331, event.kpscatpCalcTPC[id], MissMass);
	    HFProf(8332, event.kpscatpCalcTPC[id], MissMassCorr);
	    HFProf(8333, event.kpscatpCalcDETPC[id], MissMassCorrDE);
	  }
	}
      }
    }
  }
  HF1( 1, event.status++ );

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
  HF1(1, event.status++);
  event.nhTpc = nh_Tpc;

#if RawCluster
  Int_t nclTpc = 0;
  Int_t remain_nclTpc = 0;
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
      TPCHit* meanHit = cl->GetMeanHit();
      Int_t houghflag = meanHit->GetHoughFlag();
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
      event.cluster_houghflag.push_back(houghflag);
      ++nclTpc;

      if(houghflag!=100&&houghflag!=200) ++remain_nclTpc; //Clusters without track
    }
  }
  event.nclTpc = nclTpc;
  event.remain_nclTpc = remain_nclTpc;
  HF1( 1, event.status++ );
#endif

  HF1( 10, ntTpc );
  //if( event.ntTpc == 0 ) return true;

  HF1( 1, event.status++ );
  event.nhtrack.resize( ntTpc );
  event.nhtrackEff.resize( ntTpc );
  event.flag.resize( ntTpc );
  event.trackid.resize( ntTpc );
  event.isBeam.resize( ntTpc );
  event.isXi.resize( ntTpc );
  event.isKurama.resize( ntTpc );
  event.isK18.resize( ntTpc );
  event.isAccidental.resize( ntTpc );
  event.isMultiloop.resize( ntTpc );
  event.fittime.resize( ntTpc );
  event.searchtime.resize( ntTpc );
  event.niteration.resize( ntTpc );
  event.chisqr.resize( ntTpc );
  event.pval.resize( ntTpc );
  event.distTgt.resize( ntTpc );

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
#endif
  event.dz_factor.resize( ntTpc );
  event.charge.resize( ntTpc );
  event.path.resize( ntTpc );

  event.pid.resize( ntTpc );
  event.hitlayer.resize( ntTpc );
  event.hitpos_x.resize( ntTpc );
  event.hitpos_y.resize( ntTpc );
  event.hitpos_z.resize( ntTpc );
  event.calpos_x.resize( ntTpc );
  event.calpos_y.resize( ntTpc );
  event.calpos_z.resize( ntTpc );
  event.mom_x.resize( ntTpc );
  event.mom_y.resize( ntTpc );
  event.mom_z.resize( ntTpc );
  event.residual.resize( ntTpc );
  event.residual_t.resize( ntTpc );
  event.residual_x.resize( ntTpc );
  event.residual_y.resize( ntTpc );
  event.residual_z.resize( ntTpc );
  event.resolution_x.resize( ntTpc );
  event.resolution_y.resize( ntTpc );
  event.resolution_z.resize( ntTpc );
  event.pull.resize( ntTpc );
  event.helix_t.resize( ntTpc );
  event.pathhit.resize(ntTpc);
  event.alpha.resize(ntTpc);
  event.houghflag.resize(ntTpc);
  event.track_cluster_de.resize(ntTpc);
  event.track_cluster_size.resize(ntTpc);
  event.track_cluster_mrow.resize(ntTpc);
  event.track_cluster_de_center.resize(ntTpc);
  event.track_cluster_x_center.resize(ntTpc);
  event.track_cluster_y_center.resize(ntTpc);
  event.track_cluster_z_center.resize(ntTpc);
  event.track_cluster_row_center.resize(ntTpc);

  event.exresidual_t.resize( ntTpc );
  event.exresidual_x.resize( ntTpc );
  event.exresidual_y.resize( ntTpc );
  event.exresidual_z.resize( ntTpc );
  event.intrinsic_residual_t.resize( ntTpc );
  event.intrinsic_residual_x.resize( ntTpc );
  event.intrinsic_residual_y.resize( ntTpc );
  event.intrinsic_residual_z.resize( ntTpc );

  event.chargeIndistinguishable.resize( ntTpc );
  event.chisqr_inverted.resize( ntTpc );
  event.pval_inverted.resize( ntTpc );
  event.helix_cx_inverted.resize( ntTpc );
  event.helix_cy_inverted.resize( ntTpc );
  event.helix_z0_inverted.resize( ntTpc );
  event.helix_r_inverted.resize( ntTpc );
  event.helix_dz_inverted.resize( ntTpc );
  event.mom0_x_inverted.resize( ntTpc );
  event.mom0_y_inverted.resize( ntTpc );
  event.mom0_z_inverted.resize( ntTpc );
  event.mom0_inverted.resize( ntTpc );
  event.pid_inverted.resize( ntTpc );

  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( it );
    if( !tp ) continue;
    Int_t nh = tp->GetNHit();
    Int_t nhEff = tp->GetNHitsEffective();
    Double_t chisqr = tp->GetChiSquare();
    Double_t pval = 1-ROOT::Math::chisquared_cdf(chisqr*(2*nhEff-5), 2*nhEff-5);
    Double_t helix_cx = tp->Getcx(), helix_cy = tp->Getcy();
    Double_t helix_z0 = tp->Getz0(), helix_r = tp->Getr();
    Double_t helix_dz = tp->Getdz();
    TVector3 mom0 = tp->GetMom0();
    Int_t flag = tp->GetFitFlag();
    Int_t trackid = tp->GetTrackID();
    Int_t isbeam = tp->GetIsBeam();
    Int_t isxi = tp->GetIsXi();
    Int_t iskurama = tp->GetIsKurama();
    Int_t isk18 = tp->GetIsK18();
    Int_t isaccidental = tp->GetIsAccidental();
    Int_t ismultiloop = tp->GetIsMultiloop();
    Int_t charge = tp->GetCharge();
    Int_t pid = tp->GetPid();
    Int_t iteration = tp->GetNIteration();
    Double_t fittime = tp->GetFitTime();
    Double_t searchtime = tp->GetSearchTime();
    Double_t pathlen = tp->GetPath();
    Double_t distTgt = tp->GetClosestDist();

    HF1(11, nh);
    HF1(12, chisqr);

#if TruncatedMean
    event.dEdx_0[it]=tp->GetdEdx(1.0);
    event.dEdx_10[it]=tp->GetdEdx(0.9);
    event.dEdx_20[it]=tp->GetdEdx(0.8);
    event.dEdx_30[it]=tp->GetdEdx(0.7);
    event.dEdx_40[it]=tp->GetdEdx(0.6);
    event.dEdx_50[it]=tp->GetdEdx(0.5);
    event.dEdx_60[it]=tp->GetdEdx(0.4);
#endif

    event.nhtrack[it] = nh;
    event.nhtrackEff[it] = nhEff;
    event.flag[it] = flag;
    event.trackid[it] = trackid;
    event.isBeam[it] = isbeam;
    event.isXi[it] = isxi;
    event.isKurama[it] = iskurama;
    event.isK18[it] = isk18;
    event.isAccidental[it] = isaccidental;
    event.isMultiloop[it] = ismultiloop;
    event.fittime[it] = fittime;
    event.charge[it] = charge;
    event.path[it] = pathlen;
    event.distTgt[it] = distTgt;
    event.chisqr[it] = chisqr;
    event.pval[it] = pval;
    event.niteration[it] = iteration;
    event.searchtime[it] = searchtime;
    event.helix_cx[it] = helix_cx;
    event.helix_cy[it] = helix_cy;
    event.helix_z0[it] = helix_z0;
    event.helix_r[it] = helix_r ;
    event.helix_dz[it] = helix_dz;
    event.mom0_x[it] = mom0.x();
    event.mom0_y[it] = mom0.y();
    event.mom0_z[it] = mom0.z();
    event.mom0[it] = mom0.Mag();
    event.pid[it] = pid;
    event.dE[it] = tp->GetTrackdE();
    event.dEdx[it] = tp->GetdEdx(truncatedMean);
    event.dz_factor[it] = sqrt(1.+(pow(helix_dz,2)));

    HF2(20, event.mom0[it]*event.charge[it], event.dEdx[it]);
    if(event.charge[it]<0){
      HF2(21, -event.mom0[it]*event.charge[it], event.dEdx[it]);
    } else {
      HF2(22, event.mom0[it]*event.charge[it], event.dEdx[it]);
    }

    HF1(15, event.mom0[it]);
    if(src.ntK18==1) HF1(16, event.mom0[it]-src.pHS[0]);

    event.hitlayer[it].resize( nh );
    event.hitpos_x[it].resize( nh );
    event.hitpos_y[it].resize( nh );
    event.hitpos_z[it].resize( nh );
    event.calpos_x[it].resize( nh );
    event.calpos_y[it].resize( nh );
    event.calpos_z[it].resize( nh );
    event.mom_x[it].resize( nh );
    event.mom_y[it].resize( nh );
    event.mom_z[it].resize( nh );
    event.residual[it].resize( nh );
    event.residual_t[it].resize( nh );
    event.residual_x[it].resize( nh );
    event.residual_y[it].resize( nh );
    event.residual_z[it].resize( nh );
    event.resolution_x[it].resize( nh );
    event.resolution_y[it].resize( nh );
    event.resolution_z[it].resize( nh );
    event.pull[it].resize( nh );
    event.helix_t[it].resize( nh );
    event.pathhit[it].resize(nh);
    event.alpha[it].resize(nh);
    event.houghflag[it].resize(nh);
    event.track_cluster_de[it].resize(nh);
    event.track_cluster_size[it].resize(nh);
    event.track_cluster_mrow[it].resize(nh);
    event.track_cluster_de_center[it].resize(nh);
    event.track_cluster_x_center[it].resize(nh);
    event.track_cluster_y_center[it].resize(nh);
    event.track_cluster_z_center[it].resize(nh);
    event.track_cluster_row_center[it].resize(nh);
    event.exresidual_t[it].resize( nh );
    event.exresidual_x[it].resize( nh );
    event.exresidual_y[it].resize( nh );
    event.exresidual_z[it].resize( nh );
    event.intrinsic_residual_t[it].resize( nh );
    event.intrinsic_residual_x[it].resize( nh );
    event.intrinsic_residual_y[it].resize( nh );
    event.intrinsic_residual_z[it].resize( nh );

    for( int ih=0; ih<nh; ++ih ){
      TPCLTrackHit *hit = tp->GetHitInOrder( ih );
      if( !hit ) continue;

      HF1( 2, hit->GetHoughDist());
      HF1( 3, hit->GetHoughDistY());
      Int_t layer = hit->GetLayer();
      Int_t houghflag = hit->GetHoughFlag();
      Double_t residual = hit->GetResidual();
      const TVector3& resi_vect = hit->GetResidualVect();
      const TVector3& res_vect = hit->GetResolutionVect();
      const TVector3& hitpos = hit->GetLocalHitPos();
      const TVector3& calpos = hit->GetLocalCalPosHelix();
      const TVector3& mom = hit->GetMomentumHelix(charge);

      HF1(13, layer);

      Double_t clde = hit->GetDe();
      Double_t mrow = hit->GetMRow();
      event.track_cluster_de[it][ih] = clde;
      event.track_cluster_mrow[it][ih] = mrow;
      event.alpha[it][ih] = tp->GetAlpha(ih);

      TPCHit *clhit = hit->GetHit();
      TPCCluster *cl = clhit->GetParentCluster();
      Int_t clsize = cl->GetClusterSize();
      //Double_t mrow = cl->MeanRow(); // same
      TPCHit* centerHit = cl->GetCenterHit();
      const TVector3& centerPos = centerHit->GetPosition();
      Double_t centerDe = centerHit->GetCDe();
      Int_t centerRow = centerHit->GetRow();

      HF1(TPCClHid, clsize);
      HF1(TPCClHid+(layer+1)*1000, clsize);
      HF1(TPCClHid+1, clde);
      HF1(TPCClHid+(layer+1)*1000+1, clde);
      const TPCHitContainer& hc = cl -> GetHitContainer();
      for(const auto& hits : hc){
	if(!hits || !hits->IsGood()) continue;
	const TVector3& pos = hits->GetPosition();
	Double_t de = hits->GetCDe();
	Double_t transDist = TranseverseDistance(hitpos.x(), hitpos.z(), pos.x(), pos.z());
	Double_t ratio = de/clde;
	HF2(TPCClHid+2, transDist, ratio);
	HF2(TPCClHid+(layer+1)*1000+2, transDist, ratio);
      }
      event.track_cluster_size[it][ih] = clsize;
      event.track_cluster_de_center[it][ih] = centerDe;
      event.track_cluster_x_center[it][ih] = centerPos.X();
      event.track_cluster_y_center[it][ih] = centerPos.Y();
      event.track_cluster_z_center[it][ih] = centerPos.Z();
      event.track_cluster_row_center[it][ih] = centerRow;
      event.pathhit[it][ih] = hit->GetPathHelix();

      event.hitlayer[it][ih] = (double)layer;
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
      event.residual_x[it][ih] = resi_vect.x();
      event.residual_y[it][ih] = resi_vect.y();
      event.residual_z[it][ih] = resi_vect.z();
      event.resolution_x[it][ih] = res_vect.x();
      event.resolution_y[it][ih] = res_vect.y();
      event.resolution_z[it][ih] = res_vect.z();
      event.houghflag[it][ih] = houghflag;
      event.helix_t[it][ih] = hit->GetTheta();

      Double_t resi_theta = TMath::ATan2(resi_vect.z(), -resi_vect.x());
      TVector3 dir_hit(TMath::Cos(hit->GetTheta()), TMath::Sin(hit->GetTheta()), 0);
      TVector3 dir_resi(TMath::Cos(resi_theta), TMath::Sin(resi_theta), 0);
      Double_t sign = 1.;
      if(dir_hit*dir_resi<0) sign = -1.;
      Double_t resi_t = sign*TMath::Hypot(resi_vect.x(), resi_vect.z());
      event.residual_t[it][ih] = resi_t;
      Double_t res_t = TMath::Hypot(res_vect.x(), res_vect.z());
      event.pull[it][ih] = hypot(resi_t/res_t, resi_vect.y()/resi_vect.y());

      HF1(TPCInclusiveHid+layer,resi_t);
      HF1(TPCInclusiveHid+100+layer,resi_vect.x());
      HF1(TPCInclusiveHid+200+layer,resi_vect.y());
      HF1(TPCInclusiveHid+300+layer,resi_vect.z());
      HF1(TPCInclusiveHid+1000+layer,resi_t/res_t);
      HF1(TPCInclusiveHid+1100+layer,resi_vect.x()/res_vect.x());
      HF1(TPCInclusiveHid+1200+layer,resi_vect.y()/res_vect.y());
      HF1(TPCInclusiveHid+1300+layer,resi_vect.z()/res_vect.z());

      HF1(TPCInclusiveHid+32,resi_t);
      HF1(TPCInclusiveHid+100+32,resi_vect.x());
      HF1(TPCInclusiveHid+200+32,resi_vect.y());
      HF1(TPCInclusiveHid+300+32,resi_vect.z());
      HF1(TPCInclusiveHid+1000+32,resi_t/res_t);
      HF1(TPCInclusiveHid+1100+32,resi_vect.x()/res_vect.x());
      HF1(TPCInclusiveHid+1200+32,resi_vect.y()/res_vect.y());
      HF1(TPCInclusiveHid+1300+32,resi_vect.z()/res_vect.z());

      HF2(TPCInclusiveHid+10000,layer,resi_t);
      HF2(TPCInclusiveHid+100+10000,layer,resi_vect.x());
      HF2(TPCInclusiveHid+200+10000,layer,resi_vect.y());
      HF2(TPCInclusiveHid+300+10000,layer,resi_vect.z());
      HF2(TPCInclusiveHid+1000+10000,layer,resi_t/res_t);
      HF2(TPCInclusiveHid+1100+10000,layer,resi_vect.x()/res_vect.x());
      HF2(TPCInclusiveHid+1200+10000,layer,resi_vect.y()/res_vect.y());
      HF2(TPCInclusiveHid+1300+10000,layer,resi_vect.z()/res_vect.z());
      if(ExclusiveTracking){
	const TVector3& exresi_vect = hit->GetResidualVectExclusive();
	Double_t exresi_x = exresi_vect.x();
	Double_t exresi_y = exresi_vect.y();
	Double_t exresi_z = exresi_vect.z();
	Double_t exresi_t = sign*hypot(exresi_x,exresi_z);
	event.exresidual_t[it][ih] = exresi_t;
	event.exresidual_x[it][ih] = exresi_x;
	event.exresidual_y[it][ih] = exresi_y;
	event.exresidual_z[it][ih] = exresi_z;

	HF1(TPCExclusiveHid+layer,exresi_t);
	HF1(TPCExclusiveHid+100+layer,exresi_x);
	HF1(TPCExclusiveHid+200+layer,exresi_y);
	HF1(TPCExclusiveHid+300+layer,exresi_z);
	HF1(TPCExclusiveHid+1000+layer,exresi_t/res_t);
	HF1(TPCExclusiveHid+1100+layer,exresi_x/res_vect.x());
	HF1(TPCExclusiveHid+1200+layer,exresi_y/res_vect.y());
	HF1(TPCExclusiveHid+1300+layer,exresi_z/res_vect.z());

	HF1(TPCExclusiveHid+32,exresi_t);
	HF1(TPCExclusiveHid+100+32,exresi_x);
	HF1(TPCExclusiveHid+200+32,exresi_y);
	HF1(TPCExclusiveHid+300+32,exresi_z);
	HF1(TPCExclusiveHid+1000+32,exresi_t/res_t);
	HF1(TPCExclusiveHid+1100+32,exresi_x/res_vect.x());
	HF1(TPCExclusiveHid+1200+32,exresi_y/res_vect.y());
	HF1(TPCExclusiveHid+1300+32,exresi_z/res_vect.z());

	HF2(TPCExclusiveHid+10000,layer,exresi_t);
	HF2(TPCExclusiveHid+100+10000,layer,exresi_x);
	HF2(TPCExclusiveHid+200+10000,layer,exresi_y);
	HF2(TPCExclusiveHid+300+10000,layer,exresi_z);
	HF2(TPCExclusiveHid+1000+10000,layer,exresi_t/res_t);
	HF2(TPCExclusiveHid+1100+10000,layer,exresi_x/res_vect.x());
	HF2(TPCExclusiveHid+1200+10000,layer,exresi_y/res_vect.y());
	HF2(TPCExclusiveHid+1300+10000,layer,exresi_z/res_vect.z());

	Double_t intrinsic_resi_t = sqrt(abs(resi_t*exresi_t));
	Double_t intrinsic_resi_x = sqrt(abs(resi_vect.x()*exresi_x));
	Double_t intrinsic_resi_y = sqrt(abs(resi_vect.y()*exresi_y));
	Double_t intrinsic_resi_z = sqrt(abs(resi_vect.z()*exresi_z));
	if(resi_t<0) intrinsic_resi_t*=-1;
	if(resi_vect.x()<0) intrinsic_resi_x*=-1;
	if(resi_vect.y()<0) intrinsic_resi_y*=-1;
	if(resi_vect.z()<0) intrinsic_resi_z*=-1;
	event.intrinsic_residual_t[it][ih] = intrinsic_resi_t;
	event.intrinsic_residual_x[it][ih] = intrinsic_resi_x;
	event.intrinsic_residual_y[it][ih] = intrinsic_resi_y;
	event.intrinsic_residual_z[it][ih] = intrinsic_resi_z;

	HF1(TPCIntrinsicHid+layer,intrinsic_resi_t);
	HF1(TPCIntrinsicHid+100+layer,intrinsic_resi_x);
	HF1(TPCIntrinsicHid+200+layer,intrinsic_resi_y);
	HF1(TPCIntrinsicHid+300+layer,intrinsic_resi_z);
	HF1(TPCIntrinsicHid+1000+layer,intrinsic_resi_t/res_t);
	HF1(TPCIntrinsicHid+1100+layer,intrinsic_resi_x/res_vect.x());
	HF1(TPCIntrinsicHid+1200+layer,intrinsic_resi_y/res_vect.y());
	HF1(TPCIntrinsicHid+1300+layer,intrinsic_resi_z/res_vect.z());

	HF1(TPCIntrinsicHid+32,intrinsic_resi_t);
	HF1(TPCIntrinsicHid+100+32,intrinsic_resi_x);
	HF1(TPCIntrinsicHid+200+32,intrinsic_resi_y);
	HF1(TPCIntrinsicHid+300+32,intrinsic_resi_z);
	HF1(TPCIntrinsicHid+1000+32,intrinsic_resi_t/res_t);
	HF1(TPCIntrinsicHid+1100+32,intrinsic_resi_x/res_vect.x());
	HF1(TPCIntrinsicHid+1200+32,intrinsic_resi_y/res_vect.y());
	HF1(TPCIntrinsicHid+1300+32,intrinsic_resi_z/res_vect.z());

	HF2(TPCIntrinsicHid+10000,layer,intrinsic_resi_t);
	HF2(TPCIntrinsicHid+100+10000,layer,intrinsic_resi_x);
	HF2(TPCIntrinsicHid+200+10000,layer,intrinsic_resi_y);
	HF2(TPCIntrinsicHid+300+10000,layer,intrinsic_resi_z);
	HF2(TPCIntrinsicHid+1000+10000,layer,intrinsic_resi_t/res_t);
	HF2(TPCIntrinsicHid+1100+10000,layer,intrinsic_resi_x/res_vect.x());
	HF2(TPCIntrinsicHid+1200+10000,layer,intrinsic_resi_y/res_vect.y());
	HF2(TPCIntrinsicHid+1300+10000,layer,intrinsic_resi_z/res_vect.z());
      }
    }

    //Inverted charge tracks
    TPCLocalTrackHelix *tp_inverted = TPCAna.GetTrackTPCHelixChargeInverted( it );
    if( !tp_inverted ) event.chargeIndistinguishable[it] = 0;
    else{
      Double_t chisqr = tp_inverted->GetChiSquare();
      Double_t pval = 1-ROOT::Math::chisquared_cdf(chisqr*(2*nhEff-5), 2*nhEff-5);
      Double_t helix_cx = tp_inverted->Getcx(), helix_cy = tp_inverted->Getcy();
      Double_t helix_z0 = tp_inverted->Getz0(), helix_r = tp_inverted->Getr();
      Double_t helix_dz = tp_inverted->Getdz();
      TVector3 mom0 = tp_inverted->GetMom0();
      Int_t charge = tp_inverted->GetCharge();
      Int_t pid = tp_inverted->GetPid();

      event.chargeIndistinguishable[it] = 1;
      event.chisqr_inverted[it] = chisqr;
      event.pval_inverted[it] = pval;
      event.helix_cx_inverted[it] = helix_cx;
      event.helix_cy_inverted[it] = helix_cy;
      event.helix_z0_inverted[it] = helix_z0;
      event.helix_r_inverted[it] = helix_r ;
      event.helix_dz_inverted[it] = helix_dz;
      event.mom0_x_inverted[it] = mom0.x();
      event.mom0_y_inverted[it] = mom0.y();
      event.mom0_z_inverted[it] = mom0.z();
      event.mom0_inverted[it] = mom0.Mag();
      event.pid_inverted[it] = pid;
      continue;
    }
  }

  Int_t nvtxTpc = TPCAna.GetNVerticesTPC();
  event.nvtxTpc = nvtxTpc;
  event.vtx_x.resize(nvtxTpc);
  event.vtx_y.resize(nvtxTpc);
  event.vtx_z.resize(nvtxTpc);
  event.vtx_dist.resize(nvtxTpc);
  event.vtx_angle.resize(nvtxTpc);
  event.vtxid.resize(nvtxTpc);
  event.vtxmom_theta.resize(nvtxTpc);
  event.vtxpos_x.resize(nvtxTpc);
  event.vtxpos_y.resize(nvtxTpc);
  event.vtxpos_z.resize(nvtxTpc);
  event.vtxmom_x.resize(nvtxTpc);
  event.vtxmom_y.resize(nvtxTpc);
  event.vtxmom_z.resize(nvtxTpc);

  event.isLambda.resize(nvtxTpc);
  event.ncombiLambda.resize(nvtxTpc);
  event.distLambda.resize(nvtxTpc);
  event.angleLambda.resize(nvtxTpc);
  event.bestmassLambda.resize(nvtxTpc);
  event.massLambda.resize(nvtxTpc);
  event.vtxLambda_x.resize(nvtxTpc);
  event.vtxLambda_y.resize(nvtxTpc);
  event.vtxLambda_z.resize(nvtxTpc);
  event.momLambda.resize(nvtxTpc);
  event.momLambda_x.resize(nvtxTpc);
  event.momLambda_y.resize(nvtxTpc);
  event.momLambda_z.resize(nvtxTpc);
  event.decaysidLambda.resize(nvtxTpc);
  event.decaysmomLambda.resize(nvtxTpc);
  event.decaysmomLambda_x.resize(nvtxTpc);
  event.decaysmomLambda_y.resize(nvtxTpc);
  event.decaysmomLambda_z.resize(nvtxTpc);
  for( Int_t it=0; it<nvtxTpc; ++it ){
    TPCVertex *vp = TPCAna.GetTPCVertex( it );
    if( !vp ) continue;
    event.vtx_x[it] = vp -> GetVertex().x();
    event.vtx_y[it] = vp -> GetVertex().y();
    event.vtx_z[it] = vp -> GetVertex().z();
    event.vtx_dist[it] = vp -> GetClosestDist();
    event.vtx_angle[it] = vp -> GetOpeningAngle();

    event.vtxid[it].resize(2);
    event.vtxmom_theta[it].resize(2);
    event.vtxpos_x[it].resize(2);
    event.vtxpos_y[it].resize(2);
    event.vtxpos_z[it].resize(2);
    event.vtxmom_x[it].resize(2);
    event.vtxmom_y[it].resize(2);
    event.vtxmom_z[it].resize(2);

    event.vtxid[it][0] = vp -> GetTrackId(0);
    event.vtxmom_theta[it][0] = vp -> GetTrackTheta(0);
    event.vtxpos_x[it][0] = vp -> GetTrackPos(0).x();
    event.vtxpos_y[it][0] = vp -> GetTrackPos(0).y();
    event.vtxpos_z[it][0] = vp -> GetTrackPos(0).z();
    event.vtxmom_x[it][0] = vp -> GetTrackMom(0).x();
    event.vtxmom_y[it][0] = vp -> GetTrackMom(0).y();
    event.vtxmom_z[it][0] = vp -> GetTrackMom(0).z();

    event.vtxid[it][1] = vp -> GetTrackId(1);
    event.vtxmom_theta[it][1] = vp -> GetTrackTheta(1);
    event.vtxpos_x[it][1] = vp -> GetTrackPos(1).x();
    event.vtxpos_y[it][1] = vp -> GetTrackPos(1).y();
    event.vtxpos_z[it][1] = vp -> GetTrackPos(1).z();
    event.vtxmom_x[it][1] = vp -> GetTrackMom(1).x();
    event.vtxmom_y[it][1] = vp -> GetTrackMom(1).y();
    event.vtxmom_z[it][1] = vp -> GetTrackMom(1).z();

    event.isLambda[it] = vp -> GetIsLambda();
    event.distLambda[it] = vp -> GetClosestDistLambda();
    event.angleLambda[it] = vp -> GetOpeningAngleLambda();
    if(event.isLambda[it]){
      Int_t ncombi = event.ncombiLambda[it] = vp -> GetNcombiLambda();
      Double_t best_lmass = 9999;
      for( Int_t combi=0; combi<ncombi; ++combi ){

	Double_t lmass = vp -> GetMassLambda(combi);
	event.massLambda[it].push_back(lmass);
	Double_t diff = TMath::Abs(lmass - LambdaMass);
	Double_t best_diff = TMath::Abs(best_lmass - LambdaMass);
	if(diff < best_diff) best_lmass = lmass;

	TVector3 vtx = vp -> GetVertexLambda(combi);
	event.vtxLambda_x[it].push_back(vtx.x());
	event.vtxLambda_y[it].push_back(vtx.y());
	event.vtxLambda_z[it].push_back(vtx.z());

	TVector3 lmom = vp -> GetMomLambda(combi);
	event.momLambda[it].push_back(lmom.Mag());
	event.momLambda_x[it].push_back(lmom.x());
	event.momLambda_y[it].push_back(lmom.y());
	event.momLambda_z[it].push_back(lmom.z());

	Int_t pid = vp -> GetProtonIdLambda(combi);
	TVector3 pmom = vp -> GetProtonMomLambda(combi);
	event.decaysidLambda[it].push_back(pid);
	event.decaysmomLambda[it].push_back(pmom.Mag());
	event.decaysmomLambda_x[it].push_back(pmom.x());
	event.decaysmomLambda_y[it].push_back(pmom.y());
	event.decaysmomLambda_z[it].push_back(pmom.z());

	Int_t piid = vp -> GetPionIdLambda(combi);
	TVector3 pimom = vp -> GetPionMomLambda(combi);
	event.decaysidLambda[it].push_back(piid);
	event.decaysmomLambda[it].push_back(pimom.Mag());
	event.decaysmomLambda_x[it].push_back(pimom.x());
	event.decaysmomLambda_y[it].push_back(pimom.y());
	event.decaysmomLambda_z[it].push_back(pimom.z());
      }
      event.bestmassLambda[it] = best_lmass;
    }
  }

  Int_t nvtxTpcClustered = TPCAna.GetNVerticesTPCClustered();
  event.nvtxTpcClustered = nvtxTpcClustered;
  event.Clusteredvtx_x.resize(nvtxTpcClustered);
  event.Clusteredvtx_y.resize(nvtxTpcClustered);
  event.Clusteredvtx_z.resize(nvtxTpcClustered);
  event.Clusteredvtxid.resize(nvtxTpcClustered);
  for( Int_t ivtx=0; ivtx<nvtxTpcClustered; ++ivtx ){
    TPCVertex *vp = TPCAna.GetTPCVertexClustered( ivtx );
    if( !vp ) continue;
    event.Clusteredvtx_x[ivtx] = vp -> GetVertex().x();
    event.Clusteredvtx_y[ivtx] = vp -> GetVertex().y();
    event.Clusteredvtx_z[ivtx] = vp -> GetVertex().z();

    Int_t ntracks = vp -> GetNTracks();
    event.Clusteredvtxid[ivtx].resize(ntracks);
    for( Int_t it=0; it<ntracks; ++it ){
      event.Clusteredvtxid[ivtx][it] = vp -> GetTrackId(it);
    }
  }

#if TrackSearchFailed
  Int_t failed_ntTpc = TPCAna.GetNTracksTPCHelixFailed();
  event.failed_ntTpc = failed_ntTpc;
  event.failed_nhtrack.resize( failed_ntTpc );
  event.failed_flag.resize( failed_ntTpc );
  event.failed_trackid.resize( failed_ntTpc );
  event.failed_isBeam.resize( failed_ntTpc );
  event.failed_isKurama.resize( failed_ntTpc );
  event.failed_isK18.resize( failed_ntTpc );
  event.failed_nclbeforetgt.resize( failed_ntTpc );
  event.failed_isAccidental.resize( failed_ntTpc );
  event.failed_fittime.resize( failed_ntTpc );
  event.failed_searchtime.resize( failed_ntTpc );
  event.failed_niteration.resize( failed_ntTpc );

  event.failed_helix_cx.resize( failed_ntTpc );
  event.failed_helix_cy.resize( failed_ntTpc );
  event.failed_helix_z0.resize( failed_ntTpc );
  event.failed_helix_r.resize( failed_ntTpc );
  event.failed_helix_dz.resize( failed_ntTpc );
  event.failed_mom0.resize( failed_ntTpc );
  event.failed_charge.resize( failed_ntTpc );

  event.failed_hitlayer.resize( failed_ntTpc );
  event.failed_hitpos_x.resize( failed_ntTpc );
  event.failed_hitpos_y.resize( failed_ntTpc );
  event.failed_hitpos_z.resize( failed_ntTpc );
  event.failed_calpos_x.resize( failed_ntTpc );
  event.failed_calpos_y.resize( failed_ntTpc );
  event.failed_calpos_z.resize( failed_ntTpc );
  event.failed_helix_t.resize( failed_ntTpc );
  event.failed_residual.resize( failed_ntTpc );
  event.failed_residual_x.resize( failed_ntTpc );
  event.failed_residual_y.resize( failed_ntTpc );
  event.failed_residual_z.resize( failed_ntTpc );
  event.failed_track_cluster_de.resize( failed_ntTpc );
  event.failed_track_cluster_size.resize( failed_ntTpc );
  event.failed_track_cluster_mrow.resize( failed_ntTpc );

  for( Int_t it=0; it<failed_ntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelixFailed( it );
    if( !tp ) continue;
    Int_t nh = tp->GetNHit();
    Double_t helix_cx=tp->Getcx(), helix_cy=tp->Getcy();
    Double_t helix_z0=tp->Getz0(), helix_r=tp->Getr();
    Double_t helix_dz=tp->Getdz();
    TVector3 mom0 = tp->GetMom0();
    Int_t flag = tp->GetFitFlag();
    Int_t trackid = tp->GetTrackID();
    Int_t isbeam = tp->GetIsBeam();
    Int_t iskurama = tp->GetIsKurama();
    Int_t isk18 = tp->GetIsK18();
    Int_t nclbeforetgt = tp->GetNclBeforeTgt();
    Int_t isaccidental = tp->GetIsAccidental();
    Int_t fittime = tp->GetFitTime();
    Int_t charge = tp->GetCharge();
    Int_t iteration = tp->GetNIteration();
    event.failed_nhtrack[it] = nh;
    event.failed_flag[it] = flag;
    event.failed_trackid[it] = trackid;
    event.failed_isBeam[it] = isbeam;
    event.failed_isKurama[it] = iskurama;
    event.failed_isK18[it] = isk18;
    event.failed_nclbeforetgt[it] = nclbeforetgt;
    event.failed_isAccidental[it] = isaccidental;
    event.failed_fittime[it] = fittime;
    event.failed_searchtime[it] = fittime;
    event.failed_niteration[it] = iteration;

    event.failed_helix_cx[it] = helix_cx;
    event.failed_helix_cy[it] = helix_cy;
    event.failed_helix_z0[it] = helix_z0;
    event.failed_helix_r[it] = helix_r ;
    event.failed_helix_dz[it] = helix_dz;
    event.failed_mom0[it] = mom0.Mag();
    event.failed_charge[it] = charge;

    event.failed_hitlayer[it].resize( nh );
    event.failed_hitpos_x[it].resize( nh );
    event.failed_hitpos_y[it].resize( nh );
    event.failed_hitpos_z[it].resize( nh );
    event.failed_calpos_x[it].resize( nh );
    event.failed_calpos_y[it].resize( nh );
    event.failed_calpos_z[it].resize( nh );
    event.failed_helix_t[it].resize( nh );
    event.failed_residual[it].resize( nh );
    event.failed_residual_x[it].resize( nh );
    event.failed_residual_y[it].resize( nh );
    event.failed_residual_z[it].resize( nh );
    event.failed_track_cluster_de[it].resize( nh );
    event.failed_track_cluster_size[it].resize( nh );
    event.failed_track_cluster_mrow[it].resize( nh );

    for( int ih=0; ih<nh; ++ih ){
      TPCLTrackHit *hit = tp->GetHit( ih );
      if( !hit ) continue;
      const TVector3& hitpos = hit->GetLocalHitPos();
      const TVector3& calpos = hit->GetLocalCalPosHelix();
      const TVector3& res_vect = hit->GetResidualVect();
      Int_t layer = hit->GetLayer();
      Double_t residual = hit->GetResidual();
      Double_t clde = hit->GetDe();
      Double_t mrow = hit->GetMRow();
      TPCCluster *cl = hit->GetHit()->GetParentCluster();
      Int_t clsize = cl->GetClusterSize();

      event.failed_hitlayer[it][ih] = (double)layer;
      event.failed_hitpos_x[it][ih] = hitpos.x();
      event.failed_hitpos_y[it][ih] = hitpos.y();
      event.failed_hitpos_z[it][ih] = hitpos.z();
      event.failed_calpos_x[it][ih] = calpos.x();
      event.failed_calpos_y[it][ih] = calpos.y();
      event.failed_calpos_z[it][ih] = calpos.z();

      event.failed_residual[it][ih] = residual;
      event.failed_residual_x[it][ih] = res_vect.x();
      event.failed_residual_y[it][ih] = res_vect.y();
      event.failed_residual_z[it][ih] = res_vect.z();
      event.failed_track_cluster_de[it][ih] = clde;
      event.failed_track_cluster_size[it][ih] = clsize;
      event.failed_track_cluster_mrow[it][ih] = mrow;
    }
  }
#endif
  HF1( 1, event.status++ );

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::seconds sec = std::chrono::duration_cast<std::chrono::seconds>(end - start);

  HF1( 4, sec.count() );

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
  HB1(1, "Status", 21, 0., 21. );
  HB1(2, "Hough Dist [mm]", 500, 0., 50 );
  HB1(3, "Hough DistY [mm]", 1000, 0., 100 );
  HB1(4, "Process Time [sec]", 1000, 0., 100 );

  HB1(10, "NTrack TPC", 40, 0., 40. );
  HB1(11, "#Hits of Track TPC", 50, 0., 50.);
  HB1(12, "Chisqr TPC", 500, 0., 500.);
  HB1(13, "LayerId TPC", 35, 0., 35.);
  HB1(14, "pHS", 1000, 0., 2.5);
  HB1(15, "mom0", 1000, 0., 2.5);
  HB1(16, "mom0-pHS", 1000, -1.25, 1.25);

  const Int_t nbinpoq = 1000;
  const Double_t minpoq = -2.0;
  const Double_t maxpoq = 2.0;
  const Int_t nbindedx = 1000;
  const Double_t mindedx = 0.;
  const Double_t maxdedx = 350.;

  HB2(20, "<dE/dx>;p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2(21, "<dE/dx>;-p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2(22, "<dE/dx>;+p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HBProf(100, ";#beta#gamma;<-dE/dx> [MeVg^{-1}cm^{2}]", 1000000, 0.1, 10000, 1., 10. );
  HBProf(101, "#pi, CH_{2}; Momentum [GeV/#font[12]{c}]; Stopping Power [MeVcm^{-1}]", 100000, 0.1, maxpoq, 1., 200. );
  HBProf(102, "K, CH_{2}; Momentum [GeV/#font[12]{c}]; Stopping Power [MeVcm^{-1}]", 100000, 0.1, maxpoq, 1., 200. );
  HBProf(103, "p, CH_{2}; Momentum [GeV/#font[12]{c}]; Stopping Power [MeVcm^{-1}]", 100000, 0.1, maxpoq, 1., 200. );
  HBProf(104, "#Xi^{-}, CH_{2}; Momentum [GeV/#font[12]{c}]; Stopping Power [MeVcm^{-1}]", 100000, 0.1, maxpoq, 1., 200. );
  HBProf(105, "#pi, Carbon; Momentum [GeV/#font[12]{c}]; Stopping Power [MeVcm^{-1}]", 100000, 0.1, maxpoq, 1., 200. );
  HBProf(106, "K, Carbon; Momentum [GeV/#font[12]{c}]; Stopping Power [MeVcm^{-1}]", 100000, 0.1, maxpoq, 1., 200. );
  HBProf(107, "p, Carbon; Momentum [GeV/#font[12]{c}]; Stopping Power [MeVcm^{-1}]", 100000, 0.1, maxpoq, 1., 200. );
  HBProf(108, "#Xi^{-}, Carbon; Momentum [GeV/#font[12]{c}]; Stopping Power [MeVcm^{-1}]", 100000, 0.1, maxpoq, 1., 200. );
  const Double_t npoints = 1000000;
  for(Int_t i=0; i<npoints; ++i){
    Double_t x = (Double_t) i/npoints;
    Double_t x1 = TMath::Power(10., 5.*x - 1.);
    Double_t beta = TMath::Sqrt(x1*x1/(x1*x1+1.));
    Double_t Carbon = 3.223;
    HFProf(100, x1, Kinematics::HypTPCdEdx(2, 1000.*pdg::KaonMass(), beta)/Carbon);

    Double_t x2 = TMath::Power(10., (TMath::Log10(maxpoq) - TMath::Log10(0.1))*x + TMath::Log10(0.1));
    HFProf(101, x2, Kinematics::HypTPCdEdx(1, 1000.*pdg::PionMass(), x2/TMath::Sqrt(x2*x2 + pdg::PionMass()*pdg::PionMass())));
    HFProf(102, x2, Kinematics::HypTPCdEdx(1, 1000.*pdg::KaonMass(), x2/TMath::Sqrt(x2*x2 + pdg::KaonMass()*pdg::KaonMass())));
    HFProf(103, x2, Kinematics::HypTPCdEdx(1, 1000.*pdg::ProtonMass(), x2/TMath::Sqrt(x2*x2 + pdg::ProtonMass()*pdg::ProtonMass())));
    HFProf(104, x2, Kinematics::HypTPCdEdx(1, 1000.*pdg::XiMinusMass(), x2/TMath::Sqrt(x2*x2 + pdg::XiMinusMass()*pdg::XiMinusMass())));
    HFProf(105, x2, Kinematics::HypTPCdEdx(2, 1000.*pdg::PionMass(), x2/TMath::Sqrt(x2*x2 + pdg::PionMass()*pdg::PionMass())));
    HFProf(106, x2, Kinematics::HypTPCdEdx(2, 1000.*pdg::KaonMass(), x2/TMath::Sqrt(x2*x2 + pdg::KaonMass()*pdg::KaonMass())));
    HFProf(107, x2, Kinematics::HypTPCdEdx(2, 1000.*pdg::ProtonMass(), x2/TMath::Sqrt(x2*x2 + pdg::ProtonMass()*pdg::ProtonMass())));
    HFProf(108, x2, Kinematics::HypTPCdEdx(2, 1000.*pdg::XiMinusMass(), x2/TMath::Sqrt(x2*x2 + pdg::XiMinusMass()*pdg::XiMinusMass())));
  }

#if TrackClusterHist
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

  for(Int_t layer=0; layer<NumOfLayersBcOut; ++layer){
    HB1(TPCK18RKHid+100+layer, Form("BcOut RK pull, layer%d;pull;Counts",layer), 500, -10, 10);
    HB1(TPCK18RKHid+200+layer, Form("BcOut RK residual, layer%d;X residual [mm];Counts",layer), 500, -10, 10);
  }
  for(Int_t layer=0; layer<10; ++layer){
    HB1(TPCK18RKHid+300+layer, Form("X TPCRK pull, layer%d;X pull;Counts",layer), 500, -10, 10);
    HB1(TPCK18RKHid+400+layer, Form("Y TPCRK pull, layer%d;X pull;Counts",layer), 500, -10, 10);
    HB1(TPCK18RKHid+500+layer, Form("X TPCRK residual, layer%d;Residual [mm];Counts",layer), 500, -10, 10);
    HB1(TPCK18RKHid+600+layer, Form("Y TPCRK residual, layer%d;Residual [mm];Counts",layer), 500, -10, 10);
    HB1(TPCK18RKHid+700+layer, Form("TPC X residual, layer%d;X Residual [mm];Counts",layer), 500, -10, 10);
    HB1(TPCK18RKHid+800+layer, Form("TPC Y residual, layer%d;Y Residual [mm];Counts",layer), 500, -10, 10);
    HB1(TPCK18RKHid+900+layer, Form("TPC Z residual, layer%d;Z Residual [mm];Counts",layer), 500, -10, 10);

    HB1(TPCK18VPHid+layer, Form("TPC - K18 XZ, layer%d;XZ residual [mm];Counts",layer), 500, 0, 100);
    HB1(TPCK18VPHid+100+layer, Form("TPC - K18 X, layer%d;X residual [mm];Counts",layer), 500, -100, 100);
    HB1(TPCK18VPHid+200+layer, Form("TPC - K18 Y, layer%d;Y residual [mm];Counts",layer), 500, -100, 100);
    HB1(TPCK18VPHid+300+layer, Form("TPC - K18 Z, layer%d;Z residual [mm];Counts",layer), 500, -100, 100);
  }
  for(Int_t layer=0; layer<NumOfLayersSdcIn; ++layer){
    HB1(TPCKuramaRKHid+100+layer, Form("Q>0, SdcIn RK pull, layer%d;pull;Counts",layer), 500, -10, 10);
    HB1(TPCKuramaRKHid+200+layer, Form("Q>0, SdcIn RK residual, layer%d;Residual [mm];Counts",layer), 500, -10, 10);

    HB1(TPCKuramaRKHid+1100+layer, Form("Q<0, SdcIn RK pull, layer%d;pull;Counts",layer), 500, -10, 10);
    HB1(TPCKuramaRKHid+1200+layer, Form("Q<0, SdcIn RK residual, layer%d;Residual [mm];Counts",layer), 500, -10, 10);
  }
  for(Int_t layer=0; layer<NumOfLayersSdcOut; ++layer){
    HB1(TPCKuramaRKHid+300+layer, Form("Q>0, SdcOut RK pull, layer%d;pull;Counts",layer), 500, -10, 10);
    HB1(TPCKuramaRKHid+400+layer, Form("Q>0, SdcOut RK residual, layer%d;residual [mm];Counts",layer), 500, -10, 10);

    HB1(TPCKuramaRKHid+1300+layer, Form("Q<0, SdcOut RK pull, layer%d;pull;Counts",layer), 500, -10, 10);
    HB1(TPCKuramaRKHid+1400+layer, Form("Q<0, SdcOut RK residual, layer%d;residual [mm];Counts",layer), 500, -10, 10);
  }
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    HB1(TPCKuramaRKHid+500+layer, Form("Q>0, X TPCRK pull, layer%d;X pull;Counts",layer), 500, -10, 10);
    HB1(TPCKuramaRKHid+600+layer, Form("Q>0, Y TPCRK pull, layer%d;X pull;Counts",layer), 500, -10, 10);
    HB1(TPCKuramaRKHid+700+layer, Form("Q>0, X TPCRK residual, layer%d;X residual [mm];Counts",layer), 500, -10, 10);
    HB1(TPCKuramaRKHid+800+layer, Form("Q>0, Y TPCRK residual, layer%d;Y residual [mm];Counts",layer), 500, -10, 10);
    HB1(TPCKuramaRKHid+1500+layer, Form("Q<0, X TPCRK pull, layer%d;X pull;Counts",layer), 500, -10, 10);
    HB1(TPCKuramaRKHid+1600+layer, Form("Q<0, Y TPCRK pull, layer%d;X pull;Counts",layer), 500, -10, 10);
    HB1(TPCKuramaRKHid+1700+layer, Form("Q<0, X TPCRK residual, layer%d;X residual [mm];Counts",layer), 500, -10, 10);
    HB1(TPCKuramaRKHid+1800+layer, Form("Q<0, Y TPCRK residual, layer%d;Y residual [mm];Counts",layer), 500, -10, 10);

    HB1(TPCKuramaRKHid+4100+layer, Form("Q>0, TPC X residual, layer%d;X Residual [mm];Counts",layer), 500, -10, 10);
    HB1(TPCKuramaRKHid+4200+layer, Form("Q>0, TPC Y residual, layer%d;Y Residual [mm];Counts",layer), 500, -10, 10);
    HB1(TPCKuramaRKHid+4300+layer, Form("Q>0, TPC Z residual, layer%d;Z Residual [mm];Counts",layer), 500, -10, 10);
    HB1(TPCKuramaRKHid+5100+layer, Form("Q<0, TPC X residual, layer%d;X Residual [mm];Counts",layer), 500, -10, 10);
    HB1(TPCKuramaRKHid+5200+layer, Form("Q<0, TPC Y residual, layer%d;Y Residual [mm];Counts",layer), 500, -10, 10);
    HB1(TPCKuramaRKHid+5300+layer, Form("Q<0, TPC Z residual, layer%d;Z Residual [mm];Counts",layer), 500, -10, 10);
  }
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    HB1(TPCKuramaVPHid+layer, Form("Q>0, TPC - Kurama XZ, layer%d;XZ residual [mm];Counts",layer), 500, 0, 100);
    HB1(TPCKuramaVPHid+100+layer, Form("Q>0, TPC - Kurama X, layer%d;X residual [mm];Counts",layer), 500, -100, 100);
    HB1(TPCKuramaVPHid+200+layer, Form("Q>0, TPC - Kurama Y, layer%d;Y residual [mm];Counts",layer), 500, -100, 100);
    HB1(TPCKuramaVPHid+300+layer, Form("Q>0, TPC - Kurama Z, layer%d;Z residual [mm];Counts",layer), 500, -100, 100);
    HB1(TPCKuramaVPHid+1000+layer, Form("Q<0, TPC - Kurama XZ, layer%d;XZ residual [mm];Counts",layer), 500, 0, 100);
    HB1(TPCKuramaVPHid+1100+layer, Form("Q<0, TPC - Kurama X, layer%d;X residual [mm];Counts",layer), 500, -100, 100);
    HB1(TPCKuramaVPHid+1200+layer, Form("Q<0, TPC - Kurama Y, layer%d;Y residual [mm];Counts",layer), 500, -100, 100);
    HB1(TPCKuramaVPHid+1300+layer, Form("Q<0, TPC - Kurama Z, layer%d;Z residual [mm];Counts",layer), 500, -100, 100);
  }

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    HB1(TPCExclusiveHid+layer,	Form("TPC ExclusiveResidual T[mm];layer%d",layer),1000,-10,10);
    HB1(TPCExclusiveHid+100+layer,	Form("TPC ExclusiveResidual X[mm];layer%d",layer),1000,-10,10);
    HB1(TPCExclusiveHid+200+layer,	Form("TPC ExclusiveResidual Y[mm];layer%d",layer),1000,-10,10);
    HB1(TPCExclusiveHid+300+layer,	Form("TPC ExclusiveResidual Z[mm];layer%d",layer),1000,-10,10);
    HB1(TPCExclusiveHid+1000+layer,	Form("TPC ExclusivePull T;layer%d",layer),1000,-10,10);
    HB1(TPCExclusiveHid+1000+100+layer,	Form("TPC ExclusivePull X;layer%d",layer),1000,-10,10);
    HB1(TPCExclusiveHid+1000+200+layer,	Form("TPC ExclusivePull Y;layer%d",layer),1000,-10,10);
    HB1(TPCExclusiveHid+1000+300+layer,	Form("TPC ExclusivePull Z;layer%d",layer),1000,-10,10);

    HB1(TPCIntrinsicHid+layer,	Form("TPC IntrinsicResidual T[mm];layer%d",layer),1000,-10,10);
    HB1(TPCIntrinsicHid+100+layer,	Form("TPC IntrinsicResidual X[mm];layer%d",layer),1000,-10,10);
    HB1(TPCIntrinsicHid+200+layer,	Form("TPC IntrinsicResidual Y[mm];layer%d",layer),1000,-10,10);
    HB1(TPCIntrinsicHid+300+layer,	Form("TPC IntrinsicResidual Z[mm];layer%d",layer),1000,-10,10);
    HB1(TPCIntrinsicHid+1000+layer,	Form("TPC IntrinsicPull T;layer%d",layer),1000,-10,10);
    HB1(TPCIntrinsicHid+1000+100+layer,	Form("TPC IntrinsicPull X;layer%d",layer),1000,-10,10);
    HB1(TPCIntrinsicHid+1000+200+layer,	Form("TPC IntrinsicPull Y;layer%d",layer),1000,-10,10);
    HB1(TPCIntrinsicHid+1000+300+layer,	Form("TPC IntrinsicPull Z;layer%d",layer),1000,-10,10);
  }

  HB1(TPCInclusiveHid+32,	Form("TPC InclusiveResidual T[mm];AllLayer"),1000,-10,10);
  HB1(TPCInclusiveHid+100+32,	Form("TPC InclusiveResidual X[mm];AllLayer"),1000,-10,10);
  HB1(TPCInclusiveHid+200+32,	Form("TPC InclusiveResidual Y[mm];AllLayer"),1000,-10,10);
  HB1(TPCInclusiveHid+300+32,	Form("TPC InclusiveResidual Z[mm];AllLayer"),1000,-10,10);
  HB1(TPCInclusiveHid+1000+32,	Form("TPC InclusivePull T;AllLayer"),1000,-10,10);
  HB1(TPCInclusiveHid+1000+100+32,	Form("TPC InclusivePull X;AllLayer"),1000,-10,10);
  HB1(TPCInclusiveHid+1000+200+32,	Form("TPC InclusivePull Y;AllLayer"),1000,-10,10);
  HB1(TPCInclusiveHid+1000+300+32,	Form("TPC InclusivePull Z;AllLayer"),1000,-10,10);

  HB1(TPCExclusiveHid+32,	Form("TPC ExclusiveResidual T[mm];AllLayer"),1000,-10,10);
  HB1(TPCExclusiveHid+100+32,	Form("TPC ExclusiveResidual X[mm];AllLayer"),1000,-10,10);
  HB1(TPCExclusiveHid+200+32,	Form("TPC ExclusiveResidual Y[mm];AllLayer"),1000,-10,10);
  HB1(TPCExclusiveHid+300+32,	Form("TPC ExclusiveResidual Z[mm];AllLayer"),1000,-10,10);
  HB1(TPCExclusiveHid+1000+32,	Form("TPC ExclusivePull T;AllLayer"),1000,-10,10);
  HB1(TPCExclusiveHid+1000+100+32,	Form("TPC ExclusivePull X;AllLayer"),1000,-10,10);
  HB1(TPCExclusiveHid+1000+200+32,	Form("TPC ExclusivePull Y;AllLayer"),1000,-10,10);
  HB1(TPCExclusiveHid+1000+300+32,	Form("TPC ExclusivePull Z;AllLayer"),1000,-10,10);

  HB2(TPCInclusiveHid+10000,	Form("TPC InclusiveResidual T[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCInclusiveHid+100+10000,	Form("TPC InclusiveResidual X[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCInclusiveHid+200+10000,	Form("TPC InclusiveResidual Y[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCInclusiveHid+300+10000,	Form("TPC InclusiveResidual Z[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCInclusiveHid+1000+10000,	Form("TPC InclusivePull T:Layer"),31,0,31,1000,-10,10);
  HB2(TPCInclusiveHid+1000+100+10000,	Form("TPC InclusivePull X:Layer"),31,0,31,1000,-10,10);
  HB2(TPCInclusiveHid+1000+200+10000,	Form("TPC InclusivePull Y:Layer"),31,0,31,1000,-10,10);
  HB2(TPCInclusiveHid+1000+300+10000,	Form("TPC InclusivePull Z:Layer"),31,0,31,1000,-10,10);

  HB2(TPCExclusiveHid+10000,	Form("TPC ExclusiveResidual T[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCExclusiveHid+100+10000,	Form("TPC ExclusiveResidual X[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCExclusiveHid+200+10000,	Form("TPC ExclusiveResidual Y[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCExclusiveHid+300+10000,	Form("TPC ExclusiveResidual Z[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCExclusiveHid+1000+10000,	Form("TPC ExclusivePull T:Layer"),31,0,31,1000,-10,10);
  HB2(TPCExclusiveHid+1000+100+10000,	Form("TPC ExclusivePull X:Layer"),31,0,31,1000,-10,10);
  HB2(TPCExclusiveHid+1000+200+10000,	Form("TPC ExclusivePull Y:Layer"),31,0,31,1000,-10,10);
  HB2(TPCExclusiveHid+1000+300+10000,	Form("TPC ExclusivePull Z:Layer"),31,0,31,1000,-10,10);

  HB2(TPCIntrinsicHid+10000,	Form("TPC IntrinsicResidual T[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCIntrinsicHid+100+10000,	Form("TPC IntrinsicResidual X[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCIntrinsicHid+200+10000,	Form("TPC IntrinsicResidual Y[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCIntrinsicHid+300+10000,	Form("TPC IntrinsicResidual Z[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCIntrinsicHid+1000+10000,	Form("TPC IntrinsicPull T:Layer"),31,0,31,1000,-10,10);
  HB2(TPCIntrinsicHid+1000+100+10000,	Form("TPC IntrinsicPull X:Layer"),31,0,31,1000,-10,10);
  HB2(TPCIntrinsicHid+1000+200+10000,	Form("TPC IntrinsicPull Y:Layer"),31,0,31,1000,-10,10);
  HB2(TPCIntrinsicHid+1000+300+10000,	Form("TPC IntrinsicPull Z:Layer"),31,0,31,1000,-10,10);

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    HB2(TPCKuramaRKHid+layer+6100, Form("Q>0, X residual, layer%d;(Kurama track - TPC cluster) X [mm];(TPCKurama track - TPC cluster) X [mm]",layer), 1000, -20, 20, 1000, -20, 20);
    HB2(TPCKuramaRKHid+layer+6200, Form("Q>0, X residual, layer%d;(Kurama track - TPC cluster) X [mm];(TPC track - TPC cluster) X [mm]",layer), 1000, -20, 20, 1000, -20, 20);
    HB2(TPCKuramaRKHid+layer+6300, Form("Q>0, X residual, layer%d;(TPC track - TPC cluster) X [mm];(TPCKurama track - TPC cluster) X [mm]",layer), 1000, -20, 20, 1000, -20, 20);
    HB2(TPCKuramaRKHid+layer+6400, Form("Q>0, Y residual, layer%d;(Kurama track - TPC cluster) Y [mm];(TPCKurama track - TPC cluster) Y [mm]",layer), 1000, -20, 20, 1000, -20, 20);
    HB2(TPCKuramaRKHid+layer+6500, Form("Q>0, Y residual, layer%d;(Kurama track - TPC cluster) Y [mm];(TPC track - TPC cluster) Y [mm]",layer), 1000, -20, 20, 1000, -20, 20);
    HB2(TPCKuramaRKHid+layer+6600, Form("Q>0, Y residual, layer%d;(TPC track - TPC cluster) Y [mm];(TPCKurama track - TPC cluster) Y [mm]",layer), 1000, -20, 20, 1000, -20, 20);

    HB2(TPCKuramaRKHid+layer+7100, Form("Q<0, X residual, layer%d;(Kurama track - TPC cluster) X [mm];(TPCKurama track - TPC cluster) X [mm]",layer), 1000, -20, 20, 1000, -20, 20);
    HB2(TPCKuramaRKHid+layer+7200, Form("Q<0, X residual, layer%d;(Kurama track - TPC cluster) X [mm];(TPC track - TPC cluster) X [mm]",layer), 1000, -20, 20, 1000, -20, 20);
    HB2(TPCKuramaRKHid+layer+7300, Form("Q<0, X residual, layer%d;(TPC track - TPC cluster) X [mm];(TPCKurama track - TPC cluster) X [mm]",layer), 1000, -20, 20, 1000, -20, 20);
    HB2(TPCKuramaRKHid+layer+7400, Form("Q<0, Y residual, layer%d;(Kurama track - TPC cluster) Y [mm];(TPCKurama track - TPC cluster) Y [mm]",layer), 1000, -20, 20, 1000, -20, 20);
    HB2(TPCKuramaRKHid+layer+7500, Form("Q<0, Y residual, layer%d;(Kurama track - TPC cluster) Y [mm];(TPC track - TPC cluster) Y [mm]",layer), 1000, -20, 20, 1000, -20, 20);
    HB2(TPCKuramaRKHid+layer+7600, Form("Q<0, Y residual, layer%d;(TPC track - TPC cluster) Y [mm];(TPCKurama track - TPC cluster) Y [mm]",layer), 1000, -20, 20, 1000, -20, 20);
  }

  HB1(5001, "P K18", 800, 1.4, 2.2);
  HB1(5002, "P Kurama", 600, 0, 3);
  HB1(5003, "Charge*Mass", 600, -1., 2.);
  HB2(5004, "P Kurama%Charge*Mass", 400, -0.8, 2.2, 400, 0, 2.5);
  HB2(5010, "Vertex", 400, -500, 500, 400, -150, 150);
  HB1(5011, "Closest distance", 200, 0., 200.);
  HB1(5012, "Vertex X", 400, -200, 200);
  HB1(5013, "Vertex Y", 400, -200, 200);
  HB1(5014, "Vertex Z", 500, -500, 500);

  HB1(6001, "P K18", 800, 1.4, 2.2);
  HB1(6002, "P Kurama", 600, 0, 3);
  HB1(6003, "Charge*Mass", 600, -1., 2.);
  HB2(6004, "P Kurama%Charge*Mass", 400, -0.8, 2.2, 400, 0, 2.5);
  HB2(6010, "Vertex", 400, -500, 500, 400, -150, 150);
  HB1(6011, "Closest distance", 50, 0., 50.);
  HB1(6012, "Vertex X", 200, -100, 100);
  HB1(6013, "Vertex Y", 200, -100, 100);
  HB1(6014, "Vertex Z", 200, -200, 200);

  HB1(7001, "P K18", 800, 1.4, 2.2);
  HB1(7002, "P Kurama", 600, 0, 3);
  HB1(7003, "Charge*Mass", 600, -1., 2.);
  HB2(7004, "P Kurama%Charge*Mass", 400, -0.8, 2.2, 400, 0, 2.5);
  HB1(7010, "Closest distance", 50, 0., 50.);
  HB1(7011, "Vertex X", 200, -100, 100);
  HB1(7012, "Vertex Y", 200, -100, 100);
  HB1(7013, "Vertex Z", 200, -200, 200);
  HB1(7014, "MissingMass K^{+}", 140, 1.1, 1.8);
  HB1(7015, "MissingMass K^{+}", 140, 1.1, 1.8);
  HB2(7016, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(7017, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(7018, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(7019, "MissingMass%U ", 160, -0.4, 0.4, 140, 1.1, 1.8);
  HB2(7020, "MissingMass%U ", 160, -0.4, 0.4, 140, 1.1, 1.8);
  HB2(7021, "MissingMass%U ", 160, -0.4, 0.4, 160, 1.1, 1.8);
  HB2(7022, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(7023, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(7024, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(7025, "MissingMass%V ", 160, -0.4, 0.4, 140, 1.1, 1.8);
  HB2(7026, "MissingMass%V ", 160, -0.4, 0.4, 140, 1.1, 1.8);
  HB2(7027, "MissingMass%V ", 160, -0.4, 0.4, 140, 1.1, 1.8);
  HB2(7028, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(7029, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(7030, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(7031, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 140, 1.1, 1.8);
  HB2(7032, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 140, 1.1, 1.8);
  HB2(7033, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 140, 1.1, 1.8);
  HB2(7036, "#DeltaP%P_{calc.} ", 160, 0.3, 1.1, 200, -1., 1.);
  HB2(7039, "MissingMass%P_{calc.} ", 160, 0.3, 1.1, 140, 1.1, 1.8);

  HBProf(7116, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(7117, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(7118, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(7119, "MissingMass%U Prof ", 160, -0.4, 0.4, 1.1, 1.8);
  HBProf(7120, "MissingMass%U Prof ", 160, -0.4, 0.4, 1.1, 1.8);
  HBProf(7121, "MissingMass%U Prof ", 160, -0.4, 0.4, 1.1, 1.8);
  HBProf(7122, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(7123, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(7124, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(7125, "MissingMass%V Prof ", 160, -0.4, 0.4, 1.1, 1.8);
  HBProf(7126, "MissingMass%V Prof ", 160, -0.4, 0.4, 1.1, 1.8);
  HBProf(7127, "MissingMass%V Prof ", 160, -0.4, 0.4, 1.1, 1.8);
  HBProf(7128, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(7129, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(7130, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(7131, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1.1, 1.8);
  HBProf(7132, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1.1, 1.8);
  HBProf(7133, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1.1, 1.8);
  HBProf(7136, "#DeltaP%P_{calc.} Prof ", 160, 0.3, 1.1, -1., 1.);
  HBProf(7139, "MissingMass%P_{calc.} Prof ", 160, 0.3, 1.1, 1.1, 1.8);

  HB1(8001, "P K18", 800, 1.4, 2.2);
  HB1(8002, "P Kurama", 600, 0, 3);
  HB1(8003, "Charge*Mass", 600, -1., 2.);
  HB2(8004, "P Kurama%Charge*Mass", 400, -0.8, 2.2, 400, 0, 2.5);
  HB1(8010, "Closest distance", 50, 0., 50.);
  HB1(8011, "Vertex X", 200, -100, 100);
  HB1(8012, "Vertex Y", 200, -100, 100);
  HB1(8013, "Vertex Z", 200, -200, 200);
  HB1(8014, "MissingMass K^{+}", 140, 1.1, 1.8);
  HB1(8015, "MissingMass K^{+}", 140, 1.1, 1.8);
  HB2(8016, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8017, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8018, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8019, "MissingMass%U ", 160, -0.4, 0.4, 140, 1.1, 1.8);
  HB2(8020, "MissingMass%U ", 160, -0.4, 0.4, 140, 1.1, 1.8);
  HB2(8021, "MissingMass%U ", 160, -0.4, 0.4, 140, 1.1, 1.8);
  HB2(8022, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8023, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8024, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8025, "MissingMass%V ", 160, -0.4, 0.4, 140, 1.1, 1.8);
  HB2(8026, "MissingMass%V ", 160, -0.4, 0.4, 140, 1.1, 1.8);
  HB2(8027, "MissingMass%V ", 160, -0.4, 0.4, 140, 1.1, 1.8);
  HB2(8028, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(8029, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(8030, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(8031, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 140, 1.1, 1.8);
  HB2(8032, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 140, 1.1, 1.8);
  HB2(8033, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 140, 1.1, 1.8);
  HB2(8036, "#DeltaP%P_{calc.} ", 160, 0.3, 1.1, 200, -1., 1.);
  HB2(8039, "MissingMass%P_{calc.} ", 160, 0.3, 1.1, 140, 1.1, 1.8);

  HBProf(8116, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8117, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8118, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8119, "MissingMass%U Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(8120, "MissingMass%U Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(8121, "MissingMass%U Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(8122, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(8123, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(8124, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(8125, "MissingMass%V Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(8126, "MissingMass%V Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(8127, "MissingMass%V Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(8128, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(8129, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(8130, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(8131, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1., 1.8);
  HBProf(8132, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1., 1.8);
  HBProf(8133, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1., 1.8);
  HBProf(8136, "#DeltaP%P_{calc.} Prof ", 160, 0.3, 1.1, -1., 1.);
  HBProf(8139, "MissingMass%P_{calc.} Prof ", 160, 0.3, 1.1, 1.1, 1.8);

  HB1(8201, "P K18", 800, 1.4, 2.2);
  HB1(8202, "P Kurama", 600, 0, 3);
  HB1(8203, "Charge*Mass", 600, -1., 2.);
  HB2(8204, "P Kurama%Charge*Mass", 400, -0.8, 2.2, 400, 0, 2.5);
  HB1(8210, "Closest distance", 50, 0., 50.);
  HB1(8211, "Vertex X", 200, -100, 100);
  HB1(8212, "Vertex Y", 200, -100, 100);
  HB1(8213, "Vertex Z", 200, -200, 200);
  HB1(8214, "MissingMass p", 200, 0.3, 1.3);
  HB1(8215, "MissingMass p", 200, 0.3, 1.3);
  HB2(8216, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8217, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8218, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8219, "MissingMass%U ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(8220, "MissingMass%U ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(8221, "MissingMass%U ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(8222, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8223, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8224, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8225, "MissingMass%V ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(8226, "MissingMass%V ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(8227, "MissingMass%V ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(8228, "#DeltaP%P_{calc.} ", 440, 0.3, 2.5, 200, -1., 1.);
  HB2(8229, "#DeltaP%P_{calc.} ", 440, 0.3, 2.5, 200, -1., 1.);
  HB2(8230, "#DeltaP%P_{calc.} ", 440, 0.3, 2.5, 200, -1., 1.);
  HB2(8231, "MissingMass%P_{calc.} ", 440, 0.3, 2.5, 200, 0.3, 1.3);
  HB2(8232, "MissingMass%P_{calc.} ", 440, 0.3, 2.5, 200, 0.3, 1.3);
  HB2(8233, "MissingMass%P_{calc.} ", 440, 0.3, 2.5, 200, 0.3, 1.3);

  HBProf(8316, "#DeltaP%U ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8317, "#DeltaP%U ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8318, "#DeltaP%U ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8319, "MissingMass%U ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(8320, "MissingMass%U ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(8321, "MissingMass%U ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(8322, "#DeltaP%V ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8323, "#DeltaP%V ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8324, "#DeltaP%V ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8325, "MissingMass%V ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(8326, "MissingMass%V ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(8327, "MissingMass%V ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(8328, "#DeltaP%P_{calc.} ", 440, 0.3, 2.5, -1., 1.);
  HBProf(8329, "#DeltaP%P_{calc.} ", 440, 0.3, 2.5, -1., 1.);
  HBProf(8330, "#DeltaP%P_{calc.} ", 440, 0.3, 2.5, -1., 1.);
  HBProf(8331, "MissingMass%P_{calc.} ", 440, 0.3, 2.5, 0.3, 1.3);
  HBProf(8332, "MissingMass%P_{calc.} ", 440, 0.3, 2.5, 0.3, 1.3);
  HBProf(8333, "MissingMass%P_{calc.} ", 440, 0.3, 2.5, 0.3, 1.3);

  TString name[3] = { "Pion", "Kaon", "Proton"};
  for(Int_t i=0; i<NumOfSegTOF; ++i){
    for(Int_t ip=0; ip<3; ++ip){
      HB1(10000+(i+1)*100+ip+1,
	  Form("%s, TofTime-%sTime", name[ip].Data(), name[ip].Data()),
	  200, -3., 3.);
    }
  }
#endif

  HBTree( "tpc", "tree of DstTPCKuramaK18Tracking" );
  tree->Branch( "status", &event.status );
  tree->Branch( "runnum", &event.runnum );
  tree->Branch( "evnum", &event.evnum );
  tree->Branch( "trigpat", &event.trigpat );
  tree->Branch( "trigflag", &event.trigflag );
  tree->Branch( "clkTpc", &event.clkTpc);

  tree->Branch( "ntK18",      &event.ntK18);
  tree->Branch( "chisqrK18",  &event.chisqrK18);
  tree->Branch( "p_3rd",       &event.p_3rd);
  tree->Branch( "delta_3rd",  &event.delta_3rd);
  tree->Branch( "pHS",        &event.pHS);

  tree->Branch( "xoutK18",    &event.xoutK18);
  tree->Branch( "youtK18",    &event.youtK18);
  tree->Branch( "uoutK18",    &event.uoutK18);
  tree->Branch( "voutK18",    &event.voutK18);
  tree->Branch( "initmomHS", &event.initmomHS);
  tree->Branch( "xtgtHS",    &event.xtgtHS);
  tree->Branch( "ytgtHS",    &event.ytgtHS);
  tree->Branch( "ztgtHS",    &event.ztgtHS);
  tree->Branch( "utgtHS",    &event.utgtHS);
  tree->Branch( "vtgtHS",    &event.vtgtHS);
  tree->Branch( "xbcHS",    &event.xbcHS);
  tree->Branch( "ybcHS",    &event.ybcHS);
  tree->Branch( "zbcHS",    &event.zbcHS);
  tree->Branch( "ubcHS",    &event.ubcHS);
  tree->Branch( "vbcHS",    &event.vbcHS);
  tree->Branch( "xbh2HS",    &event.xbh2HS);
  tree->Branch( "ybh2HS",    &event.ybh2HS);
  tree->Branch( "zbh2HS",    &event.zbh2HS);
  tree->Branch( "ubh2HS",    &event.ubh2HS);
  tree->Branch( "vbh2HS",    &event.vbh2HS);
  tree->Branch( "pbh2HS",    &event.pbh2HS);
  tree->Branch( "xgasvesselHS",    &event.xgasvesselHS);
  tree->Branch( "ygasvesselHS",    &event.ygasvesselHS);
  tree->Branch( "zgasvesselHS",    &event.zgasvesselHS);
  tree->Branch( "ugasvesselHS",    &event.ugasvesselHS);
  tree->Branch( "vgasvesselHS",    &event.vgasvesselHS);
  tree->Branch( "pgasvesselHS",    &event.pgasvesselHS);
  tree->Branch( "xhtofHS",    &event.xhtofHS);
  tree->Branch( "yhtofHS",    &event.yhtofHS);
  tree->Branch( "zhtofHS",    &event.zhtofHS);
  tree->Branch( "uhtofHS",    &event.uhtofHS);
  tree->Branch( "vhtofHS",    &event.vhtofHS);
  tree->Branch( "xvpHS",    &event.xvpHS);
  tree->Branch( "yvpHS",    &event.yvpHS);
  tree->Branch( "zvpHS",    &event.zvpHS);
  tree->Branch( "uvpHS",    &event.uvpHS);
  tree->Branch( "vvpHS",    &event.vvpHS);
  tree->Branch("layerK18", &event.layerK18);
  tree->Branch("wireK18", &event.wireK18);
  tree->Branch("localhitposK18", &event.localhitposK18);
  tree->Branch("wposK18", &event.wposK18);

  tree->Branch( "tpcidTPCK18", &event.tpcidTPCK18);
  tree->Branch( "isgoodTPCK18", &event.isgoodTPCK18);
  tree->Branch( "niterationTPCK18", &event.niterationTPCK18);
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

  tree->Branch( "ntKurama",      &event.ntKurama);
  tree->Branch( "chisqrKurama",  &event.chisqrKurama);
  tree->Branch( "pKurama",       &event.pKurama);
  tree->Branch( "qKurama",       &event.qKurama);
  tree->Branch( "xout",          &event.xout);
  tree->Branch( "yout",          &event.yout);
  tree->Branch( "zout",          &event.zout);
  tree->Branch( "pxout",         &event.pxout);
  tree->Branch( "pyout",         &event.pyout);
  tree->Branch( "pzout",         &event.pzout);
  tree->Branch( "layer",         &event.layer);
  tree->Branch( "wire",          &event.wire);
  tree->Branch( "localhitpos",   &event.localhitpos);
  tree->Branch( "wpos",          &event.wpos);

  tree->Branch( "xtgtKurama",    &event.xtgtKurama);
  tree->Branch( "ytgtKurama",    &event.ytgtKurama);
  tree->Branch( "utgtKurama",    &event.utgtKurama);
  tree->Branch( "vtgtKurama",    &event.vtgtKurama);
  tree->Branch( "thetaKurama",   &event.thetaKurama);
  tree->Branch( "phiKurama",     &event.phiKurama);
  tree->Branch( "xvpKurama",     &event.xvpKurama);
  tree->Branch( "yvpKurama",     &event.yvpKurama);
  tree->Branch( "zvpKurama",     &event.zvpKurama);
  tree->Branch( "xhtofKurama",   &event.xhtofKurama);
  tree->Branch( "yhtofKurama",   &event.yhtofKurama);
  tree->Branch( "zhtofKurama",   &event.zhtofKurama);
  tree->Branch( "tpcidTPCKurama", &event.tpcidTPCKurama);
  tree->Branch( "isgoodTPCKurama", &event.isgoodTPCKurama);
  tree->Branch( "niterationTPCKurama", &event.niterationTPCKurama);
  tree->Branch( "kflagTPCKurama", &event.kflagTPCKurama);
  tree->Branch( "pflagTPCKurama", &event.pflagTPCKurama);
  tree->Branch( "chisqrTPCKurama", &event.chisqrTPCKurama);
  tree->Branch( "pTPCKurama", &event.pTPCKurama);
  tree->Branch( "qTPCKurama", &event.qTPCKurama);
  tree->Branch( "m2TPCKurama", &event.m2TPCKurama);
  tree->Branch( "m2OrgTPCKurama", &event.m2OrgTPCKurama);
  tree->Branch( "xtgtTPCKurama", &event.xtgtTPCKurama);
  tree->Branch( "ytgtTPCKurama", &event.ytgtTPCKurama);
  tree->Branch( "utgtTPCKurama", &event.utgtTPCKurama);
  tree->Branch( "vtgtTPCKurama", &event.vtgtTPCKurama);
  tree->Branch( "thetaTPCKurama", &event.thetaTPCKurama);
  tree->Branch( "pathTPCKurama", &event.pathTPCKurama);
  tree->Branch( "lhtofTPCKurama", &event.lhtofTPCKurama);
  tree->Branch( "xhtofTPCKurama", &event.xhtofTPCKurama);
  tree->Branch( "yhtofTPCKurama", &event.yhtofTPCKurama);
  tree->Branch( "phtofTPCKurama", &event.phtofTPCKurama);
  tree->Branch( "lgasvesselTPCKurama", &event.lgasvesselTPCKurama);
  tree->Branch( "xgasvesselTPCKurama", &event.xgasvesselTPCKurama);
  tree->Branch( "ygasvesselTPCKurama", &event.ygasvesselTPCKurama);
  tree->Branch( "pgasvesselTPCKurama", &event.pgasvesselTPCKurama);
  tree->Branch( "lvpTPCKurama", &event.lvpTPCKurama);
  tree->Branch( "xvpTPCKurama", &event.xvpTPCKurama);
  tree->Branch( "yvpTPCKurama", &event.yvpTPCKurama);

  tree->Branch( "isgoodTPC", &event.isgoodTPC);
  tree->Branch( "insideTPC", &event.insideTPC);
  tree->Branch( "vtxTPC", &event.vtxTPC);
  tree->Branch( "vtyTPC", &event.vtyTPC);
  tree->Branch( "vtzTPC", &event.vtzTPC);
  tree->Branch( "pxKmTPC", &event.pxKmTPC);
  tree->Branch( "pyKmTPC", &event.pyKmTPC);
  tree->Branch( "pzKmTPC", &event.pzKmTPC);
  tree->Branch( "pxScatTPC", &event.pxScatTPC);
  tree->Branch( "pyScatTPC", &event.pyScatTPC);
  tree->Branch( "pzScatTPC", &event.pzScatTPC);
  tree->Branch( "closeDistTPC", &event.closeDistTPC);
  tree->Branch( "MissMassTPC", &event.MissMassTPC);
  tree->Branch( "MissMassCorrTPC", &event.MissMassCorrTPC);
  tree->Branch( "MissMassCorrDETPC", &event.MissMassCorrDETPC);
  tree->Branch( "MissMassNuclTPC", &event.MissMassNuclTPC);
  tree->Branch( "MissMassNuclCorrTPC", &event.MissMassNuclCorrTPC);
  tree->Branch( "MissMassNuclCorrDETPC", &event.MissMassNuclCorrDETPC);
  tree->Branch( "pOrgTPC", &event.pOrgTPC);
  tree->Branch( "pCorrTPC", &event.pCorrTPC);
  tree->Branch( "pCorrDETPC", &event.pCorrDETPC);
  tree->Branch( "pCalcTPC", &event.pCalcTPC);
  tree->Branch( "thetaCMTPC", &event.thetaCMTPC);
  tree->Branch( "costCMTPC", &event.costCMTPC);
  tree->Branch( "pCalcDETPC", &event.pCalcDETPC);
  tree->Branch( "thetaCMDETPC", &event.thetaCMDETPC);
  tree->Branch( "costCMDETPC", &event.costCMDETPC);
  tree->Branch( "xistarpCalcDETPC", &event.xistarpCalcDETPC);
  tree->Branch( "xistarthetaCMDETPC", &event.xistarthetaCMDETPC);
  tree->Branch( "xistarcostCMDETPC", &event.xistarcostCMDETPC);
  tree->Branch( "kpscatpCalcTPC", &event.kpscatpCalcTPC);
  tree->Branch( "kpscatthetaCMTPC", &event.kpscatthetaCMTPC);
  tree->Branch( "kpscatcostCMTPC", &event.kpscatcostCMTPC);
  tree->Branch( "kpscatpCalcDETPC", &event.kpscatpCalcDETPC);
  tree->Branch( "kpscatthetaCMDETPC", &event.kpscatthetaCMDETPC);
  tree->Branch( "kpscatcostCMDETPC", &event.kpscatcostCMDETPC);
  tree->Branch( "thetaTPC", &event.thetaTPC);
  tree->Branch( "ubTPC", &event.ubTPC);
  tree->Branch( "vbTPC", &event.vbTPC);
  tree->Branch( "usTPC", &event.usTPC);
  tree->Branch( "vsTPC", &event.vsTPC);

  tree->Branch( "nhTpc", &event.nhTpc );
#if RawHit
  tree->Branch( "raw_hitpos_x", &event.raw_hitpos_x );
  tree->Branch( "raw_hitpos_y", &event.raw_hitpos_y );
  tree->Branch( "raw_hitpos_z", &event.raw_hitpos_z );
  tree->Branch( "raw_de", &event.raw_de );
  tree->Branch( "raw_padid", &event.raw_padid );
  tree->Branch( "raw_layer", &event.raw_layer );
  tree->Branch( "raw_row", &event.raw_row );
#endif

  tree->Branch( "nclTpc", &event.nclTpc );
  tree->Branch( "remain_nclTpc", &event.remain_nclTpc );
#if RawCluster
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
  tree->Branch( "cluster_houghflag", &event.cluster_houghflag );
#endif

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "ntKuramaCandidate", &event.ntKuramaCandidate );
  tree->Branch( "isKuramaCandidate", &event.isKuramaCandidate );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "nhtrackEff", &event.nhtrackEff );
  tree->Branch( "trackid", &event.trackid );
  tree->Branch( "isBeam", &event.isBeam );
  tree->Branch( "isXi", &event.isXi );
  tree->Branch( "isKurama", &event.isKurama );
  tree->Branch( "isK18", &event.isK18 );
  tree->Branch( "isAccidental", &event.isAccidental );
  tree->Branch( "isMultiloop", &event.isMultiloop );
  tree->Branch( "flag", &event.flag );
  tree->Branch( "fittime", &event.fittime );
  tree->Branch( "searchtime", &event.searchtime );
  tree->Branch( "chisqr", &event.chisqr );
  tree->Branch( "pval", &event.pval );
  tree->Branch( "distTgt", &event.distTgt );
  tree->Branch( "niteration", &event.niteration );
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
#endif
  tree->Branch( "dz_factor", &event.dz_factor );
  tree->Branch( "charge", &event.charge );
  tree->Branch( "path", &event.path );

  tree->Branch( "pid", &event.pid );
  tree->Branch( "hitlayer", &event.hitlayer );
  tree->Branch( "hitpos_x", &event.hitpos_x );
  tree->Branch( "hitpos_y", &event.hitpos_y );
  tree->Branch( "hitpos_z", &event.hitpos_z );
  tree->Branch( "calpos_x", &event.calpos_x );
  tree->Branch( "calpos_y", &event.calpos_y );
  tree->Branch( "calpos_z", &event.calpos_z );
  tree->Branch( "mom_x", &event.mom_x );
  tree->Branch( "mom_y", &event.mom_y );
  tree->Branch( "mom_z", &event.mom_z );
  tree->Branch( "residual", &event.residual );
  tree->Branch( "residual_t", &event.residual_t );
  tree->Branch( "residual_x", &event.residual_x );
  tree->Branch( "residual_y", &event.residual_y );
  tree->Branch( "residual_z", &event.residual_z );
  tree->Branch( "resolution_x", &event.resolution_x);
  tree->Branch( "resolution_y", &event.resolution_y);
  tree->Branch( "resolution_z", &event.resolution_z);
  tree->Branch( "pull", &event.pull);
  tree->Branch( "helix_t", &event.helix_t );
  tree->Branch( "alpha", &event.alpha);
  tree->Branch( "pathhit", &event.pathhit);
  tree->Branch( "houghflag", &event.houghflag );
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_size", &event.track_cluster_size);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch( "track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch( "track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch( "track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch( "track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch( "track_cluster_row_center", &event.track_cluster_row_center);
  tree->Branch( "exresidual_t", &event.exresidual_t );
  tree->Branch( "exresidual_x", &event.exresidual_x );
  tree->Branch( "exresidual_y", &event.exresidual_y );
  tree->Branch( "exresidual_z", &event.exresidual_z );
  tree->Branch( "intrinsic_residual_t", &event.intrinsic_residual_t );
  tree->Branch( "intrinsic_residual_x", &event.intrinsic_residual_x );
  tree->Branch( "intrinsic_residual_y", &event.intrinsic_residual_y );
  tree->Branch( "intrinsic_residual_z", &event.intrinsic_residual_z );

  tree->Branch( "chargeIndistinguishable", &event.chargeIndistinguishable );
  tree->Branch( "chisqr_inverted", &event.chisqr_inverted );
  tree->Branch( "pval_inverted", &event.pval_inverted );
  tree->Branch( "helix_cx_inverted", &event.helix_cx_inverted );
  tree->Branch( "helix_cy_inverted", &event.helix_cy_inverted );
  tree->Branch( "helix_z0_inverted", &event.helix_z0_inverted );
  tree->Branch( "helix_r_inverted", &event.helix_r_inverted );
  tree->Branch( "helix_dz_inverted", &event.helix_dz_inverted );
  tree->Branch( "mom0_x_inverted", &event.mom0_x_inverted );
  tree->Branch( "mom0_y_inverted", &event.mom0_y_inverted );
  tree->Branch( "mom0_z_inverted", &event.mom0_z_inverted );
  tree->Branch( "mom0_inverted", &event.mom0_inverted );
  tree->Branch( "pid_inverted", &event.pid_inverted );

  tree->Branch( "vpntTpc", &event.vpntTpc );
  tree->Branch( "vpnhtrack", &event.vpnhtrack );
  tree->Branch( "vptrackid", &event.vptrackid );
  tree->Branch( "vpisKurama", &event.vpisKurama );
  tree->Branch( "vpisK18", &event.vpisK18 );
  tree->Branch( "vphelix_cx", &event.vphelix_cx );
  tree->Branch( "vphelix_cy", &event.vphelix_cy );
  tree->Branch( "vphelix_z0", &event.vphelix_z0 );
  tree->Branch( "vphelix_r", &event.vphelix_r );
  tree->Branch( "vphelix_dz", &event.vphelix_dz );
  tree->Branch( "vpmom0", &event.vpmom0 );
  tree->Branch( "vpcharge", &event.vpcharge );
  tree->Branch( "vphelix_t", &event.vphelix_t );
  tree->Branch( "vppos_x", &event.vppos_x);
  tree->Branch( "vppos_y", &event.vppos_y);
  tree->Branch( "vppos_z", &event.vppos_z);
  tree->Branch( "residual_vppos_x", &event.residual_vppos_x);
  tree->Branch( "residual_vppos_y", &event.residual_vppos_y);
  tree->Branch( "residual_vppos_z", &event.residual_vppos_z);

#if TrackSearchFailed
  tree->Branch( "failed_ntTpc", &event.failed_ntTpc );
  tree->Branch( "failed_nhtrack", &event.failed_nhtrack );
  tree->Branch( "failed_flag", &event.failed_flag );
  tree->Branch( "failed_trackid", &event.failed_trackid );
  tree->Branch( "failed_isBeam", &event.failed_isBeam );
  tree->Branch( "failed_isKurama", &event.failed_isKurama );
  tree->Branch( "failed_isK18", &event.failed_isK18 );
  tree->Branch( "failed_nclbeforetgt", &event.failed_nclbeforetgt );
  tree->Branch( "failed_isAccidental", &event.failed_isAccidental );
  tree->Branch( "failed_fittime", &event.failed_fittime );
  tree->Branch( "failed_searchtime", &event.failed_searchtime );
  tree->Branch( "failed_niteration", &event.failed_niteration );

  tree->Branch( "failed_helix_cx", &event.failed_helix_cx );
  tree->Branch( "failed_helix_cy", &event.failed_helix_cy );
  tree->Branch( "failed_helix_z0", &event.failed_helix_z0 );
  tree->Branch( "failed_helix_r", &event.failed_helix_r );
  tree->Branch( "failed_helix_dz", &event.failed_helix_dz );
  tree->Branch( "failed_mom0", &event.failed_mom0 );
  tree->Branch( "failed_charge", &event.failed_charge );

  tree->Branch( "failed_hitlayer", &event.failed_hitlayer );
  tree->Branch( "failed_hitpos_x", &event.failed_hitpos_x );
  tree->Branch( "failed_hitpos_y", &event.failed_hitpos_y );
  tree->Branch( "failed_hitpos_z", &event.failed_hitpos_z );
  tree->Branch( "failed_calpos_x", &event.failed_calpos_x );
  tree->Branch( "failed_calpos_y", &event.failed_calpos_y );
  tree->Branch( "failed_calpos_z", &event.failed_calpos_z );
  tree->Branch( "failed_helix_t", &event.failed_helix_t );
  tree->Branch( "failed_residual", &event.failed_residual );
  tree->Branch( "failed_residual_x", &event.failed_residual_x );
  tree->Branch( "failed_residual_y", &event.failed_residual_y );
  tree->Branch( "failed_residual_z", &event.failed_residual_z );
  tree->Branch( "failed_track_cluster_de", &event.failed_track_cluster_de);
  tree->Branch( "failed_track_cluster_size", &event.failed_track_cluster_size);
  tree->Branch( "failed_track_cluster_mrow", &event.failed_track_cluster_mrow);
#endif

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

  tree->Branch( "nvtxTpcClustered", &event.nvtxTpcClustered );
  tree->Branch( "clusteredVtx_x", &event.Clusteredvtx_x );
  tree->Branch( "clusteredVtx_y", &event.Clusteredvtx_y );
  tree->Branch( "clusteredVtx_z", &event.Clusteredvtx_z );
  tree->Branch( "clusteredVtxid", &event.Clusteredvtxid );

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

  TTreeCont[kK18HSTracking]->SetBranchStatus("*", 0);
  TTreeCont[kK18HSTracking]->SetBranchStatus("ntK18",       1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("nhK18"  ,     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("chisqrK18",   1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("pHS"    ,     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("p_3rd",       1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("delta_3rd",   1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("xout",        1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("yout",        1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("uout",        1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("vout",        1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("layerK18",    1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("wireK18",     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("localhitposK18", 1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("wposK18",     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("xtgtHS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("ytgtHS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("ztgtHS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("utgtHS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("vtgtHS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("initmomHS",   1);

  TTreeCont[kK18HSTracking]->SetBranchStatus("xvp1HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("yvp1HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("zvp1HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("uvp1HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("vvp1HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("xvp2HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("yvp2HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("zvp2HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("uvp2HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("vvp2HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("xvp3HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("yvp3HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("zvp3HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("uvp3HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("vvp3HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("xvp4HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("yvp4HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("zvp4HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("uvp4HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("vvp4HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("xbh2HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("ybh2HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("zbh2HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("ubh2HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("vbh2HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("pbh2HS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("xgasvesselHS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("ygasvesselHS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("zgasvesselHS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("ugasvesselHS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("vgasvesselHS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("pgasvesselHS",      1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("xhtofHS",     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("yhtofHS",     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("zhtofHS",     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("uhtofHS",     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("vhtofHS",     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("xbcHS"  ,     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("ybcHS"  ,     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("zbcHS"  ,     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("ubcHS"  ,     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("vbcHS"  ,     1);
  TTreeCont[kK18HSTracking]->SetBranchStatus("initmomHS",   1);

  TTreeCont[kK18HSTracking]->SetBranchAddress("ntK18",     &src.ntK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("nhK18",     src.nhK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("chisqrK18", src.chisqrK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("pHS",       src.pHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("p_3rd",     src.p_3rd);
  TTreeCont[kK18HSTracking]->SetBranchAddress("delta_3rd", src.delta_3rd);
  TTreeCont[kK18HSTracking]->SetBranchAddress("xout" ,    src.xoutK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("yout" ,    src.youtK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("uout" ,    src.uoutK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("vout" ,    src.voutK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("layerK18", src.layerK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("wireK18",  src.wireK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("localhitposK18", src.localhitposK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("wposK18",  src.wposK18);
  TTreeCont[kK18HSTracking]->SetBranchAddress("xtgtHS",   src.xtgtHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("ytgtHS",   src.ytgtHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("ztgtHS",   src.ztgtHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("utgtHS",   src.utgtHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("vtgtHS",   src.vtgtHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("initmomHS",  src.initmomHS);

  TTreeCont[kK18HSTracking]->SetBranchAddress("xvp1HS",   src.xvp1HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("yvp1HS",   src.yvp1HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("zvp1HS",   src.zvp1HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("uvp1HS",   src.uvp1HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("vvp1HS",   src.vvp1HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("xvp2HS",   src.xvp2HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("yvp2HS",   src.yvp2HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("zvp2HS",   src.zvp2HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("uvp2HS",   src.uvp2HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("vvp2HS",   src.vvp2HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("xvp3HS",   src.xvp3HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("yvp3HS",   src.yvp3HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("zvp3HS",   src.zvp3HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("uvp3HS",   src.uvp3HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("vvp3HS",   src.vvp3HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("xvp4HS",   src.xvp4HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("yvp4HS",   src.yvp4HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("zvp4HS",   src.zvp4HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("uvp4HS",   src.uvp4HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("vvp4HS",   src.vvp4HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("xbh2HS",   src.xbh2HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("ybh2HS",   src.ybh2HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("zbh2HS",   src.zbh2HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("ubh2HS",   src.ubh2HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("vbh2HS",   src.vbh2HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("pbh2HS",   src.pbh2HS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("xgasvesselHS",   src.xgasvesselHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("ygasvesselHS",   src.ygasvesselHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("zgasvesselHS",   src.zgasvesselHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("ugasvesselHS",   src.ugasvesselHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("vgasvesselHS",   src.vgasvesselHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("pgasvesselHS",   src.pgasvesselHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("xhtofHS",  src.xhtofHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("yhtofHS",  src.yhtofHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("zhtofHS",  src.zhtofHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("uhtofHS",  src.uhtofHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("vhtofHS",  src.vhtofHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("xbcHS"  ,  src.xbcHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("ybcHS"  ,  src.ybcHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("zbcHS"  ,  src.zbcHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("ubcHS"  ,  src.ubcHS);
  TTreeCont[kK18HSTracking]->SetBranchAddress("vbcHS"  ,  src.vbcHS);

  TTreeCont[kKuramaTracking]->SetBranchStatus("*", 0);
  TTreeCont[kKuramaTracking]->SetBranchStatus("ntKurama",       1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("chisqrKurama",   1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("pKurama",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("qKurama",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("m2",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("xtgtKurama",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("ytgtKurama",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("utgtKurama",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("vtgtKurama",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("xtofKurama",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("ytofKurama",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("utofKurama",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("vtofKurama",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("tofsegKurama",   1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("thetaKurama",    1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("phiKurama",      1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("path",           1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("stof",           1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("vpxtpc",         1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("vpytpc",         1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("vpztpc",         1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("vpxhtof",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("vpyhtof",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("vpzhtof",        1);

  TTreeCont[kKuramaTracking]->SetBranchStatus("nh",          1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("xout",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("yout",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("zout",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("pxout",       1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("pyout",       1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("pzout",       1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("xtof",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("ytof",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("ztof",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("pxtof",       1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("pytof",       1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("pztof",       1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("layer",       1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("wire",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("localhitpos", 1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("wire",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("wpos",        1);

  TTreeCont[kKuramaTracking]->SetBranchAddress("ntKurama",       &src.ntKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("chisqrKurama",   src.chisqrKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("pKurama",        src.pKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("qKurama",        src.qKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("m2",             src.m2);
  TTreeCont[kKuramaTracking]->SetBranchAddress("xtgtKurama",     src.xtgtKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("ytgtKurama",     src.ytgtKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("utgtKurama",     src.utgtKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("vtgtKurama",     src.vtgtKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("xtofKurama",     src.xtofKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("ytofKurama",     src.ytofKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("utofKurama",     src.utofKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("vtofKurama",     src.vtofKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("tofsegKurama",   src.tofsegKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("thetaKurama",    src.thetaKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("phiKurama",      src.phiKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("path",           src.path);
  TTreeCont[kKuramaTracking]->SetBranchAddress("stof",           src.stof);
  TTreeCont[kKuramaTracking]->SetBranchAddress("vpxtpc",         src.vpxtpc);
  TTreeCont[kKuramaTracking]->SetBranchAddress("vpytpc",         src.vpytpc);
  TTreeCont[kKuramaTracking]->SetBranchAddress("vpztpc",         src.vpztpc);
  TTreeCont[kKuramaTracking]->SetBranchAddress("vpxhtof",        src.vpxhtof);
  TTreeCont[kKuramaTracking]->SetBranchAddress("vpyhtof",        src.vpyhtof);
  TTreeCont[kKuramaTracking]->SetBranchAddress("vpzhtof",        src.vpzhtof);

  TTreeCont[kKuramaTracking]->SetBranchAddress("nh",    src.nh);
  TTreeCont[kKuramaTracking]->SetBranchAddress("xout",  src.xout);
  TTreeCont[kKuramaTracking]->SetBranchAddress("yout",  src.yout);
  TTreeCont[kKuramaTracking]->SetBranchAddress("zout",  src.zout);
  TTreeCont[kKuramaTracking]->SetBranchAddress("pxout", src.pxout);
  TTreeCont[kKuramaTracking]->SetBranchAddress("pyout", src.pyout);
  TTreeCont[kKuramaTracking]->SetBranchAddress("pzout", src.pzout);
  TTreeCont[kKuramaTracking]->SetBranchAddress("xtof",  src.xtof);
  TTreeCont[kKuramaTracking]->SetBranchAddress("ytof",  src.ytof);
  TTreeCont[kKuramaTracking]->SetBranchAddress("ztof",  src.ztof);
  TTreeCont[kKuramaTracking]->SetBranchAddress("pxtof", src.pxtof);
  TTreeCont[kKuramaTracking]->SetBranchAddress("pytof", src.pytof);
  TTreeCont[kKuramaTracking]->SetBranchAddress("pztof", src.pztof);
  TTreeCont[kKuramaTracking]->SetBranchAddress("layer", src.layer);
  TTreeCont[kKuramaTracking]->SetBranchAddress("wire",  src.wire);
  TTreeCont[kKuramaTracking]->SetBranchAddress("localhitpos", src.localhitpos);
  TTreeCont[kKuramaTracking]->SetBranchAddress("wpos",  src.wpos);

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

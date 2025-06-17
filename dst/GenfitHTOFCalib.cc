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
#include "TPCPadHelper.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCVertex.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#include "HypTPCFitter.hh"
#include "HypTPCTask.hh"

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
const auto& zHSCenter = gGeom.LocalZ("HS");
const Double_t truncatedMean = 0.8; //80%

//For GenFit Setting
const Bool_t Const_field = false; //Must be false for linear tracking
const Int_t verbosity = 0;//3;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kHTOFCaib, kHodoscope1, kHodoscope2, kOutFile, nArgc
};
std::vector<TString> ArgName =
  { "[Process]", "[ConfFile]", "[HTOFCaib]"/* or DstE42*/, "[Hodoscope]", "[Hodoscope]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc", "hodo", "tree", "" };
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
  std::vector<Int_t> cluster_houghflag;

  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Int_t> isBeam;
  std::vector<Int_t> isKurama;
  std::vector<Int_t> isK18;
  std::vector<Int_t> isAccidental;
  std::vector<Double_t> chisqr;
  std::vector<Double_t> helix_cx;
  std::vector<Double_t> helix_cy;
  std::vector<Double_t> helix_z0;
  std::vector<Double_t> helix_r;
  std::vector<Double_t> helix_dz;
  std::vector<Double_t> dE;
  std::vector<Double_t> dEdx;
  std::vector<Double_t> mom0;//Helix momentum at Y = 0
  std::vector<Int_t> charge;//Helix charge
  std::vector<Double_t> path;//Helix path
  std::vector<Int_t> pid;
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

  std::vector<Int_t> isgoodTPCKurama;
  std::vector<Int_t> insideTPC;
  std::vector<Double_t> pTPCKurama;
  std::vector<Double_t> qTPCKurama;
  std::vector<Double_t> m2TPCKurama;
  std::vector<Double_t> xsTPC;
  std::vector<Double_t> ysTPC;
  std::vector<Double_t> usTPC;
  std::vector<Double_t> vsTPC;

  std::vector<Double_t> pK18;
  std::vector<Double_t> xbTPC;
  std::vector<Double_t> ybTPC;
  std::vector<Double_t> ubTPC;
  std::vector<Double_t> vbTPC;

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

  Int_t nhHtof;
  std::vector<Double_t> HtofSeg;
  std::vector<Double_t> tHtof;
  std::vector<Double_t> dtHtof;
  std::vector<Double_t> deHtof;
  std::vector<Double_t> posHtof;

  std::vector<Double_t> utimeHtof;
  std::vector<Double_t> dtimeHtof;
  std::vector<Double_t> uctimeHtof;
  std::vector<Double_t> dctimeHtof;
  std::vector<Double_t> udeHtof;
  std::vector<Double_t> ddeHtof;

  Int_t GFstatus;
  Int_t GFntTpc;
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

  std::vector<Double_t> GFmom_p;
  std::vector<Double_t> GFtracklen_p;
  std::vector<Double_t> GFtrack2vtxdist_p;
  std::vector<Double_t> GFcalctof_p;
  std::vector<Double_t> GFsegHtof_p;
  std::vector<Double_t> GFtofHtof_p;
  std::vector<Double_t> GFtdiffHtof_p;
  std::vector<Double_t> GFposHtof_p;
  std::vector<Double_t> GFposx_p;
  std::vector<Double_t> GFposy_p;
  std::vector<Double_t> GFposz_p;
  std::vector<Double_t> GFinvbeta_p;
  std::vector<Double_t> GFm2_p;

  std::vector<Double_t> GFmom_pi;
  std::vector<Double_t> GFtracklen_pi;
  std::vector<Double_t> GFtrack2vtxdist_pi;
  std::vector<Double_t> GFcalctof_pi;
  std::vector<Double_t> GFsegHtof_pi;
  std::vector<Double_t> GFtofHtof_pi;
  std::vector<Double_t> GFtdiffHtof_pi;
  std::vector<Double_t> GFposHtof_pi;
  std::vector<Double_t> GFposx_pi;
  std::vector<Double_t> GFposy_pi;
  std::vector<Double_t> GFposz_pi;
  std::vector<Double_t> GFinvbeta_pi;
  std::vector<Double_t> GFm2_pi;

  void clear( void )
  {
    runnum = 0;
    evnum = 0;
    status = 0;
    trigpat.clear();
    trigflag.clear();

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
    cluster_houghflag.clear();

    ntTpc = 0;
    nhtrack.clear();
    isBeam.clear();
    isKurama.clear();
    isK18.clear();
    isAccidental.clear();
    chisqr.clear();
    helix_cx.clear();
    helix_cy.clear();
    helix_z0.clear();
    helix_r.clear();
    helix_dz.clear();
    dE.clear();
    dEdx.clear();

    mom0.clear();
    charge.clear();
    path.clear();
    pid.clear();
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

    isgoodTPCKurama.clear();
    insideTPC.clear();
    pTPCKurama.clear();
    qTPCKurama.clear();
    m2TPCKurama.clear();
    xsTPC.clear();
    ysTPC.clear();
    usTPC.clear();
    vsTPC.clear();

    pK18.clear();
    xbTPC.clear();
    ybTPC.clear();
    ubTPC.clear();
    vbTPC.clear();

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

    nhHtof = 0;
    HtofSeg.clear();
    tHtof.clear();
    dtHtof.clear();
    deHtof.clear();
    posHtof.clear();

    utimeHtof.clear();
    dtimeHtof.clear();
    uctimeHtof.clear();
    dctimeHtof.clear();
    udeHtof.clear();
    ddeHtof.clear();

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

    GFmom_p.clear();
    GFtracklen_p.clear();
    GFtrack2vtxdist_p.clear();
    GFcalctof_p.clear();
    GFsegHtof_p.clear();
    GFtofHtof_p.clear();
    GFtdiffHtof_p.clear();
    GFposHtof_p.clear();
    GFposx_p.clear();
    GFposy_p.clear();
    GFposz_p.clear();
    GFinvbeta_p.clear();
    GFm2_p.clear();

    GFmom_pi.clear();
    GFtracklen_pi.clear();
    GFtrack2vtxdist_pi.clear();
    GFcalctof_pi.clear();
    GFsegHtof_pi.clear();
    GFtofHtof_pi.clear();
    GFtdiffHtof_pi.clear();
    GFposHtof_pi.clear();
    GFposx_pi.clear();
    GFposy_pi.clear();
    GFposz_pi.clear();
    GFinvbeta_pi.clear();
    GFm2_pi.clear();
  }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* runnum;
  TTreeReaderValue<Int_t>* evnum;
  TTreeReaderValue<std::vector<Int_t>>* trigpat;
  TTreeReaderValue<std::vector<Int_t>>* trigflag;

  TTreeReaderValue<Int_t>* nclTpc; // Number of clusters
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
  TTreeReaderValue<std::vector<Int_t>>* isBeam;
  TTreeReaderValue<std::vector<Int_t>>* isKurama;
  TTreeReaderValue<std::vector<Int_t>>* isK18;
  TTreeReaderValue<std::vector<Int_t>>* isAccidental;
  TTreeReaderValue<std::vector<Double_t>>* chisqr;
  TTreeReaderValue<std::vector<Double_t>>* helix_cx;
  TTreeReaderValue<std::vector<Double_t>>* helix_cy;
  TTreeReaderValue<std::vector<Double_t>>* helix_z0;
  TTreeReaderValue<std::vector<Double_t>>* helix_r;
  TTreeReaderValue<std::vector<Double_t>>* helix_dz;
  TTreeReaderValue<std::vector<Double_t>>* dE;
  TTreeReaderValue<std::vector<Double_t>>* dEdx; //reference dedx
  TTreeReaderValue<std::vector<Double_t>>* mom0;//Helix momentum at Y = 0
  TTreeReaderValue<std::vector<Int_t>>* charge;//Helix charge
  TTreeReaderValue<std::vector<Double_t>>* path;//Helix path
  TTreeReaderValue<std::vector<Int_t>>* pid;

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
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* alpha;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* pathhit;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_size;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_mrow;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_x_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_y_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_z_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_row_center;

  TTreeReaderValue<std::vector<Int_t>>* isgoodTPCKurama;
  TTreeReaderValue<std::vector<Int_t>>* insideTPC;
  TTreeReaderValue<std::vector<Double_t>>* pTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* qTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* m2TPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* xsTPC;
  TTreeReaderValue<std::vector<Double_t>>* ysTPC;
  TTreeReaderValue<std::vector<Double_t>>* usTPC;
  TTreeReaderValue<std::vector<Double_t>>* vsTPC;

  TTreeReaderValue<std::vector<Double_t>>* pK18;
  TTreeReaderValue<std::vector<Double_t>>* xbTPC;
  TTreeReaderValue<std::vector<Double_t>>* ybTPC;
  TTreeReaderValue<std::vector<Double_t>>* ubTPC;
  TTreeReaderValue<std::vector<Double_t>>* vbTPC;

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

  Double_t Time0;
  Double_t CTime0;

  Int_t    nhHtof;
  Int_t    csHtof[NumOfSegHTOF*MaxDepth];
  Double_t HtofSeg[NumOfSegHTOF*MaxDepth];
  Double_t tHtof[NumOfSegHTOF*MaxDepth];
  Double_t dtHtof[NumOfSegHTOF*MaxDepth];
  Double_t deHtof[NumOfSegHTOF*MaxDepth];
  Double_t posHtof[NumOfSegHTOF*MaxDepth];

  Double_t htofmt[NumOfSegHTOF][MaxDepth];
  Double_t htofde[NumOfSegHTOF];
  Double_t htofutime[NumOfSegHTOF][MaxDepth];
  Double_t htofuctime[NumOfSegHTOF][MaxDepth];
  Double_t htofdtime[NumOfSegHTOF][MaxDepth];
  Double_t htofdctime[NumOfSegHTOF][MaxDepth];
  Double_t htofhitpos[NumOfSegHTOF][MaxDepth];
  Double_t htofude[NumOfSegHTOF];
  Double_t htofdde[NumOfSegHTOF];

};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
  enum eDetHid {
    PadHid    = 100000,
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
  HypTPCFitter* fitter = new HypTPCFitter(tpcGeo.Data(),Const_field);
  //Initiallize the genfit track container
  HypTPCTask& GFtracks = HypTPCTask::GetInstance();
  GFtracks.SetVerbosity(verbosity);
  std::cout<<"GenFit verbosity = "<<"-1: Silent, 0: Minimum, 1: Errors only, 2: Errors and Warnings, 3: Verbose mode, long term debugging(default)"<<std::endl;
  std::cout<<"Current verbosity = "<<GFtracks.GetVerbosity()<<std::endl;

#if 0
  GFtracks.DebugMode();
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

  Double_t vtx_scan_range = gUser.GetParameter("VertexScanRange");

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

  event.nclTpc = **src.nclTpc;
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

  int ntTpc = **src.ntTpc;
  event.ntTpc = ntTpc;
  event.nhtrack = **src.nhtrack;
  event.isBeam = **src.isBeam;
  event.isKurama = **src.isKurama;
  event.isK18 = **src.isK18;
  event.isAccidental = **src.isAccidental;
  event.chisqr = **src.chisqr;
  event.helix_cx = **src.helix_cx;
  event.helix_cy = **src.helix_cy;
  event.helix_z0 = **src.helix_z0;
  event.helix_r = **src.helix_r;
  event.helix_dz = **src.helix_dz;

  event.dE = **src.dE;
  event.dEdx = **src.dEdx;
  event.mom0 = **src.mom0;
  event.charge = **src.charge;
  event.path = **src.path;
  event.pid = **src.pid;

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
  event.alpha = **src.alpha;
  event.pathhit = **src.pathhit;
  event.track_cluster_de = **src.track_cluster_de;
  event.track_cluster_size = **src.track_cluster_size;
  event.track_cluster_mrow = **src.track_cluster_mrow;
  event.track_cluster_de_center = **src.track_cluster_de_center;
  event.track_cluster_x_center = **src.track_cluster_x_center;
  event.track_cluster_y_center = **src.track_cluster_y_center;
  event.track_cluster_z_center = **src.track_cluster_z_center;
  event.track_cluster_row_center = **src.track_cluster_row_center;

  event.isgoodTPCKurama = **src.isgoodTPCKurama;
  event.insideTPC = **src.insideTPC;
  event.pTPCKurama = **src.pTPCKurama;
  event.qTPCKurama = **src.qTPCKurama;
  event.m2TPCKurama = **src.m2TPCKurama;
  event.xsTPC = **src.xsTPC;
  event.ysTPC = **src.ysTPC;
  event.usTPC = **src.usTPC;
  event.vsTPC = **src.vsTPC;

  event.pK18 = **src.pK18;
  event.xbTPC = **src.xbTPC;
  event.ybTPC = **src.ybTPC;
  event.ubTPC = **src.ubTPC;
  event.vbTPC = **src.vbTPC;

  event.isElectron.resize(ntTpc);
  event.nsigma_triton.resize(ntTpc);
  event.nsigma_deutron.resize(ntTpc);
  event.nsigma_proton.resize(ntTpc);
  event.nsigma_kaon.resize(ntTpc);
  event.nsigma_pion.resize(ntTpc);
  event.nsigma_electron.resize(ntTpc);
  for(int it=0; it<ntTpc; ++it){
    event.isElectron[it] = Kinematics::HypTPCdEdxElectron(event.dEdx[it], event.mom0[it]);
    event.nsigma_triton[it] = Kinematics::HypTPCdEdxNsigmaTriton(event.dEdx[it], event.mom0[it]);
    event.nsigma_deutron[it] = Kinematics::HypTPCdEdxNsigmaDeutron(event.dEdx[it], event.mom0[it]);
    event.nsigma_proton[it] = Kinematics::HypTPCdEdxNsigmaProton(event.dEdx[it], event.mom0[it]);
    event.nsigma_kaon[it]  = Kinematics::HypTPCdEdxNsigmaKaon(event.dEdx[it], event.mom0[it]);
    event.nsigma_pion[it] = Kinematics::HypTPCdEdxNsigmaPion(event.dEdx[it], event.mom0[it]);
    event.nsigma_electron[it] = Kinematics::HypTPCdEdxNsigmaElectron(event.dEdx[it], event.mom0[it]);
  }

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

  event.nhHtof = src.nhHtof;
  for(Int_t i=0; i<event.nhHtof; i++){
    event.HtofSeg.push_back(src.HtofSeg[i]);
    event.tHtof.push_back(src.tHtof[i]);
    event.dtHtof.push_back(src.dtHtof[i]);
    event.deHtof.push_back(src.deHtof[i]);
    event.posHtof.push_back(src.posHtof[i]);

    double utime = qnan; double dtime = qnan;
    double uctime = qnan; double dctime = qnan;
    double ude = qnan; double dde = qnan;

    int j = src.HtofSeg[i] - 1;
    for(Int_t m=0; m<MaxDepth; ++m){
      if(TMath::IsNaN(src.htofutime[j][m])) continue;
      double cmeantime = 0.5*(src.htofuctime[j][m] + src.htofdctime[j][m]);
      if(TMath::Abs(src.deHtof[i] - src.htofde[j])<0.001 &&
	 TMath::Abs(src.tHtof[i] - cmeantime)<0.001){

	utime = src.htofutime[j][m];
	dtime = src.htofdtime[j][m];
	uctime = src.htofuctime[j][m];
	dctime = src.htofdctime[j][m];
	ude = src.htofude[j];
	dde = src.htofdde[j];
      }
    }

    event.utimeHtof.push_back(utime);
    event.dtimeHtof.push_back(dtime);
    event.uctimeHtof.push_back(uctime);
    event.dctimeHtof.push_back(dctime);
    event.udeHtof.push_back(ude);
    event.ddeHtof.push_back(dde);
  }

  HF1( 1, event.status++ );

  HF1( 10, ntTpc );
  if( event.ntTpc == 0 )
    return true;

  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCTracks(**src.ntTpc, **src.isK18, **src.isKurama,
			 **src.charge, **src.nhtrack, **src.helix_cx,
			 **src.helix_cy, **src.helix_z0, **src.helix_r,
			 **src.helix_dz, **src.hitlayer, **src.track_cluster_mrow,
			 **src.helix_t, **src.track_cluster_de, **src.resolution_x,
			 **src.resolution_y, **src.resolution_z, **src.hitpos_x,
			 **src.hitpos_y, **src.hitpos_z);

  HypTPCTask& GFtracks = HypTPCTask::GetInstance();
  //GFtracks.Init();

  HF1( 1, event.status++ );
  for(int it=0; it<event.ntTpc; ++it){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix(it);
    if( !tp ) continue;
    std::vector<Int_t> pdgcode;
    Kinematics::HypTPCPID_PDGCode(event.charge[it], event.pid[it], pdgcode);
    GFtracks.AddHelixTrack(pdgcode, tp);
  }
  GFtracks.FitTracks();

  HF1( 2, event.GFstatus++ );

  int GFntTpc = GFtracks.GetNTrack();
  if(GFntTpc!=event.ntTpc){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    return true;
  }
  HF1( 2, event.GFstatus++ );

  HF1( genfitHid, GFntTpc);
  event.GFntTpc = GFntTpc;
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
  event.GFsegHtof.resize(GFntTpc);
  event.GFtofHtof.resize(GFntTpc);
  event.GFtdiffHtof.resize(GFntTpc);
  event.GFposHtof.resize(GFntTpc);
  event.GFposx.resize(GFntTpc);
  event.GFposy.resize(GFntTpc);
  event.GFposz.resize(GFntTpc);
  event.GFinvbeta.resize(GFntTpc);
  event.GFm2.resize(GFntTpc);
  event.nsigma_tritonHtof.resize(ntTpc);
  event.nsigma_deutronHtof.resize(ntTpc);
  event.nsigma_protonHtof.resize(ntTpc);
  event.nsigma_kaonHtof.resize(ntTpc);
  event.nsigma_pionHtof.resize(ntTpc);
  event.nsigma_electronHtof.resize(ntTpc);

  event.GFmom_p.resize(GFntTpc);
  event.GFtracklen_p.resize(GFntTpc);
  event.GFtrack2vtxdist_p.resize(GFntTpc);
  event.GFcalctof_p.resize(GFntTpc);
  event.GFsegHtof_p.resize(GFntTpc);
  event.GFtofHtof_p.resize(GFntTpc);
  event.GFtdiffHtof_p.resize(GFntTpc);
  event.GFposHtof_p.resize(GFntTpc);
  event.GFposx_p.resize(GFntTpc);
  event.GFposy_p.resize(GFntTpc);
  event.GFposz_p.resize(GFntTpc);
  event.GFinvbeta_p.resize(GFntTpc);
  event.GFm2_p.resize(GFntTpc);

  event.GFmom_pi.resize(GFntTpc);
  event.GFtracklen_pi.resize(GFntTpc);
  event.GFtrack2vtxdist_pi.resize(GFntTpc);
  event.GFcalctof_pi.resize(GFntTpc);
  event.GFsegHtof_pi.resize(GFntTpc);
  event.GFtofHtof_pi.resize(GFntTpc);
  event.GFtdiffHtof_pi.resize(GFntTpc);
  event.GFposHtof_pi.resize(GFntTpc);
  event.GFposx_pi.resize(GFntTpc);
  event.GFposy_pi.resize(GFntTpc);
  event.GFposz_pi.resize(GFntTpc);
  event.GFinvbeta_pi.resize(GFntTpc);
  event.GFm2_pi.resize(GFntTpc);

  Int_t ntrack_intarget = 0;
  Double_t x0[100] = {0};
  Double_t y0[100] = {0};
  Double_t u0[100] = {0};
  Double_t v0[100] = {0};
  for( Int_t igf=0; igf<GFntTpc; ++igf ){
    event.GFfitstatus[igf] = (int)GFtracks.TrackCheck(igf);
    HF1( 3, event.GFfitstatus[igf]);
    if(!GFtracks.TrackCheck(igf)) continue;
    int nh = GFtracks.GetNHits(igf);
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

    event.GFchisqr[igf] = GFtracks.GetChi2NDF(igf);
    event.GFcharge[igf] = GFtracks.GetCharge(igf);
    event.GFtof[igf] = GFtracks.GetTrackTOF(igf, 0, -1);
    event.GFpval[igf] = GFtracks.GetPvalue(igf);
    event.GFnhtrack[igf] = GFtracks.GetNHits(igf);
    event.GFpdgcode[igf] = GFtracks.GetPDGcode(igf);

    HF1( genfitHid+1, event.GFchisqr[igf]);
    HF1( genfitHid+2, event.GFpval[igf]);
    HF1( genfitHid+3, event.GFcharge[igf]);
    HF1( genfitHid+4, event.GFnhtrack[igf]);
    HF1( genfitHid+5, event.GFtracklen[igf]);
    HF1( genfitHid+6, event.GFtof[igf]);
    for( Int_t ihit=0; ihit<nh; ++ihit ){
      TVector3 hit = GFtracks.GetPos(igf, ihit);
      TVector3 mom = GFtracks.GetMom(igf, ihit);
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
      if(ihit==0) HF1( genfitHid+7, event.GFmom[igf][0]);
      HF1( genfitHid+8, event.GFlayer[igf][ihit]);
      HF1( genfitHid+10, event.GFresidual_x[igf][ihit]);
      HF1( genfitHid+11, event.GFresidual_y[igf][ihit]);
      HF1( genfitHid+12, event.GFresidual_z[igf][ihit]);
      HF1( genfitHid+13, event.GFresidual_p[igf][ihit]);
      HF1( genfitHid+14, event.GFresidual_px[igf][ihit]);
      HF1( genfitHid+15, event.GFresidual_py[igf][ihit]);
      HF1( genfitHid+16, event.GFresidual_pz[igf][ihit]);
      HF1( genfitHid+1000*(layer+1), event.GFresidual_x[igf][ihit]);
      HF1( genfitHid+1000*(layer+1)+1, event.GFresidual_y[igf][ihit]);
      HF1( genfitHid+1000*(layer+1)+2, event.GFresidual_z[igf][ihit]);
      HF1( genfitHid+1000*(layer+1)+3, event.GFresidual_p[igf][ihit]);
      HF1( genfitHid+1000*(layer+1)+4, event.GFresidual_px[igf][ihit]);
      HF1( genfitHid+1000*(layer+1)+5, event.GFresidual_py[igf][ihit]);
      HF1( genfitHid+1000*(layer+1)+6, event.GFresidual_pz[igf][ihit]);
    } //ihit

    //Extrapolation
    if(event.isBeam[igf]==1) continue;
    if(event.isK18[igf]==1) continue;
    if(event.isAccidental[igf]==1) continue;
    if(GFtracks.IsInsideTarget(igf)){
      event.GFinside[igf] = 1;
      TVector3 posv; TVector3 momv; double len; double tof;
      if(GFtracks.ExtrapolateToTargetCenter(igf, posv, momv, len, tof)){
	x0[ntrack_intarget] = posv.x();
	y0[ntrack_intarget] = posv.y();
	u0[ntrack_intarget] = momv.x()/momv.z();
	v0[ntrack_intarget] = momv.y()/momv.z();
	ntrack_intarget++;
      }
    }
    else event.GFinside[igf] = 0;
  } //igf

  TVector3 vertex = Kinematics::MultitrackVertex(ntrack_intarget, x0, y0, u0, v0);
  event.GFntTpc_inside = ntrack_intarget;
  event.GFprodvtx_x = vertex.x();
  event.GFprodvtx_y = vertex.y();
  event.GFprodvtx_z = vertex.z();

  for( Int_t igf=0; igf<GFntTpc; ++igf ){
    if((event.pid[igf]&1)==1){
      Int_t repid = 0; //pion

      Int_t hitid_htof; Double_t tof; Double_t len;
      TVector3 pos_htof; Double_t track2tgt_dist;
      Bool_t htofextrapolation =
	GFtracks.TPCHTOFTrackMatching(igf, repid, vertex,
				      event.HtofSeg, event.posHtof,
				      hitid_htof, tof,
				      len, pos_htof, track2tgt_dist);
      if(htofextrapolation){
	event.GFmom_pi[igf] = GFtracks.GetMom(igf, 0, repid).Mag();
	event.GFtracklen_pi[igf] = len;
	event.GFtrack2vtxdist_pi[igf] = track2tgt_dist;
	event.GFcalctof_pi[igf] = tof;
	event.GFposx_pi[igf] = pos_htof.x();
	event.GFposy_pi[igf] = pos_htof.y();
	event.GFposz_pi[igf] = pos_htof.z();
	event.GFsegHtof_pi[igf] = event.HtofSeg[hitid_htof];
	event.GFtofHtof_pi[igf] = event.tHtof[hitid_htof];
	event.GFtdiffHtof_pi[igf] = event.dtHtof[hitid_htof];
	event.GFposHtof_pi[igf] = event.posHtof[hitid_htof];

	Double_t beta = len/event.tHtof[hitid_htof]/MathTools::C();
	event.GFinvbeta_pi[igf] = 1./beta;
	Double_t mass2 = Kinematics::MassSquare(event.GFmom_pi[igf], len, event.tHtof[hitid_htof]);
	event.GFm2_pi[igf] = mass2;
      }
    } //pion

    if(event.charge[igf]==1 && (event.pid[igf]&4)==4){ //proton
      Int_t repid = 0;
      Int_t flag = 1;
      for(Int_t i=0;i<2;i++){
	Int_t temp = flag&event.pid[igf];
	if(temp==flag) repid += 1;
	flag*=2;
      }

      Int_t hitid_htof; Double_t tof; Double_t len;
      TVector3 pos_htof; Double_t track2tgt_dist;
      Bool_t htofextrapolation =
	GFtracks.TPCHTOFTrackMatching(igf, repid, vertex,
				      event.HtofSeg, event.posHtof,
				      hitid_htof, tof,
				      len, pos_htof, track2tgt_dist);
      if(htofextrapolation){
	event.GFmom_p[igf] = GFtracks.GetMom(igf, 0, repid).Mag();
	event.GFtracklen_p[igf] = len;
	event.GFtrack2vtxdist_p[igf] = track2tgt_dist;
	event.GFcalctof_p[igf] = tof;
	event.GFposx_p[igf] = pos_htof.x();
	event.GFposy_p[igf] = pos_htof.y();
	event.GFposz_p[igf] = pos_htof.z();
	event.GFsegHtof_p[igf] = event.HtofSeg[hitid_htof];
	event.GFtofHtof_p[igf] = event.tHtof[hitid_htof];
	event.GFtdiffHtof_p[igf] = event.dtHtof[hitid_htof];
	event.GFposHtof_p[igf] = event.posHtof[hitid_htof];

	Double_t beta = len/event.tHtof[hitid_htof]/MathTools::C();
	event.GFinvbeta_p[igf] = 1./beta;
	Double_t mass2 = Kinematics::MassSquare(event.GFmom_p[igf], len, event.tHtof[hitid_htof]);
	event.GFm2_p[igf] = mass2;
      }
    } //proton

    {
      Int_t repid = -1;
      Int_t hitid_htof; Double_t tof; Double_t len;
      TVector3 pos_htof; Double_t track2tgt_dist;
      Bool_t htofextrapolation =
	GFtracks.TPCHTOFTrackMatching(igf, repid, vertex,
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
	event.GFtdiffHtof[igf] = event.dtHtof[hitid_htof];
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
  } //igf

  HF1( 2, event.GFstatus++);
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
  /*
    HB1( 1, "Status", 21, 0., 21. );
    HB1( 2, "Genfit Status", 20, 0., 20. );
    HB1( 3, "Genfit Fit Status", 2, 0., 2. );
    HB1( 10, "NTrack TPC", 20, 0., 20. );

    HB1(genfitHid, "[GenFit] #Track TPC; #Track; Counts", 20, 0., 20. );
    HB1(genfitHid+1, "[GenFit] Chisqr/ndf; ; Counts", 200, 0, 10 );
    HB1(genfitHid+2, "[GenFit] p-value; p-value; Counts", 100, -0.05, 1.05);
    HB1(genfitHid+3, "[GenFit] Charge;", 6, -3, 3 );
    HB1(genfitHid+4, "[GenFit] #Hits of Track", 50, 0., 50. );
    HB1(genfitHid+5, "[GenFit] Track Length; Length [mm]; Counts [/1 mm]", 500, 0, 500 );
    HB1(genfitHid+6, "[GenFit] Tof; Tof [ns]; Counts [/0.01 ns]", 500, 0, 5 );
    HB1(genfitHid+7, "[GenFit] Reconstructed P; P [GeV/c]; Counts [/0.001 GeV/c]", 1500, 0., 1.5 );
    HB1(genfitHid+8, "[GenFit] LayerID", 33, 0., 33. );

    //Residuals
    HB1(genfitHid+10, "[GenFit] Residual x; Residual x [mm]; Counts [/0.1 mm]", 400, -20., 20.);
    HB1(genfitHid+11, "[GenFit] Residual y; Residual y [mm]; Counts [/0.1 mm]", 400, -20., 20.);
    HB1(genfitHid+12, "[GenFit] Residual z; Residual z [mm]; Counts [/0.1 mm]", 400, -20., 20.);
    HB1(genfitHid+13, "[GenFit] Residual P; Residual P [GeV/c]; Counts [/0.001 GeV/c]", 400, -0.2, 0.2 );
    HB1(genfitHid+14, "[GenFit] Residual Px; Residual Px [GeV/c]; Counts [/0.001 GeV/c]", 400, -0.2, 0.2 );
  HB1(genfitHid+15, "[GenFit] Residual Py; Residual Py [GeV/c]; Counts [/0.001 GeV/c]", 400, -0.2, 0.2 );
  HB1(genfitHid+16, "[GenFit] Residual Pz; Residual Pz [GeV/c]; Counts [/0.001 GeV/c]", 400, -0.2, 0.2 );
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    HB1(genfitHid+(layer+1)*1000, Form("[GenFit] Residual x Layer%d; Residual x [mm]; Counts [/0.1 mm]", layer), 400, -20., 20.);
    HB1(genfitHid+(layer+1)*1000+1, Form("[GenFit] Residual y Layer%d; Residual y [mm]; Counts [/0.1 mm]", layer), 400, -20., 20.);
    HB1(genfitHid+(layer+1)*1000+2, Form("[GenFit] Residual z Layer%d; Residual z [mm]; Counts [/0.1 mm]", layer), 400, -20., 20.);
    HB1(genfitHid+(layer+1)*1000+3, Form("[GenFit] Residual P Layer%d; Residual P [GeV/c]; Counts [/0.001 GeV/c]", layer), 400, -0.2, 0.2 );
    HB1(genfitHid+(layer+1)*1000+4, Form("[GenFit] Residual Px Layer%d; Residual Px [GeV/c]; Counts [/0.001 GeV/c]", layer), 400, -0.2, 0.2 );
    HB1(genfitHid+(layer+1)*1000+5, Form("[GenFit] Residual Py Layer%d; Residual Py [GeV/c]; Counts [/0.001 GeV/c]", layer), 400, -0.2, 0.2 );
    HB1(genfitHid+(layer+1)*1000+6, Form("[GenFit] Residual Pz Layer%d; Residual Pz [GeV/c]; Counts [/0.001 GeV/c]", layer), 400, -0.2, 0.2 );
  }

  //Extrapolation
  HB1(genfitHid+20, "[GenFit] Vertex X (extrapolated to the Target); X [mm]; Counts [/0.1 mm]", 500, -25, 25 );
  HB1(genfitHid+21, "[GenFit] Vertex Y (extrapolated to the Target); Vertex Y [mm]; Counts [/0.1 mm]", 500, -25, 25 );
  HB1(genfitHid+22, "[GenFit] Vertex Z (extrapolated to the Target); Vertex Z [mm]; Counts [/0.1 mm]", 500, -143 -25, -143 + 25 );
  HB1(genfitHid+23, "[GenFit] #Hits (extrapolated to the HTOF); #Hits; Counts", 10, 0, 10 );
  HB1(genfitHid+24, "[GenFit] Track Length (extrapolated to the HTOF); Length [mm]; Counts [/1 mm]", 1000, -500, 500 );
  HB1(genfitHid+25, "[GenFit] Tof (extrapolated to the HTOF); Tof [ns]; Counts [/0.01 ns]", 500, 0, 5 );
  HB1(genfitHid+26, "[GenFit] X (extrapolated to the HTOF); X [mm]; Counts [/0.1 mm]", 10000, -500, 500 );
  HB1(genfitHid+27, "[GenFit] Y (extrapolated to the HTOF); Y [mm]; Counts [/0.1 mm]", 8000, -400, 400 );
  HB1(genfitHid+28, "[GenFit] Z (extrapolated to the HTOF); Z [mm]; Counts [/0.1 mm]", 10000, -500, 500 );
  HB1(genfitHid+29, "[GenFit] HTOF ID; #ID ; Counts", 36, 0, 36 );
*/
  HBTree( "tpc", "tree of GenfitHTOFCalib" );

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

  tree->Branch( "utimeHtof", &event.utimeHtof );
  tree->Branch( "dtimeHtof", &event.dtimeHtof );
  tree->Branch( "uctimeHtof", &event.uctimeHtof );
  tree->Branch( "dctimeHtof", &event.dctimeHtof );
  tree->Branch( "udeHtof", &event.udeHtof );
  tree->Branch( "ddeHtof", &event.ddeHtof );

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
  tree->Branch( "cluster_houghflag", &event.cluster_houghflag );

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "isBeam", &event.isBeam );
  tree->Branch( "isK18", &event.isK18 );
  tree->Branch( "isKurama", &event.isKurama );
  tree->Branch( "isAccidental", &event.isAccidental );
  tree->Branch( "chisqr", &event.chisqr );
  tree->Branch( "helix_cx", &event.helix_cx );
  tree->Branch( "helix_cy", &event.helix_cy );
  tree->Branch( "helix_z0", &event.helix_z0 );
  tree->Branch( "helix_r", &event.helix_r );
  tree->Branch( "helix_dz", &event.helix_dz );
  tree->Branch( "mom0", &event.mom0 );
  tree->Branch( "dE", &event.dE );
  tree->Branch( "dEdx", &event.dEdx );
  tree->Branch( "isElectron", &event.isElectron );
  tree->Branch( "nsigma_triton", &event.nsigma_triton );
  tree->Branch( "nsigma_deutron", &event.nsigma_deutron );
  tree->Branch( "nsigma_proton", &event.nsigma_proton );
  tree->Branch( "nsigma_kaon", &event.nsigma_kaon );
  tree->Branch( "nsigma_pion", &event.nsigma_pion );
  tree->Branch( "nsigma_electron", &event.nsigma_electron );
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
  tree->Branch( "residual_x", &event.residual_x );
  tree->Branch( "residual_y", &event.residual_y );
  tree->Branch( "residual_z", &event.residual_z );
  tree->Branch( "resolution_x", &event.resolution_x );
  tree->Branch( "resolution_y", &event.resolution_y );
  tree->Branch( "resolution_z", &event.resolution_z );
  tree->Branch( "helix_t", &event.helix_t );
  tree->Branch( "alpha", &event.alpha);
  tree->Branch( "pathhit", &event.pathhit);
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_size", &event.track_cluster_size);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch( "track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch( "track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch( "track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch( "track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch( "track_cluster_row_center", &event.track_cluster_row_center);

  tree->Branch( "isgoodTPCKurama", &event.isgoodTPCKurama);
  tree->Branch( "insideTPC", &event.insideTPC);
  tree->Branch( "pTPCKurama", &event.pTPCKurama);
  tree->Branch( "qTPCKurama", &event.qTPCKurama);
  tree->Branch( "m2TPCKurama", &event.m2TPCKurama);
  tree->Branch( "xsTPC", &event.xsTPC);
  tree->Branch( "ysTPC", &event.ysTPC);
  tree->Branch( "usTPC", &event.usTPC);
  tree->Branch( "vsTPC", &event.vsTPC);
  tree->Branch( "pK18", &event.pK18);
  tree->Branch( "xbTPC", &event.xbTPC);
  tree->Branch( "ybTPC", &event.ybTPC);
  tree->Branch( "ubTPC", &event.ubTPC);
  tree->Branch( "vbTPC", &event.vbTPC);

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

  //track fitting results
  tree->Branch("GFstatus", &event.GFstatus);
  tree->Branch("GFntTpc", &event.GFntTpc);
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

  tree->Branch("GFmom_p", &event.GFmom_p);
  tree->Branch("GFtracklen_p", &event.GFtracklen_p);
  tree->Branch("GFtrack2vtxdist_p", &event.GFtrack2vtxdist_p);
  tree->Branch("GFcalctof_p", &event.GFcalctof_p);
  tree->Branch("GFsegHtof_p", &event.GFsegHtof_p);
  tree->Branch("GFtofHtof_p", &event.GFtofHtof_p);
  tree->Branch("GFtdiffHtof_p", &event.GFtdiffHtof_p);
  tree->Branch("GFposHtof_p", &event.GFposHtof_p);
  tree->Branch("GFposx_p", &event.GFposx_p);
  tree->Branch("GFposy_p", &event.GFposy_p);
  tree->Branch("GFposz_p", &event.GFposz_p);
  tree->Branch("GFinvbeta_p", &event.GFinvbeta_p);
  tree->Branch("GFm2_p", &event.GFm2_p);

  tree->Branch("GFmom_pi", &event.GFmom_pi);
  tree->Branch("GFtracklen_pi", &event.GFtracklen_pi);
  tree->Branch("GFtrack2vtxdist_pi", &event.GFtrack2vtxdist_pi);
  tree->Branch("GFcalctof_pi", &event.GFcalctof_pi);
  tree->Branch("GFsegHtof_pi", &event.GFsegHtof_pi);
  tree->Branch("GFtofHtof_pi", &event.GFtofHtof_pi);
  tree->Branch("GFtdiffHtof_pi", &event.GFtdiffHtof_pi);
  tree->Branch("GFposHtof_pi", &event.GFposHtof_pi);
  tree->Branch("GFposx_pi", &event.GFposx_pi);
  tree->Branch("GFposy_pi", &event.GFposy_pi);
  tree->Branch("GFposz_pi", &event.GFposz_pi);
  tree->Branch("GFinvbeta_pi", &event.GFinvbeta_pi);
  tree->Branch("GFm2_pi", &event.GFm2_pi);

  TTreeReaderCont[kHTOFCaib] = new TTreeReader( "tpc", TFileCont[kHTOFCaib] );
  const auto& reader = TTreeReaderCont[kHTOFCaib];
  src.runnum = new TTreeReaderValue<Int_t>( *reader, "runnum" );
  src.evnum = new TTreeReaderValue<Int_t>( *reader, "evnum" );
  src.trigpat = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trigpat" );
  src.trigflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trigflag" );

  src.nclTpc = new TTreeReaderValue<Int_t>( *reader, "nclTpc" );
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

  src.ntTpc = new TTreeReaderValue<Int_t>( *reader, "ntTpc" );
  src.nhtrack = new TTreeReaderValue<std::vector<Int_t>>( *reader, "nhtrack" );
  src.isBeam = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isBeam" );
  src.isK18 = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isK18" );
  src.isKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isKurama" );
  src.isAccidental = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isAccidental" );
  src.chisqr = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqr" );
  src.helix_cx = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_cx" );
  src.helix_cy = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_cy" );
  src.helix_z0 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_z0" );
  src.helix_r = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_r" );
  src.helix_dz = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_dz" );
  src.dE = new TTreeReaderValue<std::vector<Double_t>>( *reader, "dE" );
  src.dEdx = new TTreeReaderValue<std::vector<Double_t>>( *reader, "dEdx" );
  src.mom0 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "mom0" );
  src.charge = new TTreeReaderValue<std::vector<Int_t>>( *reader, "charge" );
  src.path = new TTreeReaderValue<std::vector<Double_t>>( *reader, "path" );
  src.pid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "pid" );

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
  src.alpha = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "alpha" );
  src.pathhit = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "pathhit" );
  src.track_cluster_de = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_de" );
  src.track_cluster_size = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_size" );
  src.track_cluster_mrow = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_mrow" );
  src.track_cluster_de_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_de_center" );
  src.track_cluster_x_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_x_center" );
  src.track_cluster_y_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_y_center" );
  src.track_cluster_z_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_z_center" );
  src.track_cluster_row_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_row_center" );

  src.isgoodTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isgoodTPCKurama" );
  src.insideTPC = new TTreeReaderValue<std::vector<Int_t>>( *reader, "insideTPC" );
  src.pTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pTPCKurama" );
  src.qTPCKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "qTPCKurama" );
  src.m2TPCKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "m2TPCKurama" );
  src.xsTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xsTPC" );
  src.ysTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ysTPC" );
  src.usTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "usTPC" );
  src.vsTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vsTPC" );

  src.pK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pK18" );
  src.xbTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xbTPC" );
  src.ybTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ybTPC" );
  src.ubTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ubTPC" );
  src.vbTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vbTPC" );

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

  ////////// Bring Address From Dst
  TTreeCont[kHodoscope1]->SetBranchStatus("*", 0);
  TTreeCont[kHodoscope1]->SetBranchStatus("CTime0",   1);
  TTreeCont[kHodoscope1]->SetBranchStatus("nhHtof",   1);
  TTreeCont[kHodoscope1]->SetBranchStatus("csHtof",   1);
  TTreeCont[kHodoscope1]->SetBranchStatus("HtofSeg",  1);
  TTreeCont[kHodoscope1]->SetBranchStatus("tHtof",    1);
  TTreeCont[kHodoscope1]->SetBranchStatus("dtHtof",   1);
  TTreeCont[kHodoscope1]->SetBranchStatus("deHtof",   1);
  TTreeCont[kHodoscope1]->SetBranchStatus("posHtof",   1);

  TTreeCont[kHodoscope1]->SetBranchAddress("CTime0",   &src.CTime0);
  TTreeCont[kHodoscope1]->SetBranchAddress("nhHtof",   &src.nhHtof);
  TTreeCont[kHodoscope1]->SetBranchAddress("HtofSeg",   src.HtofSeg);
  TTreeCont[kHodoscope1]->SetBranchAddress("tHtof",     src.tHtof);
  TTreeCont[kHodoscope1]->SetBranchAddress("dtHtof",    src.dtHtof);
  TTreeCont[kHodoscope1]->SetBranchAddress("deHtof",    src.deHtof);
  TTreeCont[kHodoscope1]->SetBranchAddress("posHtof",    src.posHtof);

  TTreeCont[kHodoscope2]->SetBranchStatus("*", 0);
  TTreeCont[kHodoscope2]->SetBranchStatus("htofmt",   1);
  TTreeCont[kHodoscope2]->SetBranchStatus("htofde",   1);
  TTreeCont[kHodoscope2]->SetBranchStatus("htofutime",   1);
  TTreeCont[kHodoscope2]->SetBranchStatus("htofuctime",   1);
  TTreeCont[kHodoscope2]->SetBranchStatus("htofdtime",   1);
  TTreeCont[kHodoscope2]->SetBranchStatus("htofdctime",   1);
  TTreeCont[kHodoscope2]->SetBranchStatus("htofhitpos",   1);
  TTreeCont[kHodoscope2]->SetBranchStatus("htofude",   1);
  TTreeCont[kHodoscope2]->SetBranchStatus("htofdde",   1);

  TTreeCont[kHodoscope2]->SetBranchAddress("htofmt", src.htofmt);
  TTreeCont[kHodoscope2]->SetBranchAddress("htofde", src.htofde);
  TTreeCont[kHodoscope2]->SetBranchAddress("htofutime", src.htofutime);
  TTreeCont[kHodoscope2]->SetBranchAddress("htofuctime", src.htofuctime);
  TTreeCont[kHodoscope2]->SetBranchAddress("htofdtime", src.htofdtime);
  TTreeCont[kHodoscope2]->SetBranchAddress("htofdctime", src.htofdctime);
  TTreeCont[kHodoscope2]->SetBranchAddress("htofhitpos", src.htofhitpos);
  TTreeCont[kHodoscope2]->SetBranchAddress("htofude", src.htofude);
  TTreeCont[kHodoscope2]->SetBranchAddress("htofdde", src.htofdde);

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

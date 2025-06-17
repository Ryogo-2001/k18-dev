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

  const Double_t vtx_scan_range = 150.; //ref
  const Double_t vtx_scan_rangeInsideL = 50.;
  const Double_t vtx_scan_rangeInsidePi = 50.;

  const Double_t lambda_masscut = 0.07;
  const Double_t p_vtx_distcut = 300;
  const Double_t pi_vtx_distcut = 300;
  const Double_t pipi_distcut = 10.; //ref
  const Double_t ppi_distcut = 10.; //ref
  const Double_t ltarget_distcut = 25.;

  const Double_t GFppi_distcut = 10.;
  const Double_t GFltarget_distcut = 25.;
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

  Bool_t lflag;
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
  Double_t GFltracklen;
  Double_t GFltof;

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
  //GFTrackCont.Init();

  HF1( 1, event.status++ );
  for(int it=0; it<event.ntTpc; ++it){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix(it);
    if( !tp ) continue;
    std::vector<Int_t> pdgcode;
    Kinematics::HypTPCPID_PDGCode(event.charge[it], event.pid[it], pdgcode);
    GFTrackCont.AddHelixTrack(pdgcode, tp);
  }
  GFTrackCont.FitTracks();

  HF1( 2, event.GFstatus++ );

  int GFntTpc = GFTrackCont.GetNTrack();
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
    event.GFfitstatus[igf] = (int)GFTrackCont.TrackCheck(igf);
    HF1( 3, event.GFfitstatus[igf]);
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

    HF1( genfitHid+1, event.GFchisqr[igf]);
    HF1( genfitHid+2, event.GFpval[igf]);
    HF1( genfitHid+3, event.GFcharge[igf]);
    HF1( genfitHid+4, event.GFnhtrack[igf]);
    HF1( genfitHid+5, event.GFtracklen[igf]);
    HF1( genfitHid+6, event.GFtof[igf]);
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
	GFTrackCont.TPCHTOFTrackMatching(igf, repid, vertex,
				      event.HtofSeg, event.posHtof,
				      hitid_htof, tof,
				      len, pos_htof, track2tgt_dist);
      if(htofextrapolation){
	event.GFmom_pi[igf] = GFTrackCont.GetMom(igf, 0, repid).Mag();
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
	GFTrackCont.TPCHTOFTrackMatching(igf, repid, vertex,
				      event.HtofSeg, event.posHtof,
				      hitid_htof, tof,
				      len, pos_htof, track2tgt_dist);
      if(htofextrapolation){
	event.GFmom_p[igf] = GFTrackCont.GetMom(igf, 0, repid).Mag();
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

  //for L, LL
  //Xi/L candidates searching
#if DebugDisp
  std::cout<<"1. L candidate searching"<<std::endl;
#endif
  Int_t l_candidates = 0;

  std::vector<TVector3> l_mom_container, l_vert_container;
  std::vector<Double_t> ppi_closedist;
  std::vector<Int_t> L_p_id_container, L_pi_id_container;
  std::vector<Int_t> L_p_repid_container, L_pi_repid_container;
  std::vector<TVector3> L_p_mom_container, L_pi_mom_container;
  std::vector<Double_t> L_mass_container;
  std::vector<Double_t> L_ppidist_container;
  std::vector<Double_t> L_targetdist_container;
  std::vector<TVector3> L_mom_container, L_vtx_container, L_targetvtx_container;

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

  {
    for(Int_t it1=0;it1<ntTpc;it1++){
      if(event.isElectron[it1]==1) continue;
      if(event.isBeam[it1]==1) continue;
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
	if(event.isBeam[it2]==1) continue;
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
      }
    }
  }

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
      if(event.isBeam[it1]==1) continue;
      if(event.isAccidental[it1]==1) continue;
      if((event.pid[it1]&1)!=1) continue; //pi+ like
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
	if(event.isBeam[it2]==1) continue;
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
    event.pipidistPipair.push_back(pipair_pipidist_container[icombi]);
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

  tree->Branch("Lflag", &event.lflag);

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

  tree->Branch("LambdaMass", &event.lmass);
  tree->Branch("LambdaDecayVtx_x", &event.ldecayvtx_x);
  tree->Branch("LambdaDecayVtx_y", &event.ldecayvtx_y);
  tree->Branch("LambdaDecayVtx_z", &event.ldecayvtx_z);
  tree->Branch("LambdaMom", &event.lmom);
  tree->Branch("LambdaMom_x", &event.lmom_x);
  tree->Branch("LambdaMom_y", &event.lmom_y);
  tree->Branch("LambdaMom_z", &event.lmom_z);
  tree->Branch("LambdaVtxCloseDist", &event.ppi_dist);
  tree->Branch("LambdaTargetCloseDist", &event.ltarget_dist);
  tree->Branch("LambdaTargetCloseVtx_x", &event.ltargetvtx_x);
  tree->Branch("LambdaTargetCloseVtx_y", &event.ltargetvtx_y);
  tree->Branch("LambdaTargetCloseVtx_z", &event.ltargetvtx_z);

  tree->Branch("DecaysTrackId", &event.decays_id);
  tree->Branch("DecaysMom", &event.decays_mom);
  tree->Branch("DecaysMom_x", &event.decays_mom_x);
  tree->Branch("DecaysMom_y", &event.decays_mom_y);
  tree->Branch("DecaysMom_z", &event.decays_mom_z);

  tree->Branch("GFDecaysHtofId", &event.GFdecays_htofid);
  tree->Branch("GFDecaysTrackLen", &event.GFdecays_tracklen);
  tree->Branch("GFDecaysTrackTof", &event.GFdecays_tof);
  tree->Branch("GFDecaysMassSquare", &event.GFdecays_mass2);
  tree->Branch("GFDecaysInvbeta", &event.GFdecays_invbeta);
  tree->Branch("GFDecaysMom", &event.GFdecays_mom);
  tree->Branch("GFDecaysMom_x", &event.GFdecays_mom_x);
  tree->Branch("GFDecaysMom_y", &event.GFdecays_mom_y);
  tree->Branch("GFDecaysMom_z", &event.GFdecays_mom_z);
  tree->Branch("GFDecaysMomLoss", &event.GFdecays_momloss);
  tree->Branch("GFDecaysELoss", &event.GFdecays_eloss);

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

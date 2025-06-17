// -*- C++ -*-

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <TLorentzVector.h>

#include <TGeoPhysicalConstants.h>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "FieldMan.hh"
#include "DatabasePDG.hh"
#include "DetectorID.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "HodoPHCMan.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "TPCAnalyzer.hh"
#include "TPCParamMan.hh"
#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCVertex.hh"
#include "TPCRKTrack.hh"
#include "TPCHit.hh"
#include "TPCPadHelper.hh"
#include <Math/ProbFunc.h>

#include "DstHelper.hh"
#include <TRandom.h>
#include "DebugCounter.hh"

#define ExclusiveResidual 0
#define KuramaK18 0

namespace
{
  using namespace root;
  using namespace dst;
  using namespace std;
  const std::string& class_name("DstTPCTrackingHelixGeant4");
  using hddaq::unpacker::GUnpacker;
  const auto qnan = TMath::QuietNaN();
  const auto& gUnpacker = GUnpacker::get_instance();
  ConfMan&            gConf = ConfMan::GetInstance();
  const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const HodoPHCMan&   gPHC  = HodoPHCMan::GetInstance();
  const auto& gTPC  = TPCParamMan::GetInstance();
  debug::ObjectCounter& gCounter  = debug::ObjectCounter::GetInstance();

  const Int_t MaxTPCHits = 10000;
  const Int_t MaxTPCTracks = 100;
  const Int_t MaxTPCnHits = 50;
  const Int_t MaxG4Hits = 500;
  const double truncatedMean = 0.8; //80%
  const bool IsWithRes = true;
  //const bool IsWithRes = true;
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
  std::vector<TTreeReader*> TTreeReaderCont;
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
  Int_t GetHitG4Tid(TVector3 TPCHit, int nhG4,int* tidG4, double* xG4,double* yG4,double* zG4){
    std::vector<TVector3> G4Hits;
    for(int ih=0;ih<nhG4;++ih){
      TVector3 G4Hit(xG4[ih],yG4[ih],zG4[ih]);
      G4Hits.push_back(G4Hit);
    }
    double dl = 5000;
    int G4ID = -1;
    int ih=0;
    for(auto g4hit:G4Hits){
      double dist = (g4hit - TPCHit).Mag();
      if(dist < dl){
	dl = dist;
	G4ID = tidG4[ih];
      }
      ih++;
    }
    return G4ID;
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

//_____________________________________________________________________
struct Event
{

  Int_t status;
  Int_t evnum;
  Int_t nhittpc;                 // Number of Hits
  Int_t max_ititpc;
  Int_t ititpc[MaxTPCHits];
  Int_t nhittpc_iti[MaxTPCHits];



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
  std::vector<Int_t> cluster_G4tid;


  Int_t ntTpc;                   // Number of Tracks
  Int_t ntKuramaCandidate; //Numer of tracks which are kurama track candidates(before TPCKurama tracking)
  std::vector<Int_t> isKuramaCandidate;
  vector<Int_t> nhtrack;
  vector<Int_t> nhtrackEff;
  vector<Int_t> trackid; //for Kurama & K1.8 tracks
  vector<Int_t> isBeam;
  vector<Int_t> isXi;
  vector<Int_t> isKurama;
  vector<Int_t> isK18;
  vector<Int_t> isAccidental;
  vector<Int_t> isMultiloop;
  vector<Int_t> flag;
  vector<Double_t> purity;
  vector<Double_t> efficiency;
  vector<Int_t> G4tid;
  vector<Double_t> chisqr;
  vector<Double_t> pval;
  vector<Double_t> distTgt;
  vector<Double_t> helix_cx;
  vector<Double_t> helix_cy;
  vector<Double_t> helix_z0;
  vector<Double_t> helix_r;
  vector<Double_t> helix_dz;
  vector<Double_t> dE;
  vector<Double_t> dEdx; //reference dedx
  vector<Double_t> dz_factor;
  vector<Double_t> mom0;//Helix momentum at Y = 0
  vector<Double_t> mom0_x;//Helix momentum at Y = 0
  vector<Double_t> mom0_y;
  vector<Double_t> mom0_z;
  vector<Int_t> charge;
  vector<Int_t> pid;
  vector<Double_t> path;
  std::vector<Int_t> isElectron;
  std::vector<Double_t> nsigma_triton;
  std::vector<Double_t> nsigma_deutron;
  std::vector<Double_t> nsigma_proton;
  std::vector<Double_t> nsigma_kaon;
  std::vector<Double_t> nsigma_pion;
  std::vector<Double_t> nsigma_electron;

  vector<vector<Double_t>> hitlayer;
  vector<vector<Double_t>> hitpos_x;
  vector<vector<Double_t>> hitpos_y;
  vector<vector<Double_t>> hitpos_z;
  vector<vector<Double_t>> calpos_x;
  vector<vector<Double_t>> calpos_y;
  vector<vector<Double_t>> calpos_z;
  vector<vector<Double_t>> mom_x;
  vector<vector<Double_t>> mom_y;
  vector<vector<Double_t>> mom_z;
  vector<vector<Double_t>> residual;
  vector<vector<Double_t>> residual_t;
  vector<vector<Double_t>> residual_x;
  vector<vector<Double_t>> residual_y;
  vector<vector<Double_t>> residual_z;
  vector<vector<Double_t>> resolution;
  vector<vector<Double_t>> resolution_x;
  vector<vector<Double_t>> resolution_y;
  vector<vector<Double_t>> resolution_z;
  vector<vector<Double_t>> pull;
  vector<vector<Double_t>> helix_t;
  vector<vector<Double_t>> pathhit;
  vector<vector<Double_t>> alpha;
  vector<vector<Double_t>> houghflag;
  vector<vector<Double_t>> track_cluster_size;
  vector<vector<Double_t>> track_cluster_de;
  vector<vector<Double_t>> track_cluster_mrow;
  vector<vector<Double_t>> pull_t;
  vector<vector<Double_t>> pull_x;
  vector<vector<Double_t>> pull_y;
  vector<vector<Double_t>> pull_z;

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

  //Geant4
  Int_t iti_g[MaxTPCTracks][MaxTPCnHits];
  Int_t idtpc[MaxTPCHits];
  Int_t ID[MaxTPCHits];
  Int_t PID[MaxTPCHits];
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
  Double_t momg_x[MaxTPCTracks][MaxTPCnHits];
  Double_t momg_y[MaxTPCTracks][MaxTPCnHits];
  Double_t momg_z[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_p[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_px[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_py[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_pz[MaxTPCTracks][MaxTPCnHits];

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

  int G4idKm,G4idKp,G4idP,G4idPi1,G4idPi2;;

  Double_t MomXi_x,MomXi_y,MomXi_z;//At decay vtx
  Double_t SpinXi_x,SpinXi_y,SpinXi_z;
  Double_t ThXi_CM;
  Double_t MomLd_x,MomLd_y,MomLd_z;
  Double_t SpinLd_x,SpinLd_y,SpinLd_z;
  Double_t ThLd_CM;

  //KpXi Kinematics, Momentum at production vtx
  double PKm,PKm_x,PKm_y,PKm_z;
  double PKp,PKp_x,PKp_y,PKp_z;
  double PXi,PXi_x,PXi_y,PXi_z;
  double PPi2,PPi2_x,PPi2_y,PPi2_z;//Pi from Xi

  double PLd,PLd_x,PLd_y,PLd_z;
  double PPi1,PPi1_x,PPi1_y,PPi1_z;//Pi from Ld
  double PP,PP_x,PP_y,PP_z;

  //HToF
  Int_t nhHtof;
  std::vector<Double_t> HtofSeg;
  std::vector<Double_t> tHtof;
  std::vector<Double_t> dtHtof;
  std::vector<Double_t> deHtof;
  std::vector<Double_t> posHtof;
  std::vector<Int_t> G4tidHtof;

  Int_t nhSch;
  Int_t tidSch[MaxG4Hits];
  Int_t pidSch[MaxG4Hits];
  Int_t didSch[MaxG4Hits];
  Int_t prtSch[MaxG4Hits];
  Int_t qSch[MaxG4Hits];
  Double_t xSch[MaxG4Hits];
  Double_t ySch[MaxG4Hits];
  Double_t zSch[MaxG4Hits];
  Double_t pxSch[MaxG4Hits];
  Double_t pySch[MaxG4Hits];
  Double_t pzSch[MaxG4Hits];
  Double_t ppSch[MaxG4Hits];
  Double_t deSch[MaxG4Hits];
  Double_t tSch[MaxG4Hits];

  //Lac
  Int_t nhLac;
  Int_t tidLac[MaxG4Hits];
  Int_t pidLac[MaxG4Hits];
  Int_t didLac[MaxG4Hits];
  Int_t prtLac[MaxG4Hits];
  Int_t qLac[MaxG4Hits];
  Double_t xLac[MaxG4Hits];
  Double_t yLac[MaxG4Hits];
  Double_t zLac[MaxG4Hits];
  Double_t pxLac[MaxG4Hits];
  Double_t pyLac[MaxG4Hits];
  Double_t pzLac[MaxG4Hits];
  Double_t ppLac[MaxG4Hits];
  Double_t deLac[MaxG4Hits];
  Double_t tLac[MaxG4Hits];

  //FToF
  Int_t nhFtof;
  Int_t tidFtof[MaxG4Hits];
  Int_t pidFtof[MaxG4Hits];
  Int_t didFtof[MaxG4Hits];
  Int_t prtFtof[MaxG4Hits];
  Int_t qFtof[MaxG4Hits];
  Double_t xFtof[MaxG4Hits];
  Double_t yFtof[MaxG4Hits];
  Double_t zFtof[MaxG4Hits];
  Double_t pxFtof[MaxG4Hits];
  Double_t pyFtof[MaxG4Hits];
  Double_t pzFtof[MaxG4Hits];
  Double_t ppFtof[MaxG4Hits];
  Double_t deFtof[MaxG4Hits];
  Double_t tFtof[MaxG4Hits];

  //SDC
  Int_t nhSdc;
  Int_t tidSdc[MaxG4Hits];
  Int_t pidSdc[MaxG4Hits];
  Int_t didSdc[MaxG4Hits];
  Int_t prtSdc[MaxG4Hits];
  Int_t qSdc[MaxG4Hits];
  Double_t xSdc[MaxG4Hits];
  Double_t ySdc[MaxG4Hits];
  Double_t zSdc[MaxG4Hits];
  Double_t pxSdc[MaxG4Hits];
  Double_t pySdc[MaxG4Hits];
  Double_t pzSdc[MaxG4Hits];
  Double_t ppSdc[MaxG4Hits];
  Double_t deSdc[MaxG4Hits];
  Double_t tSdc[MaxG4Hits];

  //WC
  Int_t nhWc;
  Int_t tidWc[MaxG4Hits];
  Int_t pidWc[MaxG4Hits];
  Int_t didWc[MaxG4Hits];
  Int_t prtWc[MaxG4Hits];
  Int_t qWc[MaxG4Hits];
  Double_t xWc[MaxG4Hits];
  Double_t yWc[MaxG4Hits];
  Double_t zWc[MaxG4Hits];
  Double_t pxWc[MaxG4Hits];
  Double_t pyWc[MaxG4Hits];
  Double_t pzWc[MaxG4Hits];
  Double_t ppWc[MaxG4Hits];
  Double_t deWc[MaxG4Hits];
  Double_t tWc[MaxG4Hits];

  //BVH
  Int_t nhBvh;
  Int_t tidBvh[MaxG4Hits];
  Int_t pidBvh[MaxG4Hits];
  Int_t didBvh[MaxG4Hits];
  Int_t prtBvh[MaxG4Hits];
  Int_t qBvh[MaxG4Hits];
  Double_t xBvh[MaxG4Hits];
  Double_t yBvh[MaxG4Hits];
  Double_t zBvh[MaxG4Hits];
  Double_t pxBvh[MaxG4Hits];
  Double_t pyBvh[MaxG4Hits];
  Double_t pzBvh[MaxG4Hits];
  Double_t ppBvh[MaxG4Hits];
  Double_t deBvh[MaxG4Hits];
  Double_t tBvh[MaxG4Hits];

  //RealData
  int ntK18;
  vector<vector<double>>xvpHS;
  vector<vector<double>>yvpHS;
  vector<vector<double>>zvpHS;
  vector<double>xtgtHS;
  vector<double>ytgtHS;
  vector<double>ztgtHS;
  vector<double>xoutK18;
  vector<double>youtK18;
  vector<double>uoutK18;
  vector<double>voutK18;
  vector<double>p_3rd;
  vector<vector<double>> layerK18;
  vector<vector<double>> wireK18;
  vector<vector<double>> localhitposK18;
  int ntKurama;
  vector<vector<double>>xvpKurama;
  vector<vector<double>>yvpKurama;
  vector<vector<double>>zvpKurama;
  vector<double>xtgtKurama;
  vector<double>ytgtKurama;
  vector<double>xout;
  vector<double>yout;
  vector<double>zout;
  vector<double>pxout;
  vector<double>pyout;
  vector<double>pzout;
  vector<vector<double>> wire;
  vector<vector<double>> layer;
  vector<vector<double>> localhitpos;



  vector<Int_t> isgoodTPCK18;
  vector<Int_t> tpcidTPCK18;
  vector<Int_t> niterationTPCK18;
  vector<Double_t> chisqrTPCK18;

  vector<Int_t> isgoodTPCKurama;
  vector<Int_t> tpcidTPCKurama;
  vector<Int_t> niterationTPCKurama;
  vector<Int_t> kflagTPCKurama;
  vector<Double_t> chisqrTPCKurama;



  void Clear(){

    HtofSeg.clear();
    tHtof.clear();
    dtHtof.clear();
    deHtof.clear();
    posHtof.clear();
    G4tidHtof.clear();

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
    cluster_G4tid.clear();


    PKm = qnan;PKm_x=qnan;PKm_y=qnan;PKm_z=qnan;
    PKp = qnan;PKp_x=qnan;PKp_y=qnan;PKp_z=qnan;
    PXi = qnan;PXi_x=qnan;PXi_y=qnan;PXi_z=qnan;
    PPi2 = qnan;PPi2_x=qnan;PPi2_y=qnan;PPi2_z=qnan;
    PLd = qnan;PLd_x=qnan;PLd_y=qnan;PLd_z=qnan;
    PPi1 = qnan;PPi1_x=qnan;PPi1_y=qnan;PPi1_z=qnan;
    PP = qnan;PP_x=qnan;PP_y=qnan;PP_z=qnan;
    nhtrack.clear();
    nhtrackEff.clear();
    trackid.clear();
    isBeam.clear();
    isKurama.clear();
    isK18.clear();
    isAccidental.clear();
    isMultiloop.clear();
    flag.clear();
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
    dz_factor.clear();
    dE.clear();
    dEdx.clear();
    mom0_x.clear();
    mom0_y.clear();
    mom0_z.clear();
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
    helix_t.clear();
    calpos_x.clear();
    calpos_y.clear();
    calpos_z.clear();
    mom_x.clear();
    mom_y.clear();
    mom_z.clear();

    resolution.clear();
    resolution_x.clear();
    resolution_y.clear();
    resolution_z.clear();
    residual.clear();
    residual_t.clear();
    residual_x.clear();
    residual_y.clear();
    residual_z.clear();
    pull_t.clear();
    pull_x.clear();
    pull_y.clear();
    pull_z.clear();
    pathhit.clear();
    alpha.clear();
    houghflag.clear();

    track_cluster_size.clear();
    track_cluster_de.clear();
    track_cluster_mrow.clear();


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

    ntK18 = 0;
    xvpHS.clear();
    yvpHS.clear();
    zvpHS.clear();
    xoutK18.clear();
    youtK18.clear();
    uoutK18.clear();
    voutK18.clear();
    p_3rd.clear();
    layerK18.clear();
    wireK18.clear();
    localhitposK18.clear();
    ntKurama = 0;
    xvpKurama.clear();
    yvpKurama.clear();
    zvpKurama.clear();
    xout.clear();
    yout.clear();
    zout.clear();
    pxout.clear();
    pyout.clear();
    pzout.clear();
    layer.clear();
    wire.clear();
    localhitpos.clear();

    isgoodTPCK18.clear();
    tpcidTPCK18.clear();
    niterationTPCK18.clear();
    chisqrTPCK18.clear();

    isgoodTPCKurama.clear();
    tpcidTPCKurama.clear();
    niterationTPCKurama.clear();
    kflagTPCKurama.clear();
    chisqrTPCKurama.clear();

  };
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
  Double_t ppPrm[MaxTPCTracks];

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
  Double_t MomXi_x,MomXi_y,MomXi_z;
  Double_t SpinXi_x,SpinXi_y,SpinXi_z;
  Double_t ThXi_CM;
  Double_t MomLd_x,MomLd_y,MomLd_z;
  Double_t SpinLd_x,SpinLd_y,SpinLd_z;
  Double_t ThLd_CM;


  Int_t ititpc[MaxTPCHits];
  Int_t idtpc[MaxTPCHits];
  Int_t parentID[MaxTPCHits];

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
  Double_t pptpc[MaxTPCHits];   // total mometum

  // Double_t masstpc[MaxTPCHits];   // mass TPC
  // Double_t betatpc[MaxTPCHits];
  Double_t edeptpc[MaxTPCHits];
  // Int_t laytpc[MaxTPCHits];
  // Int_t rowtpc[MaxTPCHits];
  // Int_t iPadtpc[MaxTPCHits];//Pad number (0 origin)

  //HToF
  Int_t nhHtof;
  Int_t tidHtof[MaxG4Hits];
  Int_t pidHtof[MaxG4Hits];
  Int_t didHtof[MaxG4Hits];
  Int_t prtHtof[MaxG4Hits];
  Int_t qHtof[MaxG4Hits];
  Double_t xHtof[MaxG4Hits];
  Double_t yHtof[MaxG4Hits];
  Double_t zHtof[MaxG4Hits];
  Double_t pxHtof[MaxG4Hits];
  Double_t pyHtof[MaxG4Hits];
  Double_t pzHtof[MaxG4Hits];
  Double_t ppHtof[MaxG4Hits];
  Double_t deHtof[MaxG4Hits];
  Double_t tHtof[MaxG4Hits];

  //SCH
  Int_t nhSch;
  Int_t tidSch[MaxG4Hits];
  Int_t pidSch[MaxG4Hits];
  Int_t didSch[MaxG4Hits];
  Int_t prtSch[MaxG4Hits];
  Int_t qSch[MaxG4Hits];
  Double_t xSch[MaxG4Hits];
  Double_t ySch[MaxG4Hits];
  Double_t zSch[MaxG4Hits];
  Double_t pxSch[MaxG4Hits];
  Double_t pySch[MaxG4Hits];
  Double_t pzSch[MaxG4Hits];
  Double_t ppSch[MaxG4Hits];
  Double_t deSch[MaxG4Hits];
  Double_t tSch[MaxG4Hits];

  //Lac
  Int_t nhLac;
  Int_t tidLac[MaxG4Hits];
  Int_t pidLac[MaxG4Hits];
  Int_t didLac[MaxG4Hits];
  Int_t prtLac[MaxG4Hits];
  Int_t qLac[MaxG4Hits];
  Double_t xLac[MaxG4Hits];
  Double_t yLac[MaxG4Hits];
  Double_t zLac[MaxG4Hits];
  Double_t pxLac[MaxG4Hits];
  Double_t pyLac[MaxG4Hits];
  Double_t pzLac[MaxG4Hits];
  Double_t ppLac[MaxG4Hits];
  Double_t deLac[MaxG4Hits];
  Double_t tLac[MaxG4Hits];

  //FToF
  Int_t nhFtof;
  Int_t tidFtof[MaxG4Hits];
  Int_t pidFtof[MaxG4Hits];
  Int_t didFtof[MaxG4Hits];
  Int_t prtFtof[MaxG4Hits];
  Int_t qFtof[MaxG4Hits];
  Double_t xFtof[MaxG4Hits];
  Double_t yFtof[MaxG4Hits];
  Double_t zFtof[MaxG4Hits];
  Double_t pxFtof[MaxG4Hits];
  Double_t pyFtof[MaxG4Hits];
  Double_t pzFtof[MaxG4Hits];
  Double_t ppFtof[MaxG4Hits];
  Double_t deFtof[MaxG4Hits];
  Double_t tFtof[MaxG4Hits];

  //SDC
  Int_t nhSdc;
  Int_t tidSdc[MaxG4Hits];
  Int_t pidSdc[MaxG4Hits];
  Int_t didSdc[MaxG4Hits];
  Int_t prtSdc[MaxG4Hits];
  Int_t qSdc[MaxG4Hits];
  Double_t xSdc[MaxG4Hits];
  Double_t ySdc[MaxG4Hits];
  Double_t zSdc[MaxG4Hits];
  Double_t pxSdc[MaxG4Hits];
  Double_t pySdc[MaxG4Hits];
  Double_t pzSdc[MaxG4Hits];
  Double_t ppSdc[MaxG4Hits];
  Double_t deSdc[MaxG4Hits];
  Double_t tSdc[MaxG4Hits];

  //WC
  Int_t nhWc;
  Int_t tidWc[MaxG4Hits];
  Int_t pidWc[MaxG4Hits];
  Int_t didWc[MaxG4Hits];
  Int_t prtWc[MaxG4Hits];
  Int_t qWc[MaxG4Hits];
  Double_t xWc[MaxG4Hits];
  Double_t yWc[MaxG4Hits];
  Double_t zWc[MaxG4Hits];
  Double_t pxWc[MaxG4Hits];
  Double_t pyWc[MaxG4Hits];
  Double_t pzWc[MaxG4Hits];
  Double_t ppWc[MaxG4Hits];
  Double_t deWc[MaxG4Hits];
  Double_t tWc[MaxG4Hits];

  //BVH
  Int_t nhBvh;
  Int_t tidBvh[MaxG4Hits];
  Int_t pidBvh[MaxG4Hits];
  Int_t didBvh[MaxG4Hits];
  Int_t prtBvh[MaxG4Hits];
  Int_t qBvh[MaxG4Hits];
  Double_t xBvh[MaxG4Hits];
  Double_t yBvh[MaxG4Hits];
  Double_t zBvh[MaxG4Hits];
  Double_t pxBvh[MaxG4Hits];
  Double_t pyBvh[MaxG4Hits];
  Double_t pzBvh[MaxG4Hits];
  Double_t ppBvh[MaxG4Hits];
  Double_t deBvh[MaxG4Hits];
  Double_t tBvh[MaxG4Hits];


  //RealData
  TTreeReaderValue<int>* ntK18=nullptr;
  TTreeReaderValue<vector<vector<double>>>* xvpHS=nullptr;
  TTreeReaderValue<vector<vector<double>>>* yvpHS=nullptr;
  TTreeReaderValue<vector<vector<double>>>* zvpHS=nullptr;
  TTreeReaderValue<vector<double>>* xtgtHS=nullptr;
  TTreeReaderValue<vector<double>>* ytgtHS=nullptr;
  TTreeReaderValue<vector<double>>* ztgtHS=nullptr;
  TTreeReaderValue<vector<double>>* xoutK18=nullptr;
  TTreeReaderValue<vector<double>>* youtK18=nullptr;
  TTreeReaderValue<vector<double>>* uoutK18=nullptr;
  TTreeReaderValue<vector<double>>* voutK18=nullptr;
  TTreeReaderValue<vector<double>>* p_3rd=nullptr;
  TTreeReaderValue<vector<vector<double>>>* layerK18=nullptr;
  TTreeReaderValue<vector<vector<double>>>* wireK18=nullptr;
  TTreeReaderValue<vector<vector<double>>>* localhitposK18=nullptr;
  TTreeReaderValue<int>* ntKurama=nullptr;
  TTreeReaderValue<vector<vector<double>>>* xvpKurama=nullptr;
  TTreeReaderValue<vector<vector<double>>>* yvpKurama=nullptr;
  TTreeReaderValue<vector<vector<double>>>* zvpKurama=nullptr;
  TTreeReaderValue<vector<double>>* xtgtKurama=nullptr;
  TTreeReaderValue<vector<double>>* ytgtKurama=nullptr;
  TTreeReaderValue<vector<double>>* xout=nullptr;
  TTreeReaderValue<vector<double>>* yout=nullptr;
  TTreeReaderValue<vector<double>>* zout=nullptr;
  TTreeReaderValue<vector<double>>* pxout=nullptr;
  TTreeReaderValue<vector<double>>* pyout=nullptr;
  TTreeReaderValue<vector<double>>* pzout=nullptr;
  TTreeReaderValue<vector<vector<double>>>* layer=nullptr;
  TTreeReaderValue<vector<vector<double>>>* wire=nullptr;
  TTreeReaderValue<vector<vector<double>>>* localhitpos=nullptr;


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
  if( !gConf.InitializeUnpacker() )
    return EXIT_FAILURE;
  // int nevent = GetEntries( TTreeCont );

  // CatchSignal::Set();

  // int ievent = 0;

  // // for(int ii=0; ii<100; ++ii){
  // //   std::cout<<"ii="<<ii<<std::endl;
  // //   ievent = 0;
  // for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
  //   gCounter.check();
  //   InitializeEvent();
  //   if( DstRead( ievent ) ) tree->Fill();
  // }
  // //  }

  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries( TTreeCont );
  if (max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();
  gRandom->SetSeed(7);
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

//_____________________________________________________________________
bool
dst::InitializeEvent( void )
{
  event.Clear();
  event.status   = 0;
  event.evnum = 0;
  event.nhittpc = 0;
  event.ntTpc = 0;
  event.max_ititpc = 0;

  for(int i=0; i<MaxTPCTracks; ++i){
    event.ititpc[i] =0;
    event.nhittpc_iti[i] =0;
    for(int j=0; j<MaxTPCnHits; ++j){
      event.momg_x[i][j] =-9999.;
      event.momg_y[i][j] =-9999.;
      event.momg_z[i][j] =-9999.;
      event.residual_px[i][j] =-9999.;
      event.residual_py[i][j] =-9999.;
      event.residual_pz[i][j] =-9999.;
      event.residual_p[i][j] =-9999.;
      event.iti_g[i][j] =0;
    }
  }
  event.NumberOfTracks = -1;
  for(int i=0;i<1000;++i){
    event.PIDOfTrack[i]=qnan;
    event.ParentIDOfTrack[i]=qnan;
    event.VertexOfTrack_x[i]=qnan;
    event.VertexOfTrack_y[i]=qnan;
    event.VertexOfTrack_z[i]=qnan;
    event.MomentumOfTrack[i]=qnan;
    event.MomentumOfTrack_x[i]=qnan;
    event.MomentumOfTrack_y[i]=qnan;
    event.MomentumOfTrack_z[i]=qnan;
  }

  event.G4idKm = -1;
  event.G4idKp = -1;
  event.G4idP = -1;
  event.G4idPi1 = -1;
  event.G4idPi2 = -1;
  event.MomXi_x = 0;
  event.MomXi_y = 0;
  event.MomXi_z = 0;
  event.SpinXi_x = 0;
  event.SpinXi_y = 0;
  event.SpinXi_z = 0;
  event.ThXi_CM = 0;

  event.MomLd_x = 0;
  event.MomLd_y = 0;
  event.MomLd_z = 0;
  event.SpinLd_x = 0;
  event.SpinLd_y = 0;
  event.SpinLd_z = 0;
  event.ThLd_CM = 0;

  event.nhHtof=-1;
  event.nhSch = -1;
  event.nhLac = -1;
  event.nhFtof=-1;
  event.nhSdc=-1;
  event.nhWc =-1;
  event.nhBvh=-1;
  for(int i=0;i<500;++i){
    event.ititpc[i]=qnan;
    event.idtpc[i]=qnan;
    event.ID[i]=qnan;
    event.PID[i]=qnan;
    event.xtpc[i]=qnan;
    event.ytpc[i]=qnan;
    event.ztpc[i]=qnan;
    event.x0tpc[i]=qnan;
    event.y0tpc[i]=qnan;
    event.z0tpc[i]=qnan;
    event.pxtpc[i]=qnan;
    event.pytpc[i]=qnan;
    event.pztpc[i]=qnan;


    event.tidSch[i]=qnan;
    event.pidSch[i]=qnan;
    event.didSch[i]=qnan;
    event.prtSch[i]=qnan;
    event.qSch[i]=qnan;
    event.xSch[i]=qnan;
    event.ySch[i]=qnan;
    event.zSch[i]=qnan;
    event.pxSch[i]=qnan;
    event.pySch[i]=qnan;
    event.pzSch[i]=qnan;
    event.ppSch[i]=qnan;
    event.deSch[i]=qnan;
    event.tSch[i]=qnan;

    event.tidLac[i]=qnan;
    event.pidLac[i]=qnan;
    event.didLac[i]=qnan;
    event.prtLac[i]=qnan;
    event.qLac[i]=qnan;
    event.xLac[i]=qnan;
    event.yLac[i]=qnan;
    event.zLac[i]=qnan;
    event.pxLac[i]=qnan;
    event.pyLac[i]=qnan;
    event.pzLac[i]=qnan;
    event.ppLac[i]=qnan;
    event.deLac[i]=qnan;
    event.tLac[i]=qnan;

    event.tidFtof[i]=qnan;
    event.pidFtof[i]=qnan;
    event.didFtof[i]=qnan;
    event.prtFtof[i]=qnan;
    event.qFtof[i]=qnan;
    event.xFtof[i]=qnan;
    event.yFtof[i]=qnan;
    event.zFtof[i]=qnan;
    event.pxFtof[i]=qnan;
    event.pyFtof[i]=qnan;
    event.pzFtof[i]=qnan;
    event.ppFtof[i]=qnan;
    event.deFtof[i]=qnan;
    event.tFtof[i]=qnan;

    event.tidSdc[i]=qnan;
    event.pidSdc[i]=qnan;
    event.didSdc[i]=qnan;
    event.prtSdc[i]=qnan;
    event.qSdc[i]=qnan;
    event.xSdc[i]=qnan;
    event.ySdc[i]=qnan;
    event.zSdc[i]=qnan;
    event.pxSdc[i]=qnan;
    event.pySdc[i]=qnan;
    event.pzSdc[i]=qnan;
    event.ppSdc[i]=qnan;
    event.deSdc[i]=qnan;
    event.tSdc[i]=qnan;

    event.tidWc[i]=qnan;
    event.pidWc[i]=qnan;
    event.didWc[i]=qnan;
    event.prtWc[i]=qnan;
    event.qWc[i]=qnan;
    event.xWc[i]=qnan;
    event.yWc[i]=qnan;
    event.zWc[i]=qnan;
    event.pxWc[i]=qnan;
    event.pyWc[i]=qnan;
    event.pzWc[i]=qnan;
    event.ppWc[i]=qnan;
    event.deWc[i]=qnan;
    event.tWc[i]=qnan;

    event.tidBvh[i]=qnan;
    event.pidBvh[i]=qnan;
    event.didBvh[i]=qnan;
    event.prtBvh[i]=qnan;
    event.qBvh[i]=qnan;
    event.xBvh[i]=qnan;
    event.yBvh[i]=qnan;
    event.zBvh[i]=qnan;
    event.pxBvh[i]=qnan;
    event.pyBvh[i]=qnan;
    event.pzBvh[i]=qnan;
    event.ppBvh[i]=qnan;
    event.deBvh[i]=qnan;
    event.tBvh[i]=qnan;
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
  static const auto KaonMass    = pdg::KaonMass();
  static const auto PionMass    = pdg::PionMass();
  static const auto ProtonMass  = pdg::ProtonMass();
  static const auto LambdaMass  = pdg::LambdaMass();
  static const auto XiMass      = pdg::XiMinusMass();
  static const auto XiStarMass  = 1.5350;
  static const auto ElectronMass = pdg::ElectronMass();
  static const Double_t Carbon12Mass = 12.*TGeoUnit::amu_c2 - 6.*ElectronMass;
  static const auto xGlobalBcOut = gGeom.GetGlobalPosition("BC3-X1").X();
  static const auto yGlobalBcOut = gGeom.GetGlobalPosition("BC3-X1").Y();
  static const auto zGlobalBcOut = gGeom.GetGlobalPosition("BC3-X1").Z();
  static const auto zLocalBcOut = gGeom.GetLocalZ("BC3-X1");
  static const auto xGlobalSdcOut = gGeom.GetGlobalPosition("SDC4-X2").X();
  static const auto yGlobalSdcOut = gGeom.GetGlobalPosition("SDC4-X2").Y();




  if( ievent%1==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  HF1( 1, event.status++ );


  vector<TVector3> G4Hits;
  vector<TVector3> G4Moms;

  for(int ih=0;ih<src.nhittpc;++ih){
    TVector3 G4Hit(src.xtpc[ih],src.ytpc[ih],src.ztpc[ih]);
    TVector3 G4Mom(src.pxtpc[ih],src.pytpc[ih],src.pztpc[ih]);
    G4Hits.push_back(G4Hit);
    G4Moms.push_back(G4Mom);
  }

  int G4idKm = -1;
  int G4idKp = -1;
  event.evnum = src.evnum -1;
  event.nhittpc = src.nhittpc;
  //debug  std::cout<<"DstTPCTracking Helix Geant4, nhit="<<event.nhittpc<<std::endl;
  event.NumberOfTracks = src.NumberOfTracks;
  for(int it=0;it<1000;++it){
    event.PIDOfTrack[it]=src.PIDOfTrack[it];
    event.ParentIDOfTrack[it]=src.ParentIDOfTrack[it];
    event.VertexOfTrack_x[it]=src.VertexOfTrack_x[it];
    event.VertexOfTrack_y[it]=src.VertexOfTrack_y[it];
    event.VertexOfTrack_z[it]=src.VertexOfTrack_z[it];
    event.MomentumOfTrack[it]=src.MomentumOfTrack[it];
    event.MomentumOfTrack_x[it]=src.MomentumOfTrack_x[it];
    event.MomentumOfTrack_y[it]=src.MomentumOfTrack_y[it];
    event.MomentumOfTrack_z[it]=src.MomentumOfTrack_z[it];
  }
  for(int it=0;it<event.NumberOfTracks;++it){
    int tid = it + 1;
    int pid = event.PIDOfTrack[tid];
    int parent = event.ParentIDOfTrack[tid];
    int pid_parent = 0;
    if(parent!= -1) pid_parent = event.PIDOfTrack[parent];
    if(pid == 11) continue;
    double px = event.MomentumOfTrack_x[tid]/1000;//MeV to GeV
    double py = event.MomentumOfTrack_y[tid]/1000;
    double pz = event.MomentumOfTrack_z[tid]/1000;
    double p = hypot(pz,hypot(px,py));
    double vert_x = src.VertexOfTrack_x[it];
    double vert_y = src.VertexOfTrack_y[it];
    double vert_z = src.VertexOfTrack_z[it];
    if(abs(pid) == 321 and parent==0){
      if(pz < 0){
        G4idKm = tid;
        event.G4idKm = tid;
	event.PKm = p;
	event.PKm_x = -px;
	event.PKm_y = -py;
	event.PKm_z = -pz;
      }
      else{
        G4idKp = tid;
	event.G4idKp = tid;
        event.PKp = p;
	event.PKp_x = px;
	event.PKp_y = py;
	event.PKp_z = pz;
      }
    }
    if(pid == 3312){
      event.PXi = p;
      event.PXi_x = px;
      event.PXi_y = py;
      event.PXi_z = pz;
    }
    if(pid == 3122 and pid_parent == 3312){
      event.PLd = p;
      event.PLd_x = px;
      event.PLd_y = py;
      event.PLd_z = pz;
    }
    if(pid == 2212 and pid_parent == 3122){
      event.G4idP = tid;
    }
    if(pid == -211){
      if (pid_parent == 3122){
        event.G4idPi1 = tid;
      }
      else if(pid_parent == 3312){
        event.G4idPi2 = tid;
      }
    }
  }
  event.MomXi_x = src.MomXi_x;
  event.MomXi_y = src.MomXi_y;
  event.MomXi_z = src.MomXi_z;
  event.SpinXi_x = src.SpinXi_x;
  event.SpinXi_y = src.SpinXi_y;
  event.SpinXi_z = src.SpinXi_z;
  event.ThXi_CM = src.ThXi_CM;
  event.MomLd_x = src.MomLd_x;
  event.MomLd_y = src.MomLd_y;
  event.MomLd_z = src.MomLd_z;
  event.SpinLd_x = src.SpinLd_x;
  event.SpinLd_y = src.SpinLd_y;
  event.SpinLd_z = src.SpinLd_z;
  event.ThLd_CM = src.ThLd_CM;
  event.nhHtof = src.nhHtof;
  for(int ih=0;ih<event.nhHtof;++ih){
    event.HtofSeg.push_back(src.didHtof[ih]+1);
    event.tHtof.push_back(src.tHtof[ih]);
    event.deHtof.push_back(src.deHtof[ih]);
    event.posHtof.push_back(src.yHtof[ih]);
    event.G4tidHtof.push_back(src.tidHtof[ih]);
  }

  event.nhSch = src.nhSch;
  for(int ihit=0;ihit<event.nhSch;++ihit){
    event.tidSch[ihit]=src.tidSch[ihit];
    event.pidSch[ihit]=src.pidSch[ihit];
    event.didSch[ihit]=src.didSch[ihit];
    event.prtSch[ihit]=src.prtSch[ihit];
    event.qSch[ihit]=src.qSch[ihit];
    event.xSch[ihit]=src.xSch[ihit];
    event.ySch[ihit]=src.ySch[ihit];
    event.zSch[ihit]=src.zSch[ihit];
    event.pxSch[ihit]=src.pxSch[ihit];
    event.pySch[ihit]=src.pySch[ihit];
    event.pzSch[ihit]=src.pzSch[ihit];
    event.ppSch[ihit]=src.ppSch[ihit];
    event.deSch[ihit]=src.deSch[ihit];
    event.tSch[ihit]=src.tSch[ihit];
  }
  event.nhLac = src.nhLac;
  for(int ihit=0;ihit<event.nhLac;++ihit){
    event.tidLac[ihit]=src.tidLac[ihit];
    event.pidLac[ihit]=src.pidLac[ihit];
    event.didLac[ihit]=src.didLac[ihit];
    event.prtLac[ihit]=src.prtLac[ihit];
    event.qLac[ihit]=src.qLac[ihit];
    event.xLac[ihit]=src.xLac[ihit];
    event.yLac[ihit]=src.yLac[ihit];
    event.zLac[ihit]=src.zLac[ihit];
    event.pxLac[ihit]=src.pxLac[ihit];
    event.pyLac[ihit]=src.pyLac[ihit];
    event.pzLac[ihit]=src.pzLac[ihit];
    event.ppLac[ihit]=src.ppLac[ihit];
    event.deLac[ihit]=src.deLac[ihit];
    event.tLac[ihit]=src.tLac[ihit];
  }
  event.nhFtof = src.nhFtof;
  for(int ihit=0;ihit<event.nhFtof;++ihit){
    event.tidFtof[ihit]=src.tidFtof[ihit];
    event.pidFtof[ihit]=src.pidFtof[ihit];
    event.didFtof[ihit]=src.didFtof[ihit];
    event.prtFtof[ihit]=src.prtFtof[ihit];
    event.qFtof[ihit]=src.qFtof[ihit];
    event.xFtof[ihit]=src.xFtof[ihit];
    event.yFtof[ihit]=src.yFtof[ihit];
    event.zFtof[ihit]=src.zFtof[ihit];
    event.pxFtof[ihit]=src.pxFtof[ihit];
    event.pyFtof[ihit]=src.pyFtof[ihit];
    event.pzFtof[ihit]=src.pzFtof[ihit];
    event.ppFtof[ihit]=src.ppFtof[ihit];
    event.deFtof[ihit]=src.deFtof[ihit];
    event.tFtof[ihit]=src.tFtof[ihit];
  }
  event.nhSdc = src.nhSdc;
  for(int ihit=0;ihit<event.nhSdc;++ihit){
    event.tidSdc[ihit]=src.tidSdc[ihit];
    event.pidSdc[ihit]=src.pidSdc[ihit];
    event.didSdc[ihit]=src.didSdc[ihit];
    event.prtSdc[ihit]=src.prtSdc[ihit];
    event.qSdc[ihit]=src.qSdc[ihit];
    event.xSdc[ihit]=src.xSdc[ihit];
    event.ySdc[ihit]=src.ySdc[ihit];
    event.zSdc[ihit]=src.zSdc[ihit];
    event.pxSdc[ihit]=src.pxSdc[ihit];
    event.pySdc[ihit]=src.pySdc[ihit];
    event.pzSdc[ihit]=src.pzSdc[ihit];
    event.ppSdc[ihit]=src.ppSdc[ihit];
    event.deSdc[ihit]=src.deSdc[ihit];
    event.tSdc[ihit]=src.tSdc[ihit];
  }

  event.nhBvh = src.nhBvh;
  for(int ihit=0;ihit<event.nhBvh;++ihit){
    event.tidBvh[ihit]=src.tidBvh[ihit];
    event.pidBvh[ihit]=src.pidBvh[ihit];
    event.didBvh[ihit]=src.didBvh[ihit];
    event.prtBvh[ihit]=src.prtBvh[ihit];
    event.qBvh[ihit]=src.qBvh[ihit];
    event.xBvh[ihit]=src.xBvh[ihit];
    event.yBvh[ihit]=src.yBvh[ihit];
    event.zBvh[ihit]=src.zBvh[ihit];
    event.pxBvh[ihit]=src.pxBvh[ihit];
    event.pyBvh[ihit]=src.pyBvh[ihit];
    event.pzBvh[ihit]=src.pzBvh[ihit];
    event.ppBvh[ihit]=src.ppBvh[ihit];
    event.deBvh[ihit]=src.deBvh[ihit];
    event.tBvh[ihit]=src.tBvh[ihit];
  }

  for(int ihit=0; ihit<event.nhittpc; ++ihit){
    event.ititpc[ihit] = src.ititpc[ihit];

    //for debug
    //debug    std::cout<<"iti:"<<src.ititpc[ihit]<<", pos:"
    //debug     <<"("<<src.xtpc[ihit]<<", "
    //debug     <<src.ytpc[ihit]<<", "
    //debug     <<src.ztpc[ihit]<<")"<<std::endl;

    if(event.max_ititpc<src.ititpc[ihit])
      event.max_ititpc = src.ititpc[ihit];
    ++event.nhittpc_iti[src.ititpc[ihit]];
  }
  if(src.nhittpc<5)
    return true;

  TPCAnalyzer *TPCAna = new TPCAnalyzer();
  if(IsWithRes){
    TPCAna->DecodeTPCHitsGeant4(src.nhittpc,
				src.xtpc, src.ytpc, src.ztpc, src.edeptpc, src.idtpc, G4Moms);
  }
  else{
    TPCAna->DecodeTPCHitsGeant4(src.nhittpc,
				src.x0tpc, src.y0tpc, src.z0tpc, src.edeptpc, src.idtpc, G4Moms);
  }
  event.xtgtHS = **src.xtgtHS;
  event.ytgtHS = **src.ytgtHS;
  event.xtgtKurama = **src.xtgtKurama;
  event.ytgtKurama = **src.ytgtKurama;
#if KuramaK18
  event.ntK18 = **src.ntK18;
  event.xvpHS = **src.xvpHS;
  event.yvpHS = **src.yvpHS;
  event.zvpHS = **src.zvpHS;
  event.ztgtHS = **src.ztgtHS;
  event.p_3rd = **src.p_3rd;
  event.xoutK18 = **src.xoutK18;
  event.youtK18 = **src.youtK18;
  event.uoutK18 = **src.uoutK18;
  event.voutK18 = **src.voutK18;
  event.layerK18 = **src.layerK18;
  event.wireK18 = **src.wireK18;
  event.localhitposK18 = **src.localhitposK18;

  vector<vector<TVector3>>vpK18;
  vpK18.resize(event.ntK18);
  vector<TVector3> initPosK18;
  vector<TVector3> initMomK18;
  vector<vector<int>> intlayerK18;
  intlayerK18.resize(event.ntK18);
  for(int it=0; it<event.ntK18; ++it){
    for(Int_t il=0; il<event.xvpHS[it].size(); ++il){
      vpK18[it].push_back(TVector3(event.xvpHS[it][il], event.yvpHS[it][il], event.zvpHS[it][il]));
      intlayerK18[it].push_back(event.layerK18[it][il]);
    }
    vpK18[it].push_back(TVector3(event.xtgtHS[it], event.ytgtHS[it], event.ztgtHS[it]));
    TVector3 posOut(xGlobalBcOut + event.xoutK18[it],
		    yGlobalBcOut + event.youtK18[it],
		    zGlobalBcOut - zLocalBcOut);
    TVector3 momOut(event.uoutK18[it], event.voutK18[it], 1.);
    momOut *= 1./momOut.Mag();
    momOut *= event.p_3rd[it];

    initPosK18.push_back(posOut);
    initMomK18.push_back(momOut);

  }
  event.ntKurama = **src.ntKurama;
  event.xvpKurama = **src.xvpKurama;
  event.yvpKurama = **src.yvpKurama;
  event.zvpKurama = **src.zvpKurama;
  event.xout = **src.xout;
  event.yout = **src.yout;
  event.zout = **src.zout;
  event.pxout = **src.pxout;
  event.pyout = **src.pyout;
  event.pzout = **src.pzout;
  event.layer = **src.layer;
  event.wire = **src.wire;
  event.localhitpos = **src.localhitpos;
  vector<Int_t> pidKurama;
  vector<TVector3> initPosKurama;
  vector<TVector3> initMomKurama;
  vector<vector<TVector3>> vpKurama;
  vector<vector<int>> intlayer;
  intlayer.resize(event.ntKurama);
  vpKurama.resize( event.ntKurama );
  for(Int_t it=0; it<event.ntKurama; ++it){
    vpKurama[it].push_back(TVector3(event.xtgtKurama[it], event.ytgtKurama[it], tpc::ZTarget));
    cout<<event.ntKurama<<", "<<event.xvpKurama.size()<<endl;
    for(Int_t il=0; il<event.xvpKurama[it].size(); ++il){
      vpKurama[it].push_back(TVector3(event.xvpKurama[it][il], event.yvpKurama[it][il], event.zvpKurama[it][il]));
      intlayer[it].push_back(event.layer[it][il]);
    }
    pidKurama.push_back(1);
    TVector3 posOut(xGlobalSdcOut + event.xout[it], yGlobalSdcOut + event.yout[it], event.zout[it]);
    TVector3 momOut(event.pxout[it], event.pyout[it], event.pzout[it]);
    initPosKurama.push_back(posOut);
    initMomKurama.push_back(momOut);
  }
#endif
#if KuramaK18
#if ExclusiveResidual
  TPCAna->TrackSearchTPCHelix(vpK18, vpKurama,true);
#else
  TPCAna->TrackSearchTPCHelix(vpK18, vpKurama);
#endif
  TPCAna->TrackSearchTPCKurama(pidKurama, initPosKurama, initMomKurama, intlayer, event.wire, event.localhitpos);
  TPCAna->TrackSearchTPCK18(initPosK18, initMomK18, intlayerK18, event.wireK18, event.localhitposK18);
#else
#if ExclusiveResidual
  TPCAna->TrackSearchTPCHelix(true);
#else
  TPCAna->TrackSearchTPCHelix();
#endif
#endif
  Int_t nclTpc = 0;
  Int_t remain_nclTpc = 0;
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    auto hc = TPCAna->GetTPCClCont( layer );
    int ic=0;
    for( const auto& cl : hc ){
      ic++;
      if( !cl || !cl->IsGood() )
        continue;
      Double_t x = cl->GetX();
      Double_t y = cl->GetY();
      Double_t z = cl->GetZ();
      int pad = tpc::findPadID(z,x);
      int layer = tpc::getLayerID(pad);
      int row = tpc::getRowID(pad);
      double val = 0;
      gTPC.GetCDe(layer,row,1,val);
      if(val == 0){
        continue;
      }
      Double_t de = cl->GetDe();
      Int_t cl_size = cl->GetClusterSizeG4();
      Double_t mrow = cl->MeanRow();
      TPCHit* meanHit = cl->GetMeanHit();
      Int_t houghflag = meanHit->GetHoughFlag();
      /*
	TPCHit* centerHit = cl->GetCenterHit();
	const TVector3& centerPos = centerHit->GetPosition();
	Double_t centerDe = centerHit->GetCDe();
	Int_t centerRow = centerHit->GetRow();
      */
      TVector3 hitpos(x,y,z);
      event.cluster_x.push_back(x);
      event.cluster_y.push_back(y);
      event.cluster_z.push_back(z);
      event.cluster_de.push_back(de);
      event.cluster_size.push_back(cl_size);
      event.cluster_layer.push_back(layer);
      event.cluster_mrow.push_back(mrow);
      /*
	event.cluster_de_center.push_back(centerDe);
	event.cluster_x_center.push_back(centerPos.X());
	event.cluster_y_center.push_back(centerPos.Y());
	event.cluster_z_center.push_back(centerPos.Z());
	event.cluster_row_center.push_back(centerRow);
      */
      event.cluster_houghflag.push_back(houghflag);
      event.cluster_G4tid.push_back(GetHitG4Tid(hitpos,src.nhittpc,src.ititpc,src.xtpc,src.ytpc,src.ztpc));
      ++nclTpc;

      if(houghflag!=100&&houghflag!=200) ++remain_nclTpc; //Clusters without track
    }
    event.nclTpc = nclTpc;
    event.remain_nclTpc = remain_nclTpc;
  }


  int ntTpc = TPCAna->GetNTracksTPCHelix();
  if( 1000<ntTpc ){
    std::cout << "#W " << func_name << " "
      	      << "too many ntTpc " << ntTpc << "/" << 1000 << std::endl;
    ntTpc = 1000;
  }
  event.ntTpc = ntTpc;
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
  event.purity.resize( ntTpc );
  event.efficiency.resize( ntTpc );
  event.G4tid.resize( ntTpc );
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

  event.dE.resize( ntTpc );
  event.dEdx.resize( ntTpc );
  event.dz_factor.resize( ntTpc );
  event.mom0.resize( ntTpc );
  event.charge.resize( ntTpc );
  event.path.resize( ntTpc );
  event.pid.resize( ntTpc );
  event.isElectron.resize(ntTpc);
  event.nsigma_triton.resize(ntTpc);
  event.nsigma_deutron.resize(ntTpc);
  event.nsigma_proton.resize(ntTpc);
  event.nsigma_kaon.resize(ntTpc);
  event.nsigma_pion.resize(ntTpc);
  event.nsigma_electron.resize(ntTpc);

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
  event.resolution.resize( ntTpc );
  event.resolution_x.resize( ntTpc );
  event.resolution_y.resize( ntTpc );
  event.resolution_z.resize( ntTpc );
  event.residual.resize( ntTpc );
  event.residual_t.resize( ntTpc );
  event.residual_x.resize( ntTpc );
  event.residual_y.resize( ntTpc );
  event.residual_z.resize( ntTpc );
  event.pull.resize( ntTpc );
  event.pull_t.resize( ntTpc );
  event.pull_x.resize( ntTpc );
  event.pull_y.resize( ntTpc );
  event.pull_z.resize( ntTpc );
  event.helix_t.resize( ntTpc );
  event.pathhit.resize(ntTpc);
  event.alpha.resize(ntTpc);
  event.houghflag.resize(ntTpc);
  event.track_cluster_size.resize(ntTpc);
  event.track_cluster_de.resize(ntTpc);
  event.track_cluster_mrow.resize(ntTpc);

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


  event.isgoodTPCK18.resize(event.ntK18);
  event.tpcidTPCK18.resize(event.ntK18);
  event.niterationTPCK18.resize(event.ntK18);
  event.chisqrTPCK18.resize(event.ntK18);
  for(Int_t itk18=0; itk18<TPCAna->GetNTracksTPCK18(); ++itk18){
    TPCRKTrack* tr_km = TPCAna->GetTPCK18Track(itk18);
    Int_t niteration = tr_km -> Niteration();
    Double_t chisqr = tr_km -> GetChiSquare();
    Int_t idtpc = tr_km -> GetTPCTrackID();
    Int_t idk18 = tr_km -> GetTrackID();
    event.isgoodTPCK18[idk18] = 1;
    event.tpcidTPCK18[idk18] = idtpc;
    event.niterationTPCK18[idk18] = niteration;
    event.chisqrTPCK18[idk18] = chisqr;
  }
  event.isgoodTPCKurama.resize( event.ntKurama );
  event.tpcidTPCKurama.resize( event.ntKurama );
  event.niterationTPCKurama.resize( event.ntKurama );
  event.kflagTPCKurama.resize( event.ntKurama );
  event.chisqrTPCKurama.resize( event.ntKurama );
  for(Int_t ittpckurama=0; ittpckurama<TPCAna->GetNTracksTPCKurama(); ++ittpckurama){
    TPCRKTrack* tr_kp = TPCAna->GetTPCKuramaTrack(ittpckurama);
    Int_t niteration = tr_kp -> Niteration();
    Double_t chisqr = tr_kp -> GetChiSquare();
    Int_t idtpc = tr_kp -> GetTPCTrackID();
    Int_t idkurama = tr_kp -> GetTrackID();
    event.isgoodTPCKurama[idkurama] = 1;
    event.tpcidTPCKurama[idkurama] = idtpc;
    event.niterationTPCKurama[idkurama] = niteration;
    event.chisqrTPCKurama[idkurama] = chisqr;
    event.kflagTPCKurama[idkurama] = 1;

  }

  vector<int> G4TrackID;
  vector<int> PureHits;

  for( int it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp= TPCAna->GetTrackTPCHelix(it);
    if(!tp) continue;
    int nh=tp->GetNHit();
    int nhEff=tp->GetNHitsEffective();
    double chisqr    = tp->GetChiSquare();
    Double_t pval = 1-ROOT::Math::chisquared_cdf( chisqr*(2*nhEff-5),2*nhEff-5);
    double helix_cx=tp->Getcx(), helix_cy=tp->Getcy();
    double helix_z0=tp->Getz0(), helix_r=tp->Getr();
    double helix_dz = tp->Getdz();
    TVector3 Mom0 = tp->GetMom0();
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
    {
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
      event.charge[it] = charge;
      event.path[it] = pathlen;
      event.distTgt[it] = distTgt;
      event.chisqr[it] = chisqr;
      event.pval[it] = pval;
      event.helix_cx[it] = helix_cx;
      event.helix_cy[it] = helix_cy;
      event.helix_z0[it] = helix_z0;
      event.helix_r[it]  = helix_r;
      event.helix_dz[it] = helix_dz;
      event.mom0_x[it] = Mom0.X();
      event.mom0_y[it] = Mom0.Y();
      event.mom0_z[it] = Mom0.Z();
      event.mom0[it] = Mom0.Mag();
      event.pid[it] = pid;
      event.dE[it] = tp->GetTrackdE();
      event.dEdx[it] = tp->GetdEdx(truncatedMean);
      event.dz_factor[it] = sqrt(1.+(pow(helix_dz,2)));
      event.isElectron[it] = Kinematics::HypTPCdEdxElectron(event.dEdx[it], event.mom0[it]);
      event.nsigma_triton[it] = Kinematics::HypTPCdEdxNsigmaTriton(event.dEdx[it], event.mom0[it]);
      event.nsigma_deutron[it] = Kinematics::HypTPCdEdxNsigmaDeutron(event.dEdx[it], event.mom0[it]);
      event.nsigma_proton[it] = Kinematics::HypTPCdEdxNsigmaProton(event.dEdx[it], event.mom0[it]);
      event.nsigma_kaon[it]  = Kinematics::HypTPCdEdxNsigmaKaon(event.dEdx[it], event.mom0[it]);
      event.nsigma_pion[it] = Kinematics::HypTPCdEdxNsigmaPion(event.dEdx[it], event.mom0[it]);
      event.nsigma_electron[it] = Kinematics::HypTPCdEdxNsigmaElectron(event.dEdx[it], event.mom0[it]);

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
      event.pull[it].resize( nh );
      event.pull_t[it].resize( nh );
      event.pull_x[it].resize( nh );
      event.pull_y[it].resize( nh );
      event.pull_z[it].resize( nh );
      event.resolution_x[it].resize( nh );
      event.resolution_y[it].resize( nh );
      event.resolution_z[it].resize( nh );
      event.helix_t[it].resize( nh );
      event.pathhit[it].resize(nh);
      event.alpha[it].resize(nh);
      event.houghflag[it].resize(nh);
      event.track_cluster_size[it].resize(nh);
      event.track_cluster_de[it].resize(nh);
      event.track_cluster_mrow[it].resize(nh);
    }

    //debug     std::cout<<"nh:"<<nh<<std::endl;
    double min_t = 9999,max_t = -9999;
    vector<TVector3> TPCHits;
    for( int ih=0; ih<nh; ++ih ){
      TPCLTrackHit *hit = tp -> GetHitInOrder(ih);
      if(!hit) continue;
      TPCHits.push_back(hit->GetLocalHitPos());
      int order = tp->GetOrder(ih);
      int layerId = layerId = hit->GetLayer();
      Int_t houghflag = hit->GetHoughFlag();
      const TVector3& resi_vect = hit->GetResidualVect();
      const TVector3& res_vect = hit->GetResolutionVect();
      const TVector3& hitpos = hit->GetLocalHitPos();
      TPCHit *clhit = hit->GetHit();
      TPCCluster *cl = clhit->GetParentCluster();
      Double_t clsize = cl->GetClusterSizeG4();
      Double_t clde = hit->GetDe();
      Double_t mrow = hit->GetMRow();
      const TVector3& calpos = hit->GetLocalCalPosHelix();
      const TVector3& mom = hit->GetMomentumHelix(tp->GetCharge());
      double residual = hit->GetResidual();
      const TVector3& excl_resi_vect = hit->GetResidualVectExclusive();
#if ExclusiveResidual
      const TVector3& ex_calpos	= hit->GetLocalCalPosHelixExclusive();
      double excl_cx=tp->GetcxExclusive(order), excl_cy=tp->GetcyExclusive(order);
      double excl_z0=tp->Getz0Exclusive(order), excl_r=tp->GetrExclusive(order);
      double excl_dz = tp->GetdzExclusive(order);
      //      const TVector3& excl_resi_vect = hit->GetResidualVectExclusive();
#endif
      for( int ih2=0; ih2<src.nhittpc; ++ih2 ){
	TVector3 setpos;
	if(IsWithRes)
	  setpos = TVector3(src.xtpc[ih2], src.ytpc[ih2], src.ztpc[ih2]);
	else
	  setpos = TVector3(src.x0tpc[ih2], src.y0tpc[ih2], src.z0tpc[ih2]);

	TVector3 d = setpos - hitpos;
	if(fabs(d.Mag()<0.1)){
	  event.momg_x[it][ih] = src.pxtpc[ih2]*1000.;
	  event.momg_y[it][ih] = src.pytpc[ih2]*1000.;
	  event.momg_z[it][ih] = src.pztpc[ih2]*1000.;
	  double momg_mag = sqrt(src.pxtpc[ih2]*src.pxtpc[ih2]
				 +src.pytpc[ih2]*src.pytpc[ih2]
				 +src.pztpc[ih2]*src.pztpc[ih2])*1000.;
	  event.residual_px[it][ih] = mom.x() - src.pxtpc[ih2]*1000.;//MeV/c
	  event.residual_py[it][ih] = mom.y() - src.pytpc[ih2]*1000.;//MeV/c
	  event.residual_pz[it][ih] = mom.z() - src.pztpc[ih2]*1000.;//MeV/c
	  event.residual_p[it][ih] = mom.Mag() - momg_mag;//MeV/c

	  event.iti_g[it][ih] = src.ititpc[ih2];
	  break;
	}
      }
      event.track_cluster_size[it][ih] = clsize;
      event.track_cluster_de[it][ih] = clde;
      event.track_cluster_mrow[it][ih] = mrow;
      event.pathhit[it][ih] = hit->GetPathHelix();
      event.alpha[it][ih] = tp->GetAlpha(ih);

      event.hitlayer[it][ih] = (double)layerId;
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
      double res_t = hypot(res_vect.x(),res_vect.z());
      double res_x = res_vect.x();
      double res_y = res_vect.y();
      double res_z = res_vect.z();

      double resi_x = resi_vect.x();
      double resi_y = resi_vect.y();
      double resi_z = resi_vect.z();
      if(res_x>1e6 or res_y>1e6 or res_z > 1e6) continue;
      double resi_theta = atan2(resi_z,-resi_x);
      TVector3 dir_hit(cos(hit->GetTheta()),sin(hit->GetTheta()),0);
      TVector3 dir_resi(cos(resi_theta),sin(resi_theta),0);
      double sign = 1.;
      if(dir_hit*dir_resi<0) sign = -1.;
      //			double hitpos_r = hypot(hitpos.z()-(-143) -cy,-hitpos.x()-helix_cx);
      //			double resi_t = hitpos_r - r;
      double resi_t = sign*hypot(resi_x,resi_z);
      event.residual_t[it][ih] = resi_t;
      event.pull_t[it][ih] = resi_t/res_t;
      event.pull_x[it][ih] = resi_x/res_x;
      event.pull_y[it][ih] = resi_y/res_y;
      event.pull_z[it][ih] = resi_z/res_z;

      if(pval>0.01) continue;
      HF1(TPCInclusiveHid+layerId,resi_t);
      HF1(TPCInclusiveHid+100+layerId,resi_x);
      HF1(TPCInclusiveHid+200+layerId,resi_y);
      HF1(TPCInclusiveHid+300+layerId,resi_z);
      HF1(TPCInclusiveHid+1000+layerId,resi_t/res_t);
      HF1(TPCInclusiveHid+1100+layerId,resi_x/res_x);
      HF1(TPCInclusiveHid+1200+layerId,resi_y/res_y);
      HF1(TPCInclusiveHid+1300+layerId,resi_z/res_z);
      HF1(TPCInclusiveHid+32,resi_t);
      HF1(TPCInclusiveHid+100+32,resi_x);
      HF1(TPCInclusiveHid+200+32,resi_y);
      HF1(TPCInclusiveHid+300+32,resi_z);
      HF1(TPCInclusiveHid+1000+32,resi_t/res_t);
      HF1(TPCInclusiveHid+1100+32,resi_x/res_x);
      HF1(TPCInclusiveHid+1200+32,resi_y/res_y);
      HF1(TPCInclusiveHid+1300+32,resi_z/res_z);

      HF2(TPCInclusiveHid+10000,layerId,resi_t);
      HF2(TPCInclusiveHid+100+10000,layerId,resi_x);
      HF2(TPCInclusiveHid+200+10000,layerId,resi_y);
      HF2(TPCInclusiveHid+300+10000,layerId,resi_z);
      HF2(TPCInclusiveHid+1000+10000,layerId,resi_t/res_t);
      HF2(TPCInclusiveHid+1100+10000,layerId,resi_x/res_x);
      HF2(TPCInclusiveHid+1200+10000,layerId,resi_y/res_y);
      HF2(TPCInclusiveHid+1300+10000,layerId,resi_z/res_z);
      double alpha = abs(event.alpha[it][ih]);
      int bin_alpha = (50*alpha / acos(-1));
      HF1(TPCInclusiveHid + 20000 + bin_alpha,resi_t);
      HF1(TPCInclusiveHid +100+ 20000 + bin_alpha,resi_x);
      HF1(TPCInclusiveHid +300+ 20000 + bin_alpha,resi_z);
#if ExclusiveResidual
      //			double excl_resi_t = excl_hitpos_r - excl_r;
      double excl_resi_x = excl_resi_vect.x();
      double excl_resi_y = excl_resi_vect.y();
      double excl_resi_z = excl_resi_vect.z();
      double excl_resi_t = sign*hypot(excl_resi_x,excl_resi_z);

      HF1(TPCExclusiveHid+layerId,excl_resi_t);
      HF1(TPCExclusiveHid+100+layerId,excl_resi_x);
      HF1(TPCExclusiveHid+200+layerId,excl_resi_y);
      HF1(TPCExclusiveHid+300+layerId,excl_resi_z);
      HF1(TPCExclusiveHid+1000+layerId,excl_resi_t/res_t);
      HF1(TPCExclusiveHid+1100+layerId,excl_resi_x/res_x);
      HF1(TPCExclusiveHid+1200+layerId,excl_resi_y/res_y);
      HF1(TPCExclusiveHid+1300+layerId,excl_resi_z/res_z);
      HF1(TPCExclusiveHid+32,excl_resi_t);
      HF1(TPCExclusiveHid+100+32,excl_resi_x);
      HF1(TPCExclusiveHid+200+32,excl_resi_y);
      HF1(TPCExclusiveHid+300+32,excl_resi_z);
      HF1(TPCExclusiveHid+1000+32,excl_resi_t/res_t);
      HF1(TPCExclusiveHid+1100+32,excl_resi_x/res_x);
      HF1(TPCExclusiveHid+1200+32,excl_resi_y/res_y);
      HF1(TPCExclusiveHid+1300+32,excl_resi_z/res_z);

      HF2(TPCExclusiveHid+10000,layerId,excl_resi_t);
      HF2(TPCExclusiveHid+100+10000,layerId,excl_resi_x);
      HF2(TPCExclusiveHid+200+10000,layerId,excl_resi_y);
      HF2(TPCExclusiveHid+300+10000,layerId,excl_resi_z);
      HF2(TPCExclusiveHid+1000+10000,layerId,excl_resi_t/res_t);
      HF2(TPCExclusiveHid+1100+10000,layerId,excl_resi_x/res_x);
      HF2(TPCExclusiveHid+1200+10000,layerId,excl_resi_y/res_y);
      HF2(TPCExclusiveHid+1300+10000,layerId,excl_resi_z/res_z);

      double intr_resi_t = sqrt(abs(resi_t*excl_resi_t));
      double intr_resi_x = sqrt(abs(resi_x*excl_resi_x));
      double intr_resi_y = sqrt(abs(resi_y*excl_resi_y));
      double intr_resi_z = sqrt(abs(resi_z*excl_resi_z));
      if(resi_t<0) intr_resi_t*=-1;
      if(resi_x<0) intr_resi_x*=-1;
      if(resi_y<0) intr_resi_y*=-1;
      if(resi_z<0) intr_resi_z*=-1;

      HF1(TPCIntrinsicHid+layerId,intr_resi_t);
      HF1(TPCIntrinsicHid+100+layerId,intr_resi_x);
      HF1(TPCIntrinsicHid+200+layerId,intr_resi_y);
      HF1(TPCIntrinsicHid+300+layerId,intr_resi_z);
      HF1(TPCIntrinsicHid+1000+layerId,intr_resi_t/res_t);
      HF1(TPCIntrinsicHid+1100+layerId,intr_resi_x/res_x);
      HF1(TPCIntrinsicHid+1200+layerId,intr_resi_y/res_y);
      HF1(TPCIntrinsicHid+1300+layerId,intr_resi_z/res_z);
      HF1(TPCIntrinsicHid+32,intr_resi_t);
      HF1(TPCIntrinsicHid+100+32,intr_resi_x);
      HF1(TPCIntrinsicHid+200+32,intr_resi_y);
      HF1(TPCIntrinsicHid+300+32,intr_resi_z);
      HF1(TPCIntrinsicHid+1000+32,intr_resi_t/res_t);
      HF1(TPCIntrinsicHid+1100+32,intr_resi_x/res_x);
      HF1(TPCIntrinsicHid+1200+32,intr_resi_y/res_y);
      HF1(TPCIntrinsicHid+1300+32,intr_resi_z/res_z);

      HF2(TPCIntrinsicHid+10000,layerId,intr_resi_t);
      HF2(TPCIntrinsicHid+100+10000,layerId,intr_resi_x);
      HF2(TPCIntrinsicHid+200+10000,layerId,intr_resi_y);
      HF2(TPCIntrinsicHid+300+10000,layerId,intr_resi_z);
      HF2(TPCIntrinsicHid+1000+10000,layerId,intr_resi_t/res_t);
      HF2(TPCIntrinsicHid+1100+10000,layerId,intr_resi_x/res_x);
      HF2(TPCIntrinsicHid+1200+10000,layerId,intr_resi_y/res_y);
      HF2(TPCIntrinsicHid+1300+10000,layerId,intr_resi_z/res_z);

      HF1(TPCExclusiveHid + 20000 + bin_alpha,resi_t);
      HF1(TPCExclusiveHid +100+ 20000 + bin_alpha,resi_x);
      HF1(TPCExclusiveHid +300+ 20000 + bin_alpha,resi_z);
      HF1(TPCIntrinsicHid + 20000 + bin_alpha,resi_t);
      HF1(TPCIntrinsicHid +100+ 20000 + bin_alpha,resi_x);
      HF1(TPCIntrinsicHid +300+ 20000 + bin_alpha,resi_z);

#endif
    }//ih
    int nPureHits;
    int G4tid = TPCToG4TrackID(TPCHits,src.nhittpc,src.ititpc,src.xtpc,src.ytpc,src.ztpc,nPureHits);
    if(G4tid == G4idKm)event.isK18[it]=1;
    if(G4tid == G4idKp)event.isKurama[it]=1;

    event.G4tid[it] = G4tid;
    event.purity[it] = (double)nPureHits/nh;
    int nG4Hits = event.nhittpc_iti[G4tid];
    event.efficiency[it] = (double)nPureHits/nG4Hits;

    //Inverted charge tracks
    TPCLocalTrackHelix *tp_inverted = TPCAna->GetTrackTPCHelixChargeInverted( it );
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
  }//it

  Int_t nvtxTpc = TPCAna->GetNVerticesTPC();
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
    TPCVertex *vp = TPCAna->GetTPCVertex( it );
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

  Int_t nvtxTpcClustered = TPCAna->GetNVerticesTPCClustered();
  event.nvtxTpcClustered = nvtxTpcClustered;
  event.Clusteredvtx_x.resize(nvtxTpcClustered);
  event.Clusteredvtx_y.resize(nvtxTpcClustered);
  event.Clusteredvtx_z.resize(nvtxTpcClustered);
  event.Clusteredvtxid.resize(nvtxTpcClustered);
  for( Int_t ivtx=0; ivtx<nvtxTpcClustered; ++ivtx ){
    TPCVertex *vp = TPCAna->GetTPCVertexClustered( ivtx );
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

  delete TPCAna;
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
    if( TTreeReaderCont[i] ) delete TTreeReaderCont[i];
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
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    HB1(TPCInclusiveHid+layer,	Form("TPC InclusiveResidual T[mm];layer%d",layer),1000,-10,10);
    HB1(TPCInclusiveHid+100+layer,	Form("TPC InclusiveResidual X[mm];layer%d",layer),1000,-10,10);
    HB1(TPCInclusiveHid+200+layer,	Form("TPC InclusiveResidual Y[mm];layer%d",layer),1000,-10,10);
    HB1(TPCInclusiveHid+300+layer,	Form("TPC InclusiveResidual Z[mm];layer%d",layer),1000,-10,10);
    HB1(TPCInclusiveHid+1000+layer,	Form("TPC InclusivePull T;layer%d",layer),1000,-10,10);
    HB1(TPCInclusiveHid+1000+100+layer,	Form("TPC InclusivePull X;layer%d",layer),1000,-10,10);
    HB1(TPCInclusiveHid+1000+200+layer,	Form("TPC InclusivePull Y;layer%d",layer),1000,-10,10);
    HB1(TPCInclusiveHid+1000+300+layer,	Form("TPC InclusivePull Z;layer%d",layer),1000,-10,10);

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

  HB1(TPCIntrinsicHid+32,	Form("TPC IntrinsicResidual T[mm];AllLayer"),1000,-10,10);
  HB1(TPCIntrinsicHid+100+32,	Form("TPC IntrinsicResidual X[mm];AllLayer"),1000,-10,10);
  HB1(TPCIntrinsicHid+200+32,	Form("TPC IntrinsicResidual Y[mm];AllLayer"),1000,-10,10);
  HB1(TPCIntrinsicHid+300+32,	Form("TPC IntrinsicResidual Z[mm];AllLayer"),1000,-10,10);
  HB1(TPCIntrinsicHid+1000+32,	Form("TPC IntrinsicPull T;AllLayer"),1000,-10,10);
  HB1(TPCIntrinsicHid+1000+100+32,	Form("TPC IntrinsicPull X;AllLayer"),1000,-10,10);
  HB1(TPCIntrinsicHid+1000+200+32,	Form("TPC IntrinsicPull Y;AllLayer"),1000,-10,10);
  HB1(TPCIntrinsicHid+1000+300+32,	Form("TPC IntrinsicPull Z;AllLayer"),1000,-10,10);


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

  for(int ih=0;ih<50;++ih){
    double dt = acos(-1)*1./50;
    double alpha_l = dt*ih;
    double alpha_h = dt*(ih+1);
    HB1(TPCInclusiveHid + 20000 + ih,Form("TPCInclusiveResidualPt_#alpha=[%g,%g)",alpha_l,alpha_h),1000,-10,10);
    HB1(TPCInclusiveHid +100+ 20000 + ih,Form("TPCInclusiveResidualPx_#alpha=[%g,%g)",alpha_l,alpha_h),1000,-10,10);
    HB1(TPCInclusiveHid +300+ 20000 + ih,Form("TPCInclusiveResidualPz_#alpha=[%g,%g)",alpha_l,alpha_h),1000,-10,10);
    HB1(TPCExclusiveHid + 20000 + ih,Form("TPCExclusiveResidualPt_#alpha=[%g,%g)",alpha_l,alpha_h),1000,-10,10);
    HB1(TPCExclusiveHid +100+ 20000 + ih,Form("TPCExclusiveResidualPx_#alpha=[%g,%g)",alpha_l,alpha_h),1000,-10,10);
    HB1(TPCExclusiveHid +300+ 20000 + ih,Form("TPCExclusiveResidualPz_#alpha=[%g,%g)",alpha_l,alpha_h),1000,-10,10);
    HB1(TPCIntrinsicHid + 20000 + ih,Form("TPCIntrinsicResidualPt_#alpha=[%g,%g)",alpha_l,alpha_h),1000,-10,10);
    HB1(TPCIntrinsicHid +100+ 20000 + ih,Form("TPCIntrinsicResidualPx_#alpha=[%g,%g)",alpha_l,alpha_h),1000,-10,10);
    HB1(TPCIntrinsicHid +300+ 20000 + ih,Form("TPCIntrinsicResidualPz_#alpha=[%g,%g)",alpha_l,alpha_h),1000,-10,10);
  }

  HBTree( "tpc", "tree of DstTPC_g" );

  tree->Branch("status", &event.status, "status/I" );
  tree->Branch("evnum", &event.evnum, "evnum/I" );
  tree->Branch("nhittpc",&event.nhittpc,"nhittpc/I");

  tree->Branch("max_ititpc",&event.max_ititpc,"max_ititpc/I");
  tree->Branch("ititpc",event.ititpc,"ititpc[nhittpc]/I");
  tree->Branch("nhittpc_iti",event.nhittpc_iti,"nhittpc_iti[max_ititpc]/I");

  tree->Branch("xtgtHS",&event.xtgtHS);
  tree->Branch("ytgtHS",&event.ytgtHS);
  tree->Branch("xtgtKurama",&event.xtgtKurama);
  tree->Branch("ytgtKurama",&event.ytgtKurama);

  tree->Branch( "nclTpc", &event.nclTpc );
  tree->Branch( "remain_nclTpc", &event.remain_nclTpc );
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
  tree->Branch( "cluster_G4tid", &event.cluster_G4tid );
  tree->Branch( "cluster_houghflag", &event.cluster_houghflag );

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "nhtrackEff", &event.nhtrackEff );
  tree->Branch( "trackid", &event.trackid );
  tree->Branch( "isBeam", &event.isBeam );
  tree->Branch( "isXi", &event.isXi );
  tree->Branch( "isKurama", &event.isKurama );
  tree->Branch( "isK18", &event.isK18 );
  tree->Branch( "isAccidental", &event.isAccidental );
  tree->Branch( "isMultiloop", &event.isMultiloop );
  tree->Branch( "purity", &event.purity );
  tree->Branch( "efficiency", &event.efficiency );
  tree->Branch( "G4tid", &event.G4tid );
  tree->Branch( "flag", &event.flag );
  tree->Branch( "chisqr", &event.chisqr );
  tree->Branch( "pval", &event.pval );
  tree->Branch( "helix_cx", &event.helix_cx );
  tree->Branch( "helix_cy", &event.helix_cy );
  tree->Branch( "helix_z0", &event.helix_z0 );
  tree->Branch( "helix_r", &event.helix_r );
  tree->Branch( "helix_dz", &event.helix_dz );
  // Momentum at Y = 0
  tree->Branch( "mom0_x", &event.mom0_x );
  tree->Branch( "mom0_y", &event.mom0_y );
  tree->Branch( "mom0_z", &event.mom0_z );
  tree->Branch( "mom0", &event.mom0 );
  tree->Branch( "dE", &event.dE );
  tree->Branch( "dEdx", &event.dEdx );

  tree->Branch( "dz_factor", &event.dz_factor );
  tree->Branch( "charge", &event.charge );
  tree->Branch( "path", &event.path );
  tree->Branch( "isElectron", &event.isElectron );
  tree->Branch( "nsigma_triton", &event.nsigma_triton );
  tree->Branch( "nsigma_deutron", &event.nsigma_deutron );
  tree->Branch( "nsigma_proton", &event.nsigma_proton );
  tree->Branch( "nsigma_kaon", &event.nsigma_kaon );
  tree->Branch( "nsigma_pion", &event.nsigma_pion );
  tree->Branch( "nsigma_electron", &event.nsigma_electron );

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

  tree->Branch("momg_x",event.momg_x,"momg_x[nttpc][64]/D");
  tree->Branch("momg_y",event.momg_y,"momg_y[nttpc][64]/D");
  tree->Branch("momg_z",event.momg_z,"momg_z[nttpc][64]/D");
  tree->Branch("iti_g",event.iti_g,"iti_g[nttpc][64]/I");
  tree->Branch("residual_p",event.residual_p,"residual_p[nttpc][64]/D");

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

  tree->Branch("nPrm",&src.nhPrm,"nPrm/I");
  tree->Branch("xPrm",src.xPrm,"xPrm[nPrm]/D");
  tree->Branch("yPrm",src.yPrm,"yPrm[nPrm]/D");
  tree->Branch("zPrm",src.zPrm,"zPrm[nPrm]/D");
  tree->Branch("pxPrm",src.pxPrm,"pxPrm[nPrm]/D");
  tree->Branch("pyPrm",src.pyPrm,"pyPrm[nPrm]/D");
  tree->Branch("pzPrm",src.pzPrm,"pzPrm[nPrm]/D");
  tree->Branch("ppPrm",src.ppPrm,"ppPrm[nPrm]/D");

  tree->Branch("xtpc",src.xtpc,"xtpc[nhittpc]/D");
  tree->Branch("ytpc",src.ytpc,"ytpc[nhittpc]/D");
  tree->Branch("ztpc",src.ztpc,"ztpc[nhittpc]/D");

  tree->Branch("x0tpc",src.x0tpc,"x0tpc[nhittpc]/D");
  tree->Branch("y0tpc",src.y0tpc,"y0tpc[nhittpc]/D");
  tree->Branch("z0tpc",src.z0tpc,"z0tpc[nhittpc]/D");

  tree->Branch("pxtpc",src.pxtpc,"pxtpc[nhittpc]/D");
  tree->Branch("pytpc",src.pytpc,"pytpc[nhittpc]/D");
  tree->Branch("pztpc",src.pztpc,"pztpc[nhittpc]/D");
  tree->Branch("pptpc",src.pptpc,"pptpc[nhittpc]/D");
  tree->Branch("idtpc",src.idtpc,"idtpc[nhittpc]/I");
  tree->Branch("parentid",src.parentID,"parentID[nhittpc]/I");

  tree->Branch("NumberOfTracks",&event.NumberOfTracks,"NumberOfTracks/I");
  tree->Branch("PIDOfTrack",event.PIDOfTrack,"PIDOfTrack[1000]/I");
  tree->Branch("ParentIDOfTrack",event.ParentIDOfTrack,"ParentIDOfTrack[1000]/I");
  tree->Branch("VertexOfTrack_x",event.VertexOfTrack_x,"VertexOfTrack_x[1000]/D");
  tree->Branch("VertexOfTrack_y",event.VertexOfTrack_y,"VertexOfTrack_y[1000]/D");
  tree->Branch("VertexOfTrack_z",event.VertexOfTrack_z,"VertexOfTrack_z[1000]/D");
  tree->Branch("MomentumOfTrack",event.MomentumOfTrack,"MomentumOfTrack[1000]/D");
  tree->Branch("MomentumOfTrack_x",event.MomentumOfTrack_x,"MomentumOfTrack_x[1000]/D");
  tree->Branch("MomentumOfTrack_y",event.MomentumOfTrack_y,"MomentumOfTrack_y[1000]/D");
  tree->Branch("MomentumOfTrack_z",event.MomentumOfTrack_z,"MomentumOfTrack_z[1000]/D");
  tree->Branch("G4TrackIDKm",event.G4idKm);
  tree->Branch("G4TrackIDKp",event.G4idKp);
  tree->Branch("G4TrackIDP",event.G4idP);
  tree->Branch("G4TrackIDPi1",event.G4idPi1);
  tree->Branch("G4TrackIDPi2",event.G4idPi2);

  tree->Branch("MomXi_x",&event.MomXi_x,"MomXi_x/D");
  tree->Branch("MomXi_y",&event.MomXi_y,"MomXi_y/D");
  tree->Branch("MomXi_z",&event.MomXi_z,"MomXi_z/D");
  tree->Branch("SpinXi_x",&event.SpinXi_x,"SpinXi_x/D");
  tree->Branch("SpinXi_y",&event.SpinXi_y,"SpinXi_y/D");
  tree->Branch("SpinXi_z",&event.SpinXi_z,"SpinXi_z/D");
  tree->Branch("ThXi_CM",&event.ThXi_CM,"ThXi_CM/D");

  tree->Branch("MomLd_x",&event.MomLd_x,"MomLd_x/D");
  tree->Branch("MomLd_y",&event.MomLd_y,"MomLd_y/D");
  tree->Branch("MomLd_z",&event.MomLd_z,"MomLd_z/D");
  tree->Branch("SpinLd_x",&event.SpinLd_x,"SpinLd_x/D");
  tree->Branch("SpinLd_y",&event.SpinLd_y,"SpinLd_y/D");
  tree->Branch("SpinLd_z",&event.SpinLd_z,"SpinLd_z/D");
  tree->Branch("ThLd_CM",&event.ThLd_CM,"ThLd_CM/D");

  tree->Branch( "nhHtof", &event.nhHtof );
  tree->Branch( "HtofSeg", &event.HtofSeg );
  tree->Branch( "tHtof", &event.tHtof );
  tree->Branch( "dtHtof", &event.dtHtof );
  tree->Branch( "deHtof", &event.deHtof );
  tree->Branch( "posHtof", &event.posHtof );
  tree->Branch( "G4tidHtof", &event.G4tidHtof );
  tree->Branch("nhSch", &event.nhSch,"nhSch/I");
  tree->Branch("tidSch", event.tidSch,"tidSch[500]/I");
  tree->Branch("pidSch", event.pidSch,"pidSch[500]/I");
  tree->Branch("didSch", event.didSch,"didSch[500]/I");
  tree->Branch("prtSch", event.prtSch,"prtSch[500]/I");
  tree->Branch("qSch", event.qSch,"qSch[500]/I");
  tree->Branch("xSch", event.xSch,"xSch[500]/D");
  tree->Branch("ySch", event.ySch,"ySch[500]/D");
  tree->Branch("zSch", event.zSch,"zSch[500]/D");
  tree->Branch("pxSch", event.pxSch,"pxSch[500]/D");
  tree->Branch("pySch", event.pySch,"pySch[500]/D");
  tree->Branch("pzSch", event.pzSch,"pzSch[500]/D");
  tree->Branch("ppSch", event.ppSch,"ppSch[500]/D");
  tree->Branch("deSch", event.deSch,"deSch[500]/D");
  tree->Branch("tSch", event.tSch,"tSch[500]/D");

  tree->Branch("nhLac", &event.nhLac,"nhLac/I");
  tree->Branch("tidLac", event.tidLac,"tidLac[500]/I");
  tree->Branch("pidLac", event.pidLac,"pidLac[500]/I");
  tree->Branch("didLac", event.didLac,"didLac[500]/I");
  tree->Branch("prtLac", event.prtLac,"prtLac[500]/I");
  tree->Branch("qLac", event.qLac,"qLac[500]/I");
  tree->Branch("xLac", event.xLac,"xLac[500]/D");
  tree->Branch("yLac", event.yLac,"yLac[500]/D");
  tree->Branch("zLac", event.zLac,"zLac[500]/D");
  tree->Branch("pxLac", event.pxLac,"pxLac[500]/D");
  tree->Branch("pyLac", event.pyLac,"pyLac[500]/D");
  tree->Branch("pzLac", event.pzLac,"pzLac[500]/D");
  tree->Branch("ppLac", event.ppLac,"ppLac[500]/D");
  tree->Branch("deLac", event.deLac,"deLac[500]/D");
  tree->Branch("tLac", event.tLac,"tLac[500]/D");

  tree->Branch("nhFtof", &event.nhFtof,"nhFtof/I");
  tree->Branch("tidFtof", event.tidFtof,"tidFtof[500]/I");
  tree->Branch("pidFtof", event.pidFtof,"pidFtof[500]/I");
  tree->Branch("didFtof", event.didFtof,"didFtof[500]/I");
  tree->Branch("prtFtof", event.prtFtof,"prtFtof[500]/I");
  tree->Branch("qFtof", event.qFtof,"qFtof[500]/I");
  tree->Branch("xFtof", event.xFtof,"xFtof[500]/D");
  tree->Branch("yFtof", event.yFtof,"yFtof[500]/D");
  tree->Branch("zFtof", event.zFtof,"zFtof[500]/D");
  tree->Branch("pxFtof", event.pxFtof,"pxFtof[500]/D");
  tree->Branch("pyFtof", event.pyFtof,"pyFtof[500]/D");
  tree->Branch("pzFtof", event.pzFtof,"pzFtof[500]/D");
  tree->Branch("ppFtof", event.ppFtof,"ppFtof[500]/D");
  tree->Branch("deFtof", event.deFtof,"deFtof[500]/D");
  tree->Branch("tFtof", event.tFtof,"tFtof[500]/D");

  tree->Branch("nhSdc", &event.nhSdc,"nhSdc/I");
  tree->Branch("tidSdc", event.tidSdc,"tidSdc[500]/I");
  tree->Branch("pidSdc", event.pidSdc,"pidSdc[500]/I");
  tree->Branch("didSdc", event.didSdc,"didSdc[500]/I");
  tree->Branch("prtSdc", event.prtSdc,"prtSdc[500]/I");
  tree->Branch("qSdc", event.qSdc,"qSdc[500]/I");
  tree->Branch("xSdc", event.xSdc,"xSdc[500]/D");
  tree->Branch("ySdc", event.ySdc,"ySdc[500]/D");
  tree->Branch("zSdc", event.zSdc,"zSdc[500]/D");
  tree->Branch("pxSdc", event.pxSdc,"pxSdc[500]/D");
  tree->Branch("pySdc", event.pySdc,"pySdc[500]/D");
  tree->Branch("pzSdc", event.pzSdc,"pzSdc[500]/D");
  tree->Branch("ppSdc", event.ppSdc,"ppSdc[500]/D");
  tree->Branch("deSdc", event.deSdc,"deSdc[500]/D");
  tree->Branch("tSdc", event.tSdc,"tSdc[500]/D");


  tree->Branch("nhBvh", &event.nhBvh,"nhBvh/I");
  tree->Branch("tidBvh", event.tidBvh,"tidBvh[500]/I");
  tree->Branch("pidBvh", event.pidBvh,"pidBvh[500]/I");
  tree->Branch("didBvh", event.didBvh,"didBvh[500]/I");
  tree->Branch("prtBvh", event.prtBvh,"prtBvh[500]/I");
  tree->Branch("qBvh", event.qBvh,"qBvh[500]/I");
  tree->Branch("xBvh", event.xBvh,"xBvh[500]/D");
  tree->Branch("yBvh", event.yBvh,"yBvh[500]/D");
  tree->Branch("zBvh", event.zBvh,"zBvh[500]/D");
  tree->Branch("pxBvh", event.pxBvh,"pxBvh[500]/D");
  tree->Branch("pyBvh", event.pyBvh,"pyBvh[500]/D");
  tree->Branch("pzBvh", event.pzBvh,"pzBvh[500]/D");
  tree->Branch("ppBvh", event.ppBvh,"ppBvh[500]/D");
  tree->Branch("deBvh", event.deBvh,"deBvh[500]/D");
  tree->Branch("tBvh", event.tBvh,"tBvh[500]/D");

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


  ////////// Bring Address From Dst
  TTreeReaderCont[kTPCGeant] = new TTreeReader( "TPC_g", TFileCont[kTPCGeant] );
  const auto& reader = TTreeReaderCont[kTPCGeant];

  TTreeCont[kTPCGeant]->SetBranchAddress("evnum", &src.evnum);
  TTreeCont[kTPCGeant]->SetBranchAddress("nhittpc", &src.nhittpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("nhPrm", &src.nhPrm);

  TTreeCont[kTPCGeant]->SetBranchAddress("NumberOfTracks",&src.NumberOfTracks);
  TTreeCont[kTPCGeant]->SetBranchAddress("PIDOfTrack",src.PIDOfTrack);
  TTreeCont[kTPCGeant]->SetBranchAddress("ParentIDOfTrack",src.ParentIDOfTrack);
  TTreeCont[kTPCGeant]->SetBranchAddress("VertexOfTrack_x",src.VertexOfTrack_x);
  TTreeCont[kTPCGeant]->SetBranchAddress("VertexOfTrack_y",src.VertexOfTrack_y);
  TTreeCont[kTPCGeant]->SetBranchAddress("VertexOfTrack_z",src.VertexOfTrack_z);
  TTreeCont[kTPCGeant]->SetBranchAddress("MomentumOfTrack",src.MomentumOfTrack);
  TTreeCont[kTPCGeant]->SetBranchAddress("MomentumOfTrack_x",src.MomentumOfTrack_x);
  TTreeCont[kTPCGeant]->SetBranchAddress("MomentumOfTrack_y",src.MomentumOfTrack_y);
  TTreeCont[kTPCGeant]->SetBranchAddress("MomentumOfTrack_z",src.MomentumOfTrack_z);
  TTreeCont[kTPCGeant]->SetBranchAddress("MomXi_x",&src.MomXi_x);
  TTreeCont[kTPCGeant]->SetBranchAddress("MomXi_y",&src.MomXi_y);
  TTreeCont[kTPCGeant]->SetBranchAddress("MomXi_z",&src.MomXi_z);
  TTreeCont[kTPCGeant]->SetBranchAddress("SpinXi_x",&src.SpinXi_x);
  TTreeCont[kTPCGeant]->SetBranchAddress("SpinXi_y",&src.SpinXi_y);
  TTreeCont[kTPCGeant]->SetBranchAddress("SpinXi_z",&src.SpinXi_z);
  TTreeCont[kTPCGeant]->SetBranchAddress("ThXi_CM",&src.ThXi_CM);
  TTreeCont[kTPCGeant]->SetBranchAddress("MomLd_x",&src.MomLd_x);
  TTreeCont[kTPCGeant]->SetBranchAddress("MomLd_y",&src.MomLd_y);
  TTreeCont[kTPCGeant]->SetBranchAddress("MomLd_z",&src.MomLd_z);
  TTreeCont[kTPCGeant]->SetBranchAddress("SpinLd_x",&src.SpinLd_x);
  TTreeCont[kTPCGeant]->SetBranchAddress("SpinLd_y",&src.SpinLd_y);
  TTreeCont[kTPCGeant]->SetBranchAddress("SpinLd_z",&src.SpinLd_z);
  TTreeCont[kTPCGeant]->SetBranchAddress("ThLd_CM",&src.ThLd_CM);


  TTreeCont[kTPCGeant]->SetBranchAddress("xPrm", src.xPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("yPrm", src.yPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("zPrm", src.zPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxPrm", src.pxPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("pyPrm", src.pyPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzPrm", src.pzPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("ppPrm", src.ppPrm);

  TTreeCont[kTPCGeant]->SetBranchAddress("ititpc", src.ititpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("idtpc", src.idtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("parentID", src.parentID);
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
  TTreeCont[kTPCGeant]->SetBranchAddress("pptpc", src.pptpc);
  //TTreeCont[kTPCGeant]->SetBranchAddress("masstpc", src.masstpc);
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

  TTreeCont[kTPCGeant]->SetBranchAddress("nhHtof", &src.nhHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("tidHtof", src.tidHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pidHtof", src.pidHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("didHtof", src.didHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("prtHtof", src.prtHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("qHtof", src.qHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("xHtof", src.xHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("yHtof", src.yHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("zHtof", src.zHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxHtof", src.pxHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pyHtof", src.pyHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzHtof", src.pzHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("ppHtof", src.ppHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("deHtof", src.deHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("tHtof", src.tHtof);



  TTreeCont[kTPCGeant]->SetBranchAddress("nhSch", &src.nhSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("tidSch", src.tidSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("pidSch", src.pidSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("didSch", src.didSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("prtSch", src.prtSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("qSch", src.qSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("xSch", src.xSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("ySch", src.ySch);
  TTreeCont[kTPCGeant]->SetBranchAddress("zSch", src.zSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxSch", src.pxSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("pySch", src.pySch);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzSch", src.pzSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("ppSch", src.ppSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("deSch", src.deSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("tSch", src.tSch);

  TTreeCont[kTPCGeant]->SetBranchAddress("nhLac", &src.nhLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("tidLac", src.tidLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("pidLac", src.pidLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("didLac", src.didLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("prtLac", src.prtLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("qLac", src.qLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("xLac", src.xLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("yLac", src.yLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("zLac", src.zLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxLac", src.pxLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("pyLac", src.pyLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzLac", src.pzLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("ppLac", src.ppLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("deLac", src.deLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("tLac", src.tLac);

  TTreeCont[kTPCGeant]->SetBranchAddress("nhFtof", &src.nhFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("tidFtof", src.tidFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pidFtof", src.pidFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("didFtof", src.didFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("prtFtof", src.prtFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("qFtof", src.qFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("xFtof", src.xFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("yFtof", src.yFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("zFtof", src.zFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxFtof", src.pxFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pyFtof", src.pyFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzFtof", src.pzFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("ppFtof", src.ppFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("deFtof", src.deFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("tFtof", src.tFtof);

  TTreeCont[kTPCGeant]->SetBranchAddress("nhSdc", &src.nhSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("tidSdc", src.tidSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pidSdc", src.pidSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("didSdc", src.didSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("prtSdc", src.prtSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("qSdc", src.qSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("xSdc", src.xSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("ySdc", src.ySdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("zSdc", src.zSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxSdc", src.pxSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pySdc", src.pySdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzSdc", src.pzSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("ppSdc", src.ppSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("deSdc", src.deSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("tSdc", src.tSdc);

  TTreeCont[kTPCGeant]->SetBranchAddress("nhBvh", &src.nhBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("tidBvh", src.tidBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("pidBvh", src.pidBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("didBvh", src.didBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("prtBvh", src.prtBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("qBvh", src.qBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("xBvh", src.xBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("yBvh", src.yBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("zBvh", src.zBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxBvh", src.pxBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("pyBvh", src.pyBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzBvh", src.pzBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("ppBvh", src.ppBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("deBvh", src.deBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("tBvh", src.tBvh);

  src.ntK18 = new TTreeReaderValue<Int_t>( *reader, "ntK18" );
  src.xvpHS = new TTreeReaderValue<vector<vector<double>>>( *reader, "xvpHS" );
  src.yvpHS = new TTreeReaderValue<vector<vector<Double_t>>>( *reader, "yvpHS" );
  src.zvpHS = new TTreeReaderValue<vector<vector<Double_t>>>( *reader, "zvpHS" );
  src.p_3rd = new TTreeReaderValue<vector<Double_t>>( *reader, "p_3rd" );
  src.xtgtHS = new TTreeReaderValue<vector<Double_t>>( *reader, "xtgtHS" );
  src.ytgtHS = new TTreeReaderValue<vector<Double_t>>( *reader, "ytgtHS" );
  src.ztgtHS = new TTreeReaderValue<vector<Double_t>>( *reader, "ztgtHS" );
  src.xoutK18 = new TTreeReaderValue<vector<Double_t>>( *reader, "xoutK18" );
  src.youtK18 = new TTreeReaderValue<vector<Double_t>>( *reader, "youtK18" );
  src.uoutK18 = new TTreeReaderValue<vector<Double_t>>( *reader, "uoutK18" );
  src.voutK18 = new TTreeReaderValue<vector<Double_t>>( *reader, "voutK18" );
  src.layerK18 = new TTreeReaderValue<vector<vector<Double_t>>>( *reader, "layerK18" );
  src.wireK18 = new TTreeReaderValue<vector<vector<Double_t>>>( *reader, "wireK18" );
  src.localhitposK18 = new TTreeReaderValue<vector<vector<Double_t>>>( *reader, "localhitposK18" );

  src.ntKurama = new TTreeReaderValue<Int_t>( *reader, "ntKurama" );
  src.xvpKurama = new TTreeReaderValue<vector<vector<Double_t>>>( *reader, "xvpKurama" );
  src.yvpKurama = new TTreeReaderValue<vector<vector<Double_t>>>( *reader, "yvpKurama" );
  src.zvpKurama = new TTreeReaderValue<vector<vector<Double_t>>>( *reader, "zvpKurama" );
  src.xtgtKurama = new TTreeReaderValue<vector<Double_t>>( *reader, "xtgtKurama" );
  src.ytgtKurama = new TTreeReaderValue<vector<Double_t>>( *reader, "ytgtKurama" );
  src.xout = new TTreeReaderValue<vector<Double_t>>( *reader, "xout" );
  src.yout = new TTreeReaderValue<vector<Double_t>>( *reader, "yout" );
  src.zout = new TTreeReaderValue<vector<Double_t>>( *reader, "zout" );
  src.pxout = new TTreeReaderValue<vector<Double_t>>( *reader, "pxout" );
  src.pyout = new TTreeReaderValue<vector<Double_t>>( *reader, "pyout" );
  src.pzout = new TTreeReaderValue<vector<Double_t>>( *reader, "pzout" );
  src.layer = new TTreeReaderValue<vector<vector<Double_t>>>( *reader, "layer" );
  src.wire = new TTreeReaderValue<vector<vector<Double_t>>>( *reader, "wire" );
  src.localhitpos = new TTreeReaderValue<vector<vector<Double_t>>>( *reader, "localhitpos" );


  return true;
}

//_____________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")   &&
      InitializeParameter<UserParamMan>("USER") &&
      InitializeParameter<TPCParamMan>("TPCPRM") &&
      InitializeParameter<HodoPHCMan>("HDPHC") );
}

//_____________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}

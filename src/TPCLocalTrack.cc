// -*- C++ -*-

//Comment by Ichikawa
//TPCLocalTrack.cc is for linear fit

//Comment by Wooseung
//Please see the discription in the TPCTrackSearch.cc
//The track coordinate origin is the target center, ***NOT TPC center***
//Track parameters are x0, y0, u0, v0 and free parameters for minuit are x0, atan2_u0, y0, atan2_v0
//x = m_x0 + m_u0*(z-tpc::ZTarget) y = m_y0 *(z-tpc::ZTarget)
//Definitions of residual
//1. ResidualVect() : Track <-> Closest point on the track (pos - calpos)
//2. ResidualVectRow() : Perpendicular to the pad direction (along the row direction)

#include "TPCLocalTrack.hh"

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <stdexcept>
#include <sstream>
#include <std_ostream.hh>

#include <TF2.h>
#include <TGraph2D.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TROOT.h>
#include <Math/Functor.h>
#include <Math/Vector3D.h>
#include <TPolyLine3D.h>
#include <Fit/Fitter.h>

#include "DetectorID.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "PrintHelper.hh"
#include "HoughTransform.hh"
#include "TPCLTrackHit.hh"
#include "TPCPadHelper.hh"
#include "UserParamMan.hh"
#include "ConfMan.hh"

#define DebugDisp 0
#define IterativeResolution 1

namespace
{
  const auto& gUser = UserParamMan::GetInstance();

  const Int_t ReservedNumOfHits = 32*10;
  //const Double_t MaxGapBtwClusters = 150.; //ref
  const Double_t MaxGapBtwClusters = 100.;
  //const Int_t MaxLayerdiffBtwClusters = 6;
  const Int_t MaxLayerdiffBtwClusters = 32;
  const Int_t MaxIteration = 50;
  const Double_t MaxChisqr = 300.;

  //for Minimization
  static Int_t gNumOfHits;
  static std::vector<TVector3> gHitPos;
  static std::vector<TVector3> gRes;
  static std::vector<Int_t> gLayer;
  static std::vector<Double_t> gPadTheta;
  static std::vector<std::vector<Double_t>> gResParam;
  static Double_t gPar[4] = {0};
  static Double_t gChisqr = 1.e+10;
  static Int_t gBadHits;
  static Int_t gMinuitStatus;
  static std::vector<Double_t> gComp;

  //x0, atan_u0, y0, atan_v0
  const Double_t FitStep[4] = { 1.0e-3, 1.0e-6, 1.0e-3, 1.0e-6};
  const Double_t LowLimit[4] = { -400., -0.5*TMath::Pi(), -400, -0.5*TMath::Pi() };
  const Double_t UpLimit[4] = { 400., 0.5*TMath::Pi(), 400, 0.5*TMath::Pi() };

  //Horizontal resolution function
  //x : alpha(track-pad angle), y : y pos of cluster (y+300 : Drift length)
  //[0] : Intrinsic XZ resolution, [1] : Attenuation term, [2] : Diffusion coefficient, [3] : Effective # of signal electrons, [4] : Pad length, [5] : Effective # of electron clusters
  static TString eq_horizontal="TMath::Sqrt(TMath::Power([0],2)+TMath::Power([2],2)*(y+300.)/([3]*TMath::Exp(-[1]*(y+300.)))+TMath::Power([4]*TMath::Tan(x),2)/(12.*[5]))";
  static TF2 *f_horizontal = new TF2("f_horizontal", eq_horizontal.Data(), -4., 4., -300., 300.);

  //Vertical resolution function
  //x : x pos of cluster (x+300 : Drift length)
  //[0] : Intrinsic Y resolution, [1] : Attenuation term, [2] : Diffusion coefficient, [3] : Effective # of signal electrons
  static TString eq_vertical="TMath::Sqrt(TMath::Power([0],2)+TMath::Power([2],2)*(x+300.)/([3]*TMath::Exp(-[1]*(x+300.))))";
  static TF1 *f_drift = new TF1("f_drift", eq_vertical.Data(), -300., 300.);

  //Window
  //Add hits into the track within the window. (ResidualWindowPull > residual/resolution)
  const Double_t ResidualWindowPullXZ = 10.;
  //const Double_t ResidualWindowPullXZ = 6.;
  const Double_t ResidualWindowPullY = 6.;
  const Double_t ResidualWindowInXZ = 5; //[mm]
  //const Double_t ResidualWindowInXZ = 10; //[mm]
  //const Double_t ResidualWindowOutXZ = 7; //[mm]
  const Double_t ResidualWindowOutXZ = 10; //[mm]
  //const Double_t ResidualWindowOutXZ = 5; //[mm]

}

//______________________________________________________________________________
static inline Bool_t CompareDist(const Int_t a, const Int_t b){

  return gComp[a] < gComp[b];
}

//______________________________________________________________________________
static inline TVector3 ResidualVect(Double_t par[4], TVector3 pos){ //Closest distance

  TVector3 x0(par[0], par[1], tpc::ZTarget);
  TVector3 x1(par[0] + par[2], par[1] + par[3], tpc::ZTarget+1.);
  TVector3 u = (x1-x0).Unit();
  TVector3 AP = pos - x0;
  Double_t dist_AX = u.Dot(AP);
  TVector3 AI(x0.x()+(u.x()*dist_AX),
	      x0.y()+(u.y()*dist_AX),
	      x0.z()+(u.z()*dist_AX));
  TVector3 d = pos-AI;
  return d;
}

//______________________________________________________________________________
static inline TVector3 ResidualVectXZ(Double_t par[4], TVector3 pos){ //Closest distance on the y=pos.y() plane

  TVector3 x0(par[0], pos.y(), tpc::ZTarget);
  TVector3 x1(par[0] + par[2], pos.y(), tpc::ZTarget+1.);
  TVector3 u = (x1-x0).Unit();
  TVector3 AP = pos - x0;
  Double_t dist_AX = u.Dot(AP);
  TVector3 AI(x0.x()+(u.x()*dist_AX),
	      x0.y()+(u.y()*dist_AX),
	      x0.z()+(u.z()*dist_AX));
  TVector3 d = pos-AI;
  return d;
}

//______________________________________________________________________________
static inline TVector3 ResidualVectRow(Double_t par[4], TVector3 pos){ //pad horizontal residual vector

  Double_t z = pos.x()*pos.x()+(pos.z()-tpc::ZTarget)*(pos.z()-tpc::ZTarget)-pos.x()*par[0];
  z /= (pos.z()-tpc::ZTarget + pos.x()*par[2]);
  TVector3 calpos(par[0] + par[2]*z,
		  par[1] + par[3]*z,
		  z + tpc::ZTarget);
  TVector3 d = pos - calpos;
  return d;
}

//_____________________________________________________________________________
static inline TVector3 CalcResolution(Double_t par[4], Int_t layer, TVector3 pos, Double_t padTheta, std::vector<Double_t> resparam, Bool_t vetoBadClusters){

  Double_t cosPad = TMath::Cos(padTheta);
  Double_t sinPad = TMath::Sin(padTheta);
  Double_t tanPad = TMath::Tan(padTheta);
  Double_t padL = tpc::padParameter[layer][tpc::kLength];
  Double_t padRadius = tpc::padParameter[layer][tpc::kRadius];
  TVector3 closestDist2TrackXZ = ResidualVectXZ(par, TVector3(0., 0., tpc::ZTarget));

  //check whether the track is crossing the layer or not
  if(vetoBadClusters && closestDist2TrackXZ.Mag() > padRadius - 0.5*padL) return TVector3(1.e+10, 1.e+10, 1.e+10);

  //Calculate resolution
  //horizontal resolution
  Double_t tanDiff = (tanPad-par[2])/(1.+tanPad*par[2]);
  Double_t alpha = TMath::ATan(tanDiff); //alpha : pad - track angle
  Double_t param_horizontal[6] = {resparam[0], resparam[1], resparam[2], resparam[3], resparam[4], resparam[5]};
  f_horizontal -> SetParameters(param_horizontal);
  Double_t res_horizontal = f_horizontal -> Eval(alpha, pos.y());

  //vertical resolution
  Double_t param_y[4] = {resparam[6], resparam[1], resparam[7], resparam[8]};
  f_drift -> SetParameters(param_y);
  Double_t res_drift = f_drift -> Eval(pos.y());

  //check whether the cluster is diffused over the layers.
  // 1. Project the cluster into the track on the XZ plane.
  // 2. For the projected point on the track, calculate a distance from the center (radius of projected position)
  // 3. By using this distance, check whether the projected point and the cluster are in the same layer or not. If not, the cluster is duffused over the layers.
  TVector3 residual_xz = ResidualVectXZ(par, pos);
  TVector3 point_projected_onTheTrack = pos - residual_xz;
  Double_t radius_projected_point = TMath::Hypot(point_projected_onTheTrack.x(), point_projected_onTheTrack.z() - tpc::ZTarget);
  //if(vetoBadClusters && radius_projected_point < padRadius - 0.5*padL) return TVector3(1.e+10, 1.e+10, 1.e+10);
  if(vetoBadClusters && radius_projected_point < padRadius - 0.5*padL) return TVector3(2.e+10, 2.e+10, 2.e+10);
  if(vetoBadClusters && radius_projected_point > padRadius + 0.5*padL) return TVector3(2.e+10, 2.e+10, 2.e+10);

  //Convert resolution along the row direction into x, y, z resolutions
  TVector3 res_row(res_horizontal*TMath::Abs(cosPad), res_drift, res_horizontal*TMath::Abs(sinPad));
  return res_row;
}

//______________________________________________________________________________
static inline Double_t CalcChi2(Double_t *par, Int_t &ndf, Bool_t vetoBadClusters)
{

  if(gHitPos.size()!=gNumOfHits || gLayer.size()!=gNumOfHits || gPadTheta.size()!=gNumOfHits || gResParam.size()!=gNumOfHits){
    hddaq::cerr << "TPCLocalTrack CalcChi2() "
		<< "#hits != # params"
		<< " # hits " << gHitPos.size()
		<< " # layer ids " << gLayer.size()
		<< " # pad thetas " << gPadTheta.size()
		<< " # res params " << gResParam.size()
		<< std::endl;
    return -1;
  }

  ndf = 0; Double_t chisqr = 0.;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 d = ResidualVectRow(par, gHitPos[i]);
    TVector3 res = CalcResolution(par, gLayer[i], gHitPos[i], gPadTheta[i], gResParam[i], vetoBadClusters);

    if(res.x() > 0.9e+10 && res.y() > 0.9e+10 && res.z() > 0.9e+10) continue; // exclude bad clusters
    chisqr += TMath::Power(d.x()/res.x(), 2) + TMath::Power(d.y()/res.y(), 2) + TMath::Power(d.z()/res.z(), 2);
    ndf++;
  }
  if(ndf < 5) return 2.*MaxChisqr;
  return chisqr/(Double_t)(ndf - 4);
}

//_____________________________________________________________________________
static void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  Double_t chisqr=0.; Int_t dof = 0;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 res = gRes[i];
    Double_t param[4] = {par[0], par[2], TMath::Tan(par[1]), TMath::Tan(par[3])};
    TVector3 d = ResidualVectRow(param, gHitPos[i]);
    if(res.x() > 0.9e+10 && res.y() > 0.9e+10 && res.z() > 0.9e+10) continue; // exclude dummy hits
    chisqr += TVector3(d.x()/res.x(), d.y()/res.y(), d.z()/res.z()).Mag2();
    dof++;
  }
  f = chisqr/(dof - 4);
}

//______________________________________________________________________________
static inline Bool_t StraightLineFit(Bool_t vetoBadClusters)
{

  if(gHitPos.size()!=gNumOfHits || gLayer.size()!=gNumOfHits || gPadTheta.size()!=gNumOfHits|| gResParam.size()!=gNumOfHits){
    hddaq::cerr << " TPCLocalTrack StraightLineFit() "
		<< "#hits != # params"
		<< " # hits " << gHitPos.size()
		<< " # layer ids " << gLayer.size()
		<< " # pad thetas " << gPadTheta.size()
		<< " # res params " << gResParam.size()
		<< std::endl;
    return false;
  }

  Double_t par[4] = {gPar[0], TMath::ATan2(gPar[2], 1.), gPar[1], TMath::ATan2(gPar[3], 1.)};
  Double_t err[4] = {-999., -999., -999., -999.};

  gRes.clear();
  for(Int_t i=0; i<gNumOfHits; i++){
    TVector3 res = CalcResolution(gPar, gLayer[i], gHitPos[i], gPadTheta[i], gResParam[i], vetoBadClusters);
    gRes.push_back(res);
  }

  //NDF value becomes different. Initialize gChisqr and fit again.
  Int_t ndf;
  if(vetoBadClusters){
    gChisqr = CalcChi2(gPar, ndf, vetoBadClusters);
    if(gChisqr>MaxChisqr) return false;
  }

  TMinuit *minuit = new TMinuit(4);
  minuit->SetPrintLevel(-1);
  minuit->SetFCN(fcn);

  Int_t ierflg = 0;
  ierflg = 1;
  Double_t arglist[10];
  arglist[0] = 1;
  minuit->mnexcm("SET NOW", arglist, 1, ierflg); // No warnings
  TString name[4] = {"x0", "atan_u0", "y0", "atan_v0"};
  for(Int_t i=0; i<4; i++){
    minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }

  minuit->Command("SET STRategy 0");
  arglist[0] = 5000.;
  arglist[1] = 0.01;
  //arglist[0] = 10000.;
  //arglist[1] = 0.001;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  //minuit->mnexcm("MINOS", arglist, 0, ierflg);
  //minuit->mnexcm("SET ERR", arglist, 2, ierflg);

  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

  Int_t parid;
  Double_t bnd1, bnd2;
  for(Int_t i=0; i<4; i++){
    minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, parid);
  }
  delete minuit;
  if(icstat==0) return false; //not calculated at all
#if DebugDisp
  if(icstat==0) std::cout<<"StraightLineFit() icstat=="<<icstat<<std::endl;
#endif

  Bool_t status = false;
  //convert x0, atan_u0, y0, atan_v0 -> x0, y0, u0, v0
  Double_t par_linear[4] = {par[0], par[2], TMath::Tan(par[1]), TMath::Tan(par[3])};
  Double_t Chisqr = CalcChi2(par_linear, ndf, vetoBadClusters);
  if(gChisqr>Chisqr){
    gBadHits = gNumOfHits - ndf;
    gChisqr = Chisqr;
    gPar[0] = par_linear[0];
    gPar[1] = par_linear[1];
    gPar[2] = par_linear[2];
    gPar[3] = par_linear[3];
    gMinuitStatus = icstat;
    status = true;
  }

  return status;
}

//_____________________________________________________________________________
TPCLocalTrack::TPCLocalTrack()
  : m_is_fitted(false),
    m_is_calculated(false),
    m_hit_array(),
    m_x0(0.), m_y0(0.), m_u0(0.), m_v0(0.),
    m_isAccidental(0),
    m_closedist(1.e+10,1.e+10, 1.e+10),
    m_edgepoint(0.,0.,-143.),
    m_chisqr(1.e+10),
    m_minuit(0),
    m_n_iteration(0),
    m_fitflag(0), m_vtxflag(0),
    m_searchtime(0), m_fittime(0),
    m_trackid(-1),
    m_x0_exclusive(), m_y0_exclusive(),
    m_u0_exclusive(), m_v0_exclusive(),
    m_chisqr_exclusive()
{
  m_hit_array.reserve(ReservedNumOfHits);
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCLocalTrack::~TPCLocalTrack()
{
  debug::ObjectCounter::decrease(ClassName());
}

//______________________________________________________________________________
TPCLocalTrack::TPCLocalTrack(TPCLocalTrack *init){

  this -> m_is_fitted = false;
  this -> m_is_calculated = false;
  for(Int_t i=0;i<init -> m_hit_array.size();i++){
    this -> m_hit_array.push_back(new TPCLTrackHit(init -> m_hit_array[i] -> GetHit()));
  }

  this -> m_x0 = init -> m_x0 ;
  this -> m_y0 = init -> m_y0 ;
  this -> m_u0 = init -> m_u0 ;
  this -> m_v0 = init -> m_v0 ;
  this -> m_isAccidental = init -> m_isAccidental ;
  this -> m_closedist = init -> m_closedist ;
  this -> m_edgepoint = init -> m_edgepoint ;
  this -> m_chisqr = init -> m_chisqr ;
  this -> m_minuit = init -> m_minuit ;
  this -> m_n_iteration = init -> m_n_iteration ;
  this -> m_fitflag = init -> m_fitflag ;
  this -> m_vtxflag = init -> m_vtxflag ;
  this -> m_searchtime = init -> m_searchtime ; //millisec
  this -> m_fittime = init -> m_fittime ; //millisec
  this -> m_trackid = init -> m_trackid ; //for k18, kurama tracks
  for(Int_t i=0;i<init -> m_x0_exclusive.size();i++){
    this -> m_x0_exclusive.push_back(init -> m_x0_exclusive[i]);
    this -> m_y0_exclusive.push_back(init -> m_y0_exclusive[i]);
    this -> m_u0_exclusive.push_back(init -> m_u0_exclusive[i]);
    this -> m_v0_exclusive.push_back(init -> m_v0_exclusive[i]);
    this -> m_chisqr_exclusive.push_back(init -> m_chisqr_exclusive[i]);
  }
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
void
TPCLocalTrack::AddTPCHit(TPCLTrackHit *hit)
{

  if(hit) m_hit_array.push_back(hit);
}

//_____________________________________________________________________________
void
TPCLocalTrack::Calculate()
{
  if(IsCalculated()){
    hddaq::cerr << FUNC_NAME << " already called" << std::endl;
    return;
  }

  const std::size_t n = m_hit_array.size();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hit = m_hit_array[i];
    hit->SetCalX0Y0(m_x0, m_y0);
    hit->SetCalUV(m_u0, m_v0);
    hit->SetCalPosition(hit->GetLocalCalPos());
    hit->SetResolution(GetResolutionVect(i, true));
  }
  m_is_calculated = true;
}

//_____________________________________________________________________________
void
TPCLocalTrack::CalculateExclusive()
{

  if(!IsCalculated()){
    hddaq::cerr << "#W " << FUNC_NAME << " No inclusive calculation" << std::endl;
    return;
  }

  const std::size_t n = m_hit_array.size();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hit = m_hit_array[i];
    hit->SetCalX0Y0Exclusive(m_x0_exclusive[i], m_y0_exclusive[i]);
    hit->SetCalUVExclusive(m_u0_exclusive[i], m_v0_exclusive[i]);
    hit->SetCalPositionExclusive(hit->GetLocalCalPosExclusive());
  }
}

//_____________________________________________________________________________
void
TPCLocalTrack::ClearHits()
{
  m_hit_array.clear();
}

//_____________________________________________________________________________
void
TPCLocalTrack::DeleteNullHit()
{

  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hit = m_hit_array[i];
    if(!hit->IsGood()){
      hddaq::cout << FUNC_NAME << " "
		  << "null hit has been deleted" << std::endl;
      m_hit_array.erase(m_hit_array.begin()+i);
      --i;
      gHitPos.clear();
      gLayer.clear();
      gPadTheta.clear();
      gResParam.clear();
      gRes.clear();
    }
  }
}

//______________________________________________________________________________
Bool_t
TPCLocalTrack::IsGoodForTracking()
{

  const std::size_t n = m_hit_array.size();
  if(GetNDF()<1){
#if DebugDisp
    hddaq::cerr << "#W " << FUNC_NAME << " "
		<< "Min layer should be > NDF" << std::endl;
#endif
    return false;
  }

  if(n>ReservedNumOfHits){
#if DebugDisp
    hddaq::cerr << "#W " << FUNC_NAME << " "
		<< "n > ReservedNumOfHits" << std::endl;
#endif
    return false;
  }

  Int_t npad = GetNPad();
  if(npad<=1) return false; //frequently many noise hits appear in a single pad

  return true;
}

//______________________________________________________________________________
Int_t
TPCLocalTrack::GetNPad() const
{

  const std::size_t n = m_hit_array.size();
  std::vector<Int_t> pads;
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hit = m_hit_array[i];
    pads.push_back(hit -> GetHit() -> GetPad());
  }
  std::sort(pads.begin(),pads.end());
  pads.erase(std::unique(pads.begin(),pads.end()),pads.end());
  return pads.size();
}

//______________________________________________________________________________
Int_t
TPCLocalTrack::GetNDF() const
{
  Int_t nhit = GetNHit();
  return nhit - 4;
}

//______________________________________________________________________________
TVector3
TPCLocalTrack::GetPosition(Double_t z) const
{

  TVector3 pos(m_x0 + m_u0*(z-tpc::ZTarget), m_y0 + m_v0*(z-tpc::ZTarget), z);
  return pos;
}

//_____________________________________________________________________________
void
TPCLocalTrack::CalcClosestDistTgt()
{

  TVector3 x0(m_x0, m_y0, tpc::ZTarget);
  TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, tpc::ZTarget+1.);
  TVector3 tgt(0., 0., tpc::ZTarget);
  TVector3 unit_v = (x1-x0).Unit();
  TVector3 diff_v = tgt - x0;
  TVector3 dist = diff_v.Cross(unit_v);
  m_closedist = dist;

}

//______________________________________________________________________________
Bool_t
TPCLocalTrack::VertexAtTarget()
{

  Bool_t status = false;
  if(m_closedist.Mag() < tpc::TargetVtxWindow) status = true;
  return status;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrack::IsBackward()
{

  if(!VertexAtTarget()) return false;
  if(m_edgepoint.z() > tpc::ZTarget) return false;
  TVector3 pos = GetPosition(-250.); //At z=-250.
  //if(TMath::Abs(pos.X()) > 50.) return false;
  if(TMath::Abs(pos.X()) > 75.) return false;
  return true;
}

//______________________________________________________________________________
void
TPCLocalTrack::SetClustersHoughFlag(Int_t hough_flag)
{
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TPCHit *hit = hitp->GetHit();
    if( !hit ) continue;
    hit->SetHoughFlag(hough_flag);
  }
}

//_____________________________________________________________________________
Bool_t
TPCLocalTrack::DoFit(Int_t MinHits)
{
  if(!IsGoodForTracking()) return false;

  Bool_t status = DoStraightTrackFit();
  m_is_fitted = status;
  if(!status || m_chisqr > MaxChisqr) return false;

  //Vtx in the target, need to check whether two tracks are merged or not
  if(!SeparateTracksAtTarget()) return false;
  //If track splitting is performed, m_is_fitted flag is set to false
  if(!m_is_fitted) return DoFit(MinHits); //Do chisqr minimization again after separation

  //Minimum # of clusters
#if 1
  Int_t nhit = GetNHit();
  if(nhit<MinHits) return false;
#else
  Int_t nbadhit = gBadHits;
  if(nhit-nbadhit<MinHits) return false;
#endif
  return status;
}

//_____________________________________________________________________________
void
TPCLocalTrack::EraseHits(std::vector<Int_t> delete_hits)
{

  std::sort(delete_hits.begin(), delete_hits.end());
  for(Int_t i=0; i<delete_hits.size(); ++i){
    //Reset houghflag
    TPCLTrackHit *hitp = m_hit_array[delete_hits[i]-i];
    TPCHit *hit = hitp->GetHit();
    hit->SetHoughFlag(0);

    m_hit_array.erase(m_hit_array.begin()+delete_hits[i]-i);
  }
  gHitPos.clear();
  gLayer.clear();
  gPadTheta.clear();
  gResParam.clear();
  gRes.clear();

}

//_____________________________________________________________________________
void
TPCLocalTrack::EraseHit(Int_t delete_hit)
{

  //Reset houghflag
  TPCLTrackHit *hitp = m_hit_array[delete_hit];
  TPCHit *hit = hitp->GetHit();
  hit->SetHoughFlag(0);

  m_hit_array.erase(m_hit_array.begin()+delete_hit);
  gHitPos.clear();
  gLayer.clear();
  gPadTheta.clear();
  gResParam.clear();
  gRes.clear();
}

//_____________________________________________________________________________
Bool_t
TPCLocalTrack::DoStraightTrackFit()
{
  Bool_t vetoBadClusters = false;

  DeleteNullHit();
  const Int_t n = m_hit_array.size();
  if(!IsGoodForTracking()) return false;
  gNumOfHits = n;
  gHitPos.clear();
  gLayer.clear();
  gPadTheta.clear();
  gResParam.clear();
  gBadHits = 0;
  gMinuitStatus = 0;

  //Initialization params with Hough-transform result or previous tracking result
  gPar[0] = m_x0;
  gPar[1] = m_y0;
  gPar[2] = m_u0;
  gPar[3] = m_v0;

  TVector3 tgt = TVector3(0., 0., tpc::ZTarget);
  m_edgepoint = tgt;
  for(Int_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    gHitPos.push_back(pos);
    Int_t layer = hitp->GetLayer();
    gLayer.push_back(layer);
    Double_t padTheta = hitp->GetPadTheta();
    gPadTheta.push_back(padTheta);
    std::vector<Double_t> resparam = hitp->GetResolutionParams();
    gResParam.push_back(resparam);
    TVector3 hit2tgt = pos - tgt;
    TVector3 dist2tgt = m_edgepoint - tgt;
    if(dist2tgt.Mag() < hit2tgt.Mag()) m_edgepoint = pos;
  }

  Int_t ndf;
  gChisqr = CalcChi2(gPar, ndf, vetoBadClusters);

#if DebugDisp
  std::cout<<FUNC_NAME+" Before linear fitting"
	   <<" n hits: "<<gNumOfHits
	   <<" bad hits: "<<gBadHits
	   <<" n_iteration: "<<m_n_iteration
	   <<" chisqr: "<<gChisqr<<std::endl
	   <<", x0: "<<m_x0
	   <<", y0: "<<m_y0
	   <<", u0: "<<m_u0
	   <<", v0: "<<m_v0
	   <<std::endl;
#endif

  if(!StraightLineFit(vetoBadClusters)) return false;
  m_chisqr = gChisqr;
  m_x0 = gPar[0];
  m_y0 = gPar[1];
  m_u0 = gPar[2];
  m_v0 = gPar[3];
  m_minuit = gMinuitStatus;

#if DebugDisp
  std::cout<<FUNC_NAME+" After linear fitting"
	   <<" n hits: "<<gNumOfHits
	   <<" bad hits: "<<gBadHits
	   <<" n_iteration: "<<m_n_iteration
	   <<" chisqr: "<<m_chisqr<<std::endl
	   <<", x0: "<<m_x0
	   <<", y0: "<<m_y0
	   <<", u0: "<<m_u0
	   <<", v0: "<<m_v0
	   <<std::endl;
#endif

  Int_t delete_hit = -1;
  Int_t false_layer = 0;
  Double_t Max_residual = -100.;
  for(Int_t i=0; i<m_hit_array.size(); ++i){
    Double_t resi = 0;
    if(!ResidualCheck(i, resi)){
      if(Max_residual<resi){
	Max_residual = resi;
	delete_hit = i;
      }
      ++false_layer;
    }
  }

  if(false_layer>0){
#if DebugDisp
    TPCLTrackHit *hitp = m_hit_array[delete_hit];
    TVector3 pos = hitp->GetLocalHitPos();
    std::cout<<"delete hits ["<<delete_hit<<"]=("
    	     <<pos.x()<<", "
      	     <<pos.y()<<", "
      	     <<pos.z()<<")"<<std::endl;
#endif
    EraseHit(delete_hit);
  }

  m_n_iteration++;
  if(m_n_iteration > MaxIteration) return false;
  if(false_layer == 0){ //Tracking is over.
#if IterativeResolution

    //Now excluding bad clusters and fitting again.
    vetoBadClusters = true;
    Bool_t update = StraightLineFit(vetoBadClusters);
    if(update){
      m_chisqr = gChisqr;
      m_x0 = gPar[0];
      m_y0 = gPar[1];
      m_u0 = gPar[2];
      m_v0 = gPar[3];
      m_minuit = gMinuitStatus;
      //iteration process
      update = StraightLineFit(vetoBadClusters);
      if(update){
	m_chisqr = gChisqr;
	m_x0 = gPar[0];
	m_y0 = gPar[1];
	m_u0 = gPar[2];
	m_v0 = gPar[3];
	m_minuit = gMinuitStatus;
      }
    }
    if(TMath::Abs(m_chisqr-MaxChisqr)<0.1) return false;

#if DebugDisp
    std::cout<<FUNC_NAME+" After excluding bad hits"
	     <<" n hits: "<<gNumOfHits
	     <<" bad hits: "<<gBadHits
	     <<" n_iteration: "<<m_n_iteration
	     <<" chisqr: "<<m_chisqr<<std::endl
	     <<", x0: "<<m_x0
	     <<", y0: "<<m_y0
	     <<", u0: "<<m_u0
	     <<", v0: "<<m_v0
	     <<std::endl;
#endif
#endif
    return true;
  }
  else return DoStraightTrackFit();
}

//_____________________________________________________________________________
TVector3
TPCLocalTrack::GetResolutionVect(Int_t i, Bool_t vetoBadClusters){

  auto hit = m_hit_array[i];
  if(vetoBadClusters && hit->IsOnTheFrame()) return TVector3(3.e+10, 3.e+10, 3.e+10);
  TVector3 pos = hit->GetLocalHitPos();
  Int_t layer = hit->GetLayer();
  Double_t padTheta = hit->GetPadTheta();
  std::vector<Double_t> resParam = hit->GetResolutionParams();
  Double_t par[4] = {m_x0, m_y0, m_u0, m_v0};

  //Convert resolution along the row direction into closets point's x, y, z resolutions
  TVector3 residual_row = ResidualVectRow(par, pos);
  TVector3 residual = ResidualVect(par, pos);
  TVector3 res_row = CalcResolution(par, layer, pos, padTheta, resParam, vetoBadClusters);
  double res_x = TMath::Abs(res_row.x()*residual.x()/residual_row.x());
  double res_y = TMath::Abs(res_row.y()*residual.y()/residual_row.y());
  double res_z = TMath::Abs(res_row.z()*residual.z()/residual_row.z());
  return TVector3(res_x, res_y, res_z);
}

//_____________________________________________________________________________
TVector3
TPCLocalTrack::CalcResidual(TVector3 position)
{

  Double_t par[4] = {m_x0, m_y0, m_u0, m_v0};
  return ResidualVect(par, position);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrack::GetHorizontalResolution(TPCHit *hit)
{

  TVector3 pos = hit->GetPosition();
  std::vector<Double_t> resparam = hit->GetResolutionParams();
  Double_t param_horizontal[6] = {resparam[0], resparam[1], resparam[2], resparam[3], resparam[4], resparam[5]};
  f_horizontal -> SetParameters(param_horizontal);
  Double_t res_horizontal = f_horizontal -> Eval(GetAlpha(hit), pos.y());

  return res_horizontal;
}

//_____________________________________________________________________________
Double_t
TPCLocalTrack::GetHorizontalResolution(Int_t i)
{
  TPCHit *hit = m_hit_array[i] -> GetHit();
  return GetHorizontalResolution(hit);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrack::GetHorizontalResidual(TPCHit *hit)
{

  //alpha : pad - track angle
  TVector3 pos = hit->GetPosition();
  Double_t par[4] = {m_x0, m_y0, m_u0, m_v0};
  Double_t padTheta = hit->GetPadTheta();
  Double_t cosTheta = TMath::Cos(padTheta);
  Double_t sinTheta = TMath::Sin(padTheta);
  TVector3 rotvecx_global2pad(cosTheta, 0., -sinTheta);
  //Local x coordinate
  TVector3 residual_row = ResidualVectRow(par, pos);
  double residual_horizontal = rotvecx_global2pad*residual_row;

  return residual_horizontal;
}

//_____________________________________________________________________________
Double_t
TPCLocalTrack::GetHorizontalResidual(Int_t i)
{
  TPCHit *hit = m_hit_array[i] -> GetHit();
  return GetHorizontalResidual(hit);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrack::GetHorizontalResidualExclusive(Int_t i)
{

  //alpha : pad - track angle
  TPCLTrackHit *hit = m_hit_array[i];
  TVector3 pos = hit->GetLocalHitPos();
  Double_t par[4] = {m_x0_exclusive[i], m_y0_exclusive[i], m_u0_exclusive[i], m_v0_exclusive[i]};
  Double_t padTheta = hit->GetPadTheta();
  Double_t cosTheta = TMath::Cos(padTheta);
  Double_t sinTheta = TMath::Sin(padTheta);
  TVector3 rotvecx_global2pad(cosTheta, 0., -sinTheta);

  //Local x coordinate
  TVector3 residual_row = ResidualVectRow(par, pos);
  double residual_horizontal = rotvecx_global2pad*residual_row;

  return residual_horizontal;
}

//_____________________________________________________________________________
Double_t
TPCLocalTrack::GetVerticalResolution(TPCHit *hit)
{

  TVector3 pos = hit->GetPosition();
  std::vector<Double_t> resparam = hit->GetResolutionParams();
  Double_t param_y[4] = {resparam[6], resparam[1], resparam[7], resparam[8]};
  f_drift -> SetParameters(param_y);
  Double_t res_drift = f_drift -> Eval(pos.y());

  return res_drift;
}

//_____________________________________________________________________________
Double_t
TPCLocalTrack::GetVerticalResolution(Int_t i)
{
  TPCHit *hit = m_hit_array[i] -> GetHit();
  return GetVerticalResolution(hit);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrack::GetVerticalResidual(TPCHit *hit)
{

  TVector3 pos = hit->GetPosition();
  Double_t par[4] = {m_x0, m_y0, m_u0, m_v0};
  TVector3 residual_row = ResidualVectRow(par, pos);

  return residual_row.y();
}

//_____________________________________________________________________________
Double_t
TPCLocalTrack::GetVerticalResidual(Int_t i)
{
  TPCHit *hit = m_hit_array[i] -> GetHit();
  return GetVerticalResidual(hit);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrack::GetVerticalResidualExclusive(Int_t i)
{

  TPCLTrackHit *hit = m_hit_array[i];
  TVector3 pos = hit->GetLocalHitPos();
  Double_t par[4] = {m_x0_exclusive[i], m_y0_exclusive[i], m_u0_exclusive[i], m_v0_exclusive[i]};
  TVector3 residual_row = ResidualVectRow(par, pos);

  return residual_row.y();
}

//_____________________________________________________________________________
Bool_t
TPCLocalTrack::ResidualCheck(Int_t i, Double_t &residual)
{

  TPCHit *hit = m_hit_array[i] -> GetHit();
  //return ResidualCheck(hit, residual);
  return IsGoodHitToAdd(hit, residual);
}

//_____________________________________________________________________________
Bool_t
TPCLocalTrack::IsGoodHitToAdd(TPCHit *hit, Double_t &residual)
{

  Int_t layer = hit->GetLayer();
  Double_t par[4] = {m_x0, m_y0, m_u0, m_v0};
  TVector3 position = hit->GetPosition();
  TVector3 resi = ResidualVect(par, position); //Closest distance
  residual = resi.Mag();

  Double_t closest_gap = 9999;
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    if( !hitp ) continue;
    TVector3 pos = hitp -> GetLocalHitPos();
    TVector3 gap = pos - position;
    if(closest_gap > gap.Mag()) closest_gap = gap.Mag();
  }
  if(closest_gap > MaxGapBtwClusters) return false;

  //Residual/resolution < window
  Double_t residual_vertical = GetVerticalResidual(hit);
  Double_t residual_horizontal = GetHorizontalResidual(hit);
  Double_t resolution_vertical = GetVerticalResolution(hit);
  Double_t resolution_horizontal = GetHorizontalResolution(hit);
  if(TMath::Abs(residual_vertical) > resolution_vertical*ResidualWindowPullY) return false;
  if(TMath::Abs(residual_horizontal) > resolution_horizontal*ResidualWindowPullXZ) return false;

  //XZ residual < window
  TVector3 residualXZ = ResidualVectXZ(par, position);
  if(layer < 10 && residualXZ.Mag() > ResidualWindowInXZ) return false;
  else if(layer >=10 && residualXZ.Mag() > ResidualWindowOutXZ) return false;

  return true;
}

//_____________________________________________________________________________
Int_t
TPCLocalTrack::Side(TVector3 pos)
{

  Int_t flag = -1;
  TVector3 tgt(0., 0., tpc::ZTarget);
  TVector3 pos_ = pos - tgt;
  TVector3 division(1./m_u0, 0., -1.);
  TVector3 norm = pos_.Cross(division);
  if(norm.Y()>0.) flag = 1;

  return flag;
}

//_____________________________________________________________________________
Double_t
TPCLocalTrack::GetTheta() const
{
  Double_t cost = 1./TMath::Sqrt(1.+m_u0*m_u0+m_v0*m_v0);
  return TMath::ACos(cost)*TMath::RadToDeg();
}

//_____________________________________________________________________________
Double_t
TPCLocalTrack::GetAlpha(TPCHit* hit) const
{
  //alpha : pad - track angle
  Double_t padTheta = hit->GetPadTheta();
  Double_t tanPad = TMath::Tan(padTheta);
  Double_t tanDiff = (tanPad-m_u0)/(1.+tanPad*m_u0);
  Double_t alpha = TMath::ATan(tanDiff);
  return alpha;
}

//_____________________________________________________________________________
Double_t
TPCLocalTrack::GetAlpha(Int_t i) const
{

  //alpha : pad - track angle
  TPCHit *hit = m_hit_array[i] -> GetHit();
  return GetAlpha(hit);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrack::GetHitLength(Int_t i) const
{

  //alpha : pad - track angle
  TPCLTrackHit *hit = m_hit_array[i];
  Int_t layer = hit->GetLayer();
  Double_t padTheta = hit->GetPadTheta();
  Double_t tanPad = TMath::Tan(padTheta);
  Double_t tanDiff = (tanPad-m_u0)/(1.+tanPad*m_u0);
  Double_t path = tpc::padParameter[layer][5];
  path *= TMath::Sqrt(1.+tanDiff*tanDiff)*TMath::Sqrt(1.+m_v0*m_v0);

  return path;
}

//_____________________________________________________________________________
void
TPCLocalTrack::Print(const TString& arg, Bool_t print_allhits) const
{

  TString tracksize = Form(" #clusters = %d", (Int_t)m_hit_array.size());
  TString trackpar = Form(", params x0:%f y0:%f u0:%f v0:%f", m_x0, m_y0, m_u0, m_v0);
  std::cout<<arg.Data()<<std::endl;
  std::cout<<"Track info : "<<tracksize.Data()<<trackpar.Data()<<std::endl;
  if(print_allhits){
    for(std::size_t i=0; i<m_hit_array.size(); ++i){
      TPCLTrackHit *hitp = m_hit_array[i];
      if( !hitp ) continue;
      hitp->Print();
    }
  }
}

//______________________________________________________________________________
Double_t
TPCLocalTrack::GetTrackdE()
{

  Double_t dE = 0;
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    dE += hitp->GetDe();
  }
  return dE;
}

//_____________________________________________________________________________
void
TPCLocalTrack::DoFitExclusive()
{
  const Int_t n = m_hit_array.size();
  m_x0_exclusive.resize(n);
  m_y0_exclusive.resize(n);
  m_u0_exclusive.resize(n);
  m_v0_exclusive.resize(n);
  m_chisqr_exclusive.resize(n);

  gNumOfHits = n-1;
  for(Int_t i=0; i<n; ++i){
    gHitPos.clear();
    gLayer.clear();
    gPadTheta.clear();
    gResParam.clear();
    gHitPos.resize(n-1);
    gLayer.resize(n-1);
    gPadTheta.resize(n-1);
    gResParam.resize(n-1);

    gChisqr = 1.e+10;
    gPar[0] = m_x0;
    gPar[1] = m_y0;
    gPar[2] = m_u0;
    gPar[3] = m_v0;

    Int_t flag=0;
    for(Int_t j=0; j<n; ++j){
      if(j==i) continue; //exclude ith hit

      TPCLTrackHit *hitp = m_hit_array[j];
      TVector3 pos = hitp->GetLocalHitPos();
      Int_t layer = hitp->GetLayer();
      Double_t padTheta = hitp->GetPadTheta();
      std::vector<Double_t> resparam = hitp->GetResolutionParams();
      gHitPos[flag] = pos;
      gLayer[flag] = layer;
      gPadTheta[flag] = padTheta;
      gResParam[flag] = resparam;
      flag++;
    } //j

    //Track fitting
    Bool_t vetoBadClusters = true;
    StraightLineFit(vetoBadClusters);

    m_chisqr_exclusive[i] = gChisqr;
    m_x0_exclusive[i] = gPar[0];
    m_y0_exclusive[i] = gPar[1];
    m_u0_exclusive[i] = gPar[2];
    m_v0_exclusive[i] = gPar[3];
  } //i
}

//______________________________________________________________________________
Bool_t
TPCLocalTrack::SeparateTracksAtTarget()
{

  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);
  m_vtxflag = 0;

  Bool_t status = true;
  if(BeamThroughTPC){
    return status;
  }

  std::vector<Int_t> side1_hits; std::vector<Int_t> side2_hits;
  CalcClosestDistTgt(); //Distance between the target & the track
  if(VertexAtTarget()){

    const std::size_t n = m_hit_array.size();
    for(std::size_t i=0; i<n; ++i){
      TPCLTrackHit *hitp = m_hit_array[i];
      TVector3 pos = hitp -> GetLocalHitPos();
      if(Side(pos)==1) side1_hits.push_back(i);
      else if(Side(pos)==-1) side2_hits.push_back(i);
    }
  }
  else{ //If track is not crossing the target.
    gComp.clear();
    std::vector<Int_t> m_hit_order;
    const std::size_t n = m_hit_array.size();
    for(std::size_t i=0; i<n; ++i){
      TPCLTrackHit *hitp = m_hit_array[i];
      TVector3 pos = hitp -> GetLocalHitPos();
      TVector3 tgt(0., 0., tpc::ZTarget);
      TVector3 pos_ = pos - tgt;
      TVector3 division(1./m_u0, 0., -1.);
      TVector3 norm = pos_.Cross(division);
      gComp.push_back(norm.Y());
      m_hit_order.push_back(i);
    }

    for(std::size_t i=0; i<m_hit_order.size(); ++i){
      std::sort(m_hit_order.begin(), m_hit_order.end(), CompareDist);
    }

    //Exclude the beam hit from the scattered track.
    Bool_t flag = false;
    TVector3 prev_pos;
    Bool_t prev_isBeamHit = false; Bool_t isBeamHit;
    for(Int_t i=0; i<n; ++i){
      Int_t id = m_hit_order[i];
      TPCLTrackHit *hitp = m_hit_array[id];
      TVector3 pos = hitp -> GetLocalHitPos();
      if(TMath::Abs(pos.x()) < 25. &&
	 TMath::Abs(pos.y()) < 30. &&
	 pos.z() < tpc::ZTarget) isBeamHit = true;
      else isBeamHit = false;

      TVector3 gap = pos - prev_pos;
      //std::cout<<i<<" gap "<<gap.Mag()<<" pos "<<pos<<std::endl;
      if(m_v0 > 0.01 && i!=0 && gap.Mag() > 30. && prev_isBeamHit!=isBeamHit) flag = true;
      if(!flag) side1_hits.push_back(id);
      else side2_hits.push_back(id);
      prev_pos = pos;
      prev_isBeamHit = isBeamHit;
    }
  }

  if(side1_hits.size()==0 || side2_hits.size()==0) return status;
  else{ //exclude the shorter side
    m_is_fitted = false; //Need to do minimization again
    m_n_iteration = 0;

    if(side1_hits.size() >= side2_hits.size()) EraseHits(side2_hits);
    else EraseHits(side1_hits);

    //After separation, recalculate track params.
    status = IsGoodForTracking();
    if(status){

      const Int_t n_remain = m_hit_array.size();
      gNumOfHits = n_remain;
      gHitPos.clear();
      gLayer.clear();
      gPadTheta.clear();
      gResParam.clear();
      for(Int_t i=0; i<n_remain; ++i){
	TPCLTrackHit *hitp = m_hit_array[i];
	TVector3 pos = hitp->GetLocalHitPos();
	gHitPos.push_back(pos);
	Int_t layer = hitp->GetLayer();
	gLayer.push_back(layer);
	Double_t padTheta = hitp->GetPadTheta();
	gPadTheta.push_back(padTheta);
	std::vector<Double_t> resparam = hitp->GetResolutionParams();
	gResParam.push_back(resparam);
      }

      gPar[0] = m_x0;
      gPar[1] = m_y0;
      gPar[2] = m_u0;
      gPar[3] = m_v0;
      Int_t ndf;
      gChisqr = CalcChi2(gPar, ndf, false);

      //Track fitting
      StraightLineFit(false);
      m_x0 = gPar[0];
      m_y0 = gPar[1];
      m_u0 = gPar[2];
      m_v0 = gPar[3];
      m_chisqr = gChisqr;

      CalcClosestDistTgt();
      TVector3 pos = m_hit_array[0] -> GetLocalHitPos();
      if(VertexAtTarget()) m_vtxflag = Side(pos);
      else m_vtxflag = 0;
    }
  }

#if DebugDisp
  if(!status) std::cout<<FUNC_NAME+" Separated track is not good for tracking"<<std::endl;
  std::cout<<FUNC_NAME+" Close distance from the target to the track : "<<m_closedist.Mag()<<std::endl;
#endif

  return status;
}

//______________________________________________________________________________
//If two tracks are merged at the target, separate them and recalculate params.
Bool_t
TPCLocalTrack::SeparateClustersWithGap()
{

  Bool_t status = true;

  gComp.clear();
  std::vector<Int_t> m_hit_order;
  const std::size_t n = m_hit_array.size();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp -> GetLocalHitPos();
    TVector3 tgt(0., 0., tpc::ZTarget);
    TVector3 pos_ = pos - tgt;
    TVector3 division(1./m_u0, 0., -1.);
    TVector3 norm = pos_.Cross(division);
    gComp.push_back(norm.Y());
    m_hit_order.push_back(i);
  }

  for(std::size_t i=0; i<m_hit_order.size(); ++i){
    std::sort(m_hit_order.begin(), m_hit_order.end(), CompareDist);
  }

  std::vector<Int_t> side1_hits; std::vector<Int_t> side2_hits;
  Int_t prev_layer; TVector3 prev_pos(0, 0, 0); Bool_t flip = false;
  for(std::size_t i=0; i<n; ++i){
    Int_t id = m_hit_order[i];
    TPCLTrackHit *hitp = m_hit_array[id];
    Int_t layer = hitp -> GetLayer();
    TVector3 pos = hitp -> GetLocalHitPos();
    TVector3 gap = pos - prev_pos;
    if(i!=0 && (TMath::Abs(layer - prev_layer) > MaxLayerdiffBtwClusters || gap.Mag() > MaxGapBtwClusters)) flip = true;
    if(!flip) side1_hits.push_back(id);
    else side2_hits.push_back(id);
    prev_layer = layer;
    prev_pos = pos;
  }

  //Check a gap between clusters.
  if(side1_hits.size()==0 || side2_hits.size()==0) return status;
  else{ //exclude the shorter side
    m_is_fitted = false; //Need to do minimization again
    m_n_iteration = 0;
    if(side1_hits.size() >= side2_hits.size()) EraseHits(side2_hits);
    else EraseHits(side1_hits);

    //After separation, recalculate track params.
    status = IsGoodForTracking();
    if(status){

      const Int_t n_remain = m_hit_array.size();
      gNumOfHits = n_remain;
      gHitPos.clear();
      gLayer.clear();
      gPadTheta.clear();
      gResParam.clear();
      for(Int_t i=0; i<n_remain; ++i){
	TPCLTrackHit *hitp = m_hit_array[i];
	TVector3 pos = hitp->GetLocalHitPos();
	gHitPos.push_back(pos);
	Double_t padTheta = hitp->GetPadTheta();
	gPadTheta.push_back(padTheta);
	Int_t layer = hitp->GetLayer();
	gLayer.push_back(layer);
	std::vector<Double_t> resparam = hitp->GetResolutionParams();
	gResParam.push_back(resparam);
      }

      gPar[0] = m_x0;
      gPar[1] = m_y0;
      gPar[2] = m_u0;
      gPar[3] = m_v0;
      Int_t ndf;
      gChisqr = CalcChi2(gPar, ndf, false);

      //Track fitting
      StraightLineFit(false);
      m_x0 = gPar[0];
      m_y0 = gPar[1];
      m_u0 = gPar[2];
      m_v0 = gPar[3];
      m_chisqr = gChisqr;
    }
  }

#if DebugDisp
  if(!status) std::cout<<FUNC_NAME+" Separated clusters are not good for tracking"<<std::endl;
#endif

  return status;
}

//_____________________________________________________________________________
void
TPCLocalTrack::RecalcTrack()
{
  m_is_calculated = true;
  m_is_fitted = true;
  m_is_calculated = true;
  CalcClosestDistTgt();
  const std::size_t n = m_hit_array.size();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hit = m_hit_array[i];
    hit->SetCalX0Y0(m_x0, m_y0);
    hit->SetCalUV(m_u0, m_v0);
    hit->SetCalPosition(hit->GetLocalCalPos());
    hit->SetResolution(GetResolutionVect(i, true));
  }

}

//_____________________________________________________________________________
void
TPCLocalTrack::CheckIsAccidental()
{

  if(!m_is_fitted) return;

  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);
  if(BeamThroughTPC || m_isAccidental==1) return;
  if(m_u0>0.02 || m_v0>0.02 || m_x0 > 40.) return;

  Int_t nhit_upstream_tgt = 0; Int_t nhit_downstream_tgt = 0;
  const std::size_t n = m_hit_array.size();
  for(Int_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp -> GetLocalHitPos();
    if(TMath::Abs(pos.x()) < 40. && pos.z() < tpc::ZTarget) nhit_upstream_tgt++;
    if(TMath::Abs(pos.x()) < 40. && pos.z() > tpc::ZTarget) nhit_downstream_tgt++;
  }
  if(nhit_upstream_tgt>1 && nhit_downstream_tgt>=5) m_isAccidental=1;

}

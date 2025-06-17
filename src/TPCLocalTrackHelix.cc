// -*- C++ -*-

/*
//Comment by Ichikawa
//TPCLocalTrackHelix.cc is for Helix fit
//Pre circle fit (TMinuit -> reduced chi2 method)

//Comment by Wooseung
Please see the discription in the TPCTrackSearch.cc
The track coordinate origin is the target center, ***NOT TPC center***

//Equation of helix
x = -X, y = Z - tpc::ZTarget, z = Y;
x = p[0] + p[3]*cos(theta);
y = p[1] + p[3]*sin(theta);
z = p[2] + p[4]*p[3]*(theta);

//Containers
1. m_hit_array : cluster(hit) container
2. m_hit_t : theta value of the cluster
3. m_hit_order : order of clusters along the track (ascending order in theta)

//FCN functions for chisqr2/ndf minimization by Minuit
1. fcn_circle
2. fcn_line
3. fcn_helix

//for K1.8, Kurama tracking, fix the momentum parameter and do mimization process.

//Fitting process
1. Preliminary fitting (Get helix initial parameters)
1-1. Circle fiting on the horizontal plane

+Optional step Line searching Hough-Transform on the vertical plane (>= 10 ms per track)
1-2. Straight-line fitting on the vertical plane

2. Helix fitting (track chisqr minimization in the 3-D)

*/

//#include <chrono>
#include "TPCLocalTrackHelix.hh"

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <stdexcept>
#include <sstream>
#include <std_ostream.hh>

#include "TH2D.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TF2.h"
#include "TMath.h"
#include "TROOT.h"

#include "TPCPadHelper.hh"
#include "Kinematics.hh"
#include "HoughTransform.hh"
#include "ConfMan.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "PrintHelper.hh"
#include "UserParamMan.hh"

#define DebugDisp 0
#define IterativeResolution 1
//Byungmin's new method for residual definition
#define PtPlane 0
#define CircCross 0

namespace
{
  const auto& gUser = UserParamMan::GetInstance();

  // B-field
  const Double_t& HS_field_0 = ConfMan::Get<Double_t>("HSFLDCALIB");
  const Double_t& HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
  const Double_t& HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");

  const Int_t ReservedNumOfHits = 32*10;
  const Double_t MaxGapBtwClusters = 100.;
  const Int_t MaxIteration = 100;
  //const Int_t MaxIteration = 500; //ref
  const Double_t TruncatedMean = 0.8; //80%

  //for minimization
  static Int_t gNumOfHits;
  static std::vector<TVector3> gHitPos;
  static std::vector<TVector3> gRes;
  static std::vector<Double_t> gHelixTheta;
  static std::vector<Int_t> gLayer;
  static std::vector<Double_t> gPadTheta;
  static std::vector<std::vector<Double_t>> gResParam;
  static TVector3 gVertex;
  static TVector3 gVertexRes;
  static Double_t gPar[5] = {0};
  static Double_t gChisqr = 1.e+10;
  static Bool_t gMomConstraint = false; //w or w/o momentum constraint
  static Bool_t gMultiLoop = false;
  static Int_t gBadHits;
  static Int_t gMinuitStatus;
  const Double_t MaxChisqr = 400.;
  //const Double_t MaxChisqr = 1000.;
  //const Int_t MaxTryMinuit = 0;
  const Int_t MaxTryMinuit = 3;

  //cx, cy, z0, r, dz
  //const Double_t FitStep[5] = { 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-5 };
  const Double_t FitStep[5] = {0.1, 0.1, 0.1, 0.1, 0.0001};
  const Double_t LowLimitBeam[5] = {-50000., -50000., -1000., 5666., -15. };// about 1.7 GeV/c
  const Double_t LowLimit[5] = { -50000., -50000., -15000., 0., -15. };
  //const Double_t UpLimit[5] = { 50000., 50000., 15000., 12000., 15. }; //3.6 GeV/c
  const Double_t UpLimit[5] = { 50000., 50000., 15000., 40000., 15. }; //12.0 GeV/c

  //Window
  //Add hits into the track within the window.
  //(ResidualWindowPull > residual/resolution), (ResidualWindowXZ > residualXZ)
  //const Double_t ResidualWindowPullXZ = 10.; //bad
  const Double_t ResidualWindowPullXZ = 6.; //ref
  const Double_t ResidualWindowPullY = 6.; //ref

  const Double_t ResidualWindowUnderTgtXZ = 5; //[mm]
  const Double_t ResidualWindowInXZ = 10; //[mm]
  const Double_t ResidualWindowOutXZ = 10; //[mm]

  //for good hit selection, hypot(pull_t,pull_y) < PullWindow
  const Double_t PullWindow = 3;

  //For theta calculation
  const Double_t ThetaWindow = 10; //[mm]
  const Double_t ThetaNSigma = 5;

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

  //for Helix tracking
  //[0]~[4] are the Helix parameters,
  //([5],[6],[7]) = (x, y, z)
  static std::string s_tmp="pow([5]-([0]+([3]*cos(x))),2)+pow([6]-([1]+([3]*sin(x))),2)+pow([7]-([2]+([3]*[4]*x)),2)";
  static TF1 fint("fint", s_tmp.c_str(), -10.*TMath::Pi(), 10.*TMath::Pi());
  static std::string s_tmpXZ="pow([2]-([0]+([3]*cos(x))),2)+pow([4]-([1]+([3]*sin(x))),2)";
  static TF1 fintXZ("fintXZ", s_tmpXZ.c_str(), -10.*TMath::Pi(), 10.*TMath::Pi());

  static std::string circ_cross="pow([2]*[2]+[2]*([0]*cos(x)+[1]*sin(x)) +[0]*[0]+[1]*[1]-[3]*[3] ,2)";
  static TF1 fcir_cross("fcir_cross", circ_cross.c_str(), -10.*TMath::Pi(), 10.*TMath::Pi());

}

//______________________________________________________________________________
static inline Bool_t CompareY(const Int_t a, const Int_t b){
  return gHitPos[a].y() < gHitPos[b].y();
}

//______________________________________________________________________________
static inline Bool_t CompareTheta(const Int_t a, const Int_t b){
  return gHelixTheta[a] < gHelixTheta[b];
}

//______________________________________________________________________________
static inline TVector3 GlobalToLocal(TVector3 pos){
  return TVector3(-pos.X(), pos.Z() - tpc::ZTarget, pos.Y());
}

//______________________________________________________________________________
static inline TVector3 LocalToGlobal(TVector3 pos){
  return TVector3(-pos.X(), pos.Z(), pos.Y() + tpc::ZTarget);
}

//______________________________________________________________________________
static inline TVector3 LocalPosition(Double_t par[5], Double_t t){

  //TPC local coordinate
  //This is the eqation of Helix
  Double_t x = par[0] + par[3]*cos(t);
  Double_t y = par[1] + par[3]*sin(t);
  Double_t z = par[2] + (par[4]*par[3]*t);
  return TVector3(x, y, z);
}

//______________________________________________________________________________
static inline TVector3 GlobalPosition(Double_t par[5], Double_t t){

  TVector3 pos = LocalPosition(par, t);
  return LocalToGlobal(pos);
}

//______________________________________________________________________________
static inline Double_t EvalTheta(Double_t par[5], TVector3 pos, Double_t window_low, Double_t window_up){

  fint.SetRange(window_low, window_up);
  Double_t fpar[8];
  TVector3 localpos = GlobalToLocal(pos);
  for(Int_t ip=0; ip<5; ++ip){
    fpar[ip] = par[ip];
  }
  fpar[5] = localpos.X();
  fpar[6] = localpos.Y();
  fpar[7] = localpos.Z();

  Int_t steps = 1440;
  //Int_t steps = TMath::Min((window_up-window_low)*par[3]*10., 200.);
  fint.SetParameters(fpar);
  fint.SetNpx(steps);
  Double_t min_t = fint.GetMinimumX();

#if CircCross
  if(hypot(localpos.X(),localpos.Z())<105) return min_t;
#if PtPlane
  min_t = (window_low+window_up)/2;
#endif
  fcir_cross.SetRange(window_low,window_up);
  double cpar[4];
  cpar[0] = par[0];
  cpar[1] = par[1];
  cpar[2] = par[3];
  cpar[3] = hypot(localpos.X(),localpos.Y());
  fcir_cross.SetParameters(cpar);
  fcir_cross.SetNpx(steps);
  double min_t_temp = fcir_cross.GetMinimumX();
  if(min_t_temp == window_low or min_t_temp == window_up) return min_t;
  else min_t = min_t_temp;
#endif

  return min_t;
}

//______________________________________________________________________________
static inline Double_t EvalThetaXZ(Double_t par[5], TVector3 pos, Double_t window_low, Double_t window_up){

  Double_t fpar[8];
  TVector3 localpos = GlobalToLocal(pos);
  fpar[0] = par[0];
  fpar[1] = par[1];
  fpar[3] = par[3];
  fpar[2] = localpos.X();
  fpar[4] = localpos.Y();

  Int_t steps = 1440;
  //Int_t steps = TMath::Min((window_up-window_low)*par[3]*10., 200.);
  fintXZ.SetRange(window_low, window_up);
  fintXZ.SetParameters(fpar);
  fintXZ.SetNpx(steps);
  Double_t min_t = fintXZ.GetMinimumX();
  return min_t;
}

/*
//______________________________________________________________________________
static inline TVector3 ExpectedPosRow(Double_t par[5], TVector3 pos){ //pad horizontal residual vector

  if(gHelixTheta.size()!=gNumOfHits){
  hddaq::cerr<< "TPCLocalTrackHelix ExpectedPosRow() "<<" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
  }

  TVector3 localpos = GlobalToLocal(pos);
  Double_t x = localpos.x(); Double_t y = localpos.y();
  Double_t theta0 = gHelixTheta[0];
  Double_t k = x*x + y*y - x*par[0] - y*par[1];
  k /= par[3];

  //Two candidates exist
  Double_t cost1 = x*k + y*TMath::Sqrt(TMath::Hypot(x, y) - k*k);
  cost1 /= TMath::Hypot(x, y);
  Double_t sint1 = y*k - x*TMath::Sqrt(TMath::Hypot(x, y) - k*k);
  sint1 /= TMath::Hypot(x, y);
  Double_t x1 = par[0] + par[3]*cost1;
  Double_t y1 = par[1] + par[3]*sint1;
  Double_t theta1 = TMath::ATan2(y1 - par[1], x1 - par[0]);
  if(TMath::Abs(theta1 - theta0) > TMath::Pi()){
  if(theta0>0 && theta1<0) theta1 += 2.*TMath::Pi();
  if(theta0<0 && theta1>0) theta1 -= 2.*TMath::Pi();
  }
  Double_t z1 = par[2] + par[3]*par[4]*theta1;

  Double_t cost2 = x*k - y*TMath::Sqrt(TMath::Hypot(x, y) - k*k);
  cost2 /= TMath::Hypot(x, y);
  Double_t sint2 = y*k + x*TMath::Sqrt(TMath::Hypot(x, y) - k*k);
  sint2 /= TMath::Hypot(x, y);
  Double_t x2 = par[0] + par[3]*cost2;
  Double_t y2 = par[1] + par[3]*sint2;
  Double_t theta2 = TMath::ATan2(y2 - par[1], x2 - par[0]);
  if(TMath::Abs(theta2 - theta0) > TMath::Pi()){
  if(theta0>0 && theta2<0) theta2 += 2.*TMath::Pi();
  if(theta0<0 && theta2>0) theta2 -= 2.*TMath::Pi();
  }
  Double_t z2 = par[2] + par[3]*par[4]*theta2;

  TVector3 calpos_row = TMath::Hypot(x - x1, y - y1) < TMath::Hypot(x - x2, y - y2) ? TVector3(x1, y1, z1) : TVector3(x2, y2, z2);
  return LocalToGlobal(calpos_row);
  }

  //______________________________________________________________________________
  static inline TVector3 ResidualVectRow(Double_t par[5], TVector3 pos){ //pad horizontal residual vector

  TVector3 calpos_row = ExpectedPosRow(par, pos);
  return pos - calpos_row;
  }
*/

//______________________________________________________________________________
static inline TVector3 ResidualVect(Double_t par[5], TVector3 pos, double theta){ //Closest distance on

  //for Helix tracking
  //[0]~[4] are the Helix parameters,
  //([5],[6],[7]) = (x, y, z)
  Double_t window = ThetaWindow/TMath::Min(par[3], 7000.); //p < 2.1 GeV/c
  Double_t evaltheta = EvalTheta(par, pos, theta - 0.5*window, theta + 0.5*window);
  TVector3 calpos = GlobalPosition(par, evaltheta);
  TVector3 d = pos - calpos;
  return d;
}

//______________________________________________________________________________
static inline TVector3 ResidualVect(Double_t par[5], TVector3 pos, Double_t theta_min, Double_t theta_max){ //Closest distance on

  //for Helix tracking
  //[0]~[4] are the Helix parameters,
  //([5],[6],[7]) = (x, y, z)
  Double_t theta = EvalTheta(par, pos, theta_min, theta_max);
  TVector3 calpos = GlobalPosition(par, theta);
  TVector3 d = pos - calpos;
  //std::cout<<"reidual pos "<<pos<<" calpos "<<calpos<<std::endl;
  return d;
}

//______________________________________________________________________________
static inline TVector3 ResidualVectXZ(Double_t par[5], TVector3 pos){ //Closest distance on the y=pos.y() plane

  TVector3 localpos = GlobalToLocal(pos);
  Double_t x = localpos.x(); Double_t y = localpos.y();
  Double_t magnitude = TMath::Hypot(x - par[0], y - par[1]) - par[3];
  TVector3 direction(x - par[0], y - par[1], 0.);
  TVector3 resi_local = magnitude*direction.Unit();
  TVector3 resi_global(-resi_local.x(), 0., resi_local.y());
  return resi_global;
}

//______________________________________________________________________________
static inline TVector3 CalcResolution(Double_t par[5], Int_t layer, TVector3 pos, Double_t padTheta, Double_t theta, std::vector<Double_t> resparam, Bool_t vetoBadClusters){

  Double_t cosPad = TMath::Cos(padTheta);
  Double_t sinPad = TMath::Sin(padTheta);
  Double_t tanPad = TMath::Tan(padTheta);
  Double_t padL = tpc::padParameter[layer][tpc::kLength];
  Double_t padRadius = tpc::padParameter[layer][tpc::kRadius];
  TVector3 closestDist2TrackXZ = ResidualVectXZ(par, TVector3(0., 0., tpc::ZTarget));

  //check whether the track is crossing the layer or not
  if(vetoBadClusters && closestDist2TrackXZ.Mag() > padRadius - 0.5*padL &&
     closestDist2TrackXZ.Mag() < padRadius + 0.5*padL) return TVector3(1.e+10, 1.e+10, 1.e+10);

  //alpha : pad - track angle
  TVector3 localpos = GlobalToLocal(pos);
  Double_t tanTrack = (localpos.y()-par[1])/(localpos.x()-par[0]);
  Double_t tanDiff = (tanPad-tanTrack)/(1.+tanPad*tanTrack);
  Double_t alpha = TMath::ATan(tanDiff);

  //Calculate resolution
  //horizontal resolution
  Double_t param_horizontal[6] = {resparam[0], resparam[1], resparam[2], resparam[3], resparam[4], resparam[5]};
  f_horizontal -> SetParameters(param_horizontal);
  Double_t res_horizontal = f_horizontal -> Eval(alpha, pos.y());

  //vertical resolution
  Double_t param_y[4] = {resparam[6], resparam[1], resparam[7], resparam[8]};
  f_drift -> SetParameters(param_y);
  Double_t res_drift = f_drift -> Eval(pos.y());

  TVector3 res(res_horizontal*TMath::Abs(cosPad), res_drift, res_horizontal*TMath::Abs(sinPad));

  //Residual/resolution < window
  if(vetoBadClusters){
    TVector3 resi = ResidualVect(par, pos, theta);
    Double_t residual_vertical = resi.y();
    Double_t residual_horizontal = TMath::Hypot(resi.x(), resi.z());
    Double_t resolution_vertical = res.y();
    Double_t resolution_horizontal = TMath::Hypot(res.x(), res.z());
    Double_t pull_t = residual_horizontal/resolution_horizontal;
    Double_t pull_y = residual_vertical/resolution_vertical;
    if(TMath::Hypot(pull_t, pull_y) > PullWindow) return TVector3(2.e+10, 2.e+10, 2.e+10);
  }

  return res;
}

//______________________________________________________________________________
static inline void fcn_helix(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  Double_t chisqr=0.; Int_t dof = 0;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 d = ResidualVect(par, gHitPos[i], gHelixTheta[i]);
    if(gRes[i].x() > 0.9e+10 && gRes[i].y() > 0.9e+10 && gRes[i].z() > 0.9e+10) continue; // exclude dummy hits in calculation
    chisqr += TMath::Power(TMath::Hypot(d.x(), d.z())/TMath::Hypot(gRes[i].x(), gRes[i].z()), 2) + TMath::Power(d.y()/gRes[i].y(), 2);
    dof++;
    dof++;
  }
  if(gMomConstraint) dof += 1; //if there is a momentum constraint
  f = chisqr/(Double_t)(dof - 5);
}

//______________________________________________________________________________
static inline void fcn_helixwVertex(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  Double_t chisqr=0.; Int_t dof = 0;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 d = ResidualVect(par, gHitPos[i], gHelixTheta[i]);
    if(gRes[i].x() > 0.9e+10 && gRes[i].y() > 0.9e+10 && gRes[i].z() > 0.9e+10) continue; // exclude dummy hits in calculation
    chisqr += TMath::Power(TMath::Hypot(d.x(), d.z())/TMath::Hypot(gRes[i].x(), gRes[i].z()), 2) + TMath::Power(d.y()/gRes[i].y(), 2);
    dof++;
    dof++;
  }
  TVector3 d = ResidualVect(par, gVertex, -2.*TMath::Pi(), 2.*TMath::Pi());
  TVector3 res = gVertexRes;
  chisqr += TMath::Power(TMath::Hypot(d.x(), d.z())/TMath::Hypot(res.x(), res.z()), 2) + TMath::Power(d.y()/res.y(), 2);
  dof++;
  dof++;

  f = chisqr/(Double_t)(dof - 5);
}

//______________________________________________________________________________
static inline void fcn_line(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  Double_t chisqr = 0.;
  Int_t dof = 0;

  Bool_t flipcheck = false;
  Double_t prev_theta = 0; Double_t theta0 = 0;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 pos = GlobalToLocal(gHitPos[i]);
    //Check ATan2 function's theta flip (-pi ~ pi)
    Double_t tmp_theta = TMath::ATan2(pos.Y() - gPar[1], pos.X() - gPar[0]);
    if(i==0) theta0 = tmp_theta;
    if(TMath::Abs(prev_theta - tmp_theta) > TMath::Pi()) flipcheck = true;
    if(flipcheck){
      if(theta0>0 && tmp_theta<0) tmp_theta += 2.*TMath::Pi();
      if(theta0<0 && tmp_theta>0) tmp_theta -= 2.*TMath::Pi();
    }
    prev_theta = tmp_theta;
    Double_t diff_z = TMath::Power(pos.Z()-(par[0]+(par[1]*gPar[3]*tmp_theta)), 2);
    /*
      Double_t pitch = 2.*TMath::Pi()*par[1]*gPar[3];
      Double_t turns = TMath::Nint((gHitPos[i].y() - ref_ypos)/pitch);
      TVector3 pos = GlobalToLocal(gHitPos[i]);
      //Check ATan2 function's theta flip (-pi ~ pi)
      Double_t tmp_theta = TMath::ATan2(pos.Y() - gPar[1], pos.X() - gPar[0]);
      Double_t param[5] = {gPar[0], gPar[1], par[0], gPar[3], par[1]};
      Double_t ref_ypos = par[0]+(par[1]*gPar[3]*tmp_theta);
      tmp_theta += 2.*TMath::Pi()*turns;
      Double_t diff_z = TMath::Power(pos.Z()-(par[0]+(par[1]*gPar[3]*tmp_theta)), 2);
    */
    /*
      Double_t param[5] = {gPar[0], gPar[1], par[0], gPar[3], par[1]};
      Double_t diff_z = ResidualVect(param, gHitPos[i], gHelixTheta[i]).y();
    */
    chisqr += diff_z/pow(gRes[i].y(),2);
    dof++;
  }
  if(gMomConstraint) dof += 1; //if there is a momentum constraint
  f = chisqr/(Double_t)(dof - 2);
}

//______________________________________________________________________________
static inline Double_t CalcChi2(Double_t *HelixPar, Int_t &ndf, Bool_t vetoBadClusters)
{

  if(gHitPos.size()!=gNumOfHits || gHelixTheta.size()!=gNumOfHits || gLayer.size()!=gNumOfHits || gPadTheta.size()!=gNumOfHits || gResParam.size()!=gNumOfHits){
    hddaq::cerr << "TPCLocalTrackHelix CalcChi2() "
		<< "#hits != # params"
		<< " # hits " << gHitPos.size()
		<< " # layer ids " << gLayer.size()
		<< " # pad angles " << gPadTheta.size()
		<< " # helix thetas " << gHelixTheta.size()
		<< " # res params " << gResParam.size()
		<< std::endl;
    return TMath::QuietNaN();
  }
  //std::cout<<"circle fit"<<std::endl;
  ndf = 0; Double_t chisqr = 0.;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 d = ResidualVect(HelixPar, gHitPos[i], gHelixTheta[i]);
    TVector3 res = CalcResolution(HelixPar, gLayer[i], gHitPos[i], gPadTheta[i], gHelixTheta[i], gResParam[i], vetoBadClusters);
    if(res.x() > 0.9e+10 && res.y() > 0.9e+10 && res.z() > 0.9e+10) continue; // exclude bad clusters
    chisqr += TMath::Power(TMath::Hypot(d.x(), d.z())/TMath::Hypot(res.x(), res.z()), 2) + TMath::Power(d.y()/res.y(), 2);
    ndf++;
    ndf++;
  }
  if(gMomConstraint) ndf += 1; //if there is a momentum constraint
  if(ndf < 6) return 1.e+10;
  return chisqr/(Double_t)(ndf-5);
}

//______________________________________________________________________________
static inline Double_t CircleFit(const Double_t *mX, const Double_t *mY, const Int_t npoints, Double_t* mXCenter, Double_t* mYCenter, Double_t* mRadius)
{

  if(npoints <= 3){
    fprintf(stderr,"CircleFit: npoints %d <= 3\n",npoints);
    return -1;
  }
  else if(npoints > 499){
    fprintf(stderr,"CircleFit: npoints %d > 499\n",npoints);
    return -1;
  }

  //Compute the center of gravity(centroid)
  Double_t xgravity = 0.0;
  Double_t ygravity = 0.0;
  for(Int_t i=0; i<npoints; i++){
    xgravity += mX[i];
    ygravity += mY[i];
  }
  xgravity /= npoints;
  ygravity /= npoints;

  //data matrix components
  Double_t x2 = 0., y2 = 0., xy = 0.;
  Double_t xz = 0., yz = 0., zz = 0.;
  for(Int_t i=0; i<npoints; i++) {
    //shift the coordinates(centroid is at the origin)
    Double_t X = mX[i] - xgravity, Y = mY[i] - ygravity;
    Double_t X2 = X*X, Y2 = Y*Y;
    Double_t Z = X2 + Y2;

    x2 += X2;
    y2 += Y2;
    xy += X*Y;
    xz += X*Z;
    yz += Y*Z;
    zz += Z*Z;
  }

  if(TMath::Abs(x2)<0.0001 || TMath::Abs(y2)<0.0001 || TMath::Abs(xy)<0.0001){
    fprintf(stderr, "CircleFit: x2 = %f, y2 = %f, xy = %f, grav=%f, %f\n", x2, y2, xy, xgravity, ygravity);
    return -1;
  }

  Double_t f, g, h, p, q, t, g0, g02, a, b, c, d;
  Double_t xroot, ff, fp, xd, yd, g1;
  Double_t dx, dy, dradius2, xnom;
  f = (3.*x2+y2)/npoints;
  g = (x2+3.*y2)/npoints;
  h = 2*xy/npoints;
  p = xz/npoints;
  q = yz/npoints;
  t = zz/npoints;
  g0 = (x2+y2)/npoints;
  g02 = g0*g0;
  a = -4.0;
  b = (f*g-t-h*h)/g02;
  c = (t*(f+g)-2.*(p*p+q*q))/(g02*g0);
  d = (t*(h*h-f*g)+2.*(p*p*g+q*q*f)-4.*p*q*h)/(g02*g02);
  xroot = 1.0;
  for(Int_t i=0; i<5; i++){
    ff = (((xroot+a)*xroot+b)*xroot+c)*xroot+d;
    fp = ((4.*xroot+3.*a)*xroot+2.*b)*xroot+c;
    xroot -= ff/fp;
  }
  g1 = xroot*g0;
  xnom = (g-g1)*(f-g1)-h*h;
  if(TMath::Abs(xnom)<0.0001 || TMath::IsNaN(xnom)){
    fprintf(stderr,"CircleFit: xnom1 = %f\n",xnom);
    return -1;
  }

  yd = (q*(f-g1)-h*p)/xnom;
  xnom = f-g1;
  if(TMath::Abs(xnom)<0.0001 || TMath::IsNaN(xnom)){
    fprintf(stderr,"CircleFit: xnom2 = %f\n",xnom);
    return -1;
  }

  xd = (p-h*yd )/xnom;

  Double_t radius2 = xd*xd+yd*yd+g1;
  *mXCenter = xd+xgravity;
  *mYCenter = yd+ygravity;

  Double_t mVariance = 0.0;
  for(Int_t i=0; i<npoints; i++){
    dx = mX[i]-(*mXCenter);
    dy = mY[i]-(*mYCenter);
    dradius2 = dx*dx+dy*dy;
    mVariance += dradius2+radius2-2.*sqrt(dradius2*radius2);
  }

  *mRadius = (Double_t) sqrt(radius2);
  if(mVariance<0){
    fprintf(stderr,"CircleFit: variance = %f\n", mVariance);
    return -1;
  }

#if DebugDisp
  std::cout<<"CircleFit fitting :"
	   <<" variance: "<<mVariance
	   <<", radius: "<<*mRadius
	   <<std::endl;
#endif

  return mVariance;
}

//______________________________________________________________________________
static inline Bool_t StraightLineFit()
{

  if(gHitPos.size()!=gNumOfHits || gHelixTheta.size()!=gNumOfHits || gLayer.size()!=gNumOfHits || gPadTheta.size()!=gNumOfHits || gResParam.size()!=gNumOfHits){
    hddaq::cerr << "TPCLocalTrackHelix StraightLineFit() "
		<< "#hits != # params"
		<< " # hits " << gHitPos.size()
		<< " # layer ids " << gLayer.size()
		<< " # pad angles " << gPadTheta.size()
		<< " # helix thetas " << gHelixTheta.size()
		<< " # res params " << gResParam.size()
		<< std::endl;
    return false;
  }

  //pre t-y fit
  Double_t par_li[2] = {gPar[2], gPar[4]};
  Double_t err_li[2] = {-999., -999.};
  TMinuit *minuit = new TMinuit(2);
  minuit->SetPrintLevel(-1);
  minuit->SetFCN(fcn_line);

  Int_t ierflg_li = 0;
  Double_t arglist_li[10];
  arglist_li[0] = 2.3;
  minuit->mnexcm("SET ERR", arglist_li,1,ierflg_li); // No warnings
  arglist_li[0] = 1;
  minuit->mnexcm("SET NOW", arglist_li,1,ierflg_li); // No warnings

  TString name_li[2] = {"z0", "dz"};
  minuit->mnparm(0, name_li[0], par_li[0], FitStep[2], LowLimit[2], UpLimit[2], ierflg_li);
  minuit->mnparm(1, name_li[1], par_li[1], FitStep[4], LowLimit[4], UpLimit[4], ierflg_li);

  minuit->Command("SET STRategy 0");
  arglist_li[0] = 1000*5*5*5;
  arglist_li[1] = 0.1/(10.*10.*10.);
  minuit->mnexcm("MIGRAD", arglist_li, 2, ierflg_li);

  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

  Int_t Err_li;
  Double_t bnd1_li, bnd2_li;
  for(Int_t i=0; i<2; i++){
    minuit->mnpout(i, name_li[i], par_li[i], err_li[i], bnd1_li, bnd2_li, Err_li);
  }

  //Double_t grad[5]; Double_t Chisqr;
  //minuit -> Eval(5, grad, Chisqr, par_li, 0);
  delete minuit;
#if DebugDisp
  if(icstat==0) std::cout<<"StraightLineFit() icstat=="<<icstat<<std::endl;
#endif

  gPar[2] = par_li[0];
  gPar[4] = par_li[1];
  gMinuitStatus = icstat;
  return true;
}

//______________________________________________________________________________
static inline Bool_t HelixFit(Int_t IsBeam, Bool_t vetoBadClusters, Bool_t ExclusiveFlag = false){

  if(gHitPos.size()!=gNumOfHits || gHelixTheta.size()!=gNumOfHits || gLayer.size()!=gNumOfHits || gPadTheta.size()!=gNumOfHits || gResParam.size()!=gNumOfHits){
    hddaq::cerr << " TPCLocalTrackHelix HelixFit() "
		<< "#hits != # params"
		<< " # hits " << gHitPos.size()
		<< " # layer ids " << gLayer.size()
		<< " # pad angles " << gPadTheta.size()
      		<< " # helix thetas " << gHelixTheta.size()
		<< " # res params " << gResParam.size()
		<< std::endl;
    return false;
  }

  Double_t par[5] = {gPar[0], gPar[1], gPar[2], gPar[3], gPar[4]};
  Double_t err[5] = {-999., -999., -999., -999., -999.};

  gRes.clear();
  for(Int_t i=0; i<gNumOfHits; i++){
    TVector3 res = CalcResolution(gPar, gLayer[i], gHitPos[i], gPadTheta[i], gHelixTheta[i], gResParam[i], vetoBadClusters);
    gRes.push_back(res);
  }

  //NDF value becomes different. Initialize gChisqr and fit again.
  Int_t ndf;
  if(vetoBadClusters){
    gChisqr = CalcChi2(gPar, ndf, vetoBadClusters);
    gBadHits = gNumOfHits - 0.5*ndf;
    if(TMath::Abs(gChisqr-1.e+10)<0.1) return false;
  }

  //for exclusive tracking, previous fitting result should not affect current fitting
  if(ExclusiveFlag) gChisqr = 1e+10;

  TMinuit *minuit = new TMinuit(5);
  minuit->SetPrintLevel(-1);
  minuit->SetFCN(fcn_helix);

  Int_t ierflg = 0;
  Double_t arglist[10];
  arglist[0] = 5.89;
  minuit->mnexcm("SET ERR", arglist, 1, ierflg); //Num of parameter
  arglist[0] = 1;
  minuit->mnexcm("SET NOW", arglist, 1, ierflg); // No warnings

  TString name[5] = {"cx", "cy", "z0", "r", "dz"};
  for(Int_t i=0; i<5; i++){
    if(IsBeam==1) minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimitBeam[i], UpLimit[i], ierflg);
    else minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }

  if(gMomConstraint) minuit -> FixParameter(3);

  minuit->Command("SET STRategy 0");
  //arglist[0] = 5000.;
  //arglist[1] = 0.01;
  arglist[0] = 1000.;
  arglist[1] = 0.1;

  Int_t Err;
  Double_t bnd1, bnd2;

  Bool_t status = false;
  Int_t itry=0; gMinuitStatus = 0;
  Double_t good_chisqr = 1.5;
  while(!status){
    if(itry>MaxTryMinuit) break;
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    minuit->mnimpr();
    //minuit->mnexcm("MINOS", arglist, 0, ierflg);
    //minuit->mnexcm("SET ERR", arglist, 2, ierflg);

    Double_t amin, edm, errdef;
    Int_t nvpar, nparx, icstat;
    minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

    //minuit->mnprin(4, amin);
    for(Int_t i=0; i<5; i++){
      minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
    }

    //Double_t grad[5];
    //minuit -> Eval(5, grad, Chisqr, par, 0);
    Double_t Chisqr = CalcChi2(par, ndf, vetoBadClusters);
    if(gChisqr>=Chisqr || TMath::Abs(gChisqr-Chisqr) < 0.01){
      gChisqr = Chisqr;
      gPar[0] = par[0];
      gPar[1] = par[1];
      gPar[2] = par[2];
      gPar[3] = par[3];
      gPar[4] = par[4];
      gMinuitStatus = icstat;
      gBadHits = gNumOfHits - 0.5*ndf;
      status = true;
    }
    arglist[0] = arglist[0]*5;
    arglist[1] = arglist[1]*0.1;
    ++itry;
  }
  delete minuit;
  if(IsBeam==1 && !status && gChisqr < 4.5) status = true;
  if(gChisqr < good_chisqr) status = true;

#if DebugDisp
  std::cout<<"HelixFit() status="<<status<<" gChisqr "<<gChisqr<<" isBeam "<<IsBeam<<std::endl;
  if(gMinuitStatus==0) std::cout<<"HelixFit() icstat==0"<<std::endl;
#endif
  return status;
}

//______________________________________________________________________________
static inline Bool_t HelixFitInvertCharge(){

  if(gHitPos.size()!=gNumOfHits || gHelixTheta.size()!=gNumOfHits || gLayer.size()!=gNumOfHits || gPadTheta.size()!=gNumOfHits || gResParam.size()!=gNumOfHits || gRes.size()!=gNumOfHits){
    hddaq::cerr << " TPCLocalTrackHelix HelixFit() "
		<< "#hits != # params"
		<< " # hits " << gHitPos.size()
		<< " # layer ids " << gLayer.size()
		<< " # pad angles " << gPadTheta.size()
      		<< " # helix thetas " << gHelixTheta.size()
		<< " # res params " << gResParam.size()
		<< " # res " << gRes.size()
		<< std::endl;
    return false;
  }

  Double_t par[5] = {gPar[0], gPar[1], gPar[2], gPar[3], gPar[4]};
  Double_t err[5] = {-999., -999., -999., -999., -999.};

  Double_t lowLimit[5] = { gPar[0] - 0.5*gPar[3], gPar[1] - 0.5*gPar[3], gPar[2] - 5000., 0, gPar[4] - 5.};
  Double_t upLimit[5] = { gPar[0] + 0.5*gPar[3], gPar[1] + 0.5*gPar[3], gPar[2] + 5000., 12000, gPar[4] + 5.}; //12.0 GeV/c

  TMinuit *minuit = new TMinuit(5);
  minuit->SetPrintLevel(-1);
  minuit->SetFCN(fcn_helix);

  Int_t ierflg = 0;
  Double_t arglist[10];
  arglist[0] = 5.89;
  minuit->mnexcm("SET ERR", arglist, 1, ierflg); //Num of parameter
  arglist[0] = 1;
  minuit->mnexcm("SET NOW", arglist, 1, ierflg); // No warnings

  TString name[5] = {"cx", "cy", "z0", "r", "dz"};
  for(Int_t i=0; i<5; i++){
    minuit->mnparm(i, name[i], par[i], FitStep[i], lowLimit[i], upLimit[i], ierflg);
  }

  minuit->Command("SET STRategy 0");
  //arglist[0] = 5000.;
  //arglist[1] = 0.01;
  arglist[0] = 1000.;
  arglist[1] = 0.1;

  Int_t Err;
  Double_t bnd1, bnd2;
  Int_t ndf; Int_t gMinuitStatus = 0;

  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  minuit->mnimpr();
  //minuit->mnexcm("MINOS", arglist, 0, ierflg);
  //minuit->mnexcm("SET ERR", arglist, 2, ierflg);

  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  //minuit->mnprin(4, amin);
  for(Int_t i=0; i<5; i++){
    minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
  }

  gPar[0] = par[0];
  gPar[1] = par[1];
  gPar[2] = par[2];
  gPar[3] = par[3];
  gPar[4] = par[4];
  gChisqr = CalcChi2(gPar, ndf, false);

  Int_t itry=0;
  Double_t good_chisqr = 1.5;
  Bool_t status = false;
  while(!status){
    if(itry>MaxTryMinuit) break;
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    minuit->mnimpr();
    //minuit->mnexcm("MINOS", arglist, 0, ierflg);
    //minuit->mnexcm("SET ERR", arglist, 2, ierflg);

    Double_t amin, edm, errdef;
    Int_t nvpar, nparx, icstat;
    minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

    //minuit->mnprin(4, amin);
    for(Int_t i=0; i<5; i++){
      minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
    }

    //Double_t grad[5];
    //minuit -> Eval(5, grad, Chisqr, par, 0);
    Double_t Chisqr = CalcChi2(par, ndf, false);
    if(gChisqr>=Chisqr || TMath::Abs(gChisqr-Chisqr) < 0.01){
      gChisqr = Chisqr;
      gPar[0] = par[0];
      gPar[1] = par[1];
      gPar[2] = par[2];
      gPar[3] = par[3];
      gPar[4] = par[4];
      gMinuitStatus = icstat;
      gBadHits = gNumOfHits - 0.5*ndf;
      status = true;
    }
    arglist[0] = arglist[0]*5;
    arglist[1] = arglist[1]*0.1;
    ++itry;
  }
  delete minuit;

#if DebugDisp
  std::cout<<"HelixFitInvertCharge() status="<<status<<" gChisqr "<<gChisqr<<std::endl;
  if(gMinuitStatus==0) std::cout<<"HelixFit() icstat==0"<<std::endl;
#endif
  return true;
}

//______________________________________________________________________________
static inline Double_t CalcChi2wVertex(Double_t *HelixPar, Int_t &ndf)
{

  if(gHitPos.size()!=gNumOfHits || gHelixTheta.size()!=gNumOfHits || gLayer.size()!=gNumOfHits || gPadTheta.size()!=gNumOfHits || gResParam.size()!=gNumOfHits){
    hddaq::cerr << "TPCLocalTrackHelix CalcChi2wVertex() "
		<< "#hits != # params"
		<< " # hits " << gHitPos.size()
		<< " # layer ids " << gLayer.size()
		<< " # pad angles " << gPadTheta.size()
		<< " # helix thetas " << gHelixTheta.size()
		<< " # res params " << gResParam.size()
		<< std::endl;
    return TMath::QuietNaN();
  }

  ndf = 0; Double_t chisqr = 0.;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 d = ResidualVect(HelixPar, gHitPos[i], gHelixTheta[i]);
    TVector3 res = gRes[i];
    if(res.x() > 0.9e+10 && res.y() > 0.9e+10 && res.z() > 0.9e+10) continue; // exclude bad clusters
    chisqr += TMath::Power(TMath::Hypot(d.x(), d.z())/TMath::Hypot(res.x(), res.z()), 2) + TMath::Power(d.y()/res.y(), 2);
    ndf++;
    ndf++;
  }
  TVector3 d = ResidualVect(HelixPar, gVertex, -2.*TMath::Pi(), 2.*TMath::Pi());
  TVector3 res = gVertexRes;
  chisqr += TMath::Power(TMath::Hypot(d.x(), d.z())/TMath::Hypot(res.x(), res.z()), 2) + TMath::Power(d.y()/res.y(), 2);
  ndf++;
  ndf++;

  if(ndf < 6) return 1.e+10;
  return chisqr/(Double_t)(ndf-5);
}

//______________________________________________________________________________
static inline Bool_t HelixFitwVertex(){

  if(gHitPos.size()!=gNumOfHits || gHelixTheta.size()!=gNumOfHits || gLayer.size()!=gNumOfHits || gPadTheta.size()!=gNumOfHits || gResParam.size()!=gNumOfHits){
    hddaq::cerr << " TPCLocalTrackHelix HelixFitwVertex() "
		<< "#hits != # params"
		<< " # hits " << gHitPos.size()
		<< " # layer ids " << gLayer.size()
		<< " # pad angles " << gPadTheta.size()
      		<< " # helix thetas " << gHelixTheta.size()
		<< " # res params " << gResParam.size()
		<< std::endl;
    return false;
  }

  Double_t par[5] = {gPar[0], gPar[1], gPar[2], gPar[3], gPar[4]};
  Double_t err[5] = {-999., -999., -999., -999., -999.};

  gChisqr = 1e+10;

  TMinuit *minuit = new TMinuit(5);
  minuit->SetPrintLevel(-1);
  minuit->SetFCN(fcn_helixwVertex);

  Int_t ierflg = 0;
  Double_t arglist[10];
  arglist[0] = 5.89;
  minuit->mnexcm("SET ERR", arglist, 1, ierflg); //Num of parameter
  arglist[0] = 1;
  minuit->mnexcm("SET NOW", arglist, 1, ierflg); // No warnings

  TString name[5] = {"cx", "cy", "z0", "r", "dz"};
  for(Int_t i=0; i<5; i++){
    minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }

  minuit->Command("SET STRategy 0");
  //arglist[0] = 5000.;
  //arglist[1] = 0.01;
  arglist[0] = 1000.;
  arglist[1] = 0.1;

  Int_t Err;
  Double_t bnd1, bnd2;

  Bool_t status = false;
  Int_t ndf=0; Int_t itry=0; gMinuitStatus = 0;
  while(true){
    if(itry>MaxTryMinuit) break;
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    minuit->mnimpr();
    //minuit->mnexcm("MINOS", arglist, 0, ierflg);
    //minuit->mnexcm("SET ERR", arglist, 2, ierflg);

    Double_t amin, edm, errdef;
    Int_t nvpar, nparx, icstat;
    minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

    //minuit->mnprin(4, amin);
    for(Int_t i=0; i<5; i++){
      minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
    }

    //Double_t grad[5];
    //minuit -> Eval(5, grad, Chisqr, par, 0);
    Double_t Chisqr = CalcChi2wVertex(par, ndf);
    if(gChisqr>=Chisqr || TMath::Abs(gChisqr-Chisqr) < 0.01){
      gChisqr = Chisqr;
      gPar[0] = par[0];
      gPar[1] = par[1];
      gPar[2] = par[2];
      gPar[3] = par[3];
      gPar[4] = par[4];
      gMinuitStatus = icstat;
      //gBadHits = gNumOfHits - 0.5*ndf;
      gBadHits = gNumOfHits - 0.5*ndf + 1; // +1 by considering the vertex point
    }
    arglist[0] = arglist[0]*5;
    arglist[1] = arglist[1]*0.1;
    ++itry;
  }
  delete minuit;

#if DebugDisp
  std::cout<<"HelixFitwVertex() status="<<status<<" gChisqr "<<gChisqr<<std::endl;
  if(gMinuitStatus==0) std::cout<<"HelixFit() icstat==0"<<std::endl;
#endif

  return true;
}

//Only for momentum constraint fitting
//______________________________________________________________________________
static inline void fcn_circle(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  Double_t temp = 0.;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 pos = GlobalToLocal(gHitPos[i]);

    TVector3 localpos(pos.x(), pos.y(), 0.);
    TVector3 center(par[0], par[1], 0.); //Circle center
    TVector3 shifted_pos = localpos - center;
    //temp += TMath::Power(shifted_pos.Mag2() - par[2]*par[2], 2); //(x^2 + y^2 - R^2)^2
    temp += TMath::Power(shifted_pos.Mag() - par[2], 2); //(x^2 + y^2 - R^2)
  }
  f = temp;
}

//______________________________________________________________________________
static inline Double_t CircleFit(Int_t IsBeam)
{

  if(gHitPos.size()!=gNumOfHits){
    hddaq::cerr << " TPCLocalTrackHelix CircleFit() "
		<< "#hits != # params"
		<< " # hits " << gHitPos.size()
		<< std::endl;
    return -1;
  }

  Double_t par_circ[3] = {gPar[0], gPar[1], gPar[3]};
  Double_t err_circ[3] = {-999., -999., -999.};
  //std::cout<<"circlefit "<<par_circ[0]<<" "<<par_circ[1]<<" "<<par_circ[2]<<std::endl;
  //cx, cy, r
  Double_t fitStep[3] = {1.0e-4, 1.0e-4, 1.0e-4};
  Double_t lowLimitBeam[3] = {gPar[0] - 500., gPar[1] - 500., 5666.};// about 1.7 GeV/c
  Double_t upLimitBeam[3] = {gPar[0] + 500., gPar[1] + 500., 6333.};// about 1.9 GeV/c
  //Double_t lowLimit[3] = {gPar[0] - 500., gPar[1] - 500., 0.};
  //Double_t upLimit[3] = {gPar[0] + 500., gPar[1] + 500., 30000.}; // ~10 GeV/c
  Double_t lowLimit[3] = {gPar[0] - 500., gPar[1] - 500., gPar[3] - 3000};
  Double_t upLimit[3] = {gPar[0] + 500., gPar[1] + 500., gPar[3] + 3000};

  TMinuit *minuit = new TMinuit(3);
  minuit->SetPrintLevel(-1);
  minuit->SetFCN(fcn_circle);

  Int_t ierflg_circ = 0;
  Double_t arglist_circ[10];
  arglist_circ[0] = 3.52; //for 3 parameters
  minuit->mnexcm("SET ERR", arglist_circ, 1, ierflg_circ);
  arglist_circ[0] = 1;
  minuit->mnexcm("SET NOW", arglist_circ, 1, ierflg_circ); // No warnings

  TString name[3] = {"cx", "cy", "r"};
  for(Int_t i=0; i<3; i++){
    if(IsBeam==1) minuit->mnparm(i, name[i], par_circ[i], fitStep[i], lowLimitBeam[i], upLimitBeam[i], ierflg_circ);
    else minuit->mnparm(i, name[i], par_circ[i], fitStep[i], lowLimit[i], upLimit[i], ierflg_circ);
  }
  if(gMomConstraint) minuit -> FixParameter(2);

  minuit->Command("SET STRategy 0");

  arglist_circ[0] = 1000*5*5*5;
  arglist_circ[1] = 0.1/(10.*10.*10.);
  Int_t Err_circ;
  Double_t bnd1_circ, bnd2_circ;

  minuit->mnexcm("MIGRAD", arglist_circ, 2, ierflg_circ);
  for(Int_t i=0; i<3; i++){
    minuit->mnpout(i, name[i], par_circ[i], err_circ[i], bnd1_circ, bnd2_circ, Err_circ);
  }
  //std::cout<<"circlefit "<<par_circ[0]<<" "<<par_circ[1]<<" "<<par_circ[2]<<std::endl;
  Double_t grad[3];
  Double_t chisqr = 0.;
  minuit -> Eval(3, grad, chisqr, par_circ, 0);
  delete minuit;

  gPar[0] = par_circ[0];
  gPar[1] = par_circ[1];
  gPar[3] = par_circ[2];
  //std::cout<<"circlefit "<<par_circ[0]<<" "<<par_circ[1]<<" "<<par_circ[2]<<std::endl;
  return chisqr;
}

//______________________________________________________________________________
TPCLocalTrackHelix::TPCLocalTrackHelix()
  : m_is_fitted(false),
    m_is_calculated(false),
    m_is_theta_calculated(false),
    m_is_fitted_exclusive(false),
    m_is_multiloop(false),
    m_hit_order(), m_hit_t(), m_kuramaid_candidate(),
    m_cx(0.), m_cy(0.), m_z0(0.), m_r(0.), m_dz(0.),
    m_pid(0),
    m_closedist(1.e+10, 1.e+10, 1.e+10),
    m_closedistXZ(1.e+10, 1.e+10, 1.e+10),
    m_chisqr(1.e+10),
    m_minuit(0),
    m_n_iteration(0),
    m_mom0(0.,0.,0.),
    m_edgepoint(0.,0.,-143.),
    m_min_t(0.), m_max_t(0.),
    m_path(0.), m_transverse_path(0.),
    m_charge(0), m_fitflag(0), m_vtxflag(0),
    m_isBeam(0), m_isK18(0), m_isKurama(0), m_isAccidental(0), m_isXi(0),
    m_trackid(-1),
    m_ncl_beforetgt(-1),
    m_searchtime(0), m_fittime(0),
    m_MomResScale(-1),
    m_dZResScale(-1),
    m_PhResScale(-1),
    m_cx_exclusive(), m_cy_exclusive(), m_z0_exclusive(),
    m_r_exclusive(), m_dz_exclusive(),
    m_chisqr_exclusive(),
    m_t_exclusive(),
    m_vp()
{

  static const Double_t MomResScale = gUser.GetParameter("MomResScale") ;
  static const Double_t dZResScale = gUser.GetParameter("dZResScale") ;
  static const Double_t PhiResScale = gUser.GetParameter("PhiResScale") ;
  m_MomResScale = MomResScale;
  m_dZResScale = dZResScale;
  m_PhResScale = PhiResScale;

  m_hit_array.reserve(ReservedNumOfHits);
  debug::ObjectCounter::increase(ClassName());
}

//______________________________________________________________________________
TPCLocalTrackHelix::~TPCLocalTrackHelix()
{
  debug::ObjectCounter::decrease(ClassName());
}

//______________________________________________________________________________
TPCLocalTrackHelix::TPCLocalTrackHelix(TPCLocalTrackHelix *init){

  this -> m_is_fitted = false;
  this -> m_is_calculated = false;
  this -> m_is_fitted_exclusive = false;

  for(Int_t i=0;i<init -> m_hit_array.size();i++){
    this -> m_hit_array.push_back(new TPCLTrackHit(init -> m_hit_array[i] -> GetHit()));
  }

  this -> m_cx = init -> m_cx ;
  this -> m_cy = init -> m_cy ;
  this -> m_z0 = init -> m_z0 ;
  this -> m_r = init -> m_r ;
  this -> m_dz = init -> m_dz ;
  this -> m_pid = init -> m_pid ;
  this -> m_closedist = init -> m_closedist ;
  this -> m_closedistXZ = init -> m_closedistXZ ;
  this -> m_chisqr = init -> m_chisqr ;
  this -> m_minuit = init -> m_minuit ;
  this -> m_n_iteration = init -> m_n_iteration ;
  this -> m_mom0 = init -> m_mom0 ;
  this -> m_edgepoint = init -> m_edgepoint ;
  this -> m_min_t = init -> m_min_t ;
  this -> m_max_t = init -> m_max_t ;
  this -> m_path = init -> m_path ;
  this -> m_transverse_path = init -> m_transverse_path ;
  this -> m_charge = init -> m_charge ;
  this -> m_fitflag = init -> m_fitflag ;
  this -> m_vtxflag = init -> m_vtxflag ;
  this -> m_isBeam = init -> m_isBeam ;
  this -> m_isK18 = init -> m_isK18 ;
  this -> m_isKurama = init -> m_isKurama ;
  this -> m_isAccidental = init -> m_isAccidental ;
  this -> m_isXi = init -> m_isXi ;
  this -> m_trackid = init -> m_trackid ;
  this -> m_ncl_beforetgt = init -> m_ncl_beforetgt ;
  this -> m_searchtime = init -> m_searchtime ; //millisec
  this -> m_fittime = init -> m_fittime ; //millisec
  this -> m_is_multiloop = init -> m_is_multiloop ;

  for(Int_t i=0;i<init -> m_hit_order.size();i++){
    this -> m_hit_order.push_back(init -> m_hit_order[i]);
  }

  this -> m_is_theta_calculated = init -> m_is_theta_calculated;

  for(Int_t i=0;i<init -> m_hit_t.size();i++){
    this -> m_hit_t.push_back(init -> m_hit_t[i]);
  }

  for(Int_t i=0;i<init -> m_kuramaid_candidate.size();i++){
    this -> m_kuramaid_candidate.push_back(init -> m_kuramaid_candidate[i]);
  }

  for(Int_t i=0;i<init -> m_cx_exclusive.size();i++){
    this -> m_cx_exclusive.push_back(init -> m_cx_exclusive[i]);
    this -> m_cy_exclusive.push_back(init -> m_cy_exclusive[i]);
    this -> m_z0_exclusive.push_back(init -> m_z0_exclusive[i]);
    this -> m_r_exclusive.push_back(init -> m_r_exclusive[i]);
    this -> m_dz_exclusive.push_back(init -> m_dz_exclusive[i]);
    this -> m_chisqr_exclusive.push_back(init -> m_chisqr_exclusive[i]);
    this -> m_t_exclusive.push_back(init -> m_t_exclusive[i]);

    this -> m_hit_array[i] -> SetCalHelixExclusive(init -> m_cx_exclusive[i],
						   init -> m_cy_exclusive[i],
						   init -> m_z0_exclusive[i],
						   init -> m_r_exclusive[i],
						   init -> m_dz_exclusive[i]);
    this -> m_hit_array[i] -> SetThetaExclusive(init -> m_t_exclusive[i]);
    this -> m_hit_array[i] -> SetCalPositionExclusive(init -> m_hit_array[i] -> GetLocalCalPosHelixExclusive());
  }

  for(Int_t i=0;i<init -> m_vp.size();i++){
    this -> m_vp.push_back(init -> m_vp[i]);
  }

  for(Int_t i=0;i<this -> m_hit_array.size();i++){
    this -> m_hit_array[i] -> SetCalHelix(init -> m_cx, init -> m_cy, init -> m_z0, init -> m_r, init -> m_dz);
    this -> m_hit_array[i] -> SetTheta(init -> m_hit_array[i] -> GetTheta());
    this -> m_hit_array[i] -> SetCalPosition(init -> m_hit_array[i] -> GetLocalCalPosHelix());
    this -> m_hit_array[i] -> SetResolution(init -> m_hit_array[i] -> GetResolutionVect());
  }

  static const Double_t MomResScale = gUser.GetParameter("MomResScale") ;
  static const Double_t dZResScale = gUser.GetParameter("dZResScale") ;
  static const Double_t PhiResScale = gUser.GetParameter("PhiResScale") ;
  m_MomResScale = MomResScale;
  m_dZResScale = dZResScale;
  m_PhResScale = PhiResScale;

  debug::ObjectCounter::increase(ClassName());
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::ClearHits()
{
  m_hit_array.clear();
  m_hit_order.clear();
  m_hit_t.clear();
}

//______________________________________________________________________________
Int_t
TPCLocalTrackHelix::GetNHitsEffective() const
{

  Int_t n = 0;
  for(auto h:m_hit_array){
    if(!h->IsGoodForTracking()) continue;
    ++n;
  }
  return n;
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::AddTPCHit(TPCLTrackHit *hit)
{
  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};

  if(hit->IsGood()){
    m_hit_order.push_back(m_hit_array.size());
    m_hit_array.push_back(hit);
    if(m_is_theta_calculated){
      //Calculate Atan2(y,x) and turns of helix. And convert them into the theta of helix
      TVector3 pos = hit->GetLocalHitPos();
      TVector3 localpos = GlobalToLocal(pos);
      Double_t theta = TMath::ATan2(localpos.y() - m_cy, localpos.x() - m_cx);
      Double_t pitch = 2.*TMath::Pi()*m_r*m_dz;
      Double_t ref_ypos = GetPosition(par, theta).y();
      Double_t turns = TMath::Nint((pos.y() - ref_ypos)/pitch);
      theta += 2.*TMath::Pi()*turns;
      m_hit_t.push_back(theta);
    }
    else m_hit_t.push_back(hit->GetTheta());
  }
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::Calculate()
{
  if(IsCalculated()){
    hddaq::cerr << "#W " << FUNC_NAME << " "
		<< "already called" << std::endl;
    return;
  }

  const std::size_t n = m_hit_array.size();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    hitp->SetCalHelix(m_cx, m_cy, m_z0, m_r, m_dz);
    hitp->SetTheta(m_hit_t[i]);
    hitp->SetCalPosition(hitp->GetLocalCalPosHelix());
    if(m_is_fitted){
      hitp->SetResolution(GetResolutionVect(i, true));
    }
  }

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  m_mom0 = CalcHelixMom(par, 0.);
  m_pid = Kinematics::HypTPCdEdxPID(GetdEdx(TruncatedMean), (Double_t) m_charge*m_mom0.Mag());
  m_is_calculated = true;

}

//______________________________________________________________________________
void
TPCLocalTrackHelix::CalculateExclusive()
{

  if(!IsCalculated()){
    hddaq::cerr << "#W " << FUNC_NAME
		<< "No inclusive calculation" << std::endl;
    return;
  }

  if(!m_is_fitted_exclusive){
    hddaq::cerr << "#W " << FUNC_NAME
		<< "No exclusive fitting" << std::endl;
    return;
  }

  const std::size_t n = m_hit_array.size();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    hitp->SetCalHelixExclusive(m_cx_exclusive[i], m_cy_exclusive[i],
			       m_z0_exclusive[i], m_r_exclusive[i],
			       m_dz_exclusive[i]);
    hitp->SetThetaExclusive(m_t_exclusive[i]);
    hitp->SetCalPositionExclusive(hitp->GetLocalCalPosHelixExclusive());
  }
}

//______________________________________________________________________________
Int_t
TPCLocalTrackHelix::GetNPad() const
{

  // #hits < MinHits
  const std::size_t n = m_hit_array.size();
  std::vector<Int_t> pads;
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hit = m_hit_array[i];
    pads.push_back(hit -> GetHit() -> GetPad());
  }
  std::sort(pads.begin(), pads.end());
  pads.erase(std::unique(pads.begin(), pads.end()), pads.end());
  return pads.size();
}

//______________________________________________________________________________
Int_t
TPCLocalTrackHelix::GetNDF() const
{

  Int_t nhit = GetNHit();
  Int_t ndf = 2*nhit - 5;
  if(gMomConstraint) ndf += 1;
  return ndf;
}

//______________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetPosition(Double_t par[5], Double_t t) const
{
  return GlobalPosition(par, t);
}

//_____________________________________________________________________________
void
TPCLocalTrackHelix::CalcClosestDistTgt()
{

  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return;
  }

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 tgt(0., 0., tpc::ZTarget);
  Double_t scantheta = 0.5*TMath::Pi();
  Double_t theta_tgt = EvalTheta(par, tgt, m_min_t - scantheta, m_max_t + scantheta);
  TVector3 closest_point = GlobalPosition(par, theta_tgt);
  TVector3 dist = closest_point - tgt;
  m_closedist = dist;

  TVector3 residualXZ = ResidualVectXZ(par, tgt);
  TVector3 horizontal_closest_point = residualXZ + tgt;
  Double_t theta_tgtXZ = EvalThetaXZ(par, horizontal_closest_point, m_min_t - scantheta, m_max_t + scantheta);
  TVector3 closest_pointXZ = GlobalPosition(par, theta_tgtXZ);
  TVector3 distXZ = closest_pointXZ - tgt;
  m_closedistXZ = distXZ;
}

//_____________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetResolutionVect(TPCLTrackHit* hit, Bool_t vetoBadClusters){

  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return TVector3(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN());
  }

  if(vetoBadClusters && hit->GetHit()->GetParentCluster()->IsOnTheFrame()) return TVector3(3.e+10, 3.e+10, 3.e+10);

  TVector3 pos = hit->GetLocalHitPos();
  Int_t layer = hit->GetLayer();
  Double_t padTheta = hit->GetPadTheta();
  Double_t theta = hit->GetTheta();
  std::vector<Double_t> resParam = hit->GetResolutionParams();
  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};

  //Convert resolution along the row direction into closets point's x, y, z resolutions
  TVector3 res = CalcResolution(par, layer, pos, padTheta, theta, resParam, vetoBadClusters);
  return res;
}

//_____________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetResolutionVect(Int_t i, Bool_t vetoBadClusters){

  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return TVector3(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN());
  }

  return GetResolutionVect(m_hit_array[i], vetoBadClusters);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetResolutionY(TPCHit* hit){

  //vertical resolution
  std::vector<Double_t> resparam = hit->GetResolutionParams();
  Double_t param_y[4] = {resparam[6], resparam[1], resparam[7], resparam[8]};
  f_drift -> SetParameters(param_y);
  TVector3 pos = hit->GetPosition();
  Double_t res_drift = f_drift -> Eval(pos.y());
  return res_drift;
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetResolutionY(Int_t i){

  TPCHit *hit = m_hit_array[i] -> GetHit();
  return GetResolutionY(hit);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetAlpha(Int_t i) const //for a point on the track
{

  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return -1;
  }

  TPCHit *hit = m_hit_array[i] -> GetHit();
  return GetAlpha(hit);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetAlpha(TPCHit* hit) const
{

  //find expected position on the track
  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 pos = hit->GetPosition();
  TVector3 localpos = GlobalToLocal(pos);
  Double_t tanTrack = (localpos.y()-par[1])/(localpos.x()-par[0]);
  Double_t padTheta = hit->GetPadTheta();
  Double_t tanPad = TMath::Tan(padTheta);

  //alpha : pad - track angle
  Double_t tanDiff = (tanPad-tanTrack)/(1.+tanPad*tanTrack);
  Double_t alpha = TMath::ATan(tanDiff);
  return alpha;
}

//______________________________________________________________________________
TVector3
TPCLocalTrackHelix::CalcResidual(TVector3 pos)
{

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  return ResidualVect(par, pos, -2.*TMath::Pi(), 2.*TMath::Pi());
}

//______________________________________________________________________________
TVector3
TPCLocalTrackHelix::CalcHelixMomCenter(Double_t par[5]) const
{

  Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  Double_t pt = fabs(par[3])*(tpc::ConstC*dMagneticField); // GeV/c
  Double_t theta = -par[2]/(par[3]*par[4]); //at y = 0

  Double_t tmp_px = pt*(-1.*sin(theta));
  Double_t tmp_py = pt*(cos(theta));
  Double_t tmp_pz = pt*(par[4]);
  Double_t px = -tmp_px*0.001;
  Double_t py = tmp_pz*0.001;
  Double_t pz = tmp_py*0.001;

  TVector3 p = TVector3(px,py,pz);
  if(m_charge < 0) p *= -1.;
  return p;
}

//______________________________________________________________________________
TVector3
TPCLocalTrackHelix::CalcHelixMom(Double_t par[5], Double_t theta) const
{

  Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  Double_t pt = fabs(par[3])*(tpc::ConstC*dMagneticField); // GeV/c

  Double_t tmp_px = pt*(-1.*sin(theta));
  Double_t tmp_py = pt*(cos(theta));
  Double_t tmp_pz = pt*(par[4]);
  Double_t px = -tmp_px*0.001;
  Double_t py = tmp_pz*0.001;
  Double_t pz = tmp_py*0.001;

  TVector3 p = TVector3(px,py,pz);
  if(m_charge < 0) p *= -1.;
  return p;
}

//______________________________________________________________________________
TPCLTrackHit*
TPCLocalTrackHelix::GetHit(std::size_t nth) const
{
  if(nth<m_hit_array.size() && nth>=0)
    return m_hit_array[nth];
  else
    return 0;
}

//______________________________________________________________________________
Int_t
TPCLocalTrackHelix::GetOrder(Int_t i) const
{
  Int_t size = m_hit_order.size();
  if(i>=size || i<0) return -1;
  Int_t id = i;
  if(m_charge<0) id = size - i - 1;
  Int_t order = m_hit_order[id];
  return order;
}

//______________________________________________________________________________
TPCLTrackHit*
TPCLocalTrackHelix::GetHitInOrder(std::size_t nth) const
{
  Int_t size = m_hit_order.size();
  Int_t id = nth;
  if(m_charge<0) id = size - nth - 1;
  Int_t order = m_hit_order[id];
  if(nth>=m_hit_array.size() || nth<0) return 0;
  return m_hit_array[order];
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::DeleteNullHit()
{

  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hit = m_hit_array[i];
    if(!hit->IsGood()){
      m_is_theta_calculated = false;
      hddaq::cout << FUNC_NAME << " "
		  << "null hit has been deleted" << std::endl;
      m_hit_array.erase(m_hit_array.begin()+i);
      m_hit_order.erase(m_hit_order.begin()+i);
      m_hit_t.erase(m_hit_t.begin()+i);
      --i;
    }
  }
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::SetClustersHoughFlag(Int_t hough_flag)
{
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TPCHit *hit = hitp->GetHit();
    if( !hit ) continue;
    hit->SetHoughFlag(hough_flag);
  }
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::SetFlag(Int_t flag)
{

  if(flag==0){ m_isBeam=0; m_isK18=0; m_isKurama=0; m_isAccidental=0; } //reset
  if((flag&1)==1) m_isBeam=1;
  if((flag&2)==2) m_isK18=1;
  if((flag&4)==4) m_isKurama=1;
  if((flag&8)==8) m_isAccidental=1;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoFit(Int_t MinHits)
{

  gMomConstraint = false; //No momentum constraint
  if(!IsGoodForTracking()) return false;

  Bool_t status = DoHelixTrackFit(); //track chisqr minimization
  m_is_fitted = status;

  if(!status || m_chisqr > MaxChisqr) return false;

  SeparateTracksAtTarget();

  //If track speration is performed, m_is_fitted flag is set to false
  if(!m_is_fitted) return DoFit(MinHits); //Do chisqr minimization again after separation

  //Minimum # of clusters
#if 1
  if(IsBackward() && GetIsBeam()==1) MinHits = 3;
  Int_t nhit = GetNHit();
  if(nhit<MinHits) return false;
#else
  Int_t nbadhit = gBadHits;
  if(nhit-nbadhit<MinHits) return false;
#endif

  return status;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoHelixTrackFit(Double_t RKpar[5]) //for beam track fitting
{

  if(!gMomConstraint){
    std::cout<<FUNC_NAME+" Fatal error : No momentum constraint"<<std::endl;
    return false;
  }

  Bool_t vetoBadClusters = false;

  DeleteNullHit();
  const std::size_t n = m_hit_array.size();
  if(!IsGoodForTracking()) return false;
  gMultiLoop = m_is_multiloop;
  gNumOfHits = n;
  gHitPos.clear();
  gLayer.clear();
  gPadTheta.clear();
  gResParam.clear();
  gHelixTheta.clear();
  gBadHits = 0;
  gMinuitStatus = 0;
  gChisqr = 1.e+10;
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
  }

#if DebugDisp
  std::cout<<FUNC_NAME+" Inital helix params before fitting"<<std::endl
	   <<" cx: "<<m_cx
	   <<", cy: "<<m_cy
	   <<", z0: "<<m_z0
	   <<", r: "<<m_r
	   <<", dz: "<<m_dz<<std::endl;
#endif

  //1. calculate chisqr by using RK parameters
  SetParam(RKpar);
  CalcHelixTheta();
  Int_t ndf;
  Double_t RKChisqr = CalcChi2(RKpar, ndf, vetoBadClusters);
  if(RKChisqr>10.){
    //2. calculate chisqr with pre-fitting parameters (gChisqr)
    if(!DoPreFit(RKpar)) return false;

#if DebugDisp
    std::cout<<FUNC_NAME+" Helix params after pre-fitting (w/ constraints from RK)"<<std::endl
	     <<", chisqr: "<<gChisqr<<std::endl
	     <<" cx: "<<gPar[0]
	     <<", cy: "<<gPar[1]
	     <<", z0: "<<gPar[2]
	     <<", r: "<<gPar[3]
	     <<", dz: "<<gPar[4]<<std::endl;
#endif
  }

  //Compare 1 & 2 and choose the better one
  if(gChisqr > RKChisqr){
    gChisqr = RKChisqr;
    gPar[0] = RKpar[0];
    gPar[1] = RKpar[1];
    gPar[2] = RKpar[2];
    gPar[3] = RKpar[3];
    gPar[4] = RKpar[4];
  }
  SetParam(gPar);
  CalcHelixTheta();

  //Helix fitting
  if(!DoHelixFit(gPar, vetoBadClusters)) return false;

#if DebugDisp
  std::cout<<FUNC_NAME+" After helix fitting"<<std::endl
	   <<" n_iteration: "<<m_n_iteration
	   <<", chisqr: "<<gChisqr<<std::endl
	   <<", cx: "<<gPar[0]
	   <<", cy: "<<gPar[1]
	   <<", z0: "<<gPar[2]
	   <<", r: "<<gPar[3]
	   <<", dz: "<<gPar[4]<<std::endl;
#endif

  Int_t delete_hit = -1;
  Int_t false_layer = FinalizeTrack(delete_hit);
  Bool_t goodtrack = true;
  if(false_layer > 0 || m_chisqr > MaxChisqr){
    goodtrack = false;

#if DebugDisp
    TPCLTrackHit *hitp = m_hit_array[delete_hit];
    TVector3 pos = hitp->GetLocalHitPos();
    std::cout<<"delete hits ["<<delete_hit<<"]=("
    	     <<pos.x()<<", "
      	     <<pos.y()<<", "
      	     <<pos.z()<<")"<<std::endl;
#endif

    EraseHit(delete_hit);
    gHelixTheta.erase(gHelixTheta.begin()+delete_hit);
  }

  m_n_iteration++;
  if(m_n_iteration > MaxIteration) return false;
  if(goodtrack){ //Tracking is over.
#if IterativeResolution
    vetoBadClusters = true;

    //Now excluding bad clusters and fitting again.
    DoHelixFit(gPar, vetoBadClusters);
#if DebugDisp
    std::cout<<FUNC_NAME+" After excluding bad hits"<<std::endl
	     <<" n_iteration: "<<m_n_iteration
	     <<", chisqr: "<<gChisqr<<std::endl
	     <<", cx: "<<gPar[0]
	     <<", cy: "<<gPar[1]
	     <<", z0: "<<gPar[2]
	     <<", r: "<<gPar[3]
	     <<", dz: "<<gPar[4]<<std::endl;
#endif
#endif
    return true;
  }
  else return DoHelixTrackFit(RKpar);
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoFit(Double_t RKpar[5], Int_t MinHits)
{

  gMomConstraint = true; //fit with the momentum constraint
  if(!IsGoodForTracking()) return false;

  Bool_t status = DoHelixTrackFit(RKpar); //track chisqr minimization
  m_is_fitted = status;
  if(!status || m_chisqr > MaxChisqr) return false;

  //SeparateTracksAtTarget();

  //Minimum # of clusters
#if 1
  if(IsBackward() && GetIsBeam()==1) MinHits = 3;
  Int_t nhit = GetNHit();
  if(nhit<MinHits) return false;
#else
  Int_t nbadhit = gBadHits;
  if(nhit-nbadhit<MinHits) return false;
#endif
  return status;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoFit(Double_t RKCharge, Double_t RKpar[5], Int_t MinHits)
{

  gMomConstraint = true; //fit with a momentum constraint
  Bool_t status = DoFit(RKpar, MinHits);
  Int_t charge = (int) RKCharge;
  if(m_charge != charge) return false; //fit with a charge constraint
  return status;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoFit(Double_t RKCharge, Int_t MinHits)
{

  Bool_t status = DoFit(MinHits);
  Int_t charge = (int) RKCharge;
  if(m_charge != charge) return false; //fit with a charge constraint
  return status;
}

//_____________________________________________________________________________
void
TPCLocalTrackHelix::EraseHits(std::vector<Int_t> delete_hits)
{

  std::sort(delete_hits.begin(), delete_hits.end());
  for(Int_t i=0; i<delete_hits.size(); ++i){
    //Reset houghflag
    EraseHit(delete_hits[i]-i);
  }
  m_is_theta_calculated = false;

}

//_____________________________________________________________________________
void
TPCLocalTrackHelix::EraseHit(Int_t delete_hit)
{

  //Reset houghflag
  TPCLTrackHit *hitp = m_hit_array[delete_hit];
  TPCHit *hit = hitp->GetHit();
  Int_t order = m_hit_order[delete_hit];
  hit->SetHoughFlag(0);
  m_hit_array.erase(m_hit_array.begin()+delete_hit);
  m_hit_t.erase(m_hit_t.begin()+delete_hit);
  m_hit_order.erase(m_hit_order.begin()+delete_hit);
  for(Int_t i=0; i<m_hit_order.size(); ++i){
    if(m_hit_order[i] > order) m_hit_order[i] -= 1;
  }
  m_is_theta_calculated = false;

}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::IsGoodForTracking()
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
Bool_t
TPCLocalTrackHelix::DoCircleFit(Double_t *par)
{

  DeleteNullHit();
  const std::size_t n = m_hit_array.size();
  if(!IsGoodForTracking()) return false;
  gNumOfHits = n;

  Double_t xp[n]; Double_t yp[n];
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 localpos = GlobalToLocal(pos);
    xp[i] = localpos.X();
    yp[i] = localpos.Y();
  }

  Double_t par_circ[3]={0};
  Int_t npad = GetNPad();
  if(npad==3) return true;
  if(npad<3 || CircleFit(xp, yp, n, &par_circ[0], &par_circ[1], &par_circ[2])<0) return false;

  par[0] = par_circ[0];
  par[1] = par_circ[1];
  par[3] = par_circ[2];

#if DebugDisp
  std::cout<<FUNC_NAME+" Circlefit results"<<std::endl
	   <<" cx: "<<par[0]
	   <<", cy: "<<par[1]
	   <<", r: "<<par[3]<<std::endl;
#endif

  return true;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoCircleFitwMomConstraint(Double_t *par)
{

  if(!gMomConstraint){
    std::cout<<FUNC_NAME+" Fatal error : No momentum constraint"<<std::endl;
    return false;
  }

#if DebugDisp
  std::cout<<FUNC_NAME+" Helix params before circle fitting w/ the momentum constraint"<<std::endl
	   <<" cx: "<<par[0]
	   <<", cy: "<<par[1]
	   <<", r: "<<par[3]<<std::endl;
#endif

  DeleteNullHit();
  const std::size_t n = m_hit_array.size();
  if(!IsGoodForTracking()) return false;
  gNumOfHits = n;
  gHitPos.clear();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    gHitPos.push_back(pos);
  }

  gPar[0] = par[0];
  gPar[1] = par[1];
  gPar[3] = par[3];
  CircleFit(m_isBeam);

  par[0] = gPar[0];
  par[1] = gPar[1];
  par[3] = gPar[3];

#if DebugDisp
  std::cout<<FUNC_NAME+" Circlefit results"<<std::endl
	   <<" cx: "<<par[0]
	   <<", cy: "<<par[1]
	   <<", r: "<<par[3]<<std::endl;
#endif

  return true;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoStraightLineFit(Double_t *par)
{
  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return false;
  }

  DeleteNullHit();
  const std::size_t n = m_hit_array.size();
  if(!IsGoodForTracking()) return false;
  gMultiLoop = m_is_multiloop;
  gNumOfHits = n;
  gHitPos.clear();
  gLayer.clear();
  gPadTheta.clear();
  gResParam.clear();
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
  }

  Bool_t vetoBadClusters = false;
  gRes.clear();
  for(Int_t i=0; i<gNumOfHits; i++){
    Double_t dummy_theta = 0;
    TVector3 res = CalcResolution(gPar, gLayer[i], gHitPos[i], gPadTheta[i], dummy_theta, gResParam[i], vetoBadClusters);
    gRes.push_back(res);
  }

#if DebugDisp
  std::cout<<FUNC_NAME+" Helix params before straight-line fitting"<<std::endl
	   <<" z0: "<<gPar[2]
	   <<", dz: "<<gPar[4]<<std::endl;
#endif

  gPar[0] = par[0];
  gPar[1] = par[1];
  gPar[2] = par[2];
  gPar[3] = par[3];
  gPar[4] = par[4];

#if 1 //Spiral-like track with low p_T
  if(gMultiLoop) return true;
#endif
  return StraightLineFit();
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoPreFit(Double_t par[5])
{

  Bool_t vetoBadClusters = false;

  gPar[0] = par[0];
  gPar[1] = par[1];
  gPar[2] = par[2];
  gPar[3] = par[3];
  gPar[4] = par[4];

  //pre circle fit
  Bool_t pass = false;
  pass = DoCircleFit(gPar);
  if(!pass && gMomConstraint) pass = DoCircleFitwMomConstraint(gPar);
  if(!pass) return false;

#if 1
  const Int_t n = m_hit_array.size();
  gNumOfHits = n;
  gHitPos.clear();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    gHitPos.push_back(pos);
  }

  Int_t MaxBin[3];
  if(pass && (gPar[0] < LowLimit[0] ||
	      gPar[0] > UpLimit[0] ||
	      gPar[1] < LowLimit[1] ||
	      gPar[1] > UpLimit[1] ||
	      gPar[3] < LowLimit[3] ||
	      gPar[3] > UpLimit[3]))
    pass = tpc::HoughTransformCircleXZ(gHitPos, MaxBin, gPar, 3);
  if(!pass) return false; //For very high momentum tracks
#endif
  SetParam(gPar);
  CalcHelixTheta();

#if 1 //Optional
  Int_t MaxBinY[3];
  if(pass && (gPar[2] < LowLimit[2] ||
	      gPar[2] > UpLimit[2] ||
	      gPar[4] < LowLimit[4] ||
	      gPar[4] > UpLimit[4])) tpc::HoughTransformLineYTheta(gHitPos, MaxBinY, gPar, 1000.);
#endif
  if(!DoStraightLineFit(gPar)) return false;
  SetParam(gPar);

  Int_t ndf;
  gChisqr = CalcChi2(gPar, ndf, vetoBadClusters);
  return true;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoHelixFit(Double_t *par, Bool_t vetoBadClusters)
{
  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return false;
  }

  gPar[0] = par[0];
  gPar[1] = par[1];
  gPar[2] = par[2];
  gPar[3] = par[3];
  gPar[4] = par[4];

  IsBackward();

  Bool_t status = false;
  if(HelixFit(m_isBeam, vetoBadClusters)){
    status = true;
    SetParam(gPar);
    CalcHelixTheta();
    m_chisqr = gChisqr;
    m_minuit = gMinuitStatus;
    Double_t window = ThetaWindow/gPar[3];
    for(Int_t i=0; i<gNumOfHits; ++i){
      Double_t theta = EvalTheta(gPar, gHitPos[i], gHelixTheta[i] - 0.5*window, gHelixTheta[i] + 0.5*window);
      m_hit_t[i] = theta;
      gHelixTheta[i] = theta;
    }
  }
  return status;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoHelixTrackFit()
{

  if(gMomConstraint){
    std::cout<<FUNC_NAME+" Fatal error : Momentum constraint is applied"<<std::endl;
    return false;
  }

#if DebugDisp
  std::cout<<FUNC_NAME+" Helix params (Hough-Transform)"<<std::endl
	   <<" cx: "<<m_cx
	   <<", cy: "<<m_cy
	   <<", z0: "<<m_z0
	   <<", r: "<<m_r
	   <<", dz: "<<m_dz<<std::endl;
#endif

  Bool_t vetoBadClusters = false;

  DeleteNullHit();
  const Int_t n = m_hit_array.size();
  if(!IsGoodForTracking()) return false;
  gMultiLoop = m_is_multiloop;
  gNumOfHits = n;
  gHitPos.clear();
  gLayer.clear();
  gPadTheta.clear();
  gResParam.clear();
  gBadHits = 0;
  gMinuitStatus = 0;

  //Initialization params with Hough-transform result or previous tracking result
  gPar[0] = m_cx;
  gPar[1] = m_cy;
  gPar[2] = m_z0;
  gPar[3] = m_r;
  gPar[4] = m_dz;

  //pre fitting
  if(!DoPreFit(gPar)) return false;
#if DebugDisp
  std::cout<<FUNC_NAME+" After pre fitting"<<std::endl
	   <<" n_iteration: "<<m_n_iteration
	   <<", chisqr: "<<gChisqr<<std::endl
	   <<", cx: "<<gPar[0]
	   <<", cy: "<<gPar[1]
	   <<", z0: "<<gPar[2]
	   <<", r: "<<gPar[3]
	   <<", dz: "<<gPar[4]<<std::endl;
#endif

  //Helix fitting
  if(!DoHelixFit(gPar, vetoBadClusters)) return false;
#if DebugDisp
  std::cout<<FUNC_NAME+" After helix fitting"<<std::endl
	   <<" n_iteration: "<<m_n_iteration
	   <<", chisqr: "<<gChisqr<<std::endl
	   <<", cx: "<<gPar[0]
	   <<", cy: "<<gPar[1]
	   <<", z0: "<<gPar[2]
	   <<", r: "<<gPar[3]
	   <<", dz: "<<gPar[4]<<std::endl;
#endif

  Int_t delete_hit = -1;
  Int_t false_layer = FinalizeTrack(delete_hit);
  Bool_t goodtrack = true;
  if(false_layer > 0 || m_chisqr > MaxChisqr){
    goodtrack = false;
#if DebugDisp
    TPCLTrackHit *hitp = m_hit_array[delete_hit];
    TVector3 pos = hitp->GetLocalHitPos();
    std::cout<<"delete hit #"<<delete_hit<<"=("
	     <<pos.x()<<", "
	     <<pos.y()<<", "
	     <<pos.z()<<")"<<std::endl;
#endif

    EraseHit(delete_hit);
    gHelixTheta.erase(gHelixTheta.begin()+delete_hit);
  }

  m_n_iteration++;
  if(m_n_iteration > MaxIteration) return false;
  if(goodtrack){ //Tracking is over.
#if IterativeResolution
    vetoBadClusters = true;

    //Now excluding bad clusters and fitting again.
    DoHelixFit(gPar, vetoBadClusters);
#if DebugDisp
    std::cout<<FUNC_NAME+" with precise resolution calculation"<<std::endl
	     <<" n_iteration: "<<m_n_iteration
	     <<", chisqr: "<<gChisqr<<std::endl
	     <<", cx: "<<gPar[0]
	     <<", cy: "<<gPar[1]
	     <<", z0: "<<gPar[2]
	     <<", r: "<<gPar[3]
	     <<", dz: "<<gPar[4]<<std::endl;
#endif

#endif
    return true;
  }
  else return DoHelixTrackFit();
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::ResidualCheck(Int_t i, Double_t &residual)
{

  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return false;
  }

  TPCHit *hit = m_hit_array[i] -> GetHit();
  Int_t layer = hit->GetLayer();
  Double_t padTheta = hit->GetPadTheta();
  std::vector<Double_t> resparam = hit->GetResolutionParams();

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 position = m_hit_array[i] -> GetLocalHitPos();
  TVector3 res = CalcResolution(par, layer, position, padTheta, gHelixTheta[i], resparam, false);
  TVector3 resi = ResidualVect(par, position, gHelixTheta[i]); //Closest distance
  residual = resi.Mag();

  //XZ residual < window
  TVector3 residualXZ = ResidualVectXZ(par, position);
  if(m_is_multiloop){
    if(residualXZ.Mag() > ResidualWindowOutXZ) return false;
  }
  else if(layer < 3 && TMath::Abs(position.y()) < 30.){
    if(residualXZ.Mag() > ResidualWindowUnderTgtXZ) return false;
  }
  else if(layer < 10 && residualXZ.Mag() > ResidualWindowInXZ) return false;
  else if(layer >=10 && residualXZ.Mag() > ResidualWindowOutXZ) return false;

  //Residual/resolution < window
  Double_t residual_vertical = resi.y();
  Double_t residual_horizontal = TMath::Hypot(resi.x(), resi.z());
  Double_t resolution_vertical = res.y();
  Double_t resolution_horizontal = TMath::Hypot(res.x(), res.z());
  if(TMath::Abs(residual_vertical) > resolution_vertical*ResidualWindowPullY) return false;
  if(TMath::Abs(residual_horizontal) > resolution_horizontal*ResidualWindowPullXZ) return false;
  return true;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::IsGoodHitToAdd(TPCHit *hit, Double_t &residual, Bool_t nolimitation)
{

  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return false;
  }

  Int_t layer = hit->GetLayer();
  Double_t padTheta = hit->GetPadTheta();
  std::vector<Double_t> resparam = hit->GetResolutionParams();

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 position = hit->GetPosition();
  Double_t dummy_theta = 0;
  TVector3 res = CalcResolution(par, layer, position, padTheta, dummy_theta, resparam, false);

  Int_t upstream_tgt = -1;
  if(TMath::Abs(position.x()) < 25. &&
     TMath::Abs(position.y()) < 10. &&
     position.z() < tpc::ZTarget) upstream_tgt = 0; //Beam section
  else if(TMath::Abs(position.x()) < 25. &&
	  TMath::Abs(position.y()) > 10. &&
	  position.z() < tpc::ZTarget) upstream_tgt = 1; //above or below the beam section

  Int_t section = -1;
  if((position.z() + position.x()) < 0 && (position.z() - position.x()) < 0) section = 1;
  if((position.z() + position.x()) < 0 && (position.z() - position.x()) > 0) section = 2;
  if((position.z() + position.x()) > 0 && (position.z() - position.x()) > 0) section = 3;
  if((position.z() + position.x()) > 0 && (position.z() - position.x()) < 0) section = 4;

  Double_t factor = 1.;
  Int_t nhit_upstream_tgt = 0;
  std::vector<Int_t> gSection;
  const std::size_t n = m_hit_array.size();
  for(Int_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    Int_t section = hitp->GetSection();
    gSection.push_back(section);

    TVector3 pos = hitp -> GetLocalHitPos();
    if(TMath::Abs(pos.x()) < 25. && pos.z() < tpc::ZTarget) nhit_upstream_tgt++;
  }

  Double_t max_scanrange = 40.;
  if(find(gSection.begin(), gSection.end(), section) == gSection.end()){ //Crossing to the other sections
    if(TMath::Abs(position.y()) > 50) factor = 2.5;
    else factor = 2.;

    Bool_t include_section1 = (find(gSection.begin(), gSection.end(), 1) == gSection.end());
    Bool_t include_section2 = (find(gSection.begin(), gSection.end(), 2) == gSection.end());
    Bool_t include_section3 = (find(gSection.begin(), gSection.end(), 3) == gSection.end());
    Bool_t include_section4 = (find(gSection.begin(), gSection.end(), 4) == gSection.end());
    if(section==1){
      if(include_section2 || include_section4) max_scanrange += 60.;
      else max_scanrange += 120.;
    }
    else if(section==2){
      if(include_section1 || include_section3) max_scanrange += 60.;
      else max_scanrange += 120.;
    }
    else if(section==3){
      if(include_section2 || include_section4) max_scanrange += 60.;
      else max_scanrange += 120.;
    }
    else if(section==4){
      if(include_section1 || include_section3) max_scanrange += 60.;
      else max_scanrange += 120.;
    }
  }
  else if(nhit_upstream_tgt==0 && upstream_tgt==1 &&
	  (TMath::Hypot(m_closedist.x(), m_closedist.z())<25. &&
	   TMath::Abs(m_closedist.y())>10.)){ //Under or over the target
    factor = 1.;
    max_scanrange += 40.;
  }
  else if(nhit_upstream_tgt==0 && upstream_tgt==0){ //Beam section
    factor = 1.;
    max_scanrange = 20.;
  }
  if(nolimitation){
    max_scanrange = 300;
    factor = 2;
  }
  max_scanrange /= m_r;

  Double_t theta_range[2];
  Double_t max_scanrange_theta = m_is_multiloop ? TMath::Pi() : 0.5*TMath::Pi();
  theta_range[0] = m_min_t - TMath::Min(max_scanrange_theta, max_scanrange);
  theta_range[1] = m_max_t + TMath::Min(max_scanrange_theta, max_scanrange);

  //XZ residual < window
  TVector3 residualXZ = ResidualVectXZ(par, position);
#if 1
  if(m_r < 250.){
    if(residualXZ.Mag() > 3.*ResidualWindowOutXZ) return false;
  }
#else
  if(m_is_multiloop){
    if(residualXZ.Mag() > factor*ResidualWindowOutXZ) return false;
  }
#endif
  else if(layer < 3 && TMath::Abs(position.y()) < 30.){
    if(residualXZ.Mag() > factor*ResidualWindowUnderTgtXZ) return false;
  }
  else if(layer < 10 && residualXZ.Mag() > factor*ResidualWindowInXZ) return false;
  else if(layer >= 10 && residualXZ.Mag() > factor*ResidualWindowOutXZ) return false;

  TVector3 resi = ResidualVect(par, position, theta_range[0], theta_range[1]);
  residual = resi.Mag();

  //Residual/resolution < window
  Double_t residual_vertical = resi.y();
  Double_t residual_horizontal = TMath::Hypot(resi.x(), resi.z());
  Double_t resolution_vertical = res.y();
  Double_t resolution_horizontal = TMath::Hypot(res.x(), res.z());

#if 1
  if(nolimitation){
    if(TMath::Abs(residual_horizontal) > 5.*resolution_horizontal*ResidualWindowPullXZ) return false;
    if(TMath::Abs(residual_vertical) > 3.*resolution_vertical*ResidualWindowPullY) return false;
  }
  else if(m_r < 250.){ //wider cut condition for low-momentum tracks (currently not supported)
    if(TMath::Abs(residual_horizontal) > factor*resolution_horizontal*ResidualWindowPullXZ) return false;
    if(TMath::Abs(residual_vertical) > factor*resolution_vertical*ResidualWindowPullY) return false;
  }
  else{
    if(TMath::Abs(residual_horizontal) > factor*resolution_horizontal*ResidualWindowPullXZ) return false;
    if(TMath::Abs(residual_vertical) > factor*resolution_vertical*ResidualWindowPullY) return false;
  }
#else
  if(TMath::Abs(residual_horizontal) > factor*resolution_horizontal*ResidualWindowPullXZ) return false;
  if(TMath::Abs(residual_vertical) > factor*resolution_vertical*ResidualWindowPullY) return false;
#endif
  return true;
}

//_____________________________________________________________________________
Int_t
TPCLocalTrackHelix::Side(TVector3 hitpos)
{

  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return 0;
  }

  if(m_hit_t.size()==0){
    std::cout<<FUNC_NAME+" Fatal error : Empty track!"<<std::endl;
    return 0;
  }

  //TPC local coordinate
  Int_t flag = -1;
  TVector3 Vect1(-hitpos.X() - m_cx, hitpos.Z() - tpc::ZTarget - m_cy, 0.); //Vect1(Hit - Helix center)
  TVector3 Vect2(-m_cx, -m_cy, 0.); //Vec2(Tgt - Helix center)

  TVector3 norm = Vect1.Cross(Vect2); //Vect1 X Vec2
  if(norm.Z()>0.) flag = 1;
  return flag;
}

//______________________________________________________________________________
Double_t
TPCLocalTrackHelix::CalcThetaExclusive(Int_t ith) const
{

  if(!m_is_fitted_exclusive){
    std::cout<<FUNC_NAME+" Fatal error : Do exclusive fitting first!"<<std::endl;
    return TMath::QuietNaN();
  }

  TPCLTrackHit *hitp = m_hit_array[ith];
  TVector3 pos = hitp -> GetLocalHitPos();
  TVector3 localpos = GlobalToLocal(pos);
  Double_t par[5] = {m_cx_exclusive[ith], m_cy_exclusive[ith], m_z0_exclusive[ith], m_r_exclusive[ith], m_dz_exclusive[ith]};
  Double_t window = ThetaWindow/par[3];
  Double_t theta = EvalTheta(par, pos, m_hit_t[ith] - 0.5*window, m_hit_t[ith] + 0.5*window);

  return theta;
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::CalcHelixTheta()
{

  if(m_hit_array.size()==0){
    std::cout<<FUNC_NAME+" Fatal error : Empty track!"<<std::endl;
    return;
  }

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
#if DebugDisp
  std::cout<<FUNC_NAME<<std::endl
	   <<", cx: "<<m_cx
	   <<", cy: "<<m_cy
	   <<", z0: "<<m_z0
	   <<", r: "<<m_r
	   <<", dz: "<<m_dz<<std::endl;
#endif

  gHelixTheta.clear();
  std::vector<Double_t> helix_theta;

  //compare a variance of y positions to deterimine theta with Atan2 or Atan2 + pi, Atan2 - pi
  Double_t variance_group0 = 0.; Double_t variance_group1 = 0.;

  Bool_t thetaflip = false;
  Double_t prev_theta = 0.; Double_t theta0 = 0;
  m_min_t = 9999; m_max_t = -9999;
  const std::size_t n = m_hit_array.size();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp -> GetLocalHitPos();

    //tmp_theta : Atan2 output with range (-pi, pi)
    TVector3 localpos = GlobalToLocal(pos);
    Double_t tmp_theta = TMath::ATan2(localpos.y() - m_cy, localpos.x() - m_cx);
    if(m_is_multiloop){
      //Calculate Atan2(y,x) and turns of helix. And convert them into the theta of helix.
      Double_t pitch = 2.*TMath::Pi()*m_r*m_dz;
      Double_t ref_ypos = GetPosition(par, tmp_theta).y();
      Double_t turns = TMath::Nint((pos.y() - ref_ypos)/pitch);
      tmp_theta += 2.*TMath::Pi()*turns;
    }
    else{
      //Check ATan2 function's theta flip (-pi ~ pi) within a loop.
      if(i==0) theta0 = tmp_theta;
      else if(TMath::Abs(prev_theta - tmp_theta) > TMath::Pi()) thetaflip = true;
      if(thetaflip){
	if(theta0>0 && tmp_theta<0) tmp_theta += 2.*TMath::Pi();
	if(theta0<0 && tmp_theta>0) tmp_theta -= 2.*TMath::Pi();
      }
      variance_group0 += TMath::Power(GetPosition(par, tmp_theta).y() - pos.y(), 2.);
      if(theta0>0) variance_group1 += TMath::Power(GetPosition(par, tmp_theta - 2.*TMath::Pi()).y() - pos.y(), 2.);
      else variance_group1 += TMath::Power(GetPosition(par, tmp_theta + 2.*TMath::Pi()).y() - pos.y(), 2.);
    }
    prev_theta = tmp_theta;
    helix_theta.push_back(tmp_theta);
  }

  for(std::size_t i=0; i<n; ++i){
    Double_t tmp_theta = helix_theta[i];
    if(!m_is_multiloop && variance_group0 > variance_group1
       && TMath::Abs(variance_group0 - variance_group1) > 100){ //almost flat track(pT~0) is excluded.
      if(helix_theta[0] > 0) tmp_theta -= 2.*TMath::Pi();
      else tmp_theta += 2.*TMath::Pi();
    }
    if(tmp_theta < m_min_t) m_min_t = tmp_theta;
    if(tmp_theta > m_max_t) m_max_t = tmp_theta;
    gHelixTheta.push_back(tmp_theta);
  }

#if DebugDisp
  std::cout<<"helix theta calculated: theta min "<<m_min_t<<" max "<<m_max_t<<std::endl;
#endif
  m_is_theta_calculated = true;

}

//______________________________________________________________________________
void
TPCLocalTrackHelix::DoFitExclusive()
{
  const std::size_t n = m_hit_array.size();

  if(m_hit_t.size() <= 1 || !m_is_fitted || !m_is_calculated) std::cout<<FUNC_NAME+" Fatal error : Wrong track. Please check"<<std::endl;

  gMultiLoop = m_is_multiloop;

  m_is_fitted_exclusive = true;
  m_cx_exclusive.resize(n);
  m_cy_exclusive.resize(n);
  m_z0_exclusive.resize(n);
  m_r_exclusive.resize(n);
  m_dz_exclusive.resize(n);
  m_t_exclusive.resize(n);
  m_chisqr_exclusive.resize(n);

  gNumOfHits = n-1;
  for(Int_t ihit=0; ihit<n; ++ihit){ //exclude the ith hit
    gHitPos.clear();
    gLayer.clear();
    gPadTheta.clear();
    gResParam.clear();
    gHelixTheta.clear();
    gHitPos.resize(n-1);
    gLayer.resize(n-1);
    gPadTheta.resize(n-1);
    gResParam.resize(n-1);
    gHelixTheta.resize(n-1);

    gChisqr = 1.e+10;
    gPar[0] = m_cx;
    gPar[1] = m_cy;
    gPar[2] = m_z0;
    gPar[3] = m_r;
    gPar[4] = m_dz;

    Int_t flag=0;
    for(Int_t j=0; j<n; ++j){
      if(j==ihit) continue; //exclude the ith hit
      TPCLTrackHit *hitp = m_hit_array[j];
      TVector3 pos = hitp->GetLocalHitPos();
      Int_t layer = hitp->GetLayer();
      Double_t padTheta = hitp->GetPadTheta();
      std::vector<Double_t> resparam = hitp->GetResolutionParams();
      gHitPos[flag] = pos;
      gLayer[flag] = layer;
      gPadTheta[flag] = padTheta;
      gHelixTheta[flag] = hitp->GetTheta();
      gResParam[flag] = resparam;
      flag++;
    } //j

    //Helix fitting
    Bool_t vetoBadClusters = true;
    Bool_t exclusive = true;
    HelixFit(m_isBeam, vetoBadClusters, exclusive);
    m_chisqr_exclusive[ihit] = gChisqr;
    m_cx_exclusive[ihit] = gPar[0];
    m_cy_exclusive[ihit] = gPar[1];
    m_z0_exclusive[ihit] = gPar[2];
    m_r_exclusive[ihit]  = gPar[3];
    m_dz_exclusive[ihit] = gPar[4];
    m_t_exclusive[ihit] = CalcThetaExclusive(ihit);

#if DebugDisp
    std::cout<<"exclusive "<<ihit<<" th chisqr : "<<m_chisqr_exclusive[ihit]<<
      " par : "<<m_cx_exclusive[ihit]<<
      " "<<m_cy_exclusive[ihit]<<
      " "<<m_z0_exclusive[ihit]<<
      " "<<m_r_exclusive[ihit]<<
      " "<<m_dz_exclusive[ihit]<<
      " "<<m_t_exclusive[ihit]<<std::endl;
#endif
  } //ihit
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::Print(const TString& arg, Bool_t print_allhits) const
{
  TString tracksize = Form(" #clusters = %d", (Int_t)m_hit_array.size());
  TString trackpar = Form(", params cx:%f cy:%f z0:%f r:%f dz:%f", m_cx, m_cy, m_z0, m_r, m_dz);
  std::cout<<arg.Data()<<std::endl;
  std::cout<<"Track info : "<<tracksize.Data()<<trackpar.Data()<<std::endl;
  if(print_allhits){
    for(std::size_t i=0; i<m_hit_array.size(); ++i){
      TPCLTrackHit *hitp = m_hit_array[i];
      if( !hitp ) continue;
      hitp->Print();
      std::cout<<"theta : "<<m_hit_t[i]<<std::endl;
    }
  }
}

//______________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetdEdx(Double_t truncatedMean)
{

  std::vector<Double_t> dEdx_vect;
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    Double_t clde = hitp->GetDe();
    Double_t pathHit = hitp->GetPathHelix();
    Double_t dEdx_cor = clde/pathHit;
    dEdx_vect.push_back(dEdx_cor);
  }

  Double_t dEdx = 0.;
  std::sort(dEdx_vect.begin(), dEdx_vect.end());
  Int_t n_truncated = (Int_t)(dEdx_vect.size()*truncatedMean);
  for( Int_t ih=0; ih<dEdx_vect.size(); ++ih ){
    if(ih<n_truncated) dEdx += dEdx_vect[ih];
  }
  dEdx /= (Double_t)n_truncated;

  return dEdx;
}

//______________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetTrackdE()
{

  Double_t dE = 0;
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    dE += hitp->GetDe();
  }
  return dE;
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::AddVPHit(TVector3 vp)
{
  m_hit_order.push_back(m_vp.size());
  m_vp.push_back(vp);
  m_hit_t.push_back(TMath::QuietNaN());

}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoVPFit()
{

  const std::size_t n = m_vp.size();
  gMultiLoop = false;
  gNumOfHits = n;
  gHitPos.clear();
  gHitPos.resize(n);
  gRes.clear();
  gRes.resize(n);
  gLayer.clear();
  gLayer.resize(n);
  gPadTheta.clear();
  gPadTheta.resize(n);
  gResParam.clear();
  gResParam.resize(n);
  gHelixTheta.clear();
  gHelixTheta.resize(n);
  gBadHits = 0;
  gMinuitStatus = 0;

  //for pre circle fitting
  Double_t xp[n]; Double_t yp[n];
  for(std::size_t i=0; i<n; ++i){
    gHitPos[i] = m_vp[i];
    xp[i] = -m_vp[i].X();
    yp[i] = m_vp[i].Z() - tpc::ZTarget;
    gRes[i] = TVector3(1., 1., 1.); //dummy values
  }

  Double_t par_circ[3]={0};
  //1st circle fitting to get initial parameters
  if(CircleFit(xp, yp, n, &par_circ[0], &par_circ[1], &par_circ[2])<0) return false;
  gPar[0] = par_circ[0];
  gPar[1] = par_circ[1];
  gPar[2] = 0.;
  gPar[3] = par_circ[2];
  gPar[4] = 0.;
  SetParam(gPar);

  Bool_t thetaflip = false;
  Double_t prev_theta = 0.; Double_t theta0 = 0;
  m_min_t = 9999; m_max_t = -9999;
  for(std::size_t i=0; i<n; ++i){
    Double_t theta = TMath::ATan2(m_vp[i].Z() - tpc::ZTarget - gPar[1], -m_vp[i].X() - gPar[0]);
    //Check ATan2 function's theta flip (-pi ~ pi)
    if(i==0) theta0 = theta;
    if(TMath::Abs(prev_theta - theta) > TMath::Pi()) thetaflip = true;
    if(thetaflip){
      if(theta0>0 && theta<0) theta += 2.*TMath::Pi();
      if(theta0<0 && theta>0) theta -= 2.*TMath::Pi();
    }
    prev_theta = theta;
    if(theta < m_min_t) m_min_t = theta;
    if(theta > m_max_t) m_max_t = theta;
    gHelixTheta[i] = theta;
    m_hit_t[i] = theta;
  }
  m_is_theta_calculated = true;

#if 0 //Optional
  Int_t MaxBinY[3];
  tpc::HoughTransformLineYTheta(gHitPos, MaxBinY, gPar, 1000.);
#endif
  if(!StraightLineFit()) return false;
  SetParam(gPar);
  m_minuit = gMinuitStatus;

#if DebugDisp
  std::cout<<FUNC_NAME+" RK helix params : "<<m_cx<<" "<<m_cy<<" "<<m_z0<<" "<<m_r<<" "<<m_dz<<std::endl;
  for(std::size_t i=0; i<n; ++i){
    TVector3 tmp = GlobalPosition(gPar, m_hit_t[i]);
    TVector3 diff = m_vp[i] - tmp;
    std::cout<<FUNC_NAME+" VP - reconstructed VP : "<<diff.Mag()<<" mm"<<std::endl;
  }
#endif

  return true;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::ResidualCheck(TVector3 pos, Double_t xzwindow, Double_t ywindow, Double_t &resi)
{

  Bool_t status = false;
  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  Double_t theta = EvalTheta(par, pos, m_min_t - 0.5*TMath::Pi(), m_max_t + 0.5*TMath::Pi());
  TVector3 fittmp = GlobalPosition(par, theta);
  TVector3 d = pos - fittmp;
  resi = d.Mag();
  Double_t xz_resi = TMath::Sqrt(d.x()*d.x()+d.z()*d.z());
  Double_t y_resi = TMath::Sqrt(d.y()*d.y());

  if(xz_resi<xzwindow && y_resi<ywindow) status = true;
  return status;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::ResidualCheck(TVector3 pos, Double_t xzwindow, Double_t ywindow)
{
  Double_t resi;
  return ResidualCheck(pos, xzwindow, ywindow, resi);
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::SortHitOrder()
{

  if(!m_is_theta_calculated) std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;

  m_hit_order.clear();
  for(std::size_t i=0; i<m_hit_t.size(); ++i){
    m_hit_order.push_back(i);
  }
  for(std::size_t i=0; i<m_hit_order.size(); ++i){
    std::sort(m_hit_order.begin(), m_hit_order.end(), CompareTheta);
  }
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DetermineCharge()
{

  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return false;
  }

  Double_t minlayer_t = 0., maxlayer_t = 0.;
  Int_t minlayer = 33, maxlayer = -1;
  if(m_hit_array.size()==0) return false;
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    Int_t layer = hitp->GetLayer();
    if(layer<minlayer){
      minlayer = layer;
      minlayer_t = m_hit_t[i];
    }
    if(layer>maxlayer){
      maxlayer = layer;
      maxlayer_t = m_hit_t[i];
    }
  }
  if(minlayer_t<maxlayer_t) m_charge = 1;
  else m_charge = -1;

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  m_edgepoint = GlobalPosition(par, maxlayer_t);

  if(m_isK18) m_charge = -1; //for K1.8 tracking.
  return true;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::VertexAtTarget()
{

  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return false;
  }

  Bool_t status = false;
  if(m_closedist.Mag() < tpc::TargetVtxWindow) status = true;
  return status;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::IsBackward()
{

  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return false;
  }

  //track is starting from the target
  if(TMath::Abs(TMath::Hypot(m_cx, m_cy) - m_r) > tpc::TargetVtxWindow) return false;

  //Track exist before the target position
  //if(m_edgepoint.z() > tpc::ZTarget) return false;
  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 start_point = GlobalPosition(par, m_min_t);
  TVector3 end_point = GlobalPosition(par, m_max_t);
  TVector3 middle_point = 0.5*(start_point + end_point);
  if(start_point.z() > tpc::ZTarget || end_point.z() > tpc::ZTarget) return false;

  //upstream end of the track is within window
  Double_t temp = m_r*m_r - TMath::Power(m_cy + tpc::ZTarget + 250. , 2);
  if(temp<0) return false;
  Double_t extrap_pointX = TMath::Min(TMath::Sqrt(temp) - m_cx, -TMath::Sqrt(temp) - m_cx); // Abs(X) of extrapolated point at Z=-250.
  if(extrap_pointX > 75.) return false; //Abs(X) < 75 at Z=-250.

  start_point -= middle_point;
  end_point -= middle_point;
  if(start_point.Mag() > 100.) return false;

  //If the backward track is accidental beam, set isbeam=1 to make fitting easier.
  if(TMath::Abs(start_point.x()) < 25. && TMath::Abs(end_point.x()) < 25. && TMath::Abs(middle_point.x()) < 25.){
    if(extrap_pointX < 25. && TMath::Abs(m_dz)<0.15) SetIsBeam();
    else if(extrap_pointX < 25. && TMath::Abs(start_point.y() - end_point.y())<15) SetIsBeam();
  }

  return true;
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::IsMultiLoop()
{

  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return ;
  }

  //High pT spiral-like track
  //Pitch > NSigma * y_resolution
  //Radius < 250.mm (TPC radius)
  //loop is larger then half circle
  Double_t pitch = TMath::Abs(2.*TMath::Pi()*m_r*m_dz);
  if(!m_is_multiloop){
    if(pitch > ThetaNSigma*GetResolutionY(0) &&
       m_r < 250. &&
       (m_max_t - m_min_t) > TMath::Pi())
      m_is_multiloop = true;
#if DebugDisp
    if(m_is_multiloop) std::cout<< " Multi-loop track!!"<<std::endl;
#endif
  }

  if(m_is_multiloop){
    if(pitch < ThetaNSigma*GetResolutionY(0) || m_r > 250.) m_is_multiloop = false;
#if DebugDisp
    if(!m_is_multiloop) std::cout<< " Not Multi-loop track!!"<<std::endl;
#endif
  }
  gMultiLoop = m_is_multiloop;

}

//______________________________________________________________________________
Int_t
TPCLocalTrackHelix::FinalizeTrack(Int_t &delete_hit)
{

  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return -1;
  }

  Int_t false_layer = 0;
  Double_t Max_residual = -100.;
  if(m_hit_array.size()==0) return -1;
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    Double_t resi = 0.;
    if(!ResidualCheck(i, resi)) ++false_layer;
    if(Max_residual<resi){
      Max_residual = resi;
      delete_hit = i;
    }
  }

  SortHitOrder();
  m_path = (m_max_t - m_min_t)*sqrt(m_r*m_r*(1. + m_dz*m_dz));
  m_transverse_path = (m_max_t - m_min_t)*m_r;
  m_mom0 = CalcHelixMom(gPar, 0.);
  if(false_layer!=0 || m_chisqr > MaxChisqr) return false_layer;

  if(!DetermineCharge()) return -1;
  IsMultiLoop();

#if DebugDisp
  std::cout<<FUNC_NAME+" chisqr: "<<m_chisqr<<std::endl
	   <<", cx: "<<m_cx
	   <<", cy: "<<m_cy
	   <<", z0: "<<m_z0
	   <<", r: "<<m_r
	   <<", dz: "<<m_dz<<std::endl
	   <<", charge: "<<m_charge
	   <<", path: "<<m_path
	   <<", fabs(m_min_t-m_max_t): "<<fabs(m_min_t-m_max_t)<<std::endl;

  //warnings
  if(m_path>550.) std::cout<<FUNC_NAME+" too long track!!! : m_path="<<m_path<<std::endl;
#endif
  if(m_hit_array.size()!=m_hit_order.size()) std::cout<<FUNC_NAME+" m_hit_array.size()!=m_hit_order.size() !!!"<<std::endl;
  if(m_hit_array.size()!=m_hit_t.size()) std::cout<<FUNC_NAME+" m_hit_array.size()!=m_hit_t.size() !!!"<<std::endl;

  return false_layer;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::ConvertParam(Double_t *linear_par)
{

  gPar[0] = 0;
  gPar[1] = 0;
  gPar[2] = linear_par[1];
  gPar[3] = 0;
  gPar[4] = linear_par[3];

  gMultiLoop = false;
  gMomConstraint = false; //No momentum constraint

  if(!DoPreFit(gPar)) return false;

  //Vtx in the target, need to check whether two tracks are merged or not
  if(!SeparateTracksAtTarget()) return false;

#if DebugDisp
  std::cout<<FUNC_NAME+" Converted track"<<std::endl
	   <<" chisqr: "<<gChisqr<<std::endl
	   <<" cx: "<<gPar[0]
	   <<", cy: "<<gPar[1]
	   <<", z0: "<<gPar[2]
	   <<", r: "<<gPar[3]
	   <<", dz: "<<gPar[4]<<std::endl;
  std::cout<<"converted track's theta min: "<<m_min_t<<" max: "<<m_max_t<<" size: "<<m_hit_t.size()<<std::endl;
#endif

  return true;
}

//______________________________________________________________________________
//If two tracks are merged at the target, separate them and recalculate params.
Bool_t
TPCLocalTrackHelix::SeparateTracksAtTarget()
{

  Bool_t status = true;
  m_vtxflag = 0;
  //If the target is at the middle of the inintial track, that track is separated with different side flag. (Two tracks are reconized as a single track)
  //Please see the Side() function.

  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);
  if(BeamThroughTPC || m_isAccidental==1){
    return status;
  }

  if(!m_is_theta_calculated) std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;

  //High pT spiral-like track
  if(m_is_multiloop && (m_max_t - m_min_t) > 2.*TMath::Pi()) return status;

  CalcClosestDistTgt(); //Distance between the target & the track
  std::vector<Int_t> side1_hits; std::vector<Int_t> side2_hits;
  if(VertexAtTarget()){ //If track is crossing the target.
    //Exclude the beam hit from the scattered track of the other scattered track's hit.
    const std::size_t n = m_hit_array.size();
    for(std::size_t i=0; i<n; ++i){
      TPCLTrackHit *hitp = m_hit_array[i];
      TVector3 pos = hitp -> GetLocalHitPos();
      if(Side(pos)==1) side1_hits.push_back(i);
      else if(Side(pos)==-1) side2_hits.push_back(i);
    }
  }
  else{ //If track is not crossing the target.
    SortHitOrder();

    //Exclude the beam hit from the scattered track.
    Bool_t flag = false;
    TVector3 prev_pos;
    Bool_t prev_isBeamHit = false; Bool_t isBeamHit;
    const std::size_t n = m_hit_array.size();
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
      if((m_dz > 0.01 || m_r < 3300.) && i!=0 && gap.Mag() > 30. && prev_isBeamHit!=isBeamHit) flag = true;
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

    Double_t par[5]={0};
    status = DoPreFit(par);
    if(status){
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
TPCLocalTrackHelix::SeparateClustersWithGap()
{

  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return false;
  }
  SortHitOrder();

  //Check a gap between clusters.
  Bool_t flag = false;
  TVector3 prev_pos;
  std::vector<Int_t> side1_hits; std::vector<Int_t> side2_hits;
  const std::size_t n = m_hit_array.size();
  for(Int_t i=0; i<n; ++i){
    Int_t id = m_hit_order[i];
    TPCLTrackHit *hitp = m_hit_array[id];
    TVector3 pos = hitp -> GetLocalHitPos();
    TVector3 gap = pos - prev_pos;
    //std::cout<<i<<" gap "<<gap.Mag()<<" pos "<<pos<<std::endl;
    if(i!=0 && gap.Mag() > MaxGapBtwClusters) flag = true;
    if(!flag) side1_hits.push_back(id);
    else side2_hits.push_back(id);
    prev_pos = pos;
  }

  if(side1_hits.size()!=0 && side2_hits.size()!=0){
    /*
    std::cout<<"size "<< side1_hits.size()<<" "<<side2_hits.size()<<std::endl;
    for(Int_t i=0; i<side1_hits.size(); ++i){
      std::cout<<i<<"/"<<side1_hits.size()<<" "<<m_hit_array[side1_hits[i]] -> GetLocalHitPos()<<std::endl;
    }
    for(Int_t i=0; i<side2_hits.size(); ++i){
      std::cout<<i<<"/"<<side2_hits.size()<<" "<<m_hit_array[side2_hits[i]] -> GetLocalHitPos()<<std::endl;
    }
    */
    if(side1_hits.size() >= side2_hits.size()) EraseHits(side2_hits);
    else EraseHits(side1_hits);
  }

  gHitPos.clear();
  return flag;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::IsKuramaTrackCandidate(Int_t id){

  Bool_t status = true;
  if(find(m_kuramaid_candidate.begin(), m_kuramaid_candidate.end(), id) == m_kuramaid_candidate.end()) status = false;
  return status;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::TestMergedTrack()
{

#if DebugDisp
  std::cout<<FUNC_NAME+" longer track's Helix params"<<std::endl
	   <<" cx: "<<m_cx
	   <<", cy: "<<m_cy
	   <<", z0: "<<m_z0
	   <<", r: "<<m_r
	   <<", dz: "<<m_dz<<std::endl;
#endif

  //Initialization of helix params
  gPar[0] = m_cx;
  gPar[1] = m_cy;
  gPar[2] = m_z0;
  gPar[3] = m_r;
  gPar[4] = m_dz;

  Bool_t vetoBadClusters = false;

  DeleteNullHit();
  const Int_t n = m_hit_array.size();
  if(!IsGoodForTracking()) return false;
  gMultiLoop = m_is_multiloop;
  gNumOfHits = n;
  gBadHits = 0;
  gMinuitStatus = 0;
  gHitPos.clear();
  gLayer.clear();
  gPadTheta.clear();
  gResParam.clear();
  gRes.clear();
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
    Double_t dummy_theta = 0;
    TVector3 res = CalcResolution(gPar, layer, pos, padTheta, dummy_theta, resparam, vetoBadClusters);
    gRes.push_back(res);
  }
  CalcHelixTheta();

  Int_t ndf;
  gChisqr = CalcChi2(gPar, ndf, vetoBadClusters);

  //pre fitting
  if(!DoPreFit(gPar)) return false;

  //Helix fitting
  if(!DoHelixFit(gPar, vetoBadClusters)) return false;
#if DebugDisp
  std::cout<<FUNC_NAME+" Helix fitting of the merged track"<<std::endl
	   <<" n_iteration: "<<m_n_iteration
	   <<", chisqr: "<<gChisqr<<std::endl
	   <<" # of clusters: "<<GetNHit()
	   <<", cx: "<<gPar[0]
	   <<", cy: "<<gPar[1]
	   <<", z0: "<<gPar[2]
	   <<", r: "<<gPar[3]
	   <<", dz: "<<gPar[4]<<std::endl;
#endif

  Int_t delete_hit = -1;
  Int_t false_layer = FinalizeTrack(delete_hit);
#if DebugDisp
  std::cout<<FUNC_NAME+" # of bad clusters : "<<false_layer<<std::endl;
#endif

  //if(false_layer > 3) return false;
  if(false_layer > 3 && (Double_t) false_layer/GetNHit() > 0.2) return false;

  //Check whether the track passing the target or not
  CalcClosestDistTgt(); //Distance between the target & the track
  CheckIsAccidental(); //check whether it is accidental beam

  if(VertexAtTarget()){
    //case1. The merged track is not crossing the target. (no problem)
    //case2. The merged track is an accidental beam crossing the target. (no problem)
    //case3. The merged track is wrongly merged two individual scattered tracks near the target <- veto

    std::vector<Int_t> side1_hits; std::vector<Int_t> side2_hits;
    for(Int_t i=0; i<n; ++i){
      Double_t resi;
      if(!ResidualCheck(i, resi)) continue;
      TPCLTrackHit *hitp = m_hit_array[i];
      TVector3 pos = hitp -> GetLocalHitPos();
      if(Side(pos)==1) side1_hits.push_back(i);
      else if(Side(pos)==-1) side2_hits.push_back(i);
    }

    if(side1_hits.size() > 0 && side2_hits.size() > 0){
      if(m_isAccidental!=1 || m_isBeam!=1) return false; //(case3) if it is not accidental beam
      /*
      if(TMath::Abs(m_closedist.x())<15. &&
	 TMath::Abs(m_closedist.y())<10. &&
	 TMath::Abs(m_closedist.z())<10.){
	 if(m_isAccidental!=1 || m_isBeam!=1) return false; //(case3) if it is not accidental beam
      }
      */
    }
  }

  //Set the vtxflag 0 to make SeparateTracksAtTarget() off.
  m_vtxflag=0;

#if IterativeResolution
  vetoBadClusters = true;

  //Now excluding bad clusters and fitting again.
  DoHelixFit(gPar, vetoBadClusters);
#endif
  return true;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoFitTrackwVertex(TVector3 vertex_pos, TVector3 vertex_res)
{

  Bool_t status = false;
  if(!IsGoodForTracking()) return status;

#if DebugDisp
  std::cout<<FUNC_NAME+" Before vertex constraint fitting, track's Helix params"<<std::endl
	   <<" cx: "<<m_cx
	   <<", cy: "<<m_cy
	   <<", z0: "<<m_z0
	   <<", r: "<<m_r
	   <<", dz: "<<m_dz<<std::endl;
#endif

  //Initialization of helix params
  gPar[0] = m_cx;
  gPar[1] = m_cy;
  gPar[2] = m_z0;
  gPar[3] = m_r;
  gPar[4] = m_dz;

  DeleteNullHit();
  const Int_t n = m_hit_array.size();

  //Vertex information
  gVertex = vertex_pos;
  gVertexRes = vertex_res;

  gMultiLoop = m_is_multiloop;
  gNumOfHits = 0;
  gBadHits = 0;
  gMinuitStatus = 0;
  gHitPos.clear();
  gLayer.clear();
  gPadTheta.clear();
  gResParam.clear();
  gRes.clear();
  for(Int_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    gNumOfHits++;
    TVector3 pos = hitp->GetLocalHitPos();
    gHitPos.push_back(pos);
    Int_t layer = hitp->GetLayer();
    gLayer.push_back(layer);
    Double_t padTheta = hitp->GetPadTheta();
    gPadTheta.push_back(padTheta);
    std::vector<Double_t> resparam = hitp->GetResolutionParams();
    gResParam.push_back(resparam);
    const TVector3& res = hitp->GetResolutionVect();
    gRes.push_back(res);
  }

  CalcHelixTheta();
  if(HelixFitwVertex()){
    status = true;
    SetParam(gPar);
    CalcHelixTheta();
    Bool_t vetoBadClusters = true;
    Int_t ndf;
    m_chisqr = CalcChi2(gPar, ndf, vetoBadClusters);
    m_minuit = gMinuitStatus;
    Double_t window = ThetaWindow/gPar[3];
    for(Int_t i=0; i<gNumOfHits; ++i){
      Double_t theta = EvalTheta(gPar, gHitPos[i], gHelixTheta[i] - 0.5*window, gHelixTheta[i] + 0.5*window);
      m_hit_t[i] = theta;
      gHelixTheta[i] = theta;
    }

    Int_t delete_hit = -1;
    FinalizeTrack(delete_hit);
  }

  return status;
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::RecalcTrack()
{

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  m_is_theta_calculated = true;
  m_is_calculated = false;
  m_mom0 = CalcHelixMom(par, 0.);
  Calculate();
  m_is_fitted = true;
  m_min_t = m_hit_t[0];
  m_max_t = m_hit_t[m_hit_t.size() - 1];
  CalcClosestDistTgt();
  m_path = (m_max_t - m_min_t)*sqrt(m_r*m_r*(1. + m_dz*m_dz));
  m_transverse_path = (m_max_t - m_min_t)*m_r;
  IsMultiLoop();

}

//_____________________________________________________________________________
void
TPCLocalTrackHelix::CheckIsAccidental()
{

  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return;
  }
  //if(!m_is_fitted) return;
  if(m_is_multiloop) return; //A multiloop track is obviously not accidental beam

  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);
  //if(BeamThroughTPC || m_isAccidental==1) return;
  if(BeamThroughTPC) return;
  m_isAccidental=0;

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 start_point = GlobalPosition(par, m_min_t);
  TVector3 end_point = GlobalPosition(par, m_max_t);
  TVector3 middle_point = 0.5*(start_point + end_point);
  //check the track is going straight twoward along the Z direction or not.
  if(TMath::Abs(start_point.x() - end_point.x()) > 200.) return;
  if(TMath::Abs(start_point.z() - middle_point.z()) < 100. ||
     TMath::Abs(middle_point.z() - end_point.z()) < 100.) return;
  if(TMath::Abs(start_point.z() - end_point.z()) < 250.) return;

  Int_t nhit_beamsection = 0;
  Int_t nhit_upstream_tgt = 0; Int_t nhit_downstream_tgt = 0;
  const std::size_t n = m_hit_array.size();
  for(Int_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp -> GetLocalHitPos();
    if(TMath::Abs(pos.x()) < 25. && pos.z() < tpc::ZTarget) nhit_beamsection++;
    if(pos.z() < tpc::ZTarget) nhit_upstream_tgt++;
    if(pos.z() > tpc::ZTarget) nhit_downstream_tgt++;
  }

  //if(nhit_upstream_tgt>=2 && nhit_downstream_tgt>=5 && m_r>3000){ //accidental beam
  if(nhit_upstream_tgt>=1 && nhit_downstream_tgt>=5 && m_r>3000){ //accidental beam
    if(TMath::Abs(m_dz)<0.05) m_isBeam=1;
    //if m_dz is large, maybe it's scattered accidental beam
    m_isAccidental=1;
  }
  else if(nhit_beamsection>=2 && nhit_downstream_tgt>=5 && m_r>1200){
    m_isAccidental=1; //scattered accidental beam crossing near the target
  }
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::TestInvertCharge()
{

  if(m_is_multiloop || m_isBeam==1 || m_isK18==1 || m_isAccidental==1 || m_isXi==1) return false;

  if(!m_is_theta_calculated){
    std::cout<<FUNC_NAME+" Fatal error : No helix theta information!!! CalcHelixTheta() should be run in front of this"<<std::endl;
    return false;
  }

#if DebugDisp
  std::cout<<FUNC_NAME+" track for testing"<<std::endl
	   <<" cx: "<<m_cx
	   <<", cy: "<<m_cy
	   <<", z0: "<<m_z0
	   <<", r: "<<m_r
	   <<", dz: "<<m_dz<<std::endl;
#endif

  //Initialization of helix params
  gPar[0] = m_cx;
  gPar[1] = m_cy;
  gPar[2] = m_z0;
  gPar[3] = m_r;
  gPar[4] = m_dz;

  Bool_t vetoBadClusters = false;
  Int_t prev_charge = m_charge;

  DeleteNullHit();
  const Int_t n = m_hit_array.size();
  if(!IsGoodForTracking()) return false;
  gMultiLoop = m_is_multiloop;
  gNumOfHits = n;
  gBadHits = 0;
  gMinuitStatus = 0;
  gHitPos.clear();
  gLayer.clear();
  gPadTheta.clear();
  gResParam.clear();
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
  }

  TVector3 localpos_start = GlobalToLocal(gHitPos[0]);
  TVector3 localpos_end = GlobalToLocal(gHitPos[n-1]);
  TVector3 localpos_avg = 0.5*(localpos_start + localpos_end);
  Double_t x_avg = localpos_avg.x();
  Double_t y_avg = localpos_avg.y();
  Double_t z_avg = localpos_avg.z();
  gPar[0] = 2.*x_avg - m_cx;
  gPar[1] = 2.*y_avg - m_cy;
  Double_t tmp_theta = TMath::ATan2(y_avg - gPar[1], x_avg - gPar[0]);
  gPar[2] = z_avg - tmp_theta*m_r*(-m_dz);
  gPar[3] = m_r;
  gPar[4] = -m_dz;
  SetParam(gPar);
  CalcHelixTheta();

  //if(!DoStraightLineFit(gPar)) return false;
  //SetParam(gPar);

#if DebugDisp
  std::cout<<FUNC_NAME+" before helixfit for inverting "<<std::endl
	   <<" cx: "<<m_cx
	   <<", cy: "<<m_cy
	   <<", z0: "<<m_z0
	   <<", r: "<<m_r
	   <<", dz: "<<m_dz<<std::endl;
#endif

  gRes.clear();
  for(Int_t i=0; i<n; ++i){
    TVector3 res = CalcResolution(gPar, gLayer[i], gHitPos[i], gPadTheta[i], gHelixTheta[i], gResParam[i], vetoBadClusters);
    gRes.push_back(res);
  }

  Bool_t status = false;
  if(HelixFitInvertCharge()){
    status = true;
    SetParam(gPar);
    CalcHelixTheta();

    Double_t window = ThetaWindow/gPar[3];
    for(Int_t i=0; i<gNumOfHits; ++i){
      Double_t theta = EvalTheta(gPar, gHitPos[i], gHelixTheta[i] - 0.5*window, gHelixTheta[i] + 0.5*window);
      m_hit_t[i] = theta;
      gHelixTheta[i] = theta;
      gRes[i] = CalcResolution(gPar, gLayer[i], gHitPos[i], gPadTheta[i], gHelixTheta[i], gResParam[i], vetoBadClusters);
    }

    Int_t ndf;
    m_chisqr = CalcChi2(gPar, ndf, vetoBadClusters);
    m_minuit = gMinuitStatus;
  }
  else return false;

#if DebugDisp
  std::cout<<FUNC_NAME+" after helixfit for inverting "<<std::endl
	   <<" cx: "<<m_cx
	   <<", cy: "<<m_cy
	   <<", z0: "<<m_z0
	   <<", r: "<<m_r
	   <<", dz: "<<m_dz<<std::endl;
#endif

  Int_t delete_hit = -1;
  Int_t false_layer = FinalizeTrack(delete_hit);
  if(!DetermineCharge() || prev_charge == m_charge) return false; //charge is not converted

#if DebugDisp
  std::cout<<FUNC_NAME+" # of bad clusters : "<<false_layer<<std::endl;
#endif

  if(false_layer > 0){
    Int_t maxloop = 0;
    while(false_layer > 0 && maxloop < 30){
      if(!HelixFitInvertCharge()) return false;
      else{
	Double_t window = ThetaWindow/gPar[3];
	for(Int_t i=0; i<gNumOfHits; ++i){
	  Double_t theta = EvalTheta(gPar, gHitPos[i], gHelixTheta[i] - 0.5*window, gHelixTheta[i] + 0.5*window);
	  m_hit_t[i] = theta;
	  gHelixTheta[i] = theta;
	  gRes[i] = CalcResolution(gPar, gLayer[i], gHitPos[i], gPadTheta[i], gHelixTheta[i], gResParam[i], vetoBadClusters);
	}

	Int_t ndf;
	Double_t chisqr = CalcChi2(gPar, ndf, vetoBadClusters);
	if(chisqr >= m_chisqr) break;
	SetParam(gPar);
	CalcHelixTheta();
	m_chisqr = CalcChi2(gPar, ndf, vetoBadClusters);
	m_minuit = gMinuitStatus;
	false_layer = FinalizeTrack(delete_hit);
      }
      maxloop++;
    } //while

    if(false_layer > 2) return false;
    if(!DetermineCharge() || prev_charge == m_charge) return false;
  }

  //Check whether the track passing the target or not
  CalcClosestDistTgt(); //Distance between the target & the track
  CheckIsAccidental(); //check whether it is accidental beam
  VertexAtTarget();

#if IterativeResolution
  vetoBadClusters = true;

  for(Int_t i=0; i<gNumOfHits; ++i) gRes[i] = CalcResolution(gPar, gLayer[i], gHitPos[i], gPadTheta[i], gHelixTheta[i], gResParam[i], vetoBadClusters);

  //Now excluding bad clusters and fitting again.
  HelixFitInvertCharge();
  SetParam(gPar);
  CalcHelixTheta();
  m_chisqr = gChisqr;
  m_minuit = gMinuitStatus;
  Double_t window = ThetaWindow/gPar[3];
  for(Int_t i=0; i<gNumOfHits; ++i){
    Double_t theta = EvalTheta(gPar, gHitPos[i], gHelixTheta[i] - 0.5*window, gHelixTheta[i] + 0.5*window);
    m_hit_t[i] = theta;
    gHelixTheta[i] = theta;
  }

#if DebugDisp
  std::cout<<FUNC_NAME+" after refitting helixfit for inverting "<<std::endl
	   <<" cx: "<<m_cx
	   <<", cy: "<<m_cy
	   <<", z0: "<<m_z0
	   <<", r: "<<m_r
	   <<", dz: "<<m_dz<<std::endl;
  std::cout<<"prev current charge "<<prev_charge<<" "<<m_charge<<std::endl;
#endif
#endif

  return true;
}

//______________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetClosestPositionTgt()
{
  TVector3 tgt(0., 0., tpc::ZTarget);
  TVector3 pos = m_closedist + tgt;
  return pos;
}

//______________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetClosestPositionTgtXZ()
{
  TVector3 tgt(0., 0., tpc::ZTarget);
  TVector3 pos = m_closedistXZ + tgt;
  return pos;
}

//Functions for momentum resolution
//_____________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetMomentumResolutionVectT(Double_t t, Double_t MomScale, Double_t PhiScale, Double_t dZScale){

  /*
    For a Multivaraible function F(M), the Covariance is given as V(F) = JV(M)J^T, where J = dFdM is a Jacobian matrix.
    We will calculate the momentum resolution from pt, theta, and dZ resolution.
  */

  Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  Double_t p_t = m_r*(tpc::ConstC*dMagneticField)*0.001;
  //Double_t pz = p_t*(cos(t));
  //Double_t py = p_t*m_dz;
  //Double_t px = p_t*(sin(t));

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 calcmom = CalcHelixMom(par, t);
  Double_t ddZ = dZScale * GetdZResolution();
  Double_t dt = PhiScale * GetTransverseAngularResolution(t);
  Double_t dp_t = MomScale * GetTransverseMomentumResolution();
  Double_t Vp_t = dp_t*dp_t;
  Double_t Vt = dt*dt;
  Double_t VdZ = ddZ*ddZ;
  /*
    p_x = p_t*sin(t);
    p_y = p_t*dZ;
    p_z = p_t*cos(t);

    Cov(px,pz,py) = J V(p_t,t,dZ)J^T
    **Note! Cov(pt,t) = res_pt*res_t, because res_t ~ res_pt!

             f
		 J = x| dx/df|

		    p_t    t     dZ
    px| sin  pt*cos   0 |;
  J=py|  dZ    0     pt |;
    pz| cos -pt*sin   0 |;
 */

  Double_t cov_pt_t = MomScale * PhiScale * GetTransverseMomentumAngularCovariance(t);
  Double_t V[3][3] =
    { Vp_t,   cov_pt_t,    0,
      cov_pt_t, Vt,        0,
      0,         0,       VdZ};
  Double_t J[3][3] =
    { sin(t),   p_t*cos(t), 0,
      m_dz,     0,          p_t,
      cos(t), -p_t*sin(t), 0};

  Double_t JT[3][3]={0};
  for(Int_t row=0;row<3;++row){
    for(Int_t col=0;col<3;++col){
      JT[row][col] = J[col][row];
    }
  }
  Double_t VJT[3][3]={0};
  for(Int_t row=0;row<3;++row){
    for(Int_t col=0;col<3;++col){
      for(Int_t itr=0;itr<3;++itr){
	VJT[row][col]+=V[row][itr]*JT[itr][col];
      }
    }
  }
  Double_t JVJT[3][3]={0};
  for(Int_t row=0;row<3;++row){
    for(Int_t col=0;col<3;++col){
      for(Int_t itr=0;itr<3;++itr){
	JVJT[row][col]+=J[row][itr]*VJT[itr][col];
      }
    }
  }
  Double_t Vpx = JVJT[0][0];
  Double_t Vpy = JVJT[1][1];
  Double_t Vpz = JVJT[2][2];
  if(Vpx<0 or Vpy<0 or Vpz < 0 or std::isnan(Vpx) or std::isnan(Vpy) or std::isnan(Vpz)){
    std::cout<<Form("MomVar = (%g,%g,%g)",Vpx,Vpy,Vpz)<<std::endl;
    std::cout<<Form("dPt,dt,ddZ = (%g,%g,%g)",dp_t,dt,ddZ)<<std::endl;
  }
  return TVector3(sqrt(Vpx), sqrt(Vpy), sqrt(Vpz));
}

//_____________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetMomentumResolutionVect(Int_t i, Double_t MomScale, Double_t PhiScale, Double_t dZScale){
  Double_t t = GetHitInOrder(i) -> GetTheta();
  return GetMomentumResolutionVectT(t, MomScale, PhiScale, dZScale);
}

//_____________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetMomentumResolutionVect(){

  Double_t t0 = GetHitInOrder(0) -> GetTheta();
  return GetMomentumResolutionVectT(t0, m_MomResScale, m_PhResScale, m_dZResScale);
}

//_____________________________________________________________________________
double
TPCLocalTrackHelix::GetTransverseMomentumAngularCovariance(Double_t t){


  if(t == -9999){
    t = GetHitInOrder(0)->GetTheta();
  }

  Double_t sign = 1;
  if(m_charge>0)sign = -1;
  Double_t dp_t =  GetTransverseMomentumResolution();
  Double_t dt = GetTransverseAngularResolution(t,0);
  return sign * dp_t * dt;
}

//_____________________________________________________________________________
double
TPCLocalTrackHelix::GetMomentumPitchAngleCovariance(){
  Double_t B = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  Double_t p_t = m_r*(tpc::ConstC*B)*0.001;
  Double_t pitch = atan2(1,m_dz);//dYdZ angle, dZ= 0 -> should return pi/2
  Double_t res_pitch = GetThetaResolution();
  return p_t*cos(pitch)/sin(pitch)/sin(pitch)*res_pitch*res_pitch;

}

//_____________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetMomentumCovarianceVectT(Double_t t, Double_t MomScale, Double_t PhiScale, Double_t dZScale){

  /*
    For a Multivaraible function F(M), the Covariance is given as V(F) = JV(M)J^T, where J = dFdM is a Jacobian matrix.
    We will calculate the momentum resolution from pt, theta, and dZ resolution.
  */
  Double_t B = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  Double_t p_t = m_r*(tpc::ConstC*B)*0.001;
  //Double_t pz = p_t*(cos(t));
  //Double_t py = p_t*m_dz;
  //Double_t px = p_t*(sin(t));

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 calcmom = CalcHelixMom(par, t);
  Double_t ddZ = dZScale * GetdZResolution();
  Double_t dt = PhiScale * GetTransverseAngularResolution(t);
  Double_t dp_t = MomScale * GetTransverseMomentumResolution();
  Double_t Vp_t = dp_t*dp_t;
  Double_t Vt = dt*dt;
  Double_t VdZ = ddZ*ddZ;
  /*
    p_x = p_t*sin(t);
    p_y = p_t*dZ;
    p_z = p_t*cos(t);

    Cov(px,pz,py) = J V(p_t,t,dZ)J^T
    **Note! Cov(pt,t) = res_pt*res_t, because res_t ~ res_pt!
  */
  Double_t cov_pt_t = MomScale*PhiScale*GetTransverseMomentumAngularCovariance(t);
  Double_t V[3][3] =
    { Vp_t,     cov_pt_t,     0,
      cov_pt_t, Vt,           0,
      0,        0,            VdZ};
  Double_t J[3][3] =
    { sin(t),   p_t*cos(t),   0,
      m_dz,     0,            p_t,
      cos(t),   -p_t*sin(t),  0};

  Double_t JT[3][3]={0};
  for(Int_t row=0;row<3;++row){
    for(Int_t col=0;col<3;++col){
      JT[row][col] = J[col][row];
    }
  }
  Double_t VJT[3][3]={0};
  for(Int_t row=0;row<3;++row){
    for(Int_t col=0;col<3;++col){
      for(Int_t itr=0;itr<3;++itr){
	VJT[row][col]+=V[row][itr]*JT[itr][col];
      }
    }
  }
  Double_t JVJT[3][3]={0};
  for(Int_t row=0;row<3;++row){
    for(Int_t col=0;col<3;++col){
      for(Int_t itr=0;itr<3;++itr){
	JVJT[row][col]+=J[row][itr]*VJT[itr][col];
      }
    }
  }
  Double_t Cxy = JVJT[0][1];
  Double_t Cyz = JVJT[1][2];
  Double_t Czx = JVJT[2][0];
  return TVector3((Cxy), (Cyz), (Czx));
}

//_____________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetMomentumCovarianceVect(Int_t i, Double_t MomScale, Double_t PhiScale, Double_t dZScale){

  Double_t t = GetHitInOrder(i) -> GetTheta();
  return GetMomentumCovarianceVectT(t, MomScale, PhiScale, dZScale);
}

//_____________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetMomentumCovarianceVect(){

  Double_t t0 = GetHitInOrder(0) -> GetTheta();
  return GetMomentumCovarianceVectT(t0, m_MomResScale, m_PhResScale, m_dZResScale);

}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetTransverseMomentumResolution(){
  Double_t B = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  Double_t dt = abs(m_max_t - m_min_t);
  if(dt > 2*acos(-1) )dt = 2*acos(-1);

  Double_t L = 2 * sin(0.5*dt)*m_r;//String length, not Arc length
  L*=0.001;//mm->m;
  Double_t res = 0;
  Int_t nh = m_hit_array.size();
  for(Int_t ih=0;ih<m_hit_array.size();++ih){
    Int_t id = m_hit_order[ih];
    TPCLTrackHit *hitp = m_hit_array[id];
    auto ResV = hitp -> GetResolutionVect();
    Double_t res_T = hypot(ResV.X(),ResV.Z());
    if(!hitp -> IsGoodForTracking()){
      nh--;
      continue;
    }
    res+=res_T*res_T;
  }
  Double_t pt = m_r*(tpc::ConstC*B)*0.001;
  if(nh<4) return pt*0.1;
  res = sqrt(3./2) * sqrt(res / nh)* 0.001;//mm-> m
  Double_t dPOverP = pt / (0.3*L*L*B)*sqrt(720./(nh+4))*res;
  if(std::isnan(dPOverP) || dPOverP < 0){
    std::cout<<Form("dPt error! nh = %d, res = %g",nh,res)<<std::endl;
  }
  return m_MomResScale*pt*dPOverP;
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetTransverseAngularResolution(Double_t t, Double_t sig0){
  //Transverse Angle Definition: atan2(pz,px);
  Double_t B = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  Double_t pt = m_r*(tpc::ConstC*B)*0.001;
  Double_t dp = GetTransverseMomentumResolution();
  Double_t dr = m_r * dp/pt;
  Double_t t_avg = 0.5*(m_max_t + m_min_t);
  Double_t dt = (t-t_avg);
  if(dt >acos(-1))dt = acos(-1);
  Double_t path = m_r * abs(dt);

  return hypot(m_PhResScale*path * dr / m_r / m_r,sig0);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetTransverseAngularResolution(){

  Double_t t0 = GetHitInOrder(0) -> GetTheta();
  Double_t sig0 = 0.001;
  return GetTransverseAngularResolution(t0, sig0);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetdZResolution(){
  Double_t res2 = 0;
  Int_t nh = m_hit_array.size();
  for(Int_t ih=0;ih<m_hit_array.size();++ih){
    Int_t id = m_hit_order[ih];
    TPCLTrackHit *hitp = m_hit_array[id];
    auto ResV = hitp -> GetResolutionVect();
    Double_t res_Y = ResV.Y();
    if(!hitp -> IsGoodForTracking()){
      nh--;
      continue;
    }
    res2 += res_Y*res_Y;
  }
  Double_t dt = abs(m_max_t - m_min_t);
  //  if(dt > 2*acos(-1)) dt = 2*acos(-1);
  Double_t path = m_r * dt;
  Double_t path_dev = path*path*nh/12;
  if(nh < 3) return 0.01;
  Double_t d_slope = 1./(nh-2)*res2/path_dev;
  if(std::isnan(d_slope) || std::isinf(d_slope)){
    std::cout<<Form("Nan || inf dZ resol! nh = %d, dt = %g, res = %g", nh,dt,res2)<<std::endl;
  }
  return m_dZResScale*sqrt(d_slope);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetThetaResolution(){
  Double_t d_slope = GetdZResolution();
  return d_slope / (1+m_dz*m_dz);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetMomentumResolution(){
  Double_t d_slope = GetdZResolution();
  return GetTransverseMomentumResolution()*hypot(1,d_slope);
}
TMatrixD
TPCLocalTrackHelix::GetCovarianceMatrix(){
  double Elements[3*3]={0};
  double cov_mom_th = GetMomentumPitchAngleCovariance();
  double cov_mom_ph = GetTransverseMomentumAngularCovariance();
  double res_mom = GetMomentumResolution();
  double res_th = GetThetaResolution();
  double res_ph = GetTransverseAngularResolution();
  Elements[0+3*0] = res_mom*res_mom;
  Elements[1+3*1] = res_th*res_th;
  Elements[2+3*2] = res_ph*res_ph;
  Elements[0*3+1] = cov_mom_th;
  Elements[1*3+0] = cov_mom_th;
  Elements[0*3+2] = cov_mom_ph;
  Elements[2*3+0] = cov_mom_ph;
  TMatrixD CovMat(3,3,Elements);
  return CovMat;
}

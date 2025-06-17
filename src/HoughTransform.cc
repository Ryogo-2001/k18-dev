// -*- C++ -*-

#include "TPCTrackSearch.hh"

#include <chrono>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <TH2D.h>
#include <TH3D.h>

#include "DebugTimer.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "DeleteUtility.hh"
#include "ConfMan.hh"
#include "TPCPadHelper.hh"
#include "TPCCluster.hh"

#include "RootHelper.hh"

#define DebugDisp      0
#define DebugEvDisp    0

#if DebugEvDisp
#include <TApplication.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
namespace { TApplication app("DebugApp", nullptr, nullptr); }
static Int_t cannum = 0;
#endif
namespace
{
  const auto qnan = TMath::QuietNaN();
  const Double_t& HSfield_Calib = ConfMan::Get<Double_t>("HSFLDCALIB");
  const Double_t& HSfield_Calc = ConfMan::Get<Double_t>("HSFLDCALC");
  const Double_t& HSfield_Hall = ConfMan::Get<Double_t>("HSFLDHALL");
  const Double_t Const = 0.299792458;

  // HoughXZ binning
  //static const Int_t    Li_theta_ndiv = 180;
  static const Int_t    Li_theta_ndiv = 360;
  static const Double_t Li_theta_min  =   0;
  static const Double_t Li_theta_max  = 180;
  //static const Int_t    Li_rho_ndiv =  180;
  static const Int_t    Li_rho_ndiv =  720;
  static const Double_t Li_rho_min  = -720;
  static const Double_t Li_rho_max  =  720;
  TH2D *histLinear = new TH2D("histLinear", "; #theta (deg.); #rho (mm)",
			      Li_theta_ndiv, Li_theta_min, Li_theta_max,
			      Li_rho_ndiv, Li_rho_min, Li_rho_max);

  // HoughCircleXZ binning
  const Int_t    nBin_rdiff = 11;
  const Double_t rdiff_min = -110.;
  const Double_t rdiff_max = 110.;
  const Int_t    nBin_theta = 720;
  const Double_t theta_min = -1.*acos(-1);
  const Double_t theta_max = acos(-1);
  const Int_t    nBin_p = 900;
  //const Double_t pmin = 50.;
  const Double_t pmin = 20.;
  const Double_t pmax = 1550.; //MeV/c

  /*
  const Int_t    nBin_rdiff = 24;
  const Double_t rdiff_min = -120.;
  const Double_t rdiff_max = 120.;
  const Int_t    nBin_theta = 36;
  const Double_t theta_min = -1.*acos(-1);
  const Double_t theta_max = acos(-1);
  const Int_t    nBin_p = 500;
  const Double_t pmin = 50.;
  const Double_t pmax = 1550.; //MeV/c
  */
  TH3D *histCircle = new TH3D("histCircle",";rd (mm); theta (rad); p(MeV/c)",
			       nBin_rdiff, rdiff_min,  rdiff_max,
			       nBin_theta, theta_min, theta_max,
			       nBin_p, pmin, pmax);

  // HoughY binning
  //const Ynt_t    thetaY_ndiv = 360; //previous value
  const Int_t    thetaY_ndiv = 1080;
  const Double_t thetaY_min  =   0.;
  const Double_t thetaY_max  = 180.;
  const Int_t    r_ndiv =  2000;
  const Double_t r_min  = -5000.;
  const Double_t r_max  =  5000.;
  TH2D *histY = new TH2D("histY",";theta (deg.);r (mm)",
			 thetaY_ndiv, thetaY_min, thetaY_max,
			 r_ndiv, r_min, r_max);
}

//_____________________________________________________________________________
static inline void
HoughFlagCheck(const std::vector<TPCClusterContainer>& ClCont, std::vector<TVector3>& gHitPos){

  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
      auto cl = ClCont[layer][ci];
      TPCHit* hit = cl->GetMeanHit();
      if(hit->GetHoughFlag()>0) continue;
      TVector3 pos = cl->GetPosition();
      gHitPos.push_back(pos);
    } //ci
  } //layer
}

//_____________________________________________________________________________
namespace tpc
{

//_____________________________________________________________________________
Bool_t
HoughTransformLineXZ(std::vector<TVector3> gHitPos, Int_t *MaxBin,
		     Double_t *LinearPar, Int_t MinNumOfHits /*=8*/)
{

#if DebugEvDisp
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  static TCanvas c1("c"+FUNC_NAME, FUNC_NAME, 900, 900);
  c1.cd();
#endif

  Bool_t status = true;

  //Hough space
  //for TPC linear track
  //rho = z*cos(theta) + x*sin(theta)
  //x = p[0] + p[2]*z
  histLinear->Reset();

  for(Int_t i=0; i<gHitPos.size(); ++i){
    for(Int_t ti=0; ti<histY -> GetNbinsX(); ti++){
      Double_t theta = histLinear->GetXaxis()->GetBinCenter(ti+1);
      Double_t rho = TMath::Cos(theta*TMath::DegToRad())*(gHitPos[i].Z() - tpc::ZTarget)
	+ TMath::Sin(theta*TMath::DegToRad())*gHitPos[i].X();
      histLinear->Fill(theta, rho);
    } //ti
  } //i

#if DebugEvDisp
  histLinear->Draw("colz");
  gPad->Modified();
  gPad->Update();
  c1.Print(Form("c%d.pdf",cannum));
  getchar();
  cannum++;
#endif

  Int_t maxbin = histLinear->GetMaximumBin();
  Int_t mx, my, mz;
  histLinear->GetBinXYZ(maxbin, mx, my, mz);
  MaxBin[0] = mx; MaxBin[1] = my;

  Double_t mtheta = histLinear->GetXaxis()->GetBinCenter(mx)*TMath::DegToRad();
  Double_t mr = histLinear->GetYaxis()->GetBinCenter(my);

  LinearPar[0] = mr/TMath::Sin(mtheta); //m_x
  LinearPar[2]  = -TMath::Cos(mtheta)/TMath::Sin(mtheta); //m_u

#if DebugDisp
  std::cout<<FUNC_NAME+" (Z, X) hough vote : "<<histLinear->GetBinContent(maxbin)<<"/"<<gHitPos.size()<<std::endl;
  std::cout<<FUNC_NAME+" Hough (theta, r)=("<<mtheta<<", "<<mr<<")"<<std::endl;
  std::cout<<FUNC_NAME+" m_x="<<LinearPar[0]<<" m_u="<<LinearPar[2]<<std::endl;
#endif

  if(histLinear->GetMaximum() < 0.5*MinNumOfHits) status = false;
  return status;

}

//_____________________________________________________________________________
Bool_t
HoughTransformLineXZ(const std::vector<TPCClusterContainer>& ClCont,
		     Int_t *MaxBin, Double_t *LinearPar,
		     Int_t MinNumOfHits /*=8*/)
{

  std::vector<TVector3> gHitPos;
  HoughFlagCheck(ClCont, gHitPos);
  return HoughTransformLineXZ(gHitPos, MaxBin, LinearPar, MinNumOfHits);

}

//_____________________________________________________________________________
void
HoughTransformLineYZ(std::vector<TVector3> gHitPos, Int_t *MaxBin,
		     Double_t *LinearPar, Double_t MaxHoughWindowY)
{

#if DebugEvDisp
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  static TCanvas c1("c"+FUNC_NAME, FUNC_NAME, 900, 900);
  c1.cd();
#endif

  //Hough space
  //for TPC linear track
  //rho = z*cos(theta) + y*sin(theta)
  //y = p[1] + p[3]*z
  histLinear->Reset();
  for(Int_t i=0; i<gHitPos.size(); ++i){
    Double_t dist = TMath::Abs(LinearPar[2]*(gHitPos[i].Z() - tpc::ZTarget) - gHitPos[i].X() +
			       LinearPar[0])/TMath::Sqrt(TMath::Sq(LinearPar[2])+1.);
    if(dist < MaxHoughWindowY){
      for(Int_t ti=0; ti<histY -> GetNbinsX(); ti++){
	Double_t theta = histLinear->GetXaxis()->GetBinCenter(ti+1);
	Double_t rho = TMath::Cos(theta*TMath::DegToRad())*(gHitPos[i].Z() - tpc::ZTarget)
	  + TMath::Sin(theta*TMath::DegToRad())*gHitPos[i].Y();
	histLinear->Fill(theta, rho);
      } //ti
    } //dist
  } //i

#if DebugEvDisp
  histLinear->Draw("colz");
  gPad->Modified();
  gPad->Update();
  c1.Print(Form("c%d.pdf",cannum));
  getchar();
  cannum++;
#endif

  Int_t maxbin = histLinear->GetMaximumBin();
  Int_t mx, my, mz;
  histLinear->GetBinXYZ(maxbin, mx, my, mz);
  MaxBin[0] = mx; MaxBin[1] = my;

  Double_t mtheta = histLinear->GetXaxis()->GetBinCenter(mx)*TMath::DegToRad();
  Double_t mr = histLinear->GetYaxis()->GetBinCenter(my);

  LinearPar[1] = mr/TMath::Sin(mtheta); //m_y
  LinearPar[3] = -TMath::Cos(mtheta)/TMath::Sin(mtheta); //m_v

#if DebugDisp
  std::cout<<FUNC_NAME+" (Z, Y) hough vote : "<<histLinear->GetBinContent(maxbin)<<"/"<<gHitPos.size()<<std::endl;
  std::cout<<FUNC_NAME+" Hough (theta, r)=("<<mtheta<<", "<<mr<<")"<<std::endl;
  std::cout<<FUNC_NAME+" m_y="<<LinearPar[1]<<" m_v="<<LinearPar[3]<<std::endl;
#endif

}

//_____________________________________________________________________________
void
HoughTransformLineYZ(const std::vector<TPCClusterContainer>& ClCont,
		     Int_t *MaxBin, Double_t *LinearPar,
		     Double_t MaxHoughWindowY)
{

  std::vector<TVector3> gHitPos;
  HoughFlagCheck(ClCont, gHitPos);
  HoughTransformLineYZ(gHitPos, MaxBin, LinearPar, MaxHoughWindowY);

}

//_____________________________________________________________________________
void
HoughTransformLineYX(std::vector<TVector3> gHitPos, Int_t *MaxBin,
		     Double_t *LinearPar, Double_t MaxHoughWindowY)
{

#if DebugEvDisp
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  static TCanvas c1("c"+FUNC_NAME, FUNC_NAME, 900, 900);
  c1.cd();
#endif

  //Hough space
  //for TPC linear track
  //rho = x*cos(theta) + y*sin(theta)
  //y = tmp0 + tmp1*x
  histLinear->Reset();
  for(Int_t i=0; i<gHitPos.size(); ++i){
    Double_t dist = TMath::Abs(LinearPar[2]*(gHitPos[i].Z() - tpc::ZTarget) - gHitPos[i].X() +
			       LinearPar[0])/TMath::Sqrt(TMath::Sq(LinearPar[2])+1.);
    if(dist < MaxHoughWindowY){
      for(Int_t ti=0; ti<histY -> GetNbinsX(); ti++){
	Double_t theta = histLinear->GetXaxis()->GetBinCenter(ti+1);
	Double_t rho = TMath::Cos(theta*TMath::DegToRad())*gHitPos[i].X()
	  + TMath::Sin(theta*TMath::DegToRad())*gHitPos[i].Y();
	histLinear->Fill(theta, rho);
      } //ti
    } //dist
  } //i

#if DebugEvDisp
  histLinear->Draw("colz");
  gPad->Modified();
  gPad->Update();
  c1.Print(Form("c%d.pdf",cannum));
  getchar();
  cannum++;
#endif

  Int_t maxbin = histLinear->GetMaximumBin();
  Int_t mx, my, mz;
  histLinear->GetBinXYZ(maxbin, mx, my, mz);
  MaxBin[0] = mx; MaxBin[1] = my;

  Double_t mtheta = histLinear->GetXaxis()->GetBinCenter(mx)*TMath::DegToRad();
  Double_t mr = histLinear->GetYaxis()->GetBinCenter(my);
  Double_t tmp0 = mr/TMath::Sin(mtheta);
  Double_t tmp1 = -TMath::Cos(mtheta)/TMath::Sin(mtheta);
  LinearPar[1] = tmp0 + LinearPar[0]*tmp1;
  LinearPar[3] = LinearPar[2]*tmp1;

#if DebugDisp
  std::cout<<FUNC_NAME+" (X, Y) hough vote : "<<histLinear->GetBinContent(maxbin)<<"/"<<gHitPos.size()<<std::endl;
  std::cout<<FUNC_NAME+" Hough (theta, r)=("<<mtheta<<", "<<mr<<")"<<std::endl;
  std::cout<<FUNC_NAME+" m_y="<<LinearPar[1]<<" m_v="<<LinearPar[3]<<std::endl;
#endif

}

//_____________________________________________________________________________
void
HoughTransformLineYX(const std::vector<TPCClusterContainer>& ClCont,
		     Int_t *MaxBin, Double_t *LinearPar,
		     Double_t MaxHoughWindowY)
{

  std::vector<TVector3> gHitPos;
  HoughFlagCheck(ClCont, gHitPos);
  HoughTransformLineYX(gHitPos, MaxBin, LinearPar, MaxHoughWindowY);

}

//_____________________________________________________________________________
void
HoughTransformLineYTheta(std::vector<TVector3> gHitPos, Int_t *MaxBin,
		       Double_t *HelixPar, Double_t MaxHoughWindowY)
{

#if DebugEvDisp
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  static TCanvas c1("c"+FUNC_NAME, FUNC_NAME, 900, 900);
  c1.cd();
#endif

  //Hough space
  //mr = tmp_xval*cos(mtheta) + tmpz*sin(mtheta)
  //tmpz = p[2] + p[4]*tmp_xval
  histY -> Reset();
  for(Int_t i=0; i<gHitPos.size(); ++i){
    Double_t tmpx = -gHitPos[i].X();
    Double_t tmpy = gHitPos[i].Z() - tpc::ZTarget;
    Double_t tmpz = gHitPos[i].Y();
    Double_t r_cal = sqrt(pow(tmpx - HelixPar[0], 2) + pow(tmpy - HelixPar[1], 2));
    Double_t dist = fabs(r_cal - HelixPar[3]);
    if(dist < MaxHoughWindowY){
      for(Int_t ti=0; ti<histY -> GetNbinsX(); ti++){
	Double_t theta = histY->GetXaxis()->GetBinCenter(ti+1);
	Double_t mtheta = theta*acos(-1.)/180.;
	Double_t tmpt = TMath::ATan2(tmpy - HelixPar[1], tmpx - HelixPar[0]);
	Double_t tmp_xval = HelixPar[3]*tmpt;
	Double_t mr = cos(mtheta)*tmp_xval + sin(mtheta)*tmpz;

	Double_t p2 = mr/sin(mtheta);
	Double_t p4 = -cos(mtheta)/sin(mtheta);
	if(TMath::Abs(p2) < 15000. && TMath::Abs(p4) < 15.) histY->Fill(theta, mr);
    } //ti
  } //dist
}

#if DebugEvDisp
  histY->Draw("colz");
  gPad->Modified();
  gPad->Update();
  c1.Print(Form("c%d.pdf",cannum));
  getchar();
  cannum++;
#endif

  Int_t maxbin = histY->GetMaximumBin();
  Int_t mx, my, mz;
  histY->GetBinXYZ(maxbin, mx, my, mz);
  MaxBin[0] = mx; MaxBin[1] = my;

  Double_t mtheta = histY->GetXaxis()->GetBinCenter(mx)*acos(-1.)/180.;
  Double_t mr = histY->GetYaxis()->GetBinCenter(my);
  HelixPar[2] = mr/sin(mtheta);
  HelixPar[4] = -cos(mtheta)/sin(mtheta);

#if DebugDisp
  std::cout<<FUNC_NAME+" (Y, phi) hough vote : "<<histY->GetBinContent(maxbin)<<"/"<<gHitPos.size()<<std::endl;
  std::cout<<FUNC_NAME+" Hough (theta, r)=("<<mtheta<<", "<<mr<<")"<<std::endl;
#endif

}

//_____________________________________________________________________________
void
HoughTransformLineYTheta(const std::vector<TPCClusterContainer>& ClCont,
		       Int_t *MaxBin, Double_t *HelixPar, Double_t MaxHoughWindowY)
{

  std::vector<TVector3> gHitPos;
  HoughFlagCheck(ClCont, gHitPos);
  HoughTransformLineYTheta(gHitPos, MaxBin, HelixPar, MaxHoughWindowY);

}

//_____________________________________________________________________________
Bool_t
HoughTransformCircleXZ(std::vector<TVector3> gHitPos,
		       Int_t *MaxBin, Double_t *HelixPar,
		       Int_t MinNumOfHits /*=8*/)
{

  // Equation
  // (x - (r + rd)*cos(theta))^2 + (y - (r + rd)*sin(theta))^2 = r^2
  // p = r * Const * dMagneticField;
  Double_t dMagneticField = HSfield_Calib*(HSfield_Hall/HSfield_Calc);
  Bool_t status = true;

#if DebugEvDisp
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  static TCanvas c1("c"+FUNC_NAME, FUNC_NAME, 900, 900);
  c1.cd();
#endif

  //for TPC circle track
  //Hough-transform
  histCircle -> Reset();
  for(Int_t i=0; i<gHitPos.size(); ++i){
    Double_t x = -gHitPos[i].X();
    Double_t y = gHitPos[i].Z() - tpc::ZTarget;
    for(Int_t ird=0; ird<histCircle -> GetNbinsX(); ++ird){
      Double_t rd = histCircle->GetXaxis()->GetBinCenter(ird+1);
      for(Int_t ip=0; ip<histCircle -> GetNbinsZ(); ++ip){
	Double_t p = histCircle->GetZaxis()->GetBinCenter(ip+1);
	Double_t r = p/(Const*dMagneticField);

	//a*sin(theta) + b*cos(theta) + c = 0
	Double_t a = 2.*(r+rd)*y;
	Double_t b = 2.*(r+rd)*x;
	Double_t c = -1.*(rd*rd + 2.*r*rd + x*x + y*y);

	Double_t r0 = sqrt(a*a + b*b);
	if(fabs(-1.*c/r0)>1.) continue;
	Double_t theta1_alpha = asin(-1.*c/r0);

	Double_t theta2_alpha;
	if(theta1_alpha>0.) theta2_alpha = acos(-1.) - theta1_alpha;
	else theta2_alpha = -1.*acos(-1.) - theta1_alpha;

	Double_t theta_alpha = TMath::ATan2(b, a);

	Double_t xcenter1 = (r+rd)*cos(theta1_alpha - theta_alpha);
	Double_t ycenter1 = (r+rd)*sin(theta1_alpha - theta_alpha);
	Double_t r_re1 = sqrt(pow(x-xcenter1,2) + pow(y-ycenter1,2));

	Double_t xcenter2 = (r+rd)*cos(theta2_alpha - theta_alpha);
	Double_t ycenter2 = (r+rd)*sin(theta2_alpha - theta_alpha);
	Double_t r_re2 = sqrt(pow(x-xcenter1,2) + pow(y-ycenter1,2));

	Double_t theta1 = TMath::ATan2(ycenter1, xcenter1);
	Double_t theta2 = TMath::ATan2(ycenter2, xcenter2);
	if(TMath::IsNaN(theta1)){
	  std::cout<<"theta1="<<theta1<<", x="<<x<<", y="<<y
		   <<"rd="<<rd<<", r"<<r<<std::endl;
	}
	if(fabs(r-r_re1)>0.01||fabs(r-r_re2)>0.01){
	  std::cout<<"r="<<r<<", r_re1="<<r_re1<<", r_re1="<<r_re2<<std::endl;
	  std::cout<<"x:"<<x<<", y:"<<y
		   <<", theta1:"<<theta1<<", theta2:"<<theta2
		   <<", theta_alpha:"<<theta_alpha<<std::endl;
	}
	histCircle->Fill(rd, theta1, p);
	histCircle->Fill(rd, theta2, p);
      } // p bins
    } // rdiff bins
  } //i

#if DebugEvDisp
  histCircle->Draw("colz");
  gPad->Modified();
  gPad->Update();
  c1.Print(Form("c%d.pdf",cannum));
  getchar();
  cannum++;
#endif

  Int_t maxbin = histCircle->GetMaximumBin();
  Int_t mx, my, mz;
  histCircle->GetBinXYZ(maxbin, mx, my, mz);
  MaxBin[0] = mx; MaxBin[1] = my; MaxBin[2] = mz;

  Double_t hough_rd = histCircle->GetXaxis()->GetBinCenter(mx);
  Double_t hough_theta = histCircle->GetYaxis()->GetBinCenter(my);
  Double_t hough_p = histCircle->GetZaxis()->GetBinCenter(mz);
  HelixPar[3] = hough_p/(Const*dMagneticField); //helix r
  HelixPar[0] = (HelixPar[3] + hough_rd)*cos(hough_theta); //helix cx
  HelixPar[1] = (HelixPar[3] + hough_rd)*sin(hough_theta); //helix cy

#if DebugDisp
  std::cout<<FUNC_NAME+" XZ hough vote : "<<histCircle->GetMaximum()<<"/"<<gHitPos.size()<<std::endl;
  std::cout<<FUNC_NAME+" Circle hough maxbin : "<<histCircle->GetBinContent(maxbin)<<std::endl;
  std::cout<<FUNC_NAME+" Hough(rd, theta, p)=("<<hough_rd<<", "
	   <<hough_theta<<", "<<hough_p<<")"<<std::endl;
#endif

  if(histCircle->GetMaximum() < 0.5*MinNumOfHits) status = false;
  return status;
}

//_____________________________________________________________________________
Bool_t
HoughTransformCircleXZ(const std::vector<TPCClusterContainer>& ClCont,
		       Int_t *MaxBin, Double_t *HelixPar,
		       Int_t MinNumOfHits /*=8*/)
{

  std::vector<TVector3> gHitPos;
  HoughFlagCheck(ClCont, gHitPos);
  return HoughTransformCircleXZ(gHitPos, MaxBin, HelixPar, MinNumOfHits);

}

} //namespace tpc

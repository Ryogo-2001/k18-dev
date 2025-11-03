/*
  CFTLocalTrack.cc

  2012/1/24
*/

#include "CFTLocalTrack.hh"
//#include "DCLTrackHit.hh"
#include "DCGeomMan.hh"
#include "DCAnalyzer.hh"
#include "DetectorID.hh"
//#include "EventDisplay.hh"

#include "MathTools.hh"

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <stdexcept>
#include <sstream>

namespace
{
  const auto& gGeom = DCGeomMan::GetInstance();
  const double Deg2Rad = TMath::DegToRad();
  const double Rad2Deg = TMath::RadToDeg();

  const int ReservedNumOfHits  = 4;
  const int ReservedNumOfHitsU  = 4;
  const int DCLocalMinNHits    = 3;

  const double offsetCATCH = 155;
  //EventDisplay&    gEvDisp = EventDisplay::GetInstance();
}


CFTLocalTrack::CFTLocalTrack()
  : status_(false), x0_(0.0), y0_(0.0), u0_(0.0), v0_(0.0),a_(0.0),b_(0.0),
    gftstatus_(true), xyFitFlag_(-1), Axy_(0.), Bxy_(0.), Az_(0.), Bz_(0.),
    zTrackFlag_(-1), TotalDEHiGain_(0.), TotalDELowGain_(0.),
    vtx_z_(-100000.),
    chisqr_dE_(1.e+10),
    chisqrXY_(1.e+10),
    chisqrZ_(1.e+10),
    nth_ztracking_(0),
    MaxDEHiGain_(0.),
    MaxDELowGain_(0.),
    MaxDEHiGain_phi_(0.),
    MaxDELowGain_phi_(0.),
    MaxDEHiGain_uv_(0.),
    MaxDELowGain_uv_(0.),
    TotalDELowGain_phi_(0.),
    TotalDELowGain_uv_(0.),
    pathlength_(0.),
    theta_(-999.),
    OffsetPos0_(2.5, -1.5, 0 )
{
  hitArray.reserve( ReservedNumOfHits );
  hitArrayU.reserve( ReservedNumOfHitsU );
}

CFTLocalTrack::~CFTLocalTrack()
{
}

CFTFiberCluster *CFTLocalTrack::GetHit( std::size_t nth ) const
{
  if( nth<hitArray.size() )
    return hitArray[nth];
  else
    return 0;
}

CFTFiberCluster *CFTLocalTrack::GetHitU( std::size_t nth ) const
{
  if( nth<hitArrayU.size() )
    return hitArrayU[nth];
  else
    return 0;
}

CFTFiberCluster *CFTLocalTrack::GetHitOfLayerNumber( int lnum ) const
{
  for( std::size_t i=0; i<hitArray.size(); ++i ) {
    CFTFiberCluster *hitp = hitArray[i];
    Int_t layer = gGeom.GetDetectorId(hitp->DetectorName()+"-"+hitp->PlaneName());
    if( layer == lnum )
      return hitArray[i];
  }

  for( std::size_t i=0; i<hitArrayU.size(); ++i ) {
    CFTFiberCluster *hitp = hitArrayU[i];
    Int_t layer = gGeom.GetDetectorId(hitp->DetectorName()+"-"+hitp->PlaneName());
    if( layer == lnum )
      return hitArrayU[i];
  }

  return 0;
}


int CFTLocalTrack::GetFirstLayerPhi() const
{
  if (hitArray.size() <= 0)
    return -1;

  int firstLayer = 999;

  for( std::size_t i=0; i<hitArray.size(); ++i ) {
    CFTFiberCluster *hitp = hitArray[i];
    Int_t layer = gGeom.GetDetectorId(hitp->DetectorName()+"-"+hitp->PlaneName());
    if( layer < firstLayer )
      firstLayer = layer;
  }

  return firstLayer;

}

int CFTLocalTrack::GetFirstLayerUV() const
{
  if (hitArrayU.size() <= 0)
    return -1;

  int firstLayer = 999;

  for( std::size_t i=0; i<hitArrayU.size(); ++i ) {
    CFTFiberCluster *hitp = hitArrayU[i];
    Int_t layer = gGeom.GetDetectorId(hitp->DetectorName()+"-"+hitp->PlaneName());
    if( layer < firstLayer )
      firstLayer = layer;
  }

  return firstLayer;

}


bool CFTLocalTrack::DoFitXY( void )
{
  const std::string funcname = "[CFTLocalTrack::DoFit()]";
  //  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  std::size_t n = hitArray.size();

  if(n < DCLocalMinNHits ) return status_ = false;

  std::vector <double> x, y, s;
  x.reserve(n); y.reserve(n); s.reserve(n);

  for( std::size_t i=0; i<n; ++i ){
    CFTFiberCluster *hitp = hitArray[i];
    if( hitp ){
      Int_t lnum = gGeom.GetDetectorId(hitp->DetectorName()+"-"+hitp->PlaneName());

      double xx = hitp->MeanX();
      double yy = hitp->MeanY();
      double ss = gGeom.GetResolution( lnum );


      x.push_back( xx ); y.push_back( yy ); s.push_back(ss);

#if 0
      double phi = hitp->MeanPhi();
      std::cout << std::setw(10) << "layer = " << lnum
		<< std::setw(10) << "x  = " << xx << ", y  =  " << yy
		<< ", phi = " << phi << ", segment " << hitp->MeanSeg()
		<< std::endl;
#endif
    }
  }

  double A=0., B=0., C=0., D=0., E=0., F=0.;

  // y = ax + b
  for (int i=0; i<n; i++) {
    A += x[i]/(s[i]*s[i]);
    B += 1./(s[i]*s[i]);
    C += y[i]/(s[i]*s[i]);
    D += x[i]*x[i]/(s[i]*s[i]);
    E += x[i]*y[i]/(s[i]*s[i]);
    F += y[i]*y[i]/(s[i]*s[i]);
  }

  Axy_=(E*B-C*A)/(D*B-A*A);
  Bxy_=(D*C-E*A)/(D*B-A*A);

  if (fabs(Axy_) <= 1.) {
    //y  = ax + b

    xyFitFlag_=0;

    chisqrXY_ = 0.;
    for (int i=0; i<n; i++) {
      double ycal;
      ycal = Axy_*x[i] + Bxy_;

      chisqrXY_ += (y[i]-ycal)*(y[i]-ycal)/(s[i]*s[i]);

      CFTFiberCluster *hitp = hitArray[i];
      if( hitp ){
	hitp->SetXYFitFlag(xyFitFlag_);
	hitp->SetCalPosition(x[i], ycal);
      }
    }
    chisqrXY_ /= n;

  } else {
    // x = ay + b

    xyFitFlag_=1;

    A=0.; B=0.; C=0.; D=0.; E=0.; F=0.;
    for (int i=0; i<n; i++) {
      A += y[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += x[i]/(s[i]*s[i]);
      D += y[i]*y[i]/(s[i]*s[i]);
      E += x[i]*y[i]/(s[i]*s[i]);
      F += x[i]*x[i]/(s[i]*s[i]);
    }

    Axy_=(E*B-C*A)/(D*B-A*A);
    Bxy_=(D*C-E*A)/(D*B-A*A);

    chisqrXY_ = 0.;
    for (int i=0; i<n; i++) {
      double xcal;
      xcal = Axy_*y[i] + Bxy_;

      chisqrXY_ += (x[i]-xcal)*(x[i]-xcal)/(s[i]*s[i]);

    }
    chisqrXY_ /= n;
  }

  return status_=true;
}

bool CFTLocalTrack::DoFitXY_2nd( void )
{
  const std::string funcname = "[CFTLocalTrack::DoFit()]";
  //const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  std::size_t n = hitArray.size();

  if(n < DCLocalMinNHits ) return status_ = false;

  std::vector <double> x, y, s;
  x.reserve(n); y.reserve(n); s.reserve(n);

  for( std::size_t i=0; i<n; ++i ){
    CFTFiberCluster *hitp = hitArray[i];
    if( hitp ){
      Int_t lnum = gGeom.GetDetectorId(hitp->DetectorName()+"-"+hitp->PlaneName());

      double xx = hitp->MeanXCor();
      double yy = hitp->MeanYCor();
      double ss = gGeom.GetResolution( lnum );


      x.push_back( xx ); y.push_back( yy ); s.push_back(ss);

#if 0
      double phi = hitp->GetPhiCor();
      std::cout << std::setw(10) << "layer = " << lnum
		<< std::setw(10) << "x  = " << xx << ", y  =  " << yy
		<< ", phi = " << phi << ", segment " << hitp->MeanPairId()
		<< std::endl;
#endif
    }
  }

  double A=0., B=0., C=0., D=0., E=0., F=0.;

  // y = ax + b
  for (int i=0; i<n; i++) {
    A += x[i]/(s[i]*s[i]);
    B += 1./(s[i]*s[i]);
    C += y[i]/(s[i]*s[i]);
    D += x[i]*x[i]/(s[i]*s[i]);
    E += x[i]*y[i]/(s[i]*s[i]);
    F += y[i]*y[i]/(s[i]*s[i]);
  }

  Axy_=(E*B-C*A)/(D*B-A*A);
  Bxy_=(D*C-E*A)/(D*B-A*A);

  if (fabs(Axy_) <= 1.) {
    //y  = ax + b

    xyFitFlag_=0;

    chisqrXY_ = 0.;
    for (int i=0; i<n; i++) {
      double ycal;
      ycal = Axy_*x[i] + Bxy_;

      chisqrXY_ += (y[i]-ycal)*(y[i]-ycal)/(s[i]*s[i]);

      CFTFiberCluster *hitp = hitArray[i];
      if( hitp ){
	hitp->SetXYFitFlag(xyFitFlag_);
	hitp->SetCalPosition(x[i], ycal);
      }
    }
    chisqrXY_ /= n;

  } else {
    // x = ay + b

    xyFitFlag_=1;

    A=0.; B=0.; C=0.; D=0.; E=0.; F=0.;
    for (int i=0; i<n; i++) {
      A += y[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += x[i]/(s[i]*s[i]);
      D += y[i]*y[i]/(s[i]*s[i]);
      E += x[i]*y[i]/(s[i]*s[i]);
      F += x[i]*x[i]/(s[i]*s[i]);
    }

    Axy_=(E*B-C*A)/(D*B-A*A);
    Bxy_=(D*C-E*A)/(D*B-A*A);

    chisqrXY_ = 0.;
    for (int i=0; i<n; i++) {
      double xcal;
      xcal = Axy_*y[i] + Bxy_;

      chisqrXY_ += (x[i]-xcal)*(x[i]-xcal)/(s[i]*s[i]);

    }
    chisqrXY_ /= n;
  }

  return status_=true;
}

bool CFTLocalTrack::DoFitZTrack()
{
  const std::string funcname = "[CFTLocalTrack::DoFitZTrack()]";
  //const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  nth_ztracking_++;

  std::size_t n = hitArrayU.size();

  if(n < DCLocalMinNHits ) return status_ = false;

  std::vector <double> z, xy, s;
  z.reserve(n); xy.reserve(n); s.reserve(n);

  int n_hitU=0;
  double mean_z = 0;
  for( std::size_t i=0; i<n; ++i ){
    CFTFiberCluster *hitp = hitArrayU[i];
    if( hitp ){
      Int_t lnum = gGeom.GetDetectorId(hitp->DetectorName()+"-"+hitp->PlaneName());
      double r = gGeom.GetLocalZ(lnum);
      double phi;
      bool ret = GetCrossPointR(r, &phi);
      if (!ret) {
	//std::cout << funcname << " return at GetCrossPointR" << std::endl;
	return false;
      }

      //double z1 = CalculateZpos2(phi, hitp);
      double z1 = CalculateZpos(phi, hitp);
      mean_z += z1;
      n_hitU++;
    }
  }
  mean_z /= n_hitU;

  for( std::size_t i=0; i<n; ++i ){
    CFTFiberCluster *hitp = hitArrayU[i];
    if( hitp ){
      Int_t lnum = gGeom.GetDetectorId(hitp->DetectorName()+"-"+hitp->PlaneName());
      double r = gGeom.GetLocalZ(lnum);
      double phi;
      bool ret = GetCrossPointR(r, &phi);
      if (!ret) {
	//std::cout << funcname << " return at GetCrossPointR" << std::endl;
	return false;
      }

      // modosu 2->1
      double z1 = CalculateZpos2(phi, hitp, mean_z);
      // double z1 = CalculateZpos(phi, hitp, mean_z);

      double xy1=-999.;
      if (xyFitFlag_==0)
	xy1 = r*cos(phi*Deg2Rad);
      else
	xy1 = r*sin(phi*Deg2Rad);

      double ss = gGeom.GetResolution( lnum );

      z.push_back( z1 ); xy.push_back( xy1 ); s.push_back(ss);

#if 0
      std::cout << std::setw(10) << "layer = " << lnum
		<< std::setw(10) << "z  = " << z1 + offsetCATCH << ", xy  =  " << xy1
		<< ", phi = " << phi << ", Z0 = " << hitp->GetZ0()
		<< ", segment = " << hitp->MeanPairId()
		<< std::endl;
#endif
    }
  }

  double slope = (xy[n-1]-xy[0])/(z[n-1]-z[0]);
  /*
  if (fabs(slope)<=1) {
    //if (0) {

    double A=0., B=0., C=0., D=0., E=0., F=0.;

    // x = az + b
    for (int i=0; i<n; i++) {
      A += z[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += xy[i]/(s[i]*s[i]);
      D += z[i]*z[i]/(s[i]*s[i]);
      E += z[i]*xy[i]/(s[i]*s[i]);
      F += xy[i]*xy[i]/(s[i]*s[i]);
    }

    Az_=(E*B-C*A)/(D*B-A*A);
    Bz_=(D*C-E*A)/(D*B-A*A);

    if (xyFitFlag_==0) {
      u0_ = Az_;
      x0_ = Bz_;
      v0_ = Axy_*Az_;
      y0_ = Axy_*Bz_+Bxy_;
    } else {
      u0_ = Axy_*Az_;
      x0_ = Axy_*Bz_+Bxy_;
      v0_ = Az_;
      y0_ = Bz_;
    }

    chisqrZ_ = 0.;
    for (int i=0; i<n; i++) {
      double xycal = Az_*z[i] + Bz_;
      double xcal = u0_*z[i] + x0_;
      double ycal = v0_*z[i] + y0_;

      chisqrZ_ += (xy[i]-xycal)*(xy[i]-xycal)/(s[i]*s[i]);

      double xy2 = xy[i];
      double zcal = (xy2-Bz_)/Az_;
      //chisqrZ_ += (z[i]-zcal)*(z[i]-zcal)/(s[i]*s[i]);

    }
    //chisqrZ_ /= 2*n;
    chisqrZ_ /= n;

    vtx_z_ = -Bz_/Az_;
    //std::cout << "chisqrZ : " << chisqrZ_ << std::endl;

    if (chisqrZ_ < 10) {
      std::cout << "chisqrZ : " << chisqrZ_ << std::endl;
      for (int i=0; i<n; i++) {
	double xycal = Az_*z[i] + Bz_;
	double xcal = u0_*z[i] + x0_;
	double ycal = v0_*z[i] + y0_;
	std::cout << n << " : z = " << z[i];
	if (xyFitFlag_==0) {
	  std::cout << " , x = " << xy[i] << " (" << xy[i]-xycal << " )" << std::endl;
	} else
	  std::cout << " , y = " << xy[i] << " (" << xy[i]-xycal << " )" << std::endl;
      }
    }

    zTrackFlag_ = 0;
  }   */
  //else {
  if (1) {
    // z = ax + b or z = ay + b
    double A=0., B=0., C=0., D=0., E=0., F=0.;
    for (int i=0; i<n; i++) {
      A += xy[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += z[i]/(s[i]*s[i]);
      D += xy[i]*xy[i]/(s[i]*s[i]);
      E += z[i]*xy[i]/(s[i]*s[i]);
      F += z[i]*z[i]/(s[i]*s[i]);
    }

    Az_=(E*B-C*A)/(D*B-A*A);
    Bz_=(D*C-E*A)/(D*B-A*A);

    if (xyFitFlag_==0) {
      u0_ = 1./Az_;
      x0_ = -Bz_/Az_;
      v0_ = Axy_/Az_;
      y0_ = -Axy_*Bz_/Az_+Bxy_;
    } else {
      u0_ = Axy_/Az_;
      x0_ = -Axy_*Bz_/Az_+Bxy_;
      v0_ = 1./Az_;
      y0_ = -Bz_/Az_;
    }

    chisqrZ_ = 0.;
    for (int i=0; i<n; i++) {
      double zcal = Az_*xy[i] + Bz_;

      chisqrZ_ += (z[i]-zcal)*(z[i]-zcal)/(s[i]*s[i]);

      double z2 = z[i];
      double xycal = 1./Az_*z2-Bz_/Az_;
      //chisqrZ_ += (xy[i]-xycal)*(xy[i]-xycal)/(s[i]*s[i]);

    }
    chisqrZ_ /= n;
    //chisqrZ_ /= 2*n;

    vtx_z_ = Bz_;

    //std::cout << "chisqrZ : " << chisqrZ_ << std::endl;
    /*
    if (chisqrZ_ < 100 && n==4) {
      std::cout << "chisqrZ : " << chisqrZ_ << std::endl;
      for (int i=0; i<n; i++) {
	double zcal = Az_*xy[i] + Bz_;

	if (xyFitFlag_==0) {
	  std::cout << "x :" << xy[i] << " , z = " << z[i] << " (" << z[i]-zcal << " )" << std::endl;
	} else
	  std::cout << "y : " << xy[i] << " , z = " << z[i] << " (" << z[i]-zcal << " )" << std::endl;
      }
    }
    */

    zTrackFlag_ = 1;

  }
#if 0
  if(nth_ztracking_ == 2 && gEvDisp.IsReady() && chisqrZ_<20){
    bool flagZX = false;
    if (xyFitFlag_==0)
      flagZX = true;

    gEvDisp.DrawCFTTrackCand(xy, z, flagZX, Axy_ );
  }
#endif
  return true;
}

bool CFTLocalTrack::DoFitZTrack_16layer()
{
  const std::string funcname = "[CFTLocalTrack::DoFitZTrack_16layer()]";
  //const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  std::size_t n = hitArrayU.size();

  if(n < DCLocalMinNHits ) return status_ = false;

  std::vector <double> z, xy, s;
  z.reserve(n); xy.reserve(n); s.reserve(n);

  for( std::size_t i=0; i<n; ++i ){
    CFTFiberCluster *hitp = hitArrayU[i];
    if( hitp ){
      Int_t lnum = gGeom.GetDetectorId(hitp->DetectorName()+"-"+hitp->PlaneName());
      double r = gGeom.GetLocalZ(lnum);
      double phi, phi1, phi2;

      bool ret = GetCrossPointR(r, &phi1, &phi2);
      if (!ret) {
	//std::cout << funcname << " return at GetCrossPointR" << std::endl;
	return false;
      }

      double phiCal = hitp->GetCalPhi();
      if (fabs(phi1-phiCal) <= fabs(phi2-phiCal) && fabs(phi1-phiCal) < 10.)
	phi = phi1;
      else if (fabs(phi1-phiCal) >= fabs(phi2-phiCal) && fabs(phi2-phiCal) < 10.)
	phi = phi2;
      else if (fabs(phi1+phiCal-360.)<=10. && (phi1<=15 || phi1 >=345))
	phi = phi1;
      else
	phi = phi2;
      /*
      if (fabs(phi1-phiCal) <= fabs(phi2-phiCal))
	phi = phi1;
      else
	phi = phi2;
      */

      //double z1 = CalculateZpos2(phi, hitp);
      double z1 = hitp->GetZIni();
      //double z1 = CalculateZpos(phi, hitp);

      double xy1=-999.;
      if (xyFitFlag_==0)
	xy1 = r*cos(phi*Deg2Rad);
      else
	xy1 = r*sin(phi*Deg2Rad);

      double ss = gGeom.GetResolution( lnum );

      z.push_back( z1 ); xy.push_back( xy1 ); s.push_back(ss);

#if 0
      std::cout << std::setw(10) << "layer = " << lnum
		<< std::setw(10) << "z  = " << z1 + offsetCATCH << ", xy  =  " << xy1
		<< ", phi = " << phi << ", Z0 = " << hitp->GetZ0()
		<< ", segment = " << hitp->MeanPairId()
		<< std::endl;
#endif
    }
  }

  double slope = (xy[n-1]-xy[0])/(z[n-1]-z[0]);
  /*
  if (fabs(slope)<=1) {

    double A=0., B=0., C=0., D=0., E=0., F=0.;

    // x = az + b
    for (int i=0; i<n; i++) {
      A += z[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += xy[i]/(s[i]*s[i]);
      D += z[i]*z[i]/(s[i]*s[i]);
      E += z[i]*xy[i]/(s[i]*s[i]);
      F += xy[i]*xy[i]/(s[i]*s[i]);
    }

    Az_=(E*B-C*A)/(D*B-A*A);
    Bz_=(D*C-E*A)/(D*B-A*A);

    if (xyFitFlag_==0) {
      u0_ = Az_;
      x0_ = Bz_;
      v0_ = Axy_*Az_;
      y0_ = Axy_*Bz_+Bxy_;
    } else {
      u0_ = Axy_*Az_;
      x0_ = Axy_*Bz_+Bxy_;
      v0_ = Az_;
      y0_ = Bz_;
    }

    chisqrZ_ = 0.;
    for (int i=0; i<n; i++) {
      double xycal = Az_*z[i] + Bz_;
      double xcal = u0_*z[i] + x0_;
      double ycal = v0_*z[i] + y0_;

      chisqrZ_ += (xy[i]-xycal)*(xy[i]-xycal)/(s[i]*s[i]);
    }
    chisqrZ_ /= n;

    //std::cout << "chisqrZ : " << chisqrZ_ << std::endl;

    zTrackFlag_ = 0;
  } else {
  */
  if (1) {
    // z = ax + b or z = ay + b
    double A=0., B=0., C=0., D=0., E=0., F=0.;
    for (int i=0; i<n; i++) {
      A += xy[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += z[i]/(s[i]*s[i]);
      D += xy[i]*xy[i]/(s[i]*s[i]);
      E += z[i]*xy[i]/(s[i]*s[i]);
      F += z[i]*z[i]/(s[i]*s[i]);
    }

    Az_=(E*B-C*A)/(D*B-A*A);
    Bz_=(D*C-E*A)/(D*B-A*A);

    if (xyFitFlag_==0) {
      u0_ = 1./Az_;
      x0_ = -Bz_/Az_;
      v0_ = Axy_/Az_;
      y0_ = -Axy_*Bz_/Az_+Bxy_;
    } else {
      u0_ = Axy_/Az_;
      x0_ = -Axy_*Bz_/Az_+Bxy_;
      v0_ = 1./Az_;
      y0_ = -Bz_/Az_;
    }

    chisqrZ_ = 0.;
    for (int i=0; i<n; i++) {
      double zcal = Az_*xy[i] + Bz_;

      chisqrZ_ += (z[i]-zcal)*(z[i]-zcal)/(s[i]*s[i]);
    }
    chisqrZ_ /= n;

    //std::cout << "chisqrZ : " << chisqrZ_ << std::endl;

    zTrackFlag_ = 1;

  }

  return true;
}


bool CFTLocalTrack::DoFitXY_wVtx( void )
{
  const std::string funcname = "[CFTLocalTrack::DoFitXY_wVtx()]";
  //const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  std::size_t n = hitArray.size() + 1; // +1 is Vtx

  if(n < DCLocalMinNHits ) return status_ = false;

  std::vector <double> x, y, s;
  x.reserve(n); y.reserve(n); s.reserve(n);

  x.push_back( CFTVtx_.x() ); y.push_back( CFTVtx_.x() ); s.push_back(5.);

  for( std::size_t i=0; i<n-1; ++i ){
    CFTFiberCluster *hitp = hitArray[i];
    if( hitp ){
      Int_t lnum = gGeom.GetDetectorId(hitp->DetectorName()+"-"+hitp->PlaneName());

      double xx = hitp->MeanX();
      double yy = hitp->MeanY();
      double ss = gGeom.GetResolution( lnum );


      x.push_back( xx ); y.push_back( yy ); s.push_back(ss);

#if 0
      double phi = hitp->GetPhi();
      std::cout << std::setw(10) << "layer = " << lnum
		<< std::setw(10) << "x  = " << xx << ", y  =  " << yy
		<< ", phi = " << phi
		<< std::endl;
#endif
    }
  }

  double A=0., B=0., C=0., D=0., E=0., F=0.;

  // y = ax + b
  for (int i=0; i<n; i++) {
    A += x[i]/(s[i]*s[i]);
    B += 1./(s[i]*s[i]);
    C += y[i]/(s[i]*s[i]);
    D += x[i]*x[i]/(s[i]*s[i]);
    E += x[i]*y[i]/(s[i]*s[i]);
    F += y[i]*y[i]/(s[i]*s[i]);
  }

  Axy_=(E*B-C*A)/(D*B-A*A);
  Bxy_=(D*C-E*A)/(D*B-A*A);

  if (fabs(Axy_) <= 1.) {
    //y  = ax + b

    xyFitFlag_=0;

    chisqrXY_ = 0.;
    for (int i=0; i<n; i++) {
      double ycal;
      ycal = Axy_*x[i] + Bxy_;

      chisqrXY_ += (y[i]-ycal)*(y[i]-ycal)/(s[i]*s[i]);

      if (i>0) {
	CFTFiberCluster *hitp = hitArray[i-1];
	if( hitp ){
	  hitp->SetXYFitFlag(xyFitFlag_);
	  hitp->SetCalPosition(x[i], ycal);
	}
      }
    }
    chisqrXY_ /= n;

  } else {
    // x = ay + b

    xyFitFlag_=1;

    A=0.; B=0.; C=0.; D=0.; E=0.; F=0.;
    for (int i=0; i<n; i++) {
      A += y[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += x[i]/(s[i]*s[i]);
      D += y[i]*y[i]/(s[i]*s[i]);
      E += x[i]*y[i]/(s[i]*s[i]);
      F += x[i]*x[i]/(s[i]*s[i]);
    }

    Axy_=(E*B-C*A)/(D*B-A*A);
    Bxy_=(D*C-E*A)/(D*B-A*A);

    chisqrXY_ = 0.;
    for (int i=0; i<n; i++) {
      double xcal;
      xcal = Axy_*y[i] + Bxy_;

      chisqrXY_ += (x[i]-xcal)*(x[i]-xcal)/(s[i]*s[i]);

    }
    chisqrXY_ /= n;
  }

  return status_=true;
}

bool CFTLocalTrack::CheckPhi( void )
{
  const std::string funcname = "[CFTLocalTrack::CheckPhi()]";
  //const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  std::size_t n = hitArray.size();

  if(n < 2 ) return false;

  double DeltaPhi = 5.;

  double phi0 = -999.;
  for( std::size_t i=0; i<n; ++i ){
    CFTFiberCluster *hitp = hitArray[i];
    if( hitp ){
      if (i==0)
	//phi0 = hitp->GetPhi();
	phi0 = hitp->MeanPhi();
      else {
	//double phi = hitp->GetPhi();
	double phi = hitp->MeanPhi();
	if (std::abs(phi-phi0)<DeltaPhi || std::abs(phi-phi0-360.)<DeltaPhi)
	  ;
	else
	  return false;
      }
    }
  }

  return true;
}

bool CFTLocalTrack::CheckPhi_1st( void )
{
  const std::string funcname = "[CFTLocalTrack::CheckPhi_1st()]";
  //const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  std::size_t n = hitArray.size();

  if(n < 2 ) return false;

  double DeltaPhi = 30.;

  double phi0 = -999.;
  for( std::size_t i=0; i<n; ++i ){
    CFTFiberCluster *hitp = hitArray[i];
    if( hitp ){
      if (i==0)
	//phi0 = hitp->GetPhi();
	phi0 = hitp->MeanPhi();
      else {
	//double phi = hitp->GetPhi();
	double phi = hitp->MeanPhi();
	if (std::abs(phi-phi0)<DeltaPhi || std::abs(phi-phi0-360.)<DeltaPhi)
	  ;
	else
	  return false;
      }
    }
  }

  return true;
}

bool CFTLocalTrack::CheckEdep( void )
{
  const std::string funcname = "[CFTLocalTrack::CheckPhi()]";
  //const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  std::size_t n = hitArray.size();

  if(n < 2 ) return false;

  bool flagLowDE = true;
  bool flagHighDE = true;

  for( std::size_t i=0; i<n; ++i ){
    CFTFiberCluster *hitp = hitArray[i];
    if (hitp) {
      double de = hitp->TotalDeltaE_Hi();

      if (de>1.5)
	flagLowDE = false;

      if (de<0.5)
	flagHighDE = false;
    }
  }


  std::size_t nU = hitArrayU.size();

  if(nU < 2 ) return false;

  for( std::size_t i=0; i<nU; ++i ){
    CFTFiberCluster *hitp = hitArrayU[i];
    if( hitp ){
      double de = hitp->TotalDeltaE_Hi();

      if (de>1.5)
	flagLowDE = false;

      if (de<1.0)
	flagHighDE = false;
    }
  }

  if (flagHighDE || flagLowDE)
    return true;
  else
    return false;
}


bool CFTLocalTrack::DoFitZTrack_wVtx()
{
  const std::string funcname = "[CFTLocalTrack:DoF:itZTrack_vVtx()]";
  //const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  std::size_t n = 0;
  std::size_t nHit = hitArrayU.size();

  if (nHit>=3)
    n=nHit;
  else
    n=nHit + 1;

  if(n < DCLocalMinNHits ) return status_ = false;

  std::vector <double> z, xy, s;
  z.reserve(n); xy.reserve(n); s.reserve(n);

  if (hitArrayU.size()<=2) {
    if (xyFitFlag_==0 ) {
      z.push_back( CFTVtx_.z() ); xy.push_back( CFTVtx_.x() ); s.push_back(5);
    } else {
      z.push_back( CFTVtx_.z() ); xy.push_back( CFTVtx_.y() ); s.push_back(10);
    }

  }

  for( std::size_t i=0; i<nHit; ++i ){
    CFTFiberCluster *hitp = hitArrayU[i];
    if( hitp ){
      Int_t lnum = gGeom.GetDetectorId(hitp->DetectorName()+"-"+hitp->PlaneName());
      double r = gGeom.GetLocalZ(lnum);
      double phi;
      bool ret = GetCrossPointR(r, &phi);
      if (!ret) {
	//std::cout << funcname << " return at GetCrossPointR" << std::endl;
	return false;
      }

      double z1 = CalculateZpos(phi, hitp);

      double xy1=-999.;
      if (xyFitFlag_==0)
	xy1 = r*cos(phi*Deg2Rad);
      else
	xy1 = r*sin(phi*Deg2Rad);

      double ss = gGeom.GetResolution( lnum );

      z.push_back( z1 ); xy.push_back( xy1 ); s.push_back(ss);

#if 0
      std::cout << std::setw(10) << "layer = " << lnum
		<< std::setw(10) << "z  = " << z1 << ", xy  =  " << xy1
		<< ", phi = " << phi
		<< std::endl;
#endif
    }
  }

  double slope = (xy[n-1]-xy[0])/(z[n-1]-z[0]);

  if (std::abs(slope)<=1) {

    double A=0., B=0., C=0., D=0., E=0., F=0.;

    // x = az + b
    for (int i=0; i<n; i++) {
      A += z[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += xy[i]/(s[i]*s[i]);
      D += z[i]*z[i]/(s[i]*s[i]);
      E += z[i]*xy[i]/(s[i]*s[i]);
      F += xy[i]*xy[i]/(s[i]*s[i]);
    }

    Az_=(E*B-C*A)/(D*B-A*A);
    Bz_=(D*C-E*A)/(D*B-A*A);

    if (xyFitFlag_==0) {
      u0_ = Az_;
      x0_ = Bz_;
      v0_ = Axy_*Az_;
      y0_ = Axy_*Bz_+Bxy_;
    } else {
      u0_ = Axy_*Az_;
      x0_ = Axy_*Bz_+Bxy_;
      v0_ = Az_;
      y0_ = Bz_;
    }

    chisqrZ_ = 0.;
    for (int i=0; i<n; i++) {
      double xycal = Az_*z[i] + Bz_;
      double xcal = u0_*z[i] + x0_;
      double ycal = v0_*z[i] + y0_;

      chisqrZ_ += (xy[i]-xycal)*(xy[i]-xycal)/(s[i]*s[i]);
    }
    chisqrZ_ /= n;

    //std::cout << "chisqrZ : " << chisqrZ_ << std::endl;

    zTrackFlag_ = 0;
  } else {

    // z = ax + b or z = ay + b
    double A=0., B=0., C=0., D=0., E=0., F=0.;
    for (int i=0; i<n; i++) {
      A += xy[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += z[i]/(s[i]*s[i]);
      D += xy[i]*xy[i]/(s[i]*s[i]);
      E += z[i]*xy[i]/(s[i]*s[i]);
      F += z[i]*z[i]/(s[i]*s[i]);
    }

    Az_=(E*B-C*A)/(D*B-A*A);
    Bz_=(D*C-E*A)/(D*B-A*A);

    if (xyFitFlag_==0) {
      u0_ = 1./Az_;
      x0_ = -Bz_/Az_;
      v0_ = Axy_/Az_;
      y0_ = -Axy_*Bz_/Az_+Bxy_;
    } else {
      u0_ = Axy_/Az_;
      x0_ = -Axy_*Bz_/Az_+Bxy_;
      v0_ = 1./Az_;
      y0_ = -Bz_/Az_;
    }

    chisqrZ_ = 0.;
    for (int i=0; i<n; i++) {
      double zcal = Az_*xy[i] + Bz_;

      chisqrZ_ += (z[i]-zcal)*(z[i]-zcal)/(s[i]*s[i]);
    }
    chisqrZ_ /= n;

    //std::cout << "chisqrZ : " << chisqrZ_ << std::endl;

    zTrackFlag_ = 1;

  }

  return true;
}

#if 0

bool CFTLocalTrack::DoFit( void )
{
  const std::string funcname = "[DCLocalTrack::DoFit()]";
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  std::size_t nPhi = hitArray.size();
  std::size_t nU   = hitArrayU.size();
  std::size_t n = nPhi + nU;

  if(n < DCLocalMinNHits ) return status_ = false;

  std::vector <double> z, w, s, ct, st;
  z.reserve(n); w.reserve(n); s.reserve(n);
  ct.reserve(n); st.reserve(n);

  for( std::size_t i=0; i<nPhi; ++i ){
    CFTFiberCluster *hitp = hitArray[i];
    if( hitp ){
      int lnum = hitp->GetTrackingLayer();
      double ww = geomMan.GetResolution( lnum );
      double zz = hitp->GetR();
      double aa = 0.;

      z.push_back( zz ); w.push_back( 1./(ww*ww) );
      s.push_back( hitp->GetPhi() );
      ct.push_back( cos(aa) ); st.push_back( sin(aa) );

#if 0
      std::cout << std::setw(10) << "layer = " << lnum
		<< std::setw(10) << "wire  = " << hitp->GetWire() << " "
		<< std::setw(20) << "WirePosition = "<<hitp->GetWirePosition() << " "
		<< std::setw(20) << "DriftLength = "<<hitp->GetDriftLength() << " "
		<< std::setw(20) << "hit position = "<<hitp->GetLocalHitPos()<< " "
		<< std::endl;
#endif
    }
  }

  for( std::size_t i=0; i<nU; ++i ){
    CFTFiberCluster *hitp = hitArrayU[i];
    if( hitp ){
      int lnum = hitp->GetTrackingLayer();
      double ww = geomMan.GetResolution( lnum );
      double zz = hitp->GetR();
      double aa = 90.+atan(hitp->GetSlopeU())*Deg2Rad;

      double r = geomMan.GetLocalZ(lnum);
      double phi;
      bool ret = GetCrossPointR(r, &phi);
      if (!ret) {
	//std::cout << funcname << " return at GetCrossPointR" << std::endl;
	return false;
      }

      double uv = CalculateUVpos(phi, hitp);

      z.push_back( zz ); w.push_back( 1./(ww*ww) );
      s.push_back( uv );
      ct.push_back( cos(aa) ); st.push_back( sin(aa) );

#if 0
      std::cout << std::setw(10) << "layer = " << lnum
		<< std::setw(10) << "wire  = " << hitp->GetWire() << " "
		<< std::setw(20) << "WirePosition = "<<hitp->GetWirePosition() << " "
		<< std::setw(20) << "DriftLength = "<<hitp->GetDriftLength() << " "
		<< std::setw(20) << "hit position = "<<hitp->GetLocalHitPos()<< " "
		<< std::endl;
#endif
    }
  }

  std::size_t nn = z.size();
  double matrx[16], *mtp[4], fitp[4];
  mtp[0]=&matrx[0]; mtp[1]=&matrx[4]; mtp[2]=&matrx[8]; mtp[3]=&matrx[12];

  for( int i=0; i<4; ++i ){
    fitp[i]=0.0;
    for( int j=0; j<4; ++j ){
      mtp[i][j]=0.0;
    }
  }

  for( std::size_t i=0; i<nn; ++i ){
    double ww=w[i], zz=z[i], ss=s[i], ctt=ct[i], stt=st[i];
    mtp[0][0] += ww*ctt*ctt;
    mtp[0][1] += ww*zz*ctt*ctt;
    mtp[0][2] += ww*ctt*stt;
    mtp[0][3] += ww*zz*ctt*stt;
    mtp[1][1] += ww*zz*zz*ctt*ctt;
    mtp[1][2] += ww*zz*ctt*stt;
    mtp[1][3] += ww*zz*zz*ctt*stt;
    mtp[2][2] += ww*stt*stt;
    mtp[2][3] += ww*zz*stt*stt;
    mtp[3][3] += ww*zz*zz*stt*stt;

    fitp[0] += ww*ss*ctt;
    fitp[1] += ww*zz*ss*ctt;
    fitp[2] += ww*ss*stt;
    fitp[3] += ww*zz*ss*stt;
  }
  mtp[1][0]=mtp[0][1]; mtp[2][0]=mtp[0][2]; mtp[3][0]=mtp[0][3];
  mtp[2][1]=mtp[1][2]; mtp[3][1]=mtp[1][3]; mtp[3][2]=mtp[2][3];

  std::vector<int> indxc(nn), indxd(nn), ipiv(nn);

#if 0
  std::cout<<"             Vector:  Q_i               "<<std::endl;
  std::cout<<"   A1=  "<<fitp[0]<<"     A2=  "<<fitp[1]<<"     A3=  "<<fitp[2]
	   <<"     A4=  "<<fitp[3]<<std::endl;

  std::cout<<"             original matrix: M_ij        "<<std::endl;
  std::cout<<"    A11="<<mtp[0][0]<<"    A12="<<mtp[0][1]<<"    A13="<<mtp[0][2]
      	   <<"    A14="<<mtp[0][3]<<std::endl;
  std::cout<<"    A21="<<mtp[1][0]<<"    A22="<<mtp[1][1]<<"    A23="<<mtp[1][2]
	   <<"    A24="<<mtp[1][3]<<std::endl;
  std::cout<<"    A31="<<mtp[2][0]<<"    A32="<<mtp[2][1]<<"    A33="<<mtp[2][2]
	   <<"    A34="<<mtp[2][3]<<std::endl;
  std::cout<<"    A41="<<mtp[3][0]<<"    A42="<<mtp[3][1]<<"    A43="<<mtp[3][2]
	   <<"    A44="<<mtp[3][3]<<std::endl;
#endif

  double Org[4][4]={0},Red[4][4]={0},Final[4][4]={0};
  double Org_vec[4]={0}, Solution_vec[4]={0};
  for(int l=0; l<4;l++){
    for(int m=0; m<4; m++){
      Org[l][m]=mtp[l][m];
    }
  }

  for (int i=0; i<4; i++){
    Org_vec[i]=fitp[i];
  }

  if( MathTools::GaussJordan(mtp,4,fitp,&indxc[0],
			     &indxd[0],&ipiv[0])==false ){
    std::cerr << funcname << ": Fitting fails" << std::endl;
    return status_=false;
  }
  f0_=fitp[0]; z0_=fitp[2]; ur0_=fitp[1]; vr0_=fitp[3];

#if 0
  std::cout<<"             reduced matrix        "<<std::endl;
  std::cout<<"     A11="<<mtp[0][0]<<"    A12="<<mtp[0][1]<<"    A13="<<mtp[0][2]
	   <<"     A14="<<mtp[0][3]<<std::endl;
  std::cout<<"     A21="<<mtp[1][0]<<"    A22="<<mtp[1][1]<<"    A23="<<mtp[1][2]
	   <<"     A24="<<mtp[1][3]<<std::endl;
  std::cout<<"     A31="<<mtp[2][0]<<"    A32="<<mtp[2][1]<<"    A33="<<mtp[2][2]
	   <<"     A34="<<mtp[2][3]<<std::endl;
  std::cout<<"     A41="<<mtp[3][0]<<"    A42="<<mtp[3][1]<<"    A43="<<mtp[3][2]
	   <<"     A44="<<mtp[3][3]<<std::endl;
#endif

#if 0
  std::cout<<"             Solution from Gauss-Jordan             "<<std::endl;
  std::cout<<"     A1="<<fitp[0]<<"    A2="<<fitp[1]<<"    A3="<<fitp[2]
	   <<"     A4="<<fitp[3]<<std::endl;
#endif


  for(int l=0; l<4;l++){
    for(int m=0; m<4; m++){
      Red[l][m]=mtp[l][m];
    }
  }

  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      Solution_vec[i] +=Red[i][j]*Org_vec[j];
    }
  }

#if 0
  std::cout<<"             Solution from explicite calculation       "<<std::endl;
  std::cout<<"    A1="<<Solution_vec[0]<<"    A2="<<Solution_vec[1]<<"    A3="<<Solution_vec[2]
	   <<"    A4="<<Solution_vec[3]<<std::endl;
#endif

  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      for(int k=0; k<4; k++){
	Final[i][j] += Red[i][k]*Org[k][j];
      }
    }
  }
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      if(Final[i][j]<1.0e-10) Final[i][j]=0.0;
    }
  }

#if 0
  std::cout<< ""<<std::endl;
  std::cout<<"            final matrix        "<<std::endl;
  std::cout<<"       A11="<<std::setw(10)<<Final[0][0]<<"       A12="<<std::setw(10)<<Final[0][1]
	   <<"       A13="<<std::setw(10)<<Final[0][2]<<"       A14="<<std::setw(10)<<Final[0][3]
	   <<std::endl;
  std::cout<<"       A21="<<std::setw(10)<<Final[1][0]<<"       A22="<<std::setw(10)<<Final[1][1]
	   <<"       A23="<<std::setw(10)<<Final[1][2]<<"       A24="<<std::setw(10)<<Final[1][3]
	   <<std::endl;
  std::cout<<"       A31="<<std::setw(10)<<Final[2][0]<<"       A32="<<std::setw(10)<<Final[2][1]
	   <<"       A33="<<std::setw(10)<<Final[2][2]<<"       A34="<<std::setw(10)<<Final[2][3]
	   <<std::endl;
  std::cout<<"       A41="<<std::setw(10)<<Final[3][0]<<"       A42="<<std::setw(10)<<Final[3][1]
	   <<"       A43="<<std::setw(10)<<Final[3][2]<<"       A44="<<std::setw(10)<<Final[3][3]
	   <<std::endl;

#endif

  double chisqr=0.0;
  for( std::size_t i=0; i<nn; ++i ){
    double ww=w[i], rr=z[i];
    double scal=GetPhiFromR(rr)*ct[i]+GetZFromR(rr)*st[i];
    chisqr += ww*(s[i]-scal)*(s[i]-scal);

#if 0
    if(1){
      //      std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
      //     std::cout<<"x coordinate = "<<(x0_+u0_*zz)<<std::endl;
      //     std::cout<<"y coordinate = "<<(y0_+v0_*zz)<<std::endl;
      std::cout<<std::setw(10)<<"layer = "<<i<<
	std::setw(10)<<"scal = "<<scal<<
	std::setw(10)<<"sdata = "<<s[i]<<std::endl;
      std::cout<<std::setw(10)<<"Res = "<<s[i]-scal<<std::endl;
      std::cout<<std::setw(10)<<"X = "<< GetX(zz)<<" Y = "<< GetY(zz)<<std::endl;
      std::cout<<std::setw(10)<<"chisqr = "<<chisqr<<std::endl;
    }
#endif
  }

  chisqr /= nn-4.;

  /*
  //if(chisqr<2){
  if(1){
    std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
    std::cout << "chisqr = " << chisqr << " nn-4 = " << nn-4 << std::endl;
    std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
  }
  */
  chisqr_=chisqr;

  for( std::size_t i=0; i<nPhi; ++i ){
    CFTFiberCluster *hitp = hitArray[i];
    if( hitp ){
      double r = hitp->GetR();
      /*
	  if(chisqr<2){
	  std::cout<<std::setw(10)<<"lnum = "<< lnum <<std::endl;
	  std::cout<<std::setw(10)<<"X = "<< GetX(zz)<<" Y = "<< GetY(zz)<<std::endl;
	  }
      */
      //hitp->SetCalPosition( GetX(zz), GetY(zz) );
      hitp->SetCalPhi(GetPhiFromR(r));
      hitp->SetCalZ(GetZFromR(r));
    }
  }

  // std::cout << "***********************************************************" << std::endl;

  return status_=true;
}

#endif

bool CFTLocalTrack::SetCalculatedValue()
{
  const std::string funcname = "[CFTLocalTrack::SetCalculatedValue()]";
  //const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  std::size_t n = hitArray.size();

  ThreeVector pos1, pos2;

  for( std::size_t i=0; i<n; ++i ){
    CFTFiberCluster *hitp = hitArray[i];
    if( hitp ){
      Int_t lnum = gGeom.GetDetectorId(hitp->DetectorName()+"-"+hitp->PlaneName());
      double r = gGeom.GetLocalZ(lnum);
      double phi;
      bool ret = GetCrossPointR(r, &phi);
      if (!ret) {
	//std::cout << funcname << " return at GetCrossPointR" << std::endl;
	return false;
      }

      double xcal = r*cos(phi*Deg2Rad);
      double ycal = r*sin(phi*Deg2Rad);
      double zcal;
      if (zTrackFlag_==0) {
	if (xyFitFlag_==0) {
	  zcal = (xcal- Bz_)/Az_;
	} else {
	  zcal = (ycal- Bz_)/Az_;
	}

	if (i==0) {
	  pos1 = ThreeVector(u0_*zcal+x0_, v0_*zcal+y0_, zcal);
	} else if (i==n-1) {
	  pos2 = ThreeVector(u0_*zcal+x0_, v0_*zcal+y0_, zcal);
	}

      } else if (zTrackFlag_==1){
	if (xyFitFlag_==0) {
	  zcal = Az_*xcal + Bz_;
	  ycal = Axy_*xcal + Bxy_;
	} else {
	  zcal = Az_*ycal + Bz_;
	  xcal = Axy_*ycal + Bxy_;
	}
	if (i==0) {
	  pos1 = ThreeVector(xcal, ycal, zcal);
	} else if (i==n-1) {
	  pos2 = ThreeVector(xcal, ycal, zcal);
	}
      }

      hitp->SetCalPhi(phi);
      hitp->SetCalZ(zcal);
      hitp->SetXYFitFlag(xyFitFlag_);
      hitp->SetCalPosition(xcal, ycal);
    }
  }

  Pos0_ = pos1;
  Dir_ = pos2-pos1;

  std::size_t nU = hitArrayU.size();

  for( std::size_t i=0; i<nU; ++i ){
    CFTFiberCluster *hitp = hitArrayU[i];
    if( hitp ){
      Int_t lnum = gGeom.GetDetectorId(hitp->DetectorName()+"-"+hitp->PlaneName());
      double r = gGeom.GetLocalZ(lnum);
      double phi;
      bool ret = GetCrossPointR(r, &phi);
      if (!ret) {
	//std::cout << funcname << " return at GetCrossPointR" << std::endl;
	return false;
      }

      // modosu 2->1
      double z1 = CalculateZpos2(phi, hitp);
      // double z1 = CalculateZpos(phi, hitp);
      double xycal = Az_*z1 + Bz_;
      double xcal = u0_*z1 + x0_;
      double ycal = v0_*z1 + y0_;

      hitp->SetZ(z1);
      hitp->SetCalPhi(phi);
      hitp->SetXYFitFlag(xyFitFlag_);
      hitp->SetCalPosition(xcal, ycal);

      double xcal2 = r*cos(phi*Deg2Rad);
      double ycal2 = r*sin(phi*Deg2Rad);
      double zcal;
      if (zTrackFlag_==0) {
	if (xyFitFlag_==0) {
	  zcal = (xcal2- Bz_)/Az_;
	} else {
	  zcal = (ycal2- Bz_)/Az_;
	}
      } else if (zTrackFlag_==1){
	if (xyFitFlag_==0) {
	  zcal = Az_*xcal2 + Bz_;
	  //ycal = Axy_*xcal + Bxy_;
	} else {
	  zcal = Az_*ycal2 + Bz_;
	  //xcal = Axy_*ycal + Bxy_;
	}
      }

      hitp->SetCalZ(zcal);

    }
  }

  double mean_dE = 0;
  double theta_cft = GetThetaCFT();
  for( std::size_t i=0; i<n; ++i ){
    CFTFiberCluster *hitp = hitArray[i];
    if( hitp ){
      mean_dE += hitp->MaxDeltaE()*sin(theta_cft*TMath::DegToRad());
    }
  }
  for( std::size_t i=0; i<nU; ++i ){
    CFTFiberCluster *hitp = hitArrayU[i];
    if( hitp ){
      mean_dE += hitp->MaxDeltaE();
    }
  }

  mean_dE /= (double)(n+nU);

  double sigma_dE = mean_dE/4.;
  double chisqr_dE = 0;
  for( std::size_t i=0; i<n; ++i ){
    CFTFiberCluster *hitp = hitArray[i];
    if( hitp ){
      double dE_max = hitp->MaxDeltaE();
      double theta = theta_cft;

      chisqr_dE += (dE_max*sin(theta*Deg2Rad) - mean_dE)*(dE_max*sin(theta*Deg2Rad) - mean_dE)/(sigma_dE*sigma_dE);

    }
  }
  for( std::size_t i=0; i<nU; ++i ){
    CFTFiberCluster *hitp = hitArrayU[i];
    if( hitp ){
      double dE_max = hitp->MaxDeltaE();
      chisqr_dE += (dE_max - mean_dE)*(dE_max - mean_dE)/(sigma_dE*sigma_dE);
    }
  }

  chisqr_dE /= (double)(n+nU);
  chisqr_dE_ = chisqr_dE;

  return true;

}

bool CFTLocalTrack::SetCorPhi()
{
  const std::string funcname = "[CFTLocalTrack::SetCorPhi()]";

  std::size_t n = hitArray.size();
  for( std::size_t i=0; i<n; ++i ){
    CFTFiberCluster *hitp = hitArray[i];

    if( hitp ){
      hitp->SetCorPhi();
    }
  }

  return true;
}

bool CFTLocalTrack::SetCalculatedValue_16layer()
{
  const std::string funcname = "[CFTLocalTrack::SetCalculatedValue_16layer()]";
  //const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  std::size_t n = hitArray.size();

  ThreeVector pos1, pos2;

  for( std::size_t i=0; i<n; ++i ){
    CFTFiberCluster *hitp = hitArray[i];
    if( hitp ){
      Int_t lnum = gGeom.GetDetectorId(hitp->DetectorName()+"-"+hitp->PlaneName());
      double r = gGeom.GetLocalZ(lnum);
      double phi, phi1, phi2;

      bool ret = GetCrossPointR(r, &phi1, &phi2);
      if (!ret) {
	//std::cout << funcname << " return at GetCrossPointR" << std::endl;
	return false;
      }
      double phiCal = hitp->GetCalPhi();
      if (std::abs(phi1-phiCal) <= std::abs(phi2-phiCal) && std::abs(phi1-phiCal) < 10.)
	phi = phi1;
      else if (std::abs(phi1-phiCal) >= std::abs(phi2-phiCal) && std::abs(phi2-phiCal) < 10.)
	phi = phi2;
      else if (std::abs(phi1+phiCal-360.)<=10. && (phi1<=15 || phi1 >=345))
	phi = phi1;
      else
	phi = phi2;

      double xcal = r*cos(phi*Deg2Rad);
      double ycal = r*sin(phi*Deg2Rad);
      double zcal;
      if (zTrackFlag_==0) {
	if (xyFitFlag_==0) {
	  zcal = (xcal- Bz_)/Az_;
	} else {
	  zcal = (ycal- Bz_)/Az_;
	}

	if (i==0) {
	  pos1 = ThreeVector(u0_*zcal+x0_, v0_*zcal+y0_, zcal);
	} else if (i==n-1) {
	  pos2 = ThreeVector(u0_*zcal+x0_, v0_*zcal+y0_, zcal);
	}

      } else if (zTrackFlag_==1){
	if (xyFitFlag_==0) {
	  zcal = Az_*xcal + Bz_;
	  ycal = Axy_*xcal + Bxy_;
	} else {
	  zcal = Az_*ycal + Bz_;
	  xcal = Axy_*ycal + Bxy_;
	}
	if (i==0) {
	  pos1 = ThreeVector(xcal, ycal, zcal);
	} else if (i==n-1) {
	  pos2 = ThreeVector(xcal, ycal, zcal);
	}
      }

      hitp->SetCalPhi(phi);
      hitp->SetCalZ(zcal);
      hitp->SetXYFitFlag(xyFitFlag_);
      hitp->SetCalPosition(xcal, ycal);
    }
  }

  Pos0_ = pos1;
  Dir_ = pos2-pos1;

  std::size_t nU = hitArrayU.size();

  for( std::size_t i=0; i<nU; ++i ){
    CFTFiberCluster *hitp = hitArrayU[i];
    if( hitp ){
      Int_t lnum = gGeom.GetDetectorId(hitp->DetectorName()+"-"+hitp->PlaneName());
      double r = gGeom.GetLocalZ(lnum);
      double phi, phi1, phi2;

      bool ret = GetCrossPointR(r, &phi1, &phi2);
      if (!ret) {
	//std::cout << funcname << " return at GetCrossPointR" << std::endl;
	return false;
      }

      double phiCal = hitp->GetCalPhi();
      if (std::abs(phi1-phiCal) <= std::abs(phi2-phiCal) && std::abs(phi1-phiCal) < 10.)
	phi = phi1;
      else if (std::abs(phi1-phiCal) >= std::abs(phi2-phiCal) && std::abs(phi2-phiCal) < 10.)
	phi = phi2;
      else if (std::abs(phi1+phiCal-360.)<=10. && (phi1<=15 || phi1 >=345))
	phi = phi1;
      else
	phi = phi2;
      /*
      if (fabs(phi1-phiCal) <= fabs(phi2-phiCal))
	phi = phi1;
      else
	phi = phi2;
      */

      // modosu 2->1
      double z1 = CalculateZpos2(phi, hitp);
      // double z1 = CalculateZpos(phi, hitp);
      double xycal = Az_*z1 + Bz_;
      double xcal = u0_*z1 + x0_;
      double ycal = v0_*z1 + y0_;

      // use 4 layer tracking value with theta correction
      // comment out here
      //hitp->SetZ(z1);
      hitp->SetCalPhi(phi);
      //hitp->SetCalZ(z1);
      hitp->SetXYFitFlag(xyFitFlag_);
      hitp->SetCalPosition(xcal, ycal);

      double xcal2 = r*cos(phi*Deg2Rad);
      double ycal2 = r*sin(phi*Deg2Rad);
      double zcal;
      if (zTrackFlag_==0) {
	if (xyFitFlag_==0) {
	  zcal = (xcal2- Bz_)/Az_;
	} else {
	  zcal = (ycal2- Bz_)/Az_;
	}
      } else if (zTrackFlag_==1){
	if (xyFitFlag_==0) {
	  zcal = Az_*xcal2 + Bz_;
	  //ycal = Axy_*xcal + Bxy_;
	} else {
	  zcal = Az_*ycal2 + Bz_;
	  //xcal = Axy_*ycal + Bxy_;
	}
      }

      hitp->SetCalZ(zcal);

    }
  }

  return true;

}

#if 0
bool CFTLocalTrack::FindZTrack()
{
  if (hitArrayU.size() != 2)
    return false;

  const double trigPosZ = 220.;

  //const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  int lnum = geomMan.GetDetectorId("CFT_U");
  double r = geomMan.GetLocalZ(lnum);

  double tmp_phi1, tmp_phi2;
  bool ret = GetCrossPointR(r, &tmp_phi1, &tmp_phi2);
  if (!ret) {
    //std::cout << "return at GetCrossPointR (1)" << std::endl;
    return false;
  }

  double phi1, phi2;
  if ((tmp_phi1>=0 && tmp_phi1 <=90) || (tmp_phi1>=270 && tmp_phi1<=360)) {
    phi1 = tmp_phi1;
    phi2 = tmp_phi2;
  } else {
    phi1 = tmp_phi2;
    phi2 = tmp_phi1;
  }

  CFTFiberCluster *cl1, *cl2;
  cl1 = hitArrayU[0];
  cl2 = hitArrayU[1];

  //std::cout << "phi1 : " << phi1 << ", phi2 : " << phi2 << std::endl;
  double z1 = CalculateZpos(phi1, cl1);
  if (z1<0) {
    std::cout << "return at CalculateZpos(phi1, cl1)" << z1 << std::endl;
    return false;
  }
  double z2 = CalculateZpos(phi1, cl2);
  if (z2<0) {
    std::cout << "return at CalculateZpos(phi1, cl2)" << z2 << std::endl;
    return false;
  }

  ThreeVector posU1, posU2;
  if (fabs(z1-trigPosZ) <= fabs(z2-trigPosZ)) {
    posU1 = ThreeVector(r*cos(phi1*Deg2Rad), r*sin(phi1*Deg2Rad), z1);

    cl1->SetCalPosition(r*cos(phi1*Deg2Rad), r*sin(phi1*Deg2Rad));
    cl1->SetCalPhi(phi1);
    cl1->SetCalZ(z1);

    double zz = CalculateZpos(phi2, cl2);
    if (zz<0) {
      std::cout << "return at CalculateZpos(phi2, cl2)" << zz << std::endl;
      return false;
    }
    posU2 = ThreeVector(r*cos(phi2*Deg2Rad), r*sin(phi2*Deg2Rad), zz);

    cl2->SetCalPosition(r*cos(phi2*Deg2Rad), r*sin(phi2*Deg2Rad));
    cl2->SetCalPhi(phi2);
    cl2->SetCalZ(zz);

  } else {
    posU1 = ThreeVector(r*cos(phi1*Deg2Rad), r*sin(phi1*Deg2Rad), z2);
    cl2->SetCalPosition(r*cos(phi1*Deg2Rad), r*sin(phi1*Deg2Rad));
    cl2->SetCalPhi(phi1);
    cl2->SetCalZ(z2);

    double zz = CalculateZpos(phi2, cl1);
    if (zz<0) {
      std::cout << "return at CalculateZpos(phi2, cl1)" << zz << std::endl;
      return false;
    }
    posU2 = ThreeVector(r*cos(phi2*Deg2Rad), r*sin(phi2*Deg2Rad), zz);

    cl1->SetCalPosition(r*cos(phi2*Deg2Rad), r*sin(phi2*Deg2Rad));
    cl1->SetCalPhi(phi2);
    cl1->SetCalZ(zz);
  }


  ThreeVector dirVec = posU1-posU2;
  // x = u0 * z + x0
  // y = v0 * z + y0

  u0_ = dirVec.x()/dirVec.z();
  v0_ = dirVec.y()/dirVec.z();
  x0_ = posU1.x()-u0_*posU1.z();
  y0_ = posU1.y()-v0_*posU1.z();

  zTrackFlag_ = 1;

  return true;
}
#endif

bool CFTLocalTrack::CalcNormalizedDE( void )
{
  int nphi = hitArray.size();
  for (int i=0; i<nphi; i++) {
    CFTFiberCluster *cl = hitArray[i];
    double r = cl->MeanR();
    double phi = cl->MeanPhi();

    double len;
    if (!GetPathLengthInFiber(r, phi, &len)) {
      //std::cout << "CFTLocalTrack::CalcNormalizedDE , error in GetPathLengthInFiber r " << r
      //<< ", phi " << phi << std::endl;
      return false;
    }
    /*
    double meanId = cl->MeanSeg();
    std::cout << "PHI : " << i << ", r " << r << ", phi " << phi
	      << ", MeanId " << meanId << ", len " << len << std::endl;
    */
    cl->SetPathLength(len);
    pathlength_ += len;
    TotalDEHiGain_  += cl->TotalDeltaE_Hi();
    TotalDELowGain_ += cl->TotalDeltaE();
    MaxDEHiGain_  += cl->MaxDeltaE_Hi();
    MaxDELowGain_ += cl->MaxDeltaE();
    MaxDEHiGain_phi_  += cl->MaxDeltaE_Hi();
    MaxDELowGain_phi_ += cl->MaxDeltaE();

    TotalDELowGain_phi_ += cl->TotalDeltaE();
  }

  int nu = hitArrayU.size();
  for (int i=0; i<nu; i++) {
    CFTFiberCluster *cl = hitArrayU[i];
    double r = cl->MeanR();
    double phi = cl->GetCalPhi();

    double len;
    if (!GetPathLengthInFiber( r, phi, &len)) {
      //std::cout << "CFTLocalTrack::CalcNormalizedDE , error in GetPathLengthInFiber r " << r
      //<< ", phi " << phi << std::endl;
      return false;
    }
    /*
    double meanId = cl->MeanSeg();
    std::cout << "U : " << i << ", r " << r << ", phi " << phi
	      << ", MeanId " << meanId << ", len " << len << std::endl;
    */
    cl->SetPathLength(len);
    pathlength_ += len;
    TotalDEHiGain_  += cl->TotalDeltaE_Hi();
    TotalDELowGain_ += cl->TotalDeltaE();
    MaxDEHiGain_  += cl->MaxDeltaE_Hi();
    MaxDELowGain_ += cl->MaxDeltaE();
    MaxDEHiGain_uv_  += cl->MaxDeltaE_Hi();
    MaxDELowGain_uv_ += cl->MaxDeltaE();

    TotalDELowGain_uv_ += cl->TotalDeltaE();
  }

  //std::cout << "Total PathLength " << pathlength_ << std::endl;

  return true;

}


bool CFTLocalTrack::GetPathLengthInFiber(double r, double phi, double *len)
{
  double x0 = r*cos(phi*Deg2Rad);
  double y0 = r*sin(phi*Deg2Rad);

  double phi1_in1, phi1_in2;
  double phi1_out1, phi1_out2;
  if (!GetCrossPointR(r-0.75/2, &phi1_in1, &phi1_in2)) {
    //std::cout << "CFTLocalTrack::GetPathLengthInFiber Error in GetCrossPointR(r-0.75/2, &phi1_in1, &phi1_in2)"
    //	      << std::endl;
    return false;
  }

  if (!GetCrossPointR(r+0.75/2, &phi1_out1, &phi1_out2)) {
    //std::cout << "CFTLocalTrack::GetPathLengthInFiber Error in GetCrossPointR(r+0.75/2, &phi1_out1, &phi1_out2)"
    //<< std::endl;
    return false;

  }

  double x1 = r*cos(phi1_in1*Deg2Rad);
  double y1 = r*sin(phi1_in1*Deg2Rad);

  double x2 = r*cos(phi1_in2*Deg2Rad);
  double y2 = r*sin(phi1_in2*Deg2Rad);

  double dist1 = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));
  double dist2 = sqrt((x0-x2)*(x0-x2) + (y0-y2)*(y0-y2));

  /*
  if (dist1 > 2 && dist2 > 2)
    std::cout << "dist1 : " << dist1 << ", dist2 : " << dist2 << std::endl;
  */
  if (dist1<dist2) {
    double x1_in1 = (r-0.75/2)*cos(phi1_in1*Deg2Rad);
    double y1_in1 = (r-0.75/2)*sin(phi1_in1*Deg2Rad);
    double z1_in1  = (x1_in1-x0_)/u0_;
    double z1_in1_y  = (y1_in1-y0_)/v0_;
    double x1_out1 = (r+0.75/2)*cos(phi1_out1*Deg2Rad);
    double y1_out1 = (r+0.75/2)*sin(phi1_out1*Deg2Rad);
    double z1_out1  = (x1_out1-x0_)/u0_;
    double z1_out1_y  = (y1_out1-y0_)/v0_;

    if (std::abs(z1_in1-z1_in1_y)>1 || std::abs(z1_out1-z1_out1_y)>1 ) {
      std::cout << "z1_in1 : " << z1_in1 << ", z1_in1_y : " << z1_in1_y << std::endl;
      std::cout << "z1_out1 : " << z1_out1 << ", z1_out1_y : " << z1_out1_y << std::endl;
    }

    *len = sqrt((x1_in1-x1_out1)*(x1_in1-x1_out1) +
		(y1_in1-y1_out1)*(y1_in1-y1_out1) +
		(z1_in1-z1_out1)*(z1_in1-z1_out1));

    /*
    if (r>30 && r<45) {
      std::cout << "(x_in1, y_in1, z_in1) = ( " << x1_in1 << ", " << y1_in1
		<< ", " << z1_in1 << ")" << std::endl;
      std::cout << "(x_out1, y_out1, z_out1) = ( " << x1_out1 << ", " << y1_out1
		<< ", " << z1_out1 << ")" << std::endl;
      std::cout << "Pathlength : " << *len << std::endl;
    }
    */

    return true;

  } else {
    double x1_in2 = (r-0.75/2)*cos(phi1_in2*Deg2Rad);
    double y1_in2 = (r-0.75/2)*sin(phi1_in2*Deg2Rad);
    double z1_in2  = (x1_in2-x0_)/u0_;
    double z1_in2_y  = (y1_in2-y0_)/v0_;
    double x1_out2 = (r+0.75/2)*cos(phi1_out2*Deg2Rad);
    double y1_out2 = (r+0.75/2)*sin(phi1_out2*Deg2Rad);
    double z1_out2  = (x1_out2-x0_)/u0_;
    double z1_out2_y  = (y1_out2-y0_)/v0_;

    if (std::abs(z1_in2-z1_in2_y)>1 || std::abs(z1_out2-z1_out2_y)>1 ) {
      std::cout << "z1_in2 : " << z1_in2 << ", z1_in2_y : " << z1_in2_y << std::endl;
      std::cout << "z1_out2 : " << z1_out2 << ", z1_out2_y : " << z1_out2_y << std::endl;
    }

    *len = sqrt((x1_in2-x1_out2)*(x1_in2-x1_out2) +
		(y1_in2-y1_out2)*(y1_in2-y1_out2) +
		(z1_in2-z1_out2)*(z1_in2-z1_out2));

    /*
    if (r>30 && r<45) {
      std::cout << "(x_in2, y_in2, z_in2) = ( " << x1_in2 << ", " << y1_in2
		<< ", " << z1_in2 << ")" << std::endl;
      std::cout << "(x_out2, y_out2, z_out2) = ( " << x1_out2 << ", " << y1_out2
		<< ", " << z1_out2 << ")" << std::endl;
      std::cout << "Pathlength : " << *len << std::endl;
    }
    */

    return true;
  }


  //std::cout << "CFTLocalTrack::GetPathLengthInFiber Strange end"
  //<< std::endl;
  return false;

}

bool CFTLocalTrack::GetCrossPointR(double r, double *phi1, double *phi2)
{
  double a = Axy_;
  double b = Bxy_;
  double x1, y1, x2, y2;

  if (xyFitFlag_==0) {
    double hanbetsu = a*a*b*b-(1.+a*a)*(b*b-r*r);
    if (hanbetsu<0)
      return false;
    else {
      x1 = (-a*b+sqrt(hanbetsu))/(1+a*a);
      x2 = (-a*b-sqrt(hanbetsu))/(1+a*a);
      y1 = Axy_*x1+Bxy_;
      y2 = Axy_*x2+Bxy_;
    }
  } else if (xyFitFlag_==1) {
    double hanbetsu = a*a*b*b-(1.+a*a)*(b*b-r*r);
    if (hanbetsu<0)
      return false;
    else {
      y1 = (-a*b+sqrt(hanbetsu))/(1+a*a);
      y2 = (-a*b-sqrt(hanbetsu))/(1+a*a);
      x1 = Axy_*y1+Bxy_;
      x2 = Axy_*y2+Bxy_;
    }
  }

  *phi1 = calcPhi(x1, y1);
  *phi2 = calcPhi(x2, y2);

  return true;
}

bool CFTLocalTrack::GetCrossPointR(double r, double *phi)
{
  double a = Axy_;
  double b = Bxy_;
  double x1, y1, x2, y2;

  if (xyFitFlag_==0) {
    double hanbetsu = a*a*b*b-(1.+a*a)*(b*b-r*r);
    if (hanbetsu<0)
      return false;
    else {
      x1 = (-a*b+sqrt(hanbetsu))/(1+a*a);
      x2 = (-a*b-sqrt(hanbetsu))/(1+a*a);
      y1 = Axy_*x1+Bxy_;
      y2 = Axy_*x2+Bxy_;
    }
  } else if (xyFitFlag_==1) {
    double hanbetsu = a*a*b*b-(1.+a*a)*(b*b-r*r);
    if (hanbetsu<0)
      return false;
    else {
      y1 = (-a*b+sqrt(hanbetsu))/(1+a*a);
      y2 = (-a*b-sqrt(hanbetsu))/(1+a*a);
      x1 = Axy_*y1+Bxy_;
      x2 = Axy_*y2+Bxy_;
    }
  }

  double phi1 = calcPhi(x1, y1);
  double phi2 = calcPhi(x2, y2);

  int nPhiHit = GetNHit();

  bool flagZeroCross = false;
  bool flagZone1 = false;
  bool flagZone4 = false;
  for (int i=0; i<nPhiHit; i++) {
    double phiVal = hitArray[i]->MeanPhi();
    if (phiVal>=0 && phiVal<=90)
      flagZone1 = true;
    else if (phiVal>=270 && phiVal<=360)
      flagZone4 = true;
  }

  if (flagZone1 && flagZone4)
    flagZeroCross = true;

  double meanPhi=0;
  for (int i=0; i<nPhiHit; i++) {
    double phiVal = hitArray[i]->MeanPhi();
    //std::cout << phiVal << "-" << phiCal << "  ";
    if (flagZeroCross && phiVal >=270 && phiVal <=360)
      phiVal -= 360.;

    meanPhi += phiVal;
  }
  meanPhi /= nPhiHit;
  //std::cout << std::endl;

  if (!flagZeroCross) {
    if (std::abs(phi1-meanPhi)<std::abs(phi2-meanPhi)) {
      *phi = phi1;
    } else {
      *phi = phi2;
    }
  } else {
    if (phi1>=270&&phi1<=360)
      phi1 -= 360.;
    if (phi2>=270&&phi2<=360)
      phi2 -= 360.;

    if (std::abs(phi1-meanPhi)<std::abs(phi2-meanPhi)) {
      if (phi1<0)
	phi1 += 360.;

      *phi = phi1;
    } else {
      if (phi2<0)
	phi2 += 360.;

      *phi = phi2;
    }
  }
  /*
  if (fabs(*phi-meanPhi)>10)
    std::cout << "MeanPhi " << meanPhi
	      << "phi1 " << phi1 << ", phi2 " << phi2 << "phi " << *phi << std::endl;
  */

  return true;
}

double CFTLocalTrack::calcPhi(double x, double y)
{
  if (x>=0 && y>=0)
    return atan(y/x)*Rad2Deg;
  else if (x<0 && y>=0)
    return 180. + atan(y/x)*Rad2Deg;
  else if (x<0 && y<0)
    return 180. + atan(y/x)*Rad2Deg;
  else if (x>=0 && y<0)
    return 360. + atan(y/x)*Rad2Deg;

}

double CFTLocalTrack::CalculateUVpos(double phi, CFTFiberCluster *cl)
{
  double z0 = cl->MeanZ0();
  double slope = cl->Slope();

  double a1 = z0+slope*phi;
  double a2 = (z0-400.)+slope*phi;
  if (slope<0)
    a2 = (z0+400.)+slope*phi;

  double theta = atan(slope)*Rad2Deg;

  if (a1>=0 && a1<=400)
    return (z0 - offsetCATCH)*cos(theta*Deg2Rad);
  else if (a2>=0 && a2<=400)
    if (slope>0)
      return (z0-400. - offsetCATCH)*cos(theta*Deg2Rad);
    else
      return (z0+400. - offsetCATCH)*cos(theta*Deg2Rad);
  else
    return -100.;
}

double CFTLocalTrack::CalculateZpos(double phi, CFTFiberCluster *cl, double mean_z)
{
  double z0 = cl->MeanZ0();
  double slope = cl->Slope();
  double time = cl->CMeanTimeMaxCluster();

  double a1 = z0+slope*phi;

  double a2 = (z0-400.)+slope*phi;
  if (slope<0)
    a2 = (z0+400.)+slope*phi;


  if (mean_z>-300 && mean_z <300) {
    // mean_z to CATCH frame
    mean_z += offsetCATCH;

    double z_glo = -99999.;
    if (std::abs(a1-mean_z) < std::abs(a2-mean_z)) {
      // std::cout << "a1 selected ( " << a1 << "), mean_z ( "
      // 		<< mean_z << " ), a2 ( " << a2 << " )" << std::endl;
      z_glo = (a1 - offsetCATCH);
    } else {
      // std::cout << "a2 selected ( " << a2 << "), mean_z ( "
      // 		<< mean_z << " ), a1 ( " << a1 << " )" << std::endl;
      z_glo =  (a2 - offsetCATCH);
    }

    return z_glo;
  }

  if (a1>=0 && a1<=400) {

    //cosmic ray no use
    if (a1>=360 && time<=-1)
      a1 = a1-400;
    else if (a1<=20 && time>=0.5)
      a1 = a1+400;

    return (a1 - offsetCATCH);
  } else if (a2>=0 && a2<=400) {

    //cosmic ray no use
    if (a2>=360 && time<=-1)
      a2 = a2-400;
    else if (a2<=20 && time>=0.5)
      a2 = a2+400;

    return (a2 - offsetCATCH);
  } else
    return -100.;
}

double CFTLocalTrack::CalculateZpos2(double phi, CFTFiberCluster *cl, double mean_z)
{
  // mean_z to CATCH frame
  mean_z += offsetCATCH;

  double z = cl->GetCorZ(phi, mean_z, theta_);
  double time = cl->CMeanTimeMaxCluster();

  //cosmic ray no use
  if (z>=360 && time<=-1)
    z = z-400;
  else if (z<=20 && time>=0.5)
    z = z+400;

  return z - offsetCATCH;
}

/*
bool CFTLocalTrack::ReCalc( bool applyRecursively )
{
  static const std::string funcname = "[CFTLocalTrack::ReCalc]";

  std::size_t n = hitArray.size();
  for( std::size_t i=0; i<n; ++i ){
    CFTFiberCluster *hitp = hitArray[i];
    if( hitp ) hitp->ReCalc( applyRecursively );
  }

  bool ret=DoFit();
  if( !ret ){
    std::cerr << funcname << ": Recalculation fails" << std::endl;
  }
  return ret;
}
*/

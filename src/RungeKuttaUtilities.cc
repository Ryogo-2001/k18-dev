// -*- C++ -*-

// Runge-Kutta Routines
// Ref.) NIM 160 (1979) 43 - 48

#include "RungeKuttaUtilities.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include <std_ostream.hh>

#include "ConfMan.hh"
#include "DatabasePDG.hh"
#include "Kinematics.hh"
#include "DCGeomMan.hh"
#include "DCGeomRecord.hh"
#include "EventDisplay.hh"
#include "Exception.hh"
#include "FieldMan.hh"
#include "FuncName.hh"
#include "KuramaTrack.hh"
#include "TPCRKTrack.hh"
#include "HSTrack.hh"
#include "PrintHelper.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCPadHelper.hh"

namespace
{
const auto& gGeom   = DCGeomMan::GetInstance();
auto&       gEvDisp = EventDisplay::GetInstance();
const auto& gField  = FieldMan::GetInstance();
const Int_t& IdHS = gGeom.DetectorId("HS");
// const Int_t& IdTOF    = gGeom.DetectorId("TOF");
const Int_t& IdTOF_UX = gGeom.DetectorId("TOF-UX");
const Int_t& IdTOF_UY = gGeom.DetectorId("TOF-UY");
const Int_t& IdTOF_DX = gGeom.DetectorId("TOF-DX");
const Int_t& IdTOF_DY = gGeom.DetectorId("TOF-DY");
const Int_t& IdTarget = gGeom.DetectorId("Target");
const Int_t& IdVPHTOF = gGeom.DetectorId("VPHTOF");
const Int_t& IdRKINIT = gGeom.DetectorId("RKINIT");
const Int_t& IdBAC = gGeom.DetectorId("BAC");
const Int_t& IdBH2 = gGeom.DetectorId("BH2");
const Int_t& IdTPCGasVessel_U = gGeom.DetectorId("VesselU");
const Int_t& IdTPCGasVessel_D = gGeom.DetectorId("VesselD");
const Int_t& IdK18TPCGasVessel_U = gGeom.DetectorId("K18VesselU");
const Int_t& IdK18TPCGasVessel_D = gGeom.DetectorId("K18VesselD");
const Int_t& IdK18Target = gGeom.DetectorId("K18Target");
const Int_t& IdHTOF = gGeom.DetectorId("HTOF");
const Int_t& IdTgtVP1 = gGeom.DetectorId("TargetVP1");
const Int_t& IdTgtVP2 = gGeom.DetectorId("TargetVP2");
const Int_t& IdTgtVP3 = gGeom.DetectorId("TargetVP3");
const Int_t& IdTgtVP4 = gGeom.DetectorId("TargetVP4");
const Int_t& IdTgtVP5 = gGeom.DetectorId("TargetVP5");
const Int_t& IdTgtVP6 = gGeom.DetectorId("TargetVP6");
const Double_t CHLB     = 2.99792458E-4;
const Double_t Polarity = 1.;
}

#define WARNOUT 0
#define ExactFFTreat 1
#define dEdxCorrection 1 //dEdx Correction (BH2 && HTOF & TPC)

//_____________________________________________________________________________
void
RKFieldIntegral::Print(std::ostream &ost) const
{
  PrintHelper helper(3, std::ios::scientific, ost);
  ost << std::setw(9) << kx << " "
      << std::setw(9) << ky << std::endl;
  ost << std::setw(9) << axu << " "
      << std::setw(9) << axv << " "
      << std::setw(9) << ayu << " "
      << std::setw(9) << ayv << std::endl;
  ost << std::setw(9) << cxx << " "
      << std::setw(9) << cxy << " "
      << std::setw(9) << cyx << " "
      << std::setw(9) << cyy << std::endl;
}

//_____________________________________________________________________________
void
RKDeltaFieldIntegral::Print(std::ostream &ost) const
{
  PrintHelper helper(3, std::ios::scientific, ost);
  ost << std::setw(9) << dkxx << " "
      << std::setw(9) << dkxy << " "
      << std::setw(9) << dkxu << " "
      << std::setw(9) << dkxv << " "
      << std::setw(9) << dkxq << std::endl;
  ost << std::setw(9) << dkyx << " "
      << std::setw(9) << dkyy << " "
      << std::setw(9) << dkyu << " "
      << std::setw(9) << dkyv << " "
      << std::setw(9) << dkyq << std::endl;
}

//_____________________________________________________________________________
void
RKCordParameter::Print(std::ostream &ost) const
{
  PrintHelper helper(3, std::ios::scientific, ost);
}

//_____________________________________________________________________________
RKCordParameter::RKCordParameter(const ThreeVector &pos,
                                 const ThreeVector &mom)
  : x(pos.x()), y(pos.y()), z(pos.z()),
    u(mom.x()/mom.z()), v(mom.y()/mom.z())
{
  Double_t p = mom.Mag();
  q = -Polarity/p;
}

//_____________________________________________________________________________
ThreeVector
RKCordParameter::MomentumInGlobal() const
{
  Double_t p  = -Polarity/q;
  Double_t pz = -std::abs(p)/std::sqrt(1.+u*u+v*v);
  return ThreeVector(pz*u, pz*v, pz);
}

//_____________________________________________________________________________
void
RKCordParameter::AddMomentum(Double_t add)
{
  Double_t p = -Polarity/q;
  Double_t mag = (std::abs(p) + add);
  q *= std::abs(p)/mag;
}

//_____________________________________________________________________________
void
RKTrajectoryPoint::Print(std::ostream &ost) const
{
  PrintHelper helper(3, std::ios::scientific, ost);
}

//_____________________________________________________________________________
RKFieldIntegral
RK::CalcFieldIntegral(Double_t U, Double_t V, Double_t Q, const ThreeVector &B)
{
  Double_t fac = std::sqrt(1.+U*U+V*V);
  Double_t f1  = U*V*B.x() - (1.+U*U)*B.y() + V*B.z();
  Double_t f2  = (1.+V*V)*B.x() - U*V*B.y() - U*B.z();

  Double_t axu = U/fac*f1 + fac*(V*B.x()-2.*U*B.y());
  Double_t axv = V/fac*f1 + fac*(U*B.x()+B.z());
  Double_t ayu = U/fac*f2 - fac*(V*B.y()+B.z());
  Double_t ayv = V/fac*f2 + fac*(2.*V*B.x()-U*B.y());

  Double_t qfac = Q*CHLB;

  return RKFieldIntegral(fac*f1*qfac, fac*f2*qfac,
                         axu*qfac, axv*qfac, ayu*qfac, ayv*qfac);

}

//_____________________________________________________________________________
RKFieldIntegral
RK::CalcFieldIntegral(Double_t U, Double_t V, Double_t Q, const ThreeVector &B,
                      const ThreeVector &dBdX, const ThreeVector &dBdY)
{
  Double_t fac = std::sqrt(1.+U*U+V*V);
  Double_t f1  = U*V*(B.x()) - (1.+U*U)*(B.y()) + V*(B.z());
  Double_t f2  = (1.+V*V)*(B.x()) - U*V*(B.y()) - U*(B.z());

  Double_t axu = U/fac*f1 + fac*(V*B.x()-2.*U*B.y());
  Double_t axv = V/fac*f1 + fac*(U*B.x()+B.z());
  Double_t ayu = U/fac*f2 - fac*(V*B.y()+B.z());
  Double_t ayv = V/fac*f2 + fac*(2.*V*B.x()-U*B.y());

  Double_t cxx = U*V*(dBdX.x()) - (1.+U*U)*(dBdX.y()) + V*(dBdX.z());
  Double_t cxy = U*V*(dBdY.x()) - (1.+U*U)*(dBdY.y()) + V*(dBdY.z());
  Double_t cyx = (1.+V*V)*(dBdX.x()) - U*V*(dBdX.y()) - U*(dBdX.z());
  Double_t cyy = (1.+V*V)*(dBdY.x()) - U*V*(dBdY.y()) - U*(dBdY.z());

  Double_t qfac = Q*CHLB;

  return RKFieldIntegral(fac*f1*qfac, fac*f2*qfac,
                         axu*qfac, axv*qfac, ayu*qfac, ayv*qfac,
                         fac*cxx*qfac, fac*cxy*qfac,
                         fac*cyx*qfac, fac*cyy*qfac);
}

//_____________________________________________________________________________
RKDeltaFieldIntegral
RK::CalcDeltaFieldIntegral(const RKTrajectoryPoint &prevPoint,
                           const RKFieldIntegral &intg)
{
  Double_t dkxx = intg.axu*prevPoint.dudx + intg.axv*prevPoint.dvdx
    + intg.cxx*prevPoint.dxdx + intg.cxy*prevPoint.dydx;
  Double_t dkxy = intg.axu*prevPoint.dudy + intg.axv*prevPoint.dvdy
    + intg.cxx*prevPoint.dxdy + intg.cxy*prevPoint.dydy;
  Double_t dkxu = intg.axu*prevPoint.dudu + intg.axv*prevPoint.dvdu
    + intg.cxx*prevPoint.dxdu + intg.cxy*prevPoint.dydu;
  Double_t dkxv = intg.axu*prevPoint.dudv + intg.axv*prevPoint.dvdv
    + intg.cxx*prevPoint.dxdv + intg.cxy*prevPoint.dydv;
  Double_t dkxq = intg.kx/prevPoint.r.q
    + intg.axu*prevPoint.dudq + intg.axv*prevPoint.dvdq
    + intg.cxx*prevPoint.dxdq + intg.cxy*prevPoint.dydq;

  Double_t dkyx = intg.ayu*prevPoint.dudx + intg.ayv*prevPoint.dvdx
    + intg.cyx*prevPoint.dxdx + intg.cyy*prevPoint.dydx;
  Double_t dkyy = intg.ayu*prevPoint.dudy + intg.ayv*prevPoint.dvdy
    + intg.cyx*prevPoint.dxdy + intg.cyy*prevPoint.dydy;
  Double_t dkyu = intg.ayu*prevPoint.dudu + intg.ayv*prevPoint.dvdu
    + intg.cyx*prevPoint.dxdu + intg.cyy*prevPoint.dydu;
  Double_t dkyv = intg.ayu*prevPoint.dudv + intg.ayv*prevPoint.dvdv
    + intg.cyx*prevPoint.dxdv + intg.cyy*prevPoint.dydv;
  Double_t dkyq = intg.ky/prevPoint.r.q
    + intg.ayu*prevPoint.dudq + intg.ayv*prevPoint.dvdq
    + intg.cyx*prevPoint.dxdq + intg.cyy*prevPoint.dydq;

  return RKDeltaFieldIntegral(dkxx, dkxy, dkxu, dkxv, dkxq,
                              dkyx, dkyy, dkyu, dkyv, dkyq);
}

//_____________________________________________________________________________
RKDeltaFieldIntegral
RK::CalcDeltaFieldIntegral(const RKTrajectoryPoint    &prevPoint,
                           const RKFieldIntegral      &intg,
                           const RKDeltaFieldIntegral &dIntg1,
                           const RKDeltaFieldIntegral &dIntg2,
                           Double_t StepSize)
{
  Double_t h  = StepSize;
  Double_t h2 = StepSize*StepSize;

  Double_t dkxx
    = intg.axu*(prevPoint.dudx + h*dIntg1.dkxx)
    + intg.axv*(prevPoint.dvdx + h*dIntg1.dkyx)
    + intg.cxx*(prevPoint.dxdx + h*prevPoint.dudx + 0.5*h2*dIntg2.dkxx)
    + intg.cxy*(prevPoint.dydx + h*prevPoint.dvdx + 0.5*h2*dIntg2.dkyx);
  Double_t dkxy
    = intg.axu*(prevPoint.dudy + h*dIntg1.dkxy)
    + intg.axv*(prevPoint.dvdy + h*dIntg1.dkyy)
    + intg.cxx*(prevPoint.dxdy + h*prevPoint.dudy + 0.5*h2*dIntg2.dkxy)
    + intg.cxy*(prevPoint.dydy + h*prevPoint.dvdy + 0.5*h2*dIntg2.dkyy);
  Double_t dkxu
    = intg.axu*(prevPoint.dudu + h*dIntg1.dkxu)
    + intg.axv*(prevPoint.dvdu + h*dIntg1.dkyu)
    + intg.cxx*(prevPoint.dxdu + h*prevPoint.dudu + 0.5*h2*dIntg2.dkxu)
    + intg.cxy*(prevPoint.dydu + h*prevPoint.dvdu + 0.5*h2*dIntg2.dkyu);
  Double_t dkxv
    = intg.axu*(prevPoint.dudv + h*dIntg1.dkxv)
    + intg.axv*(prevPoint.dvdv + h*dIntg1.dkyv)
    + intg.cxx*(prevPoint.dxdv + h*prevPoint.dudv + 0.5*h2*dIntg2.dkxv)
    + intg.cxy*(prevPoint.dydv + h*prevPoint.dvdv + 0.5*h2*dIntg2.dkyv);
  Double_t dkxq = intg.kx/prevPoint.r.q
    + intg.axu*(prevPoint.dudq + h*dIntg1.dkxq)
    + intg.axv*(prevPoint.dvdq + h*dIntg1.dkyq)
    + intg.cxx*(prevPoint.dxdq + h*prevPoint.dudq + 0.5*h2*dIntg2.dkxq)
    + intg.cxy*(prevPoint.dydq + h*prevPoint.dvdq + 0.5*h2*dIntg2.dkyq);

  Double_t dkyx
    = intg.ayu*(prevPoint.dudx + h*dIntg1.dkxx)
    + intg.ayv*(prevPoint.dvdx + h*dIntg1.dkyx)
    + intg.cyx*(prevPoint.dxdx + h*prevPoint.dudx + 0.5*h2*dIntg2.dkxx)
    + intg.cyy*(prevPoint.dydx + h*prevPoint.dvdx + 0.5*h2*dIntg2.dkyx);
  Double_t dkyy
    = intg.ayu*(prevPoint.dudy + h*dIntg1.dkxy)
    + intg.ayv*(prevPoint.dvdy + h*dIntg1.dkyy)
    + intg.cyx*(prevPoint.dxdy + h*prevPoint.dudy + 0.5*h2*dIntg2.dkxy)
    + intg.cyy*(prevPoint.dydy + h*prevPoint.dvdy + 0.5*h2*dIntg2.dkyy);
  Double_t dkyu
    = intg.ayu*(prevPoint.dudu + h*dIntg1.dkxu)
    + intg.ayv*(prevPoint.dvdu + h*dIntg1.dkyu)
    + intg.cyx*(prevPoint.dxdu + h*prevPoint.dudu + 0.5*h2*dIntg2.dkxu)
    + intg.cyy*(prevPoint.dydu + h*prevPoint.dvdu + 0.5*h2*dIntg2.dkyu);
  Double_t dkyv
    = intg.ayu*(prevPoint.dudv + h*dIntg1.dkxv)
    + intg.ayv*(prevPoint.dvdv + h*dIntg1.dkyv)
    + intg.cyx*(prevPoint.dxdv + h*prevPoint.dudv + 0.5*h2*dIntg2.dkxv)
    + intg.cyy*(prevPoint.dydv + h*prevPoint.dvdv + 0.5*h2*dIntg2.dkyv);
  Double_t dkyq = intg.ky/prevPoint.r.q
    + intg.ayu*(prevPoint.dudq + h*dIntg1.dkxq)
    + intg.ayv*(prevPoint.dvdq + h*dIntg1.dkyq)
    + intg.cyx*(prevPoint.dxdq + h*prevPoint.dudq + 0.5*h2*dIntg2.dkxq)
    + intg.cyy*(prevPoint.dydq + h*prevPoint.dvdq + 0.5*h2*dIntg2.dkyq);

  return RKDeltaFieldIntegral(dkxx, dkxy, dkxu, dkxv, dkxq,
                              dkyx, dkyy, dkyu, dkyv, dkyq);
}

//____________________________________________________________________________
RKTrajectoryPoint
RK::TraceOneStep(Double_t StepSize, const RKTrajectoryPoint &prevPoint)
{
  Double_t pre_x = prevPoint.r.x;
  Double_t pre_y = prevPoint.r.y;
  Double_t pre_z = prevPoint.r.z;
  Double_t pre_u = prevPoint.r.u;
  Double_t pre_v = prevPoint.r.v;
  Double_t pre_q = prevPoint.r.q;
  Double_t dr    = StepSize/std::sqrt(1.+pre_u*pre_u+pre_v*pre_v);

  ThreeVector Z1 = prevPoint.PositionInGlobal();
  ThreeVector B1 = gField.GetField(Z1);
#ifdef ExactFFTreat
  ThreeVector dBdX1 = gField.GetdBdX(Z1);
  ThreeVector dBdY1 = gField.GetdBdY(Z1);
  RKFieldIntegral f1 =
    RK::CalcFieldIntegral(pre_u, pre_v, pre_q,
                          B1, dBdX1, dBdY1);
#else
  RKFieldIntegral f1 =
    RK::CalcFieldIntegral(pre_u, pre_v, pre_q, B1);
#endif
  RKDeltaFieldIntegral df1 =
    RK::CalcDeltaFieldIntegral(prevPoint, f1);

  ThreeVector Z2 = Z1 +
    ThreeVector(0.5*dr,
                0.5*dr*pre_u + 0.125*dr*dr*f1.kx,
                0.5*dr*pre_v + 0.125*dr*dr*f1.ky);
  ThreeVector B2 = gField.GetField(Z2);
#ifdef ExactFFTreat
  ThreeVector dBdX2 = gField.GetdBdX(Z2);
  ThreeVector dBdY2 = gField.GetdBdY(Z2);
  RKFieldIntegral f2 =
    RK::CalcFieldIntegral(pre_u + 0.5*dr*f1.kx,
                          pre_v + 0.5*dr*f1.ky,
                          pre_q, B2, dBdX2, dBdY2);
#else
  RKFieldIntegral f2 =
    RK::CalcFieldIntegral(pre_u + 0.5*dr*f1.kx,
                          pre_v + 0.5*dr*f1.ky,
                          pre_q, B2);
#endif
  RKDeltaFieldIntegral df2 =
    RK::CalcDeltaFieldIntegral(prevPoint, f2, df1, df1, 0.5*dr);

#ifdef ExactFFTreat
  RKFieldIntegral f3 =
    RK::CalcFieldIntegral(pre_u + 0.5*dr*f2.kx,
                          pre_v + 0.5*dr*f2.ky,
                          pre_q, B2, dBdX2, dBdY2);
#else
  RKFieldIntegral f3 =
    RK::CalcFieldIntegral(pre_u + 0.5*dr*f2.kx,
                          pre_v + 0.5*dr*f2.ky,
                          pre_q, B2);
#endif
  RKDeltaFieldIntegral df3 =
    RK::CalcDeltaFieldIntegral(prevPoint, f3, df2, df1, 0.5*dr);

  ThreeVector Z4 = Z1 +
    ThreeVector(dr,
                dr*pre_u + 0.5*dr*dr*f3.kx,
                dr*pre_v + 0.5*dr*dr*f3.ky);
  ThreeVector B4 = gField.GetField(Z4);
#ifdef ExactFFTreat
  ThreeVector dBdX4 = gField.GetdBdX(Z4);
  ThreeVector dBdY4 = gField.GetdBdY(Z4);
  RKFieldIntegral f4 =
    RK::CalcFieldIntegral(pre_u + dr*f3.kx,
                          pre_v + dr*f3.ky,
                          pre_q, B4, dBdX4, dBdY4);
#else
  RKFieldIntegral f4 =
    RK::CalcFieldIntegral(pre_u + dr*f3.kx,
                          pre_v + dr*f3.ky,
                          pre_q, B4);
#endif
  RKDeltaFieldIntegral df4 =
    RK::CalcDeltaFieldIntegral(prevPoint, f4, df3, df3, dr);

  Double_t z = pre_z + dr;
  Double_t x = pre_x + dr*pre_u
    + 1./6.*dr*dr*(f1.kx+f2.kx+f3.kx);
  Double_t y = pre_y + dr*pre_v
    + 1./6.*dr*dr*(f1.ky+f2.ky+f3.ky);
  Double_t u = pre_u + 1./6.*dr*(f1.kx+2.*(f2.kx+f3.kx)+f4.kx);
  Double_t v = pre_v + 1./6.*dr*(f1.ky+2.*(f2.ky+f3.ky)+f4.ky);

  Double_t dxdx = prevPoint.dxdx + dr*prevPoint.dudx
    + 1./6.*dr*dr*(df1.dkxx+df2.dkxx+df3.dkxx);
  Double_t dxdy = prevPoint.dxdy + dr*prevPoint.dudy
    + 1./6.*dr*dr*(df1.dkxy+df2.dkxy+df3.dkxy);
  Double_t dxdu = prevPoint.dxdu + dr*prevPoint.dudu
    + 1./6.*dr*dr*(df1.dkxu+df2.dkxu+df3.dkxu);
  Double_t dxdv = prevPoint.dxdv + dr*prevPoint.dudv
    + 1./6.*dr*dr*(df1.dkxv+df2.dkxv+df3.dkxv);
  Double_t dxdq = prevPoint.dxdq + dr*prevPoint.dudq
    + 1./6.*dr*dr*(df1.dkxq+df2.dkxq+df3.dkxq);

  Double_t dydx = prevPoint.dydx + dr*prevPoint.dvdx
    + 1./6.*dr*dr*(df1.dkyx+df2.dkyx+df3.dkyx);
  Double_t dydy = prevPoint.dydy + dr*prevPoint.dvdy
    + 1./6.*dr*dr*(df1.dkyy+df2.dkyy+df3.dkyy);
  Double_t dydu = prevPoint.dydu + dr*prevPoint.dvdu
    + 1./6.*dr*dr*(df1.dkyu+df2.dkyu+df3.dkyu);
  Double_t dydv = prevPoint.dydv + dr*prevPoint.dvdv
    + 1./6.*dr*dr*(df1.dkyv+df2.dkyv+df3.dkyv);
  Double_t dydq = prevPoint.dydq + dr*prevPoint.dvdq
    + 1./6.*dr*dr*(df1.dkyq+df2.dkyq+df3.dkyq);

  Double_t dudx = prevPoint.dudx
    + 1./6.*dr*(df1.dkxx+2.*(df2.dkxx+df3.dkxx)+df4.dkxx);
  Double_t dudy = prevPoint.dudy
    + 1./6.*dr*(df1.dkxy+2.*(df2.dkxy+df3.dkxy)+df4.dkxy);
  Double_t dudu = prevPoint.dudu
    + 1./6.*dr*(df1.dkxu+2.*(df2.dkxu+df3.dkxu)+df4.dkxu);
  Double_t dudv = prevPoint.dudv
    + 1./6.*dr*(df1.dkxv+2.*(df2.dkxv+df3.dkxv)+df4.dkxv);
  Double_t dudq = prevPoint.dudq
    + 1./6.*dr*(df1.dkxq+2.*(df2.dkxq+df3.dkxq)+df4.dkxq);

  Double_t dvdx = prevPoint.dvdx
    + 1./6.*dr*(df1.dkyx+2.*(df2.dkyx+df3.dkyx)+df4.dkyx);
  Double_t dvdy = prevPoint.dvdy
    + 1./6.*dr*(df1.dkyy+2.*(df2.dkyy+df3.dkyy)+df4.dkyy);
  Double_t dvdu = prevPoint.dvdu
    + 1./6.*dr*(df1.dkyu+2.*(df2.dkyu+df3.dkyu)+df4.dkyu);
  Double_t dvdv = prevPoint.dvdv
    + 1./6.*dr*(df1.dkyv+2.*(df2.dkyv+df3.dkyv)+df4.dkyv);
  Double_t dvdq = prevPoint.dvdq
    + 1./6.*dr*(df1.dkyq+2.*(df2.dkyq+df3.dkyq)+df4.dkyq);

  Double_t dl = (ThreeVector(x,y,z)-Z1).Mag()*StepSize/std::abs(StepSize);

#if 0
  {
    PrintHelper helper(2, std::ios::fixed);
    hddaq::cout << FUNC_NAME << ": "
		<< std::setw(9) << x
		<< std::setw(9) << y
		<< std::setw(9) << z;
    hddaq::cout.precision(4);
    hddaq::cout << std::setw(10) << u
		<< std::setw(10) << v
		<< std::setw(10) << prevPoint.r.q;
    hddaq::cout.precision(2);
    hddaq::cout << std::setw(10) << prevPoint.l+dl;
    hddaq::cout.precision(4);
    hddaq::cout << std::setw(10) << B1.z()
		<< std::endl;


    helper.set(3, std::ios::scientific);
    hddaq::cout << std::setw(10) << dxdx << " "
		<< std::setw(10) << dxdy << " "
		<< std::setw(10) << dxdu << " "
		<< std::setw(10) << dxdv << " "
		<< std::setw(10) << dxdq << std::endl;
    hddaq::cout << std::setw(10) << dydx << " "
		<< std::setw(10) << dydy << " "
		<< std::setw(10) << dydu << " "
		<< std::setw(10) << dydv << " "
		<< std::setw(10) << dydq << std::endl;
    hddaq::cout << std::setw(10) << dudx << " "
		<< std::setw(10) << dudy << " "
		<< std::setw(10) << dudu << " "
		<< std::setw(10) << dudv << " "
		<< std::setw(10) << dudq << std::endl;
    hddaq::cout << std::setw(10) << dvdx << " "
		<< std::setw(10) << dvdy << " "
		<< std::setw(10) << dvdu << " "
		<< std::setw(10) << dvdv << " "
		<< std::setw(10) << dvdq << std::endl;
  }
#endif

  return RKTrajectoryPoint(x, y, z, u, v, prevPoint.r.q,
                           dxdx, dxdy, dxdu, dxdv, dxdq,
                           dydx, dydy, dydu, dydv, dydq,
                           dudx, dudy, dudu, dudv, dudq,
                           dvdx, dvdy, dvdu, dvdv, dvdq,
                           prevPoint.l+dl);
}

//_____________________________________________________________________________
Bool_t
RK::CheckCrossing(Int_t lnum, const RKTrajectoryPoint &startPoint,
                  const RKTrajectoryPoint &endPoint,
                  RKcalcHitPoint &crossPoint)
{
  const auto geom_record = gGeom.GetRecord(lnum);
  ThreeVector posVector   = geom_record->Position();
  ThreeVector normalVector = geom_record->NormalVector();

  ThreeVector startVector = startPoint.PositionInGlobal();
  ThreeVector endVector   = endPoint.PositionInGlobal();
  // move to origin
  startVector -= posVector;
  endVector   -= posVector;

  // inner product
  Double_t ip1 = normalVector * startVector;
  Double_t ip2 = normalVector * endVector;

  // judge whether start/end points are same side
  if(ip1*ip2 > 0.) return false;

  Double_t x = (ip1*endPoint.r.x - ip2*startPoint.r.x)/(ip1-ip2);
  Double_t y = (ip1*endPoint.r.y - ip2*startPoint.r.y)/(ip1-ip2);
  Double_t z = (ip1*endPoint.r.z - ip2*startPoint.r.z)/(ip1-ip2);
  Double_t u = (ip1*endPoint.r.u - ip2*startPoint.r.u)/(ip1-ip2);
  Double_t v = (ip1*endPoint.r.v - ip2*startPoint.r.v)/(ip1-ip2);
  Double_t q = (ip1*endPoint.r.q - ip2*startPoint.r.q)/(ip1-ip2);
  Double_t l = (ip1*endPoint.l   - ip2*startPoint.l)/(ip1-ip2);

  Double_t dxdx = (ip1*endPoint.dxdx - ip2*startPoint.dxdx)/(ip1-ip2);
  Double_t dxdy = (ip1*endPoint.dxdy - ip2*startPoint.dxdy)/(ip1-ip2);
  Double_t dxdu = (ip1*endPoint.dxdu - ip2*startPoint.dxdu)/(ip1-ip2);
  Double_t dxdv = (ip1*endPoint.dxdv - ip2*startPoint.dxdv)/(ip1-ip2);
  Double_t dxdq = (ip1*endPoint.dxdq - ip2*startPoint.dxdq)/(ip1-ip2);

  Double_t dydx = (ip1*endPoint.dydx - ip2*startPoint.dydx)/(ip1-ip2);
  Double_t dydy = (ip1*endPoint.dydy - ip2*startPoint.dydy)/(ip1-ip2);
  Double_t dydu = (ip1*endPoint.dydu - ip2*startPoint.dydu)/(ip1-ip2);
  Double_t dydv = (ip1*endPoint.dydv - ip2*startPoint.dydv)/(ip1-ip2);
  Double_t dydq = (ip1*endPoint.dydq - ip2*startPoint.dydq)/(ip1-ip2);

  Double_t dudx = (ip1*endPoint.dudx - ip2*startPoint.dudx)/(ip1-ip2);
  Double_t dudy = (ip1*endPoint.dudy - ip2*startPoint.dudy)/(ip1-ip2);
  Double_t dudu = (ip1*endPoint.dudu - ip2*startPoint.dudu)/(ip1-ip2);
  Double_t dudv = (ip1*endPoint.dudv - ip2*startPoint.dudv)/(ip1-ip2);
  Double_t dudq = (ip1*endPoint.dudq - ip2*startPoint.dudq)/(ip1-ip2);

  Double_t dvdx = (ip1*endPoint.dvdx - ip2*startPoint.dvdx)/(ip1-ip2);
  Double_t dvdy = (ip1*endPoint.dvdy - ip2*startPoint.dvdy)/(ip1-ip2);
  Double_t dvdu = (ip1*endPoint.dvdu - ip2*startPoint.dvdu)/(ip1-ip2);
  Double_t dvdv = (ip1*endPoint.dvdv - ip2*startPoint.dvdv)/(ip1-ip2);
  Double_t dvdq = (ip1*endPoint.dvdq - ip2*startPoint.dvdq)/(ip1-ip2);

  Double_t pz = Polarity/(std::sqrt(1.+u*u+v*v)*q);

  crossPoint.posG = ThreeVector(x, y, z);
  crossPoint.momG = ThreeVector(pz*u, pz*v, pz);

  if(lnum==IdTOF_UX || lnum==IdTOF_DX)
    crossPoint.s = crossPoint.posG.x();
  else if(lnum==IdTOF_UY || lnum==IdTOF_DY)
    crossPoint.s = crossPoint.posG.y();
  else
    crossPoint.s = gGeom.Global2LocalPos(lnum, crossPoint.posG).x();

  crossPoint.l = l;

  Double_t sx = geom_record->dsdx();
  Double_t sy = geom_record->dsdy();
  Double_t sz = geom_record->dsdz();
  Double_t ux = geom_record->dudx();
  Double_t uy = geom_record->dudy();
  Double_t uz = geom_record->dudz();

  if(uz==0.){
    crossPoint.dsdx = crossPoint.dsdy =
      crossPoint.dsdu = crossPoint.dsdv = crossPoint.dsdq = 0.;
    crossPoint.dsdxx = crossPoint.dsdxy =
      crossPoint.dsdxu = crossPoint.dsdxv = crossPoint.dsdxq = 0.;
    crossPoint.dsdyx = crossPoint.dsdyy =
      crossPoint.dsdyu = crossPoint.dsdyv = crossPoint.dsdyq = 0.;
    crossPoint.dsdux = crossPoint.dsduy =
      crossPoint.dsduu = crossPoint.dsduv = crossPoint.dsduq = 0.;
    crossPoint.dsdvx = crossPoint.dsdvy =
      crossPoint.dsdvu = crossPoint.dsdvv = crossPoint.dsdvq = 0.;
    crossPoint.dsdqx = crossPoint.dsdqy =
      crossPoint.dsdqu = crossPoint.dsdqv = crossPoint.dsdqq = 0.;
  }
  else {
    Double_t ffx = ux/uz, ffy = uy/uz;
    Double_t dzdx = -ffx*dxdx - ffy*dydx;
    Double_t dzdy = -ffx*dxdy - ffy*dydy;
    Double_t dzdu = -ffx*dxdu - ffy*dydu;
    Double_t dzdv = -ffx*dxdv - ffy*dydv;
    Double_t dzdq = -ffx*dxdq - ffy*dydq;

    crossPoint.dsdx = sx*dxdx + sy*dydx + sz*dzdx;
    crossPoint.dsdy = sx*dxdy + sy*dydy + sz*dzdy;
    crossPoint.dsdu = sx*dxdu + sy*dydu + sz*dzdu;
    crossPoint.dsdv = sx*dxdv + sy*dydv + sz*dzdv;
    crossPoint.dsdq = sx*dxdq + sy*dydq + sz*dzdq;

    crossPoint.dsdxx = sx*dzdx*dudx + sy*dzdx*dvdx;
    crossPoint.dsdxy = sx*dzdx*dudy + sy*dzdx*dvdy;
    crossPoint.dsdxu = sx*dzdx*dudu + sy*dzdx*dvdu;
    crossPoint.dsdxv = sx*dzdx*dudv + sy*dzdx*dvdv;
    crossPoint.dsdxq = sx*dzdx*dudq + sy*dzdx*dvdq;

    crossPoint.dsdyx = sx*dzdy*dudx + sy*dzdy*dvdx;
    crossPoint.dsdyy = sx*dzdy*dudy + sy*dzdy*dvdy;
    crossPoint.dsdyu = sx*dzdy*dudu + sy*dzdy*dvdu;
    crossPoint.dsdyv = sx*dzdy*dudv + sy*dzdy*dvdv;
    crossPoint.dsdyq = sx*dzdy*dudq + sy*dzdy*dvdq;

    crossPoint.dsdux = sx*dzdu*dudx + sy*dzdu*dvdx;
    crossPoint.dsduy = sx*dzdu*dudy + sy*dzdu*dvdy;
    crossPoint.dsduu = sx*dzdu*dudu + sy*dzdu*dvdu;
    crossPoint.dsduv = sx*dzdu*dudv + sy*dzdu*dvdv;
    crossPoint.dsduq = sx*dzdu*dudq + sy*dzdu*dvdq;

    crossPoint.dsdvx = sx*dzdv*dudx + sy*dzdv*dvdx;
    crossPoint.dsdvy = sx*dzdv*dudy + sy*dzdv*dvdy;
    crossPoint.dsdvu = sx*dzdv*dudu + sy*dzdv*dvdu;
    crossPoint.dsdvv = sx*dzdv*dudv + sy*dzdv*dvdv;
    crossPoint.dsdvq = sx*dzdv*dudq + sy*dzdv*dvdq;

    crossPoint.dsdqx = sx*dzdq*dudx + sy*dzdq*dvdx;
    crossPoint.dsdqy = sx*dzdq*dudy + sy*dzdq*dvdy;
    crossPoint.dsdqu = sx*dzdq*dudu + sy*dzdq*dvdu;
    crossPoint.dsdqv = sx*dzdq*dudv + sy*dzdq*dvdv;
    crossPoint.dsdqq = sx*dzdq*dudq + sy*dzdq*dvdq;
  }

  crossPoint.dxdx=dxdx; crossPoint.dxdy=dxdy; crossPoint.dxdu=dxdu;
  crossPoint.dxdv=dxdv; crossPoint.dxdq=dxdq;
  crossPoint.dydx=dydx; crossPoint.dydy=dydy; crossPoint.dydu=dydu;
  crossPoint.dydv=dydv; crossPoint.dydq=dydq;
  crossPoint.dudx=dudx; crossPoint.dudy=dudy; crossPoint.dudu=dudu;
  crossPoint.dudv=dudv; crossPoint.dudq=dudq;
  crossPoint.dvdx=dvdx; crossPoint.dvdy=dvdy; crossPoint.dvdu=dvdu;
  crossPoint.dvdv=dvdv; crossPoint.dvdq=dvdq;

#if 0
  {
    PrintHelper helper(5, std::ios::fixed);
    hddaq::cout << FUNC_NAME << ": Layer#"
		<< std::setw(3) << lnum << std::endl;
    hddaq::cout << " " << std::setw(12) << dxdx
		<< " " << std::setw(12) << dxdy
		<< " " << std::setw(12) << dxdu
		<< " " << std::setw(12) << dxdv
		<< " " << std::setw(12) << dxdq << std::endl;
    hddaq::cout << " " << std::setw(12) << dydx
		<< " " << std::setw(12) << dydy
		<< " " << std::setw(12) << dydu
		<< " " << std::setw(12) << dydv
		<< " " << std::setw(12) << dydq << std::endl;
    hddaq::cout << " " << std::setw(12) << dudx
		<< " " << std::setw(12) << dudy
		<< " " << std::setw(12) << dudu
		<< " " << std::setw(12) << dudv
		<< " " << std::setw(12) << dudq << std::endl;
    hddaq::cout << " " << std::setw(12) << dvdx
		<< " " << std::setw(12) << dvdy
		<< " " << std::setw(12) << dvdu
		<< " " << std::setw(12) << dvdv
		<< " " << std::setw(12) << dvdq << std::endl;
  }
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
RK::CheckCrossingHS(Int_t lnum, const RKTrajectoryPoint &startPoint,
		    const RKTrajectoryPoint &endPoint,
		    RKcalcHitPoint &crossPoint)
{
  const auto geom_record = gGeom.GetRecord(lnum);
  ThreeVector posVector   = geom_record->Position();
  ThreeVector normalVector = geom_record->NormalVector();

  ThreeVector startVector(startPoint.r.x,startPoint.r.y,startPoint.r.z);
  ThreeVector endVector(endPoint.r.x,endPoint.r.y,endPoint.r.z);
  // move to origin
  startVector -= posVector;
  endVector   -= posVector;

  // inner product
  Double_t ip1 = normalVector * startVector;
  Double_t ip2 = normalVector * endVector;

  // judge whether start/end points are same side
  if(ip1*ip2 > 0.) return false;

  Double_t x = (ip1*endPoint.r.x - ip2*startPoint.r.x)/(ip1-ip2);
  Double_t y = (ip1*endPoint.r.y - ip2*startPoint.r.y)/(ip1-ip2);
  Double_t z = (ip1*endPoint.r.z - ip2*startPoint.r.z)/(ip1-ip2);
  Double_t u = (ip1*endPoint.r.u - ip2*startPoint.r.u)/(ip1-ip2);
  Double_t v = (ip1*endPoint.r.v - ip2*startPoint.r.v)/(ip1-ip2);
  Double_t q = (ip1*endPoint.r.q - ip2*startPoint.r.q)/(ip1-ip2);
  Double_t l = (ip1*endPoint.l   - ip2*startPoint.l)/(ip1-ip2);

  Double_t dxdx = (ip1*endPoint.dxdx - ip2*startPoint.dxdx)/(ip1-ip2);
  Double_t dxdy = (ip1*endPoint.dxdy - ip2*startPoint.dxdy)/(ip1-ip2);
  Double_t dxdu = (ip1*endPoint.dxdu - ip2*startPoint.dxdu)/(ip1-ip2);
  Double_t dxdv = (ip1*endPoint.dxdv - ip2*startPoint.dxdv)/(ip1-ip2);
  Double_t dxdq = (ip1*endPoint.dxdq - ip2*startPoint.dxdq)/(ip1-ip2);

  Double_t dydx = (ip1*endPoint.dydx - ip2*startPoint.dydx)/(ip1-ip2);
  Double_t dydy = (ip1*endPoint.dydy - ip2*startPoint.dydy)/(ip1-ip2);
  Double_t dydu = (ip1*endPoint.dydu - ip2*startPoint.dydu)/(ip1-ip2);
  Double_t dydv = (ip1*endPoint.dydv - ip2*startPoint.dydv)/(ip1-ip2);
  Double_t dydq = (ip1*endPoint.dydq - ip2*startPoint.dydq)/(ip1-ip2);

  Double_t dudx = (ip1*endPoint.dudx - ip2*startPoint.dudx)/(ip1-ip2);
  Double_t dudy = (ip1*endPoint.dudy - ip2*startPoint.dudy)/(ip1-ip2);
  Double_t dudu = (ip1*endPoint.dudu - ip2*startPoint.dudu)/(ip1-ip2);
  Double_t dudv = (ip1*endPoint.dudv - ip2*startPoint.dudv)/(ip1-ip2);
  Double_t dudq = (ip1*endPoint.dudq - ip2*startPoint.dudq)/(ip1-ip2);

  Double_t dvdx = (ip1*endPoint.dvdx - ip2*startPoint.dvdx)/(ip1-ip2);
  Double_t dvdy = (ip1*endPoint.dvdy - ip2*startPoint.dvdy)/(ip1-ip2);
  Double_t dvdu = (ip1*endPoint.dvdu - ip2*startPoint.dvdu)/(ip1-ip2);
  Double_t dvdv = (ip1*endPoint.dvdv - ip2*startPoint.dvdv)/(ip1-ip2);
  Double_t dvdq = (ip1*endPoint.dvdq - ip2*startPoint.dvdq)/(ip1-ip2);

  Double_t pz = Polarity/(std::sqrt(1.+u*u+v*v)*q);
  crossPoint.posG = ThreeVector(x, y, z);
  crossPoint.momG = ThreeVector(pz*u, pz*v, pz);

  crossPoint.s = gGeom.Global2LocalPos(lnum, crossPoint.posG).x();

  crossPoint.l = l;

  Double_t sx = geom_record->dsdx();
  Double_t sy = geom_record->dsdy();
  Double_t sz = geom_record->dsdz();
  Double_t ux = geom_record->dudx();
  Double_t uy = geom_record->dudy();
  Double_t uz = geom_record->dudz();

  if(uz==0.){
    crossPoint.dsdx = crossPoint.dsdy =
      crossPoint.dsdu = crossPoint.dsdv = crossPoint.dsdq = 0.;
    crossPoint.dsdxx = crossPoint.dsdxy =
      crossPoint.dsdxu = crossPoint.dsdxv = crossPoint.dsdxq = 0.;
    crossPoint.dsdyx = crossPoint.dsdyy =
      crossPoint.dsdyu = crossPoint.dsdyv = crossPoint.dsdyq = 0.;
    crossPoint.dsdux = crossPoint.dsduy =
      crossPoint.dsduu = crossPoint.dsduv = crossPoint.dsduq = 0.;
    crossPoint.dsdvx = crossPoint.dsdvy =
      crossPoint.dsdvu = crossPoint.dsdvv = crossPoint.dsdvq = 0.;
    crossPoint.dsdqx = crossPoint.dsdqy =
      crossPoint.dsdqu = crossPoint.dsdqv = crossPoint.dsdqq = 0.;
  }
  else {
    Double_t ffx = ux/uz, ffy = uy/uz;
    Double_t dzdx = -ffx*dxdx - ffy*dydx;
    Double_t dzdy = -ffx*dxdy - ffy*dydy;
    Double_t dzdu = -ffx*dxdu - ffy*dydu;
    Double_t dzdv = -ffx*dxdv - ffy*dydv;
    Double_t dzdq = -ffx*dxdq - ffy*dydq;

    crossPoint.dsdx = sx*dxdx + sy*dydx + sz*dzdx;
    crossPoint.dsdy = sx*dxdy + sy*dydy + sz*dzdy;
    crossPoint.dsdu = sx*dxdu + sy*dydu + sz*dzdu;
    crossPoint.dsdv = sx*dxdv + sy*dydv + sz*dzdv;
    crossPoint.dsdq = sx*dxdq + sy*dydq + sz*dzdq;

    crossPoint.dsdxx = sx*dzdx*dudx + sy*dzdx*dvdx;
    crossPoint.dsdxy = sx*dzdx*dudy + sy*dzdx*dvdy;
    crossPoint.dsdxu = sx*dzdx*dudu + sy*dzdx*dvdu;
    crossPoint.dsdxv = sx*dzdx*dudv + sy*dzdx*dvdv;
    crossPoint.dsdxq = sx*dzdx*dudq + sy*dzdx*dvdq;

    crossPoint.dsdyx = sx*dzdy*dudx + sy*dzdy*dvdx;
    crossPoint.dsdyy = sx*dzdy*dudy + sy*dzdy*dvdy;
    crossPoint.dsdyu = sx*dzdy*dudu + sy*dzdy*dvdu;
    crossPoint.dsdyv = sx*dzdy*dudv + sy*dzdy*dvdv;
    crossPoint.dsdyq = sx*dzdy*dudq + sy*dzdy*dvdq;

    crossPoint.dsdux = sx*dzdu*dudx + sy*dzdu*dvdx;
    crossPoint.dsduy = sx*dzdu*dudy + sy*dzdu*dvdy;
    crossPoint.dsduu = sx*dzdu*dudu + sy*dzdu*dvdu;
    crossPoint.dsduv = sx*dzdu*dudv + sy*dzdu*dvdv;
    crossPoint.dsduq = sx*dzdu*dudq + sy*dzdu*dvdq;

    crossPoint.dsdvx = sx*dzdv*dudx + sy*dzdv*dvdx;
    crossPoint.dsdvy = sx*dzdv*dudy + sy*dzdv*dvdy;
    crossPoint.dsdvu = sx*dzdv*dudu + sy*dzdv*dvdu;
    crossPoint.dsdvv = sx*dzdv*dudv + sy*dzdv*dvdv;
    crossPoint.dsdvq = sx*dzdv*dudq + sy*dzdv*dvdq;

    crossPoint.dsdqx = sx*dzdq*dudx + sy*dzdq*dvdx;
    crossPoint.dsdqy = sx*dzdq*dudy + sy*dzdq*dvdy;
    crossPoint.dsdqu = sx*dzdq*dudu + sy*dzdq*dvdu;
    crossPoint.dsdqv = sx*dzdq*dudv + sy*dzdq*dvdv;
    crossPoint.dsdqq = sx*dzdq*dudq + sy*dzdq*dvdq;
  }

  crossPoint.dxdx=dxdx; crossPoint.dxdy=dxdy; crossPoint.dxdu=dxdu;
  crossPoint.dxdv=dxdv; crossPoint.dxdq=dxdq;
  crossPoint.dydx=dydx; crossPoint.dydy=dydy; crossPoint.dydu=dydu;
  crossPoint.dydv=dydv; crossPoint.dydq=dydq;
  crossPoint.dudx=dudx; crossPoint.dudy=dudy; crossPoint.dudu=dudu;
  crossPoint.dudv=dudv; crossPoint.dudq=dudq;
  crossPoint.dvdx=dvdx; crossPoint.dvdy=dvdy; crossPoint.dvdu=dvdu;
  crossPoint.dvdv=dvdv; crossPoint.dvdq=dvdq;

  return true;
}

//_____________________________________________________________________________
Bool_t
RK::CheckCrossingTPC(TPCLocalTrackHelix *tpctrack, Int_t lnum,
		     const RKTrajectoryPoint &startPoint,
		     const RKTrajectoryPoint &endPoint,
		     RKcalcHitPoint &crossPoint_x,
		     RKcalcHitPoint &crossPoint_y)
{

  Double_t qnan = TMath::QuietNaN();
  crossPoint_x.posG = ThreeVector(qnan, qnan, qnan);
  crossPoint_x.momG = ThreeVector(qnan, qnan, qnan);
  crossPoint_x.s = qnan;
  crossPoint_x.l = qnan;
  crossPoint_y.posG = ThreeVector(qnan, qnan, qnan);
  crossPoint_y.momG = ThreeVector(qnan, qnan, qnan);
  crossPoint_y.s = qnan;
  crossPoint_y.l = qnan;

  // Start and end point to find a crossing point
  ThreeVector startVector = startPoint.PositionInGlobal();
  ThreeVector endVector = endPoint.PositionInGlobal();

  // TPC Cluster
  Int_t ClusterId = lnum - PlOffsTPCHit - 1;
  TPCLTrackHit *hit = tpctrack -> GetHitInOrder(ClusterId);
  const TVector3& localhitpos = hit -> GetLocalHitPos();

  Double_t theta_pad = hit -> GetPadTheta();
  Double_t RA2_pad = TMath::ATan2(TMath::Sin(theta_pad), TMath::Cos(theta_pad))*TMath::RadToDeg();

  const auto geom_record = gGeom.GetRecord(IdHS); //for align between HS and Kurama or K1.8 line
  const Double_t tilt = geom_record -> TiltAngle();
  const Double_t RA1 = geom_record -> RotationAngle1();
  const Double_t RA2 = geom_record -> RotationAngle2() + RA2_pad;
  Double_t ct1 = TMath::Cos(RA1*TMath::DegToRad());
  Double_t st1 = TMath::Sin(RA1*TMath::DegToRad());
  Double_t ct2 = TMath::Cos(RA2*TMath::DegToRad());
  Double_t st2 = TMath::Sin(RA2*TMath::DegToRad());

  // Normal vector of the layer
  Double_t xu = ct1*st2; Double_t yu = -st1; Double_t zu = ct1*ct2;
  ThreeVector normalVector = TVector3(xu, yu, zu);

  // Global position of the TPC cluster on the Kurama frame
  ThreeVector posVector = gGeom.Local2GlobalPos(IdHS, localhitpos);

  // move to origin
  ThreeVector startVector_shifted = startVector - posVector;
  ThreeVector endVector_shifted = endVector - posVector;

  // inner product
  Double_t ip1 = normalVector * startVector_shifted;
  Double_t ip2 = normalVector * endVector_shifted;
  if(ip1*ip2 > 0.) return false;

  Double_t x = (ip1*endPoint.r.x - ip2*startPoint.r.x)/(ip1-ip2);
  Double_t y = (ip1*endPoint.r.y - ip2*startPoint.r.y)/(ip1-ip2);
  Double_t z = (ip1*endPoint.r.z - ip2*startPoint.r.z)/(ip1-ip2);
  Double_t u = (ip1*endPoint.r.u - ip2*startPoint.r.u)/(ip1-ip2);
  Double_t v = (ip1*endPoint.r.v - ip2*startPoint.r.v)/(ip1-ip2);
  Double_t q = (ip1*endPoint.r.q - ip2*startPoint.r.q)/(ip1-ip2);
  Double_t l = (ip1*endPoint.l   - ip2*startPoint.l)/(ip1-ip2);

  Double_t dxdx = (ip1*endPoint.dxdx - ip2*startPoint.dxdx)/(ip1-ip2);
  Double_t dxdy = (ip1*endPoint.dxdy - ip2*startPoint.dxdy)/(ip1-ip2);
  Double_t dxdu = (ip1*endPoint.dxdu - ip2*startPoint.dxdu)/(ip1-ip2);
  Double_t dxdv = (ip1*endPoint.dxdv - ip2*startPoint.dxdv)/(ip1-ip2);
  Double_t dxdq = (ip1*endPoint.dxdq - ip2*startPoint.dxdq)/(ip1-ip2);

  Double_t dydx = (ip1*endPoint.dydx - ip2*startPoint.dydx)/(ip1-ip2);
  Double_t dydy = (ip1*endPoint.dydy - ip2*startPoint.dydy)/(ip1-ip2);
  Double_t dydu = (ip1*endPoint.dydu - ip2*startPoint.dydu)/(ip1-ip2);
  Double_t dydv = (ip1*endPoint.dydv - ip2*startPoint.dydv)/(ip1-ip2);
  Double_t dydq = (ip1*endPoint.dydq - ip2*startPoint.dydq)/(ip1-ip2);

  Double_t dudx = (ip1*endPoint.dudx - ip2*startPoint.dudx)/(ip1-ip2);
  Double_t dudy = (ip1*endPoint.dudy - ip2*startPoint.dudy)/(ip1-ip2);
  Double_t dudu = (ip1*endPoint.dudu - ip2*startPoint.dudu)/(ip1-ip2);
  Double_t dudv = (ip1*endPoint.dudv - ip2*startPoint.dudv)/(ip1-ip2);
  Double_t dudq = (ip1*endPoint.dudq - ip2*startPoint.dudq)/(ip1-ip2);

  Double_t dvdx = (ip1*endPoint.dvdx - ip2*startPoint.dvdx)/(ip1-ip2);
  Double_t dvdy = (ip1*endPoint.dvdy - ip2*startPoint.dvdy)/(ip1-ip2);
  Double_t dvdu = (ip1*endPoint.dvdu - ip2*startPoint.dvdu)/(ip1-ip2);
  Double_t dvdv = (ip1*endPoint.dvdv - ip2*startPoint.dvdv)/(ip1-ip2);
  Double_t dvdq = (ip1*endPoint.dvdq - ip2*startPoint.dvdq)/(ip1-ip2);

  Double_t pz = Polarity/(std::sqrt(1.+u*u+v*v)*q);

  //Horizontal component
  crossPoint_x.posG = ThreeVector(x, y, z);
  crossPoint_x.momG = ThreeVector(pz*u, pz*v, pz);
  crossPoint_x.s = gGeom.Global2LocalPos(posVector, tilt, RA1, RA2, crossPoint_x.posG).x();
  crossPoint_x.l = l;

  Double_t ct0 = TMath::Cos(tilt*TMath::DegToRad());
  Double_t st0 = TMath::Sin(tilt*TMath::DegToRad());

  //Unit vector of local position
  Double_t sx = ct0*ct2+st0*st1*st2; Double_t sy = st0*ct1; Double_t sz = -ct0*st2+st0*st1*ct2;
  Double_t ux = ct1*st2; Double_t uy = -st1; Double_t uz = ct1*ct2;

  if(uz==0.){
    crossPoint_x.dsdx = crossPoint_x.dsdy =
      crossPoint_x.dsdu = crossPoint_x.dsdv = crossPoint_x.dsdq = 0.;
    crossPoint_x.dsdxx = crossPoint_x.dsdxy =
      crossPoint_x.dsdxu = crossPoint_x.dsdxv = crossPoint_x.dsdxq = 0.;
    crossPoint_x.dsdyx = crossPoint_x.dsdyy =
      crossPoint_x.dsdyu = crossPoint_x.dsdyv = crossPoint_x.dsdyq = 0.;
    crossPoint_x.dsdux = crossPoint_x.dsduy =
      crossPoint_x.dsduu = crossPoint_x.dsduv = crossPoint_x.dsduq = 0.;
    crossPoint_x.dsdvx = crossPoint_x.dsdvy =
      crossPoint_x.dsdvu = crossPoint_x.dsdvv = crossPoint_x.dsdvq = 0.;
    crossPoint_x.dsdqx = crossPoint_x.dsdqy =
      crossPoint_x.dsdqu = crossPoint_x.dsdqv = crossPoint_x.dsdqq = 0.;
  }
  else {
    Double_t ffx = ux/uz, ffy = uy/uz;
    Double_t dzdx = -ffx*dxdx - ffy*dydx;
    Double_t dzdy = -ffx*dxdy - ffy*dydy;
    Double_t dzdu = -ffx*dxdu - ffy*dydu;
    Double_t dzdv = -ffx*dxdv - ffy*dydv;
    Double_t dzdq = -ffx*dxdq - ffy*dydq;

    crossPoint_x.dsdx = sx*dxdx + sy*dydx + sz*dzdx;
    crossPoint_x.dsdy = sx*dxdy + sy*dydy + sz*dzdy;
    crossPoint_x.dsdu = sx*dxdu + sy*dydu + sz*dzdu;
    crossPoint_x.dsdv = sx*dxdv + sy*dydv + sz*dzdv;
    crossPoint_x.dsdq = sx*dxdq + sy*dydq + sz*dzdq;

    crossPoint_x.dsdxx = sx*dzdx*dudx + sy*dzdx*dvdx;
    crossPoint_x.dsdxy = sx*dzdx*dudy + sy*dzdx*dvdy;
    crossPoint_x.dsdxu = sx*dzdx*dudu + sy*dzdx*dvdu;
    crossPoint_x.dsdxv = sx*dzdx*dudv + sy*dzdx*dvdv;
    crossPoint_x.dsdxq = sx*dzdx*dudq + sy*dzdx*dvdq;

    crossPoint_x.dsdyx = sx*dzdy*dudx + sy*dzdy*dvdx;
    crossPoint_x.dsdyy = sx*dzdy*dudy + sy*dzdy*dvdy;
    crossPoint_x.dsdyu = sx*dzdy*dudu + sy*dzdy*dvdu;
    crossPoint_x.dsdyv = sx*dzdy*dudv + sy*dzdy*dvdv;
    crossPoint_x.dsdyq = sx*dzdy*dudq + sy*dzdy*dvdq;

    crossPoint_x.dsdux = sx*dzdu*dudx + sy*dzdu*dvdx;
    crossPoint_x.dsduy = sx*dzdu*dudy + sy*dzdu*dvdy;
    crossPoint_x.dsduu = sx*dzdu*dudu + sy*dzdu*dvdu;
    crossPoint_x.dsduv = sx*dzdu*dudv + sy*dzdu*dvdv;
    crossPoint_x.dsduq = sx*dzdu*dudq + sy*dzdu*dvdq;

    crossPoint_x.dsdvx = sx*dzdv*dudx + sy*dzdv*dvdx;
    crossPoint_x.dsdvy = sx*dzdv*dudy + sy*dzdv*dvdy;
    crossPoint_x.dsdvu = sx*dzdv*dudu + sy*dzdv*dvdu;
    crossPoint_x.dsdvv = sx*dzdv*dudv + sy*dzdv*dvdv;
    crossPoint_x.dsdvq = sx*dzdv*dudq + sy*dzdv*dvdq;

    crossPoint_x.dsdqx = sx*dzdq*dudx + sy*dzdq*dvdx;
    crossPoint_x.dsdqy = sx*dzdq*dudy + sy*dzdq*dvdy;
    crossPoint_x.dsdqu = sx*dzdq*dudu + sy*dzdq*dvdu;
    crossPoint_x.dsdqv = sx*dzdq*dudv + sy*dzdq*dvdv;
    crossPoint_x.dsdqq = sx*dzdq*dudq + sy*dzdq*dvdq;
  }

  crossPoint_x.dxdx=dxdx; crossPoint_x.dxdy=dxdy; crossPoint_x.dxdu=dxdu;
  crossPoint_x.dxdv=dxdv; crossPoint_x.dxdq=dxdq;
  crossPoint_x.dydx=dydx; crossPoint_x.dydy=dydy; crossPoint_x.dydu=dydu;
  crossPoint_x.dydv=dydv; crossPoint_x.dydq=dydq;
  crossPoint_x.dudx=dudx; crossPoint_x.dudy=dudy; crossPoint_x.dudu=dudu;
  crossPoint_x.dudv=dudv; crossPoint_x.dudq=dudq;
  crossPoint_x.dvdx=dvdx; crossPoint_x.dvdy=dvdy; crossPoint_x.dvdu=dvdu;
  crossPoint_x.dvdv=dvdv; crossPoint_x.dvdq=dvdq;

#if 0
  {
    PrintHelper helper(5, std::ios::fixed);
    hddaq::cout << FUNC_NAME << ": Layer#"
		<< std::setw(3) << lnum << std::endl;
    hddaq::cout << " " << std::setw(12) << dxdx
		<< " " << std::setw(12) << dxdy
		<< " " << std::setw(12) << dxdu
		<< " " << std::setw(12) << dxdv
		<< " " << std::setw(12) << dxdq << std::endl;
    hddaq::cout << " " << std::setw(12) << dydx
		<< " " << std::setw(12) << dydy
		<< " " << std::setw(12) << dydu
		<< " " << std::setw(12) << dydv
		<< " " << std::setw(12) << dydq << std::endl;
    hddaq::cout << " " << std::setw(12) << dudx
		<< " " << std::setw(12) << dudy
		<< " " << std::setw(12) << dudu
		<< " " << std::setw(12) << dudv
		<< " " << std::setw(12) << dudq << std::endl;
    hddaq::cout << " " << std::setw(12) << dvdx
		<< " " << std::setw(12) << dvdy
		<< " " << std::setw(12) << dvdu
		<< " " << std::setw(12) << dvdv
		<< " " << std::setw(12) << dvdq << std::endl;
  }
#endif

  //Vertical component
  ct0 = TMath::Cos((tilt + 90.)*TMath::DegToRad());
  st0 = TMath::Sin((tilt + 90.)*TMath::DegToRad());

  crossPoint_y.posG = ThreeVector(x, y, z);
  crossPoint_y.momG = ThreeVector(pz*u, pz*v, pz);
  crossPoint_y.s = gGeom.Global2LocalPos(posVector, tilt + 90., RA1, RA2, crossPoint_y.posG).x();
  crossPoint_y.l = l;

  //Unit vector of local position
  sx = ct0*ct2+st0*st1*st2; sy = st0*ct1; sz = -ct0*st2+st0*st1*ct2;
  if(uz==0.){
    crossPoint_y.dsdx = crossPoint_y.dsdy =
      crossPoint_y.dsdu = crossPoint_y.dsdv = crossPoint_y.dsdq = 0.;
    crossPoint_y.dsdxx = crossPoint_y.dsdxy =
      crossPoint_y.dsdxu = crossPoint_y.dsdxv = crossPoint_y.dsdxq = 0.;
    crossPoint_y.dsdyx = crossPoint_y.dsdyy =
      crossPoint_y.dsdyu = crossPoint_y.dsdyv = crossPoint_y.dsdyq = 0.;
    crossPoint_y.dsdux = crossPoint_y.dsduy =
      crossPoint_y.dsduu = crossPoint_y.dsduv = crossPoint_y.dsduq = 0.;
    crossPoint_y.dsdvx = crossPoint_y.dsdvy =
      crossPoint_y.dsdvu = crossPoint_y.dsdvv = crossPoint_y.dsdvq = 0.;
    crossPoint_y.dsdqx = crossPoint_y.dsdqy =
      crossPoint_y.dsdqu = crossPoint_y.dsdqv = crossPoint_y.dsdqq = 0.;
  }
  else {
    Double_t ffx = ux/uz, ffy = uy/uz;
    Double_t dzdx = -ffx*dxdx - ffy*dydx;
    Double_t dzdy = -ffx*dxdy - ffy*dydy;
    Double_t dzdu = -ffx*dxdu - ffy*dydu;
    Double_t dzdv = -ffx*dxdv - ffy*dydv;
    Double_t dzdq = -ffx*dxdq - ffy*dydq;

    crossPoint_y.dsdx = sx*dxdx + sy*dydx + sz*dzdx;
    crossPoint_y.dsdy = sx*dxdy + sy*dydy + sz*dzdy;
    crossPoint_y.dsdu = sx*dxdu + sy*dydu + sz*dzdu;
    crossPoint_y.dsdv = sx*dxdv + sy*dydv + sz*dzdv;
    crossPoint_y.dsdq = sx*dxdq + sy*dydq + sz*dzdq;

    crossPoint_y.dsdxx = sx*dzdx*dudx + sy*dzdx*dvdx;
    crossPoint_y.dsdxy = sx*dzdx*dudy + sy*dzdx*dvdy;
    crossPoint_y.dsdxu = sx*dzdx*dudu + sy*dzdx*dvdu;
    crossPoint_y.dsdxv = sx*dzdx*dudv + sy*dzdx*dvdv;
    crossPoint_y.dsdxq = sx*dzdx*dudq + sy*dzdx*dvdq;

    crossPoint_y.dsdyx = sx*dzdy*dudx + sy*dzdy*dvdx;
    crossPoint_y.dsdyy = sx*dzdy*dudy + sy*dzdy*dvdy;
    crossPoint_y.dsdyu = sx*dzdy*dudu + sy*dzdy*dvdu;
    crossPoint_y.dsdyv = sx*dzdy*dudv + sy*dzdy*dvdv;
    crossPoint_y.dsdyq = sx*dzdy*dudq + sy*dzdy*dvdq;

    crossPoint_y.dsdux = sx*dzdu*dudx + sy*dzdu*dvdx;
    crossPoint_y.dsduy = sx*dzdu*dudy + sy*dzdu*dvdy;
    crossPoint_y.dsduu = sx*dzdu*dudu + sy*dzdu*dvdu;
    crossPoint_y.dsduv = sx*dzdu*dudv + sy*dzdu*dvdv;
    crossPoint_y.dsduq = sx*dzdu*dudq + sy*dzdu*dvdq;

    crossPoint_y.dsdvx = sx*dzdv*dudx + sy*dzdv*dvdx;
    crossPoint_y.dsdvy = sx*dzdv*dudy + sy*dzdv*dvdy;
    crossPoint_y.dsdvu = sx*dzdv*dudu + sy*dzdv*dvdu;
    crossPoint_y.dsdvv = sx*dzdv*dudv + sy*dzdv*dvdv;
    crossPoint_y.dsdvq = sx*dzdv*dudq + sy*dzdv*dvdq;

    crossPoint_y.dsdqx = sx*dzdq*dudx + sy*dzdq*dvdx;
    crossPoint_y.dsdqy = sx*dzdq*dudy + sy*dzdq*dvdy;
    crossPoint_y.dsdqu = sx*dzdq*dudu + sy*dzdq*dvdu;
    crossPoint_y.dsdqv = sx*dzdq*dudv + sy*dzdq*dvdv;
    crossPoint_y.dsdqq = sx*dzdq*dudq + sy*dzdq*dvdq;
  }

  crossPoint_y.dxdx=dxdx; crossPoint_y.dxdy=dxdy; crossPoint_y.dxdu=dxdu;
  crossPoint_y.dxdv=dxdv; crossPoint_y.dxdq=dxdq;
  crossPoint_y.dydx=dydx; crossPoint_y.dydy=dydy; crossPoint_y.dydu=dydu;
  crossPoint_y.dydv=dydv; crossPoint_y.dydq=dydq;
  crossPoint_y.dudx=dudx; crossPoint_y.dudy=dudy; crossPoint_y.dudu=dudu;
  crossPoint_y.dudv=dudv; crossPoint_y.dudq=dudq;
  crossPoint_y.dvdx=dvdx; crossPoint_y.dvdy=dvdy; crossPoint_y.dvdu=dvdu;
  crossPoint_y.dvdv=dvdv; crossPoint_y.dvdq=dvdq;

  return true;
}

//_____________________________________________________________________________
void
RK::ELossCorrection(Int_t lnum, const RKTrajectoryPoint &prevPoint,
		    RKTrajectoryPoint &nextPoint,
		    RKcalcHitPoint &crossPoint,
		    Int_t pikp)
{

  Double_t mass = TMath::QuietNaN();
  if(pikp==0) mass = pdg::PionMass();
  else if(pikp==1) mass = pdg::KaonMass();
  else if(pikp==2) mass = pdg::ProtonMass();
  else return ;

  Double_t U2D = prevPoint.r.z < nextPoint.r.z ? 1. : -1.;
  ThreeVector momvect = crossPoint.momG;
  ThreeVector unit = momvect.Unit();
  Double_t initmom = momvect.Mag();
  Double_t E = TMath::Sqrt(initmom*initmom + mass*mass);
  Double_t u = momvect.x()/momvect.z();
  Double_t v = momvect.y()/momvect.z();

  Double_t path;
  Int_t materialid;
  if(lnum==IdBH2){
    path = 0.5*TMath::Sqrt(1.+u*u+v*v);
    materialid = 3; //Scintillator
  }
  else if(lnum==IdBAC){
    path = 0.66*TMath::Sqrt(1.+u*u+v*v);
    materialid = 6; //Aerogel
  }
  else if(lnum==IdVPHTOF || lnum==IdHTOF){
    path = 1.*TMath::Sqrt(1.+u*u+v*v);
    materialid = 3; //Scintillator
  }
  else if(lnum==IdTPCGasVessel_U || lnum==IdTPCGasVessel_D ||
	  lnum==IdK18TPCGasVessel_U || lnum==IdK18TPCGasVessel_D){
    ThreeVector lpos = gGeom.Global2LocalPos(lnum, crossPoint.posG);
    if(TMath::Abs(lpos.x())<70. && TMath::Abs(lpos.y())<255.){
      path = 0.2*TMath::Sqrt(1.+u*u+v*v);
      materialid = 4; //Mylar
    }
    else if(TMath::Abs(lpos.x())<100. && TMath::Abs(lpos.y())<285.){
      path = 0.8*TMath::Sqrt(1.+u*u+v*v);
      materialid = 5; //Al
    }
    else{
      path = 0.3*TMath::Sqrt(1.+u*u+v*v);
      materialid = 5; //Al
    }
  }
  else return;

  Double_t beta = initmom/TMath::Sqrt(initmom*initmom + mass*mass);
  Double_t dedx = 0.001*Kinematics::HypTPCdEdx(materialid, 1000.*mass, beta);
  Double_t eloss = dedx*path*U2D;
  E = E - eloss;
  Double_t pcorr = E > mass ? TMath::Sqrt(E*E - mass*mass) : 0.;

  //std::cout<<"U2D "<<U2D<<" prev step "<<prevPoint.r.z<<" mm next step "<<nextPoint.r.z<<" mm"<<std::endl;
  //std::cout<<"before : pos "<<nextPoint.PositionInGlobal()<<" mom "<<nextPoint.MomentumInGlobal().Mag()<<std::endl;
  nextPoint.r.AddMomentum(pcorr - initmom);
  //std::cout<<"after : pos "<<nextPoint.PositionInGlobal()<<" mom "<<nextPoint.MomentumInGlobal().Mag()<<std::endl;

  //std::cout<<"before : pos "<<crossPoint.PositionInGlobal()<<" mom "<<crossPoint.MomentumInGlobal().Mag()<<std::endl;
  crossPoint.momG = pcorr*unit;
  //std::cout<<"after : pos "<<crossPoint.PositionInGlobal()<<" mom "<<crossPoint.MomentumInGlobal().Mag()<<std::endl;

#if 0
  std::cout<<"mass "<<mass<<" [GeV/c2] material id "<<materialid<<" path "<<path<<" [cm] lnum "<<lnum<<std::endl;
  std::cout<<"beta "<<beta<<" dedx "<<1000.*dedx<<" [MeV/cm] eloss "<<1000.*eloss<<" [MeV] ploss "<<1000.*(initmom - pcorr)<<" [MeV/c]"<<std::endl;
  std::cout<<"E init. "<<TMath::Sqrt(initmom*initmom + mass*mass)<<" Corr "<<E<<" [GeV/c] P init. "<<initmom<<" corr. "<<pcorr<<" [GeV/c]"<<std::endl;
#endif

}

//_____________________________________________________________________________
Int_t
RK::Trace(const RKCordParameter &initial, RKHitPointContainer &hitContainer)
{
  const Int_t nPlane = hitContainer.size();
  Int_t iPlane = nPlane-1;

  RKTrajectoryPoint prevPoint(initial,
                              1., 0., 0., 0., 0.,
                              0., 1., 0., 0., 0.,
                              0., 0., 1., 0., 0.,
                              0., 0., 0., 1., 0.,
                              0.0);
  Int_t MaxStep = 10000;
  static const Double_t MaxPathLength  = 10000.; // mm
  static const Double_t NormalStepSize = -10.;   // mm
  Double_t MinStepSize = 2.;     // mm
  /*for EventDisplay*/
  std::vector<TVector3> StepPoint(MaxStep);

  Int_t iStep = 0;

  while(++iStep < MaxStep){
    Double_t StepSize = gField.StepSize(prevPoint.PositionInGlobal(),
                                        NormalStepSize, MinStepSize);
    RKTrajectoryPoint nextPoint = RK::TraceOneStep(StepSize, prevPoint);

    /* for EventDisplay */
    StepPoint[iStep-1] = nextPoint.PositionInGlobal();

    while(RK::CheckCrossing(hitContainer[iPlane].first,
                            prevPoint, nextPoint,
                            hitContainer[iPlane].second)){
#if 0
      {
	PrintHelper helper(1, std::ios::fixed);

	hddaq::cout << std::flush;
        Int_t plnum = hitContainer[iPlane].first;
        const auto& chp = hitContainer[iPlane].second;
        const auto& gpos = chp.PositionInGlobal();
        const auto& gmom = chp.MomentumInGlobal();

	hddaq::cout << FUNC_NAME << " " << iPlane << " PL#="
		    << std::setw(2) << plnum  << " X="
		    << std::setw(7) << chp.PositionInLocal()
		    << " ("  << std::setw(8) << gpos.x()
		    << ","   << std::setw(8) << gpos.y()
		    << ","   << std::setw(8) << gpos.z()
		    << ")";
	helper.precision(3);
	hddaq::cout << " P=" << std::setw(8) << gmom.Mag()
		    << "  (" << std::setw(8) << gmom.x()
		    << ","   << std::setw(8) << gmom.y()
		    << ","   << std::setw(8) << gmom.z()
		    << ")"   << std::endl;

	helper.precision(2);
	hddaq::cout << "  PL=" << std::setw(10) << chp.PathLength();
	helper.set(3, std::ios::scientific);
	hddaq::cout << "  Coeff.s = "
		    << std::setw(10) << chp.coefX() << " "
		    << std::setw(10) << chp.coefY() << " "
		    << std::setw(10) << chp.coefU() << " "
		    << std::setw(10) << chp.coefV() << " "
		    << std::setw(10) << chp.coefQ() << std::endl;
      }
#endif

      --iPlane;
      if(iPlane<0) {
	if(gEvDisp.IsReady()){
	  Double_t q = hitContainer[0].second.MomentumInGlobal().z();
	  gEvDisp.DrawKuramaTrack(iStep, StepPoint, q);
	}
	return TPCRKTrack::kPassed;
      } // if(iPlane<0)
    } // while(RKcheckCrossing())

    if(nextPoint.PathLength() > MaxPathLength){
#if WARNOUT
      hddaq::cerr << FUNC_NAME << ": Exceed MaxPathLength. "
		  << " PL=" << std::dec << nextPoint.PathLength()
		  << " Step=" << std::dec << iStep
		  << " iPlane=" << std::dec << hitContainer[iPlane+1].first
		  << std::endl;
#endif
      return KuramaTrack::kExceedMaxPathLength;
    }
    prevPoint = nextPoint;
  }// while(++iStep)

#if WARNOUT
  hddaq::cerr << FUNC_NAME << ": Exceed MaxStep. "
	      << " PL=" << std::dec << prevPoint.PathLength()
	      << " Step=" << std::dec << iStep
	      << " iPlane=" << std::dec << hitContainer[iPlane+1].first
	      << std::endl;
#endif

  return KuramaTrack::kExceedMaxStep;
}

//_____________________________________________________________________________
Int_t
RK::TraceTPC(TPCLocalTrackHelix *tpctrack, const RKCordParameter &initial, RKHitPointContainer &hitContainer, Int_t pikp)
{
  const ThreeVector& HSpos = gGeom.GetGlobalPosition("HS");

  std::vector<Int_t> checkList;
  Int_t nPlane = hitContainer.size();
  Int_t nPlaneTPC = 0;
  for(Int_t i=0;i<hitContainer.size();i++){
    if(hitContainer[i].first > PlOffsTPCHit){
      if(i%2==0) checkList.push_back(0);
      nPlaneTPC++;
    }
  }

  Int_t iPlane = nPlane - 1;
  Int_t iPlaneTPC = nPlaneTPC - 1;

  RKTrajectoryPoint prevPoint(initial,
                              1., 0., 0., 0., 0.,
                              0., 1., 0., 0., 0.,
                              0., 0., 1., 0., 0.,
                              0., 0., 0., 1., 0.,
                              0.0);
  Int_t MaxStep = 10000;
  static const Double_t MaxPathLength  = 10000.; // mm
  static const Double_t NormalStepSize = -10.;   // mm
  Double_t MinStepSize = 2.;     // mm

  /*for EventDisplay*/
  std::vector<TVector3> StepPoint(MaxStep);

  Int_t iStep = 0;
  while(++iStep < MaxStep){
    Double_t StepSize = gField.StepSize(prevPoint.PositionInGlobal(),
                                        NormalStepSize, MinStepSize);
    RKTrajectoryPoint nextPoint = RK::TraceOneStep(StepSize, prevPoint);

    /* for EventDisplay */
    StepPoint[iStep-1] = nextPoint.PositionInGlobal();

    if(iPlaneTPC >= 0 && TMath::Abs(prevPoint.PositionInGlobal().z() - HSpos.z()) < 300.){
      for(Int_t tpcid=0; tpcid<checkList.size(); ++tpcid){
	if(checkList[tpcid]!=0) continue;
	if(RK::CheckCrossingTPC(tpctrack, hitContainer[2*tpcid].first,
				prevPoint, nextPoint,
				hitContainer[2*tpcid].second,
				hitContainer[2*tpcid+1].second)){
	  checkList[tpcid] = 1;
	  --iPlaneTPC;
	  --iPlaneTPC;
	}
      }
    }
    if(iPlaneTPC >= 0 && (prevPoint.PositionInGlobal().z() - HSpos.z()) < -300.){
      iPlaneTPC = -1;
    }

    if(iPlane >= nPlaneTPC){
      while(RK::CheckCrossing(hitContainer[iPlane].first,
			      prevPoint, nextPoint,
			      hitContainer[iPlane].second)){

#if dEdxCorrection
	RK::ELossCorrection(hitContainer[iPlane].first, prevPoint,
			    nextPoint, hitContainer[iPlane].second, pikp);
#endif

#if 0
	{
	  PrintHelper helper(1, std::ios::fixed);

	  hddaq::cout << std::flush;
	  Int_t plnum = hitContainer[iPlane].first;
	  const auto& chp = hitContainer[iPlane].second;
	  const auto& gpos = chp.PositionInGlobal();
	  const auto& gmom = chp.MomentumInGlobal();

	  hddaq::cout << FUNC_NAME << " " << iPlane << " PL#="
		      << std::setw(2) << plnum  << " X="
		      << std::setw(7) << chp.PositionInLocal()
		      << " ("  << std::setw(8) << gpos.x()
		      << ","   << std::setw(8) << gpos.y()
		      << ","   << std::setw(8) << gpos.z()
		      << ")";
	  helper.precision(3);
	  hddaq::cout << " P=" << std::setw(8) << gmom.Mag()
		      << "  (" << std::setw(8) << gmom.x()
		      << ","   << std::setw(8) << gmom.y()
		      << ","   << std::setw(8) << gmom.z()
		      << ")"   << std::endl;

	  helper.precision(2);
	  hddaq::cout << "  PL=" << std::setw(10) << chp.PathLength();
	  helper.set(3, std::ios::scientific);
	  hddaq::cout << "  Coeff.s = "
		      << std::setw(10) << chp.coefX() << " "
		      << std::setw(10) << chp.coefY() << " "
		      << std::setw(10) << chp.coefU() << " "
		      << std::setw(10) << chp.coefV() << " "
		      << std::setw(10) << chp.coefQ() << std::endl;
	}
#endif

	--iPlane;
	if(iPlane < nPlaneTPC) break;
      } // while(RKcheckCrossing())
    } // if(iPlane >= 0)

    if(iPlane < nPlaneTPC && iPlaneTPC < 0) {
      if(gEvDisp.IsReady()){
	Double_t q = hitContainer[0].second.MomentumInGlobal().z();
	gEvDisp.DrawKuramaTrack(iStep, StepPoint, q);
      }
      return TPCRKTrack::kPassed;
    }

    if(nextPoint.PathLength() > MaxPathLength){
#if WARNOUT
      hddaq::cerr << FUNC_NAME << ": Exceed MaxPathLength. "
		  << " PL=" << std::dec << nextPoint.PathLength()
		  << " Step=" << std::dec << iStep
		  << " iPlane=" << std::dec << hitContainer[iPlane+1].first
		  << std::endl;
#endif
      return TPCRKTrack::kExceedMaxPathLength;
    }
    prevPoint = nextPoint;
  }// while(++iStep)

#if WARNOUT
  hddaq::cerr << FUNC_NAME << ": Exceed MaxStep. "
	      << " PL=" << std::dec << prevPoint.PathLength()
	      << " Step=" << std::dec << iStep
	      << " iPlane=" << std::dec << hitContainer[iPlane+1].first
	      << std::endl;
#endif

  return TPCRKTrack::kExceedMaxStep;
}

//_____________________________________________________________________________
Int_t
RK::Extrap(const RKCordParameter &initial, RKHitPointContainer &hitContainer, Int_t pikp)
{
  const Int_t nPlane = hitContainer.size();
  Int_t iPlane = 0;

  RKTrajectoryPoint prevPoint(initial,
                              1., 0., 0., 0., 0.,
                              0., 1., 0., 0., 0.,
                              0., 0., 1., 0., 0.,
                              0., 0., 0., 1., 0.,
                              0.0);

  Int_t MaxStep = 5000;
  static const Double_t MaxPathLength  = 5000.; // mm
  static const Double_t NormalStepSize = 5.;   // mm
  Double_t MinStepSize = 2.;     // mm
  std::vector<TVector3> StepPoint(MaxStep);

  Int_t iStep = 0;
  while(++iStep < MaxStep){
    Double_t StepSize = gField.StepSize(prevPoint.PositionInGlobal(),
                                        NormalStepSize, MinStepSize);
    RKTrajectoryPoint nextPoint = RK::TraceOneStep(StepSize, prevPoint);

    StepPoint[iStep-1] = nextPoint.PositionInGlobal();

    while(RK::CheckCrossingHS(hitContainer[iPlane].first,
			      prevPoint, nextPoint,
			      hitContainer[iPlane].second)){

#if dEdxCorrection
      RK::ELossCorrection(hitContainer[iPlane].first, prevPoint,
			  nextPoint, hitContainer[iPlane].second, pikp);
#endif


      hddaq::cout << std::flush;
      const auto& chp = hitContainer[iPlane].second;
      const auto& gPos = chp.PositionInGlobal();
      const auto& gMom = chp.MomentumInGlobal();
      const auto &lPos = gGeom.Global2LocalPos(IdTarget,gPos);
      const auto &lMom = gGeom.Global2LocalDir(IdTarget,gMom);
      //Int_t plnum = hitContainer[iPlane].first;
      //std::cout << plnum << " lPos " << lPos.X() << " " << lPos.Y() << " " << lPos.Z()<<std::endl;
      //std::cout << plnum << " gPos " << gPos.X() << " " << gPos.Y() << " " << gPos.Z()<< " path length "<<nextPoint.PathLength()<<std::endl;
      ++iPlane;
      if(iPlane == nPlane) {
	if(gEvDisp.IsReady()){
	  Double_t q = hitContainer[nPlane-1].second.MomentumInGlobal().z();
	  gEvDisp.DrawHSTrack(iStep, StepPoint, q);
	}
	return HSTrack::kPassed;
      }
    } // while(RKcheckCrossing())
    if(nextPoint.PathLength() > MaxPathLength){
      //std::cout<<"path length "<<nextPoint.PathLength()<<std::endl;
      return HSTrack::kExceedMaxPathLength;
    }
    prevPoint = nextPoint;
  }// while(++iStep)

  return HSTrack::kExceedMaxStep;
}

//_____________________________________________________________________________
Int_t
RK::ExtrapTPC(TPCLocalTrackHelix *tpctrack, const RKCordParameter &initial, RKHitPointContainer &hitContainer, Int_t pikp)
{
  const ThreeVector& HSpos = gGeom.GetGlobalPosition("HS");

  std::vector<Int_t> checkList;
  Int_t nPlane = hitContainer.size();
  Int_t nPlaneTPC = 0;
  for(Int_t i=0;i<hitContainer.size();i++){
    if(hitContainer[i].first > PlOffsTPCHit){
      if(i%2) checkList.push_back(0);
      nPlaneTPC++;
    }
  }

  Int_t iPlane = nPlaneTPC;
  Int_t iPlaneTPC = 0;

  RKTrajectoryPoint prevPoint(initial,
                              1., 0., 0., 0., 0.,
                              0., 1., 0., 0., 0.,
                              0., 0., 1., 0., 0.,
                              0., 0., 0., 1., 0.,
                              0.0);
  Int_t MaxStep = 10000;
  static const Double_t MaxPathLength  = 10000.; // mm
  static const Double_t NormalStepSize = 10.;   // mm
  Double_t MinStepSize = 2.;     // mm

  /*for EventDisplay*/
  std::vector<TVector3> StepPoint(MaxStep);

  Int_t iStep = 0;
  while(++iStep < MaxStep){
    Double_t StepSize = gField.StepSize(prevPoint.PositionInGlobal(),
                                        NormalStepSize, MinStepSize);
    RKTrajectoryPoint nextPoint = RK::TraceOneStep(StepSize, prevPoint);

    /* for EventDisplay */
    StepPoint[iStep-1] = nextPoint.PositionInGlobal();

    if(iPlaneTPC < nPlaneTPC && TMath::Abs(prevPoint.PositionInGlobal().z() - HSpos.z()) < 300.){
      for(Int_t tpcid=0; tpcid<checkList.size(); ++tpcid){
	if(checkList[tpcid]!=0) continue;
	if(RK::CheckCrossingTPC(tpctrack, hitContainer[2*tpcid].first,
				prevPoint, nextPoint,
				hitContainer[2*tpcid].second,
				hitContainer[2*tpcid+1].second)){
	  checkList[tpcid] = 1;
	  ++iPlaneTPC;
	  ++iPlaneTPC;
	}
      }
    }
    if(iPlaneTPC < nPlaneTPC && (prevPoint.PositionInGlobal().z() - HSpos.z()) > 300.) iPlaneTPC = nPlaneTPC;

    if(iPlane < nPlane){
      while(RK::CheckCrossing(hitContainer[iPlane].first,
			      prevPoint, nextPoint,
			      hitContainer[iPlane].second)){

#if dEdxCorrection
	RK::ELossCorrection(hitContainer[iPlane].first, prevPoint,
			    nextPoint, hitContainer[iPlane].second, pikp);
#endif

#if 0
	{
	  PrintHelper helper(1, std::ios::fixed);

	  hddaq::cout << std::flush;
	  Int_t plnum = hitContainer[iPlane].first;
	  const auto& chp = hitContainer[iPlane].second;
	  const auto& gpos = chp.PositionInGlobal();
	  const auto& gmom = chp.MomentumInGlobal();

	  hddaq::cout << FUNC_NAME << " " << iPlane << " PL#="
		      << std::setw(2) << plnum  << " X="
		      << std::setw(7) << chp.PositionInLocal()
		      << " ("  << std::setw(8) << gpos.x()
		      << ","   << std::setw(8) << gpos.y()
		      << ","   << std::setw(8) << gpos.z()
		      << ")";
	  helper.precision(3);
	  hddaq::cout << " P=" << std::setw(8) << gmom.Mag()
		      << "  (" << std::setw(8) << gmom.x()
		      << ","   << std::setw(8) << gmom.y()
		      << ","   << std::setw(8) << gmom.z()
		      << ")"   << std::endl;

	  helper.precision(2);
	  hddaq::cout << "  PL=" << std::setw(10) << chp.PathLength();
	  helper.set(3, std::ios::scientific);
	  hddaq::cout << "  Coeff.s = "
		      << std::setw(10) << chp.coefX() << " "
		      << std::setw(10) << chp.coefY() << " "
		      << std::setw(10) << chp.coefU() << " "
		      << std::setw(10) << chp.coefV() << " "
		      << std::setw(10) << chp.coefQ() << std::endl;
	}
#endif
	++iPlane;
	if(iPlane >= nPlaneTPC) break;
      } // while(RKcheckCrossing())
    } // if(iPlane<nPlane)

    if(iPlane >= nPlane && iPlaneTPC >= nPlaneTPC) {
      if(gEvDisp.IsReady()){
	Double_t q = hitContainer[0].second.MomentumInGlobal().z();
	gEvDisp.DrawKuramaTrack(iStep, StepPoint, q);
      }
      return TPCRKTrack::kPassed;
    }

    if(nextPoint.PathLength() > MaxPathLength){
#if WARNOUT
      hddaq::cerr << FUNC_NAME << ": Exceed MaxPathLength. "
		  << " PL=" << std::dec << nextPoint.PathLength()
		  << " Step=" << std::dec << iStep
		  << " iPlane=" << std::dec << hitContainer[iPlane+1].first
		  << std::endl;
#endif
      return TPCRKTrack::kExceedMaxPathLength;
    }
    prevPoint = nextPoint;
  }// while(++iStep)

#if WARNOUT
  hddaq::cerr << FUNC_NAME << ": Exceed MaxStep. "
	      << " PL=" << std::dec << prevPoint.PathLength()
	      << " Step=" << std::dec << iStep
	      << " iPlane=" << std::dec << hitContainer[iPlane+1].first
	      << std::endl;
#endif

  return TPCRKTrack::kExceedMaxStep;
}

//_____________________________________________________________________________
Bool_t
RK::TraceToLast(RKHitPointContainer& hitContainer)
{
  Int_t nPlane = hitContainer.size();
  Int_t iPlane = nPlane-1;
  Int_t idini  = hitContainer[iPlane].first;

  // add downstream detectors from the TOF
  for(const auto& lid: gGeom.GetDetectorIDList()){
    if(idini < lid && lid <= IdRKINIT){
      hitContainer.push_back(std::make_pair(lid, RKcalcHitPoint()));
    }
  }

  nPlane = hitContainer.size();

  const RKcalcHitPoint& hpini  = hitContainer[iPlane].second;
  // const ThreeVector&    iniPos = hpini.PositionInGlobal();
  // const ThreeVector&    iniMom = hpini.MomentumInGlobal();

  RKTrajectoryPoint prevPoint(hpini);

  const Int_t    MaxStep  = 80000;
  const Double_t StepSize =  -20.;  // mm

  iPlane += 1;

  ThreeVector StepPoint[MaxStep];
  Int_t iStep = 0;
  while(++iStep < MaxStep){
    RKTrajectoryPoint nextPoint = RK::TraceOneStep(-StepSize, prevPoint);

    StepPoint[iStep-1] = nextPoint.PositionInGlobal();

    while(RK::CheckCrossing(hitContainer[iPlane].first,
                            prevPoint, nextPoint,
                            hitContainer[iPlane].second)){
      if(++iPlane>=nPlane){
	//	if(gEvDisp.IsReady()){
	//Double_t q = hitContainer[0].second.MomentumInGlobal().z();
	//gEvDisp.DrawKuramaTrack(iStep, StepPoint, q);
	//    }
	return true;
      }
    }
    prevPoint = nextPoint;
  }

  return false;
}

//_____________________________________________________________________________
bool
RK::CheckTrackTargetCrossing(const RKTrajectoryPoint &startPoint,
			     const RKTrajectoryPoint &endPoint,
			     RKcalcHitPoint &crossPoint,
			     Int_t &crossPlaneId,
			     bool &InToOut)
{

  //Target geometry (thickness, height, width)
  Double_t tTgt = TMath::Abs(gGeom.GetGlobalPosition(IdTgtVP1).z() -
			     gGeom.GetGlobalPosition(IdTgtVP2).z());
  Double_t hTgt = TMath::Abs(gGeom.GetGlobalPosition(IdTgtVP3).y() -
			     gGeom.GetGlobalPosition(IdTgtVP4).y());
  Double_t wTgt = TMath::Abs(gGeom.GetGlobalPosition(IdTgtVP5).x() -
			     gGeom.GetGlobalPosition(IdTgtVP6).x());

  //Define virtual planes corresponding to six sides of the target
  RKHitPointContainer hitContainer;
  hitContainer.push_back(std::make_pair(IdTgtVP1, RKcalcHitPoint())); //+z
  hitContainer.push_back(std::make_pair(IdTgtVP2, RKcalcHitPoint())); //-z
  hitContainer.push_back(std::make_pair(IdTgtVP3, RKcalcHitPoint())); //+y
  hitContainer.push_back(std::make_pair(IdTgtVP4, RKcalcHitPoint())); //-y
  hitContainer.push_back(std::make_pair(IdTgtVP5, RKcalcHitPoint())); //+x
  hitContainer.push_back(std::make_pair(IdTgtVP6, RKcalcHitPoint())); //-x

  Bool_t status = false;
  RKcalcHitPoint checkPoint;
  if(RK::CheckCrossing(IdTgtVP1, startPoint, endPoint, checkPoint) &&
     TMath::Abs(gGeom.Global2LocalPos(IdTgtVP1, checkPoint.PositionInGlobal()).x()) < 0.5*wTgt &&
     TMath::Abs(gGeom.Global2LocalPos(IdTgtVP1, checkPoint.PositionInGlobal()).y()) < 0.5*hTgt){
    status = true;
    crossPoint = checkPoint;
    crossPlaneId = IdTgtVP1;
    //std::cout<<"plnum "<<crossPlaneId<<" lpos "<<gGeom.Global2LocalPos(IdK18Target, checkPoint.PositionInGlobal())<<std::endl;
  }
  else if(RK::CheckCrossing(IdTgtVP2, startPoint, endPoint, checkPoint) &&
     TMath::Abs(gGeom.Global2LocalPos(IdTgtVP2, checkPoint.PositionInGlobal()).x()) < 0.5*wTgt &&
     TMath::Abs(gGeom.Global2LocalPos(IdTgtVP2, checkPoint.PositionInGlobal()).y()) < 0.5*hTgt){
    status = true;
    crossPoint = checkPoint;
    crossPlaneId = IdTgtVP2;
    //std::cout<<"plnum "<<crossPlaneId<<" lpos "<<gGeom.Global2LocalPos(IdK18Target, checkPoint.PositionInGlobal())<<std::endl;
  }
  else if(RK::CheckCrossing(IdTgtVP3, startPoint, endPoint, checkPoint) &&
	  TMath::Abs(gGeom.Global2LocalPos(IdTgtVP3, checkPoint.PositionInGlobal()).x()) < 0.5*wTgt &&
	  TMath::Abs(gGeom.Global2LocalPos(IdTgtVP3, checkPoint.PositionInGlobal()).y()) < 0.5*tTgt){
    status = true;
    crossPoint = checkPoint;
    crossPlaneId = IdTgtVP3;
    //std::cout<<"plnum "<<crossPlaneId<<" lpos "<<gGeom.Global2LocalPos(IdK18Target, checkPoint.PositionInGlobal())<<std::endl;
  }
  else if(RK::CheckCrossing(IdTgtVP4, startPoint, endPoint, checkPoint) &&
	  TMath::Abs(gGeom.Global2LocalPos(IdTgtVP4, checkPoint.PositionInGlobal()).x()) < 0.5*wTgt &&
	  TMath::Abs(gGeom.Global2LocalPos(IdTgtVP4, checkPoint.PositionInGlobal()).y()) < 0.5*tTgt){
    status = true;
    crossPoint = checkPoint;
    crossPlaneId = IdTgtVP4;
    //std::cout<<"plnum "<<crossPlaneId<<" lpos "<<gGeom.Global2LocalPos(IdK18Target, checkPoint.PositionInGlobal())<<std::endl;
  }
  else if(RK::CheckCrossing(IdTgtVP5, startPoint, endPoint, checkPoint) &&
	  TMath::Abs(gGeom.Global2LocalPos(IdTgtVP5, checkPoint.PositionInGlobal()).x()) < 0.5*tTgt &&
	  TMath::Abs(gGeom.Global2LocalPos(IdTgtVP5, checkPoint.PositionInGlobal()).y()) < 0.5*hTgt){
    status = true;
    crossPoint = checkPoint;
    crossPlaneId = IdTgtVP5;
    //std::cout<<"plnum "<<crossPlaneId<<" lpos "<<gGeom.Global2LocalPos(IdK18Target, checkPoint.PositionInGlobal())<<std::endl;
  }
  else if(RK::CheckCrossing(IdTgtVP6, startPoint, endPoint, checkPoint) &&
	  TMath::Abs(gGeom.Global2LocalPos(IdTgtVP6, checkPoint.PositionInGlobal()).x()) < 0.5*tTgt &&
	  TMath::Abs(gGeom.Global2LocalPos(IdTgtVP6, checkPoint.PositionInGlobal()).y()) < 0.5*hTgt){
    status = true;
    crossPoint = checkPoint;
    crossPlaneId = IdTgtVP6;
    //std::cout<<"plnum "<<crossPlaneId<<" lpos "<<gGeom.Global2LocalPos(IdK18Target, checkPoint.PositionInGlobal())<<std::endl;
  }

  if(status)
    InToOut = (TMath::Abs(gGeom.Global2LocalPos(IdK18Target, startPoint.PositionInGlobal()).x()) < 0.5*wTgt &&
	       TMath::Abs(gGeom.Global2LocalPos(IdK18Target, startPoint.PositionInGlobal()).y()) < 0.5*hTgt &&
	       TMath::Abs(gGeom.Global2LocalPos(IdK18Target, startPoint.PositionInGlobal()).z()) < 0.5*tTgt);

  return status;
}

//_____________________________________________________________________________
Bool_t
RK::FindVertex(const ThreeVector XtgtKm,
	       const ThreeVector PtgtKm,
	       const ThreeVector XtgtKp,
	       const ThreeVector PtgtKp,
	       ThreeVector &Vertex,
	       Double_t &closeDist,
	       Double_t &pathKm,
	       Double_t &pathKp,
	       ThreeVector &momVertexKm,
	       ThreeVector &momVertexKp,
	       RKcalcHitPoint &inPointKm,
	       RKcalcHitPoint &outPointKp)
{

  //Target geometry (thickness, height, width)
  const ThreeVector& tgtPos = gGeom.GetGlobalPosition(IdTarget);
  Double_t tTgt = TMath::Abs(gGeom.GetGlobalPosition(IdTgtVP1).z() -
			     gGeom.GetGlobalPosition(IdTgtVP2).z());
  Double_t hTgt = TMath::Abs(gGeom.GetGlobalPosition(IdTgtVP3).y() -
			     gGeom.GetGlobalPosition(IdTgtVP4).y());
  Double_t wTgt = TMath::Abs(gGeom.GetGlobalPosition(IdTgtVP5).x() -
			     gGeom.GetGlobalPosition(IdTgtVP6).x());
  //std::cout<<"target geometry "<<tTgt<<" "<<hTgt<<" "<<wTgt<<std::endl;

  const Int_t    MaxStep  = 80000;
  const Double_t ScanWindow = 1000.;  // mm
  Double_t StepSize = 1;  // mm (1st naive scanning)
  Double_t fineStepSize = StepSize*0.01; // (2nd fine scanning)
  Int_t iStep = 0;

  //K-, K+ track at the target center
  ThreeVector XtgtGlobalKm = XtgtKm + tgtPos;
  RKCordParameter centerCordKm(XtgtGlobalKm, PtgtKm);
  RKTrajectoryPoint centerPointKm(centerCordKm,
				  1., 0., 0., 0., 0.,
				  0., 1., 0., 0., 0.,
				  0., 0., 1., 0., 0.,
				  0., 0., 0., 1., 0.,
				  0.0);
  ThreeVector XtgtGlobalKp = XtgtKp + tgtPos;
  RKCordParameter centerCordKp(XtgtGlobalKp, PtgtKp);
  RKTrajectoryPoint centerPointKp(centerCordKp,
				  1., 0., 0., 0., 0.,
				  0., 1., 0., 0., 0.,
				  0., 0., 1., 0., 0.,
				  0., 0., 0., 1., 0.,
				  0.0);

  //Scanning a point where the K- track enters the target
  std::vector<RKcalcHitPoint> crossPointKm;
  std::vector<Bool_t> InToOutKm;
  //Scanning a point where the K+ track is out of the target
  std::vector<RKcalcHitPoint> crossPointKp;
  std::vector<Bool_t> InToOutKp;

  //Scanning from the target center to downstream
  closeDist = 10000.;

  ThreeVector posVertexKm; ThreeVector posVertexKp;
  RKTrajectoryPoint prevPointKm = centerPointKm;
  RKTrajectoryPoint prevPointKp = centerPointKp;
  Double_t pathVertexKm = 0; Double_t pathVertexKp = 0;
  while(++iStep < MaxStep){

    //Search the closest point(Vertex)
    TVector3 diff = prevPointKm.PositionInGlobal() - prevPointKp.PositionInGlobal();
    if(closeDist > diff.Mag()){
      closeDist = diff.Mag();
      posVertexKm = prevPointKm.PositionInGlobal();
      posVertexKp = prevPointKp.PositionInGlobal();
      momVertexKm = prevPointKm.MomentumInGlobal();
      momVertexKp = prevPointKp.MomentumInGlobal();
      pathVertexKm = prevPointKm.PathLength();
      pathVertexKp = prevPointKp.PathLength();
    }

    //Extrapolate the K- track
    ThreeVector momKm = prevPointKm.MomentumInGlobal();
    Double_t dzKm = StepSize*std::sqrt(1. + momKm.x()*momKm.x()/(momKm.z()*momKm.z()) +
				       momKm.y()*momKm.y()/(momKm.z()*momKm.z()));
    RKTrajectoryPoint nextPointKm = RK::TraceOneStep(dzKm, prevPointKm);
    //std::cout<<"forward K- pos : "<<gGeom.Global2LocalPos(IdK18Target, prevPointKm.PositionInGlobal())<<" mom "<<momKm<<" path "<<prevPointKm.PathLength()<<std::endl;
    RKcalcHitPoint crossPoint;
    Int_t crossPlaneId; bool InToOut;
    //if(CheckTrackTargetCrossing(prevPointKm, nextPointKm, crossPoint, crossPlaneId, InToOut) && !InToOut){
    if(CheckTrackTargetCrossing(prevPointKm, nextPointKm, crossPoint, crossPlaneId, InToOut)){
      crossPointKm.push_back(crossPoint);
      InToOutKm.push_back(InToOut);
    }

    //Extrapolate the K+ track
    ThreeVector momKp = prevPointKp.MomentumInGlobal();
    Double_t dzKp = StepSize*std::sqrt(1. + momKp.x()*momKp.x()/(momKp.z()*momKp.z()) +
				       momKp.y()*momKp.y()/(momKp.z()*momKp.z()));
    RKTrajectoryPoint nextPointKp = RK::TraceOneStep(dzKp, prevPointKp);
    //std::cout<<"forawrd K+ pos : "<<gGeom.Global2LocalPos(IdK18Target, prevPointKp.PositionInGlobal())<<" mom "<<momKp<<" path "<<prevPointKp.PathLength()<<std::endl;
    //if(CheckTrackTargetCrossing(prevPointKp, nextPointKp, crossPoint, crossPlaneId, InToOut) && InToOut) crossPointKp.push_back(crossPoint);
    if(CheckTrackTargetCrossing(prevPointKp, nextPointKp, crossPoint, crossPlaneId, InToOut)){
      crossPointKp.push_back(crossPoint);
      InToOutKp.push_back(InToOut);
    }

    prevPointKm = nextPointKm;
    prevPointKp = nextPointKp;
    if(gGeom.Global2LocalPos(IdK18Target, nextPointKp.PositionInGlobal()).z() > 0.5*ScanWindow) break;
  }

  //Scanning from the target center to upstream
  iStep = 0;
  prevPointKm = centerPointKm; prevPointKp = centerPointKp;
  while(++iStep < MaxStep){

    //Search the closest point(Vertex)
    TVector3 diff = prevPointKm.PositionInGlobal() - prevPointKp.PositionInGlobal();
    if(closeDist > diff.Mag()){
      closeDist = diff.Mag();
      posVertexKm = prevPointKm.PositionInGlobal();
      posVertexKp = prevPointKp.PositionInGlobal();
      momVertexKm = prevPointKm.MomentumInGlobal();
      momVertexKp = prevPointKp.MomentumInGlobal();
      pathVertexKm = prevPointKp.PathLength();
      pathVertexKp = prevPointKp.PathLength();
    }

    //Extrapolate the K- track
    ThreeVector momKm = prevPointKm.MomentumInGlobal();
    Double_t dzKm = StepSize*std::sqrt(1. + momKm.x()*momKm.x()/(momKm.z()*momKm.z()) +
				       momKm.y()*momKm.y()/(momKm.z()*momKm.z()));
    RKTrajectoryPoint nextPointKm = RK::TraceOneStep(-dzKm, prevPointKm);
    //std::cout<<"backward K- pos : "<<gGeom.Global2LocalPos(IdK18Target, prevPointKm.PositionInGlobal())<<" mom "<<momKm<<" path "<<prevPointKm.PathLength()<<std::endl;
    RKcalcHitPoint crossPoint;
    Int_t crossPlaneId; bool OutToIn;
    //if(CheckTrackTargetCrossing(prevPointKm, nextPointKm, crossPoint, crossPlaneId, InToOut) && InToOut) crossPointKm.push_back(crossPoint);
    if(CheckTrackTargetCrossing(prevPointKm, nextPointKm, crossPoint, crossPlaneId, OutToIn)){
      crossPointKm.push_back(crossPoint);
      InToOutKm.push_back(!OutToIn);
    }

    //Extrapolate the K+ track
    ThreeVector momKp = prevPointKp.MomentumInGlobal();
    Double_t dzKp = StepSize*std::sqrt(1. + momKp.x()*momKp.x()/(momKp.z()*momKp.z()) +
				       momKp.y()*momKp.y()/(momKp.z()*momKp.z()));
    RKTrajectoryPoint nextPointKp = RK::TraceOneStep(-dzKp, prevPointKp);
    //std::cout<<"backward K+ pos : "<<gGeom.Global2LocalPos(IdK18Target, prevPointKp.PositionInGlobal())<<" mom "<<momKp<<" path "<<prevPointKp.PathLength()<<std::endl;
    //if(CheckTrackTargetCrossing(prevPointKp, nextPointKp, crossPoint, crossPlaneId, InToOut) && !InToOut) crossPointKp.push_back(crossPoint);
    if(CheckTrackTargetCrossing(prevPointKp, nextPointKp, crossPoint, crossPlaneId, OutToIn)){
      crossPointKp.push_back(crossPoint);
      InToOutKp.push_back(!OutToIn);
    }

    prevPointKm = nextPointKm;
    prevPointKp = nextPointKp;
    if(gGeom.Global2LocalPos(IdK18Target, nextPointKm.PositionInGlobal()).z() < -0.5*ScanWindow) break;
  }
  //1st scanning ends
  //std::cout<<"naive scanning result vertex : "<<gGeom.Global2LocalPos(IdK18Target, 0.5*(posVertexKm+posVertexKp))<<" Km : "<<gGeom.Global2LocalPos(IdK18Target, posVertexKm)<<" Kp : "<<gGeom.Global2LocalPos(IdK18Target, posVertexKp)<<std::endl;

  //2nd scanning starts (fine steps)
  RKCordParameter naiveVertexCordKm(posVertexKm, momVertexKm);
  RKTrajectoryPoint naiveVertexPointKm(naiveVertexCordKm,
				       1., 0., 0., 0., 0.,
				       0., 1., 0., 0., 0.,
				       0., 0., 1., 0., 0.,
				       0., 0., 0., 1., 0.,
				       pathVertexKm);
  RKCordParameter naiveVertexCordKp(posVertexKp, momVertexKp);
  RKTrajectoryPoint naiveVertexPointKp(naiveVertexCordKp,
				       1., 0., 0., 0., 0.,
				       0., 1., 0., 0., 0.,
				       0., 0., 1., 0., 0.,
				       0., 0., 0., 1., 0.,
				       pathVertexKp);

  //Scanning again to downstream
  iStep = 0;
  prevPointKm = naiveVertexPointKm; prevPointKp = naiveVertexPointKp;
  while(++iStep < MaxStep){
    //Search the closest point(Vertex)
    TVector3 diff = prevPointKm.PositionInGlobal() - prevPointKp.PositionInGlobal();
    if(closeDist > diff.Mag()){
      closeDist = diff.Mag();
      posVertexKm = prevPointKm.PositionInGlobal();
      posVertexKp = prevPointKp.PositionInGlobal();
      momVertexKm = prevPointKm.MomentumInGlobal();
      momVertexKp = prevPointKp.MomentumInGlobal();
      pathVertexKm = prevPointKm.PathLength();
      pathVertexKp = prevPointKp.PathLength();
    }

    //Extrapolate the K- track
    ThreeVector momKm = prevPointKm.MomentumInGlobal();
    Double_t dzKm = fineStepSize*std::sqrt(1. + momKm.x()*momKm.x()/(momKm.z()*momKm.z()) +
					   momKm.y()*momKm.y()/(momKm.z()*momKm.z()));
    RKTrajectoryPoint nextPointKm = RK::TraceOneStep(dzKm, prevPointKm);
    //std::cout<<"fine forward K- pos : "<<gGeom.Global2LocalPos(IdK18Target, prevPointKm.PositionInGlobal())<<" mom "<<momKm<<" path "<<prevPointKm.PathLength()<<" dzkm "<<dzKm<<std::endl;

    //Extrapolate the K+ track
    ThreeVector momKp = prevPointKp.MomentumInGlobal();
    Double_t dzKp = fineStepSize*std::sqrt(1. + momKp.x()*momKp.x()/(momKp.z()*momKp.z()) +
					   momKp.y()*momKp.y()/(momKp.z()*momKp.z()));
    RKTrajectoryPoint nextPointKp = RK::TraceOneStep(dzKp, prevPointKp);
    //std::cout<<"fine forawrd K+ pos : "<<gGeom.Global2LocalPos(IdK18Target, prevPointKp.PositionInGlobal())<<" mom "<<momKp<<" path "<<prevPointKp.PathLength()<<std::endl;

    prevPointKm = nextPointKm;
    prevPointKp = nextPointKp;
    if(nextPointKp.PositionInGlobal().z() - naiveVertexPointKp.PositionInGlobal().z() >= StepSize) break;
  }

  //Scanning again to upstream
  iStep = 0;
  prevPointKm = naiveVertexPointKm; prevPointKp = naiveVertexPointKp;
  while(++iStep < MaxStep){

    //Search the closest point(Vertex)
    TVector3 diff = prevPointKm.PositionInGlobal() - prevPointKp.PositionInGlobal();
    if(closeDist > diff.Mag()){
      closeDist = diff.Mag();
      posVertexKm = prevPointKm.PositionInGlobal();
      posVertexKp = prevPointKp.PositionInGlobal();
      momVertexKm = prevPointKm.MomentumInGlobal();
      momVertexKp = prevPointKp.MomentumInGlobal();
      pathVertexKm = prevPointKp.PathLength();
      pathVertexKp = prevPointKp.PathLength();
    }

    //Extrapolate the K- track
    ThreeVector momKm = prevPointKm.MomentumInGlobal();
    Double_t dzKm = fineStepSize*std::sqrt(1. + momKm.x()*momKm.x()/(momKm.z()*momKm.z()) +
					   momKm.y()*momKm.y()/(momKm.z()*momKm.z()));
    RKTrajectoryPoint nextPointKm = RK::TraceOneStep(-dzKm, prevPointKm);
    //std::cout<<"fine forward K- pos : "<<gGeom.Global2LocalPos(IdK18Target, prevPointKm.PositionInGlobal())<<" mom "<<momKm<<" path "<<prevPointKm.PathLength()<<" dzkm "<<dzKm<<std::endl;

    //Extrapolate the K+ track
    ThreeVector momKp = prevPointKp.MomentumInGlobal();
    Double_t dzKp = fineStepSize*std::sqrt(1. + momKp.x()*momKp.x()/(momKp.z()*momKp.z()) +
					   momKp.y()*momKp.y()/(momKp.z()*momKp.z()));
    RKTrajectoryPoint nextPointKp = RK::TraceOneStep(-dzKp, prevPointKp);
    //std::cout<<"fine forawrd K+ pos : "<<gGeom.Global2LocalPos(IdK18Target, prevPointKp.PositionInGlobal())<<" mom "<<momKp<<" path "<<prevPointKp.PathLength()<<std::endl;

    prevPointKm = nextPointKm;
    prevPointKp = nextPointKp;
    if(nextPointKp.PositionInGlobal().z() - naiveVertexPointKp.PositionInGlobal().z() <= -StepSize) break;
  }
  //2nd scanning ends
  //std::cout<<"fine scanning result vertex : "<<gGeom.Global2LocalPos(IdK18Target, 0.5*(posVertexKm+posVertexKp))<<" Km : "<<gGeom.Global2LocalPos(IdK18Target, posVertexKm)<<" Kp : "<<gGeom.Global2LocalPos(IdK18Target, posVertexKp)<<std::endl;

  ThreeVector VertexGlobal;
  VertexGlobal = 0.5*(posVertexKm + posVertexKp);
  Vertex = gGeom.Global2LocalPos(IdK18Target, VertexGlobal);
  //Check whether the vertex in the target
  pathKm = 0.; pathKp = 0.;
  if(TMath::Abs(Vertex.x()) > 0.5*wTgt || TMath::Abs(Vertex.y()) > 0.5*hTgt || TMath::Abs(Vertex.z()) > 0.5*tTgt ||
     crossPointKm.size()!=2 || crossPointKp.size()!=2){
    if(Vertex.z() > 0.5*tTgt){
      if(crossPointKm.size()==2) pathKm = TMath::Abs(crossPointKm[0].PathLength()-crossPointKm[1].PathLength());
    }
    else if(Vertex.z() < -0.5*tTgt){
      if(crossPointKp.size()==2) pathKp = TMath::Abs(crossPointKp[0].PathLength()-crossPointKp[1].PathLength());
    }
    else{
      if(crossPointKm.size()==2) pathKm = TMath::Abs(crossPointKm[0].PathLength()-crossPointKm[1].PathLength());
      if(crossPointKp.size()==2) pathKp = TMath::Abs(crossPointKp[0].PathLength()-crossPointKp[1].PathLength());
    }
    //std::cout<<"crosspoint "<<crossPointKm.size()<<" "<<crossPointKm.size()<<std::endl;
    return false;
  }
  else{
    //std::cout<<"else "<<std::endl;
    //K-, K+ target crossing points
    Int_t km_in = 0; Int_t kp_out = 0;
    //std::cout<<"km size "<<crossPointKm.size()<<" "<<InToOutKm.size()<<std::endl;
    for(Int_t i=0; i<crossPointKm.size(); ++i){
      if(!InToOutKm[i] && VertexGlobal.z() > crossPointKm[i].PositionInGlobal().z()){
	inPointKm = crossPointKm[i];
	pathKm = TMath::Abs(crossPointKm[i].PathLength()-pathVertexKm);
	//std::cout<<"K- path vertex "<<pathVertexKm<<" in "<<crossPointKm[i].PathLength()<<" "<<crossPointKm[i].PathLength()-pathVertexKm<<std::endl;
	km_in++;
      }
    }
    //std::cout<<"kp size "<<crossPointKp.size()<<" "<<InToOutKp.size()<<std::endl;
    for(Int_t i=0; i<crossPointKp.size(); ++i){
      if(InToOutKm[i] && VertexGlobal.z() < crossPointKp[i].PositionInGlobal().z()){
	outPointKp = crossPointKp[i];
	pathKp = TMath::Abs(crossPointKp[i].PathLength()-pathVertexKp);
	//std::cout<<"K+ path vertex "<<pathVertexKp<<" in "<<crossPointKp[i].PathLength()<<" "<<crossPointKp[i].PathLength()-pathVertexKp<<std::endl;
	kp_out++;
      }
    }
    if(km_in!=1 || kp_out!=1) return false;
  }
  return true;

}

//_____________________________________________________________________________
//For Kurama Tracking
RKHitPointContainer
RK::MakeHPContainer()
{
  static const auto& IdList = gGeom.GetDetectorIDList();
  RKHitPointContainer container;
  // for(auto& id: IdList){
  //   if(id <= IdTOF_DY){
  //     container.push_back(std::make_pair(id, RKcalcHitPoint()));
  //   }
  // }

  // /*** From Upstream ***/
  container.push_back(std::make_pair(IdTarget, RKcalcHitPoint()));
  for(int i=0; i<NumOfLayersVPTPC; ++i){
    container.push_back(std::make_pair(i+PlOffsVPTPC+1, RKcalcHitPoint()));
  }
  container.push_back(std::make_pair(IdVPHTOF, RKcalcHitPoint()));

  for(Int_t i=0; i<NumOfLayersSdcIn; ++i){
    container.push_back(std::make_pair(i+PlOffsSdcIn+1, RKcalcHitPoint()));
  }
  for(Int_t i=0; i<NumOfLayersVP; ++i){
    container.push_back(std::make_pair(i+PlOffsVP+1, RKcalcHitPoint()));
  }
  for(Int_t i=0; i<NumOfLayersSdcOut; ++i){
    container.push_back(std::make_pair(i+PlOffsSdcOut+1, RKcalcHitPoint()));
  }
  container.push_back(std::make_pair(IdTOF_UX, RKcalcHitPoint()));
  container.push_back(std::make_pair(IdTOF_UY, RKcalcHitPoint()));
  container.push_back(std::make_pair(IdTOF_DX, RKcalcHitPoint()));
  container.push_back(std::make_pair(IdTOF_DY, RKcalcHitPoint()));

  return container;
}

//_____________________________________________________________________________
//For TPC RK tracking
RKHitPointContainer
RK::MakeHPContainer(std::vector<Int_t> lnum)
{
  //static const auto& IdList = gGeom.GetDetectorIDList();
  RKHitPointContainer container;
  // /*** From Upstream ***/
  for(int i=0; i<lnum.size(); ++i){
    container.push_back(std::make_pair(lnum[i], RKcalcHitPoint()));
  }

  return container;
}

//_____________________________________________________________________________
//For K1.8 beam track extrapolation to the TPC(HTOF)
RKHitPointContainer
RK::MakeHSHPContainer()
{
  static const auto& IdList = gGeom.GetDetectorIDList();
  RKHitPointContainer cont;

  // /*** From Upstream ***/
  for(int i=0; i< 6; i++){
    cont.push_back(std::make_pair(i+PlOffsBcOut+1, RKcalcHitPoint()));
  }
  cont.push_back(std::make_pair(IdBAC, RKcalcHitPoint()));
  for(int i=6; i< NumOfLayersBcOut; i++){
    cont.push_back(std::make_pair(i+PlOffsBcOut+1, RKcalcHitPoint()));
  }
  cont.push_back(std::make_pair(IdBH2,RKcalcHitPoint()));
  cont.push_back(std::make_pair(IdK18TPCGasVessel_U, RKcalcHitPoint()));
  for(int i=0; i<NumOfLayersVPHS; i++){
    cont.push_back(std::make_pair(i+PlOffsVPHS+1, RKcalcHitPoint()));
  }
  cont.push_back(std::make_pair(IdK18Target,RKcalcHitPoint()));
  cont.push_back(std::make_pair(IdK18TPCGasVessel_D, RKcalcHitPoint()));
  cont.push_back(std::make_pair(IdHTOF,RKcalcHitPoint()));

  return cont;
}

//_____________________________________________________________________________
const RKcalcHitPoint&
RKHitPointContainer::HitPointOfLayer(int lnum) const
{
  for(auto itr=this->begin(), end=this->end();
      itr!=end; ++itr){
    if(lnum == itr->first) return itr->second;
  }
  throw Exception(FUNC_NAME + Form(" no record #%d", lnum));
}

//_____________________________________________________________________________
RKcalcHitPoint&
RKHitPointContainer::HitPointOfLayer(int lnum)
{
  for(auto itr=this->begin(), end=this->end();
      itr!=end; ++itr){
    if(lnum == itr->first) return itr->second;
  }
  throw Exception(FUNC_NAME + Form(" no record #%d", lnum));
}

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
#include "DCGeomMan.hh"
#include "DCGeomRecord.hh"
#include "EventDisplay.hh"
#include "Exception.hh"
#include "FieldMan.hh"
#include "FuncName.hh"
#include "KuramaTrack.hh"
#include "PrintHelper.hh"

namespace
{
const auto& gGeom   = DCGeomMan::GetInstance();
auto&       gEvDisp = EventDisplay::GetInstance();
const auto& gField  = FieldMan::GetInstance();
// const Int_t& IdTOF    = gGeom.DetectorId("TOF");
const Int_t& IdTOF_UX = gGeom.DetectorId("TOF-UX");
const Int_t& IdTOF_UY = gGeom.DetectorId("TOF-UY");
const Int_t& IdTOF_DX = gGeom.DetectorId("TOF-DX");
const Int_t& IdTOF_DY = gGeom.DetectorId("TOF-DY");
const Int_t& IdTarget = gGeom.DetectorId("Target");
const Int_t& IdHTOF = gGeom.DetectorId("HTOF");
const Int_t& IdRKINIT = gGeom.DetectorId("RKINIT");

const Int_t& IdBH2    = gGeom.DetectorId("BH2");
const Int_t& IdK18Target = gGeom.DetectorId("K18Target");

const Double_t CHLB     = 2.99792458E-4;
const Double_t Polarity = 1.;
}

#define WARNOUT 0
#define ExactFFTreat 1

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


//____________________________________________________________________________
RKTrajectoryPoint
RK::PropagateOnce(Double_t StepSize, const RKTrajectoryPoint &prevPoint)
{
  Double_t pre_x = prevPoint.r.x;
  Double_t pre_y = prevPoint.r.y;
  Double_t pre_z = prevPoint.r.z;
  Double_t pre_u = prevPoint.r.u;
  Double_t pre_v = prevPoint.r.v;
  Double_t pre_q = prevPoint.r.q;
  Double_t dr    = StepSize/std::sqrt(1.+pre_u*pre_u+pre_v*pre_v);

  //  std::cout << pre_x << "\t" << pre_y << "\t" << pre_z
  //	    << pre_u << "\t" << pre_v << "\t" << pre_q << std::endl;
  //  std::cout << StepSize << std::endl;

  ThreeVector Z1 = prevPoint.PositionInGlobal();

  ThreeVector B1 = gField.GetField(Z1);
  //  std::cout << B1.Y() << std::endl;

  ThreeVector dBdX1 = gField.GetdBdX(Z1);
  ThreeVector dBdY1 = gField.GetdBdY(Z1);
  RKFieldIntegral f1 =
    RK::CalcFieldIntegral(pre_u, pre_v, pre_q,
                          B1, dBdX1, dBdY1);
  RKDeltaFieldIntegral df1 =
    RK::CalcDeltaFieldIntegral(prevPoint, f1);

  ThreeVector Z2 = Z1 +
    ThreeVector(0.5*dr,
                0.5*dr*pre_u + 0.125*dr*dr*f1.kx,
                0.5*dr*pre_v + 0.125*dr*dr*f1.ky);
  ThreeVector B2 = gField.GetField(Z2);
  ThreeVector dBdX2 = gField.GetdBdX(Z2);
  ThreeVector dBdY2 = gField.GetdBdY(Z2);
  RKFieldIntegral f2 =
    RK::CalcFieldIntegral(pre_u + 0.5*dr*f1.kx,
                          pre_v + 0.5*dr*f1.ky,
                          pre_q, B2, dBdX2, dBdY2);
  RKDeltaFieldIntegral df2 =
    RK::CalcDeltaFieldIntegral(prevPoint, f2, df1, df1, 0.5*dr);

  RKFieldIntegral f3 =
    RK::CalcFieldIntegral(pre_u + 0.5*dr*f2.kx,
                          pre_v + 0.5*dr*f2.ky,
                          pre_q, B2, dBdX2, dBdY2);
  RKDeltaFieldIntegral df3 =
    RK::CalcDeltaFieldIntegral(prevPoint, f3, df2, df1, 0.5*dr);

  ThreeVector Z4 = Z1 +
    ThreeVector(dr,
                dr*pre_u + 0.5*dr*dr*f3.kx,
                dr*pre_v + 0.5*dr*dr*f3.ky);
  ThreeVector B4 = gField.GetField(Z4);

  ThreeVector dBdX4 = gField.GetdBdX(Z4);
  ThreeVector dBdY4 = gField.GetdBdY(Z4);
  RKFieldIntegral f4 =
    RK::CalcFieldIntegral(pre_u + dr*f3.kx,
                          pre_v + dr*f3.ky,
                          pre_q, B4, dBdX4, dBdY4);
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
  ThreeVector nomalVector = geom_record->NormalVector();

  ThreeVector startVector = startPoint.PositionInGlobal();
  ThreeVector endVector   = endPoint.PositionInGlobal();
  // move to origin
  startVector -= posVector;
  endVector   -= posVector;

  // inner product
  Double_t ip1 = nomalVector * startVector;
  Double_t ip2 = nomalVector * endVector;

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
	return KuramaTrack::kPassed;
      }
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
Bool_t
RK::TraceToLast(RKHitPointContainer& hitContainer)
{
  Int_t nPlane = hitContainer.size();
  Int_t iPlane = nPlane-1;
  Int_t idini  = hitContainer[iPlane].first;

  // add downstream detectors from TOF
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

  container.push_back(std::make_pair(IdHTOF, RKcalcHitPoint()));

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

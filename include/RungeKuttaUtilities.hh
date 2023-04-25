// -*- C++ -*-

#ifndef RUNGE_KUTTA_UTILITIES_HH
#define RUNGE_KUTTA_UTILITIES_HH

#include "ThreeVector.hh"

#include <vector>
#include <utility>
#include <iosfwd>

#include <TString.h>

class RKFieldIntegral;
class RKDeltaFieldIntegral;
class RKTrajectoryPoint;
class RKcalcHitPoint;
class RKCordParameter;
class RKHitPointContainer;

//_____________________________________________________________________________
namespace RK
{
//_____________________________________________________________________________
RKFieldIntegral
CalcFieldIntegral(Double_t U, Double_t V, Double_t Q, const ThreeVector &B);
//_____________________________________________________________________________
RKFieldIntegral
CalcFieldIntegral(Double_t U, Double_t V, Double_t Q,
                  const ThreeVector &B,
                  const ThreeVector &dBdX, const ThreeVector &dBdY);
//_____________________________________________________________________________
RKDeltaFieldIntegral
CalcDeltaFieldIntegral(const RKTrajectoryPoint &prevPoint,
                       const RKFieldIntegral &intg);
//_____________________________________________________________________________
RKDeltaFieldIntegral
CalcDeltaFieldIntegral(const RKTrajectoryPoint &prevPoint,
                       const RKFieldIntegral &intg,
                       const RKDeltaFieldIntegral &dIntg1,
                       const RKDeltaFieldIntegral &dIntg2, Double_t StepSize);
//_____________________________________________________________________________
bool
CheckCrossing(Int_t lnum, const RKTrajectoryPoint &startPoint,
              const RKTrajectoryPoint &endPoint, RKcalcHitPoint &crossPoint);
//_____________________________________________________________________________
Int_t
Trace(const RKCordParameter &initial, RKHitPointContainer &hitContainer);
//_____________________________________________________________________________
RKTrajectoryPoint
TraceOneStep(Double_t StepSize, const RKTrajectoryPoint &prevPoint);
//_____________________________________________________________________________
RKTrajectoryPoint
PropagateOnce(Double_t StepSize, const RKTrajectoryPoint &prevPoint);
//_____________________________________________________________________________
bool
TraceToLast(RKHitPointContainer &hitContainer);
//_____________________________________________________________________________
RKHitPointContainer
MakeHPContainer();

inline TString
ClassName() { static TString s_name("RK"); return s_name; }
}

//_____________________________________________________________________________
class RKFieldIntegral
{
public:
  RKFieldIntegral(Double_t Kx, Double_t Ky,
                  Double_t Axu, Double_t Axv, Double_t Ayu, Double_t Ayv,
                  Double_t Cxx=0., Double_t Cxy=0.,
                  Double_t Cyx=0., Double_t Cyy=0.)
    : kx(Kx), ky(Ky), axu(Axu), axv(Axv), ayu(Ayu), ayv(Ayv),
      cxx(Cxx), cxy(Cxy), cyx(Cyx), cyy(Cyy)
    {
    }

private:
  Double_t kx, ky;
  Double_t axu, axv, ayu, ayv;
  Double_t cxx, cxy, cyx, cyy;

public:
  void Print(std::ostream &ost) const;

  friend RKTrajectoryPoint RK::TraceOneStep(Double_t, const RKTrajectoryPoint &);
  friend RKTrajectoryPoint RK::PropagateOnce(Double_t, const RKTrajectoryPoint &);
  friend RKDeltaFieldIntegral
  RK::CalcDeltaFieldIntegral(const RKTrajectoryPoint &,
                             const RKFieldIntegral &,
                             const RKDeltaFieldIntegral &,
                             const RKDeltaFieldIntegral &, Double_t);
  friend RKDeltaFieldIntegral
  RK::CalcDeltaFieldIntegral(const RKTrajectoryPoint &,
                             const RKFieldIntegral &);
};

//_____________________________________________________________________________
class RKDeltaFieldIntegral
{
public:
  RKDeltaFieldIntegral(Double_t dKxx, Double_t dKxy, Double_t dKxu,
                       Double_t dKxv, Double_t dKxq,
                       Double_t dKyx, Double_t dKyy, Double_t dKyu,
                       Double_t dKyv, Double_t dKyq)
    : dkxx(dKxx), dkxy(dKxy), dkxu(dKxu), dkxv(dKxv), dkxq(dKxq),
      dkyx(dKyx), dkyy(dKyy), dkyu(dKyu), dkyv(dKyv), dkyq(dKyq)
    {}

private:
  Double_t dkxx, dkxy, dkxu, dkxv, dkxq;
  Double_t dkyx, dkyy, dkyu, dkyv, dkyq;
public:
  void Print(std::ostream &ost) const;
  friend RKTrajectoryPoint RK::TraceOneStep(Double_t, const RKTrajectoryPoint &);
  friend RKTrajectoryPoint RK::PropagateOnce(Double_t, const RKTrajectoryPoint &);
  friend RKDeltaFieldIntegral
  RK::CalcDeltaFieldIntegral(const RKTrajectoryPoint &,
                             const RKFieldIntegral &,
                             const RKDeltaFieldIntegral &,
                             const RKDeltaFieldIntegral &, Double_t);
  friend RKDeltaFieldIntegral
  RK::CalcDeltaFieldIntegral(const RKTrajectoryPoint &,
                             const RKFieldIntegral &);
};

//_____________________________________________________________________________
class RKCordParameter
{
public:
  RKCordParameter()
    : x(0.), y(0.), z(0.), u(0.), v(0.), q(0.)
    {}
  RKCordParameter(Double_t X, Double_t Y, Double_t Z,
                  Double_t U, Double_t V, Double_t Q)
    : x(X), y(Y), z(Z), u(U), v(V), q(Q)
    {}

  RKCordParameter(const ThreeVector &pos,
                  Double_t U, Double_t V, Double_t Q)
    : x(pos.x()), y(pos.y()), z(pos.z()),
      u(U), v(V), q(Q)
    {}

  RKCordParameter(const ThreeVector &pos,
                  const ThreeVector &mom);
private:
  Double_t x, y, z, u, v, q;
public:
  ThreeVector PositionInGlobal() const
    { return ThreeVector(x, y, z); }
  ThreeVector MomentumInGlobal() const;
  void Print(std::ostream &ost) const;

  Double_t X() const { return x; }
  Double_t Y() const { return y; }
  Double_t Z() const { return z; }
  Double_t U() const { return u; }
  Double_t V() const { return v; }
  Double_t Q() const { return q; }

  friend class RKTrajectoryPoint;
  friend RKTrajectoryPoint
  RK::TraceOneStep(Double_t, const RKTrajectoryPoint &);
  friend RKTrajectoryPoint
  RK::PropagateOnce(Double_t, const RKTrajectoryPoint &);
  friend RKDeltaFieldIntegral
  RK::CalcDeltaFieldIntegral(const RKTrajectoryPoint &,
                             const RKFieldIntegral &,
                             const RKDeltaFieldIntegral &,
                             const RKDeltaFieldIntegral &, Double_t);
  friend RKDeltaFieldIntegral
  RK::CalcDeltaFieldIntegral(const RKTrajectoryPoint &,
                             const RKFieldIntegral &);
  friend bool
  RK::CheckCrossing(int, const RKTrajectoryPoint &,
                    const RKTrajectoryPoint &, RKcalcHitPoint &);
};

//_____________________________________________________________________________
class RKcalcHitPoint
{
public:
  RKcalcHitPoint(){}
  RKcalcHitPoint(const ThreeVector &pos, const ThreeVector &mom,
                 Double_t S, Double_t L,
                 Double_t Dsdx,  Double_t Dsdy,  Double_t Dsdu,  Double_t Dsdv,  Double_t Dsdq,
                 Double_t Dsdxx, Double_t Dsdxy, Double_t Dsdxu, Double_t Dsdxv, Double_t Dsdxq,
                 Double_t Dsdyx, Double_t Dsdyy, Double_t Dsdyu, Double_t Dsdyv, Double_t Dsdyq,
                 Double_t Dsdux, Double_t Dsduy, Double_t Dsduu, Double_t Dsduv, Double_t Dsduq,
                 Double_t Dsdvx, Double_t Dsdvy, Double_t Dsdvu, Double_t Dsdvv, Double_t Dsdvq,
                 Double_t Dsdqx, Double_t Dsdqy, Double_t Dsdqu, Double_t Dsdqv, Double_t Dsdqq,
                 Double_t Dxdx,  Double_t Dxdy,  Double_t Dxdu,  Double_t Dxdv,  Double_t Dxdq,
                 Double_t Dydx,  Double_t Dydy,  Double_t Dydu,  Double_t Dydv,  Double_t Dydq,
                 Double_t Dudx,  Double_t Dudy,  Double_t Dudu,  Double_t Dudv,  Double_t Dudq,
                 Double_t Dvdx,  Double_t Dvdy,  Double_t Dvdu,  Double_t Dvdv,  Double_t Dvdq)
  : posG(pos), momG(mom), s(S), l(L),
    dsdx(Dsdx),   dsdy(Dsdy),   dsdu(Dsdu),   dsdv(Dsdv),   dsdq(Dsdq),
    dsdxx(Dsdxx), dsdxy(Dsdxy), dsdxu(Dsdxu), dsdxv(Dsdxv), dsdxq(Dsdxq),
    dsdyx(Dsdyx), dsdyy(Dsdyy), dsdyu(Dsdyu), dsdyv(Dsdyv), dsdyq(Dsdyq),
    dsdux(Dsdux), dsduy(Dsduy), dsduu(Dsduu), dsduv(Dsduv), dsduq(Dsduq),
    dsdvx(Dsdvx), dsdvy(Dsdvy), dsdvu(Dsdvu), dsdvv(Dsdvv), dsdvq(Dsdvq),
    dsdqx(Dsdqx), dsdqy(Dsdqy), dsdqu(Dsdqu), dsdqv(Dsdqv), dsdqq(Dsdqq),
    dxdx(Dxdx),   dxdy(Dxdy),   dxdu(Dxdu),   dxdv(Dxdv),   dxdq(Dxdq),
    dydx(Dydx),   dydy(Dydy),   dydu(Dydu),   dydv(Dydv),   dydq(Dydq),
    dudx(Dudx),   dudy(Dudy),   dudu(Dudu),   dudv(Dudv),   dudq(Dudq),
    dvdx(Dvdx),   dvdy(Dvdy),   dvdu(Dvdu),   dvdv(Dvdv),   dvdq(Dvdq)
    {}

private:
  ThreeVector posG, momG;
  Double_t s; // local X
  Double_t l;
  Double_t dsdx,  dsdy,  dsdu,  dsdv,  dsdq;
  Double_t dsdxx, dsdxy, dsdxu, dsdxv, dsdxq;
  Double_t dsdyx, dsdyy, dsdyu, dsdyv, dsdyq;
  Double_t dsdux, dsduy, dsduu, dsduv, dsduq;
  Double_t dsdvx, dsdvy, dsdvu, dsdvv, dsdvq;
  Double_t dsdqx, dsdqy, dsdqu, dsdqv, dsdqq;
  Double_t dxdx,  dxdy,  dxdu,  dxdv,  dxdq;
  Double_t dydx,  dydy,  dydu,  dydv,  dydq;
  Double_t dudx,  dudy,  dudu,  dudv,  dudq;
  Double_t dvdx,  dvdy,  dvdu,  dvdv,  dvdq;

public:
  const ThreeVector& PositionInGlobal() const { return posG; }
  const ThreeVector& MomentumInGlobal() const { return momG; }
  Double_t PositionInLocal() const { return s; }
  Double_t PathLength() const { return l; }
  Double_t coefX() const { return dsdx; }
  Double_t coefY() const { return dsdy; }
  Double_t coefU() const { return dsdu; }
  Double_t coefV() const { return dsdv; }
  Double_t coefQ() const { return dsdq; }
  Double_t coefXX() const { return dsdxx; }
  Double_t coefXY() const { return dsdxy; }
  Double_t coefXU() const { return dsdxu; }
  Double_t coefXV() const { return dsdxv; }
  Double_t coefXQ() const { return dsdxq; }
  Double_t coefYX() const { return dsdyx; }
  Double_t coefYY() const { return dsdyy; }
  Double_t coefYU() const { return dsdyu; }
  Double_t coefYV() const { return dsdyv; }
  Double_t coefYQ() const { return dsdyq; }
  Double_t coefUX() const { return dsdux; }
  Double_t coefUY() const { return dsduy; }
  Double_t coefUU() const { return dsduu; }
  Double_t coefUV() const { return dsduv; }
  Double_t coefUQ() const { return dsduq; }
  Double_t coefVX() const { return dsdvx; }
  Double_t coefVY() const { return dsdvy; }
  Double_t coefVU() const { return dsdvu; }
  Double_t coefVV() const { return dsdvv; }
  Double_t coefVQ() const { return dsdvq; }
  Double_t coefQX() const { return dsdqx; }
  Double_t coefQY() const { return dsdqy; }
  Double_t coefQU() const { return dsdqu; }
  Double_t coefQV() const { return dsdqv; }
  Double_t coefQQ() const { return dsdqq; }

  Double_t dXdX() const { return dxdx; }
  Double_t dXdY() const { return dxdy; }
  Double_t dXdU() const { return dxdu; }
  Double_t dXdV() const { return dxdv; }
  Double_t dXdQ() const { return dxdq; }
  Double_t dYdX() const { return dydx; }
  Double_t dYdY() const { return dydy; }
  Double_t dYdU() const { return dydu; }
  Double_t dYdV() const { return dydv; }
  Double_t dYdQ() const { return dydq; }
  Double_t dUdX() const { return dudx; }
  Double_t dUdY() const { return dudy; }
  Double_t dUdU() const { return dudu; }
  Double_t dUdV() const { return dudv; }
  Double_t dUdQ() const { return dudq; }
  Double_t dVdX() const { return dvdx; }
  Double_t dVdY() const { return dvdy; }
  Double_t dVdU() const { return dvdu; }
  Double_t dVdV() const { return dvdv; }
  Double_t dVdQ() const { return dvdq; }

  friend bool
  RK::CheckCrossing(Int_t, const RKTrajectoryPoint &,
                    const RKTrajectoryPoint &, RKcalcHitPoint &);
};

//_____________________________________________________________________________
class RKTrajectoryPoint
{
public:
  RKTrajectoryPoint(Double_t X, Double_t Y, Double_t Z,
                    Double_t U, Double_t V, Double_t Q,
                    Double_t Dxdx, Double_t Dxdy, Double_t Dxdu,
                    Double_t Dxdv, Double_t Dxdq,
                    Double_t Dydx, Double_t Dydy, Double_t Dydu,
                    Double_t Dydv, Double_t Dydq,
                    Double_t Dudx, Double_t Dudy, Double_t Dudu,
                    Double_t Dudv, Double_t Dudq,
                    Double_t Dvdx, Double_t Dvdy, Double_t Dvdu,
                    Double_t Dvdv, Double_t Dvdq,
                    Double_t L)
  : r(X,Y,Z,U,V,Q),
    dxdx(Dxdx), dxdy(Dxdy), dxdu(Dxdu), dxdv(Dxdv), dxdq(Dxdq),
    dydx(Dydx), dydy(Dydy), dydu(Dydu), dydv(Dydv), dydq(Dydq),
    dudx(Dudx), dudy(Dudy), dudu(Dudu), dudv(Dudv), dudq(Dudq),
    dvdx(Dvdx), dvdy(Dvdy), dvdu(Dvdu), dvdv(Dvdv), dvdq(Dvdq),
    l(L)
    {}

  RKTrajectoryPoint(const RKCordParameter &R,
                    Double_t Dxdx, Double_t Dxdy, Double_t Dxdu,
                    Double_t Dxdv, Double_t Dxdq,
                    Double_t Dydx, Double_t Dydy, Double_t Dydu,
                    Double_t Dydv, Double_t Dydq,
                    Double_t Dudx, Double_t Dudy, Double_t Dudu,
                    Double_t Dudv, Double_t Dudq,
                    Double_t Dvdx, Double_t Dvdy, Double_t Dvdu,
                    Double_t Dvdv, Double_t Dvdq,
                    Double_t L)
  : r(R),
    dxdx(Dxdx), dxdy(Dxdy), dxdu(Dxdu), dxdv(Dxdv), dxdq(Dxdq),
    dydx(Dydx), dydy(Dydy), dydu(Dydu), dydv(Dydv), dydq(Dydq),
    dudx(Dudx), dudy(Dudy), dudu(Dudu), dudv(Dudv), dudq(Dudq),
    dvdx(Dvdx), dvdy(Dvdy), dvdu(Dvdu), dvdv(Dvdv), dvdq(Dvdq),
    l(L)
    {}

  RKTrajectoryPoint(const ThreeVector &pos,
                    Double_t U, Double_t V, Double_t Q,
                    Double_t Dxdx, Double_t Dxdy, Double_t Dxdu,
                    Double_t Dxdv, Double_t Dxdq,
                    Double_t Dydx, Double_t Dydy, Double_t Dydu,
                    Double_t Dydv, Double_t Dydq,
                    Double_t Dudx, Double_t Dudy, Double_t Dudu,
                    Double_t Dudv, Double_t Dudq,
                    Double_t Dvdx, Double_t Dvdy, Double_t Dvdu,
                    Double_t Dvdv, Double_t Dvdq,
                    Double_t L)
  : r(pos,U,V,Q),
    dxdx(Dxdx), dxdy(Dxdy), dxdu(Dxdu), dxdv(Dxdv), dxdq(Dxdq),
    dydx(Dydx), dydy(Dydy), dydu(Dydu), dydv(Dydv), dydq(Dydq),
    dudx(Dudx), dudy(Dudy), dudu(Dudu), dudv(Dudv), dudq(Dudq),
    dvdx(Dvdx), dvdy(Dvdy), dvdu(Dvdu), dvdv(Dvdv), dvdq(Dvdq),
    l(L)
    {}

  RKTrajectoryPoint(const ThreeVector &pos,
                    const ThreeVector &mom,
                    Double_t Dxdx, Double_t Dxdy, Double_t Dxdu,
                    Double_t Dxdv, Double_t Dxdq,
                    Double_t Dydx, Double_t Dydy, Double_t Dydu,
                    Double_t Dydv, Double_t Dydq,
                    Double_t Dudx, Double_t Dudy, Double_t Dudu,
                    Double_t Dudv, Double_t Dudq,
                    Double_t Dvdx, Double_t Dvdy, Double_t Dvdu,
                    Double_t Dvdv, Double_t Dvdq,
                    Double_t L)
  : r(pos,mom),
    dxdx(Dxdx), dxdy(Dxdy), dxdu(Dxdu), dxdv(Dxdv), dxdq(Dxdq),
    dydx(Dydx), dydy(Dydy), dydu(Dydu), dydv(Dydv), dydq(Dydq),
    dudx(Dudx), dudy(Dudy), dudu(Dudu), dudv(Dudv), dudq(Dudq),
    dvdx(Dvdx), dvdy(Dvdy), dvdu(Dvdu), dvdv(Dvdv), dvdq(Dvdq),
    l(L)
    {}

  RKTrajectoryPoint(const RKcalcHitPoint& hp)
    : r(hp.PositionInGlobal(), hp.MomentumInGlobal()),
      dxdx(hp.dXdX()), dxdy(hp.dXdY()),
      dxdu(hp.dXdU()), dxdv(hp.dXdV()), dxdq(hp.dXdQ()),
      dydx(hp.dYdX()), dydy(hp.dYdY()),
      dydu(hp.dYdU()), dydv(hp.dYdV()), dydq(hp.dYdQ()),
      dudx(hp.dUdX()), dudy(hp.dUdY()),
      dudu(hp.dUdU()), dudv(hp.dUdV()), dudq(hp.dUdQ()),
      dvdx(hp.dVdX()), dvdy(hp.dVdY()),
      dvdu(hp.dVdU()), dvdv(hp.dVdV()), dvdq(hp.dVdQ()),
      l(hp.PathLength())
    {}

private:
  RKCordParameter r;
  Double_t dxdx, dxdy, dxdu, dxdv, dxdq;
  Double_t dydx, dydy, dydu, dydv, dydq;
  Double_t dudx, dudy, dudu, dudv, dudq;
  Double_t dvdx, dvdy, dvdu, dvdv, dvdq;
  Double_t l;

public:
  ThreeVector PositionInGlobal() const { return r.PositionInGlobal(); }
  ThreeVector MomentumInGlobal() const { return r.MomentumInGlobal(); }
  Double_t    PathLength() const { return l; }
  void        Print(std::ostream &ost) const;

  friend RKTrajectoryPoint
  RK::TraceOneStep(Double_t, const RKTrajectoryPoint &);
  friend RKTrajectoryPoint
  RK::PropagateOnce(Double_t, const RKTrajectoryPoint &);
  friend bool
  RK::CheckCrossing(Int_t, const RKTrajectoryPoint &,
                    const RKTrajectoryPoint &, RKcalcHitPoint &);
  friend RKDeltaFieldIntegral
  RK::CalcDeltaFieldIntegral(const RKTrajectoryPoint &,
                             const RKFieldIntegral &,
                             const RKDeltaFieldIntegral &,
                             const RKDeltaFieldIntegral &, Double_t);
  friend RKDeltaFieldIntegral
  RK::CalcDeltaFieldIntegral(const RKTrajectoryPoint &,
                             const RKFieldIntegral &);
};

//_____________________________________________________________________________
class RKHitPointContainer
  : public std::vector< std::pair<int,RKcalcHitPoint> >
{
public:
  static TString ClassName()
    { static TString s_name("RKHitPointContainer"); return s_name; }
  const RKcalcHitPoint& HitPointOfLayer(Int_t lnum) const;
  RKcalcHitPoint& HitPointOfLayer(Int_t lnum);

  typedef std::vector<std::pair<int,RKcalcHitPoint> >
  ::const_iterator RKHpCIterator;

  typedef std::vector<std::pair<int,RKcalcHitPoint> >
  ::iterator RKHpIterator;

};

//_____________________________________________________________________________
inline std::ostream&
operator <<(std::ostream &ost, const RKFieldIntegral &obj)
{
  obj.Print(ost);
  return ost;
}

//_____________________________________________________________________________
inline std::ostream&
operator <<(std::ostream &ost, const RKDeltaFieldIntegral &obj)
{
  obj.Print(ost);
  return ost;
}

//_____________________________________________________________________________
inline std::ostream&
operator <<(std::ostream &ost, const RKCordParameter &obj)
{
  obj.Print(ost);
  return ost;
}

//_____________________________________________________________________________
inline std::ostream&
operator <<(std::ostream &ost, const RKTrajectoryPoint &obj)
{
  obj.Print(ost);
  return ost;
}

#endif

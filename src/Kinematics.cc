// -*- C++ -*-

#include "Kinematics.hh"

#include <cmath>
#include <cstdio>

#include <TF2.h>
#include <TMath.h>

#include <std_ostream.hh>

#include "FuncName.hh"
#include "MathTools.hh"

namespace
{
const Double_t qnan = TMath::QuietNaN();
const Double_t z_offset      = 0.0;
const Double_t TARGETcenter  = 0.0;//Z position
const Double_t TARGEThw      = 15.0/2.0;
const Double_t TARGETsizeX   = 25.0/2.0;
const Double_t TARGETsizeY   = 15.0/2.0;
const Double_t TARGETsizeZ   = 15.0/2.0;
const Double_t TARGETradius  = 67.3/2.0;
const Double_t TARGETcenterX = 0.0;
const Double_t TARGETcenterY = 0.0;
}

//_____________________________________________________________________________
namespace Kinematics
{
//_____________________________________________________________________________
Double_t
MassSquare(Double_t p, Double_t path, Double_t time)
{
  if(time < 0) return qnan;
  Double_t beta = path/time/MathTools::C();
  return p*p*(1.-beta*beta)/beta/beta;
}

//_____________________________________________________________________________
Double_t
CalcTimeOfFlight(Double_t p, Double_t path, Double_t mass)
{
  return path*TMath::Sqrt(mass*mass+p*p)/p/MathTools::C();
}

//_____________________________________________________________________________
TVector3
VertexPoint(const TVector3& Xin, const TVector3& Xout,
            const TVector3& Pin, const TVector3& Pout)
{
  Double_t xi=Xin.x(), yi=Xin.y(), xo=Xout.x(), yo=Xout.y();
  Double_t ui=Pin.x()/Pin.z(), vi=Pin.y()/Pin.z();
  Double_t uo=Pout.x()/Pout.z(), vo=Pout.y()/Pout.z();

  Double_t z=((xi-xo)*(uo-ui)+(yi-yo)*(vo-vi))/
    ((uo-ui)*(uo-ui)+(vo-vi)*(vo-vi));
  Double_t x1=xi+ui*z, y1=yi+vi*z;
  Double_t x2=xo+uo*z, y2=yo+vo*z;
  Double_t x = 0.5*(x1+x2);
  Double_t y = 0.5*(y1+y2);
  if(std::isnan(x) || std::isnan(y) || std::isnan(z))
    return TVector3(qnan, qnan, qnan);

  return TVector3(x, y, z+z_offset);

}

//_____________________________________________________________________________
TVector3
VertexPointByHonly(const TVector3& Xin,
                   const TVector3& Xout,
                   const TVector3& Pin,
                   const TVector3& Pout)
{
  Double_t xi=Xin.x(), yi=Xin.y(), xo=Xout.x();//, yo=Xout.y();
  Double_t ui=Pin.x()/Pin.z(), vi=Pin.y()/Pin.z();
  Double_t uo=Pout.x()/Pout.z();//, vo=Pout.y()/Pout.z();

  Double_t z=(xi-xo)/(uo-ui);
  return TVector3(xi+ui*z, yi+vi*z, z+z_offset);
}

//_____________________________________________________________________________
TVector3
VertexPoint(const TVector3& Xin, const TVector3& Xout,
            const TVector3& Pin, const TVector3& Pout,
            Double_t& dist)
{
  Double_t xi=Xin.x(), yi=Xin.y(), xo=Xout.x(), yo=Xout.y();
  Double_t ui=Pin.x()/Pin.z(), vi=Pin.y()/Pin.z();
  Double_t uo=Pout.x()/Pout.z(), vo=Pout.y()/Pout.z();

  Double_t z=((xi-xo)*(uo-ui)+(yi-yo)*(vo-vi))/
    ((uo-ui)*(uo-ui)+(vo-vi)*(vo-vi));
  Double_t x1=xi+ui*z, y1=yi+vi*z;
  Double_t x2=xo+uo*z, y2=yo+vo*z;
  dist=TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

  return TVector3(0.5*(x1+x2), 0.5*(y1+y2), z+z_offset);
}

//_____________________________________________________________________________
TVector3 VertexPoint3D(const TVector3& Xin, const TVector3& Xout,
                       const TVector3& Pin, const TVector3& Pout,
                       Double_t& dist)
{
  // Beam Xin+Pin*s, Scat Xout+Pout*t
  Double_t s, t;
  TVector3 a  = Pin;
  TVector3 A0 = Xin;
  TVector3 b  = Pout;
  TVector3 B0 = Xout;
  TVector3 AB = A0-B0;

  s = (-(a*b)*(b*AB)+(b.Mag2())*(a*AB))/((a*b)*(a*b)-(a.Mag2())*(b.Mag2()));
  t = ((a*b)*(a*AB)-(a.Mag2())*(b*AB))/((a*b)*(a*b)-(a.Mag2())*(b.Mag2()));

  TVector3 pos1 = A0+a*s;
  TVector3 pos2 = B0+b*t;

  dist = (pos1-pos2).Mag();

  return (pos1+pos2)*0.5;
}

//_____________________________________________________________________________
TVector3
VertexPointReal(const TVector3& Xin, const TVector3& Xout,
                const TVector3& Pin, const TVector3& Pout,
                Double_t& dist)
{
  Double_t a1_2  = Pin.Mag()*Pin.Mag();
  Double_t a2_2  = Pout.Mag()*Pout.Mag();
  Double_t a1a2  = Pin.Dot(Pout);
  Double_t bunbo = a1_2*a2_2 -a1a2*a1a2;
  TVector3 b1_b2;
  b1_b2 = (Xin - Xout);

  Double_t t[2];
  t[0] = (b1_b2.Dot(Pout)*a1a2 - b1_b2.Dot(Pin)*a2_2)
    /bunbo;
  t[1] = (b1_b2.Dot(Pout)*a1_2 - b1_b2.Dot(Pin)*a1a2)
    /bunbo;

  TVector3 tempVector[2];
  tempVector[0] = t[0]*Pin + Xin;
  tempVector[1] = t[1]*Pout + Xout;

  Double_t vertex[3];
  vertex[0] = (tempVector[0].X() + tempVector[1].X())/2.;
  vertex[1] = (tempVector[0].Y() + tempVector[1].Y())/2.;
  vertex[2] = (tempVector[0].Z() + tempVector[1].Z())/2.;

  TVector3 DistanceVector;
  DistanceVector  = tempVector[0] - tempVector[1];
  dist = DistanceVector.Mag();

  return TVector3(vertex[0], vertex[1], vertex[2]);
}

//_____________________________________________________________________________
TVector3
VertexPointTF2(const TVector3& Xin, const TVector3& Xout,
               const TVector3& Pin, const TVector3& Pout,
               Double_t& dist)
{
  Double_t x_in = Xin.x(), y_in = Xin.y(), z_in = Xin.z();
  Double_t x_out = Xout.x(), y_out = Xout.y(), z_out = Xout.z();
  Double_t u0in = Pin.x()/Pin.z(), v0in = Pin.y()/Pin.z();
  Double_t u0out = Pout.x()/Pout.z(), v0out = Pout.y()/Pout.z();

  Double_t x0in = x_in - u0in*z_in, y0in = y_in - v0in*z_in;
  Double_t x0out = x_out - u0out*z_out, y0out = y_out - v0out*z_out;
  // in function
  // x = [0] + [1]*t
  // y = [2] + [3]*t
  // z = t

  // out function
  // x = [4] + [5]*t
  // y = [6] + [7]*t
  // z = t

  static TF2 fvert("fvert",
                   "pow(([0]+[1]*x)-([4]+[5]*y),2)+pow(([2]+[3]*x)-([6]+[7]*y),2)+pow(x-y,2)",
                   -1000.,1000.,-1000.,1000.);

  fvert.SetParameter(0, x0in);
  fvert.SetParameter(1, u0in);
  fvert.SetParameter(2, y0in);
  fvert.SetParameter(3, v0in);
  fvert.SetParameter(4, x0out);
  fvert.SetParameter(5, u0out);
  fvert.SetParameter(6, y0out);
  fvert.SetParameter(7, v0out);

  Double_t close_zin, close_zout;
  fvert.GetMinimumXY(close_zin, close_zout);
  dist = TMath::Sqrt(fvert.GetMinimum());

  Double_t vertx = ((x0in+u0in*close_zin)+(x0out+u0out*close_zout))/2.;
  Double_t verty = ((y0in+v0in*close_zin)+(y0out+v0out*close_zout))/2.;
  Double_t vertz = (close_zin+close_zout)/2.;

  return TVector3(vertx, verty, vertz);
}

//_____________________________________________________________________________
Double_t
CloseDist(const TVector3& Xin, const TVector3& Xout,
          const TVector3& Pin, const TVector3& Pout)
{
  Double_t xi=Xin.x(), yi=Xin.y(), xo=Xout.x(), yo=Xout.y();
  Double_t ui=Pin.x()/Pin.z(), vi=Pin.y()/Pin.z();
  Double_t uo=Pout.x()/Pout.z(), vo=Pout.y()/Pout.z();

  Double_t z=((xi-xo)*(uo-ui)+(yi-yo)*(vo-vi))/
    ((uo-ui)*(uo-ui)+(vo-vi)*(vo-vi));
  Double_t x1=xi+ui*z, y1=yi+vi*z;
  Double_t x2=xo+uo*z, y2=yo+vo*z;

  return TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

//_____________________________________________________________________________
TVector3
CorrElossIn(const TVector3& Pin, const TVector3& Xin,
            const TVector3& vtx, Double_t mass)
{
  Double_t length = 0.;//, dE;
  Double_t mom_new, energy_new;

  TVector3 CorPin = Pin;
  Double_t mom = Pin.Mag();

  mom_new = mom;
  energy_new = TMath::Sqrt(mass*mass+mom*mom);

  if (std::abs(vtx.z()-TARGETcenter) <= TARGETsizeZ)
    length = CalcLengthBeam(Pin, Xin, vtx);
  else if (vtx.z() > TARGETsizeZ+TARGETcenter){
    Double_t u=Pin.x()/Pin.z();
    Double_t v=Pin.y()/Pin.z();
    Double_t ztgtOut=TARGETcenter+TARGETsizeZ;
    TVector3 tgtOutPos(u*(ztgtOut)+Xin.x(), v*(ztgtOut)+Xin.y(), ztgtOut);
    length = CalcLengthBeam(Pin, Xin, tgtOutPos);
  }
  else if (vtx.z() < -TARGETsizeZ+TARGETcenter)
    return CorPin;

  if (CalcDe(mom, mass, length,&mom_new,&energy_new)) {
    CorPin = mom_new/mom*Pin;
    // hddaq::cout << "CorrElossIn:: mom = " << mom << ", mom_new = " << mom_new
    // 	      << std::endl;
    return CorPin;
  } else {
    return Pin;
  }
}

//_____________________________________________________________________________
Double_t
CalcLengthBeam(const TVector3& Pin, const TVector3& Xin,
               const TVector3& vtx)
{
  Double_t u=Pin.x()/Pin.z();
  Double_t v=Pin.y()/Pin.z();

  // hddaq::cout << "CalcLengthBeam:: vtx (x,y,z)=(" << vtx.x() << ", "
  // 	    << vtx.y() << ", " << vtx.z() << ")" << std::endl;

  TVector3 point1, point2;
  /* vertex point check */
  if (IsInsideTarget(vtx))
    point1=vtx;
  else {
    Double_t z1,z2;
    Double_t z;
    if (CalcCrossingPoint(u, v, vtx,&z1,&z2)) {
      if (std::abs(vtx.z()-z1)<std::abs(vtx.z()-z2))
        z=z1;
      else
        z=z2;

      point1 = TVector3(u*(z-vtx.z())+vtx.x(), v*(z-vtx.z())+vtx.y(), z);
    } else {
      return 0.0;
    }
  }
  // hddaq::cout << "CalcLengthBeam:: Point1 (x,y,z)=(" << point1.x() << ", "
  // 	    << point1.y() << ", " << point1.z() << ")" << std::endl;
  /* target entrance point check */
  Double_t ztgtIn=TARGETcenter-TARGETsizeZ;

  TVector3 tgtInPos(u*(ztgtIn)+Xin.x(), v*(ztgtIn)+Xin.y(), ztgtIn);
  if (IsInsideTarget(tgtInPos))
    point2=tgtInPos;
  else {
    Double_t z1,z2;
    Double_t z;
    if (CalcCrossingPoint(u, v, tgtInPos,&z1,&z2)) {
      if (std::abs(tgtInPos.z()-z1)<std::abs(tgtInPos.z()-z2))
        z=z1;
      else
        z=z2;

      point2 = TVector3(u*(z)+Xin.x(), v*(z)+Xin.y(), z);
    } else {
      return 0.0;
    }
  }

  // hddaq::cout << "CalcLengthBeam:: Point2 (x,y,z)=(" << point2.x() << ", "
  // 	    << point2.y() << ", " << point2.z() << ")" << std::endl;
  // hddaq::cout << "CalcLengthBeam:: length=" << (point1-point2).Mag() << std::endl;

  return (point1-point2).Mag();
}

/*
  Correct Energy loss (Downstream part of Vertex)
*/
//_____________________________________________________________________________
TVector3
CorrElossOut(const TVector3& Pout, const TVector3& Xout,
             const TVector3& vtx, Double_t mass)
{
  Double_t FL,FH,FTMP;
  Double_t Elow, Ehigh, Elast, EPS;
  Double_t E = 0., length = 0.;

  TVector3 CorPout = Pout;
  Double_t mom = Pout.Mag();

  if (std::abs(vtx.z()-TARGETcenter) < TARGETsizeZ)
    length =  CalcLengthScat(Pout, Xout, vtx);
  else if (vtx.z() < -TARGETsizeZ-TARGETcenter){
    Double_t u=Pout.x()/Pout.z();
    Double_t v=Pout.y()/Pout.z();
    Double_t ztgtIn=TARGETcenter-TARGETsizeZ;
    TVector3 tgtInPos(u*(ztgtIn)+Xout.x(), v*(ztgtIn)+Xout.y(), ztgtIn);
    length = CalcLengthScat(Pout, Xout, tgtInPos);
  }
  else if (vtx.z() > TARGETsizeZ+TARGETcenter)
    return CorPout;

  Elow  = mass;
  Ehigh = 10.;
  Elast = TMath::Sqrt(mass*mass+mom*mom);
  EPS = 0.001;
  FL  = DiffE(mass, Elow, length, Elast);
  FH  = DiffE(mass, Ehigh, length, Elast);
  while (std::abs((Ehigh-Elow)/Ehigh) > EPS) {
    E = Ehigh-(Ehigh-Elow)*FH/(FH-FL);
    //printf("-------E=%f, Elow=%f, Ehigh=%f, FL=%f, FH=%f----------\n",E, Elow, Ehigh, FL, FH);
    if (std::abs(FTMP=DiffE(mass, E, length, Elast)) < 0.0000001){
      Elow = E;
      Ehigh = E;
      FH = FTMP;
    } if ((FTMP=DiffE(mass, E, length, Elast)) < 0) {
      Elow = E;
      FL = FTMP;
    } else if ((FTMP=DiffE(mass, E, length, Elast)) > 0){
      Ehigh = E;
      FH = FTMP;
    }
  }

  Double_t energy_new = E;
  Double_t mom_new = TMath::Sqrt(energy_new*energy_new-mass*mass);

  CorPout = mom_new/mom*Pout;
  // hddaq::cout << "CorrElossOut:: mom = " << mom << ", mom_new = " << mom_new
  // 	    << std::endl;

  return CorPout;
}

//_____________________________________________________________________________
Double_t
CalcLengthScat(const TVector3& Pout, const TVector3& Xout,
               const TVector3& vtx)
{
  Double_t u=Pout.x()/Pout.z();
  Double_t v=Pout.y()/Pout.z();

  // hddaq::cout << "CalcLengthScat:: vtx (x,y,z)=(" << vtx.x() << ", "
  // 	    << vtx.y() << ", " << vtx.z() << ")" << std::endl;

  TVector3 point1, point2;
  /* vertex point check */
  if (IsInsideTarget(vtx))
    point1=vtx;
  else {
    Double_t z1,z2;
    Double_t z;
    if (CalcCrossingPoint(u, v, vtx,&z1,&z2)) {
      if (std::abs(vtx.z()-z1)<std::abs(vtx.z()-z2))
        z=z1;
      else
        z=z2;

      point1 = TVector3(u*(z-vtx.z())+vtx.x(), v*(z-vtx.z())+vtx.y(), z);
    } else {
      return 0.0;
    }
  }
  // hddaq::cout << "CalcLengthScat:: Point1 (x,y,z)=(" << point1.x() << ", "
  // 	    << point1.y() << ", " << point1.z() << ")" << std::endl;

  /* target exit point check */
  Double_t ztgtOut=TARGETcenter+TARGETsizeZ;

  TVector3 tgtOutPos(u*(ztgtOut)+Xout.x(), v*(ztgtOut)+Xout.y(), ztgtOut);
  if (IsInsideTarget(tgtOutPos))
    point2=tgtOutPos;
  else {
    Double_t z1,z2;
    Double_t z;
    if (CalcCrossingPoint(u, v, tgtOutPos,&z1,&z2)) {
      if (std::abs(tgtOutPos.z()-z1)<std::abs(tgtOutPos.z()-z2))
        z=z1;
      else
        z=z2;

      point2 = TVector3(u*(z)+Xout.x(), v*(z)+Xout.y(), z);
    } else {
      return 0.0;
    }
  }

  // hddaq::cout << "CalcLengthScat:: Point2 (x,y,z)=(" << point2.x() << ", "
  // 	    << point2.y() << ", " << point2.z() << ")" << std::endl;
  // hddaq::cout << "CalcLengthScata:: length=" << (point1-point2).Mag() << std::endl;

  return (point1-point2).Mag();

}

//_____________________________________________________________________________
TVector3
CorrElossOutCheck(const TVector3& Pout, const TVector3& Xout,
                  const TVector3& vtx, Double_t mass)
{
  Double_t length = 0.;//, dE;
  Double_t energy_new, mom_new;

  Double_t mom=Pout.Mag();

  mom_new = mom;
  energy_new = TMath::Sqrt(mass*mass+mom*mom);

  TVector3 CorPout = Pout;

  if (std::abs(vtx.z()-TARGETcenter) <= TARGETsizeZ)
    length =  CalcLengthBeam(Pout, Xout, vtx);
  else if (vtx.z() < -TARGETsizeZ+TARGETcenter){
    Double_t u=Pout.x()/Pout.z();
    Double_t v=Pout.y()/Pout.z();
    Double_t ztgtIn=TARGETcenter-TARGETsizeZ;
    TVector3 tgtInPos(u*(ztgtIn)+Xout.x(), v*(ztgtIn)+Xout.y(), ztgtIn);
    length = CalcLengthBeam(Pout, Xout, tgtInPos);
  }
  else if (vtx.z() > TARGETsizeZ+TARGETcenter)
    return CorPout;

  if (CalcDe(mom, mass, length,&mom_new,&energy_new)) {
    CorPout = mom_new/mom*Pout;
    return CorPout;
  } else {
    return Pout;
  }
}

//_____________________________________________________________________________
Bool_t
IsInsideTarget(const TVector3&point)
{
#if 1
  return ((std::abs(point.x() - TARGETcenterX) < TARGETsizeX)&&
          (std::abs(point.y() - TARGETcenterY) < TARGETsizeY));
#else
  return (pow((point.x()-TARGETcenterX), 2.) +
          pow((point.y()-TARGETcenterY), 2.) <=
          pow(TARGETradius, 2.));
#endif
  return false;
}

//_____________________________________________________________________________
Bool_t
CalcCrossingPoint(Double_t u, Double_t v, const TVector3& Point,
                  Double_t* z1, Double_t* z2)
{
  Double_t x0=TARGETcenterX, y0=TARGETcenterY;
  Double_t a=Point.x()-u*Point.z()-x0;
  Double_t b=Point.y()-v*Point.z()-y0;
  Double_t r=TARGETradius;
  Double_t c=(a*u+b*v)*(a*u+b*v)-(u*u+v*v)*(a*a+b*b-r*r);

  if(c < 0){
    // hddaq::cerr << FUNC_NAME << " "
    // 	      << "this track does not cross target" << std::endl;
    return false;
  }else{
    // hddaq::cout << FUNC_NAME << " "
    // 	      << "this track cross target!" << std::endl;
    *z1 = (-(a*u+b*v)+TMath::Sqrt(c))/(u*u+v*v);
    *z2 = (-(a*u+b*v)-TMath::Sqrt(c))/(u*u+v*v);
    return true;
  }
}

//_____________________________________________________________________________
Double_t
DiffE(Double_t mass, Double_t E, Double_t length, Double_t Elast)
{
  Double_t p = TMath::Sqrt(E*E-mass*mass);
  Double_t mom_new, energy_new;
  CalcDe(p, mass, length, &mom_new, &energy_new);
  return (energy_new - Elast);
  //return (TMath::Sqrt(mass*mass+p*p)-CalcDe(particle,p,length)-Elast);
}

//_____________________________________________________________________________
Int_t
CalcDe(Double_t momentum, Double_t mass, Double_t distance,
       Double_t *momentum_cor, Double_t *energy_cor)
{
  Double_t dE_dx; /*MeV/cm*/
  Double_t eloss; /*MeV*/

  Double_t beta;
  Double_t thickness = distance/10.0; /*cm*/
  Double_t m = mass*1000.0;  /*mass of incident particle(Mev)*/

  Double_t E_0;
  Double_t E;
  Double_t p = momentum*1000.0;
  //Double_t delta=0.01; /*cm*/
  Double_t delta=0.1; /*cm*/
  Int_t i;
  Double_t length=0.0;
  Double_t total_eloss=0.0;

  //E_0=Gamma(beta)*m;
  //p=Gamma(beta)*beta*m;

  E_0= TMath::Sqrt(m*m + p*p);
  E=E_0;
  beta = Beta(E,p);
  //printf("beta=%f, E=%f, p=%f\n",beta, E, p);
  if (beta<=0.0) {
    *momentum_cor = p/1000.0;
    *energy_cor = E/1000.0;
    return 1;
  }
  dE_dx=0.0;
  eloss=0.0;
  for(i=0;i<=thickness/delta;i++){
    dE_dx=CalcDedx(beta);
    eloss=dE_dx*delta;
    E=E-eloss;
    if(E<m){
      fprintf(stderr,"particle stops in material at %5.3fcm\n",length);
      *momentum_cor = p/1000.0;
      *energy_cor = E/1000.0;
      return 0;
      //break;
    }
    p=TMath::Sqrt(pow(E,2.0)-pow(m,2.0));
    beta=Beta(E,p);
    length=length+delta;
    total_eloss=total_eloss+eloss;
    /*
      printf("beta:%5.3f\n",beta);
      printf("dE_dx:%5.3f\teloss:%5.3f\n",dE_dx,eloss);
      printf("E:%5.3f(MeV)\tp:%5.3f(MeV/c)\n",E,p);
      printf("length:%5.3f(cm)\n",length);
      printf("total energy loss:%5.3f(MeV)\n",total_eloss);
    */
    //getchar();
  }
  *momentum_cor = p/1000.0;
  *energy_cor = E/1000.0;

  return 1;
}

//_____________________________________________________________________________
Double_t
CalcDedx(Double_t beta)
{
  Double_t value;
  const Double_t C=0.1535; /*MeVcm^2/g*/
  const Double_t m_e=0.511;
  Double_t logterm;
  Double_t gamma_2;
  Double_t W_max;
  Double_t gamma;
  Double_t X;
  Double_t delta = 0.;

  Double_t rho=0.0709;   /*g/cm^3 (C)*/
  Double_t I=21.8;     /*eV*/
  Double_t Z_A=0.99216;
  Int_t z=1;
  Double_t C0=-3.2632;
  Double_t X0=0.4759;
  Double_t X1=1.9215;
  Double_t a=0.13483;
  Double_t M=5.6249;

  gamma = Gamma(beta);
  X = log10(beta*gamma);
  if (X<=X0)
    delta=0.0;
  else if (X0<X&& X<X1)
    delta=4.6052*X+C0+a*pow((X1-X),M);
  else if (X>=X1)
    delta=4.6052*X+C0;

  gamma_2=1/(1-pow(beta,2.0));
  //printf("TMath::Sqrt(gamma_2):%f\n",TMath::Sqrt(gamma_2));

  W_max=2.0*m_e*pow(beta,2.0)*gamma_2;

  logterm=log(2.0*m_e*gamma_2*pow(beta,2.0)*W_max*pow(10.0,12.0)/pow(I,2.0))-2.0*pow(beta,2.0)-delta/2.0;

  value=C*rho*Z_A*pow((Double_t)z,2.0)*logterm/pow(beta,2.0);

  return value;

}

//_____________________________________________________________________________
Double_t
Gamma(Double_t beta)
{
  return 1./TMath::Sqrt(1.-beta*beta);
}

//_____________________________________________________________________________
Double_t
Beta(Double_t energy,Double_t mormentum)
{
  return mormentum/energy;
}

}

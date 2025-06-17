// -*- C++ -*-
/*
//Comment by Wooseung
Material lookup table for HypTPC dEdx calculation
0:P10, 1:Polyethylene(Target) 2:Diamond(Target) 3:Polyvinyltoluene (old TPC gas-vessel)
4:Mylar (gas vessel window) 5:Al (gas-vessel frame)
6:Silica aerogel (AC)
*/

#include "Kinematics.hh"

#include <cmath>
#include <cstdio>

#include <TF2.h>
#include <TMath.h>
#include <TMinuit.h>

#include <std_ostream.hh>

#include "FuncName.hh"
#include "MathTools.hh"
#include "TPCPadHelper.hh"

namespace
{
const Double_t qnan = TMath::QuietNaN();
const Double_t z_offset      = 0.0;
const Double_t TARGETcenter  = 0.0;//Z position
const Double_t TARGEThw      = 15.0/2.0;
const Double_t TARGETsizeX   = 25.0/2.0;
const Double_t TARGETsizeY   = 15.0/2.0;
const Double_t TARGETsizeZ   = 15.0/2.0;
const Double_t TARGETcenterX = 0.0;
const Double_t TARGETcenterY = 0.0;
const Double_t TARGETradius  = 67.3/2.0; //Legacy (Not E42 target)

const Double_t conversion_factor = 12171.3; //HypTPC's ADC to <dE/dx>
const Double_t sigma_dedx_pi[5] = {3.94842, 0.0138502, -0.110281, 12.6065, -10.9347};
const Double_t sigma_dedx_k[5] = {6.24543, -3.21037, 1.52683, 127.099, -9.1004};
const Double_t sigma_dedx_p[5] = {12.9717, -8.43799, 3.10608, 166.494, -6.56123};
const Double_t sigma_dedx_d[5] = {29.0934, -15.2229, 2.6452, 500.909, -5.07974};
const Double_t sigma_dedx_t[5] = {28.3878, 3.69034, -4.49067, 0, 0};
const Double_t offset_dedx_d = 5.84462;
const Double_t offset_dedx_t = 14.4277;
const Double_t sigma_dedx_ep = 7.8511;
const Double_t sigma_dedx_em = 8.45029;

const Double_t sigma_tof_pi[5] = {0.450116, 0.0918938, 2.46676e-06, -0.294671, 0.33648};
const Double_t sigma_tof_k[5] = {-0.959259, 1.35794, -0.427667, 1.53197, -1.82723};
const Double_t sigma_tof_p[5] = {0.150481, -0.0433547, 0.0221499, 1.90992, -10.1773};
const Double_t sigma_tof_d[5] = {0.376134, -0.193213, 0.0527238, 0, 0};
const Double_t sigma_tof_t[5] = {0.492751, -0.305449, 0.10159, 0, 0};
const Double_t sigma_tof_ep = 0.129776;
const Double_t sigma_tof_em = 0.14181;
const Double_t offset_tof_ep = -0.237853;
}

static Int_t gNumOfTracks;
static std::vector<Double_t> gX0;
static std::vector<Double_t> gY0;
static std::vector<Double_t> gU0;
static std::vector<Double_t> gV0;
static std::vector<Double_t> gSX0;
static std::vector<Double_t> gSY0;
static std::vector<Double_t> gSU0;
static std::vector<Double_t> gSV0;
static void fcn_vertex(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){

  Double_t chisqr=0.;
  for(Int_t i=0; i<gNumOfTracks; ++i){
    chisqr += TMath::Power(par[0]-gX0[i]-gU0[i]*par[2], 2) / (gSX0[i]*gSX0[i]+par[2]*par[2]*gSU0[i]*gSU0[i]);
    chisqr += TMath::Power(par[1]-gY0[i]-gV0[i]*par[2], 2) / (gSY0[i]*gSY0[i]+par[2]*par[2]*gSV0[i]*gSV0[i]);
  };
  f = chisqr;
};

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
TVector3
VertexPointHelix(const Double_t par1[5], const Double_t par2[5],
                 Double_t& dist, Double_t& t1, Double_t& t2)
{
  //helix function 1
  //x = [0] + [3]*cos(t);
  //y = [1] + [3]*sin(t);
  //z = [2] + [3]*[4]*t;

  //helix function 2
  //x = [5] + [8]*cos(t);
  //y = [6] + [8]*sin(t);
  //z = [7] + [8]*[9]*t;

  static TF2 fvert_helix("fvert_helix",
                         "pow(([0]+[3]*cos(x))-([5]+[8]*cos(y)),2)+pow(([1]+[3]*sin(x))-([6]+[8]*sin(y)),2)+pow(([2]+[3]*[4]*x)-([7]+[8]*[9]*y),2)",
                         -5.,5.,-5.,5.);

  fvert_helix.SetParameter(0, par1[0]);
  fvert_helix.SetParameter(1, par1[1]);
  fvert_helix.SetParameter(2, par1[2]);
  fvert_helix.SetParameter(3, par1[3]);
  fvert_helix.SetParameter(4, par1[4]);
  fvert_helix.SetParameter(5, par2[0]);
  fvert_helix.SetParameter(6, par2[1]);
  fvert_helix.SetParameter(7, par2[2]);
  fvert_helix.SetParameter(8, par2[3]);
  fvert_helix.SetParameter(9, par2[4]);

  Double_t close_zin, close_zout;
  fvert_helix.GetMinimumXY(close_zin, close_zout);
  t1 = close_zin;
  t2 = close_zout;
  dist = TMath::Sqrt(fvert_helix.GetMinimum());

  Double_t xin = par1[0]+par1[3]*cos(close_zin);
  Double_t xout = par2[0]+par2[3]*cos(close_zout);
  Double_t yin =  par1[1]+par1[3]*sin(close_zin);
  Double_t yout = par2[1]+par2[3]*sin(close_zout);
  Double_t zin = par1[2]+par1[3]*par1[4]*close_zin;
  Double_t zout =  par2[2]+par2[3]*par2[4]*close_zout;

  // Double_t vx = (par1[0]+par1[3]*cos(close_zin) + par2[0]+par2[3]*cos(close_zout))/2.;
  // Double_t vy = (par1[1]+par1[3]*sin(close_zin) + par2[1]+par2[3]*sin(close_zout))/2.;
  // Double_t vz = (par1[2]+par1[3]*par1[4]*close_zin + par2[2]+par2[3]*par2[4]*close_zout)/2.;
  Double_t vx = (xin+xout)/2.;
  Double_t vy = (yin+yout)/2.;
  Double_t vz = (zin+zout)/2.;

  Double_t dist2 = sqrt(pow(xin-xout,2)
			+pow(yin-yout,2)
			+pow(zin-zout,2));
  // std::cout<<"dist ="<<dist<<", dist2="<<dist2<<std::endl;
  // std::cout<<"close_zin="<<close_zin<<", close_zout="<<close_zout<<std::endl;
  dist = dist2;
  Double_t vertx = -1.*vx;
  Double_t verty = vz;
  Double_t vertz = vy + tpc::ZTarget;
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
Double_t CalcHelixCloseDist(TVector3 point, Double_t par[5], Double_t t1_start, Double_t t1_end){

  Double_t xi = -1.*point.x();
  Double_t yi = point.z() - tpc::ZTarget;
  Double_t zi = point.y();

  //helix function
  //x = [0] + [3]*cos(t);
  //y = [1] + [3]*sin(t);
  //z = [2] + [3]*[4]*t;

  TF1 fvertex("fvertex", "pow([0]+[3]*cos(x)-[5], 2)+pow([1]+[3]*sin(x)-[6], 2)+pow([2]+[3]*[4]*x-[7], 2)", t1_start, t1_end);
  fvertex.SetParameter(0, par[0]);
  fvertex.SetParameter(1, par[1]);
  fvertex.SetParameter(2, par[2]);
  fvertex.SetParameter(3, par[3]);
  fvertex.SetParameter(4, par[4]);
  fvertex.SetParameter(5, xi);
  fvertex.SetParameter(6, yi);
  fvertex.SetParameter(7, zi);
  Double_t helix_t = fvertex.GetMinimumX(t1_start, t1_end);
  return TMath::Sqrt(fvertex.GetMinimum());
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
    length = CalcLengthBeam(Pout, Xout, vtx);
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
#else //Legacy (Not E42 target)
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

//Legacy
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

/*
  Correct Energy loss in the Target
  inputs : 1. path length through the target  2.momentum vector at the target crossing point
*/

//_____________________________________________________________________________
Double_t DensityEffectCorrection(Double_t betagamma, Double_t *par){

  //Sternheimer's parameterization from the PDG
  //notation : par[0] : a, par[1] : k, par[2] : x0, par[3] : x1, par[4] : _C, par[5] : delta0
    Double_t constant = 2*TMath::Log(10);
    Double_t delta = 0.;
    Double_t X = log10(betagamma);
    if(X<=par[2]) delta = par[5]*TMath::Power(10., 2*(X - par[2]));
    else if(par[2]<X && X<par[3]) delta = constant*X - par[4] + par[0]*pow((par[3] - X), par[1]);
    else if(X>=par[3]) delta = constant*X - par[4];

  return delta;

}

//_____________________________________________________________________________
Double_t HypTPCdEdx(Int_t materialid, Double_t mass/*MeV/c2*/, Double_t beta){

  Double_t rho=0.; //[g cm-3]
  Double_t ZoverA=0.; //[mol g-1]
  Double_t I=0.; //[eV]
  Double_t density_effect_par[6]={0.}; //Sternheimer's parameterization form the PDG
  if(materialid==0){  //P10
    rho = TMath::Power(10.,-3)*(0.9*1.662 + 0.1*0.6672);
    ZoverA = 17.2/37.6;
    I = 0.9*188.0 + 0.1*41.7;
    density_effect_par[0] = 0.9*0.19714 + 0.1*0.09253;
    density_effect_par[1] = 0.9*2.9618 + 0.1*3.6257;
    density_effect_par[2] = 0.9*1.7635 + 0.1*1.6263;
    density_effect_par[3] = 0.9*4.4855 + 0.1*3.9716;
    density_effect_par[4] = 0.9*11.9480 + 0.1*9.5243;
    density_effect_par[5] = 0.;
  }
  else if(materialid==1){ //Polyethylene(Target)
    rho = 1.13;
    ZoverA = 0.57034;
    I = 57.4;
    density_effect_par[0] = 0.12108;
    density_effect_par[1] = 3.4292;
    density_effect_par[2] = 0.1489;
    density_effect_par[3] = 2.5296;
    density_effect_par[4] = 3.0563;
    density_effect_par[5] = 0.;
  }
  else if(materialid==2){ //Diamond(Target)
    rho = 3.223;
    ZoverA = 6./12.0107;
    I = 78.0;
    density_effect_par[0] = 0.26142;
    density_effect_par[1] = 2.8697;
    density_effect_par[2] = -0.1135;
    density_effect_par[3] = 2.2458;
    density_effect_par[4] = 2.4271;
    density_effect_par[5] = 0.12;
  }
  else if(materialid==3){ //Polyvinyltoluene
    rho = 1.032;
    ZoverA = 0.54141;
    I = 64.7;
    density_effect_par[0] = 0.16101;
    density_effect_par[1] = 3.2393;
    density_effect_par[2] = 0.1464;
    density_effect_par[3] = 2.4855;
    density_effect_par[4] = 3.1997;
    density_effect_par[5] = 0.00;
  }
  else if(materialid==4){ //Mylar (gas-vessel window)
    rho = 1.400;
    ZoverA = 0.52037;
    I = 78.7;
    density_effect_par[0] = 0.12679;
    density_effect_par[1] = 3.3076;
    density_effect_par[2] = 0.1562;
    density_effect_par[3] = 2.6507;
    density_effect_par[4] = 3.3262;
    density_effect_par[5] = 0.00;
  }
  else if(materialid==5){ //Al (gas-vessel frame)
    rho = 2.699;
    ZoverA = 0.481811;
    I = 166.0;
    density_effect_par[0] = 0.08024;
    density_effect_par[1] = 3.6345;
    density_effect_par[2] = 0.1708;
    density_effect_par[3] = 3.0127;
    density_effect_par[4] = 4.2395;
    density_effect_par[5] = 0.12;
  }
  else if(materialid==6){ //Slical aerogel
    rho = 0.50093;
    ZoverA = 0.200;
    I = 139.2;
    density_effect_par[0] = 0.26675;
    density_effect_par[1] = 3.0000;
    density_effect_par[2] = 0.6029;
    density_effect_par[3] = 3.0000;
    density_effect_par[4] = 6.4507;
    density_effect_par[5] = 0.00;
  }

  Double_t Z = 1.;
  Double_t me = 0.5109989461; //[MeV]
  Double_t K = 0.307075; //[MeV cm2 mol-1]
  Double_t constant = rho*K*ZoverA; //[MeV cm-1]
  Double_t I2 = I*I; //Mean excitaion energy [eV]
  Double_t beta2 = beta*beta;
  Double_t gamma2 = 1./(1.-beta2);
  Double_t MeVToeV = TMath::Power(10.,6);
  Double_t Wmax = 2*me*beta2*gamma2/((me/mass+1.)*(me/mass+1.)+2*(me/mass)*(TMath::Sqrt(gamma2)-1));
  Double_t delta = DensityEffectCorrection(TMath::Sqrt(beta2*gamma2), density_effect_par);
  Double_t dedx = constant*Z*Z/beta2*(0.5*TMath::Log(2*me*beta2*gamma2*Wmax*MeVToeV*MeVToeV/I2) - beta2 - 0.5*delta);
  return dedx;

}

//_____________________________________________________________________________
Int_t
HypTPCCalcDe(Int_t materialid, Double_t momentum, Double_t mass, Double_t distance,
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
  Double_t delta=0.01; /*cm*/
  //Double_t delta=0.1; /*cm*/
  Int_t i;
  Double_t length=0.0;
  Double_t total_eloss=0.0;

  //E_0=Gamma(beta)*m;
  //p=Gamma(beta)*beta*m;

  E_0= TMath::Sqrt(m*m + p*p);
  E=E_0;
  beta = Beta(E,p);
  //printf("beta=%f, E=%f, p=%f\n",beta, E, p);
  if(beta<=0.0) {
    *momentum_cor = p/1000.0;
    *energy_cor = E/1000.0;
    return 1;
  }
  dE_dx=0.0;
  eloss=0.0;
  for(i=0;i<=thickness/delta;i++){
    dE_dx=HypTPCdEdx(materialid, m, beta);
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
HypTPCDiffE(Int_t materialid, Double_t mass, Double_t E, Double_t length, Double_t Elast)
{

  Double_t p = TMath::Sqrt(E*E-mass*mass);
  Double_t mom_new, energy_new;
  HypTPCCalcDe(materialid, p, mass, length, &mom_new, &energy_new);
  return (energy_new - Elast);

}

//_____________________________________________________________________________
TVector3
HypTPCCorrElossOut(Int_t materialid, const TVector3& Pout, Double_t length, Double_t mass)
{

  TVector3 CorPout = Pout;
  if(length>0. && !TMath::IsNaN(mass)){
    Double_t FL,FH,FTMP;
    Double_t Elow, Ehigh, Elast, EPS;
    Double_t E = 0.;
    Double_t mom = Pout.Mag();
    Elow  = mass;
    Ehigh = 10.;
    Elast = TMath::Sqrt(mass*mass+mom*mom);
    EPS = 0.001;
    FL  = HypTPCDiffE(materialid, mass, Elow, length, Elast);
    FH  = HypTPCDiffE(materialid, mass, Ehigh, length, Elast);
    while (std::abs((Ehigh-Elow)/Ehigh) > EPS) {
      E = Ehigh-(Ehigh-Elow)*FH/(FH-FL);
      //printf("-------E=%f, Elow=%f, Ehigh=%f, FL=%f, FH=%f----------\n",E, Elow, Ehigh, FL, FH);
      if (std::abs(FTMP=HypTPCDiffE(materialid, mass, E, length, Elast)) < 0.0000001){
	Elow = E;
	Ehigh = E;
	FH = FTMP;
      } if ((FTMP=HypTPCDiffE(materialid, mass, E, length, Elast)) < 0) {
	Elow = E;
	FL = FTMP;
      } else if ((FTMP=HypTPCDiffE(materialid, mass, E, length, Elast)) > 0){
	Ehigh = E;
	FH = FTMP;
      }
    }

    Double_t energy_new = E;
    Double_t mom_new = TMath::Sqrt(energy_new*energy_new-mass*mass);
    CorPout = mom_new/mom*Pout;
    //hddaq::cout<<"CorrElossOut:: length = "<<length<<" mm"<<std::endl;
    //hddaq::cout<<"CorrElossOut:: mom = "<<mom<<", mom_new = "<<mom_new<<" dP = "<<mom_new - mom<<std::endl;
    //hddaq::cout << "CorrElossOut:: E = " << TMath::Sqrt(mass*mass+mom*mom) << ", E_new = " << energy_new << " dE = "<< energy_new - TMath::Sqrt(mass*mass+mom*mom)<<std::endl;
  }

  return CorPout;

}

//_____________________________________________________________________________
TVector3
HypTPCCorrElossIn(Int_t materialid, const TVector3& Pin, Double_t length, Double_t mass)
{

  Double_t mom = Pin.Mag();
  Double_t mom_new, energy_new;
  TVector3 CorPin = Pin;
  if(length>0. && HypTPCCalcDe(materialid, mom, mass, length, &mom_new, &energy_new)){
    CorPin = mom_new/mom*Pin;
    //hddaq::cout<<"CorrElossIn:: length = "<<length<<" [mm]"<<std::endl;
    //hddaq::cout<<"CorrElossIn:: mom = "<<mom<<", mom_new = "<<mom_new<<" dP = "<<mom_new - mom<<std::endl;
    //hddaq::cout << "CorrElossIn:: E = " << TMath::Sqrt(mass*mass+mom*mom) << ", mom_new = " << energy_new << " dE = "<< energy_new - TMath::Sqrt(mass*mass+mom*mom)<<std::endl;
  }

  return CorPin;
}

//_____________________________________________________________________________
Double_t HypTPCBethe(Double_t *x, Double_t *p){

  //x[0] : poq [GeV/c]
  //p[0] : converting factor
  //p[1] : mass [MeV/c2]
  Double_t momentum = 1000.*TMath::Abs(x[0]); /*MeV/c2*/
  Double_t beta = Beta(TMath::Sqrt(p[1]*p[1] + momentum*momentum), momentum);
  Double_t dedx = p[0]*HypTPCdEdx(0, p[1], beta); //P10

  return dedx;
}
/* Legacy
//_____________________________________________________________________________
Int_t HypTPCdEdxPID(Double_t dedx, Double_t poq){

  Double_t mpi = 139.57039;
  //Double_t mk  = 493.677;
  Double_t mp  = 938.2720813;
  //Double_t md  = 1875.612762;

  // 1 sigma of <dE/dx>_pi
  Double_t par_pi[2] = {conversion_factor, mpi};
  Double_t dedx_pi = HypTPCBethe(&poq, par_pi); //P10's <dE/dx>_pi
  Double_t sigma_pi = (sigma_dedx_pi[0] + sigma_dedx_pi[1]*TMath::Abs(poq) + sigma_dedx_pi[2]*poq*poq);

  // 1 sigma of <dE/dx>_p
  Double_t par_p[2] = {conversion_factor, mp};
  Double_t dedx_p = HypTPCBethe(&poq, par_p); //P10's <dE/dx>_p
  Double_t sigma_p = (sigma_dedx_p[0] + sigma_dedx_p[1]*TMath::Abs(poq) + sigma_dedx_p[2]*poq*poq);

  //p/pi separation power calculation
  Double_t avg_sigma = 0.5*(sigma_pi + sigma_p);
  Double_t dedx_diff = TMath::Abs(dedx_pi - dedx_p);
  Double_t separation_power = dedx_diff/avg_sigma;
  Double_t ppi_separation_cut = dedx_pi < dedx_p ? dedx_pi + separation_power*sigma_pi : dedx_p + separation_power*sigma_p;

  Int_t pid[3] = {0};
  if(separation_power < 3.){ //3 sigma separation limit
    pid[0]=1; pid[1]=1; pid[2]=1;
  }
  else{
    //p pi separation
    if(dedx > ppi_separation_cut) pid[2]=1;
    else pid[0]=1;
    //for Kaon candidate selection
    if(dedx_p > dedx) pid[1]=1;
  }

  Int_t output = pid[0] + pid[1]*2 + pid[2]*4;
  return output;
}
*/

//Legacy
//_____________________________________________________________________________
Bool_t HypTPCdEdxPID_IsKaonTemp(Double_t dedx, Double_t poq){

  Double_t mk  = 493.677;
  Double_t par_k[2] = {conversion_factor, mk};
  Double_t dedx_k = HypTPCBethe(&poq, par_k); //P10's <dE/dx>_k
  // 1 sigma of <dE/dx>_pi
  Double_t sigma_pi = (sigma_dedx_pi[0] + sigma_dedx_pi[1]*TMath::Abs(poq) + sigma_dedx_pi[2]*poq*poq); //temporary using pi's sigma value

  Bool_t flag = false;
  if(TMath::Abs((dedx-dedx_k)/sigma_pi) < 3) flag = true;
  return flag;
}

//_____________________________________________________________________________
Double_t HypTPCdEdxNsigmaProton(Double_t dedx, Double_t poq){

  Double_t mp = 938.2720813;
  Double_t par_p[2] = {conversion_factor, mp};
  Double_t dedx_p = HypTPCBethe(&poq, par_p); //P10's <dE/dx>_p
  // 1 sigma of <dE/dx>_p
  Double_t sigma_p = (sigma_dedx_p[0] + sigma_dedx_p[1]*TMath::Abs(poq) +
		      sigma_dedx_p[2]*poq*poq + sigma_dedx_p[3]*TMath::Exp(sigma_dedx_p[4]*TMath::Abs(poq)));

  Double_t nsigma = (dedx-dedx_p)/sigma_p;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCdEdxNsigmaDeutron(Double_t dedx, Double_t poq){

  Double_t md = 1875.612762;
  Double_t par_d[2] = {conversion_factor, md};
  Double_t dedx_d = HypTPCBethe(&poq, par_d) + offset_dedx_d; //P10's <dE/dx>_d
  // 1 sigma of <dE/dx>_d
  Double_t sigma_d = (sigma_dedx_d[0] + sigma_dedx_d[1]*TMath::Abs(poq) +
		      sigma_dedx_d[2]*poq*poq + sigma_dedx_d[3]*TMath::Exp(sigma_dedx_d[4]*TMath::Abs(poq)));

  Double_t nsigma = (dedx-dedx_d)/sigma_d;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCdEdxNsigmaTriton(Double_t dedx, Double_t poq){

  Double_t mt = 2808.921112;
  Double_t par_t[2] = {conversion_factor, mt};
  Double_t dedx_t = HypTPCBethe(&poq, par_t) + offset_dedx_t; //P10's <dE/dx>_t
  // 1 sigma of <dE/dx>_t
  Double_t sigma_t = (sigma_dedx_t[0] + sigma_dedx_t[1]*TMath::Abs(poq) +
		      sigma_dedx_t[2]*poq*poq + sigma_dedx_t[3]*TMath::Exp(sigma_dedx_t[4]*TMath::Abs(poq)));

  Double_t nsigma = (dedx-dedx_t)/sigma_t;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCdEdxNsigmaKaon(Double_t dedx, Double_t poq){

  Double_t mk = 493.677;
  Double_t par_k[2] = {conversion_factor, mk};
  Double_t dedx_k = HypTPCBethe(&poq, par_k); //P10's <dE/dx>_k
  // 1 sigma of <dE/dx>_k
  Double_t sigma_k = (sigma_dedx_k[0] + sigma_dedx_k[1]*TMath::Abs(poq) +
		      sigma_dedx_k[2]*poq*poq + sigma_dedx_k[3]*TMath::Exp(sigma_dedx_k[4]*TMath::Abs(poq)));

  Double_t nsigma = (dedx-dedx_k)/sigma_k;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCdEdxNsigmaPion(Double_t dedx, Double_t poq){

  Double_t mpi = 139.57039;
  Double_t par_pi[2] = {conversion_factor, mpi};
  Double_t dedx_pi = HypTPCBethe(&poq, par_pi); //P10's <dE/dx>_pi
  // 1 sigma of <dE/dx>_pi
  Double_t sigma_pi = (sigma_dedx_pi[0] + sigma_dedx_pi[1]*TMath::Abs(poq) +
		       sigma_dedx_pi[2]*poq*poq + sigma_dedx_pi[3]*TMath::Exp(sigma_dedx_pi[4]*TMath::Abs(poq)));

  Double_t nsigma = (dedx-dedx_pi)/sigma_pi;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCdEdxNsigmaElectron(Double_t dedx, Double_t poq){

  Double_t me = 0.5109989461; //[MeV]
  Double_t par_e[2] = {conversion_factor, me};
  Double_t dedx_e = HypTPCBethe(&poq, par_e); //P10's <dE/dx>_e
  // 1 sigma of <dE/dx>_e
  Double_t sigma_dedx_e = poq>0 ? sigma_dedx_ep : sigma_dedx_em;

  Double_t nsigma = (dedx-dedx_e)/sigma_dedx_e;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCHTOFNsigmaProton(Double_t poq, Double_t tracklength, Double_t tof){

  Double_t mom = TMath::Abs(poq);
  Double_t mp = 0.001*938.2720813;
  Double_t inverse_beta = tof*MathTools::C()/tracklength;
  Double_t inverse_beta_p = TMath::Sqrt(mp*mp + mom*mom)/mom;
  Double_t sigma_p = (sigma_tof_p[0] + sigma_tof_p[1]*TMath::Abs(mom) +
		      sigma_tof_p[2]*mom*mom + sigma_tof_p[3]*TMath::Exp(sigma_tof_p[4]*TMath::Abs(mom)));

  double nsigma = (inverse_beta-inverse_beta_p)/sigma_p;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCHTOFNsigmaDeutron(Double_t poq, Double_t tracklength, Double_t tof){

  Double_t mom = TMath::Abs(poq);
  Double_t md = 0.001*1875.612762;
  Double_t inverse_beta = tof*MathTools::C()/tracklength;
  Double_t inverse_beta_d = TMath::Sqrt(md*md + mom*mom)/mom;
  Double_t sigma_d = (sigma_tof_d[0] + sigma_tof_d[1]*TMath::Abs(mom) +
		      sigma_tof_d[2]*mom*mom + sigma_tof_d[3]*TMath::Exp(sigma_tof_d[4]*TMath::Abs(mom)));

  double nsigma = (inverse_beta-inverse_beta_d)/sigma_d;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCHTOFNsigmaTriton(Double_t poq, Double_t tracklength, Double_t tof){

  Double_t mom = TMath::Abs(poq);
  Double_t mt = 0.001*2808.921112;
  Double_t inverse_beta = tof*MathTools::C()/tracklength;
  Double_t inverse_beta_t = TMath::Sqrt(mt*mt + mom*mom)/mom;
  Double_t sigma_t = (sigma_tof_t[0] + sigma_tof_t[1]*TMath::Abs(mom) +
		      sigma_tof_t[2]*mom*mom + sigma_tof_t[3]*TMath::Exp(sigma_tof_t[4]*TMath::Abs(mom)));

  double nsigma = (inverse_beta-inverse_beta_t)/sigma_t;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCHTOFNsigmaKaon(Double_t poq, Double_t tracklength, Double_t tof){

  Double_t mom = TMath::Abs(poq);
  Double_t mk = 0.001*493.677;
  Double_t inverse_beta = tof*MathTools::C()/tracklength;
  Double_t inverse_beta_k = TMath::Sqrt(mk*mk + mom*mom)/mom;
  Double_t sigma_k = (sigma_tof_k[0] + sigma_tof_k[1]*TMath::Abs(mom) +
		      sigma_tof_k[2]*mom*mom + sigma_tof_k[3]*TMath::Exp(sigma_tof_k[4]*TMath::Abs(mom)));

  double nsigma = (inverse_beta-inverse_beta_k)/sigma_k;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCHTOFNsigmaPion(Double_t poq, Double_t tracklength, Double_t tof){

  Double_t mom = TMath::Abs(poq);
  Double_t mpi = 0.001*139.57039;
  Double_t inverse_beta = tof*MathTools::C()/tracklength;
  Double_t inverse_beta_pi = TMath::Sqrt(mpi*mpi + mom*mom)/mom;
  Double_t sigma_pi = (sigma_tof_pi[0] + sigma_tof_pi[1]*TMath::Abs(mom) +
		       sigma_tof_pi[2]*mom*mom + sigma_tof_pi[3]*TMath::Exp(sigma_tof_pi[4]*TMath::Abs(mom)));

  double nsigma = (inverse_beta-inverse_beta_pi)/sigma_pi;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCHTOFNsigmaElectron(Double_t poq, Double_t tracklength, Double_t tof){

  Double_t mom = TMath::Abs(poq);
  Double_t me = 0.001*0.5109989461; //[MeV]
  Double_t inverse_beta = tof*MathTools::C()/tracklength;
  Double_t inverse_beta_e = TMath::Sqrt(me*me + mom*mom)/mom;
  Double_t sigma_e = poq>0 ? sigma_tof_ep : sigma_tof_em;

  double nsigma = (inverse_beta-inverse_beta_e)/sigma_e;
  if(poq>0) nsigma = (inverse_beta-inverse_beta_e-offset_tof_ep)/sigma_e;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCHTOFNsigmaProton(Double_t poq, Double_t inverse_beta){

  Double_t mom = TMath::Abs(poq);
  Double_t mp = 0.001*938.2720813;
  Double_t inverse_beta_p = TMath::Sqrt(mp*mp + mom*mom)/mom;
  Double_t sigma_p = (sigma_tof_p[0] + sigma_tof_p[1]*TMath::Abs(mom) +
		      sigma_tof_p[2]*mom*mom + sigma_tof_p[3]*TMath::Exp(sigma_tof_p[4]*TMath::Abs(mom)));

  double nsigma = (inverse_beta-inverse_beta_p)/sigma_p;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCHTOFNsigmaDeutron(Double_t poq, Double_t inverse_beta){

  Double_t mom = TMath::Abs(poq);
  Double_t md = 0.001*1875.612762;
  Double_t inverse_beta_d = TMath::Sqrt(md*md + mom*mom)/mom;
  Double_t sigma_d = (sigma_tof_d[0] + sigma_tof_d[1]*TMath::Abs(mom) +
		      sigma_tof_d[2]*mom*mom + sigma_tof_d[3]*TMath::Exp(sigma_tof_d[4]*TMath::Abs(mom)));

  double nsigma = (inverse_beta-inverse_beta_d)/sigma_d;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCHTOFNsigmaTriton(Double_t poq, Double_t inverse_beta){

  Double_t mom = TMath::Abs(poq);
  Double_t mt = 0.001*2808.921112;
  Double_t inverse_beta_t = TMath::Sqrt(mt*mt + mom*mom)/mom;
  Double_t sigma_t = (sigma_tof_t[0] + sigma_tof_t[1]*TMath::Abs(mom) +
		      sigma_tof_t[2]*mom*mom + sigma_tof_t[3]*TMath::Exp(sigma_tof_t[4]*TMath::Abs(mom)));

  double nsigma = (inverse_beta-inverse_beta_t)/sigma_t;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCHTOFNsigmaKaon(Double_t poq, Double_t inverse_beta){

  Double_t mom = TMath::Abs(poq);
  Double_t mk = 0.001*493.677;
  Double_t inverse_beta_k = TMath::Sqrt(mk*mk + mom*mom)/mom;
  Double_t sigma_k = (sigma_tof_k[0] + sigma_tof_k[1]*TMath::Abs(mom) +
		      sigma_tof_k[2]*mom*mom + sigma_tof_k[3]*TMath::Exp(sigma_tof_k[4]*TMath::Abs(mom)));

  double nsigma = (inverse_beta-inverse_beta_k)/sigma_k;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCHTOFNsigmaPion(Double_t poq, Double_t inverse_beta){

  Double_t mom = TMath::Abs(poq);
  Double_t mpi = 0.001*139.57039;
  Double_t inverse_beta_pi = TMath::Sqrt(mpi*mpi + mom*mom)/mom;
  Double_t sigma_pi = (sigma_tof_pi[0] + sigma_tof_pi[1]*TMath::Abs(mom) +
		       sigma_tof_pi[2]*mom*mom + sigma_tof_pi[3]*TMath::Exp(sigma_tof_pi[4]*TMath::Abs(mom)));

  double nsigma = (inverse_beta-inverse_beta_pi)/sigma_pi;
  return nsigma;
}

//_____________________________________________________________________________
Double_t HypTPCHTOFNsigmaElectron(Double_t poq, Double_t inverse_beta){

  Double_t mom = TMath::Abs(poq);
  Double_t me = 0.001*0.5109989461; //[MeV]
  Double_t inverse_beta_e = TMath::Sqrt(me*me + mom*mom)/mom;
  Double_t sigma_e = poq>0 ? sigma_tof_ep : sigma_tof_em;

  double nsigma = (inverse_beta-inverse_beta_e)/sigma_e;
  if(poq>0) nsigma = (inverse_beta-inverse_beta_e-offset_tof_ep)/sigma_e;
  return nsigma;
}

//_____________________________________________________________________________
Bool_t HypTPCdEdxElectron(Double_t dedx, Double_t poq){

  Double_t nsigma_pi = HypTPCdEdxNsigmaPion(dedx, poq);
  Double_t nsigma_e = HypTPCdEdxNsigmaElectron(dedx, poq);

  Double_t me = 0.5109989461; //[MeV]
  Double_t par_e[2] = {conversion_factor, me};
  Double_t dedx_e = HypTPCBethe(&poq, par_e); //P10's <dE/dx>_e

  Bool_t flag = (nsigma_pi < -4 &&
		 TMath::Abs(nsigma_e) < 3.5 &&
		 TMath::Abs(poq) < 0.1);
  return flag;
}

//_____________________________________________________________________________
Bool_t HypTPCdEdxKaon(Double_t dedx, Double_t poq){

  Double_t nsigma = HypTPCdEdxNsigmaKaon(dedx, poq);
  Double_t window_k[2] = {-3., 3.};
  Bool_t flag = (nsigma < window_k[1] && nsigma > window_k[0]);
  return flag;
}

//_____________________________________________________________________________
Int_t HypTPCdEdxPID(Double_t dedx, Double_t poq){

  Double_t mpi = 139.57039;
  //Double_t mk  = 493.677;
  Double_t mp  = 938.2720813;
  //Double_t md  = 1875.612762;
  if(HypTPCdEdxElectron(dedx, poq)) return 0; //electron

  // 1 sigma of <dE/dx>_pi
  Double_t par_pi[2] = {conversion_factor, mpi};
  Double_t dedx_pi = HypTPCBethe(&poq, par_pi); //P10's <dE/dx>_pi
  Double_t sigma_pi = (sigma_dedx_pi[0] + sigma_dedx_pi[1]*TMath::Abs(poq) + sigma_dedx_pi[2]*poq*poq + sigma_dedx_pi[3]*TMath::Exp(sigma_dedx_pi[4]*TMath::Abs(poq)));
  Double_t nsigma_pi = HypTPCdEdxNsigmaPion(dedx, poq);
  Double_t window_pi[2] = {-3., 3.};

  // 1 sigma of <dE/dx>_p
  Double_t par_p[2] = {conversion_factor, mp};
  Double_t dedx_p = HypTPCBethe(&poq, par_p); //P10's <dE/dx>_p
  Double_t sigma_p = (sigma_dedx_p[0] + sigma_dedx_p[1]*TMath::Abs(poq) + sigma_dedx_p[2]*poq*poq + sigma_dedx_p[3]*TMath::Exp(sigma_dedx_p[4]*TMath::Abs(poq)));
  Double_t nsigma_p = HypTPCdEdxNsigmaProton(dedx, poq);
  Double_t window_p[2] = {-3., 6.};

  //p/pi separation power calculation
  Double_t avg_sigma = 0.5*(sigma_pi + sigma_p);
  Double_t dedx_diff = TMath::Abs(dedx_pi - dedx_p);
  Double_t separation_power = dedx_diff/avg_sigma;
  Double_t ppi_separation_cut = dedx_pi < dedx_p ? dedx_pi + 0.5*separation_power*sigma_pi : dedx_p + 0.5*separation_power*sigma_p;

  Int_t pid[3] = {0};
  if(separation_power < 6.){ //3 sigma separation limit
    if(nsigma_pi > window_pi[0] && nsigma_pi < window_pi[1]) pid[0]=1;
    if(nsigma_p > window_p[0]) pid[2]=1;
  }
  else{
    //p pi separation
    if(dedx > ppi_separation_cut) pid[2]=1;
    else pid[0]=1;
  }

  //for Kaon candidate selection
  Double_t nsigma_k = HypTPCdEdxNsigmaKaon(dedx, poq);
  Double_t window_k[2] = {-3., 3.};
  if(nsigma_k>window_k[0] && nsigma_k<window_k[1]) pid[1] = 1;

  Int_t output = pid[0] + pid[1]*2 + pid[2]*4;
  return output;
}

//_____________________________________________________________________________
void HypTPCPID_PDGCode(Int_t charge, Int_t pid, std::vector<Int_t>& pdg){

  const Int_t particles = 3;
  Int_t pdgcode[particles] = {211, 321, 2212};
  Int_t flag = 1;
  for(Int_t i=0;i<particles;i++){
    Int_t temp = flag&pid;
    if(temp==flag) pdg.push_back(charge*pdgcode[i]);
    if(pid==0) pdg.push_back(charge*pdgcode[i]);
    flag*=2;
  }
}

//_____________________________________________________________________________
TVector3
CalcHelixMom(Double_t Bfield, Int_t charge, Double_t par[5], Double_t t){

  Double_t pt = fabs(par[3])*tpc::ConstC*Bfield;
  Double_t tmp_px = pt*(-1.*sin(t));
  Double_t tmp_py = pt*(cos(t));
  Double_t tmp_pz = pt*(par[4]);
  Double_t px = -tmp_px*0.001;
  Double_t py = tmp_pz*0.001;
  Double_t pz = tmp_py*0.001;

  TVector3 mom(px, py, pz);
  if(charge < 0) mom *= -1.;
  return mom;
}

//_____________________________________________________________________________
void
CalcHelixParam(Double_t Bfield, Int_t charge, TVector3 mom, TVector3 pos, Double_t *par){

  double q = charge;
  Double_t tmp_px = -q*1000.*mom.x();
  Double_t tmp_py = q*1000.*mom.z();
  Double_t tmp_pz = q*1000.*mom.y();
  TVector3 pT(tmp_px, tmp_py, 0.);
  par[3] = pT.Mag()/(tpc::ConstC*Bfield);
  par[4] = tmp_pz/pT.Mag();

  TVector3 pos_(-pos.x(), pos.z() - tpc::ZTarget, pos.y());
  TVector3 norm(0., 0., 1.);
  TVector3 dummy = norm.Cross(pT);
  dummy.SetMag(par[3]);
  TVector3 center = pos_ + dummy;
  par[0] = center.x();
  par[1] = center.y();
  double theta = TMath::ATan2(pos_.y() - par[1], pos_.x() - par[0]);
  par[2] = pos_.z() - (par[4]*par[3]*theta);

}

//______________________________________________________________________________
TVector3
CalcHelixPosition(double par[5], double t)
{
  //This is the eqation of Helix
  double x = par[0] + par[3]*cos(t);
  double y = par[1] + par[3]*sin(t);
  double z = par[2] + (par[4]*par[3]*t);
  TVector3 calpos(-x, z, y + tpc::ZTarget); //local to global coordinate

  return calpos;
}

//_____________________________________________________________________________
TVector3
VertexPointHelix(Double_t par1[5], Double_t par2[5], Double_t t1_start, Double_t t1_end, Double_t t2_start, Double_t t2_end, Double_t& t1, Double_t& t2, Double_t& dist){

  //helix function 1
  //x = [0] + [3]*cos(t);
  //y = [1] + [3]*sin(t);
  //z = [2] + [3]*[4]*t;

  //helix function 2
  //x = [5] + [8]*cos(t);
  //y = [6] + [8]*sin(t);
  //z = [7] + [8]*[9]*t;

  TF2 fvertex_helix("fvertex_helix", "pow(([0]+[3]*cos(x))-([5]+[8]*cos(y)),2)+pow(([1]+[3]*sin(x))-([6]+[8]*sin(y)),2)+pow(([2]+[3]*[4]*x)-([7]+[8]*[9]*y),2)", t1_start, t1_end, t2_start, t2_end);

  fvertex_helix.SetParameter(0, par1[0]);
  fvertex_helix.SetParameter(1, par1[1]);
  fvertex_helix.SetParameter(2, par1[2]);
  fvertex_helix.SetParameter(3, par1[3]);
  fvertex_helix.SetParameter(4, par1[4]);
  fvertex_helix.SetParameter(5, par2[0]);
  fvertex_helix.SetParameter(6, par2[1]);
  fvertex_helix.SetParameter(7, par2[2]);
  fvertex_helix.SetParameter(8, par2[3]);
  fvertex_helix.SetParameter(9, par2[4]);

  Double_t close_zin, close_zout;
  fvertex_helix.GetMinimumXY(close_zin, close_zout);
  t1 = close_zin;
  t2 = close_zout;

  Double_t xin = par1[0]+par1[3]*cos(close_zin);
  Double_t xout = par2[0]+par2[3]*cos(close_zout);
  Double_t yin =  par1[1]+par1[3]*sin(close_zin);
  Double_t yout = par2[1]+par2[3]*sin(close_zout);
  Double_t zin = par1[2]+par1[3]*par1[4]*close_zin;
  Double_t zout = par2[2]+par2[3]*par2[4]*close_zout;

  Double_t vx = (xin+xout)/2.;
  Double_t vy = (yin+yout)/2.;
  Double_t vz = (zin+zout)/2.;

  dist = sqrt(pow(xin-xout,2)
	      +pow(yin-yout,2)
	      +pow(zin-zout,2));
  Double_t vertx = -1.*vx;
  Double_t verty = vz;
  Double_t vertz = vy + tpc::ZTarget;
  return TVector3(vertx, verty, vertz);
}

//_____________________________________________________________________________
TVector3
LambdaVertex(Double_t Bfield, Double_t p_par[5], Double_t pi_par[5],
	     Double_t p_theta_min, Double_t p_theta_max,
	     Double_t pi_theta_min, Double_t pi_theta_max,
	     TVector3 &p_mom, TVector3 &pi_mom, TVector3 &lambda_mom,
	     Double_t& dist){

  Double_t t1, t2;
  TVector3 vertex = VertexPointHelix(p_par, pi_par, p_theta_min, p_theta_max, pi_theta_min, pi_theta_max, t1, t2, dist);
  p_mom = CalcHelixMom(Bfield, 1, p_par, t1);
  pi_mom = CalcHelixMom(Bfield, -1, pi_par, t2);
  lambda_mom = p_mom + pi_mom;

  return vertex;
}

//_____________________________________________________________________________
TVector3 XiVertex(Double_t Bfield, Double_t pi_par[5],
		  Double_t theta_min, Double_t theta_max,
		  TVector3 Xlambda, TVector3 Plambda,
		  TVector3 &Ppi, Double_t &lambdapi_dist){

  Double_t lambdavtx_xivtx_cut = 0.;

  Double_t xi = -1.*Xlambda.x();
  Double_t yi = Xlambda.z() - tpc::ZTarget;
  Double_t zi = Xlambda.y();
  Double_t pxi = -1.*Plambda.x();
  Double_t pyi = Plambda.z();
  Double_t pzi = Plambda.y();
  Double_t ui = -pxi/pyi, vi = pzi/pyi;

  TVector3 p_L = TVector3(pxi, pyi, pzi);
  TVector3 p_unit = p_L.Unit();

  //helix function
  //x = [0] + [3]*cos(t);
  //y = [1] + [3]*sin(t);
  //z = [2] + [3]*[4]*t;

  //straight function
  //x = [5] + [6]*y;
  //z = [7] + [8]*y;

  //TF2 fvertex_helix_linear("fvertex_helix_linear", "pow(([0]+[3]*cos(x))-([5]+[6]*y), 2)+pow(([1]+[3]*sin(x))-y, 2)+pow(([2]+[3]*[4]*x)-([7]+[8]*y), 2)", theta_min, theta_max, -250.-tpc::ZTarget, 250.-tpc::ZTarget);

  Double_t scan_range[2] ={-250. - tpc::ZTarget, 250. - tpc::ZTarget};
  if(pyi>0) scan_range[1] = yi + lambdavtx_xivtx_cut/(p_unit.y());
  else scan_range[0] = yi - lambdavtx_xivtx_cut/(p_unit.y());
  TF2 fvertex_helix_linear("fvertex_helix_linear", "pow(([0]+[3]*cos(x))-([5]+[6]*y), 2)+pow(([1]+[3]*sin(x))-y, 2)+pow(([2]+[3]*[4]*x)-([7]+[8]*y), 2)", theta_min, theta_max, scan_range[0], scan_range[1]);
  //TF2 fvertex_helix_linear("fvertex_helix_linear", "pow(([0]+[3]*cos(x))-([5]+[6]*y), 2)+pow(([1]+[3]*sin(x))-y, 2)+pow(([2]+[3]*[4]*x)-([7]+[8]*y), 2)", theta_min, theta_max, -250.-tpc::ZTarget, 250.-tpc::ZTarget);

  fvertex_helix_linear.SetParameter(0, pi_par[0]);
  fvertex_helix_linear.SetParameter(1, pi_par[1]);
  fvertex_helix_linear.SetParameter(2, pi_par[2]);
  fvertex_helix_linear.SetParameter(3, pi_par[3]);
  fvertex_helix_linear.SetParameter(4, pi_par[4]);
  fvertex_helix_linear.SetParameter(5, xi + ui*yi);
  fvertex_helix_linear.SetParameter(6, -ui);
  fvertex_helix_linear.SetParameter(7, zi - vi*yi);
  fvertex_helix_linear.SetParameter(8, vi);

  Double_t helix_t, close_y;
  fvertex_helix_linear.GetMinimumXY(helix_t, close_y);
  //lambdapi_dist = TMath::Sqrt(fvertex_helix_linear.GetMinimum());
  lambdapi_dist = TMath::Sqrt(fvertex_helix_linear.Eval(helix_t, close_y));

  Ppi = CalcHelixMom(Bfield, -1, pi_par, helix_t);

  Double_t xPi = pi_par[0]+pi_par[3]*cos(helix_t);
  Double_t yPi = pi_par[1]+pi_par[3]*sin(helix_t);
  Double_t zPi = pi_par[2]+pi_par[3]*pi_par[4]*helix_t;
  Double_t xL = xi - ui*(close_y-yi);
  Double_t yL = close_y;
  Double_t zL = zi + vi*(close_y-yi);
  Double_t vx = 0.5*(xPi+xL);
  Double_t vy = 0.5*(yPi+yL);
  Double_t vz = 0.5*(zPi+zL);

  Double_t vertx = -1.*vx;
  Double_t verty = vz;
  Double_t vertz = vy + tpc::ZTarget;
  return TVector3(vertx, verty, vertz);
}

//_____________________________________________________________________________
TVector3 CalcCloseDistXi(TVector3 point, Double_t Bfield,
			 TVector3 xi_decayvtx, TVector3 xi_mom_decayvtx,
			 TVector3 &xi_mom, Double_t &dist){

  Double_t xi_par[5];
  CalcHelixParam(Bfield, -1, xi_mom_decayvtx, xi_decayvtx, xi_par);

  Double_t xi = -1.*point.x();
  Double_t yi = point.z() - tpc::ZTarget;
  Double_t zi = point.y();

  //helix function
  //x = [0] + [3]*cos(t);
  //y = [1] + [3]*sin(t);
  //z = [2] + [3]*[4]*t;

  TF1 fvertex_xi("fvertex_xi", "pow([0]+[3]*cos(x)-[5], 2)+pow([1]+[3]*sin(x)-[6], 2)+pow([2]+[3]*[4]*x-[7], 2)", -TMath::Pi(), TMath::Pi());
  fvertex_xi.SetParameter(0, xi_par[0]);
  fvertex_xi.SetParameter(1, xi_par[1]);
  fvertex_xi.SetParameter(2, xi_par[2]);
  fvertex_xi.SetParameter(3, xi_par[3]);
  fvertex_xi.SetParameter(4, xi_par[4]);
  fvertex_xi.SetParameter(5, xi);
  fvertex_xi.SetParameter(6, yi);
  fvertex_xi.SetParameter(7, zi);

  Double_t helix_t = fvertex_xi.GetMinimumX(-TMath::Pi(), TMath::Pi());
  dist = TMath::Sqrt(fvertex_xi.GetMinimum());
  xi_mom = CalcHelixMom(Bfield, -1, xi_par, helix_t);

  Double_t vx = xi_par[0]+xi_par[3]*cos(helix_t);
  Double_t vy = xi_par[1]+xi_par[3]*sin(helix_t);
  Double_t vz = xi_par[2]+xi_par[3]*xi_par[4]*helix_t;

  Double_t vertx = -1.*vx;
  Double_t verty = vz;
  Double_t vertz = vy + tpc::ZTarget;

  return TVector3(vertx, verty, vertz);
}

//_____________________________________________________________________________
TVector3
CloseDistLambda(TVector3 point, TVector3 Xlambda,
		TVector3 Plambda, Double_t &dist){

  Double_t lambdavtx_xivtx_cut = 0.;

  Double_t local_xi = -1.*point.x();
  Double_t local_yi = point.z() - tpc::ZTarget;
  Double_t local_zi = point.y();

  Double_t xi = -1.*Xlambda.x();
  Double_t yi = Xlambda.z() - tpc::ZTarget;
  Double_t zi = Xlambda.y();
  Double_t pxi = -1.*Plambda.x();
  Double_t pyi = Plambda.z();
  Double_t pzi = Plambda.y();
  Double_t ui = -pxi/pyi, vi = pzi/pyi;

  TVector3 p_L = TVector3(pxi, pyi, pzi);
  TVector3 p_unit = p_L.Unit();

  Double_t scan_range[2] ={-250. - tpc::ZTarget, 250. - tpc::ZTarget};
  if(pyi>0) scan_range[1] = yi + lambdavtx_xivtx_cut/(p_unit.y());
  else scan_range[0] = yi - lambdavtx_xivtx_cut/(p_unit.y());

  TF1 fpoint_linear("fpoint_linear", "pow([0]-([3]+[4]*x), 2)+pow([1]-x, 2)+pow([2]-([5]+[6]*x), 2)", scan_range[0], scan_range[1]);
  fpoint_linear.SetParameter(0, local_xi);
  fpoint_linear.SetParameter(1, local_yi);
  fpoint_linear.SetParameter(2, local_zi);
  fpoint_linear.SetParameter(3, xi + ui*yi);
  fpoint_linear.SetParameter(4, -ui);
  fpoint_linear.SetParameter(5, zi - vi*yi);
  fpoint_linear.SetParameter(6, vi);

  //straight function
  //x = [3] + [4]*y;
  //z = [5] + [6]*y;

  Double_t close_y = fpoint_linear.GetMinimumX(scan_range[0], scan_range[1]);
  dist = TMath::Sqrt(fpoint_linear.Eval(close_y));

  Double_t xL = xi - ui*(close_y-yi);
  Double_t yL = close_y;
  Double_t zL = zi + vi*(close_y-yi);
  Double_t vertx = -xL;
  Double_t verty = zL;
  Double_t vertz = yL + tpc::ZTarget;
  return TVector3(vertx, verty, vertz);
}

//_____________________________________________________________________________
TVector3
CalcCloseDistLambda(TVector3 point, TVector3 Xlambda,
		    TVector3 Plambda, Double_t &dist){

  Double_t lambdavtx_xivtx_cut = 0.;

  Double_t local_xi = -1.*point.x();
  Double_t local_yi = point.z() - tpc::ZTarget;
  Double_t local_zi = point.y();

  Double_t xi = -1.*Xlambda.x();
  Double_t yi = Xlambda.z() - tpc::ZTarget;
  Double_t zi = Xlambda.y();
  Double_t pxi = -1.*Plambda.x();
  Double_t pyi = Plambda.z();
  Double_t pzi = Plambda.y();
  Double_t ui = -pxi/pyi, vi = pzi/pyi;

  TVector3 p_L = TVector3(pxi, pyi, pzi);
  TVector3 p_unit = p_L.Unit();

  Double_t scan_range[2] ={-250. - tpc::ZTarget, 250. - tpc::ZTarget};
  if(pyi>0) scan_range[1] = yi + lambdavtx_xivtx_cut/(p_unit.y());
  else scan_range[0] = yi - lambdavtx_xivtx_cut/(p_unit.y());

  TF1 fpoint_linear("fpoint_linear", "pow([0]-([3]+[4]*x), 2)+pow([1]-x, 2)+pow([2]-([5]+[6]*x), 2)", scan_range[0], scan_range[1]);
  fpoint_linear.SetParameter(0, local_xi);
  fpoint_linear.SetParameter(1, local_yi);
  fpoint_linear.SetParameter(2, local_zi);
  fpoint_linear.SetParameter(3, xi + ui*yi);
  fpoint_linear.SetParameter(4, -ui);
  fpoint_linear.SetParameter(5, zi - vi*yi);
  fpoint_linear.SetParameter(6, vi);

  //straight function
  //x = [3] + [4]*y;
  //z = [5] + [6]*y;

  Double_t close_y = fpoint_linear.GetMinimumX(scan_range[0], scan_range[1]);
  dist = TMath::Sqrt(fpoint_linear.Eval(close_y));

  Double_t xL = xi - ui*(close_y-yi);
  Double_t yL = close_y;
  Double_t zL = zi + vi*(close_y-yi);
  Double_t vertx = -xL;
  Double_t verty = zL;
  Double_t vertz = yL + tpc::ZTarget;
  return TVector3(vertx, verty, vertz);
}

//_____________________________________________________________________________
TVector3
LambdaTargetCenter(TVector3 Xlambda, TVector3 Plambda,
		   Double_t &dist){

  Double_t u = Plambda.x()/Plambda.z();
  Double_t v = Plambda.y()/Plambda.z();
  Double_t zdiff = tpc::ZTarget - Xlambda.z();
  TVector3 track(u*zdiff, v*zdiff, zdiff);
  TVector3 lpos = track + Xlambda;
  TVector3 diff(lpos.x(), lpos.y(), lpos.z() - tpc::ZTarget);
  dist = diff.Mag();
  return TVector3(lpos.x(), lpos.y(), lpos.z());
}

//_____________________________________________________________________________
TVector3
LambdaLambdaVertex(TVector3 Xlambda1, TVector3 Plambda1,
		   TVector3 Xlambda2, TVector3 Plambda2,
		   TVector3 &vtxlambda1, TVector3 &vtxlambda2,
		   Double_t &dist){

  Double_t lambdavtx_xivtx_cut = 0.;

  Double_t xi1 = -1.*Xlambda1.x();
  Double_t yi1 = Xlambda1.z() - tpc::ZTarget;
  Double_t zi1 = Xlambda1.y();
  Double_t pxi1 = -1.*Plambda1.x();
  Double_t pyi1 = Plambda1.z();
  Double_t pzi1 = Plambda1.y();
  Double_t ui1 = -pxi1/pyi1, vi1 = pzi1/pyi1;

  TVector3 p_L1 = TVector3(pxi1, pyi1, pzi1);
  TVector3 p_unit1 = p_L1.Unit();

  Double_t scan_range1[2] ={-250. - tpc::ZTarget, 250. - tpc::ZTarget};
  if(pyi1>0) scan_range1[1] = yi1 + lambdavtx_xivtx_cut/(p_unit1.y());
  else scan_range1[0] = yi1 - lambdavtx_xivtx_cut/(p_unit1.y());

  Double_t xi2 = -2.*Xlambda2.x();
  Double_t yi2 = Xlambda2.z() - tpc::ZTarget;
  Double_t zi2 = Xlambda2.y();
  Double_t pxi2 = -2.*Plambda2.x();
  Double_t pyi2 = Plambda2.z();
  Double_t pzi2 = Plambda2.y();
  Double_t ui2 = -pxi2/pyi2, vi2 = pzi2/pyi2;

  TVector3 p_L2 = TVector3(pxi2, pyi2, pzi2);
  TVector3 p_unit2 = p_L2.Unit();

  Double_t scan_range2[2] ={-250. - tpc::ZTarget, 250. - tpc::ZTarget};
  if(pyi2>0) scan_range2[1] = yi2 + lambdavtx_xivtx_cut/(p_unit2.y());
  else scan_range2[0] = yi2 - lambdavtx_xivtx_cut/(p_unit2.y());

  TF2 fpoint_ll("fpoint_ll", "pow(([0]+[1]*x)-([4]+[5]*y), 2)+pow((x)-(y), 2)+pow(([2]+[3]*x)-([6]+[7]*y), 2)", scan_range1[0], scan_range1[1], scan_range2[0], scan_range2[1]);
  fpoint_ll.SetParameter(0, xi1 + ui1*yi1);
  fpoint_ll.SetParameter(1, -ui1);
  fpoint_ll.SetParameter(2, zi1 - vi1*yi1);
  fpoint_ll.SetParameter(3, vi1);
  fpoint_ll.SetParameter(4, xi2 + ui2*yi2);
  fpoint_ll.SetParameter(5, -ui2);
  fpoint_ll.SetParameter(6, zi2 - vi2*yi2);
  fpoint_ll.SetParameter(7, vi2);

  Double_t close_y1, close_y2;
  fpoint_ll.GetMinimumXY(close_y1, close_y2);
  dist = TMath::Sqrt(fpoint_ll.Eval(close_y1, close_y2));

  Double_t xL1 = xi1 - ui1*(close_y1-yi1);
  Double_t yL1 = close_y1;
  Double_t zL1 = zi1 + vi1*(close_y1-yi1);
  vtxlambda1 = TVector3(-xL1, zL1, yL1 + tpc::ZTarget);

  Double_t xL2 = xi2 - ui2*(close_y2-yi2);
  Double_t yL2 = close_y2;
  Double_t zL2 = zi2 + vi2*(close_y2-yi2);
  vtxlambda2 = TVector3(-xL2, zL2, yL2 + tpc::ZTarget);

  Double_t vertx = -0.5*(xL1 + xL2);
  Double_t verty = 0.5*(zL1 + zL2);
  Double_t vertz = 0.5*(yL1 + yL2) + tpc::ZTarget;
  return TVector3(vertx, verty, vertz);
}

//_____________________________________________________________________________
Bool_t
HelixDirection(TVector3 vertex, TVector3 start, TVector3 end, Double_t &dist){

  Bool_t status = false;
  TVector3 dist1 = vertex - start;
  TVector3 dist2 = vertex - end;
  dist = dist1.Mag();
  if(dist1.Mag() < dist2.Mag()) status = true;
  return status;
}

//_____________________________________________________________________________
TVector3
MultitrackVertex(Int_t ntrack, Double_t *x0, Double_t *y0, Double_t *u0, Double_t *v0,
		 std::vector<Double_t> Res_x0, std::vector<Double_t> Res_y0,
		 std::vector<Double_t> Res_u0, std::vector<Double_t> Res_v0){

  gNumOfTracks = ntrack;
  gX0.clear();
  gY0.clear();
  gU0.clear();
  gV0.clear();
  gSX0.clear();
  gSY0.clear();
  gSU0.clear();
  gSV0.clear();
  for(Int_t i=0; i<gNumOfTracks; ++i){
    gX0.push_back(x0[i]);
    gY0.push_back(y0[i]);
    gU0.push_back(u0[i]);
    gV0.push_back(v0[i]);
  }
  bool ResFlag = true;
  if(Res_x0.size() == 0){
    ResFlag = false;
  }
  else if(Res_x0.size() != ntrack || Res_y0.size() != ntrack
	  || Res_u0.size() != ntrack || Res_v0.size() != ntrack){
    std::cout<<"MultitrackVertex:: Resolution not set"<<std::endl;
    std::cout<<"ntracks = "<<ntrack<<std::endl;
    std::cout<<"Res_x0.size() = "<<Res_x0.size()<<std::endl;
    std::cout<<"Res_y0.size() = "<<Res_y0.size()<<std::endl;
    std::cout<<"Res_u0.size() = "<<Res_u0.size()<<std::endl;
    std::cout<<"Res_v0.size() = "<<Res_v0.size()<<std::endl;
    ResFlag = false;
  }
  if(ResFlag){
    for(Int_t i=0; i<gNumOfTracks; ++i){
      gSX0.push_back(Res_x0[i]);
      gSY0.push_back(Res_y0[i]);
      gSU0.push_back(Res_u0[i]);
      gSV0.push_back(Res_v0[i]);
    }
  }
  else{
    for(Int_t i=0;i<gNumOfTracks;++i){
      gSX0.push_back(1.);
      gSY0.push_back(1.);
      gSU0.push_back(0.);
      gSV0.push_back(0.);
    }
  }

  Double_t par[3] = {0, 0, 0};
  Double_t err[3] = {999., 999., 999.};

  TMinuit *minuit = new TMinuit(3);
  minuit->SetPrintLevel(-1);
  minuit->SetFCN(fcn_vertex);

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1; //error level for ch2 minimization
  minuit->mnexcm("SET ERR", arglist, 1, ierflg);
  minuit->mnexcm("SET NOW", arglist, 1, ierflg);

  TString name[3] = {"x", "y" ,"z"};
  const Double_t FitStep[3] = {0.001, 0.001, 0.001};
  const Double_t LowLimit[3] = {-50., -50., -100};
  const Double_t UpLimit[3] = {50., 50., 100};
  for(Int_t i=0; i<3; i++){
    minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }
  minuit->Command("SET STRategy 0");

  //arglist[0] = 500.;
  //arglist[1] = 1;
  arglist[0] = 1000.;
  arglist[1] = 0.1;

  Int_t Err;
  Double_t bnd1, bnd2;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  //minuit->mnimpr();
  //minuit->mnexcm("MINOS", arglist, 0, ierflg);
  //minuit->mnexcm("SET ERR", arglist, 2, ierflg);

  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  for(Int_t i=0; i<3; i++){
    minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
  }
  delete minuit;

  return TVector3(par[0], par[1], par[2] + tpc::ZTarget);
}

//_____________________________________________________________________________
TVector3
MultitrackVertex(Int_t ntrack, Double_t *x0, Double_t *y0, Double_t *u0, Double_t *v0,
		 std::vector<Double_t> Res_x0, std::vector<Double_t> Res_y0,
		 std::vector<Double_t> Res_u0, std::vector<Double_t> Res_v0, Double_t &chisqr){

  gNumOfTracks = ntrack;
  gX0.clear();
  gY0.clear();
  gU0.clear();
  gV0.clear();
  gSX0.clear();
  gSY0.clear();
  gSU0.clear();
  gSV0.clear();
  for(Int_t i=0; i<gNumOfTracks; ++i){
    gX0.push_back(x0[i]);
    gY0.push_back(y0[i]);
    gU0.push_back(u0[i]);
    gV0.push_back(v0[i]);
  }
  bool ResFlag = true;
  if(Res_x0.size() == 0){
    ResFlag = false;
  }
  else if(Res_x0.size() != ntrack || Res_y0.size() != ntrack
	  || Res_u0.size() != ntrack || Res_v0.size() != ntrack){
    std::cout<<"MultitrackVertex:: Resolution not set"<<std::endl;
    std::cout<<"ntracks = "<<ntrack<<std::endl;
    std::cout<<"Res_x0.size() = "<<Res_x0.size()<<std::endl;
    std::cout<<"Res_y0.size() = "<<Res_y0.size()<<std::endl;
    std::cout<<"Res_u0.size() = "<<Res_u0.size()<<std::endl;
    std::cout<<"Res_v0.size() = "<<Res_v0.size()<<std::endl;
    ResFlag = false;
  }
  if(ResFlag){
    for(Int_t i=0; i<gNumOfTracks; ++i){
      gSX0.push_back(Res_x0[i]);
      gSY0.push_back(Res_y0[i]);
      gSU0.push_back(Res_u0[i]);
      gSV0.push_back(Res_v0[i]);
    }
  }
  else{
    for(Int_t i=0;i<gNumOfTracks;++i){
      gSX0.push_back(1.);
      gSY0.push_back(1.);
      gSU0.push_back(0.);
      gSV0.push_back(0.);
    }
  }

  Double_t par[3] = {0, 0, 0};
  Double_t err[3] = {999., 999., 999.};

  TMinuit *minuit = new TMinuit(3);
  minuit->SetPrintLevel(-1);
  minuit->SetFCN(fcn_vertex);

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1; //error level for ch2 minimization
  minuit->mnexcm("SET ERR", arglist, 1, ierflg);
  minuit->mnexcm("SET NOW", arglist, 1, ierflg);

  TString name[3] = {"x", "y" ,"z"};
  const Double_t FitStep[3] = {0.001, 0.001, 0.001};
  const Double_t LowLimit[3] = {-50., -50., -100};
  const Double_t UpLimit[3] = {50., 50., 100};
  for(Int_t i=0; i<3; i++){
    minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }
  minuit->Command("SET STRategy 0");

  //arglist[0] = 500.;
  //arglist[1] = 1;
  arglist[0] = 1000.;
  arglist[1] = 0.1;

  Int_t Err;
  Double_t bnd1, bnd2;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  //minuit->mnimpr();
  //minuit->mnexcm("MINOS", arglist, 0, ierflg);
  //minuit->mnexcm("SET ERR", arglist, 2, ierflg);

  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  for(Int_t i=0; i<3; i++){
    minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
  }
  Double_t grad[3]; Double_t chi;
  minuit -> Eval(3, grad, chi, par, 0);
  chisqr = chi;

  delete minuit;

  return TVector3(par[0], par[1], par[2] + tpc::ZTarget);
}

}

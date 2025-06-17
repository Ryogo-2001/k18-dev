// -*- C++ -*-

#include "TPCLTrackHit.hh"

#include <string>
#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <TString.h>
#include <TF1.h>

#include <std_ostream.hh>

#include "ConfMan.hh"
#include "FuncName.hh"
#include "TPCAnalyzer.hh"
#include "MathTools.hh"
#include "TPCPadHelper.hh"
//#include "TPCParamMan.hh"

namespace
{
const Double_t& HSfield_Calib = ConfMan::Get<Double_t>("HSFLDCALIB");
const Double_t& HSfield_Calc = ConfMan::Get<Double_t>("HSFLDCALC");
const Double_t& HSfield_Hall = ConfMan::Get<Double_t>("HSFLDHALL");
}

//______________________________________________________________________________
TPCLTrackHit::TPCLTrackHit(TPCHit *hit)
  : m_hit(hit),
    m_layer(hit->GetLayer()),
    m_section(-1),
    m_mrow(hit->GetMRow()),
    m_local_hit_pos(hit->GetPosition()),
    m_cal_pos(TVector3(0.,0.,0.)),
    m_res(),
    m_res_param(hit->GetResolutionParams()),
    m_x0(TMath::QuietNaN()),
    m_y0(TMath::QuietNaN()),
    m_u0(TMath::QuietNaN()),
    m_v0(TMath::QuietNaN()),
    m_cx(TMath::QuietNaN()),
    m_cy(TMath::QuietNaN()),
    m_z0(TMath::QuietNaN()),
    m_r(TMath::QuietNaN()),
    m_dz(TMath::QuietNaN()),
    m_t(TMath::QuietNaN()),
    m_padtheta(hit->GetPadTheta()),
    m_padlength(hit->GetPadLength()),
    m_de(hit->GetCDe()),

    m_cal_pos_exclusive(TVector3(0.,0.,0.)),
    m_x0_exclusive(TMath::QuietNaN()),
    m_y0_exclusive(TMath::QuietNaN()),
    m_u0_exclusive(TMath::QuietNaN()),
    m_v0_exclusive(TMath::QuietNaN()),
    m_cx_exclusive(TMath::QuietNaN()),
    m_cy_exclusive(TMath::QuietNaN()),
    m_z0_exclusive(TMath::QuietNaN()),
    m_r_exclusive(TMath::QuietNaN()),
    m_dz_exclusive(TMath::QuietNaN()),
    m_t_exclusive(TMath::QuietNaN())
{

  if((m_local_hit_pos.z() + m_local_hit_pos.x()) < 0 && (m_local_hit_pos.z() - m_local_hit_pos.x()) < 0) m_section = 1;
  if((m_local_hit_pos.z() + m_local_hit_pos.x()) < 0 && (m_local_hit_pos.z() - m_local_hit_pos.x()) > 0) m_section = 2;
  if((m_local_hit_pos.z() + m_local_hit_pos.x()) > 0 && (m_local_hit_pos.z() - m_local_hit_pos.x()) > 0) m_section = 3;
  if((m_local_hit_pos.z() + m_local_hit_pos.x()) > 0 && (m_local_hit_pos.z() - m_local_hit_pos.x()) < 0) m_section = 4;

  m_hit->RegisterHits(this);
  debug::ObjectCounter::increase(ClassName());
}

//______________________________________________________________________________
TPCLTrackHit::TPCLTrackHit(const TPCLTrackHit& right)
  : m_hit(right.m_hit),
    m_layer(right.m_layer),
    m_section(-1),
    m_mrow(right.m_mrow),
    m_local_hit_pos(right.m_local_hit_pos),
    m_cal_pos(right.m_cal_pos),
    m_res(right.m_res),
    m_res_param(right.m_res_param),
    m_x0(right.m_x0),
    m_y0(right.m_y0),
    m_u0(right.m_u0),
    m_v0(right.m_v0),
    m_cx(right.m_cx),
    m_cy(right.m_cy),
    m_z0(right.m_z0),
    m_r(right.m_r),
    m_dz(right.m_dz),
    m_t(right.m_t),
    m_padtheta(right.m_padtheta),
    m_padlength(right.m_padlength),
    m_de(right.m_de),

    m_cal_pos_exclusive(right.m_cal_pos_exclusive),
    m_x0_exclusive(right.m_x0_exclusive),
    m_y0_exclusive(right.m_y0_exclusive),
    m_u0_exclusive(right.m_u0_exclusive),
    m_v0_exclusive(right.m_v0_exclusive),
    m_cx_exclusive(right.m_cx_exclusive),
    m_cy_exclusive(right.m_cy_exclusive),
    m_z0_exclusive(right.m_z0_exclusive),
    m_r_exclusive(right.m_r_exclusive),
    m_dz_exclusive(right.m_dz_exclusive),
    m_t_exclusive(right.m_t_exclusive)
{
  m_hit->RegisterHits(this);
  debug::ObjectCounter::increase(ClassName());
}

//______________________________________________________________________________
TPCLTrackHit::~TPCLTrackHit()
{
  debug::ObjectCounter::decrease(ClassName());
}

// For Helix fit
//______________________________________________________________________________
TVector3
TPCLTrackHit::GetHelixPosition(double par[5], double t) const
{
  //This is the eqation of Helix
  double x = par[0] + par[3]*cos(t);
  double y = par[1] + par[3]*sin(t);
  double z = par[2] + (par[4]*par[3]*t);

  return TVector3(x, y, z);
}

//______________________________________________________________________________
TVector3
TPCLTrackHit::GetLocalCalPos() const
{
  TVector3 pos = m_local_hit_pos;
  TVector3 x0(m_x0, m_y0, tpc::ZTarget);
  TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, tpc::ZTarget+1.);
  TVector3 u = (x1-x0).Unit();
  TVector3 AP = pos-x0;
  double dist_AX = u.Dot(AP);
  TVector3 AI(x0.x()+(u.x()*dist_AX),
	      x0.y()+(u.y()*dist_AX),
	      x0.z()+(u.z()*dist_AX));
  return AI;
}

//______________________________________________________________________________
TVector3
TPCLTrackHit::GetLocalCalPosExclusive() const
{
  TVector3 pos = m_local_hit_pos;
  TVector3 x0(m_x0_exclusive, m_y0_exclusive, tpc::ZTarget);
  TVector3 x1(m_x0_exclusive + m_u0_exclusive, m_y0_exclusive + m_v0_exclusive, tpc::ZTarget+1.);
  TVector3 u = (x1-x0).Unit();
  TVector3 AP = pos-x0;
  double dist_AX = u.Dot(AP);
  TVector3 AI(x0.x()+(u.x()*dist_AX),
	      x0.y()+(u.y()*dist_AX),
	      x0.z()+(u.z()*dist_AX));
  return AI;
}

//______________________________________________________________________________
TVector3
TPCLTrackHit::GetLocalCalPosHelix() const
{
  double par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 fittmp = GetHelixPosition(par, m_t);
  TVector3 calpos(-fittmp.X(),
		   fittmp.Z(),
		   fittmp.Y()+tpc::ZTarget);
  return calpos;
}

//______________________________________________________________________________
TVector3
TPCLTrackHit::GetLocalCalPosHelixExclusive() const
{
  double par[5] = {m_cx_exclusive, m_cy_exclusive, m_z0_exclusive, m_r_exclusive, m_dz_exclusive};
  TVector3 fittmp = GetHelixPosition(par, m_t_exclusive);
  TVector3 calpos(-fittmp.X(),
		   fittmp.Z(),
		   fittmp.Y()+tpc::ZTarget);
  return calpos;
}

//______________________________________________________________________________
TVector3
TPCLTrackHit::GetMomentumHelix(Double_t charge) const
{
  TVector3 pos(-m_cal_pos.X(),
   	       m_cal_pos.Z() - tpc::ZTarget,
   	       m_cal_pos.Y());

  const double Const = 0.299792458; // =c/10^9
  const double dMagneticField = HSfield_Calib*(HSfield_Hall/HSfield_Calc);

  double t = (pos.Z()-m_z0)/(m_r*m_dz);
  double pt = fabs(m_r)*(Const*dMagneticField); // MeV/c
  //From here!!!!
  double tmp_px = pt*(-1.*sin(t));
  double tmp_py = pt*(cos(t));
  double tmp_pz = pt*(m_dz);
  double px = -tmp_px*0.001;
  double py = tmp_pz*0.001;
  double pz = tmp_py*0.001;

  TVector3 p = TVector3(px,py,pz);
  if(charge<0.) p *= -1.;
  return p;
}

//______________________________________________________________________________
TVector3
TPCLTrackHit::GetResidualVect() const
{
  return m_local_hit_pos - m_cal_pos;
}

//______________________________________________________________________________
TVector3
TPCLTrackHit::GetResidualVectExclusive() const
{
  return m_local_hit_pos - m_cal_pos_exclusive;
}

//______________________________________________________________________________
double
TPCLTrackHit::GetResidual() const
{
  TVector3 Res = m_local_hit_pos - m_cal_pos;
  return Res.Mag();
}

//______________________________________________________________________________
double
TPCLTrackHit::GetResidualExclusive() const
{
  TVector3 Res = m_local_hit_pos - m_cal_pos_exclusive;
  return Res.Mag();
}

//______________________________________________________________________________
Double_t
TPCLTrackHit::GetPadTrackAngleHelix() const
{
  Double_t thetaDiff = m_t - m_padtheta;
  return thetaDiff;
}

//______________________________________________________________________________
Double_t
TPCLTrackHit::GetPathHelix() const
{
  //Approximation of pathlength
  Double_t thetaDiff = GetPadTrackAngleHelix();
  Double_t path_xz = m_padlength/TMath::Abs(cos(thetaDiff));
  Double_t factor = TMath::Sqrt(1.+(TMath::Power(m_dz,2)));
  return path_xz*factor;
}

//______________________________________________________________________________
Bool_t
TPCLTrackHit::IsGoodForTracking() const
{
  Bool_t isgood = true;
  if(m_res.x() > 0.9e+10 && m_res.y() > 0.9e+10 && m_res.z() > 0.9e+10) isgood = false;
  return isgood;
}

//______________________________________________________________________________
void
TPCLTrackHit::Print(const TString& arg) const
{
  const int w = 2;
  std::cout << arg.Data() << " Hough flag" << m_hit->GetHoughFlag()
	    << " L" << m_hit->GetLayer() << " R" << m_hit->GetMRow()
	    << std::setw(w) << std::right << " pos"
	    << std::fixed << std::setprecision(1) << std::setw(w) << m_local_hit_pos
	    << std::endl;

}

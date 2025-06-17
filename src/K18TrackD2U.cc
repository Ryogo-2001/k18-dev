// -*- C++ -*-

#include "K18TrackD2U.hh"

#include <algorithm>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <std_ostream.hh>

#include "ConfMan.hh"
#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DCLocalTrack.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "K18Parameters.hh"
#include "K18TransMatrix.hh"
#include "Minuit.hh"
#include "TrackHit.hh"

namespace
{
using namespace K18Parameter;
const auto& gGeom   = DCGeomMan::GetInstance();
const auto& gK18Mtx = K18TransMatrix::GetInstance();
const Double_t& zK18Target = gGeom.LocalZ("K18Target");
const Double_t LowBand[5] =
{ MinK18InX, MinK18InY, MinK18InU, MinK18InV, MinK18Delta };
const Double_t UpperBand[5] =
{ MaxK18InX, MaxK18InY, MaxK18InU, MaxK18InV, MaxK18Delta };
const Double_t LowBandOut[5] =
{ MinK18OutX, MinK18OutY, MinK18OutU, MinK18OutV, MinK18Delta };
const Double_t UpperBandOut[5] =
{ MaxK18OutX, MaxK18OutY, MaxK18OutU, MaxK18OutV, MaxK18Delta };
}

//_____________________________________________________________________________
K18TrackD2U::K18TrackD2U(Double_t local_x, DCLocalTrack* track_out,
                         Double_t p0)
  : m_local_x(local_x),
    m_track_out(track_out),
    m_p0(p0),
    m_Xi(0),
    m_Yi(0),
    m_Ui(0),
    m_Vi(0),
    m_Xo(0),
    m_Yo(0),
    m_Uo(0),
    m_Vo(0),
    m_status(false),
    m_delta(),
    m_delta3rd(),
    m_good_for_analysis(true)
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
K18TrackD2U::~K18TrackD2U()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
Bool_t
K18TrackD2U::CalcMomentumD2U()
{
  m_status = false;

  if(!gK18Mtx.IsReady()) return m_status;

  Double_t xi, yi, ui, vi;
  Double_t xo, yo, uo, vo;
  Double_t delta1, delta2;
  xi= -m_track_out->GetX0();
  yi=  m_track_out->GetY0();
  ui=  m_track_out->GetU0();
  vi= -m_track_out->GetV0();
  xo= -m_local_x;

  m_status = gK18Mtx.CalcDeltaD2U(xi, yi, ui, vi,
                                  xo, yo, uo, vo,
                                  delta1, delta2);

#if 0
  hddaq::cout << FUNC_NAME << ": after calculation. "
	      << " StatusD2U=" << m_status  << std::endl;
#endif

  if(m_status){
    // hddaq::cout << "delta1 = " << delta1 << ", delta2 = " << delta2 << std::endl;
    m_delta    = delta1;
    m_delta3rd = delta2;
    m_Xi = -xo;
    m_Yi =  yo;
    m_Ui =  uo;
    m_Vi = -vo;
    m_Xo = m_track_out->GetX0();
    m_Yo = m_track_out->GetY0();
    m_Uo = m_track_out->GetU0();
    m_Vo = m_track_out->GetV0();
  }

  return m_status;
}

//_____________________________________________________________________________
TVector3
K18TrackD2U::BeamMomentumD2U() const
{
  Double_t u  = m_track_out->GetU0(), v=m_track_out->GetV0();
  Double_t pz = P()/TMath::Sqrt(1.+u*u+v*v);
  return TVector3(pz*u, pz*v, pz);
}

//_____________________________________________________________________________
Double_t
K18TrackD2U::GetChiSquare() const
{
  return m_track_out->GetChiSquare();
}

//_____________________________________________________________________________
Double_t
K18TrackD2U::Xtgt() const
{
  return m_track_out->GetU0()*zK18Target+m_track_out->GetX0();
}

//_____________________________________________________________________________
Double_t
K18TrackD2U::Ytgt() const
{
  return m_track_out->GetV0()*zK18Target+m_track_out->GetY0();
}

//_____________________________________________________________________________
Double_t
K18TrackD2U::Utgt() const
{
  return m_track_out->GetU0();
}

//_____________________________________________________________________________
Double_t
K18TrackD2U::Vtgt() const
{
  return m_track_out->GetV0();
}

//_____________________________________________________________________________
Bool_t
K18TrackD2U::GoodForAnalysis(Bool_t status)
{
  Bool_t pre_status = m_good_for_analysis;
  m_good_for_analysis = status;
  return pre_status;
}

//_____________________________________________________________________________
Bool_t
K18TrackD2U::ReCalc(Bool_t applyRecursively)
{
  if(applyRecursively) m_track_out->ReCalc(applyRecursively);
  return CalcMomentumD2U();
}

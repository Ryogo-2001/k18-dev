// -*- C++ -*-

#include "TrackHit.hh"

#include <cstring>
#include <sstream>
#include <stdexcept>

#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "DebugCounter.hh"
#include "FuncName.hh"

//_____________________________________________________________________________
TrackHit::TrackHit(DCLTrackHit *hit)
  : m_dcltrack_hit(hit)
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TrackHit::~TrackHit()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
Double_t
TrackHit::GetLocalHitPos() const
{
  return m_dcltrack_hit->GetLocalHitPos();
}

//_____________________________________________________________________________
Double_t
TrackHit::GetResidual() const
{
  Double_t a = GetTiltAngle()*TMath::DegToRad();
  Double_t u = m_cal_global_mom.x()/m_cal_global_mom.z();
  Double_t v = m_cal_global_mom.y()/m_cal_global_mom.z();
  Double_t dsdz = u*TMath::Cos(a)+v*TMath::Sin(a);
  Double_t coss = IsHoneycomb() ? TMath::Cos(TMath::ATan(dsdz)) : 1.;
  Double_t wp = GetWirePosition();
  Double_t ss = wp+(GetLocalHitPos()-wp)/coss;
  return (ss-m_cal_local_pos)*coss;
}

//_____________________________________________________________________________
Double_t
TrackHit::GetTiltAngle() const
{
  return m_dcltrack_hit->GetTiltAngle();
}

//_____________________________________________________________________________
Bool_t
TrackHit::ReCalc(Bool_t applyRecursively)
{
  if(applyRecursively)
    return m_dcltrack_hit->ReCalc(applyRecursively);
  else
    return true;
}

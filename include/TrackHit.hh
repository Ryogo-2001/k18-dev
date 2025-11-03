// -*- C++ -*-

#ifndef TRACK_HIT_HH
#define TRACK_HIT_HH

#include "DCLTrackHit.hh"

#include <TVector3.h>

class DCLocalTrack;

//______________________________________________________________________________
class TrackHit
{
public:
  static const TString& ClassName();
  explicit TrackHit(DCLTrackHit *hit);
  ~TrackHit();

private:
  TrackHit(const TrackHit&);
  TrackHit& operator =(const TrackHit&);

private:
  DCLTrackHit* m_dcltrack_hit;
  TVector3     m_cal_global_mom;
  TVector3     m_cal_global_pos;
  Double_t     m_cal_local_pos;
  Double_t     m_cal_local_u;
  Double_t     m_cal_local_v;

public:
  Bool_t          IsHoneycomb() const { return m_dcltrack_hit->IsHoneycomb(); }
  void            SetCalGMom(const TVector3 &mom) { m_cal_global_mom = mom; }
  void            SetCalGPos(const TVector3 &pos) { m_cal_global_pos = pos; }
  void            SetCalLPos(Double_t pos) { m_cal_local_pos=pos; }
  void            SetCalLUV(Double_t u, Double_t v) { m_cal_local_u = u; m_cal_local_v = v; }
  DCLTrackHit*    GetHit() const { return m_dcltrack_hit; }
  Double_t        GetWirePosition() const
  { return m_dcltrack_hit->GetWirePosition(); }
  int             GetLayer() const { return m_dcltrack_hit->GetLayer(); }
  Double_t        GetLocalHitPos() const;
  const TVector3& GetCalGPos() const { return m_cal_global_pos; }
  Double_t        GetCalLPos() const { return m_cal_local_pos; }
  Double_t        GetResidual() const;
  Double_t        GetResolution() const
  { return m_dcltrack_hit->GetResolution(); }
  Double_t        GetTiltAngle() const;
  DCLTrackHit*    GetDCLTrack(){ return m_dcltrack_hit; }
  Bool_t          ReCalc(Bool_t applyRecursively=false);
};

//______________________________________________________________________________
inline const TString&
TrackHit::ClassName()
{
  static TString s_name("TrackHit");
  return s_name;
}

#endif

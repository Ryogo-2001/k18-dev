// -*- C++ -*-

#ifndef DC_LTRACK_HIT_HH
#define DC_LTRACK_HIT_HH

#include "DCHit.hh"

#include <TString.h>

class DCAnalyzer;

//_____________________________________________________________________________
class DCLTrackHit
{
public:
  static const TString& ClassName();
  DCLTrackHit(DCHit* hit, Double_t pos, Int_t nh);
  DCLTrackHit(const DCLTrackHit& right);
  ~DCLTrackHit();

private:
  DCHit*   m_hit;
  Int_t    m_nth_hit;
  Double_t m_local_hit_pos;
  Double_t m_cal_pos; //for VXU
  Double_t m_xcal;
  Double_t m_ycal;
  Double_t m_ucal;
  Double_t m_vcal;
  Bool_t   m_honeycomb;
  //exlcusive
  //double  m_cal_pos_exclusive; //for VXU, not supported
  Bool_t  m_is_fitted_exclusive;
  double  m_xcal_exclusive;
  double  m_ycal_exclusive;
  double  m_ucal_exclusive;
  double  m_vcal_exclusive;

public:
  void     SetLocalHitPos(Double_t xl) { m_local_hit_pos = xl; }
  void     SetCalPosition(Double_t x, Double_t y) { m_xcal = x; m_ycal = y; }
  void     SetCalUV(Double_t u, Double_t v) { m_ucal = u; m_vcal = v; }
  void     SetHoneycomb(Bool_t flag=true) { m_honeycomb = flag; }
  Bool_t   IsHoneycomb() const { return m_honeycomb; }
  Int_t    GetLayer() const { return m_hit->GetLayer(); }
  Int_t    GetMeanSeg() const { return m_hit->GetMeanSeg(); }
  Int_t    GetMaxSeg() const { return m_hit->GetMaxSeg(); }
  Double_t GetWire() const { return m_hit->GetWire(); }
  Int_t    GetTdcVal() const { return m_hit->GetTdcVal(m_nth_hit); }
  Int_t    GetTdcSize() const { return m_hit->GetTdcSize(); }
  Double_t GetDriftTime() const { return m_hit->GetDriftTime(m_nth_hit); }
  void     SetDriftTime(Int_t ith, Double_t dt)
  { m_hit->SetDriftTime(ith, dt); }
  Double_t GetDriftLength() const { return m_hit->GetDriftLength(m_nth_hit); }
  void     SetDriftLength(Int_t ith, Double_t dl)
  { m_hit->SetDriftLength(ith, dl); }
  Double_t GetTrailingTime() const
  { return m_hit->GetTrailingTime(m_nth_hit); }
  Double_t GetTot() const { return m_hit->GetTot(m_nth_hit); }
  DCHit*   GetHit() const { return m_hit; }
  Double_t GetTiltAngle() const { return m_hit->GetTiltAngle(); }
  Double_t GetWirePosition() const { return m_hit->GetWirePosition(); }
  Double_t GetLocalHitPos() const { return m_local_hit_pos; }
  Double_t GetLocalCalPos() const;
  Double_t GetXcal() const { return m_xcal; }
  Double_t GetYcal() const { return m_ycal; }
  Double_t GetUcal() const { return m_ucal; }
  Double_t GetVcal() const { return m_vcal; }
  Double_t GetResidual() const;
  Double_t GetResolution() const { return m_hit->GetResolution(); }

  ///// for XUV tracking
  void     SetLocalCalPosVXU(Double_t xcl) { m_cal_pos=xcl; }
  Double_t GetLocalCalPosVXU() const { return m_cal_pos; }
  Double_t GetResidualVXU() const { return m_local_hit_pos-m_cal_pos; }

  //exlcusive
  void     SetExclusiveReady(Bool_t flag=true) { m_is_fitted_exclusive = flag; }
  void     SetCalPositionExclusive(Double_t x, Double_t y) { m_xcal_exclusive = x; m_ycal_exclusive = y; }
  void     SetCalUVExclusive(Double_t u, Double_t v) { m_ucal_exclusive = u; m_vcal_exclusive = v; }
  Double_t GetLocalCalPosExclusive()  const;
  Double_t GetResidualExclusive() const;

  ///// for TOF
  Double_t GetZ() const { return m_hit->GetZ(); }

  ///// for CFT
  // Double_t GetPositionR() const { return m_hit->GetPositionR(); }
  // Double_t GetPositionPhi() const { return m_hit->GetPositionPhi(); }
  // Double_t GetPosPhi() const { return m_hit->GetPosPhi(); }
  // Double_t GetPosZ() const { return m_hit->GetPosZ(); }
  // Double_t GetPosR() const { return m_hit->GetPosR(); }
  // Double_t GetAdcLow() const { return m_hit->GetAdcLow(); }
  // Double_t GetMIPLow() const { return m_hit->GetMIPLow(); }
  // Double_t GetdELow() const { return m_hit->GetdELow(); }
  // Double_t GetMaxAdcLow() const { return m_hit->GetMaxAdcLow(); }
  // Double_t GetMaxMIPLow() const { return m_hit->GetMaxMIPLow(); }
  // Double_t GetMaxdELow() const { return m_hit->GetMaxdELow(); }
  // TVector3 GetVtx() const { return m_hit->GetVtx(); }
  // Double_t GetTime() const { return m_hit->GetTime(); }

  void   JoinTrack() { m_hit->JoinTrack(m_nth_hit); }
  void   QuitTrack() { m_hit->QuitTrack(m_nth_hit); }
  Bool_t BelongToTrack() const { return m_hit->BelongToTrack(m_nth_hit); }

  // void   JoinTrackCFT() { m_hit->JoinTrackCFT(m_nth_hit); }
  // void   QuitTrackCFT() { m_hit->QuitTrackCFT(m_nth_hit); }
  // Bool_t BelongToTrackCFT() const
  // { return m_hit->BelongToTrackCFT(m_nth_hit); }

  void Print(const TString& arg="") const;

  Bool_t ReCalc(Bool_t applyRecursively=false);

  friend class DCHit;
};

//_____________________________________________________________________________
inline const TString&
DCLTrackHit::ClassName()
{
  static TString s_name("DCLTrackHit");
  return s_name;
}

#endif

// -*- C++ -*-

#ifndef TPC_LTRACK_HIT_HH
#define TPC_LTRACK_HIT_HH

#include <TVector3.h>

#include "DCHit.hh"
#include "MathTools.hh"
#include "ThreeVector.hh"
#include "TPCHit.hh"
#include "TPCCluster.hh"

class TPCAnalyzer;

//_____________________________________________________________________________
class TPCLTrackHit
{
public:
  static const TString& ClassName();
  TPCLTrackHit(TPCHit *hit);
  TPCLTrackHit(const TPCLTrackHit& right);
  ~TPCLTrackHit();

private:
  TPCHit* m_hit;
  Int_t m_layer;
  Int_t m_section;
  Double_t m_mrow;
  TVector3 m_local_hit_pos;
  TVector3 m_cal_pos;
  TVector3 m_res;
  std::vector<Double_t> m_res_param;
  Double_t m_x0;
  Double_t m_y0;
  Double_t m_u0;
  Double_t m_v0;

  Double_t m_cx;
  Double_t m_cy;
  Double_t m_z0;
  Double_t m_r;
  Double_t m_dz;
  Double_t m_t;

  Double_t m_padtheta;
  Double_t m_padlength;
  Double_t m_de;

  TVector3 m_cal_pos_exclusive;
  Double_t m_x0_exclusive;
  Double_t m_y0_exclusive;
  Double_t m_u0_exclusive;
  Double_t m_v0_exclusive;
  Double_t m_cx_exclusive;
  Double_t m_cy_exclusive;
  Double_t m_z0_exclusive;
  Double_t m_r_exclusive;
  Double_t m_dz_exclusive;
  Double_t m_t_exclusive;

public:
  Bool_t IsGood() const { return (m_hit && m_hit->IsGood()); }
  Bool_t IsGoodForTracking() const;
  Bool_t IsOnTheFrame(){ return m_hit->GetParentCluster()->IsOnTheFrame(); };
  void   SetLocalHitPos(TVector3 xl) { m_local_hit_pos = xl; }
  void   SetCalPosition(TVector3 cal_pos) { m_cal_pos = cal_pos; }
  void   SetCalX0Y0(Double_t x0, Double_t y0) { m_x0 = x0; m_y0 = y0; }
  void   SetCalUV(Double_t u0, Double_t v0) { m_u0 = u0; m_v0 = v0; }
  void   SetCalHelix(Double_t cx, Double_t cy, Double_t z0, Double_t r, Double_t dz)
  {m_cx = cx; m_cy = cy; m_z0 = z0; m_r = r; m_dz = dz;}
  void   SetResolution(TVector3 res) { m_res = res; }
  void   SetTheta(Double_t t) { m_t = t; }

  TPCHit* GetHit() const { return m_hit; }

  const TVector3& GetLocalHitPos()  const { return m_local_hit_pos; }
  Int_t  GetLayer()       const { return m_layer; }
  Int_t  GetSection()       const { return m_section; }
  Double_t GetMRow()      const { return m_mrow; }
  Double_t GetPadTheta()  const { return m_padtheta; }
  Double_t GetPadLength() const { return m_padlength; }
  Double_t GetDe()        const { return m_de; }
  Double_t GetXcal()      const { return m_cal_pos.x(); }
  Double_t GetYcal()      const { return m_cal_pos.y(); }
  Double_t GetZcal()      const { return m_cal_pos.z(); }
  Double_t GetUcal()      const { return m_u0; }
  Double_t GetVcal()      const { return m_v0; }
  Double_t Getcx()        const { return m_cx; }
  Double_t Getcy()        const { return m_cy; }
  Double_t Getz0()        const { return m_z0; }
  Double_t Getr()         const { return m_r; }
  Double_t Getdz()        const { return m_dz; }
  Double_t GetTheta()     const { return m_t; }
  const TVector3& GetResolutionVect() const { return m_res; }
  const std::vector<Double_t>& GetResolutionParams() const { return m_res_param; }

  TVector3 GetLocalCalPos() const;
  TVector3 GetLocalCalPosHelix() const;
  TVector3 GetHelixPosition(Double_t par[5], Double_t t) const;
  TVector3 GetMomentumHelix(Double_t charge)  const;
  Double_t GetPadTrackAngleHelix() const;
  Double_t GetPathHelix() const;
  TVector3 GetResidualVect() const;
  Double_t GetResidual() const;
  Double_t GetResolution() const { return m_res.Mag(); }

  Int_t GetHoughFlag() const { return m_hit->GetHoughFlag(); }
  Double_t GetHoughDist() const { return m_hit->GetHoughDist(); }
  Double_t GetHoughDistY() const { return m_hit->GetHoughDistY(); }
  void Print(const TString& arg="") const;

  TVector3 GetLocalCalPosExclusive() const;
  TVector3 GetLocalCalPosHelixExclusive() const;
  void SetDe(Double_t de){ m_de = de; }
  void SetCalX0Y0Exclusive(Double_t x0, Double_t y0) { m_x0_exclusive = x0; m_y0_exclusive = y0; }
  void SetCalUVExclusive(Double_t u0, Double_t v0) { m_u0_exclusive = u0; m_v0_exclusive = v0; }
  void SetThetaExclusive(Double_t t) { m_t_exclusive = t; }
  void SetCalPositionExclusive(TVector3 cal_pos) { m_cal_pos_exclusive = cal_pos; }
  void SetCalHelixExclusive(Double_t cx, Double_t cy, Double_t z0, Double_t r, Double_t dz)
  {m_cx_exclusive = cx; m_cy_exclusive = cy; m_z0_exclusive = z0; m_r_exclusive = r, m_dz_exclusive = dz;}
  TVector3 GetResidualVectExclusive() const;
  Double_t GetResidualExclusive() const;

  friend class TPCHit;
};

//_____________________________________________________________________________
inline const TString&
TPCLTrackHit::ClassName()
{
  static TString s_name("TPCLTrackHit");
  return s_name;
}

#endif

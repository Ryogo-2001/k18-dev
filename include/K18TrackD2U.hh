// -*- C++ -*-

#ifndef K18_TRACK_D2U_HH
#define K18_TRACK_D2U_HH

#include <vector>
#include <functional>
#include <TString.h>

#include <TVector3.h>

class DCLocalTrack;
class TrackHit;
class DCAnalyzer;
class DCHit;

//_____________________________________________________________________________
class K18TrackD2U
{
public:
  static const TString& ClassName();
  K18TrackD2U(Double_t local_x, DCLocalTrack* track_out, Double_t p0);
  ~K18TrackD2U();

private:
  K18TrackD2U();
  K18TrackD2U(const K18TrackD2U &);
  K18TrackD2U& operator= (const K18TrackD2U &);

private:
  Double_t      m_local_x;
  DCLocalTrack* m_track_out;
  Double_t      m_p0;
  Double_t      m_Xi;
  Double_t      m_Yi;
  Double_t      m_Ui;
  Double_t      m_Vi;
  Double_t      m_Xo;
  Double_t      m_Yo;
  Double_t      m_Uo;
  Double_t      m_Vo;
  Bool_t        m_status;
  Double_t      m_delta;
  Double_t      m_delta3rd;
  Bool_t        m_good_for_analysis;

public:
  TVector3      BeamMomentumD2U() const;
  Bool_t        CalcMomentumD2U();
  Double_t      Delta() const { return m_delta; }
  Double_t      Delta3rd() const { return m_delta3rd; }
  Double_t      GetChiSquare() const;
  Bool_t        GoodForAnalysis(Bool_t status);
  Bool_t        GoodForAnalysis() const { return m_good_for_analysis; }
  Double_t      P() const { return m_p0*(1.+m_delta); }
  Double_t      P3rd() const { return m_p0*(1.+m_delta3rd); }
  Bool_t        ReCalc(Bool_t applyRecursively=false);
  Bool_t        StatusD2U() const { return m_status; }
  DCLocalTrack* TrackOut() { return m_track_out; }
  Double_t      Uin() const { return m_Ui; }
  Double_t      Uout() const { return m_Uo; }
  Double_t      Utgt() const;
  Double_t      Vin() const { return m_Vi; }
  Double_t      Vout() const { return m_Vo; }
  Double_t      Vtgt() const;
  Double_t      Xin() const { return m_Xi; }
  Double_t      Xout() const { return m_Xo; }
  Double_t      Xtgt() const;
  Double_t      Yin() const { return m_Yi; }
  Double_t      Yout() const { return m_Yo; }
  Double_t      Ytgt() const;
};

//_____________________________________________________________________________
inline const TString&
K18TrackD2U::ClassName()
{
  static TString s_name("K18TrackD2U");
  return s_name;
}

#endif

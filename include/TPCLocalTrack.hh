// -*- C++ -*-

#ifndef TPC_LOCAL_TRACK_HH
#define TPC_LOCAL_TRACK_HH

#include <vector>
#include <functional>

#include <TH2D.h>
#include <TMinuit.h>
#include <TString.h>
#include <TVector3.h>

#include <std_ostream.hh>

#include "DetectorID.hh"
#include "TPCHit.hh"
#include "TPCCluster.hh"
#include "TPCLTrackHit.hh"

class TPCHit;
class TPCCluster;

//_____________________________________________________________________________
class TPCLocalTrack
{
public:
  static const TString& ClassName();
  explicit TPCLocalTrack();
  ~TPCLocalTrack();
  TPCLocalTrack(TPCLocalTrack *init); //deep copy

private:
  TPCLocalTrack & operator =(const TPCLocalTrack &init);

private:
  Bool_t   m_is_fitted;     // flag of DoFit()
  Bool_t   m_is_calculated; // flag of Calculate()
  std::vector<TPCLTrackHit*> m_hit_array;

  //track coordinate origin is center of the target, ***NOT TPC center***
  //Track params
  Double_t m_x0;
  Double_t m_y0;
  Double_t m_u0;
  Double_t m_v0;
  Int_t m_isAccidental;
  TVector3 m_closedist; //closest distance from the track to the target
  TVector3 m_edgepoint; //the most outer point of the track
  Double_t m_chisqr;
  Int_t m_minuit; //Minuit output status 0:not calculated at all 1:approximation only, not accurate
  //2:full matrix, but forced positive-definite 3:full accurate covariance matrix
  Int_t m_n_iteration;
  Int_t m_fitflag;
  Int_t m_vtxflag; //for track seperation around the target
  Int_t m_searchtime; //millisec
  Int_t m_fittime; //millisec
  Int_t m_trackid; //for k18, kurama track

  //exclusive
  std::vector<Double_t> m_x0_exclusive;
  std::vector<Double_t> m_y0_exclusive;
  std::vector<Double_t> m_u0_exclusive;
  std::vector<Double_t> m_v0_exclusive;
  std::vector<Double_t> m_chisqr_exclusive;

public:

  void          AddTPCHit(TPCLTrackHit *hit);
  void          EraseHits(std::vector<Int_t> delete_hits);
  void          EraseHit(Int_t delete_hit);
  void          ClearHits();
  void          Calculate();
  void          CalcClosestDistTgt();
  void          DeleteNullHit();
  Bool_t        DoFit(Int_t MinHits=0);
  Bool_t        DoStraightTrackFit();
  TVector3      CalcResidual(TVector3 position);
  Double_t      GetHorizontalResidual(TPCHit *hit);
  Double_t      GetHorizontalResidual(Int_t i);
  Double_t      GetVerticalResidual(TPCHit *hit);
  Double_t      GetVerticalResidual(Int_t i);

  TVector3      GetResolutionVect(Int_t i, Bool_t vetoBadClusters);
  Double_t      GetHorizontalResolution(TPCHit *hit);
  Double_t      GetHorizontalResolution(Int_t i);
  Double_t      GetVerticalResolution(TPCHit *hit);
  Double_t      GetVerticalResolution(Int_t i);

  Bool_t        ResidualCheck(Int_t i, Double_t &residual);
  Bool_t        IsGoodHitToAdd(TPCHit *hit, Double_t &residual);
  Int_t         Side(TVector3 pos);
  void          Print(const TString& arg="", Bool_t print_allhits = false) const;
  Bool_t        IsGoodForTracking();

  TPCLTrackHit* GetHit(Int_t i) const { return m_hit_array.at(i); }
  Int_t         GetNHit() const { return m_hit_array.size(); }

  Double_t GetChiSquare() const { return m_chisqr; }
  Double_t GetClosestDist() const { return m_closedist.Mag(); }
  Int_t    GetMinuitStatus() const { return m_minuit; }
  Int_t    GetNIteration() const { return m_n_iteration; }
  Int_t    GetFitFlag(void) const {return m_fitflag;}
  Int_t    GetVtxFlag(void) const {return m_vtxflag;} //vtx is in the target or not
  Bool_t   IsFitted() const { return m_is_fitted; }
  Bool_t   IsCalculated() const { return m_is_calculated; }
  Int_t    GetSearchTime() const { return m_searchtime; }
  Int_t    GetFitTime() const { return m_fittime; }
  void     GetParam(Double_t *par) const { par[0] = m_x0; par[1] = m_y0; par[2] = m_u0; par[3] = m_v0; }

  Double_t GetX0() const { return m_x0; }
  Double_t GetY0() const { return m_y0; }
  Double_t GetU0() const { return m_u0; }
  Double_t GetV0() const { return m_v0; }

  Int_t    GetNPad() const;
  Int_t    GetNDF() const;
  TVector3 GetPosition(Double_t z) const;
  Double_t GetTheta() const;
  Double_t GetAlpha(TPCHit *hit) const; //alpha : pad - track angle
  Double_t GetAlpha(Int_t i) const; //alpha : pad - track angle
  Double_t GetHitLength(Int_t i) const;
  Double_t GetTrackdE();

  void SetParam(Double_t *par){ m_x0 = par[0]; m_y0 = par[1]; m_u0 = par[2]; m_v0 = par[3]; }
  void SetClustersHoughFlag(Int_t hough_flag);
  void SetTrackId(Int_t flag) { m_trackid = flag; } //Dummy function for the template(E42 is not using this)
  void SetFlag(Int_t flag) { return; } //Dummy function for the template(E42 is not using this)
  void FinalizeTrack() { return; } //Dummy function for the template (E42 is not using this)
  bool DetermineCharge() { return false; } //Dummy function for the template (E42 is not using this)
  Bool_t TestInvertCharge() { return false; } //Dummy function for the template (E42 is not using this)
  void SetFitFlag(Int_t flag) { m_fitflag = flag; }
  void SetVtxFlag(Int_t flag) { m_vtxflag = flag; }
  void SetSearchTime(Int_t time) { m_searchtime = time; }
  void SetFitTime(Int_t time) { m_fittime = time; }

  Bool_t VertexAtTarget();
  Bool_t IsBackward();
  void CheckIsAccidental();
  Bool_t SeparateTracksAtTarget();
  Bool_t SeparateClustersWithGap();
  void RecalcTrack();

  //exclusive
  void CalculateExclusive();
  void DoFitExclusive();
  Double_t GetHorizontalResidualExclusive(Int_t i);
  Double_t GetVerticalResidualExclusive(Int_t i);
  void GetParamExclusive(Int_t i, Double_t *par) const { par[0] = m_x0_exclusive[i]; par[1] = m_y0_exclusive[i]; par[2] = m_u0_exclusive[i]; par[3] = m_v0_exclusive[i]; }
  Double_t GetX0Exclusive(Int_t i) const { return m_x0_exclusive[i]; }
  Double_t GetY0Exclusive(Int_t i) const { return m_y0_exclusive[i]; }
  Double_t GetU0Exclusive(Int_t i) const { return m_u0_exclusive[i]; }
  Double_t GetV0Exclusive(Int_t i) const { return m_v0_exclusive[i]; }

};

//_____________________________________________________________________________
inline const TString&
TPCLocalTrack::ClassName()
{
  static TString s_name("TPCLocalTrack");
  return s_name;
}

#endif

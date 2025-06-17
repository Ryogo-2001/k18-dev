// -*- C++ -*-

#ifndef TPC_LOCAL_TRACK_HELIX_HH
#define TPC_LOCAL_TRACK_HELIX_HH

#include <vector>
#include <functional>

#include <std_ostream.hh>

#include "TVector3.h"
#include "ThreeVector.hh"
#include "DetectorID.hh"
#include "TPCHit.hh"
#include "TPCCluster.hh"
#include "TPCLTrackHit.hh"
#include "TMatrixD.h"
class TPCHit;
class TPCCluster;

//______________________________________________________________________________
class TPCLocalTrackHelix
{
public:
  static const TString ClassName();
  explicit TPCLocalTrackHelix();
  ~TPCLocalTrackHelix();
  TPCLocalTrackHelix(TPCLocalTrackHelix *init); //deep copy

private:
  TPCLocalTrackHelix & operator =(const TPCLocalTrackHelix &init);

private:
  Bool_t m_is_fitted;     // flag of DoFit()
  Bool_t m_is_calculated; // flag of Calculate()
  Bool_t m_is_theta_calculated; // flag of CalcHelixTheta()
  Bool_t m_is_fitted_exclusive; // flag of DoFitExclusive()
  //Bool_t m_is_thetaflip; //Atan2 theta flip occurs or not
  Bool_t m_is_multiloop; // Multi-loop track w/ High p_L, low p_T

  std::vector<TPCLTrackHit*> m_hit_array;
  std::vector<Int_t> m_hit_order;
  std::vector<Double_t> m_hit_t; //theta
  std::vector<Int_t> m_kuramaid_candidate; //for kurama track

  //equation of Helix
  //x = -X;
  //y = Z - Tgtz;
  //z = Y;
  //x = p[0] + p[3]*cos(theta);
  //y = p[1] + p[3]*sin(theta);
  //z = p[2] + p[4]*p[3]*(theta);
  // track coordinate origin is target, ***NOT TPC center***
  //Track params
  Double_t m_cx;
  Double_t m_cy;
  Double_t m_z0;
  Double_t m_r;
  Double_t m_dz;

  Int_t m_pid;
  TVector3 m_closedist; //closest distance from the track to the target
  TVector3 m_closedistXZ; //closest distance from the track to the target XZ position
  Double_t m_chisqr;
  Int_t m_minuit; //Minuit output status 0:not calculated at all 1:approximation only, not accurate
  //2:full matrix, but forced positive-definite 3:full accurate covariance matrix
  Int_t m_n_iteration;
  TVector3 m_mom0;
  TVector3 m_edgepoint; //the most outer point of the track
  Double_t m_min_t;
  Double_t m_max_t;
  Double_t m_path;
  Double_t m_transverse_path;
  Int_t m_charge;
  Int_t m_fitflag;
  Int_t m_vtxflag; //for track seperation around the target
  Int_t m_isBeam;
  Int_t m_isK18;
  Int_t m_isKurama;
  Int_t m_isAccidental;
  Int_t m_isXi; // xi track
  Int_t m_trackid; //for k18, kurama track
  Int_t m_ncl_beforetgt; //for k18 track
  Int_t m_searchtime; //millisec
  Int_t m_fittime; //millisec

  //scaling factor applied to momentum resolution
  //common for all tracks
  Double_t m_MomResScale;
  Double_t m_dZResScale;
  Double_t m_PhResScale;

  //exclusive
  std::vector<Double_t> m_cx_exclusive;
  std::vector<Double_t> m_cy_exclusive;
  std::vector<Double_t> m_z0_exclusive;
  std::vector<Double_t> m_r_exclusive;
  std::vector<Double_t> m_dz_exclusive;
  std::vector<Double_t> m_chisqr_exclusive;
  std::vector<Double_t> m_t_exclusive;
  std::vector<TVector3> m_vp; //RK virtual plane

public:

  void         AddTPCHit(TPCLTrackHit *hit);
  void         EraseHits(std::vector<Int_t> delete_hits);
  void         EraseHit(Int_t delete_hit);
  void         ClearHits();
  void         Calculate();
  void         CalcClosestDistTgt();
  void         CalcHelixTheta();
  void         DeleteNullHit();
  Bool_t       DoFit(Int_t MinHits=0);
  Bool_t       DoHelixTrackFit();
  Bool_t       DoStraightLineFit(Double_t *par);
  Bool_t       DoCircleFit(Double_t *par);
  Bool_t       DoPreFit(Double_t *par);
  Bool_t       DoHelixFit(Double_t *par, Bool_t vetoBadClusters);

  Bool_t       DetermineCharge();
  Int_t        FinalizeTrack(Int_t &delete_hit);
  void         SortHitOrder(); //Sort hits by theta

  Bool_t       ResidualCheck(Int_t i, Double_t &residual);
  Bool_t       IsGoodHitToAdd(TPCHit *hit, Double_t &residual, Bool_t nolimitation=false);
  Int_t        Side(TVector3 hitpos);
  void         Print(const TString& arg="", Bool_t print_allhits=false) const;
  Bool_t       IsGoodForTracking();

  TVector3 CalcHelixMom(Double_t par[5], Double_t theta) const;
  TVector3 CalcHelixMomCenter(Double_t par[5]) const;

  TPCLTrackHit* GetHit(std::size_t nth) const;
  TPCLTrackHit* GetHitInOrder(std::size_t nth) const;
  Int_t         GetNHit() const { return m_hit_array.size(); }
  Int_t         GetNHitsEffective() const;
  Int_t         GetOrder(Int_t i) const;

  Double_t   GetChiSquare() const { return m_chisqr; }
  TVector3   GetClosestDistVect() const { return m_closedist; }
  Double_t   GetClosestDist() const { return m_closedist.Mag(); }
  Double_t   GetClosestDistXZ() const { return m_closedistXZ.Mag(); }
  TVector3   GetClosestPositionTgt();
  TVector3   GetClosestPositionTgtXZ();

  Int_t      GetMinuitStatus() const { return m_minuit; }
  Int_t      GetPid() const { return m_pid; }
  Int_t      GetNIteration() const { return m_n_iteration; }
  Int_t      GetFitFlag() const { return m_fitflag; }
  Int_t      GetVtxFlag() const { return m_vtxflag; } //vtx is in the target or not
  Bool_t     IsFitted() const { return m_is_fitted; }
  Bool_t     IsCalculated() const { return m_is_calculated; }
  Bool_t     IsThetaCalculated() const { return m_is_theta_calculated; }
  Int_t      GetSearchTime() const { return m_searchtime; }
  Int_t      GetFitTime() const { return m_fittime; }

  TVector3 GetMom0() const { return m_mom0; }// Momentum at Y = 0

  //TVector3 GetResolutionVect(TPCHit *hit, Bool_t vetoBadClusters);
  TVector3 GetResolutionVect(TPCLTrackHit *hit, Bool_t vetoBadClusters);
  TVector3 GetResolutionVect(Int_t i, Bool_t vetoBadClusters);
  Double_t GetResolutionY(TPCHit *hit);
  Double_t GetResolutionY(Int_t i);

  TVector3 GetMomentumResolutionVectT(Double_t t, Double_t MomScale, Double_t PhiScale, Double_t dZScale);
  TVector3 GetMomentumResolutionVect(Int_t i, Double_t MomScale = 1, Double_t PhiScale = 1, Double_t dZScale= 1);//These scale factors are for scaling the resolution of individual track.
  TVector3 GetMomentumResolutionVect();
  TVector3 GetMomentumCovarianceVectT(Double_t t, Double_t MomScale, Double_t PhiScale, Double_t dZScale);
  TVector3 GetMomentumCovarianceVect(Int_t i, Double_t MomScale = 1, Double_t PhiScale = 1, Double_t dZScale = 1);
  TVector3 GetMomentumCovarianceVect();
  Double_t GetMomentumResolution();
  Double_t GetTransverseMomentumAngularCovariance(Double_t t = -9999);
  Double_t GetTransverseMomentumResolution();//returns dP, not dP/P;
  Double_t GetTransverseAngularResolution(Double_t t, Double_t sig0 = 0.01); //returns angular resolution on pad plane
  Double_t GetTransverseAngularResolution();
  TMatrixD GetCovarianceMatrix();

  Double_t GetdZResolution();//returns pitch resolution.
  Double_t GetThetaResolution();//returns pitch angle resolution.
  Double_t GetMomentumPitchAngleCovariance();

  TVector3 CalcResidual(TVector3 position);

  void GetParam(Double_t *par) const { par[0] = m_cx; par[1] = m_cy; par[2] = m_z0; par[3] = m_r; par[4] = m_dz; }
  Double_t Getcx() const { return m_cx; }
  Double_t Getcy() const { return m_cy; }
  Double_t Getz0() const { return m_z0; }
  Double_t Getr() const { return m_r; }
  Double_t Getdz() const { return m_dz; }
  Double_t GetTheta(Int_t i) const { return m_hit_t[i]; }
  TVector3 GetVPPos(Int_t i) const { return m_vp[i]; }
  Double_t GetMint() const { return m_min_t; }
  Double_t GetMaxt() const { return m_max_t; }
  Double_t GetPath() const { return m_path; }
  Double_t GetTransversePath() const { return m_transverse_path; }
  Int_t    GetCharge() const { return m_charge; }
  Double_t GetTrackdE();
  Double_t GetdEdx(Double_t truncatedMean = 1.0);

  Int_t    GetNPad() const;
  Int_t    GetNDF() const;
  TVector3 GetPosition(Double_t par[5], Double_t t) const;
  Double_t GetAlpha(Int_t i) const; //alpha : pad - track angle
  Double_t GetAlpha(TPCHit *hit) const; //alpha : pad - track angle

  Int_t GetIsBeam() const { return m_isBeam; }
  Int_t GetIsK18() const { return m_isK18; }
  Int_t GetIsKurama() const { return m_isKurama; }
  Int_t GetIsAccidental() const { return m_isAccidental; }
  Int_t GetIsMultiloop() const { return m_is_multiloop; }
  Int_t GetIsXi() const { return m_isXi; }

  void SetParam(Double_t *par){ m_cx = par[0]; m_cy = par[1]; m_z0 = par[2]; m_r = par[3]; m_dz = par[4]; }
  void SetClustersHoughFlag(Int_t hough_flag);

  void SetFlag(Int_t flag);
  void SetCharge(Int_t flag) { m_charge = flag; }
  void SetFitFlag(Int_t flag) { m_fitflag = flag; }
  void SetVtxFlag(Int_t flag) { m_vtxflag = flag; }
  void SetIsBeam(Int_t flag=1) { m_isBeam = flag; }
  void SetIsK18(Int_t flag=1) { m_isK18 = flag; }
  void SetIsKurama(Int_t flag=1) { m_isKurama = flag; }
  void SetIsAccidental(Int_t flag=1) { m_isAccidental = flag; }
  void SetIsXi(Int_t flag=1) { m_isXi = flag; }
  void SetIsFitted(Bool_t flag=true) { m_is_fitted = flag; }
  void SetIsCalculated(Bool_t flag=true) { m_is_calculated = flag; }
  void SetIsThetaCalculated(Bool_t flag=true) { m_is_theta_calculated = flag; }
  void SetSearchTime(Int_t time) { m_searchtime = time; }
  void SetFitTime(Int_t time) { m_fittime = time; }
  //void SetIsThetaFlip() { m_is_thetaflip = true; }
  void SetMint(Double_t t) { m_min_t = t; }
  void SetMaxt(Double_t t) { m_max_t = t; }
  Bool_t ConvertParam(Double_t *linear_par);
  void AdjustResolutionScale(Double_t MomScale=1, Double_t PhiScale=1, Double_t dZScale=1){
    m_MomResScale*=MomScale; m_PhResScale*=PhiScale; m_dZResScale*=dZScale;
  };

  Bool_t IsBackward();
  void IsMultiLoop();
  void CheckIsAccidental();
  Bool_t VertexAtTarget();
  Bool_t SeparateTracksAtTarget();
  Bool_t SeparateClustersWithGap();
  Bool_t TestMergedTrack();
  Bool_t TestInvertCharge();
  void RecalcTrack();

  //for K1.8 & Kurama tracking
  Bool_t       DoCircleFitwMomConstraint(Double_t *par);
  Bool_t       DoHelixTrackFit(Double_t RKpar[5]);
  Bool_t       DoFit(Double_t RKpar[5], Int_t MinHits=0); //momentum constraint
  Bool_t       DoFit(Double_t RKCharge, Int_t MinHits); //charge constraint
  Bool_t       DoFit(Double_t RKCharge, Double_t RKpar[5], Int_t MinHits=0); //charge&&momoentum constraint
  Bool_t       ResidualCheck(TVector3 pos, Double_t xzwindow, Double_t ywindow, Double_t &resi);
  Bool_t       ResidualCheck(TVector3 pos, Double_t xzwindow, Double_t ywindow);
  void         AddVPHit(TVector3 vp);
  Bool_t       DoVPFit();
  Int_t        GetVPNHit() const { return m_vp.size(); }

  //For the Kurama&TPC Runge-Kutta tracking
  Int_t GetNclBeforeTgt() const { return m_ncl_beforetgt; }
  void SetNclBeforeTgt(Int_t ncl) { m_ncl_beforetgt = ncl; }
  void SetTrackID(Int_t flag) { m_trackid = flag; }
  Int_t GetTrackID() { return m_trackid; }
  void AddTrackIDCandidate(Int_t flag) { m_kuramaid_candidate.push_back(flag); }
  Bool_t IsKuramaTrackCandidate(Int_t id);

  //exclusive
  void CalculateExclusive();
  Double_t CalcThetaExclusive(Int_t ith) const;
  void DoFitExclusive();
  Double_t GetHorizontalResidualExclusive(Int_t i);
  Double_t GetVerticalResidualExclusive(Int_t i);
  void GetParamExclusive(Int_t i, Double_t *par) const { par[0] = m_cx_exclusive[i]; par[1] = m_cy_exclusive[i]; par[2] = m_z0_exclusive[i]; par[3] = m_r_exclusive[i]; par[4] = m_dz_exclusive[i]; }
  Double_t GetcxExclusive(Int_t i) const { return m_cx_exclusive[i]; }
  Double_t GetcyExclusive(Int_t i) const { return m_cy_exclusive[i]; }
  Double_t Getz0Exclusive(Int_t i) const { return m_z0_exclusive[i]; }
  Double_t GetrExclusive(Int_t i) const { return m_r_exclusive[i]; }
  Double_t GetdzExclusive(Int_t i) const { return m_dz_exclusive[i]; }

  //Re-fit with the vertex constraint
  Bool_t DoFitTrackwVertex(TVector3 vertex_pos, TVector3 vertex_res);

};

//_____________________________________________________________________________
inline const TString
TPCLocalTrackHelix::ClassName()
{
  static TString s_name("TPCLocalTrackHelix");
  return s_name;
}

#endif

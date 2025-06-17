// -*- C++ -*-

#ifndef RK_TPC_TRACK_HH
#define RK_TPC_TRACK_HH

#include <vector>
#include <iosfwd>
#include <iostream>
#include <functional>

#include <TString.h>
#include <TVector3.h>

#include <std_ostream.hh>

#include "RungeKuttaUtilities.hh"

class DCLocalTrack;
class TrackHit;
class TPCLocalTrackHelix;

//_____________________________________________________________________________
class TPCRKTrack
{
public:
  static const TString& ClassName();
  TPCRKTrack(TPCLocalTrackHelix *track_tpc, std::vector<Int_t> lnum);
  TPCRKTrack(TPCLocalTrackHelix *tpctrack, DCLocalTrack* track_in,
	     DCLocalTrack* track_out, std::vector<Int_t> lnum);
  ~TPCRKTrack();

private:
  TPCRKTrack(const TPCRKTrack&);
  TPCRKTrack& operator =(const TPCRKTrack&);

public:
  enum RKstatus { kInit,
		  kPassed,
		  kExceedMaxPathLength,
		  kExceedMaxStep,
		  kFailedTraceLast,
		  kFailedGuess,
		  kFailedSave,
		  kFatal,
		  nRKstatus };

private:
  TString                 s_status[nRKstatus];
  RKstatus                m_status;
  Int_t                   m_trackid;
  DCLocalTrack*           m_track_in;
  DCLocalTrack*           m_track_out;
  TPCLocalTrackHelix*     m_track_tpc;
  Double_t                m_tof_seg;
  Double_t                m_initial_momentum;
  std::vector<TrackHit*>  m_dchit_array;
  RKHitPointContainer     m_HitPointCont;
  Int_t                   m_n_iteration;
  Int_t                   m_nef_iteration;
  Double_t                m_chisqr;
  Double_t                m_polarity;
  Int_t                   m_pid;
  TVector3                m_primary_position;
  TVector3                m_primary_momentum;
  TVector3                m_last_position;
  TVector3                m_last_momentum;
  Double_t                m_path_length_tof;
  Double_t                m_path_length_tgt;
  Double_t                m_path_length_total;
  TVector3                m_tgt_position;
  TVector3                m_tgt_momentum;
  TVector3                m_tof_pos;
  TVector3                m_tof_mom;
  RKCordParameter         m_cord_param;
  Bool_t                  m_is_good;

public:

  Double_t        ChiSquare() const { return m_chisqr; }
  Bool_t          DoFit(TVector3 initPos, TVector3 initMom, Bool_t U2D);
  Bool_t          DoFit(RKCordParameter initCord, Bool_t U2D);
  Double_t        GetChiSquare() const { return m_chisqr; }
  TrackHit*       GetHit(Int_t nth) const { return m_dchit_array.at(nth); }
  TrackHit*       GetHitOfLayerNumber(Int_t lnum) const;
  Double_t        GetInitialMomentum() const { return m_initial_momentum; }
  DCLocalTrack*   GetLocalTrackIn() { return m_track_in; }
  DCLocalTrack*   GetLocalTrackOut() { return m_track_out; }
  Int_t           GetNHits() const { return m_dchit_array.size(); }
  Bool_t          GetTrajectoryLocalPosition(Int_t layer,
                                             Double_t& path, Double_t& x, Double_t& y) const;
  Bool_t          GetTrajectoryLocalPosition(Int_t layer,
                                             Double_t& x, Double_t& y) const;
  Bool_t          GetTrajectoryResidual(Int_t layer, Double_t& resolution, Double_t& residual) const;
  Bool_t          GetTrajectoryResidualTPC(Int_t clusterId, TVector3& resolution, TVector3& residual) const;
  Bool_t          GetTrajectoryGlobalPosition(Int_t layer, TVector3& global_pos) const;
  Bool_t          GetTrajectoryGlobalPositionTPC(Int_t clusterId, TVector3& global_pos) const;
  Bool_t          GetTrajectoryMomentum(Int_t layer, TVector3& global_mom) const;
  Bool_t          GetTrajectoryMomentumTPC(Int_t clusterId, TVector3& global_mom) const;
  Bool_t          GoodForAnalysis() const { return m_is_good; }
  Bool_t          GoodForAnalysis(Bool_t status)
  { m_is_good = status; return status; }
  Int_t           Niteration() const { return m_n_iteration; }
  void            Print(const TString& arg="", std::ostream& ost=hddaq::cout);
  void            SetInitialMomentum(Double_t p) { m_initial_momentum = p; }
  Double_t        PathLengthToTOF() const { return m_path_length_tof; }
  Double_t        PathLengthToTarget() const { return m_path_length_tgt; }
  Double_t        PathLengthTotal() const { return m_path_length_total; }
  Double_t        Polarity() const { return m_polarity; }
  const TVector3& PrimaryMomentum() const { return m_primary_momentum; }
  Double_t        PrimaryMomMag() const { return m_primary_momentum.Mag(); }
  const TVector3& PrimaryPosition() const { return m_primary_position; }
  const TVector3& EndMomentum() const { return m_last_momentum; }
  Double_t        EndMomMag() const { return m_last_momentum.Mag(); }
  const TVector3& EndPosition() const { return m_last_position; }
  const TVector3& TargetMomentum() const { return m_tgt_momentum; }
  Double_t        TargetMomMag() const { return m_tgt_momentum.Mag(); }
  const TVector3& TargetPosition() const { return m_tgt_position; }
  const TVector3& TofMom() const { return m_tof_mom; }
  const TVector3& TofPos() const { return m_tof_pos; }
  Double_t        TofSeg() const { return m_tof_seg; }
  Bool_t          Status() const { return m_status; }

  void SetPID(Int_t pikp) { m_pid = pikp; }
  void SetTofSeg(Double_t id) { m_tof_seg = id; }
  void AddRKHit(TrackHit *hit) { m_dchit_array.push_back(hit); }
  TPCLocalTrackHelix* GetTPCTrack() const { return m_track_tpc; }
  Int_t GetTPCTrackID() const { return m_track_tpc -> GetTrackID(); }
  void SetTPCTrackID(Int_t id) { m_track_tpc -> SetTrackID(id); }
  Int_t GetTrackID() const { return m_trackid; }
  void SetTrackID(Int_t id) { m_trackid = id; }

private:
  void     ClearHitArray();
  Double_t CalcChiSqr(const RKHitPointContainer& hpCont) const;
  void     Initialize(std::vector<Int_t> lnum);
  void     Initialize(DCLocalTrack* track_in, DCLocalTrack* track_out,
		      std::vector<Int_t> lnum);
  Bool_t   GuessNextParameters(const RKHitPointContainer& hpCont,
                               RKCordParameter& Cord,
                               Double_t& estDeltaChisqr,
                               Double_t& lambdaCri, Double_t dmp=0.) const;
  void     PrintCalcHits(const RKHitPointContainer& hpCont,
                         std::ostream& ost=std::cout) const;
  Bool_t   SaveCalcPosition(const RKHitPointContainer& hpCont);
  Bool_t   SaveTrackParameters(const RKCordParameter& cp);

};

//_____________________________________________________________________________
inline const TString&
TPCRKTrack::ClassName()
{
  static TString s_name("TPCRKTrack");
  return s_name;
}

//_____________________________________________________________________________
struct TPCRKTrackComp
  : public std::binary_function<TPCRKTrack*, TPCRKTrack*, Bool_t>
{
  Bool_t operator()(const TPCRKTrack* const p1,
                    const TPCRKTrack* const p2) const
  {
    Int_t n1=p1->GetNHits(), n2=p2->GetNHits();
    if(n1>n2+1) return true;
    if(n2>n1+1) return false;
    return (p1->ChiSquare())<(p2->ChiSquare());
  }
};

#endif

// -*- C++ -*-

#ifndef KURAMA_TRACK_HH
#define KURAMA_TRACK_HH

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
class DCAnalyzer;
class Hodo2Hit;

//_____________________________________________________________________________
class KuramaTrack
{
public:
  static const TString& ClassName();
  KuramaTrack(DCLocalTrack* track_in, DCLocalTrack* track_out);
  ~KuramaTrack();

private:
  KuramaTrack(const KuramaTrack&);
  KuramaTrack& operator =(const KuramaTrack&);

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
  DCLocalTrack*           m_track_in;
  DCLocalTrack*           m_track_out;
  Double_t                m_tof_seg;
  Double_t                m_initial_momentum;
  std::vector<TrackHit*>  m_hit_array;
  RKHitPointContainer     m_HitPointCont;
  Int_t                   m_n_iteration;
  Int_t                   m_nef_iteration;
  Double_t                m_chisqr;
  Double_t                m_polarity;
  TVector3                m_primary_position;
  TVector3                m_primary_momentum;
  Double_t                m_path_length_tof;
  Double_t                m_path_length_total;
  TVector3                m_tof_pos;
  TVector3                m_tof_mom;
  RKCordParameter         m_cord_param;
  Bool_t                  m_is_good;

public:
  Double_t        ChiSquare() const { return m_chisqr; }
  Bool_t          DoFit();
  Bool_t          DoFit(RKCordParameter iniCord);
  Bool_t          DoFitMinuit();
  Double_t        GetChiSquare() const { return m_chisqr; }
  TrackHit*       GetHit(Int_t nth) const { return m_hit_array.at(nth); }
  TrackHit*       GetHitOfLayerNumber(Int_t lnum) const;
  Double_t        GetInitialMomentum() const { return m_initial_momentum; }
  DCLocalTrack*   GetLocalTrackIn() { return m_track_in; }
  DCLocalTrack*   GetLocalTrackOut() { return m_track_out; }
  Int_t           GetNHits() const { return m_hit_array.size(); }
  Bool_t          GetTrajectoryLocalPosition(Int_t layer,
                                             Double_t& path, Double_t& x, Double_t& y) const;
  Bool_t          GetTrajectoryLocalPosition(Int_t layer,
                                             Double_t& x, Double_t& y) const;
  Bool_t          GoodForAnalysis() const { return m_is_good; }
  Bool_t          GoodForAnalysis(Bool_t status)
  { m_is_good = status; return status; }
  Int_t           Niteration() const { return m_n_iteration; }
  Double_t        PathLengthToTOF() const { return m_path_length_tof; }
  Double_t        PathLengthTotal() const { return m_path_length_total; }
  Double_t        Polarity() const { return m_polarity; }
  const TVector3& PrimaryMomentum() const { return m_primary_momentum; }
  Double_t        PrimaryMomMag() const { return m_primary_momentum.Mag();}
  const TVector3& PrimaryPosition() const { return m_primary_position; }
  void            Print(const TString& arg="", std::ostream& ost=hddaq::cout);
  Bool_t          ReCalc(Bool_t applyRecursively=false);
  void            SetInitialMomentum(Double_t p) { m_initial_momentum = p; }
  Bool_t          Status() const { return m_status; }
  const TVector3& TofMom() const { return m_tof_mom; }
  const TVector3& TofPos() const { return m_tof_pos; }
  Double_t        TofSeg() const { return m_tof_seg; }

private:
  void     ClearHitArray();
  Double_t CalcChiSqr(const RKHitPointContainer& hpCont) const;
  void     FillHitArray();
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
KuramaTrack::ClassName()
{
  static TString s_name("KuramaTrack");
  return s_name;
}

//_____________________________________________________________________________
struct KuramaTrackComp
  : public std::binary_function<KuramaTrack*, KuramaTrack*, Bool_t>
{
  Bool_t operator()(const KuramaTrack* const p1,
                    const KuramaTrack* const p2) const
  {
    Int_t n1=p1->GetNHits(), n2=p2->GetNHits();
    if(n1>n2+1) return true;
    if(n2>n1+1) return false;
    return (p1->ChiSquare())<(p2->ChiSquare());
  }
};

#endif

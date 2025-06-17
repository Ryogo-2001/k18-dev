// -*- C++ -*-

#ifndef TPC_VERTEX_HH
#define TPC_VERTEX_HH

#include <vector>

#include "TPCLocalTrackHelix.hh"
#include "TPCLocalTrack.hh"
#include <TString.h>

class TPCLocalTrack;
class TPCLocalTrackHelix;

//_____________________________________________________________________________
class TPCVertex
{
public:
  static const TString& ClassName();
  TPCVertex(Int_t id1, Int_t id2);
  TPCVertex(TVector3 vertex, std::vector<Int_t> trackid);
  ~TPCVertex();

private:

  Bool_t m_is_calculated; // flag of Calculate() m_vertex;
  Int_t m_is_accidental;
  TVector3 m_vertex;
  Double_t m_angle;
  Double_t m_distance;
  std::vector<Int_t>  m_track_id;
  std::vector<Int_t>  m_track_pid;
  std::vector<Int_t>  m_track_charge;
  std::vector<TVector3> m_track_pos;
  std::vector<TVector3> m_track_mom;
  std::vector<Double_t> m_track_theta;

  //for Lambda reconstruction
  Bool_t m_is_lambda;
  TVector3 m_lambda_vertex;
  Double_t m_lambda_distance;
  Double_t m_lambda_angle;
  std::vector<Double_t> m_lambda_mass;
  std::vector<TVector3> m_lambda_mom;
  std::vector<Int_t> m_proton_id;
  std::vector<TVector3> m_proton_mom;
  std::vector<Int_t> m_pion_id;
  std::vector<TVector3> m_pion_mom;

public:

  void Calculate(TPCLocalTrackHelix* track1,
		 TPCLocalTrackHelix* track2);
  void Calculate(TPCLocalTrack* track1,
		 TPCLocalTrack* track2);
  Bool_t IsCalculated() const { return m_is_calculated; }
  Int_t GetIsAccidental() const { return m_is_accidental; }
  void SetIsAccidental(Int_t flag=1) { m_is_accidental=flag; }

  TVector3 GetVertex() const { return m_vertex; }
  Double_t GetOpeningAngle() const { return m_angle; }
  Double_t GetClosestDist() const { return m_distance; }

  Int_t GetTrackId(Int_t i) const { return m_track_id.at(i); }
  Int_t GetTrackCharge(Int_t i) const { return m_track_charge.at(i); }
  Int_t GetTrackPid(Int_t i) const { return m_track_pid.at(i); }
  TVector3 GetTrackPos(Int_t i) const { return m_track_pos.at(i); }
  TVector3 GetTrackMom(Int_t i) const { return m_track_mom.at(i); }
  Double_t GetTrackTheta(Int_t i) const { return m_track_theta.at(i); }

  //For particle reconstruction
  Bool_t ReconstructLambda(TVector3 vertex, TVector3 mom1, TVector3 mom2, Double_t ppi_distance);
  void ReconstructLambdaWithVertex(TPCLocalTrackHelix* track1, TPCLocalTrackHelix* track2);

  Bool_t GetIsLambda() const { return m_is_lambda; }
  Double_t GetClosestDistLambda() const { return m_lambda_distance; }
  Double_t GetOpeningAngleLambda() const { return m_lambda_angle; }
  Int_t GetNcombiLambda() const { return m_lambda_mass.size(); }
  TVector3 GetVertexLambda(Int_t i) const { return m_lambda_vertex; }
  Double_t GetMassLambda(Int_t i) const { return m_lambda_mass.at(i); }
  TVector3 GetMomLambda(Int_t i) const { return m_lambda_mom.at(i); }
  Int_t GetProtonIdLambda(Int_t i) const { return m_proton_id.at(i); }
  Int_t GetPionIdLambda(Int_t i) const { return m_pion_id.at(i); }
  TVector3 GetProtonMomLambda(Int_t i) const { return m_proton_mom.at(i); }
  TVector3 GetPionMomLambda(Int_t i) const { return m_pion_mom.at(i); }

  //For clustered tracks
  Int_t GetNTracks() const { return m_track_id.size(); }

  void Print(const TString& arg="", Bool_t print_all=false) const;
};

//_____________________________________________________________________________
inline const TString&
TPCVertex::ClassName()
{
  static TString s_name("TPCVertex");
  return s_name;
}

#endif

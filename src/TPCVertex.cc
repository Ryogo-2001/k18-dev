// -*- C++ -*-
/*
  //comment by wooseung

  This class has two purposes.
  The first purpose is to find the vertex between two tracks, which is used in Lambda reconstruction and the tracking process.
  The second purpose is for clustering in multitrack for accidental concidence evnet veto.
  These two cases are distinguished by using different constructors.
 */

#include "TPCVertex.hh"

#include <iostream>
#include <iterator>

#include "std_ostream.hh"

#include <TLorentzVector.h>
#include "DebugCounter.hh"
#include "FuncName.hh"
#include "Kinematics.hh"
#include "DatabasePDG.hh"
#include "UserParamMan.hh"

#define refit_wVertex 0

namespace
{
const auto& gUser = UserParamMan::GetInstance();
const Double_t ppi_distcut = 10.; //Closest distance for p, pi at the vertex point
}

//_____________________________________________________________________________
TPCVertex::TPCVertex(Int_t id1, Int_t id2)
  : m_is_calculated(false), m_is_accidental(),
    m_vertex(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN()),
    m_angle(TMath::QuietNaN()), m_distance(1.e+10),
    m_track_id(), m_track_pid(), m_track_charge(),
    m_track_pos(), m_track_mom(), m_track_theta(),
    m_is_lambda(false),
    m_lambda_vertex(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN()),
    m_lambda_distance(TMath::QuietNaN()), m_lambda_angle(TMath::QuietNaN()),
    m_lambda_mass(), m_lambda_mom(),
    m_proton_id(), m_proton_mom(),
    m_pion_id(), m_pion_mom()

{
  debug::ObjectCounter::increase(ClassName());
  m_track_id.push_back(id1);
  m_track_id.push_back(id2);

}

//_____________________________________________________________________________
TPCVertex::TPCVertex(TVector3 vertex, std::vector<Int_t> trackid)
  : m_is_calculated(false), m_is_accidental(),
    m_vertex(vertex.x(), vertex.y(), vertex.z()),
    m_angle(TMath::QuietNaN()), m_distance(1.e+10),
    m_track_id(), m_track_pid(), m_track_charge(),
    m_track_pos(), m_track_mom(), m_track_theta(),
    m_is_lambda(false),
    m_lambda_vertex(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN()),
    m_lambda_distance(TMath::QuietNaN()), m_lambda_angle(TMath::QuietNaN()),
    m_lambda_mass(), m_lambda_mom(),
    m_proton_id(), m_proton_mom(),
    m_pion_id(), m_pion_mom()

{
  debug::ObjectCounter::increase(ClassName());
  for(Int_t i=0;i<trackid.size();i++) m_track_id.push_back(trackid[i]);
}

//_____________________________________________________________________________
TPCVertex::~TPCVertex()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
TPCVertex::Calculate(TPCLocalTrackHelix* track1, TPCLocalTrackHelix* track2)
{

  // For vertex finding (Closest point betweepn helix tracks)
  // Scanning range (-VertexScanRange, Tracklength + VertexScanRange) for each helix track
  static const Double_t VertexScanRange = gUser.GetParameter("VertexScanRange"); //mm

  Double_t helix_par1[5];
  track1 -> GetParam(helix_par1);
  Double_t scantheta1 = VertexScanRange/helix_par1[3]; //mm -> rad.
  Double_t range_theta1[2] = {track1 -> GetMint() - scantheta1,
			      track1 -> GetMaxt() + scantheta1};
  Double_t helix_par2[5];
  track2 -> GetParam(helix_par2);
  Double_t scantheta2 = VertexScanRange/helix_par2[3]; //mm -> rad.
  Double_t range_theta2[2] = {track2 -> GetMint() - scantheta2,
			      track2 -> GetMaxt() + scantheta2};

  Double_t theta1, theta2, dist;
  TVector3 vertex = Kinematics::VertexPointHelix(helix_par1, helix_par2,
						 range_theta1[0], range_theta1[1],
						 range_theta2[0], range_theta2[1],
						 theta1, theta2, dist);
  if(!TMath::IsNaN(dist) &&
     TMath::Abs(vertex.x()) < 250. &&
     TMath::Abs(vertex.z()) < 250. &&
     TMath::Abs(vertex.y()) < 250.) m_is_calculated = true;

  if(m_is_calculated){
    m_vertex = vertex;
    m_track_charge.push_back(track1 -> GetCharge());
    m_track_charge.push_back(track2 -> GetCharge());
    m_track_pid.push_back(track1 -> GetPid());
    m_track_pid.push_back(track2 -> GetPid());
    m_track_pos.push_back(track1 -> GetPosition(helix_par1, theta1));
    m_track_pos.push_back(track2 -> GetPosition(helix_par2, theta2));
    TVector3 mom1 = track1 -> CalcHelixMom(helix_par1, theta1);
    TVector3 mom2 = track2 -> CalcHelixMom(helix_par2, theta2);
    m_track_mom.push_back(mom1);
    m_track_mom.push_back(mom2);
    m_track_theta.push_back(theta1);
    m_track_theta.push_back(theta2);

    m_angle = m_track_mom[0].Angle(m_track_mom[1]);
    m_distance = dist;

    if(track1 -> GetIsK18()==0 && track2 -> GetIsK18()==0 &&
       track1 -> GetIsBeam()==0 && track2 -> GetIsBeam()==0 &&
       track1 -> GetIsKurama()==0 && track2 -> GetIsKurama()==0 &&
       track1 -> GetIsAccidental()==0 && track2 -> GetIsAccidental()==0)
      ReconstructLambda(vertex, mom1, mom2, dist); //Check wheter it is L decay vertex or not
  }

#if refit_wVertex
  if(m_is_lambda) ReconstructLambdaWithVertex(track1, track2); //refit with decay vertex
#endif

}

//_____________________________________________________________________________
Bool_t
TPCVertex::ReconstructLambda(TVector3 vertex, TVector3 mom1, TVector3 mom2, Double_t ppi_distance)
{

  Double_t lambda_masscut = 0.1; //ref
  static const auto PionMass    = pdg::PionMass();
  static const auto ProtonMass  = pdg::ProtonMass();
  static const auto LambdaMass  = pdg::LambdaMass();

  if(!IsCalculated()) return false;
  if(ppi_distance > ppi_distcut) return false;

  m_lambda_vertex = vertex;
  m_lambda_distance = ppi_distance;
  m_lambda_angle = mom1.Angle(mom2);
  m_lambda_mass.clear();
  m_lambda_mom.clear();
  m_proton_mom.clear();
  m_pion_mom.clear();
  m_proton_id.clear();
  m_pion_id.clear();
  for(Int_t i=0;i<2;i++){ //p, pi- or pi- p combination
    Int_t p_id = m_track_id[i];
    Int_t pi_id = m_track_id[1-i];
    Int_t p_pid = m_track_pid[i];
    Int_t pi_pid = m_track_pid[1-i];
    Int_t p_charge = m_track_charge[i];
    Int_t pi_charge = m_track_charge[1-i];
    TVector3 p_mom = mom1;
    TVector3 pi_mom = mom2;
    if(i!=0){
      p_mom = mom2;
      pi_mom = mom1;
    }
    if((p_pid&4)!=4 || p_charge!=1) continue;
    if((pi_pid&1)!=1 || pi_charge!=-1) continue;

    TLorentzVector Lpi(pi_mom, TMath::Hypot(pi_mom.Mag(), PionMass));
    TLorentzVector Lp(p_mom, TMath::Hypot(p_mom.Mag(), ProtonMass));
    TLorentzVector Llambda = Lp + Lpi;
    TVector3 lambda_mom = pi_mom + p_mom;
    /*
    std::cout<<"lambda candidate "<<p_id<<" "<<pi_id;
    std::cout<<" lmass "<<TMath::Abs(Llambda.M())<<" vertex "<<vertex<<std::endl;
    std::cout<<" mom "<<p_mom.Mag()<<" "<<pi_mom.Mag();
    */
    if(TMath::Abs(Llambda.M() - LambdaMass) > lambda_masscut) continue;
    //std::cout<<" l reconstructed "<<p_id<<" "<<pi_id<<std::endl;

    m_is_lambda = true;
    m_lambda_mass.push_back(Llambda.M());
    m_lambda_mom.push_back(lambda_mom);
    m_proton_mom.push_back(p_mom);
    m_pion_mom.push_back(pi_mom);
    m_proton_id.push_back(p_id);
    m_pion_id.push_back(pi_id);
  }

  return m_is_lambda;
}

//_____________________________________________________________________________
void
TPCVertex::ReconstructLambdaWithVertex(TPCLocalTrackHelix* track1, TPCLocalTrackHelix* track2)
{

  if(!IsCalculated()) return;
  if(!m_is_lambda) return;

  static const Double_t VertexScanRange = gUser.GetParameter("VertexScanRange"); //mm
  TVector3 vertex_res(0.6, 0.6, 1.0);
  TPCLocalTrackHelix *Track1 = new TPCLocalTrackHelix(track1);
  if(!Track1 -> DoFitTrackwVertex(m_lambda_vertex, vertex_res)){
    delete Track1;
    return;
  }

  TPCLocalTrackHelix *Track2 = new TPCLocalTrackHelix(track2);
  if(!Track2 -> DoFitTrackwVertex(m_lambda_vertex, vertex_res)){
    delete Track2;
    return;
  }
  Double_t helix_par1[5];
  Track1 -> GetParam(helix_par1);
  Double_t scantheta1 = VertexScanRange/helix_par1[3]; //mm -> rad.
  Double_t range_theta1[2] = {Track1 -> GetMint() - scantheta1,
			      Track1 -> GetMaxt() + scantheta1};
  Double_t helix_par2[5];
  track2 -> GetParam(helix_par2);
  Double_t scantheta2 = VertexScanRange/helix_par2[3]; //mm -> rad.
  Double_t range_theta2[2] = {Track2 -> GetMint() - scantheta2,
			      Track2 -> GetMaxt() + scantheta2};

  Double_t theta1, theta2, dist;
  TVector3 vertex = Kinematics::VertexPointHelix(helix_par1, helix_par2,
						 range_theta1[0], range_theta1[1],
						 range_theta2[0], range_theta2[1],
						 theta1, theta2, dist);
  if(dist > ppi_distcut){
    delete Track1;
    delete Track2;
    return;
  }

  TVector3 mom1 = Track1 -> CalcHelixMom(helix_par1, theta1);
  TVector3 mom2 = Track2 -> CalcHelixMom(helix_par2, theta2);
  ReconstructLambda(vertex, mom1, mom2, dist);

  delete Track1;
  delete Track2;
}

//not supported (E42 is not using this)
//_____________________________________________________________________________
void
TPCVertex::Calculate(TPCLocalTrack* track1, TPCLocalTrack* track2)
{

  /*
  // For vertex finding (Closest point betweepn helix tracks)
  // Scanning range (-VertexScanRange, Tracklength + VertexScanRange) for each helix track
  static const Double_t VertexScanRange = gUser.GetParameter("VertexScanRange"); //mm

  Double_t par1[5];
  track1 -> GetParam(par1);
  Double_t par2[5];
  track2 -> GetParam(par2);

  Double_t dist;
  TVector3 vertex = Kinematics::VertexPoint();

  if(!TMath::IsNaN(dist) &&
     TMath::Abs(vertex.x()) < 250. &&
     TMath::Abs(vertex.z()) < 250. &&
     TMath::Abs(vertex.y()) < 250.) m_is_calculated = true;

  if(m_is_calculated){
    m_vertex = vertex;
    m_track_pos.push_back(track1 -> GetPosition(par1, theta1));
    m_track_pos.push_back(track2 -> GetPosition(par2, theta2));
    m_angle = ;
    m_distance = dist;
  }
  */
}

//_____________________________________________________________________________
void
TPCVertex::Print(const TString& arg, Bool_t print_all) const
{
  hddaq::cerr << arg << std::endl
	      << "Vertex point = " << m_vertex << std::endl
	      << "Closest distance = " << m_distance << " mm, "
	      << "Opening angle = " << m_angle << std::endl;
  if(print_all){
    hddaq::cerr << " Track1 #id : "<< m_track_id[0] << std::endl
		<< " closest point " << m_track_pos[0]
		<< " mm, theta "<<m_track_theta[0] << std::endl
		<< " mom "<<m_track_mom[0] << " GeV/c"<<std::endl;
    hddaq::cerr << " Track2 #id : "<< m_track_id[1] << std::endl
		<< " closest point " << m_track_pos[1]
		<< " mm, theta "<<m_track_theta[1] << std::endl
		<< " mom "<<m_track_mom[1] << " GeV/c"<<std::endl;
  }

}

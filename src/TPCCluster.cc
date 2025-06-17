 // -*- C++ -*-

#include <iostream>
#include <iterator>

#include <escape_sequence.hh>
#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "Exception.hh"
#include "FuncName.hh"
#include "PrintHelper.hh"
#include "TPCCluster.hh"
#include "TPCPadHelper.hh"
#include "TPCPositionCorrector.hh"
#include "ThreeVector.hh"

namespace
{
const auto& gTPCPos = TPCPositionCorrector::GetInstance();
}

//_____________________________________________________________________________
TPCCluster::TPCCluster(Int_t layer, const TPCHitContainer& HitCont)
  : m_is_good(false),
    m_is_onframe(false),
    m_layer(layer),
    m_cluster_de(),
    m_cluster_position(),
    m_hit_array(HitCont), // shallow copy
    m_mean_row(),
    m_mean_theta(),
    m_center_hitid(),
    m_mean_hit(new TPCHit(layer, TMath::QuietNaN()))
{
  auto itr = m_hit_array.begin();
  while(itr != m_hit_array.end()){
    if(!*itr || !(*itr)->IsGood()){
      itr = m_hit_array.erase(itr);
    }else{
      itr++;
    }
  }
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCCluster::~TPCCluster()
{
  ClearTPCHits();
  delete m_mean_hit;
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
TPCCluster::ClearTPCHits()
{
  m_hit_array.clear();
}

//_____________________________________________________________________________
void
TPCCluster::AddTPCHit(TPCHit* hit)
{
  if(hit) m_hit_array.push_back(hit);
  m_is_good = false;
  m_is_onframe = false;
}

//_____________________________________________________________________________
void
TPCCluster::CheckClusterOnTheFrame()
{
  m_is_onframe = false;
  if(m_layer<8) return; // no frame

  Int_t low_row = 10000; Int_t high_row = -1; //get edge pad ids of the cluster
  for(Int_t i=0; i<m_hit_array.size(); ++i){
    if(!m_hit_array[i]) continue;
    Int_t row = m_hit_array[i] -> GetRow();
    if(row<low_row) low_row = row;
    if(row>high_row) high_row = row;
  }

  Bool_t status = false;
  for(Int_t i=0; i<5; ++i){
    if(TMath::Abs(tpc::FrameHighEdge[m_layer][i] - low_row) <= tpc::MaxRowDifTPC) status = true;
    if(TMath::Abs(tpc::FrameLowEdge[m_layer][i] - high_row) <= tpc::MaxRowDifTPC) status = true;
  }

  m_is_onframe = status;
}

//_____________________________________________________________________________
Bool_t
TPCCluster::Calculate()
{
  static const TVector2 target_center(0., tpc::ZTarget); // (X, Z)
  //const Double_t R = tpc::GetRadius(m_layer);
  //int max_row = tpc::padParameter[m_layer][1];
  m_cluster_de = 0.;
  m_cluster_position.SetXYZ(0., 0., 0.);

  //w/o position correction
  TVector2 xz_vectorHS0(0., 0.);

  //w/ position correction
  m_mean_row = 0.;
  m_mean_theta = 0.;
  Double_t mean_y = 0.;
  TVector2 xz_vectorHS(0., 0.);
  for(const auto& hit: m_hit_array){
    const auto& pos = hit->GetPosition();
    TVector2 xz_vector(pos.X(), pos.Z());
    xz_vector -= target_center;

#if 0 //Legacy. After position correction, the size of xz_vector can be differ from the layer R.
    if(TMath::Abs(xz_vector.Mod() - R) > 1e-10){
      hit->Print();
      throw Exception(FUNC_NAME + Form(" found invalid radius %lf/%lf",
                                       R, xz_vector.Mod()));
    }
#endif
    const Double_t de = hit->GetCDe();
    mean_y += pos.Y() * de;
    m_cluster_de += de;
    xz_vectorHS += de*xz_vector;

    Int_t row = hit -> GetRow();
    Int_t padid = tpc::GetPadId(m_layer, row);
    auto pos0 = tpc::getPosition(padid);
    TVector2 xz_vector0(pos0.X(), pos0.Z());
    xz_vector0 -= target_center;
    xz_vectorHS0 += de*xz_vector0;
  }

  mean_y *= 1./m_cluster_de;
  xz_vectorHS *= 1./m_cluster_de;
  m_mean_theta = xz_vectorHS.Phi();
  TVector2 xz_vector = xz_vectorHS + target_center;
  m_cluster_position.SetXYZ(xz_vector.X(), mean_y, xz_vector.Y());
  m_mean_row = tpc::getMrow(m_layer, m_mean_theta*TMath::RadToDeg());
  m_mean_hit->AddHit(0., 0.);
  m_mean_hit->SetMRow(m_mean_row);
  m_mean_hit->SetPadLength(tpc::padParameter[m_layer][5]);
  m_mean_hit->SetPadTheta(tpc::getTheta(m_layer, m_mean_row)*TMath::DegToRad());
  m_mean_hit->SetDe(m_cluster_de);
  m_mean_hit->SetPosition(m_cluster_position);
  m_mean_hit->SetParentCluster(this);

  //center hit determination
  Int_t max_row = tpc::padParameter[m_layer][3];
  Double_t mean_phi0 = xz_vectorHS0.Phi();
  Double_t mean_row0 = tpc::getMrow(m_layer, mean_phi0*TMath::RadToDeg());

  Double_t rowdiff = 10; Int_t id = -1;
  for(Int_t i=0; i<m_hit_array.size(); ++i){
    if(!m_hit_array[i]) continue;
    Double_t row = (double) m_hit_array[i] -> GetRow();
    Double_t dif_row = std::min(abs(mean_row0-row), abs(mean_row0-row+max_row));
    dif_row = std::min(dif_row, abs(mean_row0-row-max_row));
    if(dif_row<rowdiff){
      rowdiff = dif_row;
      id = i;
    }
  }
  m_center_hitid = id;
  CheckClusterOnTheFrame(); //check whether the cluster on the frame or not
  m_is_good = true;
  return true;
}

//_____________________________________________________________________________
void
TPCCluster::Print(Option_t* opt) const
{
  PrintHelper helper(1, std::ios::fixed);
  const Double_t R = tpc::GetRadius(m_layer);
  hddaq::cout << FUNC_NAME << " " << std::endl
              << "is good = "  << m_is_good  << std::endl
              << "de = "  << m_cluster_de << std::endl
              << "position = " << m_cluster_position << std::endl
              << "Radius = " << R << std::endl
              << "mean row = " << m_mean_row << std::endl
              << "mean phi = " << m_mean_theta*TMath::RadToDeg() << " (in XZ plane)" << std::endl;
  hddaq::cout << "layer" << std::setw(2) << m_layer <<" size="
              << std::setw(3) << m_hit_array.size() << "  ";
  for(const auto& hit: m_hit_array){
    hddaq::cout << "(de=" << hddaq::unpacker::esc::k_purple
                << hit->GetCDe() << hddaq::unpacker::esc::k_default_color
                << ", y=" << hddaq::unpacker::esc::k_cyan
                << hit->GetDriftLength() << hddaq::unpacker::esc::k_default_color
                << ")" << " ";
    const auto& pos = hit->GetPosition();
    TVector2 xz_vector(pos.X(), pos.Z());
    TVector2 target_position(0., tpc::ZTarget);
    auto phi = (xz_vector - target_position).Phi()*TMath::RadToDeg();
    auto residual = tpc::ArcLength(m_layer, hit->GetRow(), m_mean_row);
    hddaq::cout << " pos=" << pos << ", phi=" << phi << ", res=" << residual << std::endl;
  }
  hddaq::cout << std::endl;
}

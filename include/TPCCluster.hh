// -*- C++ -*-

#ifndef TPC_CLUSTER_HH
#define TPC_CLUSTER_HH

#include <cstddef>
#include <string>
#include <vector>

#include <TString.h>
#include <TVector3.h>

#include "TPCHit.hh"

typedef std::vector<TPCHit*> TPCHitContainer;

//_____________________________________________________________________________
class TPCCluster
{
public:
  static const TString& ClassName();
  TPCCluster(Int_t layer, const TPCHitContainer& HitCont);
  ~TPCCluster();

private:
  Bool_t          m_is_good;
  Bool_t          m_is_onframe;
  Int_t           m_layer;
  Double_t        m_cluster_de;
  TVector3        m_cluster_position;
  TPCHitContainer m_hit_array;
  Double_t        m_mean_row;
  Double_t        m_mean_theta; // in XZ plane
  Int_t           m_center_hitid;
  Int_t           m_clsize_g4;
  TPCHit*         m_mean_hit; // representative hit for tracking

public:
  void     AddTPCHit(TPCHit* hit);
  Bool_t   IsOnTheFrame(){ return m_is_onframe; };
  void     CheckClusterOnTheFrame();
  Bool_t   Calculate();
  void     ClearTPCHits();
  Int_t    GetClusterSize() const { return m_hit_array.size(); }
  Int_t    GetClusterSizeG4() const { return m_clsize_g4; }
  Double_t GetDe() const { return m_cluster_de; }
  TPCHit*  GetHit(Int_t i) const { return m_hit_array.at(i); }
  Double_t GetPadTheta() { return m_mean_theta; }
  const TPCHitContainer& GetHitContainer() const { return m_hit_array; }
  Bool_t   IsGood() const { return m_is_good; }
  TPCHit*  GetMeanHit() const { return m_mean_hit; }
  Double_t MeanRow() const { return m_mean_row; }
  const TVector3& GetPosition() const { return m_cluster_position; }
  Double_t GetX() const { return m_cluster_position.X(); }
  Double_t GetY() const { return m_cluster_position.Y(); }
  Double_t GetZ() const { return m_cluster_position.Z(); }
  TPCHit*  GetCenterHit() const { return m_hit_array.at(m_center_hitid); }
  Int_t    GetHoughFlag() const { return m_mean_hit->GetHoughFlag(); }
  const std::vector<Double_t>& GetResolutionParams() const { return m_mean_hit -> GetResolutionParams(); }
  void     Print(Option_t* opt="") const;
  void     SetClusterSizeG4(int clsize) { m_clsize_g4 = clsize; }
};

//_____________________________________________________________________________
inline const TString&
TPCCluster::ClassName()
{
  static TString s_name("TPCCluster");
  return s_name;
}

#endif

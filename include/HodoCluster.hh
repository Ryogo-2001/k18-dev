// -*- C++ -*-

#ifndef HODO_CLUSTER_HH
#define HODO_CLUSTER_HH

#include <TString.h>

class HodoHit;
class HodoAnalyzer;

//_____________________________________________________________________________
class HodoCluster
{
public:
  static const TString& ClassName();
  HodoCluster(HodoHit* hitA, HodoHit* hitB=nullptr, HodoHit* hitC=nullptr);
  virtual ~HodoCluster();

private:
  HodoCluster(const HodoCluster &);
  HodoCluster& operator =(const HodoCluster &);

private:
  HodoHit* m_hitA;
  HodoHit* m_hitB;
  HodoHit* m_hitC;
  Int_t      m_indexA;
  Int_t      m_indexB;
  Int_t      m_indexC;
  Int_t       m_cluster_size;
  Double_t    m_mean_time;
  Double_t    m_cmean_time;
  Double_t    m_de;
  Double_t    m_mean_seg;
  Double_t    m_time_diff;
  Double_t    m_1st_seg;
  Double_t    m_1st_time;
  Bool_t      m_good_for_analysis;

public:
  Double_t  C1stSeg() const { return m_1st_seg; }
  Double_t  C1stTime() const { return m_1st_time; }
  void      Calculate();
  Int_t     ClusterSize() const { return m_cluster_size; }
  Double_t  CMeanTime() const { return m_cmean_time; }
  Double_t  DeltaE() const { return m_de; }
  HodoHit*  GetHit(Int_t i) const;
  Bool_t    GoodForAnalysis() const { return m_good_for_analysis; }
  Bool_t    GoodForAnalysis(Bool_t status);
  Double_t  MeanTime() const { return m_mean_time; }
  Double_t  MeanSeg() const { return m_mean_seg; }
  Bool_t    ReCalc(Bool_t applyRecusively=false);
  void      SetIndex(Int_t iA, Int_t iB=0, Int_t iC=0);
  Double_t  TimeDif() const { return m_time_diff; }
};

//_____________________________________________________________________________
inline const TString&
HodoCluster::ClassName()
{
  static TString s_name("HodoCluster");
  return s_name;
}

#endif

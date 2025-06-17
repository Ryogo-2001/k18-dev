// -*- C++ -*-

#ifndef BH2_CLUSTER_HH
#define BH2_CLUSTER_HH

class BH2Hit;
class HodoAnalyzer;

#include <TString.h>

//_____________________________________________________________________________
class BH2Cluster
{
public:
  static const TString& ClassName();
  BH2Cluster(BH2Hit *hitA, BH2Hit *hitB=nullptr, BH2Hit *hitC=nullptr);
  ~BH2Cluster();

private:
  BH2Cluster(const BH2Cluster&);
  BH2Cluster& operator =(const BH2Cluster&);

private:
  BH2Hit*  m_hitA;
  BH2Hit*  m_hitB;
  BH2Hit*  m_hitC;
  Int_t    m_indexA;
  Int_t    m_indexB;
  Int_t    m_indexC;
  Int_t    m_cluster_size;
  Double_t m_mean_time;
  Double_t m_cmean_time;
  Double_t m_de;
  Double_t m_mean_seg;
  Double_t m_time_diff;
  Double_t m_time0;
  Double_t m_ctime0;
  Bool_t   m_good_for_analysis;

public:
  void     Calculate();
  Int_t    ClusterSize() const { return m_cluster_size; }
  Double_t CMeanTime() const { return m_cmean_time; }
  Double_t CTime0() const { return m_ctime0; }
  Double_t DeltaE() const { return m_de; }
  BH2Hit*  GetHit(Int_t i) const;
  Bool_t   GoodForAnalysis() const { return m_good_for_analysis; }
  Bool_t   GoodForAnalysis(Bool_t status);
  Double_t MeanSeg() const { return m_mean_seg; }
  Double_t MeanTime() const { return m_mean_time; }
  Double_t Time0() const { return m_time0; }
  Double_t TimeDif() const { return m_time_diff; }
  Bool_t   ReCalc(Bool_t applyRecusively=false);
  void     SetIndex(Int_t iA, Int_t iB=0, Int_t iC=0);
};

//_____________________________________________________________________________
inline const TString&
BH2Cluster::ClassName()
{
  static TString s_name("BH2Cluster");
  return s_name;
}

#endif

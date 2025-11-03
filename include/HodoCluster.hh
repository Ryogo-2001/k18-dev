// -*- C++ -*-

#ifndef HODO_CLUSTER_HH
#define HODO_CLUSTER_HH

#include <TString.h>

#include "HodoHit.hh"

using index_t = std::vector<Int_t>;

//_____________________________________________________________________________
class HodoCluster
{
public:
  static const TString& ClassName();
  HodoCluster(const HodoHC& cont,
              const index_t& index);
  virtual ~HodoCluster();

private:
  HodoCluster(const HodoCluster&);
  HodoCluster& operator =(const HodoCluster&);

protected:
  Bool_t   m_is_good;
  HodoHC   m_hit_container;
  index_t  m_index;
  Int_t    m_cluster_size;
  Double_t m_mean_time;
  Double_t m_ctime;
  Double_t m_time_diff;
  Double_t m_de;
  Double_t m_tot;
  Double_t m_mean_position;
  Double_t m_segment;
  Double_t m_1st_seg;
  Double_t m_1st_time;
  Double_t m_time0;
  Double_t m_ctime0;

public:
  const TString& DetectorName() const
    { return m_hit_container.at(0)->DetectorName(); }
  const TString& PlaneName() const
    { return m_hit_container.at(0)->PlaneName(); }
  Int_t     DetectorId() const { return m_hit_container.at(0)->DetectorId(); }
  Int_t     PlaneId() const { return m_hit_container.at(0)->PlaneId(); }
  Bool_t    IsGood() const { return m_is_good; }
  Int_t     ClusterSize() const { return m_cluster_size; }
  Double_t  MeanPosition() const { return m_mean_position; }
  Double_t  MeanTime() const { return m_mean_time; }
  Double_t  CMeanTime() const { return m_ctime; }
  Double_t  CTime() const { return CMeanTime(); }
  Double_t  TimeDiff() const { return m_time_diff; }
  Double_t  Time0() const { return m_time0; }
  Double_t  CTime0() const { return m_ctime0; }
  Double_t  DeltaE() const { return m_de; }
  Double_t  TOT() const { return m_tot; }
  Double_t  MeanSeg() const { return m_segment; }
  Double_t  FirstSeg() const { return m_1st_seg; }
  Double_t  FirstTime() const { return m_1st_time; }
  HodoHit*  GetHit(Int_t i) const;
  void      Print(Option_t* opt="") const;
  Bool_t    ReCalc(Bool_t applyRecusively=false);

protected:
  void      Calculate();
};

//_____________________________________________________________________________
inline const TString&
HodoCluster::ClassName()
{
  static TString s_name("HodoCluster");
  return s_name;
}

#endif

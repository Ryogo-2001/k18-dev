// -*- C++ -*-

#ifndef FIBER_CLUSTER_HH
#define FIBER_CLUSTER_HH

#include <vector>

#include <TString.h>

class FLHit;

//_____________________________________________________________________________
class FiberCluster
{
public:
  static const TString& ClassName();
  FiberCluster();
  virtual ~FiberCluster();

private:
  FiberCluster(const FiberCluster& object);
  FiberCluster& operator =(const FiberCluster& object);

protected:
  typedef std::vector<FLHit*> HitContainer;
  enum FlagsFiber { Initialized, gfastatus, sizeFlagsFiber };
  HitContainer m_hit_container;
  Int_t          m_cluster_size;
  Int_t          m_cluster_id;
  Int_t          m_max_cluster_id;
  Double_t       m_mean_time;
  Double_t       m_max_time;
  Double_t       m_real_mean_time;
  // real mean (not a closest value of CTime)
  Double_t       m_max_width;
  Double_t       m_min_width;
  Double_t       m_mean_seg;
  Double_t       m_mean_pos;
  Double_t       m_sum_adc_lg;
  Double_t       m_sum_mip_lg;
  Double_t       m_sum_de_lg;
  Double_t       m_max_adc_hg;
  Double_t       m_max_adc_lg;
  Double_t       m_max_mip_lg;
  Double_t       m_max_de_lg;
  Double_t       m_max_seg;
  Bool_t         m_flag[sizeFlagsFiber];

public:
  Bool_t   Calculate();
  void     push_back(FLHit* hit) { m_hit_container.push_back(hit); }
  Int_t    VectorSize() const { return m_hit_container.size(); }
  Int_t    ClusterId() const { return m_cluster_id; }
  Int_t    ClusterSize() const { return m_cluster_size; }
  Int_t    GetMaxClusterId() const { return m_max_cluster_id; }
  Double_t CMeanTime() const { return m_mean_time; }
  Double_t CMaxTime() const { return m_max_time; }
  Double_t RCMeanTime() const { return m_real_mean_time; }
  Double_t Width() const { return m_max_width; }
  Double_t minWidth() const { return m_min_width; }
  Double_t Tot() const { return Width(); }
  Double_t MeanPosition() const { return m_mean_pos; }
  Double_t SumAdcLG() const { return m_sum_adc_lg; }
  Double_t SumMipLG() const { return m_sum_mip_lg; }
  Double_t SumDeLG() const { return m_sum_de_lg; }
  Double_t MeanSeg() const { return m_mean_seg; }
  Double_t MaxAdcHi() const { return m_max_adc_hg; }
  Double_t MaxAdcLG() const { return m_max_adc_lg; }
  Double_t MaxMipLG() const { return m_max_mip_lg; }
  Double_t MaxDeLG() const { return m_max_de_lg; }
  Double_t MaxSeg() const { return m_max_seg; }
  FLHit*   GetHit(Int_t i) const;
  Bool_t   GoodForAnalysis() const { return m_flag[gfastatus]; }
  Bool_t   GoodForAnalysis(Bool_t status);
  Bool_t   ReCalc(Bool_t applyRecusively=false);

private:
  void Debug();

};

//_____________________________________________________________________________
inline const TString&
FiberCluster::ClassName()
{
  static TString s_name("FiberCluster");
  return s_name;
}

#endif

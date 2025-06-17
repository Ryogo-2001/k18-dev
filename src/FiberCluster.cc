// -*- C++ -*-

#include "FiberCluster.hh"

#include <cmath>
#include <string>
#include <limits>

#include <TMath.h>

#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "FiberHit.hh"
#include "FLHit.hh"
#include "FuncName.hh"

namespace
{
const Bool_t reject_nan = false;
}

//_____________________________________________________________________________
FiberCluster::FiberCluster()
  : m_cluster_size(0),
    m_max_cluster_id(-1),
    m_max_time(TMath::QuietNaN()),
    m_max_adc_hg(0.),
    m_max_adc_lg(0.),
    m_max_mip_lg(0.),
    m_max_de_lg(0.)
{
  for(Int_t i=0; i<sizeFlagsFiber; ++i){
    m_flag[i] = false;
  }
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
FiberCluster::~FiberCluster()
{
  m_hit_container.clear();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
FLHit*
FiberCluster::GetHit(Int_t i) const
{
  return m_hit_container.at(i);
}

//_____________________________________________________________________________
Bool_t
FiberCluster::GoodForAnalysis(Bool_t status)
{
  Bool_t ret = m_flag[gfastatus];
  m_flag[gfastatus] = status;
  return ret;
}

//_____________________________________________________________________________
Bool_t
FiberCluster::Calculate()
{
  if(m_flag[Initialized]){
    hddaq::cerr << FUNC_NAME << " "
		<< "Already initialied" << std::endl;
    return false;
  }

  if(0 ==(m_cluster_size = m_hit_container.size())){
    hddaq::cerr << FUNC_NAME << " "
		<< "No FiberHit in local container" << std::endl;
    return false;
  }

  Double_t mean_seg       = 0.;
  Double_t mean_pos       = 0.;
  Double_t sum_adc_lg     = 0.;
  Double_t sum_mip_lg     = 0.;
  Double_t sum_de_lg      = 0.;
  Double_t mean_time      = TMath::QuietNaN();
  Double_t max_time       = TMath::QuietNaN();
  Double_t max_width      = 0.;
  Double_t min_width      = 0.;
  Double_t real_mean_time = 0.;
  Int_t    cluster_id     = 0;
  Double_t max_seg        = 0.;
  Double_t max_adc_hg     = 0.;
  Double_t max_adc_lg     = 0.;
  Double_t max_mip_lg     = 0.;
  Double_t max_de_lg      = 0.;

  for(Int_t i=0; i<m_cluster_size; ++i){
    if(reject_nan && std::isnan(m_hit_container.at(i)->GetWidth())){
      --m_cluster_size;
      continue;
    }

    mean_seg       += m_hit_container.at(i)->SegmentId();
    mean_pos       += m_hit_container.at(i)->GetPosition();
    sum_adc_lg     += m_hit_container.at(i)->GetAdcLG();
    sum_mip_lg     += m_hit_container.at(i)->GetMipLG();
    sum_de_lg      += m_hit_container.at(i)->GetDeLG();
    real_mean_time += m_hit_container.at(i)->GetCTime();
    Double_t time = m_hit_container.at(i)->GetCTime();

    if(std::isnan(mean_time)){
      mean_time = time;
    }else if(std::abs(time) < std::abs(mean_time)){
      mean_time = time;
    }

    if(max_width < m_hit_container.at(i)->GetWidth()){
      max_width = m_hit_container.at(i)->GetWidth();
    }
    if(min_width > m_hit_container.at(i)->GetWidth()){
      min_width = m_hit_container.at(i)->GetWidth();
    }

    if(max_adc_lg <= m_hit_container.at(i)->GetAdcLG()){
      max_adc_hg = m_hit_container.at(i)->GetAdcHG();
      max_adc_lg = m_hit_container.at(i)->GetAdcLG();
      max_seg = m_hit_container.at(i)->SegmentId();
      m_max_cluster_id = i;
      max_time = m_hit_container.at(i)->GetCTime();
    }

    if(max_mip_lg <= m_hit_container.at(i)->GetMipLG())
      max_mip_lg = m_hit_container.at(i)->GetMipLG();

    if(max_de_lg <= m_hit_container.at(i)->GetDeLG())
      max_de_lg = m_hit_container.at(i)->GetDeLG();

    cluster_id  += m_hit_container.at(i)->PairId();
  }

  mean_seg       /= Double_t(m_cluster_size);
  mean_pos       /= Double_t(m_cluster_size);
  real_mean_time /= Double_t(m_cluster_size);

  m_mean_time      = mean_time;
  m_max_width      = max_width;
  m_min_width      = min_width;
  m_mean_seg       = mean_seg;
  m_mean_pos       = mean_pos;
  m_sum_adc_lg     = sum_adc_lg;
  m_sum_mip_lg     = sum_mip_lg;
  m_sum_de_lg      = sum_de_lg;

  m_max_adc_hg      = max_adc_hg;
  m_max_adc_lg     = max_adc_lg;
  m_max_mip_lg     = max_mip_lg;
  m_max_de_lg      = max_de_lg;
  m_max_seg        = max_seg;
  m_max_time       = max_time;

  m_real_mean_time = real_mean_time;
  m_cluster_id     = (m_cluster_size == 1) ? 2*cluster_id : cluster_id;

  if(m_cluster_size == 0){
    //    std::cout << "FiberCluster is destroied" << std::endl;
    return false;
  }

  m_flag[Initialized] = true;

#ifdef DEBUG
  Debug();
#endif

  return true;
}

//_____________________________________________________________________________
void
FiberCluster::Debug()
{
  hddaq::cout << "Used hit" << std::endl;
  for(Int_t i = 0; i<m_cluster_size; ++i){
    m_hit_container.at(i)->Dump();
  }
}

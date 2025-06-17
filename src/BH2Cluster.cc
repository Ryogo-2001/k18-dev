// -*- C++ -*-

#include "BH2Cluster.hh"

#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "BH2Hit.hh"
#include "DebugCounter.hh"
#include "HodoAnalyzer.hh"

//_____________________________________________________________________________
BH2Cluster::BH2Cluster(BH2Hit *hitA, BH2Hit *hitB, BH2Hit *hitC)
  : m_hitA(hitA),
    m_hitB(hitB),
    m_hitC(hitC),
    m_indexA(0),
    m_indexB(0),
    m_indexC(0),
    m_cluster_size(0),
    m_mean_time(),
    m_cmean_time(),
    m_de(),
    m_mean_seg(),
    m_time_diff(),
    m_time0(),
    m_ctime0(),
    m_good_for_analysis(true)
{
  if(hitA) ++m_cluster_size;
  if(hitB) ++m_cluster_size;
  if(hitC) ++m_cluster_size;
  //  Calculate();
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
BH2Cluster::~BH2Cluster()
{
  debug::ObjectCounter::decrease(ClassName());
};

//_____________________________________________________________________________
void
BH2Cluster::Calculate()
{
  Double_t ms=0., mt=0., cmt=0., de=0., t0=0., ct0=0., dt=0;
  if(m_hitA){
    ms  += m_hitA->SegmentId();
    mt  += m_hitA->MeanTime(m_indexA);
    cmt += m_hitA->CMeanTime(m_indexA);
    de  += m_hitA->DeltaE();
    t0  += m_hitA->Time0(m_indexA);
    ct0 += m_hitA->CTime0(m_indexA);
    dt  += (m_hitA->GetTDown(m_indexA) - m_hitA->GetTUp(m_indexA));
  }
  if(m_hitB){
    ms  += m_hitB->SegmentId();
    mt  += m_hitB->MeanTime(m_indexB);
    cmt += m_hitB->CMeanTime(m_indexB);
    de  += m_hitB->DeltaE();
    t0  += m_hitB->Time0(m_indexB);
    ct0 += m_hitB->CTime0(m_indexB);
    dt  += (m_hitB->GetTDown(m_indexB) - m_hitB->GetTUp(m_indexB));
  }
  if(m_hitC){
    ms  += m_hitC->SegmentId();
    mt  += m_hitC->MeanTime(m_indexC);
    cmt += m_hitC->CMeanTime(m_indexC);
    de  += m_hitC->DeltaE();
    t0  += m_hitA->Time0(m_indexC);
    ct0 += m_hitA->CTime0(m_indexC);
    dt  += (m_hitC->GetTDown(m_indexC) - m_hitB->GetTUp(m_indexC));
  }

  ms /= Double_t(m_cluster_size);
  mt /= Double_t(m_cluster_size);
  cmt/= Double_t(m_cluster_size);
  t0 /= Double_t(m_cluster_size);
  ct0/= Double_t(m_cluster_size);
  dt /= Double_t(m_cluster_size);

  m_mean_seg  = ms;
  m_mean_time = mt;
  m_cmean_time= cmt;
  m_de        = de;
  m_time0     = t0;
  m_ctime0    = ct0;
  m_time_diff = dt;
}

//_____________________________________________________________________________
BH2Hit*
BH2Cluster::GetHit(Int_t i) const
{
  switch(i){
  case 0: return m_hitA;
  case 1: return m_hitB;
  case 2: return m_hitC;
  default: return nullptr;
  }
}

//_____________________________________________________________________________
Bool_t
BH2Cluster::GoodForAnalysis(Bool_t status)
{
  Bool_t pre_status = m_good_for_analysis;
  m_good_for_analysis = status;
  return pre_status;
}

//_____________________________________________________________________________
Bool_t
BH2Cluster::ReCalc(Bool_t applyRecursively)
{
  if(applyRecursively){
    if(m_hitA) m_hitA->ReCalc(applyRecursively);
    if(m_hitB) m_hitB->ReCalc(applyRecursively);
    if(m_hitC) m_hitC->ReCalc(applyRecursively);
  }
  Calculate();
  return true;
}

//_____________________________________________________________________________
void
BH2Cluster::SetIndex(Int_t iA, Int_t iB, Int_t iC)
{
  m_indexA = iA;
  m_indexB = iB;
  m_indexC = iC;
}

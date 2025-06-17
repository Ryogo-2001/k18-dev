// -*- C++ -*-

#include "HodoCluster.hh"

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>

#include "DebugCounter.hh"
#include "HodoHit.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"

//_____________________________________________________________________________
HodoCluster::HodoCluster(HodoHit* hitA, HodoHit* hitB, HodoHit* hitC)
  : m_hitA(hitA),
    m_hitB(hitB),
    m_hitC(hitC),
    m_indexA(0),
    m_indexB(0),
    m_indexC(0),
    m_cluster_size(0),
    m_good_for_analysis(true)
{
  if(hitA) ++m_cluster_size;
  if(hitB) ++m_cluster_size;
  if(hitC) ++m_cluster_size;
  //  Calculate();
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
HodoCluster::~HodoCluster()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
HodoCluster::Calculate()
{
  Double_t ms = 0.;
  Double_t mt = 0.;
  Double_t cmt= 0.;
  Double_t de = 0.;
  Double_t dt = 0.;
  Double_t s1 = 0.;
  Double_t t1 = 0.;

  if(m_hitA){
    ms += m_hitA->SegmentId();
    mt += m_hitA->MeanTime(m_indexA);
    cmt+= m_hitA->CMeanTime(m_indexA);
    de += m_hitA->DeltaE();
    dt += (m_hitA->GetTDown(m_indexA)-m_hitA->GetTUp(m_indexA));
    s1  = m_hitA->SegmentId();
    t1  = m_hitA->CMeanTime(m_indexA);
  }
  if(m_hitB){
    ms += m_hitB->SegmentId();
    mt += m_hitB->MeanTime(m_indexB);
    cmt+= m_hitB->CMeanTime(m_indexB);
    de += m_hitB->DeltaE();
    dt += (m_hitB->GetTDown(m_indexB)-m_hitB->GetTUp(m_indexB));
    if(m_hitB->CMeanTime(m_indexB)<t1){
      s1 = m_hitB->SegmentId();
      t1 = m_hitB->CMeanTime(m_indexB);
    }
  }
  if(m_hitC){
    ms += m_hitC->SegmentId();
    mt += m_hitC->MeanTime(m_indexC);
    cmt+= m_hitC->CMeanTime(m_indexC);
    de += m_hitC->DeltaE();
    dt += (m_hitC->GetTDown(m_indexC)-m_hitC->GetTUp(m_indexC));
    if(m_hitC->CMeanTime(m_indexC)<t1){
      s1 = m_hitC->SegmentId();
      t1 = m_hitC->CMeanTime(m_indexC);
    }
  }
  ms /= Double_t(m_cluster_size);
  mt /= Double_t(m_cluster_size);
  cmt/= Double_t(m_cluster_size);
  dt /= Double_t(m_cluster_size);

  m_mean_time = mt;
  m_cmean_time= cmt;
  m_de        = de;
  m_mean_seg  = ms;
  m_time_diff = dt;
  m_1st_seg   = s1;
  m_1st_time  = t1;
}

//_____________________________________________________________________________
HodoHit*
HodoCluster::GetHit(Int_t i) const
{
  if(i==0) return m_hitA;
  if(i==1) return m_hitB;
  if(i==2) return m_hitC;
  return nullptr;
}

//_____________________________________________________________________________
Bool_t
HodoCluster::GoodForAnalysis(Bool_t status)
{
  Bool_t pre_status = m_good_for_analysis;
  m_good_for_analysis = status;
  return pre_status;
}

//_____________________________________________________________________________
Bool_t
HodoCluster::ReCalc(Bool_t applyRecursively)
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
HodoCluster::SetIndex(Int_t iA, Int_t iB, Int_t iC)
{
  m_indexA = iA;
  m_indexB = iB;
  m_indexC = iC;
}

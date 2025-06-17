// -*- C++ -*-

#include "DCPairHitCluster.hh"

#include "DCLTrackHit.hh"
#include "DebugCounter.hh"
#include "FuncName.hh"

//_____________________________________________________________________________
DCPairHitCluster::DCPairHitCluster(DCLTrackHit* hitA, DCLTrackHit* hitB)
  : m_hitA(hitA),
    m_hitB(hitB),
    m_nhits(0)
{
  if(m_hitA) ++m_nhits;
  if(m_hitB) ++m_nhits;
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
DCPairHitCluster::~DCPairHitCluster()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
DCLTrackHit*
DCPairHitCluster::GetHit(Int_t i) const
{
  if(i==0) return m_hitA;
  if(i==1) return m_hitB;
  return 0;
}

//_____________________________________________________________________________
bool
DCPairHitCluster::IsHoneycomb() const
{
  if(m_hitA) return m_hitA->IsHoneycomb();
  if(m_hitB) return m_hitB->IsHoneycomb();
  return false;
}

//_____________________________________________________________________________
void
DCPairHitCluster::SetHoneycomb(bool flag)
{
  if(m_hitA) m_hitA->SetHoneycomb(flag);
  if(m_hitB) m_hitB->SetHoneycomb(flag);
}

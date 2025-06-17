// -*- C++ -*-

#include "FLHit.hh"

#include "DebugCounter.hh"
#include "FuncName.hh"

//_____________________________________________________________________________
// One side readout method (BFT, SFT, SCH)
FLHit::FLHit(FiberHit* ptr, Int_t index)
  : m_hit_u(ptr),
    m_hit_d(nullptr),
    m_nth_hit_u(index),
    m_nth_hit_d(-1)
{
  m_hit_u->RegisterHits(this);
  m_hit_u->SetJoined(index);
  Initialize();
  debug::ObjectCounter::increase(ClassName());
}


//_____________________________________________________________________________
// Both side readout method (FBH)
FLHit::FLHit(FiberHit* ptr1,  FiberHit* ptr2,
             Int_t index1, Int_t index2)
  : m_hit_u(ptr1),
    m_hit_d(ptr2),
    m_nth_hit_u(index1),
    m_nth_hit_d(index2)
{
  m_hit_u->RegisterHits(this);
  m_hit_u->SetJoined(index1);
  m_hit_d->SetJoined(index2);
  Initialize();
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
FLHit::~FLHit()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
FLHit::Dump() const
{
  std::cout << "plid " << m_hit_u->PairId() << " ";
  std::cout << "Pos "  << GetPosition() << " ";
  std::cout << "Time " << GetCTime() << std::endl;
}

//_____________________________________________________________________________
void
FLHit::Initialize()
{
  if(m_hit_d == nullptr){
    m_leading  = m_hit_u->GetLeading(m_nth_hit_u);
    // m_trailing = m_hit_u->GetTrailing(m_nth_hit_u);
    m_trailing = m_leading; 	// Temporary solution. Need to be changed some day.
    m_time     = m_hit_u->GetTime(m_nth_hit_u);
    m_width    = m_hit_u->GetWidth(m_nth_hit_u);
    m_ctime    = m_hit_u->GetCTime(m_nth_hit_u);
  }else{
    m_leading  =
      (m_hit_u->GetLeading(m_nth_hit_u) + m_hit_d->GetLeading(m_nth_hit_d))/2;
    m_trailing  =
      (m_hit_u->GetTrailing(m_nth_hit_u) + m_hit_d->GetTrailing(m_nth_hit_d))/2;
    m_time  =
      (m_hit_u->GetTime(m_nth_hit_u) + m_hit_d->GetTime(m_nth_hit_d))/2;
    m_ctime  =
      (m_hit_u->GetCTime(m_nth_hit_u) + m_hit_d->GetCTime(m_nth_hit_d))/2;
    m_width  =
      (m_hit_u->GetWidth(m_nth_hit_u) + m_hit_d->GetWidth(m_nth_hit_d))/2;
  }
  m_flag_fljoin = false;
}

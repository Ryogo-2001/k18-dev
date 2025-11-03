// -*- C++ -*-

#ifndef FL_HIT_HH
#define FL_HIT_HH

#include "FiberHit.hh"

#include <vector>

#include <TString.h>

//_____________________________________________________________________________
class FLHit
{
public:
  static const TString& ClassName();
  // One side readout method (BFT, SFT, SCH)
  FLHit(FiberHit* ptr, Int_t index);
  // Both side readout method (FBH)
  FLHit(FiberHit* ptr1, FiberHit* ptr2,
        Int_t index1, Int_t index2);
  ~FLHit();

private:
  FLHit();
  FLHit(const FLHit& object);
  FLHit& operator =(const FLHit& object);

private:
  FiberHit *m_hit_u;
  FiberHit *m_hit_d;
  Int_t       m_nth_hit_u;
  Int_t       m_nth_hit_d;
  Double_t    m_leading;
  Double_t    m_trailing;
  Double_t    m_time;
  Double_t    m_ctime;
  Double_t    m_width;
  Bool_t      m_flag_fljoin;

public:
  Double_t GetLeading() const { return m_leading; }
  Double_t GetTrailing() const { return m_trailing; }
  Double_t GetTime() const { return m_time; }
  Double_t GetCTime() const { return m_ctime; }
  Double_t GetWidth() const { return m_width; }
  Double_t GetPosition() const { return m_hit_u->GetPosition(); }
  Double_t GetAdcHG() const { return m_hit_u->GetAdcHG(); }
  Double_t GetAdcLG() const { return m_hit_u->GetAdcLG(); }
  Double_t GetMipHG() const { return m_hit_u->GetMipHG(); }
  Double_t GetMipLG() const { return m_hit_u->GetMipLG(); }
  Double_t GetDeHG() const { return m_hit_u->GetDeHG(); }
  Double_t GetDeLG() const { return m_hit_u->GetDeLG(); }
  Int_t    PairId() const { return m_hit_u->PairId(); }
  Double_t SegmentId() const { return m_hit_u->SegmentId(); }
  void     SetJoined() { m_flag_fljoin = true; }
  Bool_t   Joined() const { return m_flag_fljoin; }
  void     Dump() const;

  friend class FiberHit;

private:
  void Initialize();

};

//_____________________________________________________________________________
inline const TString&
FLHit::ClassName()
{
  static TString s_name("FLHit");
  return s_name;
}

#endif

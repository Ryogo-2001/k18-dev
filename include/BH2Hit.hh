// -*- C++ -*-

#ifndef BH2_HIT_HH
#define BH2_HIT_HH

#include <TString.h>

#include "DebugCounter.hh"
#include "Hodo2Hit.hh"

//_____________________________________________________________________________
class BH2Hit : public Hodo2Hit
{
public:
  static const TString& ClassName();
  BH2Hit(HodoRawHit *rhit, Double_t max_time_diff=10.);
  ~BH2Hit();

private:
  Double_t m_time_offset;

public:
  Bool_t   Calculate();
  Double_t UTime0(Int_t i=0) const;
  Double_t UCTime0(Int_t i=0) const;
  Double_t DTime0(Int_t i=0) const;
  Double_t DCTime0(Int_t i=0) const;
  Double_t Time0(Int_t i=0) const;
  Double_t CTime0(Int_t i=0) const;
  virtual Bool_t ReCalc(Bool_t applyRecursively=false) { return Calculate(); }
};

//_____________________________________________________________________________
inline const TString&
BH2Hit::ClassName()
{
  static TString s_name("BH2Hit");
  return s_name;
}

#endif

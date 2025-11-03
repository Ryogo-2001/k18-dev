// -*- C++ -*-

#ifndef FIBER_HIT_HH
#define FIBER_HIT_HH

#include <string>

#include <std_ostream.hh>

// #include "DCHit.hh"
#include "HodoHit.hh"

//_____________________________________________________________________________
class FiberHit : public HodoHit
{
public:
  static const TString& ClassName();
  explicit FiberHit(HodoRawHit* hit);
  virtual  ~FiberHit();

private:
  FiberHit();
  FiberHit(const FiberHit& rhit);
  FiberHit& operator =(const FiberHit& rhit);

protected:
  Double_t m_position;
  Double_t m_dxdw;

public:
  Bool_t   Calculate();
  Double_t Position() const { return m_position; }
  Double_t dXdW() const { return m_dxdw; }
  Double_t TimeOverThreshold(Int_t i, Int_t j=0) const
    { return m_time_trailing.at(i).at(j) - m_time_leading.at(i).at(j); }
  Double_t TOT(Int_t i, Int_t j=0) const
    { return TimeOverThreshold(i, j); }
  Double_t MeanTimeOverThreshold(Int_t j=0) const;
  Double_t MeanTOT(Int_t j=0) const { return MeanTimeOverThreshold(j); }

  // aliases
  Double_t UTOT(Int_t j) const { return TimeOverThreshold(HodoRawHit::kUp, j); }
  Double_t DTOT(Int_t j) const { return TimeOverThreshold(HodoRawHit::kDown, j); }


  // Double_t CTimeOverThreshold(Int_t i, Int_t j=0) const
  //   { return m_ctime_trailing.at(i).at(j) - m_ctime_leading.at(i).at(j); }
  // Double_t CTOT(Int_t i, Int_t j=0) const
  //   { return CTimeOverThreshold(i, j); }

  // Double_t GetAdcHG() const { return m_adc_hg; }
  // Double_t GetAdcLG() const { return m_adc_lg; }
  // Double_t GetMipHG() const { return m_mip_hg; }
  // Double_t GetMipLG() const { return m_mip_lg; }
  // Double_t GetDeHG() const { return m_dE_hg; }
  // Double_t GetDeLG() const { return m_dE_lg; }
  // void     SetPedestalCor(Double_t deltaHG, Double_t deltaLG)
  // { m_pedcor_hg = deltaHG; m_pedcor_lg = deltaLG; }
  // void     RegisterHits(FLHit* hit) { m_hit_container.push_back(hit); }

  virtual void   Print(Option_t* arg="") const;
  virtual Bool_t ReCalc(Bool_t allpyRecursively=false){ return Calculate(); }

  static Bool_t Compare(const FiberHit* left, const FiberHit* right);
};

//_____________________________________________________________________________
inline const TString&
FiberHit::ClassName()
{
  static TString s_name("FiberHit");
  return s_name;
}

//_____________________________________________________________________________
inline Bool_t
FiberHit::Compare(const FiberHit* left, const FiberHit* right)
{
  return left->Position() < right->Position();
}

#endif

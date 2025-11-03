// -*- C++ -*-

#ifndef HODO_HIT_HH
#define HODO_HIT_HH

#include <cmath>
#include <cstddef>

#include <TString.h>

#include "HodoRawHit.hh"
#include "ThreeVector.hh"

class HodoHit;
// class FiberHit;
using HodoHC = std::vector<HodoHit*>;
// using FiberHC = std::vector<FiberHit*>;

//_____________________________________________________________________________
class HodoHit
{
public:
  explicit HodoHit(const HodoRawHit *rhit,
                   Double_t max_time_diff=10.);
  virtual ~HodoHit();
  static TString ClassName();

private:
  HodoHit(const HodoHit&);
  HodoHit& operator =(const HodoHit&);

protected:
  const HodoRawHit* m_raw;
  Bool_t            m_is_calculated;
  Double_t          m_max_time_diff; // [ns]
  Int_t             m_n_ch;
  Double_t          m_time_offset;

  using data_t = std::vector<std::vector<Double_t>>; // [ch][mhit]

  data_t m_de_high;
  data_t m_de_low;
  data_t m_time_leading;
  data_t m_time_trailing;
  data_t m_ctime_leading;
  data_t m_ctime_trailing; // no use
  std::vector<Bool_t> m_is_clustered;

public:
  const TString&    GetName() const { return m_raw->DetectorName(); }
  const HodoRawHit* GetRawHit() const { return m_raw; }
  Int_t             DetectorId() const { return m_raw->DetectorId(); }
  const TString&    DetectorName() const { return m_raw->DetectorName(); }
  Int_t             PlaneId() const { return m_raw->PlaneId(); }
  const TString&    PlaneName() const { return m_raw->PlaneName(); }
  Int_t             SegmentId() const { return m_raw->SegmentId(); }
  Bool_t            Calculate();
  Bool_t            IsCalculated() const { return m_is_calculated; }

  Int_t NumOfChannel() const { return m_n_ch; }
  Int_t GetEntries(Int_t i) const
    { return m_ctime_leading.at(i).size(); }
  Int_t GetEntries() const
  { if(m_n_ch == HodoRawHit::kNChannel) return GetEntries(HodoRawHit::kExtra);
    else return GetEntries(HodoRawHit::kUp); }
  Double_t GetDeltaEHighGain(Int_t i, Int_t j=0) const
    { return m_de_high.at(i).at(j); }
  Double_t GetDeltaELowGain(Int_t i, Int_t j=0) const
    { return m_de_low.at(i).at(j); }
  Double_t GetTimeLeading(Int_t i, Int_t j=0) const
    { return m_time_leading.at(i).at(j); }
  Double_t GetTimeTrailing(Int_t i, Int_t j=0) const
    { return m_time_leading.at(i).at(j); }
  Double_t GetCTimeLeading(Int_t i, Int_t j=0) const
    { return m_ctime_leading.at(i).at(j); }
  Double_t GetCTimeTrailing(Int_t i, Int_t j=0) const
    { return m_ctime_leading.at(i).at(j); }

  Double_t DeltaEHighGain(Int_t j=0) const;
  Double_t DeltaELowGain(Int_t j=0) const;
  Double_t MeanTime(Int_t j=0) const;
  Double_t CMeanTime(Int_t j=0) const;
  Double_t TimeDiff(Int_t j=0) const;
  Double_t CTimeDiff(Int_t j=0) const;

  Double_t TimeOffset() const { return m_time_offset; }
  Double_t Time0(Int_t j=0) const { return MeanTime(j) + m_time_offset; }
  Double_t CTime0(Int_t j=0) const { return CMeanTime(j) + m_time_offset; }

  // aliases
  Double_t DeltaE(Int_t j=0) const { return DeltaEHighGain(j); }
  Double_t A(Int_t i, Int_t j=0) const { return GetDeltaEHighGain(i, j); }
  Double_t GetA(Int_t i, Int_t j=0) const { return A(i, j); }
  Double_t GetAUp(Int_t j=0) const { return A(HodoRawHit::kUp, j); }
  Double_t GetALeft(Int_t j=0) const { return A(HodoRawHit::kUp, j); }
  Double_t GetADown(Int_t j=0) const { return A(HodoRawHit::kDown, j); }
  Double_t GetARight(Int_t j=0) const { return A(HodoRawHit::kDown, j); }
  Double_t UDeltaE(Int_t j=0) const { return GetAUp(j); }
  Double_t DDeltaE(Int_t j=0) const { return GetADown(j); }

  Double_t T(Int_t i, Int_t j=0) const { return GetTimeLeading(i, j); }
  Double_t GetT(Int_t i, Int_t j=0) const { return T(i, j); }
  Double_t GetTUp(Int_t j=0) const { return T(HodoRawHit::kUp, j); }
  Double_t GetTLeft(Int_t j=0) const { return T(HodoRawHit::kUp, j); }
  Double_t GetTDown(Int_t j=0) const { return T(HodoRawHit::kDown, j); }
  Double_t GetTRight(Int_t j=0) const { return T(HodoRawHit::kDown, j); }

  Double_t CT(Int_t i, Int_t j=0) const { return GetCTimeLeading(i, j); }
  Double_t GetCT(Int_t i, Int_t j=0) const { return CT(i, j); }
  Double_t GetCTUp(Int_t j=0) const { return CT(HodoRawHit::kUp, j); }
  Double_t GetCTLeft(Int_t j=0) const { return CT(HodoRawHit::kUp, j); }
  Double_t GetCTDown(Int_t j=0) const { return CT(HodoRawHit::kDown, j); }
  Double_t GetCTRight(Int_t j=0) const { return CT(HodoRawHit::kDown, j); }

  void   JoinCluster(Int_t m) { m_is_clustered.at(m) = true; }
  Bool_t IsClustered(Int_t m) const { return m_is_clustered.at(m); }
  Bool_t IsClusteredAll();

  virtual void   Print(Option_t* arg="") const;
  virtual Bool_t ReCalc(Bool_t applyRecursively=false){ return Calculate(); }

  static  Bool_t Compare(const HodoHit* left, const HodoHit* right);
};

//_____________________________________________________________________________
inline TString
HodoHit::ClassName()
{
  static TString s_name("HodoHit");
  return s_name;
}

//_____________________________________________________________________________
inline Bool_t
HodoHit::Compare(const HodoHit* left, const HodoHit* right)
{
  return left->SegmentId() < right->SegmentId();
}

#endif

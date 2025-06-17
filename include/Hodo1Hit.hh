// -*- C++ -*-

#ifndef HODO1_HIT_HH
#define HODO1_HIT_HH

#include "HodoRawHit.hh"
#include "HodoHit.hh"

//_____________________________________________________________________________
class Hodo1Hit : public HodoHit
{
public:
  explicit Hodo1Hit(HodoRawHit *rhit, Int_t index=0);
  virtual ~Hodo1Hit();
  static TString ClassName();

private:
  Hodo1Hit(const Hodo1Hit&);
  Hodo1Hit& operator =(const Hodo1Hit&);

protected:
  HodoRawHit*           m_raw;
  Bool_t                m_is_calculated;
  Int_t                 m_multi_hit_l;
  Int_t                 m_multi_hit_t;
  std::vector<Double_t> m_a;
  std::vector<Double_t> m_t;
  std::vector<Double_t> m_ct;
  std::vector<Bool_t>   m_flag_join;
  Int_t                 m_index;

public:
  HodoRawHit* GetRawHit() { return m_raw; }
  Bool_t      Calculate(Bool_t tdc_flag = true);
  Bool_t      IsCalculated() const { return m_is_calculated; }
  Int_t       GetNumOfHit(Int_t sel=0) const;
  Double_t    GetA(Int_t n=0) const { return m_a.at(n); }
  Double_t    GetT(Int_t n=0) const { return m_t.at(n); }
  Double_t    GetCT(Int_t n=0) const { return m_ct.at(n); }
  Double_t    Time(Int_t n=0) const { return GetT(n); }
  Double_t    CTime(Int_t n=0) const { return GetCT(n); }
  virtual
  Double_t    DeltaE(Int_t n=0) const { return GetA(n); }
  Double_t    GetAUp(Int_t n=0) const { return m_a.at(n); }
  Double_t    GetALeft(Int_t n=0) const { return m_a.at(n); }
  Double_t    GetADown(Int_t n=0) const { return m_a.at(n); }
  Double_t    GetARight(Int_t n=0) const { return m_a.at(n); }
  virtual
  Double_t    GetTUp(Int_t n=0) const { return m_t.at(n); }
  Double_t    GetTLeft(Int_t n=0) const { return m_t.at(n); }
  virtual
  Double_t    GetTDown(Int_t n=0) const { return m_t.at(n); }
  Double_t    GetTRight(Int_t n=0) const { return m_t.at(n); }
  Double_t    GetCTUp(Int_t n=0) const { return m_ct.at(n); }
  Double_t    GetCTLeft(Int_t n=0) const { return m_ct.at(n); }
  Double_t    GetCTDown(Int_t n=0) const { return m_ct.at(n); }
  Double_t    GetCTRight(Int_t n=0) const { return m_ct.at(n); }
  virtual
  Double_t    MeanTime(Int_t n=0) const { return GetT(n); }
  virtual
  Double_t    CMeanTime(Int_t n=0) const { return GetCT(n); }
  virtual
  Int_t DetectorId() const { return m_raw->DetectorId(); }
  virtual
  Int_t PlaneId() const { return m_raw->PlaneId(); }
  virtual
  Int_t SegmentId() const { return m_raw->SegmentId(); }

  // For BGO
  void   ClearACont()  { m_a.clear(); }
  void   SetE(Double_t energy)  { m_a.push_back(energy); }
  void   SetJoined(Int_t m)           { m_flag_join.at(m) = true;         }
  Bool_t Joined(Int_t m)        const { return m_flag_join.at(m);         }
  Bool_t JoinedAllMhit();

  virtual Bool_t ReCalc(Bool_t applyRecursively=false)
  { return Calculate(); }
};

//_____________________________________________________________________________
inline TString
Hodo1Hit::ClassName()
{
  static TString s_name("Hodo1Hit");
  return s_name;
}

#endif

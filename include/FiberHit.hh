// -*- C++ -*-

#ifndef FIBER_HIT_HH
#define FIBER_HIT_HH

#include <string>

#include <std_ostream.hh>

#include "Hodo1Hit.hh"

class FLHit;

//_____________________________________________________________________________
class FiberHit : public Hodo1Hit
{
public:
  static const TString& ClassName();
  explicit FiberHit(HodoRawHit* object, const TString& name);
  virtual  ~FiberHit();

private:
  FiberHit();
  FiberHit(const FiberHit& object);
  FiberHit& operator =(const FiberHit& object);

protected:
  TString             m_detector_name;
  Int_t               m_segment;
  Int_t               m_ud;
  Double_t            m_position;
  Double_t            m_offset;
  Int_t               m_pair_id;
  Bool_t              m_status;
  std::vector<FLHit*> m_hit_container;
  Double_t            m_adc_hg;
  Double_t            m_adc_lg;
  Double_t            m_pedcor_hg;
  Double_t            m_pedcor_lg;
  Double_t            m_mip_hg;
  Double_t            m_mip_lg;
  Double_t            m_dE_hg;
  Double_t            m_dE_lg;

  struct data_pair
  {
    Double_t time_l;
    Double_t time_t;
    Double_t ctime_l;
    Double_t tot;
    Int_t    index_t;
  };

  std::vector<data_pair> m_pair_cont;

public:
  Bool_t   Calculate();
  // Call super class method
  Int_t    GetNLeading() const { return Hodo1Hit::GetNumOfHit(0); }
  Int_t    GetNTrailing() const { return Hodo1Hit::GetNumOfHit(1); }
  Double_t GetLeading(Int_t n=0) const
  { return m_ud==0? m_raw->GetTdc1(n) : m_raw->GetTdc2(n); }
  Double_t GetTrailing(Int_t n=0) const
  { return m_ud==0? m_raw->GetTdcT1(n) : m_raw->GetTdcT2(n); }
  // Call member in this class
  Int_t    GetNPair() const { return m_pair_cont.size(); }
  Double_t GetTime(Int_t n=0) const { return m_pair_cont.at(n).time_l; }
  Double_t GetCTime(Int_t n=0) const { return m_pair_cont.at(n).ctime_l; }
  Double_t GetTimeT(Int_t n=0) const { return m_pair_cont.at(n).time_t; }
  Double_t GetWidth(Int_t n=0) const { return m_pair_cont.at(n).tot; }
  Double_t GetTot(Int_t n=0) const { return m_pair_cont.at(n).tot; }
  Double_t GetPosition() const { return m_position + m_offset; }
  Int_t    PairId() const { return m_pair_id; }
  //  virtual Double_t SegmentId()    const { return m_segment; }
  Double_t GetAdcHG() const { return m_adc_hg; }
  Double_t GetAdcLG() const { return m_adc_lg; }
  Double_t GetMipHG() const { return m_mip_hg; }
  Double_t GetMipLG() const { return m_mip_lg; }
  Double_t GetDeHG() const { return m_dE_hg; }
  Double_t GetDeLG() const { return m_dE_lg; }
  void     SetDetectorName(const TString& name) { m_detector_name = name; }
  void     SetPedestalCor(Double_t deltaHG, Double_t deltaLG)
  { m_pedcor_hg = deltaHG; m_pedcor_lg = deltaLG; }
  void     Print(const TString& arg="", std::ostream& ost=hddaq::cout) const;
  void     RegisterHits(FLHit* hit) { m_hit_container.push_back(hit); }
  virtual Bool_t ReCalc(Bool_t allpyRecursively=false)
  { return FiberHit::Calculate(); }
  static Bool_t CompFiberHit(const FiberHit* left, const FiberHit* right);
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
FiberHit::CompFiberHit(const FiberHit* left, const FiberHit* right)
{
  return left->PairId() < right->PairId();
}

#endif

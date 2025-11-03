// -*- C++ -*-

#ifndef RAW_DATA_HH
#define RAW_DATA_HH

#include <map>
#include <vector>

#include <TString.h>

class HodoRawHit;
class DCRawHit;

using HodoRHC = std::vector<HodoRawHit*>;
using DCRHC = std::vector<DCRawHit*>;

//_____________________________________________________________________________
class RawData
{
public:
  static TString ClassName();
  RawData();
  ~RawData();

private:
  RawData(const RawData&);
  RawData& operator=(const RawData&);

private:
  template <typename T> using map_t = std::map<TString, T>;

  map_t<Bool_t>  m_is_decoded;
  map_t<HodoRHC> m_hodo_raw_hit_collection;
  map_t<DCRHC>   m_dc_raw_hit_collection;

public:
  void           Clear(const TString& name="");
  Bool_t         DecodeHits(const TString& name="");
  Bool_t         DecodeCalibHits();
  const HodoRHC& GetHodoRawHitContainer(const TString& name) const;
  const DCRHC&   GetDCRawHitContainer(const TString& name) const;
  void           TdcCutBCOut();
  void           TdcCutSDCIn();
  void           TdcCutSDCOut();
  void           Print(Option_t* arg="") const;

  // aliases
  const HodoRHC& GetHodoRawHC(const TString& name) const
    { return GetHodoRawHitContainer(name); }
  const DCRHC&   GetDCRawHC(const TString& name) const
    { return GetDCRawHitContainer(name); }
  const HodoRHC& GetHodoRawHits(const TString& name) const
    { return GetHodoRawHitContainer(name); }
  const DCRHC&   GetDCRawHits(const TString& name) const
    { return GetDCRawHitContainer(name); }

  // templates
  template <typename T> Int_t GetEntries(const TString& name) const;
  template <typename T> const T* Get(const TString& name, Int_t i) const;

private:
  Bool_t AddHodoRawHit(const TString& name, Int_t plane, Int_t seg,
                       Int_t UorD, Int_t data, Double_t val);
  Bool_t AddFiberRawHit(const TString& name, Int_t plane, Int_t seg,
                        Int_t UorD, Int_t data, Double_t val);
  Bool_t AddDCRawHit(const TString& name, Int_t plane, Int_t seg,
                     Int_t UorD, Int_t data, Double_t val);
  void TdcCut(const TString& name, Double_t min_tdc, Double_t max_tdc);
  void EraseEmptyHits(DCRHC& HitCont);
};

//_____________________________________________________________________________
inline TString
RawData::ClassName()
{
  static TString s_name("RawData");
  return s_name;
}

//_____________________________________________________________________________
template <>
inline Int_t
RawData::GetEntries<HodoRawHit>(const TString& name) const
{
  return m_hodo_raw_hit_collection.at(name).size();
}

//_____________________________________________________________________________
template <>
inline Int_t
RawData::GetEntries<DCRawHit>(const TString& name) const
{
  return m_dc_raw_hit_collection.at(name).size();
}

//_____________________________________________________________________________
template <>
inline const HodoRawHit*
RawData::Get<HodoRawHit>(const TString& name, Int_t i) const
{
  return m_hodo_raw_hit_collection.at(name).at(i);
}

//_____________________________________________________________________________
template <>
inline const DCRawHit*
RawData::Get<DCRawHit>(const TString& name, Int_t i) const
{
  return m_dc_raw_hit_collection.at(name).at(i);
}

#endif

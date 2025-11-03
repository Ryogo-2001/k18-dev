// -*- C++ -*-

#ifndef HODO_RAW_HIT_HH
#define HODO_RAW_HIT_HH

#include <vector>

#include <TString.h>

//_____________________________________________________________________________
class HodoRawHit
{
public:
  static const TString& ClassName();
  HodoRawHit(const TString& detector_name, Int_t plane_id, Int_t segment_id);
  ~HodoRawHit();
  enum EChannel { kUp, kDown, kExtra, kNChannel };

private:
  using data_t = std::vector<std::vector<Double_t>>; // [ch][mhit]
  TString m_detector_name;
  Int_t   m_detector_id;
  TString m_plane_name;
  Int_t   m_plane_id;
  Int_t   m_segment_id;
  data_t  m_adc_high;
  data_t  m_adc_low;
  data_t  m_tdc_leading;
  data_t  m_tdc_trailing;
  Bool_t  m_tdc_is_overflow; // module TDC overflow

public:
  const TString& DetectorName() const { return m_detector_name; }
  Int_t          DetectorId() const { return m_detector_id; }
  const TString& PlaneName() const { return m_plane_name; }
  Int_t          PlaneId() const { return m_plane_id; }
  Int_t          SegmentId() const { return m_segment_id; }

  const std::vector<Double_t>& GetArrayAdcHigh(Int_t i=0) const
    { return m_adc_high.at(i); }
  const std::vector<Double_t>& GetArrayAdcLow(Int_t i=0) const
    { return m_adc_low.at(i); }
  const std::vector<Double_t>& GetArrayTdcLeading(Int_t i=0) const
    { return m_tdc_leading.at(i); }
  const std::vector<Double_t>& GetArrayTdcTrailing(Int_t i=0) const
    { return m_tdc_trailing.at(i); }

  Double_t GetAdcHigh(Int_t i=0, Int_t j=0) const;
  Double_t GetAdcLow(Int_t i=0, Int_t j=0) const;
  Double_t GetTdcLeading(Int_t i=0, Int_t j=0) const;
  Double_t GetTdcTrailing(Int_t i=0, Int_t j=0) const;

  void SetAdcHigh(Int_t i, Double_t adc)
    { m_adc_high.at(i).push_back(adc); }
  void SetAdcLow(Int_t i, Double_t adc)
    { m_adc_low.at(i).push_back(adc); }
  void SetTdcLeading(Int_t i, Double_t tdc)
    { m_tdc_leading.at(i).push_back(tdc); }
  void SetTdcTrailing(Int_t i, Double_t tdc)
    { m_tdc_trailing.at(i).push_back(tdc); }

  Int_t GetSizeAdcHigh(Int_t i=0) const
    { return m_adc_high.at(i).size(); }
  Int_t GetSizeAdcLow(Int_t i=0) const
    { return m_adc_low.at(i).size(); }
  Int_t GetSizeTdcLeading(Int_t i=0) const
    { return m_tdc_leading.at(i).size(); }
  Int_t GetSizeTdcTrailing(Int_t i=0) const
    { return m_tdc_trailing.at(i).size(); }

  // aliases
  const std::vector<Double_t>& GetArrayAdc(Int_t i=0) const
    { return GetArrayAdcHigh(i); }
  const std::vector<Double_t>& GetArrayAdcUp() const
    { return GetArrayAdcHigh(kUp); }
  const std::vector<Double_t>& GetArrayAdcDown() const
    { return GetArrayAdcHigh(kDown); }
  const std::vector<Double_t>& GetArrayAdcExtra() const
    { return GetArrayAdcHigh(kExtra); }
  const std::vector<Double_t>& GetArrayTdc(Int_t i=0) const
    { return GetArrayTdcLeading(i); }
  const std::vector<Double_t>& GetArrayTdcUp() const
    { return GetArrayTdcLeading(kUp); }
  const std::vector<Double_t>& GetArrayTdcDown() const
    { return GetArrayTdcLeading(kDown); }
  const std::vector<Double_t>& GetArrayTdcLeft() const
    { return GetArrayTdcLeading(kUp); }
  const std::vector<Double_t>& GetArrayTdcRight() const
    { return GetArrayTdcLeading(kDown); }
  const std::vector<Double_t>& GetArrayTdcExtra() const
    { return GetArrayTdcLeading(kExtra); }
  Double_t GetAdc(Int_t i=0, Int_t j=0) const
    { return GetAdcHigh(i, j); }
  Double_t GetTdc(Int_t i=0, Int_t j=0) const
    { return GetTdcLeading(i, j); }

  Double_t GetAdcUp(Int_t i=0) const { return GetAdc(kUp, i); }
  Double_t GetAdcLeft(Int_t i=0) const { return GetAdc(kUp, i); }
  Double_t GetAdcDown(Int_t i=0) const { return GetAdc(kDown, i); }
  Double_t GetAdcRight(Int_t i=0) const { return GetAdc(kDown, i); }
  Double_t GetAdcExtra(Int_t i=0) const { return GetAdc(kExtra, i); }

  Double_t GetTdcUp(Int_t i=0) const { return GetTdc(kUp, i); }
  Double_t GetTdcLeft(Int_t i=0) const { return GetTdc(kUp, i); }
  Double_t GetTdcDown(Int_t i=0) const { return GetTdc(kDown, i); }
  Double_t GetTdcRight(Int_t i=0) const { return GetTdc(kDown, i); }
  Double_t GetTdcExtra(Int_t i=0) const { return GetTdc(kExtra, i); }

  Int_t GetSizeTdcUp() const { return GetSizeTdcLeading(kUp); }
  Int_t GetSizeTdcDown() const { return GetSizeTdcLeading(kDown); }
  Int_t GetSizeTdcExtra() const { return GetSizeTdcLeading(kExtra); }

  void SetAdc(Int_t i, Double_t adc) { SetAdcHigh(i, adc); }
  void SetTdc(Int_t i, Double_t tdc) { SetTdcLeading(i, tdc); }

  void SetAdcUp(Double_t adc){ SetAdc(kUp, adc); }
  void SetAdcLeft(Double_t adc){ SetAdc(kUp, adc); }
  void SetAdcDown(Double_t adc){ SetAdc(kDown, adc); }
  void SetAdcRight(Double_t adc){ SetAdc(kDown, adc); }
  void SetAdcExtra(Double_t adc){ SetAdc(kExtra, adc); }

  void SetTdcUp(Double_t tdc){ SetTdc(kUp, tdc); }
  void SetTdcLeft(Double_t tdc){ SetTdc(kUp, tdc); }
  void SetTdcDown(Double_t tdc){ SetTdc(kDown, tdc); }
  void SetTdcRight(Double_t tdc){ SetTdc(kDown, tdc); }
  void SetTdcExtra(Double_t tdc){ SetTdc(kExtra, tdc); }

  // for Multi-hit method
  Bool_t GetTdcOverflow() const { return m_tdc_is_overflow; }
  void   SetTdcOverflow(Int_t i, Double_t flag)
    { m_tdc_is_overflow = static_cast<Bool_t>(flag); }

  void Clear();
  void Print(const TString& arg="") const;
};

//_____________________________________________________________________________
inline const TString&
HodoRawHit::ClassName()
{
  static TString s_name("HodoRawHit");
  return s_name;
}

#endif

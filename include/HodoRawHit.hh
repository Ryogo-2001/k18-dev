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
  HodoRawHit(Int_t detector_id, Int_t plane_id, Int_t segment_id);
  ~HodoRawHit();

private:
  Int_t              m_detector_id;
  Int_t              m_plane_id;
  Int_t              m_segment_id;
  std::vector<Int_t> m_adc1;
  std::vector<Int_t> m_adc2;
  // leading
  std::vector<Int_t> m_tdc1;
  std::vector<Int_t> m_tdc2;
  //trailing
  std::vector<Int_t> m_tdc_t1;
  std::vector<Int_t> m_tdc_t2;
  Bool_t             m_oftdc; // module TDC overflow
  Int_t              m_nhtdc;

public:
  void SetAdc1(Int_t adc);
  void SetAdc2(Int_t adc);
  // leading
  void SetTdc1(Int_t tdc);
  void SetTdc2(Int_t tdc);
  // trailing
  void SetTdcT1(Int_t tdc);
  void SetTdcT2(Int_t tdc);

  void SetAdcUp(Int_t adc){ SetAdc1(adc); }
  void SetAdcLeft(Int_t adc){ SetAdc1(adc); }
  void SetAdcDown(Int_t adc){ SetAdc2(adc); }
  void SetAdcRight(Int_t adc){ SetAdc2(adc); }

  // leading
  void SetTdcUp(Int_t tdc){ SetTdc1(tdc); }
  void SetTdcLeft(Int_t tdc){ SetTdc1(tdc); }
  void SetTdcDown(Int_t tdc){ SetTdc2(tdc); }
  void SetTdcRight(Int_t tdc){ SetTdc2(tdc); }
  // trailing
  void SetTdcTUp(Int_t tdc){ SetTdcT1(tdc); }
  void SetTdcTLeft(Int_t tdc){ SetTdcT1(tdc); }
  void SetTdcTDown(Int_t tdc){ SetTdcT2(tdc); }
  void SetTdcTRight(Int_t tdc){ SetTdcT2(tdc); }

  Int_t DetectorId() const { return m_detector_id; }
  Int_t PlaneId() const { return m_plane_id; }
  Int_t SegmentId() const { return m_segment_id; }
  // for Multi-hit method
  void  SetTdcOverflow(Int_t flag) { m_oftdc = static_cast<Bool_t>(flag); }
  Int_t GetNumOfTdcHits() const { return m_nhtdc; }
  const std::vector<Int_t>& GetArrayAdc1() const { return m_adc1; }
  const std::vector<Int_t>& GetArrayAdc2() const { return m_adc2; }
  Int_t GetAdc1(Int_t i=0) const { return m_adc1.at(i); }
  Int_t GetAdc2(Int_t i=0) const { return m_adc2.at(i); }
  // leading
  const std::vector<Int_t>& GetArrayTdc1() const { return m_tdc1; }
  const std::vector<Int_t>& GetArrayTdc2() const { return m_tdc2; }
  Int_t GetTdc1(Int_t i=0) const { return m_tdc1.at(i); }
  Int_t GetTdc2(Int_t i=0) const { return m_tdc2.at(i); }
  // trailing
  const std::vector<Int_t>& GetArrayTdcT1() const { return m_tdc_t1; }
  const std::vector<Int_t>& GetArrayTdcT2() const { return m_tdc_t2; }
  Int_t GetTdcT1(Int_t i=0) const { return m_tdc_t1.at(i); }
  Int_t GetTdcT2(Int_t i=0) const { return m_tdc_t2.at(i); }

  Int_t GetAdcUp(Int_t i=0) const { return GetAdc1(i); }
  Int_t GetAdcLeft(Int_t i=0) const { return GetAdc1(i); }
  Int_t GetAdcDown(Int_t i=0) const { return GetAdc2(i); }
  Int_t GetAdcRight(Int_t i=0) const { return GetAdc2(i); }

  // leading
  Int_t GetTdcUp(Int_t i=0) const { return GetTdc1(i); }
  Int_t GetTdcLeft(Int_t i=0) const { return GetTdc1(i); }
  Int_t GetTdcDown(Int_t i=0) const { return GetTdc2(i); }
  Int_t GetTdcRight(Int_t i=0) const { return GetTdc2(i); }
  // trailing
  Int_t GetTdcTUp(Int_t i=0) const { return GetTdcT1(i); }
  Int_t GetTdcTLeft(Int_t i=0) const { return GetTdcT1(i); }
  Int_t GetTdcTDown(Int_t i=0) const { return GetTdcT2(i); }
  Int_t GetTdcTRight(Int_t i=0) const { return GetTdcT2(i); }

  Int_t SizeAdc1() const;
  Int_t SizeAdc2() const;
  Int_t SizeTdc1() const;
  Int_t SizeTdc2() const;
  Int_t SizeTdcT1() const;
  Int_t SizeTdcT2() const;

  Int_t GetSizeAdcUp() const { return SizeAdc1(); }
  Int_t GetSizeAdcLeft() const { return SizeAdc1(); }
  Int_t GetSizeAdcDown() const { return SizeAdc2(); }
  Int_t GetSizeAdcRight() const { return SizeAdc2(); }

  // leading
  Int_t GetSizeTdcUp() const { return SizeTdc1(); }
  Int_t GetSizeTdcLeft() const { return SizeTdc1(); }
  Int_t GetSizeTdcDown() const { return SizeTdc2(); }
  Int_t GetSizeTdcRight() const { return SizeTdc2(); }

  // trailing
  Int_t GetSizeTdcTUp() const { return SizeTdcT1(); }
  Int_t GetSizeTdcTLeft() const { return SizeTdcT1(); }
  Int_t GetSizeTdcTDown() const { return SizeTdcT2(); }
  Int_t GetSizeTdcTRight() const { return SizeTdcT2(); }

  Bool_t IsTdcOverflow() const { return m_oftdc; }
  void Clear();
  void Print(const TString& arg="");
};

//_____________________________________________________________________________
inline const TString&
HodoRawHit::ClassName()
{
  static TString s_name("HodoRawHit");
  return s_name;
}

#endif

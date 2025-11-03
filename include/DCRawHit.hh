// -*- C++ -*-

#ifndef DC_RAW_HIT_HH
#define DC_RAW_HIT_HH

#include <TString.h>

//_____________________________________________________________________________
class DCRawHit
{
public:
  static const TString& ClassName();
  DCRawHit(const TString& detector_name, Int_t plane_id, Int_t wire_id);
  ~DCRawHit();

private:
  using data_t = std::vector<Double_t>; // [mhit]
  TString m_detector_name;
  Int_t   m_detector_id;
  TString m_plane_name;
  Int_t   m_plane_id;
  Int_t   m_dcgeom_layer;
  Int_t   m_wire_id;
  data_t  m_tdc;
  data_t  m_trailing;
  Bool_t  m_oftdc; // module tdc over flow

public:
  const auto& DetectorName() const { return m_detector_name; }
  Int_t       DetectorId() const { return m_detector_id; }
  const auto& PlaneName() const { return m_plane_name; }
  Int_t       PlaneId() const { return m_plane_id; }
  Int_t       DCGeomLayerId() const { return m_dcgeom_layer; }
  Int_t       LayerId() const { return DCGeomLayerId(); }
  Int_t       WireId() const { return m_wire_id; }
  void        SetTdc(Int_t tdc) { m_tdc.push_back(tdc); }
  void        SetTrailing(Int_t tdc) { m_trailing.push_back(tdc); }
  void        SetTdcOverflow(Int_t fl) { m_oftdc = static_cast<Bool_t>(fl); }
  const auto& GetTdcArray() const { return m_tdc; }
  Int_t       GetTdc(Int_t nh) const { return m_tdc[nh]; }
  Int_t       GetTdcSize() const { return m_tdc.size(); }
  const auto& GetTrailingArray() const { return m_trailing; }
  Int_t       GetTrailing(Int_t nh) const { return m_trailing[nh]; }
  Int_t       GetTrailingSize() const { return m_trailing.size(); }
  Bool_t      IsTdcOverflow() const { return m_oftdc; }
  Bool_t      IsEmpty() const { return GetTdcSize() == 0; }
  void        TdcCut(Double_t min, Double_t max);
  void        Print(const TString& arg="") const;
};

//_____________________________________________________________________________
inline const TString&
DCRawHit::ClassName()
{
  static TString s_name("DCRawHit");
  return s_name;
}

#endif

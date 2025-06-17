// -*- C++ -*-

#ifndef DC_RAW_HIT_HH
#define DC_RAW_HIT_HH

#include <vector>

#include <TString.h>

typedef std::vector<Int_t> IntVec;

//_____________________________________________________________________________
class DCRawHit
{
public:
  static const TString& ClassName();
  DCRawHit(Int_t plane_id, Int_t wire_id);
  ~DCRawHit();

private:
  Int_t  m_plane_id;
  Int_t  m_wire_id;
  IntVec m_tdc;
  IntVec m_trailing;
  Bool_t m_oftdc; // module tdc over flow

public:
  Int_t  PlaneId() const { return m_plane_id; }
  Int_t  WireId() const { return m_wire_id; }
  void   SetTdc(Int_t tdc) { m_tdc.push_back(tdc); }
  void   SetTrailing(Int_t tdc) { m_trailing.push_back(tdc); }
  void   SetTdcOverflow(Int_t fl) { m_oftdc = static_cast<Bool_t>(fl); }
  Int_t  GetTdc(Int_t nh) const { return m_tdc[nh]; }
  Int_t  GetTdcSize() const { return m_tdc.size(); }
  Int_t  GetTrailing(Int_t nh) const { return m_trailing[nh]; }
  Int_t  GetTrailingSize() const { return m_trailing.size(); }
  Bool_t IsTdcOverflow() const { return m_oftdc; }
  void   Print(const TString& arg="") const;
};

//_____________________________________________________________________________
inline const TString&
DCRawHit::ClassName()
{
  static TString s_name("DCRawHit");
  return s_name;
}

#endif

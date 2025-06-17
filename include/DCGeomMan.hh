// -*- C++ -*-

#ifndef DC_GEOM_MAN_HH
#define DC_GEOM_MAN_HH

#include "ThreeVector.hh"
#include <TString.h>
#include <string>
#include <vector>
#include <map>

class DCGeomRecord;

//_____________________________________________________________________________
class DCGeomMan
{
public:
  static const TString& ClassName();
  static DCGeomMan& GetInstance();
  ~DCGeomMan();

private:
  DCGeomMan();
  DCGeomMan(const DCGeomMan&);
  DCGeomMan& operator =(const DCGeomMan&);

private:
  typedef std::map <Int_t, DCGeomRecord*> DCGeomContainer;
  typedef DCGeomContainer::const_iterator DCGeomIterator;
  typedef std::map <TString, Int_t>    IntList;
  typedef std::map <TString, Double_t> DoubleList;
  Bool_t          m_is_ready;
  TString         m_file_name;
  DCGeomContainer m_container;
  IntList         m_detector_id_map;
  DoubleList      m_global_z_map;
  DoubleList      m_local_z_map;

public:
  Bool_t              Initialize();
  Bool_t              Initialize(const TString& file_name);
  Bool_t              IsReady() const { return m_is_ready; }
  Double_t            GetLocalZ(Int_t lnum) const;
  Double_t            GetLocalZ(const TString& key) const;
  Double_t            GetResolution(Int_t lnum) const;
  Double_t            GetResolution(const TString& key) const;
  Double_t            GetTiltAngle(Int_t lnum) const;
  Double_t            GetTiltAngle(const TString& key) const;
  Double_t            GetRotAngle1(Int_t lnum) const;
  Double_t            GetRotAngle1(const TString& key) const;
  Double_t            GetRotAngle2(Int_t lnum) const;
  Double_t            GetRotAngle2(const TString& key) const;
  const ThreeVector&  GetGlobalPosition(Int_t lnum) const;
  const ThreeVector&  GetGlobalPosition(const TString& key) const;
  ThreeVector         NormalVector(Int_t lnum) const;
  ThreeVector         NormalVector(const TString& key) const;
  ThreeVector         UnitVector(Int_t lnum) const;
  ThreeVector         UnitVector(const TString& key) const;
  const DCGeomRecord* GetRecord(Int_t lnum) const;
  const DCGeomRecord* GetRecord(const TString& key) const;
  ThreeVector         Local2GlobalPos(Int_t lnum, const ThreeVector &in) const;
  ThreeVector         Local2GlobalPos(const TString& key, const ThreeVector &in) const;
  ThreeVector         Global2LocalPos(Int_t lnum, const ThreeVector &in) const;
  ThreeVector         Global2LocalPos(ThreeVector posVector, Double_t TA,
				      Double_t RA1, Double_t RA2, const ThreeVector& in) const;
  ThreeVector         Global2LocalPos(const TString& key, const ThreeVector &in) const;
  ThreeVector         Local2GlobalDir(Int_t lnum, const ThreeVector &in) const;
  ThreeVector         Local2GlobalDir(const TString& key, const ThreeVector &in) const;
  ThreeVector         Global2LocalDir(Int_t lnum, const ThreeVector &in) const;
  ThreeVector         Global2LocalDir(const TString& key, const ThreeVector &in) const;
  Double_t            CalcWirePosition(Int_t lnum, Double_t wire) const;
  Double_t            CalcWirePosition(const TString& key, Double_t wire) const;
  Int_t               CalcWireNumber(Int_t lnum, Double_t position) const;
  Int_t               CalcWireNumber(const TString& key, Double_t position) const;
  std::vector<Int_t>  GetDetectorIDList() const;
  Int_t               GetDetectorId(const TString &key) const;
  void                SetFileName(const TString &file_name);
  // Do not use this method except for special cases
  void                SetResolution(Int_t lnum, Double_t res);
  void                SetResolution(const TString& key, Double_t res);

  // Static method
  static const Int_t&    DetectorId(const TString& key);
  static const Double_t& GlobalZ(const TString& key);
  static const Double_t& LocalZ(const TString& key);

private:
  void ClearElements();
};

//______________________________________________________________________________
inline const TString&
DCGeomMan::ClassName()
{
  static TString s_name("DCGeomMan");
  return s_name;
}

//______________________________________________________________________________
inline DCGeomMan&
DCGeomMan::GetInstance()
{
  static DCGeomMan s_instance;
  return s_instance;
}

//______________________________________________________________________________
inline const Int_t&
DCGeomMan::DetectorId(const TString& key)
{
  return GetInstance().m_detector_id_map[key];
}

//______________________________________________________________________________
inline const Double_t&
DCGeomMan::GlobalZ(const TString& key)
{
  return GetInstance().m_global_z_map[key];
}

//______________________________________________________________________________
inline const Double_t&
DCGeomMan::LocalZ(const TString& key)
{
  return GetInstance().m_local_z_map[key];
}

#endif

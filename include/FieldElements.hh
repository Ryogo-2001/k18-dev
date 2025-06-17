// -*- C++ -*-

#ifndef FIELD_ELEMENTS_HH
#define FIELD_ELEMENTS_HH

#include <TString.h>
#include <TVector3.h>

#include "DCGeomRecord.hh"

//_____________________________________________________________________________
class FieldElements
{
public:
  static const TString& ClassName();
  FieldElements(const TString& name, const TVector3& pos,
                Double_t ta, Double_t ra1, Double_t ra2);
  virtual ~FieldElements();

private:
  enum FERegion { kFERSurface, kFERInside, kFEROutside };
  DCGeomRecord m_geom_record;

public:
  static FERegion FERSurface() { return kFERSurface; }
  static FERegion FERInside() { return kFERInside;  }
  static FERegion FEROutside() { return kFEROutside; }
  TVector3        Local2GlobalPos(const TVector3& in) const;
  TVector3        Local2GlobalDir(const TVector3& in) const;
  TVector3        Global2LocalPos(const TVector3& in) const;
  TVector3        Global2LocalDir(const TVector3& in) const;
  virtual TVector3 GetField(const TVector3& gPos) const = 0;
  virtual Bool_t   ExistField(const TVector3& gPos) const = 0;
  virtual FERegion CheckRegion(const TVector3& gPos,
                               Double_t Tolerance) const = 0;
};

//_____________________________________________________________________________
inline const TString&
FieldElements::ClassName()
{
  static TString s_name("FieldElements");
  return s_name;
}

#endif

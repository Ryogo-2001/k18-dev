// -*- C++ -*-

#ifndef DC_GEOM_RECORD_HH
#define DC_GEOM_RECORD_HH


#include <string>
#include <functional>

#include <TString.h>
#include <TVector3.h>

#include <std_ostream.hh>

//_____________________________________________________________________________
class DCGeomRecord
{
public:
  DCGeomRecord(Int_t id, const TString& name,
               Double_t x, Double_t y, Double_t z, Double_t ta,
               Double_t ra1, Double_t ra2, Double_t length, Double_t resol,
               Double_t w0, Double_t dd, Double_t ofs);
  DCGeomRecord(Int_t id, const TString& name,
               const TVector3& pos, Double_t ta,
               Double_t ra1, Double_t ra2, Double_t length, Double_t resol,
               Double_t w0, Double_t dd, Double_t ofs);
  ~DCGeomRecord();
  DCGeomRecord(const DCGeomRecord&);
  DCGeomRecord& operator =(const DCGeomRecord);

private:
  Int_t    m_id;
  TString  m_name;
  TVector3 m_pos;
  Double_t m_tilt_angle;
  Double_t m_rot_angle1;
  Double_t m_rot_angle2;
  Double_t m_length;
  Double_t m_resolution;
  Double_t m_w0;
  Double_t m_dd;
  Double_t m_offset;

  Double_t m_dxds, m_dxdt, m_dxdu;
  Double_t m_dyds, m_dydt, m_dydu;
  Double_t m_dzds, m_dzdt, m_dzdu;

  Double_t m_dsdx, m_dsdy, m_dsdz;
  Double_t m_dtdx, m_dtdy, m_dtdz;
  Double_t m_dudx, m_dudy, m_dudz;

public:
  const TVector3& Position() const { return m_pos; }
  TVector3        NormalVector() const;
  TVector3        UnitVector() const;
  Int_t           Id() const { return m_id; }
  TString         Name() const { return m_name; }
  const TVector3& Pos() const { return m_pos; }
  Double_t        TiltAngle() const { return m_tilt_angle; }
  Double_t        RotationAngle1() const { return m_rot_angle1; }
  Double_t        RotationAngle2() const { return m_rot_angle2; }
  Double_t        dXdW() const { return m_dd; }
  Double_t        Length() const { return m_length; }
  Double_t        Resolution() const { return m_resolution; }
  void            SetResolution(Double_t res) { m_resolution = res; }

  Double_t dsdx() const { return m_dsdx; }
  Double_t dsdy() const { return m_dsdy; }
  Double_t dsdz() const { return m_dsdz; }
  Double_t dtdx() const { return m_dtdx; }
  Double_t dtdy() const { return m_dtdy; }
  Double_t dtdz() const { return m_dtdz; }
  Double_t dudx() const { return m_dudx; }
  Double_t dudy() const { return m_dudy; }
  Double_t dudz() const { return m_dudz; }

  Double_t dxds() const { return m_dxds; }
  Double_t dxdt() const { return m_dxdt; }
  Double_t dxdu() const { return m_dxdu; }
  Double_t dyds() const { return m_dyds; }
  Double_t dydt() const { return m_dydt; }
  Double_t dydu() const { return m_dydu; }
  Double_t dzds() const { return m_dzds; }
  Double_t dzdt() const { return m_dzdt; }
  Double_t dzdu() const { return m_dzdu; }

  Double_t WirePos(Double_t wire)   const;
  Int_t    WireNumber(Double_t pos) const;
  void     Print(const TString& arg="", std::ostream& ost=hddaq::cout) const;

private:
  void CalcVectors();

};

//_____________________________________________________________________________
struct DCGeomRecordComp
  : public std::binary_function <DCGeomRecord*, DCGeomRecord*, bool>
{
  bool operator()(const DCGeomRecord* const p1,
                  const DCGeomRecord* const p2) const
    { return p1->Id() < p2->Id(); }
};

#endif

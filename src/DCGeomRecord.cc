// -*- C++ -*-

#include "DCGeomRecord.hh"

#include <cstdlib>
#include <iostream>

#include "MathTools.hh"
#include "PrintHelper.hh"

namespace
{
enum  EDefinition { kSks, kKurama };
const EDefinition GlobalCoordinate = kKurama;
}

//_____________________________________________________________________________
DCGeomRecord::DCGeomRecord(Int_t id, const TString& name,
                           Double_t x, Double_t y, Double_t z, Double_t ta,
                           Double_t ra1, Double_t ra2, Double_t length,
                           Double_t resol,
                           Double_t w0, Double_t dd, Double_t ofs)
  : m_id(id),
    m_name(name),
    m_pos(x, y, z),
    m_tilt_angle(ta),
    m_rot_angle1(ra1),
    m_rot_angle2(ra2),
    m_length(length),
    m_resolution(resol),
    m_w0(w0),
    m_dd(dd),
    m_offset(ofs)
{
  CalcVectors();
}

//_____________________________________________________________________________
DCGeomRecord::DCGeomRecord(Int_t id, const TString& name,
                           const TVector3& pos, Double_t ta,
                           Double_t ra1, Double_t ra2, Double_t length,
                           Double_t resol,
                           Double_t w0, Double_t dd, Double_t ofs)
  : m_id(id),
    m_name(name),
    m_pos(pos),
    m_tilt_angle(ta),
    m_rot_angle1(ra1),
    m_rot_angle2(ra2),
    m_length(length),
    m_resolution(resol),
    m_w0(w0),
    m_dd(dd),
    m_offset(ofs)
{
  CalcVectors();
}

//_____________________________________________________________________________
DCGeomRecord::~DCGeomRecord()
{
}

//_____________________________________________________________________________
void
DCGeomRecord::CalcVectors()
{
  Double_t ct0 = TMath::Cos(m_tilt_angle*TMath::DegToRad());
  Double_t st0 = TMath::Sin(m_tilt_angle*TMath::DegToRad());
  Double_t ct1 = TMath::Cos(m_rot_angle1*TMath::DegToRad());
  Double_t st1 = TMath::Sin(m_rot_angle1*TMath::DegToRad());
  Double_t ct2 = TMath::Cos(m_rot_angle2*TMath::DegToRad());
  Double_t st2 = TMath::Sin(m_rot_angle2*TMath::DegToRad());

  switch(GlobalCoordinate){
  case kSks: // SKS difinition
    m_dxds =  ct0*ct2-st0*ct1*st2;
    m_dxdt = -st0*ct2-ct0*ct1*st2;
    m_dxdu =  st1*st2;

    m_dyds =  ct0*st2+st0*ct1*ct2;
    m_dydt = -st0*st2+ct0*ct1*ct2;
    m_dydu = -st1*ct2;

    m_dzds =  st0*st1;
    m_dzdt =  ct0*st1;
    m_dzdu =  ct1;

    m_dsdx =  ct0*ct2-st0*ct1*st2;
    m_dsdy =  ct0*st2+st0*ct1*ct2;
    m_dsdz =  st0*st1;

    m_dtdx = -st0*ct2-ct0*ct1*st2;
    m_dtdy = -st0*st2+ct0*ct1*ct2;
    m_dtdz =  ct0*st1;

    m_dudx =  st1*st2;
    m_dudy = -st1*ct2;
    m_dudz =  ct1;
    break;
  case kKurama: // new difinition : RA2 is around Y-axis
    m_dxds =  ct0*ct2+st0*st1*st2;
    m_dxdt = -st0*ct2+ct0*st1*st2;
    m_dxdu =  ct1*st2;

    m_dyds =  st0*ct1;
    m_dydt =  ct0*ct1;
    m_dydu = -st1;

    m_dzds = -ct0*st2+st0*st1*ct2;
    m_dzdt =  st0*st2+ct0*st1*ct2;
    m_dzdu =  ct1*ct2;

    m_dsdx =  ct0*ct2+st0*st1*st2;
    m_dsdy =  st0*ct1;
    m_dsdz = -ct0*st2+st0*st1*ct2;

    m_dtdx = -st0*ct2+ct0*st1*st2;
    m_dtdy =  ct0*ct1;
    m_dtdz =  st0*st2+ct0*st1*ct2;

    m_dudx =  ct1*st2;
    m_dudy = -st1;
    m_dudz =  ct1*ct2;
    break;
  }
}

//_____________________________________________________________________________
TVector3
DCGeomRecord::NormalVector() const
{
  return TVector3(m_dxdu, m_dydu, m_dzdu);
}

//_____________________________________________________________________________
TVector3
DCGeomRecord::UnitVector() const
{
  return TVector3(m_dxds, m_dyds, m_dzds);
}

//_____________________________________________________________________________
Double_t
DCGeomRecord::WirePos(Double_t wire) const
{
  return m_dd*(wire - m_w0)+m_offset;
}

//_____________________________________________________________________________
Int_t
DCGeomRecord::WireNumber(Double_t pos) const
{
  Double_t dw = ((pos-m_offset)/m_dd) + m_w0;
  Int_t    iw = Int_t(dw);
  if((dw-Double_t(iw))>0.5)
    return iw+1;
  else
    return iw;
}

//_____________________________________________________________________________
void
DCGeomRecord::Print(const TString& arg, std::ostream& ost) const
{
  PrintHelper helper(3, std::ios::fixed, ost);
  ost << " id = "   << std::setw(3) << std::right << m_id << " "
      << " name = " << std::setw(9) << std::left << m_name << " "
      << " pos = (" << " " << std::right << std::fixed
      << std::setw(8)  << m_pos.x() << ", "
      << std::setw(8)  << m_pos.y() << ", "
      << std::setw(10) << m_pos.z() << ") "
      << " TA = "     << std::setw(7) << m_tilt_angle << " "
      << " RA1 = "    << std::setw(6) << m_rot_angle1 << " "
      << " RA2 = "    << std::setw(6) << m_rot_angle2 << " "
      << " L = "      << std::setw(9) << m_length << " "
      << " W0 = "     << std::setw(7) << m_w0 << " "
      << " dWdX = "   << std::setw(6) << m_dd << " "
      << " offset = " << std::setw(7) << m_offset << std::endl;
}

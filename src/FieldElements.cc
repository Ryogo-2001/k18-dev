// -*- C++ -*-

#include "FieldElements.hh"

//_____________________________________________________________________________
FieldElements::FieldElements(const TString& name,
                             const TVector3& pos,
                             Double_t ta, Double_t ra1, Double_t ra2)
  : m_geom_record(-1, name, pos, ta, ra1, ra2, 0., 0., 0., 0., 0.)
{
}

//_____________________________________________________________________________
FieldElements::~FieldElements()
{
}

//_____________________________________________________________________________
TVector3
FieldElements::Local2GlobalPos(const TVector3& in) const
{
  TVector3 gPos = m_geom_record.Position();
  Double_t x = gPos.x()
    + m_geom_record.dxds()*in.x()
    + m_geom_record.dxdt()*in.y()
    + m_geom_record.dxdu()*in.z();
  Double_t y = gPos.y()
    + m_geom_record.dyds()*in.x()
    + m_geom_record.dydt()*in.y()
    + m_geom_record.dydu()*in.z();
  Double_t z = gPos.z()
    + m_geom_record.dzds()*in.x()
    + m_geom_record.dzdt()*in.y()
    + m_geom_record.dzdu()*in.z();

  return TVector3(x, y, z);
}

//_____________________________________________________________________________
TVector3
FieldElements::Local2GlobalDir(const TVector3& in) const
{
  Double_t x
    = m_geom_record.dxds()*in.x()
    + m_geom_record.dxdt()*in.y()
    + m_geom_record.dxdu()*in.z();
  Double_t y
    = m_geom_record.dyds()*in.x()
    + m_geom_record.dydt()*in.y()
    + m_geom_record.dydu()*in.z();
  Double_t z
    = m_geom_record.dzds()*in.x()
    + m_geom_record.dzdt()*in.y()
    + m_geom_record.dzdu()*in.z();

  return TVector3(x, y, z);
}

//_____________________________________________________________________________
TVector3
FieldElements::Global2LocalPos(const TVector3& in) const
{
  TVector3 gPos = m_geom_record.Position();
  Double_t x
    = m_geom_record.dsdx()*(in.x()-gPos.x())
    + m_geom_record.dsdy()*(in.y()-gPos.y())
    + m_geom_record.dsdz()*(in.z()-gPos.z());
  Double_t y
    = m_geom_record.dtdx()*(in.x()-gPos.x())
    + m_geom_record.dtdy()*(in.y()-gPos.y())
    + m_geom_record.dtdz()*(in.z()-gPos.z());
  Double_t z
    = m_geom_record.dudx()*(in.x()-gPos.x())
    + m_geom_record.dudy()*(in.y()-gPos.y())
    + m_geom_record.dudz()*(in.z()-gPos.z());

  return TVector3(x, y, z);
}

//_____________________________________________________________________________
TVector3
FieldElements::Global2LocalDir(const TVector3& in) const
{
  Double_t x
    = m_geom_record.dsdx()*in.x()
    + m_geom_record.dsdy()*in.y()
    + m_geom_record.dsdz()*in.z();
  Double_t y
    = m_geom_record.dtdx()*in.x()
    + m_geom_record.dtdy()*in.y()
    + m_geom_record.dtdz()*in.z();
  Double_t z
    = m_geom_record.dudx()*in.x()
    + m_geom_record.dudy()*in.y()
    + m_geom_record.dudz()*in.z();

  return TVector3(x, y, z);
}

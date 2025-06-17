// -*- C++ -*-

#include "BH2Hit.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <std_ostream.hh>

#include "FuncName.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"

namespace
{
const auto& gHodo = HodoParamMan::GetInstance();
// const auto& gPHC = HodoPHCMan::GetInstance();
}

//_____________________________________________________________________________
BH2Hit::BH2Hit(HodoRawHit *rhit, double max_time_diff)
  : Hodo2Hit(rhit, max_time_diff),
    m_time_offset(0.)
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
BH2Hit::~BH2Hit()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
bool
BH2Hit::Calculate()
{
  if(!Hodo2Hit::Calculate()) return false;

  Int_t cid  = m_raw->DetectorId();
  Int_t plid = m_raw->PlaneId();
  Int_t seg  = m_raw->SegmentId();
  gHodo.GetTime(cid, plid, seg, 2, 0, m_time_offset);

  return true;
}

//_____________________________________________________________________________
Double_t
BH2Hit::UTime0(Int_t i) const
{
  return m_pair_cont[i].time1 + m_time_offset;
}

//_____________________________________________________________________________
Double_t
BH2Hit::UCTime0(Int_t i) const
{
  return m_pair_cont[i].ctime1 + m_time_offset;
}

//_____________________________________________________________________________
Double_t
BH2Hit::DTime0(Int_t i) const
{
  return m_pair_cont[i].time2 + m_time_offset;
}

//_____________________________________________________________________________
Double_t
BH2Hit::DCTime0(Int_t i) const
{
  return m_pair_cont[i].ctime2 + m_time_offset;
}

//_____________________________________________________________________________
Double_t
BH2Hit::Time0(Int_t i) const
{
  return 0.5*(m_pair_cont[i].time1 + m_pair_cont[i].time2) + m_time_offset;
}

//_____________________________________________________________________________
Double_t
BH2Hit::CTime0(Int_t i) const
{
  return 0.5*(m_pair_cont[i].ctime1 + m_pair_cont[i].ctime2) + m_time_offset;
}

// -*- C++ -*-

#include "DCTdcCalibMan.hh"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <std_ostream.hh>

#include "DeleteUtility.hh"
#include "FuncName.hh"

namespace
{
inline Int_t
MakeKey(Int_t plane_id, Double_t wire_id)
{
  return (plane_id<<10) | Int_t(wire_id);
}
}

//_____________________________________________________________________________
struct DCTdcCalMap
{
  Double_t p0, p1;
  DCTdcCalMap(Double_t q0, Double_t q1)
    : p0(q0), p1(q1)
    {}
};

//_____________________________________________________________________________
DCTdcCalibMan::DCTdcCalibMan()
  : m_is_ready(false),
    m_file_name(),
    m_container()
{
}

//_____________________________________________________________________________
DCTdcCalibMan::~DCTdcCalibMan()
{
  ClearElements();
}

//_____________________________________________________________________________
void
DCTdcCalibMan::ClearElements()
{
  del::ClearMap(m_container);
}

//_____________________________________________________________________________
Bool_t
DCTdcCalibMan::Initialize()
{
  if(m_is_ready){
    hddaq::cerr << FUNC_NAME << " already initialied" << std::endl;
    return false;
  }

  std::ifstream ifs(m_file_name);
  if(!ifs.is_open()){
    hddaq::cerr << FUNC_NAME << " "
		<< "file open fail : " << m_file_name << std::endl;
    return false;
  }

  ClearElements();

  TString line;
  while(ifs.good() && line.ReadLine(ifs)){
    if(line.IsNull() || line[0]=='#') continue;
    std::istringstream iss(line.Data());
    Int_t plane_id=-1, wire_id=-1;
    Double_t p0=-9999., p1=-9999.;
    if(iss >> plane_id >> wire_id >> p1 >> p0){
      Int_t key = MakeKey(plane_id, wire_id);
      DCTdcCalMap* tdc_calib = new DCTdcCalMap(p0, p1);
      if(m_container[key]) delete m_container[key];
      m_container[key] = tdc_calib;
    } else {
      hddaq::cerr << FUNC_NAME << ": Bad format => "
		  << line << std::endl;
    }
  }

  m_is_ready = true;
  return m_is_ready;
}

//_____________________________________________________________________________
Bool_t
DCTdcCalibMan::Initialize(const TString& file_name)
{
  m_file_name = file_name;
  return Initialize();
}

//_____________________________________________________________________________
DCTdcCalMap*
DCTdcCalibMan::GetMap(Int_t plane_id, Double_t wire_id) const
{
  Int_t key = MakeKey(plane_id, wire_id);
  DCTdcIterator itr = m_container.find(key);
  if(itr != m_container.end())
    return itr->second;
  else
    return nullptr;
}

//_____________________________________________________________________________
Bool_t
DCTdcCalibMan::GetTime(Int_t plane_id, Double_t wire_id,
                       Int_t tdc, Double_t& time) const
{
  DCTdcCalMap *tdc_calib = GetMap(plane_id, wire_id);
  if(tdc_calib){
    time = (tdc + (tdc_calib->p0)) * (tdc_calib->p1);
    return true;
  } else {
    hddaq::cerr << FUNC_NAME << ": No record. "
		<< " PlaneId=" << std::setw(3) << std::dec << plane_id
		<< " WireId="  << std::setw(3) << std::dec << wire_id
		<< std::endl;
    return false;
  }
}

//_____________________________________________________________________________
Bool_t
DCTdcCalibMan::GetTdc(Int_t plane_id, Double_t wire_id,
                      Double_t time, Int_t &tdc) const
{
  DCTdcCalMap *tdc_calib = GetMap(plane_id, wire_id);
  if(tdc_calib){
    tdc = Int_t((time-(tdc_calib->p0))/(tdc_calib->p1));
    return true;
  } else {
    hddaq::cerr << FUNC_NAME << ": No record. "
		<< " PlaneId=" << std::setw(3) << std::dec << plane_id
		<< " WireId="  << std::setw(3) << std::dec << wire_id
		<< std::endl;
    return false;
  }
}

//_____________________________________________________________________________
Bool_t
DCTdcCalibMan::GetParameter(Int_t plane_id, Double_t wire_id,
                            Double_t &p0, Double_t &p1) const
{
  DCTdcCalMap *tdc_calib = GetMap(plane_id, wire_id);
  if(tdc_calib){
    p0 = tdc_calib->p0;
    p1 = tdc_calib->p1;
    return true;
  } else {
    hddaq::cerr << FUNC_NAME << ": No record. "
                << " PlaneId=" << std::setw(3) << std::dec << plane_id
                << " WireId="  << std::setw(3) << std::dec << wire_id
                << std::endl;
    return false;
  }
}

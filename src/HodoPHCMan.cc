// -*- C++ -*-

#include "HodoPHCMan.hh"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <std_ostream.hh>

#include "DeleteUtility.hh"
#include "Exception.hh"
#include "FuncName.hh"
#include "MathTools.hh"

namespace
{
const Int_t SegMask  = 0x03FF;
const Int_t CidMask  = 0x00FF;
const Int_t PlidMask = 0x00FF;
const Int_t UdMask   = 0x0003;
const Int_t SegShift  =  0;
const Int_t CidShift  = 11;
const Int_t PlidShift = 19;
const Int_t UdShift   = 27;
}

//_____________________________________________________________________________
HodoPHCParam::HodoPHCParam(Int_t type, Int_t n_param,
                           const std::vector<Double_t>& parlist)
  : m_type(type),
    m_n_param(n_param),
    m_param_list(parlist)
{
}

//_____________________________________________________________________________
HodoPHCParam::~HodoPHCParam()
{
}

//_____________________________________________________________________________
inline Int_t
MakeKey(Int_t cid, Int_t pl, Int_t seg, Int_t ud)
{
  return (((cid&CidMask) << CidShift ) |
          ((pl&PlidMask) << PlidShift) |
          ((seg&SegMask) << SegShift ) |
          ((ud&UdMask)   << UdShift  ));
}

//_____________________________________________________________________________
Double_t
HodoPHCParam::DoPHC(Double_t time, Double_t de) const
{
  Double_t ctime = time;

  switch(m_type){
  case 0:
    ctime = time; break;
  case 1:
    ctime = Type1Correction(time, de); break;
  case 2:
    ctime = Type2Correction(time, de); break; // fiber
  default:
    hddaq::cerr << FUNC_NAME << ": No Correction Method. type="
		<< m_type << std::endl;
  }
  return ctime;
}

//_____________________________________________________________________________
Double_t
HodoPHCParam::DoRPHC(Double_t time, Double_t de) const
{
  Double_t ctime = time;

  switch(m_type){
  case 0:
    ctime = time; break;
  case 1:
    ctime = Type1RCorrection(time, de); break;
  default:
    hddaq::cerr << FUNC_NAME << ": No Correction Method. type="
		<< m_type << std::endl;
  }
  return ctime;
}

//_____________________________________________________________________________
Double_t
HodoPHCParam::DoSTC(Double_t stof, Double_t btof) const
{
  Double_t cstof = stof;

  switch(m_type){
  case 0:
    cstof = stof; break;
  case 1:
    cstof = Type1STCorrection(stof, btof); break;
  default:
    hddaq::cerr << FUNC_NAME << ": No Correction Method. type="
		<< m_type << std::endl;
  }
  return cstof;
}

//_____________________________________________________________________________
Double_t
HodoPHCParam::Type1Correction(Double_t time, Double_t de) const
{
  if(m_param_list.size()<3) throw Exception(FUNC_NAME+" invalid parameter");

  if(TMath::Abs(de-m_param_list[1])<MathTools::Epsilon())
    de = m_param_list[1] + MathTools::Epsilon();

  return (time
          - m_param_list[0]/TMath::Sqrt(TMath::Abs(de-m_param_list[1]))
          + m_param_list[2]);
}

//_____________________________________________________________________________
Double_t
HodoPHCParam::Type2Correction(Double_t time, Double_t w) const
{
  if(m_param_list.size()<3) throw Exception(FUNC_NAME+" invalid parameter");

  // Correction function for fiber is quadratic function
  return time - (m_param_list[0]*w*w + m_param_list[1]*w + m_param_list[2]);
}

//_____________________________________________________________________________
Double_t
HodoPHCParam::Type1RCorrection(Double_t time, Double_t de) const
{
  if(m_param_list.size()<3) throw Exception(FUNC_NAME+" invalid parameter");

  if(TMath::Abs(de-m_param_list[1]) < MathTools::Epsilon()){
    de = m_param_list[1] + MathTools::Epsilon();
  }
  return (time
          + m_param_list[0]/TMath::Sqrt(TMath::Abs(de-m_param_list[1]))
          - m_param_list[2]);
}

//_____________________________________________________________________________
Double_t
HodoPHCParam::Type1STCorrection(Double_t stof, Double_t btof) const
{
  if(m_param_list.size()<2) Exception(FUNC_NAME+" invalid parameter");

  //Correction function to eliminate the correlation between stof and btof
  return stof-(m_param_list[0]*btof + m_param_list[1]);
}

//_____________________________________________________________________________
HodoPHCMan::HodoPHCMan()
  : m_is_ready(false),
    m_file_name(),
    m_container()
{
}

//_____________________________________________________________________________
HodoPHCMan::~HodoPHCMan()
{
  ClearElements();
}

//_____________________________________________________________________________
void HodoPHCMan::ClearElements()
{
  del::ClearMap(m_container);
}

//_____________________________________________________________________________
Bool_t
HodoPHCMan::Initialize()
{
  if(m_is_ready){
    hddaq::cerr << FUNC_NAME << " already initialied" << std::endl;
    return false;
  }

  std::ifstream ifs(m_file_name);
  if(!ifs.is_open()){
    hddaq::cerr << FUNC_NAME << " file open fail : "
		<< m_file_name << std::endl;
    return false;
  }

  ClearElements();

  TString line;
  while(ifs.good() && line.ReadLine(ifs)){
    if(line.IsNull() || line[0]=='#') continue;
    std::istringstream input_line(line.Data());
    Int_t cid=-1, plid=-1, seg=-1, ud=-1, type=-1, np=-1;
    std::vector<Double_t> par;
    if(input_line >> cid >> plid >> seg >> ud >> type >> np){
      Double_t p = 0.;
      while(input_line >> p) par.push_back(p);
      Int_t key = MakeKey(cid, plid, seg, ud);
      HodoPHCParam* param = new HodoPHCParam(type, np, par);
      HodoPHCParam* pre_param = m_container[key];
      m_container[key] = param;
      if(pre_param){
	hddaq::cerr << FUNC_NAME << ": duplicated key "
		    << " following record is deleted." << std::endl
		    << " key = " << key << std::endl;
	delete pre_param;
      }
    } else {
      hddaq::cerr << FUNC_NAME << ": Invalid format" << std::endl
		  << " ===> " << line << std::endl;
    }
  }

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoPHCMan::Initialize(const TString& file_name)
{
  m_file_name = file_name;
  return Initialize();
}

//_____________________________________________________________________________
Bool_t
HodoPHCMan::DoCorrection(Int_t cid, Int_t plid, Int_t seg, Int_t ud,
                         Double_t time, Double_t de, Double_t& ctime) const
{
  ctime = time;
  HodoPHCParam* map = GetMap(cid, plid, seg, ud);
  if(!map) return false;
  ctime = map->DoPHC(time, de);
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoPHCMan::DoRCorrection(Int_t cid, Int_t plid, Int_t seg, Int_t ud,
                          Double_t time, Double_t de, Double_t& ctime) const
{
  ctime = time;
  HodoPHCParam* map = GetMap(cid, plid, seg, ud);
  if(!map) return false;
  ctime = map->DoRPHC(time, de);
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoPHCMan::DoStofCorrection(Int_t cid, Int_t plid, Int_t seg, Int_t ud,
                             Double_t stof, Double_t btof, Double_t& cstof) const
{
  cstof = stof;
  HodoPHCParam* map = GetMap(cid, plid, seg, ud);
  if(!map) return false;
  cstof = map->DoSTC(stof, btof);
  return true;
}

//_____________________________________________________________________________
HodoPHCParam*
HodoPHCMan::GetMap(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const
{
  Int_t key = MakeKey(cid, plid, seg, ud);
  PhcPIterator itr = m_container.find(key);
  if(itr != m_container.end())
    return itr->second;
  else
    return nullptr;
}

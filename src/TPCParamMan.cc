// -*- C++ -*-

#include "TPCParamMan.hh"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <std_ostream.hh>

#include "DeleteUtility.hh"
#include "FuncName.hh"
#include "TPCPadHelper.hh"

namespace
{
const Int_t RowMask = 0xff;
const Int_t LayerMask = 0xff;
const Int_t RowShift = 0;
const Int_t LayerShift = 8;

inline Int_t
MakeKey(Int_t layer, Int_t row)
{
  return (((layer & LayerMask) << LayerShift) |
          ((row & RowMask) << RowShift));
}

}

//_____________________________________________________________________________
TPCParamMan::TPCParamMan()
  : m_is_ready(false),
    m_file_name()
{
}

//_____________________________________________________________________________
TPCParamMan::~TPCParamMan()
{
  ClearACont();
  ClearTCont();
  ClearYCont();
  ClearCoboCont();
  ClearResCont();
}

//_____________________________________________________________________________
void
TPCParamMan::ClearACont()
{
  del::ClearMap(m_APContainer);
}

//_____________________________________________________________________________
void
TPCParamMan::ClearTCont()
{
  del::ClearMap(m_TPContainer);
}

//_____________________________________________________________________________
void
TPCParamMan::ClearYCont()
{
  del::ClearMap(m_YPContainer);
}

//_____________________________________________________________________________
void
TPCParamMan::ClearCoboCont()
{
  del::ClearMap(m_CoboContainer);
}

//_____________________________________________________________________________
void
TPCParamMan::ClearResCont()
{
  del::ClearMap(m_ResContainer);
}

//_____________________________________________________________________________
Bool_t
TPCParamMan::Initialize()
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

  ClearACont();
  ClearTCont();
  ClearYCont();
  ClearCoboCont();
  ClearResCont();

  Int_t line_number = 0;
  TString line;
  while(ifs.good() && line.ReadLine(ifs)){
    ++line_number;
    if(line.IsNull() || line[0]=='#') continue;
    std::istringstream input_line(line.Data());
    Int_t layer, row, aty;
    Double_t p0, p1, p2, p3, p4, p5, p6, p7, p8;
    if(input_line >> layer >> row >> aty >> p0 >> p1){
      Int_t key = MakeKey(layer, row);
      switch(aty){
      case kAdc: {
	TPCAParam* pre_param = m_APContainer[key];
	TPCAParam* param = new TPCAParam(p0, p1);
	m_APContainer[key] = param;
	if(pre_param){
	  hddaq::cerr << FUNC_NAME << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " layer = " << layer << ","
                      << " row = " << row
                      << " aty = " << aty
                      << std::endl;
	  delete pre_param;
	}
        break;
      }
      case kTdc: {
	TPCTParam *pre_param = m_TPContainer[key];
	TPCTParam *param = new TPCTParam(p0, p1);
	m_TPContainer[key] = param;
	if(pre_param){
	  hddaq::cerr << FUNC_NAME << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " layer = " << layer << ","
                      << " row = " << row
                      << " aty = " << aty
                      << std::endl;
	  delete pre_param;
	}
        break;
      }
      case kY: {
        TPCYParam *pre_param = m_YPContainer[key];
        TPCYParam *param = new TPCYParam(p0, p1);
        m_YPContainer[key] = param;
        if(pre_param){
	  hddaq::cerr << FUNC_NAME << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " layer = " << layer << ","
                      << " row = " << row
                      << " aty = " << aty
                      << std::endl;
          delete pre_param;
        }
        break;
      }
      case kCobo: {
        if(input_line >> p2){
          //layer is Cobo number (not real layer)
          Int_t cobo = layer;
          TPCCoboParam *pre_param = m_CoboContainer[cobo];
          // TPCCoboParam *param = new TPCCoboParam(p0);
          std::vector<Double_t> params{ p0, p1, p2 };
          TPCCoboParam *param = new TPCCoboParam(params);
          m_CoboContainer[cobo] = param;
          if(pre_param){
            hddaq::cerr << FUNC_NAME << ": duplicated key "
                        << " following record is deleted." << std::endl
                        << " layer = " << layer << ","
                        << " row = " << row
                        << " aty = " << aty
                        << std::endl;
            delete pre_param;
          }
        }
        break;
      }
      case kRes: {
	if(input_line >> p2 >> p3 >> p4 >> p5 >> p6 >> p7 >> p8){
	  TPCResParam *pre_param = m_ResContainer[key];
	  std::vector<Double_t> params{ p0, p1, p2, p3, p4, p5, p6, p7, p8 };
	  TPCResParam *param = new TPCResParam(params);
	  m_ResContainer[key] = param;

	  Int_t HS = layer; Int_t InOut = row;
	  if(HS==0&&InOut==0) m_Res_HSOFF_Inner = params;
	  else if(HS==0&&InOut==1) m_Res_HSOFF_Outer = params;
	  else if(HS==1&&InOut==0) m_Res_HSON_Inner = params;
	  else if(HS==1&&InOut==1) m_Res_HSON_Outer = params;
	  else hddaq::cerr << FUNC_NAME << ": Invalid TPC resolution parameter"
			   << std::endl << " p0, 1: HS On 0 : HS Off"
			   << std::endl << " p1, 1: Outer layers 0 : inner layers "
			   << std::endl;

	  if(pre_param){
	    hddaq::cerr << FUNC_NAME << ": duplicated key "
			<< " following record is deleted." << std::endl
			<< " layer = " << layer << ","
			<< " row = " << row
			<< " aty = " << aty
			<< std::endl;
	    delete pre_param;
	  }
	  break;
	}
      }
      default:
        hddaq::cerr << FUNC_NAME << ": Invalid Input" << std::endl
                    << " ===> L" << line_number << " " << line << std::endl;
        break;
      } // switch
    } else {
      hddaq::cerr << FUNC_NAME << ": Invalid Input" << std::endl
		  << " ===> L" << line_number << " " << line << std::endl;
    }
  } // while

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCParamMan::Initialize(const TString& file_name)
{
  m_file_name = file_name;
  return Initialize();
}

//_____________________________________________________________________________
Bool_t
TPCParamMan::GetCDe(Int_t layer, Int_t row, Double_t adc, Double_t& de) const
{
  TPCAParam* map = GetAmap(layer, row);
  if(map){
    de = map->CDeltaE(adc);
    return true;
  } else {
    hddaq::cerr << FUNC_NAME << ": No record for"
		<< " LayerId=" << std::setw(3) << std::dec << layer
		<< " RowId="   << std::setw(3) << std::dec << row
		<< std::endl;
    return false;
  }
}

//_____________________________________________________________________________
Bool_t
TPCParamMan::GetCTime(Int_t layer, Int_t row, Double_t tdc,
                      Double_t& time) const
{
  TPCTParam* map = GetTmap(layer, row);
  if(map){
    time = map->Time(tdc);
    return true;
  } else {
    hddaq::cerr << FUNC_NAME << ": No record for"
		<< " LayerId=" << std::setw(3) << std::dec << layer
		<< " RowId="   << std::setw(3) << std::dec << row
		<< std::endl;
    return false;
  }
}

//_____________________________________________________________________________
Bool_t
TPCParamMan::GetCClock(Int_t layer, Int_t row, Double_t clock,
                       Double_t& cclk) const
{
  Int_t cobo = tpc::GetCoBoId(layer, row);
  TPCCoboParam* map = GetCobomap(cobo);
  if(map){
    cclk = map->PhaseShift(clock);
    return true;
  } else {
    hddaq::cerr << FUNC_NAME << ": No record for"
		<< " CoBoId=" << std::setw(3) << std::dec << cobo
		<< " LayerId=" << std::setw(3) << std::dec << layer
		<< " RowId="   << std::setw(3) << std::dec << row
		<< std::endl;
    return false;
  }
}

//_____________________________________________________________________________
Bool_t
TPCParamMan::GetY(Int_t layer, Int_t row, Double_t time, Double_t& y) const
{
  TPCYParam* map = GetYmap(layer, row);
  if(map){
    y = map->Y(time);
    return true;
  } else {
    hddaq::cerr << FUNC_NAME << ": No record for"
		<< " LayerId=" << std::setw(3) << std::dec << layer
		<< " RowId="   << std::setw(3) << std::dec << row
		<< std::endl;
    return false;
  }
}

//_____________________________________________________________________________
TPCAParam*
TPCParamMan::GetAmap(Int_t layer, Int_t row) const
{
  Int_t key = MakeKey(layer, row);
  AIterator itr = m_APContainer.find(key);
  if(itr != m_APContainer.end())
    return itr->second;
  else
    return nullptr;
}

//_____________________________________________________________________________
TPCTParam*
TPCParamMan::GetTmap(Int_t layer, Int_t row) const
{
  Int_t key = MakeKey(layer, row);
  TIterator itr = m_TPContainer.find(key);
  if(itr != m_TPContainer.end())
    return itr->second;
  else
    return nullptr;
}

//_____________________________________________________________________________
TPCYParam*
TPCParamMan::GetYmap(Int_t layer, Int_t row) const
{
  Int_t key = MakeKey(layer, row);
  YIterator itr = m_YPContainer.find(key);
  if(itr != m_YPContainer.end())
    return itr->second;
  else
    return nullptr;
}

//_____________________________________________________________________________
TPCCoboParam*
TPCParamMan::GetCobomap(Int_t cobo) const
{
  CoboIterator itr = m_CoboContainer.find(cobo);
  if(itr != m_CoboContainer.end())
    return itr->second;
  else
    return nullptr;
}

//_____________________________________________________________________________
TPCResParam*
TPCParamMan::GetResmap(Int_t B, Int_t InnerOrOuter) const
{
  Int_t key = MakeKey(B, InnerOrOuter);
  ResIterator itr = m_ResContainer.find(key);
  if(itr != m_ResContainer.end())
    return itr->second;
  else
    return nullptr;
}

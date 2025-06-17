// -*- C++ -*-

#include "Hodo1Hit.hh"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "FuncName.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"

namespace
{
const auto& gHodo = HodoParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();
}

//_____________________________________________________________________________
HodoHit::HodoHit()
{
}

//_____________________________________________________________________________
HodoHit::~HodoHit()
{
}

//_____________________________________________________________________________
Hodo1Hit::Hodo1Hit(HodoRawHit *rhit, Int_t index)
  : HodoHit(), m_raw(rhit), m_is_calculated(false), m_index(index)
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
Hodo1Hit::~Hodo1Hit()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
Bool_t
//Hodo1Hit::Calculate()
Hodo1Hit::Calculate(Bool_t tdc_flag)
{
  if(m_is_calculated){
    hddaq::cout << FUNC_NAME << " already calculated" << std::endl;
    return false;
  }

  if(tdc_flag && m_raw->GetNumOfTdcHits()!=1)
    return false;

  if(!gHodo.IsReady()){
    hddaq::cout << FUNC_NAME << " HodoParamMan must be initialized" << std::endl;
    return false;
  }

  if(!gPHC.IsReady()){
    hddaq::cout << FUNC_NAME << " HodoPHCMan must be initialized" << std::endl;
    return false;
  }

  Int_t cid  = m_raw->DetectorId();
  Int_t plid = m_raw->PlaneId();
  Int_t seg  = m_raw->SegmentId();

  // hit information
  m_multi_hit_l = m_raw->SizeTdc1()  > m_raw->SizeTdc2()?  m_raw->SizeTdc1()  : m_raw->SizeTdc2();
  m_multi_hit_t = m_raw->SizeTdcT1() > m_raw->SizeTdcT2()? m_raw->SizeTdcT1() : m_raw->SizeTdcT2();
  Int_t UorD=0;
  Int_t adc = m_raw->GetAdc1(0);
  if(0 > m_raw->GetTdc1(0)){
    UorD=1;
    adc = m_raw->GetAdc2(0);
  }

  Double_t dE = 0.;
  if(adc>=0){
    if(!gHodo.GetDe(cid, plid, seg, UorD, adc, dE)){
      hddaq::cerr << FUNC_NAME
		  << " something is wrong at GetDe("
		  << cid  << ", " << plid << ", " << seg << ", "
		  << UorD << ", " << adc  << ", " << dE  << ")" << std::endl;
      return false;
    }
  }

  m_a.push_back(dE);

  Int_t mhit = m_multi_hit_l;
  for(Int_t m=0; m<mhit; ++m){
    Int_t    tdc  = -999;
    Double_t time = -999.;
    if(0 == UorD){
      tdc = m_raw->GetTdc1(m);
    }else{
      tdc = m_raw->GetTdc2(m);
    }

    if(tdc<0) continue;

    if(!gHodo.GetTime(cid, plid, seg, UorD, tdc, time)){
      hddaq::cerr << FUNC_NAME
		  << " something is wrong at GetTime("
		  << cid  << ", " << plid << ", " << seg  << ", "
		  << UorD << ", " << tdc  << ", " << time << ")" << std::endl;
      return false;
    }
    m_t.push_back(time);

    Double_t ctime = time;
    gPHC.DoCorrection(cid, plid, seg, UorD, time, dE, ctime);
    m_ct.push_back(ctime);

    m_flag_join.push_back(false);
  }

  m_is_calculated = true;
  return true;
}

// ____________________________________________________________
Int_t
Hodo1Hit::GetNumOfHit(Int_t sel) const
{
  return sel == 0 ? m_multi_hit_l : m_multi_hit_t;
}

// ____________________________________________________________
Bool_t
Hodo1Hit::JoinedAllMhit()
{
  Bool_t ret = true;
  for(Int_t i=0, n=m_flag_join.size(); i<n; ++i){
    ret = ret & m_flag_join[i];
  }
  return ret;
}

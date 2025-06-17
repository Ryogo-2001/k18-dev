// -*- C++ -*-

#include "Hodo2Hit.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "DebugCounter.hh"
#include "FuncName.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"

#include <std_ostream.hh>

namespace
{
const auto& gHodo = HodoParamMan::GetInstance();
const auto& gPHC = HodoPHCMan::GetInstance();
}

//_____________________________________________________________________________
Hodo2Hit::Hodo2Hit(HodoRawHit *rhit, Double_t max_time_diff)
  : HodoHit(),
    m_raw(rhit),
    m_is_calculated(false),
    m_is_tof(false),
    m_max_time_diff(max_time_diff)
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
Hodo2Hit::~Hodo2Hit()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
Bool_t
Hodo2Hit::Calculate()
{
  if(m_is_calculated){
    hddaq::cout << FUNC_NAME << " already calculated" << std::endl;
    return false;
  }

  if(m_raw->GetNumOfTdcHits() != 2){
    return false;
  }

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
  Int_t adc1 = m_raw->GetAdc1(), adc2=m_raw->GetAdc2();

  Double_t offset_vtof = 0;
  if(m_is_tof){
    gHodo.GetTime(cid, plid, seg, 2, 0, offset_vtof);
    offset_vtof = offset_vtof;
  }// offset vtof

  Double_t dE1 = 0., dE2 = 0.;
  if(adc1>=0){
    if(!gHodo.GetDe(cid, plid, seg, 0, adc1, dE1)){
      hddaq::cerr << FUNC_NAME
		  << " something is wrong at GetDe("
		  << cid  << ", " << plid << ", " << seg  << ", "
		  << "0, " << adc1 << ", " << dE1 << ")" << std::endl;
      return false;
    }
  }
  if(adc2>=0){
    if(!gHodo.GetDe(cid, plid, seg, 1, adc2, dE2)){
      hddaq::cerr << FUNC_NAME
		  << " something is wrong at GetDe("
		  << cid  << ", " << plid << ", " << seg  << ", "
		  << "0, " << adc2 << ", " << dE2 << ")" << std::endl;
      return false;
    }
  }

  m_a1.push_back(dE1);
  m_a2.push_back(dE2);

  // Multi-hit analysis ________________________________________________________
  Int_t n_mhit1 = m_raw->SizeTdc1();
  Int_t n_mhit2 = m_raw->SizeTdc2();

  // local container
  std::vector<Double_t> time1;
  std::vector<Double_t> time2;
  std::vector<Double_t> ctime1;
  std::vector<Double_t> ctime2;

  // Tdc1
  for(Int_t i = 0; i<n_mhit1; ++i){
    Int_t tdc = m_raw->GetTdc1(i);
    Double_t time = 0.;
    if(!gHodo.GetTime(cid, plid, seg, 0, tdc, time)){
      hddaq::cerr << "#E " << FUNC_NAME
		  << " something is wrong at GetTime("
		  << cid  << ", " << plid << ", " << seg  << ", "
		  << "U, " << tdc << ", " << "time"
		  << ")" << std::endl;
      return false;
    }
    time1.push_back(time);
    Double_t ctime = -999.;
    gPHC.DoCorrection(cid, plid, seg, 0, time, dE1, ctime);
    ctime1.push_back(ctime);
  }// for(i)

  // Tdc2
  for(Int_t i = 0; i<n_mhit2; ++i){
    Int_t tdc = m_raw->GetTdc2(i);
    Double_t time = 0.;
    if(!gHodo.GetTime(cid, plid, seg, 1, tdc, time)){
      hddaq::cerr << "#E " << FUNC_NAME
		  << " something is wrong at GetTime("
		  << cid  << ", " << plid << ", " << seg  << ", "
		  << "D, " << tdc << ", " << "time"
		  << ")" << std::endl;
      return false;
    }
    time2.push_back(time);
    Double_t ctime = -999.;
    gPHC.DoCorrection(cid, plid, seg, 1, time, dE2, ctime);
    ctime2.push_back(ctime);
  }

  std::sort(time1.begin(),   time1.end(), std::greater<Double_t>());
  std::sort(time2.begin(),   time2.end(), std::greater<Double_t>());
  std::sort(ctime1.begin(), ctime1.end(), std::greater<Double_t>());
  std::sort(ctime2.begin(), ctime2.end(), std::greater<Double_t>());

#if 0
  // Make coincidence
  for(Int_t i = 0, i_d = 0, last_d = 0; i<n_mhit1; ++i){
    while(i_d < n_mhit2){
      if(abs(time1[i]-time2[i_d]) < m_max_time_diff){
	// Coincidence
	data_pair a_pair;
	a_pair.time1   = time1[i]     + offset_vtof;
	a_pair.time2   = time2[i_d]   + offset_vtof;
	a_pair.ctime1  = ctime1[i]    + offset_vtof;
	a_pair.ctime2  = ctime2[i_d]  + offset_vtof;
	m_pair_cont.push_back(a_pair);
	last_d         = i_d+1;
	++i_d;

	m_flag_join.push_back(false);
	break;
      }
      ++i_d;
    }// go to next D hit
    if(last_d == n_mhit2) break; // no more D hit
    if(i_d == n_mhit2-1) i_d = last_d;
  }
#endif

#if 1
  // Make coincidence
  for(Int_t i = 0; i<n_mhit1; ++i){
    for (Int_t i_d = 0; i_d<n_mhit2; ++i_d) {
      if(abs(time1[i]-time2[i_d]) < m_max_time_diff){
        // Coincidence
        data_pair a_pair;
        a_pair.time1   = time1[i]     + offset_vtof;
        a_pair.time2   = time2[i_d]   + offset_vtof;
        a_pair.ctime1  = ctime1[i]    + offset_vtof;
        a_pair.ctime2  = ctime2[i_d]  + offset_vtof;
        m_pair_cont.push_back(a_pair);
        m_flag_join.push_back(false);
        break;
      }
    }
  }
#endif

  if(m_pair_cont.size() == 0) return false;

  m_is_calculated = true;

  return true;
}

// ____________________________________________________________
Bool_t
Hodo2Hit::JoinedAllMhit()
{
  Bool_t ret = true;
  for(Int_t i=0, n=m_flag_join.size(); i<n; ++i){
    ret = ret & m_flag_join[i];
  }
  return ret;
}

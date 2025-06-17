// -*- C++ -*-

#include "HodoRawHit.hh"

#include <iterator>

#include "DebugCounter.hh"
#include "FuncName.hh"

//_____________________________________________________________________________
HodoRawHit::HodoRawHit(Int_t detector_id, Int_t plane_id, Int_t segment_id)
  : m_detector_id(detector_id),
    m_plane_id(plane_id),
    m_segment_id(segment_id),
    m_adc1(1, -1),
    m_adc2(1, -1),
    m_tdc1(1, -1),
    m_tdc2(1, -1),
    m_tdc_t1(1, -1),
    m_tdc_t2(1, -1),
    m_oftdc(false),
    m_nhtdc(0)
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
HodoRawHit::~HodoRawHit()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
HodoRawHit::SetAdc1(Int_t adc)
{
  if(-1 == m_adc1.at(0))
    m_adc1.at(0) = adc;
  else
    m_adc1.push_back(adc);
}

//_____________________________________________________________________________
void
HodoRawHit::SetAdc2(Int_t adc)
{
  if(-1 == m_adc2.at(0))
    m_adc2.at(0) = adc;
  else
    m_adc2.push_back(adc);
}

//_____________________________________________________________________________
void
HodoRawHit::SetTdc1(Int_t tdc)
{
  if(-1 == m_tdc1.at(0)){
    m_tdc1.at(0) = tdc;
    ++m_nhtdc;
  }else{
    m_tdc1.push_back(tdc);
  }
}

//_____________________________________________________________________________
void
HodoRawHit::SetTdc2(Int_t tdc)
{
  if(-1 == m_tdc2.at(0)){
    m_tdc2.at(0) = tdc;
    ++m_nhtdc;
  }else{
    m_tdc2.push_back(tdc);
  }
}

//_____________________________________________________________________________
void
HodoRawHit::SetTdcT1(Int_t tdc)
{
  if(-1 == m_tdc_t1.at(0)){
    m_tdc_t1.at(0) = tdc;
  }else{
    m_tdc_t1.push_back(tdc);
  }
}

//_____________________________________________________________________________
void
HodoRawHit::SetTdcT2(Int_t tdc)
{
  if(-1 == m_tdc_t2.at(0)){
    m_tdc_t2.at(0) = tdc;
  }else{
    m_tdc_t2.push_back(tdc);
  }
}

//_____________________________________________________________________________
Int_t
HodoRawHit::SizeAdc1() const
{
  if(-1 == m_adc1.at(0))
    return 0;
  else
    return m_adc1.size();
}

//_____________________________________________________________________________
Int_t
HodoRawHit::SizeAdc2() const
{
  if(-1 == m_adc2.at(0))
    return 0;
  else
    return m_adc2.size();
}

//_____________________________________________________________________________
Int_t
HodoRawHit::SizeTdc1() const
{
  if(-1 == m_tdc1.at(0))
    return 0;
  else
    return m_tdc1.size();
}

//_____________________________________________________________________________
Int_t
HodoRawHit::SizeTdc2() const
{
  if(-1 == m_tdc2.at(0))
    return 0;
  else
    return m_tdc2.size();
}

//_____________________________________________________________________________
Int_t
HodoRawHit::SizeTdcT1() const
{
  if(-1 == m_tdc_t1.at(0))
    return 0;
  else
    return m_tdc_t1.size();
}

//_____________________________________________________________________________
Int_t
HodoRawHit::SizeTdcT2() const
{
  if(-1 == m_tdc_t2.at(0))
    return 0;
  else
    return m_tdc_t2.size();
}

//_____________________________________________________________________________
void
HodoRawHit::Clear()
{
  m_nhtdc = 0;
  m_adc1.clear();
  m_adc2.clear();
  m_tdc1.clear();
  m_tdc2.clear();
  m_tdc_t1.clear();
  m_tdc_t2.clear();
  m_adc1.push_back(-1);
  m_adc2.push_back(-1);
  m_tdc1.push_back(-1);
  m_tdc2.push_back(-1);
  m_tdc_t1.push_back(-1);
  m_tdc_t2.push_back(-1);
}

//_____________________________________________________________________________
void
HodoRawHit::Print(const TString& arg)
{
  hddaq::cerr << FUNC_NAME << " " << arg << std::endl
	      << "detector_id = " << m_detector_id << std::endl
	      << "plane_id    = " << m_plane_id    << std::endl
	      << "segment_id  = " << m_segment_id   << std::endl;
  std::vector<Int_t>::const_iterator itr, end;
  hddaq::cout << "adc1        = " << m_adc1.size() << " ";
  std::copy(m_adc1.begin(), m_adc1.end(),
	     std::ostream_iterator<Int_t>(hddaq::cout," "));
  hddaq::cout << std::endl
	      << "adc2        = " << m_adc2.size() << " ";
  std::copy(m_adc2.begin(), m_adc2.end(),
	     std::ostream_iterator<Int_t>(hddaq::cout," "));
  hddaq::cout << std::endl
	      << "tdc1        = " << m_tdc1.size() << " ";
  std::copy(m_tdc1.begin(), m_tdc1.end(),
	     std::ostream_iterator<Int_t>(hddaq::cout," "));
  hddaq::cout << std::endl
	      << "tdc2        = " << m_tdc2.size() << " ";
  std::copy(m_tdc2.begin(), m_tdc2.end(),
	     std::ostream_iterator<Int_t>(hddaq::cout," "));
  hddaq::cout << std::endl
	      << "tdc_t1      = " << m_tdc_t1.size() << " ";
  std::copy(m_tdc_t1.begin(), m_tdc_t1.end(),
	     std::ostream_iterator<Int_t>(hddaq::cout," "));
  hddaq::cout << std::endl
	      << "tdc_t2      = " << m_tdc_t2.size() << " ";
  std::copy(m_tdc2.begin(), m_tdc2.end(),
	     std::ostream_iterator<Int_t>(hddaq::cout," "));
  hddaq::cout << std::endl;
}

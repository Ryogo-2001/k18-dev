// -*- C++ -*-

#include "DCRawHit.hh"

#include <iostream>
#include <iterator>

#include "std_ostream.hh"

#include "DebugCounter.hh"
#include "FuncName.hh"

//_____________________________________________________________________________
DCRawHit::DCRawHit(Int_t plane_id, Int_t wire_id)
  : m_plane_id(plane_id),
    m_wire_id(wire_id),
    m_oftdc(false)
{
  m_tdc.clear();
  m_trailing.clear();
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
DCRawHit::~DCRawHit()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
DCRawHit::Print(const TString& arg) const
{
  hddaq::cerr << FUNC_NAME << " " << arg << std::endl
	      << "plane_id = " << m_plane_id    << std::endl
	      << "wire_id  = " << m_wire_id     << std::endl;
  std::vector<Int_t>::const_iterator itr, end;
  hddaq::cout << "tdc      = " << m_tdc.size() << " ";
  std::copy(m_tdc.begin(), m_tdc.end(),
	     std::ostream_iterator<Int_t>(hddaq::cout," "));
  hddaq::cout << std::endl
	      << "trailing = " << m_trailing.size() << " ";
  std::copy(m_trailing.begin(), m_trailing.end(),
	     std::ostream_iterator<Int_t>(hddaq::cout," "));
  hddaq::cout << std::endl;
}

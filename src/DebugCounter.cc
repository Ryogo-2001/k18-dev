// -*- C++ -*-

#include "DebugCounter.hh"

#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>

#include <std_ostream.hh>

#include "FuncName.hh"

namespace debug
{

//_____________________________________________________________________________
ObjectCounter::ObjectCounter()
  : m_map()
{
}

//_____________________________________________________________________________
ObjectCounter::~ObjectCounter()
{
}

//_____________________________________________________________________________
void
ObjectCounter::check(const TString& arg) const
{
#ifdef MemoryLeak
  Bool_t has_leak = false;
  ObjectIter itr, end=m_map.end();
  for(itr=m_map.begin(); itr!=end; ++itr){
    if(itr->second!=0) has_leak = true;
  }
  if(has_leak)
    print(arg+" "+FUNC_NAME);
#endif
}

//_____________________________________________________________________________
void
ObjectCounter::print(const TString& arg) const
{
  hddaq::cout << "#DCounter " << FUNC_NAME << " " << arg << std::endl;
  ObjectIter itr, end=m_map.end();
  for(itr=m_map.begin(); itr!=end; ++itr){
    hddaq::cout << std::setw(20) << std::left
		<< itr->first << " : " << itr->second << std::endl;
  }
}

}

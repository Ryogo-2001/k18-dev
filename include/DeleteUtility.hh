// -*- C++ -*-

#ifndef DELETE_UTILITY_HH
#define DELETE_UTILITY_HH

#include <algorithm>
#include <iterator>
#include <signal.h>
#include <vector>

namespace del
{
//_____________________________________________________________________________
struct DeleteFunction
{
  template <typename T>
  void operator ()(T*& p) const
    {
      if(p){
        delete p;
        p = nullptr;
      }
    }
};

//_____________________________________________________________________________
template <typename T>
inline void DeleteObject(T*& p)
{
  if(p){
    delete p;
    p = nullptr;
  }
}

//_____________________________________________________________________________
template <typename Container>
inline void ClearContainer(Container& c, bool ownership=true)
{
  if(ownership)
    std::for_each(c.begin(), c.end(), DeleteFunction());
  c.clear();
}

//_____________________________________________________________________________
template <typename Container>
inline void DeleteObject(Container& c, bool ownership=true)
{
  ClearContainer(c,ownership);
}

//_____________________________________________________________________________
template <typename Container>
inline void ClearContainerAll(Container& c, bool ownership=true)
{
  typename Container::iterator itr, itr_end=c.end();
  for(itr=c.begin(); itr!=itr_end; ++itr)
    ClearContainer(*itr, ownership);
}

//_____________________________________________________________________________
template <typename Map>
inline void ClearMap(Map& m, bool ownership=true)
{
  typename Map::const_iterator itr, end=m.end();
  if(ownership){
    for(itr=m.begin(); itr!=end; ++itr){
      delete itr->second;
    }
  }
  m.clear();
}

}

#endif

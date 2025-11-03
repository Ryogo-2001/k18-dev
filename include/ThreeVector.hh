// -*- C++ -*-

#ifndef THREE_VECTOR_HH
#define THREE_VECTOR_HH

#include <ostream>

#include <TVector3.h>
typedef TVector3 ThreeVector;

//_____________________________________________________________________________
inline std::ostream&
operator <<(std::ostream& ost, const TVector3& v)
{
  ost << "(" << v.x() << ", " << v.y()
      << ", " << v.z() << ")";
  return ost;
}

#endif

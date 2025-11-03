// -*- C++ -*-

#ifndef FIBER_CLUSTER_HH
#define FIBER_CLUSTER_HH

#include <vector>

#include <TString.h>

#include "FiberHit.hh"
#include "HodoCluster.hh"

//_____________________________________________________________________________
class FiberCluster : public HodoCluster
{
public:
  static const TString& ClassName();
  FiberCluster(const HodoHC& cont,
               const index_t& index);
  virtual ~FiberCluster();

private:
  FiberCluster(const FiberCluster&);
  FiberCluster& operator =(const FiberCluster&);

protected:
  void Calculate();
};

//_____________________________________________________________________________
inline const TString&
FiberCluster::ClassName()
{
  static TString s_name("FiberCluster");
  return s_name;
}

#endif

// -*- C++ -*-

#ifndef DC_PAIR_HIT_CLUSTER_HH
#define DC_PAIR_HIT_CLUSTER_HH

#include <TString.h>

class DCLTrackHit;

//_____________________________________________________________________________
class DCPairHitCluster
{
public:
  static const TString& ClassName();
  DCPairHitCluster(DCLTrackHit* hitA, DCLTrackHit* hitB=nullptr);
  ~DCPairHitCluster();

private:
  DCPairHitCluster(const DCPairHitCluster&);
  DCPairHitCluster& operator =(const DCPairHitCluster&);

private:
  DCLTrackHit *m_hitA;
  DCLTrackHit *m_hitB;
  Int_t        m_nhits;

public:
  Int_t        NumberOfHits() const { return m_nhits; }
  DCLTrackHit* GetHit(Int_t i=0) const;
  bool         IsHoneycomb() const;
  void         SetHoneycomb(bool flag);
};

//_____________________________________________________________________________
inline const TString&
DCPairHitCluster::ClassName()
{
  static TString s_name("DCPairHitCluster");
  return s_name;
}

#endif

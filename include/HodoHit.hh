// -*- C++ -*-

#ifndef HODO_HIT_HH
#define HODO_HIT_HH

#include <TString.h>

//_____________________________________________________________________________
class HodoHit
{
public:
  HodoHit();
  virtual ~HodoHit() = 0;

public:
  virtual Double_t DeltaE(Int_t n=0) const = 0;
  virtual Double_t GetTUp(Int_t n=0) const = 0;
  virtual Double_t GetTDown(Int_t n=0) const = 0;
  virtual Double_t MeanTime(Int_t n=0) const = 0;
  virtual Double_t CMeanTime(Int_t n=0) const = 0;
  virtual Int_t    DetectorId() const = 0;
  virtual Int_t    PlaneId() const = 0;
  virtual Int_t    SegmentId() const = 0;
  virtual Bool_t   ReCalc(Bool_t applyRecursively=false) = 0;
};

#endif

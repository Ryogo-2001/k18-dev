// -*- C++ -*-

#ifndef VEVENT_HH
#define VEVENT_HH

#include <TString.h>

//_____________________________________________________________________________
class VEvent
{
public:
  VEvent();
  virtual ~VEvent() = 0;
  virtual const TString& ClassName() = 0;
  virtual Bool_t ProcessingBegin() = 0;
  virtual Bool_t ProcessingEnd() = 0;
  virtual Bool_t ProcessingNormal() = 0;
};

#endif

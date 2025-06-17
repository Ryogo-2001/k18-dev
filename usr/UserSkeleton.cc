// -*- C++ -*-

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#include <UnpackerManager.hh>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "RootHelper.hh"
#include "HodoRawHit.hh"
#include "KuramaLib.hh"
#include "RawData.hh"

namespace
{
using namespace root;
auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
}

//_____________________________________________________________________________
class UserSkeleton : public VEvent
{
private:

public:
  UserSkeleton();
  ~UserSkeleton();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserSkeleton::ClassName()
{
  static TString s_name("UserSkeleton");
  return s_name;
}

//_____________________________________________________________________________
UserSkeleton::UserSkeleton()
  : VEvent()
{
}

//_____________________________________________________________________________
UserSkeleton::~UserSkeleton()
{
}

//_____________________________________________________________________________
struct Event
{
  Int_t runnum;
  Int_t evnum;
  Int_t spill;
  void clear()
    {
      runnum = -1;
      evnum = -1;
      spill = -1;
    }
};

//_____________________________________________________________________________
namespace root
{
Event  event;
TH1   *h[MaxHist];
TTree *tree;
}

//_____________________________________________________________________________
Bool_t
UserSkeleton::ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserSkeleton::ProcessingNormal()
{
  event.runnum = gUnpacker.get_run_number();
  event.evnum  = gUnpacker.get_event_number();
  // rawData = new RawData;
  // rawData->DecodeHits();

  // for(Int_t i=0; i<100; ++i){
  //   HF1(i, (double)i);
  // }

  return true;
}

//_____________________________________________________________________________
Bool_t
UserSkeleton::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserSkeleton;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  for(Int_t i=0; i<100; ++i){
    HB1(i, Form("hist %d", i), 100, 0., 100.);
  }

  HBTree("skeleton","tree of Skeleton");
  tree->Branch("runnum", &event.runnum, "runnum/I");
  tree->Branch("evnum",  &event.evnum,  "evnum/I");
  tree->Branch("spill",  &event.spill,  "spill/I");

  HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO") &&
     InitializeParameter<HodoParamMan>("HDPRM"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

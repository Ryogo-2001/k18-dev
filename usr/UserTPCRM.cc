// -*- C++ -*-

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RootHelper.hh"
#include "Unpacker.hh"
#include "UnpackerManager.hh"

namespace
{
using namespace root;
using hddaq::unpacker::GUnpacker;
const auto& gUnpacker = GUnpacker::get_instance();
}

//_____________________________________________________________________________
class UserTPCRM : public VEvent
{
public:
  UserTPCRM();
  ~UserTPCRM();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserTPCRM::ClassName()
{
  static TString s_name("UserTPCRM");
  return s_name;
}

//_____________________________________________________________________________
UserTPCRM::UserTPCRM()
  : VEvent()
{
}

//_____________________________________________________________________________
UserTPCRM::~UserTPCRM()
{
}

//_____________________________________________________________________________
struct Event
{
  Int_t runnum;
  Int_t evnum;
  Int_t e03_runnum;
  Int_t e03_event;
  Int_t e03_spill;

  void clear()
  {
    runnum = 0;
    evnum = 0;
    e03_runnum = 0;
    e03_event = 0;
    e03_spill = 0;
  }
  void print()
  {
    std::cout << TString('=', 80) << std::endl
              << "runnum     = " << runnum << std::endl
              << "evnum      = " << evnum << std::endl
              << "e03_runnum = " << e03_runnum << std::endl
              << "e03_event  = " << e03_event << std::endl
              << "e03_spill  = " << e03_spill << std::endl
              << std::endl;
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
bool
UserTPCRM::ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
bool
UserTPCRM::ProcessingNormal()
{
  const Int_t run_number   = gUnpacker.get_root()->get_run_number();
  const Int_t event_number = gUnpacker.get_event_number();
  event.runnum = run_number;
  event.evnum  = event_number;

  // VME RM
  {
    static const Int_t device_id = gUnpacker.get_device_id("VME-RM");
    static const Int_t plane_id = gUnpacker.get_plane_id("VME-RM", "vme08");
    // static const Int_t event_id = gUnpacker.get_data_id("VME-RM", "event");
    // static const Int_t spill_id = gUnpacker.get_data_id("VME-RM", "spill");
    // static const Int_t serial_id = gUnpacker.get_data_id("VME-RM", "serial");
    // static const Int_t time_id = gUnpacker.get_data_id("VME-RM", "time");
    enum EChannelType { kTag, kLock, kNChannelType };
    enum EDataType { kEvent, kSpill, kSerial, kTime, kNDataType };
    std::vector<std::vector<UInt_t>> data(kNChannelType);
    for (Int_t ch=0; ch<kNChannelType; ++ch) {
      data[ch].resize(kNDataType);
      for (Int_t type=0; type<kNDataType; ++type) {
        if (gUnpacker.get_entries(device_id, plane_id, 0, ch, type) > 0) {
          data[ch][type] = gUnpacker.get(device_id, plane_id, 0, ch, type);
        }
      }
    }
    static const UInt_t lock_shift = 31;
    if ((data[kLock][kEvent] >> lock_shift) == 1)
      event.e03_event = data[kTag][kEvent];
    if ((data[kLock][kSpill] >> lock_shift) == 1)
      event.e03_spill = data[kTag][kSpill];
    event.e03_runnum = data[kTag][kTime];
  }

  // event.print();

  for( Int_t i=0; i<100; ++i ){
    HF1( i, (double)i );
  }

  return true;
}

//_____________________________________________________________________________
bool
UserTPCRM::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserTPCRM;
}

//_____________________________________________________________________________
bool
ConfMan::InitializeHistograms()
{
  HBTree("tpcrm", "tree of TPCRM");
  tree->Branch("runnum", &event.runnum);
  tree->Branch("evnum",  &event.evnum);
  tree->Branch("e03_runnum", &event.e03_runnum);
  tree->Branch("e03_event", &event.e03_event);
  tree->Branch("e03_spill", &event.e03_spill);

  HPrint();
  return true;
}

//_____________________________________________________________________________
bool
ConfMan::InitializeParameterFiles()
{
  return true;
}

//_____________________________________________________________________________
bool
ConfMan::FinalizeProcess()
{
  return true;
}

// -*- C++ -*-

#include "VEvent.hh"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

#include <DAQNode.hh>
#include <Unpacker.hh>
#include <UnpackerManager.hh>

#include "ConfMan.hh"
#include "RootHelper.hh"

namespace
{
using namespace root;
using hddaq::unpacker::DAQNode;
using hddaq::unpacker::GUnpacker;
auto& gUnpacker = GUnpacker::get_instance();
}

//_____________________________________________________________________________
class UserBookMark : public VEvent
{
public:
  UserBookMark();
  ~UserBookMark();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserBookMark::ClassName()
{
  static TString s_name("UserBookMark");
  return s_name;
}

//_____________________________________________________________________________
UserBookMark::UserBookMark()
  : VEvent()
{
}

//_____________________________________________________________________________
UserBookMark::~UserBookMark()
{
}

//_____________________________________________________________________________
struct Event
{
  Int_t runnum;
  Int_t evnum;
  Int_t dsize; // [word=4Bytes]
  void clear()
    {
      runnum = -1;
      evnum = -1;
      dsize = -1;
    }
};

//_____________________________________________________________________________
namespace root
{
Event event;
TH1   *h[MaxHist];
TTree *tree;
}

//_____________________________________________________________________________
Bool_t
UserBookMark::ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserBookMark::ProcessingNormal()
{
  static const TString OutputDir(
    "/group/had/sks/E42/JPARC2021May/e42_2021may/bookmark");
  static std::ofstream ofs(Form("%s/run%05d_bookmark.dat",
                                OutputDir.Data(),
                                gUnpacker.get_root()->get_run_number()),
                           std::ios::binary);
  static ULong64_t first_bookmark = gUnpacker.get_istream_bookmark();
  static ULong64_t buf;
  buf = gUnpacker.get_istream_bookmark() - first_bookmark;
  ofs.write(reinterpret_cast<char*>(&buf), sizeof(buf));

  event.runnum = gUnpacker.get_run_number();
  event.evnum = gUnpacker.get_run_number();
  auto dsize = gUnpacker.get_node_header(gUnpacker.get_fe_id("k18eb"),
                                         DAQNode::k_data_size);
  event.dsize = dsize;
  HF1(1, dsize);

  return true;
}

//_____________________________________________________________________________
Bool_t
UserBookMark::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserBookMark;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  gUnpacker.disable_istream_bookmark();

  HB1(1, "EB Data Size", 1000, 0, 1e6);

  HBTree("bookmark","tree of BookMark");
  tree->Branch("runnum", &event.runnum, "runnum/I");
  tree->Branch("evnum", &event.evnum, "evnum/I");
  tree->Branch("dsize", &event.dsize, "dsize/I");

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

// -*- C++ -*-

#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TSystem.h>

#include <std_ostream.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DebugCounter.hh"
#include "DebugTimer.hh"
#include "DeleteUtility.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"

namespace
{
using hddaq::unpacker::GUnpacker;
auto& gConf     = ConfMan::GetInstance();
auto& gCounter  = debug::ObjectCounter::GetInstance();
auto& gUnpacker = GUnpacker::get_instance();
enum EArg
{
  kArgProcess,
  kArgConfFile,
  kArgInFile,
  kArgOutFile,
  kArgc
};
}

TROOT theROOT("k18analyzer", "k18analyzer");

//______________________________________________________________________________
int
main(int argc, char **argv)
{
  std::vector<TString> arg(argv, argv + argc);
  const TString& process = arg[kArgProcess];
  if(argc!=kArgc){
    hddaq::cout << "#D Usage: " << gSystem->BaseName(process)
  		<< " [analyzer config file]"
  		<< " [data input stream]"
  		<< " [output root file]"
  		<< std::endl;
    return EXIT_SUCCESS;
  }

  debug::Timer timer("[::main()] End of Analyzer");

  const TString& conf_file = arg[kArgConfFile];
  const TString& in_file   = arg[kArgInFile];
  const TString& out_file  = arg[kArgOutFile];

  hddaq::cout << "[::main()] recreate root file : " << out_file << std::endl;
  new TFile(out_file, "recreate");

  if(!gConf.Initialize(conf_file) || !gConf.InitializeUnpacker())
    return EXIT_FAILURE;

  gUnpacker.set_istream(in_file.Data());
  gUnpacker.initialize();

  CatchSignal::Set(SIGINT);

  for(; !gUnpacker.eof() && !CatchSignal::Stop(); ++gUnpacker){
    VEvent* event = gConf.EventAllocator();
    event->ProcessingBegin();
    event->ProcessingNormal();
    event->ProcessingEnd();
    delete event;
    gCounter.check();
  }
  gConf.Finalize();

  gFile->Write();
  gFile->Close();

  return EXIT_SUCCESS;
}

// -*- C++ -*-

#include "CatchSignal.hh"

#include <Unpacker.hh>
#include <UnpackerManager.hh>
#include <escape_sequence.hh>
#include <std_ostream.hh>

#include "FuncName.hh"

namespace CatchSignal
{
namespace
{
using namespace hddaq::unpacker;
const auto& gUnpacker = GUnpacker::get_instance();
Bool_t user_stop = false;
}

//_____________________________________________________________________________
Bool_t
Stop()
{
  return user_stop;
}

//_____________________________________________________________________________
void
Catch(Int_t sig)
{
  std::signal(SIGINT, SIG_IGN);
  user_stop = true;
  if(gUnpacker.get_root()->is_esc_on()){
    hddaq::cout << esc::k_yellow
                << FUNC_NAME << " exit process by signal " << sig
                << esc::k_default_color << std::endl;
  }else{
    hddaq::cout << FUNC_NAME << " exit process by signal " << sig
                << std::endl;
  }
}

//_____________________________________________________________________________
void
Set(Int_t sig)
{
  std::signal(sig, Catch);
}
}

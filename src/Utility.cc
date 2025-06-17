// -*- C++ -*-

#include <DAQNode.hh>
#include <Unpacker.hh>
#include <UnpackerManager.hh>

#include "Utility.hh"

namespace
{
using hddaq::unpacker::GUnpacker;
using hddaq::unpacker::DAQNode;
const auto& gUnpacker = GUnpacker::get_instance();
}

namespace utility
{
//_____________________________________________________________________________
UInt_t
EBDataSize( void )
{
  static const Int_t NodeIdEB = gUnpacker.get_fe_id("k18eb");
  return gUnpacker.get_node_header(NodeIdEB, DAQNode::k_data_size);
}

//_____________________________________________________________________________
UInt_t
UnixTime( void )
{
  static const Int_t NodeIdEB = gUnpacker.get_fe_id("k18eb");
  return gUnpacker.get_node_header(NodeIdEB, DAQNode::k_unix_time);
}

//_____________________________________________________________________________
TTimeStamp
TimeStamp( void )
{
  TTimeStamp s( UnixTime() );
  s.Add( -s.GetZoneOffset() );
  return s;
}
}

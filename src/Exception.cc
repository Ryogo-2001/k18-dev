// -*- C++ -*-

#include "Exception.hh"

#include <TString.h>

#include <escape_sequence.hh>
#include <std_ostream.hh>

#include "FuncName.hh"

//_____________________________________________________________________________
Exception::Exception(const TString& msg)
  : m_msg()
{
  m_msg = hddaq::unpacker::esc::k_yellow
    // + FUNC_NAME + " "
    + msg
    + hddaq::unpacker::esc::k_default_color;
}

//_____________________________________________________________________________
Exception::~Exception() throw()
{
}

//_____________________________________________________________________________
void
Exception::hoge(const TString& arg) const
{
  hddaq::cout << FUNC_NAME << " " << arg << std::endl;
}

//_____________________________________________________________________________
const Char_t*
Exception::what() const throw()
{
  return m_msg.Data();
}

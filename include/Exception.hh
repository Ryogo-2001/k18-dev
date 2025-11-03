// -*- C++ -*-

#ifndef EXCEPTION_HH
#define EXCEPTION_HH

#include <exception>
#include <stdexcept>

#include <TString.h>

//_____________________________________________________________________________
class Exception : public std::exception
{
public:
  static const TString& ClassName();
  Exception(const TString& msg);
  virtual ~Exception() throw();

private:
  TString m_msg;

public:
  virtual void          hoge(const TString& arg="") const;
  virtual const Char_t* what() const throw();
};

//_____________________________________________________________________________
inline const TString&
Exception::ClassName()
{
  static TString s_name("Exception");
  return s_name;
}

#endif

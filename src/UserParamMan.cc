// -*- C++ -*-

#include "UserParamMan.hh"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include <TFile.h>
#include <TMath.h>
#include <TNamed.h>

#include <std_ostream.hh>

#include "FuncName.hh"

namespace
{
const Double_t default_value = TMath::QuietNaN();
}

//_____________________________________________________________________________
UserParamMan::UserParamMan()
  : m_is_ready(false),
    m_use_default(false), // if no parameter, throw exception
    m_file_name(),
    m_param_map(),
    m_buf(),
    m_object()
{
}

//_____________________________________________________________________________
UserParamMan::~UserParamMan()
{
  if(m_object) delete m_object;
}

//_____________________________________________________________________________
void
UserParamMan::AddObject()
{
  if(m_object) delete m_object;
  m_object = new TNamed("user", m_buf.Data());
  m_object->Write();
}

//_____________________________________________________________________________
Bool_t
UserParamMan::Initialize()
{
  std::ifstream ifs(m_file_name);
  if(!ifs.is_open()){
    hddaq::cerr << FUNC_NAME << " "
		<< "No such parameter file : " << m_file_name << std::endl;
    return false;
  }

  m_buf = "\n";

  TString line;
  while(ifs.good() && line.ReadLine(ifs)){
    m_buf += line + "\n";
    if(line.IsNull() || line[0]=='#') continue;
    std::istringstream input_line(line.Data());
    TString key;
    input_line >> key;
    ParamArray param_array;
    Double_t   param;
    while(input_line >> param){
      param_array.push_back(param);
    }
    m_param_map[key] = param_array;
  }

  AddObject();

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
UserParamMan::Initialize(const TString& filename)
{
  m_file_name = filename;
  return Initialize();
};

//_____________________________________________________________________________
Int_t
UserParamMan::GetSize(const TString& key) const
{
  PIterator itr = m_param_map.find(key);
  if(itr == m_param_map.end()){
    Print(m_file_name);
    hddaq::cerr << FUNC_NAME << " "
		<< "No such key : " << key << std::endl;
    return 0;
  }

  return itr->second.size();
}

//_____________________________________________________________________________
Double_t
UserParamMan::GetParameter(const TString& key, Int_t i) const
{
  std::stringstream param;
  param << key << "(" << i << ")";

  PIterator itr = m_param_map.find(key);

  if(itr == m_param_map.end()){
    Print(m_file_name);
    if(m_use_default){
      hddaq::cerr << FUNC_NAME
                  << "set default value : " << param.str() << " -> "
                  << default_value << std::endl;
      return default_value;
    } else {
      throw std::out_of_range((FUNC_NAME+" No such key : "+key).Data());
    }
  }

  if(i+1 > itr->second.size()){
    Print(m_file_name);
    if(m_use_default){
      hddaq::cerr << FUNC_NAME
                  << "set default value : " << param.str() << " -> "
                  << default_value << std::endl;
      return default_value;
    } else {
      throw std::out_of_range((FUNC_NAME+" No such key : "+key).Data());
    }
  }

  return itr->second.at(i);
}

//_____________________________________________________________________________
void
UserParamMan::Print(const TString& arg) const
{
  hddaq::cout << FUNC_NAME << " " << arg << std::endl;

  const Int_t w = 20;
  for(PIterator itr=m_param_map.begin(), end=m_param_map.end();
       itr != end; ++itr){
    hddaq::cout << " key = " << std::setw(w) << std::left
		<< itr->first << itr->second.size() << " : ";
    for(Int_t i=0, n=itr->second.size(); i<n; ++i){
      hddaq::cout << std::setw(5) << std::right
		  << itr->second.at(i) << " ";
    }
    hddaq::cout << std::endl;
  }
}

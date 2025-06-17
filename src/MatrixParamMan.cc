// -*- C++ -*-

#include "MatrixParamMan.hh"

#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>

#include <std_ostream.hh>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "FuncName.hh"

//_____________________________________________________________________________
MatrixParamMan::MatrixParamMan()
  : m_is_ready(false),
    m_file_name_2d(),
    m_file_name_3d()
{
  // 2D
  m_enable_2d.resize(NumOfSegTOF);
  for(Int_t i=0; i<NumOfSegTOF; ++i){
    m_enable_2d[i].resize(NumOfSegSCH);
  }
  // 3D
  m_enable_3d.resize(NumOfSegTOF);
  for(Int_t i=0; i<NumOfSegTOF; ++i){
    m_enable_3d[i].resize(NumOfSegSCH);
    for(Int_t j=0; j<NumOfSegSCH; ++j){
      m_enable_3d[i][j].resize(NumOfSegSFT_Mtx);
    }
  }
}

//_____________________________________________________________________________
MatrixParamMan::~MatrixParamMan()
{
}

//_____________________________________________________________________________
Bool_t
MatrixParamMan::Initialize()
{
  return Initialize(m_file_name_2d, m_file_name_3d);
}

//_____________________________________________________________________________
Bool_t
MatrixParamMan::Initialize(const TString& filename_2d,
                           const TString& filename_3d)
{
  {
    std::ifstream ifs(filename_2d);
    if(!ifs.is_open()){
      hddaq::cerr << FUNC_NAME << " "
		  << "No such 2D parameter file : " << filename_2d << std::endl;
      std::exit(EXIT_FAILURE);
    }

    TString line;
    Int_t tofseg = 0;
    while(!ifs.eof() && line.ReadLine(ifs)){
      if(line.IsNull()) continue;

      std::string param[2];
      std::istringstream iss(line.Data());
      iss >> param[0] >> param[1];

      if(param[0].at(0) == '#') continue;

      if(param[0].substr(2).at(0)=='0')
	param[0] = param[0].substr(3);
      else
	param[0] = param[0].substr(2);

      Int_t channel = std::strtol(param[0].c_str(), NULL, 0);
      Int_t enable  = std::strtol(param[1].substr(0,1).c_str(), NULL, 0);

      m_enable_2d[tofseg][channel] = enable;
      if(channel==NumOfSegSCH-1) ++tofseg;
    }
  }

  {
    std::ifstream ifs(filename_3d);
    if(!ifs.is_open()){
      hddaq::cerr << FUNC_NAME << " "
		  << "No such 3D parameter file : " << filename_3d << std::endl;
      std::exit(EXIT_FAILURE);
    }

    TString line;
    Int_t tofseg = 0;
    while(!ifs.eof() && line.ReadLine(ifs)){
      if(line.IsNull()) continue;

      std::string param[2];
      std::istringstream iss(line.Data());
      iss >> param[0] >> param[1];

      if(param[0].at(0) == '#') continue;

      if(param[0].substr(2).at(0)=='0')
	param[0] = param[0].substr(3);
      else
	param[0] = param[0].substr(2);

      Int_t channel = std::strtol(param[0].c_str(), NULL, 0);
      Int_t enable[NumOfSegSFT_Mtx] = {};
      for(Int_t i=0; i<NumOfSegSFT_Mtx; ++i){
	enable[i] = std::strtol(param[1].substr(i,1).c_str(), NULL, 0);
	m_enable_3d[tofseg][channel][i] = enable[i];
      }

      if(channel==NumOfSegSCH-1) ++tofseg;
    }
  }

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
MatrixParamMan::IsAccept(Int_t detA, Int_t detB) const
{
  if(m_enable_2d.size() <= detA){
    hddaq::cerr << FUNC_NAME << " detA is too much : "
		<< detA << "/" << m_enable_2d.size() << std::endl;
    return false;
  }

  if(m_enable_2d[detA].size()<=detB){
    hddaq::cerr << FUNC_NAME << " detB is too much : "
		<< detB << "/" << m_enable_2d[detA].size() << std::endl;
    return false;
  }

  return (m_enable_2d[detA][detB] != 0.);
}

//_____________________________________________________________________________
Bool_t
MatrixParamMan::IsAccept(Int_t detA, Int_t detB, Int_t detC) const
{
  if(m_enable_3d.size() <= detA){
    hddaq::cerr << FUNC_NAME << " detA is too much : "
		<< detA << "/" << m_enable_3d.size() << std::endl;
    return false;
  }

  if(m_enable_3d[detA].size() <= detB){
    hddaq::cerr << FUNC_NAME << " detB is too much : "
		<< detB << "/" << m_enable_3d[detA].size() << std::endl;
    return false;
  }

  if(m_enable_3d[detA][detB].size() <= detC){
    hddaq::cerr << FUNC_NAME << " detC is too much : "
		<< detC << "/" << m_enable_3d[detA][detB].size() << std::endl;
    return false;
  }

  return (m_enable_3d[detA][detB][detC] != 0.);
}

//_____________________________________________________________________________
void
MatrixParamMan::Print2D(const TString& arg) const
{
  hddaq::cout << FUNC_NAME << " " << arg << std::endl;
  for(Int_t i=NumOfSegTOF-1; i>=0; --i){
    hddaq::cout << " detA = " << std::setw(2)
		<< std::right << i << " : ";
    for(Int_t j=0; j<NumOfSegSCH; ++j){
      hddaq::cout << m_enable_2d[i][j];
    }
    hddaq::cout << std::endl;
  }
}

//_____________________________________________________________________________
void
MatrixParamMan::Print3D(const TString& arg) const
{
  hddaq::cout << FUNC_NAME << " " << arg << std::endl;
  for(Int_t k=0; k<NumOfSegSFT_Mtx; ++k){
    hddaq::cout << " detC = " << std::setw(2)
		<< std::right << k << std::endl;
    for(Int_t i=NumOfSegTOF-1; i>=0; --i){
      hddaq::cout << " detA = " << std::setw(2)
		  << std::right << i << " : ";
      for(Int_t j=0; j<NumOfSegSCH; ++j){
	hddaq::cout << m_enable_3d[i][j][k];
      }
      hddaq::cout << std::endl;
    }
  }
}

//_____________________________________________________________________________
void
MatrixParamMan::SetMatrix2D(const TString& file_name)
{
  m_file_name_2d = file_name;
}

//_____________________________________________________________________________
void
MatrixParamMan::SetMatrix3D(const TString& file_name)
{
  m_file_name_3d = file_name;
}

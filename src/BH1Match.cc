// -*- C++ -*-

#include "BH1Match.hh"

#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>

#include "DetectorID.hh"
#include "FuncName.hh"
#include "Exception.hh"

//_____________________________________________________________________________
BH1Match::Param::Param()
  : m_xmin(0.),
    m_xmax(0.)
{
}

//_____________________________________________________________________________
BH1Match::Param::~Param()
{
}

//_____________________________________________________________________________
void
BH1Match::Param::Print(std::ostream& ost) const
{
  static const Int_t w = 5;
  ost << std::right << std::setw(w) << m_seg << ":"
      << " (" << std::right
      << std::setw(w) << m_xmin << " "
      << std::setw(w) << m_xmax << ")"
      << std::endl;
}

//_____________________________________________________________________________
BH1Match::BH1Match()
  : m_param(2*NumOfSegBH1-1)
{
  m_status.reset();
}

//_____________________________________________________________________________
BH1Match::~BH1Match()
{
}

//_____________________________________________________________________________
Bool_t
BH1Match::Initialize(const TString& file_name)
{
  std::ifstream ifs(file_name);
  if(!ifs.is_open()){
    std::cerr << FUNC_NAME << " No such parameter file("
              << file_name << ")" << std::endl;
    return false;
  }

  TString line;
  while(ifs.good() && line.ReadLine(ifs)){
    if(line.IsNull() || line[0] == '#') continue;
    std::istringstream iss(line.Data());
    std::istream_iterator<Double_t> itr_begin(iss);
    std::istream_iterator<Double_t> itr_end;
    std::vector<Double_t> cont(itr_begin, itr_end);
    if(cont.size() != kNParam){
      std::cerr << FUNC_NAME << "\n"
		<< " Number of parameters: " << cont.size() << "\n"
		<< " Required: " << kNParam << std::endl;
      throw Exception(FUNC_NAME + " invalid parameter.");
    }

    const Double_t bh1seg = cont[kBH1Segment];
    const Int_t i_bh1seg  = static_cast<Int_t>(2*cont[kBH1Segment]);
    const Double_t xmin   = cont[kXMin];
    const Double_t xmax   = cont[kXMax];

    m_param.at(i_bh1seg).m_seg  = bh1seg;
    m_param.at(i_bh1seg).m_xmin = xmin;
    m_param.at(i_bh1seg).m_xmax = xmax;
  }// read line

  if(m_status[kVerbose]) Print();
  m_status.set(kReady);

  return true;
}

//_____________________________________________________________________________
Bool_t
BH1Match::Judge(Double_t bft_xpos, Double_t bh1seg)
{
  if(!m_status[kReady]) throw Exception(FUNC_NAME + " is not initialized.");

  const Int_t i_bh1seg = static_cast<Int_t>(2*bh1seg);
  const Double_t xmin  = m_param.at(i_bh1seg).m_xmin;
  const Double_t xmax  = m_param.at(i_bh1seg).m_xmax;

  if(m_status[kVerbose]){
    hddaq::cout << FUNC_NAME << std::endl;
    m_param.at(i_bh1seg).Print();
    hddaq::cout << std::fixed;
    hddaq::cout << " BFT pos (" << std::setw(5)
		<< std::setprecision(2) << bft_xpos
		<< ") : ";
    if(xmin < bft_xpos && bft_xpos < xmax){
      hddaq::cout << "adopt" << std::endl;
    } else {
      hddaq::cout << "reject" << std::endl;
    }
  }

  return (xmin < bft_xpos && bft_xpos < xmax);
}

//_____________________________________________________________________________
void
BH1Match::Print(std::ostream& ost) const
{
  ost << FUNC_NAME << "\n"
      << "     (xmin xmax)" << std::endl;
  for(const auto& param : m_param) param.Print(ost);
}

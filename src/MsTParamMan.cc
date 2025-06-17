// -*- C++ -*-

#include "MsTParamMan.hh"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <std_ostream.hh>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "Exception.hh"
#include "FuncName.hh"

namespace
{
const Int_t NumOfSegDetA = NumOfSegTOF;
const Int_t NumOfSegDetB = NumOfSegSCH;
}

//_____________________________________________________________________________
MsTParamMan::MsTParamMan( void )
  : m_is_ready(false),
    m_nA(NumOfSegDetA),
    m_nB(NumOfSegDetB),
    m_low_threshold(),
    m_high_threshold()
{
  m_low_threshold.resize( NumOfSegDetA );
  m_high_threshold.resize( NumOfSegDetA );
  for( Int_t i=0; i<NumOfSegDetA; ++i ){
    m_low_threshold[i].resize( NumOfSegDetB );
    m_high_threshold[i].resize( NumOfSegDetB );
  }
}

//_____________________________________________________________________________
MsTParamMan::~MsTParamMan( void )
{
}

//_____________________________________________________________________________
Double_t
MsTParamMan::GetLowThreshold( Int_t detA, Int_t detB ) const
{
  if( detA<0 || m_nA<=detA ) return false;
  if( detB<0 || m_nB<=detB ) return false;
  return m_low_threshold[detA][detB];
}

//_____________________________________________________________________________
Double_t
MsTParamMan::GetHighThreshold( Int_t detA, Int_t detB ) const
{
  if( detA<0 || m_nA<=detA ) return false;
  if( detB<0 || m_nB<=detB ) return false;
  return m_high_threshold[detA][detB];
}

//_____________________________________________________________________________
Bool_t
MsTParamMan::Initialize( const TString& filename )
{
  std::ifstream ifs( filename );
  if( !ifs.is_open() ){
    hddaq::cerr << FUNC_NAME << " "
                << "No such parameter file : " << filename << std::endl;
    std::exit( EXIT_FAILURE );
  }

  TString line;
  Int_t tofseg = 0;
  Int_t LorH   = 0;
  while( ifs.good() && line.ReadLine( ifs ) ){
    if( line.IsNull() || line[0] == '#' ) continue;
    TString param[NumOfSegDetB];
    std::istringstream iss( line.Data() );

    for( Int_t i=0; i<NumOfSegDetB; ++i ){
      iss >> param[i];
    }
    if( param[0][0] == '#' ) continue;
    if( param[1] == "Mem_Thr_Low" ){
      tofseg =  0;
      LorH   = -1;
      continue;
    }
    if( param[1] == "Mem_Thr_Hi" ){
      tofseg =  0;
      LorH   =  1;
      continue;
    }

    if( param[NumOfSegDetB-1].IsNull() ) continue;
    for( Int_t i=0; i<NumOfSegDetB; ++i ){
      param[i].ReplaceAll( "(", " " );
      param[i].ReplaceAll( ")", " " );
      param[i].ReplaceAll( ",", " " );
      if( LorH == -1 )
        m_low_threshold[tofseg][i] = param[i].Atof();
      if( LorH ==  1 )
        m_high_threshold[tofseg][i] = param[i].Atof();
    }
    ++tofseg;
  }

#if 0
  // Low
  hddaq::cout << "Low Threshold" << std::endl;
  for( Int_t i=0; i<NumOfSegDetA; ++i ){
    hddaq::cout << i << "\t:";
    for( Int_t j=0; j<NumOfSegDetB; ++j ){
      hddaq::cout << " " << m_low_threshold[i][j];
    }
    hddaq::cout << std::endl;
  }
  // High
  hddaq::cout << "High Threshold" << std::endl;
  for( Int_t i=0; i<NumOfSegDetA; ++i ){
    hddaq::cout << i << "\t:";
    for( Int_t j=0; j<NumOfSegDetB; ++j ){
      hddaq::cout << " " << m_high_threshold[i][j];
    }
    hddaq::cout << std::endl;
  }
#endif

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
MsTParamMan::IsAccept( Int_t detA, Int_t detB, Int_t tdc ) const
{
  if( !m_is_ready ){
    throw Exception( FUNC_NAME+" "+ClassName()+" is not initialized" );
    // return false;
  }
  if( m_nA<=detA || m_nB<=detB ){
    std::ostringstream oss;
    oss << " detA/detB is out of range : "
	<< std::setw(2) << detA << " "
	<< std::setw(2) << detB;
    throw std::out_of_range(FUNC_NAME+oss.str());
    // return false;
  }

  Int_t low  = (Int_t)m_low_threshold[detA][detB];
  Int_t high = (Int_t)m_high_threshold[detA][detB];
  return ( low < tdc && tdc < high );
}

//_____________________________________________________________________________
void
MsTParamMan::Print( const TString& arg ) const
{
  hddaq::cout << FUNC_NAME << " " << arg << std::endl;

  for( Int_t iA=0; iA<m_nA; ++iA ){
    hddaq::cout << std::setw(2) << iA << ": ";
    Int_t iBok = 0;
    for( Int_t iB=0; iB<m_nB; ++iB ){
      // ofs << iA << "\t" << iB << "\t"
      // 	  << m_low_threshold[iA][iB] << "\t"
      // 	  << m_high_threshold[iA][iB] << std::endl;
      if( m_low_threshold[iA][iB]==0 &&
	  m_high_threshold[iA][iB]==0 )
	continue;

      iBok++;
      hddaq::cout << std::setw(2) << iB << " "
		  << std::setw(4) << m_low_threshold[iA][iB] << " "
		  << std::setw(4) << m_high_threshold[iA][iB] << "  ";
      if( iBok%5==0 ){
	hddaq::cout << std::endl << "    ";
      }
    }
    hddaq::cout << std::endl;
  }
}

//_____________________________________________________________________________
void
MsTParamMan::Print( Int_t detA, Int_t detB, Int_t tdc ) const
{
  if( detA<0 || m_nA<=detA ) return;
  if( detB<0 || m_nB<=detB ) return;
  hddaq::cout << " detA " << std::setw(2) << detA
	      << " detB " << std::setw(2) << detB
	      << " : " << IsAccept( detA, detB, tdc ) << std::endl;
}

/**
 *  file: CFTPedCorMan.cc
 *  date: 2017.04.10
 *
 */

#include "CFTPedCorMan.hh"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <std_ostream.hh>

#include "DeleteUtility.hh"
#include "MathTools.hh"
#include "HodoParamMan.hh"
#include "HodoRawHit.hh"
#include "RawData.hh"

namespace
{
  const std::string& class_name("CFTPedCorMan");
  const HodoParamMan& gHodo = HodoParamMan::GetInstance();
}

//______________________________________________________________________________
CFTPedCorMan::CFTPedCorMan( void )
  : m_is_ready(false), m_file_name("")
{
}

//______________________________________________________________________________
CFTPedCorMan::~CFTPedCorMan()
{

}

//______________________________________________________________________________
bool
CFTPedCorMan::Initialize( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_ready ){
    hddaq::cerr << "#W " << func_name
		<< " already initialied" << std::endl;
    return false;
  }

  std::ifstream ifs( m_file_name.c_str() );
  if( !ifs.is_open() ){
    hddaq::cerr << "#E " << func_name << " file open fail : "
		<< m_file_name << std::endl;
    return false;
  }

  std::string line;
  while( ifs.good() && std::getline(ifs,line) ){
    if( line.empty() || line[0]=='#' ) continue;
    std::istringstream input_line( line );

    int layer = -1, segment = -1;
    double sigma = 999.;
    int cor_layer = -1, cor_segment = -1;
    double thr;
    double p0, p1, p2, p3;

    if (input_line >> layer >> segment >> sigma >> cor_layer >> cor_segment >> thr >> p0 >> p1 >> p2 >> p3) {

      CFTPedCorParam res;
      res.layer = layer;
      res.segment = segment;
      res.sigma = sigma;
      res.cor_layer = cor_layer;
      res.cor_segment = cor_segment;
      res.thr = thr;
      res.p0_HG = p0;
      res.p1_HG = p1;
      res.p0_LG = p2;
      res.p1_LG = p3;

      m_container[layer][segment].push_back(res);
    }
    else{
      hddaq::cerr << func_name << ": Invalid format" << std::endl
		  << " ===> " << line << std::endl;
    }
  }

  m_is_ready = true;
  return true;
}

//______________________________________________________________________________
bool
CFTPedCorMan::Initialize( const std::string& file_name )
{
  m_file_name = file_name;
  return Initialize();
}

//______________________________________________________________________________
bool
CFTPedCorMan::PedestalCorrection(int layer, int segment, double &deltaHG, double &deltaLG, RawData* rawData ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  deltaHG = -9999.;
  deltaLG = -9999.;

  if (! (segment >= 0 && segment < NumOfSegCFT[layer]) ) {
    hddaq::cerr << "#W " << func_name
		<< " invalid CFT segment, layer " << layer 
		<< ", segment " << segment 
		<< " ( " << NumOfSegCFT[layer] << " )" << std::endl;
    return false;
  }

  
  int n = m_container[layer][segment].size();
  for (int i=0; i<n; i++) {
    int cor_layer = m_container[layer][segment][i].cor_layer;
    int cor_segment = m_container[layer][segment][i].cor_segment;
    double thr = m_container[layer][segment][i].thr;
    
    const HodoRHitContainer &cont = rawData->GetCFTRawHC(cor_layer);
    int nh=cont.size();
    for( int j=0; j<nh; ++j ){
      HodoRawHit *hit=cont[j];
      if( !hit ) continue;

      int seg = hit->SegmentId();
      if (seg != cor_segment)
	continue;

      int adcHi = hit->GetAdcUp();	
      if (adcHi > thr)
	break;
      else {
	//std::cout << "Corrected by ( " << cor_layer << ", " << cor_segment
	//		  << " ) " << adcHi << std::endl;

	double p0 = m_container[layer][segment][i].p0_HG;
	double p1 = m_container[layer][segment][i].p1_HG;
	deltaHG = -(p0 + p1*adcHi);

	double p2 = m_container[layer][segment][i].p0_LG;
	double p3 = m_container[layer][segment][i].p1_LG;
	deltaLG = -(p2 + p3*adcHi);

	return true;
      }
    }
  }

  return true;

}

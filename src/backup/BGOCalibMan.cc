#include "BGOCalibMan.hh"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <std_ostream.hh>

//______________________________________________________________________________

BGOCalibMan::BGOCalibMan()
  : m_is_ready(false),
    m_file_name("")
{
}


BGOCalibMan::~BGOCalibMan()
{

}


//______________________________________________________________________________
bool
BGOCalibMan::Initialize( void )
{
  const std::string& class_name("BGOCalibMan");

  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_ready ){
    std::cerr << "#W " << func_name
		<< " already initialied" << std::endl;
    return false;
  }


  std::ifstream ifs( m_file_name.c_str() );
  if( !ifs.is_open() ){
    std::cerr << "#E " << func_name << " file open fail : "
		<< m_file_name << std::endl;
    return false;
  }

  for (int i=0; i<NumOfSegBGO; i++) {
    m_Econt[i].clear();
    m_Hcont[i].clear();
  }

  int invalid=0;
  std::string line;
  while( ifs.good() && std::getline( ifs, line ) ){
    ++invalid;
    if( line[0]=='#' || line.empty() ) continue;
    std::istringstream input_line( line );
    int seg=-1, npoint=-1;
    double E, PH;

    if( input_line >> seg >> npoint >> E >> PH ){
      m_Econt[seg].push_back(E);
      m_Hcont[seg].push_back(PH);
    }
  }

  m_is_ready = true;
  return true;
}

//______________________________________________________________________________
bool
BGOCalibMan::Initialize( const std::string& file_name )
{
  m_file_name = file_name;
  return Initialize();
}

//______________________________________________________________________________
bool
BGOCalibMan::GetEnergy( int seg, double pulse_height, double &energy ) const
{
  if (!(seg>=0 && seg<NumOfSegBGO)) {
    return false;
  }

  int icut = 0;

  int Cut_Num = m_Hcont[seg].size();
  if(pulse_height<m_Hcont[seg][Cut_Num-1]){
    while((pulse_height>m_Hcont[seg][icut]) && (icut<Cut_Num-1)){
      icut+=1;
    }
  }else{
    icut = Cut_Num-1;
  }

  if(icut>0){
    energy  = m_Econt[seg][icut-1] + (m_Econt[seg][icut]-m_Econt[seg][icut-1])*(pulse_height-m_Hcont[seg][icut-1])/(m_Hcont[seg][icut]-m_Hcont[seg][icut-1]);
  }else{
    energy = 0;
  }

  return true;
}

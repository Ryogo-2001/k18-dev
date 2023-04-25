#ifndef BGO_CALIB_MAN_HH
#define BGO_CALIB_MAN_HH

#include <numeric>
#include <vector>
#include <TROOT.h>
#include <string>
#include <fstream>
#include <cmath>
#include <string>

#include "DetectorID.hh"

class BGOCalibMan{
 private:
  bool        m_is_ready;
  std::string m_file_name;
  
  std::vector <double> m_Econt[NumOfSegBGO];
  std::vector <double> m_Hcont[NumOfSegBGO];

 public:
  ~BGOCalibMan();
  BGOCalibMan();

  static BGOCalibMan& GetInstance( void );
  static const std::string& ClassName( void );

  bool Initialize( void );
  bool Initialize( const std::string& file_name );

  std::string GetFileName( void ) const { return m_file_name; };
  
  bool GetEnergy( int seg, double pulse_height, double &energy ) const;

};

//______________________________________________________________________________
inline BGOCalibMan&
BGOCalibMan::GetInstance( void )
{
  static BGOCalibMan g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const std::string&
BGOCalibMan::ClassName( void )
{
  static std::string g_name("BGOCalibMan");
  return g_name;
}



#endif  //INC_CTEMP

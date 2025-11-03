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
  Bool_t      m_is_ready;
  TString   m_file_name;
  
  std::vector <Double_t> m_Econt[NumOfSegBGO];
  std::vector <Double_t> m_Hcont[NumOfSegBGO];

 public:
  ~BGOCalibMan();
  BGOCalibMan();

  static BGOCalibMan& GetInstance( void );
  static const TString& ClassName( void );

  Bool_t Initialize( void );
  Bool_t Initialize( const TString& file_name );

  TString GetFileName( void ) const { return m_file_name; };
  
  Bool_t GetEnergy( Int_t seg, Double_t pulse_height, Double_t &energy ) const;

};

//______________________________________________________________________________
inline BGOCalibMan&
BGOCalibMan::GetInstance( void )
{
  static BGOCalibMan g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const TString&
BGOCalibMan::ClassName( void )
{
  static TString g_name("BGOCalibMan");
  return g_name;
}



#endif  //INC_CTEMP

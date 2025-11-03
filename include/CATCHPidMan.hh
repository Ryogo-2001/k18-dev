#ifndef CATCH_PID_MAN_HH
#define CATCH_PID_MAN_HH

#include <numeric>
#include <vector>
#include <TROOT.h>
#include <string>
#include <fstream>
#include <cmath>
#include <string>

#include "DetectorID.hh"

class CATCHPidMan{
 private:
  Bool_t   m_is_ready;
  TString  m_file_name;
  
  Double_t m_ParLow[4];
  Double_t m_ParHigh[4];

 public:
  ~CATCHPidMan();
  CATCHPidMan();

  static CATCHPidMan& GetInstance( void );
  static const TString& ClassName( void );

  Bool_t Initialize( void );
  Bool_t Initialize( const TString& file_name );

  TString GetFileName( void ) const { return m_file_name; };
  
  Bool_t CheckProton( Double_t BGO_E, Double_t CFT_dE, Double_t &delta ) const;
  Bool_t CheckPi( Double_t BGO_E, Double_t CFT_dE) const;

};

//______________________________________________________________________________
inline CATCHPidMan&
CATCHPidMan::GetInstance( void )
{
  static CATCHPidMan g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const TString&
CATCHPidMan::ClassName( void )
{
  static TString g_name("CATCHPidMan");
  return g_name;
}



#endif  //INC_CTEMP

#ifndef EVENT_SELECT_MAN_HH
#define EVENT_SELECT_MAN_HH

#include <numeric>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <string>
#include <TString.h>

class EventSelectMan{
 private:
  Bool_t        m_is_ready;
  TString m_file_name;
  
  std::vector <Int_t> m_EventNumCont;

  mutable Int_t         m_current_point;

 public:
  ~EventSelectMan();
  EventSelectMan();

  static EventSelectMan& GetInstance( void );
  static const TString& ClassName( void );

  Bool_t Initialize( void );
  Bool_t Initialize( const TString& file_name );

  TString GetFileName( void ) const { return m_file_name; };
  
  Bool_t RunNumCheck( Int_t runNum ) const;
  Bool_t IsGood( Int_t eventNum ) const;

};

//______________________________________________________________________________
inline EventSelectMan&
EventSelectMan::GetInstance( void )
{
  static EventSelectMan g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const TString&
EventSelectMan::ClassName( void )
{
  static TString g_name("EventSelectMan");
  return g_name;
}



#endif  //INC_CTEMP

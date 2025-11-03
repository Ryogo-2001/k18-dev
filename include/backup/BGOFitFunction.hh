#ifndef BGO_FIT_FUNCTION_HH
#define BGO_FIT_FUNCTION_HH

#include <numeric>
#include <vector>
#include <TROOT.h>
#include <string>
#include <fstream>
#include <cmath>
#include <string>

class BGOTemplate;

class BGOFitFunction{
 private:

 public:
  ~BGOFitFunction();
  BGOFitFunction(std::string filename);
  BGOTemplate *m_tempFunc;
  Double_t operator()(Double_t *x, Double_t *par);
  

};

class BGOTemplateManager{
 private:
  bool        m_is_ready;
  std::string m_file_name;
  BGOFitFunction *m_fitFunction;  

 public:
  ~BGOTemplateManager();
  BGOTemplateManager();

  static BGOTemplateManager& GetInstance( void );
  static const std::string& ClassName( void );

  bool Initialize( void );
  bool Initialize( const std::string& file_name );

  std::string GetFileName( void ) const { return m_file_name; };
  
  BGOFitFunction* GetFitFunction( void ) const { return m_fitFunction; };
};

//______________________________________________________________________________
inline BGOTemplateManager&
BGOTemplateManager::GetInstance( void )
{
  static BGOTemplateManager g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const std::string&
BGOTemplateManager::ClassName( void )
{
  static std::string g_name("BGOTemplateManager");
  return g_name;
}



#endif  //INC_CTEMP

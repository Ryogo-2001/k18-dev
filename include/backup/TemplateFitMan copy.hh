#ifndef TEMPLATE_FIT_MAN_HH
#define TEMPLATE_FIT_MAN_HH

#include <numeric>
#include <vector>
#include <TROOT.h>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#include "DetectorID.hh"

class TemplateFitFunction{
 private:
  std::vector<Double_t> m_tempx;
  std::vector<Double_t> m_tempy;
  Double_t m_area;
  Int_t  m_sample_num;
  Double_t m_interval;
  Int_t m_center;

 public:
  ~TemplateFitFunction();
  TemplateFitFunction(TString filename);
  Double_t operator()(Double_t *x, Double_t *par);
  Double_t GetTempX(Int_t i){return m_tempx[i];}
  Double_t GetTempY(Int_t i){return m_tempy[i];}
  Int_t    GetSampleNum(){return m_sample_num;}
  Double_t GetArea(){return m_area;}
  Double_t GetTemplateFunction(Double_t x);
  Double_t myTemp(Double_t *x, Double_t *par);
};

using TempFuncContainer  = std::vector<TemplateFitFunction*>;

class TemplateFitMan{
 private:
  bool        m_is_ready;
  TString     m_file_name;
  TemplateFitFunction *m_fitFunction[NumOfSegBGO];
  template <typename T> using map_t = std::map<TString, T>;
  map_t<TempFuncContainer> m_tempfunc_collection;
  bool        m_flag_ch14;
 public:
  ~TemplateFitMan();
  TemplateFitMan();

  static TemplateFitMan& GetInstance( void );
  static const TString& ClassName( void );

  bool Initialize( void );
  bool Initialize( const TString& file_name );
  const TempFuncContainer& GetHitContainer(const TString& name) const;

  
  TString GetFileName( void ) const { return m_file_name; };
  
  TemplateFitFunction* GetFitFunction( int seg ) const { return m_fitFunction[seg]; };
  bool GetCh14Flag( void  ) const { return m_flag_ch14;};
  std::vector<TString> split(const TString &s, char delim); 
};

//______________________________________________________________________________
inline TemplateFitMan&
TemplateFitMan::GetInstance( void )
{
  static TemplateFitMan g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const TString&
TemplateFitMan::ClassName( void )
{
  static TString g_name("TemplateFitMan");
  return g_name;
}



#endif  //INC_CTEMP

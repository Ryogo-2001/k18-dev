#ifndef BGO_TEMPLATE_HH
#define BGO_TEMPLATE_HH

#include <numeric>
#include <vector>
#include <cstdio>
#include <string>
#include <fstream>
#include <cmath>
#include <TROOT.h>

#include <std_ostream.hh>

class BGOTemplate
{
private:
  std::vector<Double_t> m_tempx;
  std::vector<Double_t> m_tempy;
  Double_t m_area;
  Int_t  m_sample_num;
  Double_t m_interval;
  Int_t m_center;
  Int_t m_TEMPSAMP;

public:
  ~BGOTemplate();
  BGOTemplate(std::string filename);
  Double_t GetTempX(Int_t i){return m_tempx[i];}
  Double_t GetTempY(Int_t i){return m_tempy[i];}
  Int_t    GetSampleNum(){return m_sample_num;}
  Double_t GetArea(){return m_area;}
  Double_t GetTemplateFunction(Double_t x);
  Double_t myTemp(Double_t *x, Double_t *par);
};

#endif  //INC_CTEMP

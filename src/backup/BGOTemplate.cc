#include "BGOTemplate.hh"

BGOTemplate::~BGOTemplate(){
  m_tempx.clear();
  m_tempy.clear();
  m_tempx.shrink_to_fit();
  m_tempy.shrink_to_fit();

}

BGOTemplate::BGOTemplate(std::string filename)
  : m_area(0), m_sample_num(0), m_interval(0), m_center(0), m_TEMPSAMP(0)
{

  std::string str; 

  std::ifstream ifile(filename.c_str());

  Double_t xxx;
  Double_t yyy;

  if(!ifile)
    std::cout<<"Template File: "<< filename <<" is not exist"<<std::endl;

  Int_t nlines=0;

  m_area=0;

  Double_t max = 0;
  while(getline(ifile,str)){    
    sscanf(str.data(),"%lf %lf",&xxx,&yyy);
    m_tempx.push_back(xxx);
    m_tempy.push_back(yyy);

    if(fabs(xxx) <0.0001)
      max = -yyy;

    m_area+=(Double_t) yyy;
    if(fabs(xxx)< 0.00001)
      m_center =nlines;

    nlines ++;

  }


  m_sample_num = nlines;

  for(int i=0;i<m_sample_num;i++)
    m_tempy[i] /= max;
  m_interval = m_tempx[1]-m_tempx[0];

    
}

Double_t BGOTemplate::myTemp(Double_t *x, Double_t *par)
{
  
  Double_t k=x[0];
  Double_t p=par[0];

  return par[1]*GetTemplateFunction(k-p);
}


Double_t BGOTemplate::GetTemplateFunction(Double_t x)
{

  Int_t xx = Int_t (x / m_interval) +m_center;
  if(x<m_tempx[0] || x >=m_tempx[m_sample_num-1])
    return 0;


  Int_t p;
  for(p= xx-5;p<xx+5;p++){
    if(x>= m_tempx[p] && x< m_tempx[p+1]){
      break;
    }
  }
  Double_t l = (x-m_tempx[p])/(m_tempx[p+1]-m_tempx[p]);

  return (m_tempy[p]+(m_tempy[p+1]-m_tempy[p])*l); 

}


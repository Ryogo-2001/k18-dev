#include "BGOFitFunction.hh"
#include "BGOTemplate.hh"

BGOFitFunction::BGOFitFunction(std::string filename)
  : m_tempFunc(0)
{
  //CT = new CTemplate*[CC.TemplateFile.size()];
  m_tempFunc = new BGOTemplate(filename);

  //  tge = (TGraphErrors*) gra -> Clone();
}

BGOFitFunction::~BGOFitFunction()
{
  if (m_tempFunc)
    delete m_tempFunc;
}

Double_t BGOFitFunction::operator()(Double_t *x, Double_t *par)
{
  Double_t mix=0;
  for(int i=0; i<par[0] ; i++)
    mix += m_tempFunc -> myTemp(&x[0],&par[2+2*i]);
  mix += par[1];

  return mix;
}

//______________________________________________________________________________

BGOTemplateManager::BGOTemplateManager()
  : m_is_ready(false),
    m_file_name("")
{
}


BGOTemplateManager::~BGOTemplateManager()
{
  if (m_fitFunction) {
    delete m_fitFunction;
    m_fitFunction = 0;
  }

}


//______________________________________________________________________________
bool
BGOTemplateManager::Initialize( void )
{
  const std::string& class_name("BGOTemplateManager");

  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_ready ){
    hddaq::cerr << "#W " << func_name
		<< " already initialied" << std::endl;
    return false;
  }

  m_fitFunction = new BGOFitFunction(m_file_name);



  m_is_ready = true;
  return true;
}

//______________________________________________________________________________
bool
BGOTemplateManager::Initialize( const std::string& file_name )
{
  m_file_name = file_name;
  return Initialize();
}

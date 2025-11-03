#include "TemplateFitMan.hh"
#include "DeleteUtility.hh"

#include <iostream>
#include <cstdio>

#include <TSystem.h>

// -------------------- TemplateFitFunction --------------------

TemplateFitFunction::TemplateFitFunction(TString filename)
  : m_area(0),
    m_sample_num(0),
    m_interval(0),
    m_center(0)
{
  std::ifstream ifile(filename.Data());
  if (!ifile) {
    std::cerr << "Template File: " << filename << " is not exist" << std::endl;
    return;
  }

  std::string line;
  Double_t x, y;
  Double_t max_abs = 0.0;
  Int_t    nlines  = 0;
  Int_t    peak_line = 0;

  while (std::getline(ifile, line)) {
    if (line.empty()) continue;
    if (std::sscanf(line.c_str(), "%lf %lf", &x, &y) != 2) continue;

    m_tempx.push_back(x);
    m_tempy.push_back(y);
    m_area += y;

    if (std::fabs(y) > max_abs) {
      max_abs   = std::fabs(y);
      peak_line = nlines;
    }
    ++nlines;
  }

  m_sample_num = nlines;
  m_center     = peak_line;

  if (m_sample_num < 2) {
    std::cerr << "!!! ERROR: Template file " << filename
              << " has less than 2 samples." << std::endl;
    return;
  }

  //normalize
  const Double_t eps = 1.e-9;
  for (int i = 0; i < m_sample_num; ++i)
    m_tempy[i] /= (max_abs + eps);

  // sample interval
  m_interval = m_tempx[1] - m_tempx[0];
}

TemplateFitFunction::~TemplateFitFunction()
{
  m_tempx.clear();
  m_tempy.clear();
  m_tempx.shrink_to_fit();
  m_tempy.shrink_to_fit();
}

// 1 pulse 
Double_t TemplateFitFunction::myTemp(Double_t *x, Double_t *par)
{
  // x[0] ... time
  // par[0] ... time shift
  // par[1] ...  amplitude scale
  const Double_t k = x[0];
  const Double_t p = par[0];
  return par[1] * GetTemplateFunction(k - p);
}


Double_t TemplateFitFunction::GetTemplateFunction(Double_t x)
{
  if (m_sample_num < 2) return 0.0;

  
  if (x < m_tempx.front() || x >= m_tempx.back())
    return 0.0;

  
  Int_t xx = static_cast<Int_t>( (x - m_tempx.front()) / m_interval );
  xx = std::max(0, std::min(xx, m_sample_num - 2));

  
  const Int_t pmin = std::max(0,         xx - 5);
  const Int_t pmax = std::min(m_sample_num - 2, xx + 5);

  Int_t p = pmin;
  for (; p <= pmax; ++p) {
    if (x >= m_tempx[p] && x < m_tempx[p+1])
      break;
  }

  if (p > pmax) {
    
    return m_tempy[xx];
  }

  const Double_t l =
    (x - m_tempx[p]) / (m_tempx[p+1] - m_tempx[p]);

  return m_tempy[p] + (m_tempy[p+1] - m_tempy[p]) * l;
}


Double_t TemplateFitFunction::operator()(Double_t *x, Double_t *par)
{
  // par[0] = nPulse
  // par[1] = pedestal
  Double_t mix = 0.0;
  const Int_t nPulse = static_cast<Int_t>(par[0]);
  for (int i = 0; i < nPulse; ++i)
    mix += myTemp(&x[0], &par[2 + 2*i]);
  mix += par[1];
  return mix;
}

// -------------------- TemplateFitMan --------------------

TemplateFitMan::TemplateFitMan()
  : m_is_ready(false),
    m_file_name(""),
    m_tempfunc_collection(),
    m_flag_ch14(false)
{
}

TemplateFitMan::~TemplateFitMan()
{
  for (auto &kv : m_tempfunc_collection) {
    del::ClearContainer(kv.second);
  }
}

bool TemplateFitMan::Initialize( const TString& file_name )
{
  m_file_name = file_name;
  return Initialize();
}

bool TemplateFitMan::Initialize( void )
{
  const std::string class_name("TemplateFitMan");

  if (m_is_ready) {
    std::cerr << "#W [" << class_name << "::Initialize()] already initialized\n";
    return false;
  }

  
  if (m_file_name == "RAYRAWTEMP" ||
      m_file_name == "RAYRAW"     ||
      m_file_name == "RAYRAW_TEMP") {
    
    // param/RAYRAWTEMP/RAYRAWTEMP.<seg>
    m_file_name = "param/RAYRAWTEMP/RAYRAWTEMP";
  }

  {
    // BGO  template don't use now
    TString name("BGO");
    auto &cont = m_tempfunc_collection[name];
    for (auto *f : cont) delete f;
    cont.clear();
  }

  {
    // RAYRAW template
    if (m_file_name.Contains("RAYRAWTEMP")) {
      TString name("RAYRAW");
      auto &cont = m_tempfunc_collection[name];
      for (auto *f : cont) delete f;
      cont.clear();

      for (int i = 0; i < NumOfSegRayraw; ++i) {
        TString fname = Form("%s.%d", m_file_name.Data(), i);
        TemplateFitFunction *tempfunc = new TemplateFitFunction(fname);
        cont.push_back(tempfunc);
      }
    }
  }

  m_is_ready = true;
  return true;
}

const TempFuncContainer&
TemplateFitMan::GetHitContainer(const TString& name) const
{
  auto itr = m_tempfunc_collection.find(name);
  if (itr == m_tempfunc_collection.end()) {
    static TempFuncContainer null_container;
    return null_container;
  }
  return itr->second;
}

// -*- C++ -*-

#include "TPCRawHit.hh"

#include <iostream>
#include <iterator>

#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "FuncName.hh"
#include "TPCPadHelper.hh"

//_____________________________________________________________________________
TPCRawHit::TPCRawHit(Int_t layer, Int_t row, Double_t* pars)
  : m_layer_id(layer),
    m_row_id(row)
{
  if(pars){
    m_pars.push_back(pars[0]);
    m_pars.push_back(pars[1]);
    m_pars.push_back(pars[2]);
  }
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCRawHit::~TPCRawHit()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
TPCRawHit::AddFadc(Double_t adc)
{
  m_fadc.push_back(adc);
}

//_____________________________________________________________________________
Double_t
TPCRawHit::LocMax(TB_t tmin, TB_t tmax) const
{
  tmax = TMath::Min(tmax, m_fadc.size());
  if(tmin >= tmax) return TMath::QuietNaN();
  FADC_t sub {m_fadc.begin()+tmin, m_fadc.begin()+tmax};
  return TMath::LocMax(sub.size(), sub.data());
}

//_____________________________________________________________________________
Double_t
TPCRawHit::MaxAdc(TB_t tmin, TB_t tmax) const
{
  tmax = TMath::Min(tmax, m_fadc.size());
  if(tmin >= tmax) return TMath::QuietNaN();
  FADC_t sub {m_fadc.begin()+tmin, m_fadc.begin()+tmax};
  return TMath::MaxElement(sub.size(), sub.data());
}

//_____________________________________________________________________________
Double_t
TPCRawHit::Mean(TB_t tmin, TB_t tmax) const
{
  tmax = TMath::Min(tmax, m_fadc.size());
  if(tmin >= tmax) return TMath::QuietNaN();
  FADC_t sub {m_fadc.begin()+tmin, m_fadc.begin()+tmax};
  return TMath::Mean(sub.size(), sub.data());
}

//_____________________________________________________________________________
Double_t
TPCRawHit::MinAdc(TB_t tmin, TB_t tmax) const
{
  tmax = TMath::Min(tmax, m_fadc.size());
  if(tmin >= tmax) return TMath::QuietNaN();
  FADC_t sub {m_fadc.begin()+tmin, m_fadc.begin()+tmax};
  return TMath::MinElement(sub.size(), sub.data());
}

//_____________________________________________________________________________
void
TPCRawHit::Print(Option_t*) const
{
  hddaq::cout << FUNC_NAME << " " << std::endl
              << "   layer  = " << m_layer_id  << std::endl
              << "   row    = " << m_row_id << std::endl
              << "   fadc.size() = " << m_fadc.size() << std::endl
              << "   maxadc = " << MaxAdc(0, 170) << std::endl
              << "   locmax = " << LocMax(0, 170) << std::endl
              << "   mean   = " << Mean(0, 170) << std::endl
              << "   rms    = " << RMS(0, 170) << std::endl
              << std::endl;
}

//_____________________________________________________________________________
Double_t
TPCRawHit::RMS(TB_t tmin, TB_t tmax) const
{
  tmax = TMath::Min(tmax, m_fadc.size());
  if(tmin >= tmax) return TMath::QuietNaN();
  FADC_t sub {m_fadc.begin()+tmin, m_fadc.begin()+tmax};
  return TMath::RMS(sub.size(), sub.data());
}

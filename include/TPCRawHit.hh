// -*- C++ -*-

#ifndef TPC_RAW_HIT_HH
#define TPC_RAW_HIT_HH

#include <vector>

#include <TMath.h>
#include <TString.h>

#include "DetectorID.hh"

typedef std::vector<Double_t> FADC_t;
typedef FADC_t::size_type     TB_t; // ULong_t
typedef std::vector<Double_t> PRM_t;

//_____________________________________________________________________________
class TPCRawHit
{
public:
  static TString ClassName();
  TPCRawHit(Int_t layer, Int_t row, Double_t* pars=nullptr);
  ~TPCRawHit();

private:
  Int_t  m_layer_id;
  Int_t  m_row_id;
  FADC_t m_fadc;
  Double_t m_raw_rms=0;
  // baseline correction
  PRM_t m_pars; // p0:adc_ofs, p1:scale, p2:tb_ofs

public:
  void          AddFadc(Double_t adc);
  const FADC_t& Fadc() const { return m_fadc; }
  TB_t          FadcSize() const { return m_fadc.size(); }
  const PRM_t&  GetParameters() const { return m_pars; }
  Int_t         LayerId() const { return m_layer_id; }
  Double_t      LocMax(TB_t tmin, TB_t tmax) const;
  void          Print(Option_t* opt="") const;
  Int_t         RowId() const { return m_row_id; }
  Double_t      MaxAdc(TB_t tmin, TB_t tmax) const;
  Double_t      MinAdc(TB_t tmin, TB_t tmax) const;
  Double_t      Mean(TB_t tmin, TB_t tmax) const;
  Double_t      RMS(TB_t tmin, TB_t tmax) const;
  Double_t      RMS_10() const
  { return TMath::RMS(10, m_fadc.data()); }
  Double_t	RawRMS() const {return m_raw_rms;}
  void          SetRawRMS(Double_t raw_rms){ m_raw_rms = raw_rms;}
};

//_____________________________________________________________________________
inline TString
TPCRawHit::ClassName()
{
  static TString s_name("TPCRawHit");
  return s_name;
}

#endif

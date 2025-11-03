// -*- C++ -*-

#ifndef HODO_WAVEFORM_HIT_HH
#define HODO_WAVEFORM_HIT_HH

#include <string>

#include <std_ostream.hh>

#include "HodoHit.hh"
#include "TGraphErrors.h"
#include "TF1.h"

struct SearchParam{
  TString  Name;
  Int_t    tgen[2];
  Double_t sbegin;
  Double_t send;
  Double_t cutbegin;
  Double_t cutend;
  Double_t threshold;
  Double_t width;
  Double_t risetime;
  std::vector<Double_t> foundx;
  std::vector<Double_t> foundy;
};


typedef struct {
  std::string Name;
  Int_t tgen;
  Int_t color;
  Double_t FitStart;
  Double_t FitEnd;
  Int_t wavenum;
  Int_t ParaNum;
  Double_t par[64];
  Double_t FitParam[64];
  Double_t Residual;
} FitParam;


//_____________________________________________________________________________
class HodoWaveformHit : public HodoHit
{
public:
  static const TString& ClassName();
  explicit HodoWaveformHit(HodoRawHit* hit);
  virtual  ~HodoWaveformHit();

private:
  HodoWaveformHit();
  HodoWaveformHit(const HodoWaveformHit& rhit);
  HodoWaveformHit& operator =(const HodoWaveformHit& rhit);

protected:
  Double_t m_position;
  Double_t m_dxdw;
  Double_t m_adc_integral;
  Double_t m_n_discri_pulse;
  Double_t m_n_discri_diffpulse;

  using waveform_t = std::vector<std::vector<std::pair<Double_t, Double_t>>>; // (time, de)[ch][mhit]
  waveform_t m_waveform;

  data_t     m_pulse_height;
  data_t     m_pulse_time;

  using TGraphC = std::vector<TGraphErrors*>;
  TGraphC m_TGraphC;
  TF1    *m_func;

  mutable Bool_t   m_JoinTrack;

public:
  Bool_t   Calculate();
  Bool_t   PulseSearch();
  Bool_t   MakeGraph();
  Bool_t   MakeDifGraph(Int_t index_org);
  Bool_t   PreSearch(struct SearchParam *sp);
  Bool_t   WidthCut(std::vector<Double_t> rise,
		    std::vector<Double_t> fall,
		    Double_t width, std::vector<Double_t> &outrise);
  void     CompareRise(std::vector<Double_t> rise1,
		       std::vector<Double_t> rise2,
		       Double_t width, std::vector<Double_t> &outrise);
  void     SetInitial(std::vector<Double_t> &v,
		      Double_t begin, Double_t end,
		      Double_t thre, Double_t rise);
  Double_t GXtoGY(Int_t index_graph, Double_t gx);
  Bool_t   SetFitParam(FitParam *fp, std::vector<Double_t> &inix,
		       std::vector<Double_t> &iniy);
  void     Fit1(FitParam *fp);
  Double_t FittedTrigX(FitParam fp, Double_t allowance);
  Double_t RisingResidual(Int_t tge_No, Double_t trig, Double_t &res_max);

  Int_t GetWaveformEntries(Int_t i) const
    { return m_waveform.at(i).size(); }
  Int_t GetWaveformEntriesEntries() const
  { if(m_n_ch == HodoRawHit::kNChannel) return GetWaveformEntries(HodoRawHit::kExtra);
    else return GetWaveformEntries(HodoRawHit::kUp); }

  std::pair<Double_t, Double_t> GetWaveform(Int_t j) const
    { return m_waveform.at(HodoRawHit::kUp).at(j); }

  std::pair<Double_t, Double_t> GetWaveform(Int_t i, Int_t j) const
    { return m_waveform.at(i).at(j); }

  Int_t GetNPulse(Int_t i) const
  { return m_pulse_height.at(i).size(); }
  Int_t GetNPulse() const
  { if(m_n_ch == HodoRawHit::kNChannel) return GetNPulse(HodoRawHit::kExtra);
    else return GetNPulse(HodoRawHit::kUp); }

  Int_t GetNGraph() const {return m_TGraphC.size(); }
  inline TGraphErrors * GetTGraph( std::size_t i ) const;
  inline TF1 * GetFitTF1() const;

  Double_t GetAdcIntegral() const
  { return m_adc_integral; }

  Double_t GetNDiscriPulse() const
  { return m_n_discri_pulse; }

  Double_t GetNDiscriDiffPulse() const
  { return m_n_discri_diffpulse; }

  Double_t GetPulseHeight(Int_t i, Int_t j) const
  { return m_pulse_height.at(i).at(j); }

  Double_t GetPulseHeight(Int_t j) const
  { return m_pulse_height.at(HodoRawHit::kUp).at(j); }

  Double_t GetPulseTime(Int_t i, Int_t j) const
  { return m_pulse_time.at(i).at(j); }

  Double_t GetPulseTime(Int_t j) const
  { return m_pulse_time.at(HodoRawHit::kUp).at(j); }

  void   SetTrackJoined() const { m_JoinTrack = true; }
  Bool_t IsTrackJoined() const  { return m_JoinTrack; }

  virtual void   Print(Option_t* arg="") const;
  virtual Bool_t ReCalc(Bool_t allpyRecursively=false){ return Calculate(); }

  //static Bool_t Compare(const HodoWaveformHit* left, const HodoWaveformHit* right);
};

//_____________________________________________________________________________
inline const TString&
HodoWaveformHit::ClassName()
{
  static TString s_name("HodoWaveformHit");
  return s_name;
}

//______________________________________________________________________________
inline TGraphErrors*
HodoWaveformHit::GetTGraph( std::size_t i ) const
{
  if( i<m_TGraphC.size() )
    return m_TGraphC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline TF1*
HodoWaveformHit::GetFitTF1() const
{
  return m_func;
}
//_____________________________________________________________________________
/*
inline Bool_t
HodoWaveformHit::Compare(const HodoWaveformHit* left, const HodoWaveformHit* right)
{
  return left->Position() < right->Position();
}
*/
#endif

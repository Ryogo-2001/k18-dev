#ifndef RAYRAWWAVEFORMHIT_HH
#define RAYRAWWAVEFORMHIT_HH 1

#include "HodoHit.hh"
#include "HodoRawHit.hh"
//#include "HodoAnalyzer.hh"

#include <vector>
#include <map>

#include <TGraphErrors.h>
#include <TF1.h>

struct FitParam {
  TString FitName;
  Int_t   tgen;
  Int_t   color;
  Double_t FitStart;
  Double_t FitEnd;
  Int_t    wavenum;
  Int_t    ParaNum;
  Double_t par[64]; // (ParaMax と同じサイズ)
  Double_t FitParam[64];
  Double_t Residual;
};

struct SearchParam {
  TString SearchName;
  Int_t tgen[2];
  Double_t sbegin;
  Double_t send;
  Double_t fbegin;
  Double_t fend;
  Double_t threshold;
  Double_t width;
  Double_t risetime;
  std::vector<Double_t> foundx;
  std::vector<Double_t> foundy;
};

class RayrawWaveformHit : public HodoHit
{
public:
  RayrawWaveformHit( HodoRawHit *rhit );
  virtual ~RayrawWaveformHit();

  virtual bool Calculate();
  virtual void Print(Option_t* arg="") const;
  const char* ClassName() const { return "RayrawWaveformHit"; }


  
  Int_t GetWaveformEntries(Int_t ch) const
  { return m_waveform.at(ch).size(); }
  std::pair<Double_t, Double_t> GetWaveform(Int_t ch, Int_t i) const
  { return m_waveform.at(ch).at(i); }
  
  Int_t GetNPulse(Int_t ch) const { return m_pulse_height.at(ch).size(); }
  Double_t GetPulseHeight(Int_t ch, Int_t i) const { return m_pulse_height.at(ch).at(i); }
  Double_t GetPulseTime(Int_t ch, Int_t i) const { return m_pulse_time.at(ch).at(i); }

protected:

  Bool_t MakeGraph();
  Bool_t MakeDifGraph(Int_t index_org);
  Bool_t PulseSearch( void );
  Bool_t PreSearch(struct SearchParam *sp);
  Bool_t SetFitParam(FitParam *fp, std::vector<Double_t> &inix,
         std::vector<Double_t> &iniy);
  void Fit1(FitParam *fp);
  Bool_t WidthCut(std::vector<Double_t> rise,
      std::vector<Double_t> fall,
      Double_t width, std::vector<Double_t> &outrise);
  void CompareRise(std::vector<Double_t> rise1,
       std::vector<Double_t> rise2,
       Double_t width, std::vector<Double_t> &outrise);
  void SetInitial(std::vector<Double_t> &v,
      Double_t begin, Double_t end,
      Double_t thre, Double_t rise);
  Double_t GXtoGY(Int_t index_graph, Double_t gx);
  Double_t FittedTrigX(FitParam fp, Double_t allowance);
  Double_t RisingResidual(Int_t tge_No, Double_t trig, Double_t &res_max);

  
  std::vector<std::vector<std::pair<Double_t, Double_t>>> m_waveform;
  std::vector<std::vector<Double_t>> m_pulse_height;
  std::vector<std::vector<Double_t>> m_pulse_time;
  Double_t m_position;
  Double_t m_adc_integral;
  Double_t m_n_discri_pulse;
  Double_t m_n_discri_diffpulse;
  TF1 *m_func;
  std::vector<TGraphErrors*> m_TGraphC;

  TF1* GetFittedFunc() const { return m_func; } 

  Bool_t m_JoinTrack;
};

#endif
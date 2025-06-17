// -*- C++ -*-

#include "TPCHit.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>

#include <TMath.h>
#include <TSpectrum.h>

#include <std_ostream.hh>

#include "DeleteUtility.hh"
#include "FuncName.hh"
#include "DebugCounter.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCLTrackHit.hh"
#include "TPCPadHelper.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "TPCRawHit.hh"
#include "UserParamMan.hh"
#include "ConfMan.hh"

//#define QuickAnalysis  1 // User EventSelectionTPCHits in RawData.cc
//#define FitPedestal    1
#define FitPedestal    0
//#define UseGaussian    1 //for gate noise fit
#define SingleFit      0
#define UseGaussian    0
#define UseLandau      0
#define UseExponential 0
#define UseGumbel      0
#define MultiFit       1
#define UseMultiGumbel 1
#define DebugEvDisp    0

#if DebugEvDisp
#include <TApplication.h>
#include <TCanvas.h>
#include <TSystem.h>
namespace { TApplication app("DebugApp", nullptr, nullptr); }
#endif

namespace
{
  const auto& gUser   = UserParamMan::GetInstance();
  const auto& gTPC    = TPCParamMan::GetInstance();
  const auto& gTPCPos = TPCPositionCorrector::GetInstance();
  const Int_t MaxADC = 4096;
  const Int_t MaxIteration = 3;
  const Int_t MaxPeaks = 20;
  const Double_t MaxChisqr = 1000.;
  const Double_t& HSfield_Hall = ConfMan::Get<Double_t>("HSFLDHALL");
  const auto& ResParamInnerLayerHSOn = gTPC.TPCResolutionParams(true, false); //B=1 T, Inner layers
  const auto& ResParamOuterLayerHSOn = gTPC.TPCResolutionParams(true, true); //B=1 T, Outer layers
  const auto& ResParamInnerLayerHSOff = gTPC.TPCResolutionParams(false, false); //B=0, Inner layers
  const auto& ResParamOuterLayerHSOff = gTPC.TPCResolutionParams(false, true); //B=0, Outer layers
}

//_____________________________________________________________________________
TPCHit::TPCHit(TPCRawHit* rhit)
  : DCHit(rhit->LayerId(), rhit->RowId()),
    m_rhit(rhit),
    m_layer(rhit->LayerId()),
    m_row(rhit->RowId()),
    m_padtheta(tpc::getTheta(m_layer, m_row)*TMath::DegToRad()),
    m_padlength(tpc::padParameter[m_layer][5]),
    m_mrow(TMath::Nint(m_row)),
    m_pad(tpc::GetPadId(m_layer, m_row)),
    m_pedestal(TMath::QuietNaN()),
    m_rms(TMath::QuietNaN()),
    m_raw_rms(TMath::QuietNaN()),
    m_de(),
    m_time(), // [time bucket]
    m_chisqr(),
    m_cde(),
    m_ctime(), // [ns]
    m_drift_length(),
    m_is_good(false),
    m_is_calculated(false),
    m_hough_flag(),
    m_houghY_num(),
    m_hough_dist(),
    m_hough_disty(),
    m_res_param()
{
  debug::ObjectCounter::increase(ClassName());

  if(HSfield_Hall<0.1&&m_layer<10) m_res_param = ResParamInnerLayerHSOff;
  else if(HSfield_Hall<0.1&&m_layer>=10) m_res_param = ResParamOuterLayerHSOff;
  else if(m_layer<10) m_res_param = ResParamInnerLayerHSOn;
  else m_res_param = ResParamOuterLayerHSOn;
}

//_____________________________________________________________________________
TPCHit::TPCHit(Int_t layer, Double_t mrow)
  : DCHit(layer, mrow),
    m_rhit(),
    m_layer(layer),
    m_row(TMath::Nint(mrow)),
    m_padtheta(tpc::getTheta(m_layer, mrow)*TMath::DegToRad()),
    m_padlength(tpc::padParameter[m_layer][5]),
    m_mrow(mrow),
    m_pad(tpc::GetPadId(layer, mrow)),
    m_pedestal(TMath::QuietNaN()),
    m_rms(TMath::QuietNaN()),
    m_de(),
    m_time(), // [time bucket]
    m_chisqr(),
    m_cde(),
    m_ctime(), // [ns]
    m_drift_length(),
    m_is_good(true),
    m_is_calculated(false),
    m_hough_flag(),
    m_houghY_num(),
    m_hough_dist(),
    m_hough_disty(),
    m_res_param()
{

  debug::ObjectCounter::increase(ClassName());

  if(HSfield_Hall<0.1&&m_layer<10) m_res_param = ResParamInnerLayerHSOff;
  else if(HSfield_Hall<0.1&&m_layer>=10) m_res_param = ResParamOuterLayerHSOff;
  else if(m_layer<10) m_res_param = ResParamInnerLayerHSOn;
  else m_res_param = ResParamOuterLayerHSOn;
}


//_____________________________________________________________________________
TPCHit::~TPCHit()
{
  ClearRegisteredHits();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
TPCHit::AddHit(Double_t de, Double_t time, Double_t sigma, Double_t chisqr)
{
  m_time.push_back(time);
  m_de.push_back(de);
  m_sigma.push_back(sigma);
  m_chisqr.push_back(chisqr);
  m_cde.push_back(TMath::QuietNaN());
  m_ctime.push_back(TMath::QuietNaN());
  m_drift_length.push_back(TMath::QuietNaN());
  m_position.push_back(TVector3());
}

//_____________________________________________________________________________
void
TPCHit::SetHoughYnum(Int_t houghY_num)
{
  m_houghY_num.push_back(houghY_num);
}

//_____________________________________________________________________________
Bool_t
TPCHit::Calculate(Double_t clock)
{
  if(m_de.size() != m_time.size()){
    hddaq::cerr << FUNC_NAME << "found size mismatch: "
                << "m_de.size()=" << m_de.size() << ", "
                << "m_time.size()=" << m_time.size() << std::endl;
    return false;
  }

  if(!gTPC.IsReady()){
    hddaq::cerr << FUNC_NAME << " TPCParamMan must be initialized" << std::endl;
    return false;
  }

  for(Int_t i=0, n=m_de.size(); i<n; ++i){
    Double_t cde, ctime, dl;
    if(!gTPC.GetCDe(m_layer, m_row, m_de[i], cde)){
      hddaq::cerr << FUNC_NAME << " something is wrong at GetCDe("
                  << m_layer << ", " << m_row << ", " << m_de[i]
                  << ", " << cde << ")" << std::endl;
    }
    if(!gTPC.GetCTime(m_layer, m_row, m_time[i], ctime)){
      hddaq::cerr << FUNC_NAME << " something is wrong at GetCTime("
                  << m_layer << ", " << m_row << ", " << m_time[i]
                  << ", " << ctime << ")" << std::endl;
    }

    Double_t cclk;
    if(!gTPC.GetCClock(m_layer, m_row, clock, cclk)){
      hddaq::cerr << FUNC_NAME << " something is wrong at GetCClock("
                  << m_layer << ", " << m_row << ", " << m_time[i]
                  << ", " << ctime << ")" << std::endl;
    }
    ctime += cclk;

    if(!gTPC.GetDriftLength(m_layer, m_row, ctime, dl)){
      hddaq::cerr << FUNC_NAME << " something is wrong at GetDriftLength("
                  << m_layer << ", " << m_row << ", " << ctime
                  << ", " << dl << ")" << std::endl;
    }
    m_cde[i] = cde;
    m_ctime[i] = ctime;
    m_drift_length[i] = dl;
    auto pos = tpc::getPosition(m_pad);
    pos.SetY(dl);
    auto cpos = gTPCPos.Correct(pos, m_layer, m_row);
    m_position[i] = cpos;
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCHit::DoFit()
{
  static const Double_t MinDe = gUser.GetParameter("MinDeTPC");
  static const Double_t MinRms = gUser.GetParameter("MinRmsTPC");
  static const Double_t MinRawRms = gUser.GetParameter("MinBaseRmsTPC");
  static const Int_t MinTimeBucket = gUser.GetParameter("TimeBucketTPC", 0);
  static const Int_t MaxTimeBucket = gUser.GetParameter("TimeBucketTPC", 1);
  static const Int_t NumOfTimeBucket = gUser.GetParameter("NumOfTimeBucket");

  if(m_is_calculated){
    hddaq::cerr << FUNC_NAME << " already calculated" << std::endl;
    return false;
  }

  if(!m_rhit){
    hddaq::cerr << FUNC_NAME << " m_rhit is nullptr" << std::endl;
    return false;
  }

#if DebugEvDisp
  //  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetOptStat(0);
  //  Option_t* option = "";
  Option_t* option = "W";//to avoid overflow events
#else
  //  Option_t* option = "Q";
  Option_t* option = "Q+W";//to avoid overflow events
#endif

  const auto level = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kError;

  {
    Double_t mean = m_rhit->Mean(0, NumOfTimeBucket);
    Double_t rms = m_rhit->RMS(0, NumOfTimeBucket);
    Double_t rawrms = m_rhit->RawRMS();
    m_pedestal = mean;
    m_rms = rms;
    m_raw_rms = rawrms;
  }
  // #if QuickAnalysis
  // {
  //   Double_t max_adc = m_rhit->MaxAdc(0, NumOfTimeBucket);
  //   //std::cout<<"max_dE = "<<max_adc - m_pedestal<<std::endl;
  //   if(max_adc-m_pedestal < MinDe){
  //     return false;
  //   }
  // }
  // #endif

#if 0 // for check
  if(m_rhit->MaxAdc(0, NumOfTimeBucket)-m_pedestal < 400
     || m_rhit->LocMax(0, NumOfTimeBucket) < 40
     || m_rhit->LocMax(0, NumOfTimeBucket) > 140){
    return false;
  }
#endif

  m_layer = m_rhit->LayerId();
  m_row = m_rhit->RowId();
  m_pad = tpc::GetPadId(m_layer, m_row);
  static TCanvas c1("c"+FUNC_NAME, FUNC_NAME, 800, 600);
  c1.cd();
  TH1D h1(FUNC_NAME+"-h1", "Pedestal",
          MaxADC, 0, MaxADC);
  TH1D h2(FUNC_NAME+"-h2",
          Form("Layer#%d Row#%d;Time Bucket;ADC",
               m_layer, m_row),
          NumOfTimeBucket, 0, NumOfTimeBucket);
  for(Int_t i=0, n=m_rhit->Fadc().size(); i<n; ++i){
    Double_t tb = i + 1;
    Double_t adc = m_rhit->Fadc().at(i);
    if(i < MinTimeBucket || MaxTimeBucket < i)
      continue;
    h1.Fill(adc);
    // if(adc < 4075){
    h2.SetBinContent(tb, adc);
    // h2.SetBinError(tb, 1);
    // }
  }
  h2.GetXaxis()->SetRangeUser(MinTimeBucket, MaxTimeBucket);
  h2.GetYaxis()->SetRangeUser(m_pedestal-100., m_rhit->MaxAdc(0, NumOfTimeBucket)+200);
  // h2.GetYaxis()->SetRangeUser(0, 0x1000);
  // h2.GetYaxis()->SetRangeUser(0, TMath::Max(m_rhit->MaxAdc(0, NumOfTimeBucket)+200, 2000.));
  //std::cout<<"pedestal:"<<m_pedestal<<", max_adc:"<<m_rhit->MaxAdc(0, NumOfTimeBucket)<<std::endl;

#if FitPedestal
  { //___ Pedestal Fitting
    Double_t mean = h1.GetXaxis()->GetBinCenter(h1.GetMaximumBin());
    Double_t rms = m_rhit->RMS(0, NumOfTimeBucket);
    if(mean > 1000.) mean = 400.;
    TF1 f1("f1", "gaus", 0, MaxADC);
    f1.SetParameter(0, NumOfTimeBucket/rms);
    f1.SetParLimits(0, 1, NumOfTimeBucket);
    f1.SetParameter(1, mean);
    f1.SetParLimits(1, mean-3*rms, mean+3*rms);
    f1.SetParameter(2, rms);
    f1.SetParLimits(2, 0.2*rms, 1.*rms);
    for(UInt_t i=0; i<MaxIteration; ++i){
      h1.Fit("f1", option, "", mean - 2*rms, mean + 2*rms);
      mean = f1.GetParameter(1);
      rms = f1.GetParameter(2);
    }
    Double_t chisqr = f1.GetChisquare() / f1.GetNDF();
    if(chisqr < 100.){
      m_pedestal = mean;
      // m_rms = rms;
    } else {
      // hddaq::cerr << FUNC_NAME << " failed to fit pedestal" << std::endl;
      // Print();
    }
#if DebugEvDisp
    hddaq::cout << FUNC_NAME << " (mean,rms,chisqr)=("
		<< mean << "," << rms << "," << chisqr << ")" << std::endl;
    h1.GetXaxis()->SetRangeUser(mean-10*rms, mean+10*rms);
    h1.Draw();
    c1.Modified();
    c1.Update();
    // getchar();
    Double_t max_adc = m_rhit->MaxAdc(0, NumOfTimeBucket);
    if(max_adc - mean < MinDe){
      return false;
    }
    if(m_raw_rms < MinRawRms){
      return false;
    }
#endif
  }
#endif

  if(m_rms < MinRms){
    return false;
  }


#if 1
  { //___ Peak Search
    Double_t sigma = 3;
    Double_t threshold = 0.4;
    TSpectrum spec(MaxPeaks);
    const Int_t n_peaks = spec.Search(&h2, sigma, "", threshold);
    if(n_peaks == 0)
      return false;
    Double_t* x_peaks = spec.GetPositionX();
    std::vector<Double_t> sorted_x_peaks(n_peaks);
    sorted_x_peaks.assign(x_peaks, x_peaks+n_peaks);
    std::stable_sort(sorted_x_peaks.begin(), sorted_x_peaks.end());

#if SingleFit
    for(Int_t i=0; i<n_peaks; ++i){
      Double_t time = sorted_x_peaks[i];
      Double_t de = h2.GetBinContent(h2.FindBin(time)) - m_pedestal;
      // hddaq::cout << FUNC_NAME << " time,de=" << time << "," << de << std::endl;
#if UseGaussian
      TF1 f1("f1", "gaus(0)+pol1(3)", 0, NumOfTimeBucket);
      sigma = 3.;
      f1.SetParameter(0, de);
      f1.SetParLimits(0, de*0.8, de*1.2);
      f1.SetParameter(1, time);
      f1.SetParLimits(1, time-2*sigma, time+2*sigma);
      f1.SetParameter(2, sigma);
      f1.SetParLimits(2, 2., 5.);
      // f1.FixParameter(2, 1.0);
      f1.SetParameter(3, m_pedestal);
      f1.SetParLimits(3, m_pedestal-3*m_rms, m_pedestal+3*m_rms);
      f1.FixParameter(4, 0);
      // f1.SetParLimits(4, -0.01, 0.01);
      for(UInt_t j=0; j<MaxIteration; ++j){
        h2.Fit("f1", option, "", time-2.5*sigma, time+2.5*sigma);
        time = f1.GetParameter(1);
        sigma = f1.GetParameter(2);
      }
      auto p = f1.GetParameters();
      time = f1.GetMaximumX();
      de = f1.GetMaximum() - p[3];
      // de = p[0]*p[2]*TMath::Sqrt(2*TMath::Pi());
      Double_t chisqr = f1.GetChisquare() / f1.GetNDF();
#elif UseLandau
      TF1 f1("f1", "landau(0)+pol1(3)", 0, NumOfTimeBucket);
      sigma = 1.;
      f1.SetParameter(0, de/0.18);
      f1.SetParLimits(0, de/0.18*0.8, de/0.18*1.2);
      f1.SetParameter(1, time);
      f1.SetParLimits(1, time-2*sigma, time+2*sigma);
      f1.SetParameter(2, sigma);
      f1.SetParLimits(2, 0.5, 1.5);
      // f1.FixParameter(2, 1.0);
      f1.SetParameter(3, m_pedestal);
      f1.SetParLimits(3, m_pedestal-3*m_rms, m_pedestal+3*m_rms);
      f1.FixParameter(4, 0);
      // f1.SetParLimits(4, -0.01, 0.01);
      for(UInt_t j=0; j<MaxIteration; ++j){
        h2.Fit("f1", option, "", time-3*sigma, time+7*sigma);
        time = f1.GetParameter(1);
        sigma = f1.GetParameter(2);
      }
      auto p = f1.GetParameters();
      time = f1.GetMaximumX();
      de = f1.GetMaximum() - p[3]; //p[0]*p[2]*TMath::Sqrt(2*TMath::Pi());
      Double_t chisqr = f1.GetChisquare() / f1.GetNDF();
#elif UseExponential
      TF1 f1("f1", "[0]*(x-[1])*exp(-(x-[1])/[2])+[3]", 0, NumOfTimeBucket);
      Double_t c = 3.; // decay constant
      f1.SetParameter(0, de*TMath::Exp(1)/c);
      f1.SetParLimits(0, TMath::Exp(1)/c, m_rhit->MaxAdc(0, NumOfTimeBucket)*TMath::Exp(1)/c);
      f1.SetParameter(1, time+c);
      f1.SetParLimits(1, time+c-12, time+c+12);
      f1.SetParameter(2, c);
      f1.SetParLimits(2, 1.5, 7.5);
      f1.FixParameter(3, m_pedestal);
      // f1.SetParameter(3, m_pedestal);
      // f1.SetParLimits(3, m_pedestal - m_rms, m_pedestal + m_rms);
      h2.Fit("f1", option, "", time-5, time+10);
      auto p = f1.GetParameters();
      time = p[1] + p[2]; // peak time
      // c = p[2];
      de = p[0]*p[2]*TMath::Exp(-1); // amplitude
      sigma = p[2];
      // de = p[0]*p[2]*p[2]; // integral
      Double_t chisqr = f1.GetChisquare() / f1.GetNDF();
#elif UseGumbel
      TF1 f1("f1", "[0]*exp(-(x-[1])/[2])*exp(-exp(-(x-[1])/[2]))+[3]",
             0, NumOfTimeBucket);
      f1.SetParameter(0, de);
      //f1.SetParLimits(0, 0, m_rhit->MaxAdc(0, NumOfTimeBucket)+m_pedestal+m_rms);
      //f1.SetParLimits(0, 0, de+m_rms);
      f1.SetParameter(1, time);
      f1.SetParLimits(1, time-10, time+10);
      f1.SetParameter(2, 3);
      f1.SetParLimits(2, 1., 10.);
      f1.SetParameter(3, m_pedestal);
      f1.SetParLimits(3, m_pedestal - 0.5*m_rms, m_pedestal + 0.5*m_rms);
      //h2.Fit("f1", option, "", time-6, time+10);
      //h2.Fit("f1", option, "",MinTimeBucket , MaxTimeBucket);
      h2.Fit("f1", option, "", time-6.-de*0.07, time+10.+de*0.11);
      auto p = f1.GetParameters();
      time = p[1]; // peak time
      de = f1.Eval(time) - p[3]; // amplitude
      sigma = p[2];
      Double_t chisqr = f1.GetChisquare()/f1.GetNDF();
      //std::cout<<"chisqr:"<<chisqr<<std::endl;
#endif
      if (chisqr > MaxChisqr || de < MinDe)
        continue;
      AddHit(de, time, sigma, chisqr);
#if DebugEvDisp
      hddaq::cout << FUNC_NAME << " (time,de,chisqr, rms)=("
		  << time << "," << de << "," << chisqr << "," << m_rms <<")" << std::endl;
#endif
    }
#endif

#if MultiFit
#if UseMultiGumbel
    Double_t time[n_peaks];
    Double_t de[n_peaks];
    TString func;
    for(Int_t i=0; i<n_peaks; ++i){
      time[i] = sorted_x_peaks[i];
      de[i] = h2.GetBinContent(h2.FindBin(time[i])) - m_pedestal;
      if(i!=0)
	func += "+";
      func += "["+std::to_string(i*3)+"]";
      func += "*exp(-(x-["+std::to_string(i*3+1)+"])";
      func += "/["+std::to_string(i*3+2)+"])";
      func += "*exp(-exp(-(x-["+std::to_string(i*3+1)+"])/["+std::to_string(i*3+2)+"]))";
    }
    func += "+["+std::to_string(n_peaks*3)+"]";
    TF1 f1("f1", func,
	   0, NumOfTimeBucket);
    for(Int_t i=0; i<n_peaks; ++i){
      f1.SetParameter(i*3, de[i]);
      f1.SetParameter(i*3+1, time[i]);
      f1.SetParLimits(i*3+1, time[i]-10, time[i]+10);
      f1.SetParameter(i*3+2, 3);
      f1.SetParLimits(i*3+2, 1., 10.);
      //f1.SetParLimits(i*3+2, 2.5, 10.);
    }
    f1.SetParameter(n_peaks*3, m_pedestal);
    f1.SetParLimits(n_peaks*3, m_pedestal-0.5*m_rms, m_pedestal+0.5*m_rms);
    Double_t FitRangeMin = MinTimeBucket;
    Double_t FitRangeMax = MaxTimeBucket;
    if(MinTimeBucket<time[0]-6.-de[0]*0.07)
      FitRangeMin = time[0]-6.-de[0]*0.07;
    if(MaxTimeBucket>time[n_peaks-1]+10.+de[n_peaks-1]*0.11)
      FitRangeMax = time[n_peaks-1]+10.+de[n_peaks-1]*0.11;
    h2.Fit("f1", option, "", FitRangeMin, FitRangeMax);
    auto p = f1.GetParameters();
    Double_t chisqr = f1.GetChisquare()/f1.GetNDF();
    if(chisqr < MaxChisqr){
      for(Int_t i=0; i<n_peaks; ++i){
	time[i] = p[i*3+1]; // peak time
	de[i] = f1.Eval(time[i]) - p[n_peaks*3]; // amplitude
	sigma = p[i*3+2];
	if(de[i] > MinDe && sigma > 2.5){
          AddHit(de[i], time[i], sigma, chisqr);
	}
      }
    }
#endif
#endif
  }
#endif

#if DebugEvDisp
  //h2.Draw("L");
  {
    Double_t maxde = 0.;
    if(m_de.size()>0)
      maxde = TMath::MaxElement(m_de.size(), m_de.data());
    if(m_time.size()>0 //&& maxde>200.
       || true
       ){
      //h2.Draw("");
      c1.Modified();
      c1.Update();
      gSystem->ProcessEvents();
      Print();
      c1.Print("c1.pdf");
      getchar();
    }
  }
#endif

  gErrorIgnoreLevel = level;
  m_is_good = true;
  return true;
}

//_____________________________________________________________________________
void
TPCHit::ClearRegisteredHits()
{
  del::ClearContainer(m_register_container);
}

//_____________________________________________________________________________
Bool_t
TPCHit::IsGood() const
{
  return m_is_good;
}

//_____________________________________________________________________________
void
TPCHit::Print(const std::string& arg, std::ostream& ost) const
{
  const int w = 16;
  ost << "#D " << FUNC_NAME << " " << arg << std::endl
      << std::setw(w) << std::left << "layer" << m_layer << std::endl
      << std::setw(w) << std::left << "row"  << m_row  << std::endl
      << std::setw(w) << std::left << "pad"  << m_pad  << std::endl;
  for(Int_t i=0, n=GetNHits(); i<n; ++i){
    ost << std::setw(3) << std::right << i << " "
	<< "(time, de, chisqr, ctime, cde, dl, pos)=("
        << std::fixed << std::setprecision(3) << std::setw(9) << m_time[i]
        << std::fixed << std::setprecision(3) << std::setw(9) << m_de[i]
        << std::fixed << std::setprecision(3) << std::setw(9) << m_chisqr[i]
        << std::fixed << std::setprecision(3) << std::setw(9) << m_ctime[i]
        << std::fixed << std::setprecision(3) << std::setw(9) << m_cde[i]
        << std::fixed << std::setprecision(3) << std::setw(9) << m_drift_length[i]
        << std::fixed << std::setprecision(3) << std::setw(9) << m_position[i]
        << ")" << std::endl;
  }
}

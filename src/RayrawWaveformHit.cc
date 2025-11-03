// -*- C++ -*-

#include "RayrawWaveformHit.hh" 
#include "UserParamMan.hh"
#include "UnpackerManager.hh" 

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <limits>
#include <utility>
#include <algorithm> // std::find, std::sort, std::greater

#include <std_ostream.hh>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "FuncName.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "PrintHelper.hh"
#include "RawData.hh"
#include "BGODiscriminator.hh"
#include "TemplateFitMan.hh"
// #include "BGOCalibMan.hh" // BGO
// #include "RayrawCalibMan.hh" // RAYRAW


typedef std::vector<std::vector<Double_t>> data_t;

namespace
{
const auto qnan = TMath::QuietNaN();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gHodo = HodoParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();
const auto& gTempFit = TemplateFitMan::GetInstance();
//show debug info
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance(); 

// param
const double graphStart = 0.0;    // [ns]
const double graphEnd   = 70.0;   // [ns]
const double y_err      = 1.0;    // FADC's error 
const double fitStart = 0.0;    // [ns]
const double fitEnd   = 70.0;   // [ns]
const int ParaMax = 64;
const double TrigTimeReso = 1.00; // [ns], timing resolution for fitting

// ped
const int PEDESTAL_BIN_START = 0;
const int PEDESTAL_BIN_END   = 10;
}

//_____________________________________________________________________________
RayrawWaveformHit::RayrawWaveformHit(HodoRawHit *rhit)
  : HodoHit(rhit),
    m_waveform(m_n_ch),
    m_pulse_height(m_n_ch),
    m_pulse_time(m_n_ch),
    m_position(qnan),
    m_adc_integral(qnan),
    m_n_discri_pulse(qnan),
    m_n_discri_diffpulse(qnan),
    m_func(0),
    m_JoinTrack(false)
{
  debug::ObjectCounter::increase(ClassName());

  
  if (DetectorName()=="RAYRAW") {
    auto &cont = gTempFit.GetHitContainer(DetectorName());
    if (cont.empty() || SegmentId() >= (int)cont.size()) { 
        std::cerr << "RayrawWaveformHit Error: Template for "
                  << DetectorName() << " seg " << SegmentId() << " not found!" << std::endl;
        return;
    }
    TemplateFitFunction *tempFunc = cont.at(SegmentId());
    m_func = new TF1(DetectorName()+"-"+std::to_string(SegmentId()),
         tempFunc,
         fitStart, fitEnd, ParaMax );
  }
}

//_____________________________________________________________________________
RayrawWaveformHit::~RayrawWaveformHit()
{
  del::ClearContainer( m_TGraphC );

  if (m_func)
    delete m_func;

  m_func = 0;
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
bool
RayrawWaveformHit::Calculate()
{
  if(!HodoHit::Calculate())
    ;

  // gUser and SAMPLING_INTERVAL_NS
  static const auto& gUser = UserParamMan::GetInstance();
  static const auto SAMPLING_INTERVAL_NS = gUser.GetParameter("SamplingInterval");

  m_is_clustered.clear();
  for(Int_t ch=0; ch<m_n_ch; ++ch){
    m_de_high.at(ch).clear();
  }
  
  Int_t seg   = m_raw->SegmentId();
  
  
  Int_t ch = HodoRawHit::kUp; // (ch=0)

  //  (GetArrayAdc(ch) -> GetArrayAdc()) 
  const auto& fadc_samples = m_raw->GetArrayAdc(); 
  if (fadc_samples.empty() || fadc_samples.size() <= PEDESTAL_BIN_END) {
      return true; 
  }

  // calculate pedestal
  double pedestal = 0;
  int ped_counts = 0;
  for (int i = PEDESTAL_BIN_START; i < PEDESTAL_BIN_END; ++i) {
    pedestal += fadc_samples[i];
    ped_counts++;
  }
  if (ped_counts > 0) pedestal /= ped_counts;
  else {
      return true; 
  }

  
  Int_t ns = 0; 
  for(const auto& adc_raw : fadc_samples){
    Double_t time = (double)ns * SAMPLING_INTERVAL_NS;
    Double_t de = (Double_t)adc_raw - pedestal; 

    std::pair<Double_t, Double_t> wf_pair(time, de);
    m_waveform.at(ch).push_back(wf_pair);
    ns++;
  }
  /*
  // debug info
  if (seg == 0 && !m_waveform.at(ch).empty()) {
      std::cout << "!!! DEBUG [Calculate]  : (Evt " << gUnpacker.get_event_number() << ", Seg 0) Waveform filled. Size = " << m_waveform.at(ch).size() << std::endl;
  }
  */
  
  Int_t Nped=0;
  Double_t event_pedestal=0;
  for(const auto& wf: m_waveform.at(ch)){
    Double_t time = wf.first;
    if (time >= 0.0 && time < 5.0) { // 0.5ns * 10bin
      event_pedestal += wf.second; 
      Nped++;
    }
  }
  if (Nped > 0) event_pedestal /= Nped;
  else event_pedestal = 0;

  // (ADC integral)
  m_adc_integral = 0.;
  for(const auto& wf: m_waveform.at(ch)){
    Double_t time = wf.first;
    if (time >= fitStart && time <= fitEnd) {
      m_adc_integral += (wf.second - event_pedestal); 
    }
  }

  PulseSearch();
  m_is_calculated = true;
  return true;
}

//_____________________________________________________________________________
Bool_t RayrawWaveformHit::PulseSearch( void )
{
  if (!MakeGraph())
    return false;

  Int_t index_original_graph = 0;
  MakeDifGraph(index_original_graph);
  Int_t index_diff_graph = 1;
 
  Double_t threshold = 10.0;  
  Double_t width     = 0.05;   
  Double_t risetime  = 0.01;   
  
  SearchParam sp1={"sp1", {index_original_graph, index_diff_graph},
    fitStart, fitEnd, fitStart, fitEnd,
    threshold, width, risetime};

  Bool_t flagPresearch = PreSearch(&sp1);
  if (!flagPresearch)
    return false;

  Int_t color = 4;
  FitParam fp1={"fp1", index_original_graph, color ,fitStart, fitEnd};

  Bool_t flagSetFitParam = SetFitParam(&fp1,sp1.foundx,sp1.foundy);
  if (!flagSetFitParam)
    return false;

  Fit1(&fp1);

  Double_t trigx = FittedTrigX(fp1,1.0);
  
  if(std::abs(trigx)<TrigTimeReso){
    Double_t max_res = 0;
    fp1.Residual=  RisingResidual(index_original_graph, trigx, max_res);
    
    if(fp1.Residual < 50 || fabs(max_res) < 50){ 
      Int_t waveNum = fp1.wavenum;
      for (Int_t nw=0; nw<waveNum; nw++) {
        Double_t time = fp1.FitParam[2*nw+2];
        Double_t height = fp1.FitParam[2*nw+3];
        m_pulse_time.at(HodoRawHit::kUp).push_back(time);
        m_pulse_height.at(HodoRawHit::kUp).push_back(height);
        Double_t energy=-999.;
        
        if (GetName() == "RAYRAW") {
          energy = height; 
          m_de_high.at(HodoRawHit::kUp).push_back(energy);
        }
      }
    }
  }
  return true;
}

//_____________________________________________________________________________
Bool_t RayrawWaveformHit::MakeGraph()
{
  Int_t nc = GetWaveformEntries(HodoRawHit::kUp);
  Int_t n_range = 0;
  for (Int_t i=0; i<nc; i++) {
    std::pair<Double_t, Double_t> fadc = GetWaveform(HodoRawHit::kUp, i);
    if ( fadc.first >= graphStart && fadc.first <= graphEnd )
      n_range++;
  }

  /*
  // debug info
  if (SegmentId() == 0 && nc > 0) {
      std::cout << "!!! DEBUG [MakeGraph]    : (Seg 0, Evt " << gUnpacker.get_event_number() << ") Found " << nc << " raw samples." << std::endl;
      std::cout << "!!! DEBUG [MakeGraph]    : (Seg 0, Evt " << gUnpacker.get_event_number() << ") Time range cut (" << graphStart << " ns to " << graphEnd << " ns) passed " << n_range << " samples." << std::endl;
      if (n_range == 0) {
          std::cout << "!!! ERROR: MakeGraph cut ALL samples! Check graphStart/End parameters in namespace." << std::endl;
      }
  }
  */


  if (n_range<=0) {
    return false;
  }
  
  TGraphErrors *gr = new TGraphErrors(n_range);
  Int_t index = 0;
  for (Int_t i=0; i<nc; i++) {
    std::pair<Double_t, Double_t> fadc = GetWaveform(HodoRawHit::kUp, i);
    if ( fadc.first >= graphStart && fadc.first <= graphEnd ) {
      if ( index<n_range ) {
        gr->SetPoint(index, fadc.first, fadc.second);
        if (GetName() == "RAYRAW") {
          gr->SetPointError(index, 0, y_err);
        } else {
          gr->SetPointError(index, 0, 0); 
        }
        index++;
      }
    }
  }

  m_TGraphC.push_back(gr);
  return true;
}

//_____________________________________________________________________________
Bool_t RayrawWaveformHit::MakeDifGraph(Int_t index_org)
{
  TGraphErrors *gr = m_TGraphC[index_org];
  Int_t n = gr->GetN() - 1;
  if (n <= 0) return false; 

  Double_t *refx = gr->GetX();
  Double_t *refy = gr->GetY();
  Double_t x[n], y[n];

  for (Int_t i=0; i<n; i++) {
    x[i] = refx[i];
    y[i] = refy[i+1] - refy[i];
  }
  TGraphErrors *gr_diff = new TGraphErrors(n, x, y);
  m_TGraphC.push_back(gr_diff);
  return true;
}

//_____________________________________________________________________________
// (BGO PreSearch
/*Bool_t RayrawWaveformHit::PreSearch(struct SearchParam *sp)
{
  ...
}*/


Bool_t RayrawWaveformHit::PreSearch(struct SearchParam *sp)
{
  static const std::string func_name(std::string("[")+ClassName()+"::"+__func__+"()]");
  
  if (m_TGraphC.empty() || !m_TGraphC[0]) return false;
  TGraphErrors *gr = m_TGraphC[0];
  if (gr->GetN() < 1) return false;
  
  //param
  const double ADC_THRESHOLD = 10.0; 

  double max_val = -1e9;
  double peak_time = TMath::QuietNaN();
  
  for (int i = 0; i < gr->GetN(); ++i) {
    if (gr->GetY()[i] > max_val) {
      max_val = gr->GetY()[i];
      peak_time = gr->GetX()[i];
    }
  }
  /*
  // debug info
  if (SegmentId() == 0) {
      std::cout << "!!! DEBUG [PreSearch]    : (Seg 0, Evt " << gUnpacker.get_event_number() << ") Max value in TGraph = " << max_val << std::endl;
      if (max_val <= ADC_THRESHOLD) {
          std::cout << "!!! INFO: PreSearch failed. Max value (" << max_val << ") <= Threshold (" << ADC_THRESHOLD << ")" << std::endl;
      } else {
          std::cout << "!!! SUCCESS: PreSearch PASSED. Max value (" << max_val << ") > Threshold (" << ADC_THRESHOLD << ")" << std::endl;
      }
  }
  */

  if (max_val > ADC_THRESHOLD) {
    sp->foundx.push_back(peak_time);
    sp->foundy.push_back(max_val);
  }
  
  if (sp->foundx.empty()) {
    return false;
  }

  return true;
}

//_____________________________________________________________________________
Bool_t RayrawWaveformHit::WidthCut(std::vector<Double_t> rise,
         std::vector<Double_t> fall,
         Double_t width, std::vector<Double_t> &outrise)
{
  // do not use 
  if(rise.size() != fall.size()){
    std::cout<<"RayrawWaveformHit::WidthCut rise num != fall num"<<std::endl;
    return false;
  }
  for(unsigned int i=0;i<rise.size();i++){
    if( fall[i]-rise[i] > width){
      outrise.push_back(rise[i]);
    }
  }
  rise.clear();
  fall.clear();
  return true;
}

//_____________________________________________________________________________
void RayrawWaveformHit::CompareRise(std::vector<Double_t> rise1,
          std::vector<Double_t> rise2,
          Double_t width, std::vector<Double_t> &outrise)
{
  // do not use
  for(int i=0; i<(int)rise2.size(); i++){
    int t=0;
    for(int j=0; j<(int)rise1.size(); j++){
      if( rise2[i]-rise1[j] <width )
        t++;
    }
    if(t==0)
      outrise.push_back(rise2[i]);
  }
  for(unsigned int j=0;j<rise1.size();j++)
    outrise.push_back(rise1[j]);
  rise1.clear();
  rise2.clear();
}

//_____________________________________________________________________________
void RayrawWaveformHit::SetInitial(std::vector<Double_t> &v,
         Double_t begin, Double_t end,
         Double_t thre, Double_t rise)
{
  // do not use
  Int_t size=v.size();
  Double_t SepaLimit = 1.0; 
  
  if(size>1)
    for(Int_t i=0;i<size;i++){
      if(v[i]<begin || v[i]>end)
        v[i]=-1;
    }

  if(size>1){
    std::sort(v.begin(),v.end());
    for(Int_t i=0;i<size-1;i++){
      if(v[i+1]-v[i]<SepaLimit){
        v[i+1] = (Double_t)((v[i+1]+v[i])/2);
        v[i]=-1;
      }
    }
  }

  Int_t index_original_graph = 0;
  for(Int_t i=0;i<size;i++)
    if(v[i]!=-1){
      Bool_t flagOverThr = false;
      for (Double_t ratio =0; ratio <= 1.0; ratio += 0.1)
        if(GXtoGY(index_original_graph, v[i]+rise*ratio) > thre) 
          flagOverThr = true;
      if (!flagOverThr)
        v[i]=-1;
    }

  std::sort(v.begin(),v.end(),std::greater<double>());
  std::vector<double>::iterator it = std::find(v.begin(),v.end(),-1);
  if(it!=v.end())
    v.erase(it,v.end());
  std::sort(v.begin(),v.end());
  size=v.size();
}

//_____________________________________________________________________________

Double_t RayrawWaveformHit::GXtoGY(Int_t index_graph, Double_t gx)
{
  Int_t point=-1;
  Double_t k,l;
  if (index_graph >= (int)m_TGraphC.size() || !m_TGraphC[index_graph]) return 0;
  Int_t sample_size = m_TGraphC[index_graph]->GetN();
  if (sample_size < 2) return 0;
  
  Double_t *GX = m_TGraphC[index_graph]->GetX();
  Double_t *GY = m_TGraphC[index_graph]->GetY();

  if (gx < GX[0]) return GY[0];
  if (gx > GX[sample_size-1]) return GY[sample_size-1];

  for(Int_t i=1;i<sample_size;i++){ 
    if(gx<=GX[i]){
      point = i ;
      break;
    }
  }

  if(point <= 0)
    return GY[0];
  else{
    k = GX[point-1];
    l = GX[point];
    
    if (std::abs(l-k) < 1e-9) return GY[point-1];

    Double_t r =(gx-k)/(l-k);
    k = GY[point-1];
    l = GY[point];
    return k + (l-k)*r;
  }
  return 0;
}

//_____________________________________________________________________________

Bool_t RayrawWaveformHit::SetFitParam(FitParam *fp, std::vector<Double_t> &inix,
            std::vector<Double_t> &iniy)
{
  Int_t wm = inix.size();
  fp->wavenum = wm;
  fp->ParaNum = wm*2+2;

  if(fp -> ParaNum > ParaMax){
    std::cout<<"RayrawWaveformHit::SetFitParam: Too many Fitting Prameter"<<std::endl;
    return false;
  }

  fp->par[0]=wm;
  fp->par[1]=0; 

  for(Int_t i=0;i<wm;i++){
    fp->par[2+2*i]=inix[i];
    fp->par[3+2*i]=iniy[i];
  }
  return true;
}

//_____________________________________________________________________________

void RayrawWaveformHit::Fit1(FitParam *fp)
{
  Int_t wavenum = fp->wavenum;
  Int_t ParaNum = fp->ParaNum;
  Double_t par[ParaNum];
  par[0]=wavenum;

  for(Int_t i=0;i<ParaNum;i++){
    m_func -> ReleaseParameter(i);
    par[i] = fp->par[i];
    
    if (i>=3 && i%2 == 1) 
      m_func -> SetParLimits(i, 0, 100000); 
  }

  m_func -> SetNpx(1000);
  m_func -> SetParameters(&par[0]);
  m_func -> FixParameter(0,par[0]); 
  m_func -> SetLineColor(fp->color);
  for(Int_t nine= ParaNum;nine<ParaMax;nine++)
    m_func -> FixParameter(nine,0);

  m_TGraphC[fp->tgen] -> Fit(m_func,"qN","",fp->FitStart,fp->FitEnd);

  for(Int_t i =0;i<ParaNum;i++){
    fp -> FitParam[i] = m_func -> GetParameter(i);
  }
}

//_____________________________________________________________________________
Double_t RayrawWaveformHit::FittedTrigX(FitParam fp, Double_t allowance)
{
 
  Int_t num = fp.wavenum; 
  Double_t reso = TrigTimeReso *allowance;
  Int_t inrange=0;
  std::vector<Double_t> xx;
  for(Int_t i=0;i<num;i++){
    Double_t x =fp.FitParam[2+2*i]; 
    if(x > - reso && x< reso ){
      inrange++;
      xx.push_back(x);
    }
  }
  if(inrange==0)
    return -9999.;
  else if(inrange==1)
    return xx[0];
  else{
    std::vector<Double_t> sub;
    for(Int_t i=0;i<inrange;i++)
      sub.push_back(std::abs(xx[i]));
    std::sort(sub.begin(),sub.end());
    for(Int_t i=0;i<inrange;i++){
      if(std::abs(-sub[0] - xx[i]) <0.0001) 
        return xx[i];
      if(std::abs(sub[0] - xx[i]) <0.0001) 
        return xx[i];
    }
    return -9999.;
  }
}

//_____________________________________________________________________________

Double_t RayrawWaveformHit::RisingResidual(Int_t tge_No, Double_t trig, Double_t &res_max)
{
  if(trig <0)
    return -1;

  Double_t Residual = 0, max_res = 0;
  Double_t a,b;

  
  for(Int_t i =0; i<10;i++){
    Double_t x = trig - 0.05 + i*0.01; 
    a = GXtoGY(tge_No, x);
    b = m_func->Eval(x);
    if(std::abs(a) > 1e-9) { 
      Residual += sqrt((a-b)*(a-b))/ y_err ;
      if (std::abs(max_res) < std::abs(a-b))
        max_res = a-b;
    }
  }

  
  if (std::abs(GXtoGY(tge_No, trig)) > 1e-9) { 
    Residual /= GXtoGY(tge_No, trig); 
    Residual *= 100;
  } else {
    Residual = 9999.;
  }

  res_max = max_res;
  return Residual;
}

//_____________________________________________________________________________

void
RayrawWaveformHit::Print(Option_t* arg) const
{
 
  PrintHelper helper(3, std::ios::fixed);
  hddaq::cout << FUNC_NAME << " " << arg << std::endl
        << "detector_name = " << m_raw->DetectorName() << std::endl
        << "detector_id   = " << m_raw->DetectorId() << std::endl
        << "plane_name    = " << m_raw->PlaneName()  << std::endl // ★★★ 修正 (:: -> ->) ★★★
        << "plane_id      = " << m_raw->PlaneId()    << std::endl // ★★★ 修正 (:: -> ->) ★★★
        << "segment_id    = " << m_raw->SegmentId()  << std::endl // ★★★ 修正 (:: -> ->) ★★★
              << "n_ch          = " << m_n_ch              << std::endl
              << "de            = " << DeltaE() << std::endl
              << "mt/cmt        = " << MeanTime()
              << " / " << CMeanTime() << std::endl
              << "tdiff/ctdiff  = " << TimeDiff()
              << " / " << CTimeDiff() << std::endl;
  for(const auto data_map: std::map<TString, data_t>
        {{"de-hi  ", m_de_high},      {"de-lo  ", m_de_low},
         {"time-l ", m_time_leading}, {"time-t ", m_time_trailing},
         {"ctime-l", m_ctime_leading}, {"ctime-t", m_ctime_trailing}
        }){
    for(const auto& cont: data_map.second){
      hddaq::cout << " " << data_map.first << ":" << cont.size()
                  << " ";
      std::copy(cont.begin(), cont.end(),
                std::ostream_iterator<Double_t>(hddaq::cout," "));
    }
    hddaq::cout << std::endl;
  }
  hddaq::cout << std::endl;
}
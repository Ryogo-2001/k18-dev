/**
 *  file: BGOAnalyzer.cc
 *  date: 2018.10.23
 *
 */


#include "BGOAnalyzer.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "RawData.hh"
#include "BGOFitFunction.hh"
#include "BGODiscriminator.hh"
#include "BGOCalibMan.hh"
#include "HodoParamMan.hh"


namespace
{
  const std::string& class_name("BGOAnalyzer");
  const double Time_SI = 1./33.33333;
  const double y_err = 30.;
  const double y_err_NoData = 300.;
  //const double fitStart = 2.7;
  //const double fitEnd   = 4.5;
  const double graphStart = -3.0;
  const double graphEnd   = 2.0;
  const double fitStart = -1.0;
  const double fitEnd   = 1.0;
  const double SepaLimit = 0.08;
  const int ParaMax = 64;
  //const double TrigTime = 3.70;
  const double TrigTimeReso = 1.00;

  const HodoParamMan& gHodo  = HodoParamMan::GetInstance();
  const BGOCalibMan&  gCalib = BGOCalibMan::GetInstance();

}

//______________________________________________________________________________
BGOAnalyzer::BGOAnalyzer( void )
  : m_fitFunction(0), m_func(0)
{
  debug::ObjectCounter::increase(class_name);

  const BGOTemplateManager& gBGOTemp = BGOTemplateManager::GetInstance(); 

  m_fitFunction = gBGOTemp.GetFitFunction();
  m_func = new TF1("m_func", m_fitFunction, fitStart, fitEnd, ParaMax );

  for (int seg=0; seg<NumOfSegBGO; seg++) {
    double tdc0    = gHodo.GetOffset(DetIdBGO, 0, seg, 1); // FADC UorD = 1
    double pedestal = gHodo.GetP0(DetIdBGO, 0, seg, 1);     // FADC UorD = 1
    
    m_Tdc0[seg] = tdc0;
    m_Pedestal[seg] = pedestal;
  }

}

//______________________________________________________________________________
BGOAnalyzer::~BGOAnalyzer( void )
{
  ClearBGOFadcHits();
  if (m_func) {
    delete m_func;
    m_func = 0;
  }

  debug::ObjectCounter::decrease(class_name);
}


//______________________________________________________________________________
void
BGOAnalyzer::ClearBGOFadcHits( void )
{

  for (int i=0; i<NumOfSegBGO; i++) {
    m_fadcContBGO[i].clear();
    m_fadcContBGO[i].shrink_to_fit();
    m_FitParamCont[i].clear();
    m_FitParamCont[i].shrink_to_fit();
    m_BGODataCont[i].clear();
    m_BGODataCont[i].shrink_to_fit();
    del::ClearContainer( m_TGraphCont[i] );
    del::ClearContainer( m_TF1Cont[i] );
  }
}


//______________________________________________________________________________
bool BGOAnalyzer::DecodeBGO( RawData *rawData )
{
  ClearBGOFadcHits();

  /*
  double Offset[NumOfSegBGO] = {
    0, 0, 16040, 16120, 16060, 
    16090, 15970, 15980, 15920, 15540,
    15530, 16290, 16370, 15350, 16350,
    16060, 16060, 16070, 16080, 15540,
    15480, 16300, 15620, 16120
  };
  */
  for (int seg=2; seg<NumOfSegBGO; seg++) {

    const FADCRHitContainer &cont = rawData->GetBGOFAdcRawHC(seg);

    int nhit = cont.size();
    for (int i=0; i<nhit; i++) {
      FAdcData fadc;
      fadc.x = ((double)i - m_Tdc0[seg])* Time_SI;

      if (cont[i] < 0xffff) {
	fadc.y = (double)(cont[i]-m_Pedestal[seg]);
	fadc.err = y_err;
      }
      else {
	fadc.y = 0.;
	fadc.err = y_err_NoData;
      }

      m_fadcContBGO[seg].push_back(fadc);

    }
  }

  return true;
}


bool BGOAnalyzer::PulseSearch( void )
{

  for (int seg = 2; seg<NumOfSegBGO; seg++) {
    if (m_fadcContBGO[seg].size() == 0)
      continue;

    MakeGraph(seg);
    // Original graph index = 0;
    int index_original_graph = 0;

    MakeDifGraph(seg, index_original_graph);
    int index_diff_graph = 1;

    double threshold = -100;
    double width = 0.05;
    double risetime = 0.1;

    SearchParam sp1={"sp1", {index_original_graph, index_diff_graph},
		     fitStart, fitEnd, fitStart, fitEnd, 
		     threshold, width, risetime};

    bool flagPresearch = PreSearch(seg, &sp1);
    if (!flagPresearch)
      continue;

    int color = 4;
    FitParam fp1={"fp1", index_original_graph, color ,fitStart, fitEnd};

    SetFitParam(&fp1,sp1.foundx,sp1.foundy);

    Fit1(seg, &fp1);

    //std::cout << "BGO Fitting result ; segment " << seg << std::endl;
    //for (int i=0; i<fp1.ParaNum; i++) {
    //std::cout << "p" << i << " : " << fp1.FitParam[i] << std::endl;
    //}

    double trigx = FittedTrigX(fp1,1.0);
    //trigx =1;/////////////////////                                             
    if(fabs(trigx)<TrigTimeReso){
      double max_res = 0;
      fp1.Residual=  RisingResidual(seg, index_original_graph, trigx, max_res);
      //fp1.Residual=1;/////////////

      m_FitParamCont[seg].push_back(fp1);
      TF1 *func = new TF1(Form("func%d_0", seg), m_fitFunction, 
			  fitStart, fitEnd, ParaMax );
      func->SetLineColor(fp1.color);
      for (int i=0; i<fp1.ParaNum; i++) {
	func->SetParameter(i, fp1.FitParam[i]);
      }
      m_TF1Cont[seg].push_back(func);

      if(fp1.Residual <5 || fabs(max_res)<50){
	int waveNum = fp1.wavenum;
	for (int nw=0; nw<waveNum; nw++) {
	  BGOData bgoData;
	  bgoData.chi2 = fp1.Residual;
	  bgoData.time = fp1.FitParam[2*nw+2];
	  bgoData.pulse_height = fp1.FitParam[2*nw+3];
	  double energy=-999.;
	  if (gCalib.GetEnergy(seg, bgoData.pulse_height, energy))
	    bgoData.energy = energy;
	  else
	    bgoData.energy = -999.;

	  m_BGODataCont[seg].push_back(bgoData);
	}


	// finish analysis of this segment
	// goto next segment
	continue;
      }
    }

  }

  return true;
}

bool BGOAnalyzer::MakeGraph(int seg)
{
  
  int nc = m_fadcContBGO[seg].size();

  int n_range = 0;
  for (int i=0; i<nc; i++) {
    FAdcData fadc = m_fadcContBGO[seg][i];
    //if ( fadc.x >= fitStart && fadc.x <= fitEnd ) 
    if ( fadc.x >= graphStart && fadc.x <= graphEnd ) 
      n_range++;
  }

  TGraphErrors *gr = new TGraphErrors(n_range);
  int index = 0;
  for (int i=0; i<nc; i++) {
    FAdcData fadc = m_fadcContBGO[seg][i];
    //if ( fadc.x >= fitStart && fadc.x <= fitEnd ) {
    if ( fadc.x >= graphStart && fadc.x <= graphEnd ) {
      if ( index<n_range ) {
	gr->SetPoint(index, fadc.x, fadc.y);
	gr->SetPointError(index, 0, fadc.err);
	index++;
      }

    }
  }
  
  m_TGraphCont[seg].push_back(gr);

  return true;

}

bool BGOAnalyzer::MakeDifGraph(int seg, int index_org)
{
  TGraphErrors *gr = m_TGraphCont[seg][index_org];

  int n = gr->GetN() - 1;

  double *refx = gr->GetX();
  double *refy = gr->GetY();

  double x[n], y[n];

  for (int i=0; i<n; i++) {
    x[i] = refx[i];
    y[i] = refy[i+1] - refy[i];
  }

  TGraphErrors *gr_diff = new TGraphErrors(n, x, y);

  m_TGraphCont[seg].push_back(gr_diff);

  return true;

}

bool BGOAnalyzer::PreSearch(int seg, struct SearchParam *sp)
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  double begin = sp ->sbegin;
  double end   = sp ->send;

  //BGODiscriminator BGODiscri();
  BGODiscriminator BGODiscri;
  BGODiscri.SetGraph(m_TGraphCont[seg][sp->tgen[0]]);

  std::vector<double> OriR,OriF;
  double threshold = -500.;
  bool flagSelectRange = BGODiscri.SelectRange(threshold, begin, end);

  if (!flagSelectRange)
    return false;

  BGODiscri.GetRisingPoint(OriR);
  BGODiscri.GetFallingPoint(OriF);
  BGODiscri.AllClear();

  //bool FF = FlagFrontWave(OriR,begin);
  bool flagFrontWave = false;
  if (OriR.size()>0)
    if (OriR[0] - begin < 0.1)
      flagFrontWave = true;

  BGODiscri.SetGraph(m_TGraphCont[seg][sp->tgen[1]]);

  std::vector<double> rise30,fall30,entry30;
  threshold = -30;

  flagSelectRange = BGODiscri.SelectRange(threshold, begin, end);
  if (!flagSelectRange)
    return false;

  BGODiscri.GetRisingPoint(rise30);
  BGODiscri.GetFallingPoint(fall30);

  double width_thr = 0.05;
  bool flagWidthCut1 = WidthCut(rise30,fall30,width_thr,entry30);
  if (!flagWidthCut1) {
    std::cout << func_name << " : WidthCut1, Number of rise and fall points does not match" << std::endl;
    return false;
  }

  ///Bool_t FF2= FlagFrontWave(rise30,begin);
  bool flagFrontWave2 = false;
  if (rise30.size()>0)
    if (rise30[0] - begin < 0.1)
      flagFrontWave2 = true;


  BGODiscri.Clear();

  int flagAroundBegin = 0;

  if(flagFrontWave && flagFrontWave2)
    flagAroundBegin = 1;
  if(flagFrontWave && !flagFrontWave2)
    flagAroundBegin = 2;
  if(!flagFrontWave && flagFrontWave2){
    //    std::cout<<"Event:"<<CEI->GetEventNum()<<" ch:"<<CEI->GetChannel()<<"Search1 ahoooooo"<<std::endl;
  }
  //  std::cout<<"Event FF:"<<FF<<std::endl;

  std::vector<double> rise100,fall100,entry100,entry30_100;
  threshold = -100;

  flagSelectRange = BGODiscri.SelectRange(threshold , begin, end);
  if (!flagSelectRange)
    return false;

  BGODiscri.GetRisingPoint(rise100);
  BGODiscri.GetFallingPoint(fall100);
  width_thr = 0.08;
  WidthCut(rise100, fall100, width_thr, entry100);

  CompareRise(entry30,entry100,0.05,entry30_100);
  //  for(int i=0;i<entry30_100.size();i++)                                                                
  //    std::cout<<"entry "<<entry30_100[i]<<std::endl;                                                    

  SetInitial(seg, entry30_100,begin,end,sp->threshold,sp->risetime);

  int index_original_graph = 0;
  for(unsigned int i=0;i<entry30_100.size();i++){
    sp->foundx.push_back(entry30_100[i]+sp->risetime);
    sp->foundy.push_back(-GXtoGY(seg, index_original_graph, 
				 entry30_100[i]+sp->risetime)) ;
  }

  if(flagAroundBegin==1){
    sp->foundx.insert(sp->foundx.begin(),begin);
    sp->foundy.insert(sp->foundy.begin(),
		      -GXtoGY(seg, index_original_graph, begin+0.03));
  }
  if(flagAroundBegin==2){
    sp->foundx.insert(sp->foundx.begin(),begin-0.5);
    sp->foundy.insert(sp->foundy.begin(),
		      -GXtoGY(seg, index_original_graph, begin+0.03));
  }

  return true;

}

bool BGOAnalyzer::WidthCut(std::vector<double> rise,
			   std::vector<double> fall,
			   double width, std::vector<double> &outrise)
{

  if(rise.size() != fall.size()){
    std::cout<<"CAlgorithm::WidthCut rise num != fall num"<<std::endl;
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

void BGOAnalyzer::CompareRise(std::vector<double> rise1,
			      std::vector<double> rise2,
			      double width, std::vector<double> &outrise)
{
  for(int i=0;i<rise2.size();i++){
    int t=0;
    for(int j=0;j<rise1.size();j++){
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

void BGOAnalyzer::SetInitial(int seg, std::vector<double> &v, 
			     double begin, double end,
			     double thre, double rise)
{
  int size=v.size();

  if(size>1)
    for(int i=0;i<size;i++){
      if(v[i]<begin || v[i]>end)
        v[i]=-1;
    }

  if(size>1){
    std::sort(v.begin(),v.end());
    for(int i=0;i<size-1;i++){
      if(v[i+1]-v[i]<SepaLimit){
        v[i+1] = (double)((v[i+1]+v[i])/2);
        v[i]=-1;
      }
    }
  }

  int index_original_graph = 0;
  for(int i=0;i<size;i++)
    if(v[i]!=-1){
      //      std::cout<<"GX "<<v[i]+0.15<<" GY "<<Original->GXtoGY(v[i]+0.15)<<std::endl;                 
      bool flagOverThr = false;
      for (double ratio =0; ratio <= 1.0; ratio += 0.1)
	if(GXtoGY(seg, index_original_graph, v[i]+rise*ratio)<thre)
	  flagOverThr = true;

      if (!flagOverThr)
	v[i]=-1;
      /*
      if(GXtoGY(0,v[i]+rise)>thre)
	v[i]=-1;
      */
    }

  std::sort(v.begin(),v.end(),std::greater<double>());

  std::vector<double>::iterator it = std::find(v.begin(),v.end(),-1);
  if(it!=v.end())
    v.erase(it,v.end());

  std::sort(v.begin(),v.end());
  size=v.size();

}


double BGOAnalyzer::GXtoGY(int seg, int index_graph, double gx)
{
  int point=-1;
  double k,l;
  int sample_size = m_TGraphCont[seg][index_graph]->GetN();
  double *GX = m_TGraphCont[seg][index_graph]->GetX();
  double *GY = m_TGraphCont[seg][index_graph]->GetY();

  for(int i=0;i<sample_size;i++){
    if(gx<=GX[i]){
      point = i ;
      break;
    }

  }
  if(point == -1)
    return 0;

  else{
    k = GX[point-1];
    if(point == sample_size-1)
      return GY[point];
    else
      l = GX[point];

    Double_t r =(gx-k)/(l-k);


    k = GY[point-1];
    l = GY[point];
    return k + (l-k)*r;
  }

  return 0;
}

void BGOAnalyzer::SetFitParam(FitParam *fp, std::vector<double> &inix,
			      std::vector<double> &iniy)
{
  int wm = inix.size();
  fp->wavenum = wm;
  fp->ParaNum = wm*2+2;

  if(fp -> ParaNum > ParaMax){
    std::cout<<"Too many Fitting Prameter"<<std::endl;
    std::exit(-1);
  }


  fp->par[0]=wm;
  fp->par[1]=0;

  for(int i=0;i<wm;i++){
    fp->par[2+2*i]=inix[i];
    fp->par[3+2*i]=iniy[i];
    //    std::cout<<"InitialY"<<InitialY[i]<<"fp y"<<fp->par[3+2*i]<<" iniy "<<iniy[i]<<std::endl;        
  }


}


void BGOAnalyzer::Fit1(int seg, FitParam *fp)
{

  int wavenum = fp->wavenum;
  int ParaNum = fp->ParaNum;
  double par[ParaNum];
  par[0]=wavenum;
  //std::cout << "wavenum = "  << wavenum << std::endl;   

  for(Int_t i=0;i<ParaNum;i++){
    m_func -> ReleaseParameter(i);
    par[i] = fp->par[i];
    //std::cout << "par[" << i << "] = "  << par[i] << std::endl;
  }

  //m_func -> SetName(Name.str().c_str());
  m_func -> SetNpx(1000);
  m_func -> SetParameters(&par[0]);
  m_func -> FixParameter(0,par[0]);
  //m_func -> FixParameter(1,0);                                                                             
  m_func -> SetLineColor(fp->color);
  for(int nine= ParaNum;nine<ParaMax;nine++)
    m_func -> FixParameter(nine,0);


  m_TGraphCont[seg][fp->tgen] -> Fit(m_func,"qN","",fp->FitStart,fp->FitEnd);
  //  std::cout<< "fit: " << fp->FitStart << "~" << fp->FitEnd << std::endl;

  for(Int_t i =0;i<ParaNum;i++){
    fp -> FitParam[i] = m_func -> GetParameter(i);
  }

  //std::cout<< "fit: " << fp->FitStart << "~" << fp->FitEnd  << ", ch" << GetChannel() << ", tgen=" << fp->tgen << std::endl;

}

double BGOAnalyzer::FittedTrigX(FitParam fp,double allowance)
{
  int num = fp.wavenum;
  double reso = TrigTimeReso *allowance;

  int inrange=0;

  std::vector<double> xx;
  for(int i=0;i<num;i++){
    double x =fp.FitParam[2+2*i];
    //std::cout<<"x "<<x<<"  ref "<<mean<<std::endl;                                                       
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
    std::vector<double> sub;
    for(int i=0;i<inrange;i++)
      sub.push_back(fabs(xx[i]));
    std::sort(sub.begin(),sub.end());
    for(int i=0;i<inrange;i++){
      if(-sub[0] - xx[i] <0.0001)
        return xx[i];
      if(sub[0] - xx[i] <0.0001)
        return xx[i];
    }
    return -9999.;
  }

}

double BGOAnalyzer::RisingResidual(int seg, int tge_No, double trig, double &res_max)
{

  if(trig <0)
    return -1;

  double Residual = 0, max_res = 0;
  double a,b;

  for(int i =0; i<10;i++){
    double x = trig - 0.1 + i*0.01;
    a = GXtoGY(seg, tge_No, x);
    b = m_func -> Eval(x);
    if(a !=0) {
      //Residual += sqrt((a-b)*(a-b))/ CEI->GetError() ;                                                   
      Residual += sqrt((a-b)*(a-b))/ 15 ; //kuso tekito                        

      if (fabs(max_res) < fabs(a-b))
        max_res = a-b;
                            
    }
    //std::cout<<"Residual x:"<<x<<"  a:"<<a<<"  b:"<<b<<std::endl; 
  }

  Residual /= -GXtoGY(seg, tge_No, trig);
  Residual *= 100;

  res_max = max_res;

  return Residual;
}


//______________________________________________________________________________
bool
BGOAnalyzer::GetBGOData0(int segment, BGOData &bgoData) 
{
  if( segment<0 || segment>NumOfSegBGO ) return false;

  int nc = m_BGODataCont[segment].size();
  if (nc == 0) return false;

  double time0 = -999.;
  int  index=0;
  for (int i=0; i<nc; i++) {
    if (std::abs(time0) > std::abs(m_BGODataCont[segment][i].time)) {
      time0 = m_BGODataCont[segment][i].time;
      index = i;
    }
  }
  
  bgoData = m_BGODataCont[segment][index];

  return  true;

}

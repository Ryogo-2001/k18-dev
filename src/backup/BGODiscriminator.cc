#include "BGODiscriminator.hh"

BGODiscriminator::~BGODiscriminator(){

  if(!m_cross.empty())
    m_cross.clear();
  if(!m_RisingPoint.empty())
    m_RisingPoint.clear();
  if(!m_FallingPoint.empty())
    m_FallingPoint.clear();
  if(!m_target.empty())
    m_target.clear();
  if(!m_FoundPeak.empty())
    m_FoundPeak.clear();
  //  std::cout<<"Deconstruct Discri"<<std::endl;
}


BGODiscriminator::BGODiscriminator( void )
  : m_sample_size(0)
{
}


void BGODiscriminator::SetGraph(TGraphErrors *graph)
{

  m_sample_size = graph->GetN();
  double *x = graph->GetX();
  double *y = graph->GetY();
  
  m_GraphX.clear();
  m_GraphY.clear();

  m_GraphX.resize(m_sample_size);
  m_GraphY.resize(m_sample_size);

  for(int i=0;i<m_sample_size;i++){
    m_GraphX[i] = x[i];
    m_GraphY[i] = y[i];
  }
}


void BGODiscriminator::Clear()
{
  if(!m_cross.empty())
    m_cross.clear();
  if(!m_RisingPoint.empty())
    m_RisingPoint.clear();
  if(!m_FallingPoint.empty())
    m_FallingPoint.clear();
  if(!m_target.empty())
    m_target.clear();
}


void BGODiscriminator::AllClear()
{
  if(!m_cross.empty())
    m_cross.clear();
  if(!m_RisingPoint.empty())
    m_RisingPoint.clear();
  if(!m_FallingPoint.empty())
    m_FallingPoint.clear();
  if(!m_target.empty())
    m_target.clear();
  m_sample_size=0;

  m_GraphX.clear();
  m_GraphY.clear();


}


bool BGODiscriminator::SelectRange(double threshold,double begin,double last)
{
  if(begin >= last){
    std::cout<<"Threshold search is failed. The search range is 'begin>=last' "<<std::endl;
    exit(-1);
  }
  
  for(int i=0;i<m_sample_size;i++){
    if(m_GraphX[i]>=begin && m_GraphX[i]<=last)
      m_target.push_back(m_GraphY[i]);

    else
      m_target.push_back(0);
  }

  if(m_target.size()==0){
    std::cout<<"Threshold search is failed. Data is nothing in the search range"<<std::endl;
    std::cout<<"begin=" << begin << "-> last=" << last <<std::endl;
    return false;
  }
  
  GetCrossPoint(threshold);
  
  return true;
  
};


bool BGODiscriminator::SelectRange(double threshold){
  for(int i=0;i<m_sample_size;i++)
    m_target.push_back(m_GraphX[i]);
  
  GetCrossPoint(threshold);

  return true;
};


void BGODiscriminator::GetRisingPoint(std::vector<double> &v){
  if(!v.empty())
    v.clear();
  /*
  if(RisingPoint.empty())
    std::cout<<"BGODiscriminator::GetRisingPoint : RisingPoint have no data "<<std::endl;
  */
  v.resize(m_RisingPoint.size());

  for(int i=0;i<m_RisingPoint.size();i++){
    //      std::cout<<"RisingPoint"<<RisingPoint[i]<<std::endl;
    //std::cout<<"v.size : "<<v.size()<<std::endl;
    v[i] =  m_GraphX[m_RisingPoint[i]] ;
  }
}
void BGODiscriminator::GetFallingPoint(std::vector<double> &v){
  if(!v.empty())
    v.clear();
  /*
  if(FallingPoint.empty())
    std::cout<<"BGODiscriminator::GeFallingPoint :FallingPoint have no data"<<std::endl;  
  */

  v.resize(m_FallingPoint.size());
  for(unsigned int i=0;i<m_FallingPoint.size();i++){
    v[i] = m_GraphX[m_FallingPoint[i]] ;
    
  }
  
}



bool BGODiscriminator::SloveBin(double input,double threshold,bool LorS){
  if(LorS){
    if(input>threshold)
      return true;
    else
      return false;
  }
  else
    if(input<threshold)
      return true;
    else
      return false;
}


int BGODiscriminator::GetPeakNum(double threshold,bool AbsFlag,bool LorS)
{
  int PeakNum=0;
  int DataSize = m_target.size();
  if(!m_cross.empty())
    m_cross.clear();
  
  bool thre[DataSize];
  if(AbsFlag)
    for(int i=0; i<DataSize;i++)
      thre[i]=SloveBin(fabs(m_target[i]),threshold,LorS);
  else
    for(int i=0; i<DataSize;i++)
      thre[i]=SloveBin(m_target[i],threshold,LorS);
      
  for(int i=0;i<DataSize;i++) {
    if(i>0) {
      if(thre[i-1] ^ thre[i]) { // XOR
	if(thre[i-1]<=thre[i])
	  m_cross.push_back(1);
	else
	  m_cross.push_back(2);
      } else {
	m_cross.push_back(0);
      }
    }
  }

  if (thre[0])
    m_cross.insert(m_cross.begin(),1);
  else
    m_cross.push_back(0);

  if (thre[DataSize-1])
    m_cross.push_back(2);
  else
    m_cross.push_back(0);
  for(int i=0;i<DataSize+1;i++)
    PeakNum += m_cross[i];
  
  if(m_target.size()!=m_cross.size()-1)
    std::cout<<"Error in BGODiscriminator::GetPeakNum"<<std::endl;
  
  
  PeakNum /= 3;
  return PeakNum;
}


int BGODiscriminator::GetCrossPoint(double threshold)
{
  int FoundNum = GetPeakNum(threshold,false,false);

  if(FoundNum>0){
    for(int i=0;i<m_cross.size();i++){
      switch(m_cross[i]){
      case 1:
	m_RisingPoint.push_back(i);
	break;
      case 2:
	m_FallingPoint.push_back(i);
	break;
      }
    }
    
  }
  return FoundNum;
}

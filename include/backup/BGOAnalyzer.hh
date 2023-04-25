/**
 *  file: BGOAnalyzer.hh
 *  date: 2017.04.10
 *
 */

#ifndef BGO_ANALYZER_HH
#define BGO_ANALYZER_HH

#include <vector>

#include "DetectorID.hh"
#include "RawData.hh"
#include "BGOFitFunction.hh"
#include "BGOCalibMan.hh"

#include "TGraphErrors.h"
#include "TF1.h"

class RawData;
class BGOFitFunction;
//______________________________________________________________________________

struct FAdcData 
{
  double x;
  double y;
  double err;
};

struct SearchParam{
  std::string Name;
  int tgen[2];
  double sbegin;
  double send;
  double cutbegin;
  double cutend;
  double threshold;
  double width;
  double risetime;
  std::vector<double> foundx;
  std::vector<double> foundy;
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

struct BGOData 
{
  double chi2;
  double time;
  double pulse_height;
  double energy;
};


typedef std::vector<FAdcData>      FAdcDataContainer;
typedef std::vector<TGraphErrors*> TGraphErrorsContainer;
typedef std::vector<TF1*>          TF1Container;
typedef std::vector<FitParam>      FitParamContainer;
typedef std::vector<BGOData>       BGODataContainer;


class BGOAnalyzer
{
public:
  BGOAnalyzer( void );
  ~BGOAnalyzer( void );

  static BGOAnalyzer& GetInstance( void );
  static const std::string& ClassName( void );

private:
  BGOAnalyzer( const BGOAnalyzer& );
  BGOAnalyzer& operator =( const BGOAnalyzer& );


  FAdcDataContainer     m_fadcContBGO[NumOfSegBGO];
  TGraphErrorsContainer m_TGraphCont[NumOfSegBGO];
  TF1Container          m_TF1Cont[NumOfSegBGO];
  FitParamContainer     m_FitParamCont[NumOfSegBGO];
  BGODataContainer      m_BGODataCont[NumOfSegBGO];

  double m_Pedestal[NumOfSegBGO];
  double m_Tdc0[NumOfSegBGO];

  BGOFitFunction *m_fitFunction;  
  TF1 *m_func;

public:
  inline int GetNHitBGO( void );
  int GetNGraph( int segment ) {return m_TGraphCont[segment].size(); };
  inline TGraphErrors * GetTGraph( int segment, std::size_t i );
  int GetNFunc( int segment ) {return m_TF1Cont[segment].size(); };
  inline TF1 * GetTF1( int segment, std::size_t i );
  int GetNPulse( int segment ) {return m_BGODataCont[segment].size(); };

  const FAdcDataContainer& GetBGOFadcCont(int segment) const;
  const BGODataContainer&  GetBGODataCont(int segment) const;

  bool GetBGOData0(int segment, BGOData &bgoData) ;
public:
  void ClearBGOFadcHits( void );
  bool DecodeBGO( RawData *rawData );
  bool PulseSearch( void );
  bool MakeGraph(int seg);
  bool MakeDifGraph(int seg, int index_org);
  bool PreSearch(int seg, struct SearchParam *sp);
  bool WidthCut(std::vector<double> rise,
		std::vector<double> fall,
		double width, std::vector<double> &outrise);
  void CompareRise(std::vector<double> rise1,
		   std::vector<double> rise2,
		   double width, std::vector<double> &outrise);
  void SetInitial(int seg, std::vector<double> &v, 
		  double begin, double end,
		  double thre, double rise);
  double GXtoGY(int seg, int index_graph, double gx);
  void SetFitParam(FitParam *fp, std::vector<double> &inix,
		   std::vector<double> &iniy);
  void Fit1(int seg, FitParam *fp);
  double FittedTrigX(FitParam fp,double allowance);
  double RisingResidual(int seg, int tge_No, double trig, double &res_max);



};

//______________________________________________________________________________
inline BGOAnalyzer&
BGOAnalyzer::GetInstance( void )
{
  static BGOAnalyzer g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const std::string&
BGOAnalyzer::ClassName( void )
{
  static std::string g_name("BGOAnalyzer");
  return g_name;
}

//______________________________________________________________________________
inline int 
BGOAnalyzer::GetNHitBGO( void )
{
  int nhit = 0;
  for (int seg=2; seg<NumOfSegBGO; seg++) {
    int ngr = GetNGraph(seg);
    if (ngr>0) 
      nhit++;
  }

  return nhit;
}

//______________________________________________________________________________
inline TGraphErrors*
BGOAnalyzer::GetTGraph( int segment, std::size_t i )
{
  if( i<m_TGraphCont[segment].size() )
    return m_TGraphCont[segment][i];
  else
    return 0;
}

//______________________________________________________________________________
inline TF1*
BGOAnalyzer::GetTF1( int segment, std::size_t i )
{
  if( i<m_TF1Cont[segment].size() )
    return m_TF1Cont[segment][i];
  else
    return 0;
}

//______________________________________________________________________________
inline const FAdcDataContainer&
BGOAnalyzer::GetBGOFadcCont(int segment) const
{
  if( segment<0 || segment>NumOfSegBGO ) segment=0;
  return m_fadcContBGO[segment];

}

//______________________________________________________________________________
inline const BGODataContainer&
BGOAnalyzer::GetBGODataCont(int segment) const
{
  if( segment<0 || segment>NumOfSegBGO ) segment=0;
  return m_BGODataCont[segment];

}


#endif

#ifndef CFTPOS_PARAM_MAN_H
#define CFTPOS_PARAM_MAN_H

#include <string>
#include "DetectorID.hh"

const int NumOfPhiLayer = 4;
const int NumOfAngle = 36;
const int NumOfPhiPar = 6;

const int NumOfULayer = 4;
const int NumOfUPar = 6;
const int NumOfUDataPar = 37;
const int DivPhi=10; //degree

const int NumOfUDataPar2 = 19;
const int NumOfTheta = 10;
const int DivPhi2=20; //degree

class CFTPosParamMan
{
private:
  bool       m_is_ready;
  TString    m_file_name;

  double PhiPosPar[NumOfPhiLayer][NumOfAngle][NumOfPhiPar];
  double UPosPar[NumOfULayer][NumOfSegCFT_UV4][NumOfUPar];
  double UPosDataPar[NumOfULayer][NumOfSegCFT_UV4][NumOfUDataPar];
  double UPosDataPar2[NumOfULayer][NumOfSegCFT_UV4][NumOfTheta][NumOfUDataPar2];
public:
  CFTPosParamMan();
  ~CFTPosParamMan();

  static CFTPosParamMan& GetInstance( void );
  static const std::string& ClassName( void );

  bool Initialize( void );
  bool Initialize( const TString& file_name );
  bool IsReady( void ) const { return m_is_ready; }

  TString GetFileName( void ) const { return m_file_name; };

  double GetPhiShift(int layer, int angle, double z) const;
  double PhiPosFunc(int layer, int angle, double z) const;
  double  GetZposU(int layer, int seg, double phi, double theta) const;
  int GetThetaId(double theta) const;

};

//______________________________________________________________________________
inline CFTPosParamMan&
CFTPosParamMan::GetInstance( void )
{
  static CFTPosParamMan g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const std::string&
CFTPosParamMan::ClassName( void )
{
  static std::string g_name("CFTPosParamMan");
  return g_name;
}



#endif

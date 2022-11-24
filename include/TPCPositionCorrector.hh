// -*- C++ -*-

#ifndef TPC_POSITION_CORRECTOR_HH
#define TPC_POSITION_CORRECTOR_HH

#include <string>
#include <vector>

#include <TString.h>
#include <TVector3.h>

//_____________________________________________________________________________
class TPCPositionCorrector
{
public:
  static const TString& ClassName( void );
  static TPCPositionCorrector& GetInstance( void );
  ~TPCPositionCorrector( void );

private:
  TPCPositionCorrector( void );
  TPCPositionCorrector( const TPCPositionCorrector& );
  TPCPositionCorrector & operator =( const TPCPositionCorrector& );

private:
  typedef std::vector<std::vector<std::vector<TVector3>>> CorrectionMap;
  Bool_t          m_is_ready;
  TString         m_file_name;
  CorrectionMap   m_correction_map;
  Int_t           m_n_x;
  Int_t           m_n_y;
  Int_t           m_n_z;
  Double_t        m_x0;
  Double_t        m_y0;
  Double_t        m_z0;
  Double_t        m_dx;
  Double_t        m_dy;
  Double_t        m_dz;

public:
  TVector3 Correct( const TVector3& pos ) const;
  TVector3 GetCorrectionVector( const TVector3& pos ) const;
  Bool_t   Initialize( void );
  Bool_t   Initialize( const TString& file_name );
  Bool_t   IsReady( void ) const { return m_is_ready; }
  void     SetFileName( const TString& file_name )
    { m_file_name = file_name; }
};

//_____________________________________________________________________________
inline const TString&
TPCPositionCorrector::ClassName( void )
{
  static TString s_name( "TPCPositionCorrector" );
  return s_name;
}

//_____________________________________________________________________________
inline TPCPositionCorrector&
TPCPositionCorrector::GetInstance( void )
{
  static TPCPositionCorrector s_instance;
  return s_instance;
}

#endif

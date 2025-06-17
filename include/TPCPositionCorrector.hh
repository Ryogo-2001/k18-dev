// -*- C++ -*-

#ifndef TPC_POSITION_CORRECTOR_HH
#define TPC_POSITION_CORRECTOR_HH

#include <string>
#include <vector>
#include <TVector3.h>
#include <TString.h>

typedef std::vector<std::vector<TVector3>> CorrectionMap;

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

  std::vector<CorrectionMap> m_map;
  Bool_t          m_is_ready;
  TString         m_file_name;
  Int_t           m_n_y;
  Double_t        m_y0;
  Double_t        m_dy;

public:
  TVector3 Correct( const TVector3& pos, Int_t layer, Int_t row ) const;
  TVector3 GetCorrectionVector( const TVector3& pos, Int_t layer, Int_t row ) const;
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

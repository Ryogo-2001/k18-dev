// -*- C++ -*-

#ifndef MST_PARAM_MAN_HH
#define MST_PARAM_MAN_HH

#include <vector>
#include <map>
#include <TString.h>

//_____________________________________________________________________________
class MsTParamMan
{
public:
  static const TString& ClassName( void );
  static MsTParamMan&   GetInstance( void );
  ~MsTParamMan( void );

private:
  MsTParamMan( void );
  MsTParamMan( const MsTParamMan& );
  MsTParamMan& operator =( const MsTParamMan& );

private:
  Bool_t                             m_is_ready;
  Int_t                              m_nA;
  Int_t                              m_nB;
  std::vector<std::vector<Double_t>> m_low_threshold;
  std::vector<std::vector<Double_t>> m_high_threshold;

public:
  Double_t GetLowThreshold( Int_t detA, Int_t detB ) const;
  Double_t GetHighThreshold( Int_t detA, Int_t detB ) const;
  Bool_t   Initialize( const TString& filename );
  Bool_t   IsAccept( Int_t detA, Int_t detB, int tdc ) const;
  Bool_t   IsReady( void ) const { return m_is_ready; }
  void     Print( const TString& arg="" ) const;
  void     Print( Int_t detA, Int_t detB, int tdc ) const;
};

//_____________________________________________________________________________
inline const TString&
MsTParamMan::ClassName( void )
{
  static TString s_name( "MsTParamMan" );
  return s_name;
}

//_____________________________________________________________________________
inline MsTParamMan&
MsTParamMan::GetInstance( void )
{
  static MsTParamMan s_instance;
  return s_instance;
}

#endif

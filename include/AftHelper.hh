


// -*- C++ -*-

#ifndef AFT_HELPER_HH
#define AFT_HELPER_HH

#include <map>
#include <TString.h>
#include <TVector3.h>

#include "DetectorID.hh"

// //_____________________________________________________________________________
// class TpcPadParam
// {
// public:
//   TpcPadParam( Int_t asad, Int_t aget, Int_t channel,
//                Int_t pad_id, Int_t layer_id, Int_t row_id );
//   ~TpcPadParam( void );

// public:
//   static const Int_t AsAdMask    = 0xff;
//   static const Int_t AGetMask    = 0xff;
//   static const Int_t ChannelMask = 0xff;
//   static const Int_t AsAdShift    =  0;
//   static const Int_t AGetShift    =  8;
//   static const Int_t ChannelShift = 16;
//   static Int_t Key( Int_t asad, Int_t aget, Int_t channel )
//     {
//       return ( ( ( asad    & AsAdMask )    << AsAdShift    ) |
//                ( ( aget    & AGetMask )    << AGetShift    ) |
//                ( ( channel & ChannelMask ) << ChannelShift ) );
//     }

// private:
//   Int_t m_asad;
//   Int_t m_aget;
//   Int_t m_channel;
//   Int_t m_pad_id;
//   Int_t m_layer_id;
//   Int_t m_row_id;
//   Int_t m_key;

// public:
//   Int_t AsAdId( void ) const { return m_asad; }
//   Int_t AGetId( void ) const { return m_aget; }
//   Int_t Channel( void ) const { return m_channel; }
//   Int_t PadId( void ) const { return m_pad_id; }
//   void  Print( void ) const;
//   Int_t LayerId( void ) const { return m_layer_id; }
//   Int_t RowId( void ) const { return m_row_id; }
// };

//_____________________________________________________________________________
class AftHelper
{
public:
  static std::string&  ClassName( void );
  static AftHelper& GetInstance( void );
  ~AftHelper( void );

private:
  AftHelper( void );
  AftHelper( const AftHelper& );
  AftHelper& operator=( const AftHelper );

private:
  // typedef std::map<Int_t,AftParam*> TpcGeomMap;
  // TpcGeomMap m_map;
  const Double_t m_distance_x_xp = 2.7; // distance between X  and X' layer (Y  and Y' layer)
  const Double_t m_distance_xp_y = 3.1; // distance between X' and Y  layer (Y' and X  layer)
  const Double_t m_distance_seg  = 3.1; // frame phi for a fiber [mm]
  const Double_t m_posZ0         = m_distance_seg/2.;
  const Double_t m_posX0Seg0     = m_distance_seg/2. - m_distance_seg*NumOfSegAFTX/2.;
  const Double_t m_posX1Seg0     = -m_distance_seg*NumOfSegAFTX/2.;
  const Double_t m_posY0Seg0     = m_distance_seg/2. - m_distance_seg*NumOfSegAFTY/2.;
  const Double_t m_posY1Seg0     = -m_distance_seg*NumOfSegAFTY/2.;
  const Double_t m_posSeg0[4]    = {m_posX0Seg0, m_posX1Seg0, m_posY0Seg0, m_posY1Seg0};
  const Double_t m_phi = 3.; // diameter of a fiber [mm]
  const Int_t    m_npoly = 360;

public:
  Double_t     GetX( Int_t plane, Int_t seg );
  Double_t     GetZ( Int_t plane, Int_t seg );
  Double_t     GetPhi( void ) const { return m_phi; }
  Double_t     GetNPoly( void ) const { return m_npoly; }
  // AftParam* GetParam( Int_t asad, Int_t aget, Int_t ch ) const;
  // AftParam* GetParam( Int_t layer, Int_t row ) const;
  // AftParam* GetParam( Int_t pad ) const;
  // TVector3     GetPoint( Int_t pad ) const;
  // Int_t        GetPadId( Int_t asad, Int_t aget, Int_t ch ) const;
  // Int_t        GetPadId( Int_t layer, Int_t row ) const;
  // Int_t        GetLayerId( Int_t asad, Int_t aget, Int_t ch ) const;
  // Int_t        GetRowId( Int_t asad, Int_t aget, Int_t ch ) const;
  // Bool_t       Initialize( const TString& file_name );

};

//_____________________________________________________________________________
inline std::string&
AftHelper::ClassName( void )
{
  static std::string s_name( "AftHelper" );
  return s_name;
}

//_____________________________________________________________________________
inline AftHelper&
AftHelper::GetInstance( void )
{
  static AftHelper s_instance;
  return s_instance;
}

#endif

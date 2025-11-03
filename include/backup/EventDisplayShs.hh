// -*- C++ -*-

#ifndef EVENT_DISPLAY_SHS_HH
#define EVENT_DISPLAY_SHS_HH

#include "ThreeVector.hh"
#include "DetectorID.hh"

class DCLocalTrack;

class TApplication;
class TArc;
class TBox;
class TBRIK;
class TCanvas;
class TFile;
class TGeometry;
class TGraph;
class TGraphErrors;
class TH1;
class TH2;
class TH2Poly;
class TF1;
class TLatex;
class TMarker;
class TMarker3DBox;
class TMixture;
class TNode;
class TPad;
class TPave;
class TPaveLabel;
class TPaveText;
class TPGON;
class TPolyLine;
class TPolyLine3D;
class TPolyMarker;
class TPolyMarker3D;
class TROOT;
class TRint;
class TRotMatrix;
class TStyle;
class TTRD1;
class TTRD2;
class TTUBS;
class TView;
class TLine;

typedef std::vector<TLine*> TLineContainer;

//______________________________________________________________________________
class EventDisplayShs
{
public:
  static EventDisplayShs& GetInstance( void );
  static const TString& ClassName( void );
  ~EventDisplayShs( void );

private:
  EventDisplayShs( void );
  EventDisplayShs( const EventDisplayShs& );
  EventDisplayShs& operator =( const EventDisplayShs& );

private:
  Bool_t                     m_is_ready;
  TApplication              *m_theApp;
  TGeometry                 *m_geometry;
  TNode                     *m_node;
  TCanvas                   *m_canvas;
  TH1*                       m_tpc_adc;
  TH1*                       m_tpc_tdc;
  TH2Poly*                   m_tpc_adc2d;
  TH2Poly*                   m_tpc_tdc2d;
  TH2Poly*                   m_htof_2d;
  TPGON*                     m_tpc_pgon;
  TNode*                     m_tpc_node;
  TPolyMarker3D*             m_tpc_mark3d;

public:
  void   CalcRotMatrix(Double_t TA, Double_t RA1, Double_t RA2, Double_t *rotMat);
  void   EndOfEvent(void);
  void   FillTPCADC(Int_t layer, Int_t row, Double_t adc);
  void   FillTPCTDC(Int_t layer, Int_t row, Double_t tdc);
  void   FillHTOF(Int_t seg);
  Int_t  GetCommand(void) const;
  Bool_t Initialize(void);
  Bool_t IsReady(void) const { return m_is_ready; }
  void   Reset(void);
  void   Run(Bool_t flag=true);
  void   SetTPCMarker(Double_t x, Double_t y, Double_t z);
  void   SetTPCMarker(const TVector3& pos);
  void   Update(void);
};

//______________________________________________________________________________
inline EventDisplayShs&
EventDisplayShs::GetInstance( void )
{
  static EventDisplayShs s_instance;
  return s_instance;
}

//______________________________________________________________________________
inline const TString&
EventDisplayShs::ClassName( void )
{
  static TString s_name("EventDisplayShs");
  return s_name;
}

#endif

// -*- C++ -*-

#ifndef EVENT_DISPLAY_HH
#define EVENT_DISPLAY_HH

#include <TROOT.h>

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
class TPolyLine;
class TPolyLine3D;
class TPolyMarker;
class TPolyMarker3D;
class TRint;
class TRotMatrix;
class TStyle;
class TTRD1;
class TTRD2;
class TTUBS;
class TView;
class TLine;

typedef std::vector <TLine*> TLineContainer;

//_____________________________________________________________________________
class EventDisplay //: public TObject
{
public:
  static const TString& ClassName();
  static EventDisplay& GetInstance();
  ~EventDisplay();

private:
  EventDisplay();
  EventDisplay(const EventDisplay&);
  EventDisplay& operator =(const EventDisplay&);

private:
  Bool_t                     m_is_ready;
  Bool_t                     m_is_save_mode;
  TApplication              *m_theApp;
  TGeometry                 *m_geometry;
  TNode                     *m_node;
  TCanvas                   *m_canvas;
  TCanvas                   *m_canvas_vertex;
  TCanvas                   *m_canvas_hist;
  TCanvas                   *m_canvas_hist2;
  TCanvas                   *m_canvas_hist3;
  TCanvas                   *m_canvas_hist4;
  TCanvas                   *m_canvas_hist5;
  TCanvas                   *m_canvas_hist6;
  TCanvas                   *m_canvas_hist7;
  TCanvas                   *m_canvas_hist8;
  TCanvas                   *m_canvas_hist9;
  TH1                       *m_hist_vertex_x;
  TH1                       *m_hist_vertex_y;
  TH1                       *m_hist_p;
  TH1                       *m_hist_m2;
  TH1                       *m_hist_missmass;
  TH2                       *m_hist_bh1;
  TH2                       *m_hist_bft;
  TH2                       *m_hist_bft_p;
  TH2                       *m_hist_bcIn;
  TH2                       *m_hist_bcOut;
  std::vector<TBox*>         m_BH1box_cont;
  std::vector<TLine*>        m_BcInTrack;
  std::vector<TBox*>         m_BH2box_cont;
  std::vector<TLine*>        m_BcOutTrack2;
  TH2                       *m_hist_bh2;
  TH2                       *m_hist_bcOut_sdcIn;
  TH2                       *m_hist_sdcIn_predict;
  TH2                       *m_hist_sdcIn_predict2;
  std::vector<TBox*>         m_BH2box_cont2;
  std::vector<TBox*>         m_SCHbox_cont;
  TBox                      *m_TargetXZ_box2;
  TBox                      *m_TargetYZ_box2;
  std::vector<TLine*>        m_BcOutTrack3;
  std::vector<TLine*>        m_SdcInTrack2;
  TH2                       *m_hist_sch;
  TH2                       *m_hist_tof;
  TH2                       *m_hist_sdc1;
  TH2                       *m_hist_sdc1p;
  TH2                       *m_hist_sdc3_l;
  TH2                       *m_hist_sdc3_t;
  TH2                       *m_hist_sdc3p_l;
  TH2                       *m_hist_sdc3p_t;
  TH2                       *m_hist_sdc3y_l;
  TH2                       *m_hist_sdc3y_t;
  TH2                       *m_hist_sdc3yp_l;
  TH2                       *m_hist_sdc3yp_t;
  TH2                       *m_hist_sdc4_l;
  TH2                       *m_hist_sdc4_t;
  TH2                       *m_hist_sdc4p_l;
  TH2                       *m_hist_sdc4p_t;
  TH2                       *m_hist_sdc4y_l;
  TH2                       *m_hist_sdc4y_t;
  TH2                       *m_hist_sdc4yp_l;
  TH2                       *m_hist_sdc4yp_t;

  TH2                       *m_hist_bc3;
  TH2                       *m_hist_bc3p;
  TH2                       *m_hist_bc3u;
  TH2                       *m_hist_bc3up;
  TH2                       *m_hist_bc3v;
  TH2                       *m_hist_bc3vp;
  TH2                       *m_hist_bc4;
  TH2                       *m_hist_bc4p;
  TH2                       *m_hist_bc4u;
  TH2                       *m_hist_bc4up;
  TH2                       *m_hist_bc4v;
  TH2                       *m_hist_bc4vp;
  TH1                       *m_hist_bc3_time;
  TH1                       *m_hist_bc3p_time;
  TH1                       *m_hist_bc4_time;
  TH1                       *m_hist_bc4p_time;

  TNode                     *m_target_node;
  TNode                     *m_kurama_inner_node;
  TNode                     *m_kurama_outer_node;
  TNode                     *m_BH2wall_node;
  std::vector<TNode*>        m_BH2seg_node;
  std::vector<TNode*>        m_BC3x1_node;
  std::vector<TNode*>        m_BC3x2_node;
  std::vector<TNode*>        m_BC3v1_node;
  std::vector<TNode*>        m_BC3v2_node;
  std::vector<TNode*>        m_BC3u1_node;
  std::vector<TNode*>        m_BC3u2_node;
  std::vector<TNode*>        m_BC4u1_node;
  std::vector<TNode*>        m_BC4u2_node;
  std::vector<TNode*>        m_BC4v1_node;
  std::vector<TNode*>        m_BC4v2_node;
  std::vector<TNode*>        m_BC4x1_node;
  std::vector<TNode*>        m_BC4x2_node;
  std::vector<TNode*>        m_SDC1v1_node;
  std::vector<TNode*>        m_SDC1v2_node;
  std::vector<TNode*>        m_SDC1x1_node;
  std::vector<TNode*>        m_SDC1x2_node;
  std::vector<TNode*>        m_SDC1u1_node;
  std::vector<TNode*>        m_SDC1u2_node;
  std::vector<TNode*>        m_SDC2x1_node;
  std::vector<TNode*>        m_SDC2x2_node;
  std::vector<TNode*>        m_SDC2y1_node;
  std::vector<TNode*>        m_SDC2y2_node;
  std::vector<TNode*>        m_SDC3x1_node;
  std::vector<TNode*>        m_SDC3x2_node;
  std::vector<TNode*>        m_SDC3y1_node;
  std::vector<TNode*>        m_SDC3y2_node;
  std::vector<TNode*>        m_SDC4y1_node;
  std::vector<TNode*>        m_SDC4y2_node;
  std::vector<TNode*>        m_SDC4x1_node;
  std::vector<TNode*>        m_SDC4x2_node;
  std::vector<TNode*>        m_SSD1y1_node;
  TNode                     *m_FBHwall_node;
  std::vector<TNode*>        m_FBHseg_node;
  TNode                     *m_SCHwall_node;
  std::vector<TNode*>        m_SCHseg_node;
  TNode                     *m_TOFwall_node;
  std::vector<TNode*>        m_TOFseg_node;
  TNode                     *m_WCwall_node;
  std::vector<TNode*>        m_WCseg_node;
  std::vector<TPolyLine3D*>  m_BcOutTrack;
  // TPolyMarker*               m_BcOutTrackShs;
  TF1*                       m_BcOutTrackShs;
  std::vector<TPolyLine3D*>  m_SdcInTrack;
  std::vector<TPolyLine3D*>  m_SdcOutTrack;
  TPolyMarker3D             *m_init_step_mark;
  TPolyMarker3D             *m_hs_step_mark;
  TPolyMarker3D             *m_kurama_step_mark;
  // vertex
  TBox                      *m_TargetXZ_box;
  TBox                      *m_TargetYZ_box;
  TBox                      *m_EmulsionXZ_box;
  TBox                      *m_EmulsionYZ_box;
  std::vector<TBox*>         m_FBHseg_box;
  TH1                       *m_SSD_x_hist;
  TH1                       *m_SSD_y_hist;
  TMarker                   *m_VertexPointXZ;
  TMarker                   *m_VertexPointYZ;
  std::vector<TPolyLine*>    m_BcOutXZ_line;
  std::vector<TPolyLine*>    m_BcOutYZ_line;
  std::vector<TPolyLine*>    m_SdcInXZ_line;
  std::vector<TPolyLine*>    m_SdcInYZ_line;
  TPolyMarker               *m_HSMarkVertexXShs;
  TPolyMarker               *m_KuramaMarkVertexXShs;
  TPolyMarker               *m_KuramaMarkVertexX;
  TPolyMarker               *m_KuramaMarkVertexY;
  TPolyLine                 *m_MissMomXZ_line;
  TPolyLine                 *m_MissMomYZ_line;

public:
  Bool_t Initialize();
  Bool_t IsReady() const { return m_is_ready; }
  Bool_t ConstructEmulsion();
  Bool_t ConstructTarget();
  Bool_t ConstructBH2();
  Bool_t ConstructKURAMA();
  Bool_t ConstructCollimator();
  Bool_t ConstructBcOut();
  Bool_t ConstructSdcIn();
  Bool_t ConstructSdcOut();
  Bool_t ConstructSSD();
  Bool_t ConstructFBH();
  Bool_t ConstructSCH();
  Bool_t ConstructTOF();
  Bool_t ConstructWC();
  void DrawInitTrack(Int_t nStep, ThreeVector *StepPoint);
  void DrawInitTrack();
  void DrawHitWire(Int_t lid, Int_t hit_wire,
                   Bool_t range_check=true, Bool_t tdc_check=true);
  void DrawRunEvent(Double_t xpos, Double_t ypos, const TString& arg);
  void DrawText(Double_t xpos, Double_t ypos, const TString& arg);
  void DrawTrackWire(Int_t lid, Int_t hit_wire, Int_t it);
  void DrawHitHodoscope(Int_t lid, Int_t seg, Int_t Tu=1, Int_t Td=1);
  void DrawBcOutLocalTrack(DCLocalTrack *tp);
  void DrawBcOutLocalTrack(Double_t x0, Double_t y0, Double_t u0, Double_t v0);
  void DrawSdcInLocalTrack(DCLocalTrack *tp);
  void DrawSdcOutLocalTrack(DCLocalTrack *tp);
  void DrawSsdHit(Int_t lid, Int_t seg, Double_t de);
  void DrawVertex(const ThreeVector& vertex);
  void DrawKuramaTrack(Int_t nStep, const std::vector<TVector3>& StepPoint,
                       Double_t q);
  void DrawHSTrack(Int_t nStep, const std::vector<TVector3>& StepPoint,
                       Double_t q);
  void DrawTarget();
  void DrawMissingMomentum(const ThreeVector& mom,
                           const ThreeVector& pos);
  void FillMomentum(Double_t momentum);
  void FillMassSquare(Double_t mass_square);
  void FillMissMass(Double_t mass_square);
  void FillBH1(Int_t seg, Int_t tdc);
  void SetCorrectTimeBH1(Int_t seg, Double_t de);
  void FillBFT(Int_t layer, Int_t seg, Int_t tdc);
  void SetCorrectTimeBFT(Double_t pos);
  void DrawBcInTrack(Double_t x0, Double_t u0);
  void FillBH2(Int_t seg, Int_t tdc);
  void SetCorrectTimeBH2(Int_t seg, Double_t de);
  void SetCorrectTimeBcOut(Int_t layer, Double_t pos);
  void DrawBcOutTrack(Double_t x0, Double_t u0, Double_t y0, Double_t v0, Bool_t flagGoodForTracking = true);
  void DrawSdcInTrack(Double_t x0, Double_t u0, Double_t y0, Double_t v0, Bool_t flagKurama = false, Bool_t flagBeam = false);
  void FillSCH(Int_t seg, Int_t tdc);
  void SetCorrectTimeSCH(Int_t seg, Double_t de);
  void SetCorrectTimeSdcIn(Int_t layer, Double_t pos);

  void FillTOF(Int_t seg, Int_t tdc);
  void FillBcOutHit(Int_t layer,  Int_t wire, Int_t tdc);
  void DrawBC3(Int_t wire, Int_t tdc);
  void DrawBC3p(Int_t wire, Int_t tdc);
  void DrawBC4(Int_t wire, Int_t tdc);
  void DrawBC4p(Int_t wire, Int_t tdc);
  void FillSDC1(Int_t wire, Int_t tdc);
  void FillSDC1p(Int_t wire, Int_t tdc);
  void FillSdcOutHit(Int_t layer,  Int_t wire, Int_t LorT, Int_t tdc);
  void FillSDC3_Leading(Int_t wire, Int_t tdc);
  void FillSDC3_Trailing(Int_t wire, Int_t tdc);
  void FillSDC3p_Leading(Int_t wire, Int_t tdc);
  void FillSDC3p_Trailing(Int_t wire, Int_t tdc);
  void FillSDC4_Leading(Int_t wire, Int_t tdc);
  void FillSDC4_Trailing(Int_t wire, Int_t tdc);
  void FillSDC4p_Leading(Int_t wire, Int_t tdc);
  void FillSDC4p_Trailing(Int_t wire, Int_t tdc);
  void ShowHitFiber(Int_t layer, Int_t segment, Double_t pe);// const;
  void Update();
  void EndOfEvent();
  void ResetVisibility();
  void CalcRotMatrix(Double_t TA, Double_t RA1, Double_t RA2, Double_t *rotMat);
  Int_t GetCommand();
  void Print(Int_t run_number, Int_t event_number);
  void Run(Bool_t flag=kTRUE);
  void SetSaveMode(){ m_is_save_mode = true; gROOT->SetBatch(); }

private:
  void ResetVisibility(TNode *& node, Color_t c=kWhite);
  void ResetVisibility(std::vector<TNode*>& node, Color_t c=kWhite);
  void ResetHist();
};

//_____________________________________________________________________________
inline EventDisplay&
EventDisplay::GetInstance()
{
  static EventDisplay s_instance;
  return s_instance;
}

//_____________________________________________________________________________
inline const TString&
EventDisplay::ClassName()
{
  static TString s_name("EventDisplay");
  return s_name;
}

#endif

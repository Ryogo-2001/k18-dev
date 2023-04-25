/*
  EvDispCFT.hh

  2012/1/24
*/

#ifndef EvDispCFT_h
#define EvDispCFT_h 1

#include "TROOT.h"
#include "TApplication.h"
#include "TRint.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGeometry.h"
#include "TMixture.h"
#include "TBRIK.h"
#include "TTRD1.h"
#include "TTRD2.h"
#include "TTUBS.h"
#include "TTUBS.h"
#include "TRotMatrix.h"
#include "TNode.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TView.h"
#include "TPad.h"
#include "TButton.h"
#include "TMarker3DBox.h"
#include "TPave.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TArc.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TBox.h"
#include "TStyle.h"

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include "ThreeVector.hh"
#include "DetectorID.hh"

#include <vector>

extern  TApplication *theApp;


class EvDispCFT : public TObject {
private:
  TPad      *tp_[10];
  TCanvas   *tc_;
  TH2F      *hbase_;
  TH2F      *hbaseU_;
  TH1F      *hp_[NumOfPlaneCFT];
  TH2F      *hbaseZX_;
  TH2F      *hbaseZY_;
  
  TArc      *Phi1_Arc_[NumOfSegCFT_PHI4];
  TArc      *Phi2_Arc_[NumOfSegCFT_PHI4];
  TArc      *Phi3_Arc_[NumOfSegCFT_PHI4];
  TArc      *Phi4_Arc_[NumOfSegCFT_PHI4];
  TArc      *UV1_Arc_[NumOfSegCFT_PHI4];
  TArc      *UV2_Arc_[NumOfSegCFT_PHI4];
  TArc      *UV3_Arc_[NumOfSegCFT_PHI4];
  TArc      *UV4_Arc_[NumOfSegCFT_PHI4];

  TArc      *U_Arc_[NumOfSegCFT_PHI4];
  TLine     *U_Line_[NumOfSegCFT_PHI4][2];

  TBox      *BGO_Box_[24];
  TLine     *BGO_Line_[24][4];

  TBox      *PiID_Box_[32];
  TLine     *PiID_Line_[32][4];

  TBox      *Tgt_Box_;
  TArc      *Tgt_Arc_;
  TArc      *TgtHead_Arc_;

  TBox      *CFRP_Box_;
  TArc      *CFRP_Arc_;

  mutable std::vector <TLine*> TrackXYCont;
  mutable std::vector <TLine*> TrackZXCont;
  mutable std::vector <TLine*> TrackZYCont;

public:
  EvDispCFT(void);
  ~EvDispCFT(void);
  static EvDispCFT & GetInstance( void );
  void Initialize(int RunNum);
  void ConstructCFT(void);
  void ShowHitFiber(int layer, int segment, double pe) const;
  void ShowHitBGO(int seg, int ADC) const;
  void ShowHitPiID(int seg) const;
  void ShowHitPos(double X, double Y, double pe) const;
  void ShowHitPosZX(double Z, double X, double pe) const;
  void ShowHitPosZY(double Z, double Y, double pe) const;
  void DrawTrackInXYPlane(double x0, double y0, double x1, double y1) const;
  void DrawTrackInZXPlane(double z0, double x0, double z1, double x1) const;
  void DrawTrackInZYPlane(double z0, double y0, double z1, double y1) const;
  void DrawTrackInXYPlane_(double x0, double y0, double x1, double y1) const;
  void DrawTrackInZXPlane_(double z0, double x0, double z1, double x1) const;
  void DrawTrackInZYPlane_(double z0, double y0, double z1, double y1) const;
  void DrawTrackInZXPlane__(double z0, double x0, double z1, double x1) const;
  void DrawTrackInZYPlane__(double z0, double y0, double z1, double y1) const;
  void UpdateCanvas() const;
  void EndOfEvent() const;
  void FiberPosPhi(int layer, int seg, double *x, double *y) const;
  void FiberPosU(int layer, int seg, double z, double *x, double *y) const;
  void get_command(void) const;

private:
  static EvDispCFT *evDisp_;
  static TApplication *theApp;
};


#endif

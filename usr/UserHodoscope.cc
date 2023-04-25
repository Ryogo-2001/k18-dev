// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iostream>
#include <sstream>

#include <UnpackerManager.hh>

#include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "KuramaLib.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "DCGeomMan.hh"

#define TimeCut    1 // in cluster analysis
#define FHitBranch 0 // make FiberHit branches (becomes heavy)
#define HodoHitPos 0

namespace
{
using namespace root;
const auto qnan = TMath::QuietNaN();
auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
auto& gRM = RMAnalyzer::GetInstance();
auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
class UserHodoscope : public VEvent
{
private:
  RawData*      rawData;
  HodoAnalyzer* hodoAna;

public:
  UserHodoscope();
  ~UserHodoscope();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserHodoscope::ClassName()
{
  static TString s_name("UserHodoscope");
  return s_name;
}

//_____________________________________________________________________________
UserHodoscope::UserHodoscope()
  : VEvent(),
    rawData(new RawData),
    hodoAna(new HodoAnalyzer)
{
}

//_____________________________________________________________________________
UserHodoscope::~UserHodoscope()
{
  if(hodoAna) delete hodoAna;
  if(rawData) delete rawData;
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;
  Int_t spill;

  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  Int_t bh1nhits;
  Int_t bh1hitpat[MaxHits];
  Double_t bh1ua[NumOfSegBH1];
  Double_t bh1ut[NumOfSegBH1][MaxDepth];
  Double_t bh1da[NumOfSegBH1];
  Double_t bh1dt[NumOfSegBH1][MaxDepth];

  Int_t bh2nhits;
  Int_t bh2hitpat[MaxHits];
  Double_t bh2ua[NumOfSegBH2];
  Double_t bh2ut[NumOfSegBH2][MaxDepth];
  Double_t bh2da[NumOfSegBH2];
  Double_t bh2dt[NumOfSegBH2][MaxDepth];

  Int_t bacnhits;
  Int_t bachitpat[MaxHits];
  Double_t baca[NumOfSegBAC];
  Double_t bact[NumOfSegBAC][MaxDepth];

  Int_t htofnhits;
  Int_t htofhitpat[MaxHits];
  Double_t htofua[NumOfSegHTOF];
  Double_t htofda[NumOfSegHTOF];
  Double_t htofut[NumOfSegHTOF][MaxDepth];
  Double_t htofdt[NumOfSegHTOF][MaxDepth];

  Int_t bvhnhits;
  Int_t bvhhitpat[MaxHits];
  Double_t bvht[NumOfSegBVH][MaxDepth];

  Int_t tofnhits;
  Int_t tofhitpat[MaxHits];
  Int_t tofnhits_3dmtx;
  Int_t tofhitpat_3dmtx[MaxHits];
  Double_t tofua[NumOfSegTOF];
  Double_t tofut[NumOfSegTOF][MaxDepth];
  Double_t tofda[NumOfSegTOF];
  Double_t tofdt[NumOfSegTOF][MaxDepth];

  Int_t lacnhits;
  Int_t lachitpat[MaxHits];
  Double_t lact[NumOfSegLAC][MaxDepth];

  Int_t wcnhits;
  Int_t wchitpat[MaxHits];
  Double_t wcua[NumOfSegWC];
  Double_t wcut[NumOfSegWC][MaxDepth];
  Double_t wcda[NumOfSegWC];
  Double_t wcdt[NumOfSegWC][MaxDepth];
  Int_t wcsumnhits;
  Int_t wcsumhitpat[MaxHits];
  Double_t wcsuma[NumOfSegWC];
  Double_t wcsumt[NumOfSegWC][MaxDepth];

  // SCH
  int    sch_nhits;
  int    sch_hitpat[NumOfSegSCH];
  double sch_tdc[NumOfSegSCH][MaxDepth];
  double sch_trailing[NumOfSegSCH][MaxDepth];
  double sch_tot[NumOfSegSCH][MaxDepth];
  int    sch_depth[NumOfSegSCH];
  int    sch_ncl;
  int    sch_clsize[NumOfSegSCH];
  double sch_ctime[NumOfSegSCH];
  double sch_ctot[NumOfSegSCH];
  double sch_clpos[NumOfSegSCH];


  ////////// Normalized
  Double_t bh1mt[NumOfSegBH1][MaxDepth];
  Double_t bh1cmt[NumOfSegBH1][MaxDepth];
  Double_t bh1utime[NumOfSegBH1][MaxDepth];
  Double_t bh1uctime[NumOfSegBH1][MaxDepth];
  Double_t bh1dtime[NumOfSegBH1][MaxDepth];
  Double_t bh1dctime[NumOfSegBH1][MaxDepth];
  Double_t bh1de[NumOfSegBH1];
  Double_t bh1ude[NumOfSegBH1];
  Double_t bh1dde[NumOfSegBH1];

  Double_t bh2mt[NumOfSegBH2][MaxDepth];
  Double_t bh2cmt[NumOfSegBH2][MaxDepth];
  Double_t bh2utime[NumOfSegBH2][MaxDepth];
  Double_t bh2uctime[NumOfSegBH2][MaxDepth];
  Double_t bh2dtime[NumOfSegBH2][MaxDepth];
  Double_t bh2dctime[NumOfSegBH2][MaxDepth];
  Double_t bh2hitpos[NumOfSegBH2][MaxDepth];
  Double_t bh2de[NumOfSegBH2];
  Double_t bh2ude[NumOfSegBH2];
  Double_t bh2dde[NumOfSegBH2];

  Double_t bacmt[NumOfSegBAC][MaxDepth];
  Double_t bacde[NumOfSegBAC];

  Double_t t0[NumOfSegBH2][MaxDepth];
  Double_t ct0[NumOfSegBH2][MaxDepth];
  Double_t btof[NumOfSegBH1][NumOfSegBH2];
  Double_t cbtof[NumOfSegBH1][NumOfSegBH2];

  Double_t tofmt[NumOfSegTOF][MaxDepth];
  Double_t tofde[NumOfSegTOF];

  Double_t htofmt[NumOfSegHTOF][MaxDepth];
  Double_t htofde[NumOfSegHTOF];
  Double_t htofutime[NumOfSegHTOF][MaxDepth];
  Double_t htofuctime[NumOfSegHTOF][MaxDepth];
  Double_t htofdtime[NumOfSegHTOF][MaxDepth];
  Double_t htofdctime[NumOfSegHTOF][MaxDepth];
  Double_t htofhitpos[NumOfSegHTOF][MaxDepth];
  Double_t htofude[NumOfSegHTOF];
  Double_t htofdde[NumOfSegHTOF];

  // Using WCSUM
  Double_t wcde[NumOfSegWC];
  Double_t wcmt[NumOfSegWC][MaxDepth];

  // Time0
  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;

  // Btof0 BH1
  Double_t Btof0Seg;
  Double_t deBtof0;
  Double_t Btof0;
  Double_t CBtof0;

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum      = 0;
  spill      = 0;
  bh1nhits   = 0;
  bacnhits   = 0;
  bh2nhits   = 0;
  bvhnhits   = 0;
  htofnhits  = 0;
  tofnhits   = 0;
  tofnhits_3dmtx   = 0;
  lacnhits   = 0;
  sch_nhits  = 0;
  sch_ncl  = 0;
  wcnhits    = 0;
  wcsumnhits = 0;
  Time0Seg = qnan;
  deTime0  = qnan;
  Time0    = qnan;
  CTime0   = qnan;
  Btof0Seg = qnan;
  deBtof0  = qnan;
  Btof0    = qnan;
  CBtof0   = qnan;

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<MaxHits; ++it){
    bh1hitpat[it]   = -1;
    bh2hitpat[it]   = -1;
    bachitpat[it]   = -1;
    htofhitpat[it]  = -1;
    bvhhitpat[it]   = -1;
    tofhitpat[it]   = -1;
    tofhitpat_3dmtx[it]   = -1;
    lachitpat[it]   = -1;
    wchitpat[it]    = -1;
    wcsumhitpat[it] = -1;
  }

  for(Int_t it=0; it<NumOfSegBH1; ++it){
    bh1ua[it] = qnan;
    bh1da[it] = qnan;
    bh1de[it] = qnan;
    bh1ude[it] = qnan;
    bh1dde[it] = qnan;
    for(Int_t that=0; that<NumOfSegBH2; ++that){
      btof[it][that]  = qnan;
      cbtof[it][that] = qnan;
    }
    for(Int_t m = 0; m<MaxDepth; ++m){
      bh1ut[it][m] = qnan;
      bh1dt[it][m] = qnan;
      bh1mt[it][m] = qnan;
      bh1utime[it][m] = qnan;
      bh1dtime[it][m] = qnan;
      bh1uctime[it][m] = qnan;
      bh1dctime[it][m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegBH2; ++it){
    bh2ua[it] = qnan;
    bh2da[it] = qnan;
    bh2de[it] = qnan;
    bh2dde[it] = qnan;
    bh2ude[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      bh2ut[it][m] = qnan;
      bh2dt[it][m] = qnan;
      bh2mt[it][m] = qnan;
      t0[it][m]    = qnan;
      ct0[it][m]   = qnan;
      bh2utime[it][m] = qnan;
      bh2dtime[it][m] = qnan;
      bh2uctime[it][m] = qnan;
      bh2dctime[it][m] = qnan;
      bh2hitpos[it][m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegBAC; it++){
    baca[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      bact[it][m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegHTOF; it++){
    htofua[it] = qnan;
    htofda[it] = qnan;
    htofde[it] = qnan;
    htofude[it] = qnan;
    htofdde[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      htofut[it][m] = qnan;
      htofdt[it][m] = qnan;
      htofmt[it][m] = qnan;
      htofutime[it][m] = qnan;
      htofdtime[it][m] = qnan;
      htofuctime[it][m] = qnan;
      htofdctime[it][m] = qnan;
      htofhitpos[it][m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegBVH; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      bvht[it][m]  = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegTOF; it++){
    tofua[it] = qnan;
    tofda[it] = qnan;
    tofde[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      tofut[it][m] = qnan;
      tofdt[it][m] = qnan;
      tofmt[it][m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegLAC; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      lact[it][m]  = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegWC; ++it){
    wcua[it]   = qnan;
    wcda[it]   = qnan;
    wcsuma[it] = qnan;
    wcde[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      wcut[it][m]   = qnan;
      wcdt[it][m]   = qnan;
      wcsumt[it][m] = qnan;
      wcmt[it][m]   = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegSCH; ++it){
    sch_hitpat[it]   = qnan;
    sch_depth[it]   = qnan;
    sch_clsize[it]   = qnan;
    sch_ctime[it]   = qnan;
    sch_ctot[it]   = qnan;
    sch_clpos[it]   = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      sch_tdc[it][m]   = qnan;
      sch_trailing[it][m] = qnan;
      sch_tot[it][m]   = qnan;
    }
  }


}

//_____________________________________________________________________________
struct Dst
{
  Int_t evnum;
  Int_t spill;

  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  Int_t    nhBh1;
  Int_t    csBh1[NumOfSegBH1*MaxDepth];
  Double_t Bh1Seg[NumOfSegBH1*MaxDepth];
  Double_t tBh1[NumOfSegBH1*MaxDepth];
  Double_t dtBh1[NumOfSegBH1*MaxDepth];
  Double_t deBh1[NumOfSegBH1*MaxDepth];
  Double_t btof[NumOfSegBH1*MaxDepth];
  Double_t cbtof[NumOfSegBH1*MaxDepth];

  // Btof0 BH1
  Double_t Btof0Seg;
  Double_t deBtof0;
  Double_t Btof0;
  Double_t CBtof0;

  Int_t    nhBh2;
  Int_t    csBh2[NumOfSegBH2*MaxDepth];
  Double_t Bh2Seg[NumOfSegBH2*MaxDepth];
  Double_t tBh2[NumOfSegBH2*MaxDepth];
  Double_t t0Bh2[NumOfSegBH2*MaxDepth];
  Double_t dtBh2[NumOfSegBH2*MaxDepth];
  Double_t deBh2[NumOfSegBH2*MaxDepth];
  Double_t posBh2[NumOfSegBH2*MaxDepth];

  Int_t    nhBac;
  Int_t    csBac[NumOfSegBAC*MaxDepth];
  Double_t BacSeg[NumOfSegBAC*MaxDepth];
  Double_t tBac[NumOfSegBAC*MaxDepth];
  Double_t deBac[NumOfSegBAC*MaxDepth];

  // Time0
  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;

  Int_t    nhTof;
  Int_t    csTof[NumOfSegTOF*MaxDepth];
  Double_t TofSeg[NumOfSegTOF*MaxDepth];
  Double_t tTof[NumOfSegTOF*MaxDepth];
  Double_t dtTof[NumOfSegTOF*MaxDepth];
  Double_t deTof[NumOfSegTOF*MaxDepth];

  Int_t    nhHtof;
  Int_t    csHtof[NumOfSegHTOF*MaxDepth];
  Double_t HtofSeg[NumOfSegHTOF*MaxDepth];
  Double_t tHtof[NumOfSegHTOF*MaxDepth];
  Double_t dtHtof[NumOfSegHTOF*MaxDepth];
  Double_t deHtof[NumOfSegHTOF*MaxDepth];
  Double_t posHtof[NumOfSegHTOF*MaxDepth];

  Int_t    nhBvh;
  Int_t    csBvh[NumOfSegBVH*MaxDepth];
  Double_t BvhSeg[NumOfSegBVH*MaxDepth];
  Double_t tBvh[NumOfSegBVH*MaxDepth];

  Int_t    nhLac;
  Int_t    csLac[NumOfSegLAC*MaxDepth];
  Double_t LacSeg[NumOfSegLAC*MaxDepth];
  Double_t tLac[NumOfSegLAC*MaxDepth];

  // for HodoParam
  Double_t tofua[NumOfSegTOF];
  Double_t tofut[NumOfSegTOF][MaxDepth];
  Double_t tofda[NumOfSegTOF];
  Double_t tofdt[NumOfSegTOF][MaxDepth];
  Double_t utTofSeg[NumOfSegTOF][MaxDepth];
  Double_t dtTofSeg[NumOfSegTOF][MaxDepth];
  Double_t udeTofSeg[NumOfSegTOF];
  Double_t ddeTofSeg[NumOfSegTOF];

  void clear();
};

//_____________________________________________________________________________
void
Dst::clear()
{
  nhBh1    = 0;
  nhBac    = 0;
  nhBh2    = 0;
  nhHtof   = 0;
  nhTof    = 0;
  nhBvh    = 0;
  nhLac    = 0;
  evnum    = 0;
  spill    = 0;
  Time0Seg = qnan;
  deTime0  = qnan;
  Time0    = qnan;
  CTime0   = qnan;
  Btof0Seg = qnan;
  deBtof0  = qnan;
  Btof0    = qnan;
  CBtof0   = qnan;

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<NumOfSegBH1; ++it){
    for(Int_t m=0; m<MaxDepth; ++m){
      csBh1[MaxDepth*it + m]  = 0;
      Bh1Seg[MaxDepth*it + m] = qnan;
      tBh1[MaxDepth*it + m]   = qnan;
      dtBh1[MaxDepth*it + m]  = qnan;
      deBh1[MaxDepth*it + m]  = qnan;
      btof[MaxDepth*it + m]   = qnan;
      cbtof[MaxDepth*it + m]  = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegBH2; ++it){
    for(Int_t m=0; m<MaxDepth; ++m){
      csBh2[MaxDepth*it + m]  = 0;
      Bh2Seg[MaxDepth*it + m] = qnan;
      tBh2[MaxDepth*it + m]   = qnan;
      t0Bh2[MaxDepth*it + m]  = qnan;
      dtBh2[MaxDepth*it + m]  = qnan;
      deBh2[MaxDepth*it + m]  = qnan;
      posBh2[MaxDepth*it + m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegBAC; it++){
    BacSeg[it] = qnan;
    deBac[it]  = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      tBac[MaxDepth*it+m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegHTOF; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      csHtof[MaxDepth*it + m]  = 0;
      HtofSeg[MaxDepth*it + m] = qnan;
      tHtof[MaxDepth*it + m]   = qnan;
      dtHtof[MaxDepth*it + m]  = qnan;
      deHtof[MaxDepth*it + m]  = qnan;
      posHtof[MaxDepth*it + m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegBVH; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      csBvh[MaxDepth*it + m]  = 0;
      BvhSeg[MaxDepth*it + m] = qnan;
      tBvh[MaxDepth*it + m]   = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegTOF; it++){
    tofua[it] = qnan;
    tofda[it] = qnan;
    udeTofSeg[it] = qnan;
    ddeTofSeg[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      tofut[it][m] = qnan;
      tofdt[it][m] = qnan;
      utTofSeg[it][m]  = qnan;
      dtTofSeg[it][m]  = qnan;
      csTof[MaxDepth*it + m]  = 0;
      TofSeg[MaxDepth*it + m] = qnan;
      tTof[MaxDepth*it + m]   = qnan;
      dtTof[MaxDepth*it + m]  = qnan;
      deTof[MaxDepth*it + m]  = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegLAC; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      csLac[MaxDepth*it + m]  = 0;
      LacSeg[MaxDepth*it + m] = qnan;
      tLac[MaxDepth*it + m]   = qnan;
    }
  }
}


//_____________________________________________________________________________
namespace root
{
Event  event;
Dst    dst;
TH1*   h[MaxHist];
TTree* tree;
TTree* hodo;
enum eDetHid {
  BH1Hid    = 10000,
  BH2Hid    = 20000,
  BACHid    = 30000,
  HTOFHid   = 40000,
  BVHHid    = 50000,
  TOFHid    = 60000,
  LACHid    = 70000,
  WCHid     = 80000,
  WCSUMHid  = 90000
};
}

//_____________________________________________________________________________
Bool_t
UserHodoscope::ProcessingBegin()
{
  event.clear();
  dst.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserHodoscope::ProcessingNormal()
{
  static const auto MinTdcBH1 = gUser.GetParameter("TdcBH1", 0);
  static const auto MaxTdcBH1 = gUser.GetParameter("TdcBH1", 1);
  static const auto MinTdcBH2 = gUser.GetParameter("TdcBH2", 0);
  static const auto MaxTdcBH2 = gUser.GetParameter("TdcBH2", 1);
  static const auto MinTdcBAC = gUser.GetParameter("TdcBAC", 0);
  static const auto MaxTdcBAC = gUser.GetParameter("TdcBAC", 1);
  static const auto MinTdcHTOF = gUser.GetParameter("TdcHTOF", 0);
  static const auto MaxTdcHTOF = gUser.GetParameter("TdcHTOF", 1);
  static const auto MinTdcBVH = gUser.GetParameter("TdcBVH", 0);
  static const auto MaxTdcBVH = gUser.GetParameter("TdcBVH", 1);
  static const auto MinTdcTOF = gUser.GetParameter("TdcTOF", 0);
  static const auto MaxTdcTOF = gUser.GetParameter("TdcTOF", 1);
  static const auto MinTdcLAC = gUser.GetParameter("TdcLAC", 0);
  static const auto MaxTdcLAC = gUser.GetParameter("TdcLAC", 1);
  static const auto MinTdcWC = gUser.GetParameter("TdcWC", 0);
  static const auto MaxTdcWC = gUser.GetParameter("TdcWC", 1);
  static const auto MinTdcSCH  = gUser.GetParameter("TdcSCH",  0);
  static const auto MaxTdcSCH  = gUser.GetParameter("TdcSCH",  1);
  static const auto MinTimeSCH = gUser.GetParameter("TimeSCH", 0);
  static const auto MaxTimeSCH = gUser.GetParameter("TimeSCH", 1);
#if HodoHitPos
  static const auto PropVelBH2 = gUser.GetParameter("PropagationBH2");
  static const auto PropVelHTOF = gUser.GetParameter("PropagationHTOF");
#endif

  rawData->DecodeHits();

  gRM.Decode();

  // event.evnum = gRM.EventNumber();
  event.evnum = gUnpacker.get_event_number();
  event.spill = gRM.SpillNumber();
  dst.evnum   = gRM.EventNumber();
  dst.spill   = gRM.SpillNumber();

  HF1(1, 0);

  //**************************************************************************
  //****************** RawData

  // Trigger Flag
  std::bitset<NumOfSegTrig> trigger_flag;
  for(const auto& hit: rawData->GetTrigRawHC()){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc1();
    if(tdc > 0){
      event.trigpat[trigger_flag.count()] = seg;
      event.trigflag[seg] = tdc;
      dst.trigpat[trigger_flag.count()] = seg;
      dst.trigflag[seg] = tdc;
      trigger_flag.set(seg);
      HF1(10, seg);
      HF1(10+seg, tdc);
    }
  }

  if(trigger_flag[trigger::kSpillEnd]) return true;

  HF1(1, 1);

  ///// BH1
  {
    Int_t bh1_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetBH1RawHC();
    Int_t nh = cont.size();
    HF1(BH1Hid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(BH1Hid +1, seg-0.5);
      // Up
      Int_t Au = hit->GetAdcUp();
      HF1(BH1Hid +100*seg +1, Au);
      event.bh1ua[seg-1] = Au;
      Bool_t is_hit_u = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        HF1(BH1Hid +100*seg +3, T);
        event.bh1ut[seg-1][m] = T;
        if(MinTdcBH1 < T && T < MaxTdcBH1) is_hit_u = true;
      }
      if(is_hit_u) HF1(BH1Hid +100*seg +5, Au);
      else         HF1(BH1Hid +100*seg +7, Au);
      // Down
      Int_t Ad = hit->GetAdcDown();
      HF1(BH1Hid +100*seg +2, Ad);
      event.bh1da[seg-1] = Ad;
      Bool_t is_hit_d = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcDown(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcDown(m);
        HF1(BH1Hid +100*seg +4, T);
        event.bh1dt[seg-1][m] = T;
        if(MinTdcBH1 < T && T < MaxTdcBH1) is_hit_d = true;
      }
      if(is_hit_d) HF1(BH1Hid +100*seg +6, Ad);
      else         HF1(BH1Hid +100*seg +8, Ad);
      // HitPat
      if(is_hit_u || is_hit_d){
        ++nh1; HF1(BH1Hid +3, seg-0.5);
      }
      if(is_hit_u && is_hit_d){
        event.bh1hitpat[bh1_nhits++] = seg;
        ++nh2; HF1(BH1Hid +5, seg-0.5);
      }
    }
    HF1(BH1Hid +2, nh1); HF1(BH1Hid +4, nh2);
    event.bh1nhits = bh1_nhits;
  }

  ///// BH2
  {
    Int_t bh2_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetBH2RawHC();
    Int_t nh = cont.size();
    HF1(BH2Hid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(BH2Hid +1, seg-0.5);
      // Up
      Int_t Au = hit->GetAdcUp();
      HF1(BH2Hid +100*seg +1, Au);
      event.bh2ua[seg-1] = Au;
      Bool_t is_hit_u = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        HF1(BH2Hid +100*seg +3, T);
        event.bh2ut[seg-1][m] = T;
        if(MinTdcBH2 < T && T < MaxTdcBH2) is_hit_u = true;
      }
      if(is_hit_u) HF1(BH2Hid +100*seg +5, Au);
      else         HF1(BH2Hid +100*seg +7, Au);
      // Down
      Int_t Ad = hit->GetAdcDown();
      HF1(BH2Hid +100*seg +2, Double_t(Ad));
      event.bh2da[seg-1] = Ad;
      Bool_t is_hit_d = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcDown(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcDown(m);
        HF1(BH2Hid +100*seg +4, Double_t(T));
        event.bh2dt[seg-1][m] = T;
        if(MinTdcBH2 < T && T < MaxTdcBH2) is_hit_d = true;
      }
      if(is_hit_d) HF1(BH2Hid +100*seg +6, Ad);
      else         HF1(BH2Hid +100*seg +8, Ad);
      // HitPat
      if(is_hit_u || is_hit_d){
        ++nh1; HF1(BH2Hid +3, seg-0.5);
      }
      if(is_hit_u && is_hit_d){
        event.bh2hitpat[bh2_nhits++] = seg;
        ++nh2; HF1(BH2Hid +5, seg-0.5);
      }
    }
    HF1(BH2Hid +2, nh1); HF1(BH2Hid +4, nh2);
    event.bh2nhits = bh2_nhits;
  }

  ///// BAC
  {
    Int_t bac_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetBACRawHC();
    Int_t nh = cont.size();
    HF1(BACHid, nh);
    Int_t nh1 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(BACHid+1, seg-0.5);
      Int_t A = hit->GetAdcUp();
      HF1(BACHid+100*seg+1, A);
      event.baca[seg-1] = A;
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        HF1(BACHid+100*seg+3, T);
        event.bact[seg-1][m] = T;
        if(MinTdcBAC < T && T < MaxTdcBAC) is_hit = true;
      }
      if(is_hit) HF1(BACHid+100*seg+5, A);
      else       HF1(BACHid+100*seg+7, A);
      // Hitpat
      if(is_hit){
        event.bachitpat[bac_nhits++] = seg;
        ++nh1; HF1(BACHid+3, seg-0.5);
      }
    }
    HF1(BACHid+2, nh1);
    event.bacnhits = bac_nhits;
  }

  ///// HTOF
  {
    Int_t htof_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetHTOFRawHC();
    Int_t nh = cont.size();
    HF1(HTOFHid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit* hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(HTOFHid+1, seg-0.5);
      // Up
      Int_t Au = hit->GetAdcUp();
      HF1(HTOFHid+100*seg+1, Au);
      event.htofua[seg-1] = Au;
      Bool_t is_hit_u = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        HF1(HTOFHid +100*seg +3, T);
        event.htofut[seg-1][m] = T;
        if(MinTdcHTOF < T && T < MaxTdcHTOF) is_hit_u = true;
      }
      if(is_hit_u) HF1(HTOFHid+100*seg+5, Au);
      else         HF1(HTOFHid+100*seg+7, Au);
      // Down
      Int_t Ad = hit->GetAdcDown();
      HF1(HTOFHid+100*seg+2, Ad);
      event.htofda[seg-1] = Ad;
      Bool_t is_hit_d = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcDown(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcDown(m);
        HF1(HTOFHid +100*seg +4, T);
        event.htofdt[seg-1][m] = T;
        if(MinTdcHTOF < T && T < MaxTdcHTOF) is_hit_d = true;
      }
      if(is_hit_d) HF1(HTOFHid+100*seg+6, Ad);
      else         HF1(HTOFHid+100*seg+8, Ad);
      // HitPat
      if(is_hit_u || is_hit_d){
        ++nh1; HF1(HTOFHid+3, seg-0.5);
      }
      if(is_hit_u && is_hit_d){
        event.htofhitpat[htof_nhits++] = seg;
        ++nh2; HF1(HTOFHid+5, seg-0.5);
      }
    }
    HF1(HTOFHid+2, nh1); HF1(HTOFHid+4, nh2);
    event.htofnhits = htof_nhits;
  }

  ///// BVH
  {
    Int_t bvh_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetBVHRawHC();
    Int_t nh = cont.size();
    HF1(BVHHid, nh);
    Int_t nh1 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(BVHHid+1, seg-0.5);
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        HF1(BVHHid+100*seg+3, T);
        event.bvht[seg-1][m] = T;
        if(MinTdcBVH < T && T < MaxTdcBVH) is_hit = true;
      }
      // Hitpat
      if(is_hit){
        event.bvhhitpat[bvh_nhits++] = seg;
        event.tofhitpat_3dmtx[bvh_nhits] = seg+NumOfSegTOF;
        ++nh1; HF1(BVHHid+3, seg-0.5);
      }
    }
    HF1(BVHHid+2, nh1);
    event.bvhnhits = bvh_nhits;
  }

  ///// TOF
  {
    Int_t tof_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetTOFRawHC();
    Int_t nh = cont.size();
    HF1(TOFHid, Double_t(nh));
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(TOFHid+1, seg-0.5);
      // Up
      Int_t Au = hit->GetAdcUp();
      HF1(TOFHid+100*seg+1, Au);
      event.tofua[seg-1] = Au;
      dst.tofua[seg-1] = Au;
      Bool_t is_hit_u = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        HF1(TOFHid +100*seg +3, T);
        event.tofut[seg-1][m] = T;
        dst.tofut[seg-1][m] = T;
        if(MinTdcTOF < T && T < MaxTdcTOF) is_hit_u = true;
      }
      if(is_hit_u) HF1(TOFHid+100*seg+5, Au);
      else         HF1(TOFHid+100*seg+7, Au);
      // Down
      Int_t Ad = hit->GetAdcDown();
      HF1(TOFHid+100*seg+2, Ad);
      event.tofda[seg-1] = Ad;
      dst.tofda[seg-1] = Ad;
      Bool_t is_hit_d = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcDown(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcDown(m);
        HF1(TOFHid +100*seg +4, Double_t(T));
        event.tofdt[seg-1][m] = T;
        dst.tofdt[seg-1][m] = T;
        if(MinTdcTOF < T && T < MaxTdcTOF)  is_hit_d = true;
      }
      if(is_hit_d) HF1(TOFHid+100*seg+6, Ad);
      else         HF1(TOFHid+100*seg+8, Ad);
      // Hitpat
      if(is_hit_u || is_hit_d){
        ++nh1; HF1(TOFHid+3, seg-0.5);
      }
      if(is_hit_u && is_hit_d){
        event.tofhitpat[tof_nhits++] = seg;
        event.tofhitpat_3dmtx[event.bvhnhits+tof_nhits] = seg;
        ++nh2; HF1(TOFHid+5, seg-0.5);
      }
    }
    HF1(TOFHid+2, nh1); HF1(TOFHid+4, nh2);
    event.tofnhits = tof_nhits;
    event.tofnhits_3dmtx = event.bvhnhits+tof_nhits;
  }

  ///// LAC
  {
    Int_t lac_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetLACRawHC();
    Int_t nh = cont.size();
    HF1(LACHid, nh);
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(LACHid +1, seg-0.5);
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        HF1(LACHid +100*seg +3, Double_t(T));
        event.lact[seg-1][m] = T;
        if(MinTdcLAC < T && T < MaxTdcLAC) is_hit = true;
      }
      //HitPat
      if(is_hit){
        HF1(LACHid +3, seg-0.5);
        HF1(LACHid +5, seg-0.5);
        event.lachitpat[lac_nhits++] = seg;
      }
    }
    HF1(LACHid +2, lac_nhits);
    HF1(LACHid +4, lac_nhits);
    event.lacnhits = lac_nhits;
  }

  ///// WC
  {
    Int_t wc_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetWCRawHC();
    Int_t nh = cont.size();
    HF1(WCHid, Double_t(nh));
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(WCHid +1, seg-0.5);
      // Up
      Int_t Au = hit->GetAdcUp();
      HF1(WCHid +100*seg +1, Au);
      event.wcua[seg-1] = Au;
      Bool_t is_hit_u = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        HF1(WCHid +100*seg +3, T);
        event.wcut[seg-1][m] = T;
        if(MinTdcWC < T && T < MaxTdcWC) is_hit_u = true;
      }
      if(is_hit_u) HF1(WCHid +100*seg +5, Au);
      else         HF1(WCHid +100*seg +7, Au);
      // Down
      Int_t Ad = hit->GetAdcDown();
      HF1(WCHid +100*seg +2, Ad);
      event.wcda[seg-1] = Ad;
      Bool_t is_hit_d = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcDown(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcDown(m);
        HF1(WCHid +100*seg +4, Double_t(T));
        event.wcdt[seg-1][m] = T;
        if(MinTdcWC < T && T < MaxTdcWC) is_hit_d = true;
      }
      if(is_hit_d) HF1(WCHid +100*seg +6, Ad);
      else         HF1(WCHid +100*seg +8, Ad);
      // HitPat
      if(is_hit_u || is_hit_d){
        ++nh1; HF1(WCHid +3, seg-0.5);
      }
      if(is_hit_u && is_hit_d){
        event.wchitpat[wc_nhits] = seg;
        wc_nhits++;
        ++nh2; HF1(WCHid +5, seg-0.5);
      }
    }
    HF1(WCHid +2, Double_t(nh1)); HF1(WCHid +4, Double_t(nh2));
    event.wcnhits = wc_nhits;
  }

  ///// WCSUM
  {
    Int_t wcsum_nhits = 0;
    const auto &cont = rawData->GetWCSUMRawHC();
    Int_t nh = cont.size();
    for( Int_t i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Int_t A = hit->GetAdcUp();
      HF1( WCSUMHid +100*seg +1, A);
      event.wcsuma[seg-1] = A;
      Bool_t is_hit = false;
      for(Int_t m=0, mh=hit->GetSizeTdcUp(); m<mh; ++m){
        Int_t T = hit->GetTdcUp(m);
        HF1(WCSUMHid +100*seg +3, T);
        event.wcsumt[seg-1][m] = T;
        if(MinTdcWC < T && T < MaxTdcWC) is_hit = true;
      }
      if(is_hit) HF1(WCSUMHid +100*seg +5, A);
      else       HF1(WCSUMHid +100*seg +7, A);
      // HitPat
      if(is_hit){
        event.wcsumhitpat[wcsum_nhits++] = seg;
        HF1(WCSUMHid +1, seg-0.5);
      }
    }
    HF1(WCSUMHid, wcsum_nhits);
    event.wcsumnhits = wcsum_nhits;
  }


  //**************************************************************************
  //****************** NormalizedData

  //BH1
  hodoAna->DecodeBH1Hits(rawData);
  //  hodoAna->TimeCutBH1(-10, 10);
  hodoAna->TimeCutBH1(-2, 2);
  {
    Int_t nh = hodoAna->GetNHitsBH1();
    HF1(BH1Hid+10, Double_t(nh));
    Int_t nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      Hodo2Hit *hit = hodoAna->GetHitBH1(i);
      if(!hit) continue;
      Int_t seg = hit->SegmentId()+1;

      Int_t n_mhit = hit->GetNumOfHit();
      for(Int_t m=0; m<n_mhit; ++m){
        HF1(BH1Hid+11, seg-0.5);

        Double_t au  = hit->GetAUp(),    ad = hit->GetADown();
        Double_t tu  = hit->GetTUp(m),   td = hit->GetTDown(m);
        Double_t ctu = hit->GetCTUp(m), ctd = hit->GetCTDown(m);
        Double_t mt  = hit->MeanTime(m),cmt = hit->CMeanTime(m);
        Double_t de  = hit->DeltaE();
        Double_t ude  = hit->UDeltaE();
        Double_t dde  = hit->DDeltaE();
        event.bh1mt[seg-1][m] = mt;
        event.bh1utime[seg-1][m] = tu;
        event.bh1dtime[seg-1][m] = td;
        event.bh1uctime[seg-1][m] = ctu;
        event.bh1dctime[seg-1][m] = ctd;
        event.bh1de[seg-1]    = de;
        event.bh1ude[seg-1]    = ude;
        event.bh1dde[seg-1]    = dde;
        HF1(BH1Hid+100*seg+11, tu);      HF1(BH1Hid+100*seg+12, td);
        HF1(BH1Hid+100*seg+13, mt);
        HF1(BH1Hid+100*seg+17, ctu);     HF1(BH1Hid+100*seg+18, ctd);
        HF1(BH1Hid+100*seg+19, cmt);     HF1(BH1Hid+100*seg+20, ctu-ctd);
        HF2(BH1Hid+100*seg+21, tu, au);  HF2(BH1Hid+100*seg+22, td, ad);
        HF2(BH1Hid+100*seg+23, ctu, au); HF2(BH1Hid+100*seg+24, ctd, ad);
        HF2(BH1Hid+100*seg+31, au, tu);  HF2(BH1Hid+100*seg+32, ad, td);
        HF2(BH1Hid+100*seg+33, au, ctu);  HF2(BH1Hid+100*seg+34, ad, ctd);
        HF1(BH1Hid+12, cmt);

        if(m == 0){
          HF1(BH1Hid+100*seg+14, au);    HF1(BH1Hid+100*seg+15, ad);
          HF1(BH1Hid+100*seg+16, de);    HF1(BH1Hid+13, de);
        }

        if(de>0.5){
          ++nh2; HF1(BH1Hid+15, seg-0.5);
          HF1(BH1Hid+16, cmt);
        }
      }
    }

    HF1(BH1Hid+14, Double_t(nh2));
    for(Int_t i1=0; i1<nh; ++i1){
      Hodo2Hit *hit1 = hodoAna->GetHitBH1(i1);
      if(!hit1 || hit1->DeltaE()<=0.5) continue;
      Int_t seg1 = hit1->SegmentId()+1;
      for(Int_t i2=0; i2<nh; ++i2){
        if(i1==i2) continue;
        Hodo2Hit *hit2=hodoAna->GetHitBH1(i2);
        if(!hit2 || hit2->DeltaE()<=0.5) continue;
        Int_t seg2 = hit2->SegmentId()+1;

        if(1 == hit1->GetNumOfHit() && 1 == hit2->GetNumOfHit()){
          Double_t ct1 = hit1->CMeanTime(), ct2 = hit2->CMeanTime();
          HF2(BH1Hid+21, seg1-0.5, seg2-0.5);
          HF2(BH1Hid+22, ct1, ct2);
          HF1(BH1Hid+23, ct2-ct1);
          if(std::abs(ct2-ct1)<2.0){
            HF2(BH1Hid+24, seg1-0.5, seg2-0.5);
          }
        }
      }//for(i2)
    }//for(i1)
  }

  // BH2
  hodoAna->DecodeBH2Hits(rawData);
  // hodoAna->TimeCutBH2(-2, 2);
  {
    Int_t nh = hodoAna->GetNHitsBH2();
    HF1(BH2Hid+10, Double_t(nh));
    Int_t nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      BH2Hit *hit = hodoAna->GetHitBH2(i);
      if(!hit) continue;
      Int_t seg = hit->SegmentId()+1;

      Int_t n_mhit = hit->GetNumOfHit();
      for(Int_t m=0; m<n_mhit; ++m){
        HF1(BH2Hid+11, seg-0.5);
        Double_t au  = hit->GetAUp(),   ad  = hit->GetADown();
        Double_t tu  = hit->GetTUp(m),  td  = hit->GetTDown(m);
        Double_t ctu = hit->GetCTUp(m), ctd = hit->GetCTDown(m);
        Double_t mt  = hit->MeanTime(m),cmt = hit->CMeanTime(m);
        Double_t de  = hit->DeltaE();
        Double_t ude  = hit->UDeltaE();
        Double_t dde  = hit->DDeltaE();
        Double_t ut0  = hit->UTime0(m), dt0  = hit->DTime0(m);
        Double_t uct0 = hit->UCTime0(m),dct0 = hit->DCTime0(m);
        Double_t t0   = hit->Time0(m),  ct0  = hit->CTime0(m);
#if HodoHitPos
	Double_t ctdiff = hit->TimeDiff(m);
	event.bh2hitpos[seg-1][m] = 0.5*PropVelBH2*ctdiff;
#endif
        event.bh2mt[seg-1][m] = mt;
        event.bh2utime[seg-1][m] = tu;
        event.bh2dtime[seg-1][m] = td;
        event.bh2uctime[seg-1][m] = ctu;
        event.bh2dctime[seg-1][m] = ctd;
        event.bh2de[seg-1]    = de;
        event.bh2ude[seg-1]    = ude;
        event.bh2dde[seg-1]    = dde;
        event.t0[seg-1][m]    = t0;
        event.ct0[seg-1][m]   = ct0;

	HF1(BH2Hid+100*seg+11, tu);      HF1(BH2Hid+100*seg+12, td);
        HF1(BH2Hid+100*seg+13, mt);
        HF1(BH2Hid+100*seg+17, ctu);     HF1(BH2Hid+100*seg+18, ctd);
        HF1(BH2Hid+100*seg+19, cmt);     HF1(BH2Hid+100*seg+20, ctu-ctd);
        HF1(BH2Hid+100*seg+21, ut0);     HF1(BH2Hid+100*seg+22, dt0);
        HF1(BH2Hid+100*seg+23, uct0);    HF1(BH2Hid+100*seg+24, dct0);
        HF1(BH2Hid+100*seg+25, t0);      HF1(BH2Hid+100*seg+26, ct0);
        HF2(BH2Hid+100*seg+27, tu, au);  HF2(BH2Hid+100*seg+28, td, ad);
        HF2(BH2Hid+100*seg+29, ctu, au); HF2(BH2Hid+100*seg+30, ctd, ad);

        // HF1(BH2Hid+100*seg+11, tu); HF1(BH2Hid+100*seg+13, mt);
        // HF1(BH2Hid+100*seg+14, au); HF1(BH2Hid+100*seg+16, de);
        // HF1(BH2Hid+100*seg+17, ctu); HF1(BH2Hid+100*seg+19, cmt);
        // HF2(BH2Hid+100*seg+21, tu, au); HF2(BH2Hid+100*seg+23, ctu, au);
        // HF1(BH2Hid+12, cmt);

        if(m == 0){
          HF1(BH2Hid+100*seg+14, au);	   HF1(BH2Hid+100*seg+15, ad);
          HF1(BH2Hid+100*seg+16, de);    HF1(BH2Hid+13, de);
        }

        if(de>0.5){
          ++nh2; HF1(BH2Hid+15, seg-0.5); HF1(BH2Hid+16, cmt);
        }
      }
    }//for(i)
    HF1(BH2Hid+14, Double_t(nh2));
    for(Int_t i1=0; i1<nh; ++i1){
      BH2Hit *hit1 = hodoAna->GetHitBH2(i1);
      if(!hit1 || hit1->DeltaE()<=0.5) continue;
      Int_t seg1 = hit1->SegmentId()+1;
      for(Int_t i2=0; i2<nh; ++i2){
        if(i1==i2) continue;
        BH2Hit *hit2 = hodoAna->GetHitBH2(i2);
        if(!hit2 || hit2->DeltaE()<=0.5) continue;
        Int_t seg2 = hit2->SegmentId()+1;

        if(1 == hit1->GetNumOfHit() && 1 == hit2->GetNumOfHit()){
          Double_t ct1 = hit1->CMeanTime(), ct2 = hit2->CMeanTime();
          HF2(BH2Hid+21, seg1-0.5, seg2-0.5);
          HF2(BH2Hid+22, ct1, ct2);
          HF1(BH2Hid+23, ct2-ct1);
          if(std::abs(ct2-ct1)<2.0){
            HF2(BH2Hid+24, seg1-0.5, seg2-0.5);
          }
        }
      }//for(i2)
    }//for(i1)

    Int_t nc=hodoAna->GetNClustersBH2();
    HF1(BH2Hid+30, Double_t(nc));
    Int_t nc2=0;

    for(Int_t i=0; i<nc; ++i){
      BH2Cluster *cluster=hodoAna->GetClusterBH2(i);
      if(!cluster) continue;
      Int_t cs=cluster->ClusterSize();
      Double_t ms = cluster->MeanSeg()+1;
      Double_t cmt= cluster->CMeanTime();
      Double_t de = cluster->DeltaE();
      // Double_t mt = cluster->MeanTime();
      HF1(BH2Hid+31, Double_t(cs));
      HF1(BH2Hid+32, ms-0.5);
      HF1(BH2Hid+33, cmt); HF1(BH2Hid+34, de);
      if(de>0.5){
        ++nc2; HF1(BH2Hid+36, cmt);
      }

      for(Int_t i2=0; i2<nc; ++i2){
        if(i2==i) continue;
        BH2Cluster *cl2=hodoAna->GetClusterBH2(i2);
        if(!cl2) continue;
        Double_t ms2=cl2->MeanSeg()+1, cmt2=cl2->CMeanTime(),
          de2=cl2->DeltaE();
        if(de<=0.5 || de2<=0.5) continue;
        HF2(BH2Hid+41, ms-0.5, ms2-0.5);
        HF2(BH2Hid+42, cmt, cmt2);
        HF1(BH2Hid+43, cmt2-cmt);
        if(std::abs(cmt2-cmt)<2.0){
          HF2(BH2Hid+44, ms-0.5, ms2-0.5);
        }
      }//for(i2)
    }//for(i)
    HF1(BH2Hid+35, Double_t(nc2));

    BH2Cluster* cl_time0 = hodoAna->GetTime0BH2Cluster();
    if(cl_time0){
      event.Time0Seg = cl_time0->MeanSeg()+1;
      event.deTime0 = cl_time0->DeltaE();
      event.Time0 = cl_time0->Time0();
      event.CTime0 = cl_time0->CTime0();
      dst.Time0Seg = cl_time0->MeanSeg()+1;
      dst.deTime0 = cl_time0->DeltaE();
      dst.Time0 = cl_time0->Time0();
      dst.CTime0 = cl_time0->CTime0();
      HF1(100, cl_time0->MeanTime());
      HF1(101, cl_time0->CTime0());
    }

    // BTOF0 segment
    HodoCluster* cl_btof0 = (dst.Time0Seg > 0)
      ? hodoAna->GetBtof0BH1Cluster(dst.CTime0) : nullptr;
    if(cl_btof0){
      event.Btof0Seg = cl_btof0->MeanSeg()+1;
      event.deBtof0 = cl_btof0->DeltaE();
      event.Btof0 = cl_btof0->MeanTime() - dst.Time0;
      event.CBtof0 = cl_btof0->CMeanTime() - dst.CTime0;
      dst.Btof0Seg = cl_btof0->MeanSeg()+1;
      dst.deBtof0 = cl_btof0->DeltaE();
      dst.Btof0 = cl_btof0->MeanTime() - dst.Time0;
      dst.CBtof0 = cl_btof0->CMeanTime() - dst.CTime0;
      HF1(102, dst.CBtof0);
      HF2(103, dst.deTime0, dst.CBtof0);
      HF2(104, dst.deBtof0, dst.CBtof0);
    }
  }

  // BH1 with BH2 gate
  {
    Int_t nhbh2 = hodoAna->GetNHitsBH2();
    if(nhbh2){
      Int_t    seg2 = hodoAna->GetHitBH2(0)->SegmentId()+1;
      Int_t n_mhit2 = hodoAna->GetHitBH2(0)->GetNumOfHit();
      for(Int_t m2 = 0; m2<n_mhit2; ++m2){
        Double_t mt2  = hodoAna->GetHitBH2(0)->CTime0(m2);
        Int_t    nh   = hodoAna->GetNHitsBH1();
        for(Int_t i=0; i<nh; ++i){
          Hodo2Hit *hit = hodoAna->GetHitBH1(i);
          if(!hit) continue;

          Int_t n_mhit1 = hit->GetNumOfHit();
          for(Int_t m1 = 0; m1<n_mhit1; ++m1){
            Int_t seg1 = hit->SegmentId()+1;
            Double_t tu1 = hit->GetTUp(m1), td1 = hit->GetTDown(m1);
            Double_t mt1 = hit->MeanTime(m1);
            HF1(BH1Hid+100*seg1+1100+21+seg2*10, tu1);
            HF1(BH1Hid+100*seg1+1100+22+seg2*10, td1);
            HF1(BH1Hid+100*seg1+1100+23+seg2*10, mt1);

            //For BH1vsBH2 Correlation
            HF2(BH1Hid+100*seg1+2200+21+seg2*10, tu1, mt2);
            HF2(BH1Hid+100*seg1+2200+22+seg2*10, td1, mt2);
            HF2(BH1Hid+100*seg1+2200+23+seg2*10, mt1, mt2);
          }// for(m1)
        }// for(bh1:seg)
      }// for(m2)
    }
    for(Int_t i2=0; i2<nhbh2; ++i2){
      Int_t    seg2 = hodoAna->GetHitBH2(i2)->SegmentId()+1;
      Int_t n_mhit2 = hodoAna->GetHitBH2(i2)->GetNumOfHit();
      for(Int_t m2 = 0; m2<n_mhit2; ++m2){
        Double_t ct0  = hodoAna->GetHitBH2(i2)->CTime0();
        Double_t t0   = hodoAna->GetHitBH2(i2)->Time0();
        Int_t nhbh1=hodoAna->GetNHitsBH1();
        for(Int_t i1=0; i1<nhbh1; ++i1){
          Hodo2Hit *hit=hodoAna->GetHitBH1(i1);
          if(!hit) continue;

          Int_t n_mhit1 = hit->GetNumOfHit();
          for(Int_t m1 = 0; m1<n_mhit1; ++m1){
            Int_t seg1=hit->SegmentId()+1;
            Double_t mt1 = hit->MeanTime(m1);
            Double_t cmt1 = hit->CMeanTime(m1);
            Double_t btof = mt1 - t0;
            Double_t cbtof = cmt1 - ct0;
            event.btof[seg1-1][seg2-1] = btof;
            event.cbtof[seg1-1][seg2-1] = cbtof;
            HF1(BH1Hid+100*seg1+1100+24+seg2*10, mt1-ct0);
            HF1(BH1Hid+100*seg1+2200+104, mt1-ct0);
            if(seg1==5){
              HF1(30+seg2, mt1-ct0);
            }
          }// for(m1)
        }// for(bh1:seg)
      }// for(m2)
    }// for(bh2:seg)
  }

  // BH1-BH2
  {
    Int_t nhbh1 = hodoAna->GetNHitsBH1();
    Int_t nhbh2 = hodoAna->GetNHitsBH2();
    for(Int_t i2=0; i2<nhbh2; ++i2){
      BH2Hit *hit2 = hodoAna->GetHitBH2(i2);
      if(!hit2) continue;

      Int_t n_mhit2 = hit2->GetNumOfHit();
      Int_t    seg2 = hit2->SegmentId()+1;
      Double_t de2  = hit2->DeltaE();
      for(Int_t m2 = 0; m2<n_mhit2; ++m2){
        Double_t t0      = hit2->Time0(m2);

        for(Int_t i1=0; i1<nhbh1; ++i1){
          Hodo2Hit *hit1 = hodoAna->GetHitBH1(i1);
          if(!hit1) continue;

          Int_t n_mhit1  = hit1->GetNumOfHit();
          Int_t    seg1  = hit1->SegmentId()+1;
          Double_t de1   = hit1->DeltaE();

          for(Int_t m1 = 0; m1<n_mhit1; ++m1){
            Double_t mt1   = hit1->MeanTime(m1);
            HF1(201, mt1-t0);
            HF2(202, seg1-0.5, seg2-0.5);
            //For BH1vsBH2 Correlation
            HF2(203, t0, mt1);
            HF1(204, mt1);
            HF1(205, t0);

            if(de1>0.5 && de2>0.5){
              HF1(211, mt1-t0);
              HF2(212, seg1-0.5, seg2-0.5);
            }
          }// for(m1)
        }// for(bh1:seg)
      }// for(m2)
    }// for(bh2:seg)
  }

#if 1
  // BH1-BH2 PHC
  {
    Int_t nh1 = hodoAna->GetNHitsBH1();
    Int_t nh2 = hodoAna->GetNHitsBH2();
    for(Int_t i2=0; i2<nh2; ++i2){
      BH2Hit *hit2 = hodoAna->GetHitBH2(i2);
      Int_t     seg2 = hit2->SegmentId()+1;
      Double_t  au2  = hit2->GetAUp(),  ad2  = hit2->GetADown();
      Int_t n_mhit2  = hit2->GetNumOfHit();
      for(Int_t m2 = 0; m2<n_mhit2; ++m2){
        Double_t  tu2  = hit2->GetTUp(m2),  td2  = hit2->GetTDown(m2);
        Double_t  ctu2 = hit2->GetCTUp(m2), ctd2 = hit2->GetCTDown(m2);
        // Double_t  t0   = hit2->Time0();
        Double_t  ct0  = hit2->CTime0(m2);
        Double_t  tofs = ct0-(ctu2+ctd2)/2.;
        for(Int_t i1=0; i1<nh1; ++i1){
          Hodo2Hit *hit1 = hodoAna->GetHitBH1(i1);
          Int_t       seg1 = hit1->SegmentId()+1;
          Double_t    au1  = hit1->GetAUp(),  ad1  = hit1->GetADown();
          Int_t    n_mhit1 = hit1->GetNumOfHit();
          for(Int_t m1 = 0; m1<n_mhit1; ++m1){
            Double_t    tu1  = hit1->GetTUp(m1),  td1  = hit1->GetTDown(m1);
            Double_t    ctu1 = hit1->GetCTUp(m1), ctd1 = hit1->GetCTDown(m1);
            Double_t    cmt1 = hit1->CMeanTime(m1);
            if(trigger_flag[trigger::kBeamA] == 0) continue;
            if(event.bacnhits > 0) continue;
            HF2(100*seg1+BH1Hid+81, au1, ct0-0.5*(ctu1+ctd1));
            HF2(100*seg1+BH1Hid+82, ad1, ct0-0.5*(ctu1+ctd1));
            HF2(100*seg1+BH1Hid+83, au1, ct0-tu1);
            HF2(100*seg1+BH1Hid+84, ad1, ct0-td1);
            HF2(100*seg2+BH2Hid+81, au2, (cmt1-tofs)-0.5*(ctu2+ctd2));
            HF2(100*seg2+BH2Hid+82, ad2, (cmt1-tofs)-0.5*(ctu2+ctd2));
            HF2(100*seg2+BH2Hid+83, au2, (cmt1-tofs)-tu2);
            HF2(100*seg2+BH2Hid+84, ad2, (cmt1-tofs)-td2);
          }// for(m1)
        }// for(bh1:seg)
      }// for(m2)
    }// for(bh2:seg)
  }
#endif

  // BAC
  {
    hodoAna->DecodeBACHits(rawData);
    Int_t nh=hodoAna->GetNHitsBAC();
    dst.nhBac = nh;
    HF1(BACHid+10, Double_t(nh));
    for(Int_t i=0; i<nh; ++i){
      Hodo1Hit *hit=hodoAna->GetHitBAC(i);
      if(!hit) continue;
      Int_t seg=hit->SegmentId()+1;
      HF1(BACHid+11, seg-0.5);
      Double_t a=hit->GetA();
      event.bacde[seg-1]  = a;
      Int_t n_mhit = hit->GetNumOfHit();
      for(Int_t m=0; m<n_mhit; ++m){
        Double_t t=hit->GetT(m), ct=hit->GetCT(m);
        event.bacmt[i][m] = ct;
        HF1(BACHid+100*seg+11, t);
        HF1(BACHid+100*seg+12, a);
        HF1(BACHid+100*seg+13, ct);
      }
    }
  }

  // HTOF
  hodoAna->DecodeHTOFHits(rawData);
  {
    Int_t nh = hodoAna->GetNHitsHTOF();
    HF1(HTOFHid+10, Double_t(nh));
    Int_t nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      Hodo2Hit *hit = hodoAna->GetHitHTOF(i);
      if(!hit) continue;
      Int_t seg = hit->SegmentId()+1;
      Int_t n_mhit = hit->GetNumOfHit();
      for(Int_t m=0; m<n_mhit; ++m){
        HF1(HTOFHid+11, seg-0.5);
        Double_t au = hit->GetAUp(), ad = hit->GetADown();
        Double_t tu = hit->GetTUp(), td = hit->GetTDown();
        Double_t ctu = hit->GetCTUp(), ctd = hit->GetCTDown();
        Double_t mt = hit->MeanTime(), cmt = hit->CMeanTime();
        Double_t de = hit->DeltaE();
        Double_t ude = hit->UDeltaE();
        Double_t dde = hit->DDeltaE();
#if HodoHitPos
	Double_t ctdiff = hit->TimeDiff(m);
	event.htofhitpos[seg-1][m] = 0.5*PropVelHTOF*ctdiff;
#endif
	event.htofmt[seg-1][m] = mt;
        event.htofutime[seg-1][m] = tu;
        event.htofdtime[seg-1][m] = td;
        event.htofuctime[seg-1][m] = ctu;
        event.htofdctime[seg-1][m] = ctd;
        event.htofde[seg-1] = de;
        event.htofude[seg-1] = ude;
        event.htofdde[seg-1] = dde;
        HF1(HTOFHid+100*seg+11, tu); HF1(HTOFHid+100*seg+12, td);
        HF1(HTOFHid+100*seg+13, mt);
        HF1(HTOFHid+100*seg+17, ctu); HF1(HTOFHid+100*seg+18, ctd);
        HF1(HTOFHid+100*seg+19, cmt); HF1(HTOFHid+100*seg+20, ctu-ctd);
        HF2(HTOFHid+100*seg+31, au, tu);  HF2(HTOFHid+100*seg+32, ad, td);
        HF2(HTOFHid+100*seg+33, au, ctu);  HF2(HTOFHid+100*seg+34, ad, ctd);
        HF1(HTOFHid+12, cmt);
        HF1(HTOFHid+12, cmt);
        if(m == 0){
          HF1(HTOFHid+100*seg+14, au); HF1(HTOFHid+100*seg+15, ad);
          HF1(HTOFHid+100*seg+16, de); HF1(HTOFHid+13, de);
        }
        if(de > 0.5){
          HF1(HTOFHid+15, seg-0.5);
          ++nh2;
        }
      }
    }
    Int_t nc = hodoAna->GetNClustersHTOF();
    HF1(HTOFHid+30, Double_t(nc));
    for(Int_t i=0; i<nc; ++i){
      HodoCluster *cluster = hodoAna->GetClusterHTOF(i);
      if(!cluster) continue;
      Int_t cs = cluster->ClusterSize();
      Double_t ms = cluster->MeanSeg()+1;
      Double_t cmt = cluster->CMeanTime();
      Double_t de = cluster->DeltaE();
      HF1(HTOFHid+31, Double_t(cs));
      HF1(HTOFHid+32, ms-0.5);
      HF1(HTOFHid+33, cmt); HF1(HTOFHid+34, de);
    }
  }

  // TOF
  hodoAna->DecodeTOFHits(rawData);
  {
    Int_t nh = hodoAna->GetNHitsTOF();
    HF1(TOFHid+10, Double_t(nh));
    Int_t nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      Hodo2Hit *hit = hodoAna->GetHitTOF(i);
      if(!hit) continue;
      Int_t seg = hit->SegmentId()+1;

      Int_t n_mhit = hit->GetNumOfHit();
      for(Int_t m=0; m<n_mhit; ++m){
        HF1(TOFHid+11, seg-0.5);

        Double_t au  = hit->GetAUp(),   ad  = hit->GetADown();
        Double_t tu  = hit->GetTUp(),   td  = hit->GetTDown();
        Double_t ctu = hit->GetCTUp(),  ctd = hit->GetCTDown();
        Double_t mt  = hit->MeanTime(), cmt = hit->CMeanTime();
        Double_t de  = hit->DeltaE();
        event.tofmt[seg-1][m] = mt;
        event.tofde[seg-1]    = de;
        HF1(TOFHid+100*seg+11, tu);      HF1(TOFHid+100*seg+12, td);
        HF1(TOFHid+100*seg+13, mt);
        HF1(TOFHid+100*seg+17, ctu);     HF1(TOFHid+100*seg+18, ctd);
        HF1(TOFHid+100*seg+19, cmt);     HF1(TOFHid+100*seg+20, ctu-ctd);
        //HF2(TOFHid+100*seg+21, tu, au);  HF2(TOFHid+100*seg+22, td, ad);
        //HF2(TOFHid+100*seg+23, ctu, au); HF2(TOFHid+100*seg+24, ctd, ad);
        HF1(TOFHid+12, cmt);

        dst.utTofSeg[seg-1][m] = tu;
        dst.dtTofSeg[seg-1][m] = td;
        dst.udeTofSeg[seg-1] = au;
        dst.ddeTofSeg[seg-1] = ad;

        if(m == 0){
          HF1(TOFHid+100*seg+14, au);    HF1(TOFHid+100*seg+15, ad);
          HF1(TOFHid+100*seg+16, de);    HF1(TOFHid+13, de);
        }

        if(de>0.5){
          HF1(TOFHid+15, seg-0.5);
          ++nh2;
        }
      }
    }

    HF1(TOFHid+14, Double_t(nh2));
    for(Int_t i1=0; i1<nh; ++i1){
      Hodo2Hit *hit1 = hodoAna->GetHitTOF(i1);
      if(!hit1 || hit1->DeltaE()<=0.5) continue;
      Int_t seg1 = hit1->SegmentId()+1;
      for(Int_t i2=0; i2<nh; ++i2){
        if(i1==i2) continue;
        Hodo2Hit *hit2 = hodoAna->GetHitTOF(i2);
        if(!hit2 || hit2->DeltaE()<=0.5) continue;
        Int_t seg2 = hit2->SegmentId()+1;

        if(1 == hit1->GetNumOfHit() && 1 == hit2->GetNumOfHit()){
          Double_t ct1 = hit1->CMeanTime(), ct2 = hit2->CMeanTime();
          HF2(TOFHid+21, seg1-0.5, seg2-0.5);
          HF2(TOFHid+22, ct1, ct2);
          HF1(TOFHid+23, ct2-ct1);
          if(std::abs(ct2-ct1)<3.0){
            HF2(TOFHid+24, seg1-0.5, seg2-0.5);
          }
        }//for(i2)
      }//for(i1)
    }

    Int_t nc = hodoAna->GetNClustersTOF();
    HF1(TOFHid+30, Double_t(nc));
    for(Int_t i=0; i<nc; ++i){
      HodoCluster *cluster = hodoAna->GetClusterTOF(i);
      if(!cluster) continue;
      Int_t cs = cluster->ClusterSize();
      Double_t ms  = cluster->MeanSeg()+1;
      Double_t cmt = cluster->CMeanTime();
      Double_t de  = cluster->DeltaE();
      HF1(TOFHid+31, Double_t(cs));
      HF1(TOFHid+32, ms-0.5);
      HF1(TOFHid+33, cmt); HF1(TOFHid+34, de);
    }
  }

  // WCSUM
  hodoAna->DecodeWCSUMHits(rawData);
  {
    Int_t nh = hodoAna->GetNHitsWCSUM();
    HF1(WCSUMHid+10, Double_t(nh));
    Int_t nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      Hodo1Hit *hit = hodoAna->GetHitWCSUM(i);
      if(!hit) continue;
      Int_t seg = hit->SegmentId()+1;
      for(Int_t m=0, n_mhit=hit->GetNumOfHit(); m<n_mhit; ++m){
        HF1(WCSUMHid+11, seg-0.5);
        Double_t au = hit->GetAUp();
        Double_t tu = hit->GetTUp();
        Double_t ctu = hit->GetCTUp();
        Double_t mt = hit->MeanTime(), cmt = hit->CMeanTime();
        Double_t de = hit->DeltaE();
        event.wcmt[seg-1][m] = mt;
        event.wcde[seg-1] = de;
        HF1(WCSUMHid+100*seg+11, tu);
        HF1(WCSUMHid+100*seg+13, mt);
        HF1(WCSUMHid+100*seg+17, ctu);
        HF1(WCSUMHid+100*seg+19, cmt);
        HF2(WCSUMHid+100*seg+21, tu, au);
        HF2(WCSUMHid+100*seg+23, ctu, au);
        HF1(WCSUMHid+12, cmt);
        if(m == 0){
          HF1(WCSUMHid+100*seg+14, au);
          HF1(WCSUMHid+100*seg+16, de); HF1(WCSUMHid+13, de);
        }
        if(de > 0.5){
          HF1(WCSUMHid+15, seg-0.5);
          ++nh2;
        }
      }
    }
    Int_t nc = hodoAna->GetNClustersWCSUM();
    HF1(WCSUMHid+30, Double_t(nc));
    for(Int_t i=0; i<nc; ++i){
      HodoCluster *cluster = hodoAna->GetClusterWCSUM(i);
      if(!cluster) continue;
      Int_t cs = cluster->ClusterSize();
      Double_t ms = cluster->MeanSeg()+1;
      Double_t cmt = cluster->CMeanTime();
      Double_t de = cluster->DeltaE();
      HF1(WCSUMHid+31, Double_t(cs));
      HF1(WCSUMHid+32, ms-0.5);
      HF1(WCSUMHid+33, cmt); HF1(WCSUMHid+34, de);
    }
  }



  ////////// SCH
  {
    hodoAna->DecodeSCHHits(rawData);
    int nh = hodoAna->GetNHitsSCH();
    int sch_nhits = 0;
    for(int i=0; i<nh; ++i){
      FiberHit* hit = hodoAna->GetHitSCH(i);
      if(!hit) continue;
      int mhit_l = hit->GetNLeading();
      int mhit_t = hit->GetNTrailing();
      int seg    = hit->SegmentId();
      event.sch_depth[seg] = mhit_l;
      int  prev     = 0;
      Bool_t hit_flag = false;

      for(int m=0; m<mhit_l; ++m){
        if(mhit_l > MaxDepth) break;
        double leading  = hit->GetLeading(m);
        if(leading==prev) continue;
        prev = leading;
        // HF1(SCHHid +3, leading);
        // HF2(SCHHid +5, seg, leading);
        // HF1(SCHHid +1000+seg+1, leading);

        event.sch_tdc[seg][m]      = leading;

        if(MinTdcSCH<leading && leading<MaxTdcSCH){
          hit_flag = true;
        }
      }// for(m)

      for(int m=0; m<mhit_t; ++m){
        if(mhit_t > MaxDepth) break;
        double trailing = hit->GetTrailing(m);
        event.sch_trailing[seg][m] = trailing;
      }// for(m)

      int mhit_pair = hit->GetNPair();
      for(int m=0; m<mhit_pair; ++m){
        if(mhit_pair > MaxDepth) break;

        // double time     = hit->GetTime(m);
        // double ctime    = hit->GetCTime(m);
        double width    = hit->GetWidth(m);

        // HF1(SCHHid +4, width);
        // HF2(SCHHid +6, seg, width);
        // HF1(SCHHid +21, time);
        // HF2(SCHHid +22, width, time);
        // HF1(SCHHid +31, ctime);
        // HF2(SCHHid +32, width, ctime);
        // HF1(SCHHid +2000+seg+1, width);
        // if(-10.<time && time<10.){
        //   HF2(SCHHid +3000+seg+1, width, time);
        //   HF2(SCHHid +4000+seg+1, width, ctime);
        // }

        event.sch_tot[seg][m]      = width;

      }//for(m)
      if(hit_flag){
        //	HF1(SCHHid +2, seg+0.5);
        event.sch_hitpat[sch_nhits++] = seg;
      }
    }
    //    HF1(SCHHid +1, sch_nhits);
    event.sch_nhits = sch_nhits;

    //Fiber Cluster
#if TimeCut
    hodoAna->TimeCutSCH(MinTimeSCH, MaxTimeSCH);
#endif
    int ncl = hodoAna->GetNClustersSCH();
    if(ncl > NumOfSegSCH){
      // std::cout << "#W SCH too much number of clusters" << std::endl;
      ncl = NumOfSegSCH;
    }
    event.sch_ncl = ncl;
    //HF1(SCHHid +101, ncl);
    for(int i=0; i<ncl; ++i){
      FiberCluster *cl = hodoAna->GetClusterSCH(i);
      if(!cl) continue;
      double clsize = cl->ClusterSize();
      double ctime  = cl->CMeanTime();
      double ctot   = cl->Width();
      double pos    = cl->MeanPosition();
      event.sch_clsize[i] = clsize;
      event.sch_ctime[i]  = ctime;
      event.sch_ctot[i]   = ctot;
      event.sch_clpos[i]  = pos;
      // HF1(SCHHid +102, clsize);
      // HF1(SCHHid +103, ctime);
      // HF1(SCHHid +104, ctot);
      // HF2(SCHHid +105, ctot, ctime);
      // HF1(SCHHid +106, pos);
    }
  }

  hodoAna->DecodeBVHHits(rawData);
  hodoAna->DecodeLACHits(rawData);

  ////////// Dst
  {
    Int_t nc = hodoAna->GetNClustersBH1();
    dst.nhBh1 = nc;
    for(Int_t i=0; i<nc; ++i){
      HodoCluster *cl = hodoAna->GetClusterBH1(i);
      if(!cl) continue;
      dst.csBh1[i]  = cl->ClusterSize();
      dst.Bh1Seg[i] = cl->MeanSeg()+1;
      dst.tBh1[i]   = cl->CMeanTime();
      dst.dtBh1[i]  = cl->TimeDif();
      dst.deBh1[i]  = cl->DeltaE();

      Int_t nc2 = hodoAna->GetNClustersBH2();
      for(Int_t i2=0; i2<nc2; ++i2){
        BH2Cluster *cl2 = hodoAna->GetClusterBH2(i2);
        if(!cl2) continue;
        Double_t btof = cl->MeanTime() - dst.Time0;
        Double_t cbtof = cl->CMeanTime() - dst.CTime0;
        dst.btof[i] = btof;
        dst.cbtof[i] = cbtof;
      }
    }
  }

  {
    Int_t nc = hodoAna->GetNClustersBH2();
    dst.nhBh2 = nc;
    for(Int_t i=0; i<nc; ++i){
      BH2Cluster *cl = hodoAna->GetClusterBH2(i);
      if(!cl) continue;
      dst.csBh2[i]  = cl->ClusterSize();
      dst.Bh2Seg[i] = cl->MeanSeg()+1;
      dst.tBh2[i]   = cl->CMeanTime();
      dst.t0Bh2[i]  = cl->CTime0();
      dst.dtBh2[i]  = cl->TimeDif();
      dst.deBh2[i]  = cl->DeltaE();
#if HodoHitPos
      dst.posBh2[i]  = 0.5*PropVelBH2*cl->TimeDif();
#endif
    }
  }

  {
    Int_t nc = hodoAna->GetNClustersBAC();
    dst.nhBac = nc;
    for(Int_t i=0; i<nc; ++i){
      HodoCluster *cl = hodoAna->GetClusterBAC(i);
      if(!cl) continue;
      dst.csBac[i]  = cl->ClusterSize();
      dst.BacSeg[i] = cl->MeanSeg()+1;
      dst.tBac[i]   = cl->CMeanTime();
      dst.deBac[i]  = cl->DeltaE();
    }
  }

  {
    Int_t nc = hodoAna->GetNClustersHTOF();
    dst.nhHtof = nc;
    for(Int_t i=0; i<nc; ++i){
      HodoCluster *cl = hodoAna->GetClusterHTOF(i);
      if(!cl) continue;
      dst.csHtof[i] = cl->ClusterSize();
      dst.HtofSeg[i] = cl->MeanSeg()+1;
      dst.tHtof[i] = cl->CMeanTime();
      dst.dtHtof[i] = cl->TimeDif();
      dst.deHtof[i] = cl->DeltaE();
#if HodoHitPos
      dst.posHtof[i]  = 0.5*PropVelHTOF*cl->TimeDif();
#endif
    }
  }

  {
    Int_t nc = hodoAna->GetNClustersBVH();
    dst.nhBvh = nc;
    for(Int_t i=0; i<nc; ++i){
      HodoCluster *cl = hodoAna->GetClusterBVH(i);
      if(!cl) continue;
      dst.csBvh[i]  = cl->ClusterSize();
      dst.BvhSeg[i] = cl->MeanSeg()+1;
      dst.tBvh[i]   = cl->CMeanTime();
    }
  }

  {
    Int_t nc = hodoAna->GetNClustersTOF();
    dst.nhTof = nc;
    for(Int_t i=0; i<nc; ++i){
      HodoCluster *cl = hodoAna->GetClusterTOF(i);
      if(!cl) continue;
      dst.csTof[i]  = cl->ClusterSize();
      dst.TofSeg[i] = cl->MeanSeg()+1;
      dst.tTof[i]   = cl->CMeanTime();
      dst.dtTof[i]  = cl->TimeDif();
      dst.deTof[i]  = cl->DeltaE();
    }
  }

  {
    Int_t nc = hodoAna->GetNClustersLAC();
    dst.nhLac = nc;
    for(Int_t i=0; i<nc; ++i){
      HodoCluster *cl = hodoAna->GetClusterLAC(i);
      if(!cl) continue;
      dst.csLac[i]  = cl->ClusterSize();
      dst.LacSeg[i] = cl->MeanSeg()+1;
      dst.tLac[i]   = cl->CMeanTime();
    }
  }

#if 0
  // BH1 (for parameter tuning)
  if(dst.Time0Seg==4){
    const HodoRHitContainer& cont = rawData->GetBH1RawHC();
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;

      //Up
      {
        Int_t n_mhit = hit->GetSizeTdcUp();
        for(Int_t m=0; m<n_mhit; ++m){
          Int_t T = hit->GetTdcUp(m);
          if(T > 0) HF1(BH1Hid +100*seg +9, Double_t(T));
        }// for(m)
      }

      //Down
      {
        Int_t n_mhit = hit->GetSizeTdcDown();
        for(Int_t m=0; m<n_mhit; ++m){
          Int_t T = hit->GetTdcDown(m);
          if(T > 0) HF1(BH1Hid +100*seg +10, Double_t(T));
        }// for(m)
      }
    }
  }// (Raw BH1 TDC for parameter tuning)
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
UserHodoscope::ProcessingEnd()
{
  tree->Fill();
  hodo->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserHodoscope;
}

//_____________________________________________________________________________
namespace
{
const Int_t    NbinAdc = 4096;
const Double_t MinAdc  =    0.;
const Double_t MaxAdc  = 4096.;

const Int_t    NbinTdc = 4096;
const Double_t MinTdc  =    0.;
const Double_t MaxTdc  = 4096.;

const Int_t    NbinTdcHr = 1e6/10;
const Double_t MinTdcHr  =  0.;
const Double_t MaxTdcHr  = 1e6;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  HB1(1, "Status", 20, 0., 20.);
  HB1(10, "Trigger HitPat", NumOfSegTrig, 0., Double_t(NumOfSegTrig));
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1(10+i+1, Form("Trigger Trig %d", i+1), 0x1000, 0, 0x1000);
  }

  // BH1
  // Rawdata
  HB1(BH1Hid +0, "#Hits BH1",        NumOfSegBH1+1, 0., Double_t(NumOfSegBH1+1));
  HB1(BH1Hid +1, "Hitpat BH1",       NumOfSegBH1,   0., Double_t(NumOfSegBH1));
  HB1(BH1Hid +2, "#Hits BH1(Tor)",   NumOfSegBH1+1, 0., Double_t(NumOfSegBH1+1));
  HB1(BH1Hid +3, "Hitpat BH1(Tor)",  NumOfSegBH1,   0., Double_t(NumOfSegBH1));
  HB1(BH1Hid +4, "#Hits BH1(Tand)",  NumOfSegBH1+1, 0., Double_t(NumOfSegBH1+1));
  HB1(BH1Hid +5, "Hitpat BH1(Tand)", NumOfSegBH1,   0., Double_t(NumOfSegBH1));

  for(Int_t i=1; i<=NumOfSegBH1; ++i){
    TString title1 = Form("BH1-%d UpAdc", i);
    TString title2 = Form("BH1-%d DownAdc", i);
    TString title3 = Form("BH1-%d UpTdc", i);
    TString title4 = Form("BH1-%d DownTdc", i);
    TString title5 = Form("BH1-%d UpAdc(w Tdc)", i);
    TString title6 = Form("BH1-%d DownAdc(w Tdc)", i);
    TString title7 = Form("BH1-%d UpAdc(w/o Tdc)", i);
    TString title8 = Form("BH1-%d DownAdc(w/o Tdc)", i);
    TString title9 = Form("BH1-%d UpTdc (Time0Seg==4)", i);
    TString title10= Form("BH1-%d DownTdc (Time0Seg==4)", i);
    HB1(BH1Hid +100*i +1, title1, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH1Hid +100*i +2, title2, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH1Hid +100*i +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(BH1Hid +100*i +4, title4, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(BH1Hid +100*i +5, title5, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH1Hid +100*i +6, title6, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH1Hid +100*i +7, title7, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH1Hid +100*i +8, title8, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH1Hid +100*i +9, title9, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(BH1Hid +100*i +10,title10,NbinTdcHr, MinTdcHr, MaxTdcHr);
  }

  //BH1 Normalized
  HB1(BH1Hid +10, "#Hits BH1[Hodo]",  NumOfSegBH1+1, 0., Double_t(NumOfSegBH1+1));
  HB1(BH1Hid +11, "Hitpat BH1[Hodo]", NumOfSegBH1,   0., Double_t(NumOfSegBH1));
  HB1(BH1Hid +12, "CMeanTime BH1", 200, -10., 10.);
  HB1(BH1Hid +13, "dE BH1", 200, -0.5, 4.5);
  HB1(BH1Hid +14, "#Hits BH1[HodoGood]",  NumOfSegBH1+1, 0., Double_t(NumOfSegBH1+1));
  HB1(BH1Hid +15, "Hitpat BH1[HodoGood]", NumOfSegBH1,   0., Double_t(NumOfSegBH1));
  HB1(BH1Hid +16, "CMeanTime BH1[HodoGood]", 200, -10., 10.);

  for(Int_t i=1; i<=NumOfSegBH1; ++i){
    TString title11 = Form("BH1-%d Up Time", i);
    TString title12 = Form("BH1-%d Down Time", i);
    TString title13 = Form("BH1-%d MeanTime", i);
    TString title14 = Form("BH1-%d Up dE", i);
    TString title15 = Form("BH1-%d Down dE", i);
    TString title16 = Form("BH1-%d dE", i);
    TString title17 = Form("BH1-%d Up CTime", i);
    TString title18 = Form("BH1-%d Down CTime", i);
    TString title19 = Form("BH1-%d CMeanTime", i);
    TString title20 = Form("BH1-%d Tup-Tdown", i);
    TString title21 = Form("BH1-%d Up dE%%Time", i);
    TString title22 = Form("BH1-%d Down dE%%Time", i);
    TString title23 = Form("BH1-%d Up dE%%CTime", i);
    TString title24 = Form("BH1-%d Down dE%%CTime", i);
    TString title31 = Form("BH1-%d Up Time%%dE", i);
    TString title32 = Form("BH1-%d Down Time%%dE", i);
    TString title33 = Form("BH1-%d Up CTime%%dE", i);
    TString title34 = Form("BH1-%d Down CTIme%%dE", i);
    HB1(BH1Hid +100*i +11, title11, 200, -10., 10.);
    HB1(BH1Hid +100*i +12, title12, 200, -10., 10.);
    HB1(BH1Hid +100*i +13, title13, 200, -10., 10.);
    HB1(BH1Hid +100*i +14, title14, 200, -0.5, 4.5);
    HB1(BH1Hid +100*i +15, title15, 200, -0.5, 4.5);
    HB1(BH1Hid +100*i +16, title16, 200, -0.5, 4.5);
    HB1(BH1Hid +100*i +17, title17, 200, -10., 10.);
    HB1(BH1Hid +100*i +18, title18, 200, -10., 10.);
    HB1(BH1Hid +100*i +19, title19, 200, -10., 10.);
    HB1(BH1Hid +100*i +20, title20, 200, -5.0, 5.0);
    HB2(BH1Hid +100*i +21, title21, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH1Hid +100*i +22, title22, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH1Hid +100*i +23, title23, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH1Hid +100*i +24, title24, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH1Hid +100*i +31, title31, 100, -0.5, 4.5,1000,-10.,10.);
    HB2(BH1Hid +100*i +32, title32, 100, -0.5, 4.5,1000,-10.,10.);
    HB2(BH1Hid +100*i +33, title33, 100, -0.5, 4.5,1000,-10.,10.);
    HB2(BH1Hid +100*i +34, title34, 100, -0.5, 4.5,1000,-10.,10.);

    //For BH1vsBH2 Correlation
    HB2(BH1Hid +100*i +2200 +21, Form("BH1-%d BH2 Up MT%%MT", i),
        100, -10., 10., 100, -10., 10.);
    HB2(BH1Hid +100*i +2200 +22, Form("BH1-%d BH2 Down MT%%MT", i),
        100, -10., 10., 100, -10., 10.);
    HB2(BH1Hid +100*i +2200 +23, Form("BH1-%d BH2 MeanTime MT%%MT", i),
        100, -10., 10., 100, -10., 10.);
  }
  for(Int_t i=1; i<=NumOfSegBH1; ++i){
    TString title1 = Form("BH1-%d Up Time [BH2]", i);
    HB1(BH1Hid +1100 +100*i +31, title1, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +41, title1, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +51, title1, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +61, title1, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +71, title1, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +81, title1, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +91, title1, 200, -10., 10.);
    TString title2 = Form("BH1-%d Down Time [BH2]", i);
    HB1(BH1Hid +1100 +100*i +32, title2, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +42, title2, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +52, title2, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +62, title2, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +72, title2, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +82, title2, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +92, title2, 200, -10., 10.);
    TString title3 = Form("BH1-%d MeanTime [BH2]", i);
    HB1(BH1Hid +1100 +100*i +33, title3, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +43, title3, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +53, title3, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +63, title3, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +73, title3, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +83, title3, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +93, title3, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +103, title3, 200, -10., 10.);
    TString title4 = Form("BH1-%d MeanTime-BH2MeamTime", i);
    HB1(BH1Hid +1100 +100*i +34, title4, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +44, title4, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +54, title4, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +64, title4, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +74, title4, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +84, title4, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +94, title4, 200, -10., 10.);
    HB1(BH1Hid +1100 +100*i +104, title4, 200, -10., 10.);
    HB1(BH1Hid +100*i +2200 +114, title4, 200, -10., 10.);
  }

  HB2(BH1Hid +21, "BH1HitPat%BH1HitPat[HodoGood]", NumOfSegBH1,   0., Double_t(NumOfSegBH1),
      NumOfSegBH1,   0., Double_t(NumOfSegBH1));
  HB2(BH1Hid +22, "CMeanTimeBH1%CMeanTimeBH1[HodoGood]",
      100, -5., 5., 100, -5., 5.);
  HB1(BH1Hid +23, "TDiff BH1[HodoGood]", 200, -10., 10.);
  HB2(BH1Hid +24, "BH1HitPat%BH1HitPat[HodoGood2]", NumOfSegBH1,   0., Double_t(NumOfSegBH1),
      NumOfSegBH1,   0., Double_t(NumOfSegBH1));
  HB1(BH1Hid +30, "#Clusters BH1", NumOfSegBH1+1, 0., Double_t(NumOfSegBH1+1));
  HB1(BH1Hid +31, "ClusterSize BH1", 5, 0., 5.);
  HB1(BH1Hid +32, "HitPat Cluster BH1", 2*NumOfSegBH1, 0., Double_t(NumOfSegBH1));
  HB1(BH1Hid +33, "CMeamTime Cluster BH1", 200, -10., 10.);
  HB1(BH1Hid +34, "DeltaE Cluster BH1", 100, -0.5, 4.5);
  HB1(BH1Hid +35, "#Clusters BH1(AdcGood)", NumOfSegBH1+1, 0., Double_t(NumOfSegBH1+1));
  HB1(BH1Hid +36, "CMeamTime Cluster BH1(AdcGood)", 200, -10., 10.);

  HB2(BH1Hid +41, "BH1ClP%BH1ClP",  NumOfSegBH1,   0., Double_t(NumOfSegBH1),
      NumOfSegBH1,   0., Double_t(NumOfSegBH1));
  HB2(BH1Hid +42, "CMeanTimeBH1%CMeanTimeBH1[Cluster]",
      100, -5., 5., 100, -5., 5.);
  HB1(BH1Hid +43, "TDiff BH1[Cluster]", 200, -10., 10.);
  HB2(BH1Hid +44, "BH1ClP%BH1ClP[AdcGood]",  NumOfSegBH1,   0., Double_t(NumOfSegBH1),
      NumOfSegBH1,   0., Double_t(NumOfSegBH1));

  // BH2
  HB1(BH2Hid +0, "#Hits BH2",        NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +1, "Hitpat BH2",       NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB1(BH2Hid +2, "#Hits BH2(Tor)",   NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +3, "Hitpat BH2(Tor)",  NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB1(BH2Hid +4, "#Hits BH2(Tand)",  NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +5, "Hitpat BH2(Tand)", NumOfSegBH2,   0., Double_t(NumOfSegBH2));

  for(Int_t i=1; i<=NumOfSegBH2; ++i){
    TString title1 = Form("BH2-%d UpAdc", i);
    TString title2 = Form("BH2-%d DownAdc", i);
    TString title3 = Form("BH2-%d UpTdc", i);
    TString title4 = Form("BH2-%d DownTdc", i);
    TString title5 = Form("BH2-%d UpAdc(w Tdc)", i);
    TString title6 = Form("BH2-%d DownAdc(w Tdc)", i);
    TString title7 = Form("BH2-%d UpAdc(w/o Tdc)", i);
    TString title8 = Form("BH2-%d DownAdc(w/o Tdc)", i);
    HB1(BH2Hid +100*i +1, title1, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH2Hid +100*i +2, title2, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH2Hid +100*i +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(BH2Hid +100*i +4, title4, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(BH2Hid +100*i +5, title5, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH2Hid +100*i +6, title6, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH2Hid +100*i +7, title7, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH2Hid +100*i +8, title8, NbinAdc,   MinAdc,   MaxAdc);
  }
  HB1(BH2Hid +10, "#Hits BH2[Hodo]",  NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +11, "Hitpat BH2[Hodo]", NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB1(BH2Hid +12, "CMeanTime BH2", 200, -10., 10.);
  HB1(BH2Hid +13, "dE BH2", 200, -0.5, 4.5);
  HB1(BH2Hid +14, "#Hits BH2[HodoGood]",  NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +15, "Hitpat BH2[HodoGood]", NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB1(BH2Hid +16, "CMeanTime BH2[HodoGood]", 200, -10., 10.);

  for(Int_t i=1; i<=NumOfSegBH2; ++i){
    TString title11 = Form("BH2-%d Up Time", i);
    TString title12 = Form("BH2-%d Down Time", i);
    TString title13 = Form("BH2-%d MeanTime", i);
    TString title14 = Form("BH2-%d Up dE", i);
    TString title15 = Form("BH2-%d Down dE", i);
    TString title16 = Form("BH2-%d dE", i);
    TString title17 = Form("BH2-%d Up CTime", i);
    TString title18 = Form("BH2-%d Down CTime", i);
    TString title19 = Form("BH2-%d CMeanTime", i);
    TString title20 = Form("BH2-%d Tup-Tdown", i);
    TString title21 = Form("BH2-%d Up Time0", i);
    TString title22 = Form("BH2-%d Down Time0", i);
    TString title23 = Form("BH2-%d Up CTime", i);
    TString title24 = Form("BH2-%d Down CTime", i);
    TString title25 = Form("BH2-%d MeanTime0", i);
    TString title26 = Form("BH2-%d CMeanTime0", i);
    TString title27 = Form("BH2-%d Up dE%%Time", i);
    TString title28 = Form("BH2-%d Down dE%%Time", i);
    TString title29 = Form("BH2-%d Up dE%%CTime", i);
    TString title30 = Form("BH2-%d Down dE%%CTime", i);
    HB1(BH2Hid +100*i +11, title11, 200, -10., 10.);
    HB1(BH2Hid +100*i +12, title12, 200, -10., 10.);
    HB1(BH2Hid +100*i +13, title13, 200, -10., 10.);
    HB1(BH2Hid +100*i +14, title14, 200, -0.5, 4.5);
    HB1(BH2Hid +100*i +15, title15, 200, -0.5, 4.5);
    HB1(BH2Hid +100*i +16, title16, 200, -0.5, 4.5);
    HB1(BH2Hid +100*i +17, title17, 200, -10., 10.);
    HB1(BH2Hid +100*i +18, title18, 200, -10., 10.);
    HB1(BH2Hid +100*i +19, title19, 200, -10., 10.);
    HB1(BH2Hid +100*i +20, title20, 200, -5.0, 5.0);
    HB1(BH2Hid +100*i +21, title21, 200, -10., 10.);
    HB1(BH2Hid +100*i +22, title22, 200, -10., 10.);
    HB1(BH2Hid +100*i +23, title23, 200, -10., 10.);
    HB1(BH2Hid +100*i +24, title24, 200, -10., 10.);
    HB1(BH2Hid +100*i +25, title25, 200, -10., 10.);
    HB1(BH2Hid +100*i +26, title26, 200, -10., 10.);
    HB2(BH2Hid +100*i +27, title27, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH2Hid +100*i +28, title28, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH2Hid +100*i +29, title29, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH2Hid +100*i +30, title30, 100, -10., 10., 100, -0.5, 4.5);
  }

  HB2(BH2Hid +21, "BH2HitPat%BH2HitPat[HodoGood]", NumOfSegBH2,   0., Double_t(NumOfSegBH2),
      NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB2(BH2Hid +22, "CMeanTimeBH2%CMeanTimeBH2[HodoGood]",
      100, -2.5, 2.5, 100, -2.5, 2.5);
  HB1(BH2Hid +23, "TDiff BH2[HodoGood]", 200, -10., 10.);
  HB2(BH2Hid +24, "BH2HitPat%BH2HitPat[HodoGood2]", NumOfSegBH2,   0., Double_t(NumOfSegBH2),
      NumOfSegBH2,   0., Double_t(NumOfSegBH2));

  HB1(BH2Hid +30, "#Clusters BH2", NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +31, "ClusterSize BH2", 5, 0., 5.);
  HB1(BH2Hid +32, "HitPat Cluster BH2", 2*NumOfSegBH2, 0., Double_t(NumOfSegBH2));
  HB1(BH2Hid +33, "CMeamTime Cluster BH2", 200, -10., 10.);
  HB1(BH2Hid +34, "DeltaE Cluster BH2", 100, -0.5, 4.5);
  HB1(BH2Hid +35, "#Clusters BH2(ADCGood)", NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +36, "CMeamTime Cluster BH2(ADCGood)", 200, -10., 10.);

  HB2(BH2Hid +41, "BH2ClP%BH2ClP", NumOfSegBH2,   0., Double_t(NumOfSegBH2),
      NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB2(BH2Hid +42, "CMeanTimeBH2%CMeanTimeBH2[Cluster]",
      100, -2.5, 2.5, 100, -2.5, 2.5);
  HB1(BH2Hid +43, "TDiff BH2[Cluster]", 200, -10., 10.);
  HB2(BH2Hid +44, "BH2ClP%BH2ClP(ADCGood)", NumOfSegBH2,   0., Double_t(NumOfSegBH2),
      NumOfSegBH2,   0., Double_t(NumOfSegBH2));

  HB1(201, "TimeDif BH1-BH2", 200, -10., 10.);
  HB2(202, "SegBH2%SegBH1",NumOfSegBH1,   0., Double_t(NumOfSegBH1),
      NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  //For BH1vsBH2 Corr
  HB2(203, "MTBH2%MTBH1", 200, -10., 10., 200, -10., 10.);
  HB1(204, "MTBH1", 200, -10., 10.);
  HB1(205, "MTBH2", 200, -10., 10.);

  HB1(211, "TimeDif BH1-BH2(GoodAdc)", 200, -10., 10.);
  HB2(212, "SegBH2%SegBH1(GoodAdc)",NumOfSegBH1,   0., Double_t(NumOfSegBH1),
      NumOfSegBH2,   0., Double_t(NumOfSegBH2));

  // BH1-BH2 PHC
  for(Int_t i=1; i<=NumOfSegBH1; ++i){
    TString title1 = Form("BH1-%dU CT-TOF%%dE", i);
    TString title2 = Form("BH1-%dD CT-TOF%%dE", i);
    TString title3 = Form("BH1-%dU  T-TOF%%dE", i);
    TString title4 = Form("BH1-%dD  T-TOF%%dE", i);
    HB2(BH1Hid +100*i +81, title1, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH1Hid +100*i +82, title2, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH1Hid +100*i +83, title3, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH1Hid +100*i +84, title4, 200, -0.5, 4.5, 200, -10., 10.);
  }
  for(Int_t i=1; i<=NumOfSegBH2; ++i){
    TString title1 = Form("BH2-%dU CT-TOF%%dE", i);
    TString title2 = Form("BH2-%dD CT-TOF%%dE", i);
    TString title3 = Form("BH2-%dU  T-TOF%%dE", i);
    TString title4 = Form("BH2-%dD  T-TOF%%dE", i);
    HB2(BH2Hid +100*i +81, title1, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH2Hid +100*i +82, title2, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH2Hid +100*i +83, title3, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH2Hid +100*i +84, title4, 200, -0.5, 4.5, 200, -10., 10.);
  }

  // BTOF
  HB1(100, "BH2 MeanTime0", 400, -4, 4);
  HB1(101, "CTime0", 400, -4, 4);
  HB1(102, "CBtof0", 400, -4, 4);
  HB2(103, "CBtof0%deTime0", 200, 0, 4, 200, -4, 4);
  HB2(104, "CBtof0%deBTof0", 200, 0, 4, 200, -4, 4);

  // BAC
  HB1(BACHid +0, "#Hits BAC",        NumOfSegBAC+1, 0., Double_t(NumOfSegBAC+1));
  HB1(BACHid +1, "Hitpat BAC",       NumOfSegBAC,   0., Double_t(NumOfSegBAC));
  HB1(BACHid +2, "#Hits BAC(Tor)",   NumOfSegBAC+1, 0., Double_t(NumOfSegBAC+1));
  HB1(BACHid +3, "Hitpat BAC(Tor)",  NumOfSegBAC,   0., Double_t(NumOfSegBAC));
  HB1(BACHid +4, "#Hits BAC(Tand)",  NumOfSegBAC+1, 0., Double_t(NumOfSegBAC+1));
  HB1(BACHid +5, "Hitpat BAC(Tand)", NumOfSegBAC,   0., Double_t(NumOfSegBAC));
  for(Int_t i=1; i<=NumOfSegBAC; ++i){
    TString title1 = Form("BAC-%d UpAdc", i);
    TString title3 = Form("BAC-%d UpTdc", i);
    TString title5 = Form("BAC-%d UpAdc(w Tdc)", i);
    TString title7 = Form("BAC-%d UpAdc(w/o Tdc)", i);
    HB1(BACHid +100*i +1, title1, NbinTdc, MinTdc, MaxTdc);
    HB1(BACHid +100*i +3, title3, NbinTdc, MinTdc, MaxTdc);
    HB1(BACHid +100*i +5, title5, NbinAdc, MinAdc, MaxAdc);
    HB1(BACHid +100*i +7, title7, NbinAdc, MinAdc, MaxAdc);
  }
  HB1(BACHid +10, "#Hits BAC[Hodo]",     NumOfSegBAC+1, 0., Double_t(NumOfSegBAC+1));
  HB1(BACHid +11, "Hitpat BAC[Hodo]",    NumOfSegBAC,   0., Double_t(NumOfSegBAC));
  for(Int_t i=1; i<=NumOfSegBAC; ++i){
    TString title1 = Form("BAC-%d Time", i);
    TString title3 = Form("BAC-%d dE", i);
    TString title5 = Form("BAC-%d CTime", i);
    HB1(BACHid +100*i +11, title1, 500, -5., 45.);
    HB1(BACHid +100*i +12, title3, 200, -0.5, 4.5);
    HB1(BACHid +100*i +13, title5, 500, -5., 45.);
  }

  // HTOF
  HB1(HTOFHid +0, "#Hits HTOF",        NumOfSegHTOF+1, 0., Double_t(NumOfSegHTOF+1));
  HB1(HTOFHid +1, "Hitpat HTOF",       NumOfSegHTOF,   0., Double_t(NumOfSegHTOF));
  HB1(HTOFHid +2, "#Hits HTOF(Tor)",   NumOfSegHTOF+1, 0., Double_t(NumOfSegHTOF+1));
  HB1(HTOFHid +3, "Hitpat HTOF(Tor)",  NumOfSegHTOF,   0., Double_t(NumOfSegHTOF));
  HB1(HTOFHid +4, "#Hits HTOF(Tand)",  NumOfSegHTOF+1, 0., Double_t(NumOfSegHTOF+1));
  HB1(HTOFHid +5, "Hitpat HTOF(Tand)", NumOfSegHTOF,   0., Double_t(NumOfSegHTOF));

  for(Int_t i=1; i<=NumOfSegHTOF; ++i){
    TString title1 = Form("HTOF-%d UpAdc", i);
    TString title2 = Form("HTOF-%d DownAdc", i);
    TString title3 = Form("HTOF-%d UpTdc", i);
    TString title4 = Form("HTOF-%d DownTdc", i);
    TString title5 = Form("HTOF-%d UpAdc(w Tdc)", i);
    TString title6 = Form("HTOF-%d DownAdc(w Tdc)", i);
    TString title7 = Form("HTOF-%d UpAdc(w/o Tdc)", i);
    TString title8 = Form("HTOF-%d DownAdc(w/o Tdc)", i);
    HB1(HTOFHid +100*i +1, title1, NbinAdc, MinAdc, MaxAdc);
    HB1(HTOFHid +100*i +2, title2, NbinAdc, MinAdc, MaxAdc);
    HB1(HTOFHid +100*i +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(HTOFHid +100*i +4, title4, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(HTOFHid +100*i +5, title5, NbinAdc, MinAdc, MaxAdc);
    HB1(HTOFHid +100*i +6, title6, NbinAdc, MinAdc, MaxAdc);
    HB1(HTOFHid +100*i +7, title7, NbinAdc, MinAdc, MaxAdc);
    HB1(HTOFHid +100*i +8, title8, NbinAdc, MinAdc, MaxAdc);
  }

  HB1(HTOFHid +10, "#Hits Htof[Hodo]",  NumOfSegHTOF+1, 0., Double_t(NumOfSegHTOF+1));
  HB1(HTOFHid +11, "Hitpat Htof[Hodo]", NumOfSegHTOF,   0., Double_t(NumOfSegHTOF));
  HB1(HTOFHid +12, "CMeanTime Htof", 500, -5., 45.);
  HB1(HTOFHid +13, "dE Htof", 200, -0.5, 4.5);
  HB1(HTOFHid +14, "#Hits Htof[HodoGood]",  NumOfSegHTOF+1, 0., Double_t(NumOfSegHTOF+1));
  HB1(HTOFHid +15, "Hitpat Htof[HodoGood]", NumOfSegHTOF,   0., Double_t(NumOfSegHTOF));

  for(Int_t i=1; i<=NumOfSegHTOF; ++i){
    TString title11 = Form("HTOF-%d Up Time", i);
    TString title12 = Form("HTOF-%d Down Time", i);
    TString title13 = Form("HTOF-%d MeanTime", i);
    TString title14 = Form("HTOF-%d Up dE", i);
    TString title15 = Form("HTOF-%d Down dE", i);
    TString title16 = Form("HTOF-%d dE", i);
    TString title17 = Form("HTOF-%d Up CTime", i);
    TString title18 = Form("HTOF-%d Down CTime", i);
    TString title19 = Form("HTOF-%d CMeanTime", i);
    TString title20 = Form("HTOF-%d Tup-Tdown", i);
    TString title21 = Form("HTOF-%d dE (w/ HTOF-HT)", i);
    TString title31 = Form("HTOF-%d Up Time%%dE", i);
    TString title32 = Form("HTOF-%d Down Time%%dE", i);
    TString title33 = Form("HTOF-%d Up CTime%%dE", i);
    TString title34 = Form("HTOF-%d Down CTIme%%dE", i);
    HB1(HTOFHid +100*i +11, title11, 500, -5., 45.);
    HB1(HTOFHid +100*i +12, title12, 500, -5., 45.);
    HB1(HTOFHid +100*i +13, title13, 500, -5., 45.);
    HB1(HTOFHid +100*i +14, title14, 200, -0.5, 4.5);
    HB1(HTOFHid +100*i +15, title15, 200, -0.5, 4.5);
    HB1(HTOFHid +100*i +16, title16, 200, -0.5, 4.5);
    HB1(HTOFHid +100*i +17, title17, 500, -5., 45.);
    HB1(HTOFHid +100*i +18, title18, 500, -5., 45.);
    HB1(HTOFHid +100*i +19, title19, 500, -5., 45.);
    HB1(HTOFHid +100*i +20, title20, 200, -10.0, 10.0);
    HB1(HTOFHid +100*i +21, title21, 200, -0.5, 4.5);
    HB2(HTOFHid +100*i +31, title31, 100, -0.5, 4.5,1000,-55.,-35.);
    HB2(HTOFHid +100*i +32, title32, 100, -0.5, 4.5,1000,-55.,-35.);
    HB2(HTOFHid +100*i +33, title33, 100, -0.5, 4.5,1000,-55,-35.);
    HB2(HTOFHid +100*i +34, title34, 100, -0.5, 4.5,1000,-55.,-35.);
  }

  HB2(HTOFHid +21, "HtofHitPat%HtofHitPat[HodoGood]", NumOfSegHTOF,   0., Double_t(NumOfSegHTOF),
      NumOfSegHTOF,   0., Double_t(NumOfSegHTOF));
  HB2(HTOFHid +22, "CMeanTimeHtof%CMeanTimeHtof[HodoGood]",
      120, 10., 40., 120, 10., 40.);
  HB1(HTOFHid +23, "TDiff Htof[HodoGood]", 200, -10., 10.);
  HB2(HTOFHid +24, "HtofHitPat%HtofHitPat[HodoGood2]", NumOfSegHTOF,   0., Double_t(NumOfSegHTOF),
      NumOfSegHTOF,   0., Double_t(NumOfSegHTOF));

  HB1(HTOFHid +30, "#Clusters Htof", NumOfSegHTOF+1, 0., Double_t(NumOfSegHTOF+1));
  HB1(HTOFHid +31, "ClusterSize Htof", 5, 0., 5.);
  HB1(HTOFHid +32, "HitPat Cluster Htof", 2*NumOfSegHTOF, 0., Double_t(NumOfSegHTOF));
  HB1(HTOFHid +33, "CMeamTime Cluster Htof", 500, -5., 45.);
  HB1(HTOFHid +34, "DeltaE Cluster Htof", 100, -0.5, 4.5);

  // TOF
  HB1(TOFHid +0, "#Hits TOF",        NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +1, "Hitpat TOF",       NumOfSegTOF,   0., Double_t(NumOfSegTOF));
  HB1(TOFHid +2, "#Hits TOF(Tor)",   NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +3, "Hitpat TOF(Tor)",  NumOfSegTOF,   0., Double_t(NumOfSegTOF));
  HB1(TOFHid +4, "#Hits TOF(Tand)",  NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +5, "Hitpat TOF(Tand)", NumOfSegTOF,   0., Double_t(NumOfSegTOF));

  for(Int_t i=1; i<=NumOfSegTOF; ++i){
    TString title1 = Form("TOF-%d UpAdc", i);
    TString title2 = Form("TOF-%d DownAdc", i);
    TString title3 = Form("TOF-%d UpTdc", i);
    TString title4 = Form("TOF-%d DownTdc", i);
    TString title5 = Form("TOF-%d UpAdc(w Tdc)", i);
    TString title6 = Form("TOF-%d DownAdc(w Tdc)", i);
    TString title7 = Form("TOF-%d UpAdc(w/o Tdc)", i);
    TString title8 = Form("TOF-%d DownAdc(w/o Tdc)", i);
    HB1(TOFHid +100*i +1, title1, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*i +2, title2, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*i +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(TOFHid +100*i +4, title4, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(TOFHid +100*i +5, title5, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*i +6, title6, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*i +7, title7, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*i +8, title8, NbinAdc, MinAdc, MaxAdc);
  }

  HB1(TOFHid +10, "#Hits Tof[Hodo]",  NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +11, "Hitpat Tof[Hodo]", NumOfSegTOF,   0., Double_t(NumOfSegTOF));
  HB1(TOFHid +12, "CMeanTime Tof", 500, -5., 45.);
  HB1(TOFHid +13, "dE Tof", 200, -0.5, 4.5);
  HB1(TOFHid +14, "#Hits Tof[HodoGood]",  NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +15, "Hitpat Tof[HodoGood]", NumOfSegTOF,   0., Double_t(NumOfSegTOF));

  for(Int_t i=1; i<=NumOfSegTOF; ++i){
    TString title11 = Form("TOF-%d Up Time", i);
    TString title12 = Form("TOF-%d Down Time", i);
    TString title13 = Form("TOF-%d MeanTime", i);
    TString title14 = Form("TOF-%d Up dE", i);
    TString title15 = Form("TOF-%d Down dE", i);
    TString title16 = Form("TOF-%d dE", i);
    TString title17 = Form("TOF-%d Up CTime", i);
    TString title18 = Form("TOF-%d Down CTime", i);
    TString title19 = Form("TOF-%d CMeanTime", i);
    TString title20 = Form("TOF-%d Tup-Tdown", i);
    TString title21 = Form("TOF-%d dE (w/ TOF-HT)", i);
    HB1(TOFHid +100*i +11, title11, 500, -5., 45.);
    HB1(TOFHid +100*i +12, title12, 500, -5., 45.);
    HB1(TOFHid +100*i +13, title13, 500, -5., 45.);
    HB1(TOFHid +100*i +14, title14, 200, -0.5, 4.5);
    HB1(TOFHid +100*i +15, title15, 200, -0.5, 4.5);
    HB1(TOFHid +100*i +16, title16, 200, -0.5, 4.5);
    HB1(TOFHid +100*i +17, title17, 500, -5., 45.);
    HB1(TOFHid +100*i +18, title18, 500, -5., 45.);
    HB1(TOFHid +100*i +19, title19, 500, -5., 45.);
    HB1(TOFHid +100*i +20, title20, 200, -10.0, 10.0);
    HB1(TOFHid +100*i +21, title21, 200, -0.5, 4.5);
  }

  HB2(TOFHid +21, "TofHitPat%TofHitPat[HodoGood]", NumOfSegTOF,   0., Double_t(NumOfSegTOF),
      NumOfSegTOF,   0., Double_t(NumOfSegTOF));
  HB2(TOFHid +22, "CMeanTimeTof%CMeanTimeTof[HodoGood]",
      120, 10., 40., 120, 10., 40.);
  HB1(TOFHid +23, "TDiff Tof[HodoGood]", 200, -10., 10.);
  HB2(TOFHid +24, "TofHitPat%TofHitPat[HodoGood2]", NumOfSegTOF,   0., Double_t(NumOfSegTOF),
      NumOfSegTOF,   0., Double_t(NumOfSegTOF));

  HB1(TOFHid +30, "#Clusters Tof", NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +31, "ClusterSize Tof", 5, 0., 5.);
  HB1(TOFHid +32, "HitPat Cluster Tof", 2*NumOfSegTOF, 0., Double_t(NumOfSegTOF));
  HB1(TOFHid +33, "CMeamTime Cluster Tof", 500, -5., 45.);
  HB1(TOFHid +34, "DeltaE Cluster Tof", 100, -0.5, 4.5);

  // BVH
  HB1(BVHHid +0, "#Hits BVH",        NumOfSegBVH+1, 0., Double_t(NumOfSegBVH+1));
  HB1(BVHHid +1, "Hitpat BVH",       NumOfSegBVH,   0., Double_t(NumOfSegBVH));
  HB1(BVHHid +2, "#Hits BVH(Tor)",   NumOfSegBVH+1, 0., Double_t(NumOfSegBVH+1));
  HB1(BVHHid +3, "Hitpat BVH(Tor)",  NumOfSegBVH,   0., Double_t(NumOfSegBVH));
  HB1(BVHHid +4, "#Hits BVH(Tand)",  NumOfSegBVH+1, 0., Double_t(NumOfSegBVH+1));
  HB1(BVHHid +5, "Hitpat BVH(Tand)", NumOfSegBVH,   0., Double_t(NumOfSegBVH));

  for(Int_t i=1; i<=NumOfSegBVH; ++i){
    TString title1 = Form("BVH-%d UpAdc", i);
    TString title2 = Form("BVH-%d DownAdc", i);
    TString title3 = Form("BVH-%d UpTdc", i);
    TString title4 = Form("BVH-%d DownTdc", i);
    TString title5 = Form("BVH-%d UpAdc(w Tdc)", i);
    TString title6 = Form("BVH-%d DownAdc(w Tdc)", i);
    TString title7 = Form("BVH-%d UpAdc(w/o Tdc)", i);
    TString title8 = Form("BVH-%d DownAdc(w/o Tdc)", i);
    HB1(BVHHid +100*i +1, title1, NbinTdc, MinTdc, MaxTdc);
    HB1(BVHHid +100*i +2, title2, NbinTdc, MinTdc, MaxTdc);
    HB1(BVHHid +100*i +3, title3, NbinTdc, MinTdc, MaxTdc);
    HB1(BVHHid +100*i +4, title4, NbinTdc, MinTdc, MaxTdc);
    HB1(BVHHid +100*i +5, title5, NbinAdc, MinAdc, MaxAdc);
    HB1(BVHHid +100*i +6, title6, NbinAdc, MinAdc, MaxAdc);
    HB1(BVHHid +100*i +7, title7, NbinAdc, MinAdc, MaxAdc);
    HB1(BVHHid +100*i +8, title8, NbinAdc, MinAdc, MaxAdc);
  }

  HB1(BVHHid +10, "#Hits Bvh[Hodo]",  NumOfSegBVH+1, 0., Double_t(NumOfSegBVH+1));
  HB1(BVHHid +11, "Hitpat Bvh[Hodo]", NumOfSegBVH,   0., Double_t(NumOfSegBVH));
  HB1(BVHHid +12, "CMeanTime Bvh", 500, -5., 45.);
  HB1(BVHHid +13, "dE Bvh", 200, -0.5, 4.5);
  HB1(BVHHid +14, "#Hits Bvh[HodoGood]",  NumOfSegBVH+1, 0., Double_t(NumOfSegBVH+1));
  HB1(BVHHid +15, "Hitpat Bvh[HodoGood]", NumOfSegBVH,   0., Double_t(NumOfSegBVH));

  for(Int_t i=1; i<=NumOfSegBVH; ++i){
    TString title11 = Form("BVH-%d Up Time", i);
    TString title12 = Form("BVH-%d Down Time", i);
    TString title13 = Form("BVH-%d MeanTime", i);
    TString title14 = Form("BVH-%d Up dE", i);
    TString title15 = Form("BVH-%d Down dE", i);
    TString title16 = Form("BVH-%d dE", i);
    TString title17 = Form("BVH-%d Up CTime", i);
    TString title18 = Form("BVH-%d Down CTime", i);
    TString title19 = Form("BVH-%d CMeanTime", i);
    TString title20 = Form("BVH-%d Tup-Tdown", i);
    HB1(BVHHid +100*i +11, title11, 500, -5., 45.);
    HB1(BVHHid +100*i +12, title12, 500, -5., 45.);
    HB1(BVHHid +100*i +13, title13, 500, -5., 45.);
    HB1(BVHHid +100*i +14, title14, 200, -0.5, 4.5);
    HB1(BVHHid +100*i +15, title15, 200, -0.5, 4.5);
    HB1(BVHHid +100*i +16, title16, 200, -0.5, 4.5);
    HB1(BVHHid +100*i +17, title17, 500, -5., 45.);
    HB1(BVHHid +100*i +18, title18, 500, -5., 45.);
    HB1(BVHHid +100*i +19, title19, 500, -5., 45.);
    HB1(BVHHid +100*i +20, title20, 200, -10.0, 10.0);
  }

  HB2(BVHHid +21, "BvhHitPat%BvhHitPat[HodoGood]", NumOfSegBVH,   0., Double_t(NumOfSegBVH),
      NumOfSegBVH,   0., Double_t(NumOfSegBVH));
  HB2(BVHHid +22, "CMeanTimeBvh%CMeanTimeBvh[HodoGood]",
      120, 10., 40., 120, 10., 40.);
  HB1(BVHHid +23, "TDiff Bvh[HodoGood]", 200, -10., 10.);
  HB2(BVHHid +24, "BvhHitPat%BvhHitPat[HodoGood2]", NumOfSegBVH,   0., Double_t(NumOfSegBVH),
      NumOfSegBVH,   0., Double_t(NumOfSegBVH));

  HB1(BVHHid +30, "#Clusters Bvh", NumOfSegBVH+1, 0., Double_t(NumOfSegBVH+1));
  HB1(BVHHid +31, "ClusterSize Bvh", 5, 0., 5.);
  HB1(BVHHid +32, "HitPat Cluster Bvh", 2*NumOfSegBVH, 0., Double_t(NumOfSegBVH));
  HB1(BVHHid +33, "CMeamTime Cluster Bvh", 500, -5., 45.);
  HB1(BVHHid +34, "DeltaE Cluster Bvh", 100, -0.5, 4.5);

  // LAC
  HB1(LACHid +0, "#Hits LAC",        NumOfSegLAC+1, 0., Double_t(NumOfSegLAC+1));
  HB1(LACHid +1, "Hitpat LAC",       NumOfSegLAC,   0., Double_t(NumOfSegLAC));
  HB1(LACHid +2, "#Hits LAC(Tor)",   NumOfSegLAC+1, 0., Double_t(NumOfSegLAC+1));
  HB1(LACHid +3, "Hitpat LAC(Tor)",  NumOfSegLAC,   0., Double_t(NumOfSegLAC));
  HB1(LACHid +4, "#Hits LAC(Tand)",  NumOfSegLAC+1, 0., Double_t(NumOfSegLAC+1));
  HB1(LACHid +5, "Hitpat LAC(Tand)", NumOfSegLAC,   0., Double_t(NumOfSegLAC));
  for(Int_t i=1; i<=NumOfSegLAC; ++i){
    TString title3 = Form("LAC-%d Tdc", i);
    HB1(LACHid +100*i +3, title3, NbinTdc, MinTdc, MaxTdc);
  }
  HB1(LACHid +10, "#Hits Lac[Hodo]",  NumOfSegLAC+1, 0., Double_t(NumOfSegLAC+1));
  HB1(LACHid +11, "Hitpat Lac[Hodo]", NumOfSegLAC,   0., Double_t(NumOfSegLAC));
  HB1(LACHid +12, "CMeanTime Lac", 500, -5., 45.);
  HB1(LACHid +14, "#Hits Lac[HodoGood]",  NumOfSegLAC+1, 0., Double_t(NumOfSegLAC+1));
  HB1(LACHid +15, "Hitpat Lac[HodoGood]", NumOfSegLAC,   0., Double_t(NumOfSegLAC));
  for(Int_t i=1; i<=NumOfSegLAC; ++i){
    TString title11 = Form("LAC-%d Time", i);
    HB1(LACHid +100*i +11, title11, 500, -5., 45.);
  }

  //WC
  HB1(WCHid +0, "#Hits WC",        NumOfSegWC+1, 0., Double_t(NumOfSegWC+1));
  HB1(WCHid +1, "Hitpat WC",       NumOfSegWC,   0., Double_t(NumOfSegWC));
  HB1(WCHid +2, "#Hits WC(Tor)",   NumOfSegWC+1, 0., Double_t(NumOfSegWC+1));
  HB1(WCHid +3, "Hitpat WC(Tor)",  NumOfSegWC,   0., Double_t(NumOfSegWC));
  HB1(WCHid +4, "#Hits WC(Tand)",  NumOfSegWC+1, 0., Double_t(NumOfSegWC+1));
  HB1(WCHid +5, "Hitpat WC(Tand)", NumOfSegWC,   0., Double_t(NumOfSegWC));

  for(Int_t i=1; i<=NumOfSegWC; ++i){
    TString title1 = Form("WC-%d UpAdc", i);
    TString title2 = Form("WC-%d DownAdc", i);
    TString title3 = Form("WC-%d UpTdc", i);
    TString title4 = Form("WC-%d DownTdc", i);
    TString title5 = Form("WC-%d UpAdc(w Tdc)", i);
    TString title6 = Form("WC-%d DownAdc(w Tdc)", i);
    TString title7 = Form("WC-%d UpAdc(w/o Tdc)", i);
    TString title8 = Form("WC-%d DownAdc(w/o Tdc)", i);
    HB1(WCHid +100*i +1, title1, NbinAdc,   MinAdc,   MaxAdc);
    HB1(WCHid +100*i +2, title2, NbinAdc,   MinAdc,   MaxAdc);
    HB1(WCHid +100*i +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(WCHid +100*i +4, title4, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(WCHid +100*i +5, title5, NbinAdc,   MinAdc,   MaxAdc);
    HB1(WCHid +100*i +6, title6, NbinAdc,   MinAdc,   MaxAdc);
    HB1(WCHid +100*i +7, title7, NbinAdc,   MinAdc,   MaxAdc);
    HB1(WCHid +100*i +8, title8, NbinAdc,   MinAdc,   MaxAdc);
  }

  //WCSUM
  HB1(WCSUMHid +0, "#Hits WC SUM", NumOfSegWC+1, 0., Double_t(NumOfSegWC+1));
  HB1(WCSUMHid +1, "Hitpat WC SUM", NumOfSegWC, 0., Double_t(NumOfSegWC));
  for(Int_t i=1; i<=NumOfSegWC; ++i){
    TString title1 = Form("WCSUM-%d Adc", i);
    TString title3 = Form("WCSUM-%d Tdc", i);
    TString title5 = Form("WCSUM-%d Adc(w Tdc)", i);
    TString title7 = Form("WCSUM-%d Adc(w/o Tdc)", i);
    HB1(WCSUMHid +100*i +1, title1, NbinAdc, MinAdc, MaxAdc);
    HB1(WCSUMHid +100*i +3, title3, NbinTdc, MinTdc, MaxTdc);
    HB1(WCSUMHid +100*i +5, title5, NbinAdc, MinAdc, MaxAdc);
    HB1(WCSUMHid +100*i +7, title7, NbinAdc, MinAdc, MaxAdc);
  }
  HB1(WCSUMHid +10, "#Hits WcSum[Hodo]",  NumOfSegWC+1, 0., Double_t(NumOfSegWC+1));
  HB1(WCSUMHid +11, "Hitpat WcSum[Hodo]", NumOfSegWC,   0., Double_t(NumOfSegWC)  );
  HB1(WCSUMHid +13, "dE WcSum", 200, -0.5, 4.5);
  HB1(WCSUMHid +14, "#Hits WcSum[HodoGood]",  NumOfSegWC+1, 0., Double_t(NumOfSegWC+1));
  HB1(WCSUMHid +15, "Hitpat WcSum[HodoGood]", NumOfSegWC,   0., Double_t(NumOfSegWC)  );
  for(Int_t i=1; i<=NumOfSegWC; ++i){
    TString title11 = Form("WCSUM-%d Time", i);
    TString title14 = Form("WCSUM-%d NPE", i);
    TString title16 = Form("WCSUM-%d NPE", i);
    TString title17 = Form("WCSUM-%d Up CTime", i);
    TString title21 = Form("WCSUM-%d dE (w/ WC-HT)", i);
    HB1(WCSUMHid +100*i +11, title11, 500, -5., 45.);
    HB1(WCSUMHid +100*i +14, title14, 110, -10., 100.);
    HB1(WCSUMHid +100*i +16, title16, 110, -10., 100.);
    HB1(WCSUMHid +100*i +17, title17, 500, -5., 45.);
    HB1(WCSUMHid +100*i +21, title21, 200, -0.5, 4.5);
  }
  HB2(WCSUMHid +21, "WcSumHitPat%WcSumHitPat[HodoGood]",
      NumOfSegWC, 0., Double_t(NumOfSegWC),
      NumOfSegWC, 0., Double_t(NumOfSegWC));
  HB2(WCSUMHid +22, "CMeanTimeWcSum%CMeanTimeWcSum[HodoGood]",
      120, 10., 40., 120, 10., 40.);
  HB1(WCSUMHid +23, "TDiff WcSum[HodoGood]", 200, -10., 10.);
  HB2(WCSUMHid +24, "WcSumHitPat%WcSumHitPat[HodoGood2]",
      NumOfSegWC, 0., Double_t(NumOfSegWC),
      NumOfSegWC, 0., Double_t(NumOfSegWC));
  HB1(WCSUMHid +30, "#Clusters WcSum", NumOfSegWC+1, 0., Double_t(NumOfSegWC+1));
  HB1(WCSUMHid +31, "ClusterSize WcSum", 5, 0., 5.);
  HB1(WCSUMHid +32, "HitPat Cluster WcSum", 2*NumOfSegWC, 0., Double_t(NumOfSegWC));
  HB1(WCSUMHid +33, "CMeamTime Cluster WcSum", 500, -5., 45.);
  HB1(WCSUMHid +34, "DeltaE Cluster WcSum", 100, -0.5, 4.5);

  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("spill",     &event.spill,     "spill/I");
  //Trig
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  //BH1
  tree->Branch("bh1nhits",   &event.bh1nhits,    "bh1nhits/I");
  tree->Branch("bh1hitpat",   event.bh1hitpat,   Form("bh1hitpat[%d]/I",NumOfSegBH1));
  tree->Branch("bh1ua",       event.bh1ua,       Form("bh1ua[%d]/D", NumOfSegBH1));
  tree->Branch("bh1ut",       event.bh1ut,       Form("bh1ut[%d][%d]/D", NumOfSegBH1, MaxDepth));
  tree->Branch("bh1da",       event.bh1da,       Form("bh1da[%d]/D", NumOfSegBH1));
  tree->Branch("bh1dt",       event.bh1dt,       Form("bh1dt[%d][%d]/D", NumOfSegBH1, MaxDepth));
  //BH2
  tree->Branch("bh2nhits",   &event.bh2nhits,    "bh2nhits/I");
  tree->Branch("bh2hitpat",   event.bh2hitpat,   Form("bh2hitpat[%d]/I", NumOfSegBH2));
  tree->Branch("bh2ua",       event.bh2ua,       Form("bh2ua[%d]/D", NumOfSegBH2));
  tree->Branch("bh2ut",       event.bh2ut,       Form("bh2ut[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2da",       event.bh2da,       Form("bh2da[%d]/D", NumOfSegBH2));
  tree->Branch("bh2dt",       event.bh2dt,       Form("bh2dt[%d][%d]/D", NumOfSegBH2, MaxDepth));
  //BAC
  tree->Branch("bacnhits",   &event.bacnhits,   "bacnhits/I");
  tree->Branch("bachitpat",   event.bachitpat,  Form("bachitpat[%d]/I", NumOfSegBAC));
  tree->Branch("baca",        event.baca,       Form("baca[%d]/D", NumOfSegBAC));
  tree->Branch("bact",        event.bact,       Form("bact[%d][%d]/D", NumOfSegBAC, MaxDepth));
  //HTOF
  tree->Branch("htofnhits",   &event.htofnhits,   "htofnhits/I");
  tree->Branch("htofhitpat",   event.htofhitpat,  Form("htofhitpat[%d]/I", NumOfSegHTOF));
  tree->Branch("htofua",       event.htofua,      Form("htofua[%d]/D", NumOfSegHTOF));
  tree->Branch("htofda",       event.htofda,      Form("htofda[%d]/D", NumOfSegHTOF));
  tree->Branch("htofut",       event.htofut,      Form("htofut[%d][%d]/D", NumOfSegHTOF, MaxDepth));
  tree->Branch("htofdt",       event.htofdt,      Form("htofdt[%d][%d]/D", NumOfSegHTOF, MaxDepth));
  //BVH
  tree->Branch("bvhnhits",   &event.bvhnhits,   "bvhnhits/I");
  tree->Branch("bvhhitpat",   event.bvhhitpat,  Form("bvhhitpat[%d]/I", NumOfSegBVH));
  tree->Branch("bvht" ,       event.bvht,       Form("bvht[%d][%d]/D", NumOfSegBVH, MaxDepth));
  //TOF
  tree->Branch("tofnhits",   &event.tofnhits,   "tofnhits/I");
  tree->Branch("tofhitpat",   event.tofhitpat,  Form("tofhitpat[%d]/I", NumOfSegTOF));
  tree->Branch("tofnhits_3dmtx",   &event.tofnhits_3dmtx,   "tofnhits_3dmtx/I");
  tree->Branch("tofhitpat_3dmtx",   event.tofhitpat_3dmtx,  "tofhitpat_3dmtx[tofnhits_3dmtx]/I");
  tree->Branch("tofua",       event.tofua,      Form("tofua[%d]/D", NumOfSegTOF));
  tree->Branch("tofut",       event.tofut,      Form("tofut[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofda",       event.tofda,      Form("tofda[%d]/D", NumOfSegTOF));
  tree->Branch("tofdt",       event.tofdt,      Form("tofdt[%d][%d]/D", NumOfSegTOF, MaxDepth));
  //LAC
  tree->Branch("lacnhits", &event.lacnhits, "lacnhits/I");
  tree->Branch("lachitpat", event.lachitpat, Form("lachitpat[%d]/I", NumOfSegLAC));
  tree->Branch("lact", event.lact, Form("lact[%d][%d]/D", NumOfSegLAC, MaxDepth));

  //WC
  tree->Branch("wcnhits",   &event.wcnhits,    "wcnhits/I");
  tree->Branch("wchitpat",   event.wchitpat,   Form("wchitpat[%d]/I", NumOfSegWC));
  tree->Branch("wcua",       event.wcua,       Form("wcua[%d]/D", NumOfSegWC));
  tree->Branch("wcut",       event.wcut,       Form("wcut[%d][%d]/D", NumOfSegWC, MaxDepth));
  tree->Branch("wcda",       event.wcda,       Form("wcda[%d]/D", NumOfSegWC));
  tree->Branch("wcdt",       event.wcdt,       Form("wcdt[%d][%d]/D", NumOfSegWC, MaxDepth));
  //WCSUM
  tree->Branch("wcsumnhits", &event.wcsumnhits, "wcsumnhits/I");
  tree->Branch("wcsumhitpat", &event.wcsumhitpat, Form("wcsumhitpat[%d]/I", NumOfSegWC));
  tree->Branch("wcsuma", &event.wcsuma, Form("wcsuma[%d]/I", NumOfSegWC));
  tree->Branch("wcsumt", &event.wcsumt, Form("wcsumt[%d][%d]/D", NumOfSegWC, MaxDepth));

#if FHitBranch
  tree->Branch("sch_tdc",        event.sch_tdc,          Form("sch_tdc[%d][%d]/D",
                                                              NumOfSegSCH, MaxDepth));
  tree->Branch("sch_trailing",   event.sch_trailing,     Form("sch_trailing[%d][%d]/D",
                                                              NumOfSegSCH, MaxDepth));
  tree->Branch("sch_tot",        event.sch_tot,          Form("sch_tot[%d][%d]/D",
                                                              NumOfSegSCH, MaxDepth));
  tree->Branch("sch_depth",      event.sch_depth,        Form("sch_depth[%d]/I", NumOfSegSCH));
#endif
  tree->Branch("sch_nhits",     &event.sch_nhits,        "sch_nhits/I");
  tree->Branch("sch_hitpat",     event.sch_hitpat,       "sch_hitpat[sch_nhits]/I");

  tree->Branch("sch_ncl",       &event.sch_ncl,          "sch_ncl/I");
  tree->Branch("sch_clsize",     event.sch_clsize,       "sch_clsize[sch_ncl]/I");
  tree->Branch("sch_ctime",      event.sch_ctime,        "sch_ctime[sch_ncl]/D");
  tree->Branch("sch_ctot",       event.sch_ctot,         "sch_ctot[sch_ncl]/D");
  tree->Branch("sch_clpos",      event.sch_clpos,        "sch_clpos[sch_ncl]/D");

  //Normalized data
  tree->Branch("bh1mt",     event.bh1mt,     Form("bh1mt[%d][%d]/D", NumOfSegBH1, MaxDepth));
  tree->Branch("bh1utime",     event.bh1utime,     Form("bh1utime[%d][%d]/D", NumOfSegBH1, MaxDepth));
  tree->Branch("bh1dtime",     event.bh1dtime,     Form("bh1dtime[%d][%d]/D", NumOfSegBH1, MaxDepth));
  tree->Branch("bh1uctime",     event.bh1uctime,     Form("bh1uctime[%d][%d]/D", NumOfSegBH1, MaxDepth));
  tree->Branch("bh1dctime",     event.bh1dctime,     Form("bh1dctime[%d][%d]/D", NumOfSegBH1, MaxDepth));
  tree->Branch("bh1de",     event.bh1de,     Form("bh1de[%d]/D", NumOfSegBH1));
  tree->Branch("bh1ude",     event.bh1ude,     Form("bh1ude[%d]/D", NumOfSegBH1));
  tree->Branch("bh1dde",     event.bh1dde,     Form("bh1dde[%d]/D", NumOfSegBH1));
  tree->Branch("bh2mt",     event.bh2mt,     Form("bh2mt[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2utime",     event.bh2utime,     Form("bh2utime[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2dtime",     event.bh2dtime,     Form("bh2dtime[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2uctime",     event.bh2uctime,     Form("bh2uctime[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2dctime",     event.bh2dctime,     Form("bh2dctime[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2de",     event.bh2de,     Form("bh2de[%d]/D", NumOfSegBH2));
  tree->Branch("bh2ude",     event.bh2ude,     Form("bh2ude[%d]/D", NumOfSegBH2));
  tree->Branch("bh2dde",     event.bh2dde,     Form("bh2dde[%d]/D", NumOfSegBH2));
#if HodoHitPos
  tree->Branch("bh2hitpos",     event.bh2hitpos,     Form("bh2hitpos[%d][%d]/D", NumOfSegBH2, MaxDepth));
#endif
  tree->Branch("bacmt",     event.bacmt,     Form("bacmt[%d][%d]/D", NumOfSegBAC, MaxDepth));
  tree->Branch("bacde",     event.bacde,     Form("bacde[%d]/D", NumOfSegBAC));
  tree->Branch("htofmt",   event.htofmt,   Form("htofmt[%d][%d]/D", NumOfSegHTOF, MaxDepth));
  tree->Branch("htofutime",     event.htofutime,     Form("htofutime[%d][%d]/D", NumOfSegHTOF, MaxDepth));
  tree->Branch("htofdtime",     event.htofdtime,     Form("htofdtime[%d][%d]/D", NumOfSegHTOF, MaxDepth));
  tree->Branch("htofuctime",     event.htofuctime,     Form("htofuctime[%d][%d]/D", NumOfSegHTOF, MaxDepth));
  tree->Branch("htofdctime",     event.htofdctime,     Form("htofdctime[%d][%d]/D", NumOfSegHTOF, MaxDepth));
  tree->Branch("htofde",   event.htofde,   Form("htofde[%d]/D", NumOfSegHTOF));
  tree->Branch("htofude",   event.htofude,   Form("htofude[%d]/D", NumOfSegHTOF));
  tree->Branch("htofdde",   event.htofdde,   Form("htofdde[%d]/D", NumOfSegHTOF));
#if HodoHitPos
  tree->Branch("htofhitpos",     event.htofhitpos,     Form("htofhitpos[%d][%d]/D", NumOfSegHTOF, MaxDepth));
#endif
  tree->Branch("tofmt",     event.tofmt,     Form("tofmt[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofde",     event.tofde,     Form("tofde[%d]/D", NumOfSegTOF));
  tree->Branch("wcmt",     event.wcmt,     Form("wcmt[%d][%d]/D", NumOfSegWC, MaxDepth));
  tree->Branch("wcde",     event.wcde,     Form("wcde[%d]/D", NumOfSegWC));

  tree->Branch("t0",        event.t0,        Form("t0[%d][%d]/D",  NumOfSegBH2, MaxDepth));
  tree->Branch("ct0",       event.ct0,       Form("ct0[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("btof",      event.btof,      Form("btof[%d][%d]/D",  NumOfSegBH1, NumOfSegBH2));
  tree->Branch("cbtof",     event.cbtof,     Form("cbtof[%d][%d]/D", NumOfSegBH1, NumOfSegBH2));

  tree->Branch("Time0Seg", &event.Time0Seg,  "Time0Seg/D");
  tree->Branch("deTime0",  &event.deTime0,   "deTime0/D");
  tree->Branch("Time0",    &event.Time0,     "Time0/D");
  tree->Branch("CTime0",   &event.CTime0,    "CTime0/D");

  tree->Branch("Btof0Seg", &event.Btof0Seg,  "Btof0Seg/D");
  tree->Branch("deBtof0",  &event.deBtof0,   "deBtof0/D");
  tree->Branch("Btof0",    &event.Btof0,     "Btof0/D");
  tree->Branch("CBtof0",   &event.CBtof0,    "CBtof0/D");

  ////////////////////////////////////////////
  //Dst
  hodo = new TTree("hodo","Data Summary Table of Hodoscope");
  hodo->Branch("evnum",     &dst.evnum,     "evnum/I");
  hodo->Branch("spill",     &dst.spill,     "spill/I");
  hodo->Branch("trigpat",    dst.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  hodo->Branch("trigflag",   dst.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));
  hodo->Branch("nhBh1",     &dst.nhBh1,     "nhBh1/I");
  hodo->Branch("csBh1",      dst.csBh1,     "csBh1[nhBh1]/I");
  hodo->Branch("Bh1Seg",     dst.Bh1Seg,    "Bh1Seg[nhBh1]/D");
  hodo->Branch("tBh1",       dst.tBh1,      "tBh1[nhBh1]/D");
  hodo->Branch("dtBh1",      dst.dtBh1,     "dtBh1[nhBh1]/D");
  hodo->Branch("deBh1",      dst.deBh1,     "deBh1[nhBh1]/D");

  hodo->Branch("nhBh2",     &dst.nhBh2,     "nhBh2/I");
  hodo->Branch("csBh2",      dst.csBh2,     "csBh2[nhBh2]/I");
  hodo->Branch("Bh2Seg",     dst.Bh2Seg,    "Bh2Seg[nhBh2]/D");
  hodo->Branch("tBh2",       dst.tBh2,      "tBh2[nhBh2]/D");
  hodo->Branch("t0Bh2",      dst.t0Bh2,     "t0Bh2[nhBh2]/D");
  hodo->Branch("dtBh2",      dst.dtBh2,     "dtBh2[nhBh2]/D");
  hodo->Branch("deBh2",      dst.deBh2,     "deBh2[nhBh2]/D");
#if HodoHitPos
  hodo->Branch("posBh2",     dst.posBh2,    "posBh2[nhBh2]/D");
#endif
  hodo->Branch("btof",       dst.btof,      "btof[nhBh1]/D");
  hodo->Branch("cbtof",      dst.cbtof,     "cbtof[nhBh1]/D");

  hodo->Branch("Btof0Seg",  &dst.Btof0Seg,  "Btof0Seg/D");
  hodo->Branch("deBtof0",   &dst.deBtof0,   "deBtof0/D");
  hodo->Branch("Btof0",     &dst.Btof0,     "Btof0/D");
  hodo->Branch("CBtof0",    &dst.CBtof0,    "CBtof0/D");

  hodo->Branch("Time0Seg",  &dst.Time0Seg,  "Time0Seg/D");
  hodo->Branch("deTime0",   &dst.deTime0,   "deTime0/D");
  hodo->Branch("Time0",     &dst.Time0,     "Time0/D");
  hodo->Branch("CTime0",    &dst.CTime0,    "CTime0/D");

  hodo->Branch("nhBac",     &dst.nhBac,     "nhBac/I");
  hodo->Branch("BacSeg",     dst.BacSeg,    "BacSeg[nhBac]/D");
  hodo->Branch("tBac",       dst.tBac,      "tBac[nhBac]/D");
  hodo->Branch("deBac",      dst.deBac,     "deBac[nhBac]/D");

  hodo->Branch("nhHtof",     &dst.nhHtof,     "nhHtof/I");
  hodo->Branch("csHtof",      dst.csHtof,     "csHtof[nhHtof]/I");
  hodo->Branch("HtofSeg",     dst.HtofSeg,    "HtofSeg[nhHtof]/D");
  hodo->Branch("tHtof",       dst.tHtof,      "tHtof[nhHtof]/D");
  hodo->Branch("dtHtof",      dst.dtHtof,     "dtHtof[nhHtof]/D");
  hodo->Branch("deHtof",      dst.deHtof,     "deHtof[nhHtof]/D");
#if HodoHitPos
  hodo->Branch("posHtof",     dst.posHtof,    "posHtof[nhHtof]/D");
#endif
  hodo->Branch("nhBvh",     &dst.nhBvh,     "nhBvh/I");
  hodo->Branch("csBvh",      dst.csBvh,     "csBvh[nhBvh]/I");
  hodo->Branch("BvhSeg",     dst.BvhSeg,    "BvhSeg[nhBvh]/D");
  hodo->Branch("tBvh",       dst.tBvh,      "tBvh[nhBvh]/D");

  hodo->Branch("nhTof",     &dst.nhTof,     "nhTof/I");
  hodo->Branch("csTof",      dst.csTof,     "csTof[nhTof]/I");
  hodo->Branch("TofSeg",     dst.TofSeg,    "TofSeg[nhTof]/D");
  hodo->Branch("tTof",       dst.tTof,      "tTof[nhTof]/D");
  hodo->Branch("dtTof",      dst.dtTof,     "dtTof[nhTof]/D");
  hodo->Branch("deTof",      dst.deTof,     "deTof[nhTof]/D");

  hodo->Branch("tofua", event.tofua, Form("tofua[%d]/D", NumOfSegTOF));
  hodo->Branch("tofut", event.tofut, Form("tofut[%d][%d]/D", NumOfSegTOF, MaxDepth));
  hodo->Branch("tofda", event.tofda, Form("tofda[%d]/D", NumOfSegTOF));
  hodo->Branch("tofdt", event.tofdt, Form("tofdt[%d][%d]/D", NumOfSegTOF, MaxDepth));
  hodo->Branch("utTofSeg", dst.utTofSeg,
               Form("utTofSeg[%d][%d]/D", NumOfSegTOF, MaxDepth));
  hodo->Branch("dtTofSeg", dst.dtTofSeg,
               Form("dtTofSeg[%d][%d]/D", NumOfSegTOF, MaxDepth));
  hodo->Branch("udeTofSeg",  dst.udeTofSeg,
               Form("udeTofSeg[%d]/D", NumOfSegTOF));
  hodo->Branch("ddeTofSeg",  dst.ddeTofSeg,
               Form("ddeTofSeg[%d]/D", NumOfSegTOF));

  hodo->Branch("nhLac",     &dst.nhLac,     "nhLac/I");
  hodo->Branch("csLac",      dst.csLac,     "csLac[nhLac]/I");
  hodo->Branch("LacSeg",     dst.LacSeg,    "LacSeg[nhLac]/D");
  hodo->Branch("tLac",       dst.tLac,      "tLac[nhLac]/D");

  // HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO") &&
     InitializeParameter<HodoParamMan>("HDPRM") &&
     InitializeParameter<HodoPHCMan>("HDPHC")   &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

// -*- C++ -*-

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "RootHelper.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "KuramaLib.hh"
#include "RawData.hh"
#include "UnpackerManager.hh"

#define HodoCut    0 // with BH1/BH2
#define TimeCut    1 // in cluster analysis
#define FHitBranch 1 // make FiberHit branches (becomes heavy)

namespace
{
using namespace root;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
class UserEasiroc : public VEvent
{
private:
  RawData*      rawData;
  DCAnalyzer*   DCAna;
  HodoAnalyzer* hodoAna;

public:
  UserEasiroc();
  ~UserEasiroc();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserEasiroc::ClassName()
{
  static TString s_name("UserEasiroc");
  return s_name;
}

//_____________________________________________________________________________
UserEasiroc::UserEasiroc()
  : VEvent(),
    rawData(new RawData),
    DCAna(new DCAnalyzer),
    hodoAna(new HodoAnalyzer)
{
}

//_____________________________________________________________________________
UserEasiroc::~UserEasiroc()
{
  if(rawData) delete rawData;
  if(hodoAna) delete hodoAna;
  if(DCAna) delete DCAna;
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;
  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];
  // BFT
  Int_t    bft_nhits;
  Int_t    bft_unhits;
  Int_t    bft_dnhits;
  Int_t    bft_uhitpat[NumOfSegBFT];
  Int_t    bft_dhitpat[NumOfSegBFT];
  Double_t bft_utdc[NumOfSegBFT][MaxDepth];
  Double_t bft_dtdc[NumOfSegBFT][MaxDepth];
  Double_t bft_utrailing[NumOfSegBFT][MaxDepth];
  Double_t bft_dtrailing[NumOfSegBFT][MaxDepth];
  Double_t bft_utot[NumOfSegBFT][MaxDepth];
  Double_t bft_dtot[NumOfSegBFT][MaxDepth];
  Int_t    bft_udepth[NumOfSegBFT];
  Int_t    bft_ddepth[NumOfSegBFT];
  Int_t    bft_ncl;
  Int_t    bft_clsize[NumOfSegBFT];
  Double_t bft_ctime[NumOfSegBFT];
  Double_t bft_ctot[NumOfSegBFT];
  Double_t bft_clpos[NumOfSegBFT];
  Double_t bft_clseg[NumOfSegBFT];
  // SCH
  Int_t    sch_nhits;
  Int_t    sch_hitpat[NumOfSegSCH];
  Double_t sch_tdc[NumOfSegSCH][MaxDepth];
  Double_t sch_trailing[NumOfSegSCH][MaxDepth];
  Double_t sch_tot[NumOfSegSCH][MaxDepth];
  Int_t    sch_depth[NumOfSegSCH];
  Int_t    sch_ncl;
  Int_t    sch_clsize[NumOfSegSCH];
  Double_t sch_ctime[NumOfSegSCH];
  Double_t sch_ctot[NumOfSegSCH];
  Double_t sch_clpos[NumOfSegSCH];
  Double_t sch_clseg[NumOfSegSCH];

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum      = 0;
  bft_nhits  = 0;
  bft_unhits = 0;
  bft_dnhits = 0;
  bft_ncl    = 0;
  sch_nhits  = 0;
  sch_ncl    = 0;

  for(Int_t it=0; it<NumOfSegTrig; it++){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<NumOfSegBFT; it++){
    bft_uhitpat[it] = -1;
    bft_dhitpat[it] = -1;
    bft_udepth[it]  = 0;
    bft_ddepth[it]  = 0;
    for(Int_t that=0; that<MaxDepth; that++){
      bft_utdc[it][that] = qnan;
      bft_dtdc[it][that] = qnan;
      bft_utrailing[it][that] = qnan;
      bft_dtrailing[it][that] = qnan;
      bft_utot[it][that] = qnan;
      bft_dtot[it][that] = qnan;
    }
    bft_clsize[it] = 0;
    bft_ctime[it]  = qnan;
    bft_ctot[it]   = qnan;
    bft_clpos[it]  = qnan;
    bft_clseg[it]  = qnan;
  }

  for(Int_t it=0; it<NumOfSegSCH; it++){
    sch_hitpat[it] = -1;
    sch_depth[it]  = 0;
    for(Int_t that=0; that<MaxDepth; that++){
      sch_tdc[it][that] = qnan;
      sch_trailing[it][that] = qnan;
      sch_tot[it][that] = qnan;
    }
    sch_clsize[it] = 0;
    sch_ctime[it]  = qnan;
    sch_ctot[it]   = qnan;
    sch_clpos[it]  = qnan;
    sch_clseg[it]  = qnan;
  }
}

//_____________________________________________________________________________
namespace root
{
Event  event;
TH1   *h[MaxHist];
TTree *tree;
enum eDetHid
{
  BFTHid  = 10000,
  SCHHid  = 20000
};
}

//_____________________________________________________________________________
Bool_t
UserEasiroc::ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserEasiroc::ProcessingNormal()
{
#if HodoCut
  static const Double_t MinDeBH2   = gUser.GetParameter("DeBH2", 0);
  static const Double_t MaxDeBH2   = gUser.GetParameter("DeBH2", 1);
  static const Double_t MinDeBH1   = gUser.GetParameter("DeBH1", 0);
  static const Double_t MaxDeBH1   = gUser.GetParameter("DeBH1", 1);
  static const Double_t MinBeamToF = gUser.GetParameter("BTOF",  0);
  static const Double_t MaxBeamToF = gUser.GetParameter("BTOF",  1);
#endif
  static const Double_t MinTdcBFT  = gUser.GetParameter("TdcBFT",  0);
  static const Double_t MaxTdcBFT  = gUser.GetParameter("TdcBFT",  1);
  static const Double_t MinTdcSCH  = gUser.GetParameter("TdcSCH",  0);
  static const Double_t MaxTdcSCH  = gUser.GetParameter("TdcSCH",  1);
#if TimeCut
  static const Double_t MinTimeBFT = gUser.GetParameter("TimeBFT", 0);
  static const Double_t MaxTimeBFT = gUser.GetParameter("TimeBFT", 1);
  static const Double_t MinTimeSCH = gUser.GetParameter("TimeSCH", 0);
  static const Double_t MaxTimeSCH = gUser.GetParameter("TimeSCH", 1);
#endif

  rawData->DecodeHits();

  event.evnum = gUnpacker.get_event_number();

  HF1(1, 0);

  ///// Trigger Flag
  std::bitset<NumOfSegTrig> trigger_flag;
  {
    for(const auto& hit: rawData->GetTrigRawHC()){
      Int_t seg = hit->SegmentId();
      Int_t tdc = hit->GetTdc1();
      if(tdc > 0){
	event.trigpat[trigger_flag.count()] = seg;
	event.trigflag[seg] = tdc;
        trigger_flag.set(seg);
	HF1(10, seg);
	HF1(10+seg, tdc);
      }
    }
  }

  if(trigger_flag[trigger::kSpillEnd]) return true;

  HF1(1, 1);

  ////////// BH2 time 0
  hodoAna->DecodeBH2Hits(rawData);
  Int_t nhBh2 = hodoAna->GetNHitsBH2();
#if HodoCut
  if (nhBh2==0) return true;
#endif
  HF1(1, 2);
  Double_t time0 = qnan;
  ////////// BH2 Analysis
  for(Int_t i=0; i<nhBh2; ++i){
    BH2Hit *hit = hodoAna->GetHitBH2(i);
    if(!hit) continue;
    Double_t cmt = hit->CMeanTime();
    Double_t ct0 = hit->CTime0();
    Double_t min_time = qnan;
#if HodoCut
    Double_t dE  = hit->DeltaE();
    if(dE<MinDeBH2 || MaxDeBH2<dE)
      continue;
#endif
    if(std::abs(cmt)<std::abs(min_time)){
      min_time = cmt;
      time0    = ct0;
    }
  }

  HF1(1, 3);

  ////////// BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData);
  Int_t nhBh1 = hodoAna->GetNHitsBH1();
#if HodoCut
  if(nhBh1==0) return true;
#endif
  HF1(1, 4);
  Double_t btof0 = qnan;
  for(Int_t i=0; i<nhBh1; ++i){
    Hodo2Hit* hit = hodoAna->GetHitBH1(i);
    if(!hit) continue;
    Double_t cmt  = hit->CMeanTime();
    Double_t btof = cmt - time0;
#if HodoCut
    Double_t dE   = hit->DeltaE();
    if(dE<MinDeBH1 || MaxDeBH1<dE) continue;
    if(btof<MinBeamToF || MaxBeamToF<btof) continue;
#endif
    if(std::abs(btof)<std::abs(btof0)){
      btof0 = btof;
    }
  }

  HF1(1, 5);

  HF1(1, 6);

  ////////// BFT
  {
    hodoAna->DecodeBFTHits(rawData);
    // Fiber Hit
    Int_t unhits = 0;
    Int_t dnhits = 0;
    for(Int_t plane = 0; plane<NumOfPlaneBFT; ++plane){
      Int_t nh = hodoAna->GetNHitsBFT(plane);
      enum { U, D };
      for(Int_t i=0; i<nh; ++i){
	const FiberHit* hit = hodoAna->GetHitBFT(plane, i);
	if(!hit) continue;
	Int_t mhit_l  = hit->GetNLeading();
	Int_t mhit_t  = hit->GetNTrailing();
	Int_t seg   = hit->SegmentId();
	if(plane==U) event.bft_udepth[seg] = mhit_l;
	if(plane==D) event.bft_ddepth[seg] = mhit_l;

	Int_t  prev = 0;
	Bool_t hit_flag = false;

	// raw leading data
	for(Int_t m=0; m<mhit_l; ++m){
	  if(mhit_l > MaxDepth) break;
	  Double_t leading  = hit->GetLeading(m);

	  if(leading==prev) continue;
	  prev = leading;
	  HF1(BFTHid +plane + 6, leading);
	  HF2(BFTHid +plane + 10, seg, leading);
	  HF1(BFTHid +1000*(1+plane)+seg+1, leading);

	  if(plane==U){
	    event.bft_utdc[seg][m]      = leading;
	  }
	  if(plane==D){
	    event.bft_dtdc[seg][m]      = leading;
	  }

	  if(MinTdcBFT<leading && leading<MaxTdcBFT){
	    hit_flag = true;
	  }
	}// for(m)

	// raw leading data
	for(Int_t m=0; m<mhit_t; ++m){
	  if(mhit_t > MaxDepth) break;
	  Double_t trailing = hit->GetTrailing(m);

	  if(plane==U){
	    event.bft_utrailing[seg][m] = trailing;
	  }
	  if(plane==D){
	    event.bft_dtrailing[seg][m] = trailing;
	  }

	}// for(m)

	Int_t mhit_pair  = hit->GetNPair();
	// pair data
	for(Int_t m=0; m<mhit_pair; ++m){
	  if(mhit_pair > MaxDepth) break;
	  Double_t time     = hit->GetTime(m);
	  Double_t ctime    = hit->GetCTime(m);
	  Double_t width    = hit->GetWidth(m);

	  HF1(BFTHid +plane+8, width);
	  HF2(BFTHid +plane+12, seg, width);
	  HF1(BFTHid +plane+21, time);
	  HF1(BFTHid +plane+31, ctime);
	  HF1(BFTHid +1000*(plane+3)+seg+1, width);
	  if(-10.<time && time<10.){
            HF2(BFTHid +plane+23, width, time);
            HF2(BFTHid +plane+33, width, ctime);
	    HF2(BFTHid +1000*(plane+5)+seg+1, width, time);
	    HF2(BFTHid +1000*(plane+7)+seg+1, width, ctime);
	  }
	  if(plane==U){
	    event.bft_utot[seg][m]      = width;
	  }
	  if(plane==D){
	    event.bft_dtot[seg][m]      = width;
	  }

	}
	if(hit_flag){
	  HF1(BFTHid +plane+4, seg+0.5);
	  if(plane==U) event.bft_uhitpat[unhits++] = seg;
	  if(plane==D) event.bft_dhitpat[dnhits++] = seg;
	}
      }
    }
    HF1(BFTHid +1, unhits);
    HF1(BFTHid +2, dnhits);
    HF1(BFTHid +3, unhits + dnhits);
    event.bft_unhits = unhits;
    event.bft_dnhits = dnhits;
    event.bft_nhits  = unhits + dnhits;

    // Fiber Cluster
#if TimeCut
    hodoAna->TimeCutBFT(MinTimeBFT, MaxTimeBFT);
#endif
    Int_t ncl = hodoAna->GetNClustersBFT();
    if(ncl > NumOfSegBFT){
      // std::cout << "#W BFT too much number of clusters" << std::endl;
      ncl = NumOfSegBFT;
    }
    event.bft_ncl = ncl;
    HF1(BFTHid +101, ncl);
    for(Int_t i=0; i<ncl; ++i){
      FiberCluster *cl = hodoAna->GetClusterBFT(i);
      if(!cl) continue;
      Double_t clsize = cl->ClusterSize();
      Double_t ctime  = cl->CMeanTime();
      Double_t ctot   = cl->Width();
      Double_t pos    = cl->MeanPosition();
      Double_t seg    = cl->MeanSeg();
      event.bft_clsize[i] = clsize;
      event.bft_ctime[i]  = ctime;
      event.bft_ctot[i]   = ctot;
      event.bft_clpos[i]  = pos;
      event.bft_clseg[i]  = seg;
      HF1(BFTHid +102, clsize);
      HF1(BFTHid +103, ctime);
      HF1(BFTHid +104, ctot);
      HF2(BFTHid +105, ctot, ctime);
      HF1(BFTHid +106, pos);
    }
  }

  ////////// SCH
  {
    hodoAna->DecodeSCHHits(rawData);
    Int_t nh = hodoAna->GetNHitsSCH();
    Int_t sch_nhits = 0;
    for(Int_t i=0; i<nh; ++i){
      FiberHit* hit = hodoAna->GetHitSCH(i);
      if(!hit) continue;
      Int_t mhit_l = hit->GetNLeading();
      Int_t mhit_t = hit->GetNTrailing();
      Int_t seg    = hit->SegmentId();
      event.sch_depth[seg] = mhit_l;
      Int_t  prev     = 0;
      Bool_t hit_flag = false;

      for(Int_t m=0; m<mhit_l; ++m){
	if(mhit_l > MaxDepth) break;
	Double_t leading  = hit->GetLeading(m);
	if(leading==prev) continue;
	prev = leading;
	HF1(SCHHid +3, leading);
	HF2(SCHHid +5, seg, leading);
	HF1(SCHHid +1000+seg+1, leading);

	event.sch_tdc[seg][m]      = leading;

	if(MinTdcSCH<leading && leading<MaxTdcSCH){
	  hit_flag = true;
	}
      }// for(m)

      for(Int_t m=0; m<mhit_t; ++m){
	if(mhit_t > MaxDepth) break;
	Double_t trailing = hit->GetTrailing(m);
	event.sch_trailing[seg][m] = trailing;
      }// for(m)

      Int_t mhit_pair = hit->GetNPair();
      for(Int_t m=0; m<mhit_pair; ++m){
	if(mhit_pair > MaxDepth) break;

	Double_t time     = hit->GetTime(m);
	Double_t ctime    = hit->GetCTime(m);
	Double_t width    = hit->GetWidth(m);

	HF1(SCHHid +4, width);
	HF2(SCHHid +6, seg, width);
	HF1(SCHHid +21, time);
	HF2(SCHHid +22, width, time);
	HF1(SCHHid +31, ctime);
	HF2(SCHHid +32, width, ctime);
	HF1(SCHHid +2000+seg+1, width);
	if(-10.<time && time<10.){
	  HF2(SCHHid +3000+seg+1, width, time);
	  HF2(SCHHid +4000+seg+1, width, ctime);
	}

	event.sch_tot[seg][m]      = width;

      }//for(m)
      if(hit_flag){
	HF1(SCHHid +2, seg+0.5);
	event.sch_hitpat[sch_nhits++] = seg;
      }
    }
    HF1(SCHHid +1, sch_nhits);
    event.sch_nhits = sch_nhits;

    // Fiber Cluster
#if TimeCut
    hodoAna->TimeCutSCH(MinTimeSCH, MaxTimeSCH);
#endif
    Int_t ncl = hodoAna->GetNClustersSCH();
    if(ncl > NumOfSegSCH){
      // std::cout << "#W SCH too much number of clusters" << std::endl;
      ncl = NumOfSegSCH;
    }
    event.sch_ncl = ncl;
    HF1(SCHHid +101, ncl);
    for(Int_t i=0; i<ncl; ++i){
      FiberCluster *cl = hodoAna->GetClusterSCH(i);
      if(!cl) continue;
      Double_t clsize = cl->ClusterSize();
      Double_t ctime  = cl->CMeanTime();
      Double_t ctot   = cl->Width();
      Double_t pos    = cl->MeanPosition();
      Double_t seg    = cl->MeanSeg();
      event.sch_clsize[i] = clsize;
      event.sch_ctime[i]  = ctime;
      event.sch_ctot[i]   = ctot;
      event.sch_clpos[i]  = pos;
      event.sch_clseg[i]  = seg;
      HF1(SCHHid +102, clsize);
      HF1(SCHHid +103, ctime);
      HF1(SCHHid +104, ctot);
      HF2(SCHHid +105, ctot, ctime);
      HF1(SCHHid +106, pos);
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
UserEasiroc::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserEasiroc;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  const Int_t    NbinTdc = 1000;
  const Double_t MinTdc  =    0.;
  const Double_t MaxTdc  = 1000.;

  const Int_t    NbinTot =  170;
  const Double_t MinTot  =  -10.;
  const Double_t MaxTot  =  160.;

  const Int_t    NbinTime = 100;
  const Double_t MinTime  = -50.;
  const Double_t MaxTime  =  50.;

  HB1( 1, "Status",  20,   0., 20.);
  HB1(10, "Trigger HitPat", NumOfSegTrig, 0., Double_t(NumOfSegTrig));
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1(10+i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000);
  }

  //BFT
  HB1(BFTHid + 1, "BFT Nhits U",   NumOfSegBFT, 0., (Double_t)NumOfSegBFT);
  HB1(BFTHid + 2, "BFT Nhits D",   NumOfSegBFT, 0., (Double_t)NumOfSegBFT);
  HB1(BFTHid + 3, "BFT Nhits",     NumOfSegBFT, 0., (Double_t)NumOfSegBFT);
  HB1(BFTHid + 4, "BFT Hitpat U",  NumOfSegBFT, 0., (Double_t)NumOfSegBFT);
  HB1(BFTHid + 5, "BFT Hitpat D",  NumOfSegBFT, 0., (Double_t)NumOfSegBFT);
  HB1(BFTHid + 6, "BFT Tdc U",      NbinTdc, MinTdc, MaxTdc);
  HB1(BFTHid + 7, "BFT Tdc D",      NbinTdc, MinTdc, MaxTdc);
  HB1(BFTHid + 8, "BFT Tot U",      NbinTot, MinTot, MaxTot);
  HB1(BFTHid + 9, "BFT Tot D",      NbinTot, MinTot, MaxTot);
  HB2(BFTHid +10, "BFT Tdc U%Seg",
      NumOfSegBFT, 0., (Double_t)NumOfSegBFT, NbinTdc, MinTdc, MaxTdc);
  HB2(BFTHid +11, "BFT Tdc D%Seg",
      NumOfSegBFT, 0., (Double_t)NumOfSegBFT, NbinTdc, MinTdc, MaxTdc);
  HB2(BFTHid +12, "BFT Tot U%Seg",
      NumOfSegBFT, 0., (Double_t)NumOfSegBFT, NbinTot, MinTot, MaxTot);
  HB2(BFTHid +13, "BFT Tot D%Seg",
      NumOfSegBFT, 0., (Double_t)NumOfSegBFT, NbinTot, MinTot, MaxTot);
  for(Int_t i=0; i<NumOfSegBFT; i++){
    HB1(BFTHid +1000+i+1, Form("BFT Tdc U-%d", i+1), NbinTdc, MinTdc, MaxTdc);
    HB1(BFTHid +2000+i+1, Form("BFT Tdc D-%d", i+1), NbinTdc, MinTdc, MaxTdc);
    HB1(BFTHid +3000+i+1, Form("BFT Tot U-%d", i+1), NbinTot, MinTot, MaxTot);
    HB1(BFTHid +4000+i+1, Form("BFT Tot D-%d", i+1), NbinTot, MinTot, MaxTot);
    HB2(BFTHid +5000+i+1, Form("BFT Time%%Tot U-%d", i+1),
        NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
    HB2(BFTHid +6000+i+1, Form("BFT Time%%Tot D-%d", i+1),
        NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
    HB2(BFTHid +7000+i+1, Form("BFT CTime%%Tot U-%d", i+1),
        NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
    HB2(BFTHid +8000+i+1, Form("BFT CTime%%Tot D-%d", i+1),
        NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  }
  HB1(BFTHid +21, "BFT Time U",     NbinTime, MinTime, MaxTime);
  HB1(BFTHid +22, "BFT Time D",     NbinTime, MinTime, MaxTime);
  HB2(BFTHid +23, "BFT Time%Tot U", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB2(BFTHid +24, "BFT Time%Tot D", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB1(BFTHid +31, "BFT CTime U",     NbinTime, MinTime, MaxTime);
  HB1(BFTHid +32, "BFT CTime D",     NbinTime, MinTime, MaxTime);
  HB2(BFTHid +33, "BFT CTime%Tot U", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB2(BFTHid +34, "BFT CTime%Tot D", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);

  HB1(BFTHid +101, "BFT NCluster", 100, 0, 100);
  HB1(BFTHid +102, "BFT Cluster Size", 5, 0, 5);
  HB1(BFTHid +103, "BFT CTime (Cluster)", 100., -20., 30.);
  HB1(BFTHid +104, "BFT Tot (Cluster)", NbinTot, MinTot, MaxTot);
  HB2(BFTHid +105, "BFT CTime%Tot (Cluster)",
      NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB1(BFTHid +106, "BFT Cluster Position",
      NumOfSegBFT, -0.5*(Double_t)NumOfSegBFT, 0.5*(Double_t)NumOfSegBFT);

  //SCH
  HB1(SCHHid +1, "SCH Nhits",     NumOfSegSCH, 0., (Double_t)NumOfSegSCH);
  HB1(SCHHid +2, "SCH Hitpat",    NumOfSegSCH, 0., (Double_t)NumOfSegSCH);
  HB1(SCHHid +3, "SCH Tdc",      NbinTdc, MinTdc, MaxTdc);
  HB1(SCHHid +4, "SCH Tot",      NbinTot, MinTot, MaxTot);
  HB2(SCHHid +5, "SCH Tdc U%Seg",
      NumOfSegSCH, 0., (Double_t)NumOfSegSCH, NbinTdc, MinTdc, MaxTdc);
  HB2(SCHHid +6, "SCH Tot U%Seg",
      NumOfSegSCH, 0., (Double_t)NumOfSegSCH, NbinTot, MinTot, MaxTot);
  HB2(SCHHid +13, "SCH Tot D%Seg",
      NumOfSegSCH, 0., (Double_t)NumOfSegSCH, NbinTot, MinTot, MaxTot);
  for(Int_t i=0; i<NumOfSegSCH; i++){
    HB1(SCHHid +1000+i+1, Form("SCH Tdc %d", i+1), NbinTdc, MinTdc, MaxTdc);
    HB1(SCHHid +2000+i+1, Form("SCH Tot %d", i+1), NbinTot, MinTot, MaxTot);
    HB2(SCHHid +3000+i+1, Form("SCH Time/Tot %d", i+1),
        NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
    HB2(SCHHid +4000+i+1, Form("SCH CTime/Tot %d", i+1),
        NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  }
  HB1(SCHHid +21, "SCH Time",     NbinTime, MinTime, MaxTime);
  HB2(SCHHid +22, "SCH Time/Tot", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB1(SCHHid +31, "SCH CTime",     NbinTime, MinTime, MaxTime);
  HB2(SCHHid +32, "SCH CTime/Tot", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);

  HB1(SCHHid +101, "SCH NCluster", 20, 0, 20);
  HB1(SCHHid +102, "SCH Cluster Size", 5, 0, 5);
  HB1(SCHHid +103, "SCH CTime (Cluster)", NbinTime, MinTime, MaxTime);
  HB1(SCHHid +104, "SCH Tot (Cluster)", NbinTot, MinTot, MaxTot);
  HB2(SCHHid +105, "SCH CTime%Tot (Cluster)",
      NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB1(SCHHid +106, "SCH Cluster Position",
      NumOfSegSCH, -0.5*(Double_t)NumOfSegSCH, 0.5*(Double_t)NumOfSegSCH);

  //Tree
  HBTree("ea0c", "tree of Easiroc");
  //Trig
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  //BFT
#if FHitBranch
  tree->Branch("bft_nhits",     &event.bft_nhits,        "bft_nhits/I");
  tree->Branch("bft_unhits",    &event.bft_unhits,       "bft_unhits/I");
  tree->Branch("bft_dnhits",    &event.bft_dnhits,       "bft_dnhits/I");
  tree->Branch("bft_uhitpat",    event.bft_uhitpat,      "bft_uhitpat[bft_unhits]/I");
  tree->Branch("bft_dhitpat",    event.bft_dhitpat,      "bft_dhitpat[bft_dnhits]/I");
  tree->Branch("bft_utdc",       event.bft_utdc,         Form("bft_utdc[%d][%d]/D",
							      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_dtdc",       event.bft_dtdc,         Form("bft_dtdc[%d][%d]/D",
							      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_utrailing",  event.bft_utrailing,    Form("bft_utrailing[%d][%d]/D",
							      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_dtrailing",  event.bft_dtrailing,    Form("bft_dtrailing[%d][%d]/D",
							      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_utot",       event.bft_utot,         Form("bft_utot[%d][%d]/D",
                                                              NumOfSegBFT, MaxDepth));
  tree->Branch("bft_dtot",       event.bft_dtot,         Form("bft_dtot[%d][%d]/D",
                                                              NumOfSegBFT, MaxDepth));
  tree->Branch("bft_udepth",     event.bft_udepth,       Form("bft_udepth[%d]/I", NumOfSegBFT));
  tree->Branch("bft_ddepth",     event.bft_ddepth,       Form("bft_ddepth[%d]/I", NumOfSegBFT));
#endif
  tree->Branch("bft_ncl",       &event.bft_ncl,          "bft_ncl/I");
  tree->Branch("bft_clsize",     event.bft_clsize,       "bft_clsize[bft_ncl]/I");
  tree->Branch("bft_ctime",      event.bft_ctime,        "bft_ctime[bft_ncl]/D");
  tree->Branch("bft_ctot",       event.bft_ctot,         "bft_ctot[bft_ncl]/D");
  tree->Branch("bft_clpos",      event.bft_clpos,        "bft_clpos[bft_ncl]/D");
  tree->Branch("bft_clseg",      event.bft_clseg,        "bft_clseg[bft_ncl]/D");

  //SCH
#if FHitBranch
  tree->Branch("sch_nhits",     &event.sch_nhits,        "sch_nhits/I");
  tree->Branch("sch_hitpat",     event.sch_hitpat,       "sch_hitpat[sch_nhits]/I");
  tree->Branch("sch_tdc",        event.sch_tdc,          Form("sch_tdc[%d][%d]/D",
							      NumOfSegSCH, MaxDepth));
  tree->Branch("sch_trailing",   event.sch_trailing,     Form("sch_trailing[%d][%d]/D",
							      NumOfSegSCH, MaxDepth));
  tree->Branch("sch_tot",        event.sch_tot,          Form("sch_tot[%d][%d]/D",
							      NumOfSegSCH, MaxDepth));
  tree->Branch("sch_depth",      event.sch_depth,        Form("sch_depth[%d]/I", NumOfSegSCH));
#endif
  tree->Branch("sch_ncl",       &event.sch_ncl,          "sch_ncl/I");
  tree->Branch("sch_clsize",     event.sch_clsize,       "sch_clsize[sch_ncl]/I");
  tree->Branch("sch_ctime",      event.sch_ctime,        "sch_ctime[sch_ncl]/D");
  tree->Branch("sch_ctot",       event.sch_ctot,         "sch_ctot[sch_ncl]/D");
  tree->Branch("sch_clpos",      event.sch_clpos,        "sch_clpos[sch_ncl]/D");
  tree->Branch("sch_clseg",      event.sch_clseg,        "sch_clseg[sch_ncl]/D");

  HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")    &&
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

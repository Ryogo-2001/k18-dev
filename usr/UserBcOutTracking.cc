// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

#include "ConfMan.hh"
#include "DCRawHit.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "KuramaLib.hh"
#include "RawData.hh"
#include "UnpackerManager.hh"

#define HodoCut 0
#define TotCut  1
#define Chi2Cut 1
#define ExclusiveTracking 0

namespace
{
using namespace root;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
class UserBcOutTracking : public VEvent
{
private:
  RawData*      rawData;
  DCAnalyzer*   DCAna;
  HodoAnalyzer* hodoAna;

public:
  UserBcOutTracking();
  ~UserBcOutTracking();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserBcOutTracking::ClassName()
{
  static TString s_name("UserBcOutTracking");
  return s_name;
}

//_____________________________________________________________________________
UserBcOutTracking::UserBcOutTracking()
  : VEvent(),
    rawData(new RawData),
    DCAna(new DCAnalyzer),
    hodoAna(new HodoAnalyzer)
{
}

//_____________________________________________________________________________
UserBcOutTracking::~UserBcOutTracking()
{
  if(rawData) delete rawData;
  if(DCAna) delete DCAna;
  if(hodoAna) delete hodoAna;
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;
  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  Int_t nhBh1;
  Double_t tBh1[MaxHits];
  Double_t deBh1[MaxHits];

  Int_t nhBh2;
  Double_t tBh2[MaxHits];
  Double_t deBh2[MaxHits];
  Double_t Bh2Seg[MaxHits];

  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;

  Int_t nhBac;

  Double_t btof;

  Int_t nhit[NumOfLayersBcOut];
  Int_t nlayer;
  Double_t pos[NumOfLayersBcOut][MaxHits];
  int tdc1st[NumOfLayersBcOut][MaxHits];
  int tot1st[NumOfLayersBcOut][MaxHits];
  int wire[NumOfLayersBcOut][MaxHits];

  Int_t ntrack;
  Double_t chisqr[MaxHits];
  Double_t x_Bh2[MaxHits];
  Double_t x0[MaxHits];
  Double_t y0[MaxHits];
  Double_t u0[MaxHits];
  Double_t v0[MaxHits];
  Double_t theta[MaxHits];

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum     =  0;
  nlayer    =  0;
  ntrack    =  0;
  nhBh2     =  0;
  nhBh1     =  0;
  nhBac     =  0;
  btof      = qnan;
  Time0Seg  = -1;
  deTime0   = -1;
  Time0     = qnan;
  CTime0    = qnan;

  for(Int_t it=0; it<MaxHits; it++){
    tBh1[it]   = qnan;
    deBh1[it]  = qnan;
    Bh2Seg[it] = -1;
    tBh2[it]   = qnan;
    deBh2[it]  = qnan;
    chisqr[it] = qnan;
    x0[it]     = qnan;
    y0[it]     = qnan;
    u0[it]     = qnan;
    v0[it]     = qnan;
    theta[it]  = qnan;
  }

  for(Int_t it=0; it<NumOfSegTrig; it++){
    trigpat[it] = -1;
    trigflag[it] = -1;
  }

  for (Int_t it=0; it<NumOfLayersBcOut; it++){
    nhit[it] = -1;
    for(Int_t that=0; that<MaxHits; that++){
      pos[it][that] = qnan;
      tdc1st[it][that] = qnan;
      tot1st[it][that] = qnan;
      wire[it][that] = qnan;
    }
  }
}

//_____________________________________________________________________________
namespace root
{
Event event;
TH1   *h[MaxHist];
TTree *tree;
enum eParticle { kKaon, kPion, nParticle };
}

//_____________________________________________________________________________
Bool_t
UserBcOutTracking::ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserBcOutTracking::ProcessingNormal()
{
#if HodoCut
  static const auto MinDeBH2 = gUser.GetParameter("DeBH2", 0);
  static const auto MaxDeBH2 = gUser.GetParameter("DeBH2", 1);
  static const auto MinDeBH1 = gUser.GetParameter("DeBH1", 0);
  static const auto MaxDeBH1 = gUser.GetParameter("DeBH1", 1);
  static const auto MinBeamToF = gUser.GetParameter("BTOF", 0);
  static const auto MaxBeamToF = gUser.GetParameter("BTOF", 1);
#endif
  static const auto MaxMultiHitBcOut = gUser.GetParameter("MaxMultiHitBcOut");

  static const auto MinDriftTimeBC34 = gUser.GetParameter("DriftTimeBC34", 0);
  static const auto MaxDriftTimeBC34 = gUser.GetParameter("DriftTimeBC34", 1);

#if TotCut
  static const auto MinTotBcOut = gUser.GetParameter("MinTotBcOut");
#endif
#if Chi2Cut
  static const auto MaxChisqrBcOut = gUser.GetParameter("MaxChisqrBcOut");
#endif

  rawData->DecodeHits();

  event.evnum = gUnpacker.get_event_number();

  // Trigger Flag
  std::bitset<NumOfSegTrig> trigger_flag;
  for(const auto& hit: rawData->GetTrigRawHC()){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc1();
    if(tdc > 0){
      event.trigpat[trigger_flag.count()] = seg;
      event.trigflag[seg] = tdc;
      trigger_flag.set(seg);
    }
  }

  HF1(1, 0.);

  if(trigger_flag[trigger::kSpillEnd]) return true;

  HF1(1, 1.);

  //////////////BH2 time 0
  hodoAna->DecodeBH2Hits(rawData);
  Int_t nhBh2 = hodoAna->GetNHitsBH2();
  event.nhBh2 = nhBh2;
#if HodoCut
  if(nhBh2==0) return true;
#endif
  HF1(1, 2);

  Double_t time0 = qnan;
  //////////////BH2 Analysis
  for(Int_t i=0; i<nhBh2; ++i){
    BH2Hit* hit = hodoAna->GetHitBH2(i);
    if(!hit) continue;
    Double_t seg = hit->SegmentId()+1;
    Double_t cmt = hit->CMeanTime();
    Double_t dE  = hit->DeltaE();

#if HodoCut
    if(dE<MinDeBH2 || MaxDeBH2<dE) continue;
#endif
    event.tBh2[i]   = cmt;
    event.deBh2[i]  = dE;
    event.Bh2Seg[i] = seg;
  }

  BH2Cluster *cl_time0 = hodoAna->GetTime0BH2Cluster();
  if(cl_time0){
    event.Time0Seg = cl_time0->MeanSeg()+1;
    event.deTime0  = cl_time0->DeltaE();
    event.Time0    = cl_time0->Time0();
    event.CTime0   = cl_time0->CTime0();
    time0          = cl_time0->CTime0();
  } else {
#if HodoCut
    return true;
#endif
  }

  HF1(1, 3.);

  //////////////BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData);
  Int_t nhBh1 = hodoAna->GetNHitsBH1();
  event.nhBh1 = nhBh1;
#if HodoCut
  if(nhBh1==0) return true;
#endif
  HF1(1, 4);

  for(Int_t i=0; i<nhBh1; ++i){
    Hodo2Hit *hit = hodoAna->GetHitBH1(i);
    if(!hit) continue;
    Double_t cmt = hit->CMeanTime();
    Double_t dE  = hit->DeltaE();
#if HodoCut
    if(dE<MinDeBH1 || MaxDeBH1<dE) continue;
    if(btof<MinBeamToF || MaxBeamToF<btof) continue;
#endif
    event.tBh1[i]  = cmt;
    event.deBh1[i] = dE;
  }

  Double_t btof0 = qnan;
  HodoCluster* cl_btof0 = event.Time0Seg > 0 ?
    hodoAna->GetBtof0BH1Cluster(event.CTime0) : nullptr;
  if(cl_btof0) btof0 = cl_btof0->CMeanTime() - time0;
  event.btof = btof0;

  HF1(1, 5.);

  // BAC
  hodoAna->DecodeBACHits(rawData);
  Int_t nhbac = hodoAna->GetNHitsBAC();
  event.nhBac = nhbac;

  //////////////BC3&4 number of hit layer
  DCAna->DecodeRawHits(rawData);
#if TotCut
  DCAna->TotCutBCOut(MinTotBcOut);
#endif
  DCAna->DriftTimeCutBC34(MinDriftTimeBC34, MaxDriftTimeBC34);
  // DCAna->MakeBH2DCHit(event.Time0Seg-1);

  Double_t multi_BcOut = 0.;
  Double_t multiplicity[] = {0, 0};
  struct hit_info
  {
    Double_t dt;
    Double_t pos;
  };
  {
    for(Int_t layer=1; layer<=NumOfLayersBcOut; ++layer){
      const DCHitContainer &contOut = DCAna->GetBcOutHC(layer);
      Int_t nhOut = contOut.size();
      event.nhit[layer-1] = nhOut;
      if(nhOut>0) event.nlayer++;
      if(layer ==  1) multiplicity[0] = nhOut;
      if(layer == 12) multiplicity[1] = nhOut;
      multi_BcOut += nhOut;
      HF1(100*layer, nhOut);
      Int_t plane_eff = (layer-1)*3;
      Bool_t fl_valid_sig = false;
      std::vector<hit_info> hit_cont;

      for(Int_t i=0; i<nhOut; ++i){
	DCHit *hit = contOut[i];
	Double_t wire = hit->GetWire();
	HF1(100*layer+1, wire-0.5);
	Int_t nhtdc = hit->GetTdcSize();
	Int_t tdc1st = -1;
	for(Int_t k=0; k<nhtdc; k++){
	  Int_t tdc = hit->GetTdcVal(k);
	  HF1(100*layer+2, tdc);
	  HF1(10000*layer+int(wire), tdc);
	  //	  HF2(1000*layer, tdc, wire-0.5);
	  if(tdc > tdc1st){
	    tdc1st = tdc;
	    fl_valid_sig = true;
	  }
	}
	event.tdc1st[layer-1][i] = tdc1st;
	event.wire[layer-1][i] = (int)wire;

	HF1(100*layer+6, tdc1st);
	for(Int_t k=0, n=hit->GetTdcTrailingSize(); k<n; ++k){
	  Int_t trailing = hit->GetTdcTrailing(k);
	  HF1(100*layer+10, trailing);
	}

	if(i<MaxHits)
	  event.pos[layer-1][i] = hit->GetWirePosition();

	Int_t nhdt = hit->GetDriftTimeSize();
	int tot1st = -1;
	for(Int_t k=0; k<nhdt; k++){
	  Double_t dt = hit->GetDriftTime(k);
	  HF1(100*layer+3, dt);
	  HF1(10000*layer+1000+(Int_t)(wire), dt);

	  Double_t tot = hit->GetTot(k);
	  HF1(100*layer+5, tot);
	  if(tot > tot1st){
	    tot1st = tot;
	  }

	  hit_info one_hit;
	  one_hit.dt  = dt;
	  one_hit.pos = wire;
	  hit_cont.push_back(one_hit);
	}
	HF1(100*layer+7, tot1st);
	event.tot1st[layer-1][i] = tot1st;
	Int_t nhdl = hit->GetDriftTimeSize();
	for(Int_t k=0; k<nhdl; k++){
	  Double_t dl = hit->GetDriftLength(k);
	  HF1(100*layer+4, dl);
	}
      }
      if(fl_valid_sig){
	++plane_eff;
      } else {
      }
      HF1(38, plane_eff);

      std::sort(hit_cont.begin(), hit_cont.end(),
                [](const hit_info& a_info, const hit_info& b_info)->Bool_t
                  { return a_info.dt < b_info.dt; });

      for(Int_t i=1; i<hit_cont.size(); ++i){
	HF1(100*layer+8, hit_cont.at(i).dt -hit_cont.at(0).dt);
	HF1(100*layer+9, hit_cont.at(i).pos-hit_cont.at(0).pos);
      }
    }// for(layer)
  }
  HF2(41, multiplicity[0], multiplicity[1]);

  if(multi_BcOut/Double_t(NumOfLayersBcOut) > MaxMultiHitBcOut){
    std::cout << "#W " << __FILE__ << " L" << __LINE__ << std::endl
	      << "multi_BcOut is too many: " << multi_BcOut << std::endl;
    return true;
  }

  HF1(1, 11.);

  // Bc Out
  //  std::cout << "==========TrackSearch BcOut============" << std::endl;

#if ExclusiveTracking
  Bool_t status_tracking = DCAna->TrackSearchBcOut(true);
#else
  Bool_t status_tracking = DCAna->TrackSearchBcOut();
#endif

#if Chi2Cut
  DCAna->ChiSqrCutBcOut(MaxChisqrBcOut);
#endif
  Int_t nt=DCAna->GetNtracksBcOut();
  event.ntrack=nt;
  HF1(10, Double_t(nt));
  HF1(40, status_tracking? Double_t(nt) : -1);
  for(Int_t it=0; it<nt; ++it){
    DCLocalTrack *tp=DCAna->GetTrackBcOut(it);
    Int_t nh=tp->GetNHit();
    Double_t chisqr=tp->GetChiSquare();
    Double_t x0=tp->GetX0(), y0=tp->GetY0();
    Double_t u0=tp->GetU0(), v0=tp->GetV0();
    Double_t theta = tp->GetTheta();
    event.chisqr[it]=chisqr;
    event.x0[it]=x0;
    event.y0[it]=y0;
    event.u0[it]=u0;
    event.v0[it]=v0;
    event.theta[it]=theta;

    HF1(11, Double_t(nh));
    HF1(12, chisqr);
    HF1(14, x0); HF1(15, y0);
    HF1(16, u0); HF1(17, v0);
    HF2(18, x0, u0); HF2(19, y0, v0);
    HF2(20, x0, y0);
    HF1(28 ,theta );
    Double_t xtgt=tp->GetX(1800.), ytgt=tp->GetY(1800.);
    Double_t utgt=u0, vtgt=v0;
    HF1(21, xtgt); HF1(22, ytgt);
    HF1(23, utgt); HF1(24, vtgt);
    HF2(25, xtgt, utgt); HF2(26, ytgt, vtgt);
    HF2(27, xtgt, ytgt);

    Double_t xbac=tp->GetX(603.), ybac=tp->GetY(885.);
    Double_t ubac=u0, vbac=v0;
    HF1(31, xbac); HF1(32, ybac);
    HF1(33, ubac); HF1(34, vbac);
    HF2(35, xbac, ubac); HF2(36, ybac, vbac);
    HF2(37, xbac, ybac);

    Double_t Xangle = -1000*std::atan(u0);

    HF2(51, -tp->GetX(245.), Xangle);
    HF2(52, -tp->GetX(600.), Xangle);
    HF2(53, -tp->GetX(1200.), Xangle);
    HF2(54, -tp->GetX(1600.), Xangle);

    HF2(61, tp->GetX(600.), tp->GetX(245.));
    HF2(62, tp->GetX(1200.), tp->GetX(245.));
    HF2(63, tp->GetX(1600.), tp->GetX(245.));

    for(Int_t ih=0; ih<nh; ++ih){
      DCLTrackHit *hit=tp->GetHit(ih);
      Int_t layerId=hit->GetLayer()-112;

      HF1(13, layerId);
      Double_t wire=hit->GetWire();
      Double_t dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1(100*layerId+11, wire-0.5);
      HF1(100*layerId+12, dt);
      HF1(100*layerId+13, dl);
      HF1(10000*layerId+ 5000 +(Int_t)wire, dt);
      Double_t xcal=hit->GetXcal(), ycal=hit->GetYcal();
      Double_t pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      Double_t exres=hit->GetResidualExclusive();
      HF1(100*layerId+14, pos);
      HF1(100*layerId+15, res);
      HF2(100*layerId+16, pos, res);
      HF2(100*layerId+17, xcal, ycal);
      //      HF1(100000*layerId+50000+wire, res);
      Double_t wp=hit->GetWirePosition();
      Double_t sign=1.;
      if(pos-wp<0.) sign=-1;
      HF2(100*layerId+18, sign*dl, res);
      Double_t xlcal=hit->GetLocalCalPos();
      HF2(100*layerId+19, dt, xlcal-wp);
      Double_t tot = hit->GetTot();
      HF1(100*layerId+21, exres);
      HF2(100*layerId+23, dt, res);
      HF2(100*layerId+24, xlcal-wp, dt);
      HF1(100*layerId+40, tot);

      if (theta>=0 && theta<15){
	HF1(100*layerId+71, res);
	HF1(100*layerId+81, exres);
      }
      else if (theta>=15 && theta<30){
	HF1(100*layerId+72, res);
	HF1(100*layerId+82, exres);
      }
      else if (theta>=30 && theta<45){
	HF1(100*layerId+73, res);
	HF1(100*layerId+83, exres);
      }
      else if (theta>=45){
	HF1(100*layerId+74, res);
	HF1(100*layerId+84, exres);
      }

      if (std::abs(dl-std::abs(xlcal-wp))<2.0) {
	HFProf(100*layerId+20, dt, std::abs(xlcal-wp));
	HF2(100*layerId+22, dt, std::abs(xlcal-wp));
	HFProf(100000*layerId+3000+(Int_t)wire, xlcal-wp,dt);
	HF2(100000*layerId+4000+(Int_t)wire, xlcal-wp,dt);
      }
    }
  }

  HF1(1, 12.);

  return true;
}

//_____________________________________________________________________________
Bool_t
UserBcOutTracking::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserBcOutTracking;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  const Int_t NbinBcOutTdc   = 1000;
  const Double_t MinBcOutTdc =  200.;
  const Double_t MaxBcOutTdc = 1200.;

  const Int_t NbinBcOutDT   =  96;
  const Double_t MinBcOutDT =  -10.;
  const Double_t MaxBcOutDT =   70.;

  const Int_t NbinBcOutDL   = 80;
  const Double_t MinBcOutDL =  -0.5;
  const Double_t MaxBcOutDL =   3.5;

  HB1(1, "Status", 20, 0., 20.);

  //***********************Chamber
  // BC34
  for(Int_t i=1; i<=NumOfLayersBcOut; ++i){
    TString tag = (i <= NumOfLayersBc) ? "BC3" : "BC4";
    Int_t nwire = (i <= NumOfLayersBc) ? MaxWireBC3 : MaxWireBC4;
    TString title0 = Form("#Hits %s#%2d", tag.Data(), i);
    TString title1 = Form("Hitpat %s#%2d", tag.Data(), i);
    TString title2 = Form("Tdc %s#%2d", tag.Data(), i);
    TString title3 = Form("Drift Time %s#%2d", tag.Data(), i);
    TString title4 = Form("Drift Length %s#%2d", tag.Data(), i);
    TString title5 = Form("TOT %s#%2d", tag.Data(), i);
    TString title6 = Form("Tdc 1st %s#%2d", tag.Data(), i);
    TString title7 = Form("TOT 1st %s#%2d", tag.Data(), i);
    TString title8 = Form("Time interval from 1st hit %s#%2d", tag.Data(), i);
    TString title9 = Form("Position interval from 1st hit %s#%2d", tag.Data(), i);
    TString title10 = Form("Trailing %s#%2d", tag.Data(), i);
    HB1(100*i+0, title0, nwire+1, 0., nwire+1);
    HB1(100*i+1, title1, nwire, 0., nwire);
    HB1(100*i+2, title2, NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc);
    HB1(100*i+3, title3, NbinBcOutDT, MinBcOutDT, MaxBcOutDT);
    HB1(100*i+4, title4, NbinBcOutDL, MinBcOutDL, MaxBcOutDL);
    HB1(100*i+5, title5, 500, 0, 500);
    HB1(100*i+6, title6, NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc);
    HB1(100*i+7, title7, 500, 0, 500);
    HB1(100*i+8, title8, 72, 0, 60);
    HB1(100*i+9, title9, 64, -32, 32);
    HB1(100*i+10, title10, NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc);
    for (Int_t wire=1; wire<=nwire; wire++) {
      TString title10 = Form("Tdc %s#%2d Wire#%d", tag.Data(), i, wire);
      TString title11 = Form("Drift Time %s#%2d Wire#%d", tag.Data(), i, wire);
      TString title12 = Form("Drift Length %s#%2d Wire#%d", tag.Data(), i, wire);
      TString title15 = Form("Drift Time %s#%2d Wire#%d [Track]", tag.Data(), i, wire);
      HB1(10000*i+wire, title10, NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc);
      HB1(10000*i+1000+wire, title11, NbinBcOutDT, MinBcOutDT, MaxBcOutDT);
      HB1(10000*i+2000+wire, title12, NbinBcOutDL, MinBcOutDL, MaxBcOutDL);
      HB1(10000*i+5000+wire, title15, NbinBcOutDT, MinBcOutDT, MaxBcOutDT);
    }

    // Tracking Histgrams
    TString title11 = Form("HitPat BcOut%2d [Track]", i);
    TString title12 = Form("DriftTime BcOut%2d [Track]", i);
    TString title13 = Form("DriftLength BcOut%2d [Track]", i);
    TString title14 = Form("Position BcOut%2d", i);
    TString title15 = Form("Residual BcOut%2d", i);
    TString title16 = Form("Resid%%Pos BcOut%2d", i);
    TString title17 = Form("Y%%Xcal BcOut%2d", i);
    TString title18 = Form("Res%%dl BcOut%2d", i);
    TString title19 = Form("HitPos%%DriftTime BcOut%2d", i);
    TString title20 = Form("DriftLength%%DriftTime BcOut%2d", i);
    TString title21 = title15 + " [w/o Self]";
    TString title22 = title20;
    TString title23 = Form("Res%%dt BcOut%2d", i);
    TString title24 = Form("DriftTime%%HitPos BcOut%2d", i);
    TString title40 = Form("TOT BcOut%2d [Track]", i);
    TString title71 = Form("Residual BcOut%2d (0<theta<15)", i);
    TString title72 = Form("Residual BcOut%2d (15<theta<30)", i);
    TString title73 = Form("Residual BcOut%2d (30<theta<45)", i);
    TString title74 = Form("Residual BcOut%2d (45<theta)", i);
    TString title81 = title71 + " [w/o Self]";
    TString title82 = title72 + " [w/o Self]";
    TString title83 = title73 + " [w/o Self]";
    TString title84 = title74 + " [w/o Self]";

    HB1(100*i+11, title11, 64, 0., 64.);
    HB1(100*i+12, title12, NbinBcOutDT, MinBcOutDT, MaxBcOutDT);
    HB1(100*i+13, title13, NbinBcOutDL, MinBcOutDL, MaxBcOutDL);
    HB1(100*i+14, title14, 100, -250., 250.);
    HB1(100*i+15, title15, 500, -5.0, 5.0);
    HB2(100*i+16, title16, 250, -250., 250., 100, -1.0, 1.0);
    HB2(100*i+17, title17, 100, -250., 250., 100, -250., 250.);
    HB2(100*i+18, title18, 40, -2., 2., 100, -1.0, 1.0);
    HB2(100*i+19, title19, NbinBcOutDT, MinBcOutDT, MaxBcOutDT,
        Int_t(MaxBcOutDL*20), -MaxBcOutDL, MaxBcOutDL);
    HBProf(100*i+20, title20, 100, -5., 50., 0., 12.);
    HB1(100*i+21, title21, 500, -5.0, 5.0);
    HB2(100*i+22, title22,
        NbinBcOutDT, MinBcOutDT, MaxBcOutDT,
        NbinBcOutDL, MinBcOutDL, MaxBcOutDL);
    HB2(100*i+23, title23,
	NbinBcOutDT, MinBcOutDT, MaxBcOutDT,
	100, -1.0, 1.0);
    HB2(100*i+24, title24,
	Int_t(MaxBcOutDL*20), -MaxBcOutDL, MaxBcOutDL,
	NbinBcOutDT, MinBcOutDT, MaxBcOutDT);
    HB1(100*i+40, title40, 300,    0, 300);
    HB1(100*i+71, title71, 500, -5.0, 5.0);
    HB1(100*i+72, title72, 500, -5.0, 5.0);
    HB1(100*i+73, title73, 500, -5.0, 5.0);
    HB1(100*i+74, title74, 500, -5.0, 5.0);
    HB1(100*i+81, title81, 500, -5.0, 5.0);
    HB1(100*i+82, title82, 500, -5.0, 5.0);
    HB1(100*i+83, title83, 500, -5.0, 5.0);
    HB1(100*i+84, title84, 500, -5.0, 5.0);

    for (Int_t j=1; j<=64; j++) {
      TString title = Form("XT of Layer %2d Wire #%4d", i, j);
      HBProf(100000*i+3000+j, title, 100, -4., 4., -5., 40.);
      HB2(100000*i+4000+j, title, 100, -4., 4., 100, -5., 40.);
    }
  }

  // Tracking Histgrams
  HB1(10, "#Tracks BcOut", 10, 0., 10.);
  HB1(11, "#Hits of Track BcOut", 15, 0., 15.);
  HB1(12, "Chisqr BcOut", 500, 0., 50.);
  HB1(13, "LayerId BcOut", 15, 0., 15.);
  HB1(14, "X0 BcOut", 400, -100., 100.);
  HB1(15, "Y0 BcOut", 400, -100., 100.);
  HB1(16, "U0 BcOut", 200, -0.20, 0.20);
  HB1(17, "V0 BcOut", 200, -0.20, 0.20);
  HB2(18, "U0%X0 BcOut", 100, -100., 100., 100, -0.20, 0.20);
  HB2(19, "V0%Y0 BcOut", 100, -100., 100., 100, -0.20, 0.20);
  HB2(20, "X0%Y0 BcOut", 100, -100., 100., 100, -100, 100);

  HB1(21, "Xtgt BcOut", 400, -100., 100.);
  HB1(22, "Ytgt BcOut", 400, -100., 100.);
  HB1(23, "Utgt BcOut", 200, -0.20, 0.20);
  HB1(24, "Vtgt BcOut", 200, -0.20, 0.20);
  HB2(25, "Utgt%Xtgt BcOut", 100, -100., 100., 100, -0.20, 0.20);
  HB2(26, "Vtgt%Ytgt BcOut", 100, -100., 100., 100, -0.20, 0.20);
  HB2(27, "Xtgt%Ytgt BcOut", 100, -100., 100., 100, -100, 100);
  HB1(28, "theta", 180, 0, 180);

  HB1(31, "Xbac BcOut", 400, -100., 100.);
  HB1(32, "Ybac BcOut", 400, -100., 100.);
  HB1(33, "Ubac BcOut", 200, -0.20, 0.20);
  HB1(34, "Vbac BcOut", 200, -0.20, 0.20);
  HB2(35, "Ubac%Xbac BcOut", 100, -100., 100., 100, -0.20, 0.20);
  HB2(36, "Vbac%Ybac BcOut", 100, -100., 100., 100, -0.20, 0.20);
  HB2(37, "Xbac%Ybac BcOut", 100, -100., 100., 100, -100, 100);
  HB1(38, "Plane Eff", 36, 0, 36);

  // HB2(1000*i, Form("Wire%%Tdc for LayerId = %d", i),
  // 	 NbinSdcOutTdc/4, MinSdcOutTdc, MaxSdcOutTdc,
  // 	 MaxWire+1, 0., double(MaxWire+1));

  HB2(51, "X-X' 245 BcOut", 400, -100., 100., 120, -60, 60);
  HB2(52, "X-X' 600 BcOut", 400, -100., 100., 120, -60, 60);
  HB2(53, "X-X' 1200 BcOut", 400, -100., 100., 120, -60, 60);
  HB2(54, "X-X' 1600 BcOut", 400, -100., 100., 120, -60, 60);

  HB2(61, "X-X 600 BcOut", 400, -100., 100., 400, -100, 100);
  HB2(62, "X-X 1200 BcOut", 400, -100., 100., 400, -100, 100);
  HB2(63, "X-X 1600 BcOut", 400, -100., 100., 400, -100, 100);

  // Analysis status
  HB1(40, "Tacking status", 11, -1., 10.);
  HB2(41, "BC3X0/BC4X1", 20, 0, 20, 20, 0, 20);

  ////////////////////////////////////////////
  //Tree
  HBTree("bcout", "tree of BcOutTracking");
  tree->Branch("evnum", &event.evnum, "evnum/I");
  tree->Branch("trigpat", event.trigpat, Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag", event.trigflag, Form("trigflag[%d]/I", NumOfSegTrig));

  //Hodoscope
  tree->Branch("nhBh1",    &event.nhBh1,   "nhBh1/I");
  tree->Branch("tBh1",      event.tBh1,    Form("tBh1[%d]/D",   MaxHits));
  tree->Branch("deBh1",     event.deBh1,   Form("deBh1[%d]/D",  MaxHits));

  tree->Branch("nhBh2",    &event.nhBh2,   "nhBh2/I");
  tree->Branch("tBh2",      event.tBh2,    Form("tBh2[%d]/D",   MaxHits));
  tree->Branch("deBh2",     event.deBh2,   Form("deBh2[%d]/D",  MaxHits));
  tree->Branch("Bh2Seg",    event.Bh2Seg,  Form("Bh2Seg[%d]/D", MaxHits));

  tree->Branch("Time0Seg", &event.Time0Seg,  "Time0Seg/D");
  tree->Branch("deTime0",  &event.deTime0,   "deTime0/D");
  tree->Branch("Time0",    &event.Time0,     "Time0/D");
  tree->Branch("CTime0",   &event.CTime0,    "CTime0/D");

  tree->Branch("nhBac",    &event.nhBac,     "nhBac/I");

  tree->Branch("btof",     &event.btof,     "btof/D");

  tree->Branch("nhit",     &event.nhit,     Form("nhit[%d]/I", NumOfLayersBcOut));
  tree->Branch("nlayer",   &event.nlayer,   "nlayer/I");
  tree->Branch("pos",      &event.pos,      Form("pos[%d][%d]/D", NumOfLayersBcOut, MaxHits));
  tree->Branch("ntrack",   &event.ntrack,   "ntrack/I");
  tree->Branch("chisqr",    event.chisqr,   "chisqr[ntrack]/D");
  tree->Branch("x0",        event.x0,       "x0[ntrack]/D");
  tree->Branch("y0",        event.y0,       "y0[ntrack]/D");
  tree->Branch("u0",        event.u0,       "u0[ntrack]/D");
  tree->Branch("v0",        event.v0,       "v0[ntrack]/D");
  tree->Branch("theta",     event.theta,    "theta[ntrack]/D");
  tree->Branch("tdc1st",    event.tdc1st,   Form("tdc1st[%d][%d]/I",
						 NumOfLayersBcOut, MaxHits));
  tree->Branch("tot1st",    event.tot1st,   Form("tot1st[%d][%d]/I",
						 NumOfLayersBcOut, MaxHits));
  tree->Branch("wire",      event.wire,     Form("wire[%d][%d]/I",
						 NumOfLayersBcOut, MaxHits));

  // HPrint();

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")        &&
     InitializeParameter<DCDriftParamMan>("DCDRFT") &&
     InitializeParameter<DCTdcCalibMan>("DCTDC")    &&
     InitializeParameter<HodoParamMan>("HDPRM")     &&
     InitializeParameter<HodoPHCMan>("HDPHC")       &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

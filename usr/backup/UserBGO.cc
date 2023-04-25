/**
 *  file: UserBGO.cc
 *  date: 2017.04.10
 *
 */

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
#include "VEvent.hh"
#include "CFTPedCorMan.hh"
#include "FLHit.hh"
#include "BGOAnalyzer.hh"

#define HodoCut 0 // with BH1/BH2
#define TimeCut 0 // in cluster analysis

namespace
{
  using namespace root;
  const std::string& classname("EventBGO");
  RMAnalyzer&         gRM   = RMAnalyzer::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
}

//______________________________________________________________________________
VEvent::VEvent( void )
{
}

//______________________________________________________________________________
VEvent::~VEvent( void )
{
}

//______________________________________________________________________________
class EventBGO : public VEvent
{
private:
  RawData      *rawData;
  BGOAnalyzer  *bgoAna;
  HodoAnalyzer *hodoAna;

public:
        EventBGO( void );
       ~EventBGO( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventBGO::EventBGO( void )
  : VEvent(),
    rawData(0),
    bgoAna( new BGOAnalyzer ),
    hodoAna( new HodoAnalyzer )
{
}

//______________________________________________________________________________
EventBGO::~EventBGO( void )
{
  if ( hodoAna ){
    delete hodoAna;
    hodoAna = 0;
  }


  if ( bgoAna ){
    delete bgoAna;
    bgoAna = 0;
  }

  if ( rawData ){
    delete rawData;
    rawData = 0;
  }
}

//______________________________________________________________________________


#ifndef MaxHits2
#define MaxHits2 60
#endif


struct Event
{
  int evnum;

  int trignhits;
  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];

  int nPulse[NumOfSegBGO];
  double chi2[NumOfSegBGO];
  double Height[NumOfSegBGO];
  double Energy[NumOfSegBGO];
  double Time[NumOfSegBGO];

  double TimeHUL[NumOfSegBGO];

};

//______________________________________________________________________________
namespace root
{
  Event  event;
  TH1   *h[MaxHist];
  TTree *tree;
}

//______________________________________________________________________________
bool
EventBGO::ProcessingBegin( void )
{
  InitializeEvent();
  return true;
}

//______________________________________________________________________________
bool
EventBGO::ProcessingNormal( void )
{
  static const std::string funcname("["+classname+"::"+__func__+"]");

#if HodoCut
  static const double MinDeBH2   = gUser.GetParameter("DeBH2", 0);
  static const double MaxDeBH2   = gUser.GetParameter("DeBH2", 1);
  static const double MinBeamToF = gUser.GetParameter("BTOF",  0);
  static const double MaxBeamToF = gUser.GetParameter("BTOF",  1);
#endif

  rawData = new RawData;
  rawData->DecodeHits();

  gRM.Decode();

  event.evnum = gRM.EventNumber();

  HF1(1, 0);

  // Trigger Flag
  {
    const HodoRHitContainer &cont=rawData->GetTrigRawHC();
    int trignhits = 0;
    int nh = cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      int seg = hit->SegmentId()+1;
      int tdc = hit->GetTdc1();
      if( tdc ){
	event.trigpat[trignhits++] = seg;
	event.trigflag[seg-1]      = tdc;
	HF1( 10, seg-1 );
	HF1( 10+seg, tdc );
      }
    }
    event.trignhits = trignhits;
  }

  // if( event.trigflag[SpillEndFlag] ) return true;

  HF1(1, 1);

  ////////// BH2 time 0
  hodoAna->DecodeBH2Hits(rawData);
  int nhBh2 = hodoAna->GetNHitsBH2();
#if HodoCut
  if ( nhBh2==0 ) return true;
#endif
  HF1(1, 2);
  double time0 = -999;
  ////////// BH2 Analysis
  for(int i=0; i<nhBh2; ++i){
    BH2Hit *hit = hodoAna->GetHitBH2(i);
    if(!hit) continue;
    double cmt = hit->CMeanTime();
    double ct0 = hit->CTime0();
    double min_time = -999;
#if HodoCut
    double dE  = hit->DeltaE();
    if( dE<MinDeBH2 || MaxDeBH2<dE )
      continue;
#endif
    if( std::abs(cmt)<std::abs(min_time) ){
      min_time = cmt;
      time0    = ct0;
    }
  }

  HF1(1, 3);

  for (int seg=0; seg<NumOfSegBGO; seg++) {

    const FADCRHitContainer &cont = rawData->GetBGOFAdcRawHC(seg);
    int nhit = cont.size();

    int hid = (seg+1)*100;

    for (int i=0; i<nhit; i++) {
      if (cont[i] < 0xffff) {
	HF2(hid, i, cont[i]);
      }
    }
  }



  bgoAna->DecodeBGO(rawData);
  bgoAna->PulseSearch();

  for (int seg=2; seg<NumOfSegBGO; seg++) {
    const FAdcDataContainer &cont = bgoAna->GetBGOFadcCont(seg);
    int hid = (seg+1)*100+10;

    for (int i=0; i<cont.size(); i++) {
      if (cont[i].err<100)
	HF2(hid, cont[i].x, cont[i].y);
    }

    const BGODataContainer &bgoCont = bgoAna->GetBGODataCont(seg);
    event.nPulse[seg] = bgoAna->GetNPulse(seg);

    BGOData bgoData;
    if (bgoAna->GetBGOData0(seg, bgoData)) {
      event.chi2[seg]   = bgoData.chi2;
      event.Height[seg] = bgoData.pulse_height;
      event.Energy[seg] = bgoData.energy;
      event.Time[seg]   = bgoData.time;
    }

    for (int i=0; i<bgoCont.size(); i++) {
      const BGOData bgoData = bgoCont[i];

      hid = (seg+1)*100+1;
      HF1(hid, bgoData.pulse_height);

      hid = (seg+1)*100+11;
      HF1(hid, bgoData.energy);

      hid = (seg+1)*100+12;
      HF1(hid, bgoData.time);
    }
  }

  hodoAna->DecodeBGOHits(rawData, bgoAna);

  int nhBGO = hodoAna->GetNHitsBGO();

  for(int i=0; i<nhBGO; ++i){
    Hodo1Hit *hit = hodoAna->GetHitBGO(i);
    if(!hit) continue;

    int seg = hit->SegmentId();
    double E = hit->DeltaE();

    int nh = hit->GetNumOfHit();
    double time0 = -9999;
    for (int j=0; j<nh; j++ ) {
      double time = hit->Time(j);
      int hid = (seg+1)*100+13;      
      HF1(hid, time);

      hid = (seg+1)*100+14;      
      HF2(hid, E, time);
      
      if (std::abs(time)<std::abs(time0))
	time0 = time;
    }
    event.TimeHUL[seg] = time0;

    HodoRawHit *rhit = hit->GetRawHit();
    int nhraw = rhit->SizeTdc1();
    for (int j=0; j<nhraw; j++ ) {
      int tdc = rhit->GetTdc1(j);
      int hid = (seg+1)*100+2;      
      HF1(hid, tdc);
    }

  }




  return true;
}

//______________________________________________________________________________
bool
EventBGO::ProcessingEnd( void )
{
  tree->Fill();
  return true;
}

//______________________________________________________________________________
void
EventBGO::InitializeEvent( void )
{
  event.evnum      = 0;

  for( int it=0; it<NumOfSegTrig; it++){
    event.trigpat[it]  = -1;
    event.trigflag[it] = -1;
  }


  for (int i=0; i<NumOfSegBGO; i++) {
    event.nPulse[i]  = 0;
    event.chi2[i]    = -999.;
    event.Height[i]  = -999.;
    event.Energy[i]  = -999.;
    event.Time[i]    = -999.;
    event.TimeHUL[i] = -999.;

  }
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventBGO;
}

//______________________________________________________________________________
namespace
{
  const int    NbinFADC_X = 200;
  const int    NbinFADC_Y = 200;
  const double MinFADC  = 0.;
  const double MaxFADC  = 18000.;

  const double MinTime  = -4;
  const double MaxTime  = 2.;
  const double MinPulseHeight  = -18000;
  const double MaxPulseHeight  = 2000.;

  const double MaxEnergy = 200.;

  const int    NbinX = 200;
  const double MinTimeHUL  = -100;
  const double MaxTimeHUL  = 100.;
  const double MinTDC  = 0;
  const double MaxTDC  = 2000.;

}

//______________________________________________________________________________
bool
ConfMan:: InitializeHistograms( void )
{
  HB1(  1, "Status",  20,   0., 20. );
  HB1( 10, "Trigger HitPat", NumOfSegTrig, 0., double(NumOfSegTrig) );
  for(int i=0; i<NumOfSegTrig; ++i){
    HB1( 10+i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000 );
  }


  char buf[100];

  for (int seg=0; seg<NumOfSegBGO; seg++) {
    int hid = (seg+1)*100;
    sprintf(buf, "BGO Raw Waveform segment %d", seg);
    HB2( hid, buf, NbinFADC_X, 0, NbinFADC_X, NbinFADC_Y, MinFADC, MaxFADC );

    hid = (seg+1)*100+1;
    sprintf(buf, "BGO pulse height segment %d", seg);
    HB1( hid, buf, NbinFADC_X, 0, MaxFADC );

    hid = (seg+1)*100+2;
    sprintf(buf, "BGO TDC (HUL) segment %d", seg);
    HB1( hid, buf, NbinX, MinTDC,  MaxTDC);

    hid = (seg+1)*100+10;
    sprintf(buf, "BGO Waveform segment %d", seg);
    HB2( hid, buf, NbinFADC_X, MinTime, MaxTime, NbinFADC_Y, MinPulseHeight, MaxPulseHeight );

    hid = (seg+1)*100+11;
    sprintf(buf, "BGO energy segment %d", seg);
    HB1( hid, buf, NbinFADC_X, 0, MaxEnergy );

    hid = (seg+1)*100+12;
    sprintf(buf, "BGO time segment %d", seg);
    HB1( hid, buf, NbinFADC_X, MinTime, MaxTime);

    hid = (seg+1)*100+13;
    sprintf(buf, "BGO time (HUL) segment %d", seg);
    HB1( hid, buf, NbinX, MinTimeHUL, MaxTimeHUL);

    hid = (seg+1)*100+14;
    sprintf(buf, "BGO time (HUL) - energy segment %d", seg);
    HB2( hid, buf,NbinFADC_X, 0, MaxEnergy,  NbinX, MinTimeHUL, MaxTimeHUL);

  }



  //Tree
  HBTree( "bgo", "tree of Easiroc" );
  //Trig
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trignhits", &event.trignhits, "trignhits/I");
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  tree->Branch("nPulse",   event.nPulse,  Form("nPulse[%d]/I", NumOfSegBGO));
  tree->Branch("chi2",     event.chi2,    Form("chi2[%d]/D", NumOfSegBGO));
  tree->Branch("Height",   event.Height,  Form("Height[%d]/D", NumOfSegBGO));
  tree->Branch("Energy",   event.Energy,  Form("Energy[%d]/D", NumOfSegBGO));
  tree->Branch("Time",     event.Time,    Form("Time[%d]/D", NumOfSegBGO));
  tree->Branch("TimeHUL",  event.TimeHUL, Form("TimeHUL[%d]/D", NumOfSegBGO));


  HPrint();
  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")    &&
      InitializeParameter<HodoParamMan>("HDPRM") &&
      InitializeParameter<HodoPHCMan>("HDPHC")   &&
      InitializeParameter<BGOTemplateManager>("BGOTEMP") &&
      InitializeParameter<BGOCalibMan>("BGOCALIB") &&
      InitializeParameter<UserParamMan>("USER")  );
}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}

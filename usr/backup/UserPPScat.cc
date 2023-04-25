/**
 *  file: UserPPScat.cc
 *  date: 2017.04.10
 *
 */

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>

#include "ConfMan.hh"
#include "DatabasePDG.hh"
#include "DCRawHit.hh"
#include "DebugTimer.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "KuramaLib.hh"
#include "MathTools.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"
#include "K18TrackD2U.hh"
#include "BH1Match.hh"
#include "BH2Filter.hh"
#include "CFTPedCorMan.hh"
#include "BGOAnalyzer.hh"
#include "CFTParticle.hh"

#define HodoCut 0
#define UseTOF  1

namespace
{
  using namespace root;
  const std::string& classname("PPScat");
  const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
  RMAnalyzer&         gRM   = RMAnalyzer::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
  BH2Filter&          gFilter = BH2Filter::GetInstance();
//BH1Match&           gBH1Mth = BH1Match::GetInstance();
  const hddaq::unpacker::UnpackerManager& gUnpacker
  = hddaq::unpacker::GUnpacker::get_instance();
  const double ProtonMass = pdg::ProtonMass();

  const double Bh2SegX[NumOfSegBH2]      = {35./2., 10./2., 7./2., 7./2., 7./2., 7./2., 10./2., 35./2.};
  const double Bh2SegXAcc[NumOfSegBH2]   = {20., 6.5, 5., 5., 5., 5., 6.5, 20.};
  const double localPosBh2X_dX           = 0.;
  const double localPosBh2X[NumOfSegBH2] = {-41.5 + localPosBh2X_dX,
					    -19.0 + localPosBh2X_dX,
					    -10.5 + localPosBh2X_dX,
					    -3.5  + localPosBh2X_dX,
					    3.5   + localPosBh2X_dX,
					    10.5  + localPosBh2X_dX,
					    19.0  + localPosBh2X_dX,
					    41.5  + localPosBh2X_dX};


  const double Bh1SegX[NumOfSegBH1] = {30./2., 20./2., 16./2., 12./2., 8./2., 8./2., 8./2., 12./2., 16./2., 20./2., 30./2.};
  const double Bh1SegXAcc[NumOfSegBH1] = {18., 12., 10., 8., 6., 6., 6., 8., 10., 12., 18.};
  //const double Bh1SegY[NumOfSegBH1] = {5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2.};
  const double& localPosBh1Z = gGeom.LocalZ("BH1");
  const double localPosBh1X_dX = 0.;
  const double localPosBh1X[NumOfSegBH1] = {-70. + localPosBh1X_dX,
					    -46. + localPosBh1X_dX,
					    -29. + localPosBh1X_dX,
					    -16. + localPosBh1X_dX,
					    -7. + localPosBh1X_dX,
					    0. + localPosBh1X_dX,
					    7. + localPosBh1X_dX,
					    16. + localPosBh1X_dX,
					    29. + localPosBh1X_dX,
					    46. + localPosBh1X_dX,
					    70. + localPosBh1X_dX};
  const double localPosBh1_dZ[NumOfSegBH1] = {4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5};
  const double& zBFT    = gGeom.LocalZ("BFT");

  const double offsetCATCH = 155;
  //const double offsetCATCH = 170;
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
class EventPPScat : public VEvent
{
private:
  RawData      *rawData;
  DCAnalyzer   *DCAna;
  HodoAnalyzer *hodoAna;
  BGOAnalyzer  *bgoAna;

public:
        EventPPScat( void );
       ~EventPPScat( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventPPScat::EventPPScat( void )
  : VEvent(),
    rawData(0),
    DCAna( new DCAnalyzer ),
    hodoAna( new HodoAnalyzer ),
    bgoAna( new BGOAnalyzer )
{
}

//______________________________________________________________________________
EventPPScat::~EventPPScat( void )
{
  if( DCAna )   delete DCAna;
  if( hodoAna ) delete hodoAna;
  if (bgoAna)  delete bgoAna;
  if( rawData ) delete rawData;
}

//______________________________________________________________________________
struct Event
{
  int evnum;

  int trignhits;
  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];

  double btof;
  double time0;

  int nhBh2;
  double Bh2Seg[MaxHits];
  double tBh2[MaxHits];
  double t0Bh2[MaxHits];
  double deBh2[MaxHits];

  int nhBh1;
  double Bh1Seg[MaxHits];
  double tBh1[MaxHits];
  double deBh1[MaxHits];

  // BFT
  int    bft_ncl;
  int    bft_ncl_bh1mth;
  int    bft_clsize[NumOfSegBFT];
  double bft_ctime[NumOfSegBFT];
  double bft_clpos[NumOfSegBFT];
  int    bft_bh1mth[NumOfSegBFT];

  // BcOut
  int nlBcOut;
  int ntBcOut;
  int nhBcOut[MaxHits];
  double chisqrBcOut[MaxHits];
  double x0BcOut[MaxHits];
  double y0BcOut[MaxHits];
  double u0BcOut[MaxHits];
  double v0BcOut[MaxHits];

  // K18
  int ntK18;
  int nhK18[MaxHits];
  double p_2nd[MaxHits];
  double p_3rd[MaxHits];
  double delta_2nd[MaxHits];
  double delta_3rd[MaxHits];

  double xinK18[MaxHits];
  double yinK18[MaxHits];
  double uinK18[MaxHits];
  double vinK18[MaxHits];

  double xoutK18[MaxHits];
  double youtK18[MaxHits];
  double uoutK18[MaxHits];
  double voutK18[MaxHits];

  double chisqrK18[MaxHits];
  double xtgtK18[MaxHits];
  double ytgtK18[MaxHits];
  double utgtK18[MaxHits];
  double vtgtK18[MaxHits];

  double thetaK18[MaxHits];
  double phiK18[MaxHits];

  // CFT
  int ntCFT;
  double phi[MaxDepth];
  double theta[MaxDepth];

  int seg[NumOfPlaneCFT][MaxDepth];
  int seg_max[NumOfPlaneCFT][MaxDepth];

  double dphi[NumOfPlaneCFT][MaxDepth];
  double phi_ini[NumOfPlaneCFT][MaxDepth];
  double phi_track[NumOfPlaneCFT][MaxDepth];

  double dz[NumOfPlaneCFT][MaxDepth];
  double z_ini[NumOfPlaneCFT][MaxDepth];
  double z_track[NumOfPlaneCFT][MaxDepth];
  ThreeVector Pos[MaxDepth];
  ThreeVector Dir[MaxDepth];

  double adc[NumOfPlaneCFT][MaxDepth], adc_max[NumOfPlaneCFT][MaxDepth];
  double mip[NumOfPlaneCFT][MaxDepth], mip_max[NumOfPlaneCFT][MaxDepth];
  double dE[NumOfPlaneCFT][MaxDepth] , dE_max[NumOfPlaneCFT][MaxDepth] ;

  double Total_mip[MaxDepth],    Total_mip_max[MaxDepth];

  double Total_dE[MaxDepth],    Total_dE_max[MaxDepth];
  double Total_dEphi[MaxDepth], Total_dEphi_max[MaxDepth];
  double Total_dEuv[MaxDepth],  Total_dEuv_max[MaxDepth];
  /*
  //CFTEx using vtx point
  int ntCFTEx;
  double phiEx[MaxDepth];
  double thetaEx[MaxDepth];
  int nhVtxEx[MaxDepth];

  double dphiEx[NumOfPlaneCFT][MaxDepth];
  double phiEx_track[NumOfPlaneCFT][MaxDepth];
  double dzEx[NumOfPlaneCFT][MaxDepth];
  double zEx_track[NumOfPlaneCFT][MaxDepth];
  */


  // BGO
  int segBGOt[NumOfSegBGO];// matched to track
  // PiID counter
  int segPiIDt[MaxDepth];// matched to track
  // CFT Particle
  double Total_E[MaxDepth],    Total_E_max[MaxDepth];
  double P[MaxDepth],    P_max[MaxDepth];

  // BGO
  int nhBGO;
  int segBGO[NumOfSegBGO];
  double adcbgo[NumOfSegBGO];
  double tdcbgo[NumOfSegBGO];
  double energybgo[NumOfSegBGO];
  double energybgot[MaxDepth]; // belongs to CFT track#
  double timebgo[NumOfSegBGO];

  double vtx_2p_x;
  double vtx_2p_y;
  double vtx_2p_z;
  double cdist_2p;
  double theta_2p;


  // Missing Mass
  int nPP;
  double vtx_pp_x[MaxHits];
  double vtx_pp_y[MaxHits];
  double vtx_pp_z[MaxHits];
  double missmass[MaxHits];
  double missmom[MaxHits];
  double theta_pp[MaxHits];
  double cdist_pp[MaxHits];

  // Beam * CFT
  double vtx_BeamCft_x[MaxHits];
  double vtx_BeamCft_y[MaxHits];
  double vtx_BeamCft_z[MaxHits];
  double cdist_BeamCft[MaxHits];
  double theta_BeamCft[MaxHits];

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
EventPPScat::ProcessingBegin( void )
{
  InitializeEvent();
  return true;
}

//______________________________________________________________________________
bool
EventPPScat::ProcessingNormal( void )
{
  static const std::string func_name("["+classname+"::"+__func__+"]");

#if HodoCut

  static const double MinDeBH2   = gUser.GetParameter("DeBH2", 0);
  static const double MaxDeBH2   = gUser.GetParameter("DeBH2", 1);
  static const double MinDeBH1   = gUser.GetParameter("DeBH1", 0);
  static const double MaxDeBH1   = gUser.GetParameter("DeBH1", 1);
#endif
  static const double MinTimeBFT = gUser.GetParameter("TimeBFT", 0);
  static const double MaxTimeBFT = gUser.GetParameter("TimeBFT", 1);
  static const double MinTotBcOut = gUser.GetParameter("MinTotBcOut", 0);

  rawData = new RawData;
  rawData->DecodeHits();

  gRM.Decode();

  event.evnum = gRM.EventNumber();

  //TrigFlag
  {
    const HodoRHitContainer &cont=rawData->GetTrigRawHC();
    int trignhits = 0;
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg = hit->SegmentId()+1;
      int tdc = hit->GetTdc1();
      if( tdc ){
	event.trigpat[trignhits++] = seg;
	event.trigflag[seg-1]      = tdc;
      }
    }
    event.trignhits = trignhits;
  }

  HF1( 1, 0. );

  // if( event.trigflag[SpillEndFlag] ) return true;

  HF1( 1, 1. );

  //////////////BH2 Analysis
  hodoAna->DecodeBH2Hits(rawData);
  int nhBh2 = hodoAna->GetNHitsBH2();
  event.nhBh2 = nhBh2;
#if HodoCut
  if( nhBh2==0 ) return true;
#endif
  double time0 = -999.;
  double time0_seg = -1;

  //////////////BH2 Analysis
  double min_time = -999.;
  for( int i=0; i<nhBh2; ++i ){
    BH2Hit *hit = hodoAna->GetHitBH2(i);
    if(!hit) continue;
    double seg = hit->SegmentId()+1;
    double de  = hit->DeltaE();
#if HodoCut
    if( de<MinDeBH2 || MaxDeBH2<de ) continue;
#endif
    event.deBh2[i]  = de;
    event.Bh2Seg[i] = seg;

    int multi = hit->GetNumOfHit();
    for (int m=0; m<multi; m++) {
      double mt  = hit->MeanTime(m);
      double cmt = hit->CMeanTime(m);
      double ct0 = hit->CTime0(m);

      if (std::abs(event.tBh2[i]) > cmt )
	event.tBh2[i]   = cmt;
      if (std::abs(event.t0Bh2[i]) > ct0 )
	event.t0Bh2[i]   = ct0;
	
      if( std::abs(mt)<std::abs(min_time) ){
	min_time = mt;
	time0    = ct0;
	time0_seg = seg;
      }
    }
  }
  event.time0 = time0;


  //////////////BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData);
  int nhBh1 = hodoAna->GetNHitsBH1();
  event.nhBh1 = nhBh1;
#if HodoCut
  if( nhBh1==0 ) return true;
#endif

  HF1( 1, 2. );

  double btof0 = -999.;
  for(int i=0; i<nhBh1; ++i){
    Hodo2Hit *hit = hodoAna->GetHitBH1(i);
    if(!hit) continue;
    int    seg  = hit->SegmentId()+1;
    double cmt  = hit->CMeanTime();
    double dE   = hit->DeltaE();
    double btof = cmt - time0;
#if HodoCut
    if( dE<MinDeBH1 || MaxDeBH1<dE ) continue;
    if( btof<MinBeamToF && MaxBeamToF<btof ) continue;
#endif
    event.Bh1Seg[i] = seg;
    event.tBh1[i]   = cmt;
    event.deBh1[i]  = dE;
    if( std::abs(btof)<std::abs(btof0) ){
      btof0 = btof;
    }
  }

  event.btof = btof0;

  HF1( 1, 20. );

  HF1( 1, 21. );
  DCAna->DecodeBcOutHits( rawData );
  DCAna->TotCutBCOut( MinTotBcOut );
  DCAna->DriftTimeCutBC34(-20, 60);

  int ntBcOut = 0;
  
  if( time0_seg >= 1  && time0_seg <= NumOfSegBH2){
    BH2Filter::FilterList cands;
    gFilter.Apply((Int_t)time0_seg-1, *DCAna, cands);
    DCAna->TrackSearchBcOut( cands, time0_seg-1 );
    //DCAna->TrackSearchBcOut(-1);
    //DCAna->TrackSearchBcOut(time0_seg-1);
    ntBcOut = DCAna->GetNtracksBcOut();

    for( int it=ntBcOut-1; it>=0 ; --it ){
      DCLocalTrack *tp = DCAna->GetTrackBcOut( it );
      int nh = tp->GetNHit();
      double chisqr = tp->GetChiSquare();
      double u0 = tp->GetU0(),  v0 = tp->GetV0();
      double x0 = tp->GetX(0.), y0 = tp->GetY(0.);

      event.nhBcOut[it] = nh;
      event.chisqrBcOut[it] = chisqr;
      event.x0BcOut[it] = x0;
      event.y0BcOut[it] = y0;
      event.u0BcOut[it] = u0;
      event.v0BcOut[it] = v0;

    }
  }

  ntBcOut = DCAna->GetNtracksBcOut();
  event.ntBcOut = ntBcOut;

  HF1( 1, 22. );

  HF1( 1, 30. );

  HF1( 1, 31. );

  std::vector<double> BftXCont;
  ////////// BFT
  {
    hodoAna->DecodeBFTHits(rawData);
    // Fiber Cluster
    hodoAna->TimeCutBFT(MinTimeBFT, MaxTimeBFT);
    int ncl = hodoAna->GetNClustersBFT();
    event.bft_ncl = ncl;
    for( int i=0; i<ncl; ++i ){
      FiberCluster *cl = hodoAna->GetClusterBFT(i);
      if( !cl ) continue;
      int    clsize = cl->ClusterSize();
      double ctime  = cl->CMeanTime();
      double pos    = cl->MeanPosition();
      BftXCont.push_back( pos );

      event.bft_clsize[i] = clsize;
      event.bft_ctime[i]  = ctime;
      event.bft_clpos[i]  = pos;
      
      /*
      if(btof0_seg > 0 && ncl != 1){
	if(gBH1Mth.Judge(pos, btof0_seg)){
	  event.bft_bh1mth[i] = 1;
	  xCand.push_back( pos );
	}
      }else{
	xCand.push_back( pos );
      }
      */

    }
  }


  std::vector<ThreeVector> BeamPCont, BeamXCont;

  // K18TrackingD2U
  DCAna->TrackSearchK18D2U( BftXCont );
  int ntK18=DCAna->GetNTracksK18D2U();
  //if( ntK18==0 ) return true;
  ThreeVector Pos0( 0, 0, 0 );
  ThreeVector Mom0( 0, 0, 0 );

  double p0=100.;
  int ntK18_Bh1=0;
  for( int i=0; i<ntK18; ++i ){
    K18TrackD2U *tp=DCAna->GetK18TrackD2U(i);
    if(!tp) continue;
    DCLocalTrack *tpOut = tp->TrackOut();
    std::size_t nh = tpOut->GetNHit();
    double chisqr = tpOut->GetChiSquare();

    double timeBh2 = time0;
    int bh1seg=-1;
    double btof=999.;

    for( int ibh1=0; ibh1<nhBh1; ++ibh1 ){
      Hodo2Hit *hit = hodoAna->GetHitBH1(ibh1);
      if(!hit) continue;
      double seg = hit->SegmentId();

      int multi = hit->GetNumOfHit();
      for (int m=0; m<multi; m++) {
	double cmt = hit->CMeanTime(m);
	//std::cout << cmt-timeBh2 << std::endl;
	if (std::abs(cmt-timeBh2)<std::abs(btof)) {
	  btof = cmt-timeBh2;
	  bh1seg = seg;
	}
      }
    }
    if (std::abs(btof)>1)
      bh1seg=-1;
      
    if (bh1seg>=0) {
      double zbh1 =  localPosBh1Z+localPosBh1_dZ[bh1seg];
      double xi = tp->Xin(), ui = tp->Uin();      
      double xbh1 = xi + ui*(zbh1-zBFT);
      //std::cout << xbh1 << " - " << xbh1-localPosBh1X[bh1seg] << " ( z = " << zbh1 << " )" << std::endl;
      //std::cout << "BH1 match ( " << bh1seg << " ) = "  << xbh1-localPosBh1X[bh1seg] 
      //<< " [ " << Bh1SegXAcc[bh1seg] << " ]" << std::endl;

      if (std::abs(xbh1-localPosBh1X[bh1seg]) > Bh1SegXAcc[bh1seg]) {
	continue;
      }
    }

    double xt=tp->Xtgt(), yt=tp->Ytgt();
    double ut=tp->Utgt(), vt=tp->Vtgt();

    double xin=tp->Xin(), yin=tp->Yin();
    double uin=tp->Uin(), vin=tp->Vin();

    double xout=tp->Xout(), yout=tp->Yout();
    double uout=tp->Uout(), vout=tp->Vout();

    double p = tp->P3rd();
    double pt = p/std::sqrt(1.+ut*ut+vt*vt);

    ThreeVector Pos( xt, yt, 0. );
    ThreeVector Mom( pt*ut, pt*vt, pt );

    //    if (std::abs(p-0.6) <= std::abs(p0-0.6)) {
    if (std::abs(p-0.65) <= std::abs(p0-0.65)) {
      p0 = p;
      Pos0 = Pos;
      Mom0 = Mom;
    }

    double p_2nd=tp->P();
    double p_3rd=tp->P3rd();
    double delta_2nd=tp->Delta();
    double delta_3rd=tp->Delta3rd();
    double cost = 1./std::sqrt(1.+ut*ut+vt*vt);
    double theta = std::acos(cost)*math::Rad2Deg();
    double phi   = atan2( ut, vt );

    event.p_2nd[ntK18_Bh1] = p_2nd;
    event.p_3rd[ntK18_Bh1] = p_3rd;
    event.delta_2nd[ntK18_Bh1] = delta_2nd;
    event.delta_3rd[ntK18_Bh1] = delta_3rd;

    event.xinK18[ntK18_Bh1] = xin;
    event.yinK18[ntK18_Bh1] = yin;
    event.uinK18[ntK18_Bh1] = uin;
    event.vinK18[ntK18_Bh1] = vin;

    event.xoutK18[ntK18_Bh1] = xout;
    event.youtK18[ntK18_Bh1] = yout;
    event.uoutK18[ntK18_Bh1] = uout;
    event.voutK18[ntK18_Bh1] = vout;

    event.nhK18[ntK18_Bh1]     = nh;
    event.chisqrK18[ntK18_Bh1] = chisqr;
    event.xtgtK18[ntK18_Bh1]   = xt;
    event.ytgtK18[ntK18_Bh1]   = yt;
    event.utgtK18[ntK18_Bh1]   = ut;
    event.vtgtK18[ntK18_Bh1]   = vt;
    event.thetaK18[ntK18_Bh1] = theta;
    event.phiK18[ntK18_Bh1]   = phi;
    ntK18_Bh1++;
  }
  event.ntK18 = ntK18_Bh1;

  if (ntK18_Bh1 >= 1) {
    BeamPCont.push_back( Mom0 );
    BeamXCont.push_back( Pos0 );
  }

  bgoAna->DecodeBGO(rawData);
  bgoAna->PulseSearch();
  hodoAna->DecodeBGOHits(rawData, bgoAna);

  int nhBGO = bgoAna->GetNHitBGO();
  event.nhBGO = nhBGO;
  for(int i=0; i<nhBGO; ++i){
    Hodo1Hit *hit = hodoAna->GetHitBGO(i);
    int seg = hit->SegmentId();
    int nh = hit->GetNumOfHit();
    //double adc = hit->GetAUp();
    double energy = hit->DeltaE();

    double adc = -999;
    HodoRawHit *rhit = hit->GetRawHit();
    if (rhit->SizeAdc2()>0)
      adc = rhit->GetAdc2();

    event.segBGO[i] = seg;
    event.adcbgo[seg] = adc;
    event.energybgo[seg] = energy;

    double time0 = -999.;
    for(int m=0; m<nh; ++m){
      double cmt = hit->CMeanTime(m);
      if (fabs(time0)>fabs(cmt))
	time0 = cmt;
    }
    event.tdcbgo[seg] = time0;
  }

  hodoAna->DecodePiIDHits(rawData);
  //int nhit = hodoAna->GetNHitsPiID();          

  // FiberHitCFT    
  hodoAna->DecodeCFTHits(rawData);

  // Fiber Cluster
  for(int p = 0; p<NumOfPlaneCFT; ++p){
    hodoAna->TimeCutCFT(p, -10, 10); // CATCH@J-PARC  
    hodoAna->AdcCutCFT(p, 10, 4000); // CATCH@J-PARC  for proton
    //hodoAna->WidthCutCFT(p, 60, 300); // pp scattering
  }

  // CFT tracking
  DCAna->DecodeCFTHits( rawData );
  DCAna->TrackSearchCFT();

  int ntCFT=DCAna->GetNtracksCFT();// vtx limit ver.
  event.ntCFT = ntCFT;
  DCLocalTrack *tpp[2];
  ThreeVector vtx[2];
  std::vector<ThreeVector> ScatPCont, ScatXCont;

  for( int i=0; i<ntCFT; ++i ){

    DCLocalTrack *tp=DCAna->GetTrackCFT(i);
    if(i==0){tpp[0]=tp;}
    else if(i==1){tpp[1]=tp;}

    int nh   = tp->GetNHit();
    int nhUV = tp->GetNHitUV();
    //double chisqrXY=tp->GetChiSquareXY();
    //double chisqrXYZ=tp->GetChiSquareZ();
    double theta =tp->GetThetaCFT();

    ThreeVector Pos0 = tp->GetPos0();
    ThreeVector Dir = tp->GetDir();

    // CATCH z origin --> target center origin
    double x0 = Pos0.x() + (offsetCATCH - Pos0.z())*Dir.x()/Dir.z();
    double y0 = Pos0.y() + (offsetCATCH - Pos0.z())*Dir.y()/Dir.z();
    ThreeVector PosT0(x0, y0, 0);

    double A=(Dir.x()*Dir.x()+Dir.y()*Dir.y());

    double phi = -999.;
    if(Dir.x()>=0 && Dir.y()>=0){
      phi = acos(Dir.x()/sqrt(A))*math::Rad2Deg();
    }//0~90
    else if (Dir.x()<0 && Dir.y()>=0){
      phi = acos(Dir.x()/sqrt(A))*math::Rad2Deg();
    }//90~180
    else if (Dir.x()<0 && Dir.y()<0){
      phi = 360. - acos(Dir.x()/sqrt(A))*math::Rad2Deg();
    }//180~270
    else if (Dir.x()>=0 && Dir.y()<0){
      phi = 360. - acos(Dir.x()/sqrt(A))*math::Rad2Deg(); ;
    }//270~360

    event.theta[i] = theta;
    event.phi[i] = phi;

    event.Pos[i] = Pos0;
    event.Dir[i] = Dir;

    event.Total_mip[i]     = tp->GetTotalSumMIP()   ;
    event.Total_mip_max[i] = tp->GetTotalMaxMIP()   ;

    event.Total_dE[i]   = tp->GetTotalSumdE()   ;
    event.Total_dEphi[i]= tp->GetTotalSumdEphi();
    event.Total_dEuv[i] = tp->GetTotalSumdEuv ();
    event.Total_dE_max[i]   = tp->GetTotalMaxdE()   ;
    event.Total_dEphi_max[i]= tp->GetTotalMaxdEphi();
    event.Total_dEuv_max[i] = tp->GetTotalMaxdEuv ();

    // straight layer
    for(int ip=0; ip<nh; ip++){
      DCLTrackHit *hit = tp->GetHit(ip);
      int layer = hit->GetLayer();
      int seg = (int)hit->GetMeanSeg();
      int seg_max = (int)hit->GetMaxSeg();

      double phi_ini   = tp->GetPhiIni(layer);      
      double phi_track   = tp->GetPhiTrack(layer);      
      double z_track = tp->GetZTrack(layer);
      double dphi  = tp->GetdPhi(layer);

      event.seg[layer][i]     = seg;
      event.seg_max[layer][i] = seg_max;
      event.phi_ini[layer][i] = phi_ini;
      event.phi_track[layer][i] = phi_track;
      event.dphi[layer][i] = dphi;
      event.z_track[layer][i] = z_track;

      event.adc[layer][i] = tp->GetSumAdc(layer);
      event.mip[layer][i] = tp->GetSumMIP(layer);
      event.dE[layer][i]  = tp->GetSumdE (layer);
      event.adc_max[layer][i] = tp->GetMaxAdc(layer);      
      event.mip_max[layer][i] = tp->GetMaxMIP(layer);
      event.dE_max[layer][i]  = tp->GetMaxdE (layer);

      int hid = ((layer+1)*1000+seg_max)*10+1;
      HF1(hid, tp->GetMaxAdc(layer));
      hid = ((layer+1)*1000+seg_max)*10+3;
      HF1(hid, tp->GetMaxMIP(layer));

      if (theta>=60 && theta<=120) {
	hid = ((layer+1)*1000+seg_max)*10+2;
	HF1(hid, tp->GetMaxAdc(layer));
	hid = ((layer+1)*1000+seg_max)*10+4;
	HF1(hid, tp->GetMaxMIP(layer));
      } else if (theta>=20 && theta<=30) {
	hid = ((layer+1)*1000+seg_max)*10+5;
	HF1(hid, tp->GetMaxAdc(layer));
	hid = ((layer+1)*1000+seg_max)*10+6;
	HF1(hid, tp->GetMaxMIP(layer));
      } else if (theta>=30 && theta<=40) {
	hid = ((layer+1)*1000+seg_max)*10+7;
	HF1(hid, tp->GetMaxAdc(layer));
	hid = ((layer+1)*1000+seg_max)*10+8;
	HF1(hid, tp->GetMaxMIP(layer));
      }
    }

    for(int ip=0; ip<nhUV; ip++){
      DCLTrackHit *hit = tp->GetHitUV(ip);
      int layer = hit->GetLayer();
      int seg = (int)hit->GetMeanSeg();
      int seg_max = (int)hit->GetMaxSeg();
      
      double phi_track   = tp->GetPhiTrack(layer);      
      double z_track = tp->GetZTrack(layer);
      double z_ini   = tp->GetZIni(layer);      
      double dz    = tp->GetdZ(layer);

      event.seg[layer][i]     = seg;
      event.seg_max[layer][i] = seg_max;
      event.phi_track[layer][i] = phi_track;
      event.z_ini[layer][i] = z_ini;
      event.z_track[layer][i] = z_track;
      event.dz[layer][i] = dz;

      event.adc[layer][i] = tp->GetSumAdc(layer);
      event.mip[layer][i] = tp->GetSumMIP(layer);
      event.dE[layer][i]  = tp->GetSumdE (layer);
      event.adc_max[layer][i] = tp->GetMaxAdc(layer);      
      event.mip_max[layer][i] = tp->GetMaxMIP(layer);
      event.dE_max[layer][i]  = tp->GetMaxdE (layer);

      int hid = ((layer+1)*1000+seg_max)*10+1;
      HF1(hid, tp->GetMaxAdc(layer));
      hid = ((layer+1)*1000+seg_max)*10+3;
      HF1(hid, tp->GetMaxMIP(layer));

      if (theta>=60 && theta<=120) {
	hid = ((layer+1)*1000+seg_max)*10+2;
	HF1(hid, tp->GetMaxAdc(layer));
	hid = ((layer+1)*1000+seg_max)*10+4;
	HF1(hid, tp->GetMaxMIP(layer));
      } else if (theta>=20 && theta<=30) {
	hid = ((layer+1)*1000+seg_max)*10+5;
	HF1(hid, tp->GetMaxAdc(layer));
	hid = ((layer+1)*1000+seg_max)*10+6;
	HF1(hid, tp->GetMaxMIP(layer));
      } else if (theta>=30 && theta<=40) {
	hid = ((layer+1)*1000+seg_max)*10+7;
	HF1(hid, tp->GetMaxAdc(layer));
	hid = ((layer+1)*1000+seg_max)*10+8;
	HF1(hid, tp->GetMaxMIP(layer));
      } 
    }    

    CFTParticle * CFTPart = new CFTParticle(tp, rawData);  
    int segBGOt  = CFTPart->GetTrackBGOSeg(); // BGO  track segment
    int segPiIDt = CFTPart->GetTrackPiIDSeg();// PiID track segment

    event.segBGOt[i]  = segBGOt;
    event.segPiIDt[i] = segPiIDt;
    event.energybgot[i]  = event.energybgo[event.segBGOt[i]];

    double BGO_energy = 0;
    if (segBGOt>=2 && segBGOt<=23)
      BGO_energy = event.energybgo[event.segBGOt[i]];

    event.Total_E[i]    = event.Total_dE[i] + BGO_energy; 
    event.Total_E_max[i] = event.Total_dE_max[i] + BGO_energy; 

    event.P[i]    = sqrt((event.Total_E[i]*0.001 + ProtonMass)*(event.Total_E[i]*0.001 + ProtonMass) - ProtonMass*ProtonMass);
    event.P_max[i]    = sqrt((event.Total_E_max[i]*0.001 + ProtonMass)*(event.Total_E_max[i]*0.001 + ProtonMass) - ProtonMass*ProtonMass);

    ThreeVector Mom = event.P[i]/Dir.Mag()*Dir;
    ScatPCont.push_back(Mom);
    ScatXCont.push_back(PosT0);


    // vtx Beam*Cft
#if 1
    if(BeamPCont.size()>0){
      ThreeVector pBeam ;
      ThreeVector xBeam ;
      ThreeVector pCft = ScatPCont[i];
      ThreeVector xCft = ScatXCont[i];
      double cdist_bc=999, cdist_min=999; 
      int min_k = 0; // Beam of min. cdist 
      for(int k=0; k<BeamPCont.size(); k++){
	ThreeVector pBeam_ = BeamPCont[k];
	ThreeVector xBeam_ = BeamXCont[k];	
	double cdist_;
	ThreeVector vertex_ = Kinematics::VertexPoint( xBeam_, xCft, pBeam_, pCft, cdist_ );	
	if(cdist_<cdist_min){
	  min_k = k;
	  cdist_min = cdist_;
	}	

      }      
      if(cdist_min<999){
	pBeam = BeamPCont[min_k];
	xBeam = BeamXCont[min_k];
	ThreeVector vertex = Kinematics::VertexPoint( xBeam, xCft, pBeam, pCft, cdist_bc);        
	double cost_bc = pBeam*pCft/(pBeam.Mag()*pCft.Mag());
	double theta_bc = std::acos(cost_bc)*math::Rad2Deg();      
	
	event.vtx_BeamCft_x[i] = vertex.x();
	event.vtx_BeamCft_y[i] = vertex.y();
	event.vtx_BeamCft_z[i] = vertex.z();
	event.cdist_BeamCft[i] = cdist_bc;
	event.theta_BeamCft[i] = theta_bc;
	if(ntCFT==2){vtx[i] = vertex;}
      }
    }
#endif

  }
  
#if 0
  // ReTracking using vtx obtained from oposite track
  if(ntCFT==2){    
    ThreeVector vtx0_(vtx[0].x(), vtx[0].y(), vtx[0].z()+offsetCATCH);
    ThreeVector vtx1_(vtx[1].x(), vtx[1].y(), vtx[1].z()+offsetCATCH);

    for(int iEx=0; iEx<2; iEx++){

      if(iEx==0){
	DCAna->DecodeCFTExHits(rawData,tpp[0],vtx1_);
	DCAna->TrackSearchCFTEx();
      }else if(iEx==1){
	DCAna->DecodeCFTExHits(rawData,tpp[1],vtx0_);
	DCAna->TrackSearchCFTEx();
      }

      int ntCFTEx = DCAna->GetNtracksCFTEx();
      event.ntCFTEx += ntCFTEx;      
      for( int i=0; i<ntCFTEx; ++i ){      
	DCLocalTrack *tp=DCAna->GetTrackCFTEx(i);
	int nh   = tp->GetNHit();
	int nhUV = tp->GetNHitUV();
	int nhVtx = tp->GetNHitVtx();
	double chisqrXY=tp->GetChiSquareXY();
	double chisqrXYZ=tp->GetChiSquareZ();
	double theta =tp->GetThetaCFT();
	
	ThreeVector Pos0 = tp->GetPos0();
	ThreeVector Dir = tp->GetDir();      
	// CATCH z origin --> target center origin
	double x0 = Pos0.x() + (offsetCATCH - Pos0.z())*Dir.x()/Dir.z();
	double y0 = Pos0.y() + (offsetCATCH - Pos0.z())*Dir.y()/Dir.z();
	ThreeVector PosT0(x0, y0, 0);      
	double A=(Dir.x()*Dir.x()+Dir.y()*Dir.y());      
	double phi = -999.;
	if(Dir.x()>=0 && Dir.y()>=0){
	  phi = acos(Dir.x()/sqrt(A))*math::Rad2Deg();
	}//0~90
	else if (Dir.x()<0 && Dir.y()>=0){
	  phi = acos(Dir.x()/sqrt(A))*math::Rad2Deg();
	}//90~180
	else if (Dir.x()<0 && Dir.y()<0){
	  phi = 360. - acos(Dir.x()/sqrt(A))*math::Rad2Deg();
	}//180~270
	else if (Dir.x()>=0 && Dir.y()<0){
	  phi = 360. - acos(Dir.x()/sqrt(A))*math::Rad2Deg(); ;
	}//270~360

	event.thetaEx[iEx] = theta;
	event.phiEx[iEx] = phi;
	event.nhVtxEx[iEx] = nhVtx;

	// straight layer
	for(int ip=0; ip<nh; ip++){
	  DCLTrackHit *hit = tp->GetHit(ip);
	  int layer = hit->GetLayer();
	  int seg = (int)hit->GetMeanSeg();
	  int seg_max = (int)hit->GetMaxSeg();
	  
	  double dphi  = tp->GetdPhi(layer);
	  double phi_track   = tp->GetPhiTrack(layer);      
	  double z_track = tp->GetZTrack(layer);
	  
	  event.dphiEx[layer][iEx] = dphi;
	  event.phiEx_track[layer][iEx] = phi_track;
	  event.zEx_track[layer][iEx] = z_track;
	}

	// spiral layer
	for(int ip=0; ip<nhUV; ip++){
	  DCLTrackHit *hit = tp->GetHitUV(ip);
	  int layer = hit->GetLayer();
	  int seg = (int)hit->GetMeanSeg();
	  int seg_max = (int)hit->GetMaxSeg();

	  double dz    = tp->GetdZ(layer);	  
	  double phi_track   = tp->GetPhiTrack(layer);      
	  double z_track = tp->GetZTrack(layer);
	  
	  event.dzEx[layer][iEx] = dz;	  
	  event.phiEx_track[layer][iEx] = phi_track;
	  event.zEx_track[layer][iEx] = z_track;
	}
      }
    }
  }
#endif


  if (ScatPCont.size() == 2) {
    ThreeVector pScat1 = ScatPCont[0];
    ThreeVector xScat1 = ScatXCont[0];

    ThreeVector pScat2 = ScatPCont[1];
    ThreeVector xScat2 = ScatXCont[1];

    double cdist;
    ThreeVector vertex = Kinematics::VertexPoint( xScat1, xScat2, pScat1, pScat2, cdist );
    double cost = pScat1*pScat2/(pScat1.Mag()*pScat2.Mag());
    double theta = std::acos(cost)*math::Rad2Deg();


    event.vtx_2p_x = vertex.x();
    event.vtx_2p_y = vertex.y();
    event.vtx_2p_z = vertex.z();
    event.cdist_2p = cdist;
    event.theta_2p = theta;
  }


  event.nPP = ScatPCont.size()*BeamPCont.size();

  int iPP=0;
  for (int is = 0; is < ScatPCont.size(); is++ ) {
    for (int ib = 0; ib < BeamPCont.size(); ib++ ) {
      ThreeVector pScat = ScatPCont[is];
      ThreeVector pBeam = BeamPCont[ib];
      ThreeVector xScat = ScatXCont[is];
      ThreeVector xBeam = BeamXCont[ib];

      double cdist;
      ThreeVector vertex = Kinematics::VertexPoint( xBeam, xScat, pBeam, pScat, cdist );
      LorentzVector LvBeam( BeamPCont[ib], std::sqrt( ProtonMass*ProtonMass+pBeam.Mag2()) );
      LorentzVector LvScat( ScatPCont[is], std::sqrt( ProtonMass*ProtonMass+pScat.Mag2() ) );
      LorentzVector LvP( 0., 0., 0., ProtonMass );
      LorentzVector LvRp = LvBeam+LvP-LvScat;
      ThreeVector MissMom = LvRp.Vect();
      
      double MissMass = LvRp.Mag();

      double cost = pBeam*pScat/(pBeam.Mag()*pScat.Mag());
      double theta = std::acos(cost)*math::Rad2Deg();

      event.vtx_pp_x[iPP] = vertex.x();
      event.vtx_pp_y[iPP] = vertex.y();
      event.vtx_pp_z[iPP] = vertex.z();
      event.missmass[iPP] = MissMass;
      event.missmom[iPP] = MissMom.Mag();
      event.theta_pp[iPP] = theta;
      event.cdist_pp[iPP] = cdist;

      iPP++;

    }
  }

  return true;
}

//______________________________________________________________________________
void
EventPPScat::InitializeEvent( void )
{
  event.nhBh2    = 0;
  event.nhBh1    = 0;

  event.time0 = -9999.;
  event.btof  = -9999.;

  event.bft_ncl   =  0;
  event.bft_ncl_bh1mth =  0;
  event.nlBcOut   =  0;
  event.ntBcOut   =  0;
  event.ntK18     =  0;
  event.nPP      =  0;

  for( int it=0; it<NumOfSegTrig; it++){
    event.trigpat[it] = -1;
    event.trigflag[it] = -1;
  }

  for( int it=0; it<MaxHits; it++){
    event.Bh2Seg[it] = -1;
    event.tBh2[it] = -9999.;
    event.deBh2[it] = -9999.;

    event.Bh1Seg[it] = -1;
    event.tBh1[it] = -9999.;
    event.deBh1[it] = -9999.;
  }

  for( int it=0; it<NumOfSegBFT; it++){
    event.bft_clsize[it] = -999;
    event.bft_ctime[it]  = -999.;
    event.bft_clpos[it]  = -999.;
    event.bft_bh1mth[it] = 0;
  }

  for(int i = 0; i<MaxHits; ++i){
    event.nhBcOut[i] = 0;
    event.chisqrBcOut[i] = -999.;
    event.x0BcOut[i] = -999.;
    event.y0BcOut[i] = -999.;
    event.u0BcOut[i] = -999.;
    event.v0BcOut[i] = -999.;

    event.p_2nd[i] = -999.;
    event.p_3rd[i] = -999.;
    event.delta_2nd[i] = -999.;
    event.delta_3rd[i] = -999.;

    event.xinK18[i] = -999.;
    event.yinK18[i] = -999.;
    event.uinK18[i] = -999.;
    event.vinK18[i] = -999.;

    event.xoutK18[i] = -999.;
    event.youtK18[i] = -999.;
    event.uoutK18[i] = -999.;
    event.voutK18[i] = -999.;

    event.nhK18[i]   = 0;
    event.chisqrK18[i] = -999.;
    event.xtgtK18[i] = -999.;
    event.ytgtK18[i] = -999.;
    event.utgtK18[i] = -999.;
    event.vtgtK18[i] = -999.;

    event.thetaK18[i] = -999.;
    event.phiK18[i] = -999.;
  }

  event.ntCFT  = 0;  
  //event.ntCFTEx  = 0;  
  for(int i = 0; i<MaxDepth; ++i){
    event.phi[i]  = -999.;
    event.theta[i]  = -999.;
    for(int j = 0; j<3; ++j){
      event.Pos[i][j]  = -999.;
      event.Dir[i][j]  = -999.;
    }

    //event.phiEx[i]  = -999.;
    //event.thetaEx[i]  = -999.;    
    //event.nhVtxEx[i]  = -1;    
    for( int p=0; p<NumOfPlaneCFT; ++p ){      
      event.seg[p][i]        = -999;
      event.seg_max[p][i]    = -999;
      event.dphi[p][i]       = -999.;
      event.phi_ini[p][i]    = -999.;
      event.phi_track[p][i]  = -999.;
      event.dz[p][i]       = -999.;
      event.z_ini[p][i]    = -999.;
      event.z_track[p][i]  = -999.;
      /*
      event.dphiEx[p][i]       = -999.;
      event.phiEx_track[p][i]  = -999.;
      event.dzEx[p][i]       = -999.;
      event.zEx_track[p][i]  = -999.;
      */
      event.adc[p][i] = -999.; event.adc_max[p][i] = -999.;
      event.mip[p][i] = -999.; event.mip_max[p][i] = -999.;
      event.dE[p][i]  = -999.; event.dE_max[p][i]  = -999.;

    }
    event.Total_dE[i]    = -999.; event.Total_dE_max[i]    = -999.;
    event.Total_dEphi[i] = -999.; event.Total_dEphi_max[i] = -999.;
    event.Total_dEuv[i]  = -999.; event.Total_dEuv_max[i]  = -999.;

    event.Total_mip[i]    = -999.; event.Total_mip_max[i]    = -999.;
    event.Total_E[i]    = -999.; event.Total_E_max[i]    = -999.;
    event.P[i]    = -999.; event.P_max[i]    = -999.;

    event.energybgot[i] = -999; 
    event.segBGOt[i]  = -1;
    event.segPiIDt[i] = -1;
  }

  event.nhBGO  = 0; 
  for( int i=0; i<NumOfSegBGO; ++i ){      
    event.segBGO[i]    = -1;
    event.segBGOt[i]   = -1;
    event.adcbgo[i]    = -999; 
    event.energybgo[i] = -999; 
    event.tdcbgo[i]    = -999;
  }

  event.vtx_2p_x  =  -999.;
  event.vtx_2p_y  =  -999.;
  event.vtx_2p_z  =  -999.;
  event.cdist_2p  =  -999.;
  event.theta_2p  =  -999.;


  for(int i = 0; i<MaxHits; ++i){
    event.vtx_pp_x[i]  =  -999.;
    event.vtx_pp_y[i]  =  -999.;
    event.vtx_pp_z[i]  =  -999.;
    event.cdist_pp[i]  =  -999.;
    event.missmass[i]  =  -999.;
    event.missmom[i]  =  -999.;
    event.theta_pp[i]  =  -999.;

    event.vtx_BeamCft_x[i]  =  -999.;
    event.vtx_BeamCft_y[i]  =  -999.;
    event.vtx_BeamCft_z[i]  =  -999.;
    event.cdist_BeamCft[i]  =  -999.;
    event.theta_BeamCft[i]  =  -999.;
  }


}

//______________________________________________________________________________
bool
EventPPScat::ProcessingEnd( void )
{
  tree->Fill();
  return true;
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventPPScat;
}

//______________________________________________________________________________
namespace
{
  const int    NBinDTSDC3 =  400;
  const double MinDTSDC3  = -100.;
  const double MaxDTSDC3  =  300.;
  const int    NBinDLSDC3 =  100;
  const double MinDLSDC3  = -5.0;
  const double MaxDLSDC3  = 15.0;
}
//______________________________________________________________________________
bool
ConfMan:: InitializeHistograms( void )
{
  HB1( 1, "Status", 30, 0., 30. );

  // Consume too many memories 
  // use only when the each channnel study is necessary 
#if 0
  for (int l=0; l<NumOfPlaneCFT; l++) {
    int NumOfSeg = NumOfSegCFT[l];
    for (int seg=0; seg<NumOfSeg; seg++ ) {
      int hid = ((l+1)*1000+seg)*10;

      char buf[100];
      sprintf(buf, "ADC Low Gain %d-%d (track)", l, seg);
      HB1(hid+1, buf, 4092, 0, 4092);
      sprintf(buf, "ADC Low Gain %d-%d (track 60<theta<120)", l, seg);
      HB1(hid+2, buf, 4092, 0, 4092);

      sprintf(buf, "MIP calib %d-%d (track)", l, seg);
      HB1(hid+3, buf, 1000, 0, 10);
      sprintf(buf, "MIP calib %d-%d (track 60<theta<120)", l, seg);
      HB1(hid+4, buf, 1000, 0, 10);

      sprintf(buf, "ADC Low Gain %d-%d (track 20<theta<30)", l, seg);
      HB1(hid+5, buf, 4092, 0, 4092);
      sprintf(buf, "MIP calib %d-%d (track 20<theta<30)", l, seg);
      HB1(hid+6, buf, 1000, 0, 10);

      sprintf(buf, "ADC Low Gain %d-%d (track 30<theta<40)", l, seg);
      HB1(hid+7, buf, 4092, 0, 4092);
      sprintf(buf, "MIP calib %d-%d (track 30<theta<40)", l, seg);
      HB1(hid+8, buf, 1000, 0, 10);

    }
  }
#endif

  ////////////////////////////////////////////
  //Tree
  HBTree( "kurama","tree of KuramaTracking" );
  tree->Branch("evnum",     &event.evnum,    "evnum/I");
  tree->Branch("trigpat",    event.trigpat,  Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag, Form("trigflag[%d]/I", NumOfSegTrig));

  //Hodoscope
  tree->Branch("nhBh2",   &event.nhBh2,   "nhBh2/I");
  tree->Branch("Bh2Seg",   event.Bh2Seg,  "Bh2Seg[nhBh2]/D");
  tree->Branch("tBh2",     event.tBh2,    "tBh2[nhBh2]/D");
  tree->Branch("deBh2",    event.deBh2,   "deBh2[nhBh2]/D");
  tree->Branch("time0",   &event.time0,   "time0/D");

  tree->Branch("nhBh1",   &event.nhBh1,   "nhBh1/I");
  tree->Branch("Bh1Seg",   event.Bh1Seg,  "Bh1Seg[nhBh1]/D");
  tree->Branch("tBh1",     event.tBh1,    "tBh1[nhBh1]/D");
  tree->Branch("deBh1",    event.deBh1,   "deBh1[nhBh1]/D");
  tree->Branch("btof",    &event.btof,    "btof/D");

  //BFT
  tree->Branch("bft_ncl",        &event.bft_ncl,    "bft_ncl/I");
  tree->Branch("bft_ncl_bh1mth", &event.bft_ncl_bh1mth, "bft_ncl_bh1mth/I");
  tree->Branch("bft_clsize",      event.bft_clsize, "bft_clsize[bft_ncl]/I");
  tree->Branch("bft_ctime",       event.bft_ctime,  "bft_ctime[bft_ncl]/D");
  tree->Branch("bft_clpos",       event.bft_clpos,  "bft_clpos[bft_ncl]/D");
  tree->Branch("bft_bh1mth",      event.bft_bh1mth, "bft_bh1mth[bft_ncl]/I");

  // BcOut
  tree->Branch("nlBcOut",   &event.nlBcOut,     "nlBcOut/I");
  tree->Branch("ntBcOut",   &event.ntBcOut,     "ntBcOut/I");
  tree->Branch("nhBcOut",    event.nhBcOut,     "nhBcOut[ntBcOut]/I");
  tree->Branch("chisqrBcOut",event.chisqrBcOut, "chisqrIn[ntBcOut]/D");
  tree->Branch("x0BcOut",    event.x0BcOut,     "x0BcOut[ntBcOut]/D");
  tree->Branch("y0BcOut",    event.y0BcOut,     "y0BcOut[ntBcOut]/D");
  tree->Branch("u0BcOut",    event.u0BcOut,     "u0BcOut[ntBcOut]/D");
  tree->Branch("v0BcOut",    event.v0BcOut,     "v0BcOut[ntBcOut]/D");

  // K18
  tree->Branch("ntK18",      &event.ntK18,     "ntK18/I");
  tree->Branch("nhK18",       event.nhK18,     "nhK18[ntK18]/I");
  tree->Branch("chisqrK18",   event.chisqrK18, "chisqrK18[ntK18]/D");
  tree->Branch("p_2nd",       event.p_2nd,     "p_2nd[ntK18]/D");
  tree->Branch("p_3rd",       event.p_3rd,     "p_3rd[ntK18]/D");
  tree->Branch("delta_2nd",   event.delta_2nd, "delta_2nd[ntK18]/D");
  tree->Branch("delta_3rd",   event.delta_3rd, "delta_3rd[ntK18]/D");

  tree->Branch("xinK18",    event.xinK18,   "xinK18[ntK18]/D");
  tree->Branch("yinK18",    event.yinK18,   "yinK18[ntK18]/D");
  tree->Branch("uinK18",    event.uinK18,   "uinK18[ntK18]/D");
  tree->Branch("vinK18",    event.vinK18,   "vinK18[ntK18]/D");

  tree->Branch("xoutK18",    event.xoutK18,   "xoutK18[ntK18]/D");
  tree->Branch("youtK18",    event.youtK18,   "youtK18[ntK18]/D");
  tree->Branch("uoutK18",    event.uoutK18,   "uoutK18[ntK18]/D");
  tree->Branch("voutK18",    event.voutK18,   "voutK18[ntK18]/D");

  tree->Branch("xtgtK18",    event.xtgtK18,   "xtgtK18[ntK18]/D");
  tree->Branch("ytgtK18",    event.ytgtK18,   "ytgtK18[ntK18]/D");
  tree->Branch("utgtK18",    event.utgtK18,   "utgtK18[ntK18]/D");
  tree->Branch("vtgtK18",    event.vtgtK18,   "vtgtK18[ntK18]/D");

  tree->Branch("thetaK18",   event.thetaK18,  "thetaK18[ntK18]/D");
  tree->Branch("phiK18",     event.phiK18,    "phiK18[ntK18]/D");

  tree->Branch("ntCFT",     &event.ntCFT,    "ntCFT/I");
  tree->Branch("theta",     event.theta,    "theta[ntCFT]/D");
  tree->Branch("phi",       event.phi,      "phi[ntCFT]/D");
  
  //tree->Branch("seg",      event.seg,     Form("seg[%d][%d]/I", NumOfPlaneCFT, MaxDepth ) );

  tree->Branch("seg",      event.seg,     Form("seg[%d][ntCFT]/I", NumOfPlaneCFT) );
 
  tree->Branch("seg_max",  event.seg_max, Form("seg_max[%d][%d]/I", NumOfPlaneCFT, MaxDepth ) );

  tree->Branch("dphi",     event.dphi,     Form("dphi[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("phi_ini",  event.phi_ini,  Form("phi_ini[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("phi_track",event.phi_track,Form("phi_track[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );

  tree->Branch("dz",      event.dz,      Form("dz[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("z_ini",   event.z_ini,   Form("z_ini[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("z_track", event.z_track, Form("z_track[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );

  tree->Branch("adc",     event.adc,      Form("adc[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("mip",     event.mip,      Form("mip[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("dE",      event.dE,       Form("dE[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("adc_max", event.adc_max,  Form("adc_max[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("mip_max", event.mip_max,  Form("mip_max[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("dE_max",  event.dE_max,   Form("dE_max[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );

  tree->Branch("Total_mip",     event.Total_mip,    "total_mip[ntCFT]/D");
  tree->Branch("Total_mip_max",    event.Total_mip_max,    "total_mip_max[ntCFT]/D");

  tree->Branch("Total_dE",    event.Total_dE,    "totaldE[ntCFT]/D");
  tree->Branch("Total_dEphi", event.Total_dEphi, "totaldEphi[ntCFT]/D");
  tree->Branch("Total_dEuv",  event.Total_dEuv,  "totaldEuv[ntCFT]/D");
  tree->Branch("Total_dE_max",    event.Total_dE_max,    "totaldE_max[ntCFT]/D");
  tree->Branch("Total_dEphi_max", event.Total_dEphi_max, "totaldEphi_max[ntCFT]/D");
  tree->Branch("Total_dEuv_max",  event.Total_dEuv_max,  "totaldEuv_max[ntCFT]/D");
  tree->Branch("Total_E",    event.Total_E,    "total_E[ntCFT]/D");
  tree->Branch("Total_E_max",    event.Total_E_max,    "total_E_max[ntCFT]/D");
  tree->Branch("P",    event.P,    "P[ntCFT]/D");
  tree->Branch("P_max",    event.P_max,    "P_max[ntCFT]/D");

  tree->Branch("energyBGOt", event.energybgot, "energybgot[ntCFT]/D");

  // CFTEx
  /*
  tree->Branch("ntCFTEx",     &event.ntCFTEx,    "ntCFTEx/I");  
  tree->Branch("thetaEx",     event.thetaEx,    "thetaEx[ntCFTEx]/D");
  tree->Branch("phiEx",       event.phiEx,      "phiEx[ntCFTEx]/D");
  tree->Branch("nhVtxEx",     event.nhVtxEx,    "nhVtxEx[ntCFTEx]/I");
  
  tree->Branch("dphiEx",     event.dphiEx,     Form("dphiEx[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("phiEx_track",event.phi_track,Form("phi_track[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("dzEx",      event.dzEx,      Form("dzEx[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("zEx_track", event.zEx_track, Form("zEx_track[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  */

  // BGO
  tree->Branch("nhBGO",     &event.nhBGO,    "nhBGO/I");
  tree->Branch("segBGO",    event.segBGO,    "segBGO[nhBGO]/I");
  tree->Branch("segBGOt",   event.segBGOt,   "segBGOt[ntCFT]/I");
  tree->Branch("adcBGO",    event.adcbgo,    "adcbgo[24]/D");
  tree->Branch("energyBGO", event.energybgo, "energybgo[24]/D");
  tree->Branch("tdcBGO",    event.tdcbgo,    "adcbgo[24]/D");

  // PiID
  tree->Branch("segPiIDt",   event.segPiIDt,   "segPiIDt[ntCFT]/I");

  tree->Branch("vtx_2p_x",   &event.vtx_2p_x, "vtx_2p_x/D");
  tree->Branch("vtx_2p_y",   &event.vtx_2p_y, "vtx_2p_y/D");
  tree->Branch("vtx_2p_z",   &event.vtx_2p_z, "vtx_2p_z/D");
  tree->Branch("cdist_2p",   &event.cdist_2p, "cdist_2p/D");
  tree->Branch("theta_2p",   &event.theta_2p, "theta_2p/D");



  tree->Branch("nPP",      &event.nPP,     "nPP/I");
  tree->Branch("vtx_pp_x",   event.vtx_pp_x, "vtx_pp_x[nPP]/D");
  tree->Branch("vtx_pp_y",   event.vtx_pp_y, "vtx_pp_y[nPP]/D");
  tree->Branch("vtx_pp_z",   event.vtx_pp_z, "vtx_pp_z[nPP]/D");
  tree->Branch("missmass",   event.missmass, "missmass[nPP]/D");
  tree->Branch("missmom",   event.missmom, "missmom[nPP]/D");
  tree->Branch("theta_pp",   event.theta_pp, "theta_pp[nPP]/D");
  tree->Branch("cdist_pp",   event.cdist_pp, "cdist_pp[nPP]/D");

  tree->Branch("vtx_BeamCft_x",   event.vtx_BeamCft_x, "vtx_BeamCft_x[ntCFT]/D");
  tree->Branch("vtx_BeamCft_y",   event.vtx_BeamCft_y, "vtx_BeamCft_y[ntCFT]/D");
  tree->Branch("vtx_BeamCft_z",   event.vtx_BeamCft_z, "vtx_BeamCft_z[ntCFT]/D");
  tree->Branch("theta_BeamCft",   event.theta_BeamCft, "theta_BeamCft[ntCFT]/D");
  tree->Branch("cdist_BeamCft",   event.cdist_BeamCft, "cdist_BeamCft[ntCFT]/D");

  HPrint();
  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")        &&
      InitializeParameter<DCDriftParamMan>("DCDRFT") &&
      InitializeParameter<DCTdcCalibMan>("DCTDC")    &&
      InitializeParameter<HodoParamMan>("HDPRM")     &&
      InitializeParameter<HodoPHCMan>("HDPHC")       &&
      InitializeParameter<BH2Filter>("BH2FLT")       &&
      //InitializeParameter<BH1Match>("BH1MTH")        &&
      InitializeParameter<K18TransMatrix>("K18TM")   &&
      InitializeParameter<UserParamMan>("USER")      &&
      InitializeParameter<CFTPedCorMan>("CFTPED")    &&
      InitializeParameter<BGOTemplateManager>("BGOTEMP") &&
      InitializeParameter<BGOCalibMan>("BGOCALIB")   );
}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}

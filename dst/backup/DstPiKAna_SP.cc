/**
 *  file: DstPiKAna.cc
 *  date: 2017.04.10
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include <filesystem_util.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "RootHelper.hh"
#include "DatabasePDG.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "DCRawHit.hh"
#include "KuramaLib.hh"
#include "K18TrackD2U.hh"
#include "FiberCluster.hh"
#include "NuclearMass.hh"
#include "MathTools.hh"
#include "MsTParamMan.hh"
#include "UserParamMan.hh"
#include "EnergyCorrection.hh"

#include "DstHelper.hh"

#define Sigmaonly 1 //0: all event filled to tree 1:only sigma event filled to tree (check sigmaflag)

namespace
{
  using namespace root;
  using namespace dst;
  const std::string& class_name("DstPiKAna");
  ConfMan&            gConf = ConfMan::GetInstance();
  const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
  // const MsTParamMan&  gMsT  = MsTParamMan::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const HodoPHCMan&   gPHC  = HodoPHCMan::GetInstance();
}

namespace dst
{
  enum kArgc
    {
      kProcess, kConfFile,
      kKuramaTracking, kK18Tracking, kHodoscope, kCFT,
      kOutFile, nArgc
    };
  std::vector<TString> ArgName =
    { "[Process]", "[ConfFile]", "[KuramaTracking]",
      "[K18Tracking]", "[Hodoscope]", "[CFT]",
      "[OutFile]" };
  std::vector<TString> TreeName =
    { "", "", "kurama", "k18track", "hodo", "cft",  "" };
  std::vector<TFile*> TFileCont;
  std::vector<TTree*> TTreeCont;

  const double pB_offset   = 1.000;
  const double pS_offset   = 1.000;
  const double pK18_offset = 0.000;
  const double pKURAMA_offset = 0.000;
  const double x_off = 0.000;
  const double y_off = 0.000;
  const double u_off = 0.000;
  const double v_off = 0.000;
  const double offsetCATCH =155;
}

//function for catchpid
double catchpidcurve(double param1, double param2, double x);

//functions for calculation
bool calcDecayPiMom(double p1, double m1, double m2, double m3, double cost, double *momCal);
bool calcDecayProtonMom(double p1, double m1, double m2, double m3, double cost, double *momCal1, double *momCal2);
bool calc2BodyKinema(double M1, double p1, double M2, double phi,
		     double *scatMomCal, double *scatEkinCal, double *scatThetaCM);
bool calc2BodyInelastic(double m1, double p1, double m2, double m3, double m4, 
			double theta, double *pCal1, double *pCal2, double *ThetaCM);
ThreeVector Xproduct(ThreeVector vec1, ThreeVector vec2);
bool calcBeamMomFrom2BodyKinema(double M1, double M2, double theta, double pScat,
				double *beamMomCal1, double *beamMomCal2);
bool calcBeamMomFromDecay(double M1, double M2, double M3, double p2, double theta, 
			  double *p1_1, double *p1_2);


//_____________________________________________________________________
struct Event
{
  int runnum;
  int evnum;
  int spill;

  //Trigger
  int trignhits;
  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];

  //Hodoscope
  int    nhBh1;
  int    csBh1[NumOfSegBH1];
  double Bh1Seg[NumOfSegBH1];
  double tBh1[NumOfSegBH1];
  double dtBh1[NumOfSegBH1];
  double deBh1[NumOfSegBH1];
  double btof[NumOfSegBH1];

  int    nhBh2;
  int    csBh2[NumOfSegBH2];
  double Bh2Seg[NumOfSegBH2];
  double tBh2[NumOfSegBH2];
  double t0Bh2[NumOfSegBH2];
  double dtBh2[NumOfSegBH2];
  double deBh2[NumOfSegBH2];

  int    nhTof;
  int csTof[NumOfSegTOF];
  double TofSeg[NumOfSegTOF];
  double tTof[NumOfSegTOF];
  double dtTof[NumOfSegTOF];
  double deTof[NumOfSegTOF];
  /*
  //Fiber
  int    nhBft;
  int    csBft[NumOfSegBFT];
  double tBft[NumOfSegBFT];
  double wBft[NumOfSegBFT];
  double BftPos[NumOfSegBFT];
  double BftSeg[NumOfSegBFT];
  int    nhSch;
  int    csSch[NumOfSegSCH];
  double tSch[NumOfSegSCH];
  double wSch[NumOfSegSCH];
  int NSch;
  double SchPos[NumOfSegSCH];
  int SchSeg[NumOfSegSCH];
  */
  //DC Beam
  int bft_ncl;
  int ntBcOut;
  int nlBcOut;
  int nhBcOut[MaxHits];
  double chisqrBcOut[MaxHits];
  double x0BcOut[MaxHits];
  double y0BcOut[MaxHits];
  double u0BcOut[MaxHits];
  double v0BcOut[MaxHits];
  double xtgtBcOut[MaxHits];
  double ytgtBcOut[MaxHits];
  double xbh2BcOut[MaxHits];
  double ybh2BcOut[MaxHits];

  int    ntK18;
  int    nhK18[MaxHits];
  double chisqrK18[MaxHits];
  double pK18[MaxHits];
  double xtgtK18[MaxHits];
  double ytgtK18[MaxHits];
  double utgtK18[MaxHits];
  double vtgtK18[MaxHits];
  double thetaK18[MaxHits];

  //DC KURAMA
  int ntSdcIn;
  int nlSdcIn;
  int nhSdcIn[MaxHits];
  double chisqrSdcIn[MaxHits];
  double x0SdcIn[MaxHits];
  double y0SdcIn[MaxHits];
  double u0SdcIn[MaxHits];
  double v0SdcIn[MaxHits];

  int ntSdcOut;
  int nlSdcOut;
  int nhSdcOut[MaxHits];
  double chisqrSdcOut[MaxHits];
  double u0SdcOut[MaxHits];
  double v0SdcOut[MaxHits];
  double x0SdcOut[MaxHits];
  double y0SdcOut[MaxHits];

  int    ntKurama;
  int    nhKurama[MaxHits];
  double chisqrKurama[MaxHits];
  double stof[MaxHits];
  double cstof[MaxHits];
  double path[MaxHits];
  double pKurama[MaxHits];
  double qKurama[MaxHits];
  double m2[MaxHits];
  double cm2[MaxHits];
  double xtgtKurama[MaxHits];
  double ytgtKurama[MaxHits];
  double utgtKurama[MaxHits];
  double vtgtKurama[MaxHits];
  double thetaKurama[MaxHits];
  double xtofKurama[MaxHits];
  double ytofKurama[MaxHits];
  double utofKurama[MaxHits];
  double vtofKurama[MaxHits];
  double tofsegKurama[MaxHits];

  //CFT
  int ntCFT;
  int ntProton;
  int ntOther;
  int tracknoproton[MaxHits];
  int tracknoother[MaxHits];
  double theta_cft[MaxHits];
  double PosT0_x[MaxHits];
  double PosT0_y[MaxHits];
  double PosT0_z[MaxHits];
  double Dir_x[MaxHits];
  double Dir_y[MaxHits];
  double Dir_z[MaxHits];
  double phi_cft[MaxHits];
  double vtx_cft[MaxHits];
  double vty_cft[MaxHits];
  double vtz_cft[MaxHits];
  double Total_dE[MaxHits];
  double Total_dEphi[MaxHits];
  double Total_dEuv[MaxHits];
  double Total_dE_max[MaxHits];
  double Total_dEphi_max[MaxHits];
  double Total_dEuv_max[MaxHits];
  int nhit_phi[MaxHits];
  int nhit_uv[MaxHits];
  int segBGOt[MaxHits];
  double energyBGOt[MaxHits];
  int segPiIDt[MaxHits];
  int protonflagt[MaxHits];
  int BGOnohitt[MaxHits];
  int simKuramat[MaxHits];

  //Reaction
  int    nPi;
  int    nK;
  int    nPiK;
  double vtx_KURAMA[MaxHits];
  double vty_KURAMA[MaxHits];
  double vtz_KURAMA[MaxHits];
  double ccm2[MaxHits];
  double chisqrKuramapik[MaxHits];
  double pKuramapik[MaxHits];
  double closedist_KURAMA[MaxHits];
  double theta[MaxHits];
  double MissMass[MaxHits];
  double MissMassCorr[MaxHits];
  double MissMassCorrDE[MaxHits];
  double MissMassPiPi[MaxHits];
  double MissMassPiP[MaxHits];
  double thetaCM[MaxHits];
  double costCM[MaxHits];

  double MissMom[MaxHits];
  double MissMomx[MaxHits];
  double MissMomy[MaxHits];
  double MissMomz[MaxHits];
  double MissMomcal[MaxHits];
  double MissMomxcal[MaxHits];
  double MissMomycal[MaxHits];
  double MissMomzcal[MaxHits];
  double SigmaBeta[MaxHits];

  double xpi[MaxHits];
  double ypi[MaxHits];
  double upi[MaxHits];
  double vpi[MaxHits];
  double xk[MaxHits];
  double yk[MaxHits];
  double uk[MaxHits];
  double vk[MaxHits];
  double uc[MaxHits];
  double vc[MaxHits];

  double pOrg[MaxHits];
  double pCalc[MaxHits];
  double pCorr[MaxHits];
  double pCorrDE[MaxHits];

  int KURAMAPID[MaxHits];//1:K+ 2:proton 3:pi+ 4:pi-

  //combination
  int nCatch;
  int nPiKCatch;
  double vtx_K18Catch[MaxHits];
  double vty_K18Catch[MaxHits];
  double vtz_K18Catch[MaxHits];
  double closedist_K18Catch[MaxHits];
  double theta_K18Catch[MaxHits];
  double vertex_distance[MaxHits];

  double vtx_XCatch[MaxHits];
  double vty_XCatch[MaxHits];
  double vtz_XCatch[MaxHits];
  double closedist_XCatch[MaxHits];
  double theta_XCatch[MaxHits];

  //best value of dE per events
  //NPScat assumption
  double DeltaE_NPScat;
  double ProtonMom_NPScat;
  double DecayNeutronMom;
  double DecayPionMom;
  double MissMassSigmaP_NPScat;
  double vtx_Decay2np;
  double vty_Decay2np;
  double vtz_Decay2np;
  double vtx_NPScat;
  double vty_NPScat;
  double vtz_NPScat;
  double cdistDecay2np;
  double cdistNPScat;
  double theta_NPScat;
  double thetaCM_NPScat;
  int ptrno_NPScat;
  int pitrno_NPScat;
  int pikno_NPScat;
  int priority_NPScat;

  //decayp assumption
  double DeltaE_DecayP;
  double ProtonMom_DecayP;
  double MissMassSigmaP_DecayP;
  double cdistDecayP;
  double SigmaMomcor_DecayP;
  double SigmaLength;
  double SigmaLengthcor;
  double theta_DecayP;
  int ptrno_DecayP;
  int pikno_DecayP;
  int priority_DecayP;

  //PiPScat assumption
  double DeltaP_PiPScat;
  double DecayNeutronMom2;
  double DecayPionMom2;
  double vtx_Decay2pip;
  double vty_Decay2pip;
  double vtz_Decay2pip;
  double vtx_PiPScat;
  double vty_PiPScat;
  double vtz_PiPScat;
  double cdistDecay2pip;
  double cdistPiPScat;
  int ptrno_PiPScat;
  int pitrno_PiPScat;
  int pikno_PiPScat;
  int priority_PiPScat;

  //PPScat assumption
  double DeltaP_PPScat;
  double theta_PPScat;
  double vtx_Decay2pp;
  double vty_Decay2pp;
  double vtz_Decay2pp;
  double vtx_PPScat;
  double vty_PPScat;
  double vtz_PPScat;
  double cdistDecay2pp;
  double cdistPPScat;
  int ptrno_PPScat;
  int p2trno_PPScat;
  int pikno_PPScat;
  int priority_PPScat;

  //sigma p scattering assumption
  double DeltaE_SigmaPScat;
  double MissMassSigmaP_SigmaPScat;
  double vtx_SigmaPScat;
  double vty_SigmaPScat;
  double vtz_SigmaPScat;
  double cdistSigmaPScat;
  int pikno_SigmaPScat;
  int ptrno_SigmaPScat;

  //npi+ decay mode
  double DeltaE_SigmaPScat2npi;
  double ProtonMom_SigmaPScat2npi;
  double MissMassSigmaP_SigmaPScat2npi;
  double vtx_ScatSigmaDecay2npi;
  double vty_ScatSigmaDecay2npi;
  double vtz_ScatSigmaDecay2npi;
  double cdistScatSigmaDecay2npi;
  double vdistance_SigmaPScat2npi;
  double theta_SigmaPScat2npi;
  double thetaCM_SigmaPScat2npi;
  int pikno_SigmaPScat2npi;
  int ptrno_SigmaPScat2npi;
  int pitrno_SigmaPScat2npi;
  int priority_SigmaPScat2npi;

  //ppi decay mode
  double DeltaE_SigmaPScat2ppi;
  double DeltaE2_SigmaPScat2ppi;
  double ProtonMom_SigmaPScat2ppi;
  double theta2proton_SigmaPScat2ppi;
  double MissMassSigmaP_SigmaPScat2ppi;
  double MissMassSigmaP2_SigmaPScat2ppi;
  double vtx_ScatSigmaDecay2ppi;
  double vty_ScatSigmaDecay2ppi;
  double vtz_ScatSigmaDecay2ppi;
  double cdistScatSigmaDecay2ppi;
  double vdistance_SigmaPScat2ppi;
  double theta_SigmaPScat2ppi;
  double thetaCM_SigmaPScat2ppi;
  int pikno_SigmaPScat2ppi;
  int ptrno_SigmaPScat2ppi;
  int p2trno_SigmaPScat2ppi;
  int priority_SigmaPScat2ppi;

};

//_____________________________________________________________________
struct Src
{
  int runnum;
  int evnum;
  int spill;

  //Trigger
  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];

  //Hodoscope
  int    nhBh1;
  int    csBh1[NumOfSegBH1];
  double Bh1Seg[NumOfSegBH1];
  double tBh1[NumOfSegBH1];
  double dtBh1[NumOfSegBH1];
  double deBh1[NumOfSegBH1];
  double btof[NumOfSegBH1];

  int    nhBh2;
  int    csBh2[NumOfSegBH2];
  double Bh2Seg[NumOfSegBH2];
  double tBh2[NumOfSegBH2];
  double t0Bh2[NumOfSegBH2];
  double dtBh2[NumOfSegBH2];
  double deBh2[NumOfSegBH2];

  double Time0Seg;
  double deTime0;
  double Time0;
  double CTime0;

  int    nhTof;
  int    csTof[NumOfSegTOF];
  double TofSeg[NumOfSegTOF];
  double tTof[NumOfSegTOF];
  double dtTof[NumOfSegTOF];
  double deTof[NumOfSegTOF];

  /*
  //Fiber
  int    nhBft;
  int    csBft[NumOfSegBFT];
  double tBft[NumOfSegBFT];
  double wBft[NumOfSegBFT];
  double BftPos[NumOfSegBFT];
  int    nhSch;
  int    NSch;
  int    SchSeg[NumOfSegSCH];
  int    csSch[NumOfSegSCH];
  double tSch[NumOfSegSCH];
  double wSch[NumOfSegSCH];
  double SchPos[NumOfSegSCH];
  int    nhFbh;
  */

  //DC Beam
  int bft_ncl;
  int ntBcOut;
  int nlBcOut;
  int nhBcOut[MaxHits];
  double chisqrBcOut[MaxHits];
  double x0BcOut[MaxHits];
  double y0BcOut[MaxHits];
  double u0BcOut[MaxHits];
  double v0BcOut[MaxHits];
  double xtgtBcOut[MaxHits];
  double ytgtBcOut[MaxHits];
  double xbh2BcOut[MaxHits];
  double ybh2BcOut[MaxHits];

  int    ntK18;
  int    nhK18[MaxHits];
  double chisqrK18[MaxHits];
  double pK18[MaxHits];
  double xtgtK18[MaxHits];
  double ytgtK18[MaxHits];
  double utgtK18[MaxHits];
  double vtgtK18[MaxHits];
  double thetaK18[MaxHits];

  //DC KURAMA
  int much;
  int ntSdcIn;
  int nlSdcIn;
  int nhSdcIn[MaxHits];
  double chisqrSdcIn[MaxHits];
  double x0SdcIn[MaxHits];
  double y0SdcIn[MaxHits];
  double u0SdcIn[MaxHits];
  double v0SdcIn[MaxHits];

  int ntSdcOut;
  int nlSdcOut;
  int nhSdcOut[MaxHits];
  double chisqrSdcOut[MaxHits];
  double u0SdcOut[MaxHits];
  double v0SdcOut[MaxHits];
  double x0SdcOut[MaxHits];
  double y0SdcOut[MaxHits];

  int    ntKurama;
  int    nhKurama[MaxHits];
  double chisqrKurama[MaxHits];
  double stof[MaxHits];
  double path[MaxHits];
  double pKurama[MaxHits];
  double qKurama[MaxHits];
  double m2[MaxHits];
  double xtgtKurama[MaxHits];
  double ytgtKurama[MaxHits];
  double utgtKurama[MaxHits];
  double vtgtKurama[MaxHits];
  double thetaKurama[MaxHits];
  double xtofKurama[MaxHits];
  double ytofKurama[MaxHits];
  double utofKurama[MaxHits];
  double vtofKurama[MaxHits];
  double tofsegKurama[MaxHits];

  //CATCH
  int ntCFT;
  double theta_cft[MaxHits];
  double phi_cft[MaxHits];
  double vtx_cft[MaxHits];
  double vty_cft[MaxHits];
  double vtz_cft[MaxHits];
  double Pos0_x[MaxHits];
  double Pos0_y[MaxHits];
  double Pos0_z[MaxHits];
  double Dir_x[MaxHits];
  double Dir_y[MaxHits];
  double Dir_z[MaxHits];
  double Total_dE[MaxHits];
  double Total_dEphi[MaxHits];
  double Total_dEuv[MaxHits];
  double Total_dE_max[MaxHits];
  double Total_dEphi_max[MaxHits];
  double Total_dEuv_max[MaxHits];
  int nhit_phi[MaxHits];
  int nhit_uv[MaxHits];
  int segBGOt[MaxHits];
  int segPiIDt[MaxHits];
  double energyBGO[NumOfSegBGO];

};

//_____________________________________________________________________
namespace root
{
  Event  event;
  Src    src;
  TH1   *h[MaxHist];
  TTree *tree;
}

//_____________________________________________________________________
int
main( int argc, char **argv )
{
  std::vector<std::string> arg( argv, argv+argc );

  if( !CheckArg( arg ) )
    return EXIT_FAILURE;
  if( !DstOpen( arg ) )
    return EXIT_FAILURE;
  if( !gConf.Initialize( arg[kConfFile] ) )
    return EXIT_FAILURE;

  int nevent = GetEntries( TTreeCont );

  CatchSignal::Set();

  int ievent = 0;
  for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
    InitializeEvent();
    if( DstRead( ievent ) ) tree->Fill();
  }

  std::cout << "#D Event Number: " << std::setw(6)
	    << ievent << std::endl;

  DstClose();

  return EXIT_SUCCESS;
}

//_____________________________________________________________________
bool
dst::InitializeEvent( void )
{
  event.runnum     = 0;
  event.evnum      = 0;
  event.spill      = 0;

  //Trigger
  for( int it=0; it<NumOfSegTrig; ++it ){
    event.trigpat[it]  = -1;
    event.trigflag[it] = -1;
  }

  //Hodoscope
  event.nhBh2  = 0;
  event.nhBh1  = 0;
  event.nhTof  = 0;

  for( int it=0; it<NumOfSegBH1; it++ ){
    event.csBh1[it]  = 0;
    event.Bh1Seg[it] = -1;
    event.tBh1[it]   = -9999.;
    event.dtBh1[it]  = -9999.;
    event.deBh1[it]  = -9999.;
    event.btof[it]   = -9999.;
  }
  for( int it=0; it<NumOfSegBH2; it++ ){
    event.csBh2[it]  = 0;
    event.Bh2Seg[it] = -1;
    event.tBh2[it]   = -9999.;
    event.t0Bh2[it]  = -9999.;
    event.dtBh2[it]  = -9999.;
    event.deBh2[it]  = -9999.;
  }
  for( int it=0; it<NumOfSegTOF; it++ ){
    event.csTof[it]  = 0;
    event.TofSeg[it] = -1;
    event.tTof[it]   = -9999.;
    event.dtTof[it]  = -9999.;
    event.deTof[it]  = -9999.;
  }
  /*
  //Fiber
  event.nhBft = 0;
  event.nhSch = 0;
  event.NSch = 0 ;
  for( int it=0; it<NumOfSegBFT; it++ ){
    event.csBft[it] = -999;
    event.tBft[it]  = -999.;
    event.wBft[it]   = -999.;
    event.BftPos[it]  = -999.;
  }
  for( int it=0; it<NumOfSegSCH; it++ ){
    event.csSch[it] = -999;
    event.tSch[it]  = -999.;
    event.wSch[it]   = -999.;
    event.SchPos[it]  = -999.;
    event.SchSeg[it] = -1;
  }
  */
  //DC
  event.bft_ncl  = 0;
  event.nlBcOut  = 0;
  event.nlSdcIn  = 0;
  event.nlSdcOut = 0;
  event.ntBcOut  = 0;
  event.ntSdcIn  = 0;
  event.ntSdcOut = 0;
  event.ntK18    = 0;
  event.ntKurama = 0;

  //Beam DC
  for( int it=0; it<MaxHits; it++){
    event.nhBcOut[it]     = 0;
    event.chisqrBcOut[it] = -1.0;
    event.x0BcOut[it] = -9999.0;
    event.y0BcOut[it] = -9999.0;
    event.u0BcOut[it] = -9999.0;
    event.v0BcOut[it] = -9999.0;

    event.xtgtBcOut[it] = -9999.0;
    event.ytgtBcOut[it] = -9999.0;
    event.xbh2BcOut[it] = -9999.0;
    event.ybh2BcOut[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhK18[it]     = 0;
    event.chisqrK18[it] = -1.0;
    event.xtgtK18[it] = -9999.;
    event.ytgtK18[it] = -9999.;
    event.utgtK18[it] = -9999.;
    event.vtgtK18[it] = -9999.;
    event.pK18[it]    = -9999.;
    event.thetaK18[it] = -9999.;
  }

  //KURAMA DC
  for( int it=0; it<MaxHits; it++){
    event.nhSdcIn[it]     = 0;
    event.chisqrSdcIn[it] = -1.;
    event.x0SdcIn[it] = -9999.;
    event.y0SdcIn[it] = -9999.;
    event.u0SdcIn[it] = -9999.;
    event.v0SdcIn[it] = -9999.;

    event.nhSdcOut[it]     = 0;
    event.chisqrSdcOut[it] = -1.;
    event.u0SdcOut[it] = -9999.;
    event.v0SdcOut[it] = -9999.;
    event.x0SdcOut[it] = -9999.;
    event.y0SdcOut[it] = -9999.;

    event.nhKurama[it]     = 0;
    event.chisqrKurama[it] = -1.;
    event.xtgtKurama[it]   = -9999.;
    event.ytgtKurama[it]   = -9999.;
    event.utgtKurama[it]   = -9999.;
    event.vtgtKurama[it]   = -9999.;
    event.pKurama[it]      = -9999.;
    event.qKurama[it]      = -9999.;
    event.stof[it]         = -9999.;
    event.cstof[it]        = -9999.;
    event.path[it]         = -9999.;
    event.m2[it]           = -9999.;
    event.cm2[it]          = -9999.;
    event.thetaKurama[it]  = -9999.;
    event.xtofKurama[it]   = -9999.;
    event.ytofKurama[it]   = -9999.;
    event.utofKurama[it]   = -9999.;
    event.vtofKurama[it]   = -9999.;
    event.tofsegKurama[it] = -9999.;
  }

  //CFT
  event.ntCFT=0;
  event.ntProton=0;
  event.ntOther=0;
  for(int it=0; it<MaxHits; ++it){
    event.tracknoproton[it]=-1;
    event.tracknoother[it]=-1;
    event.theta_cft[it] =-9999.;
    event.PosT0_x[it]    =-9999.;
    event.PosT0_y[it]    =-9999.;
    event.PosT0_z[it]    =-9999.;
    event.Dir_x[it]    =-9999.;
    event.Dir_y[it]    =-9999.;
    event.Dir_z[it]    =-9999.;
    event.phi_cft[it]   =-9999.;
    event.vtx_cft[it]   =-9999.;
    event.vty_cft[it]   =-9999.;
    event.vtz_cft[it]   =-9999.;
    event.Total_dE[it]  =-9999.;
    event.Total_dEphi[it]=-9999.;
    event.Total_dEuv[it]=-9999.;
    event.Total_dE_max[it]=-9999.;
    event.Total_dEphi_max[it]=-9999.;
    event.Total_dEuv_max[it]=-9999.;
    event.nhit_phi[it] =0;
    event.nhit_uv[it]  =0;
    event.segBGOt[it]  =-1;
    event.energyBGOt[it] =-9999.;
    event.segPiIDt[it] =-1;
    event.protonflagt[it]=0;
    event.BGOnohitt[it]=-1;
    event.simKuramat[it]=0;
    }

  //Reaction & combination
  event.nPi = 0;
  event.nK = 0;
  event.nPiK = 0;
  event.nCatch = 0;
  event.nPiKCatch = 0;

  for( int it=0; it<MaxHits; ++it ){
    event.vtx_KURAMA[it]       = -9999.;
    event.vty_KURAMA[it]       = -9999.;
    event.vtz_KURAMA[it]       = -9999.;
    event.ccm2[it]             = -9999.;
    event.chisqrKuramapik[it]  = -9999.;
    event.pKuramapik[it]       = -9999.;
    event.closedist_KURAMA[it] = -9999.;
    event.theta[it]     = -9999.;
    event.thetaCM[it]   = -9999.;
    event.costCM[it]    = -9999.;
    event.MissMass[it]  = -9999.;
    event.MissMassCorr[it]  = -9999.;
    event.MissMassCorrDE[it]  = -9999.;
    event.MissMassPiPi[it] =-9999.;
    event.MissMassPiP[it]=-9999.;

    event.MissMom[it]  =-9999.;
    event.MissMomx[it] =-9999.;
    event.MissMomy[it] =-9999.;
    event.MissMomz[it] =-9999.;

    event.xpi[it] = -9999.0;
    event.ypi[it] = -9999.0;
    event.upi[it] = -9999.0;
    event.vpi[it] = -9999.0;
    event.xk[it] = -9999.0;
    event.yk[it] = -9999.0;
    event.uk[it] = -9999.0;
    event.vk[it] = -9999.0;
    event.uc[it] = -9999.0;
    event.vc[it] = -9999.0;
    event.pOrg[it] = -9999.0;
    event.pCalc[it] = -9999.0;
    event.pCorr[it] = -9999.0;
    event.pCorrDE[it] = -9999.0;

    event.KURAMAPID[it] =0;
    event.MissMomcal[it]  =-9999.;
    event.MissMomxcal[it] =-9999.;
    event.MissMomycal[it] =-9999.;
    event.MissMomzcal[it] =-9999.;
    event.SigmaBeta[it] =-9999.;
    

    event.vtx_K18Catch[it] = -9999.0;
    event.vty_K18Catch[it] = -9999.0;
    event.vtz_K18Catch[it] = -9999.0;
    event.closedist_K18Catch[it] = -9999.0;
    event.theta_K18Catch[it] = -9999.0;
    event.vtx_XCatch[it] = -9999.0;
    event.vty_XCatch[it] = -9999.0;
    event.vtz_XCatch[it] = -9999.0;
    event.closedist_XCatch[it] =-9999.0;
    event.theta_XCatch[it]=-9999.0;
    event.vertex_distance[it]=-9999.0;
  }

  //np scattering assumption
  event.DeltaE_NPScat = -9999.;
  event.ProtonMom_NPScat=-9999.;
  event.MissMassSigmaP_NPScat=-9999.;
  event.vtx_Decay2np= -9999.;
  event.vty_Decay2np= -9999.;
  event.vtz_Decay2np= -9999.;
  event.vtx_NPScat= -9999.;
  event.vty_NPScat= -9999.;
  event.vtz_NPScat= -9999.;
  event.cdistDecay2np= -9999.;
  event.cdistNPScat= -9999.;
  event.DecayNeutronMom =-9999.;
  event.DecayPionMom =-9999.;
  event.theta_NPScat=-9999.;
  event.ptrno_NPScat =-1;
  event.pitrno_NPScat =-1;
  event.pikno_NPScat =-1;
  event.priority_NPScat =10;

  //decayp assumption
  event.DeltaE_DecayP =-9999.;
  event.ProtonMom_DecayP =-9999.;
  event.MissMassSigmaP_DecayP=-9999.;
  event.cdistDecayP =-9999.;
  event.SigmaMomcor_DecayP =-9999.;
  event.SigmaLength =-9999.;
  event.SigmaLengthcor=-9999.;
  event.theta_DecayP=-9999.;
  event.ptrno_DecayP=-1;
  event.pikno_DecayP=-1;
  event.priority_DecayP =10;

  //pip scattering assumption
  event.DeltaP_PiPScat = -9999.;
  event.vtx_Decay2pip= -9999.;
  event.vty_Decay2pip= -9999.;
  event.vtz_Decay2pip= -9999.;
  event.vtx_PiPScat= -9999.;
  event.vty_PiPScat= -9999.;
  event.vtz_PiPScat= -9999.;
  event.cdistDecay2pip= -9999.;
  event.cdistPiPScat= -9999.;
  event.DecayNeutronMom2 =-9999;
  event.DecayPionMom2 =-9999;
  event.ptrno_PiPScat =-1;
  event.pitrno_PiPScat =-1;
  event.pikno_PiPScat =-1;
  event.priority_PiPScat =10;

  //pp scattering assumption
  event.DeltaP_PPScat = -9999.;
  event.theta_PPScat =-9999.;
  event.vtx_Decay2pp= -9999.;
  event.vty_Decay2pp= -9999.;
  event.vtz_Decay2pp= -9999.;
  event.vtx_PPScat= -9999.;
  event.vty_PPScat= -9999.;
  event.vtz_PPScat= -9999.;
  event.cdistDecay2pp= -9999.;
  event.cdistPPScat= -9999.;
  event.ptrno_PPScat =-1;
  event.p2trno_PPScat =-1;
  event.pikno_PPScat =-1;
  event.priority_PPScat =10;

 //sigma p scattering assumption
  event.DeltaE_SigmaPScat = -9999.;
  event.MissMassSigmaP_SigmaPScat =-9999.;
  event.vtx_SigmaPScat =-9999.;
  event.vty_SigmaPScat =-9999.;
  event.vtz_SigmaPScat =-9999.;
  event.cdistSigmaPScat=-9999.;
  event.pikno_SigmaPScat=-1;
  event.ptrno_SigmaPScat=-1;

  //npi+ decay mode
  event.DeltaE_SigmaPScat2npi = -9999.;
  event.ProtonMom_SigmaPScat2npi =-9999.;
  event.MissMassSigmaP_SigmaPScat2npi=-9999.;
  event.vtx_ScatSigmaDecay2npi = -9999.;
  event.vty_ScatSigmaDecay2npi = -9999.;
  event.vtz_ScatSigmaDecay2npi = -9999.;
  event.cdistScatSigmaDecay2npi = -9999.;
  event.vdistance_SigmaPScat2npi = -9999.;
  event.theta_SigmaPScat2npi=-9999.;
  event.thetaCM_SigmaPScat2npi=-9999.;
  event.pikno_SigmaPScat2npi= -1;
  event.ptrno_SigmaPScat2npi= -1;
  event.pitrno_SigmaPScat2npi= -1;
  event.priority_SigmaPScat2npi = 10;

  //ppi decay mode
  event.DeltaE_SigmaPScat2ppi =-9999.;
  event.DeltaE2_SigmaPScat2ppi =-9999.;
  event.ProtonMom_SigmaPScat2ppi =-9999.;
  event.theta2proton_SigmaPScat2ppi =-9999.;
  event.MissMassSigmaP_SigmaPScat2ppi =-9999.;
  event.MissMassSigmaP2_SigmaPScat2ppi =-9999.;
  event.vtx_ScatSigmaDecay2ppi =-9999.;
  event.vty_ScatSigmaDecay2ppi =-9999.;
  event.vtz_ScatSigmaDecay2ppi =-9999.;
  event.cdistScatSigmaDecay2ppi =-9999.;
  event.vdistance_SigmaPScat2ppi =-9999.;
  event.theta_SigmaPScat2ppi=-9999.;
  event.thetaCM_SigmaPScat2ppi=-9999.;
  event.pikno_SigmaPScat2ppi = -1;
  event.ptrno_SigmaPScat2ppi = -1;
  event.p2trno_SigmaPScat2ppi = -1;
  event.priority_SigmaPScat2ppi = 10;

  return true;
}

//_____________________________________________________________________
bool
dst::DstOpen( std::vector<std::string> arg )
{
  int open_file = 0;
  int open_tree = 0;
  for( std::size_t i=0; i<nArgc; ++i ){
    if( i==kProcess || i==kConfFile || i==kOutFile ) continue;
    open_file += OpenFile( TFileCont[i], arg[i] );
    open_tree += OpenTree( TFileCont[i], TTreeCont[i], TreeName[i] );
  }

  if( open_file!=open_tree || open_file!=nArgc-3 )
    return false;
  if( !CheckEntries( TTreeCont ) )
    return false;

  TFileCont[kOutFile] = new TFile( arg[kOutFile].c_str(), "recreate" );

  return true;
}

//_____________________________________________________________________
bool
dst::DstRead( int ievent )
{
  static const std::string func_name("["+class_name+"::"+__func__+"]");

  static const double OffsetToF  = gUser.GetParameter("OffsetToF");
  static const double CatchPIDparam1  = gUser.GetParameter("CatchPID",0);
  static const double CatchPIDparam2  = gUser.GetParameter("CatchPID",1);

  static const double KaonMass    = pdg::KaonMass();
  static const double PionMass    = pdg::PionMass();
  static const double PiZeroMass  = 0.1349770;
  static const double ProtonMass  = pdg::ProtonMass();
  static const double NeutronMass = pdg::NeutronMass();
  static const double SigmaPMass  = pdg::SigmaPMass();
  

  if( ievent%10000==0 ){
    std::cout << "#D " << func_name << " Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  event.runnum   = src.runnum;
  event.evnum    = src.evnum;
  event.spill    = src.spill;

  event.bft_ncl  = src.bft_ncl;
  event.ntBcOut  = src.ntBcOut;
  event.ntSdcIn  = src.ntSdcIn;
  event.ntSdcOut = src.ntSdcOut;
  event.ntKurama = src.ntKurama;
  event.ntK18    = src.ntK18;
  event.ntCFT    = src.ntCFT;
  //event.nhBft    = src.nhBft;
  event.nhBh1    = src.nhBh1;
  event.nhBh2    = src.nhBh2;
  //event.nhSch    = src.nhSch;
  event.nhTof    = src.nhTof;
  //event.NSch     = src.NSch;

  const int ntBcOut  = event.ntBcOut;
  const int ntSdcIn  = event.ntSdcIn;
  const int ntSdcOut = event.ntSdcOut;
  const int ntKurama = event.ntKurama;
  const int ntK18    = event.ntK18;
  const int ntCFT    = event.ntCFT;
  //const int nhBft    = event.nhBft;
  const int nhBh1    = event.nhBh1;
  const int nhBh2    = event.nhBh2;
  //const int nhSch    = event.nhSch;
  const int nhTof    = event.nhTof;
  //const int NSch     = event.NSch;

#if 0
  std::cout << "#D DebugPrint" << std::endl
	    << " event  : " << std::setw(6) << ievent   << std::endl
	    << " BcOut  : " << std::setw(6) << ntBcOut  << std::endl
	    << " SdcIn  : " << std::setw(6) << ntSdcIn  << std::endl
	    << " SdcOut : " << std::setw(6) << ntSdcOut << std::endl
	    << " Kurama : " << std::setw(6) << ntKurama << std::endl
	    << " K18    : " << std::setw(6) << ntK18    << std::endl
            << " CFT    : " << std::setw(6) << ntCFT    << std::endl
	    << " BFT    : " << std::setw(6) << nhBft    << std::endl
	    << " BH1    : " << std::setw(6) << nhBh1    << std::endl
	    << " BH2    : " << std::setw(6) << nhBh2    << std::endl
	    << " SCH    : " << std::setw(6) << nhSch    << std::endl
	    << " TOF    : " << std::setw(6) << nhTof    << std::endl
	    << std::endl;
#endif

  // TFlag
  int trignhits =0;
  for(int i=0;i<NumOfSegTrig;++i){
    int tdc = src.trigflag[i];
    if( tdc>0 ){
    event.trigpat[trignhits]  = i + 1;
    event.trigflag[i] = tdc;
    trignhits++;
    }
  }
  event.trignhits=trignhits;

  HF1( 1, 0. );

  // if( event.ntKurama==0 ) return true;
  // HF1( 1, 1. );
  // if( event.ntK18==0 ) return true;
  // HF1( 1, 2. );
  // if( event.nhBh1==0 ) return true;
  // HF1( 1, 3. );
  // if( event.nhBh2==0 ) return true;
  // HF1( 1, 4. );
  // if( event.nhTof==0 ) return true;
  // HF1( 1, 5. );
  // if( event.nhSch==0 ) return true;
  // HF1( 1, 6. );
  // if( event.nhFbh==0 ) return true;
  // HF1( 1, 7. );

  std::vector <ThreeVector> PiPCont, PiXCont;
  std::vector <ThreeVector> KPCont,  KXCont;
  std::vector <ThreeVector> CatchDCont, CatchXCont;
  /*
  // BFT
  for( int i=0; i<nhBft; ++i ){
    event.csBft[i]  = src.csBft[i];
    event.tBft[i]   = src.tBft[i];
    event.wBft[i]   = src.wBft[i];
    event.BftPos[i] = src.BftPos[i];
  }
  */
  // BH1
  double btof=-9999.9;
  for( int i=0; i<nhBh1; ++i ){
    event.csBh1[i]  = src.csBh1[i];
    event.Bh1Seg[i] = src.Bh1Seg[i];
    event.tBh1[i]   = src.tBh1[i];
    event.dtBh1[i]  = src.dtBh1[i];
    event.deBh1[i]  = src.deBh1[i];
    event.btof[i]   = src.btof[i];
  }
  btof =event.btof[0];

  // BH2
  double time0 = src.Time0;
  for( int i=0; i<nhBh2; ++i ){
    event.csBh2[i]  = src.csBh2[i];
    event.Bh2Seg[i] = src.Bh2Seg[i];
    event.tBh2[i]   = src.tBh2[i];
    event.t0Bh2[i]  = src.t0Bh2[i];
    event.dtBh2[i]  = src.dtBh2[i];
    event.deBh2[i]  = src.deBh2[i];
  }

  /*
  // SCH
  for( int i=0; i<nhSch; ++i ){
    event.csSch[i]  = src.csSch[i];
    event.tSch[i]   = src.tSch[i];
    event.wSch[i]   = src.wSch[i];
    event.SchPos[i] = src.SchPos[i];
  }

  for( int i=0; i<NSch; ++i){
    event.SchSeg[i] = src.SchSeg[i];
  }
  */

  // TOF
  for( int i=0; i<nhTof; ++i ){
    event.csTof[i]  = src.csTof[i];
    event.TofSeg[i] = src.TofSeg[i];
    event.tTof[i]   = src.tTof[i];
    event.dtTof[i]  = src.dtTof[i];
    event.deTof[i]  = src.deTof[i];
  }

  ////////// BcOut
  event.nlBcOut = src.nlBcOut;
  for( int it=0; it<ntBcOut; ++it ){
    event.nhBcOut[it]     = src.nhBcOut[it];
    event.chisqrBcOut[it] = src.chisqrBcOut[it];
    event.x0BcOut[it]     = src.x0BcOut[it];
    event.y0BcOut[it]     = src.y0BcOut[it];
    event.u0BcOut[it]     = src.u0BcOut[it];
    event.v0BcOut[it]     = src.v0BcOut[it];
  }

  ////////// SdcIn
  event.nlSdcIn = src.nlSdcIn;
  for( int it=0; it<ntSdcIn; ++it ){
    event.nhSdcIn[it]     = src.nhSdcIn[it];
    event.chisqrSdcIn[it] = src.chisqrSdcIn[it];
    event.x0SdcIn[it]     = src.x0SdcIn[it];
    event.y0SdcIn[it]     = src.y0SdcIn[it];
    event.u0SdcIn[it]     = src.u0SdcIn[it];
    event.v0SdcIn[it]     = src.v0SdcIn[it];
  }

  ////////// SdcOut
  event.nlSdcOut = src.nlSdcOut;
  for( int it=0; it<ntSdcOut; ++it ){
    event.nhSdcOut[it]     = src.nhSdcOut[it];
    event.chisqrSdcOut[it] = src.chisqrSdcOut[it];
    event.x0SdcOut[it]     = src.x0SdcOut[it];
    event.y0SdcOut[it]     = src.y0SdcOut[it];
    event.u0SdcOut[it]     = src.u0SdcOut[it];
    event.v0SdcOut[it]     = src.v0SdcOut[it];
  }

  //catch

  for( int it=0; it<ntCFT; ++it){
    event.theta_cft[it]        = src.theta_cft[it];
    event.phi_cft[it]          = src.phi_cft[it];
    event.vtx_cft[it]          = src.vtx_cft[it];
    event.vty_cft[it]          = src.vty_cft[it];
    event.vtz_cft[it]          = src.vtz_cft[it];
    event.Total_dE[it]         = src.Total_dE[it];
    event.Total_dEphi[it]      = src.Total_dEphi[it];
    event.Total_dEuv[it]       = src.Total_dEuv[it];
    event.Total_dE_max[it]     = src.Total_dE_max[it];
    event.Total_dEphi_max[it]  = src.Total_dEphi_max[it];
    event.Total_dEuv_max[it]   = src.Total_dEuv_max[it];
    event.nhit_phi[it]         = src.nhit_phi[it];
    event.nhit_uv[it]          = src.nhit_uv[it];
    event.segBGOt[it]          = src.segBGOt[it];
    if(event.segBGOt[it]>-1){
    event.energyBGOt[it]       = src.energyBGO[src.segBGOt[it]];
    }
    event.segPiIDt[it]         = src.segPiIDt[it];
    double posx=src.Pos0_x[it];
    double posy=src.Pos0_y[it];
    double posz=src.Pos0_z[it];
    double dirx=src.Dir_x[it];
    double diry=src.Dir_y[it];
    double dirz=src.Dir_z[it];

    ThreeVector Pos_CFT(posx+(offsetCATCH-posz)*dirx/dirz,posy+(offsetCATCH-posz)*diry/dirz, 0.);
    ThreeVector Dir_CFT(dirx, diry, dirz);

    event.PosT0_x[it]           = Pos_CFT.x();
    event.PosT0_y[it]           = Pos_CFT.y();
    event.PosT0_z[it]           = Pos_CFT.z();
    event.Dir_x[it]            = src.Dir_x[it];
    event.Dir_y[it]            = src.Dir_y[it];
    event.Dir_z[it]            = src.Dir_z[it];

    CatchDCont.push_back(Dir_CFT);
    CatchXCont.push_back(Pos_CFT);
    
    //treat tracks without BGOhit
    if(event.energyBGOt[it]<0){
      event.energyBGOt[it]=0.0;
      event.BGOnohitt[it]=0;
      if((dirz>0 && pow(Pos_CFT.x()+(dirx/dirz)*190,2)+pow(Pos_CFT.y()+(diry/dirz)*190,2)<10000.0) ||(dirz<0 &&pow(Pos_CFT.x()+(dirx/dirz)*-210,2)+pow(Pos_CFT.y()+(diry/dirz)*-210,2)<10000.0)){
	event.BGOnohitt[it]=1;
      }
      if(event.BGOnohitt[it]==0 && event.Total_dEphi_max[it]*sin(event.theta_cft[it]*3.14/180)/event.nhit_phi[it]+event.Total_dEuv_max[it]/event.nhit_uv[it]>0.60){
	event.protonflagt[it]=1;//stop in CFT
      }
    }

    //CATCH PID
    if(event.energyBGOt[it]>0){
      if(event.Total_dEphi_max[it]*sin(event.theta_cft[it]*3.14/180)/event.nhit_phi[it]+event.Total_dEuv_max[it]/event.nhit_uv[it]>catchpidcurve(CatchPIDparam1,CatchPIDparam2,event.energyBGOt[it])){
	event.protonflagt[it]=1;
      }
    }

  }

  ////////// K+
  for( int itKurama=0; itKurama<ntKurama; ++itKurama ){
    int nh= src.nhKurama[itKurama];
    double chisqr = src.chisqrKurama[itKurama];
    double p = src.pKurama[itKurama];
    double x = src.xtgtKurama[itKurama];
    double y = src.ytgtKurama[itKurama];
    double u = src.utgtKurama[itKurama];
    double v = src.vtgtKurama[itKurama];
    double theta = src.thetaKurama[itKurama];
    double pt = p/std::sqrt(1.+u*u+v*v);
    ThreeVector Pos( x, y, 0. );
    ThreeVector Mom( pt*u, pt*v, pt );
    if( std::isnan( Pos.Mag() ) ) continue;
    //Calibration
    ThreeVector PosCorr( x+x_off, y+y_off, 0. );
    ThreeVector MomCorr( pt*(u+u_off), pt*(v+v_off), pt );
    double path = src.path[itKurama];
    double xt = PosCorr.x();
    double yt = PosCorr.y();
    double zt = PosCorr.z();
    double pCorr = MomCorr.Mag();
    double q  = src.qKurama[itKurama];
    double ut = xt/zt;
    double vt = yt/zt;

    // SdcOut vs TOF
    double TofSegKurama = src.tofsegKurama[itKurama];
    double stof  = -9999.;
    double cstof = -9999.;
    double m2    = -9999.;
    double cm2   = -9999.;
    // w/  TOF

    bool fl_matching = false;
    for( int j=0; j<nhTof; ++j ){
      double seg = src.TofSeg[j];
      if( seg == TofSegKurama ){
	// event.tTof[i]  = src.tTof[j];
	// event.dtTof[i] = src.dtTof[j];
	// event.deTof[i] = src.deTof[j];
	stof  = src.tTof[j] - time0 + OffsetToF;
	m2 = Kinematics::MassSquare( pCorr, path, stof );
	if(btof==-9999.9){
	  cstof=stof;
	} 
	else{
	  gPHC.DoStofCorrection( 8, 0, seg-1, 2, stof, btof, cstof );  
	  cm2 = Kinematics::MassSquare( pCorr, path, cstof );
	}
	fl_matching = true;
      }
    }
    if ( !fl_matching ) {
      for( int j=0; j<nhTof; ++j ){
	double seg = src.TofSeg[j];
	if ( seg == TofSegKurama ) { continue; }
	if( std::abs( seg - TofSegKurama ) < 2 ){
	  // event.tTof[i]  = src.tTof[j];
	  // event.dtTof[i] = src.dtTof[j];
	  // event.deTof[i] = src.deTof[j];
	  stof  = src.tTof[j] - time0 + OffsetToF;
	  m2 = Kinematics::MassSquare( pCorr, path, stof );
	  if(btof==-9999.9){
	    cstof=stof;
	  } 
	  else{
	    gPHC.DoStofCorrection( 8, 0, seg-1, 2, stof, btof, cstof ); 
	    cm2 = Kinematics::MassSquare( pCorr, path, cstof );
	  }
	}
      }
    }

    // w/o TOF
    // double minres = 1.0e10;
    // for( int j=0; j<nhTof; ++j ){
    //   double seg = src.TofSeg[j];
    //   double res = TofSegKurama - seg;
    //   if( std::abs( res ) < minres && std::abs( res ) < 1. ){
    // 	minres = res;
    // 	stof = src.tTof[j] - time0 + OffsetToF;
    //   }
    // }
    event.nhKurama[itKurama]   = nh;
    event.chisqrKurama[itKurama] = chisqr;
    event.pKurama[itKurama] = p;
    event.qKurama[itKurama] = q;
    event.xtgtKurama[itKurama] = x;
    event.ytgtKurama[itKurama] = y;
    event.utgtKurama[itKurama] = u;
    event.vtgtKurama[itKurama] = v;
    event.thetaKurama[itKurama] = theta;
    event.stof[itKurama] = stof;
    event.cstof[itKurama] = cstof;
    event.path[itKurama] = path;
    //event.m2[itKurama] =src.m2[itKurama];
    event.m2[itKurama] =m2;
    event.cm2[itKurama] = cm2;
    HF1( 3202, double(nh) );
    HF1( 3203, chisqr );
    HF1( 3204, xt ); HF1( 3205, yt );
    HF1( 3206, ut ); HF1( 3207, vt );
    HF2( 3208, xt, ut ); HF2( 3209, yt, vt );
    HF2( 3210, xt, yt );
    HF1( 3211, pCorr );
    HF1( 3212, path );

    //HF1( 3213, m2 );

    //     double xTof=(clTof->MeanSeg()-7.5)*70.;
    //     double yTof=(clTof->TimeDif())*800./12.;
    //     HF2( 4011, clTof->MeanSeg()-0.5, xtof );
    //     HF2( 4013, clTof->TimeDif(), ytof );
    //     HF1( 4015, xtof-xTof ); HF1( 4017, ytof-yTof );
    //     HF2( 4019, xtof-xTof, ytof-yTof );

    //       double ttof=clTof->CMeanTime()-time0,
    //       HF1( 322, clTof->ClusterSize() );
    //       HF1( 323, clTof->MeanSeg()-0.5 );
    //       HF1( 324, ttof );
    //       HF1( 325, clTof->DeltaE() );
    //       double u0in=trIn->GetU0();
    //       HF2( 4001, u0in, ttof ); HF2( 4003, u0in, ttof+12.5*u0in );
    HF1( 4202, double(nh) );
    HF1( 4203, chisqr );
    HF1( 4204, xt ); HF1( 4205, yt );
    HF1( 4206, ut ); HF1( 4207, vt );
    HF2( 4208, xt, ut ); HF2( 4209, yt, vt );
    HF2( 4210, xt, yt );
    HF1( 4211, pCorr );
    HF1( 4212, path );
    KXCont.push_back( PosCorr );
    KPCont.push_back( MomCorr );
  }



  ////////// pi
  for( int itK18=0; itK18<ntK18; ++itK18 ){
    int nh = src.nhK18[itK18];
    double chisqr = src.chisqrK18[itK18];
    //Calibration
    double p = src.pK18[itK18]*pB_offset+pK18_offset;
    double x = src.xtgtK18[itK18];
    double y = src.ytgtK18[itK18];
    double u = src.utgtK18[itK18];
    double v = src.vtgtK18[itK18];
    event.nhK18[itK18]     = nh;
    event.chisqrK18[itK18] = chisqr;
    event.pK18[itK18]      = p;
    event.xtgtK18[itK18]   = x;
    event.ytgtK18[itK18]   = y;
    event.utgtK18[itK18]   = u;
    event.vtgtK18[itK18]   = v;
    // double loss_bh2 = 1.09392e-3;
    // p = p - loss_bh2;
    double pt=p/std::sqrt(1.+u*u+v*v);
    ThreeVector Pos( x, y, 0. );
    ThreeVector Mom( pt*u, pt*v, pt );
    //double xo=trOut->GetX0(), yo=trOut->GetY0();
    HF1( 4104, p );
    HF1( 4105, x ); HF1( 4106, y );
    //HF1( 4107, xo ); HF1( 4108, yo ); HF1( 4109, u ); HF1( 4110, v );

    PiPCont.push_back(Mom); PiXCont.push_back(Pos);
  }

  //MissingMass and combination
  int nK  = KPCont.size();
  int nPi = PiPCont.size();
  int nCatch =CatchDCont.size();
  event.nPi = nPi;
  event.nK  = nK;
  event.nPiK = nPi*nK;
  event.nCatch =nCatch;
  if(nCatch==0){//to analyze events without CFT track
    nCatch=1;
    ThreeVector null(-9999.9,-9999.9,-9999.9);
    CatchDCont.push_back(null);
    CatchXCont.push_back(null);
  }
  event.nPiKCatch =event.nPiK*nCatch;
  bool sigmaflag =false;
  if( KPCont.size()==0 ) return false;
  HF1( 1, 8. );
  if( PiPCont.size()==0 ) return false;
  HF1( 1, 9. );
  HF1( 4101, double(nPi));
  HF1( 4201, double(nK) );
  int npik=0;
  int sigmamulti=0;
  for( int icatch=0; icatch<nCatch; ++icatch){
    ThreeVector dcatch =CatchDCont[icatch], xcatch =CatchXCont[icatch];
    for( int ipi=0; ipi<nPi; ++ipi ){
      ThreeVector ppi  = PiPCont[ipi], xpi = PiXCont[ipi];
      ThreeVector vertex_k18catch = Kinematics::VertexPoint( xcatch, xpi, dcatch, ppi);
      double cdist_k18catch =Kinematics::closeDist(xcatch,xpi,dcatch,ppi);
      double cost_k18catch =dcatch*ppi/(dcatch.Mag()*ppi.Mag());
      double theta_k18catch =std::acos(cost_k18catch)*math::Rad2Deg();

      for( int ikp=0; ikp<nK; ++ikp ){
	ThreeVector pkp = KPCont[ikp], xkp = KXCont[ikp];
	ThreeVector vertex_kurama = Kinematics::VertexPoint( xpi, xkp, ppi, pkp );
	// std::cout << "vertex : " << vert << " " << vert.Mag() << std::endl;
	double cdist_kurama = Kinematics::closeDist( xpi, xkp, ppi, pkp );

	ThreeVector diffvtx=vertex_kurama-vertex_k18catch;

	double us = pkp.x()/pkp.z(), vs = pkp.y()/pkp.z();
	double uc = dcatch.x()/dcatch.z(), vc=dcatch.y()/dcatch.z();
	double ub = ppi.x()/ppi.z(), vb = ppi.y()/ppi.z();
	double cost = ppi*pkp/(ppi.Mag()*pkp.Mag());

	double pk0   = pkp.Mag();
	double pCorr = pk0;

	ThreeVector pkpCorr( pCorr*pkp.x()/pkp.Mag(),
			     pCorr*pkp.y()/pkp.Mag(),
			     pCorr*pkp.z()/pkp.Mag() );

	//      ThreeVector ppiCorrDE = Kinematics::CorrElossIn( ppi, xpi, vert, PionMass );
	//      ThreeVector pkpCorrDE = Kinematics::CorrElossOut( pkpCorr, xkp, vert, KaonMass );

	LorentzVector LvPi( ppi, std::sqrt( PionMass*PionMass+ppi.Mag2() ) );
	//      LorentzVector LvPiCorrDE( ppiCorrDE, sqrt( PionMass*PionMass+ppiCorrDE.Mag2() ) );

	LorentzVector LvKp( pkp, std::sqrt( KaonMass*KaonMass+pkp.Mag2() ) );
	LorentzVector LvKpCorr( pkpCorr, std::sqrt( KaonMass*KaonMass+pkpCorr.Mag2() ) );
	//      LorentzVector LvKpCorrDE( pkpCorrDE, std::sqrt( KaonMass*KaonMass+pkpCorrDE.Mag2() ) );

	LorentzVector LvPiScat( pkp, std::sqrt( PionMass*PionMass+pkp.Mag2()));
	LorentzVector LvPScat( pkp, std::sqrt( ProtonMass*ProtonMass+pkp.Mag2()));

	LorentzVector LvC( 0., 0., 0., ProtonMass );
	LorentzVector LvCore( 0., 0., 0., 0. );

	LorentzVector LvRc       = LvPi+LvC-LvKp;
	LorentzVector LvRcCorr   = LvPi+LvC-LvKpCorr;
	//      LorentzVector LvRcCorrDE = LvPiCorrDE+LvC-LvKpCorrDE;

	LorentzVector LvPiPi     = LvPi+LvC-LvPiScat;
	LorentzVector LvPiP      = LvPi+LvC-LvPScat;
	double MisMass       = LvRc.Mag();//-LvC.Mag();
	double MisMass2      = LvRc.Mag2();
	double MisMassCorr   = LvRcCorr.Mag();//-LvC.Mag();
	//      double MisMassCorrDE = LvRcCorrDE.Mag();//-LvC.Mag();

	double MisMassPiPi   = LvPiPi.Mag();
	double MisMassPiP    = LvPiP.Mag();

	ThreeVector MisMom =LvRc.Vect();

	ThreeVector vertex_Xcatch =Kinematics::VertexPoint( vertex_kurama, xcatch, MisMom, dcatch);
	double cdist_Xcatch =Kinematics::closeDist(vertex_kurama, xcatch, MisMom,dcatch);
	double cost_Xcatch = MisMom*dcatch/(MisMom.Mag()*MisMom.Mag());
	double theta_Xcatch = std::acos(cost_Xcatch)*math::Rad2Deg();

	//Primary frame
	LorentzVector PrimaryLv = LvPi+LvC;
	double TotalEnergyCM = PrimaryLv.Mag();
	ThreeVector beta( 1/PrimaryLv.E()*PrimaryLv.Vect() );

	//CM
	double TotalMomCM
	  = 0.5*std::sqrt(( TotalEnergyCM*TotalEnergyCM
			    -( KaonMass+SigmaPMass )*( KaonMass+SigmaPMass ))
			  *( TotalEnergyCM*TotalEnergyCM
			     -( KaonMass-SigmaPMass )*( KaonMass-SigmaPMass )))/TotalEnergyCM;

	double costLab = cost;
	double cottLab = costLab/std::sqrt(1.-costLab*costLab);
	double bt=beta.Mag(), gamma=1./std::sqrt(1.-bt*bt);
	double gbep=gamma*bt*std::sqrt(TotalMomCM*TotalMomCM+KaonMass*KaonMass)/TotalMomCM;
	double a  = gamma*gamma+cottLab*cottLab;
	double bp = gamma*gbep;
	double c  = gbep*gbep-cottLab*cottLab;
	double dd = bp*bp-a*c;

	if( dd<0. ){
	  std::cerr << "dd<0." << std::endl;
	  dd = 0.;
	}

	double costCM = (std::sqrt(dd)-bp)/a;
	if( costCM>1. || costCM<-1. ){
	  std::cerr << "costCM>1. || costCM<-1." << std::endl;
	  costCM=-1.;
	}
	double sintCM  = std::sqrt(1.-costCM*costCM);
	double KaonMom = TotalMomCM*sintCM/std::sqrt(1.-costLab*costLab);

	double dpath =vertex_kurama.Mag();
	double cpath =(vertex_kurama.z()>0)?(event.path[ikp]-dpath):(event.path[ikp]+dpath);
	double dt    =Kinematics::CalcTimeOfFlight(event.pKurama[ikp],dpath,PionMass);
	double ccstof=(vertex_kurama.z()>0)?(event.cstof[ikp]-dt):(event.cstof[ikp]+dt);
	double m2= Kinematics::MassSquare( pCorr, cpath, ccstof);
	double chisqrKurama=src.chisqrKurama[ikp];
	double pkurama=src.pKurama[ikp];
        double qKurama=event.qKurama[ikp];


	if (npik<MaxHits) {
	  event.vtx_K18Catch[npik] = vertex_k18catch.x();
	  event.vty_K18Catch[npik] = vertex_k18catch.y();
	  event.vtz_K18Catch[npik] = vertex_k18catch.z();
	  event.closedist_K18Catch[npik] = cdist_k18catch;
	  event.theta_K18Catch[npik] =theta_k18catch;
	  event.vtx_KURAMA[npik]=vertex_kurama.x();
	  event.vty_KURAMA[npik]=vertex_kurama.y();
	  event.vtz_KURAMA[npik]=vertex_kurama.z();
	  event.closedist_KURAMA[npik] = cdist_kurama;
	  event.theta[npik]     = std::acos(cost)*math::Rad2Deg();
	  event.vertex_distance[npik] =diffvtx.Mag();
	  event.thetaCM[npik]   = std::acos(costCM)*math::Rad2Deg();
	  event.costCM[npik]    = costCM;

	  event.pKuramapik[npik]=pkurama;
	  event.ccm2[npik]=m2;
	  event.chisqrKuramapik[npik]=chisqrKurama;

	  event.MissMass[npik]       = MisMass;
	  event.MissMassCorr[npik]   = MisMassCorr;
	//	event.MissMassCorrDE[npik] = MisMassCorrDE;
	  event.MissMassCorrDE[npik] = 0;

	  event.MissMassPiPi[npik]   = MisMassPiPi;
	  event.MissMassPiP[npik]    = MisMassPiP;

	  event.MissMom[npik] =MisMom.Mag();
	  event.MissMomx[npik]=MisMom.x();
	  event.MissMomy[npik]=MisMom.y();
	  event.MissMomz[npik]=MisMom.z();

	  event.xk[npik] = xkp.x();
	  event.yk[npik] = xkp.y();
	  event.uc[npik] = uc;
	  event.vc[npik] = vc;
	  event.uk[npik] = us;
	  event.vk[npik] = vs;
	  if(pow(uc-us,2)+pow(vc-vs,2)<0.01){
	    event.simKuramat[icatch]=1;
	  }

	  event.xpi[npik] = xpi.x();
	  event.ypi[npik] = xpi.y();
	  event.upi[npik] = ub;
	  event.vpi[npik] = vb;
	  event.pOrg[npik] = pk0;
	  event.pCalc[npik] = KaonMom;
	  event.pCorr[npik] = pCorr;
	  //	event.pCorrDE[npik] = pkpCorrDE.Mag();
	  event.pCorrDE[npik] = 0;

	  event.vtx_XCatch[npik] =vertex_Xcatch.x();
	  event.vty_XCatch[npik] =vertex_Xcatch.y();
	  event.vtz_XCatch[npik] =vertex_Xcatch.z();
	  event.closedist_XCatch[npik] =cdist_Xcatch;
	  event.theta_XCatch[npik] =theta_Xcatch;
	  
	  if(chisqrKurama<50 && m2 <0.35 && 0.15 <m2 && pk0<1.0 && 0.8<pk0 && qKurama>0 && KaonMom>0){
	  event.KURAMAPID[npik] = 1; //sigma
	  event.MissMomxcal[npik] = ppi.x() - pkp.x()*KaonMom/pk0;
	  event.MissMomycal[npik] = ppi.y() - pkp.y()*KaonMom/pk0;
	  event.MissMomzcal[npik] = ppi.z() - pkp.z()*KaonMom/pk0;
	  ThreeVector MissMomV(event.MissMomxcal[npik],event.MissMomycal[npik],event.MissMomzcal[npik]);
	  event.MissMomcal[npik]=MissMomV.Mag();
	  event.SigmaBeta[npik]=MissMomV.Mag()/sqrt(SigmaPMass*SigmaPMass+MissMomV.Mag()*MissMomV.Mag());
	  if(npik<event.ntKurama){
	    HF1( 5002, MisMass2 );
	    if(MisMass<1.25 && MisMass>1.15){
	      HF1(5003, MisMass2 );
	    sigmamulti++;
	    }
	  }
	  sigmaflag = true;
	  }
	  if( 0.45<m2 && pk0<1.6 && qKurama>0){
	  event.KURAMAPID[npik] = 2; //proton
	  }
	  if( -0.1<m2 && 0.1> m2 && qKurama>0 && pk0<1.2){
	  event.KURAMAPID[npik] = 3; //scattering pi+
	  }
	  if( -0.1<m2 && m2<0.1 && qKurama<0){
	  event.KURAMAPID[npik] = 4; //pi-
	  }

	  npik++;
	} else {
	  std::cout << "#W npik: "<< npik << " exceeding MaxHits: " << MaxHits << std::endl;
	}

	HF1( 5001, vertex_kurama.z() );

	
	HF2( 5011, MisMass, us );
	HF2( 5012, MisMass, vs );
	HF2( 5013, MisMass, ub );
	HF2( 5014, MisMass, vb );
      }
    }
  }

  for( int it=0; it<ntCFT; ++it){
    if(event.simKuramat[it]==1) continue;
    if(event.protonflagt[it]==1 && event.segPiIDt[it]<0){
      event.tracknoproton[event.ntProton]=it;
      event.ntProton+=1;
    }
    else{
      if(event.Total_dE[it]>0){
	event.tracknoother[event.ntOther]=it;
	event.ntOther+=1;
      }
    }
  }

  //for light DST file (ignore not-sigma events)
#if Sigmaonly
  if(sigmaflag==false) return false;
#else
  if(sigmaflag==false) return true;
#endif

  //sigma decay(npi+) to npscattering assumption
  if(event.ntOther>0 && event.ntProton>0){
    for(int ipi=0;ipi<event.ntOther;ipi++){//pion analysis
      int tracknopi =event.tracknoother[ipi];
      if(event.simKuramat[tracknopi]==1) continue;
      ThreeVector PiDecayPosT0(event.PosT0_x[tracknopi], event.PosT0_y[tracknopi], event.PosT0_z[tracknopi]);
      ThreeVector PiDecayDir(event.Dir_x[tracknopi], event.Dir_y[tracknopi], event.Dir_z[tracknopi]);
      for(int ipik=0;ipik<event.nPiK;ipik++){//reaction analysis
	if(event.KURAMAPID[tracknopi*event.nPiK+ipik]==1){
	  ThreeVector SigmaPos(event.vtx_KURAMA[tracknopi*event.nPiK +ipik],event.vty_KURAMA[tracknopi*event.nPiK +ipik], event.vtz_KURAMA[tracknopi*event.nPiK +ipik]);
	  //ThreeVector SigmaMom(event.MissMomx[tracknopi*event.nPiK +ipik], event.MissMomy[tracknopi*event.nPiK +ipik], event.MissMomz[tracknopi*event.nPiK +ipik]);//momentum measured w/KURAMA
	  ThreeVector SigmaMom(event.MissMomxcal[tracknopi*event.nPiK +ipik], event.MissMomycal[tracknopi*event.nPiK +ipik], event.MissMomzcal[tracknopi*event.nPiK +ipik]);//momentum calculated w/KURAMA angle
	
	  double cdistDecayPi;
	  ThreeVector VertSigmaDecay =Kinematics::VertexPoint3D(SigmaPos, PiDecayPosT0, SigmaMom, PiDecayDir, cdistDecayPi);
	  double costDecayPi =SigmaMom*PiDecayDir/(SigmaMom.Mag()*PiDecayDir.Mag());
	  double thetaDecayPi =std::acos(costDecayPi)*math::Rad2Deg();

	  //sigma energy correction
	  ThreeVector SigmaLength =VertSigmaDecay - SigmaPos;
	  double SigmaMomcor, SigmaEcor;
	  int LH2Tag=0;

	  caldE(SigmaMom.Mag(), SigmaPMass, SigmaLength.Mag(), &SigmaMomcor, &SigmaEcor, LH2Tag);
	  //SigmaMom =SigmaMom*(SigmaMomcor/SigmaMom.Mag()); //sigma energy correction


	  double decayPiMomCal;
	  bool flagDecayCal = calcDecayPiMom( SigmaMom.Mag(), SigmaPMass, NeutronMass, PionMass, costDecayPi, &decayPiMomCal);
	  if(!flagDecayCal) continue;
	  ThreeVector PiDecayMomCal(decayPiMomCal*PiDecayDir.x()/PiDecayDir.Mag(),decayPiMomCal*PiDecayDir.y()/PiDecayDir.Mag(),decayPiMomCal*PiDecayDir.z()/PiDecayDir.Mag());
	  ThreeVector DecayNeutronMom =SigmaMom-PiDecayMomCal;
	  for(int ip=0;ip<event.ntProton;ip++){//proton analysis
	    int tracknop =event.tracknoproton[ip];
	    if(event.simKuramat[tracknop]==1) continue;
	    ThreeVector PScatPosT0(event.PosT0_x[tracknop], event.PosT0_y[tracknop], event.PosT0_z[tracknop]);
	    ThreeVector PScatDir(event.Dir_x[tracknop], event.Dir_y[tracknop], event.Dir_z[tracknop]);
	    
	    double cdistNPScat;
	    ThreeVector VertNPScat = Kinematics::VertexPoint3D( VertSigmaDecay, PScatPosT0, DecayNeutronMom, PScatDir, cdistNPScat);
	    double costNPScat =DecayNeutronMom*PScatDir/(DecayNeutronMom.Mag()*PScatDir.Mag());
	    double thetaNPScat=std::acos(costNPScat)*math::Rad2Deg();

	    if(thetaNPScat>0. && thetaNPScat<90.){
	      double scatNPMomCal;
	      double scatNPEkinCal;
	      double thetaNPScatCM;

	      bool flagNPScat =calc2BodyKinema(NeutronMass, DecayNeutronMom.Mag(), ProtonMass, thetaNPScat, &scatNPMomCal, &scatNPEkinCal, &thetaNPScatCM);

	      if(!flagNPScat) continue;
	      //proton energy correction
	      double ScatPEkinMeas =event.Total_dE[tracknop]+event.energyBGOt[tracknop];
	      double ScatPpMeas =sqrt(ScatPEkinMeas*ScatPEkinMeas+2.0*ScatPEkinMeas*ProtonMass*1000)/1000.;

	      double ScatPEtotcor;
	      double ScatPpcor;
	      CorrElossOutWithCFRP( &ScatPpcor, &ScatPEtotcor, ScatPpMeas, PROTON, PScatDir*(1.0/PScatDir.Mag()), VertNPScat, PScatPosT0);

	      //double DeltaE=ScatPEkinMeas-1000*scatNPEkinCal;//w/o proton dE correction
	      double DeltaE=((ScatPEtotcor-ProtonMass)-scatNPEkinCal)*1000.0;//w/ proton dE correction

	      LorentzVector LvProton(PScatDir*(ScatPpcor/PScatDir.Mag()),std::sqrt(ProtonMass*ProtonMass+pow(ScatPpcor,2)));
	      LorentzVector LvSigma(SigmaMom, std::sqrt(SigmaPMass*SigmaPMass+SigmaMom.Mag2()));
	      LorentzVector LvSigmaP=LvSigma-LvProton;

	      int priority =0;
	      if(cdistDecayPi>20) priority++;
	      if(cdistNPScat>20) priority++;
	      if(VertSigmaDecay.z()-event.vtz_KURAMA[ipik]<-10) priority++;
	      if(VertNPScat.z()-VertSigmaDecay.z()<-10) priority++;
	      //if(priority>event.priority_NPScat) continue;

	      if( pow(DeltaE,2)<pow(event.DeltaE_NPScat,2)){
		event.DeltaE_NPScat =DeltaE;
		event.ProtonMom_NPScat =ScatPpcor;
		event.DecayNeutronMom =DecayNeutronMom.Mag();
		event.DecayPionMom =PiDecayMomCal.Mag();
		event.MissMassSigmaP_NPScat=LvSigmaP.Mag2();
		event.vtx_Decay2np=VertSigmaDecay.x();
		event.vty_Decay2np=VertSigmaDecay.y();
		event.vtz_Decay2np=VertSigmaDecay.z();
		event.vtx_NPScat=VertNPScat.x();
		event.vty_NPScat=VertNPScat.y();
		event.vtz_NPScat=VertNPScat.z();
		event.cdistDecay2np=cdistDecayPi;
		event.cdistNPScat=cdistNPScat;
		event.theta_NPScat=thetaNPScat;
		event.thetaCM_NPScat=thetaNPScatCM;
		event.ptrno_NPScat=tracknop;
		event.pitrno_NPScat=tracknopi;
		event.pikno_NPScat=ipik;
		event.priority_NPScat=priority;
	      }

	    }
	  }


	}
      }							    

    }
  }
  //end of sigma decay(npi+) to npscattering assumption

  //1 proton from decay(ppi0) assumption
  if(event.ntProton>0){
    for( int ip=0; ip<event.ntProton;ip++){//proton loop
      int tracknop =event.tracknoproton[ip];
      if(event.simKuramat[tracknop]==1) continue;
      ThreeVector PDecayPosT0(event.PosT0_x[tracknop], event.PosT0_y[tracknop], event.PosT0_z[tracknop]);
      ThreeVector PDecayDir(event.Dir_x[tracknop], event.Dir_y[tracknop], event.Dir_z[tracknop]);
      double DecayPEkinMeas =event.Total_dE[tracknop]+event.energyBGOt[tracknop];
      double DecayPpMeas =sqrt(DecayPEkinMeas*DecayPEkinMeas+2.*DecayPEkinMeas*ProtonMass*1000)/1000.;
      for(int ipik=0;ipik<event.nPiK;ipik++){//reaction analysis
	if(event.KURAMAPID[tracknop*event.nPiK+ipik]!=1) continue;
	ThreeVector SigmaPos(event.vtx_KURAMA[tracknop*event.nPiK +ipik],event.vty_KURAMA[tracknop*event.nPiK +ipik], event.vtz_KURAMA[tracknop*event.nPiK +ipik]);
	//ThreeVector SigmaMom(event.MissMomx[tracknop*event.nPiK +ipik], event.MissMomy[tracknop*event.nPiK +ipik], event.MissMomz[tracknop*event.nPiK +ipik]);//momentum measured w/KURAMA
	ThreeVector SigmaMom(event.MissMomxcal[tracknop*event.nPiK +ipik], event.MissMomycal[tracknop*event.nPiK +ipik], event.MissMomzcal[tracknop*event.nPiK +ipik]);//momentum calculated w/KURAMA angle

	double cdistSigmaDecay;
	ThreeVector VertSigmaDecay = Kinematics::VertexPoint3D(SigmaPos, PDecayPosT0, SigmaMom, PDecayDir, cdistSigmaDecay);
	double costSigmaDecay=SigmaMom*PDecayDir/(SigmaMom.Mag()*PDecayDir.Mag());
	double thetaSigmaDecay=std::acos(costSigmaDecay)*math::Rad2Deg();

	//proton energy correction
	double DecayPEtotcor;
	double DecayPpcor;
	CorrElossOutWithCFRP(&DecayPpcor, &DecayPEtotcor, DecayPpMeas, PROTON, PDecayDir*(1.0/PDecayDir.Mag()), VertSigmaDecay, PDecayPosT0);

	//sigma energy correction
	ThreeVector SigmaLength =VertSigmaDecay - SigmaPos;
	double SigmaMomcor, SigmaEcor;
	int LH2Tag=0;

	caldE(SigmaMom.Mag(), SigmaPMass, SigmaLength.Mag(), &SigmaMomcor, &SigmaEcor, LH2Tag);
	
	double DecayPMomcal1= -999.;
	double DecayPMomcal2= -999.;
	double DeltaE1=-999.;
	double DeltaE2=-999.;
	double DeltaE=-999.;
       	bool flagSigmaDecay = calcDecayProtonMom( SigmaMom.Mag(), SigmaPMass, PiZeroMass, ProtonMass, costSigmaDecay, &DecayPMomcal1, &DecayPMomcal2);//w/o sigma dE correction
	//bool flagSigmaDecay = calcDecayProtonMom( SigmaMomcor, SigmaPMass, PiZeroMass, ProtonMass, costSigmaDecay, &DecayPMomcal1, &DecayPMomcal2); // w/ sigma dE correction

	if(!flagSigmaDecay) continue;
	double DecayPEkincal1 =(sqrt(ProtonMass*ProtonMass+DecayPMomcal1*DecayPMomcal1)-ProtonMass)*1000;
	double DecayPEkincal2 =(sqrt(ProtonMass*ProtonMass+DecayPMomcal2*DecayPMomcal2)-ProtonMass)*1000;
	//w/o proton dE correction
	//DeltaE1 =DecayPEkinMeas-DecayPEkincal1;
	//DeltaE2 =DecayPEkinMeas-DecayPEkincal2;
	
	//w/ proton dE correction
	DeltaE1 =(DecayPEtotcor-ProtonMass)*1000-DecayPEkincal1;
	DeltaE2 =(DecayPEtotcor-ProtonMass)*1000-DecayPEkincal2;

	if(pow(DeltaE1,2)<=pow(DeltaE2,2)){
	  DeltaE=DeltaE1;
	}
	if(pow(DeltaE2,2)<pow(DeltaE1,2)){
	  DeltaE=DeltaE2;
	}

	LorentzVector LvProton(PDecayDir*(DecayPpcor/PDecayDir.Mag()),std::sqrt(ProtonMass*ProtonMass+pow(DecayPpcor,2)));
	LorentzVector LvSigma(SigmaMom*(SigmaMomcor/SigmaMom.Mag()), std::sqrt(SigmaPMass*SigmaPMass+pow(SigmaMomcor,2)));
	LorentzVector LvSigmaP=LvSigma-LvProton;
	
	int priority=0;
	if(cdistSigmaDecay>20) priority++;
	if(VertSigmaDecay.z()-event.vtz_KURAMA[ipik]<-10) priority++;
	//if(priority>event.priority_DecayP) continue;

	if( pow(DeltaE,2)<pow(event.DeltaE_DecayP,2)){
	  event.DeltaE_DecayP =DeltaE;
	  event.ProtonMom_DecayP=DecayPEtotcor;
	  event.MissMassSigmaP_DecayP=LvSigmaP.Mag2();
	  event.SigmaMomcor_DecayP =SigmaMomcor;
	  event.cdistDecayP =cdistSigmaDecay;
	  event.SigmaLength =SigmaLength.Mag();
	  event.SigmaLengthcor =event.SigmaLength*sqrt(1-event.SigmaBeta[ipik]*event.SigmaBeta[ipik]);
	  event.theta_DecayP=thetaSigmaDecay;
	  event.ptrno_DecayP=tracknop;
	  event.pikno_DecayP=ipik;
	  event.priority_DecayP=priority;
	}


      }//reaction loop end
    }//proton loop end
  }
  //1 proton from decay(ppi0) assumption

  //end of 1 proton from decay(ppi0) assumption

  //sigma decay(npi+) to pip scattering assumption
  if(event.ntProton>0 && event.ntOther>0){
    for(int ipi=0; ipi<event.ntOther;ipi++){//pion loop
      int tracknopi =event.tracknoother[ipi];
      if(event.simKuramat[tracknopi]==1) continue;
      ThreeVector PiScatPosT0(event.PosT0_x[tracknopi], event.PosT0_y[tracknopi], event.PosT0_z[tracknopi]);
      ThreeVector PiScatDir(event.Dir_x[tracknopi], event.Dir_y[tracknopi], event.Dir_z[tracknopi]);
      for(int ip=0; ip<event.ntProton; ip++){// proton loop
	int tracknop =event.tracknoproton[ip];
	if(event.simKuramat[tracknop]==1) continue;
	ThreeVector PScatPosT0(event.PosT0_x[tracknop], event.PosT0_y[tracknop], event.PosT0_z[tracknop]);
	ThreeVector PScatDir(event.Dir_x[tracknop], event.Dir_y[tracknop], event.Dir_z[tracknop]);
	double cdistPiPScat;
	ThreeVector VertPiPScat = Kinematics::VertexPoint3D( PiScatPosT0, PScatPosT0, PiScatDir, PScatDir, cdistPiPScat);
	double costPiPScat =PiScatDir*PScatDir/(PiScatDir.Mag()*PScatDir.Mag());
	double thetaPiPScat =std::acos(costPiPScat)*math::Rad2Deg();

	double PScatKin = (event.Total_dE[tracknop]+event.energyBGOt[tracknop])*0.001;//GeV unit
	double PScatMom = sqrt(PScatKin*PScatKin +2.0*PScatKin*ProtonMass);

	//proton dE correction
	double ScatPEtotcor;
	double ScatPpcor;
	CorrElossOutWithCFRP( &ScatPpcor, &ScatPEtotcor, PScatMom, PROTON, PScatDir*(1.0/PScatDir.Mag()), VertPiPScat, PScatPosT0); 

	ThreeVector NormalVector =Xproduct(PScatDir,PiScatDir);

	for(int ipik=0;ipik<event.nPiK;ipik++){//reaction analysis
	  if(event.KURAMAPID[tracknopi*event.nPiK+ipik]!=1) continue;
	  ThreeVector SigmaPos(event.vtx_KURAMA[tracknopi*event.nPiK +ipik],event.vty_KURAMA[tracknopi*event.nPiK +ipik], event.vtz_KURAMA[tracknopi*event.nPiK +ipik]);
	  //ThreeVector SigmaMom(event.MissMomx[tracknopi*event.nPiK +ipik], event.MissMomy[tracknopi*event.nPiK +ipik], event.MissMomz[tracknopi*event.nPiK +ipik]);//momentum measured w/KURAMA
	  ThreeVector SigmaMom(event.MissMomxcal[tracknopi*event.nPiK +ipik], event.MissMomycal[tracknopi*event.nPiK +ipik], event.MissMomzcal[tracknopi*event.nPiK +ipik]);//momentum calculated w/KURAMA angle
	  double DeltaP1min=-999.0;
	  double DeltaP2min=-999.0;
	  for(double theta=0.; theta<=thetaPiPScat;theta+=0.1){//theta and deltap search end
	    ThreeVector PiDecayMom = PScatDir;
	    PiDecayMom.Rotate(theta*math::Deg2Rad(), NormalVector*(1/NormalVector.Mag()));

	    double PiDecayMomCal1, PiDecayMomCal2;
	    //bool flagPiMomCalc =calcBeamMomFrom2BodyKinema(PionMass, ProtonMass, theta, PScatMom, &PiDecayMomCal1, &PiDecayMomCal2); //w/o proton dE correction
	    bool flagPiMomCalc =calcBeamMomFrom2BodyKinema(PionMass, ProtonMass, theta, ScatPpcor, &PiDecayMomCal1, &PiDecayMomCal2); //w/ proton dE correction

	    if(!flagPiMomCalc) continue;

	    double cdistDecay2pip;
	    ThreeVector VertDecay2pip =Kinematics::VertexPoint3D(VertPiPScat, SigmaPos, PiDecayMom, SigmaMom, cdistDecay2pip);

	    //sigma energycorrection
	    ThreeVector SigmaLength =VertDecay2pip - SigmaPos;
	    double SigmaMomcor, SigmaEcor;
	    int LH2Tag=0;

	    caldE(SigmaMom.Mag(), SigmaPMass, SigmaLength.Mag(), &SigmaMomcor, &SigmaEcor, LH2Tag);
	    //SigmaMom =SigmaMom*(SigmaMomcor/SigmaMom.Mag());

	    double costPiDecay = PiDecayMom*SigmaMom/(PiDecayMom.Mag()*SigmaMom.Mag());
	    double PiDecayMomCalfromSigma;
	    bool flagPiMomCalc2 =calcDecayPiMom(SigmaMom.Mag(), SigmaPMass, NeutronMass, PionMass, costPiDecay, &PiDecayMomCalfromSigma);

	    if(!flagPiMomCalc2) continue;

	    //pion energy correction
	    ThreeVector DiffVert =VertPiPScat -VertDecay2pip;
	    double piMomcor, piEcor;
	    caldE(PiDecayMomCalfromSigma, PionMass, DiffVert.Mag(), &piMomcor, &piEcor, LH2Tag);
	
	    double DeltaP1 =-999.;
	    double DeltaP2 =-999.;

	     //w/o pion dE correction
	    if(PiDecayMomCal1>0) DeltaP1=PiDecayMomCal1-PiDecayMomCalfromSigma;
	    if(PiDecayMomCal2>0) DeltaP2=PiDecayMomCal2-PiDecayMomCalfromSigma;
	    

	    /*//w/ pion dE correction
	    if(PiDecayMomCal1>0) DeltaP1=PiDecayMomCal1-piMomcor;
	    if(PiDecayMomCal2>0) DeltaP2=PiDecayMomCal2-piMomcor;
	    */

	    double DeltaP=-999.;
	    double DeltaP1Scat=-999;
	    double DeltaP2Scat=-999;
	    double PiDecayMomCalc;
	    double ScatPiMomCal1, ScatPiMomCal2;
	    double ScatPiMomCalFromE1, ScatPiMomCalFromE2;
	    if(pow(DeltaP1,2)<pow(DeltaP1min,2)){
	      DeltaP1min =DeltaP1;
	      double thetaPiScat = thetaPiPScat-theta;
	      double mom1, mom2, thetaPiCM;
	      bool flagPiScatcalc= calc2BodyInelastic(PionMass, PiDecayMomCal1,ProtonMass,PionMass, ProtonMass, thetaPiScat, &mom1, &mom2, &thetaPiCM);
	      
	      double E1 =sqrt(PiDecayMomCal1*PiDecayMomCal1+PionMass*PionMass);
	      double E3 =ScatPEtotcor;
	      double EScatPi = E1 +ProtonMass -E3;
	      ScatPiMomCalFromE1 =sqrt(EScatPi*EScatPi-PionMass*PionMass);
	      ScatPiMomCal1 =mom1;
	      DeltaP1Scat=ScatPiMomCalFromE1-ScatPiMomCal1;

	    }
	    if(pow(DeltaP2,2)<pow(DeltaP2min,2)){
	      DeltaP2min =DeltaP2;
	      double thetaPiScat = thetaPiPScat-theta;
	      double mom1, mom2, thetaPiCM;
	      bool flagPiScatcalc= calc2BodyInelastic(PionMass, PiDecayMomCal2,ProtonMass,PionMass, ProtonMass, thetaPiScat, &mom1, &mom2, &thetaPiCM);
	      
	      double E1 =sqrt(PiDecayMomCal2*PiDecayMomCal2+PionMass*PionMass);
	      double E3 =ScatPEtotcor;
	      double EScatPi = E1 +ProtonMass -E3;
	      ScatPiMomCalFromE2 =sqrt(EScatPi*EScatPi-PionMass*PionMass);
	      ScatPiMomCal2 =mom1;
	      DeltaP2Scat=ScatPiMomCalFromE2-ScatPiMomCal2;

	    }
	    
	    if(pow(DeltaP1Scat,2)<=pow(DeltaP2Scat,2)){
	      DeltaP=DeltaP1Scat;
	      PiDecayMomCalc=PiDecayMomCal1;
	    }
	    if(pow(DeltaP1Scat,2)>pow(DeltaP2Scat,2)){
	      DeltaP=DeltaP2Scat;
	      PiDecayMomCalc=PiDecayMomCal2;
	    }	    

	    PiDecayMom =PiDecayMom*(PiDecayMomCalc/PiDecayMom.Mag());
	    ThreeVector NeutronMom =PiDecayMom -SigmaMom;

	    int priority =0;
	    if(cdistDecay2pip>20) priority++;
	    if(cdistPiPScat>20) priority++;
	    if(VertDecay2pip.z()-event.vtz_KURAMA[ipik]<-10) priority++;
	    if(VertPiPScat.z()-VertDecay2pip.z()<-10) priority++;
	    //if(priority>event.priority_PiPScat) continue;

	    if(pow(DeltaP,2)<pow(event.DeltaP_PiPScat,2)){
	      event.DeltaP_PiPScat =DeltaP;
	      event.DecayNeutronMom2 = NeutronMom.Mag();
	      event.DecayPionMom2 =PiDecayMom.Mag();
	      event.vtx_Decay2pip =VertDecay2pip.x();
	      event.vty_Decay2pip =VertDecay2pip.y();
	      event.vtz_Decay2pip =VertDecay2pip.z();
	      event.vtx_PiPScat =VertPiPScat.x();
	      event.vty_PiPScat =VertPiPScat.y();
	      event.vtz_PiPScat =VertPiPScat.z();
	      event.cdistDecay2pip =cdistDecay2pip;
	      event.cdistPiPScat =cdistPiPScat;
	      event.pikno_PiPScat =ipik;
	      event.ptrno_PiPScat =tracknop;
	      event.pitrno_PiPScat =tracknopi;
	      event.priority_PiPScat =priority;
	    }
	  
	  }//theta and deltaP search end
	}//reaction analysis end
      }// proton loop end
    }//pion loop end
  }
  //end of sigma decay(npi+) to pip scattering assumption

  //decay to sigma decay(ppi0) to pp scattering assumption
  if(event.ntProton>1){
    for(int ip=0; ip<event.ntProton;ip++){//proton loop1
      int tracknop =event.tracknoproton[ip];
      if(event.simKuramat[tracknop]==1) continue;
      ThreeVector PScat1PosT0(event.PosT0_x[tracknop], event.PosT0_y[tracknop], event.PosT0_z[tracknop]);
      ThreeVector PScat1Dir(event.Dir_x[tracknop], event.Dir_y[tracknop], event.Dir_z[tracknop]);
      double PScat1Ekin =event.Total_dE[tracknop]+event.energyBGOt[tracknop];
      double PScat1Mom =sqrt(PScat1Ekin*PScat1Ekin+2*PScat1Ekin*ProtonMass*1000.)/1000.;

      for(int ip2=0; ip2<event.ntProton;ip2++){//proton loop2
	if(ip2==ip) continue;
	int tracknop2 =event.tracknoproton[ip2];
	if(event.simKuramat[tracknop2]==1) continue;
	ThreeVector PScat2PosT0(event.PosT0_x[tracknop2], event.PosT0_y[tracknop2], event.PosT0_z[tracknop2]);
	ThreeVector PScat2Dir(event.Dir_x[tracknop2], event.Dir_y[tracknop2], event.Dir_z[tracknop2]);
	double PScat2Ekin =event.Total_dE[tracknop2]+event.energyBGOt[tracknop2];
	double PScat2Mom =sqrt(PScat1Ekin*PScat1Ekin+2*PScat1Ekin*ProtonMass*1000.)/1000.;


	double cdistPPScat;
	ThreeVector VertPPScat =Kinematics::VertexPoint3D(PScat1PosT0, PScat2PosT0, PScat1Dir, PScat2Dir, cdistPPScat);
	double costPPScat =PScat1Dir*PScat2Dir/(PScat1Dir.Mag()*PScat2Dir.Mag());
	double thetaPPScat =std::acos(costPPScat)*math::Rad2Deg();

	//proton energy correction
	double PScat1Etotcor;
	double PScat1pcor;
	double PScat2Etotcor;
	double PScat2pcor;
	CorrElossOutWithCFRP(&PScat1pcor, &PScat1Etotcor, PScat1Mom, PROTON, PScat1Dir*(1.0/PScat1Dir.Mag()), VertPPScat, PScat1PosT0);
	CorrElossOutWithCFRP(&PScat2pcor, &PScat2Etotcor, PScat2Mom, PROTON, PScat2Dir*(1.0/PScat2Dir.Mag()), VertPPScat, PScat2PosT0);

	//ThreeVector PScat1MomV =PScat1Dir*(PScat1Mom/PScat1Dir.Mag()); //w/o dE correction proton1
	//ThreeVector PScat2MomV =PScat2Dir*(PScat2Mom/PScat2Dir.Mag()); // w/o dE correction proton2

	ThreeVector PScat1MomV =PScat1Dir*(PScat1pcor/PScat1Dir.Mag()); // w/ dE correction proton1
	ThreeVector PScat2MomV =PScat2Dir*(PScat2pcor/PScat2Dir.Mag()); // w/ dE correction proton2

	ThreeVector DecayPMom =PScat1MomV+PScat2MomV;
	
	for(int ipik=0;ipik<event.nPiK;ipik++){//reaction analysis
	  if(event.KURAMAPID[tracknop*event.nPiK+ipik]!=1) continue;
	  ThreeVector SigmaPos(event.vtx_KURAMA[tracknop*event.nPiK +ipik],event.vty_KURAMA[tracknop*event.nPiK +ipik], event.vtz_KURAMA[tracknop*event.nPiK +ipik]);
	  //ThreeVector SigmaMom(event.MissMomx[tracknop*event.nPiK +ipik], event.MissMomy[tracknop*event.nPiK +ipik], event.MissMomz[tracknop*event.nPiK +ipik]);//momentum measured w/KURAMA
	  ThreeVector SigmaMom(event.MissMomxcal[tracknop*event.nPiK +ipik], event.MissMomycal[tracknop*event.nPiK +ipik], event.MissMomzcal[tracknop*event.nPiK +ipik]);//momentum calculated w/KURAMA angle

	  double cdistSigmaDecay2pp;
	  ThreeVector VertSigmaDecay2pp=Kinematics::VertexPoint3D(SigmaPos, VertPPScat, SigmaMom, DecayPMom, cdistSigmaDecay2pp);
	  double costDecayP =SigmaMom*DecayPMom/(SigmaMom.Mag()*DecayPMom.Mag());
	  double thetaDecayP =std::acos(costDecayP)*math::Rad2Deg();

	  //sigma energy correction
	  ThreeVector SigmaLength =VertSigmaDecay2pp -SigmaPos;
	  double SigmaMomcor, SigmaEcor;
	  int LH2Tag=0;

	  caldE(SigmaMom.Mag(), SigmaPMass, SigmaLength.Mag(), &SigmaMomcor, &SigmaEcor, LH2Tag);

	  double momCalDecayP = -999.;
	  double momCalDecayP2 =-999.;
	  double thetaCMDecayPPi0 =-999.;

	  bool flagDecay2pp =calc2BodyInelastic(SigmaPMass, SigmaMom.Mag(), 0, ProtonMass, PiZeroMass, thetaDecayP, &momCalDecayP, &momCalDecayP2, &thetaCMDecayPPi0); //w/o sigma dE correction
	  //bool flagDecay2pp =calc2BodyInelastic(SigmaPMass, SigmaMomcor, 0, ProtonMass, PiZeroMass, thetaDecayP, &momCalDecayP, &momCalDecayP2, &thetaCMDecayPPi0); //w/ sigma dE correction

	  //proton dE correction
	  ThreeVector DiffVert =VertPPScat -VertSigmaDecay2pp;
	  double DecayPMomcor, DecayPEcor;

	  if(!flagDecay2pp) continue;
	  double DeltaP1 =-999.;
	  double DeltaP2 =-999.;
	  if(momCalDecayP>0){
	    caldE(momCalDecayP, ProtonMass, DiffVert.Mag(), &DecayPMomcor, &DecayPEcor, LH2Tag);
	    DeltaP1 =DecayPMom.Mag()-momCalDecayP; //w/o proton correction
	    //DeltaP1 =DecayPMom.Mag()-DecayPMomcor; //w/ proton dE correction
	  }
	  if(momCalDecayP2>0){
	    caldE(momCalDecayP2, ProtonMass, DiffVert.Mag(), &DecayPMomcor, &DecayPEcor, LH2Tag);
	    DeltaP2=DecayPMom.Mag()-momCalDecayP2; //w/o proton correction
	    //DeltaP2 =DecayPMom.Mag()-DecayPMomcor; //w/ proton dE correction
	  }

	  double DeltaP =-999.;
	  if(pow(DeltaP1,2)<=pow(DeltaP2,2)) DeltaP=DeltaP1;
	  if(pow(DeltaP2,2)<pow(DeltaP1,2)) DeltaP=DeltaP2;

	  int priority =0;
	  if(cdistSigmaDecay2pp>20) priority++;
	  if(cdistPPScat>20) priority++;
	  if(VertSigmaDecay2pp.z()-event.vtz_KURAMA[ipik]<-10) priority++;
	  if(VertPPScat.z()-VertSigmaDecay2pp.z()<-10) priority++;
	  //if(priority>event.priority_PPScat) continue;


	  if(pow(DeltaP,2)<pow(event.DeltaP_PPScat,2)){
	    event.DeltaP_PPScat =DeltaP;
	    event.theta_PPScat =thetaPPScat;
	    event.vtx_Decay2pp =VertSigmaDecay2pp.x();
	    event.vty_Decay2pp =VertSigmaDecay2pp.y();
	    event.vtz_Decay2pp =VertSigmaDecay2pp.z();
	    event.vtx_PPScat =VertPPScat.x();
	    event.vty_PPScat =VertPPScat.y();
	    event.vtz_PPScat =VertPPScat.z();
	    event.cdistDecay2pp = cdistSigmaDecay2pp;
	    event.cdistPPScat =cdistPPScat;
	    event.ptrno_PPScat=tracknop;
	    event.p2trno_PPScat=tracknop2;
	    event.pikno_PPScat =ipik;
	    event.priority_PPScat =priority;
	  }


	}//reaction loop end
      }//proton loop2 end
    }//proton loop end
  }
  //end of sigma decay(ppi0) to pp scattering assumption

  //SigmaPscattering assumption

  if(event.ntProton>0){
    for( int ip=0; ip<event.ntProton;ip++){//proton loop
      int tracknop =event.tracknoproton[ip];
      if(event.simKuramat[tracknop]==1) continue;
      ThreeVector PScatPosT0(event.PosT0_x[tracknop], event.PosT0_y[tracknop], event.PosT0_z[tracknop]);
      ThreeVector PScatDir(event.Dir_x[tracknop], event.Dir_y[tracknop], event.Dir_z[tracknop]);
      for(int ipik=0;ipik<event.nPiK;ipik++){//reaction analysis
	if(event.KURAMAPID[tracknop*event.nPiK+ipik]!=1) continue;
	  ThreeVector SigmaPos(event.vtx_KURAMA[tracknop*event.nPiK +ipik],event.vty_KURAMA[tracknop*event.nPiK +ipik], event.vtz_KURAMA[tracknop*event.nPiK +ipik]);
	  //ThreeVector SigmaMom(event.MissMomx[tracknop*event.nPiK +ipik], event.MissMomy[tracknop*event.nPiK +ipik], event.MissMomz[tracknop*event.nPiK +ipik]);//momentum measured w/KURAMA
	  ThreeVector SigmaMom(event.MissMomxcal[tracknop*event.nPiK +ipik], event.MissMomycal[tracknop*event.nPiK +ipik], event.MissMomzcal[tracknop*event.nPiK +ipik]);//momentum calculated w/KURAMA angle

	  double cdistSigmaPScat;
	  ThreeVector VertSigmaPScat = Kinematics::VertexPoint3D(SigmaPos, PScatPosT0, SigmaMom, PScatDir, cdistSigmaPScat);
	  double costSigmaPScat=SigmaMom*PScatDir/(SigmaMom.Mag()*PScatDir.Mag());
	  double thetaSigmaPScat=std::acos(costSigmaPScat)*math::Rad2Deg();

	  //sigma energy correction
	  ThreeVector SigmaLength =VertSigmaPScat-SigmaPos;
	  double SigmaMomcor, SigmaEcor;
	  int LH2Tag=0;
	  caldE(SigmaMom.Mag(), SigmaPMass, SigmaLength.Mag(), &SigmaMomcor, &SigmaEcor, LH2Tag);
	  //SigmaMom=SigmaMom*(SigmaMomcor/SigmaMom.Mag());
	  
	  if(thetaSigmaPScat>0. && thetaSigmaPScat<90.){
	    double ScatPMomCal;
	    double ScatPEkinCal;
	    double thetaSigmaPScatCM;

	    bool flagSigmaPScat = calc2BodyKinema(SigmaPMass, SigmaMom.Mag(), ProtonMass, thetaSigmaPScat, &ScatPMomCal, &ScatPEkinCal, &thetaSigmaPScatCM); //w/o sigma energy correction

	    if(!flagSigmaPScat) continue;
	    double ScatPEkinMeas =event.Total_dE[tracknop]+event.energyBGOt[tracknop];
	    double ScatPMomMeas =sqrt(ScatPEkinMeas*ScatPEkinMeas+2*ScatPEkinMeas*ProtonMass*1000)/1000;

	    //proton energy correction
	    double ScatPEtotcor;
	    double ScatPpcor;
	    CorrElossOutWithCFRP( &ScatPpcor, &ScatPEtotcor, ScatPMomMeas, PROTON, PScatDir*(1.0/PScatDir.Mag()), VertSigmaPScat, PScatPosT0);


	    ThreeVector PScatMomCal=PScatDir*(ScatPMomCal/PScatDir.Mag());
	    ThreeVector SigmaScatMomCal=SigmaMom-PScatMomCal;
	    //ThreeVector PScatMomMeas=PScatDir*(ScatPMomMeas/PScatDir.Mag()); //w/o proton energy correction
	    ThreeVector PScatMomMeas=PScatDir*(ScatPpcor/PScatDir.Mag()); //w/ proton energy correction
	    ThreeVector SigmaScatMomMeas=SigmaMom-PScatMomMeas;
	    LorentzVector LvProton1(PScatMomMeas,std::sqrt(ProtonMass*ProtonMass+PScatMomMeas.Mag2()));
	    LorentzVector LvTarget(0.,0.,0.,ProtonMass);
	    LorentzVector LvSigma(SigmaMom, std::sqrt(SigmaPMass*SigmaPMass+SigmaMom.Mag2()));
	    //LorentzVector LvSigmaP1=LvSigma+LvTarget-LvProton1;
	    LorentzVector LvSigmaP1=LvSigma-LvProton1;//decay assumption MissMass
	    //double DeltaE=ScatPEkinMeas-1000*ScatPEkinCal; //w/o proton energy correction
	    double DeltaE=((ScatPEtotcor-ProtonMass)-ScatPEkinCal)*1000.0; //w/ proton energy correction
	    double MisMassSigmaP1=LvSigmaP1.Mag2();
	    if(pow(DeltaE,2)<pow(event.DeltaE_SigmaPScat,2)){
	      event.DeltaE_SigmaPScat = DeltaE;
	      event.MissMassSigmaP_SigmaPScat=MisMassSigmaP1;
	      event.vtx_SigmaPScat=VertSigmaPScat.x();
	      event.vty_SigmaPScat=VertSigmaPScat.y();
	      event.vtz_SigmaPScat=VertSigmaPScat.z();
	      event.cdistSigmaPScat=cdistSigmaPScat;
	      event.pikno_SigmaPScat=ipik;
	      event.ptrno_SigmaPScat=tracknop;
	    }
	    
	    if(event.ntOther>0){//npi+ Decay mode
	      for(int ipi=0;ipi<event.ntOther;ipi++){//pion loop
		int tracknopi =event.tracknoother[ipi];
		if(event.simKuramat[tracknopi]==1) continue;
		ThreeVector PiDecayPosT0(event.PosT0_x[tracknopi], event.PosT0_y[tracknopi], event.PosT0_z[tracknopi]);
		ThreeVector PiDecayDir(event.Dir_x[tracknopi], event.Dir_y[tracknopi], event.Dir_z[tracknopi]);
		double cdistScatSigmaDecay2npi;
		ThreeVector VertScatSigmaDecay2npi =Kinematics::VertexPoint3D(VertSigmaPScat, PiDecayPosT0, SigmaScatMomCal, PiDecayDir, cdistScatSigmaDecay2npi);
		ThreeVector VertDiff2npi=VertScatSigmaDecay2npi-VertSigmaPScat;

		double costDecayPi=SigmaScatMomCal*PiDecayDir/(SigmaScatMomCal.Mag()*PiDecayDir.Mag());
		double thetaDecayPi=std::acos(costDecayPi)*math::Rad2Deg();

		double decayPiMomCal;
		bool flagDecayCal = calcDecayPiMom(SigmaScatMomCal.Mag(), SigmaPMass, NeutronMass, PionMass, costDecayPi, &decayPiMomCal);
		if(!flagDecayCal) continue;

		int priority =0;
		if(cdistSigmaPScat>20) priority++;
		if(cdistScatSigmaDecay2npi>20) priority++;
		if(VertSigmaPScat.z()-event.vtz_KURAMA[ipik]<-10) priority++;
		if(VertScatSigmaDecay2npi.z()-VertSigmaPScat.z()<-10) priority++;
		if(priority>event.priority_SigmaPScat2npi) continue;

		if(pow(DeltaE,2)<pow(event.DeltaE_SigmaPScat2npi,2)){
		  event.DeltaE_SigmaPScat2npi =DeltaE;
		  event.ProtonMom_SigmaPScat2npi =ScatPpcor;
		  event.MissMassSigmaP_SigmaPScat2npi=MisMassSigmaP1;
		  event.vtx_ScatSigmaDecay2npi=VertScatSigmaDecay2npi.x();
		  event.vty_ScatSigmaDecay2npi=VertScatSigmaDecay2npi.y();
		  event.vtz_ScatSigmaDecay2npi=VertScatSigmaDecay2npi.z();
		  event.cdistScatSigmaDecay2npi=cdistScatSigmaDecay2npi;
		  event.vdistance_SigmaPScat2npi=VertDiff2npi.Mag();
		  event.theta_SigmaPScat2npi=thetaSigmaPScat;
		  event.thetaCM_SigmaPScat2npi=thetaSigmaPScatCM;
		  event.pikno_SigmaPScat2npi=ipik;
		  event.ptrno_SigmaPScat2npi=tracknop;
		  event.pitrno_SigmaPScat2npi=tracknopi;
		  event.priority_SigmaPScat2npi=priority;
		}

	      }//pion loop end
	    }//pion condition

	    if(event.ntProton>1){//p pi0 Decay mode
	      for(int ip2=0; ip2<event.ntProton;ip2++){//proton2 loop
		if(ip2==ip) continue;
		int tracknop2 =event.tracknoproton[ip2];
		if(event.simKuramat[tracknop2]==1) continue;
		ThreeVector PDecayPosT0(event.PosT0_x[tracknop2], event.PosT0_y[tracknop2], event.PosT0_z[tracknop2]);
		ThreeVector PDecayDir(event.Dir_x[tracknop2], event.Dir_y[tracknop2], event.Dir_z[tracknop2]);
		double cdistScatSigmaDecay2ppi;
		ThreeVector VertScatSigmaDecay2ppi =Kinematics::VertexPoint3D(VertSigmaPScat, PDecayPosT0, SigmaScatMomCal, PDecayDir, cdistScatSigmaDecay2ppi);
		ThreeVector VertDiff2ppi=VertScatSigmaDecay2ppi-VertSigmaPScat;

		//sigma energy correction 2
		double SigmaMomcor2, SigmaEcor2;

		caldE(SigmaScatMomCal.Mag(), SigmaPMass, VertDiff2ppi.Mag(), &SigmaMomcor2, &SigmaEcor2, LH2Tag);
		//SigmaScatMomCal =SigmaScatMomCal*(SigmaMomcor2/SigmaScatMomCal.Mag());

		double costDecayP=SigmaScatMomCal*PDecayDir/(SigmaScatMomCal.Mag()*PDecayDir.Mag());
		double thetaDecayP=std::acos(costDecayP)*math::Rad2Deg();

		double cost2proton =PScatDir*PDecayDir/(PScatDir.Mag()*PDecayDir.Mag());
		double theta2proton =std::acos(cost2proton)*math::Rad2Deg();

		double DecayPEkinMeas =event.Total_dE[tracknop2]+event.energyBGOt[tracknop2];
		double DecayPMomMeas =sqrt(DecayPEkinMeas*DecayPEkinMeas+2*DecayPEkinMeas*ProtonMass*1000)/1000;
		//proton energy correction
		double DecayPEtotcor;
		double DecayPpcor;
		CorrElossOutWithCFRP(&DecayPpcor, &DecayPEtotcor, DecayPMomMeas, PROTON, PDecayDir*(1.0/PDecayDir.Mag()), VertScatSigmaDecay2ppi, PDecayPosT0);

		ThreeVector PDecayMomMeas=PDecayDir*(DecayPMomMeas/PDecayDir.Mag());
		LorentzVector LvProton2(PDecayMomMeas, std::sqrt(ProtonMass*ProtonMass+PDecayMomMeas.Mag2()));
		LorentzVector LvSigmaScat(SigmaScatMomCal,std::sqrt(SigmaPMass*SigmaPMass+SigmaScatMomCal.Mag2()));
		LorentzVector LvSigmaP2=LvSigmaScat-LvProton2;
		double MisMassSigmaP2=LvSigmaP2.Mag2();

		double momCalDecayP =-999.;
		double momCalDecayP2 =-999.;
		double thetaCMDecayPPi0 =-999.;

		bool flagDecayPPi0 =calc2BodyInelastic(SigmaPMass, SigmaScatMomCal.Mag(),0, ProtonMass, PiZeroMass, thetaDecayP, &momCalDecayP, &momCalDecayP2, &thetaCMDecayPPi0);

		if(!flagDecayPPi0) continue;
		double DeltaE2_1=-999.;
		double DeltaE2_2=-999.;

		if(momCalDecayP>0){
		  //DeltaE2_1=DecayPEkinMeas-(sqrt(ProtonMass*ProtonMass+momCalDecayP*momCalDecayP)-ProtonMass)*1000; //w/o proton dE correction
		  DeltaE2_1=((DecayPEtotcor-ProtonMass)-(sqrt(ProtonMass*ProtonMass+momCalDecayP*momCalDecayP)-ProtonMass))*1000.0; //w/ proton dE correction
		}
		if(momCalDecayP2>0){
		  //DeltaE2_2=DecayPEkinMeas-(sqrt(ProtonMass*ProtonMass+momCalDecayP2*momCalDecayP2)-ProtonMass)*1000;
		  DeltaE2_2=((DecayPEtotcor-ProtonMass)-(sqrt(ProtonMass*ProtonMass+momCalDecayP2*momCalDecayP2)-ProtonMass))*1000.0; //w/ proton dE correction
		}
		double DeltaE2 =-999.;
		if(pow(DeltaE2_1,2)<=pow(DeltaE2_2,2)) DeltaE2=DeltaE2_1;
		if(pow(DeltaE2_2,2)<pow(DeltaE2_1,2)) DeltaE2=DeltaE2_2;

		int priority =0;
		if(cdistSigmaPScat>20) priority++;
		if(cdistScatSigmaDecay2ppi>20) priority++;
		if(VertSigmaPScat.z()-event.vtz_KURAMA[ipik]<-10) priority++;
		if(VertScatSigmaDecay2ppi.z()-VertSigmaPScat.z()<-10) priority++;
		//if(priority>event.priority_SigmaPScat2ppi) continue;

		if(pow(DeltaE,2)<pow(event.DeltaE_SigmaPScat2ppi,2)){
		  event.DeltaE_SigmaPScat2ppi =DeltaE;
		  event.DeltaE2_SigmaPScat2ppi =DeltaE2;
		  event.ProtonMom_SigmaPScat2ppi =ScatPpcor;
		  event.theta2proton_SigmaPScat2ppi =theta2proton;
		  event.MissMassSigmaP_SigmaPScat2ppi=MisMassSigmaP1;
		  event.MissMassSigmaP2_SigmaPScat2ppi=MisMassSigmaP2;
		  event.vtx_ScatSigmaDecay2ppi=VertScatSigmaDecay2ppi.x();
		  event.vty_ScatSigmaDecay2ppi=VertScatSigmaDecay2ppi.y();
		  event.vtz_ScatSigmaDecay2ppi=VertScatSigmaDecay2ppi.z();
		  event.cdistScatSigmaDecay2ppi=cdistScatSigmaDecay2ppi;
		  event.vdistance_SigmaPScat2ppi=VertDiff2ppi.Mag();
		  event.theta_SigmaPScat2ppi=thetaSigmaPScat;
		  event.thetaCM_SigmaPScat2ppi=thetaSigmaPScatCM;
		  event.pikno_SigmaPScat2ppi=ipik;
		  event.ptrno_SigmaPScat2ppi=tracknop;
		  event.p2trno_SigmaPScat2ppi=tracknop2;
		  event.priority_SigmaPScat2ppi =priority;
		}
		
	      }//proton 2 loop end
	    }//proton condition

	  }//scattering angle condition
      }//reaction loop end
    }//proton loop end
  }

  //end of SigmaP scattering assumption


  HF1( 1, 10. );

  //Final Hodoscope histograms
  for( int i=0; i<nhBh2; ++i ){
    HF1( 152, event.csBh2[i] );
    HF1( 153, event.Bh2Seg[i] );
    HF1( 154, event.tBh2[i] );
    HF1( 155, event.deBh2[i] );
  }

  for( int i=0; i<nhBh1; ++i ){
    HF1( 252, event.csBh1[i] );
    HF1( 253, event.Bh1Seg[i] );
    HF1( 254, event.tBh1[i] );
    HF1( 255, event.deBh1[i] );
  }

  for( int i=0; i<nhTof; ++i ){
    HF1( 352, event.csTof[i] );
    HF1( 353, event.TofSeg[i] );
    HF1( 354, event.tTof[i] );
    HF1( 355, event.deTof[i] );
  }

  return true;
}

//_____________________________________________________________________
bool calcDecayPiMom(double p1, double m1, double m2, double m3, double cost, double *momCal)
/*
  p1 : Sigma mom
  m1 : Sigma mass
  m2 : Neutron mass
  m3 : Pion mass
*/
{
  double E1 = sqrt(p1*p1+m1*m1);
  double A  = (m1*m1+m3*m3-m2*m2)/2.;
  double hanbetu = (A*p1*cost)*(A*p1*cost)-(E1*E1-p1*p1*cost*cost)*(E1*E1*m3*m3-A*A);
  if (hanbetu >= 0) {
    double ans1 = (A*p1*cost+sqrt(hanbetu))/(E1*E1-p1*p1*cost*cost);
    double ans2 = (A*p1*cost-sqrt(hanbetu))/(E1*E1-p1*p1*cost*cost);
    if (ans1>=0 && ans2<0) {
      *momCal = ans1;	  
      return true;
    } else if (ans2>=0 && ans1<0) {
      *momCal = ans2;	  
      return true;
          } else if (ans1>=0 && ans2>=0) {
      std::cout << "decayPiMomCal two answers: " << ans1 << ", " << ans2
      	<< std::endl;
      return false;
    } else if (ans1<0 && ans2<0) {
      std::cout << "decayPiMomCal two negative answers: " 
		<< ans1 << ", " << ans2
		<< std::endl;
      return false;
    }
  }

  return false;
}

//_____________________________________________________________________
//_____________________________________________________________________
bool calcDecayProtonMom(double p1, double m1, double m2, double m3, double cost, double *momCal1, double *momCal2)
/*
  p1 : Sigma mom
  m1 : Sigma mass
  m2 : PiZeroMass
  m3 : ProtonMass
*/
{
  double E1 = sqrt(p1*p1+m1*m1);
  double A  = (m1*m1+m3*m3-m2*m2)/2.;
  double hanbetu = (A*p1*cost)*(A*p1*cost)-(E1*E1-p1*p1*cost*cost)*(E1*E1*m3*m3-A*A);
  if (hanbetu >= 0) {
    double ans1 = (A*p1*cost+sqrt(hanbetu))/(E1*E1-p1*p1*cost*cost);
    double ans2 = (A*p1*cost-sqrt(hanbetu))/(E1*E1-p1*p1*cost*cost);
    if (ans1>=0 && ans2<0) {
      *momCal1 = ans1;	  
      return true;
    } else if (ans2>=0 && ans1<0) {
      *momCal2 = ans2;	  
      return true;
    } else if (ans1>=0 && ans2>=0) {
      *momCal1 = ans1;	  
      *momCal2 = ans2;	  
      return true;
    } else if (ans1<0 && ans2<0) {
      std::cout << "decayPiMomCal two negative answers: " 
		<< ans1 << ", " << ans2
		<< std::endl;
      return false;
    }
  }

  return false;
}

//_____________________________________________________________________
bool calc2BodyKinema(double M1, double p1, double M2, double phi,
		     double *scatMomCal, double *scatEkinCal, double *scatThetaCM)
/*
  M1 : sigma mass
  p1 : sigma mom
  M2 : proton mass
  phi : scat proton angle

*/
{
  double M = M2;
  double E1 = sqrt(M1*M1+p1*p1);
  double p4 = (2.*M*(E1+M)*p1*cos(phi*math::Deg2Rad()))/
    ((E1+M)*(E1+M)-p1*p1*cos(phi*math::Deg2Rad())*cos(phi*math::Deg2Rad()));
  double E4 = sqrt(M*M+p4*p4);
  double Ekin4 = E4 - M;
	    
  double beta = p1/(E1+M);
  double gamma = 1./sqrt(1.-beta*beta);

  double M3 = M1;
  double Ecm = sqrt((E1+M)*(E1+M)-p1*p1);
  double Pcm = 1./(2.*Ecm)*
    sqrt((Ecm-M3+M)*(Ecm-M3-M)*(Ecm+M3-M)*(Ecm+M3+M));
  double sintCM = p4*sin(phi*math::Deg2Rad())/Pcm;
  double thetaCM = asin(sintCM)*math::Rad2Deg();

  double P4cm_para = -beta*gamma*E4 + 
    gamma*p4*cos(phi*math::Deg2Rad());
  if (P4cm_para>0)
    thetaCM = 180. - thetaCM;

  *scatMomCal = p4;
  *scatEkinCal = Ekin4;
  *scatThetaCM = thetaCM;

  return true;
}

//_____________________________________________________________________
bool calc2BodyInelastic(double m1, double p1, double m2, double m3, double m4, 
			double theta, double *pCal1, double *pCal2, double *ThetaCM)
/*
  m1 = SigmaPlusMass(injected particle);
  m2 = Zero (Target);
  m3 = ProtonMass (Producted A, we want to know momentum);
  m4 = pi0Mass (Producted B);
  theta = thetaLambdaNConv;
*/
{

  double E1 = sqrt(p1*p1+m1*m1);
  double A = m1*m1+m2*m2+m3*m3-m4*m4+2.*E1*m2;

  double cost = cos(theta*math::Deg2Rad());

  double hanbetu = 4.*m3*m3*(p1*p1*cost*cost-(E1+m2)*(E1+m2))+A*A;
  if (hanbetu>=0) {
    double ans1 = (A*p1*cost+(E1+m2)*sqrt(hanbetu))/(2.*((E1+m2)*(E1+m2)-p1*p1*cost*cost));
    double ans2 = (A*p1*cost-(E1+m2)*sqrt(hanbetu))/(2.*((E1+m2)*(E1+m2)-p1*p1*cost*cost));
    *pCal1 = ans1;
    *pCal2 = ans2;

    /*
      std::cout << "---Calculated Lambda mom" << std::endl;
      if (ans1>=0 && ans2>= 0) {
      std::cout << "1 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
      } else if (ans1>=0) {
      std::cout << "2 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
      } else if (ans2>=0) {
      std::cout << "3 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
      } else {
      std::cout << "4 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
      }
    */		    
    double beta = p1/(E1+m2);
    double gamma = 1./sqrt(1.-beta*beta);
		    
    double Ecm = sqrt((E1+m2)*(E1+m2)-p1*p1);
    double Pcm = 1./(2.*Ecm)*
      sqrt((Ecm-m3+m4)*(Ecm-m3-m4)*(Ecm+m3-m4)*(Ecm+m3+m4));
    //double p3 = event.momCalLambda;
    double p3 = ans1;
    double E3 = sqrt(p3*p3+m3*m3);
    double sintCM = p3*sin(theta*math::Deg2Rad())/Pcm;
    double thetaCM = asin(sintCM)*math::Rad2Deg();
		    
    double P3cm_para = -beta*gamma*E3 + 
      gamma*p3*cos(theta*math::Deg2Rad());
    if (P3cm_para<0)
      thetaCM = 180. - thetaCM;
    
    *ThetaCM = thetaCM;
    
    return true;
  }

  return false;
}

//_____________________________________________________________________

ThreeVector Xproduct(ThreeVector vec1, ThreeVector vec2)
{
  double vec_x = vec1.y()*vec2.z()-vec2.y()*vec1.z();
  double vec_y = vec1.z()*vec2.x()-vec2.z()*vec1.x();
  double vec_z = vec1.x()*vec2.y()-vec2.x()*vec1.y();

  return ThreeVector(vec_x, vec_y, vec_z);
}

//_____________________________________________________________________

bool calcBeamMomFrom2BodyKinema(double M1, double M2, double theta, double pScat,
				double *beamMomCal1, double *beamMomCal2)
/*
  M1 : pi mass
  M2 : proton mass
  theta : scat proton angle
  pScat : scat proton momentum
*/
{
  double p3 = pScat;
  double E3 = sqrt(M2*M2+p3*p3);

  double A = (M2-E3)*(M2-E3)-p3*p3*cos(theta*math::Deg2Rad())*cos(theta*math::Deg2Rad());
  double B = (E3*M2-M2*M2)*p3*cos(theta*math::Deg2Rad());
  double C = (M2-E3)*(M2-E3)*M1*M1-(E3*M2-M2*M2)*(E3*M2-M2*M2);

  double hanbetsu = B*B-A*C;
  if (hanbetsu>=0) {
    double p1 = (-B+sqrt(hanbetsu))/A;
    double p2 = (-B-sqrt(hanbetsu))/A;

    //std::cout << "p1 : " << p1 << ", p2 : " << p2 << std::endl;

    *beamMomCal1 = p1;
    *beamMomCal2 = p2;
    return true;
  }

  return false;
}

//_____________________________________________________________________

bool calcBeamMomFromDecay(double M1, double M2, double M3, double p2, double theta, double *p1_1, double *p1_2)
/*
  M1: SigmaPlusMass
  M2: ProtonMass
  M3: PiZeroMass
  p2: Momentum of proton
  theta : theta of proton
*/
{
  double A = M1*M1+M2*M2-M3*M3;
  double E2 = sqrt(p2*p2+M2*M2);
  double B = E2*E2-p2*p2*cos(theta*math::Deg2Rad())*cos(theta*math::Deg2Rad());
  double C = A*p2*cos(theta*math::Deg2Rad());
  double D = M1*M1*E2*E2-A*A/4.;

  double hanbetsu = C*C-4*B*D;
  if (hanbetsu>=0) {
    
    *p1_1 = (C+sqrt(hanbetsu))/(2.*B);
    *p1_2 = (C-sqrt(hanbetsu))/(2.*B);

    return true;
  }
  
  return false;
}

//_____________________________________________________________________

double catchpidcurve(double param1, double param2, double x){
  return param1*pow(x+param2,2)/(pow(x,2)+2*param2*x);
}

//_____________________________________________________________________

bool
dst::DstClose( void )
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  const std::size_t n = TFileCont.size();
  for( std::size_t i=0; i<n; ++i ){
    if( TTreeCont[i] ) delete TTreeCont[i];
    if( TFileCont[i] ) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________
bool
ConfMan::InitializeHistograms( void )
{
  HB1( 1, "Status", 60, 0., 60. );

  HB1( 101, "#Clusters BH2",   7, 0., 7. );
  HB1( 102, "ClusterSize BH2", 7, 0., 7. );
  HB1( 103, "HitPat BH2", 10, 0., 10. );
  HB1( 104, "MeanTime BH2", 200, -10., 10. );
  HB1( 105, "Delta-E BH2", 200, -0.5, 4.5 );

  HB1( 112, "ClusterSize BH2", 7, 0., 7. );
  HB1( 113, "HitPat BH2", 10, 0., 10. );
  HB1( 114, "MeanTime BH2", 200, -10., 10. );
  HB1( 115, "Delta-E BH2", 200, -0.5, 4.5 );

  HB1( 122, "ClusterSize BH2 [T0]", 7, 0., 7. );
  HB1( 123, "HitPat BH2 [T0]", 10, 0., 10. );
  HB1( 124, "MeanTime BH2 [T0]", 200, -10., 10. );
  HB1( 125, "Delta-E BH2 [T0]", 200, -0.5, 4.5 );

  HB1( 152, "ClusterSize BH2 [kk]", 7, 0., 7. );
  HB1( 153, "HitPat BH2 [kk]", 10, 0., 10. );
  HB1( 154, "MeanTime BH2 [kk]", 200, -10., 10. );
  HB1( 155, "Delta-E BH2 [kk]", 200, -0.5, 4.5 );

  HB1( 201, "#Clusters BH1",  11, 0., 11. );
  HB1( 202, "ClusterSize BH1",11, 0., 11. );
  HB1( 203, "HitPat BH1", 11, 0., 11. );
  HB1( 204, "MeanTime BH1", 200, -10., 10. );
  HB1( 205, "Delta-E BH1", 200, -0.5, 4.5 );
  HB1( 206, "Beam ToF", 200, -10., 10. );

  HB1( 211, "#Clusters BH1 [pi]",  11, 0., 11. );
  HB1( 212, "ClusterSize BH1 [pi]",11, 0., 11. );
  HB1( 213, "HitPat BH1 [pi]", 11, 0., 11. );
  HB1( 214, "MeanTime BH1 [pi]", 200, -10., 10. );
  HB1( 215, "Delta-E BH1 [pi]", 200, -0.5, 4.5 );

  HB1( 252, "ClusterSize BH1 [kk]",11, 0., 11. );
  HB1( 253, "HitPat BH1 [kk]", 11, 0., 11. );
  HB1( 254, "MeanTime BH1 [kk]", 200, -10., 10. );
  HB1( 255, "Delta-E BH1 [kk]", 200, -0.5, 4.5 );
  HB1( 256, "Beam ToF [kk]", 200, -10., 10. );

  HB1( 301, "#Clusters Tof",  32, 0., 32. );
  HB1( 302, "ClusterSize Tof",32, 0., 32. );
  HB1( 303, "HitPat Tof", 32, 0., 32. );
  HB1( 304, "TimeOfFlight Tof", 500, -50., 100. );
  HB1( 305, "Delta-E Tof", 200, -0.5, 4.5 );

  HB1( 311, "#Clusters Tof [Good]",  32, 0., 32. );
  HB1( 312, "ClusterSize Tof [Good]",32, 0., 32. );
  HB1( 313, "HitPat Tof [Good]", 32, 0., 32. );
  HB1( 314, "TimeOfFlight Tof [Good]", 500, -50., 100. );
  HB1( 315, "Delta-E Tof [Good]", 200, -0.5, 4.5 );

  HB1( 352, "ClusterSize Tof [kk]",32, 0., 32. );
  HB1( 353, "HitPat Tof [kk]", 32, 0., 32. );
  HB1( 354, "TimeOfFlight Tof [kk]", 500, -50., 100. );
  HB1( 355, "Delta-E Tof [kk]", 200, -0.5, 4.5 );
  HB2( 501, "SegLC%SegTOF", 32, 0., 32., 28, 0., 28. );

  HB1( 1001, "#Tracks SdcIn", 10, 0., 10. );
  HB1( 1002, "#Hits SdcIn", 20, 0., 20. );
  HB1( 1003, "Chisqr SdcIn", 200, 0., 100. );
  HB1( 1004, "X0 SdcIn", 500, -200., 200. );
  HB1( 1005, "Y0 SdcIn", 500, -200., 200. );
  HB1( 1006, "U0 SdcIn",  700, -0.35, 0.35 );
  HB1( 1007, "V0 SdcIn",  400, -0.20, 0.20 );
  HB2( 1008, "X0%U0 SdcIn", 120, -200., 200., 100, -0.35, 0.35 );
  HB2( 1009, "Y0%V0 SdcIn", 100, -200., 200., 100, -0.20, 0.20 );
  HB2( 1010, "X0%Y0 SdcIn", 100, -200., 200., 100, -200., 200. );

  HB1( 1101, "#Tracks BcOut", 10, 0., 10. );
  HB1( 1102, "#Hits BcOut", 20, 0., 20. );
  HB1( 1103, "Chisqr BcOut", 200, 0., 100. );
  HB1( 1104, "X0 BcOut", 500, -200., 200. );
  HB1( 1105, "Y0 BcOut", 500, -200., 200. );
  HB1( 1106, "U0 BcOut",  700, -0.35, 0.35 );
  HB1( 1107, "V0 BcOut",  400, -0.20, 0.20 );
  HB2( 1108, "X0%U0 BcOut", 120, -200., 200., 100, -0.35, 0.35 );
  HB2( 1109, "Y0%V0 BcOut", 100, -200., 200., 100, -0.20, 0.20 );
  HB2( 1110, "X0%Y0 BcOut", 100, -200., 200., 100, -200., 200. );
  HB1( 1111, "Xtgt BcOut", 500, -200., 200. );
  HB1( 1112, "Ytgt BcOut", 500, -200., 200. );
  HB2( 1113, "Xtgt%Ytgt BcOut", 100, -200., 200., 100, -200., 200. );

  HB1( 1201, "#Tracks SdcOut", 10, 0., 10. );
  HB1( 1202, "#Hits SdcOut", 20, 0., 20. );
  HB1( 1203, "Chisqr SdcOut", 200, 0., 100. );
  HB1( 1204, "X0 SdcOut", 600, -1200., 1200. );
  HB1( 1205, "Y0 SdcOut", 600, -600., 600. );
  HB1( 1206, "U0 SdcOut",  700, -0.35, 0.35 );
  HB1( 1207, "V0 SdcOut",  400, -0.20, 0.20 );
  HB2( 1208, "X0%U0 SdcOut", 120, -1200., 1200., 100, -0.35, 0.35 );
  HB2( 1209, "Y0%V0 SdcOut", 100,  -600.,  600., 100, -0.20, 0.20 );
  HB2( 1210, "X0%Y0 SdcOut", 100, -1200., 1200., 100, -600., 600. );

  HB1( 2201, "#Tracks K18", 10, 0., 10. );
  HB1( 2202, "#Hits K18", 30, 0., 30. );
  HB1( 2203, "Chisqr K18", 500, 0., 50. );
  HB1( 2204, "P K18", 1000, 0.5, 2.0 );
  HB1( 2251, "#Tracks K18 [Good]", 10, 0., 10. );

  HB1( 3001, "#Tracks Kurama", 10, 0., 10. );
  HB1( 3002, "#Hits Kurama", 30, 0., 30. );
  HB1( 3003, "Chisqr Kurama", 500, 0., 500. );
  HB1( 3004, "Xtgt Kurama", 500, -200., 200. );
  HB1( 3005, "Ytgt Kurama", 500, -100., 100. );
  HB1( 3006, "Utgt Kurama", 400, -0.35, 0.35 );
  HB1( 3007, "Vtgt Kurama", 200, -0.20, 0.20 );
  HB2( 3008, "Xtgt%U Kurama", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 3009, "Ytgt%V Kurama", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 3010, "Xtgt%Ytgt Kurama", 100, -200., 200., 100, -100., 100. );
  HB1( 3011, "P Kurama", 200, 0.0, 1.0 );
  HB1( 3012, "PathLength Kurama", 600, 3000., 6000. );
  HB1( 3013, "MassSqr", 600, -1.2, 1.2 );

  HB1( 3101, "#Tracks Kurama [Good]", 10, 0., 10. );
  HB1( 3102, "#Hits Kurama [Good]", 30, 0., 30. );
  HB1( 3103, "Chisqr Kurama [Good]", 500, 0., 500. );
  HB1( 3104, "Xtgt Kurama [Good]", 500, -200., 200. );
  HB1( 3105, "Ytgt Kurama [Good]", 500, -100., 100. );
  HB1( 3106, "Utgt Kurama [Good]", 700, -0.35, 0.35 );
  HB1( 3107, "Vtgt Kurama [Good]", 400, -0.20, 0.20 );
  HB2( 3108, "Xtgt%U Kurama [Good]", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 3109, "Ytgt%V Kurama [Good]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 3110, "Xtgt%Ytgt Kurama [Good]", 100, -200., 200., 100, -100., 100. );
  HB1( 3111, "P Kurama [Good]", 200, 0.0, 1.0 );
  HB1( 3112, "PathLength Kurama [Good]", 600, 3000., 6000. );
  HB1( 3113, "MassSqr", 600, -1.2, 1.2 );

  HB1( 3201, "#Tracks Kurama [Good2]", 10, 0., 10. );
  HB1( 3202, "#Hits Kurama [Good2]", 30, 0., 30. );
  HB1( 3203, "Chisqr Kurama [Good2]", 500, 0., 500. );
  HB1( 3204, "Xtgt Kurama [Good2]", 500, -200., 200. );
  HB1( 3205, "Ytgt Kurama [Good2]", 500, -100., 100. );
  HB1( 3206, "Utgt Kurama [Good2]", 700, -0.35, 0.35 );
  HB1( 3207, "Vtgt Kurama [Good2]", 400, -0.20, 0.20 );
  HB2( 3208, "Xtgt%U Kurama [Good2]", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 3209, "Ytgt%V Kurama [Good2]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 3210, "Xtgt%Ytgt Kurama [Good2]", 100, -200., 200., 100, -100., 100. );
  HB1( 3211, "P Kurama [Good2]", 200, 0.0, 1.0 );
  HB1( 3212, "PathLength Kurama [Good2]", 600, 3000., 6000. );
  HB1( 3213, "MassSqr", 600, -1.2, 1.2 );

  HB1( 4101, "#Tracks K18 [KuramaP]", 10, 0., 10. );
  HB1( 4102, "#Hits K18 [KuramaP]", 30, 0., 30. );
  HB1( 4103, "Chisqr K18 [KuramaP]", 500, 0., 50. );
  HB1( 4104, "P K18 [KuramaP]", 500, 0.5, 2.0 );
  HB1( 4105, "Xtgt K18 [KuramaP]", 500, -200., 200. );
  HB1( 4106, "Ytgt K18 [KuramaP]", 500, -100., 100. );
  HB1( 4107, "Xout K18 [KuramaP]", 500, -200., 200. );
  HB1( 4108, "Yout K18 [KuramaP]", 500, -100., 100. );
  HB1( 4109, "Uout K18 [KuramaP]", 400, -0.35, 0.35 );
  HB1( 4110, "Vout K18 [KuramaP]", 200, -0.20, 0.20 );
  HB1( 4111, "Xin  K18 [KuramaP]", 500, -100., 100. );
  HB1( 4112, "Yin  K18 [KuramaP]", 500, -100., 100. );
  HB1( 4113, "Uin  K18 [KuramaP]", 400, -0.35, 0.35 );
  HB1( 4114, "Vin  K18 [KuramaP]", 200, -0.20, 0.20 );

  HB1( 4201, "#Tracks Kurama [Proton]", 10, 0., 10. );
  HB1( 4202, "#Hits Kurama [Proton]", 30, 0., 30. );
  HB1( 4203, "Chisqr Kurama [Proton]", 500, 0., 500. );
  HB1( 4204, "Xtgt Kurama [Proton]", 500, -200., 200. );
  HB1( 4205, "Ytgt Kurama [Proton]", 500, -100., 100. );
  HB1( 4206, "Utgt Kurama [Proton]", 700, -0.35, 0.35 );
  HB1( 4207, "Vtgt Kurama [Proton]", 400, -0.20, 0.20 );
  HB2( 4208, "Xtgt%U Kurama [Proton]", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 4209, "Ytgt%V Kurama [Proton]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 4210, "Xtgt%Ytgt Kurama [Proton]", 100, -200., 200., 100, -100., 100. );
  HB1( 4211, "P Kurama [Proton]", 200, 0.0, 1.0 );
  HB1( 4212, "PathLength Kurama [Proton]", 600, 3000., 6000. );
  HB1( 4213, "MassSqr", 600, -1.2, 1.2 );

  HB1( 5001, "Zvert [KK]", 1000, -1000., 1000. );
  HB1( 5002, "MissingMass [PiK]", 400, 1.0, 2.0 );
  HB1( 5003, "MissingMass [PiK] 1.15 to 1.25", 400, 1.0, 2.0 );

  HB2( 5011, "MissingMass%Us", 200, 0.0, 2.50, 100, -0.40, 0.40 );
  HB2( 5012, "MissingMass%Vs", 200, 0.0, 2.50, 100, -0.20, 0.20 );
  HB2( 5013, "MissingMass%Ub", 200, 0.0, 2.50, 100, -0.30, 0.30 );
  HB2( 5014, "MissingMass%Vb", 200, 0.0, 2.50, 100, -0.10, 0.10 );

  ////////////////////////////////////////////
  //Tree
  HBTree("pik","tree of PiKAna");
  tree->Branch("runnum", &event.runnum, "runnum/I");
  tree->Branch("evnum",  &event.evnum,  "evnum/I");
  tree->Branch("spill",  &event.spill,  "spill/I");
  //Trigger
  tree->Branch("trignhits", &event.trignhits, "trignhits/I");
  tree->Branch("trigpat",   event.trigpat, "trigpat[trignhits]/I");
  tree->Branch("trigflag",  event.trigflag, Form( "trigflag[%d]/I", NumOfSegTrig ) );

  //Hodoscope
  tree->Branch("nhBh1",   &event.nhBh1,   "nhBh1/I");
  tree->Branch("csBh1",    event.csBh1,   "csBh1[nhBh1]/I");
  tree->Branch("Bh1Seg",   event.Bh1Seg,  "Bh1Seg[nhBh1]/D");
  tree->Branch("tBh1",     event.tBh1,    "tBh1[nhBh1]/D");
  tree->Branch("dtBh1",    event.dtBh1,   "dtBh1[nhBh1]/D");
  tree->Branch("deBh1",    event.deBh1,   "deBh1[nhBh1]/D");
  tree->Branch("btof",     event.btof,    "btof[nhBh1]/D");

  tree->Branch("nhBh2",   &event.nhBh2,   "nhBh2/I");
  tree->Branch("csBh2",    event.csBh2,   "csBh2[nhBh2]/I");
  tree->Branch("Bh2Seg",   event.Bh2Seg,  "Bh2Seg[nhBh2]/D");
  tree->Branch("tBh2",     event.tBh2,    "tBh2[nhBh2]/D");
  tree->Branch("t0Bh2",    event.t0Bh2,   "t0Bh2[nhBh2]/D");
  tree->Branch("dtBh2",    event.dtBh2,   "dtBh2[nhBh2]/D");
  tree->Branch("deBh2",    event.deBh2,   "deBh2[nhBh2]/D");

  tree->Branch("nhTof",   &event.nhTof,   "nhTof/I");
  tree->Branch("csTof",    event.csTof,   "csTof[nhTof]/I");
  tree->Branch("TofSeg",   event.TofSeg,  "TofSeg[nhTof]/D");
  tree->Branch("tTof",     event.tTof,    "tTof[nhTof]/D");
  tree->Branch("dtTof",    event.dtTof,   "dtTof[nhTof]/D");
  tree->Branch("deTof",    event.deTof,   "deTof[nhTof]/D");

  /*
  //Fiber
  tree->Branch("nhBft",  &event.nhBft,  "nhBft/I");
  tree->Branch("csBft",   event.csBft,  "csBft[nhBft]/I");
  tree->Branch("tBft",    event.tBft,   "tBft[nhBft]/D");
  tree->Branch("wBft",    event.wBft,   "wBft[nhBft]/D");
  tree->Branch("BftPos",  event.BftPos, "BftPos[nhBft]/D");
  tree->Branch("nhSch",  &event.nhSch,  "nhSch/I");
  tree->Branch("csSch",   event.csSch,  "csSch[nhSch]/I");
  tree->Branch("tSch",    event.tSch,   "tSch[nhSch]/D");
  tree->Branch("wSch",    event.wSch,   "wSch[nhSch]/D");
  tree->Branch("NSch",   &event.NSch,   "NSch/I");
  tree->Branch("SchSeg",  event.SchSeg, "SchSeg[NSch]/I");
  */

  //Beam DC
  tree->Branch("bft_ncl",   &event.bft_ncl,     "bft_ncl/I");
  tree->Branch("nlBcOut",   &event.nlBcOut,     "nlBcOut/I");
  tree->Branch("ntBcOut",   &event.ntBcOut,     "ntBcOut/I");
  tree->Branch("nhBcOut",    event.nhBcOut,     "nhBcOut[ntBcOut]/I");
  tree->Branch("chisqrBcOut",event.chisqrBcOut, "chisqrBcOut[ntBcOut]/D");
  tree->Branch("x0BcOut",    event.x0BcOut,     "x0BcOut[ntBcOut]/D");
  tree->Branch("y0BcOut",    event.y0BcOut,     "y0BcOut[ntBcOut]/D");
  tree->Branch("u0BcOut",    event.u0BcOut,     "u0BcOut[ntBcOut]/D");
  tree->Branch("v0BcOut",    event.v0BcOut,     "v0BcOut[ntBcOut]/D");
  tree->Branch("xtgtBcOut",  event.xtgtBcOut,   "xtgtBcOut[ntBcOut]/D");
  tree->Branch("ytgtBcOut",  event.ytgtBcOut,   "ytgtBcOut[ntBcOut]/D");
  tree->Branch("xbh2BcOut",  event.xbh2BcOut,   "xbh2BcOut[ntBcOut]/D");
  tree->Branch("ybh2BcOut",  event.ybh2BcOut,   "ybh2BcOut[ntBcOut]/D");

  tree->Branch("ntK18",      &event.ntK18,     "ntK18/I");
  tree->Branch("nhK18",       event.nhK18,     "nhK18[ntK18]/I");
  tree->Branch("chisqrK18",   event.chisqrK18, "chisqrK18[ntK18]/D");
  tree->Branch("pK18",        event.pK18,      "pK18[ntK18]/D");
  tree->Branch("xtgtK18",     event.xtgtK18,   "xtgtK18[ntK18]/D");
  tree->Branch("ytgtK18",     event.ytgtK18,   "ytgtK18[ntK18]/D");
  tree->Branch("utgtK18",     event.utgtK18,   "utgtK18[ntK18]/D");
  tree->Branch("vtgtK18",     event.vtgtK18,   "vtgtK18[ntK18]/D");
  tree->Branch("thetaK18",    event.thetaK18,  "thetaK18[ntK18]/D");

  //KURAMA
  tree->Branch("nlSdcIn",    &event.nlSdcIn,     "nlSdcIn/I");
  tree->Branch("ntSdcIn",    &event.ntSdcIn,     "ntSdcIn/I");
  tree->Branch("nhSdcIn",     event.nhSdcIn,     "nhSdcIn[ntSdcIn]/I");
  tree->Branch("chisqrSdcIn", event.chisqrSdcIn, "chisqrSdcIn[ntSdcIn]/D");
  tree->Branch("x0SdcIn",     event.x0SdcIn,     "x0SdcIn[ntSdcIn]/D");
  tree->Branch("y0SdcIn",     event.y0SdcIn,     "y0SdcIn[ntSdcIn]/D");
  tree->Branch("u0SdcIn",     event.u0SdcIn,     "u0SdcIn[ntSdcIn]/D");
  tree->Branch("v0SdcIn",     event.v0SdcIn,     "v0SdcIn[ntSdcIn]/D");

  tree->Branch("nlSdcOut",   &event.nlSdcOut,     "nlSdcOut/I");
  tree->Branch("ntSdcOut",   &event.ntSdcOut,     "ntSdcOut/I");
  tree->Branch("nhSdcOut",    event.nhSdcOut,     "nhSdcOut[ntSdcOut]/I");
  tree->Branch("chisqrSdcOut",event.chisqrSdcOut, "chisqrOut[ntSdcOut]/D");
  tree->Branch("x0SdcOut",    event.x0SdcOut,     "x0SdcOut[ntSdcOut]/D");
  tree->Branch("y0SdcOut",    event.y0SdcOut,     "y0SdcOut[ntSdcOut]/D");
  tree->Branch("u0SdcOut",    event.u0SdcOut,     "u0SdcOut[ntSdcOut]/D");
  tree->Branch("v0SdcOut",    event.v0SdcOut,     "v0SdcOut[ntSdcOut]/D");

  tree->Branch("ntKurama",     &event.ntKurama,     "ntKurama/I");
  tree->Branch("nhKurama",      event.nhKurama,     "nhKurama[ntKurama]/I");
  tree->Branch("stof",          event.stof,         "stof[ntKurama]/D");
  tree->Branch("cstof",         event.cstof,         "cstof[ntKurama]/D");
  tree->Branch("path",          event.path,         "path[ntKurama]/D");
  tree->Branch("chisqrKurama",  event.chisqrKurama, "chisqrKurama[ntKurama]/D");
  tree->Branch("pKurama",       event.pKurama,      "pKurama[ntKurama]/D");
  tree->Branch("qKurama",       event.qKurama,      "qKurama[ntKurama]/D");
  tree->Branch("m2",            event.m2,           "m2[ntKurama]/D");
  tree->Branch("cm2",           event.cm2,          "cm2[ntKurama]/D");
  tree->Branch("xtgtKurama",    event.xtgtKurama,   "xtgtKurama[ntKurama]/D");
  tree->Branch("ytgtKurama",    event.ytgtKurama,   "ytgtKurama[ntKurama]/D");
  tree->Branch("utgtKurama",    event.utgtKurama,   "utgtKurama[ntKurama]/D");
  tree->Branch("vtgtKurama",    event.vtgtKurama,   "vtgtKurama[ntKurama]/D");
  tree->Branch("thetaKurama",   event.thetaKurama,  "thetaKurama[ntKurama]/D");
  tree->Branch("xtofKurama",    event.xtofKurama,   "xtofKurama[ntKurama]/D");
  tree->Branch("ytofKurama",    event.ytofKurama,   "ytofKurama[ntKurama]/D");
  tree->Branch("utofKurama",    event.utofKurama,   "utofKurama[ntKurama]/D");
  tree->Branch("vtofKurama",    event.vtofKurama,   "vtofKurama[ntKurama]/D");
  tree->Branch("tofsegKurama",  event.tofsegKurama, "tofsegKurama[ntKurama]/D");

  //CFT
  tree->Branch("ntCFT",         &event.ntCFT,          "ntCFT/I");
  tree->Branch("ntProton",      &event.ntProton,       "ntProton/I");
  tree->Branch("ntOther",       &event.ntOther,        "ntOther/I");
  tree->Branch("theta_cft",      event.theta_cft,      "theta_cft[ntCFT]/D");
  //  tree->Branch("phi_cft",        event.phi_cft,        "phi_cft[ntCFT]/D");
  //  tree->Branch("vtx_cft",        event.vtx_cft,        "vtx_cft[ntCFT]/D");
  //  tree->Branch("vty_cft",        event.vty_cft,        "vty_cft[ntCFT]/D");
  //  tree->Branch("vtz_cft",        event.vtz_cft,        "vtz_cft[ntCFT]/D");
  tree->Branch("Total_dE",       event.Total_dE,       "Total_dE[ntCFT]/D");
  tree->Branch("Total_dEphi",    event.Total_dEphi,    "Total_dEphi[ntCFT]/D");
  tree->Branch("Total_dEuv",     event.Total_dEuv,     "Total_dEuv[ntCFT]/D");
  tree->Branch("Total_dE_max",   event.Total_dE_max,   "Total_dE_max[ntCFT]/D");
  tree->Branch("Total_dEphi_max",event.Total_dEphi_max,"Total_dEphi_max[ntCFT]/D");
  tree->Branch("Total_dEuv_max", event.Total_dEuv_max, "Total_dEuv_max[ntCFT]/D");
  tree->Branch("nhit_phi",       event.nhit_phi,       "nhit_phi[ntCFT]/I");
  tree->Branch("nhit_uv" ,       event.nhit_uv,        "nhit_uv[ntCFT]/I");
  tree->Branch("segBGOt" ,       event.segBGOt,        "segBGOt[ntCFT]/I");
  tree->Branch("energyBGOt",     event.energyBGOt,     "energyBGOt[ntCFT]/D");
  tree->Branch("segPiIDt",       event.segPiIDt,       "segPiIDt[ntCFT]/I");
  tree->Branch("protonflagt",     event.protonflagt,     "protonflagt[ntCFT]/I");
  tree->Branch("BGOnohitt",      event.BGOnohitt,      "BGOnohitt[ntCFT]/I");
  tree->Branch("simKuramat",      event.simKuramat,      "simKuramat[ntCFT]/I");
  tree->Branch("tracknoproton",       event.tracknoproton,       "tracknoproton[ntProton]/I");
  tree->Branch("tracknoother",       event.tracknoother,       "tracknoother[ntOther]/I");

  //Reaction
  tree->Branch("nPi",           &event.nPi,            "nPi/I");
  tree->Branch("nK",            &event.nK,             "nK/I");
  tree->Branch("nPiK",          &event.nPiK,           "nPiK/I");
  tree->Branch("nCatch",        &event.nCatch,         "nCatch/I");
  tree->Branch("nPiKCatch",     &event.nPiKCatch,      "nPiKCatch/I");

  tree->Branch("vtx_KURAMA",            event.vtx_KURAMA,            "vtx_KURAMA[nPiK]/D");
  tree->Branch("vty_KURAMA",            event.vty_KURAMA,            "vty_KURAMA[nPiK]/D");
  tree->Branch("vtz_KURAMA",            event.vtz_KURAMA,            "vtz_KURAMA[nPiK]/D");
  tree->Branch("ccm2",                  event.ccm2,                  "ccm2[nPiK]/D");
  tree->Branch("chisqrKuramapik",       event.chisqrKuramapik,       "chisqrKuramapik[nPiK]/D");
  tree->Branch("pKuramapik"     ,       event.pKuramapik,            "pKuramapik[nPiK]/D");
  tree->Branch("closedist_KURAMA",      event.closedist_KURAMA,      "closedist_KURAMA[nPiKCatch]/D");
  tree->Branch("theta",          event.theta,          "theta[nPiK]/D");
  tree->Branch("MissMass",       event.MissMass,       "MissMass[nPiK]/D");
  tree->Branch("MissMassCorr",   event.MissMassCorr,   "MissMassCorr[nPiK]/D");
  tree->Branch("MissMassCorrDE", event.MissMassCorrDE, "MissMassCorrDE[nPiK]/D");
  tree->Branch("thetaCM", event.thetaCM,  "thetaCM[nPiK]/D");
  tree->Branch("costCM",  event.costCM,   "costCM[nPiK]/D");

  tree->Branch("MissMassPiPi",       event.MissMassPiPi,       "MissMassPiPi[nPiK]/D");
  tree->Branch("MissMassPiP",       event.MissMassPiP,       "MissMassPiP[nPiK]/D");

  tree->Branch("KURAMAPID", event.KURAMAPID, "KURAMAPID[nPiKCatch]/I");

  tree->Branch("MissMom", event.MissMom, "MissMom[nPiKCatch]/D");
  tree->Branch("MissMomx", event.MissMomx, "MissMomx[nPiKCatch]/D");
  tree->Branch("MissMomy", event.MissMomy, "MissMomy[nPiKCatch]/D");
  tree->Branch("MissMomz", event.MissMomz, "MissMomz[nPiKCatch]/D");
  tree->Branch("MissMomcal", event.MissMomcal, "MissMomcal[nPiKCatch]/D");
  tree->Branch("SigmaBeta", event.SigmaBeta, "SigmaBeta[nPiKCatch]/D");

  tree->Branch("xpi",        event.xpi,      "xpi[nPiKCatch]/D");
  tree->Branch("ypi",        event.ypi,      "ypi[nPiKCatch]/D");
  tree->Branch("upi",        event.upi,      "upi[nPiKCatch]/D");
  tree->Branch("vpi",        event.vpi,      "vpi[nPiKCatch]/D");
  tree->Branch("xk",         event.xk,       "xk[nPiKCatch]/D");
  tree->Branch("yk",         event.yk,       "yk[nPiKCatch]/D");
  tree->Branch("uk",         event.uk,       "uk[nPiKCatch]/D");
  tree->Branch("vk",         event.vk,       "vk[nPiKCatch]/D");
  tree->Branch("uc",         event.uc,       "uc[nPiKCatch]/D");
  tree->Branch("vc",         event.vc,       "vc[nPiKCatch]/D");
  tree->Branch("pOrg",       event.pOrg,      "pOrg[nPiKCatch]/D");
  tree->Branch("pCalc",      event.pCalc,     "pCalc[nPiKCatch]/D");
  tree->Branch("pCorr",      event.pCorr,     "pCorr[nPiKCatch]/D");
  tree->Branch("pCorrDE",    event.pCorrDE,   "pCorrDE[nPiKCatch]/D");

  //  tree->Branch("vtx_K18Catch", event.vtx_K18Catch ,"vtx_K18Catch[nPiKCatch]/D");
  //  tree->Branch("vty_K18Catch", event.vty_K18Catch ,"vty_K18Catch[nPiKCatch]/D");
  //  tree->Branch("vtz_K18Catch", event.vtz_K18Catch ,"vtz_K18Catch[nPiKCatch]/D");
  //  tree->Branch("closedist_K18Catch", event.closedist_K18Catch, "closedist_K18Catch[nPiKCatch]/D");
  tree->Branch("theta_K18Catch", event.theta_K18Catch, "theta_K18Catch[nPiKCatch]/D");
  tree->Branch("vertex_distance", event.vertex_distance, "vertex_distance[nPiKCatch]/D");
  //  tree->Branch("vtx_XCatch", event.vtx_XCatch,"vtx_XCatch[nPiKCatch]/D");
  //  tree->Branch("vty_XCatch", event.vty_XCatch,"vty_XCatch[nPiKCatch]/D");
  //  tree->Branch("vtz_XCatch", event.vtz_XCatch,"vtz_XCatch[nPiKCatch]/D");
  //  tree->Branch("closedist_XCatch", event.closedist_XCatch, "closedist_XCatch[nPiKCatch]/D");
  //  tree->Branch("theta_XCatch", event.theta_XCatch, "theta_XCatch[nPiKCatch]/D");

  tree->Branch("DeltaE_NPScat", &event.DeltaE_NPScat, "DeltaE_NPScat/D");
  tree->Branch("ProtonMom_NPScat", &event.ProtonMom_NPScat, "ProtonMom_NPScat/D");
  tree->Branch("DecayNeutronMom", &event.DecayNeutronMom, "DecayNeutronMom/D");
  tree->Branch("DecayPionMom", &event.DecayPionMom, "DecayPionMom/D");
  tree->Branch("MissMassSigmaP_NPScat", &event.MissMassSigmaP_NPScat, "MissMassSigmaP_NPScat/D");
  tree->Branch("vtx_Decay2np", &event.vtx_Decay2np, "vtx_Decay2np/D");
  tree->Branch("vty_Decay2np", &event.vty_Decay2np, "vty_Decay2np/D");
  tree->Branch("vtz_Decay2np", &event.vtz_Decay2np, "vtz_Decay2np/D");
  tree->Branch("vtx_NPScat", &event.vtx_NPScat, "vtx_NPScat/D");
  tree->Branch("vty_NPScat", &event.vty_NPScat, "vty_NPScat/D");
  tree->Branch("vtz_NPScat", &event.vtz_NPScat, "vtz_NPScat/D");
  tree->Branch("cdistDecay2np", &event.cdistDecay2np, "cdistDecay2np/D");
  tree->Branch("cdistNPScat",  &event.cdistNPScat , "cdistNPScat/D");
  tree->Branch("theta_NPScat", &event.theta_NPScat, "theta_NPScat/D");
  tree->Branch("thetaCM_NPScat", &event.thetaCM_NPScat, "thetaCM_NPScat/D");
  tree->Branch("ptrno_NPScat",  &event.ptrno_NPScat , "ptrno_NPScat/I");
  tree->Branch("pitrno_NPScat",  &event.pitrno_NPScat , "pitrno_NPScat/I");
  tree->Branch("pikno_NPScat",  &event.pikno_NPScat , "pikno_NPScat/I");
  tree->Branch("priority_NPScat", &event.priority_NPScat , "priority_NPScat/I");

  tree->Branch("DeltaE_DecayP", &event.DeltaE_DecayP, "DeltaE_DecayP/D");
  tree->Branch("ProtonMom_DecayP", &event.ProtonMom_DecayP, "ProtonMom_DecayP/D");
  tree->Branch("MissMassSigmaP_DecayP", &event.MissMassSigmaP_DecayP, "MissMassSigmaP_DecayP/D");
  tree->Branch("cdistDecayP",  &event.cdistDecayP , "cdistDecayP/D");
  tree->Branch("SigmaMomcor_DecayP", &event.SigmaMomcor_DecayP, "SigmaMomcor_DecayP/D");
  tree->Branch("SigmaLength",  &event.SigmaLength , "SigmaLength/D");
  tree->Branch("SigmaLengthcor",  &event.SigmaLengthcor , "SigmaLengthcor/D");
  tree->Branch("theta_DecayP", &event.theta_DecayP, "theta_DecayP/D");
  tree->Branch("ptrno_DecayP",  &event.ptrno_DecayP , "ptrno_DecayP/I");
  tree->Branch("pikno_DecayP",  &event.pikno_DecayP , "pikno_DecayP/I");
  tree->Branch("priority_DecayP", &event.priority_DecayP , "priority_DecayP/I");

  tree->Branch("DeltaP_PiPScat", &event.DeltaP_PiPScat, "DeltaP_PiPScat/D");
  tree->Branch("DecayNeutronMom2", &event.DecayNeutronMom2, "DecayNeutronMom2/D");
  tree->Branch("DecayPionMom2", &event.DecayPionMom2, "DecayPionMom2/D");
  tree->Branch("vtx_Decay2pip", &event.vtx_Decay2pip, "vtx_Decay2pip/D");
  tree->Branch("vty_Decay2pip", &event.vty_Decay2pip, "vty_Decay2pip/D");
  tree->Branch("vtz_Decay2pip", &event.vtz_Decay2pip, "vtz_Decay2pip/D");
  tree->Branch("vtx_PiPScat", &event.vtx_PiPScat, "vtx_PiPScat/D");
  tree->Branch("vty_PiPScat", &event.vty_PiPScat, "vty_PiPScat/D");
  tree->Branch("vtz_PiPScat", &event.vtz_PiPScat, "vtz_PiPScat/D");
  tree->Branch("cdistDecay2pip", &event.cdistDecay2pip, "cdistDecay2pip/D");
  tree->Branch("cdistPiPScat",  &event.cdistPiPScat , "cdistPiPScat/D");
  tree->Branch("ptrno_PiPScat",  &event.ptrno_PiPScat , "ptrno_PiPScat/I");
  tree->Branch("pitrno_PiPScat",  &event.pitrno_PiPScat , "pitrno_PiPScat/I");
  tree->Branch("pikno_PiPScat",  &event.pikno_PiPScat , "pikno_PiPScat/I");
  tree->Branch("priority_PiPScat", &event.priority_PiPScat , "priority_PiPScat/I");

  tree->Branch("DeltaP_PPScat", &event.DeltaP_PPScat, "DeltaP_PPScat/D");
  tree->Branch("theta_PPScat", &event.theta_PPScat, "theta_PPScat/D");
  tree->Branch("vtx_Decay2pp", &event.vtx_Decay2pp, "vtx_Decay2pp/D");
  tree->Branch("vty_Decay2pp", &event.vty_Decay2pp, "vty_Decay2pp/D");
  tree->Branch("vtz_Decay2pp", &event.vtz_Decay2pp, "vtz_Decay2pp/D");
  tree->Branch("vtx_PPScat", &event.vtx_PPScat, "vtx_PPScat/D");
  tree->Branch("vty_PPScat", &event.vty_PPScat, "vty_PPScat/D");
  tree->Branch("vtz_PPScat", &event.vtz_PPScat, "vtz_PPScat/D");
  tree->Branch("cdistDecay2pp", &event.cdistDecay2pp, "cdistDecay2pp/D");
  tree->Branch("cdistPPScat",  &event.cdistPPScat , "cdistPPScat/D");
  tree->Branch("ptrno_PPScat",  &event.ptrno_PPScat , "ptrno_PPScat/I");
  tree->Branch("p2trno_PPScat",  &event.p2trno_PPScat , "p2trno_PPScat/I");
  tree->Branch("pikno_PPScat",  &event.pikno_PPScat , "pikno_PPScat/I");
  tree->Branch("priority_PPScat", &event.priority_PPScat , "priority_PPScat/I");

  tree->Branch("DeltaE_SigmaPScat", &event.DeltaE_SigmaPScat, "DeltaE_SigmaPScat/D");
  tree->Branch("MissMassSigmaP_SigmaPScat", &event.MissMassSigmaP_SigmaPScat, "MissMassSigmaP_SigmaPScat/D");
  tree->Branch("vtx_SigmaPScat", &event.vtx_SigmaPScat, "vtx_SigmaPScat/D");
  tree->Branch("vty_SigmaPScat", &event.vty_SigmaPScat, "vty_SigmaPScat/D");
  tree->Branch("vtz_SigmaPScat", &event.vtz_SigmaPScat, "vtz_SigmaPScat/D");
  tree->Branch("cdistSigmaPScat", &event.cdistSigmaPScat, "cdistSigmaPScat/D");
  tree->Branch("ptrno_SigmaPScat",  &event.ptrno_SigmaPScat , "ptrno_SigmaPScat/I");
  tree->Branch("pikno_SigmaPScat",  &event.pikno_SigmaPScat , "pikno_SigmaPScat/I");

  tree->Branch("DeltaE_SigmaPScat2npi", &event.DeltaE_SigmaPScat2npi, "DeltaE_SigmaPScat2npi/D");
  tree->Branch("ProtonMom_SigmaPScat2npi", &event.ProtonMom_SigmaPScat2npi, "ProtonMom_SigmaPScat2npi/D");
  tree->Branch("MissMassSigmaP_SigmaPScat2npi", &event.MissMassSigmaP_SigmaPScat2npi, "MissMassSigmaP_SigmaPScat2npi/D");
  tree->Branch("vtx_ScatSigmaDecay2npi", &event.vtx_ScatSigmaDecay2npi, "vtx_ScatSigmaDecay2npi/D");
  tree->Branch("vty_ScatSigmaDecay2npi", &event.vty_ScatSigmaDecay2npi, "vty_ScatSigmaDecay2npi/D");
  tree->Branch("vtz_ScatSigmaDecay2npi", &event.vtz_ScatSigmaDecay2npi, "vtz_ScatSigmaDecay2npi/D");
  tree->Branch("cdistScatSigmaDecay2npi", &event.cdistScatSigmaDecay2npi, "cdistScatSigmaDecay2npi/D");
  tree->Branch("vdistance_SigmaPScat2npi", &event.vdistance_SigmaPScat2npi, "vdistance_SigmaPScat2npi/D");
  tree->Branch("theta_SigmaPScat2npi", &event.theta_SigmaPScat2npi, "theta_SigmaPScat2npi/D");
  tree->Branch("thetaCM_SigmaPScat2npi", &event.thetaCM_SigmaPScat2npi, "thetaCM_SigmaPScat2npi/D");
  tree->Branch("ptrno_SigmaPScat2npi",  &event.ptrno_SigmaPScat2npi , "ptrno_SigmaPScat2npi/I");
  tree->Branch("pitrno_SigmaPScat2npi",  &event.pitrno_SigmaPScat2npi , "pitrno_SigmaPScat2npi/I");
  tree->Branch("pikno_SigmaPScat2npi",  &event.pikno_SigmaPScat2npi , "pikno_SigmaPScat2npi/I");
  tree->Branch("priority_SigmaPScat2npi", &event.priority_SigmaPScat2npi , "priority_SigmaPScat2npi/I");

  tree->Branch("DeltaE_SigmaPScat2ppi", &event.DeltaE_SigmaPScat2ppi, "DeltaE_SigmaPScat2ppi/D");
  tree->Branch("DeltaE2_SigmaPScat2ppi", &event.DeltaE2_SigmaPScat2ppi, "DeltaE2_SigmaPScat2ppi/D");
  tree->Branch("theta2proton_SigmaPScat2ppi", &event.theta2proton_SigmaPScat2ppi, "theta2proton_SigmaPScat2ppi/D");
  tree->Branch("ProtonMom_SigmaPScat2ppi", &event.ProtonMom_SigmaPScat2ppi, "ProtonMom_SigmaPScat2ppi/D");
  tree->Branch("MissMassSigmaP_SigmaPScat2ppi", &event.MissMassSigmaP_SigmaPScat2ppi, "MissMassSigmaP_SigmaPScat2ppi/D");
  tree->Branch("MissMassSigmaP2_SigmaPScat2ppi", &event.MissMassSigmaP2_SigmaPScat2ppi, "MissMassSigmaP2_SigmaPScat2ppi/D");
  tree->Branch("vtx_ScatSigmaDecay2ppi", &event.vtx_ScatSigmaDecay2ppi, "vtx_ScatSigmaDecay2ppi/D");
  tree->Branch("vty_ScatSigmaDecay2ppi", &event.vty_ScatSigmaDecay2ppi, "vty_ScatSigmaDecay2ppi/D");
  tree->Branch("vtz_ScatSigmaDecay2ppi", &event.vtz_ScatSigmaDecay2ppi, "vtz_ScatSigmaDecay2ppi/D");
  tree->Branch("cdistScatSigmaDecay2ppi", &event.cdistScatSigmaDecay2ppi, "cdistScatSigmaDecay2ppi/D");
  tree->Branch("vdistance_SigmaPScat2ppi", &event.vdistance_SigmaPScat2ppi, "vdistance_SigmaPScat2ppi/D");
  tree->Branch("theta_SigmaPScat2ppi", &event.theta_SigmaPScat2ppi, "theta_SigmaPScat2ppi/D");
  tree->Branch("thetaCM_SigmaPScat2ppi", &event.thetaCM_SigmaPScat2ppi, "thetaCM_SigmaPScat2ppi/D");
  tree->Branch("ptrno_SigmaPScat2ppi",  &event.ptrno_SigmaPScat2ppi , "ptrno_SigmaPScat2ppi/I");
  tree->Branch("p2trno_SigmaPScat2ppi",  &event.p2trno_SigmaPScat2ppi , "p2trno_SigmaPScat2ppi/I");
  tree->Branch("pikno_SigmaPScat2ppi",  &event.pikno_SigmaPScat2ppi , "pikno_SigmaPScat2ppi/I");
  tree->Branch("priority_SigmaPScat2ppi", &event.priority_SigmaPScat2ppi , "priority_SigmaPScat2ppi/I");

  ////////// Bring Address From Dst
  TTreeCont[kHodoscope]->SetBranchStatus("*", 0);
  TTreeCont[kHodoscope]->SetBranchStatus("evnum",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("spill",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("trigflag", 1);
  TTreeCont[kHodoscope]->SetBranchStatus("trigpat",  1);
  TTreeCont[kHodoscope]->SetBranchStatus("nhBh1",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("csBh1",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("Bh1Seg",   1);
  TTreeCont[kHodoscope]->SetBranchStatus("tBh1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtBh1",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("deBh1",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("btof",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("nhBh2",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("csBh2",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("Bh2Seg",   1);
  TTreeCont[kHodoscope]->SetBranchStatus("t0Bh2",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("tBh2",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtBh2",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("deBh2",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("Time0",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("CTime0",   1);
  TTreeCont[kHodoscope]->SetBranchStatus("Time0Seg", 1);
  TTreeCont[kHodoscope]->SetBranchStatus("deTime0",  1);
  TTreeCont[kHodoscope]->SetBranchStatus("nhTof",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("csTof",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("TofSeg",   1);
  TTreeCont[kHodoscope]->SetBranchStatus("tTof",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtTof",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("deTof",    1);

  TTreeCont[kHodoscope]->SetBranchAddress("evnum",    &src.evnum);
  TTreeCont[kHodoscope]->SetBranchAddress("spill",    &src.spill);
  TTreeCont[kHodoscope]->SetBranchAddress("trigflag",  src.trigflag);
  TTreeCont[kHodoscope]->SetBranchAddress("trigpat",   src.trigpat);
  TTreeCont[kHodoscope]->SetBranchAddress("nhBh1",    &src.nhBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("csBh1",     src.csBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("Bh1Seg",    src.Bh1Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("tBh1",      src.tBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("dtBh1",     src.dtBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("deBh1",     src.deBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("btof",      src.btof);
  TTreeCont[kHodoscope]->SetBranchAddress("nhBh2",    &src.nhBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("csBh2",     src.csBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("Bh2Seg",    src.Bh2Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("t0Bh2",     src.t0Bh2);
  TTreeCont[kHodoscope]->SetBranchAddress("tBh2",      src.tBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("dtBh2",     src.dtBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("deBh2",     src.deBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("Time0",    &src.Time0);
  TTreeCont[kHodoscope]->SetBranchAddress("CTime0",   &src.CTime0);
  TTreeCont[kHodoscope]->SetBranchAddress("Time0Seg", &src.Time0Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("deTime0",  &src.deTime0);
  TTreeCont[kHodoscope]->SetBranchAddress("nhTof",    &src.nhTof);
  TTreeCont[kHodoscope]->SetBranchAddress("csTof",     src.csTof);
  TTreeCont[kHodoscope]->SetBranchAddress("TofSeg",    src.TofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("tTof",      src.tTof);
  TTreeCont[kHodoscope]->SetBranchAddress("dtTof",     src.dtTof);
  TTreeCont[kHodoscope]->SetBranchAddress("deTof",     src.deTof);

  TTreeCont[kKuramaTracking]->SetBranchStatus("*", 0);
  TTreeCont[kKuramaTracking]->SetBranchStatus("ntSdcIn",      1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("nlSdcIn",      1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("nhSdcIn",      1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("chisqrSdcIn",  1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("x0SdcIn",      1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("y0SdcIn",      1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("u0SdcIn",      1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("v0SdcIn",      1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("ntSdcOut",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("nhSdcOut",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("chisqrSdcOut", 1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("x0SdcOut",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("y0SdcOut",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("u0SdcOut",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("v0SdcOut",     1);

  TTreeCont[kKuramaTracking]->SetBranchStatus("ntKurama",    1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("nhKurama",    1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("stof",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("path",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("pKurama",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("qKurama",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("chisqrKurama",1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("m2",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("xtgtKurama",  1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("ytgtKurama",  1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("utgtKurama",  1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("vtgtKurama",  1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("thetaKurama", 1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("xtofKurama",  1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("ytofKurama",  1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("utofKurama",  1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("vtofKurama",  1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("tofsegKurama",1);

  TTreeCont[kKuramaTracking]->SetBranchAddress("ntSdcIn",      &src.ntSdcIn      );
  TTreeCont[kKuramaTracking]->SetBranchAddress("nlSdcIn",      &src.nlSdcIn      );
  TTreeCont[kKuramaTracking]->SetBranchAddress("nhSdcIn",       src.nhSdcIn      );
  TTreeCont[kKuramaTracking]->SetBranchAddress("chisqrSdcIn",   src.chisqrSdcIn  );
  TTreeCont[kKuramaTracking]->SetBranchAddress("x0SdcIn",       src.x0SdcIn      );
  TTreeCont[kKuramaTracking]->SetBranchAddress("y0SdcIn",       src.y0SdcIn      );
  TTreeCont[kKuramaTracking]->SetBranchAddress("u0SdcIn",       src.u0SdcIn      );
  TTreeCont[kKuramaTracking]->SetBranchAddress("v0SdcIn",       src.v0SdcIn      );
  TTreeCont[kKuramaTracking]->SetBranchAddress("ntSdcOut",     &src.ntSdcOut     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("nhSdcOut",      src.nhSdcOut     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("chisqrSdcOut",  src.chisqrSdcOut );
  TTreeCont[kKuramaTracking]->SetBranchAddress("x0SdcOut",      src.x0SdcOut     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("y0SdcOut",      src.y0SdcOut     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("u0SdcOut",      src.u0SdcOut     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("v0SdcOut",      src.v0SdcOut     );

  TTreeCont[kKuramaTracking]->SetBranchAddress("ntKurama",    &src.ntKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("nhKurama",     src.nhKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("stof",         src.stof);
  TTreeCont[kKuramaTracking]->SetBranchAddress("path",         src.path);
  TTreeCont[kKuramaTracking]->SetBranchAddress("pKurama",      src.pKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("qKurama",      src.qKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("chisqrKurama", src.chisqrKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("m2",      src.m2);
  TTreeCont[kKuramaTracking]->SetBranchAddress("xtgtKurama",   src.xtgtKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("ytgtKurama",   src.ytgtKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("utgtKurama",   src.utgtKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("vtgtKurama",   src.vtgtKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("thetaKurama",  src.thetaKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("xtofKurama",   src.xtofKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("ytofKurama",   src.ytofKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("utofKurama",   src.utofKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("vtofKurama",   src.vtofKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("tofsegKurama", src.tofsegKurama);

  TTreeCont[kK18Tracking]->SetBranchStatus("*", 0);
  TTreeCont[kK18Tracking]->SetBranchStatus("bft_ncl",    1);
  TTreeCont[kK18Tracking]->SetBranchStatus("ntBcOut",     1);
  TTreeCont[kK18Tracking]->SetBranchStatus("nlBcOut",     1);
  TTreeCont[kK18Tracking]->SetBranchStatus("nhBcOut",     1);
  TTreeCont[kK18Tracking]->SetBranchStatus("chisqrBcOut", 1);
  TTreeCont[kK18Tracking]->SetBranchStatus("x0BcOut",     1);
  TTreeCont[kK18Tracking]->SetBranchStatus("y0BcOut",     1);
  TTreeCont[kK18Tracking]->SetBranchStatus("u0BcOut",     1);
  TTreeCont[kK18Tracking]->SetBranchStatus("v0BcOut",     1);

  TTreeCont[kK18Tracking]->SetBranchStatus("ntK18",       1);
  TTreeCont[kK18Tracking]->SetBranchStatus("nhK18",       1);
  TTreeCont[kK18Tracking]->SetBranchStatus("chisqrK18",   1);
  TTreeCont[kK18Tracking]->SetBranchStatus("p_3rd",       1);
  TTreeCont[kK18Tracking]->SetBranchStatus("xtgtK18",     1);
  TTreeCont[kK18Tracking]->SetBranchStatus("ytgtK18",     1);
  TTreeCont[kK18Tracking]->SetBranchStatus("utgtK18",     1);
  TTreeCont[kK18Tracking]->SetBranchStatus("vtgtK18",     1);

  TTreeCont[kK18Tracking]->SetBranchAddress("bft_ncl",     &src.bft_ncl);
  TTreeCont[kK18Tracking]->SetBranchAddress("ntBcOut",     &src.ntBcOut     );
  TTreeCont[kK18Tracking]->SetBranchAddress("nlBcOut",     &src.nlBcOut     );
  TTreeCont[kK18Tracking]->SetBranchAddress("nhBcOut",      src.nhBcOut     );
  TTreeCont[kK18Tracking]->SetBranchAddress("chisqrBcOut",  src.chisqrBcOut );
  TTreeCont[kK18Tracking]->SetBranchAddress("x0BcOut",      src.x0BcOut     );
  TTreeCont[kK18Tracking]->SetBranchAddress("y0BcOut",      src.y0BcOut     );
  TTreeCont[kK18Tracking]->SetBranchAddress("u0BcOut",      src.u0BcOut     );
  TTreeCont[kK18Tracking]->SetBranchAddress("v0BcOut",      src.v0BcOut     );

  TTreeCont[kK18Tracking]->SetBranchAddress("ntK18", &src.ntK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("nhK18", &src.nhK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("chisqrK18", &src.chisqrK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("p_3rd", &src.pK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("xtgtK18", &src.xtgtK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("ytgtK18", &src.ytgtK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("utgtK18", &src.utgtK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("vtgtK18", &src.vtgtK18);
  /*
  TTreeCont[kEasiroc]->SetBranchStatus("*",    0);
  TTreeCont[kEasiroc]->SetBranchStatus("bft_ncl",    1);
  TTreeCont[kEasiroc]->SetBranchStatus("bft_clsize", 1);
  TTreeCont[kEasiroc]->SetBranchStatus("bft_ctime",  1);
  TTreeCont[kEasiroc]->SetBranchStatus("bft_ctot",   1);
  TTreeCont[kEasiroc]->SetBranchStatus("bft_clpos",  1);
  TTreeCont[kEasiroc]->SetBranchStatus("sch_nhits",  1);
  TTreeCont[kEasiroc]->SetBranchStatus("sch_hitpat", 1);
  TTreeCont[kEasiroc]->SetBranchStatus("sch_ncl",    1);
  TTreeCont[kEasiroc]->SetBranchStatus("sch_clsize", 1);
  TTreeCont[kEasiroc]->SetBranchStatus("sch_ctime",  1);
  TTreeCont[kEasiroc]->SetBranchStatus("sch_ctot",   1);
  TTreeCont[kEasiroc]->SetBranchStatus("sch_clpos",  1);

  TTreeCont[kEasiroc]->SetBranchAddress("bft_ncl",    &src.nhBft);
  TTreeCont[kEasiroc]->SetBranchAddress("bft_clsize", &src.csBft);
  TTreeCont[kEasiroc]->SetBranchAddress("bft_ctime",  &src.tBft);
  TTreeCont[kEasiroc]->SetBranchAddress("bft_ctot",   &src.wBft);
  TTreeCont[kEasiroc]->SetBranchAddress("bft_clpos",  &src.BftPos);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_nhits",  &src.NSch);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_hitpat", &src.SchSeg);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_ncl",    &src.nhSch);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_clsize", &src.csSch);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_ctime",  &src.tSch);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_ctot",   &src.wSch);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_clpos",  &src.SchPos);
  */
  TTreeCont[kCFT]->SetBranchStatus("*",    0);
  TTreeCont[kCFT]->SetBranchStatus("ntCFT",    1);
  TTreeCont[kCFT]->SetBranchStatus("theta",    1);
  TTreeCont[kCFT]->SetBranchStatus("phi",    1);
  TTreeCont[kCFT]->SetBranchStatus("vtx_x",    1);
  TTreeCont[kCFT]->SetBranchStatus("vtx_y",    1);
  TTreeCont[kCFT]->SetBranchStatus("vtx_z",    1);
  TTreeCont[kCFT]->SetBranchStatus("Pos0_x",    1);
  TTreeCont[kCFT]->SetBranchStatus("Pos0_y",    1);
  TTreeCont[kCFT]->SetBranchStatus("Pos0_z",    1);
  TTreeCont[kCFT]->SetBranchStatus("Dir_x",    1);
  TTreeCont[kCFT]->SetBranchStatus("Dir_y",    1);
  TTreeCont[kCFT]->SetBranchStatus("Dir_z",    1);
  TTreeCont[kCFT]->SetBranchStatus("Total_dE",    1);
  TTreeCont[kCFT]->SetBranchStatus("Total_dEphi",    1);
  TTreeCont[kCFT]->SetBranchStatus("Total_dEuv",    1);
  TTreeCont[kCFT]->SetBranchStatus("Total_dE_max",    1);
  TTreeCont[kCFT]->SetBranchStatus("Total_dEphi_max",    1);
  TTreeCont[kCFT]->SetBranchStatus("Total_dEuv_max",    1);
  TTreeCont[kCFT]->SetBranchStatus("nhit_phi",    1);
  TTreeCont[kCFT]->SetBranchStatus("nhit_uv",    1);
  TTreeCont[kCFT]->SetBranchStatus("segBGOt",    1);
  TTreeCont[kCFT]->SetBranchStatus("energyBGO",    1);
  TTreeCont[kCFT]->SetBranchStatus("segPiIDt",    1);

  TTreeCont[kCFT]->SetBranchAddress("ntCFT",       &src.ntCFT);
  TTreeCont[kCFT]->SetBranchAddress("theta",       &src.theta_cft);
  TTreeCont[kCFT]->SetBranchAddress("phi",         &src.phi_cft);
  TTreeCont[kCFT]->SetBranchAddress("vtx_x",       &src.vtx_cft);
  TTreeCont[kCFT]->SetBranchAddress("vtx_y",       &src.vty_cft);
  TTreeCont[kCFT]->SetBranchAddress("vtx_z",       &src.vtz_cft);
  TTreeCont[kCFT]->SetBranchAddress("Pos0_x",      &src.Pos0_x);
  TTreeCont[kCFT]->SetBranchAddress("Pos0_y",      &src.Pos0_y);
  TTreeCont[kCFT]->SetBranchAddress("Pos0_z",      &src.Pos0_z);
  TTreeCont[kCFT]->SetBranchAddress("Dir_x",       &src.Dir_x);
  TTreeCont[kCFT]->SetBranchAddress("Dir_y",       &src.Dir_y);
  TTreeCont[kCFT]->SetBranchAddress("Dir_z",       &src.Dir_z);
  TTreeCont[kCFT]->SetBranchAddress("Total_dE",    &src.Total_dE);
  TTreeCont[kCFT]->SetBranchAddress("Total_dEphi", &src.Total_dEphi);
  TTreeCont[kCFT]->SetBranchAddress("Total_dEuv",  &src.Total_dEuv);
  TTreeCont[kCFT]->SetBranchAddress("Total_dE_max",   &src.Total_dE_max);
  TTreeCont[kCFT]->SetBranchAddress("Total_dEphi_max",&src.Total_dEphi_max);
  TTreeCont[kCFT]->SetBranchAddress("Total_dEuv_max", &src.Total_dEuv_max);
  TTreeCont[kCFT]->SetBranchAddress("nhit_phi",    &src.nhit_phi);
  TTreeCont[kCFT]->SetBranchAddress("nhit_uv",     &src.nhit_uv);
  TTreeCont[kCFT]->SetBranchAddress("segBGOt",     &src.segBGOt);
  TTreeCont[kCFT]->SetBranchAddress("energyBGO",  &src.energyBGO);
  TTreeCont[kCFT]->SetBranchAddress("segPiIDt",    &src.segPiIDt);


  return true;
}

//_____________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")   &&
      InitializeParameter<UserParamMan>("USER") &&
      InitializeParameter<HodoPHCMan>("HDPHC"));
}

//_____________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}

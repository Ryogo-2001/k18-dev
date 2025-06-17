// -*- C++ -*-

#include "DstHelper.hh"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include <filesystem_util.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DatabasePDG.hh"
#include "DCRawHit.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "FiberCluster.hh"
#include "FuncName.hh"
#include "HodoPHCMan.hh"
#include "K18TrackD2U.hh"
#include "KuramaLib.hh"
#include "MathTools.hh"
#include "NuclearMass.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"

namespace
{
using namespace root;
using namespace dst;
const auto qnan = TMath::QuietNaN();
auto& gConf = ConfMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();

//_________________________________________________________
// local functions

//_________________________________________________________
inline const TString&
ClassName()
{
  static TString s_name("DstKKAna");
  return s_name;
}

//_________________________________________________________
[[maybe_unused]] Double_t
CalcCutLineByTOF(Double_t M, Double_t mom)
{
  const Double_t K  = 0.307075;//cosfficient [MeV cm^2/g]
  const Double_t ro = 1.032;//density of TOF [g/cm^3]
  const Double_t ZA = 0.5862;//atomic number to mass number ratio of TOF
  const Int_t    z  = 1;//charge of incident particle
  const Double_t me = 0.511;//mass of electron [MeV/c^2]
  const Double_t I  = 0.0000647;//mean excitaton potential [MeV]

  const Double_t C0 = -3.20;
  const Double_t a  = 0.1610;
  const Double_t m  = 3.24;
  const Double_t X1 = 2.49;
  const Double_t X0 = 0.1464;

  Double_t p1 = ro*K*ZA*z*z/2;
  Double_t p2 = (2*me/I/M)*(2*me/I/M);
  Double_t p3 = 0.2414*K*z*z/ro;

  Double_t y = mom/M;
  Double_t X = log10(y);
  Double_t delta;
  if(X < X0)              delta = 0;
  else if(X0 < X&&X < X1) delta = 4.6052*X+C0+a*pow(X1-X,m);
  else                      delta = 4.6052*X+C0;

  Double_t C = (0.422377*pow(y,-2)+0.0304043*pow(y,-4)-0.00038106*pow(y,-6))*pow(10,-6)*pow(I,2)
    +(3.85019*pow(y,-2)-0.1667989*pow(y,-4)+0.00157955*pow(y,-6))*pow(10,-9)*pow(I,3);

  Double_t x = mom;
  Double_t dEdx = p1*((M/x)*(M/x)+1)*log(p2*x*x*x*x/(M*M+2*me*sqrt(M*M+x*x))+me*me)-2*p1-p1*delta-p3*C;

  return dEdx;
}
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kKuramaTracking, kK18Tracking, kHodoscope, kEasiroc,
  kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[KuramaTracking]",
  "[K18Tracking]", "[Hodoscope]", "[Easiroc]",
  "[OutFile]" };
std::vector<TString> TreeName =
{ "", "", "kurama", "k18track", "hodo", "ea0c", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;

const Double_t pB_offset   = 1.000;
const Double_t pS_offset   = 1.000;
const Double_t pK18_offset = 0.000;
const Double_t pKURAMA_offset = 0.000;
const Double_t x_off = 0.000;
const Double_t y_off = 0.000;
const Double_t u_off = 0.000;
const Double_t v_off = 0.000;
}

//_____________________________________________________________________________
struct Event
{
  Int_t runnum;
  Int_t evnum;
  Int_t spill;

  //Trigger
  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  //Hodoscope
  Int_t    nhBh1;
  Int_t    csBh1[NumOfSegBH1*MaxDepth];
  Double_t Bh1Seg[NumOfSegBH1*MaxDepth];
  Double_t tBh1[NumOfSegBH1*MaxDepth];
  Double_t dtBh1[NumOfSegBH1*MaxDepth];
  Double_t deBh1[NumOfSegBH1*MaxDepth];

  Int_t    nhBh2;
  Int_t    csBh2[NumOfSegBH2*MaxDepth];
  Double_t Bh2Seg[NumOfSegBH2*MaxDepth];
  Double_t tBh2[NumOfSegBH2*MaxDepth];
  Double_t t0Bh2[NumOfSegBH2*MaxDepth];
  Double_t dtBh2[NumOfSegBH2*MaxDepth];
  Double_t deBh2[NumOfSegBH2*MaxDepth];

  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;
  Double_t Btof0Seg;
  Double_t deBtof0;
  Double_t Btof0;
  Double_t CBtof0;

  Int_t    nhTof;
  Int_t    csTof[NumOfSegTOF*MaxDepth];
  Double_t TofSeg[NumOfSegTOF*MaxDepth];
  Double_t tTof[NumOfSegTOF*MaxDepth];
  Double_t dtTof[NumOfSegTOF*MaxDepth];
  Double_t deTof[NumOfSegTOF*MaxDepth];

  Int_t    nhBvh;
  Double_t BvhSeg[NumOfSegBVH];

  //Fiber
  Int_t    nhBft;
  Int_t    csBft[NumOfSegBFT];
  Double_t tBft[NumOfSegBFT];
  Double_t wBft[NumOfSegBFT];
  Double_t BftPos[NumOfSegBFT];
  Double_t BftSeg[NumOfSegBFT];
  Int_t    nhSch;
  Int_t    csSch[NumOfSegSCH];
  Double_t tSch[NumOfSegSCH];
  Double_t wSch[NumOfSegSCH];
  Double_t SchPos[NumOfSegSCH];
  Double_t SchSeg[NumOfSegSCH];

  //DC Beam
  Int_t ntBcOut;
  Int_t nlBcOut;
  Int_t nhBcOut[MaxHits];
  Double_t chisqrBcOut[MaxHits];
  Double_t x0BcOut[MaxHits];
  Double_t y0BcOut[MaxHits];
  Double_t u0BcOut[MaxHits];
  Double_t v0BcOut[MaxHits];
  Double_t xtgtBcOut[MaxHits];
  Double_t ytgtBcOut[MaxHits];
  Double_t xbh2BcOut[MaxHits];
  Double_t ybh2BcOut[MaxHits];

  Int_t    ntK18;
  Int_t    nhK18[MaxHits];
  Double_t chisqrK18[MaxHits];
  Double_t pK18[MaxHits];
  Double_t xtgtK18[MaxHits];
  Double_t ytgtK18[MaxHits];
  Double_t utgtK18[MaxHits];
  Double_t vtgtK18[MaxHits];
  Double_t thetaK18[MaxHits];

  //DC KURAMA
  Int_t ntSdcIn;
  Int_t nlSdcIn;
  Int_t nhSdcIn[MaxHits];
  Double_t chisqrSdcIn[MaxHits];
  Double_t x0SdcIn[MaxHits];
  Double_t y0SdcIn[MaxHits];
  Double_t u0SdcIn[MaxHits];
  Double_t v0SdcIn[MaxHits];

  Int_t ntSdcOut;
  Int_t nlSdcOut;
  Int_t nhSdcOut[MaxHits];
  Double_t chisqrSdcOut[MaxHits];
  Double_t u0SdcOut[MaxHits];
  Double_t v0SdcOut[MaxHits];
  Double_t x0SdcOut[MaxHits];
  Double_t y0SdcOut[MaxHits];

  Int_t    ntKurama;
  Int_t    nhKurama[MaxHits];
  Double_t chisqrKurama[MaxHits];
  Double_t stof[MaxHits];
  Double_t cstof[MaxHits];
  Double_t path[MaxHits];
  Double_t pKurama[MaxHits];
  Double_t qKurama[MaxHits];
  Double_t m2[MaxHits];
  Double_t xtgtKurama[MaxHits];
  Double_t ytgtKurama[MaxHits];
  Double_t utgtKurama[MaxHits];
  Double_t vtgtKurama[MaxHits];
  Double_t thetaKurama[MaxHits];
  Double_t xtofKurama[MaxHits];
  Double_t ytofKurama[MaxHits];
  Double_t utofKurama[MaxHits];
  Double_t vtofKurama[MaxHits];
  Double_t tofsegKurama[MaxHits];
  Double_t best_deTof[MaxHits];
  Double_t best_TofSeg[MaxHits];

  //Reaction
  Int_t    nKm;
  Int_t    nKp;
  Int_t    nKK;
  Double_t vtx[MaxHits];
  Double_t vty[MaxHits];
  Double_t vtz[MaxHits];
  Double_t closeDist[MaxHits];
  Int_t    inside[MaxHits];
  Double_t theta[MaxHits];
  Double_t MissMass[MaxHits];
  Double_t MissMassCorr[MaxHits];
  Double_t MissMassCorrDE[MaxHits];
  Double_t thetaCM[MaxHits];
  Double_t costCM[MaxHits];
  Int_t Kflag[MaxHits];

  Double_t xkm[MaxHits];
  Double_t ykm[MaxHits];
  Double_t ukm[MaxHits];
  Double_t vkm[MaxHits];
  Double_t xkp[MaxHits];
  Double_t ykp[MaxHits];
  Double_t ukp[MaxHits];
  Double_t vkp[MaxHits];

  Double_t pOrg[MaxHits];
  Double_t pCalc[MaxHits];
  Double_t pCorr[MaxHits];
  Double_t pCorrDE[MaxHits];
};

//_____________________________________________________________________________
struct Src
{
  Int_t runnum;
  Int_t evnum;
  Int_t spill;

  //Trigger
  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  //Hodoscope
  Int_t    nhBh1;
  Int_t    csBh1[NumOfSegBH1*MaxDepth];
  Double_t Bh1Seg[NumOfSegBH1*MaxDepth];
  Double_t tBh1[NumOfSegBH1*MaxDepth];
  Double_t dtBh1[NumOfSegBH1*MaxDepth];
  Double_t deBh1[NumOfSegBH1*MaxDepth];

  Int_t    nhBh2;
  Int_t    csBh2[NumOfSegBH2*MaxDepth];
  Double_t Bh2Seg[NumOfSegBH2*MaxDepth];
  Double_t tBh2[NumOfSegBH2*MaxDepth];
  Double_t t0Bh2[NumOfSegBH2*MaxDepth];
  Double_t dtBh2[NumOfSegBH2*MaxDepth];
  Double_t deBh2[NumOfSegBH2*MaxDepth];

  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;
  Double_t Btof0Seg;
  Double_t deBtof0;
  Double_t Btof0;
  Double_t CBtof0;

  Int_t    nhTof;
  Int_t    csTof[NumOfSegTOF*MaxDepth];
  Double_t TofSeg[NumOfSegTOF*MaxDepth];
  Double_t tTof[NumOfSegTOF*MaxDepth];
  Double_t dtTof[NumOfSegTOF*MaxDepth];
  Double_t deTof[NumOfSegTOF*MaxDepth];

  Int_t    nhBvh;
  Double_t BvhSeg[NumOfSegBVH];

  //Fiber
  Int_t    nhBft;
  Int_t    csBft[NumOfSegBFT];
  Double_t tBft[NumOfSegBFT];
  Double_t wBft[NumOfSegBFT];
  Double_t BftPos[NumOfSegBFT];
  Double_t BftSeg[NumOfSegBFT];
  Int_t    nhSch;
  // Int_t    sch_hitpat[NumOfSegSCH];
  Int_t    csSch[NumOfSegSCH];
  Double_t tSch[NumOfSegSCH];
  Double_t wSch[NumOfSegSCH];
  Double_t SchPos[NumOfSegSCH];
  Double_t SchSeg[NumOfSegSCH];

  //DC Beam
  Int_t ntBcOut;
  Int_t nlBcOut;
  Int_t nhBcOut[MaxHits];
  Double_t chisqrBcOut[MaxHits];
  Double_t x0BcOut[MaxHits];
  Double_t y0BcOut[MaxHits];
  Double_t u0BcOut[MaxHits];
  Double_t v0BcOut[MaxHits];
  Double_t xtgtBcOut[MaxHits];
  Double_t ytgtBcOut[MaxHits];
  Double_t xbh2BcOut[MaxHits];
  Double_t ybh2BcOut[MaxHits];

  Int_t    ntK18;
  Int_t    nhK18[MaxHits];
  Double_t chisqrK18[MaxHits];
  Double_t pK18[MaxHits];
  Double_t xtgtK18[MaxHits];
  Double_t ytgtK18[MaxHits];
  Double_t utgtK18[MaxHits];
  Double_t vtgtK18[MaxHits];
  Double_t thetaK18[MaxHits];

  //DC KURAMA
  Int_t much;
  Int_t ntSdcIn;
  Int_t nlSdcIn;
  Int_t nhSdcIn[MaxHits];
  Double_t chisqrSdcIn[MaxHits];
  Double_t x0SdcIn[MaxHits];
  Double_t y0SdcIn[MaxHits];
  Double_t u0SdcIn[MaxHits];
  Double_t v0SdcIn[MaxHits];

  Int_t ntSdcOut;
  Int_t nlSdcOut;
  Int_t nhSdcOut[MaxHits];
  Double_t chisqrSdcOut[MaxHits];
  Double_t u0SdcOut[MaxHits];
  Double_t v0SdcOut[MaxHits];
  Double_t x0SdcOut[MaxHits];
  Double_t y0SdcOut[MaxHits];

  Int_t    ntKurama;
  Int_t    nhKurama[MaxHits];
  Double_t chisqrKurama[MaxHits];
  Double_t stof[MaxHits];
  Double_t cstof[MaxHits];
  Double_t path[MaxHits];
  Double_t pKurama[MaxHits];
  Double_t qKurama[MaxHits];
  Double_t m2[MaxHits];
  Double_t xtgtKurama[MaxHits];
  Double_t ytgtKurama[MaxHits];
  Double_t utgtKurama[MaxHits];
  Double_t vtgtKurama[MaxHits];
  Double_t thetaKurama[MaxHits];
  Double_t xtofKurama[MaxHits];
  Double_t ytofKurama[MaxHits];
  Double_t utofKurama[MaxHits];
  Double_t vtofKurama[MaxHits];
  Double_t tofsegKurama[MaxHits];
};

//_____________________________________________________________________________
namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
}

//_____________________________________________________________________________
Int_t
main(Int_t argc, char **argv)
{
  std::vector<std::string> arg(argv, argv+argc);

  if(!CheckArg(arg))
    return EXIT_FAILURE;
  if(!DstOpen(arg))
    return EXIT_FAILURE;
  if(!gConf.Initialize(arg[kConfFile]))
    return EXIT_FAILURE;

  Int_t nevent = GetEntries(TTreeCont);

  CatchSignal::Set();

  Int_t ievent = 0;
  for(; ievent<nevent && !CatchSignal::Stop(); ++ievent){
    gCounter.check();
    InitializeEvent();
    if(DstRead(ievent)) tree->Fill();
  }

  std::cout << "#D Event Number: " << std::setw(6)
	    << ievent << std::endl;

  DstClose();

  return EXIT_SUCCESS;
}

//_____________________________________________________________________________
bool
dst::InitializeEvent()
{
  event.runnum     = 0;
  event.evnum      = 0;
  event.spill      = 0;
  event.Time0Seg = qnan;
  event.deTime0  = qnan;
  event.Time0    = qnan;
  event.CTime0   = qnan;
  event.Btof0Seg = qnan;
  event.deBtof0  = qnan;
  event.Btof0    = qnan;
  event.CBtof0   = qnan;

  //Trigger
  for(Int_t it=0; it<NumOfSegTrig; ++it){
    event.trigpat[it]  = -1;
    event.trigflag[it] = -1;
  }

  //Hodoscope
  event.nhBh2  = 0;
  event.nhBh1  = 0;
  event.nhTof  = 0;
  event.nhBvh  = 0;

  for(Int_t i=0; i<NumOfSegBH1; ++i){
    for(Int_t j=0; j<MaxDepth; ++j){
      event.csBh1[i*MaxDepth+j]  = 0;
      event.Bh1Seg[i*MaxDepth+j] = -1;
      event.tBh1[i*MaxDepth+j]   = qnan;
      event.dtBh1[i*MaxDepth+j]  = qnan;
      event.deBh1[i*MaxDepth+j]  = qnan;
    }
  }
  for(Int_t i=0; i<NumOfSegBH2; ++i){
    for(Int_t j=0; j<MaxDepth; ++j){
      event.csBh2[i*MaxDepth+j]  = 0;
      event.Bh2Seg[i*MaxDepth+j] = -1;
      event.tBh2[i*MaxDepth+j]   = qnan;
      event.t0Bh2[i*MaxDepth+j]  = qnan;
      event.dtBh2[i*MaxDepth+j]  = qnan;
      event.deBh2[i*MaxDepth+j]  = qnan;
    }
  }
  for(Int_t i=0; i<NumOfSegTOF; ++i){
    for(Int_t j=0; j<MaxDepth; ++j){
      event.csTof[i*MaxDepth+j]  = 0;
      event.TofSeg[i*MaxDepth+j] = -1;
      event.tTof[i*MaxDepth+j]   = qnan;
      event.dtTof[i*MaxDepth+j]  = qnan;
      event.deTof[i*MaxDepth+j]  = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegBVH; it++){
    event.BvhSeg[it] = -1;
  }

  //Fiber
  event.nhBft = 0;
  event.nhSch = 0;
  for(Int_t it=0; it<NumOfSegBFT; it++){
    event.csBft[it] = 0;
    event.tBft[it] = qnan;
    event.wBft[it] = qnan;
    event.BftPos[it] = qnan;
    event.BftSeg[it] = qnan;
  }
  for(Int_t it=0; it<NumOfSegSCH; it++){
    event.csSch[it] = 0;
    event.tSch[it] = qnan;
    event.wSch[it] = qnan;
    event.SchPos[it] = qnan;
    event.SchSeg[it] = qnan;
  }

  //DC
  event.nlBcOut  = 0;
  event.nlSdcIn  = 0;
  event.nlSdcOut = 0;
  event.ntBcOut  = 0;
  event.ntSdcIn  = 0;
  event.ntSdcOut = 0;
  event.ntK18    = 0;
  event.ntKurama = 0;

  //Beam DC
  for(Int_t it=0; it<MaxHits; it++){
    event.nhBcOut[it]     = 0;
    event.chisqrBcOut[it] = -1.0;
    event.x0BcOut[it] = qnan;
    event.y0BcOut[it] = qnan;
    event.u0BcOut[it] = qnan;
    event.v0BcOut[it] = qnan;

    event.xtgtBcOut[it] = qnan;
    event.ytgtBcOut[it] = qnan;
    event.xbh2BcOut[it] = qnan;
    event.ybh2BcOut[it] = qnan;
  }

  for(Int_t it=0; it<MaxHits; it++){
    event.nhK18[it]     = 0;
    event.chisqrK18[it] = -1.0;
    event.xtgtK18[it] = qnan;
    event.ytgtK18[it] = qnan;
    event.utgtK18[it] = qnan;
    event.vtgtK18[it] = qnan;
    event.pK18[it]    = qnan;
    event.thetaK18[it] = qnan;
  }

  //KURAMA DC
  for(Int_t it=0; it<MaxHits; it++){
    event.nhSdcIn[it]     = 0;
    event.chisqrSdcIn[it] = -1.;
    event.x0SdcIn[it] = qnan;
    event.y0SdcIn[it] = qnan;
    event.u0SdcIn[it] = qnan;
    event.v0SdcIn[it] = qnan;

    event.nhSdcOut[it]     = 0;
    event.chisqrSdcOut[it] = -1.;
    event.u0SdcOut[it] = qnan;
    event.v0SdcOut[it] = qnan;
    event.x0SdcOut[it] = qnan;
    event.y0SdcOut[it] = qnan;

    event.nhKurama[it]     = 0;
    event.chisqrKurama[it] = -1.;
    event.xtgtKurama[it]   = qnan;
    event.ytgtKurama[it]   = qnan;
    event.utgtKurama[it]   = qnan;
    event.vtgtKurama[it]   = qnan;
    event.pKurama[it]      = qnan;
    event.qKurama[it]      = qnan;
    event.stof[it]         = qnan;
    event.cstof[it]        = qnan;
    event.path[it]         = qnan;
    event.m2[it]           = qnan;
    event.thetaKurama[it]  = qnan;
    event.xtofKurama[it]   = qnan;
    event.ytofKurama[it]   = qnan;
    event.utofKurama[it]   = qnan;
    event.vtofKurama[it]   = qnan;
    event.tofsegKurama[it] = qnan;
    event.best_deTof[it]   = qnan;
    event.best_TofSeg[it]  = qnan;
  }

  //Reaction
  event.nKm = 0;
  event.nKp = 0;
  event.nKK = 0;

  for(Int_t it=0; it<MaxHits; ++it){
    event.vtx[it] = qnan;
    event.vty[it] = qnan;
    event.vtz[it] = qnan;
    event.closeDist[it] = qnan;
    event.inside[it] = 0;
    event.theta[it] = qnan;
    event.thetaCM[it] = qnan;
    event.costCM[it] = qnan;
    event.MissMass[it]  = qnan;
    event.MissMassCorr[it]  = qnan;
    event.MissMassCorrDE[it]  = qnan;
    event.Kflag[it] = 0;

    event.xkm[it] = qnan;
    event.ykm[it] = qnan;
    event.ukm[it] = qnan;
    event.vkm[it] = qnan;
    event.xkp[it] = qnan;
    event.ykp[it] = qnan;
    event.ukp[it] = qnan;
    event.vkp[it] = qnan;
    event.pOrg[it] = qnan;
    event.pCalc[it] = qnan;
    event.pCorr[it] = qnan;
    event.pCorrDE[it] = qnan;
  }

  return true;
}

//_____________________________________________________________________________
bool
dst::DstOpen(std::vector<std::string> arg)
{
  Int_t open_file = 0;
  Int_t open_tree = 0;
  for(std::size_t i=0; i<nArgc; ++i){
    if(i==kProcess || i==kConfFile || i==kOutFile) continue;
    open_file += OpenFile(TFileCont[i], arg[i]);
    open_tree += OpenTree(TFileCont[i], TTreeCont[i], TreeName[i]);
  }

  if(open_file!=open_tree || open_file!=nArgc-3)
    return false;
  if(!CheckEntries(TTreeCont))
    return false;

  TFileCont[kOutFile] = new TFile(arg[kOutFile].c_str(), "recreate");

  return true;
}

//_____________________________________________________________________________
bool
dst::DstRead(Int_t ievent)
{
  static const auto StofOffset = gUser.GetParameter("StofOffset");
  static const auto KaonMass    = pdg::KaonMass();
  // static const auto PionMass    = pdg::PionMass();
  static const auto ProtonMass  = pdg::ProtonMass();
  static const auto SigmaNMass  = pdg::SigmaNMass();

  if(ievent%10000 == 0){
    std::cout << FUNC_NAME << " Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  event.runnum   = src.runnum;
  event.evnum    = src.evnum;
  event.spill    = src.spill;

  event.ntBcOut  = src.ntBcOut;
  event.ntSdcIn  = src.ntSdcIn;
  event.ntSdcOut = src.ntSdcOut;
  event.ntKurama = src.ntKurama;
  event.ntK18    = src.ntK18;
  event.nhBft    = src.nhBft;
  event.nhBh1    = src.nhBh1;
  event.nhBh2    = src.nhBh2;
  event.nhSch    = src.nhSch;
  event.nhTof    = src.nhTof;
  event.nhBvh    = src.nhBvh;

  const Int_t ntBcOut  = event.ntBcOut;
  const Int_t ntSdcIn  = event.ntSdcIn;
  const Int_t ntSdcOut = event.ntSdcOut;
  const Int_t ntKurama = event.ntKurama;
  const Int_t ntK18    = event.ntK18;
  const Int_t nhBft    = event.nhBft;
  const Int_t nhBh1    = event.nhBh1;
  const Int_t nhBh2    = event.nhBh2;
  const Int_t nhSch    = event.nhSch;
  const Int_t nhTof    = event.nhTof;
  const Int_t nhBvh    = event.nhBvh;

#if 0
  std::cout << "#D DebugPrint" << std::endl
	    << " event  : " << std::setw(6) << ievent   << std::endl
	    << " BcOut  : " << std::setw(6) << ntBcOut  << std::endl
	    << " SdcIn  : " << std::setw(6) << ntSdcIn  << std::endl
	    << " SdcOut : " << std::setw(6) << ntSdcOut << std::endl
	    << " Kurama : " << std::setw(6) << ntKurama << std::endl
	    << " K18    : " << std::setw(6) << ntK18    << std::endl
	    << " BFT    : " << std::setw(6) << nhBft    << std::endl
	    << " BH1    : " << std::setw(6) << nhBh1    << std::endl
	    << " BH2    : " << std::setw(6) << nhBh2    << std::endl
	    << " SCH    : " << std::setw(6) << nhSch    << std::endl
	    << " TOF    : " << std::setw(6) << nhTof    << std::endl
	    << std::endl;
#endif

  // TFlag
  for(Int_t i=0;i<NumOfSegTrig;++i){
    Int_t tdc = src.trigflag[i];
    if(tdc<=0) continue;
    event.trigpat[i]  = i + 1;
    event.trigflag[i] = tdc;
  }

  HF1(1, 0.);

  // if(event.ntKurama==0) return true;
  // HF1(1, 1.);
  // if(event.ntK18==0) return true;
  // HF1(1, 2.);
  // if(event.nhBh1==0) return true;
  // HF1(1, 3.);
  // if(event.nhBh2==0) return true;
  // HF1(1, 4.);
  // if(event.nhTof==0) return true;
  // HF1(1, 5.);
  // if(event.nhSch==0) return true;
  // HF1(1, 6.);
  // if(event.nhFbh==0) return true;
  // HF1(1, 7.);

  std::vector<ThreeVector> KmPCont, KmXCont;
  std::vector<ThreeVector> KpPCont, KpXCont;

  // BFT
  for(Int_t i=0; i<nhBft; ++i){
    event.csBft[i]  = src.csBft[i];
    event.tBft[i]   = src.tBft[i];
    event.wBft[i]   = src.wBft[i];
    event.BftPos[i] = src.BftPos[i];
    event.BftSeg[i] = src.BftSeg[i];
  }

  // BH1
  for(Int_t i=0; i<nhBh1; ++i){
    event.csBh1[i]  = src.csBh1[i];
    event.Bh1Seg[i] = src.Bh1Seg[i];
    event.tBh1[i]   = src.tBh1[i];
    event.dtBh1[i]  = src.dtBh1[i];
    event.deBh1[i]  = src.deBh1[i];
  }

  // BH2
  for(Int_t i=0; i<nhBh2; ++i){
    event.csBh2[i]  = src.csBh2[i];
    event.Bh2Seg[i] = src.Bh2Seg[i];
    event.tBh2[i]   = src.tBh2[i];
    event.t0Bh2[i]  = src.t0Bh2[i];
    event.dtBh2[i]  = src.dtBh2[i];
    event.deBh2[i]  = src.deBh2[i];
  }

  event.Time0Seg = src.Time0Seg;
  event.deTime0  = src.deTime0;
  event.Time0 = src.Time0;
  event.CTime0 = src.CTime0;
  event.Btof0Seg = src.Btof0Seg;
  event.deBtof0 = src.deBtof0;
  event.Btof0 = src.Btof0;
  event.CBtof0 = src.CBtof0;

  // SCH
  for(Int_t i=0; i<nhSch; ++i){
    event.csSch[i]  = src.csSch[i];
    event.tSch[i]   = src.tSch[i];
    event.wSch[i]   = src.wSch[i];
    event.SchPos[i] = src.SchPos[i];
    event.SchSeg[i] = src.SchSeg[i];
  }

  // TOF
  for(Int_t i=0; i<nhTof; ++i){
    event.csTof[i]  = src.csTof[i];
    event.TofSeg[i] = src.TofSeg[i];
    event.tTof[i]   = src.tTof[i];
    event.dtTof[i]  = src.dtTof[i];
    event.deTof[i]  = src.deTof[i];
  }

  // BVH
  for(Int_t i=0; i<nhBvh; ++i){
    event.BvhSeg[i] = src.BvhSeg[i];
  }

  ////////// BcOut
  event.nlBcOut = src.nlBcOut;
  for(Int_t it=0; it<ntBcOut; ++it){
    event.nhBcOut[it]     = src.nhBcOut[it];
    event.chisqrBcOut[it] = src.chisqrBcOut[it];
    event.x0BcOut[it]     = src.x0BcOut[it];
    event.y0BcOut[it]     = src.y0BcOut[it];
    event.u0BcOut[it]     = src.u0BcOut[it];
    event.v0BcOut[it]     = src.v0BcOut[it];
  }

  ////////// SdcIn
  event.nlSdcIn = src.nlSdcIn;
  for(Int_t it=0; it<ntSdcIn; ++it){
    event.nhSdcIn[it]     = src.nhSdcIn[it];
    event.chisqrSdcIn[it] = src.chisqrSdcIn[it];
    event.x0SdcIn[it]     = src.x0SdcIn[it];
    event.y0SdcIn[it]     = src.y0SdcIn[it];
    event.u0SdcIn[it]     = src.u0SdcIn[it];
    event.v0SdcIn[it]     = src.v0SdcIn[it];
  }

  ////////// SdcOut
  event.nlSdcOut = src.nlSdcOut;
  for(Int_t it=0; it<ntSdcOut; ++it){
    event.nhSdcOut[it]     = src.nhSdcOut[it];
    event.chisqrSdcOut[it] = src.chisqrSdcOut[it];
    event.x0SdcOut[it]     = src.x0SdcOut[it];
    event.y0SdcOut[it]     = src.y0SdcOut[it];
    event.u0SdcOut[it]     = src.u0SdcOut[it];
    event.v0SdcOut[it]     = src.v0SdcOut[it];
  }

  ////////// K+
  for(Int_t itKurama=0; itKurama<ntKurama; ++itKurama){
    Int_t nh= src.nhKurama[itKurama];
    Double_t chisqr = src.chisqrKurama[itKurama];
    Double_t p = src.pKurama[itKurama];
    Double_t x = src.xtgtKurama[itKurama];
    Double_t y = src.ytgtKurama[itKurama];
    Double_t u = src.utgtKurama[itKurama];
    Double_t v = src.vtgtKurama[itKurama];
    // Double_t utof = src.utofKurama[itKurama];
    // Double_t vtof = src.vtofKurama[itKurama];
    Double_t theta = src.thetaKurama[itKurama];
    Double_t pt = p/TMath::Sqrt(1.+u*u+v*v);
    ThreeVector Pos(x, y, 0.);
    ThreeVector Mom(pt*u, pt*v, pt);
    // if(TMath::IsNaN(Pos.Mag())) continue;
    //Calibration
    ThreeVector PosCorr(x+x_off, y+y_off, 0.);
    ThreeVector MomCorr(pt*(u+u_off), pt*(v+v_off), pt);
    Double_t path = src.path[itKurama];
    Double_t xt = PosCorr.x();
    Double_t yt = PosCorr.y();
    Double_t zt = PosCorr.z();
    Double_t pCorr = MomCorr.Mag();
    Double_t q  = src.qKurama[itKurama];
    Double_t ut = xt/zt;
    Double_t vt = yt/zt;

    Double_t tofsegKurama = src.tofsegKurama[itKurama];
    Double_t stof = qnan;
    Double_t cstof = qnan;
    Double_t m2   = qnan;

    // bool correct_hit[src.nhTof];
    // for(Int_t j=0; j<src.nhTof; ++j) correct_hit[j]=true;
    // for(Int_t j=0; j<src.nhTof; ++j){
    //   Double_t seg1=src.TofSeg[j]; Double_t de1=src.deTof[j];
    //   for(Int_t k=j+1; k<src.nhTof; ++k){
    //     Double_t seg2=src.TofSeg[k]; Double_t de2=src.deTof[k];
    //     if(std::abs(seg1 - seg2) < 2 && de1 > de2) correct_hit[k]=false;
    //     if(std::abs(seg1 - seg2) < 2 && de2 > de1) correct_hit[j]=false;
    //   }
    // }
    // Double_t best_de=-9999;
    // Int_t correct_num=0;
    // Int_t Dif=9999;
    // for(Int_t j=0; j<src.nhTof; ++j){
    //   if(correct_hit[j]==false) continue;
    //   Int_t dif = std::abs(src.TofSeg[j] - tofsegKurama);
    //   if((dif < Dif) || (dif==Dif&&src.deTof[j]>best_de)){
    //     correct_num=j;
    //     best_de=src.deTof[j];
    //     Dif=dif;
    //   }
    // }

    for(Int_t j=0; j<src.nhTof; ++j){
      if(tofsegKurama == src.TofSeg[j]){
        stof = event.tTof[j] - event.CTime0 + StofOffset;
        break;
      }
    }
    m2 = Kinematics::MassSquare(pCorr, path, stof);
    // m2 = Kinematics::MassSquare(pCorr, path, cstof);
    // if(TMath::IsNaN(btof)){
    //   cstof=stof;
    // }else{
    //   gPHC.DoStofCorrection(8, 0, src.TofSeg[correct_num]-1, 2, stof, btof, cstof);
    //   m2 = Kinematics::MassSquare(pCorr, path, cstof);
    // }
    // event.best_deTof[itKurama] = best_de;
    // event.best_TofSeg[itKurama] = src.TofSeg[correct_num];

    ///for Kflag///
    Int_t Kflag=0;
    // Double_t dEdx = Mip2MeV*best_de/sqrt(1+utof*utof+vtof*vtof);
    // if(CalcCutLineByTOF(PionCutMass, 1000*p) < dEdx &&
    //    dEdx < CalcCutLineByTOF(ProtonCutMass, 1000*p)){
    //   Kflag=1;
    // }
    // w/o TOF
    // Double_t minres = 1.0e10;
    // for(Int_t j=0; j<nhTof; ++j){
    //   Double_t seg = src.TofSeg[j];
    //   Double_t res = TofSegKurama - seg;
    //   if(std::abs(res) < minres && std::abs(res) < 1.){
    // 	minres = res;
    // 	stof = src.tTof[j] - time0 + OffsetToF;
    //   }
    // }
    event.nhKurama[itKurama] = nh;
    event.chisqrKurama[itKurama] = chisqr;
    event.pKurama[itKurama] = p;
    event.qKurama[itKurama] = q;
    event.xtgtKurama[itKurama] = x;
    event.ytgtKurama[itKurama] = y;
    event.utgtKurama[itKurama] = u;
    event.vtgtKurama[itKurama] = v;
    event.tofsegKurama[itKurama] = tofsegKurama;
    event.thetaKurama[itKurama] = theta;
    event.stof[itKurama] = stof;
    event.cstof[itKurama] = cstof;
    event.path[itKurama] = path;
    event.m2[itKurama] = m2;
    event.Kflag[itKurama] = Kflag;
    HF1(3202, Double_t(nh));
    HF1(3203, chisqr);
    HF1(3204, xt); HF1(3205, yt);
    HF1(3206, ut); HF1(3207, vt);
    HF2(3208, xt, ut); HF2(3209, yt, vt);
    HF2(3210, xt, yt);
    HF1(3211, pCorr);
    HF1(3212, path);

    //HF1(3213, m2);

    //     Double_t xTof=(clTof->MeanSeg()-7.5)*70.;
    //     Double_t yTof=(clTof->TimeDif())*800./12.;
    //     HF2(4011, clTof->MeanSeg()-0.5, xtof);
    //     HF2(4013, clTof->TimeDif(), ytof);
    //     HF1(4015, xtof-xTof); HF1(4017, ytof-yTof);
    //     HF2(4019, xtof-xTof, ytof-yTof);

    //       Double_t ttof=clTof->CMeanTime()-time0,
    //       HF1(322, clTof->ClusterSize());
    //       HF1(323, clTof->MeanSeg()-0.5);
    //       HF1(324, ttof);
    //       HF1(325, clTof->DeltaE());
    //       Double_t u0in=trIn->GetU0();
    //       HF2(4001, u0in, ttof); HF2(4003, u0in, ttof+12.5*u0in);
    HF1(4202, Double_t(nh));
    HF1(4203, chisqr);
    HF1(4204, xt); HF1(4205, yt);
    HF1(4206, ut); HF1(4207, vt);
    HF2(4208, xt, ut); HF2(4209, yt, vt);
    HF2(4210, xt, yt);
    HF1(4211, pCorr);
    HF1(4212, path);
    KpXCont.push_back(PosCorr);
    KpPCont.push_back(MomCorr);
  }

  if(KpPCont.size()==0) return true;

  HF1(1, 8.);

  ////////// km
  for(Int_t itK18=0; itK18<ntK18; ++itK18){
    Int_t nh = src.nhK18[itK18];
    Double_t chisqr = src.chisqrK18[itK18];
    //Calibration
    Double_t p = src.pK18[itK18]*pB_offset+pK18_offset;
    Double_t x = src.xtgtK18[itK18];
    Double_t y = src.ytgtK18[itK18];
    Double_t u = src.utgtK18[itK18];
    Double_t v = src.vtgtK18[itK18];
    event.nhK18[itK18]     = nh;
    event.chisqrK18[itK18] = chisqr;
    event.pK18[itK18]      = p;
    event.xtgtK18[itK18]   = x;
    event.ytgtK18[itK18]   = y;
    event.utgtK18[itK18]   = u;
    event.vtgtK18[itK18]   = v;
    // Double_t loss_bh2 = 1.09392e-3;
    // p = p - loss_bh2;
    Double_t pt=p/std::sqrt(1.+u*u+v*v);
    ThreeVector Pos(x, y, 0.);
    ThreeVector Mom(pt*u, pt*v, pt);
    //Double_t xo=trOut->GetX0(), yo=trOut->GetY0();
    HF1(4104, p);
    HF1(4105, x); HF1(4106, y);
    //HF1(4107, xo); HF1(4108, yo); HF1(4109, u); HF1(4110, v);

    KmPCont.push_back(Mom); KmXCont.push_back(Pos);
  }

  if(KmPCont.size()==0) return true;

  HF1(1, 9.);

  //MissingMass
  Int_t nKm = KmPCont.size();
  Int_t nKp = KpPCont.size();
  event.nKm = nKm;
  event.nKp = nKp;
  event.nKK = nKm*nKp;
  HF1(4101, Double_t(nKm));
  HF1(4201, Double_t(nKp));
  Int_t nkk=0;
  for(Int_t ikp=0; ikp<nKp; ++ikp){
    ThreeVector pkp = KpPCont[ikp], xkp = KpXCont[ikp];
    for(Int_t ikm=0; ikm<nKm; ++ikm){
      ThreeVector pkm  = KmPCont[ikm], xkm = KmXCont[ikm];
      ThreeVector vert = Kinematics::VertexPoint(xkm, xkp, pkm, pkp);
      // std::cout << "vertex : " << vert << " " << vert.Mag() << std::endl;
      Double_t closedist = Kinematics::CloseDist(xkm, xkp, pkm, pkp);

      Double_t us = pkp.x()/pkp.z(), vs = pkp.y()/pkp.z();
      Double_t ub = pkm.x()/pkm.z(), vb = pkm.y()/pkm.z();
      Double_t cost = pkm*pkp/(pkm.Mag()*pkp.Mag());
      Double_t theta = TMath::ACos(cost)*TMath::RadToDeg();
      Double_t pk0   = pkp.Mag();
      Double_t pCorr = pk0;

      ThreeVector pkpCorr(pCorr*pkp.x()/pkp.Mag(),
                          pCorr*pkp.y()/pkp.Mag(),
                          pCorr*pkp.z()/pkp.Mag());

      //      ThreeVector pkmCorrDE = Kinematics::CorrElossIn(pkm, xkm, vert, KmonMass);
      //      ThreeVector pkpCorrDE = Kinematics::CorrElossOut(pkpCorr, xkp, vert, KaonMass);

      LorentzVector LvKm(pkm, std::sqrt(KaonMass*KaonMass+pkm.Mag2()));
      //      LorentzVector LvKmCorrDE(pkmCorrDE, sqrt(KaonMass*KaonMass+pkmCorrDE.Mag2()));

      LorentzVector LvKp(pkp, std::sqrt(KaonMass*KaonMass+pkp.Mag2()));
      LorentzVector LvKpCorr(pkpCorr, std::sqrt(KaonMass*KaonMass+pkpCorr.Mag2()));
      //      LorentzVector LvKpCorrDE(pkpCorrDE, std::sqrt(KaonMass*KaonMass+pkpCorrDE.Mag2()));

      LorentzVector LvC(0., 0., 0., ProtonMass);
      LorentzVector LvCore(0., 0., 0., 0.);

      LorentzVector LvRc       = LvKm+LvC-LvKp;
      LorentzVector LvRcCorr   = LvKm+LvC-LvKpCorr;
      //      LorentzVector LvRcCorrDE = LvKmCorrDE+LvC-LvKpCorrDE;
      Double_t MissMass       = LvRc.Mag();//-LvC.Mag();
      Double_t MissMassCorr   = LvRcCorr.Mag();//-LvC.Mag();
      //      Double_t MissMassCorrDE = LvRcCorrDE.Mag();//-LvC.Mag();

      //Primary frame
      LorentzVector PrimaryLv = LvKm+LvC;
      Double_t TotalEnergyCM = PrimaryLv.Mag();
      ThreeVector beta(1/PrimaryLv.E()*PrimaryLv.Vect());

      //CM
      Double_t TotalMomCM
	= 0.5*std::sqrt((TotalEnergyCM*TotalEnergyCM
                         -(KaonMass+SigmaNMass)*(KaonMass+SigmaNMass))
			*(TotalEnergyCM*TotalEnergyCM
                          -(KaonMass-SigmaNMass)*(KaonMass-SigmaNMass)))/TotalEnergyCM;

      Double_t costLab = cost;
      Double_t cottLab = costLab/std::sqrt(1.-costLab*costLab);
      Double_t bt=beta.Mag(), gamma=1./std::sqrt(1.-bt*bt);
      Double_t gbep=gamma*bt*std::sqrt(TotalMomCM*TotalMomCM+KaonMass*KaonMass)/TotalMomCM;
      Double_t a  = gamma*gamma+cottLab*cottLab;
      Double_t bp = gamma*gbep;
      Double_t c  = gbep*gbep-cottLab*cottLab;
      Double_t dd = bp*bp-a*c;

      if(dd<0.){
	std::cerr << "dd<0." << std::endl;
	dd = 0.;
      }

      Double_t costCM = (std::sqrt(dd)-bp)/a;
      if(costCM>1. || costCM<-1.){
	std::cerr << "costCM>1. || costCM<-1." << std::endl;
	costCM=-1.;
      }
      Double_t sintCM  = std::sqrt(1.-costCM*costCM);
      Double_t KaonMom = TotalMomCM*sintCM/std::sqrt(1.-costLab*costLab);

      Int_t inside = 0;
      if(true
         && TMath::Abs(vert.x()-9) < 30.
         && TMath::Abs(vert.y()) < 25.
         && TMath::Abs(vert.z()+50) < 100
         && closedist < 100.
        ){
        inside = 1;
      }

      if(nkk < MaxHits){
	event.vtx[nkk]=vert.x();
	event.vty[nkk]=vert.y();
	event.vtz[nkk]=vert.z();
	event.closeDist[nkk] = closedist;
        event.inside[nkk]  = inside;
	event.theta[nkk]     = theta;
	event.thetaCM[nkk]   = TMath::ACos(costCM)*TMath::RadToDeg();
	event.costCM[nkk]    = costCM;

	event.MissMass[nkk]       = MissMass;
	event.MissMassCorr[nkk]   = MissMassCorr;
	//	event.MissMassCorrDE[nkk] = MissMassCorrDE;
	event.MissMassCorrDE[nkk] = 0;

	event.xkp[nkk] = xkp.x();
	event.ykp[nkk] = xkp.y();
	event.ukp[nkk] = us;
	event.vkp[nkk] = vs;
	event.xkm[nkk] = xkm.x();
	event.ykm[nkk] = xkm.y();
	event.ukm[nkk] = ub;
	event.vkm[nkk] = vb;
	event.pOrg[nkk] = pk0;
	event.pCalc[nkk] = KaonMom;
	event.pCorr[nkk] = pCorr;
	//	event.pCorrDE[nkk] = pkpCorrDE.Mag();
	event.pCorrDE[nkk] = 0;
	nkk++;
      }else{
	std::cout << FUNC_NAME << " nkk("
                  << nkk << ") exceeding MaxHits: " << MaxHits << std::endl;
      }

      // Good event histograms
      if(event.chisqrK18[ikm] < 20.
         && event.chisqrKurama[ikp] < 200.){
        HF1(6001, event.pK18[ikm]);
        HF1(6002, event.pKurama[ikp]);
        HF1(6003, event.qKurama[ikp]*event.m2[ikp]);
        HF2(6004, event.qKurama[ikp]*event.m2[ikp], event.pKurama[ikp]);
        HF1(6011, closedist);
        HF1(6012, vert.x());
        HF1(6013, vert.y());
        HF1(6014, vert.z());
        HF1(6015, MissMassCorr);
        if(inside == 1){
          HF1(6101, event.pK18[ikm]);
          HF1(6102, event.pKurama[ikp]);
          HF1(6103, event.qKurama[ikp]*event.m2[ikp]);
          HF2(6104, event.qKurama[ikp]*event.m2[ikp], event.pKurama[ikp]);
          HF1(6111, closedist);
          HF1(6112, vert.x());
          HF1(6113, vert.y());
          HF1(6114, vert.z());
          HF1(6115, MissMassCorr);
          if(event.pKurama[ikp] < 1.4 && event.qKurama[ikp] > 0){
            HF1(6201, event.pK18[ikm]);
            HF1(6202, event.pKurama[ikp]);
            HF1(6203, event.qKurama[ikp]*event.m2[ikp]);
            HF2(6204, event.qKurama[ikp]*event.m2[ikp], event.pKurama[ikp]);
            HF1(6211, closedist);
            HF1(6212, vert.x());
            HF1(6213, vert.y());
            HF1(6214, vert.z());
            HF1(6215, MissMassCorr);
            if(event.m2[ikp] > 0.15
               && event.m2[ikp] < 0.40
               // && event.trigflag[20]>0
              ){
              HF1(6301, event.pK18[ikm]);
              HF1(6302, event.pKurama[ikp]);
              HF1(6303, event.qKurama[ikp]*event.m2[ikp]);
              HF2(6304, event.qKurama[ikp]*event.m2[ikp], event.pKurama[ikp]);
              HF1(6311, closedist);
              HF1(6312, vert.x());
              HF1(6313, vert.y());
              HF1(6314, vert.z());
              HF1(6315, MissMassCorr);
            }
          }
        }
      }

      ///// Matrix
      if(event.chisqrKurama[ikp] < 200
         && event.nhTof == 1
         && event.nhSch == 1
         && event.nKK == 1
         && inside == 1){
        for(Int_t isch=0; isch<nhSch; ++isch){
          for(Int_t itof=0; itof<nhTof; ++itof){
            HF2(501, event.SchSeg[isch]-1, event.TofSeg[itof]-1);
            if(event.qKurama[ikp] > 0){
              if(event.m2[ikp] < 0.12){
                HF2(502, event.SchSeg[isch]-1, event.TofSeg[itof]-1);
              }
              if(event.pKurama[ikp] < 1.4
                 && event.m2[ikp] > 0.15
                 && event.m2[ikp] < 0.4){
                HF2(503, event.SchSeg[isch]-1, event.TofSeg[itof]-1);
              }
              if(event.pKurama[ikp] > 1.8
                 && event.m2[ikp] > 0.5
                 && theta < 10.){
                HF2(504, event.SchSeg[isch]-1, event.TofSeg[itof]-1);
              }
            }
          }
        }
      }
    }
  }

  HF1(1, 10.);

  //Final Hodoscope histograms
  for(Int_t i=0; i<nhBh2; ++i){
    HF1(152, event.csBh2[i]);
    HF1(153, event.Bh2Seg[i]);
    HF1(154, event.tBh2[i]);
    HF1(155, event.deBh2[i]);
  }

  for(Int_t i=0; i<nhBh1; ++i){
    HF1(252, event.csBh1[i]);
    HF1(253, event.Bh1Seg[i]);
    HF1(254, event.tBh1[i]);
    HF1(255, event.deBh1[i]);
  }

  for(Int_t i=0; i<nhTof; ++i){
    HF1(352, event.csTof[i]);
    HF1(353, event.TofSeg[i]);
    HF1(354, event.tTof[i]);
    HF1(355, event.deTof[i]);
  }

  return true;
}

//_____________________________________________________________________________
bool
dst::DstClose()
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  const std::size_t n = TFileCont.size();
  for(std::size_t i=0; i<n; ++i){
    if(TTreeCont[i]) delete TTreeCont[i];
    if(TFileCont[i]) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________________
bool
ConfMan::InitializeHistograms()
{
  HB1(1, "Status", 60, 0., 60.);

  HB1(101, "#Clusters BH2",   7, 0., 7.);
  HB1(102, "ClusterSize BH2", 7, 0., 7.);
  HB1(103, "HitPat BH2", 10, 0., 10.);
  HB1(104, "MeanTime BH2", 200, -10., 10.);
  HB1(105, "Delta-E BH2", 200, -0.5, 4.5);

  HB1(112, "ClusterSize BH2", 7, 0., 7.);
  HB1(113, "HitPat BH2", 10, 0., 10.);
  HB1(114, "MeanTime BH2", 200, -10., 10.);
  HB1(115, "Delta-E BH2", 200, -0.5, 4.5);

  HB1(122, "ClusterSize BH2 [T0]", 7, 0., 7.);
  HB1(123, "HitPat BH2 [T0]", 10, 0., 10.);
  HB1(124, "MeanTime BH2 [T0]", 200, -10., 10.);
  HB1(125, "Delta-E BH2 [T0]", 200, -0.5, 4.5);

  HB1(152, "ClusterSize BH2 [kk]", 7, 0., 7.);
  HB1(153, "HitPat BH2 [kk]", 10, 0., 10.);
  HB1(154, "MeanTime BH2 [kk]", 200, -10., 10.);
  HB1(155, "Delta-E BH2 [kk]", 200, -0.5, 4.5);

  HB1(201, "#Clusters BH1",  11, 0., 11.);
  HB1(202, "ClusterSize BH1",11, 0., 11.);
  HB1(203, "HitPat BH1", 11, 0., 11.);
  HB1(204, "MeanTime BH1", 200, -10., 10.);
  HB1(205, "Delta-E BH1", 200, -0.5, 4.5);
  HB1(206, "Beam ToF", 200, -10., 10.);

  HB1(211, "#Clusters BH1 [km]",  11, 0., 11.);
  HB1(212, "ClusterSize BH1 [km]",11, 0., 11.);
  HB1(213, "HitPat BH1 [km]", 11, 0., 11.);
  HB1(214, "MeanTime BH1 [km]", 200, -10., 10.);
  HB1(215, "Delta-E BH1 [km]", 200, -0.5, 4.5);

  HB1(252, "ClusterSize BH1 [kk]",11, 0., 11.);
  HB1(253, "HitPat BH1 [kk]", 11, 0., 11.);
  HB1(254, "MeanTime BH1 [kk]", 200, -10., 10.);
  HB1(255, "Delta-E BH1 [kk]", 200, -0.5, 4.5);
  HB1(256, "Beam ToF [kk]", 200, -10., 10.);

  HB1(301, "#Clusters Tof",  32, 0., 32.);
  HB1(302, "ClusterSize Tof",32, 0., 32.);
  HB1(303, "HitPat Tof", 32, 0., 32.);
  HB1(304, "TimeOfFlight Tof", 500, -50., 100.);
  HB1(305, "Delta-E Tof", 200, -0.5, 4.5);

  HB1(311, "#Clusters Tof [Good]",  32, 0., 32.);
  HB1(312, "ClusterSize Tof [Good]",32, 0., 32.);
  HB1(313, "HitPat Tof [Good]", 32, 0., 32.);
  HB1(314, "TimeOfFlight Tof [Good]", 500, -50., 100.);
  HB1(315, "Delta-E Tof [Good]", 200, -0.5, 4.5);

  HB1(352, "ClusterSize Tof [kk]",32, 0., 32.);
  HB1(353, "HitPat Tof [kk]", 32, 0., 32.);
  HB1(354, "TimeOfFlight Tof [kk]", 500, -50., 100.);
  HB1(355, "Delta-E Tof [kk]", 200, -0.5, 4.5);

  HB2(501, "HitPat TOF%SCH [All]", NumOfSegSCH, 0., NumOfSegSCH, NumOfSegTOF, 0., NumOfSegTOF);
  HB2(502, "HitPat TOF%SCH [Pion]", NumOfSegSCH, 0., NumOfSegSCH, NumOfSegTOF, 0., NumOfSegTOF);
  HB2(503, "HitPat TOF%SCH [Kaon]", NumOfSegSCH, 0., NumOfSegSCH, NumOfSegTOF, 0., NumOfSegTOF);
  HB2(504, "HitPat TOF%SCH [Proton]", NumOfSegSCH, 0., NumOfSegSCH, NumOfSegTOF, 0., NumOfSegTOF);

  HB1(1001, "#Tracks SdcIn", 10, 0., 10.);
  HB1(1002, "#Hits SdcIn", 20, 0., 20.);
  HB1(1003, "Chisqr SdcIn", 200, 0., 100.);
  HB1(1004, "X0 SdcIn", 500, -200., 200.);
  HB1(1005, "Y0 SdcIn", 500, -200., 200.);
  HB1(1006, "U0 SdcIn",  700, -0.35, 0.35);
  HB1(1007, "V0 SdcIn",  400, -0.20, 0.20);
  HB2(1008, "X0%U0 SdcIn", 120, -200., 200., 100, -0.35, 0.35);
  HB2(1009, "Y0%V0 SdcIn", 100, -200., 200., 100, -0.20, 0.20);
  HB2(1010, "X0%Y0 SdcIn", 100, -200., 200., 100, -200., 200.);

  HB1(1101, "#Tracks BcOut", 10, 0., 10.);
  HB1(1102, "#Hits BcOut", 20, 0., 20.);
  HB1(1103, "Chisqr BcOut", 200, 0., 100.);
  HB1(1104, "X0 BcOut", 500, -200., 200.);
  HB1(1105, "Y0 BcOut", 500, -200., 200.);
  HB1(1106, "U0 BcOut",  700, -0.35, 0.35);
  HB1(1107, "V0 BcOut",  400, -0.20, 0.20);
  HB2(1108, "X0%U0 BcOut", 120, -200., 200., 100, -0.35, 0.35);
  HB2(1109, "Y0%V0 BcOut", 100, -200., 200., 100, -0.20, 0.20);
  HB2(1110, "X0%Y0 BcOut", 100, -200., 200., 100, -200., 200.);
  HB1(1111, "Xtgt BcOut", 500, -200., 200.);
  HB1(1112, "Ytgt BcOut", 500, -200., 200.);
  HB2(1113, "Xtgt%Ytgt BcOut", 100, -200., 200., 100, -200., 200.);

  HB1(1201, "#Tracks SdcOut", 10, 0., 10.);
  HB1(1202, "#Hits SdcOut", 20, 0., 20.);
  HB1(1203, "Chisqr SdcOut", 200, 0., 100.);
  HB1(1204, "X0 SdcOut", 600, -1200., 1200.);
  HB1(1205, "Y0 SdcOut", 600, -600., 600.);
  HB1(1206, "U0 SdcOut",  700, -0.35, 0.35);
  HB1(1207, "V0 SdcOut",  400, -0.20, 0.20);
  HB2(1208, "X0%U0 SdcOut", 120, -1200., 1200., 100, -0.35, 0.35);
  HB2(1209, "Y0%V0 SdcOut", 100,  -600.,  600., 100, -0.20, 0.20);
  HB2(1210, "X0%Y0 SdcOut", 100, -1200., 1200., 100, -600., 600.);

  HB1(2201, "#Tracks K18", 10, 0., 10.);
  HB1(2202, "#Hits K18", 30, 0., 30.);
  HB1(2203, "Chisqr K18", 500, 0., 50.);
  HB1(2204, "P K18", 1000, 0.5, 2.0);
  HB1(2251, "#Tracks K18 [Good]", 10, 0., 10.);

  HB1(3001, "#Tracks Kurama", 10, 0., 10.);
  HB1(3002, "#Hits Kurama", 30, 0., 30.);
  HB1(3003, "Chisqr Kurama", 500, 0., 500.);
  HB1(3004, "Xtgt Kurama", 500, -200., 200.);
  HB1(3005, "Ytgt Kurama", 500, -100., 100.);
  HB1(3006, "Utgt Kurama", 400, -0.35, 0.35);
  HB1(3007, "Vtgt Kurama", 200, -0.20, 0.20);
  HB2(3008, "Xtgt%U Kurama", 100, -200., 200., 100, -0.35, 0.35);
  HB2(3009, "Ytgt%V Kurama", 100, -100., 100., 100, -0.20, 0.20);
  HB2(3010, "Xtgt%Ytgt Kurama", 100, -200., 200., 100, -100., 100.);
  HB1(3011, "P Kurama", 200, 0.0, 1.0);
  HB1(3012, "PathLength Kurama", 600, 3000., 6000.);
  HB1(3013, "MassSqr", 600, -1.2, 1.2);

  HB1(3101, "#Tracks Kurama [Good]", 10, 0., 10.);
  HB1(3102, "#Hits Kurama [Good]", 30, 0., 30.);
  HB1(3103, "Chisqr Kurama [Good]", 500, 0., 500.);
  HB1(3104, "Xtgt Kurama [Good]", 500, -200., 200.);
  HB1(3105, "Ytgt Kurama [Good]", 500, -100., 100.);
  HB1(3106, "Utgt Kurama [Good]", 700, -0.35, 0.35);
  HB1(3107, "Vtgt Kurama [Good]", 400, -0.20, 0.20);
  HB2(3108, "Xtgt%U Kurama [Good]", 100, -200., 200., 100, -0.35, 0.35);
  HB2(3109, "Ytgt%V Kurama [Good]", 100, -100., 100., 100, -0.20, 0.20);
  HB2(3110, "Xtgt%Ytgt Kurama [Good]", 100, -200., 200., 100, -100., 100.);
  HB1(3111, "P Kurama [Good]", 200, 0.0, 1.0);
  HB1(3112, "PathLength Kurama [Good]", 600, 3000., 6000.);
  HB1(3113, "MassSqr", 600, -1.2, 1.2);

  HB1(3201, "#Tracks Kurama [Good2]", 10, 0., 10.);
  HB1(3202, "#Hits Kurama [Good2]", 30, 0., 30.);
  HB1(3203, "Chisqr Kurama [Good2]", 500, 0., 500.);
  HB1(3204, "Xtgt Kurama [Good2]", 500, -200., 200.);
  HB1(3205, "Ytgt Kurama [Good2]", 500, -100., 100.);
  HB1(3206, "Utgt Kurama [Good2]", 700, -0.35, 0.35);
  HB1(3207, "Vtgt Kurama [Good2]", 400, -0.20, 0.20);
  HB2(3208, "Xtgt%U Kurama [Good2]", 100, -200., 200., 100, -0.35, 0.35);
  HB2(3209, "Ytgt%V Kurama [Good2]", 100, -100., 100., 100, -0.20, 0.20);
  HB2(3210, "Xtgt%Ytgt Kurama [Good2]", 100, -200., 200., 100, -100., 100.);
  HB1(3211, "P Kurama [Good2]", 200, 0.0, 1.0);
  HB1(3212, "PathLength Kurama [Good2]", 600, 3000., 6000.);
  HB1(3213, "MassSqr", 600, -1.2, 1.2);

  HB1(4101, "#Tracks K18 [KuramaP]", 10, 0., 10.);
  HB1(4102, "#Hits K18 [KuramaP]", 30, 0., 30.);
  HB1(4103, "Chisqr K18 [KuramaP]", 500, 0., 50.);
  HB1(4104, "P K18 [KuramaP]", 500, 0.5, 2.0);
  HB1(4105, "Xtgt K18 [KuramaP]", 500, -200., 200.);
  HB1(4106, "Ytgt K18 [KuramaP]", 500, -100., 100.);
  HB1(4107, "Xout K18 [KuramaP]", 500, -200., 200.);
  HB1(4108, "Yout K18 [KuramaP]", 500, -100., 100.);
  HB1(4109, "Uout K18 [KuramaP]", 400, -0.35, 0.35);
  HB1(4110, "Vout K18 [KuramaP]", 200, -0.20, 0.20);
  HB1(4111, "Xin  K18 [KuramaP]", 500, -100., 100.);
  HB1(4112, "Yin  K18 [KuramaP]", 500, -100., 100.);
  HB1(4113, "Uin  K18 [KuramaP]", 400, -0.35, 0.35);
  HB1(4114, "Vin  K18 [KuramaP]", 200, -0.20, 0.20);

  HB1(4201, "#Tracks Kurama [Proton]", 10, 0., 10.);
  HB1(4202, "#Hits Kurama [Proton]", 30, 0., 30.);
  HB1(4203, "Chisqr Kurama [Proton]", 500, 0., 500.);
  HB1(4204, "Xtgt Kurama [Proton]", 500, -200., 200.);
  HB1(4205, "Ytgt Kurama [Proton]", 500, -100., 100.);
  HB1(4206, "Utgt Kurama [Proton]", 700, -0.35, 0.35);
  HB1(4207, "Vtgt Kurama [Proton]", 400, -0.20, 0.20);
  HB2(4208, "Xtgt%U Kurama [Proton]", 100, -200., 200., 100, -0.35, 0.35);
  HB2(4209, "Ytgt%V Kurama [Proton]", 100, -100., 100., 100, -0.20, 0.20);
  HB2(4210, "Xtgt%Ytgt Kurama [Proton]", 100, -200., 200., 100, -100., 100.);
  HB1(4211, "P Kurama [Proton]", 200, 0.0, 1.0);
  HB1(4212, "PathLength Kurama [Proton]", 600, 3000., 6000.);
  HB1(4213, "MassSqr", 600, -1.2, 1.2);

  std::vector<TString> labels = {
    "[GoodChisqr]", "[GoodVertex]", "[GoodPKurama]", "[GoodM2]"
  };
  for(Int_t i=0, n=labels.size(); i<n; ++i){
    Int_t hofs = 6000 + i*100;
    HB1(hofs+ 1, "P K18 "+labels[i], 800, 1.4, 2.2);
    HB1(hofs+ 2, "P Kurama "+labels[i], 600, 0, 3);
    HB1(hofs+ 3, "Charge*MassSquared "+labels[i], 600, -1., 2.);
    HB2(hofs+ 4, "P Kurama%Charge*MassSquared "+labels[i],
        100, -1., 2., 100, 0, 2.5);
    HB1(hofs+11, "Closest distance "+labels[i], 400, 0., 200.);
    HB1(hofs+12, "Vertex X "+labels[i], 400, -200, 200);
    HB1(hofs+13, "Vertex Y "+labels[i], 400, -200, 200);
    HB1(hofs+14, "Vertex Z "+labels[i], 400, -1000, 1000);
    HB1(hofs+15, "MissingMass "+labels[i], 200, 1.0, 2.0);
  }

  ////////////////////////////////////////////
  //Tree
  HBTree("kk","tree of KkAna");
  tree->Branch("runnum", &event.runnum, "runnum/I");
  tree->Branch("evnum",  &event.evnum,  "evnum/I");
  tree->Branch("spill",  &event.spill,  "spill/I");
  //Trigger
  tree->Branch("trigpat",   event.trigpat,  Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",  event.trigflag, Form("trigflag[%d]/I", NumOfSegTrig));

  //Hodoscope
  tree->Branch("nhBh1",   &event.nhBh1,   "nhBh1/I");
  tree->Branch("csBh1",    event.csBh1,   "csBh1[nhBh1]/I");
  tree->Branch("Bh1Seg",   event.Bh1Seg,  "Bh1Seg[nhBh1]/D");
  tree->Branch("tBh1",     event.tBh1,    "tBh1[nhBh1]/D");
  tree->Branch("dtBh1",    event.dtBh1,   "dtBh1[nhBh1]/D");
  tree->Branch("deBh1",    event.deBh1,   "deBh1[nhBh1]/D");

  tree->Branch("nhBh2",   &event.nhBh2,   "nhBh2/I");
  tree->Branch("csBh2",    event.csBh2,   "csBh2[nhBh2]/I");
  tree->Branch("Bh2Seg",   event.Bh2Seg,  "Bh2Seg[nhBh2]/D");
  tree->Branch("tBh2",     event.tBh2,    "tBh2[nhBh2]/D");
  tree->Branch("t0Bh2",    event.t0Bh2,   "t0Bh2[nhBh2]/D");
  tree->Branch("dtBh2",    event.dtBh2,   "dtBh2[nhBh2]/D");
  tree->Branch("deBh2",    event.deBh2,   "deBh2[nhBh2]/D");

  tree->Branch("Time0Seg", &event.Time0Seg, "Time0Seg/D");
  tree->Branch("deTime0", &event.deTime0, "deTime0/D");
  tree->Branch("Time0", &event.Time0, "Time0/D");
  tree->Branch("CTime0", &event.CTime0, "CTime0/D");
  tree->Branch("Btof0Seg", &event.Btof0Seg, "Btof0Seg/D");
  tree->Branch("deBtof0", &event.deBtof0, "deBtof0/D");
  tree->Branch("Btof0", &event.Btof0, "Btof0/D");
  tree->Branch("CBtof0", &event.CBtof0, "CBtof0/D");

  tree->Branch("nhTof",   &event.nhTof,   "nhTof/I");
  tree->Branch("csTof",    event.csTof,   "csTof[nhTof]/I");
  tree->Branch("TofSeg",   event.TofSeg,  "TofSeg[nhTof]/D");
  tree->Branch("tTof",     event.tTof,    "tTof[nhTof]/D");
  tree->Branch("dtTof",    event.dtTof,   "dtTof[nhTof]/D");
  tree->Branch("deTof",    event.deTof,   "deTof[nhTof]/D");

  tree->Branch("nhBvh",   &event.nhBvh,   "nhBvh/I");
  tree->Branch("BvhSeg",   event.BvhSeg,  "BvhSeg[nhBvh]/D");

  //Fiber
  tree->Branch("nhBft",  &event.nhBft,  "nhBft/I");
  tree->Branch("csBft",   event.csBft,  "csBft[nhBft]/I");
  tree->Branch("tBft",    event.tBft,   "tBft[nhBft]/D");
  tree->Branch("wBft",    event.wBft,   "wBft[nhBft]/D");
  tree->Branch("BftPos",  event.BftPos, "BftPos[nhBft]/D");
  tree->Branch("BftSeg",  event.BftSeg, "BftSeg[nhBft]/D");
  tree->Branch("nhSch",  &event.nhSch,  "nhSch/I");
  tree->Branch("csSch",   event.csSch,  "csSch[nhSch]/I");
  tree->Branch("tSch",    event.tSch,   "tSch[nhSch]/D");
  tree->Branch("wSch",    event.wSch,   "wSch[nhSch]/D");
  tree->Branch("SchPos",  event.SchPos, "SchPos[nhSch]/D");
  tree->Branch("SchSeg",  event.SchSeg, "SchSeg[nhSch]/D");

  //Beam DC
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
  tree->Branch("chisqrKurama",  event.chisqrKurama, "chisqrKurama[ntKurama]/D");
  tree->Branch("stof",          event.stof,         "stof[ntKurama]/D");
  tree->Branch("cstof",         event.cstof,        "cstof[ntKurama]/D");
  tree->Branch("path",          event.path,         "path[ntKurama]/D");
  tree->Branch("pKurama",       event.pKurama,      "pKurama[ntKurama]/D");
  tree->Branch("qKurama",       event.qKurama,      "qKurama[ntKurama]/D");
  tree->Branch("m2",            event.m2,           "m2[ntKurama]/D");
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
  tree->Branch("best_deTof",    event.best_deTof,   "best_deTof[ntKurama]/D");
  tree->Branch("best_TofSeg",   event.best_TofSeg,  "best_TofSeg[ntKurama]/D");

  //Reaction
  tree->Branch("nKm",           &event.nKm,            "nKm/I");
  tree->Branch("nKp",           &event.nKp,            "nKp/I");
  tree->Branch("nKK",           &event.nKK,            "nKK/I");
  tree->Branch("vtx",            event.vtx,            "vtx[nKK]/D");
  tree->Branch("vty",            event.vty,            "vty[nKK]/D");
  tree->Branch("vtz",            event.vtz,            "vtz[nKK]/D");
  tree->Branch("closeDist",      event.closeDist,      "closeDist[nKK]/D");
  tree->Branch("inside",         event.inside,         "inside[nKK]/I");
  tree->Branch("theta",          event.theta,          "theta[nKK]/D");
  tree->Branch("MissMass",       event.MissMass,       "MissMass[nKK]/D");
  tree->Branch("MissMassCorr",   event.MissMassCorr,   "MissMassCorr[nKK]/D");
  tree->Branch("MissMassCorrDE", event.MissMassCorrDE, "MissMassCorrDE[nKK]/D");
  tree->Branch("thetaCM", event.thetaCM,  "thetaCM[nKK]/D");
  tree->Branch("costCM",  event.costCM,   "costCM[nKK]/D");
  tree->Branch("Kflag" ,  event.Kflag,    "Kflag[nKK]/I");

  tree->Branch("xkm",        event.xkm,      "xkm[nKK]/D");
  tree->Branch("ykm",        event.ykm,      "ykm[nKK]/D");
  tree->Branch("ukm",        event.ukm,      "ukm[nKK]/D");
  tree->Branch("vkm",        event.vkm,      "vkm[nKK]/D");
  tree->Branch("xkp",        event.xkp,      "xkp[nKK]/D");
  tree->Branch("ykp",        event.ykp,      "ykp[nKK]/D");
  tree->Branch("ukp",        event.ukp,      "ukp[nKK]/D");
  tree->Branch("vkp",        event.vkp,      "vkp[nKK]/D");
  tree->Branch("pOrg",       event.pOrg,     "pOrg[nKK]/D");
  tree->Branch("pCalc",      event.pCalc,    "pCalc[nKK]/D");
  tree->Branch("pCorr",      event.pCorr,    "pCorr[nKK]/D");
  tree->Branch("pCorrDE",    event.pCorrDE,  "pCorrDE[nKK]/D");

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
  TTreeCont[kHodoscope]->SetBranchStatus("Btof0Seg", 1);
  TTreeCont[kHodoscope]->SetBranchStatus("deBtof0",  1);
  TTreeCont[kHodoscope]->SetBranchStatus("Btof0",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("CBtof0",   1);
  TTreeCont[kHodoscope]->SetBranchStatus("nhTof",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("csTof",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("TofSeg",   1);
  TTreeCont[kHodoscope]->SetBranchStatus("tTof",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtTof",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("deTof",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("nhBvh",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("BvhSeg",   1);

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
  TTreeCont[kHodoscope]->SetBranchAddress("Btof0Seg", &src.Btof0Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("deBtof0",  &src.deBtof0);
  TTreeCont[kHodoscope]->SetBranchAddress("Btof0",    &src.Btof0);
  TTreeCont[kHodoscope]->SetBranchAddress("CBtof0",   &src.CBtof0);
  TTreeCont[kHodoscope]->SetBranchAddress("nhTof",    &src.nhTof);
  TTreeCont[kHodoscope]->SetBranchAddress("csTof",     src.csTof);
  TTreeCont[kHodoscope]->SetBranchAddress("TofSeg",    src.TofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("tTof",      src.tTof);
  TTreeCont[kHodoscope]->SetBranchAddress("dtTof",     src.dtTof);
  TTreeCont[kHodoscope]->SetBranchAddress("deTof",     src.deTof);
  TTreeCont[kHodoscope]->SetBranchAddress("nhBvh",    &src.nhBvh);
  TTreeCont[kHodoscope]->SetBranchAddress("BvhSeg",    src.BvhSeg);

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
  TTreeCont[kKuramaTracking]->SetBranchStatus("cstof",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("path",        1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("pKurama",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("qKurama",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("chisqrKurama",1);
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

  TTreeCont[kKuramaTracking]->SetBranchAddress("ntSdcIn",      &src.ntSdcIn     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("nlSdcIn",      &src.nlSdcIn     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("nhSdcIn",       src.nhSdcIn     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("chisqrSdcIn",   src.chisqrSdcIn );
  TTreeCont[kKuramaTracking]->SetBranchAddress("x0SdcIn",       src.x0SdcIn     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("y0SdcIn",       src.y0SdcIn     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("u0SdcIn",       src.u0SdcIn     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("v0SdcIn",       src.v0SdcIn     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("ntSdcOut",     &src.ntSdcOut    );
  TTreeCont[kKuramaTracking]->SetBranchAddress("nhSdcOut",      src.nhSdcOut    );
  TTreeCont[kKuramaTracking]->SetBranchAddress("chisqrSdcOut",  src.chisqrSdcOut);
  TTreeCont[kKuramaTracking]->SetBranchAddress("x0SdcOut",      src.x0SdcOut    );
  TTreeCont[kKuramaTracking]->SetBranchAddress("y0SdcOut",      src.y0SdcOut    );
  TTreeCont[kKuramaTracking]->SetBranchAddress("u0SdcOut",      src.u0SdcOut    );
  TTreeCont[kKuramaTracking]->SetBranchAddress("v0SdcOut",      src.v0SdcOut    );

  TTreeCont[kKuramaTracking]->SetBranchAddress("ntKurama",    &src.ntKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("nhKurama",     src.nhKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("stof",         src.stof);
  TTreeCont[kKuramaTracking]->SetBranchAddress("cstof",        src.cstof);
  TTreeCont[kKuramaTracking]->SetBranchAddress("path",         src.path);
  TTreeCont[kKuramaTracking]->SetBranchAddress("pKurama",      src.pKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("qKurama",      src.qKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("chisqrKurama", src.chisqrKurama);
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

  TTreeCont[kK18Tracking]->SetBranchAddress("ntBcOut",     &src.ntBcOut    );
  TTreeCont[kK18Tracking]->SetBranchAddress("nlBcOut",     &src.nlBcOut    );
  TTreeCont[kK18Tracking]->SetBranchAddress("nhBcOut",      src.nhBcOut    );
  TTreeCont[kK18Tracking]->SetBranchAddress("chisqrBcOut",  src.chisqrBcOut);
  TTreeCont[kK18Tracking]->SetBranchAddress("x0BcOut",      src.x0BcOut    );
  TTreeCont[kK18Tracking]->SetBranchAddress("y0BcOut",      src.y0BcOut    );
  TTreeCont[kK18Tracking]->SetBranchAddress("u0BcOut",      src.u0BcOut    );
  TTreeCont[kK18Tracking]->SetBranchAddress("v0BcOut",      src.v0BcOut    );

  TTreeCont[kK18Tracking]->SetBranchAddress("ntK18", &src.ntK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("nhK18", &src.nhK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("chisqrK18", &src.chisqrK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("p_3rd", &src.pK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("xtgtK18", &src.xtgtK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("ytgtK18", &src.ytgtK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("utgtK18", &src.utgtK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("vtgtK18", &src.vtgtK18);

  TTreeCont[kEasiroc]->SetBranchStatus("*", 0);
  TTreeCont[kEasiroc]->SetBranchStatus("bft_ncl",    1);
  TTreeCont[kEasiroc]->SetBranchStatus("bft_clsize", 1);
  TTreeCont[kEasiroc]->SetBranchStatus("bft_ctime",  1);
  TTreeCont[kEasiroc]->SetBranchStatus("bft_ctot",   1);
  TTreeCont[kEasiroc]->SetBranchStatus("bft_clpos",  1);
  TTreeCont[kEasiroc]->SetBranchStatus("bft_clseg",  1);
  TTreeCont[kEasiroc]->SetBranchStatus("sch_ncl",    1);
  TTreeCont[kEasiroc]->SetBranchStatus("sch_clsize", 1);
  TTreeCont[kEasiroc]->SetBranchStatus("sch_ctime",  1);
  TTreeCont[kEasiroc]->SetBranchStatus("sch_ctot",   1);
  TTreeCont[kEasiroc]->SetBranchStatus("sch_clpos",  1);
  TTreeCont[kEasiroc]->SetBranchStatus("sch_clseg",  1);

  TTreeCont[kEasiroc]->SetBranchAddress("bft_ncl",    &src.nhBft);
  TTreeCont[kEasiroc]->SetBranchAddress("bft_clsize", &src.csBft);
  TTreeCont[kEasiroc]->SetBranchAddress("bft_ctime",  &src.tBft);
  TTreeCont[kEasiroc]->SetBranchAddress("bft_ctot",   &src.wBft);
  TTreeCont[kEasiroc]->SetBranchAddress("bft_clpos",  &src.BftPos);
  TTreeCont[kEasiroc]->SetBranchAddress("bft_clseg",  &src.BftSeg);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_ncl",    &src.nhSch);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_clsize", &src.csSch);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_ctime",  &src.tSch);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_ctot",   &src.wSch);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_clseg",  &src.SchSeg);

  return true;
}

//_____________________________________________________________________________
bool
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")   &&
     InitializeParameter<UserParamMan>("USER") &&
     InitializeParameter<HodoPHCMan>("HDPHC"));
}

//_____________________________________________________________________________
bool
ConfMan::FinalizeProcess()
{
  return true;
}

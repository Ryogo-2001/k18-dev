// -*- C++ -*-

#include "DstHelper.hh"

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include <TMath.h>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DatabasePDG.hh"
#include "DCGeomMan.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "HodoPHCMan.hh"

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
const Bool_t USE_M2 = true;
const Bool_t USE_XYCut = true;
TString ClassName() { return TString("DstKuramaHodoscope"); }
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kKuramaTracking, kHodoscope, kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[KuramaTracking]",
  "[Hodoscope]", "[OutFile]" };
std::vector<TString> TreeName =
{ "", "", "kurama", "hodo", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________
struct Event
{
  Int_t status;

  // SdcOut
  Int_t ntSdcOut;
  Double_t chisqrSdcOut[MaxHits];
  Double_t x0SdcOut[MaxHits];
  Double_t y0SdcOut[MaxHits];
  Double_t u0SdcOut[MaxHits];
  Double_t v0SdcOut[MaxHits];

  // KuramaTracking
  Int_t ntKurama;
  Double_t path[MaxHits];
  Double_t pKurama[MaxHits];
  Double_t qKurama[MaxHits];
  Double_t chisqrKurama[MaxHits];
  Double_t xtgtKurama[MaxHits];
  Double_t ytgtKurama[MaxHits];
  Double_t utgtKurama[MaxHits];
  Double_t vtgtKurama[MaxHits];
  Double_t thetaKurama[MaxHits];
  Double_t vpx[NumOfLayersVP*MaxHits];
  Double_t vpy[NumOfLayersVP*MaxHits];

  // Hodoscope
  Int_t trigflag[NumOfSegTrig];
  Int_t trigpat[NumOfSegTrig];

  Int_t nhBh1;
  Int_t csBh1[NumOfSegBH1*MaxDepth];
  Double_t Bh1Seg[NumOfSegBH1*MaxDepth];
  Double_t tBh1[NumOfSegBH1*MaxDepth];
  Double_t dtBh1[NumOfSegBH1*MaxDepth];
  Double_t deBh1[NumOfSegBH1*MaxDepth];

  Int_t nhBh2;
  Int_t csBh2[NumOfSegBH2*MaxDepth];
  Double_t Bh2Seg[NumOfSegBH2*MaxDepth];
  Double_t tBh2[NumOfSegBH2*MaxDepth];
  Double_t dtBh2[NumOfSegBH2*MaxDepth];
  Double_t t0Bh2[NumOfSegBH2*MaxDepth];
  Double_t deBh2[NumOfSegBH2*MaxDepth];

  Double_t Time0Seg;

  Int_t nhTof;
  Int_t csTof[NumOfSegTOF];
  Double_t TofSeg[NumOfSegTOF];
  Double_t tTof[NumOfSegTOF];
  Double_t dtTof[NumOfSegTOF];
  Double_t deTof[NumOfSegTOF];

  Double_t btof[NumOfSegBH1*MaxDepth];

  // KuramaHodoscope
  Int_t    m2Combi;
  Double_t beta[MaxHits];
  Double_t stof[MaxHits];
  Double_t cstof[MaxHits];
  Double_t m2[MaxHits];

  Double_t tofua[NumOfSegTOF];
  Double_t tofda[NumOfSegTOF];
  enum eParticle { Pion, Kaon, Proton, nParticle };
  Double_t tTofCalc[nParticle];
  Double_t utTofSeg[NumOfSegTOF][MaxDepth];
  Double_t dtTofSeg[NumOfSegTOF][MaxDepth];
  Double_t udeTofSeg[NumOfSegTOF];
  Double_t ddeTofSeg[NumOfSegTOF];

};

//_____________________________________________________________________
struct Src
{
  Int_t trigflag[NumOfSegTrig];
  Int_t trigpat[NumOfSegTrig];

  Int_t ntSdcOut;
  Double_t chisqrSdcOut[MaxHits];
  Double_t x0SdcOut[MaxHits];
  Double_t y0SdcOut[MaxHits];
  Double_t u0SdcOut[MaxHits];
  Double_t v0SdcOut[MaxHits];

  Int_t    ntKurama;
  Double_t path[MaxHits];
  Double_t pKurama[MaxHits];
  Double_t qKurama[MaxHits];
  Double_t chisqrKurama[MaxHits];
  Double_t xtgtKurama[MaxHits];
  Double_t ytgtKurama[MaxHits];
  Double_t utgtKurama[MaxHits];
  Double_t vtgtKurama[MaxHits];
  Double_t thetaKurama[MaxHits];
  Double_t vpx[NumOfLayersVP];
  Double_t vpy[NumOfLayersVP];

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
  Double_t dtBh2[NumOfSegBH2*MaxDepth];
  Double_t t0Bh2[NumOfSegBH2*MaxDepth];
  Double_t deBh2[NumOfSegBH2*MaxDepth];

  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;

  Int_t    nhTof;
  Int_t    csTof[NumOfSegTOF];
  Double_t TofSeg[NumOfSegTOF];
  Double_t tTof[NumOfSegTOF];
  Double_t dtTof[NumOfSegTOF];
  Double_t deTof[NumOfSegTOF];

  ////////// for HodoParam
  Double_t tofua[NumOfSegTOF];
  Double_t tofda[NumOfSegTOF];
  Double_t utTofSeg[NumOfSegTOF][MaxDepth];
  Double_t dtTofSeg[NumOfSegTOF][MaxDepth];
  Double_t udeTofSeg[NumOfSegTOF];
  Double_t ddeTofSeg[NumOfSegTOF];

};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
}

//_____________________________________________________________________
int
main(int argc, char **argv)
{
  std::vector<std::string> arg(argv, argv+argc);

  if(!CheckArg(arg))
    return EXIT_FAILURE;
  if(!DstOpen(arg))
    return EXIT_FAILURE;
  if(!gConf.Initialize(arg[kConfFile]))
    return EXIT_FAILURE;
  if(!gConf.InitializeUnpacker())
    return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries(TTreeCont);
  if (max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();

  Int_t ievent = skip;
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

//_____________________________________________________________________
bool
dst::InitializeEvent()
{
  event.status   = 0;
  event.ntSdcOut = 0;
  event.ntKurama = 0;
  event.nhBh1    = 0;
  event.nhBh2    = 0;
  event.nhTof    = 0;
  event.m2Combi  = 0;
  event.Time0Seg = qnan;

  for(Int_t i=0; i<NumOfSegTrig; ++i){
    event.trigflag[i] = -1;
    event.trigpat[i]  = -1;
  }

  for(Int_t i=0; i<Event::nParticle; ++i){
    event.tTofCalc[i] = qnan;
  }

  for(Int_t i=0;i<MaxHits;++i){
    event.chisqrSdcOut[i] = qnan;
    event.x0SdcOut[i]     = qnan;
    event.y0SdcOut[i]     = qnan;
    event.u0SdcOut[i]     = qnan;
    event.v0SdcOut[i]     = qnan;
    event.path[i]         = qnan;
    event.pKurama[i]      = qnan;
    event.qKurama[i]      = qnan;
    event.chisqrKurama[i] = qnan;
    event.xtgtKurama[i]  = qnan;
    event.ytgtKurama[i]  = qnan;
    event.utgtKurama[i]  = qnan;
    event.vtgtKurama[i]  = qnan;
    event.thetaKurama[i] = qnan;
  }

  for (Int_t l = 0; l < NumOfLayersVP; ++l) {
    event.vpx[l] = qnan;
    event.vpy[l] = qnan;
  }

  for(Int_t i=0;i<NumOfSegBH1*MaxDepth;++i){
    event.Bh1Seg[i] = -1;
    event.csBh1[i]  = 0;
    event.tBh1[i]   = qnan;
    event.dtBh1[i]  = qnan;
    event.deBh1[i]  = qnan;
    event.btof[i]   = qnan;
  }

  for(Int_t i=0;i<NumOfSegBH2*MaxDepth;++i){
    event.Bh2Seg[i] = -1;
    event.csBh2[i]  = 0;
    event.tBh2[i]   = qnan;
    event.dtBh2[i]  = qnan;
    event.t0Bh2[i]  = qnan;
    event.deBh2[i]  = qnan;
  }

  for(Int_t i=0;i<NumOfSegTOF;++i){
    event.TofSeg[i] = -1;
    event.csTof[i]  = 0;
    event.tTof[i]   = qnan;
    event.dtTof[i]  = qnan;
    event.deTof[i]  = qnan;
    event.tofua[i] = qnan;
    event.tofda[i] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      event.utTofSeg[i][m] = qnan;
      event.dtTofSeg[i][m] = qnan;
    }
    event.udeTofSeg[i] = qnan;
    event.ddeTofSeg[i] = qnan;
  }

  for(Int_t i=0; i<Event::nParticle; ++i){
    event.tTofCalc[i] = qnan;
  }

  for(Int_t i=0; i<MaxHits; ++i){
    event.beta[i]  = qnan;
    event.stof[i]  = qnan;
    event.cstof[i] = qnan;
    event.m2[i]    = qnan;
  }
  return true;
}

//_____________________________________________________________________
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

//_____________________________________________________________________
bool
dst::DstRead(Int_t ievent)
{
  static const auto StofOffset = gUser.GetParameter("StofOffset");

  if(ievent%10000==0){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  HF1(1, event.status++);

  event.ntSdcOut = src.ntSdcOut;
  event.ntKurama = src.ntKurama;
  event.nhBh1    = src.nhBh1;
  event.nhBh2    = src.nhBh2;
  event.nhTof    = src.nhTof;

#if 0
  std::cout<<"[event]: "<<std::setw(6)<<ievent<<" ";
  std::cout<<"[ntKurama]: "<<std::setw(2)<<src.ntKurama<<" ";
  std::cout<<"[nhBh1]: "<<std::setw(2)<<src.nhBh1<<" ";
  std::cout<<"[nhBh2]: "<<std::setw(2)<<src.nhBh2<<" ";
  std::cout<<"[nhTof]: "<<std::setw(2)<<src.nhTof<<" "<<std::endl;
#endif

  for(Int_t i=0;i<NumOfSegTrig;++i){
    Int_t tdc = src.trigflag[i];
    if(tdc<=0) continue;
    event.trigpat[i]  = i + 1;
    event.trigflag[i] = tdc;
  }

  // if(event.nhBh1<=0) return true;
  HF1(1, event.status++);

  // if(event.nhBh2<=0) return true;
  HF1(1, event.status++);

  // if(event.nhTof<=0) return true;
  HF1(1, event.status++);

  // if(event.ntKurama<=0) return true;
  // if(event.ntKurama>MaxHits)
  //   event.ntKurama = MaxHits;

  HF1(1, event.status++);

  Double_t time0 = src.CTime0;
  for(Int_t i=0; i<src.nhBh2; ++i){
    event.Bh2Seg[i] = src.Bh2Seg[i];
    event.tBh2[i]   = src.tBh2[i];
    event.t0Bh2[i]  = src.t0Bh2[i];
    event.deBh2[i]  = src.deBh2[i];
  }
  event.Time0Seg = src.Time0Seg;

  ////////// for BeamTof
  // Double_t btof = qnan;
  for(Int_t i=0; i<src.nhBh1; ++i){
    event.Bh1Seg[i] = src.Bh1Seg[i];
    event.tBh1[i]   = src.tBh1[i];
    event.deBh1[i]  = src.deBh1[i];
    event.btof[i]   = src.tBh1[i] - time0;
    // if(i==0) btof = src.tBh1[i] - time0;
  }

  for(Int_t it=0; it<src.ntSdcOut; ++it){
    event.chisqrSdcOut[it] = src.chisqrSdcOut[it];
    event.x0SdcOut[it] = src.x0SdcOut[it];
    event.y0SdcOut[it] = src.y0SdcOut[it];
    event.u0SdcOut[it] = src.u0SdcOut[it];
    event.v0SdcOut[it] = src.v0SdcOut[it];
  }

  for(Int_t it=0; it<src.nhTof; ++it){
    event.csTof[it]  = src.csTof[it];
    event.TofSeg[it] = src.TofSeg[it];
    event.tTof[it]   = src.tTof[it];
    event.dtTof[it]  = src.dtTof[it];
    event.deTof[it]  = src.deTof[it];
  }
  for(Int_t it=0; it<NumOfSegTOF; ++it){
    event.tofua[it] = src.tofua[it];
    event.tofda[it] = src.tofda[it];
    for(Int_t m=0; m<MaxDepth; ++m){
      event.utTofSeg[it][m] = src.utTofSeg[it][m];
      event.dtTofSeg[it][m] = src.dtTofSeg[it][m];
    }
    event.udeTofSeg[it] = src.udeTofSeg[it];
    event.ddeTofSeg[it] = src.ddeTofSeg[it];
    // TOF ADC-Pedestal
    Bool_t is_hit = false;
    for(Int_t i=0; i<event.nhTof; ++i){
      if(it == event.TofSeg[i]-1) is_hit = true;
    }
    if(!is_hit){
      HF1(30000+(it+1)*100+1, event.tofua[it]);
      HF1(30000+(it+1)*100+2, event.tofda[it]);
    }
  }

  Int_t m2Combi = event.nhTof*event.ntKurama;
  if(m2Combi>MaxHits || m2Combi<0){
    std::cout << FUNC_NAME << " too much m2Combi : " << m2Combi << std::endl;
    return false;
  }
  event.m2Combi = m2Combi;

  HF1(1, event.status++);
  Int_t mm=0;
  for(Int_t it=0; it<src.ntKurama; ++it){
    event.path[it]    = src.path[it];
    event.pKurama[it] = src.pKurama[it];
    event.qKurama[it] = src.qKurama[it];
    event.xtgtKurama[it] = src.xtgtKurama[it];
    event.ytgtKurama[it] = src.ytgtKurama[it];
    event.utgtKurama[it] = src.utgtKurama[it];
    event.vtgtKurama[it] = src.vtgtKurama[it];
    event.thetaKurama[it] = src.thetaKurama[it];
    event.chisqrKurama[it] = src.chisqrKurama[it];
    Double_t xtgt = event.xtgtKurama[it];
    Double_t ytgt = event.ytgtKurama[it];
    Double_t pKurama = event.pKurama[it];
    Double_t qKurama = event.qKurama[it];
    if(event.chisqrKurama[it] < 200){
      HF1(10, pKurama);
      HF1(13, event.path[it]);
    }
    if(src.ntKurama == 1){
      for (Int_t l = 0; l < NumOfLayersVP; ++l) {
	event.vpx[l] = src.vpx[l];
	event.vpy[l] = src.vpy[l];
      }
    }
    if(it == 0){
      event.tTofCalc[Event::Pion] =
        Kinematics::CalcTimeOfFlight(event.pKurama[it],
                                     event.path[it],
                                     pdg::PionMass());
      event.tTofCalc[Event::Kaon] =
        Kinematics::CalcTimeOfFlight(event.pKurama[it],
                                     event.path[it],
                                     pdg::KaonMass());
      event.tTofCalc[Event::Proton] =
        Kinematics::CalcTimeOfFlight(event.pKurama[it],
                                     event.path[it],
                                     pdg::ProtonMass());
    }

    for(Int_t itof=0; itof<src.nhTof; ++itof){
      Int_t tofseg = (Int_t)event.TofSeg[itof];

      ////////// TimeCut
      Double_t stof = event.tTof[itof] - time0 + StofOffset;
      Double_t cstof = stof;
      // gPHC.DoStofCorrection(8, 0, tofseg-1, 2, stof, btof, cstof);
      Double_t beta = event.path[it]/cstof/MathTools::C();
      event.beta[mm] = beta;
      event.stof[mm] = stof;// - event.tTofCalc[0];
      event.cstof[mm] = cstof;
      Double_t m2 = Kinematics::MassSquare(pKurama, event.path[it], cstof);
      event.m2[mm] = m2;
      // Bool_t is_pion = (// qKurama > 0
      //                   // &&
      //                   TMath::Abs(m2-pdg::PionMass()*pdg::PionMass()) < 0.1);
      // Bool_t is_kaon = (qKurama > 0
      //                   && TMath::Abs(m2-pdg::KaonMass()*pdg::KaonMass()) < 0.12);
#if 0
      std::cout << "#D DebugPrint() Event : " << ievent << std::endl
		<< std::setprecision(3) << std::fixed
		<< "   time0   : " << time0 << std::endl
		<< "   offset  : " << StofOffset << std::endl
		<< "   tTof    : " << event.tTof[itof] << std::endl
		<< "   stof    : " << stof << std::endl
		<< "   pKurama : " << pKurama << std::endl
		<< "   m2      : " << m2 << std::endl;
#endif

      if(event.chisqrKurama[it] < 200.){
        for(Int_t ip=0; ip<Event::nParticle; ++ip){
          if(USE_M2){
            if(ip == Event::Pion && TMath::Abs(m2-0.0194) > 0.1) continue;
            if(ip == Event::Proton && TMath::Abs(m2-0.88) > 0.2) continue;
          }
          HF1(33000+tofseg*100+ip+1, event.tofua[tofseg-1]);
          HF1(33000+tofseg*100+ip+1+Event::nParticle, event.tofda[tofseg-1]);
        }
        if(TMath::Abs(event.vtgtKurama[it]) < 0.01){
          HF1(20000+tofseg*100+1, event.dtTof[itof]);
          HF2(20000+1, tofseg-1, event.dtTof[itof]);
        }
      }

      Bool_t xy_ok = USE_XYCut
        ? (TMath::Abs(xtgt) < 25. && TMath::Abs(ytgt) < 20.)
        : true;
      if(event.chisqrKurama[it] < 200. && xy_ok){
        for(Int_t ip=0; ip<Event::nParticle; ++ip){
          if(USE_M2){
            if(ip == Event::Pion && m2 > 0.2) continue;
            if(ip == Event::Kaon) continue;
            if(ip == Event::Proton && (m2 < 0.5 || !xy_ok)) continue;
          }
          HF2(10000+ip+1, tofseg-1, cstof-event.tTofCalc[ip]);
          HF1(10000+tofseg*100+ip+1, cstof-event.tTofCalc[ip]);
          if(TMath::Abs(event.dtTof[itof] < 0.1)){
            for(Int_t mh=0; mh<MaxDepth; ++mh){
              HF2(40000+tofseg*100+ip+1,
                  event.udeTofSeg[tofseg-1],
                  event.tTofCalc[ip] - stof);
              HF2(40000+tofseg*100+ip+11,
                  event.ddeTofSeg[tofseg-1],
                  event.tTofCalc[ip] - stof);
              HF2(40000+tofseg*100,
                  event.udeTofSeg[tofseg-1],
                  event.tTofCalc[ip] - stof);
              HF2(40000+tofseg*100+10,
                  event.ddeTofSeg[tofseg-1],
                  event.tTofCalc[ip] - stof);
            }
          }
        }
        HF1(11, qKurama*m2);
        HF1(12, beta);
        HF1(14, stof);
        HF2(20, qKurama*m2, pKurama);
        if(!TMath::IsNaN(event.Time0Seg)){
          HF1(1100+(Int_t)event.Time0Seg, qKurama*m2);
          HF2(1200+(Int_t)event.Time0Seg, qKurama*m2, pKurama);
        }
        HF1(2100+tofseg, qKurama*m2);
        HF2(2200+tofseg, qKurama*m2, pKurama);
      }
      ++mm;
    }
  }

  HF1(1, event.status++);

  return true;
}

//_____________________________________________________________________
bool
dst::DstClose()
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  for(Int_t i=0, n=TFileCont.size(); i<n; ++i){
    if(TTreeCont[i]) delete TTreeCont[i];
    if(TFileCont[i]) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________
bool
ConfMan::InitializeHistograms()
{
  const Int_t NBinP = 220;
  const Double_t MinP =   0;
  const Double_t MaxP = 2.2;
  const Int_t NBinM2 = 340;
  const Double_t MinM2 = -1.4;
  const Double_t MaxM2 =  2.0;

  HB1(1, "Status", 21, 0., 21.);

  TString name[Event::nParticle] = { "Pion", "Kaon", "Proton" };

  HB1(10, "pKurama",    NBinP,  MinP, MaxP);
  HB1(11, "ChargexMassSquare", NBinM2, MinM2, MaxM2);
  HB1(12, "beta", 500, 0., 1.);
  HB1(13, "path", 500, 3000., 8000.);
  HB1(14, "stof", 500, 0., 100.);
  HB2(20, "pKurama % ChargexMassSquare", NBinM2, MinM2, MaxM2, NBinP, MinP, MaxP);
  for(Int_t i=0;i<NumOfSegBH2;++i){
    HB1(1100+i+1, Form("ChargexMassSquare [BH2-%d]", i+1),
        NBinM2, MinM2, MaxM2);
    HB2(1200+i+1, Form("pKurama %% ChargexMassSquare [BH2-%d]", i+1),
        NBinM2, MinM2, MaxM2, NBinP, MinP, MaxP);
  }
  for(Int_t i=0;i<NumOfSegTOF;++i){
    HB1(2100+i+1, Form("ChargexMassSquare [TOF-%d]", i+1),
        NBinM2, MinM2, MaxM2);
    HB2(2200+i+1, Form("pKurama %% ChargexMassSquare [TOF-%d]", i+1),
        NBinM2, MinM2, MaxM2, NBinP, MinP, MaxP);
  }
  // for TOF Param
  for(Int_t ip=0; ip<Event::nParticle; ++ip){
    HB2(10000+ip+1,
        Form("TofTime-%sTime %% TofSeg", name[ip].Data()),
        NumOfSegTOF, 0., (Double_t)NumOfSegTOF, 500, -25., 25.);
  }
  HB2(20001, "Tof TimeDiff U-D",
      NumOfSegTOF, 0., (Double_t)NumOfSegTOF, 500, -25., 25.);
  for(Int_t i=0; i<NumOfSegTOF; ++i){
    HB1(30000+(i+1)*100+1,
        Form("TOF ADC-Pedestal %d-U", i+1), 4000, 0., 4000.);
    HB1(30000+(i+1)*100+2,
        Form("TOF ADC-Pedestal %d-D", i+1), 4000, 0., 4000.);
    for(Int_t ip=0; ip<Event::nParticle; ++ip){
      HB1(10000+(i+1)*100+ip+1,
          Form("Tof-%d TofTime-%sTime", i+1, name[ip].Data()),
          500, -25., 25.);
      HB1(33000+(i+1)*100+ip+1,
          Form("TOF ADC-Signal %d-U [%s]", i+1, name[ip].Data()),
          4000, 0., 4000.);
      HB1(33000+(i+1)*100+ip+1+Event::nParticle,
          Form("TOF ADC-Signal %d-D [%s]", i+1, name[ip].Data()),
          4000, 0., 4000.);
      HB2(40000+(i+1)*100+ip+1,
          Form("tCalc-Time %% TOF De %d-U [%s]", i+1, name[ip].Data()),
          100, 0., 4., 100, -3., 3.);
      HB2(40000+(i+1)*100+ip+11,
          Form("tCalc-Time %% TOF De %d-D [%s]", i+1, name[ip].Data()),
          100, 0., 4., 100, -3., 3.);
    }
    HB1(20000+(i+1)*100+1, Form("Tof TimeDiff U-D %d", i+1),
        500, -25., 25.);
    HB2(40000+(i+1)*100,
        Form("tCalc-Time %% TOF De %d-U [All]", i+1),
        100, 0., 4., 100, -3., 3.);
    HB2(40000+(i+1)*100+10,
        Form("tCalc-Time %% TOF De %d-D [All]", i+1),
        100, 0., 4., 100, -3., 3.);
  }

  HBTree("khodo", "tree of DstKuramaHodoscope");
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("status",     &event.status,      "status/I");
  tree->Branch("ntSdcOut", &event.ntSdcOut,  "ntSdcOut/I");
  tree->Branch("chisqrSdcOut", event.chisqrSdcOut, "chisqrSdcOut[ntSdcOut]/D");
  tree->Branch("x0SdcOut", event.x0SdcOut, "x0SdcOut[ntSdcOut]/D");
  tree->Branch("y0SdcOut", event.y0SdcOut, "y0SdcOut[ntSdcOut]/D");
  tree->Branch("u0SdcOut", event.u0SdcOut, "u0SdcOut[ntSdcOut]/D");
  tree->Branch("v0SdcOut", event.v0SdcOut, "v0SdcOut[ntSdcOut]/D");

  tree->Branch("ntKurama", &event.ntKurama,  "ntKurama/I");
  tree->Branch("path",      event.path,      "path[ntKurama]/D");
  tree->Branch("pKurama",   event.pKurama,   "pKurama[ntKurama]/D");
  tree->Branch("qKurama",   event.qKurama,   "qKurama[ntKurama]/D");
  tree->Branch("chisqrKurama", event.chisqrKurama, "chisqrKurama[ntKurama]/D");
  tree->Branch("xtgtKurama",    event.xtgtKurama,   "xtgtKurama[ntKurama]/D");
  tree->Branch("ytgtKurama",    event.ytgtKurama,   "ytgtKurama[ntKurama]/D");
  tree->Branch("utgtKurama",    event.utgtKurama,   "utgtKurama[ntKurama]/D");
  tree->Branch("vtgtKurama",    event.vtgtKurama,   "vtgtKurama[ntKurama]/D");
  tree->Branch("thetaKurama",   event.thetaKurama,  "thetaKurama[ntKurama]/D");
  tree->Branch("vpx", event.vpx, Form("vpx[%d]/D", NumOfLayersVP));
  tree->Branch("vpy", event.vpy, Form("vpy[%d]/D", NumOfLayersVP));

  tree->Branch("tTofCalc",  event.tTofCalc,  Form("tTofCalc[%d]/D",
  						  Event::nParticle));
  tree->Branch("nhBh1", &event.nhBh1, "nhBh1/I");
  tree->Branch("csBh1",  event.csBh1, "csBh1[nhBh1]/D");
  tree->Branch("Bh1Seg", event.Bh1Seg,"Bh1Seg[nhBh1]/D");
  tree->Branch("tBh1",   event.tBh1,  "tBh1[nhBh1]/D");
  tree->Branch("dtBh1",  event.dtBh1, "dtBh1[nhBh1]/D");
  tree->Branch("deBh1",  event.deBh1, "deBh1[nhBh1]/D");
  tree->Branch("btof",   event.btof, "btof[nhBh1]/D");

  tree->Branch("nhBh2", &event.nhBh2, "nhBh2/I");
  tree->Branch("csBh2",  event.csBh2, "csBh2[nhBh2]/D");
  tree->Branch("Bh2Seg", event.Bh2Seg,"Bh2Seg[nhBh2]/D");
  tree->Branch("tBh2",   event.tBh2,  "tBh2[nhBh2]/D");
  tree->Branch("dtBh2",  event.dtBh2, "dtBh2[nhBh2]/D");
  tree->Branch("t0Bh2",  event.t0Bh2, "t0Bh2[nhBh2]/D");
  tree->Branch("deBh2",  event.deBh2, "deBh2[nhBh2]/D");
  tree->Branch("Time0Seg", &event.Time0Seg, "Time0Seg/D");

  tree->Branch("nhTof",  &event.nhTof, "nhTof/I");
  tree->Branch("csTof",   event.csTof, "csTof[nhTof]/D");
  tree->Branch("TofSeg",  event.TofSeg,"TofSeg[nhTof]/D");
  tree->Branch("tTof",    event.tTof,  "tTof[nhTof]/D");
  tree->Branch("dtTof",   event.dtTof, "dtTof[nhTof]/D");
  tree->Branch("deTof",   event.deTof, "deTof[nhTof]/D");
  tree->Branch("tofua", event.tofua, Form("tofua[%d]/D", NumOfSegTOF));
  tree->Branch("tofda", event.tofda, Form("tofda[%d]/D", NumOfSegTOF));
  tree->Branch("utTofSeg", event.utTofSeg, Form("utTofSeg[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("dtTofSeg", event.dtTofSeg, Form("dtTofSeg[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("udeTofSeg", event.udeTofSeg, Form("udeTofSeg[%d]/D", NumOfSegTOF));
  tree->Branch("ddeTofSeg", event.ddeTofSeg, Form("ddeTofSeg[%d]/D", NumOfSegTOF));

  tree->Branch("m2Combi", &event.m2Combi, "m2Combi/I");
  tree->Branch("beta",     event.beta,    "beta[m2Combi]/D");
  tree->Branch("stof",     event.stof,    "stof[m2Combi]/D");
  tree->Branch("cstof",    event.cstof,   "cstof[m2Combi]/D");
  tree->Branch("m2",       event.m2,      "m2[m2Combi]/D");

  ////////// Bring Address From Dst
  TTreeCont[kHodoscope]->SetBranchStatus("*", 0);
  TTreeCont[kHodoscope]->SetBranchStatus("trigflag",  1);
  TTreeCont[kHodoscope]->SetBranchStatus("trigpat",   1);
  TTreeCont[kHodoscope]->SetBranchStatus("nhBh1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("csBh1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("Bh1Seg",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("tBh1",      1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtBh1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("deBh1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("nhBh2",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("csBh2",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("Bh2Seg",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("tBh2",      1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtBh2",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("t0Bh2",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("deBh2",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("nhTof",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("csTof",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("TofSeg",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("tTof",      1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtTof",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("deTof",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("tofua",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("tofda",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("utTofSeg",  1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtTofSeg",  1);
  TTreeCont[kHodoscope]->SetBranchStatus("udeTofSeg", 1);
  TTreeCont[kHodoscope]->SetBranchStatus("ddeTofSeg", 1);
  TTreeCont[kHodoscope]->SetBranchStatus("Time0",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("CTime0",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("Time0Seg",  1);
  TTreeCont[kHodoscope]->SetBranchStatus("deTime0",   1);

  TTreeCont[kHodoscope]->SetBranchAddress("trigflag", src.trigflag);
  TTreeCont[kHodoscope]->SetBranchAddress("trigpat",  src.trigpat);
  TTreeCont[kHodoscope]->SetBranchAddress("nhBh1",  &src.nhBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("csBh1",  src.csBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("Bh1Seg", src.Bh1Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("tBh1",   src.tBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("dtBh1",  src.dtBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("deBh1",  src.deBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("nhBh2", &src.nhBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("csBh2",  src.csBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("Bh2Seg",src.Bh2Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("tBh2",  src.tBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("dtBh2", src.dtBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("t0Bh2", src.t0Bh2);
  TTreeCont[kHodoscope]->SetBranchAddress("deBh2", src.deBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("nhTof", &src.nhTof);
  TTreeCont[kHodoscope]->SetBranchAddress("csTof", src.csTof);
  TTreeCont[kHodoscope]->SetBranchAddress("TofSeg",src.TofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("tTof",  src.tTof);
  TTreeCont[kHodoscope]->SetBranchAddress("dtTof", src.dtTof);
  TTreeCont[kHodoscope]->SetBranchAddress("deTof", src.deTof);
  TTreeCont[kHodoscope]->SetBranchAddress("tofua", src.tofua);
  TTreeCont[kHodoscope]->SetBranchAddress("tofda", src.tofda);
  TTreeCont[kHodoscope]->SetBranchAddress("utTofSeg", src.utTofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("dtTofSeg", src.dtTofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("udeTofSeg", src.udeTofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("ddeTofSeg", src.ddeTofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("Time0", &src.Time0);
  TTreeCont[kHodoscope]->SetBranchAddress("CTime0", &src.CTime0);
  TTreeCont[kHodoscope]->SetBranchAddress("Time0Seg", &src.Time0Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("deTime0", &src.deTime0);

  TTreeCont[kKuramaTracking]->SetBranchStatus("*",      0);
  TTreeCont[kKuramaTracking]->SetBranchStatus("ntSdcOut", 1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("chisqrSdcOut", 1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("x0SdcOut", 1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("y0SdcOut", 1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("u0SdcOut", 1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("v0SdcOut", 1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("ntKurama",     1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("path",         1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("pKurama",      1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("qKurama",      1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("chisqrKurama", 1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("thetaKurama",  1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("xtgtKurama",   1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("ytgtKurama",   1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("utgtKurama",   1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("vtgtKurama",   1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("vpx",          1);
  TTreeCont[kKuramaTracking]->SetBranchStatus("vpy",          1);

  TTreeCont[kKuramaTracking]->SetBranchAddress("ntSdcOut", &src.ntSdcOut);
  TTreeCont[kKuramaTracking]->SetBranchAddress("chisqrSdcOut", src.chisqrSdcOut);
  TTreeCont[kKuramaTracking]->SetBranchAddress("x0SdcOut", src.x0SdcOut);
  TTreeCont[kKuramaTracking]->SetBranchAddress("y0SdcOut", src.y0SdcOut);
  TTreeCont[kKuramaTracking]->SetBranchAddress("u0SdcOut", src.u0SdcOut);
  TTreeCont[kKuramaTracking]->SetBranchAddress("v0SdcOut", src.v0SdcOut);
  TTreeCont[kKuramaTracking]->SetBranchAddress("ntKurama", &src.ntKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("path",     src.path);
  TTreeCont[kKuramaTracking]->SetBranchAddress("pKurama",  src.pKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("qKurama",  src.qKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("chisqrKurama", src.chisqrKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("thetaKurama", src.thetaKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("xtgtKurama", src.xtgtKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("ytgtKurama", src.ytgtKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("utgtKurama", src.utgtKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("vtgtKurama", src.vtgtKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("vpx",        src.vpx);
  TTreeCont[kKuramaTracking]->SetBranchAddress("vpy",        src.vpy);

  return true;
}

//_____________________________________________________________________
bool
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")   &&
     InitializeParameter<UserParamMan>("USER") &&
     InitializeParameter<HodoPHCMan>("HDPHC"));
}

//_____________________________________________________________________
bool
ConfMan::FinalizeProcess()
{
  return true;
}

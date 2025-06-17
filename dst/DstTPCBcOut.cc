// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include <TSystem.h>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DatabasePDG.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DstHelper.hh"
#include "HodoPHCMan.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrack.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#define TrackCluster 1
#define GainCorrection 1
#define TrackSearch 1
#define HoughYcut 0

namespace
{
using namespace root;
using namespace dst;
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
//position cut for gain histogram
const double min_ycut = -50.;//mm
const double max_ycut = 50.;//mm

const Double_t& zK18HS = gGeom.LocalZ("K18HS");
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kTpcHit,  kBcOut, kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[TPCHit]", "[BcOut]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc","bcout" ,"" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________________
struct Event
{
  Int_t status;
  Int_t runnum;
  Int_t evnum;
  std::vector<Int_t> trigpat;
  std::vector<Int_t> trigflag;
  Int_t nhTpc;
  Int_t nclTpc;
  std::vector<Double_t> raw_hitpos_x;
  std::vector<Double_t> raw_hitpos_y;
  std::vector<Double_t> raw_hitpos_z;
  std::vector<Double_t> raw_de;
  std::vector<Int_t> raw_padid;
  std::vector<Int_t> raw_layer;
  std::vector<Int_t> raw_row;
  std::vector<Double_t> cluster_x;
  std::vector<Double_t> cluster_y;
  std::vector<Double_t> cluster_z;
  std::vector<Double_t> cluster_de;
  std::vector<Int_t> cluster_size;
  std::vector<Int_t> cluster_layer;
  std::vector<Double_t> cluster_mrow;
  std::vector<Double_t> cluster_de_center;
  std::vector<Double_t> cluster_x_center;
  std::vector<Double_t> cluster_y_center;
  std::vector<Double_t> cluster_z_center;
  std::vector<Int_t> cluster_row_center;

  //BcOut info
  Int_t ntBcOut;
  std::vector<Double_t> chisqrBcOut;
  std::vector<Double_t> x0BcOut;
  std::vector<Double_t> y0BcOut;
  std::vector<Double_t> u0BcOut;
  std::vector<Double_t> v0BcOut;

  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Double_t> chisqrTpc;
  std::vector<Double_t> x0Tpc;
  std::vector<Double_t> y0Tpc;
  std::vector<Double_t> u0Tpc;
  std::vector<Double_t> v0Tpc;
  std::vector<Double_t> theta;
  std::vector<std::vector<Double_t>> hitlayer;
  std::vector<std::vector<Double_t>> hitpos_x;
  std::vector<std::vector<Double_t>> hitpos_y;
  std::vector<std::vector<Double_t>> hitpos_z;
  std::vector<std::vector<Double_t>> calpos_x;
  std::vector<std::vector<Double_t>> calpos_y;
  std::vector<std::vector<Double_t>> calpos_z;
  std::vector<std::vector<Double_t>> residual;
  std::vector<std::vector<Double_t>> residual_x;
  std::vector<std::vector<Double_t>> residual_y;
  std::vector<std::vector<Double_t>> residual_z;
  std::vector<std::vector<Double_t>> residual_wbcout_x;
  std::vector<std::vector<Double_t>> residual_wbcout_y;
  std::vector<std::vector<Double_t>> residual_trackwbcout_x;
  std::vector<std::vector<Double_t>> residual_trackwbcout_y;
  std::vector<std::vector<Double_t>> track_cluster_de;
  std::vector<std::vector<Double_t>> track_cluster_size;
  std::vector<std::vector<Double_t>> track_cluster_mrow;
  std::vector<std::vector<Double_t>> track_cluster_de_center;
  std::vector<std::vector<Double_t>> track_cluster_x_center;
  std::vector<std::vector<Double_t>> track_cluster_y_center;
  std::vector<std::vector<Double_t>> track_cluster_z_center;
  std::vector<std::vector<Double_t>> track_cluster_row_center;

  std::vector<Double_t> xCorVec; // correction vector x
  std::vector<Double_t> yCorVec; // correction vector y
  std::vector<Double_t> zCorVec; // correction vector z
  std::vector<Double_t> xCorPos; // position x
  std::vector<Double_t> yCorPos; // position y
  std::vector<Double_t> zCorPos; // position z

  std::vector<Double_t> clkTpc;

  void clear()
    {
      trigpat.clear();
      trigflag.clear();
      clkTpc.clear();

      runnum = 0;
      evnum = 0;
      status = 0;

      nhTpc = 0;
      raw_hitpos_x.clear();
      raw_hitpos_y.clear();
      raw_hitpos_z.clear();
      raw_de.clear();
      raw_padid.clear();
      raw_layer.clear();
      raw_row.clear();

      nclTpc = 0;
      cluster_x.clear();
      cluster_y.clear();
      cluster_z.clear();
      cluster_de.clear();
      cluster_size.clear();
      cluster_layer.clear();
      cluster_mrow.clear();
      cluster_de_center.clear();
      cluster_x_center.clear();
      cluster_y_center.clear();
      cluster_z_center.clear();
      cluster_row_center.clear();

      ntBcOut = 0;
      chisqrBcOut.clear();
      x0BcOut.clear();
      y0BcOut.clear();
      u0BcOut.clear();
      v0BcOut.clear();

      ntTpc = 0;
      nhtrack.clear();
      chisqrTpc.clear();
      x0Tpc.clear();
      y0Tpc.clear();
      u0Tpc.clear();
      v0Tpc.clear();
      theta.clear();
      hitlayer.clear();
      hitpos_x.clear();
      hitpos_y.clear();
      hitpos_z.clear();
      calpos_x.clear();
      calpos_y.clear();
      calpos_z.clear();
      residual.clear();
      residual_x.clear();
      residual_y.clear();
      residual_wbcout_x.clear();
      residual_wbcout_y.clear();
      residual_trackwbcout_x.clear();
      residual_trackwbcout_y.clear();
      residual_z.clear();

      track_cluster_de.clear();
      track_cluster_size.clear();
      track_cluster_mrow.clear();
      track_cluster_de_center.clear();
      track_cluster_x_center.clear();
      track_cluster_y_center.clear();
      track_cluster_z_center.clear();
      track_cluster_row_center.clear();

      xCorVec.clear();
      yCorVec.clear();
      zCorVec.clear();
      xCorPos.clear();
      yCorPos.clear();
      zCorPos.clear();
    }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* runnum;
  TTreeReaderValue<Int_t>* evnum;
  TTreeReaderValue<std::vector<Int_t>>* trigpat;
  TTreeReaderValue<std::vector<Int_t>>* trigflag;
  TTreeReaderValue<Int_t>* npadTpc;   // number of pads
  TTreeReaderValue<Int_t>* nhTpc;     // number of hits
  // vector (size=nhTpc)
  TTreeReaderValue<std::vector<Int_t>>* layerTpc;     // layer id
  TTreeReaderValue<std::vector<Int_t>>* rowTpc;       // row id
  TTreeReaderValue<std::vector<Int_t>>* padTpc;       // pad id
  TTreeReaderValue<std::vector<Double_t>>* pedTpc;    // pedestal
  TTreeReaderValue<std::vector<Double_t>>* rmsTpc;    // rms
  TTreeReaderValue<std::vector<Double_t>>* deTpc;     // dE
  TTreeReaderValue<std::vector<Double_t>>* tTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* ctTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* chisqrTpc; // chi^2 of signal fitting
  TTreeReaderValue<std::vector<Double_t>>* clkTpc;      // time

  //BcOut input
  Int_t ntBcOut;
  Double_t chisqrBcOut[MaxHits];
  Double_t x0BcOut[MaxHits];
  Double_t y0BcOut[MaxHits];
  Double_t u0BcOut[MaxHits];
  Double_t v0BcOut[MaxHits];
};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
enum eDetHid {
  TPCCalibHid    = 10000,
  TPCGainHid     = 100000,
  TPCClHid       = 200000,
  XCorrectionMapHid = 1000000,
  YCorrectionMapHid = 2000000,
  TPCPadYHid     = 3000000,
  TPCResYHid     = 4000000,
  TPCResYClkHid  = 5000000,
  TPCResXHid     = 6000000,
  TPCResXClkHid  = 7000000,
  TPCResYCoBoHid = 8000000,
  TPCResYPosHid  = 9000000,
};

const Int_t MinPosMapXZ = -300;
const Int_t MaxPosMapXZ = 300;
const Int_t MinPosMapY = -200;
const Int_t MaxPosMapY = 200;
const Int_t Meshsize = 20;
const Int_t NumOfDivXZ = ((MaxPosMapXZ - MinPosMapXZ)/Meshsize) + 1;
const Int_t NumOfDivY = ((MaxPosMapY - MinPosMapY)/Meshsize) + 1;

//_____________________________________________________________________________
Int_t
XyzToHid(Double_t x, Double_t y, Double_t z)
{
  Int_t ix = TMath::Nint((x - MinPosMapXZ)/Meshsize);
  Int_t iy = TMath::Nint((y - MinPosMapY)/Meshsize);
  Int_t iz = TMath::Nint((z - MinPosMapXZ)/Meshsize);
  return ix*NumOfDivY*NumOfDivXZ + iy*NumOfDivXZ + iz;
}

//_____________________________________________________________________________
Int_t
XyzToHid(const TVector3& pos)
{
  return XyzToHid(pos.X(), pos.Y(), pos.Z());
}

//_____________________________________________________________________________
TVector3
HidToXyz(Int_t hid)
{
  Int_t ix = TMath::FloorNint(hid/NumOfDivY/NumOfDivXZ);
  Int_t iy = TMath::FloorNint((hid - ix*NumOfDivY*NumOfDivXZ)/NumOfDivXZ);
  Int_t iz = TMath::FloorNint((hid - ix*NumOfDivY*NumOfDivXZ - iy*NumOfDivXZ));
  return TVector3(ix*Meshsize + MinPosMapXZ,
                  iy*Meshsize + MinPosMapY,
                  iz*Meshsize + MinPosMapXZ);
}

//_____________________________________________________________________________
Double_t
TranseverseDistance(Double_t x_center, Double_t z_center, Double_t x, Double_t z)
{
  Double_t dummy = TMath::Sqrt((x-x_center)*(x-x_center) + (z-z_center)*(z-z_center));
  Double_t dist;
  if(x_center-x<0) dist=-1.*dummy;
  else dist=dummy;
  return dist;
}
}

//_____________________________________________________________________________
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

//_____________________________________________________________________________
Bool_t
dst::InitializeEvent()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstOpen(std::vector<std::string> arg)
{
  int open_file = 0;
  int open_tree = 0;
  for(Int_t i=0; i<nArgc; ++i){
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
Bool_t
dst::DstRead(int ievent)
{


  if(ievent%1000==0){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;

  event.clkTpc = **src.clkTpc;

  event.ntBcOut = src.ntBcOut;

  for(int it=0; it<src.ntBcOut; ++it){
    event.chisqrBcOut.push_back(src.chisqrBcOut[it]);
    event.x0BcOut.push_back(src.x0BcOut[it]);
    event.y0BcOut.push_back(src.y0BcOut[it]);
    event.u0BcOut.push_back(src.u0BcOut[it]);
    event.v0BcOut.push_back(src.v0BcOut[it]);
  }

  HF1(1, event.status++);

  if(**src.nhTpc == 0)
    return true;

  HF1(1, event.status++);

  if(event.clkTpc.size() != 1){
    std::cerr << "something is wrong" << std::endl;
    return true;
  }
  Double_t clock = event.clkTpc.at(0);

  DCAnalyzer DCAna;
  DCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, clock);

#if HoughYcut
  DCAna.HoughYCut(min_ycut, max_ycut);
#endif

  HF1(1, event.status++);

  Int_t nhTpc = 0;
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    auto hc = DCAna.GetTPCHC(layer);
    for(const auto& hit : hc){
      if(!hit || !hit->IsGood())
        continue;
      const auto& pos = hit->GetPosition();
      Double_t x = pos.X();
      Double_t y = pos.Y();
      Double_t z = pos.Z();
      Double_t de = hit->GetCDe();
      Double_t pad = hit->GetPad();
      Int_t row = hit->GetRow();
      Int_t cobo = tpc::GetCoBoId(layer, row);
      event.raw_hitpos_x.push_back(x);
      event.raw_hitpos_y.push_back(y);
      event.raw_hitpos_z.push_back(z);
      event.raw_de.push_back(de);
      event.raw_padid.push_back(pad);
      event.raw_layer.push_back(layer);
      event.raw_row.push_back(row);

      Double_t zTPC = zK18HS + z;
      for(int it=0; it<src.ntBcOut; ++it){
	Double_t x0BcOut = src.x0BcOut[it];
     	Double_t u0BcOut = src.u0BcOut[it];
     	Double_t y0BcOut = src.y0BcOut[it];
     	Double_t v0BcOut = src.v0BcOut[it];
     	Double_t xbc = x0BcOut + zTPC*u0BcOut;
     	Double_t ybc = y0BcOut + zTPC*v0BcOut;
        Double_t resx = x - xbc;
        Double_t resy = y - ybc;
	HF2(TPCCalibHid, layer, resy);
        if(src.ntBcOut == 1 && src.chisqrBcOut[it] < 10.){
	  //	  HF2(TPCCalibHid, layer, resy);
          HF1(TPCResYHid+layer*1000, resy);
          HF1(TPCResYHid+layer*1000+row+1, resy);
          HF2(TPCResYClkHid+layer*1000, clock, resy);
          HF2(TPCResYCoBoHid+cobo*1000, clock, resy);
          HFProf(TPCResYCoBoHid+cobo*1000+1, clock, resy);
          HF2(TPCResYPosHid+layer*1000, ybc, resy);
	  HFProf(TPCResYPosHid+cobo*1000+1, ybc, resy);
          HF1(TPCResXHid+layer*1000, resx);
          HF1(TPCResXHid+layer*1000+row+1, resx);
          HF2(TPCResXClkHid+layer*1000, clock, resx);
        }
      }
      ++nhTpc;
    }
  }
  event.nhTpc = nhTpc;

  HF1(1, event.status++);
  Int_t nclTpc = 0;
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    auto hc = DCAna.GetTPCClCont(layer);
    for(const auto& cl : hc){
      if(!cl || !cl->IsGood())
        continue;
      Double_t x = cl->GetX();
      Double_t y = cl->GetY();
      Double_t z = cl->GetZ();
      Double_t clde = cl->GetDe();
      Int_t cs = cl->GetClusterSize();
      Double_t mrow = cl->MeanRow(); // same
      TPCHit* centerHit = cl->GetCenterHit();
      const TVector3& centerPos = centerHit->GetPosition();
      Double_t centerDe = centerHit->GetCDe();
      Int_t centerRow = centerHit->GetRow();

      event.cluster_x.push_back(x);
      event.cluster_y.push_back(y);
      event.cluster_z.push_back(z);
      event.cluster_de.push_back(clde);
      event.cluster_size.push_back(cs);
      event.cluster_layer.push_back(layer);
      event.cluster_mrow.push_back(mrow);
      event.cluster_de_center.push_back(centerDe);
      event.cluster_x_center.push_back(centerPos.X());
      event.cluster_y_center.push_back(centerPos.Y());
      event.cluster_z_center.push_back(centerPos.Z());
      event.cluster_row_center.push_back(centerRow);
      ++nclTpc;
    }
  }
  event.nclTpc = nclTpc;

  HF1(1, event.status++);

#if TrackSearch
  DCAna.TrackSearchTPC();
#endif

  Int_t ntTpc = DCAna.GetNTracksTPC();
  event.ntTpc = ntTpc;
  HF1(10, ntTpc);
  if(event.ntTpc == 0)
    return true;

  HF1(1, event.status++);

  event.nhtrack.resize(ntTpc);
  event.chisqrTpc.resize(ntTpc);
  event.x0Tpc.resize(ntTpc);
  event.y0Tpc.resize(ntTpc);
  event.u0Tpc.resize(ntTpc);
  event.v0Tpc.resize(ntTpc);
  event.theta.resize(ntTpc);
  event.hitlayer.resize(ntTpc);
  event.hitpos_x.resize(ntTpc);
  event.hitpos_y.resize(ntTpc);
  event.hitpos_z.resize(ntTpc);
  event.calpos_x.resize(ntTpc);
  event.calpos_y.resize(ntTpc);
  event.calpos_z.resize(ntTpc);
  event.residual.resize(ntTpc);
  event.residual_x.resize(ntTpc);
  event.residual_y.resize(ntTpc);
  event.residual_wbcout_x.resize(ntTpc);
  event.residual_wbcout_y.resize(ntTpc);
  event.residual_trackwbcout_x.resize(ntTpc);
  event.residual_trackwbcout_y.resize(ntTpc);
  event.residual_z.resize(ntTpc);
#if TrackCluster
  event.track_cluster_de.resize(ntTpc);
  event.track_cluster_size.resize(ntTpc);
  event.track_cluster_mrow.resize(ntTpc);
  event.track_cluster_de_center.resize(ntTpc);
  event.track_cluster_x_center.resize(ntTpc);
  event.track_cluster_y_center.resize(ntTpc);
  event.track_cluster_z_center.resize(ntTpc);
  event.track_cluster_row_center.resize(ntTpc);
#endif
  for(Int_t it=0; it<ntTpc; ++it){
    TPCLocalTrack *tp = DCAna.GetTrackTPC(it);
    if(!tp) continue;
    Int_t nh = tp->GetNHit();
    Double_t chisqr = tp->GetChiSquare();
    Double_t x0Tpc=tp->GetX0(), y0Tpc=tp->GetY0();
    Double_t u0Tpc=tp->GetU0(), v0Tpc=tp->GetV0();
    Double_t theta = tp->GetTheta();
    HF1(11, nh);
    HF1(12, chisqr);
    HF1(14, x0Tpc);
    HF1(15, y0Tpc);
    HF1(16, u0Tpc);
    HF1(17, v0Tpc);
    HF2(18, x0Tpc, u0Tpc);
    HF2(19, y0Tpc, v0Tpc);
    HF2(20, x0Tpc, y0Tpc);

    event.nhtrack[it] = nh;
    event.chisqrTpc[it] = chisqr;
    event.x0Tpc[it] = x0Tpc;
    event.y0Tpc[it] = y0Tpc;
    event.u0Tpc[it] = u0Tpc;
    event.v0Tpc[it] = v0Tpc;
    event.theta[it] = theta;
    event.hitlayer[it].resize(nh);
    event.hitpos_x[it].resize(nh);
    event.hitpos_y[it].resize(nh);
    event.hitpos_z[it].resize(nh);
    event.calpos_x[it].resize(nh);
    event.calpos_y[it].resize(nh);
    event.calpos_z[it].resize(nh);
    event.residual[it].resize(nh);
    event.residual_x[it].resize(nh);
    event.residual_y[it].resize(nh);
    event.residual_wbcout_x[it].resize(nh);
    event.residual_wbcout_y[it].resize(nh);
    event.residual_trackwbcout_x[it].resize(nh);
    event.residual_trackwbcout_y[it].resize(nh);
    event.residual_z[it].resize(nh);
#if TrackCluster
    event.track_cluster_de[it].resize(nh);
    event.track_cluster_size[it].resize(nh);
    event.track_cluster_mrow[it].resize(nh);
    event.track_cluster_de_center[it].resize(nh);
    event.track_cluster_x_center[it].resize(nh);
    event.track_cluster_y_center[it].resize(nh);
    event.track_cluster_z_center[it].resize(nh);
    event.track_cluster_row_center[it].resize(nh);
#endif
    for(int ih=0; ih<nh; ++ih){
      TPCLTrackHit *hit = tp->GetHit(ih);
      if(!hit) continue;
      Int_t layer = hit->GetLayer();
      const TVector3& hitpos = hit->GetLocalHitPos();
      const TVector3& calpos = hit->GetLocalCalPos();
      const TVector3& res_vect = hit->GetResidualVect();
      Double_t residual = hit->GetResidual();
      HF1(13, layer);
      event.hitlayer[it][ih] = layer;
      event.hitpos_x[it][ih] = hitpos.x();
      event.hitpos_y[it][ih] = hitpos.y();
      event.hitpos_z[it][ih] = hitpos.z();
      event.calpos_x[it][ih] = calpos.x();
      event.calpos_y[it][ih] = calpos.y();
      event.calpos_z[it][ih] = calpos.z();
      event.residual[it][ih] = residual;
      event.residual_x[it][ih] = res_vect.x();
      event.residual_y[it][ih] = res_vect.y();
      event.residual_z[it][ih] = res_vect.z();

#if TrackCluster
      TPCHit *clhit = hit -> GetHit();
      TPCCluster *cl = clhit -> GetParentCluster();
      Int_t cs = cl->GetClusterSize();
      Double_t clde = cl->GetDe();
      Double_t mrow = cl->MeanRow(); // same
      TPCHit* centerHit = cl->GetCenterHit();
      const TVector3& centerPos = centerHit->GetPosition();
      Double_t centerDe = centerHit->GetCDe();
      Int_t centerRow = centerHit->GetRow();

      HF1(TPCClHid, cs);
      HF1(TPCClHid+(layer+1)*1000, cs);
      HF1(TPCClHid+1, clde);
      HF1(TPCClHid+(layer+1)*1000+1, clde);
      if(cs<11){
	HF1(TPCClHid+9+cs, clde);
	HF1(TPCClHid+(layer+1)*1000+9+cs, clde);
      }
      const TPCHitContainer& hc = cl -> GetHitContainer();
      for(const auto& hits : hc){
	if(!hits || !hits->IsGood()) continue;
	const TVector3& pos = hits->GetPosition();
	Double_t de = hits->GetCDe();
	Double_t transDist = TranseverseDistance(hitpos.x(), hitpos.z(), pos.x(), pos.z());
	Double_t ratio = de/clde;
	HF2(TPCClHid+2, transDist, ratio);
	HF2(TPCClHid+(layer+1)*1000+2, transDist, ratio);
	if(cs<11) HF2(TPCClHid+(layer+1)*1000+19+cs, transDist, ratio);
      }
      event.track_cluster_de[it][ih] = clde;
      event.track_cluster_size[it][ih] = cs;
      event.track_cluster_mrow[it][ih] = mrow;
      event.track_cluster_de_center[it][ih] = centerDe;
      event.track_cluster_x_center[it][ih] = centerPos.X();
      event.track_cluster_y_center[it][ih] = centerPos.Y();
      event.track_cluster_z_center[it][ih] = centerPos.Z();
      event.track_cluster_row_center[it][ih] = centerRow;
#endif
#if GainCorrection
      HF1(TPCGainHid+(layer+1)*1000+centerRow, centerDe);
#endif
    }
  }
#if 0
  std::cout<<"[event]: "<<std::setw(6)<<ievent<<" "<<std::endl;
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstClose()
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  const Int_t n = TFileCont.size();
  for(Int_t i=0; i<n; ++i){
    if(TTreeReaderCont[i]) delete TTreeReaderCont[i];
    if(TTreeCont[i]) delete TTreeCont[i];
    if(TFileCont[i]) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  const Int_t    NbinDe = 1000;
  const Double_t MinDe  =    0.;
  const Double_t MaxDe  = 2000.;
  const Int_t    NbinHitDe = 500;
  const Double_t MinHitDe  =   0.;
  const Double_t MaxHitDe  = 1000.;
  const Int_t    NbinRes = 400;
  const Double_t MinRes  = -20.;
  const Double_t MaxRes  =  20.;
  const Int_t    NbinPos = 400;
  const Double_t MinPos  = -100.;
  const Double_t MaxPos  = 100.;
  // const Int_t    NbinClk = 20000;
  // const Double_t MinClk  = 100.;
  // const Double_t MaxClk  = 100.;
#if TrackCluster
  const Int_t NbinClSize = 25;
  const Double_t MinClSize = 0;
  const Double_t MaxClSize = 25;
  const Int_t NbinDist = 300;
  const Double_t MinDist = -15.;
  const Double_t MaxDist = 15.;
  const Int_t NbinRatio = 100;
  const Double_t MinRatio = 0.;
  const Double_t MaxRatio = 1.;
#endif

  HB1(1, "Status", 21, 0., 21.);
  HB1(10, "#Tracks TPC", 40, 0., 40.);
  HB1(11, "#Hits of Track TPC", 50, 0., 50.);
  HB1(12, "Chisqr TPC", 500, 0., 500.);
  HB1(13, "LayerId TPC", 35, 0., 35.);
  HB1(14, "X0 TPC", 400, -100., 100.);
  HB1(15, "Y0 TPC", 400, -100., 100.);
  HB1(16, "U0 TPC", 200, -0.20, 0.20);
  HB1(17, "V0 TPC", 200, -0.20, 0.20);
  HB2(18, "U0%X0 TPC", 100, -100., 100., 100, -0.20, 0.20);
  HB2(19, "V0%Y0 TPC", 100, -100., 100., 100, -0.20, 0.20);
  HB2(20, "X0%Y0 TPC", 100, -100., 100., 100, -100, 100);
  HB2(TPCCalibHid,"Layer%ResY", 32, 0, 32, 400, -5, 5);
#if TrackCluster
  HB1(TPCClHid, "Cluster size;Cluster size;Counts", NbinClSize, MinClSize, MaxClSize);
  HB1(TPCClHid+1, "Cluster dE;Cluster dE;Counts", NbinDe, MinDe, MaxDe);
  HB2(TPCClHid+2, "Transverse diffusion;X_{cluster_center}-X_{pad};A/A_{sum}", NbinDist, MinDist, MaxDist, NbinRatio, MinRatio, MaxRatio);
  for(Int_t i=0; i<10; ++i) HB1(TPCClHid+10+i, Form("Cluster dE, Cluster size=%d;Cluster dE;Counts",i+1), NbinDe, MinDe, MaxDe);
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    HB1(TPCClHid+(layer+1)*1000, Form("Cluster size layer%d;Cluster size;Counts",layer), NbinClSize, MinClSize, MaxClSize);
    HB1(TPCClHid+(layer+1)*1000+1, Form("Cluster dE layer%d;Cluster dE;Counts",layer), NbinDe, MinDe, MaxDe);
    HB2(TPCClHid+(layer+1)*1000+2, Form("Transverse diffusion Layer%d;X_{cluster_center}-X_{pad};A/A_{sum}",layer), NbinDist, MinDist, MaxDist, NbinRatio, MinRatio, MaxRatio);
    for(Int_t i=0; i<10; ++i){
      HB1(TPCClHid+(layer+1)*1000+10+i, Form("Cluster dE layer%d, Cluster size=%d;Cluster dE;Counts",layer,i+1), NbinDe, MinDe, MaxDe);
      HB2(TPCClHid+(layer+1)*1000+20+i, Form("Transverse diffusion, clsize=%d Layer%d;X_{cluster_center}-X_{pad};A/A_{sum}",i+1,layer), NbinDist, MinDist, MaxDist, NbinRatio, MinRatio, MaxRatio);
    }
  }
#endif

#if GainCorrection
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
    for(Int_t r=0; r<NumOfRow; ++r){
      HB1(TPCGainHid+(layer+1)*1000+r, Form("Cluster center dE layer%d row%d;Cluster center dE;Counts",layer,r), NbinHitDe, MinHitDe, MaxHitDe);
    }
  }
#endif

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
    HB1(TPCResYHid+layer*1000,
	Form("TPC Y Residual Layer%d (TPCHit);[mm];Counts", layer),
	NbinRes, MinRes, MaxRes);
    HB2(TPCResYClkHid+layer*1000,
	Form("ResY%%ClockTime Layer%d (TPCHit);[ns];[mm];Counts", layer),
	400, -50, 50, 400, -50, 50);
    HB2(TPCResYCoBoHid+layer*1000,
	Form("ResY%%ClockTime CoBo%d (TPCHit);[ns];[mm];Counts", layer),
	400, -50, 50, 400, -50, 50);
    HBProf(1+TPCResYCoBoHid+layer*1000,
	   //    	Form("TProf ResY%%ClockTime CoBo%d (TPCHit);[ns];[mm];Counts", layer),
	   Form("Tprof CoBo%d ResY%%ClockTime (TPCHit);[ns];[mm];Counts", layer),
    	400, -50, 50, -50, 50);
    HB2(TPCResYPosHid+layer*1000,
	Form("ResY%%YBcOut Layer%d (TPCHit);[mm];[mm];Counts", layer),
	200, -100, 100, 400, -50, 50);
    HBProf(1+TPCResYPosHid+layer*1000,
	   Form("Tprof Layer%d ResY%%ClockTime (TPCHit);[ns];[mm];Counts", layer),
	   400, -50, 50, -50, 50);
    HB1(TPCResXHid+layer*1000,
	Form("TPC X Residual Layer%d (TPCHit);[mm];Counts", layer),
	NbinRes, MinRes, MaxRes);
    HB2(TPCResXClkHid+layer*1000,
	Form("ResX%%ClockTime Layer%d (TPCHit);[ns];[mm];Counts", layer),
	400, -50, 50, 400, -50, 50);
    for(Int_t r=0; r<NumOfRow; ++r){
      // HB1(TPCPadYHid+layer*1000+r, "TPC Pad Y Cor", NbinPos, MinPos, MaxPos);
      HB1(TPCResYHid+layer*1000+r+1,
	  Form("TPC Y Residual L%dR%d (TPCHit);[mm];Counts", layer, r),
	  NbinRes, MinRes, MaxRes);
      // HB1(TPCPadXHid+layer*1000+r, "TPC Pad X Cor", NbinPos, MinPos, MaxPos);
      HB1(TPCResXHid+layer*1000+r+1,
	  Form("TPC X Residual L%dR%d (TPCHit);[mm];Counts", layer, r),
	  NbinRes, MinRes, MaxRes);
    }
  }

  HBTree("tpc", "tree of DstTPCTracking");
  tree->Branch("status", &event.status);
  tree->Branch("runnum", &event.runnum);
  tree->Branch("evnum", &event.evnum);
  tree->Branch("trigpat", &event.trigpat);
  tree->Branch("trigflag", &event.trigflag);

  tree->Branch("clkTpc", &event.clkTpc);

  tree->Branch("nhTpc", &event.nhTpc);
  tree->Branch("nclTpc", &event.nclTpc);
  tree->Branch("raw_hitpos_x", &event.raw_hitpos_x);
  tree->Branch("raw_hitpos_y", &event.raw_hitpos_y);
  tree->Branch("raw_hitpos_z", &event.raw_hitpos_z);
  tree->Branch("raw_de", &event.raw_de);
  tree->Branch("raw_padid", &event.raw_padid);
  tree->Branch("raw_layer", &event.raw_layer);
  tree->Branch("raw_row", &event.raw_row);
  tree->Branch("cluster_x", &event.cluster_x);
  tree->Branch("cluster_y", &event.cluster_y);
  tree->Branch("cluster_z", &event.cluster_z);
  tree->Branch("cluster_de", &event.cluster_de);
  tree->Branch("cluster_size", &event.cluster_size);
  tree->Branch("cluster_layer", &event.cluster_layer);
  tree->Branch("cluster_mrow", &event.cluster_mrow);
  tree->Branch("cluster_de_center", &event.cluster_de_center);
  tree->Branch("cluster_x_center", &event.cluster_x_center);
  tree->Branch("cluster_y_center", &event.cluster_y_center);
  tree->Branch("cluster_z_center", &event.cluster_z_center);
  tree->Branch("cluster_row_center", &event.cluster_row_center);

  tree->Branch("ntBcOut", &event.ntBcOut);
  tree->Branch("chisqrBcOut", &event.chisqrBcOut);
  tree->Branch("x0BcOut", &event.x0BcOut);
  tree->Branch("y0BcOut", &event.y0BcOut);
  tree->Branch("u0BcOut", &event.u0BcOut);
  tree->Branch("v0BcOut", &event.v0BcOut);

  tree->Branch("ntTpc", &event.ntTpc);
  tree->Branch("nhtrack", &event.nhtrack);
  tree->Branch("chisqrTpc", &event.chisqrTpc);
  tree->Branch("x0Tpc", &event.x0Tpc);
  tree->Branch("y0Tpc", &event.y0Tpc);
  tree->Branch("u0Tpc", &event.u0Tpc);
  tree->Branch("v0Tpc", &event.v0Tpc);
  tree->Branch("theta", &event.theta);
  tree->Branch("hitlayer", &event.hitlayer);
  tree->Branch("hitpos_x", &event.hitpos_x);
  tree->Branch("hitpos_y", &event.hitpos_y);
  tree->Branch("hitpos_z", &event.hitpos_z);
  tree->Branch("calpos_x", &event.calpos_x);
  tree->Branch("calpos_y", &event.calpos_y);
  tree->Branch("calpos_z", &event.calpos_z);
  tree->Branch("residual", &event.residual);
  tree->Branch("residual_x", &event.residual_x);
  tree->Branch("residual_y", &event.residual_y);
  tree->Branch("residual_wbcout_x", &event.residual_wbcout_x);
  tree->Branch("residual_wbcout_y", &event.residual_wbcout_y);
  tree->Branch("residual_trackwbcout_x", &event.residual_trackwbcout_x);
  tree->Branch("residual_trackwbcout_y", &event.residual_trackwbcout_y);
  tree->Branch("residual_z", &event.residual_z);
#if TrackCluster
  tree->Branch("track_cluster_de", &event.track_cluster_de);
  tree->Branch("track_cluster_size", &event.track_cluster_size);
  tree->Branch("track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch("track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch("track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch("track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch("track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch("track_cluster_row_center", &event.track_cluster_row_center);
#endif
#if PositionCorrection
  tree->Branch("xCorVec", &event.xCorVec);
  tree->Branch("yCorVec", &event.yCorVec);
  tree->Branch("zCorVec", &event.zCorVec);
  tree->Branch("xCorPos", &event.xCorPos);
  tree->Branch("yCorPos", &event.yCorPos);
  tree->Branch("zCorPos", &event.zCorPos);
#endif

  TTreeReaderCont[kTpcHit] = new TTreeReader("tpc", TFileCont[kTpcHit]);
  const auto& reader = TTreeReaderCont[kTpcHit];
  src.runnum = new TTreeReaderValue<Int_t>(*reader, "runnum");
  src.evnum = new TTreeReaderValue<Int_t>(*reader, "evnum");
  src.trigpat = new TTreeReaderValue<std::vector<Int_t>>(*reader, "trigpat");
  src.trigflag = new TTreeReaderValue<std::vector<Int_t>>(*reader, "trigflag");
  src.npadTpc = new TTreeReaderValue<Int_t>(*reader, "npadTpc");
  src.nhTpc = new TTreeReaderValue<Int_t>(*reader, "nhTpc");
  src.layerTpc = new TTreeReaderValue<std::vector<Int_t>>(*reader, "layerTpc");
  src.rowTpc = new TTreeReaderValue<std::vector<Int_t>>(*reader, "rowTpc");
  src.padTpc = new TTreeReaderValue<std::vector<Int_t>>(*reader, "padTpc");
  src.pedTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "pedTpc");
  src.rmsTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "rmsTpc");
  src.deTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "deTpc");
  src.tTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "tTpc");
  src.ctTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ctTpc");
  src.chisqrTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "chisqrTpc");
  src.clkTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "clkTpc");

  TTreeCont[kBcOut]->SetBranchStatus("*", 0);
  TTreeCont[kBcOut]->SetBranchStatus("ntrack",  1);
  TTreeCont[kBcOut]->SetBranchStatus("chisqr",  1);
  TTreeCont[kBcOut]->SetBranchStatus("x0",  1);
  TTreeCont[kBcOut]->SetBranchStatus("y0",  1);
  TTreeCont[kBcOut]->SetBranchStatus("u0",  1);
  TTreeCont[kBcOut]->SetBranchStatus("v0",  1);

  TTreeCont[kBcOut]->SetBranchAddress("ntrack",  &src.ntBcOut);
  TTreeCont[kBcOut]->SetBranchAddress("chisqr",  src.chisqrBcOut);
  TTreeCont[kBcOut]->SetBranchAddress("x0",  src.x0BcOut);
  TTreeCont[kBcOut]->SetBranchAddress("y0",  src.y0BcOut);
  TTreeCont[kBcOut]->SetBranchAddress("u0",  src.u0BcOut);
  TTreeCont[kBcOut]->SetBranchAddress("v0",  src.v0BcOut);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")   &&
     InitializeParameter<TPCParamMan>("TPCPRM") &&
     InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
     InitializeParameter<UserParamMan>("USER") &&
     InitializeParameter<HodoPHCMan>("HDPHC"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

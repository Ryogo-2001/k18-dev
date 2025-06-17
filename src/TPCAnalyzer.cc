// -*- C++ -*-

#include "TPCAnalyzer.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

//#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DCGeomRecord.hh"
#include "DebugCounter.hh"
#include "DebugTimer.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "RawData.hh"
#include "UserParamMan.hh"
#include "DeleteUtility.hh"
#include "RungeKuttaUtilities.hh"
#include "TPCPadHelper.hh"
#include "TPCParamMan.hh"
 //#include "TPCPositionCorrector.hh"
#include "TPCRawHit.hh"
#include "TPCHit.hh"
#include "TPCCluster.hh"
#include "TPCLocalTrack.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCTrackSearch.hh"
#include "DCHit.hh"
#include "TrackHit.hh"
#include "TPCRKTrack.hh"
#include "PrintHelper.hh"
#include "TRandom3.h"

 /* TPCTracking */
#define UseTpcCluster 1 // 1 : Common clustering method, 0 : Cluster size=1 no clustering

namespace
{
TRandom3 RandGen;

//const auto& gConf   = ConfMan::GetInstance();
const auto& gTPC  = TPCParamMan::GetInstance();
//const auto& gTPCPos = TPCPositionCorrector::GetInstance();
const auto& gUser   = UserParamMan::GetInstance();
const auto& gGeom   = DCGeomMan::GetInstance();

//Kurama
// const Int_t& IdTOF   = gGeom.DetectorId("TOF");
const Int_t& IdTgt = gGeom.DetectorId("Target");
const Int_t& IdTOFUX = gGeom.DetectorId("TOF-UX");
const Int_t& IdTOFUY = gGeom.DetectorId("TOF-UY");
const Int_t& IdTOFDX = gGeom.DetectorId("TOF-DX");
const Int_t& IdTOFDY = gGeom.DetectorId("TOF-DY");
const Int_t& IdTPCGasVessel_U = gGeom.DetectorId("VesselU");
const Int_t& IdTPCGasVessel_D = gGeom.DetectorId("VesselD");
const Int_t& IdVPHTOF = gGeom.DetectorId("VPHTOF");
const Int_t& IdRKINIT = gGeom.DetectorId("RKINIT");
//K1.8
const Int_t& IdBH2 = gGeom.DetectorId("BH2");
const Int_t& IdBAC = gGeom.DetectorId("BAC");

  /*
const Int_t& IdVP1 = gGeom.DetectorId("VPHS1");
const Int_t& IdVP2 = gGeom.DetectorId("VPHS2");
const Int_t& IdVP3 = gGeom.DetectorId("VPHS3");
const Int_t& IdVP4 = gGeom.DetectorId("VPHS4");
const Int_t& IdHTOF = gGeom.DetectorId("HTOF");
  */

const Double_t MaxChiSqrTrack = 1000.;
static std::vector<Double_t> gChisqr;
}

static inline Bool_t CompareChisqr(const Int_t a, const Int_t b){
  return gChisqr[a] < gChisqr[b];
}

//_____________________________________________________________________________
TPCAnalyzer::TPCAnalyzer()
  : m_is_decoded(n_type),
    m_TPCHitCont(NumOfLayersTPC+1),
    m_TPCClCont(NumOfLayersTPC)
{
  for(Int_t i=0; i<n_type; ++i){
    m_is_decoded[i] = false;
  }
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCAnalyzer::~TPCAnalyzer()
{
  ClearTPCHits();
  ClearTPCClusters();
  ClearTPCTracks();
  ClearTPCVertices();
  ClearTPCKuramaTracks();
  ClearTPCK18Tracks();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::MakeUpTPCClusters(const TPCHitContainer& HitCont,
			       TPCClusterContainer& ClCont,
			       Double_t maxdy)
{
  static const Double_t MinClusterDe = gUser.GetParameter("MinClusterDeTPC");
  static const Int_t MinClusterSize = gUser.GetParameter("MinClusterSizeTPC");
  static const Double_t MinClusterYPos = gUser.GetParameter("MinClusterYPosTPC");
  static const Double_t MaxClusterYPos = gUser.GetParameter("MaxClusterYPosTPC");

  const auto nh = HitCont.size();
  if(nh==0) return false;

  std::vector<Int_t> joined(nh, 0);
  for(Int_t i=0; i<nh; ++i){
    if(joined[i] > 0) continue;
    TPCHitContainer CandCont;
    TPCHit* hit = HitCont[i];
    if(!hit || !hit->IsGood()) continue;
    Int_t layer = hit->GetLayer();
    CandCont.push_back(hit);
    joined[i]++;
    Double_t padlength = hit -> GetPadLength();
    TVector3 dist2tgt = hit -> GetPosition() - TVector3(0., 0., tpc::ZTarget);
    Double_t verticalpathlength_forpad =
      padlength*dist2tgt.y()/TMath::Hypot(dist2tgt.x(), dist2tgt.z());
    maxdy = TMath::Max(maxdy, verticalpathlength_forpad);
#if UseTpcCluster
    for(Int_t j=0; j<nh; ++j){
      if(i==j || joined[j]>0) continue;
      TPCHit* thit = HitCont[j];
      if(!thit || !thit->IsGood()) continue;
      Int_t rowID = thit->GetRow();
      for(const auto& c_hit: CandCont){
        Int_t c_rowID = c_hit->GetRow();
        if(tpc::IsClusterable(layer, rowID, c_rowID)
           && TMath::Abs(thit->GetY() - c_hit->GetY()) < maxdy){
          CandCont.push_back(thit);
          joined[j]++;
          break;
        }
      }
    }
#endif

    TPCCluster* cluster = new TPCCluster(layer, CandCont);
    if(!cluster) continue;
    if(cluster->Calculate()
       && cluster->GetDe()>=MinClusterDe && cluster->GetClusterSize()>=MinClusterSize
       && cluster->GetY()>=MinClusterYPos && cluster->GetY()<=MaxClusterYPos){
      ClCont.push_back(cluster);
    }else{
      delete cluster;
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::DecodeTPCHits(RawData *rawData, Double_t clock)
{
  if(m_is_decoded[kTPC]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  ClearTPCHits();
  ClearTPCClusters();

  for(Int_t layer=0; layer<=NumOfLayersTPC; ++layer){
    for(const auto& rhit: rawData->GetTPCCorHC(layer)){
      auto hit = new TPCHit(rhit);
      if(hit->DoFit() && hit->Calculate(clock)){
        m_TPCHitCont[layer].push_back(hit);
      }else{
        delete hit;
      }
    }
  }

#if 0 // Cluster analysis will be done by RecalcTPCHits() in Dst.
  static const Double_t MaxYDif = gUser.GetParameter("MaxYDifClusterTPC");
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    MakeUpTPCClusters(m_TPCHitCont[layer], m_TPCClCont[layer], MaxYDif);
  }
#endif

  m_is_decoded[kTPC] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::ReCalcTPCHits(const Int_t nhits,
			   const std::vector<Int_t>& pad,
			   const std::vector<Double_t>& time,
			   const std::vector<Double_t>& de,
			   Double_t clock)
{
  if(m_is_decoded[kTPC]){
    hddaq::cerr << FUNC_NAME << " already decoded" << std::endl;
    return false;
  }

  static const Double_t MinCDe = gUser.GetParameter("MinCDeTPC");

  ClearTPCHits();
  ClearTPCClusters();

  if(nhits != pad.size() || nhits != time.size() || nhits != de.size()){
    hddaq::cerr << FUNC_NAME << " vector size mismatch" << std::endl;
    return false;
  }

  for(Int_t ih=0; ih<nhits; ih++){
    const Int_t layer = tpc::getLayerID(pad[ih]);
    const Double_t row = tpc::getRowID(pad[ih]);
    auto hit = new TPCHit(layer, row);
    hit->AddHit(de[ih], time[ih]);
    if(hit->Calculate(clock) && hit->GetCDe()>=MinCDe && hit->IsGood()){
      m_TPCHitCont[layer].push_back(hit);
    }else{
      delete hit;
    }
  }

#if 1
  static const Double_t MaxYDif = gUser.GetParameter("MaxYDifClusterTPC");
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    MakeUpTPCClusters(m_TPCHitCont[layer], m_TPCClCont[layer], MaxYDif);
  }
#endif

  m_is_decoded[kTPC] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::DecodeTPCHitsGeant4(const Int_t nhits,
				 const Double_t *x, const Double_t *y,
				 const Double_t *z, const Double_t *de,
				 const Int_t *pid, std::vector<TVector3> Mom)
{
  if(m_is_decoded[kTPC]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  static const Double_t MinClusterYPos = gUser.GetParameter("MinClusterYPosTPC");
  static const Double_t MaxClusterYPos = gUser.GetParameter("MaxClusterYPosTPC");
  static const Double_t MinCDe = gUser.GetParameter("MinCDeTPC");

  bool RejectKaonHits = false;
  ClearTPCHits();
  ClearTPCClusters();
  if(nhits != Mom.size()){
    hddaq::cerr << FUNC_NAME << " vector size mismatch" << std::endl;
    hddaq::cerr << FUNC_NAME << " nhits = " << nhits << " Mom.size() = " << Mom.size() << std::endl;
    Mom = std::vector<TVector3>(nhits, TVector3(0., 0., 0.));
  }
  for(Int_t i=0; i<nhits; i++){
    Int_t pad = tpc::findPadID(z[i], x[i]);
    if(pad<0) continue;
    Int_t layer = tpc::getLayerID(pad);
    Int_t row = tpc::getRowID(pad);
    if(RejectKaonHits && abs(pid[i])==321) continue;
    double CheckPad;
    gTPC.GetCDe(layer, row, 1,CheckPad);
    if(CheckPad==0) continue;
    TVector3 hitpos = TVector3(x[i], y[i], z[i]);
    double Eff = GetDetectionEfficiency(hitpos, pid[i], Mom[i], de[i]);
    double rndm = gRandom->Uniform(0., 1.);
    if(rndm > Eff) continue;
    if(de[i] == 0. || de[i] == TMath::QuietNaN() || de[i] < MinCDe) continue;
    auto hit = new TPCHit(layer, row);
    hit->AddHit(TMath::QuietNaN(), TMath::QuietNaN()); // allocate hit
    hit->SetDe(de[i]);
    // end of tentative treatment
    int cl_size = GetClusterSize(hitpos, pid[i], Mom[i], de[i]);
    if (cl_size == 1){
      auto PadPos = tpc::getPosition(pad);
      PadPos.SetY(y[i]);
      hit->SetPosition(PadPos);
//      hit->SetPosition(TVector3(x[i], y[i], z[i]));
    }
    else{
      hit->SetPosition(TVector3(x[i], y[i], z[i]));
    }

    TPCHitContainer CandCont;
    CandCont.push_back(hit);
    TPCCluster* cluster = new TPCCluster(layer, CandCont);
    if(!cluster) continue;
    if(cluster->Calculate() && cluster->GetY()>=MinClusterYPos && cluster->GetY()<=MaxClusterYPos){
      cluster->SetClusterSizeG4(cl_size);
      m_TPCClCont[layer].push_back(cluster);
    }
    else delete cluster;
    m_TPCHitCont[layer].push_back(hit);
  }

  m_is_decoded[kTPC] = true;
  return true;
}
Double_t
TPCAnalyzer::GetDetectionEfficiency(TVector3 pos, Int_t pid, TVector3 mom, Double_t de){
  return tpc::GetDetectionEfficiency(pos, pid, mom, de);
}
Int_t
TPCAnalyzer::GetClusterSize(TVector3 pos, Int_t pid, TVector3 mom, Double_t de){
  Int_t pad = tpc::findPadID(pos.z(), pos.x());
  Int_t layer = tpc::getLayerID(pad);
  Int_t row = tpc::getRowID(pad);
  bool Inner = false;
  if(layer < 10)Inner = true;
  double Mom = mom.Mag();
  double Cl1Prob = 0;
  Cl1Prob =tpc::GetClSize1Prob(Mom, pid, layer);
  if(abs(pid)==321) Cl1Prob= sqrt(tpc::GetClSize1Prob(Mom,211,layer)* tpc::GetClSize1Prob(Mom,2212,layer));
  int ncl=1;
  if(gRandom->Uniform(0., 1.) > Cl1Prob)ncl = 2;
  return ncl;
}


//_____________________________________________________________________________
Bool_t
TPCAnalyzer::TrackSearchTPC(Bool_t exclusive)
{
  if(m_is_decoded[kTPCTracking]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  static const Int_t MinLayer = gUser.GetParameter("MinLayerTPC");

  tpc::LocalTrackSearch(m_TPCClCont, m_TPCTC, m_TPCTCFailed, exclusive, MinLayer);

  m_is_decoded[kTPCTracking] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::TrackSearchTPCHelix(Bool_t exclusive)
{
  if(m_is_decoded[kTPCTracking]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  static const Int_t MinLayer = gUser.GetParameter("MinLayerTPC");
  tpc::LocalTrackSearchHelix(m_TPCClCont, m_TPCTCHelix, m_TPCTCHelixInverted, m_TPCTCHelixFailed, m_TPCVC, m_TPCVCClustered, exclusive, MinLayer);

  m_is_decoded[kTPCTracking] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::TrackSearchTPCHelix(std::vector<std::vector<TVector3>> K18VPs,
				 std::vector<std::vector<TVector3>> KuramaVPs,
				 std::vector<Double_t> KuramaCharge,
				 Bool_t exclusive)
{
  if(m_is_decoded[kTPCTracking]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  static const Int_t MinLayer = gUser.GetParameter("MinLayerTPC");
  tpc::LocalTrackSearchHelix(K18VPs, KuramaVPs, KuramaCharge, m_TPCClCont, m_TPCTCHelix, m_TPCTCHelixInverted, m_TPCTCVP, m_TPCTCHelixFailed, m_TPCVC, m_TPCVCClustered, exclusive, MinLayer);

  m_is_decoded[kTPCTracking] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::TestHoughTransform()
{
  if(m_is_decoded[kTPCTracking]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  static const Int_t MinLayer = gUser.GetParameter("MinLayerTPC");

  tpc::HoughTransformTest(m_TPCClCont, m_TPCTC, MinLayer);

  m_is_decoded[kTPCTracking] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::TestHoughTransformHelix()
{
  if(m_is_decoded[kTPCTracking]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  static const Int_t MinLayer = gUser.GetParameter("MinLayerTPC");

  tpc::HoughTransformTestHelix(m_TPCClCont, m_TPCTCHelix, MinLayer);

  m_is_decoded[kTPCTracking] = true;
  return true;
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCHits()
{
  del::ClearContainerAll(m_TPCHitCont);
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCClusters()
{
  del::ClearContainerAll(m_TPCClCont);
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCTracks()
{
  del::ClearContainer(m_TPCTC);
  del::ClearContainer(m_TPCTCFailed);
  del::ClearContainer(m_TPCTCHelix);
  del::ClearContainer(m_TPCTCHelixInverted);
  del::ClearContainer(m_TPCTCVP);
  del::ClearContainer(m_TPCTCHelixFailed);
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCVertices()
{
  del::ClearContainer(m_TPCVC);
  del::ClearContainer(m_TPCVCClustered);
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCKuramaTracks()
{
  del::ClearContainer(m_TPCKuramaTC);
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCK18Tracks()
{
  del::ClearContainer(m_TPCK18TC);
}

//_____________________________________________________________________________
//KuramaTracking
void
TPCAnalyzer::MakeTPCKuramaHPContainer(TPCLocalTrackHelix *tpctrack,
				      std::vector<Int_t> &lnum){

  //The end point of RK tracking is the target z position.
  //So, clusters before the target center are excluded in the tracking. (Z_cluster < Z_target).

  lnum.clear();

  //Virtual planes on the TPC clusters
  Int_t nhtpc = tpctrack -> GetNHit();
  for(int id=0; id<nhtpc; ++id){
    TPCLTrackHit *thp = tpctrack -> GetHitInOrder(id);
    if(!thp) continue;
    //Exclude bad hits
    const TVector3& resolution = thp->GetResolutionVect();
    if(resolution.x() > 0.9e+10 && resolution.y() > 0.9e+10 && resolution.z() > 0.9e+10) continue;
    //Exclude clusters before the target
    const TVector3& localhitpos = thp->GetLocalHitPos();
    Double_t TGTz = gGeom.GlobalZ("Target") - gGeom.GlobalZ("HS");
    if(localhitpos.z()<TGTz){
      tpctrack -> SetIsBeam();
      continue;
    }
    lnum.push_back(id + PlOffsTPCHit + 1);
    lnum.push_back(id + 2.*PlOffsTPCHit + 1);
  }

  //Virtual planes of Kurama detectors
  // /*** From Upstream ***/
  lnum.push_back(IdTgt);
  for(int j=0; j<NumOfLayersVPTPC; ++j){
    const Int_t& IdVP = gGeom.DetectorId(Form("VPTPC%d",j+1));
    lnum.push_back(IdVP);
  }
  lnum.push_back(IdTPCGasVessel_D);
  lnum.push_back(IdVPHTOF);
  for(Int_t i=0; i<NumOfLayersSdcIn; ++i) lnum.push_back(i+PlOffsSdcIn+1);
  for(Int_t i=0; i<NumOfLayersVP; ++i) lnum.push_back(i+PlOffsVP+1);
  for(Int_t i=0; i<NumOfLayersSdcOut; ++i) lnum.push_back(i+PlOffsSdcOut+1);
  lnum.push_back(IdTOFUX);
  lnum.push_back(IdTOFUY);
  lnum.push_back(IdTOFDX);
  lnum.push_back(IdTOFDY);

}

//_____________________________________________________________________________
TPCRKTrack*
TPCAnalyzer::MakeTPCKuramaTrack(TPCLocalTrackHelix *tpctrack,
				const std::vector<Int_t> layerID,
				const std::vector<Double_t> wire,
				const std::vector<Double_t> localHitPos)
{

  std::vector<Int_t> lnum;
  MakeTPCKuramaHPContainer(tpctrack, lnum);

  TPCRKTrack *track = new TPCRKTrack(tpctrack, lnum);
  for(Int_t ih=0; ih<layerID.size(); ih++){ //Add DC Hits
    auto hit = new DCHit(layerID[ih], wire[ih]);
    hit -> SetWirePosition(gGeom.CalcWirePosition(layerID[ih], wire[ih]));
    hit -> SetTiltAngle(gGeom.GetTiltAngle(layerID[ih]));
    hit -> SetZ(gGeom.GetLocalZ(layerID[ih]));
    auto dclthit = new DCLTrackHit(hit, localHitPos[ih], ih);
    dclthit -> SetHoneycomb();
    auto thit = new TrackHit(dclthit);
    track -> AddRKHit(thit);
    if(layerID[ih]==IdTOFUX || layerID[ih]==IdTOFDX) track -> SetTofSeg(wire[ih]);
  } //nhitsDC

  return track;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::DoFitTPCKuramaTrack(Int_t KuramaTrackID,
				 const Int_t PID,
				 const TVector3 initPos,
				 const TVector3 initMom,
				 const std::vector<Int_t> layerID,
				 const std::vector<Double_t> wire,
				 const std::vector<Double_t> localHitPos,
				 Bool_t U2D)
{
  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);

  Bool_t status = false;
  gChisqr.clear();
  std::vector<Int_t> id;
  Int_t ntTpc = GetNTracksTPCHelix();
  TPCRKTrackContainer tempKuramaTC;
  for(Int_t ittpc=0; ittpc<ntTpc; ittpc++){
    TPCLocalTrackHelix* tpctrack = GetTrackTPCHelix(ittpc);
    if(tpctrack -> GetTrackID()!=-1) continue; //if the track already has a Kurama track pair
    if(tpctrack -> GetIsKurama()!=1 || !tpctrack -> IsKuramaTrackCandidate(KuramaTrackID)) continue; //whether the track is pair candidate of the kurama track
    if(!BeamThroughTPC && !tpctrack -> VertexAtTarget()){
      tpctrack -> SetIsKurama(0);
      continue;
    }

    //After track pair matching
    TPCRKTrack* track = MakeTPCKuramaTrack(tpctrack, layerID, wire, localHitPos);
    track -> SetPID(PID);
    track -> SetTPCTrackID(ittpc);
    if(track->DoFit(initPos, initMom, U2D) && track->ChiSquare()<MaxChiSqrTrack){
      tempKuramaTC.push_back(track);
      gChisqr.push_back(track -> ChiSquare());
      status = true;
      id.push_back(gChisqr.size() - 1);
    }
    else delete track;
  } //nttpc

  //find the best combination of lowest ChiSqruare()
  std::sort(id.begin(), id.end(), CompareChisqr);
  if(tempKuramaTC.size()!=0){
    tempKuramaTC.at(id.at(0)) -> SetTrackID(KuramaTrackID);
    m_TPCKuramaTC.push_back(tempKuramaTC.at(id.at(0)));
    tempKuramaTC.erase(tempKuramaTC.begin()+id.at(0));
  }
  for(auto& track: tempKuramaTC) track -> GetTPCTrack() -> SetIsKurama(0);
  del::ClearContainer(tempKuramaTC);

  return status;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::TrackSearchTPCKurama(const std::vector<Int_t> PID,
				  const std::vector<TVector3> initPos,
				  const std::vector<TVector3> initMom,
				  const std::vector<std::vector<Int_t>> layerID,
				  const std::vector<std::vector<Double_t>> wire,
				  const std::vector<std::vector<Double_t>> localHitPos,
				  Bool_t U2D)
{

  if(m_is_decoded[kTPCKurama]){
    hddaq::cerr << FUNC_NAME << " already decoded" << std::endl;
    return false;
  }

  if(initPos.size()!=initMom.size() || initPos.size()!=layerID.size() ||
     initPos.size()!=wire.size() || initPos.size()!=localHitPos.size()){
    hddaq::cerr << FUNC_NAME << " sizes of vector are not consistent" << std::endl;
    return false;
  }

  ClearTPCKuramaTracks();
  Int_t ntKurama = initPos.size();
  for(Int_t itkurama=0; itkurama<ntKurama; itkurama++){
    DoFitTPCKuramaTrack(itkurama, PID[itkurama], initPos[itkurama], initMom[itkurama],
			layerID[itkurama], wire[itkurama],
			localHitPos[itkurama], U2D);
  } //ntkurama

  m_is_decoded[kTPCKurama] = true;
  return true;
}

//_____________________________________________________________________________
//K18Tracking
void
TPCAnalyzer::MakeTPCK18HPContainer(TPCLocalTrackHelix *tpctrack,
				   std::vector<Int_t> &lnum){
  lnum.clear();

  //Virtual planes on the TPC clusters
  Int_t nhtpc = tpctrack -> GetNHit();
  for(int id=0; id<nhtpc; ++id){
    TPCLTrackHit *thp = tpctrack -> GetHitInOrder(id);
    if(!thp) continue;
    //Exclude bad hits
    const TVector3& resolution = thp->GetResolutionVect();
    if(resolution.x() > 0.9e+10 && resolution.y() > 0.9e+10 && resolution.z() > 0.9e+10) continue;
    //Exclude clusters after the target
    const TVector3& localhitpos = thp->GetLocalHitPos();
    Double_t TGTz = gGeom.GlobalZ("Target") - gGeom.GlobalZ("HS");
    if(localhitpos.z()>TGTz) continue;
    lnum.push_back(id + PlOffsTPCHit + 1);
    lnum.push_back(id + 2.*PlOffsTPCHit + 1);
  }

  //Virtual planes of Kurama detectors
  // /*** From Upstream ***/
  for(Int_t i=0; i<6; ++i) lnum.push_back(i+PlOffsBcOut+1);
  lnum.push_back(IdBAC);
  for(Int_t i=6; i<NumOfLayersBcOut; ++i) lnum.push_back(i+PlOffsBcOut+1);
  lnum.push_back(IdBH2);
  lnum.push_back(IdTPCGasVessel_U);
  for(int j=0; j<NumOfLayersVPHS; ++j){
    const Int_t& IdVP = gGeom.DetectorId(Form("VPHS%d",j+1));
    lnum.push_back(IdVP);
  }
  lnum.push_back(IdTgt);
  lnum.push_back(IdTPCGasVessel_D);
  lnum.push_back(IdVPHTOF);

}

//_____________________________________________________________________________
TPCRKTrack*
TPCAnalyzer::MakeTPCK18Track(TPCLocalTrackHelix *tpctrack,
			     const std::vector<Int_t> layerID,
			     const std::vector<Double_t> wire,
			     const std::vector<Double_t> localHitPos)
{

  std::vector<Int_t> lnum;
  MakeTPCK18HPContainer(tpctrack, lnum);
  TPCRKTrack *track = new TPCRKTrack(tpctrack, lnum);
  for(Int_t ih=0; ih<layerID.size(); ih++){ //Add DC Hits
    auto hit = new DCHit(layerID[ih], wire[ih]);
    hit -> SetWirePosition(gGeom.CalcWirePosition(layerID[ih], wire[ih]));
    hit -> SetTiltAngle(gGeom.GetTiltAngle(layerID[ih]));
    hit -> SetZ(gGeom.GetLocalZ(layerID[ih]));
    auto dclthit = new DCLTrackHit(hit, localHitPos[ih], ih);
    auto thit = new TrackHit(dclthit);
    track -> AddRKHit(thit);
  } //nhitsDC

  return track;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::DoFitTPCK18Track(Int_t K18TrackID,
			      const TVector3 initPos,
			      const TVector3 initMom,
			      const std::vector<Int_t> layerID,
			      const std::vector<Double_t> wire,
			      const std::vector<Double_t> localHitPos,
			      Bool_t U2D)
{
  Bool_t status = false;
  gChisqr.clear();
  std::vector<Int_t> id;
  Int_t ntTpc = GetNTracksTPCHelix();
  TPCRKTrackContainer tempK18TC;
  for(Int_t ittpc=0; ittpc<ntTpc; ittpc++){
    TPCLocalTrackHelix* tpctrack = GetTrackTPCHelix(ittpc);
    if(K18TrackID != tpctrack -> GetTrackID() || tpctrack -> GetIsK18()!=1) continue;
    //After track matching of a K18Track and a TPC track
    TPCRKTrack* track = MakeTPCK18Track(tpctrack, layerID, wire, localHitPos);
    if(track->DoFit(initPos, initMom, U2D) && track->ChiSquare()<MaxChiSqrTrack){
      tempK18TC.push_back(track);
      gChisqr.push_back(track -> ChiSquare());
      track -> SetTPCTrackID(ittpc);
      status = true;
      id.push_back(gChisqr.size() - 1);
    }
    else delete track;
  } //nttpc

  //find the best combination of lowest ChiSqruare()
  std::sort(id.begin(), id.end(), CompareChisqr);
  if(tempK18TC.size()!=0){
    tempK18TC.at(id.at(0)) -> SetTrackID(K18TrackID);
    m_TPCK18TC.push_back(tempK18TC.at(id.at(0)));
    tempK18TC.erase(tempK18TC.begin()+id.at(0));
  }
  for(auto& track: tempK18TC) track -> GetTPCTrack() -> SetIsK18(0);
  del::ClearContainer(tempK18TC);

  return status;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::TrackSearchTPCK18(const std::vector<TVector3> initPos,
			       const std::vector<TVector3> initMom,
			       const std::vector<std::vector<Int_t>> layerID,
			       const std::vector<std::vector<Double_t>> wire,
			       const std::vector<std::vector<Double_t>> localHitPos,
			       Bool_t U2D)
{

  if(m_is_decoded[kTPCK18]){
    hddaq::cerr << FUNC_NAME << " already decoded" << std::endl;
    return false;
  }

  if(initPos.size()!=initMom.size() || initPos.size()!=layerID.size() ||
     initPos.size()!=wire.size() || initPos.size()!=localHitPos.size()){
    hddaq::cerr << FUNC_NAME << " sizes of vector are not consistent" << std::endl;
    return false;
  }

  ClearTPCK18Tracks();
  Int_t ntK18 = initPos.size();
  for(Int_t itk18=0; itk18<ntK18; itk18++){
    DoFitTPCK18Track(itk18, initPos[itk18], initMom[itk18],
		     layerID[itk18], wire[itk18],
		     localHitPos[itk18], U2D);
  } //ntK18

  m_is_decoded[kTPCK18] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::GetProductionVertex(const TVector3 XtgtKm, const TVector3 PtgtKm,
				 const TVector3 XtgtKp, const TVector3 PtgtKp,
				 TVector3 &Vertex, Double_t &closeDist,
				 Double_t &KmPathInTgt, Double_t &KpPathInTgt,
				 TVector3 &KmMomVtx, TVector3 &KpMomVtx)
{

  RKcalcHitPoint inPointKm; RKcalcHitPoint outPointKp;
  return RK::FindVertex(XtgtKm, PtgtKm, XtgtKp, PtgtKp, Vertex, closeDist, KmPathInTgt, KpPathInTgt, KmMomVtx, KpMomVtx, inPointKm, outPointKp);
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::ReCalcTPCTracks(const Int_t ntracks,
			     const std::vector<Int_t>& nhits,
			     const std::vector<Double_t>& x0,
			     const std::vector<Double_t>& y0,
			     const std::vector<Double_t>& u0,
			     const std::vector<Double_t>& v0,
			     const std::vector<std::vector<Double_t>>& layer,
			     const std::vector<std::vector<Double_t>>& mrow,
			     const std::vector<std::vector<Double_t>>& de,
			     const std::vector<std::vector<Double_t>>& res_x,
			     const std::vector<std::vector<Double_t>>& res_y,
			     const std::vector<std::vector<Double_t>>& res_z,
			     const std::vector<std::vector<Double_t>>& localpos_x,
			     const std::vector<std::vector<Double_t>>& localpos_y,
			     const std::vector<std::vector<Double_t>>& localpos_z)
{
  if(m_is_decoded[kTPC]){
    hddaq::cerr << FUNC_NAME << " already decoded" << std::endl;
    return false;
  }

  ClearTPCHits();
  //ClearTPCClusters();
  ClearTPCTracks();
  ClearTPCVertices();

  if(ntracks != nhits.size() ||
     ntracks != x0.size() ||
     ntracks != y0.size() ||
     ntracks != u0.size() ||
     ntracks != v0.size() ||
     ntracks != layer.size() ||
     ntracks != mrow.size() ||
     ntracks != res_x.size() ||
     ntracks != res_y.size() ||
     ntracks != res_z.size() ||
     ntracks != localpos_x.size() ||
     ntracks != localpos_y.size() ||
     ntracks != localpos_z.size() ||
     ntracks != de.size()){

    hddaq::cerr << FUNC_NAME << " track params vector size mismatch" << std::endl;
    return false;
  }

  for(Int_t it=0; it<ntracks; it++){
    if(layer[it].size() != nhits[it] ||
       mrow[it].size() != nhits[it] ||
       res_x[it].size() != nhits[it] ||
       res_y[it].size() != nhits[it] ||
       res_z[it].size() != nhits[it] ||
       localpos_x[it].size() != nhits[it] ||
       localpos_y[it].size() != nhits[it] ||
       localpos_z[it].size() != nhits[it] ||
       de[it].size() != nhits[it]){

      hddaq::cerr << FUNC_NAME << " track params vector size mismatch" << std::endl;
      return false;
    }

    TPCLocalTrack *track = new TPCLocalTrack();
    Double_t Par[4] = {x0[it], y0[it], u0[it], v0[it]};
    track->SetParam(Par);
    for(Int_t ih=0; ih<nhits[it]; ih++){
      Int_t l = layer[it][ih];
      auto hit = new TPCHit(l, mrow[it][ih]);
      hit -> AddHit(0, 0.);
      hit -> SetDe(de[it][ih]);
      hit -> SetPosition(TVector3(localpos_x[it][ih], localpos_y[it][ih], localpos_z[it][ih]));
      m_TPCHitCont[l].push_back(hit);

      TPCLTrackHit *hitp = new TPCLTrackHit(hit);
      track->AddTPCHit(hitp);
    }
    track -> RecalcTrack();
    m_TPCTC.push_back(track);
  }

  m_is_decoded[kTPC] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::ReCalcTPCTracks(const Int_t ntracks,
			     const std::vector<Int_t>& isK18,
			     const std::vector<Int_t>& isKurama,
			     const std::vector<Int_t>& charge,
			     const std::vector<Int_t>& nhits,
			     const std::vector<Double_t>& cx,
			     const std::vector<Double_t>& cy,
			     const std::vector<Double_t>& z0,
			     const std::vector<Double_t>& r,
			     const std::vector<Double_t>& dz,
			     const std::vector<std::vector<Double_t>>& layer,
			     const std::vector<std::vector<Double_t>>& mrow,
			     const std::vector<std::vector<Double_t>>& helix_t,
			     const std::vector<std::vector<Double_t>>& de,
			     const std::vector<std::vector<Double_t>>& res_x,
			     const std::vector<std::vector<Double_t>>& res_y,
			     const std::vector<std::vector<Double_t>>& res_z,
			     const std::vector<std::vector<Double_t>>& localpos_x,
			     const std::vector<std::vector<Double_t>>& localpos_y,
			     const std::vector<std::vector<Double_t>>& localpos_z)
{
  if(m_is_decoded[kTPC]){
    hddaq::cerr << FUNC_NAME << " already decoded" << std::endl;
    return false;
  }

  ClearTPCHits();
  //ClearTPCClusters();
  ClearTPCTracks();
  ClearTPCVertices();

  if(ntracks != nhits.size() ||
     ntracks != cx.size() ||
     ntracks != cy.size() ||
     ntracks != z0.size() ||
     ntracks != r.size() ||
     ntracks != dz.size() ||
     ntracks != isK18.size() ||
     ntracks != isKurama.size() ||
     ntracks != charge.size() ||
     ntracks != layer.size() ||
     ntracks != mrow.size() ||
     ntracks != helix_t.size() ||
     ntracks != res_x.size() ||
     ntracks != res_y.size() ||
     ntracks != res_z.size() ||
     ntracks != localpos_x.size() ||
     ntracks != localpos_y.size() ||
     ntracks != localpos_z.size() ||
     ntracks != de.size()){

    hddaq::cerr << FUNC_NAME << " track params vector size mismatch" << std::endl;
    return false;
  }

  for(Int_t it=0; it<ntracks; it++){
    if(layer[it].size() != nhits[it] ||
       mrow[it].size() != nhits[it] ||
       helix_t[it].size() != nhits[it] ||
       res_x[it].size() != nhits[it] ||
       res_y[it].size() != nhits[it] ||
       res_z[it].size() != nhits[it] ||
       localpos_x[it].size() != nhits[it] ||
       localpos_y[it].size() != nhits[it] ||
       localpos_z[it].size() != nhits[it] ||
       de[it].size() != nhits[it]){

      hddaq::cerr << FUNC_NAME << " track params vector size mismatch" << std::endl;
      return false;
    }

    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    Double_t HelixPar[5] = {cx[it], cy[it], z0[it], r[it], dz[it]};
    track -> SetParam(HelixPar);
    track -> SetIsK18(isK18[it]);
    track -> SetIsKurama(isKurama[it]);
    track -> SetCharge(charge[it]);
    for(Int_t ih=0; ih<nhits[it]; ih++){
      Int_t id = ih;
      if(charge[it] < 0) id = nhits[it] -ih -1;
      Int_t l = layer[it][id];
      auto hit = new TPCHit(l, mrow[it][id]);
      hit -> AddHit(0, 0.);
      hit -> SetIsGood();
      hit -> SetDe(de[it][id]);
      hit -> SetPosition(TVector3(localpos_x[it][id], localpos_y[it][id], localpos_z[it][id]));
      m_TPCHitCont[l].push_back(hit);

      TPCLTrackHit *hitp = new TPCLTrackHit(hit);
      hitp -> SetResolution(TVector3(res_x[it][id], res_y[it][id], res_z[it][id]));
      hitp -> SetTheta(helix_t[it][id]);
      track -> AddTPCHit(hitp);
    }
    track -> RecalcTrack();
    m_TPCTCHelix.push_back(track);
  }

  m_is_decoded[kTPC] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::ReCalcTPCTracksGeant4(const Int_t ntracks,
				   const std::vector<Int_t>& charge,
				   const std::vector<Int_t>& nhits,
				   const std::vector<Double_t>& cx,
				   const std::vector<Double_t>& cy,
				   const std::vector<Double_t>& z0,
				   const std::vector<Double_t>& r,
				   const std::vector<Double_t>& dz,
				   const std::vector<std::vector<Double_t>>& layer,
				   const std::vector<std::vector<Double_t>>& mrow,
				   const std::vector<std::vector<Double_t>>& helix_t,
				   const std::vector<std::vector<Double_t>>& de,
				   const std::vector<std::vector<Double_t>>& res_x,
				   const std::vector<std::vector<Double_t>>& res_y,
				   const std::vector<std::vector<Double_t>>& res_z,
				   const std::vector<std::vector<Double_t>>& localpos_x,
				   const std::vector<std::vector<Double_t>>& localpos_y,
				   const std::vector<std::vector<Double_t>>& localpos_z)
{
  if(m_is_decoded[kTPC]){
    hddaq::cerr << FUNC_NAME << " already decoded" << std::endl;
    return false;
  }

  ClearTPCHits();
  //ClearTPCClusters();
  ClearTPCTracks();
  ClearTPCVertices();

  if(ntracks != nhits.size() ||
     ntracks != cx.size() ||
     ntracks != cy.size() ||
     ntracks != z0.size() ||
     ntracks != r.size() ||
     ntracks != dz.size() ||
     ntracks != charge.size() ||
     ntracks != layer.size() ||
     ntracks != mrow.size() ||
     ntracks != helix_t.size() ||
     ntracks != res_x.size() ||
     ntracks != res_y.size() ||
     ntracks != res_z.size() ||
     ntracks != localpos_x.size() ||
     ntracks != localpos_y.size() ||
     ntracks != localpos_z.size() ||
     ntracks != de.size()){

    hddaq::cerr << FUNC_NAME << " track params vector size mismatch" << std::endl;
    return false;
  }

  for(Int_t it=0; it<ntracks; it++){
    if(layer[it].size() != nhits[it] ||
       mrow[it].size() != nhits[it] ||
       helix_t[it].size() != nhits[it] ||
       res_x[it].size() != nhits[it] ||
       res_y[it].size() != nhits[it] ||
       res_z[it].size() != nhits[it] ||
       localpos_x[it].size() != nhits[it] ||
       localpos_y[it].size() != nhits[it] ||
       localpos_z[it].size() != nhits[it] ||
       de[it].size() != nhits[it]){

      hddaq::cerr << FUNC_NAME << " track params vector size mismatch" << std::endl;
      return false;
    }

    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    Double_t HelixPar[5] = {cx[it], cy[it], z0[it], r[it], dz[it]};
    track -> SetParam(HelixPar);
    track -> SetCharge(charge[it]);
    for(Int_t ih=0; ih<nhits[it]; ih++){
      Int_t id = ih;
      if(charge[it] < 0) id = nhits[it] -ih -1;
      Int_t l = layer[it][id];
      auto hit = new TPCHit(l, mrow[it][id]);
      hit -> AddHit(0, 0.);
      hit -> SetIsGood();
      hit -> SetDe(de[it][id]);
      hit -> SetPosition(TVector3(localpos_x[it][id], localpos_y[it][id], localpos_z[it][id]));
      m_TPCHitCont[l].push_back(hit);

      TPCLTrackHit *hitp = new TPCLTrackHit(hit);
      hitp -> SetResolution(TVector3(res_x[it][id], res_y[it][id], res_z[it][id]));
      hitp -> SetTheta(helix_t[it][id]);
      track -> AddTPCHit(hitp);
    }
    track -> RecalcTrack();
    m_TPCTCHelix.push_back(track);
  }

  m_is_decoded[kTPC] = true;
  return true;
}

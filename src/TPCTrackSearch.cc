// -*- C++ -*-

/*
//Comment by Wooseung

Track finding process
1. Hough-transform on the XZ plane
2. Hough-transform on the vertical plane
   -> Get initial track parameters.
3. Check residual(Houghdist) in the Hough-space and
make a inital track with clusters with small residual.

Track fitting process
1. First fitting with inital track.
2. If fitting is succeeded, check the residual of other cluster with the track.
If residual is acceptable, add cluster into the track.
3. If the track is extended by added clusters, then do fitting again. (recursive)

Especially near the target, two almost parallel tracks can be reconiged as a single track. So, there is a treatment to separate them.
track -> SeparateTracksAtTraget();

LocalTrackSearch/Helix are main functions of TPC tracking.
HoughTransformTest/Helix functions have only track finding algorithm for Hough-transform performace test. (no fitting)

//Veto process for accidental coincidence event
MarkingAccidentalTracks : find accidental beams
FindAccidentalCoincidenceTracks &&  MarkingClusteredAccidentalTracks: find accidental cpincidence track cluster for the reacted accidental event(not beam-thorugh)

Detailed fitting procedures are explained in the TPCLocalTrack/Helix.
*/

#include "TPCTrackSearch.hh"

#include <chrono>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <TH2D.h>
#include <TH3D.h>
#include <TLorentzVector.h>

#include "DebugTimer.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "DatabasePDG.hh"
#include "UserParamMan.hh"
#include "DeleteUtility.hh"
#include "ConfMan.hh"
#include "HoughTransform.hh"
#include "Kinematics.hh"
#include "TPCPadHelper.hh"
#include "TPCLTrackHit.hh"
#include "TPCLocalTrack.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCVertex.hh"
#include "TPCCluster.hh"
#include "RootHelper.hh"
#include "DebugCounter.hh"
#define DebugDisp 0
#define FragmentedTrackTest 1
#define ReassignClusterTest 1
//#define RemainingClustersTest 1
#define RemainingClustersTest 0 //not helpful
#define RefitXiTrack 1

namespace
{
  const auto qnan = TMath::QuietNaN();
  const auto& gUser = UserParamMan::GetInstance();
  const auto& gCounter = debug::ObjectCounter::GetInstance();
  const Int_t    MaxNumOfTrackTPC = 30;
  //const Double_t KuramaXZWindow = 35.; //ref
  const Double_t KuramaXZWindow = 45.;
  //const Double_t KuramaYWindow = 20.; //ref
  const Double_t KuramaYWindow = 30.; //ref
  const Double_t K18XZWindow = 10.5;
  //const Double_t K18YWindow = 10.;
  const Double_t K18YWindow = 15.;

  // Closest distance cut for vertex finding
  const Double_t VertexDistCut = 30.;
  const Double_t ppi_distcut = 10.; //Closest distance for p, pi at the vertex point

  // Maximum number of fitting steps
  //const Int_t MaxFitSteps = 5;
  const Int_t MaxFitSteps = 10;

  // Houghflags
  const Int_t GoodForTracking = 100;
  const Int_t K18Tracks = 200;
  //const Int_t KuramaTracks = 300;
  const Int_t BadHoughTransform = 300;
  const Int_t BadForTracking = 400;
  const Int_t Candidate = 1000;

  // Tracks in the Hough-Space
  std::vector<Double_t> XZhough_x;
  std::vector<Double_t> XZhough_y;
  std::vector<Double_t> XZhough_z;
  std::vector<Double_t> Yhough_x;
  std::vector<Double_t> Yhough_y;

  // B-field
  const Double_t& HS_field_0 = ConfMan::Get<Double_t>("HSFLDCALIB");
  const Double_t& HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
  const Double_t& HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");

  //_____________________________________________________________________________
  // Local Functions
  //_____________________________________________________________________________
  template <typename T> void
  CalcTracks(std::vector<T*>& TrackCont)
  {
    for(auto& track: TrackCont){
      track->Calculate();
    }
  }

  //_____________________________________________________________________________
  template <typename T> void
  ExclusiveTracking(std::vector<T*>& TrackCont)
  {
    for(auto& track: TrackCont){
      track->DoFitExclusive();
      track->CalculateExclusive();
    }
  }

  template <typename T> void
  MarkingAccidentalTracks(std::vector<T*>& TrackCont)
  {
    for(auto& track: TrackCont){
      if(track->GetIsAccidental()!=1) track->CheckIsAccidental();
    }
  }

  template <typename T> void
  MarkingClusteredAccidentalTracks(std::vector<T*>& TrackCont, std::vector<TPCVertex*>& ClusteredVertexCont)
  {
    Double_t signal_section = 30; //within this section, we don't use this veto process.
    for(auto& vertex: ClusteredVertexCont){
      if(TMath::Abs(vertex -> GetVertex().y()) < signal_section) continue;
      vertex -> SetIsAccidental();

      Int_t ntracks = vertex -> GetNTracks();
      for(Int_t i=0; i<ntracks; i++){
	Int_t id = vertex -> GetTrackId(i);
	TrackCont[id] -> SetIsAccidental();
      }
    }
  }

  //_____________________________________________________________________________
  //reset houghflag of remain clusters
  void
  ResetHoughFlag(const std::vector<TPCClusterContainer>& ClCont, Int_t flagID=Candidate)
  {

    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	auto cl = ClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	if(hit->GetHoughFlag()==flagID) hit->SetHoughFlag(0);
      } //ci
    } //layer
  }

  //_____________________________________________________________________________
  template <typename T> void
  GetTrackClCont(T* track, std::vector<TVector3>& gHitPos)
  {

    Int_t n = track->GetNHit();
    for(Int_t i=0; i<n; ++i){
      TPCLTrackHit *hitp = track->GetHit(i);
      TVector3 pos = hitp->GetLocalHitPos();
      gHitPos.push_back(pos);
    }
  }

  //_____________________________________________________________________________
  template <typename T> Bool_t
  IsGood(T* track, Int_t MinNumOfHits)
  {

    Bool_t status = false;
    if(track->IsBackward()) status = true;
    else if(track->GetNHit() >= MinNumOfHits) status = true; //track w/ enough clusters
    return status;
  }

  //_____________________________________________________________________________
  //Convert a straight line track into helix track
  Bool_t
  ConvertTrack(TPCLocalTrack *LinearTrack, TPCLocalTrackHelix *HelixTrack)
  {

    Int_t n = LinearTrack->GetNHit();
    for(Int_t i=0; i<n; ++i){
      TPCLTrackHit *hitp = LinearTrack->GetHit(i);
      HelixTrack->AddTPCHit(hitp);
    }
    Double_t LinearPar[4];
    LinearTrack->GetParam(LinearPar);
    Bool_t status = HelixTrack->ConvertParam(LinearPar);
    delete LinearTrack;
    if(status &&
       HelixTrack-> Getr() < 200.) status = false;
    return status;
  }

  //_____________________________________________________________________________
  template <typename T> Bool_t
  AddClusters(T* Track, const std::vector<TPCClusterContainer>& ClCont, Int_t HoughFlag=0)
  {

    Bool_t status = false;
    //Residual check with other hits
    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	auto cl = ClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	if(hit->GetHoughFlag()!=HoughFlag) continue;
	Double_t resi=0.;
	if(Track->IsGoodHitToAdd(hit, resi)){
	  Int_t vtxflag = Track -> GetVtxFlag();
	  TVector3 pos = hit -> GetPosition();
	  Int_t side = Track -> Side(pos);
	  //Vertex inside the target : vtxflag = -1 or 1 / outside vtxflag = 0
	  //if track and new cluster are on the same side and vertex in the target : vtxflag*side = 1
	  //if vertex is outside of the target : vtxflag*side = 0
	  if(vtxflag*side >= 0){
	    Track->AddTPCHit(new TPCLTrackHit(hit));
	    status = true;
	  }
	}
	//else std::cout<<"resi "<<resi<<std::endl;
      } //ci
    } //layer
    return status;
  }
} //namespace

namespace tpc{
//_____________________________________________________________________________
template <typename T> Bool_t
FitStep(T* Track,
	const std::vector<TPCClusterContainer>& ClCont,
	std::vector<T*>& TrackContFailed,
	Int_t MinNumOfHits)
{

  Bool_t status = true;
  if(Track->DoFit(MinNumOfHits)){ // MinNumOfHits cut is not applied in inital tracking
    Track->SetClustersHoughFlag(Candidate);

#if DebugDisp
    Track->Print(FUNC_NAME+" Fitting is succeeded");
    //Track->Print(FUNC_NAME+" Fitting is succeeded", true);
#endif
  }
  else{ //Fitting is failed. (reset clusters' Hough flag)
    Track->SetClustersHoughFlag(0);
    Track->SetFlag(0);
    TrackContFailed.push_back(Track);

#if DebugDisp
    Track->Print(FUNC_NAME+" Fitting is failed");
#endif
    status = false;
  }

  return status;
}

//_____________________________________________________________________________
template <typename T> void
FitTrack(T* Track, Int_t Houghflag,
	 const std::vector<TPCClusterContainer>& ClCont,
	 std::vector<T*>& TrackCont,
	 std::vector<T*>& TrackContFailed,
	 Int_t MinNumOfHits)
{

  std::chrono::milliseconds sec;
  auto fit_start = std::chrono::high_resolution_clock::now();

#if DebugDisp
  Track->Print(FUNC_NAME+" An initial track for fitting");
  //Track->Print(FUNC_NAME+" An initial track for fitting", true);
#endif

  Int_t nstep = 0;
  while(true){
    T *ExtendedTrack = new T(Track);
    //Int_t thr_ncl = 0.5*MinNumOfHits;
    Int_t thr_ncl = 0;
    if(nstep==0) thr_ncl = 0; //The initial track can be short.
    else if(nstep > MaxFitSteps){ //Tracking is over
      delete ExtendedTrack;
      if(IsGood(Track, MinNumOfHits)){
#if DebugDisp
	Track->Print(FUNC_NAME+" Tracking is over. track is good for tracking");
	//Track->Print(FUNC_NAME+" Tracking is over. track is good for tracking", true);
#endif
	Track->SetClustersHoughFlag(Houghflag);
	TrackCont.push_back(Track);
      }
      else{
#if DebugDisp
	Track->Print(FUNC_NAME+" Tracking is over. track is not good for tracking");
	//Track->Print(FUNC_NAME+" Tracking is over. track is not good for tracking", true);
#endif
	Track->SetClustersHoughFlag(BadForTracking); //track w/ few clusters
	Track->SetFlag(0);
	TrackContFailed.push_back(Track);
      }
      break;
    }
    else if(!AddClusters(ExtendedTrack, ClCont)){
      //No more clusters for addding, then check for clusters of bad tracks.
      Bool_t add_badhough = AddClusters(ExtendedTrack, ClCont, BadHoughTransform);
      Bool_t add_badtracking = AddClusters(ExtendedTrack, ClCont, BadForTracking);
      if(!add_badhough && !add_badtracking){
	delete ExtendedTrack;
	if(IsGood(Track, MinNumOfHits)){
#if DebugDisp
	  Track->Print(FUNC_NAME+" track is good and no more cluster to add");
	  //Track->Print(FUNC_NAME+" track is good and no more cluster to add", true);
#endif
	  Track->SetClustersHoughFlag(Houghflag);
	  TrackCont.push_back(Track);
	}
	else{
#if DebugDisp
	  Track->Print(FUNC_NAME+" Track is not good and no more cluster to add");
	  //Track->Print(FUNC_NAME+" Track is not good and no more cluster to add", true);
#endif
	  Track->SetClustersHoughFlag(BadForTracking); //track w/ few clusters
	  Track->SetFlag(0);
	  TrackContFailed.push_back(Track);
	}
	break; //No more cluster to add
      }
    } // No more cluster to add

    //After adding more clusters, fitting starts.
    if(FitStep(ExtendedTrack, ClCont, TrackContFailed, thr_ncl)){
      delete Track;
      Track = ExtendedTrack; //Updates the track
#if DebugDisp
      std::cout<<"extension"<<std::endl;
#endif
    }
    else{
      if(nstep!=0 && IsGood(Track, MinNumOfHits)){
#if DebugDisp
	Track->Print(FUNC_NAME+" no more clusters to add");
	//Track->Print(FUNC_NAME+" no more clusters to add", true);
#endif
	Track->SetClustersHoughFlag(Houghflag);
	TrackCont.push_back(Track);
      }
      else{
#if DebugDisp
	std::cout<<"delete"<<std::endl;
#endif
	delete Track;
      }
      break; //Extended track's fitting is failed
    }

#if DebugDisp
    std::cout<<FUNC_NAME+" # of fitting steps for tracking: "<<nstep<<std::endl;
    std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<" #failed track : "<<TrackContFailed.size()<<std::endl;
    std::cout<<std::endl;
#endif

    nstep++;
  } //while

  auto fittingtime = std::chrono::high_resolution_clock::now();
  sec = std::chrono::duration_cast<std::chrono::milliseconds>(fittingtime - fit_start);
  if(Track) Track -> SetFitTime(sec.count());
  if(Track) Track -> SetFitFlag(nstep);

}

//_____________________________________________________________________________
Int_t
LocalTrackSearch(const std::vector<TPCClusterContainer>& ClCont,
		 std::vector<TPCLocalTrack*>& TrackCont,
		 std::vector<TPCLocalTrack*>& TrackContFailed,
		 Bool_t Exclusive,
		 Int_t MinNumOfHits)
{

  static const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

  XZhough_x.clear();
  XZhough_y.clear();
  Yhough_x.clear();
  Yhough_y.clear();

  Bool_t prev_add = true;
  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    if(!prev_add) continue;
    prev_add = false;

#if DebugDisp
    std::cout<<FUNC_NAME+" tracki : "<<tracki<<std::endl;
#endif

    std::chrono::milliseconds sec;
    auto before_hough = std::chrono::high_resolution_clock::now();

    //Line Hough-transform on the XZ plane
    Double_t LinearPar[4]; Int_t MaxBinXZ[2];
    if(!tpc::HoughTransformLineXZ(ClCont, MaxBinXZ, LinearPar, MinNumOfHits)){
#if DebugDisp
      std::cout<<FUNC_NAME+" No more track candiate! tracki : "<<tracki<<std::endl;
#endif
      break;
    }

    //Line Hough-transform on the YZ or YX plane
    Int_t MaxBinY[2];
    if(TMath::Abs(LinearPar[2]) < 1) tpc::HoughTransformLineYZ(ClCont, MaxBinY, LinearPar, MaxHoughWindowY);
    else tpc::HoughTransformLineYX(ClCont, MaxBinY, LinearPar, MaxHoughWindowY);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrack *track = new TPCLocalTrack;
    track->SetParam(LinearPar);

    //If two tracks are merged at the target, separate them and recalculate params.
    Bool_t vtx_flag;
    prev_add = MakeLinearTrack(track, vtx_flag, ClCont, LinearPar, MaxHoughWindowY);
    if(!prev_add) continue;
    if(!track -> IsGoodForTracking() || !vtx_flag){
      track->SetClustersHoughFlag(BadHoughTransform);
      TrackContFailed.push_back(track);
      continue;
    }

    //Check for duplicates
    Bool_t hough_flag = true;
    for(Int_t i=0; i<XZhough_x.size(); ++i){
      Int_t bindiffXZ = TMath::Abs(MaxBinXZ[0] - XZhough_x[i]) + TMath::Abs(MaxBinXZ[1] - XZhough_y[i]);
      Int_t bindiffY = TMath::Abs(MaxBinY[0] - Yhough_x[i]) + TMath::Abs(MaxBinY[1] - Yhough_y[i]);
      if(bindiffXZ<=1 && bindiffY<=1){
	hough_flag = false;
#if DebugDisp
	std::cout<<"Previous hough bin on the XZ plane "<<i<<"th x: "
		 <<XZhough_x[i]<<", y: "<<XZhough_y[i]<<" on the vertical plane x: "
		 <<Yhough_x[i]<<", y: "<<Yhough_y[i]<<std::endl;
	std::cout<<"Current hough bin on the XZ plane "<<i<<"th x: "
		 <<MaxBinXZ[0]<<", y: "<<MaxBinXZ[1]<<" on the vertical plane x: "
		 <<MaxBinY[0]<<", y: "<<MaxBinY[1]<<std::endl;
#endif
      }
    }
    XZhough_x.push_back(MaxBinXZ[0]);
    XZhough_y.push_back(MaxBinXZ[1]);
    Yhough_x.push_back(MaxBinY[0]);
    Yhough_y.push_back(MaxBinY[1]);

    if(!hough_flag){
#if DebugDisp
      std::cout<<FUNC_NAME+" The same track is found by Hough-Transform : tracki : "<<tracki<<" hough cont size : "<<XZhough_x.size()<<std::endl;
#endif
      track->SetFitFlag(0);
      track->SetClustersHoughFlag(BadHoughTransform);
      TrackContFailed.push_back(track);
      continue;
    }

    auto after_hough = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_hough - before_hough);
    track->SetSearchTime(sec.count());

    //Track fitting processes
    FitTrack(track, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  }// tracki
  ResetHoughFlag(ClCont);
  ResetHoughFlag(ClCont, BadHoughTransform);

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

  CalcTracks(TrackCont);
  CalcTracks(TrackContFailed);
  //MarkingAccidentalTracks(TrackCont);
  if(Exclusive) ExclusiveTracking(TrackCont);
  return TrackCont.size();
}

//_____________________________________________________________________________
Bool_t
MakeLinearTrack(TPCLocalTrack *Track, Bool_t &VtxFlag,
		const std::vector<TPCClusterContainer>& ClCont,
		Double_t *LinearPar, Double_t MaxHoughWindowY){

  Bool_t status = false;

  //Check Hough-distance and add hits
  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
      auto cl = ClCont[layer][ci];
      TPCHit* hit = cl->GetMeanHit();
      if(!hit) continue;
      if(hit->GetHoughFlag()>0) continue;
      TVector3 pos = cl->GetPosition();
      Double_t distXZ = TMath::Abs(LinearPar[2]*(pos.Z() - tpc::ZTarget) - pos.X() +
				   LinearPar[0])/TMath::Sqrt(TMath::Sq(LinearPar[2])+1.);
      Double_t distYZ = TMath::Abs(LinearPar[3]*(pos.Z() - tpc::ZTarget) - pos.Y() +
				   LinearPar[1])/TMath::Sqrt(TMath::Sq(LinearPar[3])+1.);
      if(distXZ < MaxHoughWindowY && distYZ < MaxHoughWindowY){
	hit->SetHoughDist(distXZ);
	hit->SetHoughDistY(distYZ);
	Track->AddTPCHit(new TPCLTrackHit(hit));
	status = true;
      }
    } //ci
  } //layer

  if(status){
    //Vtx in the target, need to check whether two tracks are merged or not
    Track -> SetClustersHoughFlag(Candidate);
    VtxFlag = (Track -> SeparateClustersWithGap() || Track -> SeparateTracksAtTarget());
    //VtxFlag = Track -> SeparateTracksAtTarget();
    if(VtxFlag && AddClusters(Track, ClCont)) Track -> SetClustersHoughFlag(Candidate);
  }

#if DebugDisp
  if(status) Track->Print(FUNC_NAME+" Initial track after track finding");
  //if(status) Track->Print(FUNC_NAME+" Initial track after track finding", true);
#endif

  if(!status) delete Track;
  return status;
}

//_____________________________________________________________________________
Bool_t
MakeHelixTrack(TPCLocalTrackHelix *Track, Bool_t &VtxFlag,
	       const std::vector<TPCClusterContainer>& ClCont,
	       Double_t *HelixPar, Double_t MaxHoughWindow,
	       Double_t MaxHoughWindowY)
{

  Bool_t status = false;

  //Check Hough-distance and add hits
  Int_t id = 0;
  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
      auto cl = ClCont[layer][ci];
      TPCHit* hit = cl->GetMeanHit();
      if(!hit) continue;
      if(hit->GetHoughFlag()>0) continue;
      TVector3 pos = cl->GetPosition();
      Double_t tmpx = -pos.x();
      Double_t tmpy = pos.z() - tpc::ZTarget;
      Double_t tmpz = pos.y();
      Double_t r_cal = TMath::Sqrt(pow(tmpx - HelixPar[0], 2) + pow(tmpy - HelixPar[1], 2));
      Double_t dist = TMath::Abs(r_cal - HelixPar[3]);
      if(dist < MaxHoughWindow){

	//for the inital track, scanning theta within range (-pi, pi)
	Double_t tmpt = TMath::ATan2(tmpy - HelixPar[1], tmpx - HelixPar[0]);
	Double_t tmp_xval = HelixPar[3]*tmpt;
	Double_t distY = TMath::Abs(HelixPar[4]*tmp_xval - tmpz + HelixPar[2])/TMath::Sqrt(pow(HelixPar[4], 2) + 1.);
	Double_t distY_2pi = TMath::Abs(HelixPar[4]*(tmp_xval + 2.*TMath::Pi()*HelixPar[3]) - tmpz + HelixPar[2])/TMath::Sqrt(pow(HelixPar[4], 2) + 1.);
	Double_t distY_m2pi = TMath::Abs(HelixPar[4]*(tmp_xval - 2.*TMath::Pi()*HelixPar[3]) - tmpz + HelixPar[2])/TMath::Sqrt(pow(HelixPar[4], 2) + 1.);
	if(distY < MaxHoughWindowY){
	  //Add hit into the track
	  hit->SetHoughDist(dist);
	  hit->SetHoughDistY(distY);
	  Track->AddTPCHit(new TPCLTrackHit(hit));

	  id++;
	  status = true;
	} //distY
	else if(distY_2pi < MaxHoughWindowY){
	  //Add hit into the track
	  hit->SetHoughDist(dist);
	  hit->SetHoughDistY(distY_2pi);
	  Track->AddTPCHit(new TPCLTrackHit(hit));

	  id++;
	  status = true;
	}
	else if(distY_m2pi < MaxHoughWindowY){
	  //Add hit into the track
	  hit->SetHoughDist(dist);
	  hit->SetHoughDistY(distY_m2pi);
	  Track->AddTPCHit(new TPCLTrackHit(hit));

	  id++;
	  status = true;
	}
      } //dist
    } //ci
  } //layer

  if(status){
    Track -> SetClustersHoughFlag(Candidate);
    Track -> CalcHelixTheta();

    //If track has very large gap, then spilt it.
    if(Track -> SeparateClustersWithGap()) Track -> CalcHelixTheta();

    //If the vtx is in the target, need to check whether two tracks are fragmented tracks or independent tracks.
    VtxFlag = Track -> SeparateTracksAtTarget();
    if(VtxFlag && AddClusters(Track, ClCont)){
      Double_t par[5];
      Track -> GetParam(par);
      VtxFlag = Track -> DoPreFit(par);
      Track -> SetClustersHoughFlag(Candidate);
    }
  }

#if DebugDisp
  if(status) Track->Print(FUNC_NAME+" Initial track after track finding");
  //if(status) Track->Print(FUNC_NAME+" Initial track after track finding", true);
#endif

  if(!status) delete Track;
  return status;
}

//_____________________________________________________________________________
void
HelixTrackSearch(Int_t Trackflag, Int_t Houghflag,
		 const std::vector<TPCClusterContainer>& ClCont,
		 std::vector<TPCLocalTrackHelix*>& TrackCont,
		 std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		 Int_t MinNumOfHits)
{

  // HoughTransform binning
  static const auto MaxHoughWindow = gUser.GetParameter("MaxHoughWindow");
  static const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

  Bool_t prev_add = true;
  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    if(!prev_add) continue;
    prev_add = false;

#if DebugDisp
    std::cout<<FUNC_NAME+" tracki : "<<tracki<<std::endl;
#endif

    std::chrono::milliseconds sec;
    auto before_hough = std::chrono::high_resolution_clock::now();

    //Circle Hough-transform
    Double_t HelixPar[5]; Int_t MaxBinXZ[3];
    if(!tpc::HoughTransformCircleXZ(ClCont, MaxBinXZ, HelixPar, MinNumOfHits)){
#if DebugDisp
      std::cout<<FUNC_NAME+" No more circle candiate! tracki : "<<tracki<<std::endl;
#endif
      break;
    }

    //Linear Hough-transform
    Int_t MaxBinY[2];
    tpc::HoughTransformLineYTheta(ClCont, MaxBinY, HelixPar, MaxHoughWindow);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    track->SetParam(HelixPar);
    track->SetFlag(Trackflag);

    //If two tracks are merged at the target, separate them and recalculate params.
    Bool_t vtx_flag;
    prev_add = MakeHelixTrack(track, vtx_flag, ClCont, HelixPar, MaxHoughWindow, MaxHoughWindowY);
    if(!prev_add) continue;
    if(!track -> IsGoodForTracking() || !vtx_flag){
      track->SetClustersHoughFlag(BadHoughTransform);
      TrackContFailed.push_back(track);
      continue;
    }

    //Check for duplicates
    Bool_t hough_flag = true;
    for(Int_t i=0; i<XZhough_x.size(); ++i){
      Int_t bindiffXZ = TMath::Abs(MaxBinXZ[0] - XZhough_x[i]) + TMath::Abs(MaxBinXZ[1] - XZhough_y[i]) + TMath::Abs(MaxBinXZ[2] - XZhough_z[i]);
      Int_t bindiffY = TMath::Abs(MaxBinY[0] - Yhough_x[i]) + TMath::Abs(MaxBinY[1] - Yhough_y[i]);
      if(bindiffXZ<=1 && bindiffY<=1){
	hough_flag = false;
#if DebugDisp
	std::cout<<"Previous hough bin on the XZ plane "<<i<<"th x: "
		 <<XZhough_x[i]<<", y: "<<XZhough_y[i]<<", z: "<<XZhough_z[i]<<" on the vertical plane x: "
		 <<Yhough_x[i]<<", y: "<<Yhough_y[i]<<std::endl;
	std::cout<<"Current hough bin on the XZ plane "<<i<<"th x: "
		 <<MaxBinXZ[0]<<", y: "<<MaxBinXZ[1]<<", z: "<<MaxBinXZ[2]<<" on the vertical plane x: "
		 <<MaxBinY[0]<<", y: "<<MaxBinY[1]<<std::endl;
#endif
      }
    }
    XZhough_x.push_back(MaxBinXZ[0]);
    XZhough_y.push_back(MaxBinXZ[1]);
    XZhough_z.push_back(MaxBinXZ[2]);
    Yhough_x.push_back(MaxBinY[0]);
    Yhough_y.push_back(MaxBinY[1]);

    if(!hough_flag){
#if DebugDisp
      std::cout<<FUNC_NAME+" The same track is found by Hough-Transform : tracki : "<<tracki<<" hough cont size : "<<XZhough_x.size()<<std::endl;
#endif
      track->SetFitFlag(0);
      track->SetClustersHoughFlag(BadHoughTransform);
      TrackContFailed.push_back(track);
      continue;
    }

    auto after_hough = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_hough - before_hough);
    track->SetSearchTime(sec.count());

    //Track fitting processes
    FitTrack(track, Houghflag, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  }//tracki
  ResetHoughFlag(ClCont);
  ResetHoughFlag(ClCont, BadHoughTransform);

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

}

//_____________________________________________________________________________
void
KuramaTrackSearch(std::vector<std::vector<TVector3>> VPs,
		  std::vector<Double_t> KuramaCharge,
		  const std::vector<TPCClusterContainer>& ClCont,
		  std::vector<TPCLocalTrackHelix*>& TrackCont,
		  std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		  std::vector<TPCLocalTrackHelix*>& TrackContVP,
		  std::vector<TPCVertex*>& VertexCont,
		  Bool_t Exclusive,
		  Int_t MinNumOfHits)
{
  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);

  if(VPs.size()==0) return;
  for(Int_t nt=0; nt<VPs.size(); nt++){
#if DebugDisp
    std::cout<<FUNC_NAME+" Kurama Track "<<nt<<"/"<<VPs.size()<<std::endl;
#endif

    //Reconstruct the kurama track in the HypTPC
    TPCLocalTrackHelix *trackref = new TPCLocalTrackHelix();
    for(Int_t i=0; i<VPs[nt].size(); i++){
      TVector3 pos = VPs[nt][i];
      trackref->AddVPHit(pos);
#if DebugDisp
      std::cout<<FUNC_NAME+" Kurama VP "<<i<<"th pos : "<<VPs[nt][i]<<std::endl;
#endif
    } //i
    trackref->SetIsKurama();
    trackref->SetTrackID(nt);
    if(!trackref->DoVPFit()){
      delete trackref;
      continue;
    }

    Double_t RKHelixParam[5];
    trackref->GetParam(RKHelixParam);
    TrackContVP.push_back(trackref);

#if DebugDisp
    std::cout<<FUNC_NAME+" Kurama VP Helix cx : "<<trackref->Getcx()<<" cy : "<<trackref->Getcy()<<" z0 : "<<trackref->Getz0()<<" r : "<<trackref->Getr()<<" dz : "<<trackref->Getdz()<<std::endl;
    std::cout<<FUNC_NAME+" Kurama VP Helix p : "<<trackref->Getr()*0.299792458<<std::endl;
#endif

    Int_t Trackflag = 1*0 + 2*0 + 4*1 + 8*0; // isBeam, isK18, isKurama, isAccidental
    Int_t id = 0;
    std::vector<Int_t> kurama_candidates;
    for(auto& track: TrackCont){
      if(track){
	Int_t ncl_downstream_tgt = 0; //#cluster after the target.

	//Check whether all track custers within the window along the Kurama track
	if(BeamThroughTPC ||
	   (track -> GetIsK18()!=1 && KuramaCharge[nt]*track->GetCharge()>0)){
	  Int_t nh = track->GetNHit();
	  for(Int_t ih=0; ih<nh; ++ih){
	    TPCLTrackHit *hit = track -> GetHit( ih );
	    if( !hit ) continue;
	    const TVector3& hitpos = hit->GetLocalHitPos();
	    if(hitpos.z()>tpc::ZTarget && !trackref->ResidualCheck(hitpos, KuramaXZWindow, KuramaYWindow)) break;
	    if(hitpos.Z()>tpc::ZTarget) ncl_downstream_tgt++;
	    //If all clusters are in the window, mark the track.
	    if(ih==nh-1 && ncl_downstream_tgt>0){
	      track->AddTrackIDCandidate(nt);
	      track->SetFlag(Trackflag);
	      kurama_candidates.push_back(id);
	    }
	  }
	}
	id++;
      }
    }

    //If there is two candidates, checking whether they are a single track or not.
    if(kurama_candidates.size()==2){
      Int_t trackid1 = kurama_candidates.at(0);
      Int_t trackid2 = kurama_candidates.at(1);

      Bool_t order = TrackCont[trackid1] -> GetNHit() >= TrackCont[trackid2] -> GetNHit() ?  true : false;
      if(!order){
	trackid1 = kurama_candidates.at(1);
	trackid2 = kurama_candidates.at(0);
      }

      TPCLocalTrackHelix *track1 = TrackCont[trackid1];
      TPCLocalTrackHelix *track2 = TrackCont[trackid2];
      TPCLocalTrackHelix *MergedTrack = new TPCLocalTrackHelix(track1);
      for(Int_t hit=0;hit<track2 -> GetNHit();hit++)
	MergedTrack -> AddTPCHit(new TPCLTrackHit(track2 -> GetHitInOrder(hit) -> GetHit()));
      //Check whether tracks belong to the same track or not
      if(!MergedTrack -> TestMergedTrack()){
	delete MergedTrack;
	continue;
      }

      //Two candidates are the same track.
      //Calulate the track and erase previous two fragmented tracks.
      std::vector<Int_t> candidates;
      candidates.push_back(trackid1);
      candidates.push_back(trackid2);

      Int_t prev_size = TrackCont.size();
      FitTrack(MergedTrack, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);

      Bool_t fitstatus = true;
      Int_t post_size = TrackCont.size();
      if(prev_size+1 == post_size){
	MergedTrack = TrackCont[post_size-1];
	MergedTrack -> Calculate();
	MergedTrack -> AddTrackIDCandidate(nt);
	MergedTrack -> SetFlag(Trackflag);
	if(Exclusive){
	  MergedTrack -> DoFitExclusive();
	  MergedTrack -> CalculateExclusive();
	}
      }
      else fitstatus = false;

      if(candidates.size()>0){
	std::sort(candidates.begin(), candidates.end());
	if(fitstatus){
	  for(Int_t i=0; i<candidates.size(); ++i){
	    Int_t trackID = candidates[i] - i;
	    TrackCont.erase(TrackCont.begin() + trackID);
	  }
	  delete track1;
	  delete track2;

	  //Vertex finding again with new tracks
	  del::ClearContainer(VertexCont);
	  VertexSearch(TrackCont, VertexCont);
	}
	else{
	  for(Int_t i=0; i<candidates.size(); ++i){
	    Int_t trackID = candidates[i] - i;
	    TrackCont[trackID] -> SetClustersHoughFlag(GoodForTracking);
	  }
	}
      }
    } //merging two tracks
  } //nt
}

//_____________________________________________________________________________
void
K18TrackSearch(std::vector<std::vector<TVector3>> VPs,
	       const std::vector<TPCClusterContainer>& ClCont,
	       std::vector<TPCLocalTrackHelix*>& TrackCont,
	       std::vector<TPCLocalTrackHelix*>& TrackContVP,
	       Int_t MinNumOfHits /*=2*/)
{

  if(VPs.size()==0) return;
  for(Int_t nt=0; nt<VPs.size(); nt++){
#if DebugDisp
    std::cout<<FUNC_NAME+" K18 Track "<<nt<<"/"<<VPs.size()<<std::endl;
#endif

    std::chrono::milliseconds sec;
    auto before_track_search = std::chrono::high_resolution_clock::now();

    TPCLocalTrackHelix *trackref = new TPCLocalTrackHelix();
    for(Int_t i=0; i<VPs[nt].size(); i++){
      TVector3 pos = VPs[nt][i];
      trackref->AddVPHit(pos);
#if DebugDisp
      std::cout<<FUNC_NAME+" K18 RK "<<i<<"th pos : "<<VPs[nt][i]<<std::endl;
#endif
    } //i
    trackref->SetIsK18();
    trackref->SetIsBeam();
    trackref->DoVPFit();

#if DebugDisp
    std::cout<<FUNC_NAME+" K18 VP Helix cx : "<<trackref->Getcx()<<" cy : "<<trackref->Getcy()<<" z0 : "<<trackref->Getz0()<<" r : "<<trackref->Getr()<<" dz : "<<trackref->Getdz()<<std::endl;
    std::cout<<FUNC_NAME+" K18 VP Helix p : "<<trackref->Getr()*0.299792458<<std::endl;
#endif

    Int_t BeforeTGTHits = 0;
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    //set min, max theta for checking a closest distance and residual
    track -> SetMint(trackref->GetMint());
    track -> SetMaxt(trackref->GetMaxt());

    trackref->DoVPFit();
    std::vector<TPCClusterContainer> ClContK18(NumOfLayersTPC);
    for(Int_t layer=0; layer<10; layer++){ //inner layers
      Double_t minresi = 1000.; Int_t id = -1;
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	auto cl = ClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	if(hit->GetHoughFlag()>0) continue;
	TVector3 pos = cl->GetPosition();
	//if(pos.Z()>tpc::ZTarget) continue;
	if(pos.Z()>tpc::ZTarget-10.) continue;
	Double_t resi;
	if(trackref->ResidualCheck(pos, K18XZWindow, K18YWindow, resi)){
	  ClContK18[layer].push_back(cl);
	  if(minresi > resi){
	    minresi = resi;
	    id = ci;
	  }
	}
      } //ci
      if(minresi<100.){
	TPCHit* hit = ClCont[layer][id]->GetMeanHit();
	hit->SetHoughDist(qnan);
	track->AddTPCHit(new TPCLTrackHit(hit));
	BeforeTGTHits++;
      }
    } //layer

    Double_t RKHelixParam[5];
    trackref->GetParam(RKHelixParam);
    trackref->SetNclBeforeTgt(BeforeTGTHits);
    TrackContVP.push_back(trackref);

#if DebugDisp
    std::cout<<FUNC_NAME+" Clusters upstream of the target and within the window : "<<BeforeTGTHits<<std::endl;
    track->Print(FUNC_NAME+ " K18 track candidate");
#endif

    if(BeforeTGTHits<2){
      delete track;
      continue;
    }

    Int_t Trackflag = 1*1 + 2*1 + 4*0 + 8*0; // isBeam, isK18, isKurama, isAccidental
    auto after_track_search = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_track_search - before_track_search);
    track->SetSearchTime(sec.count());
    track->SetFlag(Trackflag);
    track->SetTrackID(nt); //Marking K18 track id

    //Track fitting processes
    auto fit_start = std::chrono::high_resolution_clock::now();
    if(track->DoFit(RKHelixParam, MinNumOfHits)){ //with constraints
      auto first_fit = std::chrono::high_resolution_clock::now();
      sec = std::chrono::duration_cast<std::chrono::milliseconds>(first_fit - fit_start);
      track->SetFitTime(sec.count());
      track->SetFitFlag(1);
      track->SetClustersHoughFlag(Candidate);

#if DebugDisp
      //track->Print(FUNC_NAME+" First fitting is succeeded");
      track->Print(FUNC_NAME+" First fitting is succeeded", true);
#endif

      TPCLocalTrackHelix *ExtendedTrack = new TPCLocalTrackHelix(track);
      //Residual check with other hits
      if(!AddClusters(ExtendedTrack, ClContK18)){ //No more cluster to add
	delete ExtendedTrack;
	track->SetClustersHoughFlag(K18Tracks);
	TrackCont.push_back(track);
#if DebugDisp
	track->Print(FUNC_NAME+" No added cluster");
	//track->Print(FUNC_NAME+" No added cluster", true);
#endif
      }
      else{ //More clusters are added into the track
	auto hitadd = std::chrono::high_resolution_clock::now();
#if DebugDisp
	ExtendedTrack->Print(FUNC_NAME+" After cluster adding");
	//ExtendedTrack->Print(FUNC_NAME+" After cluster adding", true);
#endif
	if(ExtendedTrack->DoFit(RKHelixParam, MinNumOfHits)){ //2nd fitting success
	  delete track;
	  auto second_fit = std::chrono::high_resolution_clock::now();
	  sec = std::chrono::duration_cast<std::chrono::milliseconds>(second_fit - hitadd);
	  ExtendedTrack->SetFitTime(sec.count());
	  ExtendedTrack->SetFitFlag(2);
	  ExtendedTrack->SetClustersHoughFlag(K18Tracks);
	  TrackCont.push_back(ExtendedTrack);

#if DebugDisp
	  ExtendedTrack->Print(FUNC_NAME+" Fitting success after adding hits (residual check)");
	  //ExtendedTrack->Print(FUNC_NAME+" Fitting success after adding hits (residual check)", true);
#endif
	}
	else{ //2nd fitting failure
	  delete ExtendedTrack;
	  //push 1st fitting track
	  auto hitadd_fail = std::chrono::high_resolution_clock::now();
	  sec = std::chrono::duration_cast<std::chrono::milliseconds>(hitadd_fail - hitadd);
	  track->SetFitTime(sec.count());
	  track->SetFitFlag(3);
	  track->SetClustersHoughFlag(K18Tracks);
	  TrackCont.push_back(track);

#if DebugDisp
	  track->Print(FUNC_NAME+" 2nd fitting failure!");
	  //track->Print(FUNC_NAME+" 2nd fitting failure!", true);
#endif
	}
      }
    }
    else{ //First fitting is failed. Reset clusters' Hough flag
#if DebugDisp
      track->Print(FUNC_NAME+" Fitting is failed");
      //track->Print(FUNC_NAME+" Fitting is failed", true);
#endif
      delete track;
    }
    ResetHoughFlag(ClContK18);
  } //nt
}

//_____________________________________________________________________________
Int_t
LocalTrackSearchHelix(const std::vector<TPCClusterContainer>& ClCont,
		      std::vector<TPCLocalTrackHelix*>& TrackCont,
		      std::vector<TPCLocalTrackHelix*>& TrackContInvertedCharge,
		      std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		      std::vector<TPCVertex*>& VertexCont,
		      std::vector<TPCVertex*>& ClusteredVertexCont,
		      Bool_t Exclusive,
		      Int_t MinNumOfHits)
{
  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);

  //Scattered helix track searching
  HighMomHelixTrackSearch(ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  HelixTrackSearch(0, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);

#if RemainingClustersTest
  ResetHoughFlag(ClCont, BadForTracking);
  HelixTrackSearch(0, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
#endif
  CalcTracks(TrackCont); //before the VertexSearch() calculation should proceed.

  if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);

  //Vertex finding with tracks in the TrackCont.
  VertexSearch(TrackCont, VertexCont);

#if FragmentedTrackTest
  //Merged fragmented tracks
  RestoreFragmentedTracks(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif

#if ReassignClusterTest
  ReassignClustersNearTheTarget(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif

  FindAccidentalCoincidenceTracks(TrackCont, VertexCont, ClusteredVertexCont);

#if ReassignClusterTest
  ReassignClustersVertex(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif

#if RefitXiTrack
  ReassignClustersXiTrack(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif

  TestingCharge(TrackCont, TrackContInvertedCharge, VertexCont, Exclusive);

  RestoreFragmentedAccidentalTracks(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

  CalcTracks(TrackContFailed);
  if(Exclusive) ExclusiveTracking(TrackCont);
  return TrackCont.size();
}

//_____________________________________________________________________________
Int_t
LocalTrackSearchHelix(std::vector<std::vector<TVector3>> K18VPs,
		      std::vector<std::vector<TVector3>> KuramaVPs,
		      std::vector<Double_t> KuramaCharge,
		      const std::vector<TPCClusterContainer>& ClCont,
		      std::vector<TPCLocalTrackHelix*>& TrackCont,
		      std::vector<TPCLocalTrackHelix*>& TrackContInvertedCharge,
		      std::vector<TPCLocalTrackHelix*>& TrackContVP,
		      std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		      std::vector<TPCVertex*>& VertexCont,
		      std::vector<TPCVertex*>& ClusteredVertexCont,
		      Bool_t Exclusive,
		      Int_t MinNumOfHits)
{
  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);

  XZhough_x.clear();
  XZhough_y.clear();
  XZhough_z.clear();
  Yhough_x.clear();
  Yhough_y.clear();

  //Track finding and fitting
  //for K1.8 track searching
  if(!BeamThroughTPC) K18TrackSearch(K18VPs, ClCont, TrackCont, TrackContVP); //default NimNumOfHits = 2
  //Scattered helix track searching
  HighMomHelixTrackSearch(ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  HelixTrackSearch(0, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);

#if RemainingClustersTest
  ResetHoughFlag(ClCont, BadForTracking);
  HelixTrackSearch(0, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
#endif
  CalcTracks(TrackCont); //before the VertexSearch() calculation should proceed.

  if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);

  //Vertex finding with tracks in the TrackCont.
  VertexSearch(TrackCont, VertexCont);

#if FragmentedTrackTest
  //Merged fragmented tracks
  RestoreFragmentedTracks(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif

#if ReassignClusterTest
  ReassignClustersNearTheTarget(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif

  //Search Kurama track candidates from the TrackCont.
  KuramaTrackSearch(KuramaVPs, KuramaCharge, ClCont, TrackCont, TrackContFailed, TrackContVP, VertexCont, Exclusive, MinNumOfHits);

  RestoreFragmentedAccidentalTracks(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);

  FindAccidentalCoincidenceTracks(TrackCont, VertexCont, ClusteredVertexCont);

#if ReassignClusterTest
  ReassignClustersVertex(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif

#if RefitXiTrack
  ReassignClustersXiTrack(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif

  TestingCharge(TrackCont, TrackContInvertedCharge, VertexCont, Exclusive);

  CalcTracks(TrackContVP);
  CalcTracks(TrackContFailed); //Tracking failed cases
  if(Exclusive) ExclusiveTracking(TrackCont);
  return TrackCont.size();
}


//_____________________________________________________________________________
void
HoughTransformTest(const std::vector<TPCClusterContainer>& ClCont,
		   std::vector<TPCLocalTrack*>& TrackCont,
		   Int_t MinNumOfHits /*=8*/)
{
  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);
  static const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

  XZhough_x.clear();
  XZhough_y.clear();
  Yhough_x.clear();
  Yhough_y.clear();

  Bool_t prev_add = true;
  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    if(!prev_add) continue;
    prev_add = false;

#if DebugDisp
    std::cout<<FUNC_NAME+" tracki : "<<tracki<<std::endl;
#endif

    std::chrono::milliseconds sec;
    auto before_hough = std::chrono::high_resolution_clock::now();

    //Line Hough-transform on the XZ plane
    Double_t LinearPar[4]; Int_t MaxBinXZ[2];
    if(!tpc::HoughTransformLineXZ(ClCont, MaxBinXZ, LinearPar, MinNumOfHits)){
#if DebugDisp
      std::cout<<FUNC_NAME+" No more track candiate! tracki : "<<tracki<<std::endl;
#endif
      break;
    }

    Int_t MaxBinY[2];
    //Line Hough-transform on the YZ or YX plane
    if(TMath::Abs(LinearPar[2]) < 1) tpc::HoughTransformLineYZ(ClCont, MaxBinY, LinearPar, MaxHoughWindowY);
    else tpc::HoughTransformLineYX(ClCont, MaxBinY, LinearPar, MaxHoughWindowY);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrack *track = new TPCLocalTrack;
    track->SetParam(LinearPar);

    //If two tracks are merged at the target, separate them and recalculate params.
    Bool_t vtx_flag;
    prev_add = MakeLinearTrack(track, vtx_flag, ClCont, LinearPar, MaxHoughWindowY);
    if(!prev_add) continue;
    if(!track -> IsGoodForTracking() || !vtx_flag){
      track->SetClustersHoughFlag(BadHoughTransform);
      delete track;
      continue;
    }

    //Check for duplicates
    Bool_t hough_flag = true;
    for(Int_t i=0; i<XZhough_x.size(); ++i){
      Int_t bindiffXZ = TMath::Abs(MaxBinXZ[0] - XZhough_x[i]) + TMath::Abs(MaxBinXZ[1] - XZhough_y[i]);
      Int_t bindiffY = TMath::Abs(MaxBinY[0] - Yhough_x[i]) + TMath::Abs(MaxBinY[1] - Yhough_y[i]);
      if(bindiffXZ<=1 && bindiffY<=1){
	hough_flag = false;
#if DebugDisp
	std::cout<<"Previous hough bin on the XZ plane "<<i<<"th x: "
		 <<XZhough_x[i]<<", y: "<<XZhough_y[i]<<" on the vertical plane x: "
		 <<Yhough_x[i]<<", y: "<<Yhough_y[i]<<std::endl;
	std::cout<<"Current hough bin on the XZ plane "<<i<<"th x: "
		 <<MaxBinXZ[0]<<", y: "<<MaxBinXZ[1]<<" on the vertical plane x: "
		 <<MaxBinY[0]<<", y: "<<MaxBinY[1]<<std::endl;
#endif
      }
    }
    XZhough_x.push_back(MaxBinXZ[0]);
    XZhough_y.push_back(MaxBinXZ[1]);
    Yhough_x.push_back(MaxBinY[0]);
    Yhough_y.push_back(MaxBinY[1]);

    if(!hough_flag){
#if DebugDisp
      std::cout<<FUNC_NAME+" The same track is found by Hough-Transform : tracki : "<<tracki<<" hough cont size : "<<XZhough_x.size()<<std::endl;
#endif
      track->SetClustersHoughFlag(BadHoughTransform);
      delete track;
      continue;
    }

    auto after_hough = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_hough - before_hough);
    track->SetClustersHoughFlag(tracki+1);
    track->SetSearchTime(sec.count());
    TrackCont.push_back(track);
  }// tracki

  CalcTracks(TrackCont);
  //if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);
}

//_____________________________________________________________________________
void
HoughTransformTestHelix(const std::vector<TPCClusterContainer>& ClCont,
			std::vector<TPCLocalTrackHelix*>& TrackCont,
			Int_t MinNumOfHits /*=8*/)
{
  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);

  // HoughTransform binning
  static const auto MaxHoughWindow = gUser.GetParameter("MaxHoughWindow");
  static const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

  XZhough_x.clear();
  XZhough_y.clear();
  XZhough_z.clear();
  Yhough_x.clear();
  Yhough_y.clear();

  Bool_t prev_add = true;
  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    if(!prev_add) continue;
    prev_add = false;

#if DebugDisp
    std::cout<<FUNC_NAME+" tracki : "<<tracki<<std::endl;
#endif

    //Circle Hough-transform
    std::chrono::milliseconds sec;
    auto before_hough = std::chrono::high_resolution_clock::now();
    Double_t HelixPar[5]; Int_t MaxBinXZ[3];
    if(!tpc::HoughTransformCircleXZ(ClCont, MaxBinXZ, HelixPar, MinNumOfHits)){
#if DebugDisp
      std::cout<<FUNC_NAME+" No more circle candiate! tracki : "<<tracki<<std::endl;
#endif
      break;
    }

    //Linear Hough-transform
    Int_t MaxBinY[2];
    tpc::HoughTransformLineYTheta(ClCont, MaxBinY, HelixPar, MaxHoughWindow);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    track->SetParam(HelixPar);

    //If two tracks are merged at the target, separate them and recalculate params.
    Bool_t vtx_flag;
    prev_add = MakeHelixTrack(track, vtx_flag, ClCont, HelixPar, MaxHoughWindow, MaxHoughWindowY);
    if(!prev_add) continue;
    if(!track -> IsGoodForTracking() || !vtx_flag){
      track->SetClustersHoughFlag(BadHoughTransform);
      delete track;
      continue;
    }

    //Check for duplicates
    Bool_t hough_flag = true;
    for(Int_t i=0; i<XZhough_x.size(); ++i){
      Int_t bindiffXZ = TMath::Abs(MaxBinXZ[0] - XZhough_x[i]) + TMath::Abs(MaxBinXZ[1] - XZhough_y[i]) + TMath::Abs(MaxBinXZ[2] - XZhough_z[i]);
      Int_t bindiffY = TMath::Abs(MaxBinY[0] - Yhough_x[i]) + TMath::Abs(MaxBinY[1] - Yhough_y[i]);
      if(bindiffXZ<=1 && bindiffY<=1){
	hough_flag = false;
#if DebugDisp
	std::cout<<"Previous hough bin on the XZ plane "<<i<<"th x: "
		 <<XZhough_x[i]<<", y: "<<XZhough_y[i]<<", z: "<<XZhough_z[i]<<" on the vertical plane x: "
		 <<Yhough_x[i]<<", y: "<<Yhough_y[i]<<std::endl;
	std::cout<<"Current hough bin on the XZ plane "<<i<<"th x: "
		 <<MaxBinXZ[0]<<", y: "<<MaxBinXZ[1]<<", z: "<<MaxBinXZ[2]<<" on the vertical plane x: "
		 <<MaxBinY[0]<<", y: "<<MaxBinY[1]<<std::endl;
#endif
      }
    }
    XZhough_x.push_back(MaxBinXZ[0]);
    XZhough_y.push_back(MaxBinXZ[1]);
    XZhough_z.push_back(MaxBinXZ[1]);
    Yhough_x.push_back(MaxBinY[0]);
    Yhough_y.push_back(MaxBinY[1]);

    if(!hough_flag){
#if DebugDisp
      std::cout<<FUNC_NAME+" The same track is found by Hough-Transform : tracki : "<<tracki<<" hough cont size : "<<XZhough_x.size()<<std::endl;
#endif
      track->SetClustersHoughFlag(BadHoughTransform);
      delete track;
      continue;
    }

    auto after_hough = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_hough - before_hough);

    track->SetClustersHoughFlag(tracki + 1);
    track->SetSearchTime(sec.count());
    TrackCont.push_back(track);
  }//tracki

  CalcTracks(TrackCont);
  //if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);
}

//_____________________________________________________________________________
void
HighMomHelixTrackSearch(const std::vector<TPCClusterContainer>& ClCont,
			std::vector<TPCLocalTrackHelix*>& TrackCont,
			std::vector<TPCLocalTrackHelix*>& TrackContFailed,
			Int_t MinNumOfHits)
{

  static const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

  std::vector<Double_t> tempXZhough_x;
  std::vector<Double_t> tempXZhough_y;
  std::vector<Double_t> tempYhough_x;
  std::vector<Double_t> tempYhough_y;

  Bool_t prev_add = true;
  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    if(!prev_add) continue;
    prev_add = false;

#if DebugDisp
    std::cout<<FUNC_NAME+" tracki : "<<tracki<<std::endl;
#endif

    std::chrono::milliseconds sec;
    auto before_hough = std::chrono::high_resolution_clock::now();

    //Line Hough-transform on the XZ plane
    Double_t LinearPar[4]; Int_t MaxBinXZ[2];
    if(!tpc::HoughTransformLineXZ(ClCont, MaxBinXZ, LinearPar, MinNumOfHits)){
#if DebugDisp
      std::cout<<FUNC_NAME+" No more track candiate! tracki : "<<tracki<<std::endl;
#endif
      break;
    }

    //Line Hough-transform on the YZ or YX plane
    Int_t MaxBinY[2];
    if(TMath::Abs(LinearPar[2]) < 1) tpc::HoughTransformLineYZ(ClCont, MaxBinY, LinearPar, MaxHoughWindowY);
    else tpc::HoughTransformLineYX(ClCont, MaxBinY, LinearPar, MaxHoughWindowY);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrack *trackTemp = new TPCLocalTrack;
    trackTemp->SetParam(LinearPar);

    Bool_t vtx_flag;
    prev_add = MakeLinearTrack(trackTemp, vtx_flag, ClCont, LinearPar, MaxHoughWindowY);
    if(!prev_add) continue;
    if(!trackTemp -> IsGoodForTracking() || !vtx_flag){
      trackTemp->SetClustersHoughFlag(BadHoughTransform);
      delete trackTemp;
      continue;
    }

    //Check for duplicates
    Bool_t hough_flag = true;
    for(Int_t i=0; i<tempXZhough_x.size(); ++i){
      Int_t bindiffXZ = TMath::Abs(MaxBinXZ[0] - tempXZhough_x[i]) + TMath::Abs(MaxBinXZ[1] - tempXZhough_y[i]);
      Int_t bindiffY = TMath::Abs(MaxBinY[0] - tempYhough_x[i]) + TMath::Abs(MaxBinY[1] - tempYhough_y[i]);
      if(bindiffXZ<=1 && bindiffY<=1){
	hough_flag = false;
#if DebugDisp
	std::cout<<"Previous hough bin on the XZ plane "<<i<<"th x: "
		 <<tempXZhough_x[i]<<", y: "<<tempXZhough_y[i]<<" on the vertical plane x: "
		 <<tempYhough_x[i]<<", y: "<<tempYhough_y[i]<<std::endl;
	std::cout<<"Current hough bin on the XZ plane "<<i<<"th x: "
		 <<MaxBinXZ[0]<<", y: "<<MaxBinXZ[1]<<" on the vertical plane x: "
		 <<MaxBinY[0]<<", y: "<<MaxBinY[1]<<std::endl;
#endif
      }
    }
    tempXZhough_x.push_back(MaxBinXZ[0]);
    tempXZhough_y.push_back(MaxBinXZ[1]);
    tempYhough_x.push_back(MaxBinY[0]);
    tempYhough_y.push_back(MaxBinY[1]);

    if(!hough_flag){
#if DebugDisp
      std::cout<<FUNC_NAME+" The same track is found by Hough-Transform : tracki : "<<tracki<<" hough cont size : "<<tempXZhough_x.size()<<std::endl;
#endif
      trackTemp->SetClustersHoughFlag(BadHoughTransform);
      delete trackTemp;
      continue;
    }

    //Convert the temporary straight line track into helix track
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix;
    if(!ConvertTrack(trackTemp, track)){
      //track -> Print("converting is failed", true);
#if DebugDisp
      track->Print(FUNC_NAME+ " Track converting is failed");
#endif
      track -> SetClustersHoughFlag(BadForTracking); //track w/ few clusters
      delete track;
      continue;
    }

#if DebugDisp
    track->Print(FUNC_NAME+ " Track converting is succeeded");
#endif

    auto after_hough = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_hough - before_hough);
    track->SetSearchTime(sec.count());

    //Track fitting processes
    FitTrack(track, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  }// tracki
  ResetHoughFlag(ClCont);
  ResetHoughFlag(ClCont, BadHoughTransform);
  ResetHoughFlag(ClCont, BadForTracking);

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

}

//_____________________________________________________________________________
template <typename T>
void
VertexSearch(std::vector<T*>& TrackCont,
	     std::vector<TPCVertex*>& VertexCont)
{
  //pair
  for(Int_t trackid1=0; trackid1<TrackCont.size(); trackid1++){
    T *track1 = TrackCont[trackid1];
    for(Int_t trackid2=trackid1+1; trackid2<TrackCont.size(); trackid2++){
      T *track2 = TrackCont[trackid2];
      TPCVertex *vertex = new TPCVertex(trackid1, trackid2);
      vertex -> Calculate(track1, track2);

      Double_t closedist = vertex -> GetClosestDist();
      if(closedist < VertexDistCut){
	VertexCont.push_back(vertex);

#if DebugDisp
	vertex->Print(FUNC_NAME+" Vertex finding");
	//vertex->Print(FUNC_NAME+" Vertex finding", true);
#endif
      }
      else delete vertex;
    } //trackid2
  } //trackid1
}

//_____________________________________________________________________________
template <typename T> void
RestoreFragmentedTracks(const std::vector<TPCClusterContainer>& ClCont,
			std::vector<T*>& TrackCont,
			std::vector<T*>& TrackContFailed,
			std::vector<TPCVertex*>& VertexCont,
			Bool_t Exclusive,
			Int_t MinNumOfHits)
{

  std::vector<TPCVertex*> candidates_VertexCont;
  std::vector<Int_t> candidates;
  for(auto& vertex: VertexCont){
    //Threshold conditions of a distance and an angle between two tracks.
    TVector3 vtx = vertex -> GetVertex();
    Int_t trackid1 = vertex -> GetTrackId(0);
    Int_t trackid2 = vertex -> GetTrackId(1);
    if(TrackCont[trackid1] -> GetIsK18()==1 ||
       TrackCont[trackid2] -> GetIsK18()==1) continue;

    Bool_t isbeam1 = true;
    for(Int_t ihit=0;ihit<TrackCont[trackid1] -> GetNHit();ihit++){
      TPCHit *hit = TrackCont[trackid1] -> GetHitInOrder(ihit) -> GetHit();
      TVector3 pos = hit -> GetPosition();
      if(pos.z() > tpc::ZTarget || TMath::Abs(pos.x()) > 15 || TMath::Abs(pos.y()) > 10.){
	isbeam1 = false;
	break;
      }
    }
    if(isbeam1) continue; //veto beam particle

    Bool_t isbeam2 = true;
    for(Int_t ihit=0;ihit<TrackCont[trackid2] -> GetNHit();ihit++){
      TPCHit *hit = TrackCont[trackid2] -> GetHitInOrder(ihit) -> GetHit();
      TVector3 pos = hit -> GetPosition();
      if(pos.z() > tpc::ZTarget || TMath::Abs(pos.x()) > 15 || TMath::Abs(pos.y()) > 10.){
	isbeam2 = false;
	break;
      }
    }
    if(isbeam2) continue; //veto beam particle
    /*
      std::cout<<"id "<<vertex -> GetTrackId(0)<<" "<<vertex -> GetTrackId(1)<<std::endl;
      std::cout<<"vtx "<<TMath::Abs(vtx.x())<<" "<<TMath::Abs(vtx.y())<<" "<<TMath::Abs(vtx.z())<<std::endl;
      std::cout<<"angle "<<0.1*TMath::Pi()<<" "<<vertex -> GetOpeningAngle()<<" "<<(1. - 0.1)*TMath::Pi()<<std::endl;
    */

    //case1. accidental beam crossing the target is splitted into two tracks
    //case2. merging fragmentations of commom track. for case2. closest point of two parts are not in the track

    //two parts are close.
    if(vertex -> GetClosestDist() > 20. ||
       (vertex -> GetOpeningAngle() > 0.1*TMath::Pi() &&
	vertex -> GetOpeningAngle() < (1. - 0.1)*TMath::Pi())) continue;

    //closest point is not in the target.
    //without this, two scattered tracks with opposite direction frequently wrongly merged.
    if(TMath::Abs(vtx.x()) < 15. &&
       TMath::Abs(vtx.y()) < 10. &&
       TMath::Abs(vtx.z() - tpc::ZTarget) < 10.) continue;

#if DebugDisp
    vertex -> Print(FUNC_NAME+" Vertex candidate for merging tracks");
    //vertex -> Print(FUNC_NAME+" Vertex candidate for merging tracks", true);
#endif

    //Avoid duplication
    if(find(candidates.begin(), candidates.end(), trackid1) != candidates.end()) continue;
    if(find(candidates.begin(), candidates.end(), trackid2) != candidates.end()) continue;

    //Longer track(track1) is a reference.
    //Add two tracks.
    Bool_t order = TrackCont[trackid1] -> GetNHit() >= TrackCont[trackid2] -> GetNHit() ?  true : false;
    if(!order){
      trackid1 = vertex -> GetTrackId(0);
      trackid2 = vertex -> GetTrackId(1);
    }

    T *track1 = TrackCont[trackid1];
    T *track2 = TrackCont[trackid2];
    T *MergedTrack = new T(track1);
    for(Int_t hit=0;hit<track2 -> GetNHit();hit++)
      MergedTrack -> AddTPCHit(new TPCLTrackHit(track2 -> GetHitInOrder(hit) -> GetHit()));

    //Check whether tracks belong to the same track or not
    if(!MergedTrack -> TestMergedTrack()){
      candidates_VertexCont.push_back(vertex);
      delete MergedTrack;
      continue;
    }
    candidates.push_back(trackid1);
    candidates.push_back(trackid2);

#if DebugDisp
    MergedTrack -> Print(FUNC_NAME+" Before fitting the merged track");
    //MergedTrack -> Print(FUNC_NAME+" Before fitting the merged track", true);
#endif

    Int_t prev_size = TrackCont.size();
    FitTrack(MergedTrack, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);

    Int_t post_size = TrackCont.size();
    if(prev_size+1 == post_size){
      MergedTrack = TrackCont[post_size-1];
      MergedTrack -> Calculate();
      if(Exclusive){
	MergedTrack -> DoFitExclusive();
	MergedTrack -> CalculateExclusive();
      }
      MergedTrack -> CheckIsAccidental();

#if DebugDisp
      MergedTrack -> Print(FUNC_NAME+" After fitting the merged track");
      //MergedTrack -> Print(FUNC_NAME+" After fitting the merged track", true);
#endif
    }
    else{
      candidates_VertexCont.push_back(vertex);
      for(Int_t i=0; i<2; ++i){
	Int_t trackID = candidates[candidates.size() - 1];
	TrackCont[trackID] -> SetClustersHoughFlag(GoodForTracking);
	candidates.erase(candidates.begin() + candidates.size() - 1);
      }
    }
  }

  //(Iterative process) For candidates with fitting failed, trying another way for fitting.
  for(auto& vertex: candidates_VertexCont){
    TVector3 vtx = vertex -> GetVertex();
    Int_t trackid1 = vertex -> GetTrackId(0);
    Int_t trackid2 = vertex -> GetTrackId(1);

#if DebugDisp
    vertex -> Print(FUNC_NAME+" Vertex candidate for merging tracks");
    //vertex -> Print(FUNC_NAME+" Vertex candidate for merging tracks", true);
#endif

    //Avoid duplication
    if(find(candidates.begin(), candidates.end(), trackid1) != candidates.end()) continue;
    if(find(candidates.begin(), candidates.end(), trackid2) != candidates.end()) continue;

    //Longer track(track1) is a reference.
    //Add two tracks.
    Bool_t order = TrackCont[trackid1] -> GetNHit() >= TrackCont[trackid2] -> GetNHit() ?  true : false;
    if(!order){
      trackid1 = vertex -> GetTrackId(0);
      trackid2 = vertex -> GetTrackId(1);
    }

    T *track1 = TrackCont[trackid1];
    T *track2 = TrackCont[trackid2];
    T *MergedTrack = new T(track1);
    Int_t total_nhit = track1 -> GetNHit() + track2 -> GetNHit();
    std::vector<TPCClusterContainer> clusters_fortest(NumOfLayersTPC);
    for(Int_t hitid=0;hitid<track2 -> GetNHit();hitid++){
      TPCHit *hit = track2 -> GetHitInOrder(hitid) -> GetHit();
      hit -> SetHoughFlag(0);
      Int_t layer = hit -> GetLayer();
      clusters_fortest[layer].push_back(hit -> GetParentCluster());
    }

    Int_t prev_size = TrackCont.size();
    FitTrack(MergedTrack, GoodForTracking, clusters_fortest, TrackCont, TrackContFailed, MinNumOfHits);

    Int_t post_size = TrackCont.size();
    if(prev_size+1 == post_size){
      MergedTrack = TrackCont[post_size-1];

#if DebugDisp
      std::cout<<FUNC_NAME+" #of bad clusters : "<<total_nhit - MergedTrack -> GetNHit()<<std::endl;
#endif

      if(total_nhit - MergedTrack -> GetNHit() > 3){
	delete MergedTrack;
	TrackCont.erase(TrackCont.begin() + post_size - 1);
	track1 -> SetClustersHoughFlag(GoodForTracking);
	track2 -> SetClustersHoughFlag(GoodForTracking);
	continue;
      }
      else{
	MergedTrack -> Calculate();
	if(Exclusive){
	  MergedTrack -> DoFitExclusive();
	  MergedTrack -> CalculateExclusive();
	}
	MergedTrack -> CheckIsAccidental();
      }
    }
    else{
      track1 -> SetClustersHoughFlag(GoodForTracking);
      track2 -> SetClustersHoughFlag(GoodForTracking);
      continue;
    }

    candidates.push_back(trackid1);
    candidates.push_back(trackid2);
  }
  candidates_VertexCont.clear();

  if(candidates.size()>0){
    std::sort(candidates.begin(), candidates.end());
    for(Int_t i=0; i<candidates.size(); ++i){
      Int_t trackID = candidates[i] - i;
      T *prevtrack = TrackCont[trackID];
      TrackCont.erase(TrackCont.begin() + trackID);
      delete prevtrack;
    }

    //Vertex finding again with new tracks
    del::ClearContainer(VertexCont);
    VertexSearch(TrackCont, VertexCont);
  }
}

//_____________________________________________________________________________
template <typename T> void
ReassignClustersNearTheTarget(const std::vector<TPCClusterContainer>& ClCont,
			      std::vector<T*>& TrackCont,
			      std::vector<T*>& TrackContFailed,
			      std::vector<TPCVertex*>& VertexCont,
			      Bool_t Exclusive,
			      Int_t MinNumOfHits)
{

  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);
  const Int_t MostInnerlayer_Track = 6; //testing layers from 0 to "MostInnerlayer_Cluster".
  const Int_t MostInnerlayer_Cluster = 4; //testing layers from 0 to "MostInnerlayer_Cluster".

  Bool_t status = false;
  Int_t ntracks = TrackCont.size();
  Int_t k18id = -9999;
  std::vector<Int_t> candidates_trackid;
  std::vector<T*> newtracks_fortest;
  std::vector<TPCHit*> clusters_fortest;
  std::vector<Int_t> clusters_original_trackid;
  for(Int_t trackid=0; trackid<ntracks; trackid++){
    T *track = TrackCont[trackid];
    if(track -> GetIsAccidental()==1) continue;
    if(track -> GetIsK18()==1) k18id = trackid;
    if(track -> GetClosestDist() < 50.){ //tracks from the target
      Int_t n = track->GetNHit();
      TPCLTrackHit *hitp0 = track->GetHitInOrder(0); //most inner cluster
      TPCHit *hit0 = hitp0->GetHit();
      Int_t layer0 = hitp0->GetLayer();
      TVector3 pos0 = hitp0->GetLocalHitPos();
      if(TMath::Abs(pos0.y()) > 50.) continue;
      if(layer0 > MostInnerlayer_Track) continue;
      status = true;
      candidates_trackid.push_back(trackid);

      T *CopiedTrack0 = new T(track);
      newtracks_fortest.push_back(CopiedTrack0);

      if(n <= MinNumOfHits || layer0 > MostInnerlayer_Cluster) continue;
      T *CopiedTrack1 = new T(CopiedTrack0);
      Int_t hitorder0 = CopiedTrack0 -> GetOrder(0);
      CopiedTrack1 -> EraseHit(hitorder0);
      hit0 -> SetHoughFlag(GoodForTracking);
      if(CopiedTrack1->DoFit(MinNumOfHits)){
	clusters_fortest.push_back(hit0);
	clusters_original_trackid.push_back(trackid);

	newtracks_fortest[newtracks_fortest.size() - 1] = CopiedTrack1;
	delete CopiedTrack0;

	TPCLTrackHit *hitp1 = track->GetHitInOrder(1); //second most inner cluster
	TPCHit *hit1 = hitp1->GetHit();
	Int_t layer1 = hitp1->GetLayer();
	TVector3 pos1 = hitp1->GetLocalHitPos();
	if(TMath::Abs(pos1.y()) > 50.) continue;
	if(n-1 <= MinNumOfHits || layer1 > MostInnerlayer_Cluster) continue;

	//Exclude most inner cluster and dofit.
	T *CopiedTrack2 = new T(CopiedTrack1);
	Int_t hitorder1 = CopiedTrack2 -> GetOrder(0);
	CopiedTrack2 -> EraseHit(hitorder1);
	hit1 -> SetHoughFlag(GoodForTracking);
	if(CopiedTrack2->DoFit(MinNumOfHits)){
	  clusters_fortest.push_back(hit1);
	  clusters_original_trackid.push_back(trackid);

	  newtracks_fortest[newtracks_fortest.size() - 1] = CopiedTrack2;
	  delete CopiedTrack1;
	}
	else delete CopiedTrack2;
      }
      else delete CopiedTrack1;
    }
  }
  if(!status) return;

  T *copiedK18track;
  if(k18id!=-9999){
    T *K18track = TrackCont[k18id];
    copiedK18track = new T(K18track);
    candidates_trackid.push_back(k18id);
    newtracks_fortest.push_back(copiedK18track);
  }

  //After excluding clusters near the target, try to add more clusters
  //Because of clusters near the target, frequently tracking quility becomes worse.
  //So without those clusters, it is better to try to add clusters
  Bool_t reassaign = false;
  std::vector<Int_t> reassaign_candidates;

  for(Int_t i=0; i<newtracks_fortest.size(); i++){
    Int_t trackid = candidates_trackid[i];
    if(trackid == k18id) continue;

    T *track = newtracks_fortest[i];
    T *CopiedTrack = new T(track);

    Int_t prev_size = TrackCont.size();
    FitTrack(CopiedTrack, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
    Int_t post_size = TrackCont.size();
    if(prev_size+1 == post_size){ //Fitting is succeeded
      CopiedTrack = TrackCont[post_size-1];
      if(Exclusive){
	CopiedTrack -> DoFitExclusive();
	CopiedTrack -> CalculateExclusive();
      }
      TrackCont.erase(TrackCont.begin() + post_size - 1);

      if(CopiedTrack -> GetNHit() > track -> GetNHit()){
	reassaign = true;
	reassaign_candidates.push_back(trackid);
	newtracks_fortest[i] = CopiedTrack;
	delete track;
      }
      else delete CopiedTrack;
    }
  }

#if DebugDisp
  std::cout<<FUNC_NAME+" Finding the best combination between tracks and clusters"<<std::endl;
#endif

  for(Int_t i=0; i<clusters_fortest.size(); i++){
    TPCHit *hit = clusters_fortest[i]; //cluster for testing
    TVector3 pos = hit -> GetPosition();

    Int_t best = -9999;
    Int_t best_trackid = -1;
    Double_t min_residual = 9999;
    for(Int_t j=0; j<newtracks_fortest.size(); j++){
      Double_t residual = 0;
      if(newtracks_fortest[j] -> IsGoodHitToAdd(hit, residual)){
	if(min_residual > residual){
	  min_residual = residual;
	  best = j;
	  best_trackid = candidates_trackid[j];
	}
      }
    }

    Int_t id = clusters_original_trackid[i]; //cluster's original id
    hit -> SetHoughFlag(GoodForTracking);
    if(id==best_trackid){ //cluster is suitable for the original track
      newtracks_fortest[best] -> AddTPCHit(new TPCLTrackHit(hit));
    }
    else{ //cluster is not suitable for the original track
      reassaign = true;
      if(best_trackid==-1){ //cluster is not suitable for all tracks
	reassaign_candidates.push_back(id);
	hit -> SetHoughFlag(0);
      }
      else{ //cluster is suitable for the other track not the original track
	reassaign_candidates.push_back(id);
	reassaign_candidates.push_back(best_trackid);
	newtracks_fortest[best] -> AddTPCHit(new TPCLTrackHit(hit));
      }
    }
  }
  if(newtracks_fortest.size()!=candidates_trackid.size()) std::cout<<FUNC_NAME+" FATAL Error #of tracks is not matched"<<std::endl;

  std::sort(reassaign_candidates.begin(), reassaign_candidates.end());
  reassaign_candidates.erase(std::unique(reassaign_candidates.begin(), reassaign_candidates.end()), reassaign_candidates.end());

  if(!reassaign){ //No change
    if(reassaign_candidates.size()!=0) std::cout<<FUNC_NAME+" FATAL Error #of tracks is not matched"<<std::endl;
#if DebugDisp
    std::cout<<FUNC_NAME+" No reassigning happened"<<std::endl;
#endif
    for(Int_t i=0; i<newtracks_fortest.size(); i++) delete newtracks_fortest[i];
  }
  else{ //Reassaigning happens
#if DebugDisp
    std::cout<<FUNC_NAME+" Candidates for reassigning "<<reassaign_candidates.size()<<std::endl;
#endif

    for(Int_t i=0; i<newtracks_fortest.size(); i++){
      Int_t trackid = candidates_trackid[i];
      if(std::find(reassaign_candidates.begin(), reassaign_candidates.end(), trackid)
	 == reassaign_candidates.end()) delete newtracks_fortest[i];
      else{
	Int_t threshold = MinNumOfHits;
	if(trackid == k18id) threshold = 3;

	if(newtracks_fortest[i] -> DoFit(threshold)){
	  newtracks_fortest[i] -> SetClustersHoughFlag(GoodForTracking);
	  newtracks_fortest[i] -> Calculate();
	  if(Exclusive){
	    newtracks_fortest[i] -> DoFitExclusive();
	    newtracks_fortest[i] -> CalculateExclusive();
	  }

	  //update with a new track and delete the previous track
	  T *prev_track = TrackCont[trackid];
	  TrackCont[trackid] = newtracks_fortest[i];
	  delete prev_track;
	}
	else delete newtracks_fortest[i];
      }
    } //for scattered tracks
  }
  ResetHoughFlag(ClCont);

  if(reassaign){
    //Vertex finding again with new tracks
    del::ClearContainer(VertexCont);
    VertexSearch(TrackCont, VertexCont);
#if FragmentedTrackTest
    //Merged fragmented tracks
    RestoreFragmentedTracks(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif
    if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);
  } //reassaigning process
}

//_____________________________________________________________________________
template <typename T> void
ReassignClustersVertex(const std::vector<TPCClusterContainer>& ClCont,
		       std::vector<T*>& TrackCont,
		       std::vector<T*>& TrackContFailed,
		       std::vector<TPCVertex*>& VertexCont,
		       Bool_t Exclusive,
		       Int_t MinNumOfHits)
{

  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);

  Int_t ntracks = TrackCont.size();

  std::vector<std::pair<Int_t, Int_t>> candidate_pairs;
  Double_t window = 20; //selection cut : cluster <-> track dist < window
  Bool_t flag = true;
  while(flag){
    flag = false;

    for(auto& vertex: VertexCont){
      TVector3 vtx = vertex -> GetVertex();
      Int_t trackid1 = vertex -> GetTrackId(0);
      Int_t trackid2 = vertex -> GetTrackId(1);

      if(std::find(candidate_pairs.begin(), candidate_pairs.end(), std::make_pair(trackid1, trackid2)) != candidate_pairs.end()) continue;

      T *track1 = TrackCont[trackid1];
      T *track2 = TrackCont[trackid2];
      //if(vertex -> GetClosestDist() > 5) continue;
      if(vertex -> GetClosestDist() > 10) continue;
      if(track1 -> GetIsAccidental()==1 || track2 -> GetIsAccidental()==1) continue;
      if(track1 -> GetIsK18()==1 || track2 -> GetIsK18()==1) continue;

      Double_t theta1 = vertex -> GetTrackTheta(0);
      Double_t theta2 = vertex -> GetTrackTheta(1);
      //vertex point exists outside of both tracks, then it is not a candidate.
      if((theta1 < track1 -> GetMint() - 10./track1 -> Getr() ||
	  theta1 > track1 -> GetMaxt() + 10./track1 -> Getr()) ||
	 (theta2 < track2 -> GetMint() - 10./track2 -> Getr() ||
	  theta2 > track2 -> GetMaxt() + 10./track2 -> Getr())) continue;

      Int_t closeid1;
      Double_t minresi1 = 9999;
      std::vector<Int_t> close_ids1;
      TVector3 vtx1 = vertex -> GetTrackPos(0);
      for(Int_t i=0; i<track1 -> GetNHit(); ++i){
	TVector3 dist = track1 -> GetHitInOrder(i) -> GetLocalCalPosHelix() - vtx1; //distance between the vtx to the cluster
	if(dist.Mag() < window){
	  close_ids1.push_back(track1 -> GetOrder(i));
	}
	if(dist.Mag() < minresi1){
	  closeid1 = i;
	  minresi1 = dist.Mag();
	}
      }

      Int_t closeid2;
      Double_t minresi2 = 9999;
      std::vector<Int_t> close_ids2;
      TVector3 vtx2 = vertex -> GetTrackPos(1);
      for(Int_t i=0; i<track2 -> GetNHit(); ++i){
	TVector3 dist = track2 -> GetHitInOrder(i) -> GetLocalCalPosHelix() - vtx2; //distance between the vtx to the cluster
	if(dist.Mag() < window){
	  close_ids2.push_back(track2 -> GetOrder(i));
	}
	if(dist.Mag() < minresi2){
	  closeid2 = i;
	  minresi2 = dist.Mag();
	}
      }
      if(closeid1==0 && closeid2==0 && (close_ids2.size()!=1 && close_ids1.size()!=1)) continue; //vertex at the end of two tracks

      //case1
      std::vector<T*> newtracks_fortest;
      std::vector<Int_t> candidates_trackid;
      std::vector<TPCHit*> clusters_fortest;
      std::vector<Int_t> clusters_original_trackid;
      if(close_ids2.size() > 0 || close_ids1.size() > 0){
	T *CopiedTrack1 = new T(track1);
	CopiedTrack1 -> EraseHits(close_ids1);
	if(!CopiedTrack1->DoFit(MinNumOfHits)){
	  track1 -> SetClustersHoughFlag(GoodForTracking);
	  delete CopiedTrack1;
	  continue;
	}

	T *CopiedTrack2 = new T(track2);
	CopiedTrack2 -> EraseHits(close_ids2);
	if(!CopiedTrack2->DoFit(MinNumOfHits)){
	  track1 -> SetClustersHoughFlag(GoodForTracking);
	  track2 -> SetClustersHoughFlag(GoodForTracking);
	  delete CopiedTrack1;
	  delete CopiedTrack2;
	  continue;
	}

	newtracks_fortest.push_back(CopiedTrack1);
	candidates_trackid.push_back(trackid1);
	for(Int_t i=0; i<close_ids1.size(); ++i){
	  TPCHit *hit = track1 -> GetHit(close_ids1[i]) -> GetHit();
	  clusters_fortest.push_back(hit);
	  clusters_original_trackid.push_back(trackid1);
	}

	newtracks_fortest.push_back(CopiedTrack2);
	candidates_trackid.push_back(trackid2);
	for(Int_t i=0; i<close_ids2.size(); ++i){
	  TPCHit *hit = track2 -> GetHit(close_ids2[i]) -> GetHit();
	  clusters_fortest.push_back(hit);
	  clusters_original_trackid.push_back(trackid2);
	}

	for(Int_t i=0; i<clusters_fortest.size(); i++){
	  TPCHit *hit = clusters_fortest[i]; //cluster for testing
	  hit -> SetHoughFlag(Candidate);
	  TVector3 pos = hit -> GetPosition();

	  Int_t best = -9999;
	  Int_t best_trackid = -1;
	  Double_t min_residual = 9999;
	  for(Int_t j=0; j<newtracks_fortest.size(); j++){
	    Double_t residual = 0;
	    if(newtracks_fortest[j] -> IsGoodHitToAdd(hit, residual)){
	      if(min_residual > residual){
		min_residual = residual;
		best = j;
		best_trackid = candidates_trackid[j];
	      }
	    }
	  }

	  Int_t id = clusters_original_trackid[i]; //cluster's original id
	  if(id==best_trackid){ //cluster is suitable for the original track
	    hit -> SetHoughFlag(GoodForTracking);
	    newtracks_fortest[best] -> AddTPCHit(new TPCLTrackHit(hit));
	  }
	  else{ //cluster is not suitable for the original track
	    flag = true;
	    if(best_trackid==-1){ //cluster is not suitable for all tracks
	      hit -> SetHoughFlag(0);
	    }
	    else{ //cluster is suitable for the other track not the original track
	      newtracks_fortest[best] -> AddTPCHit(new TPCLTrackHit(hit));
	    }
	  }
	} //for(Int_t i=0; i<clusters_fortest.size(); i++){

	if(!flag){ //no change -> checking next track pair
	  for(Int_t i=0; i<newtracks_fortest.size(); i++) delete newtracks_fortest[i];
	  track1 -> SetClustersHoughFlag(GoodForTracking);
	  track2 -> SetClustersHoughFlag(GoodForTracking);
	}
	else{ //reassaigning happens
	  for(Int_t i=0; i<newtracks_fortest.size(); i++){
	    Int_t trackid = candidates_trackid[i];
	    T *prev_track = TrackCont[trackid];

	    if(newtracks_fortest[i] -> DoFit(MinNumOfHits)){
	      newtracks_fortest[i] -> SetClustersHoughFlag(GoodForTracking);
	      newtracks_fortest[i] -> Calculate();
	      if(Exclusive){
		newtracks_fortest[i] -> DoFitExclusive();
		newtracks_fortest[i] -> CalculateExclusive();
	      }

	      //update with a new track and delete the previous track
	      TrackCont[trackid] = newtracks_fortest[i];
	      delete prev_track;
	    }
	    else{
	      prev_track -> SetClustersHoughFlag(GoodForTracking);
	      delete newtracks_fortest[i];
	    }
	  } //for scattered tracks

	  //Vertex finding again with new tracks
	  del::ClearContainer(VertexCont);
	  VertexSearch(TrackCont, VertexCont);
	  candidate_pairs.push_back(std::make_pair(candidates_trackid[0], candidates_trackid[1]));

	  break;
	} // else //reassaigning happens

	if(newtracks_fortest.size()!=candidates_trackid.size()) std::cout<<FUNC_NAME+" FATAL Error #of tracks is not matched"<<std::endl;
      }
    } //if(close_ids2.size() > 0 || close_ids1.size() > 0) case1
  }  // while(true)
  ResetHoughFlag(ClCont);

#if FragmentedTrackTest
  //Merged fragmented tracks
  RestoreFragmentedTracks(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif
  if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);

}

//_____________________________________________________________________________
template <typename T> void
FindAccidentalCoincidenceTracks(std::vector<T*>& TrackCont,
				std::vector<TPCVertex*>& VertexCont,
				std::vector<TPCVertex*>& ClusteredVertexCont)
{

  Double_t tgtXZ_cut = 30.; //mm
  Double_t target_section = 15; // abs(y)<target_section is not counted in this function

  TVector3 tgt(0., 0., tpc::ZTarget);
  std::vector<std::vector<TPCVertex*>> clustered_vertices; //Clustered tracks id
  std::vector<std::vector<Int_t>> clustered_tracks; //Clustered tracks id
  std::vector<TVector3> accidental_vertices;

  std::vector<Int_t> candidates_trackids;

  Int_t ntracks = TrackCont.size();
  for(Int_t trackid=0; trackid<ntracks; trackid++){
    T* ref_track = TrackCont[trackid];
    if(ref_track -> GetIsK18()==1 || ref_track -> GetIsKurama()==1) continue;
    if(ref_track -> GetIsBeam()!=1) continue;

    std::vector<TPCVertex*> candidate_clustered_vertices;
    std::vector<Int_t> candidate_clustered_tracks;
    std::vector<TVector3> candidate_clustered_pos;

    TVector3 vertex_point(0, 0, 0);
    candidate_clustered_tracks.push_back(trackid);
    candidate_clustered_pos.push_back(ref_track -> GetClosestPositionTgtXZ() - tgt);
    for(auto& vertex: VertexCont){
      //if(vertex -> GetClosestDist() > 15.) continue;
      if(vertex -> GetClosestDist() > 20.) continue;
      TVector3 vtx = vertex -> GetVertex() - tgt;
      if(TMath::Hypot(vtx.x(), vtx.z()) > tgtXZ_cut) continue;

      Int_t trackid1 = vertex -> GetTrackId(0);
      Int_t trackid2 = vertex -> GetTrackId(1);
      if(trackid2 == trackid){
	trackid1 = vertex -> GetTrackId(1);
	trackid2 = vertex -> GetTrackId(0);
      }
      else if(trackid1 != trackid) continue;

      T* coin_track = TrackCont[trackid2];
      if(coin_track -> GetIsK18()==1 || coin_track -> GetIsKurama()==1) continue;
      TVector3 residual_tgt = coin_track -> GetClosestPositionTgtXZ() - tgt;
      Double_t distXZ = TMath::Hypot(residual_tgt.x(), residual_tgt.z());
      if(distXZ > tgtXZ_cut) continue;

      candidate_clustered_vertices.push_back(vertex);
      candidate_clustered_tracks.push_back(trackid2);
      candidate_clustered_pos.push_back(residual_tgt);

      vertex_point += residual_tgt;
      candidates_trackids.push_back(trackid2);
    }

    if(candidate_clustered_tracks.size() < 2) continue;
    candidates_trackids.push_back(trackid);

    Double_t nvtx = candidate_clustered_tracks.size();

    clustered_vertices.push_back(candidate_clustered_vertices);
    clustered_tracks.push_back(candidate_clustered_tracks);
    accidental_vertices.push_back(TVector3(vertex_point.x()/nvtx,
					   vertex_point.y()/nvtx,
					   vertex_point.z()/nvtx));
  } //trackid

  std::sort(candidates_trackids.begin(), candidates_trackids.end());
  candidates_trackids.erase(std::unique(candidates_trackids.begin(), candidates_trackids.end()), candidates_trackids.end());

  for(Int_t trackid=0; trackid<ntracks; trackid++){
    if(find(candidates_trackids.begin(), candidates_trackids.end(), trackid) != candidates_trackids.end()) continue;

    T* ref_track = TrackCont[trackid];
    if(ref_track -> GetIsK18()==1 || ref_track -> GetIsKurama()==1) continue;
    if(ref_track -> GetIsBeam()==1) continue;

    TVector3 ref_residual_tgt = ref_track -> GetClosestPositionTgtXZ() - tgt;
    Double_t ref_distXZ = TMath::Hypot(ref_residual_tgt.x(), ref_residual_tgt.z());
    if(ref_distXZ > tgtXZ_cut) continue;

    std::vector<TPCVertex*> candidate_clustered_vertices;
    std::vector<Int_t> candidate_clustered_tracks;
    std::vector<TVector3> candidate_clustered_pos;

    TVector3 vertex_point(0, 0, 0);
    candidate_clustered_tracks.push_back(trackid);
    candidate_clustered_pos.push_back(ref_residual_tgt);
    for(auto& vertex: VertexCont){
      if(vertex -> GetClosestDist() > 15.) continue;
      TVector3 vtx = vertex -> GetVertex() - tgt;
      if(TMath::Hypot(vtx.x(), vtx.z()) > tgtXZ_cut) continue;

      Int_t trackid1 = vertex -> GetTrackId(0);
      Int_t trackid2 = vertex -> GetTrackId(1);

      if(trackid2 == trackid){
	trackid1 = vertex -> GetTrackId(1);
	trackid2 = vertex -> GetTrackId(0);
      }
      else if(trackid1 != trackid) continue;

      T* coin_track = TrackCont[trackid2];
      if(coin_track -> GetIsK18()==1 || coin_track -> GetIsKurama()==1) continue;
      if(coin_track -> GetIsBeam()==1) continue;
      TVector3 residual_tgt = coin_track -> GetClosestPositionTgtXZ() - tgt;
      Double_t distXZ = TMath::Hypot(residual_tgt.x(), residual_tgt.z());
      if(distXZ > tgtXZ_cut) continue;

      Double_t nvtx = candidate_clustered_tracks.size();
      TVector3 vertexpos(vertex_point.x()/nvtx,
			 vertex_point.y()/nvtx,
			 vertex_point.z()/nvtx);
      TVector3 resi = vertexpos - residual_tgt;
      candidate_clustered_vertices.push_back(vertex);
      candidate_clustered_tracks.push_back(trackid2);
      candidate_clustered_pos.push_back(residual_tgt);

      candidates_trackids.push_back(trackid2);

      vertex_point += residual_tgt;
    }

    if(candidate_clustered_tracks.size() < 2) continue;
    candidates_trackids.push_back(trackid);

    Double_t nvtx = candidate_clustered_tracks.size();

    clustered_vertices.push_back(candidate_clustered_vertices);
    clustered_tracks.push_back(candidate_clustered_tracks);
    accidental_vertices.push_back(TVector3(vertex_point.x()/nvtx,
					   vertex_point.y()/nvtx,
					   vertex_point.z()/nvtx));
  } //trackid

  for(Int_t id=0; id<accidental_vertices.size(); id++){
    if(TMath::Abs(accidental_vertices[id].y()) < target_section) continue;
    TPCVertex *clustered_vtx = new TPCVertex(accidental_vertices[id], clustered_tracks[id]);
    ClusteredVertexCont.push_back(clustered_vtx);
  }

  MarkingClusteredAccidentalTracks(TrackCont, ClusteredVertexCont);

}

//_____________________________________________________________________________
template <typename T> void
TestingCharge(std::vector<T*>& TrackCont,
	      std::vector<T*>& TrackContInvertedCharge,
	      std::vector<TPCVertex*>& VertexCont,
	      Bool_t Exclusive)
{

  Bool_t flag = false;
  TrackContInvertedCharge.resize(TrackCont.size());
  for(Int_t trackid=0; trackid<TrackCont.size(); trackid++){
    T *track = TrackCont[trackid];
    if(track -> GetIsAccidental()==1 || track -> GetIsK18()==1 || track -> GetIsBeam()==1 || track -> GetIsKurama()==1) continue;

    T *InvertedTrack = new T(track);
    Bool_t test = InvertedTrack -> TestInvertCharge();
    if(!test || (track -> GetCharge() == InvertedTrack -> GetCharge())){
      delete InvertedTrack;
      continue;
    }

    flag = true;
    InvertedTrack -> Calculate();
    if(Exclusive){
      InvertedTrack -> DoFitExclusive();
      InvertedTrack -> CalculateExclusive();
    }

    Int_t pid = track -> GetPid();
    Bool_t proton = ((pid&4)==4 && (pid&1)!=1);
    Bool_t pion = ((pid&4)!=4 && (pid&1)==1);

    Bool_t invert = false;
    if(!invert && track -> GetChiSquare() > InvertedTrack -> GetChiSquare()) invert = true;
    else{
      for(Int_t i=0; i<VertexCont.size(); i++){
	TPCVertex *vertex = VertexCont[i];
	Double_t closedist = vertex -> GetClosestDist();
	if(closedist < ppi_distcut){
	  Int_t trackid1 = vertex -> GetTrackId(0);
	  Int_t trackid2 = vertex -> GetTrackId(1);
	  if(trackid1 != trackid && trackid2 != trackid) continue;

	  T *track1 = TrackCont[trackid1];
	  T *track2 = TrackCont[trackid2];
	  if(trackid1 == trackid) track1 = InvertedTrack;
	  else if(trackid2 == trackid) track2 = InvertedTrack;

	  TPCVertex *newvertex = new TPCVertex(trackid1, trackid2);
	  newvertex -> Calculate(track1, track2);
	  if(newvertex -> GetIsLambda()) invert = true;
	  delete newvertex;
	}
      }
    }

    if(invert){
      TrackCont[trackid] = InvertedTrack;
      TrackContInvertedCharge[trackid] = track;
    }
    else{
      TrackCont[trackid] = track;
      TrackContInvertedCharge[trackid] = InvertedTrack;
    }
  }

  if(flag){
    //Vertex finding again with new tracks
    del::ClearContainer(VertexCont);
    VertexSearch(TrackCont, VertexCont);
  }
}

//_____________________________________________________________________________
template <typename T> void
ReassignClustersXiTrack(const std::vector<TPCClusterContainer>& ClCont,
			std::vector<T*>& TrackCont,
			std::vector<T*>& TrackContFailed,
			std::vector<TPCVertex*>& VertexCont,
			Bool_t Exclusive,
			Int_t MinNumOfHits)
{

  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);
  static const auto PionMass    = pdg::PionMass();
  static const auto LambdaMass  = pdg::LambdaMass();
  static const auto XiMinusMass = pdg::XiMinusMass();
  Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  Double_t xi_masscut = 0.15; //ref
  Double_t lpi_distcut = 15.; //ref

  //Xi- track is merged with p or pi- track (long-lived Xi- track)
  //Separating this into two tracks

  std::vector<Int_t> adjusted_trackid;
  std::vector<Int_t> erase_trackid;

  Bool_t flag = false;
  Int_t ntracks = TrackCont.size();
  for(Int_t trackid=0; trackid<ntracks; trackid++){
    if(std::find(adjusted_trackid.begin(), adjusted_trackid.end(), trackid) != adjusted_trackid.end()) continue;

    T *track = TrackCont[trackid];
    Int_t charge = track -> GetCharge();
    Int_t nhit = track -> GetNHit();

    if(track -> GetIsAccidental()==1) continue;
    if(track -> GetIsK18()==1) continue;
    if(track -> GetIsBeam()==1) continue;
    if(track -> GetIsKurama()==1) continue;

    //Xi- track is starting from the target
    Int_t start_section = track -> GetHitInOrder(0) -> GetSection();
    if(!track -> VertexAtTarget()) continue;
    if(start_section!=1) continue;

    //Xi- track goes forward
    TVector3 pos0 = track -> GetHitInOrder(0) -> GetLocalHitPos();
    if(pos0.z() < tpc::ZTarget) continue;

    std::vector<Int_t> candidates_trackid;
    std::vector<Int_t> candidates_charge;
    std::vector<Int_t> candidates_pid;
    std::vector<Double_t> candidates_helixtheta;
    std::vector<TVector3> candidates_vtx;
    for(auto& vertex: VertexCont){
      //Threshold conditions of a distance and an angle between two tracks.
      TVector3 vtx = vertex -> GetVertex();
      Int_t trackid1 = vertex -> GetTrackId(0);
      Int_t trackid2 = vertex -> GetTrackId(1);
      Double_t theta = vertex -> GetTrackTheta(0);
      if(std::find(adjusted_trackid.begin(), adjusted_trackid.end(), trackid1) != adjusted_trackid.end()) continue;
      if(std::find(adjusted_trackid.begin(), adjusted_trackid.end(), trackid2) != adjusted_trackid.end()) continue;

      if(trackid2==trackid){
	trackid1 = vertex -> GetTrackId(1);
	trackid2 = vertex -> GetTrackId(0);
	theta = vertex -> GetTrackTheta(1);
      }
      if(vertex -> GetClosestDist() > 10) continue;
      if(vtx.z() - tpc::ZTarget < 30) continue; //Xi- goes forward

      //vertex point exists inside of the track
      if(theta < track -> GetMint() || theta > track -> GetMaxt()) continue;

      T *track1 = TrackCont[trackid1];
      T *track2 = TrackCont[trackid2];
      if(track1 -> GetIsAccidental()==1 || track2 -> GetIsAccidental()==1) continue;
      if(track1 -> GetIsK18()==1 || track2 -> GetIsK18()==1) continue;
      if(track1 -> GetIsKurama()==1 || track2 -> GetIsKurama()==1) continue;
      if(track1 -> GetIsBeam()==1 || track2 -> GetIsBeam()==1) continue;

      candidates_trackid.push_back(trackid2);
      candidates_charge.push_back(track2 -> GetCharge());
      candidates_pid.push_back(track2 -> GetPid());
      candidates_helixtheta.push_back(theta); //track1's theta
      candidates_vtx.push_back(vtx);
    }

    Int_t segmented_trackid = -1;
    T *track_segmented;
    if(candidates_trackid.size()==3){

      Double_t theta_end = track -> GetHitInOrder(track -> GetNHit() - 1) -> GetTheta(); //Xi- track's closest point to the target
      T *track = TrackCont[trackid];

      for(Int_t candi=0; candi<candidates_vtx.size(); ++candi){
	Int_t close_hit; Double_t closedist = 9999;
	for(Int_t ihit=0; ihit<nhit; ++ihit){
	  //distance between cluster and vertex
	  TVector3 pos_diff = track->GetHitInOrder(ihit) -> GetLocalHitPos() - candidates_vtx[candi];
	  if(pos_diff.Mag() < closedist){
	    close_hit = ihit; closedist = pos_diff.Mag();
	  }
	}
	if(close_hit == nhit-1){
	  segmented_trackid = candidates_trackid[candi];
	  candidates_trackid.erase(candidates_trackid.begin() + candi);
	  candidates_charge.erase(candidates_charge.begin() + candi);
	  candidates_pid.erase(candidates_pid.begin() + candi);
	  candidates_helixtheta.erase(candidates_helixtheta.begin() + candi);
	  candidates_vtx.erase(candidates_vtx.begin() + candi);
	  track_segmented = TrackCont[segmented_trackid];
	}
      }
    }

    //Xi- track is merged with other track
    if(candidates_trackid.size()==2){
      Double_t theta0 = track -> GetHitInOrder(0) -> GetTheta(); //Xi- track's closest point to the target
      Int_t inner_trackid = candidates_trackid[0];
      Int_t outer_trackid = candidates_trackid[1];
      Int_t inner_pid = candidates_pid[0];
      Int_t outer_pid = candidates_pid[1];
      Int_t inner_charge = candidates_charge[0];
      Int_t outer_charge = candidates_charge[1];
      Double_t inner_theta = candidates_helixtheta[0];
      Double_t outer_theta = candidates_helixtheta[1];
      TVector3 inner_vtx = candidates_vtx[0];
      TVector3 outer_vtx = candidates_vtx[1];

      Double_t dist_track12 = TMath::Abs(track -> Getr()*(inner_theta - outer_theta));
      if(TMath::Abs(theta0 - candidates_helixtheta[0]) > TMath::Abs(theta0 - candidates_helixtheta[1])){
	inner_trackid = candidates_trackid[1];
	outer_trackid = candidates_trackid[0];
	inner_pid = candidates_pid[1];
	outer_pid = candidates_pid[0];
	inner_theta = candidates_helixtheta[1];
	outer_theta = candidates_helixtheta[0];
	inner_charge = candidates_charge[1];
	outer_charge = candidates_charge[0];
	inner_vtx = candidates_vtx[1];
	outer_vtx = candidates_vtx[0];
      }

      //order of tracks should be pi- and L decays
      if(inner_charge>0 && outer_charge>0) continue; //Two positive tracks are not Xi- decay.
      else if(inner_charge>0 && outer_charge<0){
	if(dist_track12 > 10) continue; //order is p, pi-, pi-. it's not a Xi-

	//make track order to pi-, p, pi- (first two tracks are too close)
	if(inner_trackid == candidates_trackid[0]){
	  inner_trackid = candidates_trackid[1];
	  outer_trackid = candidates_trackid[0];
	  inner_pid = candidates_pid[1];
	  outer_pid = candidates_pid[0];
	  inner_theta = candidates_helixtheta[1];
	  outer_theta = candidates_helixtheta[0];
	  inner_charge = candidates_charge[1];
	  outer_charge = candidates_charge[0];
	  inner_vtx = candidates_vtx[1];
	  outer_vtx = candidates_vtx[0];
	}
	else{
	  inner_trackid = candidates_trackid[0];
	  outer_trackid = candidates_trackid[1];
	  inner_pid = candidates_pid[0];
	  outer_pid = candidates_pid[1];
	  inner_theta = candidates_helixtheta[0];
	  outer_theta = candidates_helixtheta[1];
	  inner_charge = candidates_charge[0];
	  outer_charge = candidates_charge[1];
	  inner_vtx = candidates_vtx[0];
	  outer_vtx = candidates_vtx[1];
	}

	Bool_t proton = (outer_pid&4)==4;
	Bool_t pion = (inner_pid&1)==1;
	if(!proton || !pion) continue;
      }
      else if(inner_charge<0 && outer_charge>0){
	Bool_t proton = (outer_pid&4)==4;
	Bool_t pion = (inner_pid&1)==1;
	if(!proton || !pion) continue;
      }
      else{
	Bool_t pion1 = (inner_pid&1)==1;
	Bool_t pion2 = (outer_pid&1)==1;
	if(!pion1 || !pion2) continue;
      }

      T *track_in = TrackCont[inner_trackid];
      T *track_out = TrackCont[outer_trackid];
      std::vector<Int_t> xitrack_clusterids;
      std::vector<Int_t> ldecay_clusterids;
      for(Int_t i=0; i<track -> GetNHit(); ++i){
	TPCHit *hit = track -> GetHitInOrder(i) -> GetHit();
	Double_t theta = track -> GetHitInOrder(i) -> GetTheta();
	Int_t clusterid = track -> GetOrder(i);

	Double_t residual = 0;
	if(charge*theta <= charge*inner_theta || track_in -> IsGoodHitToAdd(hit, residual)) xitrack_clusterids.push_back(clusterid);
	else if(charge*theta > charge*outer_theta || track_out -> IsGoodHitToAdd(hit, residual)) ldecay_clusterids.push_back(clusterid);
	else break; //if cluster exists between lambda and pi- tracks, it is not Xi-
      }
      if(xitrack_clusterids.size()==0) continue;

      //if cluster exists between lambda and xi- tracks, it is not Xi-
      if((xitrack_clusterids.size() + ldecay_clusterids.size()) != track -> GetNHit()) continue;

      T *LdecayTrack = new T(track);
      LdecayTrack -> EraseHits(xitrack_clusterids);

      if(segmented_trackid != -1){
	for(Int_t hit=0;hit<track_segmented -> GetNHit();hit++){
	  LdecayTrack -> AddTPCHit(new TPCLTrackHit(track_segmented -> GetHitInOrder(hit) -> GetHit()));
	}
      }

      if(LdecayTrack -> DoFit(MinNumOfHits)){
	LdecayTrack -> Calculate();

	if(inner_charge<0 && outer_charge<0 && (LdecayTrack -> GetCharge()>0 || LdecayTrack -> TestInvertCharge())){ //pi-, pi-, p

	  Bool_t xi_reconstructed = false;
	  TPCVertex *Lvertex1 = new TPCVertex(outer_trackid, -1);
	  Lvertex1 -> Calculate(track_out, LdecayTrack);
	  if(Lvertex1 -> GetIsLambda()){
	    TVector3 lambda_vert; TVector3 lambda_mom;
	    for(int l=0; l<Lvertex1 -> GetNcombiLambda();l++){
	      if(Lvertex1 -> GetPionIdLambda(l) == outer_trackid){
		lambda_vert = Lvertex1 -> GetVertexLambda(l);
		lambda_mom = Lvertex1 -> GetMomLambda(l);
	      }
	    }

	    Double_t pi2_par[5];
	    track_in -> GetParam(pi2_par);
	    Double_t pi2_theta_min = track_in -> GetMint() - 150./pi2_par[3];
	    Double_t pi2_theta_max = track_in -> GetMaxt() + 150./pi2_par[3];

	    TVector3 pi2_mom; Double_t lpi_dist;
	    TVector3 xi_vert = Kinematics::XiVertex(dMagneticField,
						    pi2_par,
						    pi2_theta_min,
						    pi2_theta_max,
						    lambda_vert,
						    lambda_mom,
						    pi2_mom,
						    lpi_dist);

	    //Reconstructed xi decay vertex <-> pi2 & Xi track vertex
	    TVector3 check_vtx = inner_vtx - xi_vert;

	    TLorentzVector Lpi2(pi2_mom, TMath::Hypot(pi2_mom.Mag(), PionMass));
	    TLorentzVector Llambda_fixedmass(lambda_mom, TMath::Hypot(lambda_mom.Mag(), LambdaMass));
	    TLorentzVector Lxi = Llambda_fixedmass + Lpi2;
	    if(!TMath::IsNaN(lpi_dist) &&
	       lpi_dist < lpi_distcut &&
	       check_vtx.Mag() < lpi_distcut &&
	       TMath::Abs(Lxi.M() - XiMinusMass) < xi_masscut){
	      flag = true;
	      xi_reconstructed = true;
	    }
	  } //if(Lvertex1 -> GetIsLambda())
	  delete Lvertex1;

	  if(!xi_reconstructed && dist_track12 < 10){
	    TPCVertex *Lvertex2 = new TPCVertex(inner_trackid, -1);
	    Lvertex2 -> Calculate(track_in, LdecayTrack);
	    if(Lvertex2 -> GetIsLambda()){

	      TVector3 lambda_vert; TVector3 lambda_mom;
	      for(int l=0; l<Lvertex2 -> GetNcombiLambda();l++){
		if(Lvertex2 -> GetPionIdLambda(l) == inner_trackid){
		  lambda_vert = Lvertex2 -> GetVertexLambda(l);
		  lambda_mom = Lvertex2 -> GetMomLambda(l);
		}
	      }

	      Double_t pi2_par[5];
	      track_out -> GetParam(pi2_par);
	      Double_t pi2_theta_min = track_out -> GetMint() - 150./pi2_par[3];
	      Double_t pi2_theta_max = track_out -> GetMaxt() + 150./pi2_par[3];

	      TVector3 pi2_mom; Double_t lpi_dist;
	      TVector3 xi_vert = Kinematics::XiVertex(dMagneticField,
						      pi2_par,
						      pi2_theta_min,
						      pi2_theta_max,
						      lambda_vert,
						      lambda_mom,
						      pi2_mom,
						      lpi_dist);

	      //Reconstructed xi decay vertex <-> pi2 & Xi track vertex
	      TVector3 check_vtx = outer_vtx - xi_vert;

	      TLorentzVector Lpi2(pi2_mom, TMath::Hypot(pi2_mom.Mag(), PionMass));
	      TLorentzVector Llambda_fixedmass(lambda_mom, TMath::Hypot(lambda_mom.Mag(), LambdaMass));
	      TLorentzVector Lxi = Llambda_fixedmass + Lpi2;

	      if(!TMath::IsNaN(lpi_dist) &&
		 lpi_dist < lpi_distcut &&
		 check_vtx.Mag() < lpi_distcut &&
		 TMath::Abs(Lxi.M() - XiMinusMass) < xi_masscut){
		flag = true;
		xi_reconstructed = true;
	      }
	    } //if(Lvertex2 -> GetIsLambda())
	    delete Lvertex2;
	  }
	}
	else if(inner_charge<0 && outer_charge>0 && (LdecayTrack -> GetCharge()<0 || LdecayTrack -> TestInvertCharge())){ //pi-, p, pi-
	  TPCVertex *Lvertex = new TPCVertex(outer_trackid, -1);
	  Lvertex -> Calculate(track_out, LdecayTrack);
	  if(Lvertex -> GetIsLambda()){

	    TVector3 lambda_vert; TVector3 lambda_mom;
	    for(int l=0; l<Lvertex -> GetNcombiLambda();l++){
	      if(Lvertex -> GetProtonIdLambda(l) == outer_trackid){
		lambda_vert = Lvertex -> GetVertexLambda(l);
		lambda_mom = Lvertex -> GetMomLambda(l);
	      }
	    }

	    Double_t pi2_par[5];
	    track_in -> GetParam(pi2_par);
	    Double_t pi2_theta_min = track_in -> GetMint() - 150./pi2_par[3];
	    Double_t pi2_theta_max = track_in -> GetMaxt() + 150./pi2_par[3];

	    TVector3 pi2_mom; Double_t lpi_dist;
	    TVector3 xi_vert = Kinematics::XiVertex(dMagneticField,
						    pi2_par,
						    pi2_theta_min,
						    pi2_theta_max,
						    lambda_vert,
						    lambda_mom,
						    pi2_mom,
						    lpi_dist);

	    //Reconstructed xi decay vertex <-> pi2 & Xi track vertex
	    TVector3 check_vtx = inner_vtx - xi_vert;
	    TLorentzVector Lpi2(pi2_mom, TMath::Hypot(pi2_mom.Mag(), PionMass));
	    TLorentzVector Llambda_fixedmass(lambda_mom, TMath::Hypot(lambda_mom.Mag(), LambdaMass));
	    TLorentzVector Lxi = Llambda_fixedmass + Lpi2;

	    if(!TMath::IsNaN(lpi_dist) &&
	       lpi_dist < lpi_distcut &&
	       check_vtx.Mag() < lpi_distcut &&
	       TMath::Abs(Lxi.M() - XiMinusMass) < xi_masscut) flag = true;
	    else{
	      delete LdecayTrack;
	      continue;
	    }
	  } //if(Lvertex -> GetIsLambda())
	  delete Lvertex;
	} //else if(inner_charge<0 && outer_charge>0 && (LdecayTrack -> GetCharge()<0 || LdecayTrack -> TestInvertCharge())){ //pi-, p, pi-
	else{
	  delete LdecayTrack;
	  continue;
	}
      } //if(LdecayTrack -> DoFit(MinNumOfHits))
      else{
	delete LdecayTrack;
	continue;
      }

      if(Exclusive){
	LdecayTrack -> DoFitExclusive();
	LdecayTrack -> CalculateExclusive();
      }

      adjusted_trackid.push_back(trackid);
      adjusted_trackid.push_back(inner_trackid);
      adjusted_trackid.push_back(outer_trackid);

      track -> SetClustersHoughFlag(Candidate);
      if(segmented_trackid != -1){
	adjusted_trackid.push_back(segmented_trackid);
	erase_trackid.push_back(segmented_trackid);
	track_segmented -> SetClustersHoughFlag(Candidate);
      }

      T *XiTrack = new T(track);
      delete track;
      XiTrack -> EraseHits(ldecay_clusterids);
      TrackCont[trackid] = LdecayTrack;
      TrackCont[trackid] -> SetClustersHoughFlag(GoodForTracking);

      Int_t Xidecay_padid = tpc::findPadID(inner_vtx.z(), inner_vtx.x());
      Int_t Xidecay_layerid = tpc::getLayerID(Xidecay_padid);
      std::vector<TPCClusterContainer> ClContXi(NumOfLayersTPC);
      for(Int_t layer=0; layer<Xidecay_layerid+1; layer++){ //inner layers
	for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	  auto cl = ClCont[layer][ci];
	  TPCHit* hit = cl->GetMeanHit();
	  if(hit->GetHoughFlag()!=0 && hit->GetHoughFlag()!=1000) continue;
	  TVector3 pos = cl->GetPosition();
	  if(pos.Z()<tpc::ZTarget) continue;

	  Double_t xipos_x = inner_vtx.x()/(inner_vtx.z() - tpc::ZTarget)*(pos.z() - tpc::ZTarget);
	  Double_t xipos_y = inner_vtx.y()/(inner_vtx.z() - tpc::ZTarget)*(pos.z() - tpc::ZTarget);
	  if(TMath::Abs(xipos_x - pos.x())<20 && TMath::Abs(xipos_y - pos.y())<20) ClContXi[layer].push_back(cl);
	  if(layer<3 && TMath::Abs(xipos_x - pos.x())<5 && TMath::Abs(xipos_y - pos.y())<5) XiTrack->AddTPCHit(new TPCLTrackHit(hit));
	  else if(layer>=3 && TMath::Abs(xipos_x - pos.x())<10 && TMath::Abs(xipos_y - pos.y())<10) XiTrack->AddTPCHit(new TPCLTrackHit(hit));
	} //ci
      } //layer

      Int_t prev_size = TrackCont.size();
      FitTrack(XiTrack, GoodForTracking, ClContXi, TrackCont, TrackContFailed, 3);
      Int_t post_size = TrackCont.size();
      if(prev_size+1 == post_size){ //Fitting is succeeded
	adjusted_trackid.push_back(post_size-1);
	XiTrack = TrackCont[post_size-1];
	XiTrack -> Calculate();
	if(Exclusive){
	  XiTrack -> DoFitExclusive();
	  XiTrack -> CalculateExclusive();
	}
	XiTrack -> SetIsXi();
      }
    } // if(candidates_trackid.size()==2)
    else continue;
  }

  std::sort(erase_trackid.begin(), erase_trackid.end());
  for(Int_t i=0; i<erase_trackid.size(); ++i){
    Int_t id = erase_trackid[i] - i;

    T *track = TrackCont[id];
    TrackCont.erase(TrackCont.begin() + id);
    delete track;
  }

  ResetHoughFlag(ClCont);

  if(flag){
    //Vertex finding again with new tracks
    del::ClearContainer(VertexCont);
    VertexSearch(TrackCont, VertexCont);
  }
}

//_____________________________________________________________________________
template <typename T> void
RestoreFragmentedAccidentalTracks(const std::vector<TPCClusterContainer>& ClCont,
				  std::vector<T*>& TrackCont,
				  std::vector<T*>& TrackContFailed,
				  std::vector<TPCVertex*>& VertexCont,
				  Bool_t Exclusive,
				  Int_t MinNumOfHits)
{

  Bool_t reassaign = false;

  //case 1: accidental track is divided into a very short track on the upstream of the target and a track on the downstream of the target
  std::vector<Int_t> erase_trackid;
  std::vector<Int_t> new_trackid;
  for(Int_t trackid=0; trackid<TrackCont.size(); trackid++){
    T *track = TrackCont[trackid];
    if(track -> GetIsAccidental()==1) continue;
    if(track -> GetIsBeam()==1) continue;
    if(track -> GetIsK18()==1) continue;
    if(track -> GetIsKurama()==1) continue;
    if(std::find(erase_trackid.begin(), erase_trackid.end(), trackid)
       != erase_trackid.end()) continue;
    if(std::find(new_trackid.begin(), new_trackid.end(), trackid)
       != new_trackid.end()) continue;

    //candidates : negative charge mom > 0.5 or positive charge(flipped charge) with mom > 1.0
    Double_t slope = track -> Getdz();
    Double_t helixmom = track -> GetMom0().Mag();
    Int_t charge = track -> GetCharge();
    if(TMath::Abs(slope)>0.1 || TMath::Abs(helixmom)<0.5) continue;
    if(charge>0 && TMath::Abs(helixmom)<1.0) continue;

    for(auto& vertex: VertexCont){
      //Threshold conditions of a distance and an angle between two tracks.
      TVector3 vtx = vertex -> GetVertex();
      Int_t trackid1 = vertex -> GetTrackId(0);
      Int_t trackid2 = vertex -> GetTrackId(1);
      if(trackid1!=trackid && trackid2!=trackid) continue;
      if(trackid2==trackid){
	trackid2 = vertex -> GetTrackId(0);
	trackid1 = vertex -> GetTrackId(1);
      }
      if(std::find(erase_trackid.begin(), erase_trackid.end(), trackid1)
	 != erase_trackid.end()) continue;
      if(std::find(erase_trackid.begin(), erase_trackid.end(), trackid2)
	 != erase_trackid.end()) continue;
      if(std::find(new_trackid.begin(), new_trackid.end(), trackid1)
	 != new_trackid.end()) continue;
      if(std::find(new_trackid.begin(), new_trackid.end(), trackid2)
	 != new_trackid.end()) continue;

      T *CandidateTrack = TrackCont[trackid2];
      if(CandidateTrack -> GetIsK18()==1 ||
	 CandidateTrack -> GetIsKurama()==1) continue;
      if(CandidateTrack -> GetNHit()>4) continue; //this constraint has same effect as IsBackward()

      Bool_t cl_added = false;
      CandidateTrack -> SetClustersHoughFlag(0);
      TPCLocalTrackHelix *MergedTrack = new TPCLocalTrackHelix(track);
      for(Int_t hit=0;hit<CandidateTrack -> GetNHit();hit++){
	TPCHit* cl = CandidateTrack -> GetHitInOrder(hit) -> GetHit();

	Double_t resi=0.;
	Bool_t nolimitation = true;
	if(MergedTrack->IsGoodHitToAdd(cl, resi, nolimitation)){
	  MergedTrack -> AddTPCHit(new TPCLTrackHit(cl));
	  cl_added = true;
	}
      }

      //Check whether tracks belong to the same track or not
      if(!cl_added ||
	 (CandidateTrack -> GetNHit() + track -> GetNHit() - MergedTrack -> GetNHit()>2)){
	CandidateTrack -> SetClustersHoughFlag(GoodForTracking);
	delete MergedTrack;
	continue;
      }
      else{
	MergedTrack -> SetIsAccidental(); //it's needed to trun off VertexAtTarget()

	Int_t prev_size = TrackCont.size();
	FitTrack(MergedTrack, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);

	Int_t post_size = TrackCont.size();
	if(prev_size+1 == post_size){ //Fitting is succeeded
	  MergedTrack = TrackCont[post_size-1];
	  TrackCont.erase(TrackCont.begin() + post_size - 1);
	  MergedTrack -> Calculate();
	  MergedTrack -> CheckIsAccidental();
	  if(Exclusive){
	    MergedTrack -> DoFitExclusive();
	    MergedTrack -> CalculateExclusive();
	  }

	  if(MergedTrack -> GetNHit() > track -> GetNHit() &&
	     MergedTrack -> GetIsAccidental()==1){
	    TrackCont.push_back(MergedTrack);
	    new_trackid.push_back(post_size - 1);

	    erase_trackid.push_back(trackid1);
	    erase_trackid.push_back(trackid2);
	    reassaign = true;
	  }
	  else{
	    track -> SetClustersHoughFlag(GoodForTracking);
	    CandidateTrack -> SetClustersHoughFlag(GoodForTracking);
	    delete MergedTrack;
	  }
	}
	else{ //Fitting is failed
	  track -> SetClustersHoughFlag(GoodForTracking);
	  CandidateTrack -> SetClustersHoughFlag(GoodForTracking);
	}
      }
    }
  }

  std::sort(erase_trackid.begin(), erase_trackid.end());
  for(Int_t i=0; i<erase_trackid.size(); ++i){
    Int_t id = erase_trackid[i] - i;
    T *track = TrackCont[id];
    TrackCont.erase(TrackCont.begin() + id);
    delete track;
  }

  //case 2: accidental track is divided into clusters on the upstream of the target and a track on the downstream of the target
  std::vector<TPCClusterContainer> CandidateClCont(10);
  for(Int_t layer=0; layer<10; layer++){ //inner layers
    for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
      auto cl = ClCont[layer][ci];
      TPCHit* hit = cl->GetMeanHit();
      if(hit->GetHoughFlag()==GoodForTracking ||
	 hit->GetHoughFlag()==K18Tracks) continue;
      TVector3 pos = cl->GetPosition();
      if(pos.Z() > tpc::ZTarget) continue;
      if(TMath::Abs(pos.x()) > 20) continue;
      CandidateClCont[layer].push_back(cl);
    } //ci
  } //layer

  for(Int_t trackid=0; trackid<TrackCont.size(); trackid++){
    T *track = TrackCont[trackid];
    if(track -> GetIsAccidental()==1) continue;
    if(track -> GetIsBeam()==1) continue;
    if(track -> GetIsK18()==1) continue;
    if(track -> GetIsKurama()==1) continue;
    //if(!track -> VertexAtTarget()) continue;

    //candidates : negative charge mom > 0.5 or positive charge(flipped charge) with mom > 1.0
    Double_t slope = track -> Getdz();
    Double_t helixmom = track -> GetMom0().Mag();
    Int_t charge = track -> GetCharge();

    if(TMath::Abs(slope)>0.1 || TMath::Abs(helixmom)<0.5) continue;
    if(charge>0 && TMath::Abs(helixmom)<1.0) continue;

    T *CopiedTrack = new T(track);
    std::vector<TPCHit*> clusters_fortest;

    Bool_t cl_added = false;
    for(Int_t layer=0; layer<10; layer++){ //inner layers
      for(Int_t ci=0, n=CandidateClCont[layer].size(); ci<n; ci++){
	auto cl = CandidateClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	if(hit->GetHoughFlag()==GoodForTracking) continue;

	Double_t resi=0.;
	Bool_t nolimitation = true;
	if(track->IsGoodHitToAdd(hit, resi, nolimitation)){
	  CopiedTrack -> AddTPCHit(new TPCLTrackHit(hit));
	  clusters_fortest.push_back(hit);
	  cl_added = true;
	}
      }
    } //inner layers

    if(!cl_added) delete CopiedTrack;
    else{ //clusters are added
      CopiedTrack -> SetIsAccidental(); //it's needed to trun off VertexAtTarget()

      Int_t prev_size = TrackCont.size();
      FitTrack(CopiedTrack, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);

      Int_t post_size = TrackCont.size();
      if(prev_size+1 == post_size){ //Fitting is succeeded
	CopiedTrack = TrackCont[post_size-1];
	if(Exclusive){
	  CopiedTrack -> DoFitExclusive();
	  CopiedTrack -> CalculateExclusive();
	}
	TrackCont.erase(TrackCont.begin() + post_size - 1);
	CopiedTrack -> Calculate();
	CopiedTrack -> CheckIsAccidental();
	if(CopiedTrack -> GetNHit() > track -> GetNHit() &&
	   CopiedTrack -> GetIsAccidental()==1){
	  TrackCont[trackid] = CopiedTrack;
	  delete track;
	  reassaign = true;
	}
	else{
	  for(Int_t cl=0; cl<clusters_fortest.size(); cl++) clusters_fortest[cl] -> SetHoughFlag(0);
	  track -> SetClustersHoughFlag(GoodForTracking);
	  delete CopiedTrack;
	}
      }
      else{
	for(Int_t cl=0; cl<clusters_fortest.size(); cl++) clusters_fortest[cl] -> SetHoughFlag(0);
	track -> SetClustersHoughFlag(GoodForTracking);
      }
    } //clusters are added
  }

  if(reassaign){
    //Vertex finding again with new tracks
    del::ClearContainer(VertexCont);
    VertexSearch(TrackCont, VertexCont);
  } //reassaigning process
}

} //namespace tp

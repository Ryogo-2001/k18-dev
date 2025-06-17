// -*- C++ -*-

/*
  For on each TPC cluster, Virtual plane which is normal to track is defined.
  X position and Y position of each cluster are treated as a independent two points.

  ex) For each cluster of HypTPC, following equations are used for chisqr and ndf calculation
  chisqr += ((x_tpc-x_track)/sig_x_tpc)^2 + ((y_tpc-y_track)/sig_y_tpc)^2
  Ndf += 2*N_tpcclusters

  Kurama Tracking with the HypTPC
  : The end point of RK tracking is the target z position.
  So, clusters before the target center are excluded in the tracking. (Z_cluster < Z_target)
*/

#include "TPCRKTrack.hh"

#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <TMinuit.h>

#include <std_ostream.hh>

#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "PrintHelper.hh"
#include "TrackHit.hh"
#include "TPCLocalTrackHelix.hh"
#include "DeleteUtility.hh"

namespace
{
const auto& gGeom = DCGeomMan::GetInstance();
// const Int_t& IdTOF   = gGeom.DetectorId("TOF");
const Int_t& IdHS = gGeom.DetectorId("HS");
const Int_t& IdTgt = gGeom.DetectorId("Target");
const Int_t& IdTOFUX = gGeom.DetectorId("TOF-UX");
const Int_t& IdTOFUY = gGeom.DetectorId("TOF-UY");
const Int_t& IdTOFDX = gGeom.DetectorId("TOF-DX");
const Int_t& IdTOFDY = gGeom.DetectorId("TOF-DY");

const Int_t    MaxIteration = 100;
const Double_t InitialChiSqr = 1.e+10;
const Double_t MaxChiSqr     = 1.e+2;
const Double_t MinDeltaChiSqrR = 0.0002;
}

#define WARNOUT 1

//_____________________________________________________________________________
TPCRKTrack::TPCRKTrack(TPCLocalTrackHelix *tpctrack, std::vector<Int_t> lnum)
  : m_status(kInit),
    m_trackid(-1),
    m_track_in(),
    m_track_out(),
    m_track_tpc(tpctrack),
    m_tof_seg(-1.),
    m_initial_momentum(TMath::QuietNaN()),
    m_n_iteration(-1),
    m_nef_iteration(-1),
    m_chisqr(InitialChiSqr),
    m_polarity(0.),
    m_pid(1),
    m_primary_position(ThreeVector(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN())),
    m_primary_momentum(ThreeVector(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN())),
    m_last_position(ThreeVector(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN())),
    m_last_momentum(ThreeVector(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN())),
    m_path_length_tof(0.),
    m_path_length_tgt(0.),
    m_path_length_total(0.),
    m_tgt_position(ThreeVector(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN())),
    m_tgt_momentum(ThreeVector(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN())),
    m_tof_pos(ThreeVector(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN())),
    m_tof_mom(ThreeVector(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN())),
    m_is_good(true)
{
  s_status[kInit]                = "Initialized";
  s_status[kPassed]              = "Passed";
  s_status[kExceedMaxPathLength] = "Exceed Max Path Length";
  s_status[kExceedMaxStep]       = "Exceed Max Step";
  s_status[kFailedTraceLast]     = "Failed to Trace Last";
  s_status[kFailedGuess]         = "Failed to Guess";
  s_status[kFailedSave]          = "Failed to Save";
  s_status[kFatal]               = "Fatal";
  Initialize(lnum);
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCRKTrack::TPCRKTrack(TPCLocalTrackHelix *tpctrack, DCLocalTrack* track_in,
		       DCLocalTrack* track_out, std::vector<Int_t> lnum)
  : m_status(kInit),
    m_trackid(-1),
    m_track_in(),
    m_track_out(),
    m_track_tpc(tpctrack),
    m_tof_seg(-1.),
    m_initial_momentum(TMath::QuietNaN()),
    m_n_iteration(-1),
    m_nef_iteration(-1),
    m_chisqr(InitialChiSqr),
    m_polarity(0.),
    m_pid(1),
    m_primary_position(ThreeVector(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN())),
    m_primary_momentum(ThreeVector(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN())),
    m_last_position(ThreeVector(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN())),
    m_last_momentum(ThreeVector(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN())),
    m_path_length_tof(0.),
    m_path_length_tgt(0.),
    m_path_length_total(0.),
    m_tgt_position(ThreeVector(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN())),
    m_tgt_momentum(ThreeVector(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN())),
    m_tof_pos(ThreeVector(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN())),
    m_tof_mom(ThreeVector(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN())),
    m_is_good(true)
{
  s_status[kInit]                = "Initialized";
  s_status[kPassed]              = "Passed";
  s_status[kExceedMaxPathLength] = "Exceed Max Path Length";
  s_status[kExceedMaxStep]       = "Exceed Max Step";
  s_status[kFailedTraceLast]     = "Failed to Trace Last";
  s_status[kFailedGuess]         = "Failed to Guess";
  s_status[kFailedSave]          = "Failed to Save";
  s_status[kFatal]               = "Fatal";
  Initialize(track_out, track_out, lnum);
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCRKTrack::~TPCRKTrack()
{
  ClearHitArray();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
TPCRKTrack::Initialize(std::vector<Int_t> lnum)
{
  ClearHitArray();
  m_HitPointCont = RK::MakeHPContainer(lnum);
}

//_____________________________________________________________________________
void
TPCRKTrack::Initialize(DCLocalTrack* track_in, DCLocalTrack* track_out,
		       std::vector<Int_t> lnum)
{

  m_track_in = track_in;
  m_track_out = track_out;

  ClearHitArray();
  Int_t nIn = track_in->GetNHit();
  Int_t nOut = track_out->GetNHit();
  m_dchit_array.reserve(nIn+nOut);
  for(Int_t i=0; i<nIn; ++i){
    DCLTrackHit *hit  = track_in->GetHit(i);
    TrackHit    *thit = new TrackHit(hit);
    m_dchit_array.push_back(thit);
  }
  for(Int_t i=0; i<nOut; ++i){
    DCLTrackHit *hit  = track_out->GetHit(i);
    TrackHit    *thit = new TrackHit(hit);
    m_dchit_array.push_back(thit);
    if(hit->GetLayer() == IdTOFUX ||
       hit->GetLayer() == IdTOFDX){
      m_tof_seg = hit->GetWire();
    }
  }
  m_HitPointCont = RK::MakeHPContainer(lnum);
}

//_____________________________________________________________________________
void
TPCRKTrack::ClearHitArray()
{

  Int_t nh = m_dchit_array.size();
  for(Int_t i=nh-1; i>=0; --i){
    TrackHit *thit = m_dchit_array[i];
    DCHit *hit = thit -> GetHit() -> GetHit();
    delete hit;
    delete m_dchit_array[i];
  }
  m_dchit_array.clear();

}

//_____________________________________________________________________________
Bool_t
TPCRKTrack::DoFit(TVector3 initPos, TVector3 initMom, Bool_t U2D)
{

  m_status = kInit;
  m_initial_momentum = initMom.Mag();
  if(m_initial_momentum<0){
    hddaq::cout << FUNC_NAME << " initial momentum must be positive"
  		<< m_initial_momentum << std::endl;
    m_status = kFatal;
    return false;
  }

  RKCordParameter initCord(initPos, initMom);
  return DoFit(initCord, U2D);
}

//_____________________________________________________________________________
Bool_t
TPCRKTrack::DoFit(RKCordParameter initCord, Bool_t U2D)
{

  RKCordParameter     prevCord;
  RKHitPointContainer preHPntCont;

  Double_t chiSqr     = InitialChiSqr;
  Double_t prevChiSqr = InitialChiSqr;
  Double_t estDChisqr = InitialChiSqr;
  Double_t lambdaCri  = 0.01;
  Double_t dmp = 0.;

  RKHitPointContainer prevHPntCont;
  Int_t iItr = 0, iItrEf = 1;
  while(++iItr<MaxIteration){
    if(U2D) m_status = (RKstatus)RK::ExtrapTPC(m_track_tpc, initCord, m_HitPointCont, m_pid);
    else m_status = (RKstatus)RK::TraceTPC(m_track_tpc, initCord, m_HitPointCont, m_pid);
    if(m_status != kPassed){
#ifdef WARNOUT
      hddaq::cerr << FUNC_NAME << " "
		  << "something is wrong : " << iItr << std::endl;
#endif
      break;
    }

    chiSqr = CalcChiSqr(m_HitPointCont);
    Double_t dChiSqr  = chiSqr - prevChiSqr;
    Double_t dChiSqrR = dChiSqr/prevChiSqr;
    Double_t Rchisqr  = dChiSqr/estDChisqr;

#if 0
    {
      //PrintHelper helper(3, std::ios::scientific);
      hddaq::cout << FUNC_NAME << ": #"
		  << std::setw(3) << iItr << " ("
		  << std::setw(2) << iItrEf << ")"
		  << " chi=" << std::setw(10) << chiSqr;
      hddaq::cout.precision(5);
      hddaq::cout << " (" << std::fixed << std::setw(10) << dChiSqrR << ")"
		  << " [" << std::fixed << std::setw(10) << Rchisqr << " ]";
      //helper.precision(2);
      hddaq::cout << " df=" << std::setw(8) << dmp
		  << " (" << std::setw(8) << lambdaCri << ")" << std::endl;
      hddaq::cout << "initCord x:" << initCord.X()
		  << " y:" << std::setw(3) << initCord.Y() << " z:" << initCord.Z()
		  << " u:" << initCord.U() << " v:" << initCord.V()
		  << " p:" << 1./initCord.Q() << std::endl;
    }
#endif

#if 0
    PrintCalcHits(m_HitPointCont);
#endif

    if(std::abs(dChiSqrR)<MinDeltaChiSqrR &&
       (chiSqr<MaxChiSqr || Rchisqr>1.)){
      // Converged
      m_status = kPassed;
      if(dChiSqr>0.){
	initCord       = prevCord;
	chiSqr         = prevChiSqr;
	m_HitPointCont = prevHPntCont;
      }
      break;
    }

    // Next Guess
    if(iItr==1){
      prevCord     = initCord;
      prevChiSqr   = chiSqr;
      prevHPntCont = m_HitPointCont;
      ++iItrEf;
    }
    else if(dChiSqr <= 0.0){
      prevCord     = initCord;
      prevChiSqr   = chiSqr;
      prevHPntCont = m_HitPointCont;
      ++iItrEf;
      if(Rchisqr>=0.75){
	dmp*=0.5;
	if(dmp < lambdaCri) dmp=0.;
      }
      else if(Rchisqr>0.25){
	// nothing
      }
      else{
	if(dmp==0.) dmp=lambdaCri;
	else dmp*=2.;
      }
    }
    else {
      if(dmp==0.) dmp = lambdaCri;
      else {
	Double_t uf=2.;
	if(2.-Rchisqr > 2.) uf=2.-Rchisqr;
	dmp *= uf;
      }
      initCord       = prevCord;
      m_HitPointCont = prevHPntCont;
    }

    if(!GuessNextParameters(m_HitPointCont, initCord,
                            estDChisqr, lambdaCri, dmp)){
      hddaq::cout << FUNC_NAME << " "
		  << "cannot guess next paramters" << std::endl;
      m_status = kFailedGuess;
      return false;
    }
  }  /* End of Iteration */

  m_n_iteration   = iItr;
  m_nef_iteration = iItrEf;
  m_chisqr = chiSqr;

  if(m_track_tpc -> GetIsKurama()==1 && !RK::TraceToLast(m_HitPointCont)){
    m_status = kFailedTraceLast;
  }

  if(!SaveCalcPosition(m_HitPointCont) ||
     !SaveTrackParameters(initCord)){
    m_status = kFailedSave;
  }

#if 0
  Print("in "+FUNC_NAME);
#endif

  if(m_status != kPassed)
    return false;

  return true;
}

//_____________________________________________________________________________
Double_t
TPCRKTrack::CalcChiSqr(const RKHitPointContainer &hpCont) const
{

  Double_t chisqr=0.0;
  Int_t n=0;
  Int_t nh = m_dchit_array.size();
  for(Int_t i=0; i<nh; ++i){
    TrackHit *thp = m_dchit_array[i];
    if(!thp) continue;
    Int_t lnum = thp->GetLayer();
    const RKcalcHitPoint& calhp = hpCont.HitPointOfLayer(lnum);
    const ThreeVector& mom = calhp.MomentumInGlobal();
    Double_t w = gGeom.GetResolution(lnum);
    w = 1./(w*w);
    Double_t hitpos = thp->GetLocalHitPos();
    Double_t calpos = calhp.PositionInLocal();
    Double_t a = thp->GetTiltAngle()*TMath::DegToRad();
    Double_t u = mom.x()/mom.z();
    Double_t v = mom.y()/mom.z();
    Double_t dsdz = u*TMath::Cos(a)+v*TMath::Sin(a);
    Double_t coss = thp->IsHoneycomb() ? TMath::Cos(TMath::ATan(dsdz)) : 1.;
    Double_t wp   = thp->GetWirePosition();
    Double_t ss   = wp+(hitpos-wp)/coss;
    Double_t res  = (ss-calpos)*coss;
#if 0 //Checking L/R configuration again with RK tracking result.
    Double_t ss_pair   = wp-(hitpos-wp)/coss;
    Double_t res_pair  = (ss_pair-calpos)*coss;
    if(res_pair*res_pair<res*res){
      res = res_pair;
      thp->SetLocalHitPos(wp - (hitpos-wp)/coss);
    }
#endif
    chisqr += w*res*res;
    ++n;
  }

  Int_t nhTpc = m_track_tpc -> GetNHit();
  for(Int_t i=0; i<nhTpc; ++i){
    TPCLTrackHit *thp = m_track_tpc -> GetHitInOrder(i);
    if(!thp) continue;

    const TVector3& localhitpos = thp->GetLocalHitPos();
    Double_t TGTz = gGeom.GlobalZ("Target") - gGeom.GlobalZ("HS");
    if(m_track_tpc -> GetIsKurama()==1 && localhitpos.z()<TGTz) continue; //Exclude clusters before the target
    if(m_track_tpc -> GetIsK18()==1 && localhitpos.z()>TGTz) continue; //Exclude clusters after the target
    const TVector3& resolution = thp->GetResolutionVect();
    if(!thp->IsGoodForTracking()) continue; // exclude bad hits

    Int_t lnum = i + PlOffsTPCHit + 1;
    const RKcalcHitPoint& calhp_x = hpCont.HitPointOfLayer(lnum);
    Double_t callocalpos_x = calhp_x.PositionInLocal();
    if(TMath::IsNaN(callocalpos_x)) continue;
    chisqr += TMath::Power(callocalpos_x/TMath::Hypot(resolution.x(), resolution.z()), 2);
    ++n;

    const RKcalcHitPoint& calhp_y = hpCont.HitPointOfLayer(lnum + PlOffsTPCHit);
    Double_t callocalpos_y = calhp_y.PositionInLocal();
    if(TMath::IsNaN(callocalpos_y)) continue;
    chisqr += TMath::Power(callocalpos_y/resolution.y(), 2);
    ++n;
  }

  chisqr /= Double_t(n-5);
  return chisqr;
}

//_____________________________________________________________________________
Bool_t
TPCRKTrack::GuessNextParameters(const RKHitPointContainer& hpCont,
				RKCordParameter& Cord, Double_t& estDeltaChisqr,
				Double_t& lambdaCri, Double_t dmp) const
{
  Double_t *a2[10], a2c[10*5], *v[5],  vc[5*5];
  Double_t *a3[5],  a3c[5*5],  *v3[5], v3c[5*5], w3[5];
  Double_t dm[5];

  for(Int_t i=0; i<10; ++i){
    a2[i] =& a2c[5*i];
  }
  for(Int_t i=0; i<5; ++i){
    v[i]=&vc[5*i]; a3[i]=&a3c[5*i]; v3[i]=&v3c[5*i];
  }

  Double_t cb2[10], wSvd[5], dcb[5];
  Double_t wv[5]; // working space for SVD functions

  for(Int_t i=0; i<10; ++i){
    cb2[i]=0.0;
    for(Int_t j=0; j<5; ++j) a2[i][j]=0.0;
  }

  Int_t nth = 0;
  Int_t nh = m_dchit_array.size();
  for(Int_t i=0; i<nh; ++i){
    TrackHit *thp = m_dchit_array[i];
    if(!thp) continue;
    Int_t lnum = thp->GetLayer();
    const RKcalcHitPoint &calhp = hpCont.HitPointOfLayer(lnum);
    const ThreeVector&    mom   = calhp.MomentumInGlobal();
    Double_t hitpos = thp->GetLocalHitPos();
    Double_t calpos = calhp.PositionInLocal();
    Double_t a = thp->GetTiltAngle()*TMath::DegToRad();
    Double_t u = mom.x()/mom.z();
    Double_t v = mom.y()/mom.z();
    Double_t dsdz = u*TMath::Cos(a)+v*TMath::Sin(a);
    Double_t coss = thp->IsHoneycomb() ? TMath::Cos(TMath::ATan(dsdz)) : 1.;
    Double_t wp = thp->GetWirePosition();
    Double_t ss = wp+(hitpos-wp)/coss;
    Double_t cb = ss-calpos;
    Double_t wt = gGeom.GetResolution(lnum);
    wt = 1./(wt*wt);

    Double_t cfx=calhp.coefX(), cfy=calhp.coefY();
    Double_t cfu=calhp.coefU(), cfv=calhp.coefV(), cfq=calhp.coefQ();
    ++nth;

    cb2[0] += 2.*cfx*wt*cb;  cb2[1] += 2.*cfy*wt*cb;  cb2[2] += 2.*cfu*wt*cb;
    cb2[3] += 2.*cfv*wt*cb;  cb2[4] += 2.*cfq*wt*cb;

    a2[0][0] += 2.*wt*(cfx*cfx - cb*calhp.coefXX());
    a2[0][1] += 2.*wt*(cfx*cfy - cb*calhp.coefXY());
    a2[0][2] += 2.*wt*(cfx*cfu - cb*calhp.coefXU());
    a2[0][3] += 2.*wt*(cfx*cfv - cb*calhp.coefXV());
    a2[0][4] += 2.*wt*(cfx*cfq - cb*calhp.coefXQ());

    a2[1][0] += 2.*wt*(cfy*cfx - cb*calhp.coefYX());
    a2[1][1] += 2.*wt*(cfy*cfy - cb*calhp.coefYY());
    a2[1][2] += 2.*wt*(cfy*cfu - cb*calhp.coefYU());
    a2[1][3] += 2.*wt*(cfy*cfv - cb*calhp.coefYV());
    a2[1][4] += 2.*wt*(cfy*cfq - cb*calhp.coefYQ());

    a2[2][0] += 2.*wt*(cfu*cfx - cb*calhp.coefUX());
    a2[2][1] += 2.*wt*(cfu*cfy - cb*calhp.coefUY());
    a2[2][2] += 2.*wt*(cfu*cfu - cb*calhp.coefUU());
    a2[2][3] += 2.*wt*(cfu*cfv - cb*calhp.coefUV());
    a2[2][4] += 2.*wt*(cfu*cfq - cb*calhp.coefUQ());

    a2[3][0] += 2.*wt*(cfv*cfx - cb*calhp.coefVX());
    a2[3][1] += 2.*wt*(cfv*cfy - cb*calhp.coefVY());
    a2[3][2] += 2.*wt*(cfv*cfu - cb*calhp.coefVU());
    a2[3][3] += 2.*wt*(cfv*cfv - cb*calhp.coefVV());
    a2[3][4] += 2.*wt*(cfv*cfq - cb*calhp.coefVQ());

    a2[4][0] += 2.*wt*(cfq*cfx - cb*calhp.coefQX());
    a2[4][1] += 2.*wt*(cfq*cfy - cb*calhp.coefQY());
    a2[4][2] += 2.*wt*(cfq*cfu - cb*calhp.coefQU());
    a2[4][3] += 2.*wt*(cfq*cfv - cb*calhp.coefQV());
    a2[4][4] += 2.*wt*(cfq*cfq - cb*calhp.coefQQ());
  }

  //TPC XZ
  Int_t nhtpc = m_track_tpc -> GetNHit();
  for(Int_t i=0; i<nhtpc; ++i){
    TPCLTrackHit *thp = m_track_tpc -> GetHitInOrder(i);
    if(!thp) continue;

    const TVector3& localhitpos = thp->GetLocalHitPos();
    Double_t TGTz = gGeom.GlobalZ("Target") - gGeom.GlobalZ("HS");
    if(m_track_tpc -> GetIsKurama()==1 && localhitpos.z()<TGTz) continue; //Exclude clusters before the target
    if(m_track_tpc -> GetIsK18()==1 && localhitpos.z()>TGTz) continue; //Exclude clusters after the target
    const TVector3& resolution = thp->GetResolutionVect();
    if(!thp->IsGoodForTracking()) continue; // exclude bad hits

    Int_t lnum = i + PlOffsTPCHit + 1;
    const RKcalcHitPoint& calhp = hpCont.HitPointOfLayer(lnum);
    if(TMath::IsNaN(calhp.PositionInLocal())) continue;
    Double_t cfx=calhp.coefX(), cfy=calhp.coefY();
    Double_t cfu=calhp.coefU(), cfv=calhp.coefV(), cfq=calhp.coefQ();
    Double_t wt = 1./(resolution.x()*resolution.x() + resolution.z()*resolution.z());
    Double_t cb = 0 - calhp.PositionInLocal(); //the center of the virtual plane is 0;
    ++nth;

    cb2[0] += 2.*cfx*wt*cb;  cb2[1] += 2.*cfy*wt*cb;  cb2[2] += 2.*cfu*wt*cb;
    cb2[3] += 2.*cfv*wt*cb;  cb2[4] += 2.*cfq*wt*cb;

    a2[0][0] += 2.*wt*(cfx*cfx - cb*calhp.coefXX());
    a2[0][1] += 2.*wt*(cfx*cfy - cb*calhp.coefXY());
    a2[0][2] += 2.*wt*(cfx*cfu - cb*calhp.coefXU());
    a2[0][3] += 2.*wt*(cfx*cfv - cb*calhp.coefXV());
    a2[0][4] += 2.*wt*(cfx*cfq - cb*calhp.coefXQ());

    a2[1][0] += 2.*wt*(cfy*cfx - cb*calhp.coefYX());
    a2[1][1] += 2.*wt*(cfy*cfy - cb*calhp.coefYY());
    a2[1][2] += 2.*wt*(cfy*cfu - cb*calhp.coefYU());
    a2[1][3] += 2.*wt*(cfy*cfv - cb*calhp.coefYV());
    a2[1][4] += 2.*wt*(cfy*cfq - cb*calhp.coefYQ());

    a2[2][0] += 2.*wt*(cfu*cfx - cb*calhp.coefUX());
    a2[2][1] += 2.*wt*(cfu*cfy - cb*calhp.coefUY());
    a2[2][2] += 2.*wt*(cfu*cfu - cb*calhp.coefUU());
    a2[2][3] += 2.*wt*(cfu*cfv - cb*calhp.coefUV());
    a2[2][4] += 2.*wt*(cfu*cfq - cb*calhp.coefUQ());

    a2[3][0] += 2.*wt*(cfv*cfx - cb*calhp.coefVX());
    a2[3][1] += 2.*wt*(cfv*cfy - cb*calhp.coefVY());
    a2[3][2] += 2.*wt*(cfv*cfu - cb*calhp.coefVU());
    a2[3][3] += 2.*wt*(cfv*cfv - cb*calhp.coefVV());
    a2[3][4] += 2.*wt*(cfv*cfq - cb*calhp.coefVQ());

    a2[4][0] += 2.*wt*(cfq*cfx - cb*calhp.coefQX());
    a2[4][1] += 2.*wt*(cfq*cfy - cb*calhp.coefQY());
    a2[4][2] += 2.*wt*(cfq*cfu - cb*calhp.coefQU());
    a2[4][3] += 2.*wt*(cfq*cfv - cb*calhp.coefQV());
    a2[4][4] += 2.*wt*(cfq*cfq - cb*calhp.coefQQ());
  }

  //TPC Y terms
  for(Int_t i=0; i<nhtpc; ++i){
    TPCLTrackHit *thp = m_track_tpc -> GetHitInOrder(i);
    if(!thp) continue;

    const TVector3& localhitpos = thp->GetLocalHitPos();
    Double_t TGTz = gGeom.GlobalZ("Target") - gGeom.GlobalZ("HS");
    if(m_track_tpc -> GetIsKurama()==1 && localhitpos.z()<TGTz) continue; //Exclude clusters before the target
    if(m_track_tpc -> GetIsK18()==1 && localhitpos.z()>TGTz) continue; //Exclude clusters after the target
    const TVector3& resolution = thp->GetResolutionVect();
    if(!thp->IsGoodForTracking()) continue; // exclude bad hits

    Int_t lnum = i + 2.*PlOffsTPCHit + 1;
    const RKcalcHitPoint& calhp = hpCont.HitPointOfLayer(lnum);
    if(TMath::IsNaN(calhp.PositionInLocal())) continue;
    Double_t cfx=calhp.coefX(), cfy=calhp.coefY();
    Double_t cfu=calhp.coefU(), cfv=calhp.coefV(), cfq=calhp.coefQ();
    Double_t wt = 1./(resolution.y()*resolution.y());
    Double_t cb = 0 - calhp.PositionInLocal(); //the center of the virtual plane is 0;

    ++nth;

    cb2[0] += 2.*cfx*wt*cb;  cb2[1] += 2.*cfy*wt*cb;  cb2[2] += 2.*cfu*wt*cb;
    cb2[3] += 2.*cfv*wt*cb;  cb2[4] += 2.*cfq*wt*cb;

    a2[0][0] += 2.*wt*(cfx*cfx - cb*calhp.coefXX());
    a2[0][1] += 2.*wt*(cfx*cfy - cb*calhp.coefXY());
    a2[0][2] += 2.*wt*(cfx*cfu - cb*calhp.coefXU());
    a2[0][3] += 2.*wt*(cfx*cfv - cb*calhp.coefXV());
    a2[0][4] += 2.*wt*(cfx*cfq - cb*calhp.coefXQ());

    a2[1][0] += 2.*wt*(cfy*cfx - cb*calhp.coefYX());
    a2[1][1] += 2.*wt*(cfy*cfy - cb*calhp.coefYY());
    a2[1][2] += 2.*wt*(cfy*cfu - cb*calhp.coefYU());
    a2[1][3] += 2.*wt*(cfy*cfv - cb*calhp.coefYV());
    a2[1][4] += 2.*wt*(cfy*cfq - cb*calhp.coefYQ());

    a2[2][0] += 2.*wt*(cfu*cfx - cb*calhp.coefUX());
    a2[2][1] += 2.*wt*(cfu*cfy - cb*calhp.coefUY());
    a2[2][2] += 2.*wt*(cfu*cfu - cb*calhp.coefUU());
    a2[2][3] += 2.*wt*(cfu*cfv - cb*calhp.coefUV());
    a2[2][4] += 2.*wt*(cfu*cfq - cb*calhp.coefUQ());

    a2[3][0] += 2.*wt*(cfv*cfx - cb*calhp.coefVX());
    a2[3][1] += 2.*wt*(cfv*cfy - cb*calhp.coefVY());
    a2[3][2] += 2.*wt*(cfv*cfu - cb*calhp.coefVU());
    a2[3][3] += 2.*wt*(cfv*cfv - cb*calhp.coefVV());
    a2[3][4] += 2.*wt*(cfv*cfq - cb*calhp.coefVQ());

    a2[4][0] += 2.*wt*(cfq*cfx - cb*calhp.coefQX());
    a2[4][1] += 2.*wt*(cfq*cfy - cb*calhp.coefQY());
    a2[4][2] += 2.*wt*(cfq*cfu - cb*calhp.coefQU());
    a2[4][3] += 2.*wt*(cfq*cfv - cb*calhp.coefQV());
    a2[4][4] += 2.*wt*(cfq*cfq - cb*calhp.coefQQ());
  }

  for(Int_t i=0; i<5; ++i)
    for(Int_t j=0; j<5; ++j)
      a3[i][j]=a2[i][j];

  // Levenberg-Marqardt method
  Double_t lambda = std::sqrt(dmp);
  //  a2[5][0]=a2[6][1]=a2[7][2]=a2[8][3]=a2[9][4]=lambda;

  for(Int_t ii=0; ii<5; ++ii){
    dm[ii]       = a2[ii][ii];
    a2[ii+5][ii] = lambda * std::sqrt(a2[ii][ii]);
  }

#if 0
  {
    PrintHelper helper(3, std::ios::scientific);
    hddaq::cout << FUNC_NAME << ": A2 and CB2 before SVDcmp"
		<<  std::endl;
    for(Int_t ii=0; ii<10; ++ii)
      hddaq::cout << std::setw(12) << a2[ii][0] << ","
		  << std::setw(12) << a2[ii][1] << ","
		  << std::setw(12) << a2[ii][2] << ","
		  << std::setw(12) << a2[ii][3] << ","
		  << std::setw(12) << a2[ii][4] << "  "
		  << std::setw(12) << cb2[ii] << std::endl;
  }
#endif

  // Solve the Eq. with SVD (Singular Value Decomposition) Method
  if(!MathTools::SVDcmp(a2, 10, 5, wSvd, v, wv))
    return false;

#if 0
  {
    PrintHelper helper(3, std::ios::scientific);
    hddaq::cout << FUNC_NAME << ": A2 after SVDcmp"
		<<  std::endl;
    for(Int_t ii=0; ii<10; ++ii)
      hddaq::cout << std::setw(12) << a2[ii][0] << ","
		  << std::setw(12) << a2[ii][1] << ","
		  << std::setw(12) << a2[ii][2] << ","
		  << std::setw(12) << a2[ii][3] << ","
		  << std::setw(12) << a2[ii][4] << std::endl;
  }
#endif

#if 0
  // check orthogonality of decomposted matrics
  {
    PrintHelper helper(5, std::ios::scientific);
    hddaq::cout << FUNC_NAME << ": Check V*~V" <<  std::endl;
    for(Int_t i=0; i<5; ++i){
      for(Int_t j=0; j<5; ++j){
	Double_t f=0.0;
	for(Int_t k=0; k<5; ++k)
	  f += v[i][k]*v[j][k];
	hddaq::cout << std::setw(10) << f;
      }
      hddaq::cout << std::endl;
    }
    hddaq::cout << FUNC_NAME << ": Check U*~U" <<  std::endl;
    for(Int_t i=0; i<10; ++i){
      for(Int_t j=0; j<10; ++j){
	Double_t f=0.0;
	for(Int_t k=0; k<5; ++k)
	  f += a2[i][k]*a2[j][k];
	hddaq::cout << std::setw(10) << f;
      }
      hddaq::cout << std::endl;
    }

    hddaq::cout << FUNC_NAME << ": Check ~U*U" <<  std::endl;
    for(Int_t i=0; i<5; ++i){
      for(Int_t j=0; j<5; ++j){
	Double_t f=0.0;
	for(Int_t k=0; k<10; ++k)
	  f += a2[k][i]*a2[k][j];
	hddaq::cout << std::setw(10) << f;
      }
      hddaq::cout << std::endl;
    }
  }
#endif

  Double_t wmax=0.0;
  for(Int_t i=0; i<5; ++i)
    if(wSvd[i]>wmax) wmax=wSvd[i];

  Double_t wmin=wmax*1.E-15;
  for(Int_t i=0; i<5; ++i)
    if(wSvd[i]<wmin) wSvd[i]=0.0;

#if 0
  {
    PrintHelper helper(3, std::ios::scientific);
    hddaq::cout << FUNC_NAME << ": V and Wsvd after SVDcmp"
		<<  std::endl;
    for(Int_t ii=0; ii<5; ++ii)
      hddaq::cout << std::setw(12) << v[ii][0] << ","
		  << std::setw(12) << v[ii][1] << ","
		  << std::setw(12) << v[ii][2] << ","
		  << std::setw(12) << v[ii][3] << ","
		  << std::setw(12) << v[ii][4] << "  "
		  << std::setw(12) << wSvd[ii] << std::endl;
  }
#endif

  MathTools::SVDksb(a2, wSvd, v, 10, 5, cb2, dcb, wv);

#if 0
  {
    PrintHelper helper(5, std::ios::scientific);
    hddaq::cout << FUNC_NAME << ": "
		<< std::setw(12) << Cord.Z();
    hddaq::cout << "  Dumping Factor = " << std::setw(12) << dmp << std::endl;
    helper.setf(std::ios::fixed);
    hddaq::cout << std::setw(12) << Cord.X() << "  "
		<< std::setw(12) << dcb[0] << " ==>  "
		<< std::setw(12) << Cord.X()+dcb[0] << std::endl;
    hddaq::cout << std::setw(12) << Cord.Y() << "  "
		<< std::setw(12) << dcb[1] << " ==>  "
		<< std::setw(12) << Cord.Y()+dcb[1] << std::endl;
    hddaq::cout << std::setw(12) << Cord.U() << "  "
		<< std::setw(12) << dcb[2] << " ==>  "
		<< std::setw(12) << Cord.U()+dcb[2] << std::endl;
    hddaq::cout << std::setw(12) << Cord.V() << "  "
		<< std::setw(12) << dcb[3] << " ==>  "
		<< std::setw(12) << Cord.V()+dcb[3] << std::endl;
    hddaq::cout << std::setw(12) << Cord.Q() << "  "
		<< std::setw(12) << dcb[4] << " ==>  "
		<< std::setw(12) << Cord.Q()+dcb[4] << std::endl;
  }
#endif

  Cord = RKCordParameter(Cord.X()+dcb[0],
                         Cord.Y()+dcb[1],
                         Cord.Z(),
                         Cord.U()+dcb[2],
                         Cord.V()+dcb[3],
                         Cord.Q()+dcb[4]);

  // calc. the critical dumping factor & est. delta-ChiSqr
  Double_t s1=0., s2=0.;
  for(Int_t i=0; i<5; ++i){
    s1 += dcb[i]*dcb[i]; s2 += dcb[i]*cb2[i];
  }
  estDeltaChisqr = (-s2-dmp*s1)/Double_t(nth-5);

  if(!MathTools::SVDcmp(a3, 5, 5, w3, v3, wv))
    return false;

  Double_t spur=0.;
  for(Int_t i=0; i<5; ++i){
    Double_t s=0.;
    for(Int_t j=0; j<5; ++j)
      s += v3[i][j]*a3[i][j];
    if(w3[i]!=0.0)
      spur += s/w3[i]*dm[i];
  }

  lambdaCri = 1./spur;

  return true;
}

//_____________________________________________________________________________
Bool_t
TPCRKTrack::SaveCalcPosition(const RKHitPointContainer &hpCont)
{
  for(Int_t i=0, n=m_dchit_array.size(); i<n; ++i){
    TrackHit *thp = m_dchit_array[i];
    if(!thp) continue;
    Int_t lnum = thp->GetLayer();
    const RKcalcHitPoint &calhp = hpCont.HitPointOfLayer(lnum);
    thp->SetCalGPos(calhp.PositionInGlobal());
    thp->SetCalGMom(calhp.MomentumInGlobal());
    thp->SetCalLPos(calhp.PositionInLocal());
  }

  return true;
}

//_____________________________________________________________________________
void
TPCRKTrack::Print(const TString& arg, std::ostream& ost)
{
  PrintHelper helper(5, std::ios::fixed, ost);

  ost << FUNC_NAME << " " << arg << std::endl
      << "   status : " << s_status[m_status] << std::endl
      << " in " << std::setw(3) << m_n_iteration << " ("
      << std::setw(2) << m_nef_iteration << ") Iteractions "
      << " chisqr=" << std::setw(10) << m_chisqr << " init=" << m_initial_momentum << std::endl;
  ost << " Primary X (" << std::setprecision(2)
      << std::setw(7) << m_primary_position.x() << ", "
      << std::setw(7) << m_primary_position.y() << ", "
      << std::setw(7) << m_primary_position.z() << ")"
      << " U " << std::setprecision(2) << std::setw(7)
      << m_primary_momentum.x()/m_primary_momentum.z()
      << " V " << std::setprecision(2) << std::setw(7)
      << m_primary_momentum.x()/m_primary_momentum.z() << std::endl;
  ost << "        P " << std::setprecision(5)
      << std::setw(7) << m_primary_momentum.Mag() << " ("
      << std::setw(7) << m_primary_momentum.x() << ", "
      << std::setw(7) << m_primary_momentum.y() << ", "
      << std::setw(7) << m_primary_momentum.z() << ")" << std::endl;
  ost << " End    X (" << std::setprecision(2)
      << std::setw(7) << m_last_position.x() << ", "
      << std::setw(7) << m_last_position.y() << ", "
      << std::setw(7) << m_last_position.z() << ")" << std::endl;
  ost << "        P " << std::setprecision(5)
      << std::setw(7) << m_last_momentum.Mag() << " ("
      << std::setw(7) << m_last_momentum.x() << ", "
      << std::setw(7) << m_last_momentum.y() << ", "
      << std::setw(7) << m_last_momentum.z() << ")" << std::endl;
  ost << " Target X " << std::setprecision(2)
      << std::setw(7) << m_tgt_position.x() << ", "
      << " Y " << std::setprecision(2)
      << std::setw(7) << m_tgt_position.y() << ", "
      << " U " << std::setprecision(2)
      << std::setw(7) << m_tgt_momentum.x()/m_tgt_momentum.z() << ", "
      << " V " << std::setprecision(2)
      << std::setw(7) << m_tgt_momentum.y()/m_tgt_momentum.z() << std::endl;
  ost << "        P " << std::setprecision(5)
      << std::setw(7) << m_tgt_momentum.Mag() << " ("
      << std::setw(7) << m_tgt_momentum.x() << ", "
      << std::setw(7) << m_tgt_momentum.y() << ", "
      << std::setw(7) << m_tgt_momentum.z() << ")"
      << "        PathLength  " << std::setprecision(1)
      << std::setw(7) << m_path_length_tgt << std::endl;
  if(m_tof_seg>=0){
    ost << " TOF    X (" << std::setprecision(2)
	<< std::setw(7) << m_tof_pos.x() << ", "
	<< std::setw(7) << m_tof_pos.y() << ", "
	<< std::setw(7) << m_tof_pos.z() << ")" << std::endl;
    ost << "        P " << std::setprecision(5)
	<< std::setw(7) << m_tof_mom.Mag() << " ("
	<< std::setw(7) << m_tof_mom.x() << ", "
	<< std::setw(7) << m_tof_mom.y() << ", "
	<< std::setw(7) << m_tof_mom.z() << ")"
	<< "        PathLength  " << std::setprecision(1)
	<< std::setw(7) << m_path_length_tof << " "
	<< "TOF#" << std::setw(2) << std::right
	<< (Int_t)m_tof_seg << std::endl;
  }
  ost << "total path " << std::setw(7) << m_path_length_total << std::endl;

  PrintCalcHits(m_HitPointCont);
}

//_____________________________________________________________________________
void
TPCRKTrack::PrintCalcHits(const RKHitPointContainer &hpCont, std::ostream &ost) const
{
  PrintHelper helper(2, std::ios::fixed, ost);

  const Int_t n = m_dchit_array.size();
  RKHitPointContainer::RKHpCIterator itr, end = hpCont.end();
  for(itr=hpCont.begin(); itr!=end; ++itr){
    const RKcalcHitPoint &calhp = itr->second;
    ThreeVector pos = calhp.PositionInGlobal();
    std::string h = " ";

    Int_t lnum = itr->first;
    if(lnum > PlOffsTPCHit && lnum < 2.*PlOffsTPCHit){
      Int_t clusterId = lnum - PlOffsTPCHit - 1;
      TPCLTrackHit *thp = m_track_tpc -> GetHitInOrder(clusterId);
      const TVector3& localhitpos = thp->GetLocalHitPos();
      ThreeVector residual = gGeom.Local2GlobalPos(IdHS, localhitpos);
      residual -= pos;
      if (thp){ h = "*"; }
      ost << "#"   << std::setw(2) << lnum << h
	  << " L " << std::setw(8) << calhp.PathLength()
	  << " X " << std::setw(7) << "None"
	  << " ("  << std::setw(7) << pos.x()
	  << ", "  << std::setw(7) << pos.y()
	  << ", "  << std::setw(8) << pos.z()<< ")";
      ost << " residual (x, y, z) Local " << std::setw(7) << thp->GetResidualVect().x()
	  << " "<< std::setw(7) << thp->GetResidualVect().y()
	  << " "<< std::setw(7) << thp->GetResidualVect().z()
	  << " Calc. " << std::setw(7) << residual.x()
	  << " "<< std::setw(7) << residual.y()
	  << " "<< std::setw(7) << residual.y();
    }
    else{
      TrackHit *thp = 0;
      for(Int_t i=0; i<n; ++i){
	if(m_dchit_array[i] && m_dchit_array[i]->GetLayer()==lnum){
	  thp = m_dchit_array[i];
	  if(thp) break;
	}
      }
      if (thp){ h = "-"; if (thp->IsHoneycomb()) h = "+"; }
      ost << "#"   << std::setw(2) << lnum << h
	  << " L " << std::setw(8) << calhp.PathLength()
	  << " X " << std::setw(7) << calhp.PositionInLocal()
	  << " (" << std::setw(7) << pos.x()
	  << ", "  << std::setw(7) << pos.y()
	  << ", "  << std::setw(8) << pos.z() << ")";
      if(thp){
	ost << " local pos "   << std::setw(7) << thp->GetLocalHitPos()
	    << " residual " << std::setw(7) << thp->GetLocalHitPos() - calhp.PositionInLocal();
      }
    }
    ost << std::endl;
#if 0
    {
      ThreeVector mom = calhp.MomentumInGlobal();
      ost << "   P=" << std::setw(7) << mom.Mag()
	  << " (" << std::setw(9) << mom.x()
	  << ", " << std::setw(9) << mom.y()
	  << ", " << std::setw(9) << mom.z()
	  << ")" << std::endl;
    }
#endif
  }
}

//_____________________________________________________________________________
Bool_t
TPCRKTrack::SaveTrackParameters(const RKCordParameter &cp)
{

  Int_t dcFirst = -1;
  RKHitPointContainer::RKHpCIterator itr, end = m_HitPointCont.end();
  for(itr=m_HitPointCont.begin(); itr!=end; ++itr){
    Int_t lnum = itr->first;
    if(lnum < PlOffsTPCHit){
      dcFirst = lnum;
      break;
    }
  }
  if(dcFirst < 0) return false;

  const Int_t idFirst = dcFirst;
  const RKcalcHitPoint& hpFirst = m_HitPointCont.HitPointOfLayer(idFirst);
  const ThreeVector& posFirst = hpFirst.PositionInGlobal();
  const ThreeVector& momFirst = hpFirst.MomentumInGlobal();
  m_primary_position = gGeom.Global2LocalPos(idFirst, posFirst);
  m_primary_momentum = gGeom.Global2LocalDir(idFirst, momFirst);
  m_polarity = m_primary_momentum.z()<0. ? -1. : 1.;

  const Int_t idLast = m_HitPointCont.rbegin()->first;
  const RKcalcHitPoint& hpLast = m_HitPointCont.rbegin()->second;
  const ThreeVector& posLast = hpLast.PositionInGlobal();
  const ThreeVector& momLast = hpLast.MomentumInGlobal();
  m_last_position = gGeom.Global2LocalPos(idLast, posLast);
  m_last_momentum = gGeom.Global2LocalDir(idLast, momLast);
  m_path_length_total = std::abs(hpFirst.PathLength()-hpLast.PathLength());

  const RKcalcHitPoint& hpTgt = m_HitPointCont.HitPointOfLayer(IdTgt);
  const ThreeVector& posTgt = hpTgt.PositionInGlobal();
  const ThreeVector& momTgt = hpTgt.MomentumInGlobal();
  m_tgt_position = gGeom.Global2LocalPos(IdTgt, posTgt);
  m_tgt_momentum = gGeom.Global2LocalDir(IdTgt, momTgt);
  m_path_length_tgt = std::abs(hpFirst.PathLength()-hpTgt.PathLength());

  if(m_tof_seg>=0){
    const RKcalcHitPoint& hpTofU = m_HitPointCont.HitPointOfLayer(IdTOFUX);
    const RKcalcHitPoint& hpTofD = m_HitPointCont.HitPointOfLayer(IdTOFDX);
    if((Int_t)m_tof_seg%2==0){ // upstream
      m_path_length_tof = std::abs(hpTgt.PathLength()-hpTofU.PathLength());
      m_tof_pos = hpTofU.PositionInGlobal();
      m_tof_mom = hpTofU.MomentumInGlobal();
    }
    else if((Int_t)m_tof_seg%2==1){ // downstream
      m_path_length_tof = std::abs(hpTgt.PathLength()-hpTofD.PathLength());
      m_tof_pos = hpTofD.PositionInGlobal();
      m_tof_mom = hpTofD.MomentumInGlobal();
    }
    else{
      m_path_length_tof = std::abs(hpTgt.PathLength()-(hpTofU.PathLength()+hpTofD.PathLength())/2.);
      m_tof_pos = (hpTofU.PositionInGlobal()+hpTofD.PositionInGlobal())*0.5;
      m_tof_mom = (hpTofU.MomentumInGlobal()+hpTofD.MomentumInGlobal())*0.5;
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
TPCRKTrack::GetTrajectoryLocalPosition(Int_t layer, Double_t& path, Double_t& x, Double_t& y) const
{
  try {
    ThreeVector lpos;
    const RKcalcHitPoint& hpTgt = m_HitPointCont.begin()->second;
    const RKcalcHitPoint& HP = m_HitPointCont.HitPointOfLayer(layer);
    const ThreeVector& gpos = HP.PositionInGlobal();
    lpos = gGeom.Global2LocalPos(layer, gpos);

    path = std::abs(hpTgt.PathLength()-HP.PathLength());
    x = lpos.x();
    y = lpos.y();
    return true;
  }
  catch(const std::out_of_range&) {
    return false;
  }
}

//_____________________________________________________________________________
Bool_t
TPCRKTrack::GetTrajectoryLocalPosition(Int_t layer, Double_t& x, Double_t& y) const
{

  Double_t path;
  return GetTrajectoryLocalPosition(layer, path, x, y);
}

//_____________________________________________________________________________
Bool_t
TPCRKTrack::GetTrajectoryResidualTPC(Int_t clusterId, TVector3& resolution, TVector3& residual) const
{
  try{
    TPCLTrackHit *thp = m_track_tpc -> GetHitInOrder(clusterId);
    if(!thp) return false;
    const TVector3& localhitpos = thp->GetLocalHitPos();
    Double_t TGTz = gGeom.GlobalZ("Target") - gGeom.GlobalZ("HS");
    if(m_track_tpc -> GetIsKurama()==1 && localhitpos.z()<TGTz) return false; //Exclude clusters before the target
    if(m_track_tpc -> GetIsK18()==1 && localhitpos.z()>TGTz) return false; //Exclude clusters after the target
    resolution = thp->GetResolutionVect();
    if(!thp->IsGoodForTracking()) return false;

    Int_t lnum = clusterId + PlOffsTPCHit + 1;
    const RKcalcHitPoint& calhp = m_HitPointCont.HitPointOfLayer(lnum);
    ThreeVector calpos = calhp.PositionInGlobal();
    residual = gGeom.Local2GlobalPos(IdHS, localhitpos);
    residual -= calpos;

    return true;
  }
  catch(const std::out_of_range&) {
    return false;
  }
}

//_____________________________________________________________________________
Bool_t
TPCRKTrack::GetTrajectoryResidual(Int_t layer, Double_t& resolution, Double_t& residual) const
{
  try{
    Bool_t status = false;
    Int_t nh = m_dchit_array.size();
    for(Int_t i=0; i<nh; ++i){
      TrackHit *thp = m_dchit_array[i];
      if(!thp) continue;
      Int_t lnum = thp->GetLayer();
      if(lnum!=layer) continue;
      status = true;
      const RKcalcHitPoint& calhp = m_HitPointCont.HitPointOfLayer(lnum);
      const ThreeVector& mom = calhp.MomentumInGlobal();
      Double_t hitpos = thp->GetLocalHitPos();
      Double_t calpos = calhp.PositionInLocal();
      Double_t a = thp->GetTiltAngle()*TMath::DegToRad();
      Double_t u = mom.x()/mom.z();
      Double_t v = mom.y()/mom.z();
      Double_t dsdz = u*TMath::Cos(a)+v*TMath::Sin(a);
      Double_t coss = thp->IsHoneycomb() ? TMath::Cos(TMath::ATan(dsdz)) : 1.;
      Double_t wp   = thp->GetWirePosition();
      Double_t ss   = wp+(hitpos-wp)/coss;
      resolution = gGeom.GetResolution(lnum);
      residual = (ss-calpos)*coss;
    }
    return status;
  }
  catch(const std::out_of_range&) {
    return false;
  }
}

//_____________________________________________________________________________
Bool_t
TPCRKTrack::GetTrajectoryGlobalPosition(Int_t layer, TVector3& global_pos) const
{
  try {
    const RKcalcHitPoint& HP = m_HitPointCont.HitPointOfLayer(layer);
    const ThreeVector& gpos = HP.PositionInGlobal();
    global_pos = gpos;
    return true;
  }
  catch(const std::out_of_range&) {
    return false;
  }
}

//_____________________________________________________________________________
Bool_t
TPCRKTrack::GetTrajectoryGlobalPositionTPC(Int_t clusterId, TVector3& global_pos) const
{
  try {
    TPCLTrackHit *thp = m_track_tpc -> GetHitInOrder(clusterId);
    if(!thp) return false;
    const TVector3& localhitpos = thp->GetLocalHitPos();
    Double_t TGTz = gGeom.GlobalZ("Target") - gGeom.GlobalZ("HS");
    if(m_track_tpc -> GetIsKurama()==1 && localhitpos.z()<TGTz) return false; //Exclude clusters before the target
    if(m_track_tpc -> GetIsK18()==1 && localhitpos.z()>TGTz) return false; //Exclude clusters after the target
    if(!thp->IsGoodForTracking()) return false;

    Int_t lnum = clusterId + PlOffsTPCHit + 1;
    const RKcalcHitPoint& HP = m_HitPointCont.HitPointOfLayer(lnum);
    const ThreeVector& gpos = HP.PositionInGlobal();
    global_pos = gpos;
    return true;
  }
  catch(const std::out_of_range&) {
    return false;
  }
}

//_____________________________________________________________________________
Bool_t
TPCRKTrack::GetTrajectoryMomentum(Int_t layer, TVector3& global_mom) const
{
  try {
    const RKcalcHitPoint& HP = m_HitPointCont.HitPointOfLayer(layer);
    const ThreeVector& gmom = HP.MomentumInGlobal();
    global_mom = gmom;
    return true;
  }
  catch(const std::out_of_range&) {
    return false;
  }
}

//_____________________________________________________________________________
Bool_t
TPCRKTrack::GetTrajectoryMomentumTPC(Int_t clusterId, TVector3& global_mom) const
{
  try {
    TPCLTrackHit *thp = m_track_tpc -> GetHitInOrder(clusterId);
    if(!thp) return false;
    const TVector3& localhitpos = thp->GetLocalHitPos();
    Double_t TGTz = gGeom.GlobalZ("Target") - gGeom.GlobalZ("HS");
    if(m_track_tpc -> GetIsKurama()==1 && localhitpos.z()<TGTz) return false; //Exclude clusters before the target
    if(m_track_tpc -> GetIsK18()==1 && localhitpos.z()>TGTz) return false; //Exclude clusters after the target
    if(!thp->IsGoodForTracking()) return false;

    Int_t lnum = clusterId + PlOffsTPCHit + 1;
    const RKcalcHitPoint& HP = m_HitPointCont.HitPointOfLayer(lnum);
    const ThreeVector& gmom = HP.MomentumInGlobal();
    global_mom = gmom;
    return true;
  }
  catch(const std::out_of_range&) {
    return false;
  }
}

// -*- C++ -*-

#include "KuramaTrack.hh"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <TMinuit.h>

#include <std_ostream.hh>

#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "PrintHelper.hh"
#include "TrackHit.hh"

namespace
{
const auto& gGeom = DCGeomMan::GetInstance();
// const Int_t& IdTOF   = gGeom.DetectorId("TOF");
const Int_t& IdTOFUX = gGeom.DetectorId("TOF-UX");
const Int_t& IdTOFDX = gGeom.DetectorId("TOF-DX");
const Int_t    MaxIteraction = 100;
const Double_t InitialChiSqr = 1.e+10;
const Double_t MaxChiSqr     = 1.e+2;
const Double_t MinDeltaChiSqrR = 0.0002;
}

#define WARNOUT 0

//_____________________________________________________________________________
KuramaTrack::KuramaTrack(DCLocalTrack *track_in, DCLocalTrack *track_out)
  : m_status(kInit),
    m_track_in(track_in),
    m_track_out(track_out),
    m_tof_seg(-1),
    m_initial_momentum(0.),
    m_n_iteration(-1),
    m_nef_iteration(-1),
    m_chisqr(InitialChiSqr),
    m_polarity(0.),
    m_path_length_tof(0.),
    m_path_length_total(0.),
    m_tof_pos(ThreeVector(0., 0., 0.)),
    m_tof_mom(ThreeVector(0., 0., 0.)),
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
  FillHitArray();
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
KuramaTrack::~KuramaTrack()
{
  ClearHitArray();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
TrackHit*
KuramaTrack::GetHitOfLayerNumber(Int_t lnum) const
{
  for(Int_t i=0, n=m_hit_array.size(); i<n; ++i){
    if(m_hit_array[i]->GetLayer()==lnum)
      return m_hit_array[i];
  }
  return 0;
}

//_____________________________________________________________________________
void
KuramaTrack::FillHitArray()
{
  Int_t nIn  = m_track_in->GetNHit();
  Int_t nOut = m_track_out->GetNHit();
  ClearHitArray();
  m_hit_array.reserve(nIn+nOut);

  for(Int_t i=0; i<nIn; ++i){
    DCLTrackHit *hit  = m_track_in->GetHit(i);
    TrackHit    *thit = new TrackHit(hit);
    m_hit_array.push_back(thit);
  }

  for(Int_t i=0; i<nOut; ++i){
    DCLTrackHit *hit  = m_track_out->GetHit(i);
    TrackHit    *thit = new TrackHit(hit);
    m_hit_array.push_back(thit);
    if(hit->GetLayer() == IdTOFUX ||
       hit->GetLayer() == IdTOFDX){
      m_tof_seg = hit->GetWire();
    }
  }

}

//_____________________________________________________________________________
void
KuramaTrack::ClearHitArray()
{
  Int_t nh = m_hit_array.size();
  for(Int_t i=nh-1; i>=0; --i){
    delete m_hit_array[i];
  }
  m_hit_array.clear();
}

//_____________________________________________________________________________
Bool_t
KuramaTrack::ReCalc(Bool_t applyRecursively)
{
  if(applyRecursively){
    m_track_in->ReCalc(applyRecursively);
    m_track_out->ReCalc(applyRecursively);
    //    Int_t nh=m_hit_array.size();
    //    for(Int_t i=0; i<nh; ++i){
    //      m_hit_array[i]->ReCalc(applyRecursively);
    //    }
  }

  return DoFit(m_cord_param);
  //  return DoFit();
}

//_____________________________________________________________________________
Bool_t
KuramaTrack::DoFit()
{
  m_status = kInit;

  if(m_initial_momentum<0){
    hddaq::cout << FUNC_NAME << " initial momentum must be positive"
  		<< m_initial_momentum << std::endl;
    m_status = kFatal;
    return false;
  }

  static const auto gposTOFDX = gGeom.GetGlobalPosition("TOF-DX");
  const Double_t xOut    = m_track_out->GetX(gposTOFDX.z());
  const Double_t yOut    = m_track_out->GetY(gposTOFDX.z());
  const Double_t uOut    = m_track_out->GetU0();
  const Double_t vOut    = m_track_out->GetV0();
  const Double_t pzOut   = m_initial_momentum/std::sqrt(1.+uOut*uOut+vOut*vOut);
  const ThreeVector posOut(xOut, yOut, gposTOFDX.z());
  const ThreeVector momOut(pzOut*uOut, pzOut*vOut, pzOut);

  RKCordParameter     iniCord(posOut, momOut);
  RKCordParameter     prevCord;
  RKHitPointContainer preHPntCont;

  Double_t chiSqr     = InitialChiSqr;
  Double_t prevChiSqr = InitialChiSqr;
  Double_t estDChisqr = InitialChiSqr;
  Double_t lambdaCri  = 0.01;
  Double_t dmp = 0.;

  m_HitPointCont = RK::MakeHPContainer();
  RKHitPointContainer prevHPntCont;

  Int_t iItr = 0, iItrEf = 1;

  while(++iItr<MaxIteraction){
    m_status = (RKstatus)RK::Trace(iniCord, m_HitPointCont);
    if(m_status != kPassed){
#ifdef WARNOUT
      // hddaq::cerr << FUNC_NAME << " "
      // 		<< "something is wrong : " << iItr << std::endl;
#endif
      break;
    }

    chiSqr = CalcChiSqr(m_HitPointCont);
    Double_t dChiSqr  = chiSqr - prevChiSqr;
    Double_t dChiSqrR = dChiSqr/prevChiSqr;
    Double_t Rchisqr  = dChiSqr/estDChisqr;
#if 0
    {
      PrintHelper helper(3, std::ios::scientific);
      hddaq::cout << FUNC_NAME << ": #"
		  << std::setw(3) << iItr << " ("
		  << std::setw(2) << iItrEf << ")"
		  << " chi=" << std::setw(10) << chiSqr;
      hddaq::cout.precision(5);
      hddaq::cout << " (" << std::fixed << std::setw(10) << dChiSqrR << ")"
		  << " [" << std::fixed << std::setw(10) << Rchisqr << " ]";
      helper.precision(2);
      hddaq::cout << " df=" << std::setw(8) << dmp
		  << " (" << std::setw(8) << lambdaCri << ")" << std::endl;
      hddaq::cout << "iniCord x:" << iniCord.X()
		  << " y:" << iniCord.Y() << " z:" << iniCord.Z()
		  << " u:" << iniCord.U() << " v:" << iniCord.V()
		  << " p:" << 1./iniCord.Q() << std::endl;
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
	iniCord        = prevCord;
	chiSqr         = prevChiSqr;
	m_HitPointCont = prevHPntCont;
      }
      break;
    }

    // Next Guess
    if(iItr==1){
      prevCord     = iniCord;
      prevChiSqr   = chiSqr;
      prevHPntCont = m_HitPointCont;
      ++iItrEf;
    }
    else if(dChiSqr <= 0.0){
      prevCord     = iniCord;
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
      iniCord        = prevCord;
      m_HitPointCont = prevHPntCont;
    }

    if(!GuessNextParameters(m_HitPointCont, iniCord,
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

  if(!RK::TraceToLast(m_HitPointCont)){
    m_status = kFailedTraceLast;
  }
  if(!SaveCalcPosition(m_HitPointCont) ||
     !SaveTrackParameters(iniCord)){
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
Bool_t
KuramaTrack::DoFit(RKCordParameter iniCord)
{
  //  ClearHitArray();
  //  FillHitArray();

  RKCordParameter prevCord;
  RKHitPointContainer preHPntCont;

  Double_t chiSqr     = InitialChiSqr;
  Double_t prevChiSqr = InitialChiSqr;
  Double_t estDChisqr = InitialChiSqr;
  Double_t lambdaCri  = 0.01;
  Double_t dmp = 0.;
  m_status = kInit;

  m_HitPointCont = RK::MakeHPContainer();
  RKHitPointContainer prevHPntCont;

  Int_t iItr=0, iItrEf=1;

  while(++iItr < MaxIteraction){
    if(!RK::Trace(iniCord, m_HitPointCont)){
      // Error
#ifdef WARNOUT
      //hddaq::cerr << FUNC_NAME << ": Error in RK::Trace. " << std::endl;
#endif
      break;
    }

    chiSqr = CalcChiSqr(m_HitPointCont);
    Double_t dChiSqr  = chiSqr-prevChiSqr;
    Double_t dChiSqrR = dChiSqr/prevChiSqr;
    Double_t Rchisqr  = dChiSqr/estDChisqr;
#if 0
    {
      PrintHelper helper(3, std::ios::scientific);
      hddaq::cout << FUNC_NAME << ": #"
		  << std::setw(3) << iItr << " ("
		  << std::setw(2) << iItrEf << ")"
		  << " chi=" << std::setw(10) << chiSqr;
      helper.precision(5);
      hddaq::cout << " (" << std::fixed << std::setw(10) << dChiSqrR << ")"
		  << " [" << std::fixed << std::setw(10) << Rchisqr << " ]";
      helper.precision(2);
      hddaq::cout << " df=" << std::setw(8) << dmp
		  << " (" << std::setw(8) << lambdaCri << ")" << std::endl;
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
	iniCord        = prevCord;
	chiSqr         = prevChiSqr;
	m_HitPointCont = prevHPntCont;
      }
      break;
    }

    // Next Guess
    if(iItr==1){
      prevCord     = iniCord;
      prevChiSqr   = chiSqr;
      prevHPntCont = m_HitPointCont;
      ++iItrEf;
    }
    else if(dChiSqr <= 0.0){
      prevCord     = iniCord;
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
	if(dmp==0.0) dmp=lambdaCri;
	else dmp*=2.;
      }
    }
    else {
      if(dmp==0.0) dmp=lambdaCri;
      else {
	Double_t uf=2.;
	if(2.-Rchisqr > 2.) uf=2.-Rchisqr;
	dmp *= uf;
      }
      iniCord = prevCord;
      m_HitPointCont = prevHPntCont;
    }
    if(!GuessNextParameters(m_HitPointCont, iniCord,
                            estDChisqr, lambdaCri, dmp)){
      m_status = kFailedGuess;
      return false;
    }
  }  /* End of Iteration */

  m_n_iteration   = iItr;
  m_nef_iteration = iItrEf;
  m_chisqr        = chiSqr;

  if(!RK::TraceToLast(m_HitPointCont))
    m_status = kFailedTraceLast;

  if(!SaveCalcPosition(m_HitPointCont) ||
     !SaveTrackParameters(iniCord))
    m_status = kFailedSave;

#if 0
  Print("in "+FUNC_NAME);
#endif

  if(m_status != kPassed)
    return false;

  return true;
}

//_____________________________________________________________________________
Double_t
KuramaTrack::CalcChiSqr(const RKHitPointContainer &hpCont) const
{
  Int_t nh = m_hit_array.size();

  Double_t chisqr=0.0;
  Int_t n=0;

  for(Int_t i=0; i<nh; ++i){
    TrackHit *thp = m_hit_array[i];
    if(!thp) continue;
    Int_t    lnum = thp->GetLayer();
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
#if 0 //Checking L/R ambiguity again with RK tracking result.
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
  chisqr /= Double_t(n-5);
  return chisqr;
}

//_____________________________________________________________________________
Bool_t
KuramaTrack::GuessNextParameters(const RKHitPointContainer& hpCont,
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
  Double_t wv[5];   // working space for SVD functions

  Int_t nh = m_hit_array.size();

  for(Int_t i=0; i<10; ++i){
    cb2[i]=0.0;
    for(Int_t j=0; j<5; ++j) a2[i][j]=0.0;
  }

  Int_t nth=0;
  for(Int_t i=0; i<nh; ++i){
    TrackHit *thp = m_hit_array[i];
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
    Double_t wp   = thp->GetWirePosition();
    Double_t ss   = wp+(hitpos-wp)/coss;
    Double_t cb   = ss-calpos;
    // Double_t cb = thp->GetLocalHitPos()-calhp.PositionInLocal();
    // Double_t cb = thp->GetResidual();
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
KuramaTrack::SaveCalcPosition(const RKHitPointContainer &hpCont)
{
  for(Int_t i=0, n=m_hit_array.size(); i<n; ++i){
    TrackHit *thp = m_hit_array[i];
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
KuramaTrack::Print(const TString& arg, std::ostream& ost)
{
  PrintHelper helper(5, std::ios::fixed, ost);

  ost << FUNC_NAME << " " << arg << std::endl
      << "   status : " << s_status[m_status] << std::endl
      << " in " << std::setw(3) << m_n_iteration << " ("
      << std::setw(2) << m_nef_iteration << ") Iteractions "
      << " chisqr=" << std::setw(10) << m_chisqr << std::endl;
  ost << " Target X (" << std::setprecision(2)
      << std::setw(7) << m_primary_position.x() << ", "
      << std::setw(7) << m_primary_position.y() << ", "
      << std::setw(7) << m_primary_position.z() << ")" << std::endl
      << "        P " << std::setprecision(5)
      << std::setw(7) << m_primary_momentum.Mag() << " ("
      << std::setw(7) << m_primary_momentum.x() << ", "
      << std::setw(7) << m_primary_momentum.y() << ", "
      << std::setw(7) << m_primary_momentum.z() << ")"
      << " init=" << m_initial_momentum
      << std::endl
      << "        PathLength  " << std::setprecision(1)
      << std::setw(7) << m_path_length_tof << " "
      << "TOF#" << std::setw(2) << std::right
      << (Int_t)m_tof_seg << std::endl;
  //  << std::setw(7) << m_path_length_total << std::endl;

  PrintCalcHits(m_HitPointCont);
}

//_____________________________________________________________________________
void
KuramaTrack::PrintCalcHits(const RKHitPointContainer &hpCont, std::ostream &ost) const
{
  PrintHelper helper(2, std::ios::fixed, ost);

  const Int_t n = m_hit_array.size();
  RKHitPointContainer::RKHpCIterator itr, end = hpCont.end();
  for(itr=hpCont.begin(); itr!=end; ++itr){
    Int_t lnum = itr->first;
    TrackHit *thp = 0;
    for(Int_t i=0; i<n; ++i){
      if(m_hit_array[i] && m_hit_array[i]->GetLayer()==lnum){
	thp = m_hit_array[i];
	if(thp) break;
      }
    }
    const RKcalcHitPoint &calhp = itr->second;
    ThreeVector pos = calhp.PositionInGlobal();
    std::string h = " ";
    if (thp){ h = "-"; if (thp->IsHoneycomb()) h = "+"; }
    ost << "#"   << std::setw(2) << lnum << h
	<< " L " << std::setw(8) << calhp.PathLength()
	<< " X " << std::setw(7) << calhp.PositionInLocal()
	<< " (" << std::setw(7) << pos.x()
	<< ", "  << std::setw(7) << pos.y()
	<< ", "  << std::setw(8) << pos.z() << ")";
    if(thp){
      ost << " "   << std::setw(7) << thp->GetLocalHitPos()
	  << " -> " << std::setw(7)
	  << thp->GetResidual();
      //<< thp->GetLocalHitPos()-calhp.PositionInLocal();
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
KuramaTrack::SaveTrackParameters(const RKCordParameter &cp)
{
  m_cord_param = cp;
  const Int_t TGTid = m_HitPointCont.begin()->first;

  // const RKcalcHitPoint& hpTof  = m_HitPointCont.HitPointOfLayer(IdTOF);

  const RKcalcHitPoint& hpTgt  = m_HitPointCont.begin()->second;
  const RKcalcHitPoint& hpLast = m_HitPointCont.rbegin()->second;

  const ThreeVector& pos = hpTgt.PositionInGlobal();
  const ThreeVector& mom = hpTgt.MomentumInGlobal();

  m_primary_position = gGeom.Global2LocalPos(TGTid, pos);
  m_primary_momentum = gGeom.Global2LocalDir(TGTid, mom);

  m_polarity = m_primary_momentum.z()<0. ? -1. : 1.;

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
  else {
    m_path_length_tof = std::abs(hpTgt.PathLength()-(hpTofU.PathLength()+hpTofD.PathLength())/2.);
    m_tof_pos = (hpTofU.PositionInGlobal()+hpTofD.PositionInGlobal())*0.5;
    m_tof_mom = (hpTofU.MomentumInGlobal()+hpTofD.MomentumInGlobal())*0.5;
  }

  m_path_length_total = std::abs(hpTgt.PathLength()-hpLast.PathLength());

  return true;
}

//_____________________________________________________________________________
Bool_t
KuramaTrack::GetTrajectoryLocalPosition(Int_t layer,
                                        Double_t& path, Double_t& x, Double_t& y) const
{
  try {
    const RKcalcHitPoint& hpTgt  = m_HitPointCont.begin()->second;
    const RKcalcHitPoint& HP   = m_HitPointCont.HitPointOfLayer(layer);
    const ThreeVector&    gpos = HP.PositionInGlobal();
    ThreeVector lpos = gGeom.Global2LocalPos(layer,gpos);
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
KuramaTrack::GetTrajectoryLocalPosition(Int_t layer,
                                        Double_t& x, Double_t& y) const
{
  try {
    const RKcalcHitPoint& HP   = m_HitPointCont.HitPointOfLayer(layer);
    const ThreeVector&    gpos = HP.PositionInGlobal();
    ThreeVector lpos = gGeom.Global2LocalPos(layer,gpos);
    x = lpos.x();
    y = lpos.y();
    return true;
  }
  catch(const std::out_of_range&) {
    return false;
  }
}

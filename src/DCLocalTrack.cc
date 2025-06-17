// -*- C++ -*-

#include "DCLocalTrack.hh"

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <stdexcept>
#include <sstream>

#include <std_ostream.hh>

#include "DCAnalyzer.hh"
#include "DCLTrackHit.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "PrintHelper.hh"
#include "HodoParamMan.hh"
#include "UserParamMan.hh"

namespace
{
const auto& gUser = UserParamMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gHodo = HodoParamMan::GetInstance();
const Double_t& zK18tgt = gGeom.LocalZ("K18Target");
const Double_t& zTgt    = gGeom.LocalZ("Target");
const Int_t ReservedNumOfHits  = 16;
const Int_t DCLocalMinNHits    =  4;
const Int_t DCLocalMinNHitsVXU =  2;// for SSD
const Int_t MaxIteration       = 100;// for Honeycomb
const Double_t MaxChisqrDiff   = 1.0e-3;
}

//_____________________________________________________________________________
DCLocalTrack::DCLocalTrack()
  : m_is_fitted(false),
    m_is_calculated(false),
    m_is_fitted_exclusive(false),
    m_is_bcsdc(false),
    m_Ax(0.),
    m_Ay(0.),
    m_Au(0.),
    m_Av(0.),
    m_Chix(0.),
    m_Chiy(0.),
    m_Chiu(0.),
    m_Chiv(0.),
    m_a(0.),
    m_b(0.),
    m_de(0.),
    m_x0(0.),
    m_y0(0.),
    m_u0(0.),
    m_v0(0.),
    m_chisqr(1e+10),
    m_good_for_tracking(true),
    m_chisqr1st(1e+10),
    m_n_iteration(0),
    m_x0_exclusive(),
    m_y0_exclusive(),
    m_u0_exclusive(),
    m_v0_exclusive(),
    m_de_exclusive(),
    m_chisqr_exclusive(),
    m_chisqr1st_exclusive(),
    m_n_iteration_exclusive()
{
  m_hit_array.reserve(ReservedNumOfHits);
  m_hit_arrayUV.reserve(ReservedNumOfHits);
  debug::ObjectCounter::increase(ClassName());

}

//_____________________________________________________________________________
DCLocalTrack::~DCLocalTrack()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
DCLocalTrack::AddHit(DCLTrackHit *hitp)
{
  if(hitp) m_hit_array.push_back(hitp);
}
//_____________________________________________________________________________
void
DCLocalTrack::AddHitUV(DCLTrackHit *hitp)
{
  if(hitp) m_hit_arrayUV.push_back(hitp);
}

//_____________________________________________________________________________
void
DCLocalTrack::Calculate()
{
  if(m_is_calculated){
    hddaq::cerr << FUNC_NAME << " "
		<< "already called" << std::endl;
    return;
  }

  for(auto& hit: m_hit_array){
    Int_t lnum = hit->GetLayer();
    Double_t z0 = hit->GetZ();
    if(m_is_bcsdc && lnum >= 1 && lnum <= 10){ // SdcIn
      z0 += zK18tgt - zTgt;
    }
    hit->SetCalPosition(GetX(z0), GetY(z0));
    hit->SetCalUV(m_u0, m_v0);
    if(hit->IsHoneycomb()){
      Double_t scal = hit->GetLocalCalPos();
      Double_t wp   = hit->GetWirePosition();
      Double_t dl   = hit->GetDriftLength();
      Double_t ss   = scal-wp>0 ? wp+dl : wp-dl;
      hit->SetLocalHitPos(ss);
    }
  }
  m_is_calculated = true;
}

//_____________________________________________________________________________
void
DCLocalTrack::CalculateExclusive()
{

  if(!m_is_calculated){
    hddaq::cerr << FUNC_NAME << " "
		<< "No inclusive calculatin result" << std::endl;
    return;
  }
  if(!m_is_fitted_exclusive){
    hddaq::cerr << "#W " << FUNC_NAME
		<< "No exclusive fitting" << std::endl;
    return;
  }

  const std::size_t n = m_hit_array.size();
  for(std::size_t i=0; i<n; ++i){
    DCLTrackHit *hit = m_hit_array[i];
    Int_t lnum = hit->GetLayer();
    Double_t z0 = hit->GetZ();
    if(m_is_bcsdc && lnum >= 1 && lnum <= 10){ // SdcIn
      z0 += zK18tgt - zTgt;
    }
    Double_t x_exclusive = m_x0_exclusive[i] + m_u0_exclusive[i]*z0;
    Double_t y_exclusive = m_y0_exclusive[i] + m_v0_exclusive[i]*z0;
    hit -> SetExclusiveReady(true);
    hit -> SetCalPositionExclusive(x_exclusive, y_exclusive);
    hit -> SetCalUVExclusive(m_u0_exclusive[i], m_v0_exclusive[i]);
  }

}

//_____________________________________________________________________________
Int_t
DCLocalTrack::GetNDF() const
{
  const Int_t n = m_hit_array.size();
  Int_t ndf = 0;
  for(Int_t i=0; i<n; ++i){
    if(m_hit_array[i]) ++ndf;
  }
  return ndf-4;
}

//_____________________________________________________________________________
Int_t
DCLocalTrack::GetNHitSFT() const
{
  Int_t n_sft=0;
  for(const auto& hit : m_hit_array){
    if(hit->GetLayer() > 6) ++n_sft;
  }

  return n_sft;
}

//_____________________________________________________________________________
Int_t
DCLocalTrack::GetNHitY() const
{
  Int_t n_y=0;
  for(const auto& hit : m_hit_array){
    if(hit->GetTiltAngle() > 1) ++n_y;
  }

  return n_y;
}

//_____________________________________________________________________________
DCLTrackHit*
DCLocalTrack::GetHit(Int_t nth) const
{
  if(nth<m_hit_array.size())
    return m_hit_array[nth];
  else
    return 0;
}

//_____________________________________________________________________________
DCLTrackHit*
DCLocalTrack::GetHitUV(Int_t nth) const
{
  if(nth<m_hit_arrayUV.size())
    return m_hit_arrayUV[nth];
  else
    return 0;
}

//_____________________________________________________________________________
DCLTrackHit*
DCLocalTrack::GetHitOfLayerNumber(Int_t lnum) const
{
  const Int_t n = m_hit_array.size();
  for(Int_t i=0; i<n; ++i){
    if(m_hit_array[i]->GetLayer()==lnum)
      return m_hit_array[i];
  }
  return 0;
}

//_____________________________________________________________________________
void
DCLocalTrack::DeleteNullHit()
{
  for(Int_t i=0; i<m_hit_array.size(); ++i){
    DCLTrackHit *hitp = m_hit_array[i];
    if(!hitp){
      hddaq::cout << FUNC_NAME << " "
		  << "null hit is deleted" << std::endl;
      m_hit_array.erase(m_hit_array.begin()+i);
      --i;
    }
  }
}

//_____________________________________________________________________________
Bool_t
DCLocalTrack::DoFit()
{
  if(m_is_fitted){
    hddaq::cerr << FUNC_NAME << " "
		<< "already called" << std::endl;
    return false;
  }

  DeleteNullHit();
  const Int_t n = m_hit_array.size();
  if(n < DCLocalMinNHits){
    hddaq::cout << FUNC_NAME << " "
		<< "the number of layers is too small : " << n << std::endl;
    return false;
  }

  const Int_t nItr = HasHoneycomb() ? MaxIteration : 1;
  Double_t prev_chisqr = m_chisqr;
  std::vector <Double_t> z0(n), z(n), wp(n),
    w(n), s(n), ct(n), st(n), coss(n);
  std::vector<Bool_t> honeycomb(n);
  for(Int_t iItr=0; iItr<nItr; ++iItr){
    for(Int_t i=0; i<n; ++i){
      DCLTrackHit *hitp = m_hit_array[i];
      Int_t    lnum = hitp->GetLayer();
      honeycomb[i] = hitp->IsHoneycomb();
      wp[i] = hitp->GetWirePosition();
      z0[i] = hitp->GetZ();
      Double_t ww = gGeom.GetResolution(lnum);
      w[i] = 1./(ww*ww);
      Double_t aa = hitp->GetTiltAngle()*TMath::DegToRad();
      ct[i] = TMath::Cos(aa); st[i] = TMath::Sin(aa);
      Double_t ss = hitp->GetLocalHitPos();
      Double_t dl = hitp->GetDriftLength();
      Double_t dsdz = m_u0*TMath::Cos(aa)+m_v0*TMath::Sin(aa);
      Double_t dcos = TMath::Cos(TMath::ATan(dsdz));
      coss[i] = dcos;
      Double_t dsin = TMath::Sin(TMath::ATan(dsdz));
      Double_t ds = dl * dcos;
      Double_t dz = dl * dsin;
      Double_t scal = iItr==0 ? ss : GetS(z[i],aa);
      if(honeycomb[i]){
	s[i] = scal-wp[i]>0 ? wp[i]+ds : wp[i]-ds;
	z[i] = scal-wp[i]>0 ? z0[i]-dz : z0[i]+dz;
      }else{
	s[i] = ss;
	z[i] = z0[i];
      }
    }

    Double_t x0, u0, y0, v0;
    if(!MathTools::SolveGaussJordan(z, w, s, ct, st,
                                    x0, u0, y0, v0)){
      hddaq::cerr << FUNC_NAME << " Fitting failed" << std::endl;
      return false;
    }

    Double_t chisqr = 0.;
    Double_t de     = 0.;
    for(Int_t i=0; i<n; ++i){
      Double_t scal = (x0+u0*z0[i])*ct[i]+(y0+v0*z0[i])*st[i];
      Double_t ss   = wp[i]+(s[i]-wp[i])/coss[i];
      Double_t res  = honeycomb[i] ? (ss-scal)*coss[i] : s[i]-scal;
      chisqr += w[i]*res*res;
    }
    chisqr /= GetNDF();

    if(iItr==0) m_chisqr1st = chisqr;

    // if worse, not update
    if(prev_chisqr-chisqr>0.){
      m_x0 = x0;
      m_y0 = y0;
      m_u0 = u0;
      m_v0 = v0;
      m_chisqr = chisqr;
      m_de     = de;
    }

    // judge convergence
    if(prev_chisqr-chisqr<MaxChisqrDiff){
#if 0
      // if(chisqr<200. && GetTheta()>4.)
      if(chisqr<20.)
      {
	if(iItr==0) hddaq::cout << "=============" << std::endl;
	hddaq::cout << FUNC_NAME << " NIteration : " << iItr << " "
		    << "chisqr = " << std::setw(10) << std::setprecision(4)
		    << m_chisqr << " "
		    << "diff = " << std::setw(20) << std::left
		    << m_chisqr-m_chisqr1st << " ndf = " << GetNDF() << std::endl;
      }
#endif
      m_n_iteration = iItr;
      break;
    }

    prev_chisqr = chisqr;
  }

  m_is_fitted = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
DCLocalTrack::DoFitExclusive()
{
  if(!m_is_fitted || !m_is_calculated){
    hddaq::cerr << FUNC_NAME << " "
		<< "No inclusive fitting result" << std::endl;
    return false;
  }

  const Int_t n = m_hit_array.size();
  const Int_t nItr = HasHoneycomb() ? MaxIteration : 1;

  m_x0_exclusive.resize(n);
  m_y0_exclusive.resize(n);
  m_u0_exclusive.resize(n);
  m_v0_exclusive.resize(n);
  m_de_exclusive.resize(n);
  m_chisqr_exclusive.resize(n);
  m_chisqr1st_exclusive.resize(n);
  m_n_iteration_exclusive.resize(n);
  for(Int_t iHit=0; iHit<n; ++iHit){
    m_x0_exclusive[iHit] = 0.;
    m_y0_exclusive[iHit] = 0.;
    m_u0_exclusive[iHit] = 0.;
    m_v0_exclusive[iHit] = 0.;
    m_de_exclusive[iHit] = 0.;
    m_chisqr_exclusive[iHit] = 1e+10;
    m_chisqr1st_exclusive[iHit] = 1e+10;
    m_n_iteration_exclusive[iHit] = 0;
    Double_t prev_chisqr = m_chisqr_exclusive[iHit];
    std::vector <Double_t> z0(n-1), z(n-1), wp(n-1),
      w(n-1), s(n-1), ct(n-1), st(n-1), coss(n-1);
    std::vector<Bool_t> honeycomb(n-1);
    for(Int_t iItr=0; iItr<nItr; ++iItr){
      Int_t j=0;
      for(Int_t i=0; i<n; ++i){
	if(i==iHit) continue; //exclude the ith hit
	DCLTrackHit *hitp = m_hit_array[i];
	Int_t    lnum = hitp->GetLayer();

	honeycomb[j] = hitp->IsHoneycomb();
	wp[j] = hitp->GetWirePosition();
	z0[j] = hitp->GetZ();
	Double_t ww = gGeom.GetResolution(lnum);
	w[j] = 1./(ww*ww);
	Double_t aa = hitp->GetTiltAngle()*TMath::DegToRad();
	ct[j] = TMath::Cos(aa); st[j] = TMath::Sin(aa);
	Double_t ss = hitp->GetLocalHitPos();
	Double_t dl = hitp->GetDriftLength();
	Double_t dsdz = m_u0_exclusive[iHit]*TMath::Cos(aa)+m_v0_exclusive[iHit]*TMath::Sin(aa);
	Double_t dcos = TMath::Cos(TMath::ATan(dsdz));
	coss[j] = dcos;
	Double_t dsin = TMath::Sin(TMath::ATan(dsdz));
	Double_t ds = dl * dcos;
	Double_t dz = dl * dsin;
	Double_t x_exclusive = iItr==0 ? 0. : m_x0_exclusive[iHit] + m_u0_exclusive[iHit]*z[j];
	Double_t y_exclusive = iItr==0 ? 0. : m_y0_exclusive[iHit] + m_v0_exclusive[iHit]*z[j];
	Double_t scal_exclusive = iItr==0 ? 0. : x_exclusive*TMath::Cos(aa) + y_exclusive*TMath::Sin(aa);
	Double_t scal = iItr==0 ? ss : scal_exclusive;
	if(honeycomb[j]){
	  s[j] = scal-wp[j]>0 ? wp[j]+ds : wp[j]-ds;
	  z[j] = scal-wp[j]>0 ? z0[j]-dz : z0[j]+dz;
	}else{
	  s[j] = ss;
	  z[j] = z0[j];
	}
	j++;
      } //i

      Double_t x0, u0, y0, v0;
      if(!MathTools::SolveGaussJordan(z, w, s, ct, st,
				      x0, u0, y0, v0)){
	hddaq::cerr << FUNC_NAME << " Fitting failed" << std::endl;
	return false;
      }

      Double_t chisqr = 0.;
      Double_t de     = 0.;
      for(Int_t i=0; i<n-1; ++i){
	Double_t scal = (x0+u0*z0[i])*ct[i]+(y0+v0*z0[i])*st[i];
	Double_t ss   = wp[i]+(s[i]-wp[i])/coss[i];
	Double_t res  = honeycomb[i] ? (ss-scal)*coss[i] : s[i]-scal;
	chisqr += w[i]*res*res;
      }
      chisqr /= (GetNDF()-1);

      if(iItr==0) m_chisqr1st_exclusive[iHit] = chisqr;

      // if worse, not update
      if(prev_chisqr-chisqr>0.){
	m_x0_exclusive[iHit] = x0;
	m_y0_exclusive[iHit] = y0;
	m_u0_exclusive[iHit] = u0;
	m_v0_exclusive[iHit] = v0;
	m_chisqr_exclusive[iHit] = chisqr;
	m_de_exclusive[iHit]     = de;
      }

      // judge convergence
      if(prev_chisqr-chisqr<MaxChisqrDiff){
#if 0
	// if(chisqr<200. && GetTheta()>4.)
	if(chisqr<20.)
	  {
	    if(iItr==0) hddaq::cout << "=============" << std::endl;
	    hddaq::cout << FUNC_NAME << " NIteration : " << iItr << " "
			<< "chisqr = " << std::setw(10) << std::setprecision(4)
			<< m_chisqr_exclusive[iHit] << " "
			<< "diff = " << std::setw(20) << std::left
			<< m_chisqr_exclusive[iHit]-m_chisqr1st_exclusive[iHit] << " ndf = " << GetNDF()-1 << std::endl;
	  }
#endif
	m_n_iteration_exclusive[iHit] = iItr;
	break;
      }

      prev_chisqr = chisqr;
    }//iItr
  } //iHit

#if 0
  hddaq::cout<<"nhits : "<<n<<std::endl;
  hddaq::cout<<"inclusive params x0 "<<m_x0<<" y0 "<<m_y0<<" u0 "<<m_u0<<" v0 "<<m_v0<<std::endl;
  std::cout<<"de "<<m_de<<" chisqr "<<m_chisqr<<" nIteration "<<m_n_iteration<<std::endl;
  for(Int_t i=0; i<n; ++i){
    std::cout<<"exclusive params x0 "<<m_x0_exclusive[i]<<" y0 "<<m_y0_exclusive[i]<<" u0 "<<m_u0_exclusive[i]<<" v0 "<<m_v0_exclusive[i]<<std::endl;
    std::cout<<"de "<<m_de_exclusive[i]<<" chisqr "<<m_chisqr_exclusive[i]<<" nIteration "<<m_n_iteration_exclusive[i]<<std::endl;
  }
#endif

  m_is_fitted_exclusive = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
DCLocalTrack::DoFitBcSdc()
{
  if(m_is_fitted){
    hddaq::cerr << FUNC_NAME << " "
		<< "already called" << std::endl;
    return false;
  }

  m_is_bcsdc = true;
  DeleteNullHit();

  const auto gposBcOut = gGeom.GetGlobalPosition("BC3-X1");
  const auto gposSdcIn = gGeom.GetGlobalPosition("SDC1-U2");
  Double_t BCSdcZoffset = zK18tgt - zTgt;
  Double_t BCSdcXoffset = gposSdcIn.x() - gposBcOut.x();
  Double_t BCSdcYoffset = gposSdcIn.y() - gposBcOut.y();
  Double_t BCSdcUoffset = -0.000205984; //Measured
  Double_t BCSdcVoffset = 0.000691015; //Measured

  Int_t n = m_hit_array.size();
  if(n < DCLocalMinNHits) return false;

  const Int_t nItr = HasHoneycomb() ? MaxIteration : 1;
  Double_t prev_chisqr = m_chisqr;
  std::vector <Double_t> z0(n), z(n), wp(n),
    w(n), s(n), ct(n), st(n), coss(n), bcsdcoffset(n);
  std::vector<Bool_t> honeycomb(n);
  for(Int_t iItr=0; iItr<nItr; ++iItr){
    for(Int_t i=0; i<n; ++i){
      DCLTrackHit *hitp = m_hit_array[i];
      Int_t lnum = hitp->GetLayer();
      honeycomb[i] = hitp->IsHoneycomb();
      wp[i] = hitp->GetWirePosition();
      z0[i] = hitp->GetZ();
      Double_t ww = gGeom.GetResolution(lnum);
      w[i] = 1./(ww*ww);
      Double_t aa = hitp->GetTiltAngle()*TMath::DegToRad();
      ct[i] = TMath::Cos(aa); st[i] = TMath::Sin(aa);
      Double_t ss = hitp->GetLocalHitPos();
      Double_t dl = hitp->GetDriftLength();
      Double_t dsdz = m_u0*TMath::Cos(aa)+m_v0*TMath::Sin(aa);
      Double_t dcos = TMath::Cos(TMath::ATan(dsdz));
      coss[i] = dcos;
      Double_t dsin = TMath::Sin(TMath::ATan(dsdz));
      Double_t ds = dl * dcos;
      Double_t dz = dl * dsin;
      bcsdcoffset[i] = 0.;
      if(lnum >= 1 && lnum <= 10){ // SdcIn
	z0[i] += BCSdcZoffset;
	//bcsdcoffset[i] = BCSdcXoffset*ct[i] + BCSdcYoffset*st[i];
	bcsdcoffset[i] = (BCSdcXoffset - BCSdcUoffset*(z0[i] - zK18tgt))*ct[i] + (BCSdcYoffset - BCSdcVoffset*(z0[i] - zK18tgt))*st[i];
	wp[i] += bcsdcoffset[i];
	ss += bcsdcoffset[i];
      }
      Double_t scal = iItr==0 ? ss : GetS(z[i],aa);
      if(honeycomb[i]){
	s[i] = scal-wp[i]>0 ? wp[i]+ds : wp[i]-ds;
	z[i] = scal-wp[i]>0 ? z0[i]-dz : z0[i]+dz;
      }else{
	s[i] = ss;
	z[i] = z0[i];
      }
      /*
      if(lnum > 10){ // SdcIn
	std::cout<<"  iItr "<<iItr<<" i "<<i<<" scal "<<scal;
	std::cout<<" s "<<s[i]<<" z "<<z[i]<<" wp "<<wp[i]<<" ss "<<ss<<std::endl;
      }
      else{
	std::cout<<"iItr "<<iItr<<" i "<<i<<" scal "<<scal;
	std::cout<<" s "<<s[i]<<" z "<<z[i]<<" wp "<<wp[i]<<" ss "<<ss<<std::endl;
      }
      */
    }

    Double_t x0, u0, y0, v0;
    if(!MathTools::SolveGaussJordan(z, w, s, ct, st,
                                    x0, u0, y0, v0)){
      hddaq::cerr << FUNC_NAME << " Fitting fails" << std::endl;
      return false;
    }

    Double_t chisqr = 0.;
    Double_t de     = 0.;
    for(Int_t i=0; i<n; ++i){
      Double_t scal = (x0+u0*z0[i])*ct[i]+(y0+v0*z0[i])*st[i];
      Double_t ss   = wp[i]+(s[i]-wp[i])/coss[i];
      Double_t res  = honeycomb[i] ? (ss-scal)*coss[i] : s[i]-scal;
      chisqr += w[i]*res*res;
      //if(!honeycomb[i]) std::cout<<"  scal "<<scal<<" s "<<s[i]<<" z0 "<<z0[i]<<" wp "<<wp[i]<<" ss "<<ss<<" resi "<<res<<" coss "<<coss[i]<<std::endl;
      //if(honeycomb[i]) std::cout<<"scal "<<scal<<" s "<<s[i]<<" z0 "<<z0[i]<<" wp "<<wp[i]<<" ss "<<ss<<" resi "<<res<<std::endl;
      if(honeycomb[i]) continue;
    }
    chisqr /= GetNDF();
    std::cout<<std::endl;
    if(iItr==0) m_chisqr1st = chisqr;

    // if worse, not update
    if(prev_chisqr-chisqr>0.){
      m_x0 = x0;
      m_y0 = y0;
      m_u0 = u0;
      m_v0 = v0;
      m_chisqr = chisqr;
      m_de     = de;
    }

    // judge convergence
    if(prev_chisqr-chisqr<MaxChisqrDiff){
      m_n_iteration = iItr;
      break;
    }
    prev_chisqr = chisqr;
  }

  m_is_fitted = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
DCLocalTrack::DoFitVXU()
{
  if(m_is_fitted){
    hddaq::cerr << FUNC_NAME << " "
		<< "already called" << std::endl;
    return false;
  }

  DeleteNullHit();

  const Int_t n = m_hit_array.size();
  if(n<DCLocalMinNHitsVXU) return false;

  Double_t w[n+1],z[n+1],x[n+1];

  for(Int_t i=0; i<n; ++i){
    DCLTrackHit *hitp = m_hit_array[i];
    Int_t lnum = hitp->GetLayer();
    w[i] = gGeom.GetResolution(lnum);
    z[i] = hitp->GetZ();
    x[i] = hitp->GetLocalHitPos();
#if 0
    hddaq::cout << "" << std::endl;
    hddaq::cout << "**********" << std::endl;
    hddaq::cout << std::setw(10) << "layer = " << lnum
		<< std::setw(10) << "wire  = " << hitp->GetWire() << " "
		<< std::setw(20) << "WirePosition = "<<hitp->GetWirePosition() << " "
		<< std::setw(20) << "DriftLength = "<<hitp->GetDriftLength() << " "
		<< std::setw(20) << "hit position = "<<hitp->GetLocalHitPos()<< " "
		<< std::endl;
    hddaq::cout << "**********" << std::endl;
    hddaq::cout << "" << std::endl;
#endif
  }

  Double_t A=0, B=0, C=0, D=0, E=0;// <-Add!!
  for(Int_t i=0; i<n; ++i){
    A += z[i]/(w[i]*w[i]);
    B += 1/(w[i]*w[i]);
    C += x[i]/(w[i]*w[i]);
    D += z[i]*z[i]/(w[i]*w[i]);
    E += x[i]*z[i]/(w[i]*w[i]);
  }

  m_a = (E*B-C*A)/(D*B-A*A);
  m_b = (D*C-E*A)/(D*B-A*A);

  Double_t chisqr = 0.;
  for(Int_t i=0; i<n; ++i){
    chisqr += (x[i]-m_a*z[i]-m_b)*(x[i]-m_a*z[i]-m_b)/(w[i]*w[i]);
  }

  if(n==2) chisqr  = 0.;
  else     chisqr /= n-2.;
  m_chisqr = chisqr;

  for(Int_t i=0; i<n; ++i){
    DCLTrackHit *hitp = m_hit_array[i];
    if(hitp){
      Double_t zz = hitp->GetZ();
      hitp->SetLocalCalPosVXU(m_a*zz+m_b);
    }
  }

  m_is_fitted = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
DCLocalTrack::FindLayer(Int_t layer) const
{
  const Int_t n = m_hit_array.size();
  for(Int_t i=0; i<n; ++i){
    DCLTrackHit *hitp = m_hit_array[i];
    if(!hitp)
      continue;
    if(layer == hitp->GetLayer())
      return true;
  }
  return false;
}

//_____________________________________________________________________________
Double_t
DCLocalTrack::GetWire(Int_t layer) const
{
  const Int_t n = m_hit_array.size();
  for(Int_t i=0; i<n; ++i){
    DCLTrackHit *hitp = m_hit_array[i];
    if(!hitp) continue;
    if(layer == hitp->GetLayer()) return hitp->GetWire();
  }
  return TMath::QuietNaN();
}

//_____________________________________________________________________________
Double_t
DCLocalTrack::GetDifVXU() const
{
  static const Double_t Cu = TMath::Cos( 15.*TMath::DegToRad());
  static const Double_t Cv = TMath::Cos(-15.*TMath::DegToRad());
  static const Double_t Cx = TMath::Cos(  0.*TMath::DegToRad());

  return
    pow(m_Av/Cv - m_Ax/Cx, 2) +
    pow(m_Ax/Cx - m_Au/Cu, 2) +
    pow(m_Au/Cu - m_Av/Cv, 2);
}

//_____________________________________________________________________________
Double_t
DCLocalTrack::GetDifVXUSDC34() const
{
  static const Double_t Cu = TMath::Cos( 30.*TMath::DegToRad());
  static const Double_t Cv = TMath::Cos(-30.*TMath::DegToRad());
  static const Double_t Cx = TMath::Cos(  0.*TMath::DegToRad());

  return
    pow(m_Av/Cv - m_Ax/Cx, 2) +
    pow(m_Ax/Cx - m_Au/Cu, 2) +
    pow(m_Au/Cu - m_Av/Cv, 2);
}

//_____________________________________________________________________________
Double_t
DCLocalTrack::GetPhi() const
{
  return TMath::ATan2(m_u0, m_v0);
}

//_____________________________________________________________________________
Double_t
DCLocalTrack::GetTheta() const
{
  Double_t cost = 1./TMath::Sqrt(1.+m_u0*m_u0+m_v0*m_v0);
  return TMath::ACos(cost)*TMath::RadToDeg();
}

//_____________________________________________________________________________
Bool_t
DCLocalTrack::HasHoneycomb() const
{
  for(Int_t i=0, n = m_hit_array.size(); i<n; ++i){
    DCLTrackHit *hitp = m_hit_array[i];
    if(!hitp) continue;
    if(hitp->IsHoneycomb()) return true;
  }
  return false;
}

//_____________________________________________________________________________
Bool_t
DCLocalTrack::ReCalc(Bool_t applyRecursively)
{
  Int_t n = m_hit_array.size();
  for(Int_t i=0; i<n; ++i){
    DCLTrackHit *hitp = m_hit_array[i];
    if(hitp) hitp->ReCalc(applyRecursively);
  }

  Bool_t status = DoFit();
  if(!status){
    hddaq::cerr << FUNC_NAME << " "
		<< "Recalculation fails" << std::endl;
  }

  return status;
}

//_____________________________________________________________________________
void
DCLocalTrack::Print(const TString& arg, std::ostream& ost) const
{
  PrintHelper helper(3, std::ios::fixed, ost);

  const Int_t w = 8;
  ost << FUNC_NAME << " " << arg << std::endl
      << " X0 : " << std::setw(w) << std::left << m_x0
      << " Y0 : " << std::setw(w) << std::left << m_y0
      << " U0 : " << std::setw(w) << std::left << m_u0
      << " V0 : " << std::setw(w) << std::left << m_v0;
  // helper.setf(std::ios::scientific);
  ost << " Chisqr : " << std::setw(w) << m_chisqr << std::endl;
  helper.setf(std::ios::fixed);
  const Int_t n = m_hit_array.size();
  for(Int_t i=0; i<n; ++i){
    DCLTrackHit *hitp = m_hit_array[i];
    if(!hitp) continue;
    Int_t lnum = hitp->GetLayer();
    Double_t zz = hitp->GetZ();
    Double_t s  = hitp->GetLocalHitPos();
    Double_t res = hitp->GetResidual();
    Double_t aa = hitp->GetTiltAngle()*TMath::DegToRad();
    // Double_t scal=GetX(zz)*TMath::Cos(aa)+GetY(zz)*TMath::Sin(aa);
    if(m_is_bcsdc && lnum >= 1 && lnum <= 10){
      const auto gposBcOut = gGeom.GetGlobalPosition("BC3-X1");
      const auto gposSdcIn = gGeom.GetGlobalPosition("SDC1-U2");
      Double_t BCSdcXoffset = gposSdcIn.x() - gposBcOut.x();
      Double_t BCSdcYoffset = gposSdcIn.y() - gposBcOut.y();
      res += BCSdcXoffset*TMath::Cos(aa) + BCSdcYoffset*TMath::Sin(aa);
    }
    const TString& h = hitp->IsHoneycomb() ? "+" : "-";
    ost << "[" << std::setw(2) << i << "]"
	<< " #"  << std::setw(2) << lnum << h
	<< " S " << std::setw(w) << s
	<< " (" << std::setw(w) << GetX(zz)
	<< ", "  << std::setw(w) << GetY(zz)
	<< ", "  << std::setw(w) << zz
	<< ")"
	<< " " << std::setw(w) << s
	<< " -> " << std::setw(w) << res << std::endl;
    // << " -> " << std::setw(w) << s-scal << std::endl;
  }
  ost << std::endl;
}

//_____________________________________________________________________________
void
DCLocalTrack::PrintVXU(const TString& arg) const
{
  PrintHelper helper(3, std::ios::fixed);

  const Int_t w = 10;
  hddaq::cout << FUNC_NAME << " " << arg << std::endl
	      << "(Local X = A*z + B) "
	      << " A : " << std::setw(w) << m_a
	      << " B : " << std::setw(w) << m_b;
  helper.setf(std::ios::scientific);
  hddaq::cout << " Chisqr : " << std::setw(w) << m_chisqr << std::endl;
  helper.setf(std::ios::fixed);

  const Int_t n = m_hit_array.size();
  for(Int_t i=0; i<n; ++i){
    const DCLTrackHit * const hitp = m_hit_array[i];
    if(!hitp) continue;
    Int_t    lnum = hitp->GetLayer();
    Double_t zz   = hitp->GetZ();
    Double_t s    = hitp->GetLocalHitPos();
    Double_t res  = hitp->GetResidualVXU();
    const TString& h = hitp->IsHoneycomb() ? "+" : "-";
    hddaq::cout << "[" << std::setw(2) << i << "] "
		<< " #"  << std::setw(2) << lnum << h
		<< " S " << std::setw(w) << s
		<< " (" << std::setw(w) << (m_a*zz+m_b)
		<< ", "  << std::setw(w) << zz
		<< ")"
		<< " " << std::setw(w) << s
		<< " -> " << std::setw(w) << res << std::endl;
  }
  hddaq::cout << std::endl;
}

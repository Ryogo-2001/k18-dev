// -*- C++ -*-

#include "DCAnalyzer.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <TH2D.h>

#include "ConfMan.hh"
#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DCLocalTrack.hh"
#include "DCRawHit.hh"
#include "DCTrackSearch.hh"
#include "DebugCounter.hh"
#include "DebugTimer.hh"
#include "FiberCluster.hh"
#include "FuncName.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "K18Parameters.hh"
//#include "K18TrackU2D.hh"
#include "K18TrackD2U.hh"
#include "KuramaTrack.hh"
#include "MathTools.hh"
#include "MWPCCluster.hh"
#include "RawData.hh"
#include "UserParamMan.hh"
#include "DeleteUtility.hh"

#define DefStatic
#include "DCParameters.hh"
#undef DefStatic

// Tracking routine selection __________________________________________________
/* BcInTracking */
#define UseBcIn    0 // not supported
/* BcOutTracking */
#define BcOut_XUV  0 // XUV Tracking (slow but accerate)
#define BcOut_Pair 1 // Pair plane Tracking (fast but bad for large angle track)
/* SdcInTracking */
#define SdcIn_XUV         0 // XUV Tracking (not used in KURAMA)
#define SdcIn_Pair        1 // Pair plane Tracking (fast but bad for large angle track)
#define SdcIn_Deletion    1 // Deletion method for too many combinations

namespace
{
using namespace K18Parameter;
const auto& gConf   = ConfMan::GetInstance();
const auto& gGeom   = DCGeomMan::GetInstance();
const auto& gUser   = UserParamMan::GetInstance();

//_____________________________________________________________________________
const Double_t& pK18 = ConfMan::Get<Double_t>("PK18");
const Int_t& IdTOFUX = gGeom.DetectorId("TOF-UX");
const Int_t& IdTOFUY = gGeom.DetectorId("TOF-UY");
const Int_t& IdTOFDX = gGeom.DetectorId("TOF-DX");
const Int_t& IdTOFDY = gGeom.DetectorId("TOF-DY");

//const Double_t TimeDiffToYTOF = 77.3511; // [mm/ns]
const Double_t TimeDiffToYTOF = 77.3511; // [mm/ns]
const Double_t MaxChiSqrKuramaTrack = 10000.;
const Double_t MaxTimeDifMWPC       =   100.;

const Double_t kMWPCClusteringWireExtension =  1.0; // [mm]
const Double_t kMWPCClusteringTimeExtension = 10.0; // [nsec]

//_____________________________________________________________________________
inline Bool_t /* for MWPCCluster */
isConnectable(Double_t wire1, Double_t leading1, Double_t trailing1,
              Double_t wire2, Double_t leading2, Double_t trailing2,
              Double_t wExt,  Double_t tExt)
{
  Double_t w1Min = wire1 - wExt;
  Double_t w1Max = wire1 + wExt;
  Double_t t1Min = leading1  - tExt;
  Double_t t1Max = trailing1 + tExt;
  Double_t w2Min = wire2 - wExt;
  Double_t w2Max = wire2 + wExt;
  Double_t t2Min = leading2  - tExt;
  Double_t t2Max = trailing2 + tExt;
  Bool_t isWireOk = !(w1Min>w2Max || w1Max<w2Min);
  Bool_t isTimeOk = !(t1Min>t2Max || t1Max<t2Min);
#if 0
  hddaq::cout << __func__ << std::endl
              << " w1 = " << wire1
              << " le1 = " << leading1
              << " tr1 = " << trailing1 << "\n"
              << " w2 = " << wire2
              << " le2 = " << leading2
              << " tr2 = " << trailing2 << "\n"
              << " w1(" << w1Min << " -- " << w1Max << "), t1("
              << t1Min << " -- " << t1Max << ")\n"
              << " w2(" << w2Min << " -- " << w2Max << "), t2("
              << t2Min << " -- " << t2Max << ")\n"
              << " wire : " << isWireOk
              << ", time : " << isTimeOk
              << std::endl;
#endif
  return (isWireOk && isTimeOk);
}

//_____________________________________________________________________________
inline void
printConnectionFlag(const std::vector<std::deque<Bool_t> >& flag)
{
  for(Int_t i=0, n=flag.size(); i<n; ++i){
    hddaq::cout << std::endl;
    for(Int_t j=0, m=flag[i].size(); j<m; ++j){
      hddaq::cout << " " << flag[i][j];
    }
  }
  hddaq::cout << std::endl;
}
}

//_____________________________________________________________________________
DCAnalyzer::DCAnalyzer()
  : m_max_v0diff(90.),
    m_is_decoded(n_type),
    m_much_combi(n_type),
    m_MWPCClCont(NumOfLayersBcIn+1),
    m_TempBcInHC(NumOfLayersBcIn+1),
    m_BcInHC(NumOfLayersBcIn+1),
    m_BcOutHC(NumOfLayersBcOut+2),
    m_SdcInHC(NumOfLayersSdcIn+1),
    m_SdcOutHC(NumOfLayersSdcOut+1)
{
  for(Int_t i=0; i<n_type; ++i){
    m_is_decoded[i] = false;
    m_much_combi[i] = 0;
  }
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
DCAnalyzer::~DCAnalyzer()
{
  ClearKuramaTracks();
#if UseBcIn
  ClearK18TracksU2D();
  ClearTracksBcIn();
#endif
  ClearK18TracksD2U();
  ClearTracksSdcOut();
  ClearTracksSdcIn();
  ClearTracksBcOut();
  ClearTracksBcOutSdcIn();
  ClearTracksSdcInSdcOut();
  ClearDCHits();
  ClearVtxHits();
  debug::ObjectCounter::decrease(ClassName());
}


//_____________________________________________________________________________
void
DCAnalyzer::PrintKurama(const TString& arg) const
{
  Int_t nn = m_KuramaTC.size();
  hddaq::cout << FUNC_NAME << " " << arg << std::endl
              << "   KuramaTC.size : " << nn << std::endl;
  for(const auto& track: m_KuramaTC){
    hddaq::cout << " Niter=" << std::setw(3) << track->Niteration()
                << " ChiSqr=" << track->ChiSquare()
                << " P=" << track->PrimaryMomentum().Mag()
                << " PL(TOF)=" << track->PathLengthToTOF()
                << std::endl;
  }
}

//_____________________________________________________________________________
#if UseBcIn
Bool_t
DCAnalyzer::DecodeBcInHits(RawData *rawData)
{
  if(m_is_decoded[kBcIn]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  ClearBcInHits();

  for(Int_t layer=1; layer<=NumOfLayersBcIn; ++layer){
    const DCRHitContainer &RHitCont = rawData->GetBcInRawHC(layer);
    Int_t nh = RHitCont.size();
    for(Int_t i=0; i<nh; ++i){
      DCRawHit *rhit  = RHitCont[i];
      DCHit    *thit  = new DCHit(rhit->PlaneId()+PlOffsBc, rhit->WireId());
      Int_t       nhtdc = rhit->GetTdcSize();
      if(!thit) continue;
      for(Int_t j=0; j<nhtdc; ++j){
        thit->SetTdcVal(rhit->GetTdc(j));
        thit->SetTdcTrailing(rhit->GetTrailing(j));
      }

      if(thit->CalcMWPCObservables())
        m_TempBcInHC[layer].push_back(thit);
      else
        delete thit;
    }

    // hddaq::cout<<"*************************************"<<std::endl;
    Int_t ncl = clusterizeMWPCHit(m_TempBcInHC[layer], m_MWPCClCont[layer]);
    // hddaq::cout<<"numCl="<< ncl << std::endl;
    for(Int_t i=0; i<ncl; ++i){
      MWPCCluster *p = m_MWPCClCont[layer][i];
      if(!p) continue;

      const MWPCCluster::Statistics& mean  = p->GetMean();
      const MWPCCluster::Statistics& first = p->GetFirst();
      Double_t mwire    = mean.m_wire;
      Double_t mwirepos = mean.m_wpos;
      Double_t mtime    = mean.m_leading;
      Double_t mtrail   = mean.m_trailing;

      DCHit *hit = new DCHit(layer+PlOffsBc, mwire);
      if(!hit) continue;
      hit->SetClusterSize(p->GetClusterSize());
      hit->SetMWPCFlag(true);
      hit->SetWire(mwire);
      hit->SetMeanWire(mwire);
      hit->SetMeanWirePosition(mwirepos);
      hit->SetTrailing(mtrail);
      hit->SetDummyPair();
      hit->SetTdcVal(0);
      hit->SetTdcTrailing(0);

      if(hit->CalcMWPCObservables())
        m_BcInHC[layer].push_back(hit);
      else
        delete hit;
    }
    // hddaq::cout << "nh="<< m_BcInHC[layer].size() <<std::endl;
  }

  m_is_decoded[kBcIn] = true;
  return true;
}
#endif

//_____________________________________________________________________________
Bool_t
DCAnalyzer::DecodeBcOutHits(RawData *rawData)
{
  if(m_is_decoded[kBcOut]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  ClearBcOutHits();

  for(Int_t layer=1; layer<=NumOfLayersBcOut; ++layer){
    for(const auto& rhit: rawData->GetBcOutRawHC(layer)){
      auto hit = new DCHit(rhit->PlaneId()+PlOffsBc, rhit->WireId());
      Int_t nhtdc = rhit->GetTdcSize();
      Int_t nhtrailing = rhit->GetTrailingSize();
      if(!hit) continue;
      for(Int_t j=0; j<nhtdc; ++j){
        hit->SetTdcVal(rhit->GetTdc(j));
      }
      for(Int_t j=0; j<nhtrailing; ++j){
        hit->SetTdcTrailing(rhit->GetTrailing(j));
      }
      if(hit->CalcDCObservables()){
        m_BcOutHC[layer].push_back(hit);
      }else{
        delete hit;
      }
    }
  }

  m_is_decoded[kBcOut] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::DecodeSdcInHits(RawData *rawData)
{
  if(m_is_decoded[kSdcIn]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  ClearSdcInHits();

  for(Int_t layer=1; layer<=NumOfLayersSdcIn; ++layer){
    for(const auto& rhit: rawData->GetSdcInRawHC(layer)){
      auto hit = new DCHit(rhit->PlaneId(), rhit->WireId());
      Int_t nhtdc = rhit->GetTdcSize();
      Int_t nhtrailing = rhit->GetTrailingSize();
      if(!hit) continue;
      for(Int_t j=0; j<nhtdc; ++j){
        hit->SetTdcVal(rhit->GetTdc(j));
      }
      for(Int_t j=0; j<nhtrailing; ++j){
        hit->SetTdcTrailing(rhit->GetTrailing(j));
      }
      if(hit->CalcDCObservables()){
        m_SdcInHC[layer].push_back(hit);
      }else{
        delete hit;
      }
    }
  }

  m_is_decoded[kSdcIn] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::DecodeSdcOutHits(RawData *rawData , Double_t ofs_dt/*=0.*/)
{
  if(m_is_decoded[kSdcOut]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  ClearSdcOutHits();

  for(Int_t layer=1; layer<=NumOfLayersSdcOut; ++layer){
    for(const auto& rhit: rawData->GetSdcOutRawHC(layer)){
      auto hit = new DCHit(rhit->PlaneId(), rhit->WireId());
      Int_t nhtdc = rhit->GetTdcSize();
      Int_t nhtrailing = rhit->GetTrailingSize();
      if(!hit) continue;
      if(ofs_dt != 0.) hit->SetOfsdT(ofs_dt);
      for(Int_t j=0; j<nhtdc; ++j){
        hit->SetTdcVal(rhit->GetTdc(j));
      }
      for(Int_t j=0; j<nhtrailing; ++j){
        hit->SetTdcTrailing(rhit->GetTrailing(j));
      }
      if(hit->CalcDCObservables()){
        m_SdcOutHC[layer].push_back(hit);
      }else{
        delete hit;
      }
    }
  }

  m_is_decoded[kSdcOut] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::DecodeRawHits(RawData *rawData)
{
  ClearDCHits();
#if UseBcIn
  DecodeBcInHits(rawData);
#endif
  DecodeBcOutHits(rawData);
  DecodeSdcInHits(rawData);
  DecodeSdcOutHits(rawData);
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::DecodeTOFHits(const Hodo2HitContainer& HitCont)
{
  if(m_is_decoded[kTOF]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  ClearTOFHits();

  // for the tilting plane case
  static const Double_t RA2 = gGeom.GetRotAngle2("TOF");
  static const TVector3 TOFPos[2] = {
    gGeom.GetGlobalPosition("TOF-UX"),
    gGeom.GetGlobalPosition("TOF-DX") };

  for(const auto& hodo_hit: HitCont){
    const Double_t seg = hodo_hit->SegmentId()+1;
    const Double_t dt  = hodo_hit->TimeDiff();
    Int_t layer_x = -1;
    Int_t layer_y = -1;
    if((Int_t)seg%2==0){
      layer_x  = IdTOFUX;
      layer_y  = IdTOFUY;
    }
    if((Int_t)seg%2==1){
      layer_x  = IdTOFDX;
      layer_y  = IdTOFDY;
    }
    Double_t wpos = gGeom.CalcWirePosition(layer_x, seg);
    TVector3 w(wpos, 0., 0.);
    w.RotateY(RA2*TMath::DegToRad()); // for the tilting plane case
    const TVector3 hit_pos = TOFPos[(Int_t)seg%2] + w
      + TVector3(0., dt*TimeDiffToYTOF, 0.);
    // X
    DCHit *dc_hit_x = new DCHit(layer_x, seg);
    dc_hit_x->SetWirePosition(hit_pos.x());
    dc_hit_x->SetZ(hit_pos.z());
    dc_hit_x->SetTiltAngle(0.);
    dc_hit_x->SetDummyPair();
    m_TOFHC.push_back(dc_hit_x);
    // Y
    DCHit *dc_hit_y = new DCHit(layer_y, seg);
    dc_hit_y->SetWirePosition(hit_pos.y());
    dc_hit_y->SetZ(hit_pos.z());
    dc_hit_y->SetTiltAngle(90.);
    dc_hit_y->SetDummyPair();
    m_TOFHC.push_back(dc_hit_y);
  }

  m_is_decoded[kTOF] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::DecodeTOFHits(const HodoClusterContainer& ClCont)
{
  if(m_is_decoded[kTOF]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  ClearTOFHits();

  static const Double_t RA2 = gGeom.GetRotAngle2("TOF");
  static const TVector3 TOFPos[2] = {
    gGeom.GetGlobalPosition("TOF-UX"),
    gGeom.GetGlobalPosition("TOF-DX") };

  for(const auto& hodo_cluster: ClCont){
    const Double_t seg = hodo_cluster->MeanSeg()+1;
    const Double_t dt  = hodo_cluster->TimeDif();
    Int_t layer_x = -1;
    Int_t layer_y = -1;
    if((Int_t)seg%2==0){
      layer_x  = IdTOFUX;
      layer_y  = IdTOFUY;
    }
    if((Int_t)seg%2==1){
      layer_x  = IdTOFDX;
      layer_y  = IdTOFDY;
    }
    Double_t wpos = gGeom.CalcWirePosition(layer_x, seg);
    TVector3 w(wpos, 0., 0.);
    w.RotateY(RA2*TMath::DegToRad());
    const TVector3& hit_pos = TOFPos[(Int_t)seg%2] + w
      + TVector3(0., dt*TimeDiffToYTOF, 0.);
    // X
    DCHit *dc_hit_x = new DCHit(layer_x, seg);
    dc_hit_x->SetWirePosition(hit_pos.x());
    dc_hit_x->SetZ(hit_pos.z());
    dc_hit_x->SetTiltAngle(0.);
    dc_hit_x->SetDummyPair();
    m_TOFHC.push_back(dc_hit_x);
    // Y
    DCHit *dc_hit_y = new DCHit(layer_y, seg);
    dc_hit_y->SetWirePosition(hit_pos.y());
    dc_hit_y->SetZ(hit_pos.z());
    dc_hit_y->SetTiltAngle(90.);
    dc_hit_y->SetDummyPair();
    m_TOFHC.push_back(dc_hit_y);
  }

  m_is_decoded[kTOF] = true;
  return true;
}

//_____________________________________________________________________________
#if UseBcIn
Bool_t
DCAnalyzer::TrackSearchBcIn()
{
  track::MWPCLocalTrackSearch(&(m_BcInHC[1]), m_BcInTC);
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchBcIn(const std::vector<std::vector<DCHitContainer> >& hc)
{
  track::MWPCLocalTrackSearch(hc, m_BcInTC);
  return true;
}
#endif

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchBcOut(Double_t T0Seg, Bool_t Exclusive)
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerBcOut");

#if BcOut_Pair //Pair Plane Tracking Routine for BcOut
  Int_t ntrack = track::LocalTrackSearch(m_BcOutHC, PPInfoBcOut, NPPInfoBcOut,
                                         m_BcOutTC, Exclusive, MinLayer, T0Seg);
  return ntrack == -1 ? false : true;
#endif

#if BcOut_XUV  //XUV Tracking Routine for BcOut
  Int_t ntrack = track::LocalTrackSearchVUX(m_BcOutHC, PPInfoBcOut, NPPInfoBcOut,
                                            m_BcOutTC, MinLayer);
  return ntrack == -1 ? false : true;
#endif

  return false;
}

//_____________________________________________________________________________
// Use with BH2Filter
Bool_t
DCAnalyzer::TrackSearchBcOut(const std::vector<std::vector<DCHitContainer> >& hc, Double_t T0Seg, Bool_t Exclusive)
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerBcOut");

#if BcOut_Pair //Pair Plane Tracking Routine for BcOut
  Int_t ntrack = track::LocalTrackSearch(hc, PPInfoBcOut, NPPInfoBcOut,
					 m_BcOutTC, Exclusive, MinLayer, T0Seg);
  return ntrack == -1 ? false : true;
#endif

#if BcOut_XUV  //XUV Tracking Routine for BcOut
  Int_t ntrack = track::LocalTrackSearchVUX(hc, PPInfoBcOut, NPPInfoBcOut,
					    m_BcOutTC, MinLayer);
  return ntrack == -1 ? false : true;
#endif

  return false;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchSdcIn(Bool_t Exclusive)
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerSdcIn");

  track::LocalTrackSearch(m_SdcInHC, PPInfoSdcIn, NPPInfoSdcIn,
			  m_SdcInTC, Exclusive, MinLayer);

  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchSdcOut(Bool_t Exclusive)
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerSdcOut");

  track::LocalTrackSearchSdcOut(m_SdcOutHC, PPInfoSdcOut, NPPInfoSdcOut,
				m_SdcOutTC, Exclusive, MinLayer);

  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchSdcOut(const Hodo2HitContainer& TOFCont, Bool_t Exclusive)
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerSdcOut");

  if(!DecodeTOFHits(TOFCont)) return false;

#if 0
  for(std::size_t i=0, n=m_TOFHC.size(); i<n; ++i){
    m_TOFHC[i]->Print();
  }
#endif

  track::LocalTrackSearchSdcOut(m_TOFHC, m_SdcOutHC, PPInfoSdcOut, NPPInfoSdcOut+2,
                                m_SdcOutTC, Exclusive, MinLayer);

  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchSdcOut(const HodoClusterContainer& TOFCont, Bool_t Exclusive )
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerSdcOut");

  if(!DecodeTOFHits(TOFCont)) return false;

  track::LocalTrackSearchSdcOut(m_TOFHC, m_SdcOutHC, PPInfoSdcOut, NPPInfoSdcOut+2,
				m_SdcOutTC, Exclusive, MinLayer);

  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchBcOutSdcIn(Bool_t Exclusive)
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerBcOutSdcIn");

  track::LocalTrackSearchBcOutSdcIn(m_BcOutHC, PPInfoBcOut, m_SdcInHC,
				    PPInfoSdcIn, NPPInfoBcOut, NPPInfoSdcIn,
				    m_BcOutSdcInTC, Exclusive, MinLayer);

  return true;
}

//_____________________________________________________________________________
#if UseBcIn
Bool_t
DCAnalyzer::TrackSearchK18U2D()
{
  ClearK18TracksU2D();

  Int_t nIn  = m_BcInTC.size();
  Int_t nOut = m_BcOutTC.size();

# if 0
  hddaq::cout << "**************************************" << std::endl;
  hddaq::cout << FUNC_NAME << ": #TracksIn=" << std::setw(3) << nIn
              << " #TracksOut=" << std::setw(3) << nOut << std::endl;
# endif

  if(nIn==0 || nOut==0) return true;

  for(Int_t iIn=0; iIn<nIn; ++iIn){
    DCLocalTrack *trIn = m_BcInTC[iIn];
# if 0
    hddaq::cout << "TrackIn  :" << std::setw(2) << iIn
                << " X0=" << trIn->GetX0() << " Y0=" << trIn->GetY0()
                << " U0=" << trIn->GetU0() << " V0=" << trIn->GetV0()
                << std::endl;
# endif
    if(!trIn->GoodForTracking() ||
       trIn->GetX0()<MinK18InX || trIn->GetX0()>MaxK18InX ||
       trIn->GetY0()<MinK18InY || trIn->GetY0()>MaxK18InY ||
       trIn->GetU0()<MinK18InU || trIn->GetU0()>MaxK18InU ||
       trIn->GetV0()<MinK18InV || trIn->GetV0()>MaxK18InV) continue;
    for(Int_t iOut=0; iOut<nOut; ++iOut){
      DCLocalTrack *trOut=m_BcOutTC[iOut];
# if 0
      hddaq::cout << "TrackOut :" << std::setw(2) << iOut
                  << " X0=" << trOut->GetX0() << " Y0=" << trOut->GetY0()
                  << " U0=" << trOut->GetU0() << " V0=" << trOut->GetV0()
                  << std::endl;
# endif
      if(!trOut->GoodForTracking() ||
         trOut->GetX0()<MinK18OutX || trOut->GetX0()>MaxK18OutX ||
         trOut->GetY0()<MinK18OutY || trOut->GetY0()>MaxK18OutY ||
         trOut->GetU0()<MinK18OutU || trOut->GetU0()>MaxK18OutU ||
         trOut->GetV0()<MinK18OutV || trOut->GetV0()>MaxK18OutV) continue;

# if 0
      hddaq::cout << FUNC_NAME << ": In -> " << trIn->GetChiSquare()
                  << " (" << std::setw(2) << trIn->GetNHit() << ") "
                  << "Out -> " << trOut->GetChiSquare()
                  << " (" << std::setw(2) << trOut->GetNHit() << ") "
                  << std::endl;
# endif

      K18TrackU2D *track=new K18TrackU2D(trIn, trOut, TMath::Abs(pK18));
      if(track && track->DoFit())
        m_K18U2DTC.push_back(track);
      else
        delete track;
    }
  }

# if 0
  hddaq::cout << "********************" << std::endl;
  {
    Int_t nn = m_K18U2DTC.size();
    hddaq::cout << FUNC_NAME << ": Before sorting. #Track="
                << nn << std::endl;
    for(Int_t i=0; i<nn; ++i){
      K18TrackU2D *tp = m_K18U2DTC[i];
      hddaq::cout << std::setw(3) << i
                  << " ChiSqr=" << tp->chisquare()
                  << " Delta=" << tp->Delta()
                  << " P=" << tp->P() << "\n";
      //      hddaq::cout<<"********************"<<std::endl;
      //      hddaq::cout << "In :"
      //         << " X " << tp->Xin() << "(" << tp->TrackIn()->GetX0() << ")"
      //         << " Y " << tp->Yin() << "(" << tp->TrackIn()->GetY0() << ")"
      //         << " U " << tp->Uin() << "(" << tp->TrackIn()->GetU0() << ")"
      //         << " V " << tp->Vin() << "(" << tp->TrackIn()->GetV0() << ")"
      //         << "\n";
      //      hddaq::cout << "Out:"
      //         << " X " << tp->Xout() << "(" << tp->TrackOut()->GetX0() << ")"
      //         << " Y " << tp->Yout() << "(" << tp->TrackOut()->GetY0() << ")"
      //         << " U " << tp->Uout() << "(" << tp->TrackOut()->GetU0() << ")"
      //         << " V " << tp->Vout() << "(" << tp->TrackOut()->GetV0() << ")"
      //         << std::endl;
    }
  }
# endif

  std::sort(m_K18U2DTC.begin(), m_K18U2DTC.end(), K18TrackU2DComp());

  return true;
}
#endif // UseBcIn

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchK18D2U(const std::vector<Double_t>& XinCont)
{
  ClearK18TracksD2U();

  std::size_t nIn  = XinCont.size();
  std::size_t nOut = m_BcOutTC.size();

  if(nIn==0 || nOut==0)
    return true;

  for(std::size_t iIn=0; iIn<nIn; ++iIn){
    Double_t LocalX = XinCont.at(iIn);
    for(std::size_t iOut=0; iOut<nOut; ++iOut){
      DCLocalTrack *trOut = m_BcOutTC[iOut];
#if 0
      hddaq::cout << "TrackOut :" << std::setw(2) << iOut
                  << " X0=" << trOut->GetX0() << " Y0=" << trOut->GetY0()
                  << " U0=" << trOut->GetU0() << " V0=" << trOut->GetV0()
                  << std::endl;
#endif
      if(!trOut->GoodForTracking() ||
         trOut->GetX0()<MinK18OutX || trOut->GetX0()>MaxK18OutX ||
         trOut->GetY0()<MinK18OutY || trOut->GetY0()>MaxK18OutY ||
         trOut->GetU0()<MinK18OutU || trOut->GetU0()>MaxK18OutU ||
         trOut->GetV0()<MinK18OutV || trOut->GetV0()>MaxK18OutV) continue;

#if 0
      hddaq::cout << FUNC_NAME
                  << "Out -> " << trOut->GetChiSquare()
                  << " (" << std::setw(2) << trOut->GetNHit() << ") "
                  << std::endl;
#endif

      K18TrackD2U *track = new K18TrackD2U(LocalX, trOut, TMath::Abs(pK18));
      if(track && track->CalcMomentumD2U())
	m_K18D2UTC.push_back(track);
      else
	delete track;
    }
  }

#if 0
  hddaq::cout<<"********************"<<std::endl;
  {
    Int_t nn = m_K18D2UTC.size();
    hddaq::cout << FUNC_NAME << ": Before sorting. #Track="
                << nn << std::endl;
    for(Int_t i=0; i<nn; ++i){
      K18TrackD2U *tp = m_K18D2UTC[i];
      hddaq::cout << std::setw(3) << i
                  << " ChiSqr=" << tp->chisquare()
                  << " Delta=" << tp->Delta()
                  << " P=" << tp->P() << "\n";
      //      hddaq::cout<<"********************"<<std::endl;
      //      hddaq::cout << "In :"
      //         << " X " << tp->Xin() << "(" << tp->TrackIn()->GetX0() << ")"
      //         << " Y " << tp->Yin() << "(" << tp->TrackIn()->GetY0() << ")"
      //         << " U " << tp->Uin() << "(" << tp->TrackIn()->GetU0() << ")"
      //         << " V " << tp->Vin() << "(" << tp->TrackIn()->GetV0() << ")"
      //         << "\n";
      //      hddaq::cout << "Out:"
      //         << " X " << tp->Xout() << "(" << tp->TrackOut()->GetX0() << ")"
      //         << " Y " << tp->Yout() << "(" << tp->TrackOut()->GetY0() << ")"
      //         << " U " << tp->Uout() << "(" << tp->TrackOut()->GetU0() << ")"
      //         << " V " << tp->Vout() << "(" << tp->TrackOut()->GetV0() << ")"
      //         << std::endl;
    }
  }
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchKurama()
{
  ClearKuramaTracks();

  auto nIn = m_SdcInTC.size();
  auto nOut = m_SdcOutTC.size();
  if(nIn==0 || nOut==0) return true;
  for(Int_t iIn=0; iIn<nIn; ++iIn){
    DCLocalTrack *trIn = GetTrackSdcIn(iIn);
    if(!trIn || !trIn->GoodForTracking()) continue;
    for(Int_t iOut=0; iOut<nOut; ++iOut){
      DCLocalTrack * trOut = GetTrackSdcOut(iOut);
      if(!trOut || !trOut->GoodForTracking()) continue;
      auto trKurama = new KuramaTrack(trIn, trOut);
      if(!trKurama) continue;
      Double_t u0In    = trIn->GetU0();
      Double_t u0Out   = trOut->GetU0();
      Double_t v0In    = trIn->GetV0();
      Double_t v0Out   = trOut->GetV0();
      Double_t bending = u0Out - u0In;
      Double_t p[3] = { 0.08493, 0.2227, 0.01572 };
      Double_t initial_momentum = p[0] + p[1]/(bending-p[2]);
      if(bending>0. && initial_momentum>0.){
        trKurama->SetInitialMomentum(initial_momentum);
      } else {
        trKurama->SetInitialMomentum(1.);
      }

      if(TMath::Abs(TMath::ATan(v0In) - TMath::ATan(v0Out)) < m_max_v0diff
         && trKurama->DoFit()
         && trKurama->ChiSquare()<MaxChiSqrKuramaTrack){
        // trKurama->Print("in "+FUNC_NAME);
        m_KuramaTC.push_back(trKurama);
      }
      else{
        // trKurama->Print("in "+FUNC_NAME);
        delete trKurama;
      }
    }
  }

  std::sort(m_KuramaTC.begin(), m_KuramaTC.end(), KuramaTrackComp());

#if 0
  PrintKurama("Before Deleting");
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchKurama(Double_t initial_momentum)
{
  ClearKuramaTracks();

  Int_t nIn  = GetNtracksSdcIn();
  Int_t nOut = GetNtracksSdcOut();

  if(nIn==0 || nOut==0) return true;

  for(Int_t iIn=0; iIn<nIn; ++iIn){
    DCLocalTrack *trIn = GetTrackSdcIn(iIn);
    if(!trIn->GoodForTracking()) continue;
    for(Int_t iOut=0; iOut<nOut; ++iOut){
      DCLocalTrack * trOut = GetTrackSdcOut(iOut);
      if(!trOut->GoodForTracking()) continue;
      KuramaTrack *trKurama = new KuramaTrack(trIn, trOut);
      if(!trKurama) continue;
      trKurama->SetInitialMomentum(initial_momentum);
      if(trKurama->DoFit() && trKurama->ChiSquare()<MaxChiSqrKuramaTrack){
        m_KuramaTC.push_back(trKurama);
      }
      else{
        trKurama->Print(" in "+FUNC_NAME);
        delete trKurama;
      }
    }// for(iOut)
  }// for(iIn)

  std::sort(m_KuramaTC.begin(), m_KuramaTC.end(), KuramaTrackComp());

#if 0
  PrintKurama("Before Deleting");
#endif

  return true;
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearDCHits()
{
#if UseBcIn
  ClearBcInHits();
#endif
  ClearBcOutHits();
  ClearSdcInHits();
  ClearSdcOutHits();
  ClearTOFHits();
}

//_____________________________________________________________________________
#if UseBcIn
void
DCAnalyzer::ClearBcInHits()
{
  del::ClearContainerAll(m_TempBcInHC);
  del::ClearContainerAll(m_BcInHC);
  del::ClearContainerAll(m_MWPCClCont);
}
#endif

//_____________________________________________________________________________
void
DCAnalyzer::ClearBcOutHits()
{
  del::ClearContainerAll(m_BcOutHC);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearSdcInHits()
{
  del::ClearContainerAll(m_SdcInHC);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearSdcOutHits()
{
  del::ClearContainerAll(m_SdcOutHC);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearVtxHits()
{
  del::ClearContainer(m_VtxPoint);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTOFHits()
{
  del::ClearContainer(m_TOFHC);
}

//_____________________________________________________________________________
#if UseBcIn
void
DCAnalyzer::ClearTracksBcIn()
{
  del::ClearContainer(m_BcInTC);
}
#endif

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksBcOut()
{
  del::ClearContainer(m_BcOutTC);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksSdcIn()
{
  del::ClearContainer(m_SdcInTC);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksSdcOut()
{
  del::ClearContainer(m_SdcOutTC);
}

//_____________________________________________________________________________
#if UseBcIn
void
DCAnalyzer::ClearK18TracksU2D()
{
  del::ClearContainer(m_K18U2DTC);
}
#endif

//_____________________________________________________________________________
void
DCAnalyzer::ClearK18TracksD2U()
{
  del::ClearContainer(m_K18D2UTC);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearKuramaTracks()
{
  del::ClearContainer(m_KuramaTC);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksBcOutSdcIn()
{
  del::ClearContainer(m_BcOutSdcInTC);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksSdcInSdcOut()
{
  del::ClearContainer(m_SdcInSdcOutTC);
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcMWPCHits(std::vector<DCHitContainer>& cont,
                           Bool_t applyRecursively)
{
  const std::size_t n = cont.size();
  for(std::size_t l=0; l<n; ++l){
    const std::size_t m = cont[l].size();
    for(std::size_t i=0; i<m; ++i){
      DCHit *hit = (cont[l])[i];
      if(!hit) continue;
      hit->ReCalcMWPC(applyRecursively);
    }
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcDCHits(std::vector<DCHitContainer>& cont,
                         Bool_t applyRecursively)
{
  const std::size_t n = cont.size();
  for(std::size_t l=0; l<n; ++l){
    const std::size_t m = cont[l].size();
    for(std::size_t i=0; i<m; ++i){
      DCHit *hit = (cont[l])[i];
      if(!hit) continue;
      hit->ReCalcDC(applyRecursively);
    }
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcDCHits(Bool_t applyRecursively)
{
#if UseBcIn
  ReCalcMWPCHits(m_TempBcInHC, applyRecursively);
  ReCalcMWPCHits(m_BcInHC, applyRecursively);
#endif

  ReCalcDCHits(m_BcOutHC, applyRecursively);
  ReCalcDCHits(m_SdcInHC, applyRecursively);
  ReCalcDCHits(m_SdcOutHC, applyRecursively);

  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcTrack(DCLocalTrackContainer& cont,
                        Bool_t applyRecursively)
{
  const std::size_t n = cont.size();
  for(std::size_t i=0; i<n; ++i){
    DCLocalTrack *track = cont[i];
    if(track) track->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcTrack(K18TrackD2UContainer& cont,
                        Bool_t applyRecursively)
{
  const std::size_t n = cont.size();
  for(std::size_t i=0; i<n; ++i){
    K18TrackD2U *track = cont[i];
    if(track) track->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcTrack(KuramaTrackContainer& cont,
                        Bool_t applyRecursively)
{
  const std::size_t n = cont.size();
  for(std::size_t i=0; i<n; ++i){
    KuramaTrack *track = cont[i];
    if(track) track->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
#if UseBcIn
Bool_t
DCAnalyzer::ReCalcTrackBcIn(Bool_t applyRecursively)
{
  return ReCalcTrack(m_BcInTC);
}
#endif

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcTrackBcOut(Bool_t applyRecursively)
{
  return ReCalcTrack(m_BcOutTC);
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcTrackSdcIn(Bool_t applyRecursively)
{
  return ReCalcTrack(m_SdcInTC);
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcTrackSdcOut(Bool_t applyRecursively)
{
  return ReCalcTrack(m_SdcOutTC);
}

//_____________________________________________________________________________
#if UseBcIn
Bool_t
DCAnalyzer::ReCalcK18TrackU2D(Bool_t applyRecursively)
{
  Int_t n = m_K18U2DTC.size();
  for(Int_t i=0; i<n; ++i){
    K18TrackU2D *track = m_K18U2DTC[i];
    if(track) track->ReCalc(applyRecursively);
  }
  return true;
}
#endif

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcK18TrackD2U(Bool_t applyRecursively)
{
  return ReCalcTrack(m_K18D2UTC, applyRecursively);
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcKuramaTrack(Bool_t applyRecursively)
{
  return ReCalcTrack(m_KuramaTC, applyRecursively);
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcAll()
{
  ReCalcDCHits();
#if UseBcIn
  ReCalcTrackBcIn();
  ReCalcK18TrackU2D();
#endif
  ReCalcTrackBcOut();
  ReCalcTrackSdcIn();
  ReCalcTrackSdcOut();

  //ReCalcK18TrackD2U();
  ReCalcKuramaTrack();

  return true;
}

//_____________________________________________________________________________
// Int_t
// clusterizeMWPCHit(const DCHitContainer& hits,
//     MWPCClusterContainer& clusters)
// {
//   if (!clusters.empty()){
//       std::for_each(clusters.begin(), clusters.end(), DeleteObject());
//       clusters.clear();
//   }

//   const Int_t nhits = hits.size();
//   //   hddaq::cout << __func__ << " " << nhits << std::endl;
//   if (nhits==0)
//     return 0;

//   Int_t n = 0;
//   for (Int_t i=0; i<nhits; ++i){
//     const DCHit* h = hits[i];
//     if (!h)
//       continue;
//     n += h->GetTdcSize();
//   }

//   DCHitContainer singleHits;
//   singleHits.reserve(n);
//   for (Int_t i=0; i<nhits; ++i){
//     const DCHit* h = hits[i];
//     if (!h)
//       continue;
//     Int_t nn = h->GetTdcSize();
//     for (Int_t ii=0; ii<nn; ++ii){
//       DCHit* htmp = new DCHit(h->GetLayer(), h->GetWire());
//       htmp->SetTdcVal(h->GetTdcVal());
//       htmp->SetTdcTrailing(h->GetTdcTrailing());
//       htmp->SetTrailingTime(h->GetTrailingTime());
//       htmp->SetDriftTime(h->GetDriftTime());
//       htmp->SetDriftLength(h->GetDriftLength());
//       htmp->SetTiltAngle(h->GetTiltAngle());
//       htmp->SetWirePosition(h->GetWirePosition());
//       htmp->setRangeCheckStatus(h->rangecheck(), 0);
//       singleHits.push_back(htmp);
//     }
//   }

//   std::vector<std::deque<Bool_t> > flag(n, std::deque<Bool_t>(n, false));
//   n = singleHits.size();
//   for (Int_t i=0;  i<n; ++i){
//     flag[i][i] = true;
//     const DCHit* h1 = singleHits[i];
//     //       h1->print("h1");
//     for (Int_t j=i+1; j<n; ++j){
//       const DCHit* h2 = singleHits[j];
//       //    h2->print("h2");
//       //    hddaq::cout << " (i,j) = (" << i << ", " << j << ")" << std::endl;
//       Bool_t val
//  = isConnectable(h1->GetWirePosition(),
//    h1->GetDriftTime(),
//    h1->GetTrailingTime(),
//    h2->GetWirePosition(),
//    h2->GetDriftTime(),
//    h2->GetTrailingTime(),
//    kMWPCClusteringWireExtension,
//    kMWPCClusteringTimeExtension);
//       //    hddaq::cout << "#D val = " << val << std::endl;
//       flag[i][j] = val;
//       flag[j][i] = val;
//     }
//   }

//   //   hddaq::cout << __func__ << "  before " << std::endl;
//   //   printConnectionFlag(flag);

//   const Int_t maxLoop = static_cast<Int_t>(std::log(x)/std::log(2.))+1;
//   for (Int_t loop=0; loop<maxLoop; ++loop){
//     std::vector<std::deque<Bool_t> > tmp(n, std::deque<Bool_t>(n, false));
//     for (Int_t i=0; i<n; ++i){
//       for (Int_t j=i; j<n; ++j){
//  for (Int_t k=0; k<n; ++k){
//    tmp[i][j] |= (flag[i][k] && flag[k][j]);
//    tmp[j][i] = tmp[i][j];
//  }
//       }
//     }
//     flag = tmp;
//     //       hddaq::cout << " n iteration = " << loop << std::endl;
//     //       printConnectionFlag(flag);
//   }

//   //   hddaq::cout << __func__ << "  after " << std::endl;
//   //   printConnectionFlag(flag);

//   std::set<Int_t> checked;
//   for (Int_t i=0; i<n; ++i){
//     if (checked.find(i)!=checked.end())
//       continue;
//     MWPCCluster* c = 0;
//     for (Int_t j=i; j<n; ++j){
//       if (flag[i][j]){
//  checked.insert(j);
//  if (!c) {
//    c = new MWPCCluster;
//    //     hddaq::cout << " new cluster " << std::endl;
//  }
//  //        hddaq::cout << " " << i << "---" << j << std::endl;
//  c->Add(singleHits[j]);
//       }
//     }

//     if (c){
//       c->Calculate();
//       clusters.push_back(c);
//     }
//   }

//   //   hddaq::cout << " end of " << __func__
//   //      << " : n = " << n << ", " << checked.size()
//   //      << std::endl;

//   //   hddaq::cout << __func__ << " n clusters = " << clusters.size() << std::endl;

//   return clusters.size();
// }

//_____________________________________________________________________________
void
DCAnalyzer::ChiSqrCutBcOut(Double_t chisqr)
{
  ChiSqrCut(m_BcOutTC, chisqr);
}

//_____________________________________________________________________________
void
DCAnalyzer::ChiSqrCutSdcIn(Double_t chisqr)
{
  ChiSqrCut(m_SdcInTC, chisqr);
}

//_____________________________________________________________________________
void
DCAnalyzer::ChiSqrCutSdcOut(Double_t chisqr)
{
  ChiSqrCut(m_SdcOutTC, chisqr);
}

//_____________________________________________________________________________
void
DCAnalyzer::ChiSqrCut(DCLocalTrackContainer& TrackCont,
                      Double_t chisqr)
{
  DCLocalTrackContainer DeleteCand;
  DCLocalTrackContainer ValidCand;
  for(auto& tempTrack : TrackCont){
    if(tempTrack->GetChiSquare() > chisqr){
      DeleteCand.push_back(tempTrack);
    }else{
      ValidCand.push_back(tempTrack);
    }
  }

  del::ClearContainer(DeleteCand);

  TrackCont.clear();
  TrackCont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), TrackCont.begin());
  ValidCand.clear();
}

//_____________________________________________________________________________
void
DCAnalyzer::ChiSqrCutKurama(Double_t chisqr)
{
  KuramaTrackContainer DeleteCand;
  KuramaTrackContainer ValidCand;
  for(auto& tempTrack : m_KuramaTC){
    if(tempTrack->GetChiSquare() > chisqr){
      DeleteCand.push_back(tempTrack);
    }else{
      ValidCand.push_back(tempTrack);
    }
  }

  del::ClearContainer(DeleteCand);

  m_KuramaTC.clear();
  m_KuramaTC.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), m_KuramaTC.begin());
  ValidCand.clear();
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCutBCOut(Double_t min_tot)
{
  for(Int_t i = 0; i<NumOfLayersBcOut; ++i){
    TotCut(m_BcOutHC[i + 1], min_tot, false);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCutSDC1(Double_t min_tot)
{
  for(Int_t i = 0; i<NumOfLayersSDC1; ++i){
    TotCut(m_SdcInHC[i + 1], min_tot, false);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCutSDC2(Double_t min_tot)
{
  for(Int_t i = 0; i<NumOfLayersSDC2; ++i){
    TotCut(m_SdcInHC[i + NumOfLayersSDC1 + 1], min_tot, false);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCutSDC3(Double_t min_tot)
{
  for(Int_t i = 0; i<NumOfLayersSDC3; ++i){
    TotCut(m_SdcOutHC[i + 1], min_tot, false);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCutSDC4(Double_t min_tot)
{
  for(Int_t i = 0; i<NumOfLayersSDC4; ++i){
    TotCut(m_SdcOutHC[i + NumOfLayersSDC3 + 1], min_tot, false);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCut(DCHitContainer& HitCont,
                   Double_t min_tot, Bool_t adopt_nan)
{
  DCHitContainer ValidCand;
  DCHitContainer DeleteCand;
  for(auto *ptr : HitCont){
    ptr->TotCut(min_tot, adopt_nan);
    if(0 == ptr->GetDriftTimeSize()){
      DeleteCand.push_back(ptr);
    }else{
      ValidCand.push_back(ptr);
    }
  }

  del::ClearContainer(DeleteCand);

  HitCont.clear();
  HitCont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), HitCont.begin());
  ValidCand.clear();
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCutBC34(Double_t min_dt, Double_t max_dt)
{
  for(Int_t i = 0; i<NumOfLayersBcOut; ++i){
    DriftTimeCut(m_BcOutHC[i + 1], min_dt, max_dt, true);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCutSDC1(Double_t min_dt, Double_t max_dt)
{
  for(Int_t i = 0; i<NumOfLayersSDC1; ++i){
    DriftTimeCut(m_SdcInHC[i + 1], min_dt, max_dt, true);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCutSDC2(Double_t min_dt, Double_t max_dt)
{
  for(Int_t i = 0; i<NumOfLayersSDC2; ++i){
    DriftTimeCut(m_SdcInHC[i + NumOfLayersSDC1 + 1], min_dt, max_dt, true);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCutSDC3(Double_t min_dt, Double_t max_dt)
{
  for(Int_t i = 0; i<NumOfLayersSDC3; ++i){
    DriftTimeCut(m_SdcOutHC[i + 1], min_dt, max_dt, true);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCutSDC4(Double_t min_dt, Double_t max_dt)
{
  for(Int_t i = 0; i<NumOfLayersSDC4; ++i){
    DriftTimeCut(m_SdcOutHC[i + NumOfLayersSDC3 + 1], min_dt, max_dt, true);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCut(DCHitContainer& HitCont,
                         Double_t min_dt, Double_t max_dt, Bool_t select_1st)
{
  DCHitContainer ValidCand;
  DCHitContainer DeleteCand;
  for(auto *ptr : HitCont){
    ptr->GateDriftTime(min_dt, max_dt, select_1st);
    if(0 == ptr->GetDriftTimeSize()){
      DeleteCand.push_back(ptr);
    }else{
      ValidCand.push_back(ptr);
    }
  }

  del::ClearContainer(DeleteCand);

  HitCont.clear();
  HitCont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), HitCont.begin());
  ValidCand.clear();
}

//Legacy
//_____________________________________________________________________________
Bool_t
DCAnalyzer::MakeBH2DCHit(Int_t t0seg)
{
  static const Double_t centerbh2[] = {
    -41.8, -19.3, -10.7, -3.6, 3.6, 10.7, 19.3, 41.8
  };

  Bool_t status = true;

  Double_t bh2pos = centerbh2[t0seg];
  DCHit *dchit = new DCHit(125, t0seg);
  dchit->SetTdcVal(0.);
  if(dchit->CalcFiberObservables()){
    dchit->SetWirePosition(bh2pos);
    m_BcOutHC[13].push_back(dchit);
  }else{
    delete dchit;
    status = false;
  }

  return status;
}

// -*- C++ -*-

#ifndef DC_ANALYZER_HH
#define DC_ANALYZER_HH

#include <map>
#include <vector>
#include <TString.h>
#include <TVector3.h>

#include "DetectorID.hh"

class DCHit;
class DCLocalTrack;
class K18TrackU2D;
class K18TrackD2U;
class S2sTrack;
class RawData;
class MWPCCluster;
class FiberCluster;
class HodoCluster;

class HodoHit;
class HodoAnalyzer;

using DCHC = std::vector<DCHit*>;
using DCLocalTC = std::vector<DCLocalTrack*>;
using K18TC = std::vector<K18TrackD2U*>;
using S2sTC = std::vector<S2sTrack*>;

using HodoHC = std::vector<HodoHit*>;
using HodoCC = std::vector<HodoCluster*>;

// no use
using MWPCClusterContainer = std::vector<MWPCCluster*>;
using K18TrackU2DContainer = std::vector<K18TrackU2D*>;

//_____________________________________________________________________________
class DCAnalyzer
{
public:
  static const TString& ClassName();
  DCAnalyzer();
  DCAnalyzer(const RawData& raw_data);
  ~DCAnalyzer();

private:
  DCAnalyzer(const DCAnalyzer&);
  DCAnalyzer& operator =(const DCAnalyzer&);

private:
  template <typename T> using map_t = std::map<TString, T>;
  const RawData*         m_raw_data;
  map_t<DCHC>            m_dc_hit_collection;
  Double_t               m_max_v0diff;
  std::vector<DCHC>      m_TempBcInHC;
  std::vector<DCHC>      m_BcInHC;
  std::vector<DCHC>      m_BcOutHC;
  std::vector<DCHC>      m_SdcInHC;
  std::vector<DCHC>      m_SdcOutHC;
  DCHC                   m_TOFHC;
  DCHC                   m_VtxPoint;
  DCLocalTC              m_BcInTC;
  DCLocalTC              m_BcOutTC;
  DCLocalTC              m_SdcInTC;
  DCLocalTC              m_SdcOutTC;
  K18TC                  m_K18D2UTC;
  S2sTC                  m_S2sTC;
  DCLocalTC              m_BcOutSdcInTC;
  DCLocalTC              m_SdcInSdcOutTC;
  std::vector<DCLocalTC> m_SdcInExTC;
  std::vector<DCLocalTC> m_SdcOutExTC;
  // no use
  std::vector<MWPCClusterContainer>  m_MWPCClCont;
  K18TrackU2DContainer               m_K18U2DTC;

public:
  void   DecodeHits(const TString& name);
  void   DecodeLocalHits(const TString& name);
  Bool_t DecodeRawHits();
  // Bool_t DecodeFiberHits(FiberCluster* FiberCl, Int_t layer);
  Bool_t DecodeFiberHits();
  Bool_t DecodeBcInHits();
  Bool_t DecodeBcOutHits();
  Bool_t DecodeSdcInHits();
  Bool_t DecodeSdcOutHits(Double_t ofs_dt=0.);
  Bool_t DecodeSdcInHitsGeant4(const std::vector<Int_t>& nhit,
			       const std::vector<std::vector<TVector3>>& pos,
			       const std::vector<std::vector<Double_t>>& de);
  Bool_t DecodeSdcOutHitsGeant4(const std::vector<Int_t>& nhit,
				const std::vector<std::vector<TVector3>>& pos,
				const std::vector<std::vector<Double_t>>& de);
  void DecodeHitsGeant4(const TString& name,
			const std::vector<Int_t>& planeId, const std::vector<Int_t>& layerId,
			const std::vector<TVector3>& gpos, const std::vector<Double_t>& de);
  Bool_t DecodeTOFHits(const HodoHC& HitCont);
  Bool_t DecodeTOFHits(const HodoCC& ClCont);
  // Bool_t DecodeSimuHits(SimuData *simuData);
  Int_t  ClusterizeMWPCHit(const DCHC& hits,
                           MWPCClusterContainer& clusters);
  Double_t GetMaxV0Diff() const { return m_max_v0diff; }
  const DCHC& GetTempBcInHC(Int_t l) const
    { return m_TempBcInHC.at(l); }
  const DCHC& GetBcInHC(Int_t l) const { return m_BcInHC.at(l); }
  const DCHC& GetBcOutHC(Int_t l) const { return m_BcOutHC.at(l); }
  const DCHC& GetSdcInHC(Int_t l) const { return m_SdcInHC.at(l); }
  const DCHC& GetSdcOutHC(Int_t l) const { return m_SdcOutHC.at(l); }
  const DCHC& GetTOFHC() const { return m_TOFHC; }

  Bool_t TrackSearchBcIn();
  Bool_t TrackSearchBcIn(const std::vector< std::vector<DCHC> >& hc);
  Bool_t TrackSearchBcOut(Int_t T0Seg=-1);
  Bool_t TrackSearchBcOut(const std::vector< std::vector<DCHC> >& hc, Int_t T0Seg);
  Bool_t TrackSearchSdcIn();
  Bool_t TrackSearchSdcInFiber();
  Bool_t TrackSearchSdcOut();
  Bool_t TrackSearchSdcOut(const HodoHC& HitCont);
  Bool_t TrackSearchSdcOut(const HodoCC& ClCont);

  Int_t GetNtracksBcIn() const { return m_BcInTC.size(); }
  Int_t GetNtracksBcOut() const { return m_BcOutTC.size(); }
  Int_t GetNtracksSdcIn() const { return m_SdcInTC.size(); }
  Int_t GetNtracksSdcOut() const { return m_SdcOutTC.size(); }
  // Exclusive Tracks
  Int_t GetNtracksSdcInEx(Int_t l) const { return m_SdcInExTC.at(l).size(); }
  Int_t GetNtracksSdcOutEx(Int_t l) const { return m_SdcOutExTC.at(l).size(); }

  const DCLocalTrack* GetTrackBcIn(Int_t l) const { return m_BcInTC.at(l); }
  const DCLocalTrack* GetTrackBcOut(Int_t l) const { return m_BcOutTC.at(l); }
  const DCLocalTrack* GetTrackSdcIn(Int_t l) const { return m_SdcInTC.at(l); }
  const DCLocalTrack* GetTrackSdcOut(Int_t l) const { return m_SdcOutTC.at(l); }
  // Exclusive Tracks
  const DCLocalTrack* GetTrackSdcInEx(Int_t l, Int_t i) const
    { return m_SdcInExTC.at(l).at(i); }
  const DCLocalTrack* GetTrackSdcOutEx(Int_t l, Int_t i) const
    { return m_SdcOutExTC.at(l).at(i); }

  Bool_t TrackSearchK18U2D();
  Bool_t TrackSearchK18D2U(const std::vector<Double_t>& XinCont);
  Bool_t TrackSearchS2s(Double_t initial_momentum);
  Bool_t TrackSearchS2s();

  void ChiSqrCutBcOut(Double_t chisqr);
  void ChiSqrCutSdcIn(Double_t chisqr);
  void ChiSqrCutSdcOut(Double_t chisqr);

  void TotCutBCOut(Double_t min_tot);
  void TotCutSDC1(Double_t min_tot);
  void TotCutSDC2(Double_t min_tot);
  void TotCutSDC3(Double_t min_tot);
  void TotCutSDC4(Double_t min_tot);
  void TotCutSDC5(Double_t min_tot);

  void DriftTimeCutBC34(Double_t min_dt, Double_t max_dt);
  void DriftTimeCutSDC1(Double_t min_dt, Double_t max_dt);
  void DriftTimeCutSDC2(Double_t min_dt, Double_t max_dt);
  void DriftTimeCutSDC3(Double_t min_dt, Double_t max_dt);
  void DriftTimeCutSDC4(Double_t min_dt, Double_t max_dt);
  void DriftTimeCutSDC5(Double_t min_dt, Double_t max_dt);

  Int_t GetNTracksK18U2D() const { return m_K18U2DTC.size(); }
  Int_t GetNTracksK18D2U() const { return m_K18D2UTC.size(); }
  Int_t GetNTracksS2s() const { return m_S2sTC.size(); }

  const K18TrackU2D* GetK18TrackU2D(Int_t l) const { return m_K18U2DTC.at(l); }
  const K18TrackD2U* GetK18TrackD2U(Int_t l) const { return m_K18D2UTC.at(l); }
  const S2sTrack* GetS2sTrack(Int_t l) const { return m_S2sTC.at(l); }
  const K18TC& GetK18TracksD2U() const { return m_K18D2UTC; }
  const S2sTC& GetS2sTracks() const { return m_S2sTC; }

  Int_t GetNClustersMWPC(Int_t l) const { return m_MWPCClCont.at(l).size(); };

  const MWPCClusterContainer& GetClusterMWPC(Int_t l) const
    { return m_MWPCClCont.at(l); }
  void PrintS2s(const TString& arg="") const;

  Bool_t ReCalcMWPCHits(std::vector<DCHC>& cont,
                        Bool_t applyRecursively=false);
  Bool_t ReCalcDCHits(std::vector<DCHC>& cont,
                      Bool_t applyRecursively=false);
  Bool_t ReCalcDCHits(Bool_t applyRecursively=false);
  void HoughYCut(Double_t min_y, Double_t max_y);
  Bool_t ReCalcTrack(DCLocalTC& cont, Bool_t applyRecursively=false);
  Bool_t ReCalcTrack(K18TC& cont, Bool_t applyRecursively=false);
  Bool_t ReCalcTrack(S2sTC& cont, Bool_t applyRecursively=false);

  Bool_t ReCalcTrackBcIn(Bool_t applyRecursively=false);
  Bool_t ReCalcTrackBcOut(Bool_t applyRecursively=false);
  Bool_t ReCalcTrackSdcIn(Bool_t applyRecursively=false);
  Bool_t ReCalcTrackSdcOut(Bool_t applyRecursively=false);

  Bool_t ReCalcK18TrackD2U(Bool_t applyRecursively=false);
  // Bool_t ReCalcK18TrackU2D(Bool_t applyRecursively=false);
  Bool_t ReCalcS2sTrack(Bool_t applyRecursively=false);
  Bool_t ReCalcAll();
  void   SetMaxV0Diff(Double_t deg){ m_max_v0diff = deg; }
  Bool_t TrackSearchBcOutSdcIn();
  Bool_t TrackSearchSdcInSdcOut();
  const DCLocalTC& GetTrackContainerBcOutSdcIn() const
    { return m_BcOutSdcInTC; }
  Int_t GetNtracksBcOutSdcIn() const { return m_BcOutSdcInTC.size(); }
  Int_t GetNtracksSdcInSdcOut() const { return m_SdcInSdcOutTC.size(); }
  const DCLocalTrack* GetTrackBcOutSdcIn(Int_t i) const
    { return m_BcOutSdcInTC.at(i); }
  const DCLocalTrack* GetTrackSdcInSdcOut(Int_t i) const
    { return m_SdcInSdcOutTC.at(i); }

  Bool_t MakeBH2DCHit(Int_t t0seg);

protected:
  void ClearDCHits();
  void ClearBcInHits();
  void ClearBcOutHits();
  void ClearSdcInHits();
  void ClearSdcOutHits();

  void ClearTOFHits();
  void ClearVtxHits();

  void ClearTracksBcIn();
  void ClearTracksBcOut();
  void ClearTracksSdcIn();
  void ClearTracksSdcOut();
  void ClearTracksBcOutSdcIn();
  void ClearTracksSdcInSdcOut();
  void ClearK18TracksU2D();
  void ClearK18TracksD2U();
  void ClearS2sTracks();
  void ChiSqrCut(DCLocalTC& cont, Double_t chisqr);
  void EraseEmptyHits(const TString& name);
  void EraseEmptyHits(std::vector<DCHC>& HitCont);
  void TotCut(const TString& name, Double_t min_tot, Bool_t keep_nan);
  void DriftTimeCut(const TString& name, Double_t min_dt, Double_t max_dt, Bool_t select_1st);
  static Int_t MakeUpMWPCClusters(const DCHC& HitCont,
                                  MWPCClusterContainer& ClusterCont,
                                  Double_t maxTimeDif);

public:
  void ResetTracksBcIn()        { ClearTracksBcIn();        }
  void ResetTracksBcOut()       { ClearTracksBcOut();       }
  void ResetTracksSdcIn()       { ClearTracksSdcIn();       }
  void ResetTracksSdcOut()      { ClearTracksSdcOut();      }
  void ResetTracksBcOutSdcIn()  { ClearTracksBcOutSdcIn();  }
  void ResetTracksSdcInSdcOut()  { ClearTracksSdcInSdcOut();  }
  void ApplyBh1SegmentCut(const std::vector<Double_t>& validBh1Cluster);
  void ApplyBh2SegmentCut(const Double_t Time0_Cluster);

protected:
};

//_____________________________________________________________________________
inline const TString&
DCAnalyzer::ClassName()
{
  static TString s_name("DCAnalyzer");
  return s_name;
}

#endif

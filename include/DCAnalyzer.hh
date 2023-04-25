// -*- C++ -*-

#ifndef DC_ANALYZER_HH
#define DC_ANALYZER_HH

#include <vector>
#include <TString.h>

#include "DetectorID.hh"

class DCHit;
class DCLocalTrack;
class K18TrackU2D;
class K18TrackD2U;
class KuramaTrack;
class RawData;
class MWPCCluster;
class FiberCluster;
class HodoCluster;

class Hodo1Hit;
class Hodo2Hit;
class HodoAnalyzer;

typedef std::vector<DCHit*>        DCHitContainer;
typedef std::vector<MWPCCluster*>  MWPCClusterContainer;
typedef std::vector<DCLocalTrack*> DCLocalTrackContainer;
typedef std::vector<K18TrackU2D*>  K18TrackU2DContainer;
typedef std::vector<K18TrackD2U*>  K18TrackD2UContainer;
typedef std::vector<KuramaTrack*>  KuramaTrackContainer;

typedef std::vector<Hodo1Hit*> Hodo1HitContainer;
typedef std::vector<Hodo2Hit*> Hodo2HitContainer;
typedef std::vector<HodoCluster*> HodoClusterContainer;

//_____________________________________________________________________________
class DCAnalyzer
{
public:
  static const TString& ClassName();
  DCAnalyzer();
  ~DCAnalyzer();

private:
  DCAnalyzer(const DCAnalyzer&);
  DCAnalyzer& operator =(const DCAnalyzer&);

private:
  enum e_type
  { kBcIn, kBcOut, kSdcIn, kSdcOut, kTOF, n_type };
  Double_t                           m_max_v0diff;
  std::vector<Bool_t>                m_is_decoded;
  std::vector<Int_t>                 m_much_combi;
  std::vector<MWPCClusterContainer>  m_MWPCClCont;
  std::vector<DCHitContainer>        m_TempBcInHC;
  std::vector<DCHitContainer>        m_BcInHC;
  std::vector<DCHitContainer>        m_BcOutHC;
  std::vector<DCHitContainer>        m_SdcInHC;
  std::vector<DCHitContainer>        m_SdcOutHC;
  DCHitContainer                     m_TOFHC;
  DCHitContainer                     m_VtxPoint;
  DCLocalTrackContainer              m_BcInTC;
  DCLocalTrackContainer              m_BcOutTC;
  DCLocalTrackContainer              m_SdcInTC;
  DCLocalTrackContainer              m_SdcOutTC;
  K18TrackU2DContainer               m_K18U2DTC;
  K18TrackD2UContainer               m_K18D2UTC;
  KuramaTrackContainer               m_KuramaTC;
  DCLocalTrackContainer              m_BcOutSdcInTC;
  DCLocalTrackContainer              m_SdcInSdcOutTC;
  std::vector<DCLocalTrackContainer> m_SdcInExTC;
  std::vector<DCLocalTrackContainer> m_SdcOutExTC;

public:
  Int_t  MuchCombinationSdcIn() const { return m_much_combi[kSdcIn]; }
  Bool_t DecodeRawHits(RawData* rawData);
  // Bool_t DecodeFiberHits(FiberCluster* FiberCl, Int_t layer);
  Bool_t DecodeFiberHits(RawData* rawData);
  Bool_t DecodeBcInHits(RawData* rawData);
  Bool_t DecodeBcOutHits(RawData* rawData);
  Bool_t DecodeSdcInHits(RawData* rawData);
  Bool_t DecodeSdcOutHits(RawData* rawData, Double_t ofs_dt=0.);
  Bool_t DecodeTOFHits(const Hodo2HitContainer& HitCont);
  Bool_t DecodeTOFHits(const HodoClusterContainer& ClCont);
  // Bool_t DecodeSimuHits(SimuData *simuData);
  Int_t  ClusterizeMWPCHit(const DCHitContainer& hits,
                           MWPCClusterContainer& clusters);
  Double_t GetMaxV0Diff() const { return m_max_v0diff; }
  const DCHitContainer& GetTempBcInHC(Int_t l) const
    { return m_TempBcInHC.at(l); }
  const DCHitContainer& GetBcInHC(Int_t l) const { return m_BcInHC.at(l); }
  const DCHitContainer& GetBcOutHC(Int_t l) const { return m_BcOutHC.at(l); }
  const DCHitContainer& GetSdcInHC(Int_t l) const { return m_SdcInHC.at(l); }
  const DCHitContainer& GetSdcOutHC(Int_t l) const { return m_SdcOutHC.at(l); }
  const DCHitContainer& GetTOFHC() const { return m_TOFHC; }

  Bool_t TrackSearchBcIn();
  Bool_t TrackSearchBcIn(const std::vector< std::vector<DCHitContainer> >& hc);
  Bool_t TrackSearchBcOut(Int_t T0Seg=-1);
  Bool_t TrackSearchBcOut(const std::vector< std::vector<DCHitContainer> >& hc, Int_t T0Seg);
  Bool_t TrackSearchSdcIn();
  Bool_t TrackSearchSdcInFiber();
  Bool_t TrackSearchSdcOut();
  Bool_t TrackSearchSdcOut(const Hodo2HitContainer& HitCont);
  Bool_t TrackSearchSdcOut(const HodoClusterContainer& ClCont);

  Int_t GetNtracksBcIn() const { return m_BcInTC.size(); }
  Int_t GetNtracksBcOut() const { return m_BcOutTC.size(); }
  Int_t GetNtracksSdcIn() const { return m_SdcInTC.size(); }
  Int_t GetNtracksSdcOut() const { return m_SdcOutTC.size(); }
  // Exclusive Tracks
  Int_t GetNtracksSdcInEx(Int_t l) const { return m_SdcInExTC.at(l).size(); }
  Int_t GetNtracksSdcOutEx(Int_t l) const { return m_SdcOutExTC.at(l).size(); }

  DCLocalTrack* GetTrackBcIn(Int_t l) const { return m_BcInTC.at(l); }
  DCLocalTrack* GetTrackBcOut(Int_t l) const { return m_BcOutTC.at(l); }
  DCLocalTrack* GetTrackSdcIn(Int_t l) const { return m_SdcInTC.at(l); }
  DCLocalTrack* GetTrackSdcOut(Int_t l) const { return m_SdcOutTC.at(l); }
  // Exclusive Tracks
  DCLocalTrack* GetTrackSdcInEx(Int_t l, Int_t i) const
    { return m_SdcInExTC.at(l).at(i); }
  DCLocalTrack* GetTrackSdcOutEx(Int_t l, Int_t i) const
    { return m_SdcOutExTC.at(l).at(i); }

  Bool_t TrackSearchK18U2D();
  Bool_t TrackSearchK18D2U(const std::vector<Double_t>& XinCont);
  Bool_t TrackSearchKurama(Double_t initial_momentum);
  Bool_t TrackSearchKurama();

  void ChiSqrCutBcOut(Double_t chisqr);
  void ChiSqrCutSdcIn(Double_t chisqr);
  void ChiSqrCutSdcOut(Double_t chisqr);

  void TotCutBCOut(Double_t min_tot);
  void TotCutSDC1(Double_t min_tot);
  void TotCutSDC2(Double_t min_tot);
  void TotCutSDC3(Double_t min_tot);
  void TotCutSDC4(Double_t min_tot);

  void DriftTimeCutBC34(Double_t min_dt, Double_t max_dt);
  void DriftTimeCutSDC1(Double_t min_dt, Double_t max_dt);
  void DriftTimeCutSDC2(Double_t min_dt, Double_t max_dt);
  void DriftTimeCutSDC3(Double_t min_dt, Double_t max_dt);
  void DriftTimeCutSDC4(Double_t min_dt, Double_t max_dt);

  Int_t GetNTracksK18U2D() const { return m_K18U2DTC.size(); }
  Int_t GetNTracksK18D2U() const { return m_K18D2UTC.size(); }
  Int_t GetNTracksKurama() const { return m_KuramaTC.size(); }

  K18TrackU2D* GetK18TrackU2D(Int_t l) const { return m_K18U2DTC.at(l); }
  K18TrackD2U* GetK18TrackD2U(Int_t l) const { return m_K18D2UTC.at(l); }
  KuramaTrack* GetKuramaTrack(Int_t l) const { return m_KuramaTC.at(l); }
  const K18TrackD2UContainer& GetK18TracksD2U() const { return m_K18D2UTC; }
  const KuramaTrackContainer& GetKuramaTracks() const { return m_KuramaTC; }

  Int_t GetNClustersMWPC(Int_t l) const { return m_MWPCClCont.at(l).size(); };

  const MWPCClusterContainer& GetClusterMWPC(Int_t l) const
    { return m_MWPCClCont.at(l); }
  void PrintKurama(const TString& arg="") const;

  Bool_t ReCalcMWPCHits(std::vector<DCHitContainer>& cont,
                        Bool_t applyRecursively=false);
  Bool_t ReCalcDCHits(std::vector<DCHitContainer>& cont,
                      Bool_t applyRecursively=false);
  Bool_t ReCalcDCHits(Bool_t applyRecursively=false);
  void HoughYCut(Double_t min_y, Double_t max_y);
  Bool_t ReCalcTrack(DCLocalTrackContainer& cont, Bool_t applyRecursively=false);
  Bool_t ReCalcTrack(K18TrackD2UContainer& cont, Bool_t applyRecursively=false);
  Bool_t ReCalcTrack(KuramaTrackContainer& cont, Bool_t applyRecursively=false);

  Bool_t ReCalcTrackBcIn(Bool_t applyRecursively=false);
  Bool_t ReCalcTrackBcOut(Bool_t applyRecursively=false);
  Bool_t ReCalcTrackSdcIn(Bool_t applyRecursively=false);
  Bool_t ReCalcTrackSdcOut(Bool_t applyRecursively=false);

  Bool_t ReCalcK18TrackD2U(Bool_t applyRecursively=false);
  // Bool_t ReCalcK18TrackU2D(Bool_t applyRecursively=false);
  Bool_t ReCalcKuramaTrack(Bool_t applyRecursively=false);
  Bool_t ReCalcAll();
  void   SetMaxV0Diff(Double_t deg){ m_max_v0diff = deg; }
  Bool_t TrackSearchBcOutSdcIn();
  Bool_t TrackSearchSdcInSdcOut();
  const DCLocalTrackContainer& GetTrackContainerBcOutSdcIn() const
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
  void ClearKuramaTracks();
  void ChiSqrCut(DCLocalTrackContainer& cont, Double_t chisqr);
  void TotCut(DCHitContainer& cont, Double_t min_tot, Bool_t adopt_nan);
  void DriftTimeCut(DCHitContainer& cont, Double_t min_dt, Double_t max_dt, Bool_t select_1st);
  static Int_t MakeUpMWPCClusters(const DCHitContainer& HitCont,
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

};

//_____________________________________________________________________________
inline const TString&
DCAnalyzer::ClassName()
{
  static TString s_name("DCAnalyzer");
  return s_name;
}

#endif

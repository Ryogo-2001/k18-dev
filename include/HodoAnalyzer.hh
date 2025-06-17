// -*- C++ -*-

#ifndef HODO_ANALYZER_HH
#define HODO_ANALYZER_HH

#include <vector>

#include <TString.h>

#include "DetectorID.hh"
#include "RawData.hh"

class RawData;
class Hodo1Hit;
class Hodo2Hit;
class BH2Hit;
class FiberHit;
class FLHit;
class HodoCluster;
class BH2Cluster;
class FiberCluster;

typedef std::vector<Hodo1Hit*> Hodo1HitContainer;
typedef std::vector<Hodo2Hit*> Hodo2HitContainer;
typedef std::vector<BH2Hit*>   BH2HitContainer;
typedef std::vector<FiberHit*> FiberHitContainer;
typedef std::vector<FLHit*>    FLHitContainer;
typedef std::vector< std::vector<FiberHit*> >
MultiPlaneFiberHitContainer;

typedef std::vector<HodoCluster*>  HodoClusterContainer;
typedef std::vector<BH2Cluster*>   BH2ClusterContainer;
typedef std::vector<FiberCluster*> FiberClusterContainer;
typedef std::vector< std::vector<FiberCluster*> >
MultiPlaneFiberClusterContainer;

//_____________________________________________________________________________
class HodoAnalyzer
{
public:
  static TString       ClassName();
  static HodoAnalyzer& GetInstance();
  HodoAnalyzer();
  ~HodoAnalyzer();

private:
  HodoAnalyzer(const HodoAnalyzer&);
  HodoAnalyzer& operator =(const HodoAnalyzer&);

private:
  Hodo2HitContainer           m_BH1Cont;
  BH2HitContainer             m_BH2Cont;
  Hodo1HitContainer           m_BACCont;
  Hodo2HitContainer           m_HTOFCont;
  Hodo1HitContainer           m_BVHCont;
  Hodo2HitContainer           m_TOFCont;
  Hodo1HitContainer           m_LACCont;
  Hodo2HitContainer           m_WCCont;
  Hodo1HitContainer           m_WCSUMCont;
  MultiPlaneFiberHitContainer m_BFTCont;
  FiberHitContainer           m_SCHCont;
  Hodo1Hit*                   m_TPCClock;
  HodoClusterContainer        m_BH1ClCont;
  BH2ClusterContainer         m_BH2ClCont;
  HodoClusterContainer        m_BACClCont;
  HodoClusterContainer        m_HTOFClCont;
  HodoClusterContainer        m_BVHClCont;
  HodoClusterContainer        m_TOFClCont;
  HodoClusterContainer        m_LACClCont;
  HodoClusterContainer        m_WCClCont;
  HodoClusterContainer        m_WCSUMClCont;
  FiberClusterContainer       m_BFTClCont;
  FiberClusterContainer       m_SCHClCont;

public:
  Bool_t DecodeRawHits(RawData* rawData);
  Bool_t DecodeBH1Hits(RawData* rawData);
  Bool_t DecodeBH2Hits(RawData* rawData);
  Bool_t DecodeBACHits(RawData* rawData);
  Bool_t DecodeHTOFHits(RawData* rawData);
  Bool_t DecodeBVHHits(RawData* rawData);
  Bool_t DecodeTOFHits(RawData* rawData);
  Bool_t DecodeLACHits(RawData* rawData);
  Bool_t DecodeWCHits(RawData* rawData);
  Bool_t DecodeWCSUMHits(RawData* rawData);
  Bool_t DecodeBFTHits(RawData* rawData);
  Bool_t DecodeSCHHits(RawData* rawData);
  Bool_t DecodeTPCClock(RawData* rawData);
  Int_t  GetNHitsBH1() const { return m_BH1Cont.size(); };
  Int_t  GetNHitsBH2() const { return m_BH2Cont.size(); };
  Int_t  GetNHitsBAC() const { return m_BACCont.size(); };
  Int_t  GetNHitsHTOF() const { return m_HTOFCont.size(); };
  Int_t  GetNHitsBVH() const { return m_BVHCont.size(); };
  Int_t  GetNHitsTOF() const { return m_TOFCont.size(); };
  Int_t  GetNHitsLAC() const { return m_LACCont.size(); };
  Int_t  GetNHitsWC() const { return m_WCCont.size(); };
  Int_t  GetNHitsWCSUM() const { return m_WCSUMCont.size(); };
  Int_t  GetNHitsBFT(Int_t plane) const
  { return m_BFTCont.at(plane).size(); };
  Int_t  GetNHitsSCH()  const { return m_SCHCont.size(); };

  Hodo2Hit* GetHitBH1(UInt_t i) const { return m_BH1Cont.at(i); }
  BH2Hit*   GetHitBH2(UInt_t i) const { return m_BH2Cont.at(i); }
  Hodo1Hit* GetHitBAC(UInt_t i) const { return m_BACCont.at(i); }
  Hodo2Hit* GetHitHTOF(UInt_t i) const { return m_HTOFCont.at(i); }
  Hodo1Hit* GetHitBVH(UInt_t i) const { return m_BVHCont.at(i); }
  Hodo2Hit* GetHitTOF(UInt_t i) const { return m_TOFCont.at(i); }
  Hodo1Hit* GetHitLAC(UInt_t i) const { return m_LACCont.at(i); }
  Hodo2Hit* GetHitWC(UInt_t i) const { return m_WCCont.at(i); }
  Hodo1Hit* GetHitWCSUM(UInt_t i) const { return m_WCSUMCont.at(i); }
  FiberHit* GetHitBFT(Int_t plane, UInt_t seg) const { return m_BFTCont.at(plane).at(seg); }
  FiberHit* GetHitSCH(UInt_t seg) const { return m_SCHCont.at(seg); }
  Hodo1Hit* GetHitTPCClock() const { return m_TPCClock; }

  const Hodo2HitContainer& GetHitsBH1() const { return m_BH1Cont; }
  const BH2HitContainer& GetHitsBH2() const { return m_BH2Cont; }
  const Hodo2HitContainer& GetHitsTOF() const { return m_TOFCont; }
  const FiberHitContainer& GetHitsSCH() const { return m_SCHCont; }

  Int_t GetNClustersBH1() const { return m_BH1ClCont.size(); }
  Int_t GetNClustersBH2() const { return m_BH2ClCont.size(); }
  Int_t GetNClustersBAC() const { return m_BACClCont.size(); }
  Int_t GetNClustersHTOF() const { return m_HTOFClCont.size(); }
  Int_t GetNClustersBVH() const { return m_BVHClCont.size(); }
  Int_t GetNClustersTOF() const { return m_TOFClCont.size(); }
  Int_t GetNClustersLAC() const { return m_LACClCont.size(); }
  Int_t GetNClustersWC() const { return m_WCClCont.size(); }
  Int_t GetNClustersWCSUM() const { return m_WCSUMClCont.size(); }
  Int_t GetNClustersBFT() const { return m_BFTClCont.size(); };
  Int_t GetNClustersSCH() const { return m_SCHClCont.size(); };
  HodoCluster*  GetClusterBH1(UInt_t i) const { return m_BH1ClCont.at(i); }
  BH2Cluster*   GetClusterBH2(UInt_t i) const { return m_BH2ClCont.at(i); }
  HodoCluster*  GetClusterBAC(UInt_t i) const { return m_BACClCont.at(i); }
  HodoCluster*  GetClusterHTOF(UInt_t i) const { return m_HTOFClCont.at(i); }
  HodoCluster*  GetClusterBVH(UInt_t i) const { return m_BVHClCont.at(i); }
  HodoCluster*  GetClusterTOF(UInt_t i) const { return m_TOFClCont.at(i); }
  HodoCluster*  GetClusterLAC(UInt_t i) const { return m_LACClCont.at(i); }
  HodoCluster*  GetClusterWC(UInt_t i) const { return m_WCClCont.at(i); }
  HodoCluster*  GetClusterWCSUM(UInt_t i) const { return m_WCSUMClCont.at(i); }
  FiberCluster* GetClusterBFT(UInt_t i) const { return m_BFTClCont.at(i); }
  FiberCluster* GetClusterSCH(UInt_t i) const { return m_SCHClCont.at(i); }
  const FiberClusterContainer& GetClustersBFT() const { return m_BFTClCont; }
  const FiberClusterContainer& GetClustersSCH() const { return m_SCHClCont; }

  Bool_t ReCalcBH1Hits(Bool_t applyRecursively=false);
  Bool_t ReCalcBH2Hits(Bool_t applyRecursively=false);
  Bool_t ReCalcBACHits(Bool_t applyRecursively=false);
  Bool_t ReCalcHTOFHits(Bool_t applyRecursively=false);
  Bool_t ReCalcBVHHits(Bool_t applyRecursively=false);
  Bool_t ReCalcTOFHits(Bool_t applyRecursively=false);
  Bool_t ReCalcLACHits(Bool_t applyRecursively=false);
  Bool_t ReCalcWCHits(Bool_t applyRecursively=false);
  Bool_t ReCalcWCSUMHits(Bool_t applyRecursively=false);
  Bool_t ReCalcBH1Clusters(Bool_t applyRecursively=false);
  Bool_t ReCalcBH2Clusters(Bool_t applyRecursively=false);
  Bool_t ReCalcBACClusters(Bool_t applyRecursively=false);
  Bool_t ReCalcHTOFClusters(Bool_t applyRecursively=false);
  Bool_t ReCalcBVHClusters(Bool_t applyRecursively=false);
  Bool_t ReCalcTOFClusters(Bool_t applyRecursively=false);
  Bool_t ReCalcLACClusters(Bool_t applyRecursively=false);
  Bool_t ReCalcWCClusters(Bool_t applyRecursively=false);
  Bool_t ReCalcWCSUMClusters(Bool_t applyRecursively=false);
  Bool_t ReCalcAll();

  void TimeCutBH1(Double_t tmin, Double_t tmax);
  void TimeCutBH2(Double_t tmin, Double_t tmax);
  void TimeCutHTOF(Double_t tmin, Double_t tmax);
  void TimeCutBVH(Double_t tmin, Double_t tmax);
  void TimeCutTOF(Double_t tmin, Double_t tmax);
  void TimeCutBFT(Double_t tmin, Double_t tmax);
  void TimeCutSCH(Double_t tmin, Double_t tmax);

  void WidthCutBFT(Double_t min_width, Double_t max_width);
  void WidthCutSCH(Double_t min_width, Double_t max_width);
  BH2Cluster*  GetTime0BH2Cluster();
  HodoCluster* GetBtof0BH1Cluster(Double_t time0);

private:
  void ClearBH1Hits();
  void ClearBH2Hits();
  void ClearBACHits();
  void ClearHTOFHits();
  void ClearBVHHits();
  void ClearTOFHits();
  void ClearLACHits();
  void ClearWCHits();
  void ClearWCSUMHits();
  void ClearBFTHits();
  void ClearSCHHits();
  void ClearTPCClock();

  template<typename TypeCluster>
  void TimeCut(std::vector<TypeCluster>& cont, Double_t tmin, Double_t tmax);
  template<typename TypeCluster>
  void WidthCut(std::vector<TypeCluster>& cont,
		 Double_t min_width, Double_t max_width, Bool_t adopt_nan);
  template<typename TypeCluster>
  void WidthCutR(std::vector<TypeCluster>& cont,
		  Double_t min_width, Double_t max_width, Bool_t adopt_nan);
  template<typename TypeCluster>
  void AdcCut(std::vector<TypeCluster>& cont, Double_t amin, Double_t amax);
  static Int_t MakeUpClusters(const Hodo1HitContainer& HitCont,
                              HodoClusterContainer& ClusterCont,
                              Double_t maxTimeDif);
  static Int_t MakeUpClusters(const Hodo2HitContainer& HitCont,
                              HodoClusterContainer& ClusterCont,
                              Double_t maxTimeDif);
  static Int_t MakeUpClusters(const BH2HitContainer& HitCont,
                              BH2ClusterContainer& ClusterCont,
                              Double_t maxTimeDif);
  static Int_t MakeUpClusters(const FiberHitContainer& cont,
                              FiberClusterContainer& ClusterCont,
                              Double_t maxTimeDif,
                              Int_t DifPairId);
  static Int_t MakeUpClusters(const FLHitContainer& cont,
                              FiberClusterContainer& ClusterCont,
                              Double_t maxTimeDif,
                              Int_t DifPairId);
  static Int_t MakeUpCoincidence(const FiberHitContainer& cont,
                                 FLHitContainer& CoinCont,
                                 Double_t maxTimeDif);
};

//_____________________________________________________________________________
inline TString
HodoAnalyzer::ClassName()
{
  static TString s_name("HodoAnalyzer");
  return s_name;
}

//_____________________________________________________________________________
inline HodoAnalyzer&
HodoAnalyzer::GetInstance()
{
  static HodoAnalyzer g_instance;
  return g_instance;
}

#endif

// -*- C++ -*-

#ifndef HODO_ANALYZER_HH
#define HODO_ANALYZER_HH

#include <vector>

#include <TString.h>

#include "DeleteUtility.hh"
#include "DetectorID.hh"
#include "FiberCluster.hh"
#include "FuncName.hh"
#include "HodoCluster.hh"
#include "HodoHit.hh"
#include "RawData.hh"
#include "UserParamMan.hh"

class RawData;
class BH2Hit;
class FiberHit;
// class BH2Cluster;

using HodoClusterContainer = std::vector<HodoCluster*>;
using HodoCC = std::vector<HodoCluster*>;

//_____________________________________________________________________________
class HodoAnalyzer
{
public:
  static TString ClassName();
  explicit HodoAnalyzer(const RawData& raw_data);
  ~HodoAnalyzer();

private:
  HodoAnalyzer(const HodoAnalyzer&);
  HodoAnalyzer& operator =(const HodoAnalyzer&);

private:
  template <typename T> using map_t = std::map<TString, T>;
  const RawData* m_raw_data;
  map_t<HodoHC>  m_hodo_hit_collection;
  map_t<HodoCC>  m_hodo_cluster_collection;

public:
  template <typename T=HodoHit>
  Bool_t DecodeHits(const TString& name, Double_t max_time_diff=10.);

  const HodoHC& GetHitContainer(const TString& name) const;
  const HodoCC& GetClusterContainer(const TString& name) const;

  Int_t  GetNHits(const TString& name) const
    { return GetHitContainer(name).size(); };
  Int_t  GetNClusters(const TString& name) const
    { return GetClusterContainer(name).size(); };

  template <typename T=HodoHit>
  const T* GetHit(const TString& name, Int_t i=0) const
    { return dynamic_cast<T*>(GetHitContainer(name).at(i)); }

  template <typename T=HodoCluster>
  const T* GetCluster(const TString& name, Int_t i=0) const
    { return dynamic_cast<T*>(GetClusterContainer(name).at(i)); }

  Bool_t ReCalcHit(const TString& name, Bool_t applyRecursively=false);
  Bool_t ReCalcCluster(const TString& name, Bool_t applyRecursively=false);
  Bool_t ReCalcAll();

  void TimeCut(const TString& name, Double_t min, Double_t max);
  void TotCut(const TString& name, Double_t min, Double_t max,
              Bool_t adopt_nan=true);
  void DeCut(const TString& name, Double_t min, Double_t max);

  const HodoCluster* GetTime0BH2Cluster() const;
  const HodoCluster* GetBtof0BH1Cluster() const;
  Double_t Time0() const;
  Double_t Btof0() const;
  Double_t Time0Seg() const;
  Double_t Btof0Seg() const;

private:
  void ClearBH1Hits();
  void ClearBH2Hits();
  void ClearBACHits();
  void ClearTOFHits();
  void ClearLACHits();
  void ClearWCHits();
  void ClearWCSUMHits();
  void ClearBFTHits();

  template <typename T>
  void TimeCut(std::vector<T>& cont, Double_t min, Double_t max);
  template <typename T>
  void TotCut(std::vector<T>& cont,
              Double_t min, Double_t max, Bool_t adopt_nan);
  template <typename T>
  void WidthCutR(std::vector<T>& cont,
                 Double_t min, Double_t max, Bool_t adopt_nan);
  template <typename T>
  void AdcCut(std::vector<T>& cont, Double_t min, Double_t max);
  template <typename T>
  void DeCut(std::vector<T>& cont, Double_t min, Double_t max);

  template <typename T>
  static Int_t MakeUpClusters(const std::vector<T*>& HitCont,
                              HodoCC& ClusterCont,
                              Int_t MaxClusterSize,
                              Double_t MaxTimeDiff);
  template <typename T>
  static HodoCluster* AllocateCluster(HodoHC& HitCont,
                                      index_t index);
  static Bool_t Connectable(const HodoHit* hitA, Int_t indexA,
                            const HodoHit* hitB, Int_t indexB,
                            Double_t MaxTimeDiff);
  static Bool_t Connectable(const FiberHit* hitA, Int_t indexA,
                            const FiberHit* hitB, Int_t indexB,
                            Double_t MaxTimeDiff);

  // static Int_t MakeUpClusters(const BH2HitContainer& HitCont,
  //                             BH2ClusterContainer& ClusterCont,
  //                             Double_t maxTimeDif);
  // static Int_t MakeUpClusters(const FiberHitContainer& cont,
  //                             FiberClusterContainer& ClusterCont,
  //                             Double_t maxTimeDif,
  //                             Int_t DifPairId);
};

//_____________________________________________________________________________
inline TString
HodoAnalyzer::ClassName()
{
  static TString s_name("HodoAnalyzer");
  return s_name;
}

//_____________________________________________________________________________
template <typename T>
inline Bool_t
HodoAnalyzer::DecodeHits(const TString& name, Double_t max_time_diff)
{
  std::vector<T*> CandCont;
  for(auto& rhit: m_raw_data->GetHodoRawHitContainer(name)){
    if(!rhit) continue;
    auto hit = new T(rhit);
    if(hit && hit->Calculate()){
      CandCont.push_back(hit);
    }else{
      delete hit;
    }
  }
  std::sort(CandCont.begin(), CandCont.end(), T::Compare);

  // m_hodo_hit_collection[name] = CandCont;
  auto& cont = m_hodo_hit_collection[name];
  for(auto& hit: cont)
    delete hit;
  cont.clear();
  for(const auto& hit: CandCont)
    cont.push_back(hit);

#if 1
  static const auto& gUser = UserParamMan::GetInstance();
  const auto MaxClusterSize = gUser.Get("MaxClusterSize"+name);
  const auto MaxTimeDiff    = gUser.Get("MaxTimeDiff"+name);
  MakeUpClusters<T>(CandCont, m_hodo_cluster_collection[name],
                    MaxClusterSize, MaxTimeDiff);
#endif

  return true;
}

//_____________________________________________________________________________
template <typename T>
inline Int_t
HodoAnalyzer::MakeUpClusters(const std::vector<T*>& HitCont,
                             HodoCC& ClusterCont,
                             Int_t MaxClusterSize,
                             Double_t MaxTimeDiff)
{
  del::ClearContainer(ClusterCont);

  if(HitCont.empty())
    return 0;
  for(Int_t iA=0, n=HitCont.size(); iA<n; ++iA){
    const auto& hitA = HitCont[iA];
    if(hitA->IsClusteredAll())
      continue;
    T* hitLast = hitA;
    for(Int_t jA=0, mA=hitA->GetEntries(); jA<mA; ++jA){
      if(hitA->IsClustered(jA))
	continue;
      Int_t jLast = jA;
      HodoHC CandCont;
      index_t index;
      hitA->JoinCluster(jA);
      CandCont.push_back(hitA);
      index.push_back(jA);
      for(Int_t iB=iA+1; iB<n; ++iB){
        if(CandCont.size() == MaxClusterSize)
          break;
        const auto& hitB = HitCont[iB];
        for(Int_t jB=0, mB=hitB->GetEntries(); jB<mB; ++jB){
          if(hitB->IsClustered(jB))
	    continue;
          if(Connectable(hitLast, jLast, hitB, jB, MaxTimeDiff)){
            hitB->JoinCluster(jB);
            CandCont.push_back(hitB);
            index.push_back(jB);
            hitLast = hitB;
            jLast = jB;
            break;
          }
        }
      }
      auto cluster = AllocateCluster<T>(CandCont, index);
      if(cluster && cluster->IsGood()){
        ClusterCont.push_back(cluster);
      }else{
        delete cluster;
      }
    }
  }
  return ClusterCont.size();
}

//_____________________________________________________________________________
template <typename T>
inline HodoCluster*
HodoAnalyzer::AllocateCluster(HodoHC& HitCont,
                              index_t index)
{
  return new HodoCluster(HitCont, index);
}

//_____________________________________________________________________________
template <>
inline HodoCluster*
HodoAnalyzer::AllocateCluster<FiberHit>(HodoHC& HitCont,
                                        index_t index)
{
  return new FiberCluster(HitCont, index);
}

//_____________________________________________________________________________
inline Bool_t
HodoAnalyzer::Connectable(const HodoHit* hitA, Int_t indexA,
                          const HodoHit* hitB, Int_t indexB,
                          Double_t MaxTimeDiff)
{
  Double_t cmtA = hitA->CMeanTime(indexA);
  Double_t cmtB = hitB->CMeanTime(indexB);
  Double_t segA = hitA->SegmentId();
  Double_t segB = hitB->SegmentId();
  return (true
          && TMath::Abs(segA - segB) <= 1
          && TMath::Abs(cmtA - cmtB) < MaxTimeDiff);
}

//_____________________________________________________________________________
inline Bool_t
HodoAnalyzer::Connectable(const FiberHit* hitA, Int_t indexA,
                          const FiberHit* hitB, Int_t indexB,
                          Double_t MaxTimeDiff)
{
  TString planeA = hitA->PlaneName();
  if(planeA.Contains("P", TString::kIgnoreCase))
    planeA.ReplaceAll("P", "");
  TString planeB = hitB->PlaneName();
  if(planeB.Contains("P", TString::kIgnoreCase))
    planeB.ReplaceAll("P", "");
  Double_t cmtA = hitA->CMeanTime(indexA);
  Double_t cmtB = hitB->CMeanTime(indexB);
  Double_t posA = hitA->Position();
  Double_t posB = hitB->Position();
#if 0
  hddaq::cout << FUNC_NAME << << std::endl
              << " " << planeA << " == " << planeB << std::endl
              << " " << cmtA << " - " << cmtB << " <= " << MaxTimeDiff
              << std::endl
              << " " << posA << " - " << posB << " <= " << hitA->dXdW()
              << std::endl;
#endif
  return (true
          && planeA.EqualTo(planeB, TString::kIgnoreCase)
          && TMath::Abs(cmtA - cmtB) <= MaxTimeDiff
          && TMath::Abs(posA - posB) <= hitA->dXdW()
    );
}

#endif

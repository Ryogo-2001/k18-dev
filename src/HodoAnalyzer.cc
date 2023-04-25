// -*- C++ -*-

#include "HodoAnalyzer.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "FLHit.hh"
#include "FuncName.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoCluster.hh"
#include "RawData.hh"
#include "UserParamMan.hh"

#define Cluster 1
#define REQDE   0

namespace
{
const auto& gUser = UserParamMan::GetInstance();
const Double_t MaxTimeDifBH1   =  2.0;
const Double_t MaxTimeDifBH2   =  2.0;
const Double_t MaxTimeDifBAC   = -1.0;
const Double_t MaxTimeDifHTOF  = -1.0;
const Double_t MaxTimeDifBVH   = -1.0;
const Double_t MaxTimeDifTOF   = -1.0;
const Double_t MaxTimeDifLAC   = -1.0;
const Double_t MaxTimeDifWC    = -1.0;
const Double_t MaxTimeDifWCSUM = -1.0;
const Double_t MaxTimeDifBFT   =  8.0;
const Double_t MaxTimeDifSCH   = -1.0; //10.0;
const Int_t    MaxSizeCl       = 8;
}

//_____________________________________________________________________________
HodoAnalyzer::HodoAnalyzer()
  : m_BH1Cont(),
    m_BH2Cont(),
    m_BACCont(),
    m_HTOFCont(),
    m_BVHCont(),
    m_TOFCont(),
    m_LACCont(),
    m_WCCont(),
    m_WCSUMCont(),
    m_BFTCont(),
    m_SCHCont(),
    m_BH1ClCont(),
    m_BH2ClCont(),
    m_BACClCont(),
    m_HTOFClCont(),
    m_BVHClCont(),
    m_TOFClCont(),
    m_LACClCont(),
    m_WCClCont(),
    m_WCSUMClCont(),
    m_BFTClCont(),
    m_SCHClCont()
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
HodoAnalyzer::~HodoAnalyzer()
{
  ClearBH1Hits();
  ClearBH2Hits();
  ClearBACHits();
  ClearHTOFHits();
  ClearBVHHits();
  ClearTOFHits();
  ClearLACHits();
  ClearWCHits();
  ClearWCSUMHits();
  ClearBFTHits();
  ClearSCHHits();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
HodoAnalyzer::ClearBH1Hits()
{
  del::ClearContainer(m_BH1Cont);
  del::ClearContainer(m_BH1ClCont);
}

//_____________________________________________________________________________
void
HodoAnalyzer::ClearBH2Hits()
{
  del::ClearContainer(m_BH2Cont);
  del::ClearContainer(m_BH2ClCont);
}

//_____________________________________________________________________________
void
HodoAnalyzer::ClearBACHits()
{
  del::ClearContainer(m_BACCont);
  del::ClearContainer(m_BACClCont);
}

//_____________________________________________________________________________
void
HodoAnalyzer::ClearHTOFHits()
{
  del::ClearContainer(m_HTOFCont);
  del::ClearContainer(m_HTOFClCont);
}

//_____________________________________________________________________________
void
HodoAnalyzer::ClearBVHHits()
{
  del::ClearContainer(m_BVHCont);
  del::ClearContainer(m_BVHClCont);
}

//_____________________________________________________________________________
void
HodoAnalyzer::ClearTOFHits()
{
  del::ClearContainer(m_TOFCont);
  del::ClearContainer(m_TOFClCont);
}

//_____________________________________________________________________________
void
HodoAnalyzer::ClearLACHits()
{
  del::ClearContainer(m_LACCont);
  del::ClearContainer(m_LACClCont);
}

//_____________________________________________________________________________
void
HodoAnalyzer::ClearWCHits()
{
  del::ClearContainer(m_WCCont);
  del::ClearContainer(m_WCClCont);
}

//_____________________________________________________________________________
void
HodoAnalyzer::ClearWCSUMHits()
{
  del::ClearContainer(m_WCSUMCont);
  del::ClearContainer(m_WCSUMClCont);
}

//_____________________________________________________________________________
void
HodoAnalyzer::ClearBFTHits()
{
  for(auto itr = m_BFTCont.begin(); itr != m_BFTCont.end(); ++itr){
    del::ClearContainer(*itr);
  }
  del::ClearContainer(m_BFTClCont);
}

//_____________________________________________________________________________
void
HodoAnalyzer::ClearSCHHits()
{
  del::ClearContainer(m_SCHCont);
  del::ClearContainer(m_SCHClCont);
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::DecodeRawHits(RawData *rawData)
{
  DecodeBH1Hits(rawData);
  DecodeBH2Hits(rawData);
  DecodeBACHits(rawData);
  DecodeHTOFHits(rawData);
  DecodeBVHHits(rawData);
  DecodeTOFHits(rawData);
  DecodeLACHits(rawData);
  DecodeWCHits(rawData);
  DecodeWCSUMHits(rawData);
  DecodeBFTHits(rawData);
  DecodeSCHHits(rawData);
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::DecodeBH1Hits(RawData *rawData)
{
  ClearBH1Hits();
  for(auto& hit: rawData->GetBH1RawHC()){
    if(!hit) continue;
    auto hp = new Hodo2Hit(hit);
    if(hp && hp->Calculate()){
      m_BH1Cont.push_back(hp);
    }else{
      delete hp;
    }
  }

#if Cluster
  MakeUpClusters(m_BH1Cont, m_BH1ClCont, MaxTimeDifBH1);
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::DecodeBH2Hits(RawData *rawData)
{
  ClearBH2Hits();
  for(auto& hit: rawData->GetBH2RawHC()){
    if(!hit) continue;
    auto hp = new BH2Hit(hit);
    if(hp && hp->Calculate()){
      m_BH2Cont.push_back(hp);
    }else{
      delete hp;
    }
  }

#if Cluster
  MakeUpClusters(m_BH2Cont, m_BH2ClCont, MaxTimeDifBH2);
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::DecodeBACHits(RawData *rawData)
{
  ClearBACHits();
  for(auto& hit: rawData->GetBACRawHC()){
    if(!hit) continue;
    auto hp = new Hodo1Hit(hit);
    if(hp && hp->Calculate()){
      m_BACCont.push_back(hp);
    }else{
      delete hp;
    }
  }

#if Cluster
  MakeUpClusters(m_BACCont, m_BACClCont, MaxTimeDifBAC);
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::DecodeHTOFHits(RawData *rawData)
{
  ClearHTOFHits();
  for(auto& hit: rawData->GetHTOFRawHC()){
    if(!hit) continue;
    auto hp = new Hodo2Hit(hit);
    hp->MakeAsTof();
    if(hp && hp->Calculate()){
      m_HTOFCont.push_back(hp);
    }else{
      delete hp;
    }
  }

#if Cluster
  MakeUpClusters(m_HTOFCont, m_HTOFClCont, MaxTimeDifHTOF);
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::DecodeBVHHits(RawData *rawData)
{
  ClearBVHHits();
  for(auto& hit: rawData->GetBVHRawHC()){
    if(!hit) continue;
    auto hp = new Hodo1Hit(hit);
    if(hp && hp->Calculate()){
      m_BVHCont.push_back(hp);
    }else{
      delete hp;
    }
  }

#if Cluster
  MakeUpClusters(m_BVHCont, m_BVHClCont, MaxTimeDifBVH);
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::DecodeTOFHits(RawData *rawData)
{
  ClearTOFHits();
  for(auto& hit: rawData->GetTOFRawHC()){
    if(!hit) continue;
    auto hp = new Hodo2Hit(hit, 15.);
    if(!hp) continue;
    hp->MakeAsTof();
    if(hp->Calculate()){
      m_TOFCont.push_back(hp);
    }else{
      delete hp;
    }
  }

#if Cluster
  MakeUpClusters(m_TOFCont, m_TOFClCont, MaxTimeDifTOF);
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::DecodeLACHits(RawData *rawData)
{
  ClearLACHits();
  for(auto& hit: rawData->GetLACRawHC()){
    if(!hit) continue;
    auto hp = new Hodo1Hit(hit);
    if(hp && hp->Calculate()){
      m_LACCont.push_back(hp);
    }else{
      delete hp;
    }
  }

#if Cluster
  MakeUpClusters(m_LACCont, m_LACClCont, MaxTimeDifLAC);
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::DecodeWCHits(RawData *rawData)
{
  ClearWCHits();
  for(auto& hit: rawData->GetWCRawHC()){
    if(!hit) continue;
    auto hp = new Hodo2Hit(hit);
    if(hp && hp->Calculate()){
      m_WCCont.push_back(hp);
    }else{
      delete hp;
    }
  }

#if Cluster
  MakeUpClusters(m_WCCont, m_WCClCont, MaxTimeDifWC);
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::DecodeWCSUMHits(RawData *rawData)
{
  ClearWCSUMHits();
  for(auto& hit: rawData->GetWCSUMRawHC()){
    if(!hit) continue;
    auto hp = new Hodo1Hit(hit);
    if(hp && hp->Calculate()){
      m_WCSUMCont.push_back(hp);
    }else{
      delete hp;
    }
  }

#if Cluster
  MakeUpClusters(m_WCSUMCont, m_WCSUMClCont, MaxTimeDifWCSUM);
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::DecodeBFTHits(RawData* rawData)
{
  ClearBFTHits();

  m_BFTCont.resize(NumOfPlaneBFT);
  for(Int_t p=0; p<NumOfPlaneBFT; ++p){
    for(auto& hit: rawData->GetBFTRawHC(p)){
      if(!hit) continue;
      auto hp = new FiberHit(hit, "BFT");
      if(hp && hp->Calculate()){
	m_BFTCont.at(p).push_back(hp);
      }else{
	delete hp;
      }
    }
    std::sort(m_BFTCont.at(p).begin(), m_BFTCont.at(p).end(),
              FiberHit::CompFiberHit);
  }

#if Cluster
  FiberHitContainer cont_merge(m_BFTCont.at(0));
  cont_merge.reserve(m_BFTCont.at(0).size() + m_BFTCont.at(1).size());
  cont_merge.insert(cont_merge.end(), m_BFTCont.at(1).begin(),
                    m_BFTCont.at(1).end());
  std::sort(cont_merge.begin(), cont_merge.end(), FiberHit::CompFiberHit);
  MakeUpClusters(cont_merge, m_BFTClCont, MaxTimeDifBFT, 3);
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::DecodeSCHHits(RawData* rawData)
{
  ClearSCHHits();

  for(auto& hit: rawData->GetSCHRawHC()){
    if(!hit) continue;
    auto hp = new FiberHit(hit, "SCH");
    if(hp && hp->Calculate()){
      m_SCHCont.push_back(hp);
    }else{
      delete hp;
    }
  }

  std::sort(m_SCHCont.begin(), m_SCHCont.end(), FiberHit::CompFiberHit);

#if Cluster
  MakeUpClusters(m_SCHCont, m_SCHClCont, MaxTimeDifSCH, 1);
#endif

  return true;
}

//_____________________________________________________________________________
Int_t
HodoAnalyzer::MakeUpClusters(const Hodo1HitContainer& HitCont,
                             HodoClusterContainer& ClusterCont,
                             Double_t maxTimeDif)
{
  del::ClearContainer(ClusterCont);

  for(Int_t i=0, nh=HitCont.size(); i<nh; ++i){
    Hodo1Hit* hitA = HitCont[i];
    Int_t     segA = hitA->SegmentId();
    Hodo1Hit* hitB = 0;
    Int_t       iB = -1;
    Int_t       mB = -1;
    Int_t       mC = -1;
    Int_t     segB = -1;
    if(hitA->JoinedAllMhit()) continue;
    for(Int_t ma=0, n_mhitA=hitA->GetNumOfHit(); ma<n_mhitA; ++ma){
      if(hitA->Joined(ma)) continue;
      Double_t cmtA = hitA->CMeanTime(ma);
      Double_t cmtB = -9999;
      for(Int_t j=i+1; j<nh; ++j){
	Hodo1Hit *hit = HitCont[j];
	Int_t    seg  = hit->SegmentId();
	for(Int_t mb=0, n_mhitB = hit->GetNumOfHit(); mb<n_mhitB; ++mb){
	  if(hit->Joined(mb)) continue;
	  Double_t cmt = hit->CMeanTime(mb);
	  if(std::abs(seg-segA)==1 && std::abs(cmt-cmtA)<maxTimeDif){
	    hitB = hit; iB = j; mB = mb; segB = seg; cmtB = cmt; break;
	  }
	}// for(mb:hitB)
      }// for(j:hitB)
      if(hitB){
	Hodo1Hit *hitC = nullptr;
	for(Int_t j=i+1; j<nh; ++j){
	  if(j==iB) continue;
	  Hodo1Hit *hit = HitCont[j];
	  Int_t    seg  = hit->SegmentId();
	  for(Int_t mc=0, n_mhitC=hit->GetNumOfHit(); mc<n_mhitC; ++mc){
	    if(hit->Joined(mc)) continue;
	    Double_t cmt = hit->CMeanTime(mc);
	    if((std::abs(seg-segA)==1 && std::abs(cmt-cmtA)<maxTimeDif) ||
               (std::abs(seg-segB)==1 && std::abs(cmt-cmtB)<maxTimeDif)){
	      hitC=hit; mC=mc; break;
	    }
	  }
	}
	if(hitC){
	  hitA->SetJoined(ma);
	  hitB->SetJoined(mB);
	  hitC->SetJoined(mC);
	  HodoCluster *cluster = new HodoCluster(hitA, hitB, hitC);
	  cluster->SetIndex(ma, mB, mC);
	  cluster->Calculate();
	  if(cluster) ClusterCont.push_back(cluster);
	} else {
	  hitA->SetJoined(ma);
	  hitB->SetJoined(mB);
	  HodoCluster *cluster = new HodoCluster(hitA, hitB);
	  cluster->SetIndex(ma, mB);
	  cluster->Calculate();
	  if(cluster) ClusterCont.push_back(cluster);
	}
      } else {
	hitA->SetJoined(ma);
	HodoCluster *cluster = new HodoCluster(hitA);
	cluster->SetIndex(ma);
	cluster->Calculate();
	if(cluster) ClusterCont.push_back(cluster);
      }
    }// for(ma:hitA)
  }// for(i:hitA)
  return ClusterCont.size();
}


//_____________________________________________________________________________
Int_t
HodoAnalyzer::MakeUpClusters(const Hodo2HitContainer& HitCont,
                             HodoClusterContainer& ClusterCont,
                             Double_t maxTimeDif)
{
  del::ClearContainer(ClusterCont);

  for(Int_t i=0, nh=HitCont.size(); i<nh; ++i){
    Hodo2Hit* hitA = HitCont[i];
    Int_t     segA = hitA->SegmentId();
    Hodo2Hit* hitB = 0;
    Int_t       iB = -1;
    Int_t       mB = -1;
    Int_t       mC = -1;
    Int_t     segB = -1;
    if(hitA->JoinedAllMhit()) continue;
    for(Int_t ma=0, n_mhitA=hitA->GetNumOfHit(); ma<n_mhitA; ++ma){
      if(hitA->Joined(ma)) continue;
      Double_t    cmtA = hitA->CMeanTime(ma);
      Double_t    cmtB = -9999;
      for(Int_t j=i+1; j<nh; ++j){
	Hodo2Hit* hit = HitCont[j];
	Int_t     seg = hit->SegmentId();
	for(Int_t mb=0, n_mhitB=hit->GetNumOfHit(); mb<n_mhitB; ++mb){
	  if(hit->Joined(mb)) continue;
	  Double_t cmt = hit->CMeanTime(mb);
	  if(std::abs(seg-segA)==1 && std::abs(cmt-cmtA)<maxTimeDif){
	    hitB = hit; iB = j; mB = mb; segB = seg; cmtB = cmt; break;
	  }
	}// for(mb:hitB)
      }// for(j:hitB)
      if(hitB){
	Hodo2Hit *hitC = nullptr;
	for(Int_t j=i+1; j<nh; ++j){
	  if(j==iB) continue;
	  Hodo2Hit *hit = HitCont[j];
	  Int_t    seg  = hit->SegmentId();
	  for(Int_t mc=0, n_mhitC=hit->GetNumOfHit(); mc<n_mhitC; ++mc){
	    if(hit->Joined(mc)) continue;
	    Double_t cmt = hit->CMeanTime(mc);
	    if((std::abs(seg-segA)==1 && std::abs(cmt-cmtA)<maxTimeDif) ||
               (std::abs(seg-segB)==1 && std::abs(cmt-cmtB)<maxTimeDif)){
	      hitC=hit; mC=mc; break;
	    }
	  }
	}
	if(hitC){
	  hitA->SetJoined(ma);
	  hitB->SetJoined(mB);
	  hitC->SetJoined(mC);
	  HodoCluster *cluster = new HodoCluster(hitA, hitB, hitC);
	  cluster->SetIndex(ma, mB, mC);
	  cluster->Calculate();
	  if(cluster) ClusterCont.push_back(cluster);
	} else {
	  hitA->SetJoined(ma);
	  hitB->SetJoined(mB);
	  HodoCluster *cluster = new HodoCluster(hitA, hitB);
	  cluster->SetIndex(ma, mB);
	  cluster->Calculate();
	  if(cluster) ClusterCont.push_back(cluster);
	}
      } else {
	hitA->SetJoined(ma);
	HodoCluster *cluster = new HodoCluster(hitA);
	cluster->SetIndex(ma);
	cluster->Calculate();
	if(cluster) ClusterCont.push_back(cluster);
      }
    }// for(ma:hitA)
  }// for(i:hitA)
  return ClusterCont.size();
}

//_____________________________________________________________________________
Int_t
HodoAnalyzer::MakeUpClusters(const BH2HitContainer& HitCont,
                             BH2ClusterContainer& ClusterCont,
                             Double_t maxTimeDif)
{
  del::ClearContainer(ClusterCont);

  for(Int_t i=0, nh=HitCont.size(); i<nh; ++i){
    BH2Hit* hitA = HitCont[i];
    Int_t   segA = hitA->SegmentId();
    BH2Hit* hitB = 0;
    Int_t     iB = -1;
    Int_t     mB = -1;
    Int_t     mC = -1;
    Int_t   segB = -1;
    if(hitA->JoinedAllMhit()) continue;
    for(Int_t ma=0, n_mhitA=hitA->GetNumOfHit(); ma<n_mhitA; ++ma){
      if(hitA->Joined(ma)) continue;
      Double_t cmtA = hitA->CMeanTime(ma);
      Double_t cmtB = -9999;
      for(int j=i+1; j<nh; ++j){
	BH2Hit *hit = HitCont[j];
	Int_t   seg = hit->SegmentId();
	for(Int_t mb=0, n_mhitB=hit->GetNumOfHit(); mb<n_mhitB; ++mb){
	  if(hit->Joined(mb)) continue;
	  Double_t cmt = hit->CMeanTime(mb);
	  if(std::abs(seg-segA)==1 && std::abs(cmt-cmtA)<maxTimeDif){
	    hitB = hit; iB = j; mB = mb; segB = seg; cmtB = cmt; break;
	  }
	}// for(mb:hitB)
      }// for(j:hitB)
      if(hitB){
	BH2Hit *hitC = nullptr;
	for(Int_t j=i+1; j<nh; ++j){
	  if(j==iB) continue;
	  BH2Hit *hit = HitCont[j];
	  Int_t   seg = hit->SegmentId();
	  for(Int_t mc=0, n_mhitC=hit->GetNumOfHit(); mc<n_mhitC; ++mc){
	    if(hit->Joined(mc)) continue;
	    Double_t cmt = hit->CMeanTime(mc);
	    if((std::abs(seg-segA)==1 && std::abs(cmt-cmtA)<maxTimeDif) ||
               (std::abs(seg-segB)==1 && std::abs(cmt-cmtB)<maxTimeDif)){
	      hitC = hit; mC=mc; break;
	    }
	  }
	}
	if(hitC){
	  hitA->SetJoined(ma);
	  hitB->SetJoined(mB);
	  hitC->SetJoined(mC);
	  BH2Cluster *cluster = new BH2Cluster(hitA, hitB, hitC);
	  cluster->SetIndex(ma, mB, mC);
	  cluster->Calculate();
	  if(cluster) ClusterCont.push_back(cluster);
	} else {
	  hitA->SetJoined(ma);
	  hitB->SetJoined(mB);
	  BH2Cluster *cluster = new BH2Cluster(hitA, hitB);
	  cluster->SetIndex(ma, mB);
	  cluster->Calculate();
	  if(cluster) ClusterCont.push_back(cluster);
	}
      } else {
	hitA->SetJoined(ma);
	BH2Cluster *cluster = new BH2Cluster(hitA);
	cluster->SetIndex(ma);
	cluster->Calculate();
	if(cluster) ClusterCont.push_back(cluster);
      }
    }// for(ma:hitA)
  }// for(i:hitA)

  return ClusterCont.size();
}

//_____________________________________________________________________________
Int_t
HodoAnalyzer::MakeUpClusters(const FiberHitContainer& cont,
                             FiberClusterContainer& ClusterCont,
                             Double_t maxTimeDif,
                             Int_t DifPairId)
{
  del::ClearContainer(ClusterCont);

  for(Int_t seg=0, NofSeg=cont.size(); seg<NofSeg; ++seg){
    FiberHit* HitA = cont.at(seg);
    Bool_t fl_ClCandA = false;
    if(seg != (NofSeg -1)){
      if(DifPairId > (cont.at(seg+1)->PairId() - HitA->PairId())){
	fl_ClCandA = true;
      }
    }

    for(Int_t mhitA=0, NofHitA=HitA->GetNumOfHit(); mhitA<NofHitA; ++mhitA){
      if(HitA->Joined(mhitA)) continue;
      FiberCluster *cluster = new FiberCluster;
      cluster->push_back(new FLHit(HitA, mhitA));
      if(!fl_ClCandA){
	// there is no more candidates
	if(cluster->Calculate()){
	  ClusterCont.push_back(cluster);
	} else {
	  delete cluster;
	  cluster = nullptr;
	}
	continue;
      }
      // Start Search HitB
      Double_t cmtA    = (Double_t)HitA->GetCTime(mhitA);
      Int_t    NofHitB = cont.at(seg+1)->GetNumOfHit();
      Bool_t   fl_HitB = false;
      Double_t cmtB    = -1;
      Int_t    CurrentPair = HitA->PairId();
      for(Int_t mhitB = 0; mhitB<NofHitB; ++mhitB){
	if(cont.at(seg+1)->Joined(mhitB)) continue;
	FiberHit* HitB = cont.at(seg+1);
	cmtB = (Double_t)HitB->GetCTime(mhitB);
	if(std::abs(cmtB-cmtA)<maxTimeDif){
	  cluster->push_back(new FLHit(HitB, mhitB));
	  CurrentPair = HitB->PairId();
	  fl_HitB = true;
	  break;
	}
      }
      Bool_t fl_ClCandB  = false;
      if((seg+1) != (NofSeg -1)){
	if(DifPairId > (cont.at(seg+2)->PairId() - CurrentPair)){
	  fl_ClCandB = true;
	}
      }
      if(!fl_ClCandB){
	// there is no more candidates
	if(cluster->Calculate()){
	  ClusterCont.push_back(cluster);
	} else {
	  delete cluster;
	  cluster = nullptr;
	}
	continue;
      }
      // Start Search HitC
      Bool_t   fl_HitC = false;
      Double_t cmtC    = -1;
      for(Int_t mhitC=0, NofHitC=cont.at(seg+2)->GetNumOfHit();
          mhitC<NofHitC; ++mhitC){
	if(cont.at(seg+2)->Joined(mhitC)) continue;
	FiberHit* HitC = cont.at(seg+2);
	cmtC = (Double_t)HitC->GetCTime(mhitC);
	if(true
           && std::abs(cmtC-cmtA)<maxTimeDif
           && !(fl_HitB && (std::abs(cmtC-cmtB)>maxTimeDif))
          ){
	  cluster->push_back(new FLHit(HitC, mhitC));
	  CurrentPair = HitC->PairId();
	  fl_HitC = true;
	  break;
	}
      }
      Bool_t fl_ClCandC = false;
      if((seg+2) != (NofSeg -1)){
	if(DifPairId > (cont.at(seg+3)->PairId() - CurrentPair)){
	  fl_ClCandC = true;
	}
      }
      if(!fl_ClCandC){
	// there is no more candidates
	if(cluster->Calculate()){
	  ClusterCont.push_back(cluster);
	} else {
	  delete cluster;
	  cluster = nullptr;
	}
	continue;
      }
      // Start Search HitD
      Double_t cmtD = -1;
      for(Int_t mhitD=0, NofHitD=cont.at(seg+3)->GetNumOfHit();
          mhitD<NofHitD; ++mhitD){
	if(cont.at(seg+3)->Joined(mhitD)) continue;
	FiberHit* HitD = cont.at(seg+3);
	cmtD = (Double_t)HitD->GetCTime(mhitD);
	if(true
           && std::abs(cmtD-cmtA)<maxTimeDif
           && !(fl_HitB && (std::abs(cmtD-cmtB)>maxTimeDif))
           && !(fl_HitC && (std::abs(cmtD-cmtC)>maxTimeDif))
          ){
	  cluster->push_back(new FLHit(HitD, mhitD));
	  break;
	}
      }
      // Finish
      if(cluster->Calculate()){
	ClusterCont.push_back(cluster);
      } else {
	delete cluster;
	cluster = nullptr;
      }
    }
  }
  return ClusterCont.size();
}

//_____________________________________________________________________________
// Clustering function for FBH
Int_t
HodoAnalyzer::MakeUpClusters(const FLHitContainer& cont,
                             FiberClusterContainer& ClusterCont,
                             Double_t maxTimeDif,
                             Int_t DifPairId)
{
  std::vector<FiberCluster*> DeleteCand;
  for(Int_t seg=0, NofSeg=cont.size(); seg<NofSeg; ++seg){
    FLHit* HitA = cont.at(seg);
    if(HitA->Joined()) continue;
    HitA->SetJoined();
    FiberCluster *cluster = new FiberCluster;
    cluster->push_back(HitA);
    Int_t CurrentPair = HitA->PairId();
    for(Int_t smarker = seg+1; smarker<NofSeg; ++smarker){
      FLHit* HitB = cont.at(smarker);
      if(false
         || DifPairId < (HitB->PairId() - CurrentPair)
         || HitB->PairId() == CurrentPair
        ){ continue; }
      if(HitB->Joined()) continue;
      HitB->SetJoined();
      Double_t cmtB = HitB->GetCTime();
      Bool_t   fl_all_valid = true;
      for(Int_t cindex = 0; cindex<cluster->VectorSize(); ++cindex){
	Double_t cmt = cluster->GetHit(cindex)->GetCTime();
	if(std::abs(cmt-cmtB) > maxTimeDif){
	  fl_all_valid = false; break;
	}
      }
      if(fl_all_valid){
	// we found a cluster candidate
	cluster->push_back(HitB);
	CurrentPair = HitB->PairId();
	break;
      }
    }
    // Finish
    if(cluster->Calculate()){
      ClusterCont.push_back(cluster);
    } else {
      DeleteCand.push_back(cluster);
    }
  }
  del::ClearContainer(DeleteCand);
  return ClusterCont.size();
}

//_____________________________________________________________________________
Int_t
HodoAnalyzer::MakeUpCoincidence(const FiberHitContainer& cont,
                                FLHitContainer& CoinCont,
                                Double_t maxTimeDif)
{
  for(Int_t seg=0, NofSeg=cont.size(); seg<NofSeg; ++seg){
    FiberHit* HitA = cont.at(seg);
    Int_t NofHitA = HitA->GetNumOfHit();
    for(Int_t mhitA = 0; mhitA<NofHitA; ++mhitA){
      if(HitA->Joined(mhitA)) continue;
      Double_t cmt = HitA->GetCTime(mhitA);
      Int_t CurrentPair = HitA->PairId();
      for(Int_t smarker = seg+1; smarker<NofSeg; ++smarker){
	FiberHit* HitB = cont.at(smarker);
	if(0 != (HitB->PairId() - CurrentPair)) continue;
	Int_t NofHitB = HitB->GetNumOfHit();
	for(Int_t mhitB = 0; mhitB<NofHitB; ++mhitB){
	  if(HitB->Joined(mhitB)) continue;
	  Double_t cmtB = HitB->GetCTime(mhitB);
	  Bool_t   fl_all_valid = true;
	  if(std::abs(cmt-cmtB)>maxTimeDif){
	    fl_all_valid = false;
	    break;
	  }
	  if(fl_all_valid){
	    // we found a coin candidate
	    //	    hddaq::cout << "yes" << std::endl;
	    CoinCont.push_back(new FLHit(HitA, HitB, mhitA, mhitB));
	    break;
	  }
	}// for(mhitB)
      }// for(segB)
    }// for(mhitA)
  }// for(segA)

  return CoinCont.size();
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcBH1Hits(Bool_t applyRecursively)
{
  for(Int_t i=0, n=m_BH1Cont.size(); i<n; ++i){
    Hodo2Hit *hit = m_BH1Cont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcBACHits(Bool_t applyRecursively)
{
  for(Int_t i=0, n=m_BACCont.size(); i<n; ++i){
    Hodo1Hit *hit = m_BACCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcBH2Hits(Bool_t applyRecursively)
{
  for(Int_t i=0, n=m_BH2Cont.size(); i<n; ++i){
    BH2Hit *hit = m_BH2Cont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcHTOFHits(Bool_t applyRecursively)
{
  for(auto& hit: m_HTOFCont){
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcBVHHits(Bool_t applyRecursively)
{
  for(auto& hit: m_BVHCont){
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcTOFHits(Bool_t applyRecursively)
{
  for(Int_t i=0, n=m_TOFCont.size(); i<n; ++i){
    Hodo2Hit *hit = m_TOFCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;

}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcLACHits(Bool_t applyRecursively)
{
  for(Int_t i=0, n=m_LACCont.size(); i<n; ++i){
    Hodo1Hit *hit = m_LACCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcWCHits(Bool_t applyRecursively)
{
  for(auto& hit: m_WCCont){
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcWCSUMHits( Bool_t applyRecursively )
{
  for(auto& hit: m_WCSUMCont){
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcBH1Clusters(Bool_t applyRecursively)
{
  for(Int_t i=0, n=m_BH1ClCont.size(); i<n; ++i){
    HodoCluster *cl = m_BH1ClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcBACClusters(Bool_t applyRecursively)
{
  for(Int_t i=0, n=m_BACClCont.size(); i<n; ++i){
    HodoCluster *cl = m_BACClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcBH2Clusters(Bool_t applyRecursively)
{
  for(Int_t i=0, n=m_BH2ClCont.size(); i<n; ++i){
    BH2Cluster *cl = m_BH2ClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcHTOFClusters(Bool_t applyRecursively)
{
  for(auto& cl: m_HTOFClCont){
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcBVHClusters(Bool_t applyRecursively)
{
  for(auto& cl: m_BVHClCont){
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcTOFClusters(Bool_t applyRecursively)
{
  for(Int_t i=0, n=m_TOFClCont.size(); i<n; ++i){
    HodoCluster *cl = m_TOFClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcLACClusters(Bool_t applyRecursively)
{
  for(Int_t i=0, n=m_LACClCont.size(); i<n; ++i){
    HodoCluster *cl = m_LACClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcWCSUMClusters(Bool_t applyRecursively)
{
  for(auto& cl: m_WCSUMClCont){
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcAll()
{
  ReCalcBH1Hits();
  ReCalcBACHits();
  ReCalcBH2Hits();
  ReCalcHTOFHits();
  ReCalcBVHHits();
  ReCalcTOFHits();
  ReCalcLACHits();
  ReCalcWCHits();
  ReCalcWCSUMHits();
  ReCalcBH1Clusters();
  ReCalcBACClusters();
  ReCalcBH2Clusters();
  ReCalcHTOFClusters();
  ReCalcBVHClusters();
  ReCalcTOFClusters();
  ReCalcLACClusters();
  ReCalcWCSUMClusters();
  return true;
}

//_____________________________________________________________________________
void
HodoAnalyzer::TimeCutBH1(Double_t tmin, Double_t tmax)
{
  TimeCut(m_BH1ClCont, tmin, tmax);
}

//_____________________________________________________________________________
void
HodoAnalyzer::TimeCutBH2(Double_t tmin, Double_t tmax)
{
  TimeCut(m_BH2ClCont, tmin, tmax);
}

//_____________________________________________________________________________
void
HodoAnalyzer::TimeCutHTOF(Double_t tmin, Double_t tmax)
{
  TimeCut(m_HTOFClCont, tmin, tmax);
}

//_____________________________________________________________________________
void
HodoAnalyzer::TimeCutBVH(Double_t tmin, Double_t tmax)
{
  TimeCut(m_BVHClCont, tmin, tmax);
}

//_____________________________________________________________________________
void
HodoAnalyzer::TimeCutTOF(Double_t tmin, Double_t tmax)
{
  TimeCut(m_TOFClCont, tmin, tmax);
}

//_____________________________________________________________________________
void
HodoAnalyzer::TimeCutBFT(Double_t tmin, Double_t tmax)
{
  TimeCut(m_BFTClCont, tmin, tmax);
}

//_____________________________________________________________________________
void
HodoAnalyzer::TimeCutSCH(Double_t tmin, Double_t tmax)
{
  TimeCut(m_SCHClCont, tmin, tmax);
}

//_____________________________________________________________________________
// Implementation of Time cut for the cluster container
template <typename TypeCluster>
void
HodoAnalyzer::TimeCut(std::vector<TypeCluster>& cont,
                      Double_t tmin, Double_t tmax)
{
  std::vector<TypeCluster> DeleteCand;
  std::vector<TypeCluster> ValidCand;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    Double_t ctime = cont.at(i)->CMeanTime();
    if(tmin < ctime && ctime < tmax){
      ValidCand.push_back(cont.at(i));
    }else{
      DeleteCand.push_back(cont.at(i));
    }
  }

  del::ClearContainer(DeleteCand);

  cont.clear();
  cont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
  ValidCand.clear();
}

//_____________________________________________________________________________
void
HodoAnalyzer::WidthCutBFT(Double_t min_width, Double_t max_width)
{
  WidthCut(m_BFTClCont, min_width, max_width , true);
}

//_____________________________________________________________________________
void
HodoAnalyzer::WidthCutSCH(Double_t min_width, Double_t max_width)
{
  WidthCut(m_SCHClCont, min_width, max_width , true);
}

//_____________________________________________________________________________
// Implementation of width cut for the cluster container
template <typename TypeCluster>
void
HodoAnalyzer::WidthCut(std::vector<TypeCluster>& cont,
                       Double_t min_width, Double_t max_width,
                       Bool_t adopt_nan)
{
  std::vector<TypeCluster> DeleteCand;
  std::vector<TypeCluster> ValidCand;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    Double_t width = cont.at(i)->Width();
    if(isnan(width) && adopt_nan){
      ValidCand.push_back(cont.at(i));
    }else if(min_width < width && width < max_width){
      ValidCand.push_back(cont.at(i));
    }else{
      DeleteCand.push_back(cont.at(i));
    }
  }

  del::ClearContainer(DeleteCand);

  cont.clear();
  cont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
  ValidCand.clear();
}

//_____________________________________________________________________________
// Implementation of width cut for the cluster container
template <typename TypeCluster>
void
HodoAnalyzer::WidthCutR(std::vector<TypeCluster>& cont,
                        Double_t min_width, Double_t max_width,
                        Bool_t adopt_nan)
{
  std::vector<TypeCluster> DeleteCand;
  std::vector<TypeCluster> ValidCand;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    Double_t width = cont.at(i)->Width();
    //Double_t width = -cont.at(i)->minWidth();//reverse

    if(isnan(width) && adopt_nan){
      ValidCand.push_back(cont.at(i));
    }else if(min_width < width && width < max_width){
      ValidCand.push_back(cont.at(i));
    }else{
      DeleteCand.push_back(cont.at(i));
    }
  }

  del::ClearContainer(DeleteCand);

  cont.clear();
  cont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
  ValidCand.clear();
}

//_____________________________________________________________________________
// Implementation of ADC cut for the cluster container
template <typename TypeCluster>
void
HodoAnalyzer::AdcCut(std::vector<TypeCluster>& cont,
                     Double_t amin, Double_t amax)
{
  std::vector<TypeCluster> DeleteCand;
  std::vector<TypeCluster> ValidCand;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    Double_t adc = cont.at(i)->MaxAdcLow();
    if(amin < adc && adc < amax){
      ValidCand.push_back(cont.at(i));
    }else{
      DeleteCand.push_back(cont.at(i));
    }
  }

  del::ClearContainer(DeleteCand);

  cont.clear();
  cont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
  ValidCand.clear();
}

//_____________________________________________________________________________
BH2Cluster*
HodoAnalyzer::GetTime0BH2Cluster()
{
  static const Double_t MinMt = gUser.GetParameter("MtBH2", 0);
  static const Double_t MaxMt = gUser.GetParameter("MtBH2", 1);
#if REQDE
  static const Double_t MinDe = gUser.GetParameter("DeBH2", 0);
  static const Double_t MaxDe = gUser.GetParameter("DeBH2", 1);
#endif

  BH2Cluster* time0_cluster = nullptr;
  Double_t min_mt = -9999;
  for(const auto& cluster : m_BH2ClCont){
    Double_t mt = cluster->MeanTime();
    if(true
       && std::abs(mt) < std::abs(min_mt)
       && MinMt < mt && mt < MaxMt
#if REQDE
       && (MinDe < cluster->DeltaE() && cluster->DeltaE() < MaxDe)
#endif
      ){
      min_mt        = mt;
      time0_cluster = cluster;
    }// T0 selection
  }// for

  return time0_cluster;
}

//_____________________________________________________________________________
HodoCluster*
HodoAnalyzer::GetBtof0BH1Cluster(Double_t time0)
{
  static const Double_t MinBtof = gUser.GetParameter("BTOF",  0);
  static const Double_t MaxBtof = gUser.GetParameter("BTOF",  1);
#if REQDE
  static const Double_t MinDe   = gUser.GetParameter("DeBH1", 0);
  static const Double_t MaxDe   = gUser.GetParameter("DeBH1", 1);
#endif

  HodoCluster* time0_cluster = nullptr;
  Double_t min_btof            = -9999;
  for(const auto& cluster : m_BH1ClCont){
    Double_t cmt  = cluster->CMeanTime();
    Double_t btof = cmt - time0;
    if(true
       && std::abs(btof) < std::abs(min_btof)
       && MinBtof < btof && btof < MaxBtof
#if REQDE
       && (MinDe < cluster->DeltaE() && cluster->DeltaE() < MaxDe)
#endif
      ){
      min_btof      = btof;
      time0_cluster = cluster;
    }// T0 selection
  }// for

  return time0_cluster;
}

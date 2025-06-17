// -*- C++ -*-

#ifndef TPC_ANALYZER_HH
#define TPC_ANALYZER_HH

#include <vector>
#include <TString.h>
#include <TVector3.h>

#include "DetectorID.hh"

class RawData;
class TPCHit;
class TPCCluster;
class TPCLocalTrack;
class TPCLocalTrackHelix;
class TPCVertex;
class TPCRKTrack;

typedef std::vector<TPCHit*>        TPCHitContainer;
typedef std::vector<TPCCluster*>    TPCClusterContainer;
typedef std::vector<TPCLocalTrack*> TPCLocalTrackContainer;
typedef std::vector<TPCLocalTrackHelix*> TPCLocalTrackHelixContainer;
typedef std::vector<TPCVertex*> TPCVertexContainer;
typedef std::vector<TPCRKTrack*>    TPCRKTrackContainer;

//_____________________________________________________________________________
class TPCAnalyzer
{
public:
  static const TString& ClassName();
  TPCAnalyzer();
  ~TPCAnalyzer();

private:
  TPCAnalyzer(const TPCAnalyzer&);
  TPCAnalyzer& operator =(const TPCAnalyzer&);

private:
  enum e_type
  { kTPC, kTPCTracking, kTPCK18, kTPCKurama, kKKvertex, n_type };
  std::vector<Bool_t>                m_is_decoded;
  std::vector<TPCHitContainer>       m_TPCHitCont;
  std::vector<TPCClusterContainer>   m_TPCClCont;
  TPCLocalTrackContainer             m_TPCTC;
  TPCLocalTrackContainer             m_TPCTCFailed;
  TPCLocalTrackHelixContainer        m_TPCTCHelix;
  TPCLocalTrackHelixContainer        m_TPCTCHelixInverted;
  TPCLocalTrackHelixContainer        m_TPCTCVP; //Reconstucted K1.8, Kurama RK track's hits on virtual plane
  TPCLocalTrackHelixContainer        m_TPCTCHelixFailed;
  TPCRKTrackContainer                m_TPCKuramaTC;
  TPCRKTrackContainer                m_TPCK18TC;
  TPCVertexContainer                 m_TPCVC; //vertex between two tracks
  TPCVertexContainer                 m_TPCVCClustered; //clusted position of multi-tracks

public:

  //TPC Hit&Cluster
  Bool_t DecodeTPCHitsGeant4(const Int_t nhits,
                             const Double_t *x, const Double_t *y,
                             const Double_t *z, const Double_t *de,
			     const Int_t *pid, std::vector<TVector3> Mom = std::vector<TVector3>());
  Double_t GetDetectionEfficiency(TVector3 pos, Int_t pid, TVector3 mom = TVector3(0,0,0), Double_t de=0);
  Int_t GetClusterSize(TVector3 pos, Int_t pid, TVector3 mom = TVector3(0,0,0),Double_t de=0);
  Bool_t DecodeTPCHits(RawData* rawData, Double_t clock=0.);
  Bool_t ReCalcTPCHits(const Int_t nhits,
                       const std::vector<Int_t>& pad,
                       const std::vector<Double_t>& time,
                       const std::vector<Double_t>& de,
                       Double_t clock=0.);
  Bool_t ReCalcTPCTracks(const Int_t ntracks,
			 const std::vector<Int_t>& nhits,
			 const std::vector<Double_t>& x0,
			 const std::vector<Double_t>& y0,
			 const std::vector<Double_t>& u0,
			 const std::vector<Double_t>& v0,
			 const std::vector<std::vector<Double_t>>& layer,
			 const std::vector<std::vector<Double_t>>& mrow,
			 const std::vector<std::vector<Double_t>>& de,
			 const std::vector<std::vector<Double_t>>& res_x,
			 const std::vector<std::vector<Double_t>>& res_y,
			 const std::vector<std::vector<Double_t>>& res_z,
			 const std::vector<std::vector<Double_t>>& localpos_x,
			 const std::vector<std::vector<Double_t>>& localpos_y,
			 const std::vector<std::vector<Double_t>>& localpos_z);
  Bool_t ReCalcTPCTracks(const Int_t ntracks,
			 const std::vector<Int_t>& isK18,
			 const std::vector<Int_t>& isKurama,
			 const std::vector<Int_t>& charge,
			 const std::vector<Int_t>& nhits,
			 const std::vector<Double_t>& cx,
			 const std::vector<Double_t>& cy,
			 const std::vector<Double_t>& z0,
			 const std::vector<Double_t>& r,
			 const std::vector<Double_t>& dz,
			 const std::vector<std::vector<Double_t>>& layer,
			 const std::vector<std::vector<Double_t>>& mrow,
			 const std::vector<std::vector<Double_t>>& helix_t,
			 const std::vector<std::vector<Double_t>>& de,
			 const std::vector<std::vector<Double_t>>& res_x,
			 const std::vector<std::vector<Double_t>>& res_y,
			 const std::vector<std::vector<Double_t>>& res_z,
			 const std::vector<std::vector<Double_t>>& localpos_x,
			 const std::vector<std::vector<Double_t>>& localpos_y,
			 const std::vector<std::vector<Double_t>>& localpos_z);
  Bool_t ReCalcTPCTracksGeant4(const Int_t ntracks,
			       const std::vector<Int_t>& charge,
			       const std::vector<Int_t>& nhits,
			       const std::vector<Double_t>& cx,
			       const std::vector<Double_t>& cy,
			       const std::vector<Double_t>& z0,
			       const std::vector<Double_t>& r,
			       const std::vector<Double_t>& dz,
			       const std::vector<std::vector<Double_t>>& layer,
			       const std::vector<std::vector<Double_t>>& mrow,
			       const std::vector<std::vector<Double_t>>& helix_t,
			       const std::vector<std::vector<Double_t>>& de,
			       const std::vector<std::vector<Double_t>>& res_x,
			       const std::vector<std::vector<Double_t>>& res_y,
			       const std::vector<std::vector<Double_t>>& res_z,
			       const std::vector<std::vector<Double_t>>& localpos_x,
			       const std::vector<std::vector<Double_t>>& localpos_y,
			       const std::vector<std::vector<Double_t>>& localpos_z);
  const TPCHitContainer& GetTPCHC(Int_t l) const { return m_TPCHitCont.at(l); }
  const TPCClusterContainer& GetTPCClCont(Int_t l) const { return m_TPCClCont.at(l); }

  //HS-Off
  Bool_t TrackSearchTPC(Bool_t exclusive=false);
  Int_t GetNTracksTPC() const { return m_TPCTC.size(); }
  Int_t GetNTracksTPCFailed() const { return m_TPCTCFailed.size(); }
  TPCLocalTrack* GetTrackTPC(Int_t l) const { return m_TPCTC.at(l); }
  TPCLocalTrack* GetTrackTPCFailed(Int_t l) const { return m_TPCTCFailed.at(l); }
  Bool_t TestHoughTransform();

  //HS-On
  Bool_t TrackSearchTPCHelix(Bool_t exclusive=false);
  Bool_t TrackSearchTPCHelix(std::vector<std::vector<TVector3>> K18VPs,
			     std::vector<std::vector<TVector3>> KuramaVPs,
			     std::vector<Double_t> KuramaCharge,
			     Bool_t exclusive=false);
  Int_t GetNTracksTPCHelix() const { return m_TPCTCHelix.size(); }
  Int_t GetNTracksTPCHelixChargeInverted() const { return m_TPCTCHelixInverted.size(); }
  Int_t GetNTracksTPCVP() const { return m_TPCTCVP.size(); }
  Int_t GetNTracksTPCHelixFailed() const { return m_TPCTCHelixFailed.size(); }
  TPCLocalTrackHelix* GetTrackTPCHelix(Int_t l) const { return m_TPCTCHelix.at(l); }
  TPCLocalTrackHelix* GetTrackTPCVP(Int_t l) const { return m_TPCTCVP.at(l); }
  TPCLocalTrackHelix* GetTrackTPCHelixChargeInverted(Int_t l) const { return m_TPCTCHelixInverted.at(l); }
  TPCLocalTrackHelix* GetTrackTPCHelixFailed(Int_t l) const { return m_TPCTCHelixFailed.at(l); }
  Bool_t TestHoughTransformHelix();

  //Kurama Tracking
  void MakeTPCKuramaHPContainer(TPCLocalTrackHelix *tpctrack,
				std::vector<Int_t> &lnum);
  TPCRKTrack* MakeTPCKuramaTrack(TPCLocalTrackHelix *tpctrack,
				 const std::vector<Int_t> layerID,
				 const std::vector<Double_t> wire,
				 const std::vector<Double_t> localHitPos);
  Bool_t DoFitTPCKuramaTrack(Int_t KuramaTrackID,
			     const Int_t PID,
			     const TVector3 initPos,
			     const TVector3 initMom,
			     const std::vector<Int_t> layerID,
			     const std::vector<Double_t> wire,
			     const std::vector<Double_t> localHitPos,
			     Bool_t U2D);
  Bool_t TrackSearchTPCKurama(const std::vector<Int_t> PID,
			      const std::vector<TVector3> initPos,
			      const std::vector<TVector3> initMom,
			      const std::vector<std::vector<Int_t>> layerID,
			      const std::vector<std::vector<Double_t>> wire,
			      const std::vector<std::vector<Double_t>> localHitPos,
			      Bool_t U2D=false);
  Int_t GetNTracksTPCKurama() const { return m_TPCKuramaTC.size(); }
  TPCRKTrack* GetTPCKuramaTrack(Int_t l) const { return m_TPCKuramaTC.at(l); }
  const TPCRKTrackContainer& GetTPCKuramaTracks() const { return m_TPCKuramaTC; }

  //K18Tracking
  void MakeTPCK18HPContainer(TPCLocalTrackHelix *tpctrack,
			     std::vector<Int_t> &lnum);
  TPCRKTrack* MakeTPCK18Track(TPCLocalTrackHelix *tpctrack,
			      const std::vector<Int_t> layerID,
			      const std::vector<Double_t> wire,
			      const std::vector<Double_t> localHitPos);
  Bool_t DoFitTPCK18Track(Int_t K18TrackID,
			  const TVector3 initPos,
			  const TVector3 initMom,
			  const std::vector<Int_t> layerID,
			  const std::vector<Double_t> wire,
			  const std::vector<Double_t> localHitPos,
			  Bool_t U2D);
  Bool_t TrackSearchTPCK18(const std::vector<TVector3> initPos,
			   const std::vector<TVector3> initMom,
			   const std::vector<std::vector<Int_t>> layerID,
			   const std::vector<std::vector<Double_t>> wire,
			   const std::vector<std::vector<Double_t>> localHitPos,
			   Bool_t U2D=true);
  Int_t GetNTracksTPCK18() const { return m_TPCK18TC.size(); }
  TPCRKTrack* GetTPCK18Track(Int_t l) const { return m_TPCK18TC.at(l); }
  const TPCRKTrackContainer& GetTPCK18Tracks() const { return m_TPCK18TC; }

  //Vertex in the target
  Bool_t GetProductionVertex(const TVector3 XtgtKm, const TVector3 PtgtKm,
			     const TVector3 XtgtKp, const TVector3 PtgtKp,
			     TVector3 &Vertex, Double_t &closeDist,
			     Double_t &KmPathInTgt, Double_t &KpPathInTgt,
			     TVector3 &KmMomVtx, TVector3 &KpMomVtx);

  //Vertices between tracks
  Int_t GetNVerticesTPC() const { return m_TPCVC.size(); }
  TPCVertex* GetTPCVertex(Int_t l) const { return m_TPCVC.at(l); }
  const TPCVertexContainer& GetTPCVertices() const { return m_TPCVC; }

  Int_t GetNVerticesTPCClustered() const { return m_TPCVCClustered.size(); }
  TPCVertex* GetTPCVertexClustered(Int_t l) const { return m_TPCVCClustered.at(l); }
  const TPCVertexContainer& GetTPCVerticesClustered() const { return m_TPCVCClustered; }

protected:

  void ClearTPCHits();
  void ClearTPCClusters();
  void ClearTPCTracks();
  void ClearTPCVertices();
  void ClearTPCKuramaTracks();
  void ClearTPCK18Tracks();
  static Bool_t MakeUpTPCClusters(const TPCHitContainer& HitCont,
                                  TPCClusterContainer& ClCont,
                                  Double_t maxdy);
};

//_____________________________________________________________________________
inline const TString&
TPCAnalyzer::ClassName()
{
  static TString s_name("TPCAnalyzer");
  return s_name;
}

#endif

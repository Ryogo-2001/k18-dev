#ifndef CFTParticle_h

#define CFTParticle_h 1

#include <vector>
#include "DCLocalTrack.hh"
#include "ThreeVector.hh"
//#include "RawData.hh"

class CFTLocalTrack;
class DCLocalTrack;
//class DCLTrackHit;
class Hodo1Hit;
class Hodo2Hit;
class RawData;


class CFTParticle
{
public:
  CFTParticle(DCLocalTrack *track, RawData *rawData);
  ~CFTParticle();

private:
  DCLocalTrack* Track_;
  RawData*      RawData_;
  std::vector <Hodo1Hit *> BGOCont_;
  std::vector <Hodo2Hit *> PiVCont_;
  int     m_bgo_seg;
  int     m_piid_seg;
  /*
  double  m_cft_dEsum[NumOfPlaneCFT];//sum of cluster
  double  m_cft_dEmax[NumOfPlaneCFT];//max in cluster
  double  m_cft_TotaldE;
  double  m_cft_TotaldE_phi; // only phi sum
  double  m_cft_TotaldE_uv;  // only uv  sum
  */
  double  PathLength_;
  ThreeVector CFTVtx_;
public:
  void AddBGOHit(Hodo1Hit* hit) {BGOCont_.push_back(hit);}
  void AddPiVHit(Hodo2Hit* hit) {PiVCont_.push_back(hit);}
  void BGOPos(int seg, double *x, double *y) const;
  void PiIDPos(int seg, double *x, double *y) const;
  bool Calculate();

  DCLocalTrack * GetTrack() {return Track_;}
  int NHitBGO() { return BGOCont_.size();}

  int    GetTrackBGOSeg()  { return m_bgo_seg; }
  int    GetTrackPiIDSeg() { return m_piid_seg;}
  /*
  double GetCFTdESum(int layer) { return m_cft_dEsum[layer];}
  double GetCFTdEMax(int layer) { return m_cft_dEmax[layer];}
  double GetCFTdETotal()        { return m_cft_TotaldE;}
  double GetCFTdETotalPhi()     { return m_cft_TotaldE_phi;}
  double GetCFTdETotalUV()      { return m_cft_TotaldE_uv;}
  */
  Hodo1Hit * GetBGOHit(int i); 
  int NHitPiV() { return PiVCont_.size();}
  Hodo2Hit * GetPiVHit(int i); 

  ThreeVector GetPos0 () { return Track_->GetPos0(); }
  ThreeVector GetDir () { return Track_->GetDir(); }
  bool checkProton( double BGO_E, double *delta);
  bool checkPi( double BGO_E);
  void   SetCFTVtx(ThreeVector vtx) {CFTVtx_=vtx;}
  ThreeVector GetCFTVtx() {return CFTVtx_;}
};


#endif

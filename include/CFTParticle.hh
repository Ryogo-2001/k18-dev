#ifndef CFTParticle_h

#define CFTParticle_h 1

#include <vector>
#include "DCLocalTrack.hh"
#include "CFTLocalTrack.hh"
#include "ThreeVector.hh"
//#include "RawData.hh"

class CFTLocalTrack;
class DCLocalTrack;
//class DCLTrackHit;
class HodoWaveformHit;
class HodoHit;
//class RawData;
class HodoAnalyzer;


class CFTParticle
{
public:
  //CFTParticle(DCLocalTrack *track, RawData *rawData);
  CFTParticle(const CFTLocalTrack *track, HodoAnalyzer *hodoAna);
  ~CFTParticle();

private:
  const CFTLocalTrack* m_Track;
  //RawData*      RawData_;
  HodoAnalyzer* m_hodoAna; 
  //std::vector <Hodo1Hit *> m_BGOCont;
  std::vector <const HodoWaveformHit *> m_BGOCont;
  std::vector <HodoHit *> m_PiVCont;
  Int_t     m_bgo_seg;
  Double_t  m_bgo_adc, m_bgo_adc_cor;
  Int_t     m_piid_seg;
  Double_t  m_bgo_energy;
  Double_t  m_cft_Total_dE;
  Double_t  m_cft_Total_dE_max;
  Double_t  m_norm_Total_dE_max;

  Double_t  m_TotalE;
  Double_t  m_TotalE_max;

  Double_t  m_dist_xy; // distance from (0, 0) in xy plane
  Double_t  m_bgo_z_surface;

  Double_t  m_Mass;

  Double_t  m_PathLength;
  ThreeVector m_CFTVtx;
public:
  void AddBGOHit(const HodoWaveformHit* hit) {m_BGOCont.push_back(hit);}
  void AddPiVHit(HodoHit* hit) {m_PiVCont.push_back(hit);}
  void BGOPos(Int_t seg, Double_t *x, Double_t *y) const;
  Double_t BGOZPosAtSurface(Int_t seg) const;
  void PiIDPos(Int_t seg, Double_t *x, Double_t *y) const;
  Bool_t FindBGO();
  Bool_t Calculate();
  Double_t GetBGOEnergy() { return m_bgo_energy;}

  const CFTLocalTrack * GetTrack() {return m_Track;}
  Int_t NHitBGO() { return m_BGOCont.size();}

  Int_t    GetTrackBGOSeg()  { return m_bgo_seg; }
  Int_t    GetTrackPiIDSeg() { return m_piid_seg;}
  /*
  Double_t GetCFTdESum(Int_t layer) { return m_cft_dEsum[layer];}
  Double_t GetCFTdEMax(Int_t layer) { return m_cft_dEmax[layer];}
  Double_t GetCFTdETotal()        { return m_cft_TotaldE;}
  Double_t GetCFTdETotalPhi()     { return m_cft_TotaldE_phi;}
  Double_t GetCFTdETotalUV()      { return m_cft_TotaldE_uv;}
  */
  const HodoWaveformHit * GetBGOHit(Int_t i); 
  Int_t NHitPiV() { return m_PiVCont.size();}
  HodoHit * GetPiVHit(Int_t i); 

  ThreeVector GetPos0 () { return m_Track->GetPos0(); }
  ThreeVector GetDir () { return m_Track->GetDir(); }
  void   SetCFTVtx(ThreeVector vtx) {m_CFTVtx=vtx;}
  ThreeVector GetCFTVtx() {return m_CFTVtx;}
  Double_t GetMass() { return m_Mass;}
  Double_t GetTotalE() { return m_TotalE;}
  Double_t GetTotalEmax() { return m_TotalE_max;}
  Double_t GetDistXY() { return m_dist_xy;}
  Double_t GetBGOZPos() { return m_bgo_z_surface;}
  Double_t GetNormTotalDEmax() { return m_norm_Total_dE_max;}
  Double_t GetBGOAdc() {return m_bgo_adc;}
  Double_t GetBGOAdcCor() {return m_bgo_adc_cor;}
};


#endif

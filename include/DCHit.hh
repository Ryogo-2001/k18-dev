// -*- C++ -*-

#ifndef DC_HIT_HH
#define DC_HIT_HH

#include <cmath>
#include <vector>
#include <deque>
#include <numeric>

#include <TString.h>
#include <TVector3.h>

#include <std_ostream.hh>

class DCRawHit;
class DCLTrackHit;

//_____________________________________________________________________________
class DCHit
{
public:
  static const TString& ClassName();
  DCHit(const DCRawHit* rhit);
  DCHit(Int_t layer);
  DCHit(Int_t layer, Double_t wire);
  DCHit(Int_t plane, Int_t layer, Int_t wire, Double_t wpos); // for Geant4
  ~DCHit();

private:
  DCHit(const DCHit&);
  DCHit& operator =(const DCHit&);

protected:
  const DCRawHit* m_raw_hit;
  Int_t     m_plane; // by DIGIT
  Int_t     m_layer; // by DCGEO
  Double_t  m_wire;  // 0-based (= ch)

  using data_t = std::vector<Double_t>;
  using flag_t = std::deque<Bool_t>;

  data_t m_tdc;
  data_t m_adc;
  data_t m_trailing;

  // DCData normalized
  data_t m_drift_time;
  data_t m_drift_length;
  data_t m_tot;
  flag_t m_belong_to_track;
  flag_t m_is_good;

  Double_t  m_wpos;
  Double_t  m_angle;

  ///// for MWPC
  Int_t    m_cluster_size;
  Bool_t   m_mwpc_flag;
  Double_t m_mwpc_wire;
  Double_t m_mwpc_wpos;

  ///// for TOF
  Double_t m_z;

  ///// for CFT
  // Double_t m_meanseg;
  // Double_t m_maxseg;
  // Double_t m_adc_low;
  // Double_t m_mip_low;
  // Double_t m_dE_low;
  // Double_t max_adc_low;
  // Double_t max_mip_low;
  // Double_t max_dE_low;
  // Double_t m_r;
  // Double_t m_phi;
  // TVector3 m_vtx;
  // Double_t m_pos_phi;
  // Double_t m_pos_z;
  // Double_t m_pos_r;
  // Double_t m_time;

  std::vector<DCLTrackHit*> m_register_container;

public:
  const DCRawHit* GetRawHit() const { return m_raw_hit; }
  Int_t           PlaneId() const { return m_plane; }
  Int_t           LayerId() const { return m_layer; }
  Int_t           GetLayer() const { return LayerId(); }
  Double_t        WireId() const { return m_wire; }
  Double_t        GetWire() const {
    if(m_mwpc_flag) return m_mwpc_wire;
    else return m_wire;
  }
  Bool_t CalcDCObservables();
  Bool_t SetDCObservablesGeant4(Double_t dl, Double_t tot);
  Bool_t CalcFiberObservables();
  Bool_t CalcMWPCObservables();
  Bool_t CalcCFTObservables();
  Bool_t CalcObservablesSimulation(Double_t dlength);
  Bool_t Calculate(){ return CalcDCObservables(); }

  void ClearDCData();
  void EraseDCData(Int_t i);

  Int_t GetTdcSize() const { return m_tdc.size(); }
  Int_t GetAdcSize() const { return m_adc.size(); }
  Int_t GetTdcVal(Int_t i=0) const { return m_tdc[i]; }
  Int_t GetAdcVal(Int_t i=0) const { return m_adc[i]; }
  Int_t GetTdcTrailing(Int_t i=0) const { return m_trailing[i]; }
  Int_t GetTdcTrailingSize() const { return m_trailing.size(); }
  Int_t GetTdc1st() const;
  Double_t GetTiltAngle() const { return m_angle; }
  Double_t GetWirePosition() const {
    if(m_mwpc_flag) return m_mwpc_wpos;
    else return m_wpos;
  }
  Int_t GetClusterSize() const { return m_cluster_size; }
  Double_t GetMeamWire() const { return m_mwpc_wire; }
  Double_t GetMeamWirePosition() const { return m_mwpc_wpos; }
  Double_t GetResolution() const;
  Double_t GetZ() const { return m_z; } // for TOF

  Int_t GetEntries() const { return m_drift_time.size(); }
  Double_t GetDriftTime(Int_t i) const { return m_drift_time.at(i); }
  Double_t GetDriftLength(Int_t i) const { return m_drift_length.at(i); }
  Double_t TimeOverThreshold(Int_t i=0) const { return m_tot.at(i); }
  void JoinTrack(Int_t i=0) { m_belong_to_track.at(i) = true; }
  void QuitTrack(Int_t i=0) { m_belong_to_track.at(i) = false; }
  Bool_t BelongToTrack(Int_t i=0) const { return m_belong_to_track.at(i); }
  Bool_t IsGood(Int_t i=0) const { return m_is_good.at(i); }
  Bool_t IsEmpty() const { return GetEntries() == 0; }

  void SetLayer(Int_t layer){ m_layer = layer; }
  void SetWire(Double_t wire){ m_wire  = wire; }
  void SetTdcVal(Int_t tdc){ m_tdc.push_back(tdc); }
  void SetTdc(Int_t tdc){ m_tdc.push_back(tdc); }
  void SetAdcVal(Int_t adc){ m_adc.push_back(adc); }
  void SetTrailing(Int_t tdc){ m_trailing.push_back(tdc); }
  void SetTdcTrailing(Int_t tdc){ m_trailing.push_back(tdc); }
  void SetDCData(Double_t dt=0., Double_t dl=0.,
                 Double_t tot=TMath::QuietNaN(),
                 Bool_t belong_to_track=false,
                 Bool_t is_good=true);
  void SetDriftLength(Double_t dl){ m_drift_length.push_back(dl); }
  void SetDriftTime(Double_t dt){ m_drift_time.push_back(dt); }
  void SetTiltAngle(Double_t angleDegree){ m_angle = angleDegree; }
  void SetClusterSize(Int_t size){ m_cluster_size = size; }
  void SetMWPCFlag(Bool_t flag){ m_mwpc_flag = flag; }
  void SetMeanWire(Double_t mwire){ m_mwpc_wire = mwire; }
  void SetMeanWirePosition(Double_t mwpos){ m_mwpc_wpos = mwpos; }
  void SetWirePosition(Double_t wpos){ m_wpos = wpos; }
  void SetZ(Double_t z) { m_z = z; } // for TOF

  // aliases
  Int_t GetTdc(Int_t i=0) const { return GetTdcVal(i); }
  Int_t Tdc(Int_t i=0) const { return GetTdcVal(i); }
  Int_t Tdc1st() const { return GetTdc1st(); }
  Int_t GetTrailing(Int_t i=0) const { return GetTdcTrailing(i); }
  Int_t Trailing(Int_t i=0) const { return GetTdcTrailing(i); }
  Double_t DT(Int_t i=0) const { return DriftTime(i); }
  Double_t DL(Int_t i=0) const { return DriftLength(i); }
  Double_t DriftTime(Int_t i=0) const { return GetDriftTime(i); }
  Double_t DriftLength(Int_t i=0) const { return GetDriftLength(i); }
  Double_t TOT(Int_t i=0) const { return TimeOverThreshold(i); }
  Double_t GetTot(Int_t i=0) const { return TimeOverThreshold(i); }
  Double_t TrailingTime(Int_t i=0) const { return DT(i) + TOT(i); }
  Double_t GetDriftTimeSize() const { return GetEntries(); }
  Double_t GetDriftLengthSize() const { return GetEntries(); }

  void DriftTimeCut(Double_t min, Double_t max, Bool_t select_1st);
  void TotCut(Double_t min, Bool_t keep_nan);
  void Print(Option_t* arg="") const;

  ///// for CFT
  // void SetMeanSeg(Double_t seg) { m_meanseg = seg; }
  // void SetMaxSeg(Double_t seg) { m_maxseg = seg; }
  // void SetAdcLow(Double_t adc) { m_adc_low = adc; }
  // void SetMIPLow(Double_t mip) { m_mip_low = mip; }
  // void SetdELow(Double_t dE) { m_dE_low = dE; }
  // void SetMaxAdcLow(Double_t adc) { max_adc_low = adc; }
  // void SetMaxMIPLow(Double_t mip) { max_mip_low = mip; }
  // void SetMaxdELow(Double_t dE) { max_dE_low = dE; }
  // void SetPositionR(Double_t r) { m_r = r; }
  // void SetPositionPhi(Double_t phi) { m_phi = phi; }
  // void SetPosPhi(Double_t phi) { m_pos_phi  = phi; }
  // void SetPosZ(Double_t z) { m_pos_z = z; }
  // void SetPosR(Double_t r) { m_pos_r = r; }
  // void SetTime(Double_t time) { m_time = time; }
  // void SetVtx(const TVector3& vtx) { m_vtx = vtx; }

  ///// for CFT
  // Double_t GetMeanSeg() const { return m_meanseg; }
  // Double_t GetMaxSeg() const { return m_maxseg;  }
  // Double_t SetMaxSeg() const { return m_maxseg;  }
  // Double_t GetAdcLow() const { return m_adc_low; }
  // Double_t GetMIPLow() const { return m_mip_low; }
  // Double_t GetdELow() const { return m_dE_low ; }
  // Double_t GetMaxAdcLow() const { return max_adc_low; }
  // Double_t GetMaxMIPLow() const { return max_mip_low; }
  // Double_t GetMaxdELow() const { return max_dE_low ; }

  // Double_t GetPositionR() const { return m_r;        }
  // Double_t GetPositionPhi() const { return m_phi;      }
  // Double_t GetPosPhi() const { return m_pos_phi;  }
  // Double_t GetPosZ() const { return m_pos_z;    }
  // Double_t GetPosR() const { return m_pos_r;    }
  // Double_t GetTime() const { return m_time;     }
  // TVector3 GetVtx() const { return m_vtx;      }

  void RegisterHits(DCLTrackHit *hit)
    { m_register_container.push_back(hit); }

  Bool_t ReCalcDC(Bool_t applyRecursively=false) { return CalcDCObservables(); }
  // Bool_t ReCalcMWPC(Bool_t applyRecursively=false) { return CalcMWPCObservables(); }
  static Bool_t Compare(const DCHit* left, const DCHit* right);

protected:
  void ClearRegisteredHits();
};

//_____________________________________________________________________
inline const TString&
DCHit::ClassName()
{
  static TString s_name("DCHit");
  return s_name;
}

//_____________________________________________________________________________
inline Bool_t
DCHit::Compare(const DCHit* left, const DCHit* right)
{
  if(left->LayerId() == right->LayerId())
    return left->WireId() < right->WireId();
  else
    return left->LayerId() < right->LayerId();
}

#endif

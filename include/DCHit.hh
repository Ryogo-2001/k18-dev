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

#include "DebugCounter.hh"

//typedef std::vector<Bool_t>   BoolVec;
typedef std::deque<Bool_t>    BoolVec;
typedef std::vector<Int_t>    IntVec;
typedef std::vector<Double_t> DoubleVec;

class DCLTrackHit;

//_____________________________________________________________________________
class DCHit
{
public:
  static const TString& ClassName();
  DCHit();
  DCHit(Int_t layer);
  DCHit(Int_t layer, Double_t wire);
  ~DCHit();

private:
  DCHit(const DCHit&);
  DCHit& operator =(const DCHit&);

protected:
  Int_t     m_layer;
  Double_t  m_wire;
  IntVec    m_tdc;
  IntVec    m_adc;
  IntVec    m_trailing;

  //  DoubleVec m_dt;
  //  DoubleVec m_dl;
  //  DoubleVec m_trailing_time;

  // For DC with HUL MH-TDC
  struct data_pair
  {
    Double_t drift_time;
    Double_t drift_length;
    Double_t trailing_time;
    Double_t tot;
    Int_t    index_t;
    Bool_t   belong_track;
    Bool_t   dl_range;
  };

  std::vector<data_pair> m_pair_cont;

  Double_t  m_wpos;
  Double_t  m_angle;
  //  BoolVec m_belong_track;
  //  BoolVec m_dl_range;

  ///// for MWPC
  Int_t    m_cluster_size;
  Bool_t   m_mwpc_flag;
  Double_t m_mwpc_wire;
  Double_t m_mwpc_wpos;

  ///// for TOF
  Double_t m_z;

  ///// For E40 Acrylic TOF
  Double_t m_ofs_dt;

  ///// for CFT
  Double_t m_meanseg;
  Double_t m_maxseg;
  Double_t m_adc_low;
  Double_t m_mip_low;
  Double_t m_dE_low;
  Double_t max_adc_low;
  Double_t max_mip_low;
  Double_t max_dE_low;
  Double_t m_r;
  Double_t m_phi;
  BoolVec  m_belong_track;
  TVector3 m_vtx;
  Double_t m_pos_phi;
  Double_t m_pos_z;
  Double_t m_pos_r;
  Double_t m_time;

  ///// for TPC
  Int_t    m_hitnum;

  std::vector<DCLTrackHit*> m_register_container;

public:
  Bool_t CalcDCObservables();
  Bool_t CalcMWPCObservables();
  Bool_t CalcFiberObservables();
  Bool_t CalcCFTObservables();
  //  Bool_t CalcObservablesSimulation(Double_t dlength);

  void SetLayer(Int_t layer)              { m_layer = layer;                    }
  void SetWire(Double_t wire)             { m_wire  = wire;                     }
  void SetTdcVal(Int_t tdc)               { m_tdc.push_back(tdc);               }
  void SetAdcVal(Int_t adc)               { m_adc.push_back(adc);               }
  void SetTdcTrailing(Int_t tdc)            { m_trailing.push_back(tdc);          }
  void SetDummyPair();
  void SetDriftLength(Int_t ith, Double_t dl){m_pair_cont.at(ith).drift_length = dl;}
  void SetDriftTime(Int_t ith, Double_t dt) { m_pair_cont.at(ith).drift_time = dt;  }
  void SetTiltAngle(Double_t angleDegree) { m_angle = angleDegree;              }

  void SetClusterSize(Int_t size)          { m_cluster_size   = size;           }
  void SetMWPCFlag(Bool_t flag)            { m_mwpc_flag = flag;                }
  void SetMeanWire(Double_t mwire)         { m_mwpc_wire = mwire;               }
  void SetMeanWirePosition(Double_t mwpos) { m_mwpc_wpos = mwpos;               }
  void SetWirePosition(Double_t wpos)      { m_wpos      = wpos;                }

  ///// for TOF
  void SetZ(Double_t z) { m_z = z; }

  ///// for TPC
  void SetHitNum(Int_t hitnum) { m_hitnum = hitnum; }

  ///// For E40 Acrylic TOF
  void SetOfsdT(Double_t ofs) { m_ofs_dt = ofs; }

  ///// for CFT
  void SetMeanSeg(Double_t seg) { m_meanseg = seg; }
  void SetMaxSeg(Double_t seg) { m_maxseg = seg; }
  void SetAdcLow(Double_t adc) { m_adc_low = adc; }
  void SetMIPLow(Double_t mip) { m_mip_low = mip; }
  void SetdELow(Double_t dE) { m_dE_low = dE; }
  void SetMaxAdcLow(Double_t adc) { max_adc_low = adc; }
  void SetMaxMIPLow(Double_t mip) { max_mip_low = mip; }
  void SetMaxdELow(Double_t dE) { max_dE_low = dE; }
  void SetPositionR(Double_t r) { m_r = r; }
  void SetPositionPhi(Double_t phi) { m_phi = phi; }
  void SetPosPhi(Double_t phi) { m_pos_phi  = phi; }
  void SetPosZ(Double_t z) { m_pos_z = z; }
  void SetPosR(Double_t r) { m_pos_r = r; }
  void SetTdcCFT(Int_t tdc); // not used
  void SetTime(Double_t time) { m_time = time; }
  void SetVtx(const TVector3& vtx) { m_vtx = vtx; }


  void GateDriftTime(Double_t min, Double_t max, Bool_t select_1st);

  Int_t GetLayer() const { return m_layer; }
  Double_t GetWire()  const {
    if(m_mwpc_flag) return m_mwpc_wire;
    else return Int_t(m_wire);
  }

  Int_t GetTdcSize()             const { return m_tdc.size(); }
  Int_t GetAdcSize()             const { return m_adc.size(); }
  Int_t GetDriftTimeSize()       const { return m_pair_cont.size(); }
  Int_t GetDriftLengthSize()     const { return m_pair_cont.size(); }
  Int_t GetTdcVal(Int_t nh=0)          const { return m_tdc[nh]; }
  Int_t GetAdcVal(Int_t nh=0)          const { return m_adc[nh]; }
  Int_t GetTdcTrailing(Int_t nh=0)     const { return m_trailing[nh]; }
  Int_t GetTdcTrailingSize()     const { return m_trailing.size(); }

  Double_t GetResolution()       const;

  Double_t GetDriftTime(Int_t nh=0)    const { return m_pair_cont.at(nh).drift_time; }
  Double_t GetDriftLength(Int_t nh=0)  const { return m_pair_cont.at(nh).drift_length; }
  Double_t GetTrailingTime(Int_t nh=0) const { return m_pair_cont.at(nh).trailing_time; }
  Double_t GetTot(Int_t nh=0)          const { return m_pair_cont.at(nh).tot; }

  Double_t GetTiltAngle() const { return m_angle; }
  Double_t GetWirePosition() const {
    if(m_mwpc_flag) return m_mwpc_wpos;
    else return m_wpos;
  }

  Int_t GetClusterSize() const { return m_cluster_size; }
  Double_t GetMeamWire() const { return m_mwpc_wire; }
  Double_t GetMeamWirePosition() const { return m_mwpc_wpos; }

  ///// for TOF
  Double_t GetZ() const { return m_z; }

  ///// for TPC
  Double_t GetHitNum() { return m_hitnum; }

  ///// for CFT
  Double_t GetMeanSeg() const { return m_meanseg; }
  Double_t GetMaxSeg() const { return m_maxseg;  }
  Double_t SetMaxSeg() const { return m_maxseg;  }
  Double_t GetAdcLow() const { return m_adc_low; }
  Double_t GetMIPLow() const { return m_mip_low; }
  Double_t GetdELow() const { return m_dE_low ; }
  Double_t GetMaxAdcLow() const { return max_adc_low; }
  Double_t GetMaxMIPLow() const { return max_mip_low; }
  Double_t GetMaxdELow() const { return max_dE_low ; }

  Double_t GetPositionR() const { return m_r;        }
  Double_t GetPositionPhi() const { return m_phi;      }
  Double_t GetPosPhi() const { return m_pos_phi;  }
  Double_t GetPosZ() const { return m_pos_z;    }
  Double_t GetPosR() const { return m_pos_r;    }
  Double_t GetTime() const { return m_time;     }
  TVector3 GetVtx() const { return m_vtx;      }

  void JoinTrack(Int_t nh=0) { m_pair_cont.at(nh).belong_track = true; }
  void QuitTrack(Int_t nh=0) { m_pair_cont.at(nh).belong_track = false;}
  Bool_t BelongToTrack(Int_t nh=0) const { return m_pair_cont.at(nh).belong_track; }
  Bool_t IsWithinRange(Int_t nh=0) const { return m_pair_cont.at(nh).dl_range; }

  void JoinTrackCFT(Int_t nh=0) {m_belong_track[nh] = true; }
  void QuitTrackCFT(Int_t nh=0) {m_belong_track[nh] = false; }
  Bool_t BelongToTrackCFT(Int_t nh=0) const { return m_belong_track[nh]; }

  void RegisterHits(DCLTrackHit *hit)
    { m_register_container.push_back(hit); }

  Bool_t ReCalcDC(Bool_t applyRecursively=false) { return CalcDCObservables(); }
  Bool_t ReCalcMWPC(Bool_t applyRecursively=false) { return CalcMWPCObservables(); }

  void TotCut(Double_t min_tot, Bool_t adopt_nan);

  void Print(const TString& arg="", std::ostream& ost=hddaq::cout) const;

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

//_____________________________________________________________________
inline std::ostream&
operator <<(std::ostream& ost, const DCHit& hit)
{
  hit.Print("", ost);
  return ost;
}

#endif

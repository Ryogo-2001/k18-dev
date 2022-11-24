// -*- C++ -*-

#ifndef RAW_DATA_HH
#define RAW_DATA_HH

#include <deque>
#include <vector>

#include <TString.h>

#include "DetectorID.hh"

class HodoRawHit;
class DCRawHit;
class TPCRawHit;

typedef std::vector<HodoRawHit*> HodoRHitContainer;
typedef std::vector<DCRawHit*>   DCRHitContainer;
typedef std::vector<TPCRawHit*>  TPCRHitContainer;
typedef std::vector<Int_t>       FADCRHitContainer;

//_____________________________________________________________________________
class RawData
{
public:
  static TString ClassName();
  RawData();
  ~RawData();

private:
  RawData(const RawData&);
  RawData& operator=(const RawData&);

private:
  enum EType { kTPC, kOthers, kNType };
  std::deque<Bool_t>             m_is_decoded;
  HodoRHitContainer              m_CaenV792RawHC;
  HodoRHitContainer              m_BH1RawHC;
  HodoRHitContainer              m_BH2RawHC;
  HodoRHitContainer              m_BACRawHC;
  HodoRHitContainer              m_HTOFRawHC;
  HodoRHitContainer              m_SCHRawHC;
  HodoRHitContainer              m_BVHRawHC;
  HodoRHitContainer              m_TOFRawHC;
  HodoRHitContainer              m_LACRawHC;
  HodoRHitContainer              m_WCRawHC;
  HodoRHitContainer              m_WCSUMRawHC;
  HodoRHitContainer              m_E72BACRawHC;
  HodoRHitContainer              m_E72BACSUMRawHC;
  HodoRHitContainer              m_KVCRawHC;
  HodoRHitContainer              m_KVCSUMRawHC;
  HodoRHitContainer              m_E90SACRawHC;
  HodoRHitContainer              m_E90SACSUMRawHC;
  HodoRHitContainer              m_T4RawHC;
  HodoRHitContainer              m_T5RawHC;
  HodoRHitContainer              m_T6RawHC;
  HodoRHitContainer              m_T7RawHC;
  std::vector<HodoRHitContainer> m_BFTRawHC;
  std::vector<DCRHitContainer>   m_BcInRawHC;
  std::vector<DCRHitContainer>   m_BcOutRawHC;
  std::vector<TPCRHitContainer>  m_TPCRawHC;
  std::vector<TPCRHitContainer>  m_TPCCorHC;
  std::vector<DCRHitContainer>   m_SdcInRawHC;
  std::vector<DCRHitContainer>   m_SdcOutRawHC;
  HodoRHitContainer              m_TPCClockRawHC;
  HodoRHitContainer              m_ScalerRawHC;
  HodoRHitContainer              m_TrigRawHC;
  HodoRHitContainer              m_VmeCalibRawHC;
  TPCRawHit*                     m_baseline;

public:
  void                     ClearAll();
  void                     ClearTPC();
  Bool_t                   CorrectBaselineTPC();
  Bool_t                   DecodeHits();
  Bool_t                   DecodeCalibHits();
  Bool_t                   DecodeTPCHits();
  Bool_t                   SelectTPCHits(Bool_t maxadccut,
                                         Bool_t maxadctbcut);
  const TPCRawHit* const   GetBaselineTPC() const { return m_baseline; }
  const HodoRHitContainer& GetCaenV792RawHC() const;
  const HodoRHitContainer& GetBH1RawHC() const;
  const HodoRHitContainer& GetBH2RawHC() const;
  const HodoRHitContainer& GetBACRawHC() const;
  const HodoRHitContainer& GetHTOFRawHC() const;
  const HodoRHitContainer& GetSCHRawHC() const;
  const HodoRHitContainer& GetBVHRawHC() const;
  const HodoRHitContainer& GetTOFRawHC() const;
  const HodoRHitContainer& GetLACRawHC() const;
  const HodoRHitContainer& GetWCRawHC() const;
  const HodoRHitContainer& GetWCSUMRawHC() const;
  const HodoRHitContainer& GetBFTRawHC(Int_t plane) const;
  const DCRHitContainer&   GetBcInRawHC(Int_t layer) const;
  const DCRHitContainer&   GetBcOutRawHC(Int_t layer) const;
  const DCRHitContainer&   GetSdcInRawHC(Int_t layer) const;
  const DCRHitContainer&   GetSdcOutRawHC(Int_t layer) const;
  const TPCRHitContainer&  GetTPCRawHC(Int_t layer) const;
  const TPCRHitContainer&  GetTPCCorHC(Int_t layer) const;
  const HodoRHitContainer& GetTPCClockRawHC() const;
  const HodoRHitContainer& GetScalerRawHC() const;
  const HodoRHitContainer& GetTrigRawHC() const;
  const HodoRHitContainer& GetVmeCalibRawHC() const;
  const HodoRHitContainer& GetE72BACRawHC() const;
  const HodoRHitContainer& GetE72BACSUMRawHC() const;
  const HodoRHitContainer& GetKVCRawHC() const;
  const HodoRHitContainer& GetKVCSUMRawHC() const;
  const HodoRHitContainer& GetE90SACRawHC() const;
  const HodoRHitContainer& GetE90SACSUMRawHC() const;
  const HodoRHitContainer& GetT4RawHC() const;
  const HodoRHitContainer& GetT5RawHC() const;
  const HodoRHitContainer& GetT6RawHC() const;
  const HodoRHitContainer& GetT7RawHC() const;

private:
  enum EDCDataType { kDcLeading, kDcTrailing, kDcOverflow, kDcNDataType };
  Bool_t AddHodoRawHit(HodoRHitContainer& cont,
                       Int_t id, Int_t plane, Int_t seg,
                       Int_t UorD, Int_t type, Int_t data);
  Bool_t AddDCRawHit(DCRHitContainer& cont,
                     Int_t plane, Int_t wire, Int_t data,
                     Int_t type=kDcLeading);
  Bool_t AddTPCRawHit(TPCRHitContainer& cont,
                      Int_t layer, Int_t row, Double_t adc,
                      Double_t* pars=nullptr);
  void   DecodeHodo(Int_t id, Int_t plane, Int_t nseg, Int_t nch,
                    HodoRHitContainer& cont);
  void   DecodeHodo(Int_t id, Int_t nseg, Int_t nch,
                    HodoRHitContainer& cont);
};

//_____________________________________________________________________________
inline TString
RawData::ClassName()
{
  static TString s_name("RawData");
  return s_name;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetCaenV792RawHC() const
{
  return m_CaenV792RawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetBH1RawHC() const
{
  return m_BH1RawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetBH2RawHC() const
{
  return m_BH2RawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetBACRawHC() const
{
  return m_BACRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetHTOFRawHC() const
{
  return m_HTOFRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetE72BACRawHC() const
{
  return m_E72BACRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetE72BACSUMRawHC() const
{
  return m_E72BACSUMRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetKVCRawHC() const
{
  return m_KVCRawHC;
}

inline const HodoRHitContainer&
RawData::GetKVCSUMRawHC() const
{
  return m_KVCSUMRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetE90SACRawHC() const
{
  return m_E90SACRawHC;
}

inline const HodoRHitContainer&
RawData::GetE90SACSUMRawHC() const
{
  return m_E90SACSUMRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetT4RawHC() const
{
  return m_T4RawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetT5RawHC() const
{
  return m_T5RawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetT6RawHC() const
{
  return m_T6RawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetT7RawHC() const
{
  return m_T7RawHC;
}


//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetSCHRawHC() const
{
  return m_SCHRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetBVHRawHC() const
{
  return m_BVHRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetTOFRawHC() const
{
  return m_TOFRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetLACRawHC() const
{
  return m_LACRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetWCRawHC() const
{
  return m_WCRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetWCSUMRawHC() const
{
  return m_WCSUMRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetBFTRawHC(Int_t plane) const
{
  if(plane<0 || plane>NumOfPlaneBFT-1) plane=0;
  return m_BFTRawHC[plane];
}

//_____________________________________________________________________________
inline const DCRHitContainer&
RawData::GetBcInRawHC(Int_t layer) const
{
  if(layer<0 || layer>NumOfLayersBcIn) layer = 0;
  return m_BcInRawHC[layer];
}

//_____________________________________________________________________________
inline const DCRHitContainer&
RawData::GetBcOutRawHC(Int_t layer) const
{
  if(layer<0 || layer>NumOfLayersBcOut) layer = 0;
  return m_BcOutRawHC[layer];
}

//_____________________________________________________________________________
inline const TPCRHitContainer&
RawData::GetTPCRawHC(Int_t layer) const
{
  if(layer<0 || layer>NumOfLayersTPC) layer = 0;
  return m_TPCRawHC[layer];
}

//_____________________________________________________________________________
inline const TPCRHitContainer&
RawData::GetTPCCorHC(Int_t layer) const
{
  if(layer<0 || layer>NumOfLayersTPC) layer = 0;
  return m_TPCCorHC[layer];
}

//_____________________________________________________________________________
inline const DCRHitContainer&
RawData::GetSdcInRawHC(Int_t layer) const
{
  if(layer<0 || layer>NumOfLayersSdcIn) layer = 0;
  return m_SdcInRawHC[layer];
}

//_____________________________________________________________________________
inline const DCRHitContainer&
RawData::GetSdcOutRawHC(Int_t layer) const
{
  if(layer<0 || layer>NumOfLayersSdcOut) layer = 0;
  return m_SdcOutRawHC[layer];
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetTPCClockRawHC() const
{
  return m_TPCClockRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetScalerRawHC() const
{
  return m_ScalerRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetTrigRawHC() const
{
  return m_TrigRawHC;
}

//_____________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetVmeCalibRawHC() const
{
  return m_VmeCalibRawHC;
}

#endif

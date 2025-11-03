// -*- C++ -*-

#ifndef CFT_FIBER_CLUSTER_HH
#define CFT_FIBER_CLUSTER_HH

#include <vector>

#include <TString.h>

#include "CFTFiberHit.hh"
#include "HodoCluster.hh"

//_____________________________________________________________________________
class CFTFiberCluster : public HodoCluster
{
public:
  static const TString& ClassName();
  CFTFiberCluster(const HodoHC& cont,
               const index_t& index);
  virtual ~CFTFiberCluster();

private:
  CFTFiberCluster(const CFTFiberCluster&);
  CFTFiberCluster& operator =(const CFTFiberCluster&);

protected:
  Double_t m_total_adc_hi;
  Double_t m_max_adc_hi;  
  Double_t m_total_adc_low;
  Double_t m_max_adc_low;  
  Double_t m_total_mip_hi;
  Double_t m_max_mip_hi;  
  Double_t m_total_mip_low;
  Double_t m_max_mip_low;  
  Double_t m_total_de_hi;
  Double_t m_max_de_hi;  
  Double_t m_total_de_low;
  Double_t m_max_de_low;  

  Double_t m_max_ctime;
  Int_t    m_max_segment;

  Double_t m_mean_x;
  Double_t m_mean_y;  
  Double_t m_mean_phi;
  Double_t m_mean_r;  
  Double_t m_mean_z0;
  Double_t m_slope;

  Double_t m_z;

  Double_t m_mean_x_cor;
  Double_t m_mean_y_cor;
  Double_t m_mean_phi_cor;

  Double_t m_x_cal;
  Double_t m_y_cal;
  Double_t m_z_cal;
  Double_t m_phi_cal;

  Double_t m_pathlength;

  Bool_t   m_belongTrack;
  Int_t    m_xyFitFlag;
  
public:
  Double_t MaxAdcLow() const { return m_max_adc_low; }
  Double_t MaxMipLow() const { return m_max_mip_low; }  
  Int_t    MaxSegment() const { return m_max_segment; }  
  Double_t MaxDeltaE() const { return m_max_de_low; }
  Double_t TotalDeltaE() const { return m_total_de_low; }  
  Double_t MaxDeltaE_Hi() const { return m_max_de_hi; }
  Double_t TotalDeltaE_Hi() const { return m_total_de_hi; }  
  Double_t MeanX() const  { return m_mean_x; }
  Double_t MeanY() const  { return m_mean_y; }
  Double_t MeanPhi() const  { return m_mean_phi; }
  Double_t MeanR() const  { return m_mean_r; }
  Double_t MeanZ0() const  { return m_mean_z0; }
  Double_t Slope() const  { return m_slope; }
  Double_t CMeanTimeMaxCluster() const { return m_max_ctime; }

  void     setFlags() { m_belongTrack = true; }
  void     clearFlags() { m_belongTrack = false; }
  Bool_t   showFlags() const { return m_belongTrack; }
  
  Double_t MeanXCor() const  { return m_mean_x_cor; }
  Double_t MeanYCor() const  { return m_mean_y_cor; }
  Double_t MeanPhiCor() const  { return m_mean_phi_cor; }

  void     SetZ(Double_t z) { m_z = z; }
  void     SetCalPosition(Double_t x, Double_t y) { m_x_cal = x; m_y_cal = y; }
  void     SetCalZ(Double_t z) { m_z_cal = z; }
  void     SetCalPhi(Double_t phi) { m_phi_cal = phi; }
  void     SetPathLength(Double_t len) {m_pathlength = len; }

  void      SetXYFitFlag(Int_t flag) { m_xyFitFlag = flag; }
  Double_t    GetXcal (void) const { return m_x_cal; }
  Double_t    GetYcal (void) const { return m_y_cal; }
  Double_t    GetZcal (void) const { return m_z_cal; }
  Double_t    GetCalPhi (void) const { return m_phi_cal; }
  Double_t    GetPathLength (void) const { return m_pathlength; }
  
  Double_t    NormalizedTotalDEHiGain (void) const { return m_total_de_hi/m_pathlength; }
  Double_t    NormalizedTotalDELowGain (void) const { return m_total_de_low/m_pathlength; }
  Double_t    NormalizedMaxDEHiGain (void) const { return m_max_de_hi/m_pathlength; }
  Double_t    NormalizedMaxDELowGain (void) const { return m_max_de_low/m_pathlength; }
  
  inline Double_t  GetResidual (void);
  Double_t  GetResidualPhi (void) const { return m_mean_phi - m_phi_cal; }
  Double_t  GetResidualPhiCor (void) const { return m_mean_phi_cor - m_phi_cal; }
  void      SetCorPhi();
  Double_t  GetResidualZ (void) const { return m_z - m_z_cal; }
  Double_t  GetZIni(void)  const { return m_z; };
  Double_t  GetCorZ(Double_t phi, Double_t mean_z, Double_t theta);
  
protected:
  void Calculate();
};

//_____________________________________________________________________________
inline const TString&
CFTFiberCluster::ClassName()
{
  static TString s_name("CFTFiberCluster");
  return s_name;
}

//_____________________________________________________________________________
inline Double_t CFTFiberCluster::GetResidual()
{
  if (m_xyFitFlag == 0) {
    return m_mean_y - m_y_cal;
  } else if (m_xyFitFlag == 1) {
    return m_mean_x - m_x_cal;
  } else {
    //std::cout << "CFTFiberCluster::GetResidual : invalid xyFitFlag=" << xyFitFlag_ << std::endl;
    return -999.;
  }
}

#endif

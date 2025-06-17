// -*- C++ -*-

#ifndef TPC_HIT_HH
#define TPC_HIT_HH

#include "DCHit.hh"

#include <cmath>
#include <vector>
#include <deque>
#include <string>
#include <numeric>

#include <std_ostream.hh>

#include <TVector3.h>

class TPCCluster;
class TPCLTrackHit;
class TPCRawHit;

//_____________________________________________________________________________
class TPCHit : public DCHit
{
public:
  static TString ClassName();
  TPCHit(TPCRawHit* rhit);
  TPCHit(Int_t layer, Double_t mrow); // for cluster hit
  ~TPCHit();

private:
  TPCHit();
  TPCHit(const TPCHit&);
  TPCHit& operator =(const TPCHit&);

protected:
  TPCRawHit*            m_rhit;
  Int_t                 m_layer;
  Int_t                 m_row;
  Double_t              m_padtheta;
  Double_t              m_padlength;
  Double_t              m_mrow;
  Int_t                 m_pad;
  Double_t              m_pedestal;
  Double_t              m_rms;
  Double_t 		m_raw_rms;
  std::vector<Double_t> m_de;
  std::vector<Double_t> m_sigma;
  std::vector<Double_t> m_time;
  std::vector<Double_t> m_chisqr;
  std::vector<Double_t> m_cde;
  std::vector<Double_t> m_ctime;
  std::vector<Double_t> m_drift_length; // this means Y (beam height = 0)
  std::vector<TVector3> m_position;
  Bool_t                m_is_good;
  Int_t                 m_is_calculated;
  Int_t                 m_hough_flag;
  std::vector<Int_t>    m_houghY_num;
  Double_t              m_hough_dist;
  Double_t              m_hough_disty;
  std::vector<Double_t> m_res_param;
  TPCCluster*           m_parent_cluster = nullptr;

  std::vector<TPCLTrackHit*> m_register_container;

public:
  void            AddHit(Double_t de, Double_t time, Double_t sigma=0.,
                         Double_t chisqr=0.);
  Bool_t          Calculate(Double_t clock=0.);
  Bool_t          DoFit();
  Double_t        GetCDe(Int_t i=0) const { return m_cde.at(i); }
  Int_t           GetCDeSize() const { return m_cde.size(); }
  Double_t        GetChisqr(Int_t i=0) const { return m_chisqr.at(i); }
  Int_t           GetChisqrSize() const { return m_chisqr.size(); }
  Double_t        GetCTime(Int_t i=0) const { return m_ctime.at(i); }
  Int_t           GetCTimeSize() const { return m_ctime.size(); }
  Double_t        GetDe(Int_t i=0) const { return m_de.at(i); }
  Double_t        GetSigma(Int_t i=0) const { return m_sigma.at(i); }
  Int_t           GetDeSize() const { return m_de.size(); }
  Double_t        GetDriftLength(Int_t i=0) const
  { return m_drift_length.at(i); }
  Int_t           GetDriftLengthSize() const
  { return m_drift_length.size(); }
  Int_t           GetNHits() const { return m_de.size(); }
  TPCCluster*     GetParentCluster() const { return m_parent_cluster; }
  Int_t           GetPad() const { return m_pad; }
  TPCRawHit*      GetRawHit() const { return m_rhit; }
  Int_t           GetLayer() const { return m_layer; }
  Int_t           GetRow() const { return m_row; }
  Double_t        GetPedestal() const { return m_pedestal; }
  Double_t        GetRMS() const { return m_rms; }
  Double_t	  GetRawRMS()const{return m_raw_rms;}
  Double_t        GetX(Int_t i=0) const { return m_position.at(i).X(); }
  Double_t        GetY(Int_t i=0) const { return m_position.at(i).Y(); }
  Double_t        GetZ(Int_t i=0) const { return m_position.at(i).Z(); }
  Double_t        GetPadTheta() { return m_padtheta; }
  Double_t        GetPadLength() const { return m_padlength; }
  Double_t        GetMRow() const { return m_mrow; }

  Int_t           GetHoughY_num(Int_t i) const { return m_houghY_num.at(i); }
  Int_t           GetHoughY_num_size() const { return m_houghY_num.size(); }
  const TVector3& GetPosition(Int_t i=0) const { return m_position.at(i); }
  const std::vector<Double_t>& GetResolutionParams() const { return m_res_param; }
  Double_t        GetTime(Int_t i=0) const { return m_time.at(i); }
  Int_t           GetTimeSize() const { return m_time.size(); }
  Bool_t          IsGood() const;
  Int_t           GetHoughFlag() const { return m_hough_flag; }
  Double_t        GetHoughDist() const { return m_hough_dist; }
  Double_t        GetHoughDistY() const { return m_hough_disty; }
  void            Print(const std::string& arg="", std::ostream& ost=hddaq::cout) const;
  void            RegisterHits(TPCLTrackHit *hit){ m_register_container.push_back(hit); }
  void            SetIsGood(Bool_t flag=true) { m_is_good=flag; }
  void            SetDe(Double_t de) { m_de.at(0) = de; m_cde.at(0) = de; }
  void            SetPad(Int_t pad) { m_pad = pad; }
  void            SetLayer(Int_t layer) { m_layer  = layer; }
  void            SetRow(Int_t row) { m_row  = row; }
  void            SetPosition(const TVector3& pos){ m_position.at(0) = pos; }
  void            SetPadLength(Double_t padlength) { m_padlength = padlength; }
  void            SetMRow(Double_t mrow) { m_mrow = mrow; }
  void            SetPadTheta(Double_t mpadtheta) { m_padtheta = mpadtheta; }
  void            SetParentCluster(TPCCluster* parent){ m_parent_cluster = parent; }
  void            SetHoughFlag(Int_t hough_flag) { m_hough_flag = hough_flag; }
  void            SetHoughYnum(Int_t houghY_num);
  void            SetHoughDist(Double_t hough_dist) { m_hough_dist = hough_dist; }
  void            SetHoughDistY(Double_t hough_disty) { m_hough_disty = hough_disty; }
protected:
  void ClearRegisteredHits();
};

//_____________________________________________________________________________
inline TString
TPCHit::ClassName()
{
  static TString s_name("TPCHit");
  return s_name;
}

//_____________________________________________________________________________
inline std::ostream&
operator <<(std::ostream& ost, const TPCHit& hit)
{
  hit.Print("", ost);
  return ost;
}

#endif

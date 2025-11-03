// -*- C++ -*-

#ifndef DC_LOCAL_TRACK_HH
#define DC_LOCAL_TRACK_HH

#include <vector>
#include <functional>

#include <std_ostream.hh>

#include "ThreeVector.hh"
#include "DCLTrackHit.hh"
#include "DetectorID.hh"

class DCLTrackHit;
class DCAnalyzer;

//_____________________________________________________________________________
class DCLocalTrack
{
public:
  static const TString& ClassName();
  explicit DCLocalTrack();
  ~DCLocalTrack();

private:
  DCLocalTrack(const DCLocalTrack &);
  DCLocalTrack & operator =(const DCLocalTrack &);

private:
  Bool_t  m_is_fitted;     // flag of DoFit()
  Bool_t  m_is_calculated; // flag of Calculate()
  Bool_t  m_is_bcsdc;
  std::vector<DCLTrackHit*> m_hit_array;
  std::vector<DCLTrackHit*> m_hit_arrayUV;
  Double_t m_Ax;
  Double_t m_Ay;
  Double_t m_Au;
  Double_t m_Av;
  Double_t m_Chix;
  Double_t m_Chiy;
  Double_t m_Chiu;
  Double_t m_Chiv;
  Double_t m_x0;
  Double_t m_y0;
  Double_t m_u0;
  Double_t m_v0;
  Double_t m_a;
  Double_t m_b;
  Double_t m_chisqr;
  Bool_t   m_good_for_tracking;
  // for SSD
  Double_t m_de;
  // for Honeycomb
  Double_t m_chisqr1st; // 1st iteration for honeycomb
  Double_t m_n_iteration;

public:
  void         AddHit(DCLTrackHit *hitp);
  void         AddHitUV(DCLTrackHit *hitp);
  void         Calculate();
  void         DeleteNullHit();
  Bool_t       DoFit();
  Bool_t       DoFitBcSdc();
  Bool_t       FindLayer(Int_t layer) const;
  Int_t        GetNDF() const;
  Int_t        GetNHit() const { return m_hit_array.size();  }
  Int_t        GetNHitUV()const { return m_hit_arrayUV.size();}
  Int_t        GetNHitSFT() const;
  Int_t        GetNHitY() const;
  DCLTrackHit* GetHit(Int_t nth) const;
  const std::vector<DCLTrackHit*>& GetHitArray() const { return m_hit_array; }
  DCLTrackHit* GetHitUV(Int_t nth) const;
  DCLTrackHit* GetHitOfLayerNumber(Int_t lnum) const;
  Double_t     GetWire(Int_t layer) const;
  Bool_t       HasHoneycomb() const;
  Bool_t       IsFitted() const { return m_is_fitted; }
  Bool_t       IsCalculated() const { return m_is_calculated; }

  void SetAx(Double_t Ax) { m_Ax = Ax; }
  void SetAy(Double_t Ay) { m_Ay = Ay; }
  void SetAu(Double_t Au) { m_Au = Au; }
  void SetAv(Double_t Av) { m_Av = Av; }
  void SetChix(Double_t Chix) { m_Chix = Chix; }
  void SetChiy(Double_t Chiy) { m_Chiy = Chiy; }
  void SetChiu(Double_t Chiu) { m_Chiu = Chiu; }
  void SetChiv(Double_t Chiv) { m_Chiv = Chiv; }
  void SetDe(Double_t de) { m_de = de; }

  Double_t GetX0() const { return m_x0; }
  Double_t GetY0() const { return m_y0; }
  Double_t GetU0() const { return m_u0; }
  Double_t GetV0() const { return m_v0; }

  //For XUV Tracking
  Bool_t DoFitVXU();

  Double_t GetVXU_A() const { return m_a; }
  Double_t GetVXU_B() const { return m_b; }
  Double_t GetVXU(Double_t z) const { return m_a*z+m_b; }
  Double_t GetAx() const { return m_Ax; }
  Double_t GetAy() const { return m_Ay; }
  Double_t GetAu() const { return m_Au; }
  Double_t GetAv() const { return m_Av; }

  Double_t GetDifVXU() const ;
  Double_t GetDifVXUSDC34() const;
  Double_t GetChiSquare() const { return m_chisqr; }
  Double_t GetChiSquare1st() const { return m_chisqr1st; }
  Double_t GetChiX() const { return m_Chix; }
  Double_t GetChiY() const { return m_Chiy; }
  Double_t GetChiU() const { return m_Chiu; }
  Double_t GetChiV() const { return m_Chiv; }
  Double_t GetX(Double_t z) const { return m_x0+m_u0*z; }
  Double_t GetY(Double_t z) const { return m_y0+m_v0*z; }
  Double_t GetS(Double_t z, Double_t tilt) const
  { return GetX(z)*TMath::Cos(tilt)+GetY(z)*TMath::Sin(tilt); }
  Int_t    GetNIteration() const { return m_n_iteration; }
  Double_t GetPhi() const;
  Double_t GetTheta() const;
  Bool_t   GoodForTracking() const { return m_good_for_tracking; }
  Bool_t   GoodForTracking(Bool_t status)
  { Bool_t ret = m_good_for_tracking; m_good_for_tracking = status; return ret; }
  Bool_t   ReCalc(Bool_t ApplyRecursively=false);
  Double_t GetDe() const { return m_de; }
  void   Print(const TString& arg="", std::ostream& ost=hddaq::cout) const;
  void   PrintVXU(const TString& arg="") const;
};

//_____________________________________________________________________________
inline const TString&
DCLocalTrack::ClassName()
{
  static TString s_name("DCLocalTrack");
  return s_name;
}

//_____________________________________________________________________________
inline
std::ostream&
operator <<(std::ostream& ost,
            const DCLocalTrack& track)
{
  track.Print("", ost);
  return ost;
}

//_____________________________________________________________________________
struct DCLTrackComp
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, Bool_t>
{
  Bool_t operator()(const DCLocalTrack * const p1,
                    const DCLocalTrack * const p2) const
    {
      Int_t n1=p1->GetNHit(), n2=p2->GetNHit();
      Double_t chi1=p1->GetChiSquare(), chi2=p2->GetChiSquare();
      if(n1>n2+1)
        return true;
      if(n2>n1+1)
        return false;
      if(n1<=4 || n2<=4)
        return (n1 >= n2);
      if(n1==n2)
        return (chi1 <= chi2);

      return (chi1-chi2 <= 3./(n1-4));// 3-sigma
    }
};

//_____________________________________________________________________________
struct DCLTrackComp1
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, Bool_t>
{
  Bool_t operator()(const DCLocalTrack * const p1,
                    const DCLocalTrack * const p2) const
    {
      Int_t n1=p1->GetNHit(), n2=p2->GetNHit();
      Double_t chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();
      if(n1>n2) return true;
      if(n2>n1) return false;
      return (chi1<=chi2);
    }

};

//_____________________________________________________________________________
struct DCLTrackComp2
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, Bool_t>
{
  Bool_t operator()(const DCLocalTrack * const p1,
                    const DCLocalTrack * const p2) const
    {
      Int_t n1=p1->GetNHit(), n2=p2->GetNHit();
      Double_t chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();
      if(n1<n2) return true;
      if(n2<n1) return false;
      return (chi1<=chi2);
    }

};

//_____________________________________________________________________________
struct DCLTrackComp3
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, Bool_t>
{
  Bool_t operator()(const DCLocalTrack * const p1,
                    const DCLocalTrack * const p2) const
    {
      Int_t n1=p1->GetNHit(), n2=p2->GetNHit();
      Double_t chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();
      Double_t a1= std::abs(1.-chi1), a2=std::abs(1.-chi2);
      if(a1<a2) return true;
      if(a2<a1) return false;
      return (n1<=n2);
    }

};

//_____________________________________________________________________________
struct DCLTrackComp4
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, Bool_t>
{
  Bool_t operator()(const DCLocalTrack * const p1,
                    const DCLocalTrack * const p2) const
    {
      Int_t n1=p1->GetNHit(), n2=p2->GetNHit();
      Double_t chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();
      if((n1>n2+1) && (std::abs(chi1-chi2)<2.))
        return true;
      if((n2>n1+1) && (std::abs(chi1-chi2)<2.))
        return false;

      return (chi1<=chi2);
    }
};

//_____________________________________________________________________________

struct DCLTrackComp_Nhit
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, Bool_t>
{
  Bool_t operator()(const DCLocalTrack * const p1,
                    const DCLocalTrack * const p2) const
    {
      Int_t n1=p1->GetNHit(), n2=p2->GetNHit();

      if(n1>=n2)
        return true;
      else
        return false;
    }
};

//_____________________________________________________________________________

struct DCLTrackComp_Chisqr
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, Bool_t>
{
  Bool_t operator()(const DCLocalTrack * const p1,
                    const DCLocalTrack * const p2) const
    {
      Double_t chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();

      if (chi1 <= chi2)
        return true;
      else
        return false;
    }
};


//_____________________________________________________________________________
struct DCLTrackCompSdcInFiber
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, Bool_t>
{
  Bool_t operator()(const DCLocalTrack * const p1,
                    const DCLocalTrack * const p2) const
    {
      Int_t n1=p1->GetNHit(), n2=p2->GetNHit();
      Double_t chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();
      Int_t NofFiberHit1 = 0;
      Int_t NofFiberHit2 = 0;
      for(Int_t ii=0;ii<n1;ii++){
        Int_t layer = p1->GetHit(ii)->GetLayer();
        if(layer > 6) NofFiberHit1++;
      }
      for(Int_t ii=0;ii<n2;ii++){
        Int_t layer = p2->GetHit(ii)->GetLayer();
        if(layer > 6) NofFiberHit2++;
      }

      if((n1>n2+1)){
        return true;
      }
      else if((n2>n1+1) ){
        return false;
      }
      else if(NofFiberHit1 > NofFiberHit2){
        return true;
      }
      else if(NofFiberHit2 > NofFiberHit1){
        return false;
      }
      else{
        return (chi1<=chi2);
      }

    }
};

//_____________________________________________________________________________
struct DCLTrackCompSdcOut
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, Bool_t>
{
  Bool_t operator()(const DCLocalTrack * const p1,
                    const DCLocalTrack * const p2) const
    {
      Int_t n1=p1->GetNHit(), n2=p2->GetNHit();
      Double_t chi1=p1->GetChiSquare(), chi2=p2->GetChiSquare();
      if(n1 > n2) return true;
      if(n2 > n1) return false;
      return (chi1 <= chi2);
    }
};

#endif

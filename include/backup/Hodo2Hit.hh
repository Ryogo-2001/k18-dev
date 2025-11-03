// -*- C++ -*-

#ifndef HODO2_HIT_HH
#define HODO2_HIT_HH

#include "HodoVHit.hh"

#include <cmath>
#include <cstddef>

#include "HodoRawHit.hh"
#include "ThreeVector.hh"

//_____________________________________________________________________________
class Hodo2Hit : public HodoVHit
{
public:
  explicit Hodo2Hit(HodoRawHit *rhit, Double_t max_time_diff=10.);
  virtual ~Hodo2Hit();
  static TString ClassName();

private:
  Hodo2Hit(const Hodo2Hit&);
  Hodo2Hit& operator =(const Hodo2Hit&);

protected:
  HodoRawHit*            m_raw;
  Bool_t                 m_is_calculated;
  Bool_t                 m_is_tof;
  Double_t               m_max_time_diff; // unit:ns
  std::vector<Double_t>  m_a1;
  std::vector<Double_t>  m_a2;

  struct data_pair
  {
    Double_t time1;
    Double_t time2;
    Double_t ctime1;
    Double_t ctime2;
  };

  std::vector<data_pair> m_pair_cont;
  std::vector<Bool_t>    m_flag_join;

public:
  HodoRawHit* GetRawHit()           { return m_raw; }
  virtual
  Int_t       DetectorId()    const { return m_raw->DetectorId(); }
  virtual
  Int_t       PlaneId()       const { return m_raw->PlaneId(); }
  virtual
  Int_t       SegmentId()     const { return m_raw->SegmentId(); }

  Bool_t      Calculate();
  Bool_t      IsCalculated()  const { return m_is_calculated; }
  void        MakeAsTof()           { m_is_tof = true;}

  Int_t         GetNumOfHit()   const { return m_pair_cont.size(); }
  Double_t      GetAUp(Int_t n=0)     const { return m_a1.at(n); }
  Double_t      GetALeft(Int_t n=0)   const { return m_a1.at(n); }
  Double_t      GetADown(Int_t n=0)   const { return m_a2.at(n); }
  Double_t      GetARight(Int_t n=0)  const { return m_a2.at(n); }
  virtual
  Double_t      GetTUp(Int_t n=0)     const { return m_pair_cont[n].time1; }
  Double_t      GetTLeft(Int_t n=0)   const { return m_pair_cont[n].time1; }
  virtual
  Double_t      GetTDown(Int_t n=0)   const { return m_pair_cont[n].time2; }
  Double_t      GetTRight(Int_t n=0)  const { return m_pair_cont[n].time2; }
  Double_t      GetCTUp(Int_t n=0)    const { return m_pair_cont[n].ctime1; }
  Double_t      GetCTLeft(Int_t n=0)  const { return m_pair_cont[n].ctime1; }
  Double_t      GetCTDown(Int_t n=0)  const { return m_pair_cont[n].ctime2; }
  Double_t      GetCTRight(Int_t n=0) const { return m_pair_cont[n].ctime2; }

  virtual
  Double_t      MeanTime(Int_t n=0)   const { return 0.5*(m_pair_cont[n].time1 + m_pair_cont[n].time2); }
  virtual
  Double_t      CMeanTime(Int_t n=0)  const { return 0.5*(m_pair_cont[n].ctime1 + m_pair_cont[n].ctime2); }
  virtual
  Double_t      DeltaE(Int_t n=0)     const { return std::sqrt(std::abs(m_a1.at(n)*m_a2.at(n))); }
  virtual
	Double_t      UDeltaE(Int_t n=0)     const { return m_a1.at(n); }
	virtual
	Double_t      DDeltaE(Int_t n=0)     const { return m_a2.at(n); }
  Double_t      TimeDiff(Int_t n=0)   const { return m_pair_cont[n].ctime2 - m_pair_cont[n].ctime1; }

  void          SetJoined(Int_t m)          { m_flag_join.at(m) = true; }
  Bool_t        Joined(Int_t m)       const { return m_flag_join.at(m); }
  Bool_t        JoinedAllMhit();

  virtual Bool_t ReCalc(Bool_t applyRecursively=false)
  { return Calculate(); }

};

//_____________________________________________________________________________
inline TString
Hodo2Hit::ClassName()
{
  static TString s_name("Hodo2Hit");
  return s_name;
}

#endif

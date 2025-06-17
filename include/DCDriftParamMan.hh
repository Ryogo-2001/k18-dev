// -*- C++ -*-

#ifndef DC_DRIFT_PARAM_MAN_HH
#define DC_DRIFT_PARAM_MAN_HH

#include <map>
#include <vector>
#include <TString.h>

struct DCDriftParamRecord;

//_____________________________________________________________________________
struct DCDriftParamRecord
{
  Int_t type, np;
  std::vector<Double_t> param;
  DCDriftParamRecord(Int_t t, Int_t n, const std::vector<Double_t>& p)
    : type(t), np(n), param(p)
    {}
};

//_____________________________________________________________________________
class DCDriftParamMan
{
public:
  static const TString&   ClassName();
  static DCDriftParamMan& GetInstance();
  ~DCDriftParamMan();

private:
  DCDriftParamMan();
  DCDriftParamMan(const DCDriftParamMan&);
  DCDriftParamMan& operator =(const DCDriftParamMan&);

private:
  typedef std::map<Int_t, DCDriftParamRecord*> DCDriftContainer;
  typedef DCDriftContainer::const_iterator DCDriftIterator;
  Bool_t           m_is_ready;
  TString          m_file_name;
  DCDriftContainer m_container;

public:
  Bool_t CalcDrift(Int_t PlaneId, Double_t WireId, Double_t ctime,
                   Double_t & dt, Double_t & dl) const;
  Bool_t Initialize();
  Bool_t Initialize(const TString& file_name);
  Bool_t IsReady() const { return m_is_ready; }
  void   SetFileName(const TString& file_name) { m_file_name = file_name; }

private:
  void                ClearElements();
  static Double_t     DriftLength1(Double_t dt, Double_t vel);
  static Double_t     DriftLength2(Double_t dt,
                                   Double_t p1, Double_t p2, Double_t p3,
                                   Double_t st, Double_t p5, Double_t p6);
  static Double_t     DriftLength3(Double_t dt, Double_t p1, Double_t p2,
                                   Int_t PlaneId);
  static Double_t     DriftLength4(Double_t dt,
                                   Double_t p1, Double_t p2, Double_t p3);
  static Double_t     DriftLength5(Double_t dt,
                                   Double_t p1, Double_t p2, Double_t p3,
                                   Double_t p4, Double_t p5);
  static Double_t     DriftLength6(Int_t PlaneId,
                                   Double_t dt,
                                   Double_t p1, Double_t p2, Double_t p3,
                                   Double_t p4, Double_t p5);
  static Double_t     DriftLength7(Int_t PlaneId,
                                   Double_t dt,
                                   Double_t p1, Double_t p2, Double_t p3,
                                   Double_t p4, Double_t p5, Double_t p6);
  DCDriftParamRecord* GetParameter(Int_t PlaneId, Double_t WireId) const;
};

//_____________________________________________________________________________
inline const TString&
DCDriftParamMan::ClassName()
{
  static TString s_name("DCDriftParamMan");
  return s_name;
}

//_____________________________________________________________________________
inline DCDriftParamMan&
DCDriftParamMan::GetInstance()
{
  static DCDriftParamMan s_instance;
  return s_instance;
}

#endif

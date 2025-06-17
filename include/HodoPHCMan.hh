// -*- C++ -*-

#ifndef HODO_PHC_MAN_HH
#define HODO_PHC_MAN_HH

#include <map>
#include <vector>
#include <TString.h>

//_____________________________________________________________________________
class HodoPHCParam
{
public:
  static const TString& ClassName();
  HodoPHCParam(Int_t type, Int_t n_param,
               const std::vector<Double_t>& parlist);
  ~HodoPHCParam();

private:
  HodoPHCParam(const HodoPHCParam&);
  HodoPHCParam& operator =(const HodoPHCParam&);

private:
  Int_t                 m_type;
  Int_t                 m_n_param;
  std::vector<Double_t> m_param_list;

public:
  Double_t DoPHC(Double_t time, Double_t de) const;
  Double_t DoRPHC(Double_t time, Double_t de) const;
  Double_t DoSTC(Double_t stof, Double_t btof) const; //for stof correction

private:
  Double_t Type1Correction(Double_t time, Double_t de) const;
  Double_t Type2Correction(Double_t time, Double_t de) const; // For fiber
  Double_t Type1RCorrection(Double_t time, Double_t de) const;
  Double_t Type1STCorrection(Double_t stof, Double_t btof) const;//for stof correction
};

//_____________________________________________________________________________
inline const TString&
HodoPHCParam::ClassName()
{
  static TString s_name("HodoPHCParam");
  return s_name;
}

//_____________________________________________________________________________
class HodoPHCMan
{
public:
  static const TString& ClassName();
  static HodoPHCMan&    GetInstance();
  ~HodoPHCMan();

private:
  HodoPHCMan();
  HodoPHCMan(const HodoPHCMan&);
  HodoPHCMan& operator =(const HodoPHCMan&);

private:
  typedef std::map<Int_t, HodoPHCParam*> PhcPContainer;
  typedef PhcPContainer::const_iterator  PhcPIterator;
  Bool_t        m_is_ready;
  TString       m_file_name;
  PhcPContainer m_container;

public:
  Bool_t DoCorrection(Int_t cid, Int_t plid, Int_t seg, Int_t ud,
                      Double_t time, Double_t de, Double_t& ctime) const;
  Bool_t DoRCorrection(Int_t cid, Int_t plid, Int_t seg, Int_t ud,
                       Double_t time, Double_t de, Double_t& ctime) const;
  Bool_t DoStofCorrection(Int_t cid, Int_t plid, Int_t seg, Int_t ud,
                          Double_t stof, Double_t btof, Double_t& cstof) const;
  Bool_t Initialize();
  Bool_t Initialize(const TString& file_name);
  Bool_t IsReady() const { return m_is_ready; }
  void   SetFileName(const TString& file_name) { m_file_name = file_name; }

private:
  void          ClearElements();
  HodoPHCParam* GetMap(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const;
};

//_____________________________________________________________________________
inline const TString&
HodoPHCMan::ClassName()
{
  static TString s_name("HodoPHCMan");
  return s_name;
}

//_____________________________________________________________________________
inline HodoPHCMan&
HodoPHCMan::GetInstance()
{
  static HodoPHCMan s_instance;
  return s_instance;
}

#endif

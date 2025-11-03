// -*- C++ -*-

#ifndef DC_TDC_CALIB_MAN_HH
#define DC_TDC_CALIB_MAN_HH

#include <map>
#include <vector>
#include <TString.h>

struct DCTdcCalMap;

//_____________________________________________________________________________
class DCTdcCalibMan
{
public:
  static const TString& ClassName();
  static DCTdcCalibMan& GetInstance();
  ~DCTdcCalibMan();

private:
  DCTdcCalibMan();
  DCTdcCalibMan(const DCTdcCalibMan&);
  DCTdcCalibMan& operator =(const DCTdcCalibMan&);

private:
  typedef std::map<Int_t, DCTdcCalMap*> DCTdcContainer;
  typedef DCTdcContainer::const_iterator DCTdcIterator;
  Bool_t         m_is_ready;
  TString        m_file_name;
  DCTdcContainer m_container;

public:
  Bool_t GetParameter(Int_t plane_id, Double_t wire_id,
                      Double_t &p0, Double_t &p1) const;
  Bool_t GetTime(Int_t plane_id, Double_t wire_id, Int_t tdc,
                 Double_t& time) const;
  Bool_t GetTdc(Int_t plane_id, Double_t wire_id, Double_t time,
                Int_t& tdc) const;
  Bool_t Initialize();
  Bool_t Initialize(const TString& file_name);
  Bool_t IsReady() const { return m_is_ready; }
  void   SetFileName(const TString& file_name) { m_file_name = file_name; }

private:
  void         ClearElements();
  DCTdcCalMap* GetMap(Int_t plane_id, Double_t wire_id) const;
};

//_____________________________________________________________________________
inline DCTdcCalibMan&
DCTdcCalibMan::GetInstance()
{
  static DCTdcCalibMan s_instance;
  return s_instance;
}

//_____________________________________________________________________________
inline const TString&
DCTdcCalibMan::ClassName()
{
  static TString s_name("DCTdcCalibMan");
  return s_name;
}

#endif

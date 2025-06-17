// -*- C++ -*-

#ifndef HODO_PARAM_MAN_HH
#define HODO_PARAM_MAN_HH

#include <map>
#include <TString.h>

//_____________________________________________________________________________
//Hodo TDC to Time
class HodoTParam
{
public:
  HodoTParam(Double_t offset, Double_t gain)
    : m_offset(offset), m_gain(gain)
  {
  }
  ~HodoTParam()
  {
  }

private:
  HodoTParam();
  HodoTParam(const HodoTParam&);
  HodoTParam& operator =(const HodoTParam&);

private:
  Double_t m_offset;
  Double_t m_gain;

public:
  Double_t Offset() const { return m_offset; }
  Double_t Gain() const { return m_gain; }
  Double_t Time(Int_t tdc) const
  { return ((Double_t)tdc - m_offset) * m_gain; }
  Int_t    Tdc(Double_t time) const
  { return (Int_t)(time/m_gain + m_offset); }
};

//_____________________________________________________________________________
//Hodo ADC to DeltaE
class HodoAParam
{
public:
  HodoAParam(Double_t pedestal, Double_t gain)
    : m_pedestal(pedestal), m_gain(gain)
  {
  }
  ~HodoAParam()
  {
  }

private:
  HodoAParam();
  HodoAParam(const HodoAParam&);
  HodoAParam& operator =(const HodoAParam&);

private:
  Double_t m_pedestal;
  Double_t m_gain;

public:
  Double_t Pedestal() const { return m_pedestal; }
  Double_t Gain() const { return m_gain; }
  Double_t DeltaE(Int_t adc) const
  { return ((Double_t)adc - m_pedestal) / (m_gain - m_pedestal); }
  Int_t    Adc(Double_t de) const
  { return (Int_t)(m_gain * de + m_pedestal * (1. - de)); }
};

//_____________________________________________________________________________
//correction parameter for CFT fiber position
class HodoFParam
{
public:
  HodoFParam(Double_t par0, Double_t par1, Double_t par2,
             Double_t par3, Double_t par4, Double_t par5)
    : Par0(par0), Par1(par1), Par2(par2),
      Par3(par3), Par4(par4), Par5(par5)
  {
  }
  ~HodoFParam()
  {
  }

private:
  HodoFParam();
  HodoFParam(const HodoFParam&);
  HodoFParam& operator =(const HodoFParam&);

private:
  Double_t Par0, Par1, Par2, Par3, Par4, Par5;

public:
  Double_t par0() const { return Par0; }
  Double_t par1() const { return Par1; }
  Double_t par2() const { return Par2; }
  Double_t par3() const { return Par3; }
  Double_t par4() const { return Par4; }
  Double_t par5() const { return Par5; }
};

//_____________________________________________________________________________
//HodoParam Main Class
class HodoParamMan
{
public:
  static const TString& ClassName();
  static HodoParamMan&  GetInstance();
  ~HodoParamMan();

private:
  HodoParamMan();
  HodoParamMan(const HodoParamMan&);
  HodoParamMan& operator =(const HodoParamMan&);

private:
  enum eAorT { kAdc, kTdc, kAorT };
  typedef std::map<Int_t, HodoTParam*> TContainer;
  typedef std::map<Int_t, HodoAParam*> AContainer;
  typedef std::map<Int_t, HodoFParam*> FContainer;
  typedef TContainer::const_iterator TIterator;
  typedef AContainer::const_iterator AIterator;
  typedef FContainer::const_iterator FIterator;
  Bool_t     m_is_ready;
  TString    m_file_name;
  TContainer m_TPContainer;
  AContainer m_APContainer;
  FContainer m_FPContainer;

public:
  Bool_t   Initialize();
  Bool_t   Initialize(const TString& file_name);
  Bool_t   IsReady() const { return m_is_ready; }
  Bool_t   GetTime(Int_t cid, Int_t plid, Int_t seg,
                   Int_t ud, Int_t tdc, Double_t &time) const;
  Bool_t   GetDe(Int_t cid, Int_t plid, Int_t seg,
                 Int_t ud, Int_t adc, Double_t &de) const;
  Bool_t   GetTdc(Int_t cid, Int_t plid, Int_t seg,
                  Int_t ud, Double_t time, Int_t &tdc) const;
  Bool_t   GetAdc(Int_t cid, Int_t plid, Int_t seg,
                  Int_t ud, Double_t de, Int_t &adc) const;
  Double_t GetP0(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const;
  Double_t GetP1(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const;
  Double_t GetPar(Int_t cid, Int_t plid, Int_t seg, Int_t ud ,Int_t i) const;
  Double_t GetOffset(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const;
  Double_t GetGain(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const;
  void     SetFileName(const TString& file_name) { m_file_name = file_name; }

private:
  void        ClearACont();
  void        ClearTCont();
  void        ClearFCont();
  HodoTParam* GetTmap(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const;
  HodoAParam* GetAmap(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const;
  HodoFParam* GetFmap(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const;
};

//_____________________________________________________________________________
inline const TString&
HodoParamMan::ClassName()
{
  static TString s_name("HodoParamMan");
  return s_name;
}

//_____________________________________________________________________________
inline HodoParamMan&
HodoParamMan::GetInstance()
{
  static HodoParamMan s_instance;
  return s_instance;
}

#endif

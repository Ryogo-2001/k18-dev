// -*- C++ -*-

#ifndef TAGPL_MATHC_HH
#define TAGPL_MATHC_HH

#include <bitset>
#include <map>

#include <TString.h>

//_____________________________________________________________________________
class TAGPLMatch
{
public:
  static const TString& ClassName();
  static TAGPLMatch&      GetInstance();
  virtual ~TAGPLMatch();

private:
  TAGPLMatch();
  TAGPLMatch(const TAGPLMatch&);
  TAGPLMatch& operator =(const TAGPLMatch&);

private:
  struct Param
  {
    Param();
    ~Param();
    Double_t m_seg;
    Double_t m_xmin;
    Double_t m_xmax;
    void Print() const;
  };

  enum EStatus
  {
    kReady,
    kVerbose,
    kNStatus
  };

  std::bitset<kNStatus>     m_status;
  std::map<Double_t, Param> m_param;

public:
  enum EParam
  {
    kTAGPLSegment,
    kXMin,
    kXMax,
    kNParam
  };
  Bool_t Initialize(const TString& file_name);
  Bool_t Judge(Double_t bft_xpos, Double_t bh1seg);
  void   Print() const;
  void   SetVerbose();
};

//_____________________________________________________________________________
inline const TString&
TAGPLMatch::ClassName()
{
  static TString s_name("TAGPLMatch");
  return s_name;
}

//_____________________________________________________________________________
inline TAGPLMatch&
TAGPLMatch::GetInstance()
{
  static TAGPLMatch s_instance;
  return s_instance;
}

//_____________________________________________________________________________
inline void
TAGPLMatch::SetVerbose()
{
  m_status.set(kVerbose);
}

#endif

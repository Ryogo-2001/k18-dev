// -*- C++ -*-

#ifndef BH1_MATHC_HH
#define BH1_MATHC_HH

#include <vector>
#include <bitset>
#include <TString.h>

#include <std_ostream.hh>

//_____________________________________________________________________________
class BH1Match
{
public:
  static const TString& ClassName();
  static BH1Match&      GetInstance();
  virtual ~BH1Match();

private:
  BH1Match();
  BH1Match(const BH1Match&);
  BH1Match& operator =(const BH1Match&);

private:
  struct Param
  {
    Param();
    ~Param();
    Double_t m_seg;
    Double_t m_xmin;
    Double_t m_xmax;
    void Print(std::ostream& ost=hddaq::cout) const;
  };

  enum EStatus
  {
    kReady,
    kVerbose,
    kNStatus
  };

  std::bitset<kNStatus> m_status;
  std::vector<Param>    m_param;

public:
  enum EParam
  {
    kBH1Segment,
    kXMin,
    kXMax,
    kNParam
  };
  Bool_t Initialize(const TString& file_name);
  Bool_t Judge(Double_t bft_xpos, Double_t bh1seg);
  void   Print(std::ostream& ost=hddaq::cout) const;
  void   SetVerbose();
};

//_____________________________________________________________________________
inline const TString&
BH1Match::ClassName()
{
  static TString s_name("BH1Match");
  return s_name;
}

//_____________________________________________________________________________
inline BH1Match&
BH1Match::GetInstance()
{
  static BH1Match s_instance;
  return s_instance;
}

//_____________________________________________________________________________
inline void
BH1Match::SetVerbose()
{
  m_status.set(kVerbose);
}

#endif

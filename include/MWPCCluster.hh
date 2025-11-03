// -*- C++ -*-

#ifndef MWPC_CLUSTER_HH
#define MWPC_CLUSTER_HH

#include <vector>

#include <TString.h>

class DCHit;

//_____________________________________________________________________________
class MWPCCluster
{
public:
  static const TString& ClassName();
  MWPCCluster();
  ~MWPCCluster();

private:
  MWPCCluster(const MWPCCluster&);
  MWPCCluster& operator =(const MWPCCluster&);

public:
  struct Statistics
  {
    Double_t m_wire;
    Double_t m_wpos;
    Double_t m_leading;
    Double_t m_trailing;
    Double_t m_length;
    Double_t m_totalLength;
    Int_t    m_clusterSize;

    Statistics();
    ~Statistics();

    void Print(const TString& arg="") const;
  };

private:
  std::vector<DCHit*> m_hits;
  Statistics          m_mean;
  Statistics          m_first;
  bool                m_status;

public:
  void                       Add(DCHit* h);
  void                       Calculate();
  Int_t                      GetClusterSize() const;
  const Statistics&          GetFirst() const;
  const std::vector<DCHit*>& GetHits() const;
  const Statistics&          GetMean() const;
  Double_t                   GetMeanTime() const;
  Double_t                   GetMeanWire() const;
  Double_t                   GetMeanWirePos() const;
  Int_t                      GetNumOfHits() const;
  bool                       IsGoodForAnalysis() const;
  void                       Print(const TString& arg="") const;
  void                       SetStatus(bool status);

};

//_____________________________________________________________________________
inline const TString&
MWPCCluster::ClassName()
{
  static TString s_name("MWPCCluster");
  return s_name;
}

#endif

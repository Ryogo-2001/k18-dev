// -*- C++ -*-

#include "MWPCCluster.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <set>

#include <std_ostream.hh>

#include "DCHit.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "DeleteUtility.hh"

//_____________________________________________________________________________
void
calcFirst(const std::vector<DCHit*>& hits,
          MWPCCluster::Statistics& first)
{
  Int_t nHits = hits.size();
  if (nHits==0)
    return;

  for (Int_t i=0; i<nHits; ++i)
  {
    const DCHit* h = hits[i];
    if (!h)
      continue;
    Double_t le = h->GetDriftTime(0);
    if (i==0)
    {
      first.m_wire     = h->GetWire();
      first.m_wpos     = h->GetWirePosition();
      first.m_leading  = le;
      first.m_trailing = h->GetTrailingTime(0);
      first.m_length   = first.m_trailing - le;
    }
    else if (le<first.m_leading) // update first hit
    {
      first.m_wire     = h->GetWire();
      first.m_wpos     = h->GetWirePosition();
      first.m_leading  = le;
      first.m_trailing = h->GetTrailingTime(0);
      first.m_length   = first.m_trailing - le;
    }
  }

  Double_t wire     = 0;
  Double_t wire_pos = 0;
  Double_t trailing = 0;
  Double_t length   = 0;

  std::set<Int_t> wires;
  for (Int_t i=0; i<nHits; ++i)
  {
    const DCHit* h = hits[i];
    if (!h)
      continue;
    Double_t le = h->GetDriftTime(0);
    if (MathTools::Equal(le, first.m_leading))
    {
      Double_t tr   = h->GetTrailingTime();
      Double_t len  = tr - le;
      Double_t w    = h->GetWire();
      wire       += w*len;
      wire_pos   += h->GetWirePosition()*len;
      trailing   += tr*len;
      length     += len;
      wires.insert(static_cast<Int_t>(w));
    }
  }

  if (!wires.empty())
  {
    first.m_wire        = wire/length;
    first.m_wpos        = wire_pos/length;
    first.m_trailing    = trailing/length;
    first.m_length      = length/wires.size();
    first.m_totalLength = length;
    std::set<Int_t>::iterator imin
      = std::min_element(wires.begin(), wires.end());
    std::set<Int_t>::iterator imax
      = std::max_element(wires.begin(), wires.end());
    first.m_clusterSize = *imax - *imin + 1;
  }
  else
  {
    first.m_totalLength = first.m_length;
    first.m_clusterSize = 1;
  }

  return;
}

//_____________________________________________________________________________
void
calcMean(const std::vector<DCHit*>& hits,
         MWPCCluster::Statistics& mean)
{
  Int_t nHits = hits.size();
  if (nHits==0)
    return;
  Double_t wire     = 0;
  Double_t wire_pos = 0;
  Double_t leading  = 0;
  Double_t trailing = 0;
  Double_t length   = 0;

  std::set<Int_t> wires;

  for (Int_t i=0; i<nHits; ++i)
  {
    const DCHit* h = hits[i];
    if (!h)
      continue;
    Double_t w    = h->GetWire();
    Double_t wpos = h->GetWirePosition();
    Double_t le   = h->GetDriftTime(0);
    Double_t tr   = h->GetTrailingTime(0);
    Double_t len  = tr - le;

    length   += len;
    leading  += le*len;
    trailing += tr*len;
    wire     += w*len;
    wire_pos += wpos*len;
    wires.insert(static_cast<Int_t>(w));
  }

  mean.m_wire     = wire    /length;
  mean.m_wpos     = wire_pos/length;
  mean.m_leading  = leading /length;
  mean.m_trailing = trailing/length;
  mean.m_length   = length/wires.size();
  mean.m_totalLength = length;

  std::set<Int_t>::iterator imin = std::min_element(wires.begin(), wires.end());
  std::set<Int_t>::iterator imax = std::max_element(wires.begin(), wires.end());
  mean.m_clusterSize = *imax - *imin + 1;
  return;
}

//_____________________________________________________________________________
// struct MWPCCluster::Statistics
//_____________________________________________________________________________
MWPCCluster::
Statistics::Statistics()
  : m_wire(-0xffff),
    m_wpos(-0xffff),
    m_leading(-0xffff),
    m_trailing(-0xffff),
    m_length(-0xffff),
    m_totalLength(-0xffff),
    m_clusterSize(-0xffff)
{
}

//_____________________________________________________________________________
MWPCCluster::
Statistics::~Statistics()
{
}

//_____________________________________________________________________________
void
MWPCCluster::
Statistics::Print(const TString& arg) const
{
  hddaq::cout << FUNC_NAME << " " << arg << std::endl
	      << "  wire   = " << m_wire << std::endl
	      << "  wpos   = " << m_wpos        << " [mm]" << std::endl
	      << "  dt     = " << m_leading     << " [nsec]" << std::endl
	      << "  trail  = " << m_trailing    << " [nsec]" << std::endl
	      << "  siglen = " << m_length      << " [nsec]" << std::endl
	      << "  total  = " << m_totalLength << " [nsec]" << std::endl
	      << "  csize  = " << m_clusterSize << std::endl;
}

//_____________________________________________________________________________
// class MWPCCluster
//_____________________________________________________________________________
MWPCCluster::MWPCCluster()
  : m_hits(0),
    m_mean(),
    m_first(),
    m_status(false)
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
MWPCCluster::~MWPCCluster()
{
  del::DeleteObject(m_hits);
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
MWPCCluster::Add(DCHit* h)
{
  m_hits.push_back(h);
  return;
}

//_____________________________________________________________________________
void
MWPCCluster::Calculate()
{
  if (m_hits.empty())
    return;

  calcMean(m_hits, m_mean);
  calcFirst(m_hits, m_first);

  return;
}

//_____________________________________________________________________________
Int_t
MWPCCluster::GetClusterSize() const
{
  return m_mean.m_clusterSize;
}

//_____________________________________________________________________________
const MWPCCluster::Statistics&
MWPCCluster::GetFirst() const
{
  return m_first;
}

//_____________________________________________________________________________
const std::vector<DCHit*>&
MWPCCluster::GetHits() const
{
  return m_hits;
}

//_____________________________________________________________________________
const MWPCCluster::Statistics&
MWPCCluster::GetMean() const
{
  return m_mean;
}

//_____________________________________________________________________________
Double_t
MWPCCluster::GetMeanTime() const
{
  return m_mean.m_leading;
}

//_____________________________________________________________________________
Double_t
MWPCCluster::GetMeanWire() const
{
  return m_mean.m_wire;
}

//_____________________________________________________________________________
Double_t
MWPCCluster::GetMeanWirePos() const
{
  return m_mean.m_wpos;
}

//_____________________________________________________________________________
Int_t
MWPCCluster::GetNumOfHits() const
{
  return m_hits.size();
}

//_____________________________________________________________________________
void
MWPCCluster::Print(const TString& arg) const
{
  hddaq::cout << FUNC_NAME << " " << arg << std::endl;
  hddaq::cout << " nhits = " << m_hits.size() << std::endl;
  m_mean.Print("mean");
  m_first.Print("first");
  return;
}

//_____________________________________________________________________________
bool
MWPCCluster::IsGoodForAnalysis() const
{
  return m_status;
}

//_____________________________________________________________________________
void
MWPCCluster::SetStatus(bool status)
{
  m_status = status;
  return;
}

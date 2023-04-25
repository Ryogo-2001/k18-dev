// -*- C++ -*-

#ifndef DST_HELPER_HH
#define DST_HELPER_HH

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>

#include <filesystem_util.hh>

#include "DCAnalyzer.hh"

// if event number mismatch is found, exit process.
#define CheckEventNumberMismatch 1

//______________________________________________________________________________
namespace dst
{
// implemented in each dst
extern std::vector<TString> ArgName;
extern std::vector<TString> TreeName;
extern std::vector<TFile*>  TFileCont;
extern std::vector<TTree*>  TTreeCont;
extern std::vector<TTreeReader*> TTreeReaderCont;
Bool_t InitializeEvent();
Bool_t DstOpen(std::vector<std::string> arg);
Bool_t DstRead();
Bool_t DstRead(Int_t ievent);
Bool_t DstRead(Int_t ievent, DCAnalyzer *DCAna);
Bool_t DstClose();

//______________________________________________________________________________
inline Bool_t
CheckArg(const std::vector<std::string>& arg)
{
  const Int_t n = arg.size();
  Bool_t status = (n == ArgName.size() && n == TreeName.size());

  if(!status){
    std::cout << "#D Usage : " << hddaq::basename(arg[0]);
    for(Int_t i=1; i<ArgName.size(); ++i){
      std::cout << " " << ArgName[i];
    }
    std::cout << std::endl;
    std::cout << " ArgName.size() = " << ArgName.size() << " "
              << " TreeName.size() = " << TreeName.size() << std::endl;
    return false;
  }

  for(Int_t i=0; i<n; ++i){
    std::cout << " key = " << std::setw(18) << std::left << ArgName[i]
              << " arg[" << i << "] = " << arg[i] << std::endl;
  }

  TFileCont.resize(n); TTreeCont.resize(n); TTreeReaderCont.resize(n);
  return (ArgName.size()==TreeName.size());
}

//______________________________________________________________________________
inline Bool_t
OpenFile(TFile*& file, const TString& name)
{
  file = new TFile(name);
  if(!file || !file->IsOpen()){
    std::cerr << "#E failed to open TFile : " << name << std::endl;
    return false;
  }
  return true;
}

//______________________________________________________________________________
inline Bool_t
OpenTree(TFile* file, TTree*& tree, const TString& name)
{
  if(!file || !file->IsOpen()) return false;
  tree = (TTree*)file->Get(name);
  if(!tree){
    std::cerr << "#E failed to open TTree : " << name << std::endl;
    return false;
  }
  return true;
}

//______________________________________________________________________________
inline Bool_t
CheckEntries(const std::vector<TTree*>& TTreeCont)
{
  const Int_t n = TTreeCont.size();
  Bool_t status = true;
  std::vector<Int_t> entries(n, -1);
  for(Int_t i=0; i<n; ++i){
    if(!TTreeCont[i]) continue;
    entries[i] = TTreeCont[i]->GetEntries();
    if(i>0 && entries[i]!=entries[i-1] && entries[i-1]!=-1){
      status = false;
    }
  }
  if(!status){
#if CheckEventNumberMismatch
    std::cerr << "#E Entries Mismatch" << std::endl;
#else
    std::cerr << "#W Entries Mismatch" << std::endl;
#endif
    for(Int_t i=0; i<n; ++i){
      if(!TTreeCont[i]) continue;
      std::cerr << "   " << std::setw(8) << TTreeCont[i]->GetName()
                << " " << entries[i] << std::endl;
    }
  }
#if CheckEventNumberMismatch
  return status;
#else
  return true;
#endif
}

//______________________________________________________________________________
inline Int_t
GetEntries(const std::vector<TTree*>& TTreeCont)
{
  std::vector<Int_t> nevent;
  for(Int_t i=0, n=TTreeCont.size(); i<n; ++i){
    if(TTreeCont[i]){
#if CheckEventNumberMismatch
      return TTreeCont[i]->GetEntries();
#else
      nevent.push_back(TTreeCont[i]->GetEntries());
#endif
    }
  }
#if CheckEventNumberMismatch
  return 0;
#else
  return *std::min_element(nevent.begin(), nevent.end());
#endif
}

//______________________________________________________________________________
inline Bool_t
GetEntry(Int_t ievent)
{
  for(Int_t i=0, n=TTreeCont.size(); i<n; ++i){
    if(TTreeCont[i]){
      TTreeCont[i]->GetEntry(ievent);
      if(TTreeReaderCont[i]){
        TTreeReaderCont[i]->SetEntry(ievent);
      }
    }
  }
  return true;
}
}

#endif

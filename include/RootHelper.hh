// -*- C++ -*-

#ifndef ROOT_HELPER_HH
#define ROOT_HELPER_HH

#include <stdexcept>

#include <TApplication.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TDatabasePDG.h>
#include <TError.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH2Poly.h>
#include <TList.h>
#include <TLatex.h>
#include <TMath.h>
#include <TMatrix.h>
#include <TMinuit.h>
#include <TObject.h>
#include <TPDGCode.h>
#include <TParticlePDG.h>
#include <TPad.h>
#include <TProfile.h>
#include <TRint.h>
#include <TRandom.h>
// #include <TRandom1.h>
// #include <TRandom2.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TVector3.h>

#include <std_ostream.hh>

#include "Exception.hh"

#define ThrowError 1 // if h[i] already exist, throw error
#define OverWrite  0 // if h[i] already exist, delete and renew.

namespace root
{
const Int_t MaxHits  = 500;
const Int_t MaxDepth = 16;
const Int_t MaxHist  = 100000000;
extern TH1   *h[MaxHist];
extern TTree *tree;

//_____________________________________________________________________________
inline void
HB1(Int_t i, const Char_t* title,
    Int_t nbinx, Double_t xlow, Double_t xhigh)
{
  if(i<0 || MaxHist<=i)
    throw Exception(Form("HB1() invalid HistId : %d/%d", i, MaxHist));
  if(h[i]){
#if ThrowError
    throw Exception(Form("h%d (%s) is already exist", i, title));
#endif
#if OverWrite
    delete h[i];
    h[i] = nullptr;
#endif
  }
  h[i] = new TH1D(Form("h%d", i), title, nbinx, xlow, xhigh);
}

//_____________________________________________________________________________
inline void
HB1(Int_t i, const TString& title,
    Int_t nbinx, Double_t xlow, Double_t xhigh)
{
  HB1(i, title.Data(), nbinx, xlow, xhigh);
}

//_____________________________________________________________________________
inline void
HB2(Int_t i, const Char_t* title,
    Int_t nbinx, Double_t xlow, Double_t xhigh,
    Int_t nbiny, Double_t ylow, Double_t yhigh)
{
  if(i<0 || MaxHist<=i)
    throw Exception(Form("HB2() invalid HistId : %d/%d", i, MaxHist));
  if(h[i]){
#if ThrowError
    throw Exception(Form("h%d (%s) is already exist", i, title));
#endif
#if OverWrite
    delete h[i];
    h[i] = nullptr;
#endif
  }
  h[i] = new TH2D(Form("h%d", i), title,
                  nbinx, xlow, xhigh,
                  nbiny, ylow, yhigh);
}

//_____________________________________________________________________________
inline void
HB2(Int_t i, const TString& title,
    Int_t nbinx, Double_t xlow, Double_t xhigh,
    Int_t nbiny, Double_t ylow, Double_t yhigh)
{
  HB2(i, title.Data(), nbinx, xlow, xhigh, nbiny, ylow, yhigh);
}

//_____________________________________________________________________________
inline void
HBProf(Int_t i, const Char_t* title,
       Int_t nbinx, Double_t xlow, Double_t xhigh,
       Double_t ylow, Double_t yhigh)
{
  if(i<0 || MaxHist<=i)
    throw Exception(Form("HBProf() invalid HistId : %d/%d", i, MaxHist));
  if(h[i]){
#if ThrowError
    throw Exception(Form("h%d (%s) is already exist", i, title));
#endif
#if OverWrite
    delete h[i];
    h[i] = nullptr;
#endif
  }
  h[i] = new TProfile(Form("h%d", i), title,
                      nbinx, xlow, xhigh, ylow, yhigh);
}

//_____________________________________________________________________________
inline void
HB2Poly(Int_t i, const Char_t* title,
        Double_t xmin=-300., Double_t xmax=300.,
        Double_t ymin=-300., Double_t ymax=300.)
{
  if(i<0 || MaxHist<=i)
    throw Exception(Form("HB2Poly() invalid HistId : %d/%d", i, MaxHist));
  if(h[i]){
#if ThrowError
    throw Exception(Form("h%d (%s) is already exist", i, title));
#endif
#if OverWrite
    delete h[i];
    h[i] = nullptr;
#endif
  }
  h[i] = new TH2Poly(Form("h%d", i), title,
                     xmin, xmax, ymin, ymax);
  gDirectory->Add(h[i]);
  /*
   * Bin is set by tpc::InitializeHistograms() in TPCPadHelper.hh
   */
}

//_____________________________________________________________________________
inline void
HC2Poly(Int_t i)
{
  if(i<0 || MaxHist<=i)
    throw Exception(Form("HC2Poly() invalid HistId : %d/%d", i, MaxHist));
  if(h[i]) dynamic_cast<TH2Poly*>(h[i])->Reset("");
}

//_____________________________________________________________________________
inline void
HF1(Int_t i, Double_t x)
{
  if(i<0 || MaxHist<=i)
    throw Exception(Form("HF1() invalid HistId : %d/%d", i, MaxHist));
  if(h[i]) h[i]->Fill(x);
}

//_____________________________________________________________________________
inline void
HF2(Int_t i, Double_t x, Double_t y)
{
  if(i<0 || MaxHist<=i)
    throw Exception(Form("HF2() invalid HistId : %d/%d", i, MaxHist));
  if(h[i]) h[i]->Fill(x, y);
}

//_____________________________________________________________________________
inline void
HF2Poly(Int_t i, Double_t x, Double_t y, Double_t w=1.)
{
  if(i<0 || MaxHist<=i)
    throw Exception(Form("HF2Poly() invalid HistId : %d/%d", i, MaxHist));
  if(h[i]) dynamic_cast<TH2Poly*>(h[i])->Fill(x, y, w);
}

//_____________________________________________________________________________
inline void
HF2Poly(Int_t i, Int_t bin, Double_t val)
{
  if(i<0 || MaxHist<=i)
    throw Exception(Form("HF2Poly() invalid HistId : %d/%d", i, MaxHist));
  if(h[i]) dynamic_cast<TH2Poly*>(h[i])->SetBinContent(bin, val);
}

//_____________________________________________________________________________
inline void
HFProf(Int_t i, Double_t x, Double_t y)
{
  if(i<0 || MaxHist<=i)
    throw Exception(Form("HFProf() invalid HistId : %d/%d", i, MaxHist));
  if(h[i]) h[i]->Fill(x, y);
}

//_____________________________________________________________________________
inline void
HBTree(const Char_t* name, const Char_t* title)
{
  if(tree){
#if ThrowError
    throw Exception(Form("%s (%s) is already exist", name, title));
#endif
#if OverWrite
    delete tree;
    tree = nullptr;
#endif
  }
  tree = new TTree(name, title);
}

//_____________________________________________________________________________
inline void
HPrint()
{
  hddaq::cout << "#D HPrint() " << std::endl;
  TList* list = gDirectory->GetList();
#if 1
  Int_t count = 0;
  TIter itr(list);
  while(itr.Next() && ++count){
    const TString& name((*itr)->GetName());
    hddaq::cout << " " << std::setw(8) << std::left << name;
    if(count%10==0) hddaq::cout << std::endl;
  }
  hddaq::cout << std::endl;
#endif
  hddaq::cout << " NObject : " << list->GetEntries() << std::endl;
}

}

#endif

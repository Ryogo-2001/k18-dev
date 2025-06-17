// -*- C++ -*-

#ifndef DATABASE_PDG_HH
#define DATABASE_PDG_HH

#include <TPDGCode.h>
#include <TString.h>

namespace pdg
{
inline const TString& ClassName()
{
  static TString s_name("DatabasePDG");
  return s_name;
}

// Mass [GeV/c2]
Double_t Mass(Int_t pdg_code);
Double_t ElectronMass();
Double_t KaonMass();
Double_t PionMass();
Double_t ProtonMass();
Double_t NeutronMass();
Double_t LambdaMass();
Double_t SigmaNMass();
Double_t SigmaPMass();
Double_t XiMinusMass();
void     Print(Int_t pdg_code);
void     Print();
}

#endif

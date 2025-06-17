// -*- C++ -*-

// PDG code is defined in ROOT/include/TPDGCode.h
// Mass unit [GeV/c2]

#include "DatabasePDG.hh"

#include <string>

#include <TDatabasePDG.h>
#include <TMath.h>
#include <TParticlePDG.h>

#include <iostream>

namespace pdg
{
//_____________________________________________________________________________
Double_t
Mass(Int_t pdg_code)
{
  auto particle = TDatabasePDG::Instance()->GetParticle(pdg_code);
  return (particle ? particle->Mass() : TMath::QuietNaN());
}

//_____________________________________________________________________________
Double_t
ElectronMass()
{
  return Mass(kElectron);
}

//_____________________________________________________________________________
Double_t
KaonMass()
{
  return Mass(kKMinus);
}

//_____________________________________________________________________________
Double_t
PionMass()
{
  return Mass(kPiMinus);
}

//_____________________________________________________________________________
Double_t
ProtonMass()
{
  return Mass(kProton);
}

//_____________________________________________________________________________
Double_t
NeutronMass()
{
  return Mass(kNeutron);
}

//_____________________________________________________________________________
Double_t
LambdaMass()
{
  return Mass(kLambda0);
}

//_____________________________________________________________________________
Double_t
SigmaNMass()
{
  return Mass(kSigmaMinus);
}

//_____________________________________________________________________________
Double_t
SigmaPMass()
{
  return Mass(kSigmaPlus);
}

//_____________________________________________________________________________
Double_t
XiMinusMass()
{
  return Mass(kXiMinus);
}

//_____________________________________________________________________________
void
Print(Int_t pdg_code)
{
  auto particle = TDatabasePDG::Instance()->GetParticle(pdg_code);
  if(particle) particle->Print();
}

//_____________________________________________________________________________
void
Print()
{
  TDatabasePDG::Instance()->Print();
}

}

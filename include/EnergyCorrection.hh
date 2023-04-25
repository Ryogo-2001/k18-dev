#ifndef EnergyCorrection_hh 
#define EnergyCorrection_hh  1

#define PION   1
#define KAON   2
#define PROTON 3

#include "ThreeVector.hh"

// values about LH2 target geometry
const double CyLH2TgtR = 40.; //mm //diameter
const double CyLH2TgtZ = 300.; // mm
const double TargetVessThickness = 0.25; // mm

const double VaccumChamWinR = 80.;
const double VaccumChamThickness = 1.;
const double VaccumChamWinZ = 500.;
//

int CorrElossOut(double *mom_new, double *energy_new, double mom, int particle, ThreeVector dir, ThreeVector vertex);
int CorrElossOutWithCFRP(double *mom_new, double *energy_new, double mom, int particle, ThreeVector dir, ThreeVector vertex, ThreeVector pos0);
double diffE(int particle, double E, double length1, double length2, double Elast);
double diffE2(int particle, double E, double length1, double length2, double length3, double Elast);
int caldE(double momentum, double mass, double distance, double *momentum_cor, double *energy_cor, int tag);
double calc_dE_dx(double beta);
double calc_dE_dx2(double beta);
double calc_dE_dx3(double beta);
double mygamma(double beta);
double mybeta(double energy,double mormentum);
double calcLengthInTarget(ThreeVector dir, ThreeVector vertex);
double calcLengthInWindow(ThreeVector dir, ThreeVector vertex);
double calcLengthInCFRP(ThreeVector dir, ThreeVector pos0);

#endif 

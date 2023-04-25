#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ThreeVector.hh"
#include "EnergyCorrection.hh"

const int LH2Tag=0;
const int WindowTag=1;
const int CFRPTag=1;

/*
* Correct Energy loss (Downstream part of Vertex)
*/
int CorrElossOut(double *mom_new, double *energy_new, double mom, int particle, ThreeVector dir, ThreeVector vertex)
{
  double FL,FH,FTMP;
  double Elow, Ehigh, Elast, EPS;
  double E=0., mass, length=0.;
  double length2=0.;
  double TgtZ = CyLH2TgtZ/2.;

  switch (particle) {
  case PION:
    mass = 0.1395701;
    break;
  case KAON:
    mass = 0.493677;
    break;
  case PROTON:
    mass = 0.93827200;
    break;
  default:
    return -1;
    break;
  }

  if (fabs(vertex.z()) <= TgtZ) {
    length = calcLengthInTarget(dir, vertex);
    length2 = calcLengthInWindow(dir, vertex);
    //length2 = 0.;
  } else if (vertex.z() < (-1)*TgtZ) { /* Upstream of SciFi */
    length = calcLengthInTarget(dir, vertex);
    length2 = calcLengthInWindow(dir, vertex);
    //length2 = 0.;
  } else if (vertex.z() > TgtZ) {     /* Downstream of SciFi */
    *mom_new = mom;
    *energy_new = sqrt(mass*mass+mom*mom);
    return 1;
  }

  //printf("Length = %f, Length2 = %f\n", length, length2);

  Elow  = mass;
  Ehigh = 3.;
  Elast = sqrt(mass*mass+mom*mom);
  EPS = 0.001;
  FL  = diffE(particle, Elow, length, length2, Elast);
  FH  = diffE(particle, Ehigh, length, length2, Elast);
  while (fabs((Ehigh-Elow)/Ehigh) > EPS) {
    E = Ehigh-(Ehigh-Elow)*FH/(FH-FL);
    //printf("-------E=%f, Elow=%f, Ehigh=%f, FL=%f, FH=%f----------\n",E, Elow, Ehigh, FL, FH);
    //getchar();
     if (fabs(FTMP=diffE(particle, E, length, length2, Elast)) < 0.0000001){
      Elow = E;
      Ehigh = E;
      FH = FTMP;
    } if ((FTMP=diffE(particle, E, length, length2, Elast)) < 0) {
      Elow = E;
      FL = FTMP;
    } else if ((FTMP=diffE(particle, E, length, length2, Elast)) > 0){
      Ehigh = E;
      FH = FTMP;
    }
     //printf("******E=%f, Elow=%f, Ehigh=%f, FL=%f, FH=%f*******\n",E, Elow, Ehigh, FL, FH);
    //getchar();
  }
  *energy_new = E;
  *mom_new = sqrt(E*E-mass*mass);

  return 1;
}

int CorrElossOutWithCFRP(double *mom_new, double *energy_new, double mom, int particle, ThreeVector dir, ThreeVector vertex, ThreeVector pos0)
/*
  dir : CFT tracking direction
  vertex : Sp scattering vertex for energy loss correction in LH2
  pos0 : CFT pos0 for CFRP correction
*/

{
  double FL,FH,FTMP;
  double Elow, Ehigh, Elast, EPS;
  double E=0., mass, length=0.;
  double length2=0., length3=0.;
  double TgtZ = CyLH2TgtZ/2.;

  switch (particle) {
  case PION:
    mass = 0.1395701;
    break;
  case KAON:
    mass = 0.493677;
    break;
  case PROTON:
    mass = 0.93827200;
    break;
  default:
    return -1;
    break;
  }

  if (fabs(vertex.z()) <= TgtZ) {
    length = calcLengthInTarget(dir, vertex);
    length2 = calcLengthInWindow(dir, vertex);
    //length2 = 0.;
    length3 = calcLengthInCFRP(dir, pos0);
  } else if (vertex.z() < (-1)*TgtZ) { /* Upstream of SciFi */
    length = calcLengthInTarget(dir, vertex);
    length2 = calcLengthInWindow(dir, vertex);
    //length2 = 0.;
    length3 = calcLengthInCFRP(dir, pos0);
  } else if (vertex.z() > TgtZ) {     /* Downstream of SciFi */
    *mom_new = mom;
    *energy_new = sqrt(mass*mass+mom*mom);
    return 1;
  }

  //printf("Length = %f, Length2 = %f\n", length, length2);

  Elow  = mass;
  Ehigh = 3.;
  Elast = sqrt(mass*mass+mom*mom);
  EPS = 0.001;
  FL  = diffE2(particle, Elow, length, length2, length3, Elast);
  FH  = diffE2(particle, Ehigh, length, length2, length3, Elast);
  while (fabs((Ehigh-Elow)/Ehigh) > EPS) {
    E = Ehigh-(Ehigh-Elow)*FH/(FH-FL);
    //printf("-------E=%f, Elow=%f, Ehigh=%f, FL=%f, FH=%f----------\n",E, Elow, Ehigh, FL, FH);
    //getchar();
    if (fabs(FTMP=diffE2(particle, E, length, length2, length3, Elast)) < 0.0000001){
      Elow = E;
      Ehigh = E;
      FH = FTMP;
    } if ((FTMP=diffE2(particle, E, length, length2, length3, Elast)) < 0) {
      Elow = E;
      FL = FTMP;
    } else if ((FTMP=diffE2(particle, E, length, length2, length3, Elast)) > 0){
      Ehigh = E;
      FH = FTMP;
    }
     //printf("******E=%f, Elow=%f, Ehigh=%f, FL=%f, FH=%f*******\n",E, Elow, Ehigh, FL, FH);
    //getchar();
  }
  *energy_new = E;
  *mom_new = sqrt(E*E-mass*mass);

  return 1;
}



double diffE(int particle, double E, double length1, double length2, double Elast)
{
  double p, mass;
  double mom_new, energy_new;

  switch (particle) {
  case PION:
    mass = 0.1395701;
    break;
  case KAON:
    mass = 0.493677;
    break;
  case PROTON:
    mass = 0.93827200;
    break;
  default:
    return -1;
    break;
  }
  p = sqrt(E*E-mass*mass);

  caldE(p, mass, length1, &mom_new, &energy_new, LH2Tag);

  p = mom_new;
  caldE(p, mass, length2, &mom_new, &energy_new, WindowTag);


  return (energy_new-Elast);
  //return (sqrt(mass*mass+p*p)-caldE(particle,p,length)-Elast);
}

double diffE2(int particle, double E, double length1, double length2, double length3, double Elast)
{
  double p, mass;
  double mom_new, energy_new;

  switch (particle) {
  case PION:
    mass = 0.1395701;
    break;
  case KAON:
    mass = 0.493677;
    break;
  case PROTON:
    mass = 0.93827200;
    break;
  default:
    return -1;
    break;
  }
  p = sqrt(E*E-mass*mass);

  caldE(p, mass, length1, &mom_new, &energy_new, LH2Tag);

  p = mom_new;
  caldE(p, mass, length2, &mom_new, &energy_new, WindowTag);

  p = mom_new;
  caldE(p, mass, length3, &mom_new, &energy_new, CFRPTag);


  return (energy_new-Elast);
  //return (sqrt(mass*mass+p*p)-caldE(particle,p,length)-Elast);
}

/*
static double
caldE(int particle, double momentum, double length)
{
  switch (particle) {
  case PION:
    return 0.0002175560*length;
    break;
  case KON:
    return 0.00020610672*length;
    break;
  case PROTON:
    return 0;
    break;
  default:
    return -10;
    break;
  }
}
*/

int caldE(double momentum, double mass, double distance, double *momentum_cor, double *energy_cor, int tag)
{
  double dE_dx; /*MeV/cm*/
  double eloss; /*MeV*/

  double beta;
  double thickness = distance/10.0; /*cm*/
  double m = mass*1000.0;  /*mass of incident particle(Mev)*/

  double E_0;
  double E;
  double p = momentum*1000.0;
  //double delta=0.01; /*cm*/
  double delta=0.05; /*cm*/
  int i;
  double length=0.0;
  double total_eloss=0.0;

  //E_0=mygamma(beta)*m;
  //p=mygamma(beta)*beta*m;

  E_0= sqrt(m*m + p*p);
  E=E_0;
  beta = mybeta(E,p);
  //printf("beta=%e, E=%f, p=%f\n",beta, E, p);
  if (beta<=0.000001) {
    *momentum_cor = p/1000.0;
    *energy_cor = E/1000.0;
    return 1;
  }
  dE_dx=0.0;
  eloss=0.0;
  for(i=0;i<=thickness/delta;i++){
    if (tag == LH2Tag)
      dE_dx=calc_dE_dx(beta);
    else if (tag == WindowTag)
      dE_dx=calc_dE_dx2(beta);
    else if (tag == CFRPTag)
      dE_dx=calc_dE_dx3(beta);
    else {
      fprintf(stderr,"caldE: Invalid tag %d\n",tag);
      exit(-1);
    }
    eloss=dE_dx*delta;
    E=E-eloss;
    if(E<m){
      //fprintf(stderr,"particle stops in material at %5.3fcm\n",length);
      *momentum_cor = p/1000.0;
      *energy_cor = E/1000.0;
      return 0;
      //break;
    }
    p=sqrt(pow(E,2.0)-pow(m,2.0));
    beta=mybeta(E,p);
    length=length+delta;
    total_eloss=total_eloss+eloss;
    /*
    printf("beta:%5.3f\n",beta);
    printf("dE_dx:%5.3f\teloss:%5.3f\n",dE_dx,eloss);
    printf("E:%5.3f(MeV)\tp:%5.3f(MeV/c)\n",E,p);
    printf("length:%5.3f(cm)\n",length);
    printf("total energy loss:%5.3f(MeV)\n",total_eloss);
    */
    //getchar();
  }
  *momentum_cor = p/1000.0;
  *energy_cor = E/1000.0;

  return 1;
}


double calc_dE_dx(double beta)
{
  double value;
  const double C=0.1535; /*MeVcm^2/g*/
  const double m_e=0.511;
  double logterm;
  double gamma_2;
  double W_max;

  double rho=0.0708;   /*g/cm^3 (C)*/
  double I=21.8;     /*eV*/
  double Z_A=0.99216;
  int z=1;
  //extern double rho;   /*g/cm^3*/
  //extern double I;     /*eV*/
  //extern double Z_A;
  //extern int z;

  gamma_2=1/(1-pow(beta,2.0));
  //printf("sqrt(gamma_2):%f\n",sqrt(gamma_2));

  W_max=2.0*m_e*pow(beta,2.0)*gamma_2;

  logterm=log(2.0*m_e*gamma_2*pow(beta,2.0)*W_max*pow(10.0,12.0)/pow(I,2.0))-2.0*pow(beta,2.0);

  value=C*rho*Z_A*pow((double)z,2.0)*logterm/pow(beta,2.0);

  return value;

}

double calc_dE_dx2(double beta)
{
  double value;
  const double C=0.1535; /*MeVcm^2/g*/
  const double m_e=0.511;
  double logterm;
  double gamma_2;
  double W_max;

  double rho=1.39;   /*g/cm^3 (C)*/
  double I=64.7;     /*eV*/ // temporary
  double Z_A=0.54155;       // temporary
  int z=1;
  //extern double rho;   /*g/cm^3*/
  //extern double I;     /*eV*/
  //extern double Z_A;
  //extern int z;

  gamma_2=1/(1-pow(beta,2.0));
  //printf("sqrt(gamma_2):%f\n",sqrt(gamma_2));

  W_max=2.0*m_e*pow(beta,2.0)*gamma_2;

  logterm=log(2.0*m_e*gamma_2*pow(beta,2.0)*W_max*pow(10.0,12.0)/pow(I,2.0))-2.0*pow(beta,2.0);

  value=C*rho*Z_A*pow((double)z,2.0)*logterm/pow(beta,2.0);

  return value;

}

double calc_dE_dx3(double beta)
{
  double value;
  const double C=0.1535; /*MeVcm^2/g*/
  const double m_e=0.511;
  double logterm;
  double gamma_2;
  double W_max;

  double rho=1.53;   /*g/cm^3 (C)*/
  double I=78.0;     /*eV*/
  double Z_A=0.49954;
  int z=1;

  gamma_2=1/(1-pow(beta,2.0));
  //printf("sqrt(gamma_2):%f\n",sqrt(gamma_2));

  W_max=2.0*m_e*pow(beta,2.0)*gamma_2;

  logterm=log(2.0*m_e*gamma_2*pow(beta,2.0)*W_max*pow(10.0,12.0)/pow(I,2.0))-2.0*pow(beta,2.0);

  value=C*rho*Z_A*pow((double)z,2.0)*logterm/pow(beta,2.0);

  return value;

}

double mygamma(double beta)
{
  double value;

  value=1.0/sqrt(1.0-pow(beta,2.0));

  return value;

}

double mybeta(double energy,double mormentum)
{
  double value;

  value=mormentum/energy;

  return value;
}


double calcLengthInTarget(ThreeVector dir, ThreeVector vertex)
{
  double t1, t2;

  double TgtR = CyLH2TgtR/2.;
  double TgtZ = CyLH2TgtZ/2.;

  if (dir.z()>0.)
    t1= (TgtZ - vertex.z())/dir.z();
  else if (dir.z()<0.)
    t1= (-TgtZ - vertex.z())/dir.z();
  else // dir.z() == 0.
    t1 = 100000000000000000.;

  double ux = dir.x(), uy = dir.y();
  double vx = vertex.x(), vy = vertex.y();

  double a, b, c;
  a = ux*ux+uy*uy;
  b = ux*vx+uy*vy;
  c = vx*vx+vy*vy-TgtR*TgtR;

  if (b*b-a*c>=0)
    t2 = (-b+sqrt(b*b-a*c))/a;
  else {
    //fprintf(stderr, "calcLengthInTarget: Hanbetsu-shiki is negative %f\n", b*b-a*c);
    //t2 = 10000.;
    return 0.;
  }

  //printf("t1 = %f, t2 = %f\n", t1, t2);
  ThreeVector length_vec;
  if (t1 <= t2)
    length_vec = dir*t1;
  else
    length_vec = dir*t2;

  return length_vec.Mag();

}

double calcLengthInWindow(ThreeVector dir, ThreeVector vertex)
{
  double t1, t2, t3;

  double TgtR = CyLH2TgtR/2.;
  double TgtZ = CyLH2TgtZ/2.;

  if (dir.z()>0.)
    t1= (TgtZ - vertex.z())/dir.z();
  else if (dir.z()<0.)
    t1= (-TgtZ - vertex.z())/dir.z();
  else // dir.z() == 0.
    t1 = 100000000000000000.;

  double ux = dir.x(), uy = dir.y();
  double vx = vertex.x(), vy = vertex.y();

  double a, b, c;
  a = ux*ux+uy*uy;
  b = ux*vx+uy*vy;
  c = vx*vx+vy*vy-TgtR*TgtR;

  if (b*b-a*c>=0)
    t2 = (-b+sqrt(b*b-a*c))/a;
  else {
    //fprintf(stderr, "calcLengthInWindow: Hanbetsu-shiki is negative %f\n", b*b-a*c);
    //t2 = 10000.;
    return 0.;
  }

  double VesselR = CyLH2TgtR/2.+ TargetVessThickness;
  [[maybe_unused]] double VesselZ = CyLH2TgtZ/2.;

  a = ux*ux+uy*uy;
  b = ux*vx+uy*vy;
  c = vx*vx+vy*vy-VesselR*VesselR;

  if (b*b-a*c>=0)
    t3 = (-b+sqrt(b*b-a*c))/a;
  else {
    //fprintf(stderr, "calcLengthInWindow: Hanbetsu-shiki is negative %f\n", b*b-a*c);
    //t3 = 10000.;
    return 0.;
  }

  double t_vessel;
  if (t1 <= t2)
    t_vessel = 0;
  else
    t_vessel = t3-t2;

  /*
  double t4, t5;
  double VaccumWinR = VaccumChamWinR/2.;
  double VaccumWinZ = VaccumChamWinZ/2.;

  a = ux*ux+uy*uy;
  b = ux*vx+uy*vy;
  c = vx*vx+vy*vy-VaccumWinR*VaccumWinR;

  if (b*b-a*c>=0)
    t4 = (-b+sqrt(b*b-a*c))/a;
  else {
    //fprintf(stderr, "calcLengthInWindow: Hanbetsu-shiki is negative %f\n", b*b-a*c);
    //t4 = 10000.;
    return 0.;
  }

  double VaccumWinROut = VaccumChamWinR/2.+VaccumChamThickness;
  a = ux*ux+uy*uy;
  b = ux*vx+uy*vy;
  c = vx*vx+vy*vy-VaccumWinROut*VaccumWinROut;

  if (b*b-a*c>=0)
    t5 = (-b+sqrt(b*b-a*c))/a;
  else {
    //fprintf(stderr, "calcLengthInWindow: Hanbetsu-shiki is negative %f\n", b*b-a*c);
    //t5 = 10000.;
    return 0.;
  }

  double t_window;
  t_window = t5-t4;
  */

  ThreeVector length_vec;
  //length_vec = dir*(t_vessel+t_window);
  length_vec = dir*(t_vessel);

  return length_vec.Mag();

}

double calcLengthInCFRP(ThreeVector dir, ThreeVector pos0)
{
  [[maybe_unused]] double t1, t2, t3;

  double VaccumWinR = VaccumChamWinR/2.;
  double VaccumWinZ = VaccumChamWinZ/2.;

  if (dir.z()>0.)
    t1= (VaccumWinZ - pos0.z())/dir.z();
  else if (dir.z()<0.)
    t1= (-VaccumWinZ - pos0.z())/dir.z();
  else // dir.z() == 0.
    t1 = 100000000000000000.;

  double x1 = dir.x()*t1+pos0.x();
  double y1 = dir.y()*t1+pos0.y();
  double r1 = sqrt(x1*x1+y1*y1);
  if (r1<VaccumWinR) {
    fprintf(stderr, "CFRP does not hit\n");
    return 0.;
  }

  double ux = dir.x(), uy = dir.y();
  double vx = pos0.x(), vy = pos0.y();



  double t4, t5;

  double a, b, c;
  a = ux*ux+uy*uy;
  b = ux*vx+uy*vy;
  c = vx*vx+vy*vy-VaccumWinR*VaccumWinR;

  if (b*b-a*c>=0) {
    double t4_1 = (-b+sqrt(b*b-a*c))/a;
    double t4_2 = (-b-sqrt(b*b-a*c))/a;

    if (fabs(t4_1)<fabs(t4_2))
      t4 = t4_1;
    else
      t4 = t4_2;
  } else {
    //fprintf(stderr, "calcLengthInWindow: Hanbetsu-shiki is negative %f\n", b*b-a*c);
    //t4 = 10000.;
    return 0.;
  }

  double VaccumWinROut = VaccumChamWinR/2.+VaccumChamThickness;
  a = ux*ux+uy*uy;
  b = ux*vx+uy*vy;
  c = vx*vx+vy*vy-VaccumWinROut*VaccumWinROut;

  if (b*b-a*c>=0) {
    double t5_1 = (-b+sqrt(b*b-a*c))/a;
    double t5_2 = (-b-sqrt(b*b-a*c))/a;

    if (fabs(t5_1)<fabs(t5_2))
      t5 = t5_1;
     else
      t5 = t5_2;
  } else {
    //fprintf(stderr, "calcLengthInWindow: Hanbetsu-shiki is negative %f\n", b*b-a*c);
    //t5 = 10000.;
    return 0.;
  }

  double t_window;
  t_window = t5-t4;

  ThreeVector length_vec;
  length_vec = dir*(t_window);

  return length_vec.Mag();

}

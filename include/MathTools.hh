// -*- C++ -*-

#ifndef MATH_TOOLS_HH
#define MATH_TOOLS_HH

#include <TMath.h>
#include <TString.h>
#include <TVector3.h>

//_____________________________________________________________________________
namespace MathTools
{
const TString& ClassName();

//_____________________________________________________________________________
inline double Infinity() { return TMath::Infinity(); }
inline double Epsilon() { return TMath::Limits<double>::Epsilon(); }
inline int    Round(double a) { return int(a+0.5-(a<0.)); }

//_____________________________________________________________________________
template <typename T>
inline T Sign(T a) { return (a>0) - (a<0); }

//_____________________________________________________________________________
template <typename Container>
inline double Accumulate(const Container& c)
{
  double sum = 0.;
  typename Container::const_iterator itr, end=c.end();
  for(itr=c.begin(); itr!=end; ++itr){
    sum += *itr;
  }
  return sum;
}

//_____________________________________________________________________________
template <typename Container>
inline double Average(const Container& c)
{
  return Accumulate(c)/c.size();
}

//_____________________________________________________________________________
template <typename Container>
inline double Deviation(const Container& c)
{
  double sum  = 0.;
  double sum2 = 0.;
  typename Container::const_iterator itr, end=c.end();
  for(itr=c.begin(); itr!=end; ++itr){
    sum  += *itr;
    sum2 += (*itr) * (*itr);
  }
  std::size_t n = c.size();
  return std::sqrt((sum2-(sum*sum/n)) / n);
}

//_____________________________________________________________________________
template <typename Container>
inline double MaxElement(const Container& c)
{
  return *std::max_element(c.begin(), c.end());
}

//_____________________________________________________________________________
template <typename Container>
inline double MinElement(const Container& c)
{
  return *std::min_element(c.begin(), c.end());
}

//_____________________________________________________________________________
// compare l and r
template <typename T>
inline Bool_t Equal(const T& l, const T& r);

// 3.14159265358979323846
inline double Pi() { return TMath::Pi(); }
// velocity of light : 299.792458 mm ns^-1
inline double C() { return TMath::C()*1e-6; }
// base of natural log
inline double E() { return TMath::E(); }

Bool_t SolveGaussJordan(const std::vector<double>& z,
                        const std::vector<double>& w,
                        const std::vector<double>& s,
                        const std::vector<double>& ct,
                        const std::vector<double>& st,
                        double& x0,
                        double& u0,
                        double& y0,
                        double& v0);

Bool_t GaussElim(double **a, int n, double *b, int *indx, int *ipiv);
Bool_t GaussJordan(double **a, int n, double *b,
                   int *indxc, int *indxd, int *ipiv);
Bool_t InterpolateRatio(int n, const double *xa, const double *ya,
                        double *w1, double *w2,
                        double x, double &y, double &dy);
Bool_t InterpolatePol(int n, const double *xa, const double *ya,
                      double *w1, double *w2,
                      double x, double &y, double &dy);
Bool_t SVDksb(double **u, const double *w, double **v,
              int m, int n, const double *b, double *x, double *wv);
Bool_t SVDcmp(double **a, int m, int n, double *w, double **v, double *wv);
template <typename T>
void PrintMatrix(T *mat, const std::string& arg="",
                 const int column=4, const int line=4);
template <typename T>
void PrintVector(T *vec, const std::string& arg="",
                 const std::size_t size=4);
}

//_____________________________________________________________________________
inline const TString&
MathTools::ClassName()
{
  static TString s_name("MathTools");
  return s_name;
}

//_____________________________________________________________________________
template <typename T>
inline Bool_t
MathTools::Equal(const T& l, const T& r)
{
  double al = TMath::Abs(l), ar = TMath::Abs(r);
  return
    ((al <= Epsilon()) || (ar <= Epsilon()))
    ?
    (TMath::Abs(l-r) <= Epsilon())
    :
    ((TMath::Abs(l-r)/std::max(al, ar)) <= Epsilon());
}


#endif

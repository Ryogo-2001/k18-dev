// -*- C++ -*-

#ifndef K18_TRANS_MATRIX_HH
#define K18_TRANS_MATRIX_HH

#include <TString.h>

//_____________________________________________________________________________
class K18TransMatrix
{
public:
  static const TString&  ClassName();
  static K18TransMatrix& GetInstance();
  ~K18TransMatrix();

private:
  K18TransMatrix();
  K18TransMatrix(const K18TransMatrix&);
  K18TransMatrix& operator=(const K18TransMatrix&);

private:
  enum NameX {
    X,   A,   T,   XX,  XA,  XT,  AA,  AT,  TT,  YY,
    YB,  BB,  XXX, XXA, XXT, XAA, XAT, XTT, XYY, XYB,
    XBB, AAA, AAT, ATT, AYY, AYB, ABB, TTT, TYY, TYB,
    TBB,
    size_NameX };
  enum NameY {
    Y, B, YX, YA, YT, BX, BA, BT,
    size_NameY };

  Bool_t   m_is_ready;
  TString  m_file_name;
  Double_t m_X[size_NameX];
  Double_t m_Y[size_NameY];
  Double_t m_U[size_NameX];
  Double_t m_V[size_NameY];

public:
  Bool_t CalcDeltaD2U(Double_t xin, Double_t yin, Double_t uin, Double_t vin,
                      Double_t xout, Double_t& yout,
                      Double_t& uout, Double_t& vout,
                      Double_t& delta1, Double_t& delta2) const;
  Bool_t Initialize();
  Bool_t Initialize(const TString& file_name);
  Bool_t IsReady() const { return m_is_ready; }
  void   SetFileName(const TString& file_name) { m_file_name = file_name; }
  Bool_t Transport(Double_t xin, Double_t yin, Double_t uin, Double_t vin,
                   Double_t delta,
                   Double_t& xout, Double_t& yout,
                   Double_t& uout, Double_t& vout) const;
};

//_____________________________________________________________________________
inline const TString&
K18TransMatrix::ClassName()
{
  static TString s_name("K18TransMatrix");
  return s_name;
}

//_____________________________________________________________________________
inline K18TransMatrix&
K18TransMatrix::GetInstance()
{
  static K18TransMatrix s_instance;
  return s_instance;
}

#endif

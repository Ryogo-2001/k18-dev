// -*- C++ -*-

#ifndef MATRIX_PARAM_MAN_HH
#define MATRIX_PARAM_MAN_HH

#include <vector>
#include <map>

#include <TString.h>

//_____________________________________________________________________________
class MatrixParamMan
{
public:
  static const TString&  ClassName();
  static MatrixParamMan& GetInstance();
  ~MatrixParamMan();

private:
  MatrixParamMan();
  MatrixParamMan(const MatrixParamMan&);
  MatrixParamMan& operator =(const MatrixParamMan&);

private:
  typedef std::vector<std::vector<Double_t>> Matrix2D;
  typedef std::vector<std::vector< std::vector<Double_t>>> Matrix3D;
  Bool_t   m_is_ready;
  TString  m_file_name_2d;
  TString  m_file_name_3d;
  Matrix2D m_enable_2d;
  Matrix3D m_enable_3d;

public:
  Bool_t Initialize();
  Bool_t Initialize(const TString& filename_2d,
                    const TString& filename_3d);
  Bool_t IsAccept(Int_t detA, Int_t detB) const;
  Bool_t IsAccept(Int_t detA, Int_t detB, Int_t detC) const;
  Bool_t IsReady() const { return m_is_ready; }
  void   Print2D(const TString& arg="") const;
  void   Print3D(const TString& arg="") const;
  void   SetMatrix2D(const TString& file_name);
  void   SetMatrix3D(const TString& file_name);

};

//_____________________________________________________________________________
inline MatrixParamMan&
MatrixParamMan::GetInstance()
{
  static MatrixParamMan g_instance;
  return g_instance;
}

//_____________________________________________________________________________
inline const TString&
MatrixParamMan::ClassName()
{
  static TString g_name("MatixParamMan");
  return g_name;
}

#endif

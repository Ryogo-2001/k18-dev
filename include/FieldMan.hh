// -*- C++ -*-

#ifndef FIELD_MAN_HH
#define FIELD_MAN_HH

#include <vector>

#include <TString.h>
#include <TVector3.h>

class HypsFieldMap;
class FieldElements;

typedef std::vector<FieldElements*> FEContainer;
typedef FEContainer::const_iterator FEIterator;

//_____________________________________________________________________________
class FieldMan
{
public:
  static const TString& ClassName();
  static FieldMan&      GetInstance();
  ~FieldMan();

private:
  FieldMan();
  FieldMan(const FieldMan&);
  FieldMan & operator =(const FieldMan&);

private:
  Bool_t          m_is_ready;
  TString         m_file_name_hyps;
  TString         m_file_name_shs;
  HypsFieldMap* m_hyps_map;
  HypsFieldMap* m_shs_map;
  FEContainer     m_element_list;

public:
  Bool_t   Initialize();
  Bool_t   Initialize(const TString& file_name_hyps);
  Bool_t   Initialize(const TString& file_name_hyps,
                      const TString& file_name_shs);
  Bool_t   IsReady() const { return m_is_ready; }
  TVector3 GetField(const TVector3& position) const;
  TVector3 GetdBdX(const TVector3& position) const;
  TVector3 GetdBdY(const TVector3& position) const;
  TVector3 GetdBdZ(const TVector3& position) const;
  void     ClearElementsList();
  void     AddElement(FieldElements* element);
  void     SetHypsFileName(const TString& file_name) { m_file_name_hyps = file_name; }
  void     SetShsFileName(const TString& file_name) { m_file_name_shs = file_name; }
  Double_t StepSize(const TVector3& position,
                    Double_t default_step_size,
                    Double_t min_step_size) const;
};

//_____________________________________________________________________________
inline const TString&
FieldMan::ClassName()
{
  static TString s_name("FieldMan");
  return s_name;
}

//_____________________________________________________________________________
inline FieldMan&
FieldMan::GetInstance()
{
  static FieldMan s_instance;
  return s_instance;
}

#endif

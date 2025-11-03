// -*- C++ -*-

#ifndef USER_PARAM_MAN_HH
#define USER_PARAM_MAN_HH

#include <map>
#include <vector>
#include <TString.h>

class TNamed;

//_____________________________________________________________________________
class UserParamMan
{
public:
  static const TString& ClassName();
  static UserParamMan&  GetInstance();
  ~UserParamMan();

private:
  UserParamMan();
  UserParamMan(const UserParamMan& );
  UserParamMan& operator =(const UserParamMan&);

private:
  typedef std::vector<Double_t>         ParamArray;
  typedef std::map<TString, ParamArray> ParamMap;
  typedef ParamMap::const_iterator      PIterator;
  Bool_t   m_is_ready;
  Bool_t   m_use_default;
  TString  m_file_name;
  ParamMap m_param_map;
  TString  m_buf;
  TNamed*  m_object;

public:
  void     AddObject();
  Bool_t   Initialize();
  Bool_t   Initialize(const TString& filename);
  Bool_t   IsReady() const { return m_is_ready; }
  Int_t    GetSize(const TString& key) const;
  Double_t GetParameter(const TString& key, Int_t i=0) const;
  Double_t Get(const TString& key, Int_t i=0) const
    { return GetParameter(key, i); }
  void     Print(const TString& arg="") const;
  void     UseDefaultValue(Bool_t flag=true) { m_use_default = flag; }
};

//_____________________________________________________________________________
inline const TString&
UserParamMan::ClassName()
{
  static TString s_name("UserParamMan");
  return s_name;
}

//_____________________________________________________________________________
inline UserParamMan&
UserParamMan::GetInstance()
{
  static UserParamMan s_instance;
  return s_instance;
}

#endif

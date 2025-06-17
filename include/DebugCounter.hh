// -*- C++ -*-

#ifndef DEBUG_COUNTER_HH
#define DEBUG_COUNTER_HH

#include <map>

#include <TString.h>

#include <std_ostream.hh>

//_____________________________________________________________________________
namespace debug
{
class ObjectCounter
{
public:
  static const TString& ClassName();
  static ObjectCounter& GetInstance();
  ~ObjectCounter();

private:
  ObjectCounter();
  ObjectCounter(const ObjectCounter&);
  ObjectCounter& operator =(const ObjectCounter&);

private:
  typedef std::map<TString, Int_t>  ObjectMap;
  typedef ObjectMap::const_iterator ObjectIter;
  ObjectMap m_map;

public:
  void check(const TString& arg="") const;
  void print(const TString& arg="") const;

public:
  static void decrease(const TString& key);
  static void increase(const TString& key);
};

//_____________________________________________________________________________
inline const TString&
ObjectCounter::ClassName()
{
  static TString s_name("ObjectCounter");
  return s_name;
}

//_____________________________________________________________________________
inline ObjectCounter&
ObjectCounter::GetInstance()
{
  static ObjectCounter s_instance;
  return s_instance;
}

//_____________________________________________________________________________
inline void
ObjectCounter::decrease(const TString& key)
{
#ifdef MemoryLeak
  --(GetInstance().m_map[key]);
#endif
}

//_____________________________________________________________________________
inline void
ObjectCounter::increase(const TString& key)
{
#ifdef MemoryLeak
  ++(GetInstance().m_map[key]);
#endif
}

}

#endif

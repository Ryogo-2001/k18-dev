// -*- C++ -*-

#ifndef DEBUG_TIMER_HH
#define DEBUG_TIMER_HH

#include <ostream>
#include <string>

#include <std_ostream.hh>

namespace debug
{
//_____________________________________________________________________________
class Timer
{
public:
  explicit Timer(const std::string& msg="",
                 bool verbose=true);
  ~Timer();

private:
  Timer(const Timer&);
  Timer& operator=(const Timer&);

private:
  ::timespec* m_start;
  ::timespec* m_stop;
  std::string m_cstart;
  std::string m_cstop;
  std::string m_msg;
  bool        m_verbose;

public:
  double sec() const;
  double nsec() const;
  void   stop();
  void   print(const std::string& arg="",
               std::ostream& ost=hddaq::cout) const;
  void   print_verbose(const std::string& arg="",
                       std::ostream& ost=hddaq::cout) const;
};
}
#endif

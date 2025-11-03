#ifndef CFTPedCorMan_h

#define CFTPedCorMan_h 1

#include "DetectorID.hh"
#include <vector>

class RawData;

struct CFTPedCorParam {
  int layer;
  int segment;
  double sigma;
  double thr;
  int cor_layer;
  int cor_segment;
  double p0_HG;
  double p1_HG;
  double p0_LG;
  double p1_LG;
 
  bool operator<(const CFTPedCorParam& another) const {
    return sigma < another.sigma;
  }
};


class CFTPedCorMan
{
public:
  static CFTPedCorMan&        GetInstance( void );
  static const std::string& ClassName( void );
  ~CFTPedCorMan( void );

private:
  CFTPedCorMan( void );
  CFTPedCorMan( const CFTPedCorMan& );
  CFTPedCorMan& operator =( const CFTPedCorMan& );

private:
  std::vector <CFTPedCorParam> m_container[NumOfPlaneCFT][NumOfSegCFT_PHI4];

  bool          m_is_ready;
  std::string   m_file_name;

public:
  bool Initialize( void );
  bool Initialize( const std::string& file_name );
  bool IsReady( void ) const { return m_is_ready; }
  void SetFileName( const std::string& file_name ) { m_file_name = file_name; }

  bool PedestalCorrection(int layer, int segment, double &deltaHG, double &deltaLG, RawData* rawData ) const;
private:
};

//______________________________________________________________________________
inline CFTPedCorMan&
CFTPedCorMan::GetInstance( void )
{
  static CFTPedCorMan g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const std::string&
CFTPedCorMan::ClassName( void )
{
  static std::string g_name("CFTPedCorMan");
  return g_name;
}




#endif

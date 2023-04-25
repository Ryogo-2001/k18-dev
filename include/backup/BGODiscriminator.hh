#ifndef BGO_DISCRI_HH
#define BGO_DISCRI_HH

#include <string>
#include <vector>
#include <cstdlib>
#include <unistd.h>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#define GBUF 2048

class BGODiscriminator{
 private:
  int m_sample_size;
  //double GraphX[GBUF];
  //double GraphY[GBUF];
  std::vector <double> m_GraphX;
  std::vector <double> m_GraphY;

  std::vector<double> m_target;
  std::vector<int> m_cross;
  std::vector<int> m_FoundPeak;
  std::vector<int> m_RisingPoint;
  std::vector<int> m_FallingPoint;
 public:
  ~BGODiscriminator();
  BGODiscriminator();
  BGODiscriminator(std::vector<double> &v);

  void SetGraph(TGraphErrors *graph);

  bool SelectRange(double threshold,double begin,double last);
  bool SelectRange(double threshold);
  void Clear();  
  void AllClear();
  bool SloveBin(double input,double threshold,bool LorS);
  int GetPeakNum(double threshold,bool AbsFlag,bool LorS);
  int GetCrossPoint(double threshold);
  
  void GetRisingPoint(std::vector<double> &v);
  void GetFallingPoint(std::vector<double> &v);


};





#endif //INC_CDISCRI

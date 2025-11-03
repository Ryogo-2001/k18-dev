// -*- C++ -*-

#ifndef DC_PARAMETERS_HH
#define DC_PARAMETERS_HH

#include <TString.h>

#include <std_ostream.hh>

//_____________________________________________________________________________
struct DCPairPlaneInfo
{
  Bool_t pair;
  Bool_t honeycomb;
  Bool_t fiber;
  Int_t  id1, id2;
  Double_t CellSize;
  void Print(const TString& arg="", std::ostream& ost=hddaq::cout) const
  {
    ost << "[DCPairPlaneInfo::Print()] " << arg << std::endl
	<< " pair      : " << pair      << std::endl
	<< " honeycomb : " << honeycomb << std::endl
	<< " fiber     : " << fiber     << std::endl
	<< " layer1    : " << id1       << std::endl
	<< " layer2    : " << id2       << std::endl
	<< " cell size : " << CellSize  << std::endl;
  }
};

extern const DCPairPlaneInfo PPInfoBcOut[], PPInfoSdcIn[], PPInfoSdcOut[];
extern const Int_t NPPInfoBcOut, NPPInfoSdcIn, NPPInfoSdcOut;

#ifdef DefStatic
const DCPairPlaneInfo PPInfoBcOut[] = {
  // { pair_plane, honeycomb, fiber, id1, id2, CellSize }
  { true, false, false, 0,  1,  3.0 }, { true, false, false,  2,  3,  3.0 },
  { true, false, false, 4,  5,  3.0 }, { true, false, false,  6,  7,  3.0 },
  { true, false, false, 8, 9,  3.0 }, { true, false, false, 10, 11,  3.0 },
  // { false, false, false, 12, 12, 10.0 }
};

const DCPairPlaneInfo PPInfoSdcIn[] = {
  // { pair_plane, honeycomb, fiber, id1, id2, CellSize }
  { true, false, false, 0, 1,  5.0 }, { true, false, false, 2,  3,  5.0 },//SDC0
  { false, true, false, 4, 4,  12.0 },{ false, true, false, 5, 5, 12.0},//SDC1
  { true, true, false, 6, 7, 12.0 }, { true, true, false, 8, 9, 12.0 },//SDC1
};

const DCPairPlaneInfo PPInfoSdcOut[] = {
  // { pair_plane, honeycomb, fiber, id1, id2, CellSize }
  { false, true, false, 0, 0,  20.0 }, { true, true, false, 1, 2,  20.0 },{true, true, false, 3,4, 20.0}, //SDC2
  { false, true, false, 5, 5,  20.0 }, { true, true, false, 6, 7,  20.0 },{ true, true, false, 8, 9,  20.0 },//SDC3
};

const Int_t NPPInfoBcOut  = sizeof(PPInfoBcOut)/sizeof(DCPairPlaneInfo);
const Int_t NPPInfoSdcIn  = sizeof(PPInfoSdcIn)/sizeof(DCPairPlaneInfo);
const Int_t NPPInfoSdcOut = sizeof(PPInfoSdcOut)/sizeof(DCPairPlaneInfo);

#endif

//DL Ranges (BC1&2 for Time range -5 ns <[Time gate]<75 ns)
const Double_t MinDLBc[25] = {
   0.0,
   // BC1
  -5.0, -5.0, -5.0, -5.0, -5.0, -5.0,
   // BC2
  -5.0, -5.0, -5.0, -5.0, -5.0, -5.0,
   // BC3
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
   // BC4
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5
};

const Double_t MaxDLBc[25] = {
  0.0,
  // BC1
  75.0, 75.0, 75.0, 75.0, 75.0, 75.0,
  // BC2
  75.0, 75.0, 75.0, 75.0, 75.0, 75.0,
  // BC3
  1.8, 1.8, 1.8, 1.8, 1.8, 1.8,
  // BC4
  1.8, 1.8, 1.8, 1.8, 1.8, 1.8
};

const Double_t MinDLSdc[] = {
  0.0,
  //SDC0
  -0.5, -0.5, -0.5, -0.5,
  //SDC1
  -0.5, -0.5, -0.5, -0.5,-0.5, -0.5,
  // Dummy Id 11-30
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  // SDC2
  -0.5, -0.5, -0.5, -0.5, -0.5,
  // SDC3
  -0.5, -0.5, -0.5, -0.5, -0.5,
};

const Double_t MaxDLSdc[] = {
  0.0,
  // SDC0
  3.0, 3.0, 3.0, 3.0,
  // SDC1
  6.5, 6.5, 6.5, 6.5, 6.5, 6.5,
  // Dummy Id 11-30
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  // SDC2
  10.5, 10.5, 10.5, 10.5, 10.5, 
  // SDC3
  10.5, 10.5, 10.5, 10.5, 10.5, 
};

#endif

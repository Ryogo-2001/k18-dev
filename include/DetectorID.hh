// -*- C++ -*-

#ifndef DETECTOR_ID_HH
#define DETECTOR_ID_HH

#include <iostream>
#include <vector>

#include <TString.h>

// Counters ___________________________________________________________
const Int_t DetIdBH1      =  1;
const Int_t DetIdBH2      =  2;
const Int_t DetIdBAC      =  3;
const Int_t DetIdSCH      =  6;
const Int_t DetIdTOF      =  7;
const Int_t DetIdHTOF     =  8;
const Int_t DetIdBVH      =  9;
const Int_t DetIdLAC      = 12;
const Int_t DetIdWC       = 13;
const Int_t DetIdWCSUM    = 14; // Dummy
const Int_t NumOfSegBH1   = 11;
const Int_t NumOfSegBH2   =  8;
const Int_t NumOfSegBAC   =  2;
const Int_t NumOfSegSCH   = 64;
const Int_t NumOfSegTOF   = 24;
const Int_t NumOfSegHTOF  = 34;
const Int_t NumOfSegBVH   =  4;
const Int_t NumOfSegLAC   = 30;
const Int_t NumOfSegWC    = 20;

// Misc _______________________________________________________________
const Int_t DetIdTrig       = 21;
const Int_t DetIdScaler     = 22;
const Int_t DetIdMsT        = 25;
const Int_t DetIdMtx        = 26;
const Int_t DetIdVmeRm      = 81;
const Int_t DetIdMsTRM      = 82;
const Int_t DetIdHulRM      = 83;
const Int_t NumOfSegScaler  = 96;
const Int_t NumOfPlaneVmeRm = 2;

// Trigger Flag
namespace trigger
{
  enum ETriggerFlag
  {
    kL1SpillOn,
    kL1SpillOff,
    kSpillEnd,
    kSpillOnEnd,
    kCommonStopSdcOut,
    kMatrix2D1,
    kMatrix2D2,
    kMatrix3D,
    kBeamA,
    kBeamB,
    kBeamC,
    kBeamD,
    kBeamE,
    kBeamF,
    kTrigA,
    kTrigB,
    kTrigC,
    kTrigD,
    kTrigE,
    kTrigF,
    kTrigAPS,
    kTrigBPS,
    kTrigCPS,
    kTrigDPS,
    kTrigEPS,
    kTrigFPS,
    kLevel1A,
    kLevel1B,
    kClockPS,
    kReserve2PS,
    kLevel1OR,
    kEssDischarge,
    NTriggerFlag
  };

  const std::vector<TString> STriggerFlag =
    {
     "L1SpillOn",
     "L1SpillOff",
     "SpillEnd",
     "SpillOnEnd",
     "CommonStopSdcOut",
     "Matrix2D1",
     "Matrix2D2",
     "Matrix3D",
     "BeamA",
     "BeamB",
     "BeamC",
     "BeamD",
     "BeamE",
     "BeamF",
     "TrigA",
     "TrigB",
     "TrigC",
     "TrigD",
     "TrigE",
     "TrigF",
     "TrigA-PS",
     "TrigB-PS",
     "TrigC-PS",
     "TrigD-PS",
     "TrigE-PS",
     "TrigF-PS",
     "Level1A",
     "Level1B",
     "Clock-PS",
     "Reserve2-PS",
     "Level1OR",
     "EssDischarge",
    };
}
const Int_t NumOfSegTrig = trigger::NTriggerFlag;

const Int_t DetIdVmeCalib      = 999;
const Int_t NumOfPlaneVmeCalib =   5;
const Int_t NumOfSegVmeCalib   =  32;

// Trackers ___________________________________________________________
const Int_t DetIdBC3  = 103;
const Int_t DetIdBC4  = 104;
const Int_t DetIdSDC1 = 105;
const Int_t DetIdSDC2 = 106;
const Int_t DetIdSDC3 = 107;
const Int_t DetIdSDC4 = 108;
const Int_t DetIdSDC5 = 109;
const Int_t DetIdBFT  = 110;

const Int_t PlMinBcIn        =   1;
const Int_t PlMaxBcIn        =  12;
const Int_t PlMinBcOut       =  13;
const Int_t PlMaxBcOut       =  24;
const Int_t PlMinSdcIn       =   1;
const Int_t PlMaxSdcIn       =  10;
const Int_t PlMinSdcOut      =  31;
const Int_t PlMaxSdcOut      =  42;
const Int_t PlMInTOF         =  51; // need to change
const Int_t PlMaxTOF         =  54; // need to change
const Int_t PlOffsBc         = 100;
const Int_t PlOffsSdcIn      =   0;
const Int_t PlOffsSdcOut     =  30;
const Int_t PlOffsTOF        =  50;
const Int_t PlOffsVP         =  20;
const Int_t PlOffsTPCX       = 600;
const Int_t PlOffsTPCY       = 650;

const Int_t NumOfLayersBc     = 6;
const Int_t NumOfLayersSDC1   = 6;
const Int_t NumOfLayersSDC2   = 4;
const Int_t NumOfLayersSDC3   = 4;
const Int_t NumOfLayersSDC4   = 4;
const Int_t NumOfLayersSDC5   = 4;
const Int_t NumOfLayersBcIn   = PlMaxBcIn   - PlMinBcIn   + 1;
const Int_t NumOfLayersBcOut  = PlMaxBcOut  - PlMinBcOut  + 1;
const Int_t NumOfLayersSdcIn  = PlMaxSdcIn  - PlMinSdcIn  + 1;
const Int_t NumOfLayersSdcOut = PlMaxSdcOut - PlMinSdcOut + 1;
const Int_t NumOfLayersTOF    = PlMaxTOF    - PlMInTOF    + 1;
const Int_t NumOfLayersVP     = 5;
const Int_t NumOfLayersTPC    = 32;
const Int_t NumOfPadTPC       = 5768;
const Int_t NumOfTimeBucket   = 170;

const Int_t MaxWireBC3      =  64;
const Int_t MaxWireBC4      =  64;

const Int_t MaxWireSDC1     =  64;
const Int_t MaxWireSDC2X    =  44;
const Int_t MaxWireSDC2Y    =  44;
const Int_t MaxWireSDC3     = 128;
const Int_t MaxWireSDC4X    = 128;
const Int_t MaxWireSDC4Y    = 128;
const Int_t MaxWireSDC4     = 128;
const Int_t MaxWireSDC5X    = 128;
const Int_t MaxWireSDC5Y    =  96;

// MaxDriftLength = CellSize/2
const Double_t CellSizeBC3 = 3.0;
const Double_t CellSizeBC4 = 3.0;
const Double_t CellSizeSDC1 =  6.0;
const Double_t CellSizeSDC2 = 10.0;
const Double_t CellSizeSDC3 =  9.0;
const Double_t CellSizeSDC4 =  9.0;
const Double_t CellSizeSDC5 =  9.0;

const Int_t NumOfPlaneBFT   =   2;
const Int_t NumOfSegBFT     = 160;

// HulRm -----------------------------------------------
const Int_t NumOfHulRm   = 4;

// Matrix ----------------------------------------------
const Int_t NumOfSegSFT_Mtx = 48;

// MsT -------------------------------------------------
enum TypesMst{typeHrTdc, typeLrTdc, typeFlag, NumOfTypesMst};
const Int_t NumOfMstHrTdc = 32;
const Int_t NumOfMstLrTdc = 64;
const Int_t NumOfMstFlag  = 7;
enum dTypesMst
  {
    mstClear,
    mstAccept,
    finalClear,
    cosolationAccept,
    fastClear,
    level2,
    noDecision,
    size_dTypsMsT
  };

// Scaler ----------------------------------------------
const Int_t NumOfScaler  = 2;

#endif

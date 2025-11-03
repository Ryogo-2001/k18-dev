// -*- C++ -*-

#ifndef DETECTOR_ID_HH
#define DETECTOR_ID_HH

#include <iostream>
#include <map>
#include <vector>
#include <array>
#include <TString.h>

const int NumOfSegRayraw = 32;

const std::map<TString, std::vector<TString>> DCNameList =
{
  {"BcOut", { "BC3", "BC4" }},
  {"SdcIn", { "SDC0", "SDC1" }},
  {"SdcOut", { "SDC2", "SDC3" }},
};

// Counters ___________________________________________________________
const Int_t NumOfSegBH1    = 11;
const Int_t NumOfSegBH2    =  8;
const Int_t NumOfSegTagSF = 55;//HYPS
const Int_t NumOfSegTagPL = 10;//HYPS
const Int_t NumOfSegRF     =  1;//HYPS
const Int_t NumOfSegT0     =  1;//HYPS
const Int_t NumOfSegBAC    =  2;
const Int_t NumOfSegSAC    =  5;
const Int_t NumOfSegSCH    = 64;
const Int_t NumOfSegE_Veto =  1;//HYPS
const Int_t NumOfSegTOF    = 48; // HYPS
const Int_t NumOfSegHTOF   = 34;
const Int_t NumOfSegBVH    =  4;
const Int_t NumOfSegAC1    = 30;
const Int_t NumOfSegWC     = 12;
const Int_t NumOfSegSAC3   = 2;
const Int_t NumOfSegSFV    = 6;

// AFT
const Int_t DetIdAFT      = 112;
const Int_t NumOfPlaneAFT = 36;
const Int_t NumOfSegAFTX  = 32;
const Int_t NumOfSegAFTY  = 16;
const Int_t NumOfSegAFT   = 32;
const std::vector<int> NumOfSegAFTarr = { 32, 32, 16, 16,
					  32, 32, 16, 16,
					  32, 32, 16, 16,
					  32, 32, 16, 16,
					  32, 32, 16, 16,
					  32, 32, 16, 16,
					  32, 32, 16, 16,
					  32, 32, 16, 16,
					  32, 32, 16, 16 };
//const Int_t NumOfSegAFT[4]    = {NumOfSegAFTX, NumOfSegAFTX, NumOfSegAFTY, NumOfSegAFTY};

//const Int_t NumOfSegAFT[4]    = {NumOfSegAFTX, NumOfSegAFTX, NumOfSegAFTY, NumOfSegAFTY};

// CFT
const int NumOfPlaneCFT   =   8;
enum CFT_PLANE{CFT_U1, CFT_PHI1, CFT_V2, CFT_PHI2, CFT_U3, CFT_PHI3, CFT_V4, CFT_PHI4};
enum CFT_PLANE_{CFT_UV1, CFT_PHI1_, CFT_UV2, CFT_PHI2_, CFT_UV3, CFT_PHI3_, CFT_UV4, CFT_PHI4_};
const int NumOfSegCFT_UV1   = 426;
const int NumOfSegCFT_PHI1  = 584;
const int NumOfSegCFT_UV2   = 472;
const int NumOfSegCFT_PHI2  = 692;
const int NumOfSegCFT_UV3   = 510;
const int NumOfSegCFT_PHI3  = 800;
const int NumOfSegCFT_UV4   = 538;
const int NumOfSegCFT_PHI4  = 910;
const std::vector<int> NumOfSegCFT = {426,584,472,692,510,800,538,910};

// BGO
const Double_t BGO_X = 30.;
const Double_t BGO_Y = 25.;
const Double_t BGO_Z = 400.;
const Int_t    NumOfBGOUnit = 8;
const Int_t    NumOfBGOInOneUnit = 2;//pair unit
const Double_t RadiusOfBGOSurface = 100.;
const Int_t    NumOfBGOInOneUnit2 = 1;//single unit
const Double_t RadiusOfBGOSurface2 = 120.;
const Int_t NumOfSegBGO = NumOfBGOUnit*(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);//24

// PiID counter
const Int_t NumOfSegPiID =  32;

const Int_t NumOfPiIDUnit = 8;
const Int_t NumOfPiIDInOneUnit = 3;
const Double_t PiID_X = 30.;
const Double_t PiID_Y = 10.;
const Double_t PiID_Z = 400.;
const Double_t RadiusOfPiIDSurface = 164.;

const Int_t    NumOfPiIDInOneUnit2 = 1;//single unit
const Double_t PiID2_X = 40.;
const Double_t PiID2_Y = 10.;
const Double_t PiID2_Z = 400.;
const Double_t RadiusOfPiID2Surface = 180.;

// VMEEASIROC
const Int_t DetIdVMEASIROC = 116;
const Int_t NumOfPlaneVMEEASIROC = 96;
const Int_t NumOfSegVMEEASIROC = 64;

// Misc _______________________________________________________________
const Int_t DetIdTrig       = 21;
const Int_t DetIdScaler     = 22;
const Int_t DetIdMsT        = 25;
const Int_t DetIdMtx        = 26;
//const Int_t DetIdVmeRm      = 81; // If you use DIGIT before 2024.10.27
const Int_t DetIdVmeRm      = 101; // If you use DIGIT after 2024.10.27
const Int_t DetIdMsTRM      = 82;
const Int_t DetIdHulRM      = 83; // If you use DIGIT before 2024.10.27
const Int_t DetIdUnixTime   = 200;
const Int_t NumOfSegScaler  = 96;
const Int_t NumOfPlaneVmeRm = 2;

// Trigger Flag
namespace trigger
{
  enum ETriggerFlag
  {
    kTrigAPS,
    kTrigBPS,
    kTrigCPS,
    kTrigDPS,
    kTrigEPS,
    kTrigFPS,
    kL1SpillOn,
    kL1SpillOff,
    kSpillOnEnd,
    kSpillOffEnd,
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
     "TrigA-PS",
     "TrigB-PS",
     "TrigC-PS",
     "TrigD-PS",
     "TrigE-PS",
     "TrigF-PS",
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
// const Int_t DetIdBC3  = 103;
// const Int_t DetIdBC4  = 104;
// const Int_t DetIdSDC1 = 105;
// const Int_t DetIdSDC2 = 106;
// const Int_t DetIdSDC3 = 107;
// const Int_t DetIdSDC4 = 108;
// const Int_t DetIdSDC5 = 109;
// const Int_t DetIdBFT  = 110;

const Int_t PlMinBcIn        =   1;
const Int_t PlMaxBcIn        =  12;
const Int_t PlMinBcOut       = 113;
const Int_t PlMaxBcOut       = 124;
const Int_t PlMinSdcIn       =   1;
const Int_t PlMaxSdcIn       =  10;
const Int_t PlMinSdcOut      =  31;
const Int_t PlMaxSdcOut      =  40;
const Int_t PlMinTOF         =  53;
const Int_t PlMaxTOF         =  58;
const Int_t PlMinVP          =  16;
const Int_t PlMaxVP          =  16;
const Int_t PlOffsBc         = 100;
const Int_t PlOffsSdcIn      =   0;
const Int_t PlOffsSdcOut     =  30;
const Int_t PlOffsTOF        =  52;
const Int_t PlOffsVP         =  15;
// const Int_t PlOffsTPCX       = 600;
// const Int_t PlOffsTPCY       = 650;

const Int_t NumOfLayersBc     = 6;
const Int_t NumOfLayersSDC0   = 4;
const Int_t NumOfLayersSDC1   = 6;
const Int_t NumOfLayersSDC2   = 5;
const Int_t NumOfLayersSDC3   = 5;
const Int_t NumOfLayersSDC4   = 5;
const Int_t NumOfLayersSDC5   = 5;
const Int_t NumOfLayersBcIn   = PlMaxBcIn   - PlMinBcIn   + 1;
const Int_t NumOfLayersBcOut  = PlMaxBcOut  - PlMinBcOut  + 1;
const Int_t NumOfLayersSdcIn  = PlMaxSdcIn  - PlMinSdcIn  + 1;
const Int_t NumOfLayersSdcOut = PlMaxSdcOut - PlMinSdcOut + 1;
const Int_t NumOfLayersTOF    = PlMaxTOF    - PlMinTOF    + 1;
const Int_t NumOfLayersVP     = PlMaxVP     - PlMinVP     + 1;
// const Int_t NumOfLayersTPC    = 32;
// const Int_t NumOfPadTPC       = 5768;
// const Int_t NumOfTimeBucket   = 170;

const Int_t MaxWireBC3      =  64;
const Int_t MaxWireBC4      =  64;

const Int_t MaxWireSDC0     =  112;
const Int_t MaxWireSDC1     =  48;
const Int_t MaxWireSDC2     =  104;
const Int_t MaxWireSDC3     =  104;
const Int_t MaxWireSDC4     =  104;
const Int_t MaxWireSDC5X     =  104;
const Int_t MaxWireSDC5Y     =  104;

// MaxDriftLength = CellSize/2
const Double_t CellSizeBC3 = 3.0;
const Double_t CellSizeBC4 = 3.0;
const Double_t CellSizeSDC0 =  6.0;
const Double_t CellSizeSDC1 =  6.0;
const Double_t CellSizeSDC2 =  10.0;
const Double_t CellSizeSDC3 =  10.0;
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
const Int_t NumOfScaler  = 3;

// Parasite ___________________________________________________________
const Int_t DetIdE72BAC      =  501;
const Int_t DetIdE90SAC      =  502;
const Int_t DetIdE72KVC      =  503;
const Int_t DetIdE42BH2      =  504;
const Int_t DetIdT1          =  505;
const Int_t DetIdT2          =  506;
const Int_t NumOfSegE72BAC   =  1;
const Int_t NumOfSegE90SAC   =  2;
const Int_t NumOfSegE72KVC   =  4;
const Int_t NumOfSegE42BH2   =  8;
const Int_t NumOfSegT1       =  1;
const Int_t NumOfSegT2       =  1;

#endif

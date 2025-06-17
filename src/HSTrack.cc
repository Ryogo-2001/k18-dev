#include "HSTrack.hh"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <TMinuit.h>

#include <std_ostream.hh>

#include "ConfMan.hh"
#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "PrintHelper.hh"
#include "TrackHit.hh"
#include "K18TrackD2U.hh"

namespace {
  const auto& gGeom = DCGeomMan::GetInstance();
  const auto& gConf = ConfMan::GetInstance();
  const Double_t& zLocalK18HS = gGeom.LocalZ("K18HS");
  const Double_t& zGlobalK18HS = gGeom.GlobalZ("K18HS");
  const Double_t& pK18 = ConfMan::Get<Double_t>("PK18");
}

HSTrack::HSTrack(Double_t x, Double_t y,
		 Double_t u, Double_t v, Double_t p)
  :m_status(kInit), xInit(x), yInit(y), uInit(u), vInit(v), m_initial_momentum(p),
   m_polarity(0.), m_pid(1), m_path_length_total(0.), m_is_good(true)
{
  s_status[kInit]                = "Initialized";
  s_status[kPassed]              = "Passed";
  s_status[kExceedMaxPathLength] = "Exceed Max Path Length";
  s_status[kExceedMaxStep]       = "Exceed Max Step";
  s_status[kFailedSave]          = "Failed to Save";
  s_status[kFatal]               = "Fatal";

}

HSTrack::~HSTrack(){}

Bool_t HSTrack::Propagate() {
  m_status = kInit;

  if(m_initial_momentum <= 0){
    hddaq::cout << FUNC_NAME << " initial momentum must be exist"
		<< m_initial_momentum << std::endl;
    m_status = kFatal;
    return false;
  }

  static const auto xGlobalBcOut = gGeom.GetGlobalPosition("BC3-X1").X();
  static const auto yGlobalBcOut = gGeom.GetGlobalPosition("BC3-X1").Y();
  static const auto zGlobalBcOut = gGeom.GetGlobalPosition("BC3-X1").Z();
  static const auto zLocalBcOut = gGeom.GetLocalZ("BC3-X1");

  const Double_t uIn = uInit;
  const Double_t vIn = vInit;
  const Double_t xIn = xInit;
  const Double_t yIn = yInit;
  const Double_t pz = m_initial_momentum/std::sqrt(1. + uIn*uIn+vIn*vIn);
  const ThreeVector posIn(xGlobalBcOut + xIn, yGlobalBcOut + yIn, zGlobalBcOut - zLocalBcOut);
  const ThreeVector momIn(pz*uIn , pz*vIn, pz);
  m_v0_position.SetX(xInit);
  m_v0_position.SetY(yInit);
  m_v0_position.SetZ(zGlobalBcOut - zLocalBcOut - zGlobalK18HS);
  m_v0_momentum = momIn;

  RKCordParameter iniCord(posIn, momIn);
  m_polarity = pK18/TMath::Abs(pK18);
  if(m_polarity>0) iniCord.ConvertCharge();
  RKCordParameter prevCord;
  RKHitPointContainer preHPntCont;

  m_HitPointCont = RK::MakeHSHPContainer();
  RKHitPointContainer prevHPntCont;
  m_status = (RKstatus)RK::Extrap(iniCord, m_HitPointCont, m_pid);
  //std::cout << s_status[m_status] << std::endl;
  if(m_status != kPassed || !SaveTrack()) return false;

  return true;
}

//_____________________________________________________________________________
Bool_t
HSTrack::SaveTrack()
{

  for(Int_t i=0; i<NumOfLayersBcOut; ++i){
    Int_t IdBC = i+PlOffsBcOut+1;
    Double_t zBC = gGeom.GetLocalZ(IdBC);
    const RKcalcHitPoint& hpBC = m_HitPointCont.HitPointOfLayer(IdBC);
    const ThreeVector &posBC = hpBC.PositionInGlobal();
    const ThreeVector &momBC = hpBC.MomentumInGlobal();
    m_bc_position.push_back(gGeom.Global2LocalPos(IdBC, posBC));
    m_bc_position[i].SetZ(zBC - zLocalK18HS);
    m_bc_momentum.push_back(gGeom.Global2LocalDir(IdBC, momBC));
  }

  const Int_t IdBH2 = gGeom.DetectorId("BH2");
  const Double_t& zBH2 = gGeom.LocalZ("BH2");
  const RKcalcHitPoint& hpBH2 = m_HitPointCont.HitPointOfLayer(IdBH2);
  const ThreeVector &posBH2 = hpBH2.PositionInGlobal();
  const ThreeVector &momBH2 = hpBH2.MomentumInGlobal();
  m_bh2_position = gGeom.Global2LocalPos(IdBH2, posBH2);
  m_bh2_position.SetZ(zBH2 - zLocalK18HS);
  m_bh2_momentum = gGeom.Global2LocalDir(IdBH2, momBH2);

  const Int_t IdBAC = gGeom.DetectorId("BAC");
  const Double_t& zBAC = gGeom.LocalZ("BAC");
  const RKcalcHitPoint& hpBAC = m_HitPointCont.HitPointOfLayer(IdBAC);
  const ThreeVector &posBAC = hpBAC.PositionInGlobal();
  const ThreeVector &momBAC = hpBAC.MomentumInGlobal();
  m_bac_position = gGeom.Global2LocalPos(IdBAC, posBAC);
  m_bac_position.SetZ(zBAC - zLocalK18HS);
  m_bac_momentum = gGeom.Global2LocalDir(IdBAC, momBAC);

  const Int_t IdVP1 = gGeom.DetectorId("VPHS1");
  const Double_t& zVPHS1 = gGeom.LocalZ("VPHS1");
  const RKcalcHitPoint& hpVP1 = m_HitPointCont.HitPointOfLayer(IdVP1);
  const ThreeVector &posVP1 = hpVP1.PositionInGlobal();
  const ThreeVector &momVP1 = hpVP1.MomentumInGlobal();
  m_vp1_position = gGeom.Global2LocalPos(IdVP1, posVP1);
  m_vp1_position.SetZ(zVPHS1 - zLocalK18HS);
  m_vp1_momentum = gGeom.Global2LocalDir(IdVP1, momVP1);

  const Int_t IdVP2 = gGeom.DetectorId("VPHS2");
  const Double_t& zVPHS2 = gGeom.LocalZ("VPHS2");
  const RKcalcHitPoint& hpVP2 = m_HitPointCont.HitPointOfLayer(IdVP2);
  const ThreeVector &posVP2 = hpVP2.PositionInGlobal();
  const ThreeVector &momVP2 = hpVP2.MomentumInGlobal();
  m_vp2_position = gGeom.Global2LocalPos(IdVP2, posVP2);
  m_vp2_position.SetZ(zVPHS2 - zLocalK18HS);
  m_vp2_momentum = gGeom.Global2LocalDir(IdVP2, momVP2);

  const Int_t IdVP3 = gGeom.DetectorId("VPHS3");
  const Double_t& zVPHS3 = gGeom.LocalZ("VPHS3");
  const RKcalcHitPoint& hpVP3 = m_HitPointCont.HitPointOfLayer(IdVP3);
  const ThreeVector &posVP3 = hpVP3.PositionInGlobal();
  const ThreeVector &momVP3 = hpVP3.MomentumInGlobal();
  m_vp3_position = gGeom.Global2LocalPos(IdVP3, posVP3);
  m_vp3_position.SetZ(zVPHS3 - zLocalK18HS);
  m_vp3_momentum = gGeom.Global2LocalDir(IdVP3, momVP3);

  const Int_t IdVP4 = gGeom.DetectorId("VPHS4");
  const Double_t& zVPHS4 = gGeom.LocalZ("VPHS4");
  const RKcalcHitPoint& hpVP4 = m_HitPointCont.HitPointOfLayer(IdVP4);
  const ThreeVector &posVP4 = hpVP4.PositionInGlobal();
  const ThreeVector &momVP4 = hpVP4.MomentumInGlobal();
  m_vp4_position = gGeom.Global2LocalPos(IdVP4, posVP4);
  m_vp4_position.SetZ(zVPHS4 - zLocalK18HS);
  m_vp4_momentum = gGeom.Global2LocalDir(IdVP4, momVP4);

  const Int_t TGTid = gGeom.DetectorId("K18Target");
  const Double_t& zK18Tgt = gGeom.LocalZ("K18Target");
  const RKcalcHitPoint& hpTgt  = m_HitPointCont.HitPointOfLayer(TGTid);
  const ThreeVector& posTgt = hpTgt.PositionInGlobal();
  const ThreeVector& momTgt = hpTgt.MomentumInGlobal();
  m_tgt_position = gGeom.Global2LocalPos(TGTid, posTgt);
  m_tgt_position.SetZ(zK18Tgt - zLocalK18HS);
  m_tgt_momentum = gGeom.Global2LocalDir(TGTid, momTgt);

  const Int_t IdGasVesselU = gGeom.DetectorId("K18VesselU");
  const Double_t& zGasVesselU = gGeom.LocalZ("K18VesselU");
  const RKcalcHitPoint& hpGasVesselU = m_HitPointCont.HitPointOfLayer(IdGasVesselU);
  const ThreeVector &posGasVesselU = hpGasVesselU.PositionInGlobal();
  const ThreeVector &momGasVesselU = hpGasVesselU.MomentumInGlobal();
  m_gasvesselU_position = gGeom.Global2LocalPos(IdGasVesselU, posGasVesselU);
  m_gasvesselU_position.SetZ(zGasVesselU - zLocalK18HS);
  m_gasvesselU_momentum = gGeom.Global2LocalDir(IdGasVesselU, momGasVesselU);

  const Int_t IdGasVesselD = gGeom.DetectorId("K18VesselD");
  const Double_t& zGasVesselD = gGeom.LocalZ("K18VesselD");
  const RKcalcHitPoint& hpGasVesselD = m_HitPointCont.HitPointOfLayer(IdGasVesselD);
  const ThreeVector &posGasVesselD = hpGasVesselD.PositionInGlobal();
  const ThreeVector &momGasVesselD = hpGasVesselD.MomentumInGlobal();
  m_gasvesselD_position = gGeom.Global2LocalPos(IdGasVesselD, posGasVesselD);
  m_gasvesselD_position.SetZ(zGasVesselD - zLocalK18HS);
  m_gasvesselD_momentum = gGeom.Global2LocalDir(IdGasVesselD, momGasVesselD);

  const Int_t IdHtof = gGeom.DetectorId("HTOF");
  const Double_t& zHtof = gGeom.LocalZ("HTOF");
  const RKcalcHitPoint &hpHtof = m_HitPointCont.rbegin()->second;
  const ThreeVector &posHtof = hpHtof.PositionInGlobal();
  const ThreeVector &momHtof = hpHtof.MomentumInGlobal();
  m_htof_position = gGeom.Global2LocalPos(IdHtof, posHtof);
  m_htof_position.SetZ(zHtof - zLocalK18HS);
  m_htof_momentum = gGeom.Global2LocalDir(IdHtof, momHtof);

  m_path_length_total = std::abs(hpBH2.PathLength()-hpTgt.PathLength());

  return true;
}

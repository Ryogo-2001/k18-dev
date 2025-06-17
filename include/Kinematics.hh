// -*- C++ -*-

#ifndef KINEMATICS_HH
#define KINEMATICS_HH

#include <TString.h>
#include <TVector3.h>

//_____________________________________________________________________________
namespace Kinematics
{
  const TString& ClassName();
  Double_t MassSquare(Double_t p, Double_t path, Double_t time);
  Double_t CalcTimeOfFlight(Double_t p, Double_t path, Double_t mass);
  TVector3 VertexPoint(const TVector3& Xin, const TVector3& Xout,
		       const TVector3& Pin, const TVector3& Pout);
  TVector3 VertexPointByHonly(const TVector3& Xin,
			      const TVector3& Xout,
			      const TVector3& Pin,
			      const TVector3& Pout);
  TVector3 VertexPoint(const TVector3& Xin, const TVector3& Xout,
		       const TVector3& Pin, const TVector3& Pout,
		       Double_t& dist);
  TVector3 VertexPoint3D(const TVector3& Xin, const TVector3& Xout,
			 const TVector3& Pin, const TVector3& Pout,
			 Double_t& dist);
  TVector3 VertexPointReal(const TVector3& Xin, const TVector3& Xout,
			   const TVector3& Pin, const TVector3& Pout,
			   Double_t& dist);
  TVector3 VertexPointTF2(const TVector3& Xin, const TVector3& Xout,
			  const TVector3& Pin, const TVector3& Pout,
			  Double_t& dist);
  TVector3 VertexPointHelix(const Double_t par1[5], const Double_t par2[5],
                            Double_t& dist, Double_t& t1, Double_t& t2);
  Double_t CloseDist(const TVector3& Xin, const TVector3& Xout,
		     const TVector3& Pin, const TVector3& Pout);
  TVector3 CorrElossIn(const TVector3& Pin,
		       const TVector3& Xin,
		       const TVector3& vtx, Double_t mass);
  Double_t CalcLengthBeam(const TVector3& Pin,
			  const TVector3& Xin,
			  const TVector3& vtx);
  TVector3 CorrElossOut(const TVector3& Pout,
			const TVector3& Xout,
			const TVector3& vtx, Double_t mass);
  Double_t CalcLengthScat(const TVector3& Pout,
			  const TVector3& Xout,
			  const TVector3& vtx);
  TVector3 CorrElossOutCheck(const TVector3& Pout,
			     const TVector3& Xout,
			     const TVector3& vtx, Double_t mass);
  Bool_t   IsInsideTarget(const TVector3&point);
  Bool_t   CalcCrossingPoint(Double_t u, Double_t v, const TVector3& Point,
			     Double_t *z1, Double_t *z2);
  Double_t DiffE(Double_t mass, Double_t E, Double_t length, Double_t Elast);
  Int_t    CalcDe(Double_t momentum, Double_t mass, Double_t distance,
		  Double_t *momentum_cor, Double_t *energy_cor);
  Double_t CalcDedx(Double_t beta);
  Double_t Gamma(Double_t beta);
  Double_t Beta(Double_t energy, Double_t mormentum);

  //For HypTPC dE calculation
  //Material lookup table
  //0:P10, 1:Polyethylene(Target) 2:Diamond(Target) 3:Polyvinyltoluene (old TPC gas-vessel) 4:Mylar (gas vessel window) 5:Al (gas-vessel frame)
  Double_t DensityEffectCorrection(Double_t betagamma, Double_t *par);
  Double_t HypTPCdEdx(Int_t target_material, Double_t mass, Double_t beta);
  Int_t HypTPCCalcDe(Int_t target_material, Double_t momentum,
		     Double_t mass, Double_t distance,
		     Double_t *momentum_cor, Double_t *energy_cor);
  Double_t HypTPCDiffE(Int_t materialid, Double_t mass,
		       Double_t E, Double_t length, Double_t Elast);
  TVector3 HypTPCCorrElossOut(Int_t target_material, const TVector3& Pout,
			      Double_t length, Double_t mass);
  TVector3 HypTPCCorrElossIn(Int_t target_material, const TVector3& Pin,
			     Double_t length, Double_t mass);

  //For HypTPC dE/dx pid
  Double_t HypTPCBethe(Double_t *x, Double_t *p);
  Bool_t HypTPCdEdxPID_IsKaonTemp(Double_t dedx, Double_t poq); //temporary
  Double_t HypTPCdEdxNsigmaProton(Double_t dedx, Double_t poq);
  Double_t HypTPCdEdxNsigmaDeutron(Double_t dedx, Double_t poq);
  Double_t HypTPCdEdxNsigmaTriton(Double_t dedx, Double_t poq);
  Double_t HypTPCdEdxNsigmaKaon(Double_t dedx, Double_t poq);
  Double_t HypTPCdEdxNsigmaPion(Double_t dedx, Double_t poq);
  Double_t HypTPCdEdxNsigmaElectron(Double_t dedx, Double_t poq);
  Bool_t HypTPCdEdxElectron(Double_t dedx, Double_t poq);
  Bool_t HypTPCdEdxKaon(Double_t dedx, Double_t poq);
  Int_t HypTPCdEdxPID(Double_t dedx, Double_t poq);
  void HypTPCPID_PDGCode(Int_t charge, Int_t pid, std::vector<Int_t>& pdg);

  //For HTOF pid
  Double_t HypTPCHTOFNsigmaProton(Double_t poq, Double_t tracklength, Double_t tof);
  Double_t HypTPCHTOFNsigmaDeutron(Double_t poq, Double_t tracklength, Double_t tof);
  Double_t HypTPCHTOFNsigmaTriton(Double_t poq, Double_t tracklength, Double_t tof);
  Double_t HypTPCHTOFNsigmaKaon(Double_t poq, Double_t tracklength, Double_t tof);
  Double_t HypTPCHTOFNsigmaPion(Double_t poq, Double_t tracklength, Double_t tof);
  Double_t HypTPCHTOFNsigmaElectron(Double_t poq, Double_t tracklength, Double_t tof);
  Double_t HypTPCHTOFNsigmaProton(Double_t poq, Double_t inverse_beta);
  Double_t HypTPCHTOFNsigmaDeutron(Double_t poq, Double_t inverse_beta);
  Double_t HypTPCHTOFNsigmaTriton(Double_t poq, Double_t inverse_beta);
  Double_t HypTPCHTOFNsigmaKaon(Double_t poq, Double_t inverse_beta);
  Double_t HypTPCHTOFNsigmaPion(Double_t poq, Double_t inverse_beta);
  Double_t HypTPCHTOFNsigmaElectron(Double_t poq,Double_t inverse_beta);

  //For helix tracking & vertex reconstruction
  Double_t CalcHelixCloseDist(TVector3 point, Double_t par[5], Double_t t1_start, Double_t t1_end);
  TVector3 CalcHelixMom(Double_t Bfield, Int_t charge,
			Double_t par[5], Double_t t);
  void CalcHelixParam(Double_t Bfield, Int_t charge,
		      TVector3 mom, TVector3 pos, Double_t *par);
  TVector3 CalcHelixPosition(double par[5], double t);
  TVector3 VertexPointHelix(Double_t par1[5], Double_t par2[5],
			    Double_t t1_start, Double_t t1_end,
			    Double_t t2_start, Double_t t2_end,
			    Double_t& t1, Double_t& t2,
			    Double_t& dist);
  TVector3 LambdaVertex(Double_t Bfield, Double_t p_par[5], Double_t pi_par[5],
			Double_t t1_start, Double_t t1_end,
			Double_t t2_start, Double_t t2_end,
			TVector3 &p_mom, TVector3 &pi_mom,
			TVector3 &lambda_mom, Double_t& dist);
  TVector3 XiVertex(Double_t Bfield, Double_t pi_par[5],
		    Double_t t_start, Double_t t_end,
		    TVector3 Xlambda, TVector3 Plambda,
		    TVector3 &Ppi, Double_t &lambdapi_dist);
  TVector3 CalcCloseDistXi(TVector3 point, Double_t Bfield,
			   TVector3 xi_decayvtx, TVector3 xi_mom_decayvtx,
			   TVector3 &xi_mom, Double_t &dist); //Extrapolate to the target
  TVector3 LambdaLambdaVertex(TVector3 Xlambda1, TVector3 Plambda1,
			      TVector3 Xlambda2, TVector3 Plambda2,
			      TVector3 &vtxlambda1, TVector3 &vtxlambda2,
			      Double_t &dist);
  TVector3 CalcCloseDistLambda(TVector3 point, TVector3 Xlambda,
			       TVector3 Plambda, Double_t &dist); //Extrapolate to the target
  TVector3 LambdaTargetCenter(TVector3 Xlambda,
			      TVector3 Plambda, Double_t &dist);  //Extrapolate to z=z_target
  Bool_t HelixDirection(TVector3 vertex, TVector3 start, TVector3 end,
			Double_t &dist);
  TVector3 MultitrackVertex(Int_t ntrack, Double_t *x0, Double_t *y0,
			    Double_t *u0, Double_t *v0,
			    std::vector<Double_t> Res_x0 = {}, std::vector<Double_t> Res_y0 = {},
			    std::vector<Double_t> Res_u0 = {}, std::vector<Double_t> Res_v0 = {});
  TVector3 MultitrackVertex(Int_t ntrack, Double_t *x0, Double_t *y0,
			    Double_t *u0, Double_t *v0,
			    std::vector<Double_t> Res_x0, std::vector<Double_t> Res_y0,
			    std::vector<Double_t> Res_u0, std::vector<Double_t> Res_v0,
			    Double_t &chisqr);

}

//_____________________________________________________________________________
inline const TString&
Kinematics::ClassName()
{
  static TString s_name("Kinematic");
  return s_name;
}

#endif

#include "FourVectorCartesianFitter.hh"
#ifndef FourVectorCartesianFitter_cc
#define FourVectorCartesianFitter_cc
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit

void FourVectorCartesianFitter::UseVertex(bool status,TVector3 Vert1,TVector3 Vert2){
	MeasDir = status;
	Clear();
	if(MeasDir){//Production and Decay vertex of R is known, hence the direction is known. However, the momentum magnitude is not directly measured. I used the sum of P and Q as an initial value.
		auto dir = Vert1 - Vert2;	
		auto th = dir.Theta();
		auto ph = dir.Phi();
		R.SetTheta(th);
		R.SetPhi(ph);
	}
	Initialize();
}
void FourVectorCartesianFitter::Initialize(){
	Clear();
	if(MeasDir){
		nMeas = 8;nUnkn = 1; nConst = 4; 
	}
	else{
		nMeas = 6;nUnkn = 3; nConst = 4; //For 1-C fit, when Vertex information is useless.
	}
	mP = P.Mag();
	TVector3 TV_P = P.Vect();
	double x_P = TV_P.x();	
	double y_P = TV_P.y();	
	double z_P = TV_P.z();	
	
	mQ = Q.Mag();
	TVector3 TV_Q = Q.Vect(); 
	double x_Q = TV_Q.x();	
	double y_Q = TV_Q.y();	
	double z_Q = TV_Q.z();	

	TVector3 TV_R = R.Vect(); 
	double x_R = TV_R.x();	
	double y_R = TV_R.y();	
	double z_R = TV_R.z();	
	double p_R = TV_R.Mag();	
	double th_R = TV_R.Theta();
	double ph_R = TV_R.Phi();

	double meas[8];
	double unkn[3];
	if(MeasDir){
		double temp[] = {th_R,ph_R,x_P,y_P,z_P,x_Q,y_Q,z_Q};
		for(int i=0;i<nMeas;++i)meas[i]=temp[i];
		double temp2[] = {p_R,0,0};
		for(int i=0;i<nUnkn;++i)unkn[i]=temp2[i];
	}
	else{
		double temp[] = {x_P,y_P,z_P,x_Q,y_Q,z_Q,0,0};
		for(int i=0;i<nMeas;++i)meas[i]=temp[i];
		double temp2[] = {x_R,y_R,z_R};
		for(int i=0;i<nUnkn;++i)unkn[i]=temp2[i];
	}
	TMatrixD Meas0(nMeas,1,meas);	
	TMatrixD Unkn0(nUnkn,1,unkn);
	vector<double>Pull;
	Pull.resize(nMeas);
	vector<double>UPull;
	UPull.resize(nUnkn);
	Measurements.push_back(Meas0);
	Unknowns.push_back(Unkn0);
	Pulls.push_back(Pull);
	UPulls.push_back(UPull);
	Chi2s.push_back(-1);
	MassDiffs.push_back(1e9);
}
void FourVectorCartesianFitter::SetConstraints(){
	auto Meas = Measurements.at(step); 
	auto Unkn = Unknowns.at(step);
	double p_R,th_R,ph_R,x_R,y_R,z_R,p_P,x_P,y_P,z_P,p_Q,x_Q,y_Q,z_Q;
	if(MeasDir){
		p_R= Unkn(0,0); 
		th_R= Meas(0,0); 
		ph_R= Meas(1,0); 
		x_P= Meas(2,0); 
		y_P= Meas(3,0); 
		z_P= Meas(4,0);
		x_Q= Meas(5,0); 
		y_Q= Meas(6,0); 
		z_Q= Meas(7,0); 
	
		x_R = p_R*sin(th_R)*cos(ph_R) ;
		y_R = p_R*sin(th_R)*sin(ph_R) ;
		z_R = p_R*cos(th_R) ;
		p_P = hypot(x_P,hypot(y_P,z_P));
		p_Q = hypot(x_Q,hypot(y_Q,z_Q));
	}
	else{
		x_R= Unkn(0,0); 
		y_R= Unkn(1,0); 
		z_R= Unkn(2,0); 
		x_P= Meas(0,0); 
		y_P= Meas(1,0); 
		z_P= Meas(2,0);
		x_Q= Meas(3,0); 
		y_Q= Meas(4,0); 
		z_Q= Meas(5,0); 
		p_P = hypot(x_P,hypot(y_P,z_P));
		p_Q = hypot(x_Q,hypot(y_Q,z_Q));
		p_R = hypot(x_R,hypot(y_R,z_R));
	}
	double f1 = x_P+x_Q-x_R;
	double f2 = y_P+y_Q-y_R;
	double f3 = z_P+z_Q-z_R;
	double E_P = hypot(p_P,mP);
	double E_Q = hypot(p_Q,mQ);
	double E_R = hypot(p_R,mR);
	double f4	= E_P+E_Q-E_R;

	double df1du1 = -1; 
	double df1du2 = 0; 
	double df1du3 = 0; 

	double df2du1 = 0; 
	double df2du2 = -1; 
	double df2du3 = 0; 

	double df3du1 = 0; 
	double df3du2 = 0; 
	double df3du3 = -1; 


	//E = sqrt(px^2 + ... + m^2);
	//dEdp_x = px/E
	double df4du1=-x_R/E_R;
	double df4du2=-y_R/E_R;
	double df4du3=-z_R/E_R;

	double df1dm1 = 1;
	double df1dm2 = 0;
	double df1dm3 = 0;
	double df1dm4 = 1;
	double df1dm5 = 0;
	double df1dm6 = 0;
	
	double df2dm1 = 0;
	double df2dm2 = 1;
	double df2dm3 = 0;
	double df2dm4 = 0;
	double df2dm5 = 1;
	double df2dm6 = 0;
	
	double df3dm1 = 0;
	double df3dm2 = 0;
	double df3dm3 = 1;
	double df3dm4 = 0;
	double df3dm5 = 0;
	double df3dm6 = 1;
	
	double df4dm1=x_P/E_P;
	double df4dm2=y_P/E_P;
	double df4dm3=z_P/E_P;
	double df4dm4=x_Q/E_Q;
	double df4dm5=y_Q/E_Q;
	double df4dm6=z_Q/E_Q;



	double fs[]={f1,f2,f3,f4};
	double dfdms[200] ;
	double dfdus[200] ;
	if(MeasDir){
		double temp[] = {
		};
		for(int i=0;i<nMeas*nConst;++i){
			dfdms[i]=temp[i];
		};
		double tempu[] = {
		};	
		for(int i=0;i<nUnkn*nConst;++i){
			dfdus[i]=tempu[i];
		};
		double temp1[]={
		};
		double temp2[]={
		};
		double temp3[]={
		};
		double temp4[]={
		};
		TMatrixD d2F1dU(nUnkn,nUnkn,temp1);
		TMatrixD d2F2dU(nUnkn,nUnkn,temp2);
		TMatrixD d2F3dU(nUnkn,nUnkn,temp3);
		TMatrixD d2F4dU(nUnkn,nUnkn,temp4);
		vector<TMatrixD> d2FdU = {
			d2F1dU,
			d2F2dU,
			d2F3dU,
			d2F4dU
		};
		d2Fd2Us.push_back(d2FdU);
	}
	else{
		double temp[] = {
			df1dm1,df1dm2,df1dm3,df1dm4	,df1dm5,df1dm6,
			df2dm1,df2dm2,df2dm3,df2dm4	,df2dm5,df2dm6,
			df3dm1,df3dm2,df3dm3,df3dm4	,df3dm5,df3dm6,
			df4dm1,df4dm2,df4dm3,df4dm4	,df4dm5,df4dm6
		};
		for(int i=0;i<nMeas*nConst;++i){
			dfdms[i]=temp[i];
		};
		double tempu[] = {
			df1du1,df1du2,df1du3,
			df2du1,df2du2,df2du3,
			df3du1,df3du2,df3du3,
			df4du1,df4du2,df4du3
		};
		for(int i=0;i<nUnkn*nConst;++i){
			dfdus[i]=tempu[i];
		};
		double temp1[] = {
//			df1du1du1,df1du1du2,df1du1du3,
//			df1du2du1,df1du2du2,df1du2du3,
//			df1du3du1,df1du3du2,df1du3du3
		};
		double temp2[] = {
//			df2du1du1,df2du1du2,df2du1du3,
//			df2du2du1,df2du2du2,df2du2du3,
//			df2du3du1,df2du3du2,df2du3du3
		};
		double temp3[] = {
//			df3du1du1,df3du1du2,df3du1du3,
//			df3du2du1,df3du2du2,df3du2du3,
//			df3du3du1,df3du3du2,df3du3du3
		};
		double temp4[] = {
//			df4du1du1,df4du1du2,df4du1du3,
//			df4du2du1,df4du2du2,df4du2du3,
//			df4du3du1,df4du3du2,df4du3du3
		};
		TMatrixD d2F1dU(nUnkn,nUnkn,temp1);
		TMatrixD d2F2dU(nUnkn,nUnkn,temp2);
		TMatrixD d2F3dU(nUnkn,nUnkn,temp3);
		TMatrixD d2F4dU(nUnkn,nUnkn,temp4);
		vector<TMatrixD> d2FdU = {
			d2F1dU,
			d2F2dU,
			d2F3dU,
			d2F4dU
		};
		d2Fd2Us.push_back(d2FdU);
	}

	TMatrixD FMat(nConst,1,fs);
	TMatrixD dFdM(nConst,nMeas,dfdms);
	TMatrixD dFdU(nConst,nUnkn,dfdus);
	FMats.push_back(FMat);//Constraint Matrices for Each step
	dFdMs.push_back(dFdM);//Constraint matrix differentiated by measurement params.
	dFdUs.push_back(dFdU);// same, but for unmeasured params.
}
void FourVectorCartesianFitter::SampleStepPoint(int steps){
	auto Meas = Measurements.at(steps); 
	auto Unkn = Unknowns.at(steps); 
	double p_R,th_R,ph_R,x_R,y_R,z_R,p_P,x_P,y_P,z_P,p_Q,x_Q,y_Q,z_Q;
	if(MeasDir){
		p_R= Unkn(0,0); 
		th_R= Meas(0,0); 
		ph_R= Meas(1,0); 
		x_P= Meas(2,0); 
		y_P= Meas(3,0); 
		z_P= Meas(4,0);
		x_Q= Meas(5,0); 
		y_Q= Meas(6,0); 
		z_Q= Meas(7,0); 
	
		x_R = p_R*sin(th_R)*cos(ph_R); 
		y_R = p_R*sin(th_R)*sin(ph_R); 
		z_R = p_R*cos(th_R); 
		p_P = hypot(x_P,hypot(y_P,z_P));
		p_Q = hypot(x_Q,hypot(y_Q,z_Q));
	}
	else{
		x_R= Unkn(0,0); 
		y_R= Unkn(1,0); 
		z_R= Unkn(2,0); 
		x_P= Meas(0,0); 
		y_P= Meas(1,0); 
		z_P= Meas(2,0);
		x_Q= Meas(3,0); 
		y_Q= Meas(4,0); 
		z_Q= Meas(5,0); 
		p_P = hypot(x_P,hypot(y_P,z_P));
		p_Q = hypot(x_Q,hypot(y_Q,z_Q));
		p_R = hypot(x_R,hypot(y_R,z_R));
	}
	TLorentzVector PP(x_P,y_P,z_P,hypot(mP,p_P));
	TLorentzVector QQ(x_Q,y_Q,z_Q,hypot(mQ,p_Q));
	TLorentzVector RR(x_R,y_R,z_R,hypot(mR,p_R));
	auto V = PP + QQ;
	double MassDiff = V.Mag()-mR;
	PCor = PP;
	QCor = QQ;
	RCor = RR;
	MassDiffs.push_back(MassDiff);
}
FourVectorCartesianFitter::FourVectorCartesianFitter(TLorentzVector P_,TLorentzVector Q_, TLorentzVector R_){
	//Does Kinematic fitting for R -> P+Q decay.
	P=P_;
	Q=Q_;
	R=R_;
	Initialize();
};
TMatrixD
FourVectorCartesianFitter::JacobianSphToCart(double p, double th, double ph){
	// x = p sin(th) cos(ph)
	// y = p sin(th) sin(ph)
	// z = p cos(th)
	//V_c = J^T V J|->
	//		dxdp, dxdth,dxdph
	//J	=	dydp, dydth,dydph
	//		dzdp, dzdth,dzdph

	double dxdp = sin(th)*cos(ph);
	double dydp = sin(th)*sin(ph);
	double dzdp = cos(th);

	double dxdth = p*cos(th)*cos(ph);
	double dydth = p*cos(th)*sin(ph);
	double dzdth = -p*sin(th);

	double dxdph = -p*sin(th)*sin(ph);
	double dydph = p*sin(th)*cos(ph);
	double dzdph = 0;
	double mat[9] = 
	{ dxdp, dxdth, dxdph,
		dydp, dydth, dydph,
		dzdp, dzdth, dzdph
	};
	/*
	double mat[9] = 
	{ dxdp, dydp, dzdp,
		dxdth, dydth, dzdth,
		dxdph, dydph, dzdph
	};
	*/
	return TMatrixD(3,3,mat);

}
void
FourVectorCartesianFitter::CalcVariance(int istep){
	
	TMatrixD VMat = Variancies.at(istep);
	TMatrixD dV = dVMats.at(istep);
	double El_contract[18]={
		1,0,0,1,0,0,
		0,1,0,0,1,0,
		0,0,1,0,0,1
	};
	TMatrixD ContT(3,6,El_contract);
	TMatrixD Cont = TransposeMatrix(ContT);
	TMatrixD UVMat = ContT*VMat*Cont;
//	VarianciesU.push_back(UVMat);
}















void
FourVectorCartesianFitter::Rotate(){
	auto VMat = Variancies.at(0);
	Initialize();
	Variancies.push_back(VMat);
	TMatrixD J;
	RotateVariance(J);
}
void
FourVectorCartesianFitter::ToDecayPlane(){
	auto Zaxis =(P + Q).Vect();
	auto vP = P.Vect();
	auto vQ = Q.Vect();
	auto Yaxis = vP.Cross(vQ);
//	double YNorm = 1./(Yaxis.Mag());
//	Yaxis = YNorm * Yaxis;
	double Th_F = Zaxis.Theta();
	double Ph_F = Zaxis.Phi();
	double RotZ[9] ={
		cos(-Ph_F),	-sin(-Ph_F),	0,	
		sin(-Ph_F),	cos(-Ph_F),		0,
		0,					0,						1
	};
	double RotY[9] ={
		cos(Th_F),	0,					-sin(Th_F),
		0,					-1,					0,
		sin(Th_F),	0,					cos(Th_F)
	};
	TMatrixD RZ(3,3,RotZ);
	TMatrixD RY(3,3,RotY);
	TMatrixD R_F = RY * RZ;
	Yaxis = R_F * Yaxis;
	double Th_Y = Yaxis.Theta();
	double Ph_Y = Yaxis.Phi();
	double RotX[9] ={
		1,				0,					0,
		0,				cos(-Ph_Y),	-sin(-Ph_Y),
		0,				sin(-Ph_Y),	cos(-Ph_Y)
	};
	TMatrixD RX(3,3,RotX);
	Yaxis = RX * Yaxis;

}
#endif

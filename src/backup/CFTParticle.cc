#include "CFTParticle.hh"
#include "DCLTrackHit.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "FiberHit.hh"
#include <TRandom.h>
#include "HodoAnalyzer.hh"
#include "RawData.hh"

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);
CFTParticle::CFTParticle(DCLocalTrack *track, RawData *rawData)
  : Track_(track), RawData_(rawData), 
    CFTVtx_(-999, -999, -999),
    m_bgo_seg(-1),
    m_piid_seg(-1)
{  
Calculate();
}

CFTParticle::~CFTParticle()
{
}

Hodo1Hit * CFTParticle::GetBGOHit(int i)
{
  if (i>=0 && i<BGOCont_.size())
    return BGOCont_[i];
  else
    return 0;
}

Hodo2Hit * CFTParticle::GetPiVHit(int i)
{
  if (i>=0 && i<PiVCont_.size())
    return PiVCont_[i];
  else
    return 0;
}

bool CFTParticle::Calculate()
{
  /* here!
  FiberTotal_E_ = Track_->TotalDEHiGain();
  FiberMax_E_   = Track_->MaxDEHiGain();
  PathLength_   = Track_->GetTotalPathLength();
  NormalizedFiberTotal_E_ = Track_->NormalizedTotalDEHiGain();
  NormalizedFiberMax_E_   = Track_->NormalizedMaxDEHiGain();
  */


  int xyFlag = Track_->GetCFTxyFlag();
  int zFlag  = Track_->GetCFTzFlag() ;
  double Axy = Track_->GetAxy(); double Bxy = Track_->GetBxy();
  double Az  = Track_->GetAz() ; double Bz  = Track_->GetBz();
  ThreeVector Pos0 = Track_->GetPos0();
  ThreeVector Dir  = Track_->GetDir();

  // track to BGO segment
  double distBGO=1000.;
  HodoAnalyzer hodoAna;
  hodoAna.DecodeBGOHits(RawData_);
  hodoAna.DecodePiIDHits(RawData_);

  // BGO
  int nhBGO = hodoAna.GetNHitsBGO();
  double max_adc = 0;
  for (int i=0; i<nhBGO; i++) {
    Hodo1Hit* hit = hodoAna.GetHitBGO(i);
    int seg = hit->SegmentId();
    int nh = hit->GetNumOfHit();
    //double adc = hit->GetAUp();
    double adc = hit->DeltaE();
    bool hit_flag = false;

    if(adc>0){
      for(int m=0; m<nh; ++m){
	double cmt = hit->CMeanTime(m);
	if(-50<cmt&&cmt<50){
	  hit_flag = true;
	}
      }
    }
    if(!hit_flag)continue;
    double x, y;
    BGOPos(seg, &x, &y);
    if     (xyFlag==0){distBGO = (Axy*x-y+Bxy)/sqrt(Axy*Axy+1.*1.);}
    else if(xyFlag==1){distBGO = (Axy*y-x+Bxy)/sqrt(Axy*Axy+1.*1.);}
    double u=Dir.x(), v=Dir.y();
    double x0=Pos0.x(), y0=Pos0.y();
    double t = (u*(x-x0)+v*(y-y0))/(u*u+v*v);
    if (t>=0) {	
      if (fabs(distBGO)<25) {
	AddBGOHit(hit);
	if(adc>max_adc){
	  max_adc = adc;
	  m_bgo_seg = seg;
	}
      }
    }
  }

  // track to PiID segment
  double distPiID=1000.;
  int nhPiID = hodoAna.GetNHitsPiID();
  for(int i = 0; i<nhPiID; ++i){
    const FiberHit* hit = hodoAna.GetHitPiID(i);
    int mhit = hit->GetNLeading();
    int seg = hit->PairId();
    bool hit_flag = false;
    for(int m = 0; m<mhit; ++m){	
      double ctime  = hit->GetCTime(m);
      if(ctime>-50&&ctime<50){hit_flag=true;}    
    }
    if(!hit_flag)continue;
    double x, y;
    PiIDPos(seg, &x, &y);
    if     (xyFlag==0){distPiID = (Axy*x-y+Bxy)/sqrt(Axy*Axy+1.*1.);}
    else if(xyFlag==1){distPiID = (Axy*y-x+Bxy)/sqrt(Axy*Axy+1.*1.);}
    double u=Dir.x(), v=Dir.y();
    double x0=Pos0.x(), y0=Pos0.y();
    double t = (u*(x-x0)+v*(y-y0))/(u*u+v*v);
    if (t>=0) {	
      if (fabs(distPiID)<40) {
	m_piid_seg = seg;
      }
    }
  }

#if 0    
  int nc = BGOCont_.size();
  for (int i=0; i<nc; i++) {
    Hodo1Hit *hitp = BGOCont_[i];
    //BGO_E_ += hitp->DeltaEBGO(); here!    
  }
  int ncPiV = PiVCont_.size();
  for (int i=0; i<ncPiV; i++) {
    Hodo2Hit *hitp = PiVCont_[i];
    PiV_E_ += hitp->DeltaE();
  }
#endif

  //TotalE_ = FiberTotal_E_ + BGO_E_;


#if 0    
  if (nc ==1 || nc>= 3) {
    double delta;
    if (checkProton(BGO_E_, &delta))
      Mass_ = 0.9382720;
    else if (checkPi(BGO_E_))
      Mass_ = 0.1395701;
  } else if (nc == 2) {
    //Hodo2Hit *hitp1 = BGOCont_[0];
    Hodo1Hit *hitp1 = BGOCont_[0];
    double BGO_E1 = hitp1->DeltaE();
    //Hodo2Hit *hitp2 = BGOCont_[1];
    Hodo1Hit *hitp2 = BGOCont_[1];
    double BGO_E2 = hitp2->DeltaE();

    bool flag0 = false,  flag1 = false, flag2 = false;
    double delta0, delta1, delta2;
    flag0 = checkProton(BGO_E_, &delta0);
    flag1 = checkProton(BGO_E1, &delta1);
    flag2 = checkProton(BGO_E2, &delta2);

    if (flag0) {
      Mass_ = 0.9382720;      
    } else if (flag1 && flag2) {
      Mass_ = 0.9382720;      
      if (fabs(delta1) < fabs(delta2))
	BGO_E_ = BGO_E1;
      else
	BGO_E_ = BGO_E2;
    } else if (flag1) {
      Mass_ = 0.9382720;      
      BGO_E_ = BGO_E1;
    } else if (flag2) {
      Mass_ = 0.9382720;      
      BGO_E_ = BGO_E2;
    }  else if (checkPi(BGO_E_)){
      Mass_ = 0.1395701;
    }
  } else if (nc == 0) {
    double delta;
    if (checkProton(BGO_E_, &delta))
      Mass_ = 0.9382720;
    else if (checkPi(BGO_E_) && NormalizedFiberTotal_E_<0.5)
      Mass_ = 0.1395701;

  }
#endif  

#if 0    
  double shiftE = 32.;
  if (PiV_E_>0.2 && Mass_ < 0. ) {
    double delta;
    if (checkProton(BGO_E_ - shiftE, &delta)) {
      Mass_ = 0.9382720;
      if (nc == 1 || nc >= 3) {
	BGO_E_ -= shiftE;
      } else if (nc == 2) {
	Hodo2Hit *hitp1 = BGOCont_[0];
	double BGO_E1 = hitp1->DeltaE();
	Hodo2Hit *hitp2 = BGOCont_[1];
	double BGO_E2 = hitp2->DeltaE();

	bool flag1 = false, flag2 = false;
	double delta1, delta2;
	flag1 = checkProton(BGO_E1, &delta1);
	flag2 = checkProton(BGO_E2, &delta2);

	if (flag1 && flag2) {
	  if (fabs(delta1)<fabs(delta2))
	    BGO_E_ = BGO_E1;
	  else
	    BGO_E_ = BGO_E2;
	} else if (flag1) {
	  BGO_E_ = BGO_E1;
	} else if (flag2) {
	  BGO_E_ = BGO_E2;
	} else {
	  BGO_E_ -= shiftE;
	}
      }
    }
  }
#endif  

#if 0    
  // assume proton stop before BGO and BGO was hit by pi
  if (Mass_ < 0. ) {
    double delta;
    double E=0.;
    if (checkProton(E, &delta)) {
      Mass_ = 0.9382720;
      BGO_E_ = 0;
    }
  }
#endif  

#if 0  
  if (NormalizedFiberTotal_E_>=0.5)
    Mass_ = 0.9382720;
  else if (NormalizedFiberTotal_E_>=0 && NormalizedFiberTotal_E_<0.5)
    Mass_ = 0.1395701;
#endif  

  return true;

}

#if 0
bool CFTParticle::checkProton( double BGO_E, double *delta)
{
  double cut1=0., cut2=0.;
  if (BGO_E >= 0 && BGO_E <= Eval[0]) {
    cut1 = Range1[0];
    cut2 = Range2[0];
  } else if (BGO_E >= Eval[NumOfPIDVal-1]) {
    cut1 = Range1[NumOfPIDVal-1];
    cut2 = Range2[NumOfPIDVal-1];
  } else {
    int index=-1;
    for (int i=0; i<NumOfPIDVal-1; i++) {
      if (BGO_E >= Eval[i] && BGO_E <= Eval[i+1]) {
	index = i;
	break;
      }
    }

    if (index != -1) {
      double v1 = BGO_E - Eval[index];
      double v2 = Eval[index+1] - BGO_E;
      cut1 = (v2*Range1[index]+v1*Range1[index+1])/(Eval[index+1]-Eval[index]);
      cut2 = (v2*Range2[index]+v1*Range2[index+1])/(Eval[index+1]-Eval[index]);
    }
  }

  *delta = BGO_E - (cut1+cut2)/2.;

  if (NormalizedFiberTotal_E_ >= cut1 && NormalizedFiberTotal_E_ <= cut2)
    return true;
  else 
    return false;

}

bool CFTParticle::checkPi( double BGO_E)
{
  double cut1=0.;
  if (BGO_E >= 0 && BGO_E <= Eval[0]) {
    cut1 = Range1[0];
  } else if (BGO_E >= Eval[NumOfPIDVal-1]) {
    cut1 = Range1[NumOfPIDVal-1];
  } else {
    int index=-1;
    for (int i=0; i<NumOfPIDVal-1; i++) {
      if (BGO_E >= Eval[i] && BGO_E <= Eval[i+1]) {
	index = i;
	break;
      }
    }

    if (index != -1) {
      double v1 = BGO_E - Eval[index];
      double v2 = Eval[index+1] - BGO_E;
      cut1 = (v2*Range1[index]+v1*Range1[index+1])/(Eval[index+1]-Eval[index]);
    }
  }

  if (NormalizedFiberTotal_E_ < cut1)
    return true;
  else 
    return false;

}
#endif

void CFTParticle::BGOPos(int seg, double *x, double *y) const
{
  int UnitNum = seg/(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);
  int SegInUnit = seg%(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);

  double theta = 22.5+(double)UnitNum*45.;
  double x0 = RadiusOfBGOSurface+BGO_Y/2;
  double y0 = (double)(SegInUnit-1)*BGO_X;
  double xc = 0.;
  double yc = 0.;

  if(seg%3==2){ // single BGO
    double n=(seg+1)/3;
    double angle = +22.5+45.*(n-1); // axis change
    xc = (120.+25./2.)*cos(angle*Deg2Rad);
    yc = (120.+25./2.)*sin(angle*Deg2Rad);    
#if 1 // new
  }else if(seg==0 || seg==1){
    xc = 100.0 + 25./2.;
    if(seg==0){yc = -30.0/2.;}
    else if(seg==1){yc = 30.0/2.;}
  }else if(seg==6 || seg==7){
    yc = 100.0 + 25./2.;
    if(seg==6){xc = 30.0/2.;}
    else if(seg==7){xc = -30.0/2.;}
  }else if(seg==12 || seg==13){
    xc = -100.0 - 25./2.;
    if(seg==12){yc = 30.0/2.;}
    else if(seg==13){yc = -30.0/2.;}  
  }else if(seg==18 || seg==19){
    yc = -100.0 -25./2.;
    if     (seg==18){xc = -30.0/2.;}
    else if(seg==19){xc = 30.0/2.;}
 }else if(seg==3 || seg==4){
    double angle = 45.;
    x0 = (100. + 25./2.)*cos(angle*Deg2Rad);
    y0 = (100. + 25./2.)*sin(angle*Deg2Rad);
    if(seg==4){
      xc = x0 - 30./2.*cos(-angle*Deg2Rad);
      yc = y0 - 30./2.*sin(-angle*Deg2Rad);	
    }else if(seg==3){
      xc = x0 + 30./2.*cos(-angle*Deg2Rad);
      yc = y0 + 30./2.*sin(-angle*Deg2Rad);	
    }
  }else if(seg==9 || seg==10){
    //double angle = -45.;
    double angle = 135.;
    x0 = (100 + 25./2.)*cos(angle*Deg2Rad);
    y0 = (100 + 25./2.)*sin(angle*Deg2Rad);      
    if(seg==10){
      xc = x0 + 30./2.*cos(-angle*Deg2Rad);
      yc = y0 + 30./2.*sin(-angle*Deg2Rad);	
    }else if(seg==9){
      xc = x0 - 30./2.*cos(-angle*Deg2Rad);
      yc = y0 - 30./2.*sin(-angle*Deg2Rad);	
    }
  }else if(seg==15 || seg==16){
    double angle = -135.;
    x0 = (100. + 25./2.)*cos(angle*Deg2Rad);
    y0 = (100. + 25./2.)*sin(angle*Deg2Rad);      
    if(seg==16){
      xc = x0 - 30./2.*cos(-angle*Deg2Rad);
      yc = y0 - 30./2.*sin(-angle*Deg2Rad);	
    }else if(seg==15){
      xc = x0 + 30./2.*cos(-angle*Deg2Rad);
      yc = y0 + 30./2.*sin(-angle*Deg2Rad);	
    }
  }else if(seg==21 || seg==22){
    //double angle = 135.;
    double angle = -45.;
    x0 = (100. + 25./2.)*cos(angle*Deg2Rad);
    y0 = (100. + 25./2.)*sin(angle*Deg2Rad);      
    if(seg==22){
      xc = x0 + 30./2.*cos(-angle*Deg2Rad);
      yc = y0 + 30./2.*sin(-angle*Deg2Rad);	
    }else if(seg==21){
      xc = x0 - 30./2.*cos(-angle*Deg2Rad);
      yc = y0 - 30./2.*sin(-angle*Deg2Rad);	
    }
  }
#endif
  *x = xc;
  *y = yc;
}

void CFTParticle::PiIDPos(int seg, double *x, double *y) const
{
  int UnitNum = seg/(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);
  int SegInUnit = seg%(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);
  double theta = 22.5+(double)UnitNum*45.;
  double x0 = RadiusOfBGOSurface+BGO_Y/2;
  double y0 = (double)(SegInUnit-1)*BGO_X;
  double xc = 0.;
  double yc = 0.;
  double w  = 30.,  t = 15.;//width, thickness
  double ww = 40., tt = 15.;//width, thickness for 45 deg.

  if(seg%4==3){ // single segment
    double n=(seg+1)/4;
    double angle = +22.5+45.*(n-1); // axis change
    xc = (164.1+tt/2.)*cos(angle*Deg2Rad);
    yc = (164.1+tt/2.)*sin(angle*Deg2Rad);    
  }else if(seg==0 || seg==1 || seg==2){      
    xc = 159.0 + t/2.;
    if     (seg==0){yc = -1.*w;}
    else if(seg==1){yc =  0.*w;}
    else if(seg==2){yc =  1.*w;}      

  }else if(seg==8 || seg==9 || seg==10){
    yc = 159.0 + t/2.;
    if     (seg==8) {xc = 1.*w;}
    else if(seg==9) {xc = 0.*w;}
    else if(seg==10){xc =-1.*w;}      

  }else if(seg==16 || seg==17 || seg==18){
    xc = -159.0 - t/2.;
    if     (seg==16){yc = 1.*w;}
    else if(seg==17){yc = 0.*w;}
    else if(seg==18){yc =-1.*w;}      

  }else if(seg==24 || seg==25 || seg==26){
    yc = -159.0 - t/2.;
    if     (seg==24) {xc =-1.*w;}
    else if(seg==25) {xc = 0.*w;}
    else if(seg==26) {xc = 1.*w;}      

  }else if(seg==4 || seg==5 || seg==6){ // Line
    double angle = 45.;
    x0 = (159.+t/2.)*cos(angle*Deg2Rad);
    y0 = (159.+t/2.)*sin(angle*Deg2Rad);
    if(seg==6){
      xc = x0 - w*cos(45.*Deg2Rad);
      yc = y0 + w*cos(45.*Deg2Rad);
    }if(seg==5){
      xc = x0 - 0*cos(45.*Deg2Rad);
      yc = y0 + 0*cos(45.*Deg2Rad);
    }else if(seg==4){
      xc = x0 + w*cos(45.*Deg2Rad); 
      yc = y0 - w*cos(45.*Deg2Rad); 
    }   

  }else if(seg==12 || seg==13 || seg==14){ // Line
      double angle = 135.;
      x0 = (159.+t/2.)*cos(angle*Deg2Rad);
      y0 = (159.+t/2.)*sin(angle*Deg2Rad);
      if(seg==12){
	xc = x0 + w*cos(45.*Deg2Rad);
	yc = y0 + w*cos(45.*Deg2Rad);
      }else if(seg==13){
	xc = x0 + 0*cos(45.*Deg2Rad);
	yc = y0 + 0*cos(45.*Deg2Rad);
      }else if(seg==14){
	xc = x0 - w*cos(45.*Deg2Rad);
	yc = y0 - w*cos(45.*Deg2Rad);
      }

  }else if(seg==20 || seg==21 || seg==22){ // Line
    double angle = -135.;
    x0 = (159. +t/2.)*cos(angle*Deg2Rad);
    y0 = (159. +t/2.)*sin(angle*Deg2Rad);
    if(seg==22){
      xc = x0 + 1.*w*cos(45.*Deg2Rad); 
      yc = y0 - 1.*w*cos(45.*Deg2Rad); 
    }else if(seg==21){
      xc = x0 + 0.*w*cos(45.*Deg2Rad); 
      yc = y0 - 0.*w*cos(45.*Deg2Rad); 
    }else if(seg==20){
      xc = x0 - 1.*w*cos(45.*Deg2Rad);
      yc = y0 + 1.*w*cos(45.*Deg2Rad);
    }

  }else if(seg==28 || seg==29 || seg==30){ // Line
    double angle = -45.;
    x0 = (159.+t/2.)*cos(angle*Deg2Rad);
    y0 = (159.+t/2.)*sin(angle*Deg2Rad);
    if(seg==28){
      xc = x0 - 1.*w*cos(45.*Deg2Rad);
      yc = y0 - 1.*w*cos(45.*Deg2Rad);
    }else if(seg==29){
      xc = x0 - 0.*w*cos(45.*Deg2Rad);
      yc = y0 - 0.*w*cos(45.*Deg2Rad);
    }else if(seg==30){
      xc = x0 + 1.*w*cos(45.*Deg2Rad);
      yc = y0 + 1.*w*cos(45.*Deg2Rad);
    }
  }

  *x = xc;
  *y = yc;
}

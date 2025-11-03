/*
  CFTLocalTrack.hh

  2012/1/24
*/

#ifndef CFTLocalTrack_h
#define CFTLocalTrack_h 1

#include <vector>
#include <functional>
#include "ThreeVector.hh"
#include "CFTFiberCluster.hh"
#include "MathTools.hh"

class CFTFiberCluster;
class DCAnalyzer;

class CFTLocalTrack
{
public:
  explicit CFTLocalTrack();
  ~CFTLocalTrack();
private:
  CFTLocalTrack( const CFTLocalTrack & );
  CFTLocalTrack & operator = ( const CFTLocalTrack & );

private:
  std::vector <CFTFiberCluster *> hitArray;
  std::vector <CFTFiberCluster *> hitArrayU;
  ThreeVector CFTVtx_;

  int    xyFitFlag_;
  double Axy_, Bxy_;
  double Az_, Bz_;
  double chisqrXY_;
  double chisqrZ_;
  double vtx_z_;
  double chisqr_dE_;
  int    nth_ztracking_;
  // three dementional track
  int    zTrackFlag_;
  double u0_, v0_, x0_, y0_;
  double ur0_, vr0_, f0_, z0_;
  ThreeVector Dir_, Pos0_;
  ThreeVector OffsetPos0_;
  double theta_;

  // dE informaiton
  double TotalDEHiGain_,TotalDELowGain_;
  double MaxDEHiGain_,MaxDELowGain_;
  double MaxDEHiGain_phi_,MaxDELowGain_phi_;
  double MaxDEHiGain_uv_,MaxDELowGain_uv_;
  double TotalDELowGain_phi_;
  double TotalDELowGain_uv_;
  double pathlength_;

public:
  void AddHit( CFTFiberCluster *hitp ) { hitArray.push_back( hitp ); }
  void AddCFTVtx(ThreeVector vtx) {CFTVtx_=vtx;}
  bool DoFitXY( void );
  bool DoFitXY_2nd( void );
  void AddHitU( CFTFiberCluster *hitp ) { hitArrayU.push_back( hitp ); }
  bool DoFitZTrack();
  bool DoFitZTrack_16layer();
  bool DoFitXY_wVtx( void );
  bool CheckPhi( void );
  bool CheckPhi_1st( void );
  bool CheckEdep( void );
  bool DoFitZTrack_wVtx();
  bool DoFit( void );
  bool FindZTrack();
  bool SetCalculatedValue();
  bool SetCorPhi();
  bool SetCalculatedValue_16layer();

  std::size_t GetNHit( void ) const {return hitArray.size(); }
  std::size_t GetNHitU( void ) const {return hitArrayU.size(); }
  std::size_t GetNHitUV( void ) const {return hitArrayU.size(); }
  CFTFiberCluster * GetHit( std::size_t nth ) const;
  CFTFiberCluster * GetHitOfLayerNumber( int lnum ) const;
  CFTFiberCluster * GetHitU( std::size_t nth ) const;
  CFTFiberCluster * GetHitUV( std::size_t nth ) const { return GetHitU(nth); };
  int    GetFirstLayerPhi( void ) const;
  int    GetFirstLayerUV( void ) const;
  double GetChiSquareXY( void ) const { return chisqrXY_; }
  double GetChiSquareZ( void ) const { return chisqrZ_; }
  double GetChiSquareDE( void ) const { return chisqr_dE_; }
  /*
  void SetAv( double Av) { Av_=Av; }
  void SetAx( double Ax) { Ax_=Ax; }
  void SetAu( double Au) { Au_=Au; }
  void SetChiv( double Chiv) { Chiv_=Chiv; }
  void SetChix( double Chix) { Chix_=Chix; }
  void SetChiu( double Chiu) { Chiu_=Chiu; }
  */
  double GetX0( void ) const { return x0_; }
  double GetY0( void ) const { return y0_; }
  double GetU0( void ) const { return u0_; }
  double GetV0( void ) const { return v0_; }


  double GetChiSquare( void ) const { return chisqr_; }
  double GetX( double z ) const { return x0_+u0_*z; } 
  double GetY( double z ) const { return y0_+v0_*z; } 
  int    GetXYFitFlag( void ) const { return xyFitFlag_; }
  int    GetCFTxyFlag( void ) const { return xyFitFlag_; }
  double GetAxy( void ) const { return Axy_; }
  double GetBxy( void ) const { return Bxy_; }

  ThreeVector GetPos0( void ) const { return Pos0_ + OffsetPos0_; }
  ThreeVector GetDir( void )  const { return Dir_; }
  double GetThetaCFT( void )  const { return acos(Dir_.z()/Dir_.Mag())*TMath::RadToDeg(); }
  void SetThetaCFT( void )  { theta_ =  acos(Dir_.z()/Dir_.Mag())*TMath::RadToDeg(); }
  void ResetThetaCFT( void )  { theta_ = -999.; }
  double GetPhiCFT( void )  const { return Dir_.Phi()*TMath::RadToDeg(); }
  
  double GetPhiFromR( double r ) const { return f0_+ur0_*r; } 
  double GetZFromR( double r ) const { return z0_+vr0_*r; } 

  int   GetZTrackFlag( void ) const { return zTrackFlag_; }
  int   GetCFTzFlag( void ) const { return zTrackFlag_; }
  double GetAz( void ) const { return Az_; }
  double GetBz( void ) const { return Bz_; }
  double GetVtxZ( void ) const { return vtx_z_; }

  bool   GetPathLengthInFiber(double r, double phi, double *len);
  bool   GetCrossPointR(double r, double *phi1, double *phi2);
  bool   GetCrossPointR(double r, double *phi);
  double calcPhi(double x, double y);
  double CalculateZpos(double phi, CFTFiberCluster *cl, double mean_z = -999999.);
  double CalculateZpos2(double phi, CFTFiberCluster *cl, double mean_z = -99999.);
  double CalculateUVpos(double phi, CFTFiberCluster *cl);

  bool   CalcNormalizedDE( void ); 

  double    TotalDEHiGain() const { return TotalDEHiGain_;};
  double    MaxDEHiGain() const { return MaxDEHiGain_;};
  double    TotalDELowGain() const { return TotalDELowGain_;};
  double    GetTotalSumdE() const { return TotalDELowGain_;};
  void      SetTotalSumdE(double dE) { TotalDELowGain_=dE;};
  double    MaxDELowGain() const { return MaxDELowGain_;};
  double    GetTotalMaxdE() const { return MaxDELowGain_;};

  double    GetTotalMaxdEphi() const { return MaxDELowGain_phi_;};
  double    GetTotalMaxdEuv() const { return MaxDELowGain_uv_;};
  double    GetTotalSumdEphi() const { return TotalDELowGain_phi_;};
  double    GetTotalSumdEuv() const { return TotalDELowGain_uv_;};
  
  double    GetTotalPathLength() const  { return pathlength_;};
  double    NormalizedTotalDEHiGain (void) const {return TotalDEHiGain_/pathlength_;};
  double    NormalizedTotalDELowGain (void) const {return TotalDELowGain_/pathlength_;};
  double    NormalizedMaxDEHiGain (void) const {return MaxDEHiGain_/pathlength_;};
  double    NormalizedMaxDELowGain (void) const {return MaxDELowGain_/pathlength_;};

  bool GetStatus( void ) const { return status_; } 
  bool GoodForTracking( void ) const { return gftstatus_; }
  bool GoodForTracking( bool status )
  { bool ret=gftstatus_; gftstatus_=status; return ret; } 
  //bool ReCalc( bool ApplyRecursively=false );  

private:
  bool status_;
  //double x0_, y0_, u0_, v0_;
  double a_,b_;
  double chisqr_;
  bool gftstatus_;
};

#if 0
struct CFTLTrackComp 
  : public std::binary_function <CFTLocalTrack *, CFTLocalTrack *, bool>
{
  bool operator()( const CFTLocalTrack * const p1, 
		   const CFTLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquareXY(),chi2=p2->GetChiSquareXY();
    if( (n1>n2+1) ){
      return true;
    }
    else if( (n2>n1+1)  ){
      return false;
    }
    else{
      return (chi1<=chi2);
    }
  }
};
#endif

struct CFTLTrackComp 
  : public std::binary_function <CFTLocalTrack *, CFTLocalTrack *, bool>
{
  bool operator()( const CFTLocalTrack * const p1, 
		   const CFTLocalTrack * const p2 ) const
  {
    double thr_chisqr=500.;

    int n1=p1->GetNHit()+p1->GetNHitU();
    int n2=p2->GetNHit()+p2->GetNHitU();

    double chi1=p1->GetChiSquareXY()+p1->GetChiSquareZ();
    double chi2=p2->GetChiSquareXY()+p2->GetChiSquareZ();
    /*
    if (n1>=3 && n2 >=3) {
      std::cout << "1 : n1 " << n1 << "chi_XY " << p1->GetChiSquareXY() << ", chi_Z " << p1->GetChiSquareZ() << std::endl;
      std::cout << "2 : n2 " << n2 << "chi_XY " << p2->GetChiSquareXY() << ", chi_Z " << p2->GetChiSquareZ() << std::endl;
    }
    */

    if( n1>n2 && chi1<thr_chisqr)
      return true;
    else if( n2>n1 && chi2<thr_chisqr)
      return false;
    else
      return (chi1<=chi2);
    
  }
};

#endif

// -*- C++ -*-

#include "KuramaFieldMap.hh"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <std_ostream.hh>

//#include "ConfMan.hh"
#include "FuncName.hh"

//namespace
//{
//const auto& gConf = ConfMan::GetInstance();
//const auto& valueNMR  = ConfMan::Get<Double_t>("FLDNMR");
//const auto& valueCalc = ConfMan::Get<Double_t>("FLDCALC");
//}

//_____________________________________________________________________________
KuramaFieldMap::KuramaFieldMap(const TString& file_name)
  : m_is_ready(false),
    m_file_name(file_name),
    B(),
    Nx(0),
    Ny(0),
    Nz(0),
    X0(),
    Y0(),
    Z0(),
    dX(),
    dY(),
    dZ()
{
}

//_____________________________________________________________________________
KuramaFieldMap::KuramaFieldMap(const TString& file_name, const Double_t measure, const Double_t calc)
  : m_is_ready(false),
    m_file_name(file_name),
    B(),
    Nx(0),
    Ny(0),
    Nz(0),
    X0(),
    Y0(),
    Z0(),
    dX(),
    dY(),
    dZ(),
    valueMeasure(measure),
    valueCalc(calc)
{
}


//_____________________________________________________________________________
KuramaFieldMap::~KuramaFieldMap()
{
  ClearField();
}

//_____________________________________________________________________________
Bool_t
KuramaFieldMap::Initialize()
{
  std::ifstream ifs(m_file_name);
  if(!ifs.is_open()){
    hddaq::cerr << FUNC_NAME << " file open fail : "
		<< m_file_name << std::endl;
    return false;
  }

  if(m_is_ready){
    hddaq::cerr << FUNC_NAME
		<< " already initialied" << std::endl;
    return false;
  }

  ClearField();

  if(!(ifs >> Nx >> Ny >> Nz >> X0 >> Y0 >> Z0 >> dX >> dY >> dZ)){
    hddaq::cerr << FUNC_NAME << " invalid format" << std::endl;
    return false;
  }

  if(Nx<0 || Ny<0 || Nz<0){
    hddaq::cerr << FUNC_NAME << " Nx, Ny, Nz must be positive" << std::endl;
    return false;
  }

  B.resize(Nx);
  for(Int_t ix=0; ix<Nx; ++ix){
    B[ix].resize(Ny);
    for(Int_t iy=0; iy<Ny; ++iy){
      B[ix][iy].resize(Nz);
    }
  }

  if(valueCalc==0. || !std::isfinite(valueCalc) ||
     valueMeasure==0.  || !std::isfinite(valueMeasure) ){
    hddaq::cout << FUNC_NAME << " KuramaField is zero : "
                << " Calc = " << valueCalc
                << " Measure = " << valueMeasure << std::endl
                << " -> skip reading fieldmap" << std::endl;
    return true;
  }
  const Double_t factor = valueMeasure/valueCalc;

  Double_t x, y, z, bx, by, bz;

  hddaq::cout << " reading fieldmap " << std::flush;

  Int_t line = 0;
  while(ifs.good()){
    if(line++%1000000==0)
      hddaq::cout << "." << std::flush;
    ifs >> x >> y >> z >> bx >> by >> bz;
    Int_t ix = Int_t((x-X0+0.1*dX)/dX);
    Int_t iy = Int_t((y-Y0+0.1*dY)/dY);
    Int_t iz = Int_t((z-Z0+0.1*dZ)/dZ);
    if(ix>=0 && ix<Nx && iy>=0 && iy<Ny && iz>=0 && iz<Nz){
      B[ix][iy][iz].x = bx*factor;
      B[ix][iy][iz].y = by*factor;
      B[ix][iy][iz].z = bz*factor;
    }
  }

  hddaq::cout << " done" << std::endl;
  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
KuramaFieldMap::GetFieldValue(const Double_t pointCM[3],
                              Double_t* BfieldTesla) const
{
  Double_t xt = pointCM[0];
  Double_t yt = pointCM[1];
  Double_t zt = pointCM[2];

  Int_t ix1, ix2, iy1, iy2, iz1, iz2;
  ix1 = Int_t((xt-X0)/dX);
  iy1 = Int_t((yt-Y0)/dY);
  iz1 = Int_t((zt-Z0)/dZ);

  Double_t wx1, wx2, wy1, wy2, wz1, wz2;
  if(ix1<0) { ix1=ix2=0; wx1=1.; wx2=0.; }
  else if(ix1>=Nx-1) { ix1=ix2=Nx-1; wx1=1.; wx2=0.; }
  else { ix2=ix1+1; wx1=(X0+dX*ix2-xt)/dX; wx2=1.-wx1; }

  if(iy1<0) { iy1=iy2=0; wy1=1.; wy2=0.; }
  else if(iy1>=Ny-1) { iy1=iy2=Ny-1; wy1=1.; wy2=0.; }
  else { iy2=iy1+1; wy1=(Y0+dY*iy2-yt)/dY; wy2=1.-wy1; }

  if(iz1<0) { iz1=iz2=0; wz1=0.; wz2=0.; }
  else if(iz1>=Nz-1) { iz1=iz2=Nz-1; wz1=1.; wz2=0.; }
  else { iz2=iz1+1; wz1=(Z0+dZ*iz2-zt)/dZ; wz2=1.-wz1; }

  Double_t bx1 = wx1*wy1*B[ix1][iy1][iz1].x + wx1*wy2*B[ix1][iy2][iz1].x
    + wx2*wy1*B[ix2][iy1][iz1].x + wx2*wy2*B[ix2][iy2][iz1].x;
  Double_t bx2 = wx1*wy1*B[ix1][iy1][iz2].x + wx1*wy2*B[ix1][iy2][iz2].x
    + wx2*wy1*B[ix2][iy1][iz2].x + wx2*wy2*B[ix2][iy2][iz2].x;
  Double_t bx  = wz1*bx1 + wz2*bx2;

  Double_t by1 = wx1*wy1*B[ix1][iy1][iz1].y + wx1*wy2*B[ix1][iy2][iz1].y
    + wx2*wy1*B[ix2][iy1][iz1].y + wx2*wy2*B[ix2][iy2][iz1].y;
  Double_t by2 = wx1*wy1*B[ix1][iy1][iz2].y + wx1*wy2*B[ix1][iy2][iz2].y
    + wx2*wy1*B[ix2][iy1][iz2].y + wx2*wy2*B[ix2][iy2][iz2].y;
  Double_t by  = wz1*by1 + wz2*by2;

  Double_t bz1 = wx1*wy1*B[ix1][iy1][iz1].z + wx1*wy2*B[ix1][iy2][iz1].z
    + wx2*wy1*B[ix2][iy1][iz1].z + wx2*wy2*B[ix2][iy2][iz1].z;
  Double_t bz2 = wx1*wy1*B[ix1][iy1][iz2].z + wx1*wy2*B[ix1][iy2][iz2].z
    + wx2*wy1*B[ix2][iy1][iz2].z + wx2*wy2*B[ix2][iy2][iz2].z;
  Double_t bz  = wz1*bz1 + wz2*bz2;

  //Default
  BfieldTesla[0] = bx;
  BfieldTesla[1] = by;
  BfieldTesla[2] = bz;
  return true;
}

//_____________________________________________________________________________
void
KuramaFieldMap::ClearField()
{
  for(Int_t ix=0; ix<Nx; ++ix){
    for(Int_t iy=0; iy<Ny; ++iy){
      B[ix][iy].clear();
    }
    B[ix].clear();
  }
  B.clear();
}

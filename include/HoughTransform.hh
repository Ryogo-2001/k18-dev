// -*- C++ -*-

#ifndef HOUGH_TRANSFORM_HH
#define HOUGH_TRANSFORM_HH

#include <vector>

#include <TString.h>
#include <TVector3.h>
#include "DCAnalyzer.hh"

class TPCCluster;

namespace tpc
{

//HS-OFF
//Functions for track searching
Bool_t HoughTransformLineXZ(std::vector<TVector3> gHitPos,
			    Int_t *MaxBin, Double_t *LinearPar,
			    Int_t MinNumOfHits=8);
Bool_t HoughTransformLineXZ(const std::vector<TPCClusterContainer>& ClCont,
			    Int_t *MaxBin, Double_t *LinearPar,
			    Int_t MinNumOfHits=8);
void HoughTransformLineYZ(std::vector<TVector3> gHitPos,
			  Int_t *MaxBin, Double_t *LinearPar,
			  Double_t MaxHoughWindowY);
void HoughTransformLineYZ(const std::vector<TPCClusterContainer>& ClCont,
			  Int_t *MaxBin, Double_t *LinearPar,
			  Double_t MaxHoughWindowY);
void HoughTransformLineYX(std::vector<TVector3> gHitPos,
			  Int_t *MaxBin, Double_t *LinearPar,
			  Double_t MaxHoughWindowY);
void HoughTransformLineYX(const std::vector<TPCClusterContainer>& ClCont,
			  Int_t *MaxBin, Double_t *LinearPar,
			  Double_t MaxHoughWindowY);

//HS-ON
//Functions for the searching
Bool_t HoughTransformCircleXZ(std::vector<TVector3> gHitPos,
			      Int_t *MaxBin, Double_t *HelixPar,
			      Int_t MinNumOfHits=8);
Bool_t HoughTransformCircleXZ(const std::vector<TPCClusterContainer>& ClCont,
			      Int_t *MaxBin, Double_t *HelixPar,
			      Int_t MinNumOfHits=8);
void HoughTransformLineYTheta(std::vector<TVector3> gHitPos,
			      Int_t *MaxBin, Double_t *HelixPar,
			      Double_t MaxHoughWindowY);
void HoughTransformLineYTheta(const std::vector<TPCClusterContainer>& ClCont,
			      Int_t *MaxBin, Double_t *HelixPar,
			      Double_t MaxHoughWindowY);

}

#endif

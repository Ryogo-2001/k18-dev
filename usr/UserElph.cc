// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iostream>
#include <sstream>

#include <UnpackerManager.hh>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "HodoRawHit.hh"

#define TimeCut    1 // in cluster analysis

namespace
{
using namespace root;
const auto qnan = TMath::QuietNaN();
auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
auto& gRM = RMAnalyzer::GetInstance();
auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
class UserElph : public VEvent
{
private:
  RawData*      rawData;

public:
  UserElph();
  ~UserElph();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal(); };

//_____________________________________________________________________________
inline const TString&
UserElph::ClassName()
{
  static TString s_name("UserElph");
  return s_name;
}

//_____________________________________________________________________________
UserElph::UserElph()
  : VEvent(),
    rawData(new RawData)
{
}

//_____________________________________________________________________________
UserElph::~UserElph()
{
  if(rawData) delete rawData;
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;
/*
  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  Int_t bacnhits;
  Int_t bachitpat[MaxHits];
  Double_t baca[NumOfSegBAC];
  Double_t bact[NumOfSegBAC][MaxDepth];
*/
  Int_t trigpat[NumOfSegTFlagElph];
  Int_t trigflag[NumOfSegTFlagElph];
//  Double_t v792a[NumOfChV792];// for old CMAP (elph_20220908)
  Double_t E72BACa[NumOfSegE72BAC];
  Double_t E72BACt[NumOfSegE72BAC][MaxDepth];  
  Double_t E72BACSUMa[NumOfSegE72BACSUM];
  Double_t E72BACSUMt[NumOfSegE72BACSUM][MaxDepth];  
  Double_t KVCa[NumOfSegKVC];
  Double_t KVCt[NumOfSegKVC][MaxDepth];
  Double_t KVCSUMa[NumOfSegKVCSUM];
  Double_t KVCSUMt[NumOfSegKVCSUM][MaxDepth];
  Double_t E90SACa[NumOfSegE90SAC];
  Double_t E90SACt[NumOfSegE90SAC][MaxDepth];
  Double_t E90SACSUMa[NumOfSegE90SACSUM];
  Double_t E90SACSUMt[NumOfSegE90SACSUM][MaxDepth];
  Double_t T4a[NumOfSegT4];
  Double_t T4t[NumOfSegT4][MaxDepth];
  Double_t T5a[NumOfSegT5];
  Double_t T5t[NumOfSegT5][MaxDepth];
  Double_t T6a[NumOfSegT6];
  Double_t T6t[NumOfSegT6][MaxDepth];
  Double_t T7a[NumOfSegT7];
  Double_t T7t[NumOfSegT7][MaxDepth];

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum      = 0;
 /* 
  bacnhits   = 0;

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<MaxHits; ++it){
    bachitpat[it]   = -1;
  }

  for(Int_t it=0; it<NumOfSegBAC; it++){
    baca[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      bact[it][m] = qnan;
    }
  }
 */
  for(Int_t it=0; it<NumOfSegTFlagElph; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  // for(Int_t it=0; it<NumOfChV792; it++){
  //   v792a[it] = qnan;
  // }
  for(Int_t it=0; it<NumOfSegE72BAC; it++){
    E72BACa[it] = qnan;
  }
  for(Int_t it=0; it<NumOfSegE72BACSUM; it++){
    E72BACSUMa[it] = qnan;
  }
  for(Int_t it=0; it<NumOfSegKVC; it++){
    KVCa[it] = qnan;
  }
  for(Int_t it=0; it<NumOfSegKVCSUM; it++){
    KVCSUMa[it] = qnan;
  }
  for(Int_t it=0; it<NumOfSegE90SAC; it++){
    E90SACa[it] = qnan;
  }
  for(Int_t it=0; it<NumOfSegE90SACSUM; it++){
    E90SACSUMa[it] = qnan;
  }
  for(Int_t it=0; it<NumOfSegT4; it++){
    T4a[it] = qnan;
  }
  for(Int_t it=0; it<NumOfSegT5; it++){
    T5a[it] = qnan;
  } 
  for(Int_t it=0; it<NumOfSegT6; it++){
    T6a[it] = qnan;
  } 
  for(Int_t it=0; it<NumOfSegT7; it++){
    T7a[it] = qnan;
  } 

  for(Int_t it=0; it<NumOfSegE72BAC; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      E72BACt[it][m]  = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegE72BACSUM; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      E72BACSUMt[it][m]  = qnan;
    }
  
  }
  for(Int_t it=0; it<NumOfSegKVC; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      KVCt[it][m]  = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegKVCSUM; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      KVCSUMt[it][m]  = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegE90SAC; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      E90SACt[it][m]  = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegE90SACSUM; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      E90SACSUMt[it][m]  = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegT4; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      T4t[it][m]  = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegT5; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      T5t[it][m]  = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegT6; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      T6t[it][m]  = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegT7; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      T7t[it][m]  = qnan;
    }
  }

}


//_____________________________________________________________________________
namespace root
{
Event  event;
TH1*   h[MaxHist];
TTree* tree;
enum eDetHid {
  V792Hid    = 10000,
  BACHid    = 20000
};
}

//_____________________________________________________________________________
Bool_t
UserElph::ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserElph::ProcessingNormal()
{
  //static const auto MinTdcBH1 = gUser.GetParameter("TdcBH1", 0);
  //static const auto MaxTdcBH1 = gUser.GetParameter("TdcBH1", 1);

  rawData->DecodeHits();

  gRM.Decode();

  // event.evnum = gRM.EventNumber();
  event.evnum = gUnpacker.get_event_number();

  HF1(1, 0);

  //**************************************************************************
  //****************** RawData

  ///// V792
  // {
  //   //    Int_t tof_nhits = 0;
  //   const HodoRHitContainer& cont = rawData->GetV792RawHC();
  //   Int_t nh = cont.size();
  //   HF1(V792Hid, nh);
  //   Int_t nh1 = 0, nh2 = 0;
  //   for(Int_t i=0; i<nh; ++i){
  //     HodoRawHit *hit = cont[i];
  //     Int_t seg = hit->SegmentId()+1;
  //     HF1(V792Hid+1, seg-0.5);
  //     Int_t A = hit->GetAdcUp();
  //     HF1(V792Hid+100*seg+1, A);
  //     event.v792a[seg-1] = A;
  //   }
  // }

  const Int_t NumOfSegTFlagElph = 5;
  std::bitset<NumOfSegTFlagElph> trigger_flag;
  for(const auto& hit: rawData->GetTrigRawHC()){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc1();
    if(tdc > 0){
      event.trigpat[trigger_flag.count()] = seg;
      event.trigflag[seg] = tdc;
    }
  }

  ///// E72BAC
  {
    //    Int_t tof_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetE72BACRawHC();
    Int_t nh = cont.size();
    // HF1(V792Hid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      //HF1(V792Hid+1, seg-0.5);
      Int_t A = hit->GetAdcUp();
      //HF1(V792Hid+100*seg+1, A);
      event.E72BACa[seg-1] = A;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.E72BACt[seg-1][m] = T;
      }
    }
  }


  ///// E72BACSUM
  {
    //    Int_t tof_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetE72BACSUMRawHC();
    Int_t nh = cont.size();
    // HF1(V792Hid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      //HF1(V792Hid+1, seg-0.5);
      Int_t A = hit->GetAdcUp();
      //HF1(V792Hid+100*seg+1, A);
      event.E72BACSUMa[seg-1] = A;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.E72BACSUMt[seg-1][m] = T;
      }
    }
  }
  ///// KVC
  {
    //    Int_t tof_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetKVCRawHC();
    Int_t nh = cont.size();
    // HF1(V792Hid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      //HF1(V792Hid+1, seg-0.5);
      Int_t A = hit->GetAdcUp();
      //HF1(V792Hid+100*seg+1, A);
      event.KVCa[seg-1] = A;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.KVCt[seg-1][m] = T;
      }
    }
  }


  ///// KVCSUM
  {
    //    Int_t tof_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetKVCSUMRawHC();
    Int_t nh = cont.size();
    // HF1(V792Hid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      //HF1(V792Hid+1, seg-0.5);
      Int_t A = hit->GetAdcUp();
      //HF1(V792Hid+100*seg+1, A);
      event.KVCSUMa[seg-1] = A;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.KVCSUMt[seg-1][m] = T;
      }
    }
  }
  ///// E90SAC
  {
    //    Int_t tof_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetE90SACRawHC();
    Int_t nh = cont.size();
    // HF1(V792Hid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      //HF1(V792Hid+1, seg-0.5);
      Int_t A = hit->GetAdcUp();
      //HF1(V792Hid+100*seg+1, A);
      event.E90SACa[seg-1] = A;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.E90SACt[seg-1][m] = T;
      }
    }
  }


  ///// E90SACSUM
  {
    //    Int_t tof_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetE90SACSUMRawHC();
    Int_t nh = cont.size();
    // HF1(V792Hid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      //HF1(V792Hid+1, seg-0.5);
      Int_t A = hit->GetAdcUp();
      //HF1(V792Hid+100*seg+1, A);
      event.E90SACSUMa[seg-1] = A;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.E90SACSUMt[seg-1][m] = T;
      }
    }
  }
  ///// T4
  {
    //    Int_t tof_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetT4RawHC();
    Int_t nh = cont.size();
    // HF1(V792Hid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      //HF1(V792Hid+1, seg-0.5);
      Int_t A = hit->GetAdcUp();
      //HF1(V792Hid+100*seg+1, A);
      event.T4a[seg-1] = A;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.T4t[seg-1][m] = T;
      }
    }
  }

  ///// T5
  {
    //    Int_t tof_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetT5RawHC();
    Int_t nh = cont.size();
    // HF1(V792Hid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      //HF1(V792Hid+1, seg-0.5);
      Int_t A = hit->GetAdcUp();
      //HF1(V792Hid+100*seg+1, A);
      event.T5a[seg-1] = A;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.T5t[seg-1][m] = T;
      }
    }
  }

  ///// T6
  {
    //    Int_t tof_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetT6RawHC();
    Int_t nh = cont.size();
    // HF1(V792Hid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      //HF1(V792Hid+1, seg-0.5);
      Int_t A = hit->GetAdcUp();
      //HF1(V792Hid+100*seg+1, A);
      event.T6a[seg-1] = A;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.T6t[seg-1][m] = T;
      }
    }
  }

  ///// T7
  {
    //    Int_t tof_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetT7RawHC();
    Int_t nh = cont.size();
    // HF1(V792Hid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      //HF1(V792Hid+1, seg-0.5);
      Int_t A = hit->GetAdcUp();
      //HF1(V792Hid+100*seg+1, A);
      event.T7a[seg-1] = A;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.T7t[seg-1][m] = T;
      }
    }
  }
  
  
  return true;
}

//_____________________________________________________________________________
Bool_t
UserElph::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserElph;
}

//_____________________________________________________________________________
namespace
{
const Int_t    NbinAdc = 4096;
const Double_t MinAdc  =    0.;
const Double_t MaxAdc  = 4096.;

const Int_t    NbinTdc = 4096;
const Double_t MinTdc  =    0.;
const Double_t MaxTdc  = 4096.;

const Int_t    NbinTdcHr = 1e6/10;
const Double_t MinTdcHr  =  0.;
const Double_t MaxTdcHr  = 1e6;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  HB1(1, "Status", 20, 0., 20.);

  // V792

  // for(Int_t i=1; i<=NumOfChV792; ++i){
  //   TString title1 = Form("V792-%d Adc", i);
  //   HB1(V792Hid +100*i +1, title1, NbinAdc, MinAdc, MaxAdc);
  // }


  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
    //Trig
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTFlagElph));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTFlagElph));
  //V792
  //  tree->Branch("v792a",       event.v792a,      Form("v792a[%d]/D", NumOfChV792));
  //E72BACa
  tree->Branch("E72BACa", event.E72BACa,  Form("E72BACa[%d]/D", NumOfSegE72BAC));
  tree->Branch("E72BACt", event.E72BACt,  Form("E72BACt[%d][%d]/D", NumOfSegE72BAC, MaxDepth));
  tree->Branch("E72BACSUMa", event.E72BACSUMa,  Form("E72BACSUMa[%d]/D", NumOfSegE72BACSUM));
  tree->Branch("E72BACSUMt", event.E72BACSUMt,  Form("E72BACSUMt[%d][%d]/D", NumOfSegE72BACSUM, MaxDepth));
  //KVCa
  tree->Branch("KVCa", event.KVCa,  Form("KVCa[%d]/D", NumOfSegKVC));
  tree->Branch("KVCt", event.KVCt,  Form("KVCt[%d][%d]/D", NumOfSegKVC, MaxDepth));
  tree->Branch("KVCSUMa", event.KVCSUMa,  Form("KVCSUMa[%d]/D", NumOfSegKVCSUM));
  tree->Branch("KVCSUMt", event.KVCSUMt,  Form("KVCSUMt[%d][%d]/D", NumOfSegKVCSUM, MaxDepth));
  //E90SACa
  tree->Branch("E90SACa", event.E90SACa,  Form("E90SACa[%d]/D", NumOfSegE90SAC));
  tree->Branch("E90SACt", event.E90SACt,  Form("E90SACt[%d][%d]/D", NumOfSegE90SAC, MaxDepth));
  tree->Branch("E90SACSUMa", event.E90SACSUMa,  Form("E90SACSUMa[%d]/D", NumOfSegE90SACSUM));
  tree->Branch("E90SACSUMt", event.E90SACSUMt,  Form("E90SACSUMt[%d][%d]/D", NumOfSegE90SACSUM, MaxDepth));
  //T4
  tree->Branch("T4a", event.T4a,  Form("T4a[%d]/D", NumOfSegT4));
  tree->Branch("T4t", event.T4t,  Form("T4t[%d][%d]/D", NumOfSegT4, MaxDepth));
  //T5
  tree->Branch("T5a", event.T5a,  Form("T5a[%d]/D", NumOfSegT5));
  tree->Branch("T5t", event.T5t,  Form("T5t[%d][%d]/D", NumOfSegT5, MaxDepth));
  //T6
  tree->Branch("T6a", event.T6a,  Form("T6a[%d]/D", NumOfSegT6));
  tree->Branch("T6t", event.T6t,  Form("T6t[%d][%d]/D", NumOfSegT6, MaxDepth));
  //T7
  tree->Branch("T7a", event.T7a,  Form("T7a[%d]/D", NumOfSegT7));
  tree->Branch("T7t", event.T7t,  Form("T7t[%d][%d]/D", NumOfSegT7, MaxDepth));



  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}

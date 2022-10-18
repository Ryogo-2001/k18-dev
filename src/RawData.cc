// -*- C++ -*-

#include "RawData.hh"

#include <algorithm>
#include <iostream>
#include <string>

#include <TF1.h>

#include <std_ostream.hh>
#include <UnpackerManager.hh>

#include "ConfMan.hh"
#include "DCRawHit.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "DetectorID.hh"
#include "Exception.hh"
#include "FuncName.hh"
#include "HodoRawHit.hh"
#include "MathTools.hh"
#include "TPCPadHelper.hh"
#include "TPCRawHit.hh"
#include "UserParamMan.hh"

#define OscillationCut 0
#define DebugEvDisp    0

#if DebugEvDisp
#include <TApplication.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
namespace { TApplication app("DebugApp", nullptr, nullptr); }
#endif

namespace
{
using namespace hddaq::unpacker;
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser     = UserParamMan::GetInstance();
enum EUorD { kOneSide=1, kBothSide=2 };
enum EHodoDataType { kHodoAdc, kHodoLeading, kHodoTrailing,
  kHodoOverflow, kHodoNDataType };
#if OscillationCut
const Int_t  MaxMultiHitDC  = 16;
#endif

///// for CorrectBaselineTPC()
TH1D* h_baseline = nullptr;
Double_t f_baseline(Double_t* x, Double_t* par)
{
  // par[0]: adc offset, par[1]: scale, par[2]: time offset
  if(!h_baseline){
    throw Exception("something is wrong in [RawData::CorrectBaselineTPC()]");
    // return TMath::QuietNaN();
  }
  Int_t floor = TMath::FloorNint(par[2]);
  Double_t frac = par[2] - floor;
  Int_t bin_left = h_baseline->GetXaxis()->FindBin(x[0] + floor);
  Double_t val_left = h_baseline->GetBinContent(bin_left);
  Double_t val_right = h_baseline->GetBinContent(bin_left+1);
  return
    par[0] + par[1]*((1-frac)*val_left + frac*val_right);
}
}

//_____________________________________________________________________________
RawData::RawData()
  : m_is_decoded(kNType),
    m_BH1RawHC(),
    m_BH2RawHC(),
    m_BACRawHC(),
    m_HTOFRawHC(),
    m_SCHRawHC(),
    m_BVHRawHC(),
    m_TOFRawHC(),
    m_LACRawHC(),
    m_WCRawHC(),
    m_WCSUMRawHC(),
    m_BFTRawHC(NumOfPlaneBFT),
    m_BcInRawHC(NumOfLayersBcIn+1),
    m_BcOutRawHC(NumOfLayersBcOut+1),
    m_TPCRawHC(NumOfLayersTPC+1),
    m_TPCCorHC(NumOfLayersTPC+1),
    m_SdcInRawHC(NumOfLayersSdcIn+1),
    m_SdcOutRawHC(NumOfLayersSdcOut+1),
    m_ScalerRawHC(),
    m_TrigRawHC(),
    m_VmeCalibRawHC(),
    m_baseline()
{
  for(auto& d: m_is_decoded) d = false;
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
RawData::~RawData()
{
  ClearAll();
  ClearTPC();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
RawData::ClearAll()
{
  del::ClearContainer(m_BH1RawHC);
  del::ClearContainer(m_BACRawHC);
  del::ClearContainer(m_BH2RawHC);
  del::ClearContainer(m_HTOFRawHC);
  del::ClearContainer(m_SCHRawHC);
  del::ClearContainer(m_BVHRawHC);
  del::ClearContainer(m_TOFRawHC);
  del::ClearContainer(m_LACRawHC);
  del::ClearContainer(m_WCRawHC);
  del::ClearContainer(m_WCSUMRawHC);
 del::ClearContainer(m_E72BACRawHC);
  del::ClearContainer(m_E72BACSUMRawHC);
  del::ClearContainer(m_KVCRawHC);
  del::ClearContainer(m_KVCSUMRawHC);
  del::ClearContainer(m_E90SACRawHC);
  del::ClearContainer(m_E90SACSUMRawHC);
  del::ClearContainer(m_T4RawHC);
  del::ClearContainer(m_T5RawHC);
  del::ClearContainer(m_T6RawHC);
  del::ClearContainer(m_T7RawHC);
  del::ClearContainerAll(m_BFTRawHC);
  del::ClearContainerAll(m_BcInRawHC);
  del::ClearContainerAll(m_BcOutRawHC);
  del::ClearContainerAll(m_SdcInRawHC);
  del::ClearContainerAll(m_SdcOutRawHC);
  del::ClearContainer(m_ScalerRawHC);
  del::ClearContainer(m_TrigRawHC);
  del::ClearContainer(m_VmeCalibRawHC);
}

//_____________________________________________________________________________
void
RawData::ClearTPC()
{
  del::ClearContainerAll(m_TPCRawHC);
  del::ClearContainerAll(m_TPCCorHC);
  del::ClearContainer(m_TPCClockRawHC);
}

//_____________________________________________________________________________
Bool_t
RawData::CorrectBaselineTPC()
{
  static const Int_t MinTimeBucket = gUser.GetParameter("TimeBucketTPC", 0);
  static const Int_t MaxTimeBucket = gUser.GetParameter("TimeBucketTPC", 1);

  if(!m_is_decoded[kTPC]){
    hddaq::cerr << FUNC_NAME << " DecodeTPCHits() must be done!" << std::endl;
    return false;
  }

  del::ClearContainerAll(m_TPCCorHC);

  TH1D h1("baseline", "Baseline", NumOfTimeBucket, 0, NumOfTimeBucket);
  h_baseline = &h1;

#if DebugEvDisp
  gStyle->SetOptStat(0);
  // gStyle->SetOptStat(1110);
  // gStyle->SetOptFit(1);
  static TCanvas c1("c"+FUNC_NAME, FUNC_NAME, 1200, 900);
  c1.cd();
  TH1D h2(FUNC_NAME+"-h2", "Corrected FADC",
          NumOfTimeBucket, 0, NumOfTimeBucket);
#endif

  m_baseline = nullptr;
  Double_t min_ref = 1e5;
  for(const auto& hc : m_TPCRawHC){
    for(const auto& hit : hc){
      if(hit->FadcSize() != NumOfTimeBucket)
        continue;
      ///// Minimum RMS method (unused)
      // auto ref = hit->RMS(MinTimeBucket, MaxTimeBucket);
      // if(hit->RMS(MinTimeBucket, MaxTimeBucket) < 20.) continue;
      ///// Minimum Amplitude method
      auto ref = hit->MaxAdc(MinTimeBucket, MaxTimeBucket)
        - hit->Mean(MinTimeBucket, MaxTimeBucket);
      if(ref < min_ref){
        min_ref = ref;
        m_baseline = hit;
      }
    }
  }
  if(!m_baseline || !h_baseline){
    hddaq::cerr << FUNC_NAME << " no reference was found." << std::endl;
    return false;
  }

  const auto& base_fadc = m_baseline->Fadc();
  const Double_t base_ped = m_baseline->Mean(0, MinTimeBucket);
  for(Int_t i=0; i<NumOfTimeBucket; ++i){
    h_baseline->SetBinContent(i+1, base_fadc.at(i) - base_ped);
  }
  h_baseline->SetTitle(Form("Baseline Layer#%d Row#%d",
                            m_baseline->LayerId(),
                            m_baseline->RowId()));

#if DebugEvDisp
  {
    m_baseline->Print();
  }
#endif

  for(const auto& hc : m_TPCRawHC){
    for(const auto& hit : hc){
      TH1D h_fadc("h_fadc", Form("FADC Layer#%d Row#%d;sample# ;ADC ch",
                                 hit->LayerId(), hit->RowId()),
                  NumOfTimeBucket, 0, NumOfTimeBucket);
      h_fadc.SetLineWidth(2);
      const auto& fadc = hit->Fadc();
      for(Int_t i=0, n=fadc.size(); i<n; ++i){
        h_fadc.SetBinContent(i+1, fadc.at(i));
      }
      TF1 f1("f1", f_baseline, 0, NumOfTimeBucket, 3);
      f1.SetParameter(0, hit->Mean(0, MinTimeBucket));
      f1.SetParameter(1, 1.);
      f1.SetParameter(2, 0.);
      f1.SetParLimits(0, 0, 4000.);
      f1.SetParLimits(1, -5, 5.);
      f1.SetParLimits(2, -10, 10);
      h_fadc.Fit("f1", "Q", "", 0, MinTimeBucket);
      Double_t max_cadc = -1e10;
      for(Int_t i=0, n=fadc.size(); i<n; ++i){
        Double_t cadc = fadc.at(i) - f1.Eval(i);
        AddTPCRawHit(m_TPCCorHC[hit->LayerId()], hit->LayerId(),
                     hit->RowId(), cadc, f1.GetParameters());
        if(cadc > max_cadc && i > 25){
          max_cadc = cadc;
        }
      }
#if DebugEvDisp
      if(max_cadc > 100 && max_cadc<1000)
      {
        h2.Reset();
        h2.SetLineWidth(2);
        Double_t max_adc = -1e10;
        for(Int_t i=0, n=fadc.size(); i<n; ++i){
          h2.SetBinContent(i+1, fadc.at(i) - f1.Eval(i));
          if(fadc.at(i) > max_adc) max_adc = fadc.at(i);
        }
        h_fadc.Draw();
        // h_fadc.SetMinimum(min_adc - 100);
        h_fadc.SetMinimum(-100);
        h_fadc.SetMaximum(max_adc + 100);
        // h_fadc.SetMaximum(1000);
        // h_fadc.SetMaximum(2000);
        TF1 f2("f2", f_baseline, MinTimeBucket, NumOfTimeBucket, 3);
        f2.SetParameters(f1.GetParameters());
        f2.Draw("same");
        h2.SetLineColor(kGreen+1);
        h2.SetLineWidth(2);
        h2.Draw("same");
        // h2.SetMinimum(-500);
        // h2.SetMaximum( 500);
        gPad->Modified();
        gPad->Update();
        c1.Print("c1.pdf");
        getchar();
      }
#endif
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::DecodeHits()
{
  static const Double_t MinTdcBC3  = gUser.GetParameter("TdcBC3", 0);
  static const Double_t MaxTdcBC3  = gUser.GetParameter("TdcBC3", 1);
  static const Double_t MinTdcBC4  = gUser.GetParameter("TdcBC4", 0);
  static const Double_t MaxTdcBC4  = gUser.GetParameter("TdcBC4", 1);
  static const Double_t MinTdcSDC1 = gUser.GetParameter("TdcSDC1", 0);
  static const Double_t MaxTdcSDC1 = gUser.GetParameter("TdcSDC1", 1);
  static const Double_t MinTdcSDC2 = gUser.GetParameter("TdcSDC2", 0);
  static const Double_t MaxTdcSDC2 = gUser.GetParameter("TdcSDC2", 1);
  static const Double_t MinTdcSDC3 = gUser.GetParameter("TdcSDC3", 0);
  static const Double_t MaxTdcSDC3 = gUser.GetParameter("TdcSDC3", 1);
  static const Double_t MinTdcSDC4 = gUser.GetParameter("TdcSDC4", 0);
  static const Double_t MaxTdcSDC4 = gUser.GetParameter("TdcSDC4", 1);
  static const Double_t MinTrailingSDC1 = gUser.GetParameter("TrailingSDC1", 0);
  static const Double_t MaxTrailingSDC1 = gUser.GetParameter("TrailingSDC1", 1);
  static const Double_t MinTrailingSDC2 = gUser.GetParameter("TrailingSDC2", 0);
  static const Double_t MaxTrailingSDC2 = gUser.GetParameter("TrailingSDC2", 1);
  static const Double_t MinTrailingSDC3 = gUser.GetParameter("TrailingSDC3", 0);
  static const Double_t MaxTrailingSDC3 = gUser.GetParameter("TrailingSDC3", 1);
  static const Double_t MinTrailingSDC4 = gUser.GetParameter("TrailingSDC4", 0);
  static const Double_t MaxTrailingSDC4 = gUser.GetParameter("TrailingSDC4", 1);

  if(m_is_decoded[kOthers]){
    hddaq::cout << "#D " << FUNC_NAME << " "
		<< "already decoded!" << std::endl;
    return false;
  }

  ClearAll();

  // BH1
  DecodeHodo(DetIdBH1, NumOfSegBH1, kBothSide, m_BH1RawHC);
  // BH2
  DecodeHodo(DetIdBH2, NumOfSegBH2, kBothSide, m_BH2RawHC);
  // BAC
  DecodeHodo(DetIdBAC, NumOfSegBAC, kOneSide,  m_BACRawHC);
  // HTOF
  DecodeHodo(DetIdHTOF, NumOfSegHTOF, kBothSide, m_HTOFRawHC);
  // E72BAC
  DecodeHodo(DetIdE72BAC, NumOfSegE72BAC, kOneSide, m_E72BACRawHC);
  // E72BACSUM
  DecodeHodo(DetIdE72BACSUM, NumOfSegE72BACSUM, kOneSide, m_E72BACSUMRawHC);
  // KVC
  DecodeHodo(DetIdKVC, NumOfSegKVC, kOneSide, m_KVCRawHC);
  // KVCSUM
  DecodeHodo(DetIdKVCSUM, NumOfSegKVCSUM, kOneSide, m_KVCSUMRawHC);
  // E90SAC
  DecodeHodo(DetIdE90SAC, NumOfSegE90SAC, kOneSide, m_E90SACRawHC);
  // E90SACSUM
  DecodeHodo(DetIdE90SACSUM, NumOfSegE90SACSUM, kOneSide, m_E90SACSUMRawHC);
  // T4
  DecodeHodo(DetIdT4, NumOfSegT4, kOneSide, m_T4RawHC);
  // T5
  DecodeHodo(DetIdT5, NumOfSegT5, kOneSide, m_T5RawHC);
  // T6
  DecodeHodo(DetIdT6, NumOfSegT6, kOneSide, m_T6RawHC);
  // T7
  DecodeHodo(DetIdT7, NumOfSegT7, kOneSide, m_T7RawHC);


  // SCH
  for(Int_t seg=0; seg<NumOfSegSCH; ++seg){
    UInt_t nhit = gUnpacker.get_entries(DetIdSCH, 0, seg, 0, 0);
    for(Int_t i=0; i<nhit; ++i){
      UInt_t leading  = gUnpacker.get(DetIdSCH, 0, seg, 0, 0, i);
      UInt_t trailing = gUnpacker.get(DetIdSCH, 0, seg, 0, 1, i);
      AddHodoRawHit(m_SCHRawHC, DetIdSCH, 0, seg, 0,
                    kHodoLeading,  leading);
      AddHodoRawHit(m_SCHRawHC, DetIdSCH, 0, seg, 0,
                    kHodoTrailing, trailing);
    }
  }
  // BVH
  DecodeHodo(DetIdBVH, NumOfSegBVH, kOneSide,  m_BVHRawHC);
  // TOF
  DecodeHodo(DetIdTOF, NumOfSegTOF, kBothSide, m_TOFRawHC);
  // LAC
  DecodeHodo(DetIdLAC, NumOfSegLAC, kOneSide, m_LACRawHC);
  // WC
  DecodeHodo(DetIdWC, NumOfSegWC, kBothSide, m_WCRawHC);
  // WC SUM
  for(Int_t seg=0; seg<NumOfSegWC; ++seg){
    for(Int_t AorT=0; AorT<2; ++AorT){
      for(Int_t m=0, nhit=gUnpacker.get_entries(DetIdWC, 0, seg, 2, AorT);
          m<nhit; ++m){
	UInt_t data = gUnpacker.get(DetIdWC, 0, seg, 2, AorT, m);
	AddHodoRawHit(m_WCSUMRawHC, DetIdWCSUM, 0, seg, 0, AorT, data);
      }
    }
  }

  //BFT
  for(Int_t plane=0; plane<NumOfPlaneBFT; ++plane){
    for(Int_t seg=0; seg<NumOfSegBFT; ++seg){
      UInt_t nhit = gUnpacker.get_entries(DetIdBFT, plane, 0, seg, 0);
      for(Int_t i=0; i<nhit; ++i){
	UInt_t leading  = gUnpacker.get(DetIdBFT, plane, 0, seg, 0, i);
	UInt_t trailing = gUnpacker.get(DetIdBFT, plane, 0, seg, 1, i);
	AddHodoRawHit(m_BFTRawHC[plane], DetIdBFT, plane, seg, 0,
                      kHodoLeading,  leading);
	AddHodoRawHit(m_BFTRawHC[plane], DetIdBFT, plane, seg, 0,
                      kHodoTrailing, trailing);
      }
    }
  }

  // BC3&BC4 MWDC
  for(Int_t plane=0; plane<NumOfLayersBcOut; ++plane){
    if(plane<NumOfLayersBc){
      for(Int_t wire=0; wire<MaxWireBC3; ++wire){
        for(Int_t lt=0; lt<2; ++lt){
	  UInt_t nhit = gUnpacker.get_entries(DetIdBC3, plane, 0, wire, lt);
#if OscillationCut
	  if(nhit>MaxMultiHitDC) continue;
#endif
	  for(Int_t i=0; i<nhit; i++){
	    UInt_t data = gUnpacker.get(DetIdBC3, plane, 0, wire, lt, i);
	    if (lt == 0 && (data<MinTdcBC3 || MaxTdcBC3<data)) continue;
	    if (lt == 1 && data<MinTdcBC3) continue;
	    AddDCRawHit(m_BcOutRawHC[plane+1], plane+PlMinBcOut, wire+1,
                        data, lt);
	  }
	}
      }
    } else {
      // BC3 U-U' is dead E42
      if(gUnpacker.get_run_number() >= 5361
         && (plane == NumOfLayersBc || plane == NumOfLayersBc+1)){
        continue;
      }
      //
      for(Int_t wire=0; wire<MaxWireBC4; ++wire){
        for(Int_t lt=0; lt<2; ++lt){
	  UInt_t nhit = gUnpacker.get_entries(DetIdBC4, plane-NumOfLayersBc,
                                              0, wire, lt);
#if OscillationCut
	  if(nhit>MaxMultiHitDC) continue;
#endif
	  for(Int_t i=0; i<nhit; i++){
	    UInt_t data =  gUnpacker.get(DetIdBC4, plane-NumOfLayersBc, 0,
                                         wire, lt, i);
	    if (lt == 0 && (data<MinTdcBC4 || MaxTdcBC4<data)) continue;
	    if (lt == 1 && data<MinTdcBC4) continue;
	    AddDCRawHit(m_BcOutRawHC[plane+1], plane+PlMinBcOut, wire+1,
                        data, lt);
	  }
        }
      }
    }
  }

  // SDC1
  for(Int_t plane=0; plane<NumOfLayersSDC1; ++plane){
    for(Int_t wire=0; wire<MaxWireSDC1; ++wire){
      for(Int_t lt=0; lt<2; ++lt){
	UInt_t nhit = gUnpacker.get_entries(DetIdSDC1, plane, 0, wire, lt);
#if OscillationCut
	if(nhit > MaxMultiHitDC) continue;
#endif
	for(Int_t i=0; i<nhit; ++i){
	  UInt_t data = gUnpacker.get(DetIdSDC1, plane, 0, wire, lt, i);
          if(lt == 0 && (data < MinTdcSDC1 || MaxTdcSDC1 < data)) continue;
          if(lt == 1 && (data < MinTrailingSDC1 || MaxTrailingSDC1 < data)) continue;
          AddDCRawHit(m_SdcInRawHC[plane+1], plane+PlMinSdcIn, wire+1, data , lt);
	}
      }
    }
  }
  // SDC2
  for(Int_t plane=0; plane<NumOfLayersSDC2; ++plane){
    const Int_t MaxWireSDC2 = (plane < 2) ? MaxWireSDC2X : MaxWireSDC2Y;
    for(Int_t wire=0; wire<MaxWireSDC2; ++wire){
      for(Int_t lt=0; lt<2; ++lt){
	UInt_t nhit = gUnpacker.get_entries(DetIdSDC2, plane, 0, wire, lt);
#if OscillationCut
	if(nhit > MaxMultiHitDC) continue;
#endif
	for(Int_t i=0; i<nhit; ++i){
	  UInt_t data = gUnpacker.get(DetIdSDC2, plane, 0, wire, lt, i);
	  if(lt == 0 && (data < MinTdcSDC2 || MaxTdcSDC2 < data)) continue;
          if(lt == 1 && (data < MinTrailingSDC2 || MaxTrailingSDC2 < data)) continue;
          AddDCRawHit(m_SdcInRawHC[plane+NumOfLayersSDC1+1],
                      plane+NumOfLayersSDC1+PlMinSdcIn, wire+1, data , lt);
	}
      }
    }
  }

  // SDC3
  for(Int_t plane=0; plane<NumOfLayersSDC3; ++plane){
    for(Int_t wire=0; wire<MaxWireSDC3; ++wire){
      for(Int_t lt=0; lt<2; ++lt){
        auto nhit = gUnpacker.get_entries(DetIdSDC3, plane, 0, wire, lt);
#if OscillationCut
        if(nhit > MaxMultiHitDC) continue;
#endif
        for(Int_t i=0; i<nhit; ++i){
          auto data = gUnpacker.get(DetIdSDC3, plane, 0, wire, lt, i);
          if(lt == 0 && (data < MinTdcSDC3 || MaxTdcSDC3 < data)) continue;
          if(lt == 1 && (data < MinTrailingSDC3 || MaxTrailingSDC3 < data)) continue;
          AddDCRawHit(m_SdcOutRawHC[plane+1], plane+PlMinSdcOut, wire+1,
                      data , lt);
        }
      }
    }
  }

  // SDC4
  for(Int_t plane=0; plane<NumOfLayersSDC4; ++plane){
    const Int_t MaxWireSDC4 = (plane < 2) ? MaxWireSDC4Y : MaxWireSDC4X;
    for(Int_t wire=0; wire<MaxWireSDC4; ++wire){
      for(Int_t lt=0; lt<2; ++lt){
        auto nhit = gUnpacker.get_entries(DetIdSDC4, plane, 0, wire, lt);
#if OscillationCut
        if(nhit > MaxMultiHitDC) continue;
#endif
        for(Int_t i=0; i<nhit; ++i){
          auto data = gUnpacker.get(DetIdSDC4, plane, 0, wire, lt, i);
          if(lt == 0 && (data < MinTdcSDC4 || MaxTdcSDC4 < data)) continue;
          if(lt == 1 && (data < MinTrailingSDC4 || MaxTrailingSDC4 < data)) continue;
          AddDCRawHit(m_SdcOutRawHC[plane+NumOfLayersSDC3+1],
                      plane+NumOfLayersSDC3+PlMinSdcOut, wire+1, data , lt);
        }
      }
    }
  }

  // Scaler
  for(Int_t l = 0; l<NumOfScaler; ++l){
    for(Int_t seg=0; seg<NumOfSegScaler; ++seg){
      UInt_t nhit = gUnpacker.get_entries(DetIdScaler, l, 0, seg, 0);
      if(nhit == 0) continue;
      UInt_t data = gUnpacker.get(DetIdScaler, l, 0, seg, 0);
      AddHodoRawHit(m_ScalerRawHC, DetIdScaler, l, seg, 0, 0, data);
    }
  }

  // Trigger Flag
  DecodeHodo(DetIdTrig, NumOfSegTrig, kOneSide, m_TrigRawHC);

  m_is_decoded[kOthers] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::DecodeTPCHits()
{
  static const auto k_tpc = gUnpacker.get_device_id("TPC");
  static const auto k_adc = gUnpacker.get_data_id("TPC", "adc");
  static const Bool_t BaselineCorrectionTPC
    = (gUser.GetParameter("BaselineCorrectionTPC") == 1);

  if(m_is_decoded[kTPC]){
    hddaq::cout << FUNC_NAME << " " << "already decoded!" << std::endl;
    return false;
  }

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
    for(Int_t r=0; r<NumOfRow; ++r){
      const auto nhit = gUnpacker.get_entries(k_tpc, layer, 0, r, k_adc);
      for(Int_t i=0; i<nhit; ++i){
        auto adc = gUnpacker.get(k_tpc, layer, 0, r, k_adc, i);
	AddTPCRawHit(m_TPCRawHC[layer], layer, r, adc);
      }
    }
  }

  {
    static const auto device_id = gUnpacker.get_device_id("HTOF");
    static const auto data_id = gUnpacker.get_data_id("HTOF", "fpga_leading");
    static const Int_t segment = 34;
    for(Int_t i=0, n=gUnpacker.get_entries(device_id, 0, segment, 0, data_id);
        i<n; ++i){
      auto tdc = gUnpacker.get(device_id, 0, segment, 0, data_id, i);
      AddHodoRawHit(m_TPCClockRawHC, 0, 0, 0, 0, kHodoLeading, tdc);
    }
  }

  m_is_decoded[kTPC] = true;

  if(BaselineCorrectionTPC)
    CorrectBaselineTPC();

  /*
   * if correction is skipped or null baseline is found,
   * m_TPCCorHC is deeply copied from m_TPCRawHC.
   * So, in any case, m_TPCCorHC will be used in DCAnalyzer too.
   */
  if(!m_baseline){
    del::ClearContainerAll(m_TPCCorHC);
    for(const auto& hc: m_TPCRawHC){
      for(const auto& hit: hc){
        for(const auto& adc: hit->Fadc()){
          AddTPCRawHit(m_TPCCorHC[hit->LayerId()], hit->LayerId(),
                       hit->RowId(), adc);
        }
      }
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::SelectTPCHits(Bool_t maxadccut, Bool_t maxadctbcut)
{
  if(!m_is_decoded[kTPC]){
    hddaq::cerr << FUNC_NAME << " "
		<< "Rawdata has not been decoded!" << std::endl;
    return false;
  }
  std::vector<TPCRHitContainer> ValidCand;
  ValidCand.resize(NumOfLayersTPC+1);
  std::vector<TPCRHitContainer> DeleteCand;
  DeleteCand.resize(NumOfLayersTPC+1);

  static const Double_t MinDe = gUser.GetParameter("MinDeTPC");
  static const Int_t MinTimeBucket = gUser.GetParameter("TimeBucketTPC", 0);
  static const Int_t MaxTimeBucket = gUser.GetParameter("TimeBucketTPC", 1);
  static const Bool_t BaselineCorrectionTPC
    = (gUser.GetParameter("BaselineCorrectionTPC") == 1);

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    std::size_t nh;
    if(BaselineCorrectionTPC)
      nh = m_TPCCorHC[layer].size();
    else
      nh = m_TPCRawHC[layer].size();
    if(nh==0)
      continue;
    for(std::size_t hiti =0; hiti< nh; ++hiti){
      TPCRawHit* hit;
      if(BaselineCorrectionTPC)
	hit = m_TPCCorHC[layer][hiti];
      else
	hit = m_TPCRawHC[layer][hiti];

      Double_t mean = hit->Mean();
      //Double_t max_adc = hit->MaxAdc() - mean;
      Double_t max_adc = hit->MaxAdc(MinTimeBucket, MaxTimeBucket) - mean;
      Int_t maxadc_tb = hit->LocMax();

      if(maxadccut&&maxadctbcut){
	if(max_adc>MinDe
	   && (MinTimeBucket < maxadc_tb
	       && maxadc_tb < MaxTimeBucket))
	  ValidCand[layer].push_back(hit);
	else
	  DeleteCand[layer].push_back(hit);
      }
      if(maxadccut&&!maxadctbcut){
	if(max_adc>MinDe)
	  ValidCand[layer].push_back(hit);
	else
	  DeleteCand[layer].push_back(hit);
      }
      if(!maxadccut&&maxadctbcut){
	if(MinTimeBucket < maxadc_tb
	   && maxadc_tb < MaxTimeBucket)
	  ValidCand[layer].push_back(hit);
	else
	  DeleteCand[layer].push_back(hit);
      }
      if(!maxadccut&&!maxadctbcut){
	ValidCand[layer].push_back(hit);
      }
    }
  }

  del::ClearContainerAll(DeleteCand);

  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    if(BaselineCorrectionTPC){
      m_TPCCorHC[layer].clear();
      m_TPCCorHC[layer].resize(ValidCand[layer].size());
      std::copy(ValidCand[layer].begin(), ValidCand[layer].end(), m_TPCCorHC[layer].begin());
      ValidCand[layer].clear();
    }
    else{
      m_TPCRawHC[layer].clear();
      m_TPCRawHC[layer].resize(ValidCand[layer].size());
      std::copy(ValidCand[layer].begin(), ValidCand[layer].end(), m_TPCRawHC[layer].begin());
      ValidCand[layer].clear();
    }
  }
  return true;
}


//_____________________________________________________________________________
Bool_t
RawData::DecodeCalibHits()
{
  del::ClearContainer(m_VmeCalibRawHC);

  for(Int_t plane=0; plane<NumOfPlaneVmeCalib; ++plane){
    DecodeHodo(DetIdVmeCalib, plane, NumOfSegVmeCalib,
               kOneSide, m_VmeCalibRawHC);
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::AddHodoRawHit(HodoRHitContainer& cont,
                       Int_t id, Int_t plane, Int_t seg,
                       Int_t UorD, Int_t type, Int_t data)
{
  HodoRawHit* p = nullptr;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    HodoRawHit* q = cont[i];
    if(q->DetectorId() == id &&
       q->PlaneId() == plane &&
       q->SegmentId() == seg){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit(id, plane, seg);
    cont.push_back(p);
  }

  switch(type){
  case kHodoAdc:
    if(UorD==0) p->SetAdcUp(data);
    else        p->SetAdcDown(data);
    break;
  case kHodoLeading:
    if(UorD==0) p->SetTdcUp(data);
    else        p->SetTdcDown(data);
    break;
  case kHodoTrailing:
    if(UorD==0) p->SetTdcTUp(data);
    else        p->SetTdcTDown(data);
    break;
  case kHodoOverflow:
    p->SetTdcOverflow(data);
    break;
  default:
    hddaq::cerr << FUNC_NAME << " wrong data type " << std::endl
		<< "DetectorId = " << id    << std::endl
		<< "PlaneId    = " << plane << std::endl
		<< "SegId      = " << seg   << std::endl
		<< "AorT       = " << type  << std::endl
		<< "UorD       = " << UorD  << std::endl;
    return false;
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::AddDCRawHit(DCRHitContainer& cont,
                     Int_t plane, Int_t wire, Int_t data, Int_t type)
{
  DCRawHit* p = nullptr;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    DCRawHit* q = cont[i];
    if(q->PlaneId()==plane &&
       q->WireId()==wire){
      p=q; break;
    }
  }
  if(!p){
    p = new DCRawHit(plane, wire);
    cont.push_back(p);
  }

  switch(type){
  case kDcLeading:
    p->SetTdc(data);
    break;
  case kDcTrailing:
    p->SetTrailing(data);
    break;
  case kDcOverflow:
    p->SetTdcOverflow(data);
    break;
  default:
    hddaq::cerr << FUNC_NAME << " wrong data type " << std::endl
		<< "PlaneId    = " << plane << std::endl
		<< "WireId     = " << wire  << std::endl
		<< "DataType   = " << type  << std::endl;
    return false;
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::AddTPCRawHit(TPCRHitContainer& cont,
                      Int_t layer, Int_t row, Double_t adc,
                      Double_t* pars)
{
  TPCRawHit* p = nullptr;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    TPCRawHit* q = cont[i];
    if(q->LayerId() == layer && q->RowId() == row){
      p=q; break;
    }
  }
  if(!p){
    p = new TPCRawHit(layer, row, pars);
    cont.push_back(p);
  }
  p->AddFadc(adc);
  return true;
}

//_____________________________________________________________________________
void
RawData::DecodeHodo(Int_t id, Int_t plane, Int_t nseg, Int_t nch,
                    HodoRHitContainer& cont)
{
  for(Int_t seg=0; seg<nseg; ++seg){
    for(Int_t UorD=0; UorD<nch; ++UorD){
      for(Int_t AorT=0; AorT<2; ++AorT){
	for(Int_t m=0, nhit=gUnpacker.get_entries(id, plane, seg, UorD, AorT);
            m<nhit; ++m){
	  UInt_t data = gUnpacker.get(id, plane, seg, UorD, AorT, m);
	  AddHodoRawHit(cont, id, plane, seg, UorD, AorT, data);
	}
      }
    }
  }
}

//_____________________________________________________________________________
void
RawData::DecodeHodo(Int_t id, Int_t nseg, Int_t nch, HodoRHitContainer& cont)
{
  DecodeHodo(id, 0, nseg, nch, cont);
}

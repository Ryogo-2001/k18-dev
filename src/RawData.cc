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
    m_CaenV792RawHC(),
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
  del::ClearContainer(m_CaenV792RawHC);
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
  if(m_is_decoded[kOthers]){
    hddaq::cout << "#D " << FUNC_NAME << " "
		<< "already decoded!" << std::endl;
    return false;
  }

  ClearAll();

  DecodeHodo(DetIdCaenV792, NumOfSegCaenV792, kOneSide,  m_CaenV792RawHC);

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

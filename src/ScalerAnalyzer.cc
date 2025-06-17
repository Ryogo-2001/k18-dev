// -*- C++ -*-

#include "ScalerAnalyzer.hh"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMath.h>
#include <TSystem.h>
#include <TTimeStamp.h>

#include <filesystem_util.hh>
#include <lexical_cast.hh>
#include <std_ostream.hh>
#include <Unpacker.hh>
#include <UnpackerManager.hh>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "Exception.hh"
#include "FuncName.hh"

namespace
{
  using namespace hddaq::unpacker;
  const UnpackerManager& gUnpacker = GUnpacker::get_instance();
}

//______________________________________________________________________________
ScalerAnalyzer::ScalerAnalyzer( void )
  : m_ost( std::cout ), //m_ost( hddaq::cout ),
    m_info( MaxColumn, std::vector<ScalerInfo>( MaxRow ) ),
    m_flag( nFlag, false ),
    m_spill_increment( false ),
    m_is_spill_end( false ),
    m_is_spill_on_end( false ),
    m_run_number( -1 ),
    m_canvas()
{
  for( Int_t i=0; i<MaxColumn; ++i ){
    for( Int_t j=0; j<MaxRow; ++j ){
      m_info[i][j] = ScalerInfo( "n/a", i, j, false );
    }
  }
}

//______________________________________________________________________________
ScalerAnalyzer::~ScalerAnalyzer( void )
{
}

//______________________________________________________________________________
void
ScalerAnalyzer::Clear( Option_t* opt )
{
  for( Int_t i=0; i<MaxColumn; ++i ){
    for( Int_t j=0; j<MaxRow; ++j ){
      if( !TString(opt).EqualTo("all") &&
	  m_info[i][j].name.EqualTo("Spill") ){
	continue;
      }
      m_info[i][j].data = 0;
    }
  }
}

//______________________________________________________________________________
Bool_t
ScalerAnalyzer::Decode( void )
{
  m_spill_increment = false;
  m_is_spill_end    = false;
  m_is_spill_on_end = false;

  //////////////////// Run Number
  if( m_run_number != gUnpacker.get_root()->get_run_number() ){
    m_run_number = gUnpacker.get_root()->get_run_number();
    Clear( "all" );
  }

  //////////////////// Trigger Flag
  std::bitset<NumOfSegTrig> trigger_flag;
  {
    static const auto k_device = gUnpacker.get_device_id( "TFlag" );
    static const auto k_tdc    = gUnpacker.get_data_id( "TFlag", "tdc" );
    for( Int_t seg=0; seg<NumOfSegTrig; ++seg ){
      for( Int_t i=0, n=gUnpacker.get_entries( k_device, 0, seg, 0, k_tdc );
	   i<n; ++i ){
	auto tdc = gUnpacker.get( k_device, 0, seg, 0, k_tdc, i );
	if( tdc>0 ) trigger_flag.set(seg);
	if( trigger_flag[seg] ) break;
      }
    }
  }
  if( trigger_flag[trigger::kSpillEnd] ){
    m_is_spill_end = true;
  }
  if( trigger_flag[trigger::kSpillOnEnd] )
    m_is_spill_on_end = true;
  if( m_flag[kSpillOn] ){
    if( !trigger_flag[trigger::kSpillOnEnd] &&
	!trigger_flag[trigger::kL1SpillOn] ){
      return false;
    }
  }
  auto spill_off_end = ( trigger_flag[trigger::kSpillEnd] &&
			 !trigger_flag[trigger::kSpillOnEnd] );
  if( m_flag[kSpillOff] ){
    if( !spill_off_end &&
	!trigger_flag[trigger::kL1SpillOff] ){
      return false;
    }
  }
  if( m_flag[kScalerSheet] && !m_is_spill_end ){
    for( Int_t i=0; i<MaxColumn; ++i ){
      for( Int_t j=0; j<MaxRow; ++j ){
	m_info[i][j].prev = 0;
      }
    }
    return true;
  }

  //////////////////// Scaler Data
  {
    static const Int_t device_id  = gUnpacker.get_device_id("Scaler");

    for( Int_t i=0; i<MaxColumn; ++i ){
      for( Int_t j=0; j<MaxRow; ++j ){
	if( !m_info[i][j].flag_disp )
	  continue;
	if( m_info[i][j].name.EqualTo("Spill") )
	  continue;

	Int_t module_id = m_info[i][j].module_id;
	Int_t channel   = m_info[i][j].channel;

	if( module_id < 0 || channel < 0 )
	  continue;

	Int_t nhit = gUnpacker.get_entries( device_id, module_id, 0, channel, 0 );
	if( nhit<=0 ) continue;
	Scaler val = gUnpacker.get( device_id, module_id, 0, channel, 0 );

	if( m_info[i][j].prev > val ){
	  m_spill_increment = true;
	  m_info[i][j].prev = 0;
	}

    	if( m_flag[kSpillBySpill] && m_spill_increment )
    	  m_info[i][j].data = 0;

	m_info[i][j].curr  = val;
	m_info[i][j].data += val - m_info[i][j].prev;
	m_info[i][j].prev  = m_info[i][j].curr;
      }
    }
  }

  //////////////////// for BH1 SUM
  {
    if( Has("BH1-SUM") ){
      static Channel ch = Find("BH1-SUM");
      m_info[ch.first][ch.second].data = 0;
      for( Int_t i=0; i<NumOfSegBH1; ++i ){
	m_info[ch.first][ch.second].data += Get( Form("BH1-%02d", i+1) );
      }
    }
  }

  //////////////////// for BH2 SUM
  {
    if( Has("BH2-SUM") ){
      static Channel ch = Find("BH2-SUM");
      m_info[ch.first][ch.second].data = 0;
      // static const Int_t n = GetFlag( kScalerE42 ) ? NumOfSegBH2_E42 : NumOfSegBH2;
      for( Int_t i=0; i<NumOfSegBH2; ++i ){
	m_info[ch.first][ch.second].data += Get( Form("BH2-%02d", i+1) );
      }
    }
  }

  //////////////////// for SCH SUM
  {
    if( Has("SCH-SUM") ){
      static Channel ch = Find("SCH-SUM");
      m_info[ch.first][ch.second].data = 0;
      for( Int_t i=0; i<NumOfSegSCH; ++i ){
	m_info[ch.first][ch.second].data += Get( Form("SCH-%02d", i+1) );
      }
    }
  }

  //////////////////// for LAC SUM
  {
    if( Has("LAC-SUM") ){
      static Channel ch = Find("LAC-SUM");
      m_info[ch.first][ch.second].data = 0;
      for( Int_t i=0; i<NumOfSegLAC/2; ++i ){
	m_info[ch.first][ch.second].data += Get( Form("LAC-%02d", i+1) );
      }
    }
  }

  //////////////////// Spill
  {
    if( Has("Spill") ){
      static Bool_t first = true;
      static Channel p = Find("Spill");
      if( first && !m_flag[kScalerSheet] ){
	m_info[p.first][p.second].data++;
	first = false;
      }
      if( m_spill_increment ||
	  ( m_flag[kScalerSheet] && m_is_spill_end ) ){
	m_info[p.first][p.second].data++;
      }
    }
  }
  return true;
}

//______________________________________________________________________________
Double_t
ScalerAnalyzer::Duty( void ) const
{
  Double_t daq_eff  = Fraction("L1-Acc","L1-Req");
  Double_t live_eff = Fraction("Live-Time","Real-Time");
  Double_t duty = daq_eff/(1.-daq_eff)*(1./live_eff-1.);
  if( duty > 1. || TMath::IsNaN(duty) )
    return 1.;
  else
    return duty;
}

//______________________________________________________________________________
void
ScalerAnalyzer::DrawOneBox( Double_t x, Double_t y,
			    const TString& title1, const TString& val1 )
{
  TLatex tex;
  tex.SetNDC();
  if( m_flag[kScalerHBX] )
    tex.SetTextSize(0.036);
  else
    tex.SetTextSize(0.04);
  tex.SetTextAlign(12);
  tex.DrawLatex( x, y, title1 );
  tex.SetTextAlign(32);
  if( x > 0.1 && x < 0.5 )
    tex.DrawLatex( x+0.295, y, val1 );
  else
    tex.DrawLatex( x+0.28, y, val1 );
}

//______________________________________________________________________________
void
ScalerAnalyzer::DrawOneLine( const TString& title1, const TString& val1,
			     const TString& title2, const TString& val2,
			     const TString& title3, const TString& val3 )
{
  static Int_t i = 0;
  if( m_flag[kScalerHBX] ){
    const Double_t ystep = 0.04;
    const Double_t y0 = 0.98;
    Double_t y = y0 - (i+1)*ystep;
    Double_t x[] = { 0.05, 0.35, 0.67 };
    DrawOneBox( x[0], y, title1, val1 );
    DrawOneBox( x[1], y, title2, val2 );
    DrawOneBox( x[2], y, title3, val3 );
    TLine line;
    line.SetNDC();
    line.SetLineColor(kGray);
    line.DrawLine( 0.05, y-0.5*ystep, 0.95, y-0.5*ystep );
    line.SetLineColor(kBlack);
    line.DrawLine( 0.34, y-0.5*ystep, 0.34, y+0.5*ystep );
    line.DrawLine( 0.66, y-0.5*ystep, 0.66, y+0.5*ystep );
    if( i==1 || i==17 )
      line.DrawLine( 0.05, y-0.5*ystep, 0.95, y-0.5*ystep );
  } else {
    const Double_t ystep = 0.05;
    const Double_t y0 = 0.95;
    Double_t y = y0 - (i+1)*ystep;
    Double_t x[] = { 0.05, 0.35, 0.67 };
    DrawOneBox( x[0], y, title1, val1 );
    DrawOneBox( x[1], y, title2, val2 );
    DrawOneBox( x[2], y, title3, val3 );
    TLine line;
    line.SetNDC();
    line.SetLineColor(kGray);
    line.DrawLine( 0.05, y-0.5*ystep, 0.95, y-0.5*ystep );
    line.SetLineColor(kBlack);
    line.DrawLine( 0.34, y-0.5*ystep, 0.34, y+0.5*ystep );
    line.DrawLine( 0.66, y-0.5*ystep, 0.66, y+0.5*ystep );
    if( i==0 || i==4 || i==13 )
      line.DrawLine( 0.05, y-0.5*ystep, 0.95, y-0.5*ystep );
  }
  ++i;
}

//______________________________________________________________________________
void
ScalerAnalyzer::DrawOneLine( const TString& title1,
			     const TString& title2,
			     const TString& title3 )
{
  DrawOneLine( title1, SeparateComma( Get(title1) ),
	       title2, SeparateComma( Get(title2) ),
	       title3, SeparateComma( Get(title3) ) );
}

//______________________________________________________________________________
Channel
ScalerAnalyzer::Find( const TString& name ) const
{
  for( Int_t i=0; i<MaxColumn; ++i ){
    for( Int_t j=0; j<MaxRow; ++j ){
      if( m_info[i][j].name.EqualTo( name ) ){
	return Channel( i, j );
      }
    }
  }

  m_ost << "#W " << FUNC_NAME << " "
	<< "no such name : " << name << std::endl;

  return Channel( -1, -1 );
}

//______________________________________________________________________________
Double_t
ScalerAnalyzer::Fraction( const TString& num, const TString& den ) const
{
  return (Double_t)Get( num ) / Get( den );
}

//______________________________________________________________________________
Scaler
ScalerAnalyzer::Get( const TString& name ) const
{
  for( Int_t i=0; i<MaxColumn; ++i ){
    for( Int_t j=0; j<MaxRow; ++j ){
      if( m_info[i][j].name.EqualTo( name ) ){
	return m_info[i][j].data;
      }
    }
  }

  m_ost << "#W " << FUNC_NAME << " "
	<< "no such ScalerInfo : " << name << std::endl;

  return 0;
}

//______________________________________________________________________________
Bool_t
ScalerAnalyzer::Has( const TString& name ) const
{
  for( Int_t i=0; i<MaxColumn; ++i ){
    for( Int_t j=0; j<MaxRow; ++j ){
      if( m_info[i][j].name.EqualTo( name ) ){
	return true;
      }
    }
  }

  return false;
}

//______________________________________________________________________________
Bool_t
ScalerAnalyzer::MakeScalerText( void ) const
{
  const Int_t run_number = gUnpacker.get_root()->get_run_number();
  const TString& bin_dir( hddaq::dirname(hddaq::selfpath()) );
  const TString& data_dir( hddaq::dirname(gUnpacker.get_istream()) );

  std::stringstream run_number_ss; run_number_ss << run_number;
  const TString& recorder_log( data_dir+"/recorder.log" );
  std::ifstream ifs( recorder_log );
  if( !ifs.is_open() ){
    std::cerr << "#E " << FUNC_NAME << " "
              << "cannot open recorder.log : "
              << recorder_log << std::endl;
    return false;
  }

  const TString& scaler_dir( bin_dir+"/../scaler" );
  const TString& scaler_txt = Form( "%s/scaler_%05d.txt",
                                    scaler_dir.Data(), run_number );

  std::ofstream ofs( scaler_txt );
  if( !ofs.is_open() ){
    std::cerr << "#E " << FUNC_NAME << " "
              << "cannot open scaler.txt : "
              << scaler_txt << std::endl;
    return false;
  }

  // Int_t recorder_event_number = 0;
  Bool_t found_run_number = false;
  std::string line;
  while( ifs.good() && std::getline(ifs,line) ){
    if( line.empty() ) continue;
    std::istringstream input_line( line );
    std::istream_iterator<std::string> line_begin( input_line );
    std::istream_iterator<std::string> line_end;
    std::vector<std::string> log_column( line_begin, line_end );
    if( log_column.at(0) != "RUN" ) continue;
    if( log_column.at(1) != run_number_ss.str() ) continue;
    // recorder_event_number = hddaq::a2i( log_column.at(15) );
    ofs << line << std::endl;
    found_run_number = true;
  }

  if( !found_run_number ){
    std::cerr << "#E " << FUNC_NAME << " "
              << "not found run# " << run_number
              << " in " << recorder_log << std::endl;
    return false;
  }
  return true;
}

//______________________________________________________________________________
void
ScalerAnalyzer::Print( Option_t* ) const
{
  m_ost << "\033[2J" << std::endl;

  TString end_mark = ( m_is_spill_on_end ? "Spill On End" :
		       m_is_spill_end ? "Spill Off End" : "" );

  Int_t event_number = gUnpacker.get_event_number();
  if( GetFlag( kScalerDaq ) || GetFlag( kScalerE42 ) ){
    m_ost << std::left  << std::setw(16) << "RUN"
	  << std::right << std::setw(16) << SeparateComma( m_run_number ) << std::endl
	  << std::left  << std::setw(16) << "Event Number"
	  << std::right << std::setw(16) << SeparateComma( event_number ) << std::endl
	  << std::left  << std::setw(16) << ""
	  << std::right << std::setw(16) << end_mark << std::endl;
    for( Int_t i=0; i<MaxRow; ++i ){
      if( !m_info[kLeft][i].name.Contains("n/a") ){
	m_ost << std::left  << std::setw(16) << m_info[kLeft][i].name
	      << std::right << std::setw(16) << SeparateComma( m_info[kLeft][i].data )
	      << std::endl;
      }
    }
    if( GetFlag( kScalerDaq ) ){
      m_ost << std::endl  << std::setprecision(6) << std::fixed
	    << std::left  << std::setw(16) << "Live/Real"
	    << std::right << std::setw(16) << Fraction("Live-Time","Real-Time") << std::endl
	    << std::left  << std::setw(16) << "DAQ-Eff"
	    << std::right << std::setw(16) << Fraction("L1-Acc","L1-Req") << std::endl
	    << std::left  << std::setw(16) << "L2-Eff"
	    << std::right << std::setw(16) << Fraction("L2-Acc","L1-Acc") << std::endl
	    << std::left  << std::setw(16) << "Duty-Factor"
	    << std::right << std::setw(16) << Duty() << std::endl;
    }
  } else {
    m_ost << std::left  << std::setw(16) << "RUN"
	  << std::right << std::setw(16) << SeparateComma( m_run_number ) << " : "
	  << std::left  << std::setw(16) << "Event Number"
	  << std::right << std::setw(16) << SeparateComma( event_number );
    m_ost << " : "
	  << std::left  << std::setw(16) << ""
	  << std::right << std::setw(16) << end_mark
	  << std::endl << std::endl;
    for( Int_t i=0; i<MaxRow; ++i ){
      if( m_info[kLeft][i].name.Contains("n/a") )
        continue;
      m_ost << std::left  << std::setw(16) << m_info[kLeft][i].name
	    << std::right << std::setw(16) << SeparateComma( m_info[kLeft][i].data )
	    << " : "
	    << std::left  << std::setw(16) << m_info[kCenter][i].name
	    << std::right << std::setw(16) << SeparateComma( m_info[kCenter][i].data ) << " : "
	    << std::left  << std::setw(16) << m_info[kRight][i].name
	    << std::right << std::setw(16) << SeparateComma( m_info[kRight][i].data )
	    <<std::endl;
    }
    if( !GetFlag( kScalerSch ) && !GetFlag( kScalerE42 ) && !GetFlag( kScalerHBX ) ){
      m_ost << std::endl  << std::setprecision(6) << std::fixed
	    << std::left  << std::setw(16) << "BH2/TM"
	    << std::right << std::setw(16) << Fraction("BH2", "TM") << " : "
	    << std::left  << std::setw(16) << "Live/Real"
	    << std::right << std::setw(16) << Fraction("Live-Time","Real-Time") << " : "
	    << std::left  << std::setw(16) << "DAQ-Eff"
	    << std::right << std::setw(16) << Fraction("L1-Acc","L1-Req") << std::endl
	    << std::left  << std::setw(16) << "L1Req/BH2"
	    << std::right << std::setw(16) << Fraction("L1-Req", "BH2") << " : "
	    << std::left  << std::setw(16) << "L2-Eff"
	    << std::right << std::setw(16) << Fraction("L2-Acc","L1-Acc") << " : "
	    << std::left  << std::setw(16) << "Duty-Factor"
	    << std::right << std::setw(16) << Duty() << std::endl
	    << std::endl;
    }
  }
}

//______________________________________________________________________________
TString
ScalerAnalyzer::SeparateComma( Scaler number ) const
{
  if( m_flag[kSeparateComma] ){
    std::vector<Scaler> sep_num;

    while( number/1000 ){
      sep_num.push_back( number%1000 );
      number /= 1000;
    }

    std::stringstream ss;  ss << number;
    std::vector<Scaler>::reverse_iterator
      itr, itr_end = sep_num.rend();
    for( itr=sep_num.rbegin(); itr!=itr_end; ++itr ){
      ss << "," << std::setfill('0') << std::setw(3) << *itr;
    }

    return TString( ss.str() );
  }
  else {
    return TString::LLtoa( number, 10 );
  }
}

//______________________________________________________________________________
void
ScalerAnalyzer::PrintFlags( void ) const
{
  m_ost << "#D " << FUNC_NAME << std::endl << std::left;
  for( Int_t i=0; i<nFlag; ++i ){
    m_ost << " key = " << std::setw(20) << sFlag[i]
	  << " val = " << m_flag[i] << std::endl;
  }
}

//______________________________________________________________________________
void
ScalerAnalyzer::PrintScalerSheet( void )
{
  if( !m_canvas )
    m_canvas = new TCanvas( "c1", "c1", 1200, 800 );
  else
    m_canvas->Clear();

  TTimeStamp stamp;
  stamp.Add( -stamp.GetZoneOffset() );

  DrawOneLine( stamp.AsString("s"), "",
               "Event#", SeparateComma( gUnpacker.get_event_number() ),
               Form("#color[%d]{Run#}", kRed+1), SeparateComma( m_run_number ) );

  TString mode = ( m_flag[kSpillOn] ? "Spill On" :
                   m_flag[kSpillOff] ? Form( "#color[%d]{Spill Off}", kBlue+1 ) :
                   "Unknown" );

  if( m_flag[kScalerHBX] ){
    DrawOneLine( mode, SeparateComma( Get("Spill") ),
                 "Clock", SeparateComma( Get("10M-Clock") ),
                 "Level1-PS", SeparateComma( Get("Level1-PS") ) );
    // for( Int_t i=0; i<NumOfSegGe; ++i )
    //   DrawOneLine( Form( "TFA-%02d", i ),
    //                Form( "CRM-%02d", i ),
    //                Form( "Reset-%02d", i ) );
    DrawOneLine( "LSO1", "LSO1High", "GeCoin1" );
    DrawOneLine( "LSO2", "LSO2High", "GeCoin2" );
    DrawOneLine( "LSO1xGe13", "GeHigh1", "Other3" );
    DrawOneLine( "LSO2xGe24", "GeHigh2", "Other4" );
  } else {
    DrawOneLine( mode, SeparateComma( Get("Spill") ),
                 "BH1", SeparateComma( Get("BH1") ),
                 "TM", SeparateComma( Get("TM") ) );
    DrawOneLine( "Clock", SeparateComma( Get("10M-Clock") ),
                 "BH2", SeparateComma( Get("BH2") ),
                 "SY", SeparateComma( Get("SY") ) );
    DrawOneLine( "K-Beam", SeparateComma( Get("BEAM-E") ),
                 "BH1-SUM", SeparateComma( Get("BH1-SUM") ),
                 "BH1-1/100-PS", SeparateComma( Get("BH1-1/100-PS") ) );
    DrawOneLine( "#pi-Beam", SeparateComma( Get("BEAM-F") ),
                 "BH2-SUM", SeparateComma( Get("BH2-SUM") ),
                 "TOF-24", SeparateComma( Get("TOF-24") ) );

    DrawOneLine( "BH1",  "TRIG-A",  "L1-Req"    );
    DrawOneLine( "BH2",  "TRIG-B",  "L1-Acc"    );
    DrawOneLine( "BAC",  "TRIG-C",  "Clear"     );
    DrawOneLine( "PVAC", "TRIG-D",  "L2-Req"    );
    DrawOneLine( "FAC",  "TRIG-E",  "L2-Acc"    );
    DrawOneLine( "SCH",  "TRIG-F",  "LSO1xGe13" );
    DrawOneLine( "TOF",  "Mtx2D-1", "LSO2xGe24" );
    DrawOneLine( "LAC",  "Mtx2D-2", "GeCoin1"   );
    DrawOneLine( "WC",   "Mtx3D",   "GeCoin2"   );
    //DrawOneLine( "TOF-24", "Other4",   "p-Scat"  );

    DrawOneLine( "BH2/TM",   Form("%.6lf", Fraction("BH2","TM")),
                 "Live/Real", Form("%.6lf", Fraction("Live-Time","Real-Time")),
                 "DAQ Eff",   Form("%.6lf", Fraction("L1-Acc","L1-Req")) );
    DrawOneLine( "L1Req/BH2", Form("%.6lf", Fraction("L1-Req","BH2")),
                 "L2 Eff",       Form("%.6lf", Fraction("L2-Acc","L1-Acc")),
                 "Duty Factor",  Form("%.6lf", Duty()) );
  }

  const TString& scaler_sheet_pdf( "/tmp/scaler_sheet.pdf" );
  m_canvas->Print( scaler_sheet_pdf );

  const TString& print_command( "lpr "+scaler_sheet_pdf );
  gSystem->Exec( print_command );
}

//______________________________________________________________________________
void
ScalerAnalyzer::Set( Int_t i, Int_t j, const ScalerInfo& info )
{
  if( i >= MaxColumn || j >= MaxRow ){
    m_ost << "#E " << FUNC_NAME << std::endl
	  << " * Exceed Column/Row : " << info.name << std::endl;
    throw Exception( FUNC_NAME+" "+info.name );
  } else {
    m_info[i][j] = info;
  }
}

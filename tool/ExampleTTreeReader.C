// -*- C++ -*-

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

//_____________________________________________________________________________
void
ExampleTTreeReader( void )
{
  TFile::Open( "tmp.root" );
  TTreeReader reader( "tpc", gFile );
  TTreeReaderValue<Int_t> runnum( reader, "runnum" );
  TTreeReaderValue<Int_t> evnum( reader, "evnum" );
  TTreeReaderValue<Int_t> ntTpc( reader, "ntTpc" );
  TTreeReaderValue<std::vector<Double_t>> chisqr( reader, "chisqr" );
  while( reader.Next() ){
    // Int_t ievent = reader.GetCurrentEntry();

    // basic for
    for( Int_t i=0; i<*ntTpc; ++i ){
      std::cout << (*chisqr).at(i) << std::endl;
    }

    // range-based for
    // for( const auto& val : *chisqr ){
    //   std::cout << val << std::endl;
    // }
  }
}

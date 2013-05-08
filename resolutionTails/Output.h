// $Id: $

#ifndef RESOLUTION_TAILS_OUTPUT
#define RESOLUTION_TAILS_OUTPUT

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TNamed.h"
#include "TROOT.h"
#include "TString.h"


// Write output to ROOT, eps, and/or ASCII file

namespace resolutionTails {
  class Output {
  public:
    static bool DEBUG;

    Output(const TString &prefix, bool toROOT, bool toEPS, bool archivePlots);
    ~Output();

    TString namePrefix() const { return prefix_; }

    void toFiles(TNamed* obj, const TString &name = "");
    void toFiles(TCanvas* obj, const TString &name = "");
    TString toFile(const TString &text, const TString &name, bool append = true) const;


  private:
    const TString prefix_;
    const bool toROOT_;
    const bool toEPS_;
    const bool tarPlots_;

    TFile* outFile_;
  };


  bool Output::DEBUG = false;


  
  // ------------------------------------------------------------------------------------
  Output::Output(const TString &prefix, bool toROOT, bool toEPS, bool archivePlots)
    : prefix_(prefix), toROOT_ (toROOT), toEPS_(toEPS), tarPlots_(archivePlots) {
    outFile_ = new TFile(prefix_+".root","RECREATE");
    if( !outFile_->IsOpen() ) {
      std::cerr << "Output: ERROR creating output file '" << prefix_ << ".root'" << std::endl;
      exit(-1);
    }
  }


  // ------------------------------------------------------------------------------------
  Output::~Output() {
    // Close ROOT file
    outFile_->Close();
    delete outFile_;

    // Clean up working directory
    if( tarPlots_ ) {
      std::cout << "Cleaning up working directory" << std::endl;
      if( toEPS_ ) {
	TString filesInTar = prefix_+"*.tex "+prefix_+"*.eps ";
	gROOT->ProcessLine(".! tar -zcf "+prefix_+".tar.gz "+filesInTar);
	gROOT->ProcessLine(".! rm "+filesInTar);
	std::cout << "  Plots in eps format and .tex files: "+prefix_+".tar.gz" << std::endl;
      } else {
	gROOT->ProcessLine(".! rm "+prefix_+"*.tex");
      }
      std::cout << "  Plots in ROOT format: "+prefix_+".root" << std::endl;
    }
  }


  
  // ------------------------------------------------------------------------------------
  void Output::toFiles(TNamed* obj, const TString &name) {
    if( toROOT_ ) {
      TString objName = obj->GetName();
      TString objTitle = obj->GetTitle();
      if( name != "" ) {
	obj->SetName(prefix_+"_"+name);
	obj->SetTitle(prefix_+"_"+name);
      }
      outFile_->WriteTObject(obj);
      obj->SetName(objName);
      obj->SetTitle(objTitle);
    }
  }


  // ------------------------------------------------------------------------------------
  void Output::toFiles(TCanvas* obj, const TString &name) {
    if( toEPS_ ) {
      TString objName = obj->GetName();
      TString objTitle = obj->GetTitle();
      if( name != "" ) {
	obj->SetName(prefix_+"_"+name);
	obj->SetTitle(prefix_+"_"+name);
      }
      obj->SaveAs(prefix_+"_"+name+".eps","eps");
      if( toROOT_ ) outFile_->WriteTObject(obj);
      obj->SetName(objName);
      obj->SetTitle(objTitle);
    }
  }


  // Write content of 'text' to ASCII file namePrefix()'name'.
  // 'name' should contain the file type, e.g. '.tex'.
  // 'append' defines whether 'text' is appended to an already
  // existing file or whether that file is replaced.
  //
  // Returns name of created file.
  // ------------------------------------------------------------------------------------
  TString Output::toFile(const TString &text, const TString &name, bool append) const {
    TString outFileName = prefix_;
    if( !(name.BeginsWith(".") || name == "") ) outFileName += "_";
    outFileName += name;
    ofstream file(outFileName,(append ? std::ios::app : std::ios::trunc));
    if( file.good() ) {
      file << text;
      file.close();
    } else {
      std::cerr << "\n\nERROR writing to file '" << outFileName << "'\n\n\n";
    }

    return outFileName;
  }
}
#endif

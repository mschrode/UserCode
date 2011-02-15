// $Id: $

#include "OutputManager.h"


#include <iostream>

namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  OutputManager* OutputManager::createOutputManager(OutputManager::Mode mode, const TString &fileNameBase) {
    OutputManager* out = 0;
    if( mode == OutputManager::EPSSingleFiles ) {
      out = new OutputManagerEPSSingleFiles(fileNameBase);
    } else if( mode == OutputManager::PSAllInOne ) {
      out = new OutputManagerPSAllInOne(fileNameBase);
    } else {
      std::cerr << "ERROR in OutputManager::createOutputManager(): No Mode '" << mode << "'" << std::endl;
      exit(1);
    }

    return out;
  }


  // -------------------------------------------------------------------------------------
  OutputManager::OutputManager(const TString &fileNameBase) {
    can_ = new TCanvas("OutputManager:Canvas","",500,500);
  }

  
  // -------------------------------------------------------------------------------------
  OutputManager::~OutputManager() {
    delete can_;
  }




  // -------------------------------------------------------------------------------------
  OutputManagerPSAllInOne::OutputManagerPSAllInOne(const TString &fileNameBase) 
    : OutputManager(fileNameBase) {
    ps_ = new TPostScript(fileNameBase+".ps",111);
    multiCan_ = new TCanvas("OutputManager:MultiCanvas","",450,600);
    multiCan_->Divide(3,4);
    multiCanPadIdx_ = 1;
  }


  // -------------------------------------------------------------------------------------
  OutputManagerPSAllInOne::~OutputManagerPSAllInOne() {
    ps_->NewPage();
    ps_->Close();
    delete ps_;
    delete multiCan_;
  }


  // -------------------------------------------------------------------------------------
  void OutputManagerPSAllInOne::newPage(const TString &title) {
    ps_->NewPage();
    multiCanPadIdx_ = 1;
  }


  // -------------------------------------------------------------------------------------
  void OutputManagerPSAllInOne::nextMultiPad(const TString &title) {
    if( multiCanPadIdx_ == 12 ) {
      ps_->NewPage();
      multiCanPadIdx_ = 1;
    }
    lastPad_ = multiCan_->cd(multiCanPadIdx_);
    lastPad_->SetLogx(0);
    lastPad_->SetLogy(0);
    ++multiCanPadIdx_;
  }


  // -------------------------------------------------------------------------------------
  void OutputManagerPSAllInOne::nextPad(const TString &title) {
    ps_->NewPage();
    lastPad_ = can_->cd();
    lastPad_->SetLogx(0);
    lastPad_->SetLogy(0);
  }


  // -------------------------------------------------------------------------------------
  void OutputManagerPSAllInOne::saveCurrentPad(const TString &name) {
    lastPad_->Draw();
  }




  // -------------------------------------------------------------------------------------
  OutputManagerEPSSingleFiles::OutputManagerEPSSingleFiles(const TString &fileNameBase) 
    : OutputManager(fileNameBase) { }


  // -------------------------------------------------------------------------------------
  void OutputManagerEPSSingleFiles::nextPad(const TString &title) {
    if( title.Length() ) can_->SetTitle(title);
    else can_->SetTitle("");
    can_->cd();
    can_->SetLogx(0);
    can_->SetLogy(0);
  }


  // -------------------------------------------------------------------------------------
  void OutputManagerEPSSingleFiles::saveCurrentPad(const TString &name) {
    can_->SaveAs(name,"eps");
  }
}

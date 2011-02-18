// $Id: OutputManager.cc,v 1.2 2011/02/17 13:42:32 mschrode Exp $

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
    topCan_ = util::HistOps::createRatioTopCanvas();
    bottomPad_ = 0;
    lastPad_ = can_;
  }

  
  // -------------------------------------------------------------------------------------
  OutputManager::~OutputManager() {
    delete can_;
    delete topCan_;
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
  void OutputManagerPSAllInOne::logx() {
    lastPad_->SetLogx(1);
    if( lastPad_ == topCan_ ) bottomPad_->SetLogx(1);
 }


  // -------------------------------------------------------------------------------------
  void OutputManagerPSAllInOne::logy() {
    lastPad_->SetLogy(1);
    if( lastPad_ == topCan_ ) bottomPad_->SetLogy(1);
  }


  // -------------------------------------------------------------------------------------
  void OutputManagerPSAllInOne::newPage(const TString &title) {
    ps_->NewPage();
    multiCanPadIdx_ = 1;
  }


  // -------------------------------------------------------------------------------------
  void OutputManagerPSAllInOne::nextMultiPad(const TString &title) {
    if( multiCanPadIdx_ == 13 ) {
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
    if( title.Length() ) can_->SetTitle(title);
    else can_->SetTitle("");
    can_->SetLogx(0);
    can_->SetLogy(0);
    lastPad_ = can_->cd();
  }


  // -------------------------------------------------------------------------------------
  void OutputManagerPSAllInOne::nextMainPad(const TString &title) {
    ps_->NewPage();
    if( title.Length() ) topCan_->SetTitle(title);
    else topCan_->SetTitle("");
    topCan_->SetLogx(0);
    topCan_->SetLogy(0);
    lastPad_ = topCan_->cd();
  }


  // -------------------------------------------------------------------------------------
  void OutputManagerPSAllInOne::nextRatioPad() {
    if( lastPad_ == topCan_ ) {
      bottomPad_ = util::HistOps::createRatioBottomPad();
      bottomPad_->Draw();
      bottomPad_->cd();
      bottomPad_->SetLogx(0);
      bottomPad_->SetLogy(0);
    }
  }



  // -------------------------------------------------------------------------------------
  void OutputManagerPSAllInOne::saveCurrentPad(const TString &name) {
    lastPad_->Draw();
  }




  // -------------------------------------------------------------------------------------
  OutputManagerEPSSingleFiles::OutputManagerEPSSingleFiles(const TString &fileNameBase) 
    : OutputManager(fileNameBase) { }


  // -------------------------------------------------------------------------------------
  void OutputManagerEPSSingleFiles::logx() {
    lastPad_->SetLogx(1);
    if( lastPad_ == topCan_ ) bottomPad_->SetLogx(1);
 }


  // -------------------------------------------------------------------------------------
  void OutputManagerEPSSingleFiles::logy() {
    lastPad_->SetLogy(1);
    if( lastPad_ == topCan_ ) bottomPad_->SetLogy(1);
  }


  // -------------------------------------------------------------------------------------
  void OutputManagerEPSSingleFiles::nextPad(const TString &title) {
    if( title.Length() ) can_->SetTitle(title);
    else can_->SetTitle("");
    lastPad_ = can_->cd();
    can_->SetLogx(0);
    can_->SetLogy(0);
  }


  // -------------------------------------------------------------------------------------
  void OutputManagerEPSSingleFiles::nextMainPad(const TString &title) {
    if( title.Length() ) topCan_->SetTitle(title);
    else topCan_->SetTitle("");
    lastPad_ = topCan_->cd();
    topCan_->SetLogx(0);
    topCan_->SetLogy(0);
  }


  // -------------------------------------------------------------------------------------
  void OutputManagerEPSSingleFiles::nextRatioPad() {
    if( lastPad_ == topCan_ ) {
      bottomPad_ = util::HistOps::createRatioBottomPad();
      bottomPad_->Draw();
      bottomPad_->cd();
      bottomPad_->SetLogx(0);
      bottomPad_->SetLogy(0);
    }
  }


  // -------------------------------------------------------------------------------------
  void OutputManagerEPSSingleFiles::saveCurrentPad(const TString &name) {
    lastPad_->SaveAs(name,"eps");
  }
}

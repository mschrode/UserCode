// $Id: OutputManager.cc,v 1.9 2012/06/15 23:05:55 mschrode Exp $

#include "OutputManager.h"

#include <iostream>


namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  OutputManager* OutputManager::createOutputManager(OutputManager::Mode mode, const TString &fileNameBase) {
    OutputManager* out = 0;
    if( mode == OutputManager::EPSSingleFiles ) {
      out = new OutputManagerEPSSingleFiles(fileNameBase,false);
    } else if( mode == OutputManager::EPSSingleFilesPlusROOT ) {
      out = new OutputManagerEPSSingleFiles(fileNameBase,true);
    } else if( mode == OutputManager::PSAllInOne ) {
      out = new OutputManagerPSAllInOne(fileNameBase);
    } else {
      std::cerr << "ERROR in OutputManager::createOutputManager(): No Mode '" << mode << "'" << std::endl;
      exit(1);
    }

    return out;
  }


  // -------------------------------------------------------------------------------------
  bool OutputManager::isValidMode(Mode mode) {
    bool result = true;
    if( mode != PSAllInOne && mode != EPSSingleFiles && mode != EPSSingleFilesPlusROOT ) {
      result = false;
      std::cerr << "ERROR in OutputManager::isValidMode: Invalid OutputMode" << std::endl;
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  OutputManager::OutputManager(const TString &fileNameBase) {
    can_ = new TCanvas("OutputManager:Canvas","",500,500);
    can_->SetWindowSize(500+(500-can_->GetWw()),500+(500-can_->GetWh()));
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
  TH1* OutputManager::mainFrame(const TH1* h) const {
    TH1* tmp = util::HistOps::createRatioTopFrame(h);
    double yMin = 3E-2;
    double yMax = h->GetMaximum();
    tmp->GetYaxis()->SetRangeUser(yMin,yMax);

    return tmp;
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
  OutputManagerEPSSingleFiles::OutputManagerEPSSingleFiles(const TString &fileNameBase, bool rootOutput) 
    : OutputManager(fileNameBase),
      rootOutput_(rootOutput) {
    if( rootOutput_ ) outFile_ = new TFile(fileNameBase+".root","UPDATE");
  }
  
  // -------------------------------------------------------------------------------------
  OutputManagerEPSSingleFiles::~OutputManagerEPSSingleFiles() {
    if( rootOutput_ ) {
      outFile_->Close();
      delete outFile_;
    }
  }
    

  // -------------------------------------------------------------------------------------
  void OutputManagerEPSSingleFiles::logx() {
    lastPad_->SetLogx(1);
    if( lastPad_ == topCan_ ) bottomPad_->SetLogx(1);
 }


  // -------------------------------------------------------------------------------------
  void OutputManagerEPSSingleFiles::logy() {
    lastPad_->SetLogy(1);
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
    // .eps output
    lastPad_->SaveAs(name+".eps","eps");
    
    // TCanvas into .root file
    if( rootOutput_ ) {
      TString nameTmp = lastPad_->GetName();
      lastPad_->SetName(name);
      if( rootOutput_ && outFile_->IsOpen() ) {
	outFile_->WriteTObject(lastPad_);
      }
      lastPad_->SetName(nameTmp);
    }
  }


  // -------------------------------------------------------------------------------------
  void OutputManagerEPSSingleFiles::saveCanvas(TCanvas* can, const TString &name) const {
    can->SetName(name);
    can->SaveAs(name+".eps","eps");
    if( rootOutput_ && outFile_->IsOpen() ) {
      outFile_->WriteTObject(can,name);
    }
  }

  // -------------------------------------------------------------------------------------
  void OutputManagerEPSSingleFiles::saveTObject(TObject* obj, const TString &name) const {
    TObject* c = obj->Clone(name);
    if( rootOutput_ && outFile_->IsOpen() ) {
      outFile_->WriteTObject(c);
    }
  }
}

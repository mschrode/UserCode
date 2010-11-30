// $Id: fitTailsFromAsym.C,v 1.9 2010/11/29 09:20:21 mschrode Exp $

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TPaveText.h"

#include "../sampleTools/BinningAdmin.h"
#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/FileOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"

#include "globalFunctions.h"



///////////////////////////////////////////////////////////////////////////
//
//  Type definitions
//
///////////////////////////////////////////////////////////////////////////


// +++++ Class Parameters +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Parameters {
 public:
  Parameters(const sampleTools::BinningAdmin* admin, double nSigCore, double nSigTailStart, double pt3Cut, const TString &fileDataPre, const TString &fileDataSuf, const TString &fileMCPre, const TString &fileMCSuf);
  
  const sampleTools::BinningAdmin* binningAdmin() const { return admin_; }
  double nSigCore() const { return nSigCore_; };
  double nSigTailStart() const { return nSigTailStart_; }
  double pt3Cut() const { return pt3Cut_; }
  TString inFileNameData(unsigned int etaBin) const { return inFileNameData_.at(etaBin); }
  TString inFileNameMC(unsigned int etaBin) const { return inFileNameMC_.at(etaBin); }
  TString outFileNamePrefix() const { return outFileNamePrefix_; }
  TString outFileNamePrefix(unsigned int etaBin) const { return outFileNamePrefix_+"Eta"+util::toTString(etaBin)+"_"; }
  TString outFileNamePrefix(unsigned int etaBin, unsigned int ptBin) const { return outFileNamePrefix_+"Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin)+"_"; }
  TString labelJetAlgo() const { return labelJetAlgo_; }
  TString labelDeltaPhiCut() const { return labelDeltaPhiCut_; }
  TString labelPSoftCut() const { return labelPSoftCut_; }
  int fillColor() const { return (labelJetAlgo().Contains("PF")) ? 38 : 5; }


 private:
  const sampleTools::BinningAdmin* admin_;
  const double nSigCore_;
  const double nSigTailStart_;
  const double pt3Cut_;

  std::vector<TString> inFileNameData_;
  std::vector<TString> inFileNameMC_;
  TString outFileNamePrefix_;
  TString labelJetAlgo_;
  TString labelDeltaPhiCut_;
  TString labelPSoftCut_;
};


// ------------------------------------------------------------------------------------
Parameters::Parameters(const sampleTools::BinningAdmin* admin, double nSigCore, double nSigTailStart, double pt3Cut, const TString &fileDataPre, const TString &fileDataSuf, const TString &fileMCPre, const TString &fileMCSuf)
  : admin_(admin), nSigCore_(nSigCore), nSigTailStart_(nSigTailStart), pt3Cut_(pt3Cut) {

  for(unsigned int etaBin = 0; etaBin < admin->nEtaBins(); ++etaBin) {
    inFileNameData_.push_back(fileDataPre+"_Eta"+util::toTString(etaBin)+"_"+fileDataSuf);
    inFileNameMC_.push_back(fileMCPre+"_Eta"+util::toTString(etaBin)+"_"+fileMCSuf);
  }
    
  // Prefix for output files
  outFileNamePrefix_ = "Tails_nSCore"+util::toTString(10.*nSigCore_)+"_nSTail"+util::toTString(10.*nSigTailStart_)+"_ptSoft"+util::toTString(100.*pt3Cut_)+"_";
  if( inFileNameData_.at(0).Contains("Calo") && inFileNameMC_.at(0).Contains("Calo") ) {
    outFileNamePrefix_ += "Calo_";
    labelJetAlgo_ = "AK5 Calo-Jets";
  } else if( inFileNameData_.at(0).Contains("PF") && inFileNameMC_.at(0).Contains("PF") ) {
    outFileNamePrefix_ += "PF_";  
    labelJetAlgo_ = "AK5 PF-Jets";
  } else {
    std::cerr << "WARNING in Parameters: unknown or inconsistent jet algorithms for data and MC" << std::endl;
  }

  labelPSoftCut_ = "p^{L2L3}_{T,3} < "+util::toTString(pt3Cut_)+"#upoint#bar{p^{ave}_{T}}";

  labelDeltaPhiCut_ = "|#Delta#phi| > 2.7";
}



// +++++ Class Label ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class Label {
public:
  Label(const Parameters* par);
  
  TPaveText* info(unsigned int etaBin, unsigned int ptBin) const;
  TPaveText* info(unsigned int etaBin) const;
  TPaveText* eta(unsigned int etaBin) const;
  TLegend* legend(const TH1* hData, const TH1* hMC, const TF1* tailData = 0, const TF1* tailMC = 0) const;
  

private:
  const Parameters* par_;
};


// ------------------------------------------------------------------------------------
 Label::Label(const Parameters* par) : par_(par) {}


// ------------------------------------------------------------------------------------
TPaveText* Label::info(unsigned int etaBin, unsigned int ptBin) const {
  double lumi = par_->binningAdmin()->hltLumi(etaBin,ptBin);
  double min = par_->binningAdmin()->ptMin(etaBin,ptBin);
  double max = par_->binningAdmin()->ptMax(etaBin,ptBin);

  TPaveText *txt = util::LabelFactory::createPaveText(3,-0.68);
  if( lumi > 0. ) txt->AddText("L = "+util::toTString(lumi)+" pb^{-1},  "+par_->labelJetAlgo()+",");
  else txt->AddText(par_->labelJetAlgo()+",");
  txt->AddText(par_->labelDeltaPhiCut()+",  "+par_->labelPSoftCut());
  txt->AddText(util::toTString(min)+" < p^{ave}_{T} < "+util::toTString(max)+" GeV");

  return txt;
}


// ------------------------------------------------------------------------------------
TPaveText* Label::info(unsigned int etaBin) const {
  double min = par_->binningAdmin()->ptMin(etaBin);
  double max = par_->binningAdmin()->ptMax(etaBin);

  TPaveText *txt = util::LabelFactory::createPaveText(3,-0.68);
  txt->AddText(par_->labelJetAlgo()+",");
  txt->AddText(par_->labelDeltaPhiCut()+",  "+par_->labelPSoftCut()+"#upoint#bar{p^{ave}_{T}}");
  txt->AddText(util::toTString(min)+" < p^{ave}_{T} < "+util::toTString(max)+" GeV");

  return txt;
}


// ------------------------------------------------------------------------------------
TPaveText* Label::eta(unsigned int etaBin) const {
  double min = par_->binningAdmin()->etaMin(etaBin);
  double max = par_->binningAdmin()->etaMax(etaBin);

  TPaveText *txt = util::LabelFactory::createPaveText(1,0.32);
  txt->AddText(util::toTString(min)+" < |#eta| < "+util::toTString(max));

  return txt;
}


// ------------------------------------------------------------------------------------
TLegend* Label::legend(const TH1* hData, const TH1* hMC, const TF1* tailData, const TF1* tailMC) const {
  TLegend *leg = 0;
  if( tailData && tailMC ) {
    leg = util::LabelFactory::createLegendColWithOffset(4,0.3,1);
    leg->AddEntry(hData,"Data","P");
    leg->AddEntry(tailData,"Fit Data","L");
    leg->AddEntry(hMC,"MC","F");
    leg->AddEntry(tailMC,"Fit MC","L");
  } else {
    leg = util::LabelFactory::createLegendColWithOffset(2,0.3,1);
    leg->AddEntry(hData,"Data","P");
    leg->AddEntry(hMC,"MC","F");
  }
  return leg;
}




// +++++ Class AsymmetryBin +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class AsymmetryBin {
public:
  AsymmetryBin(const Parameters* par, int type, unsigned int etaBin, unsigned int ptBin);
  ~AsymmetryBin();

  int type() const { return type_; }
  
  TH1 *hPtAsym() const { return hist(hPtAsym_); }
  TH1 *hPtAsymSmeared() const { return hist(hPtAsymSmeared_); }
  TH1 *hTailClean() const { return hist(hTailClean_); }
  double tailCleanStart() const;
  TF1 *fTail() const;
  
  double widthAsym() const { return widthAsym_; }
  double widthAsymErr() const { return widthAsymErr_; }  
  double widthAsymSmeared() const { return widthAsymSmeared_; }
  double widthAsymErrSmeared() const { return widthAsymErrSmeared_; }  

  double nTotal() const { return nTotal_; }
  double nTail() const { return nTail_; }

  void smearMCAsymmetry(const AsymmetryBin* dataAsymmetryBin);


private:
  const Parameters* par_;
  const int type_;

  TH1* hPtAsym_;
  TH1* hPtAsymSmeared_;
  TH1* hTail_;
  TH1* hTailClean_;
  TF1* fTailFit_;

  double widthAsym_;
  double widthAsymErr_;
  double widthAsymSmeared_;
  double widthAsymErrSmeared_;
  double nTotal_;
  double nTail_;

  mutable unsigned int count_;

  void fitWidth(const TH1* h, double &width, double &widthErr) const {
    func::fitCoreWidth(h,par_->nSigCore(),width,widthErr);
  }
  void getTail(const TH1* h, const TF1* fGauss = 0, double tailStart = 0);
  TH1* hist(const TH1* h) const;
  void setStyle(TH1 *h) const;
  void smearAsymmetry(double scaling, const TH1* hOrig, TH1* &hSmeared) const;
};


// ------------------------------------------------------------------------------------
AsymmetryBin::AsymmetryBin(const Parameters* par, int type, unsigned int etaBin, unsigned int ptBin)
  : par_(par), type_(type) {
  count_ = 0;

  hPtAsymSmeared_ = 0;
  hTail_ = 0;
  hTailClean_ = 0;
  fTailFit_ = 0;


  // Read histograms from file
  TString fileName;
  if( type_ == 0 ) fileName = par_->inFileNameData(etaBin);
  else if( type_ == 1 ) fileName = par_->inFileNameMC(etaBin);
  TString histName = "hPtAsym_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin);
  TString newHistName = "AsymmetryBin";
  if( type_ == 0 ) newHistName += "_Data_";
  else if( type_ == 1 ) newHistName += "_MC_";
  newHistName += histName;
  hPtAsym_ = util::FileOps::readTH1(fileName,histName,newHistName);
  if( hPtAsym_ ) {
    hPtAsym_->GetXaxis()->SetRangeUser(-1.,1);
    util::HistOps::setAxisTitles(hPtAsym_,"Asymmetry","","events",true);      
    hPtAsym_->SetTitle("");
    setStyle(hPtAsym_);
    nTotal_ = hPtAsym_->GetEntries();
  } else {
    std::cerr << "ERROR reading histogram '" << histName << "' from file '" << fileName << "'\n";
    exit(1);
  }

  // Fitting width of asymmetry
  fitWidth(hPtAsym_,widthAsym_,widthAsymErr_);
  // Extract tails from measured asymmetry
  //  getTail(hPtAsym_);


  // Initialize smeared histograms
//   smearAsymmetry(0.1,hPtAsym_,hPtAsymSmeared_);
//   fitWidth(hPtAsymSmeared_,widthAsymSmeared_,widthAsymErrSmeared_);
//   getTail(hPtAsymSmeared_);


  if( type_ == 0 ) {
    // Extract tails from measured asymmetry
    getTail(hPtAsym_);
  } else {
    // Smearing asymmetry
    smearAsymmetry(0.1,hPtAsym_,hPtAsymSmeared_);
    fitWidth(hPtAsymSmeared_,widthAsymSmeared_,widthAsymErrSmeared_);

    // Extract tails from smeared asymmetry
    getTail(hPtAsymSmeared_);
  }

}


// ------------------------------------------------------------------------------------
AsymmetryBin::~AsymmetryBin() {
  delete hPtAsym_;
  delete hPtAsymSmeared_;
  delete fTailFit_;
  delete hTail_;
  delete hTailClean_;
}


// ------------------------------------------------------------------------------------
double AsymmetryBin::tailCleanStart() const {
  double start = 1.;
  for(int bin = hTailClean_->FindBin(0); bin <= hTailClean_->GetNbinsX(); ++bin) {
    if( hTailClean_->GetBinContent(bin) > 0. ) {
      start = hTailClean_->GetBinCenter(bin);
      break;
    }
  }
  return start;
}


// ------------------------------------------------------------------------------------
TH1* AsymmetryBin::hist(const TH1* h) const {
  ++count_;
  TString name = h->GetName();
  name += count_;
  return static_cast<TH1*>(h->Clone(name));
}


// ------------------------------------------------------------------------------------
TF1* AsymmetryBin::fTail() const {
  ++count_;
  TString name = fTailFit_->GetName();
  name += count_;
  return static_cast<TF1*>(fTailFit_->Clone(name));
}


// ------------------------------------------------------------------------------------
void AsymmetryBin::smearMCAsymmetry(const AsymmetryBin* dataAsymmetryBin) {

  // Width ratio data / MC
  double scale = dataAsymmetryBin->widthAsym()/widthAsym() - 1.;

  // Smearing asymmetry
  smearAsymmetry(scale,hPtAsym_,hPtAsymSmeared_);
  fitWidth(hPtAsymSmeared_,widthAsymSmeared_,widthAsymErrSmeared_);
  
  // Extract tails from smeared asymmetry
  getTail(hPtAsymSmeared_,dataAsymmetryBin->fTail(),dataAsymmetryBin->tailCleanStart());
}


// ------------------------------------------------------------------------------------
void AsymmetryBin::smearAsymmetry(double scaling, const TH1* hOrig, TH1* &hSmeared) const {
  if( hSmeared ) delete hSmeared;
  func::smearHistogram(hOrig,hSmeared,nTotal(),widthAsym(),scaling);
  setStyle(hSmeared);
}
 
 
// ------------------------------------------------------------------------------------
void AsymmetryBin::getTail(const TH1* h, const TF1* fGauss, double tailStart) {
  double tmpResult = 0.;

  if( hTail_ ) delete hTail_;
  if( hTailClean_ ) delete hTailClean_;
  if( fTailFit_ ) delete fTailFit_;

  if( fGauss && tailStart ) {
    tmpResult = func::getTailFromGauss(h,fGauss,tailStart,par_->nSigCore(),hTail_,hTailClean_,fTailFit_);
  } else {
    tmpResult = func::getTail(h,par_->nSigCore(),par_->nSigTailStart(),hTail_,hTailClean_,fTailFit_);
  }
  if( tmpResult < 0 ) {
    std::cerr << "ERROR in AsymmetryBin::getTail(): Fitting of tail did not work.\n";
    nTail_ = 0;
  } else {
    nTail_ = tmpResult;
  }

  if( type() == 1 ) {
    fTailFit_->SetLineStyle(2);
  }
}


// ------------------------------------------------------------------------------------
void AsymmetryBin::setStyle(TH1 *h) const {
  h->SetLineWidth(1);
  if( type() == 0 ) {
    h->SetMarkerStyle(20);
  } else if( type() == 1 ) {
    h->SetMarkerStyle(1);
    h->SetFillColor(par_->fillColor());
  }
}




// +++++ Class AsymmetryAdmin +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class AsymmetryAdmin {
public:
  AsymmetryAdmin(const Parameters* par, unsigned int etaBin, unsigned int ptBin);
  ~AsymmetryAdmin();

  double widthAsymData() const { return data_->widthAsym(); }
  double widthAsymErrData() const { return data_->widthAsymErr(); }  
  double widthAsymMC() const { return mc_->widthAsym(); }
  double widthAsymErrMC() const { return mc_->widthAsymErr(); }  
  double widthSmearedAsymMC() const { return mc_->widthAsymSmeared(); }
  double widthSmearedAsymErrMC() const { return mc_->widthAsymErrSmeared(); }  
  double nTailData() const { return data_->nTail(); }
  double nTailMC() const { return mc_->nTail(); }
  double nTotalData() const { return data_->nTotal(); }
  double nTotalMC() const { return mc_->nTotal(); }

  void plotAsymmetry() const;
  void plotSmearedAsymmetry() const;
  void plotTails() const;
  
private:
  const Parameters* par_;
  const unsigned int etaBin_;
  const unsigned int ptBin_;

  Label* label_;
  AsymmetryBin* data_;
  AsymmetryBin* mc_;

  TString name_;
  TString outFileNamePrefix_;
};

typedef std::vector< std::vector<AsymmetryAdmin*> >::const_iterator AsymEtaPtIt;
typedef std::vector<AsymmetryAdmin*>::const_iterator AsymPtIt;


// ------------------------------------------------------------------------------------
AsymmetryAdmin::AsymmetryAdmin(const Parameters* par, unsigned int etaBin, unsigned int ptBin)
  : par_(par), etaBin_(etaBin), ptBin_(ptBin) {

  outFileNamePrefix_ = par_->outFileNamePrefix(etaBin,ptBin);
  std::cout << "Creating samples '" << outFileNamePrefix_ << "'" << std::endl;
  
  // Create samples holding histograms
  data_ = new AsymmetryBin(par_,0,etaBin,ptBin);
  mc_ = new AsymmetryBin(par_,1,etaBin,ptBin);
  mc_->smearMCAsymmetry(data_);

  label_ = new Label(par); 
}


// ------------------------------------------------------------------------------------
AsymmetryAdmin::~AsymmetryAdmin() {
  delete label_;
  delete data_;
  delete mc_;
}


// ------------------------------------------------------------------------------------
void AsymmetryAdmin::plotAsymmetry() const {
  TH1 *hMC = mc_->hPtAsym();
  TH1 *hData = data_->hPtAsym();

  // Log scale
  util::HistOps::setYRange(hMC,3,3E-5);
  util::HistOps::setYRange(hData,3,3E-5);
  hMC->GetXaxis()->SetRangeUser(-1.,1.);

  TString canName = outFileNamePrefix_+"PtAsym";
  TCanvas *can = new TCanvas(canName,canName,500,500);
  can->cd();
  hMC->Draw("HISTE");
  hData->Draw("PE1same");
  label_->info(etaBin_,ptBin_)->Draw("same");
  label_->eta(etaBin_)->Draw("same");
  label_->legend(hData,hMC)->Draw("same");
  can->SetLogy(1);
  can->SaveAs(outFileNamePrefix_+"PtAsymLog.eps","eps");

  // Linear scale
  util::HistOps::setYRange(hMC,3);
  hMC->GetXaxis()->SetRangeUser(-0.4,0.4);
  can->cd();
  hMC->Draw("HISTE");
  hData->Draw("PE1same");
  label_->info(etaBin_,ptBin_)->Draw("same");
  label_->eta(etaBin_)->Draw("same");
  label_->legend(hData,hMC)->Draw("same");
  can->SetLogy(0);
  can->SaveAs(outFileNamePrefix_+"PtAsym.eps","eps");
  
  // Superimpose fit
  TF1 *tailMC = mc_->fTail();
  TF1 *tailData = data_->fTail();
  can->cd();
  hMC->Draw("HISTE");
  hData->Draw("PE1same");
  tailMC->Draw("same");
  tailData->Draw("same");
  label_->info(etaBin_,ptBin_)->Draw("same");
  label_->eta(etaBin_)->Draw("same");
  label_->legend(hData,hMC,tailData,tailMC)->Draw("same");
  can->SetLogy(0);
  can->SaveAs(outFileNamePrefix_+"PtAsymFit.eps","eps");

  // Ratio data / MC
  TH1 *hRatio = util::HistOps::createRatioPlot(hData,hMC);
  TH1 *hRatioFrame = util::HistOps::createRatioFrame(hData,"Data / MC",0.,3.);
  can->cd();
  hRatioFrame->Draw();
  hRatio->Draw("PE1same");
  label_->info(etaBin_,ptBin_)->Draw("same");
  label_->eta(etaBin_)->Draw("same");
  can->SetLogy(0);
  can->SaveAs(outFileNamePrefix_+"PtAsymLogRatio.eps","eps");

  hRatioFrame->GetXaxis()->SetRangeUser(-0.4,0.4);
  can->cd();
  hRatioFrame->Draw();
  hRatio->Draw("PE1same");
  label_->info(etaBin_,ptBin_)->Draw("same");
  label_->eta(etaBin_)->Draw("same");
  can->SetLogy(0);
  can->SaveAs(outFileNamePrefix_+"PtAsymRatio.eps","eps");

  // Bottom ratio plot
  delete can;
  can = util::HistOps::createRatioTopCanvas();
  TPad *bRatioBottomPad = util::HistOps::createRatioBottomPad();
  TH1 *bRatioTopFrame = util::HistOps::createRatioTopHist(hMC);
  TH1 *bRatioBottomFrame = util::HistOps::createRatioBottomFrame(hMC,"Asymmetry","",0.61,1.39);
  can->cd();
  bRatioTopFrame->GetXaxis()->SetRangeUser(-0.4,0.4);
  bRatioTopFrame->GetYaxis()->SetRangeUser(0.,9.5);
  bRatioTopFrame->Draw("HISTE");
  hData->Draw("PE1same");
  label_->info(etaBin_,ptBin_)->Draw("same");
  label_->eta(etaBin_)->Draw("same");
  label_->legend(hData,hMC)->Draw("same");
  bRatioBottomPad->Draw();
  bRatioBottomPad->cd();
  bRatioBottomFrame->GetXaxis()->SetRangeUser(-0.4,0.4);
  bRatioBottomFrame->GetXaxis()->SetMoreLogLabels();
  bRatioBottomFrame->Draw();
  hRatio->Draw("PE1same");
  can->SaveAs(outFileNamePrefix_+"PtAsymBottomRatio.eps","eps");

  delete can;
}


// ------------------------------------------------------------------------------------
void AsymmetryAdmin::plotSmearedAsymmetry() const {
  TString canName = outFileNamePrefix_+"PtAsymSmear";
  TCanvas *can = new TCanvas(canName,canName,500,500);
  
  TH1 *hMC = mc_->hPtAsymSmeared();
  TH1 *hData = data_->hPtAsym();
  TF1 *tailMC = mc_->fTail();
  TF1 *tailData = data_->fTail();

  // Linear scale    
  util::HistOps::setYRange(hMC,3);
  hMC->GetXaxis()->SetRangeUser(-0.4,0.4);
  can->cd();
  hMC->Draw("HISTE");
  hData->Draw("PE1same");
  label_->info(etaBin_,ptBin_)->Draw("same");
  label_->eta(etaBin_)->Draw("same");
  label_->legend(hData,hMC)->Draw("same");
  can->SetLogy(0);
  can->SaveAs(outFileNamePrefix_+"PtSmearAsym.eps","eps");
  
  // Superimpose fit
  can->cd();
  hMC->Draw("HISTE");
  hData->Draw("PE1same");
  label_->info(etaBin_,ptBin_)->Draw("same");
  label_->eta(etaBin_)->Draw("same");
  tailMC->Draw("same");
  tailData->Draw("same");
  label_->legend(hData,hMC,tailData,tailMC)->Draw("same");
  can->SaveAs(outFileNamePrefix_+"TailFitLinear.eps","eps");

  // Log scale    
  util::HistOps::setYRange(hMC,3,3E-5);
  util::HistOps::setYRange(hData,3,3E-5);
  hMC->GetXaxis()->SetRangeUser(-1.,1.);
  can->cd();
  hMC->Draw("HISTE");
  hData->Draw("PE1same");
  label_->info(etaBin_,ptBin_)->Draw("same");
  label_->eta(etaBin_)->Draw("same");
  label_->legend(hData,hMC)->Draw("same");
  can->SetLogy(1);
  can->SaveAs(outFileNamePrefix_+"PtSmearAsymLog.eps","eps");

  // Superimpose fit
  hMC->GetXaxis()->SetRangeUser(0.,1.);
  can->cd();
  hMC->Draw("HISTE");
  hData->Draw("PE1same");
  label_->info(etaBin_,ptBin_)->Draw("same");
  label_->eta(etaBin_)->Draw("same");
  tailMC->Draw("same");
  tailData->Draw("same");
  label_->legend(hData,hMC,tailData,tailMC)->Draw("same");
  can->SetLogy(1);
  can->SaveAs(outFileNamePrefix_+"TailFit.eps","eps");
  
  // Ratio data / MC
  TH1 *hRatio = util::HistOps::createRatioPlot(hData,hMC);
  TH1 *hRatioFrame = util::HistOps::createRatioFrame(hData,"Data / MC",0.,3.);
  can->cd();
  hRatioFrame->Draw();
  hRatio->Draw("PE1same");
  label_->info(etaBin_,ptBin_)->Draw("same");
  label_->eta(etaBin_)->Draw("same");
  can->SetLogy(0);
  can->SaveAs(outFileNamePrefix_+"PtSmearAsymLogRatio.eps","eps");

  // Ratio data / MC
  hRatioFrame->GetXaxis()->SetRangeUser(-0.4,0.4);
  can->cd();
  hRatioFrame->Draw();
  hRatio->Draw("PE1same");
  label_->info(etaBin_,ptBin_)->Draw("same");
  label_->eta(etaBin_)->Draw("same");
  can->SetLogy(0);
  can->SaveAs(outFileNamePrefix_+"PtSmearAsymRatio.eps","eps");

  // Bottom ratio plot
  delete can;
  can = util::HistOps::createRatioTopCanvas();
  TPad *bRatioBottomPad = util::HistOps::createRatioBottomPad();
  TH1 *bRatioTopFrame = util::HistOps::createRatioTopHist(hMC);
  TH1 *bRatioBottomFrame = util::HistOps::createRatioBottomFrame(hMC,"Asymmetry","",0.61,1.39);
  can->cd();
  bRatioTopFrame->GetXaxis()->SetRangeUser(-0.4,0.4);
  bRatioTopFrame->GetYaxis()->SetRangeUser(0.,9.5);
  bRatioTopFrame->Draw("HISTE");
  hData->Draw("PE1same");
  label_->info(etaBin_,ptBin_)->Draw("same");
  label_->eta(etaBin_)->Draw("same");
  label_->legend(hData,hMC)->Draw("same");
  bRatioBottomPad->Draw();
  bRatioBottomPad->cd();
  bRatioBottomFrame->GetXaxis()->SetMoreLogLabels();
  bRatioBottomFrame->GetXaxis()->SetRangeUser(-0.4,0.4);
  bRatioBottomFrame->Draw();
  hRatio->Draw("PE1same");
  can->SaveAs(outFileNamePrefix_+"PtSmearAsymBottomRatio.eps","eps");

  delete can;
}


// ------------------------------------------------------------------------------------
void AsymmetryAdmin::plotTails() const {
  TH1 *hMC = mc_->hTailClean();
  TH1 *hData = data_->hTailClean();
  double min = 10.;
  double max = 0.;
  util::HistOps::findYRange(hData,min,max);
  hMC->GetYaxis()->SetRangeUser(min,2.*max);
  //util::HistOps::setYRange(hMC,5);
  hMC->GetXaxis()->SetRangeUser(0.,1.);
  TString canName = outFileNamePrefix_+"Tail";
  TCanvas *can = new TCanvas(canName,canName,500,500);
  can->cd();
  hMC->Draw("HISTE");
  hData->Draw("PE1same");
  label_->info(etaBin_,ptBin_)->Draw("same");
  label_->eta(etaBin_)->Draw("same");
  label_->legend(hData,hMC)->Draw("same");
  can->SaveAs(outFileNamePrefix_+"Tail.eps","eps");

  delete can;
}





///////////////////////////////////////////////////////////////////////////
//
//  Global functions
//
///////////////////////////////////////////////////////////////////////////


// ------------------------------------------------------------------------------------
void createSlides(const Parameters* par, const std::vector< std::vector<AsymmetryAdmin*> > &admins) {

  TString name = par->outFileNamePrefix()+"Slides.tex";
  std::ofstream oFile(name);

  TString info = "\\begin{small}\n  \\begin{itemize}\n";
  info += "    \\item "+par->labelJetAlgo()+", "+par->labelDeltaPhiCut()+", "+par->labelPSoftCut()+", "+par->labelPSoftCut()+"\n";
  info += "    \\item Core range $"+util::toTString(par->nSigCore())+"\\sigma$, tail start $"+util::toTString(par->nSigTailStart())+"\\sigma$\n";
  info += "  \\end{itemize}\n\\end{small}\n";
  
  // Loop over eta bins
  for(unsigned int etaBin = 0; etaBin < admins.size(); ++etaBin) {
    unsigned int nPtBins = admins[etaBin].size();
    
    TString etaRange = "$"+util::toTString(par->binningAdmin()->etaMin(etaBin))+" < |\\eta| < "+util::toTString(par->binningAdmin()->etaMax(etaBin))+"$";

    
    oFile << "\n\n\n% ----- Pt Asymmetry ---------------------------" << std::endl;
    unsigned int nSlides = nPtBins/3;
    if( nPtBins%3 > 0 ) nSlides++;
    unsigned int ptBin = 0;
    for(unsigned int slide = 0; slide < nSlides; ++slide) {
      oFile << "\n% --------------------------------------------------\n";
      oFile << "\\begin{frame}\n";
      oFile << "  \\frametitle{\\pt Asymmetry "+etaRange+" (" << slide+1 << "/" << nSlides << ")}\n";
      oFile << info;
      oFile << "  \\begin{columns}[T] \n";
      for(int col = 0; col < 3; ++col, ++ptBin) {
	oFile << "    \\begin{column}{0.3333\\textwidth} \n";
	oFile << "    \\centering\n";
	if( ptBin < nPtBins ) {
	  oFile << "      \\includegraphics[width=0.9\\textwidth]{" << par->outFileNamePrefix(etaBin,ptBin) << "PtAsym}\\\\ \n";
	  oFile << "      \\includegraphics[width=0.9\\textwidth]{" << par->outFileNamePrefix(etaBin,ptBin) << "PtAsymRatio}\\\\ \n";
	}
	oFile << "  \\end{column} \n";
      }
      oFile << "  \\end{columns} \n";
      oFile << "\\end{frame} \n";
    }

    oFile << "\n\n\n% ----- Pt Asymmetry Log ------------------------" << std::endl;
    nSlides = nPtBins/3;
    if( nPtBins%3 > 0 ) nSlides++;
    ptBin = 0;
    for(unsigned int slide = 0; slide < nSlides; ++slide) {
      oFile << "\n% --------------------------------------------------\n";
      oFile << "\\begin{frame}\n";
      oFile << "  \\frametitle{\\pt Asymmetry Log "+etaRange+" (" << slide+1 << "/" << nSlides << ")}\n";
      oFile << info;
      oFile << "  \\begin{columns}[T] \n";
      for(int col = 0; col < 3; ++col, ++ptBin) {
	oFile << "    \\begin{column}{0.3333\\textwidth} \n";
	oFile << "    \\centering\n";
	if( ptBin < nPtBins ) {
	  oFile << "      \\includegraphics[width=0.9\\textwidth]{" << par->outFileNamePrefix(etaBin,ptBin) << "PtAsymLog}\\\\ \n";
	  oFile << "      \\includegraphics[width=0.9\\textwidth]{" << par->outFileNamePrefix(etaBin,ptBin) << "PtAsymLogRatio}\\\\ \n";
	}
	oFile << "  \\end{column} \n";
      }
      oFile << "  \\end{columns} \n";
      oFile << "\\end{frame} \n";
    }

    oFile << "\n\n\n% ----- Smeared Pt Asymmetry ---------------------------" << std::endl;
    nSlides = nPtBins/3;
    if( nPtBins%3 > 0 ) nSlides++;
    ptBin = 0;
    for(unsigned int slide = 0; slide < nSlides; ++slide) {
      oFile << "\n% --------------------------------------------------\n";
      oFile << "\\begin{frame}\n";
      oFile << "  \\frametitle{Smeared \\pt Asymmetry "+etaRange+" (" << slide+1 << "/" << nSlides << ")}\n";
      oFile << info;
      oFile << "  \\begin{columns}[T] \n";
      for(int col = 0; col < 3; ++col, ++ptBin) {
	oFile << "    \\begin{column}{0.3333\\textwidth} \n";
	oFile << "    \\centering\n";
	if( ptBin < nPtBins ) {
	  oFile << "      \\includegraphics[width=0.9\\textwidth]{" << par->outFileNamePrefix(etaBin,ptBin) << "PtSmearAsym}\\\\ \n";
	  oFile << "      \\includegraphics[width=0.9\\textwidth]{" << par->outFileNamePrefix(etaBin,ptBin) << "PtSmearAsymRatio}\\\\ \n";
	}
	oFile << "  \\end{column} \n";
      }
      oFile << "  \\end{columns} \n";
      oFile << "\\end{frame} \n";
    }

    oFile << "\n\n\n% ----- Smeared Pt Asymmetry Log ------------------------" << std::endl;
    nSlides = nPtBins/3;
    if( nPtBins%3 > 0 ) nSlides++;
    ptBin = 0;
    for(unsigned int slide = 0; slide < nSlides; ++slide) {
      oFile << "\n% --------------------------------------------------\n";
      oFile << "\\begin{frame}\n";
      oFile << "  \\frametitle{Smeared \\pt Asymmetry Log "+etaRange+" (" << slide+1 << "/" << nSlides << ")}\n";
      oFile << info;
      oFile << "  \\begin{columns}[T] \n";
      for(int col = 0; col < 3; ++col, ++ptBin) {
	oFile << "    \\begin{column}{0.3333\\textwidth} \n";
	oFile << "    \\centering\n";
	if( ptBin < nPtBins ) {
	  oFile << "      \\includegraphics[width=0.9\\textwidth]{" << par->outFileNamePrefix(etaBin,ptBin) << "PtSmearAsymLog}\\\\ \n";
	  oFile << "      \\includegraphics[width=0.9\\textwidth]{" << par->outFileNamePrefix(etaBin,ptBin) << "PtSmearAsymLogRatio}\\\\ \n";
	}
	oFile << "  \\end{column} \n";
      }
      oFile << "  \\end{columns} \n";
      oFile << "\\end{frame} \n";
    }

    oFile << "\n\n\n% ----- Tails ---------------------------" << std::endl;
    nSlides = nPtBins/3;
    if( nPtBins%3 > 0 ) nSlides++;
    ptBin = 0;
    for(unsigned int slide = 0; slide < nSlides; ++slide) {
      oFile << "\n% --------------------------------------------------\n";
      oFile << "\\begin{frame}\n";
      oFile << "  \\frametitle{Tails "+etaRange+" (" << slide+1 << "/" << nSlides << ")}\n";
      oFile << info;
      oFile << "  \\begin{columns}[T] \n";
      for(int col = 0; col < 3; ++col, ++ptBin) {
	oFile << "    \\begin{column}{0.3333\\textwidth} \n";
	oFile << "    \\centering\n";
	if( ptBin < nPtBins ) {
	  oFile << "      \\includegraphics[width=\\textwidth]{" << par->outFileNamePrefix(etaBin,ptBin) << "TailFit}\\\\ \n";
	  oFile << "      \\includegraphics[width=\\textwidth]{" << par->outFileNamePrefix(etaBin,ptBin) << "Tail}\\\\ \n";
	}
	oFile << "  \\end{column} \n";
      }
      oFile << "  \\end{columns} \n";
      oFile << "\\end{frame} \n";
    }
  } // End of loop over eta bins
  
  oFile.close();
}





///////////////////////////////////////////////////////////////////////////
//
//  Main function
//
///////////////////////////////////////////////////////////////////////////


// ------------------------------------------------------------------------------------
void fitTailsFromAsym() {
  

  // +++++ Set up parameters +++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Setting parameters" << std::endl;
  gErrorIgnoreLevel = 1001;
  util::StyleSettings::presentation();
  sampleTools::BinningAdmin* binAdm = new sampleTools::BinningAdmin("BinningAdmin.cfg");
  Parameters* par = new Parameters(binAdm,2.,3.,0.1,"tails/Tails_Calo_Data_Pt3Cut","PtSoft0.root","tails/Tails_Calo_MCFall10_Pt3Cut","PtSoft0.root");
  // Global labels
  Label* label = new Label(par);



  // +++++ Reading histograms and fitting asymmetry ++++++++++++++++++++++++++++++

  std::cout << "Reading histograms and extracting tails" << std::endl;
  std::vector< std::vector<AsymmetryAdmin*> > asymAdms;
  for(unsigned int etaBin = 0; etaBin < binAdm->nEtaBins(); ++etaBin) {
    std::vector<AsymmetryAdmin*> tmp;
    for(unsigned int ptBin = 0; ptBin < binAdm->nPtBins(etaBin); ++ptBin) {
      tmp.push_back(new AsymmetryAdmin(par,etaBin,ptBin));
    }
    asymAdms.push_back(tmp);
  }




  // +++++ Plot spectra and asymmetry distributions ++++++++++++++++++++++++++++++

  std::cout << "Plotting asymmetry and tail distributions" << std::endl;
  for(AsymEtaPtIt etaIt = asymAdms.begin(); etaIt != asymAdms.end(); ++etaIt) {
    for(AsymPtIt ptIt = etaIt->begin(); ptIt != etaIt->end(); ++ptIt) {
      (*ptIt)->plotAsymmetry();
      (*ptIt)->plotSmearedAsymmetry();
      (*ptIt)->plotTails();
    }
  }



  
  // +++++ Common plotting objects +++++++++++++++++++++++++++++++++++++++++++++++

  TFile outFile(par->outFileNamePrefix()+".root","RECREATE");
  



  // +++++ Compare width of asymmetry ++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Plotting comparison of width of asymmetry" << std::endl;
  for(unsigned int etaBin = 0; etaBin < asymAdms.size(); ++etaBin) {
    TString suffix = "_Eta"+util::toTString(etaBin);
    TPaveText* labelInfo = label->info(etaBin);
    TPaveText* labelEta = label->eta(etaBin);

    // Width of asymmetry      
    TH1 *hWidthData = util::HistOps::createTH1D("hWidthData"+suffix,binAdm->ptBinEdges(etaBin),
						"p^{ave}_{T}","GeV","#sigma(Asymmetry)");
    hWidthData->GetXaxis()->SetMoreLogLabels();
    TH1 *hWidthMC = static_cast<TH1D*>(hWidthData->Clone("hWidthMC"+suffix));
    hWidthData->SetMarkerStyle(20);

    // Width of smeared asymmetry
    TH1 *hWidthSmearedMC = static_cast<TH1D*>(hWidthMC->Clone("hWidthSmearedMC"+suffix));
    hWidthMC->SetLineStyle(2);

    // Fill width
    int bin = 1;
    for(AsymPtIt ptIt = asymAdms[etaBin].begin(); ptIt != asymAdms[etaBin].end(); ++ptIt, ++bin) {
      hWidthData->SetBinContent(bin,(*ptIt)->widthAsymData());
      hWidthData->SetBinError(bin,(*ptIt)->widthAsymErrData());
      hWidthMC->SetBinContent(bin,(*ptIt)->widthAsymMC());
      hWidthMC->SetBinError(bin,(*ptIt)->widthAsymErrMC());
      hWidthSmearedMC->SetBinContent(bin,(*ptIt)->widthSmearedAsymMC());
      hWidthSmearedMC->SetBinError(bin,(*ptIt)->widthSmearedAsymErrMC());
    }

    // Plot asymmetry width
    TCanvas *can = new TCanvas("can","Width of Asymmetry",500,500);
    can->cd();
    hWidthMC->GetYaxis()->SetRangeUser(0.,0.4);
    hWidthMC->Draw("HISTE");
    hWidthData->Draw("PE1same");
    labelInfo->Draw("same");
    labelEta->Draw("same");
    TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,0.3,1);
    leg->AddEntry(hWidthData,"Data","P");
    leg->AddEntry(hWidthMC,"MC","L");
    leg->Draw("same");
    can->SetLogx();
    can->SaveAs(par->outFileNamePrefix(etaBin)+"WidthAsym.eps","eps");
    delete leg;

    TH1 *hRatioWidth = util::HistOps::createRatioPlot(hWidthData,hWidthMC);
    TH1 *hRatioWidthFrame = util::HistOps::createRatioFrame(hRatioWidth,"#sigma(Asymmetry) Data / MC",0.8,1.5);
    can->cd();
    hRatioWidthFrame->Draw();
    hRatioWidth->Draw("PE1same");
    labelInfo->Draw("same");
    labelEta->Draw("same");
    can->SetLogx();
    can->SaveAs(par->outFileNamePrefix(etaBin)+"WidthAsymRatio.eps","eps");

    can->cd();
    hWidthSmearedMC->GetYaxis()->SetRangeUser(0.,0.4);
    hWidthSmearedMC->Draw("HISTE");
    hWidthMC->Draw("HISTEsame");
    hWidthData->Draw("PE1same");
    labelInfo->Draw("same");
    labelEta->Draw("same");
    leg = util::LabelFactory::createLegendColWithOffset(3,0.3,1);
    leg->AddEntry(hWidthData,"Data","P");
    leg->AddEntry(hWidthMC,"MC","L");
    leg->AddEntry(hWidthSmearedMC,"MC (Smeared)","L");
    leg->Draw("same");
    can->SetLogx();
    can->SaveAs(par->outFileNamePrefix(etaBin)+"WidthSmearedAsym.eps","eps");
    delete leg;

    TH1 *hRatioWidthSmear = util::HistOps::createRatioPlot(hWidthData,hWidthSmearedMC);
    can->cd();
    hRatioWidthFrame->Draw();
    hRatioWidth->SetMarkerStyle(24);
    hRatioWidth->Draw("PE1same");
    hRatioWidthSmear->Draw("PE1same");
    labelInfo->Draw("same");
    labelEta->Draw("same");
    leg = util::LabelFactory::createLegendColWithOffset(2,0.3,1);
    leg->AddEntry(hRatioWidth,"Regular MC","P");
    leg->AddEntry(hRatioWidthSmear,"Smeared MC","P");
    leg->Draw("same");
    can->SetLogx();
    can->SaveAs(par->outFileNamePrefix(etaBin)+"WidthSmearedAsymRatio.eps","eps");
    delete leg;


    // Write selected histograms to file
    outFile.WriteTObject(hRatioWidth);
    outFile.WriteTObject(hRatioWidthSmear);

    
    delete labelInfo;
    delete labelEta;
    delete can;
  } // End of loop over eta bins

  
  
  // +++++ Number of tail events and scaling factors +++++++++++++++++++++++++++++++++++

  std::cout << "Plotting number of tail events" << std::endl;
  for(unsigned int etaBin = 0; etaBin < asymAdms.size(); ++etaBin) {
    TString suffix = "_Eta"+util::toTString(etaBin);
    TPaveText* labelInfo = label->info(etaBin);
    TPaveText* labelEta = label->eta(etaBin);

    // Normal pt binning    
    TH1 *hNTailData = util::HistOps::createTH1D("hNTailData"+suffix,binAdm->ptBinEdges(etaBin),
						"p^{ave}_{T}","GeV","Number Tail Events");
    hNTailData->GetXaxis()->SetMoreLogLabels();
    hNTailData->SetMarkerStyle(20);
    TH1* hNTailMC = static_cast<TH1D*>(hNTailData->Clone("hNTailMC"+suffix));
    hNTailMC->SetFillColor(par->fillColor());

    // Integrated pt binning
    std::vector<double> ptBinEdgesInt = binAdm->ptBinEdgesInt(etaBin);
    TH1 *hNTailIntData = util::HistOps::createTH1D("hNTailIntData"+suffix,ptBinEdgesInt,
						   "p^{ave}_{T}","GeV","Number Tail Events");
    hNTailIntData->GetXaxis()->SetMoreLogLabels();
    hNTailIntData->SetMarkerStyle(20);
    TH1* hNTailIntMC = static_cast<TH1D*>(hNTailIntData->Clone("hNTailIntMC"+suffix));
    hNTailIntMC->SetFillColor(par->fillColor());

    
    // Fill number of events (from data x-sec)
    std::cout << "\n\nNUMBER OF TAILS EVENTS\n";
    int bin = 1;
    int binInt = 1;
    double nTailIntData = 0.;
    double nTailIntMC = 0.;
    for(AsymPtIt ptIt = asymAdms[etaBin].begin(); ptIt != asymAdms[etaBin].end(); ++ptIt, ++bin) {
      double nTotal = (*ptIt)->nTotalData();
      double nTailData = nTotal*(*ptIt)->nTailData();
      double nTailMC = nTotal*(*ptIt)->nTailMC();
      
      hNTailData->SetBinContent(bin,nTailData);
      hNTailData->SetBinError(bin,sqrt(nTailData));
      hNTailMC->SetBinContent(bin,nTailMC);
      hNTailMC->SetBinError(bin,sqrt(nTailMC));

      nTailIntData += nTailData;
      nTailIntMC += nTailMC;
      if( hNTailMC->GetXaxis()->GetBinUpEdge(bin) == ptBinEdgesInt.at(binInt) ) {
	hNTailIntData->SetBinContent(binInt,nTailIntData);
	hNTailIntData->SetBinError(binInt,sqrt(nTailIntData));
	hNTailIntMC->SetBinContent(binInt,nTailIntMC);
	hNTailIntMC->SetBinError(binInt,sqrt(nTailIntMC));

	binInt++;
	nTailIntData = 0;
	nTailIntMC = 0;
      }

      std::cout << bin-1 << ": \t" << hNTailData->GetBinContent(bin) << " +/- " << hNTailData->GetBinError(bin) << ",  \t" << hNTailMC->GetBinContent(bin) << " +/- " << hNTailMC->GetBinError(bin) << std::endl;
    } // End of loop over pt bins


    // Plot numbers
    TCanvas *can = new TCanvas("can","Number of tail events",500,500);
    can->cd();
    hNTailMC->GetYaxis()->SetRangeUser(0,500.);
    hNTailMC->Draw("HIST");
    hNTailData->Draw("PE1same");
    labelInfo->Draw("same");
    labelEta->Draw("same");
    TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,0.3,1);
    leg->AddEntry(hNTailData,"Data","P");
    leg->AddEntry(hNTailMC,"MC","F");
    leg->Draw("same");
    can->SetLogx();
    can->SaveAs(par->outFileNamePrefix(etaBin)+"NumTailEvts.eps","eps");

    can->cd();
    hNTailIntMC->GetYaxis()->SetRangeUser(0,500.);
    hNTailIntMC->Draw("HIST");
    hNTailIntData->Draw("PE1same");
    labelInfo->Draw("same");
    labelEta->Draw("same");
    leg->Draw("same");
    can->SetLogx();
    can->SaveAs(par->outFileNamePrefix(etaBin)+"NumTailEvtsInt.eps","eps");
    delete leg;


    // Scale factors
    TH1* hScalFac = util::HistOps::createRatioPlot(hNTailData,hNTailMC,"");
    TH1* hScalFacInt = util::HistOps::createRatioPlot(hNTailIntData,hNTailIntMC,"");
    TH1 *hScalFacFrame = util::HistOps::createRatioFrame(hScalFac,"Tail Scaling Factor",2);

    TCanvas *canScalFac = new TCanvas("canScalFac","Scaling factors",500,500);
    canScalFac->cd();
    hScalFacFrame->GetYaxis()->SetRangeUser(0.,5.);
    hScalFacFrame->GetXaxis()->SetMoreLogLabels();
    hScalFacFrame->Draw();
    hScalFac->Draw("PE1same");
    labelInfo->Draw("same");
    labelEta->Draw("same");
    canScalFac->SetLogx();
    canScalFac->SaveAs(par->outFileNamePrefix(etaBin)+"ScalingFactors.eps","eps");

    canScalFac->cd();
    hScalFacFrame->Draw();
    hScalFacInt->Draw("PE1same");
    labelInfo->Draw("same");
    labelEta->Draw("same");
    canScalFac->SetLogx();
    canScalFac->SaveAs(par->outFileNamePrefix(etaBin)+"ScalingFactorsInt.eps","eps");


    // Writ to root file
    outFile.WriteTObject(hScalFac);
    outFile.WriteTObject(hScalFacInt);

    
    delete labelInfo;
    delete labelEta;
    delete canScalFac;
    delete can;
  } // End of loop over eta bins    



  // +++++ Write slides and clean up +++++++++++++++++++++++++++++++++++++++++++++++++

  outFile.Close();
  createSlides(par,asymAdms);

  delete par;
  delete binAdm;
}


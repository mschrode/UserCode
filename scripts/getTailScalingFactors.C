// $ Id: $

#include <vector>

#include <cmath>
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

#include "globalFunctions.h"
#include "../sampleTools/BinningAdmin.h"
#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/FileOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"



// ------------------------------------------------------------------------------------
class Pt3Bin {
public:
  Pt3Bin(const TString &fileNameData, const TString &fileNameMC, unsigned int etaBin, unsigned int ptBin, unsigned int pt3Bin, double pt3Thres, double nSigCore, double coreScalingFactor, double tailWindowMin, double tailWindowMax);
  ~Pt3Bin();

  
  double coreScalingFactor() const { return coreScalingFactor_; }
  double fNData() const { return fNData_; }
  double fNMCSmeared() const { return fNMCSmeared_; }
  double fNMC() const { return fNMC_; }


private:
  const unsigned int etaBin_;
  const unsigned int ptBin_;
  const unsigned int pt3Bin_;
  const double pt3Thres_;

  double coreScalingFactor_;
  double tailWindowMin_;
  double tailWindowMax_;
  double fNMCSmeared_;
  double fNMC_;
  double fNData_;

  TH1* hAsymData_;
  TH1* hAsymMC_;
  TH1* hAsymMCSmeared_;

  TH1* readHist(const TString &fileName, const TString &id);  
  void init(double nSigCore, double coreScalingFactor, double tailWindowMin, double tailWindowMax);
};



// ------------------------------------------------------------------------------------
class EtaPtBin {
public:
  EtaPtBin(unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScalingFactor, double nSigTailWindowMin, double nSigTailWindowMax);
  ~EtaPtBin();


  unsigned int nPt3Bins() const { return pt3Bins_.size(); }

  void addPt3Bin(double thres, const TString &fileNameData, const TString &fileNameMC);


private:
  const unsigned int etaBin_;
  const unsigned int ptBin_;
  const double nSigCore_;
  const double coreScalingFactor_;
  const double nSigTailWindowMin_;
  const double nSigTailWindowMax_;
  
  double tailWindowMin_;
  double tailWindowMax_;

  std::vector<Pt3Bin*> pt3Bins_;
};




// ------------------------------------------------------------------------------------
EtaPtBin::EtaPtBin(unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScalingFactor, double nSigTailWindowMin, double nSigTailWindowMax)
  : etaBin_(etaBin), ptBin_(ptBin),
    nSigCore_(nSigCore), coreScalingFactor_(coreScalingFactor),
    nSigTailWindowMin_(nSigTailWindowMin), nSigTailWindowMax_(nSigTailWindowMax) {

  tailWindowMin_ = 0.;
  tailWindowMax_ = 0.;
}


// ------------------------------------------------------------------------------------
EtaPtBin::~EtaPtBin() {
  for(std::vector<Pt3Bin*>::iterator it = pt3Bins_.begin();
      it != pt3Bins_.end(); ++it) {
    delete *it;
  }
}


// ------------------------------------------------------------------------------------
void EtaPtBin::addPt3Bin(double thres, const TString &fileNameData, const TString &fileNameMC) {
  // If has not been done yet, find asymmetry corresponding
  // to tail window borders from Gaussian fit to MC asymmetry
  // in this pt3 bin
  if( tailWindowMin_ == tailWindowMax_ ) {
    TString histName = "hPtAsym_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_);
    TH1* h = util::FileOps::readTH1(fileNameMC,histName);
    h->GetXaxis()->SetRangeUser(-1.,1);
    double width = 0.;
    double widthErr = 1000.;
    if( util::HistOps::fitCoreWidth(h,nSigCore_,width,widthErr) ) {
      tailWindowMin_ = std::abs(nSigTailWindowMin_*width);
      tailWindowMax_ = std::abs(nSigTailWindowMax_*width);
    }
    delete h;
  }

  // Add pt3 bin
  pt3Bins_.push_back(new Pt3Bin(fileNameData,fileNameMC,etaBin_,ptBin_,pt3Bins_.size(),thres,nSigCore_,coreScalingFactor_,tailWindowMin_,tailWindowMax_));
}


// ------------------------------------------------------------------------------------
Pt3Bin::Pt3Bin(const TString &fileNameData, const TString &fileNameMC, unsigned int etaBin, unsigned int ptBin, unsigned int pt3Bin, double pt3Thres, double nSigCore, double coreScalingFactor, double tailWindowMin, double tailWindowMax)
  : etaBin_(etaBin), ptBin_(ptBin), pt3Bin_(pt3Bin), pt3Thres_(pt3Thres),
    coreScalingFactor_(coreScalingFactor),
    tailWindowMin_(tailWindowMin),
    tailWindowMax_(tailWindowMax) {

  hAsymData_ = readHist(fileNameData,"Data");
  hAsymMC_ = readHist(fileNameMC,"MC");
  init(nSigCore,coreScalingFactor_,tailWindowMin_,tailWindowMax_);
}


// ------------------------------------------------------------------------------------
Pt3Bin::~Pt3Bin() {
  delete hAsymData_;
  delete hAsymMC_;
  delete hAsymMCSmeared_;
}


// Get asymmetry histograms from Kalibri input files
// ------------------------------------------------------------------------------------
TH1* Pt3Bin::readHist(const TString &fileName, const TString &id) {

  TString histName = "hPtAsym_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_);
  TH1* h = util::FileOps::readTH1(fileName,histName,id+histName+"_Pt3Bin"+util::toTString(pt3Bin_));
  h->GetXaxis()->SetRangeUser(-1.,1);
  util::HistOps::setAxisTitles(h,"Asymmetry","","events",true);      
  h->SetTitle("");

  return h;
}


// ------------------------------------------------------------------------------------
void Pt3Bin::init(double nSigCore, double coreScalingFactor, double tailWindowMin, double tailWindowMax) {

  // Smear MC asymmetry
  double width = 0.;
  double widthErr = 1000.;
  if( util::HistOps::fitCoreWidth(hAsymMC_,nSigCore,width,widthErr) ) {
    func::smearHistogram(hAsymMC_,hAsymMCSmeared_,hAsymMC_->GetEntries(),width,coreScalingFactor);
  } else {
    func::smearHistogram(hAsymMC_,hAsymMCSmeared_,hAsymMC_->GetEntries(),0.,0.);
  }

  // Get relative number of entries in tail
  int startBin = hAsymMCSmeared_->FindBin(tailWindowMin);
  int endBin = hAsymMCSmeared_->FindBin(tailWindowMax);
  fNMCSmeared_ = hAsymMCSmeared_->Integral(startBin,endBin) / hAsymMCSmeared_->Integral();
  fNMC_ = hAsymMC_->Integral(startBin,endBin) / hAsymMC_->Integral();
  fNData_ = hAsymData_->Integral(startBin,endBin) / hAsymData_->Integral();
}

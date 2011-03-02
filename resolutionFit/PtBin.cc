// $ Id: $

#include <iostream>

#include "PtBin.h"


namespace resolutionFit{

  // -------------------------------------------------------------------------------------
  PtBin::PtBin(unsigned int etaBin, unsigned int ptBin, const Parameters* par)
    : par_(par), etaBin_(etaBin), ptBin_(ptBin) {
    
    mcTruthSample_ = 0;
  }


  // -------------------------------------------------------------------------------------
  PtBin::~PtBin() {
    for(std::map<SampleLabel,Sample*>::iterator it = samples_.begin(); it != samples_.end(); ++it) {
      delete it->second;
    }
    if( hasMCTruthSample() ) delete mcTruthSample_;
  }


  // -------------------------------------------------------------------------------------
  const Sample* PtBin::findSample(SampleLabel label) const {
    SampleIt sIt = samples_.find(label);
    if( sIt == samples_.end() ) {
      std::cerr << "ERROR in PtBin::findSample(): No sample '" << label << "'" << std::endl;
      exit(1);
    }

    return sIt->second;
  }


  // -------------------------------------------------------------------------------------
  bool PtBin::addSample(Sample::Type type, const TString &label, const TString &baseFileName) {
    bool result = true;
    std::vector<TString> fileNames;
    for(unsigned int ptSoftBin = 0; ptSoftBin < par_->nPtSoftBins(); ++ptSoftBin) {
      fileNames.push_back(baseFileName+par_->fileNameSuffix(etaBin(),ptBin(),ptSoftBin));
    }
    if( type == Sample::Data ) {
      DataSample* s = new DataSample(label,fileNames,par_->ptSoft(),par_->histNameSuffix(etaBin(),ptBin()),par_->verbosity());
      dataSamples_[label] = s;
      samples_[label] = s;
    } else if( type == Sample::MC ) {
      MCSample* s = new MCSample(label,fileNames,par_->ptSoft(),par_->histNameSuffix(etaBin(),ptBin()),par_->verbosity());
      mcSamples_[label] = s;
      samples_[label] = s;
    } else {
      result = false;
      if( par_->verbosity() ) std::cerr << "ERROR in PtBin::addSample: Unknown Sample::Type '" << type << "'" << std::endl;
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  bool PtBin::addMCTruthSample(const TString &label, const TString &baseFileName) {
    std::vector<TString> fileNames;
    for(unsigned int ptSoftBin = 0; ptSoftBin < par_->nPtSoftBins(); ++ptSoftBin) {
      fileNames.push_back(baseFileName+par_->fileNameSuffix(etaBin(),ptBin(),ptSoftBin));
    }
    mcTruthSample_ = new MCTruthSample(label,fileNames,par_->ptSoft(),par_->histNameSuffix(etaBin(),ptBin()),par_->verbosity());
    mcTruthSample_->addFitResult(FitResult::PtGenAsym);
    
    return true;
  }


  // -------------------------------------------------------------------------------------
  bool PtBin::addFitResult(FitResult::Type type) {
    bool result = true;
    for(std::map<SampleLabel,Sample*>::iterator it = samples_.begin(); it != samples_.end(); ++it) {
      if( it->second->type() != Sample::MCTruth ) {
	result = result && it->second->addFitResult(type);
      }      
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  bool PtBin::setKSoftFit(const SampleLabel &label, FitResult::Type type, const TF1* fit) {
    bool result = true;
    SampleIt it = samples_.find(label);
    if( it != samples_.end() ) {
      result = it->second->setKSoftFit(type,fit);
    } else {
      result = false;
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  TString PtBin::toTString() const {
    TString str = "(Eta ";
    str += etaBin(); 
    str += ", Pt ";
    str += ptBin();
    str += ")";

    return str;
  }
}

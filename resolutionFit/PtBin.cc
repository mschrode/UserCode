// $ Id: $

#include <iostream>

#include "TH1.h"

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
  const Sample* PtBin::findSample(const SampleLabel &label) const {
    SampleIt sIt = samples_.find(label);
    if( sIt == samples_.end() ) {
      std::cerr << "ERROR in PtBin::findSample(): No sample '" << label << "'" << std::endl;
      exit(1);
    }

    return sIt->second;
  }


  // -------------------------------------------------------------------------------------
  bool PtBin::addSample(Sample::Type type, const TString &label, const TString &baseFileName, const TString &baseFileNameSpectrum) {
    bool result = true;
    std::vector<TString> fileNames;
    for(unsigned int ptSoftBin = 0; ptSoftBin < par_->nPtSoftBins(); ++ptSoftBin) {
      fileNames.push_back(baseFileName+par_->fileNameSuffix(etaBin(),ptBin(),ptSoftBin));
    }
    if( type == Sample::Data ) {
      DataSample* s = new DataSample(label,par_->ptMin(etaBin(),ptBin()),par_->ptMax(etaBin(),ptBin()),
				     fileNames,par_->ptSoft(),par_->histNameSuffix(etaBin(),ptBin()),
				     baseFileNameSpectrum+par_->fileNameSuffix(etaBin()),
				     par_->verbosity());
      dataSamples_[label] = s;
      samples_[label] = s;
    } else if( type == Sample::MC ) {
      MCSample* s = new MCSample(label,par_->ptMin(etaBin(),ptBin()),par_->ptMax(etaBin(),ptBin()),
				 fileNames,par_->ptSoft(),par_->histNameSuffix(etaBin(),ptBin()),
				 baseFileNameSpectrum+par_->fileNameSuffix(etaBin()),
				 par_->verbosity());
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
    mcTruthSample_ = new MCTruthSample(label,par_->ptMin(etaBin(),ptBin()),
				       par_->ptMax(etaBin(),ptBin()),
				       fileNames,par_->ptSoft(),
				       par_->histNameSuffix(etaBin(),ptBin()),
				       par_->verbosity());
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


  // Compute sample weight for 'sampleToBeWeighted' to weight it to 
  // 'referenceSample'. The weights are determined by comparison of
  // the ptAve distributions.
  // Note: This has to be done per ptAve bin as different bins contain
  // data from a different number triggers.
  // -------------------------------------------------------------------------------------
  void PtBin::weightSampleRelTo(const SampleLabel &sampleToBeWeighted, const SampleLabel &referenceSample) {
    const Sample* sReference = findSample(referenceSample);
    std::map<SampleLabel,Sample*>::iterator it = samples_.find(sampleToBeWeighted);
    if( it == samples_.end() ) {
      std::cerr << "ERROR in PtBin::weightSampleRelTo(): No sample '" << sampleToBeWeighted << "'" << std::endl;
      exit(1);
    }
    Sample* sToBeWeighted = it->second;
    if( Sample::type(sampleToBeWeighted) == Sample::Data ) {
      std::cerr << "\n\n**** WARNING: Weighting data sample '" << sToBeWeighted->label() << "' ****\n" << std::endl;
    }

    if( sToBeWeighted && sReference ) {
      for(unsigned int ptSoftBin = 0; ptSoftBin < sToBeWeighted->nPtSoftBins(); ++ptSoftBin) {
	// Get entries in ptAve distributions per sample
	TH1* h = sToBeWeighted->histPtAve(ptSoftBin);
	double entriesToBeWeighted = h->Integral();
	delete h;
	h = sReference->histPtAve(ptSoftBin);
	double entriesReference = h->Integral();
	delete h;
	double relWeight = entriesToBeWeighted > 0. ? entriesReference/entriesToBeWeighted : 1.;
	sToBeWeighted->setRelativeWeightTo(referenceSample,ptSoftBin,relWeight);
      }
    }
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

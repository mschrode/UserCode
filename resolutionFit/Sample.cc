#include "Sample.h"

#include <cassert>
#include <iostream>

#include "../util/HistOps.h"


namespace resolutionFit {

  std::map<TString,int> Sample::COLOR;
  std::map<TString,int> Sample::MARKER_STYLE;
  std::map<TString,Sample::Type> Sample::TYPE;
  unsigned int DataSample::N_DATA_SAMPLES = 0;
  unsigned int MCSample::N_MC_SAMPLES = 0;
  unsigned int MCTruthSample::N_MCTRUTH_SAMPLES = 0;


  // -------------------------------------------------------------------
  bool Sample::validType(Type type) {
    if( type != Data && type != MC && type != MCTruth ) {
      std::cerr << "ERROR in Sample::validType(): Unknown type '" << type << "'" << std::endl;
      exit(1);
    }
    
    return true;
  }


  // -------------------------------------------------------------------
  TString Sample::toString(Type type) {
    TString str = "";
    if( validType(type) ) {
      if( type == Data ) str = "Data";
      else if( type == MC ) str = "MC";
      else if( type == MCTruth ) str = "MC Truth";
    }

    return str;
  }


  // -------------------------------------------------------------------
  int Sample::color(const TString &label) {
    int c = 0;
    std::map<TString,int>::const_iterator it = COLOR.find(label);
    if( it != COLOR.end() ) c = it->second;
    
    return c;
  }


  // -------------------------------------------------------------------
  int Sample::markerStyle(const TString &label) {
    int s = 0;
    std::map<TString,int>::const_iterator it = MARKER_STYLE.find(label);
    if( it != MARKER_STYLE.end() ) s = it->second;
    
    return s;
  }


  // -------------------------------------------------------------------
  Sample::Type Sample::type(const TString &label) {
    Sample::Type t = NONE;
    std::map<TString,Type>::const_iterator it = TYPE.find(label);
    if( it != TYPE.end() ) t = it->second;
    
    return t;
  }
  


  // -------------------------------------------------------------------
  Sample::Sample(const TString &label, unsigned int verbosity)
    : label_(label), verbosity_(verbosity) {
  }


  // -------------------------------------------------------------------
  Sample::~Sample() {
    for(std::vector<Measurement*>::iterator it = meas_.begin(); it != meas_.end(); ++it) {
      delete *it;
    }
    for(std::map<FitResult::Type,FitResult*>::iterator it = fitResult_.begin(); it != fitResult_.end(); ++it) {
      delete it->second;
    }
  }



  // -------------------------------------------------------------------
  double Sample::relativeWeightTo(const SampleLabel &other, unsigned int ptSoftBin) const {
    assert( ptSoftBin < nPtSoftBins() );
    double relWeight = 1.;
    std::map<SampleLabel,std::vector<double> >::const_iterator it = relWeightToOtherSample_.find(other);
    if( it != relWeightToOtherSample_.end() ) relWeight = it->second.at(ptSoftBin);

    return relWeight;
  }



  // -------------------------------------------------------------------
  void Sample::setRelativeWeightTo(const SampleLabel &other, unsigned int ptSoftBin, double relWeight) {
    assert( ptSoftBin < nPtSoftBins() );
    std::map<SampleLabel,std::vector<double> >::iterator it = relWeightToOtherSample_.find(other);
    if( it != relWeightToOtherSample_.end() ) {
      it->second.at(ptSoftBin) = relWeight;
    } else {
      std::vector<double> weights(nPtSoftBins(),1.);
      weights.at(ptSoftBin) = relWeight;
      relWeightToOtherSample_[other] = weights;
    }
  }


  // -------------------------------------------------------------------
  bool Sample::addFitResult(FitResult::Type type, double minPt3) {
    bool result = true;
    FitResultMapIt it = fitResult_.find(type);
    if( it != fitResult_.end() ) {
      result = false;
      if( verbosity_ ) std::cerr << "ERROR in Sample: FitResult of type '" << FitResult::toString(type) << "' does already exist" << std::endl;
    } else {
      if( verbosity_ ) std::cout << "Sample::addFitResult(): Adding FitResult '" << FitResult::toString(type) << "'" << std::endl;
      fitResult_[type] = FitResult::createFitResult(type,meas_,minPt3,verbosity_);
    }

    return result;
  }


  // -------------------------------------------------------------------
  bool Sample::setKSoftFit(FitResult::Type type, const TF1* fit) {
    bool result = true;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) fr->setKSoftFit(fit);
    else result = false;

    return result;
  }


  // -------------------------------------------------------------------
  TF1* Sample::kSoftFit(FitResult::Type type, const TString &name) const {
    TF1* result = 0;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->kSoftFit(name);
    
    return result;
  }


  // -------------------------------------------------------------------
  void Sample::ptSoft(std::vector<double> &ptSoftVals) const {
    ptSoftVals.clear();
    for(unsigned int ptSoftBin = 0; ptSoftBin < nPtSoftBins(); ++ptSoftBin) {
      ptSoftVals.push_back(ptSoft(ptSoftBin));
    }
  }


  //! \brief Get y values of points which enter the extrapolation
  //! These might be absolute or relative resolutions depending on
  //! the chosen FitResultType. Note also that this is not necessarily
  //! the same as 'values()', which always returns the relative resolutions. 
  // -------------------------------------------------------------------
  void Sample::valuesInExtrapolation(FitResult::Type type, std::vector<double> &val, std::vector<double> &uncert) const {
    val.clear();
    uncert.clear();
    for(unsigned int ptSoftBin = 0; ptSoftBin < nPtSoftBins(); ++ptSoftBin) {
      val.push_back(valueInExtrapolation(type,ptSoftBin));
      uncert.push_back(uncertInExtrapolation(type,ptSoftBin));
    }
  }
  
  
  // -------------------------------------------------------------------
  void Sample::fittedValues(FitResult::Type type, std::vector<double> &val, std::vector<double> &uncert) const {
    val.clear();
    uncert.clear();
    for(unsigned int ptSoftBin = 0; ptSoftBin < nPtSoftBins(); ++ptSoftBin) {
      val.push_back(fittedValue(type,ptSoftBin));
      uncert.push_back(fittedUncert(type,ptSoftBin));
    }
  }


  // -------------------------------------------------------------------
  double Sample::meanPt(FitResult::Type type) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->meanPt();
    
    return result;
  }


  // -------------------------------------------------------------------
  double Sample::meanPtStatUncert(FitResult::Type type) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->meanPtUncert();
    
    return result;
  }


  // -------------------------------------------------------------------
  double Sample::valueInExtrapolation(FitResult::Type type, unsigned int ptSoftBin) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->valueInExtrapolation(ptSoftBin);
    
    return result;
  }


  // -------------------------------------------------------------------
  double Sample::uncertInExtrapolation(FitResult::Type type, unsigned int ptSoftBin) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->statUncertInExtrapolation(ptSoftBin);
    
    return result;
  }


  // -------------------------------------------------------------------
  double Sample::fittedValue(FitResult::Type type, unsigned int ptSoftBin) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->value(ptSoftBin);
    
    return result;
  }


  // -------------------------------------------------------------------
  double Sample::fittedUncert(FitResult::Type type, unsigned int ptSoftBin) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->statUncert(ptSoftBin);
    
    return result;
  }


  // -------------------------------------------------------------------
  double Sample::kSoftSlope(FitResult::Type type) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->kSoftSlope();
    
    return result;
  }


  // -------------------------------------------------------------------
  double Sample::kSoftSlopeStatUncert(FitResult::Type type) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->kSoftSlopeStatUncert();
    
    return result;
  }


  // -------------------------------------------------------------------
  TF1* Sample::extrapolationFunction(FitResult::Type type, const TString &name) const {
    TF1* result = 0;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->extrapolationFunction(name);
    
    return result;
  }


  // -------------------------------------------------------------------
  double Sample::extrapolatedValue(FitResult::Type type) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->extrapolatedValue();
    
    return result;
  }


  // -------------------------------------------------------------------
  double Sample::extrapolatedStatUncert(FitResult::Type type) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->extrapolatedStatUncert();
    
    return result;
  }


  // -------------------------------------------------------------------
  double Sample::extrapolatedSystUncert(FitResult::Type type) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->extrapolatedSystUncert();
    
    return result;
  }


  // -------------------------------------------------------------------
  TString Sample::labelQuantityInExtrapolation(FitResult::Type type) const {
    TString result = "";
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->labelQuantityInExtrapolation();
    
    return result;
  }


  // -------------------------------------------------------------------
  bool Sample::findFitResult(FitResult::Type type, FitResult* &fitResult) const {
    bool result = false;
    FitResultMapIt it = fitResult_.find(type);
    if( it != fitResult_.end() ) {
      result = true;
      fitResult = it->second;
    } else {
      std::cerr << "WARNING in Sample::findFitResult(): no fit result type '" << FitResult::toString(type) << "'" << std::endl;
    }

    return result;
  }


  // -------------------------------------------------------------------
  int Sample::color(unsigned int idx) const {
    int col[6] = { 1, 4, 8, 28, 38, 6 };
    return (idx>=0 && idx<6) ? col[idx] : 1;
  }


  //! Fit Gaussian to asymmetry histograms for the different ptSoft bins
  //! and store standard deviation and uncertainty
  // -------------------------------------------------------------------
  void Sample::fitAsymmetryWidths() {
    for(unsigned int ptSoftBin = 0; ptSoftBin < meas_.size(); ++ptSoftBin) {
      double width = 0.;
      double widthErr = 1000.;
      TH1* hPtAsym = meas_.at(ptSoftBin)->histPtAsym();
      hPtAsym->GetXaxis()->SetRangeUser(-1.,1.);
      if( hPtAsym->Integral(1,hPtAsym->GetNbinsX()) > 0. ) {
	if( util::HistOps::fitCoreWidth(hPtAsym,2.,width,widthErr) ) {
	  asymmetryWidths_.push_back(width);
	  asymmetryWidthErrs_.push_back(widthErr);
	} else {
	  asymmetryWidths_.push_back(0.);
	  asymmetryWidthErrs_.push_back(1000.);
	  std::cerr << "WARNING: error when fitting asymmetry of sample '" << label() << "' for " << meas_.at(ptSoftBin)->ptMin() << " < ptAve < " << meas_.at(ptSoftBin)->ptMax() << " GeV and ptSoft/ptAve < " << meas_.at(ptSoftBin)->ptSoft() << std::endl;
	}
      } else {
	asymmetryWidths_.push_back(0.);
	asymmetryWidthErrs_.push_back(1000.);
      }
      delete hPtAsym;
    }
  }



  // -------------------------------------------------------------------
  DataSample::DataSample(const TString &label, unsigned int etaBin, unsigned int ptBin, double ptMin, double ptMax, const std::vector<double> &ptSoft, const TString &fileName, unsigned int verbosity)
    : Sample(label,verbosity) {
    // Read measurements and fitted values from file
    for(unsigned int i = 0; i < ptSoft.size(); ++i) {
      meas_.push_back(new Measurement(fileName,etaBin,ptBin,i,ptMin,ptMax,ptSoft.at(i),verbosity));
    }
    fitAsymmetryWidths();

    // Style
    if( Sample::type(label) == Sample::NONE ) {
      TYPE[label] = Sample::Data;
      MARKER_STYLE[label] = 20 + N_DATA_SAMPLES;
      COLOR[label] = color(2*N_DATA_SAMPLES);
      ++N_DATA_SAMPLES;
    }
  }



  // -------------------------------------------------------------------
  MCSample::MCSample(const TString &label, unsigned int etaBin, unsigned int ptBin, double ptMin, double ptMax, const std::vector<double> &ptSoft, const TString &fileName, unsigned int verbosity)
    : Sample(label,verbosity) {

    // Read measurements and fitted values from file
    for(unsigned int i = 0; i < ptSoft.size(); ++i) {
      meas_.push_back(new Measurement(fileName,etaBin,ptBin,i,ptMin,ptMax,ptSoft.at(i),verbosity));
    }
    fitAsymmetryWidths();

    // Style
    if( Sample::type(label) == Sample::NONE ) {
      TYPE[label] = Sample::MC;
      MARKER_STYLE[label] = 25 + N_MC_SAMPLES;
      COLOR[label] = color(1+2*N_MC_SAMPLES);
      ++N_MC_SAMPLES;
    }
  }



  // -------------------------------------------------------------------
  MCTruthSample::MCTruthSample(const TString &label, unsigned int etaBin, unsigned int ptBin, double ptMin, double ptMax, const std::vector<double> &ptSoft, const TString &fileName, unsigned int verbosity)
    : Sample(label,verbosity) {

    // Read measurements from file
    for(unsigned int i = 0; i < ptSoft.size(); ++i) {
      meas_.push_back(new Measurement(fileName,etaBin,ptBin,i,ptMin,ptMax,ptSoft.at(i),verbosity));
    }

    // Style
    if( Sample::type(label) == Sample::NONE ) {
      TYPE[label] = Sample::MCTruth;
      MARKER_STYLE[label] = 24 + N_MCTRUTH_SAMPLES;
      COLOR[label] = 2;
      ++N_MCTRUTH_SAMPLES;
    }
  }
}

// $Id: EtaBin.h,v 1.2 2011/02/18 18:42:22 mschrode Exp $

#ifndef ETA_BIN_H
#define ETA_BIN_H

#include <vector>

#include "TString.h"

#include "FitResult.h"
#include "Parameters.h"
#include "PtBin.h"
#include "ResolutionFunction.h"
#include "Sample.h"


namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  class EtaBin {
  public:
    EtaBin(unsigned int etaBin, unsigned int nPtBins, const Parameters* par);
    ~EtaBin();

    unsigned int etaBin() const { return etaBin_; }
    TString toString() const;

    PtBinIt ptBinsBegin() const { return ptBins_.begin(); }
    PtBinIt ptBinsEnd() const { return ptBins_.end(); }
    const PtBin* ptBinsFront() const { return ptBins_.front(); }
    FitResultTypeIt fitResultTypesBegin() const { return fitResultTypes_.begin(); }
    FitResultTypeIt fitResultTypesEnd() const { return fitResultTypes_.end(); }
    unsigned int nSampleTypes() const { return sampleTypes_.size(); }
    SampleTypeIt sampleTypesBegin() const { return sampleTypes_.begin(); }
    SampleTypeIt sampleTypesEnd() const { return sampleTypes_.end(); }
    unsigned int nComparedSamples() const { return compSamples_.size(); }
    ComparedSamplesIt comparedSamplesBegin() const { return compSamples_.begin(); }
    ComparedSamplesIt comparedSamplesEnd() const { return compSamples_.end(); }

    Sample::Type sampleType(const SampleLabel &label) const;

    double meanPt(SampleLabel label, FitResult::Type type, unsigned int ptBinIdx) const;
    double meanPtStatUncert(SampleLabel label, FitResult::Type type, unsigned int ptBinIdx) const;
    double extrapolatedValue(SampleLabel label, FitResult::Type type, unsigned int ptBinIdx) const;
    double extrapolatedStatUncert(SampleLabel label, FitResult::Type type, unsigned int ptBinIdx) const;
    double correctedResolution(SampleLabel label, FitResult::Type type, unsigned int ptBinIdx) const;
    double correctedResolutionStatUncert(SampleLabel label, FitResult::Type type, unsigned int ptBinIdx) const;
    double pli(SampleLabel label, FitResult::Type type, unsigned int ptBinIdx) const;
    double mcTruthResolution(SampleLabel label, FitResult::Type type, unsigned int ptBinIdx) const;

    TF1* mcTruthResoFunc(const TString &name) const { return mcTruthReso_->func(name); }
    TF1* pliFunc(const TString &name) const { return pli_->func(name); }

    bool addDataSample(const TString &label, const TString &baseFileName);
    bool addMCSample(const TString &label, const TString &baseFileName);
    bool addFitResult(FitResult::Type type);

    bool compareSamples(const SampleLabel &label1, const SampleLabel &label2);

    void setMCTruthResolution(ResolutionFunction* mcTruthReso);
    void setPLI(ResolutionFunction* pli);
    

  private:
    const Parameters* par_;
    const unsigned int etaBin_;

    PtBins ptBins_;
    FitResultTypes fitResultTypes_;
    SampleTypes sampleTypes_;    
    ComparedSamples compSamples_;

    ResolutionFunction* mcTruthReso_;
    ResolutionFunction* pli_;
  };

  typedef std::vector<EtaBin*> EtaBins;
  typedef std::vector<EtaBin*>::const_iterator EtaBinIt;
}
#endif

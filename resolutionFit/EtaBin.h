// $Id: EtaBin.h,v 1.4 2011/02/25 19:50:21 mschrode Exp $

#ifndef ETA_BIN_H
#define ETA_BIN_H

#include <set>
#include <vector>

#include "TString.h"

#include "FitResult.h"
#include "Parameters.h"
#include "PtBin.h"
#include "ResolutionFunction.h"
#include "Sample.h"
#include "SystematicUncertainty.h"


namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  class EtaBin {
  public:
    EtaBin(unsigned int etaBin, unsigned int nPtBins, const Parameters* par);
    ~EtaBin();

    unsigned int etaBin() const { return etaBin_; }
    TString toString() const;
    Sample::Type sampleType(const SampleLabel &label) const;

    PtBinIt ptBinsBegin() const { return ptBins_.begin(); }
    PtBinIt ptBinsEnd() const { return ptBins_.end(); }
    const PtBin* ptBinsFront() const { return ptBins_.front(); }
    unsigned int nFitResultTypes() const { return fitResultTypes_.size(); }
    FitResultTypeIt fitResultTypesBegin() const { return fitResultTypes_.begin(); }
    FitResultTypeIt fitResultTypesEnd() const { return fitResultTypes_.end(); }
    unsigned int nSampleTypes() const { return sampleTypes_.size(); }
    SampleTypeIt sampleTypesBegin() const { return sampleTypes_.begin(); }
    SampleTypeIt sampleTypesEnd() const { return sampleTypes_.end(); }
    unsigned int nComparedSamples() const { return compSamples_.size(); }
    ComparedSamplesIt comparedSamplesBegin() const { return compSamples_.begin(); }
    ComparedSamplesIt comparedSamplesEnd() const { return compSamples_.end(); }
    bool hasSystematicUncertainty(const SampleLabel &label, FitResult::Type type) const;
    bool findSystematicUncertainty(const SampleLabel &label, FitResult::Type type, const SystematicUncertainty* &uncert) const;

    double meanPt(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    double meanPtStatUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    double extrapolatedValue(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    double extrapolatedStatUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    double extrapolatedSystUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    double correctedResolution(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx, double scalePLI = 1.) const;
    double correctedResolutionStatUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    double pli(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    double mcTruthResolution(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    double relativeMCClosure(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;

    TF1* mcTruthResoFunc(const TString &name) const { return mcTruthReso_->func(name); }
    TF1* pliFunc(const TString &name) const { return pli_->func(name); }

    bool addDataSample(const TString &label, const TString &baseFileName);
    bool addMCSample(const TString &label, const TString &baseFileName);
    bool addFitResult(FitResult::Type type);

    bool addExtrapolationUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color);
    bool addPLIUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color);
    bool addMCClosureUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color);
    bool addUncertaintyFromVariedSample(const TString &uncertaintyLabel, double fraction, const SampleLabel &nominalSample, FitResult::Type type, const TString &variedSampleDown, const TString &variedSampleUp, int color);

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
    std::set<SystematicUncertainty*> systUncerts_;

    ResolutionFunction* mcTruthReso_;
    ResolutionFunction* pli_;

    SystematicUncertainty* findSystematicUncertainty(const SampleLabel &label, FitResult::Type type);
  };

  typedef std::vector<EtaBin*> EtaBins;
  typedef std::vector<EtaBin*>::const_iterator EtaBinIt;
}
#endif

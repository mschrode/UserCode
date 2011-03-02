// $Id: EtaBin.h,v 1.7 2011/03/01 16:52:41 mschrode Exp $

#ifndef ETA_BIN_H
#define ETA_BIN_H

#include <map>
#include <set>
#include <vector>

#include "TF1.h"
#include "TGraphAsymmErrors.h"
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
    bool hasMCTruthSample() const { return ptBins_.front()->hasMCTruthSample(); }
    SampleLabel mcTruthSampleLabel() const { return ptBins_.front()->mcTruthSample()->label(); }

    double meanPt(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    double meanPtStatUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    double fittedValue(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx, unsigned int ptSoftIdx) const;
    double fittedStatUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx, unsigned int ptSoftIdx) const;
    TGraphAsymmErrors* kSoftSlope(const SampleLabel &label, FitResult::Type type) const;
    TF1* kSoftFit(const SampleLabel &label, FitResult::Type type, const TString &name) const;
    double extrapolatedValue(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    double extrapolatedStatUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    double extrapolatedSystUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    double correctedResolution(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx, double scalePLI = 1.) const;
    double correctedResolutionStatUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    double pli(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    double mcTruthResolution(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    double relativeMCClosure(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const;
    TGraphAsymmErrors* biasCorrectedResolution(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type) const;

    TGraphAsymmErrors* extrapolatedResolution(const SampleLabel &label, FitResult::Type type) const;
    TGraphAsymmErrors* correctedResolution(const SampleLabel &label, FitResult::Type type) const;
    TF1* mcTruthResoFunc(const TString &name) const { return mcTruthReso_->func(name); }
    TF1* pliFunc(const TString &name) const { return pli_->func(name); }

    TF1* scaledMCTruthResoFunc(const TString &name) const { return scaledMCTruthReso_->func(name); }
    TGraphAsymmErrors* scaledMCTruthUncertaintyBand() const;
    TGraphAsymmErrors* scaledMCTruthRatioBand() const;

    bool hasKValue(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type) const;
    double kValue(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type) const { return kVal_; }
    double kStat(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type) const { return kValStat_; }
    double kSystDown(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type) const { return kValSystDown_; }
    double kSystUp(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type) const { return kValSystUp_; }
    double kTotalDown(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type) const { return kValTotalDown_; }
    double kTotalUp(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type) const { return kValTotalUp_; }
    TF1* kValueLine(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type, const TString &name, double xMin, double xMax) const;
    TGraphAsymmErrors* kValueStatBand(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type, double xMin, double xMax) const;
    TGraphAsymmErrors* kValueSystBand(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type, double xMin, double xMax) const;

    TGraphAsymmErrors* ratioGraph(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type) const;

    bool addDataSample(const TString &label, const TString &baseFileName);
    bool addMCSample(const TString &label, const TString &baseFileName);
    bool addFitResult(FitResult::Type type);

    bool addExtrapolationUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color);
    bool addPLIUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color);
    bool addMCClosureUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color);
    bool addUncertaintyFromVariedSample(const TString &uncertaintyLabel, double fraction, const SampleLabel &nominalSample, FitResult::Type type, const TString &variedSampleDown, const TString &variedSampleUp, int color);

    bool compareSamples(const SampleLabel &label1, const SampleLabel &label2);
    void fitKValue(FitResult::Type type);

    void setMCTruthResolution(ResolutionFunction* mcTruthReso);
    void fitPLI(const TString &label, const TString &baseFileName, ResolutionFunction::Type type);
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
    ResolutionFunction* scaledMCTruthReso_;
    ResolutionFunction* pli_;

    // Ratios of resolution of different samples
/*     std::map<SampleLabelPair,double> kVal_;  */
/*     std::map<SampleLabelPair,double> kValStat_;  */
/*     std::map<SampleLabelPair,double> kValSyst_;  */

    FitResult::Type kValType_;
    double kVal_; 
    double kValStat_; 
    double kValSystDown_; 
    double kValSystUp_; 
    double kValTotalDown_; 
    double kValTotalUp_; 

    SystematicUncertainty* findSystematicUncertainty(const SampleLabel &label, FitResult::Type type);
    TF1* fitKSoftSlope(const TString &name, const SampleLabel &label, FitResult::Type type) const;
    ResolutionFunction* fitResolution(const TGraphAsymmErrors* g, ResolutionFunction::Type type) const;
  };

  typedef std::vector<EtaBin*> EtaBins;
  typedef std::vector<EtaBin*>::const_iterator EtaBinIt;
}
#endif

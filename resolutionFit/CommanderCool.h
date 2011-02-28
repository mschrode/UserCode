// $Id: CommanderCool.h,v 1.5 2011/02/26 17:55:50 mschrode Exp $

#ifndef COMMANDER_COOL_H
#define COMMANDER_COOL_H

#include <vector>

#include "TString.h"

#include "EtaBin.h"
#include "FitResult.h"
#include "Parameters.h"
#include "PlotMaker.h"
#include "ResolutionFunction.h"


namespace resolutionFit {

  // Main steering class controling the workflow
  //
  // Make Singleton
  class CommanderCool {
  public:
    CommanderCool(const Parameters* par);
    ~CommanderCool();

    void setMCTruthResolution(const TString &fileName, ResolutionFunction::Type type);
    void fitPLI(const TString &label, const TString &baseFileName, ResolutionFunction::Type type);
    void setPLI(const TString &fileName, ResolutionFunction::Type type);

    void addDataSample(const TString &label, const TString &baseFileName);
    void addMCSample(const TString &label, const TString &baseFileName);
    void addFitResult(FitResult::Type type);
    void addExtrapolationUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color);
    void addPLIUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color);
    void addMCClosureUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color);
    void addUncertaintyFromVariedSample(const TString &uncertaintyLabel, double fraction, const SampleLabel &nominalSample, FitResult::Type type, const TString &variedSampleDown, const TString &variedSampleUp, int color);
    void compareSamples(const SampleLabel &label1, const SampleLabel &label2);

    void printSetup() const;
    void printResult() const;
    void makeAllPlots() const { plotMaker_->makeAllPlots(); }


  private:
    const Parameters* par_;

    PlotMaker* plotMaker_;
    EtaBins etaBins_;

    bool isConsistentInputName(const TString &name) const;
  };
}
#endif

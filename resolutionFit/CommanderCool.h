// $Id: CommanderCool.h,v 1.2 2011/02/18 18:42:22 mschrode Exp $

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
    void setPLI(const TString &fileName, ResolutionFunction::Type type);

    void addDataSample(const TString &label, const TString &baseFileName);
    void addMCSample(const TString &label, const TString &baseFileName);
    void addFitResult(FitResult::Type type);
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

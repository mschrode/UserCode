// $Id: $

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

    void setMCTruthResolution(ResolutionFunction::Type type);
    void setPLI(ResolutionFunction::Type type);

    void addDataSample(const TString &label, const TString &baseFileName);
    void addMCSample(const TString &label, const TString &baseFileName);
    void addFitResult(FitResult::Type type);

    void printSetup() const;
    void printResult() const;
    void makeAllPlots() const { plotMaker_->makeAllPlots(); }


  private:
    const Parameters* par_;

    PlotMaker* plotMaker_;
    EtaBins etaBins_;
  };
}
#endif

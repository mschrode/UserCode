// $Id: PlotMaker.h,v 1.1 2011/02/15 18:22:25 mschrode Exp $

#ifndef PLOT_MAKER_H
#define PLOT_MAKER_H

#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"

#include "EtaBin.h"
#include "OutputManager.h"
#include "Parameters.h"
#include "PtBin.h"
#include "Sample.h"


namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  class PlotMaker {
  public:
    PlotMaker(const Parameters *par, const EtaBins &etaBins);
    ~PlotMaker();

    void makeAllPlots() const { 
      plotExtrapolation();
      plotPtGenSpectra();
      plotResolution();
    }


    void plotExtrapolation() const;
    void plotPtGenSpectra() const;
    void plotResolution() const;


  private:
    class LabelMaker {
    public:
      LabelMaker(const Parameters* par);

      double start() const { return 0.75; }

      TPaveText* ptSoftBin(SampleLabel label, unsigned int etaBinIdx, unsigned int ptBinIdx, unsigned int ptSoftBinIdx) const;
      TPaveText* ptBin(SampleLabel label, unsigned int etaBinIdx, unsigned int ptBinIdx) const;
      TPaveText* etaBin(SampleLabel label, unsigned int etaBinIdx, unsigned int nExtraEntries = 0) const;
      TString etaRange(unsigned int etaBin) const;
      TString ptRange(unsigned int etaBin, unsigned int ptBin) const;
      TString ptSoftRange(unsigned int ptSoftBinIdx) const;
      TString jets() const;
      TString pt() const;
      TString ptSoft() const;
      TString label(FitResult::Type type) const;
      

    private:
      const Parameters* par_;
    };

    const Parameters* par_;
    const EtaBins etaBins_;

    const double xMinPt_;
    const double xMaxPt_;
    const double yMinExtraRes_;
    const double yMaxExtraRes_;
    const double yMinResRatio_;
    const double yMaxResRatio_;

    OutputManager* out_;
    LabelMaker* labelMk_;

    TString histFileName(const TString &id, const EtaBin* etaBin, SampleLabel sampleLabel, FitResult::Type type) const;
    TString histFileName(const TString &id, const PtBin* ptBin, const Sample* sample, FitResult::Type type) const;
    TString histFileName(const TString &id, const PtBin* ptBin, const Sample* sample, unsigned int ptSoftBinIdx) const;

    int markerStyleExtrapolatedResolution(Sample::Type type) const;
    int markerStyleCorrectedResolution(Sample::Type type) const;
  };
}
#endif

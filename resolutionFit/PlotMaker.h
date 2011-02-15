// $Id: $

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
      
/*       TPaveText* binLabel(unsigned int etaBin, unsigned int ptBin) const; */
/*       TPaveText* binLabel(const Bin* bin) const { */
/* 	return binLabel(bin->etaBin(),bin->ptBin()); */
/*       } */
/*       TPaveText* binLabel(unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBinIdx) const; */
/*       TPaveText* binLabel(const Bin* bin, unsigned int ptSoftBinIdx) const { */
/* 	return binLabel(bin->etaBin(),bin->ptBin(),ptSoftBinIdx); */
/*       } */
      
/*       TString etaRange(unsigned int etaBin) const; */
/*       TString ptRange(unsigned int etaBin, unsigned int ptBin) const; */
/*       TString ptSoftRange(unsigned int ptSoftBinIdx) const; */
/*       TString jetAlgo() const; */
/*       TString pt() const; */
/*       TString ptSoft() const; */
      
      
    private:
      const Parameters* par_;
    };

    const Parameters* par_;
    const EtaBins etaBins_;

    const double xMinPt_;
    const double xMaxPt_;
    const double yMinExtraRes_;
    const double yMaxExtraRes_;

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

// $Id: Measurement.h,v 1.10 2011/08/19 08:33:56 mschrode Exp $

#ifndef MEASUREMENT_H
#define MEASUREMENT_H


#include <map>
#include <vector>

#include "TString.h"


class TH1;

namespace resolutionFit {

  //! \brief Reads and stores control histograms and fitted parameter values
  //! from Kalibri resolution fits
  //!
  //! This is the interface to the resolution fit results and control
  //! plots produced by the Kalibri tool described at
  //! https://twiki.cern.ch/twiki/bin/view/CMS/HamburgWikiAnalysisCalibration
  //! and specifically the output of http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Bromo/Calibration/CalibCore/ControlPlotsResolution.h?view=log
  //! The histograms are expected to be stored in a ROOT file as produced
  //! by this script,
  //! http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/mschrode/sampleTools/combineBins.C?view=log
  //! which collects the Kalibri output from the different eta, pt, and
  //! ptSoft bins.
  // -------------------------------------------------------------------
  class Measurement {
  public:
    Measurement(const TString &fileName, unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBin, double ptMin, double ptMax, double ptSoft, unsigned int verbosity = 0);

    ~Measurement();

    double ptMin() const { return ptMin_; }
    double ptMax() const { return ptMax_; }
    double ptSoft() const { return ptSoft_; }

    unsigned int nFittedValues() const { return values_.size(); }
    double fittedValue(unsigned int i) const { return values_.at(i); }
    double fittedUncert(unsigned int i) const { return statUncert_.at(i); }
    double pdfPtTrue(double pt) const;
    double meanPtGen() const { return meanPtGen_; }
    double meanPtGenUncert() const { return meanPtGenUncert_; }
    double meanPtAve() const { return meanPtAve_; }
    double meanPtAveUncert() const { return meanPtAveUncert_; }
    double meanPdfPtTrue() const { return meanPdfPtTrue_; }
    double meanPdfPtTrueUncert() const { return meanPtAveUncert_; }

    TH1* histPdfPtAsym() const { return getClone("hFitPtAsym"); }
    TH1* histPdfPtTrue() const { return getClone("hTruthPDF"); }
    TH1* histPtAsym() const { return getClone("hPtAbsAsym"); }
    TH1* histPtGenAsym() const { return getClone("hPtGenAsym"); }
    TH1* histPtGen() const { return getClone("hPtGen"); }
    TH1* histPtAve() const { return getClone("hPtAveAbs"); }
    TH1* histPt1() const { return getClone("hPtJet1"); }
    TH1* histPt2() const { return getClone("hPtJet2"); }
    TH1* histPt3() const { return getClone("hPtJet3"); }
    TH1* histNumVtx() const { return getClone("hNumVtx"); }
    TH1* histMCWeight() const { return getClone("hWeights"); }
    TH1* histMCNumPU() const { return getClone("hNumPU"); }
    TH1* histDeltaPt() const { return getClone("hDeltaPtJet12"); }




  private:
    static unsigned int HIST_COUNT;

    typedef std::map<TString,TH1*> HistMap;
    typedef std::map<TString,TH1*>::iterator HistMapIt;

    const double ptMin_;
    const double ptMax_;
    const double ptSoft_;
    const unsigned int verbosity_;

    TString hnSuffix_;		//! Histogram name suffix
    TString hnPrefix_;		//! Histogram name prefix
    std::vector<double> values_;     //! Fitted parameter values
    std::vector<double> statUncert_; //! Statistical uncertainty of the fitted parameter values
    HistMap hists_;		//! The stored histograms: [Name][Pointer to histogram]

    double meanPtGen_;		//! Mean value of ptGen distribution
    double meanPtGenUncert_;	//! Uncertainty on mean value of ptGen distribution
    double meanPtAve_;         	//! Mean value of ptAve distribution
    double meanPtAveUncert_;	//! Uncertainty on mean value of ptAve distribution
    double meanPdfPtTrue_;	//! Mean value of ptTrue pdf
    double meanPdfPtTrueUncert_;	//! Uncertainty on mean value of ptTrue pdf

    TH1* hSpectrum_;		//! The assumed underlying dijet spectrum, i.e. without modification for migration effects

    TH1* getClone(const TString &name) const;
    TString histNamePath(unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBin) const;
    TString histNameSuffix(unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBin) const;
    void init(const TString &fileName, unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBin);
    bool parse(const TString &fileName);
    void setMean(const TString &name, double &mean, double &uncert);
  };

  typedef std::vector<Measurement*> Meas;
  typedef std::vector<Measurement*>::const_iterator MeasIt;
}
#endif

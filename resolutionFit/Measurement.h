// $Id: Measurement.h,v 1.4 2011/02/28 10:53:15 mschrode Exp $

#ifndef MEASUREMENT_H
#define MEASUREMENT_H


#include <map>
#include <vector>

#include "TString.h"


class TH1;

namespace resolutionFit {
  class Measurement {
  public:
    Measurement(const TString &fileName, const TString &histNameSuffix, double ptSoft, unsigned int verbosity = 0);
    Measurement(const TString &fileName, const TString &histNameSuffix, double ptSoft, bool hasFittedParameters, unsigned int verbosity = 0);

    ~Measurement();

    double ptSoft() const { return ptSoft_; }

    unsigned int nFittedValues() const { return values_.size(); }
    double fittedValue(unsigned int i) const { return values_.at(i); }
    double fittedUncert(unsigned int i) const { return statUncert_.at(i); }
    double meanPtGen() const { return meanPtGen_; }
    double meanPtGenUncert() const { return meanPtGenUncert_; }
    double meanPtAve() const { return meanPtAve_; }
    double meanPtAveUncert() const { return meanPtAveUncert_; }
    double meanPdfPtTrue() const { return meanPdfPtTrue_; }
    double meanPdfPtTrueUncert() const { return meanPdfPtTrueUncert_; }

    TH1* histPdfPtAsym() const { return getClone("hFitPtAsym"); }
    TH1* histPdfPtTrue() const { return getClone("hTruthPDF"); }
    TH1* histPtAsym() const { return getClone("hPtAsym"); }
    TH1* histPtGenAsym() const { return getClone("hPtGenAsym"); }
    TH1* histPtGen() const { return getClone("hPtGen"); }
    TH1* histPtAve() const { return getClone("hPtAveAbs"); }
    TH1* histPt1() const { return getClone("hPtJet1"); }
    TH1* histPt2() const { return getClone("hPtJet2"); }
    TH1* histPt3() const { return getClone("hPtJet3"); }


  private:
    static unsigned int HIST_COUNT;

    typedef std::map<TString,TH1*> HistMap;
    typedef std::map<TString,TH1*>::iterator HistMapIt;

    const TString histNameSuffix_;
    const double ptSoft_;
    const unsigned int verbosity_;

    std::vector<double> values_;     //! Fitted parameter values
    std::vector<double> statUncert_; //! Statistical uncertainty of the fitted parameter values
    HistMap hists_;

    double meanPtGen_;		//! Mean value of ptGen distribution
    double meanPtGenUncert_;	//! Uncertainty on mean value of ptGen distribution
    double meanPtAve_;         	//! Mean value of ptAve distribution
    double meanPtAveUncert_;	//! Uncertainty on mean value of ptAve distribution
    double meanPdfPtTrue_;	//! Mean value of ptTrue pdf
    double meanPdfPtTrueUncert_;	//! Uncertainty on mean value of ptTrue pdf

    TH1* getClone(const TString &name) const;
    void init(const TString &fileName, bool hasFittedParameters);
    bool parse(const TString &fileName, bool hasFittedParameters);
    void setMean(const TString &name, double &mean, double &uncert);
  };

  typedef std::vector<Measurement*> Meas;
  typedef std::vector<Measurement*>::const_iterator MeasIt;
}
#endif

// $Id: $

#ifndef KALIBRI_FILE_PARSER_H
#define KALIBRI_FILE_PARSER_H

#include <cassert>
#include <vector>
#include <map>

#include "TString.h"

class TH1F;

namespace resolutionFit {

  //! \brief Parse Kalibri fit result files
  //!
  //! Parses the ROOT-files produced by the Kalibri
  //! \p ControlPlotsJetSmearing class. So far, the
  //! following information can be obtained from the
  //! parser
  //!  - fitted parameter values and uncertainties
  //!    (from the histograms written by
  //!     \p ControlPlotsJetSmearing::plotParameters() )
  //!  - different estimators of the mean pt
  //!    (from the histograms written by
  //!     \p ControlPlotsJetSmearing::plotResponse() )
  //! Furthermore, the following histograms (written by
  //! \p ControlPlotsJetSmearing::plotResponse() ) can be
  //! obtained from the parser
  //!  - ptGen spectrum
  //!  - ptDijet spectrum
  //!  - ptTrue pdf assumed by the fit
  //!  - response distribution pt / ptGen
  //!  - fitted resolution
  //!
  //! \author Matthias Schroeder
  //! \date 2009/03/05
  //! $Id: $
  // --------------------------------------------
  class KalibriFileParser {
  public:
    //! Constructor
    KalibriFileParser(const TString &fileName, int verbose = 1);
    //! Destructor
    ~KalibriFileParser();

    //! Returns number of fitted parameters
    int nValues() const { return static_cast<int>(values_.size()); }
    //! Returns value of the i-th fitted parameter
    double value(int i = 0) const { assert( i>=0 && i < nValues() ); return values_[i]; }
    //! Returns statistical uncertainty of the i-th parameter
    double statUncert(int i = 0) const { assert( i>=0 && i < nValues() ); return statUncert_[i]; }
    //! Returns mean value of the ptGen distribution
    double meanPtGen() const { return meanPtGen_; }
    //! Returns mean value of the ptDijet distribution
    double meanPtDijet() const { return meanPtDijet_; }
    //! Returns mean value of the ptTrue pdf assumed in the fit
    double meanPdfPtTrue() const { return meanPdfPtTrue_; }
    //! Returns the specified histogram filled by \p ControlPlotsJetSmearing
    TH1F *hist(const TString &name, const TString &newName) const;

  private:
    typedef std::map<TString,TH1F*>::const_iterator HistIt;

    const int verbose_;		//! Verbosity level

    std::vector<double> values_; //! Fitted parameter values
    std::vector<double> statUncert_; //! Statistical uncertainty of the fitted parameter values
    std::map<TString,TH1F*> hists_; //! The histograms read from ROOT file

    double meanPtGen_;		//! Mean value of ptGen distribution
    double meanPtDijet_;	//! Mean value of ptDijet distribution
    double meanPdfPtTrue_;	//! Mean value of ptTrue pdf

    //! Parse the ROOT file with name \p fileName
    int parse(const TString &fileName);
    //! Set the mean values (attributes) from the corresponding distributions
    void setMeanPt();
  };
}
#endif


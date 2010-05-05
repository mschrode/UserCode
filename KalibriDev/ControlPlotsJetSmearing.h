//
// $Id: ControlPlotsJetSmearing.h,v 1.10 2010/04/18 14:40:38 mschrode Exp $
//
#ifndef JS_CONTROLPLOTS_JETSMEARING_H
#define JS_CONTROLPLOTS_JETSMEARING_H

#include <string>
#include <vector>

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"

class Jet;
class TCanvas;
class TH1;
class TLegend;
class TObject;
class TPaveText;
class TPostScript;
class TRandom3;


//!  \brief Generates validation plots for jet-smearing method
//!  \author Matthias Schroeder
//!  \date Thu May  7 11:30:28 CEST 2009 
//!  $Id: ControlPlotsJetSmearing.h,v 1.10 2010/04/18 14:40:38 mschrode Exp $
// --------------------------------------------------
class ControlPlotsJetSmearing {
 public:
  ControlPlotsJetSmearing(const std::string& configfile,const std::vector<Event*> * data, TParameters * param, const std::string &outDir = "./controlPlots");
  ~ControlPlotsJetSmearing();

  void makePlots() const;
  void setBinningResp(int nbins, double min, double max) { respNBins_ = nbins; respMin_ = min; respMax_ = max;}


 private:
  static double spectrum(double *x, double *par);

  typedef std::vector<Event*>::const_iterator DataIt;

  const std::vector<Event*> * data_;   //!< The data which is plotted
  const ConfigFile          * config_; //!< The configuration file
  mutable TParameters       * param_;  //!< The parametrization
  
  int          respNBins_;             //!< Number of bins in response control plots \p plotResponse()
  double       respMin_;               //!< Minimum of response control plots \p plotResponse()
  double       respMax_;               //!< Maximum of response control plots \p plotResponse()
  std::string  dir_;                   //!< Directory in which the control plots are written
  TRandom3 *rand_;
  std::string parClass_;
  std::vector<double> startParJet_;
  std::vector<double> scale_;
  std::vector<double> truthPar_;
  std::string ptBinningVar_;
  std::vector<double> ptBinEdges_;
  std::vector<double> ptBinCenters_;

  void plotDijets() const;
  void plotResponse() const;
  void plotParameters() const;
  //! Plots the negative log-likelihood for different parameter values
  void plotParameterScan() const;
  //! Plots the distributions of the probability density of
  //! each event before and after the fit
  void plotLogP() const;
  void plotMeanResponseAndResolution() const;
  void plot3rdJet() const;

  double gaussianWidth(double pt) const;
  double gaussianWidthError(double pt) const;
  double gaussianWidthTruth(double pt) const;

  int nPtBins() const { return static_cast<int>(ptBinEdges_.size()-1); }
  double ptBinsMin() const { return ptBinEdges_.front(); }
  double ptBinsMax() const { return ptBinEdges_.back(); }
  int findPtBin(const Jet *jet) const;
  int findPtBin(double pt) const;
  int findBin(double x, const std::vector<double> &binEdges) const;
  bool equidistLogBins(std::vector<double>& bins, int nBins, double first, double last) const;

  TLegend *createLegend(int nEntries, double width = 1., double lineHgt = -1., double yOffset = 0.) const;
  TPaveText *createPaveText(int nEntries, double width = 1., double lineHgt = -1.) const;
  void drawPSPage(TPostScript * ps, TCanvas * can, TObject * obj, const std::string &option = "", bool log = false) const;
  void drawPSPage(TPostScript * ps, TCanvas * can, std::vector<TObject*> objs, const std::string &option = "", bool log = false) const;
  void findYRange(const TH1 * h, double& min, double& max) const;
  double lineHeight() const { return 0.06; }
  void normHist(TH1 * h, std::string option = "") const;
  void normHist(TH1 *h, double min, double max, std::string option = "") const;
  void setGStyle() const;
  void setYRange(TH1 * h, double c1 = 0.9, double c2 = 1.1, double minLimit = 0.) const;
  template <class T> std::string toString(const T& t) const;
};
#endif

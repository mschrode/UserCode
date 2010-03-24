//
// $Id: ControlPlotsJetSmearing.h,v 1.6 2010/02/25 15:28:18 mschrode Exp $
//
#ifndef JS_CONTROLPLOTS_JETSMEARING_H
#define JS_CONTROLPLOTS_JETSMEARING_H

#include <string>
#include <vector>

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"

class TH1F;
class TH2F;
class TObject;
class TPostScript;
class TCanvas;
class TRandom3;


//!  \brief Generates validation plots for jet-smearing method
//!  \author Matthias Schroeder
//!  \date Thu May  7 11:30:28 CEST 2009 
//!  $Id: ControlPlotsJetSmearing.h,v 1.6 2010/02/25 15:28:18 mschrode Exp $
// --------------------------------------------------
class ControlPlotsJetSmearing {
 public:
  ControlPlotsJetSmearing(const std::string& configfile,const std::vector<Event*> * data, TParameters * param, const std::string &outDir = "./controlPlots");
  ~ControlPlotsJetSmearing();

  void makePlots() const;
  void setBinningResp(int nbins, double min, double max) { respNBins_ = nbins; respMin_ = min; respMax_ = max;}


 private:
  typedef std::vector<Event*>::const_iterator DataIt;

  const std::vector<Event*> * data_;   //!< The data which is plotted
  const ConfigFile          * config_; //!< The configuration file
  mutable TParameters       * param_;  //!< The parametrization
  
  int          respNBins_;             //!< Number of bins in response control plots \p plotResponse()
  double       respMin_;               //!< Minimum of response control plots \p plotResponse()
  double       respMax_;               //!< Maximum of response control plots \p plotResponse()
  std::string  dir_;                   //!< Directory in which the control plots are written
  TRandom3 *rand_;

  void plotDijets() const;
  void plotResponse() const;
  void plotParameters() const;
  //! Plots the negative log-likelihood for different parameter values
  void plotParameterScan() const;
  //! Plots the distributions of the probability density of
  //! each event before and after the fit
  void plotLogP() const;
  void plotMeanResponseAndResolution() const;

  void drawPSPage(TPostScript * ps, TCanvas * can, TObject * obj, std::string option = "", bool log = false) const;
  void drawPSPage(TPostScript * ps, TCanvas * can, std::vector<TObject*> objs, std::string option = "", bool log = false) const;
  void findYRange(const TH1F * h, double& min, double& max) const;
  void normHist(TH1F * h, std::string option = "") const;
  void normHist(TH1F *h, double min, double max, std::string option = "") const;
  void setGStyle() const;
  void setYRange(TH1F * h, double c1 = 0.9, double c2 = 1.1, double minLimit = 0.) const;
  template <class T> std::string toString(const T& t) const;
};
#endif

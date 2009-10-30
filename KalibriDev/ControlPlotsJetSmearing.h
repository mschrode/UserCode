// $Id: ControlPlotsJetSmearing.h,v 1.5 2009/08/07 12:19:24 mschrode Exp $

#ifndef JS_CONTROLPLOTS_JETSMEARING_H
#define JS_CONTROLPLOTS_JETSMEARING_H

#include <string>
#include <vector>

#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TObject.h"
#include "TPostScript.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"

using std::vector;



//!  \brief Generates validation plots for jet-smearing method
//!  \author Matthias Schroeder
//!  \date Thu May  7 11:30:28 CEST 2009 
//!  $Id: ControlPlotsJetSmearing.h,v 1.5 2009/08/07 12:19:24 mschrode Exp $
// --------------------------------------------------
class ControlPlotsJetSmearing {
 public:
  ControlPlotsJetSmearing(const std::string& configfile,const std::vector<TData*> * data, TParameters * param);
  ~ControlPlotsJetSmearing() {};
  
  void plotDijets() const;
  void plotMeanResponseAndResolution() const;
  void plotResponse() const;
  void plotParameterScan(const std::vector<unsigned int>& pars) const;
  void setBinningResp(int nbins, double min, double max) { respNBins_ = nbins; respMin_ = min; respMax_ = max;}


 private:
  typedef std::vector<TData*>::const_iterator DataIt;

  const std::vector<TData*> * data_;   //!< The data which is plotted
  const ConfigFile          * config_; //!< The configuration file
  mutable TParameters       * param_;  //!< The parametrization
  
  int          respNBins_;             //!< Number of bins in response control plots \p plotResponse()
  double       respMin_;               //!< Minimum of response control plots \p plotResponse()
  double       respMax_;               //!< Maximum of response control plots \p plotResponse()
  std::string  dir_;                   //!< Directory in which the control plots are written

  void drawPSPage(TPostScript * ps, TCanvas * can, TObject * obj, std::string option = "", bool log = false) const;
  void drawPSPage(TPostScript * ps, TCanvas * can, std::vector<TObject*> objs, std::string option = "", bool log = false) const;
  void findYRange(const TH1F * h, double& min, double& max) const;
  void fitSlices(const TH2F * h2, std::vector<TH1F*>& hFit) const;
  void normHist(TH1F * h, std::string option = "") const
    { if( h->Integral(option.c_str()) ) h->Scale(1./h->Integral(option.c_str())); }
  void setGStyle() const;
  void setYRange(TH1F * h, double c1 = 0.9, double c2 = 1.1, double minLimit = 0.) const;
  template <class T> std::string toString(const T& t) const;
};
#endif

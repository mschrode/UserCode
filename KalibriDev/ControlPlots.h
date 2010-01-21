#ifndef PLOTS_H
#define PLOTS_H

#include <string>
#include <vector>

#include "ControlPlotsConfig.h"
#include "ControlPlotsFunction.h"

class ConfigFile;
class Event;


//!  \brief Create control plots
//!
//!  Creates control plots via the \p makePlots() method from
//!  different \p Event types. The output is in .eps and .root format.
//!  The attributes of the control plots are  specified via 
//!  the configuration file.
//!
//!  \author Christian Autermann
//!  \date Fri Jan 18 13:55:15 2008 UTC
//!  $Id: ControlPlots.h,v 1.34 2010/01/04 17:04:51 mschrode Exp $
// -------------------------------------------------------------
class ControlPlots {
 public:
  typedef std::vector<Event*>::const_iterator DataIt;

  ControlPlots(const ConfigFile *configFile, const std::vector<Event*> *data);
  ~ControlPlots() {};

  void makePlots() const;

 private:
  const ConfigFile *config_;                    //!< The configuration file
  const std::vector<Event*> *data_;             //!< The plotted data

  void createJetTruthEventPlots() const;
  ControlPlotsFunction::Function findJetTruthEventFunction(const std::string& varName, ControlPlotsConfig::CorrectionType type = ControlPlotsConfig::Uncorrected) const;
  void setGStyle() const;
};

#endif

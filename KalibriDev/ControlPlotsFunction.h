// $Id: ControlPlotsFunction.h,v 1.1 2010/01/04 17:04:51 mschrode Exp $

#ifndef CONTROLPLOTS_FUNCTION_H
#define CONTROLPLOTS_FUNCTION_H

#include "ControlPlotsConfig.h"

class Event;


//!  \brief Functions for the profile control plots
//!
//!  Stores the functions for the profile control plots. This class
//!  is used by \p ControlPlotsProfile and returns for a given \p Event
//!  the x value, the y value and the binning value for the profile
//!  plots. There are different y values depending on the applied
//!  correction.
//!
//!  For example, for a "Response vs GenJetPt" plot in bins of "Eta",
//!  this class stores the functions to get Response, GenJetPt and
//!  Eta from an \p Event.
//!
//!  \sa \p ControlPlotsConfig, \p ControlPlotsProfile  
//!
//!  \author Matthias Schroeder
//!  \date 2009/12/18
//!  $Id: ControlPlotsFunction.h,v 1.1 2010/01/04 17:04:51 mschrode Exp $
// ----------------------------------------------------------------   
class ControlPlotsFunction {
 public:
  //! The function's signature
  typedef double (ControlPlotsFunction::*Function)(const Event *evt) const;

  //! Constructor
  ControlPlotsFunction()
    : binFunc_(0), xFunc_(0) {};

  //! Check if all functions are initialised
  bool isInit() const { return binFunc_ && xFunc_ && yFuncs_.size(); }

  //! Return iterator to first \p CorrectionType
  ControlPlotsConfig::CorrectionTypeIt correctionTypeBegin() const { return types_.begin(); }
  //! Return iterator to last \p CorrectionType
  ControlPlotsConfig::CorrectionTypeIt correctionTypeEnd() const { return types_.end(); }

  //! Interface to the profile: return the value of the binning quantity from \p evt
  double binValue(const Event * evt) const { return (this->*binFunc_)(evt); }
  //! Interface to the profile: return the value of the x quantity from \p evt
  double xValue(const Event * evt) const { return (this->*xFunc_)(evt); }
  //! Interface to the profile: return the value of the y quantity from \p evt for the correction type \p type
  double yValue(const Event * evt, ControlPlotsConfig::CorrectionType type) const;

  //! Set the function returning the binning value from an event
  void setBinFunction(Function func) { binFunc_ = func; }
  //! Set the function returning the x value from an event
  void setXFunction(Function func) { xFunc_ = func; }
  //! Set the functions returning the y value for different corrections from an event
  void addYFunction(ControlPlotsConfig::CorrectionType type, Function func);

  double jetTruthEventJetEta(const Event *evt) const;
  double jetTruthEventTruthPt(const Event *evt) const;
  double jetTruthEventResponse(const Event * evt) const;
  double jetTruthEventResponseKalibriCorrected(const Event * evt) const;
  double jetTruthEventResponseL2L3Corrected(const Event * evt) const;
    
 private:
  //! The different correction types of the y quantity
  std::vector<ControlPlotsConfig::CorrectionType> types_;
  //! The binning function
  Function binFunc_;
  //! The x value function
  Function xFunc_;
  //! The y value functions for the different correction types
  std::map<ControlPlotsConfig::CorrectionType,Function> yFuncs_;    
};
#endif

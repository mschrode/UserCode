// $Id: CoreScaleFactors.h,v 1.1 2013/05/08 13:07:28 mschrode Exp $

#ifndef RESOLUTION_TAILS_CORE_SCALE_FACTORS
#define RESOLUTION_TAILS_CORE_SCALE_FACTORS

#include "Uncertainty.h"

// Container class to store core-scale factors with uncertainties
// for various analysis versions

namespace resolutionTails {
  class CoreScaleFactors {
  public:
    enum Version { Run2012ABCReReco53X };
    
    CoreScaleFactors(CoreScaleFactors::Version version)
      : version_(version) {};
    
    std::vector<double> operator()(Uncertainty::SystematicVariation var) const;
    

  private:
    Version version_;
  };


  // Return vector with core scale factors
  // One factor per eta bin
  std::vector<double> CoreScaleFactors::operator()(resolutionTails::Uncertainty::SystematicVariation var) const {
    std::vector<double> factors;

    if( version_ == Run2012ABCReReco53X ) {
      // Core scale factors from runs 190456-202305 (Run2012ABC ReReco 53X)
      // Version 2012/11/28
      if( var == Uncertainty::CoreDn ) {	// core down 1 sigma
	factors.push_back(0.);
	factors.push_back(0.);
	factors.push_back(0.);
	factors.push_back(0.);
	factors.push_back(0.);
      } else if( var == Uncertainty::CoreUp ) { // core up 1 sigma 
	factors.push_back(0.);
	factors.push_back(0.);
	factors.push_back(0.);
	factors.push_back(0.);
	factors.push_back(0.);
      } else {              	// nominal core scale
	factors.push_back(0.094);
	factors.push_back(0.106);
	factors.push_back(0.155);
	factors.push_back(0.272);
	factors.push_back(0.080);
      }
    }
    
    return factors;
  }
}
#endif

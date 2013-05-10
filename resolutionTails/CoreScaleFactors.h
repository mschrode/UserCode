// $Id: CoreScaleFactors.h,v 1.2 2013/05/10 13:22:19 mschrode Exp $

#ifndef RESOLUTION_TAILS_CORE_SCALE_FACTORS
#define RESOLUTION_TAILS_CORE_SCALE_FACTORS

#include "Uncertainty.h"

// Container class to store core-scale factors with uncertainties
// for various analysis versions

namespace resolutionTails {
  class CoreScaleFactors {
  public:
    enum Version { Run2011A_42X_S3, Run2011AB_42X_S6, Run2012ABCReReco_53X_S10 };
    
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

    if( version_ == Run2011A_42X_S3 ) {
      // Core scale factors from runs 163337-167151 (1/fb)
      // Version 2011/07/20: closure uncertainty removed
      if( var == Uncertainty::CoreDn ) {	// core down 1 sigma
	factors.push_back(0.000);
	factors.push_back(0.001);
	factors.push_back(0.032);
	factors.push_back(0.043);
	factors.push_back(0.090);
      } else if( var == Uncertainty::CoreUp ) {	// core up 1 sigma
	factors.push_back(0.116);
	factors.push_back(0.115);
	factors.push_back(0.162);
	factors.push_back(0.228);
	factors.push_back(0.489);
      } else {	// nominal core scale
	factors.push_back(0.052);
	factors.push_back(0.057);
	factors.push_back(0.096);
	factors.push_back(0.134);
	factors.push_back(0.288);
      }
    }
    
    // Core scale factors from runs 163337-167151
    // Version 2012/06/09 (thesis):
    // - corrected MC statistical uncertainties
    // - corrected extrapolation uncertainty determination
    // - corrected uncertainty propagation to ratio
    else if( version_ == Run2011AB_42X_S6 ) {
      if( var == Uncertainty::CoreDn ) {	// core down 1 sigma
	factors.push_back(0.0000);
	factors.push_back(0.0000);
	factors.push_back(0.0192);
	factors.push_back(0.0283);
	factors.push_back(0.0737);
      } else if( var == Uncertainty::CoreUp ) {	// core up 1 sigma
	factors.push_back(0.1232);
	factors.push_back(0.1254);
	factors.push_back(0.1659);
	factors.push_back(0.2567);
	factors.push_back(0.5174);
      } else {	// nominal core scale
	factors.push_back(0.0541);
	factors.push_back(0.0599);
	factors.push_back(0.0920);
	factors.push_back(0.1409);
	factors.push_back(0.2942);
      }
    }

    // Core scale factors from runs 190456-202305 (Run2012ABC ReReco 53X)
    // Version 2012/11/28
    if( version_ == Run2012ABCReReco_53X_S10 ) {
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

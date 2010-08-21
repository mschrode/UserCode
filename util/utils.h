// $Id: utils.h,v 1.3 2010/07/27 17:08:39 mschrode Exp $

#ifndef UTILS_H
#define UTILS_H

#include <cmath>

#include <string>
#include <sstream>

#include "TString.h"


//!  Encapsulates useful classes and methods
//!
//!  \author   Matthias Schroeder (www.desy.de/~matsch)
//!  \date     2010/03/09
//!  $Id: utils.h,v 1.3 2010/07/27 17:08:39 mschrode Exp $
// -------------------------------------------------------------------------------------
namespace util {

  //! Returns error of A / B for two uncorrelated variables
  //! A, B with error AE, BE
  // -------------------------------------------------------------------------------------
  static double ratioError(double A, double AE, double B, double BE) {
    return sqrt( AE*AE/B/B + BE*BE*A*A/B/B/B/B );
  }


  // -------------------------------------------------------------------------------------
  static double round(double d, int decPlaces) {
    d *= pow(10.,1.*decPlaces);
    d = std::floor(d+0.5);
    d /= pow(10.,1.*decPlaces);
    
    return d;
  }  


  // -------------------------------------------------------------------------------------
  static std::string toString(double d) {
    std::stringstream ss;
    ss << d;
    return ss.str();
  }
  
  
  // -------------------------------------------------------------------------------------
  static std::string toString(double d, int decPlaces) {
    std::stringstream ss;
    ss << round(d,decPlaces);
    return ss.str();
  }


  // -------------------------------------------------------------------------------------
  static TString toTString(double d) {
    return toString(d).c_str();
  }
  
  
  // -------------------------------------------------------------------------------------
  static TString toTString(double d, int decPlaces) {
    return toString(d,decPlaces).c_str();
  }
}
#endif

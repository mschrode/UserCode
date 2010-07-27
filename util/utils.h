// $Id: utils.h,v 1.1 2010/07/27 10:05:43 mschrode Exp $

#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <sstream>

#include "TString.h"


//!  Encapsulates useful classes and methods
//!
//!  \author   Matthias Schroeder (www.desy.de/~matsch)
//!  \date     2010/03/09
//!  $Id: utils.h,v 1.1 2010/07/27 10:05:43 mschrode Exp $
// -------------------------------------------------------------------------------------
namespace util {

  // -------------------------------------------------------------------------------------
  double round(double d, int decPlaces) {
    d *= pow(10.,1.*decPlaces);
    d = std::floor(d+0.5);
    d /= pow(10.,1.*decPlaces);
    
    return d;
  }  


  // -------------------------------------------------------------------------------------
  std::string toString(double d) {
    std::stringstream ss;
    ss << d;
    return ss.str();
  }
  
  
  // -------------------------------------------------------------------------------------
  std::string toString(double d, int decPlaces) {
    std::stringstream ss;
    ss << round(d,decPlaces);
    return ss.str();
  }


  // -------------------------------------------------------------------------------------
  TString toTString(double d) {
    return toString(d).c_str();
  }
  
  
  // -------------------------------------------------------------------------------------
  TString toTString(double d, int decPlaces) {
    return toString(d,decPlaces).c_str();
  }
}
#endif

// $Id: utils.h,v 1.4 2010/08/21 13:40:21 mschrode Exp $

#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <string>
#include <sstream>
#include <vector>

#include "TString.h"


//!  Encapsulates useful classes and methods
//!
//!  \author   Matthias Schroeder (www.desy.de/~matsch)
//!  \date     2010/03/09
//!  $Id: utils.h,v 1.4 2010/08/21 13:40:21 mschrode Exp $
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


  // -------------------------------------------------------------------------------------
  static bool findBin(double x, const std::vector<double> &binEdges, unsigned int &bin) {
    bin = 0;
    bool inRange = false;
    if( x >= binEdges.front() && x <= binEdges.back() ) {
      inRange = true;
      for(unsigned int i = 0; i < (binEdges.size()-1); ++i) {
	if( x > binEdges[i] ) bin = i;
	else break;
      }
    }
    
    return inRange;
  }
}
#endif

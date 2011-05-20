// $Id: utils.h,v 1.7 2011/05/17 16:36:43 mschrode Exp $

#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>

#include "TString.h"


//!  Encapsulates useful classes and methods
//!
//!  \author   Matthias Schroeder (www.desy.de/~matsch)
//!  \date     2010/03/09
//!  $Id: utils.h,v 1.7 2011/05/17 16:36:43 mschrode Exp $
// -------------------------------------------------------------------------------------
namespace util {

  //! Returns error of A / B for two uncorrelated variables
  //! A, B with error AE, BE
  // -------------------------------------------------------------------------------------
  static inline double ratioError(double A, double AE, double B, double BE) {
    return sqrt( AE*AE/B/B + BE*BE*A*A/B/B/B/B );
  }


  // -------------------------------------------------------------------------------------
  static inline double round(double d, int decPlaces) {
    d *= pow(10.,1.*decPlaces);
    d = std::floor(d+0.5);
    d /= pow(10.,1.*decPlaces);
    
    return d;
  }  


  // -------------------------------------------------------------------------------------
  static inline std::string toString(double d) {
    std::stringstream ss;
    ss << d;
    return ss.str();
  }
  
  
  // -------------------------------------------------------------------------------------
  static inline std::string toString(double d, int decPlaces) {
    std::stringstream ss;
    ss << round(d,decPlaces);
    return ss.str();
  }


  // -------------------------------------------------------------------------------------
  static inline TString toTString(double d) {
    return toString(d).c_str();
  }
  
  
  // -------------------------------------------------------------------------------------
  static inline TString toTString(double d, int decPlaces) {
    return toString(d,decPlaces).c_str();
  }


  // -------------------------------------------------------------------------------------
  static inline TString extractFileName(const TString &name) {
    TString fileName = name;
    if( fileName.Contains("/") ) {
      Ssiz_t pos = fileName.Last('/');
      fileName = fileName(pos+1,fileName.Length()-pos);
    }

    return fileName;
  }


  // -------------------------------------------------------------------------------------
  static inline bool findBin(double x, const std::vector<double> &binEdges, unsigned int &bin) {
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


  // -------------------------------------------------------------------------------------
  static inline TString jetAlgo(const TString &name) {
    TString algo = "";
    if( name.Contains("Calo") || name.Contains("calo") || name.Contains("CALO") ) algo = "Calo";
    else if( name.Contains("Jpt") || name.Contains("jpt") || name.Contains("JPT") ) algo = "JPT";
    else if( name.Contains("Pf") || name.Contains("pf") || name.Contains("PF") ) algo = "PF";

    return algo;
  }


  // For sorting jets in pt
  // --------------------------------------------------
  class JetIndexCol {
  private:
    class Jet {
    public:
      Jet(unsigned int jetIdx, double jetPt) : idx_(jetIdx), pt_(jetPt) {};
	const unsigned int idx_;
	const double pt_;
	// For sorting jets in pt
	static bool ptGreaterThan(const Jet *idx1, const Jet *idx2) {
	  // check for 0
	  if(idx1 == 0) {
	    return idx2 != 0;
	  } else if (idx2 == 0) {
	    return false;
	  } else {
	    return idx1->pt_ > idx2->pt_;
	  }
	}
    };

    std::vector<Jet*> jets_;


  public:
    JetIndexCol() {}
    ~JetIndexCol() { clear(); }
    
    unsigned int operator()(unsigned int i) { return idx(i); }
    unsigned int nJets() const { return jets_.size(); }
    unsigned int idx(unsigned int i) const { return jets_.at(i)->idx_; }
    double pt(unsigned int i) const { return jets_.at(i)->pt_; }
    
    void add(unsigned int jetIdx, double jetPt) {
      jets_.push_back(new Jet(jetIdx,jetPt));
    }
    void clear() {
      for(std::vector<Jet*>::iterator it = jets_.begin(); it != jets_.end(); ++it) {
	delete *it;
      }
      jets_.clear();
    }
    void sort() {
      std::sort(jets_.begin(),jets_.end(),Jet::ptGreaterThan);
    }
  };

}
#endif

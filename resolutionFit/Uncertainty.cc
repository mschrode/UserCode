// $Id: $

#include "Uncertainty.h"

#include <cassert>
#include <cmath>


namespace resolutionFit {
  //! Creates an uncertainty of zero with label
  //! "DefaultUncertainty"; this is a combined
  //! uncertainty.
  // --------------------------------------------
  Uncertainty::Uncertainty()
    : isCombined_(true), down_(0.), up_(0.), label_("DefaultUncertainty") {};


  //! Creates a combined uncertainty with label
  //! \p label. At construction, the upper and 
  //! lower uncertainties are zero. Sub-uncertainties 
  //! can be added by \p addUncertainty.
  // --------------------------------------------
  Uncertainty::Uncertainty(const TString &label)
    : isCombined_(true), down_(0.), up_(0.), label_(label) {};
    

  //! Creates a non-combined uncertainty with upper and lower
  //! value \p uncert and label \p label.
  // --------------------------------------------
  Uncertainty::Uncertainty(const TString &label, double uncert) 
    : isCombined_(false), down_(uncert), up_(uncert), label_(label) {};


  //! Creates a non-combined uncertainty with upper value \p up
  //! and lower value \p down and label \p label.
  // --------------------------------------------
  Uncertainty::Uncertainty(const TString &label, double up, double down)
    : isCombined_(false), down_(down), up_(up), label_(label) {};


  //! Deletes all sub-uncertainties
  // --------------------------------------------
  Uncertainty::~Uncertainty() {
    for(std::vector<Uncertainty*>::iterator it = uncerts_.begin();
	it != uncerts_.end(); it++) {
      delete *it;
    }
  }



  //! Adds \p uncert to the list of sub-uncertainties
  //! and increases \p nUncerts() by one. The upper and
  //! lower value will be added quardratically to the
  //! current values accessed by \p up() and \p down().
  //!
  //! \note This works only for combined uncertainties
  //! i.e. <tt>isCombined() == true<\tt>
  // --------------------------------------------
  void Uncertainty::addUncertainty(Uncertainty *uncert) {
    assert( isCombined() );
    uncerts_.push_back(uncert);
    down_ = sqrt( down()*down() + uncert->down()*uncert->down() );
    up_ = sqrt( up()*up() + uncert->up()*uncert->up() );
  }
}

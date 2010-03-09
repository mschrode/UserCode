#include "Uncertainty.h"

#include <cmath>

namespace resolutionFit {
  Uncertainty::Uncertainty()
    : down_(0.), up_(0.), label_("DefaultUncertainty") {};

  Uncertainty::Uncertainty(const TString &label)
    : down_(0.), up_(0.), label_(label) {};
    
  Uncertainty::Uncertainty(const TString &label, double uncert) 
    : down_(uncert), up_(uncert), label_(label) {};

  Uncertainty::Uncertainty(const TString &label, double up, double down)
    : down_(down), up_(up), label_(label) {};

  Uncertainty::~Uncertainty() {
    for(std::vector<Uncertainty*>::iterator it = uncerts_.begin();
	it != uncerts_.end(); it++) {
      delete *it;
    }
  }

  void Uncertainty::addUncertainty(Uncertainty *uncert) {
    uncerts_.push_back(uncert);
    down_ = sqrt( down()*down() + uncert->down()*uncert->down() );
    up_ = sqrt( up()*up() + uncert->up()*uncert->up() );
  }
}

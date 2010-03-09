#ifndef UNCERTAINTY_H
#define UNCERTAINTY_H

#include <cassert>
#include <vector>

#include "TString.h"

namespace resolutionFit {
  class Uncertainty {
  public: 
    Uncertainty();
    Uncertainty(const TString &label);
    Uncertainty(const TString &label, double uncert);
    Uncertainty(const TString &label, double up, double down);
      
    ~Uncertainty();
      
    TString label() const { return label_; }
    double up() const { return up_; }
    double down() const { return down_; }

    void addUncertainty(Uncertainty *uncert);

    bool isCombined() const { return nUncerts()==0 ? false : true; }
    int nUncerts() const { return static_cast<int>(uncerts_.size()); }
    const Uncertainty *uncert(int i) const { assert( i>=0 && i < nUncerts() ); return uncerts_[i]; }
    TString label(int i) const { assert( i>=0 && i < nUncerts() ); return uncerts_[i]->label(); }
    double up(int i) const { assert( i>=0 && i < nUncerts() ); return uncerts_[i]->up(); }
    double down(int i) const { assert( i>=0 && i < nUncerts() ); return uncerts_[i]->down(); }

  private:
    TString label_;
    double down_;
    double up_;

    std::vector<Uncertainty*> uncerts_;
  };
}
#endif

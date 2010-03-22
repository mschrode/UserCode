// $Id: $

#ifndef UNCERTAINTY_H
#define UNCERTAINTY_H

#include <cassert>
#include <vector>

#include "TString.h"

namespace resolutionFit {

  //! \brief An upper and lower uncertainty value
  //!
  //! An \p Uncertainty has an upper and a lower
  //! value as well as a label defining the
  //! uncertainty. It can be combined of several
  //! sub-uncertainties; in that case the upper
  //! and lower value are the quadratic sum of
  //! the respective values of the uncertainties
  //! it is composed of.
  //!
  //! \author Matthias Schroeder
  //! \date 2009/03/08
  //! $Id: $
  // --------------------------------------------
  class Uncertainty {
  public: 
    //! Default constructor
    Uncertainty();
    //! Constructor for combined uncertainty
    Uncertainty(const TString &label);
    //! Constructor for non-combined symmetric uncertainty
    Uncertainty(const TString &label, double uncert);
    //! Constructor for non-combined asymmetric uncertainty
    Uncertainty(const TString &label, double up, double down);
      
    //! Destructor
    ~Uncertainty();

    //! Returns the label      
    TString label() const { return label_; }
    //! \brief Returns the upper value of the uncertainty
    //!
    //! In case of a combined uncertainty, this is the
    //! quadratic sum of the upper values of the
    //! sub-uncertainties.
    double up() const { return up_; }
    //! \brief Returns the lower value of the uncertainty
    //!
    //! In case of a combined uncertainty, this is the
    //! quadratic sum of the lower values of the
    //! sub-uncertainties.
    double down() const { return down_; }

    //! Adds an uncertainty as sub-uncertainty
    void addUncertainty(Uncertainty *uncert);

    //! Returns \p true in case of a combined uncertainty
    bool isCombined() const { return isCombined_; }
    //! Returns number of sub-uncertainties
    int nUncerts() const { return static_cast<int>(uncerts_.size()); }
    //! Returns the i-th sub-uncertainty
    const Uncertainty *uncert(int i) const { assert( i>=0 && i < nUncerts() ); return uncerts_[i]; }
    //! Returns the label of the i-th sub-uncertainty
    TString label(int i) const { assert( i>=0 && i < nUncerts() ); return uncerts_[i]->label(); }
    //! Returns the upper value of the i-th sub-uncertainty
    double up(int i) const { assert( i>=0 && i < nUncerts() ); return uncerts_[i]->up(); }
    //! Returns the lower value of the i-th sub-uncertainty
    double down(int i) const { assert( i>=0 && i < nUncerts() ); return uncerts_[i]->down(); }


  private:
    const bool isCombined_;	//! True for combined uncertainties

    TString label_;		//! Label of the uncertainty
    double down_;		//! Lower value
    double up_;			//! Upper value

    std::vector<Uncertainty*> uncerts_; //! Sub-uncertainties
  };
}
#endif

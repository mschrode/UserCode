// $Id: SmearPhotonJet.cc,v 1.2 2010/01/21 16:49:26 mschrode Exp $

#include "SmearPhotonJet.h"



//!  \brief Get the negative log-likelihood of this event
//!  \return The negative log-likelihood of this event
// --------------------------------------------------
double SmearPhotonJet::chi2() const {
  double respPD = respPDF( GetMess()->pt / GetTruth(), GetTruth() );
  return -1. * GetWeight() * log( respPD / GetTruth() ); // Need to divide by truth to have probability (!= density)
}



//!  \brief Get the negative log-likelihood and the
//!         contributions to the 1. and 2. derivatives
//!
//!  Calculates the negative log-likelihood \f$ -\ln L \f$
//!  of this event. Moreover the contribution to the 
//!  first and second derivative ('temp_derivative1',
//!  'temp_derivative2') of the global \f$ -\ln L \f$
//!  function is calculated numerically and returned
//!  by reference, where
//!  \f[
//!   \frac{\partial (-\ln L) }{\partial p}
//!   = \sum \frac{\textrm{temp\_derivative1}}{2\epsilon}
//!  \f]
//!  and
//!  \f[
//!   \frac{\partial^{2} (-\ln L) }{\partial p^{2}}
//!   = \sum \frac{\textrm{temp\_derivative2}}{\epsilon^{2}}
//!  \f]
//!
//!  \param temp_derivative1 Pointer to first derivative contribution
//!  \param temp_derivative2 Pointer to second derivative contribution
//!  \param epsilon Step size \f$ \epsilon \f$ for derivative calculation
//!  \return The negative log-likelihood of this event
// --------------------------------------------------
double SmearPhotonJet::chi2_fast(double * temp_derivative1,
			     double * temp_derivative2,
			     double const epsilon) const {
  double f = chi2();

  double   oldpar;
  double   temp1;
  double   temp2;

  // Vary parameters of response pdf
  for(int i = 0; i < pdf_.nRespPars(); i++) {
    oldpar = pdf_.respPar()[i];

    pdf_.respPar()[i] += epsilon;
    temp1 = chi2();

    pdf_.respPar()[i] -= 2.*epsilon;
    temp2 = chi2();

    pdf_.respPar()[i] = oldpar;

    temp_derivative1[pdf_.respParIdx()+i] += temp1 - temp2;
    temp_derivative2[pdf_.respParIdx()+i] += temp1 + temp2 - 2*f;
  }

  return f;
}

//
// $Id: CalibData.cc,v 1.32 2009/11/24 17:13:17 stadie Exp $
//
#include "CalibData.h"

#include <map>
#include <cmath>

double (*Event::ScaleResidual)(double z2) = &Event::ScaleNone;

//!  Scale the normalized, squared residual
//!  \f$ z^{2} = \chi^{2}/\textrm{weight} \f$
//!  using the Cauchy-Function
//!  \f$ z^{2} \rightarrow c^{2}\ln( 1 + (z/c)^{2} ) \f$
//!
//!  \param z2 Normalized and squared residual
//!  \return Scaled residual
double Event::ScaleCauchy(double const z2)
{
  double const c = 2.3849;
  return (c*c) * log( 1 + z2*(1.0/(c*c)) );
}

//!  Scale the normalized, squared residual
//! \f$ z^{2} = \chi^{2}/\textrm{weight} \f$
//!  using the Huber-Function
//!  \f[
//!  z^{2} \rightarrow 
//!  \left\{
//!     \begin{array}{ll}
//!        z & \textrm{for } |z| <= c \\ c ( 2|z| - c ) & \textrm{for } |z| > c
//!      \end{array}
//!    \right.
//!  \f]
//!
//!  \param z2 Normalized and squared residual
//!  \return Scaled residual
double Event::ScaleHuber(double const z2)
{
  static double const c = 1.345;
  double const z = sqrt(z2);
  return (  std::abs(z) <= c  ?  z2  :  c*(2.*std::abs(z) - c)  );
}


//!  \brief Cut on residuals
//!
//!  discards events with $|residual| > sqrt(c2) \sigma$
//!
//!  \param z2 Normalized and squared residual
//!  \return Scaled residual
double Event::ScaleTukey(const double z2)
{
  const double c2 = 16;
  if(z2 > c2) return 0;
  double w = 1-z2/c2;
  return w*w * z2;
}


unsigned int TAbstractData::total_n_pars = 0;


double TData_ParLimit::chi2_fast(double* temp_derivative1, 
 				 double* temp_derivative2, double const epsilon) const {
  // Penalty term with current parameter values
  double new_chi2  = chi2();
 
  // Variation of parameters
  double oldpar    = _par[0];
  _par[0]         += epsilon;
  double temp2     = chi2();
  _par[0]          = oldpar - epsilon;
  double temp1     = chi2();
 
  // Difference of chi2 at par+epsilon and par-epsilon
  temp_derivative1[_index] += (temp2 - temp1);                // for 1st derivative
  temp_derivative2[_index] += (temp2 + temp1 - 2.*new_chi2);  // for 2nd derivative
 
  // Reset original parameter value
  _par[0]  = oldpar;
 
  return new_chi2;
}

std::vector<TAbstractData*> TData_ParLimit::_cache;

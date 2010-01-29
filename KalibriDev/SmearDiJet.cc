// $Id: SmearDiJet.cc,v 1.2 2010/01/21 16:49:24 mschrode Exp $

#include "SmearDiJet.h"

//!  \brief Constructor
//!  \param jet1 First jet
//!  \param jet2 Second jet
//!  \param jet3 Third jet (for cuts)
//!  \param weight Event weight
//!  \param respPDF Response probability density
//!  \param truthPDF True pt probability density
//!  \param min Minimum true pt in integration range
//!  \param max Maximum true pt in integration range
//!  \param eps Integration precision for convergence
//!  \param niter Maximum number of iterations in integration
// --------------------------------------------------
SmearDiJet::SmearDiJet(Jet * jet1,
		       Jet * jet2,
		       Jet * jet3,
		       double ptHat,
		       double weight,
		       const Function& respPDF,
		       const Function& truthPDF,
		       double min,
		       double max,
		       double eps,
		       int niter)
  : SmearData(TypeSmearDiJet,jet1,0,ptHat,weight,respPDF),
    kMaxNIter_(niter),
    kEps_(eps),
    kMin_(min),
    kMax_(max),
    jet2_(jet2),
    jet3_(jet3),
    truthPDF_(truthPDF) { };



// --------------------------------------------------
SmearDiJet::~SmearDiJet() { 
  delete jet2_;
  delete jet3_;
}



// --------------------------------------------------
void SmearDiJet::ChangeParAddress(double* oldpar, double* newpar) {
  respPDF_.changeParBase(oldpar,newpar);
  truthPDF_.changeParBase(oldpar,newpar);
}



//!  \brief Get the negative log-likelihood of this event
//!
//!  Calculates the probability of this event configuration
//!  \f[
//!   P(m_{1},m_{2}) = \int\;dt\;f(t)r(m_{1}/t)r(m_{2}/t),
//!  \f]
//!  where \f$ f \f$ is the truth pdf and \f$ r \f$ is the
//!  response pdf. \f$ r \f$ is normalized to unity,
//!  \f$ f \f$ is normalized such that \f$ P \f$ is 
//!  normalized to unity (see truthPDF(t)). The method
//!  returns \f$ -\ln(P) \f$.
//!
//!  The integral is calculated numerically using the
//!  Simpson's 3/8 rule.
//!
//!  \return The negative log-likelihood of this event
// --------------------------------------------------
double SmearDiJet::chi2() const
{
  double         h         = kMax_ - kMin_;     // Integration interval
  double         pint      = 0.;              // Current value of integral over response pdfs
  double         pint_old  = 1.;              // Value of integral over response pdfs from previous iteration
  int            nIter     = 0;               // Current iteration in interval splitting
  std::vector<double> pp;         // Product of current function values of response pdfs
  std::vector<double> pp_old;     // Product of function values of response pdfs from previous iteration

  // Iterate until precision or max. number iterations reached
  while(fabs((pint - pint_old) / pint_old) > kEps_ && nIter < kMaxNIter_) {
    pint_old = pint;
    pint     = 0;
    pp_old   = pp;
    pp.clear();
    h       /= 3.;    // In each iteration, split h into 3 new intervals
    
    // Loop over nodes xi i.e. interval borders
    for(int i = 0; i <= pow(3.0,nIter+1); ++i){
      double t = kMin_ + i * h;  // Pt at node
      
      // Calculate probability only at new nodes
      if(nIter == 0 || i % 3 != 0) {
	double p0 = respPDF( jet1()->pt() / t, t );
	double p1 = respPDF( jet2()->pt() / t, t );
	pp.push_back(p0 * p1 * truthPDF(t)); // Store product of pdfs and normalization
      } else {
	pp.push_back(pp_old.at(i/3));       // Store product of pdfs previously calcluated
      }
    }

    // Sum up weighted function values
    for(unsigned int i = 0; i < pp.size(); i++)	{
      double w = 1.;                       // Weight w from Simpson's rule
      if( i > 0 && i < (pp.size() - 1) ) { // w = 1 for x0 and last node
	if( i % 3 == 0 ) {                 // w = 2 for x3, x6, ...
	  w = 2.;
	} else {
	  w = 3.;
	}
      }
      pint += w * (pp.at(i));              // Sum up weighted function values
    }
    pint *= (3. * h / 8.);                 // Apply overall normalization
    nIter++;
  }

  if( pint <= 0 ) return 0.;

  return  -1. * GetWeight() * log(pint);
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
double SmearDiJet::chi2_fast(double * temp_derivative1,
			     double * temp_derivative2,
			     double const epsilon) const {
  double f = chi2();

  int      idx;
  double * par;
  double   oldpar;
  double   temp1;
  double   temp2;

  // Vary parameters of response pdf
  idx = respPDF_.parIndex();
  par = respPDF_.firstPar();
  for(int i = 0; i < respPDF_.nPars(); i++) {
    oldpar = par[i];
    par[i] += epsilon;
    temp1   = chi2();

    par[i] -= 2.*epsilon;
    temp2   = chi2();

    par[i]  = oldpar;

    temp_derivative1[idx+i] += temp1 - temp2;
    temp_derivative2[idx+i] += temp1 + temp2 - 2*f;
  }

  // Vary parameters of truth pdf
  idx = truthPDF_.parIndex();
  par = truthPDF_.firstPar();
  for(int i = 0; i < truthPDF_.nPars(); i++) {
    oldpar = par[i];

    par[i] += epsilon;
    temp1   = chi2();

    par[i] -= 2.*epsilon;
    temp2   = chi2();

    par[i]  = oldpar;

    temp_derivative1[idx+i] += temp1 - temp2;
    temp_derivative2[idx+i] += temp1 + temp2 - 2*f;
  }

  return f;
}



//!  \brief Truth pdf
//!  \note The probability is not normalized to unity but in such
//!        a way that the negative log-likelihood as returned by
//!        chi2() is properly normalized.
//!  \param t Truth, pt in GeV
//!  \return The probability density of the truth \p t,
//!          normalized such that the negative log-likelihood
//!          is normalized
// --------------------------------------------------
double SmearDiJet::truthPDF(double t) const {
  Measurement meas;
  meas.pt = t;
  return truthPDF_(&meas);
}

// --------------------------------------------------
double SmearDiJet::truthPDFSigma(double t) const {
  Measurement meas;
  meas.pt = t;
  return truthPDF_.sigma(&meas);
}


//!  \brief Print event parameters
// --------------------------------------------------
void SmearDiJet::printInitStats() const {
  std::cout << "Event type: " << GetType() << "\n";
  std::cout << "Integration parameters\n";
  std::cout << "  niter: " << kMaxNIter_ << "\n";
  std::cout << "  eps:   " << kEps_ << "\n";
  std::cout << "  range: " << kMin_ << " < truth < " << kMax_ << " (GeV)\n";
}

//
//  $Id: Parametrization.cc,v 1.1 2010/01/21 16:46:05 mschrode Exp $
//
#include "Parametrization.h"


#include "TH1D.h"
#include "TMath.h"



//!  \param tMin          Minimum of non-zero part of truth pdf
//!  \param tMax          Maximum of non-zero part of truth pdf
//!  \param rMin          Minimum of binned part of response pdf \f$ H_{inter} \f$
//!  \param rMax          Maximum of binned part of response pdf \f$ H_{inter} \f$
//!  \param ptDijetMin    Minimum of dijet pt
//!  \param ptDijetMax    Maximum of dijet pt
//!  \param rNBins        Number of bins of binned part of response pdf \f$ H_{inter} \f$
//!  \param rParScales    Jet parameter scales
//!  \param gaussPar   Parameters for mean of Gaussian
// ------------------------------------------------------------------------
SmearStepGaussInter::SmearStepGaussInter(double tMin, double tMax, double rMin, double rMax, int rNBins, double ptDijetMin, double ptDijetMax, const std::vector<double>& rParScales, const std::vector<double>& gaussPar)
  : Parametrization(0,rNBins+2,0,0),
    tMin_(tMin),
    tMax_(tMax),
    rMin_(rMin),
    rMax_(rMax),
    ptDijetMin_(ptDijetMin),
    ptDijetMax_(ptDijetMax),
    nStepPar_(rNBins),
    nGaussPar_(2),
    binWidth_((rMax_ - rMin_)/nStepPar_),
    respParScales_(rParScales),
    gaussPar_(gaussPar) {
  for(int i = 0; i < nStepPar_+1; i++) {
    binCenters_.push_back( rMin_ + (0.5+i)*binWidth_ );
  }
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= rMin_ && rMin_ < rMax_ );
  assert( 0.0 <= ptDijetMin_ && ptDijetMin_ < ptDijetMax_ );
  assert( respParScales_.size() >= nJetPars() );
  assert( gaussPar_.size() >= 1 );
  
  print();
  
  // Integral over dijet resolution for truth pdf
  ptDijetCutInt_ = new TH1D("norm","",400,tMin_,tMax_);
}



// ------------------------------------------------------------------------
SmearStepGaussInter::~SmearStepGaussInter() { 
  binCenters_.clear(); 
  delete ptDijetCutInt_; 
}



//!  The truth pdf contains an integral over the dijet response
//!  describing the cuts on ptdijet:
//!  \f[
//!   \int^{x_{1}}_{x_{0}}dx\,\frac{r(x/t)\cdot t}{\sqrt{2}}
//!  \f]
//!  The values of this integral for different truth \p t are
//!  stored in the histogram \p ptDijetCutInt_.
//!
//!  Calling this function recalculates the integrals using the current
//!  parameter values for the response function \p r. The integrals
//!  are calculated numerically by summing over 400 bins between
//!  \p ptDijetMin_ and \p ptDijetMax_.
//!  \sa correctedGlobalJetEt()
// ------------------------------------------------------------------------
void SmearStepGaussInter::update(const double * par) {
  std::cout << "'" << name() << "': Updating ptDijet cut integral... ";
  ptDijetCutInt_->Reset();
    
  for(int bin = 1; bin < ptDijetCutInt_->GetNbinsX(); bin++) {
    Measurement x;
    x.pt = ptDijetCutInt_->GetBinCenter(bin);
    double integral = 0.;
    int nSteps = 400;
    double dPtDijet = (ptDijetMax_ - ptDijetMin_) / nSteps;
    for(int i = 0; i < nSteps; i++) {
      double ptDijet = ptDijetMin_ + i*dPtDijet;
      x.E = ptDijet / x.pt;
      double prob = correctedJetEt(&x,par);
      prob *= x.pt;
      integral += prob;
    }
    integral /= sqrt(2.);
    integral *= dPtDijet;
    ptDijetCutInt_->SetBinContent(bin,integral);
  }
  std::cout << "ok\n";
}



//!  \param x   Measurement::E is the response, Measurement::pt the true pt
//!  \param par Pointer to parameters (parameters are multiplied by the corresponding scale)
//!  \return Probability density of response
// ------------------------------------------------------------------------
double SmearStepGaussInter::correctedJetEt(const Measurement *x,const double *par) const {
  // Probability density from Gaussian part
//   double c = respParScales_.at(0)*par[0];
//   double u = gaussPar_.at(0);
//   double a1 = respParScales_.at(1)*par[1];
//   double a2 = respParScales_.at(2)*par[2];
//   double a3 = respParScales_.at(3)*par[3];


  double u = gaussPar_.at(0);
  double c = respParScales_.at(0)*par[0];
  double s = respParScales_.at(1)*par[1];

  double p = c / sqrt(2* M_PI ) / s * exp(-pow((x->E - u) / s, 2) / 2);
  // Probability density from step function part
  // Copy parameter values to temp vector for normalization
  std::vector<double> stepPar(nStepPar_,1.);
  double norm = 0.;
  for(int i = 0; i < nStepPar_; i++) {
    stepPar.at(i) = respParScales_.at(nGaussPar_+i)*par[nGaussPar_+i];
    if( stepPar.at(i) < 0. ) stepPar.at(i) = 0.;
    norm += stepPar.at(i);
  }
  norm = (1.-c)/norm/binWidth_;
    
  p += norm*interpolate(x->E,stepPar);
    
  return p;
}



// ------------------------------------------------------------------------
double SmearStepGaussInter::correctedJetEtSigma(const Measurement *x,const double *par,const double *cov, const std::vector<int> &covIdx) const {
  // Store derivatives
  std::vector<double> df(nJetPars());
  for(size_t i = 0; i < nJetPars(); i++) {
    df[i] = respDerivative(x,par,i);
  }

  // Calculate variance
  double var = 0.;
  if( x->E > 0.99 && x->E < 1. ) std::cout << std::endl;
  for(int i = 0; i < static_cast<int>(nJetPars()); i++) { // Outer loop over parameters
    for(int j = 0; j < i+1; j++) { // Inner loop over parameters
      int idx = (i*i + i)/2 + j; // Index of (i,j) in covariance vector
      if( cov[idx] ) {
	if( i == j ) { // Diagonal terms
	  var += df[i]*df[i]*respParScales_[i]*respParScales_[i]*cov[idx];
	} else { // Off-diagonal terms
	  var += 2*df[i]*df[j]*respParScales_[i]*respParScales_[j]*cov[idx];
	}
      }
    } // End of inner loop over parameters
  } // End of outer loop over parameters
  df.clear();
  // Return standard deviation
  return sqrt(var);
}



//!         of dijet probability (see also \p SmearDiJet::truthPDF(t)).
//!  \param x   \p Measurement::pt is truth for which the 
//!             probability density is returned
//!  \param par Pointer to parameters (parameters are multiplied by the corresponding scale)
//!  \return Probability density of truth times normalization of dijet probability
// ------------------------------------------------------------------------
double SmearStepGaussInter::correctedGlobalJetEt(const Measurement *x,const double *par) const {
//   // Norm of probability of dijet event configuration
//   double norm = 0.;
//   for(int bin = 1; bin < ptDijetCutInt_->GetNbinsX(); bin++) {
//     double pt = ptDijetCutInt_->GetBinCenter(bin);
//     norm += pow(pt,2-par[0]) * ptDijetCutInt_->GetBinContent(bin) * ptDijetCutInt_->GetXaxis()->GetBinWidth(1);
//   }
  
//   double p = 0.;
//   if( norm ) {
//     p = truthPDF(x->pt,par[0]);
//     p /= norm;
//   }
  
//   return p;

  double norm = (pow(tMax_,3)-pow(tMin_,3))/3.;
  return x->pt > tMin_ && x->pt < tMax_ ? 1./norm : 0.;
}



// ------------------------------------------------------------------------
double SmearStepGaussInter::respDerivative(const Measurement *x, const double *par, int i) const {
  double r = x->E;

  // Probability density from Gaussian part
  double c = respParScales_.at(0)*par[0];
  double u = gaussPar_.at(0);
//   double a1 = respParScales_.at(1)*par[1];
//   double a2 = respParScales_.at(2)*par[2];
//   double a3 = respParScales_.at(3)*par[3];
//   double s = sqrt( a1*a1/x->pt/x->pt + a2*a2/x->pt + a3*a3 );
  double s = respParScales_[1]*par[1];
  double G = 1. / sqrt(2* M_PI ) / s * exp(-pow((r - u) / s, 2) / 2);
        
  // Probability density from step function part
  // Copy parameter values to temp vector for normalization
  std::vector<double> stepPar(nStepPar_,1.);
  double norm = 0.;
  for(int k = 0; k < nStepPar_; k++) {
    stepPar.at(k) = respParScales_.at(nGaussPar_+k)*par[nGaussPar_+k];
    if( stepPar.at(k) < 0. ) stepPar.at(k) = 0.;
    norm += stepPar.at(k);
  }
  norm = 1./norm/binWidth_;
  double S = norm*interpolate(r,stepPar);

  double df = 0.;
  // df/dc
  if( i == 0 ) {
    df = G - S;
  }
//   // df/da1
//   else if( i == 1 ) {
//     df = c*G*((r - u)*(r - u)/s/s - 1.)*a1/(x->pt)/(x->pt)/s/s;
//   }
//   // df/da2
//   else if( i == 2 ) {
//     df = c*G*((r - u)*(r - u)/s/s - 1.)*a2/(x->pt)/s/s;
//   }
//   // df/da3
//   else if( i == 3 ) {
//     df = c*G*((r - u)*(r - u)/s/s - 1.)*a3/s/s;
//   }
  // df/dsigma
  else if( i == 1 ) {
    df = c*G*((r - u)*(r - u)/s/s - 1.)/s;
  }
  // df/dbi
  else if( i < static_cast<int>(nJetPars()) ) {
    if( r > rMin_ - 0.5*binWidth_  &&  r < rMax_ + 0.5*binWidth_ ) {
      // Find index j of the next larger bin center
      int j = 0;                              
      while( r > binCenters_[j] ) j++;

      double dr = (binCenters_[j] - r)/binWidth_;
      if( i-nGaussPar_ == j ) {
	df = 1. - dr - binWidth_*S;
      } else if( i-nGaussPar_ == j-1) {
	df = dr - binWidth_*S;
      }
      df *= (1.-c)*norm;
    }
  }

  return df;
}



//!  Returns the linearly weighted mean value of two adjacent bin
//!  contents, where the two bins are the ones whose bin centers
//!  are closest to \p r. The contents are weighted with the distance
//!  of r from the corresponding bin center. The interpolation
//!  assumes \p nStepPar_ bins. At the left and right edges of the
//!  binned range, the interpolation is done assuming 0 bin content
//!  outside the binned range.
//!
//!  \param r   Response
//!  \param par Parameters / bin content to be interpolated
//!  \return    Interpolated bin content
// ------------------------------------------------------------------------
double SmearStepGaussInter::interpolate(double r, const std::vector<double>& par) const {
  double p = 0.;     // The interpolated parameter

  // Check that r is in range of binning +/- 0.5*mBinWidth
  if( r > rMin_ - 0.5*binWidth_  &&  r < rMax_ + 0.5*binWidth_ ) {
    // Find index i of the next larger bin center
    unsigned int i = 0;                              
    while( r > binCenters_.at(i) ) i++;
    assert( i < binCenters_.size() );

    double dx  = binCenters_.at(i) - r;  // Distance to next larger bin center
    dx        /= binWidth_;              // dx relative to bin width
    double a   = 0.;
    double b   = 0.;

    // Find bin contents to be interpolated
    if( i == 0 ) { // Left edge
      a = 0.;
      b = par.front();
    }
    else if( i == binCenters_.size()-1 ) { // Right edge
      a = par.back();
      b = 0.;
    }
    else { // Central part
      a = par.at(i-1);
      b = par.at(i);
    }

    // Interpolate contents
    p = b - (b-a)*dx;
  }

  return p;
}



//!  The probability density of \p pt if cuts
//!  \f$ \texttt{ptDijetMin\_} < p^{\textrm{dijet}}_{T} < \texttt{ptDijetMax\_} \f$
//!  are applied:
//!  \f[ t^{-n} \int^{\texttt{ptDijetMax\_}}_{\texttt{ptDijetMin\_}}
//!      dx\,\frac{r(x/t) \cdot t}{\sqrt{2}} \f]
// ------------------------------------------------------------------------
double  SmearStepGaussInter::truthPDF(double pt, double n) const {
//   double prob = 0.;
//   if( tMin_ < pt && pt < tMax_ ) {
//     prob = 1. / pow( pt, n );
//     int bin = ptDijetCutInt_->FindBin(pt);
//     prob *= ptDijetCutInt_->GetBinContent(bin);
//   }
  
//   return prob;

  return pt > tMin_ && pt < tMax_ ? 1./(tMax_-tMin_) : 0.;
}



// ------------------------------------------------------------------------
void SmearStepGaussInter::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  Probability density of ptTrue spectrum:\n";
  std::cout << "    Powerlaw\n";
  std::cout << "    " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
  std::cout << "    " << ptDijetMin_ << " < ptDijet < " << ptDijetMax_ << " GeV\n";
  std::cout << "  Probability density of response:\n";
  std::cout << "    Non-zero for " << rMin_ << " < R < " << rMax_ << "\n";
  std::cout << "    " << nStepPar_ << " step parameters\n";
  std::cout << "    " << gaussPar_.size() << " Gauss parameters: ";
  for(size_t i = 0; i < gaussPar_.size(); i++) {
    std::cout << gaussPar_.at(i);
    if( i != gaussPar_.size() - 1 ) std::cout << ", ";
  }
  std::cout << std::endl;
}




// ------------------------------------------------------------------------
SmearCrystalBall::SmearCrystalBall(double tMin, double tMax, double rMin, double rMax, double ptDijetMin, double ptDijetMax, const std::vector<double>& rParScales)
  : Parametrization(0,4,0,0),
    tMin_(tMin),
    tMax_(tMax),
    rMin_(rMin),
    rMax_(rMax),
    ptDijetMin_(ptDijetMin),
    ptDijetMax_(ptDijetMax),
    respParScales_(rParScales) {
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= rMin_ && rMin_ < rMax_ );
  assert( 0.0 <= ptDijetMin_ && ptDijetMin_ < ptDijetMax_ );
  assert( respParScales_.size() >= nJetPars() );
  
  print();
  
  // Integral over dijet resolution for truth pdf
  ptDijetCutInt_ = new TH1D("norm","",400,tMin_,tMax_);
}



// ------------------------------------------------------------------------
SmearCrystalBall::~SmearCrystalBall() { 
  delete ptDijetCutInt_; 
}



//!  The truth pdf contains an integral over the dijet response
//!  describing the cuts on ptdijet:
//!  \f[
//!   \int^{x_{1}}_{x_{0}}dx\,\frac{r(x/t)\cdot t}{\sqrt{2}}
//!  \f]
//!  The values of this integral for different truth \p t are
//!  stored in the histogram \p ptDijetCutInt_.
//!
//!  Calling this function recalculates the integrals using the current
//!  parameter values for the response function \p r. The integrals
//!  are calculated numerically by summing over 400 bins between
//!  \p ptDijetMin_ and \p ptDijetMax_.
//!  \sa correctedGlobalJetEt()
// ------------------------------------------------------------------------
void SmearCrystalBall::update(const double * par) {
  std::cout << "'" << name() << "': Updating ptDijet cut integral... ";

  ptDijetCutInt_->Reset();
  for(int bin = 1; bin < ptDijetCutInt_->GetNbinsX(); bin++) {
    Measurement x;
    x.pt = ptDijetCutInt_->GetBinCenter(bin);
    double integral = 0.;
    int nSteps = 400;
    double dPtDijet = (ptDijetMax_ - ptDijetMin_) / nSteps;
    for(int i = 0; i < nSteps; i++) {
      double ptDijet = ptDijetMin_ + i*dPtDijet;
      x.E = ptDijet / x.pt;
      double prob = correctedJetEt(&x,par);
      prob *= x.pt;
      integral += prob;
    }
    integral /= sqrt(2.);
    integral *= dPtDijet;
    ptDijetCutInt_->SetBinContent(bin,integral);
  }
  std::cout << "ok\n";
}



//! The response pdf parametrization is a crystall ball function
//! with the parameters
//! - 0: mean response
//! - 1: sigma
//! - 2: alpha
//! - 3: n
//! The pdf is 0 for \f$ R < rMin\_ \f$ and \f$ R > rMax\_ \f$.
//! \param x   Measurement::E is the response, Measurement::pt the true pt
//! \param par Pointer to parameters (parameters are multiplied by the corresponding scale)
//! \return Probability density of response
// ------------------------------------------------------------------------
double SmearCrystalBall::correctedJetEt(const Measurement *x,const double *par) const {
  double p = 0.;

  //if( x->E > rMin_ && x->E < rMax_ ) {
    double mean = respParScales_[0]*par[0];
    double sigma = respParScales_[1]*par[1];
    double alpha = respParScales_[2]*par[2];
    double n = respParScales_[3]*par[3];

    if( sigma <= 0. ) sigma = 1E-5;
    if( alpha <= 0. ) alpha = 1E-5;
    if( n <= 0. ) n = 1. + 1E-5;
    
    //if( n <= alpha*alpha ) n = alpha*alpha + 1E-5;
    
//     double norm = 0.;
//     int nSteps = 100;
//     double deltaResp = (rMax_-rMin_)/nSteps;
//     for(int i = 0; i < nSteps; i++) {
//       norm += crystallBallFunc(rMin_+deltaResp*i,mean,sigma,alpha,n);
//     }

//     if( norm > 0 ) p = crystallBallFunc(x->E,mean,sigma,alpha,n)/norm/deltaResp;

    p = crystalBallFunc(x->E,mean,sigma,alpha,n);
    //}

  return p;
}

double SmearCrystalBall::correctedJetEtSigma(const Measurement *x,const double *par,const double *cov, const std::vector<int> &covIdx) const {
  // Scale parameter values
  double mean = respParScales_[0]*par[0];
  double sigma = respParScales_[1]*par[1];
  double alpha = respParScales_[2]*par[2];
  double n = respParScales_[3]*par[3];

  // Check parameter limits
  if( sigma <= 0 ) sigma = 1E-5;
  if( alpha <= 0 ) alpha = 1E-5;
  if( n <= 1 ) n = 1. + 1E-5;
  //if( n <= alpha*alpha ) n = alpha*alpha + 1E-5;

  // Store derivatives
  std::vector<double> df(nJetPars());
  for(size_t i = 0; i < nJetPars(); i++) {
    df[i] = crystalBallDerivative(x->E,mean,sigma,alpha,n,i);
  }

  // Calculate variance
  double var = 0.;
  for(int i = 0; i < static_cast<int>(nJetPars()); i++) { // Outer loop over parameters
    for(int j = 0; j < i+1; j++) { // Inner loop over parameters
      int idx = (i*i + i)/2 + j; // Index of (i,j) in covariance vector
      if( cov[idx] ) {
	if( i == j ) { // Diagonal terms
	  var += df[i]*df[i]*respParScales_[i]*respParScales_[i]*cov[idx];
	} else { // Off-diagonal terms
	  var += 2*df[i]*df[j]*respParScales_[i]*respParScales_[j]*cov[idx];
	}
      }
    } // End of inner loop over parameters
  } // End of outer loop over parameters
  // Return standard deviation
  return sqrt(var);
}



//!  \param x   \p Measurement::pt is truth for which the 
//!             probability density is returned
//!  \param par Pointer to parameters (parameters are multiplied by the corresponding scale)
//!  \return Probability density of truth times normalization of dijet probability
// ------------------------------------------------------------------------
double SmearCrystalBall::correctedGlobalJetEt(const Measurement *x,const double *par) const {
  // Norm of probability of dijet event configuration
//   double norm = 0.;
//   for(int bin = 1; bin < ptDijetCutInt_->GetNbinsX(); bin++) {
//     double pt = ptDijetCutInt_->GetBinCenter(bin);
//     //    norm += pow(pt,2-par[0]) * ptDijetCutInt_->GetBinContent(bin) * ptDijetCutInt_->GetXaxis()->GetBinWidth(1);
//     norm += pow(pt,2-par[0]) * ptDijetCutInt_->GetXaxis()->GetBinWidth(bin);
//   }
  
//   double p = 0.;
//   if( norm ) {
//     p = truthPDF(x->pt,par[0]);
//     p /= norm;
//   }
 
//  return p;

  
  double norm = (pow(tMax_,3)-pow(tMin_,3))/3.;
  return 1./norm;
}



// ------------------------------------------------------------------------
double SmearCrystalBall::crystalBallFunc(double x, double mean, double sigma, double alpha, double n) const {
  double f = 0.;

  if( x > rMin_ && x < rMax_ ) {
    double u = (x - mean)/sigma;
    double A = pow(n/alpha,n)*exp(-0.5*alpha*alpha);
    double B = n/alpha - alpha;
    
    if( u > -alpha ) {             // Gaussian part
      f = exp(-0.5*u*u);
    } else {                       // Powerlaw part
      f =  A*pow(B-u,-n);
    }
  }

  return f/crystalBallNorm(mean,sigma,alpha,n);
}


// ------------------------------------------------------------------------
double SmearCrystalBall::crystalBallNorm(double mean, double sigma, double alpha, double n) const {
  double A = pow(n/alpha,n)*exp(-0.5*alpha*alpha);
  double B = n/alpha - alpha;
  double m = n - 1.;

  // Norm from Gaussian part
  double norm = sigma*sqrt(M_PI/2.)*( erf(alpha/sqrt(2)) - erf((mean-rMax_)/sqrt(2)/sigma) );
  // Norm from powerlaw part
  if( n == 1. ) {
    norm += A*sigma/pow(B,m)*log( (B+(mean-rMin_)/sigma)/(B+alpha) );
  } else {
    norm += A*sigma/m*(1./pow(B+alpha,m) - 1./pow(B+(mean-rMin_)/sigma,m));
  }

  return norm;
}



double SmearCrystalBall::crystalBallDerivative(double x, double mean, double sigma, double alpha, double n, int i) const {
  double d = 0.;
  
  double u = (x - mean)/sigma;
  double beta = (mean - rMin_)/sigma;
  double gamma = (mean - rMax_)/sigma;
  double A = pow(n/alpha,n)*exp(-0.5*alpha*alpha);
  double B = n/alpha - alpha;
  double m = n - 1.;
  double term = pow(B+alpha,-m) - pow(B+beta,-m);

  // df/dmean
  if( i == 0 ) {
    d = 0.;
  }
  // df/dsigma
  else if( i == 1 ) { 
    if( u > -alpha ) {
      d = u*u/sigma;
    } else {
      d = -n*u/sigma/(B-u);
    }
    d -= ( A*term/m - A*beta*sigma/pow(B+beta,n)
	   + sqrt(M_PI/2.)*(erf(alpha/sqrt(2.)) - erf(gamma/sqrt(2.)))
	   + gamma*exp(-gamma*gamma/2.) ) / crystalBallNorm(mean,sigma,alpha,n) ;
  }
  // df/alpha
  else if( i == 2 ) {
    if( u > -alpha ) {
      d = 0.;
    } else {
      d = (n*n/alpha/alpha + n)/(B - u) - (alpha + 1./alpha);
    }
    d -=( -A*term*sigma/m*(alpha + 1./alpha) 
	  - A*sigma/m*( (1.+1./alpha)*m/pow(B+alpha,n) - m/alpha/pow(B+beta,n) )
	  + exp(-alpha*alpha/2.) - sqrt(M_PI/2.)*sigma*erf(gamma/sqrt(2.))
	  ) / crystalBallNorm(mean,sigma,alpha,n);
  }
  // df/n
  else if( i == 3 ) {
    if( u > -alpha ) {
      d = 0.;
    } else {
      d = (1.+log(n/alpha)) - (log(B-u)+n/alpha/(B-u));
    }
    d -= ( A*sigma/m*(
		    1.+log(n/alpha)-term/m
		    -1./pow(B+alpha,m)*(log(B+alpha)+m/(B+alpha)*(1.+1./alpha))
		    +1./pow(B+beta,m)*(log(B+beta)+m/(B+beta)/alpha)
		    ) ) / crystalBallNorm(mean,sigma,alpha,n);
  }

  d *= crystalBallFunc(x,mean,sigma,alpha,n);

  return d;
}



//!  The probability density of \p pt if cuts
//!  \f$ \texttt{ptDijetMin\_} < p^{\textrm{dijet}}_{T} < \texttt{ptDijetMax\_} \f$
//!  are applied:
//!  \f[ t^{-n} \int^{\texttt{ptDijetMax\_}}_{\texttt{ptDijetMin\_}}
//!      dx\,\frac{r(x/t) \cdot t}{\sqrt{2}} \f]
// ------------------------------------------------------------------------
double  SmearCrystalBall::truthPDF(double pt, double n) const {
  double prob = 0.;
  if( tMin_ < pt && pt < tMax_ ) {
    prob = 1. / pow( pt, n );
    //int bin = ptDijetCutInt_->FindBin(pt);
    //prob *= ptDijetCutInt_->GetBinContent(bin);
    prob = 1.;
  }
  
  return prob;
}



// ------------------------------------------------------------------------
void SmearCrystalBall::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  Probability density of ptTrue spectrum:\n";
  std::cout << "    Powerlaw\n";
  std::cout << "    " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
  std::cout << "    " << ptDijetMin_ << " < ptDijet < " << ptDijetMax_ << " GeV\n";
  std::cout << "  Probability density of response:\n";
  std::cout << "    Non-zero for " << rMin_ << " < R < " << rMax_ << "\n";
  std::cout << std::endl;
}
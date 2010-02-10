//
//  $Id: Parametrization.cc,v 1.4 2010/02/09 10:19:23 mschrode Exp $
//
#include "Parametrization.h"


#include "TF1.h"
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
  : Parametrization(0,rNBins+2,0,1),
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
  ptDijetCutInt_ = new TH1D("norm","",100,tMin_,tMax_);
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
  for(int i = 0; i < static_cast<int>(nJetPars()); i++) { // Outer loop over parameters
    for(int j = 0; j < i+1; j++) { // Inner loop over parameters
      int idx = (i*i + i)/2 + j; // Index of (i,j) in covariance vector
      if( cov[covIdx[idx]] ) {
	if( i == j ) { // Diagonal terms
	  var += df[i]*df[i]*respParScales_[i]*respParScales_[i]*cov[covIdx[idx]];
	} else { // Off-diagonal terms
	  var += 2*df[i]*df[j]*respParScales_[i]*respParScales_[j]*cov[covIdx[idx]];
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
  // Norm of probability of dijet event configuration
//   double norm = 0.;
//   for(int bin = 1; bin < ptDijetCutInt_->GetNbinsX(); bin++) {
//     double pt = ptDijetCutInt_->GetBinCenter(bin);
//     //    norm += pow(pt,2-par[0]) * ptDijetCutInt_->GetBinContent(bin) * ptDijetCutInt_->GetXaxis()->GetBinWidth(1);
//     norm += pow(pt,2-par[0]) * ptDijetCutInt_->GetXaxis()->GetBinWidth(1);
//   }

  double norm = 0.;
  double expo = 3. - par[0];
  if( expo ) {
    norm = (pow(tMax_,expo) - pow(tMin_,expo))/expo;
  } else {
    norm = log(tMax_/tMin_);
  }
  
  double p = 0.;
  if( norm ) {
    p = truthPDF(x->pt,par[0]);
    p /= norm;
  }
  
  return p;

//   double norm = (pow(tMax_,3)-pow(tMin_,3))/3.;
//   return x->pt > tMin_ && x->pt < tMax_ ? 1./norm : 0.;
}

// ------------------------------------------------------------------------
double SmearStepGaussInter::correctedGlobalJetEtSigma(const Measurement *x,const double *par,const double *cov, const std::vector<int> &covIdx) const {
  double var = 0.;
  if( tMin_ < x->pt && x->pt < tMax_ ) {
    double n = par[0];
    double norm = 0.;
    double dnorm = 0.;
    if( n == 1. ) {
      norm = log(tMax_/tMin_);
      dnorm = 0.;
    } else {
      norm = (pow(tMax_,1.-n)-pow(tMin_,1.-n))/(1.-n);
      dnorm = norm/(1.-n) + (pow(tMax_,-n)-pow(tMin_,-n));
    }
    
    double df = -1./norm/pow(x->pt,n)*( n/x->pt + dnorm/norm );
    var = df*df*cov[covIdx[0]];
  }

  // Return standard deviation
  return sqrt(var);
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
  double prob = 0.;
  if( tMin_ < pt && pt < tMax_ ) {
    prob = 1. / pow( pt, n );
    //int bin = ptDijetCutInt_->FindBin(pt);
    //prob *= ptDijetCutInt_->GetBinContent(bin);
  }
  
  return prob;

//  return pt > tMin_ && pt < tMax_ ? 1./(tMax_-tMin_) : 0.;
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
  : Parametrization(0,4,0,1),
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
double SmearCrystalBall::resolution(double r, double pt, const double *rPar) const {
  return crystalBallFunc(r,mean(rPar),sigma(rPar),alpha(rPar),n(rPar));
}


double SmearCrystalBall::resolutionError(double r, double pt, const double *rPar, const double *cov, const std::vector<int> &covIdx) const {
  // Store derivatives
  std::vector<double> df(nJetPars());
  for(size_t i = 0; i < nJetPars(); i++) {
    df[i] = crystalBallDerivative(r,mean(rPar),sigma(rPar),alpha(rPar),n(rPar),i);
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
double SmearCrystalBall::spectrum(double pt, const double *sPar, const double *rPar) const {
  // Norm of probability of dijet event configuration
  double norm = 0.;

//   double expo = 3. - sPar[0];
//    if( expo ) {
//     norm = (pow(tMax_,expo) - pow(tMin_,expo))/expo;
//   } else {
//     norm = log(tMax_/tMin_);
//   }
  
  int nSteps = 100;
  double dt = (tMax_-tMin_)/nSteps;
  for(int i = 0; i < nSteps; i++) {
    double t = tMin_ + (i+0.5)*dt;
    norm += t*t*truthPDF(t,sPar,rPar);
  }
  norm *= dt;

  double p = 0.;
  if( norm ) {
    p = truthPDF(pt,sPar,rPar);
    p /= norm;
  }

  return p;
}


//!  The probability density of \p pt if cuts
//!  \f$ \texttt{ptDijetMin\_} < p^{\textrm{dijet}}_{T} < \texttt{ptDijetMax\_} \f$
//!  are applied:
//!  \f[ t^{-n} \int^{\texttt{ptDijetMax\_}}_{\texttt{ptDijetMin\_}}
//!      dx\,\frac{r(x/t) \cdot t}{\sqrt{2}} \f]
// ------------------------------------------------------------------------
double  SmearCrystalBall::truthPDF(double pt, const double *sPar, const double *rPar) const {
  double prob = 0.;
  if( tMin_ < pt && pt < tMax_ ) {
    prob = crystalBallInt(mean(rPar),sigma(rPar),alpha(rPar),n(rPar),ptDijetMin_/pt,ptDijetMax_/pt)
      /sqrt(2.)/pow(pt,sPar[0]);
  }
  
  return prob;
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
    f *= crystalBallNorm(mean,sigma,alpha,n);
  }

  return f;
}


// ------------------------------------------------------------------------
double SmearCrystalBall::crystalBallInt(double mean, double sigma, double alpha, double n, double min, double max) const {
  double A = pow(n/alpha,n)*exp(-0.5*alpha*alpha);
  double B = n/alpha - alpha;
  double m = n - 1.;
  double c = mean - alpha*sigma;

  double in = 0.;
  if( min > c ) {
    // Integral from Gaussian part
    in = sigma*sqrt(M_PI/2.)*( erf((mean-min)/sqrt(2)/sigma) - erf((mean-max)/sqrt(2)/sigma) );
  } else if( max < c ) {
    // Integral from powerlaw part
    if( n == 1. ) {
      in = A*sigma/pow(B,m)*log( (B+(mean-min)/sigma)/(B+(mean-max)/sigma) );
    } else {
      in = A*sigma/m*(1./pow(B+(mean-max)/sigma,m) - 1./pow(B+(mean-min)/sigma,m));
    }
  } else {
    // Integral from both parts
    in = sigma*sqrt(M_PI/2.)*( erf(alpha/sqrt(2)) - erf((mean-max)/sqrt(2)/sigma) );
    if( n == 1. ) {
      in += A*sigma/pow(B,m)*log( (B+(mean-min)/sigma)/(B+alpha) );
    } else {
      in += A*sigma/m*(1./pow(B+alpha,m) - 1./pow(B+(mean-min)/sigma,m));
    }
  }

  return in;
}


// ------------------------------------------------------------------------
double SmearCrystalBall::crystalBallNorm(double mean, double sigma, double alpha, double n) const {
  double norm = crystalBallInt(mean,sigma,alpha,n,rMin_,rMax_);
  return norm > 0 ? 1./norm : 0.;
}



double SmearCrystalBall::crystalBallDerivative(double x, double mean, double sigma, double alpha, double n, int i) const {
  double d = 0.;

  double u = (x - mean)/sigma;
  double beta = (rMax_ - mean)/sigma;
  double B = n/alpha - alpha;
  double m = n - 1.;
  double c = rMin_-1.;
  double s = n*sigma;
  double k = s - alpha*c - alpha*alpha*sigma;
  double ea = exp(-0.5*alpha*alpha);

  // df/dmu
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
    d -= (-beta*exp(-0.5*beta*beta)
      + ea/(m*alpha*pow(k,n))
      *((alpha*alpha-n)*pow(s,n)+n*n*(n-alpha*alpha)*sigma*pow(k,m)
	+n*pow(k,n)-n*k/sigma*pow(s,n)-(s*pow(k,n)-pow(s,n)*k)*n*(n-alpha*alpha)/k)
      +sqrt(M_PI/2.)*(erf(alpha/sqrt(2.))+erf(beta/sqrt(2.))))
      * crystalBallNorm(mean,sigma,alpha,n);
  }
  // df/alpha
  else if( i == 2 ) {
    if( u > -alpha ) {
      d = 0.;
    } else {
      d = -(alpha*alpha+n)*(u+alpha)/(alpha*alpha-n+alpha*u);
    }
    d -= (-ea*(n+alpha*alpha)*pow(s/k,n)/(m*alpha*alpha)
	  *(c*alpha+sigma*alpha*alpha-sigma*(1.-pow(k/s,n))))
      * crystalBallNorm(mean,sigma,alpha,n);
  }
  // df/n
  else if( i == 3 ) {
    if( u > -alpha ) {
      d = 0.;
    } else {
      d = (1.+log(n/alpha)) - (log(B-u)+n/alpha/(B-u));
    }
    d -= (ea/m/m/alpha/pow(k,n)
	  *pow(s,n)*(c*alpha*(n-2.)+sigma*(1.+alpha*alpha*(n-2.))
		     -sigma*pow(k/s,n)-m*k*log(s)+m*k*log(k)))
      * crystalBallNorm(mean,sigma,alpha,n);
  }
  d *= crystalBallFunc(x,mean,sigma,alpha,n);

  return d;
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




// ------------------------------------------------------------------------
SmearCrystalBallPt::SmearCrystalBallPt(double tMin, double tMax, double rMin, double rMax, double ptDijetMin, double ptDijetMax, const std::vector<double>& rParScales)
  : Parametrization(0,6,0,3),
    tMin_(tMin),
    tMax_(tMax),
    rMin_(rMin),
    rMax_(rMax),
    ptDijetMin_(ptDijetMin),
    ptDijetMax_(ptDijetMax),
    scale_(rParScales),
    warningPrinted_(false) {
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= rMin_ && rMin_ < rMax_ );
  assert( 0.0 <= ptDijetMin_ && ptDijetMin_ < ptDijetMax_ );
  assert( scale_.size() >= nJetPars() );
  
  print();
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
double SmearCrystalBallPt::correctedJetEt(const Measurement *x, const double *par) const {
  double rMean = mean(x,par);
  double rSigma = sigma(x,par);
  double rAlpha = alpha(x,par);
  double rN = n(x,par);

  if( rSigma <= 0. ) {
    if( !warningPrinted_ ) {
      std::cerr << "WARNING: sigma < 0\n";
      warningPrinted_ = true;
    }
    rSigma = 1E-5;
  }
  if( rAlpha <= 0. ) {
    if( !warningPrinted_ ) {
      std::cerr << "WARNING: alpha < 0\n";
      warningPrinted_ = true;
    }
    rAlpha = 1E-5;
  }
  if( rN <= 0. ) {
    if( !warningPrinted_ ) {
      std::cerr << "WARNING: n < 0\n";
      warningPrinted_ = true;
    }
    rN = 1. + 1E-5;
  }
    
  return crystalBallFunc(x->E,rMean,rSigma,rAlpha,rN);
}



//!  \param x   \p Measurement::pt is truth for which the 
//!             probability density is returned
//!  \param par Pointer to parameters (parameters are multiplied by the corresponding scale)
//!  \return Probability density of truth times normalization of dijet probability
// ------------------------------------------------------------------------
double SmearCrystalBallPt::correctedGlobalJetEt(const Measurement *x,const double *par) const {
  // Norm of probability of dijet event configuration
  double norm = 0.;

//   double expo = 3. - m(x,par);
//    if( expo ) {
//     norm = (pow(tMax_,expo) - pow(tMin_,expo))/expo;
//   } else {
//     norm = log(tMax_/tMin_);
//   }
  
//   TF1 *f = new TF1("f","x^(x[0]+[1]*x+[2]/x)",tMin_,tMax_);
//   for(unsigned int i = 0; i < nGlobalJetPars(); i++) {
//     f->SetParameter(i,par[i]);
//   }
//   norm = f->Integral(tMin_,tMax_);
//  delete f;

  int N = 200;
  double dPt = (tMax_-tMin_)/N;
  for(int i = 0; i < N; i++) {
    double pt = tMin_+(i+0.5)*dPt;
    norm += truthPDF(pt,m(pt,par));
  }
  double p = 0.;
  if( norm ) {
    p = truthPDF(x->pt,m(x->pt,par));
    p /= norm;
  }
  
  return p;
}



// ------------------------------------------------------------------------
double SmearCrystalBallPt::crystalBallFunc(double x, double mean, double sigma, double alpha, double n) const {
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

  return f*crystalBallNorm(mean,sigma,alpha,n);
}


// ------------------------------------------------------------------------
double SmearCrystalBallPt::crystalBallNorm(double mean, double sigma, double alpha, double n) const {
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

  return norm > 0 ? 1./norm : 0.;
}



double SmearCrystalBallPt::crystalBallDerivative(double x, double mean, double sigma, double alpha, double n, int i) const {
  double d = 0.;

  double u = (x - mean)/sigma;
  double beta = (rMax_ - mean)/sigma;
  double B = n/alpha - alpha;
  double m = n - 1.;
  double c = rMin_-1.;
  double s = n*sigma;
  double k = s - alpha*c - alpha*alpha*sigma;
  double ea = exp(-0.5*alpha*alpha);

  // df/dmu
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
    d -= (-beta*exp(-0.5*beta*beta)
      + ea/(m*alpha*pow(k,n))
      *((alpha*alpha-n)*pow(s,n)+n*n*(n-alpha*alpha)*sigma*pow(k,m)
	+n*pow(k,n)-n*k/sigma*pow(s,n)-(s*pow(k,n)-pow(s,n)*k)*n*(n-alpha*alpha)/k)
      +sqrt(M_PI/2.)*(erf(alpha/sqrt(2.))+erf(beta/sqrt(2.))))
      * crystalBallNorm(mean,sigma,alpha,n);
  }
  // df/alpha
  else if( i == 2 ) {
    if( u > -alpha ) {
      d = 0.;
    } else {
      d = -(alpha*alpha+n)*(u+alpha)/(alpha*alpha-n+alpha*u);
    }
    d -= (-ea*(n+alpha*alpha)*pow(s/k,n)/(m*alpha*alpha)
	  *(c*alpha+sigma*alpha*alpha-sigma*(1.-pow(k/s,n))))
      * crystalBallNorm(mean,sigma,alpha,n);
  }
  // df/n
  else if( i == 3 ) {
    if( u > -alpha ) {
      d = 0.;
    } else {
      d = (1.+log(n/alpha)) - (log(B-u)+n/alpha/(B-u));
    }
    d -= (ea/m/m/alpha/pow(k,n)
	  *pow(s,n)*(c*alpha*(n-2.)+sigma*(1.+alpha*alpha*(n-2.))
		     -sigma*pow(k/s,n)-m*k*log(s)+m*k*log(k)))
      * crystalBallNorm(mean,sigma,alpha,n);
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
double  SmearCrystalBallPt::truthPDF(double pt, double m) const {
  double prob = 0.;
  if( tMin_ < pt && pt < tMax_ ) {
    prob = 1. / pow( pt, m );
  }
  
  return prob;
}



// ------------------------------------------------------------------------
void SmearCrystalBallPt::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  Probability density of ptTrue spectrum:\n";
  std::cout << "    Powerlaw\n";
  std::cout << "    " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
  std::cout << "    " << ptDijetMin_ << " < ptDijet < " << ptDijetMax_ << " GeV\n";
  std::cout << "  Probability density of response:\n";
  std::cout << "    Non-zero for " << rMin_ << " < R < " << rMax_ << "\n";
  std::cout << std::endl;
}

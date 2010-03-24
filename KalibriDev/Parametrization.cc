//
//  $Id: Parametrization.cc,v 1.7 2010/02/25 15:28:19 mschrode Exp $
//
#include "Parametrization.h"


#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
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
  ptDijetCutInt_ = new TH1F("norm","",100,tMin_,tMax_);
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
SmearCrystalBall::SmearCrystalBall(double tMin, double tMax, double ptDijetMin, double ptDijetMax, const std::vector<double>& rParScales)
  : Parametrization(0,9,0,0),
    tMin_(tMin),
    tMax_(tMax),
    ptDijetMin_(ptDijetMin),
    ptDijetMax_(ptDijetMax),
    scale_(rParScales) {
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= ptDijetMin_ && ptDijetMin_ < ptDijetMax_ );
  assert( scale_.size() >= nJetPars() );

  hPdfPtTrue_ = new TH1F("hPdfPtTrue_","",1000,tMin,tMax);
  //Fill values of unnormalised truth pdf
  double *parSpec = new double[9];
  parSpec[0] = 0.;
  parSpec[1] = 0.;
  parSpec[2] = 0.;
  parSpec[3] = 5.01;
  parSpec[4] = 2.77;
  parSpec[5] = 8.63;
  parSpec[6] = 1.37;
  parSpec[7] = 1.188;
  parSpec[8] = 7.59;
  // For pt-dependent sigma
  double *parSigma = new double[3];
  parSigma[0] = 0.;
  parSigma[1] = 1.19;
  parSigma[2] = 0.03;
  double alpha = 2.;
  double n = 4.411;
  for(int bin = 1; bin <= hPdfPtTrue_->GetNbinsX(); bin++) {
    double ptTrue = hPdfPtTrue_->GetBinCenter(bin);
    hPdfPtTrue_->SetBinContent(bin,pdfPtTrueNotNorm(ptTrue,parSpec,parSigma,alpha,n));
  }
  // Normalise values of truth pdf
  hPdfPtTrue_->Scale(1./hPdfPtTrue_->Integral("width"));
  delete [] parSpec;
  delete [] parSigma;
  
  print();
}

SmearCrystalBall::~SmearCrystalBall() {
  delete hPdfPtTrue_;
}

// ------------------------------------------------------------------------
double SmearCrystalBall::pdfPtMeas(double ptMeas, double ptTrue, const double *par) const {
  return crystalBallFunc(ptMeas,ptTrue,sigma(par),alpha(par),n(par));
}


// ------------------------------------------------------------------------
double SmearCrystalBall::pdfPtTrue(double ptTrue, const double *par) const {
//   double pdf = 0.;
//   if( tMin_ < ptTrue && ptTrue < tMax_ ) {
// //     double m = 1.-exponent(par);
// //     double norm = ( m == 0. ? log(tMax_/tMin_) : (pow(tMax_,m)-pow(tMin_,m))/m );
// //     pdf = 1./pow(ptTrue,exponent(par))/norm;

//     double norm = exp(-specPar(par,0))*(exp(-specPar(par,1)*tMin_)-exp(-specPar(par,1)*tMax_))/specPar(par,1);
//     norm += exp(-specPar(par,2))*(exp(-specPar(par,3)*tMin_)-exp(-specPar(par,3)*tMax_))/specPar(par,3);
//     norm += exp(-specPar(par,4))*(exp(-specPar(par,5)*tMin_)-exp(-specPar(par,5)*tMax_))/specPar(par,5);

//     pdf = exp(-specPar(par,0)-specPar(par,1)*ptTrue);
//     pdf += exp(-specPar(par,2)-specPar(par,3)*ptTrue);
//     pdf += exp(-specPar(par,4)-specPar(par,5)*ptTrue);
//     pdf /= norm;
//   }
//   return pdf;

  return hPdfPtTrue_->GetBinContent(hPdfPtTrue_->FindBin(ptTrue));
}


// ------------------------------------------------------------------------
double SmearCrystalBall::pdfResponse(double r, double ptTrue, const double *par) const {
  return crystalBallFunc(r,1.,sigma(par)/ptTrue,alpha(par),n(par));
}


// ------------------------------------------------------------------------
double SmearCrystalBall::pdfResponseError(double r, double ptTrue, const double *par, const double *cov, const std::vector<int> &covIdx) const {
  // Store derivatives
  std::vector<double> df(nJetPars());
  for(size_t i = 0; i < nJetPars(); i++) {
    df[i] = crystalBallDerivative(r,1.,sigma(par)/ptTrue,alpha(par),n(par),i);
  }

  // Calculate variance
  double var = 0.;
  for(int i = 0; i < static_cast<int>(nJetPars()); i++) { // Outer loop over parameters
    for(int j = 0; j < i+1; j++) { // Inner loop over parameters
      int idx = (i*i + i)/2 + j; // Index of (i,j) in covariance vector
      if( cov[idx] ) {
	if( i == j ) { // Diagonal terms
	  var += df[i]*df[i]*scale_[i]*scale_[i]*cov[idx];
	} else { // Off-diagonal terms
	  var += 2*df[i]*df[j]*scale_[i]*scale_[j]*cov[idx];
	}
      }
    } // End of inner loop over parameters
  } // End of outer loop over parameters
  // Return standard deviation
  return sqrt(var);
}


// ------------------------------------------------------------------------
double SmearCrystalBall::pdfDijetAsym(double a, double ptTrue, const double *par) const {
  double s = sigma(par)/ptTrue/sqrt(2.);
  double u = a/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


//! The Crystal Ball function value normalized to
//! the integral from 0 to infinity
// ------------------------------------------------------------------------
double SmearCrystalBall::crystalBallFunc(double x, double mean, double sigma, double alpha, double n) const {
  double f = 0.;
  if( x > 0. ) {
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


//! The inverse of the integral over Crystal Ball function 
//! from 0 to infinity
// ------------------------------------------------------------------------
double SmearCrystalBall::crystalBallNorm(double mean, double sigma, double alpha, double n) const {
  double norm = sigma*exp(-0.5*alpha*alpha)*n/alpha/(n-1.);
  norm *= (1. - pow(1. - alpha*alpha/n + alpha*mean/n/sigma, 1.-n));
  norm += sqrt(M_PI/2.)*sigma*(1. + erf(alpha/sqrt(2.)));

  return norm == 0. ? 0. : 1./norm;
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


double SmearCrystalBall::crystalBallDerivative(double x, double mean, double sigma, double alpha, double n, int i) const {
  double d = 0.;

  // TODO: rewrite derivative for rMin = 0 and rMax --> infty
  double rMin = 0.;
  double rMax = 2.;

  double u = (x - mean)/sigma;
  double beta = (rMax - mean)/sigma;
  double B = n/alpha - alpha;
  double m = n - 1.;
  double c = rMin-1.;
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


double  SmearCrystalBall::pdfPtTrueNotNorm(double ptTrue, const double *parSpec, const double *parSigma, double alpha, double n) const {
  double pdf = 0.;

  // Underlying truth pdf, normalised from 0 to infinity
  //  pdf = exp(-ptTrue/par[3])/par[3];

  double norm = exp(-specPar(parSpec,0))*(exp(-specPar(parSpec,1)*tMin_)-exp(-specPar(parSpec,1)*tMax_))/specPar(parSpec,1);
  norm += exp(-specPar(parSpec,2))*(exp(-specPar(parSpec,3)*tMin_)-exp(-specPar(parSpec,3)*tMax_))/specPar(parSpec,3);
  norm += exp(-specPar(parSpec,4))*(exp(-specPar(parSpec,5)*tMin_)-exp(-specPar(parSpec,5)*tMax_))/specPar(parSpec,5);

  pdf = exp(-specPar(parSpec,0)-specPar(parSpec,1)*ptTrue);
  pdf += exp(-specPar(parSpec,2)-specPar(parSpec,3)*ptTrue);
  pdf += exp(-specPar(parSpec,4)-specPar(parSpec,5)*ptTrue);
  pdf /= norm;

//   // Convolution with cuts on ptdijet assuming Gaussian response
//   double s = sqrt( par[9]*par[9] + par[10]*par[10]*ptTrue + par[11]*par[11]*ptTrue*ptTrue );
//   double min = -ptTrue/sqrt(2.)/s;
//   double max = (2*ptDijetMax_ - ptTrue)/sqrt(2.)/s;
//   int nSteps = hPdfPtTrue_->GetNbinsX();
//   double dz = (max-min)/nSteps;
//   double c = 0.;
//   for(int i = 0; i < nSteps; i++) {
//     double z = min + (0.5+i)*dz;
//     c += exp(-z*z)*( erf(sqrt(2.)*(ptDijetMax_-ptTrue)/s-z) - erf(sqrt(2.)*(ptDijetMin_-ptTrue)/s-z) );
//   }
//   c *= dz/2./sqrt(M_PI) *4./pow((1.+erf(ptTrue/sqrt(2.)/s)),2.);

  // Convolution with cuts on ptdijet assuming Crystal Ball response
  double s = sqrt( parSigma[0]*parSigma[0] + parSigma[1]*parSigma[1]*ptTrue + parSigma[2]*parSigma[2]*ptTrue*ptTrue );
  int nSteps = hPdfPtTrue_->GetNbinsX();
  double dz = 2.*ptDijetMax_/nSteps;
  double c = 0.;
  for(int i = 0; i < nSteps; i++) {
    double z = (0.5+i)*dz;
    c += crystalBallFunc(z,ptTrue,s,alpha,n)*crystalBallInt(ptTrue,s,alpha,n,2.*ptDijetMin_-z,2.*ptDijetMax_-z);
  }
  c *= dz;

//   // Convolution with cuts on pt of 1. jet assuming Crystal Ball response
//   double s = sqrt( parSigma[0]*parSigma[0] + parSigma[1]*parSigma[1]*ptTrue + parSigma[2]*parSigma[2]*ptTrue*ptTrue );
//   int nSteps = hPdfPtTrue_->GetNbinsX();
//   double dz = (ptDijetMax_-ptDijetMin_)/nSteps;
//   double c = 0.;
//   for(int i = 0; i < nSteps; i++) {
//     double z = (0.5+i)*dz;
//     c += crystalBallFunc(z,ptTrue,s,alpha,n)*crystalBallInt(ptTrue,s,alpha,n,ptDijetMin_,ptDijetMax_);
//   }
//   c *= dz;
  
  return c*pdf;
}



// ------------------------------------------------------------------------
void SmearCrystalBall::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  Probability density of ptTrue spectrum:\n";
  std::cout << "    Powerlaw\n";
  std::cout << "    " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
  std::cout << "    " << ptDijetMin_ << " < ptDijet < " << ptDijetMax_ << " GeV\n";
  std::cout << std::endl;
}




// ------------------------------------------------------------------------
SmearCrystalBallPt::SmearCrystalBallPt(double tMin, double tMax, double ptDijetMin, double ptDijetMax, const std::vector<double>& rParScales)
  : Parametrization(0,10,0,0),
    tMin_(tMin),
    tMax_(tMax),
    ptDijetMin_(ptDijetMin),
    ptDijetMax_(ptDijetMax),
    scale_(rParScales) {
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= ptDijetMin_ && ptDijetMin_ < ptDijetMax_ );
  assert( scale_.size() >= nJetPars() );
  
  print();
}


// ------------------------------------------------------------------------
double SmearCrystalBallPt::pdfPtMeas(double ptMeas, double ptTrue, const double *par) const {
  return crystalBallFunc(ptMeas,ptTrue,sigma(ptTrue,par),alpha(ptTrue,par),n(par));
}


// ------------------------------------------------------------------------
double SmearCrystalBallPt::pdfPtTrue(double ptTrue, const double *par) const {
  double pdf = 0.;
  if( tMin_ < ptTrue && ptTrue < tMax_ ) {
//     // Powerlaw spectrum
//     double m = 1.-exponent(par);
//     double norm = ( m == 0. ? log(tMax_/tMin_) : (pow(tMax_,m)-pow(tMin_,m))/m );
//     pdf = 1./pow(ptTrue,exponent(par))/norm;

    double norm = exp(-specPar(par,0))*(exp(-specPar(par,1)*tMin_)-exp(-specPar(par,1)*tMax_))/specPar(par,1);
    norm += exp(-specPar(par,2))*(exp(-specPar(par,3)*tMin_)-exp(-specPar(par,3)*tMax_))/specPar(par,3);
    pdf = exp(-specPar(par,0)-specPar(par,1)*ptTrue);
    pdf += exp(-specPar(par,2)-specPar(par,3)*ptTrue);
    pdf /= norm;
  }
  return pdf;
}


// ------------------------------------------------------------------------
double SmearCrystalBallPt::pdfResponse(double r, double ptTrue, const double *par) const {
  return crystalBallFunc(r,1.,sigma(ptTrue,par)/ptTrue,alpha(ptTrue,par),n(par));
}


//! The Crystal Ball function value normalized to
//! the integral from 0 to infinity
// ------------------------------------------------------------------------
double SmearCrystalBallPt::crystalBallFunc(double x, double mean, double sigma, double alpha, double n) const {
  double f = 0.;
  if( x > 0. ) {
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


//! The inverse of the integral over Crystal Ball function 
//! from 0 to infinity
// ------------------------------------------------------------------------
double SmearCrystalBallPt::crystalBallNorm(double mean, double sigma, double alpha, double n) const {
  double norm = sigma*exp(-0.5*alpha*alpha)*n/alpha/(n-1.);
  norm *= (1. - pow(1. - alpha*alpha/n + alpha*mean/n/sigma, 1.-n));
  norm += sqrt(M_PI/2.)*sigma*(1. + erf(alpha/sqrt(2.)));

  return norm == 0. ? 0. : 1./norm;
}


// ------------------------------------------------------------------------
void SmearCrystalBallPt::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  Probability density of ptTrue spectrum:\n";
  std::cout << "    Powerlaw\n";
  std::cout << "    " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
  std::cout << "    " << ptDijetMin_ << " < ptDijet < " << ptDijetMax_ << " GeV\n";
  std::cout << std::endl;
}



// ------------------------------------------------------------------------
SmearGauss::SmearGauss(double tMin, double tMax, double ptDijetMin, double ptDijetMax, const std::vector<double>& rParScales)
  : Parametrization(0,4,0,0),
    tMin_(tMin),
    tMax_(tMax),
    ptDijetMin_(ptDijetMin),
    ptDijetMax_(ptDijetMax),
    scale_(rParScales),
    hPdfPtTrue_(0) {
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= ptDijetMin_ && ptDijetMin_ < ptDijetMax_ );
  assert( scale_.size() >= nJetPars() );

  hPdfPtTrue_ = new TH1F("hPdfPtTrue_","",10000,tMin,tMax);
  //Fill values of unnormalised truth pdf
  double *truePar = new double[nJetPars()];
  // Sigma parameters
  truePar[0] = 4.;
  truePar[1] = 1.2;
  truePar[2] = 0.05;
  // Spectrum parameters
  truePar[3] = 80.;
  for(int bin = 1; bin <= hPdfPtTrue_->GetNbinsX(); bin++) {
    double ptTrue = hPdfPtTrue_->GetBinCenter(bin);
    hPdfPtTrue_->SetBinContent(bin,pdfPtTrueNotNorm(ptTrue,truePar));
  }
  // Normalise values of truth pdf
  hPdfPtTrue_->Scale(1./hPdfPtTrue_->Integral("width"));
  delete [] truePar;

// //   TFile file("~/UserCode/mschrode/resolutionFit/jsResponse.root","READ");
// //   file.GetObject("hPtGen",hPdfPtTrue_);
// //   if( !hPdfPtTrue_ ) {
// //     std::cerr << "ERROR: No histogram found in file '" << file.GetName() << "'\n";
// //     exit(1);
// //   } else {
// //     hPdfPtTrue_->SetDirectory(0);
// //     hPdfPtTrue_->SetName("hPdfPtTrue");
// //     int binMin = hPdfPtTrue_->FindBin(tMin_);
// //     int binMax = hPdfPtTrue_->FindBin(tMax_);
// //     std::cout << "Integal: " << hPdfPtTrue_->Integral(binMin,binMax,"width") << std::endl;
// //   }
// //   file.Close();

  print();
}


SmearGauss::~SmearGauss() { if( hPdfPtTrue_ ) delete hPdfPtTrue_; }


// ------------------------------------------------------------------------
double SmearGauss::pdfPtMeas(double ptMeas, double ptTrue, const double *par) const {
  double s = sigma(ptTrue,par);
  double u = (ptMeas - ptTrue)/s;
  double norm = sqrt(M_PI/2.)*s*sqrt(( erf((ptDijetMax_-ptTrue)/sqrt(2.)/s) - erf((ptDijetMin_-ptTrue)/sqrt(2.)/s) ));
  // This should be caught more cleverly
  if( norm < 1E-10 ) norm = 1E-10;

  return exp(-0.5*u*u)/norm;
}


// ------------------------------------------------------------------------
double  SmearGauss::pdfPtTrue(double ptTrue, const double *par) const {
//   double pdf = 0.;
//   if( tMin_ < ptTrue && ptTrue < tMax_ ) {
// //     double m = 1.-exponent(par);
// //     double norm = ( m == 0. ? log(tMax_/tMin_) : (pow(tMax_,m)-pow(tMin_,m))/m );
// //     pdf = 1./pow(ptTrue,exponent(par))/norm;

//     // ToyMC parametrization of spectrum
//     double tau = exponent(par);
//     double norm = tau*(exp(-tMin_/tau)-exp(-tMax_/tau));
//     pdf = exp(-ptTrue/tau)/norm;

// //     double norm = exp(-specPar(par,0))*(exp(-specPar(par,1)*tMin_)-exp(-specPar(par,1)*tMax_))/specPar(par,1);
// //     norm += exp(-specPar(par,2))*(exp(-specPar(par,3)*tMin_)-exp(-specPar(par,3)*tMax_))/specPar(par,3);
// //     norm += exp(-specPar(par,4))*(exp(-specPar(par,5)*tMin_)-exp(-specPar(par,5)*tMax_))/specPar(par,5);

// //     pdf = exp(-specPar(par,0)-specPar(par,1)*ptTrue);
// //     pdf += exp(-specPar(par,2)-specPar(par,3)*ptTrue);
// //     pdf += exp(-specPar(par,4)-specPar(par,5)*ptTrue);
// //     pdf /= norm;

// //    pdf = 1./(tMax_-tMin_);
//   }
//   return pdf;
  
  return hPdfPtTrue_->GetBinContent(hPdfPtTrue_->FindBin(ptTrue));
}


// ------------------------------------------------------------------------
double SmearGauss::pdfResponse(double r, double ptTrue, const double *par) const {
  double s = sigma(ptTrue,par)/ptTrue;
  double u = (r - 1.)/s;
  double cut = (1.+erf(ptTrue/sqrt(2.)/s))/2.;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s/cut;
}


// ------------------------------------------------------------------------
double SmearGauss::pdfDijetAsym(double a, double ptTrue, const double *par) const {
  double s = sigma(ptTrue,par)/ptTrue/sqrt(2.);
  double u = a/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


// ------------------------------------------------------------------------
double SmearGauss::pdfResponseError(double r, double ptTrue, const double *par, const double *cov, const std::vector<int> &covIdx) const {
  // Store derivatives
  std::vector<double> df(nJetPars());
  for(size_t i = 0; i < nJetPars(); i++) {
    df[i] = pdfResponseDeriv(r,ptTrue,par,i);
  }

  // Calculate variance
  double var = 0.;
  for(int i = 0; i < static_cast<int>(nJetPars()); i++) { // Outer loop over parameters
    for(int j = 0; j < i+1; j++) { // Inner loop over parameters
      int idx = (i*i + i)/2 + j; // Index of (i,j) in covariance vector
      if( cov[idx] ) {
	if( i == j ) { // Diagonal terms
	  var += df[i]*df[i]*scale_[i]*scale_[i]*cov[idx];
	} else { // Off-diagonal terms
	  var += 2*df[i]*df[j]*scale_[i]*scale_[j]*cov[idx];
	}
      }
    } // End of inner loop over parameters
  } // End of outer loop over parameters
  // Return standard deviation
  return sqrt(var);
}


// ------------------------------------------------------------------------
double SmearGauss::pdfResponseDeriv(double r, double ptTrue, const double *par, int i) const {
  double df = 0.;
  if( i < 3 ) {
    double s = sigma(ptTrue,par);
    double u = ptTrue*(r-1.)/s;
    df = pdfResponse(r,ptTrue,par) * scale_[i]*par[i]/s/s * (u*u - 1.);
    if( i == 1 ) df *= ptTrue;
    if( i == 2 ) df *= ptTrue*ptTrue;
  }

  return df;
}


// ------------------------------------------------------------------------
double  SmearGauss::pdfPtTrueNotNorm(double ptTrue, const double *par) const {
  double pdf = 0.;

  // Underlying truth pdf

  // ToyMC parametrization of spectrum, normalized from 0 to infty
  pdf = exp(-ptTrue/par[3])/par[3];

//   double norm = exp(-specPar(par,0))*(exp(-specPar(par,1)*tMin_)-exp(-specPar(par,1)*tMax_))/specPar(par,1);
//   norm += exp(-specPar(par,2))*(exp(-specPar(par,3)*tMin_)-exp(-specPar(par,3)*tMax_))/specPar(par,3);
//   norm += exp(-specPar(par,4))*(exp(-specPar(par,5)*tMin_)-exp(-specPar(par,5)*tMax_))/specPar(par,5);

//   pdf = exp(-specPar(par,0)-specPar(par,1)*ptTrue);
//   pdf += exp(-specPar(par,2)-specPar(par,3)*ptTrue);
//   pdf += exp(-specPar(par,4)-specPar(par,5)*ptTrue);
//   pdf /= norm;

//   // Convolution with cuts on ptdijet
//   double s = sigma(ptTrue,par)/sqrt(2.);
//   double c = 0.5*( erf((ptDijetMax_-ptTrue)/s/sqrt(2.)) - erf((ptDijetMin_-ptTrue)/s/sqrt(2.)) );

  // Convolution with cuts on 1. jet (randomly assigned)
  double s = sqrt( par[0]*par[0] + par[1]*par[1]*ptTrue + par[2]*par[2]*ptTrue*ptTrue );
  double c = 0.5*( erf((ptDijetMax_-ptTrue)/s/sqrt(2.)) - erf((ptDijetMin_-ptTrue)/s/sqrt(2.)) );

  return c*pdf;
}



// ------------------------------------------------------------------------
void SmearGauss::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  Probability density of ptTrue spectrum:\n";
  std::cout << "    Powerlaw\n";
  std::cout << "    " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
  std::cout << "    " << ptDijetMin_ << " < ptDijet < " << ptDijetMax_ << " GeV\n";
  std::cout << std::endl;
}



// ------------------------------------------------------------------------
SmearGaussPtBin::SmearGaussPtBin(double tMin, double tMax, double ptDijetMin, double ptDijetMax, const std::vector<double>& rParScales)
  : Parametrization(0,7,0,0),
    tMin_(tMin),
    tMax_(tMax),
    ptDijetMin_(ptDijetMin),
    ptDijetMax_(ptDijetMax),
    scale_(rParScales) {
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= ptDijetMin_ && ptDijetMin_ < ptDijetMax_ );
  assert( scale_.size() >= nJetPars() );

  TFile file("input/jsResponse_Unweighted.root","READ");
  file.GetObject("hPtGen",hPurePdfPtTrue_);
  if( !hPurePdfPtTrue_ ) {
    std::cerr << "ERROR: No histogram found in file '" << file.GetName() << "'\n";
    exit(1);
  } else {
    hPurePdfPtTrue_->SetDirectory(0);
    hPurePdfPtTrue_->SetName("hPurePdfPtTrue_");
    int binMin = hPurePdfPtTrue_->FindBin(tMin_);
    int binMax = hPurePdfPtTrue_->FindBin(tMax_);
    if( hPurePdfPtTrue_->Integral(binMin,binMax,"width") ) {
      hPurePdfPtTrue_->Scale(1./(hPurePdfPtTrue_->Integral(binMin,binMax,"width")));
    }
  }
  file.Close();

  hPdfPtTrue_ = new TH1F("hPdfPtTrue_","",1000,tMin,tMax);
  //Fill values of unnormalised truth pdf
  double *truePar = new double[nJetPars()+3];
  truePar[0] = 0.;
  truePar[1] = 5.7;
  truePar[2] = 3.00;
  truePar[3] = 3.91;
  truePar[4] = 1.52;
  truePar[5] = 7.15;
  truePar[6] = 8.37;
  // For pt-dependent sigma
  truePar[7] = 0.;
  truePar[8] = 1.4885;//1.145;
  truePar[9] = 0.0481;//0.037;
  for(int bin = 1; bin <= hPdfPtTrue_->GetNbinsX(); bin++) {
    double ptTrue = hPdfPtTrue_->GetBinCenter(bin);
    hPdfPtTrue_->SetBinContent(bin,pdfPtTrueNotNorm(ptTrue,truePar));
  }
  // Normalise values of truth pdf
  hPdfPtTrue_->Scale(1./hPdfPtTrue_->Integral("width"));
  delete [] truePar;
  
  print();
}


SmearGaussPtBin::~SmearGaussPtBin() { 
  if( hPdfPtTrue_ ) delete hPdfPtTrue_;
  if( hPurePdfPtTrue_ ) delete hPurePdfPtTrue_;
 }


// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfPtMeas(double ptMeas, double ptTrue, const double *par) const {
  double s = sigma(par);
  double u = (ptMeas - ptTrue)/s;
  double norm = sqrt(M_PI/2.)*s*sqrt(( erf((ptDijetMax_-ptTrue)/sqrt(2.)/s) - erf((ptDijetMin_-ptTrue)/sqrt(2.)/s) ));
  // This should be caught more cleverly
  if( norm < 1E-10 ) norm = 1E-10;

  return 0.;//exp(-0.5*u*u)/norm;
}

double SmearGaussPtBin::pdfPtMeasJet1(double ptMeas, double ptTrue, const double *par) const {
  double pdf = 0.;
  if( ptDijetMin_ < ptMeas && ptMeas < ptDijetMax_ ) {
    double s = sigma(par);
    double u = (ptMeas - ptTrue)/s;
    double norm = sqrt(M_PI/2.)*s*( erf((ptDijetMax_-ptTrue)/sqrt(2.)/s) - erf((ptDijetMin_-ptTrue)/sqrt(2.)/s) );
    // This should be caught more cleverly
    if( norm < 1E-10 ) norm = 1E-10;
    
    pdf = exp(-0.5*u*u)/norm; 
  }

  return pdf;
}

double SmearGaussPtBin::pdfPtMeasJet2(double ptMeas, double ptTrue, const double *par) const {
  double s = sigma(par);
  double u = (ptMeas - ptTrue)/s;
  double norm = sqrt(M_PI/2.)*s*( 1. + erf(ptTrue/sqrt(2.)/s) );
  // This should be caught more cleverly
  if( norm < 1E-10 ) norm = 1E-10;

  return exp(-0.5*u*u)/norm;
}



// ------------------------------------------------------------------------
double  SmearGaussPtBin::pdfPtTrue(double ptTrue, const double *par) const {
//   double pdf = 0.;
//   if( tMin_ < ptTrue && ptTrue < tMax_ ) {
// //     double m = 1.-exponent(par);
// //     double norm = ( m == 0. ? log(tMax_/tMin_) : (pow(tMax_,m)-pow(tMin_,m))/m );
// //     pdf = 1./pow(ptTrue,exponent(par))/norm;

// //     double tau = par[1];
// //     double norm = tau*(exp(-tMin_/tau)-exp(-tMax_/tau));
// //     pdf = exp(-ptTrue/tau)/norm;
    
//     double norm = exp(-specPar(par,0))*(exp(-specPar(par,1)*tMin_)-exp(-specPar(par,1)*tMax_))/specPar(par,1);
//     norm += exp(-specPar(par,2))*(exp(-specPar(par,3)*tMin_)-exp(-specPar(par,3)*tMax_))/specPar(par,3);
//     norm += exp(-specPar(par,4))*(exp(-specPar(par,5)*tMin_)-exp(-specPar(par,5)*tMax_))/specPar(par,5);

//     pdf = exp(-specPar(par,0)-specPar(par,1)*ptTrue);
//     pdf += exp(-specPar(par,2)-specPar(par,3)*ptTrue);
//     pdf += exp(-specPar(par,4)-specPar(par,5)*ptTrue);
//     pdf /= norm;

// //    pdf = 1./(tMax_-tMin_);
//   }
//   return pdf;
  
  return hPdfPtTrue_->GetBinContent(hPdfPtTrue_->FindBin(ptTrue));
}


// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfResponse(double r, double ptTrue, const double *par) const {
  double s = sigma(par)/ptTrue;
  double u = (r - 1.)/s;
  //double cut = (1.+erf(ptTrue/sqrt(2.)/s))/2.;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfResponseError(double r, double ptTrue, const double *par, const double *cov, const std::vector<int> &covIdx) const {
  // Store derivatives
  double s = sigma(par)/ptTrue;
  double u = (r-1.)/s;
  double df = pdfResponse(r,ptTrue,par)*(u*u - 1.)/s;

  // Calculate variance
  double var = df*df*scale_[0]*scale_[0]*cov[0]/ptTrue/ptTrue;

  // Return standard deviation
  return sqrt(var);
}


// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfDijetAsym(double a, double ptTrue, const double *par) const {
  double s = sigma(par)/ptTrue/sqrt(2.);
  double u = a/s;
  //double cut = (1.+erf(ptTrue/sqrt(2.)/s))/2.;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


// ------------------------------------------------------------------------
double  SmearGaussPtBin::pdfPtTrueNotNorm(double ptTrue, const double *par) const {
  double pdf = 0.;

  // Underlying truth pdf
  //  pdf = exp(-ptTrue/par[1])/par[1];

  double norm = exp(-specPar(par,0))*(exp(-specPar(par,1)*tMin_)-exp(-specPar(par,1)*tMax_))/specPar(par,1);
  norm += exp(-specPar(par,2))*(exp(-specPar(par,3)*tMin_)-exp(-specPar(par,3)*tMax_))/specPar(par,3);
  norm += exp(-specPar(par,4))*(exp(-specPar(par,5)*tMin_)-exp(-specPar(par,5)*tMax_))/specPar(par,5);

  pdf = exp(-specPar(par,0)-specPar(par,1)*ptTrue);
  pdf += exp(-specPar(par,2)-specPar(par,3)*ptTrue);
  pdf += exp(-specPar(par,4)-specPar(par,5)*ptTrue);
  pdf /= norm;

//  pdf = hPurePdfPtTrue_->GetBinContent(hPurePdfPtTrue_->FindBin(ptTrue));


//   // Convolution with cuts on ptdijet
//   double s = sqrt( par[7]*par[7] + par[8]*par[8]*ptTrue + par[9]*par[9]*ptTrue*ptTrue )/sqrt(2.);
//   double c = 0.5*( erf((ptDijetMax_-ptTrue)/s/sqrt(2.)) - erf((ptDijetMin_-ptTrue)/s/sqrt(2.)) );

  // Convolution with cuts on 1. jet
  double s = sqrt( par[7]*par[7] + par[8]*par[8]*ptTrue + par[9]*par[9]*ptTrue*ptTrue );
  double c = 0.5*( erf((ptDijetMax_-ptTrue)/s/sqrt(2.)) - erf((ptDijetMin_-ptTrue)/s/sqrt(2.)) );
 
  return c*pdf;
}


// ------------------------------------------------------------------------
void SmearGaussPtBin::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  Probability density of ptTrue spectrum:\n";
  std::cout << "    Powerlaw\n";
  std::cout << "    " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
  std::cout << "    " << ptDijetMin_ << " < ptDijet < " << ptDijetMax_ << " GeV\n";
  std::cout << std::endl;
}

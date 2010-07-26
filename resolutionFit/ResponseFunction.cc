#include "ResponseFunction.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "TRandom3.h"


namespace resolutionFit {
  ResponseFunction::ResponseFunction(Type type)
    : nBinsIntegration_(5000), zMinIntegration_(-1.), zMaxIntegration_(1.),
      deltaZIntegration_((zMaxIntegration_-zMinIntegration_)/nBinsIntegration_),
      epsIntegration_(1E-4), maxNIterIntegration_(20), type_(type) {
    rand_ = new TRandom3(0);
    fCB_ = new CrystalBallFunction();
  }


  
  ResponseFunction::~ResponseFunction() {
    delete rand_;
    delete fCB_;
  }



  double ResponseFunction::randomGauss(const std::vector<double> &par) const {
    return rand_->Gaus(par[0],par[1]);
  }


  int ResponseFunction::nPars() const {
    int n = 0;
    if( type() == Gauss ) n = 2;
    else if( type() == CrystalBall ) n = 4;
    else if( type() == TruncCrystalBall ) n = 5;

    return n;
  }


  TString ResponseFunction::typeLabel() const {
    TString label = "";
    if( type() == Gauss ) label = "Gauss";
    else if( type() == CrystalBall ) label = "Crystal Ball";
    else if( type() == TruncCrystalBall ) label = "Truncated Crystal Ball";
    return label;   
  }


  double ResponseFunction::pdfGauss(double pt, const std::vector<double> &par) const {
    double u = (pt-par[0])/par[1];
    double f =  exp(-0.5*u*u)/sqrt(2.*M_PI)/par[1];
    return f;
  }


  double ResponseFunction::sigmaGauss(double pt, const std::vector<double> &par) const {
    return sqrt( par[0]*par[0]/pt/pt + par[1]*par[1]/pt + par[2]*par[2]);
  }


  double ResponseFunction::pdf(double pt, Type type, const std::vector<double> &par) const {
    double f = 0.;
    if( type == Gauss ) {
      assert( par.size() == 2 );
      f = pdfGauss(pt,par);
    } else if( type == CrystalBall ) {
      assert( par.size() == 4 );
      f = pdfCrystalBall(pt,par);
    } else if( type == TruncCrystalBall ) {
      assert( par.size() == 5 );
      f = pdfTruncCrystalBall(pt,par);
    }

    return f;
  }


  double ResponseFunction::pdfAsymmetry(double a, Type type, const std::vector<double> &par) const {
    double h = zMaxIntegration_-zMinIntegration_;
    double asym = 0.;   
    double asymOld = 1.;
    double eps = 1.;
    int nIter = 0;
    std::vector<double> sAsym;
    std::vector<double> sAsymOld;

    if( a == 0 ) {
      std::cerr << "WARNING: A = 0" << std::endl;
    } else {
      // Iterate until precision or max. number iterations reached
      while( eps > epsIntegration_ && nIter < maxNIterIntegration_ ) {
	// Init iteration
	asymOld = asym;
	asym = 0.;
	sAsymOld = sAsym;
	sAsym.clear();
	
	// In each iteration, split h into 3 new intervals
	h /= 3.;        
	
	// Loop over nodes xi i.e. interval borders
	for(int i = 0; i <= pow(3.0,nIter+1); ++i){
	  double z = zMinIntegration_ + i*h;
	  
	  // Calculate pdf only at new nodes
	  if( nIter == 0 || i % 3 != 0 ) {
	    if( z == 0 ) {
	      std::cerr << "WARNING: z = 0" << std::endl;
	      sAsym.push_back(0.);
	    } else {
	      sAsym.push_back(pdf((z+z/a)/2.,type,par)*pdf((z/a-z)/2.,type,par)*std::abs(z));
	    }
	  } else {
	    sAsym.push_back(sAsymOld.at(i/3));       // Store value from previous iteration
	  }
	}

	// Sum up weighted function values
	for(size_t i = 0; i < sAsym.size(); i++) {
	  double w = 1.;                       // Weight w from Simpson's rule
	  if( i > 0 && i < (sAsym.size() - 1) ) { // w = 1 for x0 and last node
	    if( i % 3 == 0 ) {                 // w = 2 for x3, x6, ...
	      w = 2.;
	    } else {
	      w = 3.;
	    }
	  }
	  asym += w*(sAsym.at(i));
	}
	// Apply overall normalization
	asym *= (3.*h/16./a/a);                 

	if( asymOld ) eps = std::abs((asym - asymOld) / asymOld);
	nIter++;
      }
    }

    return asym;
  }
      


  double ResponseFunction::random(Type type, const std::vector<double> &par) const {
    double f = 0.;
    if( type == Gauss ) {
      assert( par.size() == 2 );
      f = randomGauss(par);
    } else if( type == CrystalBall ) {
      assert( par.size() == 4 );
      f = randomCrystalBall(par);
    } else if( type == TruncCrystalBall ) {
      assert( par.size() == 5 );
      f = randomTruncCrystalBall(par);
    }

    return f;
  }



  //! A Crystal Ball pdf
  // ------------------------------------------------------------------------
  ResponseFunction::CrystalBallFunction::CrystalBallFunction()
    : rand_(new TRandom3(0)) {
    par_.push_back(1.);
    par_.push_back(0.1);
    par_.push_back(1.);
    par_.push_back(2.);
  }

  ResponseFunction::CrystalBallFunction::CrystalBallFunction(const std::vector<double> &par)
    : rand_(new TRandom3(0)), par_(par) {
    assert( par_.size() == 4 );
  }
	  
  ResponseFunction::CrystalBallFunction::~CrystalBallFunction() {
    delete rand_;
  }

    
    
  //! The Crystal Ball function value (not normalized)
  // ------------------------------------------------------------------------
  double ResponseFunction::CrystalBallFunction::value(double x, double mean, double sigma, double alpha, double n) const {
    double f = 0.;
    if( x > 0. ) {
      double u = (x - mean)/sigma;
      if( u > -alpha ) {             // Gaussian part
	f = exp(-0.5*u*u);
      } else {                       // Powerlaw part
	f = exp(-0.5*alpha*alpha)/pow(1.-alpha*alpha/n-alpha*u/n,n); //A*pow(B-u,-n);
      }
    }
    
    return f;
  }


  //! The inverse of the integral over Crystal Ball function 
  //! from 0 to infinity
  // ------------------------------------------------------------------------
  double ResponseFunction::CrystalBallFunction::norm(double mean, double sigma, double alpha, double n) const {
    double m = n-1.;
    double k = n*sigma*exp(-0.5*alpha*alpha)/m/alpha;
    double norm = sigma*sqrt(M_PI/2.)*( 1. + erf(alpha/sqrt(2)) );
    if( n == 1. ) {
      double B = n/alpha - alpha;
      norm += k/pow(1-alpha*alpha/n,m)*log( (B+mean/sigma)/(B+alpha) );
    } else {
      norm += k*( 1. - pow( 1 + alpha/n*( mean/sigma - alpha ),-m ) );
    }
    norm = 1./norm;
    if( norm < 1E-10 ) norm = 1E-10;
    return norm;
  }


  // ------------------------------------------------------------------------
  double ResponseFunction::CrystalBallFunction::random(double mean, double sigma, double alpha, double n) const {
    double max = pdf(mean,mean,sigma,alpha,n);
    double x = -1.;
    do {
      x = rand_->Uniform(0.,2.);
    } while( pdf(x,mean,sigma,alpha,n) < rand_->Uniform(0.,1.)*max );
    return x;
  }



  // ------------------------------------------------------------------------
  double ResponseFunction::CrystalBallFunction::integral(double min, double max, double mean, double sigma, double alpha, double n) const {
    double m = n - 1.;
    double c = mean - alpha*sigma;
    
    double in = 0.;
    if( min > c ) {
      // Integral from Gaussian part
      in = sigma*sqrt(M_PI/2.)*( erf((mean-min)/sqrt(2)/sigma) - erf((mean-max)/sqrt(2)/sigma) );
    } else if( max < c ) {
      // Integral from powerlaw part
      double k = n*sigma*exp(-0.5*alpha*alpha)/m/alpha;
      if( n == 1. ) {
	double B = n/alpha - alpha;
	in = k/pow(1-alpha*alpha/n,m)*log( (B+(mean-min)/sigma)/(B+(mean-max)/sigma) );
      } else {
	in = k*( pow( 1 + alpha/n*( (mean-max)/sigma - alpha ),-m )
		 -pow( 1 + alpha/n*( (mean-min)/sigma - alpha ),-m ) );
      }
    } else {
      // Integral from both parts
      double k = n*sigma*exp(-0.5*alpha*alpha)/m/alpha;
      in = sigma*sqrt(M_PI/2.)*( erf(alpha/sqrt(2)) - erf((mean-max)/sqrt(2)/sigma) );
      if( n == 1. ) {
	double B = n/alpha - alpha;
	in += k/pow(1-alpha*alpha/n,m)*log( (B+(mean-min)/sigma)/(B+alpha) );
      } else {
	in += k*( 1. - pow( 1 + alpha/n*( (mean-min)/sigma - alpha ),-m ) );
      }
    }

    return in;
  }


  // ------------------------------------------------------------------------
  double ResponseFunction::CrystalBallFunction::truncPdf(double x, double mean, double sigma, double alpha, double n, double min) const {
    double f = value(x,mean,sigma,alpha,n);
    if( f ) {
      double norm = integral(min,2.,mean,sigma,alpha,n);
      if( norm > 0. ) {
	f /=norm;
      }
    }
    return f;
  }


  // ------------------------------------------------------------------------
  double ResponseFunction::CrystalBallFunction::truncRandom(double mean, double sigma, double alpha, double n, double min) const {
    double max = truncPdf(mean,mean,sigma,alpha,n,min);
    double x = -1.;
    do {
      x = rand_->Uniform(0.,2.);
    } while( truncPdf(x,mean,sigma,alpha,n,min) < rand_->Uniform(0.,1.)*max );
    return x;
  }


  // ------------------------------------------------------------------------
  double ResponseFunction::CrystalBallFunction::truncValue(double x, double mean, double sigma, double alpha, double n, double min) const {
    double f = 0.;
    if( x > min ) {
      f = value(x,mean,sigma,alpha,n);
    }
    return f;
  }

}





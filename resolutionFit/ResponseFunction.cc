#include "ResponseFunction.h"

#include <cassert>
#include <cmath>


namespace resolutionFit {
  double ResponseFunction::operator()(double x, const std::vector<double> &par) const {
    double f = 0.;
    if( type() == Gauss ) {
      assert( par.size() == 1 );
      f = pdfGauss(x,1.,par[0]);
    } else if( type() == CrystalBall ) {
      assert( par.size() == 3 );
      f = pdfCrystalBall(x,1.,par[0],par[1],par[2]);
    }

    return f;
  }


  int ResponseFunction::nPars() const {
    int n = 0;
    if( type() == Gauss ) n = 1;
    else if( type() == CrystalBall ) n = 3;

    return n;
  }


  double ResponseFunction::pdfGauss(double x, double mean, double sigma) const {
    double f = 0.;
    if( x > 0. ) {
      double u = (x-mean)/sigma;
      f =  exp(-0.5*u*u)/sqrt(2.*M_PI)/sigma;
    }
    
    return f;
  }


  double ResponseFunction::pdfCrystalBall(double x, double mean, double sigma, double alpha, double n) const {
    double f = 0.;
    if( x > 0. ) {
      double u = (x - mean)/sigma;
      if( u > -alpha ) {             // Gaussian part
	f = exp(-0.5*u*u);
      } else {                       // Powerlaw part
	f = exp(-0.5*alpha*alpha)/pow(1.-alpha*alpha/n-alpha*u/n,n);
      }
      double norm = crystalBallInt(mean,sigma,alpha,n,0.,2.);
      if( norm ) f /= norm;
    }
    return f;
  }


  double ResponseFunction::crystalBallInt(double mean, double sigma, double alpha, double n, double min, double max) const {
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

}



#include <cmath>
#include <iostream>

double cb(double x, double mu, double sigma, double alpha, double n, double a, double b) {
  double f = 0.;

  if( x > a && x < b ) {
    double u = (x - mu)/sigma;
    double A = pow(n/alpha,n)*exp(-0.5*alpha*alpha);
    double B = n/alpha - alpha;
    
    if( u > -alpha ) {             // Gaussian part
      f = exp(-0.5*u*u);
    } else {                       // Powerlaw part
      f =  A*pow(B-u,-n);
    }
  }

  return f;
}


double in(double mu, double sigma, double alpha, double n, double a, double b) {
  double A = pow(n/alpha,n)*exp(-0.5*alpha*alpha);
  double B = n/alpha - alpha;
  double m = n - 1.;

  // Norm from Gaussian part
  double in = sigma*sqrt(M_PI/2.)*( erf(alpha/sqrt(2)) - erf((mu-b)/sqrt(2)/sigma) );
  // Norm from powerlaw part
  if( n == 1. ) {
    in += A*sigma/pow(B,m)*log( (B+(mu-a)/sigma)/(B+alpha) );
  } else {
    in += A*sigma/m*(1./pow(B+alpha,m) - 1./pow(B+(mu-a)/sigma,m));
  }

  return in;
}


double norm(double mu, double sigma, double alpha, double n, double a, double b) {
  return 1./in(mu,sigma,alpha,n,a,b);
}


double der(int i, double x, double mu, double sigma, double alpha, double n, double a, double b) {
  double d = 0.;
  
  double u = (x - mu)/sigma;
  double beta = (b-mu)/sigma;
  double B = n/alpha - alpha;
  double m = n - 1.;
  double c = a-1.;
  double s = n*sigma;
  double k = s - alpha*c - alpha*alpha*sigma;
  double ea = exp(-0.5*alpha*alpha);

  std::cout << "\nu " << u << std::endl;
  std::cout << "beta " << beta << std::endl;
  std::cout << "B " << B << std::endl;
  std::cout << "m " << m << std::endl;
  std::cout << "c " << c << std::endl;
  std::cout << "s " << s << std::endl;
  std::cout << "k " << k << std::endl;
  std::cout << "ea " << ea << std::endl;


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
      * norm(mu,sigma,alpha,n,a,b);
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
      * norm(mu,sigma,alpha,n,a,b);
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
      * norm(mu,sigma,alpha,n,a,b);
  }
  d *= norm(mu,sigma,alpha,n,a,b)*cb(x,mu,sigma,alpha,n,a,b);

  return d;
}


void crystalBall(double x = 1.15, double sigma = 0.07, double alpha = 1.8, double n = 5.6, double a = 0.2, double b = 1.6) {

  std::cout << "PARAMETER\n";
  std::cout << " x  " << x << std::endl;
  std::cout << " sigma  " << sigma << std::endl;
  std::cout << " alpha  " << alpha << std::endl;
  std::cout << " n  " << n << std::endl;
  std::cout << " a  " << a << std::endl;
  std::cout << " b  " << b << std::endl;

  std::cout << "\nFUNCTION\n";
  std::cout << " f  " << cb(x,1.,sigma,alpha,n,a,b) << std::endl;
  std::cout << " int  " << in(1.,sigma,alpha,n,a,b) << std::endl;
  std::cout << " N  " << norm(1.,sigma,alpha,n,a,b) << std::endl;

  std::cout << "\nDERIVATIVES\n";
  std::cout << " sigma  " << der(1,x,1.,sigma,alpha,n,a,b) << std::endl;
  std::cout << " alpha  " << der(2,x,1.,sigma,alpha,n,a,b) << std::endl;
  std::cout << " n  " << der(3,x,1.,sigma,alpha,n,a,b) << std::endl;
}

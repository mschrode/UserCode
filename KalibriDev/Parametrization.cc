//
//  $Id: Parametrization.cc,v 1.5 2010/04/13 13:38:24 mschrode Exp $
//
#include "Parametrization.h"


#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TMath.h"




// ------------------------------------------------------------------------
SmearGauss::SmearGauss(double tMin, double tMax, double xMin, double xMax, const std::vector<double>& parScales)
  : Parametrization(0,4,0,0),
    tMin_(tMin),
    tMax_(tMax),
    xMin_(xMin),
    xMax_(xMax),
    scale_(parScales),
    hPdfPtTrue_(0) {
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= xMin_ && xMin_ < xMax_ );
  assert( scale_.size() >= nJetPars() );

  //  hPdfPtTrue_ = new TH1F("hPdfPtTrue_","",10000,tMin,tMax);
//   //Fill values of unnormalised truth pdf
//   double *truePar = new double[nJetPars()];
//   // Sigma parameters
//   truePar[0] = 4.;
//   truePar[1] = 1.2;
//   truePar[2] = 0.05;
//   // Spectrum parameters
//   truePar[3] = 80.;
//   for(int bin = 1; bin <= hPdfPtTrue_->GetNbinsX(); bin++) {
//     double ptTrue = hPdfPtTrue_->GetBinCenter(bin);
//     hPdfPtTrue_->SetBinContent(bin,pdfPtTrueNotNorm(ptTrue,truePar));
//   }
//   // Normalise values of truth pdf
//   hPdfPtTrue_->Scale(1./hPdfPtTrue_->Integral("width"));
//   delete [] truePar;

//  TFile file("controlPlots/truthSpectrum.root","READ");
//   file.GetObject("hPtGenJet1",hPdfPtTrue_);
//   if( !hPdfPtTrue_ ) {
//     std::cerr << "ERROR: No histogram found in file '" << file.GetName() << "'\n";
//     exit(1);
//   } else {
//     hPdfPtTrue_->SetDirectory(0);
//     hPdfPtTrue_->SetName("hPdfPtTrue");
//     int binMin = hPdfPtTrue_->FindBin(tMin_);
//     int binMax = hPdfPtTrue_->FindBin(tMax_);
//     std::cout << "Integal: " << hPdfPtTrue_->Integral(binMin,binMax,"width") << std::endl;
//   }
//   file.Close();

  print();
}


SmearGauss::~SmearGauss() { if( hPdfPtTrue_ ) delete hPdfPtTrue_; }


// ------------------------------------------------------------------------
double SmearGauss::pdfPtMeasJet1(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
  double s = sigma(ptTrue,par);
  //  s *= corr3rdJet(ptTrue);
  double u = (ptMeas - ptTrue)/s;
  double norm = sqrt(M_PI/2.)*s*( erf((xMax_-ptTrue)/sqrt(2.)/s) - erf((xMin_-ptTrue)/sqrt(2.)/s) );
  // This should be caught more cleverly
  if( norm < 1E-10 ) norm = 1E-10;

  return exp(-0.5*u*u)/norm;
}

// ------------------------------------------------------------------------
double SmearGauss::pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
  double s = sigma(ptTrue,par);
  //  s *= corr3rdJet(ptTrue);
  double u = (ptMeas - ptTrue)/s;
  double norm = sqrt(M_PI/2.)*s*( 1. + erf(ptTrue/sqrt(2.)/s) );
  // This should be caught more cleverly
  if( norm < 1E-10 ) norm = 1E-10;

  return exp(-0.5*u*u)/norm;
}


// ------------------------------------------------------------------------
double  SmearGauss::pdfPtTrue(double ptTrue, const double *par) const {
  double pdf = 0.;
  if( tMin_ < ptTrue && ptTrue < tMax_ ) {
//     double m = 1.-exponent(par);
//     double norm = ( m == 0. ? log(tMax_/tMin_) : (pow(tMax_,m)-pow(tMin_,m))/m );
//     pdf = 1./pow(ptTrue,exponent(par))/norm;

//     // ToyMC parametrization of spectrum
//     double tau = exponent(par);
//     double norm = tau*(exp(-tMin_/tau)-exp(-tMax_/tau));
//     pdf = exp(-ptTrue/tau)/norm;

    double norm = exp(-specPar(par,0))*(exp(-specPar(par,1)*tMin_)-exp(-specPar(par,1)*tMax_))/specPar(par,1);
    norm += exp(-specPar(par,2))*(exp(-specPar(par,3)*tMin_)-exp(-specPar(par,3)*tMax_))/specPar(par,3);
    norm += exp(-specPar(par,4))*(exp(-specPar(par,5)*tMin_)-exp(-specPar(par,5)*tMax_))/specPar(par,5);

    pdf = exp(-specPar(par,0)-specPar(par,1)*ptTrue);
    pdf += exp(-specPar(par,2)-specPar(par,3)*ptTrue);
    pdf += exp(-specPar(par,4)-specPar(par,5)*ptTrue);
    pdf /= norm;

    //  pdf = 1./(tMax_-tMin_);
  }
  return pdf;
  
//  return hPdfPtTrue_->Interpolate(ptTrue);
}


// ------------------------------------------------------------------------
double SmearGauss::pdfResponse(double r, double ptTrue, const double *par) const {
  double s = sigma(ptTrue,par)/ptTrue;
  double u = (r - 1.)/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
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
double SmearGauss::corr3rdJet(double ptTrue) const {
//   std::vector<double> a(3);
//   a[0] = 1.;
//   a[1] = 8E-4;
//   a[2] = -7E-7;
//   return a[0] + a[1]*ptTrue + a[2]*ptTrue*ptTrue;

  return 0.;
}


// ------------------------------------------------------------------------
double  SmearGauss::pdfPtTrueNotNorm(double ptTrue, const double *par) const {
  double pdf = 0.;

  // Underlying truth pdf

  // ToyMC parametrization of spectrum, normalized from 0 to infty
  //pdf = exp(-ptTrue/par[3])/par[3];

//   double norm = exp(-specPar(par,0))*(exp(-specPar(par,1)*tMin_)-exp(-specPar(par,1)*tMax_))/specPar(par,1);
//   norm += exp(-specPar(par,2))*(exp(-specPar(par,3)*tMin_)-exp(-specPar(par,3)*tMax_))/specPar(par,3);
//   norm += exp(-specPar(par,4))*(exp(-specPar(par,5)*tMin_)-exp(-specPar(par,5)*tMax_))/specPar(par,5);

//   pdf = exp(-specPar(par,0)-specPar(par,1)*ptTrue);
//   pdf += exp(-specPar(par,2)-specPar(par,3)*ptTrue);
//   pdf += exp(-specPar(par,4)-specPar(par,5)*ptTrue);
//   pdf /= norm;

//   // Convolution with cuts on ptdijet
//   double s = sigma(ptTrue,par)/sqrt(2.);
//   double c = 0.5*( erf((xMax_-ptTrue)/s/sqrt(2.)) - erf((xMin_-ptTrue)/s/sqrt(2.)) );

  // Convolution with cuts on 1. jet (randomly assigned)
//   double s = sqrt( par[0]*par[0] + par[1]*par[1]*ptTrue + par[2]*par[2]*ptTrue*ptTrue );
//   double c = 0.5*( erf((xMax_-ptTrue)/s/sqrt(2.)) - erf((xMin_-ptTrue)/s/sqrt(2.)) );

  double c = 1.;

  return c*pdf;
}



// ------------------------------------------------------------------------
void SmearGauss::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  Probability density of ptTrue spectrum:\n";
  std::cout << "    " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
  std::cout << "    " << xMin_ << " < pt < " << xMax_ << " GeV\n";
  std::cout << std::endl;
}



// ------------------------------------------------------------------------
SmearGaussPtBin::SmearGaussPtBin(double tMin, double tMax, double xMin, double xMax, const std::vector<double> &parScales, const std::vector<double> &startPar)
  : Parametrization(0,10,0,0),
    tMin_(tMin),
    tMax_(tMax),
    xMin_(xMin),
    xMax_(xMax),
    scale_(parScales) {
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= xMin_ && xMin_ < xMax_ );
  assert( scale_.size() >= nJetPars() );
  assert( startPar.size() >= nJetPars() );

  print();

  TFile file("config/truthSpectrum_Spring10_Eta1.root","READ");
  file.GetObject("hPtGen",hPdfPtTrue_);
  if( !hPdfPtTrue_ ) {
    std::cerr << "ERROR: No histogram found in file '" << file.GetName() << "'\n";
    exit(1);
  } else {
    hPdfPtTrue_->SetDirectory(0);
    hPdfPtTrue_->SetName("hPdfPtTrue");
    int binMin = hPdfPtTrue_->FindBin(tMin_);
    int binMax = hPdfPtTrue_->FindBin(tMax_);
    if( hPdfPtTrue_->Integral(binMin,binMax,"width") )
      hPdfPtTrue_->Scale(1./hPdfPtTrue_->Integral(binMin,binMax,"width"));
    std::cout << "Integal: " << hPdfPtTrue_->Integral(binMin,binMax,"width") << std::endl;
  }
  file.Close();
  
  hashTablePdfPtTrue_ = new TH1D("hashTablePdfPtTrue_","",5000,tMin_,tMax_);
  hashPdfPtTrue(&(startPar.front()));
}

SmearGaussPtBin::~SmearGaussPtBin() { 
  delete hashTablePdfPtTrue_;
  if( hPdfPtTrue_ ) delete hPdfPtTrue_;
}



// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfPtMeasJet1(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
  double pdf = 0.;
  if( xMin_ < ptMeas && ptMeas < xMax_ ) {
    double s = sigma(par);
    double u = (ptMeas - ptTrue)/s;
    double norm = sqrt(M_PI/2.)*s*( erf((xMax_-ptTrue)/sqrt(2.)/s) - erf((xMin_-ptTrue)/sqrt(2.)/s) );
    // This should be caught more cleverly
    if( norm < 1E-10 ) norm = 1E-10;
    
    pdf = exp(-0.5*u*u)/norm; 
  }

  return pdf;
}



// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
  double s = sigma(par);
  double u = (ptMeas - ptTrue)/s;
  double norm = sqrt(M_PI/2.)*s*( 1. + erf(ptTrue/sqrt(2.)/s) );
  // This should be caught more cleverly
  if( norm < 1E-10 ) norm = 1E-10;

  return exp(-0.5*u*u)/norm;
}



// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfPtTrue(double ptTrue, const double *par) const {
  return hashTablePdfPtTrue_->GetBinContent(hashTablePdfPtTrue_->FindBin(ptTrue));
}


// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfResponse(double r, double ptTrue, const double *par) const {
  double s = sigma(par)/ptTrue;
  double u = (r - 1.)/s;
  double cut = (1.+erf(ptTrue/sqrt(2.)/s))/2.;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s/cut;
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
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


// ------------------------------------------------------------------------
void SmearGaussPtBin::hashPdfPtTrue(const double *par) const {
  std::cout << "  Hashing truth pdf... " << std::flush;

  // Loop over tMin_ < ptTrue < tMax_ values
  for(int bin = 1; bin <= hashTablePdfPtTrue_->GetNbinsX(); bin++) {
    double ptTrue = hashTablePdfPtTrue_->GetBinCenter(bin);
    // Store (un-normalized) truth pdf for
    // this value of ptTrue in hash table
    hashTablePdfPtTrue_->SetBinContent(bin,pdfPtTrueNotNorm(ptTrue,par));
  }
  // Normalise values of truth pdf
  hashTablePdfPtTrue_->Scale(1./hashTablePdfPtTrue_->Integral("width"));

  std::cout << "ok" << std::endl;
}



// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfPtTrueNotNorm(double ptTrue, const double *par) const {
  double pdf = underlyingPdfPtTrue(ptTrue,par);
  // Convolution with cuts on 1. jet
  double s = sqrt( specSigmaPar(par,0)*specSigmaPar(par,0) +
		   specSigmaPar(par,1)*specSigmaPar(par,1)*ptTrue +
		   specSigmaPar(par,2)*specSigmaPar(par,2)*ptTrue*ptTrue );
  double c = 0.5*( erf((xMax_-ptTrue)/s/sqrt(2.)) - erf((xMin_-ptTrue)/s/sqrt(2.)) );

  return c*pdf;
}


// ------------------------------------------------------------------------
double SmearGaussPtBin::underlyingPdfPtTrue(double ptTrue, const double *par) const {
//   // For ToyMC
//   double tau = specSlopePar(par,0);
//   double norm = tau*(exp(-tMin_/tau)-exp(-tMax_/tau));
//   double pdf = exp(-ptTrue/tau)/norm;

//   // For QCD
//   double norm = exp(-specSlopePar(par,0))
//     *(exp(-specSlopePar(par,1)*tMin_)-exp(-specSlopePar(par,1)*tMax_))/specSlopePar(par,1);
//   norm += exp(-specSlopePar(par,2))
//     *(exp(-specSlopePar(par,3)*tMin_)-exp(-specSlopePar(par,3)*tMax_))/specSlopePar(par,3);
//   norm += exp(-specSlopePar(par,4))
//     *(exp(-specSlopePar(par,5)*tMin_)-exp(-specSlopePar(par,5)*tMax_))/specSlopePar(par,5);
  
//   double pdf = exp(-specSlopePar(par,0)-specSlopePar(par,1)*ptTrue);
//   pdf += exp(-specSlopePar(par,2)-specSlopePar(par,3)*ptTrue);
//   pdf += exp(-specSlopePar(par,4)-specSlopePar(par,5)*ptTrue);

//   return pdf / norm;

  return hPdfPtTrue_->Interpolate(ptTrue);
  //return hPdfPtTrue_->GetBinContent(hPdfPtTrue_->FindBin(ptTrue));
}


// ------------------------------------------------------------------------
void SmearGaussPtBin::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
  std::cout << "  " << xMin_ << " < ptMeas < " << xMax_ << " GeV\n";
  std::cout << std::endl;
}



// ------------------------------------------------------------------------
SmearGaussExtrapolation::SmearGaussExtrapolation(double tMin, double tMax, double xMin, double xMax, const std::vector<double>& parScales)
  : Parametrization(0,6,0,0),
    tMin_(tMin),
    tMax_(tMax),
    xMin_(xMin),
    xMax_(xMax),
    scale_(parScales),
    hPdfPtTrue_(0) {
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= xMin_ && xMin_ < xMax_ );
  assert( scale_.size() >= nJetPars() );

  print();
}


SmearGaussExtrapolation::~SmearGaussExtrapolation() { if( hPdfPtTrue_ ) delete hPdfPtTrue_; }


// ------------------------------------------------------------------------
double SmearGaussExtrapolation::pdfPtMeasJet1(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
  double s = sigmaWithOffset(ptTrue,pt3Rel,par);
  double u = (ptMeas - ptTrue)/s;
  double norm = sqrt(M_PI/2.)*s*( erf((xMax_-ptTrue)/sqrt(2.)/s) - erf((xMin_-ptTrue)/sqrt(2.)/s) );
  // This should be caught more cleverly
  if( norm < 1E-10 ) norm = 1E-10;

  return exp(-0.5*u*u)/norm;
}

// ------------------------------------------------------------------------
double SmearGaussExtrapolation::pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
  double s = sigmaWithOffset(ptTrue,pt3Rel,par);
  double u = (ptMeas - ptTrue)/s;
  double norm = sqrt(M_PI/2.)*s*( 1. + erf(ptTrue/sqrt(2.)/s) );
  // This should be caught more cleverly
  if( norm < 1E-10 ) norm = 1E-10;

  return exp(-0.5*u*u)/norm;
}


// ------------------------------------------------------------------------
double  SmearGaussExtrapolation::pdfPtTrue(double ptTrue, const double *par) const {
  double pdf = 0.;
  if( tMin_ < ptTrue && ptTrue < tMax_ ) {
//     double m = 1.-exponent(par);
//     double norm = ( m == 0. ? log(tMax_/tMin_) : (pow(tMax_,m)-pow(tMin_,m))/m );
//     pdf = 1./pow(ptTrue,exponent(par))/norm;

    // ToyMC parametrization of spectrum
    double tau = exponent(par);
    double norm = tau*(exp(-tMin_/tau)-exp(-tMax_/tau));
    pdf = exp(-ptTrue/tau)/norm;

//     double norm = exp(-specPar(par,0))*(exp(-specPar(par,1)*tMin_)-exp(-specPar(par,1)*tMax_))/specPar(par,1);
//     norm += exp(-specPar(par,2))*(exp(-specPar(par,3)*tMin_)-exp(-specPar(par,3)*tMax_))/specPar(par,3);
//     norm += exp(-specPar(par,4))*(exp(-specPar(par,5)*tMin_)-exp(-specPar(par,5)*tMax_))/specPar(par,5);

//     pdf = exp(-specPar(par,0)-specPar(par,1)*ptTrue);
//     pdf += exp(-specPar(par,2)-specPar(par,3)*ptTrue);
//     pdf += exp(-specPar(par,4)-specPar(par,5)*ptTrue);
//     pdf /= norm;

    //  pdf = 1./(tMax_-tMin_);
  }
  return pdf;
  
//  return hPdfPtTrue_->Interpolate(ptTrue);
}


// ------------------------------------------------------------------------
double SmearGaussExtrapolation::pdfResponse(double r, double ptTrue, const double *par) const {
  double s = sigma(ptTrue,par)/ptTrue;
  double u = (r - 1.)/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


// ------------------------------------------------------------------------
double SmearGaussExtrapolation::pdfDijetAsym(double a, double ptTrue, const double *par) const {
  double s = sigma(ptTrue,par)/ptTrue/sqrt(2.);
  double u = a/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


// ------------------------------------------------------------------------
double SmearGaussExtrapolation::pdfResponseError(double r, double ptTrue, const double *par, const double *cov, const std::vector<int> &covIdx) const {
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
double SmearGaussExtrapolation::pdfResponseDeriv(double r, double ptTrue, const double *par, int i) const {
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
void SmearGaussExtrapolation::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  Probability density of ptTrue spectrum:\n";
  std::cout << "    " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
  std::cout << "    " << xMin_ << " < pt < " << xMax_ << " GeV\n";
  std::cout << std::endl;
}



// ------------------------------------------------------------------------
SmearGaussImbalance::SmearGaussImbalance(double tMin, double tMax, double xMin, double xMax, const std::vector<double>& parScales)
  : Parametrization(0,9,0,0),
    tMin_(tMin),
    tMax_(tMax),
    xMin_(xMin),
    xMax_(xMax),
    scale_(parScales),
    dPtTrue_(0.1),
    hPdfPtTrue_(0) {
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= xMin_ && xMin_ < xMax_ );
  assert( scale_.size() >= nJetPars() );

  parMin_ = std::vector<double>(3);
  parMin_[0] = 0.;
  parMin_[1] = 1E-3;
  parMin_[2] = 1E-3;

  print();
}


SmearGaussImbalance::~SmearGaussImbalance() { if( hPdfPtTrue_ ) delete hPdfPtTrue_; }


// ------------------------------------------------------------------------
double SmearGaussImbalance::pdfPtMeasJet1(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
  double s = sigma(ptTrue,par);
  double u = (ptMeas - ptTrue)/s;
  double norm = sqrt(M_PI/2.)*s*( erf((xMax_-ptTrue)/sqrt(2.)/s) - erf((xMin_-ptTrue)/sqrt(2.)/s) );
  // This should be caught more cleverly
  if( norm < 1E-10 ) norm = 1E-10;

  return exp(-0.5*u*u)/norm;
}

// ------------------------------------------------------------------------
double SmearGaussImbalance::pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
//   double min = ptTrue2Min(ptTrue);
//   double max = ptTrue2Max(ptTrue);
//   double sMin = sigma(min,par);
//   double sMax = sigma(max,par);
//   double norm = 2.*(max-min);
//   double pa = param(1,par);

//   double pdf = erfc((min-ptMeas)/sqrt(2.)/sMin) - erfc((max-ptMeas)/sqrt(2.)/sMax);

//   double factor1 = exp(2.*ptMeas/pa/pa);
//   double factor2 = erf((min+ptMeas)/sqrt(2.)/sMin) - erf((max+ptMeas)/sqrt(2.)/sMax);

//   if( factor1 != HUGE_VAL ) {
//     pdf -= factor1*factor2;
//   } else if( factor2 > 0. ) {
//     std::cerr << "WARNING HUGE_VAL in SmearGaussImbalance::pdfPtMeasJet2\n";
//   }
//   pdf /= norm;

//   return pdf;

  double pdf = 0.;
  double s2 = 1.02*sqrt(ptTrue);
  double t2Min = ptTrue - 5*s2;
  double t2Max = ptTrue + 5*s2;
  int nt2Bins = 10;
  double dt2 = (t2Max-t2Min)/nt2Bins;
  for(int t2Bin = 0; t2Bin < nt2Bins; ++t2Bin) {
    double t2 = t2Min + (t2Bin+0.5)*dt2;
    double sx = sigma(t2,par);
    double dt = (t2-ptTrue)/s2;
    double dx = (ptMeas-t2)/sx;
    pdf += exp(-0.5*dt*dt)*exp(-0.5*dx*dx)/sx;
  }
  pdf *= dt2/s2/2./M_PI;

  return pdf;


//   double s = sigma(ptTrue,par);
//   double u = (ptMeas - ptTrue)/s;
//   double norm = sqrt(M_PI/2.)*s*( erf((xMax_-ptTrue)/sqrt(2.)/s) - erf((xMin_-ptTrue)/sqrt(2.)/s) );
//   // This should be caught more cleverly
//   if( norm < 1E-10 ) norm = 1E-10;

//   return exp(-0.5*u*u)/norm;
}


// ------------------------------------------------------------------------
double  SmearGaussImbalance::pdfPtTrue(double ptTrue, const double *par) const {
  double pdf = 0.;
  if( tMin_ < ptTrue && ptTrue < tMax_ ) {
//     // ToyMC parametrization of spectrum
//     double tau = exponent(par);
//     double norm = tau*(exp(-tMin_/tau)-exp(-tMax_/tau));
//     pdf = exp(-ptTrue/tau)/norm;

// Triple exponential spectrum
    double norm = exp(-specPar(par,0))*(exp(-specPar(par,1)*tMin_)-exp(-specPar(par,1)*tMax_))/specPar(par,1);
    norm += exp(-specPar(par,2))*(exp(-specPar(par,3)*tMin_)-exp(-specPar(par,3)*tMax_))/specPar(par,3);
    norm += exp(-specPar(par,4))*(exp(-specPar(par,5)*tMin_)-exp(-specPar(par,5)*tMax_))/specPar(par,5);

    pdf = exp(-specPar(par,0)-specPar(par,1)*ptTrue);
    pdf += exp(-specPar(par,2)-specPar(par,3)*ptTrue);
    pdf += exp(-specPar(par,4)-specPar(par,5)*ptTrue);
    pdf /= norm;
  }
  return pdf;
}


// ------------------------------------------------------------------------
double SmearGaussImbalance::pdfResponse(double r, double ptTrue, const double *par) const {
  double s = sigma(ptTrue,par)/ptTrue;
  double u = (r - 1.)/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


// ------------------------------------------------------------------------
double SmearGaussImbalance::pdfDijetAsym(double a, double ptTrue, const double *par) const {
  double s = sigma(ptTrue,par)/ptTrue/sqrt(2.);
  double u = a/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


// ------------------------------------------------------------------------
double SmearGaussImbalance::pdfResponseError(double r, double ptTrue, const double *par, const double *cov, const std::vector<int> &covIdx) const {
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
double SmearGaussImbalance::pdfResponseDeriv(double r, double ptTrue, const double *par, int i) const {
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
void SmearGaussImbalance::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
  std::cout << "  " << xMin_ << " < pt < " << xMax_ << " GeV\n";
  std::cout << "  dPtTrue = " << dPtTrue_ << std::endl;
  std::cout << std::endl;
}



// ------------------------------------------------------------------------
SmearCrystalBallPtBin::SmearCrystalBallPtBin(double tMin, double tMax, double xMin, double xMax, const std::vector<double> &parScales, const std::vector<double> &startPar)
  : Parametrization(0,9,0,0),
    tMin_(tMin),
    tMax_(tMax),
    xMin_(xMin),
    xMax_(xMax),
    scale_(parScales) {
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= xMin_ && xMin_ < xMax_ );
  assert( scale_.size() >= nJetPars() );
  
  print();

  hashTablePdfPtTrue_ = new TH1D("hashTablePdfPtTrue_","",5000,tMin_,tMax_);
  hashPdfPtTrue(&(startPar.front()));
}

SmearCrystalBallPtBin::~SmearCrystalBallPtBin() {
  delete hashTablePdfPtTrue_;
}


// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::pdfPtMeasJet1(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
  double pdf = 0.;
  if( xMin_ < ptMeas && ptMeas < xMax_ ) {
    double norm = crystalBallInt(ptTrue,sigma(par),alpha(par),n(par),xMin_,xMax_);
    if( norm < 1E-10 ) norm = 1E-10;
    pdf = crystalBallFunc(ptMeas,ptTrue,sigma(par),alpha(par),n(par)) / norm;
  }
  return pdf;
}


// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
  return crystalBallNorm(ptTrue,sigma(par),alpha(par),n(par))
    *crystalBallFunc(ptMeas,ptTrue,sigma(par),alpha(par),n(par));
}


// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::pdfPtTrue(double ptTrue, const double *par) const {
  return hashTablePdfPtTrue_->GetBinContent(hashTablePdfPtTrue_->FindBin(ptTrue));
}


// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::pdfResponse(double r, double ptTrue, const double *par) const {
  return crystalBallFunc(r,1.,sigma(par)/ptTrue,alpha(par),n(par))*crystalBallNorm(1.,sigma(par)/ptTrue,alpha(par),n(par));
}


// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::pdfDijetAsym(double a, double ptTrue, const double *par) const {
  double s = sigma(par)/ptTrue/sqrt(2.);
  double u = a/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


//! The Crystal Ball function value (not normalized)
// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::crystalBallFunc(double x, double mean, double sigma, double alpha, double n) const {
  double f = 0.;
  if( x > 0. ) {
    double u = (x - mean)/sigma;
    if( u > -alpha ) {             // Gaussian part
      f = exp(-0.5*u*u);
    } else {                       // Powerlaw part
      f = exp(-alpha*alpha)/pow(1.-alpha*alpha/n-alpha*u/n,n); //A*pow(B-u,-n);
    }
  }

  return f;
}


//! The inverse of the integral over Crystal Ball function 
//! from 0 to infinity
// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::crystalBallNorm(double mean, double sigma, double alpha, double n) const {
  double m = n-1.;
  double k = n*sigma*exp(-0.5*alpha*alpha)/m/alpha;
  double norm = sigma*sqrt(M_PI/2.)*( 1. + erf(alpha/sqrt(2)) );
  if( n == 1. ) {
    double B = n/alpha - alpha;
    norm += k/pow(1-alpha*alpha/n,m)*log( (B+mean/sigma)/(B+alpha) );
  } else {
    norm += k*( 1. - pow( 1 + alpha/n*( mean/sigma - alpha ),-m ) );
  }
  if( norm < 1E-10 ) norm = 1E-10;
  return norm;
}


// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::crystalBallInt(double mean, double sigma, double alpha, double n, double min, double max) const {
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
void SmearCrystalBallPtBin::hashPdfPtTrue(const double *par) const {
  std::cout << "  Hashing truth pdf... " << std::flush;

  // Loop over tMin_ < ptTrue < tMax_ values
  for(int bin = 1; bin <= hashTablePdfPtTrue_->GetNbinsX(); bin++) {
    double ptTrue = hashTablePdfPtTrue_->GetBinCenter(bin);
    // Store (un-normalized) truth pdf for
    // this value of ptTrue in hash table
    hashTablePdfPtTrue_->SetBinContent(bin,pdfPtTrueNotNorm(ptTrue,par));
  }
  // Normalise values of truth pdf
  hashTablePdfPtTrue_->Scale(1./hashTablePdfPtTrue_->Integral("width"));

  std::cout << "ok" << std::endl;
}


// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::pdfPtTrueNotNorm(double ptTrue, const double *par) const {
  double pdf = underlyingPdfPtTrue(ptTrue,par);

  // Description of cuts on 1. jet pt
  double s = sqrt( specSigmaPar(par,0)*specSigmaPar(par,0) +
		   specSigmaPar(par,1)*specSigmaPar(par,1)*ptTrue +
		   specSigmaPar(par,2)*specSigmaPar(par,2)*ptTrue*ptTrue );

  // Assuming Crystal Ball resolution
  double c = crystalBallNorm(ptTrue,s,alpha(par),n(par))
    *crystalBallInt(ptTrue,s,alpha(par),n(par),xMin_,xMax_);


  return c*pdf;
}


// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::underlyingPdfPtTrue(double ptTrue, const double *par) const {
//   // For ToyMC
//   double tau = specSlopePar(par,0);
//   double norm = tau*(exp(-tMin_/tau)-exp(-tMax_/tau));
//   double pdf = exp(-ptTrue/tau)/norm;

  // For QCD
  double norm = exp(-specSlopePar(par,0))
    *(exp(-specSlopePar(par,1)*tMin_)-exp(-specSlopePar(par,1)*tMax_))/specSlopePar(par,1);
  norm += exp(-specSlopePar(par,2))
    *(exp(-specSlopePar(par,3)*tMin_)-exp(-specSlopePar(par,3)*tMax_))/specSlopePar(par,3);
  norm += exp(-specSlopePar(par,4))
    *(exp(-specSlopePar(par,5)*tMin_)-exp(-specSlopePar(par,5)*tMax_))/specSlopePar(par,5);
  
  double pdf = exp(-specSlopePar(par,0)-specSlopePar(par,1)*ptTrue);
  pdf += exp(-specSlopePar(par,2)-specSlopePar(par,3)*ptTrue);
  pdf += exp(-specSlopePar(par,4)-specSlopePar(par,5)*ptTrue);

  return pdf / norm;

//  return hPdfPtTrue_->Interpolate(ptTrue);
}


// ------------------------------------------------------------------------
void SmearCrystalBallPtBin::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
  std::cout << "  " << xMin_ << " < ptMeas < " << xMax_ << " GeV\n";
  std::cout << std::endl;
}

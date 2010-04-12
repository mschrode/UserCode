//
//  $Id: Parametrization.cc,v 1.9 2010/03/25 08:30:59 mschrode Exp $
//
#include "Parametrization.h"


#include "TF1.h"
#include "TFile.h"
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

  hPdfPtTrue_ = new TH1F("hPdfPtTrue_","",10000,tMin,tMax);
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
double SmearGauss::pdfPtMeasJet1(double ptMeas, double ptTrue, const double *par) const {
  double s = sigma(ptTrue,par);
  s *= corr3rdJet(ptTrue);
  double u = (ptMeas - ptTrue)/s;
  double norm = sqrt(M_PI/2.)*s*( erf((xMax_-ptTrue)/sqrt(2.)/s) - erf((xMin_-ptTrue)/sqrt(2.)/s) );
  // This should be caught more cleverly
  if( norm < 1E-10 ) norm = 1E-10;

  return exp(-0.5*u*u)/norm;
}

// ------------------------------------------------------------------------
double SmearGauss::pdfPtMeasJet2(double ptMeas, double ptTrue, const double *par) const {
  double s = sigma(ptTrue,par);
  s *= corr3rdJet(ptTrue);
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

//    pdf = 1./(tMax_-tMin_);
  }
  return pdf;
  
  //  return hPdfPtTrue_->GetBinContent(hPdfPtTrue_->FindBin(ptTrue));
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
  std::vector<double> a(3);
  a[0] = 1.;
  a[1] = 8E-4;
  a[2] = -7E-7;
  return a[0] + a[1]*ptTrue + a[2]*ptTrue*ptTrue;
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
//   double c = 0.5*( erf((xMax_-ptTrue)/s/sqrt(2.)) - erf((xMin_-ptTrue)/s/sqrt(2.)) );

  // Convolution with cuts on 1. jet (randomly assigned)
  double s = sqrt( par[0]*par[0] + par[1]*par[1]*ptTrue + par[2]*par[2]*ptTrue*ptTrue );
  double c = 0.5*( erf((xMax_-ptTrue)/s/sqrt(2.)) - erf((xMin_-ptTrue)/s/sqrt(2.)) );

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

  hashTablePdfPtTrue_ = new TH1D("hashTablePdfPtTrue_","",5000,tMin_,tMax_);
  hashPdfPtTrue(&(startPar.front()));
}

SmearGaussPtBin::~SmearGaussPtBin() { 
  delete hashTablePdfPtTrue_;
}



// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfPtMeasJet1(double ptMeas, double ptTrue, const double *par) const {
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
double SmearGaussPtBin::pdfPtMeasJet2(double ptMeas, double ptTrue, const double *par) const {
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
}


// ------------------------------------------------------------------------
void SmearGaussPtBin::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
  std::cout << "  " << xMin_ << " < ptMeas < " << xMax_ << " GeV\n";
  std::cout << std::endl;
}

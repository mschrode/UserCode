#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TRandom3.h"

#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"



const int nEvts_ = 50000;
const double aMin_ = -0.2;
const double aMax_ = 0.2;
const int nABins_ = 5000;
const double zMin_ = -1.;
const double zMax_ = 1.;
const double eps_ = 1E-4;
const int maxNIter_ = 20;


TRandom3 *rand_ = new TRandom3(0);


double pdfGauss(double x, const std::vector<double> &par);
double pdf(double x, const std::vector<double> &par);
double randomGauss(const std::vector<double> &par);
double random(const std::vector<double> &par);
double pdfAsymmetry(double a, const std::vector<double> &par);


void AsymmetryToyExample() {
  std::vector<double> pars;
  pars.push_back(1.);
  pars.push_back(0.1);


  // Measured asymmetry distribution
  TH1 *hMeas = util::HistOps::createTH1D("hMeas",21,aMin_,aMax_,"p_{T} asymmetry","","events");
  hMeas->SetMarkerStyle(20);
  for(int n = 0; n < nEvts_; ++n) {
    double m1 = random(pars);
    double m2 = random(pars);
    if( m1+m2 ) {
      double asym = (m1-m2)/(m1+m2);
      hMeas->Fill(asym);
    }
  }
  for(int aBin = 1; aBin <= hMeas->GetNbinsX(); ++aBin) {
    hMeas->SetBinError(aBin,0.);
  }


  // Predicted asymmetry distribution
  TH1 *hPred = static_cast<TH1D*>(hMeas->Clone("hPred"));
  hPred->Sumw2();
  hPred->Reset();
  hPred->SetMarkerStyle(0);
  hPred->SetLineColor(2);
  double deltaA = hMeas->GetNbinsX()*hMeas->GetBinWidth(1)/nABins_;
  for(int aBin = 0; aBin < nABins_; ++aBin) {
    double a = aMin_ + (aBin+0.5)*deltaA;
    double pdfA = pdfAsymmetry(a,pars);
    hPred->Fill(a,pdfA);
  }
  hPred->Scale(1.*nEvts_*hMeas->GetBinWidth(1)*hPred->GetNbinsX()/hPred->GetEntries());
  for(int aBin = 1; aBin <= hPred->GetNbinsX(); ++aBin) {
    hPred->SetBinError(aBin,sqrt(hPred->GetBinContent(aBin)));
  }


  // Draw histograms
  TCanvas *canAsym = new TCanvas("canAsym","Asymmetry",500,500);
  canAsym->cd();
  hPred->Draw("H");
  hMeas->Draw("PE1same");

  TCanvas *canRatio = new TCanvas("canRatio","Ratio",500,500);
  canRatio->cd();
  TH1 *hRatio = util::HistOps::createRatioPlot(hPred,hMeas,"Prediction / Measurement");
  TH1 *hFrame = util::HistOps::createRatioFrame(hRatio,"Prediction / Measurement");
  hFrame->Draw();
  hRatio->Draw("PE1same");


  // Chi2 goodness-of-fit test
  double chi2 = 0.;
  for(int aBin = 1; aBin <= hPred->GetNbinsX(); ++aBin) {
    chi2 += (hMeas->GetBinContent(aBin) - hPred->GetBinContent(aBin))*(hMeas->GetBinContent(aBin) - hPred->GetBinContent(aBin))/(hPred->GetBinContent(aBin));
  }
  int ndof = hPred->GetNbinsX() - 3; // Parameter mu, sigma, and normalisation (num entries)
  std::cout << "Chi2        = " << chi2 << std::endl;
  std::cout << "Chi2 / ndof = " << chi2 / ndof << std::endl;
  std::cout << "Prob        = " << TMath::Prob(chi2,ndof) << std::endl;
}




double pdfAsymmetry(double a, const std::vector<double> &par) {
  double h = zMax_-zMin_;
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
    while( eps > eps_ && nIter < maxNIter_ ) {
      // Init iteration
      asymOld = asym;
      asym = 0.;
      sAsymOld = sAsym;
      sAsym.clear();
      
      // In each iteration, split h into 3 new intervals
      h /= 3.;        
      
      // Loop over nodes xi i.e. interval borders
      for(int i = 0; i <= pow(3.0,nIter+1); ++i){
	double z = zMin_ + i*h;
	
	// Calculate pdf only at new nodes
	if( nIter == 0 || i % 3 != 0 ) {
	  if( z == 0 ) {
	    std::cerr << "WARNING: z = 0" << std::endl;
	    sAsym.push_back(0.);
	  } else {
	    sAsym.push_back(pdf((z+z/a)/2.,par)*pdf((z/a-z)/2.,par)*std::abs(z));
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



double pdf(double x, const std::vector<double> &par) {
  return pdfGauss(x,par);
}


double random(const std::vector<double> &par) {
  return randomGauss(par);
}


double pdfGauss(double x, const std::vector<double> &par) {
  assert( par.size() == 2 );
  return exp(-0.5*(x-par[0])*(x-par[0])/par[1]/par[1])/sqrt(2.*M_PI)/par[1];
}


double randomGauss(const std::vector<double> &par) {
  assert( par.size() == 2 );
  return rand_->Gaus(par[0],par[1]);
}

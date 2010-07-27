//  $Id: AsymmetryToyExample.C,v 1.5 2010/07/27 10:48:15 mschrode Exp $
//
//  Toy validation of pt asymmetry pdf
// -------------------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TRandom3.h"

#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"
#include "../util/utils.h"


const int nEvts_ = 50000;
const double aMin_ = -0.2;
const double aMax_ = 0.2;
const double zMin_ = -1.;
const double zMax_ = 1.;
const double eps_ = 1E-4;
const int maxNIter_ = 20;

TRandom3 *rand_ = new TRandom3(0);


double pdfGauss(double x, const std::vector<double> &par);
double pdf(double x, const std::vector<double> &par);
double randomGauss(double t, const std::vector<double> &par);
double random(double t, const std::vector<double> &par);
double pdfAsymmetry(double a, const std::vector<double> &par);
double sigma(double t, const std::vector<double> &par);


void AsymmetryToyExample() {
  util::StyleSettings::presentationNoTitle();

  // Measured asymmetry distribution
  std::vector<double> parSig;
  parSig.push_back(2.0);
  parSig.push_back(1.2);
  parSig.push_back(0.04);

  TH1 *hMeas = util::HistOps::createTH1D("hMeas",31,aMin_,aMax_,"p_{T} asymmetry","","events");
  hMeas->SetMarkerStyle(20);
  TH1 *hResp = util::HistOps::createTH1D("hResp",31,0.,2.,"p^{meas}_{T} / p^{true}_{T}","","jets");
  hResp->SetMarkerStyle(20);
  for(int n = 0; n < nEvts_; ++n) {
    double truth = 125.;//rand_->Uniform(500.,600.);
    double m1 = random(truth,parSig);
    double m2 = random(truth,parSig);
    hResp->Fill(m1/truth);
    hResp->Fill(m2/truth);
    if( m1+m2 ) {
      double asym = (m1-m2)/(m1+m2);
      hMeas->Fill(asym);
    }
  }
  for(int aBin = 1; aBin <= hMeas->GetNbinsX(); ++aBin) {
    hMeas->SetBinError(aBin,0.);
  }


  // Fit mean response
  hResp->Fit("gaus","0QI");
  TF1 *fResp = hResp->GetFunction("gaus");
  fResp->SetLineWidth(2);
  fResp->SetLineColor(4);
  std::vector<double> meanSig;
  meanSig.push_back(fResp->GetParameter(2));

  // Predicted asymmetry distribution
  TH1 *hPred = static_cast<TH1D*>(hMeas->Clone("hPred"));
  hPred->Sumw2();
  hPred->Reset();
  hPred->SetMarkerStyle(0);
  hPred->SetLineColor(2);
  for(int aBin = 1; aBin <= hPred->GetNbinsX(); ++aBin) {
    double aMin = hPred->GetXaxis()->GetBinLowEdge(aBin);
    double aMax = hPred->GetXaxis()->GetBinUpEdge(aBin);
    int n = 500;
    double deltaA = (aMax-aMin)/n;
    double pdfA = 0.;
    for(int i = 0; i < n; ++i) {
      double a = aMin + (i+0.5)*deltaA;
      pdfA += pdfAsymmetry(a,meanSig);
    }
    pdfA *= nEvts_*deltaA;
    hPred->SetBinContent(aBin,pdfA);
    hPred->SetBinError(aBin,sqrt(pdfA));
  }


  double aMin = -1.;
  double aMax = 1.;
  int n = 5000;
  double deltaA = (aMax-aMin)/n;
  double pdfA = 0.;
  for(int i = 0; i < n; ++i) {
    double a = aMin + (i+0.5)*deltaA;
    pdfA += pdfAsymmetry(a,meanSig);
  }
  pdfA *= deltaA;
  std::cout << "Int(pdfA) = " << pdfA << std::endl;


  // Chi2 goodness-of-fit test
  double chi2 = 0.;
  for(int aBin = 1; aBin <= hPred->GetNbinsX(); ++aBin) {
    chi2 += (hMeas->GetBinContent(aBin) - hPred->GetBinContent(aBin))*(hMeas->GetBinContent(aBin) - hPred->GetBinContent(aBin))/(hPred->GetBinContent(aBin));
  }
  int ndof = hPred->GetNbinsX() - 3; // Parameter mu, sigma, and normalisation (num entries)
  std::cout << "Chi2        = " << chi2 << std::endl;
  std::cout << "Chi2 / ndof = " << chi2 / ndof << std::endl;
  std::cout << "Prob        = " << TMath::Prob(chi2,ndof) << std::endl;



  // Draw histograms
  TString name = "ToyMC_AsymPred_ConstSigma_";

  TLegend *legResp = util::LabelFactory::createLegend(2);
  //legResp->AddEntry(hResp,"Toy MC,  #sigma(p^{true}_{T}),  p^{true}_{T} #in U(500,600)","P");
  legResp->AddEntry(hResp,"Toy MC,  #sigma = const,  p^{true}_{T} = 125","P");
  legResp->AddEntry(fResp,"Gaussian fit","L");
  util::HistOps::setYRange(hResp,2);
  TCanvas *canResp = new TCanvas("canResp","Response",500,500);
  canResp->cd();
  hResp->Draw("PE1");
  fResp->Draw("same");
  legResp->Draw("same");
  canResp->SaveAs(name+"Resp.eps","eps");

  TLegend *legAsym = util::LabelFactory::createLegend(3);
  legAsym->AddEntry(hMeas,"Toy MC measurement","P");
  //legAsym->AddEntry(hPred,"Prediction,  #bar{#sigma}","L");
  legAsym->AddEntry(hPred,"Prediction,  #sigma","L");
  util::LabelFactory::addExtraLegLine(legAsym,"#chi^{2}/ndof = "+util::toTString(chi2/ndof,3)+",  Prob = "+util::toTString(TMath::Prob(chi2,ndof),3));
  util::HistOps::setYRange(hPred,3);
  TCanvas *canAsym = new TCanvas("canAsym","Asymmetry",500,500);
  canAsym->cd();
  hPred->Draw("H");
  hMeas->Draw("PE1same");
  legAsym->Draw("same");
  canAsym->SaveAs(name+"Asym.eps","eps");

  TCanvas *canRatio = new TCanvas("canRatio","Ratio",500,500);
  canRatio->cd();
  TH1 *hRatio = util::HistOps::createRatioPlot(hPred,hMeas,"Prediction / Measurement");
  TH1 *hFrame = util::HistOps::createRatioFrame(hRatio,"Prediction / Measurement");
  TPaveText *txtRatio = util::LabelFactory::createPaveText(1);
  txtRatio->AddText("#chi^{2}/ndof = "+util::toTString(chi2/ndof,3)+",  Prob = "+util::toTString(TMath::Prob(chi2,ndof),3));
  util::HistOps::setYRange(hFrame,0.9,1.2);
  hFrame->Draw();
  hRatio->Draw("PE1same");
  txtRatio->Draw("same");
  canRatio->SaveAs(name+"Ratio.eps","eps");
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


double random(double t, const std::vector<double> &par) {
  return randomGauss(t,par);
}


double pdfGauss(double x, const std::vector<double> &par) {
  assert( par.size() == 1 );
  return exp(-0.5*(x-1.)*(x-1.)/par[0]/par[0])/sqrt(2.*M_PI)/par[0];
}


double randomGauss(double t, const std::vector<double> &par) {
  assert( par.size() == 3 );
  return rand_->Gaus(t,sigma(t,par));
}


double sigma(double t, const std::vector<double> &par) {
  assert( par.size() == 3 );
  return sqrt( par[0]*par[0] + par[1]*par[1]*t + par[2]*par[2]*t*t );
}

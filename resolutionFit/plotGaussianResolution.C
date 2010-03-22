// $Id: $

#include <cassert>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TPaveText.h"

#include "/Users/matthias/UserCode/mschrode/util/LabelFactory.h"
#include "/Users/matthias/UserCode/mschrode/util/StyleSettings.h"


std::vector< std::vector<double> > pars_;  // <Parameter < Variation > >
std::vector< std::vector<double> > errors_;
std::vector< std::vector<double> > corr_;
std::vector<TString> label_;
std::vector<int> color_;
std::vector<int> lineStyle_;
std::vector<int> markerStyle_;
TString outNamePrefix_;

int nPars();
int nVariations();
double par(int p, int v);
double error(int p, int v);
double corr(int p1, int p2, int v);
double sigma(double pt, int v);
double sigma(double pt, const std::vector<double> &par);
double sigmaError(double pt, int v);
void addVariation(const TString &fileName, const TString &label);
void plotSigma(int nBins, double ptMin, double ptMax, const std::vector<double> &truePar);
void plotRelDiff(int nBins, double ptMin, double ptMax, const std::vector<double> &truePar);
void plotResolution(const TString &fileName, const std::vector<double> &ptBinEdges);
void plotSpectrum(const TString &fileName, double ptMin = 1., double ptMax = 0.);
void print();


int nPars() {
  assert( pars_.size() == errors_.size() );
  return static_cast<int>(pars_.size());
}

int nVariations() {
  assert( nPars() > 0 );
  assert( pars_[0].size() == errors_[0].size() );
  return static_cast<int>(pars_[0].size());
}

double par(int p, int v) {
  assert( p >= 0 && p < nPars() );
  assert( v >= 0 && v < nVariations() );
  return pars_[p][v];
}

double error(int p, int v) {
  assert( p >= 0 && p < nPars() );
  assert( v >= 0 && v < nVariations() );
  return errors_[p][v];
}

double corr(int p1, int p2, int v) {
  int idx1 = p1;
  int idx2 = p2;
  if( p2 > p1 ) {
    idx1 = p2;
    idx2 = p1;
  }
  int idx = idx1 == 0 ? 0 : idx1*(idx1-1) + idx2 + 1;
  return corr_[idx][v];
}


double sigma(double pt, int v) {
  std::vector<double> pars(3);
  for(int p = 0; p < 3; p++) {
    pars[p] = par(p,v);
  }
  return sigma(pt,pars);
}


double sigma(double pt, const std::vector<double> &par) {
  assert(par.size() == 3);  
  return sqrt( par[0]*par[0] + par[1]*par[1]*pt + par[2]*par[2]*pt*pt );
}


double sigmaError(double pt, int v) {
  // Calculate derivatives
  std::vector<double> ds(3);
  double s = sigma(pt,v);
  for(int i = 0; i < 3; i++) {
    ds[i] = par(i,v)/s;
    if( i == 1 ) ds[i] *= pt;
    if( i == 2 ) ds[i] *= pt*pt;
  }

  //  std::cout << "Derivate ok" << std::endl;

  // Calculate variance
  double var = 0.;
  for(int i = 0; i < nPars(); i++) { // Outer loop over parameters
    for(int j = 0; j < i+1; j++) { // Inner loop over parameters
      if( i == j ) { // Diagonal terms
	var += ds[i]*ds[i]*corr(i,i,v)*error(i,v)*error(i,v);
      } else { // Off-diagonal terms
	var += 2*ds[i]*ds[j]*corr(i,j,v)*error(i,v)*error(j,v);
      }
    } // End of inner loop over parameters
  } // End of outer loop over parameters

  //  std::cout << "Variance ok" << std::endl;
  

  return sqrt(var);
}


void addVariation(const TString &fileName, const TString &label) {
  std::cout << "Adding parameter set from file " << fileName << std::endl;

  std::set<int> fixedPar;

  TFile file(fileName,"READ");
  TH1D *hPar = 0;
  TH2D *hCorr = 0;
  file.GetObject("hAbsoluteParameters",hPar);
  file.GetObject("hParameterCorrelations",hCorr);
  if( hPar == 0 ) {
    std::cerr << "ERROR: No histogram 'hAbsoluteParameters' found." << std::endl;
    exit(-1);
  } else if( hCorr == 0 ) {
    std::cerr << "ERROR: No histogram 'hParameterCorrelations' found." << std::endl;
    exit(-1);
  } else {
    hPar->SetDirectory(0);
    hCorr->SetDirectory(0);

    // Find fixed parameters
    for(int bin = 1; bin <= hPar->GetNbinsX(); bin++) {
      if( hPar->GetBinError(bin) == 0. ) {
	fixedPar.insert(bin-1);
      }
    }

    if( nPars() == 0 ) {
      int nPar = hPar->GetNbinsX()-fixedPar.size();
      pars_ = std::vector< std::vector<double> >(nPar);
      errors_ = std::vector< std::vector<double> >(nPar);
      corr_ = std::vector< std::vector<double> >(nPar*(nPar+1)/2);
    }
    int p = 0;
    int c = 0;
    for(int bin = 1; bin <= hPar->GetNbinsX(); bin++) {
      if( fixedPar.find(bin-1) == fixedPar.end() ) {
	pars_[p].push_back(std::abs(hPar->GetBinContent(bin)));
	errors_[p].push_back(hPar->GetBinError(bin));
	p++;

	for(int bin2 = 1; bin2 <= bin; bin2++) {
	  if( fixedPar.find(bin2-1) == fixedPar.end() ) {
	    corr_[c].push_back(hCorr->GetBinContent((hCorr->GetBin(bin,bin2))));
	    //std::cout << ">> " << c << " (" << bin << "," << bin2 << "): " << corr_[c].back() << std::endl;
	    c++; // What a nice joke...
	  }
	}
      }
    }      
    label_.push_back(label);
    lineStyle_.push_back(1);
    markerStyle_.push_back(label_.size()+19);
    if( label_.size() > color_.size() ) {
      color_.push_back(1);
    }
  }
}


void plotSigma(int nBins, double ptMin, double ptMax, const std::vector<double> &truePar) {
  assert( truePar.size() == 3 );
  
  // Initializing histograms
  TH1D *hRelSigma = new TH1D("hRelSigma",";p_{T} (GeV);#sigma(p_{T}) / p_{T}",nBins,ptMin,ptMax);
  hRelSigma->SetLineColor(2);
  hRelSigma->SetLineWidth(2);

  TH1D *hRelSigmaErr = static_cast<TH1D*>(hRelSigma->Clone("hRelSigmaErr"));
  hRelSigmaErr->SetLineColor(43);
  hRelSigmaErr->SetFillColor(43);
  hRelSigmaErr->SetLineWidth(1);

  TH1D *hRelSigmaTrue = static_cast<TH1D*>(hRelSigma->Clone("hRelSigmaTrue"));
  hRelSigmaTrue->SetLineColor(1);
  hRelSigmaTrue->SetLineStyle(2);


  // Filling histograms
  double min = 0.;
  double max = 0.;
  for(int bin = 1; bin <= hRelSigma->GetNbinsX(); bin++) {
    double pt = hRelSigma->GetBinCenter(bin);
    hRelSigma->SetBinContent( bin, sigma(pt,0)/pt );
    hRelSigmaErr->SetBinContent( bin, sigma(pt,0)/pt );
    hRelSigmaErr->SetBinError( bin, sigmaError(pt,0)/pt );
    hRelSigmaTrue->SetBinContent( bin, sigma(pt,truePar)/pt );

    // Store y range
    if( bin == 1 ) {
      max = hRelSigma->GetBinContent(bin);
      if( hRelSigmaTrue->GetBinContent(bin) > max ) max = hRelSigmaTrue->GetBinContent(bin);
    } else if( bin == hRelSigma->GetNbinsX() ) {
      min = hRelSigma->GetBinContent(bin);
      if( hRelSigmaTrue->GetBinContent(bin) < min ) min = hRelSigmaTrue->GetBinContent(bin);
    }      
  }
  hRelSigmaErr->GetYaxis()->SetRangeUser(0.5*min,1.5*max);

  //Draw histograms
  TCanvas *can1 = new TCanvas("canRelSigma","Relative sigma",500,500);
  can1->cd();
  hRelSigmaErr->Draw("LE3");
  hRelSigma->Draw("Lsame");
  hRelSigmaTrue->Draw("Lsame");

  TLegend *leg = util::LabelFactory::createLegend(3);
  leg->AddEntry(hRelSigmaTrue,"Generated resolution","L");
  leg->AddEntry(hRelSigma,"Fitted resolution","L");
  leg->AddEntry(hRelSigmaErr,"Statistical uncertainties","F");
  leg->Draw("same");

  if( outNamePrefix_ != "" ) can1->SaveAs(outNamePrefix_+"Sigma.eps","eps"); 
}


void plotRelDiff(int nBins, double ptMin, double ptMax, const std::vector<double> &truePar) {
  assert( truePar.size() == 3 );
  
  // Initializing histograms
  std::vector<TH1D*> hRelDiff(nVariations());
  for(int v = 0; v < nVariations(); v++) {
    TString name = "hRelDiff";
    name += v;
    hRelDiff[v] = new TH1D(name,";p_{T} (GeV);(#sigma - #sigma_{true}) / #sigma_{true}",nBins,ptMin,ptMax);
    hRelDiff[v]->SetLineColor(color_[v]);
    hRelDiff[v]->SetLineWidth(2);
  }
  TH1D *hFit = static_cast<TH1D*>(hRelDiff[0]->Clone("hFit"));
  hFit->SetLineColor(2);
  TH1D *hFitError = static_cast<TH1D*>(hRelDiff[0]->Clone("hFitError"));
  hFitError->SetLineColor(43);
  hFitError->SetFillColor(43);
  hFitError->SetLineWidth(1);
  TH1D *hTruth = static_cast<TH1D*>(hRelDiff[0]->Clone("hTruth"));
  hTruth->SetLineColor(1);
  hTruth->SetLineStyle(2);

  // Filling histograms
  double min = 100.;
  double max = 0.;
  for(int bin = 1; bin <= nBins; bin++) {
    hTruth->SetBinContent(bin,0.);
    double pt = hTruth->GetBinCenter(bin);
    double sigTrue = sigma(pt,truePar);
    double sigFit = sigma(pt,0);
    hFit->SetBinContent(bin,(sigFit-sigTrue)/sigTrue);
    hFitError->SetBinContent(bin,(sigFit-sigTrue)/sigTrue);
    hFitError->SetBinError(bin,sigmaError(pt,0)/sigTrue);
    for(int v = 0; v < nVariations(); v++) {
      double sig = sigma(pt,v);
      sig =  (sig - sigFit)/sigFit;
      hRelDiff[v]->SetBinContent(bin,sig);

      // Store y range
      if( sig > max ) max = sig;
      if( sig < min ) min = sig;
    }
  }
  hTruth->GetYaxis()->SetRangeUser(-0.1,0.2);
  hRelDiff[0]->GetYaxis()->SetRangeUser(-0.1,0.2);
  
  //Draw histograms
  TCanvas *can1 = new TCanvas("canUncert","Relative uncertainties",500,500);
  can1->cd();
  TLegend *leg1 = util::LabelFactory::createLegend(nVariations());
  leg1->AddEntry(hRelDiff[0],"Correct","L");
  hRelDiff[0]->SetLineStyle(2);
  hRelDiff[0]->Draw();
  for(int v = 1; v < nVariations(); v++) {
    hRelDiff[v]->Draw("Lsame");
    leg1->AddEntry(hRelDiff[v],label_[v],"L");
  }
  leg1->Draw("same");
  if( outNamePrefix_ != "" ) can1->SaveAs(outNamePrefix_+"SigmaUncertainties.eps","eps"); 

  TCanvas *can2 = new TCanvas("canRelDiff","Relative difference",500,500);
  can2->cd();
  hTruth->Draw();
  hFitError->Draw("LE3same");
  hTruth->Draw("same");
  hFit->Draw("Lsame");
  can2->RedrawAxis();

  TLegend *leg2 = util::LabelFactory::createLegend(3);
  leg2->AddEntry(hTruth,"Generated resolution","L");
  leg2->AddEntry(hFit,"Fitted resolution","L");
  leg2->AddEntry(hFitError,"Statistical uncertainties","F");
  leg2->Draw("same");

  if( outNamePrefix_ != "" ) can2->SaveAs(outNamePrefix_+"SigmaRelDifference.eps","eps"); 
}


void plotResolution(const TString &fileName, const std::vector<double> &ptBinEdges) {
  TFile file(fileName,"READ");
  std::vector<TH1F*> hTrueRes;
  std::vector<TH1F*> hPdfRes;
  for(size_t i = 0; i < ptBinEdges.size()-1; i++) {
    TH1F *hTruth = 0;
    TH1F *hPdf = 0;
    TString name = "hRespMeas_";
    name += i;
    file.GetObject(name,hTruth);
    if( hTruth == 0 ) {
      std::cerr << "ERROR: No histogram '" << name << "' found." << std::endl;
      exit(-1);
    }
    name = "hRespFit_";
    name += i;
    file.GetObject(name,hPdf);
    if( hPdf == 0 ) {
      std::cerr << "ERROR: No histogram '" << name << "' found." << std::endl;
      exit(-1);
    }
    hTruth->SetDirectory(0);
    hPdf->SetDirectory(0);
    hTrueRes.push_back(hTruth);
    hPdfRes.push_back(hPdf);
  }
  for(size_t i = 0; i < hTrueRes.size(); i++) {
    TString name = "Resolution_";
    name += i;
    TCanvas *can = new TCanvas(name,name,500,500);
    can->cd();
    hTrueRes[i]->Draw();
    hPdfRes[i]->Draw("Lsame");

    TLegend *leg = util::LabelFactory::createLegend(2);
    char entry[50];
    sprintf(entry,"Truth (%.0f < p^{gen}_{T} < %.0f GeV)",ptBinEdges[i],ptBinEdges[i+1]);
    leg->AddEntry(hTrueRes[i],entry,"L");
    leg->AddEntry(hPdfRes[i],"Fitted resolution","L");
    leg->Draw("same");
    
    if( outNamePrefix_ != "" ) {
      name = outNamePrefix_;
      name += "ResolutionBin";
      name += i;
      name += ".eps";
      can->SaveAs(name,"eps"); 
    }
  }
}


void plotSpectrum(const TString &fileName, double ptMin, double ptMax) {
  TFile file(fileName,"READ");
  TH1F *hPtGen = 0;
  TH1F *hPdfPtTrue = 0;
  file.GetObject("hPtGen",hPtGen);
  file.GetObject("hTruthPDF",hPdfPtTrue);
  if( hPtGen == 0 ) {
    std::cerr << "ERROR: No histogram 'hPtGen' found." << std::endl;
    exit(-1);
  } else if( hPdfPtTrue == 0 ) {
    std::cerr << "ERROR: No histogram 'hTruthPDF' found." << std::endl;
    exit(-1);
  } else {
    hPtGen->SetDirectory(0);
    hPdfPtTrue->SetDirectory(0);
  }

  TLegend *leg = 0;
  TPaveText *text = 0;
  if( ptMax > ptMin ) {
    char entry[50];
    sprintf(entry,"%.0f < p_{T} < %.0f GeV",ptMin,ptMax);
    text = util::LabelFactory::createPaveText(1,0.6);
    text->AddText(entry);
    leg = util::LabelFactory::createLegend(2,0.6,util::LabelFactory::lineHeight(),1.1*util::LabelFactory::lineHeight());
  } else {
    leg = util::LabelFactory::createLegend(2,0.6);
  }
  leg->AddEntry(hPtGen,"p^{gen}_{T} spectrum","L");
  leg->AddEntry(hPdfPtTrue,"Pdf  f(p^{true}_{T})","L");

  TCanvas *can1 = new TCanvas("canSpectrumLin","Spectrum linear",500,500);
  can1->cd();
  hPtGen->Draw();
  hPdfPtTrue->Draw("Lsame");
  if( ptMax > ptMin ) text->Draw("same");
  leg->Draw("same");
  can1->RedrawAxis();
  if( outNamePrefix_ != "" ) can1->SaveAs(outNamePrefix_+"SpectrumLinear.eps","eps"); 

  TCanvas *can2 = new TCanvas("canSpectrumLog","Spectrum log",500,500);
  can2->cd();
  hPtGen->SetMinimum(6E-8);
  hPtGen->Draw();
  hPdfPtTrue->Draw("Lsame");
  if( ptMax > ptMin ) text->Draw("same");
  leg->Draw("same");
  can2->SetLogy();
  can2->RedrawAxis();
  if( outNamePrefix_ != "" ) can2->SaveAs(outNamePrefix_+"SpectrumLog.eps","eps"); 
}



void print() {
  std::cout << "\n\nFitted Parameters\n\n";
  std::cout << "Var" << std::flush;
  for(int p = 0; p < nPars(); p++) {
    std::cout << "\t\tParameter " << p << std::flush;
  }
  std::cout << std::endl;
  for(int v = 0; v < nVariations(); v++) {
    std::cout << v << std::flush;
    for(int p = 0; p < nPars(); p++) {
      std::cout << "\t" << par(p,v) << " +/- " << error(p,v) << std::flush;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "\nCorrelations\n";
  std::cout << "Var" << std::flush;
  for(int p1 = 0; p1 < nPars(); p1++) {
    for(int p2 = 0; p2 <= p1; p2++) {
      std::cout << "\t(" << p1 << "," << p2 << ")" << std::flush;
    }
  }
  std::cout << std::endl;
  for(int v = 0; v < nVariations(); v++) {
    std::cout << v << std::flush;
    for(int p1 = 0; p1 < nPars(); p1++) {
      for(int p2 = 0; p2 <= p1; p2++) {
	std::cout << "\t" << corr(p1,p2,v) << std::flush;
      }
    }
    std::cout << std::endl;
  }
}


void plotGaussianResolution() {
  // 0: PtGen cuts
  // 1: PtCuts and variations
  int mode = 1;

  util::StyleSettings::presentationNoTitle();

  color_.push_back(1);
  color_.push_back(2);
  color_.push_back(4);
  color_.push_back(8);

  std::vector<double> truePar(3);
  truePar[0] = 4.;
  truePar[1] = 1.2;
  truePar[2] = 0.05;

  std::vector<double> ptBinEdges;
  ptBinEdges.push_back(50.);
  ptBinEdges.push_back(60.);
  ptBinEdges.push_back(80.);
  ptBinEdges.push_back(100.);
  ptBinEdges.push_back(120.);
  ptBinEdges.push_back(140.);
  ptBinEdges.push_back(170.);
  ptBinEdges.push_back(200.);
  ptBinEdges.push_back(250.);
  ptBinEdges.push_back(300.);
  ptBinEdges.push_back(400.);
  ptBinEdges.push_back(500.);
  ptBinEdges.push_back(700.);
  ptBinEdges.push_back(1000.);

  if( mode == 0 ) {
    // PtGen cuts
    outNamePrefix_ = "resFit_ToyMC_PtGenCuts_";
    addVariation("toyMC_gauss_PtGenCuts_0050-1000/jsResponse.root","PtGen cuts");
    plotSigma(500,50.,1000.,truePar);
    plotRelDiff(500,50.,1000.,truePar);
    plotSpectrum("toyMC_gauss_PtGenCuts_0050-1000/jsResponse.root");
    plotResolution("toyMC_gauss_PtGenCuts_0050-1000/jsResponse.root",ptBinEdges);
  } else if( mode == 1 ) {
    // Pt cuts
    outNamePrefix_ = "resFit_ToyMC_PtCuts_";
    addVariation("toyMC_gauss_PtCuts_0080-0800/jsResponse.root","Pt cuts");
    addVariation("toyMC_gauss_PtCuts_0080-0800_SpecDown30/jsResponse.root","Spectrum -30%");
    addVariation("toyMC_gauss_PtCuts_0080-0800_SpecUp30/jsResponse.root","Spectrum +30%");
    addVariation("toyMC_gauss_PtCuts_0080-0800_SigmaDown30/jsResponse.root","Resolution -30%");
    addVariation("toyMC_gauss_PtCuts_0080-0800_SigmaUp30/jsResponse.root","Resolution +30%");
    plotSigma(500,80.,800.,truePar);
    plotRelDiff(500,80.,800.,truePar);
    plotSpectrum("toyMC_gauss_PtCuts_0080-0800/jsResponse.root",80.,800.);
    plotResolution("toyMC_gauss_PtCuts_0080-0800/jsResponse.root",ptBinEdges);
  }
  print();
}

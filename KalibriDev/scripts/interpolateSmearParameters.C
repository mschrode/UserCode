#include <iomanip>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TString.h"
#include "TStyle.h"


std::vector< std::vector<double> > pars_;
std::vector< std::vector<double> > errors_;
std::vector<double> pt_;
std::vector<double> ptWidth_;
std::vector<double> ptBinEdges_;
std::vector<TString> fileNames_;
std::vector<TString> parNames_;



void init() { 
  //  ptBinEdges_.push_back(70.);
  //  ptBinEdges_.push_back(100.);
  //  ptBinEdges_.push_back(150.);
  ptBinEdges_.push_back(200.);
  ptBinEdges_.push_back(250.);
  ptBinEdges_.push_back(300.);
  ptBinEdges_.push_back(350.);
  ptBinEdges_.push_back(400.);
  ptBinEdges_.push_back(500.);
  ptBinEdges_.push_back(700.);
  ptBinEdges_.push_back(1000.);
  ptBinEdges_.push_back(1500.);
  ptBinEdges_.push_back(3000.);

  for(size_t i = 0; i < ptBinEdges_.size()-1; i++) {
    TString name = "~/results/JetSmearing/PtGenBins/PtGen_";
    if( ptBinEdges_[i] < 1000. ) name += "0";
    if( ptBinEdges_[i] < 100. ) name += "0";
    char val[10];
    sprintf(val,"%.0f",ptBinEdges_[i]);
    name += val;
    name += "-";
    if( ptBinEdges_[i+1] < 1000. ) name += "0";
    if( ptBinEdges_[i+1] < 100. ) name += "0";
    sprintf(val,"%.0f",ptBinEdges_[i+1]);
    name += val;
    name += "/jsResponse.root";
    fileNames_.push_back(name);
  }

  for(size_t i = 0; i < ptBinEdges_.size()-1; i++) {
    pt_.push_back( (ptBinEdges_[i]+ptBinEdges_[i+1])/2. );
    ptWidth_.push_back( (ptBinEdges_[i+1]-ptBinEdges_[i])/2. );
    parNames_.push_back("");
    parNames_.back() += i;
  }

  parNames_[0] = "#mu";
  parNames_[1] = "#sigma";
  parNames_[2] = "#alpha";
  parNames_[3] = "n";
  parNames_[4] = "m";  
}



void readParameters() {
  for(size_t f = 0; f < fileNames_.size(); f++) {
    std::cout << "Reading parameters from file " << fileNames_[f] << std::endl;
    TFile file(fileNames_[f],"READ");
    TH1D *h = 0;
    file.GetObject("hParameters",h);
    if( h == 0 ) {
      std::cerr << "ERROR: No histogram 'hAbsoluteParameters' found." << std::endl;
      exit(-1);
    } else {
      h->SetDirectory(0);
      if( f == 0 ) {
	pars_ = std::vector< std::vector<double> >(h->GetNbinsX());
	errors_ = std::vector< std::vector<double> >(h->GetNbinsX());
      }
      for(size_t i = 0; i < pars_.size(); i++) {
	pars_[i].push_back(h->GetBinContent(1+i));
	errors_[i].push_back(h->GetBinError(1+i));
	if( i == 1 ) {
	  if( f == 0 ) {
	    pars_[i].back() *= 0.1;
	    errors_[i].back() *= 0.1;
	  } else {
	    pars_[i].back() *= 0.01;
	    errors_[i].back() *= 0.01;
	  }
	}
      }      
    }
  }
}



void interpolate() {
  gStyle->SetOptFit(0011);

  std::vector<TGraphErrors*> gPars(pars_.size());
  std::vector<TCanvas*> cans(pars_.size());
  for(size_t i = 0; i < pars_.size(); i++) {
    TString name = "can";
    name += i;
    TString title = "Par ";
    title += i;
    cans[i] = new TCanvas(name,title,500,500);
    cans[i]->cd();

    gPars[i] = new TGraphErrors(pars_[i].size(),&(pt_.front()),&(pars_[i].front()),
				&(ptWidth_.front()),&(errors_[i].front()));
    gPars[i]->SetMarkerStyle(20);
    title = "Parameter ";
    title += i;
    title += ": ";
    title += parNames_[i];
    title += ";p_{T} (GeV)";
    gPars[i]->SetTitle(title);
    gPars[i]->Draw("AP");
  }

  // Fit graphs
  std::vector<TF1*> fits(pars_.size());
  fits[0] = new TF1("fitsMean","[0]",50,1500);

  fits[1] = new TF1("fitsSigma","sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/x/x)",50,1500);
  fits[1]->SetParameter(0,0.03);
  fits[1]->SetParameter(1,1.2);
  fits[1]->SetParameter(2,4.);

  fits[2] = new TF1("fitsAlpha","[0] + [1]/x/x",50,1500);
  fits[2]->SetParameter(0,2.);
  fits[2]->SetParameter(1,0.);

  fits[3] = new TF1("fitsN","[0] + [1]/x/x",50,1500);
  fits[3]->SetParameter(0,3.);
  fits[3]->SetParameter(1,0.);

  fits[4] = new TF1("fitsSpec","sqrt([0] + [1]*x)",50,1500);
  fits[4]->SetParameter(0,3.);
  fits[4]->SetParameter(1,0.005);

  for(size_t i = 0; i < pars_.size(); i++) {
    gPars[i]->Fit(fits[i],"0I");
    cans[i]->cd();
    fits[i]->Draw("same");
  }
}



void printParameters() {
  int w = 12;
  std::cout << std::setw(30) << "\n";
  for(size_t i = 0; i < pars_.size(); i++) {
    std::cout << std::setw(w) << i;
  }
  std::cout << std::endl;
  for(size_t p = 0; p < pars_[0].size(); p++) {
    std::cout << ptBinEdges_[p] << " - " << ptBinEdges_[p+1];
    for(size_t i = 0; i < pars_.size(); i++) {
      std::cout << std::setw(w) << pars_[i][p];
    }
    std::cout << std::endl;
  }
}



void interpolateSmearParameters() { 
  init();
  readParameters();
  printParameters();
  interpolate();
}


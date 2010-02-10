#include <cassert>
#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TString.h"
#include "TStyle.h"


double tMin_;
double tMax_;
TH1F *hPtGen_;


void init(double min, double max) {
  TFile f("jsResponse.root","READ");
  hPtGen_ = 0;
  f.GetObject("hPtGen",hPtGen_);
  if( !hPtGen_ ) {
    std::cerr << "ERROR: Histogram not found.\n";
    exit(-1);
  } else {
    hPtGen_->SetDirectory(0);
    tMin_ = min;
    tMax_ = max;
  }
}


double expo(double *x, double *par) {
  //  return sqrt( par[0] + par[1]*x[0] );
  return par[0] + par[1]*x[0] + par[2]/x[0];
}


double spec(double *x, double *par) {
  double f = 0.;
  if( tMin_ < x[0] && x[0] < tMax_ ) {
    double m = expo(x,par);
    f = pow(x[0],-m);
  }

  return f;
}


double norm(double *par) {
  TF1 *fn = new TF1("fn",spec,tMin_,tMax_,3);
  for(unsigned int i = 0; i < 3; i++) {
    fn->SetParameter(i,par[i]);
  }
  return fn->Integral(tMin_,tMax_);
}


double normSpec(double *x, double *par) {
  return spec(x,par)/norm(par);
}


double specBins(double *x, double *par) {
  double f = 0.;
  if( par[0] < x[0] && x[0] < par[1] ) {
    double norm = 1.;
    double m = par[3];
    if( m == 1 ) {
      norm = log(par[1]/par[0]);
    } else {
      norm = ( pow(par[1],1.-m) - pow(par[0],1.-m) )/(1.-m);
    }
    f = par[2] * pow(x[0],-m) / norm;
  }

  return f;
}


void fitPtGenSpectrum() {
  TF1 *fit = new TF1("fit",normSpec,tMin_,tMax_,3);
  fit->SetParameter(0,5.1);
  fit->SetParameter(1,3.6E-4);
  fit->SetParameter(2,58.5);
  //fit->FixParameter(2,0);
  //fit->SetParameter(3,1E-1);
  fit->SetLineColor(2);
  hPtGen_->Fit(fit,"0I");

  std::cout << "50: " << fit->Eval(50.) << std::endl;
  std::cout << "100: " << fit->Eval(100.) << std::endl;
  std::cout << "150: " << fit->Eval(150.) << std::endl;
  std::cout << "500: " << fit->Eval(500.) << std::endl;

  std::cout << "Int(50,100) = " << fit->Integral(50,100) << std::endl;
  std::cout << "Int(100,200) = " << fit->Integral(100,200) << std::endl;
  std::cout << "Int(20,1500) = " << fit->Integral(20,1500) << std::endl;
  
  TCanvas *can = new TCanvas("can","Spectrum",1200,600);
  can->Divide(2,1);
  can->cd(1);
  hPtGen_->Draw();
  fit->Draw("same");
  hPtGen_->Draw("same");
  gPad->SetLogy();

  can->cd(2);
  hPtGen_->Draw();
  fit->Draw("same");
  hPtGen_->Draw("same");
  gPad->SetLogx();
  gPad->SetLogy();
}


void fitBins(const std::vector<double> edges) {
  int nBins = static_cast<int>(edges.size()-1);
  for(int i = 0; i < nBins; i++) {
    assert( edges[i+1] > edges[i] );
  }

  std::vector<TH1F*> hSpecs(nBins);
  std::vector<TF1*> fSpecs(nBins);
  std::vector<double> ptMean(nBins);
  std::vector<double> ptMin(nBins);
  std::vector<double> ptMax(nBins);
  std::vector<double> par(nBins);
  std::vector<double> error(nBins);
  for(int i = 0; i < nBins; i++) {
    int binMin = hPtGen_->FindBin(edges[i]);
    double min = hPtGen_->GetBinLowEdge(binMin);
    int binMax = hPtGen_->FindBin(edges[i+1]);
    double max = hPtGen_->GetBinLowEdge(binMax) + hPtGen_->GetBinWidth(binMax);

    TString name = "hSpecBin";
    name += i;
    char title[100];
    sprintf(title,"%.0f < p_{T} < %.0f GeV",min,max);
    hSpecs[i] = new TH1F(name,title,1+binMax-binMin,min,max);
    for(int bin = 1; bin <= 1+binMax-binMin; bin++) {
      hSpecs[i]->SetBinContent(bin,hPtGen_->GetBinContent(binMin+bin-1));
      hSpecs[i]->SetBinError(bin,hPtGen_->GetBinError(binMin+bin-1));
    }
    double norm = hSpecs[i]->Integral("width");

    name = "fSpecBin";
    name += i;
    fSpecs[i] = new TF1(name,specBins,min,max,4);
    fSpecs[i]->SetParameter(0,min);
    fSpecs[i]->FixParameter(0,min);
    fSpecs[i]->SetParameter(1,max);
    fSpecs[i]->FixParameter(1,max);
    fSpecs[i]->SetParameter(2,norm);
    fSpecs[i]->FixParameter(2,norm);
    fSpecs[i]->SetParameter(3,5.);  
    hSpecs[i]->Fit(fSpecs[i],"0QRI");
    fSpecs[i]->SetLineWidth(1);

    ptMean[i] = hSpecs[i]->GetMean();
    ptMin[i] = ptMean[i] - min;
    ptMax[i] = max - ptMean[i];
    par[i] = fSpecs[i]->GetParameter(3);
    error[i] = fSpecs[i]->GetParError(3);
    std::cout << i << ": " << par[i] << " +/- " << error[i] << std::endl;
  }

  TCanvas *can = new TCanvas("canFitBins","Specs",1200,800);
  int div = 1;
  while( 6*div*div < nBins ) div++;
  can->Divide(3*div,2*div);
  for(int i = 0; i < nBins; i++) {
    can->cd(1+i);
    hSpecs[i]->Draw();
    fSpecs[i]->Draw("same");
  }

  TCanvas *can2 = new TCanvas("canFittedPars","Parameter",600,600);
  can2->cd();
  TGraphAsymmErrors *graph
    = new TGraphAsymmErrors(nBins,&(ptMean.front()),&(par.front()),
			    &(ptMin.front()),&(ptMax.front()),
			    &(error.front()),&(error.front()));
  graph->SetMarkerStyle(20);
  graph->SetTitle("Exponent of spectrum;p_{T} (GeV)");
  graph->Draw("APE1");
}


void fitSpectrum(double min = 20, double max = 1500) {
  init(min,max);

  fitPtGenSpectrum();

  std::vector<double> edges;
  edges.push_back(20.);
  edges.push_back(30.);
  edges.push_back(40.);
  edges.push_back(50.);
  edges.push_back(70.);
  edges.push_back(100.);
  edges.push_back(150.);
  edges.push_back(200.);
  edges.push_back(250.);
  edges.push_back(300.);
  edges.push_back(350.);
  edges.push_back(400.);
  edges.push_back(500.);
  edges.push_back(600.);
  edges.push_back(700.);
  edges.push_back(800.);
  edges.push_back(900.);
  edges.push_back(1000.);
  edges.push_back(1100.);
  edges.push_back(1200.);
  edges.push_back(1500.);

  //fitBins(edges);
}

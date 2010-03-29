// $Id: $

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"



//! Investigate influence of a 3rd jet on fitted jet resolution
//! and find correction function
// ------------------------------------------------------------
namepace correct3rdJetPt {

  // ----- Global parameters -----
  const int nEvts_ = 100000;
  const double tMin_ = 50.;
  const double tMax_ = 500.;
  const double dt_ = 0.9;
  const double b_[3] = { 4., 1.2, 0.05 };
  const int nBins_ = 50;


  // ----- Global variables -----
  TH1 *hBinning_ = 0;


  // ----- Type defintions and classes -----
  class Dijet {
  public:
    Dijet(double t1, double x1, double t2, double x2) : t1_(t1), t2_(t2), x1_(x1), x2_(x2) {};

    double x1() const { return x1_; }
    double x2() const { return x2_; } 
    double t1() const { return t1_; }
    double t2() const { return t2_; }

  private:
    const double t1_, t2_, x1_, x2_;
  };
  typedef std::vector<Dijet*> Data;
  typedef std::vector<Dijet*>::const_iterator DataIt;

  typedef std::vector<TH1*> Hists;
  typedef std::vector<TF1*> Functs;


  // ----- Function definitions -----

  //! Return width of Gaussian resolution
  double sigma(double t) {
    return sqrt( b_[0]*b_[0] + b_[1]*b_[1]*t + b_[2]*b_[2]*t*t );
  }

  //! Generate \p nEvts_ \p Dijet events
  Data generateData() {
    std::cout << "Generating " << nEvts_ << " dijet events... " << std::flush;
    Data data(nEvts_);
  
    TRandom3 rand(0);
    double t[2] = { 0., 0. };
    double x[2] = { 0., 0. };
    for(size_t i = 0; i < data.size(); ++i) {
      t[0] = rand.Uniform(tMin_,tMax_);
      t[1] = dt_*t[0];
      for(int j = 0; j < 2; ++j) {
	x[j] = rand.Gaus(t[j],sigma(t[j]));
      }
      data[i] = new Dijet(t[0],x[0],t[1],x[1]);
    }

    std::cout << "ok" << std::endl;
    return data;
  }

  //! Find the bin 0 < bin < \p nBins_ for a given true pt
  int bin(double t) {
    int b = static_cast<int>(hBinning_->GetBinContent(hBinning_->FindBin(t)));
    if( b < 0 ) b = 0;
    else if( b >= nBins_ ) b = nBins_-1;
    return b;
  }

  //! Fit different defintions of the resolution and 
  //! find a correction function
  void fitResolution(const Data &data, bool plotDistributions = true) {
    std::cout << "Fitting resolution... " << std::flush;
    // Spectrum
    TH1 *hSpectrum = new TH1D("hSpectrum",";p_{T} (GeV)",nBins_,0.,tMax_);
    // Resolution histograms in bins of t
    Hists hResMC(nBins_);
    Hists hResMeanT(nBins_);
    Functs fResMC(nBins_);
    Functs fResMeanT(nBins_);
    double dt = (tMax_-tMin_)/nBins_;
    for(int i = 0; i < nBins_; ++i) {
      TString title = "";
      title += tMin_+i*dt;
      title += " < p_{T} < ";
      title += tMin_+(1+i)*dt;
      title += ";Resolution";

      TString name = "hResMC";
      name += i;
      hResMC[i] = new TH1D(name,title,51,0.,2.);
      name = "fResMC";
      fResMC[i] = new TF1(name,"gaus",0.,2.);

      name = "hResMeanT";
      name += i;
      hResMeanT[i] = static_cast<TH1D*>(hResMC[i]->Clone(name));
      hResMeanT[i]->SetLineColor(2);
      name = "fResMeanT";
      fResMeanT[i] = static_cast<TF1*>(fResMC[i]->Clone(name));
    }
    // Filling with resolution
    for(DataIt it = data.begin(); it != data.end(); ++it) {
      const Dijet *dijet = *it;
      // Spectrum
      hSpectrum->Fill(dijet->t1());
      hSpectrum->Fill(dijet->t2());
      // MC truth resolution
      hResMC[bin(dijet->t1())]->Fill(dijet->x1()/dijet->t1());
      hResMC[bin(dijet->t2())]->Fill(dijet->x2()/dijet->t2());
      // Fitted resolution around mean t
      double meanT = 0.5*(dijet->t1()+dijet->t2());
      hResMeanT[bin(meanT)]->Fill(dijet->x1()/meanT);
      hResMeanT[bin(meanT)]->Fill(dijet->x2()/meanT);
    }

    // Fits and histograms of Gaussian width per bin
    TH1 *hSigmaMC = new TH1D("hSigmaMC",";p_{T} (GeV);#sigma(p_{T}) / p_{T}",nBins_,tMin_,tMax_);
    hSigmaMC->Sumw2();
    hSigmaMC->SetMarkerStyle(7);
    TH1 *hSigmaRelDiff = static_cast<TH1D*>(hSigmaMC->Clone("hSigmaRelDiff"));
    hSigmaRelDiff->SetTitle(";p_{T} (GeV);#sigma / #sigma_{true})");

    TH1 *hSigmaMeanT = static_cast<TH1D*>(hSigmaMC->Clone("hSigmaMeanT"));
    hSigmaMeanT->SetMarkerColor(2);
    hSigmaMeanT->SetLineColor(2);
    TH1 *hSigmaRelDiffMeanT = static_cast<TH1D*>(hSigmaMeanT->Clone("hSigmaRelDiffMeanT"));

    // Fitting Gaussian width per bin
    for(int i = 0; i < nBins_; ++i) {
      hResMC[i]->Fit(fResMC[i],"IL0Q");
      hSigmaMC->SetBinContent(1+i,fResMC[i]->GetParameter(2));
      hSigmaMC->SetBinError(1+i,fResMC[i]->GetParError(2));
      hSigmaRelDiff->SetBinContent(1+i,1.);

      hResMeanT[i]->Fit(fResMeanT[i],"IL0Q");
      hSigmaMeanT->SetBinContent(1+i,fResMeanT[i]->GetParameter(2));
      hSigmaMeanT->SetBinError(1+i,fResMeanT[i]->GetParError(2));
      hSigmaRelDiffMeanT->Divide(hSigmaMeanT,hSigmaMC);
    }

    // Fitting sigma(pt)
    std::vector<double> startPar(3);
    startPar[0] = 4.;
    startPar[1] = 1.2;
    startPar[2] = 0.05;
    TF1 *fSigmaMC = new TF1("fSigmaMC","sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",tMin_,tMax_);
    fSigmaMC->SetLineWidth(1);
    fSigmaMC->SetLineColor(hSigmaMC->GetLineColor());
    for(int i = 0; i < 3; ++i) fSigmaMC->SetParameter(i,startPar[i]);
    hSigmaMC->Fit(fSigmaMC,"R0Q");
    TF1 *fSigmaMeanT = new TF1("fSigmaMeanT","sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",tMin_,tMax_);
    fSigmaMeanT->SetLineWidth(1);
    fSigmaMeanT->SetLineColor(hSigmaMeanT->GetLineColor());
    for(int i = 0; i < 3; ++i) fSigmaMeanT->SetParameter(i,startPar[i]);
    hSigmaMeanT->Fit(fSigmaMeanT,"R0Q");

    // Plot resolution per bin
    if( plotDistributions ) {
      int nCans = 0;
      int nPadsX = 5;
      int nPadsY = 4;
      int nPads = nPadsX*nPadsY;
      while( nCans*nPads < nBins_ ) {
	TString name = "Resolution";
	name += nCans;
	TCanvas *can = new TCanvas(name,name,nPadsX*200,nPadsY*200);
	can->Divide(nPadsX,nPadsY);
	for(int p = 1; p <= nPads; ++p) {
	  can->cd(p);
	  int idx = nCans*nPads + p - 1;
	  if( idx < nBins_ ) {
	    hResMC[idx]->Draw();
	    hResMeanT[idx]->Draw("same");
	  }
	}
	nCans++;
      }
    }

    // Plot mean spectrum
    TCanvas *can1 = new TCanvas("canSpectrum","Spectrum",600,600);
    can1->cd();
    hSpectrum->Draw();

    // Plot mean sigmas
    TCanvas *can2 = new TCanvas("canMeanSigmas","Sigma",600,600);
    can2->cd();
    hSigmaMC->Draw("PE1");
    fSigmaMC->Draw("same");
    hSigmaMeanT->Draw("PE1same");
    fSigmaMeanT->Draw("same");

    // Plot rel difference
    TCanvas *can3 = new TCanvas("canRelDiffs","Differences",600,600);
    can3->cd();
    hSigmaRelDiff->SetLineStyle(2);
    hSigmaRelDiff->GetYaxis()->SetRangeUser(0.95,1.45);
    hSigmaRelDiff->Draw();
    hSigmaRelDiffMeanT->Draw("PE1same");
  
    std::cout << "ok" << std::endl;

    // Print fit results
    std::cout << std::endl;
    std::cout << setw(12) << "MC truth\t" << std::flush;
    for(int i = 0; i < 3; ++i) {
      std::cout << fSigmaMC->GetParameter(i) << " +/- " << fSigmaMC->GetParError(i) << "\t";
    }
    std::cout << std::endl;
    std::cout << setw(12) << "MC truth\t" << std::flush;
    for(int i = 0; i < 3; ++i) {
      std::cout << fSigmaMeanT->GetParameter(i) << " +/- " << fSigmaMeanT->GetParError(i) << "\t";
    }
    std::cout << std::endl;
  }

  //! Init global variables
  void init() {
    gStyle->SetErrorX(0);
    gStyle->SetOptStat(0);
    hBinning_ = new TH1I("hBinning_","",nBins_,tMin_,tMax_);
    for(int i = 0; i < nBins_; ++i) {
      hBinning_->SetBinContent(1+i,i);
    }
  }

  //! Run the fit
  void run() {
    init();
    Data data = generateData();
    fitResolution(data,false); 
  }
}

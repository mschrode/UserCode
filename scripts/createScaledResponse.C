#include <cassert>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"

#include "../util/utils.h"
#include "../util/FileOps.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"

#include "globalFunctions.h"



const double nSigCore_ = 2.;

TH1* hScalingFactors_ = 0;
TH1* hScalingFactorsUp_ = 0;
TH1* hScalingFactorsDown_ = 0;



std::vector<double> getSmearFactors(double pt) {
  std::vector<double> v(1,1.1);
  return v;
}



double getTailScalingFactor(const TH1* hFactors, double pt) {
  double min = hFactors->GetXaxis()->GetBinLowEdge(1);
  double max = hFactors->GetXaxis()->GetBinUpEdge(hFactors->GetNbinsX());
  double fac = hFactors->GetBinContent(1);
  if( pt > min && pt < max ) {
    fac = hFactors->GetBinContent(hFactors->FindBin(pt));
  } else if( pt >= max ) {
    fac = hFactors->GetBinContent(hFactors->GetNbinsX());
  }

  return fac;
}



double getTailScalingFactor(double pt) {
  return getTailScalingFactor(hScalingFactors_,pt);
}



void createScaledResponse() {

  // Get smearing and tail scaling factors
  std::cout << "Getting smearing and tail scaling factors\n";
  hScalingFactors_ = readTH1("","","hScalingFactors_");
  hScalingFactorsUp_ = readTH1("","","hScalingFactorsUp_");
  hScalingFactorsDown_ = readTH1("","","hScalingFactorsDown_");


  // Get MC truth response histograms
  std::cout << "Getting MC truth response histograms\n";
  std::vector<double> ptGenBinEdges;
  //  10 20 30 40 50 60 70 80 90 100 120 150 170 200 250 300 350 400 500 1000 1500
  util::HistVec hResMC = util::FileOps::readHistVec("","hRespMeas_");


  // Sanity checks
  unsigned int nPtBins = ptGenBinEdges.size()-1;

  assert( hResMC.size() >= nPtBins );

  std::vector<double> ptGenBinCenters(nPtBins);
  for(unsigned int i = 0; i < nPtBins; ++i) {
    ptGenBinCenters[i] = 0.5*(ptGenBinEdges[i]+ptGenBinEdges[1+i]);
  }


  // Smear MC truth response
  std::cout << "Smearing MC truth response\n";
  util::HistVec hResSmeared(nPtBins);
  for(unsigned int i = 0; i < nPtBins; ++i) {
    hResMC[i]->GetXaxis()->SetRangeUser(0.,2.);
    double width = 0.;
    double widthErr = 1000.;
    if( func::fitCoreWidth(hResMC[i],nSigCore_,width,widthErr) ) {
      func::smearHistogram(hResMC[i],
			   hResSmeared[i],
			   hResMC[i]->GetEntries(),
			   width,
			   getSmearFactors(ptBinCenters[i]));
    } else {
      std::cerr << "ERROR: Fit of Gaussian core failed." << std::endl;
      exit(1);
    }
  }
  
  
  // Dividing MC truth into core and tails
  std::cout << "Dividing MC truth into core and tails\n";
  util::HistVec hCore(nPtBins);
  util::HistVec hTails(nPtBins);
  util::HistVec hTailsClean(nPtBins);
  std::vector<TF1*> fGauss(nPtBins);
  for(unsigned int i = 0; i < nPtBins; ++i) {
    func::getTail(hResSmeared[i],nSigTail_,hTails[i],hTailsClean[i],fGauss[i]);
    hCore[i] = static_cast<TH1D*>(hResSmeared[i]->Clone("hCore"+util::toTString(i)));
    hCore[i]->Add(hTailsClean[i],-1.);
    hCore[i]->SetLineWidth(2);
  }

  
  // Scaled response
  std::cout << "Scaled response\n";
  util::HistVec hResScaled(nPtBins);
  util::HistVec hTailsScaled(nPtBins);
  for(unsigned int i = 0; i < nPtBins; ++i) {
    hTailsScaled[i] = static_cast<TH1D*>(hTailsClean[i]->Clone("hTailsScaled"+util::toTString(i)));
    hTailsScaled[i]->Scale(getTailScalingFactor(ptGenBinCenters[i]));
    hTailsClean[i]->SetLineStyle(2);
    hTailsScaled[i]->SetLineColor(2);
    hResScaled[i] = static_cast<TH1D*>(hResSmeared[i]->Clone("hResScaled"+util::toTString(i)));
    hResScaled[i]->Reset();
    for(int bin = 1; bin <= hResScaled[i]->GetNbinsX(); ++bin) {
      hResScaled[i]->SetBinContent(bin,hCore[i]->GetBinContent(bin)+hTailsScaled[i]->GetBinContent(bin));
    }
    hResScaled[i]->SetLineColor(2);
  }


   // Plots
   for(unsigned int i = 0; i < nPtBins; ++i) {
     TCanvas* canSmear = new TCanvas("canSmear"+util::toTString(i),"Smeared "+util::toTString(i),500,500);
     canSmear->cd();
     hResSmeared[i]->GetXaxis()->UnZoom();
     hResSmeared[i]->Draw("HIST");
     fGauss[i]->Draw("same");
     hResMC[i]->SetLineStyle(2);
     hResMC[i]->Draw("HISTsame");
     canSmear->SetLogy();

     TCanvas* canParts = new TCanvas("canParts"+util::toTString(i),"Parts "+util::toTString(i),500,500);
     canParts->cd();
     hCore[i]->Draw("HIST");
     hTailsScaled[i]->Draw("HISTsame");
     hTailsClean[i]->Draw("HISTsame");
     canParts->SetLogy();

     TCanvas* canScaled = new TCanvas("canScaled"+util::toTString(i),"Scaled Resp "+util::toTString(i),500,500);
     canScaled->cd();
     hResSmeared[i]->Draw("HIST");
     hResScaled[i]->Draw("HISTsame");
     canScaled->SetLogy();
   }
}

#include <cassert>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TPaveText.h"

#include "../util/utils.h"
#include "../util/FileOps.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"

#include "globalFunctions.h"



const double nSigCore_ = 2.;
const double nSigTail_ = 3.;
const double min_ = 3E-5;
const double max_ = 3E3;

TH1* hScalingFactors_ = 0;
TH1* hScalingFactorsUp_ = 0;
TH1* hScalingFactorsDown_ = 0;



double getSmearFactors(double pt) {
  return 0.;
}



double getTailScalingFactor(TH1* hFactors, double pt) {
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
  util::StyleSettings::presentationNoTitle();

  // Get smearing and tail scaling factors
  std::cout << "Getting smearing and tail scaling factors\n";
  TString fileScalingFactors = "Tails_PF_.root";
  hScalingFactors_ = util::FileOps::readTH1(fileScalingFactors,"hScaleFactorsInt_Eta0","hScalingFactors_");
  hScalingFactorsUp_ = util::FileOps::readTH1(fileScalingFactors,"hScaleFactorsIntUp_Eta0","hScalingFactorsUp_");
  hScalingFactorsDown_ = util::FileOps::readTH1(fileScalingFactors,"hScaleFactorsIntDown_Eta0","hScalingFactorsDown_");


  // Get MC truth response histograms
  std::cout << "Getting MC truth response histograms\n";
  std::vector<double> ptGenBinEdges;
  ptGenBinEdges.push_back(10.);
  ptGenBinEdges.push_back(20.);
  ptGenBinEdges.push_back(30.);
  ptGenBinEdges.push_back(40.);
  ptGenBinEdges.push_back(50.);
  ptGenBinEdges.push_back(60.);
  ptGenBinEdges.push_back(70.);
  ptGenBinEdges.push_back(80.);
  ptGenBinEdges.push_back(90.);
  ptGenBinEdges.push_back(100.);
  ptGenBinEdges.push_back(120.);
  ptGenBinEdges.push_back(150.);
  ptGenBinEdges.push_back(170.);
  ptGenBinEdges.push_back(200.);
  ptGenBinEdges.push_back(250.);
  ptGenBinEdges.push_back(300.);
  ptGenBinEdges.push_back(350.);
  ptGenBinEdges.push_back(400.);
  ptGenBinEdges.push_back(500.);
  ptGenBinEdges.push_back(1000.);
  ptGenBinEdges.push_back(1500.);

  //  10 20 30 40 50 60 70 80 90 100 120 150 170 200 250 300 350 400 500 1000 1500
  util::HistVec hResMC = util::FileOps::readHistVec("MCFall10_TruthResolution_PF.root","hRespMeas_");


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
    hResMC[i]->GetYaxis()->SetRangeUser(min_,max_);
    util::HistOps::setAxisTitles(hResMC[i],"Response","","jets",true);
    //    hResMC[i]->SetTitle("");    
    double width = 0.;
    double widthErr = 1000.;
    if( func::fitCoreWidth(hResMC[i],nSigCore_,width,widthErr) ) {
      double ratioDataMC = getSmearFactors(ptGenBinCenters[i]);
      std::cout << "  PtBin " << i << ": Data / MC (Core Width) = " << 1+ratioDataMC << std::endl;
      func::smearHistogram(hResMC[i],
			   hResSmeared[i],
			   hResMC[i]->GetEntries(),
			   width,
			   ratioDataMC);
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
    func::getTail(hResSmeared[i],nSigCore_,nSigTail_,hTails[i],hTailsClean[i],fGauss[i]);
    hCore[i] = static_cast<TH1D*>(hResSmeared[i]->Clone("hCore"+util::toTString(i)));
    hCore[i]->Add(hTails[i],-1.);
    hCore[i]->SetLineWidth(2);
  }

  
  // Scaled response
  std::cout << "Scaled response\n";
  util::HistVec hResScaled(nPtBins);
  util::HistVec hTailsScaled(nPtBins);
  for(unsigned int i = 0; i < nPtBins; ++i) {
    double tailScale = getTailScalingFactor(ptGenBinCenters[i]);
    std::cout << "  PtBin " << i << ": Tail Scaling Factor " << tailScale << std::endl;
    hTailsScaled[i] = static_cast<TH1D*>(hTails[i]->Clone("hTailsScaled"+util::toTString(i)));
    hTailsScaled[i]->Scale(tailScale);
    hTails[i]->SetLineStyle(2);
    hTailsScaled[i]->SetLineColor(2);
    hResScaled[i] = static_cast<TH1D*>(hResSmeared[i]->Clone("hResScaled"+util::toTString(i)));
    hResScaled[i]->Reset();
    for(int bin = 1; bin <= hResScaled[i]->GetNbinsX(); ++bin) {
      hResScaled[i]->SetBinContent(bin,hCore[i]->GetBinContent(bin)+hTailsScaled[i]->GetBinContent(bin));
    }
    hResScaled[i]->SetLineColor(2);
  }


   // Plots
   for(unsigned int i = 10; i < nPtBins; ++i) {
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
     hCore[i]->GetYaxis()->SetRangeUser(min_,max_);
     hCore[i]->Draw("HIST");
     hTailsScaled[i]->Draw("HISTsame");
     hTails[i]->Draw("HISTsame");
     canParts->SetLogy();

     TCanvas* canScaled = new TCanvas("canScaled"+util::toTString(i),"Scaled Resp "+util::toTString(i),500,500);
     canScaled->cd();
     hResSmeared[i]->GetYaxis()->SetRangeUser(min_,max_);
     hResSmeared[i]->Draw("HIST");
     hResScaled[i]->Draw("HISTsame");
     canScaled->SetLogy();
   }
}

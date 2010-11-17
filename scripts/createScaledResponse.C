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



std::vector<double> getSmearFactors() {
  std::vector<double> v(1,1.1);
  return v;
}


bool readScalingFactors(const TString &fileName, std::vector<double> &scalingFactors) {
  bool result = false;
  scalingFactors.clear();
  TH1* facs = util::FileOps::readTH1(fileName,"Ratio_hNTailData0");
  if( facs ) {
    scalingFactors = std::vector<double>(facs->GetNbinsX());
    for(unsigned int i = 0; i < scalingFactors.size(); ++i) {
      scalingFactors[i] = 2.;//facs->GetBinContent(1+i);
    }
    result = true;
  } else {
    std::cerr << "ERROR reading scaling factors" << std::endl;
  }
  
  return result;
}


void createScaledResponse() {
  // Get smearing and tail scaling factors
  std::cout << "Getting smearing and tail scaling factors\n";
  std::vector<double> scalingFactors;
  readScalingFactors("results/Tails_Calo_Eta00-11_Pp10.root",scalingFactors);
  std::vector<double> smearFactors = getSmearFactors();

  // Get MC truth response histograms
  std::cout << "Getting MC truth response histograms\n";
  util::HistVec hResMC = util::FileOps::readHistVec("Fall10QCDFlat_MCTruthReso.root","hRespMeas_");

  // Sanity checks
  unsigned int nPtBins = 1;//scalingFactors.size();
  assert( hResMC.size() >= nPtBins );
  assert( smearFactors.size() >= nPtBins );

  // Smear MC truth response
  std::cout << "Smearing MC truth response\n";
  util::HistVec hResSmeared(nPtBins);
  for(unsigned int i = 0; i < nPtBins; ++i) {
    hResMC[i]->GetXaxis()->SetRangeUser(0.,2.);
    double width = 0.;
    double widthErr = 1000.;
    func::fitCoreWidth(hResMC[i],width,widthErr);
    func::smearHistogram(hResMC[i],hResSmeared[i],hResMC[i]->GetEntries(),width,smearFactors[i]);
  }

   // Dividing MC truth into core and tails
   std::cout << "Dividing MC truth into core and tails\n";
   util::HistVec hCore(nPtBins);
   util::HistVec hTails(nPtBins);
   util::HistVec hTailsClean(nPtBins);
   std::vector<TF1*> fGauss(nPtBins);
   for(unsigned int i = 0; i < nPtBins; ++i) {
     func::getTail(hResSmeared[i],hTails[i],hTailsClean[i],fGauss[i]);
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
     hTailsScaled[i]->Scale(scalingFactors[i]);
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

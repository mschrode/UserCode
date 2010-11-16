#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"

#include "../util/utils.h"
#include "../util/FileOps.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"



void plotTailsVsPP3() {
  util::StyleSettings::presentationNoTitle();

  std::cout << "Setting up script" << std::endl;

  // File name and pp3 lists
  std::vector<TString> fileNames;
  fileNames.push_back("results/Tails_Calo_Eta00-11_Pp06.root");
  fileNames.push_back("results/Tails_Calo_Eta00-11_Pp10.root");
  fileNames.push_back("results/Tails_Calo_Eta00-11_Pp15.root");

  std::vector<double> pp3Limits;
  pp3Limits.push_back(0.06);
  pp3Limits.push_back(0.10);
  pp3Limits.push_back(0.15);

  assert( fileNames.size() == pp3Limits.size() );

  // Scaling factors
  util::HistVec hScales = util::FileOps::readTH1(fileNames,"Ratio_hNTailData0","ScaleFactors");
  TCanvas* cScales = new TCanvas("cScales","Scaling factors",500,500);
  cScales->cd();
  TLegend *legScales = util::LabelFactory::createLegendCol(fileNames.size(),0.4);
  util::HistOps::setYRange(hScales[0],fileNames.size());
  for(unsigned int i = 0; i < hScales.size(); ++i) {
    util::HistOps::setStyleColor(hScales[i],i);
    legScales->AddEntry(hScales[i],"p_{||} < "+util::toTString(pp3Limits[i]),"P");
    if( i == 0 ) {
      hScales[i]->GetXaxis()->SetMoreLogLabels();
      hScales[i]->Draw("PE1");
    } else {
      hScales[i]->Draw("PE1same");
    }
  }
  legScales->Draw("same");
  cScales->SetLogx();
  cScales->SaveAs("Tails_ScalingFactors.eps","eps");


  std::cout << "Plotting scaling factors" << std::endl;

  // Scaling factors vs pp3
  for(int ptBin = 1; ptBin <= hScales[0]->GetNbinsX(); ++ptBin) {
    std::vector<double> pp3(fileNames.size());
    std::vector<double> pp3Err(fileNames.size(),0.);
    std::vector<double> scale(fileNames.size());
    std::vector<double> scaleErr(fileNames.size());
    for(unsigned int pp3Bin = 0; pp3Bin < fileNames.size(); ++pp3Bin) {
      pp3[pp3Bin] = pp3Limits[pp3Bin];
      scale[pp3Bin] = hScales[pp3Bin]->GetBinContent(ptBin);
      scaleErr[pp3Bin] = hScales[pp3Bin]->GetBinError(ptBin);
    }
    TGraphErrors *gExtra = new TGraphErrors(pp3.size(),&(pp3.front()),&(scale.front()),&(pp3Err.front()),&(scaleErr.front()));
    gExtra->SetMarkerStyle(20);

    double max = *(std::max_element(gExtra->GetY(),gExtra->GetY()+gExtra->GetN()));
    TH1 *hExtraFrame = util::HistOps::createTH1D("hExtraFrame"+util::toTString(ptBin),fileNames.size(),0.,1.5*pp3Limits.back(),"p_{||} Threshold","","Scaling Factor");
    hExtraFrame->SetNdivisions(505);
    hExtraFrame->GetYaxis()->SetRangeUser(0.,2.5*max);

    TCanvas* cExtra = new TCanvas("cExtra"+util::toTString(ptBin),"Extrapolation "+util::toTString(ptBin),500,500);
    cExtra->cd();
    hExtraFrame->Draw();
    gExtra->Draw("PE1same");

    double ptMin = hScales[0]->GetXaxis()->GetBinLowEdge(ptBin);
    double ptMax = hScales[0]->GetXaxis()->GetBinUpEdge(ptBin);
    TPaveText *label = util::LabelFactory::createPaveText(1);
    label->AddText(util::toTString(ptMin)+" < p^{ave}_{T} < "+util::toTString(ptMax)+" GeV");
    label->Draw("same");
    cExtra->SaveAs("Tails_ScalingFactorsVsPP3_"+util::toTString(ptBin-1)+".eps","eps");
  }
}

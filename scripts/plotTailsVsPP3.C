#include <algorithm>
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

bool gPP3Bins_ = false;

std::vector<TString> fileNames(const TString &baseName, const std::vector<double> &ppLimits) {
  std::vector<TString> names( gPP3Bins_ ? ppLimits.size()-1 : ppLimits.size() );
  for(unsigned int i = 0; i < names.size(); ++i) {
    if( gPP3Bins_ ) {
      TString thresMin = util::toTString(100.*ppLimits[i]);
      TString thresMax = util::toTString(100.*ppLimits[i+1]);
      while( thresMin.Length() < 2 ) thresMin = "0"+thresMin;
      while( thresMax.Length() < 2 ) thresMax = "0"+thresMax;
      names[i] = baseName+"_Pp"+thresMin+"-"+thresMax+".root";
    } else {
      TString thres = util::toTString(100.*ppLimits[i]);
      while( thres.Length() < 2 ) thres = "0"+thres;
      names[i] = baseName+"_Pp"+thres+".root";
    }
  }

  return names;
}



void plotScalesVsGauss() {
  util::StyleSettings::presentationNoTitle();

  std::cout << "Setting up script" << std::endl;

  // File name and pp3 lists
  std::vector<TString> fileNames;
//   fileNames.push_back("results/Tails_Calo_Eta00-11_Pp10_Gauss12.root");
//   fileNames.push_back("results/Tails_Calo_Eta00-11_Pp10.root");
//   fileNames.push_back("results/Tails_Calo_Eta00-11_Pp10_Gauss22.root");

  std::vector<double> coreLimits;
  coreLimits.push_back(1.4);
  coreLimits.push_back(1.8);
  coreLimits.push_back(2.2);

  assert( fileNames.size() == coreLimits.size() );

  // Scaling factors
  util::HistVec hScales = util::FileOps::readTH1(fileNames,"Ratio_hNTailData0","ScaleFactors");
  TCanvas* cScales = new TCanvas("cScales","Scaling factors",500,500);
  cScales->cd();
  TLegend *legScales = util::LabelFactory::createLegendCol(fileNames.size(),0.6);
  util::HistOps::setYRange(hScales[0],fileNames.size());
  for(unsigned int i = 0; i < hScales.size(); ++i) {
    util::HistOps::setStyleColor(hScales[i],i);
    legScales->AddEntry(hScales[i],"Core < "+util::toTString(coreLimits[i])+" #sigma","P");
    if( i == 0 ) {
      hScales[i]->GetYaxis()->SetRangeUser(0.,6.);
      hScales[i]->GetXaxis()->SetMoreLogLabels();
      hScales[i]->SetYTitle("Scaling Factor");
      hScales[i]->Draw("PE1");
    } else {
      hScales[i]->Draw("PE1same");
    }
  }
  legScales->Draw("same");
  cScales->SetLogx();
  cScales->SaveAs("Tails_ScalingFactorsVsCoreRegion.eps","eps");
}



void plotTailsVsPP3() {
  util::StyleSettings::presentationNoTitle();

  std::cout << "Setting up script" << std::endl;

  std::vector<double> pp3Limits;
//   pp3Limits.push_back(0.06);
//   pp3Limits.push_back(0.10);
//   pp3Limits.push_back(0.15);

  pp3Limits.push_back(2.5);
  pp3Limits.push_back(3.);
  pp3Limits.push_back(3.5);
  pp3Limits.push_back(4.);



  //  std::vector<TString> names = fileNames("results/Tails_Calo_Eta00-11",pp3Limits);
  std::vector<TString> names;
  names.push_back("results/Tails_Calo_Eta00-11_PSoft10_NSigTail25.root");
  names.push_back("results/Tails_Calo_Eta00-11_PSoft10_NSigTail30.root");
  names.push_back("results/Tails_Calo_Eta00-11_PSoft10_NSigTail35.root");
  names.push_back("results/Tails_Calo_Eta00-11_PSoft10_NSigTail40.root");


  // Scaling factors
  util::HistVec hScales = util::FileOps::readTH1(names,"hScalFacTrend","ScaleFactors");
  //  util::HistVec hScales = util::FileOps::readTH1(names,"Ratio_hNTailData0","ScaleFactors");
  TCanvas* cScales = new TCanvas("cScales","Scaling factors",500,500);
  cScales->cd();
  TLegend *legScales = util::LabelFactory::createLegendCol(names.size(),0.4);
  util::HistOps::setYRange(hScales[0],names.size());
  for(unsigned int i = 0; i < hScales.size(); ++i) {
    util::HistOps::setStyleColor(hScales[i],i);
    if( gPP3Bins_ ) {
      legScales->AddEntry(hScales[i],util::toTString(pp3Limits[i])+" < p^{rel}_{||} < "+util::toTString(pp3Limits[i+1]),"P");
    } else {
      //legScales->AddEntry(hScales[i],"p^{rel}_{||} < "+util::toTString(pp3Limits[i]),"P");
      legScales->AddEntry(hScales[i],"Tail: "+util::toTString(pp3Limits[i])+" #sigma","P");
    }
    if( i == 0 ) {
      hScales[i]->GetXaxis()->SetMoreLogLabels();
      hScales[i]->Draw("PE1");
    } else {
      hScales[i]->Draw("PE1same");
    }
  }
  legScales->Draw("same");
  cScales->SetLogx();
  //cScales->SaveAs("Tails_ScalingFactors.eps","eps");
  cScales->SaveAs("Tails_ScalingFactors_Calo.eps","eps");


  std::cout << "Plotting scaling factors" << std::endl;

  // Scaling factors vs pp3
  for(int ptBin = 1; ptBin <= hScales[0]->GetNbinsX(); ++ptBin) {
    std::vector<double> pp3(names.size());
    std::vector<double> pp3Err(names.size(),0.);
    std::vector<double> scale(names.size());
    std::vector<double> scaleErr(names.size());
    for(unsigned int pp3Bin = 0; pp3Bin < names.size(); ++pp3Bin) {
      pp3[pp3Bin] = pp3Limits[pp3Bin];
      scale[pp3Bin] = hScales[pp3Bin]->GetBinContent(ptBin);
      scaleErr[pp3Bin] = hScales[pp3Bin]->GetBinError(ptBin);
    }
    TGraphErrors *gExtra = new TGraphErrors(pp3.size(),&(pp3.front()),&(scale.front()),&(pp3Err.front()),&(scaleErr.front()));
    gExtra->SetMarkerStyle(20);

    double max = *(std::max_element(gExtra->GetY(),gExtra->GetY()+gExtra->GetN()));
    //TH1 *hExtraFrame = util::HistOps::createTH1D("hExtraFrame"+util::toTString(ptBin),500,0.,1.5*pp3Limits.back(),"p_{||} Threshold","","Scaling Factor");
    TH1 *hExtraFrame = util::HistOps::createTH1D("hExtraFrame"+util::toTString(ptBin),500,0.,1.5*pp3Limits.back(),"N #sigma for Tail Definition","","Scaling Factor");
    hExtraFrame->SetNdivisions(505);
    hExtraFrame->GetYaxis()->SetRangeUser(0.,2.5*max);

    TCanvas* cExtra = new TCanvas("cExtra"+util::toTString(ptBin),"Extrapolation "+util::toTString(ptBin),500,500);
    cExtra->cd();
    hExtraFrame->Draw();

    TH1* hMeanFactor = static_cast<TH1D*>(hExtraFrame->Clone("hMeanFactor"+util::toTString(ptBin)));
    hMeanFactor->SetLineWidth(2);
    TH1* hMeanFactorErr = static_cast<TH1D*>(hExtraFrame->Clone("hMeanFactor"+util::toTString(ptBin)));
    hMeanFactorErr->SetFillColor(1);
    hMeanFactorErr->SetFillStyle(3013);
    TF1* fitMeanFactor = new TF1("fitMeanFactor"+util::toTString(ptBin),"pol0");
    if( gExtra->Fit(fitMeanFactor,"0Q") == 0 ) {
      for(int bin = 1; bin <= hMeanFactor->GetNbinsX(); ++bin) {
	hMeanFactor->SetBinContent(bin,fitMeanFactor->GetParameter(0));
	hMeanFactorErr->SetBinContent(bin,fitMeanFactor->GetParameter(0));
	hMeanFactorErr->SetBinError(bin,fitMeanFactor->GetParError(0));
      }
    }
    hMeanFactorErr->Draw("E3same");
    hMeanFactor->Draw("HISTsame");
    gExtra->Draw("PE1same");

    double ptMin = hScales[0]->GetXaxis()->GetBinLowEdge(ptBin);
    double ptMax = hScales[0]->GetXaxis()->GetBinUpEdge(ptBin);
    TPaveText *label = util::LabelFactory::createPaveText(1);
    label->AddText(util::toTString(ptMin)+" < p^{ave}_{T} < "+util::toTString(ptMax)+" GeV");
    label->Draw("same");
    //cExtra->SaveAs("Tails_ScalingFactorsVsPP3_"+util::toTString(ptBin-1)+".eps","eps");
    cExtra->SaveAs("Tails_ScalingFactorsVsNSigma_Calo_"+util::toTString(ptBin-1)+".eps","eps");
  }
}

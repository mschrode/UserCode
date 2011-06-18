// $Id: $

// Plot combined resolution from photon+jet and dijet measurements
// for 2010 SUS and JME analyses.


#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"

#include "../util/utils.h"
#include "../util/ConfigParser.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


const int COLOR_SYST = 5;
const int COLOR_STAT = 1;
const int COLOR_PHOTONJET = 46;
const int COLOR_DIJETS = 38;


void getCombinedRatio(unsigned int etaBin, double ptMin, double ptMax, TF1* &fRatio, TGraphAsymmErrors* &gRatioStat, TGraphAsymmErrors* &gUncertTotal, TGraphAsymmErrors* &gUncertJES) {

  std::cout << "Entering getCombinedRatio(etaBin = " << etaBin << ")" << std::endl;

  // combined ratios and uncertainties for different eta bins
  std::vector<double> ratio;
  std::vector<double> stat;
  std::vector<double> uncertOtherD;
  std::vector<double> uncertOtherU;
  std::vector<double> uncertJESD;
  std::vector<double> uncertJESU;
  std::vector<double> uncertTotalD;
  std::vector<double> uncertTotalU;

//   // 2010-12-31
//   ratio.push_back(1.048);
//   ratio.push_back(1.058);
//   ratio.push_back(0.983);
//   ratio.push_back(1.110);

//   stat.push_back(0.012);
//   stat.push_back(0.021);
//   stat.push_back(0.032);
//   stat.push_back(0.043);

//   uncertOtherD.push_back(0.007);
//   uncertOtherD.push_back(0.006);
//   uncertOtherD.push_back(0.010);
//   uncertOtherD.push_back(0.016);  
			 
//   uncertOtherU.push_back(0.008);
//   uncertOtherU.push_back(0.007);
//   uncertOtherU.push_back(0.011);
//   uncertOtherU.push_back(0.018);

//   uncertJESD.push_back(0.035);
//   uncertJESD.push_back(0.040);
//   uncertJESD.push_back(0.050);
//   uncertJESD.push_back(0.088);

//   uncertJESU.push_back(0.031);
//   uncertJESU.push_back(0.033);
//   uncertJESU.push_back(0.052);
//   uncertJESU.push_back(0.064);


  // 2011-03-03
  ratio.push_back(1.06621);
  ratio.push_back(1.16675);
  ratio.push_back(1.08772);
  ratio.push_back(1.17224);
  
  stat.push_back(0.00654984);
  stat.push_back(0.0161794);
  stat.push_back(0.0252353);
  stat.push_back(0.038656);
  
  uncertOtherD.push_back(0.0291134);
  uncertOtherD.push_back(0.029887);
  uncertOtherD.push_back(0.045277);
  uncertOtherD.push_back(0.0455763);  
  
  uncertOtherU.push_back(0.0275251);
  uncertOtherU.push_back(0.024335);
  uncertOtherU.push_back(0.0490106);
  uncertOtherU.push_back(0.046487);
  
  uncertJESD.push_back(0.0210994);
  uncertJESD.push_back(0.0191421);
  uncertJESD.push_back(0.0234846);
  uncertJESD.push_back(0.13327);
  
  uncertJESU.push_back(0.0217332);
  uncertJESU.push_back(0.019671);
  uncertJESU.push_back(0.0248322);
  uncertJESU.push_back(0.112174);

  for(size_t i = 0; i < uncertJESU.size(); ++i) {
    uncertTotalD.push_back( sqrt( uncertJESD.at(i)*uncertJESD.at(i) + 
				  uncertOtherD.at(i)*uncertOtherD.at(i) ) );
    uncertTotalU.push_back( sqrt( uncertJESU.at(i)*uncertJESU.at(i) + 
			          uncertOtherU.at(i)*uncertOtherU.at(i) ) );
  }

  
  // create ratio and graphs as bands for syst uncertainties
  fRatio = new TF1("fRatio_Eta"+util::toTString(etaBin),"pol0",ptMin,ptMax);
  fRatio->SetParameter(0,ratio.at(etaBin));
  fRatio->SetParError(0,stat.at(etaBin));
  fRatio->SetLineWidth(1);

  double x = 0.5*(ptMin+ptMax);
  double xe = 0.5*(ptMax-ptMin);
  double y = ratio.at(etaBin);
  double yed = uncertJESD.at(etaBin);
  double yeu = uncertJESU.at(etaBin);
  gUncertJES = new TGraphAsymmErrors(1,&x,&y,&xe,&xe,&yed,&yeu);
  gUncertJES->SetFillStyle(1001);
  gUncertJES->SetFillColor(5);
  gUncertJES->SetLineColor(gUncertJES->GetFillColor());

  yed = uncertTotalD.at(etaBin);
  yeu = uncertTotalU.at(etaBin);
  gUncertTotal = new TGraphAsymmErrors(1,&x,&y,&xe,&xe,&yed,&yeu);
  gUncertTotal->SetFillStyle(1001);
  gUncertTotal->SetFillColor(COLOR_SYST);
  gUncertTotal->SetLineColor(gUncertTotal->GetFillColor());

  yed = stat.at(etaBin);
  yeu = stat.at(etaBin);
  gRatioStat = new TGraphAsymmErrors(1,&x,&y,&xe,&xe,&yed,&yeu);
  gRatioStat->SetFillStyle(3013);
  gRatioStat->SetFillColor(COLOR_STAT);
  gRatioStat->SetLineColor(gRatioStat->GetFillColor());


  std::cout << "Leaving getCombinedRatio(etaBin = " << etaBin << ")" << std::endl;
}



TGraphAsymmErrors* getPoints(const TString &fileName, const TString &mode, double &etaMin, double &etaMax) {
  std::cout << "Entering getPoints()" << std::endl;

  // Read points and errors from file
  util::ConfigParser parser(fileName.Data());

  // Eta bins
  etaMin = parser.readDouble("EtaMin",":");
  etaMax = parser.readDouble("EtaMax",":");

  // Points
  TString prefix;
  if( mode == "Combined" ) prefix = "";
  else if( mode == "PhotonJet" ) prefix = "PhotonJet";
  else if( mode == "Dijet" ) prefix = "Dijet";
  else {
    std::cerr << "ERROR in getPoints(): Unknown mode '" << mode << "'" << std::endl;
    exit(1);
  }

  std::vector<double> pt = parser.readDoubleVec((mode+"MeanPt").Data(),":");
  std::vector<double> ptErr = parser.readDoubleVec((mode+"MeanPtError").Data(),":");
  std::vector<double> ratio = parser.readDoubleVec((mode+"Ratio").Data(),":");
  std::vector<double> ratioErr = parser.readDoubleVec((mode+"RatioError").Data(),":");
    
  // Check for sanity
  assert( pt.size() == ptErr.size() );
  assert( pt.size() == ratio.size() );
  assert( pt.size() == ratioErr.size() );
  if( pt.size() == 0 ) {
    std::cerr << "WARNING: no field with name '" << (mode+"MeanPt") << "' in file '" << fileName << "'" << std::endl;
  }

  // Create TGraph of nominal values and statistical errors
  TGraphAsymmErrors* gStat = new TGraphAsymmErrors(pt.size(),&(pt.front()),&(ratio.front()),
						   &(ptErr.front()),&(ptErr.front()),
						   &(ratioErr.front()),&(ratioErr.front()));
  gStat->SetMarkerStyle(20);

  std::cout << "Leaving getPoints()" << std::endl;

  return gStat;
}



void plotCombinedResolution(const TString fileNameBase = "results/GaussCoreCombined_2011-03-03") {

  bool showTitle = false;
  TString title = "";
  if( showTitle ) {
    util::StyleSettings::paper();
    gStyle->SetTitleAlign(13);
    gStyle->SetTitleFontSize(0.1);
    gStyle->SetTitleX(0.7);
    gStyle->SetTitleH(0.038);
    title = "CMS preliminary";
  } else {
    util::StyleSettings::paperNoTitle();
  }

  std::vector<double> etaBins;
  etaBins.push_back(0.);
  etaBins.push_back(1.1);
  etaBins.push_back(1.7);
  etaBins.push_back(2.3);
  etaBins.push_back(5.0);

  std::vector<double> etaMean;
  std::vector<double> etaErr;
  for(size_t i = 0; i < etaBins.size()-1; ++i) {
    etaMean.push_back(0.5*(etaBins.at(i)+etaBins.at(1+i)));
    etaErr.push_back(0.5*(etaBins.at(1+i)-etaBins.at(i)));
  }
  std::vector<double> ratioVsEta;
  std::vector<double> systVsEtaU;
  std::vector<double> systVsEtaD;

  TH1* hRatioVsEta = new TH1D("hRatioVsEta",title+";|#eta|;#sigma(Data) / #sigma(MC)",etaBins.size()-1,&(etaBins.front()));
  hRatioVsEta->GetYaxis()->SetRangeUser(0.71,1.89);
  hRatioVsEta->SetMarkerStyle(20);

  for(unsigned int etaBin = 0; etaBin < 4; ++etaBin) {
    std::cout << "Creating plots for etaBin " << etaBin << std::endl;

    TString fileName = fileNameBase+"_Eta"+util::toTString(etaBin)+".txt";

    std::cout << "  Reading points" << std::endl;
    double etaMin = 10.;
    double etaMax = 10.;
    std::cout << "    Reading photon jets" << std::endl;
    TGraphAsymmErrors* gPhotonJet = getPoints(fileName,"PhotonJet",etaMin,etaMax);
    gPhotonJet->SetMarkerStyle(22);
    gPhotonJet->SetMarkerColor(COLOR_PHOTONJET);
    gPhotonJet->SetLineColor(gPhotonJet->GetMarkerColor());
    std::cout << "    Reading photon jets done" << std::endl;
    std::cout << "    Reading dijets" << std::endl;
    TGraphAsymmErrors* gDijet = getPoints(fileName,"Dijet",etaMin,etaMax);
    gDijet->SetMarkerColor(COLOR_DIJETS);
    gDijet->SetLineColor(gDijet->GetMarkerColor());
    gDijet->SetMarkerStyle(20);
    std::cout << "    Reading dijets done" << std::endl;
    std::cout << ">> NPhotonJet = " << gPhotonJet->GetN() << std::endl;
    std::cout << ">> NDijet = " << gDijet->GetN() << std::endl;

    double ptMin = 0.5*gPhotonJet->GetX()[0];
    double ptMax = 2.*gDijet->GetX()[gDijet->GetN()-1];
    std::cout << "  Points read" << std::endl;
    
    std::cout << "  Filling ratio vs eta" << std::endl;
    // Bands around ratio
    TGraphAsymmErrors* gSystJes = 0;
    TGraphAsymmErrors* gSystTotal = 0;
    TF1* fRatio = 0;
    TGraphAsymmErrors* gRatioStat = 0;
    getCombinedRatio(etaBin,ptMin,ptMax,fRatio,gRatioStat,gSystTotal,gSystJes);

    hRatioVsEta->SetBinContent(etaBin+1,fRatio->GetParameter(0));
    hRatioVsEta->SetBinError(etaBin+1,fRatio->GetParError(0));
    ratioVsEta.push_back(fRatio->GetParameter(0));
    systVsEtaD.push_back(gSystTotal->GetEYlow()[0]);
    systVsEtaU.push_back(gSystTotal->GetEYhigh()[0]);

    std::cout << "&&&& " << etaBin << ": " << systVsEtaD.back() << std::endl;

    std::cout << "  Done filling ratio vs eta" << std::endl;

    // Labels
    TPaveText* label = util::LabelFactory::createPaveText(3,-0.46);
    label->AddText("#sqrt{s} = 7 TeV");
    label->AddText("Anti-k_{T} (R=0.5) PF Jets");
    label->AddText(util::toTString(etaMin)+" < |#eta| < "+util::toTString(etaMax));

    TLegend* leg = util::LabelFactory::createLegendCol(5,0.52);
    leg->AddEntry(gPhotonJet,"Photon+Jet (L = 36 pb^{-1})","P");
    leg->AddEntry(gDijet,"Dijets (L = 33 pb^{-1})","P");
    leg->AddEntry(fRatio,"Combined ratio","L");
    leg->AddEntry(gRatioStat,"Stat. uncertainty","F");
    leg->AddEntry(gSystTotal,"Syst. uncertainty","F");

    // Plot
    TH1* hFrame = new TH1D("hFrame_Eta"+util::toTString(etaBin),title+";p_{T} (GeV);#sigma(Data) / #sigma(MC)",
			   1000,ptMin,ptMax);
    for(int bin = 1; bin <= hFrame->GetNbinsX(); ++bin) {
      hFrame->SetBinContent(bin,1);
    }
    hFrame->SetLineStyle(2);
    hFrame->GetXaxis()->SetMoreLogLabels();
    hFrame->GetXaxis()->SetNoExponent();
    hFrame->GetYaxis()->SetRangeUser(0.61,1.99);

    TCanvas* can = new TCanvas("Can_Eta"+util::toTString(etaBin),"Ratio Eta"+util::toTString(etaBin),500,500);
    can->cd();
    hFrame->Draw("HIST");
    gSystTotal->Draw("E2same");
    gRatioStat->Draw("E2same");
    gPhotonJet->Draw("PE1same");
    gDijet->Draw("PE1same");
    fRatio->Draw("same");
    label->Draw("same");
    leg->Draw("same");
    can->SetLogx();
    gPad->RedrawAxis();
    can->SaveAs("CombinedResolution_Eta"+util::toTString(etaBin)+".eps","eps");
  }

//   std::cout << "ratio.at(0) = " << ratio.at(0) << std::endl;
//   std::cout << "etaMean.at(0) = " << etaMean.at(0) << std::endl;
//   std::cout << "etaErr.at(0) = " << etaErr.at(0) << std::endl;

   TGraphAsymmErrors* gUncert = new TGraphAsymmErrors(ratioVsEta.size(),&(etaMean.front()),&(ratioVsEta.front()),
 						     &(etaErr.front()),&(etaErr.front()),
 						     &(systVsEtaD.front()),&(systVsEtaU.front()));
   gUncert->SetFillStyle(1001);
   gUncert->SetFillColor(COLOR_SYST);
   gUncert->SetLineColor(gUncert->GetFillColor());

   TPaveText* label = util::LabelFactory::createPaveText(6,-0.48);
   label->AddText("#sqrt{s} = 7 TeV");
   label->AddText("Anti-k_{T} (R=0.5) PF Jets");
   label->AddText("Photon+Jet (L = 36 pb^{-1})");
   label->AddText("   22 < p_{T} < 180 GeV");
   label->AddText("Dijets (L = 33 pb^{-1})");
   label->AddText("   43 < p_{T} < 1000 GeV");

   TLegend* leg = util::LabelFactory::createLegendCol(2,0.49);
   leg->AddEntry(hRatioVsEta,"Combined ratio","P");
   leg->AddEntry(gUncert,"Syst. uncertainty","F");


  TCanvas* can = new TCanvas("Can_ResVsEta","Ratio vs Eta",500,500);
  can->cd();
  hRatioVsEta->Draw("PE1");
  gUncert->Draw("E2same");
  hRatioVsEta->Draw("PE1same");
  label->Draw("same");
  leg->Draw("same");
  gPad->RedrawAxis();
  can->SaveAs("CombinedResolutionVsEta.eps","eps");
}





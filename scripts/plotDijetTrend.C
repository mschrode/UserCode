#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"

#include "globalFunctions.h"
#include "../sampleTools/BinningAdmin.h"
#include "../util/utils.h"
#include "../util/FileOps.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"





// ------------------------------------------------------------------------------
void getNTailMCTruth(const TString &fileName, const TString &outNamePrefix, double scaling, double nSigCore, double nSigTailStart, double nSigTailEnd, TH1* &hNumFit, TH1* &hNumFromExt);



// ------------------------------------------------------------------------------
void plotDijetTrend() {
  std::cout << "Setting up parameters" << std::endl;

  util::StyleSettings::paperNoTitle();
  sampleTools::BinningAdmin* binAdm = new sampleTools::BinningAdmin("input/BinningAdminTailsMC.cfg");
  


  // +++++ Input +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // Input
  const TString id = "Nominal"; 
  //const TString fileNameAsym = "results/PF/Tails_PF_nSCore20_nSTail30";
  const TString fileNameAsym = "results/PF/Tails_PF_nSCore20_nSTail30";
  const TString nameMCTruthPoins = "results/PF/Tails_PF_MCTruth.root";
  const TString jetAlgo = "PF";
  const TString ptBinVar = "p^{ave}_{T}";
  const double lumi = 32.7;

  // Extrapolation
  const double fitMin = 0.;
  const double fitMax = 0.4;

  // MC truth closure
  const bool hasMCTruthPoint = true;
  const double mcTruthScaling = 0.041;
  const double nSigCore = 2.;
  const double nSigTailStart = 3.;
  const double nSigTailEnd = 6.;


  int startBinComb = 1;
  std::vector<int> nCombBins;

  std::vector<TString> names;
//   names.push_back(fileNameAsym+"_ptSoft05_.root");
//   names.push_back(fileNameAsym+"_ptSoft07.5_.root");
//   names.push_back(fileNameAsym+"_ptSoft10_.root");
//   names.push_back(fileNameAsym+"_ptSoft12.5_.root");
//   names.push_back(fileNameAsym+"_ptSoft15_.root");
//   names.push_back(fileNameAsym+"_ptSoft17.5_.root");
//   names.push_back(fileNameAsym+"_ptSoft20_.root");
//   names.push_back(fileNameAsym+"_ptSoft22.5_.root");
//   names.push_back(fileNameAsym+"_ptSoft25_.root");
//   names.push_back(fileNameAsym+"_ptSoft27.5_.root");
//   names.push_back(fileNameAsym+"_ptSoft30_.root");

  names.push_back(fileNameAsym+"_ptSoft050_.root");
  names.push_back(fileNameAsym+"_ptSoft075_.root");
  names.push_back(fileNameAsym+"_ptSoft100_.root");
  names.push_back(fileNameAsym+"_ptSoft125_.root");
  names.push_back(fileNameAsym+"_ptSoft150_.root");
  names.push_back(fileNameAsym+"_ptSoft175_.root");
  names.push_back(fileNameAsym+"_ptSoft200_.root");
  names.push_back(fileNameAsym+"_ptSoft225_.root");
  names.push_back(fileNameAsym+"_ptSoft250_.root");
  names.push_back(fileNameAsym+"_ptSoft275_.root");
  names.push_back(fileNameAsym+"_ptSoft300_.root");


  // +++++ Parameters +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::vector<double> pt3Cuts;
  for(unsigned int i = 0; i < binAdm->nPtSoftBins(); ++i) {
    pt3Cuts.push_back(binAdm->ptSoftMax(i));
  }
  assert( pt3Cuts.size() == names.size() );

  TString outName = "Tails_"+jetAlgo+"_"+id;
  if( hasMCTruthPoint ) outName += "_Closure";

  TFile outFile(outName+".root","RECREATE");



  // Loop over eta bins
  unsigned int nEtaBins = binAdm->nEtaBins();
  if( hasMCTruthPoint ) nEtaBins = 1;
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    

    // +++++ Style +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    TString outNamePrefix = outName+"_Eta"+util::toTString(etaBin)+"_";

    TString jetLabel = "Anti-k_{T} (d=0.5) ";
    if( jetAlgo == "Calo" ) jetLabel += "Calo Jets";
    else if( jetAlgo == "PF" ) jetLabel += "PF Jets";
    TString etaBinLabel =  util::toTString(binAdm->etaMin(etaBin))+" < |#eta| < "+util::toTString(binAdm->etaMax(etaBin));

    TPaveText* labelMC = util::LabelFactory::createPaveText(3);
    labelMC->AddText("CMS Simulation,  #sqrt{s} = 7 TeV");
    labelMC->AddText(jetLabel+", "+etaBinLabel);

    TPaveText* labelData = util::LabelFactory::createPaveText(3);
    labelData->AddText("#sqrt{s} = 7 TeV, L = "+util::toTString(lumi)+" pb^{-1}");
    labelData->AddText(jetLabel+", "+etaBinLabel);


    // +++++ Number of tail events +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // Combine these bins
    nCombBins.clear();
    if( etaBin == 0 ) {
      nCombBins.push_back(5);
      nCombBins.push_back(6);
      nCombBins.push_back(5);
    }

    // Get number of events from asymmetry from file
    std::vector<TH1*> hNumMC;
    std::vector<TH1*> hNumData;
    for(unsigned int i = 0; i < names.size(); ++i) {
      hNumData.push_back(util::FileOps::readTH1(names[i],"hNTailData_Eta"+util::toTString(etaBin),"hNTailIntData_Eta"+util::toTString(etaBin)+"_File"+util::toTString(i)));
      hNumMC.push_back(util::FileOps::readTH1(names[i],"hNTailMC_Eta"+util::toTString(etaBin),"hNTailIntMC_Eta"+util::toTString(etaBin)+"_File"+util::toTString(i)));
    }

    // Get number of events from MC truth
    TH1* hNumMCTruth = 0;
    TH1* hNumMCTruthBias = 0;  
    if( hasMCTruthPoint ) 
      getNTailMCTruth(nameMCTruthPoins,outName,mcTruthScaling,nSigCore,nSigTailStart,nSigTailEnd,hNumMCTruth,hNumMCTruthBias);


    // Pt summary plots
    TH1* hAbsDiff = new TH1D("hAbsDiff_Eta"+util::toTString(etaBin),"",binAdm->nPtBins(etaBin),&((binAdm->ptBinEdges(etaBin)).front()));
    hAbsDiff->SetXTitle(ptBinVar+" (GeV)");
    hAbsDiff->SetYTitle("f_{Tail}(Data) - f_{Tail}(MC)");
    hAbsDiff->GetYaxis()->SetRangeUser(-0.029,0.049);
    hAbsDiff->SetMarkerStyle(20);
    TH1* hRelDiff = static_cast<TH1*>(hAbsDiff->Clone("hRelDiff_Eta"+util::toTString(etaBin)));
    hRelDiff->SetYTitle("(f_{Tail}(Data) - f_{Tail}(MC)) / f_{Tail}(MC)");
    hRelDiff->GetYaxis()->SetRangeUser(-0.69,0.99);

    TH1* hChi2Closure = new TH1D("hChi2AbsDiff","",binAdm->nPtBins(etaBin),&((binAdm->ptBinEdges(etaBin)).front()));
    hChi2Closure->SetXTitle(ptBinVar+" (GeV)");
    hChi2Closure->SetYTitle("#chi^{2} / ndof");
    hChi2Closure->SetMarkerStyle(25);
    hChi2Closure->SetMarkerColor(kBlue);
    hChi2Closure->SetLineColor(kBlue);

    TH1* hClosure = new TH1D("hClosure","",binAdm->nPtBins(etaBin),&((binAdm->ptBinEdges(etaBin)).front()));
    hClosure->SetXTitle(ptBinVar+" (GeV)");
    hClosure->SetYTitle("(f_{Tail}(Asym) - f_{Tail}(MCTruth)) / f_{Tail}(MCTruth)");
    hClosure->GetYaxis()->SetRangeUser(-0.029,0.049);
    hClosure->SetMarkerStyle(26);
    hClosure->SetMarkerColor(kRed);
    hClosure->SetLineColor(kRed);

    TH1* hMCTruth = new TH1D("hBiasedMCTruthResponse","",binAdm->nPtBins(etaBin),&((binAdm->ptBinEdges(etaBin)).front()));
    hMCTruth->SetXTitle(ptBinVar+" (GeV)");
    hMCTruth->SetYTitle("Biased MCTruth");

    TH1* hExtra = new TH1D("hExtrapolatedAsymmetry","",binAdm->nPtBins(etaBin),&((binAdm->ptBinEdges(etaBin)).front()));
    hExtra->SetXTitle(ptBinVar+" (GeV)");
    hExtra->SetYTitle("Extrapolated Asymmetry");


    // Loop over pt bins and create trend and ratio plots
    for(unsigned int ptBin = 0; ptBin < binAdm->nPtBins(etaBin); ++ptBin) {
      

      // +++++++++ Number of tail events +++++++++++++++++++++++++++++++++

      // Number of tail events from asymmetry for different pt3 cuts
      std::vector<double> pt3;
      std::vector<double> pt3E;
      std::vector<double> nMC;
      std::vector<double> nMCE;
      std::vector<double> nData;
      std::vector<double> nDataE;
      for(unsigned int i = 0; i < binAdm->nPtSoftBins(); ++i) {
	pt3.push_back(pt3Cuts[i]);
	pt3E.push_back(0.);
	nMC.push_back(hNumMC[i]->GetBinContent(ptBin+1));
	nMCE.push_back(hNumMC[i]->GetBinError(ptBin+1));
	nData.push_back(hNumData[i]->GetBinContent(ptBin+1));
	nDataE.push_back(hNumData[i]->GetBinError(ptBin+1));
      }
   
      TGraphAsymmErrors* gMC = new TGraphAsymmErrors(pt3.size(),&(pt3.front()),&(nMC.front()),
						     &(pt3E.front()),&(pt3E.front()),
						     &(nMCE.front()),&(nMCE.front()));
      gMC->SetMarkerStyle(25);
      gMC->SetMarkerColor(kBlue);
      gMC->SetLineColor(kBlue);

      TGraphAsymmErrors* gData = new TGraphAsymmErrors(pt3.size(),&(pt3.front()),&(nData.front()),
						       &(pt3E.front()),&(pt3E.front()),
						       &(nDataE.front()),&(nDataE.front()));
      gData->SetMarkerStyle(20);


      // +++++++++ Extrapolation +++++++++++++++++++++++++++++++++

       // Extrapolation
//          TF1* fitMC = new TF1("fitMC","pol1",fitMin,fitMax);
//          gMC->Fit(fitMC,"0QR");
//          fitMC->SetLineWidth(1);
//          fitMC->SetLineColor(gMC->GetLineColor());
  
//          TF1* fitMC = new TF1("fitMC","sqrt( sq([0]) + sq([1]*x) )",fitMin,fitMax);
      TF1* fitMC = new TF1("fitMC","[0] + sq([1])*x + [2]*sq(x)",fitMin,fitMax);
      fitMC->SetParameter(0,0.005);
      fitMC->SetParameter(1,0.);
      fitMC->SetParameter(2,0.05);
      gMC->Fit(fitMC,"0QR");
      fitMC->SetLineWidth(1);
      fitMC->SetLineColor(gMC->GetLineColor());

//       TF1* fitMC = new TF1("fitMC","pol2",fitMin,fitMax);
//       gMC->Fit(fitMC,"0QR");
//       fitMC->SetLineWidth(1);
//       fitMC->SetLineColor(gMC->GetLineColor());

      hChi2Closure->SetBinContent(1+ptBin,fitMC->GetChisquare()/fitMC->GetNDF());

      TF1* fitMCDraw = static_cast<TF1*>(fitMC->Clone("fitMCDraw"));
      fitMCDraw->SetRange(0.,fitMax);
      fitMCDraw->SetLineStyle(2);



      // +++++ Number of tail events from MC truth +++++++++++++++++++++++++++++++
      TGraphAsymmErrors* gMCTruth = 0;
      TGraphAsymmErrors* gMCTruthBias = 0;
      if( hasMCTruthPoint ) {
	// MC truth
	pt3.clear();
	pt3E.clear();
	nMC.clear();
	nMCE.clear();

	pt3.push_back(0.);
	pt3E.push_back(0.);
	nMC.push_back(hNumMCTruth->GetBinContent(1+ptBin));
	nMCE.push_back(hNumMCTruth->GetBinError(1+ptBin));
	gMCTruth = new TGraphAsymmErrors(pt3.size(),&(pt3.front()),&(nMC.front()),
					 &(pt3E.front()),&(pt3E.front()),
					 &(nMCE.front()),&(nMCE.front()));
	gMCTruth->SetMarkerStyle(23);
	gMCTruth->SetMarkerColor(hClosure->GetLineColor());
	gMCTruth->SetLineColor(hClosure->GetLineColor());

	// Bias corrected MC truth
	gMCTruthBias = static_cast<TGraphAsymmErrors*>(gMCTruth->Clone());
	gMCTruthBias->SetPoint(0,0.,hNumMCTruthBias->GetBinContent(1+ptBin));
	gMCTruthBias->SetPointError(0,0.,0.,hNumMCTruthBias->GetBinError(1+ptBin),
				    hNumMCTruthBias->GetBinError(1+ptBin));
	gMCTruthBias->SetMarkerStyle(hClosure->GetMarkerStyle());


	// Relative difference to extrapolation
	double fExtra = fitMC->GetParameter(0);
	double fExtraErr = fitMC->GetParError(0);
	double fTrue = hNumMCTruthBias->GetBinContent(1+ptBin);
	double fTrueErr = hNumMCTruthBias->GetBinError(1+ptBin);

	hExtra->SetBinContent(1+ptBin,fExtra);
	hExtra->SetBinError(1+ptBin,fExtraErr);
	hMCTruth->SetBinContent(1+ptBin,fTrue);
	hMCTruth->SetBinError(1+ptBin,fTrueErr);

	if( fTrue > 0. ) {
	  hClosure->SetBinContent(1+ptBin,fExtra/fTrue - 1.);
	  hClosure->SetBinError(1+ptBin,sqrt(pow(fExtraErr/fTrue,2)+pow(fExtra*fTrueErr/fTrue/fTrue,2)));
	}
      }


      
      // +++++++++ Data / MC difference +++++++++++++++++++++++++++++++++
    
      // Absolute and relative difference
      TH1* hDiff = new TH1D("hDiff",";#delta_{i} = f^{Data}_{Asym,i} - f^{MC}_{Asym}(p_{T,3})",1000,-fitMC->Eval(fitMax),fitMC->Eval(fitMax));
      TGraphAsymmErrors* gDiff = static_cast<TGraphAsymmErrors*>(gData->Clone());
      TGraphAsymmErrors* gRelDiff = static_cast<TGraphAsymmErrors*>(gData->Clone());
      for(int i = 0; i < gData->GetN(); ++i) {
	double x = gMC->GetX()[i];
	double ydata = gData->GetY()[i];
	double yedata = gData->GetEYhigh()[i];
	double ymc = fitMC->Eval(x);
	double yemc = 0.;	// Replace by fit error!
	
	// Absolute difference (data - MC fit)
	double y = ydata - ymc;
	double ye = sqrt(pow(yedata,2.) + pow(yemc,2.));
	gDiff->SetPoint(i,x,y);
	//gDiff->SetPointError(i,0.,0.,ye,ye);
	gDiff->SetPointError(i,0.,0.,0.,0.);
	hDiff->Fill(y);
	
	// Relative difference (data - MC)/MC
	y = (ydata - ymc)/ymc;
	ye = sqrt(pow(yedata/ymc,2.) + pow(ydata*yemc/ymc/ymc,2.));
	gRelDiff->SetPoint(i,x,y);
	gRelDiff->SetPointError(i,0.,0.,ye,ye);
      }
      gDiff->SetMarkerStyle(20);
      gRelDiff->SetMarkerStyle(20);    

      
      // Average absolute difference and spread
      double delta = hDiff->GetMean();
      double deltaErr = hDiff->GetRMS();
      hAbsDiff->SetBinContent(1+ptBin,delta);
      hAbsDiff->SetBinError(1+ptBin,deltaErr);
      delete hDiff;

      TF1* fitData = static_cast<TF1*>(fitMC->Clone("fitData"));
      fitData->SetParameter(0,fitData->GetParameter(0)+delta);      
      fitData->SetLineStyle(2);


      


      // +++++++++ Plots +++++++++++++++++++++++++++++++++

      TString ptBinLabel = util::toTString(binAdm->ptMin(etaBin,ptBin))+" < "+ptBinVar+" < "+util::toTString(binAdm->ptMax(etaBin,ptBin))+" GeV";

      TPaveText* labelMCBin = util::LabelFactory::createPaveText(3);
      labelMCBin->AddText("CMS Simulation,  #sqrt{s} = 7 TeV");
      labelMCBin->AddText(jetLabel);
      labelMCBin->AddText(etaBinLabel+", "+ptBinLabel);
    
      TPaveText* labelDataBin = util::LabelFactory::createPaveText(3);
      labelDataBin->AddText("#sqrt{s} = 7 TeV, L = "+util::toTString(lumi)+" pb^{-1}");
      labelDataBin->AddText(jetLabel);
      labelDataBin->AddText(etaBinLabel+",  "+ptBinLabel);

      TLegend* legMC = 0;
      if( hasMCTruthPoint ) {
	legMC = util::LabelFactory::createLegendColWithOffset(4,-0.8,3);
	legMC->AddEntry(gMCTruth,"MC Truth Response","P");
	legMC->AddEntry(gMCTruthBias,"MC Truth Response (Biased)","P");
	legMC->AddEntry(gMC,"MC Asymmetry","P");
	legMC->AddEntry(fitMC,"Fit","L");
      } else {
	legMC = util::LabelFactory::createLegendColWithOffset(3,-0.65,3);
	legMC->AddEntry(gMC,"MC Asymmetry","P");
	legMC->AddEntry(fitMC,"Fit","L");
      }

      TLegend* legData = 0;
      if( hasMCTruthPoint ) {
	legData = util::LabelFactory::createLegendColWithOffset(4,-0.8,3);
	legData->AddEntry(gMCTruthBias,"MC Truth Response (Biased)","P");
	legData->AddEntry(gMC,"MC Asymmetry","P");
	legData->AddEntry(fitMC,"Fit","L");
	legData->AddEntry(gData,"Data Asymmetry","P");
      } else {
	legData = util::LabelFactory::createLegendColWithOffset(3,-0.6,3);
	legData->AddEntry(gMC,"MC Asymmetry","P");
	legData->AddEntry(fitMC,"Fit","L");
	legData->AddEntry(gData,"Data Asymmetry","P");
      }

      TLegend* legData2 = 0;
      if( hasMCTruthPoint ) {
	legData2 = util::LabelFactory::createLegendColWithOffset(5,-0.8,3);
	legData2->AddEntry(gMCTruthBias,"MC Truth Response (Biased)","P");
	legData2->AddEntry(gMC,"MC Asymmetry","P");
	legData2->AddEntry(fitMC,"Fit to MC Asymmetry","L");
	legData2->AddEntry(gData,"Data Asymmetry","P");
	legData2->AddEntry(fitData,"Shifted Fit","L");
      } else {
	legData2 = util::LabelFactory::createLegendColWithOffset(4,-0.6,3);
	legData2->AddEntry(gMC,"MC Asymmetry","P");
	legData2->AddEntry(fitMC,"Fit to MC Asymmetry","L");
	legData2->AddEntry(gData,"Data Asymmetry","P");
	legData2->AddEntry(fitData,"Shifted Fit","L");
      }

      // Extrapolation MC Closure
      TCanvas* can = new TCanvas("can","Number of events",500,500);
      can->cd();
      TH1* hFrame = new TH1D("hFrame",";p_{T,3} Threshold;f_{Asym}",1000,0.,0.44);
      hFrame->GetYaxis()->SetRangeUser(0.,2.3*gMC->GetY()[gMC->GetN()-1]);
      hFrame->Draw();
      gMC->Draw("PE1same");
      fitMC->Draw("same");
      if( hasMCTruthPoint ) {
	gMCTruth->Draw("PE1same");
	gMCTruthBias->Draw("PE1same");
      }
      labelMCBin->Draw("same");
      legMC->Draw("same");
      can->SaveAs(outNamePrefix+"PtBin"+util::toTString(ptBin)+"_ExtrapolationMCClosure.eps","eps");

      // Extrapolation MC Closure + Data
      can->cd();
      hFrame->Draw();
      gMC->Draw("PE1same");
      fitMC->Draw("same");
      if( hasMCTruthPoint ) gMCTruthBias->Draw("PE1same");
      gData->Draw("PE1same");
      labelDataBin->Draw("same");
      legData->Draw("same");
      can->SaveAs(outNamePrefix+"PtBin"+util::toTString(ptBin)+"_Extrapolation.eps","eps");

      // Extrapolation MC Closure + Data + Shifted Extrapolation
      can->cd();
      hFrame->Draw();
      gMC->Draw("PE1same");
      fitMC->Draw("same");
      if( hasMCTruthPoint ) gMCTruthBias->Draw("PE1same");
      fitData->Draw("same");
      gData->Draw("PE1same");
      labelDataBin->Draw("same");
      legData2->Draw("same");
      can->SaveAs(outNamePrefix+"PtBin"+util::toTString(ptBin)+"_Extrapolation2.eps","eps");


      // Absolute difference (Delta)
      TH1* hDelta = new TH1D("hDelta",";p_{T,3} Threshold;#delta_{i} = f^{Data}_{Asym,i} - f^{MC}_{Asym}(p_{T,3})",1000,0.,hFrame->GetXaxis()->GetBinUpEdge(hFrame->GetNbinsX()));
      hDelta->SetLineColor(kBlack);
      TH1* hDeltaErr = static_cast<TH1D*>(hDelta->Clone("hDeltaErr"));
      for(int i = 1; i < hDelta->GetNbinsX(); ++i) {
	hDelta->SetBinContent(i,delta);
	hDeltaErr->SetBinContent(i,delta);
	hDeltaErr->SetBinError(i,deltaErr);
      }
      hDeltaErr->SetFillColor(40);
      hDeltaErr->SetLineColor(hDeltaErr->GetFillColor());
      hDeltaErr->SetFillStyle(1001);
      hDeltaErr->GetYaxis()->SetRangeUser(-0.023,0.033);
     
      TLegend* legDelta = util::LabelFactory::createLegendColWithOffset(2,-0.6,3);
      legDelta->AddEntry(hDelta,"#Delta = <#delta_{i}>","L");
      legDelta->AddEntry(hDeltaErr,"#sigma(#delta_{i})","F");

      can->cd();
      hDeltaErr->Draw("E3");
      hDelta->Draw("same");
      gDiff->Draw("Psame");
      labelDataBin->Draw("same");
      legDelta->Draw("same");
      can->SaveAs(outNamePrefix+"PtBin"+util::toTString(ptBin)+"_Delta.eps","eps");


      delete hDelta;
      delete hDeltaErr;
      delete hFrame;
      delete legMC;
      delete legData;
      delete legData2;
      delete can;
    } // end of loop over pt bins


  

    // +++++ Summary plots for all pt bins ++++++++++++++++++++++++++++++++++++++++++++

    TCanvas* can = new TCanvas("can","Summary",500,500);

    TH1* hAbsDiffFrame = new TH1D("hAbsDiffFrame_Eta"+util::toTString(etaBin),"",1000,binAdm->ptMin(etaBin),binAdm->ptMax(etaBin));
    for(int i = 1; i < hAbsDiffFrame->GetNbinsX(); ++i) {
      hAbsDiffFrame->SetBinContent(i,0.);
    }
    hAbsDiffFrame->SetLineStyle(2);
    hAbsDiffFrame->GetXaxis()->SetMoreLogLabels();
    hAbsDiffFrame->SetXTitle(ptBinVar+" (GeV)");
    hAbsDiffFrame->SetYTitle("f_{Tail}(Data) - f_{Tail}(MC)");
    hAbsDiffFrame->GetYaxis()->SetRangeUser(-0.019,0.039);

    can->cd();
    hAbsDiffFrame->Draw();
    hAbsDiff->Draw("PE1same");
    labelData->Draw("same");
    can->SetLogx();
    can->SaveAs(outNamePrefix+"AbsDiff.eps","eps");


    TH1* hRelDiffFrame = static_cast<TH1*>(hAbsDiffFrame->Clone("hRelDiffFrame_Eta"+util::toTString(etaBin)));
    hRelDiffFrame->SetYTitle("(f_{Tail}(Data) - f_{Tail}(MC)) / f_{Tail}(MC)");
    hRelDiffFrame->GetYaxis()->SetRangeUser(-0.69,0.99);

    can->cd();
    hRelDiffFrame->Draw();
    hRelDiff->Draw("PE1same");
    labelData->Draw("same");
    can->SetLogx();
    can->SaveAs(outNamePrefix+"RelDiff.eps","eps");

    // Rebinned differences
    if( nCombBins.size() ) {
      TH1* hRelDiffComb = util::HistOps::combineBins(hRelDiff,nCombBins,startBinComb);
      TH1* hAbsDiffComb = util::HistOps::combineBins(hAbsDiff,nCombBins,startBinComb);

      can->cd();
      hAbsDiffFrame->Draw();
      hAbsDiffComb->Draw("PE1same");
      labelData->Draw("same");
      can->SetLogx();
      can->SaveAs(outNamePrefix+"AbsDiffComb.eps","eps");
    
      can->cd();
      hRelDiffFrame->Draw();
      hRelDiffComb->Draw("PE1same");
      labelData->Draw("same");
      can->SetLogx();
      can->SaveAs(outNamePrefix+"RelDiffComb.eps","eps");

      delete hAbsDiffComb;
      delete hRelDiffComb;
    }

    delete hAbsDiffFrame;
    delete hRelDiffFrame;

    if( hasMCTruthPoint )  {
    
      TH1* hClosureFrame = new TH1D("hClosureFrame","",1000,binAdm->ptMin(etaBin),binAdm->ptMax(etaBin));
      for(int i = 1; i < hClosureFrame->GetNbinsX(); ++i) {
	hClosureFrame->SetBinContent(i,0.);
      }
      hClosureFrame->SetLineStyle(2);
      hClosureFrame->GetXaxis()->SetMoreLogLabels();
      hClosureFrame->SetXTitle(ptBinVar+" (GeV)");
      hClosureFrame->SetYTitle("( f^{MC}_{Asym}(0) - f^{MCTruth}_{Resp} ) / f^{MCTruth}_{Resp}");
      hClosureFrame->GetYaxis()->SetRangeUser(-0.99,0.39);

      can->cd();
      hClosureFrame->Draw();
      hClosure->Draw("PE1same");
      labelMC->Draw("same");
      can->SetLogx();
      can->SaveAs(outNamePrefix+"Closure.eps","eps");

      can->cd();
      hChi2Closure->GetYaxis()->SetRangeUser(0.,199.);
      hChi2Closure->Draw("P");
      labelMC->Draw("same");
      can->SetLogx();
      can->SaveAs(outNamePrefix+"Chi2Closure.eps","eps");

      if( nCombBins.size() ) {
	TH1* hClosureComb = util::HistOps::combineBins(hClosure,nCombBins,startBinComb);

	can->cd();
	hClosureFrame->Draw();
	hClosureComb->Draw("PE1same");
	labelMC->Draw("same");
	can->SetLogx();
	can->SaveAs(outNamePrefix+"ClosureComb.eps","eps");

	delete hClosureComb;
      }

      delete hClosureFrame;

    }

    
    // +++++ Write output to ROOT-file for exchange ++++++++++++++++++++++++++++++++++
    
    outFile.WriteTObject(hAbsDiff);
    outFile.WriteTObject(hExtra);
    outFile.WriteTObject(hMCTruth);
    
    
    delete hExtra;
    delete hMCTruth;
    delete hAbsDiff;
    delete hRelDiff;
    delete hClosure;
    delete hChi2Closure;
    delete can;
  } // End of loop over eta bins

  outFile.Close();
}  



void getNTailMCTruth(const TString &fileName, const TString &outNamePrefix, double scaling, double nSigCore, double nSigTailStart, double nSigTailEnd, TH1* &hNumFit, TH1* &hNumFromExt) {
  std::cout << "Getting tails of MC truth response" << std::endl;
  
  std::cout << "  Getting TH2 from file" << std::endl;
  TH2* hSymRespVsPtGen = util::FileOps::readTH2(fileName,"MeanResp_hSymRespVsPtGen","hSymRespVsPtGen");

  // Fill response distributions per pt bin
  std::cout << "  Filling response distributions per ptGen bin" << std::endl;
  util::HistVec hSymResp;
  std::vector<double> entries;
  util::HistOps::fillSlices(hSymRespVsPtGen,hSymResp,"hSymResp",entries);
  for(util::HistItConst hIt = hSymResp.begin(); hIt != hSymResp.end(); ++hIt) {
    (*hIt)->Scale(1./(*hIt)->Integral("width"));
    util::HistOps::setAxisTitles(*hIt,"Symmetrized response","","jets",true);
    (*hIt)->SetMarkerStyle(20);

//     TF1* fit = new TF1("fit","gaus",0.72,1.28);
//     (*hIt)->Fit(fit,"I0R");
//     std::cout << "++++++ " << (*hIt)->Integral("width") << " : " << fit->Integral(0.,2.) << std::endl;
  }

  // Getting tail events
  std::cout << "  Getting relative number of tail events" << std::endl;
  hNumFit = new TH1D("hNumMCTruth","",hSymRespVsPtGen->GetNbinsX(),hSymRespVsPtGen->GetXaxis()->GetXbins()->GetArray());
  util::HistOps::setAxisTitles(hNumFit,"p^{gen}_{T}","GeV","Number of tail events");

   // BIAS PF
   std::vector<double> bias;
   bias.push_back(0.0196111);
   bias.push_back(0.0110882);
   bias.push_back(0.0144756);
   bias.push_back(0.0129489);
   bias.push_back(0.0103704);
   bias.push_back(0.00881784);
   bias.push_back(0.00409732);
   bias.push_back(0.00316566);
   bias.push_back(0.00581672);
   bias.push_back(0.00219176);
   bias.push_back(0.00337096);
   bias.push_back(0.00254131);
   bias.push_back(0.000874077);
   bias.push_back(0.000304794);
   bias.push_back(0.00227166);
   bias.push_back(0.0025912);


//   // BIAS Calo
//   std::vector<double> bias;

//   bias.push_back(0.0353731);
//   bias.push_back(0.0202273);
//   bias.push_back(0.0133078);
//   bias.push_back(0.0122776);
//   bias.push_back(0.00362893);
//   bias.push_back(0.00371523);
//   bias.push_back(0.00207232);
//   bias.push_back(-0.00273322);
//   bias.push_back(-0.000663845);
//   bias.push_back(-0.00156201);
//   bias.push_back(-0.00359953);
//   bias.push_back(-0.00153167);
//   bias.push_back(-0.0010165);
//   bias.push_back(-0.00387256);
//   bias.push_back(-0.00149045);
//   bias.push_back(-0.000380134);

  hNumFromExt = new TH1D("hNumMCTruthBias","",hSymRespVsPtGen->GetNbinsX(),hSymRespVsPtGen->GetXaxis()->GetXbins()->GetArray());
  util::HistOps::setAxisTitles(hNumFromExt,"p^{gen}_{T}","GeV","Number of tail events");

  // Loop over pt bins
  for(unsigned int i = 0; i < hSymResp.size(); ++i) {
//     std::cout << "\n\nIntegral(distribution) = " << hSymResp[i]->Integral("width") << std::endl;

    // Smear symmetrized MC truth response
    double width = 0.;
    double widthErr = 1000.;
    util::HistOps::fitCoreWidth(hSymResp.at(i),nSigCore,width,widthErr);
    TH1* hMC = 0;
    func::smearHistogram(hSymResp.at(i),hMC,hSymResp.at(i)->GetEntries(),width,scaling);

    TH1* hTail = 0;
    TH1* hTailClean = 0;
    TF1* gauss = 0;
    double nTailEvts = func::getTail(hMC,nSigCore,nSigTailStart,nSigTailEnd,hTail,hTailClean,gauss);
    hNumFit->SetBinContent(1+i,nTailEvts);
    hNumFit->SetBinError(1+i,sqrt(nTailEvts/entries.at(i)));
//     std::cout << "NumTail = " << nTailEvts << std::endl;

    TH1* hTailBias = 0;
    TH1* hTailCleanBias = 0;
    TF1* gaussBias = 0;
    TF1* gaussBiasCorr = static_cast<TF1*>(gauss->Clone("gaussBiasCorr"));
    gaussBiasCorr->SetLineStyle(2);

    double gaussInt = gauss->Integral(0.,2.);
    double sigmaBias = gauss->GetParameter(2)+bias.at(i);
//     std::cout << "Integral " << gaussInt << "  " << gaussBiasCorr->Integral(0.,2.) << std::endl;
    std::cout << "Changing resolution from " << gaussBiasCorr->GetParameter(2) << " to " << std::flush;
    gaussBiasCorr->SetParameter(0,gaussInt/sqrt(2.*M_PI)/sigmaBias);
    gaussBiasCorr->SetParameter(2,sigmaBias);
    std::cout << gaussBiasCorr->GetParameter(2) << std::endl;
//     std::cout << "Integral " << gaussInt << "  " << gaussBiasCorr->Integral(0.,2.) << std::endl;
    nTailEvts = func::getTailFromGauss(hMC,gaussBiasCorr,nSigTailStart,nSigTailEnd,nSigCore,hTailBias,hTailCleanBias,gaussBias);
    hNumFromExt->SetBinContent(1+i,nTailEvts);
    hNumFromExt->SetBinError(1+i,sqrt(nTailEvts/entries.at(i)));
//     std::cout << "NumTail = " << nTailEvts << "\n\n" << std::endl;

    TCanvas* can = new TCanvas("can","Symmetrized MC truth response "+util::toTString(i),500,500);
    can->cd();
    hMC->SetTitle(util::toTString(hSymRespVsPtGen->GetXaxis()->GetBinLowEdge(1+i))+" < p^{gen}_{T} < "+util::toTString(hSymRespVsPtGen->GetXaxis()->GetBinUpEdge(1+i)));
    hMC->Draw("PE1");
    gaussBiasCorr->Draw("same");
    gauss->Draw("same");
    can->SaveAs(outNamePrefix+"_MCTruthLinear_PtBin"+util::toTString(i)+".eps","eps");

    hMC->Draw("PE1");
    gaussBiasCorr->Draw("same");
    gauss->Draw("same");
    can->SetLogy(1);
    can->SaveAs(outNamePrefix+"_MCTruthLog_PtBin"+util::toTString(i)+".eps","eps");

    hTailClean->Draw("PE1");
    can->SetLogy(0);
    can->SaveAs(outNamePrefix+"_MCTruthTail_PtBin"+util::toTString(i)+".eps","eps");
    
    delete hTail;
    delete hTailClean;
    delete gauss;
    delete hTailBias;
    delete hTailCleanBias;
    delete gaussBias;
    delete gaussBiasCorr;
    delete can;
  }

  hNumFit->SetMarkerStyle(26);
  hNumFromExt->SetMarkerStyle(21);
  TCanvas* can = new TCanvas("can","Number of MC truth tail events",500,500);
  can->cd();
  hNumFit->DrawClone("PE1");
  hNumFromExt->DrawClone("PE1same");
  can->SetLogx();
  can->SaveAs(outNamePrefix+"_MCTruthNumEvts.eps","eps");
    
  delete can;
}



//  Combine variations and plot final result
// ----------------------------------------------------------------
void makeFinalPlots(const TString &fileNominal) {
  

  TH1* hDelta = util::FileOps::readTH1(fileNominal,"hAbsDiff_Eta0","hDeltaNominal");
  TH1* hExtra = util::FileOps::readTH1(fileNominal,"hExtrapolatedAsymmetry","hExtrapolatedAsymmetryNominal");
  TH1* hBiasedTruth = util::FileOps::readTH1(fileNominal,"hBiasedMCTruthResponse","hBiasedMCTruthResponseNominal");

  // Nominal scaling factors
  TH1* hScale = static_cast<TH1*>(hDelta->Clone("hScale"));
  hScale->Reset();
  
  // Variations
  TH1* hScaleSpreadDelta = static_cast<TH1*>(hScale->Clone("hScaleSpreadDelta"));
  TH1* hScaleNonClosure = static_cast<TH1*>(hScale->Clone("hScaleNonClosure"));

  for(int ptBin = 1; ptBin <= hScale->GetNbinsX(); ++ptBin) {

    double fExtra = hExtra->GetBinContent(ptBin);
    double fTruth = hBiasedTruth->GetBinContent(ptBin);
    double delta = hDelta->GetBinContent(ptBin);
    double deltaErr = hDelta->GetBinError(ptBin);
    double nonClosure = std::abs(fExtra-fTruth);

    // Nominal scale
    double scale = (fExtra + delta)/fExtra;
    hScale->SetBinContent(ptBin,scale);
    hScaleSpreadDelta->SetBinContent(ptBin,scale);
    hScaleNonClosure->SetBinContent(ptBin,scale);
       
    // Variations
    double scaleVarDelta = std::abs((fExtra + delta + deltaErr)/fExtra);
    hScaleSpreadDelta->SetBinError(ptBin,std::abs(scale-0.5*(scaleVarDelta)));

    double scaleVarClosure = std::abs((fExtra + nonClosure + delta)/fExtra);
    hScaleNonClosure->SetBinError(ptBin,std::abs(scale-0.5*(scaleVarClosure)));
  }

  hScale->SetMarkerStyle(20);
  hScaleSpreadDelta->SetFillStyle(1001);
  hScaleSpreadDelta->SetFillColor(8);
  hScaleNonClosure->SetFillStyle(1001);
  hScaleNonClosure->SetFillColor(40);


  // Plot

  TH1* hFrame = util::HistOps::createRatioFrame(hScale,"Tail Scaling Factors");

  TCanvas* can = new TCanvas("can","Tail Scaling Factors",500,500);
  can->cd();
  hFrame->Draw();
  //hScaleNonClosure->Draw("HISTEsame");
  hScaleSpreadDelta->Draw("HISTEsame");
  hScale->Draw("Psame");
  //  can->SaveAs();
  
}

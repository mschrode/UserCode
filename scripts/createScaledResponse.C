#include <cassert>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TPaveText.h"

#include "../util/utils.h"
#include "../util/FileOps.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"
#include "../sampleTools/BinningAdmin.h"

#include "globalFunctions.h"



const double nSigCore_ = 2.;
const double nSigTail_ = 3.;
const double min_ = 3E-5;
const double max_ = 3E4;

bool pf_ = true;

TH1* hScalingFactors_ = 0;
TH1* hScalingFactorsUp_ = 0;
TH1* hScalingFactorsDown_ = 0;

TF1* diffDataMC_ = 0;
TF1* diffDataMCUp_ = 0;
TF1* diffDataMCDown_ = 0;



void getSmearFactors(double pt, double fitMin, double fitMax, double &scale, double &scaleUp, double &scaleDown) {
  if( diffDataMC_ == 0 ) {
    diffDataMC_ = new TF1("diffDataMC_","[0] - exp(x/[1]) - 1.",0.,2000.);
    diffDataMCUp_ = new TF1("diffDataMCUp_","[0] - exp(x/[1]) - 1.",0.,2000.);
    diffDataMCDown_ = new TF1("diffDataMCDown_","[0] - exp(x/[1]) - 1.",0.,2000.);
    if( pf_ ) {
      diffDataMC_->SetParameter(0,1.10271);
      diffDataMC_->SetParameter(1,-49.7541);
      diffDataMCUp_->SetParameter(0,1.13059);
      diffDataMCUp_->SetParameter(1,-56.6096);
      diffDataMCDown_->SetParameter(0,1.07482);
      diffDataMCDown_->SetParameter(1,-42.8986);
    } else {
      diffDataMC_->SetParameter(0,1.06352);
      diffDataMC_->SetParameter(1,-33.4677);
      diffDataMCUp_->SetParameter(0,1.07392);
      diffDataMCUp_->SetParameter(1,-30.2579);
      diffDataMCDown_->SetParameter(0,1.05313);
      diffDataMCDown_->SetParameter(1,-36.6775);
    }
  }

  if( pt < fitMin ) {
    scale = diffDataMC_->Eval(fitMin);
    scaleUp = diffDataMCUp_->Eval(fitMin);    
    scaleDown = diffDataMCDown_->Eval(fitMin);    
  } else if( pt > fitMax ) {
    scale = diffDataMC_->Eval(fitMax);
    scaleUp = diffDataMCUp_->Eval(fitMax);    
    scaleDown = diffDataMCDown_->Eval(fitMax);    
  } else {
    scale = diffDataMC_->Eval(pt);
    scaleUp = diffDataMCUp_->Eval(pt);    
    scaleDown = diffDataMCDown_->Eval(pt);    
  }
  if( scale < 0. ) scale = 0.;
  if( scaleUp < 0. ) scaleUp = 0.;
  if( scaleDown < 0. ) scaleDown = 0.;
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


double getTailScalingFactorUp(double pt) {
  return getTailScalingFactor(hScalingFactorsUp_,pt);
}


double getTailScalingFactorDown(double pt) {
  return getTailScalingFactor(hScalingFactorsDown_,pt);
}


void createScaledResponse() {
  util::StyleSettings::presentationNoTitle();
  gErrorIgnoreLevel = 1001;        // Do not print ROOT message if eps file has been created

  // Get smearing and tail scaling factors
  std::cout << "Getting smearing and tail scaling factors\n";
  TString fileScalingFactors = "results/Tails_Calo_.root";
  TString outNamePrefix = "ScaledMCTruthResponse_";
  unsigned int etaBin = 0;

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
  util::HistVec hResMC = util::FileOps::readHistVec("results/MCFall10_TruthResolution_Calo.root","hRespMeas_");


  // Sanity checks
  unsigned int nPtBins = ptGenBinEdges.size()-1;

  assert( hResMC.size() >= nPtBins );

  std::vector<double> ptGenBinCenters(nPtBins);
  for(unsigned int i = 0; i < nPtBins; ++i) {
    ptGenBinCenters[i] = 0.5*(ptGenBinEdges[i]+ptGenBinEdges[1+i]);
  }

  pf_ = fileScalingFactors.Contains("PF");

  sampleTools::BinningAdmin adm("BinningAdmin.cfg");

  if( pf_ ) {
    outNamePrefix += "PF_";
  } else {
    outNamePrefix += "Calo_";
  }
  outNamePrefix += "Eta"+util::toTString(etaBin)+"_";



  // Smear MC truth response
  std::cout << "Smearing MC truth response\n";
  std::vector<TH1*> hResSmeared(nPtBins);
  std::vector<TH1*> hResSmearedUp(nPtBins);
  std::vector<TH1*> hResSmearedDown(nPtBins);
  for(unsigned int i = 0; i < nPtBins; ++i) {
    hResMC[i]->GetXaxis()->SetRangeUser(0.,2.);
    hResMC[i]->GetYaxis()->SetRangeUser(min_,max_);
    hResMC[i]->SetLineWidth(1);
    hResMC[i]->SetLineStyle(2);
    hResMC[i]->SetTitle("");
    util::HistOps::setAxisTitles(hResMC[i],"Response","","jets",true);
    double width = 0.;
    double widthErr = 1000.;
    if( func::fitCoreWidth(hResMC[i],nSigCore_,width,widthErr) ) {
      double ratioDataMC = 1.;
      double ratioDataMCUp = 1.;
      double ratioDataMCDown = 1.;
      getSmearFactors(ptGenBinCenters[i],adm.ptMin(etaBin),adm.ptMax(etaBin),ratioDataMC,ratioDataMCUp,ratioDataMCDown);
      std::cout << "  PtBin " << i << ": Data / MC (Core Width) = " << 1+ratioDataMC << " +"  << 1+ratioDataMCUp << " -" << 1+ratioDataMCDown << std::endl;
      func::smearHistogram(hResMC[i],hResSmeared[i],hResMC[i]->GetEntries(),width,ratioDataMC);
      func::smearHistogram(hResMC[i],hResSmearedUp[i],hResMC[i]->GetEntries(),width,ratioDataMCUp);
      func::smearHistogram(hResMC[i],hResSmearedDown[i],hResMC[i]->GetEntries(),width,ratioDataMCDown);
      hResSmeared[i]->SetLineColor(kRed);
      hResSmearedUp[i]->SetLineColor(kBlue);
      hResSmearedDown[i]->SetLineColor(kBlue);
    } else {
      std::cerr << "ERROR: Fit of Gaussian core failed." << std::endl;
      exit(1);
    }
  }
  
  
  // Dividing MC truth into core and tails
  std::cout << "Dividing MC truth into core and tails\n";
  util::HistVec hCore(nPtBins);
  util::HistVec hCoreUp(nPtBins);
  util::HistVec hCoreDown(nPtBins);
  util::HistVec hTails(nPtBins);
  util::HistVec hTailsUp(nPtBins);
  util::HistVec hTailsDown(nPtBins);
  util::HistVec hTailsClean(nPtBins);
  util::HistVec hTailsCleanUp(nPtBins);
  util::HistVec hTailsCleanDown(nPtBins);
  std::vector<TF1*> fGauss(nPtBins);
  std::vector<TF1*> fGaussUp(nPtBins);
  std::vector<TF1*> fGaussDown(nPtBins);
  for(unsigned int i = 0; i < nPtBins; ++i) {
    func::getTail(hResSmeared[i],nSigCore_,nSigTail_,hTails[i],hTailsClean[i],fGauss[i]);
    hCore[i] = static_cast<TH1D*>(hResSmeared[i]->Clone("hCore"+util::toTString(i)));
    hCore[i]->Add(hTails[i],-1.);
    hCore[i]->SetLineStyle(1);
    hCore[i]->SetLineColor(1);
    hResSmeared[i]->SetLineStyle(1);
    hResSmeared[i]->SetLineColor(kBlack);
    hTails[i]->SetLineStyle(1);
    hTails[i]->SetLineColor(kBlack);
    hTailsClean[i]->SetLineStyle(1);
    fGauss[i]->SetLineColor(kRed);

    func::getTail(hResSmearedUp[i],nSigCore_,nSigTail_,hTailsUp[i],hTailsCleanUp[i],fGaussUp[i]);
    hCoreUp[i] = static_cast<TH1D*>(hResSmearedUp[i]->Clone("hCoreUp"+util::toTString(i)));
    hCoreUp[i]->Add(hTails[i],-1.);
    hCoreUp[i]->SetLineStyle(1);
    hCoreUp[i]->SetLineColor(1);
    hResSmearedUp[i]->SetLineStyle(1);
    hTailsUp[i]->SetLineStyle(1);
    hTailsCleanUp[i]->SetLineStyle(1);
    fGaussUp[i]->SetLineColor(kBlue);
    fGaussUp[i]->SetLineStyle(1);

    func::getTail(hResSmearedDown[i],nSigCore_,nSigTail_,hTailsDown[i],hTailsCleanDown[i],fGaussDown[i]);
    hCoreDown[i] = static_cast<TH1D*>(hResSmearedDown[i]->Clone("hCoreDown"+util::toTString(i)));
    hCoreDown[i]->Add(hTails[i],-1.);
    hCoreDown[i]->SetLineStyle(1);
    hCoreDown[i]->SetLineColor(1);
    hResSmearedDown[i]->SetLineStyle(1);
    hTailsDown[i]->SetLineStyle(1);
    hTailsCleanDown[i]->SetLineStyle(1);
    fGaussDown[i]->SetLineColor(kBlue);
    fGaussDown[i]->SetLineStyle(1);
  }

  
  // Scaled response
  std::cout << "Scaleing response\n";
  util::HistVec hResScaled(nPtBins);
  util::HistVec hResScaledUp(nPtBins);
  util::HistVec hResScaledDown(nPtBins);
  util::HistVec hTailsScaled(nPtBins);
  util::HistVec hTailsScaledUp(nPtBins);
  util::HistVec hTailsScaledDown(nPtBins);
  for(unsigned int i = 0; i < nPtBins; ++i) {
    double tailScale = getTailScalingFactor(ptGenBinCenters[i]);
    double tailScaleUp = getTailScalingFactorUp(ptGenBinCenters[i]);
    double tailScaleDown = getTailScalingFactorDown(ptGenBinCenters[i]);
    std::cout << "  PtBin " << i << ": Tail Scaling Factor " << tailScale << " +" << tailScaleUp << " -" << tailScaleDown << std::endl;

    hTailsScaled[i] = static_cast<TH1D*>(hTails[i]->Clone("hTailsScaled"+util::toTString(i)));
    hTailsScaled[i]->SetLineColor(kRed);
    hTailsScaled[i]->Scale(tailScale);
    hTailsScaledUp[i] = static_cast<TH1D*>(hTailsUp[i]->Clone("hTailsScaledUp"+util::toTString(i)));
    hTailsScaledUp[i]->SetLineColor(kBlue);
    hTailsScaledUp[i]->Scale(tailScaleUp);
    hTailsScaledDown[i] = static_cast<TH1D*>(hTailsDown[i]->Clone("hTailsScaledDown"+util::toTString(i)));
    hTailsScaledDown[i]->SetLineColor(kBlue);
    hTailsScaledDown[i]->Scale(tailScaleDown);

    hResScaled[i] = static_cast<TH1D*>(hResSmeared[i]->Clone("hResScaled"+util::toTString(i)));
    hResScaled[i]->Reset();
    hResScaled[i]->SetLineColor(kRed);
    hResScaledUp[i] = static_cast<TH1D*>(hResSmearedUp[i]->Clone("hResScaledUp"+util::toTString(i)));
    hResScaledUp[i]->SetLineColor(kBlue);
    hResScaledDown[i] = static_cast<TH1D*>(hResSmearedDown[i]->Clone("hResScaledDown"+util::toTString(i)));
    hResScaledDown[i]->SetLineColor(kBlue);
    for(int bin = 1; bin <= hResScaled[i]->GetNbinsX(); ++bin) {
      hResScaled[i]->SetBinContent(bin,hCore[i]->GetBinContent(bin)+hTailsScaled[i]->GetBinContent(bin));
      hResScaledUp[i]->SetBinContent(bin,hCoreUp[i]->GetBinContent(bin)+hTailsScaledUp[i]->GetBinContent(bin));
      hResScaledDown[i]->SetBinContent(bin,hCoreDown[i]->GetBinContent(bin)+hTailsScaledDown[i]->GetBinContent(bin));
    }
  }


  // Plots
  TLegend* legSmear = util::LabelFactory::createLegendCol(3,0.6);
  legSmear->AddEntry(hResMC[0],"MC truth","L");
  legSmear->AddEntry(hResSmeared[0],"Smeared MC truth","L");
  legSmear->AddEntry(hResSmearedUp[0],"Smear uncertainty","L");

  TLegend* legParts = util::LabelFactory::createLegendCol(4,0.6);
  legParts->AddEntry(hCore[0],"Smeared core","L");
  legParts->AddEntry(hTails[0],"Smeared tails","L");
  legParts->AddEntry(hTailsScaled[0],"+ tail scaling","L");
  legParts->AddEntry(hTailsScaledUp[0],"Scaling uncertainty","L");

  TLegend* legScaled = util::LabelFactory::createLegendCol(3,0.6);
  legScaled->AddEntry(hResMC[0],"MC truth","L");
  legScaled->AddEntry(hResScaled[0],"+ tail scaling","L");
  legScaled->AddEntry(hResScaledUp[0],"Scaling uncertainty","L");

  for(unsigned int i = 0; i < nPtBins; ++i) {
    TPaveText* label = util::LabelFactory::createPaveText(3,-0.4);
    label->AddText(util::LabelFactory::labelJetAlgo(outNamePrefix));
    label->AddText(util::LabelFactory::labelEta(adm.etaMin(etaBin),adm.etaMax(etaBin)));
    label->AddText(util::toTString(ptGenBinEdges[i])+" < p^{gen}_{T} < "+util::toTString(ptGenBinEdges[1+i]));

    TCanvas* canSmear = new TCanvas("canSmear"+util::toTString(i),"Smeared "+util::toTString(i),500,500);
    canSmear->cd();
    hResSmeared[i]->SetLineColor(kRed);
    hResSmeared[i]->GetXaxis()->SetRangeUser(0.4,1.6);
    hResSmeared[i]->GetYaxis()->SetRangeUser(0.,12.);
    hResSmeared[i]->Draw("HIST");
    hResMC[i]->Draw("HISTsame");
    hResSmearedUp[i]->Draw("HISTsame");
    hResSmearedDown[i]->Draw("HISTsame");
    hResSmeared[i]->Draw("HISTsame");
    legSmear->Draw("same");
    label->Draw("same");
    canSmear->SetLogy(0);
    canSmear->SaveAs(outNamePrefix+"SmearedLinear_PtBin"+util::toTString(i)+".eps","eps");

    canSmear->cd();
    hResSmeared[i]->GetXaxis()->UnZoom();
    hResSmeared[i]->GetYaxis()->SetRangeUser(min_,max_);
    hResSmeared[i]->Draw("HIST");
    hResMC[i]->Draw("HISTsame");
    fGaussUp[i]->Draw("same");
    fGaussDown[i]->Draw("same");
    fGauss[i]->Draw("same");
    hResSmearedUp[i]->Draw("HISTsame");
    hResSmearedDown[i]->Draw("HISTsame");
    hResSmeared[i]->Draw("HISTsame");
    legSmear->Draw("same");
    label->Draw("same");
    canSmear->SetLogy();
    canSmear->SaveAs(outNamePrefix+"Smeared_PtBin"+util::toTString(i)+".eps","eps");
    
    TCanvas* canParts = new TCanvas("canParts"+util::toTString(i),"Parts "+util::toTString(i),500,500);
    canParts->cd();
    hCore[i]->GetYaxis()->SetRangeUser(min_,max_);
    hCore[i]->Draw("HIST");
    hTails[i]->SetLineColor(hCore[i]->GetLineColor());
    hTails[i]->SetLineStyle(2);
    hTails[i]->Draw("HISTsame");
    hTailsScaledUp[i]->Draw("HISTsame");
    hTailsScaledDown[i]->Draw("HISTsame");
    hTailsScaled[i]->Draw("HISTsame");
    legParts->Draw("same");
    label->Draw("same");
    canParts->SetLogy();
    canParts->SaveAs(outNamePrefix+"TailParts_PtBin"+util::toTString(i)+".eps","eps");

    canParts->cd();
    hCore[i]->GetYaxis()->SetRangeUser(min_,0.35);
    hCore[i]->Draw("HIST");
    hTails[i]->Draw("HISTsame");
    hTailsScaledUp[i]->Draw("HISTsame");
    hTailsScaledDown[i]->Draw("HISTsame");
    hTailsScaled[i]->Draw("HISTsame");
    legParts->Draw("same");
    label->Draw("same");
    canParts->SetLogy(0);
    canParts->SaveAs(outNamePrefix+"TailPartsLinear_PtBin"+util::toTString(i)+".eps","eps");
    
    TCanvas* canScaled = new TCanvas("canScaled"+util::toTString(i),"Scaled Resp "+util::toTString(i),500,500);
    canScaled->cd();
    hResMC[i]->GetYaxis()->SetRangeUser(min_,max_);
    hResMC[i]->Draw("HIST");
    hResScaledUp[i]->Draw("HISTsame");
    hResScaledDown[i]->Draw("HISTsame");
    hResScaled[i]->Draw("HISTsame");
    legScaled->Draw("same");
    label->Draw("same");
    canScaled->SetLogy();
    canScaled->SaveAs(outNamePrefix+"Scaled_PtBin"+util::toTString(i)+".eps","eps");
  }
}

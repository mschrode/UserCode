// $Id: makeTailScaleFactorPlots_SUS-12-011.C,v 1.1 2012/07/23 12:24:38 mschrode Exp $
//
// Produce plots for SUS-12-011 from .root files in same directory
// Requires tag 'SUS-12-011' of the 
//  UserCode/mschrode/sampleTools
//  UserCode/mschrode/util
// directories.


#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

#include "TArrow.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TString.h"

#include "TObject.h"
#include "TList.h"
#include "TCollection.h"

#include "../../sampleTools/BinningAdmin.h"
#define UTILS_AS_HEADER_FILE
#include "../../util/utils.h"
#include "../../util/HistOps.h"
#include "../../util/FileOps.h"
#include "../../util/LabelFactory.h"
#include "../../util/StyleSettings.h"


const bool SHOW_HEADER = false;
const TString LUMI = "4.90 fb^{-1}";//util::StyleSettings::luminosity(4600);

const bool DEBUG = false;
//const double BINLABEL_WIDTH = -0.48;
const double BINLABEL_WIDTH = -0.52;
const double LEG_WIDTH = 0.48;
const TString PT3RELVAR = "p_{T,3} / p^{ave}_{T} threshold";
const double PT3PLOTMAX = 0.23;
const TString FASYM = "f_{asym}";
const TString FASYMMC = "f^{mc}_{asym}";
const TString FASYMDATA = "f^{data}_{asym}";

const TString LUMI_LABEL = SHOW_HEADER ? "CMS preliminary, L = "+LUMI+",  #sqrt{s} = 7 TeV" : "#sqrt{s} = 7 TeV,  L = "+LUMI;
const TString HEADER = SHOW_HEADER ? LUMI_LABEL : "";

const int COLOR_GAUSS = 46;
const int COLOR_FILLED_ASYM = 38;
const int COLOR_FILLED_ASYM_SMEAR = 29;
const int COLOR_LINE_ASYM_SMEAR = 30;
const double MARKER_SIZE = SHOW_HEADER ? 1. : 1.4;
const int LINE_WIDTH = SHOW_HEADER ? 1 : 2;



// ------------------------------------------------------------------------------------
TGraphAsymmErrors* nomRatio(const TH1* hNom, const TH1* hPtMean) {

  std::vector<double> x;
  std::vector<double> xe;
  std::vector<double> y;
  std::vector<double> ye;
  for(int bin = 1; bin <= hNom->GetNbinsX(); ++bin) {
    x.push_back(hPtMean->GetBinContent(bin));
    xe.push_back(hPtMean->GetBinError(bin));
    y.push_back(hNom->GetBinContent(bin));
    ye.push_back(hNom->GetBinError(bin));
  }
  TGraphAsymmErrors* g = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),
					       &(xe.front()),&(xe.front()),&(ye.front()),&(ye.front()));
  g->SetMarkerStyle(20);
  g->SetMarkerSize(MARKER_SIZE);
  g->SetLineWidth(LINE_WIDTH);

  return g;
}


// ------------------------------------------------------------------------------------
TGraphAsymmErrors* relUncertainty(const TH1* hNom, int color, double weight, const TH1* hVarUp, const TH1* hVarDown = 0) {

  std::vector<double> x;
  std::vector<double> xe;
  std::vector<double> y;
  std::vector<double> ye;
  for(int bin = 1; bin <= hNom->GetNbinsX(); ++bin) {
    x.push_back(hNom->GetBinCenter(bin));
    xe.push_back(0.5*hNom->GetBinWidth(bin));
    y.push_back(0.);
    double nom = hNom->GetBinContent(bin);
    double var = hVarUp->GetBinContent(bin);
    double err = std::abs((var-nom)/nom);
    if( hVarDown ) {
      var = hVarDown->GetBinContent(bin);
      err += std::abs((var-nom)/nom);
      err /= 2.;
    }
    err *= weight;
    ye.push_back(err);
  }
  TGraphAsymmErrors* g = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),
					       &(xe.front()),&(xe.front()),&(ye.front()),&(ye.front()));
  g->SetFillStyle(1001);
  g->SetFillColor(color);
  g->SetLineColor(color);

  return g;
}


// ------------------------------------------------------------------------------------
TGraphAsymmErrors* totalUncertainty(const std::vector<TGraphAsymmErrors*> &uncerts, int color) {

  TGraphAsymmErrors* g = static_cast<TGraphAsymmErrors*>(uncerts.at(0)->Clone());
  for(int i = 0; i < g->GetN(); ++i) {
    double eu2 = 0.;
    double ed2 = 0.;
    for(unsigned int j = 0; j < uncerts.size(); ++j) {
      eu2 += pow(uncerts.at(j)->GetEYhigh()[i],2); 
      ed2 += pow(uncerts.at(j)->GetEYlow()[i],2);
    }
    g->SetPointError(i,g->GetEXlow()[i],g->GetEXhigh()[i],sqrt(ed2),sqrt(eu2));
  }
//   g->SetFillStyle(0);
//   g->SetFillColor(0);
//   g->SetLineColor(color);

  g->SetFillStyle(3444);
  g->SetFillColor(kBlack);
  g->SetLineColor(kBlack);


  return g;
}


// ------------------------------------------------------------------------------------
TGraphAsymmErrors* uncertaintyBand(const TH1* hNom, const TGraphAsymmErrors* gRelUncert, int color) {
  TGraphAsymmErrors* g = static_cast<TGraphAsymmErrors*>(gRelUncert->Clone());
  for(int i = 0; i < g->GetN(); ++i) {
    double y = hNom->GetBinContent(1+i);
    g->SetPoint(i,g->GetX()[i],y);
    double eu = (g->GetEYhigh()[i])*y;
    double ed = (g->GetEYlow()[i])*y;
    g->SetPointError(i,g->GetEXlow()[i],g->GetEXhigh()[i],ed,eu);
  }
  g->SetFillStyle(1001);
  g->SetFillColor(color);
  g->SetLineColor(color);

  return g;
}
  

// ------------------------------------------------------------------------------------
void plotFinalResult() {
  if( SHOW_HEADER ) {
    util::StyleSettings::setStylePAS();
  } else {
    util::StyleSettings::setStyleNoteNoTitle();
  }
  gErrorIgnoreLevel = 1001;

  sampleTools::BinningAdmin* binAdm = new sampleTools::BinningAdmin("../../resolutionFit/config/Analysis2011/Binning/BinningAdmin2011_Tails_mergedPtBins.cfg");  

  const TString fileNamePrefix = "Tail_163337-180252_Sig25-Inf_PF";
  const TString outNamePrefix = "Tail_163337-180252_SUS-12-011_Sig25-Inf_PF_ScaleFactors";

  //  std::vector<EtaPtBin*> etaPtBins;
  unsigned int nEtaBins = binAdm->nEtaBins();
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {

    // ****** Nominal scaling factors *****************************************************
    
    // Read mean pt
    TString histName = "hPtAveMeanData_Eta"+util::toTString(etaBin);
    TH1* hPtAveMeanData = util::FileOps::readTH1(fileNamePrefix+".root",histName,histName);

    // Read nominal scaling factors
    histName = "hScaleFrame_Eta"+util::toTString(etaBin);
    TH1* hScaleFrame = util::FileOps::readTH1(fileNamePrefix+".root",histName,histName);
    hScaleFrame->GetXaxis()->SetTitle("p^{ave}_{T} [GeV]");
    hScaleFrame->SetLineWidth(LINE_WIDTH);
    hScaleFrame->SetLineStyle(2);
    histName = "hScale_Eta"+util::toTString(etaBin);
    TH1* hScaleNom = util::FileOps::readTH1(fileNamePrefix+".root",histName,histName+"_Nominal");

    // Graph of nominal scaling factors with ptAveMean as x
    TGraphAsymmErrors* gNomRatio = nomRatio(hScaleNom,hPtAveMeanData);


    // ****** Variations for systematic uncertainties *************************************

    // Read histograms with variations
    TH1* hScaleVarCoreUp = util::FileOps::readTH1(fileNamePrefix+"_VarCoreUp.root",histName,histName+"_VarCoreUp");
    TH1* hScaleVarCoreDown = util::FileOps::readTH1(fileNamePrefix+"_VarCoreDown.root",histName,histName+"_VarCoreDown");
    TH1* hScaleVarClosure = util::FileOps::readTH1(fileNamePrefix+"_VarClosure.root",histName,histName+"_VarClosure");
    TH1* hScaleVarPUUp = util::FileOps::readTH1(fileNamePrefix+"_VarPUUp.root",histName,histName+"_VarPUUp");
    TH1* hScaleVarPUDown = util::FileOps::readTH1(fileNamePrefix+"_VarPUDown.root",histName,histName+"_VarPUDown");
    TH1* hScaleVarExtra = util::FileOps::readTH1(fileNamePrefix+"_VarExtrapolation.root",histName,histName+"_VarExtrapolation");

    // Get relative uncertainties
    std::vector<TGraphAsymmErrors*> uncerts;
    uncerts.push_back(relUncertainty(hScaleNom,46,1.,hScaleVarCoreUp,hScaleVarCoreDown));
    uncerts.push_back(relUncertainty(hScaleNom,11,1.,hScaleVarPUUp,hScaleVarPUDown));
    uncerts.push_back(relUncertainty(hScaleNom,38,1.,hScaleVarExtra));
    uncerts.push_back(relUncertainty(hScaleNom,8,0.5,hScaleVarClosure));

    // Define labels
    std::vector<TString> uncertLabels;
    uncertLabels.push_back("Core adaption");
    uncertLabels.push_back("Pile-up");
    uncertLabels.push_back("Extrapolation");
    uncertLabels.push_back("Non-closure");

    // Add up uncertainties
    TGraphAsymmErrors* gUncertRelTotal = totalUncertainty(uncerts,kBlack);
    TGraphAsymmErrors* gUncertAbs = uncertaintyBand(hScaleNom,gUncertRelTotal,5);
    TGraphAsymmErrors* gUncertAbsPos = static_cast<TGraphAsymmErrors*>(gUncertAbs->Clone());
    for(int i = 0; i < gUncertAbsPos->GetN(); ++i) {
      gUncertAbsPos->SetPointError(i,gUncertAbsPos->GetEXlow()[i],gUncertAbsPos->GetEXhigh()[i],
				   gUncertAbsPos->GetEYlow()[i]>gUncertAbsPos->GetY()[i] ? gUncertAbsPos->GetY()[i] : gUncertAbsPos->GetEYlow()[i],gUncertAbsPos->GetEYhigh()[i]);
    }

    // ****** Plotting ********************************************************************

    // Label
    TPaveText* label = 0;
    if( SHOW_HEADER ) {
      label = util::LabelFactory::createPaveText(3,BINLABEL_WIDTH);
    } else {
      label = util::LabelFactory::createPaveText(5,BINLABEL_WIDTH);
      label->AddText("CMS");
      label->AddText(LUMI_LABEL);
    }
    label->AddText(util::LabelFactory::labelJet("AK5PF"));
    label->AddText(util::LabelFactory::etaCut(binAdm->etaMin(etaBin),binAdm->etaMax(etaBin)));
    if( fileNamePrefix.Contains("Sig15") ) {
      label->AddText("Tail: A_{tail} = 1.5 #sigma_{A}");
    } else if( fileNamePrefix.Contains("Sig25") ) {
      label->AddText("Tail: A_{tail} = 2.5 #sigma_{A}");
    } else if( fileNamePrefix.Contains("Sig35") ) {
      label->AddText("Tail: A_{tail} = 3.5 #sigma_{A}");
    }

    TLegend* legScale = util::LabelFactory::createLegendCol(2,LEG_WIDTH);
    legScale->AddEntry(gNomRatio,"Stat. uncertainty","L");
    legScale->AddEntry(gUncertAbs,"Syst. uncertainty","F");

    TLegend* legUncert = util::LabelFactory::createLegendCol(uncerts.size()+1,LEG_WIDTH);
    for(unsigned int i = 0; i < uncerts.size(); ++i) {
      legUncert->AddEntry(uncerts.at(i),uncertLabels.at(i),"F");
    }
    legUncert->AddEntry(gUncertRelTotal,"Total","F");


    // Plot scaling factors and total uncertainty    
    hScaleFrame->SetTitle(HEADER);
    hScaleFrame->GetYaxis()->SetRangeUser(0.01,2.99);
    if( gNomRatio->GetY()[0] > 2. ) hScaleFrame->GetYaxis()->SetRangeUser(0.01,6.99);
    hScaleFrame->GetYaxis()->SetTitle("Scale factor");
    hScaleFrame->GetXaxis()->SetMoreLogLabels();
    hScaleFrame->GetXaxis()->SetNoExponent();
    TCanvas* canScale = new TCanvas("canScale"+util::toTString(etaBin),"Scale Factors Eta "+util::toTString(etaBin),500,500);
    canScale->cd();
    hScaleFrame->Draw("HIST");
    gUncertAbs->Draw("E2same");
    hScaleFrame->Draw("HISTsame");
    gNomRatio->Draw("PE1same");
    label->Draw("same");
    legScale->Draw("same");
    canScale->SetLogx();
    gPad->RedrawAxis();
    canScale->SaveAs(outNamePrefix+"_Eta"+util::toTString(etaBin)+".eps","eps");
    canScale->SaveAs(outNamePrefix+"_Eta"+util::toTString(etaBin)+".png");

    // Plot relative uncertainties
    TH1* hUncertsFrame = static_cast<TH1D*>(hScaleFrame->Clone("hUncertsFrame"+util::toTString(etaBin)));
    for(int bin = 1; bin <= hUncertsFrame->GetNbinsX(); ++bin) {
      hUncertsFrame->SetBinContent(bin,0.);
    }
    hUncertsFrame->GetYaxis()->SetRangeUser(-0.99,1.49);
    hUncertsFrame->GetXaxis()->SetMoreLogLabels();
    hUncertsFrame->GetXaxis()->SetNoExponent();
    hUncertsFrame->GetYaxis()->SetTitle("Relative uncertainty");
    TCanvas* canRelUncerts = new TCanvas("canRelUncerts"+util::toTString(etaBin),"Relative Uncertainties Eta "+util::toTString(etaBin),500,500);
    canRelUncerts->cd();
    hUncertsFrame->Draw("HIST");
    for(unsigned int i = 0; i < uncerts.size(); ++i) {
      uncerts.at(i)->Draw("E2same");
    }
    gUncertRelTotal->Draw("E2same");
    hUncertsFrame->Draw("HISTsame");
    label->Draw("same");
    legUncert->Draw("same");
    canRelUncerts->SetLogx();
    gPad->RedrawAxis();
    //    canRelUncerts->SaveAs(outNamePrefix+"_Uncertainties_Eta"+util::toTString(etaBin)+".eps","eps");

    
    // store factors and total uncertainties for evolution plots
    TH1* hOutStat = new TH1D("ScaleFactors_StatUncert_Eta"+util::toTString(etaBin),";PtAveBin;ScaleFactor",gNomRatio->GetN(),0,gNomRatio->GetN());
    TH1* hOutSyst = new TH1D("ScaleFactors_SystUncert_Eta"+util::toTString(etaBin),";PtAveBin;ScaleFactor",gNomRatio->GetN(),0,gNomRatio->GetN());
    for(int i = 0; i < gNomRatio->GetN(); ++i) {
      hOutStat->SetBinContent(1+i,gNomRatio->GetY()[i]);
      hOutStat->SetBinError(1+i,gNomRatio->GetEYhigh()[i]);
      hOutSyst->SetBinContent(1+i,gNomRatio->GetY()[i]);
      hOutSyst->SetBinError(1+i,gUncertAbs->GetEYhigh()[i]);
    }
    TFile* out = new TFile(outNamePrefix+".root","UPDATE");
    out->WriteTObject(hOutStat);
    out->WriteTObject(hOutSyst);
    out->Close();
    delete out;
    delete hOutStat;
    delete hOutSyst;
    

    // ****** Print results ***************************************************************

    // Print factors with total uncertainty (stat + syst)
    if( etaBin == 0 ) {
      std::cout << "\\begin{tabular}{ccccc}\n\\hline\n";
      std::cout << "$|\\eta|$ & $\\ptave\\,(\\gevnospace)$ & $\\mean{\\ptave}\\,(\\gevnospace)$ & Scale Factor ($\\pm\\text{stat}{}^{+\\text{syst}}_{-\\text{syst}}$) & Scale Factor ($\\pm\\text{total}$) \\\\\n";
      std::cout << "\\hline\n";
    }
    for(int bin = 1; bin <= hScaleNom->GetNbinsX(); ++bin) {
      cout.setf(ios::fixed,ios::floatfield);
      std::cout << std::setprecision(1) << "  $" << util::toTString(binAdm->etaMin(etaBin)) << " - " << util::toTString(binAdm->etaMax(etaBin)) << "$ & $";
      std::cout << std::setprecision(0) << util::toTString(binAdm->ptMin(etaBin,bin-1)) << " - " << util::toTString(binAdm->ptMax(etaBin,bin-1)) << "$ & $";
      std::cout << std::setprecision(1) << hPtAveMeanData->GetBinContent(bin) << " \\pm " << hPtAveMeanData->GetBinError(bin) << "$ & $";
      std::cout << std::setprecision(3) << hScaleNom->GetBinContent(bin);
      double estat = hScaleNom->GetBinError(bin);
      double esystd = gUncertAbs->GetEYlow()[bin-1];
      double esystu = gUncertAbs->GetEYhigh()[bin-1];
      double etotd = sqrt( estat*estat + esystd*esystd );
      double etotu = sqrt( estat*estat + esystu*esystu );
      std::cout << " \\pm " << estat << "^{+" << esystu << "}_{-" << esystd << "}$ & $";
      std::cout << hScaleNom->GetBinContent(bin) << "^{+" << etotu << "}_{-" << etotd << "} $ \\\\\n";
    }    
    std::cout << "\\hline\n";
    if( etaBin == binAdm->nEtaBins()-1 ) std::cout << "\\end{tabular}\n";
  }

  delete binAdm;
}

// ------------------------------------------------------------------------------------
void rescuePlots() {
  util::StyleSettings::setStyleNoteNoTitle();

  TFile f("Tail_163337-180252_Sig25-Inf_PF.root","READ");

  TCanvas* can1 = 0;
  f.GetObject("Tail_163337-180252_Sig25-Inf_PF_EtaBin0_PtBin3_Pt3Bin2_PtSmearAsymTail",can1);
  can1->Draw();

  TPaveText* label1 = util::LabelFactory::createPaveText(6,BINLABEL_WIDTH);
  label1->AddText("CMS");
  label1->AddText(LUMI_LABEL);
  label1->AddText(util::LabelFactory::labelJet("AK5PF"));
  label1->AddText("0.0 < |#eta| < 0.5");
  label1->AddText("312 < p^{ave}_{T} < 360 GeV");
  label1->AddText("p_{T,3} / p^{ave}_{T} < 0.075");
  label1->Draw("same");

  TPaveText* cover1 = util::LabelFactory::createPaveTextWithOffset(4,0.5,4);
  cover1->Draw("same");

  can1->SaveAs("Tail_163337-180252_SUS-12-011_Sig25-Inf_PF_EtaBin0_PtBin3_Pt3Bin2_PtSmearAsymTail.eps","eps");
  can1->SaveAs("Tail_163337-180252_SUS-12-011_Sig25-Inf_PF_EtaBin0_PtBin3_Pt3Bin2_PtSmearAsymTail.png");

  TCanvas* can2 = 0;
  f.GetObject("Tail_163337-180252_Sig25-Inf_PF_EtaBin0_PtBin3_Extrapolation",can2);
  can2->Draw();

  TPaveText* label2 = util::LabelFactory::createPaveText(6,BINLABEL_WIDTH);
  label2->AddText("CMS");
  label2->AddText(LUMI_LABEL);
  label2->AddText(util::LabelFactory::labelJet("AK5PF"));
  label2->AddText("0.0 < |#eta| < 0.5");
  label2->AddText("312 < p^{ave}_{T} < 360 GeV");
  label2->AddText("Tail: A_{tail} = 2.5 #sigma_{A}");
  label2->Draw("same");

  can2->SaveAs("Tail_163337-180252_SUS-12-011_Sig25-Inf_PF_EtaBin0_PtBin3_Extrapolation.eps","eps");
  can2->SaveAs("Tail_163337-180252_SUS-12-011_Sig25-Inf_PF_EtaBin0_PtBin3_Extrapolation.png");
}

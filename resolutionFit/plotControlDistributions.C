// $Id: plotControlDistributions.C,v 1.1 2012/02/02 13:20:42 mschrode Exp $

//!  Control and n-1 distributions for
//!  histograms in Kalibri skims from
//!  sampleTools::writeDijetSkims.C

#include <iomanip>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TLine.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TString.h"

#define UTILS_AS_HEADER_FILE
#include "../sampleTools/BinningAdmin.h"
#include "../util/utils.h"
#include "../util/FileOps.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"

const int COLOR_DATA = 1;
const int COLOR_MC = 38;
const int MARKER_DATA = 20;
const TString LABEL_DATA = "Data, L = 0.84 fb^{-1}";
const TString LABEL_MC = "CMS Simulation";


// ------------------------------------------------------------
TPaveText* createLabel(double etaMin, double etaMax, double ptAveMin, double ptAveMax, const TString &selectionLabel = "NONE") {
  TPaveText* label = 0;
  if( selectionLabel == "NONE" ) {
    label = util::LabelFactory::createPaveText(2,-0.9);
  } else {
    label = util::LabelFactory::createPaveText(3,-0.9);
  }
  label->AddText("#sqrt{s} = 7 TeV,  Anti-k_{T} (R=0.5) PF Jets");
  if( ptAveMax > ptAveMin ) {
    label->AddText(util::LabelFactory::labelEta(etaMin,etaMax)+",  "+util::toTString(ptAveMin)+" < p^{ave}_{T} < "+util::toTString(ptAveMax)+" GeV");
  } else {
    label->AddText(util::LabelFactory::labelEta(etaMin,etaMax));
  }
  if( selectionLabel != "NONE" ) label->AddText(selectionLabel);

  return label;
}


// ------------------------------------------------------------
TPaveText* createLabel(const sampleTools::BinningAdmin &binningAdmin, unsigned int etaBin, unsigned int ptBin, const TString &selectionLabel = "NONE") {
  return createLabel(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin),binningAdmin.ptMin(etaBin,ptBin),binningAdmin.ptMax(etaBin,ptBin),selectionLabel);
}


// ------------------------------------------------------------
TPaveText* createLabel(const sampleTools::BinningAdmin &binningAdmin, unsigned int etaBin, const TString &selectionLabel = "NONE") {
  return createLabel(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin),0.,-1.,selectionLabel);
}


// ------------------------------------------------------------
void plotNMinus1(const TString &fileNameData, const TString &fileNameMC, const TString &histNamePrefix, const sampleTools::BinningAdmin &binningAdmin, unsigned int etaBin, unsigned int ptBin, int rebin, double xMin, double xMax, bool logY, const TString &xLabel, const TString &xUnit, const TString &yLabel, const TString &selectionLabel, const TString &outNamePrefix) {

  // Get histograms
  TString binId = "_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin);
  TH1* hData = util::FileOps::readTH1(fileNameData,histNamePrefix+binId);
  hData->SetMarkerStyle(MARKER_DATA);
  hData->SetMarkerColor(COLOR_DATA);
  hData->SetLineColor(COLOR_DATA);
  hData->SetLineWidth(2);
  hData->SetTitle("");
  hData->Rebin(rebin);
  TH1* hMC = util::FileOps::readTH1(fileNameMC,histNamePrefix+binId);
  hMC->SetFillColor(COLOR_MC);
  hMC->SetLineColor(kBlack);
  hMC->SetLineWidth(2);
  hMC->SetTitle("");
  hMC->Rebin(rebin);
  hMC->Scale(hData->Integral()/hMC->Integral());
  util::HistOps::setAxisTitles(hMC,xLabel,xUnit,yLabel);
  if( xMax > xMin ) {
    hMC->GetXaxis()->SetRangeUser(xMin,xMax);
  }

  // Create labels
  TPaveText* label = createLabel(binningAdmin,etaBin,ptBin,selectionLabel);
  TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,-0.48,label->GetSize());
  leg->AddEntry(hData,LABEL_DATA,"P");
  leg->AddEntry(hMC,LABEL_MC,"F");

  // Plots
  util::HistOps::setYRange(hMC,label->GetSize()+leg->GetNRows(),(logY ? 3E-1 : -1.));
  TCanvas* can = new TCanvas("can","N-1 plot",500,500);
  can->cd();
  hMC->Draw("HIST");
  hData->Draw("PE1same");
  label->Draw("same");
  leg->Draw("same");
  if( logY ) can->SetLogy();
  can->SaveAs(outNamePrefix+binId+".eps","eps");

  // Clean up
  delete hData;
  delete hMC;
  delete label;
  delete leg;
  delete can;
}


// ------------------------------------------------------------
void plotNMinus1Eta(const TString &fileNameData, const TString &fileNameMC, const sampleTools::BinningAdmin &binningAdmin, unsigned int ptBin, int rebin, double xMin, double xMax, bool logY, const TString &xLabel, const TString &xUnit, const TString &yLabel, const TString &selectionLabel, const TString &outNamePrefix) {

  // Get histograms
  TH1* hEtaNMin1Data = 0;
  TH1* hEtaNMin1MC = 0;
  TH1* hEtaData = 0;
  TH1* hEtaMC = 0;
  for(unsigned int etaBin = 0; etaBin < binningAdmin.nEtaBins(); ++etaBin) {
    TString binId = "_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin);
    if( etaBin == 0 ) {
      hEtaNMin1Data = util::FileOps::readTH1(fileNameData,"hEtaNMin1"+binId);
      hEtaNMin1Data->SetMarkerStyle(21);
      hEtaNMin1Data->SetMarkerColor(COLOR_DATA);
      hEtaNMin1Data->SetLineColor(COLOR_DATA);
      hEtaNMin1Data->SetLineWidth(2);
      hEtaNMin1Data->SetTitle("");
      hEtaNMin1Data->Rebin(rebin);

      hEtaNMin1MC = util::FileOps::readTH1(fileNameMC,"hEtaNMin1"+binId);
      hEtaNMin1MC->SetFillStyle(3004);
      hEtaNMin1MC->SetFillColor(kBlack);
      hEtaNMin1MC->SetLineColor(kBlack);
      hEtaNMin1MC->SetLineWidth(2);
      hEtaNMin1MC->SetTitle("");
      hEtaNMin1MC->Rebin(rebin);

      hEtaData = util::FileOps::readTH1(fileNameData,"hEta"+binId);
      hEtaData->SetMarkerStyle(MARKER_DATA);
      hEtaData->SetMarkerColor(COLOR_DATA);
      hEtaData->SetLineColor(COLOR_DATA);
      hEtaData->SetLineWidth(2);
      hEtaData->SetTitle("");
      hEtaData->Rebin(rebin);

      hEtaMC = util::FileOps::readTH1(fileNameMC,"hEta"+binId);
      hEtaMC->SetFillColor(38);
      hEtaMC->SetLineColor(kBlack);
      hEtaMC->SetLineWidth(2);
      hEtaMC->SetTitle("");
      hEtaMC->Rebin(rebin);
    } else {
      TH1* h = util::FileOps::readTH1(fileNameData,"hEtaNMin1"+binId);
      h->Rebin(rebin);
      hEtaNMin1Data->Add(h);

      h = util::FileOps::readTH1(fileNameMC,"hEtaNMin1"+binId);
      h->Rebin(rebin);
      hEtaNMin1MC->Add(h);

      h = util::FileOps::readTH1(fileNameData,"hEta"+binId);
      h->Rebin(rebin);
      hEtaData->Add(h);

      h = util::FileOps::readTH1(fileNameMC,"hEta"+binId);
      h->Rebin(rebin);
      hEtaMC->Add(h);
    }
  }
  if( hEtaNMin1MC->Integral() )
    hEtaNMin1MC->Scale(hEtaNMin1Data->Integral()/hEtaNMin1MC->Integral());
  if( hEtaMC->Integral() )
    hEtaMC->Scale(hEtaData->Integral()/hEtaMC->Integral());
  util::HistOps::setAxisTitles(hEtaNMin1MC,xLabel,xUnit,yLabel);
  if( xMax > xMin ) {
    hEtaNMin1MC->GetXaxis()->SetRangeUser(xMin,xMax);
  }

  // Create labels
  TPaveText* label = util::LabelFactory::createPaveText(3,-0.9);
  label->AddText("#sqrt{s} = 7 TeV,  Anti-k_{T} (R=0.5) PF Jets");
  label->AddText(util::toTString(binningAdmin.ptMin(0,ptBin))+" < p^{ave}_{T} < "+util::toTString(binningAdmin.ptMax(0,ptBin))+" GeV");
  label->AddText(selectionLabel);

  TPaveText* labelData = util::LabelFactory::createPaveTextWithOffset(1,-0.48,label->GetSize());
  labelData->AddText(LABEL_DATA);
  TLegend* legData = util::LabelFactory::createLegendColWithOffset(2,-0.48,label->GetSize()+1);
  legData->AddEntry(hEtaNMin1Data,"All","P");
  legData->AddEntry(hEtaData,"|#eta| selection","P");

  TPaveText* labelMC = util::LabelFactory::createPaveTextWithOffset(1,0.48,label->GetSize());
  labelMC->AddText(LABEL_MC);
  TLegend* legMC = util::LabelFactory::createLegendColWithOffset(2,0.48,label->GetSize()+1);
  legMC->AddEntry(hEtaNMin1MC,"All","F");
  legMC->AddEntry(hEtaMC,"|#eta| selection","F");

  // Plots
  util::HistOps::setYRange(hEtaNMin1MC,label->GetSize()+3,logY?3E-1:-1.);
  TCanvas* can = new TCanvas("can","N-1 plot",500,500);
  can->cd();
  hEtaNMin1MC->Draw("HIST");
  hEtaMC->Draw("HISTsame");
  hEtaNMin1Data->Draw("PE1same");
  hEtaData->Draw("PE1same");
  label->Draw("same");
  labelData->Draw("same");
  legData->Draw("same");
  labelMC->Draw("same");
  legMC->Draw("same");
  gPad->RedrawAxis();
  if( logY ) can->SetLogy();
  can->SaveAs(outNamePrefix+"_PtBin"+util::toTString(ptBin)+".eps","eps");

  // Clean up
  delete hEtaNMin1MC;
  delete hEtaNMin1Data;
  delete hEtaMC;
  delete hEtaData;
  delete label;
  delete labelData;
  delete legData;
  delete labelMC;
  delete legMC;
  delete can;
}


// ------------------------------------------------------------
void plotPtSpectrum(const TString &fileNameData, const TString &fileNameMC, const TString &histNamePrefix, const sampleTools::BinningAdmin &binningAdmin, unsigned int etaBin, int rebin, double xMin, double xMax, bool logY, const TString &xLabel, const TString &xUnit, const TString &yLabel, const TString &selectionLabel, const TString &outNamePrefix) {

  // Get histograms
  TString binId = "_EtaBin"+util::toTString(etaBin);
  TH1* hData = util::FileOps::readTH1(fileNameData,histNamePrefix+binId);
  hData->SetMarkerStyle(MARKER_DATA);
  hData->SetMarkerColor(COLOR_DATA);
  hData->SetLineColor(COLOR_DATA);
  hData->SetLineWidth(2);
  hData->SetTitle("");
  TH1* hMC = util::FileOps::readTH1(fileNameMC,histNamePrefix+binId);
  hMC->SetFillColor(COLOR_MC);
  hMC->SetLineColor(kBlack);
  hMC->SetLineWidth(2);
  hMC->SetTitle("");

  // Scale MC to data considering effective pre-scales
  std::vector<TString> trigNames = binningAdmin.triggerNames();
  unsigned int minPtBin = 1;
  unsigned int maxPtBin = 1;
  for(std::vector<TString>::const_iterator trigName = trigNames.begin();
      trigName != trigNames.end(); ++trigName) {
    double minPt = binningAdmin.ptMin(etaBin,binningAdmin.hltMinPtBin(*trigName,etaBin));
    double maxPt = binningAdmin.ptMax(etaBin,binningAdmin.hltMaxPtBin(*trigName,etaBin));
    minPtBin = hData->FindBin(minPt);
    maxPtBin = hData->FindBin(maxPt)-1;
    if( hMC->Integral(minPtBin,maxPtBin) ) {
      double norm = hData->Integral(minPtBin,maxPtBin)/hMC->Integral(minPtBin,maxPtBin);
      for(unsigned int bin = minPtBin; bin <= maxPtBin; ++bin) {
	hMC->SetBinContent(bin,norm*hMC->GetBinContent(bin));
	hMC->SetBinError(bin,norm*hMC->GetBinError(bin));
      }
    }
  }
  hData->Rebin(rebin);
  hMC->Rebin(rebin);

  util::HistOps::setAxisTitles(hMC,xLabel,xUnit,yLabel);
  if( xMax > xMin )
    hMC->GetXaxis()->SetRangeUser(xMin,xMax);

  // Create labels
  TPaveText* label = createLabel(binningAdmin,etaBin,selectionLabel);
  TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,-0.48,label->GetSize());
  leg->AddEntry(hData,LABEL_DATA,"P");
  leg->AddEntry(hMC,LABEL_MC,"F");

  // Plots
  util::HistOps::setYRange(hMC,label->GetSize()+leg->GetNRows(),(logY ? 3E-1 : -1.));
  TCanvas* can = new TCanvas("can","control plot",500,500);
  can->cd();
  hMC->Draw("HIST");
  hData->Draw("PE1same");
  label->Draw("same");
  leg->Draw("same");
  if( logY ) can->SetLogy();
  can->SaveAs(outNamePrefix+binId+".eps","eps");

  // Clean up
  delete hData;
  delete hMC;
  delete label;
  delete leg;
  delete can;
}


// ------------------------------------------------------------
void plotPt3RelRecoVsGen(const TString &fileNameMC, const sampleTools::BinningAdmin &binningAdmin, unsigned int etaBin, unsigned int ptBin, int rebin, double xMin, double xMax, bool logZ, const TString &selectionLabel, const TString &outNamePrefix) {

  // Get histograms
  TString binId = "_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin);
  TH2* h2 = util::FileOps::readTH2(fileNameMC,"hPt3RelRecoVsGen"+binId);
  h2->SetTitle("");
  h2->GetXaxis()->SetTitle("p^{gen}_{T,3} / p^{gen,ave}_{T}");
  h2->GetYaxis()->SetTitle("p_{T,3} / p^{ave}_{T}");
  if( xMax > xMin ) {
    h2->GetXaxis()->SetRangeUser(xMin,xMax);
    h2->GetYaxis()->SetRangeUser(xMin,xMax);
  }
  util::HistOps::normHist(h2,"width");

  // Create labels
  TPaveText* label = createLabel(binningAdmin,etaBin,ptBin,selectionLabel);

  // Plots
  TCanvas* can = new TCanvas("can","control plot",550,500);
  can->cd();
  can->SetMargin(gStyle->GetPadLeftMargin(),gStyle->GetPadRightMargin()+(1.*can->GetWindowWidth()/can->GetWindowHeight())-1.,gStyle->GetPadBottomMargin(),gStyle->GetPadTopMargin());
  h2->Draw("COLZ");
  label->Draw("same");
  if( logZ ) can->SetLogz();
  can->SaveAs(outNamePrefix+binId+".eps","eps");

  // Clean up
  delete h2;
  delete label;
  delete can;
}



// ------------------------------------------------------------
void plotControlDistributions() {
  gErrorIgnoreLevel = 1001;
  util::StyleSettings::setStyleNoteNoTitle();

  const TString fileNameData = "Kalibri_ControlPlots_Run2011A-163337-167151.root";
  const TString fileNameMC = "Kalibri_ControlPlots_MCSummer11_AK5PF.root";
  const sampleTools::BinningAdmin binAdmin("config/Analysis2011/Binning/BinningAdmin2011_v2.cfg");
  const TString outNamePrefix = "DijetSel_163337-167151_";
  const TString ptRelVar = "p_{T,3} / p^{ave}_{T}";
  const TString ptRelLabel = ptRelVar+" < 0.14";
  const TString deltaPhiVar = "|#Delta#phi|";
  const TString deltaPhiLabel = deltaPhiVar+" > 2.7";

  plotNMinus1Eta(fileNameData,fileNameMC,binAdmin,3,3,-5.1,-5.1,false,"#eta","","jets",deltaPhiLabel+",  "+ptRelLabel,outNamePrefix+"NMin1Eta");

  int ptBin[2] = { 3,14 };
  for(int i = 0; i < 2; ++i) {
    // N-1 plots
    plotNMinus1(fileNameData,fileNameMC,"hDeltaPhiNMin1",binAdmin,0,ptBin[i],2,2.1,3.5,false,deltaPhiVar,"","events",ptRelLabel,outNamePrefix+"NMin1DeltaPhi");
    plotNMinus1(fileNameData,fileNameMC,"hPtRelNMin1",binAdmin,0,ptBin[i],3,0.,1.,false,ptRelVar,"","events",ptRelLabel,outNamePrefix+"NMin1PtRel");
    plotPt3RelRecoVsGen(fileNameMC,binAdmin,0,ptBin[i],1,0.,0.35,false,deltaPhiLabel,outNamePrefix+"Pt3RelRecoVsGen");
  }
  
  // Spectrum control plots
  plotPtSpectrum(fileNameData,fileNameMC,"hPtAve",binAdmin,0,20,1.,0.,true,"p^{ave}_{T}","GeV","events",deltaPhiLabel+",  "+ptRelLabel,outNamePrefix+"PtAve");
  plotPtSpectrum(fileNameData,fileNameMC,"hPtJet0",binAdmin,0,20,1.,0.,false,"p_{T,1}","GeV","events",deltaPhiLabel+",  "+ptRelLabel,outNamePrefix+"PtJet1");
  plotPtSpectrum(fileNameData,fileNameMC,"hPtJet1",binAdmin,0,20,1.,0.,false,"p_{T,2}","GeV","events",deltaPhiLabel+",  "+ptRelLabel,outNamePrefix+"PtJet2");
  plotPtSpectrum(fileNameData,fileNameMC,"hPtJet2",binAdmin,0,3,1.,250.,false,"p_{T,3}","GeV","events",deltaPhiLabel+",  "+ptRelLabel,outNamePrefix+"PtJet3");

}

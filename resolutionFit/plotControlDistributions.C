// $Id: plotControlDistributions.C,v 1.4 2012/05/31 20:18:37 mschrode Exp $

//!  Control and n-1 distributions for
//!  histograms in Kalibri skims from
//!  sampleTools::writeDijetSkims.C

#include <iomanip>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
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
const TString LUMI = util::StyleSettings::luminosity(855);
TFile* OUT_FILE= 0;


// ------------------------------------------------------------
TPaveText* createLabel(double etaMin, double etaMax, double ptAveMin, double ptAveMax, const TString &selectionLabel = "NONE", const TString &dataLabel = "NONE") {
  int nLines = 1;
  if( selectionLabel != "NONE" ) ++nLines;
  if( dataLabel == "DATA" || dataLabel == "MC" ) ++nLines;

  TPaveText* label = util::LabelFactory::createPaveText(nLines,-0.9);
  if( dataLabel == "DATA" ) {
    label->AddText(util::LabelFactory::data(LUMI));
  } else if( dataLabel == "MC" ) {
    label->AddText(util::LabelFactory::mc());
  }
  if( ptAveMax > ptAveMin ) {
    label->AddText(util::LabelFactory::etaCut(etaMin,etaMax)+",  "+util::LabelFactory::ptAveCut(ptAveMin,ptAveMax));
  } else {
    label->AddText(util::LabelFactory::etaCut(etaMin,etaMax));
  }
  if( selectionLabel != "NONE" ) label->AddText(selectionLabel);

  return label;
}


// ------------------------------------------------------------
TPaveText* createLabel(const sampleTools::BinningAdmin &binningAdmin, unsigned int etaBin, unsigned int ptBin, const TString &selectionLabel = "NONE", const TString &dataLabel = "NONE") {
  return createLabel(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin),binningAdmin.ptMin(etaBin,ptBin),binningAdmin.ptMax(etaBin,ptBin),selectionLabel,dataLabel);
}


// ------------------------------------------------------------
TPaveText* createLabel(const sampleTools::BinningAdmin &binningAdmin, unsigned int etaBin, const TString &selectionLabel = "NONE", const TString &dataLabel = "NONE") {
  return createLabel(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin),0.,-1.,selectionLabel, dataLabel);
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
  leg->AddEntry(hData,util::LabelFactory::data(LUMI),"P");
  leg->AddEntry(hMC,util::LabelFactory::mc(),"F");

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
  TPaveText* label = util::LabelFactory::createPaveText(2,-0.9);
  label->AddText(util::LabelFactory::ptAveCut(binningAdmin.ptMin(0,ptBin),binningAdmin.ptMax(0,ptBin)));
  label->AddText(selectionLabel);

  TPaveText* labelData = util::LabelFactory::createPaveTextWithOffset(1,-0.48,label->GetSize());
  labelData->AddText(util::LabelFactory::data(LUMI));
  TLegend* legData = util::LabelFactory::createLegendColWithOffset(2,-0.48,label->GetSize()+1);
  legData->AddEntry(hEtaNMin1Data,"All","P");
  legData->AddEntry(hEtaData,"|#eta| selection","P");

  TPaveText* labelMC = util::LabelFactory::createPaveTextWithOffset(1,0.48,label->GetSize());
  labelMC->AddText(util::LabelFactory::mc());
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
  leg->AddEntry(hData,util::LabelFactory::data(LUMI),"P");
  leg->AddEntry(hMC,util::LabelFactory::mc(),"F");

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
void plotCorrelations(bool data, const TString &fileName, const sampleTools::BinningAdmin &binningAdmin, const TString &histNamePrefix, unsigned int etaBin, unsigned int ptBin, int rebin, double xMin, double xMax, double yMin, double yMax, bool logZ, const TString &xLabel, const TString &yLabel, const TString &selectionLabel, const TString &outNamePrefix) {

  // Get histograms
  TString binId = "_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin);
  TH2* h2 = util::FileOps::readTH2(fileName,histNamePrefix+binId);
  h2->SetTitle("");
  h2->GetXaxis()->SetTitle(xLabel);//"p^{gen}_{T,3} / p^{gen,ave}_{T}");
  h2->GetYaxis()->SetTitle(yLabel);//"p_{T,3} / p^{ave}_{T}");
  if( xMax > xMin ) {
    h2->GetXaxis()->SetRangeUser(xMin,xMax);
  }
  if( yMax > yMin ) {
    h2->GetYaxis()->SetRangeUser(yMin,yMax);
  }
  util::HistOps::normHist(h2,"width");

  // Create labels
  TPaveText* label = createLabel(binningAdmin,etaBin,ptBin,selectionLabel,(data?"DATA":"MC"));

  // Plots
  TCanvas* can = new TCanvas(outNamePrefix+binId,outNamePrefix+binId,550,500);
  can->cd();
  can->SetMargin(gStyle->GetPadLeftMargin(),gStyle->GetPadRightMargin()+(1.*can->GetWindowWidth()/can->GetWindowHeight())-1.,gStyle->GetPadBottomMargin(),gStyle->GetPadTopMargin());
  h2->Draw("COLZ");
  label->Draw("same");
  if( logZ ) can->SetLogz();
  util::FileOps::toFiles(can,OUT_FILE);

  // Clean up
  delete h2;
  delete label;
  delete can;
}


// Plot the target pile-up distribution expected
// for data.
// ------------------------------------------------------------
void plotExpectedDataPileUpDistribution(const TString &fileName, const TString &histName, const TString &outNamePrefix) {
  TH1* h = util::FileOps::readTH1(fileName,histName);
  h->SetTitle("");
  h->Sumw2();
  util::HistOps::normHist(h,"width");
  util::HistOps::setAxisTitles(h,"N_{PU}","","Probability Density");
  util::HistOps::setYRange(h,1);
  h->GetXaxis()->SetRangeUser(0,18.5);
  h->SetLineWidth(2);

  TPaveText* label = util::LabelFactory::createPaveText(1);
  label->AddText("#LTN_{PU}#GT = "+util::toTString(h->GetMean(),4)+" #pm "+util::toTString(h->GetMeanError(),4));
  
  TCanvas* c = new TCanvas(outNamePrefix+"_ExpectedPileUpDistribution","Expected pile-up distribution",500,500);
  c->cd();
  h->Draw("HIST");
  label->Draw("same");
  util::FileOps::toFiles(c,OUT_FILE);
}



// Plot the underlying di-jet pt spectrum
// ------------------------------------------------------------
void plotUnderlyingSpectra(const TString &fileName, const sampleTools::BinningAdmin &binningAdmin, const TString &outNamePrefix) {
  unsigned int ptSoftBin = 6;
  std::vector<TH1*> hists;
  TLegend* leg = util::LabelFactory::createLegendCol(binningAdmin.nEtaBins(),0.5);
  for(unsigned int etaBin = 0; etaBin < binningAdmin.nEtaBins(); ++etaBin) {
    TH1* h = util::FileOps::readTH1(fileName,"hPtGen_Eta"+util::toTString(etaBin)+"_PtSoft"+util::toTString(ptSoftBin));
    h->SetTitle("");
    util::HistOps::normHist(h,"width");
    util::HistOps::setAxisTitles(h,"p^{gen}_{T}","GeV","Probability Density",true);
    h->SetLineWidth(2);
    h->SetLineColor(util::StyleSettings::color(etaBin));
    h->SetLineStyle(etaBin);
    leg->AddEntry(h,util::LabelFactory::etaGenCut(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin)),"L");
    hists.push_back(h);
  }

  TPaveText* label = util::LabelFactory::createPaveText(3,-0.5);
  label->AddText(util::LabelFactory::mc());
  label->AddText(util::LabelFactory::deltaPhiGenCut(2.7));
  label->AddText(util::LabelFactory::pt3RelGenCut(binningAdmin.ptSoftMax(ptSoftBin)));

  // linear x-scale
  TCanvas* c1 = new TCanvas(outNamePrefix+"_PtGenDijetSpectra","Underlying di-jet spectra",500,500);
  c1->cd();
  TH1* h = hists.front();
  h->GetXaxis()->SetRangeUser(20.,1900.);
  h->GetYaxis()->SetRangeUser(3E-12,109);
  h->GetXaxis()->SetNdivisions(505);
  h->Draw("HIST");
  for(std::vector<TH1*>::reverse_iterator rit = hists.rbegin(); rit != hists.rend(); ++rit) {
    (*rit)->Draw("HISTsame");
  }
  label->Draw("same");
  leg->Draw("same");
  c1->SetLogy();
  gPad->RedrawAxis();
  util::FileOps::toFiles(c1,OUT_FILE);

  // log x-scale
  TCanvas* c2 = new TCanvas(outNamePrefix+"_PtGenDijetSpectra_Logx","Underlying di-jet spectra",500,500);
  c2->cd();
  h->GetXaxis()->SetRangeUser(20.,1900.);
  h->GetXaxis()->SetNoExponent();
  h->GetXaxis()->SetMoreLogLabels();
  h->Draw("HIST");
  for(std::vector<TH1*>::reverse_iterator rit = hists.rbegin(); rit != hists.rend(); ++rit) {
    (*rit)->Draw("HISTsame");
  }
  label->Draw("same");
  leg->Draw("same");
  c2->SetLogy();
  c2->SetLogx();
  gPad->RedrawAxis();
  util::FileOps::toFiles(c2,OUT_FILE);
}


// ------------------------------------------------------------
void plotVariedSpectrum(const std::vector<TString> &fileNames, const std::vector<TString> &sampleNames, const sampleTools::BinningAdmin &binningAdmin, const TString &outNamePrefix) {

  unsigned int etaBin = 0;
  unsigned int ptSoftBin = 6;
  std::vector<TString> legEntries;
  for(unsigned int i = 0; i < sampleNames.size(); ++i) {
    legEntries.push_back(sampleNames.at(i));
  }

  // read histograms
  std::vector<TH1*> spectra;
  for(std::vector<TString>::const_iterator name = fileNames.begin();
      name != fileNames.end(); ++name) {
    spectra.push_back(util::FileOps::readTH1(*name,"hPtGen_Eta"+util::toTString(etaBin)+"_PtSoft"+util::toTString(ptSoftBin),"Spectrum_"+util::toTString(spectra.size())));
  }

  // vary nominal histogram
  TH1* hUp = static_cast<TH1*>(spectra.front()->Clone("Spectrum_"+util::toTString(spectra.size()+1)));
  TH1* hDn = static_cast<TH1*>(spectra.front()->Clone("Spectrum_"+util::toTString(spectra.size()+2)));
  for(int bin = 1; bin <= spectra.front()->GetNbinsX(); ++bin) {
    double x = spectra.front()->GetBinCenter(bin);
    double w = pow(x/30.,0.5);
    hUp->SetBinContent(bin,w*(spectra.front()->GetBinContent(bin)));
    hUp->SetBinError(bin,w*spectra.front()->GetBinError(bin));
    w = pow(x/30.,-0.5);
    hDn->SetBinContent(bin,w*(spectra.front()->GetBinContent(bin)));
    hDn->SetBinError(bin,w*spectra.front()->GetBinError(bin));
  }    
  spectra.push_back(hUp);
  spectra.push_back(hDn);
  legEntries.push_back(sampleNames.front()+" #left(#frac{p^{gen}_{T}}{30}#right)^{+1/2}");
  legEntries.push_back(sampleNames.front()+" #left(#frac{p^{gen}_{T}}{30}#right)^{-1/2}");

  // set style
  TLegend* leg = util::LabelFactory::createLegendCol(spectra.size()+5,0.5);
  for(unsigned int i = 0; i < spectra.size(); ++i) {
    spectra.at(i)->SetTitle("");
    util::HistOps::normHist(spectra.at(i),"width");
    util::HistOps::setAxisTitles(spectra.at(i),"p^{gen}_{T}","GeV","Probability Density",true);
    util::HistOps::setYRange(spectra.at(i),2,3E-13);
    spectra.at(i)->GetXaxis()->SetNdivisions(505);
    spectra.at(i)->SetLineWidth(2);
    spectra.at(i)->SetLineColor(util::StyleSettings::color(i));
    spectra.at(i)->SetLineStyle(i);
    leg->AddEntry(spectra.at(i),legEntries.at(i),"L");
  }
  TPaveText* label = util::LabelFactory::createPaveText(4,-0.5);
  label->AddText(util::LabelFactory::mc());
  label->AddText(util::LabelFactory::etaGenCut(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin)));
  label->AddText(util::LabelFactory::deltaPhiGenCut(2.7));
  label->AddText(util::LabelFactory::pt3RelGenCut(binningAdmin.ptSoftMax(ptSoftBin)));

  TCanvas* c = new TCanvas(outNamePrefix+"_VariedPtGenDijetSpectra","Varied spectra",500,500);
  c->cd();
  TH1* h = spectra.front();
  h->GetXaxis()->SetRangeUser(20.,1900.);
  h->GetYaxis()->SetRangeUser(3E-12,109);
  h->Draw("HIST");
  for(std::vector<TH1*>::reverse_iterator rit = spectra.rbegin(); rit != spectra.rend(); ++rit) {
    (*rit)->Draw("HISTsame");
  }
  label->Draw("same");
  leg->Draw("same");
  c->SetLogy();
  gPad->RedrawAxis();
  util::FileOps::toFiles(c,OUT_FILE);
}


// Plot pthat spectrum before and after weighting
// ------------------------------------------------------------
void plotPtHat(const TString &fileName, const TString &treeName, const TString &outNamePrefix) {
  TH1* hPtHat = new TH1D("hPtHat","",500,0,3000);
  util::HistOps::setAxisTitles(hPtHat,"#hat{p}_{T}","GeV","Probability Density");
  hPtHat->SetLineStyle(2);
  hPtHat->SetLineWidth(2);
  hPtHat->SetLineColor(kBlue);

  TH1* hPtHatWeighted = static_cast<TH1*>(hPtHat->Clone("hPtHatWeighted"));
  hPtHatWeighted->SetLineStyle(1);
  hPtHatWeighted->SetLineColor(kRed);

  TChain* chain = new TChain(treeName,treeName);
  chain->Add(fileName);
  chain->SetBranchStatus("*",0);
  chain->SetBranchStatus("GenEvtScale",1);
  chain->SetBranchStatus("Weight",1);
  Float_t         GenEvtScale = 0.;
  Float_t         Weight = 0.;
  chain->SetBranchAddress("GenEvtScale",&GenEvtScale);
  chain->SetBranchAddress("Weight",&Weight);

  int nEntries = chain->GetEntries();
  int nEntriesNearlyDone = static_cast<int>(0.95*nEntries);
  int nEntriesCool = static_cast<int>(0.7*nEntries);
  for(int i = 0; i < nEntries; ++i) {
    chain->GetEntry(i);

    if( GenEvtScale > 0. && Weight > 0. ) {
      hPtHat->Fill(GenEvtScale);
      hPtHatWeighted->Fill(GenEvtScale,Weight);
    }

    if( i % 500000 == 0 ) {
      if( i > nEntriesNearlyDone ) 
	std::cout << "I processed already " << i << " entries - almost finished :D" << std::endl;
      else if( i > nEntriesCool ) 
	std::cout << "I processed already " << i << " entries - am I the man or what?!" << std::endl;
      else
	std::cout << "I processed already " << i << " entries." << std::endl;
    }
  }

  TPaveText* label = util::LabelFactory::createPaveText(1,-0.6);
  label->AddText(util::LabelFactory::mc());

  TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,-0.4,1);
  leg->AddEntry(hPtHat,"Generated","L");
  leg->AddEntry(hPtHatWeighted,"QCD","L");

  util::HistOps::normHist(hPtHat,"weight");
  util::HistOps::normHist(hPtHatWeighted,"weight");
  util::HistOps::setYRange(hPtHat,1,3E-18);
  hPtHat->GetXaxis()->SetNdivisions(505);
  hPtHat->GetXaxis()->SetRangeUser(20.,2850.);
  hPtHat->GetXaxis()->SetMoreLogLabels();
  hPtHat->GetXaxis()->SetNoExponent();

  TCanvas* c = new TCanvas(outNamePrefix+"_PtHatSpectrum","PtHat spectra",500,500);
  c->cd();
  hPtHat->Draw("HIST");
  hPtHatWeighted->Draw("HISTsame");
  label->Draw("same");
  leg->Draw("same");
  c->SetLogy();
  gPad->RedrawAxis();
  util::FileOps::toFiles(c,OUT_FILE);

  TCanvas* c2 = new TCanvas(outNamePrefix+"_PtHatSpectrum_Logx","PtHat spectra",500,500);
  c2->cd();
  hPtHat->Draw("HIST");
  hPtHatWeighted->Draw("HISTsame");
  label->Draw("same");
  leg->Draw("same");
  c2->SetLogy();
  c2->SetLogx();
  gPad->RedrawAxis();
  util::FileOps::toFiles(c2,OUT_FILE);
}


// ------------------------------------------------------------
void plotControlDistributions() {
  gErrorIgnoreLevel = 1001;
  util::StyleSettings::setStyleNoteNoTitle();

  const TString outNamePrefix = "ControlPlots_163337-167151_Summer11-PythiaZ2";
  const TString fileNameData = "../results/Analysis2011/ControlPlots/Kalibri_ControlPlots_Run2011A-163337-167151.root";
  const TString fileNameMC = "../results/Analysis2011/ControlPlots/Kalibri_ControlPlots_MCSummer11_AK5PF.root";
  const sampleTools::BinningAdmin binAdmin("config/Analysis2011/Binning/BinningAdmin2011_v2.cfg");

  const TString ptRelVar = util::LabelFactory::pt3Rel();
  const TString ptGenRelVar = util::LabelFactory::pt3RelGen();
  const TString ptRelLabel = util::LabelFactory::pt3RelCut(0.14);
  const TString imbalGen = "#alpha_{imbal}";
  const TString deltaPhiVar = util::LabelFactory::deltaPhi();
  const TString deltaPhiGenVar = util::LabelFactory::deltaPhiGen();
  const TString deltaPhiLabel = util::LabelFactory::deltaPhiCut(2.7);
  const TString deltaPhiGenLabel = util::LabelFactory::deltaPhiGenCut(2.7);

  OUT_FILE = new TFile(outNamePrefix+".root","RECREATE");


  plotNMinus1Eta(fileNameData,fileNameMC,binAdmin,3,3,-5.1,-5.1,false,"#eta","","jets",deltaPhiLabel+",  "+ptRelLabel,outNamePrefix+"_NMin1Eta");

  for(int ptBin = 0; ptBin < 16; ++ptBin) {
//     // N-1 plots
    plotNMinus1(fileNameData,fileNameMC,"hDeltaPhiNMin1",binAdmin,0,ptBin,2,2.1,3.5,false,deltaPhiVar,"","events",ptRelLabel,outNamePrefix+"_NMin1DeltaPhi");
    plotNMinus1(fileNameData,fileNameMC,"hPtRelNMin1",binAdmin,0,ptBin,3,0.,1.,false,ptRelVar,"","events",deltaPhiLabel,outNamePrefix+"_NMin1PtRel");

    // Pt3 correlations
    plotCorrelations(true,fileNameData,binAdmin,"hDeltaPhiVsPt3Rel",0,ptBin,-1,0.,1.,1.95,3.5,false,ptRelVar,deltaPhiVar,deltaPhiLabel,outNamePrefix+"_DeltaPhiVsPt3Rel");
    plotCorrelations(false,fileNameMC,binAdmin,"hPt3RelGenVsReco",0,ptBin,-1,0.,0.35,0.,0.35,false,"#alpha","#alpha^{gen}",deltaPhiLabel,outNamePrefix+"_Pt3RelGenVsReco");
    plotCorrelations(false,fileNameMC,binAdmin,"hPt3RelGenVsImbalGen",0,ptBin,-1,0.,0.35,0.,0.35,false,"#alpha^{imbal}","#alpha^{gen}",deltaPhiGenLabel,outNamePrefix+"_Pt3RelGenVsImbalGen");
  }


  // Spectra
  std::vector<TString> fileNames;
  fileNames.push_back("~/Kalibri/input/Kalibri_DijetSpectra_PythiaZ2_Summer11.root");
  fileNames.push_back("~/Kalibri/input/Kalibri_DijetSpectra_Herwigpp23_Summer11.root");
  std::vector<TString> sampleNames;
  sampleNames.push_back("Pythia");
  sampleNames.push_back("Herwig++");
  plotUnderlyingSpectra(fileNames.front(),binAdmin,outNamePrefix);
  plotVariedSpectrum(fileNames,sampleNames,binAdmin,outNamePrefix);
  plotPtHat("~/lustre/mc/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6job_*_ak5FastPF.root","DiJetTree",outNamePrefix);
  plotExpectedDataPileUpDistribution("~/PileUp/Pileup_2011_to_172255.root","pileup",outNamePrefix);

  OUT_FILE->Close();
  delete OUT_FILE;
}

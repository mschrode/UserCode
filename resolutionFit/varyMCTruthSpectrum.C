// $Id: $

#define UTILS_AS_HEADER_FILE

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"

#include "../util/FileOps.h"
#include "../util/StyleSettings.h"
#include "../util/utils.h"


unsigned int N_VAR = 0;


TH1* vary(const TString &outName, const TH1* hOrig, double expo) {
  // Create clone to store in file
  TH1* hVar = static_cast<TH1*>(hOrig->Clone("hPtGen"));
  hVar->Reset();

  // Vary spectrum
  double x0 = 30.;
  for(int bin = 1; bin <= hOrig->GetNbinsX(); ++bin) {
    double x = hOrig->GetBinCenter(bin);
    double w = pow(x/x0,expo);
    hVar->SetBinContent(bin,w*(hOrig->GetBinContent(bin)));
    hVar->SetBinError(bin,hOrig->GetBinError(bin));
  }    
  if( hVar->Integral() ) hVar->Scale(1./hVar->Integral("width"));

  // Store varied spectrum in file
  TFile outFile(outName,"RECREATE");
  outFile.WriteTObject(hVar);
  outFile.Close();

  // Change name of varied spectrum and return it
  hVar->SetName("hVar"+util::toTString(++N_VAR));

  return hVar;
}

void varyMCTruthSpectrum(const TString &name) {
  util::StyleSettings::paperNoTitle();

  TH1* hSpec = util::FileOps::readTH1(name,"hPtGen","hSpecOrig");
  hSpec->UseCurrentStyle();
  hSpec->SetLineWidth(1);

  TH1* hVarUp = vary(name(0,name.Last('.'))+"_Up.root",hSpec,0.5);
  hVarUp->SetLineColor(kRed);
  TH1* hVarDown = vary(name(0,name.Last('.'))+"_Down.root",hSpec,-0.5);
  hVarDown->SetLineColor(kBlue);

  TPaveText* label = util::LabelFactory::createPaveText(2,0.7);
  label->AddText("CMS Simulation, #sqrt{s} = 7 TeV");
  label->AddText("Anti-k_{T} (d=0.5) Gen Jets");

  TLegend* leg = util::LabelFactory::createLegendColWithOffset(3,0.55,2);
  leg->AddEntry(hSpec,"Nominal","L");
  leg->AddEntry(hVarDown,"Variation down","L");
  leg->AddEntry(hVarUp,"Variation up","L");

  TCanvas* can2 = new TCanvas("can2","Varied Spectra",500,500);
  can2->cd();
  util::HistOps::setAxisTitles(hSpec,"p^{gen}_{T}","GeV","jets",true);
  hSpec->SetTitle("");
  hSpec->GetYaxis()->SetRangeUser(3E-15,1.);
  hSpec->Draw("HIST");
  hVarUp->Draw("HISTsame");  
  hVarDown->Draw("HISTsame");  
  hSpec->Draw("HISTsame");
  label->Draw("same");
  leg->Draw("same");
  can2->SetLogy();
  can2->SaveAs("VariedDijetSpectrum_Eta0.eps","eps");
}



void plotSpectra() {
  util::StyleSettings::paperNoTitle();

  std::vector<TString> fileNames;
  fileNames.push_back("~/Kalibri_new/input/Kalibri_DijetSpectrum_Pt0020-1500_Eta00-11.root");
  fileNames.push_back("~/Kalibri_new/input/Kalibri_DijetSpectrum_Pt0020-1500_Eta11-17.root");
  fileNames.push_back("~/Kalibri_new/input/Kalibri_DijetSpectrum_Pt0020-1500_Eta17-23.root");
  fileNames.push_back("~/Kalibri_new/input/Kalibri_DijetSpectrum_Pt0020-1500_Eta23-50.root");
  util::HistVec hPtGen = util::FileOps::readTH1(fileNames,"hPtGen");

  std::vector<double> etaBins;
  etaBins.push_back(0.);
  etaBins.push_back(1.1);
  etaBins.push_back(1.7);
  etaBins.push_back(2.3);
  etaBins.push_back(5.0);

  TLegend* leg = util::LabelFactory::createLegendCol(hPtGen.size(),0.4);
  for(unsigned int i = 0; i < hPtGen.size(); ++i) {
    hPtGen[i]->UseCurrentStyle();
    hPtGen[i]->SetTitle("");
    hPtGen[i]->GetXaxis()->SetRangeUser(20.,1500.);
    hPtGen[i]->GetXaxis()->SetMoreLogLabels();
    hPtGen[i]->GetXaxis()->SetTitle("p^{gen}_{T} (GeV)");
    hPtGen[i]->GetYaxis()->SetRangeUser(3E-14,8.);
    hPtGen[i]->GetYaxis()->SetTitle("Probability");
    hPtGen[i]->SetLineWidth(1);
    hPtGen[i]->SetLineColor(util::StyleSettings::color(i));

    leg->AddEntry(hPtGen[i],util::toTString(etaBins.at(i))+" < |#eta| < "+util::toTString(etaBins.at(i+1)),"L");
  }

  TPaveText* label = util::LabelFactory::createPaveText(2,-0.55);
  label->AddText("CMS Simulation, #sqrt{s} = 7 TeV");
  label->AddText("Anti-k_{T} (d=0.5) Gen Jets");

  TCanvas* can = new TCanvas("can","Dijet Spectrum",500,500);
  can->cd();
  hPtGen.back()->Draw("HISTL");
  for(int i = static_cast<int>(hPtGen.size()-2); i >= 0;  i--) {
    hPtGen[i]->Draw("HISTLsame");
  }
  label->Draw("same");
  leg->Draw("same");
  can->SetLogy();
  can->SaveAs("PtGenSpectra.eps","eps");  
}

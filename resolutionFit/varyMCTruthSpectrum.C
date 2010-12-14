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
  util::StyleSettings::presentationNoTitle();

  TH1* hSpec = util::FileOps::readTH1(name,"hPtGen","hSpecOrig");

  TH1* hVarUp = vary(name(0,name.Last('.'))+"_Up.root",hSpec,0.5);
  hVarUp->SetLineColor(kRed);
  TH1* hVarDown = vary(name(0,name.Last('.'))+"_Down.root",hSpec,-0.5);
  hVarDown->SetLineColor(kBlue);

  TCanvas* can = new TCanvas("can","Varied Spectra",500,500);
  can->cd();
  util::HistOps::setAxisTitles(hSpec,"p^{gen}_{T}","GeV","jets",true);
  hSpec->SetTitle("");
  hSpec->GetYaxis()->SetRangeUser(3E-15,1.);
  hSpec->Draw("HIST");
  hVarUp->Draw("HISTsame");  
  hVarDown->Draw("HISTsame");  
  can->SetLogy();
}

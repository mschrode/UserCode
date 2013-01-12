// Plot the lepton ptGen spectra
void plotHadTau1(const TString &fileName = "HadTau_WJetMC.root") {
  gROOT->ProcessLine(".L StyleMatters.h+");
  StyleMatters::init();

  // Get histograms from file
  TH1* hGenMuPt = 0;
  TH1* hGenTauPt = 0;
  TFile file(fileName,"READ");
  file.GetObject("hGenMuPt",hGenMuPt);
  file.GetObject("hGenTauPt",hGenTauPt);
  if( !(hGenMuPt && hGenTauPt) ) {
    std::cerr << "ERROR: Histograms not found" << std::endl;
    exit(-1);
  }
  hGenMuPt->SetDirectory(0);
  hGenMuPt->UseCurrentStyle();
  hGenTauPt->SetDirectory(0);
  hGenTauPt->UseCurrentStyle();
  file.Close();


  // Create ratio
  TH1* hRatio = static_cast<TH1*>(hGenTauPt->Clone("hRatio"));
  hRatio->Divide(hGenMuPt);
  hRatio->GetXaxis()->SetTitle("p_{T} [GeV]");
  hRatio->GetYaxis()->SetTitle("N(#tau #rightarrow had) / N(#mu)");


  // Set style
  hGenMuPt->GetXaxis()->SetTitle("p_{T} [GeV]");
  hGenMuPt->GetYaxis()->SetTitle("N");
  hGenMuPt->SetMarkerStyle(20);
  hGenMuPt->SetMarkerColor(kBlue);
  hGenMuPt->SetLineColor(hGenMuPt->GetMarkerColor());

  hGenTauPt->GetXaxis()->SetTitle("p_{T} [GeV]");
  hGenTauPt->GetYaxis()->SetTitle("N");
  hGenTauPt->SetMarkerStyle(21);
  hGenTauPt->SetMarkerColor(kRed);
  hGenTauPt->SetLineColor(hGenTauPt->GetMarkerColor());

  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetLineColor(hRatio->GetMarkerColor());

  
  // Create legend
  TLegend* leg = new TLegend(0.6,0.75,0.9,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->AddEntry(hGenMuPt," #mu","P");
  leg->AddEntry(hGenTauPt," #tau #rightarrow had","P");


  // Draw histograms
  TCanvas* canSpecs = new TCanvas("canSpecs","Pt Spectra",600,600);
  canSpecs->cd();
  hGenMuPt->Draw("PE1");
  hGenTauPt->Draw("PE1same");
  leg->Draw("same");
  gPad->SetLogy();

  TCanvas* canRatio = new TCanvas("canRatio","Ratio",600,600);
  canRatio->cd();
  hRatio->Draw("PE1");
}

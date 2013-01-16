// Plot the tau-response templates
void plotHadTau2(const TString &fileName = "HadTau_WJetMC.root") {
  gROOT->ProcessLine(".L ../utils/StyleMatters.h+");
  StyleMatters::init();

  
  // Get histograms from file
  const unsigned int kNBins = 4;
  TH1* hResp[kNBins];
  TFile file(fileName,"READ");
  for(unsigned int i = 0; i < kNBins; ++i) {
    TString name = "hTauResp_";
    name += i;
    file.GetObject(name,hResp[i]);
    if( !hResp[i] ) {
      std::cerr << "ERROR: Histograms not found" << std::endl;
      exit(-1);
    }
    hResp[i]->SetDirectory(0);
    hResp[i]->UseCurrentStyle();
  }
  file.Close();
  
  
  // Normalize the response distributions
  // to get probability density
  for(unsigned int i = 0; i < kNBins; ++i) {
    if( hResp[i]->Integral() ) hResp[i]->Scale(1./hResp[i]->Integral("width"));
  }
  

  // Set style
  for(unsigned int i = 0; i < kNBins; ++i) {
    char title[100];  
    sprintf(title,"Probability / %.2f",hResp[i]->GetXaxis()->GetBinWidth(1));
    hResp[i]->GetYaxis()->SetTitle(title);
    hResp[i]->SetLineColor(1+i);
    hResp[i]->GetYaxis()->SetRangeUser(0.,2.05);
  }


  // Create legend
  TLegend* leg = new TLegend(0.63,0.55,0.9,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader("p_{T}(#tau^{gen}) / GeV");
  leg->AddEntry(hResp[0],"20 - 30","L");
  leg->AddEntry(hResp[1],"30 - 50","L");
  leg->AddEntry(hResp[2],"50 - 100","L");
  leg->AddEntry(hResp[3],"> 100","L");



  // Draw
  TCanvas* can = new TCanvas("can","Response",600,600);
  can->cd();
  hResp[0]->Draw("HIST");
  for(unsigned int i = 1; i < kNBins; ++i) {
    hResp[i]->Draw("HISTsame");
  }
  leg->Draw("same");
}

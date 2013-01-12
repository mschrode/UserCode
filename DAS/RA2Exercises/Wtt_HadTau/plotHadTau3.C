// Plot the closure test
void plotHadTau3(double scale = 1.,
		 const TString &fileName = "HadTau_WJetMC_PredGen.root") {
  gROOT->ProcessLine(".L StyleMatters.h+");
  StyleMatters::init();

  
  // Get histograms from file
  const unsigned int kNDists = 3;
  TH1* hTrue[kNDists];
  TH1* hPred[kNDists];
  TFile file(fileName,"READ");
  for(unsigned int i = 0; i < kNDists; ++i) {
    TString name = "";
    if( i == 0 )      name = "TauJetPt";
    else if( i == 1 ) name = "Ht";
    else if( i == 2 ) name = "Mht";
    file.GetObject("hTrue"+name,hTrue[i]);
    file.GetObject("hPred"+name,hPred[i]);
    if( !(hTrue[i] && hPred[i]) ) {
      std::cerr << "ERROR: Histograms not found" << std::endl;
      exit(-1);
    }
    hTrue[i]->SetDirectory(0);
    hTrue[i]->UseCurrentStyle();
    hPred[i]->SetDirectory(0);
    hPred[i]->UseCurrentStyle();
  }
  file.Close();
  
  
  // Apply correction factors
  for(unsigned int i = 0; i < kNDists; ++i) {
    hPred[i]->Scale(scale);
  }


  // Set style
  for(unsigned int i = 0; i < kNDists; ++i) {
    TString title = "";
    if( i == 0 )      title = "p_{T}(#tau) [GeV]";
    else if( i == 1 ) title = "H_{T} [GeV]";
    else if( i == 2 ) title = "#slash{H}_{T} [GeV]";
    hTrue[i]->GetXaxis()->SetTitle(title);
    hPred[i]->GetXaxis()->SetTitle(title);

    hTrue[i]->SetLineColor(kBlue);

    hPred[i]->SetMarkerStyle(20);
    hPred[i]->SetMarkerColor(kRed);
    hPred[i]->SetLineColor(hPred[i]->GetMarkerColor());
  }


  // Create legend
  TLegend* leg = new TLegend(0.4,0.75,0.9,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->AddEntry(hTrue[0],"MC Expectation");
  if( fileName.Contains("Gen") ) leg->AddEntry(hPred[0],"Gen-Based Pred.");
  else if( fileName.Contains("Reco") ) leg->AddEntry(hPred[0],"Data-Based Pred.");
  else leg->AddEntry(hPred[0],"Prediction");


  // Draw
  for(unsigned int i = 0; i < kNDists; ++i) {
    TString name = "";
    if( i == 0 )      name = "TauJetPt";
    else if( i == 1 ) name = "HT";
    else if( i == 2 ) name = "MHT";

    TCanvas* can = new TCanvas(name,name,600,600);
    can->cd();
    hTrue[i]->Draw("HISTE");
    hPred[i]->Draw("PE1same");
    leg->Draw("same");
    gPad->SetLogy();

    if( fileName.Contains("Gen") ) name = "hGenClosure"+name;
    else if( fileName.Contains("Reco") ) name = "hRecoClosure"+name;
    else name = "hClosure"+name;
    can->SaveAs(name+".eps","eps");
  }
}


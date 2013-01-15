// Plot the closure test
void plotHadTau3(double scale = 1.,
		 const TString &fileName = "HadTau_WJetMC_PredGen.root") {
  gROOT->ProcessLine(".L ../utils/StyleMatters.h+");
  StyleMatters::init();

  bool isMCPred = false;
  if( fileName.Contains("MC") ) isMCPred = true;
  bool isGenPred = false;
  if( isMCPred && fileName.Contains("Gen") ) isGenPred = true;

  
  // Get histograms from file
  const unsigned int kNDists = 3;
  TH1* hTrue[kNDists];
  TH1* hPred[kNDists];
  TH1* hMuonPt = 0;
  TFile file(fileName,"READ");
  for(unsigned int i = 0; i < kNDists; ++i) {
    TString name = "";
    if( i == 0 )      name = "TauJetPt";
    else if( i == 1 ) name = "Ht";
    else if( i == 2 ) name = "Mht";
    file.GetObject("hPred"+name,hPred[i]);
    if( !hPred[i] ) {
      std::cerr << "ERROR: Histograms not found" << std::endl;
      exit(-1);
    }
    hPred[i]->SetDirectory(0);
    hPred[i]->UseCurrentStyle();
    if( isMCPred ) {
      file.GetObject("hTrue"+name,hTrue[i]);
      if( !hTrue[i] ) {
	std::cerr << "ERROR: Histograms not found" << std::endl;
	exit(-1);
      }
      hTrue[i]->SetDirectory(0);
      hTrue[i]->UseCurrentStyle();
    }
  }
  file.GetObject("hMuonPt",hMuonPt);
  if( !hMuonPt ) {
    std::cerr << "ERROR: Histogram not found" << std::endl;
    exit(-1);
  }
  hMuonPt->SetDirectory(0);
  hMuonPt->UseCurrentStyle();
  file.Close();
  
  
  // Apply correction factors
  for(unsigned int i = 0; i < kNDists; ++i) {
    hPred[i]->Scale(scale);
  }


  // Set style
  for(unsigned int i = 0; i < kNDists; ++i) {
    TString xTitle = "";
    if( i == 0 )      xTitle = "p_{T}(#tau) [GeV]";
    else if( i == 1 ) xTitle = "H_{T} [GeV]";
    else if( i == 2 ) xTitle = "#slash{H}_{T} [GeV]";

    TString yTitle = "N(events) / ";
    yTitle += static_cast<int>(hPred[i]->GetYaxis()->GetBinWidth(1));
    yTitle += " GeV";

    hPred[i]->GetXaxis()->SetTitle(xTitle);
    hPred[i]->GetYaxis()->SetTitle(yTitle);
    hPred[i]->SetMarkerStyle(20);
    hPred[i]->SetMarkerColor(kRed);
    hPred[i]->SetLineColor(hPred[i]->GetMarkerColor());

    if( isMCPred ) {
      hTrue[i]->GetXaxis()->SetTitle(xTitle);
      hTrue[i]->GetYaxis()->SetTitle(yTitle);
      hTrue[i]->SetLineColor(kBlue);
    }
  }
  if( isGenPred ) hMuonPt->GetXaxis()->SetTitle("p_{T}(#mu^{gen}) [GeV]");
  else            hMuonPt->GetXaxis()->SetTitle("p_{T}(#mu) [GeV]");
  hMuonPt->SetMarkerStyle(20);
  hMuonPt->SetMarkerColor(kBlack);
  hMuonPt->SetLineColor(hMuonPt->GetMarkerColor());


  // Create legend
  TLegend* leg = new TLegend(0.4,0.75,0.9,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  if( isMCPred ) {
    leg->AddEntry(hTrue[0],"MC Expectation");
    if( isGenPred ) leg->AddEntry(hPred[0],"Gen-Based Pred.");
    else            leg->AddEntry(hPred[0],"Data-Based Pred.");
  }


  // Draw
  for(unsigned int i = 0; i < kNDists; ++i) {
    TString name = "";
    if( i == 0 )      name = "TauJetPt";
    else if( i == 1 ) name = "HT";
    else if( i == 2 ) name = "MHT";

    TCanvas* can = new TCanvas(name,name,600,600);
    can->cd();
    if( isMCPred ) {
      hTrue[i]->Draw("HISTE");
      hPred[i]->Draw("PE1same");
      leg->Draw("same");
      if( isGenPred ) name = "hGenClosure"+name;
      else            name = "hRecoClosure"+name;
    } else {
      hPred[i]->Draw("PE1");
      name = "hDataPred"+name;
    }
    gPad->SetLogy();
    can->SaveAs(name+".eps","eps");
  }

  TCanvas* can = new TCanvas("can","muon pt",600,600);
  can->cd();
  hMuonPt->Draw("PE1");
  gPad->SetLogy();
  TString name = "MuonPt";
  if( isMCPred ) {
    if( isGenPred ) name = "hGenClosure"+name;
    else            name = "hRecoClosure"+name;
  } else {
    name = "hDataPred"+name;
  }
  can->SaveAs(name+".eps","eps");
}


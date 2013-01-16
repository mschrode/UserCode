void plotGeneral1(int sampleId) {
  gROOT->ProcessLine(".L ../utils/StyleMatters.h+");
  StyleMatters::init();

  std::cout << "Plotting distributions for " << sampleLabel(sampleId) << std::endl;

 
  // Get histograms from file
  TFile file("General_"+fileName(sampleId)+".root","READ");

  TH1* hNJets = 0;
  TH1* hHt = 0;
  TH1* hMht = 0;
  TH1* hMEff = 0;
  file.GetObject("hNJets",hNJets);
  file.GetObject("hHt",hHt);
  file.GetObject("hMht",hMht);
  file.GetObject("hMEff",hMEff);
  if( !(hNJets && hHt && hMht && hMEff) ) {
    std::cerr << "ERROR reading histograms" << std::endl;
    exit(-1);
  }
  hNJets->SetDirectory(0);
  hHt->SetDirectory(0);
  hMht->SetDirectory(0);
  hMEff->SetDirectory(0);
  hNJets->UseCurrentStyle();
  hHt->UseCurrentStyle();
  hMht->UseCurrentStyle();
  hMEff->UseCurrentStyle();
  if( sampleId == 0 ) {
    hNJets->SetMarkerStyle(20);
    hHt->SetMarkerStyle(20);
    hMht->SetMarkerStyle(20);
    hMEff->SetMarkerStyle(20);
  }

  const int kNJetHists = 3;
  TH1* hJetPt[kNJetHists];
  TH1* hJetPhi[kNJetHists];
  TH1* hJetEta[kNJetHists];
  for(unsigned int i = 0; i < kNJetHists; ++i) {
    hJetPt[i] = 0;
    hJetPhi[i] = 0;
    hJetEta[i] = 0;

    TString name = "hJetPt_";
    name += i;
    file.GetObject(name,hJetPt[i]);
    name = "hJetPhi_";
    name += i;
    file.GetObject(name,hJetPhi[i]);
    name = "hJetEta_";
    name += i;
    file.GetObject(name,hJetEta[i]);
    if( !(hJetPt[i] && hJetEta[i] && hJetPhi[i]) ) {
      std::cerr << "ERROR reading histograms" << std::endl;
      exit(-1);
    }
    hJetPt[i]->SetDirectory(0);
    hJetPt[i]->UseCurrentStyle();
    hJetEta[i]->SetDirectory(0);
    hJetEta[i]->UseCurrentStyle();
    hJetPhi[i]->SetDirectory(0);
    hJetPhi[i]->UseCurrentStyle();
    if( sampleId == 0 ) {
      hJetPt[i]->SetMarkerStyle(20);
      hJetEta[i]->SetMarkerStyle(20);
      hJetPhi[i]->SetMarkerStyle(20);
    }
  }
  file.Close();


  
  // Draw
  TString drawOption = "HISTE";
  if( sampleId == 0 ) drawOption = "PE1"; // Markers for data

  TString outName = "h"+fileName(sampleId)+"_";

  // Create legend
  TLegend* leg = new TLegend(0.3,0.8,0.9,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  char entry[100];
  sprintf(entry,"N(total) = %.1f",hHt->Integral());
  leg->AddEntry(hHt,entry,"");

  TCanvas* canNJets = new TCanvas("canNJets","NJets",600,600);
  canNJets->cd();
  hNJets->Draw(drawOption);
  leg->Draw("same");
  canNJets->SaveAs(outName+"NJets.eps","eps");

  TCanvas* canHt = new TCanvas("canHt","Ht",600,600);
  canHt->cd();
  hHt->Draw(drawOption);
  leg->Draw("same");
  canHt->SaveAs(outName+"Ht.eps","eps");

  TCanvas* canMht = new TCanvas("canMht","Mht",600,600);
  canMht->cd();
  hMht->Draw(drawOption);
  leg->Draw("same");
  canMht->SaveAs(outName+"Mht.eps","eps");

  for(unsigned int i = 0; i < kNJetHists; ++i) {
    TString name = "JetPt";
    name += i+1;
    TCanvas* can = new TCanvas(name,name,600,600);
    can->cd();
    hJetPt[i]->Draw(drawOption);
    can->SaveAs(outName+name+".eps","eps");
  }

  for(unsigned int i = 0; i < kNJetHists; ++i) {
    TString name = "JetPhi";
    name += i+1;
    TCanvas* can = new TCanvas(name,name,600,600);
    can->cd();
    hJetPhi[i]->Draw(drawOption);
    can->SaveAs(outName+name+".eps","eps");
  }

  for(unsigned int i = 0; i < kNJetHists; ++i) {
    TString name = "JetEta";
    name += i+1;
    TCanvas* can = new TCanvas(name,name,600,600);
    can->cd();
    hJetEta[i]->Draw(drawOption);
    can->SaveAs(outName+name+".eps","eps");
  }

}



// === Implementation of Auxiliary Functions =====================


// Return the label for a given sample
TString sampleLabel(int sampleId) {
  TString label = "";
  if( sampleId == 0 )      label += "Data";
  else if( sampleId == 1 ) label += "QCD";
  else if( sampleId == 2 ) label += "t#bar{t}+Jets";
  else if( sampleId == 3 ) label += "W(l#nu)+Jets";
  else if( sampleId == 4 ) label += "Z(#nu#bar{#nu})+Jets";
  else if( sampleId == 5 ) label += "LM6";
  else if( sampleId == 6 ) label += "LM9";
  else {
    std::cerr << "ERROR: no sample with id " << sampleId << std::endl;
    exit(-1);
  }

  return label;
}


// Return the file name for a given sample
TString fileName(int sampleId) {
  TString name = "";
  if( sampleId == 0 )      name += "Data";
  else if( sampleId == 1 ) name += "QCD";
  else if( sampleId == 2 ) name += "TTJets";
  else if( sampleId == 3 ) name += "WJets";
  else if( sampleId == 4 ) name += "ZInv";
  else if( sampleId == 5 ) name += "LM6";
  else if( sampleId == 6 ) name += "LM9";
  else {
    std::cerr << "ERROR: no sample with id " << sampleId << std::endl;
    exit(-1);
  }

  return name;
}

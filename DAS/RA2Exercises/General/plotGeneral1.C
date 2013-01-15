void plotGeneral1(int sampleId) {
  gROOT->ProcessLine(".L ../StyleMatters.h+");
  StyleMatters::init();

  std::cout << "Plotting distributions for " << sampleLabel(sampleId) << std::endl;

 
  // Get histograms from file
  TFile file(fileName(sampleId)+".root","READ");

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
  }
  file.Close();


  TString drawOption = "HISTE";
  if( sampleId == 0 ) drawOption = "PE1"; // Markers for data

  // Draw
  TCanvas* canNJets = new TCanvas("canNJets","NJets",600,600);
  canNJets->cd();
  hNJets->Draw(drawOption);

  TCanvas* canHt = new TCanvas("canHt","Ht",600,600);
  canHt->cd();
  hHt->Draw(drawOption);
  std::cout << "Int(Ht) = " << hHt->Integral() << std::endl;

  TCanvas* canMht = new TCanvas("canMht","Mht",600,600);
  canMht->cd();
  hMht->Draw(drawOption);

  for(unsigned int i = 0; i < kNJetHists; ++i) {
    TString name = "Jet Pt ";
    name += i+1;
    TCanvas* can = new TCanvas(name,name,600,600);
    can->cd();
    hJetPt[i]->Draw(drawOption);
  }

  for(unsigned int i = 0; i < kNJetHists; ++i) {
    TString name = "Jet Phi ";
    name += i+1;
    TCanvas* can = new TCanvas(name,name,600,600);
    can->cd();
    hJetPhi[i]->Draw(drawOption);
  }

  for(unsigned int i = 0; i < kNJetHists; ++i) {
    TString name = "Jet Eta ";
    name += i+1;
    TCanvas* can = new TCanvas(name,name,600,600);
    can->cd();
    hJetEta[i]->Draw(drawOption);
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
  else {
    std::cerr << "ERROR: no sample with id " << sampleId << std::endl;
    exit(-1);
  }

  return label;
}


// Return the file name for a given sample
TString fileName(int sampleId) {
  TString name = "General_";
  if( sampleId == 0 )      name += "QCDFlat_test";
  else if( sampleId == 1 ) name += "QCDFlat_test";
  else if( sampleId == 2 ) name += "TTJets";
  else if( sampleId == 3 ) name += "WJets";
  else if( sampleId == 4 ) name += "ZInv";
  else if( sampleId == 5 ) name += "LM6";
  else {
    std::cerr << "ERROR: no sample with id " << sampleId << std::endl;
    exit(-1);
  }

  return name;
}

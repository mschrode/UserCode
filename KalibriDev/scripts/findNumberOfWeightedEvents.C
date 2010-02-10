#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TChain.h"
#include "TH1F.h"



void findNumberOfWeightedEvents(const std::string& fileList, float min, float max, float expo) {
  // Opening root files in input list
  std::cout << "Opening files in list '" << fileList << "'\n";
  TChain* chain = new TChain("DiJetTree"); 
  std::ifstream file;
  file.open(fileList.c_str());
  int nOpenedFiles = 0;
  std::string name = "";
  while( !file.eof() ) {
    file >> name;
    if( file.eof() ) break;
    chain->Add( name.c_str() );
    nOpenedFiles++;
  }
  file.close();
  std::cout << "Opened " << nOpenedFiles << " files\n";

  // Creating histograms of pthat
  std::cout << "Filling histograms of ptHat\n";
  TH1F *hPtHat[2];
  hPtHat[0] = new TH1F("hPtHat0",";#hat{p}_{T} GeV;Number of events",
		       500,min,max);
  hPtHat[1] = static_cast<TH1F*>(hPtHat[0]->Clone("hPtHat1"));
  hPtHat[1]->SetLineColor(4);

  float ptHat = 0.;
  chain->SetBranchAddress("GenEvtScale",&ptHat);
  for(int n = 0; n < chain->GetEntries(); n++) {
    chain->GetEntry(n);
    hPtHat[0]->Fill(ptHat);
    hPtHat[1]->Fill(ptHat,pow(ptHat,expo));
  }
  
  // Output
  std::cout << "Processed " << chain->GetEntries() << " events.\n";
  std::cout << "Number of events\n";
  std::cout << " unweighted: " << hPtHat[0]->Integral() << std::endl;
  std::cout << " weighted:   " << hPtHat[1]->Integral() << std::endl;

  TCanvas *can = new TCanvas("can","PtHat spectra",1000,500);
  can->Divide(2,1);
  can->cd(1);
  hPtHat[0]->Draw();
  can->cd(2);
  hPtHat[1]->Draw();
}

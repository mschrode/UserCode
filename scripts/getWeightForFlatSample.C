// $Id: $

#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TH1D.h"
#include "TString.h"

#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


void plot(const TString &name, int nEvts = -1) {

  TH1 *hFlat = util::HistOps::createTH1D("hFlat",100,0.,3000.,"#hat{p}_{T}","GeV","events");
  TH1 *hWeighted = util::HistOps::createTH1D("hWeighted",100,0.,3000.,"#hat{p}_{T}","GeV","events");
  TH1 *hQCD = util::HistOps::createTH1D("hQCD",100,0.,3000.,"#hat{p}_{T}","GeV","events");

  TChain *chain = new TChain("DiJetTree");
  chain->Add(name);

  float ptHat = 0.;
  float weight = 0.;
  chain->SetBranchAddress("GenEvtScale",&ptHat);
  chain->SetBranchAddress("Weight",&weight);
  
  if( nEvts < 0 || nEvts > chain->GetEntries() ) nEvts = chain->GetEntries();
  for(int n = 0; n < nEvts; ++n) {
    if( n%50000 == 0 ) std::cout << " Entry " << n << std::endl;
    chain->GetEntry(n);

    hFlat->Fill(ptHat);
    hWeighted->Fill(ptHat,weight);
    hQCD->Fill(ptHat,pow(static_cast<double>(ptHat),-4.5));
  }

  std::cout << "Flat    " << hFlat->Integral() << "  (" << chain->GetEntries() << ")" << std::endl;
  std::cout << "Weight  " << hWeighted->Integral() << std::endl;
  std::cout << "QCD     " << hQCD->Integral() << std::endl;

  TCanvas *can1 = new TCanvas("can1","Flat",500,500);
  can1->cd();
  hFlat->Draw();
  can1->SetLogy();

  TCanvas *can2 = new TCanvas("can2","Weighted (QCD)",500,500);
  can2->cd();
  hWeighted->Draw();
  can2->SetLogy();

  TCanvas *can3 = new TCanvas("can3","Weighted (RAW)",500,500);
  can3->cd();
  hQCD->Draw();
  can3->SetLogy();
}


void getWeightForFlatSample() {
  util::StyleSettings::presentationNoTitle();

  plot("/scratch/hh/current/cms/user/mschrode/mc/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10-START38_V12-v1B/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6job_*_ak5Calo.root");
  //  plot("~/scratch/test.root",-1);
}

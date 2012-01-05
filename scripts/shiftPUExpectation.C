// $Id: $
// Original author: Christian Sander
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/SusyAnalysis/RA2/QCDSmearProd/test/PUWeights.C?view=log

{
  const TString inputFile = "/afs/naf.desy.de/user/m/mschrode/PileUp/ExpectedPileUpDist_160404_180252_SelfCombined.root";
  TFile *d = TFile::Open(inputFile,"READ");
  
  d->cd();

  TH1F* htmp;
  gDirectory->GetObject("pileup;1", htmp);
  TH1F h(*htmp);

  //h->Draw();
  double Norm = h->Integral();
  TF1* fit = new TF1("f", "[1]*TMath::PoissonI(x+0.5,[0])", -0.5, 50.5);
  TF1* fitPlus = new TF1("f", "[1]*TMath::PoissonI(x+0.5,[0])", -0.5, 50.5);
  TF1* fitMinus = new TF1("f", "[1]*TMath::PoissonI(x+0.5,[0])", -0.5, 50.5);
  fit->SetParameters((double) h->GetMean(), Norm);
  fit->FixParameter(1, Norm);
  h.Fit(fit, "LRQN");
  double Mean = fit->GetParameter(0);

  fitMinus->SetParameters(Mean-1., Norm);
  fitPlus->SetParameters(Mean+1., Norm);

  TH1F pileup(h);
  pileup.Reset();
  pileup.Add(fit);

  TH1F pileupPlus(h);
  pileupPlus.Reset();
  pileupPlus.Add(fitPlus);

  TH1F pileupMinus(h);
  pileupMinus.Reset();
  pileupMinus.Add(fitMinus);

  //   d->Close();

  TString outputFile = inputFile;
  outputFile.ReplaceAll(".root","_Fit.root");
  TFile *d1 = TFile::Open(outputFile,"NEW");
  d1->cd();
  pileup.Write();
  //   d1->Close();

  outputFile = inputFile;
  outputFile.ReplaceAll(".root","_FitMeanDnOne.root");
  TFile *d2 = TFile::Open(outputFile,"NEW");
  d2->cd();
  pileupMinus.Write();
  //   d2->Close();

  outputFile = inputFile;
  outputFile.ReplaceAll(".root","_FitMeanUpOne.root");
  TFile *d3 = TFile::Open(outputFile,"NEW");
  d3->cd();
  pileupPlus.Write();
  //   d3->Close();

  d->Close();
  d1->Close();
  d2->Close();
  d3->Close();
}

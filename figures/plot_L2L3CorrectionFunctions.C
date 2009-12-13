#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TPostScript.h"
#include "TStyle.h"

//1.044     1.131        11       4    2000         1.3      -1.137      -1.237      0.6313     0.02118      0.9658       3.452       2.834       1.673
const int nParL2_ = 5;
const int nParL3_ = 4;

std::vector< std::vector<double> > parL2_;
std::vector<double> parL3_;


void readParameters(const std::string &fileName) {
  parL2_.clear();
  parL3_.clear();

  std::ifstream file;
  file.open(fileName.c_str());
  unsigned int line = 0;
  double       val  = -1.;
  while( !file.eof() ) {
    std::vector<double> par;        // Store parameters
    for(int i = 0; i < 2; i++) {    // Loop over eta edges of bin
      file >> val;
      par.push_back(val);
    }
    int n = 0;
    file >> n;
    if( n != 0 ) {
      for(int i = 0; i < 2; i++) {
	file >> val;
      }      
      for(int i = 0; i < nParL2_; i++) {    // Loop over parameters
	file >> val;
	par.push_back(val);
      }
      parL2_.push_back(par);
      for(int i = 0; i < nParL3_; i++) {    // Loop over parameters
	file >> val;
	if( static_cast<int>(parL3_.size()) < nParL3_ )
	  parL3_.push_back(val);
      }      
    }
    line++;
  }
  file.close();
}


void printParameters() {
  std::cout << "L2 parameters:\n";
  std::vector< std::vector<double> >::const_iterator etaBinIt = parL2_.begin();
  for(; etaBinIt != parL2_.end(); etaBinIt++) {
    std::vector<double>::const_iterator parIt = etaBinIt->begin();
    for(; parIt != etaBinIt->end(); parIt++) {
      std::cout << *parIt << "  ";
    }
    std::cout << "\n";
  }

  std::cout << "\nL3 parameters:\n";
  std::vector<double>::const_iterator parIt = parL3_.begin();
  for(; parIt != parL3_.end(); parIt++) {
    std::cout << *parIt << "  ";
  }
  std::cout << "\n";
}


double funcL2(double *x, double *par) {
  double pt = (x[0] < 4.0) ? 4.0 : (x[0] > 2000.0) ? 2000.0 : x[0];
  double logpt = log10(pt);
  return par[0]+logpt*(0.1 * par[1]+logpt *(0.01* par[2]+logpt*(0.01*par[3]+logpt*(0.01*par[4]))));
}

double funcL3(double *x, double *par) {
  double pt = (x[0] < 4.0) ? 4.0 : (x[0] > 2000.0) ? 2000.0 : x[0];
  double logpt = log10(pt);
  return par[0] + par[1]/(pow(logpt,par[2]) + par[3]);
}

void draw(double min, double max) {
  gStyle->SetOptStat(0);

  TPostScript *ps = new TPostScript("L2L3CorrectionFunctions.ps",111);
  TCanvas *can = new TCanvas("can","L2L3",500,500);

  TH1D *hFrame = new TH1D("hFrame","",100,min,max);
  hFrame->GetYaxis()->SetRangeUser(0,2);
  hFrame->GetYaxis()->SetTitle("Relative correction");
  hFrame->GetXaxis()->SetTitle("p_{T} (GeV)");

  std::vector< std::vector<double> >::const_iterator etaBinIt = parL2_.begin();
  for(; etaBinIt != parL2_.end(); etaBinIt++) {
    TF1 *fL2 = new TF1("fL2",funcL2,min,max,nParL2_);
    double etaMin = 0.; 
    double etaMax = 0.;
    int i = 0;
    std::vector<double>::const_iterator parIt = etaBinIt->begin();
    for(; parIt != etaBinIt->end(); parIt++, i++) {
      if( i == 0 ) etaMin = *parIt;
      else if( i == 1 ) etaMax = *parIt;
      else fL2->SetParameter(i-2,*parIt);
    }
    
    TPaveText * etaLabel = new TPaveText(0.34,0.82,0.80,0.93,"NDC");
    etaLabel->SetFillColor(0);
    etaLabel->SetTextFont(42);
    etaLabel->SetBorderSize(0);
    char title[50];
    sprintf(title,"%.2f < #eta < %.2f",etaMin,etaMax);
    etaLabel->AddText(title);

    can->cd();
    hFrame->Draw();
    fL2->Draw("same");
    etaLabel->Draw("same");
    can->SetLogx(1);
    can->Draw();
    ps->NewPage();

    delete fL2;
    delete etaLabel;
  }

  TF1 *fL3 = new TF1("fL3",funcL3,min,max,nParL3_);
  int i = 0;
  std::vector<double>::const_iterator parIt = parL3_.begin();
  for(; parIt != parL3_.end(); parIt++, i++) {
    fL3->SetParameter(i,*parIt);
  }
    
  hFrame->GetYaxis()->SetRangeUser(0,4);
  hFrame->GetYaxis()->SetTitle("Absolute correction");
  
  can->cd();
  hFrame->Draw();
  fL3->Draw("same");
  can->SetLogx(1);
  can->Draw();
  ps->NewPage();

  ps->Close();

  //delete hFrame;
  //delete can;
  //delete ps;
}


void run(const std::string& fileName, double min = 4, double max = 2000) {
  readParameters(fileName);
  //  printParameters();
  draw(min,max);
}


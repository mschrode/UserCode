#include <iostream>
#include <vector>

#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"       
#include "Minuit2/MnHesse.h"        
#include "Minuit2/MnUserParameters.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"



// --- Global parameters ---
double ptGenMin_ = 50.;
double ptGenMax_ = 100.;
double etaMin_ = -3.;
double etaMax_ = 3.;
double stochTerm_ = 1.25;



// --- Classes and typedefs ---
class Event {
public:
  Event(double ptGen, double pt1, double pt2, double eta1, double eta2, double error1, double error2)
    : ptGen_(ptGen),
      pt1_(pt1),
      pt2_(pt2),
      eta1_(eta1),
      eta2_(eta2),
      error1_(error1),
      error2_(error2) {};
  ~Event() {};

  double ptGen() const { return ptGen_; }
  double pt1() const { return pt1_; }
  double pt2() const { return pt2_; }
  double eta1() const { return eta1_; }
  double eta2() const { return eta2_; }
  double error1() const { return error1_; }
  double error2() const { return error2_; }

  double pt1Corr(const std::vector<double>& par) const { return ptCorr(par,1); };
  double pt2Corr(const std::vector<double>& par) const { return ptCorr(par,2); };

private:
  double ptGen_;
  double pt1_;
  double pt2_;
  double eta1_;
  double eta2_;
  double error1_;
  double error2_;

  double ptCorr(const std::vector<double>& par, int n) const;
};

double Event::ptCorr(const std::vector<double>& par, int n) const {
  double pt = ( n == 1 ) ? pt1() : pt2();
  double eta = ( n == 1 ) ? eta1() : eta2();
  double c = ( eta < 0 ) ? par.at(0) : par.at(1);
  return c*pt;
}

typedef std::vector<Event> Data;
typedef std::vector<Event>::const_iterator DataIt;


class PtBalanceFitterFCN : public ROOT::Minuit2::FCNBase {
public:
  PtBalanceFitterFCN(const Data& data)
    : data_(data) {};
  ~PtBalanceFitterFCN();

  virtual double operator()(const std::vector<double>& par) const;
  virtual double Up() const { return 1.; }

private:
  Data data_;
};

double PtBalanceFitterFCN::operator()(const std::vector<double>& par) const {
  double chi2 = 0.;
  for(DataIt evt = data_.begin(); evt != data.end(); evt++) {
    double cPt1 = evt->pt1Corr(par);
    double cPt2 = evt->pt2Corr(par);
    double dPt1 = cPt1 / evt->pt1();
    double dPt2 = cPt2 / evt->pt2();

    double res = cPt1 - cPt2;    
    double dRes2 = dPt1 * dPt1 * evt->error1() * evt->error1();
    dRes2 += dPt2 * dPt2 * evt->error2() * evt->error2();

    chi2 += res * res / dRes2;
  }
  return chi2;
}



// --- Global functions ---
Data generateDijets(int nEvents) {
  Data data(nEvents);

  std::vector<double> par(2);
  par.at(0) = 0.7;
  par.at(1) = 1.3;
  TRandom3 randGen(0);
  double rand[3];
  for(int n = 0; n < nEvents; n++) {
    randGen.RndmArray(3,rand);
    double ptGen = ptGenMin_ + (ptGenMax_ - ptGenMin_)*rand[0];
    double eta1 = etaMin_ + (etaMax_ - etaMin_)*rand[1];
    double eta2 = etaMin_ + (etaMax_ - etaMin_)*rand[2];
    double pt1Mean = eta1 < 0 ? par.at(0)*ptGen : par.at(1)*ptGen;
    double pt2Mean = eta2 < 0 ? par.at(0)*ptGen : par.at(1)*ptGen;
    double sigma1 = stochTerm_ * sqrt(pt1Mean);
    double sigma2 = stochTerm_ * sqrt(pt2Mean);
    double pt1 = randGen.Gaus(pt1Mean,sigma1);
    double pt2 = randGen.Gaus(pt2Mean,sigma2);

    Event evt(ptGen,pt1,pt2,eta1,eta2,sigma1,sigma2);
    data.at(n) = evt;
  }

  return data;
}



void plotDijets(const Data& data) {
  double ptMin = 0.;
  double ptMax = 200.;
  TH2D * hPt1vsPt2 = new TH2D("hPt1vsPt2",
			      ";p^{1}_{T};p^{2}_{T}",
			      50,ptMin,ptMax,50,ptMin,ptMax);
  TH1D * hPt1 = new TH1D("hPt1",";p^{1}_{T}",
			 50,ptMin,ptMax);
  TH1D * hPt2 = new TH1D("hPt2",";p^{2}_{T}",
			 50,ptMin,ptMax);

  for(DataIt evt = data.begin(); evt != data.end(); evt++) {
    hPt1vsPt2->Fill(evt->pt1,evt->pt2);
    hPt1->Fill(evt->pt1);
    hPt2->Fill(evt->pt2);
  }

  hPt1->Fit("gaus","0I");
  hPt2->Fit("gaus","0I");
  TF1 * fitPt1 = static_cast<TF1*>(hPt1->GetFunction("gaus"));
  TF1 * fitPt2 = static_cast<TF1*>(hPt2->GetFunction("gaus"));

  TCanvas * c1 = new TCanvas("c1","pt1 vs pt2",600,600);
  c1->cd();
  hPt1vsPt2->Draw("BOX");
  
  TCanvas * c2 = new TCanvas("c2","pt1 and pt2",1200,600);
  c2->Divide(2,1);

  c2->cd(1);
  hPt1->Draw();
  fitPt1->Draw("same");

  c2->cd(2);
  hPt2->Draw();
  fitPt2->Draw("same");
  
  std::cout << "Covariance:  " << hPt1vsPt2->GetCovariance() << std::endl;
  std::cout << "Correlation: " << hPt1vsPt2->GetCorrelationFactor() << std::endl;
}


void run(int nEvents) {
  plotDijets(generateDijets(nEvents));
}

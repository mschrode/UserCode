#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TStyle.h"

#include "TROOT.h"
// #include "TVirtualFitter.h"
// #include "TMinuit.h"
#include "Minuit2/FCNBase.h"
// #include "Minuit2/FunctionMinimum.h"
// #include "Minuit2/MnMigrad.h"
// #include "Minuit2/MnHesse.h"
// #include "Minuit2/MnUserParameters.h"
// #include "Minuit2/MnPrint.h"
// #include "Minuit2/MnContours.h"
// #include "Minuit2/MnPlot.h"
// #include "Minuit2/MnParameterScan.h"



// --- Global parameters ---
double ptGenMin_ = 50.;
double ptGenMax_ = 100.;
double etaMin_ = -3.;
double etaMax_ = 3.;
double stochTerm_ = 1.25;



// --- Classes and typedefs ---
class Event {
public:
  Event()
    : ptGen_(0.),
      pt1_(0.),
      pt2_(0.),
      eta1_(0.),
      eta2_(0.),
      error1_(0.),
      error2_(0.) {};
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


class PtBalanceFCN : public ROOT::Minuit2::FCNBase {
public:
  PtBalanceFCN(const Data& data) : data_(data) {};
  ~PtBalanceFCN() {};

  virtual double operator()(const std::vector<double>& par) const;
  virtual double Up() const {return 1.;}
  
private:
  Data data_;
};

double PtBalanceFCN::operator()(const std::vector<double>& par) const {
  double chi2 = 0.;
  for(DataIt evt = data_.begin(); evt != data_.end(); evt++) {
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
Data generateDijets(int nEvents, const std::vector<double>& par) {
  std::cout << "Generating " << nEvents << " dijet events... " << std::flush;
  Data data(nEvents);

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
  std::cout << "ok\n";

  return data;
}



std::vector<double> fitDijets(const Data& data) {
  using namespace ROOT::Minuit2;

  std::cout << "Fitting parameters...\n";
  PtBalanceFCN fcn(data);
  MnUserParameters upar;
  upar.Add("par0", 1.0, 0.1);
  upar.Add("par1", 1.0, 0.1);

  MnMigrad migrad(fcn, upar);
  FunctionMinimum min = migrad();
  std::cout << min << std::endl;  

  std::vector<double> parCorr(2);
  corrPar.at(0) = 1.;
  corrPar.at(1) = 1.;
  if( min.IsValid() ) {
    MnUserParameters uparMin = min.UserParameters();
    for(size_t i = 0; i < parCorr.size(); i++){
      parCorr.at(i) = uparMin.Value(i);
    }
  }
 
  return parCorr;
}



void plotDijets(const Data& data, const std::vector<double>& par) {
  std::cout << "Plotting histograms... " << std::flush;
  gStyle->SetOptStat(0);

  double ptMin = 0.;
  double ptMax = 200.;
  TH2D * hPt1vsPt2 = new TH2D("hPt1vsPt2",";p^{1}_{T};p^{2}_{T}",
			      50,ptMin,ptMax,50,ptMin,ptMax);
  TH2D * hRespVsEta = new TH2D("hRespVsEta",";#eta;p_{T} / p^{gen}_{T}",
			       25,etaMin_,etaMax_,51,0,2);
  TH2D * hCorrRespVsEta = static_cast<TH2D*>(hRespVsEta->Clone("hCorrRespVsEta"));

  for(DataIt evt = data.begin(); evt != data.end(); evt++) {
    double ptGen = evt->ptGen();
    double pt1 = evt->pt1();
    double pt2 = evt->pt2();
    double pt1Corr = evt->pt1Corr(par);
    double pt2Corr = evt->pt2Corr(par);

    hPt1vsPt2->Fill(pt1,pt2);
    hRespVsEta->Fill(evt->eta1(),pt1/ptGen);
    hRespVsEta->Fill(evt->eta2(),pt2/ptGen);
    hCorrRespVsEta->Fill(evt->eta1(),pt1Corr/ptGen);
    hCorrRespVsEta->Fill(evt->eta2(),pt2Corr/ptGen);
  }

  TProfile * pRespVsEta = hRespVsEta->ProfileX("pRespVsEta",1,-1);
  pRespVsEta->GetYaxis()->SetRangeUser(0,2);
  TProfile * pCorrRespVsEta = hCorrRespVsEta->ProfileX("pCorrRespVsEta",1,-1);
  pCorrRespVsEta->GetYaxis()->SetRangeUser(0,2);

  TCanvas * c1 = new TCanvas("c1","pt1 vs pt2",600,600);
  c1->cd();
  hPt1vsPt2->Draw("BOX");
  
  TCanvas * c2 = new TCanvas("c2","Response",800,800);
  c2->Divide(2,2);

  c2->cd(1);
  hRespVsEta->Draw("BOX");

  c2->cd(2);
  hCorrRespVsEta->Draw("BOX");
  
  c2->cd(3);
  pRespVsEta->Draw();

  c2->cd(4);
  pCorrRespVsEta->Draw();

  std::cout << "ok\n";
}


void run(int nEvents) {
  // Generate data
  std::vector<double> parResp(2);
  parResp.at(0) = 0.7;
  parResp.at(1) = 1.3;
  Data data = generateDijets(nEvents,parResp);

  // Plot data
  std::vector<double> parCorr = fitDijets(data);
  plotDijets(data,parCorr);
}

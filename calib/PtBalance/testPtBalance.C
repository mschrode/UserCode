#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TStyle.h"

#include "TROOT.h"
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"



// === Global parameters ===
double ptGenMin_ = 50.;
double ptGenMax_ = 500.;
double etaMin_ = -3.;
double etaMax_ = 3.;
double stochTerm_ = 1.25;



// === Declaration of classes and typedefs ===
class Event;
class PtBalanceFCN;

typedef std::vector<Event> Data;
typedef std::vector<Event>::const_iterator DataIt;



// === Declaration of global functions ===
Data generateDijets(int nEvents, const std::vector<double>& par);
std::vector<double> fitDijets(const Data& data);
void plotDijets(const Data& data, const std::vector<double>& par);
double getEnergy(double pt, double eta);



// === Main function ===
void run(int nEvents) {
  std::vector<double> parResp(2);
  parResp.at(0) = 0.4;
  parResp.at(1) = 1.1;
  Data data = generateDijets(nEvents,parResp);

  std::vector<double> parCorr = fitDijets(data);
  plotDijets(data,parCorr);
}



// === Implementation of classes ===
// ---------------------------------------------------------------
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
  double c = ( eta > 0 ) ? par.at(0) : par.at(1);
  return c*pt;
}


// ---------------------------------------------------------------
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
    double cPtDijet = 0.5*(cPt1+cPt2);
    double cPtBal = (cPt1-cPt2)/(cPt1+cPt2);
    double dPt1 = cPt1 / evt->pt1();
    double dPt2 = cPt2 / evt->pt2();
    double meanPt = 0.5*(evt->pt1()+evt->pt2());//evt->ptGen();
    double error1 = stochTerm_*sqrt(getEnergy(meanPt,evt->eta1()));
    double error2 = stochTerm_*sqrt(getEnergy(meanPt,evt->eta2()));

    double res = cPt1 - cPt2;    
    double dRes2 = dPt1 * dPt1 * error1 * error1;
    dRes2 += dPt2 * dPt2 * error2 * error2;

//     double res = 2*cPtBal;
//     double dRes2 = (1 - cPtBal)*(1 - cPtBal) * dPt1*dPt1;
//     dRes2 += (1 + cPtBal)*(1 + cPtBal) * dPt2*dPt2;
//     dRes2 *= (error*error);
//     dRes2 /= (cPtDijet*cPtDijet);
    
    chi2 += res * res / dRes2;
  }
  return chi2;
}



// === Implementation of global functions ===
// ---------------------------------------------------------------
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
    double eGen1 = getEnergy(ptGen,eta1);
    double eGen2 = getEnergy(ptGen,eta2);
    double resp1 = eta1 > 0 ? par.at(0) : par.at(1);
    double resp2 = eta2 > 0 ? par.at(0) : par.at(1);
    double sigma1 = stochTerm_ * sqrt(eGen1) * resp1;
    double sigma2 = stochTerm_ * sqrt(eGen2) * resp2;
    double pt1 = randGen.Gaus(ptGen*resp1,sigma1);
    double pt2 = randGen.Gaus(ptGen*resp2,sigma2);

    Event evt(ptGen,pt1,pt2,eta1,eta2,sigma1,sigma2);
    data.at(n) = evt;
  }
  std::cout << "ok\n";

  return data;
}


// ---------------------------------------------------------------
std::vector<double> fitDijets(const Data& data) {
  using namespace ROOT::Minuit2;

  std::cout << "Fitting parameters...\n";
  PtBalanceFCN fcn(data);
  MnUserParameters upar;
  upar.Add("par0", 1.0, 0.01);
  upar.Add("par1", 1.0, 0.01);

  MnMigrad migrad(fcn, upar);
  FunctionMinimum min = migrad();
  std::cout << min << std::endl;  

  std::vector<double> parCorr(2);
  parCorr.at(0) = 1.;
  parCorr.at(1) = 1.;
  if( min.IsValid() ) {
    MnUserParameters uparMin = min.UserParameters();
    for(size_t i = 0; i < parCorr.size(); i++){
      parCorr.at(i) = uparMin.Value(i);
    }
  }
 
  return parCorr;
}


// ---------------------------------------------------------------
void plotDijets(const Data& data, const std::vector<double>& par) {
  std::cout << "Plotting histograms... " << std::flush;
  gStyle->SetOptStat(0);

  double ptMin = 0.;
  double ptMax = 800.;
  TH2D * hPt1vsPt2 = new TH2D("hPt1vsPt2",";p^{1}_{T};p^{2}_{T}",
			      50,ptMin,ptMax,50,ptMin,ptMax);
  TH2D * hRespVsEta = new TH2D("hRespVsEta",";#eta;p_{T} / p^{gen}_{T}",
			       25,etaMin_,etaMax_,51,0,2);
  TH2D * hCorrRespVsEta = static_cast<TH2D*>(hRespVsEta->Clone("hCorrRespVsEta"));

  TH2D * hBalVsEta = new TH2D("hBalVsEta",";#eta;Balance",
			       25,etaMin_,etaMax_,51,-1,1);
  TH2D * hCorrBalVsEta = static_cast<TH2D*>(hBalVsEta->Clone("hCorrBalVsEta"));

  TH2D * hRespVsPtGen = new TH2D("hRespVsPtGen",";p^{gen}_{T};p_{T} / p^{gen}_{T}",
			       25,ptGenMin_,ptGenMax_,51,0,2);
  TH2D * hCorrRespVsPtGen = static_cast<TH2D*>(hRespVsPtGen->Clone("hCorrRespVsPtGen"));


  for(DataIt evt = data.begin(); evt != data.end(); evt++) {
    double ptGen = evt->ptGen();
    double pt1 = evt->pt1();
    double pt2 = evt->pt2();
    double pt1Corr = evt->pt1Corr(par);
    double pt2Corr = evt->pt2Corr(par);

    hPt1vsPt2->Fill(pt1,pt2);
    hRespVsEta->Fill(evt->eta1(),pt1/ptGen);
    hRespVsEta->Fill(evt->eta2(),pt2/ptGen);
    hRespVsPtGen->Fill(ptGen,pt1/ptGen);
    hRespVsPtGen->Fill(ptGen,pt2/ptGen);
    hBalVsEta->Fill(evt->eta1(),2*(pt1-pt2)/(pt1+pt2));

    hCorrRespVsEta->Fill(evt->eta1(),pt1Corr/ptGen);
    hCorrRespVsEta->Fill(evt->eta2(),pt2Corr/ptGen);
    hCorrRespVsPtGen->Fill(ptGen,pt1Corr/ptGen);
    hCorrRespVsPtGen->Fill(ptGen,pt2Corr/ptGen);
    hCorrBalVsEta->Fill(evt->eta1(),2*(pt1Corr-pt2Corr)/(pt1Corr+pt2Corr));
  }

  TProfile * pRespVsEta = hRespVsEta->ProfileX("pRespVsEta");
  pRespVsEta->GetYaxis()->SetRangeUser(0,2);
  TProfile * pCorrRespVsEta = hCorrRespVsEta->ProfileX("pCorrRespVsEta");
  pCorrRespVsEta->SetLineColor(2);
  pCorrRespVsEta->SetMarkerColor(2);

  TProfile * pRespVsPtGen = hRespVsPtGen->ProfileX("pRespVsPtGen");
  pRespVsPtGen->GetYaxis()->SetRangeUser(0,2);
  TProfile * pCorrRespVsPtGen = hCorrRespVsPtGen->ProfileX("pCorrRespVsPtGen");
  pCorrRespVsPtGen->SetLineColor(2);
  pCorrRespVsPtGen->SetMarkerColor(2);

  TProfile * pBalVsEta = hBalVsEta->ProfileX("pBalVsEta");
  pBalVsEta->GetYaxis()->SetRangeUser(-1,1);
  TProfile * pCorrBalVsEta = hCorrBalVsEta->ProfileX("pCorrBalVsEta");
  pCorrBalVsEta->SetLineColor(2);
  pCorrBalVsEta->SetMarkerColor(2);

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
  hRespVsPtGen->Draw("BOX");
  c2->cd(4);
  hCorrRespVsPtGen->Draw("BOX");

  TCanvas * c3 = new TCanvas("c3","Response profile",800,800);
  c3->Divide(2,2);
  c3->cd(1);
  pRespVsEta->Draw();
  pCorrRespVsEta->Draw("same");
  c3->cd(2);
  double meanRespEta = pCorrRespVsEta->GetMean(2);
  pCorrRespVsEta->GetYaxis()->SetRangeUser(meanRespEta-0.05,meanRespEta+0.05);
  pCorrRespVsEta->Draw();
  c3->cd(3);
  pRespVsPtGen->Draw();
  pCorrRespVsPtGen->Draw("same");
  c3->cd(4);
  pCorrRespVsPtGen->Draw();

  TCanvas * c4 = new TCanvas("c4","Balance profile",800,800);
  c4->Divide(2,2);
  c4->cd(1);
  hBalVsEta->Draw("BOX");
  c4->cd(2);
  hCorrBalVsEta->Draw("BOX");
  c4->cd(3);
  pBalVsEta->Draw();
  pCorrBalVsEta->Draw("same");
  c4->cd(4);
  double meanBalEta = pCorrBalVsEta->GetMean(2);
  pCorrBalVsEta->GetYaxis()->SetRangeUser(meanBalEta-0.05,meanBalEta+0.05);
  pCorrBalVsEta->Draw();

  std::cout << "ok\n";
}



// ---------------------------------------------------------------
double getEnergy(double pt, double eta) {
  return pt / std::abs(sin(2*atan(exp(-eta))));
}

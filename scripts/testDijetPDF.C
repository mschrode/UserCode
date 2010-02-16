#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"


// ===== Type and function declarations =====

class Jet;
class DijetEvent;

typedef std::vector<DijetEvent*> Data;
typedef std::vector<DijetEvent*>::const_iterator DataIt;

Data generateData(int nEvts, double tMin, double tMax, double n);
double sigma(double pt);
double measPDF(double pt, double t);
double norm(double tMin, double tMax, double n);
double dijetPDF(double pt1, double pt2, double tMin, double tMax, double n);
double nLogL(const Data &data, double tMin, double tMax, double n);
double prob(double pt1Min, double pt1Max, double pt2Min, double pt2Max, double tMin, double tMax, double n);
double truthPDF(double *x, double *par);
double meanTruth(double tMin, double tMax, double n);
void plotData(const Data &data, const std::vector<double> &tBinEdges, double n);
void scanLkh(const Data &data, int nSteps, double eps, double tMin, double tMax, double trueN);
void compareNEvents(const Data &data, double tMin, double tMax, double n);
void testDijetPDF(int nEvts = 10000);



// ===== Class implementations =====

class Jet {
public:
  Jet(double ptGen, double pt)
    : pt_(pt), ptGen_(ptGen), resp_(pt/ptGen) {};

  Jet(const Jet *j) 
    : pt_(j->pt()), ptGen_(j->ptGen()), resp_(j->resp()) {};

  double pt() const { return pt_; }
  double ptGen() const { return ptGen_; }
  double resp() const { return resp_; }

private:
  const double pt_;
  const double ptGen_;
  const double resp_;
};


class DijetEvent {
public:
  DijetEvent(const Jet *j1, const Jet *j2) 
    : j1_(new Jet(j1)), j2_(new Jet(j2)) {};

  ~DijetEvent() {
    delete j1_;
    delete j2_;
  }

  double pt(int i) const { return i == 0 ? j1_->pt() : j2_->pt(); }
  double ptGen(int i) const { return i == 0 ? j1_->ptGen() : j2_->ptGen(); }
  double resp(int i) const { return i == 0 ? j1_->resp() : j2_->resp(); }
  
private:
  Jet *j1_;
  Jet *j2_;
};





// ===== Function implementations =====

Data generateData(int nEvts, double tMin, double tMax, double n) {
  std::cout << "Generating " << nEvts << " dijet events\n";
  
  Data data;
  TRandom3 *rand = new TRandom3(0);
  for(int i = 0; i < nEvts; i++) {
    double ptGen = tMax + 1.;
    while( ptGen > tMax ) {
      double x = rand->Uniform();
      ptGen = tMin*pow(x,-1./(n-1.));
    }
    double s = sigma(ptGen);
    Jet *jet[2];
    for(int j = 0; j < 2; j++) {
      double pt = rand->Gaus(ptGen,s);
      jet[j] = new Jet(ptGen,pt);
    }
    data.push_back(new DijetEvent(jet[0],jet[1]));
    for(int j = 0; j < 2; j++) {
      delete jet[j];
    }	
  }
  delete rand;
  
  return data;
}


double sigma(double pt) {
  double a = 4.44;
  double b = 1.11;
  double c = 0.03;
  return sqrt( a*a + b*b*pt + c*c*pt*pt );
  //  return 10.;
}


double measPDF(double pt, double t) {
  double s = sigma(t);
  double u = (pt-t)/s;
  return exp(-u*u/2.)/sqrt(2.*M_PI)/s;
}


double norm(double tMin, double tMax, double n) {
  double m = 1.-n;
  return m == 0. ? log(tMax/tMin) : ( pow(tMax,m) - pow(tMin,m) )/m;
}


double dijetPDF(double pt1, double pt2, double tMin, double tMax, double n) {
  double p = 0.;
  int nSteps = 500;
  double dt = (tMax-tMin)/nSteps;
  for(int i = 0; i < nSteps; i++) {
    double t = tMin + (i+0.5)*dt;
    p += measPDF(pt1,t)*measPDF(pt2,t)/pow(t,n);
  }
  p *= dt;

  double c = norm(tMin,tMax,n);
  if( c ) p /= c;

  assert( p == p );
  if( p <= 0 ) {
    std::cerr << "ERROR: Probability zero or negative (" << p << ")\n";
    exit(-1);
  }

  return p;
}


double nLogL(const Data &data, double tMin, double tMax, double n) {
  double lkh = 0.;
  for(DataIt it = data.begin(); it != data.end(); it++) {
    DijetEvent *evt = *it;
    lkh -= log(dijetPDF(evt->pt(0),evt->pt(1),tMin,tMax,n));
  }
  return lkh;
}


double prob(double pt1Min, double pt1Max, double pt2Min, double pt2Max, double tMin, double tMax, double n) {
  double p = 0.;
  int nSteps = 50;
  double dpt1 = (pt1Max-pt1Min)/nSteps;
  double dpt2 = (pt2Max-pt2Min)/nSteps;
  for(int i1 = 0; i1 < nSteps; i1++) {
    double pt1 = pt1Min + (i1+0.5)*dpt1;
    for(int i2 = 0; i2 < nSteps; i2++) {
      double pt2 = pt2Min + (i2+0.5)*dpt2;
      p += dijetPDF(pt1,pt2,tMin,tMax,n);
    }
  }
  p *= dpt1*dpt2;

  return p;
}


double respPDF(double *x, double *par) {
  double s = sigma(par[0])/par[0];
  double u = (x[0]-1.)/s;
  return exp(-u*u/2.)/sqrt(2.*M_PI)/s;

  double m = 1.-par[0];
  return m/(pow(par[2],m)-pow(par[1],m))/pow(x[0],par[0]);
}


double truthPDF(double *x, double *par) {
  double m = 1.-par[0];
  return m/(pow(par[2],m)-pow(par[1],m))/pow(x[0],par[0]);
}


double meanTruth(double tMin, double tMax, double n) {
  double m = 2.-n;
  return m == 0. ? log(tMax/tMin)/norm(tMin,tMax,n) : ( pow(tMax,m) - pow(tMin,m) )/m/norm(tMin,tMax,n);
}


void plotData(const Data &data, const std::vector<double> &tBinEdges, double n) {
  std::cout << "Plotting data\n";
  
  TH1D *hPtGen = new TH1D("hPtGen",";p^{gen}_{T} (GeV)",
			  100,tBinEdges.front(),tBinEdges.back());
  int nBins = static_cast<int>(tBinEdges.size()-1);
  std::vector<TH1D*> hResp(nBins);
  for(int i = 0; i < nBins; i++) {
    TString name = "hResp";
    name += i;
    char title[50];
    sprintf(title,"%.0f < p^{gen}_{T} < %.0f GeV;Response",tBinEdges[i],tBinEdges[i+1]);
    hResp[i] = new TH1D(name,title,100,0,2);
  }
  for(DataIt it = data.begin(); it != data.end(); it++) {
    DijetEvent *evt = *it;

    for(int tBin = 0; tBin < nBins; tBin++) {
      if( tBinEdges[tBin] <= evt->ptGen(0) && evt->ptGen(0) < tBinEdges[tBin+1] ) {
	for(int j = 0; j < 2; j++) {
	  hPtGen->Fill(evt->ptGen(j));
	  hResp[tBin]->Fill(evt->pt(j)/evt->ptGen(j));
	}
	break;
      }
    }
  }

  if( hPtGen->Integral("width") ) hPtGen->Scale(1./hPtGen->Integral("width"));
  TF1 *fPtGen = new TF1("fPtGen",truthPDF,tBinEdges.front(),tBinEdges.back(),3);
  fPtGen->SetLineWidth(1);
  fPtGen->SetLineColor(2);
  fPtGen->SetParameter(0,n);
  fPtGen->SetParameter(1,tBinEdges.front());
  fPtGen->SetParameter(2,tBinEdges.back());

  TCanvas *can = new TCanvas("can","Spectrum",500,500);
  can->cd();
  hPtGen->Draw();
  fPtGen->Draw("same");
  gPad->SetLogy();

  for(int i = 0; i < nBins; i++) {
    if( hResp[i]->Integral("width") ) hResp[i]->Scale(1./hResp[i]->Integral("width"));
    TString name = "fResp";
    name += i;
    TF1 *fResp = new TF1(name,respPDF,0,2,1);
    fResp->SetLineWidth(1);
    fResp->SetLineColor(2);
    fResp->SetParameter(0,meanTruth(tBinEdges[i],tBinEdges[i+1],n));
    
    std::cout << " " << tBinEdges[i] << " - " << tBinEdges[i+1] << ": " << std::flush;
    std::cout << meanTruth(tBinEdges[i],tBinEdges[i+1],n) << std::endl;

    name = "ResBin_";
    name += i;
    TCanvas *canR = new TCanvas(name,name,500,500);
    canR->cd();
    hResp[i]->Draw();
    fResp->Draw("same");
  }
}


void scanLkh(const Data &data, int nSteps, double eps, double tMin, double tMax, double trueN) {
  std::cout << "Scanning n" << std::endl;
  int nPoints = 2*nSteps + 1;
  std::vector<double> nVal(nPoints);
  std::vector<double> lkh(nPoints);
  for(int i = 0; i < nPoints; i++) {
    nVal[i] = eps*(i-nSteps)+trueN;
    lkh[i] = nLogL(data,tMin,tMax,nVal[i]);

    std::cout << " " << nVal[i] << ": " << lkh[i] << std::endl;
  }
  TGraph *gLkh = new TGraph(nPoints,&(nVal.front()),&(lkh.front()));
  gLkh->SetMarkerStyle(20);
  gLkh->SetTitle("Parameter scan;n;-ln(L)");

  TCanvas *can = new TCanvas("canScan","Scan",500,500);
  can->cd();
  gLkh->Draw("AP");
}


void compareNEvents(const Data &data, double tMin, double tMax, double n) {
  int nBins = 10;
  double min = 0.5*tMin;
  double max = 1.5*tMax;
  TH2D *h2Data = new TH2D("h2Data",";p^{1}_{T} (GeV);p^{2}_{T} (GeV)",
			  nBins,min,max,nBins,min,max);
  for(DataIt it = data.begin(); it != data.end(); it++) {
    if( (*it)->pt(0) < min || (*it)->pt(1) < min ) continue;
    if( (*it)->pt(0) > max || (*it)->pt(1) > max ) continue;
    h2Data->Fill((*it)->pt(0),(*it)->pt(1));
  }

  TH2D *h2Pred = static_cast<TH2D*>(h2Data->Clone("h2Pred"));
  std::vector<TH1D*> hData(nBins);
  std::vector<TH1D*> hPred(nBins);
  std::vector<TH1D*> hPredErr(nBins);
  for(int xBin = 1; xBin <= nBins; xBin++) {
    double xMin = h2Data->GetXaxis()->GetBinLowEdge(xBin);
    double xMax = h2Data->GetXaxis()->GetBinUpEdge(xBin);
    std::cout << "\n" << xMin << " < pt1 < " << xMax << std::endl;

    TString name = "hData";
    name += xBin-1;
    char title[50];
    sprintf(title,"%.1f < p^{1}_{T} < %.1f GeV;p^{2}_{T} (GeV)",xMin,xMax);
    hData[xBin-1] = new TH1D(name,title,nBins,
			     h2Data->GetYaxis()->GetBinLowEdge(1),
			     h2Data->GetYaxis()->GetBinLowEdge(nBins));
    hData[xBin-1]->SetMarkerStyle(20);

    name = "hPred";
    name += xBin-1;
    hPred[xBin-1] = static_cast<TH1D*>(hData[xBin-1]->Clone(name));
    hPred[xBin-1]->SetMarkerStyle(1);
    hPred[xBin-1]->SetLineWidth(2);

    name = "hPredErr";
    name += xBin-1;
    hPredErr[xBin-1] = static_cast<TH1D*>(hPred[xBin-1]->Clone(name));
    hPredErr[xBin-1]->SetFillColor(45);


    for(int yBin = 1; yBin <= nBins; yBin++) {
      double yMin = h2Data->GetYaxis()->GetBinLowEdge(yBin);
      double yMax = h2Data->GetYaxis()->GetBinUpEdge(yBin);
      std::cout << " " << yMin << " < pt2 < " << yMax << std::flush;
      
      double val = prob(xMin,xMax,yMin,yMax,tMin,tMax,n);
      std::cout << ": " << val << std::flush;
      int bin = h2Data->GetBin(xBin,yBin);
      h2Pred->SetBinContent(bin,val);
      std::cout << " (" << val*data.size() << ")" << std::flush;
      std::cout << " -- " << h2Data->GetBinContent(bin) << std::endl;

      hData[xBin-1]->SetBinContent(yBin,h2Data->GetBinContent(bin));
      hPred[xBin-1]->SetBinContent(yBin,val*data.size());
      hPredErr[xBin-1]->SetBinContent(yBin,val*data.size());
      hPredErr[xBin-1]->SetBinError(yBin,sqrt(val*data.size()));
    }
  }    


  TCanvas *cProb = new TCanvas("cProb","# events",1000,500);
  cProb->Divide(2,1);
  cProb->cd(1);
  h2Data->Draw("COLZ");
  cProb->cd(2);
  h2Pred->Draw("COLZ");

  for(int i = 0; i < nBins; i ++) {
    TString name = "NEvts_";
    name += i;
    TCanvas *c = new TCanvas(name,name,500,500);
    c->cd();
    hPredErr[i]->Draw("E3");
    hData[i]->Draw("Psame");
    hPred[i]->Draw("same");
  }
}


void testDijetPDF(int nEvts) {
  std::vector<double> tBinEdges;
  tBinEdges.push_back(100.);
  tBinEdges.push_back(120.);
  tBinEdges.push_back(150.);
  tBinEdges.push_back(200.);
  tBinEdges.push_back(500.);

  double tMin = tBinEdges.front();
  double tMax = tBinEdges.back();
  double n = 6.;
  
  Data data = generateData(nEvts,tMin,tMax,n);
  plotData(data,tBinEdges,n);
  //scanLkh(data,5,0.1,tMin,tMax,n);
  compareNEvents(data,tMin,tMax,n);
}

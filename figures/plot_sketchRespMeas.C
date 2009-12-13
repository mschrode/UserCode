#include <iostream>
#include <vector>

#include "TArrow.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"


// === Global parameters ===
double ptGenMin_ = 100.;
double ptGenMax_ = 500.;
double stochTerm_ = 2.5;
int nBinsX_ = 20;
double xMin_ = 0.;
double xMax_ = 600.;
double fitRange_ = 1.25;


// === Declaration of classes and typedefs ===
class Event;

typedef std::vector<Event> Data;
typedef std::vector<Event>::const_iterator DataIt;


// === Declaration of global functions ===
Data generateData(int nEvents, const std::vector<double>& par);
void plotRespMeasVsTruth(const Data& data);
void plotRespMeasVsMeas(const Data& data);
TString truthTitle();
TString measTitle();
TString responseTitle();
void setGStyle();


// === Main function ===
void run(int nEvents) {
  setGStyle();
  std::vector<double> parResp(2);
  parResp.at(0) = 0.;
  parResp.at(1) = 1.;
  Data data = generateData(nEvents,parResp);
  plotRespMeasVsTruth(data);
  plotRespMeasVsMeas(data);
}



// === Implementation of classes ===
// ---------------------------------------------------------------
class Event {
public:
  Event()
    : ptGen_(0.),
      pt_(0.),
      error_(0.) {};
  Event(double ptGen, double pt, double error)
    : ptGen_(ptGen),
      pt_(pt),
      error_(error) {};

  double ptGen() const { return ptGen_; }
  double pt() const { return pt_; }
  double error() const { return error_; }

private:
  double ptGen_;
  double pt_;
  double error_;
};


// === Implementation of global functions ===
// ---------------------------------------------------------------
Data generateData(int nEvents, const std::vector<double>& par) {
  std::cout << "Generating " << nEvents << " events... " << std::flush;
  Data data(nEvents);

  TRandom3 randGen(0);
  double rand[1];
  for(int n = 0; n < nEvents; n++) {
    randGen.RndmArray(1,rand);
    double ptGen = ptGenMin_ + (ptGenMax_ - ptGenMin_)*rand[0];
    double resp = 1. - par[0] / (ptGen + par[1]);
    double sigma = stochTerm_ * sqrt(ptGen) * resp;
    double pt = randGen.Gaus(ptGen*resp,sigma);
    while( pt < 0. || pt > 2*ptGen*resp ) {
      pt = randGen.Gaus(ptGen*resp,sigma);
    }

    Event evt(ptGen,pt,sigma);
    data.at(n) = evt;
  }
  std::cout << "ok\n";

  return data;
}



void plotRespMeasVsTruth(const Data& data) {
  TH2D *hRespVsTruth = new TH2D("hRespVsTruth","",
				nBinsX_,xMin_,xMax_,
				51,0.,2.);
  hRespVsTruth->SetXTitle(truthTitle());
  hRespVsTruth->SetYTitle(responseTitle());
  TH2D *hMeasVsTruth = new TH2D("hMeasVsTruth","",
				nBinsX_,xMin_,xMax_,
				nBinsX_,xMin_,xMax_);
  hMeasVsTruth->SetXTitle(truthTitle());
  hMeasVsTruth->SetYTitle(measTitle());

  // Filling histos
  for(DataIt evt = data.begin(); evt != data.end(); evt++) {
    hRespVsTruth->Fill(evt->ptGen(),evt->pt()/evt->ptGen());
    hMeasVsTruth->Fill(evt->ptGen(),evt->pt());
  }

  // Fitting profiles
  // Create profile histogram
  TString name = hRespVsTruth->GetName();
  name += "Profile";
  TH1D *hProf = new TH1D(name,"",nBinsX_,xMin_,xMax_);
  hProf->SetMarkerStyle(20);
  hProf->SetXTitle(hRespVsTruth->GetXaxis()->GetTitle());
  name = "GaussFit <";
  name += hRespVsTruth->GetYaxis()->GetTitle();
  name += ">";
  hProf->SetYTitle(name);
  hProf->GetYaxis()->SetRangeUser(0.95,1.05);

  // Fit distributions
  std::vector<TH1D*> hDists;
  std::vector<TF1*> fits;
  for(int xBin = 1; xBin <= hRespVsTruth->GetNbinsX(); xBin++) {
    name = hRespVsTruth->GetName();
    name += "Dist";
    name += xBin;
    TH1D *htemp = new TH1D(name,"",
			   hRespVsTruth->GetNbinsY(),
			   hRespVsTruth->GetYaxis()->GetXmin(),
			   hRespVsTruth->GetYaxis()->GetXmax());
    htemp->SetXTitle(hRespVsTruth->GetYaxis()->GetTitle());
    htemp->Sumw2();

    // Copy over bin contents (y-slice)
    for(int yBin = 1; yBin <= hRespVsTruth->GetNbinsY(); yBin++) {
      htemp->SetBinContent(yBin,hRespVsTruth->GetBinContent(hRespVsTruth->GetBin(xBin,yBin)));
      htemp->SetBinError(yBin,hRespVsTruth->GetBinError(xBin,yBin));
    }  
    hDists.push_back(htemp);

    name = hRespVsTruth->GetName();
    name += "Fit";
    name += xBin;
    double mean = htemp->GetMean();
    double sigma = htemp->GetRMS();
    TF1 *fitTemp = new TF1(name,"gaus",mean - fitRange_*sigma,mean + fitRange_*sigma);
    htemp->Fit(fitTemp,"0QLR");
    fits.push_back(fitTemp);
    hProf->SetBinContent(xBin,fitTemp->GetParameter(1));
    hProf->SetBinError(xBin,fitTemp->GetParError(1));
  }

  // Define bin whose distribution is drawn
  int bin = 7;
  double binMin = hProf->GetXaxis()->GetBinLowEdge(bin);
  double binMax = hProf->GetXaxis()->GetBinUpEdge(bin);
  hDists[bin]->GetXaxis()->SetRangeUser(0.6,1.4);
  hDists[bin]->GetYaxis()->SetRangeUser(0,1.2*hDists[bin]->GetMaximum());
  fits[bin]->SetLineColor(4);
  fits[bin]->SetLineWidth(2);
  
  TLine *diag = new TLine(xMin_,xMin_,xMax_,xMax_);
  diag->SetLineColor(4);
  diag->SetLineStyle(2);
  diag->SetLineWidth(2);

  TBox *box = new TBox(binMin,xMin_,binMax,xMax_);
  box->SetFillStyle(0);
  box->SetLineColor(2);
  box->SetLineWidth(3);

  double aLeft   = 0.4*(binMax/(xMax_-xMin_)+0.1);
  double aRight  = 0.38;
  double aHeight = 0.4;
  TArrow *hArrow1 = new TArrow(aLeft,aHeight,aRight,aHeight+0.1,0.03,"|>");
  hArrow1->SetAngle(30);
  hArrow1->SetLineColor(2);
  hArrow1->SetFillColor(2);
  hArrow1->SetLineWidth(3);

  aLeft   = 0.56;
  aRight  = 0.6 + 0.4*(0.1+(hProf->GetBinCenter(bin))/(xMax_ - xMin_));
  aHeight = 0.4;
  TArrow *hArrow2 = new TArrow(aLeft,aHeight,aRight,aHeight+0.12,0.03,"|>");
  hArrow2->SetAngle(30);
  hArrow2->SetLineColor(2);
  hArrow2->SetFillColor(2);
  hArrow2->SetLineWidth(3);

  TLine *line = new TLine(xMin_,1,xMax_,1);
  line->SetLineWidth(2);
  line->SetLineColor(4);

  TEllipse *ellipse = new TEllipse(hProf->GetBinCenter(bin),
				   hProf->GetBinContent(bin),
				   20,0.007);
  ellipse->SetLineColor(2);
  ellipse->SetLineWidth(3);
  ellipse->SetFillStyle(0);

  // Draw stuff
  TCanvas *can = new TCanvas("canRespVsTruth","RespVsTruth",1000,400);
  can->cd();
  TPad *lPad = new TPad("lPadRespVsTruth","",0.,0.,0.4,1.);
  lPad->SetFillStyle(1001);
  lPad->SetFrameFillColor(10);
  lPad->SetFrameBorderMode(0);
  lPad->Draw();
  lPad->cd();
  hMeasVsTruth->Draw("BOX");
  diag->Draw("same");
  box->Draw("same");

  can->cd();
  TPad *rPad = new TPad("rPadRespVsTruth","",0.6,0.,1.,1.);
  rPad->SetFillStyle(1001);
  rPad->SetFrameFillColor(10);
  rPad->SetFrameBorderMode(0);
  rPad->Draw();
  rPad->cd();
  hProf->Draw("PE1");
  line->Draw("same");
  hProf->Draw("PE1same");
  ellipse->Draw("same");

  can->cd();
  TPad *cPad = new TPad("cPadRespVsTruth","",0.32,0.2,0.59,0.95);
  cPad->SetFillStyle(1001);
  cPad->SetFrameFillColor(10);
  cPad->SetFrameBorderMode(0);
  cPad->SetLeftMargin(0.1);

  cPad->Draw();
  cPad->cd();
  hDists[bin]->Draw("HIST");
  fits[bin]->Draw("same");

  can->cd();
  TPad *pad = new TPad("padRespVsTruth","",0.,0.,1.,1.);
  pad->SetFillStyle(0);
  pad->SetFrameFillColor(10);
  pad->SetFrameBorderMode(0);
  pad->Draw();
  pad->cd();
  hArrow1->Draw();
  hArrow2->Draw();

}


void plotRespMeasVsMeas(const Data& data) {
  TH2D *hRespVsMeas = new TH2D("hRespVsMeas","",
			       nBinsX_,xMin_,xMax_,
			       51,0.,2.);
  hRespVsMeas->SetXTitle(measTitle());
  hRespVsMeas->SetYTitle(responseTitle());
  TH2D *hMeasVsTruth = new TH2D("hMeasVsTruth2","",
				nBinsX_,xMin_,xMax_,
				nBinsX_,xMin_,xMax_);
  hMeasVsTruth->SetXTitle(truthTitle());
  hMeasVsTruth->SetYTitle(measTitle());

  // Filling histos
  for(DataIt evt = data.begin(); evt != data.end(); evt++) {
    hRespVsMeas->Fill(evt->pt(),evt->pt()/evt->ptGen());
    hMeasVsTruth->Fill(evt->ptGen(),evt->pt());
  }

  // Fitting profiles
  // Create profile histogram
  TString name = hRespVsMeas->GetName();
  name += "Profile";
  TH1D *hProf = new TH1D(name,"",nBinsX_,xMin_,xMax_);
  hProf->SetMarkerStyle(20);
  hProf->SetXTitle(hRespVsMeas->GetXaxis()->GetTitle());
  name = "GaussFit <";
  name += hRespVsMeas->GetYaxis()->GetTitle();
  name += ">";
  hProf->SetYTitle(name);
  double respYMin = 0.5;
  double respYMax = 1.2;
  hProf->GetYaxis()->SetRangeUser(respYMin,respYMax);

  // Fit distributions
  std::vector<TH1D*> hDists;
  std::vector<TF1*> fits;
  for(int xBin = 1; xBin <= hRespVsMeas->GetNbinsX(); xBin++) {
    name = hRespVsMeas->GetName();
    name += "Dist";
    name += xBin;
    TH1D *htemp = new TH1D(name,"",
			   hRespVsMeas->GetNbinsY(),
			   hRespVsMeas->GetYaxis()->GetXmin(),
			   hRespVsMeas->GetYaxis()->GetXmax());
    htemp->SetXTitle(hRespVsMeas->GetYaxis()->GetTitle());
    htemp->Sumw2();

    // Copy over bin contents (y-slice)
    for(int yBin = 1; yBin <= hRespVsMeas->GetNbinsY(); yBin++) {
      htemp->SetBinContent(yBin,hRespVsMeas->GetBinContent(hRespVsMeas->GetBin(xBin,yBin)));
      htemp->SetBinError(yBin,hRespVsMeas->GetBinError(xBin,yBin));
    }  
    hDists.push_back(htemp);

    name = hRespVsMeas->GetName();
    name += "Fit";
    name += xBin;
    double mean = htemp->GetMean();
    double sigma = htemp->GetRMS();
    TF1 *fitTemp = new TF1(name,"gaus",mean - fitRange_*sigma,mean + fitRange_*sigma);
    htemp->Fit(fitTemp,"0QLR");
    fits.push_back(fitTemp);
    hProf->SetBinContent(xBin,fitTemp->GetParameter(1));
    hProf->SetBinError(xBin,fitTemp->GetParError(1));
  }

  // Define bin whose distribution is drawn
  int bin1 = 3;
  double bin1Min = hMeasVsTruth->GetYaxis()->GetBinLowEdge(bin1);
  double bin1Max = hMeasVsTruth->GetYaxis()->GetBinUpEdge(bin1);
  hDists[bin1-1]->GetXaxis()->SetRangeUser(0.5,1.5);
  hDists[bin1-1]->GetYaxis()->SetRangeUser(0,1.2*hDists[bin1-1]->GetMaximum());
  fits[bin1-1]->SetLineColor(4);
  fits[bin1-1]->SetLineWidth(2);

  int bin2 = 7;
  double bin2Min = hMeasVsTruth->GetYaxis()->GetBinLowEdge(bin2);
  double bin2Max = hMeasVsTruth->GetYaxis()->GetBinUpEdge(bin2);
  hDists[bin2-1]->GetXaxis()->SetRangeUser(0.5,1.5);
  hDists[bin2-1]->GetYaxis()->SetRangeUser(0,1.2*hDists[bin2-1]->GetMaximum());
  fits[bin2-1]->SetLineColor(4);
  fits[bin2-1]->SetLineWidth(2);
  
  TLine *diag = new TLine(xMin_,xMin_,xMax_,xMax_);
  diag->SetLineColor(4);
  diag->SetLineStyle(2);
  diag->SetLineWidth(2);

  TBox *box1 = new TBox(xMin_,bin1Min,xMax_,bin1Max);
  box1->SetFillStyle(0);
  box1->SetLineColor(2);
  box1->SetLineWidth(3);

  TBox *box2 = new TBox(xMin_,bin2Min,xMax_,bin2Max);
  box2->SetFillStyle(0);
  box2->SetLineColor(2);
  box2->SetLineWidth(3);

  double aLeft   = 0.3;
  double aRight  = 0.43;
  double aHeight = 0.1 + 0.8*(gStyle->GetPadBottomMargin() + (1-gStyle->GetPadTopMargin()-gStyle->GetPadBottomMargin())*(bin1Min/(xMax_-xMin_)));
  TArrow *hArrow11 = new TArrow(aLeft,aHeight,aRight,0.2,0.03,"|>");
  hArrow11->SetAngle(30);
  hArrow11->SetLineColor(2);
  hArrow11->SetFillColor(2);
  hArrow11->SetLineWidth(3);

  aHeight = 0.1 + 0.8*(gStyle->GetPadBottomMargin() + (1-gStyle->GetPadTopMargin()-gStyle->GetPadBottomMargin())*(bin2Max/(xMax_-xMin_)));
  TArrow *hArrow12 = new TArrow(aLeft,aHeight,aRight,0.8,0.03,"|>");
  hArrow12->SetAngle(30);
  hArrow12->SetLineColor(2);
  hArrow12->SetFillColor(2);
  hArrow12->SetLineWidth(3);

  aLeft   = 0.56;
  aRight  = 0.6 + 0.4*(gStyle->GetPadLeftMargin()+(1-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin())*((hProf->GetBinCenter(bin1))/(xMax_ - xMin_)));
  aHeight = 0.1 + 0.8*(gStyle->GetPadBottomMargin() + (1-gStyle->GetPadTopMargin()-gStyle->GetPadBottomMargin())*(hProf->GetBinContent(bin1)-respYMin)/(respYMax-respYMin));
  std::cout << hProf->GetBinContent(bin1) << std::endl;
  TArrow *hArrow21 = new TArrow(aLeft,0.2,aRight,aHeight-0.04,0.03,"|>");
  hArrow21->SetAngle(30);
  hArrow21->SetLineColor(2);
  hArrow21->SetFillColor(2);
  hArrow21->SetLineWidth(3);

  aLeft   = 0.56;
  aRight  = 0.6 + 0.4*(gStyle->GetPadLeftMargin()+(1-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin())*((hProf->GetBinCenter(bin2))/(xMax_ - xMin_)));
  aHeight = 0.1 + 0.8*(gStyle->GetPadBottomMargin() + (1-gStyle->GetPadTopMargin()-gStyle->GetPadBottomMargin())*(hProf->GetBinContent(bin2)-respYMin)/(respYMax-respYMin));
  std::cout << hProf->GetBinContent(bin1) << std::endl;
  TArrow *hArrow22 = new TArrow(aLeft,0.8,aRight,aHeight+0.04,0.03,"|>");
  hArrow22->SetAngle(30);
  hArrow22->SetLineColor(2);
  hArrow22->SetFillColor(2);
  hArrow22->SetLineWidth(3);

  TLine *line = new TLine(xMin_,1,xMax_,1);
  line->SetLineWidth(2);
  line->SetLineColor(4);
  line->SetLineStyle(2);

  TLine *line1 = new TLine(1.,hDists[bin2-1]->GetMinimum(),1.,hDists[bin2-1]->GetMaximum());
  line1->SetLineWidth(2);
  line1->SetLineColor(4);
  line1->SetLineStyle(2);

  TLine *line2 = new TLine(1.,hDists[bin1-1]->GetMinimum(),1.,hDists[bin1-1]->GetMaximum());
  line2->SetLineWidth(2);
  line2->SetLineColor(4);
  line2->SetLineStyle(2);

  TEllipse *ellipse1 = new TEllipse(hProf->GetBinCenter(bin1),
				    hProf->GetBinContent(bin1),
				    24,0.04);
  ellipse1->SetLineColor(2);
  ellipse1->SetLineWidth(3);
  ellipse1->SetFillStyle(0);

  TEllipse *ellipse2 = new TEllipse(hProf->GetBinCenter(bin2),
				    hProf->GetBinContent(bin2),
				    24,0.04);
  ellipse2->SetLineColor(2);
  ellipse2->SetLineWidth(3);
  ellipse2->SetFillStyle(0);

  // Draw stuff
  TCanvas *can = new TCanvas("canRespVsMeas","RespVsMeas",1000,500);
  can->cd();
  TPad *lPad = new TPad("lPadRespVsMeas","",0.,0.1,0.4,0.9);
  lPad->SetFillStyle(1001);
  lPad->SetFrameFillColor(10);
  lPad->SetFrameBorderMode(0);
  lPad->Draw();
  lPad->cd();
  hMeasVsTruth->Draw("BOX");
  diag->Draw("same");
  box1->Draw("same");
  box2->Draw("same");

  can->cd();
  TPad *rPad = new TPad("rPadRespVsMeas","",0.6,0.1,1.,0.9);
  rPad->SetFillStyle(1001);
  rPad->SetFrameFillColor(10);
  rPad->SetFrameBorderMode(0);
  rPad->Draw();
  rPad->cd();
  hProf->Draw("PE1");
  line->Draw("same");
  hProf->Draw("PE1same");
  ellipse1->Draw("same");
  ellipse2->Draw("same");

  can->cd();
  TPad *cPad1 = new TPad("cPad1RespVsMeas","",0.37,0.5,0.6,1.);
  cPad1->SetFillStyle(1001);
  cPad1->SetFrameFillColor(10);
  cPad1->SetFrameBorderMode(0);
  cPad1->SetLeftMargin(0.1);

  cPad1->Draw();
  cPad1->cd();
  hDists[bin2-1]->Draw("HIST");
  fits[bin2-1]->Draw("same");
  line1->Draw("same");

  can->cd();
  TPad *cPad2 = new TPad("cPad2RespVsMeas","",0.37,0.,0.6,0.5);
  cPad2->SetFillStyle(1001);
  cPad2->SetFrameFillColor(10);
  cPad2->SetFrameBorderMode(0);
  cPad2->SetLeftMargin(0.1);

  cPad2->Draw();
  cPad2->cd();
  hDists[bin1-1]->Draw("HIST");
  fits[bin1-1]->Draw("same");
  line2->Draw("same");

  can->cd();
  TPad *pad = new TPad("padRespVsMeas","",0.,0.,1.,1.);
  pad->SetFillStyle(0);
  pad->SetFrameFillColor(10);
  pad->SetFrameBorderMode(0);
  pad->Draw();
  pad->cd();
  hArrow11->Draw();
  hArrow12->Draw();
  hArrow21->Draw();
  hArrow22->Draw();
}


TString truthTitle() {
  return "p^{true}_{T} (GeV)";
}

TString measTitle() {
  return "p^{meas}_{T} (GeV)";
}

TString responseTitle() {
  return "p^{meas}_{T} / p^{true}_{T}";
}

void setGStyle() {
  //  For the canvas
  // -------------------------------------------
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(500); //Height of canvas
  gStyle->SetCanvasDefW(500); //Width of canvas
  gStyle->SetCanvasDefX(0);   //Position on screen
  gStyle->SetCanvasDefY(0);


  //  For the frame
  // -------------------------------------------
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(kBlack);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetFrameLineStyle(0);
  gStyle->SetFrameLineWidth(1);


  //  For the Pad
  // -------------------------------------------
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);


  //  Margins
  // -------------------------------------------
  gStyle->SetPadTopMargin(0.04);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.04);


  //  For the histo:
  // -------------------------------------------
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);


  //  For the statistics box:
  // -------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatX(0.92);              
  gStyle->SetStatY(0.86);              
  gStyle->SetStatH(0.16);
  gStyle->SetStatW(0.22);


  //  For the Global title:
  // -------------------------------------------
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFont(42,"");
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleX(0.58);
  gStyle->SetTitleH(0.05);
  gStyle->SetTitleBorderSize(0);


  //  For the axis
  // -------------------------------------------
  gStyle->SetAxisColor(1,"XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03,"XYZ");
  gStyle->SetNdivisions(510,"XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);


  //  For the axis labels and titles
  // -------------------------------------------
  gStyle->SetTitleColor(1,"XYZ");
  gStyle->SetLabelColor(1,"XYZ");
  // For the axis labels:
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetLabelOffset(0.007,"XYZ");
  gStyle->SetLabelSize(0.045,"XYZ");
  
  // For the axis titles:
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetTitleSize(0.08,"XYZ");
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.0);


  //  For the legend
  // -------------------------------------------
  gStyle->SetLegendBorderSize(1);
  gStyle->SetStatFontSize(0.04);
  gStyle->SetStatX(0.92);              
  gStyle->SetStatY(0.92);              
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.3);
}

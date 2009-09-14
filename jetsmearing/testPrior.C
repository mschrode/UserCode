//!  \brief Example to test the pt true prior for events cut
//!  in measured dijet pt.

#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"

namespace testPrior {
  
  //  TYPE DEFINITIONS
  //==================================================================

  //!  \brief Stores global parameter values
  //------------------------------------------------------------------
  struct Parameter {
    double ptGenMin;
    double ptGenMax;
    double ptDijetMin;
    double ptDijetMax;
    double expo;
    double stochTerm;
    double pdfMin;
    double pdfMax;
  };


  //!  \brief A dijet event
  //------------------------------------------------------------------
  class DijetEvent {
  public:
    DijetEvent() : ptJet1_(0.), ptJet2_(0.), ptGen_(0.) {};
    DijetEvent(double ptGen, double ptJet1, double ptJet2)
      : ptJet1_(ptJet1), ptJet2_(ptJet2), ptGen_(ptGen) {};
  
    double ptDijet() const { return 0.5 * ( ptJet1() + ptJet2() ); }
    double ptJet1() const { return ptJet1_; }
    double ptJet2() const { return ptJet2_; }
    double ptGen() const { return ptGen_; }

  private:
    double ptJet1_;
    double ptJet2_;
    double ptGen_;
  };
  typedef vector<DijetEvent> Data;
  typedef vector<DijetEvent>::const_iterator DataIt;



  //  FUNCTION DECLARATIONS
  //==================================================================
  TF1 * createPDF(const TString& name, const Parameter& par, bool cut);
  Data generateData(int nEvents, const Parameter& par);
  void plotData(const Data& data, const Parameter& par);
  double pdfTruth(double * x, double * par);



  //  FUNCTION IMPLEMENTATIONS
  //==================================================================



  //!  \brief Run the example
  //------------------------------------------------------------------
  void run() {
    Parameter pars;
    pars.ptGenMin = 100.;
    pars.ptGenMax = 800.;
    pars.ptDijetMin = 200.;
    pars.ptDijetMax = 400.;
    pars.expo = 7.;
    pars.stochTerm = 1.25;
    pars.pdfMin = 50.;
    pars.pdfMax = 1000.;
 
    Data data = generateData(500000,pars);
    plotData(data,pars);
  }



  //!  \brief Generate dijet events
  //!
  //!  Generate \p nEvents dijet events with true pt between
  //!  \p Parameter::ptGenMin and infinity according to a powerlaw
  //!  spectrum falling like \f$ 1/x^{n} \f$ with \p n being
  //!  \p Parameter::expo. Two measured pt, the pt of the measured
  //!  jets, are generated from the truth by smearing with a Gaussian
  //!  around the true pt with a width \p Parameter::stochTerm.
  //------------------------------------------------------------------
  Data generateData(int nEvents, const Parameter& par) {
    cout << "Generating " << nEvents << " dijet events\n";
    Data data = vector<DijetEvent>(nEvents);
    TRandom3 * rand = new TRandom3(0);
    for(int i = 0; i < nEvents; i++) {
      double ptGen = rand->Uniform(0.,1.);
      ptGen = par.ptGenMin * pow(ptGen,-1.0/(par.expo-1.));
      double ptJet1 = rand->Gaus( ptGen, par.stochTerm*sqrt(ptGen) );
      double ptJet2 = rand->Gaus( ptGen, par.stochTerm*sqrt(ptGen) );
      data.at(i) = DijetEvent(ptGen,ptJet1,ptJet2);
    }

    return data;
  }



  //!  \brief Plot the dijet pt spectrum and prediction
  //!
  //!  Plots
  //!  - Dijet pt versus true pt
  //!  - Dijet pt spectrum (normalized)
  //!  - True pt spectrum (normalized)
  //!  for all events in \p data and for all events with dijet pt
  //!  between \p Parameters::ptDijetMin and \p Parameters::ptDijetMax.
  //!  probability density of the truth spectrum is superimposed.
  //------------------------------------------------------------------
  void plotData(const Data& data, const Parameter& par) {
    cout << "Plotting data\n";
    gStyle->SetOptStat(0);

    vector<TH2D*> hPtCalVsPtGen(2);
    vector<TH1D*> hPtCal(2);
    vector<TH1D*> hPtGen(2);

    for(int i = 0; i < 2; i++) {
      TString name = "hPtCalVsPtGen";
      name += i;
      TString title = "All events";
      if( i == 1 ) {
	title = "p^{dijet}_{T} > ";
	title += par.ptDijetMin;
	title += " GeV";
      }
      hPtCalVsPtGen.at(i) = new TH2D(name,";p^{true}_{T} (GeV);p^{dijet}_{T} (GeV)",
				     50,0.5*par.ptGenMin,1.1*par.ptGenMax,
				     50,0.5*par.ptGenMin,1.1*par.ptGenMax);
      hPtCalVsPtGen.at(i)->SetTitle(title);

      name = "hPtCal";
      name += i;
      hPtCal.at(i) = new TH1D(name,";p^{dijet}_{T} (GeV)",50,0.5*par.ptGenMin,1.1*par.ptGenMax);
      hPtCal.at(i)->SetTitle(title);
      hPtCal.at(i)->SetLineWidth(2);
      hPtCal.at(i)->Sumw2();

      name = "hPtGen";
      name += i;
      hPtGen.at(i) = new TH1D(name,";p^{true}_{T} (GeV)",50,0.5*par.ptGenMin,1.1*par.ptGenMax);
      hPtGen.at(i)->SetTitle(title);
      hPtGen.at(i)->SetLineWidth(2);
      hPtGen.at(i)->Sumw2();
    }

    for(DataIt evt = data.begin(); evt != data.end(); evt++) {
      hPtCalVsPtGen.at(0)->Fill(evt->ptGen(),evt->ptDijet());
      hPtCal.at(0)->Fill(evt->ptDijet());
      hPtGen.at(0)->Fill(evt->ptGen());
      if( evt->ptDijet() > par.ptDijetMin && evt->ptDijet() < par.ptDijetMax ) {
	hPtCalVsPtGen.at(1)->Fill(evt->ptGen(),evt->ptDijet());
	hPtCal.at(1)->Fill(evt->ptDijet());
	hPtGen.at(1)->Fill(evt->ptGen());
      }
    }

    for(int i = 0; i < 2; i++) {
      if( hPtCal.at(i)->Integral("width") ) hPtCal.at(i)->Scale(1./hPtCal.at(i)->Integral("width"));
      if( hPtGen.at(i)->Integral("width") ) hPtGen.at(i)->Scale(1./hPtGen.at(i)->Integral("width"));
    }

    TF1 * pdfNoCut = createPDF("pdfNoCut",par,false);
    TF1 * pdfCut = createPDF("pdfCut",par,true);

    TCanvas * can1 = new TCanvas("can1","Meas vs truth",1000,500);
    can1->Divide(2,1);
    TCanvas * can2 = new TCanvas("can2","Meas",1000,500);
    can2->Divide(2,1);
    TCanvas * can3 = new TCanvas("can3","Truth",1000,500);
    can3->Divide(2,1);

    for(int i = 0; i < 2; i++) {
      can1->cd(1+i);
      hPtCalVsPtGen.at(i)->Draw("BOX");

      can2->cd(1+i);
      hPtCal.at(i)->Draw();
      can2->cd(1+i)->SetLogy();

      can3->cd(1+i);
      hPtGen.at(i)->Draw();
      if( i == 0 ) {
	pdfNoCut->Draw("same");
      }
      else if( i == 1 ) pdfCut->Draw("same");
      can3->cd(1+i)->SetLogy();
    }
  }



  //!  \brief Returns probability density of the truth spectrum
  //!
  //!  If \p cut is false, the probability density for the complete
  //!  truth spectrum as generated by \p generateData() is returned;
  //!  else the probability density for the truth for only the events
  //!  with dijet pt between \p Parameters::ptDijetMin and 
  //!  \p Parameters::ptDijetMax.
  //------------------------------------------------------------------
  TF1 * createPDF(const TString& name, const Parameter& par, bool cut) {
    TF1 * pdf = new TF1(name,pdfTruth,par.pdfMin,par.pdfMax,9);
    if( cut ) pdf->SetParameter(0,1.);
    else pdf->SetParameter(0,0.);
    pdf->SetParameter(1,1.);
    pdf->SetParameter(2,par.stochTerm);
    pdf->SetParameter(3,par.pdfMin);
    pdf->SetParameter(4,par.pdfMax);
    pdf->SetParameter(5,par.ptDijetMin);
    pdf->SetParameter(6,par.ptDijetMax);
    pdf->SetParameter(7,par.expo);
    pdf->SetParameter(8,par.ptGenMin);

    if( cut ) {
      double norm = pdf->Integral(par.pdfMin,par.pdfMax);
      pdf->SetParameter(1,1./norm);
    }
    pdf->SetLineWidth(2);

    return pdf;
  }



  //!  \brief The actual function used by \p createPDF()
  //------------------------------------------------------------------
  double pdfTruth(double * x, double * par) {
    double pt    = x[0];
    double stochTerm = par[2];
    double pdfMin = par[3];
    double pdfMax = par[4];
    double ptDijetMin = par[5];
    double ptDijetMax = par[6];
    double expo = par[7];
    double ptGenMin = par[8];
    double sigma = stochTerm * sqrt(pt) / sqrt(2);

    double weight = 1.;
    if( par[0] == 1. ) {
      TF1 fResolution("fResolution","gaus",pdfMin,pdfMax);
      fResolution.SetParameter(0,1./sqrt(2*M_PI)/sigma);
      fResolution.SetParameter(1,pt);
      fResolution.SetParameter(2,sigma);
      weight = fResolution.Integral(ptDijetMin,ptDijetMax);
    }

    double prob = 0.;
    if( pt > ptGenMin ) {
      prob = 1. / pow(pt,expo);
    }

    double norm = (expo-1.) / pow(ptGenMin,1.-expo);

    return par[1] * norm * weight * prob;
  }
}





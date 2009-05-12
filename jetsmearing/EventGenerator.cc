// $Id: EventGenerator.cc,v 1.6 2009/05/12 16:27:51 mschrode Exp $

#include "EventGenerator.h"

#include <cmath>
#include <iostream>

#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "Jet.h"
#include "NJetEvent.h"
#include "PhotonJetEvent.h"

namespace js
{
  // --------------------------------------------------
  EventGenerator::EventGenerator()
    : mRandom(new TRandom3(0)),
      mModelTruth("Uniform"),
      mModelResp("Gauss+Uniform"),
      mHistResp(0)
  {
    mParTruth.push_back(100.);
    mParTruth.push_back(1000.);
    mParResp.push_back(1.2);
    mParResp.push_back(0.1);
    mParResp.push_back(0.05);
  }



  // --------------------------------------------------
  EventGenerator::EventGenerator(const std::string& modelTruth, const std::vector<double>& parTruth,
				 const std::string& modelResp, const std::vector<double>& parResp)
    : mRandom(new TRandom3(0)),
      mModelTruth(modelTruth),
      mParTruth(parTruth),
      mModelResp(modelResp),
      mParResp(parResp),
      mHistResp(0)
  {
    assert( mParTruth.size() == 2 );

    assert( mModelResp == "Histogram"
	    || mModelResp == "Gauss+Uniform"
	    || mModelResp == "TwoGauss" );

    if( mModelResp == "Gauss+Uniform" ) assert( mParResp.size() == 3 );
    if( mModelResp == "TwoGauss" )
      {
	assert( mParResp.size() == 5 );

	// Set up sum of two Gaussians as pdf
	double c  = mParResp.at(0);  // Normalization
	double u0 = 1.;              // Mean of central Gaussian (scale)
	double s0 = mParResp.at(1);  // Width of central Gaussian
	double u1 = mParResp.at(2);  // Mean of second Gaussian
	double s1 = mParResp.at(3);  // Width of central Gaussian
	// mParResp.at(4) is temeperature for phi variation

	double minResp = 0.;
	double maxResp = 2.;
	TF1 * f = new TF1("f","gaus(0)+gaus(3)",minResp,maxResp);
	f->SetParameter(0,c);
	f->SetParameter(1,u0);
	f->SetParameter(2,s0);
	f->SetParameter(3,1.-c);
	f->SetParameter(4,u1);
	f->SetParameter(5,s1);

	// Fill response histogram according to f
	//	mHistResp = new TH1F("hHistResp","",200,minResp,maxResp);
	mHistResp = new TH1F("hHistResp",";p^{jet}_{T} / p^{true}_{T};1/(Nw) dN / d(p^{jet}_{T} / p^{true}_{T})",
			     200,minResp,maxResp);
	for(int bin = 1; bin <= mHistResp->GetNbinsX(); bin++)
	  {
	    double r = f->Eval(mHistResp->GetBinCenter(bin));
	    mHistResp->SetBinContent(bin,r);
	  }
	double norm = mHistResp->Integral("width");
	mHistResp->Scale(1./norm);
	delete f;
      }

  }



  // --------------------------------------------------
  EventGenerator::~EventGenerator()
  {
    delete mRandom;
  }



  // --------------------------------------------------
  Data EventGenerator::GenerateDijetEvents(int n) const
  {
    std::cout << "Generating " << n << " dijet events... " << std::flush;

    Data data;

    // Generate n dijet events      
    for(int i = 0; i < n; i++)
      {
	//if( i%1000 == 0 ) std::cout << i << std::endl;

	// Generate truth of first jet
	// between min and max
	TLorentzVector pTrue = GenerateTruth();

	// Smear truth to get measurement
	// of first jet
	TLorentzVector pMeas = SmearTruth(pTrue);

	// Set up first jet
	js::Jet jet1(pTrue,pMeas);

	// Get truth of second jet
	pTrue.SetPtEtaPhiM(pTrue.Pt(),0,M_PI+(pTrue.Phi()),0);
	  
	// Smear truth to get measurement
	// of second jet
	//delete pMeas;
	pMeas = SmearTruth(pTrue);

	// Set up second jet
	Jet jet2(pTrue,pMeas);
	//	delete pTrue;
	//	delete pMeas;

	// Construct Dijet event
	DiJetEvent * evt = new DiJetEvent();
	evt->AddJet(jet1);
	evt->AddJet(jet2);
	//	evt->Print();

	// Add event to data
	data.push_back(evt);
      }

    std::cout << "ok" << std::endl;

    return data;
  }



  // --------------------------------------------------
  Data EventGenerator::GeneratePhotonJetEvents(int n) const
  {
    std::cout << "Generating " << n << " photon-jet events... " << std::flush;

    Data data;

    // Generate n photon-jet events      
    for(int i = 0; i < n; i++)
      {
	//if( i%1000 == 0 ) std::cout << i << std::endl;

	// Generate truth of jet
	// between min and max
	TLorentzVector pTrue = GenerateTruth();

	// Smear truth to get measurement
	// of jet
	TLorentzVector pMeas = SmearTruth(pTrue);

	// Set up jet
	js::Jet jet(pTrue,pMeas);

	// Get antiparallel momentum of truth
	// to get photon momentum
	pTrue.SetPtEtaPhiM(pTrue.Pt(),0,M_PI+(pTrue.Phi()),0);
	  
	// Construct Photon-Jet event
	PhotonJetEvent * evt = new PhotonJetEvent(pTrue,jet);

	// Add event to data
	data.push_back(evt);
      }

    std::cout << "ok" << std::endl;

    return data;
  }



  // --------------------------------------------------
  TLorentzVector EventGenerator::GenerateTruth() const
  {
    double pt  = mRandom->Uniform(mParTruth.at(0),mParTruth.at(1));
    double phi = mRandom->Uniform(0.,2*M_PI);
    
    TLorentzVector pTrue;
    pTrue.SetPtEtaPhiM(pt,0,phi,0);
    
    return pTrue;
  }



  // --------------------------------------------------
  TLorentzVector EventGenerator::SmearTruth(const TLorentzVector& pTrue) const
    {
      TLorentzVector pMeas = pTrue;

      if( mModelResp == "Gauss+Uniform")
	SmearTruthGaussUniform(pMeas);
      else if( mModelResp == "Histogram")
	SmearTruthHistogram(pMeas);
      else if( mModelResp == "TwoGauss")
	SmearTruthTwoGauss(pMeas);

      return pMeas;
    }



  // --------------------------------------------------
  void EventGenerator::SmearTruthGaussUniform(TLorentzVector& p) const
  {
    // Gaussian smearing of pt
    double pt = 0.;
    do
      {
	pt = mRandom->Gaus(p.Pt(),mParResp.at(0)*sqrt(p.Pt()));
      } while(pt < 0.);

    // Non Gaussian tail
    if( mRandom->Uniform() < mParResp.at(1) )
      pt = mRandom->Uniform(0.,p.Pt());
      
    // Exponential smearing of phi
    double dPhi = mRandom->Exp(mParResp.at(2));
    double phi  = p.Phi();
    if( mRandom->Uniform() < 0.5 ) phi += dPhi;
    else                           phi -= dPhi;
    
    p.SetPtEtaPhiM(pt,0,phi,0);
  }



  // --------------------------------------------------
  void EventGenerator::SmearTruthHistogram(TLorentzVector& p) const
  {
    // Pt from response histogram
    double pt = mHistResp->GetRandom();
    assert( pt >= 0. );
    pt *= p.Pt();
      
    // Exponential smearing of phi
    double dPhi = mRandom->Exp(mParResp.at(0));
    double phi  = p.Phi();
    if( mRandom->Uniform() < 0.5 ) phi += dPhi;
    else                           phi -= dPhi;
    
    p.SetPtEtaPhiM(pt,0,phi,0);
  }



  // --------------------------------------------------
  void EventGenerator::SmearTruthTwoGauss(TLorentzVector& p) const
  {
    // Pt from response histogram
    double pt = mHistResp->GetRandom();
    assert( pt >= 0. );
    pt *= p.Pt();
      
    // Exponential smearing of phi
    double dPhi = mRandom->Exp(mParResp.at(4));
    double phi  = p.Phi();
    if( mRandom->Uniform() < 0.5 ) phi += dPhi;
    else                           phi -= dPhi;
    
    p.SetPtEtaPhiM(pt,0,phi,0);
  }



  //!  \brief Write internal response histogram to file
  // --------------------------------------------------
  void EventGenerator::WriteResponseHist(const std::string& name) const
  {
    std::cout << "Writing response histogram to file " << name << "... " << std::flush;

    // For the canvas:
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetCanvasDefH(800); //Height of canvas
    gStyle->SetCanvasDefW(800); //Width of canvas
    gStyle->SetCanvasDefX(0);   //Position on screen
    gStyle->SetCanvasDefY(0);

    // For the frame
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameBorderSize(1);
    gStyle->SetFrameFillColor(kBlack);
    gStyle->SetFrameFillStyle(0);
    gStyle->SetFrameLineColor(kBlack);
    gStyle->SetFrameLineStyle(0);
    gStyle->SetFrameLineWidth(1);

    // For the Pad:
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(kWhite);
    gStyle->SetPadGridX(false);
    gStyle->SetPadGridY(false);
    gStyle->SetGridColor(0);
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);

    // For the histo:
    gStyle->SetHistLineColor(kBlack);
    gStyle->SetHistLineStyle(0);
    gStyle->SetHistLineWidth(1);

    // For the statistics box:
    gStyle->SetOptStat(0);
    //  gStyle->SetOptStat("neMR");
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

    // For the legend
    gStyle->SetLegendBorderSize(1);

    // Margins:
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadLeftMargin(0.2);
    gStyle->SetPadRightMargin(0.04);

    // For the Global title:
    gStyle->SetOptTitle(1);
    gStyle->SetTitleFont(42,"");
    gStyle->SetTitleColor(1);
    gStyle->SetTitleTextColor(1);
    gStyle->SetTitleFillColor(0);
    gStyle->SetTitleFontSize(0.1);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleX(0.58);
    gStyle->SetTitleH(0.05);
    gStyle->SetTitleXOffset(0);
    gStyle->SetTitleYOffset(0);
    gStyle->SetTitleBorderSize(0);

    // For the axis titles:
    gStyle->SetTitleColor(1,"XYZ");
    gStyle->SetTitleFont(42,"XYZ");
    gStyle->SetTitleSize(0.04,"XYZ");
    gStyle->SetTitleXOffset(1.5);
    gStyle->SetTitleYOffset(2.0);

    // For the axis labels:
    gStyle->SetLabelColor(1,"XYZ");
    gStyle->SetLabelFont(42,"XYZ");
    gStyle->SetLabelOffset(0.007,"XYZ");
    gStyle->SetLabelSize(0.04,"XYZ");

    // For the axis:
    gStyle->SetAxisColor(1,"XYZ");
    gStyle->SetStripDecimals(kTRUE);
    gStyle->SetTickLength(0.03,"XYZ");
    gStyle->SetNdivisions(510,"XYZ");
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);

    TCanvas * c = new TCanvas("c","",500,500);
    c->cd();
    if( mHistResp )
      {
	mHistResp->UseCurrentStyle();
	mHistResp->Draw();
	c->SetLogy();
      }
    c->SaveAs(name.c_str());
    delete c;
    std::cout << "ok" << std::endl;
  }
}

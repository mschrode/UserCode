// $Id: EventGenerator.cc,v 1.5 2009/05/08 12:14:17 mschrode Exp $

#include "EventGenerator.h"

#include <cmath>
#include <iostream>

#include "TF1.h"
#include "TCanvas.h"

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
	mHistResp = new TH1F("hHistResp","",25,minResp,maxResp);
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
    TCanvas * c = new TCanvas("c","",500,500);
    c->cd();
    if( mHistResp ) mHistResp->Draw();
    c->SaveAs(name.c_str());
    delete c;
    std::cout << "ok" << std::endl;
  }
}

// $Id: $

#include "EventGenerator.h"

#include <cmath>
#include <iostream>

namespace js
{
  // --------------------------------------------------
  EventGenerator::EventGenerator()
    : mRandom(new TRandom3(0)),
      mModelTruth("Uniform"),
      mModelResp("Gauss+Uniform"),
      mHistResp(0),
      mMin(100.),
      mMax(1000.)
  {
    mParResp.push_back(1.2);
    mParResp.push_back(0.1);
    mParResp.push_back(0.05);
  }



  // --------------------------------------------------
  EventGenerator::EventGenerator(double min, double max,
				 const std::string& modelTruth, const std::vector<double>& parTruth,
				 const std::string& modelResp, const std::vector<double>& parResp)
    : mRandom(new TRandom3(0)),
      mModelTruth(modelTruth),
      mParTruth(parTruth),
      mModelResp(modelResp),
      mParResp(parResp),
      mHistResp(0),
      mMin(min),
      mMax(max) {}



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
	pTrue.SetPtEtaPhiM(-1.*(pTrue.Pt()),0,M_PI+(pTrue.Phi()),0);
	  
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
  TLorentzVector EventGenerator::GenerateTruth() const
  {
    double pt  = mRandom->Uniform(mMin,mMax);
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
    pt *= p.Pt();
      
    // Exponential smearing of phi
    double dPhi = mRandom->Exp(mParResp.at(0));
    double phi  = p.Phi();
    if( mRandom->Uniform() < 0.5 ) phi += dPhi;
    else                           phi -= dPhi;
    
    p.SetPtEtaPhiM(pt,0,phi,0);
  }

}

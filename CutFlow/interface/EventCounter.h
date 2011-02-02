#ifndef USERCODE_CUTFLOW_EVENTCOUNTER_H
#define USERCODE_CUTFLOW_EVENTCOUNTER_H

// -*- C++ -*-
//
// Package:    CutFlow
// Class:      EventCounter
// 
/**\class CutFlow EventCounter.h UserCode/CutFlow/src/EventCounter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Schroeder,,,
//         Created:  Wed Feb  2 18:03:02 CET 2011
// $Id$
//
//


// system include files
#include <memory>

#include "TH1.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

class EventCounter : public edm::EDAnalyzer {
 public:
  explicit EventCounter(const edm::ParameterSet&);
  ~EventCounter();
  
  
 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  TH1* hNumEvts_;
};

#endif

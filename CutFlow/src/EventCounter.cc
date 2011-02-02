// $Id$

#include "TH1D.h"
#include "TString.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "UserCode/CutFlow/interface/EventCounter.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EventCounter::EventCounter(const edm::ParameterSet& iConfig) {
   //now do what ever initialization is needed

  TString name = iConfig.getParameter<std::string>("FilterName");

  edm::Service<TFileService> fs;
  hNumEvts_ = fs->make<TH1D>("NumEvts"+name,"Numer of events after '"+name+"' filter",1,0.,1.);
}


EventCounter::~EventCounter() {
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
EventCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  hNumEvts_->Fill(0.);
}


// ------------ method called once each job just before starting event loop  ------------
void 
EventCounter::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EventCounter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventCounter);

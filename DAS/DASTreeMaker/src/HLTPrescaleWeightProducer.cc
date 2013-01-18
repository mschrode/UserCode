// system include files
#include <cmath>
#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

//
// class decleration
//

class HLTPrescaleWeightProducer: public edm::EDProducer {
public:
  explicit HLTPrescaleWeightProducer(const edm::ParameterSet&);
  ~HLTPrescaleWeightProducer();
  
private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(edm::Run&, edm::EventSetup const&);
  
  const std::string triggerName_;

  /// The instance of the HLTConfigProvider as a data member
  HLTConfigProvider hltConfig_;
};



HLTPrescaleWeightProducer::HLTPrescaleWeightProducer(const edm::ParameterSet& iConfig)
  : triggerName_(iConfig.getParameter<std::string>("HLTName")) {
  
  //register your products
  produces<double>();
}

HLTPrescaleWeightProducer::~HLTPrescaleWeightProducer() {
}

// ------------ method called to produce the data  ------------
void HLTPrescaleWeightProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  int prescale = 1;
  // trigger info
  edm::Handle< edm::TriggerResults > hltresults;
  if( iEvent.getByLabel(edm::InputTag("TriggerResults::HLT"),hltresults) ) {
    const edm::TriggerNames & trigNames = iEvent.triggerNames(*hltresults);
    for(unsigned int itrig = 0; itrig != hltresults->size(); ++itrig) {
      std::string trigName = trigNames.triggerName(itrig);
      if(trigName.find(triggerName_) != std::string::npos) {
	prescale = hltConfig_.prescaleValue(iEvent,iSetup,trigName);
	//std::cout << "    " << trigName << ": " << prescale << std::endl;
      }
    }
  }
  
  // put prescale into the Event
  std::auto_ptr<double> pOut(new double(prescale));
  iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
void HLTPrescaleWeightProducer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void HLTPrescaleWeightProducer::endJob() {
}

// ------------ method called at beginning for each run  ---------------------------------
void HLTPrescaleWeightProducer::beginRun(edm::Run& iRun, edm::EventSetup const& iSetup) {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTPrescaleWeightProducer);

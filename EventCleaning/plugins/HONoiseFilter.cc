/**\class HONoiseFilter HONoiseFilter.cc EventCleaning/plugins/HONoiseFilter.cc

 Description: filters events with jets that have a large HO-energy fraction

*/
//
// Original Author:  Matthias Schroeder,,,
//         Created:  Tue Feb 26 17:23:48 CET 2013
// $Id: $
//
//


#include <memory>

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"



class HONoiseFilter : public edm::EDFilter {
public:
  
  explicit HONoiseFilter(const edm::ParameterSet & iConfig);
  ~HONoiseFilter() {}
  
private:
  const bool taggingMode_;
  const edm::InputTag patJetsInputTag_;
  const double jetPtMin_;
  const double jetEtaMax_;
  const double maxHOEfrac_;
  
  virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);
};



HONoiseFilter::HONoiseFilter(const edm::ParameterSet & iConfig) :
  taggingMode_ (iConfig.getParameter<bool>("taggingMode")),
  patJetsInputTag_(iConfig.getParameter<edm::InputTag>("patJetsInputTag")),
  jetPtMin_(iConfig.getParameter<double>("jetPtMin")),
  jetEtaMax_(iConfig.getParameter<double>("jetEtaMax")),
  maxHOEfrac_(iConfig.getParameter<double>("maxHOEfrac")) {
  produces<bool>();
}


bool HONoiseFilter::filter(edm::Event & iEvent, const edm::EventSetup & iSetup) {
  // Filter decision
  bool pass = true;

  // PAT jets
  edm::Handle< edm::View<pat::Jet> > patJets;
  iEvent.getByLabel(patJetsInputTag_,patJets);
  if( patJets.isValid() ) {
    // Loop over all PAT jets
    for(edm::View<pat::Jet>::const_iterator j = patJets->begin();
	j != patJets->end(); ++j) {
      if( j->isPFJet() ) {	// HO fraction computed from PF candidates
	if( j->pt() > jetPtMin_ && std::abs(j->eta()) < jetEtaMax_ ) {
	  double hoe = 0.;
	  std::vector<reco::PFCandidatePtr> pfCands = j->getPFConstituents();
	  for(std::vector<reco::PFCandidatePtr>::const_iterator pfCandIt = pfCands.begin();
	      pfCandIt != pfCands.end(); ++pfCandIt) {
	    hoe += (*pfCandIt)->hoEnergy();
	  }
	  if( j->energy() ) {
	    if( hoe / j->energy() > maxHOEfrac_ ) {
	      pass = false;
	      break;
	    }
	  }
	}
      }
    }
  }

  // Store filter decision in event
  iEvent.put( std::auto_ptr<bool>(new bool(pass)) );
  
  // Depending on mode, return filter decision or true
  return taggingMode_ || pass; 
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(HONoiseFilter);

// -*- C++ -*-
//
// Package:    CorrJetProducer
// Class:      CorrJetProducer
// 
/**\class CorrJetProducer CorrJetProducer.cc UserCode/CorrJetProducer/src/CorrJetProducer.cc

 Description: Produce collection of corrected jets from input jet collection

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Schroeder
//         Created:  Tue Aug  9 21:13:57 CEST 2011
// $Id: CorrJetProducer.cc,v 1.1 2011/08/10 11:53:20 mschrode Exp $
//
//


// system include files
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/PtComparator.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"


//
// class declaration
//

class CorrJetProducer : public edm::EDProducer {
   public:
      explicit CorrJetProducer(const edm::ParameterSet&);
      ~CorrJetProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
  const edm::InputTag src_;
  const std::string jecLevel_;
  const std::string instanceName_;
};


//
// constructors and destructor
//
CorrJetProducer::CorrJetProducer(const edm::ParameterSet& iConfig)
  : src_(iConfig.getParameter<edm::InputTag>("src")),
    jecLevel_(iConfig.getParameter<std::string>("jecLevel")),
    instanceName_(iConfig.getParameter<std::string>("instanceName"))
{
  produces< std::vector<pat::Jet> >(instanceName_);
}


CorrJetProducer::~CorrJetProducer()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
CorrJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Get the original jet collection from the event
  edm::Handle< edm::View<pat::Jet> > origJetsHandle;
  iEvent.getByLabel(edm::InputTag(src_),origJetsHandle);
  edm::View<pat::Jet> origJets = *origJetsHandle;
  
  // This will store the corrected jets
  std::auto_ptr< std::vector<pat::Jet> > correctedJets(new std::vector<pat::Jet>);

  // Loop over original jet collection, retrieve JEC factor from the
  // jets 4-momenta, and store corrected jet. Details at
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrOnTheFly
  const JetCorrector* corrector = JetCorrector::getJetCorrector(jecLevel_,iSetup);
  for(edm::View<pat::Jet>::const_iterator j = origJets.begin(); j != origJets.end(); ++j) {
    pat::Jet correctedJet = *j;
    double corrFactor = corrector->correction(correctedJet.p4());
    //    std::cout << correctedJet.pt() << std::flush;
    correctedJet.scaleEnergy(corrFactor);
    //    std::cout << " --> " << correctedJet.pt() << std::endl;
    correctedJets->push_back(correctedJet);
  }

  // Sort corrected jets by pt
  GreaterByPt<pat::Jet> ptComparator;
  std::sort(correctedJets->begin(),correctedJets->end(),ptComparator);

  // Store collection of corrected jets in the event
  iEvent.put(correctedJets,instanceName_);
}

// ------------ method called once each job just before starting event loop  ------------
void 
CorrJetProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CorrJetProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
CorrJetProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CorrJetProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CorrJetProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CorrJetProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CorrJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CorrJetProducer);

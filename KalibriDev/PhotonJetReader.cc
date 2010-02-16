//
//  $Id: PhotonJetReader.cc,v 1.2 2010/01/21 16:49:20 mschrode Exp $
//
#include "PhotonJetReader.h"

#include "CalibData.h"
#include "SmearPhotonJet.h"
#include "JetTruthEvent.h"
#include "Jet.h"
#include "JetWithTowers.h"
#include "ConfigFile.h"
#include "ToyMC.h"
#include "Parameters.h"
#include "GammaJetSel.h"
#include "CorFactors.h"
#include "CorFactorsFactory.h"

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "TVector2.h"
#include "TLorentzVector.h"

PhotonJetReader::PhotonJetReader(const std::string& configfile, TParameters* p) :
  EventReader(configfile,p), gammaJet_(new GammaJetSel())
{
  // Maximum number of read events
  nGammaJetEvents_     = config_->read<int>("use Gamma-Jet events",-1); 
  if(nGammaJetEvents_ == 0) return;

  // Cuts
  minJetEt_           = config_->read<double>("Et min cut on jet",0.0); 
  maxJetEt_           = config_->read<double>("Et max cut on jet",0.0); 
  minGammaEt_         = config_->read<double>("Et min cut on gamma",0.0);
  maxGammaEt_         = config_->read<double>("Et max cut on gamma",0.0);
  maxRel2ndJetEt_     = config_->read<double>("Relative n+1 Jet Et Cut",0.1);
  minGenJetEt_        = config_->read<double>("Et genJet min",0.0);
  maxGenJetEt_        = config_->read<double>("Et genJet max",10000.0);
  maxDeltaR_          = config_->read<double>("DeltaR cut on jet matching",0.25);
  maxJetEta_          = config_->read<double>("Eta cut on jet",5.0);
  minJetHadFraction_  = config_->read<double>("Min had fraction",0.07);
  maxJetHadFraction_  = config_->read<double>("Max had fraction",0.95);
  // Counter for cutflow
  nMinJetEt_          = 0;
  nMaxJetEt_          = 0;
  nMinGammaEt_        = 0;
  nMaxGammaEt_        = 0;
  nMaxRel2ndJetEt_    = 0;
  nMinGenJetEt_       = 0;    
  nMaxGenJetEt_       = 0;     
  nMaxDeltaR_         = 0;
  nMaxJetEta_         = 0;  
  nMinJetHadFraction_ = 0;     
  nMaxJetHadFraction_ = 0;     

  // Data class
  dataClass_ = config_->read<int>("Gamma-Jet data class", 1);
  if( !( dataClass_ == 1 || dataClass_ == 2 || dataClass_ == 5) ) {
    std::cerr << "PhotonJetReader: Unknown data class " << dataClass_ << ". Aborting!" << std::endl;
    exit(9);
  }

  // Input files
  gammaJet_->Init(createTree("gammaJet"));
}

PhotonJetReader::~PhotonJetReader()
{
}

int PhotonJetReader::readEvents(std::vector<Event*>& data)
{
  if(nGammaJetEvents_ == 0) return 0;

  // Reset counters of rejected events
  nMinJetEt_          = 0;
  nMaxJetEt_          = 0;
  nMinGammaEt_        = 0;
  nMaxGammaEt_        = 0;
  nMaxRel2ndJetEt_    = 0;
  nMinGenJetEt_       = 0;    
  nMaxGenJetEt_       = 0;     
  nMaxDeltaR_         = 0;
  nMaxJetEta_         = 0;  
  nMinJetHadFraction_ = 0;     
  nMaxJetHadFraction_ = 0;     

  int nevent    = gammaJet_->fChain->GetEntries();  // Number of events in chain
  int nReadEvts = 0;                              // Number of read events
  int nGoodEvts = 0;                              // Number of events passing all cuts

  // Some informative output for the interested calibrator
  // Check of correct data class
  std::cout << "Reading events of type ";
  if((dataClass_ == 1)  || (dataClass_ == 2)) {
    std::cout << "'JetTruthEvent'";
  } else if(dataClass_ == 5) {
    std::cout << "'SmearData'";
    if( !correctToL3_ && !correctL2L3_ ) {
      std::cerr << "WARNING: Jets are not corrected! Aborting\n";
      exit(9);
    }
  } else {
    std::cerr << "Unknown data class " << dataClass_ << '\n';
    exit(9);
  }
  std::cout << " (data class " << dataClass_ << "):\n";

  // Read events
  for (int i=0;i<nevent;i++) {
    nReadEvts++;

    if((i+1)%10000==0) std::cout << "  " << i+1 << std::endl;
    gammaJet_->fChain->GetEvent(i); 
    if (gammaJet_->NobjTowCal>200) {
      std::cerr <<"ERROR: Increase array sizes in GammaJetSelector; NobjTowCal="
		<<gammaJet_->NobjTowCal<<"!"<< std::endl;
      exit(8);
    }
 
    // Trivial cuts
    bool goodEvent = true;
    if( gammaJet_->PhotonEt < minGammaEt_ ) {
      nMinGammaEt_++;
      goodEvent = false;
    } else if( gammaJet_->PhotonEt > maxGammaEt_ ) {
      nMaxGammaEt_++;
      goodEvent = false;
    } else if( gammaJet_->JetCalEt < minJetEt_ ) {
      nMinJetEt_++;
      goodEvent = false;
    } else if( gammaJet_->JetCalEt > maxJetEt_ ) {
      nMaxJetEt_++;
      goodEvent = false;
    } else if( gammaJet_->NonLeadingJetPt / gammaJet_->PhotonPt > maxRel2ndJetEt_) {
      nMaxRel2ndJetEt_++;
      goodEvent = false;
    } else if( std::abs(gammaJet_->JetCalEta) > maxJetEta_ ) {
      nMaxJetEta_++;
      goodEvent = false;
    } else if( gammaJet_->JetGenEt < minGenJetEt_ ) {
      nMinGenJetEt_++;
      goodEvent = false;
    } else if( gammaJet_->JetGenEt > maxGenJetEt_ ) {
      nMaxGenJetEt_++;
      goodEvent = false;
    } else if( pow(gammaJet_->JetCalEta - gammaJet_->JetGenEta,2)
	       + pow(TVector2::Phi_mpi_pi(gammaJet_->JetCalPhi - gammaJet_->JetGenPhi),2)
	       > pow(maxDeltaR_,2) ) {
      nMaxDeltaR_++;
      goodEvent = false;
    }

    if( goodEvent ) {
      
      Event*                                  ev = 0;
      if(dataClass_ == 1 || dataClass_ == 2)  ev = createJetTruthEvent();
      else if(dataClass_ == 5)                ev = createSmearEvent();
      
      if(ev) {
	data.push_back(ev); 
	nGoodEvts++;
      }
    }

    if((nGammaJetEvents_ >= 0) && (nReadEvts >= nGammaJetEvents_)) break;
  }

  // Print cut flow
  std::cout << "Read " << nReadEvts << " gamma-jet events:\n";
  std::cout << "  " << (nReadEvts-=nMinGammaEt_) << std::flush;
  std::cout << " gamma-jet events photon Et > " << minGammaEt_ << "\n";
  std::cout << "  " << (nReadEvts-=nMaxGammaEt_) << std::flush;
  std::cout << " gamma-jet events photon Et < " << maxGammaEt_ << "\n";
  std::cout << "  " << (nReadEvts-=nMinJetEt_) << std::flush;
  std::cout << " gamma-jet events jet Et > " << minJetEt_ << "\n";
  std::cout << "  " << (nReadEvts-=nMaxJetEt_) << std::flush;
  std::cout << " gamma-jet events jet Et < " << maxJetEt_ << "\n";
  std::cout << "  " << (nReadEvts-=nMaxRel2ndJetEt_) << std::flush;
  std::cout << " gamma-jet events (non-leading jet Et) / (photon Et) > " << maxRel2ndJetEt_ << "\n";
  std::cout << "  " << (nReadEvts-=nMaxJetEta_) << std::flush;
  std::cout << " gamma-jet events with |eta| < " << maxJetEta_ << "\n";
  std::cout << "  " << (nReadEvts-=nMinGenJetEt_) << std::flush;
  std::cout << " gamma-jet events with ptgen > " << minGenJetEt_ << "\n";
  std::cout << "  " << (nReadEvts-=nMaxGenJetEt_) << std::flush;
  std::cout << " gamma-jet events with ptgen < " << maxGenJetEt_ << "\n";
  std::cout << "  " << (nReadEvts-=nMaxDeltaR_) << std::flush;
  std::cout << " gamma-jet events with DeltaR < " << maxDeltaR_ << "\n";
  std::cout << "  " << (nReadEvts-=nMinJetHadFraction_) << std::flush;
  std::cout << " gamma-jet events with hadronic fraction > " << minJetHadFraction_ << "\n";
  std::cout << "  " << (nReadEvts-=nMaxJetHadFraction_) << std::flush;
  std::cout << " gamma-jet events with hadronic fraction < " << maxJetHadFraction_ << "\n";
  std::cout << "Stored " << nGoodEvts << " gamma-jet events for analysis.\n";

  return nGoodEvts;
}



// ----------------------------------------------------------------   
Event* PhotonJetReader::createJetTruthEvent()
{
  double em        = 0;
  double had       = 0;
  double out       = 0;
  double dR        = 10;
  int closestTower = 0;

  TLorentzVector LJet(0,0,0,0);
  LJet.SetPtEtaPhiE(gammaJet_->JetCalPt,gammaJet_->JetCalEta,gammaJet_->JetCalPhi,gammaJet_->JetCalE);
  TLorentzVector LGenJet(0,0,0,0);
  LGenJet.SetPtEtaPhiE(gammaJet_->JetGenPt,gammaJet_->JetGenEta,gammaJet_->JetGenPhi,gammaJet_->JetGenE);

  // Loop over towers, find closest tower to jet axis,
  // and sum up emf, hadf, outf
  double seta = 0;
  double seta2 = 0;
  double sumpt = 0;
  for(int n = 0; n < gammaJet_->NobjTowCal; ++n) {
    em          += gammaJet_->TowEm[n];
    had         +=  gammaJet_->TowHad[n];
    out         +=  gammaJet_->TowOE[n]; 
 
    double phi =  gammaJet_->TowPhi[n];
    double eta =  gammaJet_->TowEta[n];
    seta += gammaJet_->TowEt[n] * eta;
    seta2 += gammaJet_->TowEt[n] * eta * eta;
    sumpt += gammaJet_->TowEt[n];
    double dphi  = TVector2::Phi_mpi_pi(gammaJet_->JetCalPhi-phi);
    double dr    = sqrt((gammaJet_->JetCalEta-eta)*(gammaJet_->JetCalEta-eta)+
			dphi*dphi);     
    if(dr < dR) {
      dR = dr;
      closestTower = n;
    }
  } // End of loop over towers
  
  // Cuts on hadronic fraction 
  if(had/(had + em) < minJetHadFraction_) {
    nMinJetHadFraction_++;
    return 0;
  } else if(had/(had + em) > maxJetHadFraction_) { 
    nMaxJetHadFraction_++;
    return 0;
  }

  double factor = gammaJet_->JetCalEt /  gammaJet_->JetCalE;
  double etaeta = sqrt(seta2/sumpt - seta * seta /(sumpt * sumpt));
  Jet *j;
  if(dataClass_ == 2) {
    JetWithTowers *jt = 
      new JetWithTowers(gammaJet_->JetCalEt,em * factor,had * factor,
			out * factor,gammaJet_->JetCalE,gammaJet_->JetCalEta,
			gammaJet_->JetCalPhi,etaeta,Jet::uds,gammaJet_->JetGenEt,
			LJet.DeltaR(LGenJet),createCorFactors(0),
			par_->jet_function(gammaJet_->TowId_eta[closestTower],
					   gammaJet_->TowId_phi[closestTower]),
			jet_error_param,par_->global_jet_function(),minJetEt_);
    for(int i = 0; i < gammaJet_->NobjTowCal; ++i) {
      double scale = gammaJet_->TowEt[i]/gammaJet_->TowE[i];
      jt->addTower(gammaJet_->TowEt[i],gammaJet_->TowEm[i]*scale,
		   gammaJet_->TowHad[i]*scale,gammaJet_->TowOE[i]*scale,
		   gammaJet_->TowE[i],gammaJet_->TowEta[i],gammaJet_->TowPhi[i],
		   par_->tower_function(gammaJet_->TowId_eta[i],gammaJet_->TowId_phi[i]),
		   tower_error_param);
    }
    j = jt;
  }
  else { 
    j = new Jet(gammaJet_->JetCalEt,em * factor,had * factor,out * factor,
		gammaJet_->JetCalE,gammaJet_->JetCalEta,gammaJet_->JetCalPhi,
		etaeta,Jet::uds,gammaJet_->JetGenEt,LJet.DeltaR(LGenJet),
		createCorFactors(0),
		par_->jet_function(gammaJet_->TowId_eta[closestTower],
				   gammaJet_->TowId_phi[closestTower]),
		jet_error_param,par_->global_jet_function(),minJetEt_);
  }
  if(corFactorsFactory_) j->updateCorFactors(corFactorsFactory_->create(j));
  if(correctToL3_) j->correctToL3();
  JetTruthEvent * jte = new JetTruthEvent(j,gammaJet_->PhotonEt,gammaJet_->EventWeight);
  
  return jte;
}



//!  \brief Create SmearPhotonJet event for jet smearing
//!  \note Measured pt is L2L3 corrected
// ----------------------------------------------------------------   
Event* PhotonJetReader::createSmearEvent()
{
  //Find the jets eta & phi index using the nearest tower to jet axis:
  double min_tower_dr = 10.;
  double em           = 0;
  double had          = 0;
  double out          = 0;
  int    closestTower = 0;

  TLorentzVector LJet(0,0,0,0);
  LJet.SetPtEtaPhiE(gammaJet_->JetCalEt,gammaJet_->JetCalEta,gammaJet_->JetCalPhi,gammaJet_->JetCalE);
  TLorentzVector LGenJet(0,0,0,0);
  LGenJet.SetPtEtaPhiE(gammaJet_->JetGenPt,gammaJet_->JetGenEta,gammaJet_->JetGenPhi,gammaJet_->JetGenE);

  for (int n=0; n<gammaJet_->NobjTowCal; ++n) {
    em  += gammaJet_->TowEm[n];
    had += gammaJet_->TowHad[n];
    out += gammaJet_->TowOE[n];
    TLorentzVector LTower(0,0,0,0);
    LTower.SetPtEtaPhiE(gammaJet_->TowEt[n],gammaJet_->TowEta[n],gammaJet_->TowPhi[n],gammaJet_->TowE[n]);
    double dr = LTower.DeltaR(LJet);
    if (dr<min_tower_dr) {
      min_tower_dr = dr;
      closestTower = n;
    }
  }

  // Cuts on hadronic fraction 
  if(had/(had + em) < minJetHadFraction_) {
    nMinJetHadFraction_++;
    return 0;
  } else if(had/(had + em) > maxJetHadFraction_) { 
    nMaxJetHadFraction_++;
    return 0;
  }

  // Set up measurement
  double projFac   = gammaJet_->JetCalEt /  gammaJet_->JetCalE;
  Jet *jet = new Jet(gammaJet_->JetCalEt,
			  em * projFac,
			  had * projFac,
			  out * projFac,
			  gammaJet_->JetCalE,
			  gammaJet_->JetCalEta,
			  gammaJet_->JetCalPhi,
			  0.,
			  Jet::uds,
			  gammaJet_->JetGenPt,
			  LJet.DeltaR(LGenJet),
			  createCorFactors(0),
			  par_->jet_function(gammaJet_->TowId_eta[closestTower],
					     gammaJet_->TowId_phi[closestTower]),
			  jet_error_param,
			  par_->global_jet_function());

  // Create smear event
  return new SmearPhotonJet(jet,gammaJet_->PhotonEt,1.,1.,par_->resolutionFitPDF(1,1));
}


CorFactors* PhotonJetReader::createCorFactors(int jetid) const
{
  return new CorFactors(gammaJet_->JetCorrZSP, // L1
			gammaJet_->JetCorrL2,  // L2
			gammaJet_->JetCorrL3,  // L3
			1.,                  // L4
			1.,                  // L5
			gammaJet_->JetCorrJPT,
			gammaJet_->JetCorrL2L3JPT);
}

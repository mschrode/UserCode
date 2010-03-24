//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: DiJetReader.cc,v 1.7 2010/02/25 15:28:18 mschrode Exp $
//   
#include "DiJetReader.h"

#include "CalibData.h"
#include "SmearDiJet.h"
#include "ConfigFile.h"
#include "Parameters.h"
#include "JetTruthEvent.h"
#include "Jet.h"
#include "JetWithTowers.h"
#include "TwoJetsPtBalanceEvent.h"
#include "NJetSel.h"
#include "CorFactors.h"
#include "CorFactorsFactory.h"
#include "Function.h"
#include "SmearFunction.h"

#include "TVector2.h"
#include "TRandom3.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>



//!  \brief Constructor
//!
//!  Reads data from ROOT trees and stores them in an \p NJetSel selector.
//!  The data can be stored in a format derived from \p Event (as specified
//!  in the 'Di-Jet data class' field in the config file) by calling the
//!  method readEvents(<tt>std::vector<Event*>& data</tt>). Additionally,
//!  the cut thresholds are read from the configfile.
//!
//!  \param configfile Name of configfile
//!  \param p Pointer to \p TParameters object
// ----------------------------------------------------------------   
DiJetReader::DiJetReader(const std::string& configfile, TParameters* p)
  : EventReader(configfile,p), nJet_(new NJetSel()), rand_(new TRandom3(0))
{
  // Maximum number of read events
  nDijetEvents_ = config_->read<int>("use Di-Jet events",-1);
  if(nDijetEvents_ == 0) return;
  prescale_ = config_->read<int>("Di-Jet prescale",1);

  // Cuts
  minJetEt_          = config_->read<double>("Et min cut on jet",0.0);
  minDijetEt_        = config_->read<double>("Et min cut on dijet",0.0); 
  maxDijetEt_        = config_->read<double>("Et max cut on dijet",100.0); 
  max3rdJetEt_       = config_->read<double>("Et cut on n+1 Jet",10.0);
  maxRel3rdJetEt_    = config_->read<double>("Relative n+1 Jet Et Cut",0.2);
  maxJetEta_         = config_->read<double>("Eta cut on jet",5.0);
  minJetHadFraction_ = config_->read<double>("Min had fraction",0.07);
  maxJetHadFraction_ = config_->read<double>("Max had fraction",0.95);
  minDeltaPhi_       = config_->read<double>("Min Delta Phi",2.5);
  minGenJetEt_       = config_->read<double>("Et genJet min",0.0);
  maxGenJetEt_       = config_->read<double>("Et genJet max",10000.0);
  maxDeltaR_         = config_->read<double>("DeltaR cut on jet matching",0.25);
  // Counter for cutflow
  nDiJetCut_          = 0;
  nMinJetEt_          = 0;
  nMinDijetEt_        = 0;
  nMaxDijetEt_        = 0;
  nCutOn3rdJet_       = 0;
  nMaxJetEta_         = 0;  
  nMinDeltaPhi_       = 0;
  nMinJetHadFraction_ = 0;     
  nMaxJetHadFraction_ = 0;     
  nMinGenJetEt_       = 0;    
  nMaxGenJetEt_       = 0;     
  nMaxDeltaR_         = 0;
  // Integration parameter for SmearData
  maxNIter_          = config_->read<int>("DiJet integration number of iterations",5);
  eps_               = config_->read<double>("DiJet integration epsilon",1.E-5);
  min_               = config_->read<double>("DiJet integration min",0.);
  max_               = config_->read<double>("DiJet integration max",1.);
  // Data class
  dataClass_ = config_->read<int>("Di-Jet data class", 0);
  if((dataClass_ != 1)&&(dataClass_ != 11)&&(dataClass_ != 12)&&(dataClass_ != 5)) {
    std::cerr << "DiJetReader: Unknown data class " << dataClass_ << ". Aborting!" << std::endl;
    exit(9);
  }
  // Input files
  nJet_->Init(createTree("dijet"));
}

DiJetReader::~DiJetReader() {
  delete rand_;
}



//!  \brief Read dijet data and store in format as specified
//!         in the config file
//!  \param data Read data objects are appended to data
//!  \return Number of appended objects
// ----------------------------------------------------------------   
int DiJetReader::readEvents(std::vector<Event*>& data)
{
  if(nDijetEvents_ == 0) return 0;

  // Reset counters of rejected events
  nDiJetCut_          = 0;
  nMinJetEt_          = 0;
  nMinDijetEt_        = 0;
  nMaxDijetEt_        = 0;
  nCutOn3rdJet_       = 0;
  nMinGenJetEt_       = 0;    
  nMaxGenJetEt_       = 0;     
  nMaxDeltaR_         = 0;
  nMaxJetEta_         = 0;  
  nMinJetHadFraction_ = 0;     
  nMaxJetHadFraction_ = 0;     
  nMinDeltaPhi_       = 0;

  //Run jet-Jet stuff  
  int nevent    = nJet_->fChain->GetEntries();  // Number of events in chain
  int nReadEvts = 0;                          // Number of read events
  int nGoodEvts = 0;                          // Number of events passing all cuts

  // Some informative output for the interested calibrator
  // Check of correct data class
  std::cout << "Reading events of type ";
  if(dataClass_ == 1) {
    std::cout << "'TwoJetsPtBalanceEvent'";
  } else if((dataClass_ == 11)  || (dataClass_ == 12)) {
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

  // Read the events
  for (int i=0 ; i < nevent ; i+= prescale_) {
    if((i+1)%10000==0) std::cout << "  " << i+1 << std::endl;
    nJet_->fChain->GetEvent(i); 
    if (nJet_->NobjTow>10000 || nJet_->NobjJet>100) {
      std::cerr << "ERROR: Increase array sizes in NJet_Selector; NobjTow="
		<< nJet_->NobjTow<<", NobjJet="<<nJet_->NobjJet<<"!"<<std::endl;
      exit(9);
    }
    if(dataClass_ == 1) {
      nReadEvts++;
      Event* td = createTwoJetsPtBalanceEvent(); 
      if(td) {
	nGoodEvts++;    
	data.push_back(td ); 
      } 
    } else if((dataClass_ == 11)  || (dataClass_ == 12)) {
      nReadEvts++;
      int nAddedJets = createJetTruthEvents(data);
      if( nAddedJets ) nGoodEvts += nAddedJets;    
    } else if(dataClass_ == 5) {
      for(int calls = 0; calls < 2; calls++) {
	nReadEvts++;
	Event* td = createSmearEvent(calls); 
	if(td) {
	  nGoodEvts++;
	  data.push_back(td ); 
	} 
      }
    }
    if(nReadEvts>=nDijetEvents_ && nDijetEvents_>=0 ) break;
  }

  // Print cut flow
  std::cout << "Read " << nReadEvts << " dijet events:\n";
  if( dataClass_ == 1 ) {
    std::cout << "  " << (nReadEvts-=nDiJetCut_) << std::flush;
    std::cout << " events with 2 or more jets\n";
    std::cout << "  " << (nReadEvts-=nMinGenJetEt_) << std::flush;
    std::cout << " dijet events with ptgen > " << minGenJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMaxGenJetEt_) << std::flush;
    std::cout << " dijet events with ptgen < " << maxGenJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMaxDeltaR_) << std::flush;
    std::cout << " dijet events with DeltaR < " << maxDeltaR_ << "\n";
    std::cout << "  " << (nReadEvts-=nMinJetEt_) << std::flush;
    std::cout << " dijet events Et > " << minJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMaxJetEta_) << std::flush;
    std::cout << " dijet events with |eta| < " << maxJetEta_ << "\n";
    std::cout << "  " << (nReadEvts-=nMinJetHadFraction_) << std::flush;
    std::cout << " dijet events with hadronic fraction > " << minJetHadFraction_ << "\n";
    std::cout << "  " << (nReadEvts-=nMaxJetHadFraction_) << std::flush;
    std::cout << " dijet events with hadronic fraction < " << maxJetHadFraction_ << "\n";
    std::cout << "  " << (nReadEvts-=nCutOn3rdJet_) << std::flush;
    std::cout << " dijet events with pt(jet3) / pt(dijet) < " << maxRel3rdJetEt_ << " or ";
    std::cout << "pt(jet3) < " << max3rdJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMinDeltaPhi_) << std::flush;
    std::cout << " dijet events with DeltaPhi > " << minDeltaPhi_ << "\n";
  } else if( dataClass_ == 11 || dataClass_ == 12 ) {
    std::cout << "  " << (nReadEvts-=nDiJetCut_) << std::flush;
    std::cout << " events with 2 or more jets\n";
    std::cout << "  That are " << (nReadEvts*=2) << " jet-truth events:\n";
    std::cout << "    " << (nReadEvts-=nMinGenJetEt_) << std::flush;
    std::cout << " jet-truth events with ptgen > " << minGenJetEt_ << "\n";
    std::cout << "    " << (nReadEvts-=nMaxGenJetEt_) << std::flush;
    std::cout << " jet-truth events with ptgen < " << maxGenJetEt_ << "\n";
    std::cout << "    " << (nReadEvts-=nMaxDeltaR_) << std::flush;
    std::cout << " jet-truth events with DeltaR < " << maxDeltaR_ << "\n";
    std::cout << "    " << (nReadEvts-=nMinJetEt_) << std::flush;
    std::cout << " jet-truth events Et > " << minJetEt_ << "\n";
    std::cout << "    " << (nReadEvts-=nMaxJetEta_) << std::flush;
    std::cout << " jet-truth events with |eta| < " << maxJetEta_ << "\n";
    std::cout << "    " << (nReadEvts-=nMinJetHadFraction_) << std::flush;
    std::cout << " jet-truth events with hadronic fraction > " << minJetHadFraction_ << "\n";
    std::cout << "    " << (nReadEvts-=nMaxJetHadFraction_) << std::flush;
    std::cout << " jet-truth events with hadronic fraction < " << maxJetHadFraction_ << "\n";
  } else if( dataClass_ == 5 ) {
    std::cout << "  " << (nReadEvts-=nDiJetCut_) << std::flush;
    std::cout << " events with more than 3 jets\n";
    std::cout << "  " << (nReadEvts-=nMinGenJetEt_) << std::flush;
    std::cout << " dijet events with ptgen > " << minGenJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMaxGenJetEt_) << std::flush;
    std::cout << " dijet events with ptgen < " << maxGenJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMaxDeltaR_) << std::flush;
    std::cout << " dijet events with DeltaR < " << maxDeltaR_ << "\n";
    std::cout << "  " << (nReadEvts-=nMinJetEt_) << std::flush;
    std::cout << " dijet events Et > " << minJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMaxJetEta_) << std::flush;
    std::cout << " dijet events with |eta| < " << maxJetEta_ << "\n";
    std::cout << "  " << (nReadEvts-=nMinDijetEt_) << std::flush;
    std::cout << " dijet events dijet pt > " << minDijetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMaxDijetEt_) << std::flush;
    std::cout << " dijet events dijet pt < " << maxDijetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMinJetHadFraction_) << std::flush;
    std::cout << " dijet events with hadronic fraction > " << minJetHadFraction_ << "\n";
    std::cout << "  " << (nReadEvts-=nMaxJetHadFraction_) << std::flush;
    std::cout << " dijet events with hadronic fraction < " << maxJetHadFraction_ << "\n";
    std::cout << "  " << (nReadEvts-=nCutOn3rdJet_) << std::flush;
    std::cout << " dijet events with pt(jet3) / pt(dijet) < " << maxRel3rdJetEt_ << " or ";
    std::cout << "pt(jet3) < " << max3rdJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMinDeltaPhi_) << std::flush;
    std::cout << " dijet events with DeltaPhi > " << minDeltaPhi_ << "\n";
  }
  std::cout << "Stored " << nGoodEvts << " dijet events for analysis.\n";

  return nGoodEvts;
}

//!  \brief Use first two jets of dijet event as two JetTruthEvent
//!         objects where genjet Et is truth
//!  \note The dijets are ordered in genjet Et
//!  \param data Read JetTruthEvent objects are appended to data
//!  \return Number of appended JetTruthEvent objects (0 - 2)
// ----------------------------------------------------------------   
int DiJetReader::createJetTruthEvents(std::vector<Event*>& data)
{
  int injet = 2;
  if( nJet_->NobjGenJet < injet ) {
    nDiJetCut_++;
    return 0;
  }
  int     njets = 0;   // Number of stored JetTruthEvents; 0 or 2
  double * terr = new double[nJet_->NobjTow];


  // Loop over two jets with highest genjet Et
  for(int genJetIdx = 0; genJetIdx < 2; genJetIdx++) {
    int calJetIdx = nJet_->GenJetColJetIdx[genJetIdx]; // Closest (DeltaR) calo jet
    if( nJet_->NobjJet <= calJetIdx ) {
      nDiJetCut_++;
      return 0;
    }

    // Cuts
    if( nJet_->GenJetColEt[genJetIdx] < minGenJetEt_ ) {
      nMinGenJetEt_++;
      continue;
    } else if( nJet_->GenJetColEt[genJetIdx] > maxGenJetEt_ ) {
      nMaxGenJetEt_++;
      continue;
    }
    double dphi        = TVector2::Phi_mpi_pi( nJet_->JetPhi[calJetIdx] - nJet_->GenJetColPhi[genJetIdx] );
    double deta        = nJet_->JetEta[calJetIdx] - nJet_->GenJetColEta[genJetIdx];
    double drJetGenjet = sqrt( deta*deta + dphi*dphi );
    if( drJetGenjet > maxDeltaR_ ) {
      nMaxDeltaR_++;
      continue;
    } else if( nJet_->JetEt[calJetIdx] < minJetEt_ ) {
      nMinJetEt_++;
      continue;
    } else if( std::abs(nJet_->JetEta[calJetIdx]) > maxJetEta_ ) {
      nMaxJetEta_++;
      continue;
    }

    // Construct event
    double em   = 0;
    double had  = 0;
    double out  = 0;
    double err2 = 0;
    Measurement tower;
    double dR        = 10;
    int closestTower = 0;
    double seta = 0;
    double seta2 = 0;
    double sumpt = 0;
    for(int n=0; n<nJet_->NobjTow; ++n){
      if(nJet_->Tow_jetidx[n] != calJetIdx) continue;//look for ij-jet's towers

      em          += nJet_->TowEm[n];
      had         += nJet_->TowHad[n];
      out         += nJet_->TowOE[n];  
      tower.pt     = nJet_->TowEt[n];
      double scale = nJet_->TowEt[n]/nJet_->TowE[n];
      tower.EMF    = nJet_->TowEm[n]*scale;
      tower.HadF   = nJet_->TowHad[n]*scale;
      tower.OutF   = nJet_->TowOE[n]*scale;
      tower.eta    = nJet_->TowEta[n];
      tower.phi    = nJet_->TowPhi[n];
      tower.E      = nJet_->TowE[n];
      terr[n]      = tower_error_param(&tower.pt,&tower,0); 
      if(terr[n] == 0) {
	//assume toy MC???
	terr[n] = TParameters::toy_tower_error_parametrization(&tower.pt,&tower);
      }
      terr[n]  *= terr[n];
      err2     += terr[n];
      dphi      = TVector2::Phi_mpi_pi(nJet_->JetPhi[calJetIdx]-tower.phi);
      deta      = nJet_->JetEta[calJetIdx]-tower.eta;
      seta += tower.pt * tower.eta;
      seta2 += tower.pt * tower.eta * tower.eta;
      sumpt += tower.pt;
      double dr = sqrt( deta*deta + dphi*dphi );     
      if(dr < dR) {
	dR = dr;
	closestTower = n;
      }
    }
    // Cuts on hadronic fraction 
    if(had/(had + em) < minJetHadFraction_) {
      nMinJetHadFraction_++;
      continue;
    } else if(had/(had + em) > maxJetHadFraction_) { 
      nMaxJetHadFraction_++;
      continue;
    }
    double factor = nJet_->JetEt[calJetIdx] /  nJet_->JetE[calJetIdx];
    tower.pt      = nJet_->JetEt[calJetIdx];
    tower.EMF     = em * factor;
    tower.HadF    = had * factor;
    tower.OutF    = out * factor;
    tower.eta     = nJet_->JetEta[calJetIdx];
    tower.phi     = nJet_->JetPhi[calJetIdx];
    tower.E       = nJet_->JetE[calJetIdx];
    tower.etaeta  = sqrt(seta2/sumpt - seta * seta /(sumpt * sumpt));
    double err    = jet_error_param(&tower.pt,&tower,0);
    err2         += err * err;

    Jet *jet;
    if(dataClass_ == 12) {
      JetWithTowers *jt = 
	new JetWithTowers(nJet_->JetEt[calJetIdx],em * factor,had * factor,
			  out * factor,nJet_->JetE[calJetIdx],nJet_->JetEta[calJetIdx],
			  nJet_->JetPhi[calJetIdx],tower.etaeta, Jet::uds,
			  nJet_->GenJetColEt[genJetIdx],drJetGenjet,
			  createCorFactors(calJetIdx),
			  par_->jet_function(nJet_->TowId_eta[closestTower],
					     nJet_->TowId_phi[closestTower]),
			  jet_error_param,par_->global_jet_function(),minJetEt_);
      for(int j = 0 ; j < nJet_->NobjTow ; ++j) {
	if (nJet_->Tow_jetidx[j]!= calJetIdx) continue;//look for ij-jet's towers
	double scale = nJet_->TowEt[j]/nJet_->TowE[j];
	jt->addTower(nJet_->TowEt[j],nJet_->TowEm[j]*scale,
		     nJet_->TowHad[j]*scale,nJet_->TowOE[j]*scale,
		     nJet_->TowE[j],nJet_->TowEta[j],nJet_->TowPhi[j],
		     par_->tower_function(nJet_->TowId_eta[calJetIdx],nJet_->TowId_phi[calJetIdx]),
		     tower_error_param);
      }
      jet = jt;
    }
    else { 
      jet = new Jet(nJet_->JetEt[calJetIdx],em * factor,had * factor,out * factor,
		    nJet_->JetE[calJetIdx],nJet_->JetEta[calJetIdx],nJet_->JetPhi[calJetIdx],tower.etaeta,
		    Jet::uds,nJet_->GenJetColEt[genJetIdx],drJetGenjet,
		    createCorFactors(calJetIdx),
		    par_->jet_function(nJet_->TowId_eta[closestTower],
				       nJet_->TowId_phi[closestTower]),
		    jet_error_param,par_->global_jet_function(),minJetEt_);    
    }
    if(corFactorsFactory_) {
      jet->updateCorFactors(corFactorsFactory_->create(jet));
    }
    if(correctToL3_) {
      jet->correctToL3();
    }
    if(correctL2L3_) {
      jet->correctL2L3();
    }
    JetTruthEvent* jte = new JetTruthEvent(jet,nJet_->GenJetColEt[genJetIdx],1.,nJet_->GenEvtScale);//nJet_->Weight);
    data.push_back(jte);
    ++njets;
  }     
  delete [] terr;
  return njets;
}



//!  \brief Create \p SmearDiJet event for jet smearing
//!
//!  Uses the three jets with the highest corrected
//!  calo pt. For creation of a \p SmearDiJet event,
//!  the JetMET L2L3 correction is applied.
//!
//!  \return A \p SmearDiJet if all cuts are passed,
//!          else 0
// ----------------------------------------------------------------   
Event* DiJetReader::createSmearEvent(int callIdx)
{
  // There should be at least three jets in the event
  if( nJet_->NobjJet < 3 ) {
    nDiJetCut_++;
    return 0;
  }

  // Read jets and apply L3 correction
  std::vector<Jet*> jets = readCaloJets(-1);

  // Sort jets by corrected pt
  std::sort(jets.begin(),jets.end(),Jet::caloPtGreaterThan);

  
//   if( nJet_->NobjGenJet < 3 ) {
//     nDiJetCut_++;
//     return 0;
//   }
//   if( nJet_->GenJetColPt[1] < minGenJetEt_ ) {
//     nMinGenJetEt_++;
//     return 0;
//   } else if( nJet_->GenJetColPt[0] > maxGenJetEt_ ) {
//     nMaxGenJetEt_++;
//     return 0;
//   }
//   // Read jets and apply L3 correction
//   std::vector<Jet*> jets = readGenJetSortedJets(3);

  
  // Create SmearDiJet event from three jets
  // with highest corrected pt
  int jetIdx[3] = { 0, 1, 2 };
//   if( rand_->Uniform() > 0.5 ) {
//     jetIdx[0] = 1;
//     jetIdx[1] = 0;
//   }
  if( callIdx == 1 ) {
    jetIdx[0] = 1;
    jetIdx[1] = 0;
  }
  SmearDiJet * smearDijet = new SmearDiJet(jets[jetIdx[0]],                    // First jet
					   jets[jetIdx[1]],                    // Second jet
					   jets[jetIdx[2]],                    // Third jet
					   nJet_->GenEvtScale,
					   1.,                         // Weights from EventProcessor
					   par_->resolutionFitPDF(1,1),
					   min_,                       // Integration minimum
					   max_,                       // Integration maximum
					   eps_,                       // Integration step length
					   maxNIter_);                 // Integration n iterations
  // Delete other jets
  for(std::vector<Jet*>::iterator jetIt = jets.begin()+3;
      jetIt != jets.end(); jetIt++) {
    delete *jetIt;
  }
  jets.erase(jets.begin()+3,jets.end());

  // Check if event is ok and return
  bool isGoodEvt = true;
  const Jet * j1 = smearDijet->jet1();
  const Jet * j2 = smearDijet->jet2();
  const Jet * j3 = smearDijet->jet3();

  if     ( j1->genPt() < minGenJetEt_ || j2->genPt() < minGenJetEt_ ) {
    nMinGenJetEt_++;
    isGoodEvt = false;
  }
  else if( j1->genPt() > maxGenJetEt_ || j2->genPt() > maxGenJetEt_ ) {
    nMaxGenJetEt_++;
    isGoodEvt = false;
  }
  else if( j1->dR() > maxDeltaR_ || j2->dR() > maxDeltaR_ ) {
    nMaxDeltaR_++;
    isGoodEvt = false;
  }
  else if( j1->pt() < minJetEt_ || j2->pt() < minJetEt_ ) {
    nMinJetEt_++;
    isGoodEvt = false;
  }
  else if( std::abs(j1->eta()) > maxJetEta_  || std::abs(j2->eta()) > maxJetEta_ ) {
    nMaxJetEta_++;
    isGoodEvt = false;
  }

//   else if( smearDijet->dijetPt() < minDijetEt_ ) {
//     nMinDijetEt_++;
//     isGoodEvt = false;
//   }
//   else if( smearDijet->dijetPt() > maxDijetEt_ ) {
//     nMaxDijetEt_++;
//     isGoodEvt = false;
//   }

  
  else if( j1->pt() < minDijetEt_ ) {
    nMinDijetEt_++;
    isGoodEvt = false;
  }
  else if( j1->pt() > maxDijetEt_ ) {
    nMaxDijetEt_++;
    isGoodEvt = false;
  }


  else if( j1->HadEt()/(j1->HadEt() + j1->EMF) < minJetHadFraction_ ||
	   j2->HadEt()/(j2->HadEt() + j2->EMF) < minJetHadFraction_ ) {
    nMinJetHadFraction_++;
    isGoodEvt = false;
  }
  else if( j1->HadEt()/(j1->HadEt() + j1->EMF) > maxJetHadFraction_ ||
	   j2->HadEt()/(j2->HadEt() + j2->EMF) > maxJetHadFraction_ ) {
    nMaxJetHadFraction_++;
    isGoodEvt = false;
  }
  else if( j3->pt() / smearDijet->dijetPt() > maxRel3rdJetEt_ && j3->pt() > max3rdJetEt_ ) {
    nCutOn3rdJet_++;
    isGoodEvt = false;
  }
//   else if( 2*j3->genPt() / (j1->genPt()+j2->genPt()) > maxRel3rdJetEt_ && j3->genPt() > max3rdJetEt_ ) {
//     nCutOn3rdJet_++;
//     isGoodEvt = false;
//   }

  else if( std::abs(TVector2::Phi_mpi_pi(j1->phi() - j2->phi())) < minDeltaPhi_ ) {
    nMinDeltaPhi_++;
    isGoodEvt = false;
  }


  //Hack: Reject tails
//   double sigma = 1.16*sqrt(j1->genPt());
//   if( std::abs(j1->pt() - j1->genPt()) > 3.5*sigma
//     || std::abs(j2->pt() - j2->genPt()) > 3.5*sigma  ) {
//     isGoodEvt = false;
//   }
  // End hack

  if(! isGoodEvt) {
    if( smearDijet ) delete smearDijet;
    smearDijet = 0;
  }
  return smearDijet;
}



//!  \brief Create an event of type \p TwoJetsPtBalanceEvent
//!
//!  The two jets leading in uncorrected calo pt are randomly
//!  assigned jet 1 and 2.
// ----------------------------------------------------------------   
Event* DiJetReader::createTwoJetsPtBalanceEvent()
{
  std::cerr << "DiJetReader::createTwoJetsPtBalanceEvent: needs to get fixed. Construction of Jet object wrong\n";
  exit(9);
  /*
  // Number of jets in the event
  int nJets = nJet_->NobjJet;

  // There should be at least two jets
  // and at the most three jets
  if( nJets < 2 ) {
  nDiJetCut_++;
  return 0;
  } else if( nJets >= 3 ) {
  nJets = 3;
  }

  // Pointer to the three jets leading in pt calo
  Jet * jet1 = 0;
  Jet * jet2 = 0;
  Jet * jet3 = 0;

  // Loop over the two or, if existing, three jets
  // leading in calo pt; the first two jets are
  // assigned randomly to have an unbiased balance
  // measure
  int calJetIdx[3];
  if( rand_->Uniform(0.,1.) < 0.5 ) {
  calJetIdx[0] = 0;
  calJetIdx[1] = 1;
  } else {
  calJetIdx[0] = 1;
  calJetIdx[1] = 0;
  }
  calJetIdx[2] = 2;
  for(int i = 0; i < nJets; i++) {
  // Jet - GenJet matching
  double dphi         = TVector2::Phi_mpi_pi( nJet_->JetPhi[calJetIdx[i]] - nJet_->GenJetPhi[calJetIdx[i]] );
  double deta         = nJet_->JetEta[calJetIdx[i]] - nJet_->GenJetEta[calJetIdx[i]];
  double drJetGenjet  = sqrt( deta*deta + dphi*dphi );

  // Summing up of hadronic fractions
  // and finding tower closest to the
  // jet axis
  double min_tower_dr = 10.;
  double emf          = 0;
  double had          = 0;
  double out          = 0;
  int    closestTower = 0; 

  // Loop over towers
  for(int n=0; n<nJet_->NobjTow; ++n) {
  if(nJet_->Tow_jetidx[n] != calJetIdx[i]) continue;//look for i-jet's towers
  emf += nJet_->TowEm[n];
  had += nJet_->TowHad[n];
  out += nJet_->TowOE[n];
  dphi = TVector2::Phi_mpi_pi( nJet_->JetPhi[calJetIdx[i]] - nJet_->TowPhi[n] );
  deta = nJet_->JetEta[calJetIdx[i]] - nJet_->TowEta[n];
  double dr = sqrt( deta*deta + dphi*dphi );     
  if (dr < min_tower_dr) {
  min_tower_dr = dr;
  closestTower = n;
  }
  } // End of loop over towers


    // Projection factor E --> Et
    // The following is not quite correct, as this factor is different for all towers
    // These values should be in the n-tupel as well
    double projFac   = nJet_->JetEt[calJetIdx[i]] /  nJet_->JetE[calJetIdx[i]];

    // Create jet
    Jet * jet = new Jet(nJet_->JetPt[calJetIdx[i]],
    emf * projFac,
    had * projFac,
    out * projFac,
    nJet_->JetE[calJetIdx[i]],
    nJet_->JetEta[calJetIdx[i]],
    nJet_->JetPhi[calJetIdx[i]],
    Jet::uds,
    nJet_->GenJetPt[calJetIdx[i]],
    drJetGenjet,createCorFactors(calJetIdx[i]),
    par_->jet_function(nJet_->TowId_eta[closestTower],
    nJet_->TowId_phi[closestTower]),
    jet_error_param,
    par_->global_jet_function(),
    minJetEt_ );    
    
    // Store jet
    if( i == 0 ) {
    jet1 = jet;
    } else if( i == 1 ) {
    jet2 = jet;
    } else if( i == 2 ) {
    jet3 = jet;
    }
    }  // End of loop over jets

    // Create TwoJetsInvMassEvent
    TwoJetsPtBalanceEvent * evt
    = new TwoJetsPtBalanceEvent(jet1,jet2,jet3,nJet_->GenEvtScale,nJet_->Weight) ;

    // Check if event is ok and return
    bool isGoodEvt = true;
    Jet * j1 = static_cast<Jet*>(evt->getJet1());
    Jet * j2 = static_cast<Jet*>(evt->getJet2());
    Jet * j3 = static_cast<Jet*>(evt->getJet3());
    if     ( j1->genPt < minGenJetEt_ || j2->genPt < minGenJetEt_ ) {
    nMinGenJetEt_++;
    isGoodEvt = false;
    }
    else if( j1->genPt > maxGenJetEt_ || j2->genPt > maxGenJetEt_ ) {
    nMaxGenJetEt_++;
    isGoodEvt = false;
    }
    else if( j1->dR > maxDeltaR_ || j2->dR > maxDeltaR_ ) {
    nMaxDeltaR_++;
    isGoodEvt = false;
    }
    else if( j1->pt < minJetEt_ || j2->pt < minJetEt_ ) {
    nMinJetEt_++;
    isGoodEvt = false;
    }
    else if( std::abs(j1->eta) > maxJetEta_ || std::abs(j2->eta) > maxJetEta_ ) {
    nMaxJetEta_++;
    isGoodEvt = false;
    }
    else if( j1->HadF/(j1->HadF + j1->EMF) < minJetHadFraction_ ||
    j2->HadF/(j2->HadF + j2->EMF) < minJetHadFraction_ ) {
    nMinJetHadFraction_++;
    isGoodEvt = false;
    }
    else if( j1->HadF/(j1->HadF + j1->EMF) > maxJetHadFraction_ ||
    j2->HadF/(j2->HadF + j2->EMF) > maxJetHadFraction_ ) {
    nMaxJetHadFraction_++;
    isGoodEvt = false;
    }
    else if( evt->relPtJet3() > maxRel3rdJetEt_ && j3->pt > max3rdJetEt_ ) {
    nCutOn3rdJet_++;
    isGoodEvt = false;
    }
    else if( std::abs(TVector2::Phi_mpi_pi(j1->phi - j2->phi)) < minDeltaPhi_ ) {
    nMinDeltaPhi_++;
    isGoodEvt = false;
    }

    if(! isGoodEvt) {
    if( evt ) delete evt;
    evt = 0;
    }

    return evt;
  */
}

CorFactors* DiJetReader::createCorFactors(int jetid) const 
{
  return new CorFactors(nJet_->JetCorrZSP[jetid], // L1
			nJet_->JetCorrL2[jetid],  // L2
			nJet_->JetCorrL3[jetid],  // L3
			1.,                         // L4
			1.,                         // L5
			nJet_->JetCorrJPT[jetid],
			nJet_->JetCorrL2L3JPT[jetid]); //JPTL2L3
}


std::vector<Jet*> DiJetReader::readCaloJets(int nJets) const {
  int maxNJets = nJets;
  if( maxNJets < 0 ) maxNJets = nJet_->NobjJet;
  std::vector<Jet*> caloJets(maxNJets);
  for(int j = 0; j < maxNJets; j++) {
    double dphi = nJet_->JetPhi[j] - nJet_->GenJetPhi[j];
    if( std::abs(dphi) > 1E-5 && std::abs(dphi) < 1E10 ) {
      dphi = TVector2::Phi_mpi_pi( dphi );
    }
    double deta         = nJet_->JetEta[j] - nJet_->GenJetEta[j];
    double drJetGenjet  = sqrt( deta*deta + dphi*dphi );
    double min_tower_dr = 10.;
    double emf          = 0;
    double had          = 0;
    double out          = 0;
    int    closestTower = 0; 
    // Loop over towers
    for (int n=0; n<nJet_->NobjTow; ++n) {
      if (nJet_->Tow_jetidx[n]!=(int)j) continue;//look for j-jet's towers
      emf += nJet_->TowEm[n];
      had += nJet_->TowHad[n];
      out += nJet_->TowOE[n];
      dphi = TVector2::Phi_mpi_pi( nJet_->JetPhi[j] - nJet_->TowPhi[n] );
      deta = nJet_->JetEta[j] - nJet_->TowEta[n];
      double dr = sqrt( deta*deta + dphi*dphi );     
      if (dr < min_tower_dr) {
	min_tower_dr = dr;
	closestTower = n;
      }
    } // End of loop over towers

    // Projection factor E --> Et
    // The following is not quite correct, as this factor is different for all towers
    // These values should be in the n-tupel as well
    double projFac   = nJet_->JetEt[j] /  nJet_->JetE[j];

    // Set up jet
    caloJets[j] = new Jet(nJet_->JetPt[j],
			  emf * projFac,
			  had * projFac,
			  out * projFac,
			  nJet_->JetE[j],
			  nJet_->JetEta[j],
			  nJet_->JetPhi[j],
			  0.,
			  Jet::uds,
			  nJet_->GenJetPt[j],
			  drJetGenjet,
			  createCorFactors(j),
			  par_->jet_function(nJet_->TowId_eta[closestTower],
					     nJet_->TowId_phi[closestTower]),
			  jet_error_param,
			  par_->global_jet_function());
    // Read external correction factors
    if(corFactorsFactory_) {
      caloJets[j]->updateCorFactors(corFactorsFactory_->create(caloJets[j]));
    }
    // Correct measurement to L3 (L1*L2*L3)
    if(correctToL3_) {
      caloJets[j]->correctToL3();
    }
    // Correct measurement with L2*L3
    if(correctL2L3_) {
      caloJets[j]->correctL2L3();
    }
  }

  return caloJets;
}



//!  Jet will be sorted by genjet pt
std::vector<Jet*> DiJetReader::readGenJetSortedJets(int nJets) const {
  int maxNJets = nJets;
  if( maxNJets < 0 ) maxNJets = nJet_->NobjJet;

  std::vector<Jet*> jets(maxNJets);
  for(int genJetIdx = 0; genJetIdx < maxNJets; genJetIdx++) {
    int calJetIdx = nJet_->GenJetColJetIdx[genJetIdx]; // Closest (DeltaR) calo jet

    double dphi = nJet_->JetPhi[calJetIdx] - nJet_->GenJetColPhi[genJetIdx];
    if( std::abs(dphi) > 1E-5 && std::abs(dphi) < 1E10 ) {
      dphi = TVector2::Phi_mpi_pi( dphi );
    }
    double deta         = nJet_->JetEta[calJetIdx] - nJet_->GenJetColEta[genJetIdx];
    double drJetGenjet  = sqrt( deta*deta + dphi*dphi );
    double min_tower_dr = 10.;
    double emf          = 0;
    double had          = 0;
    double out          = 0;
    int    closestTower = 0; 
    // Loop over towers
    for (int n=0; n<nJet_->NobjTow; ++n) {
      if (nJet_->Tow_jetidx[n]!=(int)calJetIdx) continue;//look for j-jet's towers
      emf += nJet_->TowEm[n];
      had += nJet_->TowHad[n];
      out += nJet_->TowOE[n];
      dphi = TVector2::Phi_mpi_pi( nJet_->JetPhi[calJetIdx] - nJet_->TowPhi[n] );
      deta = nJet_->JetEta[calJetIdx] - nJet_->TowEta[n];
      double dr = sqrt( deta*deta + dphi*dphi );     
      if (dr < min_tower_dr) {
	min_tower_dr = dr;
	closestTower = n;
      }
    } // End of loop over towers

    // Projection factor E --> Et
    // The following is not quite correct, as this factor is different for all towers
    // These values should be in the n-tupel as well
    double projFac   = nJet_->JetEt[calJetIdx] /  nJet_->JetE[calJetIdx];

    // Set up jet
    jets[genJetIdx] = new Jet(nJet_->JetPt[calJetIdx],
			      emf * projFac,
			      had * projFac,
			      out * projFac,
			      nJet_->JetE[calJetIdx],
			      nJet_->JetEta[calJetIdx],
			      nJet_->JetPhi[calJetIdx],
			      0.,
			      Jet::uds,
			      nJet_->GenJetColPt[genJetIdx],
			      drJetGenjet,
			      createCorFactors(calJetIdx),
			      par_->jet_function(nJet_->TowId_eta[closestTower],
						 nJet_->TowId_phi[closestTower]),
			      jet_error_param,
			      par_->global_jet_function());
    // Read external correction factors
    if(corFactorsFactory_) {
      jets[genJetIdx]->updateCorFactors(corFactorsFactory_->create(jets[genJetIdx]));
    }
    // Correct measurement to L3 (L1*L2*L3)
    if(correctToL3_) {
      jets[genJetIdx]->correctToL3();
    }
    // Correct measurement with L2*L3
    if(correctL2L3_) {
      jets[genJetIdx]->correctL2L3();
    }
  }

  return jets;
}




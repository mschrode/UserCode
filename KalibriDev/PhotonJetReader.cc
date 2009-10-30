//  $Id: PhotonJetReader.cc,v 1.21 2009/10/26 20:56:29 mschrode Exp $

#include "PhotonJetReader.h"

#include "CalibData.h"
#include "SmearPhotonJet.h"
#include "JetTruthEvent.h"
#include "Jet.h"
#include "JetWithTowers.h"
#include "ConfigFile.h"
#include "ToyMC.h"
#include "Parameters.h"

#include <iostream>
#include <cstdlib>

#include <TVector2.h>
#include <TLorentzVector.h>



// ----------------------------------------------------------------   
PhotonJetReader::PhotonJetReader(const std::string& configfile, TParameters* p) :
  EventReader(configfile,p)
{
  // Maximum number of read events
  nGammaJetEvents_     = config_->read<int>("use Gamma-Jet events",-1); 
  if(nGammaJetEvents_ == 0) return;

  // Cuts
  minJetEt_           = config_->read<double>("Et cut on jet",0.0); 
  minGammaEt_         = config_->read<double>("Et cut on gamma",0.0);
  maxRel2ndJetEt_     = config_->read<double>("Relative Rest Jet Cut",0.2);
  minGenJetEt_        = config_->read<double>("Et genJet min",0.0);
  maxGenJetEt_        = config_->read<double>("Et genJet max",10000.0);
  maxDeltaR_          = config_->read<double>("DeltaR cut on jet matching",0.25);
  maxJetEta_          = config_->read<double>("Eta cut on jet",5.0);
  minJetHadFraction_  = config_->read<double>("Min had fraction",0.07);
  maxJetHadFraction_  = config_->read<double>("Max had fraction",0.95);
  // Counter for cutflow
  nMinJetEt_          = 0;
  nMinGammaEt_        = 0;
  nMaxRel2ndJetEt_    = 0;
  nMinGenJetEt_       = 0;    
  nMaxGenJetEt_       = 0;     
  nMaxDeltaR_         = 0;
  nMaxJetEta_         = 0;  
  nMinJetHadFraction_ = 0;     
  nMaxJetHadFraction_ = 0;     

  // Data class
  dataClass_ = config_->read<int>("Gamma-Jet data class", 0);
  if( !( dataClass_ == 0 || dataClass_ == 1 || dataClass_ == 2 || dataClass_ == 5) ) {
    std::cout << "PhotonJetReader: Unknown data class " << dataClass_ << ". Using data class 0." << std::endl;
    dataClass_ = 0;
  }

  // Input files
  string default_tree_name = config_->read<string>("Default Tree Name","CalibTree");
  string treename_gammajet = config_->read<string>("Gamma-Jet tree", default_tree_name);
  TTree* tchain_gammajet;
  vector<string> input_gammajet = 
    bag_of_string(config_->read<string>( "Gamma-Jet input file", "input/gammajet.root" ));
  if(input_gammajet[0] == "toy") {
    std::cout << "generating " << nGammaJetEvents_ << " Gamma-Jet events\n";
    ToyMC* mc = new ToyMC();
    mc->init(configfile);
    mc->print();
    tchain_gammajet = new TTree(treename_gammajet.c_str(),"Gamma Jet events");
    mc->generatePhotonJetTree(tchain_gammajet,nGammaJetEvents_);
    delete mc;
  } else {
    TChain* chain = new TChain(treename_gammajet.c_str());
    for (bag_of_string::const_iterator it = input_gammajet.begin(); it!=input_gammajet.end(); ++it){
      cout << "...opening root-file " << (*it) << " for Gamma-Jet analysis." << endl;
      chain->Add( it->c_str() );
    }
    tchain_gammajet = chain;
  }
  gammaJet_.Init( tchain_gammajet );
}



// ----------------------------------------------------------------   
int PhotonJetReader::readEvents(std::vector<TData*>& data)
{
  if(nGammaJetEvents_ == 0) return 0;

  // Reset counters of rejected events
  nMinJetEt_          = 0;
  nMinGammaEt_        = 0;
  nMaxRel2ndJetEt_    = 0;
  nMinGenJetEt_       = 0;    
  nMaxGenJetEt_       = 0;     
  nMaxDeltaR_         = 0;
  nMaxJetEta_         = 0;  
  nMinJetHadFraction_ = 0;     
  nMaxJetHadFraction_ = 0;     

  int nevent    = gammaJet_.fChain->GetEntries();  // Number of events in chain
  int nReadEvts = 0;                              // Number of read events
  int nGoodEvts = 0;                              // Number of events passing all cuts

  // Some informative output for the interested calibrator
  // Check of correct data class
  cout << "\nGammaJetReader: Reading events of type ";
  if(dataClass_ == 0) {
    std::cout << "'TruthMultMessEvent'";
  } else if((dataClass_ == 1)  || (dataClass_ == 2)) {
    std::cout << "'JetTruthEvent'";
  } else if(dataClass_ == 5) {
    std::cout << "'SmearData'";
  } else {
    std::cerr << "Unknown data class " << dataClass_ << '\n';
    exit(9);
  }
  std::cout << " (data class " << dataClass_ << "):\n";

  for (int i=0;i<nevent;i++) {
    nReadEvts++;

    if((i+1)%10000==0) cout << "  " << i+1 << endl;
    gammaJet_.fChain->GetEvent(i); 
    if (gammaJet_.NobjTowCal>200) {
      cerr<<"ERROR: Increase array sizes in GammaJetSelector; NobjTowCal="
	  <<gammaJet_.NobjTowCal<<"!"<<endl;
      exit(8);
    }
 
    // Trivial cuts
    bool goodEvent = true;
    if( gammaJet_.JetGenEt < minGenJetEt_ ) {
      nMinGenJetEt_++;
      goodEvent = false;
    } else if( gammaJet_.JetGenEt > maxGenJetEt_ ) {
      nMaxGenJetEt_++;
      goodEvent = false;
    } else if( pow(gammaJet_.JetCalEta - gammaJet_.JetGenEta,2)
	       + pow(TVector2::Phi_mpi_pi(gammaJet_.JetCalPhi - gammaJet_.JetGenPhi),2)
	       > pow(maxDeltaR_,2) ) {
      nMaxDeltaR_++;
      goodEvent = false;
    } else if( gammaJet_.PhotonEt < minGammaEt_ ) {
      nMinGammaEt_++;
      goodEvent = false;
    } else if( gammaJet_.JetCalEt < minJetEt_ ) {
      nMinJetEt_++;
      goodEvent = false;
    } else if( gammaJet_.NonLeadingJetPt / gammaJet_.PhotonPt > maxRel2ndJetEt_) {
      nMaxRel2ndJetEt_++;
      goodEvent = false;
    } else if( fabs(gammaJet_.JetCalEta) > maxJetEta_ ) {
      nMaxJetEta_++;
      goodEvent = false;
    }

    if( goodEvent ) {
      
      TData*                                     ev = 0;
      if(dataClass_ == 0)                         ev = createTruthMultMessEvent();
      else if(dataClass_ == 1 || dataClass_ == 2)  ev = createJetTruthEvent();
      else if(dataClass_ == 5)                    ev = createSmearEvent();
      
      if(ev) {
	data.push_back(ev); 
	nGoodEvts++;
      }
    }

    if((nGammaJetEvents_ >= 0) && (nReadEvts >= nGammaJetEvents_)) break;
  }

  // Print cut flow
  std::cout << "Read " << nReadEvts << " gamma-jet events:\n";
  std::cout << "  " << (nReadEvts-=nMinGenJetEt_) << std::flush;
  std::cout << " gamma-jet events with ptgen > " << minGenJetEt_ << "\n";
  std::cout << "  " << (nReadEvts-=nMaxGenJetEt_) << std::flush;
  std::cout << " gamma-jet events with ptgen < " << maxGenJetEt_ << "\n";
  std::cout << "  " << (nReadEvts-=nMaxDeltaR_) << std::flush;
  std::cout << " gamma-jet events with DeltaR < " << maxDeltaR_ << "\n";
  std::cout << "  " << (nReadEvts-=nMinGammaEt_) << std::flush;
  std::cout << " gamma-jet events photon Et > " << minGammaEt_ << "\n";
  std::cout << "  " << (nReadEvts-=nMinJetEt_) << std::flush;
  std::cout << " gamma-jet events jet Et > " << minJetEt_ << "\n";
  std::cout << "  " << (nReadEvts-=nMaxRel2ndJetEt_) << std::flush;
  std::cout << " gamma-jet events (non-leading jet Et) / (photon Et) > " << maxRel2ndJetEt_ << "\n";
  std::cout << "  " << (nReadEvts-=nMaxJetEta_) << std::flush;
  std::cout << " gamma-jet events with |eta| < " << maxJetEta_ << "\n";
  std::cout << "  " << (nReadEvts-=nMinJetHadFraction_) << std::flush;
  std::cout << " gamma-jet events with hadronic fraction > " << minJetHadFraction_ << "\n";
  std::cout << "  " << (nReadEvts-=nMaxJetHadFraction_) << std::flush;
  std::cout << " gamma-jet events with hadronic fraction < " << maxJetHadFraction_ << "\n";
  std::cout << "Stored " << nGoodEvts << " gamma-jet events for analysis.\n";

  return nGoodEvts;
}



// ----------------------------------------------------------------   
TData* PhotonJetReader::createJetTruthEvent()
{
  double em        = 0;
  double had       = 0;
  double out       = 0;
  double dR        = 10;
  int closestTower = 0;
  TMeasurement tower;

  TLorentzVector LJet(0,0,0,0);
  LJet.SetPtEtaPhiE(gammaJet_.JetCalPt,gammaJet_.JetCalEta,gammaJet_.JetCalPhi,gammaJet_.JetCalE);
  TLorentzVector LGenJet(0,0,0,0);
  LGenJet.SetPtEtaPhiE(gammaJet_.JetGenPt,gammaJet_.JetGenEta,gammaJet_.JetGenPhi,gammaJet_.JetGenE);

  // Loop over towers, find closest tower to jet axis,
  // and sum up emf, hadf, outf
  for(int n = 0; n < gammaJet_.NobjTowCal; ++n) {
    em          += gammaJet_.TowEm[n];
    had         +=  gammaJet_.TowHad[n];
    out         +=  gammaJet_.TowOE[n];  
    
    double dphi  = TVector2::Phi_mpi_pi(gammaJet_.JetCalPhi-tower.phi);
    double dr    = sqrt((gammaJet_.JetCalEta-tower.eta)*(gammaJet_.JetCalEta-tower.eta)+
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

  double factor = gammaJet_.JetCalEt /  gammaJet_.JetCalE;

  Jet *j;
  if(dataClass_ == 2) {
    JetWithTowers *jt = 
      new JetWithTowers(gammaJet_.JetCalEt,em * factor,had * factor,
			out * factor,gammaJet_.JetCalE,gammaJet_.JetCalEta,
			gammaJet_.JetCalPhi,TJet::uds,gammaJet_.JetGenEt,LJet.DeltaR(LGenJet),
			TJet::CorFactors(gammaJet_.JetCorrZSP, // L1
					 gammaJet_.JetCorrL2,  // L2
					 gammaJet_.JetCorrL3,  // L3
					 1.,                  // L4
					 1.,                  // L5
					 gammaJet_.JetCorrJPT,
					 gammaJet_.JetCorrL2L3JPT),
			par_->jet_function(gammaJet_.TowId_eta[closestTower],
					gammaJet_.TowId_phi[closestTower]),
			jet_error_param,par_->global_jet_function(),minJetEt_);
    for(int i = 0; i < gammaJet_.NobjTowCal; ++i) {
      double scale = gammaJet_.TowEt[i]/gammaJet_.TowE[i];
      jt->addTower(gammaJet_.TowEt[i],gammaJet_.TowEm[i]*scale,
		   gammaJet_.TowHad[i]*scale,gammaJet_.TowOE[i]*scale,
		   gammaJet_.TowE[i],gammaJet_.TowEta[i],gammaJet_.TowPhi[i],
		   par_->tower_function(gammaJet_.TowId_eta[i],gammaJet_.TowId_phi[i]),
		   tower_error_param);
    }
    j = jt;
  }
  else { 
    j = new Jet(gammaJet_.JetCalEt,em * factor,had * factor,out * factor,
		gammaJet_.JetCalE,gammaJet_.JetCalEta,gammaJet_.JetCalPhi,
		TJet::uds,gammaJet_.JetGenEt,LJet.DeltaR(LGenJet),
		TJet::CorFactors(gammaJet_.JetCorrZSP, // L1
				 gammaJet_.JetCorrL2,  // L2
				 gammaJet_.JetCorrL3,  // L3
				 1.,                  // L4
				 1.,                  // L5
				 gammaJet_.JetCorrJPT,
				 gammaJet_.JetCorrL2L3JPT),
		par_->jet_function(gammaJet_.TowId_eta[closestTower],
				gammaJet_.TowId_phi[closestTower]),
		jet_error_param,par_->global_jet_function(),minJetEt_);
  }
  JetTruthEvent * jte = new JetTruthEvent(j,gammaJet_.PhotonEt,gammaJet_.EventWeight);
  
  return jte;
}



//!  \brief Create SmearPhotonJet event for jet smearing
//!  \note Measured pt is L2L3 corrected
// ----------------------------------------------------------------   
TData* PhotonJetReader::createSmearEvent()
{
  //Find the jets eta & phi index using the nearest tower to jet axis:
  int    jet_index    = -1;
  double min_tower_dr = 10.;
  double em           = 0;
  double had          = 0;
  double out          = 0;
  int    closestTower = 0;

  TLorentzVector LJet(0,0,0,0);
  LJet.SetPtEtaPhiE(gammaJet_.JetCalEt,gammaJet_.JetCalEta,gammaJet_.JetCalPhi,gammaJet_.JetCalE);
  TLorentzVector LGenJet(0,0,0,0);
  LGenJet.SetPtEtaPhiE(gammaJet_.JetGenPt,gammaJet_.JetGenEta,gammaJet_.JetGenPhi,gammaJet_.JetGenE);

  for (int n=0; n<gammaJet_.NobjTowCal; ++n) {
    em  += gammaJet_.TowEm[n];
    had += gammaJet_.TowHad[n];
    out += gammaJet_.TowOE[n];
    TLorentzVector LTower(0,0,0,0);
    LTower.SetPtEtaPhiE(gammaJet_.TowEt[n],gammaJet_.TowEta[n],gammaJet_.TowPhi[n],gammaJet_.TowE[n]);
    double dr = LTower.DeltaR(LJet);
    if (dr<min_tower_dr) {
      min_tower_dr = dr;
      closestTower = n;
    }
  }

  if (jet_index<0){ cerr<<"WARNING: jet_index = " << jet_index << endl; return 0; }
  if(had/(had + em) < minJetHadFraction_) { return 0;}
  if(had/(had + em) > maxJetHadFraction_) { return 0;}

  // Set up measurement
  TJet * jet      = new TJet;
  jet->pt         = gammaJet_.JetCorrL2 * gammaJet_.JetCorrL3 * gammaJet_.JetCalEt;
  jet->eta        = gammaJet_.JetCalEta;
  jet->phi        = gammaJet_.JetCalPhi;
  jet->E          = gammaJet_.JetCalE;
  jet->genPt      = gammaJet_.JetGenPt;
  jet->dR         = LJet.DeltaR(LGenJet);
  jet->corFactors = TJet::CorFactors(gammaJet_.JetCorrZSP, // L1
				     gammaJet_.JetCorrL2,  // L2
				     gammaJet_.JetCorrL3,  // L3
				     1.,                  // L4
				     1.,                  // L5
				     gammaJet_.JetCorrJPT);
  //the following is not quite correct, as this factor is different for all towers. These values should be in the n-tupel as well
  double factor    = gammaJet_.JetCalEt /  gammaJet_.JetCalE;
  jet->HadF       = had * factor;
  jet->EMF        = em * factor;
  jet->OutF       = out * factor;

  SmearPhotonJet * pje = new SmearPhotonJet(jet,gammaJet_.PhotonEt,1.,
					    par_->jet_function(gammaJet_.TowId_eta[closestTower],
							    gammaJet_.TowId_phi[closestTower]));

  return pje;
}



TData* PhotonJetReader::createTruthMultMessEvent() 
{
    //Find the jets eta & phi index using the nearest tower to jet axis:
    int jet_index=-1;
    double min_tower_dr = 10.0;
    double em = 0;
    double had = 0;
    double out = 0;

    TLorentzVector LJet(0,0,0,0);
    LJet.SetPtEtaPhiE(gammaJet_.JetCalEt,gammaJet_.JetCalEta,gammaJet_.JetCalPhi,gammaJet_.JetCalE);
    TLorentzVector LGenJet(0,0,0,0);
    LGenJet.SetPtEtaPhiE(gammaJet_.JetGenPt,gammaJet_.JetGenEta,gammaJet_.JetGenPhi,gammaJet_.JetGenE);

    /*
    for (int n=0; n<gammaJet_.NobjTowCal; ++n){
      TLorentzVector Ltower(0,0,0,0);
      Ltower.SetPtEtaPhiE(gammaJet_.TowEt[n],gammaJet_.TowEta[n],gammaJet_.TowPhi[n],gammaJet_.TowE[n]);
      Ljet += Ltower;
    }
    */
    //Ljet.SetPtEtaPhiE(gammaJet_.JetCalEt,gammaJet_.JetCalEta,gammaJet_.JetCalPhi,gammaJet_.JetCalE);
    for (int n=0; n<gammaJet_.NobjTowCal; ++n){
      em += gammaJet_.TowEm[n];
      had +=  gammaJet_.TowHad[n];
      out +=  gammaJet_.TowOE[n];
      TLorentzVector Ltower(0,0,0,0);
      Ltower.SetPtEtaPhiE(gammaJet_.TowEt[n],gammaJet_.TowEta[n],gammaJet_.TowPhi[n],gammaJet_.TowE[n]);
      double dr = Ltower.DeltaR(LJet);
      if (dr<min_tower_dr) {
	jet_index = par_->GetJetBin(par_->GetJetEtaBin(gammaJet_.TowId_eta[n]),
				 par_->GetJetPhiBin(gammaJet_.TowId_phi[n]));
	min_tower_dr = dr;
      }
    }
    if (jet_index<0){ cerr<<"WARNING: jet_index = " << jet_index << endl; return 0; }
    if(had/(had + em) < minJetHadFraction_) { return 0;}
    if(had/(had + em) > maxJetHadFraction_) { return 0;}
    //jet_index: par_->eta_granularity*par_->phi_granularity*par_->GetNumberOfTowerParametersPerBin()
    //           has to be added for a correct reference to k[...].

    TJet* jetp  = new TJet;
    jetp->pt  = gammaJet_.JetCalEt;
    jetp->eta = gammaJet_.JetCalEta;
    jetp->phi = gammaJet_.JetCalPhi;
    jetp->E   = gammaJet_.JetCalE;
    jetp->genPt =gammaJet_.JetGenPt;
    jetp->dR    = LJet.DeltaR(LGenJet);
    jetp->corFactors = TJet::CorFactors(gammaJet_.JetCorrZSP, // L1
					gammaJet_.JetCorrL2,  // L2
					gammaJet_.JetCorrL3,  // L3
					1.,                  // L4
					1.,                  // L5
					gammaJet_.JetCorrJPT,
					gammaJet_.JetCorrL2L3JPT);
    //the following is not quite correct, as this factor is different for all towers. These values should be in the n-tupel as well
    double factor =  gammaJet_.JetCalEt /  gammaJet_.JetCalE;
    jetp->HadF = had * factor;
    jetp->EMF = em * factor;
    jetp->OutF = out * factor;

    //Create an Gamma/Jet TData event
    TData_TruthMultMess * gj_data = new TData_TruthMultMess
      (
       jet_index  * par_->GetNumberOfJetParametersPerBin() + par_->GetNumberOfTowerParameters(),
       //gammaJet_.PhotonEt,				    //truth//
       gammaJet_.JetGenPt,
       sqrt(pow(0.5,2)+pow(0.10*gammaJet_.PhotonEt,2)),   //error//
       //0.10*gammaJet_.PhotonEt.//error//				    
       gammaJet_.EventWeight,                             //weight//
       //1.0,                                            //weight//
       par_->GetJetParRef( jet_index ),                     //params
       par_->GetNumberOfJetParametersPerBin(),              //number of free jet param. p. bin
       par_->jet_parametrization,                           //function
       jet_error_param,                                  //error param. function
       jetp                                              //measurement
       );

    //double EM=0.,F=0.;
    //Add the jet's towers to "gj_data":
    for (int n=0; n<gammaJet_.NobjTowCal; ++n){
      //if (gammaJet_.TowEt[n]<0.01) continue;
   
      int index = par_->GetBin(par_->GetEtaBin(gammaJet_.TowId_eta[n]),
			    par_->GetPhiBin(gammaJet_.TowId_phi[n]));
      if (index<0){ cerr<<"WARNING: towewer_index = " << index << endl; continue; }

      //double dR = deltaR(gammaJet_.JetCalEta, gammaJet_.JetCalPhi, gammaJet_.TowEta[n], gammaJet_.TowPhi[n]);
	      
      double relativEt = gammaJet_.TowEt[n]/gammaJet_.JetCalEt;  
      //if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
      //This relativeE is used *only* for plotting! Therefore no cuts on this var!
      //create array with multidimensional measurement
      TMeasurement * mess = new TTower;
      mess->pt = double(gammaJet_.TowEt[n]);
      double scale = gammaJet_.TowEt[n]/gammaJet_.TowE[n];
      mess->EMF = double(gammaJet_.TowEm[n]*scale);
      mess->HadF = double(gammaJet_.TowHad[n]*scale);
      mess->OutF = double(gammaJet_.TowOE[n]*scale);
      mess->eta = double(gammaJet_.TowEta[n]);
      mess->phi = double(gammaJet_.TowPhi[n]);
      mess->E = double(gammaJet_.TowE[n]);
      //mess[7] = double( cos( gammaJet_.JetCalPhi-gammaJet_.TowPhi[n] ) ); // Projection factor for summing tower Pt
      //EM+=mess->EMF;
      //F+=mess->pt;
      gj_data->AddMess(new TData_TruthMess(index,
					   mess,                                                    //mess//
					   gammaJet_.PhotonEt * relativEt,                           //truth//
					   //sqrt(1.3 * 1.3/gammaJet_.TowHad[n] + 0.056 * 0.056) * mess->HadF,
					   sqrt(pow(0.5,2)+pow(0.1*gammaJet_.PhotonEt*relativEt,2)), //error//
					   1.,                                                      //weight//
					   par_->GetTowerParRef( index ),                              //parameter//
					   par_->GetNumberOfTowerParametersPerBin(),                   //number of free tower param. p. bin//
					   par_->tower_parametrization,                                //function//
					   tower_error_param                                       //error param.func.//
					   ));
    } 

    //Add the jet's tracks to "gj_data":
    double* EfficiencyMap = par_->GetEffMap();
    int track_index;
    for (int n=0; n<gammaJet_.NobjTrack; ++n){
      if((gammaJet_.TrackTowIdEta[n] == 0) || (gammaJet_.TrackTowIdPhi[n] == 0)) {
	if(gammaJet_.TrackPt[n] > 2){
	  std::cerr << "WARNING: eta or phi id of track is zero!\n";
	  continue;
	}
	else track_index = 0; //bent low momentum tracks with no HCAL hit
      }
      else
	//one trackindex for all tracks in jet = track_index
	track_index = par_->GetTrackBin(par_->GetTrackEtaBin(gammaJet_.TrackTowIdEta[n]),
					 par_->GetTrackPhiBin(gammaJet_.TrackTowIdPhi[n]));
      if (track_index<0){ cerr<<"WARNING: track_index = " << track_index << endl; continue; }
      //create array with multidimensional measurement
      //TMeasurement * Tmess = new TTrack;
      TTrack * Tmess = new TTrack;
      Tmess->TrackId = int(gammaJet_.TrackId[n]);
      Tmess->TowerId = int(gammaJet_.TrackTowId[n]);
      Tmess->pt = double(gammaJet_.TrackPt[n]);
      double scale = gammaJet_.TrackP[n]/gammaJet_.TrackPt[n];
      Tmess->EM1 = double(gammaJet_.TrackEMC1[n]*scale);
      Tmess->EMF = double(gammaJet_.TrackEMC3[n]*scale);
      Tmess->EM5 = double(gammaJet_.TrackEMC5[n]*scale);
      Tmess->Had1 = double(gammaJet_.TrackHAC1[n]*scale);
      Tmess->HadF = double(gammaJet_.TrackHAC3[n]*scale);
      Tmess->Had5 = double(gammaJet_.TrackHAC5[n]*scale);
      Tmess->OutF = 0;
      Tmess->DR = double(gammaJet_.TrackDR[n]);
      Tmess->DRout = double(gammaJet_.TrackDROut[n]);
      Tmess->eta = double(gammaJet_.TrackEta[n]);
      Tmess->etaOut = double(gammaJet_.TrackEtaOut[n]);
      Tmess->phi = double(gammaJet_.TrackPhi[n]);
      Tmess->phiOut = double(gammaJet_.TrackPhiOut[n]);
      Tmess->E = double(gammaJet_.TrackP[n]);
      Tmess->TrackChi2 = double(gammaJet_.TrackChi2[n]);
      Tmess->NValidHits = int(gammaJet_.TrackNHits[n]);
      Tmess->TrackQualityT = bool(gammaJet_.TrackQualityT[n]);
      Tmess->MuDR = double(gammaJet_.MuDR[n]);
      Tmess->MuDE = double(gammaJet_.MuDE[n]);
      int TrackEffBin = par_->GetTrackEffBin(gammaJet_.TrackPt[n],gammaJet_.TrackEta[n]);
      Tmess->Efficiency = EfficiencyMap[TrackEffBin];
      //mess[7] = double( cos( gammaJet_.JetCalPhi-gammaJet_.TowPhi[n] ) ); // Projection factor for summing tower Pt
      //EM+=mess->EMF;
      //F+=mess->pt;
      gj_data->AddTrack(new TData_TruthMess(
					    track_index  * par_->GetNumberOfTrackParametersPerBin() + par_->GetNumberOfTowerParameters() + par_->GetNumberOfJetParameters() ,
					   Tmess,                                                    //mess//
					   0,                           //truth//
					   0.05 + 0.00015 * gammaJet_.TrackPt[n], //error//
					   1.,                                                      //weight//
					   par_->GetTrackParRef( track_index ),                              //parameter//
					   par_->GetNumberOfTrackParametersPerBin(),                   //number of free tower param. p. bin//
					   par_->track_parametrization,                                //function//
					   track_error_param                                        //error param.func.//
					   ));
    }
    gj_data->UseTracks(useTracks_);   //check if track information is sufficient to use Track Parametrization
    return  gj_data;
}

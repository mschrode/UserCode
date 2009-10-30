//
//    Reader for ttbar events
//
//    This class reads events according fo the TopSel
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: TopReader.cc,v 1.16 2009/10/26 20:56:29 mschrode Exp $
//   
#include "TopReader.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"
#include "TwoJetsInvMassEvent.h"
#include "Jet.h"
#include "JetWithTowers.h"

#include "TLorentzVector.h"

#include <iostream>
#include <cstdlib>

TopReader::TopReader(const std::string& configfile, TParameters* p) :
  EventReader(configfile,p),
  minJetEt_     (0.),
  maxJetEta_    (0.),
  minJetHadFrac_(0.07),
  maxJetHadFrac_(0.95),
  useToL3CorrectedJets_(false),
  useMassConstraintW_  (false),
  useMassConstraintTop_(false),
  useGenJetInvMass_(false),
  massConstraintW_  (80.4),
  massConstraintTop_(172.4),
  dataClass_(0),
  createGenWHist_(false)
{
  nTopEvents_     = config_->read<int>("use Top events",-1); 
  if(nTopEvents_ == 0) return;


  minJetEt_  = config_->read<double>("Et cut on jet",0.0);
  maxJetEta_ = config_->read<double>("Eta cut on jet",100.0);
  minJetHadFrac_ = config_->read<double>("Min had fraction", 0.07);
  maxJetHadFrac_ = config_->read<double>("Max had fraction", 0.95);
  useToL3CorrectedJets_ = config_->read<bool>("use to L3 corrected jets",false);
  useMassConstraintW_   = config_->read<bool>("use W mass constraint",true);
  useMassConstraintTop_ = config_->read<bool>("use Top mass constraint",false);
  useGenJetInvMass_     = config_->read<bool>("use GenJet inv mass constraint",false);
  massConstraintW_   = config_->read<double>("W mass",80.4);
  massConstraintTop_ = config_->read<double>("Top mass",172.4);
  if(!useMassConstraintW_ && !useMassConstraintTop_) {
    std::cout << "W mass constraint and Top mass constraint were both turned off," << std::endl
	      << "you should enable at least one of these constraints if you want to run with Top...!" << std::endl;
  }
  dataClass_ = config_->read<int>("Top data class", 0);
  if((dataClass_ != 0)&&(dataClass_ != 1)&&(dataClass_ != 2)) dataClass_ = 0;
  
  string default_tree_name = config_->read<string>("Default Tree Name","CalibTree");
  string treename_top    = config_->read<string>( "Top tree", default_tree_name );
  TChain * tchain_top = new TChain( treename_top.c_str() );
  vector<string> input_top = bag_of_string( config_->read<string>( "Top input file", "input/top.root" ) );
  for (bag_of_string::const_iterator it = input_top.begin(); it!=input_top.end(); ++it){
    cout << "...opening root-file " << (*it) << " for Top analysis." << endl;
    tchain_top->Add( it->c_str() );
  }  
  top_.Init( tchain_top );

  createGenWHist_ = config_->read<bool>("create genW histogram",false);
  if(createGenWHist_)
    genWPtEta_ = new TH2F("genWPtEta", "genWPtEta", 100, -5., 5., 400, 0., 400.);
}

TopReader::~TopReader()
{
  if(createGenWHist_) {
    TFile file("toymcPtEtaInput.root", "recreate");
    genWPtEta_->Write();
    delete genWPtEta_;
  }
}

int TopReader::readEvents(std::vector<TData*>& data)
{
  if(nTopEvents_ == 0) return 0;
  //Run Top stuff  
  int nevent = top_.fChain->GetEntries();
  int evt=0;
  
  for (int i=0;i<nevent;i++) {
    if((i+1)%1000==0) 
      cout<<"Top Event: "<<i+1<<endl;
    top_.fChain->GetEvent(i); 
    if (top_.NobjTow>1000 || top_.NobjJet>8) {
      cerr << "ERROR: Increase array sizes in topSelector; NobjTow="
	   << top_.NobjTow<<", NobjBJet="<<top_.NobjJet<<"!"<<endl;
      exit(10);
    }
    //--------------
    
    if(useMassConstraintW_ && (dataClass_ == 0)) {
      bool goodevent = false;
      TData_InvMass2 * top_data[2]; //two W-jets
      top_data[0] = 0;
      int nstoredjets = 0;
      for (unsigned int ij = 0; ij<3; ++ij){
	if(top_.JetPt[ij] < minJetEt_ || fabs(top_.JetEta[ij]) > maxJetEta_) continue;
	if((TJet::Flavor)top_.JetFlavor[ij] != TJet::uds) continue;
	
	//Find the jets eta & phi index using the nearest tower to jet axis:
	int jet_index=-1;
	double min_tower_dr = 10.0;
	double em = 0;
	double had = 0;
	double out = 0;
	TLorentzVector Ljet(0,0,0,0);
	
	Ljet.SetPtEtaPhiE(top_.JetPt[ij],top_.JetEta[ij],top_.JetPhi[ij],top_.JetE[ij]);
	for (int n=0; n<top_.NobjTow; ++n){
	  if (top_.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	  em += top_.TowEm[n];
	  had += top_.TowHad[n];
	  out += top_.TowOE[n];
	  TLorentzVector Ltower(0,0,0,0);
	  Ltower.SetPtEtaPhiE(top_.TowEt[n],top_.TowEta[n],top_.TowPhi[n],top_.TowE[n]);
	  double dr = Ltower.DeltaR(Ljet);
	  if (dr<min_tower_dr) {
	    jet_index = par_->GetJetBin(par_->GetJetEtaBin(top_.TowId_eta[n]),
				     par_->GetJetPhiBin(top_.TowId_phi[n]));
	    min_tower_dr = dr;
	  }
	}
	if (jet_index<0){
	  cerr<<"WARNING: Top jet_index = " << jet_index << endl;
	  continue;
	}
	
	double * direction = new double[3];
	direction[0] = sin(top_.JetPhi[ij]);
	direction[1] = cos(top_.JetPhi[ij]);
	direction[2] = sqrt(top_.JetE[ij]*top_.JetE[ij] - top_.JetEt[ij]*top_.JetEt[ij]);
	if(top_.JetEta[ij]<0) direction[2] *= -1.;
	TJet* jetp  = new TJet;
	jetp->pt  = top_.JetEt[ij];
	jetp->eta = top_.JetEta[ij];
	jetp->phi = top_.JetPhi[ij];
	jetp->E   = top_.JetE[ij];
	jetp->flavor = (TJet::Flavor)top_.JetFlavor[ij];
	//the following is not quite correct, as this factor is different for all towers.
	//These values should be in the n-tupel as well
	double factor =  top_.JetEt[ij] /  top_.JetE[ij];
	jetp->HadF = had * factor;
	jetp->EMF = em * factor;
	jetp->OutF = out * factor;
	//Create an Top TData event
	top_data[nstoredjets] = 
	  new TData_InvMass2(jet_index * par_->GetNumberOfJetParametersPerBin() + par_->GetNumberOfTowerParameters(),
			     direction,                                     //p_T direction of this jet
			     massConstraintW_,                                         //truth//
			     sqrt(pow(0.5,2)+pow(0.10*top_.JetPt[ij],2)),    //error//
			     top_.Weight,                                    //weight//
			     //1.,                                          //weight//
			     par_->GetJetParRef( jet_index ),                  //params
			     par_->GetNumberOfJetParametersPerBin(),           //number of free jet param. p. bin
			     par_->jet_parametrization,                        //function
			     //par_->dummy_parametrization,
			     jet_error_param,                               //error param. function
			     jetp                                           //jet momentum for plotting and scale
			     );
	
	//Add the jet's towers to "top_data":
	for (int n=0; n<top_.NobjTow; ++n){
	  if (top_.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	  //if (top_.TowEt[n]<0.01) continue;
	  
	  int index = par_->GetBin(par_->GetEtaBin(top_.TowId_eta[n]),
				par_->GetPhiBin(top_.TowId_phi[n]));
	  //std::cout << "jet:" << ij << "bin index:" << index << "\n";
	  if (index<0){ cerr<<"WARNING: Top tower_index = " << index << endl; continue; }
	  
	  double relativEt = top_.TowEt[n]/top_.JetEt[ij];
	  //if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
	  //This relativeE is used *only* for plotting! Therefore no cuts on this var!
	  //create array with multidimensional measurement
	  TMeasurement * mess = new TTower;
	  mess->pt = double(top_.TowEt[n]);
	  double scale = 0.;
	  if (top_.TowE[n]!=0.) scale = top_.TowEt[n]/top_.TowE[n];
	  mess->EMF = double(top_.TowEm[n]*scale);
	  mess->HadF = double(top_.TowHad[n]*scale);
	  mess->OutF = double(top_.TowOE[n]*scale);
	  mess->eta = double(top_.TowEta[n]);
	  mess->phi = double(top_.TowPhi[n]);
	  mess->E = double(top_.TowE[n]);
	  //mess[7] = double( cos( top_.JetCalPhi-top_.TowPhi[n] ) ); // Projection factor for summing tower Pt
	  
	  top_data[nstoredjets]->AddMess(new TData_TruthMess(index,
							     mess,                                                   //mess//
							     top_.JetPt[ij] * relativEt,                             //truth//
							     sqrt(pow(0.5,2)+pow(0.1*top_.JetPt[ij]*relativEt,2)),   //error//
							     //1., //weight//
							     top_.Weight,                                            //weight//
							     par_->GetTowerParRef( index ),                             //parameter//
							     par_->GetNumberOfTowerParametersPerBin(),                  //number of free tower param. p. bin//
							     par_->tower_parametrization,                               //function//
							     tower_error_param                                       //error param. function//
							     ));
	}
	if(nstoredjets> 0) {
          top_data[0]->AddNewMultMess( top_data[nstoredjets] );
	}
	++nstoredjets;
	goodevent = (nstoredjets==2 ? true : false);
      }//loop over all top-jets
      
      if (goodevent) {
	data.push_back( top_data[0] );
      }
    }

    if(useMassConstraintTop_ && (dataClass_ == 0)) {
      bool goodevent = false;
      TData_InvMass2 * top_data[3]; //two W-jets and one b-jet
      top_data[0] = 0;
      int nstoredjets = 0;
      for (unsigned int ij = 0; ij<3; ++ij){
	if(top_.JetPt[ij] < minJetEt_ || fabs(top_.JetEta[ij]) > maxJetEta_) continue;
	
	//Find the jets eta & phi index using the nearest tower to jet axis:
	int jet_index=-1;
	double min_tower_dr = 10.0;
	double em = 0;
	double had = 0;
	double out = 0;
	TLorentzVector Ljet(0,0,0,0);
	
	Ljet.SetPtEtaPhiE(top_.JetPt[ij],top_.JetEta[ij],top_.JetPhi[ij],top_.JetE[ij]);
	for (int n=0; n<top_.NobjTow; ++n){
	  if (top_.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	  em += top_.TowEm[n];
	  had += top_.TowHad[n];
	  out += top_.TowOE[n];
	  TLorentzVector Ltower(0,0,0,0);
	  Ltower.SetPtEtaPhiE(top_.TowEt[n],top_.TowEta[n],top_.TowPhi[n],top_.TowE[n]);
	  double dr = Ltower.DeltaR(Ljet);
	  if (dr<min_tower_dr) {
	    jet_index = par_->GetJetBin(par_->GetJetEtaBin(top_.TowId_eta[n]),
				     par_->GetJetPhiBin(top_.TowId_phi[n]));
	    min_tower_dr = dr;
	  }
	}
	if (jet_index<0){
	  cerr<<"WARNING: Top jet_index = " << jet_index << endl;
	  continue;
	}
	
	double * direction = new double[3];
	direction[0] = sin(top_.JetPhi[ij]);
	direction[1] = cos(top_.JetPhi[ij]);
	direction[2] = sqrt(top_.JetE[ij]*top_.JetE[ij] - top_.JetEt[ij]*top_.JetEt[ij]);
	TJet* jetp  = new TJet;
	jetp->pt  = top_.JetEt[ij];
	jetp->eta = top_.JetEta[ij];
	jetp->phi = top_.JetPhi[ij];
	jetp->E   = top_.JetE[ij];
	jetp->flavor = (TJet::Flavor)top_.JetFlavor[ij];
	//the following is not quite correct, as this factor is different for all towers.
	//These values should be in the n-tupel as well
	double factor =  top_.JetEt[ij] /  top_.JetE[ij];
	jetp->HadF = had * factor;
	jetp->EMF = em * factor;
	jetp->OutF = out * factor;
	//Create an Top TData event
	top_data[nstoredjets] = 
	  new TData_InvMass2(jet_index * par_->GetNumberOfJetParametersPerBin() + par_->GetNumberOfTowerParameters(),
			     direction,                                     //p_T direction of this jet
			     massConstraintTop_,                                         //truth//
			     sqrt(pow(0.5,2)+pow(0.10*top_.JetPt[ij],2)),    //error//
			     top_.Weight,                                    //weight//
			     //1.,                                          //weight//
			     par_->GetJetParRef( jet_index ),                  //params
			     par_->GetNumberOfJetParametersPerBin(),           //number of free jet param. p. bin
			     par_->jet_parametrization,                        //function
			     //par_->dummy_parametrization,
			     jet_error_param,                               //error param. function
			     jetp                                           //jet momentum for plotting and scale
			     );
	
	//Add the jet's towers to "top_data":
	for (int n=0; n<top_.NobjTow; ++n){
	  if (top_.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	  //if (top_.TowEt[n]<0.01) continue;
	  
	  int index = par_->GetBin(par_->GetEtaBin(top_.TowId_eta[n]),
				par_->GetPhiBin(top_.TowId_phi[n]));
	  //std::cout << "jet:" << ij << "bin index:" << index << "\n";
	  if (index<0){ cerr<<"WARNING: Top tower_index = " << index << endl; continue; }
	  
	  double relativEt = top_.TowEt[n]/top_.JetEt[ij];
	  //if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
	  //This relativeE is used *only* for plotting! Therefore no cuts on this var!
	  //create array with multidimensional measurement
	  TMeasurement * mess = new TTower;
	  mess->pt = double(top_.TowEt[n]);
	  double scale = 0.;
	  if (top_.TowE[n]!=0.) scale = top_.TowEt[n]/top_.TowE[n];
	  mess->EMF = double(top_.TowEm[n]*scale);
	  mess->HadF = double(top_.TowHad[n]*scale);
	  mess->OutF = double(top_.TowOE[n]*scale);
	  mess->eta = double(top_.TowEta[n]);
	  mess->phi = double(top_.TowPhi[n]);
	  mess->E = double(top_.TowE[n]);
	  //mess[7] = double( cos( top_.JetCalPhi-top_.TowPhi[n] ) ); // Projection factor for summing tower Pt
	  
	  top_data[nstoredjets]->AddMess(new TData_TruthMess(index,
							     mess,                                                   //mess//
							     top_.JetPt[ij] * relativEt,                             //truth//
							     sqrt(pow(0.5,2)+pow(0.1*top_.JetPt[ij]*relativEt,2)),   //error//
							     //1., //weight//
							     top_.Weight,                                            //weight//
							     par_->GetTowerParRef( index ),                             //parameter//
							     par_->GetNumberOfTowerParametersPerBin(),                  //number of free tower param. p. bin//
							     par_->tower_parametrization,                               //function//
							     tower_error_param                                       //error param. function//
							     ));
	}
	if(nstoredjets> 0) {
          top_data[0]->AddNewMultMess( top_data[nstoredjets] );
	}
	++nstoredjets;
	goodevent = (nstoredjets==3 ? true : false);
      }//loop over all top-jets
      
      if (goodevent) {
	data.push_back( top_data[0] );
      }
    }
    if(dataClass_ > 0) {
      TData* ev = createTwoJetsInvMassEvents();
      if(ev) data.push_back(ev);
      else continue;
    }
    ++evt;
    if (evt>=nTopEvents_ && nTopEvents_>=0)
      break;
  }
  return evt;
}

TData* TopReader::createTwoJetsInvMassEvents()
{
  Jet *jets[2] = {0,0};
  double* terr = new double[top_.NobjTow];
  for(int i = 0; i < 3; ++i) {
    if((TJet::Flavor)top_.JetFlavor[i] != TJet::uds) continue;
    
    double em = 0;
    double had = 0;
    double out = 0;
    double err2 = 0;
    TMeasurement tower;
    double dR = 10;
    int closestTower = 0; 
    for(int n=0; n<top_.NobjTow; ++n){
      if(top_.Tow_jetidx[n] != i) continue;//look for ij-jet's towers
      
      em += top_.TowEm[n];
      had +=  top_.TowHad[n];
      out +=  top_.TowOE[n];  
      tower.pt = top_.TowEt[n];
      double scale = top_.TowEt[n]/top_.TowE[n];
      tower.EMF = top_.TowEm[n]*scale;
      tower.HadF = top_.TowHad[n]*scale;
      tower.OutF = top_.TowOE[n]*scale;
      tower.eta = top_.TowEta[n];
      tower.phi = top_.TowPhi[n];
      tower.E = top_.TowE[n];
      terr[n] = tower_error_param(&tower.pt,&tower,0); 
      if(terr[n] == 0) {
	//assume toy MC???
	terr[n] = TParameters::toy_tower_error_parametrization(&tower.pt,&tower);
      }
      terr[n] *= terr[n];
      err2 += terr[n];
      double dphi = TVector2::Phi_mpi_pi(top_.JetPhi[i]-tower.phi);
      double dr = sqrt((top_.JetEta[i]-tower.eta)*(top_.JetEta[i]-tower.eta)+
		       dphi*dphi);     
      if(dr < dR) {
	dR = dr;
	closestTower = n;
      }
    } 
    if(had/(had + em) < minJetHadFrac_) { delete jets[0]; return 0;}
    if(had/(had + em) > maxJetHadFrac_) { delete jets[0]; return 0;}
    double factor = top_.JetEt[i] / top_.JetE[i];
    tower.pt = top_.JetEt[i];
    tower.EMF = em * factor;
    tower.HadF = had * factor;
    tower.OutF = out * factor;
    tower.eta = top_.JetEta[i];
    tower.phi = top_.JetPhi[i];
    tower.E   = top_.JetE[i];
    double err =  jet_error_param(&tower.pt,&tower,0);
    err2 += err * err;
    Jet **jet = jets[0] ? &jets[1] : &jets[0];
    TJet::CorFactors corFactors = TJet::CorFactors(top_.JetCorrL1[i],
						   top_.JetCorrL2[i],
						   top_.JetCorrL3[i],
						   top_.JetCorrL4[i],            
						   top_.JetCorrL5[i]);
    // use to L3 corrected jets as input if desired
    double jecFactor = 1.;
    if(useToL3CorrectedJets_)
      jecFactor = corFactors.getToL3();
    if(dataClass_ == 2) {
      JetWithTowers *jt = 
	new JetWithTowers(top_.JetEt[i] * jecFactor,
			  em * jecFactor * factor,
			  had * jecFactor * factor,
			  out * jecFactor * factor,
			  top_.JetE[i] * jecFactor,
			  top_.JetEta[i], top_.JetPhi[i],
			  TJet::uds, top_.GenJetPt[i], 0.,
			  corFactors,
			  par_->jet_function(top_.TowId_eta[closestTower],
					  top_.TowId_phi[closestTower]),
			  jet_error_param, par_->global_jet_function(), minJetEt_);
      for(int j = 0 ; j < top_.NobjTow ; ++j) {
	if (top_.Tow_jetidx[j]!= i) continue;//look for ij-jet's towers
	double scale = top_.TowEt[j]/top_.TowE[j];
	jt->addTower(top_.TowEt[j],top_.TowEm[j]*scale,
		     top_.TowHad[j]*scale,top_.TowOE[j]*scale,
		     top_.TowE[j],top_.TowEta[j],top_.TowPhi[j],
		     par_->tower_function(top_.TowId_eta[i],top_.TowId_phi[i]),
		     tower_error_param);
      }
      *jet = jt;
    }
    else { 
      *jet = new Jet(top_.JetEt[i] * jecFactor,
		     em * jecFactor * factor,
		     had * jecFactor * factor,
		     out * jecFactor * factor,
		     top_.JetE[i] * jecFactor,
		     top_.JetEta[i], top_.JetPhi[i],
		     TJet::uds, top_.GenJetPt[i], 0.,
		     corFactors,
		     par_->jet_function(top_.TowId_eta[closestTower],
				     top_.TowId_phi[closestTower]),
		     jet_error_param, par_->global_jet_function(), minJetEt_);    
    }
  }
  delete [] terr;
  if(createGenWHist_) {
    TLorentzVector genW, genJet;
    for(int i = 0; i < 3; ++i) {
      if((TJet::Flavor)top_.JetFlavor[i] != TJet::uds) continue;
      genJet.SetPtEtaPhiE(top_.GenJetPt[i],top_.GenJetEta[i],top_.GenJetPhi[i],top_.GenJetE[i]);
      if(i==0) genW = genJet;
      else genW += genJet;
    }
    genWPtEta_->Fill(genW.Eta(), genW.Pt());
  }
  if(jets[1] &&
     jets[0]->Et() > minJetEt_ && std::abs(jets[0]->eta()) < maxJetEta_ &&
     jets[1]->Et() > minJetEt_ && std::abs(jets[1]->eta()) < maxJetEta_) {
    if(useGenJetInvMass_ ) {
      //compute inv mass of GenJets
      TLorentzVector p1,p2;
      p1.SetPtEtaPhiM(jets[0]->GenPt(),jets[0]->eta(),jets[0]->phi(),0);
      p2.SetPtEtaPhiM(jets[1]->GenPt(),jets[1]->eta(),jets[1]->phi(),0);
      return new TwoJetsInvMassEvent(jets[0],jets[1],(p1+p2).M(),1.0,par_->GetPars());
    }
    return new TwoJetsInvMassEvent(jets[0],jets[1],massConstraintW_,1.0,par_->GetPars());
  } else {
    delete jets[0];
    delete jets[1];
  }
  return 0;
}

//
//    Reader for ttbar events
//
//    This class reads events according fo the TopSel
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: TopReader.cc,v 1.21 2009/11/27 15:28:12 stadie Exp $
//   
#include "TopReader.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"
#include "TwoJetsInvMassEvent.h"
#include "Jet.h"
#include "JetWithTowers.h"
#include "TopSel.h"
#include "CorFactors.h"
#include "CorFactorsFactory.h"

#include "TLorentzVector.h"
#include "TH2F.h"

#include <iostream>
#include <cstdlib>

TopReader::TopReader(const std::string& configfile, TParameters* p) :
  EventReader(configfile,p), top_(new TopSel()),
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
  
  std::string default_tree_name = config_->read<std::string>("Default Tree Name","CalibTree");
  std::string treename_top    = config_->read<std::string>( "Top tree", default_tree_name );
  TChain * tchain_top = new TChain( treename_top.c_str() );
  std::vector<std::string> input_top = bag_of_string( config_->read<std::string>( "Top input file", "input/top.root" ) );
  for (bag_of_string::const_iterator it = input_top.begin(); it!=input_top.end(); ++it){
    std::cout << "...opening root-file " << (*it) << " for Top analysis." 
	      << std::endl;
    tchain_top->Add( it->c_str() );
  }  
  top_->Init( tchain_top );

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

int TopReader::readEvents(std::vector<Event*>& data)
{
  if(nTopEvents_ == 0) return 0;
  //Run Top stuff  
  int nevent = top_->fChain->GetEntries();
  int evt=0;
  
  for (int i=0;i<nevent;i++) {
    if((i+1)%1000==0) 
      std::cout<<"Top Event: "<<i+1<<std::endl;
    top_->fChain->GetEvent(i); 
    if (top_->NobjTow>1000 || top_->NobjJet>8) {
      std::cerr << "ERROR: Increase array sizes in topSelector; NobjTow="
		<< top_->NobjTow<<", NobjBJet="<<top_->NobjJet<<"!"<<std::endl;
      exit(10);
    }

    Event* ev = createTwoJetsInvMassEvents();
    if(ev) data.push_back(ev);
    else continue;
    
    ++evt;
    if (evt>=nTopEvents_ && nTopEvents_>=0)
      break;
  }
  return evt;
}

Event* TopReader::createTwoJetsInvMassEvents()
{
  Jet *jets[2] = {0,0};
  double* terr = new double[top_->NobjTow];
  for(int i = 0; i < 3; ++i) {
    if((Jet::Flavor)top_->JetFlavor[i] != Jet::uds) continue;
    
    double em = 0;
    double had = 0;
    double out = 0;
    double err2 = 0;
    Measurement tower;
    double dR = 10;
    int closestTower = 0; 
    double seta = 0;
    double seta2 = 0;
    double sumpt = 0;
    for(int n=0; n<top_->NobjTow; ++n){
      if(top_->Tow_jetidx[n] != i) continue;//look for ij-jet's towers
      
      em += top_->TowEm[n];
      had +=  top_->TowHad[n];
      out +=  top_->TowOE[n];  
      tower.pt = top_->TowEt[n];
      double scale = top_->TowEt[n]/top_->TowE[n];
      tower.EMF = top_->TowEm[n]*scale;
      tower.HadF = top_->TowHad[n]*scale;
      tower.OutF = top_->TowOE[n]*scale;
      tower.eta = top_->TowEta[n];
      tower.phi = top_->TowPhi[n];
      tower.E = top_->TowE[n];
      terr[n] = tower_error_param(&tower.pt,&tower,0); 
      if(terr[n] == 0) {
	//assume toy MC???
	terr[n] = TParameters::toy_tower_error_parametrization(&tower.pt,&tower);
      }
      terr[n] *= terr[n];
      err2 += terr[n];
      seta += tower.pt  * tower.eta;
      seta2 += tower.pt  * tower.eta * tower.eta;
      sumpt += tower.pt;
      double dphi = TVector2::Phi_mpi_pi(top_->JetPhi[i]-tower.phi);
      double dr = sqrt((top_->JetEta[i]-tower.eta)*(top_->JetEta[i]-tower.eta)+
		       dphi*dphi);     
      if(dr < dR) {
	dR = dr;
	closestTower = n;
      }
    } 
    if(had/(had + em) < minJetHadFrac_) { delete jets[0]; return 0;}
    if(had/(had + em) > maxJetHadFrac_) { delete jets[0]; return 0;}
    double factor = top_->JetEt[i] / top_->JetE[i];
    tower.pt = top_->JetEt[i];
    tower.EMF = em * factor;
    tower.HadF = had * factor;
    tower.OutF = out * factor;
    tower.eta = top_->JetEta[i];
    tower.phi = top_->JetPhi[i];
    tower.E   = top_->JetE[i];
    tower.etaeta = sqrt(seta2/sumpt - seta * seta /(sumpt * sumpt));
    double err =  jet_error_param(&tower.pt,&tower,0);
    err2 += err * err;
    Jet **jet = jets[0] ? &jets[1] : &jets[0];
    CorFactors* corFactors = createCorFactors(i);
    if(dataClass_ == 2) {
      JetWithTowers *jt = 
	new JetWithTowers(top_->JetEt[i],
			  em * factor,
			  had * factor,
			  out * factor,
			  top_->JetE[i],
			  top_->JetEta[i], top_->JetPhi[i],
			  tower.etaeta,
			  Jet::uds, top_->GenJetPt[i], 0.,
			  corFactors,
			  par_->jet_function(top_->TowId_eta[closestTower],
					     top_->TowId_phi[closestTower]),
			  jet_error_param, par_->global_jet_function(), minJetEt_);
      for(int j = 0 ; j < top_->NobjTow ; ++j) {
	if (top_->Tow_jetidx[j]!= i) continue;//look for ij-jet's towers
	double scale = top_->TowEt[j]/top_->TowE[j];
	jt->addTower(top_->TowEt[j],top_->TowEm[j]*scale,
		     top_->TowHad[j]*scale,top_->TowOE[j]*scale,
		     top_->TowE[j],top_->TowEta[j],top_->TowPhi[j],
		     par_->tower_function(top_->TowId_eta[i],top_->TowId_phi[i]),
		     tower_error_param);
      }
      *jet = jt;
    }
    else { 
      *jet = new Jet(top_->JetEt[i],
		     em * factor,
		     had * factor,
		     out * factor,
		     top_->JetE[i],
		     top_->JetEta[i], top_->JetPhi[i],
		     tower.etaeta,
		     Jet::uds, top_->GenJetPt[i], 0.,
		     corFactors,
		     par_->jet_function(top_->TowId_eta[closestTower],
					top_->TowId_phi[closestTower]),
		     jet_error_param, par_->global_jet_function(), minJetEt_);    
    }
    if(corFactorsFactory_) {
      (*jet)->updateCorFactors(corFactorsFactory_->create(*jet));
    }
    if(correctToL3_ || useToL3CorrectedJets_) (*jet)->correctToL3();
  } 
  delete [] terr;
  if(createGenWHist_) {
    TLorentzVector genW, genJet;
    for(int i = 0; i < 3; ++i) {
      if((Jet::Flavor)top_->JetFlavor[i] != Jet::uds) continue;
      genJet.SetPtEtaPhiE(top_->GenJetPt[i],top_->GenJetEta[i],top_->GenJetPhi[i],top_->GenJetE[i]);
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
      p1.SetPtEtaPhiM(jets[0]->genPt(),jets[0]->eta(),jets[0]->phi(),0);
      p2.SetPtEtaPhiM(jets[1]->genPt(),jets[1]->eta(),jets[1]->phi(),0);
      return new TwoJetsInvMassEvent(jets[0],jets[1],(p1+p2).M(),1.0,par_->GetPars());
    }
    return new TwoJetsInvMassEvent(jets[0],jets[1],massConstraintW_,1.0,par_->GetPars());
  } else {
    delete jets[0];
    delete jets[1];
  }
  return 0;
}

CorFactors* TopReader::createCorFactors(int jetid) const
{
  return new CorFactors(top_->JetCorrL1[jetid],
			top_->JetCorrL2[jetid],
			top_->JetCorrL3[jetid],
			top_->JetCorrL4[jetid],            
			top_->JetCorrL5[jetid]);
}

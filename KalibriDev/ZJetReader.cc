//
//    Reader for Z Jet Events
//
//    This class reads events according fo the ZJetSel
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: ZJetReader.cc,v 1.24 2010/04/13 13:44:10 mschrode Exp $
//   
#include "ZJetReader.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"
#include "TLorentzVector.h"
#include "JetTruthEvent.h"
#include "Jet.h"
#include "JetWithTowers.h"
#include "JetWithTracks.h"
#include "ZJetSel.h"
#include "CorFactors.h"
#include "CorFactorsFactory.h"

#include <cstdlib>
#include <iostream>

#include "TLorentzVector.h"


ZJetReader::ZJetReader(const std::string& configfile, TParameters* p) :
  EventReader(configfile,p),zjet(new ZJetSel()),Et_cut_on_Z(0),
  Et_cut_on_jet(0),Had_cut_min(0),Had_cut_max(1)
{
  n_zjet_events     = config_->read<int>("use Z-Jet events",-1); 
  if(n_zjet_events == 0) return ;

  Et_cut_on_Z       = config_->read<double>("Et cut on Z",0.0); 
  Et_cut_on_jet     = config_->read<double>("Et cut on jet",0.0);
  Et_cut_on_genJet  = config_->read<double>("Et genJet min",0.0);
  Eta_cut_on_jet    = config_->read<double>("Eta cut on jet",5.0);
  Had_cut_min       = config_->read<double>("Min had fraction",0.07);
  Had_cut_max       = config_->read<double>("Max had fraction",0.95);

  std::string default_tree_name = config_->read<std::string>("Default Tree Name","CalibTree");
  std::string treename_zjet     = config_->read<std::string>( "Z-Jet tree", default_tree_name );
  TChain * tchain_zjet      = new TChain( treename_zjet.c_str() );
  std::vector<std::string> input_zjet = bag_of_string( config_->read<std::string>( "Z-Jet input file", "input/zjet.root" ) );
  for (bag_of_string::const_iterator it = input_zjet.begin(); it!=input_zjet.end(); ++it){
    std::cout << "...opening root-file " << (*it) << " for Z-Jet analysis." 
	      << std::endl;
    tchain_zjet->Add( it->c_str() );
  }  
  zjet->Init( tchain_zjet );
    
  dataClass = config_->read<int>("Z-Jet data class", 0);
  if((dataClass < 1) || (dataClass > 3)) {
    std::cout << "ZJetReader: Unknown data class " << dataClass << ". Using data class 0." << std::endl;
    dataClass = 0;
  }
}

ZJetReader::~ZJetReader()
{
}

//calculates from Z energy a truth value for one calo tower of the jet.
int ZJetReader::readEvents(std::vector<Event*>& data)
{
  if(n_zjet_events == 0) return 0;
  //Run Z-Jet stuff  
  int nevent = zjet->fChain->GetEntries();
  int nevents_added = 0;
  for (int i=0;i<nevent;i++) {
    if(i%1000==0) std::cout<<"Z-Jet Event: "<<i<<std::endl;
    zjet->fChain->GetEvent(i); 
    if (zjet->NobjTowCal>200) {
      std::cerr<<"ERROR: Increase array sizes in ZJetSelector; NobjTowCal="
	       <<zjet->NobjTowCal<<"!"<<std::endl;
      exit(8);
    }
    if (zjet->NobjETowCal>200) {
      std::cerr<<"ERROR: Increase array sizes in ZJetSelector; NobjETowCal="
	       <<zjet->NobjETowCal<<"!"<<std::endl;
      exit(8);
    }   
    if (zjet->NobjTrack>200) {
      std::cerr<<"ERROR: Increase array sizes in ZJetSelector; NobjTrackCal="
	       <<zjet->NobjTrack<<"!"<<std::endl;
      exit(8);
    }
    //if (zjet->ZPt<Et_cut_on_Z || zjet->JetCalPt<Et_cut_on_jet) continue;
    //if(zjet->JetGenEt < Et_cut_on_genJet  Et_cut_on_Z || zjet->JetCalPt<Et_cut_on_jet) continue;
    if(zjet->JetGenEt < Et_cut_on_genJet || zjet->JetCalPt<Et_cut_on_jet) continue;
    //temporary solution to get rid of wrong Zs as long as muon cleaning is not applied
    if(std::abs(zjet->ZPt - zjet->GenZPt) > 80) continue;    
    //cut on Eta to avoid scarcely populated bins
    //if(std::abs(zjet->JetCalEta) > 3.0) continue;

    Event* ev = createJetTruthEvent();
     
    if(ev) {
      data.push_back(ev); 
      ++nevents_added;
      if((n_zjet_events>=0) && (nevents_added >= n_zjet_events))
	break;
    }
  }
  return nevents_added;
}

Event* ZJetReader::createJetTruthEvent()
{
  double em = 0;
  double had = 0;
  double out = 0;
  double err2 = 0;
  Measurement tower;
  double* terr = new double[zjet->NobjTowCal];
  double dR = 10;
  int closestTower = 0;
  double seta = 0;
  double seta2 = 0; 
  double sphi = 0;
  double sphi2 = 0;
  double sumpt = 0;
  for(int n = 0; n < zjet->NobjTowCal; ++n) {
    em += zjet->TowEm[n];
    had +=  zjet->TowHad[n];
    out +=  zjet->TowOE[n];  
    tower.pt = zjet->TowEt[n];
    double scale = zjet->TowEt[n]/zjet->TowE[n];
    tower.EMF = zjet->TowEm[n]*scale;
    tower.HadF = zjet->TowHad[n]*scale;
    tower.OutF = zjet->TowOE[n]*scale;
    tower.eta = zjet->TowEta[n];
    tower.phi = zjet->TowPhi[n];
    tower.E = zjet->TowE[n];
    terr[n] = tower_error_param(&tower.pt,&tower,0); 
    if(terr[n] == 0) {
      //assume toy MC???
      terr[n] = TParameters::toy_tower_error_parametrization(&tower.pt,&tower);
    }
    terr[n] *= terr[n];
    err2 += terr[n];  
    seta += tower.pt  * tower.eta;
    seta2 += tower.pt * tower.eta * tower.eta;
    sphi += tower.pt  * tower.phi;
    sphi2 += tower.pt * tower.phi * tower.phi;
    sumpt += tower.pt;
    double dphi = TVector2::Phi_mpi_pi(zjet->JetCalPhi-tower.phi);
    double dr = sqrt((zjet->JetCalEta-tower.eta)*(zjet->JetCalEta-tower.eta)+
		     dphi*dphi);     
    if(dr < dR) {
      dR = dr;
      closestTower = n;
    }
  }  //calc jet error
  if(had/(had + em) < Had_cut_min) { return 0;}
  if(had/(had + em) > Had_cut_max) { return 0;}
  double factor =  zjet->JetCalEt /  zjet->JetCalE;
  tower.pt = zjet->JetCalEt;
  tower.EMF = em * factor;
  tower.HadF = had * factor;
  tower.OutF = out * factor;
  tower.eta = zjet->JetCalEta;
  tower.phi = zjet->JetCalPhi;
  tower.E   = zjet->JetCalE;
  tower.etaeta = sqrt(seta2/sumpt - seta * seta /(sumpt * sumpt));
  tower.phiphi = sqrt(sphi2/sumpt - sphi * sphi /(sumpt * sumpt));
  double err =  jet_error_param(&tower.pt,&tower,0);
  err2 += err * err;

  TLorentzVector LJet(0,0,0,0);
  LJet.SetPtEtaPhiE(zjet->JetCalPt,zjet->JetCalEta,zjet->JetCalPhi,zjet->JetCalE);
  TLorentzVector LGenJet(0,0,0,0);
  LGenJet.SetPtEtaPhiE(zjet->JetGenPt,zjet->JetGenEta,zjet->JetCalPhi,zjet->JetGenE);

  Jet *j;
  if(dataClass == 2) {
    JetWithTowers *jt = 
      new JetWithTowers(zjet->JetCalEt,em * factor,had * factor,
			out * factor,zjet->JetCalE,zjet->JetCalEta,
			zjet->JetCalPhi,tower.phiphi,tower.etaeta,Jet::uds,
			zjet->JetGenEt,LJet.DeltaR(LGenJet),
			createCorFactors(0),
			par_->jet_function(zjet->TowId_eta[closestTower],
					   zjet->TowId_phi[closestTower]),
			jet_error_param,par_->global_jet_function());
    for(int i = 0; i < zjet->NobjTowCal; ++i) {
      double scale = zjet->TowEt[i]/zjet->TowE[i];
      jt->addTower(zjet->TowEt[i],zjet->TowEm[i]*scale,
		   zjet->TowHad[i]*scale,zjet->TowOE[i]*scale,
		   zjet->TowE[i],zjet->TowEta[i],zjet->TowPhi[i], 
		   par_->tower_function(zjet->TowId_eta[i],zjet->TowId_phi[i]),
		   tower_error_param);
    }
    j = jt;
  } else  if(dataClass == 3) {
    JetWithTracks *jt = 
      new JetWithTracks(zjet->JetCalEt,em * factor,had * factor,
			out * factor, zjet->JetCalE,zjet->JetCalEta,
			zjet->JetCalPhi,tower.phiphi,tower.etaeta,Jet::uds,
			zjet->JetGenEt,LJet.DeltaR(LGenJet),createCorFactors(0),
			par_->jet_function(zjet->TowId_eta[closestTower],
					   zjet->TowId_phi[closestTower]),
			jet_error_param,par_->global_jet_function());
    double* EfficiencyMap = par_->GetEffMap();
    for(int i = 0; i < zjet->NobjTrack; ++i) {
      double scale = zjet->TrackP[i]/zjet->TrackPt[i];  
      int TrackEffBin = par_->GetTrackEffBin(zjet->TrackPt[i],zjet->TrackEta[i]);
      double eff = EfficiencyMap[TrackEffBin];
      jt->addTrack(zjet->TrackP[i],zjet->TrackEMC3[i]*scale,zjet->TrackHAC3[i]*scale,0,zjet->TrackP[i],
		   zjet->TrackEta[i],zjet->TrackPhi[i],zjet->TrackId[i],zjet->TrackTowId[i],
		   zjet->TrackDR[i],zjet->TrackDROut[i],zjet->TrackEtaOut[i],zjet->TrackPhiOut[i],
		   zjet->TrackEMC1[i]*scale,zjet->TrackEMC5[i]*scale,zjet->TrackHAC1[i],zjet->TrackHAC5[i],
		   zjet->TrackChi2[i],zjet->TrackNHits[i],zjet->TrackQualityT[i],zjet->MuDR[i],zjet->MuDE[i],
		   eff,par_->track_function(zjet->TrackTowIdEta[i],zjet->TrackTowIdPhi[i]),track_error_param);
    }
    j = jt;
  } else { 
    j = new Jet(zjet->JetCalEt,em * factor,had * factor,out * factor,
		zjet->JetCalE,zjet->JetCalEta,zjet->JetCalPhi,tower.etaeta,
		tower.phiphi,Jet::uds,zjet->JetGenEt,LJet.DeltaR(LGenJet),
		createCorFactors(0),
		par_->jet_function(zjet->TowId_eta[closestTower],
				   zjet->TowId_phi[closestTower]),
		jet_error_param,par_->global_jet_function());
  } 
  if(corFactorsFactory_) {
      j->updateCorFactors(corFactorsFactory_->create(j));
  }
  if(correctToL3_) j->correctToL3();
  else if(correctL2L3_) j->correctL2L3();
  JetTruthEvent* jte = new JetTruthEvent(j,zjet->JetGenEt,1.0);//zjet->EventWeight);
  delete [] terr;
  return jte;
}

CorFactors* ZJetReader::createCorFactors(int jetid) const
{
  return new CorFactors(zjet->JetCorrZSP, // L1
			zjet->JetCorrL2,  // L2
			zjet->JetCorrL3,  // L3
			1.,              // L4
			1.,              // L5
			zjet->JetCorrJPT,
			zjet->JetCorrL2L3JPT);
}

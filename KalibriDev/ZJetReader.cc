//
//    Reader for Z Jet Events
//
//    This class reads events according fo the ZJetSel
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: ZJetReader.cc,v 1.17 2009/10/26 20:56:30 mschrode Exp $
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

#include <cstdlib>
#include <iostream>

#include "TLorentzVector.h"


ZJetReader::ZJetReader(const std::string& configfile, TParameters* p) :
 EventReader(configfile,p),Et_cut_on_Z(0),Et_cut_on_jet(0),Had_cut_min(0),Had_cut_max(1)
{
  n_zjet_events     = config_->read<int>("use Z-Jet events",-1); 
  if(n_zjet_events == 0) return ;

  Et_cut_on_Z       = config_->read<double>("Et cut on Z",0.0); 
  Et_cut_on_jet     = config_->read<double>("Et cut on jet",0.0);
  Et_cut_on_genJet  = config_->read<double>("Et genJet min",0.0);
  Eta_cut_on_jet    = config_->read<double>("Eta cut on jet",5.0);
  Had_cut_min       = config_->read<double>("Min had fraction",0.07);
  Had_cut_max       = config_->read<double>("Max had fraction",0.95);

  string default_tree_name = config_->read<string>("Default Tree Name","CalibTree");
  string treename_zjet      = config_->read<string>( "Z-Jet tree", default_tree_name );
  TChain * tchain_zjet      = new TChain( treename_zjet.c_str() );
  vector<string> input_zjet = bag_of_string( config_->read<string>( "Z-Jet input file", "input/zjet.root" ) );
  for (bag_of_string::const_iterator it = input_zjet.begin(); it!=input_zjet.end(); ++it){
    cout << "...opening root-file " << (*it) << " for Z-Jet analysis." << endl;
    tchain_zjet->Add( it->c_str() );
  }  
  zjet.Init( tchain_zjet );
    
  dataClass = config_->read<int>("Z-Jet data class", 0);
  if((dataClass < 0) || (dataClass > 3)) {
    std::cout << "ZJetReader: Unknown data class " << dataClass << ". Using data class 0." << std::endl;
    dataClass = 0;
  }
}

ZJetReader::~ZJetReader()
{
}

//calculates from Z energy a truth value for one calo tower of the jet.
int ZJetReader::readEvents(std::vector<TData*>& data)
{
  if(n_zjet_events == 0) return 0;
  //Run Z-Jet stuff  
  int nevent = zjet.fChain->GetEntries();
  int nevents_added = 0;
  for (int i=0;i<nevent;i++) {
    if(i%1000==0) cout<<"Z-Jet Event: "<<i<<endl;
    zjet.fChain->GetEvent(i); 
    if (zjet.NobjTowCal>200) {
      cerr<<"ERROR: Increase array sizes in ZJetSelector; NobjTowCal="
	  <<zjet.NobjTowCal<<"!"<<endl;
      exit(8);
    }
    if (zjet.NobjETowCal>200) {
      cerr<<"ERROR: Increase array sizes in ZJetSelector; NobjETowCal="
	  <<zjet.NobjETowCal<<"!"<<endl;
      exit(8);
    }   
    if (zjet.NobjTrack>200) {
      cerr<<"ERROR: Increase array sizes in ZJetSelector; NobjTrackCal="
	  <<zjet.NobjTrack<<"!"<<endl;
      exit(8);
    }
    //if (zjet.ZPt<Et_cut_on_Z || zjet.JetCalPt<Et_cut_on_jet) continue;
    //if(zjet.JetGenEt < Et_cut_on_genJet  Et_cut_on_Z || zjet.JetCalPt<Et_cut_on_jet) continue;
    if(zjet.JetGenEt < Et_cut_on_genJet || zjet.JetCalPt<Et_cut_on_jet) continue;
    //temporary solution to get rid of wrong Zs as long as muon cleaning is not applied
    if(std::abs(zjet.ZPt - zjet.GenZPt) > 80) continue;    
    //cut on Eta to avoid scarcely populated bins
    //if(std::abs(zjet.JetCalEta) > 3.0) continue;

    TData* ev = 0;
    if(dataClass == 0) ev = createTruthMultMessEvent();
    else if(dataClass > 0) ev = createJetTruthEvent();
     
    if(ev) {
      data.push_back(ev); 
      ++nevents_added;
      if((n_zjet_events>=0) && (nevents_added >= n_zjet_events))
	break;
    }
  }
  return nevents_added;
}

TData* ZJetReader::createJetTruthEvent()
{
  double em = 0;
  double had = 0;
  double out = 0;
  double err2 = 0;
  TMeasurement tower;
  double* terr = new double[zjet.NobjTowCal];
  double dR = 10;
  int closestTower = 0;
  for(int n = 0; n < zjet.NobjTowCal; ++n) {
    em += zjet.TowEm[n];
    had +=  zjet.TowHad[n];
    out +=  zjet.TowOE[n];  
    tower.pt = zjet.TowEt[n];
    double scale = zjet.TowEt[n]/zjet.TowE[n];
    tower.EMF = zjet.TowEm[n]*scale;
    tower.HadF = zjet.TowHad[n]*scale;
    tower.OutF = zjet.TowOE[n]*scale;
    tower.eta = zjet.TowEta[n];
    tower.phi = zjet.TowPhi[n];
    tower.E = zjet.TowE[n];
    terr[n] = tower_error_param(&tower.pt,&tower,0); 
    if(terr[n] == 0) {
      //assume toy MC???
      terr[n] = TParameters::toy_tower_error_parametrization(&tower.pt,&tower);
    }
    terr[n] *= terr[n];
    err2 += terr[n];
    double dphi = TVector2::Phi_mpi_pi(zjet.JetCalPhi-tower.phi);
    double dr = sqrt((zjet.JetCalEta-tower.eta)*(zjet.JetCalEta-tower.eta)+
		     dphi*dphi);     
    if(dr < dR) {
      dR = dr;
      closestTower = n;
    }
  }  //calc jet error
  if(had/(had + em) < Had_cut_min) { return 0;}
  if(had/(had + em) > Had_cut_max) { return 0;}
  double factor =  zjet.JetCalEt /  zjet.JetCalE;
  tower.pt = zjet.JetCalEt;
  tower.EMF = em * factor;
  tower.HadF = had * factor;
  tower.OutF = out * factor;
  tower.eta = zjet.JetCalEta;
  tower.phi = zjet.JetCalPhi;
  tower.E   = zjet.JetCalE;
  double err =  jet_error_param(&tower.pt,&tower,0);
  err2 += err * err;

  TLorentzVector LJet(0,0,0,0);
  LJet.SetPtEtaPhiE(zjet.JetCalPt,zjet.JetCalEta,zjet.JetCalPhi,zjet.JetCalE);
  TLorentzVector LGenJet(0,0,0,0);
  LGenJet.SetPtEtaPhiE(zjet.JetGenPt,zjet.JetGenEta,zjet.JetCalPhi,zjet.JetGenE);

  Jet *j;
  if(dataClass == 2) {
    JetWithTowers *jt = 
      new JetWithTowers(zjet.JetCalEt,em * factor,had * factor,
			out * factor,zjet.JetCalE,zjet.JetCalEta,
			zjet.JetCalPhi,TJet::uds,zjet.JetGenEt,LJet.DeltaR(LGenJet),
			TJet::CorFactors(zjet.JetCorrZSP, // L1
					 zjet.JetCorrL2,  // L2
					 zjet.JetCorrL3,  // L3
					 1.,              // L4
					 1.,              // L5
					 zjet.JetCorrJPT,
					 zjet.JetCorrL2L3JPT),
			par_->jet_function(zjet.TowId_eta[closestTower],
					zjet.TowId_phi[closestTower]),
			jet_error_param,par_->global_jet_function());
    for(int i = 0; i < zjet.NobjTowCal; ++i) {
      double scale = zjet.TowEt[i]/zjet.TowE[i];
      jt->addTower(zjet.TowEt[i],zjet.TowEm[i]*scale,
		   zjet.TowHad[i]*scale,zjet.TowOE[i]*scale,
		   zjet.TowE[i],zjet.TowEta[i],zjet.TowPhi[i], 
		   par_->tower_function(zjet.TowId_eta[i],zjet.TowId_phi[i]),
		   tower_error_param);
    }
    j = jt;
  } else  if(dataClass == 3) {
    JetWithTracks *jt = 
      new JetWithTracks(zjet.JetCalEt,em * factor,had * factor,
			out * factor, zjet.JetCalE,zjet.JetCalEta,
			zjet.JetCalPhi,TJet::uds,zjet.JetGenEt,LJet.DeltaR(LGenJet),
			TJet::CorFactors(zjet.JetCorrZSP, // L1
					 zjet.JetCorrL2,  // L2
					 zjet.JetCorrL3,  // L3
					 1.,              // L4
					 1.,              // L5
					 zjet.JetCorrJPT,
					 zjet.JetCorrL2L3JPT),
			par_->jet_function(zjet.TowId_eta[closestTower],
					zjet.TowId_phi[closestTower]),
			jet_error_param,par_->global_jet_function());
    double* EfficiencyMap = par_->GetEffMap();
    for(int i = 0; i < zjet.NobjTrack; ++i) {
      double scale = zjet.TrackP[i]/zjet.TrackPt[i];  
      int TrackEffBin = par_->GetTrackEffBin(zjet.TrackPt[i],zjet.TrackEta[i]);
      double eff = EfficiencyMap[TrackEffBin];
      jt->addTrack(zjet.TrackP[i],zjet.TrackEMC3[i]*scale,zjet.TrackHAC3[i]*scale,0,zjet.TrackP[i],
		   zjet.TrackEta[i],zjet.TrackPhi[i],zjet.TrackId[i],zjet.TrackTowId[i],
		   zjet.TrackDR[i],zjet.TrackDROut[i],zjet.TrackEtaOut[i],zjet.TrackPhiOut[i],
		   zjet.TrackEMC1[i]*scale,zjet.TrackEMC5[i]*scale,zjet.TrackHAC1[i],zjet.TrackHAC5[i],
		   zjet.TrackChi2[i],zjet.TrackNHits[i],zjet.TrackQualityT[i],zjet.MuDR[i],zjet.MuDE[i],
		   eff,par_->track_function(zjet.TrackTowIdEta[i],zjet.TrackTowIdPhi[i]),track_error_param);
    }
    j = jt;
  } else { 
    j = new Jet(zjet.JetCalEt,em * factor,had * factor,out * factor,
		zjet.JetCalE,zjet.JetCalEta,zjet.JetCalPhi,TJet::uds,
		zjet.JetGenEt,LJet.DeltaR(LGenJet),
		TJet::CorFactors(zjet.JetCorrZSP, // L1
				 zjet.JetCorrL2,  // L2
				 zjet.JetCorrL3,  // L3
				 1.,              // L4
				 1.,              // L5
				 zjet.JetCorrJPT,
				 zjet.JetCorrL2L3JPT),
		par_->jet_function(zjet.TowId_eta[closestTower],
				zjet.TowId_phi[closestTower]),
		jet_error_param,par_->global_jet_function());
  }
  JetTruthEvent* jte = new JetTruthEvent(j,zjet.JetGenEt,1.0);//zjet.EventWeight);
  delete [] terr;
  return jte;
}

TData* ZJetReader::createTruthMultMessEvent()
{
  //Find the jets eta & phi index using the nearest tower to jet axis:
    int jet_index=-1;
    double min_tower_dr = 10.0;
    double em = 0;
    double had = 0;
    double out = 0;
    TLorentzVector Ljet(0,0,0,0);
    Ljet.SetPtEtaPhiE(zjet.JetCalEt,zjet.JetCalEta,zjet.JetCalPhi,zjet.JetCalE);
    for (int n=0; n<zjet.NobjTowCal; ++n){
      em += zjet.TowEm[n];
      had +=  zjet.TowHad[n];
      out +=  zjet.TowOE[n];
      TLorentzVector Ltower(0,0,0,0);
      Ltower.SetPtEtaPhiE(zjet.TowEt[n],zjet.TowEta[n],zjet.TowPhi[n],zjet.TowE[n]);
      double dr = Ltower.DeltaR(Ljet);
      if (dr<min_tower_dr) {
	jet_index = par_->GetJetBin(par_->GetJetEtaBin(zjet.TowId_eta[n]),
				 par_->GetJetPhiBin(zjet.TowId_phi[n]));
	min_tower_dr = dr;
      }
    }

    if (jet_index<0){ cerr<<"WARNING: jet_index = " << jet_index << endl; return 0; }
    if(em == 0) { return 0;}
    //jet_index: par_->eta_granularity*par_->phi_granularity*par_->GetNumberOfTowerParametersPerBin()
    //           has to be added for a correct reference to k[...].

    TJet* jetp  = new TJet;
    jetp->pt  = zjet.JetCalEt;
    jetp->eta = zjet.JetCalEta;
    jetp->phi = zjet.JetCalPhi;
    jetp->E   = zjet.JetCalE;
    jetp->genPt =zjet.JetGenPt;
    jetp->corFactors = TJet::CorFactors(zjet.JetCorrZSP, // L1
					zjet.JetCorrL2,  // L2
					zjet.JetCorrL3,  // L3
					1.,              // L4
					1.,              // L5
					zjet.JetCorrJPT,
					zjet.JetCorrL2L3JPT);
    //the following is not quite correct, as this factor is slightly different for all towers.
    double factor =  zjet.JetCalEt / zjet.JetCalE;
    jetp->HadF = had * factor;
    jetp->EMF = em * factor;
    jetp->OutF = out * factor;

    //Create an Z/Jet TData event
    TData_TruthMultMess * gj_data = new 
      TData_TruthMultMess(jet_index  * par_->GetNumberOfJetParametersPerBin() + par_->GetNumberOfTowerParameters(),
			  //zjet.ZPt,				    //truth//
			  zjet.JetGenPt,
			  sqrt(pow(0.5,2)+pow(0.10*zjet.ZEt,2)),    //error//
			  //zjet.EventWeight,                       //weight//
			  1.0,                                      //weight//
			  par_->GetJetParRef( jet_index ),             //params
			  par_->GetNumberOfJetParametersPerBin(),      //number of free jet param. p. bin
			  par_->jet_parametrization,                   //function
			  jet_error_param,                          //error param. function
			  jetp
			  );
    
    //Add the jet's towers to "gj_data":
    for (int n=0; n<zjet.NobjTowCal; ++n){
      //if (zjet.TowEt[n]<0.01) continue;
      
      int index = par_->GetBin(par_->GetEtaBin(zjet.TowId_eta[n]),
			    par_->GetPhiBin(zjet.TowId_phi[n]));
      if (index<0){ cerr<<"WARNING: towewer_index = " << index << endl; continue; }
      
      //double dR = deltaR(zjet.JetCalEta, zjet.JetCalPhi, zjet.TowEta[n], zjet.TowPhi[n]);
      
      double relativEt = zjet.TowEt[n]/zjet.JetCalEt;  
      //if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
      //This relativeE is used *only* for plotting! Therefore no cuts on this var!
      //create array with multidimensional measurement
      TMeasurement * mess = new TTower;
      mess->pt = double(zjet.TowEt[n]);
      double scale = zjet.TowEt[n]/zjet.TowE[n];
      mess->EMF = double(zjet.TowEm[n]*scale);
      mess->HadF = double(zjet.TowHad[n]*scale);
      mess->OutF = double(zjet.TowOE[n]*scale);
      mess->eta = double(zjet.TowEta[n]);
      mess->phi = double(zjet.TowPhi[n]);
      mess->E = double(zjet.TowE[n]);
      //mess[7] = double( cos( zjet.JetCalPhi-zjet.TowPhi[n] ) ); // Projection factor for summing tower Pt
      gj_data->AddMess(new TData_TruthMess(index,
					   mess,                                           //mess//
					   zjet.ZEt * relativEt,                           //truth//
					   sqrt(pow(0.5,2)+pow(0.1*zjet.ZEt*relativEt,2)), //error//
					   1.,                                             //weight//
					   par_->GetTowerParRef( index ),                     //parameter//
					   par_->GetNumberOfTowerParametersPerBin(),          //number of free tower param. p. bin//
					   par_->tower_parametrization,                       //function//
					   tower_error_param                                        //error param.func.//
					   ));
    } 

    
    //Add the jet's tracks to "gj_data":
    int index;
    double* EfficiencyMap = par_->GetEffMap();
    for (int n=0; n<zjet.NobjTrack; ++n){
      if((zjet.TrackTowIdEta[n] == 0) || (zjet.TrackTowIdPhi[n] == 0)) {
	if(zjet.TrackPt[n] > 2){
	  std::cerr << "WARNING: eta or phi id of track is zero!\n";
	  continue;
	}
	else index = 0; //bent low momentum tracks with no HCAL hit
      }
      else
	index = par_->GetTrackBin(par_->GetTrackEtaBin(zjet.TrackTowIdEta[n]),
			       par_->GetTrackPhiBin(zjet.TrackTowIdPhi[n]));
      if (index<0){ cerr<<"WARNING: track_index = " << index << endl; continue; }
      //create array with multidimensional measurement
      //TMeasurement * Tmess = new TTrack;
      TTrack * Tmess = new TTrack;
      Tmess->TrackId = int(zjet.TrackId[n]);
      Tmess->TowerId = int(zjet.TrackTowId[n]);
      Tmess->pt = double(zjet.TrackPt[n]);
      double scale = zjet.TrackP[n]/zjet.TrackPt[n];
      Tmess->EM1 = double(zjet.TrackEMC1[n]*scale);
      Tmess->EMF = double(zjet.TrackEMC3[n]*scale);
      Tmess->EM5 = double(zjet.TrackEMC5[n]*scale);
      Tmess->Had1 = double(zjet.TrackHAC1[n]*scale);
      Tmess->HadF = double(zjet.TrackHAC3[n]*scale);
      Tmess->Had5 = double(zjet.TrackHAC5[n]*scale);
      Tmess->OutF = 0;
      Tmess->DR = double(zjet.TrackDR[n]);
      Tmess->DRout = double(zjet.TrackDROut[n]);
      Tmess->eta = double(zjet.TrackEta[n]);
      Tmess->etaOut = double(zjet.TrackEtaOut[n]);
      Tmess->phi = double(zjet.TrackPhi[n]);
      Tmess->phiOut = double(zjet.TrackPhiOut[n]);
      Tmess->E = double(zjet.TrackP[n]);
      Tmess->TrackChi2 = double(zjet.TrackChi2[n]);
      Tmess->NValidHits = int(zjet.TrackNHits[n]);
      Tmess->TrackQualityT = bool(zjet.TrackQualityT[n]);
      Tmess->MuDR = double(zjet.MuDR[n]);
      Tmess->MuDE = double(zjet.MuDE[n]);
      int TrackEffBin = par_->GetTrackEffBin(zjet.TrackPt[n],zjet.TrackEta[n]);
      Tmess->Efficiency = EfficiencyMap[TrackEffBin];
      //mess[7] = double( cos( zjet.JetCalPhi-zjet.TowPhi[n] ) ); // Projection factor for summing tower Pt
      //EM+=mess->EMF;
      //F+=mess->pt;
       gj_data->AddTrack(new TData_TruthMess(index  * par_->GetNumberOfTrackParametersPerBin() + par_->GetNumberOfTowerParameters() + par_->GetNumberOfJetParameters() ,
					   Tmess,                                                    //mess//
					   0,                           //truth//
					   0.015 * zjet.TrackPt[n], //error//
					   1.,                                                      //weight//
					   par_->GetTrackParRef( index ),                              //parameter//
					   par_->GetNumberOfTrackParametersPerBin(),                   //number of free tower param. p. bin//
					   par_->track_parametrization,                                //function//
					   track_error_param                                        //error param.func.//
					   ));
    } 
    gj_data->UseTracks(useTracks_);   //check if track information is sufficient to use Track Parametrization
    return gj_data;
}

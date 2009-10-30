//
//    Reader for Track Cluster Events
//
//    This class reads events according to the TrackClusterSel
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: TrackClusterReader.cc,v 1.2 2009/10/26 20:56:30 mschrode Exp $
//   
#include "TrackClusterReader.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"

#include "TLorentzVector.h"
#include <iostream>

TrackClusterReader::TrackClusterReader(const std::string& configfile, TParameters* p) :
  EventReader(configfile,p),Et_cut_on_track(0),Et_cut_on_cluster(0)
{
  n_trackcluster_events = config_->read<int>("use Track-Cluster events",-1);
  if(n_trackcluster_events == 0) return ;

  Et_cut_on_track = config_->read<double>("Et cut on track",0.0); 
  Et_cut_on_cluster = config_->read<double>("Et cut on cluster",0.0);

  string default_tree_name = config_->read<string>( "Default Tree Name","CalibTree");
  string treename_trackcluster    = config_->read<string>( "Track-Cluster tree", default_tree_name );
  TChain * tchain_trackcluster = new TChain( treename_trackcluster.c_str() );
  vector<string> input_trackcluster = bag_of_string( 
						    config_->read<string>( "Track-Cluster input file", "input/trackcluster.root" ) );
  for (bag_of_string::const_iterator it = input_trackcluster.begin(); it!=input_trackcluster.end(); ++it){
    cout << "...opening root-file " << (*it) << " for Track-Cluster analysis." << endl;
    tchain_trackcluster->Add( it->c_str() );
  }  
  trackcluster.Init( tchain_trackcluster );
}

TrackClusterReader::~TrackClusterReader()
{
}

int TrackClusterReader::readEvents(std::vector<TData*>& data)
{
  if(n_trackcluster_events == 0) return 0;
  //Run Track-cluster stuff
  int nevent = trackcluster.fChain->GetEntries();
  int evt=0;
  for (int i=0;i<nevent;i++) {
    trackcluster.GetEntry(i); 
    if (trackcluster.NobjTowCal>200)
      cerr<<"ERROR: Increase array sizes in TrackClusterSelector; NobjTowCal="
	  <<trackcluster.NobjTowCal<<"!"<<endl;
     
    //Calculate cluster energy (needed for plotting)
    double cluster_energy = 0.0;	   
    for (int n=0; n<trackcluster.NobjTowCal; ++n)
      cluster_energy += trackcluster.TowEt[n];

    if (trackcluster.TrackEt < Et_cut_on_track || cluster_energy < Et_cut_on_cluster) continue;
    TMeasurement* clusterp  = new TJet;
    clusterp->pt  = cluster_energy;
    clusterp->eta = 0;
    clusterp->phi = 0;
    clusterp->E   = 0;
    //Define Track-Cluster event	
    TData_TruthMultMess * tc = new TData_TruthMultMess( 0,
							trackcluster.TrackEt,  			           //truth//
							sqrt(pow(0.5,2)+pow(0.10*trackcluster.TrackEt,2)), //error//
							//trackcluster.EventWeight,                        //weight//
							1.,                                                //weight//
							0,                                                 //params
							0,                                                 //number of free jet param. p. bin
							par_->dummy_parametrization,                          //function
			                                jet_error_param,                                  //error param. function
							clusterp);
    tc->SetType( TrackCluster );
    //Add the towers to the event
    for (int n=0; n<trackcluster.NobjTowCal; ++n){
      //if (trackcluster.TrackEt[n]<Et_cut_on_track)
      //   continue;

      int index=par_->GetBin(par_->GetEtaBin(trackcluster.TowId_eta[n]),par_->GetPhiBin(trackcluster.TowId_phi[n]));
      if (index<0) {
	cerr << "INDEX = "<< index << endl;
	continue;
      }
      //create array with multidimensional measurement
      TMeasurement * mess = new TTower;
      mess->pt = double(trackcluster.TowEt[n]);
      double scale = trackcluster.TowEt[n]/trackcluster.TowE[n];
      mess->EMF = double(trackcluster.TowEm[n]*scale);
      mess->HadF = double(trackcluster.TowHad[n]*scale);
      mess->OutF = double(trackcluster.TowOE[n]*scale);
      mess->eta = double(trackcluster.TowEta[n]);
      mess->phi = double(trackcluster.TowPhi[n]);
      mess->E = double(trackcluster.TowE[n]);
      //mess[7] = double( cos( trackcluster.JetCalPhi-trackcluster.TowPhi[n] ) ); // Projection factor for summing tower Pt

      TData_TruthMess * tower = new TData_TruthMess(index,
						    mess,                                                      //mess//
						    trackcluster.TrackEt*trackcluster.TowEt[n]/cluster_energy, //"truth" for plotting only!//
						    //trackcluster.TrackEterr[n],                              //error//
						    sqrt(pow(0.5,2)+ pow(0.1*trackcluster.TrackEt ,2)),        //error//
						    1.,                                                        //weight//
						    par_->GetTowerParRef( index ),                                //parameter//
						    par_->GetNumberOfTowerParametersPerBin(),                     //number of free cluster param. p. bin//
						    par_->tower_parametrization,                                  //function//
					            tower_error_param                                        //error param.func.//
						    );
      tc->AddMess( tower );
    } 
     
    //Save event
    data.push_back( tc ); 

    if((evt++)%1000==0) cout<<"Track-Cluster Event: "<<evt<<endl;
    if (n_trackcluster_events>=0 && evt>=n_trackcluster_events)
      break;
  }
  return evt;
}

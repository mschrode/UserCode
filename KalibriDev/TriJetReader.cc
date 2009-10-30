//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: TriJetReader.cc,v 1.3 2009/10/26 20:56:30 mschrode Exp $
//   
#include "TriJetReader.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "ToyMC.h"
#include "Parameters.h"

#include <cstdlib>

TriJetReader::TriJetReader(const std::string& configfile, TParameters* p) 
  : EventReader(configfile,p),Et_cut_nplus1Jet(0),Rel_cut_on_nJet(10),n_trijet_events(0)
{
  n_trijet_events = config_->read<int>("use Tri-Jet events",-1);
  if(n_trijet_events == 0) return;

  Et_cut_nplus1Jet = config_->read<double>("Et cut on n+1 Jet",10.0);
  Rel_cut_on_nJet  =  config_->read<double>("Relative n+1 Jet Et Cut",0.2);
  
  string default_tree_name = config_->read<string>("Default Tree Name","CalibTree"); 
  string treename_trijet    = config_->read<string>( "Tri-Jet tree", default_tree_name );
  TChain * tchain_trijet = new TChain( treename_trijet.c_str() );
  vector<string> input_trijet = bag_of_string(config_->read<string>( "Tri-Jet input file", "input/trijet.root" ) );
  for (bag_of_string::const_iterator it = input_trijet.begin(); it!=input_trijet.end(); ++it){
    cout << "...opening root-file " << (*it) << " for Tri-Jet analysis." << endl;
    tchain_trijet->Add( it->c_str() );
  }  
  njet.Init( tchain_trijet );
}

TriJetReader::~TriJetReader()
{
  
}

int TriJetReader::readEvents(std::vector<TData*>& data)
{
  if(n_trijet_events == 0) return 0;
  int injet = 3;
  //Run jet-Jet stuff  
  int nevent = njet.fChain->GetEntries();
  int evt=0;
  for (int i=0;i<nevent;i++) {
    if((i+1)%10000==0) cout<<injet<<"-Jet Event: "<<i+1<<endl;
    njet.fChain->GetEvent(i); 
    if (njet.NobjTow>10000 || njet.NobjJet>100) {
      cerr << "ERROR: Increase array sizes in NJetSelector; NobjTow="
	   << njet.NobjTow<<", NobjJet="<<njet.NobjJet<<"!"<<endl;
      exit(9);
    }
    //--------------
    //  n - Jet
    //--------------
    TData_PtBalance * jj_data[njet.NobjJet];
    jj_data[0] = 0;
    //std::cout << "reading " << njet.NobjJet << " jets\n";

    int nstoredjets = 0;
    for (unsigned int ij = 0; (int)ij<njet.NobjJet; ++ij){
      if(njet.JetPt[ij] < Et_cut_nplus1Jet) continue;
      //Find the jets eta & phi index using the nearest tower to jet axis:
      int jet_index=-1;
      double min_tower_dr = 10.0;
      double em = 0;
      double had = 0;
      double out = 0;
      TLorentzVector Ljet(0,0,0,0);
      Ljet.SetPtEtaPhiE(njet.JetPt[ij],njet.JetEta[ij],njet.JetPhi[ij],njet.JetE[ij]);
      for (int n=0; n<njet.NobjTow; ++n){
        if (njet.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	em += njet.TowEm[n];
	had += njet.TowHad[n];
	out += njet.TowOE[n];
	TLorentzVector Ltower(0,0,0,0);
	Ltower.SetPtEtaPhiE(njet.TowEt[n],njet.TowEta[n],njet.TowPhi[n],njet.TowE[n]);
	double dr = Ltower.DeltaR(Ljet);
	if (dr<min_tower_dr) {
	  jet_index = par_->GetJetBin(par_->GetJetEtaBin(njet.TowId_eta[n]),
				   par_->GetJetPhiBin(njet.TowId_phi[n]));
	  min_tower_dr = dr;
	}
      }
      if (jet_index<0){ 
	cerr<<"WARNING: JJ jet_index = " << jet_index << endl; 
	continue; 
      }
      double * direction = new double[2];
      direction[0] = sin(njet.JetPhi[ij]);
      direction[1] = cos(njet.JetPhi[ij]);
      TMeasurement* jetp  = new TJet;
      jetp->pt  = njet.JetEt[ij];
      jetp->eta = njet.JetEta[ij];
      jetp->phi = njet.JetPhi[ij];
      jetp->E   = njet.JetE[ij];
    //the following is not quite correct, as this factor is different for all towers. These values should be in the n-tupel as well
      double factor =  njet.JetEt[ij] /  njet.JetE[ij];
      jetp->HadF = had * factor;
      jetp->EMF = em * factor;
      jetp->OutF = out * factor;
      //Create an jet/Jet TData event
      jj_data[nstoredjets] = new TData_PtBalance( 
          jet_index * par_->GetNumberOfJetParametersPerBin() + par_->GetNumberOfTowerParameters(),
	  direction,                                     //p_T direction of this jet
	  0.0,                                           //truth//
	  sqrt(pow(0.5,2)+pow(0.10*njet.JetPt[ij],2)),   //error//
	  njet.Weight,                                   //weight//
	  //1.,                                          //weight//
	  par_->GetJetParRef( jet_index ),                  //params
	  par_->GetNumberOfJetParametersPerBin(),           //number of free jet param. p. bin
	  par_->jet_parametrization,                        //function
	  //par_->dummy_parametrization,
          jet_error_param,                               //error param. function
	  jetp                                           //jet momentum for plotting and scale
        );
//cout << "jet "<<nstoredjets<<"'s E="<<njet.JetE[ij]
//     << ", ntower:"<<endl;
      //Add the jet's towers to "jj_data":
      for (int n=0; n<njet.NobjTow; ++n){
        if (njet.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	//if (njet.TowEt[n]<0.01) continue;

	int index = par_->GetBin(par_->GetEtaBin(njet.TowId_eta[n]),
			      par_->GetPhiBin(njet.TowId_phi[n]));
//std::cout << "jet:" << ij << ", towid=" << n << ", bin index:" << index << "\n";
	if (index<0){ cerr<<"WARNING: JJ tower_index = " << index << endl; continue; }

	double relativEt = njet.TowEt[n]/njet.JetEt[ij];  
	//if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
	//This relativeE is used *only* for plotting! Therefore no cuts on this var!
	//create array with multidimensional measurement
	TMeasurement * mess = new TTower;
	mess->pt = double(njet.TowEt[n]);
	double scale = njet.TowEt[n]/njet.TowE[n];
	mess->EMF = double(njet.TowEm[n]*scale);
	mess->HadF = double(njet.TowHad[n]*scale);
	mess->OutF = double(njet.TowOE[n]*scale);
	mess->eta = double(njet.TowEta[n]);
	mess->phi = double(njet.TowPhi[n]);
	mess->E = double(njet.TowE[n]);
	//mess[7] = double( cos( njet.JetCalPhi-njet.TowPhi[n] ) ); // Projection factor for summing tower Pt

	jj_data[nstoredjets]->AddMess(new TData_TruthMess(
	    index,
	    mess,                                                   //mess//
	    njet.JetPt[ij] * relativEt,                             //truth//
	    sqrt(pow(0.5,2)+pow(0.1*njet.JetPt[ij]*relativEt,2)),   //error//
            //1.,                                                   //weight//
	    njet.Weight,                                            //weight//
	    par_->GetTowerParRef( index ),                             //parameter//
	    par_->GetNumberOfTowerParametersPerBin(),                  //number of free tower param. p. bin//
	    par_->tower_parametrization,                               //function//
	    tower_error_param                                      //error param. function//
	  ));
      }
      //Add the jet's tracks to "gj_data":
      for (int n=0; n<njet.NobjTrack; ++n){
        if (njet.Track_jetidx[n]!=(int)ij) continue;//look for ij-jet's tracks

	int track_index = par_->GetTrackBin(par_->GetTrackEtaBin(njet.TrackTowIdEta[n]),
					 par_->GetTrackPhiBin(njet.TrackTowIdPhi[n]));
	if (track_index<0){ cerr<<"WARNING: JJ track_index = " << track_index << endl; continue; }
	//create array with multidimensional measurement
	//TMeasurement * Tmess = new TTrack;
	TTrack * Tmess = new TTrack;
	Tmess->TrackId = int(njet.TrackId[n]);
	Tmess->TowerId = int(njet.TrackTowId[n]);
	Tmess->pt = double(njet.TrackPt[n]);
	double scale = njet.TrackP[n]/njet.TrackPt[n];
	Tmess->EM1 = double(njet.TrackEMC1[n]*scale);
	Tmess->EMF = double(njet.TrackEMC3[n]*scale);
	Tmess->EM5 = double(njet.TrackEMC5[n]*scale);
	Tmess->Had1 = double(njet.TrackHAC1[n]*scale);
	Tmess->HadF = double(njet.TrackHAC3[n]*scale);
	Tmess->Had5 = double(njet.TrackHAC5[n]*scale);
	Tmess->OutF = 0;
	Tmess->DR = double(njet.TrackDR[n]);
	Tmess->DRout = double(njet.TrackDROut[n]);
	Tmess->eta = double(njet.TrackEta[n]);
	Tmess->etaOut = double(njet.TrackEtaOut[n]);
	Tmess->phi = double(njet.TrackPhi[n]);
	Tmess->phiOut = double(njet.TrackPhiOut[n]);
	Tmess->E = double(njet.TrackP[n]);
	Tmess->TrackChi2 = double(njet.TrackChi2[n]);
	Tmess->NValidHits = int(njet.TrackNHits[n]);
	Tmess->MuDR = double(njet.MuDR[n]);
	Tmess->MuDE = double(njet.MuDE[n]);
	//mess[7] = double( cos( njet.JetCalPhi-njet.TowPhi[n] ) ); // Projection factor for summing tower Pt
	//EM+=mess->EMF;
	//F+=mess->pt;
	jj_data[nstoredjets]->AddTrack(new TData_TruthMess(
					      track_index  * par_->GetNumberOfTrackParametersPerBin() + par_->GetNumberOfTowerParameters() + par_->GetNumberOfJetParameters() ,
					      Tmess,                                                    //mess//
					      0,                           //truth//
					      0.05 + 0.00015 * njet.TrackPt[n], //error//
					      1.,                                                      //weight//
					      par_->GetTrackParRef( track_index ),                              //parameter//
					      par_->GetNumberOfTrackParametersPerBin(),                   //number of free tower param. p. bin//
					      par_->track_parametrization,                                //function//
					      track_error_param                                        //error param.func.//
					      ));
      }
      jj_data[nstoredjets]->UseTracks(useTracks_);   //check if track information is sufficient to use Track Parametrization
      
      if(nstoredjets> 0)  
      	jj_data[0]->AddNewMultMess( jj_data[nstoredjets] );
      ++nstoredjets;
    }//loop over all n-jets
    bool goodevent=true;
    if (nstoredjets < injet) goodevent = false;
    if (nstoredjets > injet){
      /*
      for (int i=0; i < nstoredjets; ++i){
	cout<<i<<"-ter Jet Pt: "<<jj_data[i]->GetMess()->pt<<endl;
      }
      */
      //relative Pt cut only works if jets are Pt sorted
      double scale=0;
      for (int i=0; i < injet; ++i){
	scale += jj_data[i]->GetMess()->pt;
      }
      scale /= injet;
      //cout<<"scale: "<<scale<<endl;
      if ( jj_data[injet]->GetMess()->pt > scale*Rel_cut_on_nJet ) goodevent = false;
    }
      /*
      //sort jets. 1st is barrel, 2nd is probe
      if( nstoredjets ==  2) {
      if(std::abs(jj_data[0]->GetMess()[1]) > 1.2) {
      if(std::abs(jj_data[1]->GetMess()[1]) > 1.2) {
      delete jj_data[0];
      continue;
      } else {
      jj_data[0]->ClearMultMess();
      jj_data[1]->AddNewMultMess(jj_data[0]);
      TData_PtBalance* tmp = jj_data[1];
      jj_data[1] = jj_data[0];
      jj_data[0] = tmp; 
      }
      } else if(std::abs(jj_data[1]->GetMess()[1]) < 1.2) {
      //both jets central, roll the dice and swap
      if(rand()/(RAND_MAX+1.0) > 0.5) {
      jj_data[0]->ClearMultMess();
      jj_data[1]->AddNewMultMess(jj_data[0]);
      TData_PtBalance* tmp = jj_data[1];
      jj_data[1] = jj_data[0];
      jj_data[0] = tmp; 
      }
      }
      }    
      */
    if (goodevent) {
      ++evt;    
      data.push_back( jj_data[0] ); 
    } else {
      delete jj_data[0];
    }
    if(evt>=n_trijet_events) break;
  }
  return evt;
}

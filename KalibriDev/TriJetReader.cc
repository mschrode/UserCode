//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: TriJetReader.cc,v 1.4 2009/11/24 16:52:59 stadie Exp $
//   
#include "TriJetReader.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "ToyMC.h"
#include "Parameters.h"
#include "NJetSel.h"

#include <cstdlib>

TriJetReader::TriJetReader(const std::string& configfile, TParameters* p) 
  : EventReader(configfile,p),Et_cut_nplus1Jet(0),Rel_cut_on_nJet(10),n_trijet_events(0)
{
  n_trijet_events = config_->read<int>("use Tri-Jet events",-1);
  if(n_trijet_events == 0) return;

  Et_cut_nplus1Jet = config_->read<double>("Et cut on n+1 Jet",10.0);
  Rel_cut_on_nJet  = config_->read<double>("Relative n+1 Jet Et Cut",0.2);
  
  std::string default_tree_name = config_->read<std::string>("Default Tree Name","CalibTree"); 
  std::string treename_trijet   = config_->read<std::string>( "Tri-Jet tree", default_tree_name );
  TChain * tchain_trijet = new TChain( treename_trijet.c_str() );
  std::vector<std::string> input_trijet = bag_of_string(config_->read<std::string>( "Tri-Jet input file", "input/trijet.root" ) );
  for (bag_of_string::const_iterator it = input_trijet.begin(); it!=input_trijet.end(); ++it){
    std::cout << "...opening root-file " << (*it) << " for Tri-Jet analysis." 
	      << std::endl;
    tchain_trijet->Add( it->c_str() );
  }  
  njet->Init( tchain_trijet );
}

TriJetReader::~TriJetReader()
{
  
}

int TriJetReader::readEvents(std::vector<Event*>& data)
{
  if(n_trijet_events == 0) return 0;
  int injet = 3;
  //Run jet-Jet stuff  
  int nevent = njet->fChain->GetEntries();
  int evt=0;
  for (int i=0;i<nevent;i++) {
    if((i+1)%10000==0) std::cout << injet << "-Jet Event: " << i+1 << std::endl;
    njet->fChain->GetEvent(i); 
    if (njet->NobjTow>10000 || njet->NobjJet>100) {
      std::cerr << "ERROR: Increase array sizes in NJetSelector; NobjTow="
		<< njet->NobjTow<<", NobjJet="<<njet->NobjJet<<"!"<< std::endl;
      exit(9);
    }
    
    std::cerr << "TriJetReader::readEvents: missing routine to create events!\n";
    
    if(evt>=n_trijet_events) break;
  }
  return evt;
}

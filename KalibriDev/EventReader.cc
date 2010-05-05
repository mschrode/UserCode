//
// $Id: EventReader.cc,v 1.12 2010/04/29 13:29:41 stadie Exp $
//
#include "EventReader.h"

#include "ConfigFile.h"
#include "Parameters.h" 
#include "Parametrization.h"
#include "CorFactorsFactory.h"
#include "JetConstraintEvent.h" 
#include "TChain.h"
#include "ToyMC.h"
#include "TTree.h"

#include <dlfcn.h>

unsigned int EventReader::numberOfEventReaders_ = 0;
std::vector<JetConstraintEvent*> EventReader::constraints_;

EventReader::EventReader(const std::string& configfile, TParameters* param) 
  : config_(0),par_(param),corFactorsFactory_(0),cp_(new ConstParametrization())
{
  numberOfEventReaders_++;

  config_ = new ConfigFile(configfile.c_str());
  
  useTracks_ = config_->read<bool>("use Tracks",true);
  if(par_->GetNumberOfTrackParameters() < 1) useTracks_ = false;

  //Error Parametrization...
  //...for tracks:
  track_error_param = par_->track_error_parametrization;
  //...for tower:
  string te = config_->read<string>("tower error parametrization","standard"); 
  if (te=="standard")
    tower_error_param = par_->tower_error_parametrization;
  else if (te=="fast")
    tower_error_param = par_->fast_error_parametrization;
  else if (te=="Jans E parametrization")
    tower_error_param = par_->jans_E_tower_error_parametrization;
  else if(te=="const")
    tower_error_param = par_->const_error_parametrization;
  else if(te=="toy")
    tower_error_param = par_->toy_tower_error_parametrization;
  else if(te=="jet")
    tower_error_param = par_->jet_only_tower_error_parametrization;
  else  
    tower_error_param = par_->tower_error_parametrization;
  //...for jets:
   std::string je = config_->read<string>("jet error parametrization","standard");
  if (je=="standard")
    jet_error_param   = par_->jet_error_parametrization;
  else if (je=="fast")
    jet_error_param   = par_->fast_error_parametrization;
  else if (je=="dummy")
    jet_error_param   = par_->dummy_error_parametrization;
  else if(je=="const")
    jet_error_param = par_->const_error_parametrization;
  else if(je=="toy")
    jet_error_param   = par_->toy_jet_error_parametrization;
  else if(je=="jet et")
    jet_error_param   = par_->jet_only_jet_error_parametrization_et;
  else if(je=="jet energy")
    jet_error_param   = par_->jet_only_jet_error_parametrization_energy;
  else  
    jet_error_param   = par_->jet_error_parametrization;

  std::string jcs = config_->read<string>("jet correction source","");
  std::string jcn = config_->read<string>("jet correction name","");
  
  if(jcs != "") {
    std::string libname = "lib/lib"+jcs+".so";
    void *hndl = dlopen(libname.c_str(), RTLD_NOW);
    if(hndl == NULL){
      std::cerr << "failed to load plugin: " << dlerror() << std::endl;
      exit(-1);
    }
  }
  corFactorsFactory_ = CorFactorsFactory::map[jcn];
  if(jcn !="" && (! corFactorsFactory_)) {
    std::cerr << "Failed to apply correction " << jcn << " from " << jcs << std::endl;
    exit(-1);
  } 
  if(corFactorsFactory_) {
    std::cout << "Jet corrections will be overwritten with " << jcn << " from " << jcs << std::endl; 
  }
  correctToL3_ = config_->read<bool>("correct jets to L3",false);
  correctL2L3_ = config_->read<bool>("correct jets L2L3",false);
  if( correctToL3_ && correctL2L3_ ) {
    std::cerr << "WARNING: Jets are corrected twice (to L3 and L2L3).\n" << std::endl;
    exit(-9);
  }

  // Print info only once for all readers
  if( numberOfEventReaders_ == 1 ) {
    // Track usage
    if(useTracks_)  std::cout<<"Tracks are used to calibrate jets"<< std::endl;
    else std::cout<<"Only Calorimeter information is used"<< std::endl;
    // Correction of jets
    if(correctToL3_) {
      std::cout << "Jets will be corrected to Level 3 (i.e. with L1 * L2 * L3)" << std::endl;
    } else if(correctL2L3_) {
      std::cout << "Jets will be corrected with L2 * L3" << std::endl;
    } 
  }

  if(! constraints_.size() ) {
    std::vector<double> jet_constraint = bag_of<double>(config_->read<std::string>( "jet constraints",""));
    if(jet_constraint.size() % 5 == 0) {
      for(unsigned int i = 0 ; i < jet_constraint.size() ; i += 5) {
	constraints_.push_back(new JetConstraintEvent(jet_constraint[i],jet_constraint[i+1],jet_constraint[i+2],jet_constraint[i+3],jet_constraint[i+4]));
      } 
    } else if(jet_constraint.size() > 1) {
      std::cout << "wrong number of arguments for jet constraint:" << jet_constraint.size() << '\n';
    }
    for(unsigned int i = 0 ; i < constraints_.size() ; ++i) {
      const JetConstraintEvent* jce = constraints_[i];
      std::cout << "adding constraint for jets with " << jce->minEta() << " < |eta| <  " 
		<< jce->maxEta() << " and " << jce->minPt() << " < pt < " << jce->maxPt() 
		<< " with weight " << jce->weight() << "\n";
    }
  }
}

EventReader::~EventReader()
{
  delete config_;
  for(unsigned int i = 0 ; i < constraints_.size() ; ++i) {
    delete constraints_[i];
  }
  constraints_.clear();
}



TTree * EventReader::createTree(const std::string &dataType) const {
  std::string readerName;
  std::string treeName;
  std::vector<std::string> inputFileNames;
  int nEvts = 0;
  if( dataType == "dijet" ) {
    readerName = "DiJetReader";
    treeName = config_->read<string>("Di-Jet tree","CalibTree");
    inputFileNames = bag_of_string(config_->read<std::string>("Di-Jet input file","input/dijet.root"));  
    nEvts = config_->read<int>("use Di-Jet events",-1);
  }

  int inputMode = -1;
  if( inputFileNames[0] == "toy" ) {
    inputMode = 0;
  } else {
    std::string fileEnding = "";
    if( inputFileNames[0].size() > 5 ) {
      fileEnding = inputFileNames[0].substr(inputFileNames[0].size()-5,inputFileNames[0].size());
      if( fileEnding == ".root" ) {
	inputMode = 1;
      }
    }
    if( fileEnding != ".root" && inputFileNames.size() == 1 ) {
      inputMode = 2;
    }
  }
  
  TTree *tree = 0;

  if( inputMode == 0 ) { // Generate Toy MC sample
    std::cout << "\n" << readerName << ": generating ToyMC events\n";
    ToyMC* mc = new ToyMC();
    mc->init(config_);
    mc->print();
    tree = new TTree(treeName.c_str(),dataType.c_str());
    if( dataType == "dijet" ) {
      mc->generateDiJetTree(tree,nEvts);
    }
    delete mc;
  } else if( inputMode == 1 ) { // Open all files listed in configfile
    TChain* chain = new TChain(treeName.c_str()); 
    std::cout << "\n" << readerName << ": opening files\n";
    for(std::vector<std::string>::const_iterator it = inputFileNames.begin();
	it!=inputFileNames.end(); ++it){
      std::cout << " " << (*it) << std::endl;
      chain->Add( it->c_str() );
    }  
    tree = chain;
  } else if( inputMode == 2 ) { // Open all files listed in input file
    std::cout << "\n" << readerName << ": opening files in list '" << inputFileNames[0] << "'\n";
    TChain* chain = new TChain(treeName.c_str()); 
    std::ifstream filelist;
    filelist.open(inputFileNames[0].c_str());
    int nOpenedFiles = 0;
    if( filelist.is_open() ) {
      std::string name = "";
      while( !filelist.eof() ) {
	filelist >> name;
	if( filelist.eof() ) break;
	chain->Add( name.c_str() );
	nOpenedFiles++;
      }
    } else {
      std::cerr << "ERROR opening file '" << inputFileNames[0] << "'\n";
      exit(1);
    }
    filelist.close();
    tree = chain;
    std::cout << "Opened " << nOpenedFiles << " files\n";
  } else {
    tree = new TTree();
    std::cerr << "WARNING: Wrong input file name syntax. No files opened.\n";
  }

  return tree;
}
 
int EventReader::addConstraints(std::vector<Event*>& data) {
  unsigned int n = constraints_.size();
  for(unsigned int i = 0 ; i < n ; ++i) { 
    JetConstraintEvent* jce = constraints_[i];
    std::cout << "added constraint for jets with " << jce->minEta() << " < |eta| <  " 
	      << jce->maxEta() << " and " << jce->minPt() << " < pt < " << jce->maxPt() 
	      << " with weight " << jce->weight() << " and " << jce->nJets() 
	      << " jets " << "\n";
    data.push_back(jce);
  }
  constraints_.clear();
  return n;
}

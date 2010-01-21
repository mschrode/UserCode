//
// $Id: EventReader.cc,v 1.6 2010/01/08 18:23:28 mschrode Exp $
//
#include "EventReader.h"

#include "ConfigFile.h"
#include "Parameters.h" 
#include "CorFactorsFactory.h"
#include "TChain.h"
#include "ToyMC.h"
#include "TTree.h"

#include <dlfcn.h>

unsigned int EventReader::numberOfEventReaders_ = 0;

EventReader::EventReader(const std::string& configfile, TParameters* param) 
  : config_(0),par_(param),corFactorsFactory_(0)
{
  numberOfEventReaders_++;

  config_ = new ConfigFile(configfile.c_str());
  
  useTracks_ = config_->read<bool>("use Tracks",true);
  if(par_->GetNumberOfTrackParameters() < 1) useTracks_ = false;

  // Print info on track usage only once for all readers
  if( numberOfEventReaders_ == 1 ) {
    if(useTracks_)  std::cout<<"Tracks are used to calibrate jets"<< std::endl;
    else std::cout<<"Only Calorimeter information is used"<< std::endl;
  }
  
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

  correctToL3_ = config_->read<bool>("correct jets to L3",false);
}

EventReader::~EventReader()
{
  delete config_;
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
    std::string name = "";
    while( !filelist.eof() ) {
      filelist >> name;
      if( filelist.eof() ) break;
      chain->Add( name.c_str() );
      nOpenedFiles++;
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
 

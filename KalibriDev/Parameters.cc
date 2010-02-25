// $Id: Parameters.cc,v 1.7 2010/02/16 13:33:16 mschrode Exp $

#include <fstream>
#include <cassert>
#include <pwd.h>
#include <unistd.h>
#include <cstdlib>
#include <ctime>
#include <iomanip>


#include "Parameters.h"

#include "TMath.h"

using namespace std;

TParameters* TParameters::instance = 0;


// -----------------------------------------------------------------
Parametrization* TParameters::CreateParametrization(const std::string& name, const ConfigFile& config) {
  if(name == "StepParametrization") {
    return new StepParametrization();
  } else if(name == "StepParametrizationEnergy") {
    return new StepParametrizationEnergy();
  } else if(name == "StepEfracParametrization") {
    return new StepEfracParametrization();
  } else if(name == "StepJetParametrization") {
    return new StepJetParametrization();
  } else if(name == "MyParametrization") {
    return new MyParametrization();
  }  else if(name == "JetMETParametrization") {
    return new JetMETParametrization();
  }  else if(name == "GlobalScaleFactorParametrization") {
    return new GlobalScaleFactorParametrization();
  }  else if(name == "SimpleParametrization") {
    return new SimpleParametrization();
  }  else if(name == "ToyParametrization") {
    return new ToyParametrization();
  }  else if(name == "ToyJetParametrization") {
    return new ToyJetParametrization(); 
  }  else if(name == "ToyStepParametrization") {
    return new ToyStepParametrization();
  }  else if(name == "ToyStepJetParametrization") {
    return new ToyStepJetParametrization();
  } else if(name == "TrackParametrization") {
    return new TrackParametrization();
  } else if(name == "L2L3JetParametrization") {
    return new L2L3JetParametrization();
  } else if(name == "L2L3JetParametrization2") {
    return new L2L3JetParametrization2();
  } else if(name == "L2L3JetTrackParametrization") {
    return new L2L3JetTrackParametrization();
  } else if(name == "ToySimpleInverseParametrization") {
    return new ToySimpleInverseParametrization();
  } else if(name == "SmearFermiTail") {
    return new SmearFermiTail();
  } else if(name == "SmearStepGaussInter") {
    double rMin       = config.read<double>("Response pdf min",0.);
    double rMax       = config.read<double>("Response pdf max",1.8);
    int    rNBins     = config.read<int>("Response pdf nsteps",10);
    double tMin       = config.read<double>("DiJet integration min",0.);
    double tMax       = config.read<double>("DiJet integration max",1.);
    double ptDijetMin = config.read<double>("Et min cut on dijet",0.);
    double ptDijetMax = config.read<double>("Et max cut on dijet",1.);
    std::vector<double> scale = bag_of<double>(config.read<string>("jet parameter scales",""));
    std::vector<double> gaussPar = bag_of<double>(config.read<string>("mean response parameters","1 0"));
    return new SmearStepGaussInter(tMin,tMax,rMin,rMax,rNBins,ptDijetMin,ptDijetMax,scale,gaussPar);
  } else if(name == "SmearCrystalBall") {
    double tMin       = config.read<double>("DiJet integration min",0.);
    double tMax       = config.read<double>("DiJet integration max",1.);
    double ptDijetMin = config.read<double>("Et min cut on dijet",0.);
    double ptDijetMax = config.read<double>("Et max cut on dijet",1.);
    std::vector<double> scale = bag_of<double>(config.read<string>("jet parameter scales",""));
    return new SmearCrystalBall(tMin,tMax,ptDijetMin,ptDijetMax,scale);
  } else if(name == "SmearCrystalBallPt") {
    double tMin       = config.read<double>("DiJet integration min",0.);
    double tMax       = config.read<double>("DiJet integration max",1.);
    double ptDijetMin = config.read<double>("Et min cut on dijet",0.);
    double ptDijetMax = config.read<double>("Et max cut on dijet",1.);
    std::vector<double> scale = bag_of<double>(config.read<string>("jet parameter scales",""));
    return new SmearCrystalBallPt(tMin,tMax,ptDijetMin,ptDijetMax,scale);
  } else if(name == "SmearGauss") {
    double tMin       = config.read<double>("DiJet integration min",0.);
    double tMax       = config.read<double>("DiJet integration max",1.);
    double ptDijetMin = config.read<double>("Et min cut on dijet",0.);
    double ptDijetMax = config.read<double>("Et max cut on dijet",1.);
    std::vector<double> scale = bag_of<double>(config.read<string>("jet parameter scales",""));
    return new SmearGauss(tMin,tMax,ptDijetMin,ptDijetMax,scale);
  } else if(name == "SmearGaussPtBin") {
    double tMin       = config.read<double>("DiJet integration min",0.);
    double tMax       = config.read<double>("DiJet integration max",1.);
    double ptDijetMin = config.read<double>("Et min cut on dijet",0.);
    double ptDijetMax = config.read<double>("Et max cut on dijet",1.);
    std::vector<double> scale = bag_of<double>(config.read<string>("jet parameter scales",""));
    return new SmearGaussPtBin(tMin,tMax,ptDijetMin,ptDijetMax,scale);
  } else if(name == "GroomParametrization") {
    return new GroomParametrization();
  } else if(name == "EtaEtaParametrization") {
    return new EtaEtaParametrization();
  }
  return 0;
}
TParameters* TParameters::CreateParameters(const std::string& configfile) 
{
  static Cleaner cleanup;
  if(  instance != 0  )
  {
    delete instance; 
    instance = 0; 
  }  
 
  ConfigFile config( configfile.c_str() );
  
  string parclass = config.read<string>("Parametrization Class","");
  //create Parameters
  if(parclass == "TStepParameters") {
    parclass = "StepParametrization";
  } else if(parclass == "TMyParameters") {
    parclass = "MyParametrization";
  } else if(parclass == "TStepParametersEnergy") {
    parclass = "StepParametrizationEnergy";
  } else if(parclass == "TStepEfracParameters") {
    parclass = "StepEfracParametrization";
  } else if(parclass == "TJetMETParameters") {
    parclass = "JetMETParametrization";
  }  else if(parclass == "TGlobalScaleFactorParameters") {
    parclass = "GlobalScaleFactorParametrization";
  }  else if(parclass == "TSimpleParameters") {
    parclass = "SimpleParametrization";
  }  else if(parclass == "TToyParameters") {
    parclass = "ToyParametrization";
  }  else if(parclass == "TToyJetParameters") {
    parclass = "ToyJetParametrization";
  }  else if(parclass == "TToyStepParametersEnergy") {
    parclass = "ToyStepParametrizationEnergy";
  } else if(parclass == "StepJetParametrization") {
    parclass = "StepJetParametrization";
  } else if(parclass == "TTrackParameters") {
    parclass = "TrackParametrization";
  } else if(parclass == "SmearParametrizationFermiTail") {
    parclass = "SmearFermiTail";
  } else if(parclass == "SmearParametrizationStepGaussInter") {
    parclass = "SmearStepGaussInter";
  } else if(parclass == "SmearParametrizationCrystalBall") {
    parclass = "SmearCrystalBall";
  } else if(parclass == "SmearParametrizationCrystalBallPt") {
    parclass = "SmearCrystalBallPt";
  } else if(parclass == "SmearParametrizationGauss") {
    parclass = "SmearGauss";
  } else if(parclass == "SmearParametrizationGaussPtBin") {
    parclass = "SmearGaussPtBin";
  }

  Parametrization *param = CreateParametrization(parclass,config);
  if(! param) {
    cerr << "TParameters::CreateParameters: could not instantiate class " << parclass << '\n';
    exit(1);
  }
  instance = new TParameters(param);
  instance->Init(config);
  return instance;
}



// -----------------------------------------------------------------
void TParameters::Init(const ConfigFile& config)
{
  eta_ntwr_used   = config.read<unsigned>("maximum eta twr used",82); 
  eta_granularity = config.read<unsigned>("granularity in eta",1); 
  phi_granularity = config.read<unsigned>("granularity in phi",1); 
  eta_symmetry    = config.read<bool>("symmetry in eta",false);
  eta_granularity_jet = config.read<unsigned>("jet granularity in eta",1); 
  phi_granularity_jet = config.read<unsigned>("jet granularity in phi",1); 
  eta_granularity_track = config.read<unsigned>("track granularity in eta",1); 
  phi_granularity_track = config.read<unsigned>("track granularity in phi",1); 

  if (eta_ntwr_used%2 !=0){
    cerr << "WARNING: Use even number of eta towers! Forced exit."<< endl;    
    exit(1);
  }

  if (phi_ntwr%phi_granularity!=0) {
    cerr << "WARNING: Check phi granularity! Forced exit."<< endl;
    exit(1);
  }
      
  if (eta_symmetry && (eta_granularity!=1 && eta_granularity!=3 &&eta_granularity!=5&&eta_granularity!=11&&
      eta_granularity!=21 && eta_granularity!=41 )){
    cerr << "WARNING: Check eta granularity! Should be 1, 3, 5, 11, 21, or 41: Forced exit."<< endl;
    exit(1);
  }

  if (!eta_symmetry && (eta_granularity!=2 && eta_granularity!=6 &&eta_granularity!=10&&eta_granularity!=22&&
      eta_granularity!=42 && eta_granularity!=82 )){
    cerr << "WARNING: Check eta granularity! Should be 2, 6, 10, 22, 42 or 82: Forced exit."<< endl;
    exit(1);
  }

  start_values = bag_of<double>(config.read<string>("start values","")); 
  if ( start_values.size()< p->nTowerPars()){
    cerr<< "ERROR: Number of start values and free parameters does not match!"<<endl
        << "       There must be at least " << p->nTowerPars() << " parameters!" << endl;
    exit(2);    
  }

  jet_start_values = bag_of<double>(config.read<string>("jet start values",""));
  if ( jet_start_values.size()< p->nJetPars()){
    cerr<< "ERROR: Number of jet start values and free jet parameters does not match!"<<endl
        << "       There must be at least " << p->nJetPars() << " parameters!" << endl;
    exit(3);
  }
  track_start_values = bag_of<double>(config.read<string>("track start values","")); 
  if ( track_start_values.size()< p->nTrackPars()){
    cerr<< "ERROR: Number of track start values and free track parameters does not match!"<<endl
        << "       There must be at least " << p->nTrackPars() << " parameters!" << endl;
    exit(3);
  }
  global_jet_start_values = bag_of<double>(config.read<string>("global jet start values","")); 
  if( global_jet_start_values.size() < p->nGlobalJetPars() ) {
    cerr<< "ERROR: Number of global jet start values and free global jet parameters does not match!"<<endl
        << "       There must be at least " << p->nGlobalJetPars() << " parameters!" << endl;
    exit(3);
  }

  // Initialize storage for parameter values and errors
  k = new double[GetNumberOfParameters()];
  parErrors_ = new double[GetNumberOfParameters()];
  parGCorr_ = new double[GetNumberOfParameters()];
  parCov_ = new double[(GetNumberOfParameters()*GetNumberOfParameters()+GetNumberOfParameters())/2];
  trackEff = new double[169];

  for(int i = 0; i < GetNumberOfParameters(); i++) {
    k[i] = 0.;
    parErrors_[i] = 0.;
    parGCorr_[i] = 0.;
  }
  for(int i = 0; i < (GetNumberOfParameters()*GetNumberOfParameters()+GetNumberOfParameters())/2; i++) {
    parCov_[i] = 0.;
  }
  isFixedPar_ = std::vector<bool>(GetNumberOfParameters(),false); 

  for (unsigned int bin=0; bin<eta_granularity*phi_granularity; ++bin){
    for (unsigned int tp=0; tp < p->nTowerPars(); ++tp){
      //step[ bin*free_pars_per_bin + tp ]   = step_sizes[ tp ];
      k[ bin*p->nTowerPars() + tp ] = start_values[ tp ];
      parErrors_[ bin*p->nTowerPars() + tp ] = 0.0;
    }
  }
  for (unsigned int bin=0; bin<eta_granularity_jet*phi_granularity_jet; ++bin){
    for (unsigned int jp=0; jp < p->nJetPars(); ++jp){
      int i = GetNumberOfTowerParameters() + bin*p->nJetPars() + jp;   
      k[i] = jet_start_values[jp];
    }
  }
  for (unsigned int bin=0; bin<eta_granularity_track*phi_granularity_track; ++bin){
    for (unsigned int trp=0; trp < p->nTrackPars(); ++trp){
      int i = GetNumberOfTowerParameters() + GetNumberOfJetParameters() + bin*p->nTrackPars() + trp;   
      k[i] = track_start_values[trp];
    }
  }

  for(int etabin=0; etabin<13; ++etabin)
    {
      for(int ptbin=0; ptbin<13; ++ptbin)
	{
	  trackEff[13*etabin+ptbin] = 1;
	}
    }

  for (unsigned int gjp = 0 ; gjp < p->nGlobalJetPars() ; ++gjp){
    int i = GetNumberOfTowerParameters() + GetNumberOfJetParameters() + GetNumberOfTrackParameters() + gjp;   
    k[i] = global_jet_start_values[gjp];
  }
  
  // read predefined calibration contants from cfi
  // or txt file depending on the ending of the name
  std::vector<std::string> inputCalibration = bag_of_string(config.read<string>("input calibration",";"));

  // Check whether start values for calibration are to be read
  // from file
  if( inputCalibration.empty() ) {
    cout << "Using calibration start values from config file.\n";
  } else {
    cout << "Using calibration start values from file ";

    // Check whether calibration constants are given in one of the
    // Kalibri formats or in official JetMET format and prepare
    // vector for further processing
    std::string inputFormat = "UNKNOWN";
    if( inputCalibration.front() == "Kalibri" ) {
      inputFormat = "Kalibri";
      inputCalibration.erase(inputCalibration.begin());
    }
    else if( inputCalibration.front() == "JetMET" ) {
      inputFormat = "JetMET";
      inputCalibration.erase(inputCalibration.begin());
    }
    else if( inputCalibration.size() == 1 ) { // For backward compatibility
      inputFormat = "Kalibri";
    }

    // Read calibration constants
    if( inputCalibration.size() == 0 ) {
      std::cerr << "\nWARNING: No file name specified to read start values from.\n"; 
      std::cerr << "         Using start values from config file.\n";
    }
    else {
      if( inputFormat == "Kalibri" ) {
	// Read predefined calibration contants from cfi
	// or txt file depending on the ending of the name
	std::string inputFileName = inputCalibration.front();
	if( !inputFileName.substr(inputFileName.rfind(".")+1).compare("cfi") ) {
	  cout << inputFileName << std::endl;
	  readCalibrationCfi(inputFileName); 
	} 
	else if( !inputFileName.substr(inputFileName.rfind(".")+1).compare("txt") ) {
	  cout << inputFileName << std::endl;
	  readCalibrationTxt(inputFileName); 
	}
	else {
	  cerr << "\nERROR: Unknown file format: '" ;
	  cerr << inputFileName.substr(inputFileName.rfind(".")) << "'\n";
	  cerr << "       Using start values from config file.\n";
	}
      }
      else if( inputFormat == "JetMET" ) {
	cout << ":\n";
	readCalibrationJetMET(inputCalibration); 
      }
      else {
	std::cerr << "\nWARNING: Unknown input format.\n"; 
	std::cerr << "         Using start values from config file.\n";
      }
    }
  }

  std::string trackEffFileName = config.read<string>("track efficiency","");
  if(!trackEffFileName.empty()){
  cout << "Reading Track Efficiency from file '" << trackEffFileName << endl;
  readTrackEffTxt(trackEffFileName);
  }

  parNames_ = std::vector<std::string>(GetNumberOfParameters(),"");
  if( p->name() == "SmearCrystalBall" ) {
    parNames_[0] = "#sigma";
    parNames_[1] = "#alpha";
    parNames_[2] = "n";
    parNames_[3] = "m";
  }
}


TParameters::~TParameters() {
  delete p;
  delete [] k;
  delete [] parErrors_;
  delete [] parGCorr_;
  delete [] parCov_;
  delete [] trackEff;
}



// -----------------------------------------------------------------
std::string TParameters::trim(std::string const& source, char const* delims) 
{
  std::string result(source);
  std::string::size_type index = result.find_last_not_of(delims);
  if(index != std::string::npos)
    result.erase(++index);

  index = result.find_first_not_of(delims);
  if(index != std::string::npos)
    result.erase(0, index);
  else
    result.erase();
  
  //replace all "," by " "  :
  std::string::size_type  pos = result.find(",");
  while(pos != string::npos) {
    result.replace(pos,1," ");
    pos = result.find(",",pos);
  }
    
  return result;
}



//!  \brief Read predefined calibration constants from txt file 
//!
//! fills start parameters for fit when read from txt file; expects 
//! 72 lines for 72 bins in phi for each eta bin ranging from -41
//! to 41 (skipping the 0) and the following parameter format:
//! maxEta minEta nPar towerParameters jetParameters separated by
//! blanks
// ---------------------------------------------------------------
void TParameters::readCalibrationTxt(std::string const& configFile)
{
  std::ifstream file(configFile.c_str());
  std::string line; // buffer line

  int      etaBin=-42;
  unsigned phiBin=  0;
  //  unsigned iLines=  0;
  while( std::getline(file,line) ){
    phiBin = 1; // Only 1 phibin is written to the file
    // determine phi bin on the fly
    //    phiBin=(iLines%72)+1;      // phi counts from 1...72 for each eta bin
//     ++iLines;
    // determine eta bin on the fly
    if(phiBin==1) ++etaBin; // increas etaValue by for the first phi bin
    if(etaBin==0) ++etaBin; // and don't forget to skip the 0

    //    cout << "etaBin: " << etaBin << " :: " << "phiBin: " << phiBin << endl;

    // buffers for input parameters
    unsigned nPar=0; //this is not needed but read out for control reasons 
    double etaMax=0; //this is not needed but read out for control reasons  
    double etaMin=0; //this is not needed but read out for control reasons 
    double etMin=0;
    double etMax=0;
    std::vector<double> twrPars, jetPars, trkPars,globaljetPars;
    unsigned entry=0; // controls which parameter is to filled
    while( line.length()>line.substr(0, line.find(" ")).size() ){
      if( 0<line.find(" ")){
   	// extract value
	switch(++entry){
	case 1 : etaMin = std::atof( line.substr(0, line.find(" ")).c_str() ); 
	  break;
	case 2 : etaMax = std::atof( line.substr(0, line.find(" ")).c_str() ); 
	  break; 
	case 3 : nPar   = std::atoi( line.substr(0, line.find(" ")).c_str() ); 
	  break; 
	case 4 : etMin  = std::atoi( line.substr(0, line.find(" ")).c_str() ); 
	  break; 
	case 5 : etMax  = std::atoi( line.substr(0, line.find(" ")).c_str() ); 
	  break; 
	default:
	  if((entry-5)<=p->nTowerPars()){
	    twrPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));
	  }
	  else if((entry-5)<=p->nTowerPars()+p->nJetPars()){
	    jetPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));
	  }
	  else if((entry-5)<=p->nTowerPars()+p->nJetPars()+p->nTrackPars()) {
	    trkPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));
	  } else {
	    globaljetPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));
	  }
	  break;
	}
	// cut string
	line = line.substr(line.find(" "));
      }
      else{
	//cut string
	if(line.find(" ")<std::string::npos){
	  line = line.substr(line.find(" ")+1);
	}
      }
    }
    // catch last character
    //trkPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));
    globaljetPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));

    // fill parameters
    for(phiBin = 1; phiBin <= 72; phiBin++) {
      int towerIdx = GetBin(GetEtaBin(etaBin),GetPhiBin(phiBin));
      if( towerIdx<0 ) continue;
      for (unsigned n=0; n< p->nTowerPars(); ++n) {
	k[towerIdx*p->nTowerPars()+n] = twrPars[n];
	//e[towerIdx*p->nTowerPars()+n] = NOT_READ_OUT;
      }
      int jetIdx = GetJetBin(GetJetEtaBin(etaBin),GetJetPhiBin(phiBin));
      if( jetIdx<0 ) continue;
      for (unsigned n=0; n<p->nJetPars(); ++n) {
	k[GetNumberOfTowerParameters()+jetIdx*p->nJetPars()+n] = jetPars[n];
	//e[GetNumberOfTowerParameters()+jetIdx*p->nJetPars()+n] = NOT_READ_OUT;
      }
      int trackIdx = GetTrackBin(GetTrackEtaBin(etaBin),GetTrackPhiBin(phiBin));
      if( trackIdx<0 ) continue;
      for (unsigned n=0; n<p->nTrackPars(); ++n) {
	k[GetNumberOfTowerParameters()+GetNumberOfJetParameters()+trackIdx*p->nTrackPars()+n] = trkPars[n];
	//e[GetNumberOfTowerParameters()+jetIdx*p->nJetPars()+n] = NOT_READ_OUT;
      }
      for (unsigned n=0; n<p->nGlobalJetPars(); ++n) {
	k[GetNumberOfTowerParameters()+GetNumberOfJetParameters()+GetNumberOfTrackParameters() +n] = globaljetPars[n];
      }
    }
  }
}



//!  \brief Read predefined calibration constants from cfi file 
// -----------------------------------------------------------------
void TParameters::readCalibrationCfi(std::string const& configFile)
{
  std::ifstream file(configFile.c_str());

  std::string line, name;
  char * dummy = new char[28];
  std::vector<int> eta, phi;
  std::vector<double> param[p->nTowerPars()], error[p->nTowerPars()];
  std::vector<int> eta_jet, phi_jet;
  std::vector<double> param_jet[p->nJetPars()], error_jet[p->nJetPars()];
  std::vector<int> eta_track, phi_track;
  std::vector<double> param_track[p->nTrackPars()], error_track[p->nTrackPars()];
  std::vector<double> param_globaljet, error_globaljet;
  int posEqual;
  while (std::getline(file,line)) {
    if (! line.length()) continue;
    if( line.find("#") != string::npos) continue;
    
    //Read Tower Calibration: ---------------------------------------------------
    //if ( line.find("module ccctm = CalibratedCaloTowerMaker") != string::npos ) {
    if ( line.find("block TowerCalibConstants = {") != string::npos ) {
      while( std::getline(file,line) ) {
	if( line.find("block") != string::npos) break;
	posEqual=line.find('=');
	name  = line.substr(0,posEqual);
	if( name.find("TowMapEta") != string::npos) 
  	  eta = bag_of<int>(trim(line.substr(posEqual+1)));
	if( name.find("TowMapPhi") != string::npos) 
  	  phi = bag_of<int>(trim(line.substr(posEqual+1)));
	for (unsigned i=0; i < p->nTowerPars() ; ++i) {
	  sprintf(dummy,"TowerParam%d ",i);
	  if( name.find(dummy) != string::npos) 
  	    param[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	  sprintf(dummy,"TowerError%d ",i);
	  if( name.find(dummy) != string::npos) 
  	    error[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	}
      }
    }
    //Read Jet Calibration --------------------------------------------------------
    if ( line.find("block JetCalibConstants = {") != string::npos ) {
      while (std::getline(file,line)) {
	if( line.find("block") != string::npos) break;
	posEqual=line.find('=');
	name  = line.substr(0,posEqual);
	//std::cout << name << ".\n";
	if( name.find("JetMapEta") != string::npos) 
	  eta_jet = bag_of<int>(trim(line.substr(posEqual+1)));
	if( name.find("JetMapPhi") != string::npos) 
  	  phi_jet = bag_of<int>(trim(line.substr(posEqual+1)));
	for (unsigned i=0; i<p->nJetPars(); ++i) {
	  sprintf(dummy,"JetParam%d ",i);
	  if( name.find(dummy) != string::npos) 
  	    param_jet[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	  sprintf(dummy,"JetError%d ",i);
	  if( name.find(dummy) != string::npos) 
  	    error_jet[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	}
	if( name.find("GlobalJetParams") != string::npos) 
	  param_globaljet = bag_of<double>(trim(line.substr(posEqual+1)));
	if( name.find("GlobalJetErrors") != string::npos) 
	  error_globaljet = bag_of<double>(trim(line.substr(posEqual+1)));
      }
    }
    //Read Track Calibration --------------------------------------------------------
    if ( line.find("block TrackCalibConstants = {") != string::npos ) {
      while (std::getline(file,line)) {
	if( line.find("block") != string::npos) break;
	posEqual=line.find('=');
	name  = line.substr(0,posEqual);
	//std::cout << name << ".\n";
	if( name.find("TrackMapEta") != string::npos) 
	  eta_track = bag_of<int>(trim(line.substr(posEqual+1)));
	if( name.find("TrackMapPhi") != string::npos) 
  	  phi_track = bag_of<int>(trim(line.substr(posEqual+1)));
	for (unsigned i=0; i<p->nTrackPars(); ++i) {
	  sprintf(dummy,"TrackParam%d ",i);
	  if( name.find(dummy) != string::npos) 
  	    param_track[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	  sprintf(dummy,"TrackError%d ",i);
	  if( name.find(dummy) != string::npos) 
  	    error_track[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	}
      }
    }
  }
  //check if the read calibration is ok:
  bool ok=eta.size()==phi.size();
  for (unsigned i=0; i < p->nTowerPars(); ++i){
    ok &= eta.size()==param[i].size();
    ok &= eta.size()==error[i].size();
  }
  //fill tower parameters and errors:  
  if (ok) {
    for (unsigned i=0; i<eta.size(); ++i){
      int index = GetBin(GetEtaBin(eta[i]),GetPhiBin(phi[i]));
      if (index<0) continue;
      for (unsigned n=0; n< p->nTowerPars(); ++n) {
        k[index*p->nTowerPars()+n] = param[n][i];
        parErrors_[index*p->nTowerPars()+n] = error[n][i];
      }
    }
  }

  //check if the read calibration is ok:
  ok=eta_jet.size()==phi_jet.size();
  for (unsigned i=0; i<p->nJetPars(); ++i){
    ok &= eta_jet.size()==param_jet[i].size();
    ok &= eta_jet.size()==error_jet[i].size();
  } 
  //fill Jet parameters and errors:  
  if (ok) {
    for (unsigned i=0; i<eta_jet.size(); ++i){
      int index = GetJetBin(GetJetEtaBin(eta_jet[i]),GetJetPhiBin(phi_jet[i]));
      if (index<0) continue;
      for (unsigned n=0; n<p->nJetPars(); ++n) {
        k[GetNumberOfTowerParameters() + index*p->nJetPars()+n] = param_jet[n][i];
        parErrors_[GetNumberOfTowerParameters() + index*p->nJetPars()+n] = error_jet[n][i];
      }
    }
  }

  //check if the read calibration is ok:
  ok=eta_track.size()==phi_track.size();
  for (unsigned i=0; i<p->nTrackPars(); ++i){
    ok &= eta_track.size()==param_track[i].size();
    ok &= eta_track.size()==error_track[i].size();
  } 
  //fill Track parameters and errors:  
  if (ok) {
    for (unsigned i=0; i<eta_track.size(); ++i){
      int index = GetTrackBin(GetTrackEtaBin(eta_track[i]),GetTrackPhiBin(phi_track[i]));
      if (index<0) continue;
      for (unsigned n=0; n<p->nTrackPars(); ++n) {
        k[GetNumberOfTowerParameters() + GetNumberOfJetParameters () + index*p->nTrackPars()+n] = param_track[n][i];
        parErrors_[GetNumberOfTowerParameters() + GetNumberOfJetParameters () + index*p->nTrackPars()+n] = error_track[n][i];
      }
    }
  }
  ok =  (param_globaljet.size() == p->nGlobalJetPars());
  //fill global Jet parameters and errors:  
  if (ok) {
    for (unsigned n=0; n<p->nGlobalJetPars(); ++n) {
      k[GetNumberOfTowerParameters() + GetNumberOfJetParameters () + GetNumberOfTrackParameters()+n] = param_globaljet[n];
      parErrors_[GetNumberOfTowerParameters() + GetNumberOfJetParameters () + GetNumberOfTrackParameters()+n] = error_globaljet[n];
    }
  }
  delete[] dummy;
}



//!  \brief Read correction factors in CondDB format
// -----------------------------------------------------------------
void TParameters::readCalibrationJetMET(const std::vector<std::string>& inputFileNames) {
  std::string corrL2FileName;
  std::string corrL3FileName;
  for(size_t i = 0; i < inputFileNames.size(); i++) {
    if( inputFileNames.at(i).find("L2") != std::string::npos ) 
      corrL2FileName = inputFileNames.at(i);
    else if( inputFileNames.at(i).find("L3") != std::string::npos ) 
      corrL3FileName = inputFileNames.at(i);
  }

  if( !corrL2FileName.empty() ) {
    std::cout << "  L2: " << corrL2FileName << std::endl;
    readCalibrationJetMETL2(corrL2FileName);
  }
  if( !corrL3FileName.empty() ) {
    std::cout << "  L3: " << corrL3FileName << std::endl;
    readCalibrationJetMETL3(corrL3FileName);
  }
}



//!  \brief Read L2 correction factors in CondDB format
//!
//!  Read parameters of L2 correction from
//!  txt file in CondDB format i.e.
//!  <tt>etaMin etaMax nPar EtMin EtMax Par1 Par2 Par3 Par4 Par5 Par6</tt>.
//!  If there are more than 6 L2 parameters in the file,
//!  they are ignored.
//!
//!  The pt ranges of validity are not considered.
//!
//!  In case some eta bins are missing, default parameter
//!  values 1 0 0 are assumed.
//!
//!  \note There are scaling factors for the L2 parameters
//!  in the \p L2L3JetParametrization . The parameter
//!  values read by this method are scaled accordingly.
// -----------------------------------------------------------------
void TParameters::readCalibrationJetMETL2(const std::string& inputFileName) {
  // There are scaling factors in "L2L3JetParametrization"
  std::vector<double> scale(6,1.);
  scale.at(1) = 10.;
  scale.at(2) = 100.;
  scale.at(3) = 100.;
  scale.at(4) = 100.;

  std::vector< std::vector<double> > parL2;

  std::ifstream file;
  file.open(inputFileName.c_str());

  int    etaBin = -41;
  float  etaMin = 0.;
  float  etaMax = 0.;
  int      nPar = 0;
  float     val = -1.;  // Needs float precision as etaEdge() has
                        // float precision; do we need it there?

  if( file.is_open() ) {
    while( !file.eof() && etaBin < 42 ) {
      if( etaBin == 0 ) etaBin++;          // No bin index 0
      file >> etaMin;                      // Eta min
      file >> etaMax;                      // Eta max
      val = 0.;
      file >> val;                         // Number of values following
      if( val != 0 ) {                     // Avoid reading of empty last line
	nPar = static_cast<int>(val - 2);  // Number of L2 parameters in file
	std::vector<double> par(GetNumberOfJetParametersPerBin(),0.);  // Storage of the L2 parameters in this bin
	file >> val;                       // Et min
	file >> val;                       // Et max
	// Store L2 parameters
	for(int i = 0; i < GetNumberOfJetParametersPerBin(); i++) {
	  file >> val;
	  par.at(i) = scale.at(i) * val;
	}
	// In case of different numbers of parameters in
	// JetMET and Kalibri L2 parametrization
	for(int i = 0; i < nPar - GetNumberOfJetParametersPerBin(); i++) {
	  file >> val;
	}      
	// In case some eta bin is missing,
	// add default parameters
	while( etaBin < 42 && 
	       etaMin != etaLowerEdge(etaBin) && 
	       etaMax != etaUpperEdge(etaBin)    ) {
	  std::cout << "    WARNING: No parameters for eta bin " << etaBin;
	  std::cout << "; using default parameters instead.\n";

	  std::vector<double> defaultPar(GetNumberOfJetParametersPerBin(),0.);
	  for(int i = 0; i < GetNumberOfJetParametersPerBin(); i++) {
	    if( i == 0 ) defaultPar.at(i) = 1.;
	    else         defaultPar.at(i) = 0.;
	  }
	  parL2.push_back(defaultPar);
	  etaBin++;
	}
      
	if( etaBin < 42 ) {
	  parL2.push_back(par);
	}
	etaBin++;
      }
    }
  }
  file.close();

  // In case last eta bins are missing,
  // add default parameters  
  while( etaBin < 42 ) {
    if( etaBin == 0 ) etaBin++;          // No bin index 0

    std::cout << "    WARNING: No parameters for eta bin " << etaBin;
    std::cout << "; using default parameters instead.\n";

    std::vector<double> defaultPar(GetNumberOfJetParametersPerBin(),0.);
    for(int i = 0; i < GetNumberOfJetParametersPerBin(); i++) {
      if( i == 0 ) defaultPar.at(i) = 1.;
      else         defaultPar.at(i) = 0.;
    }
    parL2.push_back(defaultPar);
    etaBin++;
  }

  // Write read constants to array
  int etaIdx = 0;
  for(etaBin = -41; etaBin <= 41; etaBin++, etaIdx++) {
    if( etaBin == 0 ) etaBin++;
    for(int phiBin = 1; phiBin <= 72; phiBin++) {
      int jetIdx = GetJetBin(GetJetEtaBin(etaBin),GetJetPhiBin(phiBin));
      if( jetIdx<0 ) continue;
      for(int i = 0; i < GetNumberOfJetParametersPerBin(); i++) {
	k[GetNumberOfTowerParameters() +
	  jetIdx*GetNumberOfJetParametersPerBin() + i] = parL2.at(etaIdx).at(i);
      }
    }
  }
}



//!  \brief Read L3 correction factors in CondDB format
//!
//!  Read parameters of L3 correction from
//!  txt file in CondDB format i.e.
//!  <tt>etaMin etaMax nPar EtMin EtMax Par1 Par2 Par3 Par4</tt>
//!
//!  The pt ranges of validity are not considered.
// -----------------------------------------------------------------
void TParameters::readCalibrationJetMETL3(const std::string& inputFileName) {
  std::ifstream file;
  file.open(inputFileName.c_str());

  std::vector<double> parL3;   // Storage for the L3 parameters
  double val  = -1.;
  int n = 0;

  if( file.is_open() ) {
    file >> val;                       // Eta min
    file >> val;                       // Eta max
    file >> val;                       // Number of values following
    n = static_cast<int>(val - 2);     // Number of L3 parameters
    file >> val;                       // Et min
    file >> val;                       // Et max
    for(int i = 0; i < n; i++) {       // Store L3 parameters
      file >> val;
      parL3.push_back(val);
    }
  }
  file.close();

  if( n == GetNumberOfGlobalJetParameters() ) {
    for(int i = 0; i < GetNumberOfGlobalJetParameters(); i++) {
      k[GetNumberOfTowerParameters() + 
	GetNumberOfJetParameters() + 
	GetNumberOfTrackParameters() + i] = parL3[i];
    }
    
  } else {
    std::cerr << "ERROR: Number of read global jet parameters too small.\n";
    std::cerr << "       Using start values from config file.\n";
  }
}



//!  \brief Read track efficiency from txt file
// -----------------------------------------------------------------
void TParameters::readTrackEffTxt(std::string const& configFile)
{
  // ---------------------------------------------------------------
  //read Track Efficiency as used in JPT Algorithm
  // ---------------------------------------------------------------
  std::ifstream file(configFile.c_str());
  std::string line; // buffer line
  int      etaBin=  -1;
  int      ptBin=  0;
  unsigned iLines=  0; 
  if(! file) cout<<configFile.c_str()<<" does not exist"<<endl;
  while( std::getline(file,line) ){
    // determine pt bin on the fly
    ptBin=(iLines%13);      // pt counts from 0..12 for each eta bin
    // determine eta bin on the fly
    if(iLines%13==0) ++etaBin; // increas etaValue by 1 each 13 lines...
    ++iLines;

    //cout << "etaBin: " << etaBin << " :: " << "ptBin: " << ptBin << endl;

    unsigned entry=0; // controls which parameter is to filled
    while( line.length()>line.substr(0, line.find(" ")).size() ){
      if( 0<line.find(" ")){
   	// extract value
	switch(++entry){
	case 1 : if( std::atof( line.substr(0, line.find(" ")).c_str() ) != etaBin) cout<<"error in Track efficiency file"<<endl;
	  break;
	case 2 : if( std::atof( line.substr(0, line.find(" ")).c_str() ) != ptBin) cout<<"error in Track efficiency file"<<endl;
	  break; 
	case 3 :  break; 
	default:
	  trackEff[iLines-1] = std::atof( line.substr(line.find(" ")).c_str() );
	  //cout<<iLines-1<<"  "<<etaBin<<"  "<<ptBin<<"  :  "<< trackEff[iLines-1]<<endl;
	  break;
	}
	// cut string
	line = line.substr(line.find(" "));
      }
      else{
	//cut string
	if(line.find(" ")<std::string::npos){
	  line = line.substr(line.find(" ")+1);
	}
      }
    }
  }
}



// -----------------------------------------------------------------
int TParameters::GetEtaBin(int eta_id, int etagranu, int phigranu, bool etasym) const
{  
  assert(eta_id != 0);
  assert(eta_id <= 41);
  assert(eta_id >= -41);
  //This function knows the number of wanted eta-bins and returns 
  //in which eta-bin the tower with eta-ID "eta_id" is located.
  //Case 1 bin:
  //cout << "eta="<<eta_id<<", etagranu:"<< etagranu<< ", eta_ntwr_used:"<< eta_ntwr_used<< "eta sym:"
  //     << etasym << endl;
  if(etagranu<=1) return 0;
  if(etagranu==2) return (eta_id < 0) ? 0 : 1;

  //check if tower is within wanted etarange:
  if( etasym && std::abs(eta_id)*2>(int)eta_ntwr_used)   return -2; 
  //calculate an index:
  unsigned index=(unsigned)(41+eta_id);
  if (eta_id>0) --index;
  if (index<0 || index>81) return -3;
  //std::cout << "eta id:" << eta_id << " index:" << index << '\n';
  //Lookup tables for some binning scenarios:                      eta 1   ->|<-  eta 2                                                |<- eta 3                         
  static const unsigned ta_42[82]={ 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,20,20,21,21,22,23,23,24,24,25,25,26,26,27,27,28,28,29,29,30,30,31,31,32,32,33,33,34,34,35,35,36,36,37,37,38,38,39,39,40,40,41,41};
  static const unsigned ta_22[82]={ 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9,10,10,10,10,11,11,11,11,12,12,12,12,13,13,13,13,14,14,14,14,15,15,15,15,16,16,16,16,17,17,17,17,18,18,18,18,19,19,19,19,20,20,20,21,21};
  static const unsigned ta_10[82]={ 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9};
  static const unsigned ta_6[82]= { 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5};
  
  static const unsigned ts_21[41]={ 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,20,20};
  static const unsigned ts_11[41]={ 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9,10,10};
  static const unsigned ts_5[41]= { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4};
  static const unsigned ts_3[41]= { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2};
  if (!etasym){
    if (etagranu==82) return index;
    else if (etagranu==42) return ta_42[index];
    else if (etagranu==22) return ta_22[index];
    else if (etagranu==10) return ta_10[index];
    else if (etagranu==6) return ta_6[index];
  } else {
    if (etagranu==41) return std::abs(eta_id)-1;
    else if (etagranu==21) return ts_21[std::abs(eta_id)-1];
    else if (etagranu==11) return ts_11[abs(eta_id)-1];
    else if (etagranu== 5) return ts_5[ abs(eta_id)-1];
    else if (etagranu== 3) return ts_3[ abs(eta_id)-1];
  }
  
  //Default value, should never be returned!
  return -4;
}



//!  \brief Return upper or lower eta edge
//!  \param etaBin Index of eta bin,
//!         \f$ -41 \leq \texttt{etaBin} \leq 41 \f$
//!  \param lowerEdge If true, return value of lower
//!                   edge else of upper edge
// -----------------------------------------------------------------
float TParameters::etaEdge(int const etaBin, bool lowerEdge)
{
  // return eta bin - eta edge mappting
  switch(etaBin){
  case -41: return (lowerEdge ? -5.191 : -4.889); break;
  case -40: return (lowerEdge ? -4.889 : -4.716); break;
  case -39: return (lowerEdge ? -4.716 : -4.538); break;
  case -38: return (lowerEdge ? -4.538 : -4.363); break;
  case -37: return (lowerEdge ? -4.363 : -4.191); break;
  case -36: return (lowerEdge ? -4.191 : -4.013); break;
  case -35: return (lowerEdge ? -4.013 : -3.839); break;
  case -34: return (lowerEdge ? -3.839 : -3.664); break;
  case -33: return (lowerEdge ? -3.664 : -3.489); break;
  case -32: return (lowerEdge ? -3.489 : -3.314); break;
  case -31: return (lowerEdge ? -3.314 : -3.139); break;
  case -30: return (lowerEdge ? -3.139 : -2.964); break;
  case -29: return (lowerEdge ? -2.964 : -2.853); break; 
  case -28: return (lowerEdge ? -2.853 :  -2.65); break;
  case -27: return (lowerEdge ?  -2.65 :   -2.5); break;
  case -26: return (lowerEdge ?   -2.5 : -2.322); break;
  case -25: return (lowerEdge ? -2.322 : -2.172); break;
  case -24: return (lowerEdge ? -2.172 : -2.043); break;
  case -23: return (lowerEdge ? -2.043 :  -1.93); break;
  case -22: return (lowerEdge ?  -1.93 :  -1.83); break;
  case -21: return (lowerEdge ?  -1.83 :  -1.74); break;
  case -20: return (lowerEdge ?  -1.74 : -1.653); break;
  case -19: return (lowerEdge ? -1.653 : -1.566); break;
  case -18: return (lowerEdge ? -1.566 : -1.479); break;
  case -17: return (lowerEdge ? -1.479 : -1.392); break;
  case -16: return (lowerEdge ? -1.392 : -1.305); break;
  case -15: return (lowerEdge ? -1.305 : -1.218); break;
  case -14: return (lowerEdge ? -1.218 : -1.131); break;
  case -13: return (lowerEdge ? -1.131 : -1.044); break;
  case -12: return (lowerEdge ? -1.044 : -0.957); break;
  case -11: return (lowerEdge ? -0.957 : -0.879); break;
  case -10: return (lowerEdge ? -0.879 : -0.783); break;
  case  -9: return (lowerEdge ? -0.783 : -0.696); break;
  case  -8: return (lowerEdge ? -0.696 : -0.609); break;
  case  -7: return (lowerEdge ? -0.609 : -0.522); break;
  case  -6: return (lowerEdge ? -0.522 : -0.435); break;
  case  -5: return (lowerEdge ? -0.435 : -0.348); break;
  case  -4: return (lowerEdge ? -0.348 : -0.261); break;
  case  -3: return (lowerEdge ? -0.261 : -0.174); break;
  case  -2: return (lowerEdge ? -0.174 : -0.087); break;
  case  -1: return (lowerEdge ? -0.087 :      0); break;
  case  +1: return (lowerEdge ?      0 :  0.087); break;
  case  +2: return (lowerEdge ?  0.087 :  0.174); break;
  case  +3: return (lowerEdge ?  0.174 :  0.261); break;
  case  +4: return (lowerEdge ?  0.261 :  0.348); break;
  case  +5: return (lowerEdge ?  0.348 :  0.435); break;
  case  +6: return (lowerEdge ?  0.435 :  0.522); break;
  case  +7: return (lowerEdge ?  0.522 :  0.609); break;
  case  +8: return (lowerEdge ?  0.609 :  0.696); break;
  case  +9: return (lowerEdge ?  0.696 :  0.783); break;
  case +10: return (lowerEdge ?  0.783 :  0.879); break;
  case +11: return (lowerEdge ?  0.879 :  0.957); break;
  case +12: return (lowerEdge ?  0.957 :  1.044); break;
  case +13: return (lowerEdge ?  1.044 :  1.131); break;
  case +14: return (lowerEdge ?  1.131 :  1.218); break;
  case +15: return (lowerEdge ?  1.218 :  1.305); break;
  case +16: return (lowerEdge ?  1.305 :  1.392); break;
  case +17: return (lowerEdge ?  1.392 :  1.479); break;
  case +18: return (lowerEdge ?  1.479 :  1.566); break;
  case +19: return (lowerEdge ?  1.566 :  1.653); break;
  case +20: return (lowerEdge ?  1.653 :   1.74); break;
  case +21: return (lowerEdge ?   1.74 :   1.83); break;
  case +22: return (lowerEdge ?   1.83 :   1.93); break;
  case +23: return (lowerEdge ?   1.93 :  2.043); break;
  case +24: return (lowerEdge ?  2.043 :  2.172); break;
  case +25: return (lowerEdge ?  2.172 :  2.322); break;
  case +26: return (lowerEdge ?  2.322 :    2.5); break;
  case +27: return (lowerEdge ?    2.5 :   2.65); break;
  case +28: return (lowerEdge ?   2.65 :  2.853); break;
  case +29: return (lowerEdge ?  2.853 :  2.964); break;
  case +30: return (lowerEdge ?  2.964 :  3.139); break;
  case +31: return (lowerEdge ?  3.139 :  3.314); break;
  case +32: return (lowerEdge ?  3.314 :  3.489); break;
  case +33: return (lowerEdge ?  3.489 :  3.664); break;
  case +34: return (lowerEdge ?  3.664 :  3.839); break;
  case +35: return (lowerEdge ?  3.839 :  4.013); break;
  case +36: return (lowerEdge ?  4.013 :  4.191); break;
  case +37: return (lowerEdge ?  4.191 :  4.363); break;
  case +38: return (lowerEdge ?  4.363 :  4.538); break;
  case +39: return (lowerEdge ?  4.538 :  4.716); break;
  case +40: return (lowerEdge ?  4.716 :  4.889); break;
  case +41: return (lowerEdge ?  4.889 :  5.191); break;
    //something went wrong;
  default : return -1; break;
  }
}



// -----------------------------------------------------------------
int TParameters::GetPhiBin(int phi_id, int phigranu) const
//This function knows the number of wanted phi-bins and returns 
//in which phi-bin the tower with eta-ID "phi_id" is located.
{
  assert(phi_id >0);
  assert(phi_id <= 72);
  return (phi_id-1)*phigranu/phi_ntwr;
}



// -----------------------------------------------------------------
void TParameters::print() const
{
  std::cout  << p->name() << " resulting in:\n "
	     << eta_granularity << " x " << phi_granularity << " tower bins with " 
	     << GetNumberOfTowerParametersPerBin() << " free parameters each, or " 
	     << GetNumberOfTowerParameters() << " in total, and\n "
	     << eta_granularity_jet << " x " << phi_granularity_jet << " JES bins with " 
	     << GetNumberOfJetParametersPerBin() << " free parameters each, or " 
	     << GetNumberOfJetParameters() << " in total \n "
	     << eta_granularity_track << " x " << phi_granularity_track << " track bins with " 
	     << GetNumberOfTrackParametersPerBin() << " free parameters each, or " 
	     << GetNumberOfTrackParameters() << " in total \n "
	     << "and " << GetNumberOfGlobalJetParameters() << " global jet parameters\n"; 
}



// -----------------------------------------------------------------
void TParameters::writeCalibrationTxt(const char* name)
{
  cout << "Writing calibration to file '" << name << "'" << endl;

  ofstream file(name, ofstream::binary);
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned iphi=1; iphi<=1; ++iphi){ // No phi binning
      int towerIdx = GetBin(GetEtaBin(ieta),GetPhiBin(iphi));
      int jetIdx = GetJetBin(GetJetEtaBin(ieta),GetJetPhiBin(iphi));
      int trackIdx = GetTrackBin(GetTrackEtaBin(ieta),GetTrackPhiBin(iphi));
      if(towerIdx<0 || jetIdx<0 || trackIdx<0) continue;
      // write: lower eta | upper eta | nparameters, for
      // each eta id of the tower and n times for n phi bins
      file << std::setw(10) << etaLowerEdge(ieta) 
	   << std::setw(10) << etaUpperEdge(ieta)  
	   << std::setw(10) << 2 + p->nTowerPars()+p->nJetPars()+p->nTrackPars()+p->nGlobalJetPars();
      // Dummy: pt range of validity
      file << std::setw(8) << std::setprecision(4) << 4;
      file << std::setw(8) << std::setprecision(4) << 2000;
      // write: each tower parameter
      for(unsigned itower=0; itower<p->nTowerPars(); ++itower){
	file << std::setw(12) << std::setprecision(4) << k[towerIdx*p->nTowerPars()+itower];
      }
      // write: each jet parameter
      for(unsigned ijet=0; ijet<p->nJetPars(); ++ijet){
	file << std::setw(12) << std::setprecision(4) << k[GetNumberOfTowerParameters()+jetIdx*p->nJetPars()+ijet];
      }
      // write: each track parameter
      for(unsigned itrack=0; itrack<p->nTrackPars(); ++itrack){
	file << std::setw(12) << std::setprecision(4) << k[GetNumberOfTowerParameters()+GetNumberOfJetParameters()+trackIdx*p->nTrackPars()+itrack];
      }
      for(unsigned igjet=0; igjet<p->nGlobalJetPars(); ++igjet){
	file << std::setw(12) << std::setprecision(4) << k[GetNumberOfTowerParameters()+GetNumberOfJetParameters()+GetNumberOfTrackParameters()+igjet];
      }
      // complete line
      file << std::endl;
    }
  }
  file.close();
}



// -----------------------------------------------------------------
void TParameters::writeCalibrationCfi(const char* name)
{
  cout << "Writing calibration to file '" << name << "'" << endl;
  ofstream file(name, ofstream::binary);
  
  time_t rawtime = time(0);
  struct tm * timeinfo;
  char buffer [80];
  timeinfo = localtime ( &rawtime );
  strftime (buffer,80," # Hamburg Calorimeter Calibration Tool, created %c",timeinfo);
  struct passwd* pw = getpwuid(getuid());	
  
  //double aux[10000], fsum;
  //int npar = GetNumberOfParameters(), iflag=0;
  //fitfunction(npar, aux, fsum, k, iflag);  

  file << buffer << " by " << pw->pw_name << "." << endl 
       << " block CalibParameters = {" << endl
       << "    untracked string  Parametrization    = " << '\"' << p->name() << '\"' <<  endl
       << "    untracked int32  NTowerParamsPerBin = " << GetNumberOfTowerParametersPerBin() << endl
       << "    untracked int32  NJetParamsPerBin   = " << GetNumberOfJetParametersPerBin() << endl
       << "    untracked int32  NTrackParamsPerBin   = " << GetNumberOfTrackParametersPerBin() << endl
       << "    untracked int32  NGlobalJetParams   = " << GetNumberOfGlobalJetParameters() << endl
       << "    untracked int32  NEtaBins           = " << eta_granularity << endl
       << "    untracked int32  NPhiBins           = " << phi_granularity << endl
       << "    untracked bool   EtaSymmetryUsed    = " << eta_symmetry << endl
       << "    untracked double FitChi2            = " << GetFitChi2() << endl
       << " }";
  //--------------------------------------------------------------------
  file << endl
       << " block TowerCalibConstants = {" << endl
       << "    untracked vint32 TowMapEta       = { ";
  //1. ieta
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
      if (ieta!=-41 || iphi!=1)
        file << ", " << ieta;
      else
        file << ieta;	
    }
  }
  file << " }" << endl << "    untracked vint32 TowMapPhi       = { ";

  //2. iphi
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
      if (ieta!= -41 || iphi!=1)
        file << ", " << iphi;
      else
        file << iphi;
    }
  }
  file << " }"<< endl;
  file << endl;

  //3. tower calibration constants
  file << "    untracked  int32  TowerParam  = " << GetNumberOfTowerParametersPerBin() << endl;
  for (unsigned int n=0; n<p->nTowerPars(); ++n) {
    file << "    untracked vdouble TowerParam"<< n <<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
        int index = GetBin(GetEtaBin(ieta),GetPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          file << ", " << k[index*p->nTowerPars()+n];
	else
          file << k[index*p->nTowerPars()+n];;
      }
    }
    file << " }" << endl; 
  }

  //4. calibration constants errors
  for (unsigned int n=0; n<p->nTowerPars(); ++n) {
    file << "    untracked vdouble TowerError"<<n<<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
        int index = GetBin(GetEtaBin(ieta),GetPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          file << ", " << parErrors_[index*p->nTowerPars()+n];
	else
          file << parErrors_[index*p->nTowerPars()+n];;
      }
    }
    file << " }" << endl; 
  }
  file << " }" << endl; 
  //--------------------------------------------------------------------
  file << endl
       << " block JetCalibConstants = {" << endl
       << "    InputTag Jets    = MyFavoriteJetAlgorithm" << endl
       << "    string CalibJets = \"\" " << endl
       << endl
       << "    untracked vint32 JetMapEta     = { ";
  
  //5. ieta
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
      if (ieta!=-41 || iphi!=1)
        file << ", " << ieta;
      else
        file << ieta;	
    }
  }
  file << " }" << endl << "    untracked vint32 JetMapPhi     = { ";
  
  //6. iphi
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
      if (ieta!= -41 || iphi!=1)
        file << ", " << iphi;
      else
        file << iphi;
    }
  }
  file << " }" << endl;
  file << endl;

  //7. jet calibration constants
  file << "    untracked  int32  JetParam  = " << GetNumberOfJetParametersPerBin() << endl;
  for (unsigned int n=0; n<p->nJetPars(); ++n) {
    file << "    untracked vdouble JetParam"<<n<<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
	int index = GetJetBin(GetJetEtaBin(ieta),GetJetPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          file << ", " << k[GetNumberOfTowerParameters()+index*p->nJetPars()+n];
	else
          file << k[GetNumberOfTowerParameters()+index*p->nJetPars()+n];
      }
    }
    file << " }" << endl; 
  }
  //8. calibration constants errors
  for (unsigned int n=0; n<p->nJetPars(); ++n) {
    file << "    untracked vdouble JetError"<<n<<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
	int index = GetBin(GetJetEtaBin(ieta),GetJetPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          file << ", " << parErrors_[GetNumberOfTowerParameters() + index*p->nJetPars()+n];
	else
          file << parErrors_[GetNumberOfTowerParameters() + index*p->nJetPars()+n];
      }
    }
    file << " }" << endl; 
  }
  //8a. global jet constants
  file << "    untracked vdouble GlobalJetParams = { ";
  for(unsigned int n = 0 ; n < p->nGlobalJetPars() ; ++n) {
    if(n != 0) file << ", " << k[GetNumberOfTowerParameters()+GetNumberOfJetParameters()+GetNumberOfTrackParameters()+n];
    else file << k[GetNumberOfTowerParameters()+GetNumberOfJetParameters()+GetNumberOfTrackParameters()+n];
  } 
  file << " }" << endl; 
  file << "    untracked vdouble GlobalJetErrors = { ";
  for(unsigned int n = 0 ; n < p->nGlobalJetPars() ; ++n) {
    if(n != 0) file << ", " << parErrors_[GetNumberOfTowerParameters()+GetNumberOfJetParameters()+GetNumberOfTrackParameters()+n];
    else file << parErrors_[GetNumberOfTowerParameters()+GetNumberOfJetParameters()+GetNumberOfTrackParameters()+n];
  } 
  file << " }" << endl; 
  
  file << " }" << endl; 
  //--------------------------------------------------------------------
  file << endl
       << " block TrackCalibConstants = {" << endl
       << "    InputTag Tracks    = MyFavoriteTrackAlgorithm" << endl
       << "    string CalibTracks = \"\" " << endl
       << endl
       << "    untracked vint32 TrackMapEta     = { ";
  
  //9. ieta
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
      if (ieta!=-41 || iphi!=1)
        file << ", " << ieta;
      else
        file << ieta;	
    }
  }
  file << " }" << endl << "    untracked vint32 TrackMapPhi     = { ";
  
  //10. iphi
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
      if (ieta!= -41 || iphi!=1)
        file << ", " << iphi;
      else
        file << iphi;
    }
  }
  file << " }" << endl;
  file << endl;

  //11. track calibration constants
  file << "    untracked  int32  TrackParam  = " << GetNumberOfTrackParametersPerBin() << endl;
  for (unsigned int n=0; n<p->nTrackPars(); ++n) {
    file << "    untracked vdouble TrackParam"<<n<<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
	int index = GetTrackBin(GetTrackEtaBin(ieta),GetTrackPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          file << ", " << k[GetNumberOfTowerParameters()+GetNumberOfJetParameters()+index*p->nTrackPars()+n];
	else
          file << k[GetNumberOfTowerParameters()+GetNumberOfJetParameters()+index*p->nTrackPars()+n];
      }
    }
    file << " }" << endl; 
  }
  //12. calibration constants errors
  for (unsigned int n=0; n<p->nTrackPars(); ++n) {
    file << "    untracked vdouble TrackError"<<n<<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
	int index = GetBin(GetTrackEtaBin(ieta),GetTrackPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          file << ", " << parErrors_[GetNumberOfTowerParameters() + GetNumberOfJetParameters() + index*p->nTrackPars()+n];
	else
          file << parErrors_[GetNumberOfTowerParameters() + GetNumberOfJetParameters() + index*p->nTrackPars()+n];
      }
    }
    file << " }" << endl; 
  }
  file << " }" << endl;  file << endl
       << " block TrackCalibConstants = {" << endl
       << "    InputTag Tracks    = MyFavoriteTrackAlgorithm" << endl
       << "    string CalibTracks = \"\" " << endl
       << endl
       << "    untracked vint32 TrackMapEta     = { ";
  
  //9. ieta
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
      if (ieta!=-41 || iphi!=1)
        file << ", " << ieta;
      else
        file << ieta;	
    }
  }
  file << " }" << endl << "    untracked vint32 TrackMapPhi     = { ";
  
  //10. iphi
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
      if (ieta!= -41 || iphi!=1)
        file << ", " << iphi;
      else
        file << iphi;
    }
  }
  file << " }" << endl;
  file << endl;

  //11. track calibration constants
  file << "    untracked  int32  TrackParam  = " << GetNumberOfTrackParametersPerBin() << endl;
  for (unsigned int n=0; n<p->nTrackPars(); ++n) {
    file << "    untracked vdouble TrackParam"<<n<<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
	int index = GetTrackBin(GetTrackEtaBin(ieta),GetTrackPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          file << ", " << k[GetNumberOfTowerParameters()+GetNumberOfJetParameters()+index*p->nTrackPars()+n];
	else
          file << k[GetNumberOfTowerParameters()+GetNumberOfJetParameters()+index*p->nTrackPars()+n];
      }
    }
    file << " }" << endl; 
  }
  //12. calibration constants errors
  for (unsigned int n=0; n<p->nTrackPars(); ++n) {
    file << "    untracked vdouble TrackError"<<n<<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
	int index = GetBin(GetTrackEtaBin(ieta),GetTrackPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          file << ", " << parErrors_[GetNumberOfTowerParameters() + GetNumberOfJetParameters() + index*p->nTrackPars()+n];
	else
          file << parErrors_[GetNumberOfTowerParameters() + GetNumberOfJetParameters() + index*p->nTrackPars()+n];
      }
    }
    file << " }" << endl; 
  }
  file << " }" << endl;  
  // complete line
  file << std::endl;
  file.close();
}



//-----------------------------------------------------
void TParameters::writeCalibrationTex(const char* name, const ConfigFile& config)
{
  cout << "Writing calibration to file '" << name << "'" << endl;

  // Getting fitted jet parameter values from TParameters
  std::vector<double> pJetFit(GetNumberOfJetParameters());
  std::vector<double> pJetError(GetNumberOfJetParameters());
  std::vector<double> pJetGCorr(GetNumberOfJetParameters());
  for(int i = 0; i < GetNumberOfJetParameters(); i++) {
    int jetbin = 0;
    pJetFit[i] = GetJetParRef(jetbin)[i];
    pJetError[i] = GetJetParErrorRef(jetbin)[i];
    pJetGCorr[i] = GetJetParGlobalCorrCoeffRef(jetbin)[i];
  }
  
  // Getting scales from config file
  std::vector<double> pJetScale = bag_of<double>(config.read<string>("jet parameter scales","")); 

  // Getting fitted global parameter values from TParameters
  std::vector<double> pGlobalJetFit(GetNumberOfGlobalJetParameters());
  std::vector<double> pGlobalJetParError(GetNumberOfGlobalJetParameters());
  std::vector<double> pGlobalJetParGCorr(GetNumberOfGlobalJetParameters());
  for(int i = 0; i < GetNumberOfGlobalJetParameters(); i++) {
    pGlobalJetFit[i] = GetGlobalJetParRef()[i];
    pGlobalJetParError[i] = GetGlobalJetParErrorRef()[i];
    pGlobalJetParGCorr[i] = GetGlobalJetParGlobalCorrCoeffRef()[i];
  }

  // Write to tex file
  ofstream outfile(name, ofstream::binary);
  if(outfile.is_open()) {
    outfile << "\\documentclass{article}\n";
    outfile << "\\usepackage{color}\n";
    outfile << "\\definecolor{gray}{rgb}{0.5,0.5,0.5}\n";

    outfile << "\\begin{document}\n";

    // Parametriaztion
    outfile << "Parametrization class: \\texttt{";
    outfile << config.read<string>("Parametrization Class","");
    outfile << "}\n";

    // BFGS parameters
    outfile << "\\begin{flushleft}\n\\begin{tabular}{lcl}\n";
    outfile << texTabularLine<double>(config,"BFGS derivative step");
    outfile << texTabularLine<int>(config,"BFGS mvec");
    outfile << texTabularLine<int>(config,"BFGS niter");
    outfile << texTabularLine<double>(config,"BFGS eps");
    outfile << texTabularLine<double>(config,"BFGS 1st wolfe parameter");
    outfile << texTabularLine<double>(config,"BFGS 2nd wolfe parameter");
    outfile << "\\end{tabular}\\end{flushleft}\n";

    // Event selection cuts
    outfile << "\\begin{flushleft}\n\\begin{tabular}{lcl}\n";
    outfile << texTabularLine<double>(config,"Et genJet min");
    outfile << texTabularLine<double>(config,"Et genJet max");
    outfile << texTabularLine<double>(config,"DeltaR cut on jet matching");
    outfile << texTabularLine<int>(config,"Et cut on jet");
    outfile << texTabularLine<int>(config,"Eta cut on jet");
    outfile << texTabularLine<double>(config,"Min had fraction");
    outfile << texTabularLine<double>(config,"Max had fraction");
    outfile << texTabularLine<double>(config,"Relative n+1 Jet Et Cut");
    outfile << "\\end{tabular}\\end{flushleft}\n";

    // Dijet integration parameters
    outfile << "\\begin{flushleft}\n\\begin{tabular}{lcl}\n";
    outfile << texTabularLine<int>(config,"DiJet integration number of iterations");
    outfile << texTabularLine<double>(config,"DiJet integration epsilon");
    outfile << texTabularLine<double>(config,"DiJet integration min");
    outfile << texTabularLine<double>(config,"DiJet integration max");
    outfile << texTabularLine<double>(config,"DiJet integration pt bin edges");
    outfile << "\\end{tabular}\\end{flushleft}\n";

    // Start and fitted parameters
    outfile << "\\begin{center}\n";
    outfile << "\\begin{tabular}{ccccc}\n";
    outfile << "\\hline\n\\hline\n";
    outfile << "Index & Scale & Start value & Fitted value & Global correlation \\\\ \n\\hline \n";
    for(unsigned int i = 0; i < pJetScale.size() && i < jet_start_values.size() && i < pJetFit.size(); i++) {
      if( isFixedPar(i) ) {
	outfile << "\\textcolor{gray}{$" << i << "$} & \\textcolor{gray}{$ ";
	outfile << pJetScale.at(i) << " $} & \\textcolor{gray}{$ ";
 	outfile << jet_start_values.at(i) << "$ } & \\textcolor{gray}{$ ";
 	outfile << pJetFit[i] << " \\pm ";
 	outfile << pJetError[i] << " $} & \\textcolor{gray}{$ ";
	outfile << pJetGCorr[i] << " $} \\\\ \n";
      } else {
	outfile << "$" << i << "$ & $";
	outfile << pJetScale.at(i) << "$ & $";
	outfile << jet_start_values.at(i) << "$ & $";
	outfile << pJetFit[i] << " \\pm ";
	outfile << pJetError[i] << "$ & $";
	outfile << pJetGCorr[i] << "$ \\\\ \n";
      }
    }
    for(unsigned int i = 0; i < pGlobalJetFit.size(); i++) {
      if( i == 0 ) {
	outfile << "\\hline\n";
      }
      outfile << "$" << i << "$ & $";
      outfile << 1. << "$ & $";
      outfile << global_jet_start_values.at(i) << "$ & $";
      outfile << pGlobalJetFit[i] << " \\pm ";
      outfile << pGlobalJetParError[i] << "$ & $";
      outfile << pGlobalJetParGCorr[i] << "$ \\\\ \n";
    }
    outfile << "\\hline\n\\hline\n";
    outfile << "\\end{tabular}\n";
    outfile << "\\end{center}\n";
    outfile << "\\end{document}\n";
    outfile.close();
  }
}



// -----------------------------------------------------------------
int TParameters::GetTrackEffBin(double pt, double eta)
{
  int bin, etabin, ptbin;
  etabin = (int)(fabs(eta)/0.2);
  if (etabin > 12) etabin = 12;
  if (pt < 5) ptbin = (int)(pt); //bin 0-4
  else{
    if(pt < 30) ptbin = (int)(5+(pt-5)/5);  //bin 5-9
    else{
      if(pt < 40) ptbin = 10;
      else{
	if(pt<50) ptbin = 11;
	else ptbin = 12;
      }
    }
  }
  bin = 13 * etabin + ptbin;
  return bin;
}

Function TParameters::tower_function(int etaid, int phiid) {
  int id = GetBin(GetEtaBin(etaid),GetPhiBin(phiid));
  if (id <0) { 
    std::cerr<<"WARNING: TParameters::tower_function::index = " << id << endl; 
    exit(-2);  
  }
  int parIndex = id*GetNumberOfTowerParametersPerBin();
  return Function(&Parametrization::correctedTowerEt,0,parIndex,
		  GetNumberOfTowerParametersPerBin(),
		  GetTowerParRef(id),p);
}

Function TParameters::jet_function(int etaid, int phiid) {
  int id = GetJetBin(GetJetEtaBin(etaid),GetJetPhiBin(phiid));
  if (id <0) { 
    std::cerr<<"WARNING: TParameters::jet_function::index = " << id << endl; 
    exit(-2);  
  }
  int parIndex = id * GetNumberOfJetParametersPerBin() + GetNumberOfTowerParameters();
  return Function(&Parametrization::correctedJetEt,
		  p->hasInvertedCorrection() ? &Parametrization::inverseJetCorrection : 0,
		  parIndex,GetNumberOfJetParametersPerBin(),
		  GetJetParRef(id),p);
}

Function TParameters::track_function(int etaid, int phiid) {
  int id = (etaid == 0) && (phiid == 0) ? 0: GetTrackBin(GetTrackEtaBin(etaid),GetTrackPhiBin(phiid));
  if (id <0) { 
    std::cerr<<"WARNING: TParameters::track_function::index = " << id << endl; 
    exit(-2);  
  }
  int parIndex = id * GetNumberOfTrackParametersPerBin() +
    GetNumberOfTowerParameters() + GetNumberOfJetParameters();
  return Function(&Parametrization::GetExpectedResponse,0,
		  parIndex,GetNumberOfTrackParametersPerBin(),
		  GetTrackParRef(id),p);
}

Function TParameters::global_jet_function() {
  int parIndex = GetNumberOfTowerParameters()+GetNumberOfJetParameters()+GetNumberOfTrackParameters();
  return Function(&Parametrization::correctedGlobalJetEt,0,
		  parIndex,GetNumberOfGlobalJetParameters(),
		  GetGlobalJetParRef(),p);
}

SmearFunction TParameters::resolutionFitPDF(int etaid, int phiid) {
  int id = GetJetBin(GetJetEtaBin(etaid),GetJetPhiBin(phiid));
  if (id <0) { 
    std::cerr<<"WARNING: TParameters::resolutionFitPDF::index = " << id << endl; 
    exit(-2);  
  }
  int jetIdx = id * GetNumberOfJetParametersPerBin() + GetNumberOfTowerParameters();
  return SmearFunction(&Parametrization::pdfPtMeas,
		       &Parametrization::pdfPtTrue,
		       &Parametrization::pdfPtTrueError,
		       &Parametrization::pdfResponse,
		       &Parametrization::pdfResponseError,
		       &Parametrization::pdfDijetAsym,
		       jetIdx,GetNumberOfJetParametersPerBin(),GetJetParRef(id),
		       findParStatus(jetIdx,GetNumberOfJetParameters()),
		       findCovIndices(jetIdx,GetNumberOfJetParameters()),
		       GetCovCoeff(),p);
}

double TParameters::toy_tower_error_parametrization(const double *x, const Measurement *xorig, double errorig)
{
  double hadet = x[0] - xorig->EMF - xorig->OutF;
  if(hadet < 0.001) hadet = 0.001;
  double hade = hadet * xorig->E / xorig->pt; 
  //std::cout << "had Et:" << hadet << " , " << "had E:" << hade << '\n';
  
  double a = 4.44;
  double b = 1.11;
  double c = 0.03;
  
  double var = a*a/hade/hade + b*b/hade + c*c;
  //truncate variance accordingly
  double truncvar = - sqrt(var) * exp(-0.5/var) * sqrt(2/M_PI) + var * TMath::Erf(1/(sqrt(2 * var)));
  return sqrt(truncvar) * hadet;
}

double TParameters::toy_jet_error_parametrization(const double *x, const Measurement *xorig, double errorig)
{
  double a = 4.44;
  double b = 1.11;
  double c = 0.03;
  
  //return sqrt(a*a/x[0]/x[0] + b*b/x[0] + c*c)*x[0];
  
  double e   = x[0] * xorig->E / xorig->pt;
  double var = a*a/e/e + b*b/e + c*c;
  //truncate variance accordingly
  double truncvar = - sqrt(var) * exp(-0.5/var) * sqrt(2/M_PI) + var * TMath::Erf(1/(sqrt(2 * var)));
  return sqrt(truncvar) * x[0];
}

//!  \brief Return one line of LaTeX tabular containing the
//!         name and value of a given parameter from config file
//!
//!  The line is the following: "\texttt{<fieldname>} & = & <value> \\"
//!
//!  \param config Configfile
//!  \param fieldname Name of parameter in config file
//!  \return The line for the LaTeX tabular
// --------------------------------------------------
template<class T> std::string TParameters::texTabularLine(const ConfigFile& config, const std::string& fieldname) const {
  std::stringstream line;
  line << "\\texttt{" << fieldname << "} & = & $";
  line << config.read<T>(fieldname.c_str(),-1) << "$ \\\\ \n";
  
  return line.str();
}


std::vector<int> TParameters::findCovIndices(int firstPar, int nPar) const {
  // Dimension of the submatrix
  int nCov = (nPar*nPar+nPar)/2;
  std::vector<int> indices(nCov);

  int idx = ((firstPar+1)*(firstPar+1)+firstPar+1)/2 - 1;
  int rowIdx = firstPar;
  int i = 0;
  while( i < nCov ) {
    indices[i] = idx;

    int max = ((rowIdx+1)*(rowIdx+1)+rowIdx+1)/2 - 1;
    if( idx == max ) {
      rowIdx++;
      idx += firstPar;
    }
    i++;
    idx++;
  }

  int maxNCov = GetNumberOfCovCoeffs();
  assert( indices.back() < maxNCov );

  return indices;
}


std::vector<bool> TParameters::findParStatus(int firstPar, int nPar) const {
  std::vector<bool> isFixed(nPar);
  for(int i = 0; i < nPar; i++) {
    isFixed[i] = isFixedPar(firstPar+i);
  }

  return isFixed;
}

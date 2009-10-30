//  $Id: caliber.cc,v 1.95 2009/10/30 08:14:24 mschrode Exp $

#include "caliber.h"

#include <algorithm>
#include <iomanip>

//Boost
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
boost::mutex io_mutex;
// User
#include "ConfigFile.h"
#include "Parameters.h"
#include "ControlPlots.h"
#include "ControlPlotsJetSmearing.h"
#include "CalibMath.h"
#include "external.h"
#include "ToyMC.h"
#include "PhotonJetReader.h"
#include "DiJetReader.h"
#include "TriJetReader.h"
#include "ZJetReader.h"
#include "TopReader.h"
#include "TrackClusterReader.h"
#include "ParameterLimitsReader.h"
#include "TowerConstraintsReader.h"
#include "JetConstraintsReader.h"
#include "EventProcessor.h"
#include "EventWeightProcessor.h"
#include "Jet.h"
#include "JetTruthEvent.h"

using namespace std;

typedef std::vector<TData*>::iterator DataIter;
typedef std::vector<TData*>::const_iterator DataConstIter;



//!  \brief Outlier Rejection
// -----------------------------------------------------------------
struct OutlierRejection {
  OutlierRejection(double cut):_cut(cut){};
  bool operator()(TData *d){
    if(d->GetType()==typeTowerConstraint) return true;
    return (d->chi2()/d->GetWeight())<_cut;
  }
  double _cut;
};



// -----------------------------------------------------------------
class ComputeThread {
private:
  int npar;
  double chi2;
  double * td1;
  double * td2;
  double *parorig, *mypar;
  double epsilon;
  std::vector<TData*> data;
  bool data_changed;
  struct calc_chi2_on
  {
  private:
    ComputeThread *parent;
  public:
    calc_chi2_on(ComputeThread *parent) : parent(parent) {}
    void operator()()
    {
      //      {
      // 	boost::mutex::scoped_lock lock(io_mutex);
      // 	std::cout << "start Thread for " << parent << std::endl; 
      //       }   
      if(parent->data_changed) {
	for (DataIter it=parent->data.begin() ; it!= parent->data.end() ; ++it)
	  (*it)->ChangeParAddress(parent->parorig,parent->mypar); 
	parent->data_changed = false;
      }
      for (int param=0; param< parent->npar ; ++param) {
	parent->td1[param]= 0.0;
	parent->td2[param]= 0.0;
	parent->mypar[param] = parent->parorig[param];
      }
      parent->chi2 =0.0;   
      for (DataIter it=parent->data.begin() ; it!= parent->data.end() ; ++it) { 
	//boost::mutex::scoped_lock lock(io_mutex);
	parent->chi2 += (*it)->chi2_fast(parent->td1, parent->td2, parent->epsilon);
      }
    }
  };
  boost::thread *thread;
  friend class calc_chi2_on;
public:
  ComputeThread(int npar,double *par, double epsilon) 
    : npar(npar), td1(new double[npar]), td2(new double[npar]), parorig(par),
      mypar(new double[npar]), epsilon(epsilon), data_changed(false) {
    //std::cout << "threads par array:" << mypar << '\n';
  }
  ~ComputeThread() {
    ClearData();
    delete [] td1;
    delete [] td2;
    delete [] mypar;
  }
  void AddData(TData* d) { 
    //d->ChangeParAddress(parorig, mypar);
    data_changed = true;
    data.push_back(d);
  }
  void ClearData() {   
    for (DataIter it= data.begin() ; it!= data.end() ; ++it)  
      (*it)->ChangeParAddress(mypar,parorig);
    data.clear();
  }
  void Start() { thread = new boost::thread(calc_chi2_on(this)); }
  bool IsDone() { thread->join(); delete thread; return true;}
  void SyncParameters() {
    for (int param=0; param< npar ; ++param) mypar[param] = parorig[param];
  }
  double Chi2() const { return chi2;}
  double TempDeriv1(int i) const { return td1[i];}
  double TempDeriv2(int i) const { return td2[i];}
};



//--------------------------------------------------------------------------------------------
void TCaliber::run()
{
  if (fitMethod_!=3){

    time_t start = time(0);
    
    //calls FlattenSpectra and BalanceSpectra if enabled in config
    EventProcessor ep(configFile_,par_);
    ep.process(data_);

    // Apply event weights if enabled
    EventWeightProcessor ewp(configFile_,par_);
    ewp.process(data_);

    if (fitMethod_==1) {
      run_Lvmini();
      time_t end = time(0);
      cout << "Done, fitted " << par_->GetNumberOfParameters() << " parameters in " << difftime(end,start) << " sec." << endl;
    }
  } 
  //Dummy Configuration: Nothing to be done, start-values are written to file
}



// -----------------------------------------------------------------
void TCaliber::run_Lvmini()
{ 
  int naux = 3000000, iret=0;
  
  int npar = par_->GetNumberOfParameters();

  naux = lvmdim_(npar,mvec_);
  cout<<"array of size "<<naux<<" needed."<<endl;

  double* aux = new double[naux], fsum = 0;

  double *temp_derivative1 = new double[npar];
  double *temp_derivative2 = new double[npar];

  cout << "\nFitting " << npar << " parameters; \n";
  par_->print();
  cout << " with LVMINI.\n" << "Using " << data_.size() << " total events and ";
  cout << nThreads_ << " threads.\n";

  // Fixed pars
  if( fixedJetPars_.size() > 0 ) cout << "Fixed jet parameters:\n";
  for(unsigned int i = 0; i < fixedJetPars_.size(); i++) {
    int idx = fixedJetPars_.at(i);
    cout << "  " << idx+1 << ": " << par_->GetPars()[idx] << endl;
  }
  if( fixedGlobalJetPars_.size() > 0 ) cout << "Fixed global jet parameters:\n";
  for(unsigned int i = 0; i < fixedGlobalJetPars_.size(); i++) {
    int idx = fixedGlobalJetPars_.at(i);
    cout << "  " << idx+1 << ": " << par_->GetPars()[idx] << endl;
  }
  
  ComputeThread *t[nThreads_];
  for (int ithreads=0; ithreads<nThreads_; ++ithreads){
    t[ithreads] = new ComputeThread(npar, par_->GetPars(),derivStep_);
  }

  lvmeps_(data_.size()*eps_,wlf1_,wlf2_);
  lvmeps_(eps_,wlf1_,wlf2_);

  //Set errors per default to 0 //@@ doesn't seem to work...
  int error_index=2;
  error_index = lvmind_(error_index);
  par_->FillErrors(aux+error_index);

  for( unsigned int loop = 0; loop < residualScalingScheme_.size() ; ++loop ) {
    cout<<"Updating Di-Jet Errors"<<endl;
    for(DataIter it = data_.begin()  ; it < data_.end() ; ++it) {
      (*it)->updateError();
    }

    if( par_->needsUpdate() ) par_->update();

    // Setting function to scale residuals in chi2 calculation
    cout << loop+1 << flush;
    if(  loop+1 == 1  ) cout << "st" << flush;
    else if(  loop+1 == 2  ) cout << "nd" << flush;
    else if(  loop+1 == 3  ) cout << "rd" << flush;
    else cout << "th" << flush;
    cout << " of " << residualScalingScheme_.size() <<" iteration(s): " << flush;
    if(  residualScalingScheme_.at(loop) == 0  ) {
	TData::ScaleResidual = &TData::ScaleNone;	
	cout << "no scaling of residuals." << endl;

	cout << "Rejecting outliers " << flush;
	DataIter beg = partition(data_.begin(), data_.end(), OutlierRejection(outlierChi2Cut_));
	for(DataIter i = beg ; i != data_.end() ; ++i) {
	  delete *i;
	}
	data_.erase(beg,data_.end());
	cout << "and using " << data_.size() << " events." << endl;
      }
    else if(  residualScalingScheme_.at(loop) == 1  ) {
	TData::ScaleResidual = &TData::ScaleCauchy;	
	cout << "scaling of residuals with Cauchy-Function." << endl;
      }
    else if(  residualScalingScheme_.at(loop) == 2  ) {
	TData::ScaleResidual = &TData::ScaleHuber;	
	cout << "scaling of residuals with Huber-Function." << endl;
      }
    else if(  residualScalingScheme_.at(loop) == 3  ) {
      TData::ScaleResidual = &TData::ScaleTukey;	
      cout << "scaling of residuals a la Tukey." << endl;
    }
    else {
      cerr << "ERROR: " << residualScalingScheme_.at(loop) << " is not a valid scheme for resdiual scaling! Breaking iteration!" << endl;
      break;
    }
    if(lvmdim_(npar,mvec_) > naux)
      cout<<"Aux field too small. "<<lvmdim_(npar,mvec_)<<" enntires needed."<<endl;
    if (npar>0) npar*=-1; //Show output
    //initialization
    lvmini_( npar, mvec_, nIter_, aux);
    npar=std::abs(npar);
    
    int n = 0;
    
    for(DataIter it = data_.begin()  ; it < data_.end() ; ++it) {
      t[n]->AddData(*it);
      n++;
      if(n == nThreads_) n = 0;
    }
    
    do {
      //set storage for temporary derivative storage to zero
      for (int param=0; param< npar ; ++param) {
	temp_derivative1[param]=0.0;
	temp_derivative2[param]=0.0;
      } 
      //set local parameters to global value
      for( std::vector<int>::const_iterator iter = globalJetPars_.begin();
	   iter != globalJetPars_.end() ; ++ iter) {
	double val = par_->GetPars()[*iter];
	for(int id = *iter + par_->GetNumberOfJetParametersPerBin(); 
	    id < par_->GetNumberOfJetParameters() ; 
	    id += par_->GetNumberOfJetParametersPerBin()) {
	  par_->GetPars()[id] = val;
	}
      }
      fsum = 0;
      for (int ithreads=0; ithreads<nThreads_; ++ithreads) t[ithreads]->Start();
      
      for (int ithreads=0; ithreads<nThreads_; ++ithreads){
	if(t[ithreads]->IsDone()) {
	  fsum += t[ithreads]->Chi2();
	  for (int param=0 ; param < npar ; ++param) {
	    temp_derivative1[param] += t[ithreads]->TempDeriv1(param);
	    temp_derivative2[param] += t[ithreads]->TempDeriv2(param);
	  }
	}
      }
      //sum up derivative results for global par
      for( std::vector<int>::const_iterator iter = globalJetPars_.begin();
	   iter != globalJetPars_.end() ; ++ iter) {
	int gid = *iter;
	for(int id = *iter + par_->GetNumberOfJetParametersPerBin(); 
	    id < par_->GetNumberOfJetParameters() ; 
	    id += par_->GetNumberOfJetParametersPerBin()) {
	  temp_derivative1[gid] += temp_derivative1[id];
	  temp_derivative2[gid] += temp_derivative2[id];
	  temp_derivative1[id] = 0;
	  temp_derivative2[id] = 0;
	}
      }
      //zero derivative of fixed pars
      for( std::vector<int>::const_iterator iter = fixedJetPars_.begin();
	   iter != fixedJetPars_.end() ; ++ iter) {
	temp_derivative1[*iter] = 0;
	temp_derivative2[*iter] = 0;
      }
      for( std::vector<int>::const_iterator iter = fixedGlobalJetPars_.begin();
	   iter != fixedGlobalJetPars_.end() ; ++ iter) {
	temp_derivative1[*iter] = 0;
	temp_derivative2[*iter] = 0;
      }
      //fast derivative calculation:
      for( int param = 0 ; param < npar ; ++param ) {
	aux[param]      = temp_derivative1[param]/(2.0*derivStep_);
	aux[param+npar] = temp_derivative2[param]/(derivStep_*derivStep_);
	assert(aux[param] == aux[param]);
	assert(aux[param+npar] == aux[param+npar]);
      }
      //print derivatives:
      if(printParNDeriv_) {
	std::cout << std::setw(5) << "\npar";
	std::cout << std::setw(15) << "p";
	std::cout << std::setw(15) << "dp/dx";
	std::cout << std::setw(15) << "d^2p/dx^2\n";
	for( int param = 0 ; param < npar ; ++param ) {
	  std::cout << std::setw(5) << param;
	  std::cout << std::setw(15) << par_->GetPars()[param];
	  std::cout << std::setw(15) << aux[param];
	  std::cout << std::setw(15) << aux[param+npar] << std::endl;
	}
      }
      lvmfun_(par_->GetPars(),fsum,iret,aux);
      //par_->SetParameters(aux + par_index); 
      lvmprt_(2,aux,2); //print out
    } while (iret<0); 

    lvmprt_(2,aux,2); //print out
    for (int ithreads=0; ithreads<nThreads_; ++ithreads){
      t[ithreads]->ClearData();
    }  
    int par_index = 1;
    par_index = lvmind_(par_index);
    par_->SetParameters(aux + par_index);
  }
  //Copy Parameter errors from aux array to the TParameter::e array
  error_index=2;
  error_index = lvmind_(error_index);
  par_->SetErrors(aux+error_index);
  for( std::vector<int>::const_iterator iter = globalJetPars_.begin();
       iter != globalJetPars_.end() ; ++ iter) {
    double val =  par_->GetPars()[*iter];
    double err = par_->GetErrors()[*iter];
    for(int id = *iter + par_->GetNumberOfJetParametersPerBin(); 
	id < par_->GetNumberOfJetParameters() ; 
	id += par_->GetNumberOfJetParametersPerBin()) {
      par_->GetPars()[id] = val;
      par_->GetErrors()[id] = err;
    }
  }
  par_->SetFitChi2(fsum);
  
  for (int ithreads=0; ithreads<nThreads_; ++ithreads){
    delete t[ithreads];
  }
  delete [] aux;  
  delete [] temp_derivative1;
  delete [] temp_derivative2;
}



//--------------------------------------------------------------------------------------------
void TCaliber::done()
{
  ConfigFile config( configFile_.c_str() );

  // write calibration to cfi output file if ending is cfi
  bool cfi=false;
  bool txt=false;
  bool tex=false;
  std::string fileName(getOutputFile());
  if( fileName.find(".cfi")!=std::string::npos ){
    if( fileName.substr(fileName.find(".cfi")).compare(".cfi")==0 ){
      par_->writeCalibrationCfi( fileName.c_str() );
      cfi=true; // file has a real .cfi ending
    }
  }
  // write calibration to cfi output file if ending is txt
  if( fileName.find(".txt")!=std::string::npos ){
    if( fileName.substr(fileName.find(".txt")).compare(".txt")==0 ){
      par_->writeCalibrationTxt( fileName.c_str() );
      txt=true; // file has a real .txt ending
    }
  }
  // write calibration to cfi output file if ending is txt
  if( fileName.find(".tex")!=std::string::npos ){
    if( fileName.substr(fileName.find(".tex")).compare(".tex")==0 ){
      par_->writeCalibrationTex( fileName.c_str(), config );
      tex=true; // file has a real .txt ending
    }
  }

  // write calibration to cfi & txt output file if w/o ending
  if( !cfi && !txt && !tex ){
    par_->writeCalibrationCfi( (fileName+".cfi").c_str() );
    par_->writeCalibrationTxt( (fileName+".txt").c_str() );
    par_->writeCalibrationTex( (fileName+".tex").c_str(), config );
  }


  // Make control plots
  if( config.read<bool>("create plots",0) ) {
    int mode = config.read<int>("Mode",0);
    if( mode == 0 ) {  // Control plots for calibration
      TControlPlots * plots = new TControlPlots(configFile_,&data_,par_);
      plots->makePlots();
      delete plots;
    } else if( mode == 1 ) {  // Control plots for jetsmearing
      ControlPlotsJetSmearing * plotsjs = new ControlPlotsJetSmearing(configFile_,&data_,par_);
      plotsjs->plotResponse();
      //plotsjs->plotMeanResponseAndResolution();
      plotsjs->plotDijets();
      delete plotsjs;
    }
  }
  
  // Clean-up
  cout << endl << "Cleaning up... " << flush;
  for(DataIter i = data_.begin() ; i != data_.end() ; ++i) {
    delete *i;
  }
  data_.clear();
  cout << "Done" << endl;
}



//--------------------------------------------------------------------------------------------
void TCaliber::init()
{
  ConfigFile config(configFile_.c_str() );

  par_ = TParameters::CreateParameters(configFile_);

  //initialize temp arrays for fast derivative calculation
  TAbstractData::total_n_pars     = par_->GetNumberOfParameters();
  //--------------------------------------------------------------------------
  //read config file
  fitMethod_ = config.read<int>("Fit method",1);
  nThreads_ = config.read<int>("Number of Threads",1);
  
  // Residual scaling
  const char* resScheme = ( config.read<string>("Residual Scaling Scheme","221").c_str() );
  while(  *resScheme != 0  )
    {
      int scheme = static_cast<int>(*resScheme - '0');
      if(  scheme < 0  ||  scheme > 3  )
	{
	  cerr << "ERROR: " << scheme << " is not a valid scheme for resdiual scaling! Using default scheme 221." << endl << endl;
	  residualScalingScheme_.clear();
	  residualScalingScheme_.push_back(2);
	  residualScalingScheme_.push_back(2);
	  residualScalingScheme_.push_back(1);
	  break;
	}

      residualScalingScheme_.push_back( static_cast<int>(*resScheme - '0') );
      resScheme++;
    }
  outlierChi2Cut_        = config.read<double>("Outlier Cut on Chi2",100.0);

  //BFGS fit parameters
  derivStep_ = config.read<double>("BFGS derivative step",1e-03);
  mvec_       = config.read<int>("BFGS mvec",6);
  nIter_      = config.read<int>("BFGS niter",100);
  eps_        = config.read<double>("BFGS eps",1e-02);
  wlf1_       = config.read<double>("BFGS 1st wolfe parameter",1e-04);
  wlf2_       = config.read<double>("BFGS 2nd wolfe parameter",0.9);
  printParNDeriv_ = config.read<bool>("BFGS print derivatives",false);
  //global parameters ?
  globalJetPars_ = bag_of<int>(config.read<string>("global jet parameters","")); 

  //fixed jet parameters
  std::vector<int> fixJetPars = bag_of<int>(config.read<string>("fixed jet parameters",""));
  if(fixJetPars.size() % 3 == 0) {
    // Fix the specified parameters
    for(unsigned int i = 0 ; i < fixJetPars.size() ; i += 3) {
      int etaid = fixJetPars[i];
      int phiid = fixJetPars[i+1];
      int parid = fixJetPars[i+2];
      if(parid >= par_->GetNumberOfJetParametersPerBin()) continue;
      int jetbin = par_->GetJetBin(par_->GetJetEtaBin(etaid),par_->GetJetPhiBin(phiid));
      if(jetbin < 0) {
	std::cerr<<"WARNING: fixed jet parameter bin index = " << jetbin << endl; 
	exit(-2);  
      }
      //std::cout << "jetbin:" << jetbin << '\n';
      fixedJetPars_.push_back(jetbin * par_->GetNumberOfJetParametersPerBin() + par_->GetNumberOfTowerParameters() + parid);
    }
  }
  else if( fixJetPars.size() % 2 == 0 ) {
    // Fix all parameters in the specified jet bins
    for(unsigned int i = 0 ; i < fixJetPars.size() ; i += 2) {
      int etaid = fixJetPars[i];
      int phiid = fixJetPars[i+1];
      int jetbin = par_->GetJetBin(par_->GetJetEtaBin(etaid),par_->GetJetPhiBin(phiid));
      if(jetbin < 0) {
	std::cerr<<"WARNING: fixed jet parameter bin index = " << jetbin << endl; 
	exit(-2);  
      }

      for(int parIdx = 0; parIdx < par_->GetNumberOfJetParametersPerBin(); parIdx++) {
	fixedJetPars_.push_back( par_->GetNumberOfTowerParameters() + 
				 jetbin * par_->GetNumberOfJetParametersPerBin() +
				 parIdx );
      }
    }
  } else {
    cerr << "ERROR: Syntax error for fixed jet parameters. Syntax is:\n";
    cerr << "       'fixed jet parameter = { <eta_id> <phi_id> <par_id> }' or\n"; 
    cerr << "       'fixed jet parameter = { <eta_id> <phi_id> }'\n"; 
  }

  //fixed global parameters
  std::vector<int> fixGlobalJetPars = bag_of<int>(config.read<string>("fixed global jet parameters",""));
  for(size_t globalJetBin = 0; globalJetBin < fixGlobalJetPars.size(); globalJetBin++) {
      if( globalJetBin < 0 ) {
	std::cerr << "ERROR: fixed global jet parameter bin index = " << globalJetBin << std::endl; 
	exit(-2);  
      } else if( static_cast<int>(globalJetBin) > par_->GetNumberOfGlobalJetParameters() ) {
	std::cerr << "ERROR: fixed global jet parameter bin index = " << globalJetBin;
	std::cerr << " which is larger than the max number ";
	std::cerr << par_->GetNumberOfGlobalJetParameters() << " of global parameters." << std::endl;
	exit(-2);  
      } else {
	fixedGlobalJetPars_.push_back( par_->GetNumberOfTowerParameters() +
				       par_->GetNumberOfJetParameters()   +
				       par_->GetNumberOfTrackParameters() +
				       globalJetBin );
      }
  }

  outputFile_ = config.read<string>( "Output file", "calibration_k.cfi" );

  //fill data vector
  PhotonJetReader pjr(configFile_,par_);
  nGammajetEvents_ = pjr.readEvents(data_);
  
  DiJetReader djr(configFile_,par_);
  nDijetEvents_ = djr.readEvents(data_);

  TriJetReader tjr(configFile_,par_);
  nTrijetEvents_ = tjr.readEvents(data_);

  ZJetReader zjr(configFile_,par_);
  nZjetEvents_ = zjr.readEvents(data_);

  TopReader tr(configFile_,par_);
  nTopEvents_ = tr.readEvents(data_);
  
  TrackClusterReader tcr(configFile_,par_);
  nTrackClusterEvents_ = tcr.readEvents(data_);

  ParameterLimitsReader plr(configFile_,par_);
  plr.readEvents(data_);

  TowerConstraintsReader cr(configFile_,par_);
  cr.readEvents(data_);  

  JetConstraintsReader jcr(configFile_,par_);
  jcr.readEvents(data_);
}
//--^-TCaliber class-^------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
int caliber(int argc, char *argv[])
{
  std::cout << "The University Hamburg Calorimeter Calibration Tool, 2007/08/15." << std::endl;
  
  TCaliber * Calibration;
  if (argc>1)
    Calibration = new TCaliber( argv[1] );
  else  
    Calibration = new TCaliber("config/calibration.cfg"); //Read input defined in config file
  
  Calibration->init();
  Calibration->run();  //Run Fit
  Calibration->done(); //Do Plots & Write Calibration to file
  JetTruthEvent::printStats();
  Jet::printInversionStats();
  delete Calibration;    

  return 0;
}



//--------------------------------------------------------------------------------------------
void printUsage()
{
  std::cerr << "ERROR: You did something wrong! Better fix it." << std::endl;
}



//--------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  if (argc>2) {
    printUsage();
    exit(EXIT_FAILURE);
  }
  return caliber(argc, argv);
}

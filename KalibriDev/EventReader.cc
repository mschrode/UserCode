//
//    Base class for all event readers
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: EventReader.cc,v 1.2 2009/10/26 20:56:29 mschrode Exp $
//   
#include "EventReader.h"

#include "ConfigFile.h"
#include "Parameters.h" 


unsigned int EventReader::numberOfEventReaders_ = 0;

EventReader::EventReader(const std::string& configfile, TParameters* param) 
  : config_(0),par_(param)
{
  numberOfEventReaders_++;

  config_ = new ConfigFile(configfile.c_str());
  
  useTracks_ = config_->read<bool>("use Tracks",true);
  if(par_->GetNumberOfTrackParameters() < 1) useTracks_ = false;

  // Print info on track usage only once for all readers
  if( numberOfEventReaders_ == 1 ) {
    if(useTracks_) cout<<"Tracks are used to calibrate jets"<<endl;
    else cout<<"Only Calorimeter information is used"<<endl;
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
  string je = config_->read<string>("jet error parametrization","standard");
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
}

EventReader::~EventReader()
{
  delete config_;
}

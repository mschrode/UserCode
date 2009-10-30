//!
//!  \brief Reader for Parameter Limits
//!
//!  This class add user defined parameter limits
//!
//!  \author Hartmut Stadie
//!  \date  2008/12/12
//!  $Id: ParameterLimitsReader.cc,v 1.8 2009/10/26 20:56:29 mschrode Exp $
//!   
#include "ParameterLimitsReader.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"

#include <iostream>

ParameterLimitsReader::ParameterLimitsReader(const std::string& configfile, TParameters* p) :
  EventReader(configfile,p)
{
  vector<double> limits = bag_of<double>(config_->read<string>( "Jet Parameter Limits",""));
  
  // In case limits are explicitly set via config file
  if( limits.size() > 0 && limits.size() % 4 == 0 ) {
    std::cout << "Using user defined parameter limits:" << std::endl;
    for(unsigned int i = 0 ; i < limits.size() ; i += 4) {
      int index = (int)limits[i];
      for(int j = par_->GetNumberOfTowerParameters() + index; 
	  j <  par_->GetNumberOfParameters() ; 
	  j += par_->GetNumberOfJetParametersPerBin()) {
	par_limits.push_back(ParameterLimit(j,limits[i+1],limits[i+2],
					    limits[i+3]));
      }
    }
  }
  // In case default limits are to be used
  else if( limits.size() == 1 ) {
    string parclass = config_->read<string>("Parametrization Class","");
    std::cout << "Using default parameter limits for '" << parclass << "':" << std::endl;

    // For SmearHistGauss
    if( parclass == "SmearParametrizationStepGauss" ) {
      // Loop over parameters in one bin
      for(int i = 0; i < par_->GetNumberOfJetParametersPerBin(); i++) {
	double min = 0.;
	double max = 10000.;   // Sigma and histogrammed parameters have no upper limit
	if( i == 0 ) max = 1.; // Normalization constant between 0 and 1

	// Loop over eta and phi bins
	for(int j = par_->GetNumberOfTowerParameters() + i; 
	    j <  par_->GetNumberOfParameters(); 
	    j += par_->GetNumberOfJetParametersPerBin()) {
	  if( j < par_->GetNumberOfParameters()-1 )
	    par_limits.push_back(ParameterLimit(j,min,max,limits.at(0)));
	} // End of loop over eta and phi bins
      } // End of loop over parameters in one bin
    }
    // For SmearHistGaussInter
    else if( parclass == "SmearParametrizationStepGaussInter" ) {
      // Loop over parameters in one bin
      for(int i = 0; i < par_->GetNumberOfJetParametersPerBin(); i++) {
	double min = 0.;
	double max = 10000.;   // Sigma and histogrammed parameters have no upper limit
	if( i == 0 ) max = 1.; // Normalization constant between 0 and 1

	// Loop over eta and phi bins
	for(int j = par_->GetNumberOfTowerParameters() + i; 
	    j <  par_->GetNumberOfParameters(); 
	    j += par_->GetNumberOfJetParametersPerBin()) {
	  if( j < par_->GetNumberOfParameters()-1 )
	    par_limits.push_back(ParameterLimit(j,min,max,limits.at(0)));
	} // End of loop over eta and phi bins
      } // End of loop over parameters in one bin
    }
    // For SmearTwoGauss
    else if( parclass == "SmearParametrizationTwoGauss" ) {
      // Loop over parameters in one bin
      for(int i = 0; i < par_->GetNumberOfJetParametersPerBin(); i++) {
	double min = 0.;
	double max = 10000.;   // Sigma and histogrammed parameters have no upper limit
	if( i > 2 && (i-3)%3 == 0 ) max = 1.; // Normalization constant between 0 and 1

	// Loop over eta and phi bins
	for(int j = par_->GetNumberOfTowerParameters() + i; 
	    j <  par_->GetNumberOfParameters(); 
	    j += par_->GetNumberOfJetParametersPerBin()) {
	  if( j < par_->GetNumberOfParameters()-1 )
	    par_limits.push_back(ParameterLimit(j,min,max,limits.at(0)));
	} // End of loop over eta and phi bins
      } // End of loop over parameters in one bin
    }
    // For SmearStepGaussInterPtBinned
    else if( parclass == "SmearParametrizationStepGaussInterPtBinned" ) {
      int rNBins = config_->read<int>("Response pdf nsteps",10);
      // Loop over parameters in one bin
      for(int i = 0; i < par_->GetNumberOfJetParametersPerBin(); i++) {
	double min = 0.;
	double max = 10000.;   // Step function parameters have no upper limit
	if( i > 1 && (i-2)%(rNBins+1) == 0 ) max = 1.; // Normalization constant between 0 and 1

	// Loop over eta and phi bins
	for(int j = par_->GetNumberOfTowerParameters() + i; 
	    j <  par_->GetNumberOfParameters(); 
	    j += par_->GetNumberOfJetParametersPerBin()) {
	  if( j < par_->GetNumberOfParameters()-1 )
	    par_limits.push_back(ParameterLimit(j,min,max,limits.at(0)));
	} // End of loop over eta and phi bins
      } // End of loop over parameters in one bin
    }

  }
  // Catch wrong syntax in config file
  else if(limits.size() > 1) {
    std::cout << "wrong number of arguments for Jet Parameter Limits:" 
	      << limits.size() << '\n';
  }
  limits.clear();
  limits = bag_of<double>(config_->read<string>( "Tower Parameter Limits",""));
  
  if(limits.size() % 4 == 0) {
    for(unsigned int i = 0 ; i < limits.size() ; i += 4) {
      int index = (int)limits[i];
      for(int j = index; 
	  j <  par_->GetNumberOfTowerParameters() ; 
	  j += par_->GetNumberOfTowerParametersPerBin()) {
	par_limits.push_back(ParameterLimit(j,limits[i+1],limits[i+2],
					    limits[i+3]));
      }
    }
  } else if(limits.size() > 1) {
    std::cout << "wrong number of arguments for Tower Parameter Limits:" 
	      << limits.size() << '\n';
  }
}

ParameterLimitsReader::~ParameterLimitsReader()
{
}

int ParameterLimitsReader::readEvents(std::vector<TData*>& data)
{
  for(std::vector<ParameterLimit>::const_iterator pl = par_limits.begin() ;
      pl != par_limits.end() ; ++pl) {
    std::cout << " Adding limit for parameter " << std::flush;
    if     ( pl->index+1 < 10   ) std::cout << pl->index+1 << ":  " << std::flush;
    else if( pl->index+1 < 100  ) std::cout << pl->index+1 << ": "  << std::flush;
    else if( pl->index+1 < 1000 ) std::cout << pl->index+1 << ":"   << std::flush;
    std::cout << "  min: " << pl->min << "\t" << std::flush;
    std::cout <<  " max: " << pl->max << "\t" << std::flush;
    std::cout <<  " k: "   << pl->k << std::endl;

    TMeasurement* limitp  = new TMeasurement;
    limitp->pt  = pl->min;
    limitp->EMF = pl->max;
    TData_ParLimit * parlim = new TData_ParLimit(pl->index,limitp,pl->k,par_->GetPars()+pl->index,par_->parameter_limit);
    data.push_back(parlim);
  }
  return par_limits.size();
}
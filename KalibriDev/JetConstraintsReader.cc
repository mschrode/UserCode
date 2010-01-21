//
//    Reader for Jet Constraints
//
//    This class add user defined jet constraints
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: JetConstraintsReader.cc,v 1.8 2009/11/26 18:24:41 stadie Exp $
//   

#include "JetConstraintsReader.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"
#include "JetConstraintEvent.h"
#include "Parametrization.h"
#include "CorFactors.h"

#include <iostream>

JetConstraintsReader::JetConstraintsReader(const std::string& configfile, TParameters* p) :
  EventReader(configfile,p),cp(new ConstParametrization())
{
  //specify constraints
  std::vector<double> jet_constraint = bag_of<double>(config_->read<std::string>( "jet constraint",""));
  if(jet_constraint.size() % 4 == 0) {
    for(unsigned int i = 0 ; i < jet_constraint.size() ; i += 4) {
      jet_constraints.push_back(JetConstraint((int)jet_constraint[i],(int)jet_constraint[i+1],
					      jet_constraint[i+2],jet_constraint[i+3]));
    } 
  } else if(jet_constraint.size() > 1) {
    std::cout << "wrong number of arguments for jet constraint:" << jet_constraint.size() << '\n';
  }
}

JetConstraintsReader::~JetConstraintsReader()
{
}

int JetConstraintsReader::readEvents(std::vector<Event*>& data)
{
  for(std::vector<JetConstraint>::const_iterator ic = jet_constraints.begin() ;
      ic != jet_constraints.end() ; ++ic) {
    std::cout << "adding constraint for jets " << ic->mineta << " to " << ic->maxeta
	      << " for Et=" << ic->Et << " with weight w=" << ic->weight << "\n";
    //constrain average jet response
    int njets= (ic->maxeta - ic->mineta + 1) * 72;
    if((ic->maxeta  > 0) && (ic->mineta < 0)) njets -= 72;
    JetConstraintEvent* jce = new JetConstraintEvent(ic->Et,ic->weight);
    for(int ideta = ic->mineta ; ideta <= ic->maxeta  ; ++ideta) {
      if(ideta == 0) ideta = 1;
      jce->addJet(new  Jet(ic->Et,0,ic->Et,0,ic->Et,0,0,0.13,Jet::uds,ic->Et,0,
			   new CorFactors(1,1,1,1,1,1,1),par_->jet_function(ideta,1),
			   jet_error_param,Function(&Parametrization::correctedJetEt,0,0,0,0,cp),0));
    }
    data.push_back(jce);
  }
  return jet_constraints.size();
}

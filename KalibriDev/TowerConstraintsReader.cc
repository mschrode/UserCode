//
//    Reader for Tower Constraints
//
//    This class add user defined tower constraints
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: TowerConstraintsReader.cc,v 1.2 2009/10/26 20:56:30 mschrode Exp $
//   

#include "TowerConstraintsReader.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"

#include <iostream>

TowerConstraintsReader::TowerConstraintsReader(const std::string& configfile, TParameters* p) :
  EventReader(configfile,p)
{
  //specify constraints
  vector<double> tower_constraint = bag_of<double>(config_->read<string>( "Tower Constraint",""));
  if(tower_constraint.size() % 5 == 0) {
    for(unsigned int i = 0 ; i < tower_constraint.size() ; i += 5) {
      tower_constraints.push_back(TowerConstraint((int)tower_constraint[i],(int)tower_constraint[i+1],
						  tower_constraint[i+2],tower_constraint[i+3],
						  tower_constraint[i+4]));
    } 
  } else if(tower_constraint.size() > 1) {
    std::cout << "wrong number of arguments for tower constraint:" << tower_constraint.size() << '\n';
  }
}

TowerConstraintsReader::~TowerConstraintsReader()
{
}

int TowerConstraintsReader::readEvents(std::vector<TData*>& data)
{
  for(std::vector<TowerConstraint>::const_iterator ic = tower_constraints.begin() ;
      ic != tower_constraints.end() ; ++ic) {
    std::cout << "adding constraint for towers " << ic->mineta << " to " << ic->maxeta
	      << " for em Et=" << ic->emEt << " and had Et=" << ic->hadEt 
	      << " with weight w=" << ic->weight << "\n";
    //constrain average tower response
    int ntowers= (ic->maxeta - ic->mineta + 1) * 72;
    if((ic->maxeta  > 0) && (ic->mineta < 0)) ntowers -= 72;
    double etsum = ic->hadEt + ic->emEt;
    TMeasurement* constraintp  = new TMeasurement;
    constraintp->pt  = etsum;
    constraintp->eta = 0;
    constraintp->phi = 0;
    constraintp->E   = etsum;
   
    TData_TruthMultMess * tc = new TData_TruthMultMess(0,
						       etsum * ntowers, //truth
						       sqrt(pow(0.5,2)+ pow(0.1*etsum * ntowers,2)), //error
						       ic->weight, //weight
						       0, //params
						       0, //number of free jet param. p. bin
						       par_->dummy_parametrization, // function
						       par_->const_error<10000>, // function
						       constraintp);
    tc->SetType( typeTowerConstraint );
    //Add the towers to the event
    for(int ideta = ic->mineta ; ideta <= ic->maxeta  ; ++ideta) {
      if(ideta == 0) ideta = 1;
      for(int idphi = 1 ; idphi <= 72  ; ++idphi) {
	int index=par_->GetBin(par_->GetEtaBin(ideta),par_->GetPhiBin(idphi));
	if (index<0) {
	  cerr << "INDEX = "<< index << endl;
	  continue;
	}
	//create array with multidimensional measurement
	TMeasurement * mess = new TMeasurement;
	mess->pt = etsum;
	mess->EMF = ic->emEt;
	mess->HadF = ic->hadEt;
	mess->OutF = 0;
	mess->eta = 0;
	mess->phi = 0;
	mess->E = etsum;
	//mess[7] = double( cos( trackcluster.JetCalPhi-trackcluster.TowPhi[n] ) ); // Projection factor for summing tower Pt
	TData_TruthMess *tower = new TData_TruthMess(index,
						     mess, //mess
						     etsum, //"truth" for plotting only
						     sqrt(pow(0.5,2)+ pow(0.1*etsum,2)), //error
						     1.0, //weight ???
						     par_->GetTowerParRef(index), //parameter
						     par_->GetNumberOfTowerParametersPerBin(), //number of free cluster param. p. bin
						     par_->tower_parametrization, //function
						     par_->const_error<10> //error param.
						     );
	tc->AddMess(tower);
      } 
    } 
    data.push_back(tc);
  }
  return tower_constraints.size();
}

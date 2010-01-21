//
//!  \brief   Container class for jet correction factors
//
//    $Id: CorFactorsFactory.cc,v 1.2 2009/11/26 13:03:28 stadie Exp $
// 
#include "CorFactorsFactory.h"  
#include "CorFactors.h"

#include <iostream>

std::map<std::string,CorFactorsFactory*> CorFactorsFactory::map;


CorFactorsFactory::CorFactorsFactory(const std::string& name) : name_(name)
{
  static Cleaner cleaner;
  std::cout << "creating CorFactorsFactory: " << name_ << '\n';

  map[name] = this;
}

CorFactorsFactory::~CorFactorsFactory()
{
  std::cout << "deleting CorFactorsFactory: " << name_ << '\n';
  map[name_] = 0;
}

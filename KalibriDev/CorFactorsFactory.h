//!  \brief   Container class for jet correction factors
//
//    $Id: CorFactorsFactory.h,v 1.1 2009/11/26 10:27:48 stadie Exp $
//   
#ifndef CORFACTORSFACTORY_H
#define CORFACTORSFACTORY_H

#include <map>
#include <string>

class Jet;
class CorFactors;

class CorFactorsFactory
{
 public:
  CorFactorsFactory(const std::string& name);
  virtual ~CorFactorsFactory();
  virtual CorFactors* create(const Jet* j) = 0;

  static std::map<std::string,CorFactorsFactory*> map;

 private:
  CorFactorsFactory(const CorFactorsFactory&) {} 
  std::string name_;

  class Cleaner
  {
  public:
    Cleaner() {}
    ~Cleaner()
      {
	for( std::map<std::string,CorFactorsFactory*>::iterator i = map.begin();
	     i != map.end() ; ++i) delete i->second;
	map.clear();
      }
  };
  friend class Cleaner;
};

#endif

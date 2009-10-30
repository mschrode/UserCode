//!  \brief Run toy Monte Carlo
//!
//!  \author Hartmut Stadie
//!  \date   Mon Jun 30 11:00:00 CEST 2008
//!  $Id: toy.cc,v 1.9 2009/07/27 13:51:59 stadie Exp $
//!  \sa ToyMC
#include "ToyMC.h"

#include <iostream>

#include "ConfigFile.h"

int main(int argc, char* argv[]) {
  if( argc != 2 ) {
    std::cout << "Usage:\n    " << argv[0] << " <config.file>\n";
    return 1;
  }

  // Initialize ToyMC
  ToyMC * mc = new ToyMC();
  mc->init(argv[1]);
  mc->print();

  // Generate specified number of events
  ConfigFile config(argv[1]);
  std::string file = config.read<std::string>("ToyMC output file","input/toy.root");
  int nevents      = config.read<int>("ToyMC events",100);
  int type         = config.read<int>("ToyMC type",2);
  if(type == 1) 
    mc->makePhotonJet(file.c_str(),nevents);
  else if(type == 2) 
    mc->makeDiJet(file.c_str(),nevents);
  else if(type == 3)
    mc->makeTop(file.c_str(),nevents);    
  else {
    std::cerr << "unknown ToyMC type:" << type << '\n';
    delete mc;
    return 1;
  }
  delete mc;

  return 0;
}


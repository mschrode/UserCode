#include "Kalibri.h"
#include "JetTruthEvent.h"
#include "Jet.h"

#include <iostream>


int caliber(int argc, char *argv[])
{
  std::cout << "The University Hamburg Calorimeter Calibration Tool, 2007/08/15." << std::endl;
  
  Kalibri* Calibration;
  if (argc>1)
    Calibration = new Kalibri( argv[1] );
  else  
    Calibration = new Kalibri("config/calibration.cfg"); //Read input defined in config file
  
  Calibration->init();
  Calibration->run();  //Run Fit
  Calibration->done(); //Do Plots & Write Calibration to file
  JetTruthEvent::printStats();
  Jet::printInversionStats();
  delete Calibration;    

  return 0;
}


void printUsage()
{
  std::cerr << "ERROR: You did something wrong! Better fix it." << std::endl;
}



int main(int argc, char *argv[])
{
  if (argc>2) {
    printUsage();
    exit(EXIT_FAILURE);
  }
  return caliber(argc, argv);
}

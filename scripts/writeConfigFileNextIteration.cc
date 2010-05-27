#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"

bool goodLine(const std::string &line) {
  bool result = false;
  if( line.size() ) {
    if( line[0] != '#' ) {
      result = true;
    }
  }
  return result;
}
  

std::string id(const std::string &line, size_t &pos) {
  std::string result = "";
  pos = line.find("=");
  if( pos != std::string::npos ) {
    result = line.substr(0,pos);
    while( result[result.size()-1] == ' ' ) {
      result.erase(result.size()-1,result.size());
    }
  }
  return result;
}


void readFittedParameters(const std::string &fileName, std::vector<double> &pars) {
  TFile file(fileName.c_str(),"READ");
  assert( !file.IsZombie() );

  TH1 *h = 0;
  file.GetObject("hParameters",h);
  assert( h );
  pars.clear();
  for(int i = 0; i < h->GetNbinsX(); i++) {
    pars.push_back(h->GetBinContent(1+i));
  }
  file.Close();
}


void writeConfig(const std::string &iFileName, const std::string &oFilePath, const std::vector<double> &pars) {
  // Dummy sstream
  std::stringstream strstr;

  // Find out current iteration
  size_t pos = iFileName.rfind("_It");
  std::string strCurIt = iFileName.substr(pos+3,pos+4);
  strstr << strCurIt;
  int curIt = -1;
  strstr >> curIt;
  strstr.str("");

  // Open config of previous iteration
  std::ifstream iFile;
  iFile.open(iFileName.c_str());
  
  // Create config of next iteration
  curIt++;
  assert( curIt > 0 );
  strstr << curIt;
  std::string oFileName = iFileName;
  oFileName.replace(pos+3,1,strstr.str());
  strstr.str("");
  std::string path = oFilePath;
  oFileName = path.append("/"+oFileName);
  assert( oFileName != iFileName );
  std::ofstream oFile(oFileName.c_str());

  // Write new config
  if( iFile.is_open() && oFile.is_open() ) {
    std::string line = "";
    while( !iFile.eof() ) {
      getline(iFile,line);
      if( goodLine(line) ) {
	size_t pos = 0;
	std::string idt = id(line,pos);

	if( idt == "fixed jet parameters" ) {
	  oFile << idt << " = " << std::flush;
	  if( curIt % 2 == 0 ) {
	    oFile << "1 1 2   1 1 3   1 1 4   1 1 5" << std::endl;
	  } else {
	    oFile << "1 1 1   1 1 3   1 1 4   1 1 5" << std::endl;
	  }
	} else if( idt == "jet start values" ) {
	  oFile << idt << " = " << std::flush;
	  for(std::vector<double>::const_iterator it = pars.begin();
	      it != pars.end(); ++it) {
	    oFile << *it << "  " << std::flush;
	  }
	  oFile << std::endl;
	} else {
	  oFile << line << std::endl;
	}
      } else {
	oFile << line << std::endl;
      }
    }
  }  
  
}


int main(int argc, char *argv[]) {
  int result = 0;

  gSystem->Load("libHistPainter");

  if( argc > 3 ) {
    std::vector<double> pars;
    readFittedParameters(argv[2],pars);
    writeConfig(argv[1],argv[3],pars);
  } else {
    std::cerr << "ERROR: too few arguments.\n";
    result = -1;
  }

  return result;
}

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>


typedef std::map< std::string, std::vector<std::string> > ValueMap;
typedef std::map< std::string, std::vector<std::string> >::const_iterator ValueMapIt;

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


ValueMap readValues(const std::string &fileName) {

  ValueMap vMap;
  size_t nIterations = 0;

  std::ifstream file;
  file.open(fileName.c_str());
  if( file.is_open() ) {
    std::string line = "";
    while( !file.eof() ) {
      getline(file,line);
      if( goodLine(line) ) {
	size_t pos = 0;
	std::string idt = id(line,pos);
	pos++;
	std::vector<std::string> vals;
	for(; pos < line.size(); ++pos) {
	  std::string val = line.substr(pos,line.find("|",pos)-pos);
	  vals.push_back(val);
	  pos += val.size();
	}
	if( nIterations > 0 ) {
	  if( vals.size() != nIterations ) {
	    std::cerr << "ERROR: different number of iterations\n";
	    exit(-1);
	  }
	} else {
	  nIterations = vals.size();
	}
	vMap[idt] = vals;
      }
    }
  } else {
    std::cerr << "ERROR opening file '" << fileName << "'\n";
    exit(1);
  }
  file.close();


  return vMap;
}


void writeConfig(const std::string &iFileName, const ValueMap &vMap, int iteration) {
  std::ifstream iFile;
  iFile.open(iFileName.c_str());
  
  std::string oFileName = iFileName.substr(0,iFileName.find("."));
  oFileName.append("_");
  std::stringstream toString;
  toString << iteration;
  oFileName.append(toString.str());
  oFileName.append(".cfg");
  std::ofstream oFile(oFileName.c_str());

  if( iFile.is_open() && oFile.is_open() ) {
    std::string line = "";
    while( !iFile.eof() ) {
      getline(iFile,line);
      if( goodLine(line) ) {
	size_t pos = 0;
	std::string idt = id(line,pos);

	ValueMapIt vMapIt = vMap.find(idt);
	if( vMapIt != vMap.end() ) {
	  oFile << idt << " = " << vMapIt->second.at(iteration) << std::endl;
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

  if( argc > 2 ) {
    ValueMap vMap = readValues(argv[2]);
    if( vMap.size() ) {
      for(size_t i = 0; i < vMap.begin()->second.size(); ++i) {
	writeConfig(argv[1],vMap,i);
      }
    }
  } else {
    std::cerr << "ERROR: too few arguments.\n";
    result = -1;
  }

  return result;
}

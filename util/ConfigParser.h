// $Id: $

#ifndef CONFIG_PARSER_H
#define CONFIG_PARSER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace util {
  //! Parse text files for lines of the form
  //!   'tag <delimiter> values'
  // -------------------------------------------------------------------------------------
  class ConfigParser {
  public:
    ConfigParser(const std::string &name) {
      file_ = new std::ifstream(name.c_str());
    }

    ~ConfigParser() {
      file_->close();
      delete file_;
    }

    std::string readString(const std::string &tag, const std::string &delim = ":") const;
    std::vector<std::string> readStringVec(const std::string &tag, const std::string &delim = ":") const;
    std::vector<double> readDoubleVec(const std::string &tag, const std::string &delim = ":") const;


  private:
    std::ifstream* file_;

    double stringToDouble(const std::string &str) const;
    std::string getTag(const std::string &line, const std::string &delim, size_t &pos) const;
    bool noComment(const std::string &line) const;
    void trim(std::string &str) const;
    void trim(std::string &str, const std::string delims) const;
  };



  // -------------------------------------------------------------------------------------
  std::string ConfigParser::readString(const std::string &tag, const std::string &delim) const {
    std::string result;
    if( file_->is_open() ) {
      std::string line = "";
      while( !file_->eof() ) {
	std::getline(*file_,line);
	if( noComment(line) ) {
	  size_t pos = 0;
	  if( getTag(line,delim,pos) == tag) {
	    pos += delim.size();
	    result = line.substr(pos);
	    break;
	  }
	}
      }
      trim(result);
    } else {
      std::cerr << "ConfigParser: ERROR reading from file\n";
      exit(1);
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  std::vector<std::string> ConfigParser::readStringVec(const std::string &tag, const std::string &delim) const {
    std::vector<std::string> result;
    std::string line = readString(tag,delim);
    std::string delims = " ;|";
    while( line.size() ) {
      trim(line,delims);
      size_t end = line.find_first_of(delims);
      result.push_back(line.substr(0,end));
      line.erase(0,end);
    }    
    return result;
  }


  // -------------------------------------------------------------------------------------
  std::vector<double> ConfigParser::readDoubleVec(const std::string &tag, const std::string &delim) const {
    std::vector<double> result;
    std::vector<std::string> strVec = readStringVec(tag,delim);
    for(std::vector<std::string>::const_iterator it = strVec.begin();
	it != strVec.end(); ++it) {
      result.push_back(stringToDouble(*it));
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  bool ConfigParser::noComment(const std::string &line) const {
    bool result = false;
    if( line.size() ) {
      if( line[0] != '#' ) {
	result = true;
      }
    }
    return result;
  }
  

  // -------------------------------------------------------------------------------------
  std::string ConfigParser::getTag(const std::string &line, const std::string &delim, size_t &pos) const {
    std::string result = "";
    pos = line.find(delim);
    if( pos != std::string::npos ) {
      result = line.substr(0,pos);
      trim(result);
    }
    return result;
  }


  // -------------------------------------------------------------------------------------
  void ConfigParser::trim(std::string &str) const {
    while( str.size() && str[0] == ' ' ) str.erase(0,1);
    while( str.size() && str[str.size()-1] == ' ' ) str.erase(str.size()-1,str.size());    
  }


  // -------------------------------------------------------------------------------------
  void ConfigParser::trim(std::string &str, const std::string delims) const {
    while( str.size() && str.find_first_of(delims) == 0 ) str.erase(0,1);
    while( str.size() && str.find_last_of(delims) == str.size()-1 ) str.erase(str.size()-1,str.size());    
  }


  // -------------------------------------------------------------------------------------
  double ConfigParser::stringToDouble(const std::string &str) const {
    std::istringstream stm;
    stm.str(str);
    double d;
    stm >> d;
    return d;
  }
}
#endif

// $Id: ControlPlotsConfig.cc,v 1.1 2010/01/04 17:04:51 mschrode Exp $

#include "ControlPlotsConfig.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include "TFile.h"

#include "ConfigFile.h"


//! \brief Constructor
//!
//! Reads the parameters for a profile control plot
//! of name \p name from the configuration file
//! \p configFile.
// --------------------------------------------------
ControlPlotsConfig::ControlPlotsConfig(const ConfigFile *configFile, const std::string &name)
  : config_(configFile), name_(name) {
  init();
}



//! This is the name used e.g. for histogram names and
//! consists of "binVariable()+binIdx".
// --------------------------------------------------
std::string ControlPlotsConfig::binName(int binIdx) const { 
  std::string name = binVariable();
  name += toString(binIdx);

  return name;
}



//! This is used for the profile histogram title
//! and consists of "min < varTitle(binning) < max (unit)".
// --------------------------------------------------
std::string ControlPlotsConfig::binTitle(double min, double max) const {
  std::string title = toString(min);
  title += " < ";
  title += varTitle(binVariable());
  title += " < ";
  title += toString(max);

  std::string unit = unitTitle(binVariable());
  if( unit != "" ) {
    title += " " + unit;
  }

  return title;
}



//! This is the name used e.g. for histogram names and
//! consists of "xVariable()+xBinIdx".
// --------------------------------------------------
std::string ControlPlotsConfig::xBinName(int xBinIdx) const { 
  std::string name = xVariable();
  name += toString(xBinIdx);

  return name;
}



//! This is a label specifying the y distribution
//! in x bin \p xBinIdx and with a value of the
//! binning variable between \p binMin and \p binMax.
//! It is drawn on the histogram and consists of
//! "min < varTitle() < max (unit), xBinMin < varTitle(x) < xBinMax (unit)"
// --------------------------------------------------
std::string ControlPlotsConfig::xBinTitle(int xBinIdx, double binMin, double binMax) const {
  std::string title = binTitle(binMin,binMax);
  if( xBinIdx >= 0 && xBinIdx < nXBins() ) {
    title += ",  ";
    title += toString(round(xBinEdges_.at(xBinIdx)));
    title += " < " + varTitle(xVariable()) + " < ";
    title += toString(round(xBinEdges_.at(xBinIdx+1)));
    std::string unit = unitTitle(xVariable());
    if( unit != "" ) {
      title += " " + unit;
    }
  }

  return title;
}



// --------------------------------------------------
std::string ControlPlotsConfig::yProfileTitle(ProfileType type) const {
  std::string title = "";

  if( type == Mean )
    title = "< " + yTitle() + ">";
  else if( type == StandardDeviation )
    title = "#sigma( " + yTitle() + ") / <" + yTitle() + ">";
  else if( type == GaussFitMean )
    title = "GaussFit < " + yTitle() + ">";
  else if( type == GaussFitWidth )
    title = "GaussFit #sigma( " + yTitle() + ") / <" + yTitle() + ">";
  else if( type == Median )
    title = "Median " + yTitle();
  else if( type == Chi2 )
    title = "#chi^{2} / ndof " + yTitle();
  else if( type == Probability )
    title = "Probability " + yTitle();
  else if( type == Quantiles )
    title = "Quantiles " + yTitle();
  else
    std::cerr << "WARNING: Undefined ProfileType '" << type << "'\n";    

  return title;
}



// --------------------------------------------------
int ControlPlotsConfig::color(CorrectionType type) const {
  int color = 1;
  std::map<CorrectionType,int>::const_iterator it = colors_.find(type);
  if( it != colors_.end() ) color = it->second;

  return color;
}



// --------------------------------------------------
int ControlPlotsConfig::markerStyle(CorrectionType type) const {
  int markerStyle = 7;
  std::map<CorrectionType,int>::const_iterator it = markerStyles_.find(type);
  if( it != markerStyles_.end() ) markerStyle = it->second;

  return markerStyle; 
}



// --------------------------------------------------
std::string ControlPlotsConfig::legendLabel(CorrectionType type) const {
  std::string label = "DEFAULT";
  std::map<CorrectionType,std::string>::const_iterator it = legendLabels_.find(type);
  if( it != legendLabels_.end() ) label = it->second;

  return label;
}



//! Possible names \p typeName are
//! - "Uncorrected": \p CorrectionType::Uncorrected
//! - "Kalibri": \p CorrectionType::Kalibri
//! - "L2L3": \p CorrectionType::L2L3
// --------------------------------------------------
ControlPlotsConfig::CorrectionType ControlPlotsConfig::correctionType(const std::string &typeName) const {
  CorrectionType type = Uncorrected;

  if( typeName == "Uncorrected" )
    type = Uncorrected;
  else if( typeName == "Kalibri" )
    type = Kalibri;
  else if( typeName == "L2L3" )
    type = L2L3;
  else
    std::cerr << "WARNING: Undefined CorrectionType '" << typeName << "'\n";

  return type;
}



//! Possible correction types \p corrType are
//! - \p CorrectionType::Uncorrected: "Uncorrected"
//! - \p CorrectionType::Kalibri: "Kalibri"
//! - \p CorrectionType::L2L3: "L2L3"
// --------------------------------------------------
std::string ControlPlotsConfig::correctionTypeName(CorrectionType corrType) const {
  std::string name = "corrTypeName";

  if( corrType == Uncorrected )
    name = "Uncorrected";
  else if( corrType == Kalibri )
    name = "Kalibri";
  else if( corrType == L2L3 )
    name = "L2L3";
  else
    std::cerr << "WARNING: Undefined CorrectionType '" << corrType << "'\n";

  return name;
}



//! Possible names \p typeName are
//! - "Mean": \p ProfileType::Mean
//! - "StandardDeviation": \p ProfileType::StandardDeviation
//! - "GaussFitMean": \p ProfileType::GaussFitMean
//! - "GaussFitWidth": \p ProfileType::GaussFitWidth
//! - "Median": \p ProfileType::Median
//! - "Chi2": \p ProfileType::Chi2
//! - "Probability": \p ProfileType::Probability
//! - "Quantiles": \p ProfileType::Quantiles
// --------------------------------------------------
ControlPlotsConfig::ProfileType ControlPlotsConfig::profileType(const std::string &typeName) const {
  ProfileType type = GaussFitMean;

  if( typeName == "Mean" )
    type = Mean;
  else if( typeName == "StandardDeviation" )
    type = StandardDeviation;
  else if( typeName == "GaussFitMean" )
    type = GaussFitMean;
  else if( typeName == "GaussFitWidth" )
    type = GaussFitWidth;
  else if( typeName == "Median" )
    type = Median;
  else if( typeName == "Chi2" )
    type = Chi2;
  else if( typeName == "Probability" )
    type = Probability;
  else if( typeName == "Quantiles" )
    type = Quantiles;
  else
    std::cerr << "WARNING: Undefined ProfileType '" << typeName << "'\n";

  return type;
}



//! Possible profile types \p profType are
//! - \p ProfileType::Mean: "Mean" 
//! - \p ProfileType::StandardDeviation: "StandardDeviation" 
//! - \p ProfileType::GaussFitMean: "GaussFitMean" 
//! - \p ProfileType::GaussFitWidth: "GaussFitWidth" 
//! - \p ProfileType::Median: "Median" 
//! - \p ProfileType::Chi2: "Chi2" 
//! - \p ProfileType::Probability: "Probability" 
//! - \p ProfileType::Quantiles: "Quantiles" 
// --------------------------------------------------
std::string ControlPlotsConfig::profileTypeName(ProfileType profType) const {
  std::string name = "ProfileTypeName";

  if( profType == Mean )
    name = "Mean";
  else if( profType == StandardDeviation )
    name = "StandardDeviation";
  else if( profType == GaussFitMean )
    name = "GaussFitMean";
  else if( profType == GaussFitWidth )
    name = "GaussFitWidth";
  else if( profType == Median )
    name = "Median";
  else if( profType == Chi2 )
    name = "Chi2";
  else if( profType == Probability )
    name = "Probability";
  else if( profType == Quantiles )
    name = "Quantiles";
  else
    std::cerr << "WARNING: Undefined ProfileType '" << profType << "'\n";    

  return name;
}



//!  Write \p obj to the file \p outFileName_ into the
//!  directory \p outDirName(). Inside the ROOT file,
//!  \p obj is written into the directory \p outFile_:/name().
//!  If it does not exist, it is created first.
//---------------------------------------------------------------
void ControlPlotsConfig::toRootFile(TObject *obj) const {
  // Create / open ROOT file for output
  TFile *outFile = new TFile((outDirName()+"/KalibriPlots.root").c_str()
			     ,"UPDATE","Kalibri control plots");

  std::string directory = outFile->GetName();
  directory += ":";
  gDirectory->cd(directory.c_str());
  directory += "/" + name();
  bool dirExists = gDirectory->GetDirectory(directory.c_str());
  if( !dirExists ) {
    gDirectory->mkdir(name().c_str());
  }
  gDirectory->cd(directory.c_str());
  if( !(gDirectory->WriteTObject(obj)) ) {
    std::cerr << "ERROR writing object '" << obj->GetName() << "' to ROOT file." << std::endl;
  }

  outFile->Close();
  delete outFile;
}



// --------------------------------------------------
void ControlPlotsConfig::init() {
  // Read binning
  std::vector<std::string> strVar = bag_of_string(config_->read<std::string>(name_+" bin variable","Eta"));
  binVar_ = strVar.at(0);
  binEdges_ = bag_of<double>(config_->read<std::string>(name_+" bin edges","0. 1."));
  assert( binEdges_.size() > 1 );
  for(size_t i = 1; i < binEdges_.size(); i++) {
    assert( binEdges_.at(i) > binEdges_.at(i-1) );
  }
  nBins_ = static_cast<int>(binEdges_.size())-1;

  // Read x axis
  strVar = bag_of_string(config_->read<std::string>(name_+" x variable","GenJetPt"));
  xVar_ = strVar.at(0);
  logX_ = false;
  if( strVar.size() == 2 ) {
    if( strVar.at(1) == "log" ) {
      logX_ = true;
    }
  }
  double min = 20.;
  double max = 100.;
  std::vector<double> var = bag_of<double>(config_->read<std::string>(name_+" x edges","15 10 1000"));
  if( var.size() == 3 ) {
    nXBins_ = static_cast<int>(var.at(0));
    min = var.at(1);
    max = var.at(2);
  } else {
    std::cerr << "WARNING: Wrong number of arguments in config line '" << name_ << " x edges'\n";
    nXBins_ = 5;
  }
  xBinEdges_ = std::vector<double>(nXBins_+1);
  if( logX_ ) {
    if( !equidistLogBins(xBinEdges_,nXBins_,min,max) )
      std::cerr << "ERROR creating equidistant logarithmic binning.\n";
  } else {
    double width = (max - min) / nXBins_;
    for(int i = 0; i < nXBins_+1; i++) {
      xBinEdges_.at(i) = min + width*i;
    }
  }
  for(int i = 0; i < nXBins_; i++) {
    assert( xBinEdges_.at(i) < xBinEdges_.at(i+1) );
  }

  // Read y axis
  strVar = bag_of_string(config_->read<std::string>(name_+" y variable","GenJetResponse"));
  yVar_ = strVar.at(0);
  min = 0.;
  max = 2.;
  var = bag_of<double>(config_->read<std::string>(name_+" y edges","51 0 2 0.5 1.5"));
  if( var.size() == 5 ) {
    nYBins_ = static_cast<int>(var.at(0));
    min = var.at(1);
    max = var.at(2);
    yMinZoom_ = var.at(3);
    yMaxZoom_ = var.at(4);
  } else {
    std::cerr << "WARNING: Wrong number of arguments in config line '" << name_ << " y edges'\n";
    nYBins_ = 51;
    min = 0.;
    max = 2.;
    yMinZoom_ = 0.5;
    yMaxZoom_ = 1.5;
  }
  yBinEdges_ = std::vector<double>(nYBins_+1);
  double width = (max - min) / nYBins_;
  for(int i = 0; i < nYBins_+1; i++) {
    yBinEdges_.at(i) = min + width*i;
  }
  for(int i = 0; i < nYBins_; i++) {
    assert( yBinEdges_.at(i) < yBinEdges_.at(i+1) );
    assert( yMinZoom_ < yMaxZoom_ );
  }

  // Store which correction types are to be drawn
  // in the profile plots
  std::vector<std::string> corrTypesStr = bag_of_string(config_->read<std::string>(name_+" correction types","Uncorrected"));
  for(std::vector<std::string>::const_iterator corrTypesIt = corrTypesStr.begin();
      corrTypesIt != corrTypesStr.end(); corrTypesIt++) {
    corrTypes_.push_back(correctionType(*corrTypesIt));
  }

  // Store which correction types are to be drawn
  // in the distributions
  corrTypesStr = bag_of_string(config_->read<std::string>(name_+" distributions",";"));
  for(std::vector<std::string>::const_iterator corrTypesIt = corrTypesStr.begin();
      corrTypesIt != corrTypesStr.end(); corrTypesIt++) {
    corrTypesDistributions_.push_back(correctionType(*corrTypesIt));
  }

  // Store which profile types are to be drawn
  std::vector<std::string> profTypesStr = bag_of_string(config_->read<std::string>(name_+" profile types","Uncorrected"));
  for(std::vector<std::string>::const_iterator profTypesIt = profTypesStr.begin();
      profTypesIt != profTypesStr.end(); profTypesIt++) {
    profTypes_.push_back(profileType(*profTypesIt));
  }

  // Store directory name for output
  outDirName_ = config_->read<std::string>("plots output directory","controlPlots");

  // Define style for different correction types
  // This should become configurable via config file
  colors_[Uncorrected] = 1;
  colors_[Kalibri] = 2;
  colors_[L2L3] = 4;

  markerStyles_[Uncorrected] = 20;
  markerStyles_[Kalibri] = 21;
  markerStyles_[L2L3] = 23;

  // Define default legend labels for the different corrections
  legendLabels_[Uncorrected] = "Uncorrected";
  legendLabels_[Kalibri] = "Kalibri";
  legendLabels_[L2L3] = "L2L3";
  // Read optional legend labels
  std::vector<std::string> legLabelStr = bag_of_string(config_->read<std::string>(name_+" legend label",";"));
  for(std::vector<std::string>::const_iterator legLabelIt = legLabelStr.begin();
      legLabelIt != legLabelStr.end(); legLabelIt++) {
    size_t pos = legLabelIt->find(":");
    if( pos != std::string::npos ) {
      legendLabels_[correctionType(legLabelIt->substr(0,pos))] = legLabelIt->substr(pos+1);
    }
  }
}



//!  Filling \p bins with borders of \p nBins bins between \p first
//!  and \p last that are equidistant when viewed in log scale,
//!  so \p bins must have length \p nBins+1. If \p first, \p last
//!  or \p nBins are not positive, failure is reported.
// -------------------------------------------------------------
bool ControlPlotsConfig::equidistLogBins(std::vector<double>& bins, int nBins, double first, double last) const {
  if( nBins < 1 || first <= 0. || last <= 0. || first >= last ) return false;

  bins[0]     = first;
  bins[nBins] = last;
  const double firstLog = log10(bins[0]);
  const double lastLog  = log10(bins[nBins]);
  for (int i = 1; i < nBins; ++i) {
    bins[i] = pow(10., firstLog + i*(lastLog-firstLog)/(nBins));
  }

  return true;
}



//! The title consists of "varTitle(varName) (unit)"
// --------------------------------------------------
std::string ControlPlotsConfig::axisTitle(const std::string &varName) const {
  std::string title = varTitle(varName);
  std::string unit = unitTitle(varName);
  if( unit != "" ) {
    title += " (" + unit + ")";
  }

  return title;
}



// --------------------------------------------------
std::string ControlPlotsConfig::unitTitle(const std::string &varName) const {
  std::string title = "";

  if( varName == "GenJetPt" )
    title = "GeV";

  return title;
}



// --------------------------------------------------
std::string ControlPlotsConfig::varTitle(const std::string &varName) const {
  std::string title = "";
  
  if( varName == "Eta" )
    title = "#eta";
  else if( varName == "GenJetPt" )
    title = "p^{gen}_{T}";
  else if( varName == "GenJetResponse" )
    title = "p_{T} / p^{gen}_{T}";

  return title;
}


// --------------------------------------------------
template <class T> std::string ControlPlotsConfig::toString(const T& t) const {
  std::stringstream ss;
  ss << t;
  return ss.str();
}

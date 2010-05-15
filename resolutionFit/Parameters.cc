// $Id:  $

#include "Parameters.h"

#include <cassert>


namespace resolutionFit {
  Parameters::Parameters(double etaMin, double etaMax, const TString &fileBaseNameStdSel, const std::vector<double> &ptBinEdges, int startIdx, int endIdx, const TString &outNamePrefix, int verbosity)
    : etaMin_(etaMin), etaMax_(etaMax),
      startIdx_(startIdx), endIdx_(endIdx),
      outNamePrefix_(outNamePrefix), verbosity_(verbosity) {

    ptBinEdges_ = ptBinEdges;
    trueResPar_ = std::vector<double>(3,0.);

    writeFileNames(namesStdSel_,fileBaseNameStdSel);

    assert( nPtBins() > 0 );
    assert( (endIdx-startIdx) == nPtBins()-1 );
  }

  Parameters::~Parameters() {
    for(std::list<PtBinParameters*>::iterator it = listOfPtBinParameters_.begin();
	it != listOfPtBinParameters_.end(); ++it) {
      delete *it;
    }
  }


  const Parameters::PtBinParameters* Parameters::createPtBinParameters(int ptBinIdx) const {
    PtBinParameters *par = new PtBinParameters(ptBinIdx,this);
    listOfPtBinParameters_.push_back(par);
    return par;
  }


  void Parameters::addPt3Cut(double pt3RelCutValue, const TString &fileBaseName) {
    pt3RelCutValues_.push_back(pt3RelCutValue);
    std::vector<TString> names;
    writeFileNames(names,fileBaseName);
    namesCutVars_.push_back(names);
  }


  void Parameters::addSystUncert(const TString &label, const TString &fileBaseNameDown, const TString &fileBaseNameUp) {
    labelSyst_.push_back(label);
    std::vector<TString> names;
    writeFileNames(names,fileBaseNameDown);
    namesSystDown_.push_back(names);
    writeFileNames(names,fileBaseNameUp);
    namesSystUp_.push_back(names);
  }


  void Parameters::writeFileNames(std::vector<TString> &names, const TString &baseName) const {
    names.clear();
    for(int bin = startIdx_; bin <= endIdx_; ++bin) {
      names.push_back(baseName);
      names.back() += bin;
      names.back() += "/jsResponse.root";
    }
  }
}

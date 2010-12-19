// $Id: FileOps.h,v 1.5 2010/12/04 16:18:05 mschrode Exp $

#ifndef FileOps_h
#define FileOps_h

#include <iostream>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TObject.h"
#include "TString.h"

#include "utils.h"
#include "HistOps.h"

namespace util
{
  class FileOps
  {
  public:
    static TF1* readTF1(const TString &fileName, const TString &fName, const TString &newFName = "");
    static TH2* readTH2(const TString &fileName, const TString &hName, const TString &newHName = "");
    static TH1* readTH1(const TString &fileName, const TString &histName, const TString &newHistName = "");
    static util::HistVec readTH1(const std::vector<TString> &fileNames, const TString &histName, const TString &newHistName = "");
    static util::HistVec readHistVec(const TString &fileName, const TString &histName, const TString &newHistName = "");
  };


  // -------------------------------------------------------------------------------------
  TH2* FileOps::readTH2(const TString &fileName, const TString &hName, const TString &newHName) {
    TFile file(fileName,"READ");
    TH2* h = 0;
    file.GetObject(hName,h);
    if( h ) {
      h->SetDirectory(0);
      if( newHName.Length() ) h->SetName(newHName);
    } else {
      std::cerr << "ERROR in FileOps::readTH2: No TH2 with name '" << hName << "' in file '" << fileName << "'\n.";
    }
    file.Close();
    
    return h;
  }


  //! Read TH1 histogram from file
  // -------------------------------------------------------------------------------------
  TH1* FileOps::readTH1(const TString &fileName, const TString &histName, const TString &newHistName) {
    TFile file(fileName,"READ");
    TH1 *h = 0;
    file.GetObject(histName,h);
    if( h ) {
      h->SetDirectory(0);
      if( newHistName.Length() ) h->SetName(newHistName);
    } else {
      std::cerr << "ERROR in FileOps::readTH1: Histogram with name '" << histName << "' does not exist in file '" << fileName << "'\n.";
    }
    file.Close();
    
    return h;
  }
  
  
  //! Read TH1 histograms from different files
  // -------------------------------------------------------------------------------------
  util::HistVec FileOps::readTH1(const std::vector<TString> &fileNames, const TString &histName, const TString &newHistName) {
    util::HistVec v(fileNames.size());
    for(unsigned int i = 0; i < fileNames.size(); ++i) {
      v[i] = readTH1(fileNames[i],histName,newHistName+util::toTString(i));
    }
    std::cout << "Done\n";
    
    return v;
  }
  
  
  //! Read TH1 histograms from one file
  // -------------------------------------------------------------------------------------
  util::HistVec FileOps::readHistVec(const TString &fileName, const TString &histName, const TString &newHistName) {
    util::HistVec v;
    TFile file(fileName,"READ");
    bool binExists = true;
    unsigned int bin = 0;
    while( binExists ) {
      TH1 *h = 0;
      file.GetObject(histName+util::toTString(bin),h);
      if( h ) {
	h->SetDirectory(0);
	if( newHistName.Length() ) h->SetName(newHistName+util::toTString(bin));
	v.push_back(h);
	++bin;
      } else {
	binExists = false;
      }
    }
    file.Close();
    
    if( v.size() == 0 ) std::cerr << "WARNING in util::FileOps::readHistVec(): No histogram read!\n";
    
    return v;
  }


  // -------------------------------------------------------------------------------------
  TF1* FileOps::readTF1(const TString &fileName, const TString &fName, const TString &newFName) {
    TFile file(fileName,"READ");
    TF1 *f = 0;
    file.GetObject(fName,f);
    if( f ) {
      //f->SetDirectory(0);
      if( newFName.Length() ) f->SetName(newFName);
    } else {
      std::cerr << "ERROR in FileOps::readTF1: TF1 with name '" << fName << "' does not exist in file '" << fileName << "'\n.";
    }
    file.Close();
    
    return f;
  }
}
#endif

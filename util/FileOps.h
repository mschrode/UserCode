// $Id: FileOps.h,v 1.9 2011/05/17 16:36:43 mschrode Exp $

#ifndef FileOps_h
#define FileOps_h

#include <fstream>
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
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
    static TGraph* readTGraph(const TString &fileName, const TString &gName, const TString &newGName = "");
    static TH2* readTH2(const TString &fileName, const TString &hName, const TString &newHName = "");
    static TH1* readTH1(const TString &fileName, const TString &histName, const TString &newHistName = "");
    static util::HistVec readTH1(const std::vector<TString> &fileNames, const TString &histName, const TString &newHistName = "");
    static util::HistVec readHistVec(const TString &fileName, const TString &histName, const TString &newHistName = "");
    static std::vector<TF1*> readTF1Vec(const TString &fileName, const TString &fName, const TString &newFName = "");
    static TChain* createTChain(const TString &fileName, const TString &treeName = "", unsigned int verbosity = 1);
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


  //! Read TGraph from file
  // -------------------------------------------------------------------------------------
  TGraph* FileOps::readTGraph(const TString &fileName, const TString &gName, const TString &newGName) {
    TFile file(fileName,"READ");
    TGraph *g = 0;
    file.GetObject(gName,g);
    if( g ) {
      if( newGName.Length() ) g->SetName(newGName);
    } else {
      std::cerr << "ERROR in FileOps::readTGraph: TGraph with name '" << gName << "' does not exist in file '" << fileName << "'\n.";
    }
    file.Close();
    
    return g;
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


  // -------------------------------------------------------------------------------------
  std::vector<TF1*> FileOps::readTF1Vec(const TString &fileName, const TString &fName, const TString &newFName) {
    std::vector<TF1*> v;
    TFile file(fileName,"READ");
    bool binExists = true;
    unsigned int bin = 0;
    while( binExists ) {
      TF1* f = 0;
      file.GetObject(fName+util::toTString(bin),f);
      if( f ) {
	//f->SetDirectory(0);
	if( newFName.Length() ) f->SetName(newFName+util::toTString(bin));
	v.push_back(f);
	++bin;
      } else {
	binExists = false;
      }
    }
    file.Close();
    
    if( v.size() == 0 ) std::cerr << "WARNING in util::FileOps::readTF1Vec(): No TF1 read!\n";
    
    return v;
  }



  //! Create TChain from input root files. The root
  //! files are expected to contain a TTree "DiJetTree".
  //! There are two possible input options:
  //!
  //! 1) 'fileName' specifies a single root file; it ends
  //!    with '.root';
  //! 2) 'fileName' contains a list of root file names.
  // --------------------------------------------------
  TChain* FileOps::createTChain(const TString &fileName, const TString &treeName, unsigned int verbosity) {
    if( verbosity >= 1 ) std::cout << "Creating TChain" << std::endl;

    TString tree = (treeName=="" ? "DiJetTree" : treeName);
    TChain* chain = new TChain(tree); 
    
    // Option 1: single root file
    if( fileName.EndsWith(".root") ) {
      if( verbosity >= 1 ) std::cout << "  Adding '" << tree << "' from file '" << fileName << "'" << std::endl;
      chain->Add(fileName);
    }
    // Option 2: list of root files
    else {
      if( verbosity >= 1 ) std::cout << "  Opening files from list '" << fileName << "'" << std::endl;
      std::ifstream filelist;
      filelist.open(fileName);
      int nOpenedFiles = 0;
      if( filelist.is_open() ) {
	TString name = "";
	while( !filelist.eof() ) {
	  filelist >> name;
	  if( filelist.eof() ) break;
	  if( verbosity >= 1 ) std::cout << "  Adding '" << tree << "' from file '" << name << "'" << std::endl;
	  chain->Add(name);
	  nOpenedFiles++;
	}
      } else {
	std::cerr << "ERROR opening file '" << fileName << "'\n";
	exit(1);
      }
      filelist.close();
    }
    
    return chain;
  }
}
#endif

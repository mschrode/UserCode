// $Id:  $

#ifndef FileOps_h
#define FileOps_h

#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TString.h"

#include "utils.h"
#include "HistOps.h"

namespace util
{
  class FileOps
  {
  public:

    //! Read TH1 histogram from file
    // -------------------------------------------------------------------------------------
    static TH1* readTH1(const TString &fileName, const TString &histName, const TString &newHistName = "") {
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
    static util::HistVec readTH1(const std::vector<TString> &fileNames, const TString &histName, const TString &newHistName = "") {
      util::HistVec v(fileNames.size());
      for(unsigned int i = 0; i < fileNames.size(); ++i) {
	v[i] = readTH1(fileNames[i],histName,newHistName+util::toTString(i));
      }
      std::cout << "Done\n";

      return v;
    }



    // -------------------------------------------------------------------------------------
    static util::HistVec readHistVec(const TString &fileName, const TString &histName, const TString &newHistName = "") {
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

      return v;
    }
  };
}
#endif

// $Id: $

#ifndef HistOps_h
#define HistOps_h

#include <string>
#include <vector>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <cmath>

#include <TH2F.h>
#include <TObject.h>
#include <TString.h>

namespace util
{
  //!  \brief    Collection of operations on histograms
  //!
  //!  HistOps bundles some useful operations on histograms.
  //!
  //!  \note     So far all operations are static.
  //!  
  //!  \author   Matthias Schroeder (www.desy.de/~matsch)
  //!  \date     2009/03/20
  //!  $Id: $
  class HistOps
  {
  public:
    static void HistOps::DrawRatioPlot(TCanvas& can, const TH1F& h1, const TH1F& h2, TString drawOption, bool drawRatioError) { HistOps::DrawRatioPlotGeneric(can,h1,h2,drawOption,drawRatioError); } //!< Draw ratio plot of two histograms 
    static void HistOps::DrawRatioPlot(TCanvas& can, const TH1F& h, const TF1& f, TString drawOption, bool drawRatioError) { HistOps::DrawRatioPlotGeneric(can,h,f,drawOption,drawRatioError); } //!< Draw ratio plot of a histogram and a function
    static std::vector<TH1F*> HistOps::FitMean(const TH2F *h2, const std::string namePrefix, const std::string nameSuffix);
    static TH1F* GetRelDiff(const TH1F *hBase, const TH1F* hist);
    static std::vector<TH1F*> GetRelDiff(const TH1F *hBase, const std::vector<TH1F*> hists);
    static int WriteToRootFile(const std::vector<TH1F*>& hists, TString fileName);

    HistOps() {;}   //!< Default constructor
    ~HistOps() {;}  //!< Destructor


  private:
    static void DrawRatioPlotGeneric(TCanvas& can, const TH1F& h1, const TObject& h2, TString drawOption, bool drawRatioError);
    static int WriteToRootFileGeneric(const std::vector<TObject*>& obj, TString fileName);
  };
}
#endif

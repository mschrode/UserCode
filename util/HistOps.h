// $Id: HistOps.h,v 1.1 2009/05/04 14:22:16 mschrode Exp $

#ifndef HistOps_h
#define HistOps_h

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "TH1.h"
#include "TH1D.h"
#include "TString.h"
#include "TStyle.h"

#include "LabelFactory.h"

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
  //!  $Id: HistOps.h,v 1.1 2009/05/04 14:22:16 mschrode Exp $
  class HistOps
  {
  public:
    // -------------------------------------------------------------------------------------
    static TH1D *createTH1D(const TString &name, int n, double xMin, double xMax, const TString &title) {
      return new TH1D(name,title,n,xMin,xMax);
    }


    // -------------------------------------------------------------------------------------
    static TH1D *createTH1D(const TString &name, int n, double xMin, double xMax, const TString &xTitle, const TString &xUnit, const TString &yTitle, bool norm = false) {
      // create histogram without axis label
      TH1D *h = createTH1D(name,n,xMin,xMax,"");

      // x axis label
      h->SetXTitle(xUnit.Length()==0 ? xTitle : xTitle+" ("+xUnit+")");

      // y axis label
      TString yAxisTitle;
      if( norm ) yAxisTitle += "1 / N  ";
      if( yTitle == "jets" || yTitle == "events" ) {
	if( yTitle == "jets" )        yAxisTitle += "Number of Jets / ";
	else if( yTitle == "events" ) yAxisTitle += "Number of Events / ";
	yAxisTitle += doubleToTString(h->GetBinWidth(1));
	if( xUnit.Length() ) yAxisTitle += " "+xUnit;
      } else {
	yAxisTitle = yTitle;
      }
      h->SetYTitle(yAxisTitle);

      return h;
    }


    // -------------------------------------------------------------------------------------
    static void findYRange(const TH1 *h, double& min, double& max) {
      min = h->GetMinimum();
      max = h->GetMaximum();
    }


    // -------------------------------------------------------------------------------------
    static void setYRange(TH1 *h, double cMin = 1., double cMax = 1., double minLimit = 3E-8) {
      double min = 0.;
      double max = 0.;
      findYRange(h,min,max);
      min *= cMin;
      max *= cMax;
      if( min < minLimit ) min = minLimit;
      h->GetYaxis()->SetRangeUser(min,max);
    }


    // -------------------------------------------------------------------------------------
    static void setYRange(TH1 *h, int nLabelLines) {
      double min = 0.;
      double max = 0.;
      findYRange(h,min,max);
      double padHeight = 1. - gStyle->GetPadTopMargin() - gStyle->GetPadBottomMargin();
      double labelHeight = util::LabelFactory::lineHeight()*(0.5+nLabelLines);
      max = max + labelHeight/padHeight*(max-min)/(1.-labelHeight/padHeight);
      h->GetYaxis()->SetRangeUser(min,max);
    }


    // -------------------------------------------------------------------------------------
    static void normHist(TH1 *h, const TString &option = "") { 
      if( h->Integral(option) ) h->Scale(1./h->Integral(option)); 
    }


    // -------------------------------------------------------------------------------------
    static void normHist(TH1 *h, double min, double max, const TString &option = "") { 
      double norm = h->Integral(h->FindBin(min),h->FindBin(max),option);
      if( norm ) h->Scale(1./norm);
    }




    //    static createRatioPlot();


  private:
    // -------------------------------------------------------------------------------------
    static std::string toString(double d) {
      std::stringstream ss;
      ss << d;
      return ss.str();
    }


    // -------------------------------------------------------------------------------------
    static TString doubleToTString(double d) {
      return toString(d).c_str();
    }
  };
}
#endif

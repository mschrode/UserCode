// $Id: HistOps.h,v 1.2 2010/07/25 20:45:56 mschrode Exp $

#ifndef HistOps_h
#define HistOps_h

#include <cassert>
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
  //!  $Id: HistOps.h,v 1.2 2010/07/25 20:45:56 mschrode Exp $
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
      setAxisTitles(h,xTitle,xUnit,yTitle,norm);

      return h;
    }


    // -------------------------------------------------------------------------------------
    static void findYRange(const TH1 *h, double& min, double& max) {
      min = h->GetMinimum();
      max = h->GetMaximum();
    }



    // -------------------------------------------------------------------------------------
    static void findYRange(const TH1 *h, int nLabelLines, double& min, double& max) {
      findYRange(h,min,max);
      double padHeight = 1. - gStyle->GetPadTopMargin() - gStyle->GetPadBottomMargin();
      double labelHeight = util::LabelFactory::lineHeight()*(0.5+nLabelLines);
      max = max + labelHeight/padHeight*(max-min)/(1.-labelHeight/padHeight);
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
      findYRange(h,nLabelLines,min,max);
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



    // -------------------------------------------------------------------------------------
    static void *setAxisTitles(TH1 *h, const TString &xTitle, const TString &xUnit, const TString &yTitle, bool norm = false) {
      // x axis label
      h->SetXTitle(xUnit.Length()==0 ? xTitle : xTitle+" ("+xUnit+")");

      // y axis label
      TString yAxisTitle;
      if( norm ) yAxisTitle += "1 / N  ";
      if( yTitle == "jets" || yTitle == "events" ) {
	if( yTitle == "jets" )        yAxisTitle += "Jets / ";
	else if( yTitle == "events" ) yAxisTitle += "Events / ";
	yAxisTitle += doubleToTString(h->GetBinWidth(1));
	if( xUnit.Length() ) yAxisTitle += " "+xUnit;
      } else {
	yAxisTitle = yTitle;
      }
      h->SetYTitle(yAxisTitle);

      return h;
    }



    // -------------------------------------------------------------------------------------
    static TH1D *createRatioFrame(const TH1 *h, const TString &yTitle) {
      double min = 0.;
      double max = 0.;
      findYRange(h,min,max);
      return createRatioFrame(h,yTitle,min,max);
    }

    // -------------------------------------------------------------------------------------
    static TH1D *createRatioFrame(const TH1 *h, const TString &yTitle, int nLabelLines) {
      double min = 0.;
      double max = 0.;
      findYRange(h,nLabelLines,min,max);
      return createRatioFrame(h,yTitle,min,max);
    }

    // -------------------------------------------------------------------------------------
    static TH1D *createRatioFrame(const TH1 *h, const TString &yTitle, double yMin, double yMax) {
      TString name = "Frame_";
      name += h->GetName();
      TH1D *hFrame = new TH1D(name,"",10*h->GetNbinsX(),
			      h->GetXaxis()->GetBinLowEdge(1),
			      h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
      hFrame->SetLineStyle(2);
      hFrame->GetYaxis()->SetRangeUser(yMin,yMax);
      hFrame->GetYaxis()->SetTitle(yTitle);
      for(int bin = 1; bin <= hFrame->GetNbinsX(); ++bin) {
	hFrame->SetBinContent(bin,1.);
	hFrame->SetBinError(bin,0.);
      }

      return hFrame;
    }



    // -------------------------------------------------------------------------------------
    static TH1D *createRatioPlot(const TH1 *h1, const TH1 *h2, const TString &yTitle, double yMin, double yMax) {
      TH1D *hRatio = createRatioPlot(h1,h2,yTitle);
      hRatio->GetYaxis()->SetRangeUser(yMin,yMax);

      return hRatio;
    }

    // -------------------------------------------------------------------------------------
    static TH1D *createRatioPlot(const TH1 *h1, const TH1 *h2, const TString &yTitle, int nLabelLines) {
      TH1D *hRatio = createRatioPlot(h1,h2,yTitle);
      setYRange(hRatio,nLabelLines);

      return hRatio;
    }

    // -------------------------------------------------------------------------------------
    static TH1D *createRatioPlot(const TH1 *h1, const TH1 *h2, const TString &yTitle) {
      assert( h1->GetNbinsX() == h2->GetNbinsX() );
      TString name = "Ratio_";
      name += h1->GetName();
      TH1D *hRatio = static_cast<TH1D*>(h1->Clone(name));
      hRatio->SetMarkerStyle(20);
      hRatio->SetMarkerColor(h1->GetLineColor());
      hRatio->SetYTitle(yTitle);
      hRatio->Divide(h2);

      return hRatio;
    }


  private:
    // -------------------------------------------------------------------------------------
    static double round(double d, int decPlaces) {
      d *= pow(10.,1.*decPlaces);
      d = std::floor(d+0.5);
      d /= pow(10.,1.*decPlaces);

      return d;
    }



    // -------------------------------------------------------------------------------------
    static std::string toString(double d) {
      std::stringstream ss;
      ss << round(d,3);
      return ss.str();
    }


    // -------------------------------------------------------------------------------------
    static TString doubleToTString(double d) {
      return toString(d).c_str();
    }
  };
}
#endif

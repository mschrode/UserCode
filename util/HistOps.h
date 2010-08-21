// $Id: HistOps.h,v 1.10 2010/08/19 09:24:09 mschrode Exp $

#ifndef HistOps_h
#define HistOps_h

#include <cassert>
#include <cmath>
#include <vector>

#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TString.h"
#include "TStyle.h"

#include "LabelFactory.h"
#include "StyleSettings.h"
#include "utils.h"

namespace util
{
  //!  \brief    Collection of encapsulated methods on histograms
  //!
  //!  HistOps encapsulates some useful operations on histograms.
  //!  
  //!  \author   Matthias Schroeder (www.desy.de/~matsch)
  //!  \date     2009/03/20
  //!  $Id: HistOps.h,v 1.10 2010/08/19 09:24:09 mschrode Exp $
  class HistOps
  {
  public:
    // -------------------------------------------------------------------------------------
    static TH1D *createTH1D(const TString &name, int n, double xMin, double xMax, const TString &title) {
      return new TH1D(name,title,n,xMin,xMax);
    }


    // -------------------------------------------------------------------------------------
    static TH1D *createTH1D(const TString &name, int n, const double *xBinEdges, const TString &title) {
      return new TH1D(name,title,n,xBinEdges);
    }


    // -------------------------------------------------------------------------------------
    static TH1D *createTH1D(const TString &name, int n, double xMin, double xMax, const TString &xTitle, const TString &xUnit, const TString &yTitle, bool norm = false) {
      // create histogram without axis label
      TH1D *h = createTH1D(name,n,xMin,xMax,"");
      setAxisTitles(h,xTitle,xUnit,yTitle,norm);

      return h;
    }


    // -------------------------------------------------------------------------------------
    static TH1D *createTH1D(const TString &name, int n, const double *xBinEdges, const TString &xTitle, const TString &xUnit, const TString &yTitle, bool norm = false) {
      // create histogram without axis label
      TH1D *h = createTH1D(name,n,xBinEdges,"");
      setAxisTitles(h,xTitle,xUnit,yTitle,norm);

      return h;
    }


    // -------------------------------------------------------------------------------------
    static void findYRange(const TH1 *h, double& min, double& max) {
      min = 1E10;
      max = 0.;
      for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
	double val = h->GetBinContent(bin);
	if( val < min ) min = val;
	if( val > max ) max = val;
      }
      if( min > max ) {
	min = 1E-3;
	max = 1;
      }
    }



    // -------------------------------------------------------------------------------------
    static void findYRange(const TH1 *h, int nLabelLines, double& min, double& max, bool log = false) {
      findYRange(h,min,max);
      double padHeight = 1. - gStyle->GetPadTopMargin() - gStyle->GetPadBottomMargin();
      double labelHeight = util::LabelFactory::lineHeight()*(1+nLabelLines) + util::LabelFactory::labelTopOffset();
      if( log ) {
	if( min <= 0. ) min = 3E-10;
	max = exp((1.-labelHeight)*(log10(max/pow(min,labelHeight))));
      } else {
	max = max + labelHeight/padHeight*(max-min)/(1.-labelHeight/padHeight);
      }
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
    static void setYRange(TH1 *h, int nLabelLines, bool log = false) {
      double min = 0.;
      double max = 0.;
      findYRange(h,nLabelLines,min,max,log);
      h->GetYaxis()->SetRangeUser(min,max);
    }


    // -------------------------------------------------------------------------------------
    static void setStyleColor(TH1 *h, int styleColor) {
      setColor(h,util::StyleSettings::color(styleColor));
    }


    // -------------------------------------------------------------------------------------
    static void setColor(TH1 *h, int color) {
      h->SetMarkerColor(color);
      h->SetLineColor(color);
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
    static void setAxisTitles(TH1 *h, const TString &xTitle, const TString &xUnit, const TString &yTitle, bool norm = false) {
      // x axis label
      h->SetXTitle(xUnit.Length()==0 ? xTitle : xTitle+" ("+xUnit+")");

      // y axis label
      TString yAxisTitle;
      if( norm ) yAxisTitle += "1 / N  ";
      if( yTitle == "jets" || yTitle == "events" ) {
	if( yTitle == "jets" )        yAxisTitle += "Jets / ";
	else if( yTitle == "events" ) yAxisTitle += "Events / ";
	if( norm ) yAxisTitle += util::toTString(h->GetBinWidth(1),3);
	if( xUnit.Length() ) yAxisTitle += " "+xUnit;
      } else {
	yAxisTitle = yTitle;
      }
      h->SetYTitle(yAxisTitle);
    }


    // -------------------------------------------------------------------------------------
    static TH1D *createRatioFrame(double xMin, double xMax, double yMin, double yMax, const TString &xTitle, const TString &yTitle) {
      TString name = "util::HistOps::hRatioFrame";
      name += nFrames_;
      ++nFrames_;
      TH1D *h = new TH1D(name,"",1000,xMin,xMax);
      for(int xBin = 1; xBin <= h->GetNbinsX(); ++xBin) {
	h->SetBinContent(xBin,1.);
      }
      h->GetXaxis()->SetTitle(xTitle);
      h->GetYaxis()->SetTitle(yTitle);
      h->GetYaxis()->SetRangeUser(yMin,yMax);
      h->SetLineStyle(2);
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
      TH1D *hFrame = 0;
      bool hasConstBinWidth = true;
      for(int bin = 2; bin <= h->GetNbinsX(); ++bin) {
	if( h->GetBinWidth(bin) != h->GetBinWidth(1) ) {
	  hasConstBinWidth = false;
	  break;
	}
      }
      if( hasConstBinWidth ) {
	hFrame = new TH1D(name,"",10*h->GetNbinsX(),
			  h->GetXaxis()->GetBinLowEdge(1),
			  h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
      } else {
	hFrame = new TH1D(name,"",h->GetNbinsX(),
			  h->GetXaxis()->GetXbins()->GetArray());
      }
      hFrame->SetLineStyle(2);
      hFrame->GetYaxis()->SetRangeUser(yMin,yMax);
      hFrame->GetYaxis()->SetTitle(yTitle);
      hFrame->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
      for(int bin = 1; bin <= hFrame->GetNbinsX(); ++bin) {
	hFrame->SetBinContent(bin,1.);
	hFrame->SetBinError(bin,0.);
      }

      return hFrame;
    }


    // -------------------------------------------------------------------------------------
    static TH1D *createFrame(const TH1 *h, const TString &yTitle, double yMin, double yMax) {
      TString name = "Frame_";
      name += h->GetName();
      TH1D *hFrame = new TH1D(name,"",10*h->GetNbinsX(),
			      h->GetXaxis()->GetBinLowEdge(1),
			      h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
      hFrame->GetYaxis()->SetRangeUser(yMin,yMax);
      hFrame->GetYaxis()->SetTitle(yTitle);
      hFrame->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());

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
    static TH1D *createRatioPlot(const TH1 *h1, const TH1 *h2, const TString &yTitle = "") {
      assert( h1->GetNbinsX() == h2->GetNbinsX() );
      TString name = "Ratio_";
      name += h1->GetName();
      TH1D *hRatio = static_cast<TH1D*>(h1->Clone(name));
      hRatio->SetMarkerStyle(20);
      hRatio->SetMarkerColor(h1->GetLineColor());
      hRatio->SetYTitle(yTitle);
      hRatio->Divide(h2);
      setYRange(hRatio,1.,1.,0.);

      return hRatio;
    }

    // -------------------------------------------------------------------------------------
    static void fillSlices(const TH2 *h2, std::vector<TH1*> &hSlices, const TString &namePrefix) {
      for(int xBin = 1; xBin <= h2->GetNbinsX(); ++xBin) {
	TH1 *h = new TH1D(namePrefix+toTString(xBin-1),"",h2->GetNbinsY(),
			  h2->GetYaxis()->GetXmin(),h2->GetYaxis()->GetXmax());
	h->Sumw2();
	h->SetXTitle(h2->GetYaxis()->GetTitle());
	for(int yBin = 1; yBin <= h2->GetNbinsY(); yBin++) {
	  h->SetBinContent(yBin,h2->GetBinContent(h2->GetBin(xBin,yBin)));
	  h->SetBinError(yBin,h2->GetBinError(h2->GetBin(xBin,yBin)));
	}  
	hSlices.push_back(h);
      }
    }

    // -------------------------------------------------------------------------------------
    static bool equidistLogBins(std::vector<double>& binEdges, unsigned int nBins, double min, double max) {
      if( binEdges.size() != nBins+1 ) return false;
      if( min <= 0. || max <= 0. || min >= max ) return false;
      
      binEdges[0]     = min;
      binEdges[nBins] = max;
      const double minLog = log10(binEdges[0]);
      const double maxLog  = log10(binEdges[nBins]);
      for (unsigned int i = 1; i < nBins; ++i) {
	binEdges[i] = pow(10., minLog + i*(maxLog-minLog)/(nBins));
      }
      
      return true;
    }
 

  private:
    static unsigned int nFrames_;
  };

  unsigned int HistOps::nFrames_ = 0;
}
#endif

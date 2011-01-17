// $Id: HistOps.h,v 1.28 2010/12/27 13:31:17 mschrode Exp $

#ifndef HistOps_h
#define HistOps_h

#include <cassert>
#include <cmath>
#include <vector>

#include "TAttPad.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

#include "LabelFactory.h"
#include "StyleSettings.h"
#include "utils.h"

namespace util
{
  typedef std::vector<TH1*> HistVec;
  typedef std::vector<TH1*>::iterator HistIt;
  typedef std::vector<TH1*>::const_iterator HistItConst;


  //!  \brief    Collection of encapsulated methods on histograms
  //!
  //!  HistOps encapsulates some useful operations on histograms.
  //!  
  //!  \author   Matthias Schroeder (www.desy.de/~matsch)
  //!  \date     2009/03/20
  //!  $Id: HistOps.h,v 1.28 2010/12/27 13:31:17 mschrode Exp $
  class HistOps
  {
  public:
    // -------------------------------------------------------------------------------------
    static TH1D *createTH1D(const TString &name, int n, double xMin, double xMax, const TString &title) {
      return new TH1D(name,title,n,xMin,xMax);
    }


    // -------------------------------------------------------------------------------------
    static TH1D *createTH1D(const TString &name, const std::vector<double> &xBinEdges, const TString &title) {
      return new TH1D(name,title,xBinEdges.size()-1,&(xBinEdges.front()));
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
    static TH1D *createTH1D(const TString &name, const std::vector<double> &xBinEdges, const TString &xTitle, const TString &xUnit, const TString &yTitle, bool norm = false) {
      // create histogram without axis label
      TH1D *h = createTH1D(name,xBinEdges,"");
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
    static void findXRange(const TH1 *h, int& min, int& max) {
      min = 1;
      for(int i = 1; i <= h->GetNbinsX(); ++i) {
	if( h->GetBinContent(i) > 0. ) {
	  min = i;
	  break;
	}
      }
      max = h->GetNbinsX();
      for(int i = h->GetNbinsX(); i > 0; --i) {
	if( h->GetBinContent(i) > 0. ) {
	  max = i;
	  break;
	}
      }
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
    static void findYRange(const TH1 *h, int nLabelLines, double& min, double& max, double logMin = -1.) {
      findYRange(h,min,max);
      double padHeight = 1. - gStyle->GetPadTopMargin() - gStyle->GetPadBottomMargin();
      double labelHeight = util::LabelFactory::lineHeight()*(1+nLabelLines) + util::LabelFactory::labelTopOffset();
      if( logMin > 0. ) {
	min = logMin;
	max = pow((log10(max) - log10(min)*0.7*labelHeight/padHeight)/(1.-0.7*labelHeight/padHeight),10.);
      } else {
	max = (max - min*labelHeight/padHeight)/(1.-labelHeight/padHeight);
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
    static void setYRange(TH1 *h, size_t nLabelLines, double logMin = -1.) {
      return setYRange(h,static_cast<int>(nLabelLines),logMin);
    }


    // -------------------------------------------------------------------------------------
    static void setYRange(TH1 *h, int nLabelLines, double logMin = -1.) {
      double min = 0.;
      double max = 0.;
      findYRange(h,nLabelLines,min,max,logMin);
      h->GetYaxis()->SetRangeUser(min,max);
    }


    // -------------------------------------------------------------------------------------
    static void setYRange(HistVec &hists, int nLabelLines, double logMin = -1.) {
      double min = 0.;
      double max = 0.;
      for(HistIt it = hists.begin(); it != hists.end(); ++it) {
	double tmpMin = 0.;
	double tmpMax = 0.;
	findYRange(*it,nLabelLines,tmpMin,tmpMax,logMin);
	if( tmpMin < min ) min = tmpMin;
	if( tmpMax > max ) max = tmpMax;
      }
      for(HistIt it = hists.begin(); it != hists.end(); ++it) {
	(*it)->GetYaxis()->SetRangeUser(min,max);
      }

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
    static void getBinBorders(const TH1* h, double &min, double &max, int bin = -1) { 
      min = h->GetXaxis()->GetBinLowEdge(( bin>0 ? bin : 1 ));
      max = h->GetXaxis()->GetBinUpEdge(( bin>0 ? bin : h->GetNbinsX() ));
    }


    // -------------------------------------------------------------------------------------
    static void setAxisTitles(TH1 *h, const TString &xTitle, const TString &xUnit, const TString &yTitle, bool norm = false) {
      // x axis label
      h->SetXTitle(xUnit.Length()==0 ? xTitle : xTitle+" ("+xUnit+")");

      // y axis label
      TString yAxisTitle;
      if( norm ) yAxisTitle += "1 / N  ";
      if( yTitle == "jets" || yTitle == "events" || yTitle == "Jets" || yTitle == "Events" ) {
	if( yTitle == "jets" || yTitle == "Jets" )        yAxisTitle += "Jets";
	else if( yTitle == "events" || yTitle == "Events" ) yAxisTitle += "Events";
	if( norm ) {
	  yAxisTitle += (" / "+util::toTString(h->GetBinWidth(1),3));
	  if( xUnit.Length() ) yAxisTitle += " "+xUnit;
	}
      } else {
	yAxisTitle = yTitle;
      }
      h->SetYTitle(yAxisTitle);
    }


    // -------------------------------------------------------------------------------------
    static TGraphAsymmErrors* combineTGraphs(const TGraphAsymmErrors* g1, const TGraphAsymmErrors* g2) {
      std::vector<double> x;
      std::vector<double> xeh;
      std::vector<double> xel;
      std::vector<double> y;
      std::vector<double> yeh;
      std::vector<double> yel;
      for(int i = 0; i < g1->GetN(); ++i) {
	x.push_back(g1->GetX()[i]);
	xeh.push_back(g1->GetEXhigh()[i]);
	xel.push_back(g1->GetEXlow()[i]);
	y.push_back(g1->GetY()[i]);
	yeh.push_back(g1->GetEYhigh()[i]);
	yel.push_back(g1->GetEYlow()[i]);
      }
      for(int i = 0; i < g2->GetN(); ++i) {
	x.push_back(g2->GetX()[i]);
	xeh.push_back(g2->GetEXhigh()[i]);
	xel.push_back(g2->GetEXlow()[i]);
	y.push_back(g2->GetY()[i]);
	yeh.push_back(g2->GetEYhigh()[i]);
	yel.push_back(g2->GetEYlow()[i]);
      }
      return new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),&(xel.front()),&(xeh.front()),&(yel.front()),&(yeh.front()));
    }



    // -------------------------------------------------------------------------------------
    static TH1D *createRatioFrame(double xMin, double xMax, double yMin, double yMax, const TString &xTitle, const TString &yTitle) {
      ++COUNT_;
      TH1D *h = new TH1D("util::HistOps::hRatioFrame"+toTString(COUNT_),"",1000,xMin,xMax);
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
      ++COUNT_;
      TString name = "Frame_";
      name += h->GetName();
      name += COUNT_;
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
      hRatio->SetMarkerStyle(h1->GetMarkerStyle());
      hRatio->SetMarkerColor(h1->GetLineColor());
      hRatio->SetYTitle(yTitle);
      hRatio->Divide(h2);
      setYRange(hRatio,1.,1.,0.);

      return hRatio;
    }

    // -------------------------------------------------------------------------------------
    static TH1D *createRatioPlot(const TH1 *h, const TF1 *f, const TString &yTitle = "") {
      TString name = "Ratio_";
      name += h->GetName();
      TH1D *hRatio = static_cast<TH1D*>(h->Clone(name));
      hRatio->SetMarkerStyle(h->GetMarkerStyle());
      hRatio->SetMarkerColor(h->GetLineColor());
      hRatio->SetYTitle(yTitle);
      hRatio->Reset();
      for(int bin = 1; bin <= hRatio->GetNbinsX(); ++bin) {
	double ref = f->Eval(hRatio->GetBinCenter(bin));
	if( ref != 0. ) {
	  hRatio->SetBinContent(bin,h->GetBinContent(bin)/ref);
	  hRatio->SetBinError(bin,h->GetBinError(bin)/ref);
	}
      }
      setYRange(hRatio,1.,1.,0.);

      return hRatio;
    }

    // -------------------------------------------------------------------------------------
    static TGraphAsymmErrors *createRatioGraph(const TGraphAsymmErrors *g, const TF1 *f) {
      TGraphAsymmErrors *gRatio = new TGraphAsymmErrors(g->GetN());
      for(int i = 0; i < g->GetN(); ++i) {
	double yTrue = f->Eval(g->GetX()[i]);
	if( yTrue != 0. ) {
	  gRatio->SetPoint(i,g->GetX()[i],g->GetY()[i]/yTrue);
	  gRatio->SetPointError(i,g->GetEXhigh()[i],g->GetEXlow()[i],g->GetEYhigh()[i]/yTrue,g->GetEYlow()[i]/yTrue);
	} else {
	  gRatio->SetPoint(i,g->GetX()[i],0.);
	  gRatio->SetPointError(i,g->GetEXhigh()[i],g->GetEXlow()[i],0.,0.);
	}
      }
      gRatio->SetMarkerStyle(g->GetMarkerStyle());
      gRatio->SetMarkerColor(g->GetMarkerColor());
      gRatio->SetLineColor(g->GetLineColor());
      
      return gRatio;
    }


    // -------------------------------------------------------------------------------------
    static TGraphErrors *createRatioGraph(const TGraphErrors *g, const TF1 *f) {
      TGraphErrors *gRatio = new TGraphErrors(g->GetN());
      for(int i = 0; i < g->GetN(); ++i) {
	double yTrue = f->Eval(g->GetX()[i]);
	if( yTrue != 0. ) {
	  gRatio->SetPoint(i,g->GetX()[i],g->GetY()[i]/yTrue);
	  gRatio->SetPointError(i,g->GetEX()[i],g->GetEY()[i]/yTrue);
	} else {
	  gRatio->SetPoint(i,g->GetX()[i],0.);
	  gRatio->SetPointError(i,g->GetEXhigh()[i],0.);
	}
      }
      gRatio->SetMarkerStyle(g->GetMarkerStyle());
      gRatio->SetMarkerColor(g->GetMarkerColor());
      gRatio->SetLineColor(g->GetLineColor());
      
      return gRatio;
    }


    // -------------------------------------------------------------------------------------
    static TGraphErrors *createRatioGraph(const TGraphErrors *g1, const TGraphErrors *g2) {
      TGraphErrors *gRatio = new TGraphErrors(g1->GetN());
      for(int i = 0; i < g1->GetN() && i < g2->GetN(); ++i) {
	double denom = g2->GetY()[i];
	double x = 0.5*(g1->GetX()[i]+g2->GetX()[i]); // Mean value of x1 and x2
	double xe = sqrt(g1->GetEX()[i]*g1->GetEX()[i]+g2->GetEX()[i]*g2->GetEX()[i]); // Quadratic mean of x1e and x2e
	if( denom != 0. ) {
	  double y = g1->GetY()[i]/denom;
	  double ye = sqrt(pow(g1->GetEY()[i]/denom,2.)+pow(g1->GetY()[i]*g2->GetEY()[i]/denom/denom,2.));
	  gRatio->SetPoint(i,x,y);
	  gRatio->SetPointError(i,xe,ye);
	} else {
	  gRatio->SetPoint(i,x,0.);
	  gRatio->SetPointError(i,xe,0.);
	}
      }
      gRatio->SetMarkerStyle(g1->GetMarkerStyle());
      gRatio->SetMarkerColor(g1->GetMarkerColor());
      gRatio->SetLineColor(g1->GetLineColor());
      
      return gRatio;
    }


    // -------------------------------------------------------------------------------------
    static TGraphAsymmErrors *createRatioGraph(const TGraphAsymmErrors *g1, const TGraphAsymmErrors *g2) {
      TGraphAsymmErrors *gRatio = new TGraphAsymmErrors(g1->GetN());
      for(int i = 0; i < g1->GetN() && i < g2->GetN(); ++i) {
	double denom = g2->GetY()[i];
	double x = 0.5*(g1->GetX()[i]+g2->GetX()[i]); // Mean value of x1 and x2
	double e1 = g1->GetEXhigh()[i];
	double e2 = g2->GetEXhigh()[i];
	double xe = sqrt( e1*e1 + e2*e2 ); // Quadratic mean of x1e and x2e
	if( denom != 0. ) {
	  double y = g1->GetY()[i]/denom;
	  e1 = g1->GetEYhigh()[i];
	  e2 = g2->GetEYhigh()[i];
	  double ye = sqrt(pow(e1/denom,2.)+pow(g1->GetY()[i]*e2/denom/denom,2.));
	  gRatio->SetPoint(i,x,y);
	  gRatio->SetPointError(i,xe,xe,ye,ye);
	} else {
	  gRatio->SetPoint(i,x,0.);
	  gRatio->SetPointError(i,xe,0.,0.);
	}
      }
      gRatio->SetMarkerStyle(g1->GetMarkerStyle());
      gRatio->SetMarkerColor(g1->GetMarkerColor());
      gRatio->SetLineColor(g1->GetLineColor());
      
      return gRatio;
    }


    // -------------------------------------------------------------------------------------
    static TH1D *createRatio(const TF1 *f1, const TF1 *f2, double xMin, double xMax, const TString &xTitle, const TString &yTitle) {
      TString name = f1->GetName();
      name += "Ratio";
      TH1D *h = new TH1D(name,"",1000,xMin,xMax);
      h->SetXTitle(xTitle);
      h->SetYTitle(yTitle);
      for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
	if( f2->Eval(h->GetBinCenter(bin)) ) {
	  h->SetBinContent(bin,f1->Eval(h->GetBinCenter(bin))/f2->Eval(h->GetBinCenter(bin)));
	  h->SetBinError(bin,0.);
	}
      }
      h->SetLineWidth(f1->GetLineWidth());
      h->SetLineColor(f1->GetLineColor());
      return h;
    }

    // -------------------------------------------------------------------------------------
    static TCanvas *createRatioTopCanvas() {
      COUNT_++;      
      TCanvas *topCan = new TCanvas("TopRatio"+toTString(COUNT_),"",500,500);
      topCan->SetBottomMargin(0.2 + 0.8*topCan->GetBottomMargin()-0.2*topCan->GetTopMargin());
      return topCan;
    }

    // -------------------------------------------------------------------------------------
    static TPad *createRatioBottomPad() {
      COUNT_++;
      TPad *ratioPad = new TPad("BottomPad"+toTString(COUNT_),"",0,0,1,1);
      ratioPad->SetTopMargin(0.8 - 0.8*ratioPad->GetBottomMargin()+0.2*ratioPad->GetTopMargin());
      ratioPad->SetFillStyle(0);
      ratioPad->SetFrameFillColor(10);
      ratioPad->SetFrameBorderMode(0);
      return ratioPad;
    }

    // -------------------------------------------------------------------------------------
    static TH1 *createRatioTopFrame(const TH1 *h) {
      TH1 *hFrame = createRatioTopHist(h);
      hFrame->Reset();
      return hFrame;
    }

    // -------------------------------------------------------------------------------------
    static TH1 *createRatioTopFrame(double xMin, double xMax, double yMin, double yMax, const TString &yTitle) {
      COUNT_++;
      TH1D *h = new TH1D("util::HistOps::hRatioTopFrame"+toTString(COUNT_),"",1000,xMin,xMax);
      setAxisTitles(h,"","",yTitle);
      h->GetYaxis()->SetRangeUser(yMin,yMax);
      h->GetXaxis()->SetLabelSize(0);
      h->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/0.8);
      return h;
    }

    // -------------------------------------------------------------------------------------
    static TH1 *createRatioTopHist(const TH1 *h) {
      COUNT_++;
      TH1D *hFrame = static_cast<TH1D*>(h->Clone("util::HistOps::hRatioTopHist"+toTString(COUNT_)));
      double min = 1.;
      double max = 0.;
      findYRange(h,min,max);
      hFrame->GetYaxis()->SetRangeUser(min,max);
      hFrame->GetXaxis()->SetTitle("");
      hFrame->GetXaxis()->SetLabelSize(0);
      hFrame->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/0.8);
      return hFrame;
    }

    // -------------------------------------------------------------------------------------
    static TH1 *createRatioBottomFrame(const TH1 *h, const TString &xTitle, const TString &xUnit, double yMin = 0.81, double yMax = 1.19) {
      COUNT_++;
      TH1D *hFrame = createRatioFrame(h,"",yMin,yMax);
      setAxisTitles(hFrame,xTitle,xUnit,"");
      hFrame->GetYaxis()->SetNdivisions(205);
      hFrame->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/0.2);
      return hFrame;
    }


    // -------------------------------------------------------------------------------------
    static void fillSlices(const TH2 *h2, std::vector<TH1*> &hSlices, const TString &namePrefix) {
      std::vector<double> entries;
      fillSlices(h2,hSlices,namePrefix,entries);
      entries.clear();
    }


    // -------------------------------------------------------------------------------------
    static void fillSlices(const TH2 *h2, std::vector<TH1*> &hSlices, const TString &namePrefix, std::vector<double> &entries) {
      entries.clear();
      for(int xBin = 1; xBin <= h2->GetNbinsX(); ++xBin) {
	TH1 *h = new TH1D(namePrefix+toTString(xBin-1),"",h2->GetNbinsY(),
			  h2->GetYaxis()->GetXmin(),h2->GetYaxis()->GetXmax());
	h->Sumw2();
	h->SetXTitle(h2->GetYaxis()->GetTitle());
	for(int yBin = 1; yBin <= h2->GetNbinsY(); yBin++) {
	  h->SetBinContent(yBin,h2->GetBinContent(h2->GetBin(xBin,yBin)));
	  h->SetBinError(yBin,h2->GetBinError(h2->GetBin(xBin,yBin)));
	}  
	entries.push_back(h2->Integral(xBin,xBin,1,h2->GetNbinsY()));
	hSlices.push_back(h);
      }
    }


    // --------------------------------------------------
    static bool fitCoreWidth(const TH1* hist, double nSig, double &width, double &widthErr) {
      double rms = 0.;
      double rmsErr = 0.;
      bool result = fitCoreWidth(hist,nSig,width,widthErr,rms,rmsErr);
      return result;
    }
    
    
    // --------------------------------------------------
    static bool fitCoreWidth(const TH1* hist, double nSig, double &width, double &widthErr, double &rms, double &rmsErr) {
      TF1* dummy = 0;
      bool result = fitCoreWidth(hist,nSig,dummy,width,widthErr,rms,rmsErr);
      delete dummy;
      return result;
    }


    // --------------------------------------------------
    static bool fitCoreWidth(const TH1* hist, double nSig, TF1* &gauss, double &width, double &widthErr) {
      double rms = 0.;
      double rmsErr = 0.;
      bool result = fitCoreWidth(hist,nSig,gauss,width,widthErr,rms,rmsErr);
      return result;
    }


    // --------------------------------------------------
    static bool fitCoreWidth(const TH1* hist, double nSig, TF1* &gauss, double &width, double &widthErr, double &rms, double &rmsErr) {
      bool result  = false;
      
      TString name = hist->GetName();
      name += "_GaussFit";
      gauss = new TF1(name,"gaus",hist->GetXaxis()->GetBinLowEdge(1),hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX()));
      gauss->SetLineWidth(1);
      gauss->SetLineColor(kRed);
      
      TH1* h = static_cast<TH1*>(hist->Clone("util::HistOps::fitCoreWidth::h"));
      
      double mean = h->GetMean();
      rms = h->GetRMS();
      rmsErr = h->GetRMSError();
      double sig = 2.5*rms;
      if( h->Fit(gauss,"0QIR","",mean-sig,mean+sig) == 0 ) {
	mean = gauss->GetParameter(1);
	sig = nSig*gauss->GetParameter(2);
	if( h->Fit(gauss,"0QIR","",mean-sig,mean+sig) == 0 ) {
	  result = true;
	  width = gauss->GetParameter(2);
	  widthErr = gauss->GetParError(2);
	}
      } else {
	std::cerr << "WARNING in util::HistOps::fitCoreWidth: No convergence when fitting width of '" << h->GetName() << "'\n";
	width = 0.;
	widthErr = 10000.;
      }
      delete h;
      
      return result;
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

    // -------------------------------------------------------------------------------------
    static void setMarginsColz(TAttPad *pad) {
      pad->SetRightMargin(0.1+pad->GetRightMargin());
      pad->SetBottomMargin(0.1+pad->GetBottomMargin());
    }


    // -------------------------------------------------------------------------------------
    static TH1D *getUncertaintyBand(const TF1 *mean, double relUncert, double xMin, double xMax, int color = -1, int fillStyle = -1) {
      TString name = mean->GetName();
      name += toTString(COUNT_);
      TH1D *hBand = new TH1D(name,"",1000,xMin,xMax);
      hBand->SetMarkerStyle(1);
      hBand->SetFillColor((color < 0) ? mean->GetLineColor()+1 : color);
      if( fillStyle > 0 ) hBand->SetFillStyle(fillStyle);
      hBand->SetMarkerColor(hBand->GetFillColor());
      for(int bin = 1; bin <= hBand->GetNbinsX(); ++bin) {
	double val = mean->Eval(hBand->GetBinCenter(bin));
	hBand->SetBinContent(bin,val);
	hBand->SetBinError(bin,relUncert*val);
      }

      return hBand;
    }


    // -------------------------------------------------------------------------------------
    static TH1D *getUncertaintyBand(const TH1 *mean, double relUncert, double xMin, double xMax, int color = -1, int fillStyle = -1) {
      TString name = mean->GetName();
      name += toTString(COUNT_);
      TH1D *hBand = new TH1D(name,"",1000,xMin,xMax);
      hBand->SetMarkerStyle(1);
      hBand->SetFillColor((color < 0) ? mean->GetLineColor()+1 : color);
      if( fillStyle > 0 ) hBand->SetFillStyle(fillStyle);
      hBand->SetMarkerColor(hBand->GetFillColor());
      for(int bin = 1; bin <= hBand->GetNbinsX(); ++bin) {
	double val = mean->GetBinContent(bin);
	hBand->SetBinContent(bin,val);
	hBand->SetBinError(bin,relUncert*val);
      }

      return hBand;
    }


    // -------------------------------------------------------------------------------------
    static TH1 *getUncertaintyBandSym(const TH1* hNominal, const TH1* hUp, const TH1* hDown, int color = -1, int fillStyle = -1) {
      // Clone nominal hist
      TString name = hNominal->GetName();
      name += "Uncertainty";
      name += toTString(COUNT_);
      TH1 *hBand = static_cast<TH1*>(hNominal->Clone(name));

      if( hUp->GetNbinsX() == hBand->GetNbinsX() && hDown->GetNbinsX() == hBand->GetNbinsX() ) {
	// Set symmetrized entry and uncertainties
	for(int bin = 1; bin <= hBand->GetNbinsX(); ++bin) {
	  double y1 = hUp->GetBinContent(bin);
	  double y2 = hDown->GetBinContent(bin);
	  double e = 0.5*(y1-y2);
	  double yn = e > 0 ? y2 + 0.5*e : y1 - 0.5*e;
	  hBand->SetBinContent(bin,yn);
	  hBand->SetBinError(bin,std::abs(e));
	}
	hBand->SetMarkerStyle(0);
	hBand->SetFillColor((color < 0) ? hBand->GetLineColor()+50 : color);
	if( fillStyle > 0 ) hBand->SetFillStyle(fillStyle);
	hBand->SetMarkerColor(hBand->GetFillColor());
      } else {
	std::cerr << "WARNING in util::HistOps::getUncertaintyBand(): different binning of uncertainty and nominal histograms!" << std::endl;
      }
      
      return hBand;
    }


    // -------------------------------------------------------------------------------------
    static TGraphAsymmErrors *getUncertaintyBand(const TH1* hNominal, const TH1* hUp, const TH1* hDown, int color = -1, int fillStyle = -1) {
      TGraphAsymmErrors* g = 0;
      if( hUp->GetNbinsX() == hDown->GetNbinsX() ) {
	std::vector<double> x;
	std::vector<double> xe;
	std::vector<double> y;
	std::vector<double> yel;
	std::vector<double> yeh;

	// Set symmetrized entry and uncertainties
	for(int bin = 1; bin <= hNominal->GetNbinsX(); ++bin) {
	  x.push_back(hNominal->GetBinCenter(bin));
	  xe.push_back(0.5*hNominal->GetBinWidth(bin));
	  double nom = hNominal->GetBinContent(bin);
	  double y1 = std::abs(hDown->GetBinContent(bin)-nom);
	  double y2 = std::abs(hUp->GetBinContent(bin)-nom);
	  y.push_back(nom);
	  yel.push_back(y1);
	  yeh.push_back(y2);
	}
	g = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),&(xe.front()),&(xe.front()),&(yel.front()),&(yeh.front()));
	
	g->SetMarkerStyle(0);
	g->SetFillColor((color < 0) ? g->GetLineColor()+50 : color);
	if( fillStyle > 0 ) g->SetFillStyle(fillStyle);
	g->SetMarkerColor(g->GetFillColor());
      } else {
	std::cerr << "WARNING in util::HistOps::getUncertaintyBand(): different binning of uncertainty histograms!" << std::endl;
      }
      
      return g;
    }


    // -------------------------------------------------------------------------------------
    static TH1 *getUncertaintyBand(const TH1* hNominal, const std::vector<TH1*> &hUp, const std::vector<TH1*> &hDown, std::vector<TH1*> &hBands, int color = -1, int fillStyle = -1) {

      // Total uncertainty band
      TH1* hTotalUp = 0;
      TH1* hTotalDown = 0;
      getTotalUncertainty(hUp,hDown,hTotalUp,hTotalDown);
      
      // Individual uncertainty bands
      for(unsigned int i = 0; i < hUp.size(); ++i) {
	hBands.push_back(getUncertaintyBandSym(hNominal,hUp[i],hDown[i],color,fillStyle));
      }

      return getUncertaintyBandSym(hNominal,hTotalUp,hTotalDown,color,fillStyle);
    }


    // -------------------------------------------------------------------------------------
    static TGraphAsymmErrors *getTotalUncertainty(const std::vector<TGraphAsymmErrors*> &uncerts, int color = -1, int fillStyle = -1) {
      TGraphAsymmErrors* g = 0;
      bool isSane = true;
      for(unsigned int i = 1; i < uncerts.size(); ++i) {
	if( uncerts[i]->GetN() != uncerts[0]->GetN() ) {
	  std::cerr << "ERROR in util::HistOps::getTotalUncertainty(): different number of points in uncertatinties" << std::endl;
	  isSane = false;
	}
      }
      if( isSane ) {
	std::vector<double> x;
	std::vector<double> xel;
	std::vector<double> xeh;
	std::vector<double> y;
	std::vector<double> yel;
	std::vector<double> yeh;
	for(int i = 0; i < uncerts[0]->GetN(); ++i) {
	  double u2 = 0.;
	  double d2 = 0.;
	  for(unsigned int k = 0; k < uncerts.size(); ++k) {
	    d2 += pow(uncerts[k]->GetEYlow()[i],2.);
	    u2 += pow(uncerts[k]->GetEYhigh()[i],2.);
	  }
	  x.push_back(uncerts[0]->GetX()[i]);
	  xel.push_back(uncerts[0]->GetEXlow()[i]);
	  xeh.push_back(uncerts[0]->GetEXhigh()[i]);
	  y.push_back(uncerts[0]->GetY()[i]);
	  yel.push_back(sqrt(d2));
	  yeh.push_back(sqrt(u2));
	}
	g = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),&(xel.front()),&(xeh.front()),&(yel.front()),&(yeh.front()));	
	g->SetMarkerStyle(0);
	g->SetFillColor((color < 0) ? g->GetLineColor()+50 : color);
	if( fillStyle > 0 ) g->SetFillStyle(fillStyle);
	g->SetMarkerColor(g->GetFillColor());
      }

      return g;
    }



    // -------------------------------------------------------------------------------------
    static void getTotalUncertainty(const std::vector<TH1*> &hUp, const std::vector<TH1*> &hDown, TH1* &hTotalUp, TH1* &hTotalDown) {
      if( hUp.size() == hDown.size() ) {
	bool sameBinning = true;
	for(unsigned int i = 0; i < hUp.size(); ++i) {
	  if( hUp[i]->GetNbinsX() != hDown[i]->GetNbinsX() ) {
	    sameBinning = false;
	    break;
	  }
	}
	if( sameBinning ) {
	  // Create total uncertainties
	  hTotalUp = static_cast<TH1*>(hUp.at(0)->Clone("util::HistOps::TotalUncertainty"+toTString(COUNT_)));
	  hTotalDown = static_cast<TH1*>(hDown.at(0)->Clone("util::HistOps::TotalUncertainty"+toTString(COUNT_)));

	  // Add up individual uncertainties in quadrature
	  for(int bin = 1; bin <= hTotalUp->GetNbinsX(); ++bin) {
	    double u2 = 0.;
	    double d2 = 0.;
	    for(unsigned int i = 0; i < hUp.size(); ++i) {
	      u2 += pow(hUp[i]->GetBinContent(bin),2);
	      d2 += pow(hDown[i]->GetBinContent(bin),2);
	    }
	    hTotalUp->SetBinContent(bin,sqrt(u2));
	    hTotalDown->SetBinContent(bin,sqrt(d2));
	  }
	} else {
	  std::cerr << "WARNING in util::HistOps::getTotalUncertainty(): different binning of of up and down variations." << std::endl;
	  exit(1);
	}
      } else {
	std::cerr << "WARNING in util::HistOps::getTotalUncertainty(): different number of up and down variations." << std::endl;
	exit(1);
      }
    }



    // -------------------------------------------------------------------------------------
    static TH1* combineBins(const TH1* hOrig, const std::vector<int> &nCombBins, int startBin) {

      // Edges of combined bins
      std::vector<double> newBinEdges;
      newBinEdges.push_back(hOrig->GetXaxis()->GetBinLowEdge(startBin));
      int bin = startBin-1;
      for(unsigned int i = 0; i < nCombBins.size(); ++i) {
	bin += nCombBins.at(i);
	if( bin > hOrig->GetNbinsX() ) break;
	newBinEdges.push_back(hOrig->GetXaxis()->GetBinUpEdge(bin));
      }
      TString name = hOrig->GetName();
      TH1* hNew = new TH1D(name+"CombinedBins",hOrig->GetTitle(),newBinEdges.size()-1,&(newBinEdges.front()));
      hNew->SetXTitle(hOrig->GetXaxis()->GetTitle());
      hNew->SetYTitle(hOrig->GetYaxis()->GetTitle());
      hNew->SetMarkerStyle(hOrig->GetMarkerStyle());
      hNew->SetMarkerColor(hOrig->GetMarkerColor());
      hNew->SetLineColor(hOrig->GetLineColor());


      // Combine bin content using error weighted mean
      bin = 1;
      double y = 0.;
      double ye = 0.;
      for(int binOrig = startBin; binOrig <= hOrig->GetNbinsX(); ++binOrig) {
	double yorig = hOrig->GetBinContent(binOrig);
	double yeorig = hOrig->GetBinError(binOrig);

	if( yeorig ) {
	  y  += yorig/yeorig/yeorig;
	  ye += 1./yeorig/yeorig;
	}
	
	if( hOrig->GetXaxis()->GetBinUpEdge(binOrig) == hNew->GetXaxis()->GetBinUpEdge(bin) ) {
	  hNew->SetBinContent(bin,y/ye);
	  hNew->SetBinError(bin,1./sqrt(ye));

	  bin++;
	  y = 0.;
	  ye = 0.;
	}
      }


      return hNew;
    }
    

    
  private:
    static unsigned int COUNT_;
  };

  unsigned int HistOps::COUNT_ = 0;
}
#endif

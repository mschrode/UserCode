//  $Id: $

#ifndef LABEL_FACTORY_H
#define LABEL_FACTORY_H

#include <cmath>

#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"

namespace util {

  //!  Factory for TLegend and TPaveText objects
  //!
  //!  \author   Matthias Schroeder (www.desy.de/~matsch)
  //!  \date     2010/03/09
  //!  $Id: $
  // -------------------------------------------------------------------------------------
  class LabelFactory {
  public:
    static double lineHeight() {
      double height = 0.044;
      TString mode = "Presentation";
      if( mode.CompareTo(gStyle->GetTitle()) == 0 ) {
	height = 0.05;
      }
      return height;
    }

    static double labelTopOffset() { return 0.06; }
    static double labelSideOffset() { return 0.04; }
      

    static TLegend *createLegend(int nEntries, double width, double yOffset, double lineHgt) {
      double x0 = 0.;
      double y0 = 0.;
      double x1 = 0.;
      double y1 = 0.;
      cornerCoordinates(nEntries,width,1.1*lineHgt,x0,y0,x1,y1);
      TLegend *leg = new TLegend(x0,y0-yOffset,x1,y1-yOffset);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextFont(42);
      return leg;
    }
    static TLegend *createLegend(int nEntries) {
      return createLegend(nEntries,1.,0.,lineHeight());
    }
    static TLegend *createLegend(int nEntries, double lineHgt) {
      return createLegend(nEntries,1.,0.,lineHgt);
    }
    static TLegend *createLegendWithOffset(int nEntries, int yOffset) {
      return createLegend(nEntries,1.,1.16*lineHeight()*yOffset,lineHeight());
    }
    static TLegend *createLegendWithOffset(int nEntries, double yOffset) {
      return createLegend(nEntries,1.,yOffset,lineHeight());
    }
    static TLegend *createLegendWithOffset(int nEntries, int yOffset, double lineHgt) {
      return createLegend(nEntries,1.,1.16*lineHeight()*yOffset,lineHgt);
    }
    static TLegend *createLegendWithOffset(int nEntries, double yOffset, double lineHgt) {
      return createLegend(nEntries,1.,yOffset,lineHgt);
    }
    static TLegend *createLegendCol(int nEntries, double width) {
      return createLegend(nEntries,width,0.,lineHeight());
    }
    static TLegend *createLegendCol(int nEntries, double width, double lineHgt) {
      return createLegend(nEntries,width,0.,lineHgt);
    }
    static TLegend *createLegendColWithOffset(int nEntries, double width, int yOffset) {
      return createLegend(nEntries,width,1.16*lineHeight()*yOffset,lineHeight());
    }

    static TH1 *addExtraLegLine(TLegend *leg, const TString &entry) {
      ++nDummies_;
      TString name = "hDummy";
      name += nDummies_;
      TH1 *hDummy = new TH1D(name,"",1,0,1);
      hDummy->SetLineColor(0);
      hDummy->SetMarkerColor(0);
      hDummy->SetMarkerStyle(0);
      leg->AddEntry(hDummy,entry,"L");
      
      return hDummy;
    }


    static TPaveText *createPaveText(int nEntries, double width = 1., double lineHgt = -1) {
      double x0 = 0.;
      double y0 = 0.;
      double x1 = 0.;
      double y1 = 0.;
      cornerCoordinates(nEntries,width,lineHgt>0 ? lineHgt : lineHeight(),x0,y0,x1,y1);
      TPaveText *txt = new TPaveText(x0,y0,x1,y1,"NDC");
      txt->SetBorderSize(0);
      txt->SetFillColor(0);
      txt->SetTextFont(42);
      txt->SetTextAlign(12);
      return txt;
    }

  private:
    static int nDummies_;

    static void cornerCoordinates(int nEntries, double width, double lineHgt, double &x0, double &y0, double &x1, double &y1) {
      x0 = gStyle->GetPadLeftMargin()+labelSideOffset();
      x1 = 1.-(gStyle->GetPadRightMargin()+labelSideOffset());
      if( width > 0. ) x0 = x0 + (1.-width)*(x1-x0);
      else x1 = x1 - (1.+width)*(x1-x0);
      y1 = 1.-(gStyle->GetPadTopMargin()+labelTopOffset());
      double height = lineHgt;
      y0 = y1 - nEntries*height;
    }
  };

  int LabelFactory::nDummies_ = 0;
}
#endif

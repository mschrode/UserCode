#ifndef LABEL_FACTORY_H
#define LABEL_FACTORY_H

#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"

namespace util {
  class LabelFactory {
  public:
    static double lineHeight() { return 0.05; }
    static TLegend *createLegend(int nEntries, double width = 1., double lineHgt = -1, double yOffset = 0.);
    static TPaveText *createPaveText(int nEntries, double width = 1., double lineHgt = -1);
  };


  TLegend *LabelFactory::createLegend(int nEntries, double width, double lineHgt, double yOffset) {
    double margin = 0.04;
    double x0 = gStyle->GetPadLeftMargin()+margin;
    double x1 = 1.-(gStyle->GetPadRightMargin()+margin);
    x0 = x0 + (1.-width)*(x1-x0);
    double y1 = 1.-(gStyle->GetPadTopMargin()+margin+0.02+yOffset);
    double height = lineHeight();
    if( lineHgt > 0 ) height = lineHgt;
    double y0 = y1 - nEntries*height;
    TLegend *leg = new TLegend(x0,y0,x1,y1);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    return leg;
  }


  TPaveText *LabelFactory::createPaveText(int nEntries, double width, double lineHgt) {
    double margin = 0.04;
    double x0 = gStyle->GetPadLeftMargin()+margin;
    double x1 = 1.-(gStyle->GetPadRightMargin()+margin);
    x0 = x0 + (1.-width)*(x1-x0);
    double y1 = 1.-(gStyle->GetPadTopMargin()+margin+0.02);
    double height = lineHeight();
    if( lineHgt > 0 ) height = lineHgt;
    double y0 = y1 - nEntries*height;
    TPaveText *txt = new TPaveText(x0,y0,x1,y1,"NDC");
    txt->SetBorderSize(0);
    txt->SetFillColor(0);
    txt->SetTextFont(42);
    txt->SetTextAlign(12);
    return txt;
  }
}
#endif

//  $Id: LabelFactory.h,v 1.21 2011/08/26 10:21:56 mschrode Exp $

#ifndef LABEL_FACTORY_H
#define LABEL_FACTORY_H

#include <cmath>
#include <vector>

#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"

#include "utils.h"


namespace util {

  //!  Factory for TLegend and TPaveText objects
  //!
  //!  \author   Matthias Schroeder (www.desy.de/~matsch)
  //!  \date     2010/03/09
  //!  $Id: LabelFactory.h,v 1.21 2011/08/26 10:21:56 mschrode Exp $
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

    static double labelTopOffset() { 
      double offset = 0.03;
      TString mode = "Presentation";
      if( mode.CompareTo(gStyle->GetTitle()) == 0 ) {
	offset = 0.04;
      }
      return offset;
    }

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
    static TLegend *createLegendColWithOffset(int nEntries, double width, double yOffset) {
      return createLegend(nEntries,width,yOffset,lineHeight());
    }

    static void addExtraLegLine(TLegend *leg, const TString &entry) {
      TString name = "hDummy";
      name += hDummies_.size();
      TH1 *hDummy = new TH1D(name,"",1,0,1);
      hDummy->SetLineColor(0);
      hDummy->SetMarkerColor(0);
      hDummy->SetMarkerStyle(0);
      leg->AddEntry(hDummy,entry,"L");
      hDummies_.push_back(hDummy);
    }

    static void deleteDummies() {
      for(std::vector<TH1*>::iterator it = hDummies_.begin();
	  it != hDummies_.end(); ++it) {
	delete *it;
      }
      hDummies_.clear();
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


    static TPaveText *createPaveTextWithOffset(int nEntries, double width, int yOffset) {
      double x0 = 0.;
      double y0 = 0.;
      double x1 = 0.;
      double y1 = 0.;
      cornerCoordinates(nEntries,width,lineHeight(),x0,y0,x1,y1);
      TPaveText *txt = new TPaveText(x0,y0-1.16*lineHeight()*yOffset,x1,y1-1.16*lineHeight()*yOffset,"NDC");
      txt->SetBorderSize(0);
      txt->SetFillColor(0);
      txt->SetTextFont(42);
      txt->SetTextAlign(12);
      return txt;
    }


    static TPaveText *createPaveTextWithOffset(int nEntries, double width, double yOffset) {
      double x0 = 0.;
      double y0 = 0.;
      double x1 = 0.;
      double y1 = 0.;
      cornerCoordinates(nEntries,width,lineHeight(),x0,y0,x1,y1);
      TPaveText *txt = new TPaveText(x0,y0-yOffset,x1,y1-yOffset,"NDC");
      txt->SetBorderSize(0);
      txt->SetFillColor(0);
      txt->SetTextFont(42);
      txt->SetTextAlign(12);
      return txt;
    }


    // -------------------------------------------------------------------------------------
    static TString jetAlgo(const TString &name) {
      TString algo = "DefaultJets";
      if( name.Contains("AK5") || name.Contains("ak5") ) {
	algo = "AK5";
      } else {
	std::cerr << "WARNING in LabelFactory::jetAlgo(): unknown jet algorithm in name '" << name << "'" << std::endl;
      }

      return algo;
    }


    // -------------------------------------------------------------------------------------
    static TString labelJetAlgo(const TString &name) {
      TString algo = "default jets";
      if( jetAlgo(name) == "AK5" ) {
	algo = "Anti k_{T} (R=0.5)";
      } else {
	std::cerr << "WARNING in LabelFactory::jetAlgo(): unknown jet algorithm in name '" << name << "'" << std::endl;
      }

      return algo;
    }


    // -------------------------------------------------------------------------------------
    static TString jetType(const TString &name) {
      TString algo = "DefaultJets";
      if( name.Contains("Calo") || name.Contains("calo") ) {
	algo = "Calo";
      } else if( name.Contains("PF") || name.Contains("pf") ) {
	algo = "PF";
      } else {
	std::cerr << "WARNING in LabelFactory::jetType(): unknown jet type in name '" << name << "'" << std::endl;
      }

      return algo;
    }


    // -------------------------------------------------------------------------------------
    static TString labelJetType(const TString &name) {
      TString algo = "default jets";
      if( jetType(name) == "Calo" ) {
	algo = "Calo-Jets";
      } else if( jetType(name) == "PF" ) {
	algo = "PF-Jets";
      } else {
	std::cerr << "WARNING in LabelFactory::jetType(): unknown jet type in name '" << name << "'" << std::endl;
      }

      return algo;
    }


    // -------------------------------------------------------------------------------------
    static TString labelJet(const TString &name) {
      return labelJetAlgo(name)+" "+labelJetType(name);
    }


    // -------------------------------------------------------------------------------------
    static TString labelJetAlgo(const TString &name1, const TString &name2) {
      TString algo = labelJetAlgo(name1);
      if( algo != labelJetAlgo(name2) ) {
	std::cerr << "WARNING in LabelFactory::jetAlgo(): inconsistent jet algorithms in files '" << name1 << "' and '" << name2 << "'" << std::endl;
	algo = "default jets";
      }
      return algo;
    }


    // -------------------------------------------------------------------------------------
    static TString labelEta(double etaMin, double etaMax, const TString &etaExp = "") {
      TString min = util::toTString(etaMin,1);
      TString max = util::toTString(etaMax,1);
      if( min.Length() == 1 ) min += ".0";
      if( max.Length() == 1 ) max += ".0";

      return etaExp == "" ? min+" < |#eta| < "+max :  min+" < |#eta^{"+etaExp+"}| < "+max;
    }


    // -------------------------------------------------------------------------------------
    static TString labelEtaGen(double etaMin, double etaMax) {
      return labelEta(etaMin,etaMax,"gen");
    }


    // -------------------------------------------------------------------------------------
    static TString labelPt3(double pt3) {
      return "p_{T,3} < "+util::toTString(pt3)+" #upoint p^{ave}_{T}";
    }
    
    
    
  private:
    static std::vector<TH1*> hDummies_;

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
}


#ifdef UTILS_AS_HEADER_FILE
std::vector<TH1*> util::LabelFactory::hDummies_;
#endif

#endif

// $Id: StyleSettings.h,v 1.13 2011/07/18 08:58:40 mschrode Exp $

#ifndef STYLE_SETTINGS_H
#define STYLE_SETTINGS_H

#include <cmath>
#include <iostream>

#include "TString.h"
#include "TStyle.h"

#include "../util/utils.h"


namespace util {

  //!  Encapsulates different pad and histogram styles
  //!
  //!  \author   Matthias Schroeder (www.desy.de/~matsch)
  //!  \date     2010/03/09
  //!  $Id: StyleSettings.h,v 1.13 2011/07/18 08:58:40 mschrode Exp $
  // -------------------------------------------------------------------------------------
  class StyleSettings {
  public:
    // Different style settings
    enum Style { Screen, Presentation, Note, PAS };

    // Convert between style representations
    static Style style() {
      Style st = Screen;
      TString mode = gStyle->GetTitle();
      if( mode == "Presentation" ) st = Presentation;
      else if( mode == "Note" ) st = Note;
      else if( mode == "PAS" ) st = PAS;
      
      return st;
    }
    static TString toString(Style mode) {
      TString st = "Screen";
      if( mode == Presentation ) st = "Presentation";
      else if( mode == Note ) st = "Note";
      else if( mode == PAS ) st = "PAS";
      
      return st;
    }

    // Set different styles
    static void setStyleScreen() { setStyle(Screen,true); }
    static void setStyleScreenNoTitle() { setStyle(Screen,false); }
    static void setStylePresentation() { setStyle(Presentation,true); }
    static void setStylePresentationNoTitle() { setStyle(Presentation,false); }
    static void setStyleNote() { setStyle(Note,true); }
    static void setStyleNoteNoTitle() { setStyle(Note,false); }
    static void setStylePAS() { setStyle(PAS,true); }

    // Retrieve some style attributes to be used in
    // user code. Their value depends on the current
    // global style settings:

    // Get a histogram title string for CMS PAS style
    static TString title(double lumi, bool isPreliminary = false) {
      TString title = "CMS";
      if( isPreliminary ) title += " preliminary";
      title += ", L = "+luminosity(lumi)+",  #sqrt{s} = 7 TeV";
      
      return title;
    }

    // This needs the lumi in 1/pb and returns
    // a nicely formatted TString
    static TString luminosity(double lumi) {
      TString lab = "";
      if( lumi > 800. ) {
	lumi /= 1000.;
	lab = toTString(lumi,2)+" fb^{-1}";
      } else {
	lab = toTString(lumi,3)+" pb^{-1}";
      }

      return lab;
    }

    // Return a readable color; useful for loops
    static int color(int i) {
      int col[5] = { 1, 2, 4, 7, 8 };
      int idx = i%5;
      return (idx>=0 && idx<5) ? col[idx] : 1;
    }

    // Line width for histograms
    static int lineWidth() {
      int width = 1;
      if( style() == Presentation ) {
	width = 3;
      } else if( style() == Note || style() == PAS ) {
	width = 2;
      }

      return width;
    }

    
  private:
    // Helper function: set the actual gStyle attributes
    static void setStyle(Style mode, bool spaceForTitle) {
      // Set title of current style object
      gStyle->SetTitle(toString(mode));

      // Zero horizontal error bars
      gStyle->SetErrorX(0);

      //  For 'colz' TH2
      gStyle->SetPalette(1);
    
      //  For the canvas
      gStyle->SetCanvasBorderMode(0);
      gStyle->SetCanvasColor(kWhite);
      gStyle->SetCanvasDefH(800); //Height of canvas
      gStyle->SetCanvasDefW(800); //Width of canvas
      gStyle->SetCanvasDefX(0);   //Position on screen
      gStyle->SetCanvasDefY(0);
    
      //  For the frame
      gStyle->SetFrameBorderMode(0);
      gStyle->SetFrameBorderSize(1);
      gStyle->SetFrameFillColor(kBlack);
      gStyle->SetFrameFillStyle(0);
      gStyle->SetFrameLineColor(kBlack);
      gStyle->SetFrameLineStyle(0);
      gStyle->SetFrameLineWidth(1);
    
      //  For the Pad
      gStyle->SetPadBorderMode(0);
      gStyle->SetPadColor(kWhite);
      gStyle->SetPadGridX(false);
      gStyle->SetPadGridY(false);
      gStyle->SetGridColor(0);
      gStyle->SetGridStyle(3);
      gStyle->SetGridWidth(1);

      //  Margins
      if( mode == Presentation ) {
	if( spaceForTitle ) {
	  gStyle->SetPadTopMargin(0.11);
	  gStyle->SetPadBottomMargin(0.18);
	  gStyle->SetPadLeftMargin(0.25);
	  gStyle->SetPadRightMargin(0.04);
	} else {
	  gStyle->SetPadTopMargin(0.05);
	  gStyle->SetPadBottomMargin(0.18);
	  gStyle->SetPadLeftMargin(0.19);
	  gStyle->SetPadRightMargin(0.04);
	}
      } else if( mode == Note || mode == PAS  ) {
	if( spaceForTitle ) {
	  gStyle->SetPadTopMargin(0.06);
	  gStyle->SetPadBottomMargin(0.18);
	  gStyle->SetPadLeftMargin(0.2);
	  gStyle->SetPadRightMargin(0.04);
	} else {
	  gStyle->SetPadTopMargin(0.05);
	  gStyle->SetPadBottomMargin(0.17);
	  gStyle->SetPadLeftMargin(0.18);
	  gStyle->SetPadRightMargin(0.04);
	}
      } else {
	if( spaceForTitle ) {
	  gStyle->SetPadTopMargin(0.10);
	  gStyle->SetPadBottomMargin(0.14);
	  gStyle->SetPadLeftMargin(0.18);
	  gStyle->SetPadRightMargin(0.04);
	} else {
	  gStyle->SetPadTopMargin(0.08);
	  gStyle->SetPadBottomMargin(0.14);
	  gStyle->SetPadLeftMargin(0.18);
	  gStyle->SetPadRightMargin(0.04);
	}
      }

      //  For the histo:
      gStyle->SetHistLineColor(kBlack);
      gStyle->SetHistLineStyle(0);
      gStyle->SetHistLineWidth(1);
    
      //  For the statistics box:
      if( mode == Screen ) {
	gStyle->SetOptStat("eMR");
	gStyle->SetStatColor(kWhite);
	gStyle->SetStatFont(42);
	gStyle->SetStatFontSize(0.03);
	gStyle->SetStatTextColor(1);
	gStyle->SetStatFormat("6.4g");
	gStyle->SetStatBorderSize(1);
	gStyle->SetStatX(0.94);              
	gStyle->SetStatY(0.86);              
	gStyle->SetStatH(0.16);
	gStyle->SetStatW(0.22);
      } else {
	gStyle->SetOptStat(0);
      }
    
      //  For the Global title:
      gStyle->SetOptTitle(1);
      gStyle->SetTitleFont(42,"");
      gStyle->SetTitleColor(1);
      gStyle->SetTitleTextColor(1);
      gStyle->SetTitleFillColor(0);
      gStyle->SetTitleFontSize(0.1);
      gStyle->SetTitleAlign(13);
      gStyle->SetTitleX(0.6);
      gStyle->SetTitleH(0.05);
      gStyle->SetTitleBorderSize(0);
      if( mode == PAS  ) {
	if( spaceForTitle ) {
	  gStyle->SetTitleAlign(13);
	  gStyle->SetTitleX(0.19);
	  gStyle->SetTitleH(0.038);
	}
      }	  

      //  For the axis
      gStyle->SetAxisColor(1,"XYZ");
      gStyle->SetTickLength(0.03,"XYZ");
      gStyle->SetNdivisions(510,"XYZ");
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);
      gStyle->SetStripDecimals(kFALSE);
    
      //  For the axis labels and titles
      gStyle->SetTitleColor(1,"XYZ");
      gStyle->SetLabelColor(1,"XYZ");
      if( mode == Presentation ) {
	// For the axis labels:
	gStyle->SetLabelFont(42,"XYZ");
	gStyle->SetLabelOffset(0.007,"XYZ");
	gStyle->SetLabelSize(0.045,"XYZ");
      
	// For the axis titles:
	gStyle->SetTitleFont(42,"XYZ");
	gStyle->SetTitleSize(0.06,"XYZ");
	gStyle->SetTitleXOffset(1.2);
	if( spaceForTitle ) gStyle->SetTitleYOffset(2.0);
	else                gStyle->SetTitleYOffset(1.5);
      } else if ( mode == Note || mode == PAS ) {
	// For the axis labels:
	gStyle->SetLabelFont(42,"XYZ");
	gStyle->SetLabelOffset(0.007,"XYZ");
	gStyle->SetLabelSize(0.04,"XYZ");
      
	// For the axis titles:
	gStyle->SetTitleFont(42,"XYZ");
	gStyle->SetTitleSize(0.045,"XYZ");
	gStyle->SetTitleXOffset(1.5);
	if( spaceForTitle ) gStyle->SetTitleYOffset(2.1);
	else                gStyle->SetTitleYOffset(1.8);
      } else {
	// For the axis labels:
	gStyle->SetLabelFont(42,"XYZ");
	gStyle->SetLabelOffset(0.007,"XYZ");
	gStyle->SetLabelSize(0.035,"XYZ");
      
	// For the axis titles:
	gStyle->SetTitleFont(42,"XYZ");
	gStyle->SetTitleSize(0.04,"XYZ");
	gStyle->SetTitleXOffset(1.5);
	if( spaceForTitle ) gStyle->SetTitleYOffset(2.1);
	else                gStyle->SetTitleYOffset(1.8);
      }


      //  For the legend
      gStyle->SetLegendBorderSize(0);

      //  For the statistics box
      if( mode == Presentation ) {
	if( spaceForTitle ) {
	  gStyle->SetStatFontSize(0.04);
	  gStyle->SetStatX(0.92);              
	  gStyle->SetStatY(0.86);              
	  gStyle->SetStatH(0.2);
	  gStyle->SetStatW(0.3);
	} else {
	  gStyle->SetStatFontSize(0.04);
	  gStyle->SetStatX(0.92);              
	  gStyle->SetStatY(0.92);              
	  gStyle->SetStatH(0.2);
	  gStyle->SetStatW(0.3);
	}
      } else {
	if( spaceForTitle ) {
	  gStyle->SetStatFontSize(0.03);
	  gStyle->SetStatX(0.92);              
	  gStyle->SetStatY(0.86);              
	  gStyle->SetStatH(0.16);
	  gStyle->SetStatW(0.22);
	} else {
	  gStyle->SetStatFontSize(0.03);
	  gStyle->SetStatX(0.92);              
	  gStyle->SetStatY(0.92);              
	  gStyle->SetStatH(0.16);
	  gStyle->SetStatW(0.22);
	}
      }

      std::cout << "Adjusted gStyle for " << std::flush;
      if( mode == Screen ) std::cout << "screen viewing" << std::flush;
      else if( mode == Note ) std::cout << "CMS Analysis Notes" << std::flush;
      else if( mode == PAS ) std::cout << "CMS PAS" << std::flush;
      else std::cout << "presentations" << std::flush;
      std::cout << " and " << std::flush;
      if( spaceForTitle ) std::cout << "histograms with title." << std::endl;
      else std::cout << "histograms without title." << std::endl;
    }
  };  
}
#endif

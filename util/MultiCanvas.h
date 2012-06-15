// $Id: $

#ifndef MULTI_CANVAS_H
#define MULTI_CANVAS_H

#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TH1.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

namespace util {

  class MultiCanvas {
  public:
    MultiCanvas(const TString& name, unsigned int nRows, unsigned int nCols, bool hasBottomRatio);
    ~MultiCanvas();

    unsigned int nPads() const { return nPads_; }
    unsigned int col(unsigned int p) const { return p%nCols_; }
    unsigned int row(unsigned int p) const { return p/nCols_; }

    void cd() { can_->cd(); };
    void setLogx();
    void setLogy();
    TCanvas* canvas() { return can_; }
    TPad* mainPad(unsigned int p) { return padMain_.at(p); }

    TH1* mainFrame(const TH1* h, unsigned int p);
    void noExponentX();
    void noExponentY();
    void moreLogLabelsX();
    void moreLogLabelsY();
    void deleteFrames();
    void clearFrames() {
      framesMain_.clear();
      if( hasRatio_ ) framesRatio_.clear();
    }

  private:
    const bool hasRatio_;
    const unsigned int nCols_;
    const unsigned int nRows_;
    const unsigned int nPads_;

    TCanvas* can_;
    std::vector<TPad*> padMain_;
    std::vector<TPad*> padRatio_;

    std::vector<TH1*> framesMain_;
    std::vector<TH1*> framesRatio_;
  };


  MultiCanvas::MultiCanvas(const TString& name, unsigned int nRows, unsigned int nCols, bool hasBottomRatio)
    : hasRatio_(hasBottomRatio), nCols_(nCols), nRows_(nRows), nPads_(nCols*nRows) {

    // Dimensions
    const double relMarginTop = gStyle->GetPadTopMargin()/(1.-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin());
    const double relMarginBot = gStyle->GetPadBottomMargin()/(1.-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin());
    const double relMarginLeft = gStyle->GetPadLeftMargin()/(1.-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin());
    const double relMarginRight = gStyle->GetPadRightMargin()/(1.-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin());

//     const double relMarginTop = 0.02;
//     const double relMarginBot = 0.3;
//     const double relMarginLeft = 0.3;
//     const double relMarginRight = 0.02;


    const double padHeight = 1./(nRows_+relMarginTop+relMarginBot);
    const double marginTop = relMarginTop*padHeight;
    const double marginBot = relMarginBot*padHeight;
    const double padWidth = 1./(nCols_+relMarginTop+relMarginBot);
    const double marginLeft = relMarginLeft*padWidth;
    const double marginRight = relMarginRight*padWidth;

//     std::cout << "  padHeight = " << padHeight << std::endl;
//     std::cout << "   padWidth = " << padWidth << std::endl;
//     std::cout << "  marginBot = " << marginBot << std::endl;
//     std::cout << "  marginTop = " << marginTop << std::endl;
//     std::cout << " marginLeft = " << marginLeft << std::endl;
//     std::cout << "marginRight = " << marginRight << std::endl;

    // Main canvas
    const int canWidth = 700;
    const int canHeight = static_cast<int>(padWidth/padHeight*canWidth);
    can_ = new TCanvas(name,name,canWidth,canHeight);
    can_->SetWindowSize(canWidth + (canWidth - can_->GetWw()), canHeight + (canHeight - can_->GetWh())); // No clue why I need this line: ROOT sucks!

    //    std::cout << "\n can: " << canWidth << " x " << canHeight << std::endl;

    // Pads
    for(unsigned int r = 0; r < nRows_; ++r) {
      for(unsigned int c = 0; c < nCols_; ++c) {
	unsigned int p = nCols_*r+c;

// 	double left = 0.;
// 	if( c > 0 ) left += marginLeft+c*padWidth;
// 	double right = marginLeft+(1+c)*padWidth;
// 	if( c == nCols_-1 ) right = 1.;
// 	double top = 1.;
// 	if( r > 0 ) top -= marginTop+r*padHeight;
// 	double bottom = 1. - (marginTop+(1.+r)*padHeight);
// 	if( r == nRows_-1 ) bottom = 0.;

	double left = c*padWidth;
	double right = marginLeft+(1+c)*padWidth+marginRight;
	double top = 1. - r*padHeight;
	double bottom = 1. - (marginTop+(1+r)*padHeight+marginBot);
	
// 	std::cout << "\n\nCreating pad " << p << ":" << std::endl;
// 	std::cout << "   left  = " << left << std::endl;
// 	std::cout << "   right = " << right << std::endl;
// 	std::cout << "  bottom = " << bottom << std::endl;
// 	std::cout << "     top = " << top << std::endl;

	TString id = "util::MultiCanvas::MainPad";
	id += p;
	TPad* pad = new TPad(id,id,left,bottom,right,top);
	pad->SetFillStyle(0);
	pad->SetFrameFillColor(0);
	pad->SetFrameBorderMode(0);
//  	if( r == 0 && nRows_ == 1 ) {
//  	  pad->SetTopMargin(marginTop/(marginTop+padHeight+marginBot));
//  	  pad->SetBottomMargin(marginBot/(marginTop+padHeight+marginBot));
//  	} else if( r == 0 ) {
//  	  pad->SetTopMargin(marginTop/(marginTop+padHeight));
//  	  pad->SetBottomMargin(0.);
//  	} else if( r == nRows_-1 ) {
//  	  pad->SetTopMargin(0.);
//  	  pad->SetBottomMargin(marginBot/(padHeight+marginBot));
//  	} else {
//  	  pad->SetTopMargin(0.);
//  	  pad->SetBottomMargin(0.);
//  	}
//  	if( c == 0 && nCols_ == 1 ) {
//  	  pad->SetLeftMargin(marginLeft/(marginLeft+padWidth+marginRight));
//  	  pad->SetRightMargin(marginRight/(marginLeft+padWidth+marginRight));
//  	} else if( c == 0 ) {
//  	  pad->SetLeftMargin(marginLeft/(marginLeft+padWidth));
//  	  pad->SetRightMargin(0.);
//  	} else if( c == nCols_-1 ) {
//  	  pad->SetLeftMargin(0.);
//  	  pad->SetRightMargin(marginRight/(padWidth+marginRight));
//  	} else {
//  	  pad->SetLeftMargin(0.);
//  	  pad->SetRightMargin(0.);
//  	}

	pad->SetTopMargin(marginTop/(marginTop+padHeight+marginBot));
	pad->SetBottomMargin(marginBot/(marginTop+padHeight+marginBot));
	pad->SetLeftMargin(marginLeft/(marginLeft+padWidth+marginRight));
	pad->SetRightMargin(marginRight/(marginLeft+padWidth+marginRight));

	padMain_.push_back(pad);
	if( hasRatio_ ) {
	  std::cerr << "ERROR: ratio pad not yet implemented." << std::endl;
	} else {
	  id = "util::MultiCanvas::RatioPadDummy";
	  id += p;
	  padRatio_.push_back(static_cast<TPad*>(pad->Clone(id)));
	}
      }	// End of loop over columns
    } // End of loop over rows
  }

  MultiCanvas::~MultiCanvas() {
    for(unsigned int i = 0; i < nPads_; ++i) {
      delete padMain_.at(i);
      delete padRatio_.at(i);
    }
    delete can_;
  }


  TH1* MultiCanvas::mainFrame(const TH1* h, unsigned int p) {
    if( p > nPads_ ) {
      std::cerr << "ERROR in MultiCanvas::mainFrame(): parameter value p = " << p << " > nPads()" << std::endl;
      exit(1);
    }
    TString name = h->GetName();
    name += "_MainFrame";
    name += p;
    TH1* hFrame = static_cast<TH1*>(h->Clone(name));
    if( col(p) > 0 ) {
      hFrame->GetYaxis()->SetTitle("");
      hFrame->GetYaxis()->SetLabelSize(0);
    }
    if( row(p) != nRows_-1 ) {
      hFrame->GetXaxis()->SetTitle("");
      hFrame->GetXaxis()->SetLabelSize(0);
    }

    framesMain_.push_back(hFrame);

    return hFrame;
  }

  void MultiCanvas::deleteFrames() {
    for(unsigned int i = 0; i < framesMain_.size(); ++i) {
      delete framesMain_.at(i);
      if( hasRatio_ ) delete framesRatio_.at(i);
    }
  }

  void MultiCanvas::moreLogLabelsX() {
    for(unsigned int i = 0; i < framesMain_.size(); ++i) {
      framesMain_.at(i)->GetXaxis()->SetMoreLogLabels();
      if( hasRatio_ ) framesRatio_.at(i)->GetXaxis()->SetMoreLogLabels();
    }
  }

  void MultiCanvas::moreLogLabelsY() {
    for(unsigned int i = 0; i < framesMain_.size(); ++i) {
      framesMain_.at(i)->GetYaxis()->SetMoreLogLabels();
    }
  }

  void MultiCanvas::noExponentX() {
    for(unsigned int i = 0; i < framesMain_.size(); ++i) {
      framesMain_.at(i)->GetXaxis()->SetNoExponent();
      if( hasRatio_ ) framesRatio_.at(i)->GetXaxis()->SetNoExponent();
    }
  }

  void MultiCanvas::noExponentY() {
    for(unsigned int i = 0; i < framesMain_.size(); ++i) {
      framesMain_.at(i)->GetYaxis()->SetNoExponent();
    }
  }

  void MultiCanvas::setLogx() {
    for(unsigned int i = 0; i < nPads_; ++i) {
      padMain_.at(i)->SetLogx();
      padRatio_.at(i)->SetLogx();
    }
  }

  void MultiCanvas::setLogy() {
    for(unsigned int i = 0; i < nPads_; ++i) {
      padMain_.at(i)->SetLogy();
      padRatio_.at(i)->SetLogy();
    }
  }
}
#endif

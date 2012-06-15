// $Id: MultiCanvas.h,v 1.1 2012/06/15 00:14:40 mschrode Exp $

#ifndef MULTI_CANVAS_H
#define MULTI_CANVAS_H

#include <cassert>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TObject.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"

namespace util {

  class MultiCanvas {
  public:
    MultiCanvas(const TString& name, unsigned int nRows, unsigned int nCols, bool hasBottomRatio);
    MultiCanvas(const TString& name, unsigned int nRows, unsigned int nCols, unsigned int nUsedPads, bool hasBottomRatio);
    ~MultiCanvas();

    unsigned int nPads() const { return nPads_; }
    unsigned int col(unsigned int p) const { return p%nCols_; }
    unsigned int row(unsigned int p) const { return p/nCols_; }

    void cd() { can_->cd(); };
    void setLogx();
    void setLogy();
    void noExponentX();
    void noExponentY();
    void moreLogLabelsX();
    void moreLogLabelsY();

    TCanvas* canvas() { return can_; }
    TPad* mainPad(unsigned int p) { return padMain_.at(p); }
    TPad* ratioPad(unsigned int p) { return padRatio_.at(p); }
    TH1* mainFrame(const TH1* h, unsigned int p);
    TH1* ratioFrame(const TH1 *h, const TString &yTitle, double yMin, double yMax, unsigned int p);

    void setName(const TString &name) { can_->SetName(name); }

    void adjustLegend(TLegend* leg) const;
    void adjustPaveText(TPaveText* txt) const;
    void markForDeletion(TObject* obj) const { objForDeletion_.push_back(obj); }

    void reset();
    void deleteObjects() const;
    void deleteFrames();
    void clearFrames() {
      framesMain_.clear();
      if( hasRatio_ ) framesRatio_.clear();
    }


  private:
    const double topPadFrac_;
    const bool hasRatio_;
    const unsigned int nCols_;
    const unsigned int nRows_;
    const unsigned int nPads_;
    const unsigned int nUsedPads_;

    TCanvas* can_;
    std::vector<TPad*> padMain_;
    std::vector<TPad*> padRatio_;

    std::vector<TH1*> framesMain_;
    std::vector<TH1*> framesRatio_;
    mutable std::vector<TObject*> objForDeletion_; // Garbage collector

    void init(const TString &name);
  };


  MultiCanvas::MultiCanvas(const TString& name, unsigned int nRows, unsigned int nCols, bool hasBottomRatio)
    : topPadFrac_(0.8), hasRatio_(hasBottomRatio), nCols_(nCols), nRows_(nRows), nPads_(nCols*nRows), nUsedPads_(nPads_) {
    init(name);
  }


  MultiCanvas::MultiCanvas(const TString& name, unsigned int nRows, unsigned int nCols, unsigned int nUsedPads, bool hasBottomRatio)
    : topPadFrac_(0.8), hasRatio_(hasBottomRatio), nCols_(nCols), nRows_(nRows), nPads_(nCols*nRows), nUsedPads_(nUsedPads) {
    assert( nUsedPads_ <= nPads_ );
    init(name);
  }


  void MultiCanvas::init(const TString &name) {
    // Dimensions
    const double gt = gStyle->GetPadTopMargin();
    const double gb = gStyle->GetPadBottomMargin();
    const double gl = gStyle->GetPadLeftMargin();
    const double gr = gStyle->GetPadRightMargin();
    const double relMarginTop = gt/(1.-gb-gt);
    const double relMarginBot = gb/(1.-gb-gt);
    const double relMarginLeft = gl/(1.-gl-gr);
    const double relMarginRight = gr/(1.-gl-gr);

    const double padHeight = 1./(nRows_+relMarginTop+relMarginBot);
    const double marginTop = relMarginTop*padHeight;
    const double marginBot = relMarginBot*padHeight;
    const double padWidth = 1./(nCols_+relMarginTop+relMarginBot);
    const double marginLeft = relMarginLeft*padWidth;
    const double marginRight = relMarginRight*padWidth;

//      std::cout << "  padHeight = " << padHeight << std::endl;
//      std::cout << "   padWidth = " << padWidth << std::endl;
//      std::cout << "  marginBot = " << marginBot << std::endl;
//      std::cout << "  marginTop = " << marginTop << std::endl;
//      std::cout << " marginLeft = " << marginLeft << std::endl;
//      std::cout << "marginRight = " << marginRight << std::endl;

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

	double left = c*padWidth;
	double right = marginLeft+(1+c)*padWidth+marginRight;
	double top = 1. - r*padHeight;
	double bottom = 1. - (marginTop+(1+r)*padHeight+marginBot);
	
//  	std::cout << "\n\nCreating pad " << p << ":" << std::endl;
//  	std::cout << "   left  = " << left << std::endl;
//  	std::cout << "   right = " << right << std::endl;
//  	std::cout << "  bottom = " << bottom << std::endl;
//  	std::cout << "     top = " << top << std::endl;

	can_->cd();
	TString id = name+"::MainPad";
	id += p;
	TPad* pad = new TPad(id,id,left,bottom,right,top);
	pad->SetFillStyle(0);
	pad->SetFrameFillColor(0);
	pad->SetFrameBorderMode(0);
	pad->SetTopMargin(marginTop/(marginTop+padHeight+marginBot));
	pad->SetBottomMargin(marginBot/(marginTop+padHeight+marginBot));
	pad->SetLeftMargin(marginLeft/(marginLeft+padWidth+marginRight));
	pad->SetRightMargin(marginRight/(marginLeft+padWidth+marginRight));
	padMain_.push_back(pad);

	if( hasRatio_ ) {
	  pad->SetBottomMargin((1.-topPadFrac_) + topPadFrac_*pad->GetBottomMargin() - (1.-topPadFrac_)*pad->GetTopMargin());

	  can_->cd();
	  id = name+"::RatioPad";
	  id += p;
	  TPad* ratio = new TPad(id,id,0.,0.,1.,1.);
	  ratio->SetFillStyle(0);
	  ratio->SetFrameFillColor(0);
	  ratio->SetFrameBorderMode(0);
	  ratio->SetTopMargin(marginTop/(marginTop+padHeight+marginBot));
	  ratio->SetBottomMargin(marginBot/(marginTop+padHeight+marginBot));
	  ratio->SetLeftMargin(marginLeft/(marginLeft+padWidth+marginRight));
	  ratio->SetRightMargin(marginRight/(marginLeft+padWidth+marginRight));
	  ratio->SetTopMargin(topPadFrac_ - topPadFrac_*ratio->GetBottomMargin() + (1.-topPadFrac_)*ratio->GetTopMargin());
	  padRatio_.push_back(ratio);
	} else {
	  id = "util::MultiCanvas::RatioPadDummy";
	  id += p;
	  padRatio_.push_back(static_cast<TPad*>(pad->Clone(id)));
	}
      }	// End of loop over columns
    } // End of loop over rows
  }

  MultiCanvas::~MultiCanvas() {
//     for(unsigned int i = 0; i < nPads_; ++i) {
//       delete padMain_.at(i);
//     }
    delete can_;
    deleteFrames();
    deleteObjects();
  }

  void MultiCanvas::reset() {
    TString name = can_->GetName();
    for(unsigned int i = 0; i < nPads_; ++i) {
      delete padMain_.at(i);
    }
    delete can_;
    deleteObjects();
    deleteFrames();
    init(name);
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
    if( hasRatio_ || (row(p) != nRows_-1 && p+nCols_ < nUsedPads_) ) {
      hFrame->GetXaxis()->SetTitle("");
      hFrame->GetXaxis()->SetLabelSize(0);
    }
    if( hasRatio_ ) {
      hFrame->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/topPadFrac_);
      //hFrame->GetYaxis()->CenterTitle();
    }

    framesMain_.push_back(hFrame);

    return hFrame;
  }

  TH1* MultiCanvas::ratioFrame(const TH1 *h, const TString &yTitle, double yMin, double yMax, unsigned int p) {
    if( p > nPads_ ) {
      std::cerr << "ERROR in MultiCanvas::ratioFrame(): parameter value p = " << p << " > nPads()" << std::endl;
      exit(1);
    }
    TString name = h->GetName();
    name += "_RatioFrame";
    name += p;
    TH1* hFrame = static_cast<TH1*>(h->Clone(name));
    for(int bin = 1; bin <= hFrame->GetNbinsX(); ++bin) {
      hFrame->SetBinContent(bin,1.);
    }
    hFrame->SetLineStyle(2);
    hFrame->SetLineWidth(2);
    hFrame->GetYaxis()->SetTitle(yTitle);
    hFrame->GetYaxis()->CenterTitle();
    hFrame->GetYaxis()->SetTitleOffset(0.9*hFrame->GetYaxis()->GetTitleOffset());
    hFrame->GetYaxis()->SetNdivisions(205);
    hFrame->GetYaxis()->SetRangeUser(yMin,yMax);
    hFrame->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/(1.-topPadFrac_));
    hFrame->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X"));
    if( col(p) > 0 ) {
      hFrame->GetYaxis()->SetTitle("");
      hFrame->GetYaxis()->SetLabelSize(0);
    }
    if( row(p) != nRows_-1 && p+nCols_ < nUsedPads_ ) {
      hFrame->GetXaxis()->SetTitle("");
      hFrame->GetXaxis()->SetLabelSize(0);
    }

    framesRatio_.push_back(hFrame);

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

  void MultiCanvas::adjustLegend(TLegend* leg) const {
    leg->SetTextSize(1.2*leg->GetTextSize());
    leg->SetY1(leg->GetY1()-0.02);
    leg->SetY2(leg->GetY2()-0.02);
  }

  void MultiCanvas::adjustPaveText(TPaveText* txt) const {
    txt->SetTextSize(1.2*txt->GetTextSize());
    txt->SetY1(txt->GetY1()-0.02);
    txt->SetY2(txt->GetY2()-0.02);
  }

  void MultiCanvas::deleteObjects() const {
    for(std::vector<TObject*>::iterator it = objForDeletion_.begin();
	it != objForDeletion_.end(); ++it) {
      delete *it;
    }
  }
}
#endif

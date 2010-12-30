#ifndef FUNC
#define FUNC

#include "TH1.h"
#include "TF1.h"
#include "TString.h"

#include "../util/HistOps.h"


namespace func {
  double gaussInt(double mean, double sigma, double min, double max);
  double getTail(const TH1* hAsym, double nSigCore, double nSigTailStart, double nSigTailEnd, TH1* &hTail, TH1* &hTailClean, TF1* &gauss);
  double getTailCut(const TH1* hAsym, double cut, TH1* &hTail, TH1* &hTailClean);
  double getTailFromGauss(const TH1* hAsym, const TF1* extGauss, double nSigTailStart, double nSigTailEnd, double nSig, TH1* &hTail, TH1* &hTailClean, TF1* &gauss);
  void smearHistogram(const TH1* hOrig, TH1* &hSmeared, double nTotal, double width, double scaling);


  // --------------------------------------------------
  double gaussInt(double mean, double sigma, double min, double max) {
    return 0.5*(erf((max-mean)/sqrt(2.)/sigma) - erf((min-mean)/sqrt(2.)/sigma));
  }


  // --------------------------------------------------
  double getTail(const TH1* hAsym, double nSigCore, double nSigTailStart, double nSigTailEnd, TH1* &hTail, TH1* &hTailClean, TF1* &gauss) {
    double numTail = -1.;
  
    TString name = hAsym->GetName();
    name += "_Tails";
    hTail = static_cast<TH1D*>(hAsym->Clone(name));
    hTail->SetLineColor(kBlue);
    hTail->SetMarkerColor(kBlue);

    name = hAsym->GetName();
    name += "_TailsClean";
    hTailClean = static_cast<TH1D*>(hAsym->Clone(name));
    hTailClean->Reset();

    double width = 0.;
    double widthErr = 1000.;
    if( util::HistOps::fitCoreWidth(hAsym,nSigCore,gauss,width,widthErr) ) {
      double mean = gauss->GetParameter(1);
      int optLeftTailStartBin = hTail->FindBin(mean-1.*std::abs(nSigTailStart*width));
      int optLeftTailEndBin = hTail->FindBin(mean-1.*std::abs(nSigTailEnd*width));
      int optRightTailStartBin = hTail->FindBin(mean+std::abs(nSigTailStart*width));
      int optRightTailEndBin = hTail->FindBin(mean+std::abs(nSigTailEnd*width));
      for(int bin = 1; bin <= hTail->GetNbinsX(); ++bin) {
	double min = hTail->GetXaxis()->GetBinLowEdge(bin);
	double max = hTail->GetXaxis()->GetBinUpEdge(bin);
	double gaussPdf = gauss->Integral(min,max)/hTail->GetBinWidth(1);
	double tailPdf = hTail->GetBinContent(bin) - gaussPdf;
	if( tailPdf < 0. ) tailPdf = 0.;
	hTail->SetBinContent(bin,tailPdf);
	if( (bin >= optLeftTailEndBin && bin <= optLeftTailStartBin) ||
	    (bin >= optRightTailStartBin && bin <= optRightTailEndBin ) ) {
	  hTailClean->SetBinContent(bin,hTail->GetBinContent(bin));
	  hTailClean->SetBinError(bin,hTail->GetBinError(bin));
	}
      }
      numTail = hTailClean->Integral("width");  // Assumes normalized (area) distrbutions i.e. pdfs
    }

    return numTail;
  }


  // --------------------------------------------------
  double getTailFromGauss(const TH1* hAsym, const TF1* extGauss, double nSigTailStart, double nSigTailEnd, double nSig, TH1* &hTail, TH1* &hTailClean, TF1* &gauss) {
    double numTail = -1.;
  
    TString name = hAsym->GetName();
    name += "_Tails";
    hTail = static_cast<TH1D*>(hAsym->Clone(name));
    hTail->SetLineColor(kBlue);
    hTail->SetMarkerColor(kBlue);
  
    name = hAsym->GetName();
    name += "_TailsClean";
    hTailClean = static_cast<TH1D*>(hAsym->Clone(name));
    hTailClean->Reset();
  
    // Fit Gaussian
    double width = 0.;
    double widthErr = 1000.;
    util::HistOps::fitCoreWidth(hAsym,nSig,gauss,width,widthErr);

    // Extract tail using external Gaussian and tail start
    TF1 *extGaussTmp = static_cast<TF1*>(extGauss->Clone("func::getTailFromGauss::extGaussTmp"));
    double mean = extGaussTmp->GetParameter(1);
    width = extGaussTmp->GetParameter(2);
    int optLeftTailStartBin = hTail->FindBin(mean-1.*std::abs(nSigTailStart*width));
    int optLeftTailEndBin = hTail->FindBin(mean-1.*std::abs(nSigTailEnd*width));
    int optRightTailStartBin = hTail->FindBin(mean+std::abs(nSigTailStart*width));
    int optRightTailEndBin = hTail->FindBin(mean+std::abs(nSigTailEnd*width));
    for(int bin = 1; bin <= hTail->GetNbinsX(); ++bin) {
      double min = hTail->GetXaxis()->GetBinLowEdge(bin);
      double max = hTail->GetXaxis()->GetBinUpEdge(bin);
      double gaussPdf = extGaussTmp->Integral(min,max)/hTail->GetBinWidth(1);
      double tailPdf = hTail->GetBinContent(bin) - gaussPdf;
      if( tailPdf < 0. ) tailPdf = 0.;
      hTail->SetBinContent(bin,tailPdf);
      if( (bin >= optLeftTailEndBin && bin <= optLeftTailStartBin) ||
	  (bin >= optRightTailStartBin && bin <= optRightTailEndBin ) ) {
	hTailClean->SetBinContent(bin,hTail->GetBinContent(bin));
	hTailClean->SetBinError(bin,hTail->GetBinError(bin));
      }
      numTail = hTailClean->Integral("width");  // Assumes normalized (area) distrbutions i.e. pdfs
    }
    delete extGaussTmp;
  
    return numTail;
  }


  // --------------------------------------------------
  double getTailCut(const TH1* hAsym, double cut, TH1* &hTail, TH1* &hTailClean) {
    double numTail = -1.;
  
    TString name = hAsym->GetName();
    name += "_Tails";
    hTail = static_cast<TH1D*>(hAsym->Clone(name));
    hTail->SetLineColor(kBlue);
    hTail->SetMarkerColor(kBlue);
  
    name = hAsym->GetName();
    name += "_TailsClean";
    hTailClean = static_cast<TH1D*>(hAsym->Clone(name));
    hTailClean->Reset();
  
    // Core range = non-tail
    int binMin = hTail->FindBin(-1.*std::abs(cut));
    int binMax = hTail->FindBin(std::abs(cut));
    for(int bin = binMin+1; bin < binMax; ++bin) {
      hTail->SetBinContent(bin,0.);
      hTail->SetBinError(bin,0.);
    }
    for(int bin = 1; bin < hTail->GetNbinsX(); ++bin) {
      if( hTail->GetBinContent(bin) ) {
	hTailClean->SetBinContent(bin,hTail->GetBinContent(bin));
	hTailClean->SetBinError(bin,hTail->GetBinError(bin));
      }
    }
    numTail = hTailClean->Integral("width");
  
    return numTail;
  }


  // --------------------------------------------------
  void smearHistogram(const TH1* hOrig, TH1* &hSmeared, double nTotal, double width, double scaling) {
    TString name = hOrig->GetName();
    name += "Smeared";
    hSmeared = static_cast<TH1D*>(hOrig->Clone(name));
    scaling += 1.;
    if( scaling > 1. ) {
      scaling = sqrt( scaling*scaling - 1. )*width;
      hSmeared->Reset();
      for(int bin = 1; bin <= hOrig->GetNbinsX(); ++bin) {
	double entries = hOrig->GetBinContent(bin);
	if( entries ) {
	  double mean = hOrig->GetBinCenter(bin);
	  for(int i = 1; i <= hSmeared->GetNbinsX(); ++i) {
	    double min = hSmeared->GetXaxis()->GetBinLowEdge(i);
	    double max = hSmeared->GetXaxis()->GetBinUpEdge(i);
	    double weight = gaussInt(mean,scaling,min,max)*entries;
	    hSmeared->Fill(hSmeared->GetBinCenter(i),weight);
	  }
	}
      }
      for(int bin = 1; bin <= hSmeared->GetNbinsX(); ++bin) {
	hSmeared->SetBinError(bin,sqrt(hSmeared->GetBinContent(bin)/nTotal));
      }
    } else {
      std::cerr << "WARNING in func::smearHistogram(): scaling = " << scaling << std::endl;
    }
  }
}
#endif

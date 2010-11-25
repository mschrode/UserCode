#ifndef FUNC
#define FUNC

#include "TH1.h"
#include "TF1.h"
#include "TString.h"

namespace func {
  bool fitCoreWidth(const TH1* hist, double nSig, double &width, double &widthErr);
  bool fitCoreWidth(const TH1* hist, double nSig, TF1* &gauss, double &width, double &widthErr);
  double gaussInt(double mean, double sigma, double min, double max);
  double getTail(const TH1* hAsym, double nSig, TH1* &hTail, TH1* &hTailClean, TF1* &gauss);
  double getTailCut(const TH1* hAsym, double cut, TH1* &hTail, TH1* &hTailClean);
  void smearHistogram(const TH1* hOrig, TH1* &hSmeared, double nTotal, double width, double scaling);



  // --------------------------------------------------
  double gaussInt(double mean, double sigma, double min, double max) {
    return 0.5*(erf((max-mean)/sqrt(2.)/sigma) - erf((min-mean)/sqrt(2.)/sigma));
  }


  // --------------------------------------------------
  bool fitCoreWidth(const TH1* hist, double nSig, double &width, double &widthErr) {
    TF1* dummy = 0;
    bool result = fitCoreWidth(hist,nSig,dummy,width,widthErr);
    delete dummy;
    return result;
  }


  // --------------------------------------------------
  bool fitCoreWidth(const TH1* hist, double nSig, TF1* &gauss, double &width, double &widthErr) {
    bool result  = false;

    TString name = hist->GetName();
    name += "_GaussFit";
    gauss = new TF1(name,"gaus",hist->GetXaxis()->GetBinLowEdge(1),hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX()));
    gauss->SetLineWidth(1);
    gauss->SetLineColor(kRed);

    TH1* h = static_cast<TH1*>(hist->Clone("func::fitCoreWidth::h"));

    double mean = h->GetMean();
    double sig = 1.8*h->GetRMS();
    if( h->Fit(gauss,"0QIR","",mean-sig,mean+sig) == 0 ) {
      mean = gauss->GetParameter(1);
      sig = nSig*gauss->GetParameter(2);
      if( h->Fit(gauss,"0QIR","",mean-sig,mean+sig) == 0 ) {
	result = true;
	width = gauss->GetParameter(2);
	widthErr = gauss->GetParError(2);
      }
    } else {
      std::cerr << "WARNING in func::fitCoreWidth: No convergence when fitting width of '" << h->GetName() << "'\n";
      width = 0.;
      widthErr = 10000.;
    }
    delete h;

    return result;
  }


  // --------------------------------------------------
  double getTail(const TH1* hAsym, double nSig, TH1* &hTail, TH1* &hTailClean, TF1* &gauss) {
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
    if( fitCoreWidth(hAsym,nSig,gauss,width,widthErr) ) {
      for(int bin = 1; bin <= hTail->GetNbinsX(); ++bin) {
	double min = hTail->GetXaxis()->GetBinLowEdge(bin);
	double max = hTail->GetXaxis()->GetBinUpEdge(bin);
	double gaussPdf = gauss->Integral(min,max)/hTail->GetBinWidth(1);
	double tailPdf = hTail->GetBinContent(bin) - gaussPdf;
	if( tailPdf < 0. ) tailPdf = 0.;
	hTail->SetBinContent(bin,tailPdf);
	
	if( tailPdf > gaussPdf ) {
	  hTailClean->SetBinContent(bin,hTail->GetBinContent(bin));
	  hTailClean->SetBinError(bin,hTail->GetBinError(bin));
	}
      }
      numTail = hTailClean->Integral("width");
    }
  
    return numTail;
  }


  // --------------------------------------------------
  double getTailFromGauss(const TH1* hAsym, const TF1* extGauss, double tailStart, double nSig, TH1* &hTail, TH1* &hTailClean, TF1* &gauss) {
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
    fitCoreWidth(hAsym,nSig,gauss,width,widthErr);

    // Extract tail using external Gaussian and tail start
    TF1 *extGaussTmp = static_cast<TF1*>(extGauss->Clone("func::getTailFromGauss::extGaussTmp"));
    int optLeftTailStartBin = hTail->FindBin(-1.*std::abs(tailStart));
    int optRightTailStartBin = hTail->FindBin(std::abs(tailStart));
    for(int bin = 1; bin <= hTail->GetNbinsX(); ++bin) {
      double min = hTail->GetXaxis()->GetBinLowEdge(bin);
      double max = hTail->GetXaxis()->GetBinUpEdge(bin);
      double gaussPdf = extGaussTmp->Integral(min,max)/hTail->GetBinWidth(1);
      double tailPdf = hTail->GetBinContent(bin) - gaussPdf;
      if( tailPdf < 0. ) tailPdf = 0.;
      hTail->SetBinContent(bin,tailPdf);
	
      if( bin <= optLeftTailStartBin || bin >= optRightTailStartBin ) {
	hTailClean->SetBinContent(bin,hTail->GetBinContent(bin));
	hTailClean->SetBinError(bin,hTail->GetBinError(bin));
      }
      numTail = hTailClean->Integral("width");
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
    hSmeared->Reset();

    scaling += 1.;
    if( scaling > 1. ) {
      scaling = sqrt( scaling*scaling - 1. )*width;
    } else {
      std::cerr << "WARNING in func::smearHistogram(): scaling = " << scaling << std::endl;
    }
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
  }

}
#endif

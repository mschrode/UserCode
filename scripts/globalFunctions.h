#ifndef FUNC
#define FUNC

#include "TH1.h"
#include "TF1.h"
#include "TString.h"

namespace func {
  // --------------------------------------------------
  double getTail(const TH1* hAsym, TH1* &hTail, TH1* &hTailClean, TF1* &gauss) {
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
  
    name = hAsym->GetName();
    name += "_GaussFit";
    gauss = new TF1(name,"gaus",hAsym->GetXaxis()->GetBinLowEdge(1),hAsym->GetXaxis()->GetBinUpEdge(hAsym->GetNbinsX()));
    gauss->SetLineWidth(1);
    gauss->SetLineColor(kRed);
    if( hAsym->GetFillColor() > 0 ) gauss->SetLineStyle(2);
    double mean = hTail->GetMean();  
    double sigma = 2.5*hTail->GetRMS();
    if( hTail->Fit(gauss,"I0QR","",mean-sigma,mean+sigma) == 0 ) {
      mean = gauss->GetParameter(1);
      sigma = 1.8*gauss->GetParameter(2);
      if( hTail->Fit(gauss,"I0QR","",mean-sigma,mean+sigma) == 0) {
	sigma = gauss->GetParameter(2);
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
    }
  
    return numTail;
  }


  // --------------------------------------------------
  void fitCoreWidth(const TH1 *hist, double &width, double &widthErr) {
    TH1* h = static_cast<TH1*>(hist->Clone("func::fitCoreWidth::h"));
    double mean = h->GetMean();
    double sig = 1.5*h->GetRMS();
    if( h->Fit("gaus","0QIR","",mean-sig,mean+sig) == 0 ) {
      mean = h->GetFunction("gaus")->GetParameter(1);
      sig = 1.8*h->GetFunction("gaus")->GetParameter(2);
      if( h->Fit("gaus","0QIR","",mean-sig,mean+sig) == 0 ) {
	width = h->GetFunction("gaus")->GetParameter(2);
	widthErr = h->GetFunction("gaus")->GetParError(2);
      }
    } else {
      std::cerr << "WARNING in func::fitCoreWidth: No convergence when fitting width of '" << h->GetName() << "'\n";
      width = 0.;
      widthErr = 10000.;
    }
    delete h;
  }


  // --------------------------------------------------
  double gaussInt(double mean, double sigma, double min, double max) {
    return 0.5*(erf((max-mean)/sqrt(2.)/sigma) - erf((min-mean)/sqrt(2.)/sigma));
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

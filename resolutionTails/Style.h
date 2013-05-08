// $Id: $

#ifndef RESOLUTION_TAILS_STYLE
#define RESOLUTION_TAILS_STYLE

#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TString.h"


// Collection of plot-style parameters
// and helper methods to apply style to
// histograms, graphs etc.

namespace resolutionTails {
  class Style {
  public:
    Style();
    ~Style() {};

    // Dimensions
    double binLabelWidth() const { return binLabelWidth_; }
    double legWidth() const { return legWidth_; }
    double pt3PlotMax() const { return pt3PlotMax_; }

    // Plot style parameters
    int colorGauss() const { return colorGauss_; }
    int colorFilledAsym() const { return colorFilledAsym_; }
    int colorFilledAsymSmear() const { return colorFilledAsymSmear_; }
    int colorLineAsymSmear() const { return colorLineAsymSmear_; }
    double markerSize() const { return markerSize_; }
    int lineWidth() const { return lineWidth_; }
    int hatchStyle() const { return hatchStyle_; }

    // Plot title
    bool showTitle() const { return showTitle_; }
    TString title() const { return showTitle_ ? title_ : ""; }

    // Axis and legend labels
    TString labelLumi() const { return labelLumi_; }
    TString labelFAsym() const { return labelFAsym_; }
    TString labelFAsymMC() const { return labelFAsymMC_; }
    TString labelFAsymData() const { return labelFAsymData_; }
    TString labelFAsymGauss() const { return labelFAsymGauss_; }
    TString labelFAsymToy() const { return labelFAsymToy_; }
    TString labelScaleFactor() const { return labelScaleFactor_; }
    TString labelAsymTailStart() const { return labelAsymTailStart_; }
    TString labelAsymTailStartEff() const { return labelAsymTailStartEff_; }
    TString labelCommonSigma() const { return labelCommonSigma_; }
    TString labelData() const { return labelData_; }
    TString labelMC() const { return labelMC_; }
    TString labelMCSmear() const { return labelMCSmear_; }
    TString labelMCGauss() const { return labelMCGauss_; }
    TString labelWindow(double nSigMin, double nSigMax) const;

    // Apply style attributes
    void applyToDataMarker(TH1* h) const;
    void applyToMCFilled(TH1* h) const;
    void applyToData(TGraphAsymmErrors* g) const;
    void applyToMC(TGraphAsymmErrors* g) const;

    // Strings for file names
    TString nameWindow(double nSigMin, double nSigMax) const;

    // Setters
    void showTitle(bool val) { showTitle_ = val; }
    void setTitle(const TString &val) { title_ = val; }
    void setLabelLumi(const TString &val) { labelLumi_ = val; }
    void setLabelData(const TString &val) { labelData_ = val; }
    void setLabelMC(const TString &val) { labelMC_ = val; }
    void setLabelMCSmear(const TString &val) { labelMCSmear_ = val; }


  private:
    // Dimensions
    double binLabelWidth_;
    double legWidth_;
    double pt3PlotMax_;

    // Plot style parameters
    int colorGauss_;
    int colorFilledAsym_;
    int colorFilledAsymSmear_;
    int colorLineAsymSmear_;
    double markerSize_;
    int lineWidth_;
    int hatchStyle_;

    // Plot title
    bool showTitle_;
    TString title_;

    // Axis and legend labels
    TString labelLumi_;
    TString labelFAsym_;
    TString labelFAsymMC_;
    TString labelFAsymData_;
    TString labelFAsymGauss_;
    TString labelFAsymToy_;
    TString labelScaleFactor_;
    TString labelAsymTailStart_;
    TString labelAsymTailStartEff_;
    TString labelCommonSigma_;
    TString labelData_;
    TString labelMC_;
    TString labelMCSmear_;
    TString labelMCGauss_;
  };


  // ------------------------------------------------------------------------------------
  Style::Style() {
    binLabelWidth_ = -0.52;
    legWidth_ = 0.48;
    pt3PlotMax_ = 0.23;

    // Plot style parameters
    colorGauss_ = 46;
    colorFilledAsym_ = 38;
    colorFilledAsymSmear_ = 29;
    colorLineAsymSmear_ = 30;
    markerSize_ = 1.4;
    lineWidth_ = 2;
    hatchStyle_ = 3354;

    // Plot title
    showTitle_ = false;
    title_ = "";

    // Axis and legend labels
    labelLumi_ = "";
    labelFAsym_ = "f_{asym}";
    labelFAsymMC_ = "f^{mc}_{asym}";
    labelFAsymData_ = "f^{data}_{asym}";
    labelFAsymGauss_ = "f^{gauss}_{asym}";
    labelFAsymToy_ = "f^{toy}_{asym}";
    labelScaleFactor_ = "#rho_{tail}";
    labelAsymTailStart_ = "A_{tail}";
    labelAsymTailStartEff_ = "#hat{A}_{tail}";
    labelCommonSigma_ = "#sigma_{c}";
    labelData_ = "Data";
    labelMC_ = "MC";
    labelMCSmear_ = "MCSmear";
    labelMCGauss_ = "Gaussian Asymmetry";
  }


  // Label for window in asymmetry from 'nSigMin' to 'nSigMax'
  // to be used in legends etc.
  // ------------------------------------------------------------------------------------
  TString Style::labelWindow(double nSigMin, double nSigMax) const {
    TString label = "";
    char tmp[100];
    if( nSigMax < 50. ) {
      sprintf(tmp,"%.1f - %.1f",nSigMin,nSigMax);
      label += tmp;
      label += " "+labelCommonSigma();
    } else {
      label = labelAsymTailStart()+" = ";
      sprintf(tmp,"%.1f",nSigMin);
      label += " "+labelCommonSigma();
    }

    return label;
  }


  // Name for window in asymmetry from 'nSigMin' to 'nSigMax'
  // to be used in file names, plot names etc.
  // ------------------------------------------------------------------------------------
  TString Style::nameWindow(double nSigMin, double nSigMax) const {
    TString name = "";
    char tmp[50];
    sprintf(tmp,"%.1f-",nSigMin);
    name += tmp;
    if( nSigMax > 10. ) {
      name += "Inf";
    } else {
      sprintf(tmp,"%.1f",nSigMax);
      name += tmp;
    }
    name.ReplaceAll(".","");

    return name;
  }



  // ------------------------------------------------------------------------------------
  void Style::applyToDataMarker(TH1* h) const {
    h->SetMarkerSize(markerSize());
    h->SetLineWidth(lineWidth());
    h->SetMarkerStyle(20);
    h->SetLineColor(kBlack);
    h->SetMarkerColor(kBlack);
  }

  
  // ------------------------------------------------------------------------------------
  void Style::applyToMCFilled(TH1* h) const {
    h->SetLineWidth(1);
    h->SetMarkerStyle(1);
    TString name = h->GetName();
    if( name.Contains("Smear") )
      h->SetFillColor(colorFilledAsymSmear());
    else 
      h->SetFillColor(colorFilledAsym());
  }

  
  // ------------------------------------------------------------------------------------
  void Style::applyToData(TGraphAsymmErrors* g) const {
    g->SetMarkerSize(markerSize());
    g->SetLineWidth(lineWidth());
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlack);
    g->SetLineColor(g->GetMarkerColor());
    g->SetLineStyle(2);
  }

  
  // ------------------------------------------------------------------------------------
  void Style::applyToMC(TGraphAsymmErrors* g) const {
    g->SetMarkerSize(markerSize());
    g->SetLineWidth(lineWidth());
    g->SetMarkerStyle(21);
    g->SetMarkerColor(colorFilledAsymSmear());
    g->SetLineColor(g->GetMarkerColor());
  }
}
#endif

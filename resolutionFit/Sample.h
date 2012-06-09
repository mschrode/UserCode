// $Id: Sample.h,v 1.19 2012/06/08 21:14:44 mschrode Exp $

#ifndef SAMPLE_H
#define SAMPLE_H

#include <map>
#include <set>
#include <vector>

#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TString.h"

#include "FitResult.h"
#include "Measurement.h"


namespace resolutionFit {

  typedef TString SampleLabel;

  // -------------------------------------------------------------------------------------
  class Sample {
  public:
    enum Type { Data, MC, MCTruth, NONE };

    static bool validType(Type type);
    static TString toString(Type type);
    static int color(const SampleLabel &label);
    static int markerStyle(const SampleLabel &label);
    static Sample::Type type(const SampleLabel &label);

    static void setOnlySolidMarkerStyles();


    Sample(const SampleLabel &label, unsigned int verbosity = 0);
    virtual ~Sample();

    virtual Sample::Type type() const = 0;
    virtual double relativeWeightTo(const SampleLabel &other, unsigned int ptSoftBin) const;

    TString label() const { return label_; }
    TString printLabel() const { return printLabel_; }
    int color() const { return color(label()); }
    int markerStyle() const { return markerStyle(label()); }

    TH1* histPtAsym(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histPtAsym(); }
    TH1* histPtGenAsym(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histPtGenAsym(); }
    TH1* histPtGen(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histPtGen(); }
    TH1* histPdfPtTrue(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histPdfPtTrue(); }
    TH1* histPtAve(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histPtAve(); }
    TH1* histPt1(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histPt1(); }
    TH1* histPt2(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histPt2(); }
    TH1* histPt3(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histPt3(); }
    TH1* histNumVtx(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histNumVtx(); }
    TH1* histMCWeight(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histMCWeight(); }
    TH1* histMCNumPU(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histMCNumPU(); }
    TH1* histDeltaPt(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histDeltaPt(); }

    void setPrintLabel(const TString &printLabel) { printLabel_ = printLabel; }
    bool addFitResult(FitResult::Type type, double minPt3, int wpIdx);
    bool setKSoftFit(FitResult::Type type, const TF1* fit);
    void addSystematicUncertainty(FitResult::Type type, const TString &label, double variedValue, double fraction);
    void addSystematicUncertainty(FitResult::Type type, const TString &label, double variedValueDown, double variedValueUp, double fraction);
    void setRelativeWeightTo(const SampleLabel &other, unsigned int ptSoftBin, double relWeight);


    unsigned int nPtSoftBins() const { return meas_.size(); }
    double ptSoft(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->ptSoft(); }
    void ptSoft(std::vector<double> &ptSoftVals) const;

    double asymmetryWidth(unsigned int ptSoftBin) const { return asymmetryWidths_.at(ptSoftBin); }
    double asymmetryWidthStatUncert(unsigned int ptSoftBin) const { return asymmetryWidthErrs_.at(ptSoftBin); }
    double meanPtAve(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->meanPtAve(); }
    double meanPtAveStatUncert(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->meanPtAveUncert(); }

    double meanPt(FitResult::Type type) const;
    double meanPtStatUncert(FitResult::Type type) const;
    void valuesInExtrapolation(FitResult::Type type, std::vector<double> &val, std::vector<double> &uncert) const;
    void fittedValues(FitResult::Type type, std::vector<double> &val, std::vector<double> &uncert) const;
    unsigned int firstPointInExtrapolation(FitResult::Type type) const;
    double valueInExtrapolation(FitResult::Type type, unsigned int ptSoftBin) const;
    double uncertInExtrapolation(FitResult::Type type, unsigned int ptSoftBin) const;
    double fittedValue(FitResult::Type type, unsigned int ptSoftBin) const;
    double fittedUncert(FitResult::Type type, unsigned int ptSoftBin) const;
    double kSoftSlope(FitResult::Type type) const;
    double kSoftSlopeStatUncert(FitResult::Type type) const;
    TF1* kSoftFit(FitResult::Type type, const TString &name) const;
    TF1* extrapolationFunction(FitResult::Type type, const TString &name) const;
    double extrapolatedValue(FitResult::Type type) const;
    double extrapolatedStatUncert(FitResult::Type type) const;
    double extrapolatedSystUncert(FitResult::Type type) const;
    TString labelQuantityInExtrapolation(FitResult::Type type) const;


  protected:
    static std::map<TString,int> COLOR;
    static std::map<TString,int> MARKER_STYLE;
    static std::map<TString,Sample::Type> TYPE;

    typedef std::map<FitResult::Type,FitResult*> FitResultMap;
    typedef std::map<FitResult::Type,FitResult*>::const_iterator FitResultMapIt;

    const TString label_;
    const unsigned int verbosity_;

    TString printLabel_;
    Meas meas_;
    FitResultMap fitResult_;
    std::map<SampleLabel, std::vector<double> > relWeightToOtherSample_;

    int color(unsigned int idx) const;
    bool findFitResult(FitResult::Type type, FitResult* &fitResult) const;
    void fitAsymmetryWidths();
    void cureStatUncerts();


  private:
    std::vector<double> asymmetryWidths_; //! Standard deviation of Gaussian fitted to asymmetry histogram for the different ptSoft bins
    std::vector<double> asymmetryWidthErrs_;  //! Uncertainty on standard deviation of Gaussian fitted to asymmetry histogram for the different ptSoft bins
  };

  typedef std::map<SampleLabel,Sample::Type> SampleTypes;
  typedef std::map<SampleLabel,Sample::Type>::const_iterator SampleTypeIt;
  typedef std::map<SampleLabel,Sample*> Samples;
  typedef std::map<SampleLabel,Sample*>::const_iterator SampleIt;
  

  class SampleLabelPair {
  public:
    SampleLabelPair(const SampleLabel &label1, const SampleLabel &label2)
      : label1_(label1), label2_(label2) {}    

    bool operator==(const SampleLabelPair &other) const {
      return ( (this->label1() == other.label1()) && (this->label2() == other.label2()) );
    }
    bool operator<(const SampleLabelPair &other) const {
      TString thisLabel = this->label1()+this->label2();
      TString otherLabel = other.label1()+other.label2();
      return thisLabel < otherLabel;
    }

    SampleLabel label1() const { return label1_; }
    SampleLabel label2() const { return label2_; }
    bool contains(const SampleLabel &label) const {
      return (label == label1()) || (label == label2());
    }

  private:
    const SampleLabel label1_;
    const SampleLabel label2_;
  };

  typedef std::set<SampleLabelPair*> ComparedSamples;
  typedef std::set<SampleLabelPair*>::const_iterator ComparedSamplesIt;  



  // -------------------------------------------------------------------------------------
  class DataSample : public Sample {
  public:
    DataSample(const TString &label, unsigned int etaBin, double etaMin, double etaMax, unsigned int ptBin, double ptMin, double ptMax, const std::vector<double> &ptSoft, const TString &fileName, unsigned int verbosity);
    
    Sample::Type type() const { return Data; }
    double relativeWeightTo(const SampleLabel &other) const { return 1.; }


  private:
    static unsigned int N_DATA_SAMPLES;
  };

  typedef std::map<SampleLabel,DataSample*> DataSamples;
  typedef std::map<SampleLabel,DataSample*>::const_iterator DataSampleIt;




  // -------------------------------------------------------------------------------------
  class MCSample : public Sample {
  public:
    MCSample(const TString &label, unsigned int etaBin, double etaMin, double etaMax, unsigned int ptBin, double ptMin, double ptMax, const std::vector<double> &ptSoft, const TString &fileName, unsigned int verbosity);
    
    Sample::Type type() const { return MC; }


  private:
    static unsigned int N_MC_SAMPLES;
  };

  typedef std::map<SampleLabel,MCSample*> MCSamples;
  typedef std::map<SampleLabel,MCSample*>::const_iterator MCSampleIt;



  // -------------------------------------------------------------------------------------
  class MCTruthSample : public Sample {
  public:
    MCTruthSample(const TString &label, unsigned int etaBin, double etaMin, double etaMax, unsigned int ptBin, double ptMin, double ptMax, const std::vector<double> &ptSoft, const TString &fileName, unsigned int verbosity);
    
    Sample::Type type() const { return MCTruth; }


  private:
    static unsigned int N_MCTRUTH_SAMPLES;
  };
}
#endif

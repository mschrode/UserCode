// $Id: Sample.h,v 1.7 2011/02/28 10:53:15 mschrode Exp $

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

  // -------------------------------------------------------------------------------------
  class Sample {
  public:
    enum Type { Data, MC, MCTruth };

    static bool validType(Type type);
    static TString toString(Type type);

    Sample(const TString &label, const std::vector<TString> &fileNames, const std::vector<double> &ptSoft, const TString &histNameSuffix, unsigned int verbosity = 0);
    virtual ~Sample();

    virtual Sample::Type type() const = 0;

    TString label() const { return label_; }

    TH1* histPtAsym(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histPtAsym(); }
    TH1* histPtGenAsym(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histPtGenAsym(); }
    TH1* histPtGen(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histPtGen(); }
    TH1* histPdfPtTrue(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histPdfPtTrue(); }
    TH1* histPtAve(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histPtAve(); }
    TH1* histPt1(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histPt1(); }
    TH1* histPt2(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histPt2(); }
    TH1* histPt3(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->histPt3(); }

    bool addFitResult(FitResult::Type type);
    void addSystematicUncertainty(FitResult::Type type, const TString &label, double variedValue, double fraction);
    void addSystematicUncertainty(FitResult::Type type, const TString &label, double variedValueDown, double variedValueUp, double fraction);

    unsigned int nPtSoftBins() const { return meas_.size(); }
    double ptSoft(unsigned int ptSoftBin) const { return meas_.at(ptSoftBin)->ptSoft(); }
    void ptSoft(std::vector<double> &ptSoftVals) const;
    double meanPt(FitResult::Type type) const;
    double meanPtStatUncert(FitResult::Type type) const;
    void values(FitResult::Type type, std::vector<double> &val, std::vector<double> &uncert) const;
    double fittedValue(FitResult::Type type, unsigned int ptSoftBin) const;
    double fittedUncert(FitResult::Type type, unsigned int ptSoftBin) const;
    TF1* extrapolationFunction(FitResult::Type type, const TString &name) const;
    double extrapolatedValue(FitResult::Type type) const;
    double extrapolatedStatUncert(FitResult::Type type) const;
    double extrapolatedSystUncert(FitResult::Type type) const;


  protected:
    typedef std::map<FitResult::Type,FitResult*> FitResultMap;
    typedef std::map<FitResult::Type,FitResult*>::const_iterator FitResultMapIt;

    const TString label_;
    const unsigned int verbosity_;

    Meas meas_;
    FitResultMap fitResult_;

    bool findFitResult(FitResult::Type type, FitResult* &fitResult) const;
  };

  typedef TString SampleLabel;
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

  private:
    const SampleLabel label1_;
    const SampleLabel label2_;
  };

  typedef std::set<SampleLabelPair*> ComparedSamples;
  typedef std::set<SampleLabelPair*>::const_iterator ComparedSamplesIt;  



  // -------------------------------------------------------------------------------------
  class DataSample : public Sample {
  public:
    DataSample(const TString &label, const std::vector<TString> &fileNames, const std::vector<double> &ptSoft, const TString &histNameSuffix, unsigned int verbosity = 0);
    
    Sample::Type type() const { return Data; }
  };

  typedef std::map<SampleLabel,DataSample*> DataSamples;
  typedef std::map<SampleLabel,DataSample*>::const_iterator DataSampleIt;




  // -------------------------------------------------------------------------------------
  class MCSample : public Sample {
  public:
    MCSample(const TString &label, const std::vector<TString> &fileNames, const std::vector<double> &ptSoft, const TString &histNameSuffix, unsigned int verbosity = 0);
    
    Sample::Type type() const { return MC; }
  };

  typedef std::map<SampleLabel,MCSample*> MCSamples;
  typedef std::map<SampleLabel,MCSample*>::const_iterator MCSampleIt;



  // -------------------------------------------------------------------------------------
  class MCTruthSample : public Sample {
  public:
    MCTruthSample(const TString &label, const std::vector<TString> &fileNames, const std::vector<double> &ptSoft, const TString &histNameSuffix, unsigned int verbosity = 0);
    
    Sample::Type type() const { return MCTruth; }
  };
}
#endif

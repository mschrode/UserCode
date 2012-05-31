// $ Id: $

#ifndef PT_BIN_H
#define PT_BIN_H

#include <vector>

#include "TString.h"

#include "FitResult.h"
#include "Parameters.h"
#include "Sample.h"


namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  class PtBin {
  public:
    PtBin(unsigned int etaBin, unsigned int ptBin, const Parameters* par);
    ~PtBin();

    unsigned int etaBin() const { return etaBin_; }
    unsigned int ptBin() const { return ptBin_; }
    TString toTString() const;
    double min() const { return par_->ptMin(etaBin(),ptBin()); }
    double max() const { return par_->ptMax(etaBin(),ptBin()); }

    SampleIt samplesBegin() const { return samples_.begin(); }
    SampleIt samplesEnd() const { return samples_.end(); }
    DataSampleIt dataSamplesBegin() const { return dataSamples_.begin(); }
    DataSampleIt dataSamplesEnd() const { return dataSamples_.end(); }
    MCSampleIt mcSamplesBegin() const { return mcSamples_.begin(); }
    MCSampleIt mcSamplesEnd() const { return mcSamples_.end(); }
    const Sample* findSample(const SampleLabel &label) const;
    bool hasMCTruthSample() const { return mcTruthSample_; }
    const Sample* mcTruthSample() const { return mcTruthSample_; }
    
    bool addDataSample(const TString &label, const TString &fileName, const TString &printLabel) {
      return addSample(Sample::Data,label,fileName,printLabel);
    }
    bool addMCSample(const TString &label, const TString &fileName, const TString &printLabel) {
      return addSample(Sample::MC,label,fileName,printLabel);
    }
    bool addMCTruthSample(const TString &label, const TString &fileName);
    bool addFitResult(FitResult::Type type);
    bool setKSoftFit(const SampleLabel &label, FitResult::Type type, const TF1* fit);

    void weightSampleRelTo(const SampleLabel &sampleToBeWeighted, const SampleLabel &referenceSample);
    

  private:
    const Parameters* par_;
    const unsigned int etaBin_;
    const unsigned int ptBin_;

    DataSamples dataSamples_;
    MCSamples mcSamples_;
    Sample* mcTruthSample_;
    Samples samples_;

    bool addSample(Sample::Type type, const TString &label, const TString &fileName, const TString &printLabel);
  };

  typedef std::vector<PtBin*> PtBins;
  typedef std::vector<PtBin*>::const_iterator PtBinIt;
}
#endif

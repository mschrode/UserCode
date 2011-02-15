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

    SampleIt samplesBegin() const { return samples_.begin(); }
    SampleIt samplesEnd() const { return samples_.end(); }
    DataSampleIt dataSamplesBegin() const { return dataSamples_.begin(); }
    DataSampleIt dataSamplesEnd() const { return dataSamples_.end(); }
    MCSampleIt mcSamplesBegin() const { return mcSamples_.begin(); }
    MCSampleIt mcSamplesEnd() const { return mcSamples_.end(); }
    const Sample* findSample(SampleLabel label) const;
    
    bool addDataSample(const TString &label, const TString &baseFileName) {
      return addSample(Sample::Data,label,baseFileName);
    }
    bool addMCSample(const TString &label, const TString &baseFileName) {
      return addSample(Sample::MC,label,baseFileName);
    }
    bool addFitResult(FitResult::Type type);
    

  private:
    const Parameters* par_;
    const unsigned int etaBin_;
    const unsigned int ptBin_;

    DataSamples dataSamples_;
    MCSamples mcSamples_;
    Samples samples_;

    bool addSample(Sample::Type type, const TString &label, const TString &baseFileName);
  };

  typedef std::vector<PtBin*> PtBins;
  typedef std::vector<PtBin*>::const_iterator PtBinIt;
}
#endif

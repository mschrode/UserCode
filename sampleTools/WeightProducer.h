// $Id: $

#ifndef WEIGHT_PRODUCER_H
#define WEIGHT_PRODUCER_H

#include <iostream>
#include <vector>

#include "TH1.h"

// Produce weights for spectrum and PU reweighting
//
// PU scenarios from: https://twiki.cern.ch/twiki/bin/view/CMS/Pileup_MC_Gen_Scenarios
// Code adapted from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting


namespace sampleTools {

class WeightProducer {
public:
  // Different QCD samples
  // - Flat: Generic flat sample, supposed to be flattened with (ptHat/15)^4.5
  // - Summer11PtHatBinned: /QCD_Pt-*to*_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM
  enum QCDSample { Flat, Summer11PtHatBinned };
  enum PUScenario { Flat10, Fall11, Summer12S7, Summer12S10 };

  static void init(QCDSample sa, PUScenario sc, const TH1* data_npu_estimated);
  static double puWeight(unsigned int nPUVtx);
  static double qcdSpectrumWeight(double ptHat);

private:
  static QCDSample qcdSample_;
  static PUScenario puScenario_;
  static std::vector<double> puWeights_;
};

std::vector<double> WeightProducer::puWeights_ = std::vector<double>(0);
WeightProducer::QCDSample WeightProducer::qcdSample_ = WeightProducer::Flat;
WeightProducer::PUScenario WeightProducer::puScenario_ = WeightProducer::Flat10;



// --------------------------------------------------
void WeightProducer::init(QCDSample sa, PUScenario sc, const TH1* data_npu_estimated) {
  WeightProducer::qcdSample_ = sa;
  WeightProducer::puScenario_ = sc;

  unsigned int nMaxPU = 0;
  double *npuProbs = 0;

  if( WeightProducer::puScenario_ == Flat10 ) {
    nMaxPU = 25;
    // see SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py; copy and paste from there:
    double npu_probs[25] = { 0.0698146584,
			     0.0698146584,
			     0.0698146584,
			     0.0698146584,
			     0.0698146584,
			     0.0698146584,
			     0.0698146584,
			     0.0698146584,
			     0.0698146584,
			     0.0698146584,
			     0.0698146584, /* <-- 10*/
			     0.0630151648,
			     0.0526654164,
			     0.0402754482,
			     0.0292988928,
			     0.0194384503,
			     0.0122016783,
			     0.007207042,
			     0.004003637,
			     0.0020278322,
			     0.0010739954,
			     0.0004595759,
			     0.0002229748,
			     0.0001028162,
			     4.58337152809607E-05 /* <-- 24 */};
    npuProbs = npu_probs;
  } else if( WeightProducer::puScenario_ == Fall11 ) {
    nMaxPU = 50;
    double npu_probs[50] = {
      0.003388501,
      0.010357558,
      0.024724258,
      0.042348605,
      0.058279812,
      0.068851751,
      0.072914824,
      0.071579609,
      0.066811668,
      0.060672356,
      0.054528356,
      0.04919354,
      0.044886042,
      0.041341896,
      0.0384679,
      0.035871463,
      0.03341952,
      0.030915649,
      0.028395374,
      0.025798107,
      0.023237445,
      0.020602754,
      0.0180688,
      0.015559693,
      0.013211063,
      0.010964293,
      0.008920993,
      0.007080504,
      0.005499239,
      0.004187022,
      0.003096474,
      0.002237361,
      0.001566428,
      0.001074149,
      0.000721755,
      0.000470838,
      0.00030268,
      0.000184665,
      0.000112883,
      6.74043E-05,
      3.82178E-05,
      2.22847E-05,
      1.20933E-05,
      6.96173E-06,
      3.4689E-06,
      1.96172E-06,
      8.49283E-07,
      5.02393E-07,
      2.15311E-07,
      9.56938E-08
    };
    npuProbs = npu_probs;
  } else if( WeightProducer::puScenario_ == Summer12S7 ) {
    nMaxPU = 50;
    double npu_probs[50] = {
      0.003388501,
      0.010357558,
      0.024724258,
      0.042348605,
      0.058279812,
      0.068851751,
      0.072914824,
      0.071579609,
      0.066811668,
      0.060672356,
      0.054528356,
      0.04919354,
      0.044886042,
      0.041341896,
      0.0384679,
      0.035871463,
      0.03341952,
      0.030915649,
      0.028395374,
      0.025798107,
      0.023237445,
      0.020602754,
      0.0180688,
      0.015559693,
      0.013211063,
      0.010964293,
      0.008920993,
      0.007080504,
      0.005499239,
      0.004187022,
      0.003096474,
      0.002237361,
      0.001566428,
      0.001074149,
      0.000721755,
      0.000470838,
      0.00030268,
      0.000184665,
      0.000112883,
      6.74043E-05,
      3.82178E-05,
      2.22847E-05,
      1.20933E-05,
      6.96173E-06,
      3.4689E-06,
      1.96172E-06,
      8.49283E-07,
      5.02393E-07,
      2.15311E-07,
      9.56938E-08
    };
    npuProbs = npu_probs;
  } else if( WeightProducer::puScenario_ == Summer12S10 ) {
    nMaxPU = 60;
    double npuSummer12_S10[60] = {
      2.560E-06,
      5.239E-06,
      1.420E-05,
      5.005E-05,
      1.001E-04,
      2.705E-04,
      1.999E-03,
      6.097E-03,
      1.046E-02,
      1.383E-02,
      1.685E-02,
      2.055E-02,
      2.572E-02,
      3.262E-02,
      4.121E-02,
      4.977E-02,
      5.539E-02,
      5.725E-02,
      5.607E-02,
      5.312E-02,
      5.008E-02,
      4.763E-02,
      4.558E-02,
      4.363E-02,
      4.159E-02,
      3.933E-02,
      3.681E-02,
      3.406E-02,
      3.116E-02,
      2.818E-02,
      2.519E-02,
      2.226E-02,
      1.946E-02,
      1.682E-02,
      1.437E-02,
      1.215E-02,
      1.016E-02,
      8.400E-03,
      6.873E-03,
      5.564E-03,
      4.457E-03,
      3.533E-03,
      2.772E-03,
      2.154E-03,
      1.656E-03,
      1.261E-03,
      9.513E-04,
      7.107E-04,
      5.259E-04,
      3.856E-04,
      2.801E-04,
      2.017E-04,
      1.439E-04,
      1.017E-04,
      7.126E-05,
      4.948E-05,
      3.405E-05,
      2.322E-05,
      1.570E-05,
      5.005E-06};
    npuProbs = npuSummer12_S10;
  }

  WeightProducer::puWeights_ = std::vector<double>(nMaxPU);
  double s = 0.0;
  for(unsigned int npu = 0; npu < nMaxPU; ++npu) {
    double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));
    WeightProducer::puWeights_[npu] = npu_estimated / npuProbs[npu];
    s += npu_estimated;
  }
  // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
  for (unsigned int npu = 0; npu < nMaxPU; ++npu) {
    WeightProducer::puWeights_[npu] /= s;
  }
}

  
double WeightProducer::puWeight(unsigned int nPUVtx) {
  double w = 1.;
  if( nPUVtx < WeightProducer::puWeights_.size() ) {
    w = WeightProducer::puWeights_.at(nPUVtx);
  } else {
    std::cerr << "WARNING in WeightProducer: no event weights for nPUVtx = " << nPUVtx << " - using weights for nPUVtx = " << WeightProducer::puWeights_.size()-1 << " instead" <<  std::endl;
    w = WeightProducer::puWeights_.back();
  }

  return w;
}



// Get weight for ptHat binned MC samples
//
// The weights returned are hard-coded for the
// Summer11 PYTHIA QCD samples Tune Z2,
// /QCD_Pt-*to*_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM
// --------------------------------------------------
double WeightProducer::qcdSpectrumWeight(double ptHat) {
  double weight = 1.;

  if( WeightProducer::qcdSample_ == WeightProducer::Flat ) {
    weight = pow((ptHat/15.),-4.5);
  } else if( WeightProducer::qcdSample_ == WeightProducer::Summer11PtHatBinned ) {
    double lumi = 1.;
    unsigned int numEvts = 0;
    double xs = 0.;
    double scale = 1.;		// To take into account missing jobs
    if( ptHat < 15. ) {
      numEvts   = 1650000;
      xs        = 3.675e+10;
    } else if( ptHat < 30. ) {
      numEvts   = 11000000;
      xs        = 8.159e+08;
      scale = 382./381.;
    } else if( ptHat < 50. ) {
      numEvts   = 6583068;
      xs        = 5.312e+07;
    } else if( ptHat < 80. ) {
      numEvts   = 6600000;
      xs        = 6.359e+06;
      scale = 233./232.;
    } else if( ptHat < 120. ) {
      numEvts   = 6589956;
      xs        = 7.843e+05;
    } else if( ptHat < 170. ) {
      numEvts   = 6127528;
      xs        = 1.151e+05;
    } else if( ptHat < 300. ) {
      numEvts   = 6220160;
      xs        = 2.426e+04;
      scale = 210./209.;
    } else if( ptHat < 470. ) {
      numEvts   = 6432669;
      xs        = 1.168e+03;
      scale = 217./215.;
    } else if( ptHat < 600. ) {
      numEvts   = 3990085;
      xs        = 7.022e+01;
      scale     = 134./129.;
    } else if( ptHat < 800. ) {
      numEvts   = 4245695;
      xs        = 1.555e+01;
    } else if( ptHat < 1000. ) {
      numEvts   = 4053888;
      xs        = 1.844e+00;
      scale     = 137./134.;
    } else if( ptHat < 1400. ) {
      numEvts   = 2093222;
      xs        = 3.321e-01;
    } else if( ptHat < 1800. ) {
      numEvts   = 2196200;
      xs        = 1.087e-02;
      scale     = 74./73.;
    } else {
      numEvts   = 293139;
      xs        = 3.575e-04;
    }
    weight = numEvts > 0 ? scale * lumi * xs / numEvts : 0.;
  }     

  return weight;
}

}
#endif

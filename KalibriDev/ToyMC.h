#ifndef TOYMC_H
#define TOYMC_H

#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TRandom.h"

class TTree;


//!  \brief Generate toy MC data
//!
//!  \author Hartmut Stadie
//!  \date   Mon Jun 30 11:00:00 CEST 2008
//!  $Id: ToyMC.h,v 1.1 2009/10/30 08:59:45 mschrode Exp $
// ----------------------------------------------------------------  
class ToyMC {

 private:
  //!  \brief Resolution model
  //! 
  //!  - For calibration:
  //!    - 'Gauss': An energy dependent Gaussian resolution
  //!      \f[
  //!       \frac{\sigma}{E} = \frac{b}{E} \oplus \frac{b}{\sqrt{E}} \oplus b
  //!      \f]
  //!      
  //!    - 'Landau': An energy dependent Landau resolution
  //!      
  //!    - 'Dirac': A Dirac delta function as perfect resolution
  //!
  //!  - For jet smearing:
  //!    - 'GaussUniform': A Gaussian resolution with a uniform
  //!      low energy tail.
  //!      The appropriate parametrizations are
  //!      - SmearParametrizationFermiTail
  //!  
  //!    - 'TwoGauss': Resolution given by a central Gaussian around
  //!      1 and a second Gaussian to model tails.
  //!      The appropriate parametrizations are
  //!      - SmearParametrizationStepGauss
  // ----------------------------------------------------------------  
  enum ResolutionModel { Gauss, Landau, GaussUniform, TwoGauss, Dirac };


  //!  \brief Response for generation
  //!
  //!  The generated response
  //!  \f[
  //!   R(E^{\textrm{true}}_{T}) = \frac{E^{\textrm{meas}}_{T}}{E^{\textrm{true}}_{T}}
  //!  \f]
  //!  depends on the parameters \f$ A_{i} \f$ as defined in the config file
  //!  via the field 'ToyMC response parameters'. The response model is definded
  //!  via the field 'ToyMC response model'.
  //!
  //!  Possible response models are:
  //!    - 'Constant': Apply a constant response factor to the
  //!      hadronic part of the jet Et:
  //!      \f[ R = \frac{1}{A_{0}} \f]
  //!      This is the default.
  //!      The appropriate correction functions are
  //!      - ToyParametrization
  //!      - ToyJetParametrization
  //!      - ToyStepParametrization
  //!      - ToyStepJetParametrization
  //!
  //!    - 'L3': Level 3 response from JetMET group. In each event,
  //!      the response factor is calculated from the total photon pt
  //!      (i.e. \f$ E^{\textrm{true}}_{T} = P^{\gamma}_{T} \f$)
  //!      but applied only to the hadronic part of the jet Et:
  //!      \f[ R = A_{0}
  //!            - \frac{A_{1}}{\log^{A_{2}}_{10}(E^{\textrm{true}}_{T}) + A_{3}}
  //!            +\frac{A_{4}}{E^{\textrm{true}}_{T}}  \f]
  //!      The appropriate correction functions are
  //!      - L2L3JetParametrization
  //!      - L2L3JetParametrization2
  //!   
  //!    - 'SimpleInverse': In each event,
  //!      the response factor is calculated from the total photon pt
  //!      (i.e. \f$ E^{\textrm{true}}_{T} = P^{\gamma}_{T} \f$)
  //!      but applied only to the hadronic part of the jet Et:
  //!      \f[ R = 1 - \frac{A_{0}}{E^{\textrm{true}}_{T} + A_{1}} \f]
  //!      The appropriate correction function is
  //!      - ToySimpleInverseParametrization
  //!
  //!    - 'Flat'
  //!    - 'Exp'
  //!    - 'Slope'
  //!
  //!    - 'StepEta': The jet response is given by
  //!      \f[ R(\eta) = A_{i},
  //!          \qquad i = 0 \textrm{ for } \eta < 0,
  //!          i = 1 \textrm{ else} \f]
  //!      The appropriate correction function is
  //!      - ToyJetParametrization
  //!    - 'SinusEta': The jet response is given by
  //!      \f[ R(\eta) = 1 + A_{0} \sin(A_{1} \eta) \f]
  //!      The appropriate correction function is
  //!      - ToyJetParametrization
  //!    - 'SinusEtaSimpleInversePt': The jet response is given by
  //!      the product of the response models 'SinusEta' and
  //!      'SimpleInverse'
  // ----------------------------------------------------------------  
  enum ResponseModel { Constant, L3, SimpleInverse, Flat, Slope, Exp, StepEta, SinusEta, SinusEtaSimpleInversePt };


  //!  \brief Truth pt spectrum
  // ----------------------------------------------------------------  
  enum TruthSpectrum { Uniform, PowerLaw, PtEtaHistogram };


  // Global variables
  TRandom*        random_;             //!< Random generator
  int             type_;               //!< Event type: Photonjet (1), Dijet (2)

  // Parameters for truth
  double          minEta_;             //!< Minimum truth eta
  double          maxEta_;             //!< Maximum truth eta
  double          minPt_;              //!< Minimum truth pt
  double          maxPt_;              //!< Maximum truth pt
  TruthSpectrum   ptSpectrum_;         //!< Truth pt spectrum
  std::vector<double> parTruth_;       //!< Parameters for truth spectrum
  TLorentzVector  pInput_;             //!< Stores the truth lorentz vector of the current event

  TH2F          * histPtEta_;          //!< For histogramed truth spectrum

  // Parameters for measurement 
  int             chunks_;
  double          jetSpreadA_;
  double          jetSpreadB_;
  bool            useTowerCenterEtaPhi_;
  bool            noOutOfCone_;
  double          maxPi0Frac_;
  double          maxEmf_;

  ResponseModel responseModel_;        //!< Response models
  std::vector<double> parResp_;        //!< Parameters for Response
  TH1F          * histResp_;           //!< For histogramed response

  ResolutionModel resolutionModel_;    //!< Resolution model
  std::vector<double> parReso_;        //!< Parameters for Respolution

  double          smearFactor_;        //!< Combined smear factor from response and resolution
  bool            smearTowersIndividually_;  //!< If true, mSmearTowersIndividually is determined individually for each tower, else for each jet


  void genInput();
  void calIds(float& eta, float &phi, int& ieta, int& iphi);
  void smearTower(const TLorentzVector& jet, double e, bool calcSmearFactor, float& te, float& tem, float& thad, 
		  float& tout, float& temtrue, float& thadtrue, float& touttrue);  
  int  splitJet(const TLorentzVector& jet ,float* et,float* eta,float * phi, int* ieta,int* iphi);
  void calculateSmearFactor(const TLorentzVector& jet, double pt);


 public:
  ToyMC();
  ~ToyMC() {
    delete random_;
    if( histResp_ ) delete histResp_;
    if( histPtEta_ ) delete histPtEta_;
  }
  int etaBin(float eta) const;
  float etaBinCenter(int etaBin) const { return 0.5*(etaLowerBinEdge(etaBin)+etaUpperBinEdge(etaBin)); }
  float etaUpperBinEdge(int etaBin) const { return etaBinEdge(etaBin, false); }
  float etaLowerBinEdge(int etaBin) const { return etaBinEdge(etaBin, true ); }
  float etaBinEdge(int etaBin, bool lowerEdge) const;
  int generatePhotonJetTree(TTree *tree, int nevents);
  int generateTrackClusterTree(TTree *tree, int nevents);
  int generateDiJetTree(TTree* CalibTree, int nevents);
  int generateTopTree(TTree* CalibTree, int nevents);
  int makeTrackCluster(const char* filename, int nevents);
  int makePhotonJet(const char* filename, int nevents);
  int makeDiJet(const char* filename, int nevents);
  int makeTop(const char* filename, int nevents);
  void init(const std::string& configfile);
  void print() const;
};

#endif

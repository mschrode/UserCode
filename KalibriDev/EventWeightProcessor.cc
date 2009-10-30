#include "EventWeightProcessor.h"


//!  \param configfile Name of configuration file
//!  \param param The parameters
// -----------------------------------------------------------------
EventWeightProcessor::EventWeightProcessor(const std::string& configfile, TParameters* param)
  : EventProcessor(configfile,param) {
  
  ConfigFile config(configfile.c_str());
  weightEvents_ = config.read<bool>("weight events",false);
  if( weightEvents_ ) {
    std::string weightConfigFileName  = config.read<string>("event weights config","");

    std::cout << "Weighting events:\n";
    std::cout << "  Config file " << weightConfigFileName << "\n";

    ConfigFile weightConfig(weightConfigFileName.c_str());
    minPtHat_                    = bag_of<double>(weightConfig.read<string>("Min pthat","0."));
    std::vector<double> xSection = bag_of<double>(weightConfig.read<string>("XS","0."));
    std::vector<int> nEvents     = bag_of<int>(weightConfig.read<string>("Number evts","0"));
    double lumi                  = weightConfig.read<double>("Lumi",0.);
    int refPtHatBin              = weightConfig.read<int>("Reference pthat bin",-1);

    assert( xSection.size() == minPtHat_.size() );
    assert( nEvents.size() == minPtHat_.size() );
    calculateWeights(xSection,nEvents,lumi,refPtHatBin);
    assert( weights_.size() == minPtHat_.size() );

    std::cout << "  Calculating weights ";
    if( refPtHatBin < 0 ) std::cout << "for a luminosity of " << lumi << " pb-1:\n";
    else                  std::cout << "relative to pthat bin " << refPtHatBin << ":\n";
    for(size_t i = 0; i < minPtHat_.size(); i++) {
      std::cout << "    " << i << ": " << minPtHat_.at(i) << " < pthat < ";
      if( i < minPtHat_.size() - 1 ) std::cout << minPtHat_.at(i+1);
      else                           std::cout << " infty";
      std::cout << ":  " << weights_.at(i) << std::endl;
    }
  }
}



// -----------------------------------------------------------------
EventWeightProcessor::~EventWeightProcessor() {
  minPtHat_.clear();
  weights_.clear();
}



//!  \brief Apply weights to the events
//!
//!  Find the \f$ \hat{p}_{T} \f$ bin of each event and apply
//!  the corresponding weight, if \p weightEvents_ \p
//!  is true.
//!
//!  \param data The data
//!  \return Number of weighted events
// -----------------------------------------------------------------
int EventWeightProcessor::process(std::vector<TData*>& data) {
  int nProcEvts = 0; // Number of processed events

  if( weightEvents_ ) {
    std::vector<TData*>::iterator evt = data.begin();
    for(; evt != data.end(); evt++) {
      if( (*evt)->GetType() == ParLimit ) continue;

      double weight   = 0.;
      for(size_t ptHatBin = 0; ptHatBin < minPtHat_.size(); ptHatBin++) {
	if( (*evt)->ptHat() > minPtHat_.at(ptHatBin) )
	  weight = weights_.at(ptHatBin);
        else
	  break;
      }
      (*evt)->setWeight( weight );

      nProcEvts++;
    }

    std::cout << "  Applied weights for " << nProcEvts << " events\n";
  }

  return nProcEvts;
}



//!  \brief Calculate weights
//!
//!  A weighting factor \f$ w \f$ is calculated per
//!  \f$ \hat{p}_{T} \f$ bin \f$ i \f$ as
//!  \f[
//!   w = N_{ref} \frac{\sigma_{i}}{N_{i}}
//!  \f]
//!  where \f$ \sigma_{i} \f$ is the cross section and
//!  \f$ N_{i} \f$ the number of events in that bin.
//!  \f$ N_{ref} \f$ is the reference to which the weight
//!  is normalized. This is either, if \p refPtHatBin \p
//!  is negative, the luminosity (\f$ N_{ref} = \mathcal{L} \f$)
//!  or else the number of events in bin \p k = refPtHatBin \p
//!  (\f$ N_{ref} = N_{k} / \sigma_{k}  \f$).
//!
//!  \note \p xSection \p and \p nEvents \p must have the same size.
//!
//!  \param xSection Cross sections \f$ \sigma_{i} \f$ for the different
//!                  \f$ \hat{p}_{T} \f$ bins \f$ i \f$
//!  \param nEvents Number of events \f$ N_{i} \f$ in the different
//!                 \f$ \hat{p}_{T} \f$ bins \f$ i \f$
//!  \param lumi Luminosity \f$ \mathcal{L} \f$
//!  \param refPtHatBin If \f$ \ge 0 \f$, the \f$ \hat{p}_{T} \f$
//!                     bin to which the events are weighted
// -----------------------------------------------------------------
void EventWeightProcessor::calculateWeights(const std::vector<double>& xSection,
					    const std::vector<int>& nEvents,
					    double lumi, int refPtHatBin) {
  // Set up new empty vector of weights
  weights_.clear();
  weights_ = std::vector<double>(xSection.size());

  // The reference event number
  double refNum = lumi;      // Weight events relative to luminosity
  if( refPtHatBin >= 0 ) {   // or weight events relative to reference bin
    refNum = nEvents.at(refPtHatBin) / xSection.at(refPtHatBin);
  }

  // Calculate weights
  for( size_t i = 0; i < xSection.size(); i++ ) {
    weights_.at(i) = refNum * xSection.at(i) / nEvents.at(i);
  }
}


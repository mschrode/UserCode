#ifndef DIJETREADER_H
#define DIJETREADER_H
//!
//!  \brief Reader for dijet events
//!
//!  This class reads dijet events from a ROOT-tree. They are read
//!  by calling the method readEvents(std::vector<Event*>& data).
//!  The data is stored in a format derived from Event depending on the
//!  'Di-Jet data class' field in the config file:
//!    - class 0: Event_PtBalance
//!
//!  It is also possible to store the first two jets of the dijet
//!  event as a JetTruthEvent, where the \f$ p^{\textrm{gen}}_{T} \f$
//!  is used as truth. This is useful for Monte Carlo based calibration:
//!    - class 11: JetTruthEvent of Jet
//!    - class 12: JetTruthEvent of JetWithTowers
//!  In this case, the following cuts as defined in the config file
//!  are applied on each jet in this order:
//!    -# Number of jets greater than 2
//!    -# \f$ p^{\textrm{gen}}_{T} > \texttt{Et genJet min both Jets} \f$
//!    -# \f$ p^{\textrm{gen}}_{T} < \texttt{Et genJet max both Jets} \f$
//!    -# \f$ \Delta R < \texttt{DeltaR cut on jet matching} \f$
//!    -# \f$ p^{\textrm{jet}}_{T} < \texttt{Et cut on jet} \f$
//!    -# \f$ |\eta| > \texttt{Eta cut on jet} \f$
//!    -# \f$ \textrm{Hadronic fraction} < \texttt{Max had fraction} \f$
//!    -# \f$ \textrm{Hadronic fraction} > \texttt{Min had fraction} \f$
//!  
//!
//!  \author Hartmut Stadie
//!  \date 2008/12/12
//!  $Id: DiJetReader.h,v 1.2 2010/01/21 16:48:47 mschrode Exp $
// ----------------------------------------------------------------   



#include "EventReader.h"

#include <string>
#include <memory>


class Jet;
class NJetSel;
class TRandom;

class DiJetReader : public EventReader{
 public:
  DiJetReader(const std::string& configfile, TParameters *p);
  ~DiJetReader();
  int readEvents(std::vector<Event*>& data);


 private:
  Event* createTwoJetsPtBalanceEvent();
  Event* createSmearEvent();
  int createJetTruthEvents(std::vector<Event*>& data);
  CorFactors* createCorFactors(int jetid) const;
  std::vector<Jet*> readCaloJets(int nJets) const;
  std::vector<Jet*> readGenJetSortedJets(int nJets) const;


  std::auto_ptr<NJetSel> nJet_;                //!< Njet Selector
  TRandom* rand_;             //!< Random number generator

  int    dataClass_;            //!< Data class, see also Event
  int    nDijetEvents_;         //!< Maximum number of read dijet events
  int    prescale_;             //!< only read every nth event

  double minJetEt_;             //!< Minimum pt of jet
  double minDijetEt_;           //!< Minimum dijet pt
  double maxDijetEt_;           //!< Maximum dijet pt
  double max3rdJetEt_;          //!< Maximum pt of 3rd jet in dijet event
  double maxRel3rdJetEt_;       //!< Maximum relative pt of 3rd jet in dijet event
  double minDeltaPhi_;          //!< Minimum DeltaPhi for 0 < DeltaPhi < Pi
  double maxJetEta_;            //!< Maximum absolute jet eta
  double minJetHadFraction_;    //!< Minimum jet Had/(Had+EMF)
  double maxJetHadFraction_;    //!< Maximum jet Had/(Had+EMF)
  double minGenJetEt_;          //!< Minimum pt of genJets of dijets
  double maxGenJetEt_;          //!< Maximum pt of genJets of dijets
  double maxDeltaR_;            //!< Maximum DeltaR

  int    nDiJetCut_;            //!< Number of events with less than 2 jets
  int    nMinJetEt_;            //!< Number of events rejected by minJetEt_ cut
  int    nMinDijetEt_;          //!< Number of events rejected by minDijetEt_ cut
  int    nMaxDijetEt_;          //!< Number of events rejected by maxDijetEt_ cut
  int    nCutOn3rdJet_;         //!< Number of events rejected by max3rdJetEt_ or maxRelJetEt_ cut
  int    nMinDeltaPhi_;         //!< Number of events rejected by maxDeltaPhi_ cut
  int    nMaxJetEta_;           //!< Number of events rejected by maxJetEta_ cut
  int    nMinJetHadFraction_;   //!< Number of events rejected by minJetHadFraction_ cut
  int    nMaxJetHadFraction_;   //!< Number of events rejected by maxJetHadFraction_ cut
  int    nMaxGenJetEt_;         //!< Number of events rejected by maxGenJetEt_ cut
  int    nMinGenJetEt_;         //!< Number of events rejected by minGenJetEt_ cut
  int    nMaxDeltaR_;           //!< Number of events rejected by maxDeltaR_ cut

  int    maxNIter_;             //!< Max number of iterations in integration
  double eps_;                  //!< Integration precision for convergence
  double min_;                  //!< Minimum of truth spectrum in integration
  double max_;                  //!< Maximum of truth spectrum in integration
  double truthSpecExp_;         //!< Exponent of truth spectrum
};


#endif

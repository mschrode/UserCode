#ifndef PHOTONJETREADER_H
#define PHOTONJETREADER_H
//!
//!  \brief Reader for Photon Jet Events
//!
//!  This class reads events according to the GammaJetSel.
//!
//!  \author Hartmut Stadie
//!  \date 2008/12/12
//!  $Id: PhotonJetReader.h,v 1.2 2010/01/21 16:49:21 mschrode Exp $
// ----------------------------------------------------------------   


#include "EventReader.h"

#include <string>
#include <memory>

class CorFactors;
class GammaJetSel;

class PhotonJetReader : public EventReader{
 public:
  PhotonJetReader(const std::string& configfile, TParameters *p);
  ~PhotonJetReader();
  int readEvents(std::vector<Event*>& data);

 private:
  Event* createJetTruthEvent();
  Event* createSmearEvent();
  CorFactors* createCorFactors(int jetid) const;

  std::auto_ptr<GammaJetSel> gammaJet_;        //!< Gamma-jet Selector

  int    dataClass_;            //!< Data class, see also Event
  int    nGammaJetEvents_;    //!< Maximum number of read photon jet events

  double minJetEt_;            //!< Minimum pt of jet
  double maxJetEt_;            //!< Maximum pt of jet
  double minGammaEt_;          //!< Minimum pt of photon
  double maxGammaEt_;          //!< Maximum pt of photon
  double maxRel2ndJetEt_;      //!< Maximum relative pt of non-leading jet
  double minGenJetEt_;         //!< Minimum pt of genJet
  double maxGenJetEt_;         //!< Maximum pt of genJet
  double maxDeltaR_;           //!< Maximum DeltaR
  double maxJetEta_;           //!< Maximum absolute jet eta
  double minJetHadFraction_;   //!< Minimum jet Had/(Had+EMF)
  double maxJetHadFraction_;   //!< Maximum jet Had/(Had+EMF)

  int    nMinJetEt_;           //!< Number of events rejected by \p minJetEt_ cut
  int    nMaxJetEt_;           //!< Number of events rejected by \p maxJetEt_ cut
  int    nMinGammaEt_;         //!< Number of events rejected by \p minGammaEt_ cut
  int    nMaxGammaEt_;         //!< Number of events rejected by \p maxGammaEt_ cut
  int    nMaxRel2ndJetEt_;     //!< Number of events rejected by \p maxRel2ndJetEt_ cut
  int    nMinGenJetEt_;        //!< Number of events rejected by \p minGenJetEt_ cut
  int    nMaxGenJetEt_;        //!< Number of events rejected by \p maxGenJetEt_ cut
  int    nMaxDeltaR_;          //!< Number of events rejected by \p maxDeltaR_ cut
  int    nMaxJetEta_;          //!< Number of events rejected by \p maxJetEta_ cut
  int    nMinJetHadFraction_;  //!< Number of events rejected by \p minJetHadFraction_ cut
  int    nMaxJetHadFraction_;  //!< Number of events rejected by \p maxJetHadFraction_ cut
};


#endif

#ifndef DASTreeMaker_h
#define DASTreeMaker_h

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TString.h>


class DASTreeMaker : public edm::EDAnalyzer {
 public:
  explicit DASTreeMaker(const edm::ParameterSet&);
  ~DASTreeMaker();


 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  void GetMCObjects(const reco::GenJetCollection& GenJets, const reco::GenMETCollection& GenMet, const reco::GenParticleCollection& Gen);
  void GetWdecay(const reco::GenParticleCollection& Gen);
  void GetZdecay(const reco::GenParticleCollection& Gen);
  void GetTTdecay(const reco::GenParticleCollection& Gen);
  void GetGenPhoton(const reco::GenParticleCollection& Gen);
  void GetRecoObjects(const edm::View<reco::Candidate>& jets, const pat::METCollection& patMet, const edm::View<reco::Candidate>& muons, const edm::View<reco::Candidate>& eles, const pat::PhotonCollection& patPhotons, const reco::VertexCollection& Vtx);
  void GetSUSYs(const LHEEventProduct& lhep, const GenEventInfoProduct& genProd);
  int GetProcID(int procID);
  int hadronicTauFlag(const reco::Candidate &cand) const;

  const int maxColSize_;  // Maximum number of elements of a collection stored in the ntuple

  double pfEventRho_;
  float evtWgt_;
  int sampleID_;
  bool isMCdata_, isSUSY_;

  edm::InputTag genJetsTag_, genMetsTag_, vertexTag_, jetsTag_, patMetsTag_, muonsTag_, elesTag_, patPhotonsTag_, pfRhoTag_, evtWgtTag_;

  std::string outFileName_;

  TFile* outFile_;
  TTree* dasTree_;
      
  int nvtx_, njet_, nmu_, nele_, nphot_, ngjet_, runnr_, evtnr_, lumib_, flgW_, flgZ_, flgTT_, flgSUSY_, flgTauHad_;
  float M0, M12, A0, tanB, sgnMu;
  float bosM, bosID, bosPx, bosPy, bosPz, bosE, bosQ;
  float *lepM, *lepID, *lepPx, *lepPy, *lepPz, *lepE, *lepQ;
  float gphoM,gphoID,gphoPx,gphoPy,gphoPz,gphoE;
  float genmet, genmetphi, recomet, recometphi;
      
  float *vtxz;
  float *jpx, *jpy, *jpz, *jen, *jpt, *jphi, *jeta;
  float *phpx, *phpy, *phpz, *phen, *phpt, *phphi, *pheta;
  float *muq, *mupx, *mupy, *mupz, *muen, *mupt, *muphi, *mueta;
  float *eleq, *elepx, *elepy, *elepz, *eleen, *elept, *elephi, *eleeta;

  float *gjpx, *gjpy, *gjpz, *gjen, *gjpt, *gjphi, *gjeta;
};

#endif


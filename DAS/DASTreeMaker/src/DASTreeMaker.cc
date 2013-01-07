// system include files
#include <memory>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <map>
#include "DAS/DASTreeMaker/interface/DASTreeMaker.h"



DASTreeMaker::DASTreeMaker(const edm::ParameterSet& conf)
{
  isMCdata_    = conf.getParameter< bool > ("MCdata");
  isSUSY_      = conf.getParameter< bool > ("isSUSY");
  sampleID_    = conf.getParameter< int > ("sampleID");
  evtWgt_      = conf.getParameter< double > ("evtWgt");
  evtWgtTag_   = conf.getParameter< edm::InputTag > ("evtWgtTag");
  genJetsTag_  = conf.getParameter< edm::InputTag > ("genjets");
  genMetsTag_  = conf.getParameter< edm::InputTag > ("genmet");
  vertexTag_   = conf.getParameter< edm::InputTag > ("vertex");
  patJetsTag_  = conf.getParameter< edm::InputTag > ("PATjets");
  patMetsTag_  = conf.getParameter< edm::InputTag > ("PATmet");
  patMuonsTag_ = conf.getParameter< edm::InputTag > ("PATmuons");
  muID_        = conf.getParameter<std::string>("muID" );
  patElesTag_  = conf.getParameter< edm::InputTag > ("PATelectrons");
  eleID_       =  conf.getParameter<std::string>("eleID" );
  patPhotonsTag_ = conf.getParameter< edm::InputTag > ("PATphotons");
  pfRhoTag_    = conf.getParameter<edm::InputTag>("PFRhoTag");
  outFileName_       =  conf.getParameter<std::string>("OutFile" );
}


DASTreeMaker::~DASTreeMaker()
{ 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

// ------------ method called to for each event  ------------
void DASTreeMaker::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  // Event provenance information
  runnr_ = evt.id().run();
  evtnr_ = evt.id().event();
  lumib_ = evt.luminosityBlock();

  // Event weight
  if (evtWgt_<0) {
    edm::Handle<double> evtWgtHandle;
    evt.getByLabel(evtWgtTag_,evtWgtHandle);
    evtWgt_ = *evtWgtHandle;
  }

  //---get generator level objects
  if (isMCdata_) {
    edm::Handle<reco::GenJetCollection> genJetsHandle;
    edm::Handle<reco::GenMETCollection> genMetsHandle;
    edm::Handle<reco::GenParticleCollection> genPartsHandle;

    evt.getByLabel (genJetsTag_,genJetsHandle);
    evt.getByLabel (genMetsTag_,genMetsHandle);
    evt.getByLabel ("genParticles",genPartsHandle);
    
    if (sampleID_==24) GetWdecay(*genPartsHandle);
    if (sampleID_==23) GetZdecay(*genPartsHandle);
    if (sampleID_==6) GetTTdecay(*genPartsHandle);
    if (sampleID_==22) GetGenPhoton(*genPartsHandle);

    GetMCObjects(*genJetsHandle, *genMetsHandle, *genPartsHandle);
    //---get SUSY parameters
    if (isSUSY_) {
      edm::Handle<GenEventInfoProduct> genInfProdHandle;
      edm::Handle<LHEEventProduct> lheEvtProdHandle;

      evt.getByLabel("source", lheEvtProdHandle);
      evt.getByLabel("generator",genInfProdHandle);
    
      if (sampleID_==0) GetSUSYs(*lheEvtProdHandle,*genInfProdHandle);
      if (sampleID_==1) flgSUSY_=GetProcID((*genInfProdHandle).signalProcessID());
    }
  }
  
  //---get reco level objects
  edm::Handle<reco::VertexCollection> VtxHandle;
  edm::Handle<pat::JetCollection> patJetsHandle;
  edm::Handle<pat::METCollection> patMETHandle;
  edm::Handle<pat::MuonCollection> patMuonsHandle;
  edm::Handle<pat::ElectronCollection> patElectronsHandle;
  edm::Handle<pat::PhotonCollection> patPhotonsHandle;
  edm::Handle<double> pfRhoHandle;

  evt.getByLabel (vertexTag_,VtxHandle);
  evt.getByLabel (patJetsTag_,patJetsHandle);
  evt.getByLabel (patMetsTag_,patMETHandle);
  evt.getByLabel (patMuonsTag_,patMuonsHandle);
  evt.getByLabel (patElesTag_,patElectronsHandle);
  evt.getByLabel (patPhotonsTag_,patPhotonsHandle);

  // get pfRho for pileup correction to isolation
  evt.getByLabel(pfRhoTag_, pfRhoHandle);
  pfEventRho_ = *pfRhoHandle;

  GetRecoObjects(*patJetsHandle, *patMETHandle, *patMuonsHandle, *patElectronsHandle, *patPhotonsHandle, *VtxHandle);
  
  dasTree_->Fill();
}

void DASTreeMaker::GetMCObjects(const reco::GenJetCollection& GenJets, const reco::GenMETCollection& GenMet, const reco::GenParticleCollection& Gen) 
{

  int jgen(0);

  if (&GenMet){
    for (reco::GenMETCollection::const_iterator gmt=GenMet.begin(); gmt!=GenMet.end();gmt++){
      genmet = gmt->pt();
      genmetphi = gmt->phi();
    }
  }
  
  if (&GenJets){
    for (reco::GenJetCollection::const_iterator gjt=GenJets.begin(); gjt!=GenJets.end(); gjt++){
      gjpx[jgen] = gjt->px();
      gjpy[jgen] = gjt->py();
      gjpz[jgen] = gjt->pz();
      gjen[jgen] = gjt->energy();
      gjpt[jgen] = gjt->pt();
      gjphi[jgen] = gjt->phi();
      gjeta[jgen] = gjt->eta();
      jgen++;
    }
    //printf("Njgen=%i\n",jgen);
  }
  ngjet_ = jgen;
  
  
}

void DASTreeMaker::GetWdecay(const reco::GenParticleCollection& Gen){
  int dcy_lep = 0;
  flgW_=0;
  flgTauHad_=0;
  //----flag W+jets events according to the W decays
  if (&Gen) {
    for (reco::GenParticleCollection::const_iterator igp = Gen.begin(); igp!=Gen.end();igp++) {
      const reco::GenParticle& genP = *igp;
      if (genP.status()!=3) continue;
      if (std::abs(genP.pdgId())>25) continue;
      //cout <<genP.pdgId()<<" "<<genP.status()<<" "<<genP.px()<<" "<<genP.py()<<" "<<genP.pz()<<" "<<genP.energy()<<" "<<genP.mass()<<endl;
      //---pick-up the W boson 4-vector
      if (std::abs(genP.pdgId())==24){
        //cout <<genP.pdgId()<<" "<<genP.status()<<" "<<genP.px()<<" "<<genP.py()<<" "<<genP.pz()<<" "<<genP.energy()<<" "<<genP.mass()<<endl;
        bosM = genP.mass();
        bosID = genP.pdgId();
        bosPx = genP.px();
        bosPy = genP.py();
        bosPz = genP.pz();
        bosE = genP.energy();
        bosQ = genP.charge();
      }
      
      //----pick-up the leptonic W decays 4-vectors
      if (std::abs(genP.pdgId())<17 && std::abs(genP.pdgId())>10) {
	if (std::abs(genP.mother()->pdgId())==24) {	
	  if (dcy_lep<2) {
	    lepM[dcy_lep] = genP.mass();
	    lepID[dcy_lep] = genP.pdgId();
	    lepPx[dcy_lep] = genP.px();
	    lepPy[dcy_lep] = genP.py();
	    lepPz[dcy_lep] = genP.pz();
	    lepE[dcy_lep] = genP.energy();
	    lepQ[dcy_lep] = genP.charge(); 
	    //---flag is equal to the lepton ID
	    if (std::abs(genP.pdgId())%2==1) flgW_ = std::abs(genP.pdgId());
	    //---flag the tau decay
	    if( std::abs(genP.pdgId()) == 15 ) flgTauHad_ = hadronicTauFlag(genP);
	  }
	  dcy_lep++;
	}
      }
    }
  }
  //cout <<"W-decay ID="<<flgW_<<endl;
}
void DASTreeMaker::GetZdecay(const reco::GenParticleCollection& Gen){
  //----flag Z+jets events according to the Z decays
  flgZ_=0;
  int dcy_lep = 0;
  if (&Gen) {
    for (reco::GenParticleCollection::const_iterator igp = Gen.begin(); igp!=Gen.end();igp++) {
      const reco::GenParticle& genP = *igp;

      if (genP.status()!=3) continue;
      if (std::abs(genP.pdgId())>25) continue;
      //---pick-up the Z boson 4-vector
      if (genP.pdgId()==23){
        //cout <<genP.pdgId()<<" "<<genP.status()<<" "<<genP.px()<<" "<<genP.py()<<" "<<genP.pz()<<" "<<genP.energy()<<" "<<genP.mass()<<endl;
        bosM = genP.mass();
        bosID = genP.pdgId();
        bosPx = genP.px();
        bosPy = genP.py();
        bosPz = genP.pz();
        bosE = genP.energy();
        bosQ = genP.charge();
      }
      //----pick-up the WZ decays 4-vectors
      if (std::abs(genP.pdgId())<17 && std::abs(genP.pdgId())>10) {
	if (std::abs(genP.mother()->pdgId())==23) {
	  if (dcy_lep<2) {
	    lepM[dcy_lep] = genP.mass();
	    lepID[dcy_lep] = genP.pdgId();
	    lepPx[dcy_lep] = genP.px();
	    lepPy[dcy_lep] = genP.py();
	    lepPz[dcy_lep] = genP.pz();
	    lepE[dcy_lep] = genP.energy();
	    lepQ[dcy_lep] = genP.charge(); 
	  }
	  dcy_lep++;
	  //---flag is equal to the lepton ID
	  flgZ_ = std::abs(genP.pdgId());
	}
      }
    }
  }
}
void DASTreeMaker::GetGenPhoton(const reco::GenParticleCollection& Gen){
  if (&Gen) {
    for (reco::GenParticleCollection::const_iterator igp = Gen.begin(); igp!=Gen.end();igp++) {
      const reco::GenParticle& genP = *igp;

      if (genP.status()!=3) continue;
      if (std::abs(genP.pdgId())>22) continue;
      //---pick-up the photon 4-vector
      if (genP.pdgId()==22){
        //cout <<genP.pdgId()<<" "<<genP.status()<<" "<<genP.px()<<" "<<genP.py()<<" "<<genP.pz()<<" "<<genP.energy()<<" "<<genP.mass()<<endl;
        gphoM  = genP.mass();
        gphoID = genP.pdgId();
        gphoPx = genP.px();
        gphoPy = genP.py();
        gphoPz = genP.pz();
        gphoE  = genP.energy();
      }
    }
  }
}
void DASTreeMaker::GetTTdecay(const reco::GenParticleCollection& Gen){
  //----flag TT+jets events according to the TT decays
  int dcy_lep(0), qcnt(0), taucnt(0), lepcnt(0);
  flgTT_=0;
  if (&Gen) {
    for (reco::GenParticleCollection::const_iterator igp = Gen.begin(); igp!=Gen.end();igp++) {
      const reco::GenParticle& genP = *igp;

      if (genP.status()!=3) continue;
      if (std::abs(genP.pdgId())>25) continue;
      if (std::abs(genP.pdgId())<17 && std::abs(genP.pdgId())>0) {
	if (std::abs(genP.mother()->pdgId())==6 || std::abs(genP.mother()->pdgId())==24) {
	  if (dcy_lep<6) {
	    lepM[dcy_lep] = genP.mass();
	    lepID[dcy_lep] = genP.pdgId();
	    lepPx[dcy_lep] = genP.px();
	    lepPy[dcy_lep] = genP.py();
	    lepPz[dcy_lep] = genP.pz();
	    lepE[dcy_lep] = genP.energy();
	    lepQ[dcy_lep] = genP.charge(); 
	  }
	  dcy_lep++;
	}
      }
    }
  }
  //---label the decays: 0-allhad, 1-tau+had, 2-lep+had, 3-ditau, 4-dilep, 5-tau+lep
  for (int ip=0;ip<6;ip++) {
    if (std::abs(lepID[ip])<5) qcnt++;
    if (std::abs(lepID[ip])==15) taucnt++;
    if (std::abs(lepID[ip])==11 || std::abs(lepID[ip])==13) lepcnt++;
  }
  if (qcnt==4) flgTT_=0;
  if (qcnt==2 && taucnt==1) flgTT_=1;
  if (qcnt==2 && lepcnt==1) flgTT_=2;
  if (taucnt==2) flgTT_=3;
  if (lepcnt==2) flgTT_=4;
  if (lepcnt==1 && taucnt==1) flgTT_=5;
}


// Returns 1 for fully-hadronic decay of tau
// 0 otherwise
int DASTreeMaker::hadronicTauFlag(const reco::Candidate &cand) const {
  int flag = 1;			// Assume hadronically decaying tau
  for(reco::Candidate::const_iterator itc = cand.begin();
      itc != cand.end(); ++itc) {
    int pdgIdDaughter = std::abs(itc->pdgId());
    // tau might not be the finally decaying tau but have
    // an intermediate tau or w in the decay
    if( pdgIdDaughter == 15 || pdgIdDaughter == 24) {
      flag = hadronicTauFlag(*itc);
      if( flag == 0 ) break;      
    } else if( pdgIdDaughter == 11 || pdgIdDaughter == 13 ) { // Leptonic decay
      flag = 0;
      break;
    }
  }

  return flag;
}



void DASTreeMaker::GetRecoObjects(const pat::JetCollection& patJets, const pat::METCollection& patMet, const pat::MuonCollection& patMuons, const pat::ElectronCollection& patEles, const pat::PhotonCollection& patPhotons, const reco::VertexCollection& Vtx) {

  int jcal(0), mucnt(0), elecnt(0), phcnt(0), vtxcnt(0);
  float mudz(10), mudxy(10);
  float eledz(10), eledxy(10);
  math::XYZPoint primVtx;
  double trackerIso, ecalIso, hcalIso, sigEtaEta;

  //---pat met
  if (&patMet){
    for (pat::METCollection::const_iterator mt=patMet.begin(); mt!=patMet.end();mt++){
      //---corr MET (JES+mu)
      recomet = mt->pt();
      recometphi = mt->phi();
    }
  }
  
  //--vertex
  if (&Vtx){
    for (reco::VertexCollection::const_iterator itv=Vtx.begin(); itv!=Vtx.end(); itv++){
      if (vtxcnt>=20) break;
      if (vtxcnt==0) primVtx = itv->position();
      vtxz[vtxcnt] = itv->z();     
      vtxcnt++;
    }
  }
  nvtx_=vtxcnt;

  //----pat muons
  if (&patMuons){
    //cout <<"Nmus:"<<patMuons.size()<<endl;
    //nmu_ = patMuons.size();
    for(pat::MuonCollection::const_iterator mu = patMuons.begin(); mu != patMuons.end(); mu++) {
      //cout << mu->pt()<<" "<<mu->isGood("GlobalMuonPromptTight")<<" "<<mu->dB()<<" "<<mu->numberOfValidHits()<<" "<<mu->normChi2()<<" "<<mu->isolationR03().sumPt<<" "<<mu->isolationR03().emEt<<" "<<mu->isolationR03().hadEt<<" "<<mu->isolationR03().trackerVetoPt<<" "<<mu->isolationR03().emVetoEt<<" "<<mu->isolationR03().hadVetoEt<<" "<<mu->trackIso()<<" "<<mu->ecalIso()<<" "<<mu->hcalIso()<<endl;
      //if (mu->isGood("GlobalMuonPromptTight")==0) continue;
      if (mu->isGood(muID_)==0) continue;
      if ((mu->chargedHadronIso()+mu->photonIso()+mu->neutralHadronIso())/(mu->pt()) >= 0.2) continue;
      if (mu->innerTrack()->numberOfValidHits() <11) continue;
      if (vtxcnt!=0) mudz = mu->innerTrack()->dz(primVtx);
      if (vtxcnt!=0) mudxy = mu->innerTrack()->dxy(primVtx);
      if (mudxy>=0.02) continue;
      if (mudz>=1) continue;
      
      muq[mucnt] = mu->charge();
      mupx[mucnt] = mu->px();
      mupy[mucnt] = mu->py();
      mupz[mucnt] = mu->pz();
      muen[mucnt] = mu->energy();
      mupt[mucnt] = mu->pt();
      muphi[mucnt] = mu->phi();
      mueta[mucnt] = mu->eta();
      mucnt++;
    }
    //printf("Nmu=%i\n",mucnt);
  }
  nmu_ = mucnt;
  
  //---pat electrons
  if (&patEles){
    for(pat::ElectronCollection::const_iterator ele = patEles.begin(); ele != patEles.end(); ele++) {
      //if (!isMCdata_ && ele->electronID("eidTight")==0) continue;
      //if (isMCdata_ && ele->electronID("eidTightMC")==0) continue;
      if (ele->electronID(eleID_)==0) continue;
      if (std::abs(ele->eta()) >= 2.5 ) continue;
      if (std::abs(ele->eta()) >= 1.4442 && std::abs(ele->eta()) < 1.566 ) continue;
      if ((ele->chargedHadronIso()+ele->photonIso()+ele->neutralHadronIso())/(ele->pt()) >= 0.2) continue;
      if (vtxcnt!=0) eledz = std::abs(ele->vz() - primVtx.Z());
      if (vtxcnt!=0) eledxy = ele->gsfTrack()->dxy(primVtx);
      if (eledxy>=0.04) continue;
      if (eledz>=1) continue;

      eleq[elecnt] = ele->charge();
      elepx[elecnt] = ele->px();
      elepy[elecnt] = ele->py();
      elepz[elecnt] = ele->pz();
      eleen[elecnt] = ele->energy();
      elept[elecnt] = ele->pt();
      elephi[elecnt] = ele->phi();
      eleeta[elecnt] = ele->eta();
      elecnt++;
    }
  }
  nele_ = elecnt;

  //--pat photons
  double areaR04       = 0.4 * 0.4 * 4.0*atan(1.0);
  if (&patPhotons) {
    for(pat::PhotonCollection::const_iterator pht = patPhotons.begin(); pht != patPhotons.end(); pht++) {
      if (pht->hasPixelSeed()) continue;
      if (pht->hadronicOverEm()>0.05) continue;
      sigEtaEta   = pht->sigmaIetaIeta();
      bool   showerShape = ( (std::abs(pht->eta())<1.379 && sigEtaEta<0.01) || (std::abs(pht->eta())>1.579 && sigEtaEta<0.028) );
      if (!showerShape) continue;
      trackerIso  = pht->trkSumPtHollowConeDR04();
      ecalIso     = pht->ecalRecHitSumEtConeDR04();
      hcalIso     = pht->hcalTowerSumEtConeDR04();
      bool   isPhotonIso = ( (trackerIso+ecalIso+hcalIso - areaR04*pfEventRho_)<5.0 );
      if (!isPhotonIso) continue;

      phpx[phcnt] = pht->px();
      phpy[phcnt] = pht->py();
      phpz[phcnt] = pht->pz();
      phen[phcnt] = pht->energy();
      phpt[phcnt] = pht->pt();
      phphi[phcnt] = pht->phi();
      pheta[phcnt] = pht->eta();
      phcnt++;
    }
  }
  nphot_ = phcnt;

  //---TO DO: work with the jets after leptons and photons so that they can be cleaned of them
  //---pat jets
  if (&patJets) {    
    for(pat::JetCollection::const_iterator jt = patJets.begin(); jt != patJets.end(); jt++) {
      //cout <<" pt:"<<jt->pt()<<" "<<jt->correctedJet("RAW").pt()<<" eta:"<<jt->eta()<<" "<<jt->correctedJet("RAW").eta()<<" phi:"<<jt->phi()<<" "<<jt->correctedJet("RAW").phi()<<" emf:"<<jt->emEnergyFraction()<<" "<<jt->correctedJet("RAW").emEnergyFraction()<<endl;
      jpx[jcal] = jt->px();
      jpy[jcal] = jt->py();
      jpz[jcal] = jt->pz();
      jen[jcal] = jt->energy();
      jpt[jcal] = jt->pt();
      jphi[jcal] = jt->phi();
      jeta[jcal] = jt->eta();
      jcal++;
    }
    //printf("Njcal=%i\n",jcal);
  }
  njet_ = jcal;


}
      
void DASTreeMaker::GetSUSYs(const LHEEventProduct& lhep, const GenEventInfoProduct& genProd) {
  typedef std::vector<std::string>::const_iterator comments_const_iterator;
  comments_const_iterator c_begin = lhep.comments_begin();
  comments_const_iterator c_end = lhep.comments_end();
  
  double mzero, mhalf, tanb, azero, mu=1.0;
  double signMu;
  for( comments_const_iterator cit=c_begin; cit!=c_end; ++cit) {
    size_t found = (*cit).find("model");
    if( found != std::string::npos)   {    
      //         std::cout << *cit << std::endl;  
      size_t foundLength = (*cit).size();
      found = (*cit).find("=");
      std::string smaller = (*cit).substr(found+1,foundLength);
      found = smaller.find("_");
      smaller = smaller.substr(found+1,smaller.size());
      //
      std::istringstream iss(smaller);
      iss >> mzero;
      iss.clear();
      //
      found = smaller.find("_");
      smaller = smaller.substr(found+1,smaller.size());
      iss.str(smaller);
      iss >> mhalf;
      iss.clear();
      //
      found = smaller.find("_");
      smaller = smaller.substr(found+1,smaller.size());
      iss.str(smaller);
      iss >> tanb;
      iss.clear();
      //
      found = smaller.find("_");
      smaller = smaller.substr(found+1,smaller.size());
      iss.str(smaller);
      iss >> azero;
      iss.clear();
      found = smaller.find("_");
      smaller = smaller.substr(found+1,smaller.size());
      iss.str(smaller);
      iss >> signMu;
      iss.clear();
      mu *= signMu;
    }
  }
  //char buffer[100];
  //int n =sprintf(buffer,"mSugra model with parameters m0=%6.2f m12=%6.2f tanb=%6.2f A0=%6.2f mu=%6.2f\n",mzero,mhalf,tanb,azero,mu);
  //std::cout << buffer ;
  
  M0 = (float)mzero;
  M12 = (float)mhalf;
  A0 = (float)azero;
  tanB = (float)tanb;
  sgnMu = (float)mu;

}
int DASTreeMaker::GetProcID(int procID){
  
  int susy_proc_labelID = 0;

  if (procID<=214 && procID>=201) susy_proc_labelID = 4;//---ll
  if (procID<=236 && procID>=216) susy_proc_labelID = 3;//---nn
  if (procID<=242 && procID>=237) susy_proc_labelID = 1;//---ng
  if (procID<=244 && procID>=243) susy_proc_labelID = 9;//---gg
  if (procID<=256 && procID>=246) susy_proc_labelID = 2;//---ns
  if (procID<=259 && procID>=258) susy_proc_labelID = 10;//---sg
  if (procID<=265 && procID>=261) susy_proc_labelID = 11;//---tb
  if (procID<=273 && procID>=271) susy_proc_labelID = 6;//---ss
  if (procID<=280 && procID>=274) susy_proc_labelID = 5;//---sb
  if (procID<=283 && procID>=281) susy_proc_labelID = 6;//---ss
  if (procID<=286 && procID>=284) susy_proc_labelID = 5;//---sb
  if (procID<=290 && procID>=287) susy_proc_labelID = 8;//---bb
  if (procID<=293 && procID>=291) susy_proc_labelID = 6;//---ss
  if (procID<=295 && procID>=294) susy_proc_labelID = 7;//---bg
  if (procID==296) susy_proc_labelID = 8;//---bb

  //---so far 12 possible labels 
  return susy_proc_labelID;
}

void DASTreeMaker::beginJob() {
  std::cout << " Beginning Analysis " << std::endl;
  
  outFile_ = new TFile(outFileName_.c_str(),"RECREATE");
  dasTree_ = new TTree("AnaTree","");
   
  //--gen jets
  gjpx = new float[200];
  gjpy = new float[200];
  gjpz = new float[200];
  gjen = new float[200];
  gjpt = new float[200];
  gjphi = new float[200];
  gjeta = new float[200];
  //---gen muons
  lepM = new float[6];
  lepID = new float[6];
  lepQ = new float[6];
  lepPx = new float[6];
  lepPy = new float[6];
  lepPz = new float[6];
  lepE = new float[6];

  //---gen photons
  gphpx = new float[200];
  gphpy = new float[200];
  gphpz = new float[200];
  gphen = new float[200];
  gphpt = new float[200];
  gphphi = new float[200];
  gpheta = new float[200];

  //--reco jets
  jpx = new float[200];
  jpy = new float[200];
  jpz = new float[200];
  jen = new float[200];  
  jpt = new float[200];
  jphi = new float[200];
  jeta = new float[200];
  //---muons
  muq = new float[200];
  mupx = new float[200];
  mupy = new float[200];
  mupz = new float[200];
  muen = new float[200];
  mupt = new float[200];
  muphi = new float[200];
  mueta = new float[200];
  //---electrons
  eleq = new float[200];
  elepx = new float[200];
  elepy = new float[200];
  elepz = new float[200];
  eleen = new float[200];
  elept = new float[200];
  elephi = new float[200];
  eleeta = new float[200];
  //---photons
  phpx = new float[200];
  phpy = new float[200];
  phpz = new float[200];
  phen = new float[200];
  phpt = new float[200];
  phphi = new float[200];
  pheta = new float[200];

  //---vertices
  vtxz = new float[50];

  //---default values
  nvtx_    = -1;
  njet_    = -1;
  nmu_     = -1;
  nele_    = -1;
  nphot_   = -1;
  ngjet_   = -1;
  runnr_   = -1;
  evtnr_   = -1;
  lumib_   = -1;
  flgW_    = -1;
  flgZ_    = -1;
  flgTT_   = -1;
  flgSUSY_ = -1;
  flgTauHad_ = 0;
  

  //--------Tree Branches---------
  
  //----event info
  dasTree_->Branch("RunNr",&runnr_,"RunNr/I");
  dasTree_->Branch("EvtNr",&evtnr_,"EvtNr/I");
  dasTree_->Branch("LumiB",&lumib_,"LumiB/I");
  dasTree_->Branch("EvtWgt",&evtWgt_,"EvtWgt/F");
  
  if (isMCdata_) {
    //---generator level jets
    dasTree_->Branch("NrecoJetGen",&ngjet_,"NrecoJetGen/I");
    dasTree_->Branch("recoJetGenPx",gjpx,"recoJetGenPx[NrecoJetGen]/F");
    dasTree_->Branch("recoJetGenPy",gjpy,"recoJetGenPy[NrecoJetGen]/F");
    dasTree_->Branch("recoJetGenPz",gjpz,"recoJetGenPz[NrecoJetGen]/F");
    dasTree_->Branch("recoJetGenE",gjen,"recoJetGenE[NrecoJetGen]/F");
    dasTree_->Branch("recoJetGenPt",gjpt,"recoJetGenPt[NrecoJetGen]/F");
    dasTree_->Branch("recoJetGenPhi",gjphi,"recoJetGenPhi[NrecoJetGen]/F");
    dasTree_->Branch("recoJetGenEta",gjeta,"recoJetGenEta[NrecoJetGen]/F");
    //---generator level MET
    dasTree_->Branch("recoMetGen",&genmet,"recoMetGen/F");
    dasTree_->Branch("recoMetGenPhi",&genmetphi,"recoMetGenPhi/F");
    //----gen level flags
    dasTree_->Branch("flgW",&flgW_,"flgW/I");
    dasTree_->Branch("flgZ",&flgZ_,"flgZ/I");
    dasTree_->Branch("flgTT",&flgTT_,"flgTT/I");
    dasTree_->Branch("flgSUSY",&flgSUSY_,"flgSUSY/I");
    dasTree_->Branch("flgTauHad",&flgTauHad_,"flgTauHad/I");
    //----generator level boson 
    dasTree_->Branch("bosM",&bosM,"bosM/F");
    dasTree_->Branch("bosID",&bosID,"bosID/F");
    dasTree_->Branch("bosPx",&bosPx,"bosPx/F");
    dasTree_->Branch("bosPy",&bosPy,"bosPy/F");
    dasTree_->Branch("bosPz",&bosPz,"bosPz/F");
    dasTree_->Branch("bosE",&bosE,"bosE/F");
    dasTree_->Branch("bosQ",&bosQ,"bosQ/F");
    //----generator level boson decays
    dasTree_->Branch("lepM",lepM,"lepM[6]/F");
    dasTree_->Branch("lepID",lepID,"lepID[6]/F");
    dasTree_->Branch("lepPx",lepPx,"lepPx[6]/F");
    dasTree_->Branch("lepPy",lepPy,"lepPy[6]/F");
    dasTree_->Branch("lepPz",lepPz,"lepPz[6]/F");
    dasTree_->Branch("lepE",lepE,"lepE[6]/F");
    dasTree_->Branch("lepQ",lepQ,"lepQ[6]/F");
    //----generator level photon
    dasTree_->Branch("gphoM", &gphoM, "gphoM/F");
    dasTree_->Branch("gphoID",&gphoID,"gphoID/F");
    dasTree_->Branch("gphoPx",&gphoPx,"gphoPx/F");
    dasTree_->Branch("gphoPy",&gphoPy,"gphoPy/F");
    dasTree_->Branch("gphoPz",&gphoPz,"gphoPz/F");
    dasTree_->Branch("gphoE", &gphoE, "gphoE/F");
  }
 
  //---reco level jets
  dasTree_->Branch("NrecoJet",&njet_,"NrecoJet/I");
  dasTree_->Branch("recoJetPx",jpx,"recoJetPx[NrecoJet]/F");
  dasTree_->Branch("recoJetPy",jpy,"recoJetPy[NrecoJet]/F");
  dasTree_->Branch("recoJetPz",jpz,"recoJetPz[NrecoJet]/F");
  dasTree_->Branch("recoJetE",jen,"recoJetE[NrecoJet]/F");
  dasTree_->Branch("recoJetPt",jpt,"recoJetPt[NrecoJet]/F");
  dasTree_->Branch("recoJetPhi",jphi,"recoJetPhi[NrecoJet]/F");
  dasTree_->Branch("recoJetEta",jeta,"recoJetEta[NrecoJet]/F");
  //---MET
  dasTree_->Branch("recoMetCal",&recomet,"recoMetCal/F");
  dasTree_->Branch("recoMetCalPhi",&recometphi,"recoMetCalPhi/F");
  
  //---vertex
  dasTree_->Branch("NVtx",&nvtx_,"NVtx/I");
  dasTree_->Branch("VtxZ",vtxz,"VtxZ[NVtx]/F");
  
  //----muons
  dasTree_->Branch("NrecoMu",&nmu_,"NrecoMu/I");    
  dasTree_->Branch("recoMuQ",muq,"recoMuQ[NrecoMu]/F");   
  dasTree_->Branch("recoMuPx",mupx,"recoMuPx[NrecoMu]/F");
  dasTree_->Branch("recoMuPy",mupy,"recoMuPy[NrecoMu]/F");
  dasTree_->Branch("recoMuPz",mupz,"recoMuPz[NrecoMu]/F");
  dasTree_->Branch("recoMuEn",muen,"recoMuEn[NrecoMu]/F");
  dasTree_->Branch("recoMuPt",mupt,"recoMuPt[NrecoMu]/F");
  dasTree_->Branch("recoMuPhi",muphi,"recoMuPhi[NrecoMu]/F");
  dasTree_->Branch("recoMuEta",mueta,"recoMuEta[NrecoMu]/F");
  //---electrons
  dasTree_->Branch("NrecoEle",&nele_,"NrecoEle/I");    
  dasTree_->Branch("recoEleQ",eleq,"recoEleQ[NrecoEle]/F");   
  dasTree_->Branch("recoElePx",elepx,"recoElePx[NrecoEle]/F");
  dasTree_->Branch("recoElePy",elepy,"recoElePy[NrecoEle]/F");
  dasTree_->Branch("recoElePz",elepz,"recoElePz[NrecoEle]/F");
  dasTree_->Branch("recoEleEn",eleen,"recoEleEn[NrecoEle]/F");
  dasTree_->Branch("recoElePt",elept,"recoElePt[NrecoEle]/F");
  dasTree_->Branch("recoElePhi",elephi,"recoElePhi[NrecoEle]/F");
  dasTree_->Branch("recoEleEta",eleeta,"recoEleEta[NrecoEle]/F");
  //---photons
  dasTree_->Branch("NrecoPho",&nphot_,"NrecoPho/I");    
  dasTree_->Branch("recoPhoPx",phpx,"recoPhoPx[NrecoPho]/F");
  dasTree_->Branch("recoPhoPy",phpy,"recoPhoPy[NrecoPho]/F");
  dasTree_->Branch("recoPhoPz",phpz,"recoPhoPz[NrecoPho]/F");
  dasTree_->Branch("recoPhoEn",phen,"recoPhoEn[NrecoPho]/F");
  dasTree_->Branch("recoPhoPt",phpt,"recoPhoPt[NrecoPho]/F");
  dasTree_->Branch("recoPhoPhi",phphi,"recoPhoPhi[NrecoPho]/F");
  dasTree_->Branch("recoPhoEta",pheta,"recoPhoEta[NrecoPho]/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void DASTreeMaker::endJob() {
  outFile_->cd();
  dasTree_->Write();
  delete dasTree_;
  dasTree_=0;

  if (outFile_!=0){
    outFile_->Write();
    delete outFile_;
    outFile_=0;
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(DASTreeMaker);


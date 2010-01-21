// $Id: ToyMC.cc,v 1.41 2010/01/08 18:14:36 mschrode Exp $

#include "ToyMC.h"

#include <cmath>
#include <iostream>
#include <map>
#include <cassert> 
#include <ext/hash_map>

#include "TF1.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"


#include "ConfigFile.h"


//!  \brief Default constructor
// -----------------------------------------------------------------
ToyMC::ToyMC() : type_(1), minEta_(-2.5),maxEta_(2.5),minPt_(30), maxPt_(400),ptSpectrum_(Uniform),
		 histPtEta_(0),chunks_(200),jetSpreadA_(0.5),jetSpreadB_(0),noOutOfCone_(true),
		 maxPi0Frac_(0.5),maxEmf_(0.5),responseModel_(Constant),histResp_(0),
		 resolutionModel_(Gauss), smearFactor_(1.), smearTowersIndividually_(true)
{
  pInput_ = new TLorentzVector();
  random_ = new TRandom3();
  random_->SetSeed(0);
}

ToyMC::~ToyMC() 
{
  delete random_;
  delete histResp_;
  delete histPtEta_;
  delete pInput_;
}

//!  \brief Generates the input (truth) of an event
//!
//!  Generates a random input lorentz vector of an event 
//!  (e.g. the lorentz vector of the photon for photon-jet
//!  events) which is then smeared due to detector response
//!  and  resolution effects etc.
//!  The randomly generated components are
//!  - pt between \p minPt_ and \p maxPt_ according to \p ptSpectrum_
//!  - eta uniformly between \p minEta_ and \p maxEta_
//!  - phi uniformly between 0 and 2Pi
//!  - a zero mass
// -----------------------------------------------------------------
void ToyMC::genInput() { 
  static double rand[3];
  random_->RndmArray(3,rand);
  double pt  = 0;
  double eta = 0;
  if(ptSpectrum_ == Uniform) {
    pt  = rand[0]*(maxPt_ - minPt_)+minPt_;
    eta = rand[1]*(maxEta_ - minEta_)+minEta_;
  } else if(ptSpectrum_ == PowerLaw) {
    pt  = minPt_ * pow(rand[0],-1.0/(parTruth_.at(0) - 1.0));
    eta = rand[1]*(maxEta_ - minEta_)+minEta_;
  } else if(ptSpectrum_ == PtEtaHistogram) {
    histPtEta_->GetRandom2(eta, pt);
  }
  pInput_->SetPtEtaPhiM(pt, eta, (rand[2]*2 * M_PI - M_PI), 0);
}


//!  \brief Transform eta and phi values of the center of
//!         the corresponding tower and find the tower indices
//!
//!  \param eta Eta, is transformed to eta of corresponding tower center
//!  \param phi Phi, is transformed to phi of corresponding tower center
//!  \param ieta Index of tower covering eta
//!  \param iphi Index of tower covering phi
// -----------------------------------------------------------------
void ToyMC::calIds(float& eta, float &phi, int& ieta, int& iphi) 
{
  ieta = etaBin(eta);
  if(useTowerCenterEtaPhi_) eta = etaBinCenter(ieta);

  assert(ieta != 0);
  assert(ieta <= 41);
  assert(ieta >= -41);

  const static float dPhi =  2 * M_PI / 72;
  if(phi < 0) phi += 2 * M_PI;
  iphi = (int)(phi / dPhi);
  iphi++;
  if(iphi > 72) iphi = 72; 
  if(useTowerCenterEtaPhi_) {
    phi = (iphi-0.5) * dPhi;
  }

  assert(phi < 2 * M_PI);
  assert(iphi > 0);
  assert(iphi <= 72);
}



//!  \brief Calculate pt smear factor due to
//!         response and resolution
//!  \param jet the four vector of the jet
//!  \param E  Energy for calculation of some smear factors
//----------------------------------------------------------
void ToyMC::calculateSmearFactor(const TLorentzVector* jet, double E) {
  // Reset smear factor
  smearFactor_ = 1.;

  // Pt
  double pt    = E * jet->Pt()/jet->E();

  // Apply response
  if( responseModel_ == Constant
      || (responseModel_ == Flat)
      || (responseModel_ == Exp)
      || (responseModel_ == Slope)) {
    smearFactor_ *= parResp_.at(0);
  }
  else if( responseModel_ == L3 ) {
    if( pt < 0.1 ) {
      smearFactor_ *= parResp_.at(0) - parResp_.at(1)/parResp_.at(3) +  parResp_.at(4);
    } else {
      smearFactor_ *= parResp_.at(0) - parResp_.at(1)/(pow(log10(pt),parResp_.at(2)) + parResp_.at(3)) + parResp_.at(4)/pt;
    }
  }
  else if( responseModel_ == SimpleInverse ) { 
    smearFactor_ *= 1. - parResp_.at(0)/(pt + parResp_.at(1));
  }
  else if( responseModel_ == StepEta ) {
    smearFactor_ *= jet->Eta() < 0 ? parResp_.at(0) : parResp_.at(1);
  }
  else if( responseModel_ == SinusEta ) {
    smearFactor_ *= 1. + parResp_.at(0)*sin( parResp_.at(1)*jet->Eta() );
  }
  else if( responseModel_ == SinusEtaSimpleInversePt ) {
    smearFactor_ *= 1. + parResp_.at(0)*sin( parResp_.at(1)*jet->Eta() );
    smearFactor_ *= 1. - parResp_.at(2)/(pt + parResp_.at(3));
  }

  // Apply resolution
  double smear = 1.;
  if( resolutionModel_ == Landau) {
    do {
      smear =  random_->Landau(1,sqrt(parReso_.at(0)*parReso_.at(0)/E/E +
				      parReso_.at(1)*parReso_.at(1)/E   +
				      parReso_.at(2)*parReso_.at(2))      );
    } while((smear < 0) || (smear > 2));
    smear = 2 - smear;
  }
  else if ( (resolutionModel_ == Gauss) ) {
    do {
      smear = random_->Gaus(1.0,sqrt(parReso_.at(0)*parReso_.at(0)/E/E +
				     parReso_.at(1)*parReso_.at(1)/E   +
				     parReso_.at(2)*parReso_.at(2))      );
    } while((smear < 0) || (smear > 2));
  }
  else if( resolutionModel_ == GaussUniform ) {
    do{
      smear = random_->Gaus( 1.0, parReso_.at(0) );
    } while( smear < 0 || smear > 2 );
    if( random_->Uniform() < parReso_.at(1) )
      smear = random_->Uniform();
  }
  else if( resolutionModel_ == TwoGauss ) {
    do {
      smear = histResp_->GetRandom();
    } while ( smear < 0.5 || smear > 1.2 );
  }

  smearFactor_ *= smear;
}




//!  \brief Calculate emf, scale hadronic tower energy with response factor,
//!         and smear with resolution
//!
//!  \param jet the four vector of the jet
//!  \param e True tower energy without pi0 part
//!  \param calcSmearFactor If true, the smear factor is calculated newly during this function call
//!  \param te Tower energy after scaling
//!  \param tem Measured em part of tower energy after scaling and smearing
//!  \param thad Measured had part of tower energy after scaling and smearing
//!  \param tout Measured HO part of tower energy after scaling and smearing
//!  \param temtrue True em part of tower energy after scaling
//!  \param thadtrue True had part of tower energy after scaling
//!  \param touttrue True HO part of tower energy after scaling
// -----------------------------------------------------------------
void ToyMC::smearTower(const TLorentzVector* jet, double e, bool calcSmearFactor, float& te, float& tem, float& thad, 
		       float& tout, float& temtrue, float& thadtrue, float& touttrue) 
{
  // Generate emf and set electromagnetic
  // and hadronic fraction; smear hadronic
  // fraction with response
  touttrue   = 0.;
  tout       = touttrue;
  float emf  = random_->Uniform(maxEmf_);
  temtrue    = emf * e;
  tem        = temtrue;
  thadtrue   = (1-emf) * e;
  thad       = thadtrue;

  // Apply response and resolution to hadronic fraction
  if( calcSmearFactor ) {
    if( smearTowersIndividually_ ) calculateSmearFactor(jet, thad);
    else                           calculateSmearFactor(jet, jet->E());
  }
  thad      *= smearFactor_;

  // Add up tower parts to total tower energy
  te = tem + thad + tout;
}



//!  \brief Split generated jet into towers
//!
//!  Splits generated true jet pt into 'chunks_' portions
//!  (particles), spreads them in eta and phi and sums up
//!  tower pt.
// -----------------------------------------------------------------
int ToyMC::splitJet(const TLorentzVector* jet ,float* et,float* eta,float * phi, int* ieta,int* iphi) {
  typedef __gnu_cxx::hash_map<int,int> TowerMap;
  TowerMap towers;
  double jphi = jet->Phi();
  if(jphi < 0) jphi += 2 * M_PI;
  //std::cout << "jet: Pt:" << jet.Pt() << " Phi:" << jet.Phi() << " Eta:" << jet.Eta() << '\n';
  //double de = jet.E() / chunks_;
  int ntowers = 0;
  TLorentzVector rec(0,0,0,0);
  TLorentzVector tow;
  double lostPt = 0;
  double dpt = jet->Pt() / chunks_;
  if(chunks_ < 0) {
    dpt = 0.3;
    chunks_ = (int)std::ceil(jet->Pt() / dpt);
  }
  for(int i = 0 ; i < chunks_ ; ++i) {
    //float teta = random_->Gaus(jet.Eta(), jetspread);
    //float tphi = random_->Gaus(jet.Phi(), jetspread);
    float R = 0.0;
    float PHI = 0.0;
    if(jetSpreadA_ != 0 || jetSpreadB_ != 0) {
      R = random_->Exp(1/(jetSpreadA_ +jetSpreadB_ * jet->E()));
      PHI = random_->Uniform(2 * M_PI);
    }
    //std::cout << "E:" << jet.E() << "  R:" << R << '\n';
    float teta = jet->Eta() + R * cos(PHI);
    float tphi = jet->Phi() + R * sin(PHI);
    tphi = TVector2::Phi_0_2pi(tphi);
    if(std::abs(teta) > 3.33333) {
      //std::cout << "chunk outside simulated calo\n";
      if(noOutOfCone_) --i;
      else lostPt += dpt;
      continue;
    }
    double de = dpt/cos(tphi-jphi);
    int ie, ip;
    calIds(teta, tphi, ie, ip); 
    //std::cout << "vorher:" << teta << ", " << tphi << ", " << ie << ", " 
    //	      << ip << "\n";
    float dphi = TVector2::Phi_mpi_pi(tphi-jphi);
    float deta = teta-jet->Eta();
    if(sqrt(deta*deta + dphi*dphi) > 0.5) {
      //std::cout << "Out of cone:" << teta << ":" << jet.Eta() << " , " << dphi << '\n';
      if(noOutOfCone_) --i;
      else lostPt += dpt;
      continue;
    }
    int id = 0;
    assert(ie*1000 + ip != 0);
    TowerMap::const_iterator towit = towers.find(ie*1000 + ip);
    if(towit != towers.end()) {
      id = towit->second;
    } else {
      ++ntowers;
      id = ntowers;
      towers[ie*1000 + ip] = id;
      et[id-1] = 0;
    }
    --id;
    tow.SetPtEtaPhiM(de,teta,tphi,0);
    //tow *= de/tow.E();
    et[id] += tow.Pt();
    eta[id] = teta;
    phi[id] = TVector2::Phi_mpi_pi(tphi);
    ieta[id] = ie;
    iphi[id] = ip;
    rec += tow;
  }
  //std::cout  << "Eta:" << jet.Eta() <<  "       : " << rec.Pt() << "," << rec.E() << "  == " << jet.Pt() << "," << jet.E() << '\n';
  //std::cout << "lost energy:" << lostE/jet.E() << '\n';
  if(noOutOfCone_) {
    double scale = jet->Pt()/rec.Pt();
    assert(scale < 1.1); 
    TLorentzVector rec2(0,0,0,0);
    for(int i = 0 ; i < ntowers ; ++i) {
      et[i] *= scale; 
      tow.SetPtEtaPhiM(et[i],eta[i],phi[i],0);
      rec2 += tow;
    }
    assert(std::abs((rec2.Pt()-jet->Pt())/jet->Pt()) < 0.001); 
    //std::cout << " vorher:" << scale << "  nachher:" << rec2.E()/jet.E() << "\n";
  }
  return ntowers;
}



// -----------------------------------------------------------------
int ToyMC::generateTrackClusterTree(TTree* CalibTree, int nevents) 
{
  //make tree
  const int kMaxTower = 1;
  int NobjTowCal;
  float towet[kMaxTower];
  float toweta[kMaxTower];
  float towphi[kMaxTower];
  float towen[kMaxTower];
  float towem[kMaxTower];
  float towhd[kMaxTower];
  float towoe[kMaxTower];
  float towemtrue[kMaxTower];
  float towhdtrue[kMaxTower];
  float towoetrue[kMaxTower];
  int towid_phi[kMaxTower];
  int towid_eta[kMaxTower];
  int towid [kMaxTower];
  float tracket;
  float tracketerr;
  float tracketa;
  float trackphi;
  float tracken;
  CalibTree->Branch("NobjTowCal",&NobjTowCal,"NobjTowCal/I");
  CalibTree->Branch("TowId",towid,"TowId[NobjTowCal]/I");
  CalibTree->Branch("TowId_phi",towid_phi,"TowId_phi[NobjTowCal]/I");
  CalibTree->Branch("TowId_eta",towid_eta,"TowId_eta[NobjTowCal]/I");
  CalibTree->Branch("TowEt",towet,"TowEt[NobjTowCal]/F");
  CalibTree->Branch("TowEta",toweta,"TowEta[NobjTowCal]/F");
  CalibTree->Branch("TowPhi",towphi,"TowPhi[NobjTowCal]/F");
  CalibTree->Branch("TowE",towen,"TowE[NobjTowCal]/F");
  CalibTree->Branch("TowEm",towem,"TowEm[NobjTowCal]/F");
  CalibTree->Branch("TowHad",towhd,"TowHad[NobjTowCal]/F");
  CalibTree->Branch("TowOE",towoe,"TowOE[NobjTowCal]/F");
  CalibTree->Branch("TowEmTrue",towemtrue,"TowEmTrue[NobjTowCal]/F");
  CalibTree->Branch("TowHadTrue",towhdtrue,"TowHadTrue[NobjTowCal]/F");
  CalibTree->Branch("TowOETrue",towoetrue,"TowOETrue[NobjTowCal]/F");
  CalibTree->Branch("TrackEt",&tracket,"TrackEt/F");
  CalibTree->Branch("TrackEterr",&tracketerr,"TrackEterr/F");
  CalibTree->Branch("TrackEta",&tracketa,"TrackEta/F");
  CalibTree->Branch("TrackPhi",&trackphi,"TrackPhi/F");
  CalibTree->Branch("TrackE",&tracken,"TrackE/F");

  for(int i = 0; i < nevents ; ++i) {
    genInput();
    tracket = pInput_->Pt();
    tracketerr = 0;
    tracketa = pInput_->Eta();
    trackphi = pInput_->Phi();
    tracken = pInput_->E();
    NobjTowCal = 1;
    towphi[0] = pInput_->Phi(); 
    toweta[0] = pInput_->Eta();
    towid[0] = 0;
    //std::cout << "vorher:" << toweta[0] << ", " << towphi[0] << ", " << towid_eta[0] << ", " 
    //	      << towid_phi[0] << "\n";
    calIds(toweta[0],towphi[0],towid_eta[0],towid_phi[0]);
    //std::cout << "nachher:" << toweta[0] << ", " << towphi[0] << ", " << towid_eta[0] << ", " 
    //	      << towid_phi[0] << "\n";

    // Calculate smear factor
    bool calcSmearFactor = false;
    if( i == 0 || smearTowersIndividually_ ) calcSmearFactor = true;

    smearTower(pInput_,pInput_->E(),calcSmearFactor,towen[0],towem[0],towhd[0],towoe[0],towemtrue[0],
	       towhdtrue[0],towoetrue[0]);
    towet[0] = towen[0]/pInput_->E() * pInput_->Pt();
    CalibTree->Fill();
    if(i % 1000 == 0) std::cout << "generated event " << i << '\n';
  }
  return CalibTree->GetEntriesFast(); 
}



// -----------------------------------------------------------------
int ToyMC::makeTrackCluster(const char* filename, int nevents) {
  TFile* file = new TFile(filename,"recreate");
  TTree* CalibTree = new TTree("TrackClusterTree","TrackClusterTre");
 
  nevents = generateTrackClusterTree(CalibTree, nevents);
  file->Write();
  file->Close();
  return  nevents;
}



//!  \brief Generate photon-jet events and write into tree
//!  \param CalibTree ROOT tree
//!  \param nevents Number of photon-jet events
//!  \return Number of generated events
// -----------------------------------------------------------------
int ToyMC::generatePhotonJetTree(TTree* CalibTree, int nevents)
{
  //make tree 
  const int kMaxTower = 1000;
  int NobjTowCal,NobjTrack = 0;
  float towet[kMaxTower];
  float toweta[kMaxTower];
  float towphi[kMaxTower];
  float towen[kMaxTower];
  float towem[kMaxTower];
  float towhd[kMaxTower];
  float towoe[kMaxTower]; 
  float towemtrue[kMaxTower];
  float towhdtrue[kMaxTower];
  float towoetrue[kMaxTower];
  int towid_phi[kMaxTower];
  int towid_eta[kMaxTower];
  int towid [kMaxTower];  
  const int kMaxTracks = 1;
  float trackpt[kMaxTracks];
  float tracketa[kMaxTracks];
  float trackphi[kMaxTracks];
  float trackp[kMaxTracks];
  float trackdr[kMaxTracks];
  float tracketaout[kMaxTracks];
  float trackphiout[kMaxTracks];
  float trackdrout[kMaxTracks];
  float trackemc1[kMaxTracks];
  float trackemc3[kMaxTracks];
  float trackemc5[kMaxTracks];
  float trackhac1[kMaxTracks];
  float trackhac3[kMaxTracks];
  float trackhac5[kMaxTracks];
  int tracktowid[kMaxTracks];
  int tracktowidphi[kMaxTracks];
  int tracktowideta[kMaxTracks];
  int trackid[kMaxTracks];
  int tracknhits[kMaxTracks];
  int trackQualityL[kMaxTracks];
  int trackQualityT[kMaxTracks];
  int trackQualityHP[kMaxTracks];
  float trackchi2[kMaxTracks];
  float muDR[kMaxTracks];
  float muDE[kMaxTracks];

  float jcalpt,jcalphi,jcaleta,jcalet,jcale;

  // All correction factors are 1 in ToyMC
  float jscaleZSP    = 1.;
  float jscalel2     = 1.;
  float jscalel3     = 1.;
  float jscalel23    = 1.;
  float jscaleJPT    = 1.;
  float jscalel23JPT = 1.;

  float jgenpt,jgenphi,jgeneta,jgenet,jgene;
  float mcalmet,mcalphi,mcalsum;
  float photonpt,photonphi,photoneta,photonet,photone; 
  float gphotonpt,gphotonphi,gphotoneta,gphotonet,gphotone;
  float nonleadingjetspt = 0.0;
  float weight = 1.0;
  CalibTree->Branch("NobjTowCal",&NobjTowCal,"NobjTowCal/I");   
  CalibTree->Branch("NobjTrack", &NobjTrack, "NobjTowCal/I");
  //tower branches
  CalibTree->Branch("TowId",towid,"TowId[NobjTowCal]/I");
  CalibTree->Branch("TowId_phi",towid_phi,"TowId_phi[NobjTowCal]/I");
  CalibTree->Branch("TowId_eta",towid_eta,"TowId_eta[NobjTowCal]/I");
  CalibTree->Branch("TowEt",towet,"TowEt[NobjTowCal]/F");
  CalibTree->Branch("TowEta",toweta,"TowEta[NobjTowCal]/F");
  CalibTree->Branch("TowPhi",towphi,"TowPhi[NobjTowCal]/F");
  CalibTree->Branch("TowE",towen,"TowE[NobjTowCal]/F");
  CalibTree->Branch("TowEm",towem,"TowEm[NobjTowCal]/F");
  CalibTree->Branch("TowHad",towhd,"TowHad[NobjTowCal]/F");
  CalibTree->Branch("TowOE",towoe,"TowOE[NobjTowCal]/F");
  CalibTree->Branch("TowEmTrue",towemtrue,"TowEmTrue[NobjTowCal]/F");
  CalibTree->Branch("TowHadTrue",towhdtrue,"TowHadTrue[NobjTowCal]/F");
  CalibTree->Branch("TowOETrue",towoetrue,"TowOETrue[NobjTowCal]/F");
  //track branches
  CalibTree->Branch( "NobjTrack",  &NobjTrack, "NobjTrack/I"             );
  CalibTree->Branch( "TrackTowId", tracktowid, "TrackTowId[NobjTrack]/I" );
  CalibTree->Branch( "TrackTowIdPhi", tracktowidphi, "TrackTowIdPhi[NobjTrack]/I" );
  CalibTree->Branch( "TrackTowIdEta", tracktowideta, "TrackTowIdEta[NobjTrack]/I" );
  CalibTree->Branch( "TrackId",    trackid,    "TrackId[NobjTrack]/I"    );
  CalibTree->Branch( "TrackNHits", tracknhits, "TrackNHits[NobjTrack]/I" );
  CalibTree->Branch( "TrackQualityL",trackQualityL,"TrackQualityL[NobjTrack]/O");
  CalibTree->Branch( "TrackQualityT",trackQualityT,"TrackQualityT[NobjTrack]/O");
  CalibTree->Branch( "TrackQualityHP",trackQualityHP,"TrackQualityHP[NobjTrack]/O");
  CalibTree->Branch( "TrackChi2",  trackchi2,  "TrackChi2[NobjTrack]/F"  );
  CalibTree->Branch( "TrackPt",    trackpt,    "TrackPt[NobjTrack]/F"    );
  CalibTree->Branch( "TrackEta",   tracketa,   "TrackEta[NobjTrack]/F"   );
  CalibTree->Branch( "TrackPhi",   trackphi,   "TrackPhi[NobjTrack]/F"   );
  CalibTree->Branch( "TrackP" ,    trackp,     "TrackP[NobjTrack]/F"     );
  CalibTree->Branch( "TrackDR" ,   trackdr,    "TrackDR[NobjTrack]/F"    );
  CalibTree->Branch( "TrackPhiOut",trackphiout,"TrackPhiout[NobjTrack]/F");
  CalibTree->Branch( "TrackEtaOut",tracketaout,"TrackEtaout[NobjTrack]/F");
  CalibTree->Branch( "TrackDROut", trackdrout, "TrackDRout[NobjTrack]/F" );
  CalibTree->Branch( "TrackEMC1",  trackemc1,  "TrackEMC1[NobjTrack]/F"  );
  CalibTree->Branch( "TrackEMC3",  trackemc3,  "TrackEMC3[NobjTrack]/F"  );
  CalibTree->Branch( "TrackEMC5",  trackemc5,  "TrackEMC5[NobjTrack]/F"  );
  CalibTree->Branch( "TrackHAC1",  trackhac1,  "TrackHAC1[NobjTrack]/F"  );
  CalibTree->Branch( "TrackHAC3",  trackhac3,  "TrackHAC3[NobjTrack]/F"  );
  CalibTree->Branch( "TrackHAC5",  trackhac5,  "TrackHAC5[NobjTrack]/F"  );
  CalibTree->Branch( "MuDR", muDR,  "MuDR[NobjTrack]/F"  );
  CalibTree->Branch( "MuDE", muDE,  "MuDE[NobjTrack]/F"  );

  // Jet- MEt-specific branches of the tree
  CalibTree->Branch("JetCalPt",&jcalpt,"JetCalPt/F");
  CalibTree->Branch("JetCalPhi",&jcalphi,"JetCalPhi/F");
  CalibTree->Branch("JetCalEta",&jcaleta,"JetCalEta/F");
  CalibTree->Branch("JetCalEt",&jcalet,"JetCalEt/F");
  CalibTree->Branch("JetCalE",&jcale,"JetCalE/F");  
  CalibTree->Branch( "JetCorrZSP",&jscaleZSP, "JetCorrZSP/F" );
  CalibTree->Branch( "JetCorrL2", &jscalel2,  "JetCorrL2/F" );
  CalibTree->Branch( "JetCorrL3", &jscalel3,  "JetCorrL3/F" );
  CalibTree->Branch( "JetCorrL2L3", &jscalel23,  "JetCorrL2L3/F" );
  CalibTree->Branch( "JetCorrJPT",&jscaleJPT, "JetCorrJPT/F" );
  CalibTree->Branch( "JetCorrL2L3JPT", &jscalel23JPT,  "JetCorrL2L3JPT/F" );

  // Gen- Jet- branches of the tree
  CalibTree->Branch("JetGenPt",&jgenpt,"JetGenPt/F");
  CalibTree->Branch("JetGenPhi",&jgenphi,"JetGenPhi/F");
  CalibTree->Branch("JetGenEta",&jgeneta,"JetGenEta/F");
  CalibTree->Branch("JetGenEt",&jgenet,"JetGenEt/F");
  CalibTree->Branch("JetGenE",&jgene,"JetGenE/F");
  //met
  CalibTree->Branch("MetCal",&mcalmet,"MetCal/F");
  CalibTree->Branch("MetCalPhi",&mcalphi,"MetCalPhi/F");
  CalibTree->Branch("MetCalSum",&mcalsum,"MetCalSum/F");
  //photons
  CalibTree->Branch("PhotonPt",&photonpt,"PhotonPt/F");
  CalibTree->Branch("PhotonPhi",&photonphi,"PhotonPhi/F");
  CalibTree->Branch("PhotonEta",&photoneta,"PhotonEta/F");
  CalibTree->Branch("PhotonEt",&photonet,"PhtonEt/F");
  CalibTree->Branch("PhotonE",&photone,"PhotonE/F");  
  // GenPhotons branches
  CalibTree->Branch( "GenPhotonPt",  &gphotonpt,  "GenPhotonPt/F"  );
  CalibTree->Branch( "GenPhotonPhi", &gphotonphi, "GenPhotonPhi/F" );
  CalibTree->Branch( "GenPhotonEta", &gphotoneta, "GenPhotonEta/F" );
  CalibTree->Branch( "GenPhotonEt",  &gphotonet,  "GenPhotonEt/F"   );
  CalibTree->Branch( "GenPhotonE",   &gphotone,   "GenPhotonE/F"   );
  // NonLeadingJetPt branch
  CalibTree->Branch( "NonLeadingJetPt", &nonleadingjetspt,   "NonLeadingJetPt/F"   );


  CalibTree->Branch("EventWeight",&weight,"EventWeight/F");
  

  // Generate events
  TLorentzVector jet;
  TLorentzVector genjet;
  TLorentzVector tower;

  // Loop over events
  for(int i = 0; i < nevents ; ++i) {
    // Generate truth 4-momentum
    genInput();

    // Assign it's variables to genphoton
    // and genphoton (perfectly measured)
    photonpt   = pInput_->Pt();
    photoneta  = pInput_->Eta();
    photonphi  = pInput_->Phi();
    photonet   = pInput_->Pt();
    photone    = pInput_->E();
    gphotonpt  = photonpt;
    gphotoneta = photoneta;
    gphotonphi = photonphi;
    gphotonet  = photonet;
    gphotone   = photone;

    // Gen jet
    jgenpt     = gphotonpt;
    jgeneta    = random_->Gaus(photoneta,1.0);
    if((jgeneta > 3.0) || (jgeneta < -3.0)) {
      --i;
      continue;
    }
    jgenphi    = photonphi + M_PI;

    // Create 4-momentum of genjet (massless)
    genjet.SetPtEtaPhiM(jgenpt,jgeneta,jgenphi,0);
    *pInput_= genjet;
    // Split it into towers and set tower truth
    NobjTowCal = splitJet(&genjet,towet,toweta,towphi,towid_eta,towid_phi);

    // Reset jet and genjet 4-momenta. They will
    // be repopulated by sum of tower 4-momenta
    jet.SetPtEtaPhiM(0,0,0,0);
    genjet.SetPtEtaPhiM(0,0,0,0);

    // Set random response parameters for this event
    // if required by response model
    if      (responseModel_ == Flat) parResp_.at(0) = random_->Uniform(1.5);
    else if (responseModel_ == Exp)  parResp_.at(0) = random_->Exp(0.5);
    else if (responseModel_ == Slope) {
      double u1 = random_->Uniform(2);
      double u2 = random_->Uniform(2);
      parResp_.at(0) = 2. - std::max(u1,u2);
    }

    // Generate pi0 fraction i.e. non-hadronic fraction
    double p0frac = random_->Uniform(maxPi0Frac_);

    // Loop over towers and smear truth with response factor
    for(int j = 0; j < NobjTowCal ; ++j) {
      // Set tower truth 4-momentum
      tower.SetPtEtaPhiM(towet[j],toweta[j],towphi[j],0);
      towen[j]  =  tower.E();

      // Add tower truth to genjet
      genjet   += tower;

      // Calculate smear factor
      bool calcSmearFactor = false;
      if( j == 0 || smearTowersIndividually_ ) calcSmearFactor = true;

      // Smear hadronic par of tower energy
      smearTower(pInput_, (1 - p0frac) * tower.E(),calcSmearFactor,towen[j],towem[j],towhd[j],towoe[j],
		 towemtrue[j],towhdtrue[j],towoetrue[j]); 

      // Add remaining em part
      towen[j]     += p0frac * tower.E();
      towem[j]     += p0frac * tower.E();
      towemtrue[j] += p0frac * tower.E();

      // Multiply tower 4-momentum with response
      tower        *= towen[j]/tower.E();
      towet[j]      = tower.Pt();

      // Add tower to jet
      jet += tower;
    } // End loop over towers

    // Set jet and genjet variables
    jcalpt   = jet.Pt();
    jcaleta  = jet.Eta();
    jcalphi  = jet.Phi();
    jcalet   = jet.Pt();
    jcale    = jet.E(); 
    jgenphi  = genjet.Phi();
    jgenet   = genjet.Pt();
    jgenpt   = genjet.Pt();
    jgene    = genjet.E();
    jgeneta  = genjet.Eta();
    mcalmet  = 0;
    mcalphi  = 0;
    mcalsum  = 0;

    if( !smearTowersIndividually_ ) {
      jscalel2  = 1. / smearFactor_;
      jscalel23 = jscalel2;
    }

    // Fill tree
    CalibTree->Fill();
    if(i % 1000 == 0) std::cout << "generating event " << i << '\n';
  } 
  return CalibTree->GetEntriesFast();
}



//!  \brief Generate dijet events and write into tree
//!
//!  The dijet event contains two real jets and a third
//!  dummy jet with Et = 0 GeV to fit to the DiJetReader
//!  which reads a third jet for potential cuts.
//!
//!  \param CalibTree ROOT tree
//!  \param nevents Number of dijet events
//!  \return Number of generated events
//----------------------------------------------------------
int ToyMC::generateDiJetTree(TTree* CalibTree, int nevents)
{
  //make tree 
  const int kMAX = 1000;
  int NobjTow;
  float towet[kMAX];
  float toweta[kMAX];
  float towphi[kMAX];
  float towen[kMAX];
  float towem[kMAX];
  float towhd[kMAX];
  float towoe[kMAX];  
  float towemtrue[kMAX];
  float towhdtrue[kMAX];
  float towoetrue[kMAX];
  int towid_phi[kMAX];
  int towid_eta[kMAX];
  int towid[kMAX];
  int tow_jetidx[kMAX];
  
  float ttowet[kMAX];
  float ttoweta[kMAX];
  float ttowphi[kMAX];
  int ttowid_phi[kMAX];
  int ttowid_eta[kMAX]; 

  int NobjTrack = 0;
  float trackpt[kMAX];
  float tracketa[kMAX];
  float trackphi[kMAX];
  float trackp[kMAX];
  float trackdr[kMAX];
  float tracketaout[kMAX];
  float trackphiout[kMAX];
  float trackdrout[kMAX];
  float trackemc1[kMAX];
  float trackemc3[kMAX];
  float trackemc5[kMAX];
  float trackhac1[kMAX];
  float trackhac3[kMAX];
  float trackhac5[kMAX];
  int tracktowid[kMAX];
  int tracktowidphi[kMAX];
  int tracktowideta[kMAX];
  int track_jetidx[kMAX];
  int trackid[kMAX]; // abs(PiD) if available, guess: muons only; =0: unknown
  int tracknhits[kMAX]; 
  bool trackQualityL[kMAX];
  bool trackQualityT[kMAX];
  bool trackQualityHP[kMAX];
  float trackchi2[ kMAX];
  float muDR[kMAX];
  float muDE[kMAX];

  const int kjMAX = 3;
  int NobjJet = 3;
  float jetpt[kjMAX];
  float jetphi[kjMAX];
  float jeteta[kjMAX];
  float jetet[kjMAX];
  float jete[kjMAX]; 

  float genevtscale;
  float jetgenpt[kjMAX];
  float jetgenphi[kjMAX];
  float jetgeneta[kjMAX];
  float jetgenet[kjMAX];
  float jetgene[kjMAX];

  float genJetColPt[kjMAX];
  float genJetColPhi[kjMAX];
  float genJetColEta[kjMAX];
  float genJetColEt[kjMAX];
  float genJetColE[kjMAX];
  int   genJetColJetIdx[kjMAX];

  float weight = 1; 

  // All correction factors are 1 in ToyMC
  float jscaleZSP[3]    = { 1., 1., 1. };
  float jscalel2[3]     = { 1., 1., 1. };
  float jscalel3[3]     = { 1., 1., 1. };
  float jscalel23[3]    = { 1., 1., 1. };
  float jscaleJPT[3]    = { 1., 1., 1. };
  float jscalel23JPT[3] = { 1., 1., 1. };

  // fill zero MET which is wrong....
  float mmet = 0;
  float mphi = 0;
  float msum = 0;

  // CaloTower branches
  CalibTree->Branch("NobjTow",&NobjTow,"NobjTow/I");
  CalibTree->Branch("TowId",towid,"TowId[NobjTow]/I");
  CalibTree->Branch("TowId_phi",towid_phi,"TowId_phi[NobjTow]/I");
  CalibTree->Branch("TowId_eta",towid_eta,"TowId_eta[NobjTow]/I");
  CalibTree->Branch("TowEt",towet,"TowEt[NobjTow]/F");
  CalibTree->Branch("TowEta",toweta,"TowEta[NobjTow]/F");
  CalibTree->Branch("TowPhi",towphi,"TowPhi[NobjTow]/F");
  CalibTree->Branch("TowE",towen,"TowE[NobjTow]/F");
  CalibTree->Branch("TowEm",towem,"TowEm[NobjTow]/F");
  CalibTree->Branch("TowHad",towhd,"TowHad[NobjTow]/F");
  CalibTree->Branch("TowOE",towoe,"TowOE[NobjTow]/F"); 
  CalibTree->Branch("TowEmTrue",towemtrue,"TowEmTrue[NobjTow]/F");
  CalibTree->Branch("TowHadTrue",towhdtrue,"TowHadTrue[NobjTow]/F");
  CalibTree->Branch("TowOETrue",towoetrue,"TowOETrue[NobjTow]/F");
  CalibTree->Branch("Tow_jetidx",tow_jetidx,"Tow_jetidx[NobjTow]/I");

  // track branches
  CalibTree->Branch("NobjTrack",&NobjTrack,"NobjTrack/I");
  CalibTree->Branch("TrackTowId",tracktowid,"TrackTowId[NobjTrack]/I");
  CalibTree->Branch("TrackTowIdPhi",tracktowidphi,"TrackTowIdPhi[NobjTrack]/I");
  CalibTree->Branch("TrackTowIdEta",tracktowideta,"TrackTowIdEta[NobjTrack]/I");
  CalibTree->Branch("TrackId",trackid,"TrackId[NobjTrack]/I");
  CalibTree->Branch("TrackNHits",tracknhits,"TrackNHits[NobjTrack]/I");
  CalibTree->Branch("TrackQualityL",trackQualityL,"TrackQualityL[NobjTrack]/O");
  CalibTree->Branch("TrackQualityT",trackQualityT,"TrackQualityT[NobjTrack]/O");
  CalibTree->Branch("TrackQualityHP",trackQualityHP,"TrackQualityHP[NobjTrack]/O");
  CalibTree->Branch("TrackChi2",trackchi2,"TrackChi2[NobjTrack]/F");
  CalibTree->Branch("TrackPt",trackpt,"TrackPt[NobjTrack]/F");
  CalibTree->Branch("TrackEta",tracketa,"TrackEta[NobjTrack]/F");
  CalibTree->Branch("TrackPhi",trackphi,"TrackPhi[NobjTrack]/F");
  CalibTree->Branch("TrackP",trackp,"TrackP[NobjTrack]/F");
  CalibTree->Branch("TrackDR",trackdr,"TrackDR[NobjTrack]/F");
  CalibTree->Branch("TrackPhiOut",trackphiout,"TrackPhiout[NobjTrack]/F");
  CalibTree->Branch("TrackEtaOut",tracketaout,"TrackEtaout[NobjTrack]/F");
  CalibTree->Branch("TrackDROut",trackdrout,"TrackDRout[NobjTrack]/F");
  CalibTree->Branch("TrackEMC1",trackemc1,"TrackEMC1[NobjTrack]/F");
  CalibTree->Branch("TrackEMC3",trackemc3,"TrackEMC3[NobjTrack]/F");
  CalibTree->Branch("TrackEMC5",trackemc5,"TrackEMC5[NobjTrack]/F");
  CalibTree->Branch("TrackHAC1",trackhac1,"TrackHAC1[NobjTrack]/F");
  CalibTree->Branch("TrackHAC3",trackhac3,"TrackHAC3[NobjTrack]/F");
  CalibTree->Branch("TrackHAC5",trackhac5,"TrackHAC5[NobjTrack]/F");
  CalibTree->Branch("Track_jetidx",track_jetidx,"Track_jetidx[NobjTrack]/I");
  CalibTree->Branch("MuDR",muDR,"MuDR[NobjTrack]/F");
  CalibTree->Branch("MuDE",muDE,"MuDE[NobjTrack]/F");
  // Jet-specific branches of the tree
  CalibTree->Branch("NobjJet",&NobjJet,"NobjJet/I"             );
  CalibTree->Branch("JetPt",jetpt,"JetPt[NobjJet]/F" );
  CalibTree->Branch("JetPhi",jetphi,"JetPhi[NobjJet]/F");
  CalibTree->Branch("JetEta",jeteta,"JetEta[NobjJet]/F");
  CalibTree->Branch("JetEt",jetet,"JetEt[NobjJet]/F" );
  CalibTree->Branch("JetE",jete,"JetE[NobjJet]/F"  );
  CalibTree->Branch("GenJetEt",jetgenet,"GenJetEt[NobjJet]/F" );
  CalibTree->Branch("GenJetPt",jetgenpt,"GenJetPt[NobjJet]/F" );
  CalibTree->Branch("GenJetE",jetgene,"GenJetE[NobjJet]/F"  );
  CalibTree->Branch("GenJetEta",jetgeneta,"GenJetEta[NobjJet]/F" );
  CalibTree->Branch("GenJetPhi",jetgenphi,"GenJetPhi[NobjJet]/F" );

  // Correction factors
  CalibTree->Branch( "JetCorrZSP",     jscaleZSP, "JetCorrZSP[NobjJet]/F" );
  CalibTree->Branch( "JetCorrL2",      jscalel2,  "JetCorrL2[NobjJet]/F" );
  CalibTree->Branch( "JetCorrL3",      jscalel3,  "JetCorrL3[NobjJet]/F" );
  CalibTree->Branch( "JetCorrL2L3",    jscalel23,  "JetCorrL2L3[NobjJet]/F" );
  CalibTree->Branch( "JetCorrJPT",     jscaleJPT, "JetCorrJPT[NobjJet]/F" );
  CalibTree->Branch( "JetCorrL2L3JPT", jscalel23JPT,  "JetCorrL2L3JPT[NobjJet]/F" );

  // Genjet collection
  CalibTree->Branch("NobjGenJet",&NobjJet,"NobjGenJet/I");
  CalibTree->Branch("GenJetColEt",genJetColEt,"GenJetColEt[NobjGenJet]/F" );
  CalibTree->Branch("GenJetColE",genJetColE,"GenJetColE[NobjGenJet]/F"     );
  CalibTree->Branch("GenJetColPt",genJetColPt,"GenJetColPt[NobjGenJet]/F" );
  CalibTree->Branch("GenJetColEta",genJetColEta,"GenJetColEta[NobjGenJet]/F" );
  CalibTree->Branch("GenJetColPhi",genJetColPhi,"GenJetColPhi[NobjGenJet]/F" );
  CalibTree->Branch("GenJetColJetIdx",genJetColJetIdx,"GenJetColJetIdx[NobjJet]/I");

  // Event branches
  CalibTree->Branch("GenEvtScale",&genevtscale,"GenEvtScale/F" );
  CalibTree->Branch("Weight",&weight,"Weight/F"  );

  // MET
  CalibTree->Branch("Met",&mmet,"Met/F");
  CalibTree->Branch("MetPhi",&mphi,"MetPhi/F");
  CalibTree->Branch("MetSum",&msum,"MetSum/F");


  // Generate events
  TLorentzVector jet[2];
  TLorentzVector orijet;
  TLorentzVector genjet[2];
  TLorentzVector tower;

  // Loop over events
  for(int n = 0; n < nevents ; ++n) {
    // Generate truth 4-momentum
    genInput();
    genevtscale  = pInput_->Pt();

    // Assign it's variables to first genjet
    jetgenpt[0]  = pInput_->Pt();
    jetgeneta[0] = pInput_->Eta();
    jetgenphi[0] = pInput_->Phi();
    jetgenet[0]  = pInput_->Pt();
    jetgene[0]   = pInput_->E();

    // Second genjet gets random eta
    // between -3 and 3, and phi+PI
    jetgenpt[1]  = jetgenpt[0];
    //    jetgeneta[1] = random_->Gaus(jetgeneta[0],1.0);
    jetgeneta[1] = random_->Uniform(minEta_,maxEta_);
//     if((jetgeneta[1] > 3.0) || (jetgeneta[1] < -3.0)) {
//       --n;
//       continue;
//     }
    jetgenphi[1] = jetgenphi[0] + M_PI;
    // Set random response paramters for this event
    // if required by response model
    if     (responseModel_ == Flat) parResp_.at(0) = random_->Uniform(1.5);
    else if(responseModel_ == Exp)  parResp_.at(0) = random_->Exp(0.5);
    else if(responseModel_ == Slope) {
      double u1 = random_->Uniform(2);
      double u2 = random_->Uniform(2);
      parResp_.at(0) = 2. - std::max(u1,u2);
    }

    // Loop over first two jets and
    // split genjets into towers
    NobjTow = 0;
    for(int i = 0 ; i < 2 ; ++i) {
      // Create 4-momentum of genjet (massless)
      orijet.SetPtEtaPhiM(jetgenpt[i],jetgeneta[i],jetgenphi[i],0);
      // Split it into towers and set truth of towers
      int ntow = splitJet(&orijet,ttowet,ttoweta,ttowphi,ttowid_eta,ttowid_phi);  
      // Reset jet and genjet 4-momenta. They will
      // be repopulated by sum of tower 4-momenta
      jet[i].SetPtEtaPhiM(0,0,0,0);
      genjet[i].SetPtEtaPhiM(0,0,0,0);

      // Generate pi0 fraction i.e. non-hadronic fraction
      double p0frac = random_->Uniform(maxPi0Frac_);

      // Loop over towers and smear truth with response
      // and resolution
      for(int j = 0 ; j < ntow ; ++j) {
	int k         = NobjTow + j;

	// Copy tower truth to tower
	towet[k]      = ttowet[j];
	toweta[k]     = ttoweta[j];
	towphi[k]     = ttowphi[j];	
	towid[k]      = 0;
	towid_eta[k]  = ttowid_eta[j];
	towid_phi[k]  = ttowid_phi[j];
	tow_jetidx[k] = i;
	tower.SetPtEtaPhiM(towet[k],toweta[k],towphi[k],0);
	towen[k]      = tower.E();

	// Add unsmeared tower to genjet
	genjet[i]    += tower;

	// Calculate smear factor
	bool calcSmearFactor = false;
	if( j == 0 || smearTowersIndividually_ ) calcSmearFactor = true;

	// Smear hadronic par of tower energy
	smearTower(&orijet, (1 - p0frac) * tower.E(),calcSmearFactor,towen[k],towem[k],towhd[k],
		   towoe[k],towemtrue[k],towhdtrue[k],towoetrue[k]); 

	// Add remaining em part
	towen[k]     += p0frac * tower.E();
	towem[k]     += p0frac * tower.E();
	towemtrue[k] += p0frac * tower.E();

	// Multiply tower 4-momentum with response
	tower        *= towen[k]/tower.E();
	towet[k]      = tower.Pt();

	// Add tower to jet
	jet[i]       += tower;
      } // End loop over towers

      // Set jet and genjet variables
      NobjTow        += ntow; 
      jetpt[i]        = jet[i].Pt();
      jetphi[i]       = jet[i].Phi();
      jeteta[i]       = jet[i].Eta();
      jetet[i]        = jet[i].Pt();
      jete[i]         = jet[i].E(); 
      jetgenpt[i]     = genjet[i].Pt();
      jetgenphi[i]    = genjet[i].Phi();
      jetgeneta[i]    = genjet[i].Eta();
      jetgenet[i]     = genjet[i].Pt();
      jetgene[i]      = genjet[i].E();

      if( !smearTowersIndividually_ ) {
	jscalel2[i]  = 1. / smearFactor_;
	jscalel23[i] = jscalel2[i];
      }

      genJetColPt[i]  = genjet[i].Pt();
      genJetColPhi[i] = genjet[i].Phi();
      genJetColEta[i] = genjet[i].Eta();
      genJetColEt[i]  = genjet[i].Pt();
      genJetColE[i]   = genjet[i].E();
      genJetColJetIdx[i] = i;
    } // End of loop over jets

    // Check generated eta measurement
    // of first jet
    if((jeteta[0] < minEta_) || (jeteta[0] > maxEta_)) {
      --n;
      continue;
    }

    // Order jets by pt
    //    if((jeteta[1] > minEta_) && (jeteta[1] < maxEta_)
    // && (jetpt[1] > jetpt[0])) {

    if( jetpt[1] > jetpt[0]) {

      //swap jets
      jetpt[0]     = jet[1].Pt();
      jetphi[0]    = jet[1].Phi();
      jeteta[0]    = jet[1].Eta();
      jetet[0]     = jet[1].Pt();
      jete[0]      = jet[1].E(); 
      jetgenpt[0]  = genjet[1].Pt();
      jetgenphi[0] = genjet[1].Phi();
      jetgeneta[0] = genjet[1].Eta();
      jetgenet[0]  = genjet[1].Pt();
      jetgene[0]   = genjet[1].E();
      genJetColJetIdx[0] = 1;

      double l2CorrTmp = jscalel2[0];
      jscalel2[0]   = jscalel2[1];
      jscalel23[0] = jscalel2[0];

      jetpt[1]     = jet[0].Pt();
      jetphi[1]    = jet[0].Phi();
      jeteta[1]    = jet[0].Eta();
      jetet[1]     = jet[0].Pt();
      jete[1]      = jet[0].E(); 
      jetgenpt[1]  = genjet[0].Pt();
      jetgenphi[1] = genjet[0].Phi();
      jetgeneta[1] = genjet[0].Eta();
      jetgenet[1]  = genjet[0].Pt();
      jetgene[1]   = genjet[0].E();
      genJetColJetIdx[1] = 0;

      jscalel2[1]   = l2CorrTmp;
      jscalel23[1] = jscalel2[1];

      // Swap towers
      for(int j = 0 ; j < NobjTow ; ++j) {
	tow_jetidx[j] = (tow_jetidx[j] == 0) ? 1 : 0;
      }

      // Set third dummy jet to zero
      jetpt[2]        = 0.;
      jetphi[2]       = 0.;
      jeteta[2]       = 0.;
      jetet[2]        = 0.;
      jete[2]         = 0.; 
      jetgenpt[2]     = 0.;
      jetgenphi[2]    = 0.;
      jetgeneta[2]    = 0.;
      jetgenet[2]     = 0.;
      jetgene[2]      = 0.;

      genJetColPt[2]  = 0.;
      genJetColPhi[2] = 0.;
      genJetColEta[2] = 0.;
      genJetColEt[2]  = 0.; 
      genJetColE[2]   = 0.;
      genJetColJetIdx[2] = 2;
    }

    // Fill tree
    CalibTree->Fill();

    if(n % 1000 == 0) std::cout << "generating event " << n << '\n';
  }  // End of loop over events

  return CalibTree->GetEntriesFast();
}

//!  \brief Generate w decay events and write into tree
//!
//!  The w decay event contains two real jets. and a third
//!  dummy jet with Et = 0 GeV to fit to the DiJetReader
//!  which reads a third jet for potential cuts.
//!
//!  \param CalibTree ROOT tree
//!  \param nevents Number of dijet events
//!  \return Number of generated events
//----------------------------------------------------------
int ToyMC::generateTopTree(TTree* CalibTree, int nevents)
{
  //make tree 
 // CaloTowers
  const int kMAX = 10000;
  int NobjTow;
  float towet[kMAX];
  float toweta[kMAX];
  float towphi[kMAX];
  float towen[kMAX];
  float towem[kMAX];
  float towhd[kMAX];
  float towoe[kMAX];  
  float towemtrue[kMAX];
  float towhdtrue[kMAX];
  float towoetrue[kMAX];
  int towid_phi[kMAX];
  int towid_eta[kMAX];
  int towid[kMAX];
  int tow_jetidx[kMAX];

  CalibTree->Branch( "NobjTow", &NobjTow,"NobjTow/I"             );
  CalibTree->Branch( "TowId",     towid,      "TowId[NobjTow]/I"    );
  CalibTree->Branch( "TowId_phi", towid_phi,  "TowId_phi[NobjTow]/I");
  CalibTree->Branch( "TowId_eta", towid_eta,  "TowId_eta[NobjTow]/I");
  CalibTree->Branch( "TowEt",     towet,      "TowEt[NobjTow]/F"    );
  CalibTree->Branch( "TowEta",    toweta,     "TowEta[NobjTow]/F"   );
  CalibTree->Branch( "TowPhi",    towphi,     "TowPhi[NobjTow]/F"   );
  CalibTree->Branch( "TowE",      towen,      "TowE[NobjTow]/F"     );
  CalibTree->Branch( "TowEm",     towem,      "TowEm[NobjTow]/F"    );
  CalibTree->Branch( "TowHad",    towhd,      "TowHad[NobjTow]/F"   );
  CalibTree->Branch( "TowOE",     towoe,      "TowOE[NobjTow]/F"    );
  CalibTree->Branch( "Tow_jetidx",tow_jetidx, "Tow_jetidx[NobjTow]/I");
  // RecoJets
  const int kjMAX = 8;
  int NobjJet = 2;;
  float jetpt[kjMAX];
  float jetphi[kjMAX];
  float jeteta[kjMAX];
  float jetet[kjMAX];
  float jete[kjMAX];
  int jetflavor[kjMAX];
  int jettopid[kjMAX];
  float jscaleL1[kjMAX];
  float jscaleL2[kjMAX];
  float jscaleL3[kjMAX];
  float jscaleL4[kjMAX];
  float jscaleL5[kjMAX];
  for(int i = 0; i < kjMAX; i++) {
    jscaleL1[i] = 1.;
    jscaleL2[i] = 1.;
    jscaleL3[i] = 1.;
    jscaleL4[i] = 1.;
    jscaleL5[i] = 1.;
  }
  CalibTree->Branch( "NobjJet",&NobjJet,"NobjJet/I"     );
  CalibTree->Branch( "JetPt", jetpt, "JetPt[NobjJet]/F" );
  CalibTree->Branch( "JetPhi",jetphi,"JetPhi[NobjJet]/F");
  CalibTree->Branch( "JetEta",jeteta,"JetEta[NobjJet]/F");
  CalibTree->Branch( "JetEt", jetet, "JetEt[NobjJet]/F" );
  CalibTree->Branch( "JetE",  jete,  "JetE[NobjJet]/F"  );
  CalibTree->Branch( "JetFlavor", jetflavor, "JetFlavor[NobjJet]/I" );
  CalibTree->Branch( "JetTopID",  jettopid,  "JetTopID[NobjJet]/I"  );
  CalibTree->Branch( "JetCorrL1", jscaleL1,  "JetCorrL1[NobjJet]/F" );
  CalibTree->Branch( "JetCorrL2", jscaleL2,  "JetCorrL2[NobjJet]/F" );
  CalibTree->Branch( "JetCorrL3", jscaleL3,  "JetCorrL3[NobjJet]/F" );
  CalibTree->Branch( "JetCorrL4", jscaleL4,  "JetCorrL4[NobjJet]/F" );
  CalibTree->Branch( "JetCorrL5", jscaleL5,  "JetCorrL5[NobjJet]/F" );
  // GenJets
  float genjetpt[kjMAX];
  float genjetphi[kjMAX];
  float genjeteta[kjMAX];
  float genjetet[kjMAX];
  float genjete[kjMAX];
  CalibTree->Branch( "GenJetPt", genjetpt, "GenJetPt[NobjJet]/F" );
  CalibTree->Branch( "GenJetPhi",genjetphi,"GenJetPhi[NobjJet]/F");
  CalibTree->Branch( "GenJetEta",genjeteta,"GenJetEta[NobjJet]/F");
  CalibTree->Branch( "GenJetEt", genjetet, "GenJetEt[NobjJet]/F" );
  CalibTree->Branch( "GenJetE",  genjete,  "GenJetE[NobjJet]/F"  );

  //EventWeight
  float weight = 1;
  CalibTree->Branch( "Weight",&weight,"Weight/F"   );

  // Loop over events
  TLorentzVector wjet[2];
  const double Wmass = 80.403; 
  float ttowet[kMAX];
  float ttoweta[kMAX];
  float ttowphi[kMAX];
  int ttowid_phi[kMAX];
  int ttowid_eta[kMAX]; 
  TLorentzVector jet, tower;
  for(int n = 0; n < nevents ; ++n) {
    // Generate truth 4-momentum of W boson
    genInput();	
    pInput_->SetVectM(pInput_->Vect(),Wmass);
    // Generate two jets
    wjet[0].SetPtEtaPhiM(0,0,0,0);
    wjet[1].SetPtEtaPhiM(0,0,0,0);    
    double theta = random_->Rndm() * M_PI;
    double phi = random_->Rndm() * M_PI * 2.0 - M_PI;
    double p = Wmass / 2;
    double pt = p * sin(theta);
    wjet[0].SetXYZM(pt * cos(phi),pt * sin(phi),p * cos(theta),0);
    wjet[1].SetVectM(-wjet[0].Vect(),0);
    // Boost jets according to momentum of W boson
    TVector3 b = pInput_->BoostVector();
    wjet[0].Boost(b);
    wjet[1].Boost(b);    

    // set gen jets from W
    for(int i = 0 ; i < 2 ; ++i) {
      genjetpt[i]  = wjet[i].Pt();
      genjeteta[i] = wjet[i].Eta();
      genjetphi[i] = wjet[i].Phi();
      genjetet[i]  = wjet[i].Pt();
      genjete[i]   = wjet[i].E();
    }
    if((std::abs(wjet[0].Eta()) > 3) || (std::abs(wjet[1].Eta()) > 3)) {
      n--;
      continue;
    }
  
    // Set random response paramters for this event
    // if required by response model
    if     (responseModel_ == Flat) parResp_.at(0) = random_->Uniform(1.5);
    else if(responseModel_ == Exp)  parResp_.at(0) = random_->Exp(0.5);
    else if(responseModel_ == Slope) {
      double u1 = random_->Uniform(2);
      double u2 = random_->Uniform(2);
      parResp_.at(0) = 2. - std::max(u1,u2);
    }

    // Loop over first two jets and
    // split genjets into towers
    NobjTow = 0;
    for(int i = 0 ; i < 2 ; ++i) {
      // Split it into towers and set truth of towers
      int ntow = splitJet(&wjet[i],ttowet,ttoweta,ttowphi,ttowid_eta,ttowid_phi);  

      // Reset jet and genjet 4-momenta. They will
      // be repopulated by sum of tower 4-momenta
      *pInput_ = wjet[i];
      jet.SetPtEtaPhiM(0,0,0,0);
      wjet[i].SetPtEtaPhiM(0,0,0,0);

      // Generate pi0 fraction i.e. non-hadronic fraction
      double p0frac = random_->Uniform(maxPi0Frac_);

      // Loop over towers and smear truth with response
      // and resolution
      for(int j = 0 ; j < ntow ; ++j) {
	int k         = NobjTow + j;

	// Copy tower truth to tower
	towet[k]      = ttowet[j];
	toweta[k]     = ttoweta[j];
	towphi[k]     = ttowphi[j];	
	towid[k]      = 0;
	towid_eta[k]  = ttowid_eta[j];
	towid_phi[k]  = ttowid_phi[j];
	tow_jetidx[k] = i;
	tower.SetPtEtaPhiM(towet[k],toweta[k],towphi[k],0);
	towen[k]      = tower.E();

	// Add unsmeared tower to genjet
	wjet[i]    += tower;

	// Calculate smear factor
	bool calcSmearFactor = false;
	if( j == 0 || smearTowersIndividually_ ) calcSmearFactor = true;

	// Smear hadronic par of tower energy
	smearTower(pInput_,(1 - p0frac) * tower.E(),calcSmearFactor,towen[k],towem[k],towhd[k],
		   towoe[k],towemtrue[k],towhdtrue[k],towoetrue[k]); 

	// Add remaining em part
	towen[k]     += p0frac * tower.E();
	towem[k]     += p0frac * tower.E();
	towemtrue[k] += p0frac * tower.E();

	// Multiply tower 4-momentum with response
	tower        *= towen[k]/tower.E();
	towet[k]      = tower.Pt();

	// Add tower to jet
	jet += tower;
      } // End loop over towers

      // Set jet and genjet variables
      NobjTow        += ntow; 
      jetpt[i]        = jet.Pt();
      jetphi[i]       = jet.Phi();
      jeteta[i]       = jet.Eta();
      jetet[i]        = jet.Pt();
      jete[i]         = jet.E();
      jetflavor[i]    = 1;//uds
      jettopid[i]      = 0;
      /*
	jetgenpt[i]     = genjet[i].Pt();
	jetgenphi[i]    = genjet[i].Phi();
	jetgeneta[i]    = genjet[i].Eta();
	jetgenet[i]     = genjet[i].Pt();
	jetgene[i]      = genjet[i].E();
      */
    } // End of loop over jets

    // Fill tree
    CalibTree->Fill();
    if(n % 1000 == 0) std::cout << "generating event " << n << '\n';
  }  // End of loop over events

  return CalibTree->GetEntriesFast();
}

//!  \brief Generate photon-jet events and write into tree
//!  \param filename Name of output file
//!  \param nevents Number of photon-jet events
//!  \return Number of generated events
//----------------------------------------------------------
int ToyMC::makePhotonJet(const char* filename, int nevents) {
  TFile* file = new TFile(filename,"recreate");
  TTree* CalibTree = new TTree("GammaJetTree","GammaJetTree");

  nevents = generatePhotonJetTree(CalibTree,nevents);
  //file->ls();
  file->Write();
  file->Close();
  return nevents;
}



//!  \brief Generate dijet events
//!  \param filename Name of output file
//!  \param nevents Number of dijet events
//!  \return Number of generated events
//----------------------------------------------------------
int ToyMC::makeDiJet(const char* filename, int nevents) {
  TFile* file = new TFile(filename,"recreate");
  TTree* CalibTree = new TTree("DiJetTree","DiJetTree");

  nevents = generateDiJetTree(CalibTree,nevents);
  //file->ls();
  file->Write();
  file->Close();
  return nevents;
}

//!  \brief Generate top events
//!  \param filename Name of output file
//!  \param nevents Number of dijet events
//!  \return Number of generated events
//----------------------------------------------------------
int ToyMC::makeTop(const char* filename, int nevents) {
  TFile* file = new TFile(filename,"recreate");
  TTree* CalibTree = new TTree("TopTree","TopTree");

  nevents = generateTopTree(CalibTree,nevents);
  //file->ls();
  file->Write();
  file->Close();
  return nevents;
}



//----------------------------------------------------------
void ToyMC::init(const std::string& configfile) {
  ConfigFile config(configfile.c_str());
  init(&config);
}

//----------------------------------------------------------
void ToyMC::init(const ConfigFile* config) {

  // Ranges
  minEta_ = config->read<double>("ToyMC min eta",-2.5);
  maxEta_ = config->read<double>("ToyMC max eta",2.5);
  minPt_  = config->read<double>("ToyMC min pt",30);
  maxPt_  = config->read<double>("ToyMC max pt",400);

  // Truth spectrum
  std::string spectrum = config->read<std::string>("ToyMC pt spectrum","uniform");
  parTruth_            = bag_of<double>(config->read<string>("ToyMC pt spectrum parameters",";"));
  if(spectrum == "powerlaw") {
    ptSpectrum_ = PowerLaw; 
    assert( parTruth_.size() > 0 );
  } else if(spectrum == "uniform") {
    ptSpectrum_ = Uniform;
  } else if(spectrum == "PtEtaHistogram") {
    ptSpectrum_ = PtEtaHistogram;
    TFile file("toymcPtEtaInput.root");
    histPtEta_ = (TH2F*) file.Get("genWPtEta");
    histPtEta_->SetDirectory(0);
  } else {
    std::cerr << "unknown ToyMC pt spectrum:" << spectrum << '\n';
    exit(1);
  }

  // Response model
  parResp_ = bag_of<double>(config->read<string>("ToyMC response parameters","1"));
  std::string response = config->read<std::string>("ToyMC response model","Constant");
  if        ( response == "Constant" ) {
    responseModel_ = Constant;
    assert( parResp_.size() >= 1 );
  } else if ( response == "Flat" ) {
    responseModel_ = Flat;
    assert( parResp_.size() >= 1 );
  } else if ( response == "Slope" ) {
    responseModel_ = Slope;
    assert( parResp_.size() >= 1 );
  } else if ( response == "Exp" ) {
    responseModel_ = Exp;
    assert( parResp_.size() >= 1 );
  } else if ( response == "L3" ) {
    responseModel_ = L3;
    assert( parResp_.size() >= 4 );
  } else if( response == "SimpleInverse" ) {
    responseModel_ = SimpleInverse;
    assert( parResp_.size() >= 2 );
  } else if( response == "StepEta" ) {
    responseModel_ = StepEta;
    assert( parResp_.size() >= 1 );
  } else if( response == "SinusEta" ) {
    responseModel_ = SinusEta;
    assert( parResp_.size() >= 1 );
  } else if( response == "SinusEtaSimpleInversePt" ) {
    responseModel_ = SinusEtaSimpleInversePt;
    assert( parResp_.size() >= 3 );
  } else {
    std::cerr << "unknown ToyMC response model: " << response << '\n';
    exit(1);
  }

  // Resolution model
  parReso_               = bag_of<double>(config->read<string>("ToyMC resolution parameters","4.44 1.11 0.03"));
  std::string resolution = config->read<std::string>("ToyMC resolution model","Gauss");
  if(resolution == "Gauss") {
    resolutionModel_ = Gauss;
    assert( parReso_.size() >= 3 );
  } else if(resolution  == "Landau") {
    resolutionModel_ = Landau;
    assert( parReso_.size() >= 3 );
  } else if(resolution  == "Dirac") {
    resolutionModel_ = Dirac;
  } else if( resolution == "GaussUniform" ) {
    resolutionModel_ = GaussUniform; 
    assert( parReso_.size() >= 3 );
  } else if( resolution == "TwoGauss" ) {
    resolutionModel_ = TwoGauss;
    assert( parReso_.size() >= 4 );
    
    // Set up sum of two Gaussians as pdf
    double c  = parReso_.at(0);  // Normalization
    double u0 = 1.;              // Mean of central Gaussian (scale)
    double s0 = parReso_.at(1);  // Width of central Gaussian
    double u1 = parReso_.at(2);  // Mean of second Gaussian
    double s1 = parReso_.at(3);  // Width of central Gaussian

    double minResp = 0.;
    double maxResp = 2.;

    TF1 * f = new TF1("f","gaus(0)+gaus(3)",minResp,maxResp);
    f->SetParameter(0,c/sqrt(2*M_PI)/s0);
    f->SetParameter(1,u0);
    f->SetParameter(2,s0);
    f->SetParameter(3,(1.-c)/sqrt(2*M_PI)/s1);
    f->SetParameter(4,u1);
    f->SetParameter(5,s1);

    // Fill response histogram according to f
    histResp_ = new TH1F("hHistResp",";p^{jet}_{T} / p^{true}_{T};1/(Nw) dN / d(p^{jet}_{T} / p^{true}_{T})",			     1000,minResp,maxResp);
    for(int bin = 1; bin <= histResp_->GetNbinsX(); bin++) {
      double r = f->Eval(histResp_->GetBinCenter(bin));
      histResp_->SetBinContent(bin,r);
    }

    double norm = histResp_->Integral("width");
    histResp_->Scale(1./norm);
    delete f;
   } else {
     std::cerr << "unknown ToyMC resolution model: " << resolution << '\n';
     exit(1);
   }

  // Calculate smear factor for each tower or each jet
  smearTowersIndividually_ = config->read<bool>("ToyMC smear towers individually",false);

  // Jets
  jetSpreadA_  = config->read<double>("ToyMC jet spread A",0.5);
  jetSpreadB_  = config->read<double>("ToyMC jet spread B",0);
  noOutOfCone_ = config->read<bool>("ToyMC avoid out-of-cone",true);
  useTowerCenterEtaPhi_ = config->read<bool>("ToyMC use eta and phi at tower center",true);
  chunks_      = config->read<int>("ToyMC chunks",200);
  maxPi0Frac_  = config->read<double>("ToyMC max pi0 fraction",0.5);
  maxEmf_      = config->read<double>("ToyMC tower max EMF",0.5);
  
  // General
  int seed = config->read<int>("ToyMC seed",0); 
  random_->SetSeed(seed);
  type_ = config->read<int>("ToyMC type",1);
  if( !( type_ == 1 || type_ == 2 || type_ == 3 ) ) {
    std::cout << "unknown ToyMC event type " << type_ << std::endl;
    exit(1);
  }
}



//----------------------------------------------------------
void ToyMC::print() const {
  std::cout << "\n  ToyMC configuration:\n";
  std::cout << " -----------------------------------------\n";
  std::cout << "  primary:      " << minEta_ << " < eta < " << maxEta_ << '\n';
  std::cout << "                " << minPt_ << " < pt < " << maxPt_ << "\n";

  std::cout << "  spectrum:     ";
  if( ptSpectrum_ == Uniform )
    std::cout << "Uniform\n";
  else if( ptSpectrum_ == PowerLaw )
    std::cout << "PowerLaw: 1/pt^" << parTruth_.at(0) << std::endl;
  else if( ptSpectrum_ == PtEtaHistogram )
    std::cout << "PtEtaHistogram\n";

  std::cout << "  response:     ";
  if( responseModel_ == Constant )
    std::cout << "'Constant'";
  else if( responseModel_ == Flat )
    std::cout << "'Flat'";
  else if( responseModel_ == Exp )
    std::cout << "'Exp'";
  else if( responseModel_ == Slope )
    std::cout << "'Slope'";
  else if( responseModel_ == L3 )
    std::cout << "'L3'";
  else if( responseModel_ == SimpleInverse )
    std::cout << "'SimpleInverse'";
  else if( responseModel_ == StepEta )
    std::cout << "'StepEta'";
  else if( responseModel_ == SinusEta )
    std::cout << "'SinusEta'";
  else if( responseModel_ == SinusEtaSimpleInversePt )
    std::cout << "'SinusEtaSimpleInversePt'";

  std::cout << " with parameters ";
  for(unsigned int i = 0; i < parResp_.size(); i++)
    std::cout << parResp_.at(i) << ",  ";
  std::cout << "\n";

  std::cout << "  resolution:   ";
  if( resolutionModel_ == Gauss )
    std::cout << "'Gauss'";
  else if( resolutionModel_ == Landau )
    std::cout << "'Landau'";
  else if( resolutionModel_ == Dirac )
    std::cout << "'Dirac'";
  else if( resolutionModel_ == GaussUniform )
    std::cout << "'GaussUniform'";
  else if( resolutionModel_ == TwoGauss )
    std::cout << "'TwoGauss'";

  std::cout << " with parameters ";
  for(unsigned int i = 0; i < parReso_.size(); i++)
    std::cout << parReso_.at(i) << ",  ";
  std::cout << "\n";

  std::cout << "  Smear factor is calculated for each ";
  if( smearTowersIndividually_ ) std::cout << "tower\n";
  else                           std::cout << "jet\n";

  std::cout << "  max EMF:      " << maxEmf_ << "\n";
  std::cout << "  max pi0:      " << maxPi0Frac_ << "\n";

  std::cout << "  jet spread:   ";
  std::cout << "A = " << jetSpreadA_ << ",  B = " << jetSpreadB_ << "\n";
  std::cout << "  n chunks:     " << chunks_ << "\n";

  std::cout << "  use eta and phi at tower center:  ";
  if( useTowerCenterEtaPhi_ )
    std::cout << "yes\n";    
  else
    std::cout << "no\n";

  std::cout << "  out-of-cone:  ";
  if( noOutOfCone_ )
    std::cout << "no\n";
  else 
    std::cout << "yes\n";

  std::cout << "  type:         ";
  if( type_ == 1 )
    std::cout << "Photon-jet events\n";
  else if( type_ == 2 )
    std::cout << "Dijet events\n";
  else if( type_ == 3 ) 
    std::cout << "Top events\n";

  std::cout << "\n";
}



// -----------------------------------------------------------------
int ToyMC::etaBin(float eta) const {
  assert( eta > -5.191 && eta < 5.191 ); 

  int etaBin = 10000;

  if( eta < -4.889 )      etaBin = -41;
  else if( eta < -4.716 ) etaBin = -40;
  else if( eta < -4.538 ) etaBin = -39;
  else if( eta < -4.363 ) etaBin = -38;
  else if( eta < -4.191 ) etaBin = -37;
  else if( eta < -4.013 ) etaBin = -36;
  else if( eta < -3.839 ) etaBin = -35;
  else if( eta < -3.664 ) etaBin = -34;
  else if( eta < -3.489 ) etaBin = -33;
  else if( eta < -3.314 ) etaBin = -32;
  else if( eta < -3.139 ) etaBin = -31;
  else if( eta < -2.964 ) etaBin = -30;
  else if( eta < -2.853 ) etaBin = -29; 
  else if( eta <  -2.65 ) etaBin = -28;
  else if( eta <   -2.5 ) etaBin = -27;
  else if( eta < -2.322 ) etaBin = -26;
  else if( eta < -2.172 ) etaBin = -25;
  else if( eta < -2.043 ) etaBin = -24;
  else if( eta <  -1.93 ) etaBin = -23;
  else if( eta <  -1.83 ) etaBin = -22;
  else if( eta <  -1.74 ) etaBin = -21;
  else if( eta < -1.653 ) etaBin = -20;
  else if( eta < -1.566 ) etaBin = -19;
  else if( eta < -1.479 ) etaBin = -18;
  else if( eta < -1.392 ) etaBin = -17;
  else if( eta < -1.305 ) etaBin = -16;
  else if( eta < -1.218 ) etaBin = -15;
  else if( eta < -1.131 ) etaBin = -14;
  else if( eta < -1.044 ) etaBin = -13;
  else if( eta < -0.957 ) etaBin = -12;
  else if( eta < -0.879 ) etaBin = -11;
  else if( eta < -0.783 ) etaBin = -10;
  else if( eta < -0.696 ) etaBin = -9;
  else if( eta < -0.609 ) etaBin = -8;
  else if( eta < -0.522 ) etaBin = -7;
  else if( eta < -0.435 ) etaBin = -6;
  else if( eta < -0.348 ) etaBin = -5;
  else if( eta < -0.261 ) etaBin = -4;
  else if( eta < -0.174 ) etaBin = -3;
  else if( eta < -0.087 ) etaBin = -2;
  else if( eta <      0 ) etaBin = -1;
  else if( eta <  0.087 ) etaBin =  1;
  else if( eta <  0.174 ) etaBin =  2;
  else if( eta <  0.261 ) etaBin =  3;
  else if( eta <  0.348 ) etaBin =  4;
  else if( eta <  0.435 ) etaBin =  5;
  else if( eta <  0.522 ) etaBin =  6;
  else if( eta <  0.609 ) etaBin =  7;
  else if( eta <  0.696 ) etaBin =  8;
  else if( eta <  0.783 ) etaBin =  9;
  else if( eta <  0.879 ) etaBin =  10;
  else if( eta <  0.957 ) etaBin =  11;
  else if( eta <  1.044 ) etaBin =  12;
  else if( eta <  1.131 ) etaBin =  13;
  else if( eta <  1.218 ) etaBin =  14;
  else if( eta <  1.305 ) etaBin =  15;
  else if( eta <  1.392 ) etaBin =  16;
  else if( eta <  1.479 ) etaBin =  17;
  else if( eta <  1.566 ) etaBin =  18;
  else if( eta <  1.653 ) etaBin =  19;
  else if( eta <   1.74 ) etaBin =  20;
  else if( eta <   1.83 ) etaBin =  21;
  else if( eta <   1.93 ) etaBin =  22;
  else if( eta <  2.043 ) etaBin =  23;
  else if( eta <  2.172 ) etaBin =  24;
  else if( eta <  2.322 ) etaBin =  25;
  else if( eta <    2.5 ) etaBin =  26;
  else if( eta <   2.65 ) etaBin =  27;
  else if( eta <  2.853 ) etaBin =  28;
  else if( eta <  2.964 ) etaBin =  29;
  else if( eta <  3.139 ) etaBin =  30;
  else if( eta <  3.314 ) etaBin =  31;
  else if( eta <  3.489 ) etaBin =  32;
  else if( eta <  3.664 ) etaBin =  33;
  else if( eta <  3.839 ) etaBin =  34;
  else if( eta <  4.013 ) etaBin =  35;
  else if( eta <  4.191 ) etaBin =  36;
  else if( eta <  4.363 ) etaBin =  37;
  else if( eta <  4.538 ) etaBin =  38;
  else if( eta <  4.716 ) etaBin =  39;
  else if( eta <  4.889 ) etaBin =  40;
  else if( eta <  5.191 ) etaBin =  41;

  return etaBin;
}



// -----------------------------------------------------------------
float ToyMC::etaBinEdge(int etaBin, bool lowerEdge) const {
  assert( etaBin >= -41 && etaBin <= 41 );
  // return eta bin - eta edge mappting
  switch(etaBin){
  case -41: return (lowerEdge ? -5.191 : -4.889); break;
  case -40: return (lowerEdge ? -4.889 : -4.716); break;
  case -39: return (lowerEdge ? -4.716 : -4.538); break;
  case -38: return (lowerEdge ? -4.538 : -4.363); break;
  case -37: return (lowerEdge ? -4.363 : -4.191); break;
  case -36: return (lowerEdge ? -4.191 : -4.013); break;
  case -35: return (lowerEdge ? -4.013 : -3.839); break;
  case -34: return (lowerEdge ? -3.839 : -3.664); break;
  case -33: return (lowerEdge ? -3.664 : -3.489); break;
  case -32: return (lowerEdge ? -3.489 : -3.314); break;
  case -31: return (lowerEdge ? -3.314 : -3.139); break;
  case -30: return (lowerEdge ? -3.139 : -2.964); break;
  case -29: return (lowerEdge ? -2.964 : -2.853); break; 
  case -28: return (lowerEdge ? -2.853 :  -2.65); break;
  case -27: return (lowerEdge ?  -2.65 :   -2.5); break;
  case -26: return (lowerEdge ?   -2.5 : -2.322); break;
  case -25: return (lowerEdge ? -2.322 : -2.172); break;
  case -24: return (lowerEdge ? -2.172 : -2.043); break;
  case -23: return (lowerEdge ? -2.043 :  -1.93); break;
  case -22: return (lowerEdge ?  -1.93 :  -1.83); break;
  case -21: return (lowerEdge ?  -1.83 :  -1.74); break;
  case -20: return (lowerEdge ?  -1.74 : -1.653); break;
  case -19: return (lowerEdge ? -1.653 : -1.566); break;
  case -18: return (lowerEdge ? -1.566 : -1.479); break;
  case -17: return (lowerEdge ? -1.479 : -1.392); break;
  case -16: return (lowerEdge ? -1.392 : -1.305); break;
  case -15: return (lowerEdge ? -1.305 : -1.218); break;
  case -14: return (lowerEdge ? -1.218 : -1.131); break;
  case -13: return (lowerEdge ? -1.131 : -1.044); break;
  case -12: return (lowerEdge ? -1.044 : -0.957); break;
  case -11: return (lowerEdge ? -0.957 : -0.879); break;
  case -10: return (lowerEdge ? -0.879 : -0.783); break;
  case  -9: return (lowerEdge ? -0.783 : -0.696); break;
  case  -8: return (lowerEdge ? -0.696 : -0.609); break;
  case  -7: return (lowerEdge ? -0.609 : -0.522); break;
  case  -6: return (lowerEdge ? -0.522 : -0.435); break;
  case  -5: return (lowerEdge ? -0.435 : -0.348); break;
  case  -4: return (lowerEdge ? -0.348 : -0.261); break;
  case  -3: return (lowerEdge ? -0.261 : -0.174); break;
  case  -2: return (lowerEdge ? -0.174 : -0.087); break;
  case  -1: return (lowerEdge ? -0.087 :      0); break;
  case  +1: return (lowerEdge ?      0 :  0.087); break;
  case  +2: return (lowerEdge ?  0.087 :  0.174); break;
  case  +3: return (lowerEdge ?  0.174 :  0.261); break;
  case  +4: return (lowerEdge ?  0.261 :  0.348); break;
  case  +5: return (lowerEdge ?  0.348 :  0.435); break;
  case  +6: return (lowerEdge ?  0.435 :  0.522); break;
  case  +7: return (lowerEdge ?  0.522 :  0.609); break;
  case  +8: return (lowerEdge ?  0.609 :  0.696); break;
  case  +9: return (lowerEdge ?  0.696 :  0.783); break;
  case +10: return (lowerEdge ?  0.783 :  0.879); break;
  case +11: return (lowerEdge ?  0.879 :  0.957); break;
  case +12: return (lowerEdge ?  0.957 :  1.044); break;
  case +13: return (lowerEdge ?  1.044 :  1.131); break;
  case +14: return (lowerEdge ?  1.131 :  1.218); break;
  case +15: return (lowerEdge ?  1.218 :  1.305); break;
  case +16: return (lowerEdge ?  1.305 :  1.392); break;
  case +17: return (lowerEdge ?  1.392 :  1.479); break;
  case +18: return (lowerEdge ?  1.479 :  1.566); break;
  case +19: return (lowerEdge ?  1.566 :  1.653); break;
  case +20: return (lowerEdge ?  1.653 :   1.74); break;
  case +21: return (lowerEdge ?   1.74 :   1.83); break;
  case +22: return (lowerEdge ?   1.83 :   1.93); break;
  case +23: return (lowerEdge ?   1.93 :  2.043); break;
  case +24: return (lowerEdge ?  2.043 :  2.172); break;
  case +25: return (lowerEdge ?  2.172 :  2.322); break;
  case +26: return (lowerEdge ?  2.322 :    2.5); break;
  case +27: return (lowerEdge ?    2.5 :   2.65); break;
  case +28: return (lowerEdge ?   2.65 :  2.853); break;
  case +29: return (lowerEdge ?  2.853 :  2.964); break;
  case +30: return (lowerEdge ?  2.964 :  3.139); break;
  case +31: return (lowerEdge ?  3.139 :  3.314); break;
  case +32: return (lowerEdge ?  3.314 :  3.489); break;
  case +33: return (lowerEdge ?  3.489 :  3.664); break;
  case +34: return (lowerEdge ?  3.664 :  3.839); break;
  case +35: return (lowerEdge ?  3.839 :  4.013); break;
  case +36: return (lowerEdge ?  4.013 :  4.191); break;
  case +37: return (lowerEdge ?  4.191 :  4.363); break;
  case +38: return (lowerEdge ?  4.363 :  4.538); break;
  case +39: return (lowerEdge ?  4.538 :  4.716); break;
  case +40: return (lowerEdge ?  4.716 :  4.889); break;
  case +41: return (lowerEdge ?  4.889 :  5.191); break;
    //something went wrong;
  default : return -1; break;
  }
}

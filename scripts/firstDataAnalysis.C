#include <algorithm>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TH1D.h"
#include "TString.h"



// ===== Type declarations =====
// =============================

class Jet;
class Event;
class Sample;

typedef std::vector<Event*> Data;
typedef std::vector<Event*>::iterator DataIt;

typedef std::vector<TH1D*> Distributions;
typedef std::vector<TH1D*>::iterator DistIt;




// ===== Class implementations =====
// =================================

// ----- Jet -----
class Jet {
 public:
  static bool ptGreaterThan(const Jet *j1, const Jet *j2);
  static bool corrPtGreaterThan(const Jet *j1, const Jet *j2);

  Jet(double pt,
      double eta,
      double phi,
      double emf,
      int n90Hits,
      double fHPD,
      double fRBX,
      double corrL2L3);

  double pt() const { return pt_; }
  double eta() const { return eta_; }
  double phi() const { return phi_; }
  double emf() const { return emf_; }
  int n90Hits() const { return n90Hits_; }
  double fHPD() const { return fHPD_; }
  double fRBX() const { return fRBX_; }
  double corrL2L3Pt() const { return corrL2L3Pt_; }

 private:
  const double pt_;
  const double eta_;
  const double phi_;
  const double emf_;
  const int n90Hits_;
  const double fHPD_;
  const double fRBX_;
  const double corrL2L3Pt_;
};

Jet::Jet(double pt,
	 double eta,
	 double phi,
	 double emf,
	 int n90Hits,
	 double fHPD,
	 double fRBX,
	 double corrL2L3)
  : pt_(pt),
    eta_(eta),
    phi_(phi),
    emf_(emf),
    n90Hits_(n90Hits),
    fHPD_(fHPD),
    fRBX_(fRBX),
    corrL2L3Pt_(corrL2L3*pt) {};

//! For sorting jets in calo pt
bool Jet::ptGreaterThan(const Jet *j1, const Jet *j2) {
  // check for 0
  if (j1 == 0) {
    return j2 != 0;
  } else if (j2 == 0) {
    return false;
  } else {
    return j1->pt() > j2->pt();
  }
}

//! For sorting jets in corrected pt
bool Jet::corrPtGreaterThan(const Jet *j1, const Jet *j2) {
  // check for 0
  if (j1 == 0) {
    return j2 != 0;
  } else if (j2 == 0) {
    return false;
  } else {
    return j1->corrL2L3Pt() > j2->corrL2L3Pt();
  }
}


// ----- Event -----
class Event {
 public:
  Event(unsigned int runNumber,
	unsigned int lumiBlockNumber,
	unsigned int vtxNTracks,
	double vtxPosZ,
	double sumEt,
	double met,
	const std::vector<Jet*>& jets);
  ~Event();

  unsigned int runNumber() const { return runNumber_; }
  unsigned int lumiBlockNumber() const { return lumiBlockNumber_; }
  unsigned int vtxNTracks() const { return vtxNTracks_; }
  double vtxPosZ() const { return vtxPosZ_; }
  double sumEt() const { return sumEt_; }
  double met() const { return met_; }
  int nJets() const { return static_cast<int>(jets_.size()); }
  const Jet *jet(int idx) const { return idx >= 0 && idx < nJets() ? jets_[idx] : 0; }
  const Jet *corrJet(int idx) const { return idx >= 0 && idx < nJets() ? corrJets_[idx] : 0; }

 private:
  const unsigned int runNumber_;
  const unsigned int lumiBlockNumber_;
  const unsigned int vtxNTracks_;
  const double vtxPosZ_;
  const double sumEt_;
  const double met_;
  std::vector<Jet*> jets_;
  std::vector<Jet*> corrJets_;
};

Event::Event(unsigned int runNumber,
	     unsigned int lumiBlockNumber,
	     unsigned int vtxNTracks,
	     double vtxPosZ,
	     double sumEt,
	     double met,
	     const std::vector<Jet*>& jets)
  : runNumber_(runNumber),
    lumiBlockNumber_(lumiBlockNumber),
    vtxNTracks_(vtxNTracks),
    vtxPosZ_(vtxPosZ),
    sumEt_(sumEt),
    met_(met),
    jets_(jets),
    corrJets_(jets) {
  // Sort jets by pt
  std::sort(jets_.begin(),jets_.end(),Jet::ptGreaterThan);
  // Sort corrected jets by L2L3 corrected pt
  std::sort(corrJets_.begin(),corrJets_.end(),Jet::corrPtGreaterThan);  
}

Event::~Event() {
  std::vector<Jet*>::iterator jetIt = jets_.begin();
  for(; jetIt != jets_.end(); jetIt++) {
    delete *jetIt;
  }
  jets_.clear();
  corrJets_.clear();
}


// ----- Sample -----
class Sample {
public:
  Sample(const TString &sampleName, const TString &treeName);
  ~Sample();

  TString name() const { return name_; }
  TString drawOptions() const { return drawOpt_; }

  int nDists() const { return static_cast<int>(dists_.size()); }
  TString distName(int idx) const;
  void normaliseDists(const Sample *sample);
  double distIntegral(int idx) const { return isValidDistIdx(idx) ? dists_[idx]->Integral() : 0; }
  TH1D *dist(int idx) const { return isValidDistIdx(idx) ? dists_[idx]: 0; }

  void addFile(const TString &name);
  void fillDistributions();
  void readData(int nMaxEvts = -1);
  void setDrawOptions(const TString &opt) { drawOpt_ = opt; }

private:
  const TString name_;
  const int maxNJet_;

  TChain *chain_;
  Data data_;
  Distributions dists_;
  std::vector<unsigned int> runs_;
  TString drawOpt_;

  bool isValidDistIdx(int idx) const { return (idx >= 0 && idx < nDists() ); }
};

Sample::Sample(const TString &sampleName,const TString &treeName)
  : name_(sampleName), maxNJet_(50) {
  chain_ = new TChain(treeName);
  dists_ = Distributions(1);

  TString name = name_;
  name += ":Pt";
  dists_[0] = new TH1D(name,"Selected dijet events",15,0,60);
  dists_[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
  dists_[0]->GetYaxis()->SetTitle("Number of jets");

  for(int d = 0; d < nDists(); d++) {
    dists_[d]->SetMarkerStyle(20);
  }
}

Sample::~Sample() {
  delete chain_;
  for(DataIt it = data_.begin(); it != data_.end(); it++) {
    delete *it;
  }
  data_.clear();
  for(DistIt it = dists_.begin(); it != dists_.end(); it++) {
    delete *it;
  }
  dists_.clear();
  runs_.clear();
}

TString Sample::distName(int idx) const {
  TString name;
  if( isValidDistIdx(idx) ) {
    name = dists_[idx]->GetName();
    size_t pos = name.First(':');
    name = name(pos+1,name.Length()-pos);
  }
  return name;
}

void Sample::addFile(const TString &name) {
  std::cout << "Sample '" << name_ << "': Adding file '" << name << "... " << std::flush;
  chain_->AddFile(name);
  std::cout << "ok\n";
}

void Sample::fillDistributions() {
  std::cout << "Sample '" << name_ << "': Filling distributions... " << std::flush;

  // Reset distributions
  for(DistIt it = dists_.begin(); it != dists_.end(); it++) {
    (*it)->Reset();
  }

  // Loop over data
  for(DataIt it = data_.begin(); it != data_.end(); it++) {
    Event *evt = *it;

    // Loop over two jets leading in pt
    for(int j = 0; j < 2; j++) {
      dists_[0]->Fill(evt->jet(j)->pt());
// 	hPtCorr_[i]->Fill(jetCorrL2L3[j]*jetPt[j]);
// 	hEta_[i]->Fill(jetEta[j]);
// 	hPhi_[i]->Fill(jetPhi[j]);
    } // End of loop over dijets
  } // End of loop over data

      // Fill distributions of selected events
//       double diff = jetPt[0]-jetPt[1];
//       if( rand_->Uniform() > 0.5 ) diff = jetPt[1]-jetPt[0];
//       hPtAsym_[i]->Fill(diff/(jetPt[0]+jetPt[1]));
//       if( deltaPhi < 0 ) deltaPhi += 2*M_PI;
//       hDeltaPhi_[i]->Fill(deltaPhi);
  std::cout << "ok" << std::endl;
}

void Sample::normaliseDists(const Sample *sample) {
  for(int i = 0; i < nDists(); i++) {
    double norm = sample->distIntegral(i);
    if( distIntegral(i) ) {
      norm /= distIntegral(i);
      dists_[i]->Scale(norm);
    }
  }
}

void Sample::readData(int nMaxEvts) {
  std::cout << "Sample '" << name_ << "': Reading data... " << std::flush;

    // Reset data
    for(DataIt it = data_.begin(); it != data_.end(); it++) {
      delete *it;
    }
    data_.clear();

    // Init read quantities
    unsigned int runNumber = 0;
    unsigned int lumiBlockNumber = 0;
    int vtxNTracks = 0;
    float vtxPosZ = 0.;
    float sumEt = 0.;
    float met = 0.;
    int nObjJet = 0;
    float jetPt[maxNJet_];
    float jetEta[maxNJet_];
    float jetPhi[maxNJet_];
    float jetEMF[maxNJet_];
    int jetN90Hits[maxNJet_];
    float jetFHPD[maxNJet_];
    float jetFRBX[maxNJet_];
    float jetCorrL2L3[maxNJet_];

    // Set branch addresses
    chain_->SetBranchAddress("RunNumber",&runNumber);
    chain_->SetBranchAddress("LumiBlockNumber",&lumiBlockNumber);
    chain_->SetBranchAddress("VtxNTracks",&vtxNTracks);
    chain_->SetBranchAddress("VtxPosZ",&vtxPosZ);
    chain_->SetBranchAddress("Met",&met);
    chain_->SetBranchAddress("MetSum",&sumEt);
    chain_->SetBranchAddress("NobjJet",&nObjJet);
    chain_->SetBranchAddress("JetPt",jetPt);
    chain_->SetBranchAddress("JetEta",jetEta);
    chain_->SetBranchAddress("JetPhi",jetPhi);
    chain_->SetBranchAddress("JetEMF",jetEMF);
    chain_->SetBranchAddress("JetN90Hits",jetN90Hits);
    chain_->SetBranchAddress("JetFHPD",jetFHPD);
    chain_->SetBranchAddress("JetFRBX",jetFRBX);
    chain_->SetBranchAddress("JetCorrL2L3",jetCorrL2L3);
  
    // Loop over tree entries
    int nEntries = chain_->GetEntries();
    if( nMaxEvts > 0 && nEntries > nMaxEvts ) nEntries = nMaxEvts;
    for(int n = 0; n < nEntries; n++) {
      chain_->GetEntry(n);

      if( nObjJet > maxNJet_ ) {
	std::cerr << "WARNING: nObjJet = " << nObjJet << " > maxNJet_. Skipping event.\n";
	continue;
      }

      if( nObjJet < 2 ) continue;

      // Create event
      std::vector<Jet*> jets(nObjJet);
      for(int i = 0; i < nObjJet; i++) {
	jets[i] = new Jet(jetPt[i],jetEta[i],jetPhi[i],jetEMF[i],jetN90Hits[i],
			  jetFHPD[i],jetFRBX[i],jetCorrL2L3[i]);
      }
      Event *evt = new Event(runNumber,lumiBlockNumber,vtxNTracks,vtxPosZ,sumEt,met,jets);
      data_.push_back(evt);

      // Run number
      bool isNewRun = true;
      for(std::vector<unsigned int>::const_iterator r = runs_.begin();
	  r != runs_.end(); r++) {
	if( runNumber == *r ) {
	  isNewRun = false;
	  break;
	}
      }
      if( isNewRun ) {
	runs_.push_back(runNumber);
      }
    } // End of loop over tree entries
    std::sort(runs_.begin(),runs_.end());
    std::cout << "ok\n";
}




// ===== Main script =====
// =======================

void draw(std::vector<Sample*> samples, const TString &prefix) {
  std::cout << "Drawing distributions... " << std::flush;
  // Normalize distributions to sample 0
  for(size_t i = 1; i < samples.size(); i++) {
    samples[i]->normaliseDists(samples[0]);
  }

  // Draw and save distributions
  for(int d = 0; d < samples[0]->nDists(); d++) {
    std::cout << d << " " << std::flush;
    TString canName = "can";
    canName += samples[0]->distName(d);
    TCanvas *can = new TCanvas(canName,samples[0]->distName(d),500,500);
    can->cd();
    for(size_t s = 0; s < samples.size(); s++) {
      TString opt = samples[s]->drawOptions();
      if( s > 0 ) opt += "same";
      TH1D *h = samples[s]->dist(d);
      if( h ) h->Draw(opt);
    }
    TString fileName = prefix;
    fileName += "_";
    fileName += samples[0]->distName(d);
    fileName += ".eps";
    can->SaveAs(fileName,"eps");
  }
  std::cout << "ok" << std::endl;
}

void run() {
  int nSamples = 2;
  std::vector<Sample*> samples(nSamples);

  // The data sample
  samples[0] = new Sample("Data","DiJetTree");
  samples[0]->addFile("/Users/matthias/CMS/MinBias-BeamCommissioning09-2360GeV-Dec19thReReco_336p3_v2.root");
  samples[0]->setDrawOptions("P");

  // The MC sample
  samples[1] = new Sample("MC","DiJetTree");
  samples[1]->addFile("/Users/matthias/CMS/MinBias-BeamCommissioning09-2360GeV-Dec19thReReco_336p3_v2.root");

  for(int s = 0; s < nSamples; s++) {
    samples[s]->readData();
    samples[s]->fillDistributions();
  }

  draw(samples,"plots/test");
}

#ifndef JETRESOLUTION_H
#define JETRESOLUTION_H


#include <vector>
#include <string>

#include "TH1F.h"
#include "TF1.h"

//tool to provide jet resolutions
class JetResolution{
 public:
  JetResolution();
  virtual ~JetResolution(){};
  
  //Get a random number distributed according to the full non-Gaussian jet energy
  //resolution for a given jet pt, eta, and rank w.r.t. pt in the event.
  float GetRandom(float pt, float eta, int i_th);

  //Sigma^2 of the Gaussian Core jet energy resolution for a given jet w.r.t. pt, eta and phi
  double JetResolution_Pt2(const double& pt, const double& eta, const int& i_th);
  double JetResolution_Ptrel(const double& pt, const double& eta, const int& i_th);
  double JetResolution_Eta2(const double& energy, const double& eta);
  double JetResolution_Phi2(const double& energy, const double& eta);

 private:
  
  //Copied from SmearProducer:
  void LoadResponseHistograms();
  double GetAdditionalSmearing(const double&, const double&);
  double GetLowerTailScaling(const double&, const double&);
  double GetUpperTailScaling(const double&, const double&);
  int GetIndex(const double&, const std::vector<double>*);
  void FoldWithGaussian(const TH1&, TH1&, const double&);
  void StretchHisto(const TH1&, TH1&, const double&);

  std::vector<double> PtBinEdges_scaling_;
  std::vector<double> EtaBinEdges_scaling_;
  std::vector<double> AdditionalSmearing_;
  std::vector<double> LowerTailScaling_;
  std::vector<double> UpperTailScaling_;
  double AdditionalSmearing_variation_;
  double LowerTailScaling_variation_;
  double UpperTailScaling_variation_; 
  std::vector<std::vector<std::vector<TH1F*> > > smearFunc;
  std::vector<std::vector<std::vector<TH1F*> > > smearFunc_Core;
  std::vector<std::vector<std::vector<TH1F*> > > smearFunc_LowerTail;
  std::vector<std::vector<std::vector<TH1F*> > > smearFunc_UpperTail;
  std::vector<std::vector<std::vector<TH1F*> > > smearFunc_scaled;
  std::vector<std::vector<TH1F*> > SigmaPtHist;
  std::vector<std::vector<TF1*> > SigmaPt;
  std::vector<std::vector<TH1F*> > SigmaPtHist_scaled;
  std::vector<std::vector<TF1*> > SigmaPt_scaled;
  std::vector<double> PtBinEdges_;
  std::vector<double> EtaBinEdges_;
  std::string inputhist1_;
  std::string inputhist2_;
  std::string inputhist3p_;
  std::string smearingfile_;
  std::string outputfile_;
  int NRebin_;
  double A0RMS_;
  double A1RMS_;
  double probExtreme_;
  bool absoluteTailScaling_;

};



#endif

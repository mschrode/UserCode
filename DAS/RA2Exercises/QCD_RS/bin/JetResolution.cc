#include "JetResolution.h"

#include <cmath>
#include <iostream>

#include "TFile.h"
#include "TMath.h"

JetResolution::JetResolution()
{
   PtBinEdges_scaling_.push_back(0);PtBinEdges_scaling_.push_back(7000); 
   EtaBinEdges_scaling_.push_back(0);EtaBinEdges_scaling_.push_back(5);
   AdditionalSmearing_.push_back(0.5);
   LowerTailScaling_.push_back(1.0);
   UpperTailScaling_.push_back(1.0);
   AdditionalSmearing_variation_= 1.0;
   LowerTailScaling_variation_ 	= 1.0;
   UpperTailScaling_variation_ 	= 1.0;
   inputhist1_ 			= "h_tot_Jet1_ResponsePt";
   inputhist2_ 			= "h_tot_Jet2_ResponsePt";
   inputhist3p_ 		= "h_tot_Jet3p_ResponsePt";
   smearingfile_ 		= "JetResolution.root";
   outputfile_ 			= "SmearedAndScaledResponse";
   NRebin_ 			= 1;
   absoluteTailScaling_ 	= true;
   A0RMS_ 			= 2.5;
   A1RMS_ 			= 10.;
   probExtreme_ 		= 0.;
   PtBinEdges_.push_back(0); 
   PtBinEdges_.push_back(20); 
   PtBinEdges_.push_back(30); 
   PtBinEdges_.push_back(50); 
   PtBinEdges_.push_back(80); 
   PtBinEdges_.push_back(120); 
   PtBinEdges_.push_back(170); 
   PtBinEdges_.push_back(230); 
   PtBinEdges_.push_back(300); 
   PtBinEdges_.push_back(380); 
   PtBinEdges_.push_back(470); 
   PtBinEdges_.push_back(570); 
   PtBinEdges_.push_back(680); 
   PtBinEdges_.push_back(800); 
   PtBinEdges_.push_back(1000); 
   PtBinEdges_.push_back(1300); 
   PtBinEdges_.push_back(1700); 
   PtBinEdges_.push_back(2200); 
   PtBinEdges_.push_back(2800); 
   PtBinEdges_.push_back(3500);
   EtaBinEdges_.push_back(0); 
   EtaBinEdges_.push_back(0.3); 
   EtaBinEdges_.push_back(0.5); 
   EtaBinEdges_.push_back(0.8); 
   EtaBinEdges_.push_back(1.1); 
   EtaBinEdges_.push_back(1.4); 
   EtaBinEdges_.push_back(1.7); 
   EtaBinEdges_.push_back(2.0); 
   EtaBinEdges_.push_back(2.3); 
   EtaBinEdges_.push_back(2.8); 
   EtaBinEdges_.push_back(3.2); 
   EtaBinEdges_.push_back(4.1); 
   EtaBinEdges_.push_back(5.0);

   LoadResponseHistograms();
}


//--------------------------------------------------------------------------
// pt resolution for KinFitter
double JetResolution::JetResolution_Pt2(const double& pt, const double& eta, const int& i) {
   int i_jet;
   i < 2 ? i_jet = i : i_jet = 2;
   int i_eta = GetIndex(eta, &EtaBinEdges_);
   //return pow(pt, 2) * pow(SigmaPt_scaled.at(i_jet).at(i_eta)->Eval(pt), 2);
   return pow(pt, 2) * pow(SigmaPt.at(i_jet).at(i_eta)->Eval(pt), 2);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// relative pt resolution for KinFitter
double JetResolution::JetResolution_Ptrel(const double& pt, const double& eta, const int& i) {
   int i_jet;
   i < 2 ? i_jet = i : i_jet = 2;
   int i_eta = GetIndex(eta, &EtaBinEdges_);
   return SigmaPt_scaled.at(i_jet).at(i_eta)->Eval(pt);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// eta resolution for KinFitter
double JetResolution::JetResolution_Eta2(const double& e, const double& eta) {
   //may be artifically reduced (no angular fit)
   return (pow(0.05 / TMath::Sqrt(e), 2) + pow(0.005, 2)) / 1.e6;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// phi resolution for KinFitter
double JetResolution::JetResolution_Phi2(const double& e, const double& eta) {
   //may be artifically reduced (no angular fit)
   return (pow(0.05 / TMath::Sqrt(e), 2) + pow(0.005, 2)) / 1.e6;
}
//--------------------------------------------------------------------------

float JetResolution::GetRandom(float pt, float eta, int i_th)
{
   if (i_th>2) i_th = 2;
   int i_Pt  = GetIndex(pt,  &PtBinEdges_);
   int i_eta = GetIndex(eta, &EtaBinEdges_);

   return smearFunc_scaled.at(i_th).at(i_eta).at(i_Pt)->GetRandom();
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////


void JetResolution::LoadResponseHistograms()
{

   //// open root file/tree and create SmearingFunction histo
   TFile *f1 = new TFile(smearingfile_.c_str(), "READ", "", 0);

   smearFunc.resize(3); //// three bins for jet rank
   for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc.begin(); it != smearFunc.end(); ++it) {
      it->resize(EtaBinEdges_.size() - 1);
      for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
         jt->resize(PtBinEdges_.size() - 1);
      }
   }

   smearFunc_Core.resize(3); //// three bins for jet rank
   for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_Core.begin(); it
         != smearFunc_Core.end(); ++it) {
      it->resize(EtaBinEdges_.size() - 1);
      for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
         jt->resize(PtBinEdges_.size() - 1);
      }
   }

   smearFunc_LowerTail.resize(3); //// three bins for jet rank
   for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_LowerTail.begin(); it
         != smearFunc_LowerTail.end(); ++it) {
      it->resize(EtaBinEdges_.size() - 1);
      for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
         jt->resize(PtBinEdges_.size() - 1);
      }
   }

   smearFunc_UpperTail.resize(3); //// three bins for jet rank
   for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_UpperTail.begin(); it
         != smearFunc_UpperTail.end(); ++it) {
      it->resize(EtaBinEdges_.size() - 1);
      for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
         jt->resize(PtBinEdges_.size() - 1);
      }
   }

   smearFunc_scaled.resize(3); //// three bins for jet rank
   for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_scaled.begin(); it
         != smearFunc_scaled.end(); ++it) {
      it->resize(EtaBinEdges_.size() - 1);
      for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
         jt->resize(PtBinEdges_.size() - 1);
      }
   }

   SigmaPtHist_scaled.resize(3);
   for (int i_jet = 0; i_jet < 3; ++i_jet) {
      SigmaPtHist_scaled.at(i_jet).resize(EtaBinEdges_.size() - 1);
      for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
         char hname[100];
         sprintf(hname, "SigmaPtHist_scaled_Eta%i_Jet%i", i_eta, i_jet + 1);
         SigmaPtHist_scaled.at(i_jet).at(i_eta) = new TH1F (hname, hname, PtBinEdges_.size() - 1,
               &(PtBinEdges_.at(0)));
         SigmaPtHist_scaled.at(i_jet).at(i_eta)->Sumw2();
      }
   }

   SigmaPtHist.resize(3);
   for (int i_jet = 0; i_jet < 3; ++i_jet) {
      SigmaPtHist.at(i_jet).resize(EtaBinEdges_.size() - 1);
      for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
         char hname[100];
         sprintf(hname, "SigmaPtHist_Eta%i_Jet%i", i_eta, i_jet + 1);
         SigmaPtHist.at(i_jet).at(i_eta) = new TH1F (hname, hname, PtBinEdges_.size() - 1, &(PtBinEdges_.at(0)));
         SigmaPtHist.at(i_jet).at(i_eta)->Sumw2();
      }
   }

   SigmaPt_scaled.resize(3);
   for (unsigned int i_jet = 0; i_jet < 3; ++i_jet) {
      SigmaPt_scaled.at(i_jet).resize(EtaBinEdges_.size() - 1);
   }

   SigmaPt.resize(3);
   for (unsigned int i_jet = 0; i_jet < 3; ++i_jet) {
      SigmaPt.at(i_jet).resize(EtaBinEdges_.size() - 1);
   }

   //// Fetch histos and fit gaussian core
   for (unsigned int i_Pt = 0; i_Pt < PtBinEdges_.size() - 1; ++i_Pt) {
      for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
         //// Get the histos
         char hname[100];
         sprintf(hname, "%s_Pt%i_Eta%i", inputhist1_.c_str(), i_Pt, i_eta);
         smearFunc.at(0).at(i_eta).at(i_Pt) = (TH1F*) f1->FindObjectAny(hname);
         sprintf(hname, "%s_Pt%i_Eta%i", inputhist2_.c_str(), i_Pt, i_eta);
         smearFunc.at(1).at(i_eta).at(i_Pt) = (TH1F*) f1->FindObjectAny(hname);
         sprintf(hname, "%s_Pt%i_Eta%i", inputhist3p_.c_str(), i_Pt, i_eta);
         smearFunc.at(2).at(i_eta).at(i_Pt) = (TH1F*) f1->FindObjectAny(hname);
         for (unsigned int i_jet = 0; i_jet < 3; ++i_jet) {
            smearFunc.at(i_jet).at(i_eta).at(i_Pt)->Rebin(NRebin_);
            if (probExtreme_ > 0) {
               double p = probExtreme_ * smearFunc.at(i_jet).at(i_eta).at(i_Pt)->Integral();
               smearFunc.at(i_jet).at(i_eta).at(i_Pt)->SetBinContent(1, p);
            }
            //// Get width of gaussian core
            if (smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetEntries() > 50) {
               double RMS = smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetRMS();
               double MEAN = smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetMean();
               TF1* fitfunction = new TF1("f", "gaus(0)", MEAN - 1 * RMS, MEAN + 1 * RMS);
               fitfunction->SetParameters(smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetMaximum(), MEAN, RMS);
               smearFunc.at(i_jet).at(i_eta).at(i_Pt)->Fit(fitfunction, "LLRQN");
               double Pt = SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetBinCenter(i_Pt);
               double eta = (EtaBinEdges_.at(i_eta) + EtaBinEdges_.at(i_eta + 1)) / 2;
               double f = GetAdditionalSmearing(Pt, eta);
               SigmaPtHist.at(i_jet).at(i_eta)->SetBinContent(i_Pt + 1, std::abs(fitfunction->GetParameter(2)));
               SigmaPtHist.at(i_jet).at(i_eta)->SetBinError(i_Pt + 1, fitfunction->GetParError(2));
               SigmaPtHist_scaled.at(i_jet).at(i_eta)->SetBinContent(i_Pt + 1, std::abs(fitfunction->GetParameter(2))
                     * f);
               SigmaPtHist_scaled.at(i_jet).at(i_eta)->SetBinError(i_Pt + 1, fitfunction->GetParError(2) * f);

               //// Split smearFunc in Core and Tail
               TH1F* hResponseFit = new TH1F(*smearFunc.at(i_jet).at(i_eta).at(i_Pt));
               hResponseFit->Reset();
               for (int i = 0; i < hResponseFit->GetNbinsX(); ++i) {
                  hResponseFit->SetBinContent(i, fitfunction->Eval(hResponseFit->GetBinCenter(i)));
               }

               //// Split lower tail
               smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_jet).at(i_eta).at(i_Pt));
               smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->Add(hResponseFit, -1.);
               for (int i = 0; i <= smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->GetNbinsX(); ++i) {
                  double tmp_v = smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinContent(i);
                  double tmp_e = smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinError(i);
                  //// by definition a tail has positive entries
                  if (tmp_v < 0) {
                     tmp_v = 0;
                     tmp_e = 0;
                  } else {
                     //// suppress everything except for low response tail
                     double scale = 1;
                     double x = smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinCenter(i);
                     if (x > MEAN - 1 * RMS)
                        scale = 0.;
                     tmp_v = scale * smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinContent(i);
                     tmp_e = scale * smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinError(i);
                  }
                  smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->SetBinContent(i, tmp_v);
                  smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->SetBinError(i, tmp_e);
               }

               //// Split upper tail
               smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_jet).at(i_eta).at(i_Pt));
               smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->Add(hResponseFit, -1.);
               for (int i = 0; i <= smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->GetNbinsX(); ++i) {
                  double tmp_v = smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinContent(i);
                  double tmp_e = smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinError(i);
                  //// by definition a tail has positive entries
                  if (tmp_v < 0) {
                     tmp_v = 0;
                     tmp_e = 0;
                  } else {
                     //// suppress everything except for low response tail
                     double scale = 1;
                     double x = smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinCenter(i);
                     if (x < MEAN + 1 * RMS)
                        scale = 0.;
                     tmp_v = scale * smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinContent(i);
                     tmp_e = scale * smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinError(i);
                  }
                  smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->SetBinContent(i, tmp_v);
                  smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->SetBinError(i, tmp_e);
               }

               smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_jet).at(i_eta).at(i_Pt));
               smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt)->Add(smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt), -1.);
               smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt)->Add(smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt), -1.);

            } else {
               //// Set core and tail if needed
               smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_jet).at(i_eta).at(i_Pt));
               smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_jet).at(i_eta).at(i_Pt));
               smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->Reset();
               smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_jet).at(i_eta).at(i_Pt));
               smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->Reset();
            }
         }
      }
   }

   //// Fit scaled gaussian sigma as function of pt
   for (unsigned int i_jet = 0; i_jet < 3; ++i_jet) {
      for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
         char fname[100];
         sprintf(fname, "SigmaPtScaled_Eta%i_Jet%i", i_eta, i_jet + 1);
         bool first = false;
         int FirstBin = 1;
         int LastBin = SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetNbinsX();
         for (int j = 1; j <= SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetNbinsX(); ++j) {
            if (!first && SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetBinContent(j) > 0) {
               first = true;
               FirstBin = j;
            }
            if (first && SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetBinContent(j) < 0.001) {
               LastBin = j - 1;
               break;
            }
         }
         SigmaPt_scaled.at(i_jet).at(i_eta) = new TF1(fname,
               "TMath::Sqrt(sign([0])*pow([0]/x,2)+pow([1],2)*pow(x,[2]-1.)+pow([3],2))",
               SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetBinCenter(FirstBin),
               SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetBinCenter(LastBin));
         SigmaPt_scaled.at(i_jet).at(i_eta)->SetParameters(1.2, 0., 0.03);
         SigmaPtHist_scaled.at(i_jet).at(i_eta)->Fit(SigmaPt_scaled.at(i_jet).at(i_eta), "LLRQ");
      }
   }

   //// Fit gaussian sigma as function of pt
   for (unsigned int i_jet = 0; i_jet < 3; ++i_jet) {
      for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
         char fname[100];
         sprintf(fname, "SigmaPt_Eta%i_Jet%i", i_eta, i_jet + 1);
         bool first = false;
         int FirstBin = 1;
         int LastBin = SigmaPtHist.at(i_jet).at(i_eta)->GetNbinsX();
         for (int j = 1; j <= SigmaPtHist.at(i_jet).at(i_eta)->GetNbinsX(); ++j) {
            if (!first && SigmaPtHist.at(i_jet).at(i_eta)->GetBinContent(j) > 0) {
               first = true;
               FirstBin = j;
            }
            if (first && SigmaPtHist.at(i_jet).at(i_eta)->GetBinContent(j) < 0.001) {
               LastBin = j - 1;
               break;
            }
         }
         SigmaPt.at(i_jet).at(i_eta) = new TF1(fname,
               "TMath::Sqrt(sign([0])*pow([0]/x,2)+pow([1],2)*pow(x,[2]-1.)+pow([3],2))", SigmaPtHist.at(i_jet).at(
                     i_eta)->GetBinCenter(FirstBin), SigmaPtHist.at(i_jet).at(i_eta)->GetBinCenter(LastBin));
         SigmaPt.at(i_jet).at(i_eta)->SetParameters(1.2, 0., 0.03);
         SigmaPtHist.at(i_jet).at(i_eta)->Fit(SigmaPt.at(i_jet).at(i_eta), "LLRQ");
      }
   }

   //// Book and fill histograms for smeared and scaled response functions
   for (unsigned int i_jet = 0; i_jet < 3; ++i_jet) {
      for (unsigned int i_Pt = 0; i_Pt < PtBinEdges_.size() - 1; ++i_Pt) {
         for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
            char hname[100];
            sprintf(hname, "SmearedAndScaledResolution_Pt%i_Eta%i_Jet%i", i_Pt, i_eta, i_jet + 1);
            smearFunc_scaled.at(i_jet).at(i_eta).at(i_Pt) = new TH1F (hname, hname,
                  smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetNbinsX(),
                  smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetXaxis()->GetXmin(),
                  smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetXaxis()->GetXmax());

            if (smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetEntries() > 50) {
               //// fold core and tail with additional gaussian
               TH1F smearFunc_Core_tmp(*smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt));
               TH1F smearFunc_LowerTail_tmp(*smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt));
               TH1F smearFunc_UpperTail_tmp(*smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt));
               double AddSmear = GetAdditionalSmearing(PtBinEdges_.at(i_Pt), EtaBinEdges_.at(i_eta));
               if (AddSmear > 1) {
                  //// additional sigma from (1+x)*sigma = sqrt(sigma^2+add_sigma^2)
                  //// or from sigma' = sqrt((sigma'/(1+x))^2+add_sigma^2)
                  double sigma = SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetBinContent(i_Pt + 1);
                  //// if no sigma was fitted use the extrapolation
                  if (sigma == 0)
                     sigma = SigmaPt_scaled.at(i_jet).at(i_eta)->Eval(
                           SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetBinCenter(i_Pt + 1));
                  double AdditionalSigma = TMath::Sqrt(1 - 1 / pow(AddSmear, 2)) * sigma;
                  smearFunc_Core_tmp.Reset();
                  smearFunc_LowerTail_tmp.Reset();
                  smearFunc_UpperTail_tmp.Reset();
                  FoldWithGaussian(*smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt), smearFunc_Core_tmp, AdditionalSigma);
                  FoldWithGaussian(*smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt), smearFunc_LowerTail_tmp,
                        AdditionalSigma);
                  FoldWithGaussian(*smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt), smearFunc_UpperTail_tmp,
                        AdditionalSigma);
               } else if (AddSmear < 1) {
                  smearFunc_Core_tmp.Reset();
                  StretchHisto(*smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt), smearFunc_Core_tmp, AddSmear);
               }

               //// Scale tails
               double LowerTailScale = GetLowerTailScaling(PtBinEdges_.at(i_Pt), EtaBinEdges_.at(i_eta));
               double UpperTailScale = GetUpperTailScaling(PtBinEdges_.at(i_Pt), EtaBinEdges_.at(i_eta));
               //cout << "absolute scaling factor: " << TailScale << endl;
               if (absoluteTailScaling_) {
                  double RMS = smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetRMS();
                  //cout << "Integral from " << 1-A1RMS_*RMS << " to " << 1-A0RMS_*RMS << endl;

                  //// get integral of tails
                  int i_min = smearFunc.at(i_jet).at(i_eta).at(i_Pt)->FindBin(1 - A1RMS_ * RMS);
                  int i_max = smearFunc.at(i_jet).at(i_eta).at(i_Pt)->FindBin(1 - A0RMS_ * RMS);
                  double RLowerTail = smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->Integral(i_min, i_max);
                  double Rcore = smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt)->Integral(i_min, i_max);
                  if (RLowerTail > 0)
                     LowerTailScale = (LowerTailScale * (RLowerTail + Rcore) - Rcore) / RLowerTail;

                  i_min = smearFunc.at(i_jet).at(i_eta).at(i_Pt)->FindBin(1 + A0RMS_ * RMS);
                  i_max = smearFunc.at(i_jet).at(i_eta).at(i_Pt)->FindBin(1 + A1RMS_ * RMS);
                  double RUpperTail = smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->Integral(i_min, i_max);
                  Rcore = smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt)->Integral(i_min, i_max);
                  if (RUpperTail > 0)
                     UpperTailScale = (UpperTailScale * (RUpperTail + Rcore) - Rcore) / RUpperTail;

               }
               //cout << "tail scaling factor: " << TailScale << endl;
               smearFunc_scaled.at(i_jet).at(i_eta).at(i_Pt)->Add(&smearFunc_Core_tmp);
               smearFunc_scaled.at(i_jet).at(i_eta).at(i_Pt)->Add(&smearFunc_LowerTail_tmp, LowerTailScale);
               smearFunc_scaled.at(i_jet).at(i_eta).at(i_Pt)->Add(&smearFunc_UpperTail_tmp, UpperTailScale);
            } else {
               //// Replace Histograms with only few entries by gaussians
               double N = 1;
               if (smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetEntries() > 0) {
                  N = smearFunc.at(i_jet).at(i_eta).at(i_Pt)->Integral();
               }
               //std::cout << "Too few entries for (i_Pt, i_eta, i_jet): " << i_Pt << ", " << i_eta << ", " << i_jet
               //      << ", entries = " << smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetEntries() << std::endl;
               for (int j = 1; j <= smearFunc_scaled.at(i_jet).at(i_eta).at(i_Pt)->GetNbinsX(); ++j) {
                  double pt = (PtBinEdges_.at(i_Pt) + PtBinEdges_.at(i_Pt + 1)) / 2;
                  double g = N * TMath::Gaus(smearFunc_scaled.at(i_jet).at(i_eta).at(i_Pt)->GetBinCenter(j), 1.,
                        SigmaPt_scaled.at(i_jet).at(i_eta)->Eval(pt));
                  smearFunc_scaled.at(i_jet).at(i_eta).at(i_Pt)->SetBinContent(j, g);
                  smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt)->SetBinContent(j, g);
                  smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->SetBinContent(j, 0.);
                  smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->SetBinContent(j, 0.);
               }
            }
         }
      }
   }
}

//--------------------------------------------------------------------------
double JetResolution::GetAdditionalSmearing(const double& pt, const double& eta) {
   int i_Pt = GetIndex(pt, &PtBinEdges_scaling_);
   int i_eta = GetIndex(eta, &EtaBinEdges_scaling_);
   double result = AdditionalSmearing_.at(i_eta * (PtBinEdges_scaling_.size() - 1) + i_Pt)
         * AdditionalSmearing_variation_;
   //if (result < 1) result = 1;
   return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double JetResolution::GetLowerTailScaling(const double& pt, const double& eta) {
   int i_Pt = GetIndex(pt, &PtBinEdges_scaling_);
   int i_eta = GetIndex(eta, &EtaBinEdges_scaling_);
   double result = LowerTailScaling_.at(i_eta * (PtBinEdges_scaling_.size() - 1) + i_Pt) * LowerTailScaling_variation_;
   return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double JetResolution::GetUpperTailScaling(const double& pt, const double& eta) {
   int i_Pt = GetIndex(pt, &PtBinEdges_scaling_);
   int i_eta = GetIndex(eta, &EtaBinEdges_scaling_);
   double result = UpperTailScaling_.at(i_eta * (PtBinEdges_scaling_.size() - 1) + i_Pt) * UpperTailScaling_variation_;
   return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Fold histogram input with gaussian resolution
void JetResolution::FoldWithGaussian(const TH1& input, TH1& output, const double& sigma) {

   double min = input.GetXaxis()->GetXmin();
   double max = input.GetXaxis()->GetXmax();
   for (int i = 0; i < input.GetNbinsX(); ++i) {
      double weight = input.GetBinContent(i);
      double mean = input.GetBinCenter(i);
      TF1 gauss("gauss", "gaus(0)", min, max);
      gauss.SetParameters(weight * 1 / sigma / sqrt(2 * TMath::Pi()), mean, sigma);
      for (int j = 0; j < output.GetNbinsX(); ++j) {
         double xmin = output.GetBinLowEdge(j);
         double xmax = output.GetBinLowEdge(j) + output.GetBinWidth(j);
         output.AddBinContent(j, gauss.Integral(xmin, xmax));
      }
   }

}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Fold histogram input with gaussian resolution
void JetResolution::StretchHisto(const TH1& input, TH1& output, const double& f) {

   if (input.Integral() > 0) {
      double mean = input.GetMean();
      for (int i = 0; i < 1000000; ++i) {
         double r = input.GetRandom();
         double rprime = mean + (r - mean) * f;
         output.Fill(rprime);
      }
      output.Scale(input.Integral() / output.Integral());
   }

}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
int JetResolution::GetIndex(const double& x, const std::vector<double>* vec) {
   int index = -1;
   // this is a check
   //int index = 0;
   for (std::vector<double>::const_iterator it = vec->begin(); it != vec->end(); ++it) {
      if ((*it) > fabs(x))
         break;
      ++index;
   }
   if (index < 0)
      index = 0;
   if (index > (int) vec->size() - 2)
      index = vec->size() - 2;

   return index;
}
//--------------------------------------------------------------------------

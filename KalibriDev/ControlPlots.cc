#include "ControlPlots.h"

#include <algorithm>
#include <iostream>
using namespace std;

#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1I.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TPostScript.h"
#include "TString.h"
#include "TROOT.h"

#include "CalibMath.h"

#include "JetTruthEvent.h"
#include "TwoJetsInvMassEvent.h"
#include "TwoJetsPtBalanceEvent.h"


//!  \brief Constructor
//! 
//!  \param configfile  Path to config file
//!  \param data        TData objects to create controlplots from
//!  \param par         Parameters
// -------------------------------------------------------------
TControlPlots::TControlPlots(const std::string& configfile, const std::vector<TData*> *data, TParameters *par)
  : config_(new ConfigFile(configfile.c_str())), data_(data), outFile_(0), par_(par)
{ 
  // Open ROOT file for writing
  if(config_->read<bool>("plot output format",0)) {
    outputROOT_ = false;
  } else {
    outputROOT_ = true;
    outFile_ = new TFile("controlplots.root","RECREATE","Kalibri control plots");
  }

  // Copy fitted parameter values in case the values in
  // par_->GetPars() are changed to do plots with different
  // constants
  fittedPar_ = std::vector<double>(par_->GetNumberOfParameters());
  for(int i = 0; i < par_->GetNumberOfParameters(); i++) {
    fittedPar_.at(i) = par_->GetPars()[i];
  }

  setGStyle();
}



// -------------------------------------------------------------
TControlPlots::~TControlPlots()
{
  if( outFile_ !=0 )
    {
      if( outFile_->IsOpen() ) outFile_->Close();
      delete outFile_;
    }
}



//!  \brief Create all plots as specified in the config file
// -------------------------------------------------------------
void TControlPlots::makePlots()
{
  cout << endl << "Writing control plots in .ps " << flush;
  if( outputROOT_ ) cout << "and .root " << flush;
  cout << "format:" << endl;
  if( config_->read<bool>("create binned response plots",false) ) {
    cout << "Creating binned control plots\n";
    makeControlPlotsBinnedResponse();
  }
  if( config_->read<bool>("create chi2 plots",false) ) {
    cout << "Creating chi2 control plots\n";
    makeControlPlotsChi2();
  }
  if( config_->read<bool>("create jet-truth event plots",false) ) {
    cout << "Creating jet-truth event control plots\n";
    makeControlPlotsJetTruthEventResponse();
  }
  if( config_->read<bool>("create L2L3 MC truth plots",false) ) {
    cout << "Creating L2L3 MC truth control plots\n";
    makeControlPlotsL2L3MCTruth();
  }
  if( config_->read<bool>("create parameter scan plots" ,false) ) {
    cout << "Creating parameter scan control plots\n";
    makeControlPlotsParameterScan();
  }
  if( config_->read<bool>("create top plots",false) ) {
    cout << "Creating top control plots\n";
    makeControlPlotsTop();
  }
  if( config_->read<bool>("create dijets pt-balance plots",false) ) {
    cout << "Creating dijet pt-balance control plots\n";
    makeControlPlotsTwoJetsPtBalance();
  }
}



//!  \brief Response distribution in bins of truth Et and eta
//!
//!  Creates distributions of the response
//!  \f$ E^{jet}_{T} / E^{true}_{T} \f$
//!  in bins of true Et \f$ E^{true}_{T} \f$ and measured
//!  pseudorapidity \f$ \eta \f$, where \f$ E^{jet}_{T} \f$
//!  is the corrected jet Et.
//!
//!  There are two sets of distributions:
//!   - Correction from global fit
//!   - Comparison of correction from global fit
//!     and JetMET L2L3 correction
//!
//!  The plots are written to the file
//!  "controlplotsBinnedResponse.ps" and (if enabled)
//!  to the directory "BinnedResponse" in the file
//!  "controlplots.root".
// -------------------------------------------------------------
void TControlPlots::makeControlPlotsBinnedResponse()
{
  std::vector<TObject*>                objToBeWritten;
  std::vector<TData*>::const_iterator  data_it;


  // Create a Et and eta binning
  // Et = x, eta = y
  std::vector<double> binEdgesEt
    = bag_of<double>(config_->read<std::string>("Control plots pt bin edges","0 100"));
  std::vector<double> binEdgesEta
    = bag_of<double>(config_->read<std::string>("Control plots eta bin edges","-5 5"));
  Binning bins(binEdgesEt,binEdgesEta);


  // Create histograms of response in
  // Et and eta bins
  std::vector<TH1F*> hResp;
  std::vector<TH1F*> hRespJetMet;
  for(int i = 0; i < bins.nBins(); i++)
    {
      char name[50];
      sprintf(name,"hResponse_pt%i_eta%i",bins.iX(i),bins.iY(i));
      TH1F * h = new TH1F(name,"",51,0,2);
      h->SetLineColor(2);
      h->GetXaxis()->SetTitleSize(0.05);
      h->GetXaxis()->SetTitleOffset(1.2);
      h->GetYaxis()->SetTitleSize(0.05);
      h->GetYaxis()->SetTitleOffset(1.6);
      hResp.push_back(h);
      objToBeWritten.push_back(h);

      sprintf(name,"hResponseJetMET_pt%i_eta%i",bins.iX(i),bins.iY(i));
      h = static_cast<TH1F*>(h->Clone(name));
      h->SetLineColor(1);
      hRespJetMet.push_back(h);
      objToBeWritten.push_back(h);
    }


  // Loop over data and fill response into
  // the corresponding histogram
  for(  data_it = data_->begin();  data_it != data_->end();  data_it++ )
    {
      double weight = (*data_it)->GetWeight();
      double etmeas = (*data_it)->GetMess()->pt;
      double etcorr = (*data_it)->GetParametrizedMess();
      double eta    = (*data_it)->GetMess()->eta;
      double ettrue = (*data_it)->GetTruth();

      // Discard events flagged bad
      JetTruthEvent *jte = dynamic_cast<JetTruthEvent*>(*data_it);
      if( jte )
	{
	  if( jte->FlaggedBad() ) continue;
	}

      int bin = bins.bin(ettrue,eta);
      if( 0 <= bin && bin < bins.nBins() )
	{
	  hResp.at(bin)->Fill(etcorr/ettrue,weight);

	  // For comparison reasons, correct with JetMET
	  // correction; this works only for data class > 0
	  TJet *jet = dynamic_cast<TJet*>((*data_it)->GetMess());
	  if( jet )
	    {
	      double cjetmet = jet->corFactors.getL2L3();
	      hRespJetMet.at(bin)->Fill(cjetmet*etmeas/ettrue,weight);
	    }
	}
    }


  // Draw histograms into ps file, put 
  // at the most 6 histograms per page
  TPostScript * const ps = new TPostScript("controlplotsBinnedResponse.ps",112);
  ps->Range(25,1);

  double cw  = 300;
  double mw  = 70;

  double w   = 3*cw + mw;
  double h   = 2*cw + mw;
  
  double crw = cw/w;
  double mrw = mw/w;
  double crh = cw/h;
  double mrh = mw/h;

  TCanvas     * const c1 = new TCanvas("c1","",(int)w,(int)h);

  // Pads for the histograms
  std::vector<TPad*> cPads;
  cPads.push_back(new TPad("cPad0","",mrw,mrh+crh,mrw+crw,1.));
  cPads.push_back(new TPad("cPad1","",mrw+crw,mrh+crh,mrw+2*crw,1.));
  cPads.push_back(new TPad("cPad2","",mrw+2*crw,mrh+crh,1.,1.));
  cPads.push_back(new TPad("cPad3","",mrw,mrh,mrw+crw,mrh+crh));
  cPads.push_back(new TPad("cPad4","",mrw+crw,mrh,mrw+2*crw,mrh+crh));
  cPads.push_back(new TPad("cPad5","",mrw+2*crw,mrh,1.,mrh+crh));
  for(unsigned int i = 0; i < cPads.size(); i++)
    {
      cPads.at(i)->SetFillStyle(1001);
      cPads.at(i)->SetFrameFillColor(10);
      cPads.at(i)->SetFrameBorderMode(0);
      cPads.at(i)->SetTopMargin(0.04);
      cPads.at(i)->SetBottomMargin(0.1);
      cPads.at(i)->SetLeftMargin(0.1);
      cPads.at(i)->SetRightMargin(0.04);

      c1->cd();
      cPads.at(i)->Draw();
    }

  // Pads for margins holding titles
  TPad * mbPad = new TPad("mbPad","",0.,0.,1.,mrh);
  mbPad->SetFillStyle(1001);
  mbPad->SetFrameFillColor(10);
  mbPad->SetFrameBorderMode(0);
  TPaveText * bLabel = new TPaveText(0.8,0.5,0.95,1.,"NDC");
  bLabel->SetFillColor(0);
  bLabel->SetTextFont(42);
  bLabel->SetTextSize(0.7);
  bLabel->SetBorderSize(0);
  bLabel->AddText("p^{jet}_{T} / p^{true}_{T}");
  c1->cd();
  mbPad->Draw();
  mbPad->cd();
  bLabel->Draw();

  //   TPad * mlPad = new TPad("mlPad","",0.,mrh,mrw,1.);
  //   mlPad->SetFillStyle(1001);
  //   mlPad->SetFrameFillColor(10);
  //   mlPad->SetFrameBorderMode(0);
  //   TPaveText * lLabel = new TPaveText(0.1,0.5,0.95,1.,"NDC");
  //   lLabel->SetFillColor(0);
  //   lLabel->SetTextFont(42);
  //   lLabel->SetTextSize(0.7);
  //   lLabel->SetBorderSize(0);
  //   lLabel->SetTextAngle(90);
  //   lLabel->AddText("dN / d( E^{jet}_{T} / E^{true}_{T} )");
  //   c1->cd();
  //   mlPad->Draw();
  //   mlPad->cd();
  //   lLabel->Draw();

  // For some reason, this prevents the first ps-page
  // from looking weird...
  c1->Draw();
  ps->NewPage();
  for(unsigned int i = 0; i < hResp.size(); i++)
    {
      // Plot histogram
      int c = ( i % 6 );
      cPads.at(c)->cd();
      TH1F *h = hResp.at(i);
      h->GetYaxis()->SetRangeUser(0,1.6*(h->GetMaximum()));
      h->Draw();

      // Fit a gaussian in central part
      h->Fit("gaus","0QIL","",h->GetMean() - 1.5*(h->GetRMS()),h->GetMean() + 1.5*(h->GetRMS()));
      TF1 *fit = h->GetFunction("gaus");
      fit->Draw("same");

      // Draw fit parameters
      char label[100];
      TPaveText * fitstat = new TPaveText(0.13,0.66,0.93,0.93,"NDC");
      fitstat->SetFillColor(0);
      fitstat->SetTextFont(42);
      fitstat->SetTextAlign(12);
      sprintf(label,"%.0f < p^{true}_{T} < %.0f GeV, %.1f < #eta < %.1f",
	      bins.xLow(i),bins.xUp(i),bins.yLow(i),bins.yUp(i));
      fitstat->AddText(label);
      sprintf(label,"#mu = %.3f #pm %.3f",fit->GetParameter(1),fit->GetParError(1));
      fitstat->AddText(label);
      sprintf(label,"#sigma = %.3f #pm %.3f",fit->GetParameter(2),fit->GetParError(2));
      fitstat->AddText(label);
      fitstat->Draw("same");

      if( c == 5 )
	{
	  c1->Draw();
	  ps->NewPage();
	}
    }

  // Now the comparison plots to JetMET 
  TLegend *legComp = new TLegend(0.13,0.66,0.6,0.84);
  legComp->SetBorderSize(0);
  legComp->SetFillColor(0);
  legComp->SetTextFont(42);
  legComp->AddEntry(hResp.at(0),"Global fit","L");
  legComp->AddEntry(hRespJetMet.at(0),"JetMET L2L3","L");
  for(unsigned int i = 0; i < hResp.size(); i++)
    {
      int c = ( i % 6 );
      c1->cd();
      cPads.at(c)->cd();

      // Set maximum from histo with larger bin content
      double max = (hResp.at(i)->GetMaximum())/1.6; // Compensate for already setting this maximum above!
      if( hRespJetMet.at(i)->GetMaximum() > max ) max = hRespJetMet.at(i)->GetMaximum();
      hResp.at(i)->GetYaxis()->SetRangeUser(0,1.6*max);
      hResp.at(i)->Draw();
      hRespJetMet.at(i)->Draw("same");

      // Draw a vertical line at 1
      TLine *line = new TLine(1.,0.,1.,max);
      line->SetLineWidth(1);
      line->SetLineColor(4);
      line->SetLineStyle(2);
      line->Draw("same");

      // Label bin
      char label[100];
      TPaveText * fitstat = new TPaveText(0.13,0.84,0.93,0.93,"NDC");
      fitstat->SetFillColor(0);
      fitstat->SetTextFont(42);
      fitstat->SetTextAlign(12);
      sprintf(label,"%.0f < p^{true}_{T} < %.0f GeV, %.1f < #eta < %.1f",
	      bins.xLow(i),bins.xUp(i),bins.yLow(i),bins.yUp(i));
      fitstat->AddText(label);
      fitstat->Draw("same");

      if( c == 5 )
	{
	  legComp->Draw("same");
	  c1->Draw();
	  ps->NewPage();
	}
    }

  if( outputROOT_ ) writeToRootFile(objToBeWritten,"BinnedResponse");


  // Clean up
  ps->Close();
  for(std::vector<TH1F*>::iterator it = hResp.begin();
      it != hResp.end(); it++)
    {
      delete *it;
    }
  hResp.clear();
  for(std::vector<TH1F*>::iterator it = hRespJetMet.begin();
      it != hRespJetMet.end(); it++)
    {
      delete *it;
    }
  hRespJetMet.clear();
  objToBeWritten.clear();
  delete c1;
  delete ps;
}



//!  \brief Chi2 distribution
//!
//!  Creates the following distributions of
//!  \f$ \chi^{2}_{i} \f$ from the sum
//!  \f[
//!   \chi^{2} = \sum_{i} w_{i}\chi^{2}_{i}
//!  \f]
//!  where \f$ w_{i} \f$ is the event weight
//!  factor:
//!   - The residual scaling scheme from the last 
//!     iteration is used.
//!   - All three residual scaling schemes are
//!     compared:
//!      - Distributions of scaled residuals
//!      - Scaled residuals vs residuals
//!   - The \f$ \chi^{2}_{i} \f$ probability
//!     prob(chi2,1).
//!
//!  The plots are written to the file
//!  "controlplotsChi2.ps" and (if enabled)
//!  to the directory "Chi2" in the file
//!  "controlplots.root".
// -------------------------------------------------------------
void TControlPlots::makeControlPlotsChi2()
{
  std::vector<TObject*>                objToBeWritten;
  std::vector<TData*>::const_iterator  data_it;

  // Distribution of chi2 summands with last scaling
  // configuration scheme
  TH1F *h_chi2_last = new TH1F("h_chi2last",";#chi^{2}_{i};dN / d#chi^{2}_{i}",50,0,20);
  h_chi2_last->SetLineColor(1);
  h_chi2_last->SetLineStyle(1);
  objToBeWritten.push_back(h_chi2_last);

  // Distribution of chi2 probability with last scaling
  // configuration scheme
  TH1F *h_prob = new TH1F("h_prob",";prob(#chi^{2}_{i},1)",50,0,1);
  h_prob->SetMarkerStyle(7);
  objToBeWritten.push_back(h_prob);

  // Distribution of chi2 summands (normalized residuals)
  TH1F *h_chi2 = new TH1F("h_chi2","Scaled residuals z^{2} = 1/w #chi^{2};f(z^{2});dN / df(z^{2})",50,0,20);
  h_chi2->SetLineColor(1);
  h_chi2->SetLineStyle(1);
  objToBeWritten.push_back(h_chi2);

  // Distribution of Cauchy scaled normalized residuals
  TH1F *h_cauchy = static_cast<TH1F*>(h_chi2->Clone("h_cauchy"));
  h_cauchy->SetLineColor(2);
  h_cauchy->SetLineStyle(2);
  objToBeWritten.push_back(h_cauchy);

  // Distribution of Huber scaled normalized residuals
  TH1F *h_huber = static_cast<TH1F*>(h_cauchy->Clone("h_huber"));
  h_huber->SetLineColor(4);
  h_huber->SetLineStyle(1);
  objToBeWritten.push_back(h_huber);

  // Cauchy-scaled versus no scaling
  TH2F *h_none_cauchy = new TH2F("h_none_cauchy","Residuals z^{2};z^{2};f(z^{2})",50,0,10,50,0,10);
  h_none_cauchy->SetLineColor(2);
  objToBeWritten.push_back(h_none_cauchy);

  // Huber-scaled versus no scaling
  TH2F *h_none_huber = static_cast<TH2F*>(h_none_cauchy->Clone("h_none_huber"));
  h_none_huber->SetLineColor(4);
  objToBeWritten.push_back(h_none_huber);


  // Loop over data and fill histograms
  for(  data_it = data_->begin();  data_it != data_->end();  data_it++ )
    {
      double weight   = (*data_it)->GetWeight();
      double lastchi2 = ((*data_it)->chi2_plots())/weight;

      h_chi2_last->Fill(lastchi2);
      h_prob->Fill(TMath::Prob(lastchi2,1));

      TData::ScaleResidual = &TData::ScaleNone;
      double res = ( (*data_it)->chi2() / weight );
      TData::ScaleResidual = &TData::ScaleCauchy;
      double res_cauchy = ( (*data_it)->chi2() / weight );
      TData::ScaleResidual = &TData::ScaleHuber;
      double res_huber = ( (*data_it)->chi2() / weight );

      h_chi2->Fill(res);
      h_cauchy->Fill(res_cauchy);
      h_huber->Fill(res_huber);
      h_none_cauchy->Fill(res,res_cauchy);
      h_none_huber->Fill(res,res_huber);
    } // End loop over data


  // Draw histograms
  TPostScript * const ps = new TPostScript("controlplotsChi2.ps",111);
  TCanvas     * const c1 = new TCanvas("c1","",500,500);
  c1->cd();

  h_chi2_last->Draw();
  c1->SetLogy(1);
  c1->Draw();
  ps->NewPage();

  h_prob->GetYaxis()->SetRangeUser(0.,1.2*(h_prob->GetMaximum()));
  h_prob->Draw("PE1");
  c1->SetLogy(0);
  c1->Draw();
  ps->NewPage();

  h_cauchy->Draw();
  h_huber->Draw("same");
  h_chi2->Draw("same");
  c1->SetLogy(1);

  TLegend *l_res = new TLegend(0.35,0.68,0.7,0.88);
  l_res->SetFillColor(0);
  l_res->SetBorderSize(0);
  l_res->SetTextFont(42);
  l_res->SetHeader("Scaling function f");
  l_res->AddEntry(h_chi2,"None","L");
  l_res->AddEntry(h_huber,"Huber","L");
  l_res->AddEntry(h_cauchy,"Cauchy","L");
  l_res->Draw("same");

  c1->Draw();
  ps->NewPage();


  h_none_cauchy->Draw("box");
  h_none_huber->Draw("boxsame");
  c1->SetLogy(0);

  TLine *line1 = new TLine(0,0,7,7);
  line1->SetLineStyle(2);
  line1->SetLineColor(1);
  line1->SetLineWidth(1);
  line1->Draw("same");

  TLegend *l_res2 = new TLegend(0.35,0.68,0.7,0.88);
  l_res2->SetFillColor(0);
  l_res2->SetBorderSize(0);
  l_res2->SetTextFont(42);
  l_res2->SetHeader("Scaling function f");
  l_res2->AddEntry(line1,"None","L");
  l_res2->AddEntry(h_none_huber,"Huber","L");
  l_res2->AddEntry(h_none_cauchy,"Cauchy","L");
  l_res2->Draw("same");

  c1->Draw();

  if( outputROOT_ ) writeToRootFile(objToBeWritten, "Chi2");
  
  ps->Close();


  // Clean up
  delete l_res;
  delete line1;
  delete l_res2;
  delete c1;
  delete ps;
  for(std::vector<TObject*>::iterator it = objToBeWritten.begin();
      it != objToBeWritten.end(); it++)
    {
      delete *it;
    }
  objToBeWritten.clear();
}



//!  \brief Response and resolution distribution in bins of
//!         truth Pt and eta
//!
//!  The plots are written to the file
//!  "controlplotsJetTruthEventResponse.ps",
//!  "controlplotsJetTruthEventResolution.ps", and (if enabled)
//!  to the directory "JetTruthEventResponse" in the file
//!  "controlplots.root".
// -------------------------------------------------------------
void TControlPlots::makeControlPlotsJetTruthEventResponse() {
  std::vector<TObject*>                objToBeWritten;
  std::vector<TData*>::const_iterator  data_it;

  // -- Create 2D histograms of response vs eta, pt --------------

  // Create a pttrue and eta binning
  // pttrue = x, eta = y
  std::vector<double> binEdgesPt = bag_of<double>(config_->read<std::string>("Control plots pt bin edges",""));

  // Absolute eta bins
  std::vector<double> binEdgesEta = bag_of<double>(config_->read<std::string>("Control plots eta bin edges",""));
  Binning bins(binEdgesPt,binEdgesEta);
  //  bins.Print();


  // Create 2D histograms of response vs eta in
  // Pt and eta bins
  std::vector<TH2F*> h2EtaUncorr;    // Uncorrected response vs eta
  std::vector<TH2F*> h2EtaCorr;      // Response corrected by kalibri fit vs eta
  std::vector<TH2F*> h2EtaCorrL2L3;  // Response corrected by JetMET L2L3 correction vs eta
  for(int ptbin = 0; ptbin < bins.nBinsX(); ptbin++) { // Loop over pttrue bins
    char name[50];
    sprintf(name,"h2EtaUncorr_pttrue%i",ptbin);
    TH2F * h2 = new TH2F(name,";#eta;< p^{jet}_{T} / p^{true}_{T} >",20,-5,5,51,0,2);
    h2EtaUncorr.push_back(h2);
    
    sprintf(name,"h2EtaCorr_pttrue%i",ptbin);
    h2EtaCorr.push_back(static_cast<TH2F*>(h2->Clone(name)));
    
    sprintf(name,"h2EtaCorrL2L3_pttrue%i",ptbin);
    h2EtaCorrL2L3.push_back(static_cast<TH2F*>(h2->Clone(name)));
  } // End of loop over pttrue bins

  // Create 2D histograms of response vs pttrue in
  // eta bins
  std::vector<TH2F*> h2PttrueUncorr;    // Uncorrected response vs pttrue
  std::vector<TH2F*> h2PttrueCorr;      // Response corrected by kalibri fit vs pttrue
  std::vector<TH2F*> h2PttrueCorrL2L3;  // Response corrected by JetMET L2L3 correction vs pttru
  // Logarithmic binning
  const int nLogBins = 15;
  double logBins[nLogBins+1];
  equidistLogBins(logBins,nLogBins,bins.xLow(0),bins.xUp(bins.nBinsX()-1));
  for(int etabin = 0; etabin < bins.nBinsY(); etabin++) { // Loop over eta bins
    char name[50];
    sprintf(name,"h2PttrueUncorr_pttrue%i",etabin);
    TH2F * h2 = new TH2F(name,";p^{true}_{T} (GeV);< p^{jet}_{T} / p^{true}_{T} >",nLogBins,logBins,51,0,2);
    h2PttrueUncorr.push_back(h2);
    
    sprintf(name,"h2PttrueCorr_eta%i",etabin);
    h2PttrueCorr.push_back(static_cast<TH2F*>(h2->Clone(name)));
    
    sprintf(name,"h2PttrueCorrL2L3_eta%i",etabin);
    h2PttrueCorrL2L3.push_back(static_cast<TH2F*>(h2->Clone(name)));
  } // End of loop over eta bins

  // Other control quantities
  TH1F * hPttrueSpectrum = new TH1F("hPttrueSpectrum",";p^{true}_{T} (GeV);dN / dp^{true}_{T}  1 / (GeV)",50,bins.xLow(0),bins.xUp(bins.nBinsX()-1));
  hPttrueSpectrum->Sumw2();
  hPttrueSpectrum->SetLineWidth(2);
  objToBeWritten.push_back(hPttrueSpectrum);


  // First loop over data and fill response of uncorrected
  // and Kalibri corrected data into the corresponding
  // 2D histogram
  for( data_it = data_->begin(); data_it != data_->end(); data_it++ ) {
    JetTruthEvent *jte = dynamic_cast<JetTruthEvent*>(*data_it);
    if( jte ) {
      if( jte->FlaggedBad() ) continue;   // Discard events flagged bad
      
      Jet * jet         = static_cast<Jet*>(jte->GetMess());
      double weight     = jte->GetWeight();
      double ptmeas     = jet->pt;
      double ptcorr     = jte->GetParametrizedMess();
      double eta        = jte->GetMess()->eta;
      double pttrue     = jte->GetTruth();

      // Find pttrue bin for response vs eta plot
      int ptbin = bins.iX(pttrue);
      if( 0 <= ptbin && ptbin < bins.nBinsX() ) {
	h2EtaUncorr.at(ptbin)->Fill(eta,ptmeas/pttrue,weight);
	h2EtaCorr.at(ptbin)->Fill(eta,ptcorr/pttrue,weight);
      }

      // Find pttrue bin for response vs eta plot
      int etabin = bins.iY(eta);
      if( 0 <= etabin && etabin < bins.nBinsY() ) {
	h2PttrueUncorr.at(etabin)->Fill(pttrue,ptmeas/pttrue,weight);
	h2PttrueCorr.at(etabin)->Fill(pttrue,ptcorr/pttrue,weight);
      }

      // Pttrue spectrum
      hPttrueSpectrum->Fill(pttrue,weight);
    }
  } // End of first loop over data

  // Either use L2L3 correction stored in TData or read from file
  bool useReadConstants = readJetMETParameters();

  // Second loop over data and fill response of uncorrected
  // and Kalibri corrected data into the corresponding
  // 2D histogram
  for( data_it = data_->begin(); data_it != data_->end(); data_it++ ) {
    JetTruthEvent *jte = dynamic_cast<JetTruthEvent*>(*data_it);
    if( jte ) {
      if( jte->FlaggedBad() ) continue;   // Discard events flagged bad

      Jet * jet         = static_cast<Jet*>(jte->GetMess());
      double weight     = jte->GetWeight();
      double ptcorr     = 0.;
      if( useReadConstants ) ptcorr = jte->GetParametrizedMess();
      else                   ptcorr = jet->corFactors.getL2L3() * jet->pt;
      double eta        = jte->GetMess()->eta;
      double pttrue     = jte->GetTruth();
      // Find pttrue bin for response vs eta plot
      int ptbin = bins.iX(pttrue);
      if( 0 <= ptbin && ptbin < bins.nBinsX() ) {
	h2EtaCorrL2L3.at(ptbin)->Fill(eta,ptcorr/pttrue,weight);
      }

      // Find pttrue bin for response vs eta plot
      int etabin = bins.iY(eta);
      if( 0 <= etabin && etabin < bins.nBinsY() ) {
	h2PttrueCorrL2L3.at(etabin)->Fill(pttrue,ptcorr/pttrue,weight);
      }
    }
  } // End of second loop over data

  // If parameter values were changed, reset them
  // to fitted values
  if( useReadConstants ) resetFittedParameters();


  // -- Mean response vs eta, pt ---------------------------------
  // Indices:
  //  0 - Mean
  //  1 - Sigma
  //  2 - Mean of Gauss
  //  3 - Width of Gauss
  // Suffix:
  //  Uncorr   - Uncorrected mean response
  //  Corr     - Mean response corrected by kalibri fit
  //  CorrL2L3 - Mean response corrected by JetMET L2L3 correction
  std::vector< std::vector<TH1F*> > hRespEtaUncorr(4,std::vector<TH1F*>(bins.nBinsX()));
  std::vector< std::vector<TH1F*> > hRespEtaCorr(4,std::vector<TH1F*>(bins.nBinsX()));
  std::vector< std::vector<TH1F*> > hRespEtaCorrL2L3(4,std::vector<TH1F*>(bins.nBinsX()));

  std::vector< std::vector<TH1F*> > hRespPttrueUncorr(4,std::vector<TH1F*>(bins.nBinsY()));
  std::vector< std::vector<TH1F*> > hRespPttrueCorr(4,std::vector<TH1F*>(bins.nBinsY()));
  std::vector< std::vector<TH1F*> > hRespPttrueCorrL2L3(4,std::vector<TH1F*>(bins.nBinsY()));

  for(int ptbin = 0; ptbin < bins.nBinsX(); ptbin++) { // Loop over pt bins
    char title[50];
    sprintf(title,"%.1f < p^{true}_{T} < %.1f GeV",bins.xLow(bins.bin(ptbin,0)),bins.xUp(bins.bin(ptbin,0)));
    TH1F * hProjection[8];
    TH1F * gaussplots[4];
    TF1 * gaussfits[4];
    // Uncorrected response
    fit2D(h2EtaUncorr.at(ptbin),hProjection,gaussplots,gaussfits);
    for(int i = 0; i < 8; i++) {
      if( i < 4 ) {
	hProjection[i]->SetTitle(title);
	hProjection[i]->SetMarkerStyle(20);
	hProjection[i]->SetMarkerColor(1);
	hProjection[i]->SetLineColor(1);
	hRespEtaUncorr.at(i).at(ptbin) = hProjection[i];
	objToBeWritten.push_back(hProjection[i]);
	delete gaussplots[i];
	delete gaussfits[i];
      } else {
	delete hProjection[i];
      }
    }
    
    // Kalibri corrected response
    fit2D(h2EtaCorr.at(ptbin),hProjection,gaussplots,gaussfits);
    for(int i = 0; i < 8; i++) {
      if( i < 4 ) {
	hProjection[i]->SetTitle(title);
	hProjection[i]->SetMarkerStyle(21);
	hProjection[i]->SetMarkerColor(2);
	hProjection[i]->SetLineColor(2);
	hRespEtaCorr.at(i).at(ptbin) = hProjection[i];
	objToBeWritten.push_back(hProjection[i]);
	delete gaussplots[i];
	delete gaussfits[i];
      } else {
	delete hProjection[i];
      }
    }

    // L2L3 corrected response
    fit2D(h2EtaCorrL2L3.at(ptbin),hProjection,gaussplots,gaussfits);
    for(int i = 0; i < 8; i++) {
      if( i < 4 ) {
	hProjection[i]->SetTitle(title);
	hProjection[i]->SetMarkerStyle(23);
	hProjection[i]->SetMarkerColor(4);
	hProjection[i]->SetLineColor(4);
	hRespEtaCorrL2L3.at(i).at(ptbin) = hProjection[i];
	objToBeWritten.push_back(hProjection[i]);
	delete gaussplots[i];
	delete gaussfits[i];
      } else {
	delete hProjection[i];
      }
    }
  } // End of loop over pt bins

  //Delete 2D histograms response vs eta
  for(int i = 0; i < bins.nBinsX(); i++) {
    delete h2EtaUncorr.at(i);
    delete h2EtaCorr.at(i);
    delete h2EtaCorrL2L3.at(i);
  }

  for(int etabin = 0; etabin < bins.nBinsY(); etabin++) { // Loop over eta bins
    char title[50];
    sprintf(title,"%.1f <  #eta < %.1f",bins.yLow(bins.bin(0,etabin)),bins.yUp(bins.bin(0,etabin)));
    TH1F * hProjection[8];
    TH1F * gaussplots[4];
    TF1 * gaussfits[4];
    // Uncorrected response
    fit2D(h2PttrueUncorr.at(etabin),hProjection,gaussplots,gaussfits);
    for(int i = 0; i < 8; i++) {
      if( i < 4 ) {
	hProjection[i]->SetTitle(title);
	hProjection[i]->SetMarkerStyle(20);
	hProjection[i]->SetMarkerColor(1);
	hProjection[i]->SetLineColor(1);
	hRespPttrueUncorr.at(i).at(etabin) = hProjection[i];
	objToBeWritten.push_back(hProjection[i]);
	delete gaussplots[i];
	delete gaussfits[i];
      } else {
	delete hProjection[i];
      }
    }
    
    // Kalibri corrected response
    fit2D(h2PttrueCorr.at(etabin),hProjection,gaussplots,gaussfits);
    for(int i = 0; i < 8; i++) {
      if( i < 4 ) {
	hProjection[i]->SetTitle(title);
	hProjection[i]->SetMarkerStyle(21);
	hProjection[i]->SetMarkerColor(2);
	hProjection[i]->SetLineColor(2);
	hRespPttrueCorr.at(i).at(etabin) = hProjection[i];
	objToBeWritten.push_back(hProjection[i]);
	delete gaussplots[i];
	delete gaussfits[i];
      } else {
	delete hProjection[i];
      }
    }

    // L2L3 corrected response
    fit2D(h2PttrueCorrL2L3.at(etabin),hProjection,gaussplots,gaussfits);
    for(int i = 0; i < 8; i++) {
      if( i < 4 ) {
	hProjection[i]->SetTitle(title);
	hProjection[i]->SetMarkerStyle(23);
	hProjection[i]->SetMarkerColor(4);
	hProjection[i]->SetLineColor(4);
	hRespPttrueCorrL2L3.at(i).at(etabin) = hProjection[i];
	objToBeWritten.push_back(hProjection[i]);
	delete gaussplots[i];
	delete gaussfits[i];
      } else {
	delete hProjection[i];
      }
    }
  } // End of loop over pt bins


  //Delete 2D histograms response vs pttrue
  for(int i = 0; i < bins.nBinsY(); i++) {
    delete h2PttrueUncorr.at(i);
    delete h2PttrueCorr.at(i);
    delete h2PttrueCorrL2L3.at(i);
  }


  // -- Draw histograms into ps file ----------------------------
  TLegend * leg = new TLegend(0.4,0.65,0.93,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(hRespEtaUncorr.at(0).at(0),"Uncorrected","P");
  leg->AddEntry(hRespEtaCorr.at(0).at(0),"Kalibri","P");
  leg->AddEntry(hRespEtaCorrL2L3.at(0).at(0),(config_->read<string>("Control plots JetMET L2L3 label","NO LABEL")).c_str(),"P");

  // Mean response vs eta per pttrue bin
  TPostScript       * ps = new TPostScript("controlplotsJetTruthEventResponse.ps",111);
  TCanvas     * const c1 = new TCanvas("c1","",500,500);
  c1->cd();

  TH1F * h = hRespEtaUncorr.at(0).at(0);
  TLine *leta = new TLine(h->GetXaxis()->GetBinLowEdge(1),1,
			  h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()),1);
  leta->SetLineColor(1);
  leta->SetLineWidth(2);
  leta->SetLineStyle(2);

  for(int ptbin = 0; ptbin < bins.nBinsX(); ptbin++) {
    // Mean response vs eta
    hRespEtaUncorr.at(0).at(ptbin)->GetYaxis()->SetRangeUser(0,2);
    hRespEtaUncorr.at(0).at(ptbin)->GetYaxis()->SetTitle("< p^{jet}_{T} / p^{true}_{T} >");
    hRespEtaUncorr.at(0).at(ptbin)->Draw("PE1");
    leta->Draw("same");
    hRespEtaCorr.at(0).at(ptbin)->Draw("PE1same");
    hRespEtaCorrL2L3.at(0).at(ptbin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();

    // Mean of Gauss fit on response vs eta
    hRespEtaUncorr.at(2).at(ptbin)->GetYaxis()->SetRangeUser(0,2);
    hRespEtaUncorr.at(2).at(ptbin)->GetYaxis()->SetTitle("Gauss fit < p^{jet}_{T} / p^{true}_{T} >");
    hRespEtaUncorr.at(2).at(ptbin)->Draw("PE1");
    leta->Draw("same");
    hRespEtaCorr.at(2).at(ptbin)->Draw("PE1same");
    hRespEtaCorrL2L3.at(2).at(ptbin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
  }
  for(int ptbin = 0; ptbin < bins.nBinsX(); ptbin++) {
    // Zoom: Mean response vs eta
    hRespEtaUncorr.at(0).at(ptbin)->GetYaxis()->SetRangeUser(0.8,1.2);
    hRespEtaUncorr.at(0).at(ptbin)->Draw("PE1");
    leta->Draw("same");
    hRespEtaCorr.at(0).at(ptbin)->Draw("PE1same");
    hRespEtaCorrL2L3.at(0).at(ptbin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
    // Zoom: Mean of Gauss fit on response vs eta
    hRespEtaUncorr.at(2).at(ptbin)->GetYaxis()->SetRangeUser(0.8,1.2);
    hRespEtaUncorr.at(2).at(ptbin)->GetYaxis()->SetTitle("Gauss fit < p^{jet}_{T} / p^{true}_{T} >");
    hRespEtaUncorr.at(2).at(ptbin)->Draw("PE1");
    leta->Draw("same");
    hRespEtaCorr.at(2).at(ptbin)->Draw("PE1same");
    hRespEtaCorrL2L3.at(2).at(ptbin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
  }

  // Mean response vs pttrue per eta bin
  h = hRespPttrueUncorr.at(0).at(0);
  TLine *lpttrue = new TLine(h->GetXaxis()->GetBinLowEdge(1),1,
			  h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()),1);
  lpttrue->SetLineColor(1);
  lpttrue->SetLineWidth(2);
  lpttrue->SetLineStyle(2);
  for(int etabin = 0; etabin < bins.nBinsY(); etabin++) {
    // Mean response vs pttrue
    hRespPttrueUncorr.at(0).at(etabin)->GetYaxis()->SetRangeUser(0,2);
    hRespPttrueUncorr.at(0).at(etabin)->GetYaxis()->SetTitle("< p^{jet}_{T} / p^{true}_{T} >");
    hRespPttrueUncorr.at(0).at(etabin)->Draw("PE1");
    lpttrue->Draw("same");
    hRespPttrueCorr.at(0).at(etabin)->Draw("PE1same");
    hRespPttrueCorrL2L3.at(0).at(etabin)->Draw("PE1same");
    leg->Draw("same");
    c1->SetLogx(1);
    c1->Draw();
    ps->NewPage();
    // Mean response of Gauss fit vs pttrue
    hRespPttrueUncorr.at(2).at(etabin)->GetYaxis()->SetRangeUser(0,2);
    hRespPttrueUncorr.at(2).at(etabin)->GetYaxis()->SetTitle("Gauss fit < p^{jet}_{T} / p^{true}_{T} >");
    hRespPttrueUncorr.at(2).at(etabin)->Draw("PE1");
    lpttrue->Draw("same");
    hRespPttrueCorr.at(2).at(etabin)->Draw("PE1same");
    hRespPttrueCorrL2L3.at(2).at(etabin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
  }
  for(int etabin = 0; etabin < bins.nBinsY(); etabin++) {
    // Zoom: Mean response vs pttrue
    hRespPttrueUncorr.at(0).at(etabin)->GetYaxis()->SetRangeUser(0.8,1.2);
    hRespPttrueUncorr.at(0).at(etabin)->Draw("PE1");
    lpttrue->Draw("same");
    hRespPttrueCorr.at(0).at(etabin)->Draw("PE1same");
    hRespPttrueCorrL2L3.at(0).at(etabin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
    // Zoom: Mean response of Gauss fit vs pttrue
    hRespPttrueUncorr.at(2).at(etabin)->GetYaxis()->SetRangeUser(0.8,1.2);
    hRespPttrueUncorr.at(2).at(etabin)->Draw("PE1");
    lpttrue->Draw("same");
    hRespPttrueCorr.at(2).at(etabin)->Draw("PE1same");
    hRespPttrueCorrL2L3.at(2).at(etabin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
  }
  c1->cd();
  c1->SetLogx(0);
  c1->SetLogy(1);
  hPttrueSpectrum->Draw("E");
  c1->Draw();
  ps->NewPage();

  ps->Close();
  delete ps;
  delete leta;
  delete lpttrue;

  // Resolution vs eta, pttrue
  ps = new TPostScript("controlplotsJetTruthEventResolution.ps",111);
  c1->cd();
  c1->SetLogx(0);
  c1->SetLogy(0);

  // Resolution vs eta per pttrue bin
  for(int ptbin = 0; ptbin < bins.nBinsX(); ptbin++) {
    // Sigma
    hRespEtaUncorr.at(1).at(ptbin)->GetYaxis()->SetRangeUser(0,0.4);
    hRespEtaUncorr.at(1).at(ptbin)->GetYaxis()->SetTitle("#sigma( p^{jet}_{T} / p^{true}_{T} )  /  < p^{jet}_{T} / p^{true}_{T} >");
    hRespEtaUncorr.at(1).at(ptbin)->Draw("PE1");
    hRespEtaCorr.at(1).at(ptbin)->Draw("PE1same");
    hRespEtaCorrL2L3.at(1).at(ptbin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
    // Width of Gauss fit
    hRespEtaUncorr.at(3).at(ptbin)->GetYaxis()->SetRangeUser(0,0.4);
    hRespEtaUncorr.at(3).at(ptbin)->GetYaxis()->SetTitle("Gauss fit #sigma( p^{jet}_{T} / p^{true}_{T} )  /  < p^{jet}_{T} / p^{true}_{T} >");
    hRespEtaUncorr.at(3).at(ptbin)->Draw("PE1");
    hRespEtaCorr.at(3).at(ptbin)->Draw("PE1same");
    hRespEtaCorrL2L3.at(3).at(ptbin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
  }
  // Resolution vs pttrue per eta bin
  for(int etabin = 0; etabin < bins.nBinsY(); etabin++) {
    // Sigma
    hRespPttrueUncorr.at(1).at(etabin)->GetYaxis()->SetRangeUser(0,0.4);
    hRespPttrueUncorr.at(1).at(etabin)->GetYaxis()->SetTitle("#sigma( p^{jet}_{T} / p^{true}_{T} )  /  < p^{jet}_{T} / p^{true}_{T} >");
    hRespPttrueUncorr.at(1).at(etabin)->Draw("PE1");
    hRespPttrueCorr.at(1).at(etabin)->Draw("PE1same");
    hRespPttrueCorrL2L3.at(1).at(etabin)->Draw("PE1same");
    leg->Draw("same");
    c1->SetLogx(1);
    c1->Draw();
    ps->NewPage();
    // Width of Gauss fit
    hRespPttrueUncorr.at(3).at(etabin)->GetYaxis()->SetRangeUser(0,0.4);
    hRespPttrueUncorr.at(3).at(etabin)->GetYaxis()->SetTitle("Gauss fit #sigma( p^{jet}_{T} / p^{true}_{T} )  /  < p^{jet}_{T} / p^{true}_{T} >");
    hRespPttrueUncorr.at(3).at(etabin)->Draw("PE1");
    hRespPttrueCorr.at(3).at(etabin)->Draw("PE1same");
    hRespPttrueCorrL2L3.at(3).at(etabin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
  }

  if( outputROOT_ ) writeToRootFile(objToBeWritten,"JetTruthEventResponse");

  // Clean up
  ps->Close();
  std::vector<TObject*>::iterator it = objToBeWritten.begin();
  for(; it != objToBeWritten.end(); it++) {
    delete *it;
  }
  delete c1;
  delete ps;
  delete leg;
}



//!  \brief Control plots for comparison with L2L3 MC truth correction
//!
//!  The plots are written to the file
//!  "controlplotsL2L3MCTruth.ps" and (if enabled)
//!  to the directory "L2L3MCTruth" in the file
//!  "controlplots.root".
// -------------------------------------------------------------
void TControlPlots::makeControlPlotsL2L3MCTruth() {
  std::vector<TObject*>                objToBeWritten;
  std::vector<TData*>::const_iterator  data_it;


  // Absolute eta bins
  std::vector<double> etaBinEdge = bag_of<double>(config_->read<std::string>("Control plots eta bin edges",""));

  std::vector<double>::iterator itNeg = etaBinEdge.end();
  std::vector<double>::iterator it = etaBinEdge.begin();
  for(; it != etaBinEdge.end(); it++) {
    if( *it < 0 ) itNeg = it;
  }
  if( itNeg != etaBinEdge.end() ) {
    etaBinEdge.erase(etaBinEdge.begin(),itNeg+1); // Want symmetry in eta
    etaBinEdge.insert(etaBinEdge.begin(),0.);
  }

  // Find histogram ranges
  double ptGenMin = 10000.;
  double ptGenMax = 0.;
  for( data_it = data_->begin(); data_it != data_->end(); data_it++ ) {
    if( (*data_it)->GetTruth() < ptGenMin ) ptGenMin = (*data_it)->GetTruth();
    if( (*data_it)->GetTruth() > ptGenMax ) ptGenMax = (*data_it)->GetTruth();
  }
  ptGenMax *= 1.1;
  ptGenMin *= 0.9;


  // Create histograms:
  // Gen Pt in different eta bins
  std::vector<TH1F*> hGenPtSpec;
  for(unsigned int i = 1; i < etaBinEdge.size(); i++) {
    char name[50];
    sprintf(name,"hGenPtSpec_eta%i",i-1);
    TH1F * h = new TH1F(name,";p^{gen}_{T} (GeV);dN_{jet} / dp^{gen}_{T}  1 / (GeV)",
			50,ptGenMin,ptGenMax);
    h->GetXaxis()->SetNdivisions(505);
    h->SetMarkerStyle(19+i);
    int color = i;
    if( i >= 3 ) color++; // Avoid yellow
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    hGenPtSpec.push_back(h);
    objToBeWritten.push_back(h);
  }


  // Loop over data and fill histograms
  for( data_it = data_->begin(); data_it != data_->end(); data_it++ ) {
    JetTruthEvent *jte = dynamic_cast<JetTruthEvent*>(*data_it);
    if( jte ) {
      if( jte->FlaggedBad() ) continue;    // Discard events flagged bad

      double weight = jte->GetWeight();
      double ptgen  = jte->GetTruth();
      double eta    = jte->GetMess()->eta;

      unsigned int etaBin = 0;
      if( std::abs(eta) > etaBinEdge.back() ) continue;
      while( std::abs(eta) > etaBinEdge.at(etaBin+1) ) etaBin++;
      hGenPtSpec.at(etaBin)->Fill( weight * ptgen );
    }
  }


  // Draw histograms
  TPostScript * const ps = new TPostScript("controlplotsL2L3MCTruth.ps",111);
  TCanvas     * const c1 = new TCanvas("c1","",500,500);
  c1->cd();

  double legminx = 0.6;
  double legminy = 0.65;
  TLegend * leg = new TLegend(legminx,legminy,0.9,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  for(unsigned int i = 1; i < etaBinEdge.size(); i++) {
    char entry[50];
    sprintf(entry,"%.1f < |#eta| < %.1f",etaBinEdge.at(i-1),etaBinEdge.at(i));
    leg->AddEntry(hGenPtSpec.at(i-1),entry,"P");
  }

  hGenPtSpec.at(0)->GetYaxis()->SetRangeUser(0.1,15*(hGenPtSpec.at(0)->GetMaximum()));

  for(unsigned int i = 0; i < hGenPtSpec.size(); i++) {
    if( i == 0 ) hGenPtSpec.at(i)->Draw("PE1");
    else         hGenPtSpec.at(i)->Draw("PE1same");
  }
  c1->SetLogy(1);
  leg->Draw("same");
  c1->Draw();
  ps->NewPage();


  if( outputROOT_ ) writeToRootFile(objToBeWritten, "L2L3MCTruth");
  
  ps->Close();


  // Clean up
  delete c1;
  delete ps;
  delete leg;
  std::vector<TObject*>::iterator objit = objToBeWritten.begin();
  for(; objit != objToBeWritten.end(); objit++) {
    delete *objit;
  }
  objToBeWritten.clear();
  hGenPtSpec.clear();
}



//!  \brief Top Control Histograms
// -------------------------------------------------------------
void TControlPlots::makeControlPlotsTop()
{
  std::vector<TObject*> objToBeWritten;

  TCanvas * const c = new TCanvas("c","",600,600);

  TPostScript * const ps = new TPostScript("controlplotsTop.ps",111);

  bool individualPdf = config_->read<bool>("create individual pdf files", false);

  // book hists

  //  double binningLogPt[20] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
  //			     120, 150, 200, 250, 300, 400, 500, 700, 1000, 2000};
  double binningLogPt[15] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
			     120, 150, 200, 250, 300};

  unsigned int nTowersPhi = 72;
  double edgePhiHist = TMath::Pi()+2*TMath::Pi()/nTowersPhi;

  TH1F* scale  = new TH1F("scale" , "", 200,    0, 200);
  TH1F* weight = new TH1F("weight", "", 200,    0, 200);
  TH1F* truth  = new TH1F("truth" , "", 200,    0, 200);
  TH1F* pt     = new TH1F("pt"    , "", 200,    0, 200);
  TH1F* eta    = new TH1F("eta"   , "",  96,  -4.,  4.);
  TH1F* phi    = new TH1F("phi"   , "",  nTowersPhi+2, -edgePhiHist, edgePhiHist);
  TH1F* genPt  = new TH1F("genPt" , "", 200,    0, 200);

  TH2F* recPt12 = new TH2F("recPt12", "", 100,  0., 200., 100,  0., 200.);
  TH2F* genPt12 = new TH2F("genPt12", "", 100,  0., 200., 100,  0., 200.);
  TH2F* eta12   = new TH2F("eta12"  , "",  80, -4.,   4.,  80, -4.,   4.);

  TH1F* deltaRecPt = new TH1F("deltaRecPt", "",  75,  0., 150.);
  TH1F* deltaGenPt = new TH1F("deltaGenPt", "",  75,  0., 150.);
  TH1F* deltaEta   = new TH1F("deltaEta"  , "",  24,  0.,   4.);
  TH2F* deltaRecPtRecPt1 = new TH2F("deltaRecPtRecPt1", "", 100,  0., 200.,  75,  0., 150.);
  TH2F* deltaRecPtEta1   = new TH2F("deltaRecPtEta1"  , "",  80, -4.,   4.,  75,  0., 150.);

  TH1F* meanPt  = new TH1F("meanPt" , "", 75,  0., 150.);
  TH1F* meanEta = new TH1F("meanEta", "", 80, -4.,   4.);

  TH1F* wPt  = new TH1F("wPt" , "", 150,  0., 300.);
  TH1F* wEta = new TH1F("wEta", "",  80, -4.,   4.);

  TH1F* dR      = new TH1F("deltaR"     , "",  40,  0.,   4.);
  TH2F* dRrecPt = new TH2F("deltaRrecPt", "", 100,  0., 200., 40,   0.,  4.);
  TH2F* dRgenPt = new TH2F("deltaRgenPt", "", 100,  0., 200., 40,   0.,  4.);
  TH2F* dReta   = new TH2F("deltaReta"  , "",  80, -4.,   4., 40,   0.,  4.);
  TH2F* dRwPt   = new TH2F("deltaRwPt"  , "", 100,  0., 200., 40,   0.,  4.);

  TH1F* invMass       [2];
  TH1F* messTruth     [2];
  TH2F* messTruthPtGen[2];
  TH2F* messTruthPtRec[2];
  TH2F* messTruthEta  [2];
  TH2F* responsePtGen [2];
  TH2F* responsePtRec [2];
  TH2F* responseEta   [2];
  TString suffix[2] = { "Before", "After" };
  for(unsigned a=0; a<2; a++){
    invMass       [a] = new TH1F("invMass"        +suffix[a], "", 100, 0., 200.);
    messTruth     [a] = new TH1F("messTruth"      +suffix[a], "",  40, 0.,   2.);
    messTruthPtGen[a] = new TH2F("messTruthPtGen" +suffix[a], "", 14, binningLogPt, 51, 0., 2.);
    messTruthPtRec[a] = new TH2F("messTruthPtRec" +suffix[a], "", 14, binningLogPt, 51, 0., 2.);
    messTruthEta  [a] = new TH2F("messTruthEta"   +suffix[a], "", 40,  -4.,  4.   , 51, 0., 2.);
    responsePtGen [a] = new TH2F("responsePtGen"  +suffix[a], "", 14, binningLogPt, 51, 0., 2.);
    responsePtRec [a] = new TH2F("responsePtRec"  +suffix[a], "", 14, binningLogPt, 51, 0., 2.);
    responseEta   [a] = new TH2F("responseEta"    +suffix[a], "", 40,  -4.,  4.   , 51, 0., 2.);
  }
  TH1F* messTruthPtGen_proj[2][4];
  TH1F* messTruthPtRec_proj[2][4];
  TH1F* messTruthEta_proj  [2][4];
  TH1F* responsePtGen_proj [2][4];
  TH1F* responsePtRec_proj [2][4];
  TH1F* responseEta_proj   [2][4];

  TH2F* corrFacsPt [5];
  TH2F* corrFacsEta[5];
  TH2F* correctedResponsePt [6];
  TH2F* correctedResponseEta[6];
  TString corrLevels[6] = { "Raw", "L1", "L2", "L3", "L4", "L5" };
  for(unsigned i=0; i<5; i++){
    corrFacsPt [i] = new TH2F("corrFacsPt" +corrLevels[i+1], "", 14, binningLogPt, 51, 0., 2.);
    corrFacsEta[i] = new TH2F("corrFacsEta"+corrLevels[i+1], "", 40,  -4.,  4.   , 51, 0., 2.);
  }
  for(unsigned i=0; i<6; i++){
    correctedResponsePt [i] = new TH2F("correctedResponsePt" +corrLevels[i], "", 14, binningLogPt, 51, 0., 2.);
    correctedResponseEta[i] = new TH2F("correctedResponseEta"+corrLevels[i], "", 40,  -4.,  4.   , 51, 0., 2.);
  }
  TH1F* corrFacsPt_mgf [5];
  TH1F* corrFacsEta_mgf[5];
  TH1F* correctedResponsePt_mgf [6];
  TH1F* correctedResponseEta_mgf[6];

  //loop over all fit-events and fill hists

  for( std::vector<TData*>::const_iterator i = data_->begin() ; i != data_->end() ; ++i )  
    {

      TData* data = *i;
      if(data->GetType() != InvMass) continue;

      std::vector<TJet*> jets;
      
      TData_InvMass2 *invM2 = dynamic_cast<TData_InvMass2*>(data);	
      TwoJetsInvMassEvent* ev = (TwoJetsInvMassEvent*)data;
      if(invM2) {	
	jets.push_back( (TJet*) invM2->GetMess() );
	for(unsigned j=0; j<invM2->MultMessSize(); j++)
	  jets.push_back( (TJet*) (*invM2->GetSecondaryJets())[j]->GetMess() );
      } 
      else if(ev) {
	jets.push_back( ev->GetJet1() );
	jets.push_back( ev->GetJet2() );
      } else { 
      continue;
      }

      TLorentzVector jet4Vec;
      TLorentzVector recW;
      for(unsigned j=0; j<jets.size(); j++) {
	jet4Vec.SetPtEtaPhiE( jets[j]->pt, jets[j]->eta, jets[j]->phi, jets[j]->E );
	if(j==0) recW = jet4Vec;
	else recW += jet4Vec;
      }

      double s = 0.;
      double t = 0.;
      double w = 0.;      
      if(invM2) {
	s = invM2->GetScale();
	w = invM2->GetWeight();
	t = invM2->GetTruth();
      }
      if(ev) {
	s = ev->GetTruth();
	w = ev->GetWeight();
	t = ev->GetTruth();
      }
      scale ->Fill(s);
      weight->Fill(w);
      truth ->Fill(t);

      if(jets[1]) {

	recPt12->Fill(jets[0]->pt   , jets[1]->pt   );
	genPt12->Fill(jets[0]->genPt, jets[1]->genPt);
	eta12  ->Fill(jets[0]->eta  , jets[1]->eta  );

	deltaRecPt->Fill( std::abs(jets[0]->pt    - jets[1]->pt    ) );
	deltaGenPt->Fill( std::abs(jets[0]->genPt - jets[1]->genPt ) );
	deltaEta  ->Fill( std::abs(jets[0]->eta   - jets[1]->eta   ) );

	deltaRecPtRecPt1->Fill( jets[0]->pt , std::abs(jets[0]->pt - jets[1]->pt) );
	deltaRecPtEta1  ->Fill( jets[0]->eta, std::abs(jets[0]->pt - jets[1]->pt) );

	double dr = deltaR(jets[0]->eta, jets[0]->phi,
			   jets[1]->eta, jets[1]->phi);
	dR->Fill(dr);
	dRrecPt->Fill(jets[0]->pt      , dr);
	dRgenPt->Fill(jets[0]->genPt   , dr);
	dReta  ->Fill(jets[0]->eta     , dr);
	dRwPt  ->Fill(recW.Pt(), dr);
      }

      wPt ->Fill( recW.Pt()  );
      wEta->Fill( recW.Eta() );
      
      invMass  [0]->Fill( recW.M()   );
      messTruth[0]->Fill( recW.M()/t );

      invMass  [1]->Fill( invM2 ? invM2->GetMessCombination()   : ev->correctedMass()   );
      messTruth[1]->Fill( invM2 ? invM2->GetMessCombination()/t : ev->correctedMass()/t );

      double mPt  = 0.;
      double mEta = 0.;

      //      messTruthPtGen [0]->Fill( recW.Pt() , recW.M()/t );
      //      messTruthPtGen [1]->Fill( recW.Pt() , invM2 ? invM2->GetMessCombination()/t : ev->correctedMass()/t );

      for(unsigned j=0; j<jets.size(); j++) {
	if(jets[j]->flavor != TJet::uds) continue;

	pt ->Fill( jets[j]->pt  );
	eta->Fill( jets[j]->eta );
	phi->Fill( jets[j]->phi );

	genPt->Fill( jets[j]->genPt );

	mPt  += jets[j]->pt;
	mEta += jets[j]->eta;

	messTruthPtGen[0]->Fill( jets[j]->genPt , recW.M()/t );
	messTruthPtRec[0]->Fill( jets[j]->pt    , recW.M()/t );
	messTruthEta  [0]->Fill( jets[j]->eta   , recW.M()/t );
	
	messTruthPtGen[1]->Fill( jets[j]->genPt , invM2 ? invM2->GetMessCombination()/t : ev->correctedMass()/t );
	messTruthPtRec[1]->Fill( jets[j]->pt    , invM2 ? invM2->GetMessCombination()/t : ev->correctedMass()/t );
	messTruthEta  [1]->Fill( jets[j]->eta   , invM2 ? invM2->GetMessCombination()/t : ev->correctedMass()/t );

	double response[2]; // before and after Kalibri fit
	if(j==0) {
	  response[0] = invM2 ? invM2->GetMess()->pt / jets[j]->genPt : 
	                           ev->GetMess()->pt / jets[j]->genPt;
	  response[1] = invM2 ? invM2->GetParametrizedMess() / jets[j]->genPt : 
	                           ev->GetParametrizedMess() / jets[j]->genPt;
	}
	else {
	  if(invM2) {
	    const TData_MessMess* mm = (*invM2->GetSecondaryJets())[j-1];
	    response[0] = mm->GetMess()->pt         / jets[j]->genPt;
	    response[1] = mm->GetParametrizedMess() / jets[j]->genPt;
	  } else {
	    response[0] = ev->GetMess2()->pt         / jets[j]->genPt;
	    response[1] = ev->GetParametrizedMess2() / jets[j]->genPt;
	  }
	}
	for(unsigned a=0; a<2; a++){
	  responsePtGen[a]->Fill( jets[j]->genPt , response[a] );
	  responsePtRec[a]->Fill( jets[j]->pt    , response[a] );
	  responseEta  [a]->Fill( jets[j]->eta   , response[a] );
	}

	corrFacsPt [0]->Fill( jets[j]->genPt  , jets[j]->corFactors.getL1() );
	corrFacsPt [1]->Fill( jets[j]->genPt  , jets[j]->corFactors.getL2() );
	corrFacsPt [2]->Fill( jets[j]->genPt  , jets[j]->corFactors.getL3() );
	corrFacsPt [3]->Fill( jets[j]->genPt  , jets[j]->corFactors.getL4() );
	corrFacsPt [4]->Fill( jets[j]->genPt  , jets[j]->corFactors.getL5() );

	corrFacsEta[0]->Fill( jets[j]->eta , jets[j]->corFactors.getL1() );
	corrFacsEta[1]->Fill( jets[j]->eta , jets[j]->corFactors.getL2() );
	corrFacsEta[2]->Fill( jets[j]->eta , jets[j]->corFactors.getL3() );
	corrFacsEta[3]->Fill( jets[j]->eta , jets[j]->corFactors.getL4() );
	corrFacsEta[4]->Fill( jets[j]->eta , jets[j]->corFactors.getL5() );

	correctedResponsePt [0]->Fill( jets[j]->genPt  ,                               jets[j]->pt/jets[j]->genPt );
	correctedResponsePt [1]->Fill( jets[j]->genPt  , jets[j]->corFactors.getL1()  *jets[j]->pt/jets[j]->genPt );
	correctedResponsePt [2]->Fill( jets[j]->genPt  , jets[j]->corFactors.getToL2()*jets[j]->pt/jets[j]->genPt );
	correctedResponsePt [3]->Fill( jets[j]->genPt  , jets[j]->corFactors.getToL3()*jets[j]->pt/jets[j]->genPt );
	correctedResponsePt [4]->Fill( jets[j]->genPt  , jets[j]->corFactors.getToL4()*jets[j]->pt/jets[j]->genPt );
	correctedResponsePt [5]->Fill( jets[j]->genPt  , jets[j]->corFactors.getToL5()*jets[j]->pt/jets[j]->genPt );

	correctedResponseEta[0]->Fill( jets[j]->eta ,                               jets[j]->pt/jets[j]->genPt );
	correctedResponseEta[1]->Fill( jets[j]->eta , jets[j]->corFactors.getL1()  *jets[j]->pt/jets[j]->genPt );
	correctedResponseEta[2]->Fill( jets[j]->eta , jets[j]->corFactors.getToL2()*jets[j]->pt/jets[j]->genPt );
	correctedResponseEta[3]->Fill( jets[j]->eta , jets[j]->corFactors.getToL3()*jets[j]->pt/jets[j]->genPt );
	correctedResponseEta[4]->Fill( jets[j]->eta , jets[j]->corFactors.getToL4()*jets[j]->pt/jets[j]->genPt );
	correctedResponseEta[5]->Fill( jets[j]->eta , jets[j]->corFactors.getToL5()*jets[j]->pt/jets[j]->genPt );
      }

      mPt  /= jets.size();
      mEta /= jets.size();

      meanPt ->Fill( mPt  );
      meanEta->Fill( mEta );

    }  //End of loop over all fit-events

  // perform Gauss fits for bins in 2-dim. hists

  for(unsigned i=0; i<2; i++){

    std::vector<TH1F*> results_ptGen;
    std::vector<TH1F*> results_ptRec;
    std::vector<TH1F*> results_eta;

    fit2D(messTruthPtGen[i], results_ptGen);
    fit2D(messTruthPtRec[i], results_ptRec);
    fit2D(messTruthEta  [i], results_eta);
    for(unsigned j=0; j<8; j++) {
      if(j<4) {
	messTruthPtGen_proj[i][j] = results_ptGen[j];
	messTruthPtRec_proj[i][j] = results_ptRec[j];
	messTruthEta_proj  [i][j] = results_eta  [j];
      }
      else {
	delete results_ptGen[j];
	delete results_ptRec[j];
	delete results_eta  [j];
      }
    }

  }

  for(unsigned i=0; i<2; i++){

    std::vector<TH1F*> results_ptGen;
    std::vector<TH1F*> results_ptRec;
    std::vector<TH1F*> results_eta;

    fit2D(responsePtGen[i], results_ptGen);
    fit2D(responsePtRec[i], results_ptRec);
    fit2D(responseEta  [i], results_eta);
    for(unsigned j=0; j<8; j++) {
      if(j<4) {
	responsePtGen_proj[i][j] = results_ptGen[j];
	responsePtRec_proj[i][j] = results_ptRec[j];
	responseEta_proj  [i][j] = results_eta  [j];
      }
      else {
	delete results_ptGen[j];
	delete results_ptRec[j];
	delete results_eta  [j];
      }
    }

  }

  for(unsigned i=0; i<5; i++){

    std::vector<TH1F*> results_pt;
    std::vector<TH1F*> results_eta;

    fit2D(corrFacsPt [i], results_pt);
    fit2D(corrFacsEta[i], results_eta);
    for(unsigned j=0; j<8; j++) {
      if(j==2) { //mean of Gauss fits
	corrFacsPt_mgf [i] = results_pt [j];
	corrFacsEta_mgf[i] = results_eta[j];
      }
      else {
	delete results_pt [j];
	delete results_eta[j];
      }
    }

  }

  for(unsigned i=0; i<6; i++){

    std::vector<TH1F*> results_pt;
    std::vector<TH1F*> results_eta;

    fit2D(correctedResponsePt [i], results_pt);
    fit2D(correctedResponseEta[i], results_eta);
    for(unsigned j=0; j<8; j++) {
      if(j==2) { //mean of Gauss fits
	correctedResponsePt_mgf [i] = results_pt [j];
	correctedResponseEta_mgf[i] = results_eta[j];
      }
      else {
	delete results_pt [j];
	delete results_eta[j];
      }
    }

  }

  // configure hists

  scale ->SetXTitle( "scale"  );
  weight->SetXTitle( "weight" );
  truth ->SetXTitle( "truth"  );
  pt   ->SetXTitle( "p_{T} (rec) [GeV]" );
  eta  ->SetXTitle( "#eta"              );
  phi  ->SetXTitle( "#phi"              );
  genPt->SetXTitle( "p_{T} (gen) [GeV]" );
  dR   ->SetXTitle( "#DeltaR_{jj}" );

  scale ->SetYTitle( "events" );
  weight->SetYTitle( "events" );
  truth ->SetYTitle( "events" );
  pt    ->SetYTitle( "events" );
  eta   ->SetYTitle( "events" );
  phi   ->SetYTitle( "events" );
  genPt ->SetYTitle( "events" );
  dR    ->SetYTitle( "events" );

  recPt12->SetXTitle( "p_{T,1} (rec) [GeV]" );
  genPt12->SetXTitle( "p_{T,1} (gen) [GeV]" );
  eta12  ->SetXTitle( "#eta_{1}" );

  recPt12->SetYTitle( "p_{T,2} (rec) [GeV]" );
  genPt12->SetYTitle( "p_{T,2} (gen) [GeV]" );
  eta12  ->SetYTitle( "#eta_{2}" );

  deltaRecPt->SetXTitle( "#Deltap_{T,jj} (rec) [GeV]" );
  deltaGenPt->SetXTitle( "#Deltap_{T,jj} (gen) [GeV]" );
  deltaEta  ->SetXTitle( "#Delta#eta_{jj}"            );

  deltaRecPt->SetYTitle( "events" );
  deltaGenPt->SetYTitle( "events" );
  deltaEta  ->SetYTitle( "events" );

  deltaRecPtRecPt1->SetXTitle( "p_{T,1} (rec) [GeV]" );
  deltaRecPtEta1  ->SetXTitle( "#eta_{1}" );

  deltaRecPtRecPt1->SetYTitle( "#Deltap_{T,jj} (rec) [GeV]" );
  deltaRecPtEta1  ->SetYTitle( "#Deltap_{T,jj} (rec) [GeV]" );

  meanPt ->SetXTitle( "(p_{T,1}+p_{T,2})/2 [GeV]"   );
  meanEta->SetXTitle( "(#eta_{1}+#eta_{2})/2" );
  wPt    ->SetXTitle( "p_{T,W} [GeV]"  );
  wEta   ->SetXTitle( "#eta_{W}" );

  meanPt ->SetYTitle( "events" );
  meanEta->SetYTitle( "events" );
  wPt    ->SetYTitle( "events" );
  wEta   ->SetYTitle( "events" );

  dRrecPt->SetXTitle( "p_{T,1} (rec) [GeV]" );
  dRgenPt->SetXTitle( "p_{T,1} (gen) [GeV]" );
  dReta  ->SetXTitle( "#eta_{1}" );
  dRwPt  ->SetXTitle( "p_{T,W} [GeV]"  );

  dRrecPt->SetYTitle( "#DeltaR_{jj}" );
  dRgenPt->SetYTitle( "#DeltaR_{jj}" );
  dReta  ->SetYTitle( "#DeltaR_{jj}" );
  dRwPt  ->SetYTitle( "#DeltaR_{jj}" );

  int markerColor[2] = { 2, 1 };
  int markerStyle[2] = { 22, 20 };

  TPaveText *paveText[2];
  paveText[0] = new TPaveText(0.7, 0.6 , 0.91, 0.75, "NDC");
  paveText[1] = new TPaveText(0.7, 0.45, 0.91, 0.6 , "NDC");

  TF1 *fInvMass  [2];
  TF1 *fMessTruth[2];

  for(unsigned a=0; a<2; a++){

    invMass     [a]->Fit( "gaus", "Q0" );
    messTruth   [a]->Fit( "gaus", "Q0" );

    double invMassMu      = invMass[a]  ->GetFunction("gaus")->GetParameter(1);
    double invMassSigma   = invMass[a]  ->GetFunction("gaus")->GetParameter(2);
    double messTruthMu    = messTruth[a]->GetFunction("gaus")->GetParameter(1);
    double messTruthSigma = messTruth[a]->GetFunction("gaus")->GetParameter(2);

    invMass  [a]->Fit( "gaus", "Q0", "", invMassMu  -1.5*invMassSigma  , invMassMu  +1.5*invMassSigma   );
    messTruth[a]->Fit( "gaus", "Q0", "", messTruthMu-1.5*messTruthSigma, messTruthMu+1.5*messTruthSigma );

    invMass  [a]->SetLineColor( markerColor[a] );
    messTruth[a]->SetLineColor( markerColor[a] );

    invMass  [a]->SetMarkerColor( markerColor[a] );
    messTruth[a]->SetMarkerColor( markerColor[a] );

    invMass  [a]->SetMarkerStyle( markerStyle[a] );
    messTruth[a]->SetMarkerStyle( markerStyle[a] );

    invMass  [a]->SetMarkerSize( 1.5 );
    messTruth[a]->SetMarkerSize( 1.5 );

    invMass  [a]->SetXTitle( "m_{jj} [GeV]" );
    messTruth[a]->SetXTitle( "m_{jj} / 80.4 GeV" );

    invMass  [a]->SetYTitle( "events" );
    messTruth[a]->SetYTitle( "events" );

    fInvMass  [a] = invMass  [a]->GetFunction("gaus");
    fMessTruth[a] = messTruth[a]->GetFunction("gaus");

    fInvMass  [a]->SetLineColor( markerColor[a] );
    fMessTruth[a]->SetLineColor( markerColor[a] );

    paveText[a]->SetFillColor( 0 );
    paveText[a]->SetBorderSize( 1 );
    paveText[a]->SetTextColor( markerColor[a] );
    paveText[a]->SetTextAlign( 12 );

    double mu    = fInvMass[a]->GetParameter(1);
    double sigma = fInvMass[a]->GetParameter(2);
    double relSigma = sigma / mu;

    char *tmpTxt = new char[100];

    sprintf(tmpTxt, "#mu = %4.1f GeV", mu);
    paveText[a]->AddText(tmpTxt);
    sprintf(tmpTxt, "#sigma = %4.1f GeV", sigma);
    paveText[a]->AddText(tmpTxt);
    sprintf(tmpTxt, "#sigma/#mu = %4.2f", relSigma);
    paveText[a]->AddText(tmpTxt);

    TString yTitle_messTruth_proj[4] = {"< m_{jj} / 80.4 GeV >",
					"#sigma( m_{jj} / 80.4 GeV )",
					"Gauss fit < m_{jj} / 80.4 GeV >",
					"Gauss fit #sigma( m_{jj} / 80.4 GeV )"};

    TString yTitle_reponse_proj[4] = {"< p_{T} (rec) / p_{T} (gen) >",
				      "#sigma( p_{T} (rec) / p_{T} (gen) )",
				      "Gauss fit < p_{T} (rec) / p_{T} (gen) >",
				      "Gauss fit #sigma( p_{T} (rec) / p_{T} (gen) )"};

    double yMin_proj[4] = {0.4, -0.01, 0.4, -0.01};
    double yMax_proj[4] = {1.6, 0.56, 1.6, 0.56};

    for(unsigned j=0; j<4; j++) {

      // messTruth projections

      messTruthPtGen_proj[a][j]->SetMarkerColor( markerColor[a] );
      messTruthPtRec_proj[a][j]->SetMarkerColor( markerColor[a] );
      messTruthEta_proj  [a][j]->SetMarkerColor( markerColor[a] );
      
      messTruthPtGen_proj[a][j]->SetMarkerStyle( markerStyle[a] );
      messTruthPtRec_proj[a][j]->SetMarkerStyle( markerStyle[a] );
      messTruthEta_proj  [a][j]->SetMarkerStyle( markerStyle[a] );
      
      messTruthPtGen_proj[a][j]->SetMarkerSize( 1.5 );
      messTruthPtRec_proj[a][j]->SetMarkerSize( 1.5 );
      messTruthEta_proj  [a][j]->SetMarkerSize( 1.5 );
      
      messTruthPtGen_proj[a][j]->SetTitle( "" );
      messTruthPtRec_proj[a][j]->SetTitle( "" );
      messTruthEta_proj  [a][j]->SetTitle( "" );
      
      messTruthPtGen_proj[a][j]->SetXTitle( "p_{T} (gen) [GeV]" );
      messTruthPtRec_proj[a][j]->SetXTitle( "p_{T} (rec) [GeV]" );
      messTruthEta_proj  [a][j]->SetXTitle( "#eta" );
      
      messTruthPtGen_proj[a][j]->SetYTitle( yTitle_messTruth_proj[j] );
      messTruthPtRec_proj[a][j]->SetYTitle( yTitle_messTruth_proj[j] );
      messTruthEta_proj  [a][j]->SetYTitle( yTitle_messTruth_proj[j] );
      
      messTruthPtGen_proj[a][j]->SetMinimum( yMin_proj[j] );
      messTruthPtRec_proj[a][j]->SetMinimum( yMin_proj[j] );
      messTruthEta_proj  [a][j]->SetMinimum( yMin_proj[j] );
      
      messTruthPtGen_proj[a][j]->SetMaximum( yMax_proj[j] );
      messTruthPtRec_proj[a][j]->SetMaximum( yMax_proj[j] );
      messTruthEta_proj  [a][j]->SetMaximum( yMax_proj[j] );

      // response projections
      
      responsePtGen_proj[a][j]->SetMarkerColor( markerColor[a] );
      responsePtRec_proj[a][j]->SetMarkerColor( markerColor[a] );
      responseEta_proj  [a][j]->SetMarkerColor( markerColor[a] );
      
      responsePtGen_proj[a][j]->SetMarkerStyle( markerStyle[a] );
      responsePtRec_proj[a][j]->SetMarkerStyle( markerStyle[a] );
      responseEta_proj  [a][j]->SetMarkerStyle( markerStyle[a] );
      
      responsePtGen_proj[a][j]->SetMarkerSize( 1.5 );
      responsePtRec_proj[a][j]->SetMarkerSize( 1.5 );
      responseEta_proj  [a][j]->SetMarkerSize( 1.5 );
      
      responsePtGen_proj[a][j]->SetTitle( "" );
      responsePtRec_proj[a][j]->SetTitle( "" );
      responseEta_proj  [a][j]->SetTitle( "" );
      
      responsePtGen_proj[a][j]->SetXTitle( "p_{T} (gen) [GeV]" );
      responsePtRec_proj[a][j]->SetXTitle( "p_{T} (rec) [GeV]" );
      responseEta_proj  [a][j]->SetXTitle( "#eta" );
      
      responsePtGen_proj[a][j]->SetYTitle( yTitle_reponse_proj[j] );
      responsePtRec_proj[a][j]->SetYTitle( yTitle_reponse_proj[j] );
      responseEta_proj  [a][j]->SetYTitle( yTitle_reponse_proj[j] );
      
      responsePtGen_proj[a][j]->SetMinimum( yMin_proj[j] );
      responsePtRec_proj[a][j]->SetMinimum( yMin_proj[j] );
      responseEta_proj  [a][j]->SetMinimum( yMin_proj[j] );
      
      responsePtGen_proj[a][j]->SetMaximum( yMax_proj[j] );
      responsePtRec_proj[a][j]->SetMaximum( yMax_proj[j] );
      responseEta_proj  [a][j]->SetMaximum( yMax_proj[j] );

    }

  }

  int markerColorCorr[6] = { 2, 3, 4, 6, 28, 1 };
  int markerStyleCorr[6] = { 22, 24, 21, 23, 25, 20 };

  for(unsigned i=0; i<5; i++){
    corrFacsPt_mgf [i]->SetMarkerColor( markerColorCorr[i+1] );
    corrFacsEta_mgf[i]->SetMarkerColor( markerColorCorr[i+1] );

    corrFacsPt_mgf [i]->SetMarkerStyle( markerStyleCorr[i+1] );
    corrFacsEta_mgf[i]->SetMarkerStyle( markerStyleCorr[i+1] );

    corrFacsPt_mgf [i]->SetMarkerSize( 1.5 );
    corrFacsEta_mgf[i]->SetMarkerSize( 1.5 );

    corrFacsPt_mgf [i]->SetTitle( "" );
    corrFacsEta_mgf[i]->SetTitle( "" );

    corrFacsPt_mgf [i]->SetXTitle( "p_{T} (gen) [GeV]" );
    corrFacsEta_mgf[i]->SetXTitle( "#eta" );

    corrFacsPt_mgf [i]->SetYTitle( "JEC factor" );
    corrFacsEta_mgf[i]->SetYTitle( "JEC factor" );

    corrFacsPt_mgf [i]->SetMinimum( 0.7 );
    corrFacsEta_mgf[i]->SetMinimum( 0.7 );

    corrFacsPt_mgf [i]->SetMaximum( 2.2 );
    corrFacsEta_mgf[i]->SetMaximum( 2.2 );
  }

  for(unsigned i=0; i<6; i++){
    correctedResponsePt_mgf [i]->SetMarkerColor( markerColorCorr[i] );
    correctedResponseEta_mgf[i]->SetMarkerColor( markerColorCorr[i] );

    correctedResponsePt_mgf [i]->SetMarkerStyle( markerStyleCorr[i] );
    correctedResponseEta_mgf[i]->SetMarkerStyle( markerStyleCorr[i] );

    correctedResponsePt_mgf [i]->SetMarkerSize( 1.5 );
    correctedResponseEta_mgf[i]->SetMarkerSize( 1.5 );

    correctedResponsePt_mgf [i]->SetTitle( "" );
    correctedResponseEta_mgf[i]->SetTitle( "" );

    correctedResponsePt_mgf [i]->SetXTitle( "p_{T} (gen) [GeV]" );
    correctedResponseEta_mgf[i]->SetXTitle( "#eta" );

    correctedResponsePt_mgf [i]->SetYTitle( "p_{T} (rec) / p_{T} (gen)" );
    correctedResponseEta_mgf[i]->SetYTitle( "p_{T} (rec) / p_{T} (gen)" );

    correctedResponsePt_mgf [i]->SetMinimum( 0.4 );
    correctedResponseEta_mgf[i]->SetMinimum( 0.4 );

    correctedResponsePt_mgf [i]->SetMaximum( 1.6 );
    correctedResponseEta_mgf[i]->SetMaximum( 1.6 );
  }

  TF1 *fL5 = new TF1("fL5", "[0]+[1]*log10(x)+[2]*log10(x)*log10(x)", 10., 2000.);
  // paramters from CondFormats/JetMETObjects/data/L5Flavor_IC5.txt in V01-08-04
  //fL5->SetParameters(0.692, 0.176, -0.026); // dijet
  fL5->SetParameters(0.599, 0.176, -0.006); // ttbar
  fL5->SetLineWidth(1);

  // create legends

  TLegend* legend = new TLegend(0.7, 0.75, 0.91, 0.85);
  legend->SetFillColor(0);
  legend->AddEntry(invMass[0],"before fit");
  legend->AddEntry(invMass[1],"after fit");

  TLegend* legendCorrFacs = new TLegend(0.8, 0.76, 0.91, 0.85);
  legendCorrFacs->SetFillColor(0);
  for(unsigned i=0; i<5; i++) {
    if(i!=0 && i!=3)
      legendCorrFacs->AddEntry(corrFacsPt_mgf[i], corrLevels[i+1]);
  }

  TLegend* legendCorrResponse = new TLegend(0.8, 0.73, 0.91, 0.85);
  legendCorrResponse->SetFillColor(0);
  for(unsigned i=0; i<6; i++) {
    if(i!=1 && i!=4)
      legendCorrResponse->AddEntry(correctedResponsePt_mgf[i], corrLevels[i]);
  }

  // create a line

  TLine* line = new TLine();
  line->SetLineStyle(2);

  // draw hists

  c->cd(1);

  scale->Draw();
  c->Draw();
  if(individualPdf) c->Print("top_scale.pdf");
  ps->NewPage();

  weight->Draw();
  c->Draw();
  if(individualPdf) c->Print("top_weight.pdf");
  ps->NewPage();

  truth->Draw();
  c->Draw();
  if(individualPdf) c->Print("top_truth.pdf");
  ps->NewPage();

  pt->Draw();
  c->Draw();
  if(individualPdf) c->Print("top_pt.pdf");
  ps->NewPage();

  eta->Draw();
  c->Draw();
  if(individualPdf) c->Print("top_eta.pdf");
  ps->NewPage();

  phi->Draw();
  c->Draw();
  if(individualPdf) c->Print("top_phi.pdf");
  ps->NewPage();

  genPt->Draw();
  c->Draw();
  if(individualPdf) c->Print("top_genPt.pdf");
  ps->NewPage();

  recPt12->Draw("box");
  c->Draw();
  if(individualPdf) c->Print("top_recPt12.pdf");
  ps->NewPage();

  genPt12->Draw("box");
  c->Draw();
  if(individualPdf) c->Print("top_genPt12.pdf");
  ps->NewPage();

  eta12->Draw("box");
  c->Draw();
  if(individualPdf) c->Print("top_eta12.pdf");
  ps->NewPage();

  deltaRecPt->Draw();
  c->Draw();
  if(individualPdf) c->Print("top_deltaRecPt.pdf");
  ps->NewPage();

  deltaGenPt->Draw();
  c->Draw();
  if(individualPdf) c->Print("top_deltaGenPt.pdf");
  ps->NewPage();

  deltaEta->Draw();
  c->Draw();
  if(individualPdf) c->Print("top_deltaEta.pdf");
  ps->NewPage();

  deltaRecPtRecPt1->Draw("box");
  c->Draw();
  if(individualPdf) c->Print("top_deltaRecPtRecPt1.pdf");
  ps->NewPage();

  deltaRecPtEta1->Draw("box");
  c->Draw();
  if(individualPdf) c->Print("top_deltaRecPtEta1.pdf");
  ps->NewPage();

  meanPt->Draw();
  c->Draw();
  if(individualPdf) c->Print("top_meanPt.pdf");
  ps->NewPage();

  meanEta->Draw();
  c->Draw();
  if(individualPdf) c->Print("top_meanEta.pdf");
  ps->NewPage();

  wPt->Draw();
  c->Draw();
  if(individualPdf) c->Print("top_wPt.pdf");
  ps->NewPage();

  wEta->Draw();
  c->Draw();
  if(individualPdf) c->Print("top_wEta.pdf");
  ps->NewPage();

  dR->Draw();
  c->Draw();
  if(individualPdf) c->Print("top_dR.pdf");
  ps->NewPage();

  dRrecPt->Draw("box");
  c->Draw();
  if(individualPdf) c->Print("top_dRrecPt.pdf");
  ps->NewPage();

  dRgenPt->Draw("box");
  c->Draw();
  if(individualPdf) c->Print("top_dRgenPt.pdf");
  ps->NewPage();

  dReta->Draw("box");
  c->Draw();
  if(individualPdf) c->Print("top_dReta.pdf");
  ps->NewPage();

  dRwPt->Draw("box");
  c->Draw();
  if(individualPdf) c->Print("top_dRwPt.pdf");
  ps->NewPage();

  invMass [0]->Draw("p");
  fInvMass[0]->Draw("same");
  invMass [1]->Draw("p same");
  fInvMass[1]->Draw("same");  
  legend->Draw("same");
  paveText[0]->Draw();
  paveText[1]->Draw();
  c->Draw();
  if(individualPdf) c->Print("top_invMass.pdf");
  ps->NewPage();

  messTruth [0]->Draw("p");
  fMessTruth[0]->Draw("same");
  messTruth [1]->Draw("p same");
  fMessTruth[1]->Draw("same");  
  legend->Draw("same");
  c->Draw();
  if(individualPdf) c->Print("top_messTruth.pdf");
  ps->NewPage();

  TString appendix_proj[4] = {"mean","sigma","mgf","sgf"};

  double yLine[4] = {1., 0., 1., 0.};

  for(unsigned j=0; j<4; j++) {
    c->SetLogx(1);
    messTruthPtGen_proj[0][j]->Draw("p");
    messTruthPtGen_proj[1][j]->Draw("p same");
    legend->Draw("same");
    line->DrawLine(messTruthPtGen_proj[0][j]->GetXaxis()->GetXmin(), yLine[j],
		   messTruthPtGen_proj[0][j]->GetXaxis()->GetXmax(), yLine[j]);
    c->Draw();
    if(individualPdf) c->Print("top_messTruthPtGen_"+appendix_proj[j]+".pdf");
    ps->NewPage();
    c->SetLogx(0);
  }

  for(unsigned j=0; j<4; j++) {
    c->SetLogx(1);
    messTruthPtRec_proj[0][j]->Draw("p");
    messTruthPtRec_proj[1][j]->Draw("p same");
    legend->Draw("same");
    line->DrawLine(messTruthPtRec_proj[0][j]->GetXaxis()->GetXmin(), yLine[j],
		   messTruthPtRec_proj[0][j]->GetXaxis()->GetXmax(), yLine[j]);
    c->Draw();
    if(individualPdf) c->Print("top_messTruthPtRec_"+appendix_proj[j]+".pdf");
    ps->NewPage();
    c->SetLogx(0);
  }

  for(unsigned j=0; j<4; j++) {
    messTruthEta_proj[0][j]->Draw("p");
    messTruthEta_proj[1][j]->Draw("p same");
    legend->Draw("same");
    line->DrawLine(messTruthEta_proj[0][j]->GetXaxis()->GetXmin(), yLine[j],
		   messTruthEta_proj[0][j]->GetXaxis()->GetXmax(), yLine[j]);
    c->Draw();
    if(individualPdf) c->Print("top_messTruthEta_"+appendix_proj[j]+".pdf");
    ps->NewPage();
  }

  for(unsigned j=0; j<4; j++) {
    c->SetLogx(1);
    responsePtGen_proj[0][j]->Draw("p");
    responsePtGen_proj[1][j]->Draw("p same");
    legend->Draw("same");
    line->DrawLine(responsePtGen_proj[0][j]->GetXaxis()->GetXmin(), yLine[j],
		   responsePtGen_proj[0][j]->GetXaxis()->GetXmax(), yLine[j]);
    c->Draw();
    if(individualPdf) c->Print("top_responsePtGen_"+appendix_proj[j]+".pdf");
    ps->NewPage();
    c->SetLogx(0);
  }

  for(unsigned j=0; j<4; j++) {
    c->SetLogx(1);
    responsePtRec_proj[0][j]->Draw("p");
    responsePtRec_proj[1][j]->Draw("p same");
    legend->Draw("same");
    line->DrawLine(responsePtRec_proj[0][j]->GetXaxis()->GetXmin(), yLine[j],
		   responsePtRec_proj[0][j]->GetXaxis()->GetXmax(), yLine[j]);
    c->Draw();
    if(individualPdf) c->Print("top_responsePtRec_"+appendix_proj[j]+".pdf");
    ps->NewPage();
    c->SetLogx(0);
  }

  for(unsigned j=0; j<4; j++) {
    responseEta_proj[0][j]->Draw("p");
    responseEta_proj[1][j]->Draw("p same");
    legend->Draw("same");
    line->DrawLine(responseEta_proj[0][j]->GetXaxis()->GetXmin(), yLine[j],
		   responseEta_proj[0][j]->GetXaxis()->GetXmax(), yLine[j]);
    c->Draw();
    if(individualPdf) c->Print("top_responseEta_"+appendix_proj[j]+".pdf");
    ps->NewPage();
  }

  c->SetLogx(1);
  corrFacsPt_mgf[1]->Draw("p");
  line->DrawLine(corrFacsPt_mgf[0]->GetXaxis()->GetXmin(), 1.,
		 corrFacsPt_mgf[0]->GetXaxis()->GetXmax(), 1.);
  fL5->Draw("same");
  for(unsigned int i=1; i<5; i++) {
    if(i!=0 && i!=3)
      corrFacsPt_mgf[i]->Draw("p same");
  }
  legendCorrFacs->Draw("same");
  c->Draw();
  if(individualPdf) c->Print("top_corrFacsPt.pdf");
  ps->NewPage();

  c->SetLogx(0);
  corrFacsEta_mgf[1]->Draw("p");
  for(unsigned int i=1; i<5; i++) {
    if(i!=0 && i!=3)
      corrFacsEta_mgf[i]->Draw("p same");
  }
  legendCorrFacs->Draw("same");
  line->DrawLine(corrFacsEta_mgf[0]->GetXaxis()->GetXmin(), 1.,
		 corrFacsEta_mgf[0]->GetXaxis()->GetXmax(), 1.);
  c->Draw();
  if(individualPdf) c->Print("top_corrFacsEta.pdf");
  ps->NewPage();

  c->SetLogx(1);
  correctedResponsePt_mgf[0]->Draw("p");
  for(unsigned int i=1; i<6; i++) {
    if(i!=1 && i!=4)
      correctedResponsePt_mgf[i]->Draw("p same");
  }
  legendCorrResponse->Draw("same");
  line->DrawLine(correctedResponsePt_mgf[0]->GetXaxis()->GetXmin(), 1.,
		 correctedResponsePt_mgf[0]->GetXaxis()->GetXmax(), 1.);
  c->Draw();
  if(individualPdf) c->Print("top_correctedResponsePt.pdf");
  ps->NewPage();

  c->SetLogx(0);
  correctedResponseEta_mgf[0]->Draw("p");
  for(unsigned int i=1; i<6; i++) {
    if(i!=1 && i!=4)
      correctedResponseEta_mgf[i]->Draw("p same");
  }
  legendCorrResponse->Draw("same");
  line->DrawLine(correctedResponseEta_mgf[0]->GetXaxis()->GetXmin(), 1.,
		 correctedResponseEta_mgf[0]->GetXaxis()->GetXmax(), 1.);
  c->Draw();
  if(individualPdf) c->Print("top_correctedResponseEta.pdf");

  ps->Close();

  // write to root-file

  objToBeWritten.push_back( scale  );
  objToBeWritten.push_back( weight );
  objToBeWritten.push_back( truth  );
  objToBeWritten.push_back( pt  );
  objToBeWritten.push_back( eta );
  objToBeWritten.push_back( phi );
  objToBeWritten.push_back( genPt  );
  objToBeWritten.push_back( recPt12 );
  objToBeWritten.push_back( genPt12 );
  objToBeWritten.push_back( eta12 );
  objToBeWritten.push_back( deltaRecPt );
  objToBeWritten.push_back( deltaGenPt );
  objToBeWritten.push_back( deltaEta );
  objToBeWritten.push_back( deltaRecPtRecPt1 );
  objToBeWritten.push_back( deltaRecPtEta1 );
  objToBeWritten.push_back( meanPt  );
  objToBeWritten.push_back( meanEta );
  objToBeWritten.push_back( wPt  );
  objToBeWritten.push_back( wEta );
  objToBeWritten.push_back( dR   );
  objToBeWritten.push_back( dRrecPt );
  objToBeWritten.push_back( dRgenPt );
  objToBeWritten.push_back( dReta   );
  objToBeWritten.push_back( dRwPt   );

  for(unsigned a=0; a<2; a++){
    objToBeWritten.push_back( invMass  [a] );
    objToBeWritten.push_back( messTruth[a] );
    objToBeWritten.push_back( messTruthPtGen[a] );
    objToBeWritten.push_back( messTruthPtRec[a] );
    objToBeWritten.push_back( messTruthEta  [a] );
    objToBeWritten.push_back( responsePtGen [a] );
    objToBeWritten.push_back( responsePtRec [a] );
    objToBeWritten.push_back( responseEta   [a] );
    for(unsigned j=0; j<4; j++) {
      objToBeWritten.push_back( messTruthPtGen_proj[a][j] );
      objToBeWritten.push_back( messTruthPtRec_proj[a][j] );
      objToBeWritten.push_back( messTruthEta_proj  [a][j] );
      objToBeWritten.push_back( responsePtGen_proj [a][j] );
      objToBeWritten.push_back( responsePtRec_proj [a][j] );
      objToBeWritten.push_back( responseEta_proj   [a][j] );
    }
  }
  for(unsigned i=0; i<5; i++){
    objToBeWritten.push_back( corrFacsPt [i] );
    objToBeWritten.push_back( corrFacsEta[i] );
    objToBeWritten.push_back( corrFacsPt_mgf [i] );
    objToBeWritten.push_back( corrFacsEta_mgf[i] );
  }
  for(unsigned i=0; i<6; i++){
    objToBeWritten.push_back( correctedResponsePt [i] );
    objToBeWritten.push_back( correctedResponseEta[i] );
    objToBeWritten.push_back( correctedResponsePt_mgf [i] );
    objToBeWritten.push_back( correctedResponseEta_mgf[i] );
  }

  if( outputROOT_ ) writeToRootFile( objToBeWritten, "Top" );

  // clean up

  delete scale;
  delete weight;
  delete truth;
  delete pt;
  delete eta;
  delete phi;
  delete genPt;
  delete recPt12;
  delete genPt12;
  delete eta12;
  delete deltaRecPt;
  delete deltaGenPt;
  delete deltaEta;
  delete deltaRecPtRecPt1;
  delete deltaRecPtEta1;
  delete meanPt;
  delete meanEta;
  delete wPt;
  delete wEta;
  delete dR;
  delete dRrecPt;
  delete dRgenPt;
  delete dReta;
  delete dRwPt;
  for(unsigned a=0; a<2; a++){
    delete paveText [a];
    delete invMass  [a];
    delete messTruth[a];
    delete messTruthPtGen[a];
    delete messTruthPtRec[a];
    delete messTruthEta  [a];
    delete responsePtGen [a];
    delete responsePtRec [a];
    delete responseEta   [a];
    for(unsigned j=0; j<4; j++) {
      delete messTruthPtGen_proj[a][j];
      delete messTruthPtRec_proj[a][j];
      delete messTruthEta_proj  [a][j];
      delete responsePtGen_proj [a][j];
      delete responsePtRec_proj [a][j];
      delete responseEta_proj   [a][j];
    }
  }
  for(unsigned i=0; i<5; i++){
    delete corrFacsPt [i];
    delete corrFacsEta[i];
    delete corrFacsPt_mgf [i];
    delete corrFacsEta_mgf[i];
  }
  for(unsigned i=0; i<6; i++){
    delete correctedResponsePt [i];
    delete correctedResponseEta[i];
    delete correctedResponsePt_mgf [i];
    delete correctedResponseEta_mgf[i];
  }

  delete fL5;

  delete legend;
  delete line;

}



//---------------------------------------------------------------
//!   \brief Vary parameters around fitted value and plot
//!          chi2 profile.
//!
//!   The method generates one a 1-dim plot
//!   for each parameter in which that parameter is varied
//!   and all other parameters are left as in \p TParameters::GetPars().
//!   Additionally, correlations between parameters are
//!   shown in 2-dim plots, where two parameters are varied.
//!   So far, only correlations between the 1st and 2nd, the
//!   3rd and 4th and so on are plotted. This has to be extended
//!   in a clever way...
//!
//!   \note Chi2 is calculated with that scheme for scaling of
//!   residuals that was used in the last fit iteration
//---------------------------------------------------------------
void TControlPlots::makeControlPlotsParameterScan()
{
  std::vector<TObject*> objToBeWritten;

  TCanvas * const c1 = new TCanvas("c1","",600,600);
  TPostScript * const ps = new TPostScript("controlplotsParameterScan.ps",111);


  // Loop over parameters and vary one parameter
  // while the others are left as in TParameter::GetPars()
  TGraph *gParScan[par_->GetNumberOfParameters()];
  TH1F *hFrame[par_->GetNumberOfParameters()];
  TLine *line[par_->GetNumberOfParameters()];
  for(int i = 0; i < par_->GetNumberOfParameters(); i++)
    {
      // Store original value of parameter i
      double origPar = par_->GetPars()[i];

      // Vary parameter and get chi2
      double x_param[21];
      double y_chi2[21];
      for(int a = 0; a < 21; a++)
	{
	  x_param[a] = 0.;
	  y_chi2[a] = 0.;
	}
      for(int a = 0; a < 21; a++)
	{
	  double variedPar = origPar - 0.1 + 0.01*a;
	  x_param[a] = variedPar;
	  par_->GetPars()[i] = variedPar;
	  for( std::vector<TData*>::const_iterator it = data_->begin();  it < data_->end();  ++it )
	    {
	      y_chi2[a] += (*it)->chi2();
	    }
	}

      // Reset original parameter
      par_->GetPars()[i] = origPar;


      // Draw graphs
      TString name = "gParScan";
      name += i;
      TString title = "Parameter ";
      title += i;
      if( i < par_->GetNumberOfTowerParameters() )
	{
	  title += " (tower parameter ";
	  title += i;
	}
      else
	{
	  title += " (jet parameter ";
	  title += i - (par_->GetNumberOfTowerParameters());
	}
      title += ");Parameter p_{";
      title += i;
      title += "};#chi^{2}";
      gParScan[i] = new TGraph(21,x_param,y_chi2);
      gParScan[i]->SetName(name);
      gParScan[i]->SetTitle(title);
      gParScan[i]->SetMarkerStyle(20);
      objToBeWritten.push_back(gParScan[i]);

      name = "hParScanFrame";
      name += i;
      hFrame[i] = new TH1F(name,title,1,x_param[0]-0.01,x_param[20]+0.01);
      double y_min = *min_element(y_chi2,y_chi2+21);
      double y_max = *max_element(y_chi2,y_chi2+21);
      double range = y_max - y_min;
      hFrame[i]->GetYaxis()->SetRangeUser(y_min - range/20, y_max + range/20);

      line[i] = new TLine(origPar,hFrame[i]->GetMinimum(),origPar,hFrame[i]->GetMaximum());
      line[i]->SetLineStyle(2);
      line[i]->SetLineWidth(1);

      c1->cd();
      hFrame[i]->Draw();
      gParScan[i]->Draw("Psame");
      line[i]->Draw("same");
      c1->Draw();   
    }



  // Correlations between 2 parameters
  int nPlots = par_->GetNumberOfParameters()/2;
  TGraph2D *gParScan2D[nPlots];
  for(int i = 0; i < nPlots; i++)
    {
      int parIdx = 2*i;

      // Store original value of parameters i and i+1
      double origParX = par_->GetPars()[parIdx];
      double origParY = par_->GetPars()[parIdx+1];

      // Vary parameter and get chi2
      double x_param[441];
      double y_param[441];
      double z_chi2[441];
      for(int a = 0; a < 441; a++)
	{
	  x_param[a] = 0.;
	  y_param[a] = 0.;
	  z_chi2[a] = 0.;
	}
      for(int a = 0; a < 21; a++)
	{
	  double variedParX = origParX - 0.1 + 0.01*a;
	  par_->GetPars()[parIdx] = variedParX;
	  for(int b = 0; b < 21; b++)
	    {
	      double variedParY = origParY - 0.1 + 0.01*b;
	      par_->GetPars()[parIdx+1] = variedParY;

	      int pointIdx = b + 21*a;
	      x_param[pointIdx] = variedParX;
	      y_param[pointIdx] = variedParY;

	      for( std::vector<TData*>::const_iterator it = data_->begin();  it < data_->end();  ++it )
		{
		  z_chi2[pointIdx] += (*it)->chi2();
		}
	    }
	}

      // Reset original parameter
      par_->GetPars()[parIdx] = origParX;
      par_->GetPars()[parIdx+1] = origParY;


      // Draw graphs
      TString name = "gParScan2D";
      name += i;
      TString title = "Parameter ";
      title += parIdx;
      title += " and ";
      title += parIdx+1;
      title += ";Parameter p_{";
      title += parIdx;
      title += "};Parameter p_{";
      title += parIdx+1;
      title += "};#chi^{2}";
      gParScan2D[i] = new TGraph2D(441,x_param,y_param,z_chi2);
      gParScan2D[i]->SetName(name);
      gParScan2D[i]->SetTitle(title);
      objToBeWritten.push_back(gParScan2D[i]);

      c1->cd();
      gParScan2D[i]->Draw("cont3");
      c1->SetGrid();
      c1->Draw();   
      if( i < nPlots-1 ) ps->NewPage();
    }




  // Clean up
  ps->Close();

  if( outputROOT_ ) writeToRootFile( objToBeWritten, "ParScan" );

  delete c1;
  delete ps;

  for(int i = 0; i < par_->GetNumberOfParameters(); i++)
    {
      delete gParScan[i];
      delete hFrame[i];
      delete line[i];
    }
  for(int i = 0; i < nPlots; i++)
    {
      delete gParScan2D[i];
    }
}



//!  \brief Control plots for events of type \p TwoJetsPtBalanceEvent
// -------------------------------------------------------------
void TControlPlots::makeControlPlotsTwoJetsPtBalance() {
  bool debug = false;

  if( debug ) std::cout << "Entering makeControlPlotsTwoJetsPtBalance()\n";

  std::vector<TObject*>                objToBeWritten;
  std::vector<TData*>::const_iterator  data_it;

  // -- Create 2D histograms of balance etc vs eta, ptdijet --------------

  // Create a ptdijet and eta binning
  // ptdijet = x, eta = y
  std::vector<double> binEdgesPtDijet = bag_of<double>(config_->read<std::string>("Control plots pt bin edges",""));

  // Eta bins
  std::vector<double> binEdgesEta = bag_of<double>(config_->read<std::string>("Control plots eta bin edges",""));
  Binning bins(binEdgesPtDijet,binEdgesEta);
  //  bins.Print();

  if( debug ) std::cout << "Creating 2D histograms\n";

  // Create 2D histograms of balance etc vs eta,
  // ptDijet, ptGenDijet in bins of eta, ptDijet,
  // ptGenDijet
  // Outer index: x-axis quantity
  // Inner index: binning
  //  0: Balance vs eta, bins of ptDijet
  //  1: Balance vs eta, bins of ptDijet(gen)
  //  2: Balance vs ptDijet, bins of eta
  //  3: Balance vs ptDijet(gen), bins of eta
  std::vector< std::vector<TH2F*> > h2BUncorr(4);    // Uncorrected dijet balance
  std::vector< std::vector<TH2F*> > h2BCorr(4);      // Dijet balance corrected by kalibri fit
  std::vector< std::vector<TH2F*> > h2BCorrL2L3(4);  // Dijet balance corrected by JetMET L2L3 correction

  // Create 2D histograms of response vs eta in
  // bins of ptGen
  // 0: Response vs eta, bins of ptGen
  // 1: Response vs ptGen, bins of eta
  std::vector< std::vector<TH2F*> > h2RUncorr(2);    // Uncorrected response
  std::vector< std::vector<TH2F*> > h2RCorr(2);      // Response corrected by kalibri fit
  std::vector< std::vector<TH2F*> > h2RCorrL2L3(2);  // Response corrected by JetMET L2L3 correction

  TString ptSumLabel = "2 #upoint (p^{1}_{T} - p^{2}_{T}) / (p^{1}_{T} + p^{2}_{T})";
  TString respLabel = "p_{T} / p_{T,gen}";

  // Loop over ptDijet bins
  for(int ptbin = 0; ptbin < bins.nBinsX(); ptbin++) {
    char name[50];
    sprintf(name,"h2BUncorrVsEta_ptDijet%i",ptbin);
    TH2F * h2 = new TH2F(name,";#eta;"+ptSumLabel,21,-5,5,51,-2,2);
    h2BUncorr.at(0).push_back(h2);
    
    sprintf(name,"h2BCorrVsEta_ptDijet%i",ptbin);
    h2BCorr.at(0).push_back(static_cast<TH2F*>(h2->Clone(name)));
    
    sprintf(name,"h2BCorrL2L3VsEta_ptDijet%i",ptbin);
    h2BCorrL2L3.at(0).push_back(static_cast<TH2F*>(h2->Clone(name)));


    sprintf(name,"h2BUncorrVsEta_ptDijetGen%i",ptbin);
    h2 = new TH2F(name,";#eta;"+ptSumLabel,21,-5,5,51,-2,2);
    h2BUncorr.at(1).push_back(h2);
    
    sprintf(name,"h2BCorrVsEta_ptDijetGen%i",ptbin);
    h2BCorr.at(1).push_back(static_cast<TH2F*>(h2->Clone(name)));
    
    sprintf(name,"h2BCorrL2L3VsEta_ptDijetGen%i",ptbin);
    h2BCorrL2L3.at(1).push_back(static_cast<TH2F*>(h2->Clone(name)));


    sprintf(name,"h2RUncorrVsEta_ptGen%i",ptbin);
    h2 = new TH2F(name,";#eta;"+respLabel,21,-5,5,51,0,2);
    h2RUncorr.at(0).push_back(h2);
    
    sprintf(name,"h2RCorrVsEta_ptGen%i",ptbin);
    h2RCorr.at(0).push_back(static_cast<TH2F*>(h2->Clone(name)));
    
    sprintf(name,"h2RCorrL2L3VsEta_ptGen%i",ptbin);
    h2RCorrL2L3.at(0).push_back(static_cast<TH2F*>(h2->Clone(name)));
  } // End of loop over ptDijet bins


  // Logarithmic binning
  const int nLogBins = 15;
  double logBins[nLogBins+1];
  equidistLogBins(logBins,nLogBins,bins.xLow(0),bins.xUp(bins.nBinsX()-1));
  // Loop over eta bins
  for(int etaBin = 0; etaBin < bins.nBinsY(); etaBin++) {
    char name[50];
    sprintf(name,"h2BUncorrVsPtDijet_eta%i",etaBin);
    TH2F * h2 = new TH2F(name,";p^{dijet}_{T} (GeV);"+ptSumLabel,nLogBins,logBins,51,-2,2);
    h2BUncorr.at(2).push_back(h2);
    
    sprintf(name,"h2BCorrVsPtDijet_eta%i",etaBin);
    h2BCorr.at(2).push_back(static_cast<TH2F*>(h2->Clone(name)));
    
    sprintf(name,"h2BCorrL2L3VsPtDijet_eta%i",etaBin);
    h2BCorrL2L3.at(2).push_back(static_cast<TH2F*>(h2->Clone(name)));


    sprintf(name,"h2BUncorrVsPtDijetGen_eta%i",etaBin);
    h2 = new TH2F(name,";p^{dijet}_{T,gen} (GeV);"+ptSumLabel,nLogBins,logBins,51,-2,2);
    h2BUncorr.at(3).push_back(h2);
    
    sprintf(name,"h2BCorrVsPtDijetGen_eta%i",etaBin);
    h2BCorr.at(3).push_back(static_cast<TH2F*>(h2->Clone(name)));
    
    sprintf(name,"h2BCorrL2L3VsPtDijetGen_eta%i",etaBin);
    h2BCorrL2L3.at(3).push_back(static_cast<TH2F*>(h2->Clone(name)));


    sprintf(name,"h2RUncorrVsPtGen_eta%i",etaBin);
    h2 = new TH2F(name,";p^{gen}_{T} (GeV);"+respLabel,nLogBins,logBins,51,0,2);
    h2RUncorr.at(1).push_back(h2);
    
    sprintf(name,"h2RCorrVsPtGen_eta%i",etaBin);
    h2RCorr.at(1).push_back(static_cast<TH2F*>(h2->Clone(name)));
    
    sprintf(name,"h2RCorrL2L3VsPtGen_eta%i",etaBin);
    h2RCorrL2L3.at(1).push_back(static_cast<TH2F*>(h2->Clone(name)));
  } // End of loop over eta bins

  // Other control quantities
  TH2F * h2Eta = new TH2F("h2Eta",";#eta^{1};#eta^{2}",21,-5.,5.,20,-5.,5.);
  objToBeWritten.push_back(h2Eta);

  TH1F * hPtDijetSpectrum = new TH1F("hPtDijetSpectrum",";p^{dijet}_{T} (GeV);dN / dp^{dijet}_{T}  1 / (GeV)",50,bins.xLow(0),bins.xUp(bins.nBinsX()-1));
  hPtDijetSpectrum->Sumw2();
  hPtDijetSpectrum->SetLineWidth(2);
  objToBeWritten.push_back(hPtDijetSpectrum);

  
  if( debug ) std::cout << "Filling 2D histograms\n";

  // First loop over data and fill balance uncorrected
  // and Kalibri corrected data into the corresponding
  // 2D histogram
  for( data_it = data_->begin(); data_it != data_->end(); data_it++ ) {
    TwoJetsPtBalanceEvent *evt = dynamic_cast<TwoJetsPtBalanceEvent*>(*data_it);
    if( evt ) {
      if( evt->flaggedBad() ) continue;   // Discard events flagged bad
      double weight = evt->GetWeight();

      h2Eta->Fill(evt->getJet1()->eta(),evt->getJet2()->eta(),evt->GetWeight());
      hPtDijetSpectrum->Fill(evt->ptDijet(),evt->GetWeight());

      // Find ptDijet bins
      int ptDijetBin = bins.iX(evt->ptDijet());
      int ptDijetGenBin = bins.iX(evt->ptDijetGen());
      
      for(int i = 0; i < 2; i++) {
	Jet * jet = evt->getJet1();
	if( i == 1 ) evt->getJet2();

	if( 0 <= ptDijetBin && ptDijetBin < bins.nBinsX() ) {
	  h2BUncorr.at(0).at(ptDijetBin)->Fill(jet->eta(),evt->ptBalance(),weight);
	  h2BCorr.at(0).at(ptDijetBin)->Fill(jet->eta(),evt->ptBalanceCorr(),weight);
	  h2BCorrL2L3.at(0).at(ptDijetBin)->Fill(jet->eta(),evt->ptBalanceCorrL2L3(),weight);
	}

	if( 0 <= ptDijetGenBin && ptDijetGenBin < bins.nBinsX() ) {
	  h2BUncorr.at(1).at(ptDijetGenBin)->Fill(jet->eta(),evt->ptBalance(),weight);
	  h2BCorr.at(1).at(ptDijetGenBin)->Fill(jet->eta(),evt->ptBalanceCorr(),weight);
	  h2BCorrL2L3.at(1).at(ptDijetGenBin)->Fill(jet->eta(),evt->ptBalanceCorrL2L3(),weight);
	}

	// Find ptGen bin
	int ptGenBin = bins.iX(jet->GenPt());
	if( 0 <= ptGenBin && ptGenBin < bins.nBinsX() ) {
	  h2RUncorr.at(0).at(ptGenBin)->Fill(jet->eta(),jet->Et()/jet->GenPt(),weight);
	  h2RCorr.at(0).at(ptGenBin)->Fill(jet->eta(),jet->correctedEt()/jet->GenPt(),weight);
	  h2RCorrL2L3.at(0).at(ptGenBin)->Fill(jet->eta(),jet->corFactors.getL2L3()*jet->Et()/jet->GenPt(),weight);
	}
	
	// Find eta bin
	int etaBin = bins.iY(jet->eta());
	if( 0 <= etaBin && etaBin < bins.nBinsY() ) {
	  h2BUncorr.at(2).at(etaBin)->Fill(evt->ptDijet(),evt->ptBalance(),weight);
	  h2BCorr.at(2).at(etaBin)->Fill(evt->ptDijet(),evt->ptBalanceCorr(),weight);
	  h2BCorrL2L3.at(2).at(etaBin)->Fill(evt->ptDijet(),evt->ptBalanceCorrL2L3(),weight);

	  h2BUncorr.at(3).at(etaBin)->Fill(evt->ptDijetGen(),evt->ptBalance(),weight);
	  h2BCorr.at(3).at(etaBin)->Fill(evt->ptDijetGen(),evt->ptBalanceCorr(),weight);
	  h2BCorrL2L3.at(3).at(etaBin)->Fill(evt->ptDijetGen(),evt->ptBalanceCorrL2L3(),weight);

	  h2RUncorr.at(1).at(etaBin)->Fill(jet->GenPt(),jet->Et()/jet->GenPt(),weight);
	  h2RCorr.at(1).at(etaBin)->Fill(jet->GenPt(),jet->correctedEt()/jet->GenPt(),weight);
	  h2RCorrL2L3.at(1).at(etaBin)->Fill(jet->GenPt(),jet->corFactors.getL2L3()*jet->Et()/jet->GenPt(),weight);
	}
      }
    }
  } // End of first loop over data

  if( debug ) std::cout << "Projecting 1D histograms\n";

  // -- Balance vs eta, pt ---------------------------------
  // Get 1D projections from 2D histograms
  //  0: Balance vs eta, bins of ptDijet
  //  1: Balance vs eta, bins of ptDijet(gen)
  //  2: Balance vs ptDijet, bins of eta
  //  3: Balance vs ptDijet(gen), bins of eta
  std::vector< std::vector<TH1F*> > hBUncorr(4);
  std::vector< std::vector<TH1F*> > hBCorr(4);
  std::vector< std::vector<TH1F*> > hBCorrL2L3(4);

  std::vector< std::vector< std::vector<TH1F*> > > hBDistUncorr(4);
  std::vector< std::vector< std::vector<TH1F*> > > hBDistCorr(4);
  std::vector< std::vector< std::vector<TH1F*> > > hBDistCorrL2L3(4);

  // 0: Response vs eta, bins of ptGen
  // 1: Response vs ptGen, bins of eta
  std::vector< std::vector<TH1F*> > hRUncorr(2);
  std::vector< std::vector<TH1F*> > hRCorr(2);
  std::vector< std::vector<TH1F*> > hRCorrL2L3(2);

  for(int i = 0; i < 4; i++) {
    int nBins = bins.nBinsX();
    if( i > 1 ) nBins = bins.nBinsY();

    hBUncorr.at(i) = std::vector<TH1F*>(nBins);
    hBCorr.at(i) = std::vector<TH1F*>(nBins);
    hBCorrL2L3.at(i) = std::vector<TH1F*>(nBins);
    
    hBDistUncorr.at(i) = std::vector< std::vector<TH1F*> >(nBins);
    hBDistCorr.at(i) = std::vector< std::vector<TH1F*> >(nBins);
    hBDistCorrL2L3.at(i) = std::vector< std::vector<TH1F*> >(nBins);

    if( i < 2 ) {
      if( i == 1 ) nBins = bins.nBinsY();
      hRUncorr.at(i) = std::vector<TH1F*>(nBins);
      hRCorr.at(i) = std::vector<TH1F*>(nBins);
      hRCorrL2L3.at(i) = std::vector<TH1F*>(nBins);
    }
  }
  int colorUncorr = 1;
  int colorCorr = 2;
  int colorGen = 4;

  if( debug ) std::cout << " Loop over ptDijet bins\n";

  for(int ptDijetBin = 0; ptDijetBin < bins.nBinsX(); ptDijetBin++) { // Loop over pt bins
    if( debug ) {
      std::cout << "  Bin " << ptDijetBin << ": ";
      std::cout << bins.xLow(bins.bin(ptDijetBin,0)) << " - ";
      std::cout << bins.xUp(bins.bin(ptDijetBin,0)) << std::endl;
    }
    char title[100];
    sprintf(title,"%.0f < p^{dijet}_{T} < %.0f GeV",bins.xLow(bins.bin(ptDijetBin,0)),bins.xUp(bins.bin(ptDijetBin,0)));

    if( debug ) std::cout << "   Uncorrected balance vs eta, bins of ptDijet\n";
    // Uncorrected balance vs eta, bins of ptDijet
    fit2DMean(h2BUncorr.at(0).at(ptDijetBin),
		   hBUncorr.at(0).at(ptDijetBin),
		   hBDistUncorr.at(0).at(ptDijetBin),colorUncorr);
    hBUncorr.at(0).at(ptDijetBin)->SetTitle(title);    
    objToBeWritten.push_back(hBUncorr.at(0).at(ptDijetBin));

    if( debug ) std::cout << "   Corrected balance vs eta, bins of ptDijet\n";
    // Corrected balance vs eta, bins of ptDijet
    fit2DMean(h2BCorr.at(0).at(ptDijetBin),
		   hBCorr.at(0).at(ptDijetBin),
		   hBDistCorr.at(0).at(ptDijetBin),colorCorr);
    hBCorr.at(0).at(ptDijetBin)->SetTitle(title);    
    objToBeWritten.push_back(hBCorr.at(0).at(ptDijetBin));

    if( debug ) std::cout << "   Genjet balance vs eta, bins of ptDijet\n";
    // Genjet balance vs eta, bins of ptDijet
    fit2DMean(h2BCorrL2L3.at(0).at(ptDijetBin),
		   hBCorrL2L3.at(0).at(ptDijetBin),
		   hBDistCorrL2L3.at(0).at(ptDijetBin),colorGen);
    hBCorrL2L3.at(0).at(ptDijetBin)->SetTitle(title);    
    objToBeWritten.push_back(hBCorrL2L3.at(0).at(ptDijetBin));

    sprintf(title,"%.0f < p^{dijet}_{T,gen} < %.0f GeV",bins.xLow(bins.bin(ptDijetBin,0)),bins.xUp(bins.bin(ptDijetBin,0)));
    if( debug ) std::cout << "   Uncorrected balance vs eta, bins of ptDijet\n";
    // Uncorrected balance vs eta, bins of ptDijetGen
    fit2DMean(h2BUncorr.at(1).at(ptDijetBin),
		   hBUncorr.at(1).at(ptDijetBin),
		   hBDistUncorr.at(1).at(ptDijetBin),colorUncorr);
    hBUncorr.at(1).at(ptDijetBin)->SetTitle(title);    
    objToBeWritten.push_back(hBUncorr.at(1).at(ptDijetBin));

    if( debug ) std::cout << "   Corrected balance vs eta, bins of ptDijet\n";
    // Corrected balance vs eta, bins of ptDijetGen
    fit2DMean(h2BCorr.at(1).at(ptDijetBin),
		   hBCorr.at(1).at(ptDijetBin),
		   hBDistCorr.at(1).at(ptDijetBin),colorCorr);
    hBCorr.at(1).at(ptDijetBin)->SetTitle(title);    
    objToBeWritten.push_back(hBCorr.at(1).at(ptDijetBin));

    if( debug ) std::cout << "   Genjet balance vs eta, bins of ptDijet\n";
    // Genjet balance vs eta, bins of ptDijetGen
    fit2DMean(h2BCorrL2L3.at(1).at(ptDijetBin),
		   hBCorrL2L3.at(1).at(ptDijetBin),
		   hBDistCorrL2L3.at(1).at(ptDijetBin),colorGen);
    hBCorrL2L3.at(1).at(ptDijetBin)->SetTitle(title);    
    objToBeWritten.push_back(hBCorrL2L3.at(1).at(ptDijetBin));


    if( debug ) std::cout << "   Uncorrected response vs eta, bins of ptGen\n";
    sprintf(title,"%.0f < p_{T,gen} < %.0f GeV",bins.xLow(bins.bin(ptDijetBin,0)),bins.xUp(bins.bin(ptDijetBin,0)));
    // Uncorrected response vs eta, bins of ptGen
    fit2DGaussMean(h2RUncorr.at(0).at(ptDijetBin),
		   hRUncorr.at(0).at(ptDijetBin),colorUncorr);
    hRUncorr.at(0).at(ptDijetBin)->SetTitle(title);    
    objToBeWritten.push_back(hRUncorr.at(0).at(ptDijetBin));

    if( debug ) std::cout << "   Corrected response vs eta, bins of ptGen\n";
    // Corrected response vs eta, bins of ptGen
    fit2DGaussMean(h2RCorr.at(0).at(ptDijetBin),
		   hRCorr.at(0).at(ptDijetBin),colorCorr);
    hRCorr.at(0).at(ptDijetBin)->SetTitle(title);    
    objToBeWritten.push_back(hRCorr.at(0).at(ptDijetBin));

    if( debug ) std::cout << "   L2L3 corrected response vs eta, bins of ptGen\n";
    // L2L3 corrected response vs eta, bins of ptGen
    fit2DGaussMean(h2RCorrL2L3.at(0).at(ptDijetBin),
		   hRCorrL2L3.at(0).at(ptDijetBin),colorGen);
    hRCorrL2L3.at(0).at(ptDijetBin)->SetTitle(title);    
    objToBeWritten.push_back(hRCorrL2L3.at(0).at(ptDijetBin));
  } // End of loop over pt bins

  if( debug ) std::cout << " Loop over eta bins\n";

  for(int etaBin = 0; etaBin < bins.nBinsY(); etaBin++) { // Loop over eta bins
    char title[50];
    sprintf(title,"%.1f < #eta < %.1f",bins.yLow(bins.bin(0,etaBin)),bins.yUp(bins.bin(0,etaBin)));
    // Uncorrected balance vs ptDijet, bins of eta
    fit2DMean(h2BUncorr.at(2).at(etaBin),
	      hBUncorr.at(2).at(etaBin),
	      hBDistUncorr.at(2).at(etaBin),colorUncorr);
    hBUncorr.at(2).at(etaBin)->SetTitle(title);    
    objToBeWritten.push_back(hBUncorr.at(2).at(etaBin));

    // Corrected balance vs ptDijet, bins of eta
    fit2DMean(h2BCorr.at(2).at(etaBin),
	      hBCorr.at(2).at(etaBin),
	      hBDistCorr.at(2).at(etaBin),colorCorr);
    hBCorr.at(2).at(etaBin)->SetTitle(title);    
    objToBeWritten.push_back(hBCorr.at(2).at(etaBin));

    // Genjet balance vs ptDijet, bins of eta
    fit2DMean(h2BCorrL2L3.at(2).at(etaBin),
	      hBCorrL2L3.at(2).at(etaBin),
	      hBDistCorrL2L3.at(2).at(etaBin),colorGen);
    hBCorrL2L3.at(2).at(etaBin)->SetTitle(title);    
    objToBeWritten.push_back(hBCorrL2L3.at(2).at(etaBin));

    // Uncorrected balance vs ptDijetGen, bins of eta
    fit2DMean(h2BUncorr.at(3).at(etaBin),
	      hBUncorr.at(3).at(etaBin),
	      hBDistUncorr.at(3).at(etaBin),colorUncorr);
    hBUncorr.at(3).at(etaBin)->SetTitle(title);    
    objToBeWritten.push_back(hBUncorr.at(3).at(etaBin));

    // Corrected balance vs ptDijetGen, bins of eta
    fit2DMean(h2BCorr.at(3).at(etaBin),
	      hBCorr.at(3).at(etaBin),
	      hBDistCorr.at(3).at(etaBin),colorCorr);
    hBCorr.at(3).at(etaBin)->SetTitle(title);    
    objToBeWritten.push_back(hBCorr.at(3).at(etaBin));

    // Genjet balance vs ptDijetGen, bins of eta
    fit2DMean(h2BCorrL2L3.at(3).at(etaBin),
	      hBCorrL2L3.at(3).at(etaBin),
	      hBDistCorrL2L3.at(3).at(etaBin),colorGen);
    hBCorrL2L3.at(3).at(etaBin)->SetTitle(title);    
    objToBeWritten.push_back(hBCorrL2L3.at(3).at(etaBin));


    // Uncorrected response vs ptGen, bins of eta
    fit2DGaussMean(h2RUncorr.at(1).at(etaBin),
		   hRUncorr.at(1).at(etaBin),colorUncorr);
    hRUncorr.at(1).at(etaBin)->SetTitle(title);    
    objToBeWritten.push_back(hRUncorr.at(1).at(etaBin));

    // Corrected response vs ptGen, bins of eta
    fit2DGaussMean(h2RCorr.at(1).at(etaBin),
		   hRCorr.at(1).at(etaBin),colorCorr);
    hRCorr.at(1).at(etaBin)->SetTitle(title);    
    objToBeWritten.push_back(hRCorr.at(1).at(etaBin));

    // L2L3 corrected response vs ptGen, bins of eta
    fit2DGaussMean(h2RCorrL2L3.at(1).at(etaBin),
		   hRCorrL2L3.at(1).at(etaBin),colorGen);
    hRCorrL2L3.at(1).at(etaBin)->SetTitle(title);    
    objToBeWritten.push_back(hRCorrL2L3.at(1).at(etaBin));

  } // End of loop over eta bins

  if( debug ) std::cout << "Drawing histograms\n";


  // -- Draw histograms into ps file ----------------------------
  // Draw mean balance
  TPostScript * const ps1 = new TPostScript("controlplotsTwoJetsPtBalanceMean.ps",111);
  TCanvas * const c1 = new TCanvas("c1","",500,500); // Mean values

  c1->Draw();
  ps1->NewPage();
  h2Eta->Draw();
  c1->Draw();
  ps1->NewPage();

  hPtDijetSpectrum->Draw("E");
  c1->SetLogy(1);
  c1->Draw();
  ps1->NewPage();
  c1->SetLogy(0);

  TLegend * legBal = new TLegend(0.4,0.65,0.93,0.85);
  legBal->SetBorderSize(0);
  legBal->SetFillColor(0);
  legBal->SetTextFont(42);
  legBal->AddEntry(hBUncorr.at(0).at(0),"Uncorrected","P");
  legBal->AddEntry(hBCorr.at(0).at(0),"Kalibri","P");
  //legBal->AddEntry(hBCorrL2L3.at(0).at(0),"L2L3 correction","P");

  TH1F * hist = hBUncorr.at(0).at(0);
  TLine *lEtaBal = new TLine(hist->GetXaxis()->GetBinLowEdge(1),0,
			  hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX()),0);
  lEtaBal->SetLineColor(1);
  lEtaBal->SetLineWidth(2);
  lEtaBal->SetLineStyle(2);

  double bMin = -1.;
  double bMax = 2.;
  double bMinZoom = -0.1;
  double bMaxZoom = 0.2;

  for(int ptDijetBin = 0; ptDijetBin < bins.nBinsX(); ptDijetBin++) {
    // Balance vs eta
    c1->cd();
    hBUncorr.at(0).at(ptDijetBin)->GetYaxis()->SetRangeUser(bMin,bMax);
    hBUncorr.at(0).at(ptDijetBin)->GetYaxis()->SetTitle(ptSumLabel);
    hBUncorr.at(0).at(ptDijetBin)->Draw("PE1");
    lEtaBal->Draw("same");
    hBCorr.at(0).at(ptDijetBin)->Draw("PE1same");
    //hBCorrL2L3.at(0).at(ptDijetBin)->Draw("PE1same");
    legBal->Draw("same");
    c1->Draw();
    ps1->NewPage();
  }
  for(int ptDijetBin = 0; ptDijetBin < bins.nBinsX(); ptDijetBin++) {
    // Zoom: Balance vs eta
    c1->cd();
    hBUncorr.at(0).at(ptDijetBin)->GetYaxis()->SetRangeUser(bMinZoom,bMaxZoom);
    hBUncorr.at(0).at(ptDijetBin)->Draw("PE1");
    lEtaBal->Draw("same");
    hBCorr.at(0).at(ptDijetBin)->Draw("PE1same");
    //hBCorrL2L3.at(0).at(ptDijetBin)->Draw("PE1same");
    legBal->Draw("same");
    c1->Draw();
    ps1->NewPage();
  }

  for(int ptDijetBin = 0; ptDijetBin < bins.nBinsX(); ptDijetBin++) {
    // Balance vs eta
    c1->cd();
    hBUncorr.at(1).at(ptDijetBin)->GetYaxis()->SetRangeUser(bMin,bMax);
    hBUncorr.at(1).at(ptDijetBin)->GetYaxis()->SetTitle(ptSumLabel);
    hBUncorr.at(1).at(ptDijetBin)->Draw("PE1");
    lEtaBal->Draw("same");
    hBCorr.at(1).at(ptDijetBin)->Draw("PE1same");
    //hBCorrL2L3.at(1).at(ptDijetBin)->Draw("PE1same");
    legBal->Draw("same");
    c1->Draw();
    ps1->NewPage();
  }
  for(int ptDijetBin = 0; ptDijetBin < bins.nBinsX(); ptDijetBin++) {
    // Zoom: Balance vs eta
    c1->cd();
    hBUncorr.at(1).at(ptDijetBin)->GetYaxis()->SetRangeUser(bMinZoom,bMaxZoom);
    hBUncorr.at(1).at(ptDijetBin)->Draw("PE1");
    lEtaBal->Draw("same");
    hBCorr.at(1).at(ptDijetBin)->Draw("PE1same");
    //hBCorrL2L3.at(1).at(ptDijetBin)->Draw("PE1same");
    legBal->Draw("same");
    c1->Draw();
    ps1->NewPage();
  }

  hist = hBUncorr.at(2).at(0);
  TLine *lPtDijet = new TLine(hist->GetXaxis()->GetBinLowEdge(1),0,
			  hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX()),0);
  lPtDijet->SetLineColor(1);
  lPtDijet->SetLineWidth(2);
  lPtDijet->SetLineStyle(2);

  for(int etaBin = 0; etaBin < bins.nBinsY(); etaBin++) {
    // Balance vs ptDijet
    c1->cd()->SetLogx(1);
    hBUncorr.at(2).at(etaBin)->GetYaxis()->SetRangeUser(bMin,bMax);
    hBUncorr.at(2).at(etaBin)->GetYaxis()->SetTitle(ptSumLabel);
    hBUncorr.at(2).at(etaBin)->Draw("PE1");
    lPtDijet->Draw("same");
    hBCorr.at(2).at(etaBin)->Draw("PE1same");
    //hBCorrL2L3.at(2).at(etaBin)->Draw("PE1same");
    legBal->Draw("same");
    c1->Draw();
    ps1->NewPage();
  }
  for(int etaBin = 0; etaBin < bins.nBinsY(); etaBin++) {
    // Zoom: Balance vs ptDijet
    c1->cd()->SetLogx(1);
    hBUncorr.at(2).at(etaBin)->GetYaxis()->SetRangeUser(bMinZoom,bMaxZoom);
    hBUncorr.at(2).at(etaBin)->Draw("PE1");
    lPtDijet->Draw("same");
    hBCorr.at(2).at(etaBin)->Draw("PE1same");
    //hBCorrL2L3.at(2).at(etaBin)->Draw("PE1same");
    legBal->Draw("same");
    c1->Draw();
    ps1->NewPage();
  }

  for(int etaBin = 0; etaBin < bins.nBinsY(); etaBin++) {
    // Balance vs ptDijetGen
    c1->cd()->SetLogx(1);
    hBUncorr.at(3).at(etaBin)->GetYaxis()->SetRangeUser(bMin,bMax);
    hBUncorr.at(3).at(etaBin)->GetYaxis()->SetTitle(ptSumLabel);
    hBUncorr.at(3).at(etaBin)->Draw("PE1");
    lPtDijet->Draw("same");
    hBCorr.at(3).at(etaBin)->Draw("PE1same");
    //hBCorrL2L3.at(3).at(etaBin)->Draw("PE1same");
    legBal->Draw("same");
    c1->Draw();
    ps1->NewPage();
  }
  for(int etaBin = 0; etaBin < bins.nBinsY(); etaBin++) {
    // Zoom: Balance vs ptDijetGen
    c1->cd()->SetLogx(1);
    hBUncorr.at(3).at(etaBin)->GetYaxis()->SetRangeUser(bMinZoom,bMaxZoom);
    hBUncorr.at(3).at(etaBin)->Draw("PE1");
    lPtDijet->Draw("same");
    hBCorr.at(3).at(etaBin)->Draw("PE1same");
    //hBCorrL2L3.at(3).at(etaBin)->Draw("PE1same");
    legBal->Draw("same");
    c1->Draw();
    ps1->NewPage();
  }


  // Draw response
  double rMin = 0.;
  double rMax = 2.5;
  double rMinZoom = 0.8;
  double rMaxZoom = 1.4;

  TLegend * legResp = new TLegend(0.4,0.65,0.93,0.85);
  legResp->SetBorderSize(0);
  legResp->SetFillColor(0);
  legResp->SetTextFont(42);
  legResp->AddEntry(hRUncorr.at(0).at(0),"Uncorrected","P");
  legResp->AddEntry(hRCorr.at(0).at(0),"Kalibri","P");
  //legResp->AddEntry(hRCorrL2L3.at(0).at(0),"L2L3 correction","P");

  hist = hRUncorr.at(0).at(0);
  TLine *lEtaResp = new TLine(hist->GetXaxis()->GetBinLowEdge(1),1,
			  hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX()),1);
  lEtaResp->SetLineColor(1);
  lEtaResp->SetLineWidth(2);
  lEtaResp->SetLineStyle(2);

  hist = hRUncorr.at(1).at(0);
  TLine *lPtGenResp = new TLine(hist->GetXaxis()->GetBinLowEdge(1),1,
			  hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX()),1);
  lPtGenResp->SetLineColor(1);
  lPtGenResp->SetLineWidth(2);
  lPtGenResp->SetLineStyle(2);

  for(int ptGenBin = 0; ptGenBin < bins.nBinsX(); ptGenBin++) {
    // Response vs eta
    c1->cd()->SetLogx(0);
    hRUncorr.at(0).at(ptGenBin)->GetYaxis()->SetRangeUser(rMin,rMax);
    hRUncorr.at(0).at(ptGenBin)->GetYaxis()->SetTitle(respLabel);
    hRUncorr.at(0).at(ptGenBin)->Draw("PE1");
    lEtaResp->Draw("same");
    hRCorr.at(0).at(ptGenBin)->Draw("PE1same");
    //hRCorrL2L3.at(0).at(ptGenBin)->Draw("PE1same");
    legResp->Draw("same");
    c1->Draw();
    ps1->NewPage();
  }
  for(int ptGenBin = 0; ptGenBin < bins.nBinsX(); ptGenBin++) {
    // Zoom: Response vs eta
    c1->cd();
    hRUncorr.at(0).at(ptGenBin)->GetYaxis()->SetRangeUser(rMinZoom,rMaxZoom);
    hRUncorr.at(0).at(ptGenBin)->Draw("PE1");
    hRCorr.at(0).at(ptGenBin)->Draw("PE1same");
    //hRCorrL2L3.at(0).at(ptGenBin)->Draw("PE1same");
    legResp->Draw("same");
    c1->Draw();
    ps1->NewPage();
  }
  for(int etaBin = 0; etaBin < bins.nBinsY(); etaBin++) {
    // Response vs ptGen
    c1->cd()->SetLogx(1);
    hRUncorr.at(1).at(etaBin)->GetYaxis()->SetRangeUser(rMin,rMax);
    hRUncorr.at(1).at(etaBin)->GetYaxis()->SetTitle(respLabel);
    hRUncorr.at(1).at(etaBin)->Draw("PE1");
    lPtGenResp->Draw("same");
    hRCorr.at(1).at(etaBin)->Draw("PE1same");
    //hRCorrL2L3.at(1).at(etaBin)->Draw("PE1same");
    legResp->Draw("same");
    c1->Draw();
    ps1->NewPage();
  }
  for(int etaBin = 0; etaBin < bins.nBinsY(); etaBin++) {
    // Zoom: Response vs ptGen
    c1->cd();
    hRUncorr.at(1).at(etaBin)->GetYaxis()->SetRangeUser(rMinZoom,rMaxZoom);
    hRUncorr.at(1).at(etaBin)->Draw("PE1");
    hRCorr.at(1).at(etaBin)->Draw("PE1same");
    //hRCorrL2L3.at(1).at(etaBin)->Draw("PE1same");
    legResp->Draw("same");
    c1->Draw();
    ps1->NewPage();
  }

  ps1->Close();


  // Draw distributions
  for(int i = 0; i < 4; i++) {
    int nBins = bins.nBinsX();
    if( i > 1 ) nBins = bins.nBinsY();
    for(int bin = 0; bin < nBins; bin++) {
      for(size_t k = 0; k < hBDistUncorr.at(i).at(bin).size(); k++) {
	double min = 0.;
	double maxUncorr = 0.;
	double maxCorr = 0.;
	double maxCorrL2L3 = 0.;
	findYRange(hBDistUncorr.at(i).at(bin).at(k),min,maxUncorr);
	findYRange(hBDistCorr.at(i).at(bin).at(k),min,maxCorr);
	findYRange(hBDistCorrL2L3.at(i).at(bin).at(k),min,maxCorrL2L3);

	//	if( maxCorrL2L3 > maxCorr ) maxCorr = maxCorrL2L3;
	if( maxUncorr > maxCorr ) maxCorr = maxUncorr;
	maxCorr *= 1.3;
	hBDistUncorr.at(i).at(bin).at(k)->GetYaxis()->SetRangeUser(0.,maxCorr);
	hBDistCorr.at(i).at(bin).at(k)->GetYaxis()->SetRangeUser(0.,maxCorr);
	hBDistCorrL2L3.at(i).at(bin).at(k)->GetYaxis()->SetRangeUser(0.,maxCorr);
      }
    }
  }

  TPostScript * const ps2 = new TPostScript("controlplotsTwoJetsPtBalanceDistributions.ps",112);
  ps2->Range(25,1);

  double cw  = 300;
  double mw  = 70;

  double w   = 3*cw + mw;
  double h   = 2*cw + 2*mw;
  
  double crw = cw/w;
  double mrw = mw/w;
  double crh = cw/h;
  double mrh = mw/h;

  TCanvas * const c2 = new TCanvas("c2","",(int)w,(int)h);

  // Pads for the histograms
  std::vector<TPad*> cPads;
  cPads.push_back(new TPad("cPad0","",mrw,mrh+crh,mrw+crw,1.-mrw));
  cPads.push_back(new TPad("cPad1","",mrw+crw,mrh+crh,mrw+2*crw,1.-mrw));
  cPads.push_back(new TPad("cPad2","",mrw+2*crw,mrh+crh,1.,1.-mrw));
  cPads.push_back(new TPad("cPad3","",mrw,mrh,mrw+crw,mrh+crh));
  cPads.push_back(new TPad("cPad4","",mrw+crw,mrh,mrw+2*crw,mrh+crh));
  cPads.push_back(new TPad("cPad5","",mrw+2*crw,mrh,1.,mrh+crh));
  for(unsigned int i = 0; i < cPads.size(); i++) {
    cPads.at(i)->SetFillStyle(1001); 
    cPads.at(i)->SetFrameFillColor(10 );
    cPads.at(i)->SetFrameBorderMode(0);
    cPads.at(i)->SetTopMargin(0.04);
    cPads.at(i)->SetBottomMargin(0.1);
    cPads.at(i)->SetLeftMargin(0.1);
    cPads.at(i)->SetRightMargin(0.04);
    
    c2->cd();
    cPads.at(i)->Draw();
  }

  // Pads for margins holding titles
  TPad * mbPad = new TPad("mbPad","",0.,0.,1.,mrh);
  mbPad->SetFillStyle(1001);
  mbPad->SetFrameFillColor(10);
  mbPad->SetFrameBorderMode(0);
 
  TPaveText * bLabel = new TPaveText(0.6,0.5,0.95,1.,"NDC");
  bLabel->SetFillColor(0);
  bLabel->SetTextFont(42);
  bLabel->SetTextSize(0.7);
  bLabel->SetBorderSize(0);
  bLabel->AddText(ptSumLabel);
  c2->cd();
  mbPad->Draw();
  mbPad->cd();
  bLabel->Draw();

  TPad * mtPad = new TPad("mtPad","",0.,mrh+2.*crh+0.01,1.,1.);
  mtPad->SetFillStyle(1001);
  mtPad->SetFrameFillColor(10);
  mtPad->SetFrameBorderMode(0);
  c2->cd();
  mtPad->Draw();
 
  std::vector<TPaveText*> tLabels(2*(bins.nBinsX()+bins.nBinsY()));
  for(int i = 0; i < 2*(bins.nBinsX()+bins.nBinsY()); i++) {
    tLabels.at(i) = new TPaveText(0.3,0.6,0.7,1.,"NDC");
    tLabels.at(i)->SetFillColor(0);
    tLabels.at(i)->SetTextFont(42);
    tLabels.at(i)->SetTextSize(0.5);
    tLabels.at(i)->SetBorderSize(0);
  }

  // For some reason, this prevents the first ps-page
  // from looking weird...
  c2->Draw();
  ps2->NewPage();
  // Loop over ptDijet bins
  for(int ptDijetBin = 0; ptDijetBin < bins.nBinsX(); ptDijetBin++) {
    char label[100];
    sprintf(label,"%.0f < p^{dijet}_{T} < %.0f GeV",
	    bins.xLow(bins.bin(ptDijetBin,0)),bins.xUp(bins.bin(ptDijetBin,0)));
    tLabels.at(ptDijetBin)->AddText(label);
    mtPad->Clear();
    mtPad->cd();
    tLabels.at(ptDijetBin)->Draw();

    for(size_t k = 0; k < hBDistUncorr.at(0).at(ptDijetBin).size(); k++) {
      int padIdx = k % 6;
      cPads.at(padIdx)->cd();

      TString etaBinText = hBDistUncorr.at(0).at(ptDijetBin).at(k)->GetTitle();
      hBDistUncorr.at(0).at(ptDijetBin).at(k)->SetTitle("");

      hBDistUncorr.at(0).at(ptDijetBin).at(k)->Draw("H");
      hBDistCorr.at(0).at(ptDijetBin).at(k)->Draw("Hsame");
      //      hBDistCorrL2L3.at(0).at(ptDijetBin).at(k)->Draw("Hsame");

      // Label bin
      TPaveText * etaBinLabel = new TPaveText(0.13,0.84,0.93,0.93,"NDC");
      etaBinLabel->SetFillColor(0);
      etaBinLabel->SetTextFont(42);
      etaBinLabel->SetTextAlign(12);
      etaBinLabel->AddText(etaBinText);
      etaBinLabel->Draw("same");

      if( padIdx == 5 || k == hBDistUncorr.at(0).at(ptDijetBin).size()-1 ) {
	c2->Draw();
	ps2->NewPage();
	for(std::vector<TPad*>::iterator p = cPads.begin(); p != cPads.end(); p++) {
	  (*p)->Clear();
	}
      }
    }
  }
  // Loop over ptDijetGen bins
  for(int ptDijetBin = 0; ptDijetBin < bins.nBinsX(); ptDijetBin++) {
    char label[100];
    sprintf(label,"%.0f < p^{dijet}_{T,gen} < %.0f GeV",
	    bins.xLow(bins.bin(ptDijetBin,0)),bins.xUp(bins.bin(ptDijetBin,0)));
    tLabels.at(bins.nBinsX()+ptDijetBin)->AddText(label);
    mtPad->Clear();
    mtPad->cd();
    tLabels.at(bins.nBinsX()+ptDijetBin)->Draw();

    for(size_t k = 0; k < hBDistUncorr.at(1).at(ptDijetBin).size(); k++) {
      int padIdx = k % 6;
      cPads.at(padIdx)->cd();

      TString etaBinText = hBDistUncorr.at(1).at(ptDijetBin).at(k)->GetTitle();
      hBDistUncorr.at(1).at(ptDijetBin).at(k)->SetTitle("");

      hBDistUncorr.at(1).at(ptDijetBin).at(k)->Draw("H");
      hBDistCorr.at(1).at(ptDijetBin).at(k)->Draw("Hsame");
      //      hBDistCorrL2L3.at(1).at(ptDijetBin).at(k)->Draw("Hsame");

      // Label bin
      TPaveText * etaBinLabel = new TPaveText(0.13,0.84,0.93,0.93,"NDC");
      etaBinLabel->SetFillColor(0);
      etaBinLabel->SetTextFont(42);
      etaBinLabel->SetTextAlign(12);
      etaBinLabel->AddText(etaBinText);
      etaBinLabel->Draw("same");

      if( padIdx == 5 || k == hBDistUncorr.at(1).at(ptDijetBin).size()-1 ) {
	c2->Draw();
	ps2->NewPage();
	for(std::vector<TPad*>::iterator p = cPads.begin(); p != cPads.end(); p++) {
	  (*p)->Clear();
	}
      }
    }
  }
  // Loop over eta bins
  for(int etaBin = 0; etaBin < bins.nBinsY(); etaBin++) {
    char label[100];
    sprintf(label,"%.2f < #eta < %.2f",
	    bins.yLow(bins.bin(0,etaBin)),bins.yUp(bins.bin(0,etaBin)));
    tLabels.at(2*bins.nBinsX()+etaBin)->AddText(label);
    mtPad->Clear();
    mtPad->cd();
    tLabels.at(2*bins.nBinsX()+etaBin)->Draw();

    for(size_t k = 0; k < hBDistUncorr.at(2).at(etaBin).size(); k++) {
      int padIdx = k % 6;
      cPads.at(padIdx)->cd();

      TString ptDijetBinText = hBDistUncorr.at(2).at(etaBin).at(k)->GetTitle();
      hBDistUncorr.at(2).at(etaBin).at(k)->SetTitle("");

      hBDistUncorr.at(2).at(etaBin).at(k)->Draw("H");
      hBDistCorr.at(2).at(etaBin).at(k)->Draw("Hsame");
      //      hBDistCorrL2L3.at(2).at(etaBin).at(k)->Draw("Hsame");

      // Label bin
      TPaveText * ptDijetBinLabel = new TPaveText(0.13,0.84,0.93,0.93,"NDC");
      ptDijetBinLabel->SetFillColor(0);
      ptDijetBinLabel->SetTextFont(42);
      ptDijetBinLabel->SetTextAlign(12);
      ptDijetBinLabel->AddText(ptDijetBinText);
      ptDijetBinLabel->Draw("same");

      if( padIdx == 5 || k == hBDistUncorr.at(2).at(etaBin).size()-1 ) {
	c2->Draw();
	ps2->NewPage();
	for(std::vector<TPad*>::iterator p = cPads.begin(); p != cPads.end(); p++) {
	  (*p)->Clear();
	}
      }
    }
  }
  // Loop over eta bins
  for(int etaBin = 0; etaBin < bins.nBinsY(); etaBin++) {
    char label[100];
    sprintf(label,"%.2f < #eta < %.2f",
	    bins.yLow(bins.bin(0,etaBin)),bins.yUp(bins.bin(0,etaBin)));
    tLabels.at(2*bins.nBinsX()+bins.nBinsY()+etaBin)->AddText(label);
    mtPad->Clear();
    mtPad->cd();
    tLabels.at(2*bins.nBinsX()+bins.nBinsY()+etaBin)->Draw();

    for(size_t k = 0; k < hBDistUncorr.at(3).at(etaBin).size(); k++) {
      int padIdx = k % 6;
      cPads.at(padIdx)->cd();

      TString ptDijetGenBinText = hBDistUncorr.at(3).at(etaBin).at(k)->GetTitle();
      hBDistUncorr.at(3).at(etaBin).at(k)->SetTitle("");

      hBDistUncorr.at(3).at(etaBin).at(k)->Draw("H");
      hBDistCorr.at(3).at(etaBin).at(k)->Draw("Hsame");
      //      hBDistCorrL2L3.at(3).at(etaBin).at(k)->Draw("Hsame");

      // Label bin
      TPaveText * ptDijetGenBinLabel = new TPaveText(0.13,0.84,0.93,0.93,"NDC");
      ptDijetGenBinLabel->SetFillColor(0);
      ptDijetGenBinLabel->SetTextFont(42);
      ptDijetGenBinLabel->SetTextAlign(12);
      ptDijetGenBinLabel->AddText(ptDijetGenBinText);
      ptDijetGenBinLabel->Draw("same");

      if( padIdx == 5 || k == hBDistUncorr.at(3).at(etaBin).size()-1 ) {
	c2->Draw();
	ps2->NewPage();
	for(std::vector<TPad*>::iterator p = cPads.begin(); p != cPads.end(); p++) {
	  (*p)->Clear();
	}
      }
    }
  }
  ps2->Close();

  if( outputROOT_ ) writeToRootFile(objToBeWritten,"TwoJetsPtBalanceEvent");

  if( debug ) std::cout << "Cleaning up\n";

  // Clean up
  for(size_t i = 0; i < objToBeWritten.size(); i++) {
    delete objToBeWritten.at(i);
  }
  for(int i = 0; i < 4; i++) {
    int nBins = bins.nBinsX();
    if( i > 1 ) nBins = bins.nBinsY();
    for(int bin = 0; bin < nBins; bin++) {
      delete h2BUncorr.at(i).at(bin);
      delete h2BCorr.at(i).at(bin);
      delete h2BCorrL2L3.at(i).at(bin);
      for(size_t k = 0; k < hBDistUncorr.size(); k++) {
	delete hBDistUncorr.at(i).at(bin).at(k);
      }
      for(size_t k = 0; k < hBDistCorr.size(); k++) {
	delete hBDistCorr.at(i).at(bin).at(k);
      }
      for(size_t k = 0; k < hBDistCorrL2L3.size(); k++) {
	delete hBDistCorrL2L3.at(i).at(bin).at(k);
      }
    }
  }
  for(int i = 0; i < 2; i++) {
    int nBins = bins.nBinsX();
    if( i > 0 ) nBins = bins.nBinsY();
    for(int bin = 0; bin < nBins; bin++) {
      delete h2RUncorr.at(i).at(bin);
      delete h2RCorr.at(i).at(bin);
      delete h2RCorrL2L3.at(i).at(bin);
    }
  }

  delete c1;
  delete c2;
  delete ps1;
  delete ps2;
  delete legBal;
  delete legResp;
  delete lPtDijet;
  delete lEtaBal;
  delete lEtaResp;
  delete lPtGenResp;

  if( debug ) std::cout << "Leaving makeControlPlotsTwoJetsPtBalance()\n";
}



//!  Filling \p bins with borders of \p nBins bins between \p first
//!  and \p last that are equidistant when viewed in log scale,
//!  so \p bins must have length \p nBins+1. If \p first, \p last
//!  or \p nBins are not positive, failure is reported.
// -------------------------------------------------------------
bool TControlPlots::equidistLogBins(double * bins, int nBins, double first, double last) const {
  if( nBins < 1 || first <= 0. || last <= 0. ) return false;
  bins[0]     = first;
  bins[nBins] = last;
  const double firstLog = log10(bins[0]);
  const double lastLog  = log10(bins[nBins]);
  for (int i = 1; i < nBins; ++i) {
    bins[i] = pow(10., firstLog + i*(lastLog-firstLog)/(nBins));
  }

  return true;
}


//!  \brief Find y-axis range
//!
//!  Sets \p min and \p max to the minimum (non-zero) and
//!  maximum bin content of \p h, respectively.
// --------------------------------------------------
void TControlPlots::findYRange(const TH1F * h, double& min, double& max) const {
  min = 10000.;
  max = 0.;
  for(int bin = 1; bin <= h->GetNbinsX(); bin++) {
    double val = h->GetBinContent(bin);
    if( val > 0. && val < min ) min = val;
    if( val > 0. && val > max ) max = val;
  }
  if( min > max ) {
    min = 1E-3;
    max = 1;
  }
}



//!  \brief Interface to fit2D(const TH2F* hist, std::vector<TH1F*>& hresults,
//!  std::vector<TH1F*>& distributions, std::vector<TF1*>& gaussFits, const bool plotgauss)
//!  using the old signature of \p fit2D.
//---------------------------------------------------------------
void TControlPlots::fit2D(const TH2F* hist, TH1F* hresults[8], TH1F* gaussplots[4], TF1* gf[4], const bool plotgauss) const {
  std::vector<TH1F*> vHresults;
  std::vector<TH1F*> vDistributions;
  std::vector<TF1*> vGaussFits;
  fit2D(hist,vHresults,vDistributions,vGaussFits,plotgauss);
  for(int i = 0; i < 8; i++) {
    hresults[i] = vHresults.at(i);
  }
  for(int i = 0; i < 4; i++) {
    int idx = i * vGaussFits.size() / 4;
    gaussplots[i] = vDistributions.at(idx);
    gf[i] = vGaussFits.at(idx);
  }
}



//!  \brief Get mean of x-profiles of a 2D histogram
//!
//!  Calculates the 1D projections of \p hist along the
//!  x-axis for different x bins. The means versus x are filled
//!  in the 1D histogram \p hresult.
//!
//!  \p distributions contains the actual 1D distributions
//!  per x-bin of \p hist.
//!
//!  The marker color of \p hresult is set to \p color.
//!
//!  \sa fit2D(const TH2F* hist, std::vector<TH1F*>& hresults,
//!  std::vector<TH1F*>& distributions, std::vector<TF1*>& gaussFits, const bool plotgauss)
//---------------------------------------------------------------
void TControlPlots::fit2DMean(const TH2F* hist, TH1F*& hresult,
			      std::vector<TH1F*>& distributions,
			      int color) const {
  std::vector<TH1F*> hResults;
  std::vector<TF1*> gaussFits;
  fit2D(hist,hResults,distributions,gaussFits,true);
  for(size_t i = 0; i < hResults.size(); i++) {
    if( i == 0 ) {
      hresult = hResults.at(i);
      hresult->SetMarkerStyle(20);
      hresult->SetMarkerColor(color);
      hresult->SetLineColor(color);
    } else {
      delete hResults.at(i);
    }
  }
  for(size_t k = 0; k < distributions.size(); k++) {
    distributions.at(k)->SetLineColor(color);
    delete gaussFits.at(k);
  }
}



//!  \brief Get mean of x-profiles of a 2D histogram
//!
//!  Calculates the 1D projections of \p hist along the
//!  x-axis for different x bins. The means versus x are filled
//!  in the 1D histogram \p hresult.
//!
//!  The marker color of \p hresult is set to \p color.
//!
//!  \sa fit2D(const TH2F* hist, std::vector<TH1F*>& hresults,
//!  std::vector<TH1F*>& distributions, std::vector<TF1*>& gaussFits, const bool plotgauss)
//---------------------------------------------------------------
void TControlPlots::fit2DMean(const TH2F* hist, TH1F*& hresult, int color) const {
  std::vector<TH1F*> hResults;
  fit2D(hist,hResults);
  for(size_t i = 0; i < hResults.size(); i++) {
    if( i == 0 ) {
      hresult = hResults.at(i);
      hresult->SetMarkerStyle(20);
      hresult->SetMarkerColor(color);
      hresult->SetLineColor(color);
    } else {
      delete hResults.at(i);
    }
  }
}



//!  \brief Get mean of x-profiles of a 2D histogram from Gauss fit
//!
//!  Calculates the 1D projections of \p hist along the
//!  x-axis for different x bins. The means of central
//!  (\f$ \pm3\sigma \f$) Gauss fits versus x are filled
//!  in the 1D histogram \p hresult.
//!
//!  \p distributions contains the actual 1D distributions
//!  per x-bin of \p hist. \p gaussFits contains the Gaussian
//!  fits on these distributions.
//!
//!  The marker color of \p hresult is set to \p color.
//!
//!  \sa fit2D(const TH2F* hist, std::vector<TH1F*>& hresults,
//!  std::vector<TH1F*>& distributions, std::vector<TF1*>& gaussFits, const bool plotgauss)
//---------------------------------------------------------------
void TControlPlots::fit2DGaussMean(const TH2F* hist, TH1F*& hresult,
				   std::vector<TH1F*>& distributions,
				   std::vector<TF1*>& gaussFits, int color) const {
  std::vector<TH1F*> hResults;
  fit2D(hist,hResults,distributions,gaussFits,true);
  for(size_t i = 0; i < hResults.size(); i++) {
    if( i == 2 ) {
      hresult = hResults.at(i);
      hresult->SetMarkerStyle(20);
      hresult->SetMarkerColor(color);
      hresult->SetLineColor(color);
    } else {
      delete hResults.at(i);
    }
  }
  for(size_t k = 0; k < distributions.size(); k++) {
    distributions.at(k)->SetLineColor(color);
    gaussFits.at(k)->SetLineColor(color);
    gaussFits.at(k)->SetLineWidth(1);
  }
}



//!  \brief Get mean of x-profiles of a 2D histogram from Gauss fit
//!
//!  Calculates the 1D projections of \p hist along the
//!  x-axis for different x bins. The means of central
//!  (\f$ \pm3\sigma \f$) Gauss fits versus x are filled
//!  in the 1D histogram \p hresult.
//!
//!  The marker color of \p hresult is set to \p color.
//!
//!  \sa fit2D(const TH2F* hist, std::vector<TH1F*>& hresults,
//!  std::vector<TH1F*>& distributions, std::vector<TF1*>& gaussFits, const bool plotgauss)
//---------------------------------------------------------------
void TControlPlots::fit2DGaussMean(const TH2F* hist, TH1F*& hresult, int color) const {
  std::vector<TH1F*> hResults;
  fit2D(hist,hResults);
  for(size_t i = 0; i < hResults.size(); i++) {
    if( i == 2 ) {
      hresult = hResults.at(i);
      hresult->SetMarkerStyle(20);
      hresult->SetMarkerColor(color);
      hresult->SetLineColor(color);
    } else {
      delete hResults.at(i);
    }
  }
}



//!  \brief Get different x-profiles of a 2D histogram
//!
//!  Calculates the 1D projections of \p hist along the
//!  x-axis for different x bins. Different quantities
//!  of these projections are calculated and shown
//!  versus x in the 1D histograms \p hresults
//!   - 0: Mean
//!   - 1: RMS
//!   - 2: Mean of central (\f$ \pm3\sigma \f$) Gauss fit
//!   - 3: Width of central Gauss fit
//!   - 4: Median
//!   - 5: \f$ \chi^{2} / ndof\f$
//!   - 6: Fit probability
//!   - 7: Quantiles
//!
//!  The binning and the x-axis title of the histograms
//!  in \p hresults is the same as the x binning of \p hist, the
//!  y-axis title and the histogram title are adjusted
//!  to the displayed quantity.
//!
//!  \p distributions contains the actual 1D distributions
//!  per x-bin of \p hist. \p gaussFits contains the Gaussian
//!  fits on these distributions. The latter two vectors are
//!  only filled, if \p plotgauss is true.
//!
//!  \p hresults are newly created (take care of deleting them!);
//!  their object names are set to:
//!     "(hist-name)_result(X)",
//!  where (hist-name) = hist->GetName() and (X) is the index of
//!  the above specified property (i.e. 0 for "mean").
//---------------------------------------------------------------
void TControlPlots::fit2D(const TH2F* hist,
			  std::vector<TH1F*>& hresults,
			  std::vector<TH1F*>& distributions,
			  std::vector<TF1*>& gaussFits,
			  const bool plotgauss) const
{
  TString quantityName[8];
  quantityName[0] = "mean"; 
  quantityName[1] = "standard deviation"; 
  quantityName[2] = "mean of Gauss fit"; 
  quantityName[3] = "width of Gauss fit"; 
  quantityName[4] = "median"; 
  quantityName[5] = "#chi^{2} / n.d.f."; 
  quantityName[6] = "probability"; 
  quantityName[7] = "quantiles"; 

  //book hists
  hresults = std::vector<TH1F*>(8);
  TString s = hist->GetName();
  s += "_result0";
  if( hist->GetXaxis()->GetXbins()->GetSize() == hist->GetNbinsX() +1)
    {
      hresults[0] = new TH1F(s,"",hist->GetNbinsX(),hist->GetXaxis()->GetXbins()->GetArray());
    }
  else
    {
      hresults[0] = new TH1F(s,"",hist->GetNbinsX(),hist->GetXaxis()->GetXmin(),
			   hist->GetXaxis()->GetXmax());
    }
  hresults[0]->SetXTitle(hist->GetXaxis()->GetTitle());
  hresults[0]->SetMarkerStyle(hist->GetMarkerStyle());
  hresults[0]->SetMarkerColor(hist->GetMarkerColor());
  hresults[0]->SetLineColor(hist->GetLineColor());
  hresults[0]->SetMarkerSize(hist->GetMarkerSize());
  for(int i = 1; i < 8 ; ++i)
    {
      s = hist->GetName();
      s += "_result";
      s += i;
      hresults[i] = (TH1F*) hresults[0]->Clone(s);
      s = hist->GetTitle();
      hresults[i]->SetTitle(s + ",  " + quantityName[i]); 
    }
  s = hist->GetTitle();
  hresults[0]->SetTitle(s + ",  " + quantityName[0]); 

  hresults[5]->SetMinimum(0.0);
  hresults[5]->SetMaximum(100);
  hresults[6]->SetMinimum(0.0);
  hresults[6]->SetMaximum(1.05);

  TH1F* htemp = new TH1F("htemp","",hist->GetNbinsY(),hist->GetYaxis()->GetXmin(),
			   hist->GetYaxis()->GetXmax());
  htemp->Sumw2();

  if(plotgauss) {
    distributions = std::vector<TH1F*>(hist->GetNbinsX());
    gaussFits = std::vector<TF1*>(hist->GetNbinsX());
  }

  const int nq = 2;
  double yq[2],xq[2];
  xq[0] = 0.5;
  xq[1] = 0.90;
  for(int i = 1 ; i <= hist->GetNbinsX() ; ++i)
    {
      htemp->Reset();
      for(int j = 1 ; j <= hist->GetNbinsY() ; ++j)
	{
	  htemp->SetBinContent(j,hist->GetBinContent(hist->GetBin(i,j)));
	  htemp->SetBinError(j,hist->GetBinError(i,j));
	}  
      if(plotgauss) {
	int idx = i - 1;
	char title[100];
	sprintf(title,
		"%.2f < %s < %.2f",
		hist->GetXaxis()->GetBinLowEdge(i),
		hist->GetXaxis()->GetTitle(),
		hist->GetXaxis()->GetBinUpEdge(i));
	s = hist->GetName();
	s += "_distribution";
	s += idx;
	distributions.at(idx) = static_cast<TH1F*>(htemp->Clone(s));
	distributions.at(idx)->SetTitle(title);
      }
      double mean = htemp->GetMean(); 
      double meanerror = htemp->GetMeanError();
      double width = htemp->GetRMS();
      if(width < 0.1) width = 0.1;
      if(htemp->GetSumOfWeights() <= 0) {
	if(plotgauss) {
	  int idx = i - 1;
	  s = hist->GetName();
	  s += "_gaussFit";
	  s += idx;
	  gaussFits.at(idx) = new TF1(s,"gaus");//static_cast<TF1*>(f->Clone(s));
	}
	continue; 
      } else {
	htemp->Fit("gaus","QNO","", mean - 3 * width,mean + 3 * width);
	TF1 *f = (TF1*)gROOT->GetFunction("gaus")->Clone();
	mean = f->GetParameter(1);
	meanerror = f->GetParError(1);
	width = f->GetParameter(2);
	if(width < 0.05) width = 0.05;
	if( (htemp->Fit(f,"LLQNO","goff",mean - 1.5 * width, mean + 1.5 * width) == 0) ) {
	  mean = f->GetParameter(1);
	  meanerror = f->GetParError(1);
	  width = f->GetParameter(2);
	  
	  hresults[2]->SetBinContent(i,mean);
	  hresults[2]->SetBinError(i,meanerror);
	  hresults[3]->SetBinContent(i,width/mean);
	  hresults[3]->SetBinError(i, f->GetParError(2)/mean);
	}
	hresults[5]->SetBinContent(i, f->GetChisquare() / f->GetNumberFreeParameters());
	hresults[5]->SetBinError(i, 0.01);
	hresults[6]->SetBinContent(i, f->GetProb());
	hresults[6]->SetBinError(i, 0.01);
	if(plotgauss) {
	  int idx = i - 1;
	  s = hist->GetName();
	  s += "_gaussFit";
	  s += idx;
	  gaussFits.at(idx) = static_cast<TF1*>(f->Clone(s));
	}
	
	mean = htemp->GetMean();
	meanerror = htemp->GetMeanError();
	width = htemp->GetRMS();
	hresults[0]->SetBinContent(i,mean);
	hresults[0]->SetBinError(i,meanerror);
	hresults[1]->SetBinContent(i,width/mean); 
	hresults[1]->SetBinError(i,htemp->GetRMSError()/mean);
	htemp->GetQuantiles(nq,yq,xq);
	hresults[4]->SetBinContent(i,yq[0]);
	hresults[4]->SetBinError(i,0.0001);
	hresults[7]->SetBinContent(i,yq[1]/yq[0]-1);
	hresults[7]->SetBinError(i,0.0001);
	delete f;
      }
    }
  delete htemp;
}



//!  \brief Read JetMET parameters from file
//!
//!  Reads the constants for L2 and L3 correction from
//!  an ascii file in \p CondDB format and stores them
//!  in \p TParameters::GetPars(). Afterwards 
//!  \p TData::GetParametrizedMess() will return the
//!  correction pt using these constants.
//!  The file names have to be specified in the
//!  'Control plots JetMET L2L3 constants' string in
//!  the config file. If no file name is given, the
//!  parameters are left unchanged.
//!
//!  \note This works only with the
//!        \p L2L3JetParametrization
//!
//!  \note The fitted parameter values can be reset
//!        using \p resetFittedParameters()
//!
//!  \return True, if parameter values where changed,
//!          false otherwise.
//---------------------------------------------------------------
bool TControlPlots::readJetMETParameters() {
  bool setReadConstants = false;

  std::vector<std::string> inputFileNames
    = bag_of_string( config_->read<string>("Control plots JetMET L2L3 constants",";") );

  if( inputFileNames.size() ) {
    setReadConstants = true;

    std::cout << "  Using constants for comparison to JetMET L2L3 corrections from files:\n";

    // Read JetMET parameter values
    std::string corrL2FileName;
    std::string corrL3FileName;
    for(size_t i = 0; i < inputFileNames.size(); i++) {
      if( inputFileNames.at(i).find("L2") != std::string::npos ) 
	corrL2FileName = inputFileNames.at(i);
      else if( inputFileNames.at(i).find("L3") != std::string::npos ) 
	corrL3FileName = inputFileNames.at(i);
    }
    
    if( !corrL2FileName.empty() || !corrL3FileName.empty() ) {
      if( !corrL2FileName.empty() ) {
	std::cout << "    L2: " << corrL2FileName << std::endl;
	par_->readCalibrationJetMETL2(corrL2FileName);
      }
      if( !corrL3FileName.empty() ) {
	std::cout << "    L3: " << corrL3FileName << std::endl;
	par_->readCalibrationJetMETL3(corrL3FileName);
      }
    } else {
      std::cerr << "    ERROR: Could not find specified file(s).\n";
      std::cerr << "           Using stored JetMET L2L3 corrections.\n";
    }
  } else {
    std::cout << "  Using stored JetMET L2L3 corrections for comparison.\n";
  }

  return setReadConstants;
}



//!  \brief Reset the parameter values in \p TParameters 
//!         to the fitted values
//---------------------------------------------------------------
void TControlPlots::resetFittedParameters() {
  for(int i = 0; i < par_->GetNumberOfParameters(); i++) {
    par_->GetPars()[i] = fittedPar_.at(i);
  }
}



//!  Set style option for ps output.
//---------------------------------------------------------------
void TControlPlots::setGStyle() const
{
  // For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(800); //Height of canvas
  gStyle->SetCanvasDefW(800); //Width of canvas
  gStyle->SetCanvasDefX(0);   //Position on screen
  gStyle->SetCanvasDefY(0);

  // For the frame
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(kBlack);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetFrameLineStyle(0);
  gStyle->SetFrameLineWidth(1);

  // For the Pad:
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  // For the histo:
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);

  // For the statistics box:
  gStyle->SetOptStat(0);
  //  gStyle->SetOptStat("neMR");
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatX(0.92);              
  gStyle->SetStatY(0.86);              
  gStyle->SetStatH(0.16);
  gStyle->SetStatW(0.22);

  // For the leegnd
  gStyle->SetLegendBorderSize(1);

  // Margins:
  gStyle->SetPadTopMargin(0.11);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadLeftMargin(0.25);
  gStyle->SetPadRightMargin(0.04);

  // For the Global title:
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFont(42,"");
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleX(0.6);
  gStyle->SetTitleH(0.05);
  gStyle->SetTitleBorderSize(0);

  // For the axis titles:
  gStyle->SetTitleColor(1,"XYZ");
  gStyle->SetLabelColor(1,"XYZ");
  // For the axis labels:
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetLabelOffset(0.007,"XYZ");
  gStyle->SetLabelSize(0.045,"XYZ");
  // For the axis titles:
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(2.0);

  // For the axis:
  gStyle->SetAxisColor(1,"XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03,"XYZ");
  gStyle->SetNdivisions(510,"XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);
}


//!  \brief Adjust y-axis range
//!
//!  Sets the y-axis range of \p h from
//!  <tt> c1 * min</tt> to <tt> c2 * max</tt>,
//!  where \p min and \p max are the minimal and
//!  the maximal bin non-zero content, respectively.
//!  If <tt>min < minLimit</tt>, \p minLimit is used
//!  instead as minimum.
// --------------------------------------------------
void TControlPlots::setYRange(TH1F * h, double c1, double c2, double minLimit) const {
  double min = 0.;
  double max = 0.;
  findYRange(h,min,max);
  min *= c1;
  max *= c2;
  if( min < minLimit ) min = minLimit;
  h->GetYaxis()->SetRangeUser( min, max );
}



//!  Write TObjects to ROOT file
//!
//!  Write all TObjects in \p obj to the file \p outFile_ into the
//!  directory 'outFile_:/dir'. If 'outFile_:/dir' does not exist,
//!  it is created first.
//---------------------------------------------------------------
void TControlPlots::writeToRootFile(std::vector<TObject*> obj, std::string dir)
{
  std::string directory = outFile_->GetName();
  directory += ":";
  gDirectory->cd(directory.c_str());
  directory += "/";
  directory += dir;
  bool dirExists = gDirectory->GetDirectory(directory.c_str());
  if( !dirExists )
    {
      gDirectory->mkdir(dir.c_str());
    }
  gDirectory->cd(directory.c_str());
  for(std::vector<TObject*>::const_iterator it = obj.begin(); it < obj.end(); it++)
    {
      int ok = gDirectory->WriteTObject( *it );
      if( !ok ) std::cerr << "Error writing object '" << (*it)->GetName() << "' to file." << std::endl;
    }
}



// -------------------------------------------------------------
TControlPlots::Binning::Binning(const std::vector<double>& binEdgesX, const std::vector<double>& binEdgesY)
{
  assert( binEdgesX.size() > 1 );
  assert( binEdgesY.size() > 1 );

  for(unsigned int i = 0; i < binEdgesX.size(); i++)
    {
      if( i > 0 ) assert( binEdgesX.at(i) > edgesX_.at(i-1) );
      edgesX_.push_back(binEdgesX.at(i));
    }
  for(unsigned int i = 0; i < binEdgesY.size(); i++)
    {
      if( i > 0 ) assert( binEdgesY.at(i) > edgesY_.at(i-1) );
      edgesY_.push_back(binEdgesY.at(i));
    }
}



// -------------------------------------------------------------
int TControlPlots::Binning::bin(int ix, int iy) const
{
  int bin = -1;
  if( ix >= 0 && ix < nBinsX() && iy >= 0 && iy < nBinsY() )
    {
      bin = ix + iy*nBinsX();
    }
  return bin;
}



// -------------------------------------------------------------
int TControlPlots::Binning::iX(double x) const
{
  int ix = -1;                     // Underflow
  if( x > edgesX_.at(nBinsX()) )   // Overflow
    {
      ix = nBinsX();
    }
  else if( x > edgesX_.at(0) )
    {
      ix = 0;
      while( x > edgesX_.at(ix+1) ) ix++;
    }
  return ix;
}



// -------------------------------------------------------------
int TControlPlots::Binning::iY(double y) const
{
  int iy = -1;                     // Underflow
  if( y > edgesY_.at(nBinsY()) )   // Overflow
    {
      iy = nBinsY();
    }
  else if( y > edgesY_.at(0) )
    {
      iy = 0;
      while( y > edgesY_.at(iy+1) ) iy++;
    }
  return iy;
}



// -------------------------------------------------------------
void TControlPlots::Binning::print() const
{
  for(int i = 0; i < nBins(); i++)
    {
      std::cout << i << ":  " << iX(i) << " (" << xLow(i) << ", " << xUp(i) << "),  " << iY(i) << " (" << yLow(i) << ", " << yUp(i) << ")" << std::endl;
    }
}


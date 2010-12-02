// $Id: plotTailUncertainties.C,v 1.2 2010/11/30 14:34:52 mschrode Exp $

#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"

#include "../sampleTools/BinningAdmin.h"
#include "../util/utils.h"
#include "../util/FileOps.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"



class Variation {
public:
  Variation(const TString &fileName, const TString &id, const TString &label);

  TString label() const { return label_; }
  double val(unsigned int etaBin, unsigned int bin) const { return hScales_.at(etaBin)->GetBinContent(1+bin); }
  double err(unsigned int etaBin, unsigned int bin) const { return hScales_.at(etaBin)->GetBinError(1+bin); }

private:
  const TString label_;
  util::HistVec hScales_;
};

Variation::Variation(const TString &fileName, const TString &id, const TString &label) 
  : label_(label) {
    hScales_ = util::FileOps::readHistVec(fileName,"Ratio_hNTailIntData_Eta",id+":hScaleFactorsInt_Eta");
}



void plotTailUncertainties() {
  std::cout << "Setting up parameters" << std::endl;

  util::StyleSettings::presentationNoTitle();
  sampleTools::BinningAdmin* binAdm = new sampleTools::BinningAdmin("BinningAdmin.cfg");
  
  TString nomName = "results/Tails_nSCore20_nSTail30_ptSoft20_PF_.root";


  std::vector<Variation*> vars;
  vars.push_back(new Variation("results/Tails_nSCore20_nSTail25_ptSoft20_PF_.root","TailStart25","Tail Start 2.5#sigma"));
  vars.push_back(new Variation("results/Tails_nSCore20_nSTail35_ptSoft20_PF_.root","TailStart35","Tail Start 3.5#sigma"));
  vars.push_back(new Variation("results/Tails_nSCore15_nSTail30_ptSoft20_PF_.root","TailStart35","Core 1.8#sigma"));
  vars.push_back(new Variation("results/Tails_nSCore25_nSTail30_ptSoft20_PF_.root","TailStart35","Core 2.2#sigma"));
  vars.push_back(new Variation("results/Tails_nSCore20_nSTail30_ptSoft10_PF_.root","TailStart35","p_{T,3} < 0.1#upoint#bar{p^{ave}_{T}}"));
  vars.push_back(new Variation("results/Tails_nSCore20_nSTail30_ptSoft30_PF_.root","TailStart35","p_{T,3} < 0.3#upoint#bar{p^{ave}_{T}}"));


  TString outName = "Tails_";
  TString jetAlgo;
  if( nomName.Contains("Calo") ) {
    jetAlgo = "AK5 Calo-Jets";
    outName += "Calo_";
  } else if( nomName.Contains("PF") ) {
    jetAlgo = "AK5 PF-Jets";
    outName += "PF_";
  }

  TFile outFile(outName+".root","RECREATE");

  util::HistVec hScales = util::FileOps::readHistVec(nomName,"Ratio_hNTailIntData_Eta","hScaleFactorsInt_Eta");


  // Plot variations per eta and pt (integrated) bin
  for(unsigned int etaBin = 0; etaBin < 1; ++etaBin) {//binAdm->nEtaBins(); ++etaBin) {
    std::vector<double> errLow;
    std::vector<double> errHigh;
    for(int ptBin = 1; ptBin <= hScales[etaBin]->GetNbinsX(); ++ptBin) {

      TH1* hNom = new TH1D("hNom",";;Scaling Factor",vars.size()+1,-0.5,vars.size()+0.5);
      hNom->SetNdivisions(-1.*vars.size());
      TH1* hVars = static_cast<TH1D*>(hNom->Clone("hVars"));
      hVars->SetMarkerStyle(20);
      hNom->SetLineStyle(2);

      double minVar = hScales[etaBin]->GetBinContent(ptBin);
      double maxVar = hScales[etaBin]->GetBinContent(ptBin);
      hVars->SetBinContent(1,hScales[etaBin]->GetBinContent(ptBin));
      hVars->SetBinError(1,hScales[etaBin]->GetBinError(ptBin));
      hNom->SetBinContent(1,hScales[etaBin]->GetBinContent(ptBin));
      hNom->SetBinError(1,0.);
      hNom->GetXaxis()->SetBinLabel(1,"Nominal");
      for(unsigned int i = 0; i < vars.size(); ++i) {
	double var = vars[i]->val(etaBin,ptBin-1);
	hVars->SetBinContent(2+i,var);
	hVars->SetBinError(2+i,vars[i]->err(etaBin,ptBin-1));
	hNom->SetBinContent(2+i,hScales[etaBin]->GetBinContent(ptBin));
	hNom->SetBinError(2+i,0.);
	hNom->GetXaxis()->SetBinLabel(2+i,vars[i]->label());

	if( var < minVar ) minVar = var;
	if( var > maxVar ) maxVar = var;
      }
      hNom->LabelsOption("v");
      errLow.push_back(minVar);
      errHigh.push_back(maxVar);

      //std::cout << "PtBin " << ptBin << ": " << minVar << " - " << maxVar << std::endl;

      hNom->GetYaxis()->SetRangeUser(0.5,2.);
      double min = hScales[etaBin]->GetXaxis()->GetBinLowEdge(ptBin);
      double max = hScales[etaBin]->GetXaxis()->GetBinUpEdge(ptBin);
      TPaveText *txt = util::LabelFactory::createPaveText(3,-0.6);
      txt->AddText(jetAlgo);
      txt->AddText(util::toTString(binAdm->etaMin(etaBin))+" < |#eta| < "+util::toTString(binAdm->etaMax(etaBin)));
      txt->AddText(util::toTString(min)+" < p^{ave}_{T} < "+util::toTString(max)+" GeV");

      gStyle->SetPadTickX(0);
      TCanvas *can = new TCanvas("can","Variation",500,500);
      can->cd();
      hNom->Draw("HIST");
      hVars->Draw("PE1same");
      txt->Draw("same");
      can->SaveAs(outName+"Variations_Eta"+util::toTString(etaBin)+"_PtInt"+util::toTString(ptBin-1)+".eps");
      delete txt;
      delete hNom;
      delete hVars;
      delete can;
      gStyle->SetPadTickX(1);
    } // End of loop over integrated pt bins

    // Plot nominal scaling factors + uncertainty band
    TCanvas *can = new TCanvas("can","Uncertainty",500,500);
    can->cd();
    std::vector<double> pt;
    std::vector<double> ptEl;
    std::vector<double> ptEr;
    std::vector<double> fac;
    for(int bin = 1; bin <= hScales[etaBin]->GetNbinsX(); ++bin) {
      pt.push_back(hScales[etaBin]->GetBinCenter(bin));
      ptEl.push_back(hScales[etaBin]->GetBinCenter(bin)-hScales[etaBin]->GetXaxis()->GetBinLowEdge(bin));
      ptEr.push_back(hScales[etaBin]->GetXaxis()->GetBinUpEdge(bin)-hScales[etaBin]->GetBinCenter(bin));
      fac.push_back(hScales[etaBin]->GetBinContent(bin));
      errLow.at(bin-1) = fac.back()-errLow.at(bin-1);
      errHigh.at(bin-1) = errHigh.at(bin-1)-fac.back();
    }
    TGraphAsymmErrors *band = new TGraphAsymmErrors(pt.size(),&(pt.front()),&(fac.front()),&(ptEl.front()),&(ptEr.front()),&(errLow.front()),&(errHigh.front()));
    band->SetFillColor(kGray);
    band->SetLineColor(0);
//     band->SetFillColor(1);
//     band->SetFillStyle(3013);
    
    TPaveText *txt = util::LabelFactory::createPaveText(2,-0.5);
    txt->AddText(jetAlgo);
    txt->AddText(util::toTString(binAdm->etaMin(etaBin))+" < |#eta| < "+util::toTString(binAdm->etaMax(etaBin)));

    hScales[etaBin]->UseCurrentStyle();
    hScales[etaBin]->SetMarkerStyle(20);
    hScales[etaBin]->SetYTitle("Scaling Factor");
    hScales[etaBin]->GetYaxis()->SetRangeUser(0.5,2.);
    hScales[etaBin]->Draw();
    band->Draw("E2same");
    hScales[etaBin]->Draw("PE1same");
    txt->Draw("same");
    TLegend* leg = util::LabelFactory::createLegendCol(2,0.4);
    leg->AddEntry(hScales[etaBin],"Stat. uncert","L");
    leg->AddEntry(band,"Syst. uncert","F");
    leg->Draw("same");
    can->SetLogx();
    can->SaveAs(outName+"Variations_Eta"+util::toTString(etaBin)+".eps");

    // Sum up systematic and statistical errors
    TH1* hUp = static_cast<TH1D*>(hScales[etaBin]->Clone("hScaleFactorsIntUp_Eta"+util::toTString(etaBin)));
    hUp->SetMarkerStyle(1);
    hUp->Reset();
    TH1* hDown = static_cast<TH1D*>(hUp->Clone("hScaleFactorsIntDown_Eta"+util::toTString(etaBin)));
    for(int bin = 1; bin <= hScales[etaBin]->GetNbinsX(); ++bin) {
      double nom = hScales[etaBin]->GetBinContent(bin);
      double stat = hScales[etaBin]->GetBinError(bin);
      double eup = band->GetEYhigh()[bin-1];
      double edo = band->GetEYlow()[bin-1];
      hUp->SetBinContent(bin,nom+sqrt(stat*stat+eup*eup));
      hDown->SetBinContent(bin,nom-sqrt(stat*stat+edo*edo));
    }
    can->cd();
    hScales[etaBin]->Draw("P");
    hUp->Draw("HISTsame");
    hDown->Draw("HISTsame");    
    can->SaveAs(outName+"SumErr_Eta"+util::toTString(etaBin)+".eps");

    outFile.WriteTObject(hScales[etaBin]);
    outFile.WriteTObject(hUp);
    outFile.WriteTObject(hDown);

    delete leg;
    delete txt;
    delete band;
    delete hUp;
    delete hDown;
    delete can;
  } // End of loop over eta bins

  outFile.Close();
}

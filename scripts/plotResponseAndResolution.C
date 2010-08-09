// $Id: plotResponseAndResolution.C,v 1.2 2010/08/04 09:29:51 mschrode Exp $

//!  Fit mean response and resolution from
//!  Kalibri::ControlPlotsJetSmearing

#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TString.h"

#include "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/util/HistOps.h"
#include "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/util/LabelFactory.h"
#include "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/util/StyleSettings.h"


unsigned int numHists_ = 0;


std::vector<TH1*> readHistograms(const std::vector<TString> &fileNames, const TString &histName, const TString &xTitle, const TString &xUnit, const TString &yTitle) {
  std::vector<TH1*> hists(fileNames.size(),0);
  for(size_t i = 0; i < fileNames.size(); ++i) {
    TFile file(fileNames[i],"READ");
    file.GetObject(histName,hists[i]);
    if( hists[i] == 0 ) {
      std::cerr << "ERROR getting histogram '" << histName << "' from file '" << fileNames[i] << "'\n";
      exit(1);
    }
    hists[i]->SetDirectory(0);
    TString name = histName;
    name += i;
    name += "_";
    name += numHists_;
    hists[i]->SetName(name);
    util::HistOps::setAxisTitles(hists[i],xTitle,xUnit,yTitle);
      util::HistOps::setStyleColor(hists[i],1+i/2);
    //hists[i]->SetMarkerStyle(20+i);
    if( i%2 == 0 ) hists[i]->SetMarkerStyle(20+i/2);
    else hists[i]->SetMarkerStyle(24+i/2); 
    util::HistOps::setYRange(hists[i],(int)(fileNames.size()+1));
    hists[i]->GetXaxis()->SetRangeUser(20.,3500.);
    numHists_++;
  }
  return hists;
}

  
void plotResponse(const std::vector<TString> &fileNames, const std::vector<TString> &labels, const TString &outName = "") {
  assert( fileNames.size() == labels.size() );

  std::vector<TH1*> hResp = readHistograms(fileNames,"MeanResp_hRespGauss","p^{gen}_{T}","GeV","Mean response");
  TLegend *leg = util::LabelFactory::createLegendCol(fileNames.size(),0.7);
  TH1 *hFrame = util::HistOps::createRatioFrame(hResp[0],"Mean response",0.8,1.4);
  TCanvas *can = new TCanvas("canResp","Mean Response",500,500);
  can->cd();
  hFrame->Draw();
  for(size_t i = 0; i < hResp.size(); ++i) {
    hResp[i]->Draw("PE1same");
    TString entry = "";
    entry += i;
    entry += ": "+labels[i];
    leg->AddEntry(hResp[i],entry,"P");
  }
  leg->Draw("same");
  can->SetLogx();

  if( outName != "" ) {
    can->SaveAs(outName+".eps","eps");
  }
}


void plotResolution(const std::vector<TString> &fileNames, const std::vector<TString> &labels, bool fit, const TString &outName = "") {
  assert( fileNames.size() == labels.size() );

  std::vector<TH1*> hReso = readHistograms(fileNames,"MeanResp_hResoGauss","p^{gen}_{T}","GeV","Resolution");
  TLegend *leg = util::LabelFactory::createLegendCol(fileNames.size(),0.7);
  TH1 *hFrame = util::HistOps::createFrame(hReso[0],"Resolution",0.,0.4);
  std::vector<TF1*> fits(hReso.size(),0);
  for(size_t i = 0; i < hReso.size(); ++i) {
    TString name = "Fit";
    name += i;
    fits[i] = new TF1(name,"sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",20.,3500.);
// 		      hReso[i]->GetXaxis()->GetBinLowEdge(1),
// 		      hReso[i]->GetXaxis()->GetBinUpEdge(hReso[i]->GetNbinsX()));
    fits[i]->SetLineWidth(1);
    fits[i]->SetLineColor(hReso[i]->GetLineColor());
    fits[i]->SetParameter(0,2.);
    fits[i]->SetParameter(1,1.3);
    fits[i]->SetParameter(2,0.04);
    if( fit ) {
      for(int bin = 1; bin <= hReso[i]->GetNbinsX(); ++bin) {
	hReso[i]->SetBinError(bin,0.01);
	if( hReso[i]->GetBinContent(bin) > 0.38 ) hReso[i]->SetBinError(bin,100.);
      }
      hReso[i]->Fit(fits[i],"INQR");
    }
  }

  if( fit ) {
    std::cout << "\n\n********** INTERPOLATED RESOLUTION **********\n";
    std::cout << "\\begin{tabular}{lccc}\n\\toprule\n";
    std::cout << " & $b_{0}\\,(\\text{Ge}\\kern-0.06667em\\text{V})$ & $b_{1}\\,(\\sqrt{\\text{Ge}\\kern-0.06667em\\text{V}})$ & $b_{2}$ \\\\\n\\midrule\n";
    for(size_t i = 0; i < fits.size(); ++i) {
      std::cout << "$" << i << "$" << std::flush;
      for(int p = 0; p < 3; ++p) {
	std::cout << " & $" << fits[i]->GetParameter(p) << " \\pm " << fits[i]->GetParError(p) << "$" << std::flush;
      }
      std::cout << " \\\\ \n";
    }
    std::cout << "\\bottomrule\n\\end{tabular}\n";
    std::cout << "\n\n";
  }

  TCanvas *can = new TCanvas("canReso","Resolution",500,500);
  can->cd();
  hFrame->Draw();
  for(size_t i = 0; i < hReso.size(); ++i) {
    hReso[i]->Draw("PE1same");
    TString entry = "";
    entry += i;
    entry += ": "+labels[i];
    leg->AddEntry(hReso[i],entry,"P");
    if( fit ) fits[i]->Draw("same");
  }
  leg->Draw("same");
  can->SetLogx();
  if( outName != "" ) {
    TString name = outName;
    if( fit ) name += "_Fits";
    name += ".eps";
    can->SaveAs(name,"eps");
  }


  if( fit ) {
    for(size_t i = 0; i < hReso.size(); ++i) {
      TString name = "canReso";
      name += i;
      TString title = "Resolution Fit ";
      title += i;
      TCanvas *can2 = new TCanvas(name,title,500,500);
      can2->cd();
      hFrame->Draw();
      hReso[i]->Draw("PE1same");
      fits[i]->Draw("same");
      TLegend *leg2 = util::LabelFactory::createLegendCol(1,0.7);
      TString entry = "";
      entry += i;
      entry += ": "+labels[i];
      leg2->AddEntry(hReso[i],entry,"P");
      leg2->Draw("same");
      can2->SetLogx();  
      if( outName != "" ) {
	name = outName+"_Fit";
	name += i;
	name += ".eps";
	can2->SaveAs(name,"eps");
      }
    }
  }
}


void plotResolutionScale1(const std::vector<TString> &fileNames, const std::vector<TString> &labels) {
  assert( fileNames.size() == labels.size() );

  std::vector<TH1*> hReso = readHistograms(fileNames,"MeanResp_hResoGauss","p^{gen}_{T}","GeV","Resolution");
  std::vector<TH1*> hResp = readHistograms(fileNames,"MeanResp_hRespGauss","p^{gen}_{T}","GeV","Mean response");
  for(size_t i = 0; i < hReso.size(); ++i) {
    assert( hReso[i]->GetNbinsX() == hResp[i]->GetNbinsX() );
    for(int ptBin = 1; ptBin <= hReso[i]->GetNbinsX(); ++ptBin) {
      hReso[i]->SetBinContent(ptBin,hReso[i]->GetBinContent(ptBin)*hResp[i]->GetBinContent(ptBin));
    }
  }

  TLegend *leg = util::LabelFactory::createLegendCol(fileNames.size(),0.7);
  TH1 *hFrame = util::HistOps::createFrame(hReso[0],"Resolution",0.,0.4);
  TCanvas *can = new TCanvas("canResoScale1","Resolution Scale 1",500,500);
  can->cd();
  hFrame->Draw();
  for(size_t i = 0; i < hReso.size(); ++i) {
    hReso[i]->Draw("PE1same");
    leg->AddEntry(hReso[i],labels[i],"P");
  }
  leg->Draw("same");
  can->SetLogx();
}



void plotResponseAndResolution() {
  util::StyleSettings::presentationNoTitle();
  
  TString base = "~/results/ResolutionFit/Spring10QCDDiJet_MCTruthResponse0010-3500_";

  std::vector<TString> fileNames;
  std::vector<TString> labels;

//   fileNames.push_back(base+"Sigma_FreeMu_CaloOrdered_Rel3rdJet100.root");
//   labels.push_back("Core");
//   fileNames.push_back(base+"Core_Sigma_FreeMu_CaloOrdered_3rdJet010.root");
//   labels.push_back("Core,  p_{T,3} < 0.1 <p^{gen}_{T;1,2}>");
//   fileNames.push_back(base+"All_Sigma_FreeMu_CaloOrdered_Rel3rdJet010.root");
//   labels.push_back("All,  p_{T,3} < 0.1 <p^{gen}_{T;1,2}>");

  fileNames.push_back(base+"Core_Sigma_FreeMu_CaloOrdered_Rel3rdJet100.root");
  labels.push_back("Calo ordered");
  fileNames.push_back(base+"Core_Sigma_FreeMu_GenOrdered_Rel3rdJet100.root");
  labels.push_back("Gen ordered");

  fileNames.push_back(base+"Core_Sigma_FreeMu_CaloOrdered_Rel3rdJet010.root");
  labels.push_back("Calo ordered, p^{rel}_{T,3} < 0.1");
  fileNames.push_back(base+"Core_Sigma_FreeMu_GenOrdered_Rel3rdJet010.root");
  labels.push_back("Gen ordered, p^{rel}_{T,3} < 0.1");

  fileNames.push_back(base+"Core_Sigma_FreeMu_CaloOrdered_3rdJet010.root");
  labels.push_back("Calo ordered, p_{T,3} < 0.1 <p^{gen}_{T;1,2}>");
  fileNames.push_back(base+"Core_Sigma_FreeMu_GenOrdered_3rdJet010.root");
  labels.push_back("Gen ordered, p_{T,3} < 0.1 <p^{gen}_{T;1,2}>");


  plotResponse(fileNames,labels,"Spring10_QCDDiJet_MeanResp_DijetCuts_Core");
  plotResolution(fileNames,labels,false,"Spring10_QCDDiJet_Resolution_DijetCuts_Core");
  //  plotResolutionScale1(fileNames,labels);
}

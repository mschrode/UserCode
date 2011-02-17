// $Id: fitMCTruth.C,v 1.11 2011/01/30 19:36:27 mschrode Exp $

//!  Fit mean response and resolution from
//!  Kalibri::ControlPlotsJetSmearing

#define UTILS_AS_HEADER_FILE

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TString.h"

#include "../util/FileOps.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


unsigned int gCOLOR_GAUSS = kRed;
unsigned int gMARKER_GAUSS = 20;
unsigned int gCOLOR_RMS = kBlue;
unsigned int gMARKER_RMS = 24;


TH1* createTH1FromTH2(const TH2* h2, const TString &name, const TString &xTitle, const TString &xUnit, const TString &yTitle) {
  TH1* h = new TH1D(name,"",h2->GetNbinsX(),h2->GetXaxis()->GetXbins()->GetArray());
  util::HistOps::setAxisTitles(h,xTitle,xUnit,yTitle);

  return h;
}



void fitProfile(const TString &fileName, double nSigCore, TH1* &hMean, TH1* &hMeanGauss, TH1* &hRMS, TH1* &hSigmaGauss, TH1* &hChi2Ndof, const TString &outNamePrefix, double etaMin, double etaMax, const TString &rootOutFileName, const TString &histNameSuffix) {
  std::cout << "Fitting profiles" << std::endl;

  // Get 2D histogram resp vs pt from file
  std::cout << "  Getting 2D histogram 'response vs ptGen' from file '" << fileName << "'" << std::endl;
  TH2* hRespVsPt = util::FileOps::readTH2(fileName,"MeanResp_hRespVsPtGen","hRespVsPt_File"+histNameSuffix);


  // Creating products
  std::cout << "  Creating result histograms" << std::endl;
  hMean = createTH1FromTH2(hRespVsPt,"hMean"+histNameSuffix,"p^{gen}_{T}","GeV","Response Mean");
  hMean->SetMarkerStyle(gMARKER_RMS);
  hMean->SetMarkerColor(gCOLOR_RMS);
  hMean->SetLineColor(gCOLOR_RMS);
  hMeanGauss = createTH1FromTH2(hRespVsPt,"hMeanGauss"+histNameSuffix,"p^{gen}_{T}","GeV","Response Gauss Mean");
  hMeanGauss->SetMarkerStyle(gMARKER_GAUSS);
  hMeanGauss->SetMarkerColor(gCOLOR_GAUSS);
  hMeanGauss->SetLineColor(gCOLOR_GAUSS);
  hRMS = createTH1FromTH2(hRespVsPt,"hRMS"+histNameSuffix,"p^{gen}_{T}","GeV","Response Standard Deviation");
  hRMS->SetMarkerStyle(gMARKER_RMS);
  hRMS->SetMarkerColor(gCOLOR_RMS);
  hRMS->SetLineColor(gCOLOR_RMS);
  hSigmaGauss = createTH1FromTH2(hRespVsPt,"hSigmaGauss"+histNameSuffix,"p^{gen}_{T}","GeV","Response Gaussian Width");
  hSigmaGauss->SetMarkerStyle(gMARKER_GAUSS);
  hSigmaGauss->SetMarkerColor(gCOLOR_GAUSS);
  hSigmaGauss->SetLineColor(gCOLOR_GAUSS);

  hChi2Ndof = static_cast<TH1D*>(hMean->Clone("hChi2Ndof"));
  hChi2Ndof->SetYTitle("#chi^{2} / ndof");
  hChi2Ndof->SetMarkerColor(kBlack);


  // Fill response distributions per pt bin
  std::cout << "  Filling response distributions per ptGen bin" << std::endl;
  util::HistVec hResp;
  util::HistOps::fillSlices(hRespVsPt,hResp,"hMCTruthResp"+histNameSuffix);
  for(util::HistItConst hIt = hResp.begin(); hIt != hResp.end(); ++hIt) {
    util::HistOps::setAxisTitles(*hIt,"Response","","jets");
    (*hIt)->UseCurrentStyle();
    (*hIt)->SetMarkerStyle(20);
  }

  std::cout << "  Writing response distributions to file" << std::endl;
  TFile outFile(rootOutFileName,"UPDATE");
  for(util::HistItConst hIt = hResp.begin(); hIt != hResp.end(); ++hIt) {
    outFile.WriteTObject(*hIt);
  }
  outFile.Close();


  // Fitting response distributions
  std::cout << "  Fitting response distributions per ptGen bin" << std::endl;
  std::vector<TF1*> fGauss;
  util::HistVec hRatio;
  int bin = 1;
  for(util::HistItConst hIt = hResp.begin(); hIt != hResp.end(); ++hIt, ++bin) {
    double width = 0.;
    double widthErr = 0.;
    double rms = 0.;
    double rmsErr = 0.;
    TF1* fit = 0;
    if( util::HistOps::fitCoreWidth(*hIt,nSigCore,fit,width,widthErr,rms,rmsErr) ) {
      hMean->SetBinContent(bin,(*hIt)->GetMean());
      hMean->SetBinError(bin,(*hIt)->GetMeanError());
      hMeanGauss->SetBinContent(bin,fit->GetParameter(1));
      hMeanGauss->SetBinError(bin,fit->GetParError(1));
      hRMS->SetBinContent(bin,rms);
      hRMS->SetBinError(bin,rmsErr);
      hSigmaGauss->SetBinContent(bin,width);
      hSigmaGauss->SetBinError(bin,widthErr);
      
      hChi2Ndof->SetBinContent(bin,fit->GetChisquare() / fit->GetNDF());

    } else {
      fit = new TF1("fGauss"+histNameSuffix+"_Bin"+util::toTString(bin),"gaus");
    }
    fit->SetLineWidth(1);
    fit->SetLineColor(gCOLOR_GAUSS);    
    fGauss.push_back(fit);

    // Ratio to Gaussian fit
    hRatio.push_back(util::HistOps::createRatioPlot(*hIt,fit));    
  }


  // Plotting resonse distributions and fits
  std::cout << "  Plotting resonse distributions and fits" << std::endl;
  bin = 1;
  for(util::HistItConst hIt = hResp.begin(); hIt != hResp.end(); ++hIt, ++bin) {
    double ptMin = hRespVsPt->GetXaxis()->GetBinLowEdge(bin);
    double ptMax = hRespVsPt->GetXaxis()->GetBinUpEdge(bin);
    TPaveText* label = util::LabelFactory::createPaveText(3);
    label->AddText("CMS Simulation,  #sqrt{s} = 7 TeV");
    TString jetLabel = "Anti-k_{T} (R=0.5) ";
    if( outNamePrefix.Contains("Calo") ) jetLabel += "Calo Jets";
    else if( outNamePrefix.Contains("PF") ) jetLabel += "PF Jets";
    else if( outNamePrefix.Contains("JPT") ) jetLabel += "JPT Jets";
    label->AddText(jetLabel);
    label->AddText(util::toTString(etaMin)+" < |#eta| < "+util::toTString(etaMax)+",  "+util::toTString(ptMin,0)+" < p^{gen}_{T} < "+util::toTString(ptMax,0)+" GeV");


    TCanvas* can = new TCanvas("can","Response PtBin "+util::toTString(bin),500,500);

    can->cd();
    util::HistOps::setYRange(*hIt,2,true);
    //(*hIt)->GetXaxis()->SetRangeUser(0.,2.);
    (*hIt)->Draw("PE1same");
    fGauss.at(bin-1)->Draw("same");
    label->Draw("same");
    can->SetLogy();
    can->SaveAs(outNamePrefix+"_Log_PtBin"+util::toTString(bin-1)+".eps","eps");

    delete can;
    can = util::HistOps::createRatioTopCanvas();
    TPad *bPad = util::HistOps::createRatioBottomPad();
    TH1 *tFrame = util::HistOps::createRatioTopHist((*hIt));
    TH1 *bFrame = util::HistOps::createRatioBottomFrame((*hIt),"Response","",0.91,1.09);
    can->cd();
    double min = 10.;
    double max = 0.;
    util::HistOps::findYRange((*hIt),4,min,max);
    tFrame->GetYaxis()->SetRangeUser(0.,max);
    std::vector<TLine*> lines;
    for(int i = 1; i <= 4; ++i) {
      double mean = fGauss.at(bin-1)->GetParameter(1);
      double sig = fGauss.at(bin-1)->GetParameter(2);
      lines.push_back(new TLine(mean+i*sig,0.,mean+i*sig,max));
      lines.push_back(new TLine(mean-i*sig,0.,mean-i*sig,max));
    }

    tFrame->GetXaxis()->SetRange(tFrame->FindBin(0.71),tFrame->FindBin(1.31));
    bFrame->GetXaxis()->SetRange(bFrame->FindBin(0.71),bFrame->FindBin(1.31));
    tFrame->Draw("PE1");
    fGauss.at(bin-1)->Draw("same");
    for(size_t i = 0; i < lines.size(); ++i) {
      lines.at(i)->SetLineColor(kBlue);
      //lines.at(i)->Draw("same");
    }
    label->Draw("same");
    bPad->Draw();
    bPad->cd();
    bFrame->Draw();
    hRatio.at(bin-1)->Draw("PE1same");
    can->SaveAs(outNamePrefix+"_Linear_PtBin"+util::toTString(bin-1)+".eps","eps");

    delete label;
    delete can;
  }


  // Creating tex slides
  std::cout << "  Creating LaTeX slides" << std::endl;

  std::ofstream oFile(outNamePrefix+".tex");

  unsigned int nPtBins = hResp.size();
  
  oFile << "\n\n\n% ----- MC Truth Response ---------------------------" << std::endl;
  unsigned int nSlides = nPtBins/3;
  if( nPtBins%3 > 0 ) nSlides++;
  unsigned int ptBin = 0;
  for(unsigned int slide = 0; slide < nSlides; ++slide) {
    oFile << "\n% --------------------------------------------------\n";
    oFile << "\\begin{frame}\n";
    oFile << "  \\frametitle{MC Truth Response (" << slide+1 << "/" << nSlides << ")}\n";
    oFile << "  \\begin{columns}[T] \n";
    for(int col = 0; col < 3; ++col, ++ptBin) {
      oFile << "    \\begin{column}{0.3333\\textwidth} \n";
      oFile << "    \\centering\n";
      if( ptBin < nPtBins ) {
	oFile << "      \\includegraphics[width=\\textwidth]{" << outNamePrefix+"_Linear_PtBin"+util::toTString(ptBin)+"}\\\\ \n";
	oFile << "      \\includegraphics[width=\\textwidth]{" << outNamePrefix+"_Log_PtBin"+util::toTString(ptBin)+"}\\\\ \n";
      }
      oFile << "    \\end{column} \n";
    }
    oFile << "  \\end{columns} \n";
    oFile << "\\end{frame} \n";
  }
  oFile.close();
}



TF1* fitResolution(TH1* h, const TString &name, double min, double max) {
//   TF1* fit = new TF1(name,"sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",min,max);
//   fit->SetParameter(0,4.);
//   fit->SetParameter(1,1.);
//   fit->SetParameter(2,0.01);

  TF1* fit = new TF1(name,"sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))",min,max);
  fit->SetParameter(0,-0.34);
  fit->SetParameter(1,0.3);
  fit->FixParameter(2,0.);
  fit->SetParameter(3,0.47);

  fit->SetLineWidth(1);
  fit->SetLineColor(kRed);
  h->Fit(fit,"BINR");


  return fit;
}



void fitMCTruth(const TString &fileName, double nSigCore, double minPt) {
  util::StyleSettings::paperNoTitle();

  TString jetAlgo;
  if( fileName.Contains("PF") ) jetAlgo = "PF";
  else if( fileName.Contains("Calo") ) jetAlgo = "Calo";
  else if( fileName.Contains("JPT") ) jetAlgo = "JPT";
  
  TString outNamePrefix = "MCTruthResponse_"+jetAlgo;
  double etaMin = 0.;
  double etaMax = 0.;
  int etaBin = 0;
  if( fileName.Contains("Eta0") ) {
    outNamePrefix += "_Eta0";
    etaMin = 0.;
    etaMax = 1.1;
    etaBin = 0;
  } else if( fileName.Contains("Eta1") ) {
    outNamePrefix += "_Eta1";
    etaMin = 1.1;
    etaMax = 1.7;
    etaBin = 1;
  } else if( fileName.Contains("Eta2") ) {
    outNamePrefix += "_Eta2";
    etaMin = 1.7;
    etaMax = 2.3;
    etaBin = 2;
  } else if( fileName.Contains("Eta3") ) {
    outNamePrefix += "_Eta3";
    etaMin = 2.3;
    etaMax = 5.;
    etaBin = 3;
  }
  TString histNameSuffix = "";
  if( fileName.Contains("Dijets") ) {
    outNamePrefix += "_Dijets";
    histNameSuffix += "_Dijets";
  } else if( fileName.Contains("DeltaR10") ) {
    outNamePrefix += "_DeltaR10";
    histNameSuffix += "_DeltaR10"; 
  } else if( fileName.Contains("DeltaR20") ) {
    outNamePrefix += "_DeltaR20";
    histNameSuffix += "_DeltaR20"; 
  } else if( fileName.Contains("DeltaR25") ) {
    outNamePrefix += "_DeltaR25";
    histNameSuffix += "_DeltaR25";
  }

  TString jetLabel = "Anti-k_{T} (R=0.5) ";
  if( jetAlgo == "Calo" ) jetLabel += "Calo Jets";
  else if( jetAlgo == "PF" ) jetLabel += "PF Jets";
  else if( jetAlgo == "JPT" ) jetLabel += "JPT Jets";
  jetLabel += ", "+util::toTString(etaMin)+" < |#eta| < "+util::toTString(etaMax);


  // Root out file name
  TString rootOutFileName = "MCTruthResponse_"+jetAlgo+".root";



  // Fit distributions
  TH1* hMean = 0;
  TH1* hMeanGauss = 0;
  TH1* hRMS = 0;
  TH1* hSigmaGauss = 0;
  TH1* hChi2Ndof = 0;
  fitProfile(fileName,nSigCore,hMean,hMeanGauss,hRMS,hSigmaGauss,hChi2Ndof,outNamePrefix,etaMin,etaMax,rootOutFileName,("_Eta"+util::toTString(etaBin)+histNameSuffix));

  // Fit resolution
  TF1* fit = fitResolution(hSigmaGauss,"Res_Eta"+util::toTString(etaBin)+histNameSuffix,minPt,2000.);

  std::cout << std::endl;
  std::cout << "par->setTrueGaussResPar(" << std::flush;
  for(int i = 0; i < fit->GetNpar(); ++i) {
    if( i > 0 ) std::cout << "," << std::flush;
    std::cout << fit->GetParameter(i) << std::flush;
  }
  std::cout << ");" << std::endl;

  std::cout << "\n\n$" << etaMin << " - " << etaMax;
  for(int i = 0; i < fit->GetNpar(); ++i) {
    std::cout << "$ & $" << std::flush;
    std::cout << std::setprecision(4) << fit->GetParameter(i) << " \\pm " << fit->GetParError(i) << std::flush;
  }
  std::cout << "$ \\\\\n\n";

  
  // Plot response and resolution
  TPaveText* label = util::LabelFactory::createPaveText(2);
  label->AddText("CMS Simulation");
  label->AddText(jetLabel);

  TLegend* legMean = util::LabelFactory::createLegendColWithOffset(2,0.7,2);
  legMean->AddEntry(hMean,"Arithmetic Mean","P");
  legMean->AddEntry(hMeanGauss,"Gaussian Mean","P");

  TCanvas* can = new TCanvas("can","",500,500);
  can->cd();
  TH1* hFrame = util::HistOps::createRatioFrame(hMeanGauss,"Mean Response",0.83,1.27);
  hFrame->Draw();
  hMean->Draw("PE1same");
  hMeanGauss->Draw("PE1same");
  label->Draw("same");
  legMean->Draw("same");
  can->SetLogx();
  can->SaveAs(outNamePrefix+"_MeanResponse.eps","eps");

  can->cd();
  hChi2Ndof->GetYaxis()->SetRangeUser(0.,13.);
  hChi2Ndof->GetXaxis()->SetMoreLogLabels();
  hChi2Ndof->SetMarkerStyle(20);
  hChi2Ndof->Draw("P");
  label->Draw("same");
  can->SetLogx();
  can->SaveAs(outNamePrefix+"_Chi2Ndof.eps","eps");

  TLegend* legRes = util::LabelFactory::createLegendColWithOffset(3,0.75,2);
  legRes->AddEntry(hRMS,"Arithmetic Mean","P");
  legRes->AddEntry(hSigmaGauss,"Gaussian Mean","P");
  legRes->AddEntry(fit,"Fit to Gaussian Mean","L");

  TCanvas* canRes = util::HistOps::createRatioTopCanvas();
  TPad *bottomPadRes = util::HistOps::createRatioBottomPad();
  TH1 *topFrameRes = util::HistOps::createRatioTopHist(hRMS);
  TH1 *bottomFrameRes = util::HistOps::createRatioBottomFrame(hRMS,"p^{gen}_{T}","GeV",0.91,1.09);
  canRes->cd();
  topFrameRes->GetXaxis()->SetRange(topFrameRes->FindBin(minPt),topFrameRes->GetNbinsX());
  hSigmaGauss->GetXaxis()->SetRange(hSigmaGauss->FindBin(minPt),hSigmaGauss->GetNbinsX());
  topFrameRes->GetYaxis()->SetRangeUser(1E-3,0.43);
  topFrameRes->GetYaxis()->SetTitle("Resolution");
  topFrameRes->Draw("PE1");
  hSigmaGauss->Draw("PE1same");
  fit->Draw("same");
  label->Draw("same");
  legRes->Draw("same");
  canRes->SetLogx();
  bottomPadRes->Draw();
  bottomPadRes->cd();
  bottomFrameRes->GetXaxis()->SetRange(bottomFrameRes->FindBin(minPt),bottomFrameRes->GetNbinsX());
  bottomFrameRes->GetXaxis()->SetMoreLogLabels();
  bottomFrameRes->Draw();
  TH1* hSigmaGaussRatio = static_cast<TH1D*>(hSigmaGauss->Clone("hSigmaGaussRatio"));
  hSigmaGaussRatio->Divide(fit);
  hSigmaGaussRatio->Draw("PE1same");
  bottomPadRes->SetLogx();
  canRes->SaveAs(outNamePrefix+"_Resolution.eps","eps");

  hSigmaGauss->SetName(("hSigmaGauss_Eta"+util::toTString(etaBin)+histNameSuffix));
  fit->SetName(("fit_Eta"+util::toTString(etaBin)+histNameSuffix));

  TFile outFile("MCTruthResponse_"+jetAlgo+".root","UPDATE");
  outFile.WriteTObject(hSigmaGauss);
  outFile.WriteTObject(fit);
  outFile.Close();
  
  delete label;
  delete legMean;
  delete legRes;
  delete can;
}



// Compare MCTruth of different eta in one plot
void plotMCTruthForDifferentEta(const TString &file, const TString &jetAlgo, double minPt) {
  util::StyleSettings::paperNoTitle();

  util::HistVec hReso = util::FileOps::readHistVec(file,"hSigmaGauss_Eta");
  std::vector<TF1*> fReso = util::FileOps::readTF1Vec(file,"fit_Eta");

  std::vector<double> etaBins;
  etaBins.push_back(0.);
  etaBins.push_back(1.1);
  etaBins.push_back(1.7);
  etaBins.push_back(2.3);
  etaBins.push_back(5.0);

  TLegend* leg = util::LabelFactory::createLegendCol(hReso.size(),0.4);
  for(unsigned int i = 0; i < hReso.size(); ++i) {
    for(int bin = 1; bin <= hReso[i]->GetNbinsX(); ++bin) {
      if( hReso[i]->GetBinError(bin) > 0.1 ) {
	hReso[i]->SetBinContent(bin,-1.);
	hReso[i]->SetBinError(bin,0.);
      }
      if( hReso[i]->GetXaxis()->GetBinLowEdge(bin) < minPt ) {
	hReso[i]->SetBinContent(bin,-1.);
	hReso[i]->SetBinError(bin,0.);
      }
      if( hReso[i]->GetXaxis()->GetBinUpEdge(bin) > 2000. ) {
	hReso[i]->SetBinContent(bin,-1.);
	hReso[i]->SetBinError(bin,0.);
      }
    }
    double maxPt = 0.;
    for(int bin = hReso[i]->GetNbinsX(); bin > 0; --bin) {
      if( hReso[i]->GetBinContent(bin) > 0. ) {
	maxPt = hReso[i]->GetXaxis()->GetBinUpEdge(bin);
	break;
      }
    }
    hReso[i]->UseCurrentStyle();
    hReso[i]->SetMarkerStyle(20+i);
    hReso[i]->SetMarkerColor(util::StyleSettings::color(i));
    hReso[i]->SetLineColor(util::StyleSettings::color(i));
    fReso[i]->SetRange(minPt,maxPt);
    fReso[i]->SetLineColor(util::StyleSettings::color(i));
    fReso[i]->SetLineWidth(1);
    
    if( i < etaBins.size() ) {
      leg->AddEntry(hReso[i],util::toTString(etaBins.at(i))+" < |#eta| < "+util::toTString(etaBins.at(i+1)));
    }
  }

  TH1* hFrame = new TH1D("hFrame",";p^{gen}_{T} (GeV);Resolution",10000,9.,2500.);
  hFrame->GetXaxis()->SetMoreLogLabels();
  hFrame->GetXaxis()->SetNoExponent();
  hFrame->GetYaxis()->SetRangeUser(1E-3,0.38);

  TString jetLabel = "Anti-k_{T} (R=0.5) ";
  if( jetAlgo == "Calo" ) jetLabel += "Calo Jets";
  else if( jetAlgo == "PF" ) jetLabel += "PF Jets";
  else if( jetAlgo == "JPT" ) jetLabel += "JPT Jets";

  TPaveText* label = util::LabelFactory::createPaveText(2,-0.55);
  label->AddText("CMS Simulation, #sqrt{s} = 7 TeV");
  label->AddText(jetLabel);


  TCanvas* can = new TCanvas("can","MC Truth Resolution",500,500);
  can->cd();
  hFrame->Draw();
  for(unsigned int i = 0; i < hReso.size(); ++i) {
    hReso[i]->Draw("PE1same");
  }
  for(unsigned int i = 0; i < fReso.size(); ++i) {
    fReso[i]->Draw("same");
  }
  label->Draw("same");
  leg->Draw("same");
  can->SetLogx();
  can->SaveAs("MCTruthReso_"+jetAlgo+".eps","eps");  
}



void plots(const TString &id, const std::vector<TH1*> &reso, const std::vector<TF1*> &fits, const std::vector<TString> &labels, const TString &jetAlgo, double minPt) {
  util::StyleSettings::paperNoTitle();

  assert( reso.size() == fits.size() );
  assert( reso.size() == labels.size() );

  TLegend* leg = util::LabelFactory::createLegendCol(reso.size(),0.48);
  for(unsigned int i = 0; i < reso.size(); ++i) {
    for(int bin = 1; bin <= reso[i]->GetNbinsX(); ++bin) {
      if( reso[i]->GetBinError(bin) > 0.1 ) {
	reso[i]->SetBinContent(bin,-1.);
	reso[i]->SetBinError(bin,0.);
      }
      if( reso[i]->GetXaxis()->GetBinLowEdge(bin) < minPt ) {
	reso[i]->SetBinContent(bin,-1.);
	reso[i]->SetBinError(bin,0.);
      }
      if( reso[i]->GetXaxis()->GetBinUpEdge(bin) > 2000. ) {
	reso[i]->SetBinContent(bin,-1.);
	reso[i]->SetBinError(bin,0.);
      }
    }
    double maxPt = 0.;
    for(int bin = reso[i]->GetNbinsX(); bin > 0; --bin) {
      if( reso[i]->GetBinContent(bin) > 0. ) {
	maxPt = reso[i]->GetXaxis()->GetBinUpEdge(bin);
	break;
      }
    }
    reso[i]->UseCurrentStyle();
    reso[i]->SetMarkerStyle(20+i);
    reso[i]->SetMarkerColor(util::StyleSettings::color(i));
    reso[i]->SetLineColor(util::StyleSettings::color(i));
    fits[i]->SetRange(minPt,maxPt);
    fits[i]->SetLineColor(util::StyleSettings::color(i));
    fits[i]->SetLineWidth(1);
    
    leg->AddEntry(reso[i],labels[i],"PL");
  }

  TH1* hFrame = new TH1D("hFrame"+id,";p^{gen}_{T} (GeV);Resolution",10000,9.,2500.);
  hFrame->GetXaxis()->SetMoreLogLabels();
  hFrame->GetXaxis()->SetNoExponent();
  hFrame->GetYaxis()->SetRangeUser(1E-3,0.38);

  TString jetLabel = "Anti-k_{T} (R=0.5) ";
  if( jetAlgo == "Calo" ) jetLabel += "Calo Jets";
  else if( jetAlgo == "PF" ) jetLabel += "PF Jets";
  else if( jetAlgo == "JPT" ) jetLabel += "JPT Jets";

  TPaveText* label = util::LabelFactory::createPaveText(3,-0.5);
  label->AddText("CMS Simulation");
  label->AddText("#sqrt{s} = 7 TeV,  |#eta| < 1.1");
  label->AddText(jetLabel);

  TCanvas* can = new TCanvas("can"+id,"MC Truth Resolution "+id,500,500);
  can->cd();
  hFrame->Draw();
  for(unsigned int i = 0; i < reso.size(); ++i) {
    reso[i]->Draw("PE1same");
  }
  for(unsigned int i = 0; i < fits.size(); ++i) {
    fits[i]->Draw("same");
  }
  label->Draw("same");
  leg->Draw("same");
  can->SetLogx();
  can->SaveAs("MCTruthReso_"+id+"_"+jetAlgo+".eps","eps"); 
}



// Compare MCTruth resolutions for different DeltaRMax and dijet selection
// vs no selection
void plotMCTruthForSelections(const TString &file, const TString &jetAlgo, double minPt) {
  std::vector<TH1*> reso;
  std::vector<TF1*> fits;
  std::vector<TString> labels;
  
  // Nominal vs dijet selection
  if( jetAlgo == "PF" || jetAlgo == "JPT" ) {
    reso.push_back(util::FileOps::readTH1(file,"hSigmaGauss_Eta0_DeltaR10"));
    fits.push_back(util::FileOps::readTF1(file,"fit_Eta0_DeltaR10"));
  } else if( jetAlgo == "Calo" ) {
    reso.push_back(util::FileOps::readTH1(file,"hSigmaGauss_Eta0_DeltaR25"));
    fits.push_back(util::FileOps::readTF1(file,"fit_Eta0_DeltaR25"));
  }    
  reso.push_back(util::FileOps::readTH1(file,"hSigmaGauss_Eta0_Dijets"));
  fits.push_back(util::FileOps::readTF1(file,"fit_Eta0_Dijets"));
  labels.push_back("All events");
  labels.push_back("Dijet events");

  plots("DijetSelection",reso,fits,labels,jetAlgo,minPt);
  
  reso.clear();
  fits.clear();
  labels.clear();


  // Different DeltaR
  reso.push_back(util::FileOps::readTH1(file,"hSigmaGauss_Eta0_DeltaR10"));
  reso.push_back(util::FileOps::readTH1(file,"hSigmaGauss_Eta0_DeltaR20"));
  reso.push_back(util::FileOps::readTH1(file,"hSigmaGauss_Eta0_DeltaR25"));
  fits.push_back(util::FileOps::readTF1(file,"fit_Eta0_DeltaR10"));
  fits.push_back(util::FileOps::readTF1(file,"fit_Eta0_DeltaR20"));
  fits.push_back(util::FileOps::readTF1(file,"fit_Eta0_DeltaR25"));
  labels.push_back("#DeltaR_{max} < 0.10    ");
  labels.push_back("#DeltaR_{max} < 0.20    ");
  labels.push_back("#DeltaR_{max} < 0.25    ");

  plots("DeltaRMax",reso,fits,labels,jetAlgo,minPt);
  
  reso.clear();
  fits.clear();
  labels.clear();
}


// Compare MCTruth distributions for dijet vs all selection
void plotMCTruthDistributionsForSelection(const TString &fileNameAll, const TString &fileNameDijets) {
  util::StyleSettings::paperNoTitle();

  util::HistVec hRespAll = util::FileOps::readHistVec(fileNameAll,"hMCTruthResp_Eta0_DeltaR10","hRespAll");
  util::HistVec hRespDijets = util::FileOps::readHistVec(fileNameDijets,"hMCTruthResp_Eta0_Dijets","hRespDijets");
  assert( hRespAll.size() == hRespDijets.size() );

  for(unsigned int i = 0; i < hRespAll.size(); ++i) {
    hRespAll[i]->UseCurrentStyle();
    if( hRespAll[i]->Integral("width") ) hRespAll[i]->Scale(1./hRespAll[i]->Integral("width"));
    hRespAll[i]->SetMarkerStyle(1);
    hRespAll[i]->SetLineColor(kBlack);

    hRespDijets[i]->UseCurrentStyle();
    if( hRespDijets[i]->Integral("width") ) hRespDijets[i]->Scale(1./hRespDijets[i]->Integral("width"));
    hRespDijets[i]->SetMarkerStyle(1);
    hRespDijets[i]->SetLineColor(kRed);
  }

  TCanvas* can = new TCanvas("can","MCTruthResponse",500,500);
  for(unsigned int i = 0; i < hRespAll.size(); ++i) {
    can->cd();
    hRespAll[i]->Draw("HIST");
    hRespDijets[i]->Draw("HISTsame");
    can->SetLogy();
    can->SaveAs("MCTruthResponseDistribution_DijetSelection_PtBin"+util::toTString(i)+".eps","eps");
  }
  delete can;
}


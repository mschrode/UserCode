// $Id: proofOfPrinciple.C,v 1.1 2012/01/05 11:35:41 mschrode Exp $

#include <cassert>
#include <iostream>
#include <vector>

#include "TArrow.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TString.h"

#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


const double sigma_ = 0.1;
const double gaussWidthNSig_ = 2.5;
const double gaussStart_ = 1. - gaussWidthNSig_*sigma_;
const double gaussEnd_ = 1. + gaussWidthNSig_*sigma_;
const double tailStart_ = 0.3;
const double gaussStartAsym_ = -gaussWidthNSig_*sigma_/sqrt(2.);
const double gaussEndAsym_ = -gaussStartAsym_;

TRandom* rand_ = new TRandom3(0);

double drawResponse(double tailFrac) {
  double res = 0.;
  if( rand_->Uniform() < tailFrac ) {
    res = rand_->Uniform(tailStart_,gaussStart_);
  } else {
    res = rand_->Gaus(0.,sigma_);
    while( res < gaussStart_ || res > gaussEnd_ ) res = rand_->Gaus(1.,sigma_);
  }
  assert( res > 0. );

  return res;
}


void generateEvents(int id, double tailFrac, unsigned int nEntries, TH1* &hResp, TH1* &hAsym) {
  TString name = "hResp";
  name += id;
  hResp = new TH1D(name,";Response;Probability density",100,0.,2.);
  name = "hAsym";
  name += id;
  hAsym = new TH1D(name,";Asymmetry;Probability density",100,-1.,1.);

  for(unsigned int n = 0; n < nEntries; ++n) {
    double r1 = drawResponse(tailFrac);
    double r2 = drawResponse(tailFrac);
    double a = (r1 - r2) / (r1 + r2);
    
    hResp->Fill(r1);
    hResp->Fill(r2);
    hAsym->Fill(a);
    hAsym->Fill(-a);
  }
}


void getTails(const TH1* hResp, const TH1* hAsym, TH1* &hRespTail, double &fTailResp, TH1* &hAsymTail, double &fTailAsym) {

  TString name = hResp->GetName();
  hRespTail = static_cast<TH1*>(hResp->Clone(name+"Tail"));

  name = hAsym->GetName();
  hAsymTail = static_cast<TH1*>(hAsym->Clone(name+"Tail"));
  for(int bin = hRespTail->FindBin(gaussStart_)+1; bin <= hRespTail->GetNbinsX(); ++bin) {
    hRespTail->SetBinContent(bin,0.);
    hRespTail->SetBinError(bin,0.);
  }

//   double asymCoreStart = -gaussWidthNSig_*sigma_/sqrt(2.);
//   double asymCoreEnd = -asymCoreStart;
  for(int bin = hAsymTail->FindBin(gaussStartAsym_); bin <= hAsymTail->FindBin(gaussEndAsym_); ++bin) {
    hAsymTail->SetBinContent(bin,0.);
    hAsymTail->SetBinError(bin,0.);
  }

  fTailResp = hRespTail->Integral()/hResp->Integral();
  fTailAsym = hAsymTail->Integral()/hAsym->Integral();
}



void proofOfPrinciple(unsigned int nEntries = 500000) {
  util::StyleSettings::setStyleNoteNoTitle();

  int nSets = 9;
  int nVar = nSets/2;
  std::vector<TH1*> hResp(nSets);
  std::vector<TH1*> hRespTail(nSets);
  std::vector<TH1*> hAsym(nSets);
  std::vector<TH1*> hAsymTail(nSets);

  std::vector<double> fTailResp(nSets,0.);
  std::vector<double> fTailAsym(nSets,0.);

  std::vector<double> varTailResp(nSets,0.);
  std::vector<double> varTailAsym(nSets,0.);
  std::vector<double> relVarTailAsym(nSets,0.);

  for(int i = 0; i < nSets; ++i) {
    double scale = (1.+0.5/nVar*(i-nVar)); // Vary fraction from 50% to 150%
    double tailFrac = 0.1*scale;
    generateEvents(i,tailFrac,nEntries,hResp.at(i),hAsym.at(i));
    getTails(hResp.at(i),hAsym.at(i),hRespTail.at(i),fTailResp.at(i),hAsymTail.at(i),fTailAsym.at(i));
    
    std::cout << "f(Tail) = " << tailFrac << ": f(TailAsym) = " << fTailAsym.at(i) << std::endl;
  }

  std::cout << std::endl;
  for(int i = 0; i < nSets; ++i) {
    varTailResp.at(i) = fTailResp.at(i)/fTailResp.at(nVar);
    varTailAsym.at(i) = fTailAsym.at(i)/fTailAsym.at(nVar);
    std::cout << i << "/" << nVar << ": " << varTailResp.at(i) << " --> " << varTailAsym.at(i) << std::endl;
    relVarTailAsym.at(i) = varTailAsym.at(i)/varTailResp.at(i);
  }


  double yMinResp = 3E-2;
  double yMinAsym = 3E-2;
  for(unsigned int i = 0; i < hResp.size(); ++i) {
    util::HistOps::normHist(hResp.at(i),"width");
    util::HistOps::normHist(hAsym.at(i),"width");
    util::HistOps::setYRange(hResp.at(i),2,yMinResp);
    util::HistOps::setYRange(hAsym.at(i),2,yMinAsym);
    hResp.at(i)->SetLineWidth(2);
    hAsym.at(i)->SetLineWidth(2);
  }
  hResp.front()->SetLineColor(kRed);
  hResp.front()->SetLineStyle(2);
  hResp.back()->SetLineColor(kBlue);
  hResp.back()->SetLineStyle(4);

  hAsym.front()->SetLineColor(kRed);
  hAsym.front()->SetLineStyle(2);
  hAsym.back()->SetLineColor(kBlue);
  hAsym.back()->SetLineStyle(4);

  TLine* linTailResp = new TLine(gaussStart_,yMinResp,gaussStart_,4.);
  linTailResp->SetLineWidth(2);
  linTailResp->SetLineColor(28);
  //  linTailResp->SetLineStyle(2);
  TArrow* arrResp = new TArrow(gaussStart_,1.5,tailStart_+0.09,1.5);
  arrResp->SetLineWidth(linTailResp->GetLineWidth());
  arrResp->SetLineColor(linTailResp->GetLineColor());
  arrResp->SetArrowSize(0.05);
  arrResp->SetAngle(30);
  TLatex* arrRespLabel = new TLatex(tailStart_+0.19,1.85,"Tail");
  arrRespLabel->SetTextFont(42);
  arrRespLabel->SetTextColor(arrResp->GetLineColor());
  TArrow* arrRespCore = new TArrow(gaussStart_,0.06,1.,0.06,0.05,"<-|");
  arrRespCore->SetLineWidth(linTailResp->GetLineWidth());
  arrRespCore->SetLineColor(linTailResp->GetLineColor());
  arrRespCore->SetAngle(40);
  TLatex* arrRespCoreLabel = new TLatex(gaussStart_+0.02,0.09,"2.5#sigma_{toy}");
  arrRespCoreLabel->SetTextFont(42);
  arrRespCoreLabel->SetTextColor(arrResp->GetLineColor());


  TLine* linTailAsymLeft = static_cast<TLine*>(linTailResp->Clone());
  linTailAsymLeft->SetX1(gaussStartAsym_);
  linTailAsymLeft->SetY1(yMinAsym);
  linTailAsymLeft->SetX2(gaussStartAsym_);
  linTailAsymLeft->SetY2(4.);
  TLine* linTailAsymRight = static_cast<TLine*>(linTailAsymLeft->Clone());
  linTailAsymRight->SetX1(gaussEndAsym_);
  linTailAsymRight->SetX2(gaussEndAsym_);
  TArrow* arrAsymLeft = static_cast<TArrow*>(arrResp->Clone());
  arrAsymLeft->SetX1(gaussStartAsym_);
  arrAsymLeft->SetY1(1.5);
  arrAsymLeft->SetX2(-0.6);
  arrAsymLeft->SetY2(arrAsymLeft->GetY1());
  TArrow* arrAsymRight = static_cast<TArrow*>(arrAsymLeft->Clone());
  arrAsymRight->SetX1(gaussEndAsym_);
  arrAsymRight->SetX2(0.6);
  TLatex* arrAsymLeftLabel = static_cast<TLatex*>(arrRespLabel->Clone());
  arrAsymLeftLabel->SetX(-0.56);
  arrAsymLeftLabel->SetY(1.85);
  TArrow* arrAsymCore = static_cast<TArrow*>(arrRespCore->Clone());
  arrAsymCore->SetX1(gaussStartAsym_);
  arrAsymCore->SetX2(0.);
  TLatex* arrAsymCoreLabel = new TLatex(gaussStartAsym_+0.01,0.12,"#frac{2.5}{#sqrt{2}}#sigma_{toy}");
  arrAsymCoreLabel->SetTextFont(42);
  arrAsymCoreLabel->SetTextColor(arrResp->GetLineColor());

  TPaveText* label = util::LabelFactory::createPaveText(1,-0.5);
  label->AddText("Toy simulation");

  TLegend* leg = util::LabelFactory::createLegendCol(5,0.3);
  leg->AddEntry(hResp.front(),"f^{toy}_{resp} #times "+util::toTString(varTailResp.front(),1),"L");
  leg->AddEntry(hResp.at(nVar),"f^{toy}_{resp} #times 1.0","L");
  leg->AddEntry(hResp.back(),"f^{toy}_{resp} #times "+util::toTString(varTailResp.back(),1),"L");

  TCanvas* cResp = new TCanvas("cResp","Response",500,500);
  cResp->cd();
  hResp.front()->Draw("L");
  hResp.back()->Draw("Lsame");
  hResp.at(nVar)->Draw("Lsame");
  leg->Draw("same");
  label->Draw("same");
  linTailResp->Draw("same");
  arrRespLabel->Draw("same");
  arrResp->Draw();
  arrRespCore->Draw();
  arrRespCoreLabel->Draw("same");
  cResp->SetLogy();
  cResp->SaveAs("TailStudy_ProofOfPrinciples_Response.eps","eps");

  TCanvas* cAsym = new TCanvas("cAsym","Asymmetry",500,500);
  cAsym->cd();
  hAsym.front()->Draw("L");
  hAsym.back()->Draw("Lsame");
  hAsym.at(nVar)->Draw("Lsame");
  leg->Draw("same");
  label->Draw("same");
  linTailAsymLeft->Draw("same");
  arrAsymLeftLabel->Draw("same");
  arrAsymLeft->Draw();
  linTailAsymRight->Draw("same");
  arrAsymRight->Draw();
  arrAsymCore->Draw();
  arrAsymCoreLabel->Draw("same");
  cAsym->SetLogy();
  cAsym->SaveAs("TailStudy_ProofOfPrinciples_Asymmetry.eps","eps");

  TGraph* gFAsymVsFResp = new TGraph(varTailResp.size(),&(varTailResp.front()),&(varTailAsym.front()));
  gFAsymVsFResp->SetMarkerStyle(20);
  gFAsymVsFResp->SetMarkerSize(1.5);

  TCanvas* cFvsF = new TCanvas("cFvsF","FAsym vs FResp",500,500);
  cFvsF->cd();
  TH1* hFrameFAsymVsFResp = new TH1D("hFrameFAsymVsFResp",";Scale of f^{toy}_{resp};Scale of f^{toy}_{asym}",1000,0.4,1.6);
  for(int bin = 1; bin <= hFrameFAsymVsFResp->GetNbinsX(); ++bin) {
    hFrameFAsymVsFResp->SetBinContent(bin,hFrameFAsymVsFResp->GetXaxis()->GetBinCenter(bin));
  }
  hFrameFAsymVsFResp->SetLineStyle(2);
  hFrameFAsymVsFResp->GetYaxis()->SetRangeUser(0.4,1.6);
  hFrameFAsymVsFResp->Draw();
  gFAsymVsFResp->Draw("Psame");
  label->Draw("same");
  cFvsF->SaveAs("TailStudy_ProofOfPrinciples_Variation.eps","eps");

  TGraph* gRelScale = new TGraph(varTailResp.size(),&(varTailResp.front()),&(relVarTailAsym.front()));
  gRelScale->SetMarkerStyle(20);
  gRelScale->SetMarkerSize(1.5);

  TCanvas* cRelVar = new TCanvas("cRelVar","Relative variation",500,500);
  cRelVar->cd();
  TH1* hFrameRelScale = new TH1D("hFrameRelScale",";Scale of f^{toy}_{resp};Scale of f^{toy}_{asym} / Scale of f^{toy}_{resp}",1000,0.4,1.6);
  for(int bin = 1; bin <= hFrameRelScale->GetNbinsX(); ++bin) {
    hFrameRelScale->SetBinContent(bin,1.);
  }
  hFrameRelScale->SetLineStyle(2);
  hFrameRelScale->GetYaxis()->SetRangeUser(0.55,1.45);
  hFrameRelScale->Draw();
  gRelScale->Draw("Psame");
  label->Draw("same");
  cRelVar->SaveAs("TailStudy_ProofOfPrinciples_RelVariation.eps","eps");
}



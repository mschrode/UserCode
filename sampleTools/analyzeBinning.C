// $Id: $
//
// Plot ptAve spectra from Kalibri ntuples for
//  - all selected events (inclusive trigger paths)
//  - events selected from >99% efficient trigger paths
// Determine ptAve bins such that
//  1) primary bin edges correspond to the trigger thresholds; and
//  2) in between, bins are further divided if the number of events
//     remaining in the new bins is greater than a certain limit

#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TString.h"
#include "TVector2.h"

#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/FileOps.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"



unsigned int COUNT = 0;



void readPtAveSpectra(const TString &fileName, const std::vector<double> &etaBinEdges, double minDeltaPhi, double maxPt3, int nEvts, const std::vector<double> &trigThres, std::vector<TH1*> &hPtAveIncl, std::vector<TH1*> &hPtAveFullEff) {
  
  unsigned int nEtaBins = etaBinEdges.size()-1;
  unsigned int nEvtsSel = 0;

  std::cout << "Creating histograms...\n";
  for(size_t etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    TH1* h = util::HistOps::createTH1D("hPtAveIncl_EtaBin"+util::toTString(etaBin),500,0.,1500.,"p^{ave}_{T}","GeV","events");
    h->Sumw2();
    h->SetMarkerStyle(20);
    hPtAveIncl.push_back(h);

    hPtAveFullEff.push_back(static_cast<TH1*>(h->Clone("hPtAveFullEff_EtaBin"+util::toTString(etaBin))));
  }

  std::cout << "Adding files to chain...\n";
  // Set up chain
  TChain* chain = util::FileOps::createTChain(fileName);

  // Set branch addresses
  const int maxNJet = 50;

  // For jet ordering
  util::JetIndexCol corrJets;

  int nObjJet = 0;
  float jetPt[maxNJet];
  float jetCorrL1[maxNJet];
  float jetCorrL2L3[maxNJet];
  float jetEta[maxNJet];
  float jetPhi[maxNJet];
  bool jetId[maxNJet];

  bool hlt30 = false;
  bool hlt60 = false;
  bool hlt80 = false;
  bool hlt110 = false;
  bool hlt150 = false;
  bool hlt190 = false;
  bool hlt240 = false;
  bool hlt300 = false;
  bool hlt370 = false;

  chain->SetBranchAddress("NobjJet",&nObjJet);
  chain->SetBranchAddress("JetPt",jetPt);
  chain->SetBranchAddress("JetCorrL1",jetCorrL1);
  chain->SetBranchAddress("JetCorrL2L3",jetCorrL2L3);
  chain->SetBranchAddress("JetEta",jetEta);
  chain->SetBranchAddress("JetPhi",jetPhi);
  chain->SetBranchAddress("JetIDLoose",jetId);

  chain->SetBranchAddress("HltDiJetAve30",&hlt30);
  chain->SetBranchAddress("HltDiJetAve60",&hlt60);
  chain->SetBranchAddress("HltDiJetAve80",&hlt80);
  chain->SetBranchAddress("HltDiJetAve110",&hlt110);
  chain->SetBranchAddress("HltDiJetAve150",&hlt150);
  chain->SetBranchAddress("HltDiJetAve190",&hlt190);
  chain->SetBranchAddress("HltDiJetAve240",&hlt240);
  chain->SetBranchAddress("HltDiJetAve300",&hlt300);
  chain->SetBranchAddress("HltDiJetAve370",&hlt370);


  // Loop over tree entries and fill histograms
  if( nEvts < 0 || nEvts > chain->GetEntries() ) nEvts = chain->GetEntries();
  std::cout << "Reading " << nEvts << " entries...\n";
  for(int n = 0; n < nEvts; ++n) {
    if( n%50000 == 0 ) std::cout << "  " << n << std::endl;
    chain->GetEntry(n);

    if( nObjJet > maxNJet ) {
      std::cerr << "WARNING: nObjJet = " << nObjJet << " > " << maxNJet << ". Skipping event.\n";
      continue;
    }

    if( nObjJet > 1 ) {
      
      // Order L2L3 corrected jets
      corrJets.clear();
      for(int j = 0; j < nObjJet; ++j) {
	corrJets.add(j,jetCorrL1[j]*jetCorrL2L3[j]*jetPt[j]);
      }
      corrJets.sort();
      
      if( std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets(0)]-jetPhi[corrJets(1)])) < minDeltaPhi ) continue;
      else if( !jetId[corrJets(0)] || !jetId[corrJets(1)] ) continue;
      double ptAve = 0.5*(corrJets.pt(0)+corrJets.pt(1));
      if( nObjJet > 2 && corrJets.pt(2) > maxPt3*ptAve ) continue;
	       
      // Find eta  bin
      unsigned int etaBin0 = 1000;
      unsigned int etaBin1 = 1000;
      if( util::findBin(std::abs(jetEta[corrJets(0)]),etaBinEdges,etaBin0) && util::findBin(std::abs(jetEta[corrJets(1)]),etaBinEdges,etaBin1) ) {
	if( etaBin0 == etaBin1 ) {
	  ++nEvtsSel;
	  hPtAveIncl.at(etaBin0)->Fill(ptAve);

	  unsigned int trigBin = 1000;
	  if( util::findBin(ptAve,trigThres,trigBin) ) {
	    bool triggered = false;
	    if( trigBin == 0 && hlt30 )
	      triggered = true;
	    else if( trigBin == 1 && ( hlt30 || hlt60 ) )
	      triggered = true;
	    else if( trigBin == 2 && ( hlt30 || hlt60 || hlt80 ) )
	      triggered = true;
	    else if( trigBin == 3 && ( hlt30 || hlt60 || hlt80 || hlt110 ) )
	      triggered = true;
	    else if( trigBin == 4 && ( hlt30 || hlt60 || hlt80 || hlt110 || hlt150 ) )
	      triggered = true;
	    else if( trigBin == 5 && ( hlt30 || hlt60 || hlt80 || hlt110 || hlt150 || hlt190 ) )
	      triggered = true;
	    else if( trigBin == 6 && ( hlt30 || hlt60 || hlt80 || hlt110 || hlt150 || hlt190 || hlt240 ) )
	      triggered = true;
	    else if( trigBin == 7 && ( hlt30 || hlt60 || hlt80 || hlt110 || hlt150 || hlt190 || hlt240 || hlt300 ) )
	      triggered = true;
	    else if( trigBin == 8 && ( hlt30 || hlt60 || hlt80 || hlt110 || hlt150 || hlt190 || hlt240 || hlt300 || hlt370 ) )
	      triggered = true;

	    if( triggered ) hPtAveFullEff.at(etaBin0)->Fill(ptAve);
	  }
	}
      }
    } // End if nObjJet > 1
  } // End loop over entries

  std::cout << "\n***** " << nEvtsSel << " dijet events selected *****\n" << std::endl;
}


void findBins(TH1* h, double low, double up, int &binMin, int &binMax) {
  binMin = h->FindBin(low);
  binMax = std::min(h->FindBin(up),h->FindBin(up-0.001));
}


double getNEvts(TH1* h, double low, double up) {
  int binMin = 0;
  int binMax = 0;
  findBins(h,low,up,binMin,binMax);
  
  return h->Integral(binMin,binMax);
}



std::vector<double> splitBin(TH1* h, double min, double max, double minBinContent) {
  ++COUNT;
  // The bin edges
  std::vector<double> edges;

  // Find corresponding bins and bin content in this range
  int binMin = 0;
  int binMax = 0;
  findBins(h,min,max,binMin,binMax);
  double content = h->Integral(binMin,binMax);
  if( content > 2.*minBinContent && binMax > binMin ) {
    // Copy bin content of selected range into new histogram
    TH1* hTmp = static_cast<TH1*>(h->Clone("hTmp_"+util::toTString(binMin)));
    hTmp->Reset();
    for(int b = binMin; b <= binMax; ++b) {
      hTmp->SetBinContent(b,h->GetBinContent(b));
      hTmp->SetBinError(b,h->GetBinError(b));
    }
    double mean = hTmp->GetBinLowEdge(hTmp->FindBin(hTmp->GetMean()));
    std::vector<double> edgesLeft = splitBin(hTmp,min,mean,minBinContent);
    std::vector<double> edgesRight = splitBin(hTmp,mean,max,minBinContent);
    delete hTmp;
    // Concat bin edges from left and right bin and take
    // into account double counted edge at 'mean'
    edges = edgesLeft;
    edges.pop_back();
    edges.insert(edges.end(),edgesRight.begin(),edgesRight.end());
  } else {
    edges.push_back(min);
    edges.push_back(max);
  }

  return edges;
}



std::vector<double> findNumEvtsPerTrigBin(TH1* hPtAve, const std::vector<double> &trigThres, int nMinEvts, TH1* &hNumTrigThres, TH1* &hNum) {
  TH1* hPtAveTmp = static_cast<TH1*>(hPtAve->Clone("hPtAveTmp"));
  TString name = hPtAve->GetName();

  // Number of events in each trigger threshold bin
  hNumTrigThres = static_cast<TH1*>(hPtAveTmp->Clone(name+"_NumTrigThres"));
  hNumTrigThres->Reset();
  std::cout << "\n\n  BINS FROM TRIGGER THRESHOLDS" << std::endl;
  for(unsigned int i = 1; i < trigThres.size(); ++i) {
    int binMin = 0;
    int binMax = 0;
    findBins(hPtAveTmp,trigThres.at(i-1),trigThres.at(i),binMin,binMax);
    double nEvts = hPtAveTmp->Integral(binMin,binMax);
    for(int bin = binMin; bin <= binMax; ++bin) {
      hNumTrigThres->SetBinContent(bin,nEvts);
    }
    std::cout << "    " << i << ": " << trigThres.at(i-1) << " - " << trigThres.at(i) << std::endl;
  }
  
  // Number of events in each finer bin
  std::vector<double> fineEdges;
  hNum = static_cast<TH1*>(hNumTrigThres->Clone(name+"_NumThresRebinned"));
  hNum->Reset();
  std::cout << "\n  FINE BINS" << std::endl;
  for(unsigned int i = 1; i < trigThres.size(); ++i) {
    COUNT = 0;
    std::vector<double> edges = splitBin(hPtAveTmp,trigThres.at(i-1),trigThres.at(i),nMinEvts);
    for(unsigned int j = 1; j < edges.size(); ++j) {
      std::cout << "    " << j << ": " << edges.at(j-1) << " - " << edges.at(j) << " --> " << std::flush;
      int binMin = 0;
      int binMax = 0;
      findBins(hPtAveTmp,edges.at(j-1),edges.at(j),binMin,binMax);
      double nEvts = hPtAveTmp->Integral(binMin,binMax);
      std::cout << nEvts << std::flush;
      for(int bin = binMin; bin <= binMax; ++bin) {
	hNum->SetBinContent(bin,nEvts);
      }
      if( j == 1 || j == edges.size()-1 ) std::cout << " *" << std::flush;
      std::cout << std::endl;
    }
    if( fineEdges.size() == 0 ) {
      fineEdges = edges;
    } else {
      fineEdges.pop_back();
      fineEdges.insert(fineEdges.end(),edges.begin(),edges.end());
    }
  }

  delete hPtAveTmp;

  return fineEdges;
}



void plotPtAveSpectra(std::vector<TH1*> &hPtAveIncl, std::vector<TH1*> &hPtAveFullEff, const std::vector<double> &etaBinEdges, const std::vector<double> &trigThres, double minBinContent, const TString &jetLabel, const TString &lumiLabel) {

  TString outNamePrefix = "ResolutionBinning_2011_AK5";
  if( jetLabel.Contains("Calo") ) outNamePrefix += "Calo";
  else if( jetLabel.Contains("PF") ) outNamePrefix += "PF";

  std::vector< std::vector<double> > binEdges;
  for(unsigned int etaBin = 0; etaBin < etaBinEdges.size()-1; ++etaBin) {
    TString outName = outNamePrefix + "_Eta" + util::toTString(etaBin);
    std::cout << "\n\n\nPROCESSING ETA BIN " << etaBin << ": " << etaBinEdges.at(etaBin) << " - " << etaBinEdges.at(etaBin+1) << std::endl;
    TH1* hIncl = hPtAveIncl.at(etaBin);
    hIncl->GetXaxis()->SetMoreLogLabels();
    hIncl->GetXaxis()->SetNoExponent();
//     hIncl->SetMarkerSize(1.5);
//     hIncl->SetLineWidth(3);
    util::HistOps::setYRange(hIncl,5,3E-1);

    TH1* hFullEff = hPtAveFullEff.at(etaBin);
    hFullEff->GetXaxis()->SetMoreLogLabels();
    hFullEff->GetXaxis()->SetNoExponent();
    util::HistOps::setYRange(hFullEff,5,3E-1);

    std::vector<TLine*> lines;
    for(unsigned int i = 0; i < trigThres.size(); ++i) {
      TLine* line = new TLine(trigThres.at(i),0.,trigThres.at(i),hIncl->GetBinContent(hIncl->GetMaximumBin()));
      line->SetLineColor(kRed);
      line->SetLineStyle(2);
      line->SetLineWidth(hIncl->GetLineWidth());
      lines.push_back(line);
    }
    
    TPaveText* lab = util::LabelFactory::createPaveText(1);
    lab->AddText(jetLabel+",  "+lumiLabel);
    
    TLegend* legIncl = util::LabelFactory::createLegendWithOffset(2,1);
    legIncl->AddEntry(hIncl,"Data: incl. trigger selection","P");
    legIncl->AddEntry(lines.front(),"Trigger thresholds","L");

    TLegend* legFullEff = util::LabelFactory::createLegendWithOffset(2,1);
    legFullEff->AddEntry(hFullEff,"Data: >99% efficient triggers","P");
    legFullEff->AddEntry(lines.front(),"Trigger thresholds","L");


    TCanvas* canIncl = new TCanvas("canPtAveInclEta"+util::toTString(etaBin),"PtAve Incl Eta"+util::toTString(etaBin),500,500);
    canIncl->cd();
    hIncl->GetXaxis()->SetRange(hIncl->FindBin(trigThres.at(0))-1,hIncl->GetNbinsX());
    hIncl->Draw("PE");
    for(unsigned int i = 0; i < trigThres.size(); ++i) {
      lines.at(i)->Draw("same");
    }
    lab->Draw("same");
    legIncl->Draw("same");
    canIncl->SetLogx();
    canIncl->SetLogy();
    canIncl->SaveAs(outName+"_PtAveInclusive.eps","eps");


    TCanvas* canFullEff = new TCanvas("canPtAveFullEffEta"+util::toTString(etaBin),"PtAve FullEff Eta"+util::toTString(etaBin),500,500);
    canFullEff->cd();
    hFullEff->GetXaxis()->SetRange(hFullEff->FindBin(trigThres.at(0))-1,hFullEff->GetNbinsX());
    hFullEff->Draw("PE");
    for(unsigned int i = 0; i < trigThres.size(); ++i) {
      lines.at(i)->Draw("same");
    }
    lab->Draw("same");
    legFullEff->Draw("same");
    canFullEff->SetLogx();
    canFullEff->SetLogy();
    canFullEff->SaveAs(outName+"_PtAveFullyEfficientTriggers.eps","eps");
    
    TH1* hNumTrigThres = 0;
    TH1* hNum = 0;
    std::vector<double> fineEdges = findNumEvtsPerTrigBin(hFullEff,trigThres,minBinContent,hNumTrigThres,hNum);
    binEdges.push_back(fineEdges);
    hNumTrigThres->GetYaxis()->SetTitle("Events");
    hNum->GetYaxis()->SetTitle("Events");
    util::HistOps::setYRange(hNumTrigThres,1);
    util::HistOps::setYRange(hNum,1);

    TCanvas* canNumTrigThres = new TCanvas("canNumTrigThres"+util::toTString(etaBin),"Num Evts Eta"+util::toTString(etaBin),500,500);
    canNumTrigThres->cd();
    hNumTrigThres->GetXaxis()->SetRange(hNumTrigThres->FindBin(trigThres.at(0))-1,hNumTrigThres->GetNbinsX());
    hNumTrigThres->GetXaxis()->SetMoreLogLabels();
    hNumTrigThres->GetXaxis()->SetNoExponent();
    hNumTrigThres->Draw("HIST");
    lab->Draw("same");
    canNumTrigThres->SetLogx();
    canNumTrigThres->SaveAs(outName+"_NumEvtsInTrigBins.eps","eps");

    TCanvas* canNum = new TCanvas("canNum"+util::toTString(etaBin),"Num Evts Eta"+util::toTString(etaBin),500,500);
    canNum->cd();
    hNum->GetXaxis()->SetRange(hNum->FindBin(trigThres.at(0))-1,hNum->GetNbinsX());
    hNum->GetXaxis()->SetMoreLogLabels();
    hNum->GetXaxis()->SetNoExponent();
    hNum->Draw("HIST");
    lab->Draw("same");
    canNum->SetLogx();
    canNum->SaveAs(outName+"_NumEvtsFlattenedBins.eps","eps");
  }

  std::cout << std::endl << std::endl;
  std::cout << "\\begin{tabular}{ll}\n\\toprule" << std::endl;
  std::cout << "$|\\eta|$ bins & \\ptave bins (\\gevnospace) \\\\" << std::endl;
  std::cout << "\\midrule" << std::endl;
  for(unsigned int etaBin = 0; etaBin < etaBinEdges.size()-1; ++etaBin) {
    std::cout << "$" << etaBinEdges.at(etaBin) << " -- " << std::flush;
    std::cout << etaBinEdges.at(etaBin+1) << "$ & $" << std::flush;
    for(unsigned int i = 0; i < binEdges.at(etaBin).size(); ++i ) {
      std::cout << binEdges.at(etaBin).at(i) << "$  $" << std::flush;
    }
    std::cout << "$ \\\\" << std::endl;
  }
  std::cout << "\\bottomrule\n\\end{tabular}" << std::endl;

  std::cout << std::endl << std::endl;
  for(unsigned int etaBin = 0; etaBin < etaBinEdges.size()-1; ++etaBin) {
    std::cout << etaBinEdges.at(etaBin) << "  " << std::flush;
    std::cout << etaBinEdges.at(etaBin+1) << "\t" << std::flush;
    std::cout << binEdges.at(etaBin).size() << "  " << std::flush;
    for(unsigned int i = 0; i < binEdges.at(etaBin).size(); ++i ) {
      std::cout << binEdges.at(etaBin).at(i) << " " << std::flush;
    }
    std::cout << std::endl;
  }
}





void analyzeBinning(int nEvts = -1) {
  util::StyleSettings::presentationNoTitle();

  std::vector<double> etaBinEdges;
  etaBinEdges.push_back(0.0);
  etaBinEdges.push_back(0.5);
  etaBinEdges.push_back(1.1);
  etaBinEdges.push_back(1.7);
  etaBinEdges.push_back(2.3);
  etaBinEdges.push_back(5.0);

  std::vector<double> trigThres;
  trigThres.push_back(45.);
  trigThres.push_back(75.);
  trigThres.push_back(100.);
  trigThres.push_back(135.);
  trigThres.push_back(175.);
  trigThres.push_back(220.);
  trigThres.push_back(270.);
  trigThres.push_back(335.);
  trigThres.push_back(405.);
  trigThres.push_back(1500.);

  const TString fileName = "~/lustre/data/Jet_Run2011A-PromptReco-v2_163337-163757/Jet_Run2011A-PromptReco-v2_163337-163757_ak5Calo.root";

  std::vector<TH1*> hPtAveIncl;
  std::vector<TH1*> hPtAveFullEff;
  readPtAveSpectra(fileName,etaBinEdges,2.7,0.15,nEvts,trigThres,hPtAveIncl,hPtAveFullEff);
  plotPtAveSpectra(hPtAveIncl,hPtAveFullEff,etaBinEdges,trigThres,1300,util::LabelFactory::labelJetAlgo(fileName),"L = 62 / pb");
}



void test() {
  TF1* s1 = new TF1("s1","[0]*exp(x*[1])",0,10);
  TF1* s2 = static_cast<TF1*>(s1->Clone("s2"));
  TF1* s3 = new TF1("s1","[0]*exp(x*[1])+0.5*(TMath::Erf([4]*(x-[5]))+1)*[2]*exp(x*[3])",0,10);
  TF1* ef = new TF1("ef","0.5*[2]*(TMath::Erf([0]*(x-[1]))+1)",0,10);

  s1->SetLineStyle(2);
  s1->SetParameter(0,0.5);
  s1->SetParameter(1,-0.5);

  s2->SetParameter(0,1.);
  s2->SetParameter(1,s1->GetParameter(1));

  ef->SetLineColor(kRed);
  ef->SetParameter(0,1.85);
  ef->SetParameter(1,4.);

  s3->SetLineColor(kBlue);
  s3->SetParameter(0,s1->GetParameter(0));
  s3->SetParameter(1,s1->GetParameter(1));
  s3->SetParameter(2,s2->GetParameter(0));
  s3->SetParameter(3,s2->GetParameter(1));
  s3->SetParameter(4,ef->GetParameter(0));
  s3->SetParameter(5,ef->GetParameter(1));

  ef->SetParameter(2,s2->Eval(3.));



  TH1* hFrame = new TH1D("hFrame","",500,0,10);
  hFrame->GetYaxis()->SetRangeUser(0.,2.);

  TCanvas* can = new TCanvas("can","",500,500);
  can->cd();
  hFrame->Draw();
  s1->Draw("same");
  s2->Draw("same");
  ef->Draw("same");
  s3->Draw("same");
}

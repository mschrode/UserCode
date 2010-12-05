#include <vector>

#include "TCanvas.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TString.h"

#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


bool pf_ = false;


void fill(std::vector<double> &ptData, std::vector<double> &ptDataErr,
	  std::vector<double> &resData, std::vector<double> &resDataErr,
	  std::vector<double> &ptMC, std::vector<double> &ptMCErr,
	  std::vector<double> &resMC, std::vector<double> &resMCErr) {


  if( pf_ ) {

    // Data
    ptData.push_back(91.8021);
    ptDataErr.push_back(0.104491);
    resData.push_back(0.0954454);
    resDataErr.push_back(0.0025649);
    ptData.push_back(105.351);
    ptDataErr.push_back(0.185798);
    resData.push_back(0.0976386);
    resDataErr.push_back(0.00235333);
    ptData.push_back(128.456);
    ptDataErr.push_back(0.232014);
    resData.push_back(0.0962946);
    resDataErr.push_back(0.0020285);
    ptData.push_back(155.155);
    ptDataErr.push_back(0.195605);
    resData.push_back(0.0911605);
    resDataErr.push_back(0.00257226);
    ptData.push_back(178.654);
    ptDataErr.push_back(0.319159);
    resData.push_back(0.0863274);
    resDataErr.push_back(0.00263398);
    ptData.push_back(215.677);
    ptDataErr.push_back(0.457658);
    resData.push_back(0.0799562);
    resDataErr.push_back(0.00208855);
    ptData.push_back(266.112);
    ptDataErr.push_back(0.702986);
    resData.push_back(0.0761535);
    resDataErr.push_back(0.00310977);
    ptData.push_back(317.615);
    ptDataErr.push_back(0.966533);
    resData.push_back(0.0596831);
    resDataErr.push_back(0.0035009);
    ptData.push_back(366.062);
    ptDataErr.push_back(1.20303);
    resData.push_back(0.0687347);
    resDataErr.push_back(0.00535838);
    ptData.push_back(431.08);
    ptDataErr.push_back(2.77011);
    resData.push_back(0.0700223);
    resDataErr.push_back(0.0057261);
    ptData.push_back(584.213);
    ptDataErr.push_back(9.44555);
    resData.push_back(0.0559435);
    resDataErr.push_back(0.00431718);


    // MC
    ptMC.push_back(91.5326);
    ptMCErr.push_back(0.0732383);
    resMC.push_back(0.104627);
    resMCErr.push_back(0.000440823);
    ptMC.push_back(105.19);
    ptMCErr.push_back(0.105585);
    resMC.push_back(0.099957);
    resMCErr.push_back(0.000385801);
    ptMC.push_back(128.855);
    ptMCErr.push_back(0.124863);
    resMC.push_back(0.0906682);
    resMCErr.push_back(0.000446479);
    ptMC.push_back(155.579);
    ptMCErr.push_back(0.0956096);
    resMC.push_back(0.0820148);
    resMCErr.push_back(0.000770148);
    ptMC.push_back(179.319);
    ptMCErr.push_back(0.120627);
    resMC.push_back(0.0796547);
    resMCErr.push_back(0.000856976);
    ptMC.push_back(216.225);
    ptMCErr.push_back(0.154511);
    resMC.push_back(0.0733606);
    resMCErr.push_back(0.000952259);
    ptMC.push_back(266.655);
    ptMCErr.push_back(0.153905);
    resMC.push_back(0.0688778);
    resMCErr.push_back(0.00140933);
    ptMC.push_back(317.139);
    ptMCErr.push_back(0.156735);
    resMC.push_back(0.0636228);
    resMCErr.push_back(0.00194805);
    ptMC.push_back(367.087);
    ptMCErr.push_back(0.155769);
    resMC.push_back(0.0601754);
    resMCErr.push_back(0.00264833);
    ptMC.push_back(433.519);
    ptMCErr.push_back(0.21931);
    resMC.push_back(0.0590016);
    resMCErr.push_back(0.00269661);
    ptMC.push_back(584.954);
    ptMCErr.push_back(0.451223);
    resMC.push_back(0.0541064);
    resMCErr.push_back(0.00286609);

  } else {

    // Data
    ptData.push_back(89.2315);
    ptDataErr.push_back(0.0423283);
    resData.push_back(0.144381);
    resDataErr.push_back(0.00160808);
    ptData.push_back(102.899);
    ptDataErr.push_back(0.0819252);
    resData.push_back(0.134734);
    resDataErr.push_back(0.00139581);
    ptData.push_back(126.47);
    ptDataErr.push_back(0.110565);
    resData.push_back(0.121451);
    resDataErr.push_back(0.0010883);
    ptData.push_back(153.032);
    ptDataErr.push_back(0.0999986);
    resData.push_back(0.113812);
    resDataErr.push_back(0.00143938);
    ptData.push_back(177.159);
    ptDataErr.push_back(0.163714);
    resData.push_back(0.102011);
    resDataErr.push_back(0.00140342);
    ptData.push_back(214.015);
    ptDataErr.push_back(0.227726);
    resData.push_back(0.0952304);
    resDataErr.push_back(0.00111573);
    ptData.push_back(264.051);
    ptDataErr.push_back(0.34904);
    resData.push_back(0.0881422);
    resDataErr.push_back(0.00163901);
    ptData.push_back(314.356);
    ptDataErr.push_back(0.509634);
    resData.push_back(0.0821156);
    resDataErr.push_back(0.00227603);
    ptData.push_back(363.886);
    ptDataErr.push_back(0.666311);
    resData.push_back(0.0837059);
    resDataErr.push_back(0.00309157);
    ptData.push_back(431.569);
    ptDataErr.push_back(1.41495);
    resData.push_back(0.0705763);
    resDataErr.push_back(0.00294646);
    ptData.push_back(582.175);
    ptDataErr.push_back(4.90939);
    resData.push_back(0.061655);
    resDataErr.push_back(0.00296379);


    // MC
    ptMC.push_back(89.0548);
    ptMCErr.push_back(0.0460599);
    resMC.push_back(0.147168);
    resMCErr.push_back(0.000340937);
    ptMC.push_back(103.108);
    ptMCErr.push_back(0.0684027);
    resMC.push_back(0.130479);
    resMCErr.push_back(0.000290813);
    ptMC.push_back(126.759);
    ptMCErr.push_back(0.0869041);
    resMC.push_back(0.117156);
    resMCErr.push_back(0.000341095);
    ptMC.push_back(153.795);
    ptMCErr.push_back(0.0705495);
    resMC.push_back(0.103921);
    resMCErr.push_back(0.000605248);
    ptMC.push_back(177.5);
    ptMCErr.push_back(0.0871607);
    resMC.push_back(0.0984676);
    resMCErr.push_back(0.000648818);
    ptMC.push_back(214.431);
    ptMCErr.push_back(0.110263);
    resMC.push_back(0.0901789);
    resMCErr.push_back(0.000701328);
    ptMC.push_back(264.826);
    ptMCErr.push_back(0.107452);
    resMC.push_back(0.0835493);
    resMCErr.push_back(0.00107623);
    ptMC.push_back(315.131);
    ptMCErr.push_back(0.108093);
    resMC.push_back(0.0782712);
    resMCErr.push_back(0.00154342);
    ptMC.push_back(365.499);
    ptMCErr.push_back(0.109099);
    resMC.push_back(0.0708297);
    resMCErr.push_back(0.00203336);
    ptMC.push_back(431.78);
    ptMCErr.push_back(0.15682);
    resMC.push_back(0.0683832);
    resMCErr.push_back(0.00212362);
    ptMC.push_back(582.871);
    ptMCErr.push_back(0.35049);
    resMC.push_back(0.0616795);
    resMCErr.push_back(0.00250518);
  }
}


TF1* fitResolution(TGraphErrors* g, const TString &name, double min, double max) {
  TF1* fit = new TF1(name,"sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",min,max);
  fit->FixParameter(0,0.);
  fit->SetParameter(1,1.);
  fit->SetParameter(2,0.01);
  fit->SetLineWidth(1);
  fit->SetLineColor(kRed);
  g->Fit(fit,"BNR");

  return fit;
}


void plotDataMCComparison() {

  util::StyleSettings::presentationNoTitle();

  // Input
  double ptMin = 80.;
  double ptMax = 1000.;

  TString outNamePrefix = "ResFit_DataMCComparison_";
  if( pf_ ) outNamePrefix += "PF_";
  else outNamePrefix += "Calo_";
  outNamePrefix += "Eta0_";

  std::vector<double> ptData;
  std::vector<double> ptDataErr;
  std::vector<double> resData;
  std::vector<double> resDataErr;
  std::vector<double> ptMC;
  std::vector<double> ptMCErr;
  std::vector<double> resMC;
  std::vector<double> resMCErr;

  fill(ptData,ptDataErr,resData,resDataErr,ptMC,ptMCErr,resMC,resMCErr);



  // Resolution
  TGraphErrors* gData = new TGraphErrors(ptData.size(),&(ptData.front()),&(resData.front()),
					 &(ptDataErr.front()),&(resDataErr.front()));
  gData->SetMarkerStyle(25);
  gData->SetMarkerColor(2);
  gData->SetLineColor(2);

  TGraphErrors* gMC = new TGraphErrors(ptMC.size(),&(ptMC.front()),&(resMC.front()),
				       &(ptMCErr.front()),&(resMCErr.front()));
  gMC->SetMarkerStyle(21);

  // Fit resolution
  TF1* fitData = fitResolution(gData,"fitData",ptMin,ptMax);
  TF1* fitMC = fitResolution(gMC,"fitData",ptMin,ptMax);
  fitMC->SetLineStyle(2);
  fitMC->SetLineColor(kBlack);


  // Ratio
  TGraphErrors* gRatio = util::HistOps::createRatioGraph(gData,gMC);
  gRatio->SetMarkerColor(1);
  gRatio->SetLineColor(1);

  TF1* fitRatio = new TF1("fitRatio","sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2]) / sqrt([3]*[3]/x/x + [4]*[4]/x + [5]*[5])",ptMin,ptMax);
  for(int i = 0; i < 3; ++i) {
    fitRatio->SetParameter(i,fitData->GetParameter(i));
    fitRatio->SetParError(i,fitData->GetParError(i));
    fitRatio->SetParameter(3+i,fitMC->GetParameter(i));
    fitRatio->SetParError(3+i,fitMC->GetParError(i));
  }
  fitRatio->SetLineWidth(1);
  fitRatio->SetLineColor(kBlack);


  TF1* fitRatio1 = new TF1("fitRatio1","pol0",ptMin,ptMax);
  fitRatio1->SetLineWidth(1);
  fitRatio1->SetLineColor(2);
  fitRatio1->SetLineStyle(2);

  TF1* fitRatio2 = new TF1("fitRatio2","[0] - exp(x/[1])",ptMin,ptMax);
  fitRatio2->SetParameter(0,1.2);
  fitRatio2->SetParameter(1,-100.);
  fitRatio2->SetLineWidth(1);
  fitRatio2->SetLineColor(2);
  fitRatio2->SetLineStyle(1);  

  std::cout << "\nfitRatio1: " << std::endl;
  gRatio->Fit(fitRatio1,"R0");
  std::cout << "\n\nfitRatio2: " << std::endl;
  gRatio->Fit(fitRatio2,"R0");

  std::cout << "\nPar0 = " << fitRatio2->GetParameter(0) << " (" << fitRatio2->GetParameter(0)-fitRatio2->GetParError(0) << ", " << fitRatio2->GetParameter(0)+fitRatio2->GetParError(0) << ")" << std::endl;
  std::cout << "Par1 = " << fitRatio2->GetParameter(1) << " (" << fitRatio2->GetParameter(1)-fitRatio2->GetParError(1) << ", " << fitRatio2->GetParameter(1)+fitRatio2->GetParError(1) << ")" << std::endl;

  

  // Resolution plots
  TPaveText* label = util::LabelFactory::createPaveText(1);
  if( pf_ ) label->AddText("AK5 PF-Jets,  |#eta| < 1.1");
  else label->AddText("AK5 Calo-Jets,  |#eta| < 1.1");

  TLegend* legRes = util::LabelFactory::createLegendColWithOffset(4,-0.8,1);
  legRes->AddEntry(gData,"Measurement (Data)","P");
  legRes->AddEntry(fitData,"Interpolation (Data)","L");  
  legRes->AddEntry(gMC,"Measurement (MC)","P");  
  legRes->AddEntry(fitMC,"Interpolation (MC)","L");  

  TLegend* legRatio = util::LabelFactory::createLegendColWithOffset(3,-0.8,1);
  legRatio->AddEntry(fitRatio,"Ratio of interpolations","L");
  legRatio->AddEntry(fitRatio1,"Fit: pol0","L");
  legRatio->AddEntry(fitRatio2,"Fit: [0]-exp(x/[1])","L");



  TCanvas* canRes = new TCanvas("canRes","Resolution",500,500);
  canRes->cd();
  TH1* hFrame = util::HistOps::createTH1D("hFrame",1000,ptMin,ptMax,"p^{ref}_{T}","GeV","#sigma / p^{ref}_{T}");
  hFrame->GetYaxis()->SetRangeUser(1E-3,0.29);
  hFrame->GetXaxis()->SetMoreLogLabels();
  hFrame->Draw();
  gMC->Draw("PE1same");
  gData->Draw("PE1same");
  label->Draw("same");
  legRes->Draw("same");
  canRes->SetLogx();
  canRes->SaveAs(outNamePrefix+"Resolution.eps","eps");

  canRes->cd();
  hFrame->Draw();
  gMC->Draw("PE1same");
  gData->Draw("PE1same");
  fitData->Draw("same");
  fitMC->Draw("same");
  label->Draw("same");
  legRes->Draw("same");
  canRes->SetLogx();
  canRes->SaveAs(outNamePrefix+"ResolutionFit.eps","eps");


  TCanvas* canRatio = new TCanvas("canRatio","Ratio",500,500);
  canRatio->cd();
  TH1* hRatioFrame = util::HistOps::createRatioFrame(hFrame,"#sigma/p^{ref}_{T}(Data)  /  #sigma/p^{ref}_{T}(MC)",0.9,1.55);
  hRatioFrame->GetXaxis()->SetMoreLogLabels();
  hRatioFrame->Draw();
  fitRatio->Draw("Esame");
  fitRatio1->Draw("same");
  fitRatio2->Draw("same");
  gRatio->Draw("PE1same");
  label->Draw("same");
  legRatio->Draw("same");
  canRatio->SetLogx();
  canRatio->SaveAs(outNamePrefix+"Ratio.eps","eps");
}




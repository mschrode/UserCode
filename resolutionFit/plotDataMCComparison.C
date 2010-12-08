#include <vector>

#include "TCanvas.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TString.h"

#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


const bool pf_ = true;
const unsigned int etaBin = 0;
const unsigned int interpolationMode = 1;



void fill(std::vector<double> &ptData, std::vector<double> &ptDataErr,
	  std::vector<double> &resData, std::vector<double> &resDataErr,
	  std::vector<double> &ptMC, std::vector<double> &ptMCErr,
	  std::vector<double> &resMC, std::vector<double> &resMCErr) {


  if( pf_ ) {
    if( etaBin == 0 ) {
      ptData.push_back(91.6322);
      ptDataErr.push_back(0.0829401);
      resData.push_back(0.0968704);
      resDataErr.push_back(0.0032796);
      ptData.push_back(104.991);
      ptDataErr.push_back(0.154839);
      resData.push_back(0.101744);
      resDataErr.push_back(0.00306573);
      ptData.push_back(128.397);
      ptDataErr.push_back(0.202448);
      resData.push_back(0.0957896);
      resDataErr.push_back(0.00263518);
      ptData.push_back(155.084);
      ptDataErr.push_back(0.173374);
      resData.push_back(0.0912976);
      resDataErr.push_back(0.00339121);
      ptData.push_back(178.499);
      ptDataErr.push_back(0.292372);
      resData.push_back(0.088414);
      resDataErr.push_back(0.00357359);
      ptData.push_back(215.716);
      ptDataErr.push_back(0.421154);
      resData.push_back(0.0797679);
      resDataErr.push_back(0.00274093);
      ptData.push_back(265.782);
      ptDataErr.push_back(0.635974);
      resData.push_back(0.0774815);
      resDataErr.push_back(0.0041593);
      ptData.push_back(317.07);
      ptDataErr.push_back(0.887554);
      resData.push_back(0.0628395);
      resDataErr.push_back(0.00479996);
      ptData.push_back(364.25);
      ptDataErr.push_back(1.16013);
      resData.push_back(0.0776843);
      resDataErr.push_back(0.00773629);
      ptData.push_back(431.036);
      ptDataErr.push_back(2.61106);
      resData.push_back(0.0684853);
      resDataErr.push_back(0.00952402);
      ptData.push_back(584.61);
      ptDataErr.push_back(9.10935);
      resData.push_back(0.0575793);
      resDataErr.push_back(0.006751);
      
      
      ptMC.push_back(91.4504);
      ptMCErr.push_back(0.0610422);
      resMC.push_back(0.104515);
      resMCErr.push_back(0.000560147);
      ptMC.push_back(105.058);
      ptMCErr.push_back(0.0890344);
      resMC.push_back(0.100225);
      resMCErr.push_back(0.000497395);
      ptMC.push_back(128.722);
      ptMCErr.push_back(0.110365);
      resMC.push_back(0.0908618);
      resMCErr.push_back(0.000583306);
      ptMC.push_back(155.544);
      ptMCErr.push_back(0.0861829);
      resMC.push_back(0.08232);
      resMCErr.push_back(0.00103103);
      ptMC.push_back(179.163);
      ptMCErr.push_back(0.109961);
      resMC.push_back(0.0803683);
      resMCErr.push_back(0.0011492);
      ptMC.push_back(216.035);
      ptMCErr.push_back(0.142244);
      resMC.push_back(0.0738747);
      resMCErr.push_back(0.00128465);
      ptMC.push_back(266.414);
      ptMCErr.push_back(0.142446);
      resMC.push_back(0.0696178);
      resMCErr.push_back(0.00190508);
      ptMC.push_back(317.006);
      ptMCErr.push_back(0.146337);
      resMC.push_back(0.0639395);
      resMCErr.push_back(0.00264645);
      ptMC.push_back(366.921);
      ptMCErr.push_back(0.145815);
      resMC.push_back(0.0606798);
      resMCErr.push_back(0.00361633);
      ptMC.push_back(433.335);
      ptMCErr.push_back(0.207263);
      resMC.push_back(0.0592833);
      resMCErr.push_back(0.00367353);
      ptMC.push_back(584.588);
      ptMCErr.push_back(0.430757);
      resMC.push_back(0.0547319);
      resMCErr.push_back(0.0039373);
    }



    // Data
//     ptData.push_back(91.8021);
//     ptDataErr.push_back(0.104491);
//     resData.push_back(0.0954454);
//     resDataErr.push_back(0.0025649);
//     ptData.push_back(105.351);
//     ptDataErr.push_back(0.185798);
//     resData.push_back(0.0976386);
//     resDataErr.push_back(0.00235333);
//     ptData.push_back(128.456);
//     ptDataErr.push_back(0.232014);
//     resData.push_back(0.0962946);
//     resDataErr.push_back(0.0020285);
//     ptData.push_back(155.155);
//     ptDataErr.push_back(0.195605);
//     resData.push_back(0.0911605);
//     resDataErr.push_back(0.00257226);
//     ptData.push_back(178.654);
//     ptDataErr.push_back(0.319159);
//     resData.push_back(0.0863274);
//     resDataErr.push_back(0.00263398);
//     ptData.push_back(215.677);
//     ptDataErr.push_back(0.457658);
//     resData.push_back(0.0799562);
//     resDataErr.push_back(0.00208855);
//     ptData.push_back(266.112);
//     ptDataErr.push_back(0.702986);
//     resData.push_back(0.0761535);
//     resDataErr.push_back(0.00310977);
//     ptData.push_back(317.615);
//     ptDataErr.push_back(0.966533);
//     resData.push_back(0.0596831);
//     resDataErr.push_back(0.0035009);
//     ptData.push_back(366.062);
//     ptDataErr.push_back(1.20303);
//     resData.push_back(0.0687347);
//     resDataErr.push_back(0.00535838);
//     ptData.push_back(431.08);
//     ptDataErr.push_back(2.77011);
//     resData.push_back(0.0700223);
//     resDataErr.push_back(0.0057261);
//     ptData.push_back(584.213);
//     ptDataErr.push_back(9.44555);
//     resData.push_back(0.0559435);
//     resDataErr.push_back(0.00431718);


    // MC
//     ptMC.push_back(91.5326);
//     ptMCErr.push_back(0.0732383);
//     resMC.push_back(0.104627);
//     resMCErr.push_back(0.000440823);
//     ptMC.push_back(105.19);
//     ptMCErr.push_back(0.105585);
//     resMC.push_back(0.099957);
//     resMCErr.push_back(0.000385801);
//     ptMC.push_back(128.855);
//     ptMCErr.push_back(0.124863);
//     resMC.push_back(0.0906682);
//     resMCErr.push_back(0.000446479);
//     ptMC.push_back(155.579);
//     ptMCErr.push_back(0.0956096);
//     resMC.push_back(0.0820148);
//     resMCErr.push_back(0.000770148);
//     ptMC.push_back(179.319);
//     ptMCErr.push_back(0.120627);
//     resMC.push_back(0.0796547);
//     resMCErr.push_back(0.000856976);
//     ptMC.push_back(216.225);
//     ptMCErr.push_back(0.154511);
//     resMC.push_back(0.0733606);
//     resMCErr.push_back(0.000952259);
//     ptMC.push_back(266.655);
//     ptMCErr.push_back(0.153905);
//     resMC.push_back(0.0688778);
//     resMCErr.push_back(0.00140933);
//     ptMC.push_back(317.139);
//     ptMCErr.push_back(0.156735);
//     resMC.push_back(0.0636228);
//     resMCErr.push_back(0.00194805);
//     ptMC.push_back(367.087);
//     ptMCErr.push_back(0.155769);
//     resMC.push_back(0.0601754);
//     resMCErr.push_back(0.00264833);
//     ptMC.push_back(433.519);
//     ptMCErr.push_back(0.21931);
//     resMC.push_back(0.0590016);
//     resMCErr.push_back(0.00269661);
//     ptMC.push_back(584.954);
//     ptMCErr.push_back(0.451223);
//     resMC.push_back(0.0541064);
//     resMCErr.push_back(0.00286609);

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
  TF1* fit = 0;
  if( interpolationMode == 0 ) {
    fit = new TF1(name,"sqrt(sq([0]/x) + sq([1])/x + sq([2]))",min,max);
    fit->FixParameter(0,0.);
    fit->SetParameter(1,1.);
    fit->SetParameter(2,0.01);
  } else if( interpolationMode == 1 ) {
    fit = new TF1(name,"sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))",min,max);
    fit->SetParameter(0,-0.34);
    fit->SetParameter(1,0.3);
    fit->FixParameter(2,0.);
    fit->SetParameter(3,0.47);
  }
  fit->SetLineWidth(1);
  g->Fit(fit,"BNR");

  return fit;
}


TF1* fitAddConst(TGraphErrors* g, const TString &name, const TF1* mcFit) {
  TF1* fit = 0;
  if( interpolationMode == 0 ) {
    fit = new TF1(name,"sqrt(sq([0]/x) + sq([1])/x + sq([2]) + sq([3]))",mcFit->GetXmin(),mcFit->GetXmax());
    fit->FixParameter(0,mcFit->GetParameter(0));
    fit->FixParameter(1,mcFit->GetParameter(1));
    fit->FixParameter(2,mcFit->GetParameter(2));
    fit->SetParameter(3,mcFit->GetParameter(3));
  } else if( interpolationMode == 1 ) {
    fit = new TF1(name,"sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))",mcFit->GetXmin(),mcFit->GetXmax());
    fit->FixParameter(0,mcFit->GetParameter(0));
    fit->FixParameter(1,mcFit->GetParameter(1));
    fit->SetParameter(2,0.);
    fit->FixParameter(3,mcFit->GetParameter(3));
  }
  fit->SetLineWidth(1);
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
  TF1* fitMC = fitResolution(gMC,"fitData",ptMin,ptMax);
  fitMC->SetLineStyle(2);
  fitMC->SetLineColor(kBlack);

  TF1* fitData1 = fitResolution(gData,"fitData1",ptMin,ptMax);
  fitData1->SetLineColor(kRed);
  TF1* fitData2 = fitAddConst(gData,"fitData2",fitMC);
  fitData2->SetLineColor(kBlue);


  // Ratio
  TGraphErrors* gRatio = util::HistOps::createRatioGraph(gData,gMC);
  gRatio->SetMarkerColor(1);
  gRatio->SetLineColor(1);

  TF1* fitRatio1 = 0;
  if( interpolationMode == 0 ) {
    fitRatio1 = new TF1("fitRatio1","sqrt(sq([0]/x) + sq([1])/x + sq([2])) / sqrt(sq([3]/x) + sq([4])/x + sq([5]))",ptMin,ptMax);
  } else if( interpolationMode == 1 ) {
    fitRatio1 = new TF1("fitRatio1","sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2])) / sqrt(((TMath::Sign(1,[4])*sq([4]/x))+(sq([5])*(x^([7]-1))))+sq([6]))",ptMin,ptMax);
  }
  int nPar = fitData1->GetNpar();
  for(int i = 0; i < nPar; ++i) {
    fitRatio1->SetParameter(i,fitData1->GetParameter(i));
    fitRatio1->SetParError(i,fitData1->GetParError(i));
    fitRatio1->SetParameter(nPar+i,fitMC->GetParameter(i));
    fitRatio1->SetParError(nPar+i,fitMC->GetParError(i));
  }
  fitRatio1->SetLineWidth(1);
  fitRatio1->SetLineColor(fitData1->GetLineColor());

  TF1* fitRatio2 = 0;
  if( interpolationMode == 0 ) {
    fitRatio2 = new TF1("fitRatio2","sqrt(sq([0]/x) + sq([1])/x + sq([2]) + sq([3])) / sqrt(sq([4]/x) + sq([5])/x + sq([6]))",ptMin,ptMax);
    for(int i = 0; i < 3; ++i) {
      fitRatio2->SetParameter(i,fitMC->GetParameter(i));
      fitRatio2->SetParError(i,fitMC->GetParError(i));
      fitRatio2->SetParameter(4+i,fitMC->GetParameter(i));
      fitRatio2->SetParError(4+i,fitMC->GetParError(i));
    }
    fitRatio2->SetParameter(3,fitData2->GetParameter(3));
    fitRatio2->SetParError(3,fitData2->GetParError(3));
  } else if( interpolationMode == 1 ) {
    fitRatio2 = new TF1("fitRatio1","sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2])) / sqrt(((TMath::Sign(1,[4])*sq([4]/x))+(sq([5])*(x^([7]-1))))+sq([6]))",ptMin,ptMax);
    for(int i = 0; i < fitData2->GetNpar(); ++i) {
      fitRatio2->SetParameter(i,fitData2->GetParameter(i));
      fitRatio2->SetParError(i,fitData2->GetParError(i));
      fitRatio2->SetParameter(fitData2->GetNpar()+i,fitMC->GetParameter(i));
      fitRatio2->SetParError(fitData2->GetNpar()+i,fitMC->GetParError(i));
    }
  } 
  fitRatio2->SetLineWidth(1);
  fitRatio2->SetLineColor(fitData2->GetLineColor());
  
  


//   TF1* fitRatio1 = new TF1("fitRatio1","pol0",ptMin,ptMax);
//   fitRatio1->SetLineWidth(1);
//   fitRatio1->SetLineColor(2);
//   fitRatio1->SetLineStyle(2);

//   TF1* fitRatio2 = new TF1("fitRatio2","[0] - exp(x/[1])",ptMin,ptMax);
//   fitRatio2->SetParameter(0,1.2);
//   fitRatio2->SetParameter(1,-100.);
//   fitRatio2->SetLineWidth(1);
//   fitRatio2->SetLineColor(2);
//   fitRatio2->SetLineStyle(1);  

//   std::cout << "\nfitRatio1: " << std::endl;
//   gRatio->Fit(fitRatio1,"R0");
//   std::cout << "\n\nfitRatio2: " << std::endl;
//   gRatio->Fit(fitRatio2,"R0");

//   std::cout << "\nPar0 = " << fitRatio2->GetParameter(0) << " (" << fitRatio2->GetParameter(0)-fitRatio2->GetParError(0) << ", " << fitRatio2->GetParameter(0)+fitRatio2->GetParError(0) << ")" << std::endl;
//   std::cout << "Par1 = " << fitRatio2->GetParameter(1) << " (" << fitRatio2->GetParameter(1)-fitRatio2->GetParError(1) << ", " << fitRatio2->GetParameter(1)+fitRatio2->GetParError(1) << ")" << std::endl;

  

  // Resolution plots
  TPaveText* label = util::LabelFactory::createPaveText(1);
  if( pf_ ) label->AddText("AK5 PF-Jets,  |#eta| < 1.1");
  else label->AddText("AK5 Calo-Jets,  |#eta| < 1.1");

  TLegend* legRes = util::LabelFactory::createLegendColWithOffset(4,-0.8,1);
  legRes->AddEntry(gData,"Measurement (Data)","P");
  legRes->AddEntry(fitData1,"Interpolation (Data)","L");  
  legRes->AddEntry(gMC,"Measurement (MC)","P");  
  legRes->AddEntry(fitMC,"Interpolation (MC)","L");  

  TLegend* legRatio = util::LabelFactory::createLegendColWithOffset(3,-0.8,1);
  legRatio->AddEntry(fitRatio1,"Ratio of interpolations","L");
  legRatio->AddEntry(fitRatio2,"Fit: pol0","L");



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
  //  canRes->SaveAs(outNamePrefix+"Resolution.eps","eps");

  canRes->cd();
  hFrame->Draw();
  gMC->Draw("PE1same");
  gData->Draw("PE1same");
  fitMC->Draw("same");
  fitData1->Draw("same");
  fitData2->Draw("same");
  label->Draw("same");
  legRes->Draw("same");
  canRes->SetLogx();
  //  canRes->SaveAs(outNamePrefix+"ResolutionFit.eps","eps");


  TCanvas* canRatio = new TCanvas("canRatio","Ratio",500,500);
  canRatio->cd();
  TH1* hRatioFrame = util::HistOps::createRatioFrame(hFrame,"#sigma/p^{ref}_{T}(Data)  /  #sigma/p^{ref}_{T}(MC)",0.9,1.55);
  hRatioFrame->GetXaxis()->SetMoreLogLabels();
  hRatioFrame->Draw();
  fitRatio1->Draw("same");
  fitRatio2->Draw("same");
  gRatio->Draw("PE1same");
  label->Draw("same");
  legRatio->Draw("same");
  canRatio->SetLogx();
  //  canRatio->SaveAs(outNamePrefix+"Ratio.eps","eps");
}




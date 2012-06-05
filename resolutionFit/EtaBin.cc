// $Id: EtaBin.cc,v 1.16 2012/05/31 20:17:43 mschrode Exp $

#include <algorithm>
#include <iostream>

#include "../util/HistOps.h"

#include "EtaBin.h"


namespace resolutionFit{

  // -------------------------------------------------------------------------------------
  EtaBin::EtaBin(unsigned int etaBin, unsigned int nPtBins, const Parameters* par)
    : par_(par), etaBin_(etaBin) {
    for(unsigned int i = 0; i < nPtBins; ++i) {
      ptBins_.push_back(new PtBin(etaBin_,i,par));
    }
    mcTruthReso_ = new ResolutionFunctionModifiedNSC();
    scaledMCTruthReso_ = new ResolutionFunctionModifiedNSC();
    pli_ = new ResolutionFunctionModifiedNSC();

    kVal_ = 1.;
    kValStat_ = 0.;
    kValSystDown_ = 0.;
    kValSystUp_ = 0.;
  }


  // -------------------------------------------------------------------------------------
  EtaBin::~EtaBin() {
    for(std::vector<PtBin*>::iterator it = ptBins_.begin(); it != ptBins_.end(); ++it) {
      delete *it;
    }
    for(ComparedSamples::iterator it = compSamples_.begin(); it != compSamples_.end(); ++it) {
      delete *it;
    }
    delete pli_;
    delete mcTruthReso_;
    delete scaledMCTruthReso_;
  }


  // -------------------------------------------------------------------------------------
  TString EtaBin::toString() const {
    TString str = "Eta";
    str += etaBin(); 

    return str;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::hasDataSample() const {
    bool hasData = false;
    for(SampleTypeIt it = sampleTypesBegin(); it != sampleTypesEnd(); ++it) {
      if( it->second == Sample::Data ) {
	hasData = true;
	break;
      }
    }

    return hasData;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::addDataSample(const TString &label, const TString &fileName, const TString &printLabel) {
    bool result = true;
    for(std::vector<PtBin*>::iterator it = ptBins_.begin(); it != ptBins_.end(); ++it) {
      bool status = (*it)->addDataSample(label,fileName,printLabel);
      result = result && status;
    }
    sampleTypes_[label] = Sample::Data;

    return result;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::addMCSample(const TString &label, const TString &fileName, const TString &printLabel) {
    bool result = true;
    for(std::vector<PtBin*>::iterator it = ptBins_.begin(); it != ptBins_.end(); ++it) {
      bool status = (*it)->addMCSample(label,fileName,printLabel);
      result = result && status;
    }
    sampleTypes_[label] = Sample::MC;

    return result;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::addFitResult(FitResult::Type type) {
    bool result = true;

    // Add FitResult object to samples and perform extrapolation
    for(std::vector<PtBin*>::iterator it = ptBins_.begin(); it != ptBins_.end(); ++it) {
      bool status = (*it)->addFitResult(type);
      result = result && status;
    }

    // Perfom kSoft fits
    std::map<SampleLabel,TF1*> kSoftFits;
    for(SampleTypeIt sTIt = sampleTypesBegin(); sTIt != sampleTypesEnd(); ++sTIt) {
      TString name = "kSoftFit"+sTIt->first;
      kSoftFits[sTIt->first] = fitKSoftSlope(name,sTIt->first,type);
    }

    // Add kSoft fits to FitResult for each sample
    for(std::map<SampleLabel,TF1*>::const_iterator sIt = kSoftFits.begin();
	sIt != kSoftFits.end(); ++sIt) {
      for(std::vector<PtBin*>::iterator ptIt = ptBins_.begin(); ptIt != ptBins_.end(); ++ptIt) {
	bool status = (*ptIt)->setKSoftFit(sIt->first,type,sIt->second);
	result = result && status;
      }
    }
    for(std::map<SampleLabel,TF1*>::iterator it = kSoftFits.begin();
	it != kSoftFits.end(); ++it) {
      delete it->second;
    }
    
    fitResultTypes_.insert(type);

    return result;
  }


  // -------------------------------------------------------------------------------------
  TF1* EtaBin::fitKSoftSlope(const TString &name, const SampleLabel &label, FitResult::Type type) const {
    TGraphAsymmErrors* g = kSoftSlope(label,type);
    double min = *std::min_element(g->GetX(),g->GetX()+g->GetN());
    double max = *std::max_element(g->GetX(),g->GetX()+g->GetN());

    TF1* fit = new TF1(name,"[0] + [1]*log(x/100.) + [2]*sq(log(x/100.)) + [3]*log(x/100.)*log(x/100.)*log(x/100.)",min,max);
    fit->SetParameter(0,1.);
    for(int i = 1; i < fit->GetNpar(); ++i) {
      fit->SetParameter(i,0.);
    }
    fit->SetLineWidth(1);
    g->Fit(fit,"0QR");
    
    return fit;
  }


  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors* EtaBin::kSoftSlope(const SampleLabel &label, FitResult::Type type) const {

    std::vector<double> pt;
    std::vector<double> ptErr;
    std::vector<double> slope;
    std::vector<double> slopeErr;
    for(PtBinIt ptBinIt = ptBins_.begin(); ptBinIt != ptBins_.end(); ++ptBinIt) {
      const Sample* sample = (*ptBinIt)->findSample(label);
      pt.push_back(sample->meanPt(type));
      ptErr.push_back(0.);
      slope.push_back(sample->kSoftSlope(type));
      slopeErr.push_back(sample->kSoftSlopeStatUncert(type));
    }
	  
    return new TGraphAsymmErrors(pt.size(),&(pt.front()),&(slope.front()),
				 &(ptErr.front()),&(ptErr.front()),
				 &(slopeErr.front()),&(slopeErr.front()));
  }


  // -------------------------------------------------------------------------------------
  TF1* EtaBin::kSoftFit(const SampleLabel &label, FitResult::Type type, const TString &name) const {
    return ptBins_.front()->findSample(label)->kSoftFit(type,name);
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::hasSystematicUncertainty(const SampleLabel &label, FitResult::Type type) const {
    bool result = false;
    for(std::set<SystematicUncertainty*>::const_iterator it = systUncerts_.begin();
	it != systUncerts_.end(); ++it) {
      if( (*it)->nominalSampleLabel() == label && (*it)->fitResultType() == type ) {
	result = true;
	break;
      }
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::findSystematicUncertainty(const SampleLabel &label, FitResult::Type type, const SystematicUncertainty* &uncert) const {
    uncert = 0;
    bool result = false;
    for(std::set<SystematicUncertainty*>::const_iterator it = systUncerts_.begin();
	it != systUncerts_.end(); ++it) {
      if( (*it)->nominalSampleLabel() == label && (*it)->fitResultType() == type ) {
	result = true;
	uncert = *it;
	break;
      }
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  SystematicUncertainty* EtaBin::findSystematicUncertainty(const SampleLabel &label, FitResult::Type type) {
    SystematicUncertainty* uncert = 0;

    // Check whether there are uncertainties already associated to this Sample and FitResult
    for(std::set<SystematicUncertainty*>::const_iterator it = systUncerts_.begin();
	it != systUncerts_.end(); ++it) {
      if( (*it)->nominalSampleLabel() == label && (*it)->fitResultType() == type ) {
	uncert = *it;
	break;
      }
    }

    // If not, create new combined uncertainty and add to others
    if( uncert == 0 ) {
      uncert = new SystematicUncertainty("Total",5,label,type);
      systUncerts_.insert(uncert);
    }

    return uncert;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::addExtrapolationUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color) {
    SystematicUncertainty* uncert = findSystematicUncertainty(nominalSample,type);

    TString label = "Extrapolation";
    if( !(uncert->hasComponent(label)) ) {
      std::vector<double> ptMean;
      std::vector<double> ptMin;
      std::vector<double> ptMax;
      std::vector<double> relUncerts;
      for(unsigned int i = 0; i < par_->nPtBins(etaBin()); ++i) {
	ptMean.push_back(meanPt(nominalSample,type,i));
	ptMin.push_back(par_->ptMin(etaBin(),i));
	ptMax.push_back(par_->ptMax(etaBin(),i));
	relUncerts.push_back(std::abs(extrapolatedSystUncert(nominalSample,type,i)/extrapolatedValue(nominalSample,type,i)));
      }
      uncert->addComponent(label,color,nominalSample,type,ptMean,ptMin,ptMax,relUncerts,relUncerts);
    }
    
    return true;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::addPLIUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color) {
    SystematicUncertainty* uncert = findSystematicUncertainty(nominalSample,type);
    
    TString label = "PLI";
    if( !(uncert->hasComponent(label)) ) {
      std::vector<double> ptMean;
      std::vector<double> ptMin;
      std::vector<double> ptMax;
      std::vector<double> relUncertDown;
      std::vector<double> relUncertUp;
      for(unsigned int i = 0; i < par_->nPtBins(etaBin()); ++i) {
 	ptMean.push_back(meanPt(nominalSample,type,i));
 	ptMin.push_back(par_->ptMin(etaBin(),i));
 	ptMax.push_back(par_->ptMax(etaBin(),i));
	double res = correctedResolution(nominalSample,type,i);
	relUncertDown.push_back(std::abs((correctedResolution(nominalSample,type,i,0.75)-res)/res));
	relUncertUp.push_back(std::abs((correctedResolution(nominalSample,type,i,1.25)-res)/res));
      }
      uncert->addComponent(label,color,nominalSample,type,ptMean,ptMin,ptMax,relUncertDown,relUncertUp);
    }
    
    return true;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::addMCClosureUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color) {
    SystematicUncertainty* uncert = findSystematicUncertainty(nominalSample,type);
    
    TString label = "MC Closure";
    if( !(uncert->hasComponent(label)) ) {
      std::vector<double> ptMean;
      std::vector<double> ptMin;
      std::vector<double> ptMax;
      std::vector<double> relUncert;
      for(unsigned int i = 0; i < par_->nPtBins(etaBin()); ++i) {
 	ptMean.push_back(meanPt(nominalSample,type,i));
 	ptMin.push_back(par_->ptMin(etaBin(),i));
 	ptMax.push_back(par_->ptMax(etaBin(),i));
	relUncert.push_back(0.5*std::abs(1.-relativeMCClosure(nominalSample,type,i)));
      }
      uncert->addComponent(label,color,nominalSample,type,ptMean,ptMin,ptMax,relUncert,relUncert);
    }
    
    return true;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::addUncertaintyFromVariedSample(const TString &uncertaintyLabel, double fraction, const SampleLabel &nominalSample, FitResult::Type type, const TString &variedSample, int color) {
    SystematicUncertainty* uncert = findSystematicUncertainty(nominalSample,type);
    
    if( !(uncert->hasComponent(uncertaintyLabel)) ) {
      std::vector<double> ptMean;
      std::vector<double> ptMin;
      std::vector<double> ptMax;
      std::vector<double> relUncert;
      for(unsigned int i = 0; i < par_->nPtBins(etaBin()); ++i) {
 	ptMean.push_back(meanPt(nominalSample,type,i));
 	ptMin.push_back(par_->ptMin(etaBin(),i));
 	ptMax.push_back(par_->ptMax(etaBin(),i));
	double res = correctedResolution(nominalSample,type,i);
	double resVar = correctedResolution(variedSample,type,i);
	relUncert.push_back(fraction*std::abs((resVar-res)/res));
      }
      uncert->addComponent(uncertaintyLabel,color,nominalSample,type,ptMean,ptMin,ptMax,relUncert,relUncert);
    }
    
    return true;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::addUncertaintyFromVariedSample(const TString &uncertaintyLabel, double fraction, const SampleLabel &nominalSample, FitResult::Type type, const TString &variedSampleDown, const TString &variedSampleUp, int color) {
    SystematicUncertainty* uncert = findSystematicUncertainty(nominalSample,type);
    
    if( !(uncert->hasComponent(uncertaintyLabel)) ) {
      std::vector<double> ptMean;
      std::vector<double> ptMin;
      std::vector<double> ptMax;
      std::vector<double> relUncert;
      for(unsigned int i = 0; i < par_->nPtBins(etaBin()); ++i) {
 	ptMean.push_back(meanPt(nominalSample,type,i));
 	ptMin.push_back(par_->ptMin(etaBin(),i));
 	ptMax.push_back(par_->ptMax(etaBin(),i));
	double res = correctedResolution(nominalSample,type,i);
	double resDown = correctedResolution(variedSampleDown,type,i);
	double resUp = correctedResolution(variedSampleUp,type,i);
	relUncert.push_back(fraction*0.5*(std::abs((resDown-res)/res)+std::abs((resUp-res)/res)));
      }
      uncert->addComponent(uncertaintyLabel,color,nominalSample,type,ptMean,ptMin,ptMax,relUncert,relUncert);
    }
    
    return true;
  }



  // Prepare two samples for comparison:
  // - store sample labels of to be compared samples. They can be
  //   referenced in further steps e.g. when calculating scale factors
  // - compute relative sample weights.
  // -------------------------------------------------------------------------------------
  bool EtaBin::compareSamples(const SampleLabel &label1, const SampleLabel &label2) {
    if( sampleTypes_.find(label1) != sampleTypes_.end() &&
	sampleTypes_.find(label2) != sampleTypes_.end() ) {
      compSamples_.insert(new SampleLabelPair(label1,label2));
      for(PtBins::iterator it = ptBins_.begin(); it != ptBins_.end(); ++it) {
	(*it)->weightSampleRelTo(label2,label1);
      }
    } else {
      std::cerr << "ERROR in EtaBin::compareSamples(): no samples with label '" << label1 << "' and '" << label2 << "'" << std::endl;
      exit(1);
    }

    return true;
  }


  // -------------------------------------------------------------------------------------
  void EtaBin::fitKValue(FitResult::Type type) {
    if( nComparedSamples() == 1 ) {
      const SampleLabelPair* slp = *(compSamples_.begin());
      TGraphAsymmErrors* g1 = correctedResolution(slp->label1(),type);
      TGraphAsymmErrors* g2 = correctedResolution(slp->label2(),type);
      TGraphAsymmErrors* g = util::HistOps::createRatioGraph(g1,g2);

      TF1* fit = new TF1("fit","pol0",(*std::min_element(g->GetX(),g->GetX()+g->GetN()))-1.,(*std::max_element(g->GetX(),g->GetX()+g->GetN()))+1);

      // HACK for Eta4 bin, remove first point
      // Should be replace by proper chi2 criterion in Extrapolation.cc
      if( etaBin() == 4 ) {
	std::cout << "\n\n************* HACK IN FITKVALUE *********************\n\n" << std::endl;
	fit->SetRange(g->GetX()[1]-1.,(*std::max_element(g->GetX(),g->GetX()+g->GetN()))+1);
      }

      g->Fit(fit,"0QR");
      kVal_ = fit->GetParameter(0);
      kValStat_ = fit->GetParError(0);

      const SystematicUncertainty* systUncert = 0;
      if( findSystematicUncertainty(slp->label2(),type,systUncert) ) {
	TGraphAsymmErrors* gDown = static_cast<TGraphAsymmErrors*>(g->Clone());
	TGraphAsymmErrors* gUp = static_cast<TGraphAsymmErrors*>(g->Clone());
	for(int i = 0; i < g->GetN(); ++i) {
	  gDown->SetPoint(i,gDown->GetX()[i],gDown->GetY()[i]-systUncert->relUncertDown(i));
	  gUp->SetPoint(i,gUp->GetX()[i],gUp->GetY()[i]+systUncert->relUncertUp(i));
	}
	gDown->Fit(fit,"0QR");
	kValSystDown_ = std::abs(kVal_-fit->GetParameter(0));
	gUp->Fit(fit,"0QR");
	kValSystUp_ = std::abs(kVal_-fit->GetParameter(0));
	delete gDown;
	delete gUp;
      } else {
	kValSystDown_ = 0.;
	kValSystUp_ = 0.;
      }

      kValTotalDown_ = sqrt( kValStat_*kValStat_ + kValSystDown_*kValSystDown_ );
      kValTotalUp_ = sqrt( kValStat_*kValStat_ + kValSystUp_*kValSystUp_ );
      scaledMCTruthReso_ = ResolutionFunction::createScaledResolutionFunction(mcTruthReso_,kVal_,kValStat_,kValTotalDown_,kValTotalUp_);
      kValType_ = type;

      delete g1;
      delete g2;
      delete g;
      delete fit;
    }
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::hasKValue(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type) const {
    const SampleLabelPair* slp = *compSamples_.begin();
    return ( nComparedSamples() == 1 && slp->label1() == label1 && slp->label2() == label2 && type == kValType_ );
  }


  // Return ratio graph
  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors* EtaBin::ratioGraph(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type) const {
    TGraphAsymmErrors* g1 = correctedResolution(label1,type);
    TGraphAsymmErrors* g2 = correctedResolution(label2,type);
    TGraphAsymmErrors* g = util::HistOps::createRatioGraph(g1,g2);
    delete g1;
    delete g2;

    return g;
  }


  // Return TF1 of fitted kValue
  // -------------------------------------------------------------------------------------
  TF1* EtaBin::kValueLine(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type, const TString &name, double xMin, double xMax) const {
    double kVal = kValue(label1,label2,type);
    double kStatErr = kStat(label1,label2,type);
    TF1* kValueLine = new TF1(name,"pol0",xMin,xMax);
    kValueLine->SetLineWidth(1);
    kValueLine->SetLineColor(kBlack);
    kValueLine->SetParameter(0,kVal);
    kValueLine->SetParError(0,kStatErr);

    return kValueLine;
  }


  // Return band of statistical uncertainty of fitted kValue
  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors* EtaBin::kValueStatBand(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type, double xMin, double xMax) const {
    double kVal = kValue(label1,label2,type);
    double kStatErr = kStat(label1,label2,type);
    double x = 0.5*(xMin+xMax);
    double xe = std::abs(x-xMin);
    TGraphAsymmErrors* kStatBand = new TGraphAsymmErrors(1,&x,&kVal,&xe,&xe,&kStatErr,&kStatErr);
    //kStatBand->SetFillStyle(3013);
    kStatBand->SetFillStyle(1001);
    kStatBand->SetFillColor(40);
    kStatBand->SetLineColor(kStatBand->GetFillColor());

    return kStatBand;
  }


  // Return band of systematic of fitted kValue
  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors* EtaBin::kValueSystBand(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type, double xMin, double xMax) const {
    double kVal = kValue(label1,label2,type);
    double kSystErrDown = kSystDown(label1,label2,type);
    double kSystErrUp = kSystUp(label1,label2,type);
    double x = 0.5*(xMin+xMax);
    double xe = std::abs(x-xMin);
    TGraphAsymmErrors* kSystBand = new TGraphAsymmErrors(1,&x,&kVal,&xe,&xe,&kSystErrDown,&kSystErrUp);
    kSystBand->SetFillStyle(1001);
    kSystBand->SetFillColor(kYellow);
    const SystematicUncertainty* uncert = 0;
    if( findSystematicUncertainty(label2,type,uncert) ) {
      kSystBand->SetFillColor(uncert->color());
    }
    kSystBand->SetLineColor(kSystBand->GetFillColor());

    return kSystBand;
  }


  // -------------------------------------------------------------------------------------
  Sample::Type EtaBin::sampleType(const SampleLabel &label) const {
    SampleTypeIt it = sampleTypes_.find(label);
    if( it == sampleTypes_.end() ) {
      std::cerr << "ERROR in EtaBin::sampleType(): No sample with label '" << label << "'" << std::endl;
      exit(1);
    }

    return it->second;
  }


  //! \brief Fit the intrinsic particle level imbalance
  //!
  //! This will add a Sample of type MCTruth, which has the correct
  //! FitResult::Type PtGenAsym. The fitted ptGen resolution is used
  //! as PLI correction to compensate for jet hadronisation and
  //! fragmentation effects.
  //!
  //! The value of the fitted ptGen resolution can can be computed
  //! via EtaBin::pli, the fitted function via EtaBin::pliFunc.
  // -------------------------------------------------------------------------------------
  void EtaBin::fitPLI(const TString &label, const TString &fileName, ResolutionFunction::Type type) {

    // Add Samples containing ptGen asymmetry histograms
    // They will per default have the correct FitResult::Type PtGenAsym
    for(std::vector<PtBin*>::iterator it = ptBins_.begin(); it != ptBins_.end(); ++it) {
      if( !((*it)->addMCTruthSample(label,fileName)) ) {
	std::cerr << "ERROR in EtaBin::fitPLI(): Error adding MCTruth Sample '" << label << "'" << std::endl;
	exit(1);
      }
    }

    // Create graph of extrapolated ptGen asymmetry width
    std::vector<double> pt;
    std::vector<double> ptErr;
    std::vector<double> res;
    std::vector<double> resStatErr;
    for(PtBinIt ptBinIt = ptBinsBegin(); ptBinIt != ptBinsEnd(); ++ptBinIt) {
      const Sample* sample = (*ptBinIt)->mcTruthSample();
      pt.push_back(sample->meanPt(FitResult::PtGenAsym));
      ptErr.push_back(0.);
      res.push_back(sample->extrapolatedValue(FitResult::PtGenAsym));
      resStatErr.push_back(sample->extrapolatedStatUncert(FitResult::PtGenAsym));
    }
    TGraphAsymmErrors* g = new TGraphAsymmErrors(pt.size(),&(pt.front()),&(res.front()),
						 &(ptErr.front()),&(ptErr.front()),
						 &(resStatErr.front()),&(resStatErr.front()));

    // Fit new PLI from extrapolated ptGen asymmetry width
    delete pli_;
    std::cout << "Fitting PLI for EtaBin " << etaBin() << std::endl;
    pli_ = ResolutionFunction::fitTGraph(g,type);
  }


  // -------------------------------------------------------------------------------------
  void EtaBin::setPLI(ResolutionFunction* pli) { 
    delete pli_;
    pli_ = pli;
  }


  // -------------------------------------------------------------------------------------
  void EtaBin::setMCTruthResolution(ResolutionFunction* mcTruthReso) { 
    delete mcTruthReso_;
    mcTruthReso_ = mcTruthReso;
  }


  // -------------------------------------------------------------------------------------
  double EtaBin::meanPt(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
    double result = 0.;
    for(PtBinIt ptBinIt = ptBins_.begin(); ptBinIt != ptBins_.end(); ++ptBinIt) {
      if( (*ptBinIt)->ptBin() == ptBinIdx ) {
	result = (*ptBinIt)->findSample(label)->meanPt(type);
	break;
      }
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  double EtaBin::meanPtStatUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
    double result = 0.;
    for(PtBinIt ptBinIt = ptBins_.begin(); ptBinIt != ptBins_.end(); ++ptBinIt) {
      if( (*ptBinIt)->ptBin() == ptBinIdx ) {
	result = (*ptBinIt)->findSample(label)->meanPtStatUncert(type);
	break;
      }
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  double EtaBin::fittedValue(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx, unsigned int ptSoftIdx) const {
    double result = 0.;
    for(PtBinIt ptBinIt = ptBins_.begin(); ptBinIt != ptBins_.end(); ++ptBinIt) {
      if( (*ptBinIt)->ptBin() == ptBinIdx ) {
	result = (*ptBinIt)->findSample(label)->fittedValue(type,ptSoftIdx);
	break;
      }
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  double EtaBin::fittedStatUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx, unsigned int ptSoftIdx) const {
    double result = 0.;
    for(PtBinIt ptBinIt = ptBins_.begin(); ptBinIt != ptBins_.end(); ++ptBinIt) {
      if( (*ptBinIt)->ptBin() == ptBinIdx ) {
	result = (*ptBinIt)->findSample(label)->fittedUncert(type,ptSoftIdx);
	break;
      }
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  double EtaBin::extrapolatedValue(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
    double result = 0.;
    for(PtBinIt ptBinIt = ptBins_.begin(); ptBinIt != ptBins_.end(); ++ptBinIt) {
      if( (*ptBinIt)->ptBin() == ptBinIdx ) {
	result = (*ptBinIt)->findSample(label)->extrapolatedValue(type);
	break;
      }
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  double EtaBin::extrapolatedStatUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
    double result = 0.;
    for(PtBinIt ptBinIt = ptBins_.begin(); ptBinIt != ptBins_.end(); ++ptBinIt) {
      if( (*ptBinIt)->ptBin() == ptBinIdx ) {
	result = (*ptBinIt)->findSample(label)->extrapolatedStatUncert(type);
	break;
      }
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  double EtaBin::extrapolatedSystUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
    double result = 0.;
    for(PtBinIt ptBinIt = ptBins_.begin(); ptBinIt != ptBins_.end(); ++ptBinIt) {
      if( (*ptBinIt)->ptBin() == ptBinIdx ) {
	result = (*ptBinIt)->findSample(label)->extrapolatedSystUncert(type);
	break;
      }
    }

    return result;
  }


  // Should be integral over bin!
  // -------------------------------------------------------------------------------------
  double EtaBin::mcTruthResolution(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
    return mcTruthReso_->val(meanPt(label,type,ptBinIdx));
  }


  // Should be integral over bin!
  // -------------------------------------------------------------------------------------
  double EtaBin::pli(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
    return pli_->val(meanPt(label,type,ptBinIdx));
  }


  // -------------------------------------------------------------------------------------
  double EtaBin::correctedResolution(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx, double scalePLI) const {
    double result = extrapolatedValue(label,type,ptBinIdx);
    double plim = scalePLI*pli(label,type,ptBinIdx);
    double sqDiff = result*result - plim*plim;
    return sqDiff < 0. ? 0.01 : sqrt(sqDiff);
  }

  
  // -------------------------------------------------------------------------------------
  double EtaBin::correctedResolutionStatUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
    return extrapolatedStatUncert(label,type,ptBinIdx);
  }


  //! Returns MCTruth / Measurement
  // -------------------------------------------------------------------------------------
  double EtaBin::relativeMCClosure(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
    double meas = correctedResolution(label,type,ptBinIdx);
    double truth = mcTruthResolution(label,type,ptBinIdx);
    double result = 0.;
    if( meas > 0. ) result = truth/meas;

    return result;
  }


  //! Create graph of extrapolated resolution
  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors* EtaBin::extrapolatedResolution(const SampleLabel &label, FitResult::Type type) const {
    std::vector<double> pt;
    std::vector<double> ptErr;
    std::vector<double> res;
    std::vector<double> resStatErr;
    for(PtBinIt ptBinIt = ptBinsBegin(); ptBinIt != ptBinsEnd(); ++ptBinIt) {
      unsigned int ptBinIdx = (*ptBinIt)->ptBin();
      pt.push_back(meanPt(label,type,ptBinIdx));
      ptErr.push_back(meanPtStatUncert(label,type,ptBinIdx));
      res.push_back(extrapolatedValue(label,type,ptBinIdx));
      resStatErr.push_back(extrapolatedStatUncert(label,type,ptBinIdx));
    }

    return new TGraphAsymmErrors(pt.size(),&(pt.front()),&(res.front()),&(ptErr.front()),&(ptErr.front()),
				 &(resStatErr.front()),&(resStatErr.front()));
  }



  //! Create graph of extrapolated and PLI-corrected resolution
  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors* EtaBin::correctedResolution(const SampleLabel &label, FitResult::Type type) const {
    std::vector<double> pt;
    std::vector<double> ptErr;
    std::vector<double> res;
    std::vector<double> resStatErr;
    for(PtBinIt ptBinIt = ptBinsBegin(); ptBinIt != ptBinsEnd(); ++ptBinIt) {
      unsigned int ptBinIdx = (*ptBinIt)->ptBin();
      pt.push_back(meanPt(label,type,ptBinIdx));
      ptErr.push_back(meanPtStatUncert(label,type,ptBinIdx));
      res.push_back(correctedResolution(label,type,ptBinIdx));
      resStatErr.push_back(correctedResolutionStatUncert(label,type,ptBinIdx));
    }

    return new TGraphAsymmErrors(pt.size(),&(pt.front()),&(res.front()),&(ptErr.front()),&(ptErr.front()),
				 &(resStatErr.front()),&(resStatErr.front()));
  }


  //! Return graph of bias corrected measurement
  //! - label1: Sample which is corrected
  //! - label2: Sample from which the bias is determined
  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors* EtaBin::biasCorrectedResolution(const SampleLabel &label1, const SampleLabel &label2, FitResult::Type type) const {
    TGraphAsymmErrors* g = correctedResolution(label1,type);
    for(unsigned int i = 0; i < g->GetN(); ++i) {
      g->SetPoint(i,g->GetX()[i],g->GetY()[i]*relativeMCClosure(label2,type,i));
    }

    return g;
  }


  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors* EtaBin::scaledMCTruthUncertaintyBand() const {
    int n = 1000;
    double dx = (scaledMCTruthReso_->ptMax()-scaledMCTruthReso_->ptMin())/n;
    std::vector<double> x;
    std::vector<double> xe;
    std::vector<double> y;
    std::vector<double> yed;
    std::vector<double> yeu;
    for(int i = 0; i < n; ++i) {
      x.push_back(scaledMCTruthReso_->ptMin()+(i+0.5)*dx);
      xe.push_back(0.);
      y.push_back(scaledMCTruthReso_->val(x.back()));
      yed.push_back(scaledMCTruthReso_->uncertDown(x.back()));
      yeu.push_back(scaledMCTruthReso_->uncertUp(x.back()));
    }
    TGraphAsymmErrors* g = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),
						 &(xe.front()),&(xe.front()),&(yed.front()),&(yeu.front()));
    g->SetFillStyle(1001);
    g->SetFillColor(5);
    g->SetLineColor(g->GetFillColor());

    return g;
  }


  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors* EtaBin::scaledMCTruthRatioBand() const {
    int n = 1000;
    double dx = (scaledMCTruthReso_->ptMax()-scaledMCTruthReso_->ptMin())/n;
    std::vector<double> x;
    std::vector<double> xe;
    std::vector<double> y;
    std::vector<double> yed;
    std::vector<double> yeu;
    for(int i = 0; i < n; ++i) {
      x.push_back(scaledMCTruthReso_->ptMin()+(i+0.5)*dx);
      xe.push_back(0.);
      y.push_back(1.);
      yed.push_back(scaledMCTruthReso_->relUncertDown(x.back()));
      yeu.push_back(scaledMCTruthReso_->relUncertUp(x.back()));
    }
    TGraphAsymmErrors* g = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),
						 &(xe.front()),&(xe.front()),&(yed.front()),&(yeu.front()));
    g->SetFillStyle(1001);
    g->SetFillColor(5);
    g->SetLineColor(g->GetFillColor());

    return g;
  }


  //! Return TGraphAsymmErrors of the standard deviations of a Gaussian fit to 
  //! the asymmetry histogram in each ptBin, for the specified ptSoftBin
  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors* EtaBin::asymmetryWidths(const SampleLabel &label, unsigned int ptSoftBin) const {
    std::vector<double> pt;
    std::vector<double> ptErr;
    std::vector<double> width;
    std::vector<double> widthStatErr;
    for(PtBinIt ptBinIt = ptBinsBegin(); ptBinIt != ptBinsEnd(); ++ptBinIt) {
      const Sample* s = (*ptBinIt)->findSample(label);
      pt.push_back(s->meanPtAve(ptSoftBin));
      ptErr.push_back(s->meanPtAveStatUncert(ptSoftBin));
      width.push_back(s->asymmetryWidth(ptSoftBin));
      widthStatErr.push_back(s->asymmetryWidthStatUncert(ptSoftBin));
    }

    return new TGraphAsymmErrors(pt.size(),&(pt.front()),&(width.front()),&(ptErr.front()),&(ptErr.front()),
				 &(widthStatErr.front()),&(widthStatErr.front()));
  }
}

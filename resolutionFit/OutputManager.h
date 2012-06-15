// $Id: OutputManager.h,v 1.6 2012/06/01 18:32:55 mschrode Exp $

#ifndef OUTPUT_MANAGER_H
#define OUTPUT_MANAGER_H

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TPad.h"
#include "TPostScript.h"
#include "TString.h"
#include "TVirtualPad.h"

#include "../util/HistOps.h"


namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  class OutputManager {
  public:
    enum Mode { PSAllInOne, EPSSingleFiles, EPSSingleFilesPlusROOT };

    static OutputManager* createOutputManager(OutputManager::Mode mode, const TString &fileNameBase);
    static bool isValidMode(Mode mode);

    OutputManager(const TString &fileNameBase);
    virtual ~OutputManager();

    virtual void logx() = 0;
    virtual void logy() = 0;
    virtual void newPage(const TString &title) = 0;
    virtual void nextMultiPad(const TString &title) = 0;
    virtual void nextPad(const TString &title) = 0;
    virtual void nextMainPad(const TString &title) = 0;
    virtual void nextRatioPad() = 0;
    virtual void saveCurrentPad(const TString &name) = 0;
    virtual void saveCanvas(TCanvas* can, const TString &name) const = 0;

    TH1* mainFrame(double xMin, double xMax, double yMin, double yMax, const TString &yTitle) const {
      return util::HistOps::createRatioTopFrame(xMin,xMax,yMin,yMax,yTitle);
    }
    TH1* mainFrame(const TH1* h) const;
    TH1* ratioFrame(const TH1 *h, const TString &xTitle, const TString &xUnit, const TString &yTitle, double yMin, double yMax) const {
      return util::HistOps::createRatioBottomFrame(h,xTitle,xUnit,yTitle,yMin,yMax);
    }
    TH1* ratioFrame(const TH1 *h, const TString &xTitle, const TString &xUnit, double yMin, double yMax) const {
      return util::HistOps::createRatioBottomFrame(h,xTitle,xUnit,yMin,yMax);
    }
    TH1* ratioFrame(const TH1 *h, double yMin, double yMax) const {
      return util::HistOps::createRatioBottomFrame(h,yMin,yMax); 
    }


  protected:
    TCanvas* can_;
    TCanvas* topCan_;
    TPad* bottomPad_;
    TVirtualPad* lastPad_;
  };


  //! Plots are stored in one common .ps file
  //! Several multipads are put on one page
  // -------------------------------------------------------------------------------------
  class OutputManagerPSAllInOne : public OutputManager {
  public:
    OutputManagerPSAllInOne(const TString &fileNameBase);
    ~OutputManagerPSAllInOne();

    void logx();
    void logy();
    void newPage(const TString &title);
    void nextMultiPad(const TString &title);
    void nextPad(const TString &title);
    void nextMainPad(const TString &title);
    void nextRatioPad();
    void saveCurrentPad(const TString &name);
    void saveCanvas(TCanvas* can, const TString &name) const {}


  private:
    TCanvas* multiCan_;
    int multiCanPadIdx_;
    TPostScript* ps_;
  };


  //! Each plot is stored in a separate .eps file; all
  //! plots are stored as TCanvas objects in one common
  //! .root file
  // -------------------------------------------------------------------------------------
  class OutputManagerEPSSingleFiles : public OutputManager {
  public:
    OutputManagerEPSSingleFiles(const TString &fileNameBase, bool rootOutput);
    ~OutputManagerEPSSingleFiles();

    void logx();
    void logy();
    void newPage(const TString &title) {};
    void nextMultiPad(const TString &title) { nextPad(title); }
    void nextPad(const TString &title);
    void nextMainPad(const TString &title);
    void nextRatioPad();
    void saveCurrentPad(const TString &name);
    void saveCanvas(TCanvas* can, const TString &name) const;

  private:
    const bool rootOutput_;
    TFile* outFile_;
  };
}
#endif

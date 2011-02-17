// $Id: OutputManager.h,v 1.1 2011/02/15 18:22:25 mschrode Exp $

#ifndef OUTPUT_MANAGER_H
#define OUTPUT_MANAGER_H

#include "TCanvas.h"
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
    enum Mode { PSAllInOne, EPSSingleFiles };

    static OutputManager* createOutputManager(OutputManager::Mode mode, const TString &fileNameBase);

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

    TH1* mainFrame(double xMin, double xMax, double yMin, double yMax, const TString &yTitle) const {
      return util::HistOps::createRatioTopFrame(xMin,xMax,yMin,yMax,yTitle); }
    TH1* ratioFrame(const TH1 *h, const TString &xTitle, const TString &xUnit, double yMin, double yMax) const {
      return util::HistOps::createRatioBottomFrame(h,xTitle,xUnit,yMin,yMax); }


  protected:
    TCanvas* can_;
    TCanvas* topCan_;
    TPad* bottomPad_;
    TVirtualPad* lastPad_;
  };


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


  private:
    TCanvas* multiCan_;
    int multiCanPadIdx_;
    TPostScript* ps_;
  };


  // -------------------------------------------------------------------------------------
  class OutputManagerEPSSingleFiles : public OutputManager {
  public:
    OutputManagerEPSSingleFiles(const TString &fileNameBase);

    void logx();
    void logy();
    void newPage(const TString &title) {};
    void nextMultiPad(const TString &title) { nextPad(title); }
    void nextPad(const TString &title);
    void nextMainPad(const TString &title);
    void nextRatioPad();
    void saveCurrentPad(const TString &name);
  };
}
#endif

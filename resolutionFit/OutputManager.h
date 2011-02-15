// $Id: $

#ifndef OUTPUT_MANAGER_H
#define OUTPUT_MANAGER_H

#include "TCanvas.h"
#include "TPad.h"
#include "TPostScript.h"
#include "TString.h"
#include "TVirtualPad.h"


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
    virtual void saveCurrentPad(const TString &name) = 0;


  protected:
    TCanvas* can_;
  };


  // -------------------------------------------------------------------------------------
  class OutputManagerPSAllInOne : public OutputManager {
  public:
    OutputManagerPSAllInOne(const TString &fileNameBase);
    ~OutputManagerPSAllInOne();

    void logx() { lastPad_->SetLogx(1); }
    void logy() { lastPad_->SetLogy(1); }
    void newPage(const TString &title);
    void nextMultiPad(const TString &title);
    void nextPad(const TString &title);
    void saveCurrentPad(const TString &name);


  private:
    TCanvas* multiCan_;
    TVirtualPad* lastPad_;
    int multiCanPadIdx_;
    TPostScript* ps_;
  };


  // -------------------------------------------------------------------------------------
  class OutputManagerEPSSingleFiles : public OutputManager {
  public:
    OutputManagerEPSSingleFiles(const TString &fileNameBase);

    void logx() { can_->SetLogx(1); }
    void logy() { can_->SetLogy(1); }
    void newPage(const TString &title) {};
    void nextMultiPad(const TString &title) { nextPad(title); }
    void nextPad(const TString &title);
    void saveCurrentPad(const TString &name);
  };
}
#endif

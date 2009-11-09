#ifndef TControlPlots_h
#define TControlPlots_h

#include <set>
#include <string>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TObject.h"
#include "TStyle.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"




//!  \brief Create control plots
//!
//!  Objects of \p TControlPlots can create control plots via
//!  the \p makeControlPlots() method from several TData
//!  objects. The output is in .ps or both .ps and .root format.
//!  The kind of the control plots and the output format is
//!  specified via the config file.
//!
//!  \author Christian Autermann
//!  \date Fri Jan 18 13:55:15 2008 UTC
//!  $Id: ControlPlots.h,v 1.1 2009/10/30 08:59:42 mschrode Exp $
// -------------------------------------------------------------
class TControlPlots
{
public:
  TControlPlots(const std::string& configfile, const std::vector<TData*> *data, TParameters *par);
  ~TControlPlots();

  void makePlots();
  bool outputFormatRoot() const { return outputROOT_; }

 private:  
  void makeControlPlotsBinnedResponse();
  void makeControlPlotsChi2();
  void makeControlPlotsJetTruthEventResponse();
  void makeControlPlotsL2L3MCTruth();
  void makeControlPlotsParameterScan();
  void makeControlPlotsTop();
  void makeControlPlotsTwoJetsPtBalance();

  bool equidistLogBins(double * bins, int nBins, double first, double last) const;
  void findYRange(const TH1F * h, double& min, double& max) const;
  void fit2D(const TH2F* hist, std::vector<TH1F*>& hresults, std::vector<TH1F*>& distributions,
	     std::vector<TF1*>& gaussFits, const bool plotgauss = true) const;
  void fit2D(const TH2F* hist, TH1F* hresults[8], TH1F* gaussplots[4], TF1* gf[4], const bool plotgauss=true) const;
  void fit2D(const TH2F* hist, std::vector<TH1F*>& hresults) const { std::vector<TH1F*> v1; std::vector<TF1*> v2; fit2D(hist, hresults, v1, v2, false); }
  void fit2DMean(const TH2F* hist, TH1F*& hresult,
		 std::vector<TH1F*>& distributions, int color) const;
  void fit2DMean(const TH2F* hist, TH1F*& hresult, int color) const;
  void fit2DGaussMean(const TH2F* hist, TH1F*& hresult,
		      std::vector<TH1F*>& distributions,
		      std::vector<TF1*>& gaussFits, int color) const;
  void fit2DGaussMean(const TH2F* hist, TH1F*& hresult, int color) const;
  std::vector<double> getEtaBinEdges(int binningModel = 0) const;
  bool readJetMETParameters();
  void resetFittedParameters();
  void setGStyle() const;
  void setYRange(TH1F * h, double c1 = 0.9, double c2 = 1.1, double minLimit = 0.) const;
  void writeToRootFile(std::vector<TObject*> obj, std::string dir);

  ConfigFile  *config_;                         //!< Pointer to config file
  const std::vector<TData*> *data_;             //!< Pointer to data
  std::vector<double> fittedPar_;               //!< Stores fitted parameters
  TFile       *outFile_;                        //!< Pointer to root output file
  bool         outputROOT_;                     //!< If true, histograms are written to ROOT file
  TParameters *par_;                            //!< Pointer to parameter values



  //!  \brief A two-dimensional grid
  //!
  //!  The two dimensions of the grid are named 'x' and 'y'.
  //!  They are divided into bins whose borders can be specified
  //!  i.e. the binsize can vary. The index of the bins in x
  //!  direction is named 'ix' and ranges from 0 to NBinsX() - 1
  //!  and likewise for the y direction. Additionally, there is
  //!  a global bin numbering scheme 'bin' from 0 to NBins() - 1,
  //!  counting the bins first in x and then in y direction i.e.
  //!  
  //!     bin    ix  iy
  //!  -------------------
  //!      0      0   0
  //!      1      1   0
  //!     ...
  //!   NBinsX()  0   1
  //!  
  //!  \author Matthias Schroeder
  //!  \date Thu Apr 23 13:05:54 CEST 2009
  // -------------------------------------------------------------
  class Binning
    {
    public:
      //!  \brief Creates a binning in (x,y)
      //!
      //!  The bins are created from the given bin edges.
      //!  \note The vectors of bin edges must have at least two entries
      //!        and must be ordered from lowest to highest value.
      //!  \param binEdgesX Bin edges in x direction (binEdgesX.size() == NBinsX()+1)
      //!  \param binEdgesY Bin edges in y direction (binEdgesY.size() == NBinsY()+1)
      // -------------------------------------------------------------
      Binning(const std::vector<double>& binEdgesX, const std::vector<double>& binEdgesY);
      ~Binning() {};

      //!  \brief Lower bin edge in x direction
      //!  \param bin Global bin index
      //!  \return Lower bin edge in x direction
      // -------------------------------------------------------------
      double xLow(int bin) const { return edgesX_.at(iX(bin)); }

      //!  \brief Upper bin edge in x direction
      //!  \param bin Global bin index
      //!  \return Upper bin edge in x direction
      // -------------------------------------------------------------
      double xUp(int bin) const { return edgesX_.at(iX(bin)+1); }

      //!  \brief Lower bin edge in y direction
      //!  \param bin Global bin index
      //!  \return Lower bin edge in y direction
      // -------------------------------------------------------------
      double yLow(int bin) const { return edgesY_.at(iY(bin)); }

      //!  \brief Upper bin edge in y direction
      //!  \param bin Global bin index
      //!  \return Upper bin edge in y direction
      // -------------------------------------------------------------
      double yUp(int bin) const { return edgesY_.at(iY(bin)+1); }

      //!  \brief X bin index ix of global bin
      //!
      //!  Finds the index ix of the bin in x direction corresponding
      //!  to a global bin.
      //!  \param bin Global bin index
      //!  \return ix of global bin
      // -------------------------------------------------------------
      int iX(int bin) const { return bin % nBinsX(); }

      //!  \brief X bin index ix of value x
      //!
      //!  Finds the index ix of the bin in x direction that contains
      //!  the value x.
      //!  \param x x value
      //!  \return ix of the bin containing x
      //!          ( ix == -1 for x < XLow(0): Underflow, 
      //!            ix == NBinsX() for x > XUp(NBins()-1): Overflow )
      // -------------------------------------------------------------
      int iX(double x) const;

      //!  \brief Y bin index iy of global bin
      //!
      //!  Finds the index iy of the bin in y direction corresponding
      //!  to a global bin.
      //!  \param bin Global bin index
      //!  \return iy of global bin
      // -------------------------------------------------------------
      int iY(int bin) const { return bin / nBinsX(); }

      //!  \brief Y bin index iy of value y
      //!
      //!  Finds the index iy of the bin in y direction that contains
      //!  the value y.
      //!  \param y y value
      //!  \return iy of the bin containing y
      //!          ( iy == -1 for y < YLow(0): Underflow, 
      //!            iy == NBinsyY() for y > YUp(NBins()-1): Overflow )
      // -------------------------------------------------------------
      int iY(double y) const;

      //!  \brief Global bin index of bin with x and y indices ix and iy
      //!  \param ix Bin index in x direction
      //!  \param iy Bin index in y direction
      //!  \return Global bin index
      // -------------------------------------------------------------
      int bin(int ix, int iy) const;

      //!  \brief Global bin index of bin containing values x and y
      //!  \param x x value
      //!  \param y y value
      //!  \return Global bin index
      // -------------------------------------------------------------
      int bin(double x, double y) const { return bin(iX(x),iY(y)); }

      //!  \brief Number of global bins
      //!
      //!  NBins() == NBinsX() * NBinsY()
      //!  \return Number of global bins
      // -------------------------------------------------------------
      int nBins() const { return  nBinsX()*nBinsY(); }

      //!  \brief Number of x bins
      //!  \return Number of x bins
      // -------------------------------------------------------------
      int nBinsX() const { return static_cast<int>(edgesX_.size()) - 1; }

      //!  \brief Number of y bins
      //!  \return Number of y bins
      // -------------------------------------------------------------
      int nBinsY() const { return static_cast<int>(edgesY_.size()) - 1; }

      //!  \brief Print the binning to std
      // -------------------------------------------------------------
      void print() const;

    private:
      std::vector<double> edgesX_;   //!< Bin edges in x direction
      std::vector<double> edgesY_;   //!< Bin edges in y direction
    };
};
#endif

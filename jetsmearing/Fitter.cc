// $Id: Fitter.cc,v 1.8 2009/05/08 15:41:42 mschrode Exp $

#include "Fitter.h"

#include <cmath>
#include <iostream>

namespace js
{
  //!  \brief Constructor
  //!
  //!  Takes the assumed model for the response pdf (see
  //!  class description). If the number of specified 
  //!  start parameters is not equal to the number of
  //!  parameters of that model, the start parameters
  //!  are adjusted in the following way:
  //!   - If par.size() > GetNPar(model):
  //!     The first GetNPar(model) parameters are taken
  //!     and the rest is ignored
  //!   - If par.size() < GetNPar(model):
  //!     The par.size() parameters are taken and the
  //!     rest of necessary parameters is initialized
  //!     to 1.     
  //!
  //!  \param data   Input Data for the fit
  //!  \param min    Minimum of considered \f$ E^{true}_{T}\f$ spectrum
  //!  \param max    Maximum of considered \f$ E^{true}_{T}\f$ spectrum
  //!  \param model  The assumed response model
  //!  \param par    Start values of parameters
  //----------------------------------------------------------
  Fitter::Fitter(Data& data, double min, double max,
	 const std::string& model, const std::vector<double>& par)
    : mData(data),
      mMin(min),
      mMax(max),
      mModel(model),
      mPar(par),
      mPDFHistMin(0.),
      mPDFHistMax(2.),
      mNBadEvts(0)
  {
    // Adjust number of parameters to fit model
    int n = GetNPar(model);
    if( GetNPar() > n )
      {
	mPar.erase(mPar.begin()+n,mPar.end());
      }
    else if( GetNPar() < n )
      {
	while( GetNPar() < n )
	  {
	    mPar.push_back(1.);
	  }
      }

    assert( GetNPar() == GetNPar(model) );

    // Set up histogram for models with 
    // histogramed pdfs
    if( mModel == "Hist" )
      {
	mPDFHist = new TH1D("mPDFHist","",GetNPar(),mPDFHistMin,mPDFHistMax);
	int minbin = mPDFHist->FindBin(0.9);
	int maxbin = mPDFHist->FindBin(1.1);
	for(int i = 0; i < GetNPar(); i++)
	  {
	    int bin = 1+i;
	    if( bin >= minbin && bin <= maxbin ) mPar.at(i) = 1.;
	    else                                 mPar.at(i) = 0.01;
	  }	
      }
    else if( mModel == "HistGauss" )
      {
	mPDFHist = new TH1D("mPDFHistGauss","",GetNPar()-2,mPDFHistMin,mPDFHistMax);
      }
    else
      {
	mPDFHist = 0;
      }

    std::cout << "Initializing Fitter with the following start parameters:" << std::endl;
    for(int i = 0; i < GetNPar(); i++)
      {
	std::cout << " " << i << "   " << mPar.at(i) << std::endl;
      }
  }



  //----------------------------------------------------------
  Fitter::~Fitter()
  {
    if( mPDFHist ) delete mPDFHist;
  }
  

  //!  \brief Do the fit
  //!
  //!  This function has to be called to do the fit.
  //!  It calls the LVMINI program.
  //----------------------------------------------------------
  void Fitter::Fit()
  {
    ////// LVMINI //////
    int NPAR =  GetNPar();
    int NITER = 1000;
    int MVEC;
    NPAR < 29 ? MVEC = NPAR : MVEC = 29;
    int NAUX = 10000;
    
    std::vector<double> AUX(NAUX, 0.);
    double FSUM;
    double FOPT, FEDM, DUMMY;

    float eps  = float(1.E-3);
    float wlf1 = 1.E-4;
    float wlf2 = 0.9;
    
    lvmeps_(eps, wlf1, wlf2);
    int NPARNEG = -abs(NPAR);
    if(lvmdim_(NPAR, MVEC) > NAUX)
      {
	std::cout << "Aux field too small: " << NAUX << "<" << lvmdim_(NPAR, MVEC) << " needed entries" << std::endl;
      }
    lvmini_(NPARNEG, MVEC, NITER, &AUX.front());
    int IFLAG;
    int IRET = 0;
    int ITER = 0;
    do {
      ++ITER;
      FSUM = NLogLSum(AUX);
      lvmfun_(&mPar.front(), FSUM, IRET, &AUX.front());
      lvmprt_(2, &AUX.front(), 2); //print out

    } while (IRET < 0 && ITER < NITER);
    int error_index = 2;
    error_index = lvmind_(error_index);

    std::cout << "\n" << mNBadEvts << " events flagged as bad." << std::endl;
  }


  //!  \brief Return the negative Log-Likelihood
  //----------------------------------------------------------
  double Fitter::NLogLSum(std::vector<double>& grad)
  {
    double L = 0.;
    int nFlaggedBad = 0;

    for(DataIt datait = mData.begin(); datait != mData.end(); datait++)
      {
	Event * evt = *datait;

	// Check status of event
	if( evt->Status() != 0 )
	  {
	    nFlaggedBad++;
	    continue;
	  }

	// Calculate probability
	double prob = NLogL(evt);
	L += prob;
	    
	// calculate gradients
	double eps = 1.e-10;
	for(int i = 0; i < GetNPar(); i++)
	  {
	    double oldPar = mPar.at(i);
	    mPar.at(i)    = oldPar + eps;

	    double prob1  = NLogL(evt);
	    mPar.at(i)    = oldPar - eps;

	    double prob2  = NLogL(evt);
	    mPar.at(i)    = oldPar;

	    grad.at(i)             += (prob1 - prob2) / 2. / eps;
	    grad.at(i + GetNPar()) += (prob1 + prob2 - 2*prob) / 4. / eps / eps;
	  }
      }

    mNBadEvts += (nFlaggedBad - mNBadEvts);

    return L;
  }



  //!  \brief Returns the negative logarithm of the probability
  //!         for an event
  //!
  //!  Depending on the event type, NLogLDiJet() or 
  //!  NLogLPhotonJet() are called.
  //!
  //!  Parameter limits are enforced by introducing a
  //!  penalty term if the parameter values exceed the
  //!  limits.
  //!
  //!  \param evt Pointer to the event
  //!  \return The negative logarithm of the probability
  //!          associated with this event
  //----------------------------------------------------------
  double Fitter::NLogL(Event * evt) const
  {
    double L = 0.;

    if( evt->Type() == "DiJetEvent" )
      {
	DiJetEvent * dijet = static_cast<DiJetEvent*>(evt);
	L = NLogLDiJet(dijet);
      }
    else if( evt->Type() == "PhotonJetEvent" )
      {
	PhotonJetEvent * photonjet = static_cast<PhotonJetEvent*>(evt);
	L = NLogLPhotonJet(photonjet);
      }

    // Check parameter ranges and set penalty term
    double min = 0.;
    double max = 1.;
    for(int i = 0; i < GetNPar(); i++)
      {
	if( mPar.at(i) < 0. ) L += (min - mPar.at(i));
      }
    if( mPar.at(0) > 1. ) L+= ( mPar.at(0) - max );

    return L;
  }



  //----------------------------------------------------------
  double Fitter::NLogLDiJet(DiJetEvent * dijet) const
  {
    int    maxNIter  = 5;       // Max number of iterations in interval splitting
    double eps = 1.e-5;

    // Normalize integral
    double norm = ( pow(mMax,3) - pow(mMin,3) ) / 3.;

    double         pint      = 0.;      // Current value of integral over response pdfs
    double         pint_old  = 1.;      // Value of integral over response pdfs from previous iteration
    int            nIter     = 0;       // Current iteration in interval splitting
    std::vector<double> pp;             // Product of current function values of response pdfs
    std::vector<double> pp_old;         // Product of function values of response pdfs from previous iteration
    double         h = mMax - mMin;     // Integration interval


    // Iterate until precision or max. number iterations reached
    while(fabs((pint - pint_old) / pint_old) > eps && nIter < maxNIter)
      {
	pint_old = pint;
	pint     = 0;
	pp_old   = pp;
	pp.clear();
	h       /= 3;    // In each iteration, split h into 3 new intervals

	// Loop over nodes xi i.e. interval borders
	for(int i = 0; i <= pow(3,nIter+1); i++)
	  {
	    double t = mMin + i * h;  // Pt at node

	    // Calculate probability only at new nodes
	    if(nIter == 0 || i % 3 != 0)
	      {
		double r = dijet->PtMeas(0) / t;
		double p0 = PDF(r);
		r =  dijet->PtMeas(1) / t;
		double p1 = PDF(r);
		pp.push_back(p0 * p1); // Store product of pdfs
	      }
	    else
	      {
		pp.push_back(pp_old.at(i/3));       // Store product of pdfs previously calcluated
	      }
	  }

	// Sum up weighted function values
	for(unsigned int i = 0; i < pp.size(); i++)
	  {
	    double w = 1;                      // Weight w from Simpson's rule
	    if( i > 0 && i < (pp.size() - 1) ) // w = 1 for x0 and last node
	      {
		if( i % 3 == 0 )               // w = 2 for x3, x6, ...
		  {
		    w = 2;
		  }
		else
		  {
		    w = 3;
		  }
	      }
	    pint += w * (pp.at(i));            // Sum up weighted function values
	  }
	pint *= (3 * h / 8);                   // Apply overall normalization
	nIter++;
      }
    pint /= norm;

    if( pint < 0. || pint > 1. ) dijet->SetStatus(1);

    return -1.*log(pint);
  }



  //----------------------------------------------------------
  double Fitter::NLogLPhotonJet(PhotonJetEvent * photonjet) const
  {
    double t = photonjet->PtPhoton();
    double m = photonjet->PtMeas();
    double p = PDF(m/t);
    p /= t;

    if( p < 0. || p > 1. ) photonjet->SetStatus(1);

    return -1.*log(p);
  }



  //----------------------------------------------------------
  double Fitter::PDF(double x) const
  {
    double p = 0.;
    if( mModel == "FermiTail" )
      p = PDFFermiTail(x,mPar);
    else if( mModel == "TwoGauss" )
      p = PDFTwoGauss(x,mPar);
    else if( mModel == "ThreeGauss" )
      p = PDFThreeGauss(x,mPar);
    else if( mModel == "ExpTail" )
      p = PDFExpTail(x,mPar);
    else if( mModel == "Hist" )
      p = PDFHist(x,mPar);
    else if( mModel == "HistGauss" )
      p = PDFHistGauss(x,mPar);
    
    return p;
  }



  //----------------------------------------------------------
  double Fitter::PDFFermiTail(double x, const std::vector<double>& par) const
  {
    double c = par.at(0);
    if (c < 0) {
      c = 0;
    }
    if (c > 1) {
      c = 1;
    }
    double mu = 1.;
    double sigma = par.at(1);
    double T = par.at(2);
    if (T < 0) {
      T = 0;
    }

    double p = c / TMath::Sqrt(2* M_PI ) / sigma * exp(-pow((x - mu) / sigma, 2) / 2) + (1 - c)
      / (T * log(1 + exp(mu / T))) / (exp((x - mu) / T) + 1);
    return p;
  }



  // --------------------------------------------------
  double Fitter::FermiTail(double * x, double * par)
  {
    double c = par[0];
    if (c < 0) {
      c = 0;
    }
    if (c > 1) {
      c = 1;
    }
    double mu = 1.;
    double sigma = par[1];
    double T = par[2];
    if (T < 0) {
      T = 0;
    }

    double p = c / sqrt(2* M_PI ) / sigma * exp(-pow((x[0] - mu) / sigma, 2) / 2) + (1 - c)
      / (T * log(1 + exp(mu / T))) / (exp((x[0] - mu) / T) + 1);
    return p;
  }



  //----------------------------------------------------------
  double Fitter::PDFThreeGauss(double x, const std::vector<double>& par) const
  {
    double u0 = 1.;
    double s0 = par.at(0);  // Sigma of main Gaussian
    double c1 = par.at(1);  // Norm of 1. Gaussian tail
    double u1 = par.at(2);  // Mean of 1. Gaussian tail
    double s1 = par.at(3);  // Sigma of 1. Gaussian tail
    double c2 = par.at(4);  // Norm of 2. Gaussian tail
    double u2 = par.at(5);  // Mean of 2. Gaussian tail
    double s2 = par.at(6);  // Sigma of 2. Gaussian tail

    //Take care of proper normalization
    if(c1 < 0.)  c1 = 0.;
    if(c1 > 0.2) c1 = 0.2;
    if(c2 < 0.)  c2 = 0.;
    if(c2 > 0.2) c2 = 0.2;

    double p  = (1 - (c1+c2)) / sqrt(2* M_PI) / s0 * exp(-pow((x - u0) / s0, 2) / 2);
    p        +=            c1 / sqrt(2* M_PI) / s1 * exp(-pow((x - u1) / s1, 2) / 2);
    p        +=            c2 / sqrt(2* M_PI) / s2 * exp(-pow((x - u2) / s2, 2) / 2);

    if( p < 0 ) p = 1.E-10;
  
    return p;
  }



  //----------------------------------------------------------
  double Fitter::ThreeGauss(double *x, double *par)
  {
    double u0 = 1.;
    double s0 = par[0];  // Sigma of main Gaussian
    double c1 = par[1];  // Norm of 1. Gaussian tail
    double u1 = par[2];  // Mean of 1. Gaussian tail
    double s1 = par[3];  // Sigma of 1. Gaussian tail
    double c2 = par[4];  // Norm of 2. Gaussian tail
    double u2 = par[5];  // Mean of 2. Gaussian tail
    double s2 = par[6];  // Sigma of 2. Gaussian tail

    // Take care of proper normalization
    if(c1 < 0.) c1 = 0.;
    if(c1 > 0.2) c1 = 0.2;
    if(c2 < 0.) c2 = 0.;
    if(c2 > 0.2) c2 = 0.2;

    double p  = (1 - (c1+c2)) / sqrt(2* M_PI) / s0 * exp(-pow((x[0] - u0) / s0, 2) / 2);
    p        +=            c1 / sqrt(2* M_PI) / s1 * exp(-pow((x[0] - u1) / s1, 2) / 2);
    p        +=            c2 / sqrt(2* M_PI) / s2 * exp(-pow((x[0] - u2) / s2, 2) / 2);

    if( p < 0 ) p = 1.E-10;

    return p;
  }



  //----------------------------------------------------------
  double Fitter::PDFTwoGauss(double x, const std::vector<double>& par) const
  {
    double u0 = 1.;
    double s0 = par.at(0);  // Sigma of main Gaussian
    double c1 = par.at(1);  // Norm of Gaussian tail
    double u1 = par.at(2);  // Mean of Gaussian tail
    double s1 = par.at(3);  // Sigma of  Gaussian tail

    //Take care of proper normalization
    if(c1 < 0.) c1 = 0.;
    if(c1 > 1.) c1 = 1.;


    double p  = (1. - c1) / sqrt(2* M_PI) / s0 * exp(-pow((x - u0) / s0, 2) / 2);
    p        +=       c1 / sqrt(2* M_PI) / s1 * exp(-pow((x - u1) / s1, 2) / 2);

    return p;
  }



  //----------------------------------------------------------
  double Fitter::TwoGauss(double *x, double *par)
  {
    double u0 = 1.;
    double s0 = par[0];  // Sigma of main Gaussian
    double c1 = par[1];  // Norm of Gaussian tail
    double u1 = par[2];  // Mean of Gaussian tail
    double s1 = par[3];  // Sigma of Gaussian tail

    // Take care of proper normalization
    if(c1 < 0.) c1 = 0.;
    if(c1 > 1.) c1 = 1.;

    double p  = (1. - c1) / sqrt(2* M_PI) / s0 * exp(-pow((x[0] - u0) / s0, 2) / 2);
    p        +=       c1 / sqrt(2* M_PI) / s1 * exp(-pow((x[0] - u1) / s1, 2) / 2);

    return p;
  }



  //----------------------------------------------------------
  double Fitter::PDFExpTail(double x, const std::vector<double>& par) const
  {
    double c = par.at(0);   // Normalization
    double u = 1.;          // Mean of central Gaussian
    double s = par.at(1);   // Sigma of central Gaussian
    double k = par.at(2);   // Decay constant of exponential
    double m = par.at(3);   // Cutoff of exponential

    if(c < 0.) c = 0.;
    if(c > 1.) c = 1.;

    double p  = c / sqrt(2* M_PI) / s * exp(-pow((x - u) / s, 2) / 2.);
    if( x > m ) p += (1. - c) * ( exp( k*(m - x) ) ) / k;
	   
    return p;
  }



  // --------------------------------------------------
  double Fitter::ExpTail(double * x, double * par)
  {
    double c = par[0];   // Normalization
    double u = 1.;       // Mean of central Gaussian
    double s = par[1];   // Sigma of central Gaussian
    double k = par[2];   // Decay constant of exponential
    double m = par[3];   // Cutoff of exponential

    if (c < 0.) c = 0.;
    if (c > 1.) c = 1.;

    double p  = c / sqrt(2* M_PI) / s * exp(-pow((x[0] - u) / s, 2) / 2.);
    if( x[0] > m ) p += (1. - c) * ( exp( k*(m - x[0]) ) ) / k;

    return p;
  }



  // --------------------------------------------------
  double Fitter::PDFHist(double x, const std::vector<double>& par) const
  {
    // Set histogram bin content according
    // to parameters
    for(int i = 0; i < GetNPar(); i++)
      {
	double p = par.at(i);
	if( p < 0. ) p = 0.;
	mPDFHist->SetBinContent(1+i,p);
      }

    // Normalize histogrammed pdf
    mPDFHist->Scale( 1. / ( mPDFHist->Integral("width") ) );

    // Return probability density for given x
    double p = mPDFHist->Interpolate(x);
    // Interpolation into underflow/overflow bins
    // keeps value of last bin in range
    if( x < mPDFHistMin || x > mPDFHistMax ) p = 0.;

    return p;
  }



  // --------------------------------------------------
  double Fitter::PDFHistGauss(double x, const std::vector<double>& par) const
  {
    double c = par.at(0);
    double u = 1.;
    double s = par.at(1);

    if( c < 0. ) c = 0.;
    if( c > 1. ) c = 1.;

    // Set histogram bin content according
    // to parameters
    for(int i = 2; i < GetNPar(); i++)
      {
	double p = par.at(i);
	if( p < 0. ) p = 0.;
	mPDFHist->SetBinContent(i-1,p);
      }

    // Normalize histogrammed pdf
    mPDFHist->Scale( 1. / ( mPDFHist->Integral("width") ) );

    // Return probability density for given x
    double p = c / sqrt(2* M_PI ) / s * exp(-pow((x - u) / s, 2) / 2)
      //  + (1. - c) * (mPDFHist->Interpolate(x));
    + (1. - c) * (mPDFHist->GetBinContent(mPDFHist->FindBin(x)));

    return p;
  }



  // --------------------------------------------------
  TF1 * Fitter::GetTF1pdf() const
  {
    assert( mModel != "Hist" );

    TF1 * f = 0;
    if( mModel == "FermiTail" )
      {
	f = new TF1("fPDFFermiTail",&Fitter::FermiTail,0,6,3);
	for(int i = 0; i < GetNPar(mModel); i++)
	  {
	    f->SetParameter(i,mPar.at(i));
	    std::cout << "PAR " << i << "  " << mPar.at(i) << std::endl;
	  }
      }
    else if( mModel == "TwoGauss" )
      {
	f = new TF1("fPDFTwoGauss",&Fitter::TwoGauss,0,6,GetNPar(mModel));
	for(int i = 0; i < GetNPar(mModel); i++)
	  {
	    f->SetParameter(i,mPar.at(i));
	    std::cout << "PAR " << i << "  " << mPar.at(i) << std::endl;
	  }
      }
    else if( mModel == "ThreeGauss" )
      {
	f = new TF1("fPDFThreeGauss",&Fitter::ThreeGauss,0,6,GetNPar(mModel));
	for(int i = 0; i < GetNPar(mModel); i++)
	  {
	    f->SetParameter(i,mPar.at(i));
	    std::cout << "PAR " << i << "  " << mPar.at(i) << std::endl;
	  }
      }
    else if( mModel == "ExpTail" )
      {
	f = new TF1("fPDFExpTail",&Fitter::ExpTail,0,6,GetNPar(mModel));
	for(int i = 0; i < GetNPar(mModel); i++)
	  {
	    f->SetParameter(i,mPar.at(i));
	    std::cout << "PAR " << i << "  " << mPar.at(i) << std::endl;
	  }
      }
    return f;
  }



  // --------------------------------------------------
  TH1F * Fitter::GetTH1Fpdf() const
  {
    assert( mModel == "Hist" || mModel == "HistGauss" );

    TH1F * h = 0;

    if( mModel == "Hist" )
      {
	h = new TH1F("fPDFHist","",GetNPar(),mPDFHistMin,mPDFHistMax);
	for(int bin = 1; bin <= h->GetNbinsX(); bin++)
	  {
	    double r = h->GetBinCenter(bin);
	    h->SetBinContent(bin,PDFHist(r,mPar));
	  }
	double norm = h->Integral("width");
	h->Scale(1./norm);
      }
    else if( mModel == "HistGauss" )
      {
	h = new TH1F("fPDFHistGauss","",10*(GetNPar()-2),mPDFHistMin,mPDFHistMax);
	for(int bin = 1; bin <= h->GetNbinsX(); bin++)
	  {
	    double r = h->GetBinCenter(bin);
	    h->SetBinContent(bin,PDFHistGauss(r,mPar));
	  }
	double norm = h->Integral("width");
	h->Scale(1./norm);
      }

    return h;
  }



  // --------------------------------------------------
  int Fitter::GetNPar(const std::string& model) const
  {
    int n = 0;
    if( model == "FermiTail" )        n = 3;
    else if( model == "TwoGauss" )    n = 4;
    else if( model == "ThreeGauss" )  n = 7;
    else if( model == "ExpTail" )     n = 4;
    else if( model == "Hist" )        n = 15;
    else if( model == "HistGauss" )   n = 12;
    
    return n;
  }
}

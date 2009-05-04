#include "Fitter.h"

#include <cmath>
#include <iostream>

#include "NJetEvent.h"

namespace js
{
  void Fitter::Fit()
  {
    ////// LVMINI //////
    int NPAR =  GetNPar();
    int NITER = 100;
    int MVEC;
    NPAR < 29 ? MVEC = NPAR : MVEC = 29;
    int NAUX = 10000;
    
    std::vector<double> AUX(NAUX, 0.);
    double FSUM;
    double FOPT, FEDM, DUMMY;

//     std::vector<double> par;
//     par.push_back(0.9);
//     par.push_back(0.1);
//     par.push_back(0.05);


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
      FSUM = EvalFDiJetData(mPar,AUX);
      lvmfun_(&mPar.front(), FSUM, IRET, &AUX.front());
      lvmprt_(2, &AUX.front(), 2); //print out
    } while (IRET < 0 && ITER < NITER);
    int error_index = 2;
    error_index = lvmind_(error_index);
  }



  double Fitter::NLogLDiJetData(const std::vector<double>& par) const
  {
    int    maxNIter  = 5;       // Max number of iterations in interval splitting
    double eps = 1.e-5;

    // Normalize integral
    double norm = ( pow(mMax,3) - pow(mMin,3) ) / 3.;
    //std::cout << "norm " << logNorm << std::endl;

    double L = 0.;
    for(DataIt datait = mData.begin(); datait != mData.end(); datait++)
      {
	if( (*datait)->Type() == "DiJetEvent" )
	  {
	    const DiJetEvent * dijet = static_cast<const DiJetEvent*>(*datait);

	    double         pint      = 0.;      // Current value of integral over response pdfs
	    double         pint_old  = 1.;      // Value of integral over response pdfs from previous iteration
	    int            nIter     = 0;       // Current iteration in interval splitting
	    std::vector<double> pp;             // Product of current function values of response pdfs
	    std::vector<double> pp_old;         // Product of function values of response pdfs from previous iteration
	    double         h = mMax - mMin;       // Integration interval


	    // Iterate until precision or max. number iterations reached
	    while(fabs((pint - pint_old) / pint_old) > eps && nIter < maxNIter)
	      {
		//std::cout << "Iteration: " << nIter << std::endl;
		
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
			double p0 = Prob(r,par);
			r =  dijet->PtMeas(1) / t;
			double p1 = Prob(r,par);
			pp.push_back(p0 * p1); // Store product of pdfs
			//std::cout << "    pp.back(): " << pp.back() << std::endl;
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
	    L -= log(pint);
	  }
      }

    //std::cout << "NLogL is " << L << std::endl;

    return L;
  }



  //----------------------------------------------------------
  double Fitter::EvalFDiJetData(const std::vector<double>& par, std::vector<double>& grad) const
  {
    double FSUM = NLogLDiJetData(par);
    //std::cout << FSUM << std::flush;

    // calculate gradients
    double eps = 1.e-10;
    size_t npar = par.size();
    for (size_t i = 0; i < npar; i++)
      {
	std::vector<double> par1 = par;
	std::vector<double> par2 = par;
	par1.at(i) += eps;
	par2.at(i) -= eps;
	grad.at(i)  = (NLogLDiJetData(par1) - NLogLDiJetData(par2)) / 2. / eps;
	par1.at(i) += eps;
	par2.at(i) -= eps;
	grad.at(i + npar) = ( NLogLDiJetData(par1) - 2* FSUM + NLogLDiJetData(par2) ) / 4. / eps / eps;

	//	std::cout << "  " << grad.at(i) << std::flush;
      }

    std::cout << std::endl;

    return FSUM;
  }


  //----------------------------------------------------------
  double Fitter::Prob(double x, const std::vector<double>& par) const
  {
    double p = 0.;
    if( mModel == "FermiTail" )
      p = ProbFermiTail(x,par);
    else if( mModel == "TwoGauss" )
      p = ProbTwoGauss(x,par);
    else if( mModel == "ThreeGauss" )
      p = ProbThreeGauss(x,par);
    
    return p;
  }



  //----------------------------------------------------------
  double Fitter::ProbFermiTail(double x, const std::vector<double>& par) const
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

    double p = c / TMath::Sqrt(2* M_PI ) / sigma * exp(-pow((x[0] - mu) / sigma, 2) / 2) + (1 - c)
      / (T * log(1 + exp(mu / T))) / (exp((x[0] - mu) / T) + 1);
    return p;
  }



  //----------------------------------------------------------
  double Fitter::ProbThreeGauss(double x, const std::vector<double>& par) const
  {
    double u0 = 1.;
    double s0 = par.at(0);  // Sigma of main Gaussian
    double c1 = par.at(1);  // Norm of 1. Gaussian tail
    double c2 = par.at(2);  // Norm of 2. Gaussian tail
    double u1 = par.at(3);  // Mean of 1. Gaussian tail
    double s1 = par.at(4);  // Sigma of 1. Gaussian tail
    double u2 = par.at(5);  // Mean of 2. Gaussian tail
    double s2 = par.at(6);  // Sigma of 2. Gaussian tail

    //Take care of proper normalization
    if(c1 < 0.) c1 = 0.;
    if(c1 > 1.) c1 = 1.;
    if(u1 < 1.) u1 = 1.;
    if(c2 < 0.) c2 = 0.;
    if(c2 > 0.25) c2 = 0.25;
    if(u2 > 1.) u2 = 1.;


    double p  = (1 - (c1+c2)) / sqrt(2* M_PI) / s0 * exp(-pow((x - u0) / s0, 2) / 2);
    p        +=            c1 / sqrt(2* M_PI) / s1 * exp(-pow((x - u1) / s1, 2) / 2);
    p        +=            c2 / sqrt(2* M_PI) / s2 * exp(-pow((x - u2) / s2, 2) / 2);

    return p;
  }



  //----------------------------------------------------------
  double Fitter::ThreeGauss(double *x, double *par)
  {
    double u0 = 1.;
    double s0 = par[0];  // Sigma of main Gaussian
    double c1 = par[1];  // Norm of 1. Gaussian tail
    double c2 = par[2];  // Norm of 2. Gaussian tail
    double u1 = par[3];  // Mean of 1. Gaussian tail
    double s1 = par[4];  // Sigma of 1. Gaussian tail
    double u2 = par[5];  // Mean of 2. Gaussian tail
    double s2 = par[6];  // Sigma of 2. Gaussian tail

    // Take care of proper normalization
    if(c1 < 0.) c1 = 0.;
    if(c1 > 0.25) c1 = 0.25;
    if(u1 < 1.) u1 = 1.;
    if(c2 < 0.) c2 = 0.;
    if(c2 > 0.25) c2 = 0.25;

    double p  = (1 - (c1+c2)) / sqrt(2* M_PI) / s0 * exp(-pow((x[0] - u0) / s0, 2) / 2);
    p        +=            c1 / sqrt(2* M_PI) / s1 * exp(-pow((x[0] - u1) / s1, 2) / 2);
    p        +=            c2 / sqrt(2* M_PI) / s2 * exp(-pow((x[0] - u2) / s2, 2) / 2);

    return p;
  }



  //----------------------------------------------------------
  double Fitter::ProbTwoGauss(double x, const std::vector<double>& par) const
  {
    double u0 = 1.;
    double s0 = par.at(0);  // Sigma of main Gaussian
    double c1 = par.at(1);  // Norm of Gaussian tail
    double u1 = par.at(2);  // Mean of Gaussian tail
    double s1 = par.at(3);  // Sigma of  Gaussian tail

    //Take care of proper normalization
    if(c1 < 0.) c1 = 0.;
    if(c1 > 1.) c1 = 1.;


    double p  = (1 - c1) / sqrt(2* M_PI) / s0 * exp(-pow((x - u0) / s0, 2) / 2);
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

    double p  = (1 - c1) / sqrt(2* M_PI) / s0 * exp(-pow((x[0] - u0) / s0, 2) / 2);
    p        +=       c1 / sqrt(2* M_PI) / s1 * exp(-pow((x[0] - u1) / s1, 2) / 2);

    return p;
  }



  // --------------------------------------------------
  TF1 * Fitter::GetTF1() const
  {
    TF1 * f = 0;
    if( mModel == "FermiTail" )
      {
	f = new TF1("fPDFFermiTail",&Fitter::FermiTail,0,3,3);
	for(int i = 0; i < GetNPar(mModel); i++)
	  {
	    f->SetParameter(i,mPar.at(i));
	    std::cout << "PAR " << i << "  " << mPar.at(i) << std::endl;
	  }
      }
    else if( mModel == "TwoGauss" )
      {
	f = new TF1("fPDFTwoGauss",&Fitter::TwoGauss,0,3,GetNPar(mModel));
	for(int i = 0; i < GetNPar(mModel); i++)
	  {
	    f->SetParameter(i,mPar.at(i));
	    std::cout << "PAR " << i << "  " << mPar.at(i) << std::endl;
	  }
      }
    else if( mModel == "ThreeGauss" )
      {
	f = new TF1("fPDFThreeGauss",&Fitter::ThreeGauss,0,3,GetNPar(mModel));
	for(int i = 0; i < GetNPar(mModel); i++)
	  {
	    f->SetParameter(i,mPar.at(i));
	    std::cout << "PAR " << i << "  " << mPar.at(i) << std::endl;
	  }
      }

    return f;
  }



  // --------------------------------------------------
  int Fitter::GetNPar(const std::string& model) const
  {
    int n = 0;
    if( model == "FermiTail" ) n = 3;
    else if( model == "TwoGauss" ) n = 4;
    else if( model == "ThreeGauss" ) n = 7;
    
    return n;
  }
}

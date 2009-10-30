
// List of all modules not written in C/C++ language... 
extern "C" { 
  //Initialization of Limited-memory Variable-metrik Minimization
  void lvmini_( int &npar, int &mvec, int const &nfcn, double *aux);
   
  //test printout
  void lvmout_(int &n, int &m, double * aux);
   
  //print on unit LUN
  void lvmprt_(int const &lun, double *aux, int const &jarg);
   
  //pair update
  void lvmupd_(double &fac, double *svec, double *yvecp, double *yvecm, double *aux);

  //multiply, addition
  void dbaxpy_(int &n, double &da, double * dx, double * dy); 
 
  //matrix inversion  
  //          V = symmetric N-by-N matrix in symmetric storage mode
  //              V(1) = V11, V(2) = V12, V(3) = V22, V(4) = V13, . . .
  //              replaced by inverse matrix  
  //    DIAG(N) =  double precision scratch array
  //    NEXT(N) =  aux array (double precision instead of integer)
  void lvminv_(double * v, int &n, int &nrank, double &diag, double &next);

  //minimization of [Minuit] function 
  void lvmidi_(int &npar, double *par, void (*funct_)(int&,double*,double&,double*,int), double *aux);

  //Linesearch Algorithms wih Guaranteed Sufficient Decrease"
  //by JJ More and D Thuente, ACM Transactions on Mathenatical
  //Software 20 (1994), pp 286-307 - modified by V.B.
  void mtline_(int &n,double *x,double &f,double *g,double *s,double &stp,
               double *wa, double *dgout,int &info);
  //void cstepm_(...) //probably internal function, only used by mtline_
  
  //pointer and total mem calculation
  int lvmdim_(int &npar, int &mvec);

  //return index for ...
  //     IARG = 0   function value
  //          = 1   parameter vector
  //          = 2   approximate errors
  //          = 3   accurate errors
  //          = 4   global correlations
  //          = 5   covariance matrix  
  int lvmind_(int &iarg);
  
  //Dot product of two vectors: DX * DY
  double dbdot_(int &n, double *dx, double *dy);
  
  //change fit constants 
  void lvmeps_(float const &eps, float const &wlf1, float const &wlf2);
  
  //called with new F and gradient
  void lvmfun_(double *x,double &f, int &jret, double *aux);
}


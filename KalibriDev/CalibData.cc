//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: CalibData.cc,v 1.30 2009/06/11 17:34:45 mschrode Exp $
//
#include "CalibData.h"

#include <map>
#include <cmath>

unsigned int TAbstractData::total_n_pars = 0;
double (*TData::ScaleResidual)(double z2) = &TData::ScaleNone;
std::vector<TAbstractData*> TData_TruthMess::resultcache = std::vector<TAbstractData*>(1);

double TData_TruthMess::chi2_fast(double *temp_derivative1, double *temp_derivative2, double const epsilon) const
{ 
  double new_mess  = GetParametrizedMess();
  double new_error = GetParametrizedErr(&new_mess);
  double weight = GetWeight();
  double new_chi2 = weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(new_error*new_error) );

  double temp1 = 0.;		// Value of chi2 at par+epsilon
  double temp2 = 0.;		// Value of chi2 at par-epsilon
  
  double dmess_dp, derror_dp;   
  unsigned idx = _index*_n_par; //_index==bin; idx==bin*Free_parameters_per_bin
  for (unsigned i=idx; i<idx+_n_par; ++i){
    temp1 = 0.;
    temp2 = 0.;
    double oldpar = _par[i-idx];
    _par[i-idx]  += epsilon;
    dmess_dp  = GetParametrizedMess();
    derror_dp = GetParametrizedErr(&dmess_dp);
    temp2 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/(derror_dp*derror_dp) );


    _par[i-idx]  = oldpar - epsilon;
    dmess_dp  = GetParametrizedMess();
    derror_dp = GetParametrizedErr(&dmess_dp);
    temp1 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/(derror_dp*derror_dp) );

    // Difference of chi2 at par+epsilon and par-epsilon
    temp_derivative1[i] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i] += (temp2 + temp1 - 2.*new_chi2); // for 2nd derivative

    _par[i-idx]  = oldpar;
  }
  return new_chi2;
};


double TData_TruthMultMess::chi2_fast(double* temp_derivative1, double* temp_derivative2, double const epsilon) const {
  double sum_mess=0.0, sum_error2=0.0, new_error, new_error2, new_mess, new_chi2, dmess_dp, derror_dp;   
  double sm1[total_n_pars],sm2[total_n_pars],se1[total_n_pars],se2[total_n_pars];//store sum_m & sum_e2 for df/dp_i
  double weight = GetWeight();
  for (unsigned i=0; i<total_n_pars; ++i){
    sm1[i] = 0.0;
    sm2[i] = 0.0;
    se1[i] = 0.0;
    se2[i] = 0.0;
  }
  unsigned idx;
  double temp1=0.;
  double temp2=0.;

  if(trackuse){
    //idx = (*_vectrack.begin())->GetIndex()*(*_vectrack.begin())->GetNumberOfPars();
    idx = (*_vectrack.begin())->GetIndex();
    new_mess = GetParametrizedTrackMess();
    new_error = GetParametrizedErr();
    new_chi2  = weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/( new_error*new_error) );

    //Get Track Parameters
    for (unsigned i=idx; i<idx+(*_vectrack.begin())->GetNumberOfPars() ; ++i){
      std::vector<TAbstractData*>::const_iterator it=_vectrack.begin();
      double oldpar = (*it)->GetPar()[i-idx];
	  
      (*it)->GetPar()[i-idx]  += epsilon;
      dmess_dp  =GetParametrizedTrackMess();
      new_error = GetParametrizedErr();
      temp2 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/ (new_error * new_error) );
      
      (*it)->GetPar()[i-idx]  = oldpar - epsilon;
      dmess_dp  =GetParametrizedTrackMess();
      new_error = GetParametrizedErr();
      temp1 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/ (new_error * new_error ));
      // Difference of chi2 at par+epsilon and par-epsilon
      //int n_JetAndTowerPars = _n_par + (*_vecmess.begin())->GetNumberOfPars();
      temp_derivative1[i] += (temp2 - temp1); // for 1st derivative
      temp_derivative2[i] += (temp2 + temp1 - 2.*new_chi2); // for 2nd derivative
      (*it)->GetPar()[i-idx] = oldpar;      
    }
  }

  else{
    for (std::vector<TAbstractData*>::const_iterator it=_vecmess.begin();
	 it!=_vecmess.end(); ++it) {
      new_mess    = (*it)->GetParametrizedMess();	 
      sum_mess   += new_mess; 
      new_error   = (*it)->GetParametrizedErr(&new_mess);
      new_error2  = new_error * new_error;
      sum_error2 +=  new_error2;

      idx = (*it)->GetIndex()*(*it)->GetNumberOfPars();
      for (unsigned i=0; i<total_n_pars; ++i){
	if (i>=idx && i<idx+(*it)->GetNumberOfPars()) {   
	  double oldpar = (*it)->GetPar()[i-idx];
	  
	  (*it)->GetPar()[i-idx]  += epsilon;
	  dmess_dp  =(*it)->GetParametrizedMess();
	  sm2[i]   += dmess_dp;
	  new_error = (*it)->GetParametrizedErr(&dmess_dp);
	  se2[i]   += new_error * new_error;
	  
	  (*it)->GetPar()[i-idx]  = oldpar - epsilon;
	  dmess_dp  =(*it)->GetParametrizedMess();
	  sm1[i]   += dmess_dp;
	  new_error = (*it)->GetParametrizedErr(&dmess_dp);
	  se1[i]   += new_error * new_error;
	  (*it)->GetPar()[i-idx] = oldpar;
	  
	} else {
	  sm1[i] += new_mess;
	  sm2[i] += new_mess;
	  se1[i] += new_error2;
	  se2[i] += new_error2;
	}
      }
    } 
    new_mess  = GetParametrizedMess(&sum_mess);
    new_error = _err(&new_mess, _mess, _error);    //saves time not to use total GetParametrizedError(incl. tower & tracks)
    new_chi2  = weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error) );
    

    temp1 = 0.;		// Value of chi2 at par+epsilon
    temp2 = 0.;		// Value of chi2 at par-epsilon
    idx = _index; //@@to be fixed -> introduce a eta-phi binning for JES
    for (unsigned i=0; i<total_n_pars; ++i){
      if ((i>=idx && i<idx+_n_par)) continue;//considered below
      temp1 = 0.;
      temp2 = 0.;
      new_mess  = sm2[i];  //the measurement with modified parameter "i"
      dmess_dp  = GetParametrizedMess(&new_mess);
      derror_dp = _err(&dmess_dp, _mess, _error);
      temp2 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/(se2[i] + derror_dp*derror_dp) );
      
      // same for p_i-epsilon:
      new_mess  = sm1[i];  
      dmess_dp  = GetParametrizedMess(&new_mess);
      derror_dp = _err(&dmess_dp, _mess, _error);
      temp1 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/(se1[i] + derror_dp*derror_dp) );
      
      // Difference of chi2 at par+epsilon and par-epsilon
      temp_derivative1[i] += (temp2 - temp1); // for 1st derivative
      temp_derivative2[i] += (temp2 + temp1 - 2.*new_chi2); // for 2nd derivative
    }
  }
  
  idx = _index; //@@to be fixed -> introduce a eta-phi binning for JES	
  for (unsigned i=idx; i<idx+_n_par; ++i){
    temp1 = 0.;
    temp2 = 0.;
    //ok, we have to change the jet's parametrization:
    double oldpar =  _par[i-idx];
    _par[i-idx]  += epsilon;
    if(trackuse){
      dmess_dp  =GetParametrizedTrackMess();
      new_error = GetParametrizedErr();
      temp2 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/ (new_error * new_error ));
    }
    else{
      dmess_dp  = GetParametrizedMess(&sum_mess);
      derror_dp = _err(&dmess_dp, _mess, _error);
      temp2 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/(se2[i] + derror_dp*derror_dp) );
    }
    _par[i-idx]  = oldpar - epsilon;
    if(trackuse){
      dmess_dp  =GetParametrizedTrackMess();
      new_error = GetParametrizedErr();
      temp1 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/ (new_error * new_error ));
    }
    else{
      dmess_dp  = GetParametrizedMess(&sum_mess);  
      derror_dp = _err(&dmess_dp, _mess, _error);
      temp1 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/(se1[i] + derror_dp*derror_dp) );
    }
    // Difference of chi2 at par+epsilon and par-epsilon
    temp_derivative1[i] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i] += (temp2 + temp1 - 2.*new_chi2); // for 2nd derivative
    _par[i-idx]  = oldpar;
  }

  return new_chi2;
};

double TData_MessMess::chi2_fast(double * temp_derivative1, double*  temp_derivative2, 
 			         double const epsilon) const
{
  double new_chi2, new_mess, new_error, sum_error2=0.0;
  double weight = GetWeight();

  //Get all tower parameter used in this event:  
  std::map<int,double*> tpars;
  for (std::vector<TAbstractData*>::const_iterator it=_vecmess.begin();
       it!=_vecmess.end(); ++it) {
    for (unsigned i= 0 ; i < (*it)->GetNumberOfPars(); ++i) {
      tpars[ (*it)->GetIndex() * (*it)->GetNumberOfPars() + i ] = &((*it)->GetPar()[i]);
    }
  }
  for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
       mit!=_m2.end();++mit) {
    std::vector<TAbstractData*>::const_iterator mitend=(*mit)->GetRef().end();  
    for (std::vector<TAbstractData*>::const_iterator it=(*mit)->GetRef().begin();
	 it!=mitend; ++it) { 
      for (unsigned i= 0 ; i < (*it)->GetNumberOfPars(); ++i) {
	tpars[ (*it)->GetIndex() * (*it)->GetNumberOfPars() + i ] = &((*it)->GetPar()[i]);
      }
    }  
  }




  //Out of Cone can hardly be determined by MessMess (but Punch throughs can)


  

  //Add all jet parameters in this event:
  for (unsigned i=0; i<_n_par; ++i)
    tpars[ i+_index ]= &(_par[i]);
  for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
       mit!=_m2.end();++mit)
    for (unsigned i=0; i<(*mit)->GetNumberOfPars(); ++i) {
      tpars[ i+(*mit)->GetIndex() ]= &((*mit)->GetPar()[i]);
    }

  

  
  //This event's chi^2 for the current (unchanged) parameters:
  new_chi2 = chi2();

  //Calc. & Cache derivatives w.r.t. tower parameters:
  double temp1 = 0.;		//Value of chi2 at par+epsilon
  double temp2 = 0.;		//Value of chi2 at par-epsilon
  for (unsigned i=0; i<total_n_pars; ++i){
     // in case the ith parameter is used for this event, calculate 
     // derivative dchi/dpar[i] and cache it:
     if ( tpars.find(i)!=tpars.end() ) {
       double oldpar =  *tpars[i];
       *tpars[i] += epsilon;
       sum_error2 = 0.0;
       for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
            mit!=_m2.end(); ++mit) {
	 new_error   = (*mit)->GetParametrizedErr();
	 sum_error2 += new_error * new_error;
       }
       new_error = GetParametrizedErr();
       new_mess  = GetMessCombination();  
       temp2 = weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error) );

       sum_error2 = 0.0;

       *tpars[i] = oldpar - epsilon;
       for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
            mit!=_m2.end(); ++mit) {
  	  new_error   = (*mit)->GetParametrizedErr();
	  sum_error2 += new_error * new_error;
       }
       new_error = GetParametrizedErr();
       new_mess  = GetMessCombination();  
       temp1 = weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error) );

       // Difference of chi2 at par+epsilon and par-epsilon
       temp_derivative1[i] += (temp2 - temp1); // for 1st derivative
       temp_derivative2[i] += (temp2 + temp1 - 2.*new_chi2); // for 2nd derivative
       *tpars[i] = oldpar;
    }
  }
  //cout << "chi2 = " <<new_chi2<<", temp1="<<temp1<<", temp2="<<temp2<<endl;
  return new_chi2;
}





double TData_PtBalance::combine() const{
  double x, y, dummy = GetParametrizedMess();
  //here _direction[0,1] is a normalized vector in eta-phi plane in pT direction.
  
  x = dummy * _direction[0];
  y = dummy * _direction[1];
  for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
       it!=_m2.end();++it){
    dummy = (*it)->GetParametrizedMess(); 
    x += dummy * (*it)->GetDirection()[0];
    y += dummy * (*it)->GetDirection()[1];  
  }
  return sqrt(x*x+y*y);
};


double TData_PtBalance::chi2_fast(double * temp_derivative1, double*  temp_derivative2, 
 			         double const epsilon) const
{
  double new_chi2, new_mess;
  double weight = GetWeight();

  //Get all tower parameter used in this event:  
  std::map<int,double*> tpars;
  for (std::vector<TAbstractData*>::const_iterator it=_vecmess.begin();
       it!=_vecmess.end(); ++it) {
    for (unsigned i= 0 ; i < (*it)->GetNumberOfPars(); ++i) {
      tpars[ (*it)->GetIndex() * (*it)->GetNumberOfPars() + i ] = &((*it)->GetPar()[i]);
    }
  }
  for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
       mit!=_m2.end();++mit) {
    std::vector<TAbstractData*>::const_iterator mitend=(*mit)->GetRef().end();  
    for (std::vector<TAbstractData*>::const_iterator it=(*mit)->GetRef().begin();
	 it!=mitend; ++it) { 
      for (unsigned i= 0 ; i < (*it)->GetNumberOfPars(); ++i) {
	tpars[ (*it)->GetIndex() * (*it)->GetNumberOfPars() + i ] = &((*it)->GetPar()[i]);
      }
    }  
  }

  //Out of Cone can hardly be determined by Di-Jet (but Punch throughs can)


  

  //Add all jet parameters in this event:
  for (unsigned i=0; i<_n_par; ++i)
    tpars[ i+_index ]= &(_par[i]);
  for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
       mit!=_m2.end();++mit)
    for (unsigned i=0; i<(*mit)->GetNumberOfPars(); ++i) {
      tpars[ i+(*mit)->GetIndex() ]= &((*mit)->GetPar()[i]);
    }

  
     
    double parascale;
      
  //This event's chi^2 for the current (unchanged) parameters:
  new_chi2 = chi2();

  //Calc. & Cache derivatives w.r.t. tower parameters:
  double temp1 = 0.;		//Value of chi2 at par+epsilon
  double temp2 = 0.;		//Value of chi2 at par-epsilon
  for (unsigned i=0; i<total_n_pars; ++i){
     // in case the ith parameter is used for this event, calculate 
     // derivative dchi/dpar[i] and cache it:
     if ( tpars.find(i)!=tpars.end() ) {
       double oldpar =  *tpars[i];
       *tpars[i] += epsilon;
       parascale = 0;
       for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin(); mit!=_m2.end(); ++mit) {
	 new_mess    = (*mit)->GetParametrizedMess();	//aufraeumen, wenn richtig 
	 parascale += new_mess;
       }
       new_mess  = GetParametrizedMess();
       parascale += new_mess;
       parascale /= 2.;
    //should be changed to ptsum (scalar) of the part projected to leading jet axis devided by 2 for n>2 n-jets 
       new_mess  = combine() / parascale;
       temp2 = weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(_error * _error));

       parascale = 0.0;

       *tpars[i] = oldpar - epsilon;
       for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin(); mit!=_m2.end(); ++mit) {
	  new_mess    = (*mit)->GetParametrizedMess();	 	 
	  parascale += new_mess;
       }
       new_mess  = GetParametrizedMess();
       parascale += new_mess;
       parascale /= 2;
    //should be changed to ptsum (scalar) of the part projected to leading jet axis devided by 2 for n>2 n-jets 
       new_mess  = combine() / parascale;
       temp1 = weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(_error * _error) );

       // Difference of chi2 at par+epsilon and par-epsilon
       temp_derivative1[i] += (temp2 - temp1); // for 1st derivative
       temp_derivative2[i] += (temp2 + temp1 - 2.*new_chi2); // for 2nd derivative

       //if((fabs(temp1-new_chi2) > 50* weight) || (fabs(temp2 - new_chi2) > 50*weight))
       // cout<<new_chi2<<"  T2: "<<temp2<<"    T1: "<<temp1<<"   W: "<<weight<<"      sum: "<<sum<<"  PS: "<<parascale<<"   C: "<<combine()<<endl;

       *tpars[i] = oldpar;
    }
  }

  return new_chi2;
}

double TData_InvMass2::combine() const{
  double x, y, z=0., e, tx, ty, tz, dummy = GetParametrizedMess();
  //here _direction is three dimensional. dimension[0,1] is a normalized vector
  //in eta-phi plane in pT direction. _direction[2] is the p_z component of 
  //the measurement.
  
  x = dummy * _direction[0];
  y = dummy * _direction[1];
  if (GetMess()->pt!=0.) z = dummy/GetMess()->pt * _direction[2];
  e = sqrt( x*x + y*y + z*z );
  for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
       it!=_m2.end();++it){
    dummy = (*it)->GetParametrizedMess();
    tx = dummy * (*it)->GetDirection()[0];
    ty = dummy * (*it)->GetDirection()[1];  
    if ((*it)->GetMess()->pt!=0.) tz = dummy/(*it)->GetMess()->pt * (*it)->GetDirection()[2];  
    else tz = 0.;
    x += tx;
    y += ty;
    z += tz;
    e += sqrt( tx*tx + ty*ty + tz*tz );
  }
  //cout << "m = "<<sqrt(e*e - x*x - y*y - z*z)<<endl<<endl;
  return sqrt(e*e - x*x - y*y - z*z);     
};

std::vector<TAbstractData*> TData_ParLimit::_cache = std::vector<TAbstractData*>(1);

double TData_ParLimit::chi2_fast(double* temp_derivative1, 
				 double* temp_derivative2, double const epsilon) const {
  // Penalty term with current parameter values
  double new_chi2  = chi2();

  // Variation of parameters
  double oldpar    = _par[0];
  _par[0]         += epsilon;
  double temp2     = chi2();
  _par[0]          = oldpar - epsilon;
  double temp1     = chi2();

  // Difference of chi2 at par+epsilon and par-epsilon
  temp_derivative1[_index] += (temp2 - temp1);                // for 1st derivative
  temp_derivative2[_index] += (temp2 + temp1 - 2.*new_chi2);  // for 2nd derivative

  // Reset original parameter value
  _par[0]  = oldpar;

  return new_chi2;
}
    

//!  Scale the normalized, squared residual
//!  \f$ z^{2} = \chi^{2}/\textrm{weight} \f$
//!  using the Cauchy-Function
//!  \f$ z^{2} \rightarrow c^{2}\ln( 1 + (z/c)^{2} ) \f$
//!
//!  \param z2 Normalized and squared residual
//!  \return Scaled residual
double TData::ScaleCauchy(double const z2)
{
  double const c = 2.3849;
  return (c*c) * log( 1 + z2*(1.0/(c*c)) );
}

//!  Scale the normalized, squared residual
//! \f$ z^{2} = \chi^{2}/\textrm{weight} \f$
//!  using the Huber-Function
//!  \f[
//!  z^{2} \rightarrow 
//!  \left\{
//!     \begin{array}{ll}
//!        z & \textrm{for } |z| <= c \\ c ( 2|z| - c ) & \textrm{for } |z| > c
//!      \end{array}
//!    \right.
//!  \f]
//!
//!  \param z2 Normalized and squared residual
//!  \return Scaled residual
double TData::ScaleHuber(double const z2)
{
  static double const c = 1.345;
  double const z = sqrt(z2);
  return (  std::abs(z) <= c  ?  z2  :  c*(2.*std::abs(z) - c)  );
}

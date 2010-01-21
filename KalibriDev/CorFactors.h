//!  \brief   Container class for jet correction factors
//
//    $Id: CorFactors.h,v 1.2 2009/11/26 10:27:48 stadie Exp $
//   
#ifndef CORFACTORS_H
#define CORFACTORS_H

class CorFactors
{
 public :
  CorFactors(double L1=1.0, double L2=1.0, double L3=1.0, double L4=1.0, 
	     double L5=1.0, double JPT=1.0, double JPTL2L3=1.0) 
    : l1_(L1),l2_(L2),l3_(L3),l4_(L4),l5_(L5),jpt_(JPT),jptL2L3_(JPTL2L3) {};
  double getL1()  const { return l1_; }    //!< Return L1 correction factor (zero-suppression)
  double getL2()  const { return l2_; }    //!< Return L2 correction factor (relative, in eta)
  double getL3()  const { return l3_; }    //!< Return L3 correction factor (absolute, in pt)
  double getL4()  const { return l4_; }    //!< Return L4 correction factor (electromagnetic fraction)
  double getL5()  const { return l5_; }    //!< Return L5 correction factor (flavor)
  double getJPT() const { return jpt_; }   //!< Return Jet+Track correction factor
  double getL2L3() const { return l2_*l3_; }   //!< Return product of L2 and L3 correction factors
  double getJPTL2L3() const { return jptL2L3_; }   //!< Return product of L2 and L3 correction factors for Jet+Track
  double getToL2() const { return l1_*l2_; }         //!< Return factor needed to get L2 corrected from raw jets: L1*L2
  double getToL3() const { return getToL2()*l3_; }   //!< Return factor needed to get L3 corrected from raw jets: L1*L2*L3
  double getToL4() const { return getToL3()*l4_; }   //!< Return factor needed to get L4 corrected from raw jets: L1*L2*L3*L4
  double getToL5() const { return getToL4()*l5_; }   //!< Return factor needed to get L5 corrected from raw jets: L1*L2*L3*L4*L5
  double getToJPTL3() const { return jpt_*l1_*jptL2L3_; }   //!< Return factor needed to get L3 corrected from raw jets for JPT: JPT*L1*JPTL2L3
 private :
  double l1_;      //!< Level 1 correction factor (zero-suppression)
  double l2_;      //!< Level 2 correction factor (relative, in eta)
  double l3_;      //!< Level 3 correction factor (absolute, in pt)
  double l4_;      //!< Level 4 correction factor (electromagnetic fraction)
  double l5_;      //!< Level 5 correction factor (flavor)
  double jpt_;     //!< Jet+Track correction factor
  double jptL2L3_; //!< Product of level 2 and level 3 correction factors for Jet+Track
}; 

#endif

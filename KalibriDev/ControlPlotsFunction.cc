// $Id: ControlPlotsFunction.cc,v 1.5 2010/04/01 16:29:17 stadie Exp $

#include "ControlPlotsFunction.h"

#include "CorFactors.h"
#include "Jet.h"
#include "JetTruthEvent.h"



//!  Will return 0 if \p type does not exist in \p types_
// ----------------------------------------------------------------   
double ControlPlotsFunction::yValue(const Event * evt, ControlPlotsConfig::CorrectionType type) const {
  double yValue = 0.;
  std::map<ControlPlotsConfig::CorrectionType,Function>::const_iterator it = yFuncs_.find(type);
  if( it != yFuncs_.end() ) yValue = (this->*(it->second))(evt); 

  return yValue;
}



// ----------------------------------------------------------------   
void ControlPlotsFunction::addYFunction(ControlPlotsConfig::CorrectionType type, Function func) {
  // Store correction types in a vector for easy access
  types_.push_back(type);
  // Store correction function
  std::map<ControlPlotsConfig::CorrectionType,Function>::const_iterator it = yFuncs_.find(type);
  if( it == yFuncs_.end() ) {
    yFuncs_[type] = func;
  } else {
    std::cerr << "WARNING: Adding function for already existing CorrectionType '" << type << "'\n";
  }
}



//!  \brief Returns #eta of the jet
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetEta(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->jet()->eta();
}


//!  \brief Returns p_{T} of the jet
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetPt(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->jet()->pt();
}

//!  \brief Returns ECal fraction of the jet
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetEMF(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->jet()->EmEt() / (jte->jet()->EmEt() + jte->jet()->HadEt());
}


//!  \brief Returns the #phi #phi moment of the jet
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetMomentPhiPhi(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->jet()->momentPhiPhi();
}


//!  \brief Returns the #eta #eta moment of the jet
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetMomentEtaEta(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->jet()->momentEtaEta();
}


//!  \brief Returns the scaled #eta #eta moment of the jet
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetScaledEtaEta(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->jet()->scaledEtaWidth();
}


//!  \brief Returns the scaled #phi #phi moment of the jet
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetScaledPhiPhi(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->jet()->scaledPhiWidth();
}


//!  \brief Returns truth pt
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventTruthPt(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->GetTruth();
}

//!  \brief Returns flavor of jet
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetFlavor(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->jet()->flavor();
}


//!  \brief Returns the jet response
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  The response is defined as
//!  \f[  p^{jet}_{T} / p^{true}_{T}\f].
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventResponse(const Event * evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  Jet * jet = static_cast<Jet*>(jte->GetMess());

  return jet->pt() / jte->GetTruth();
}



//!  \brief Returns the corrected jet response
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  The response is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the Kalibri JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventResponseKalibriCorrected(const Event * evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  Jet * jet = static_cast<Jet*>(jte->GetMess());

  return jet->correctedEt() / jte->GetTruth();
}



//!  \brief Returns the corrected jet response
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  The response is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3 JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventResponseL2L3Corrected(const Event * evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  Jet * jet = static_cast<Jet*>(jte->GetMess());

  return jet->corFactors().getL2L3() * jet->pt() / jte->GetTruth();
}



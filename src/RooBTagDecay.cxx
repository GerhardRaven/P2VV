/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl           *
 *   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl           *
 *                                                                           *
 * Copyright (c) 2011, Nikhef. All rights reserved.                          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/


//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Description of a B decay time distribution with flavour tags, effects
// of CP violation, mixing and life time differences. This function can 
// be analytically convolved with any RooResolutionModel implementation.
// END_HTML
//


#include "RooFit.h"
#include "RooMsgService.h"

#include "TMath.h"
#include "RooBTagDecay.h"
#include "RooRealVar.h"
#include "RooRandom.h"

ClassImp(RooBTagDecay);

//_____________________________________________________________________________
RooBTagDecay::RooBTagDecay(const char *name, const char* title, 
    RooRealVar& time, RooAbsReal& tau, RooAbsReal& dGamma, RooAbsReal& dm,
    RooAbsReal& coshCoef, RooAbsReal& sinhCoef, RooAbsReal& cosCoef,
    RooAbsReal& sinCoef, const RooResolutionModel& model, DecayType type,
    Bool_t checkVars) :
  RooAbsAnaConvPdf(name, title, model, time),
  _time("time", "proper time", this, time),
  _iTag("iTag", "initial state tag", this),
  _fTag("fTag", "final state tag", this),
  _tau("tau", "average decay time", this, tau),
  _dGamma("dGamma", "Delta Gamma", this, dGamma),
  _dm("dm", "Delta mass", this, dm),
  _dilution("dilution", "mis-tag dilution", this),
  _ADilWTag("ADilWTag", "dilution/wrong tag asymmetry", this),
  _ANorm("ANorm", "normalization asymmetry", this),
  _avgCEven("avgCEven", "CP average even coefficients", this),
  _avgCOdd("avgCOdd", "CP average odd coefficients", this),
  _coshCoef("coshCoef", "cosh coefficient", this, coshCoef),
  _sinhCoef("sinhCoef", "sinh coefficient", this, sinhCoef),
  _cosCoef("cosCoef", "cos coefficient", this, cosCoef),
  _sinCoef("sinCoef", "sin coefficient", this, sinCoef),
  _type(type),
  _tags(0),
  _checkVars(checkVars)
{
  // constructor without flavour tags (behaves like RooBDecay)
  if (!checkVarDep(time, kTRUE)) assert(0);

  declareBases();
}

//_____________________________________________________________________________
RooBTagDecay::RooBTagDecay(const char *name, const char* title,
    RooRealVar& time, RooAbsCategory& iTag, RooAbsCategory& fTag,
    RooAbsReal& tau, RooAbsReal& dGamma, RooAbsReal& dm, RooAbsReal& dilution,
    RooAbsReal& ADilWTag, RooAbsReal& ANorm, RooAbsReal& avgCEven,
    RooAbsReal& avgCOdd, RooAbsReal& cosCoef, const RooResolutionModel& model,
    DecayType type, Bool_t checkVars) :
  RooAbsAnaConvPdf(name, title, model, time),
  _time("time", "proper time", this, time),
  _iTag("iTag", "initial state tag", this, iTag),
  _fTag("fTag", "final state tag", this, fTag),
  _tau("tau", "average decay time", this, tau),
  _dGamma("dGamma", "Delta Gamma", this, dGamma),
  _dm("dm", "Delta mass", this, dm),
  _dilution("dilution", "mis-tag dilution", this, dilution),
  _ADilWTag("ADilWTag", "dilution/wrong tag asymmetry", this, ADilWTag),
  _ANorm("ANorm", "normalization asymmetry", this, ANorm),
  _avgCEven("avgCEven", "CP average even coefficients", this, avgCEven),
  _avgCOdd("avgCOdd", "CP average odd coefficients", this, avgCOdd),
  _coshCoef("coshCoef", "cosh coefficient", this, cosCoef),
  _sinhCoef("sinhCoef", "sinh coefficient", this),
  _cosCoef("cosCoef", "cos coefficient", this, cosCoef),
  _sinCoef("sinCoef", "sin coefficient", this),
  _type(type),
  _tags(3),
  _checkVars(checkVars)
{
  // constructor with both initial and final state flavour tags
  // (decay into flavour specific final state)

  if (!checkVarDep(time, kTRUE)) assert(0);
  if (!checkVarDep(iTag, kTRUE)) assert(0);
  if (!checkVarDep(fTag, kTRUE)) assert(0);
  if (!checkTag(iTag, kTRUE)) assert(0);
  if (!checkTag(fTag, kTRUE)) assert(0);

  declareBases();
}

//_____________________________________________________________________________
RooBTagDecay::RooBTagDecay(const char *name, const char* title,
    RooRealVar& time, RooAbsCategory& iTag, RooAbsReal& tau,
    RooAbsReal& dGamma, RooAbsReal& dm, RooAbsReal& dilution,
    RooAbsReal& ADilWTag, RooAbsReal& avgCEven, RooAbsReal& avgCOdd,
    RooAbsReal& coshCoef, RooAbsReal& sinhCoef, RooAbsReal& cosCoef,
    RooAbsReal& sinCoef, const RooResolutionModel& model, DecayType type,
    Bool_t checkVars) :
  RooAbsAnaConvPdf(name, title, model, time),
  _time("time", "proper time", this, time),
  _iTag("iTag", "initial state tag", this, iTag),
  _fTag("fTag", "final state tag", this),
  _tau("tau", "average decay time", this, tau),
  _dGamma("dGamma", "Delta Gamma", this, dGamma),
  _dm("dm", "Delta mass", this, dm),
  _dilution("dilution", "mis-tag dilution", this, dilution),
  _ADilWTag("ADilWTag", "dilution/wrong tag asymmetry", this, ADilWTag),
  _ANorm("ANorm", "normalization asymmetry", this),
  _avgCEven("avgCEven", "CP average even coefficients", this, avgCEven),
  _avgCOdd("avgCOdd", "CP average odd coefficients", this, avgCOdd),
  _coshCoef("coshCoef", "cosh coefficient", this, coshCoef),
  _sinhCoef("sinhCoef", "sinh coefficient", this, sinhCoef),
  _cosCoef("cosCoef", "cos coefficient", this, cosCoef),
  _sinCoef("sinCoef", "sin coefficient", this, sinCoef),
  _type(type),
  _tags(1),
  _checkVars(checkVars)
{
  // constructor with only an initial state flavour tag
  // (decay into CP self-conjugate state)

  if (!checkVarDep(time, kTRUE)) assert(0);
  if (!checkVarDep(iTag, kTRUE)) assert(0);
  if (!checkTag(iTag, kTRUE)) assert(0);

  declareBases();
}

//_____________________________________________________________________________
RooBTagDecay::RooBTagDecay(const RooBTagDecay& other, const char* name) :
  RooAbsAnaConvPdf(other, name),
  _time("time", this, other._time),
  _iTag("iTag", this, other._iTag),
  _fTag("fTag", this, other._fTag),
  _tau("tau", this, other._tau),
  _dGamma("dGamma", this, other._dGamma),
  _dm("dm", this, other._dm),
  _dilution("dilution", this, other._dilution),
  _ADilWTag("ADilWTag", this, other._ADilWTag),
  _ANorm("ANorm", this, other._ANorm),
  _avgCEven("avgCEven", this, other._avgCEven),
  _avgCOdd("avgCOdd", this, other._avgCOdd),
  _coshCoef("coshCoef", this, other._coshCoef),
  _sinhCoef("sinhCoef", this, other._sinhCoef),
  _cosCoef("cosCoef", this, other._cosCoef),
  _sinCoef("sinCoef", this, other._sinCoef),
  _coshBasis(other._coshBasis),
  _sinhBasis(other._sinhBasis),
  _cosBasis(other._cosBasis),
  _sinBasis(other._sinBasis),
  _type(other._type),
  _tags(other._tags),
  _checkVars(other._checkVars)
{
  // copy constructor
}

//_____________________________________________________________________________
RooBTagDecay::~RooBTagDecay()
{
  // destructor
}

//_____________________________________________________________________________
Double_t RooBTagDecay::coefficient(Int_t basisIndex) const
{
  // return the value of the specified basis function's coefficient

  // get coefficient's value
  Double_t coefVal = 0.;
  if (basisIndex == _coshBasis) {
    coefVal = _coshCoef;
  } else if (basisIndex == _sinhBasis) {
    if (_tags > 1) return 0.;
    else coefVal = _sinhCoef;
  } else if (basisIndex == _cosBasis) {
    coefVal = _cosCoef;
  } else if (basisIndex == _sinBasis) {
    if (_tags > 1) return 0.;
    else coefVal = _sinCoef;
  }

  // return the coefficient if there are no explicit tags
  if (_tags == 0) return coefVal;

  // terms that are even in initial state tag
  if (basisIndex == _coshBasis || basisIndex == _sinhBasis) {
    if (_tags == 1)
      return (_avgCEven + _iTag * _dilution
          * (_avgCOdd - _avgCEven * _ADilWTag)) * coefVal;
    else
      return (1. - _fTag * _ANorm) * (_avgCEven + _iTag * _dilution
          * (_avgCOdd - _avgCEven * _ADilWTag)) * coefVal;
  }

  // terms that are odd in initial state tag
  if (basisIndex == _cosBasis || basisIndex == _sinBasis) {
    if (_tags == 1)
      return (_avgCOdd + _iTag * _dilution
          * (_avgCEven - _avgCOdd * _ADilWTag)) * coefVal;
    else
      return (_fTag - _ANorm) * (_avgCOdd + _iTag * _dilution
          * (_avgCEven - _avgCOdd * _ADilWTag)) * coefVal;
  }

  return 0.;
}

//_____________________________________________________________________________
RooArgSet* RooBTagDecay::coefVars(Int_t basisIndex) const 
{
  // Return the set of variables used exclusively by the coefficient of the
  // specified basis function. The caller of this function is responsible for
  // deleting the returned argset.

  RooArgSet* coefVars = 0;

  // get variables from specified coefficient
  if (basisIndex == _coshBasis)
    coefVars = _coshCoef.arg().getParameters((RooArgSet*)0);
  if (basisIndex == _sinhBasis && _tags < 2)
    coefVars = _sinhCoef.arg().getParameters((RooArgSet*)0);
  if (basisIndex == _cosBasis)
    coefVars = _cosCoef.arg().getParameters((RooArgSet*)0);
  if (basisIndex == _sinBasis && _tags < 2)
    coefVars = _sinCoef.arg().getParameters((RooArgSet*)0);

  if (coefVars == 0) {
    coefVars = new RooArgSet("parameters");
  } else {
    // add tagging variables
    if (_tags > 0) {
      RooArgSet* tempSet = _iTag.arg().getParameters((RooArgSet*)0);
      if (tempSet->getSize() > 0)
        coefVars->add(*tempSet, kTRUE);
      else
        coefVars->add(_iTag.arg());
      delete tempSet;

      tempSet = _dilution.arg().getParameters((RooArgSet*)0);
      if (tempSet->getSize() > 0)
        coefVars->add(*tempSet, kTRUE);
      else
        coefVars->add(_dilution.arg());
      delete tempSet;

      tempSet = _ADilWTag.arg().getParameters((RooArgSet*)0);
      if (tempSet->getSize() > 0)
        coefVars->add(*tempSet, kTRUE);
      else
        coefVars->add(_ADilWTag.arg());
      delete tempSet;

      tempSet = _avgCEven.arg().getParameters((RooArgSet*)0);
      if (tempSet->getSize() > 0)
        coefVars->add(*tempSet, kTRUE);
      else
        coefVars->add(_avgCEven.arg());
      delete tempSet;

      tempSet = _avgCOdd.arg().getParameters((RooArgSet*)0);
      if (tempSet->getSize() > 0)
        coefVars->add(*tempSet, kTRUE);
      else
        coefVars->add(_avgCOdd.arg());
      delete tempSet;

      if (_tags > 1) {
        tempSet = _fTag.arg().getParameters((RooArgSet*)0);
        if (tempSet->getSize() > 0)
          coefVars->add(*tempSet, kTRUE);
        else
          coefVars->add(_fTag.arg());
        delete tempSet;

        tempSet = _ANorm.arg().getParameters((RooArgSet*)0);
        if (tempSet->getSize() > 0)
          coefVars->add(*tempSet, kTRUE);
        else
          coefVars->add(_ANorm.arg());
        delete tempSet;
      }
    }

    // remove variables that are in basis functions
    TIterator* coefVarsIter = coefVars->createIterator();
    RooAbsArg* coefArg;
    while ((coefArg = (RooAbsArg*)coefVarsIter->Next()) != 0) {
      for (Int_t basisIter = 0; basisIter < _convSet.getSize(); ++basisIter) {
        if (_convSet.at(basisIter)->dependsOn(*coefArg))
          coefVars->remove(*coefArg, kTRUE);
      }
    }
    delete coefVarsIter;
  }

  return coefVars;
}

//_____________________________________________________________________________
Int_t RooBTagDecay::getCoefAnalyticalIntegral(Int_t coef, RooArgSet& allVars,
    RooArgSet& analVars, const char* rangeName) const 
{
  // copy variables that can be integrated over analytically from allVars to
  // analVars for the specified basis function's coefficient and return
  // integration code

  // integrate numerically if variables are unchecked
  if (!_checkVars || !checkVarDep(_time.arg())
      || (_tags > 0 && (!checkVarDep(_iTag.arg()) || !checkTag(_iTag.arg())))
      || (_tags > 1 && (!checkVarDep(_iTag.arg()) || !checkTag(_fTag.arg()))))
    return 0;

  // get integration code
  Int_t  intCode = 0;
  Bool_t intITag = kFALSE;
  Bool_t intFTag = kFALSE;
  RooArgSet intVars = allVars;
  if (_tags > 0) intITag = intVars.remove(_iTag.arg(), kTRUE, kTRUE);
  if (_tags > 1) intFTag = intVars.remove(_fTag.arg(), kTRUE, kTRUE);

  if (intVars.getSize() > 0) {
    if (coef == _coshBasis) {
      intCode = 4 * _coshCoef.arg().getAnalyticalIntegral(intVars, analVars,
          rangeName);
    } else if (coef == _sinhBasis) {
      if (_tags > 1)
        return 0;
      else
        intCode = 4 * _sinhCoef.arg().getAnalyticalIntegral(intVars, analVars,
            rangeName);
    } else if (coef == _cosBasis) {
      intCode = 4 * _cosCoef.arg().getAnalyticalIntegral(intVars, analVars,
          rangeName);
    } else if (coef == _sinBasis) {
      if (_tags > 1)
        return 0;
      else
        intCode = 4 * _sinCoef.arg().getAnalyticalIntegral(intVars, analVars,
            rangeName);
    }
  }

  // return the integration code if there are no explicit tags
  if (_tags == 0) return intCode;

  // check that tagging parameters don't depend on integration variables
  TIterator* analVarIter = analVars.createIterator();
  RooAbsArg* analVar = 0;
  while ((analVar = (RooAbsArg*)analVarIter->Next()) != 0) {
    if (_dilution.arg().dependsOn(*analVar)
        || _ADilWTag.arg().dependsOn(*analVar)
        || (_tags > 1 && _ANorm.arg().dependsOn(*analVar))
        || _avgCEven.arg().dependsOn(*analVar)
        || _avgCOdd.arg().dependsOn(*analVar)) {
      analVars.removeAll();
      return 0;
    }
  }
  delete analVarIter;

  // add the initial state tag to integration variables
  if (intITag) {
    intCode += 1;
    analVars.add(_iTag.arg());
  }

  // add the final state tag to integration variables
  if (intFTag) {
    intCode += 2;
    analVars.add(_fTag.arg());
  }

  return intCode;
}

//_____________________________________________________________________________
Double_t RooBTagDecay::coefAnalyticalIntegral(Int_t coef, Int_t code,
    const char* rangeName) const 
{
  // return analytical integral for basis function's coefficient

  // get integration code for coefficient
  Int_t  coefCode = TMath::FloorNint(code / 4.);
  Bool_t intITag = (Bool_t)(code & 1);
  Bool_t intFTag = (Bool_t)(code & 2);

  // get coefficient's integral
  Double_t coefInt = 0.;
  if (coef == _coshBasis) {
    if (coefCode != 0)
      coefInt = _coshCoef.arg().analyticalIntegral(coefCode, rangeName);
    else
      coefInt = _coshCoef.arg().getVal();
  } else if (coef == _sinhBasis) {
    if (_tags > 1)
      return 0.;
    else if (coefCode != 0)
      coefInt = _sinhCoef.arg().analyticalIntegral(coefCode, rangeName);
    else
      coefInt = _sinhCoef.arg().getVal();
  } else if (coef == _cosBasis) {
    if (coefCode != 0)
      coefInt = _cosCoef.arg().analyticalIntegral(coefCode, rangeName);
    else
      coefInt = _cosCoef.arg().getVal();
  } else if (coef == _sinBasis) {
    if (_tags > 1)
      return 0.;
    else if (coefCode != 0)
      coefInt = _sinCoef.arg().analyticalIntegral(coefCode, rangeName);
    else
      coefInt = _sinCoef.arg().getVal();
  }

  // return the integral if there are no explicit tags
  if (_tags == 0) return coefInt;

  // terms that are even in initial state tag
  if (coef == _coshBasis || coef == _sinhBasis) {
    if (_tags == 1) {
      if (intITag)
        return 2. * _avgCEven * coefInt;
      else
        return (_avgCEven + _iTag * _dilution
            * (_avgCOdd - _avgCEven * _ADilWTag)) * coefInt;
    } else {
      if (intITag && intFTag)
        return 4. * _avgCEven * coefInt;
      else if (intITag)
        return 2. * (1. - _fTag * _ANorm) * _avgCEven * coefInt;
      else if (intFTag)
        return 2. * (_avgCEven + _iTag * _dilution
            * (_avgCOdd - _avgCEven * _ADilWTag)) * coefInt;
      else
        return (1. - _fTag * _ANorm) * (_avgCEven + _iTag * _dilution
            * (_avgCOdd - _avgCEven * _ADilWTag)) * coefInt;
    }
  }

  // terms that are odd in initial state tag
  if (coef == _cosBasis || coef == _sinBasis) {
    if (_tags == 1) {
      if (intITag)
        return 2. * _avgCOdd * coefInt;
      else
        return (_avgCOdd + _iTag * _dilution
            * (_avgCEven - _avgCOdd * _ADilWTag)) * coefInt;
    } else {
      if (intITag && intFTag)
        return -4. * _ANorm * _avgCOdd * coefInt;
      else if (intITag)
        return 2. * (_fTag - _ANorm) * _avgCOdd * coefInt;
      else if (intFTag)
        return -2. * _ANorm * (_avgCOdd + _iTag * _dilution
            * (_avgCEven - _avgCOdd * _ADilWTag)) * coefInt;
      else
        return (_fTag - _ANorm) * (_avgCOdd + _iTag * _dilution
            * (_avgCEven - _avgCOdd * _ADilWTag)) * coefInt;
    }
  }

  return 0.;
}

//_____________________________________________________________________________
Int_t RooBTagDecay::getGenerator(const RooArgSet& directVars,
    RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  // copy variables that can be generated directly by RooBTagDecay from
  // directVars to generateVars and return generation code

  // use accept/reject if the flavour tags are unchecked
  if (!_checkVars || !checkVarDep(_time.arg())
      || (_tags > 0 && (!checkVarDep(_iTag.arg()) || !checkTag(_iTag.arg())))
      || (_tags > 1 && (!checkVarDep(_iTag.arg()) || !checkTag(_fTag.arg()))))
    return 0;

  Int_t genCode = 0;

  // search for time variable
  RooAbsArg* arg = directVars.find(_time.arg().GetName());
  if (arg == 0) return genCode;

  genCode += 1;
  generateVars.add(*arg);

  if (_tags == 0) return genCode;

  // search for initial state tag variable
  arg = directVars.find(_iTag.arg().GetName());
  if (arg == 0) return genCode;

  if (_tags == 1) {
    genCode += 2;
    generateVars.add(*arg);
  } else {
    // search for final state tag variable
    RooAbsArg* arg1 = directVars.find(_fTag.arg().GetName());
    if (arg1 == 0) return genCode;

    genCode += 2 + 4;
    generateVars.add(*arg);
    generateVars.add(*arg1);
  }

  return genCode;
}

//_____________________________________________________________________________
void RooBTagDecay::generateEvent(Int_t code)
{
  // generate values for the variables corresponding to the generation code

  // check generation code
  if (code == 7) {
    if (_tags != 3) {
      coutE(InputArguments) << "RooBTagDecay::generateEvent(" << GetName()
          << ") generation code 7 requires a PDF with both an initial state and a final state explicit flavour tag"
          << endl;
      assert(0);
    }
  } else if (code == 3) {
    if (_tags != 1) {
      coutE(InputArguments) << "RooBTagDecay::generateEvent(" << GetName()
          << ") generation code 3 requires a PDF with (only) an explicit initial state flavour tag"
          << endl;
      assert(0);
    }
  } else if (code != 1) {
    coutE(InputArguments) << "RooBTagDecay::generateEvent(" << GetName()
        << ") unrecognized generation code: " << code << endl;
    assert(0);
  }

  // set parameters
  Double_t gammaMin  = 1. / _tau - TMath::Abs(_dGamma) / 2.;

  // generate event
  while(true) {
    // generate time variable with the exponential envelope function
    Double_t timeGen = -log(RooRandom::uniform()) / gammaMin;
    if (_type == Flipped || (_type == DoubleSided
        && RooRandom::uniform() > 0.5))
      timeGen *= -1.;
    if (timeGen < _time.min() || timeGen > _time.max()) continue;

    // get (unnormalized) PDF and envelope values for generated time
    Double_t tAbs         = TMath::Abs(timeGen);
    Double_t pdfVal       = exp(-tAbs / _tau);
    Double_t envVal       = exp(-gammaMin * tAbs);
    Double_t evenTerms    = _coshCoef * cosh(_dGamma * timeGen / 2.);
    Double_t oddTerms     = _cosCoef * cos(_dm * timeGen);
    Double_t evenTermsEnv = 0.;
    Double_t oddTermsEnv  = 0.;
    if (_tags < 2) {
        evenTerms    += _sinhCoef * sinh(_dGamma * timeGen / 2.);
        oddTerms     += _sinCoef * sin(_dm * timeGen);
        evenTermsEnv  = TMath::Abs(_coshCoef) + TMath::Abs(_sinhCoef) / 2.;
        oddTermsEnv   = sqrt(_cosCoef * _cosCoef + _sinCoef * _sinCoef);
    } else {
        evenTermsEnv  = TMath::Abs(_coshCoef);
        oddTermsEnv   = TMath::Abs(_cosCoef);
    }

    if (code == 1) {
      // generate a value for the time variable only
      if (_tags == 0) {
        // use PDF without flavour tags
        pdfVal *= oddTerms + evenTerms;
        envVal *= evenTermsEnv + oddTermsEnv;

      } else {
        // use PDF with flavour tag(s)
        Double_t cEven = _avgCEven + _iTag * _dilution
            * (_avgCOdd - _avgCEven * _ADilWTag);
        Double_t cOdd = _avgCOdd + _iTag * _dilution
            * (_avgCEven - _avgCOdd * _ADilWTag);
        if (_tags > 1) {
          cEven *= 1. - _fTag * _ANorm;
          cOdd  *= _fTag - _ANorm;
        }

        pdfVal *= cEven * evenTerms + cOdd * oddTerms;
        envVal *= TMath::Abs(cEven) * evenTermsEnv
            + TMath::Abs(cOdd) * oddTermsEnv;
      }

    } else {
      // generate both time and flavour tag(s): use PDF integrated over tag(s)
      if (_tags == 1) {
        pdfVal *= _avgCEven * evenTerms + _avgCOdd * oddTerms;
        envVal *= TMath::Abs(_avgCEven) * evenTermsEnv
            + TMath::Abs(_avgCOdd) * oddTermsEnv;
      } else {
        pdfVal *= _avgCEven * evenTerms - _ANorm * _avgCOdd * oddTerms;
        envVal *= TMath::Abs(_avgCEven) * evenTermsEnv
            + TMath::Abs(_ANorm * _avgCOdd) * oddTermsEnv;
      }
    }

    // check if PDF value is positive
    if (pdfVal < 0.) {
      coutE(Generation) << "RooBTagDecay::generateEvent(" << GetName()
          << ") PDF is negative for generated time value" << endl;
      assert(0);
    }

    // accept/reject generated time value, using the envelope as maximum
    if (pdfVal < envVal * RooRandom::uniform()) continue;
    _time = timeGen;

    if (code != 1) {
      // generate value for initial state tag
      Int_t iTagGen = 1;
      Double_t ACP = 0.;
      Double_t cEvenDil = _dilution * (_avgCOdd - _avgCEven * _ADilWTag);
      Double_t cOddDil  = _dilution * (_avgCEven - _avgCOdd * _ADilWTag);
      if (_tags == 1)
        ACP += (cEvenDil * evenTerms + cOddDil * oddTerms)
            / (_avgCEven * evenTerms + _avgCOdd * oddTerms);
      else
        ACP += (cEvenDil * evenTerms - _ANorm * cOddDil * oddTerms)
            / (_avgCEven * evenTerms - _ANorm * _avgCOdd * oddTerms);

      if (2. * RooRandom::uniform() > 1. + ACP) iTagGen = -1;

      _iTag = iTagGen;

      if (_tags > 1) {
        Int_t fTagGen = 1;
        Double_t cEven = _avgCEven + iTagGen * cEvenDil;
        Double_t cOdd  = _avgCOdd  + iTagGen * cOddDil;
        Double_t Af = (-_ANorm * cEven * evenTerms + cOdd * oddTerms)
            / (cEven * evenTerms - _ANorm * cOdd * oddTerms);

        if (2. * RooRandom::uniform() > 1. + Af) fTagGen = -1;

        _fTag = fTagGen;
      }
    }

    // exit generation loop
    break;
  }
}

//_____________________________________________________________________________
Bool_t RooBTagDecay::checkVarDep(const RooAbsArg& var, Bool_t warn) const
{
  if (_checkVars && (_tau.arg().dependsOn(var)
      || _dGamma.arg().dependsOn(var) || _dm.arg().dependsOn(var)
      || (_tags > 0 && (_dilution.arg().dependsOn(var)
      || _ADilWTag.arg().dependsOn(var)
      || (_tags > 1 && _ANorm.arg().dependsOn(var))
      || _avgCEven.arg().dependsOn(var) || _avgCOdd.arg().dependsOn(var)))
      || _coshCoef.arg().dependsOn(var) || _cosCoef.arg().dependsOn(var)
      || (_tags < 2 && (_sinhCoef.arg().dependsOn(var)
      || _sinCoef.arg().dependsOn(var))))) {
    coutE(InputArguments) << "RooBTagDecay::checkVarDep(" << GetName()
        << ") parameters depend on " << var.GetName()
        << ": use \"checkVars = false\" if you insist on using this dependence"
        << endl;
    return kFALSE;
  } else if (warn && !_checkVars) {
    coutW(InputArguments) << "RooBTagDecay::checkVarDep(" << GetName()
        << ") parameters dependence on "
        << var.GetName()
        << " is unchecked (use of \"checkVars = false\" is not recommended): integrals will be calculated numerically and event generation will use accept/reject"
        << endl;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t RooBTagDecay::checkTag(const RooAbsCategory& tag, Bool_t warn) const
{
  if (_checkVars && (tag.numTypes() != 2 || !tag.isValidIndex(+1)
      || !tag.isValidIndex(-1))) {
    coutE(InputArguments) << "RooBTagDecay::checkTag(" << GetName()
        << ") flavour tag " << tag.GetName()
        << " can only have values +1 and -1: use \"checkVars = false\" if you insist on using other/additional values"
        << endl;
    return kFALSE;
  } else if (warn && !_checkVars) {
    coutW(InputArguments) << "RooBTagDecay::checkTag(" << GetName()
        << ") unchecked flavour tag "
        << tag.GetName()
        << " (use of \"checkVars = false\" is not recommended): integrals will be calculated numerically and event generation will use accept/reject"
        << endl;
  }

  return kTRUE;
}

//_____________________________________________________________________________
void RooBTagDecay::declareBases()
{
  // create basis functions for time variable
  if (_type == SingleSided) {
    // create basis functions for positive time values
    _coshBasis = declareBasis("exp(-@0/@1)*cosh(@0*@2/2)",
        RooArgList(_tau.arg(), _dGamma.arg()));
    _cosBasis = declareBasis("exp(-@0/@1)*cos(@0*@2)",
        RooArgList(_tau.arg(), _dm.arg()));

    if (_tags < 2) {
      _sinhBasis = declareBasis("exp(-@0/@1)*sinh(@0*@2/2)",
          RooArgList(_tau.arg(), _dGamma.arg()));
      _sinBasis = declareBasis("exp(-@0/@1)*sin(@0*@2)",
          RooArgList(_tau.arg(), _dm.arg()));
    }

  } else if (_type == Flipped) {
    // create basis functions for negative time values
    _coshBasis = declareBasis("exp(@0/@1)*cosh(@0*@2/2)",
        RooArgList(_tau.arg(), _dGamma.arg()));
    _cosBasis = declareBasis("exp(@0/@1)*cos(@0*@2)",
        RooArgList(_tau.arg(), _dm.arg()));

    if (_tags < 2) {
      _sinhBasis = declareBasis("exp(@0/@1)*sinh(@0*@2/2)",
          RooArgList(_tau.arg(), _dGamma.arg()));
      _sinBasis = declareBasis("exp(@0/@1)*sin(@0*@2)",
          RooArgList(_tau.arg(), _dm.arg()));
    }

  } else if (_type == DoubleSided) {
    // create basis functions for both positive and negative time values
    _coshBasis = declareBasis("exp(-TMath::Abs(@0)/@1)*cosh(@0*@2/2)",
        RooArgList(_tau.arg(), _dGamma.arg()));
    _cosBasis = declareBasis("exp(-TMath::Abs(@0)/@1)*cos(@0*@2)",
        RooArgList(_tau.arg(), _dm.arg()));

    if (_tags < 2) {
      _sinhBasis = declareBasis("exp(-TMath::Abs(@0)/@1)*sinh(@0*@2/2)",
          RooArgList(_tau.arg(), _dGamma.arg()));
      _sinBasis = declareBasis("exp(-TMath::Abs(@0)/@1)*sin(@0*@2)",
          RooArgList(_tau.arg(), _dm.arg()));
    }
  }
}


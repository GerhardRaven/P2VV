/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   JvL, Jeroen van Leerdam, Nikhef,           j.van.leerdam@nikhef.nl      *
 *   GR,  Gerhard Raven,      Nikhef & VU,      Gerhard.Raven@nikhef.nl      *
 *   PL,  Parker C Lund,      UC Irvine                                      *
 *   DK,  David Kirkby,       UC Irvine,        dkirkby@uci.edu              *
 *   WV,  Wouter Verkerke,    UC Santa Barbara, verkerke@slac.stanford.edu   *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
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
    RooAbsReal& sinCoef, const RooResolutionModel& model, DecayType type) :
  RooAbsAnaConvPdf(name, title, model, time),
  _time("time", "proper time", this, time),
  _iTag("iTag", "initial state tag", this),
  _fTag("fTag", "final state tag", this),
  _tau("tau", "average decay time", this, tau),
  _dGamma("dGamma", "Delta Gamma", this, dGamma),
  _dm("dm", "Delta mass", this, dm),
  _dilution("dilution", "mis-tag dilution", this),
  _ANorm("ANorm", "normalization asymmetry", this),
  _ATagEff("ATagEff", "tagging efficiency asymmetry", this),
  _ADilMisTag("ADilMisTag", "dilution/mis-tag asymmetry", this),
  _coshCoef("coshCoef", "cosh coefficient", this, coshCoef),
  _sinhCoef("sinhCoef", "sinh coefficient", this, sinhCoef),
  _cosCoef("cosCoef", "cos coefficient", this, cosCoef),
  _sinCoef("sinCoef", "sin coefficient", this, sinCoef),
  _type(type),
  _tags(0)
{
  // constructor without flavour tags (behaves like RooBDecay)

  declareBases();
}

//_____________________________________________________________________________
RooBTagDecay::RooBTagDecay(const char *name, const char* title,
    RooRealVar& time, RooAbsCategory& iTag, RooAbsReal& tau,
    RooAbsReal& dGamma, RooAbsReal& dm, RooAbsReal& dilution,
    RooAbsReal& ANorm, RooAbsReal& ATagEff, RooAbsReal& ADilMisTag,
    RooAbsReal& coshCoef, RooAbsReal& sinhCoef, RooAbsReal& cosCoef,
    RooAbsReal& sinCoef, const RooResolutionModel& model, DecayType type) :
  RooAbsAnaConvPdf(name, title, model, time),
  _time("time", "proper time", this, time),
  _iTag("iTag", "initial state tag", this, iTag),
  _fTag("fTag", "final state tag", this),
  _tau("tau", "average decay time", this, tau),
  _dGamma("dGamma", "Delta Gamma", this, dGamma),
  _dm("dm", "Delta mass", this, dm),
  _dilution("dilution", "mis-tag dilution", this, dilution),
  _ANorm("ANorm", "normalization asymmetry", this, ANorm),
  _ATagEff("ATagEff", "tagging efficiency asymmetry", this, ATagEff),
  _ADilMisTag("ADilMisTag", "dilution/mis-tag asymmetry", this, ADilMisTag),
  _coshCoef("coshCoef", "cosh coefficient", this, coshCoef),
  _sinhCoef("sinhCoef", "sinh coefficient", this, sinhCoef),
  _cosCoef("cosCoef", "cos coefficient", this, cosCoef),
  _sinCoef("sinCoef", "sin coefficient", this, sinCoef),
  _type(type),
  _tags(1)
{
  // constructor with initial state flavour tag (decay into CP eigenstate)

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
  _ANorm("ANorm", this, other._ANorm),
  _ATagEff("ATagEff", this, other._ATagEff),
  _ADilMisTag("ADilMisTag", this, other._ADilMisTag),
  _coshCoef("coshCoef", this, other._coshCoef),
  _sinhCoef("sinhCoef", this, other._sinhCoef),
  _cosCoef("cosCoef", this, other._cosCoef),
  _sinCoef("sinCoef", this, other._sinCoef),
  _coshBasis(other._coshBasis),
  _sinhBasis(other._sinhBasis),
  _cosBasis(other._cosBasis),
  _sinBasis(other._sinBasis),
  _type(other._type),
  _tags(other._tags)
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
  // return the coefficient of the specified basis function

  if (_tags == 0) {
    if (basisIndex == _coshBasis) return _coshCoef;
    if (basisIndex == _sinhBasis) return _sinhCoef;
    if (basisIndex == _cosBasis)  return _cosCoef;
    if (basisIndex == _sinBasis)  return _sinCoef;

  } else {
    Double_t cA = 1. + _ANorm * _ATagEff;
    Double_t cB = _ANorm + _ATagEff;

    if (basisIndex == _coshBasis) 
      return (cA + _iTag * _dilution * (cB - cA * _ADilMisTag)) * _coshCoef;
    if (basisIndex == _sinhBasis)
      return (cA + _iTag * _dilution * (cB - cA * _ADilMisTag)) * _sinhCoef;
    if (basisIndex == _cosBasis)
      return (cB + _iTag * _dilution * (cA - cB * _ADilMisTag)) * _cosCoef;
    if (basisIndex == _sinBasis)
      return (cB + _iTag * _dilution * (cA - cB * _ADilMisTag)) * _sinCoef;
  }

  return 0.;
}

//_____________________________________________________________________________
RooArgSet* RooBTagDecay::coefVars(Int_t basisIndex) const 
{
  // return the set of variables used exclusively by the coefficient of the
  // specified basis function

  RooArgSet* coefVars = 0;

  // get variables from specified coefficient
  if (basisIndex == _coshBasis)
    coefVars = _coshCoef.arg().getParameters((RooArgSet*)0);
  if (basisIndex == _sinhBasis)
    coefVars = _sinhCoef.arg().getParameters((RooArgSet*)0);
  if (basisIndex == _cosBasis)
    coefVars = _cosCoef.arg().getParameters((RooArgSet*)0);
  if (basisIndex == _sinBasis)
    coefVars = _sinCoef.arg().getParameters((RooArgSet*)0);

  // add tagging variables
  if (_tags > 0) {
    coefVars->add(_iTag.arg());
    coefVars->add(_dilution.arg());
    coefVars->add(_ANorm.arg());
    coefVars->add(_ATagEff.arg());
    coefVars->add(_ADilMisTag.arg());
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

  return coefVars;
}

//_____________________________________________________________________________
Int_t RooBTagDecay::getCoefAnalyticalIntegral(Int_t coef, RooArgSet& allVars,
    RooArgSet& analVars, const char* rangeName) const 
{
  // copy variables that can be integrated over analytically from allVars to
  // analVars for the specified basis function's coefficient and return
  // integration code

  // get coefficient
  const RooAbsReal* coefArg = 0;
  if (coef == _coshBasis) coefArg = &_coshCoef.arg();
  if (coef == _sinhBasis) coefArg = &_sinhCoef.arg();
  if (coef == _cosBasis)  coefArg = &_cosCoef.arg();
  if (coef == _sinBasis)  coefArg = &_sinCoef.arg();

  // get code from coefficient if it doesn't depend on tagging variables
  if (_tags == 0 || !(coefArg->dependsOn(_iTag.arg())
      || coefArg->dependsOn(_dilution.arg())
      || coefArg->dependsOn(_ANorm.arg())
      || coefArg->dependsOn(_ATagEff.arg())
      || coefArg->dependsOn(_ADilMisTag.arg())))
    return coefArg->getAnalyticalIntegral(allVars, analVars, rangeName);

  return 0;
}

//_____________________________________________________________________________
Double_t RooBTagDecay::coefAnalyticalIntegral(Int_t coef, Int_t code,
    const char* rangeName) const 
{
  // return analytical integral for basis function's coefficient

  if (_tags == 0) {
    if (coef == _coshBasis)
      return _coshCoef.arg().analyticalIntegral(code, rangeName);
    if (coef == _sinhBasis)
      return _sinhCoef.arg().analyticalIntegral(code, rangeName);
    if (coef == _cosBasis)
      return _cosCoef.arg().analyticalIntegral(code, rangeName);
    if (coef == _sinBasis)
      return _sinCoef.arg().analyticalIntegral(code, rangeName);

  } else {
    Double_t cA = 1. + _ANorm * _ATagEff;
    Double_t cB = _ANorm + _ATagEff;

    if (coef == _coshBasis) 
      return (cA + _iTag * _dilution * (cB - cA * _ADilMisTag))
          * _coshCoef.arg().analyticalIntegral(code, rangeName);
    if (coef == _sinhBasis)
      return (cA + _iTag * _dilution * (cB - cA * _ADilMisTag))
          * _sinhCoef.arg().analyticalIntegral(code, rangeName);
    if (coef == _cosBasis)
      return (cB + _iTag * _dilution * (cA - cB * _ADilMisTag))
          * _cosCoef.arg().analyticalIntegral(code, rangeName);
    if (coef == _sinBasis)
      return (cB + _iTag * _dilution * (cA - cB * _ADilMisTag))
          * _sinCoef.arg().analyticalIntegral(code, rangeName);
  }

  return 0.;
}

//_____________________________________________________________________________
Int_t RooBTagDecay::getGenerator(const RooArgSet& directVars,
    RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  // copy variables that can be generated directly by RooBTagDecay from
  // directVars to generateVars and return generation code

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

  genCode += 2;
  generateVars.add(*arg);

  if (_tags == 1) return genCode;

  // search for final state tag variable
  arg = directVars.find(_fTag.arg().GetName());
  if (arg == 0) return genCode;

  genCode += 4;
  generateVars.add(*arg);

  return genCode;
}

//_____________________________________________________________________________
void RooBTagDecay::generateEvent(Int_t code)
{
  // generate values for the variables corresponding to the generation code

  // check generation code
  if (code == 3) {
    if (_tags < 1) {
      coutE(InputArguments) << "RooBTagDecay::generateEvent(" << GetName()
          << ") generation code 3 requires a PDF with explicit tagging variables: abort" << endl;
      assert(0);
    }
  } else if (code != 1) {
    coutE(InputArguments) << "RooBTagDecay::generateEvent(" << GetName()
        << ") unrecognized generation code: " << code << endl;
    assert(0);
  }

  // check if Delta Gamma is positive
  if (_dGamma < 0.) {
    coutE(Generation) << "RooBTagDecay::generateEvent(" << GetName()
        << ") Delta Gamma has a negative value: abort" << endl;
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
    Double_t tAbs     = TMath::Abs(timeGen);
    Double_t pdfVal   = exp(-tAbs / _tau);
    Double_t envVal   = exp(-gammaMin * tAbs);
    Double_t cEvenAvg = 1.;
    Double_t cOddAvg  = 0.;
    Double_t evenTerms = _coshCoef * cosh(_dGamma * timeGen / 2.)
        + _sinhCoef * sinh(_dGamma * timeGen / 2.);
    Double_t oddTerms  = _cosCoef * cos(_dm * timeGen)
        + _sinCoef * sin(_dm * timeGen);

    if (code == 1) {
      // generate a value for the time variable only
      if (_tags == 0) {
        // use PDF without tags
        pdfVal *= oddTerms + evenTerms;
        envVal *= TMath::Abs(_coshCoef) + TMath::Abs(_sinhCoef) / 2.
            + sqrt(_cosCoef * _cosCoef + _sinCoef * _sinCoef);

      } else if (_tags == 1) {
        // use PDF with initial state tag
        cEvenAvg += _ANorm * _ATagEff;
        cOddAvg  += _ANorm + _ATagEff;
        Double_t cEven = cEvenAvg
            + _iTag * _dilution * (cOddAvg - cEvenAvg * _ADilMisTag);
        Double_t cOdd = cOddAvg
            + _iTag * _dilution * (cEvenAvg - cOddAvg * _ADilMisTag);

        pdfVal *= cEven * evenTerms + cOdd * oddTerms;
        envVal *= TMath::Abs(cEven)
            * (TMath::Abs(_coshCoef) + TMath::Abs(_sinhCoef) / 2.)
            + TMath::Abs(cOdd)
            * sqrt(_cosCoef * _cosCoef + _sinCoef * _sinCoef);
      }

    } else {
      // generate both time and initial state tag: use CP averaged PDF
      cEvenAvg += _ANorm * _ATagEff;
      cOddAvg  += _ANorm + _ATagEff;
      pdfVal *= cEvenAvg * evenTerms + cOddAvg * oddTerms;
      envVal *= TMath::Abs(cEvenAvg)
          * (TMath::Abs(_coshCoef) + TMath::Abs(_sinhCoef) / 2.)
          + TMath::Abs(cOddAvg)
          * sqrt(_cosCoef * _cosCoef + _sinCoef * _sinCoef);
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
      Double_t ACP = _dilution
          * ((cOddAvg - cEvenAvg * _ADilMisTag) * evenTerms
          + (cEvenAvg - cOddAvg * _ADilMisTag) * oddTerms)
          / (cEvenAvg * evenTerms + cOddAvg * oddTerms);

      if (2. * RooRandom::uniform() > 1. + ACP) iTagGen = -1;

      _iTag = iTagGen;
    }

    // exit generation loop
    break;
  }
}

//_____________________________________________________________________________
void RooBTagDecay::declareBases()
{
  // create basis functions for time variable
  if (_type == SingleSided) {
    // create basis functions for positive time values
    _coshBasis = declareBasis("exp(-@0/@1)*cosh(@0*@2/2)",
        RooArgList(_tau.arg(), _dGamma.arg()));
    _sinhBasis = declareBasis("exp(-@0/@1)*sinh(@0*@2/2)",
        RooArgList(_tau.arg(), _dGamma.arg()));
    _cosBasis = declareBasis("exp(-@0/@1)*cos(@0*@2)",
        RooArgList(_tau.arg(), _dm.arg()));
    _sinBasis = declareBasis("exp(-@0/@1)*sin(@0*@2)",
        RooArgList(_tau.arg(), _dm.arg()));
  } else if (_type == Flipped) {
    // create basis functions for negative time values
    _coshBasis = declareBasis("exp(@0/@1)*cosh(@0*@2/2)",
        RooArgList(_tau.arg(), _dGamma.arg()));
    _sinhBasis = declareBasis("exp(@0/@1)*sinh(@0*@2/2)",
        RooArgList(_tau.arg(), _dGamma.arg()));
    _cosBasis = declareBasis("exp(@0/@1)*cos(@0*@2)",
        RooArgList(_tau.arg(), _dm.arg()));
    _sinBasis = declareBasis("exp(@0/@1)*sin(@0*@2)",
        RooArgList(_tau.arg(), _dm.arg()));
  } else if (_type == DoubleSided) {
    // create basis functions for both positive and negative time values
    _coshBasis = declareBasis("exp(-TMath::Abs(@0)/@1)*cosh(@0*@2/2)",
        RooArgList(_tau.arg(), _dGamma.arg()));
    _sinhBasis = declareBasis("exp(-TMath::Abs(@0)/@1)*sinh(@0*@2/2)",
        RooArgList(_tau.arg(), _dGamma.arg()));
    _cosBasis = declareBasis("exp(-TMath::Abs(@0)/@1)*cos(@0*@2)",
        RooArgList(_tau.arg(), _dm.arg()));
    _sinBasis = declareBasis("exp(-TMath::Abs(@0)/@1)*sin(@0*@2)",
        RooArgList(_tau.arg(), _dm.arg()));
  }
}


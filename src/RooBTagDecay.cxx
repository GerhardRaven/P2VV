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
// Description of a B decay time distribution with flavour tagging, effects
// of CP violation, mixing and life time differences. This function can 
// be analytically convolved with any RooResolutionModel implementation.
// END_HTML
//


#include <memory>

#include "Riostream.h"
#include "RooFit.h"
#include "RooMsgService.h"

#include "TMath.h"
#include "RooBTagDecay.h"
#include "RooCategory.h"
#include "RooRandom.h"
#include "RooRealVar.h"

ClassImp(RooBTagDecay);

//_____________________________________________________________________________
RooBTagDecay::RooBTagDecay(const char *name, const char* title, 
    RooRealVar& time, RooAbsReal& tau, RooAbsReal& dGamma, RooAbsReal& dm,
    RooAbsReal& coshCoef, RooAbsReal& sinhCoef, RooAbsReal& cosCoef,
    RooAbsReal& sinCoef, const RooResolutionModel& model, DecayType type,
    Bool_t checkVars) :
  RooAbsAnaConvPdf(name, title, model, time),
  _time("time", "B lifetime", this, time),
  _iTag("iTag", "initial state tag", this),
  _fTag("fTag", "final state tag", this),
  _tagCat("tagCat", "tagging category", this),
  _tau("tau", "mean B lifetime", this, tau),
  _dGamma("dGamma", "Delta Gamma", this, dGamma),
  _dm("dm", "Delta mass", this, dm),
  _dilutions("dilutions", 0, this),
  _ADilWTags("ADilWTags", 0, this),
  _ANorm("ANorm", "normalization asymmetry", this),
  _avgCEvenSum("avgCEvenSum", "weighted sum of average even coeffs", this),
  _avgCOddSum("avgCOddSum", "weighted sum of average odd coeffs", this),
  _avgCEvens("avgCEvens", 0, this),
  _avgCOdds("avgCOdds", 0, this),
  _tagCatCoefs("tagCatCoefs", 0, this),
  _coshCoef("coshCoef", "cosh coefficient", this, coshCoef),
  _sinhCoef("sinhCoef", "sinh coefficient", this, sinhCoef),
  _cosCoef("cosCoef", "cos coefficient", this, cosCoef),
  _sinCoef("sinCoef", "sin coefficient", this, sinCoef),
  _createdVars("createdVars", 0, this),
  _decayType(type),
  _tagCatType(0),
  _tags(0),
  _iTagVal(2),
  _fTagVal(2),
  _checkVars(checkVars)
{
  // constructor without flavour tags (behaves like RooBDecay)
  if (!checkVarDep(time, kTRUE)) assert(0);

  declareBases();
}

//_____________________________________________________________________________
RooBTagDecay::RooBTagDecay(const char *name, const char* title,
    RooRealVar& time, RooCategory& iTag, RooCategory& fTag,
    RooAbsReal& tau, RooAbsReal& dGamma, RooAbsReal& dm, RooAbsReal& dilution,
    RooAbsReal& ADilWTag, RooAbsReal& ANorm, RooAbsReal& avgCEven,
    RooAbsReal& avgCOdd, RooAbsReal& cosCoef, const RooResolutionModel& model,
    DecayType type, Bool_t checkVars) :
  RooAbsAnaConvPdf(name, title, model, time),
  _time("time", "B lifetime", this, time),
  _iTag("iTag", "initial state tag", this, iTag),
  _fTag("fTag", "final state tag", this, fTag),
  _tagCat("tagCat", "tagging category", this),
  _tau("tau", "mean B lifetime", this, tau),
  _dGamma("dGamma", "Delta Gamma", this, dGamma),
  _dm("dm", "Delta mass", this, dm),
  _dilutions("dilutions", 0, this),
  _ADilWTags("ADilWTags", 0, this),
  _ANorm("ANorm", "normalization asymmetry", this, ANorm),
  _avgCEvenSum("avgCEvenSum", "weighted sum of average even coeffs", this),
  _avgCOddSum("avgCOddSum", "weighted sum of average odd coeffs", this),
  _avgCEvens("avgCEvens", 0, this),
  _avgCOdds("avgCOdds", 0, this),
  _tagCatCoefs("tagCatCoefs", 0, this),
  _coshCoef("coshCoef", "cosh coefficient", this, cosCoef),
  _sinhCoef("sinhCoef", "sinh coefficient", this),
  _cosCoef("cosCoef", "cos coefficient", this, cosCoef),
  _sinCoef("sinCoef", "sin coefficient", this),
  _createdVars("createdVars", 0, this),
  _decayType(type),
  _tagCatType(1),
  _tags(3),
  _iTagVal(2),
  _fTagVal(2),
  _checkVars(checkVars)
{
  // constructor with both initial and final state flavour tags
  // (decay into a flavour specific final state)

  initTag(kTRUE);
  initTag(kFALSE);

  RooArgList tagCatCoefs;
  RooArgList dilutions(dilution);
  RooArgList ADilWTags(ADilWTag);
  RooArgList avgCEvens(avgCEven);
  RooArgList avgCOdds(avgCOdd);
  initTaggingCats(tagCatCoefs, dilutions, ADilWTags, avgCEvens, avgCOdds);

  if (!checkVarDep(time, kTRUE)) assert(0);
  if (!checkVarDep(iTag, kTRUE)) assert(0);
  if (!checkVarDep(fTag, kTRUE)) assert(0);

  declareBases();
}

RooBTagDecay::RooBTagDecay(const char *name, const char* title,
    RooRealVar& time, RooCategory& iTag, RooAbsReal& tau,
    RooAbsReal& dGamma, RooAbsReal& dm, RooAbsReal& dilution,
    RooAbsReal& ADilWTag, RooAbsReal& avgCEven, RooAbsReal& avgCOdd,
    RooAbsReal& coshCoef, RooAbsReal& sinhCoef, RooAbsReal& cosCoef,
    RooAbsReal& sinCoef, const RooResolutionModel& model, DecayType type,
    Bool_t checkVars) :
  RooAbsAnaConvPdf(name, title, model, time),
  _time("time", "B lifetime", this, time),
  _iTag("iTag", "initial state tag", this, iTag),
  _fTag("fTag", "final state tag", this),
  _tagCat("tagCat", "tagging category", this),
  _tau("tau", "mean B lifetime", this, tau),
  _dGamma("dGamma", "Delta Gamma", this, dGamma),
  _dm("dm", "Delta mass", this, dm),
  _dilutions("dilutions", 0, this),
  _ADilWTags("ADilWTags", 0, this),
  _ANorm("ANorm", "normalization asymmetry", this),
  _avgCEvenSum("avgCEvenSum", "weighted sum of average even coeffs", this),
  _avgCOddSum("avgCOddSum", "weighted sum of average odd coeffs", this),
  _avgCEvens("avgCEvens", 0, this),
  _avgCOdds("avgCOdds", 0, this),
  _tagCatCoefs("tagCatCoefs", 0, this),
  _coshCoef("coshCoef", "cosh coefficient", this, coshCoef),
  _sinhCoef("sinhCoef", "sinh coefficient", this, sinhCoef),
  _cosCoef("cosCoef", "cos coefficient", this, cosCoef),
  _sinCoef("sinCoef", "sin coefficient", this, sinCoef),
  _createdVars("createdVars", 0, this),
  _decayType(type),
  _tagCatType(1),
  _tags(1),
  _iTagVal(2),
  _fTagVal(2),
  _checkVars(checkVars)
{
  // constructor with only an initial state flavour tag
  // (decay into CP self-conjugate state)

  initTag(kTRUE);

  RooArgList tagCatCoefs;
  RooArgList dilutions(dilution);
  RooArgList ADilWTags(ADilWTag);
  RooArgList avgCEvens(avgCEven);
  RooArgList avgCOdds(avgCOdd);
  initTaggingCats(tagCatCoefs, dilutions, ADilWTags, avgCEvens, avgCOdds);

  if (!checkVarDep(time, kTRUE)) assert(0);
  if (!checkVarDep(iTag, kTRUE)) assert(0);

  declareBases();
}

//_____________________________________________________________________________
RooBTagDecay::RooBTagDecay(const char *name, const char* title,
    RooRealVar& time, RooCategory& iTag, RooCategory& fTag,
    RooCategory& tagCat, RooAbsReal& tau, RooAbsReal& dGamma, RooAbsReal& dm,
    const RooArgList& dilutions, const RooArgList& ADilWTags,
    RooAbsReal& ANorm, const RooArgList& avgCEvens, const RooArgList& avgCOdds,
    const RooArgList& tagCatCoefs, RooAbsReal& cosCoef,
    const RooResolutionModel& model, DecayType type, Bool_t checkVars) :
  RooAbsAnaConvPdf(name, title, model, time),
  _time("time", "B lifetime", this, time),
  _iTag("iTag", "initial state tag", this, iTag),
  _fTag("fTag", "final state tag", this, fTag),
  _tagCat("tagCat", "tagging category", this, tagCat),
  _tau("tau", "mean B lifetime", this, tau),
  _dGamma("dGamma", "Delta Gamma", this, dGamma),
  _dm("dm", "Delta mass", this, dm),
  _dilutions("dilutions", 0, this),
  _ADilWTags("ADilWTags", 0, this),
  _ANorm("ANorm", "normalization asymmetry", this, ANorm),
  _avgCEvenSum("avgCEvenSum", "weighted sum of average even coeffs", this),
  _avgCOddSum("avgCOddSum", "weighted sum of average odd coeffs", this),
  _avgCEvens("avgCEvens", 0, this),
  _avgCOdds("avgCOdds", 0, this),
  _tagCatCoefs("tagCatCoefs", 0, this),
  _coshCoef("coshCoef", "cosh coefficient", this, cosCoef),
  _sinhCoef("sinhCoef", "sinh coefficient", this),
  _cosCoef("cosCoef", "cos coefficient", this, cosCoef),
  _sinCoef("sinCoef", "sin coefficient", this),
  _createdVars("createdVars", 0, this),
  _decayType(type),
  _tagCatType(2),
  _tags(3),
  _iTagVal(2),
  _fTagVal(2),
  _checkVars(checkVars)
{
  // constructor with both initial and final state flavour tags (decay into
  // a flavour specific final state) and with tagging categories

  initTag(kTRUE);
  initTag(kFALSE);

  initTaggingCats(tagCatCoefs, dilutions, ADilWTags, avgCEvens, avgCOdds);

  if (!checkVarDep(time, kTRUE))   assert(0);
  if (!checkVarDep(tagCat, kTRUE)) assert(0);
  if (!checkVarDep(iTag, kTRUE))   assert(0);
  if (!checkVarDep(fTag, kTRUE))   assert(0);

  declareBases();
}

//_____________________________________________________________________________
RooBTagDecay::RooBTagDecay(const char *name, const char* title,
    RooRealVar& time, RooCategory& iTag, RooCategory& tagCat, RooAbsReal& tau,
    RooAbsReal& dGamma, RooAbsReal& dm, const RooArgList& dilutions,
    const RooArgList& ADilWTags, const RooArgList& avgCEvens,
    const RooArgList& avgCOdds, const RooArgList& tagCatCoefs,
    RooAbsReal& coshCoef, RooAbsReal& sinhCoef, RooAbsReal& cosCoef,
    RooAbsReal& sinCoef, const RooResolutionModel& model, DecayType type,
    Bool_t checkVars) :
  RooAbsAnaConvPdf(name, title, model, time),
  _time("time", "B lifetime", this, time),
  _iTag("iTag", "initial state tag", this, iTag),
  _fTag("fTag", "final state tag", this),
  _tagCat("tagCat", "tagging category", this, tagCat),
  _tau("tau", "mean B lifetime", this, tau),
  _dGamma("dGamma", "Delta Gamma", this, dGamma),
  _dm("dm", "Delta mass", this, dm),
  _dilutions("dilutions", 0, this),
  _ADilWTags("ADilWTags", 0, this),
  _ANorm("ANorm", "normalization asymmetry", this),
  _avgCEvenSum("avgCEvenSum", "weighted sum of average even coeffs", this),
  _avgCOddSum("avgCOddSum", "weighted sum of average odd coeffs", this),
  _avgCEvens("avgCEvens", 0, this),
  _avgCOdds("avgCOdds", 0, this),
  _tagCatCoefs("tagCatCoefs", 0, this),
  _coshCoef("coshCoef", "cosh coefficient", this, coshCoef),
  _sinhCoef("sinhCoef", "sinh coefficient", this, sinhCoef),
  _cosCoef("cosCoef", "cos coefficient", this, cosCoef),
  _sinCoef("sinCoef", "sin coefficient", this, sinCoef),
  _createdVars("createdVars", 0, this),
  _decayType(type),
  _tagCatType(2),
  _tags(1),
  _iTagVal(2),
  _fTagVal(2),
  _checkVars(checkVars)
{
  // constructor with only an initial state flavour tag
  // (decay into CP self-conjugate state) and with tagging categories

  initTag(kTRUE);

  initTaggingCats(tagCatCoefs, dilutions, ADilWTags, avgCEvens, avgCOdds);

  if (!checkVarDep(time, kTRUE)) assert(0);
  if (!checkVarDep(tagCat, kTRUE)) assert(0);
  if (!checkVarDep(iTag, kTRUE)) assert(0);

  declareBases();
}

//_____________________________________________________________________________
RooBTagDecay::RooBTagDecay(const RooBTagDecay& other, const char* name) :
  RooAbsAnaConvPdf(other, name),
  _time("time", this, other._time),
  _iTag("iTag", this, other._iTag),
  _fTag("fTag", this, other._fTag),
  _tagCat("tagCat", this, other._tagCat),
  _tau("tau", this, other._tau),
  _dGamma("dGamma", this, other._dGamma),
  _dm("dm", this, other._dm),
  _dilutions("dilutions", this, other._dilutions),
  _ADilWTags("ADilWTags", this, other._ADilWTags),
  _ANorm("ANorm", this, other._ANorm),
  _avgCEvenSum("avgCEvenSum", this, other._avgCEvenSum),
  _avgCOddSum("avgCOddSum", this, other._avgCOddSum),
  _avgCEvens("avgCEvens", this, other._avgCEvens),
  _avgCOdds("avgCOdds", this, other._avgCOdds),
  _tagCatCoefs("tagCatCoefs", this, other._tagCatCoefs),
  _coshCoef("coshCoef", this, other._coshCoef),
  _sinhCoef("sinhCoef", this, other._sinhCoef),
  _cosCoef("cosCoef", this, other._cosCoef),
  _sinCoef("sinCoef", this, other._sinCoef),
  _createdVars("createdVars", 0, this),
  _tagCatPositions(other._tagCatPositions),
  _tagCatIndices(other._tagCatIndices),
  _coshBasis(other._coshBasis),
  _sinhBasis(other._sinhBasis),
  _cosBasis(other._cosBasis),
  _sinBasis(other._sinBasis),
  _decayType(other._decayType),
  _tagCatType(other._tagCatType),
  _tags(other._tags),
  _iTagVal(other._iTagVal),
  _fTagVal(other._fTagVal),
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

  // determine in which tagging category we are
  Int_t tagCatPos = -1;
  coefVal *= tagCatCoef(tagCatPos);

  // get value of initial state flavour tag
  Int_t iTagValue = _iTagVal;
  if (_iTagVal == 2 || _iTagVal == 3) iTagValue = _iTag;

  if (basisIndex == _coshBasis || basisIndex == _sinhBasis) {
    // terms that are even in the initial state tag
    if (iTagValue == 0) {
      // "untagged": dilution = 0 / integrate over initial state tag
      coefVal *= ((RooAbsReal*)_avgCEvens.at(tagCatPos))->getVal();
      if (_iTagVal == 0) coefVal *= 2.;
    } else {
      // calculate even coefficient with initial state tag
      coefVal *= ((RooAbsReal*)_avgCEvens.at(tagCatPos))->getVal()
          + iTagValue * ((RooAbsReal*)_dilutions.at(tagCatPos))->getVal()
          * (((RooAbsReal*)_avgCOdds.at(tagCatPos))->getVal()
          - ((RooAbsReal*)_ADilWTags.at(tagCatPos))->getVal()
          * ((RooAbsReal*)_avgCEvens.at(tagCatPos))->getVal());
    }

    if (_tags == 1) {
      // no final state tag
      return coefVal;
    } else {
      if (_fTagVal == 0) {
        // integrate over final state tag
        return 2. * coefVal;
      } else if (_fTagVal == 2) {
        // calculate even coefficient with varying final state tag
        return (1. - _fTag * _ANorm) * coefVal;
      } else {
        // calculate even coefficient with fixed final state tag
        return (1. - _fTagVal * _ANorm) * coefVal;
      }
    }

  } else if (basisIndex == _cosBasis || basisIndex == _sinBasis) {
    // terms that are odd in the initial state tag
    if (iTagValue == 0) {
      // "untagged": dilution = 0 / integrate over initial state tag
      coefVal *= ((RooAbsReal*)_avgCOdds.at(tagCatPos))->getVal();
      if (_iTagVal == 0) coefVal *= 2.;
    } else {
      // calculate odd coefficient with initial state tag
      coefVal *= ((RooAbsReal*)_avgCOdds.at(tagCatPos))->getVal()
          + iTagValue * ((RooAbsReal*)_dilutions.at(tagCatPos))->getVal()
          * (((RooAbsReal*)_avgCEvens.at(tagCatPos))->getVal()
          - ((RooAbsReal*)_ADilWTags.at(tagCatPos))->getVal()
          * ((RooAbsReal*)_avgCOdds.at(tagCatPos))->getVal());
    }

    if (_tags == 1) {
      // no final state tag
      return coefVal;
    } else {
      if (_fTagVal == 0) {
        // integrate over final state tag
        return -2. * _ANorm * coefVal;
      } else if (_fTagVal == 2) {
        // calculate even coefficient with varying final state tag
        return (_fTag - _ANorm) * coefVal;
      } else {
        // calculate even coefficient with fixed final state tag
        return (_fTagVal - _ANorm) * coefVal;
      }
    }
  }

  return 0.;
}

//_____________________________________________________________________________
Double_t RooBTagDecay::tagCatCoef() const
{
  Int_t category = -1;
  return tagCatCoef(category);
}

//_____________________________________________________________________________
Double_t RooBTagDecay::tagCatCoef(Int_t& category) const
{
  // check if we have tagging categories
  if (_tagCatType < 1) {
    coutE(Eval) << "RooBTagDecay::tagCatCoef(" << GetName()
        << ") this PDF does not contain tagging categories" << endl;
    if (category < 0) category = -1;
    return -1.;
  } else if (_tagCatType == 1) {
    if (category < 1) {
      category = 0;
      return 1.;
    } else {
      coutE(InputArguments) << "RooBTagDecay::tagCatCoef(" << GetName()
          << ") category " << category << "does not exist" << endl;
      return -1.;
    }
  }

  // get category
  if (category < 0) {
    category = getTagCatPosition(_tagCat.arg().getIndex());
  } else if (category >= _tagCat.arg().numTypes()) {
    coutE(InputArguments) << "RooBTagDecay::tagCatCoef(" << GetName()
        << ") category " << category << "does not exist" << endl;
    return -1.;
  }

  return tagCatCoefUnsafe(category);
}

//_____________________________________________________________________________
Double_t RooBTagDecay::tagCatCoefUnsafe(Int_t category) const
{
  // get tagging category coefficient
  Double_t tagCatCoef = ((RooAbsReal*)_tagCatCoefs.at(category))->getVal();

  // check coefficient
  if (tagCatCoef < 0. || (_tagCatType > 2 && tagCatCoef > 1.)) {
    coutE(Eval) << "RooBTagDecay::tagCatCoefUnsafe(" << GetName()
        << ") tagging category coefficient " << category
        << " is not in range (" << tagCatCoef << "): returning -1"
        << endl;
    return -1.;
  }

  return tagCatCoef;
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
    coefVars = _coshCoef.arg().getVariables();
  if (basisIndex == _sinhBasis && _tags < 2)
    coefVars = _sinhCoef.arg().getVariables();
  if (basisIndex == _cosBasis)
    coefVars = _cosCoef .arg().getVariables();
  if (basisIndex == _sinBasis && _tags < 2)
    coefVars = _sinCoef .arg().getVariables();

  if (coefVars == 0) {
    coefVars = new RooArgSet("parameters");
  } else {
    // add tagging variables
    std::auto_ptr<RooArgSet> tempSet;
    if (_tags > 0) {
      // add tagging category
      if (_tagCatType > 1) coefVars->add(_tagCat.arg());

      // add initial state tag
      coefVars->add(_iTag.arg());

      if (_tags > 1) {
        // add final state tag
        coefVars->add(_fTag.arg());
        tempSet.reset(_ANorm.arg().getVariables());
        coefVars->add((tempSet->getSize() > 0) ? *tempSet : _ANorm.arg());
      }

      // add average tagging category coefficient sums
      tempSet.reset(_avgCEvenSum.arg().getVariables());
      coefVars->add((tempSet->getSize() > 0) ? *tempSet: _avgCEvenSum.arg());
      tempSet.reset(_avgCOddSum.arg().getVariables());
      coefVars->add((tempSet->getSize() > 0) ? *tempSet : _avgCOddSum.arg());

      // add tagging category parameters
      Int_t  numTagCats = (_tagCatType > 1 ? _tagCat.arg().numTypes() : 1);
      for (Int_t tagCatIter = 0; tagCatIter < numTagCats; ++tagCatIter) {
        tempSet.reset( _dilutions.at(tagCatIter)->getVariables());
        coefVars->add((tempSet->getSize() > 0)
            ? *tempSet : *_dilutions.at(tagCatIter));

        tempSet.reset( _ADilWTags.at(tagCatIter)->getVariables());
        coefVars->add((tempSet->getSize() > 0)
            ? *tempSet :*_ADilWTags.at(tagCatIter));

        tempSet.reset(_avgCEvens.at(tagCatIter)->getVariables());
        coefVars->add(tempSet->getSize() > 0
            ? *tempSet : *_avgCEvens.at(tagCatIter));

        tempSet.reset(_avgCOdds.at(tagCatIter)->getVariables());
        coefVars->add( tempSet->getSize() > 0
            ? *tempSet : *_avgCOdds.at(tagCatIter));

        if (_tagCatType > 1) {
          tempSet.reset(_tagCatCoefs.at(tagCatIter)->getVariables());
          coefVars->add((tempSet->getSize() > 0)
              ? *tempSet : *_tagCatCoefs.at(tagCatIter));
        }
      }
    }

    // remove variables that are in basis functions
    TIterator* coefVarsIter = coefVars->createIterator();
    RooAbsArg* coefArg;
    while ((coefArg = (RooAbsArg*)coefVarsIter->Next()) != 0) {
      for (Int_t basisIter = 0; basisIter < _convSet.getSize(); ++basisIter) {
        if (_convSet.at(basisIter)->dependsOn(*coefArg))
          coefVars->remove(*coefArg);
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
      || (_tagCatType > 1 && !checkVarDep(_tagCat.arg()))
      || (_tags > 0 && (!checkVarDep(_iTag.arg()) || !checkTag(kTRUE)))
      || (_tags > 1 && (!checkVarDep(_fTag.arg()) || !checkTag(kFALSE))))
    return 0;

  // initialize integration code variables
  Int_t  intCode   = 0;
  Bool_t intTagCat = kFALSE;
  Bool_t intITag   = kFALSE;
  Bool_t intFTag   = kFALSE;
  RooArgSet intVars = allVars;

  // integrate over tagging category?
  if (_tagCatType > 1) intTagCat = intVars.remove(_tagCat.arg(), kTRUE, kTRUE);

  // integrate over initial state tag?
  if (_tags > 0) intITag = intVars.remove(_iTag.arg(), kTRUE, kTRUE);

  // integrate over final state tag?
  if (_tags > 1) intFTag = intVars.remove(_fTag.arg(), kTRUE, kTRUE);

  // get integration code
  if (intVars.getSize() > 0) {
    if (coef == _coshBasis) {
      intCode = 8 * _coshCoef.arg().getAnalyticalIntegral(intVars, analVars,
          rangeName);
    } else if (coef == _sinhBasis) {
      if (_tags > 1)
        return 0;
      else
        intCode = 8 * _sinhCoef.arg().getAnalyticalIntegral(intVars, analVars,
            rangeName);
    } else if (coef == _cosBasis) {
      intCode = 8 * _cosCoef.arg().getAnalyticalIntegral(intVars, analVars,
          rangeName);
    } else if (coef == _sinBasis) {
      if (_tags > 1)
        return 0;
      else
        intCode = 8 * _sinCoef.arg().getAnalyticalIntegral(intVars, analVars,
            rangeName);
    }
  }

  // return the integration code if there are no explicit tags
  if (_tags == 0) return intCode;

  // check that tagging parameters don't depend on integration variables
  TIterator* analVarIter = analVars.createIterator();
  RooAbsArg* analVar = 0;
  while ((analVar = (RooAbsArg*)analVarIter->Next()) != 0) {
    if (!checkVarDep(*analVar, kFALSE, kTRUE)) return 0;
  }
  delete analVarIter;

  if (intTagCat && (rangeName == 0 || !_tagCat.hasRange(rangeName))) {
    // add the tagging category to integration variables
    intCode += 4;
    analVars.add(_tagCat.arg());

    // initialize tagging category maps
    if (_tagCatIndices.empty()) initTagCatMaps();
  }

  if (intITag && (rangeName == 0 || !_iTag.hasRange(rangeName))) {
    // add the initial state tag to integration variables
    intCode += 2;
    analVars.add(_iTag.arg());
  }

  if (intFTag && (rangeName == 0 || !_fTag.hasRange(rangeName))) {
    // add the final state tag to integration variables
    intCode += 1;
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
  Int_t  coefCode = TMath::FloorNint(code / 8.);
  Bool_t intTagCat = (Bool_t)(code & 4);
  Bool_t intITag = (Bool_t)(code & 2);
  Bool_t intFTag = (Bool_t)(code & 1);

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

  // return the integral if we don't have to evaluate explicit tags
  if (_tags == 0 || coefInt == 0.) return coefInt;

  // get value of initial state flavour tag
  Int_t iTagValue = _iTagVal;
  if (_iTagVal == 2 || _iTagVal == 3) iTagValue = _iTag;

  // determine in which tagging category we are
  Int_t tagCatPos = -1;
  if (_tagCatType > 1 && !intTagCat) coefInt *= tagCatCoef(tagCatPos);
  else if (_tagCatType == 1) tagCatPos = 0;

  if (coef == _coshBasis || coef == _sinhBasis) {
    // terms that are even in the initial state tag
    if ((intITag && (_iTagVal == 2 || _iTagVal == 3)) || iTagValue == 0) {
      // integrate over the initial state tag
      // (or dilution = 0 for a tag with three states)
      if (intTagCat) coefInt *= _avgCEvenSum;
      else coefInt *= ((RooAbsReal*)_avgCEvens.at(tagCatPos))->getVal();

      if (intITag || _iTagVal == 0) coefInt *= 2.;
    } else if (intTagCat) {
      // don't integrate over initial state tag, but over tagging category
      Double_t catInt = _avgCEvenSum;
      std::map<Int_t, Int_t>::const_iterator catIt;
      for (catIt = _tagCatIndices.begin(); catIt != _tagCatIndices.end();
          ++catIt) {
        catInt += tagCatCoefUnsafe(catIt->first)
            * (iTagValue * ((RooAbsReal*)_dilutions.at(catIt->first))->getVal()
            * (((RooAbsReal*)_avgCOdds.at(catIt->first))->getVal()
            - ((RooAbsReal*)_ADilWTags.at(catIt->first))->getVal()
            * ((RooAbsReal*)_avgCEvens.at(catIt->first))->getVal()));
      }

      coefInt *= catInt;
    } else {
      // don't integrate over the initial state tag or the tagging category
      coefInt *= ((RooAbsReal*)_avgCEvens.at(tagCatPos))->getVal()
          + iTagValue * ((RooAbsReal*)_dilutions.at(tagCatPos))->getVal()
          * (((RooAbsReal*)_avgCOdds.at(tagCatPos))->getVal()
          - ((RooAbsReal*)_ADilWTags.at(tagCatPos))->getVal()
          * ((RooAbsReal*)_avgCEvens.at(tagCatPos))->getVal());
    }

    // integrate over final state tag if required
    if (_tags == 1)
      return coefInt;
    else if ((intFTag && _fTagVal == 2) || _fTagVal == 0)
      return 2. * coefInt;
    else if (_fTagVal == 2)
      return (1. - _fTag * _ANorm) * coefInt;
    else
      return (1. - _fTagVal * _ANorm) * coefInt;

  } else if (coef == _cosBasis || coef == _sinBasis) {
    // terms that are odd in the initial state tag
    if ((intITag && (_iTagVal == 2 || _iTagVal == 3)) || iTagValue == 0) {
      // integrate over the initial state tag
      // (or dilution = 0 for a tag with three states)
      if (intTagCat) coefInt *= _avgCOddSum;
      else coefInt *= ((RooAbsReal*)_avgCOdds.at(tagCatPos))->getVal();

      if (intITag || _iTagVal == 0) coefInt *= 2.;
    } else if (intTagCat) {
      // don't integrate over initial state tag, but over tagging category
      Double_t catInt = _avgCOddSum;
      std::map<Int_t, Int_t>::const_iterator catIt;
      for (catIt = _tagCatIndices.begin(); catIt != _tagCatIndices.end();
          ++catIt) {
        catInt += tagCatCoefUnsafe(catIt->first)
            * (iTagValue * ((RooAbsReal*)_dilutions.at(catIt->first))->getVal()
            * (((RooAbsReal*)_avgCEvens.at(catIt->first))->getVal()
            - ((RooAbsReal*)_ADilWTags.at(catIt->first))->getVal()
            * ((RooAbsReal*)_avgCOdds.at(catIt->first))->getVal()));
      }

      coefInt *= catInt;
    } else {
      // don't integrate over the initial state tag or the tagging category
      coefInt *= ((RooAbsReal*)_avgCOdds.at(tagCatPos))->getVal()
          + iTagValue * ((RooAbsReal*)_dilutions.at(tagCatPos))->getVal()
          * (((RooAbsReal*)_avgCEvens.at(tagCatPos))->getVal()
          - ((RooAbsReal*)_ADilWTags.at(tagCatPos))->getVal()
          * ((RooAbsReal*)_avgCOdds.at(tagCatPos))->getVal());
    }

    // integrate over final state tag if required
    if (_tags == 1)
      return coefInt;
    else if ((intFTag && _fTagVal == 2) || _fTagVal == 0)
      return -2. * _ANorm * coefInt;
    else if (_fTagVal == 2)
      return (_fTag - _ANorm) * coefInt;
    else
      return (_fTagVal - _ANorm) * coefInt;
  }

  return 0.;
}

//_____________________________________________________________________________
Int_t RooBTagDecay::getGenerator(const RooArgSet& directVars,
    RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  // copy variables that can be generated directly by RooBTagDecay from
  // directVars to generateVars and return generation code

  // use accept/reject in certain cases
  if (!_checkVars || !checkVarDep(_time.arg())
      || (_tagCatType > 1 && !checkVarDep(_tagCat.arg()))
      || (_tags > 0 && ((_iTagVal != 2 && _iTagVal != 3)
          || !checkVarDep(_iTag.arg()) || !checkTag(kTRUE)))
      || (_tags > 1 && (_fTagVal != 2 || !checkVarDep(_fTag.arg())
          || !checkTag(kFALSE))))
    return 0;

  Int_t genCode = 0;

  // find time variable
  RooAbsArg* arg = directVars.find(_time.arg().GetName());
  if (arg == 0) return genCode;

  genCode += 1;
  generateVars.add(*arg);

  if (_tags == 0) return genCode;

  // check which tagging variables to generate
  arg = 0;
  if (_tagCatType > 1) {
    // find the tagging category variable
    arg = directVars.find(_tagCat.arg().GetName());
    if (arg == 0) return genCode;
  }

  // find initial state tag variable
  RooAbsArg* arg1 = directVars.find(_iTag.arg().GetName());
  if (arg1 == 0) return genCode;

  RooAbsArg* arg2 = 0;
  if (_tags > 1) {
    // find final state tag variable
    arg2 = directVars.find(_fTag.arg().GetName());
    if (arg2 == 0) return genCode;
  }

  // build generation code and add variables to generation variables set
  if (arg != 0) {
    // generate tagging category
    genCode += 2;
    generateVars.add(*arg);
  }

  // generate initial state tag
  genCode += 4;
  generateVars.add(*arg1);

  if (arg2 != 0) {
    // generate final state tag
    genCode += 8;
    generateVars.add(*arg2);
  }

  return genCode;
}

//_____________________________________________________________________________
void RooBTagDecay::generateEvent(Int_t code)
{
  // generate values for the variables corresponding to the generation code

  // set minimum Gamma for time envelope: Gamma - |DeltaGamma| / 2
  Double_t gammaMin  = 1. / _tau - TMath::Abs(_dGamma) / 2.;

  // generate event
  while(true) {
    // generate time variable with the exponential envelope function
    Double_t timeGen = -log(RooRandom::uniform()) / gammaMin;
    if (_decayType == Flipped || (_decayType == DoubleSided
        && RooRandom::uniform() > 0.5))
      timeGen *= -1.;
    if (timeGen < _time.min() || timeGen > _time.max()) continue;

    // get (unnormalized) PDF and envelope values for generated time
    Double_t tAbs         = TMath::Abs(timeGen);
    Double_t pdfVal       = 0.;
    Double_t envVal       = 0.;
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
        pdfVal = oddTerms + evenTerms;
        envVal = evenTermsEnv + oddTermsEnv;

      } else {
        // use PDF with flavour tag(s) and tagging categories
        Int_t tagCatPos = -1;
        Double_t tagCoef = tagCatCoef(tagCatPos);

        // check tagging category coefficient
        if (tagCoef < 0.) {
          coutF(Generation) << "RooBTagDecay::generateEvent(" << GetName()
              << ") not a valid tagging category coefficient"
              << endl;
          assert(0);
        }

        Double_t cEven = tagCoef
          * (((RooAbsReal*)_avgCEvens.at(tagCatPos))->getVal()
          + _iTag * ((RooAbsReal*)_dilutions.at(tagCatPos))->getVal()
          * (((RooAbsReal*)_avgCOdds.at(tagCatPos))->getVal()
          - ((RooAbsReal*)_ADilWTags.at(tagCatPos))->getVal()
          * ((RooAbsReal*)_avgCEvens.at(tagCatPos))->getVal()));

        Double_t cOdd = tagCoef
          * (((RooAbsReal*)_avgCOdds.at(tagCatPos))->getVal()
          + _iTag * ((RooAbsReal*)_dilutions.at(tagCatPos))->getVal()
          * (((RooAbsReal*)_avgCEvens.at(tagCatPos))->getVal()
          - ((RooAbsReal*)_ADilWTags.at(tagCatPos))->getVal()
          * ((RooAbsReal*)_avgCOdds.at(tagCatPos))->getVal()));

        if (_tags > 1) {
          cEven *= 1. - _fTag * _ANorm;
          cOdd  *= _fTag - _ANorm;
        }

        pdfVal = cEven * evenTerms + cOdd * oddTerms;
        envVal = TMath::Abs(cEven) * evenTermsEnv
            + TMath::Abs(cOdd) * oddTermsEnv;
      }

    } else {
      // generate time, tagging category and flavour tag(s):
      // use PDF integrated over tagging variables
      if (_tags == 1) {
        pdfVal = _avgCEvenSum * evenTerms + _avgCOddSum * oddTerms;
        envVal = TMath::Abs(_avgCEvenSum) * evenTermsEnv
            + TMath::Abs(_avgCOddSum) * oddTermsEnv;
      } else {
        pdfVal = _avgCEvenSum * evenTerms - _ANorm * _avgCOddSum * oddTerms;
        envVal = TMath::Abs(_avgCEvenSum) * evenTermsEnv
            + TMath::Abs(_ANorm * _avgCOddSum) * oddTermsEnv;
      }
    }

    // check if PDF value is positive
    if (pdfVal < 0.) {
      coutF(Generation) << "RooBTagDecay::generateEvent(" << GetName()
          << ") PDF is negative for generated time value" << endl;
      assert(0);
    }

    // accept/reject generated time value, using the envelope as maximum
    if (exp(-tAbs / _tau) * pdfVal
        < RooRandom::uniform() * exp(-gammaMin * tAbs) * envVal)
      continue;

    // set time value
    _time = timeGen;

    // exit generation loop if we don't generate tagging variables
    if (code == 1) break;

    Int_t tagCatGen = 0;
    Double_t catPdfSum = 0.;
    if (_tagCatType > 1) {
      // generate value for the tagging category
      Double_t rand = RooRandom::uniform();
      while (tagCatGen < _tagCat.arg().numTypes()) {
        // get tagging category coefficient
        Double_t tagCoef = tagCatCoefUnsafe(tagCatGen);

        // check coefficient
        if (tagCoef < 0.) {
          coutF(Generation) << "RooBTagDecay::generateEvent(" << GetName()
              << ") not a valid tagging category coefficient"
              << endl;
          assert(0);
        }

        // add the value of the (tags-integrated) PDF to the categories sum
        if (_tags == 1)
          catPdfSum += tagCoef *
              (((RooAbsReal*)_avgCEvens.at(tagCatGen))->getVal() * evenTerms
              + ((RooAbsReal*)_avgCOdds.at(tagCatGen))->getVal() * oddTerms);
        else
          catPdfSum += tagCoef *
              (((RooAbsReal*)_avgCEvens.at(tagCatGen))->getVal() * evenTerms
              - _ANorm * ((RooAbsReal*)_avgCOdds.at(tagCatGen))->getVal()
              * oddTerms);

        if (catPdfSum < rand * pdfVal) ++tagCatGen;
        else break;
      }

      // check generated value
      if (tagCatGen >= _tagCat.arg().numTypes()) {
        coutW(Generation) << "RooBTagDecay::generateEvent(" << GetName()
            << ") generation of event aborted due to a numerical problem in generation of the tagging category"
            << endl;
        continue;
      }

      // set tagging category value
      _tagCat = getTagCatIndex(tagCatGen);
    }

    // calculate dilution * (coef_(O/E) - ADilWTag * coef_(E/O))
    Double_t cEvenDil = ((RooAbsReal*)_dilutions.at(tagCatGen))->getVal()
        * (((RooAbsReal*)_avgCOdds.at(tagCatGen))->getVal()
        - ((RooAbsReal*)_ADilWTags.at(tagCatGen))->getVal()
        * ((RooAbsReal*)_avgCEvens.at(tagCatGen))->getVal());
    Double_t cOddDil = ((RooAbsReal*)_dilutions.at(tagCatGen))->getVal()
        * (((RooAbsReal*)_avgCEvens.at(tagCatGen))->getVal()
        - ((RooAbsReal*)_ADilWTags.at(tagCatGen))->getVal()
        * ((RooAbsReal*)_avgCOdds.at(tagCatGen))->getVal());

    // generate value for initial state tag
    Int_t iTagGen = 1;
    if (_iTagVal == 3 && tagCatGen == 0) {
      // tagging category 0 is interpreted as the "untagged" category
      // if the initial state tag can assume values +1, -1 and 0!!!!!!!!!
      iTagGen = 0;
    } else {
      // use CP asymmetry if initial state tag only contains values +1 and -1
      Double_t ACP = 0.;
      if (_tags == 1)
        ACP = (cEvenDil * evenTerms + cOddDil * oddTerms)
        / (((RooAbsReal*)_avgCEvens.at(tagCatGen))->getVal() * evenTerms
        + ((RooAbsReal*)_avgCOdds.at(tagCatGen))->getVal() * oddTerms);
      else
        ACP = (cEvenDil * evenTerms - _ANorm * cOddDil * oddTerms)
        / (((RooAbsReal*)_avgCEvens.at(tagCatGen))->getVal() * evenTerms
        - _ANorm * ((RooAbsReal*)_avgCOdds.at(tagCatGen))->getVal()*oddTerms);

      if (2. * RooRandom::uniform() > 1. + ACP) iTagGen = -1;
    }

    // set initial state tag value
    _iTag = iTagGen;

    if (_tags > 1) {
      // generate value for final state tag
      Int_t fTagGen = 1;
      Double_t cEven = ((RooAbsReal*)_avgCEvens.at(tagCatGen))->getVal()
          + iTagGen * cEvenDil;
      Double_t cOdd  = ((RooAbsReal*)_avgCOdds.at(tagCatGen))->getVal()
          + iTagGen * cOddDil;
      Double_t Af = (-_ANorm * cEven * evenTerms + cOdd * oddTerms)
          / (cEven * evenTerms - _ANorm * cOdd * oddTerms);

      if (2. * RooRandom::uniform() > 1. + Af) fTagGen = -1;

      // set final state tag value
      _fTag = fTagGen;
    }

    // exit generation loop
    break;
  }
}

//_____________________________________________________________________________
Int_t RooBTagDecay::getTagCatPosition(Int_t tagCatIndex) const
{
  if (_tagCatPositions.empty()) initTagCatMaps();

  return _tagCatPositions.find(tagCatIndex)->second;
}

//_____________________________________________________________________________
Int_t RooBTagDecay::getTagCatIndex(Int_t tagCatPosition) const
{
  if (_tagCatIndices.empty()) initTagCatMaps();

  return _tagCatIndices.find(tagCatPosition)->second;
}

//_____________________________________________________________________________
void RooBTagDecay::initTaggingCats(const RooArgList& tagCatCoefs,
    const RooArgList& dilutions, const RooArgList& ADilWTags,
    const RooArgList& avgCEvens, const RooArgList& avgCOdds)
{
  // set proxies for tagging category coefficients, dilutions,
  // dilution/wrong-tag asymmetries, average even coefficients and
  // average odd coefficients

  // set tagging categories type
  Int_t numTagCats = 0;
  if (_tagCatType == 1) numTagCats = 1;
  else if (_tagCatType > 1) numTagCats = _tagCat.arg().numTypes();

  if (numTagCats == 0) {
    coutW(InputArguments) << "RooBTagDecay::initTaggingCats(" << GetName()
        << ") there are no tagging categories" << endl;
    _tagCatType = 0;
    return;
  } else if (numTagCats == 1) {
    _tagCatType = 1;
    if (tagCatCoefs.getSize() != 0) 
      coutW(InputArguments) << "RooBTagDecay::initTaggingCats(" << GetName()
          << ") there is only one tagging category: category coefficient(s) will not be used"
          << endl;
  } else if (tagCatCoefs.getSize() == numTagCats) {
    _tagCatType = 2;
  } else if (tagCatCoefs.getSize() == numTagCats - 1) {
    _tagCatType = 3;
  } else {
    coutF(InputArguments) << "RooBTagDecay::initTaggingCats(" << GetName()
        << ") number of tagging category coefficients does not match number of category types"
        << endl;
    assert(0);
  }

  // check number of dilutions
  if (dilutions.getSize() != numTagCats) {
    coutF(InputArguments) << "RooBTagDecay::initTaggingCats(" << GetName()
        << ") number of dilutions does not match number of category types"
        << endl;
    assert(0);
  }

  // check number of dilution/wrong-tag asymmetries
  if (ADilWTags.getSize() != numTagCats) {
    coutF(InputArguments) << "RooBTagDecay::initTaggingCats(" << GetName()
        << ") number of dilution/wrong-tag asymmetries does not match number of category types"
        << endl;
    assert(0);
  }

  // check number of average even coefficients
  if (avgCEvens.getSize() != numTagCats) {
    coutF(InputArguments) << "RooBTagDecay::initTaggingCats(" << GetName()
        << ") number of average even coefficients does not match number of category types"
        << endl;
    assert(0);
  }

  // check number of average odd coefficients
  if (avgCOdds.getSize() != numTagCats) {
    coutF(InputArguments) << "RooBTagDecay::initTaggingCats(" << GetName()
        << ") number of average odd coefficients does not match number of category types"
        << endl;
    assert(0);
  }

  // initialize tagging category maps
  initTagCatMaps();

  // loop over tagging categories
  RooArgList avgCEvenList;
  RooArgList avgCOddList;
  RooArgList tagCatCoefList;
  TString avgCoefFormString;
  TString tagCatCoefFormString("1.");
  for (Int_t tagCatIter = 0; tagCatIter < numTagCats; ++tagCatIter) {
    // get parameters
    RooAbsReal* tagCatCoef = 0;
    if (_tagCatType == 2)
      tagCatCoef = dynamic_cast<RooAbsReal*>(tagCatCoefs.at(tagCatIter));
    else if (_tagCatType > 2 && tagCatIter > 0)
      tagCatCoef = dynamic_cast<RooAbsReal*>(tagCatCoefs.at(tagCatIter - 1));
    RooAbsReal* dilution = dynamic_cast<RooAbsReal*>(dilutions.at(tagCatIter));
    RooAbsReal* ADilWTag = dynamic_cast<RooAbsReal*>(ADilWTags.at(tagCatIter));
    RooAbsReal* avgCEven = dynamic_cast<RooAbsReal*>(avgCEvens.at(tagCatIter));
    RooAbsReal* avgCOdd  = dynamic_cast<RooAbsReal*>(avgCOdds.at(tagCatIter));

    // get tagging category coefficient
    if (tagCatCoef != 0) {
      tagCatCoefList.add(*tagCatCoef);
    } else if (_tagCatType == 2 || tagCatIter > 0) {
      coutF(InputArguments) << "RooBTagDecay::initTaggingCats("
          << GetName() << ") tagging category coefficient '"
          << tagCatCoefs.at(tagCatIter)->GetName()
          << "' is not a real-valued variable" << endl;
      assert(0);
    }

    // get dilution
    if (dilution != 0) _dilutions.add(*dilution);
    else {
      coutF(InputArguments) << "RooBTagDecay::initTaggingCats("
          << GetName() << ") dilution '"
          << dilutions.at(tagCatIter)->GetName()
          << "' is not a real-valued variable" << endl;
      assert(0);
    }

    // get dilution/wrong-tag asymmetry
    if (ADilWTag != 0) _ADilWTags.add(*ADilWTag);
    else {
      coutF(InputArguments) << "RooBTagDecay::initTaggingCats("
          << GetName() << ") dilution/wrong-tag asymmetry '"
          << ADilWTags.at(tagCatIter)->GetName()
          << "' is not a real-valued variable" << endl;
      assert(0);
    }

    // get average even coefficient
    if (avgCEven != 0) {
      if (_tagCatType == 2 || (_tagCatType == 3 && tagCatIter > 0)) {
        avgCEvenList.add(*avgCEven);

        Int_t formStringIndex = 0;
        if (avgCoefFormString == "") {
          avgCoefFormString += "@";
        } else {
          avgCoefFormString += " + @";
          if (_tagCatType == 2) formStringIndex = tagCatIter;
          else formStringIndex = tagCatIter - 1;
        }

        avgCoefFormString += formStringIndex + numTagCats;
        avgCoefFormString += "* @";
        avgCoefFormString += formStringIndex;

        tagCatCoefFormString += " - @";
        tagCatCoefFormString += formStringIndex;
      } else
        _avgCEvenSum.setArg(*avgCEven);
    } else {
      coutF(InputArguments) << "RooBTagDecay::initTaggingCats("
          << GetName() << ") average even coefficient '"
          << dilutions.at(tagCatIter)->GetName()
          << "' is not a real-valued variable" << endl;
      assert(0);
    }

    // get average odd coefficient
    if (avgCOdd != 0) {
      if (_tagCatType == 2 || (_tagCatType == 3 && tagCatIter > 0))
        avgCOddList.add(*avgCOdd);
      else
        _avgCOddSum.setArg(*avgCOdd);
    } else {
      coutF(InputArguments) << "RooBTagDecay::initTaggingCats("
          << GetName() << ") average odd coefficient '"
          << dilutions.at(tagCatIter)->GetName()
          << "' is not a real-valued variable" << endl;
      assert(0);
    }
  }

  // set average even and average odd coefficients
  if (_tagCatType == 1) {
    _avgCEvens.add(*_avgCEvenSum.absArg());
    _avgCOdds.add(*_avgCOddSum.absArg());
  } else if (_tagCatType == 2) {
    // set average even, average odd and category coefficients lists
    _avgCEvens.add(avgCEvenList);
    _avgCOdds.add(avgCOddList);
    _tagCatCoefs.add(tagCatCoefList);

    // set the sums of the average even and average odd coefficients
    avgCEvenList.add(tagCatCoefList);
    avgCOddList.add(tagCatCoefList);
    RooFormulaVar* avgCEvenSum = new RooFormulaVar("avgCEvenSum",
        "avgCEvenSum", avgCoefFormString, avgCEvenList);
    RooFormulaVar* avgCOddSum = new RooFormulaVar("avgCOddSum",
        "avgCOddSum", avgCoefFormString, avgCOddList);
    _createdVars.addOwned(*avgCEvenSum);
    _createdVars.addOwned(*avgCOddSum);
    _avgCEvenSum.setArg(*avgCEvenSum);
    _avgCOddSum.setArg(*avgCOddSum);

  } else if (_tagCatType == 3) {
    // create tagging category coefficient for complementary category
    RooFormulaVar* tagCatCoef0 = new RooFormulaVar("tagCatCoef0",
        "tagCatCoef0", tagCatCoefFormString, tagCatCoefList);
    _createdVars.addOwned(*tagCatCoef0);
    _tagCatCoefs.add(*tagCatCoef0);
    _tagCatCoefs.add(tagCatCoefList);

    // create even and odd coefficients for complementary category
    TString formString("1. / @");
    formString += avgCEvenList.getSize() + tagCatCoefList.getSize() + 1;
    formString += " * (@";
    formString += avgCEvenList.getSize();
    avgCoefFormString = formString + " - (" + avgCoefFormString + "))";

    RooArgList tempEvenList(avgCEvenList);
    tempEvenList.add(*_avgCEvenSum.absArg());
    tempEvenList.add(tagCatCoefList);
    tempEvenList.add(*tagCatCoef0);
    RooFormulaVar* avgCEven0 = new RooFormulaVar("avgCEven0",
        "avgCEven0", avgCoefFormString, tempEvenList);
    _createdVars.addOwned(*avgCEven0);
    _avgCEvens.add(*avgCEven0);
    _avgCEvens.add(avgCEvenList);

    RooArgList tempOddList(avgCOddList);
    tempOddList.add(*_avgCOddSum.absArg());
    tempOddList.add(tagCatCoefList);
    tempOddList.add(*tagCatCoef0);
    RooFormulaVar* avgCOdd0 = new RooFormulaVar("avgCOdd0",
        "avgCOdd0", avgCoefFormString, tempOddList);
    _createdVars.addOwned(*avgCOdd0);
    _avgCOdds.add(*avgCOdd0);
    _avgCOdds.add(avgCOddList);
  }
}

//_____________________________________________________________________________
void RooBTagDecay::initTag(Bool_t iTag)
{
    // initialize flavour tag

    // get tag parameters
    const char* tagName;
    Int_t  tagVal = 2, numTypes = 0;
    Bool_t BIndex= kFALSE, BbarIndex = kFALSE, noTagIndex = kFALSE;
    if (iTag) {
      tagName    = _iTag.GetName();
      numTypes   = _iTag.arg().numTypes();
      BIndex     = _iTag.arg().isValidIndex(+1);
      BbarIndex  = _iTag.arg().isValidIndex(-1);
      noTagIndex = _iTag.arg().isValidIndex(0);
    } else {
      tagName   = _fTag.GetName();
      numTypes  = _fTag.arg().numTypes();
      BIndex    = _fTag.arg().isValidIndex(+1);
      BbarIndex = _fTag.arg().isValidIndex(-1);
    }

    if (!_checkVars) {
      // print warning message if flavour tags are not to be checked
      coutW(InputArguments) << "RooBTagDecay::initTag(" << GetName()
          << ") unchecked flavour tag "
          << tagName
          << " (use of \"checkVars = false\" is not recommended): integrals will be calculated numerically and event generation will use accept/reject"
          << endl;
    } else if (numTypes == 2 && BIndex && BbarIndex) {
      // both B and Bbar
      tagVal = 2;
    } else if (numTypes == 1 && BIndex) {
      // only B
      tagVal = +1;
    } else if (numTypes == 1 && BbarIndex) {
      // only Bbar
      tagVal = -1;
    } else if (iTag && numTypes == 1 && noTagIndex) {
      // sum of B and Bbar
      tagVal = 0;
    } else if (iTag && numTypes == 3  && BIndex && BbarIndex && noTagIndex) {
      // B, Bbar and sum of B and Bbar
      tagVal = 3;
      coutW(InputArguments) << "RooBTagDecay::initTag(" << GetName()
          << ") initial state tag can assume three values (-1, 0, +1): please make sure you know what you are doing:"
          << endl
          << "    * the value of the decay rate with tag = 0 is equal to the value with dilution = 0"
          << endl
          << "    * the integral of the decay rate is equal to the sum of B (+1) and Bbar (-1)"
          << endl
          << "    * events with tag = 0 will be generated in tagging category 0"
          << endl;
    } else {
      // not a valid configuration
      coutE(InputArguments) << "RooBTagDecay::initTag(" << GetName()
          << ") flavour tag \"" << tagName
          << "\" can only have values +1, -1 and 0: use \"checkVars = false\" if you insist on using other/additional values"
          << endl;
      assert(0);
    }

    // set tag value
    if (iTag) _iTagVal = tagVal;
    else _fTagVal = tagVal;
}

//_____________________________________________________________________________
void RooBTagDecay::initTagCatMaps() const
{
  // check number of categories
  if (_tagCatType < 2) return;

  // loop over tagging categories
  TIterator*  tagCatTypeIter = _tagCat.arg().typeIterator();
  RooCatType* tagCatIndex    = 0;
  Int_t       position       = 0;
  while ((tagCatIndex = (RooCatType*)tagCatTypeIter->Next()) != 0) {
    // set tagging category position and index
    Int_t index = tagCatIndex->getVal();
    _tagCatPositions[index] = position;
    _tagCatIndices[position] = index;

    // update position
    ++position;
  }

  delete tagCatTypeIter;
}

//_____________________________________________________________________________
void RooBTagDecay::declareBases()
{
  // create basis functions for time variable

  if (_decayType == SingleSided) {
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

  } else if (_decayType == Flipped) {
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

  } else if (_decayType == DoubleSided) {
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

//_____________________________________________________________________________
Bool_t RooBTagDecay::checkVarDep(const RooAbsArg& var, Bool_t warn,
    Bool_t onlyTagPars) const
{
  // check if parameters depend on specified variable

  if (_checkVars) {
    Bool_t checks = kTRUE;
    if ((!onlyTagPars && (_tau.arg().dependsOn(var)
        || _dGamma.arg().dependsOn(var) || _dm.arg().dependsOn(var)
        || _coshCoef.arg().dependsOn(var) || _cosCoef.arg().dependsOn(var)
        || (_tags < 2 && (_sinhCoef.arg().dependsOn(var)
        || _sinCoef.arg().dependsOn(var)))))
        || (_tags > 1 && _ANorm.arg().dependsOn(var))
        || (_tags > 0 && (_avgCEvenSum.arg().dependsOn(var)
        || _avgCOddSum.arg().dependsOn(var))))
      checks = kFALSE;

    if (checks && _tags > 0) {
      Int_t numTagCats = 1;
      if (_tagCatType > 1) numTagCats = _tagCat.arg().numTypes();
      for (Int_t tagCatIter = 0; tagCatIter < numTagCats; ++tagCatIter) {
        if (_dilutions.at(tagCatIter)->dependsOn(var)
            || _ADilWTags.at(tagCatIter)->dependsOn(var)
            || _avgCEvens.at(tagCatIter)->dependsOn(var)
            || _avgCOdds.at(tagCatIter)->dependsOn(var)
            || _tagCatCoefs.at(tagCatIter)->dependsOn(var))
          checks = kFALSE;
      }
    }

    if (!checks) 
      coutW(InputArguments) << "RooBTagDecay::checkVarDep(" << GetName()
          << ") parameters depend on " << var.GetName()
          << ": use \"checkVars = false\" if you insist on using this dependence"
          << endl;

    return checks;

  } else if (warn) {
    coutW(InputArguments) << "RooBTagDecay::checkVarDep(" << GetName()
        << ") parameters dependence on "
        << var.GetName()
        << " is unchecked (use of \"checkVars = false\" is not recommended): integrals will be calculated numerically and event generation will use accept/reject"
        << endl;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t RooBTagDecay::checkTag(Bool_t iTag) const
{
  // check flavour tag state
  if (_checkVars) {
    // get tag parameters
    Int_t  tagVal = 2, nTypes = 0;
    Bool_t hasB = kFALSE, hasBbar = kFALSE, hasNoTag = kFALSE;
    if (iTag) {
      tagVal   = _iTagVal;
      nTypes   = _iTag.arg().numTypes();
      hasB     = _iTag.arg().isValidIndex(+1);
      hasBbar  = _iTag.arg().isValidIndex(-1);
      hasNoTag = _iTag.arg().isValidIndex(0);
    } else {
      tagVal  = _fTagVal;
      nTypes  = _fTag.arg().numTypes();
      hasB    = _fTag.arg().isValidIndex(+1);
      hasBbar = _fTag.arg().isValidIndex(-1);
    }

    // check tag
    if ((tagVal == 2 && (nTypes != 2 || !hasB || !hasBbar))
        || (tagVal == +1 && (nTypes != 1 || !hasB))
        || (tagVal == -1 && (nTypes != 1 || !hasBbar))
        || (tagVal ==  3 && (nTypes != 3 || !hasB || !hasBbar || !hasNoTag))
        || (tagVal ==  0 && (nTypes != 1 || !hasNoTag)))
      return kFALSE;
  }

  return kTRUE;
}


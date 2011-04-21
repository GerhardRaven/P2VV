/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                *
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
// <p>RooMultiMultinomial is a binned function (not a PDF) with an arbitrary
// number of binnings. Optionally, the bins are associated with values of
// continuous variables (which results in a step function).</p>
//
// <p>A coefficient is specified by the user for each bin. The function's value
// is determined with the value of a particular bin's coefficient. In case of
// multiple binnings, the function's value is the product of the values of the
// binnings. The coefficients are either interpreted as bin heights or as bin
// integrals.</p>
//
// <p>The user may either specify coefficients for all the bins of a binning or
// leave it to the function to calculate the value of the last bin. In the
// former case, the coefficients may have any value. In the latter case, the
// value of the last bin's coefficient is equal to one minus the sum of the
// other bin coefficients. Negative bin values are then considered to be zero.
// If the sum of the specified bin coefficients is larger than one, the
// coefficient of the last bin is zero and the other coefficients are scaled
// such that their sum yields one.</p>
// END_HTML
//

#include "RooAbsCategory.h"
#include "RooAbsRealLValue.h"
#include "RooBinningCategory.h"
#include "RooArgSet.h"
#include "RooMsgService.h"
#include "RooMultiMultinomial.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"

ClassImp(RooMultiMultinomial);

//_____________________________________________________________________________
RooMultiMultinomial::RooMultiMultinomial(const char* name, const char* title,
    RooAbsCategory& baseCat, RooArgList& coefList, Bool_t ignoreFirstBin) :
  RooAbsReal(name, title),
  _numCats(1),
  _baseCatsList(TString(name) + "_baseCatsList", 0, this),
  _baseVarsList(TString(name) + "_baseVarsList", 0, this),
  _coefLists(1, 0),
  _continuousBase(kFALSE),
  _binIntegralCoefs(kTRUE),
  _ignoreFirstBin(ignoreFirstBin),
  _calcLastCoefs(kFALSE)
{
  // constructor with one binning, which depends on the value of a category
  //
  // The "current bin" is given by the value of the RooAbsCategory "baseCat".
  // The number of bins is given by the number of category types (N).
  //
  // Bin coefficients are specified in the list "coefList" and inherit from
  // RooAbsReal. One may either specify N or N - 1 coefficients (see also the
  // class description). The order of the coefficients is given by the indices
  // of "baseCat".
  //
  // If "ignoreFirstBin" is true, evaluation of the value of the first bin
  // always yields zero.

  // add category to list of base categories
  _baseCatsList.add(baseCat);

  // initialize coefficients
  TObjArray coefLists(1, 0);
  coefLists.Add(&coefList);
  initCoefs(coefLists);
}

//_____________________________________________________________________________
RooMultiMultinomial::RooMultiMultinomial(const char* name, const char* title,
    const RooArgList& baseCats, const TObjArray& coefLists,
    Bool_t ignoreFirstBin) :
  RooAbsReal(name, title),
  _numCats(0),
  _baseCatsList(TString(name) + "_baseCatsList", 0, this),
  _baseVarsList(TString(name) + "_baseVarsList", 0, this),
  _coefLists(1, 0),
  _continuousBase(kFALSE),
  _binIntegralCoefs(kTRUE),
  _ignoreFirstBin(ignoreFirstBin),
  _calcLastCoefs(kFALSE)
{
  // constructor with an arbitrary number of binnings, which depend on the
  // values of an equal number of categories
  //
  // The "current bins" are given by the values of the RooAbsCategories
  // contained by "baseCats". Exactly one category is specified for each
  // binning. The number of bins in a binning is given by the number of
  // category types (N_i).
  //
  // Bin coefficients are specified in the RooArgLists contained by coefLists.
  // The coefficients inherit from RooAbsReal. One may either specify N_i or
  // N_i - 1 coefficients (see also the class description). The order of the
  // coefficients is given by the indices of the "baseCats".
  //
  // If "ignoreFirstBin" is true, evaluation of the function value always
  // yields zero in the case that all the "current bins" are the first bin of
  // the binnings.

  // loop over input categories
  TIterator* catIter = baseCats.createIterator();
  RooAbsArg* cat = 0;
  while ((cat = (RooAbsArg*)catIter->Next()) != 0) {
    // check if this is a category
    if (dynamic_cast<RooAbsCategory*>(cat) == 0) {
      coutF(InputArguments) << "RooMultiMultinomial::RooMultiMultinomial("
          << GetName() << ") base category '" << cat->GetName()
          << "' is not a RooAbsCategory" << endl;
      reset();
      break;
    }

    // add category to list of base categories
    ++_numCats;
    _baseCatsList.add(*cat);
  }

  delete catIter;

  // initialize coefficients
  if (_numCats > 0) initCoefs(coefLists);
}

//_____________________________________________________________________________
RooMultiMultinomial::RooMultiMultinomial(const char* name, const char* title,
    RooAbsRealLValue& baseVar, const char* binningName,
    RooArgList& coefList, Bool_t binIntegralCoefs, Bool_t ignoreFirstBin) :
  RooAbsReal(name, title),
  _numCats(0),
  _baseCatsList(TString(name) + "_baseCatsList", 0, this),
  _baseVarsList(TString(name) + "_baseVarsList", 0, this),
  _coefLists(1, 0),
  _continuousBase(kTRUE),
  _binIntegralCoefs(binIntegralCoefs),
  _ignoreFirstBin(ignoreFirstBin),
  _calcLastCoefs(kFALSE)
{
  // constructor with one binning, which depends on the value of a continuous
  // variable with a binning
  //
  // The "current bin" is given by the value of the RooAbsRealLValue "baseVar"
  // and its binning with N bins, given by "binningName".
  //
  // Bin coefficients are specified in the list "coefList" and inherit from
  // RooAbsReal. One may either specify N or N - 1 coefficients (see also the
  // class description). The order of the coefficients is given by the order of
  // the "baseVar" bins. The coefficients may either be interpreted as bin
  // heights or as bin integrals. This is controlled by "binIntegralCoefs".
  //
  // If "ignoreFirstBin" is true, evaluation of the value of the first bin
  // always yields zero.

  // put the base variable in a list
  RooArgList baseVars(baseVar);

  // put the binning name in a list
  TObjArray binningNames(1, 0);
  TObjString strName(binningName);
  binningNames.Add(&strName);

  // create base category and initialize coefficients
  if (createBaseCats(baseVars, binningNames) > 0) {
    TObjArray coefLists(1, 0);
    coefLists.Add(&coefList);
    initCoefs(coefLists);
  }
}

//_____________________________________________________________________________
RooMultiMultinomial::RooMultiMultinomial(const char* name, const char* title,
    const RooArgList& baseVars, const TObjArray& binningNames,
    const TObjArray& coefLists, Bool_t binIntegralCoefs,
    Bool_t ignoreFirstBin) :
  RooAbsReal(name, title),
  _numCats(0),
  _baseCatsList(TString(name) + "_baseCatsList", 0, this),
  _baseVarsList(TString(name) + "_baseVarsList", 0, this),
  _coefLists(1, 0),
  _continuousBase(kTRUE),
  _binIntegralCoefs(binIntegralCoefs),
  _ignoreFirstBin(ignoreFirstBin),
  _calcLastCoefs(kFALSE)
{
  // constructor with an arbitrary number of binnings, which depend on the
  // values of continuous variables with binnings
  //
  // The "current bins" are given by the values of the RooAbsRealLValues
  // contained by "baseVars" and their binnings with N_i bins, given by
  // "binningNames". Exactly one variable is specified for each binning. The
  // "binningNames" are TObjStrings.
  //
  // Bin coefficients are specified in the lists contained by "coefLists". The
  // coefficients inherit from RooAbsReal. One may either specify N_i or
  // N_i - 1 coefficients (see also the class description). The order of the
  // coefficients is given by the order of the "baseVars" bins. The
  // coefficients may either be interpreted as bin heights or as bin integrals.
  // This is controlled by "binIntegralCoefs".
  //
  // If "ignoreFirstBin" is true, evaluation of the function value always
  // yields zero in the case that all the "current bins" are the first bin of
  // the binnings.

  // create base categories and initialize coefficients
  if (createBaseCats(baseVars, binningNames) > 0) initCoefs(coefLists);
}

//_____________________________________________________________________________
RooMultiMultinomial::RooMultiMultinomial(const RooMultiMultinomial& other,
    const char* name) :
  RooAbsReal(other, name),
  _numCats(other._numCats),
  _baseCatsList(TString(name) + "_baseCatsList", this, other._baseCatsList),
  _baseVarsList(TString(name) + "_baseVarsList", this, other._baseVarsList),
  _coefLists(_numCats, 0),
  _missingCoefs(other._missingCoefs),
  _indexPositions(other._indexPositions),
  _binningNames(other._binningNames),
  _continuousBase(other._continuousBase),
  _binIntegralCoefs(other._binIntegralCoefs),
  _ignoreFirstBin(other._ignoreFirstBin),
  _calcLastCoefs(other._calcLastCoefs)
{
  // copy constructor

  // make coefficient lists array owner of its lists
  _coefLists.SetOwner(kTRUE);

  // copy coefficient lists
  for (Int_t cListIter = 0; cListIter < _numCats; ++cListIter) {
    // get other's coefficients list
    RooListProxy* cListOther
        = (RooListProxy*)other._coefLists.UncheckedAt(cListIter);

    // create new coefficients list
    TString cListName(GetName());
    cListName += "_coefList";
    cListName += cListIter;
    RooListProxy* cListProxy = new RooListProxy(cListName, this, *cListOther);
    _coefLists.Add(cListProxy);
  }
}

//_____________________________________________________________________________
RooMultiMultinomial::~RooMultiMultinomial()
{
  // destructor
}

//_____________________________________________________________________________
Double_t RooMultiMultinomial::evaluate() const
{
  // evaluates and returns the current value of the function

  // temporary result of evaluation
  Double_t value = 1.;

  // temporary array of coefficient positions
  Int_t* cPos = new Int_t [_numCats];

  // loop over base categories
  Bool_t ignoreBin = _ignoreFirstBin;
  for (Int_t catIter = 0; catIter < _numCats; ++catIter) {
    // get coefficient position
    std::map<Int_t, Int_t> indexMap = _indexPositions[catIter];
    cPos[catIter] = indexMap[((RooAbsCategory*)_baseCatsList.at(catIter))
        ->getIndex()];
    if (cPos[catIter] != 0) ignoreBin = kFALSE;

    // multiply temporary result with value of coefficient
    if (_missingCoefs[catIter] != cPos[catIter]) {
      value *= ((RooAbsReal*)((RooArgList*)_coefLists.UncheckedAt(catIter))
          ->at(cPos[catIter]))->getVal();

      // divide by bin width if coefficients are bin integrals
      if (_continuousBase && _binIntegralCoefs) {
        value /= ((RooAbsRealLValue*)_baseVarsList.at(catIter))
            ->getBinning(_binningNames[catIter]).binWidth(cPos[catIter]);
      }
    }
  }

  // return value if we don't have to calculate last coefficients
  if (ignoreBin) return 0.;
  if (!_calcLastCoefs || !(TMath::Abs(value) > 0.)) return value;

  // loop over base categories again
  for (Int_t catIter = 0; catIter < _numCats; ++catIter) {
    if (_missingCoefs[catIter] >= 0) {
      // last coefficient is not explicitely specified for this category
      // * coefficients have values between zero and one
      // * sum of coefficients is equal to one
      // * last coefficient is given by one minus the sum of the other coefs.

      // loop over specified coefficients (#category types - 1)
      Double_t coefSum = 0.;
      RooArgList* coefList = (RooArgList*)_coefLists.UncheckedAt(catIter);
      for (Int_t coefIter = 0; coefIter < coefList->getSize(); ++coefIter) {
        // get coefficient's value
        Double_t coefVal = ((RooAbsReal*)coefList->at(coefIter))->getVal();

        // check for negative coefficients
        if (coefVal < 0) {
          // make negative values equal to zero
          if (coefIter == cPos[catIter]) return 0.;
        } else {
          // add positive values to sum
          coefSum += coefVal;
        }
      }

      if (_missingCoefs[catIter] == cPos[catIter]) {
        // multiply result with the calculated value of last coefficient
        if (coefSum > 1.) {
          return 0.;
        } else {
          value *= 1. - coefSum;

          // divide by bin width if coefficients are bin integrals
          if (_continuousBase && _binIntegralCoefs) {
            value /= ((RooAbsRealLValue*)_baseVarsList.at(catIter))
                ->getBinning(_binningNames[catIter]).binWidth(cPos[catIter]);
          }
        }
      } else if (coefSum > 1.) {
        // scale these coefficients to make their sum equal to one
        value /= coefSum;
      }
    }
  }

  delete[] cPos;

  return value;
}

//_____________________________________________________________________________
Int_t RooMultiMultinomial::createBaseCats(const RooArgList& baseVars,
    const TObjArray& binningNames)
{
  // creates a base category for each continuous input base variable

  _numCats = baseVars.getSize();

  // check the number of specified binnings
  if (binningNames.GetEntries() != _numCats) {
    coutF(InputArguments) << "RooMultiMultinomial::createBaseCats("
        << GetName() << ") number of binnings (" << binningNames.GetEntries()
        << ") does not match the number of base variables (" << _numCats << ")"
        << endl;
    reset();
    return -1;
  }

  // loop over base variables
  for (Int_t varIter = 0; varIter < _numCats; ++varIter) {
    // get input variable
    TString catName(baseVars.at(varIter)->GetName());
    catName += "_cat";
    RooAbsRealLValue* var
        = dynamic_cast<RooAbsRealLValue*>(baseVars.at(varIter));

    // check type of input variable
    if (var == 0) {
      coutF(InputArguments) << "RooMultiMultinomial::createBaseCats("
          << GetName() << ") base variable '" << var->GetName()
          << "' is not a RooAbsRealLValue" << endl;
      reset();
      return -2;
    }

    // check if list of base variables already contains this variable
    if (_baseVarsList.find(var->GetName()) != 0) {
      coutF(InputArguments) << "RooMultiMultinomial::createBaseCats("
          << GetName() << ") variable '" << var->GetName()
          << "' is not unique in set of base variables" << endl;
      reset();
      return -3;
    }

    // add variable to list of base variables
    _baseVarsList.add(*var);

    // create binning of input variable if necessary
    TString bins = ((const TObjString*)binningNames.At(varIter))->GetString();
    if (!var->hasBinning(bins))
      var->getBinning(bins, kTRUE, kTRUE);

    // add binning name to list
    _binningNames.push_back(bins);

    // create category
    RooBinningCategory* cat = new RooBinningCategory(catName, catName, *var,
        bins);
    _baseCatsList.addOwned(*cat);
  }

  return _baseCatsList.getSize();
}

//_____________________________________________________________________________
Int_t RooMultiMultinomial::initCoefs(const TObjArray& coefLists)
{
  // initializes bin coefficients

  // set size of coefficient lists array
  _coefLists.Expand(_numCats);

  // make array owner of the coefficient lists
  _coefLists.SetOwner(kTRUE);

  // loop over coefficient lists
  Int_t catPos = 0;
  TIterator* cListIter = coefLists.MakeIterator();
  RooArgList* cList = 0;
  while ((cList = dynamic_cast<RooArgList*>(cListIter->Next())) != 0) {
    // create new coefficients list
    TString cListName(GetName());
    cListName += "_coefList";
    cListName += catPos;
    RooListProxy* cListProxy = new RooListProxy(cListName, 0, this);
    _coefLists.Add(cListProxy);

    // loop over coefficients
    TIterator* coefIter = cList->createIterator();
    RooAbsReal* coef = 0;
    while ((coef = dynamic_cast<RooAbsReal*>(coefIter->Next())) != 0)
      cListProxy->add(*coef);

    delete coefIter;

    // get category corresponding to this list
    RooAbsCategory* cat = (RooAbsCategory*)_baseCatsList.at(catPos);

    // set/check number of coefficients
    if (cListProxy->getSize() == cat->numTypes()) {
      _missingCoefs.push_back(-1);
    } else if (cListProxy->getSize() == cat->numTypes() - 1) {
      _calcLastCoefs = kTRUE;
      _missingCoefs.push_back(cListProxy->getSize());
    } else {
      coutF(InputArguments) << "RooMultiMultinomial::initCoefs("
          << GetName() << ") number of coefficients (" << cListProxy->getSize()
          << ") does not match number of base variable types ("
          << cat->numTypes() << ")" << endl;
      reset();
      return -2;
    }

    // get category indices
    vector<Int_t> catIndices;
    RooCatType* catType = 0;
    TIterator* catTypeIter = cat->typeIterator();
    while ((catType = (RooCatType*)catTypeIter->Next()) != 0)
      catIndices.push_back(catType->getVal());
    delete catTypeIter;
    sort(catIndices.begin(), catIndices.end());

    // set position of category indices
    map<Int_t, Int_t> indexPosMap;
    for (Int_t indexPos = 0; indexPos < (Int_t)catIndices.size(); ++indexPos)
      indexPosMap[catIndices.at(indexPos)] = indexPos;
    _indexPositions.push_back(indexPosMap);

    ++catPos;
  }

  delete cListIter;

  // check number of coefficient lists
  if (_coefLists.GetEntries() != _numCats) {
    coutF(InputArguments) << "RooMultiMultinomial::initCoefs("
        << GetName() << ") number of coefficient lists ("
        << _coefLists.GetEntries()
        << ") does not match number of base variables (" << _numCats << ")"
        << endl;
    reset();
    return -2;
  }

  return _numCats;
}

//_____________________________________________________________________________
RooArgList* RooMultiMultinomial::baseVariables()
{
  // Returns the list of base variables (continuous variables or categories).
  // The caller of this function is responsible for deleting the returned
  // arglist.

  if (_continuousBase) return new RooArgList(_baseVarsList, "base_variables");
  else                 return new RooArgList(_baseCatsList, "base_variables");
}

//_____________________________________________________________________________
void RooMultiMultinomial::reset()
{
  _numCats = 0;
  _baseCatsList.removeAll();
  _baseVarsList.removeAll();
  _coefLists.Clear();
  _missingCoefs.clear();
  _indexPositions.clear();
  _binningNames.clear();
  _calcLastCoefs = kFALSE;
}


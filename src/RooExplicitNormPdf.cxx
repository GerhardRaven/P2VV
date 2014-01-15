/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                *
 *                                                                           *
 * Copyright (c) 2012, Nikhef. All rights reserved.                          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/


//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// PDF with explicitly specified normalization function. The function and the
// normalization function of the PDF are specified separately.
// END_HTML
//

#include <iostream>
using std::endl;

#include "Riostream.h"
#include "RooMsgService.h"
#include "TMath.h"
#include "P2VV/RooExplicitNormPdf.h"
#include "RooArgSet.h"
#include "RooCustomizer.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"

//_____________________________________________________________________________
RooExplicitNormPdf::RooExplicitNormPdf(const char *name, const char* title,
      const RooArgSet& obsSet, const RooAbsReal& function,
      const RooAbsReal& normFunc, Double_t normFactor) :
  RooAbsReal(name, title),
  _obsSet("obsSet", "observables set", this),
  _parSet("parSet", "parameters set", this),
  _functionOrig(&function),
  _normFuncOrig(&normFunc),
  _function(0),
  _normFunc(0),
  _functionInteg(0),
  _normFuncInteg(0),
  _normFactor(normFactor),
  _projData(0)
{
  initVariables(&obsSet, 0);
}

//_____________________________________________________________________________
RooExplicitNormPdf::RooExplicitNormPdf(const char *name, const char* title,
      const RooArgSet& obsSet, const RooAbsReal& function,
      const RooAbsReal& normFunc, Double_t normFactor,
      const RooAbsData& projectionData) :
  RooAbsReal(name, title),
  _obsSet("obsSet", "observables set", this),
  _parSet("parSet", "parameters set", this),
  _functionOrig(&function),
  _normFuncOrig(&normFunc),
  _function(0),
  _normFunc(0),
  _functionInteg(0),
  _normFuncInteg(0),
  _normFactor(normFactor),
  _projData(&projectionData)
{
  initVariables(&obsSet, 0);
}

//_____________________________________________________________________________
RooExplicitNormPdf::RooExplicitNormPdf(const char *name, const char* title,
      const RooArgSet& obsSet, const RooArgSet& intObsSet,
      const RooAbsReal& function, const RooAbsReal& normFunc,
      Double_t normFactor, const char* intRangeFunc,
      const char* intRangeNorm) :
  RooAbsReal(name, title),
  _obsSet("obsSet", "observables set", this),
  _parSet("parSet", "parameters set", this),
  _functionOrig(&function),
  _normFuncOrig(&normFunc),
  _function(0),
  _normFunc(0),
  _functionInteg(0),
  _normFuncInteg(0),
  _normFactor(normFactor),
  _projData(0),
  _intRangeFunc(intRangeFunc),
  _intRangeNorm(intRangeNorm)
{
  initVariables(&obsSet, &intObsSet);
}

//_____________________________________________________________________________
RooExplicitNormPdf::RooExplicitNormPdf(const char *name, const char* title,
      const RooArgSet& obsSet, const RooArgSet& intObsSet,
      const RooAbsReal& function, const RooAbsReal& normFunc,
      Double_t normFactor, const RooAbsData& projectionData,
      const char* intRangeFunc, const char* intRangeNorm) :
  RooAbsReal(name, title),
  _obsSet("obsSet", "observables set", this),
  _parSet("parSet", "parameters set", this),
  _functionOrig(&function),
  _normFuncOrig(&normFunc),
  _function(0),
  _normFunc(0),
  _functionInteg(0),
  _normFuncInteg(0),
  _normFactor(normFactor),
  _projData(&projectionData),
  _intRangeFunc(intRangeFunc),
  _intRangeNorm(intRangeNorm)
{
  initVariables(&obsSet, &intObsSet);
}

//_____________________________________________________________________________
RooExplicitNormPdf::RooExplicitNormPdf(const RooExplicitNormPdf& other,
      const char* name) :
  RooAbsReal(other, name),
  _obsSet("obsSet", this, other._obsSet),
  _intObsSet(other._intObsSet),
  _parSet("parSet", this, other._parSet),
  _functionOrig(other._functionOrig),
  _normFuncOrig(other._normFuncOrig),
  _function(0),
  _normFunc(0),
  _functionInteg(0),
  _normFuncInteg(0),
  _normFactor(other._normFactor),
  _projData(other._projData),
  _intRangeFunc(other._intRangeFunc),
  _intRangeNorm(other._intRangeNorm)
{}

//_____________________________________________________________________________
RooExplicitNormPdf::~RooExplicitNormPdf()
{
  if (_function != 0 && _function != _functionOrig) delete _function;
  if (_normFunc != 0 && _normFunc != _normFuncOrig) delete _normFunc;
}

//_____________________________________________________________________________
Double_t RooExplicitNormPdf::evaluate() const
{
  if (_function == 0) initFunctions();

  Double_t intVal = 0.;
  if (_projData == 0) {
    // get PDF integral value
    intVal = _function->getVal() / _normFunc->getVal();
  } else {
    // get data-weighted average for PDF integral with projection data set
    for (Int_t dataIt = 0; dataIt < _projData->numEntries(); ++dataIt) {
      _projData->get(dataIt);
      intVal += _projData->weight() * _function->getVal() / _normFunc->getVal();
    }
    intVal /= _projData->sumEntries();
    std::cout << "." << std::flush;
  }

  return _normFactor * intVal;
}

//_____________________________________________________________________________
void RooExplicitNormPdf::initVariables(const RooArgSet* obsSet,
    const RooArgSet* intObsSet)
{
  // get variables in functions
  RooArgSet* functionVars = _functionOrig->getVariables();
  RooArgSet* normFuncVars = _normFuncOrig->getVariables();

  // set observables
  RooAbsArg* var = 0;
  RooFIter varIter = obsSet->fwdIterator();
  while ((var = (RooAbsArg*)varIter.next()) != 0) {
    if (functionVars->find(var->GetName()) == 0) {
      coutW(InputArguments) << "RooExplicitNormPdf::initVariables("
          << GetName()
          << "): specified function does not depend on observable \""
          << var->GetName() << "\"" << std::endl;
    }

    _obsSet.add(*var);
  }

  if (intObsSet != 0) {
    // set integration observables
    varIter = intObsSet->fwdIterator();
    while ((var = (RooAbsArg*)varIter.next()) != 0) {
      if (_obsSet.find(var->GetName()) != 0) {
        coutW(InputArguments) << "RooExplicitNormPdf::initVariables("
            << GetName() << "): integration observable \"" << var->GetName()
            << "\" is in the set of observables; not integrating over \""
            << var->GetName() << "\"" << std::endl;
        continue;
      }

      _intObsSet.add(*var);
    }
  }

  // set parameters
  RooArgSet parSet(*functionVars);
  parSet.add(*normFuncVars, kTRUE);
  parSet.remove(_obsSet, kFALSE, kTRUE);
  parSet.remove(_intObsSet, kFALSE, kTRUE);
  if (_projData) parSet.remove(*_projData->get(), kFALSE, kTRUE);
  _parSet.add(parSet);
}

//_____________________________________________________________________________
void RooExplicitNormPdf::initFunctions() const
{
  _functionOrig->getVal();
  _normFuncOrig->getVal();
  if (_function != 0 && _function != _functionOrig) {
    delete _function;
    _function = 0;
  }
  if (_normFunc != 0 && _normFunc != _normFuncOrig) {
    delete _normFunc;
    _normFunc = 0;
  }
  if (_functionInteg != 0 && _functionInteg != _functionOrig) {
    delete _functionInteg;
    _functionInteg = 0;
  }
  if (_normFuncInteg != 0 && _normFuncInteg != _normFuncOrig) {
    delete _normFuncInteg;
    _normFuncInteg = 0;
  }

  // create customizers for the function and normalization function
  RooCustomizer funcCust(*_functionOrig, "function");
  RooCustomizer normCust(*_normFuncOrig, "normFunc");

  RooArgSet* funcObsSet = _functionOrig->getVariables();
  RooArgSet* normObsSet = _normFuncOrig->getVariables();
  RooAbsArg* obs = 0;

  // replace observables in functions by observables in our set
  RooFIter obsIter = _obsSet.fwdIterator();
  while ((obs = (RooAbsArg*)obsIter.next()) != 0) {
    RooAbsArg* obsFunc = funcObsSet->find(obs->GetName());
    if (obsFunc != 0) funcCust.replaceArg(*obsFunc, *obs);
    obsFunc = normObsSet->find(obs->GetName());
    if (obsFunc != 0) normCust.replaceArg(*obsFunc, *obs);
  }

  if (_projData) {
    // replace observables in functions by observables in projection data set
    obsIter = _projData->get()->fwdIterator();
    while ((obs = (RooAbsArg*)obsIter.next()) != 0) {
      if (_obsSet.find(obs->GetName()) != 0) continue;

      RooAbsArg* obsFunc = funcObsSet->find(obs->GetName());
      if (obsFunc != 0) funcCust.replaceArg(*obsFunc, *obs);
      obsFunc = normObsSet->find(obs->GetName());
      if (obsFunc != 0) normCust.replaceArg(*obsFunc, *obs);
    }

    coutI(Eval) << "RooExplicitNormPdf::initFunctions(" << GetName()
        << "): evaluating integral as a data-weighted average with data set \""
        << _projData->GetName() << "\": " << _projData->numEntries() << " ("
        << _projData->sumEntries() << ") entries" << std::endl;
    coutI(Eval) << "conditional observables: ";
    _projData->get()->Print();
  }

  delete funcObsSet;
  delete normObsSet;

  // build functions with replaced observables
  _functionInteg = (RooAbsReal*)funcCust.build();
  _normFuncInteg = (RooAbsReal*)normCust.build();
  funcObsSet = _functionInteg->getVariables();
  normObsSet = _normFuncInteg->getVariables();

  // integrate function
  RooArgSet funcIntSet;
  obsIter = _intObsSet.fwdIterator();
  while ((obs = (RooAbsArg*)obsIter.next()) != 0) {
    RooAbsArg* obsFunc = funcObsSet->find(obs->GetName());
    if (obsFunc != 0) funcIntSet.add(*obsFunc);
  }
  if (funcIntSet.getSize() > 0) {
    _function = _functionInteg->createIntegral(funcIntSet, 0, 0,
        (_intRangeFunc.Length() > 0 ? _intRangeFunc.Data() : 0));
  } else {
    _function = _functionInteg;
  }

  // integrate normalization function
  RooArgSet normIntSet;
  obsIter = _obsSet.fwdIterator();
  while ((obs = (RooAbsArg*)obsIter.next()) != 0) {
    RooAbsArg* obsFunc = normObsSet->find(obs->GetName());
    if (obsFunc != 0) normIntSet.add(*obsFunc);
  }
  obsIter = _intObsSet.fwdIterator();
  while ((obs = (RooAbsArg*)obsIter.next()) != 0) {
    RooAbsArg* obsFunc = normObsSet->find(obs->GetName());
    if (obsFunc != 0) normIntSet.add(*obsFunc);
  }
  if (normIntSet.getSize() > 0) {
    _normFunc = _normFuncInteg->createIntegral(normIntSet, 0, 0,
        (_intRangeNorm.Length() > 0 ? _intRangeNorm.Data() : 0));
  } else {
    _normFunc = _normFuncInteg;
  }

  delete funcObsSet;
  delete normObsSet;

  // print functions
  _function->getVal();
  _normFunc->getVal();
  coutI(Eval) << "RooExplicitNormPdf::initFunctions(" << GetName()
      << "): function:" << std::endl;
  _function->Print();
  coutI(Eval) << "RooExplicitNormPdf::initFunctions(" << GetName()
      << "): normalization function:" << std::endl;
  _normFunc->Print();
}

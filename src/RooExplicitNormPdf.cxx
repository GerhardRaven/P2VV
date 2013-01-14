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
#include "RooExplicitNormPdf.h"
#include "RooArgSet.h"
#include "RooCustomizer.h"
#include "RooRealVar.h"

ClassImp(RooExplicitNormPdf);

//_____________________________________________________________________________
RooExplicitNormPdf::RooExplicitNormPdf(const char *name, const char* title,
      const RooArgSet& obsSet, const RooAbsReal& function,
      const RooAbsReal& normFunc, Double_t normFactor) :
  RooAbsReal(name, title),
  _obsSet("obsSet", "observables set", this),
  _functionOrig(&function),
  _normFuncOrig(&normFunc),
  _function(0),
  _normFunc(0),
  _normFactor(normFactor),
  _projData(0)
{
  _obsSet.add(obsSet);
}

//_____________________________________________________________________________
RooExplicitNormPdf::RooExplicitNormPdf(const char *name, const char* title,
      const RooArgSet& obsSet, const RooAbsReal& function,
      const RooAbsReal& normFunc, Double_t normFactor,
      const RooAbsData& projectionData) :
  RooAbsReal(name, title),
  _obsSet("obsSet", "observables set", this),
  _functionOrig(&function),
  _normFuncOrig(&normFunc),
  _function(0),
  _normFunc(0),
  _normFactor(normFactor),
  _projData(&projectionData)
{
  _obsSet.add(obsSet);
}

//_____________________________________________________________________________
RooExplicitNormPdf::RooExplicitNormPdf(const RooExplicitNormPdf& other,
      const char* name) :
  RooAbsReal(other, name),
  _obsSet("obsSet", this, other._obsSet),
  _functionOrig(other._functionOrig),
  _normFuncOrig(other._normFuncOrig),
  _function(0),
  _normFunc(0),
  _normFactor(other._normFactor),
  _projData(other._projData)
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
    cout << "." << flush;
  }

  return _normFactor * intVal;
}

//_____________________________________________________________________________
void RooExplicitNormPdf::initFunctions() const
{
  if (_function != 0 && _function != _functionOrig) delete _function;
  if (_normFunc != 0 && _normFunc != _normFuncOrig) delete _normFunc;
  _functionOrig->getVal();
  _normFuncOrig->getVal();

  // create customizers for the function and normalization function
  RooCustomizer funcCust(*_functionOrig, "function");
  RooCustomizer normCust(*_normFuncOrig, "normFunc");

  RooArgSet* funcObsSet = _functionOrig->getVariables();
  RooArgSet* normObsSet = _normFuncOrig->getVariables();
  RooAbsArg* obs = 0;

  // replace observables in functions by observables in our set
  TIterator* obsIter = _obsSet.createIterator();
  while ((obs = (RooAbsArg*)obsIter->Next()) != 0) {
    RooAbsArg* obsFunc = funcObsSet->find(obs->GetName());
    if (obsFunc != 0) funcCust.replaceArg(*obsFunc, *obs);
    obsFunc = normObsSet->find(obs->GetName());
    if (obsFunc != 0) normCust.replaceArg(*obsFunc, *obs);
  }
  delete obsIter;

  if (_projData) {
    // replace observables in functions by observables in projection data set
    obsIter = _projData->get()->createIterator();
    while ((obs = (RooAbsArg*)obsIter->Next()) != 0) {
      if (_obsSet.find(obs->GetName()) != 0) continue;

      RooAbsArg* obsFunc = funcObsSet->find(obs->GetName());
      if (obsFunc != 0) funcCust.replaceArg(*obsFunc, *obs);
      obsFunc = normObsSet->find(obs->GetName());
      if (obsFunc != 0) normCust.replaceArg(*obsFunc, *obs);
    }
    delete obsIter;

    coutI(Eval) << "RooExplicitNormPdf::initFunctions(" << GetName()
        << ") evaluating integral as a data-weighted average with data set \""
        << _projData->GetName() << "\": " << setprecision(1) << fixed
        << _projData->numEntries() << " (" << _projData->sumEntries()
        << ") entries" << endl;
    coutI(Eval) << "conditional observables:" << endl;
    _projData->get()->Print();
  }

  delete funcObsSet;
  delete normObsSet;

  // build functions with replaced observables
  _function = (RooAbsReal*)funcCust.build();
  _normFunc = (RooAbsReal*)normCust.build();
}

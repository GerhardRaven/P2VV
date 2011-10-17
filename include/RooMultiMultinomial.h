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

#ifndef ROO_MULTI_MULTINOMIAL
#define ROO_MULTI_MULTINOMIAL

#include <map>
#include <vector>

#include "RooAbsReal.h"
#include "RooListProxy.h"

class RooAbsCategory;
class RooAbsRealLValue;
class RooArgSet;
class TObjArray;

class RooMultiMultinomial : public RooAbsReal
{

public:
  inline RooMultiMultinomial() {}

  RooMultiMultinomial(const char *name, const char *title,
      RooAbsCategory& baseCat, RooArgList& coefList,
      Bool_t ignoreFirstBin = kFALSE);

  RooMultiMultinomial(const char *name, const char *title,
      const RooArgList& baseCats, const TObjArray& coefLists,
      Bool_t ignoreFirstBin = kFALSE);

  RooMultiMultinomial(const char *name, const char *title,
      RooAbsRealLValue& baseVar, const char* binningName,
      RooArgList& coefList, Bool_t binIntegralCoefs = kFALSE,
      Bool_t ignoreFirstBin = kFALSE);

  RooMultiMultinomial(const char *name, const char *title,
      const RooArgList& baseVars, const TObjArray& binningNames,
      const TObjArray& coefLists, Bool_t binIntegralCoefs = kFALSE,
      Bool_t ignoreFirstBin = kFALSE);

  RooMultiMultinomial(const RooMultiMultinomial& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const 
  { 
    return new RooMultiMultinomial(*this, newname);
  }

  virtual ~RooMultiMultinomial();

  RooArgList* baseVariables();

  Bool_t continuousBase() {return _continuousBase;}

  Bool_t binIntegralCoefs() {return _binIntegralCoefs;}
  void   setBinIntegralCoefs(Bool_t integralCoefs = kTRUE)
  {
    _binIntegralCoefs = integralCoefs;
  }

  Bool_t ignoreFirstBin() {return _ignoreFirstBin;}

private:
  virtual Double_t evaluate() const;

  Int_t createBaseCats(const RooArgList& baseVars,
      const TObjArray& binningNames);
  Int_t initCoefs(const TObjArray& coefLists);

  void reset();

  Int_t _numCats;

  RooListProxy _baseCatsList;
  RooListProxy _baseVarsList;
  TObjArray    _coefLists;

  std::vector<Int_t>                    _missingCoefs;
  std::vector< std::map<Int_t, Int_t> > _indexPositions;
  std::vector<TString>                  _binningNames;

  Bool_t _continuousBase;
  Bool_t _binIntegralCoefs;
  Bool_t _ignoreFirstBin;
  Bool_t _calcLastCoefs;

  ClassDef(RooMultiMultinomial, 1) // multi-multinomial function
};

#endif


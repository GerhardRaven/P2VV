/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                *
 *   RA,  Roel Aaij,          Nikhef                                         *
 *                                                                           *
 * Copyright (c) 2011, Nikhef. All rights reserved.                          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROO_BINNED_PDF
#define ROO_BINNED_PDF

#include <map>
#include <vector>

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooAbsCategory;
class RooAbsRealLValue;
class RooArgSet;
class TObjArray;

class RooBinnedPdf : public RooAbsPdf
{

public:
  inline RooBinnedPdf() {}

  RooBinnedPdf(const char *name, const char *title,
      RooAbsCategory& baseCat, const RooArgList& coefList,
      Bool_t ignoreFirstBin = kFALSE);

  RooBinnedPdf(const char *name, const char *title,
      const RooArgList& baseCats, const TObjArray& coefLists,
      Bool_t ignoreFirstBin = kFALSE);

  RooBinnedPdf(const char *name, const char *title,
      RooAbsRealLValue& baseVar, const char* binningName,
      const RooArgList& coefList, Bool_t binIntegralCoefs = kFALSE,
      Bool_t ignoreFirstBin = kFALSE);

  RooBinnedPdf(const char *name, const char *title,
      const RooArgList& baseVars, const TObjArray& binningNames,
      const TObjArray& coefLists, Bool_t binIntegralCoefs = kFALSE,
      Bool_t ignoreFirstBin = kFALSE);

  RooBinnedPdf(const char* name, const char* title,
      const RooArgList& baseVars, const TObjArray& binningNames,
      RooAbsReal* function);

  RooBinnedPdf(const char* name, const char* title,
      const RooAbsArg& baseVar, const char* binning,
      RooAbsReal* function);

  RooBinnedPdf(const RooBinnedPdf& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const 
  { 
    return new RooBinnedPdf(*this, newname);
  }

  virtual ~RooBinnedPdf();

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
   
  const RooAbsRealLValue& function() {
    return dynamic_cast<const RooAbsRealLValue&>(_function.arg());
  }

  Double_t evaluateCoef() const;
  Double_t evaluateFunction() const;

  Int_t _numCats;

  RooListProxy _baseCatsList;
  RooListProxy _baseVarsList;
  TObjArray    _coefLists;
  RooRealProxy _function;

  std::vector< std::map<Int_t, Int_t> > _indexPositions;
  std::vector<TString>                  _binningNames;
  std::vector<Bool_t>                   _calcCoefZeros;

  Bool_t _continuousBase;
  Bool_t _binIntegralCoefs;
  Bool_t _ignoreFirstBin;

  ClassDef(RooBinnedPdf, 1) // multi-multinomial function
};

#endif


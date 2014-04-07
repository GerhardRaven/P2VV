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

#ifndef ROO_EXPLICIT_NORM_PDF
#define ROO_EXPLICIT_NORM_PDF

#include "TString.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooSetProxy.h"
#include "RooAbsData.h"

class RooArgSet;

class RooExplicitNormPdf : public RooAbsReal
{

public:
  RooExplicitNormPdf() {};

  RooExplicitNormPdf(const char *name, const char *title,
      const RooArgSet& obsSet, const RooAbsReal& function,
      const RooAbsReal& normFunc, Double_t normFactor);

  RooExplicitNormPdf(const char *name, const char *title,
      const RooArgSet& obsSet, const RooAbsReal& function,
      const RooAbsReal& normFunc, Double_t normFactor,
      const RooAbsData& projectionData);

  RooExplicitNormPdf(const char *name, const char *title,
      const RooArgSet& obsSet, const RooArgSet& intObsSet,
      const RooAbsReal& function, const RooAbsReal& normFunc,
      Double_t normFactor, const char* intRangeFunc = 0,
      const char* intRangeNorm = 0);

  RooExplicitNormPdf(const char *name, const char *title,
      const RooArgSet& obsSet, const RooArgSet& intObsSet,
      const RooAbsReal& function, const RooAbsReal& normFunc,
      Double_t normFactor, const RooAbsData& projectionData,
      const char* intRangeFunc = 0, const char* intRangeNorm = 0);

  RooExplicitNormPdf(const RooExplicitNormPdf& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const 
  { 
    return new RooExplicitNormPdf(*this, newname);
  }

  virtual ~RooExplicitNormPdf();

  RooArgSet* obsSet() const
  {
    return new RooArgSet(_obsSet, "observables");
  }
  RooArgSet* intObsSet() const
  {
    return new RooArgSet(_intObsSet, "integrationObservables");
  }
  RooArgSet* parSet() const
  {
    return new RooArgSet(_parSet, "parameters");
  }

  const RooAbsReal* function() const {return _function;}
  const RooAbsReal* normFunc() const {return _normFunc;}
  Double_t normFactor()        const {return _normFactor;}
  const RooAbsData* projData() const {return _projData;}

  void initFunctions() const;

protected:
  RooSetProxy _obsSet;
  RooArgSet   _intObsSet;
  RooSetProxy _parSet;
  const RooAbsReal* _functionOrig;
  const RooAbsReal* _normFuncOrig;
  mutable RooArgSet* _functionClones; //!
  mutable RooArgSet* _normFuncClones; //!
  mutable RooAbsReal* _function; //!
  mutable RooAbsReal* _normFunc; //!
  Double_t _normFactor;
  const RooAbsData* _projData;
  TString _intRangeFunc;
  TString _intRangeNorm;

  Double_t evaluate() const;
  void initVariables(const RooArgSet* obsSet, const RooArgSet* intObsSet);

private:
  ClassDef(RooExplicitNormPdf, 1) // PDF with explicitly specified normalization function
};

#endif


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

#include "RooAbsReal.h"
#include "RooListProxy.h"
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

  RooExplicitNormPdf(const RooExplicitNormPdf& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const 
  { 
    return new RooExplicitNormPdf(*this, newname);
  }

  virtual ~RooExplicitNormPdf();

protected:
  RooListProxy _obsSet;
  const RooAbsReal* _functionOrig;
  const RooAbsReal* _normFuncOrig;
  mutable RooAbsReal* _function; //!
  mutable RooAbsReal* _normFunc; //!
  Double_t _normFactor;
  const RooAbsData* _projData;

  Double_t evaluate() const;
  void initFunctions() const;

private:
  ClassDef(RooExplicitNormPdf, 1) // PDF with explicitly specified normalization function
};

#endif


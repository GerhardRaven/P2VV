/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
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

#include "RooRealProxy.h"
#include "RooCategoryProxy.h"

class RooMultiMultinomial : public RooAbsReal
{

public:
  inline RooMultiMultinomial() {}

  RooMultiMultinomial(const char *name, const char *title);

  RooMultiMultinomial(const RooMultiMultinomial& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const 
  { 
    return new RooMultiMultinomial(*this, newname);
  }

  virtual ~RooMultiMultinomial();

  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
      const char* rangeName=0) const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName = 0)
      const;

protected:
  virtual Double_t evaluate() const;

  ClassDef(RooMultiMultinomial, 1) // 
};

#endif


/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   GR, Gerhard Raven,   Nikhef & VU, Gerhard.Raven@nikhef.nl
 *                                                                           *
 * Copyright (c) 2010, Nikhef & VU. All rights reserved.
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_P2VVAngleBasis
#define ROO_P2VVAngleBasis

#include "RooRealVar.h"
#include "RooProduct.h"

class RooP2VVAngleBasis : public RooProduct {
public:
  RooP2VVAngleBasis();
  RooP2VVAngleBasis(const char *name, const char *title,
      RooAbsReal& cpsi, RooAbsReal& ctheta, RooAbsReal& phi,
      Int_t i, Int_t j, Int_t l, Int_t m, Double_t c = 1. );
  RooP2VVAngleBasis(const char *name, const char *title,
      RooAbsReal& cpsi, RooAbsReal& ctheta, RooAbsReal& phi,
      Int_t i1, Int_t j1, Int_t l1, Int_t m1,
      Int_t i2, Int_t j2, Int_t l2, Int_t m2, Double_t c = 1.);
  RooP2VVAngleBasis(const RooP2VVAngleBasis& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const
  {
    return new RooP2VVAngleBasis(*this, newname);
  }

  inline virtual ~RooP2VVAngleBasis() { }

  // create a new RooAbsReal which is 'us' multiplied by an efficiency factor
  // Note: we can only multiply once... 
  RooP2VVAngleBasis* createProduct(Int_t i, Int_t j, Int_t l, Int_t m,
      Double_t c) const;
  RooP2VVAngleBasis* createProduct(const RooP2VVAngleBasis& basis,
      Double_t c = 1.) const
  { 
    // TODO: verify that basis and *this use the same angular observables!!!
    return !basis._prod ? createProduct(basis._i, basis._j, basis._l,
        basis._m, c * basis._c) : 0;
  }

  Int_t i() const { return _i; }
  Int_t j() const { return _j; }
  Int_t l() const { return _l; }
  Int_t m() const { return _m; }
  Double_t c() const { return _c; }
  bool prod() const { return _prod; }

private: 
  Double_t _c;
  Int_t _i,_j,_l,_m;
  bool _prod;

  ClassDef(RooP2VVAngleBasis, 1) // 
};

#endif

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
  RooP2VVAngleBasis() ;
  RooP2VVAngleBasis(const char *name, const char *title, RooAbsReal& cpsi, RooAbsReal& ctheta, RooAbsReal& phi, int i, int j, int l, int m, double c = 1. );
  RooP2VVAngleBasis(const char *name, const char *title, RooAbsReal& cpsi, RooAbsReal& ctheta, RooAbsReal& phi, int i1, int j1, int l1, int m1
                                                                                                              , int i2, int j2, int l2, int m2
                                                                                                              , double c = 1. );

  RooP2VVAngleBasis(const RooP2VVAngleBasis& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooP2VVAngleBasis(*this, newname); }
  inline virtual ~RooP2VVAngleBasis() { }

  // create a new RooAbsReal which is 'us' multiplied by an efficiency factor
  // Note: we can only multiply once... 
  RooP2VVAngleBasis* createProduct(int i, int j, int l, int m, double c) const;
  RooP2VVAngleBasis* createProduct(const RooP2VVAngleBasis& basis, double c=1.) const
  { 
      // TODO: verify that basis and *this use the same angular observables!!!
      return !basis._prod ? createProduct( basis._i, basis._j, basis._l, basis._m, c*basis._c ): 0;
  }

  int i() const { return _i; }
  int j() const { return _j; }
  int l() const { return _l; }
  int m() const { return _m; }

private: 
  double _c;
  int _i,_j,_l,_m;
  bool _prod;

  ClassDef(RooP2VVAngleBasis,1) // 
};

#endif

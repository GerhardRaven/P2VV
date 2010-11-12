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
  RooP2VVAngleBasis(const char *name, const char *title, RooRealVar& cpsi, RooRealVar& ctheta, RooRealVar& phi, int i, int l, int m, double c);

  RooP2VVAngleBasis(const P2VVAngleBasis& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooP2VVAngleBasis(*this, newname); }
  inline virtual ~RooP2VVAngleBasis() { }

  // create a new RooAbsReal which is 'us' multiplied by an efficiency factor
  RooP2VVAngleBasis* createProduct(int i, int j, int k, double eps_ijk) {
        return 0; //TODO: implement
  }

private: 
  int _i,_j,_l,_m;

  ClassDef(RooP2VVAngleBasis,1) // 
};

#endif

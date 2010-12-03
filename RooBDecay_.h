/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   GR, Gerhard Raven,  Vrije Universiteit Amsterdam & Nikhef               *
 *                                                                           *
 * Copyright (c) 2010, Vrije Universiteit Amsterdam & Nikhef                 *
 *                     All rights reserved.                                  *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROO_BDECAY_
#define ROO_BDECAY_

#include "RooBDecay.h"
#include "RooMultiCatIter.h"

class RooBDecay_ : public RooBDecay
{

public:
  inline RooBDecay_() { }
  RooBDecay_(const char *name, const char *title, RooRealVar& t,
			RooAbsReal& tau, RooAbsReal& dgamma,
			RooAbsReal& f0,
			RooAbsReal& f1, RooAbsReal& f2, 
			RooAbsReal& f3, RooAbsReal& dm, 
			const RooResolutionModel& model,
			DecayType type);
  RooBDecay_(const RooBDecay_& other, const char* name=0);
  virtual TObject* clone(const char* newname) const 
  { 
    return new RooBDecay_(*this,newname);
  }
  virtual ~RooBDecay_();

  Int_t getMaxVal(const RooArgSet& vars) const;
  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  void initGenerator(Int_t code);
  void generateEvent(Int_t code);

protected:

  class CacheElem : public RooAbsCacheElement {
  public :
        virtual ~CacheElem();
        // payload
        double *n;
        RooMultiCatIter* iter;
        virtual RooArgList containedArgs(Action);
  };

  ClassDef(RooBDecay_, 1) // P.d.f of general description of B decay time distribution
};

#endif


/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   RA,  Roel Aaij,          Nikhef, roel.aaij@nikhef.nl                    *
 *   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                *
 *                                                                           *
 * Copyright (c) 2014, Nikhef. All rights reserved.                          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROO_COMB_EFF_CONSTRAINT
#define ROO_COMB_EFF_CONSTRAINT

#include <memory>
#include <vector>

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooAbsArg;
class RooArgSet;
class RooArgList;

class RooCombEffConstraint : public RooAbsPdf
{
public:

  RooCombEffConstraint();

  RooCombEffConstraint(const char *name, const char *title, Int_t numBins,
      const RooArgList& nu, const RooArgList& eps1A, const RooArgList& eps1B,
      const RooArgList& eps2A, const RooArgList& eps2B,
      const std::vector< std::vector<Double_t> >& sumW,
      const std::vector< std::vector<Double_t> >& sumWSq, Int_t strategy = 0);

  RooCombEffConstraint(const RooCombEffConstraint& other, const char* name=0);

  virtual ~RooCombEffConstraint();

  virtual TObject* clone(const char* newname) const
  {
    return new RooCombEffConstraint(*this, newname);
  }
  
  virtual Bool_t forceAnalyticalInt(const RooAbsArg& /*dep*/) const;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
      const char* rangeName=0) const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

  virtual Double_t getLogVal(const RooArgSet* nset = 0) const;

  virtual Int_t strategy() const {return _strat;}

  virtual const RooArgList& nu()    const {return _nu;   }
  virtual const RooArgList& eps1A() const {return _eps1A;}
  virtual const RooArgList& eps1B() const {return _eps1B;}
  virtual const RooArgList& eps2A() const {return _eps2A;}
  virtual const RooArgList& eps2B() const {return _eps2B;}
  
protected:
   
  virtual Double_t evaluate() const;
   
private:

  inline Double_t nuNorm(const Double_t nu, const Double_t sumW,
      std::vector<Double_t>& normFac, const Int_t ind) const;

  Int_t _strat;
  RooListProxy _nu;
  RooListProxy _eps1A;
  RooListProxy _eps1B;
  RooListProxy _eps2A;
  RooListProxy _eps2B;
  std::vector< std::vector<Double_t> > _sumW;
  std::vector< std::vector<Double_t> > _sumWSq;
  mutable std::vector< std::vector<Double_t> > _normFac;
  Int_t _nBins;

  ClassDef(RooCombEffConstraint, 1) // RooCombEffConstraint
};

#endif

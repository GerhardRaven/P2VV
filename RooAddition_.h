/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooAddition.h,v 1.3 2007/05/11 09:11:30 verkerke Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_ADDITION
#define ROO_ADDITION

#include "RooAbsReal.h"
#include "RooListProxy.h"
#include "RooAICRegistry.h"

class RooRealVar;
class RooArgList ;

class RooAddition_ : public RooAbsReal {
public:

  RooAddition_() ;
  RooAddition_(const char *name, const char *title, const RooArgSet& sumSet, Bool_t takeOwnerShip=kFALSE) ;
  RooAddition_(const char *name, const char *title, const RooArgList& sumSet1, const RooArgList& sumSet2, Bool_t takeOwnerShip=kFALSE) ;
  virtual ~RooAddition_() ;

  RooAddition_(const RooAddition_& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooAddition_(*this, newname); }

  virtual Double_t defaultErrorLevel() const ;

  void printMetaArgs(ostream& os) const ;

  const RooArgList& list() const { return _set ; }

  virtual Bool_t forceAnalyticalInt(const RooAbsArg& /*dep*/) const {
      // Force RooRealIntegral to offer all observables for internal integration
      return kTRUE ;
  }
  Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& numVars, const RooArgSet* normSet, const char* rangeName=0) const;
  Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName=0) const ;
     

protected:

  RooArgList   _ownedList ;      // List of owned components
  RooListProxy _set ;            // set of terms to be summed
  mutable TIterator* _setIter ;  //! Iterator over set

  mutable RooAICRegistry _codeReg;

  Double_t evaluate() const;

  ClassDef(RooAddition_,1) // Sum of RooAbsReal objects
};

#endif

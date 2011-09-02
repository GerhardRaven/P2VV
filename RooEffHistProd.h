/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooEffProd.h,v 1.2 2007/05/11 10:14:56 verkerke Exp $
 * Authors:                                                                  *
 *   GR, Gerhard Raven, NIKHEF/VU                                            *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_EFFHIST_PROD
#define ROO_EFFHIST_PROD

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooSetProxy.h"
#include "RooObjCacheManager.h"

class RooEffHistProd: public RooAbsPdf {
public:
  // Constructors, assignment etc
  inline RooEffHistProd() : _nset(0), _fixedNset(0) { };
  virtual ~RooEffHistProd();
  RooEffHistProd(const char *name, const char *title, RooAbsPdf& pdf, RooAbsReal& efficiency);
  RooEffHistProd(const RooEffHistProd& other, const char* name=0);

  virtual TObject* clone(const char* newname) const { return new RooEffHistProd(*this,newname); }

  virtual RooAbsGenContext* genContext(const RooArgSet &vars, const RooDataSet *prototype,
                                       const RooArgSet* auxProto, Bool_t verbose) const;

  virtual Double_t getVal(const RooArgSet* set=0) const ;

  virtual Bool_t forceAnalyticalInt(const RooAbsArg& /*dep*/) const { 
    // Return kTRUE to force RooRealIntegral to offer all observables for internal integration
    return kTRUE ; 
  }
  Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& numVars, const RooArgSet* normSet, const char* rangeName=0) const ;
  Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName=0) const ;
  
protected:
  
  const RooAbsPdf* pdf() const { 
    // Return pointer to pdf in product
    return (RooAbsPdf*) _pdf.absArg() ; 
  }
  const RooAbsReal* eff() const { 
    // Return pointer to efficiency function in product
    return (RooAbsReal*) _eff.absArg() ; 
  }
  RooRealVar& x() const { return *static_cast<RooRealVar*>(  _observables.first() ) ; }
  
  // Function evaluation
  virtual Double_t evaluate() const ;

   // the real stuff...
  RooRealProxy _pdf ;     // Probability Density function
  RooRealProxy _eff;      // Efficiency function
  RooSetProxy  _observables ; // Observables in the efficiency histogram

  mutable std::map<Int_t, RooArgList> _integralmap ;
  typedef std::vector<Double_t> BinBoundaries ;
  BinBoundaries _binboundaries ;
  
  mutable const RooArgSet* _nset  ; //! Normalization set to be used in evaluation

  RooArgSet* _fixedNset ; //! Fixed normalization set overriding default normalization set (if provided)

  ClassDef(RooEffHistProd,2) // Product operator p.d.f of (PDF x efficiency) implementing optimized generator context
};

#endif

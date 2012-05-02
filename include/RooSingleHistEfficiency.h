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
#ifndef ROO_SINGLEHIST_EFFICIENCY
#define ROO_SINGLEHIST_EFFICIENCY

#include "RooEffHistProd.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooSetProxy.h"
#include "RooObjCacheManager.h"

class RooSingleHistEfficiency : public RooEffHistProd {
public:

   // Constructors, assignment etc
   inline RooSingleHistEfficiency(){};
   virtual ~RooSingleHistEfficiency();
   RooSingleHistEfficiency(const char *name, const char *title, RooAbsPdf& pdf, RooAbsReal& efficiency);
   RooSingleHistEfficiency(const RooSingleHistEfficiency& other, const char* name = 0);

   virtual TObject* clone(const char* newname) const { return new RooSingleHistEfficiency(*this,newname); }

   virtual RooAbsGenContext* genContext(const RooArgSet &vars, const RooDataSet *prototype,
                                        const RooArgSet* auxProto, Bool_t verbose) const;

protected:
   
   virtual const RooArgSet* effObservables() const;
   virtual double effVal() const;

private:
  
   const RooAbsReal* eff() const { 
      // Return pointer to efficiency function in product
      return static_cast<RooAbsReal*>(_eff.absArg());
   }

   // the real stuff...
   RooRealProxy _eff;      // Efficiency function
   RooSetProxy  _observables ; // Observables in the efficiency histogram

   ClassDef(RooSingleHistEfficiency, 1) // Product operator p.d.f of (PDF x efficiency) implementing optimized generator context
};

#endif

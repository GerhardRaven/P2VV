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
   inline RooEffHistProd()   { };
   virtual ~RooEffHistProd();
   RooEffHistProd(const char *name, const char *title, RooAbsPdf& pdf, RooAbsReal& efficiency);
   RooEffHistProd(const RooEffHistProd& other, const char* name=0);

   virtual TObject* clone(const char* newname) const { return new RooEffHistProd(*this,newname); }

   virtual RooAbsGenContext* genContext(const RooArgSet &vars, const RooDataSet *prototype,
                                        const RooArgSet* auxProto, Bool_t verbose) const;


   virtual Bool_t forceAnalyticalInt(const RooAbsArg& /*dep*/) const { 
      // Return kTRUE to force RooRealIntegral to offer all observables for internal integration
      return kTRUE ; 
   }
   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
   Double_t analyticalIntegral(Int_t code,const char* rangeName=0) const ;

   virtual void selectNormalization(const RooArgSet*,Bool_t);
   virtual ExtendMode extendMode() const;

   virtual Double_t expectedEvents(const RooArgSet* nset) const;

private:
  
   const char* makeFPName(const TString& prefix, const RooArgSet& iset, const RooArgSet *nset,
                          const TString& postfix) const;

   const RooAbsPdf* pdf() const { 
      // Return pointer to pdf in product
      return (RooAbsPdf*) _pdf.absArg() ; 
   }
   const RooAbsReal* eff() const { 
      // Return pointer to efficiency function in product
      return (RooAbsReal*) _eff.absArg() ; 
   }
   RooRealVar& x() const { return *static_cast<RooRealVar*>(  _observables.first() ) ; }

   const RooArgSet& observables() const { return static_cast<const RooArgSet&>(_observables); }
  
   // Function evaluation
   virtual Double_t evaluate() const ;

   // the real stuff...
   RooRealProxy _pdf ;     // Probability Density function
   RooRealProxy _eff;      // Efficiency function
   RooSetProxy  _observables ; // Observables in the efficiency histogram

   typedef std::vector<double> BinBoundaries ;

   class CacheElem : public RooAbsCacheElement {
   public:
      CacheElem(const RooEffHistProd* parent,const RooArgSet& iset, const RooArgSet* nset,
                const char *rangeName);
      CacheElem() {}
      virtual ~CacheElem();

      virtual RooArgList containedArgs(Action) ;
      Double_t getVal(Int_t i = 0) { 
         return _I[i]->getVal();
      }
      
      bool trivial() const { return _trivial; }

   private:
      // Payload
      std::vector<RooAbsReal*> _I;
      bool _trivial;
   };
   
   friend class CacheElem;
   CacheElem *getCache(const RooArgSet* nset, const RooArgSet* iset, 
                       const char* rangeName = 0) const;

   mutable RooObjCacheManager _cacheMgr;
   BinBoundaries _binboundaries;
  
   ClassDef(RooEffHistProd, 4) // Product operator p.d.f of (PDF x efficiency) implementing optimized generator context
};

#endif

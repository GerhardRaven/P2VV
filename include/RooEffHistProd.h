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
   RooEffHistProd(const char *name, const char *title, RooAbsPdf& pdf);
   RooEffHistProd(const RooEffHistProd& other, const char* name = 0);

   virtual Bool_t forceAnalyticalInt(const RooAbsArg& /*dep*/) const { 
      // Return kTRUE to force RooRealIntegral to offer all observables for internal integration
      return kTRUE ; 
   }

   Int_t getAnalyticalIntegralWN(RooArgSet& allDeps, RooArgSet& analDeps, 
                                 const RooArgSet* normSet, const char* rangeName) const;
   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code,const char* rangeName=0) const ;

   virtual void selectNormalization(const RooArgSet*,Bool_t);
   virtual ExtendMode extendMode() const;

   virtual Double_t expectedEvents(const RooArgSet* nset) const;

   // Function evaluation
   virtual Double_t getValV(const RooArgSet* set = 0) const;

protected:

   virtual Double_t evaluate() const;

   const RooAbsPdf* pdf() const { 
      // Return pointer to pdf in product
      return (RooAbsPdf*) _pdf.absArg() ; 
   }

   virtual const RooArgSet* effObservables() const = 0;
   virtual double effVal() const = 0;

   typedef std::vector<double> BinBoundaries;
   BinBoundaries _binboundaries;

private:

   const char* makeFPName(const TString& prefix, const RooArgSet& iset, const RooArgSet *nset,
                          const TString& postfix) const;

   // the real stuff...
   RooRealProxy _pdf ;     // Probability Density function
   mutable const RooArgSet* _pdfNormSet;
   mutable const RooArgSet* _fixedNormSet;

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

      void setClone(RooEffHistProd* cl) { _clone = cl; }

      RooEffHistProd* clone() { return _clone; }

      RooArgSet& intObs() { return _intObs; }

   private:
      friend class RooEffHistProd;
      // Payload
      RooArgSet _intObs ;
      RooEffHistProd* _clone;
      std::vector<RooAbsReal*> _I;
      bool _trivial;
   };
   
   friend class CacheElem;
   CacheElem *getCache(const RooArgSet* nset, const RooArgSet* iset, 
                       const char* rangeName = 0, const bool makeClone = false) const;

   mutable RooObjCacheManager _cacheMgr;
  
   ClassDef(RooEffHistProd, 1) // Product operator p.d.f of (PDF x efficiency) implementing optimized generator context
};

#endif

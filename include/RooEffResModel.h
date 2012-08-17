/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooGaussModel.h,v 1.21 2007/05/11 09:13:07 verkerke Exp $
 * Authors:                                                                  *
 *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_EFF_RES_MODEL
#define ROO_EFF_RES_MODEL

#include "RooResolutionModel.h"
#include "RooRealProxy.h"
#include "RooObjCacheManager.h"

class RooCustomizer;
class RooEffResModel : public RooResolutionModel {
public:

  // Constructors, assignment etc
  inline RooEffResModel()  { }
  RooEffResModel(const char *name, const char *title, RooResolutionModel& , RooAbsReal& );
  RooEffResModel(const RooEffResModel& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooEffResModel(*this,newname) ; }
  virtual ~RooEffResModel();
  
  virtual Int_t basisCode(const char* name) const ;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName) const ;
  virtual Bool_t forceAnalyticalInt(const RooAbsArg& /*dep*/) const;

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  //void generateEvent(Int_t code);

protected:
   class CacheElem : public RooAbsCacheElement {
   public:
      CacheElem(const RooEffResModel& parent,const RooArgSet& iset, const TNamed *rangeName);
      virtual ~CacheElem();

      virtual RooArgList containedArgs(Action) ;
      Double_t getVal() { return _I->getVal(); }

      const RooAbsReal* integral() const { return _I; }

   private:
      friend class RooEffResModel;
      // Payload
      RooAbsReal* _I;
      std::vector<RooCustomizer*> _customizers;
   };

  virtual Double_t evaluate() const ;

  // Pointers to our underlying components
  RooRealProxy _model;     // RooResolutionModel
  RooRealProxy _eff;       // RooAbsReal 

  RooResolutionModel& model() const { return (RooResolutionModel&)_model.arg(); }
  RooAbsReal& eff() const { return (RooAbsReal&)*_eff.absArg(); }

   CacheElem *getCache(const RooArgSet* iset, const TNamed* rangeName = 0 ) const;
   mutable RooObjCacheManager _cacheMgr;

  ClassDef(RooEffResModel,1) // EffResian Resolution Model
};

#endif

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
#ifndef ROO_ABS_EFF_RES_MODEL
#define ROO_ABS_EFF_RES_MODEL

#include "RooResolutionModel.h"
#include "RooRealProxy.h"
#include "RooObjCacheManager.h"
#include "RooSetProxy.h"

class RooAbsEffResModel {
public:

   // Constructors, assignment etc
   inline RooAbsEffResModel() { }

   virtual ~RooAbsEffResModel() {}
  
   /** 
    * Get a RooArgSet of all observables
    * @return RooArgSet of observables
    */
   virtual RooArgSet observables() const = 0;

   virtual const RooAbsReal* efficiency() const = 0;

   //TODO: can we do without this one??
   virtual std::vector<const RooAbsReal*> efficiencies() const = 0;

};

#endif

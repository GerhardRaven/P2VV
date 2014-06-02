/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * Authors:                                                                  *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_TPDECAY
#define ROO_TPDECAY

#include "RooAbsAnaConvPdf.h"
#include "RooDecay.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooConstVar.h"
#include "RooRealVar.h"

class RooTPDecay : public RooDecay {
public:

   // Constructors, assignment etc
   inline RooTPDecay() { }
   RooTPDecay(const char *name, const char *title, RooRealVar& t, 
              RooAbsReal& tau, const RooResolutionModel& model,
              DecayType type, const char* normRange);
   RooTPDecay(const RooTPDecay& other, const char* name=0);
   virtual TObject* clone(const char* newname) const { return new RooTPDecay(*this, newname); }
   virtual ~RooTPDecay();

   virtual RooAbsGenContext* genContext(const RooArgSet &vars, const RooDataSet *prototype=0, 
                                        const RooArgSet* auxProto=0, Bool_t verbose= kFALSE) const ;
  
protected:

   RooListProxy _tps;

   ClassDef(RooTPDecay,1) // General decay function p.d.f 
};

#endif

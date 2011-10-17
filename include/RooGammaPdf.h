/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   GR, Gerhard Raven,   Nikhef & VU, Gerhard.Raven@nikhef.nl
 *                                                                           *
 * Copyright (c) 2010, Nikhef & VU. All rights reserved.
 *           
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_GAMMAPDF
#define ROO_GAMMAPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooBinning.h"

class RooRealVar;
class RooArgList ;

class RooGammaPdf : public RooAbsPdf {
public:

  RooGammaPdf(const char *name, const char *title, RooAbsReal& x, RooAbsReal& gamma, RooAbsReal& beta, RooAbsReal& mu);

  RooGammaPdf(const RooGammaPdf& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooGammaPdf(*this, newname); }
  virtual ~RooGammaPdf() ;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:


  RooRealProxy _x;
  RooRealProxy _gamma;
  RooRealProxy _beta;
  RooRealProxy _mu;

  Double_t evaluate() const;

  ClassDef(RooGammaPdf,1) // GammaPdf Pdf
};

#endif

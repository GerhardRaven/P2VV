/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
 * Authors:                                                                  *
 *   GR, Gerhard Raven,   Nikhef & VU, Gerhard.Raven@nikhef.nl
 *                                                                           *
 * Copyright (c) 2010, Nikhef & VU. All rights reserved.
 *           
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////

#include "RooFit.h"

#include "Riostream.h"

#include "RooGammaPdf.h"
#include "RooAbsReal.h"
#include "TMath.h"

ClassImp(RooGammaPdf)
;


//_____________________________________________________________________________
RooGammaPdf::RooGammaPdf(const char* name, const char* title, RooAbsReal& x
                        , RooAbsReal& gamma, RooAbsReal& beta, RooAbsReal& mu)
  : RooAbsPdf(name, title)
  , _x("x", "Dependent", this, x)
  , _gamma("gamma", "gamma", this, gamma)
  , _beta("beta", "beta", this, beta)
  , _mu("mu", "mu", this, mu)
{
}



//_____________________________________________________________________________
RooGammaPdf::RooGammaPdf(const RooGammaPdf& other, const char* name) 
 : RooAbsPdf(other, name)
 , _x("x", this, other._x)
 , _gamma("gamma", this, other._gamma)
 , _beta("beta", this, other._beta)
 , _mu("mu", this, other._mu)
{
  // Copy constructor
}



//_____________________________________________________________________________
RooGammaPdf::~RooGammaPdf()
{
  // Destructor
}

//_____________________________________________________________________________
Int_t RooGammaPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
  if (matchArgs(allVars, analVars, _x)) return 1;
  return 0;
}

//_____________________________________________________________________________
Double_t RooGammaPdf::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  assert(code==1) ;
  Double_t zmax = (_x.max(rangeName)-_mu)/_beta ;
  // 1/Gamma(a) Int_0^x t^(a-1) exp(-t) dt
  double I = TMath::Gamma(_gamma,zmax);
  Double_t zmin = (_x.min(rangeName)-_mu)/_beta ;
  if (zmin!=0) I -= TMath::Gamma(_gamma,zmin);
  return I;
}

//_____________________________________________________________________________
Double_t RooGammaPdf::evaluate() const 
{
    // z = (x-mu)/beta
    // z^(gamma-1) exp(-z) / ( beta Gamma(gamma) )
    return TMath::GammaDist(_x,_gamma,_mu,_beta);
}

/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
 * Authors:                                                                  *
 *   GR, Gerhard Raven,  Vrije Universiteit Amsterdam & Nikhef               *
 *                                                                           *
 * Copyright (c) 2010, Vrije Universiteit Amsterdam & Nikhef                 *
 *                     All rights reserved.                                  *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/


//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Extend RooBDecay to make toy generation to allow the usage of FOAM
// in the presence of RooCategories
// END_HTML
//


#include "RooFit.h"

#include "RooBDecay_.h"
#include "RooRandom.h"
#include "RooAbsCategory.h"
#include "TIterator.h"

ClassImp(RooBDecay_);


//_____________________________________________________________________________
RooBDecay_::RooBDecay_(const char *name, const char* title, 
	       RooRealVar& t, RooAbsReal& tau, RooAbsReal& dgamma,
	       RooAbsReal& f0, RooAbsReal& f1, RooAbsReal& f2, RooAbsReal& f3, 
	       RooAbsReal& dm, const RooResolutionModel& model, DecayType type) :
  RooBDecay(name,title,t,tau,dgamma,f0,f1,f2,f3,dm,model,type)
{
}

//_____________________________________________________________________________
RooBDecay_::RooBDecay_(const RooBDecay_& other, const char* name) :
  RooBDecay(other, name)
{
  //Copy constructor
}

//_____________________________________________________________________________
RooBDecay_::~RooBDecay_()
{
  //Destructor
}

//_____________________________________________________________________________
Int_t RooBDecay_::getMaxVal(const RooArgSet& vars) const {
  cout << "RooBDecay_("<<GetName()<<")::getMaxVal got request for max over " << vars << endl;
  return 0;
}

//_____________________________________________________________________________
Int_t RooBDecay_::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const
{
  cout << "RooBDecay_("<<GetName()<<")::getGenerator got request to generate " << directVars << endl;
  Int_t bcode = RooBDecay::getGenerator(directVars,generateVars,staticInitOK);
  cout << "RooBDecay_("<<GetName()<<")::getGenerator parent returns " << bcode << " and offers " << generateVars << endl;
  return bcode;
}

//_____________________________________________________________________________
void RooBDecay_::initGenerator(Int_t code)
{
    return RooBDecay::initGenerator(code);
}


//_____________________________________________________________________________
void RooBDecay_::generateEvent(Int_t code)
{
  return RooBDecay::generateEvent(code);
}

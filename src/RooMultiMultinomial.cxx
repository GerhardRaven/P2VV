/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id$
 * Authors:                                                                  *
 *   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                *
 *                                                                           *
 * Copyright (c) 2011, Nikhef. All rights reserved.                          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/


//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// END_HTML
//


#include "RooFit.h"
#include "RooMsgService.h"
#include "RooMultiMultinomial.h"

ClassImp(RooMultiMultinomial);

//_____________________________________________________________________________
RooMultiMultinomial::RooMultiMultinomial(const char *name, const char* title) :
  RooAbsReal(name, title)
{
}

//_____________________________________________________________________________
RooMultiMultinomial::RooMultiMultinomial(const RooMultiMultinomial& other,
    const char* name) :
  RooAbsReal(other, name)
{
  // copy constructor
}

//_____________________________________________________________________________
RooMultiMultinomial::~RooMultiMultinomial()
{
  // destructor
}

//_____________________________________________________________________________
Double_t RooMultiMultinomial::evaluate() const
{
  return 0.;
}

//_____________________________________________________________________________
Int_t RooMultiMultinomial::getAnalyticalIntegral(RooArgSet& /*integSet*/,
    RooArgSet& /*anaIntSet*/, const char* /*rangeName*/) const
{
  return 0;
}

//_____________________________________________________________________________
Double_t RooMultiMultinomial::analyticalIntegral(Int_t code,
const char* /*rangeName*/) const
{
  coutF(Eval) << "RooMultiMultinomial::analyticalIntegral(" << GetName()
      << ") code " << code << " not implemented" << endl;

  return 0.;
}


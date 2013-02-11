/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                *
 *                                                                           *
 * Copyright (c) 2012, Nikhef. All rights reserved.                          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/


//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// PDF with a value 1/(2*width) in a window 2*width around its mean and a value
// equal to zero everywhere else.
// END_HTML
//

#include "TMath.h"
#include "P2VV/RooBoxPdf.h"
#include "RooRealVar.h"
#include "RooArgSet.h"

//_____________________________________________________________________________
RooBoxPdf::RooBoxPdf(const char *name, const char* title, RooRealVar& var,
    RooAbsReal& mean, RooAbsReal& width) :
  RooAbsPdf(name, title),
  _var("var", "variable", this, var),
  _mean("mean", "mean of var", this, mean),
  _width("width", "width in var", this, width)
{}

RooBoxPdf::RooBoxPdf(const RooBoxPdf& other, const char* name) :
  RooAbsPdf(other, name),
  _var("var", this, other._var),
  _mean("mean", this, other._mean),
  _width("width", this, other._width)
{}

Double_t RooBoxPdf::evaluate() const
{
  return _var > _mean + _width || _var < _mean - _width ? 0. : 1.;
}

Int_t RooBoxPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
    const char* rangeName) const
{
  // only integrate analytically if the integral is over the full range of _var
  if (rangeName != 0 && _var.absArg()->hasRange(rangeName)) return 0;
  if (matchArgs(allVars, analVars, _var)) return 1;

  return 0;
}

Double_t RooBoxPdf::analyticalIntegral(Int_t /*code*/,
    const char* /*rangeName*/) const
{
  return 2. * _width;
}

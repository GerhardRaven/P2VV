/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id: RooEfficiency.cxx 25184 2008-08-20 13:59:55Z wouter $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
// 
// BEGIN_HTML
// RooEfficiency is a PDF helper class to fit efficiencies parameterized
// by a supplied function F.
// 
// Given a dataset with a category C that determines if a given
// event is accepted or rejected for the efficiency to be measured,
// this class evaluates as F if C is 'accept' and as (1-F) if
// C is 'reject'. Values of F below 0 and above 1 are clipped.
// F may have an arbitrary number of dependents and parameters
// END_HTML
//

#include "RooFit.h"

#include "RooEfficiency.h"
#include "RooEfficiency.h"
#include "RooStreamParser.h"
#include "RooArgList.h"

ClassImp(RooEfficiency)
  ;


//_____________________________________________________________________________
RooEfficiency::RooEfficiency(const char *name, const char *title, const RooAbsReal& effFunc, const RooAbsCategory& cat, const char* sigCatName) :
  RooAbsPdf(name,title),
  _cat("cat","Signal/Background category",this,(RooAbsCategory&)cat),
  _effFunc("effFunc","Efficiency modeling function",this,(RooAbsReal&)effFunc),
  _sigCatName(sigCatName)
{  
  // Construct an N+1 dimensional efficiency p.d.f from an N-dimensional efficiency
  // function and a category cat with M states (0..M=-1) that indicate if a given
  // event should be counted as in the i-th category
}



//_____________________________________________________________________________
RooEfficiency::RooEfficiency(const RooEfficiency& other, const char* name) : 
  RooAbsPdf(other, name),
  _cat("cat",this,other._cat),
  _effFunc("effFunc",this,other._effFunc),
  _sigCatName(other._sigCatName)
{
  // Copy constructor
}



//_____________________________________________________________________________
RooEfficiency::~RooEfficiency() 
{
  // Destructor
}



//_____________________________________________________________________________
Double_t RooEfficiency::evaluate() const
{
  // Calculate the raw value of this p.d.f which is the effFunc
  // value if cat==1 and it is (1-effFunc) if cat==0

  Double_t effFuncVal = _effFunc ;

  // Truncate efficiency function in range 0.0-1.0
  if (_effFunc>1) {
    effFuncVal = 1.0 ;
  } else if (_effFunc<0) {
    effFuncVal = 0.0 ;
  }

  if (_cat.label() == _sigCatName) {
    // Accept case
    return effFuncVal ;
  } else {
    // Reject case
    return 1 - effFuncVal ;
  }
}



//_____________________________________________________________________________
Int_t RooEfficiency::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const 
{
  if (matchArgs(allVars,analVars,_cat)) return 1 ;
  return 0 ;
}



//_____________________________________________________________________________
Double_t RooEfficiency::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  assert(code==1) ;
  return 1.0 ;
}







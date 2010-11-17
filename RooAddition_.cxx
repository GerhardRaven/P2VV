/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id: RooAddition.cxx 28259 2009-04-16 16:21:16Z wouter $
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
// RooAddition calculates the sum of a set of RooAbsReal terms, or
// when constructed with two sets, it sums the product of the terms
// in the two sets. This class does not (yet) do any smart handling of integrals, 
// i.e. all integrals of the product are handled numerically
// END_HTML
//


#include "RooFit.h"

#include "Riostream.h"
#include <math.h>
#include <memory>

#include "RooAddition_.h"
#include "RooProduct.h"
#include "RooAbsReal.h"
#include "RooErrorHandler.h"
#include "RooArgSet.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMsgService.h"

ClassImp(RooAddition_)
;


//_____________________________________________________________________________
RooAddition_::RooAddition_()
    : _codeReg(10)
{
  _setIter = _set.createIterator() ;
}



//_____________________________________________________________________________
RooAddition_::RooAddition_(const char* name, const char* title, const RooArgSet& sumSet, Bool_t takeOwnership) 
  : RooAbsReal(name, title)
  , _set("!set","set of components",this)
  , _codeReg(10)
{
  // Constructor with a single set of RooAbsReals. The value of the function will be
  // the sum of the values in sumSet. If takeOwnership is true the RooAddition object
  // will take ownership of the arguments in sumSet

  _setIter = _set.createIterator() ;

  std::auto_ptr<TIterator> inputIter( sumSet.createIterator() );
  RooAbsArg* comp ;
  while((comp = (RooAbsArg*)inputIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(comp)) {
      coutE(InputArguments) << "RooAddition::ctor(" << GetName() << ") ERROR: component " << comp->GetName() 
			    << " is not of type RooAbsReal" << endl ;
      RooErrorHandler::softAbort() ;
    }
    _set.add(*comp) ;
    if (takeOwnership) _ownedList.addOwned(*comp) ;
  }

}



//_____________________________________________________________________________
RooAddition_::RooAddition_(const char* name, const char* title, const RooArgList& sumSet1, const RooArgList& sumSet2, Bool_t takeOwnership) 
    : RooAbsReal(name, title)
    , _set("!set","First set of components",this)
    , _codeReg(10)
{
  // Constructor with two set of RooAbsReals. The value of the function will be
  //
  //  A = sum_i sumSet1(i)*sumSet2(i) 
  //
  // If takeOwnership is true the RooAddition object will take ownership of the arguments in sumSet

  _setIter = _set.createIterator() ;

  if (sumSet1.getSize() != sumSet2.getSize()) {
    coutE(InputArguments) << "RooAddition::ctor(" << GetName() << ") ERROR: input lists should be of equal length" << endl ;
    RooErrorHandler::softAbort() ;    
  }

  std::auto_ptr<TIterator> inputIter1( sumSet1.createIterator() );
  std::auto_ptr<TIterator> inputIter2( sumSet2.createIterator() );
  RooAbsArg *comp1(0),*comp2(0) ;
  while((comp1 = (RooAbsArg*)inputIter1->Next())) {
    if (!dynamic_cast<RooAbsReal*>(comp1)) {
      coutE(InputArguments) << "RooAddition::ctor(" << GetName() << ") ERROR: component " << comp1->GetName() 
			    << " in first list is not of type RooAbsReal" << endl ;
      RooErrorHandler::softAbort() ;
    }
    comp2 = (RooAbsArg*)inputIter2->Next();
    if (!dynamic_cast<RooAbsReal*>(comp2)) {
      coutE(InputArguments) << "RooAddition::ctor(" << GetName() << ") ERROR: component " << comp2->GetName() 
			    << " in first list is not of type RooAbsReal" << endl ;
      RooErrorHandler::softAbort() ;
    }
    // TODO: add flag to RooProduct c'tor to make it assume ownership...
    TString name(comp1->GetName());
    name.Append( "_x_");
    name.Append(comp2->GetName());
    RooProduct  *prod = new RooProduct( name, name , RooArgSet(*comp1, *comp2) /*, takeOwnership */ ) ;
    _set.add(*prod);
    _ownedList.addOwned(*prod) ;
    if (takeOwnership) {
        _ownedList.addOwned(*comp1) ;
        _ownedList.addOwned(*comp2) ;
    }
  }
}



//_____________________________________________________________________________
RooAddition_::RooAddition_(const RooAddition_& other, const char* name) 
    : RooAbsReal(other, name)
    , _set("!set",this,other._set)
    , _codeReg(other._codeReg)
{
  // Copy constructor
  _setIter = _set.createIterator() ;
  
  // Member _ownedList is intentionally not copy-constructed -- ownership is not transferred
}


//_____________________________________________________________________________
RooAddition_::~RooAddition_() 
{
  // Destructor
  delete _setIter ;
}

//_____________________________________________________________________________
Double_t RooAddition_::evaluate() const 
{
  // Calculate and return current value of self
  Double_t sum(0);
  const RooArgSet* nset = _set.nset() ;

  _setIter->Reset() ;
  RooAbsReal* comp ;
  while((comp=(RooAbsReal*)_setIter->Next())) {
      sum += comp->getVal(nset) ;
  }
  return sum ;
}

//_____________________________________________________________________________
Double_t RooAddition_::defaultErrorLevel() const 
{
  // Return the default error level for MINUIT error analysis
  // If the addition contains one or more RooNLLVars and 
  // no RooChi2Vars, return the defaultErrorLevel() of
  // RooNLLVar. If the addition contains one ore more RooChi2Vars
  // and no RooNLLVars, return the defaultErrorLevel() of
  // RooChi2Var. If the addition contains neither or both
  // issue a warning message and return a value of 1

  RooAbsReal* nllArg(0) ;
  RooAbsReal* chi2Arg(0) ;

  RooAbsArg* arg ;

  _setIter->Reset() ;
  while((arg=(RooAbsArg*)_setIter->Next())) {
    if (dynamic_cast<RooNLLVar*>(arg)) {
      nllArg = (RooAbsReal*)arg ;
    }
    if (dynamic_cast<RooChi2Var*>(arg)) {
      chi2Arg = (RooAbsReal*)arg ;
    }
  }

  if (nllArg && !chi2Arg) {
    coutI(Fitting) << "RooAddition::defaultErrorLevel(" << GetName() 
		   << ") Summation contains a RooNLLVar, using its error level" << endl ;
    return nllArg->defaultErrorLevel() ;
  } else if (chi2Arg && !nllArg) {
    coutI(Fitting) << "RooAddition::defaultErrorLevel(" << GetName() 
		   << ") Summation contains a RooChi2Var, using its error level" << endl ;
    return chi2Arg->defaultErrorLevel() ;
  } else if (!nllArg && !chi2Arg) {
    coutI(Fitting) << "RooAddition::defaultErrorLevel(" << GetName() << ") WARNING: "
		   << "Summation contains neither RooNLLVar nor RooChi2Var server, using default level of 1.0" << endl ;
  } else {
    coutI(Fitting) << "RooAddition::defaultErrorLevel(" << GetName() << ") WARNING: "
		   << "Summation contains BOTH RooNLLVar and RooChi2Var server, using default level of 1.0" << endl ;
  }

  return 1.0 ;
}


//_____________________________________________________________________________
void RooAddition_::printMetaArgs(ostream& os) const 
{
  _setIter->Reset() ;

  Bool_t first(kTRUE) ;
  RooAbsArg* arg;
  while((arg=(RooAbsArg*)_setIter->Next())) {
    if (!first) { os << " + " ;
    } else { first = kFALSE ; 
    }
    os << arg->GetName() ; 
  }  
  os << " " ;    
}

Int_t RooAddition_::getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars, const RooArgSet* normSet, const char* rangeName) const
{
    std::auto_ptr<RooArgSet> allDepVars( getObservables(allVars) );
    RooArgSet allAnalVars(*allDepVars) ;

    // check which variables can be analytically integrated over...
    _setIter->Reset();
    std::auto_ptr<TIterator> avIter( allVars.createIterator() );
    RooAbsReal* arg;
    while((arg=(RooAbsReal*)_setIter->Next())) { // checked in c'tor to be RooAbsReal
        RooArgSet subAnalVars;
        arg->getAnalyticalIntegralWN(allVars,subAnalVars,normSet,rangeName);
        avIter->Reset();
        RooAbsArg *var;
        while ( (var=(RooAbsArg*)avIter->Next())) {
            if (!subAnalVars.find(var->GetName())&& arg->dependsOn(*arg) ) allAnalVars.remove(*arg,kTRUE,kTRUE);
        }
    }
    if (allAnalVars.getSize()==0) return 0;

// Now retrieve codes for integration over common set of analytically integrable observables for each component
  _setIter->Reset() ;
  Int_t n=0 ;
  Int_t* subCode = new Int_t[_set.getSize()];
  Bool_t allOK(kTRUE) ;
  while((arg=(RooAbsReal*)_setIter->Next())) {
    RooArgSet subAnalVars ;
    std::auto_ptr<RooArgSet> allAnalVars2( arg->getObservables(allAnalVars) );
    subCode[n] = arg->getAnalyticalIntegralWN(*allAnalVars2,subAnalVars,normSet,rangeName) ;
    if (subCode[n]==0 && allAnalVars2->getSize()>0) {
      coutE(InputArguments) << "RooAddition::getAnalyticalIntegral(" << GetName() << ") WARNING: component RooAbsReal " << arg->GetName()
                << "   advertises inconsistent set of integrals (e.g. (X,Y) but not X or Y individually."
                << "   Distributed analytical integration disabled. Please fix" << endl ;
      allOK = kFALSE ;
    }
    n++ ;
  }
  if (!allOK) return 0 ;

  // Mare all analytically integrated observables as such
  analVars.add(allAnalVars) ;

  // Store set of variables analytically integrated
  RooArgSet* intSet = new RooArgSet(allAnalVars) ;
  Int_t masterCode = _codeReg.store(subCode,_set.getSize(),intSet)+1 ;
  delete[] subCode;

  return masterCode;
}

Double_t 
RooAddition_::analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName) const 
{
    if (code==0) return getVal(normSet);
    // Retrieve analytical integration subCodes and set of observabels integrated over
    RooArgSet* intSet ;
    const Int_t* subCode = _codeReg.retrieve(code-1,intSet) ;
    if (!subCode) {
      coutE(InputArguments) << "RooAddPdf::analyticalIntegral(" << GetName() << "): ERROR unrecognized integration code, " << code << endl ;
      assert(0) ;
    }
    _setIter->Reset() ;
    Int_t n(0);
    double result(0);
    RooAbsReal *arg;
    while((arg=(RooAbsReal*)_setIter->Next())) {
        result += arg->analyticalIntegralWN( subCode[n], normSet, rangeName );
        ++n;
    }
    return result;
}

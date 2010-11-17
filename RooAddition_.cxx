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
#include "RooNameReg.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMsgService.h"

ClassImp(RooAddition_)
;


//_____________________________________________________________________________
RooAddition_::RooAddition_()
{
  _setIter = _set.createIterator() ;
}



//_____________________________________________________________________________
RooAddition_::RooAddition_(const char* name, const char* title, const RooArgSet& sumSet, Bool_t takeOwnership) 
  : RooAbsReal(name, title)
  , _set("!set","set of components",this)
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
{
  // Copy constructor
  _setIter = _set.createIterator() ;
  
  // Member _ownedList is intentionally not copy-constructed -- ownership is not transferred
}


//_____________________________________________________________________________
RooAddition_::~RooAddition_() 
{ // Destructor
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

Int_t RooAddition_::getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars, const RooArgSet* nSet, const char* rangeName) const
{
  // we always do things ourselves -- actually, always delegate further down the line ;-)
  analVars.add(allVars);

  // check if we already have integrals for this combination of factors
  Int_t sterileIndex(-1);
  CacheElem* cache = (CacheElem*) _cacheMgr.getObj(&allVars,nSet,&sterileIndex,RooNameReg::ptr(rangeName));
  if (cache!=0) {
    Int_t code = _cacheMgr.lastIndex();
    return code;
  }

  // we don't, so we make it right here....

  cache = new CacheElem;
  _setIter->Reset();
  RooAbsReal *arg(0);
  while( (arg=(RooAbsReal*)_setIter->Next())!=0 ) {  // checked in c'tor that this will work...
      RooAbsReal *I = arg->createIntegral(allVars,*nSet,rangeName);
      cache->_I.addOwned(*I);
  }

  Int_t code = _cacheMgr.setObj(&allVars,nSet,(RooAbsCacheElement*)cache,RooNameReg::ptr(rangeName));
  cout << " RooAddition_("<< GetName()<<")::getAnalyticalIntegral -- returning " << 1+code << " and analVars: " <<  endl;
  analVars.Print("V");
  return 1+code;
}

Double_t 
RooAddition_::analyticalIntegralWN(Int_t code, const RooArgSet* nSet, const char* rangeName) const 
{
  // Calculate integral internally from appropriate integral cache

  // note: rangeName implicit encoded in code: see _cacheMgr.setObj in getPartIntList...
  CacheElem *cache = (CacheElem*) _cacheMgr.getObjByIndex(code-1);
  if (cache==0) {
    // cache got sterilized, trigger repopulation of this slot, then try again...
    std::auto_ptr<RooArgSet> vars( getParameters(RooArgSet()) );
    std::auto_ptr<RooArgSet> iset(  _cacheMgr.nameSet2ByIndex(code-1)->select(*vars) );
    RooArgSet dummy;
    Int_t code2 = getAnalyticalIntegralWN(*iset,dummy,nSet,rangeName);
    assert(code==code2); // must have revived the right (sterilized) slot...
    return analyticalIntegralWN(code2,nSet,rangeName);
  }
  assert(cache!=0);

  // loop over cache, and sum...
  std::auto_ptr<TIterator> iter( cache->_I.createIterator() );
  RooAbsReal *I;
  double result(0);
  while ( ( I=(RooAbsReal*)iter->Next() ) != 0 ) {
        result += I->getVal();
  }
  return result;

}

//_____________________________________________________________________________
RooArgList RooAddition_::CacheElem::containedArgs(Action)
{
  // Return list of all RooAbsArgs in cache element
  RooArgList ret(_I) ;
  return ret ;
}

RooAddition_::CacheElem::~CacheElem()
{
  // Destructor
}



Bool_t RooAddition_::isJacobianOK(const RooArgSet& ) const 
{
    return kTRUE; 
}

Double_t RooAddition_::jacobian() const
{
    return 1.0; 
}

    


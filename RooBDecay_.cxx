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
  RooBDecay_(other, name)
{
  //Copy constructor
}

//_____________________________________________________________________________
RooBDecay_::~RooBDecay_()
{
  //Destructor
}

//_____________________________________________________________________________
Int_t RooBDecay_::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const
{
  if (staticInitOK) {
      // we will generate _all_ categories, plus any real ones supported by our parent...
      Int_t sterileIndex(-1);
      allObs = this->getObservables( directVars );
      catObs = grab list of categories
      otherObs( allObs )
      otherObs.remove( catObs )
      parentCode = RooBDecay::getGenerator( otherObs, generateVars, staticInitOK)
      generateVars.add( catObs );
      // make sure catObs & otherObs stay alive!

      // first time, we will fail ;-)
      CacheElem *cache = (CacheElem*) _cacheMgr.getObj( &catObs, &generateVars, &sterileIndex)
      if (cache!=0) {
         Int_t code = _cacheMgr.lastIndex();
         return code*100+RooBDecay::getGenerator(directVars,generateVars,staticInitOK);
      }
      cache = new CacheElem;

      cache._genObs = &generateVars;
      cache._catObs = &catObs;
      cache->iter = new RooMultCatIter( *cache->catObs );

      code = _cacheMgr.setObj( cache._allObs, cache._catObs, cache );
      return code*100+RooBDecay::getGenerator(directVars,generateVars,staticInitOK);
  } 
  return RooBDecay::getGenerator(directVars,generateVars,staticInitOK);
}

//_____________________________________________________________________________
Int_t RooBDecay_::initGenerator(Int_t code)
{
    CacheElem *cache = (CacheElem*) _cacheMgr.getObjByIndex(code/100);
    if (cache == 0) { // (re)generate it and try again...
        // resurrect the two relevant sets...
        RooArgSet s = _cacheMgr.nameSet1ByIndex(code/100).select( getObservables() );
        s.add(        _cacheMgr.nameSet2ByIndex(code/100).select( getObservables() ) )
        getGenerator( s, dummy, true  );
        return initGenerator(code);
    }
    // we have a cache... but get getGenerator can not actually compute the limits, we
    // have to wait until we've arrived here.....
    if (cache->n==0) {

        int ncomb = 0;
        const TCollection *collection = cache->iter->GetCollection()
        if (collection!=0) {  // is a TIterator always backed up by a collection??
            ncomb = collection.GetSize();
        } else {
            cache->iter.Reset();
            while ( iter.Next() != 0 ) { ++ncomb; }
        }
        cache->n = new double[ncomb];

        std::auto_ptr<RooAbsPdf> marginal( createProjection( cache->otherObs ) );
        RooArgSet nset( otherObs )
        nset.add( catObs );
        std::auto_ptr<RooAbsReal> normInt( createIntegral( nset ) );
        double norm = normInt->getVal()
        int i=0;
        cache->iter->Reset();
        while ( cache->iter->Next() ) {
            cache->n[i] = marginal->getVal()/norm; // fraction of events in this combination 
            if (i>0) cache->n[i]+=cache->n[i-1];
            ++i;
        }
    }
}


//_____________________________________________________________________________
void RooBDecay_::generateEvent(Int_t code)
{
  if (code>=100) {
     CacheElem *cache = (CacheElem*) _cacheMgr.getObjByIndex( code/100 );
     if (cache==0) { // repopulate and try again...
            initGenerator(code);
            return generateEvent(code);
     }
     Double_t r = RooRandom::uniform();
     _cache->iter.Reset();
     int i=0;
     while ( _cache->iter.Next()!=0 ) { if (r < _cache->n[i++]) break; }
     return RooBDecay::generateEvent(code%100);
  }
  return RooBDecay::generateEvent(code);
}

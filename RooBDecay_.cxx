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
#if 0
  if (staticInitOK) {
      // we will generate _all_ categories, on top of whatever our parent is capable off...
      Int_t sterileIndex(-1);
      // we do generateVars + all categories...
      RooArgSet marginalObs( directVars );  // we need to figure out the marginalization set, which is directVars - generateVars 
      marginalObs.remove( generateVars );
      RooArgSet catObs;
      std::auto_ptr<TIterator> iter( marginalObs.createIterator() );
      TObject *obj(0);
      while ( ( obj = iter->Next() ) !=0 ) {
        RooAbsCategory *cat = dynamic_cast<RooAbsCategory*>(obj);
        if (cat!=0) catObs.add(*cat);
      }
      marginalObs.remove( catObs );
      generateVars.add( catObs ); 

      // first time, we will fail ;-)
      CacheElem *cache = (CacheElem*) _cacheMgr.getObj( &generateVars, &generateVars, &sterileIndex);
      if (cache!=0) {
         Int_t code = _cacheMgr.lastIndex();
         return code*100+bcode;
      }
      cache = new CacheElem;

      cache->_genObs = &generateVars;
      cache->_catObs = &catObs;
      cache->iter = new RooMultiCatIter( *cache->catObs );

      code = _cacheMgr.setObj( cache._allObs, cache._catObs, cache );
      return code*100+bcode;
  }  else {
      return bcode;
  }
#endif
}

//_____________________________________________________________________________
void RooBDecay_::initGenerator(Int_t code)
{
    return;
#if 0
    if (code/100==0) return;
    CacheElem *cache = (CacheElem*) _cacheMgr.getObjByIndex(code/100);
    if (cache == 0) { // (re)generate it and try again...
        // resurrect the two relevant sets...
        RooArgSet s = _cacheMgr.nameSet1ByIndex(code/100).select( getObservables() );
        s.add(        _cacheMgr.nameSet2ByIndex(code/100).select( getObservables() ) )
        getGenerator( s, dummy, true  ); // we can put true, because we got a code/100 != 0
        return initGenerator(code);
    }
    // we have a cache... but get getGenerator can not actually compute the limits, we
    // have to wait until we've arrived here.....
    if (cache->n==0) {

        int ncomb = 0;
        const TCollection *collection = cache->iter->GetCollection();
        if (collection!=0) {  // is a TIterator always backed up by a collection??
            ncomb = collection.GetSize();
        } else {
            cache->iter.Reset();
            while ( cache->iter.Next() != 0 ) { ++ncomb; }
        }
        cache->n = new double[ncomb];

        std::auto_ptr<RooAbsPdf> marginal( createProjection( cache->otherObs ) );
        RooArgSet nset( otherObs );
        nset.add( catObs );
        std::auto_ptr<RooAbsReal> normInt( createIntegral( nset ) );
        double norm = normInt->getVal();
        int i=0;
        cache->iter->Reset();
        while ( cache->iter->Next() ) {
            cache->n[i] = marginal->getVal()/norm; // fraction of events in this combination 
            if (i>0) cache->n[i]+=cache->n[i-1];   // cumulative
            ++i;
        }
    }
#endif
}


//_____________________________________________________________________________
void RooBDecay_::generateEvent(Int_t code)
{
#if 0
  if (code/100!=0) {
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
#endif
  return RooBDecay::generateEvent(code);
}

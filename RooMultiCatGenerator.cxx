/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
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
// Class RooMultiCatGenerator is a generic toy monte carlo generator that 
// generates discrete categories, and dispatches any real variables to the
// underlying PDF.
// END_HTML
//


#include "RooFit.h"
#include "Riostream.h"

#include "RooMultiCatGenerator.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooCategory.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooRandom.h"
#include "RooErrorHandler.h"

#include "TString.h"
#include "TIterator.h"
#include "RooMsgService.h"
#include "TClass.h"
#include "RooNumGenFactory.h"
#include "RooNumGenConfig.h"

#include <assert.h>

ClassImp(RooMultiCatGenerator)
  ;


//_____________________________________________________________________________
void RooMultiCatGenerator::registerSampler(RooNumGenFactory& fact)
{
  // Register RooIntegrator1D, is parameters and capabilities with RooNumIntFactory
  fact.storeProtoSampler(new RooMultiCatGenerator,RooArgSet()) ;
}


const char *makeFPName(const RooAbsReal& func, const RooArgSet& terms ) {
    static TString pname;
    pname = func.GetName();
    pname.Append("_");
    std::auto_ptr<TIterator> i( terms.createIterator() );
    RooAbsArg *arg;
    Bool_t first(kTRUE);
    while((arg=(RooAbsArg*)i->Next())) {
         if (first) { first=kFALSE;}
         else pname.Append("_X_");
         pname.Append(arg->GetName());
    }
    return pname.Data();
}

RooSuperCategory makeSuper(const RooAbsReal& func, const RooArgSet& _catVars  ) {
  const char *catName =  makeFPName(func, _catVars );
  return RooSuperCategory(catName,catName,_catVars);
}

//_____________________________________________________________________________
RooMultiCatGenerator::RooMultiCatGenerator(const RooAbsReal &func, const RooArgSet &genVars, const RooNumGenConfig& config, Bool_t verbose, const RooAbsReal* maxFuncVal) 
  : RooAbsNumGenerator(func,genVars,verbose,maxFuncVal)
  , _super( makeSuper( func, _catVars ) )
{
  // could do createIntegral... but for now do brute-force
  RooAbsPdf& pdf = dynamic_cast<RooAbsPdf&>(*_funcClone);
  std::auto_ptr<RooAbsPdf> marginal( pdf.createProjection( _realVars ) );
  std::auto_ptr<TIterator> superIter( _super.MakeIterator() );
  while ( superIter->Next() ) {
            double n = marginal->getVal(); // fraction of events in this combination 
            if (!_realGenerators.empty()) n += _realGenerators.back().first;   // cumulative
            _realGenerators.push_back(make_pair(n, RooNumGenFactory::instance().createSampler(pdf,_realVars,RooArgSet(),*pdf.getGeneratorConfig(), verbose ))); 
  }
  for (Generators::iterator i=_realGenerators.begin();i!=_realGenerators.end();++i) i->first /= _realGenerators.back().first; // normalize
}


//_____________________________________________________________________________
RooMultiCatGenerator::~RooMultiCatGenerator() 
{
  // Destructor
  for (Generators::iterator i=_realGenerators.begin();i!=_realGenerators.end();++i) delete i->second;
}



//_____________________________________________________________________________
const RooArgSet *RooMultiCatGenerator::generateEvent(UInt_t remaining, Double_t& resampleRatio) 
{
  // are we actually generating anything? (the cache always contains at least our function value)
  const RooArgSet *event= _cache->get();
  if(event->getSize() == 1) return event;

   Double_t r = RooRandom::uniform();
   std::auto_ptr<TIterator> superIter( _super.MakeIterator() );
   // find the right generator, and generate categories at the same time...
   Generators::iterator gen = _realGenerators.begin();
   while ( superIter->Next()!=0 && gen->first < r)  ++gen;
   _super.setLabel( dynamic_cast<TObjString&>(*superIter).String() ); // this should assign _catVars...

   cout << "RooMultCatGenerator("<<GetName()<<")::generateEvent: " << _catVars << endl;

   // and generate the real vars for this combination of categories
  gen->second->generateEvent( remaining, resampleRatio);

  // calculate and store our function value at this new point
  Double_t val= _funcClone->getVal();
  _funcValPtr->setVal(val);

   cout << "RooMultCatGenerator("<<GetName()<<")::generateEvent: " << _cloneSet << endl;

  // Transfer contents to dataset
  return _cloneSet;
}

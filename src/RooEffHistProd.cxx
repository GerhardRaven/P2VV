/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id: RooEffHistProd.cxx 25184 2008-08-20 13:59:55Z wouter $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, NIKHEF
 *   GR, Gerhard Raven, NIKHEF/VU                                            *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/


/////////////////////////////////////////////////////////////////////////////////////
// BEGIN_HTML
// The class RooEffHistProd implements the product of a PDF with an efficiency function.
// The normalization integral of the product is calculated numerically, but the
// event generation is handled by a specialized generator context that implements
// the event generation in a more efficient for cases where the PDF has an internal
// generator that is smarter than accept reject. 
// END_HTML
//

#include "RooFit.h"
#include "RooEffHistProd.h"
#include "RooEffGenContext.h"
#include "RooNameReg.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooCategory.h"
#include "RooCatType.h"
#include "TString.h"
#include <RooCustomizer.h>
#include <RooConstVar.h>
#include <RooAddition.h>

#include <memory>
#include <algorithm>
#include <exception>

ClassImp(RooEffHistProd);

namespace EffHistProd {
//_____________________________________________________________________________
void cloneRanges(const RooArgSet& observables, const RooArgSet& iset,
            const RooArgSet* nset, const char* rangeName,
            const char* newName)
{
   assert(rangeName != 0);
   assert(newName != rangeName);
   // Create a set which is the union of iset and nset to find other variables with the same range
   RooArgSet uni(iset);
   uni.add(*nset, kTRUE);

   // Skip the observables for which we are explicitely modifying the range
   uni.remove(observables);

   std::auto_ptr<TIterator> iter(uni.createIterator());
   RooAbsArg* arg = 0;
   while ((arg = static_cast<RooAbsArg*>(iter->Next()))) {
      if (arg->hasRange(rangeName)) {
         if (RooRealVar* real = dynamic_cast<RooRealVar*>(arg)) {
            const RooAbsBinning& range = real->getBinning(rangeName);
            if (range.isParameterized()) {
               // Deal with parametrized ranges
               real->setRange(newName, *range.lowBoundFunc(), *range.highBoundFunc());
            } else {
               // Regular range
               real->setRange(newName, range.lowBound(), range.highBound());
            }
         } else if (RooCategory* cat = dynamic_cast<RooCategory*>(arg)) {
            std::auto_ptr<TIterator> it(cat->typeIterator());
            TString states;
            const RooCatType* ct = 0;
            while ((ct = static_cast<const RooCatType*>(it->Next()))) {
               if (cat->isStateInRange(rangeName, ct->GetName())) {
                  states.Append(ct->GetName());
                  states.Append(",");
               }
            }
            states.Remove(TString::kTrailing, ',');
            cat->setRange(newName, states.Data());
         } else {
            throw EffHistProd::Exception("Got type for which a range cannot be cloned");
         }
      }
   }
}
}

//_____________________________________________________________________________
RooEffHistProd::CacheElem::CacheElem(const RooEffHistProd* parent, const RooArgSet& iset,
                                     const RooArgSet* nset, const char *rangeName)
   : _clone(0), _I(0)
{
   // create a function object for the corresponding integral of the underlying PDF.
   // i.e. if we have  eps(x)f(x,y)  and we get (x,y) as allVars, 
   // construct I_i,j = int_xmin,i^xmax,j dx int dy f(x,y)
   // so that later we can do sum_i eps( (x_i+x_i+1)/2 ) * I(x_i,x_i+1)
   RooArgSet x_;
   const RooArgSet *effObs = parent->eff()->getObservables(iset); // the subset of iset on which _eff depends
   RooFIter iter = iset.fwdIterator();
   while(RooAbsArg* arg = iter.next()) {
      if (effObs->contains(*arg)) {
         x_.add(*arg);
      }
   }
   delete effObs;

   _intObs.addClone(*nset);

   const BinBoundaries& bounds = parent->_binboundaries;

   assert(x_.getSize() < 2); // for now, we only do 1D efficiency histograms...
   _trivial = x_.getSize()==0;

   if (!trivial()) {
      RooRealVar* x = static_cast<RooRealVar*>(x_.first());
      Double_t xmin = x->getMin(rangeName);
      Double_t xmax = x->getMax(rangeName);

      RooArgList effList;
      RooArgList intList;

      for(unsigned int i = 0; i < bounds.size() - 1; ++i) {
         const char* newRange = rangeName;
         Double_t low = bounds[i];
         Double_t high = bounds[i + 1];
         if (high < xmin) continue; // not there yet...
         if (low  > xmax) break;    // past the requested interval...

         Double_t thisxmin = std::max(low,  xmin);
         Double_t thisxmax = std::min(high, xmax);

         TString range = TString::Format("_I_Range_%d", i);
         // Add original rangeName if there is one
         if (0 != newRange) range.Replace(3, 8, newRange);
      
         // Create a new name for the range
         newRange = parent->makeFPName(parent->GetName(), iset, nset, range);
         x->setRange(newRange, thisxmin, thisxmax);

         if (0 != rangeName) {
            EffHistProd::cloneRanges(x_, iset, nset, rangeName, newRange);
         }
         TString suffix = "bin_"; suffix += i;
         TString name = x->GetName(); name += newRange; name += suffix;
         TString binName = x->GetName(); binName += "_"; binName += suffix;
         RooCustomizer customizer(*parent->eff(), name.Data());
         customizer.replaceArg(*x, RooConstVar(binName.Data(), binName.Data(),
                                              0.5 * (thisxmin + thisxmax)));
         RooAbsReal* I = parent->pdf()->createIntegral(iset, nset, parent->getIntegratorConfig(),
                                                       newRange);
         effList.add(*customizer.build(kTRUE));
         intList.add(*I);
      }
      _I = new RooAddition("integral", "integral", effList, intList, kTRUE);
   } else {
      _I = parent->pdf()->createIntegral(iset, nset, parent->getIntegratorConfig(), rangeName);
   }
}

//_____________________________________________________________________________
RooArgList RooEffHistProd::CacheElem::containedArgs(Action) 
{
   // Return list of all RooAbsArgs in cache element
   RooArgList l(_intObs);
   l.add(*_I);
   l.add(*_clone);
   return l;
}

//_____________________________________________________________________________
RooEffHistProd::CacheElem::~CacheElem() 
{
   if (_I) delete _I;
   if (_clone) delete _clone;
}

//_____________________________________________________________________________
RooEffHistProd::RooEffHistProd(const char *name, const char *title, 
                               RooAbsPdf& inPdf, RooAbsReal& eff)
   : RooAbsPdf(name, title),
     _pdf("pdf", "pre-efficiency pdf", this, inPdf),
     _eff("efficiency", "efficiency", this, eff),
     _pdfNormSet(0),
     _fixedNormSet(0),
     _cacheMgr(this, 10)
{  

}

//_____________________________________________________________________________
RooEffHistProd::RooEffHistProd(const RooEffHistProd& other, const char* name):
   RooAbsPdf(other, name),
   _binboundaries(other._binboundaries),
   _pdf("pdf", this, other._pdf),
   _eff("eff", this, other._eff),
   _pdfNormSet(other._pdfNormSet),
   _fixedNormSet(other._fixedNormSet),
   _cacheMgr(other._cacheMgr, this)
{
   // Copy constructor
}

//_____________________________________________________________________________
RooEffHistProd::~RooEffHistProd() 
{
   // Destructor
}

//_____________________________________________________________________________
Double_t RooEffHistProd::getValV(const RooArgSet* normSet) const 
{  
   // Return p.d.f. value normalized over given set of observables
   // cout << "RooEffHistProd::getValV " << (normSet ? *normSet : RooArgSet()) << endl;
   // FIXME: memory leak!!
   _pdfNormSet = _fixedNormSet ? _fixedNormSet : pdf()->getObservables(*normSet);
   return RooAbsPdf::getValV(normSet);
}

//_____________________________________________________________________________
Double_t RooEffHistProd::evaluate() const
{
   // Calculate and return 'raw' unnormalized value of p.d.f
   double pdfVal = pdf()->getVal(_pdfNormSet);
   double effVal = eff()->getVal();
   cout << "RooEffHistProd::evaluate " << effVal << " " << pdfVal << " " << (_pdfNormSet ? *_pdfNormSet : RooArgSet()) << endl;
   return effVal * pdfVal;
}

//_____________________________________________________________________________
const char* RooEffHistProd::makeFPName(const TString& prefix,const RooArgSet& iset, 
                                       const RooArgSet* nset, const TString& postfix) const
{
   static TString pname;
   pname = prefix;
   if (prefix.Sizeof()) pname.Append("_");
   std::auto_ptr<TIterator> i(iset.createIterator());
   RooAbsArg *arg;
   Bool_t first(kTRUE);
   pname.Append("I_");
   while((arg=(RooAbsArg*)i->Next())) {
      if (first) { 
         first = kFALSE;
      } else {
         pname.Append("_X_");
      }
      pname.Append(arg->GetName());
   }
   if (nset) {
      pname.Append("_N_");
      std::auto_ptr<TIterator> it(nset->createIterator());
      first = kTRUE;
      while ((arg = (RooAbsArg*)it->Next())) {
         if (first) {
            first=kFALSE;
         } else {
            pname.Append("_X_");
         }
         pname.Append(arg->GetName());
      }
   }
   pname.Append(postfix);
   return RooNameReg::str(RooNameReg::ptr(pname.Data()));
}

//_____________________________________________________________________________
void RooEffHistProd::selectNormalization(const RooArgSet* nset,Bool_t force) {
   RooAbsPdf::selectNormalization(nset,force);
}

//_____________________________________________________________________________
RooAbsPdf::ExtendMode RooEffHistProd::extendMode() const
{
   // If this product contains exactly one extendable p.d.f return the extension abilities of
   // that p.d.f, otherwise return CanNotBeExtended
   return pdf()->extendMode();
}

//_____________________________________________________________________________
Double_t RooEffHistProd::expectedEvents(const RooArgSet* nset) const
{
   // Return the expected number of events associated with the extentable input p.d.f
   // in the product. If there is no extendable term, return zero and issue and error
   return pdf()->expectedEvents(nset);
}

//_____________________________________________________________________________
RooAbsGenContext* RooEffHistProd::genContext(const RooArgSet &vars, const RooDataSet *prototype,
                                             const RooArgSet* auxProto, Bool_t verbose) const
{
   std::auto_ptr<RooArgSet> pdfObs(pdf()->getObservables(vars));
   // Return specialized generator context for RooEffHistProds that implements generation
   // in a more efficient way than can be done for generic correlated products
   assert(pdf() != 0);
   assert(eff() != 0);
   return new RooEffGenContext(*this, *pdf(), *eff(), *pdfObs.get(), prototype, auxProto, verbose);
}


//_____________________________________________________________________________
RooEffHistProd::CacheElem *RooEffHistProd::getCache(const RooArgSet* nset,
                                                    const RooArgSet* iset,
                                                    const char* rangeName,
                                                    const bool makeClone) const 
{
   Int_t sterileIndex(-1);
   CacheElem* cache = (CacheElem*) _cacheMgr.getObj(nset, iset, &sterileIndex,
                                                    RooNameReg::ptr(rangeName));
   if (cache) return cache;

   const RooEffHistProd* parent = makeClone ? static_cast<const RooEffHistProd*>
      (clone(Form("%s_clone", GetName()))) : this;
   
   cache = new CacheElem(parent, *iset, nset, rangeName);
   _cacheMgr.setObj(nset, iset, cache, RooNameReg::ptr(rangeName));

   if (makeClone) {
      cache->setClone(const_cast<RooEffHistProd*>(parent));
      cache->clone()->_fixedNormSet = &cache->intObs();
   }
   return getCache(nset, iset, rangeName, makeClone);
}

//_____________________________________________________________________________
Int_t RooEffHistProd::getAnalyticalIntegralWN(RooArgSet& allDeps, RooArgSet& analDeps, 
                                              const RooArgSet* normSet, const char* rangeName) const
{
   // Variant of getAnalyticalIntegral that is also passed the normalization set
   // that should be applied to the integrand of which the integral is request.
   // For certain operator p.d.f it is useful to overload this function rather
   // than analyticalIntegralWN() as the additional normalization information
   // may be useful in determining a more efficient decomposition of the
   // requested integral
   
   // No special handling required if a normalization set is given
   if (normSet && normSet->getSize() > 0) {
      // FIXME: here's another possible memory leak
      // if(_pdfNormSet) delete _pdfNormSet;
      _pdfNormSet = pdf()->getObservables(*normSet);
      Int_t code = _forceNumInt ? 0 : getAnalyticalIntegral(allDeps, analDeps, rangeName);
      cout << "RooEffHistProd::getAnalyticalIntegralWN " << allDeps << " " << analDeps << " "
           << (normSet ? *normSet : RooArgSet()) << " " << rangeName << endl;
      return code;
   } else if (_fixedNormSet) {    
      _pdfNormSet = _fixedNormSet;
      Int_t code = _forceNumInt ? 0 : getAnalyticalIntegral(allDeps, analDeps, rangeName);
      cout << "RooEffHistProd::getAnalyticalIntegralWN " << allDeps << " " << analDeps << " "
           << (normSet ? *normSet : RooArgSet()) << " " << rangeName << endl;
      return code;
   } else {
      // No normSet passed

      // Declare that we can analytically integrate all requested observables
      std::auto_ptr<RooArgSet> pdfObs(pdf()->getObservables(allDeps));
      analDeps.add(*pdfObs.get());
      // Construct cache with clone of p.d.f that has fixed normalization set that is passed to input pdf
      cout << "RooEffHistProd::getAnalyticalIntegralWN " << allDeps << " " << *pdfObs.get() 
           << " " << analDeps << " "
           << (normSet ? *normSet : RooArgSet()) << " " << rangeName << endl;

      getCache(pdfObs.get(), pdfObs.get(), rangeName, true);
      Int_t code = _cacheMgr.lastIndex();
      return 1 + code;
   }
}

//_____________________________________________________________________________
Int_t RooEffHistProd::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& iset,
                                            const char* rangeName) const 
{
   if (_forceNumInt) return 0;

   std::auto_ptr<RooArgSet> pdfObs(pdf()->getObservables(allVars));
   if (pdfObs->getSize() == 0) {
      return 0;
   } else {
      iset.add(*pdfObs.get());
   }

   getCache(_pdfNormSet, &iset, rangeName);
   Int_t code = _cacheMgr.lastIndex();
   return 1 + code;
}

//_____________________________________________________________________________
Double_t RooEffHistProd::analyticalIntegral(Int_t code, const char* rangeName) const 
{
   assert(code > 0);

   CacheElem* cache = static_cast<CacheElem*>(_cacheMgr.getObjByIndex(code - 1));
   if (!cache) {
      std::auto_ptr<RooArgSet> vars(getParameters(RooArgSet()));
      std::auto_ptr<RooArgSet> nset( _cacheMgr.nameSet1ByIndex(code - 1)->select(*vars));
      std::auto_ptr<RooArgSet> iset( _cacheMgr.nameSet2ByIndex(code - 1)->select(*vars));
      // const RooArgSet* normSet = _pdfNormSet ? _pdfNormSet : vars.get();
      cache = getCache(nset.get(), iset.get(), rangeName, (_pdfNormSet == 0));
   }
   return cache->getVal();
}

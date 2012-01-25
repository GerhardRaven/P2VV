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

#include <memory>
#include <algorithm>
#include <exception>

ClassImp(RooEffHistProd);

namespace {

class Exception : public exception
{
public:
   Exception(const std::string& message)
      : m_message(message)
   {

   }

   virtual ~Exception() throw()
   {
      
   }

   virtual const char* what() const throw()
   {
      return m_message.c_str();
   }

private:
   const std::string m_message;
   
};

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
            throw Exception("Got type for which a range cannot be cloned");
         }
      }
   }
}
}

//_____________________________________________________________________________
RooEffHistProd::CacheElem::CacheElem(const RooEffHistProd* parent, const RooArgSet& iset,
                                     const RooArgSet* nset, const char *rangeName)
   : _I(parent->_binboundaries.size() - 1, 0)
{
   // create a function object for the corresponding integral of the underlying PDF.
   // i.e. if we have  eps(x)f(x,y)  and we get (x,y) as allVars, 
   // construct I_i,j = int_xmin,i^xmax,j dx int dy f(x,y)
   // so that later we can do sum_i eps( (x_i+x_i+1)/2 ) * I(x_i,x_i+1)
   RooArgSet *x_ = parent->eff()->getObservables(&iset); // the subset of iset on which _eff depends
   assert(x_->getSize() < 2); // for now, we only do 1D efficiency histograms...
   const BinBoundaries& bounds = parent->_binboundaries;

   if (x_->getSize() == 1) {
      assert( *x_->first() == parent->x() );
      _trivial = false;
   } else {
      _trivial = true;
   }

   Double_t xmin = parent->x().getMin(rangeName);
   Double_t xmax = parent->x().getMax(rangeName);

   if (!trivial()) {
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
         parent->x().setRange(newRange, thisxmin, thisxmax);
         if (0 != rangeName) {
            cloneRanges(parent->observables(), iset, nset, rangeName, newRange);
         }
         _I[i] = parent->pdf()->createIntegral(iset, nset, parent->getIntegratorConfig(), newRange);
         _I[i]->printMultiline(cout, 1, kFALSE, "");
      }
   } else {
      _I[0] = parent->pdf()->createIntegral(iset, nset, parent->getIntegratorConfig(), rangeName);
      _I[0]->printMultiline(cout, 1, kFALSE, "");
   }
}

//_____________________________________________________________________________
RooArgList RooEffHistProd::CacheElem::containedArgs(Action) 
{
   // Return list of all RooAbsArgs in cache element
   RooArgList l;
   for (std::vector<RooAbsReal*>::const_iterator it = _I.begin(), end = _I.end();
        it != end; ++it) {
      l.add(**it);
   }
   return l;
}

//_____________________________________________________________________________
RooEffHistProd::CacheElem::~CacheElem() 
{
   for(std::vector<RooAbsReal*>::const_iterator it = _I.begin(), end = _I.end();
       it != end; ++it) {
      if (*it) delete *it;
   }
}

//_____________________________________________________________________________
RooEffHistProd::RooEffHistProd(const char *name, const char *title, 
                               RooAbsPdf& inPdf, RooAbsReal& inEff) :
   RooAbsPdf(name, title),
   _pdf("pdf", "pre-efficiency pdf", this, inPdf),
   _eff("eff", "efficiency function", this, inEff),
   _observables("obs", "observables in efficiency function", this),
   _cacheMgr(this, 10)
{  
   // Constructor of a a production of p.d.f inPdf with efficiency
   // function inEff.

   // to figure out what the observable is, we look at the overlap of
   // the variables of efficiency function and pdf.
  
   std::auto_ptr<RooArgSet> pdfpars(inPdf.getVariables());
   std::auto_ptr<RooArgSet> effpars(inEff.getVariables());
   std::auto_ptr<TIterator> iter(effpars->createIterator());

   RooAbsArg *effelem(0);
   while((effelem = (RooAbsArg*)iter->Next())) {
      if(pdfpars->find(effelem->GetName()))  _observables.add(*effelem);
   }

   if(_observables.getSize() == 0) {
      throw std::string("WARNING: RooEffHistProd: PDF and Efficiency function "
                        "factorise. Please use RooProd");
   } else if(_observables.getSize() > 1) {
      throw std::string("WARNING: RooEffHistProd not yet implemented for more than 1D efficiency" );
   }
  
   // an interesting hack. need to discuss with wouter. one idea: let
   // every function that is 'binned' (discrete,quantized,..) add a 
   // special binning object to its dependents.
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,32,0)
   std::list<Double_t>* binboundaries = inEff.binBoundaries(x(), x().getMin(), x().getMax());
   std::copy(binboundaries->begin(), binboundaries->end(), std::back_inserter(_binboundaries));
#else
   std::list<Double_t>* binboundaries = inEff.plotSamplingHint(x(), x().getMin(), x().getMax());
   if (binboundaries) {
      for(std::list<Double_t>::const_iterator it = binboundaries->begin();
           it != binboundaries->end(); ++it) {
         double x1 = *it++, x2 = *it;
         _binboundaries.push_back(0.5 * (x1 + x2));
      }
   }
#endif
   delete binboundaries;
}

//_____________________________________________________________________________
RooEffHistProd::RooEffHistProd(const RooEffHistProd& other, const char* name):
   RooAbsPdf(other, name),
   _pdf("pdf", this, other._pdf),
   _eff("acc", this, other._eff),
   _observables("obs", this, other._observables), 
   _cacheMgr(other._cacheMgr, this),
   _binboundaries(other._binboundaries)
{
   // Copy constructor
}

//_____________________________________________________________________________
RooEffHistProd::~RooEffHistProd() 
{
   // Destructor
}

//_____________________________________________________________________________
Double_t RooEffHistProd::evaluate() const
{
   // Calculate and return 'raw' unnormalized value of p.d.f
   return eff()->getVal() * pdf()->getVal(_normSet);
}

//_____________________________________________________________________________
RooAbsGenContext* RooEffHistProd::genContext(const RooArgSet &vars, const RooDataSet *prototype,
                                             const RooArgSet* auxProto, Bool_t verbose) const
{
   // Return specialized generator context for RooEffHistProds that implements generation
   // in a more efficient way than can be done for generic correlated products
   assert(pdf() != 0);
   assert(eff() != 0);
   return new RooEffGenContext(*this, *pdf(), *eff(), vars, prototype, auxProto, verbose);
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
RooEffHistProd::CacheElem *RooEffHistProd::getCache(const RooArgSet* nset,
                                                    const RooArgSet* iset,
                                                    const char* rangeName) const 
{
   Int_t sterileIndex(-1);
   CacheElem* cache = (CacheElem*) _cacheMgr.getObj(nset, iset, &sterileIndex,
                                                    RooNameReg::ptr(rangeName));
   if (cache) return cache;
   _cacheMgr.setObj(nset, iset, new CacheElem(this, *iset, nset, rangeName),
                    RooNameReg::ptr(rangeName));
   return getCache(nset, iset, rangeName);
}

//_____________________________________________________________________________
Int_t RooEffHistProd::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& iset,
                                            const char* rangeName) const 
{
   if (_forceNumInt) return 0;
   if (allVars.getSize() == 0) return 0;
   iset.add(allVars);

   RooArgSet *nset = &iset;
   CacheElem* cache = getCache(nset, &iset, rangeName);
   Int_t code = _cacheMgr.lastIndex();
   return 1 + code;
}

//_____________________________________________________________________________
Double_t RooEffHistProd::analyticalIntegral(Int_t code, const char* rangeName) const 
{
   assert(code > 0);

   CacheElem *cache = (CacheElem*)_cacheMgr.getObjByIndex(code-1);
   if (!cache) {
        std::auto_ptr<RooArgSet> vars(getParameters(RooArgSet()));
        std::auto_ptr<RooArgSet> iset( _cacheMgr.nameSet2ByIndex(code - 1)->select(*vars));
        CacheElem* cache = getCache(_normSet, iset.get(), rangeName);
   }

   Double_t xmin = x().getMin(rangeName), xmax = x().getMax(rangeName);

   // make sure the range is contained within the binboundaries...
   assert(_binboundaries.size() > 1);
   assert(xmin <= xmax);
   assert(xmin >= _binboundaries.front());
   assert(_binboundaries.back() - xmax > -1e-10);

   if (cache->trivial()) {// no integral over efficiency dependant...
      return eff()->getVal() * cache->getVal();
   }

   double xorig = x().getVal();
   Bool_t origState = inhibitDirty();

   double sum(0);
   for(unsigned int i = 0; i < _binboundaries.size() - 1; ++i) {
      Double_t low = _binboundaries[i];
      Double_t high = _binboundaries[i + 1];
      if (high < xmin) continue; // not there yet...
      if (low  > xmax) break;    // past the requested interval...
      Double_t thisxmin = std::max(low,       xmin);
      Double_t thisxmax = std::min(high, xmax);
      if (thisxmin >= thisxmax) continue;

      setDirtyInhibit(kTRUE);
      x().setVal(0.5 * (thisxmin + thisxmax)); // get the efficiency for this bin
      double eps = eff()->getVal();
      x().setVal(xorig);  // and restore the original value ASAP...
      eff()->getVal() ; // restore original state
      setDirtyInhibit(origState) ;	

      // Passing the bin index is not very nice, but I prefer not having to loop twice.
      sum += eps * cache->getVal(i);
   }
   eff()->getVal() ; // restore original state
   setDirtyInhibit(origState) ;	
   return  sum;
}

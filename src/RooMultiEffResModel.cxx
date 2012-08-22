/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooMultiEffResModel.cxx 44982 2012-07-10 08:36:13Z moneta $
 * Authors:                                                                  *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Class RooMultiEffResModel implements a RooResolutionModel that models a EffResian
// distribution. Object of class RooMultiEffResModel can be used
// for analytical convolutions with classes inheriting from RooAbsAnaConvPdf
// END_HTML
//

#include <memory>
#include <sstream>

#include "RooFit.h"
#include "Riostream.h"
#include "RooMultiEffResModel.h"
#include "RooRealConstant.h"
#include "RooCustomizer.h"
#include "RooAddition.h"
#include "RooSuperCategory.h"

namespace {
   TString makeName(const char* name, const RooArgSet& terms ) {
      TString pname;
      pname = name;
      pname.Append("_");
      RooFIter i = terms.fwdIterator();
      RooAbsArg *arg;
      bool first = true;
      while((arg = i.next())) {
         if (first) {
            first= false;
         } else {
            pname.Append("_X_");
         }
         pname.Append(arg->GetName());
      }
      return pname;
   }

   RooSuperCategory* makeSuper(const char* name, const RooArgSet& _catVars) {
      TString catName = makeName(name, _catVars );
      return new RooSuperCategory(catName.Data(), catName, _catVars);
   }

   using std::string;
   using std::stringstream;
   using std::vector;
   using std::list;
   using std::map;
   using std::pair;
}

ClassImp(RooMultiEffResModel) 
;

//_____________________________________________________________________________
RooMultiEffResModel::CacheElem::~CacheElem()
{
   for (Integrals::const_iterator it = _I.begin(), end = _I.end(); it != end;
        ++it) {
      delete it->second;
   }
}

//_____________________________________________________________________________
RooArgList RooMultiEffResModel::CacheElem::containedArgs(Action) 
{
   // Return list of all RooAbsArgs in cache element
   RooArgList l;
   for (Integrals::const_iterator it = _I.begin(), end = _I.end(); it != end;
        ++it) {
      l.add(*(it->second));
   }
   l.add(_iset);
   return l;
}

//_____________________________________________________________________________
RooMultiEffResModel::CacheElem::CacheElem(const HistEntries& entries,
                                          const RooArgSet& iset, const TNamed* rangeName)
   : _iset(iset)
{
   bool cats = false;
   RooArgSet categories = entries.begin()->second->categories();
   RooArgSet observables(iset);
   RooAbsArg* cat = categories.first();
   if (iset.contains(*cat)) {
      // Integral over only categories
      cats = true;
      observables.remove(categories);
   }

   for(HistEntries::const_iterator it = entries.begin(), end = entries.end();
       it != end; ++it) {
      Int_t index = it->first;
      const HistEntry* entry = it->second;
      if (cats && iset.getSize() == categories.getSize()) {
         stringstream s;
         string name;
         s << "integral_" << index;
         s >> name;
         _I.insert(make_pair(index, new RooConstVar(name.c_str(), name.c_str(), 1.)));
      } else {
         RooAbsReal* I = entry->efficiency()->createIntegral(observables, RooNameReg::str(rangeName));
         _I.insert(make_pair(index, I));
      }
   }
}

//_____________________________________________________________________________
Double_t RooMultiEffResModel::CacheElem::getVal(const Int_t index) const
{
   Integrals::const_iterator it = _I.find(index);
   assert(it != _I.end());
   double val = it->second->getVal(_iset);
   // cout << "CacheElem::getVal " << it->second->GetName() << " = " << val << endl;
   return val;
}

//_____________________________________________________________________________
RooMultiEffResModel::RooMultiEffResModel(const char *name, const char *title,
                                         std::vector<HistEntry*> entries)
   : RooResolutionModel(name,title, (*entries.begin())->efficiency()->convVar()),
     _binboundaries(0),
     _prodGenCode(0),
     _super(0),
     _cacheMgr(this, 10)
{  
   RooArgSet observables;
   RooArgSet categories;
   for(vector<HistEntry*>::const_iterator it = entries.begin(),
          end = entries.end(); it != end; ++it) {
      HistEntry* entry = *it;
      if (observables.getSize() == 0) {
         observables.add(*(entry->efficiency()->observables()));
         categories.add(entry->categories());
      } else {
         assert(observables.equals(*(entry->efficiency()->observables())));
         assert(categories.equals(entry->categories()));
      }
   }
   
   // Observable
   RooAbsRealLValue* x = dynamic_cast<RooAbsRealLValue*>(observables.first());

   // For the binnings, we just assume that they are "matched" by the user.
   unsigned int most = 0;
   for(vector<HistEntry*>::const_iterator it = entries.begin(),
          end = entries.end(); it != end; ++it) {
      RooAbsReal* eff = (*it)->efficiency()->efficiency();
      auto_ptr<BinBoundaries> bounds(eff->binBoundaries(*x, x->getMin(), x->getMax()));
      if (!bounds.get()) {
         continue;
      } else if (bounds->size() > most) {
         if (_binboundaries) delete _binboundaries;
         _binboundaries = bounds.release();
      }
   }

   cout << "RooMultHistEfficiency::ctor(): bins [";
   for(list<Double_t>::const_iterator it = _binboundaries->begin();
       it != _binboundaries->end(); ++it) {
      if (it != _binboundaries->begin()) cout << " ";
      cout << *it;
   }
   cout << "]" << endl;

   // Build entries.
   _super = makeSuper(GetName(), categories);

   typedef map<RooAbsCategory*, string> categories_t;
   categories_t signal;

   vector<HistEntry*> ownedEntries;
   for(vector<HistEntry*>::const_iterator it = entries.begin(),
          end = entries.end(); it != end; ++it) {
      HistEntry* entry = new HistEntry(**it);
      entry->setParent(this);
      ownedEntries.push_back(entry);
   }
   TString current = _super->getLabel();

   for(vector<HistEntry*>::const_iterator it = ownedEntries.begin(),
          end = ownedEntries.end(); it != end; ++it) {
      HistEntry* entry = *it;
      entry->select();
      Int_t index = _super->getIndex();
      pair<HistEntries::iterator, bool> r = _entries.insert(make_pair(index, entry));
      assert(r.second);
      _intVals.insert(make_pair(index, 0.));
   }
   _super->setLabel(current.Data());
}

//_____________________________________________________________________________
RooMultiEffResModel::RooMultiEffResModel(const RooMultiEffResModel& other, const char* name) 
   : RooResolutionModel(other,name),
     _intVals(other._intVals),
     _prodGenObs(other._prodGenObs),
     _prodGenCode(other._prodGenCode),
     _levels(other._levels),
     _cacheMgr(other._cacheMgr,this)
{
   // Copy constructor
   _binboundaries = new BinBoundaries(*other._binboundaries);
   _super = new RooSuperCategory(*other._super);

   for (HistEntries::const_iterator it = other._entries.begin(), end = other._entries.end();
        it != end; ++it) {
      HistEntry* entry = new HistEntry(*(it->second), this);
      _entries.insert(make_pair(it->first, entry));
   }
}

//_____________________________________________________________________________
RooMultiEffResModel::~RooMultiEffResModel()
{
   // Destructor
   if (_binboundaries) delete _binboundaries;
   if (_super) delete _super;
}

//_____________________________________________________________________________
RooMultiEffResModel* 
RooMultiEffResModel::convolution(RooFormulaVar* inBasis, RooAbsArg* owner) const
{
   // Instantiate a clone of this resolution model representing a convolution with given
   // basis function. The owners object name is incorporated in the clones name
   // to avoid multiple convolution objects with the same name in complex PDF structures.
   // 
   // Note: The 'inBasis' formula expression must be a RooFormulaVar that encodes the formula
   // in the title of the object and this expression must be an exact match against the
   // implemented basis function strings (see derived class implementation of method basisCode()
   // for those strings

   //IMPLEMENT!!!

   // Check that primary variable of basis functions is our convolution variable  
   if (inBasis->getParameter(0) != x.absArg()) {
      coutE(InputArguments) << "RooMultiEffResModel::convolution(" << GetName() << "," << this
                            << ") convolution parameter of basis function and PDF don't match" << endl
                            << "basis->findServer(0) = " << inBasis->findServer(0) << endl
                            << "x.absArg()           = " << x.absArg() << endl ;
      return 0 ;
   }

   if (basisCode(inBasis->GetTitle())==0) {
      coutE(InputArguments) << "RooMultiEffResModel::convolution(" << GetName() << "," << this
                            << ") basis function '" << inBasis->GetTitle() << "' is not supported." << endl ;
      return 0 ;
   }

   vector<HistEntry*> entries;
   vector<RooResolutionModel*> models;
   for (HistEntries::const_iterator it = _entries.begin(), end = _entries.end();
        it != end; ++it) {
      RooEffResModel *conv = it->second->efficiency()->convolution(inBasis, owner);
      models.push_back(conv);

      HistEntry* entry = new HistEntry(*(it->second), const_cast<RooMultiEffResModel*>(this));
      entry->setEfficiency(conv);
      entries.push_back(entry);
   }

   TString newName(GetName());
   newName.Append("_conv_") ;
   newName.Append(inBasis->GetName()) ;
   newName.Append("_[") ;
   newName.Append(owner->GetName()) ;
   newName.Append("]") ;

   TString newTitle(GetTitle()) ;
   newTitle.Append(" convoluted with basis function ") ;
   newTitle.Append(inBasis->GetName()) ;

   RooMultiEffResModel *effConv = new RooMultiEffResModel(newName, newTitle, entries);
   for (vector<RooResolutionModel*>::iterator it = models.begin(), end = models.end();
        it != end; ++it) {
      effConv->addOwnedComponents(**it);
   }
   effConv->changeBasis(inBasis) ;

   for (vector<HistEntry*>::const_iterator it = entries.begin(),
           end = entries.end(); it != end; ++it) {
      delete *it;
   }

   return effConv ;
}

//_____________________________________________________________________________
Int_t RooMultiEffResModel::basisCode(const char* name) const 
{ 
   return _entries.begin()->second->efficiency()->basisCode(name);
} 


//_____________________________________________________________________________
Double_t RooMultiEffResModel::evaluate() const 
{  
   bool onlyCats = false;
   bool onlyVars = false;

   // Calculate the raw value of this p.d.f
   RooArgSet categories(_entries.begin()->second->categories());
   RooAbsArg* cat = categories.first();
   if (_normSet && _normSet->getSize() == categories.getSize() && 
       _normSet->contains(*cat)) {
      onlyCats = true;
   } else if (_normSet && !_normSet->contains(*cat)) {
      onlyVars = true;
   }

   double val = 0;
   for (HistEntries::const_iterator it = _entries.begin(), end = _entries.end();
        it != end; ++it) {
      if (it->second->thisEntry()) {
         if (onlyCats) {
            val = it->second->relative()->getVal();
         } else if (onlyVars) {
            val = it->second->relative()->getVal() * it->second->efficiency()->getVal();
         } else {
            val = it->second->relative()->getVal() * it->second->efficiency()->getVal();
         }
         // cout << "RooMultiHistEfficiency::evaluate " 
         //      << it->second->efficiency()->GetName() << " case = " << code << " " 
         //      << " norm " << (_normSet ? *_normSet : RooArgSet()) << " = " << val << endl;
      }
   }
   return val;
}


//_____________________________________________________________________________
Bool_t RooMultiEffResModel::forceAnalyticalInt(const RooAbsArg& /*dep*/) const
{
   // Force RooRealIntegral to offer all observables for internal integration
   return kTRUE ;
}

//_____________________________________________________________________________
Int_t RooMultiEffResModel::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
   if (_forceNumInt) return 0;

   RooArgSet categories;
   for (HistEntries::const_iterator it = _entries.begin(), end = _entries.end();
        it != end; ++it) {
      categories.add(it->second->categories());
   }

   bool all = true;
   bool none = true;
   RooFIter iter = categories.fwdIterator();
   while (RooAbsArg* cat = iter.next()) {
      if (!allVars.find(*cat)) {
         all &= false;
      } else {
         none &= false;
      }
   }
   if (!(all ^ none)) {
      return 0;
   }

   // we always do things ourselves -- actually, always delegate further down the line ;-)
   analVars.add(allVars);
   // The underlying pdf does not depend on the categories
   // analVars.remove(categories);

   bool vars = (!all && analVars.getSize() != 0) || (all && analVars.getSize() > categories.getSize());

   // check if we already have integrals for this combination of factors
   Int_t sterileIndex(-1);
   CacheElem* cache = (CacheElem*) _cacheMgr.getObj(&analVars, &analVars, &sterileIndex,
                                                    RooNameReg::ptr(rangeName));
   if (cache != 0) {
      Int_t code = _cacheMgr.lastIndex();
      code = (code + 1) << 1;
      code |= vars;
      code = code << 1;
      code |= all;
      return code;
   }

   // we don't, so we make it right here....
   RooArgSet pdfVars(analVars);
   pdfVars.remove(categories);
   cache = new CacheElem(_entries, analVars, RooNameReg::ptr(rangeName));
   
   Int_t code = _cacheMgr.setObj(&analVars, &analVars, cache, RooNameReg::ptr(rangeName));
   // cout << "getAnalyticalIntegral cache code = " << code << endl;
   code = (code + 1) << 1;
   code |= vars;
   code = code << 1;
   code |= all;
   // cout << "getAnalyticalIntegral return code = " << code << endl;
   return code;
}

//_____________________________________________________________________________
Double_t RooMultiEffResModel::analyticalIntegral(Int_t code, const char* rangeName) const 
{
   assert(code > 0);
   // Calculate integral internally from appropriate integral cache
   // note: rangeName implicit encoded in code: see _cacheMgr.setObj in getPartIntList...
   bool cats = code & 0x1;
   bool vars = code & 0x2;
   Int_t cacheCode = (code >> 2);

   // cout << "analyticalIntegral code = " << code << " cachecode = " << cacheCode << endl;
 
   CacheElem *cache = static_cast<CacheElem*>(_cacheMgr.getObjByIndex(cacheCode - 1));
   
   if (cache==0) {
      // cache got sterilized, trigger repopulation of this slot, then try again...
      auto_ptr<RooArgSet> params(getParameters(RooArgSet()));
      auto_ptr<RooArgSet> iset(_cacheMgr.nameSet2ByIndex(cacheCode - 1)->select(*params));
      RooArgSet dummy;
      Int_t code2 = getAnalyticalIntegral(*iset, dummy, rangeName);
      // must have revived the right (sterilized) slot...
      // cout << "code = " << code << " code2 = " << code2 << endl;
      assert(code == code2);
      return analyticalIntegral(code2, rangeName);
   }
   assert(cache != 0);

   double sum = 0;
   for (HistEntries::const_iterator it = _entries.begin(), end = _entries.end();
        it != end; ++it) {
      if (!cats && vars && it->second->thisEntry()) {
         double entry = it->second->relative()->getVal() * cache->getVal(it->first);
         // cout << "AnalyticalIntegral: entry = " << entry << endl;
         return entry;
      } else if (cats) {
         sum += it->second->relative()->getVal() * cache->getVal(it->first);
      }
   }
   // cout << "AnalyticalIntegral: sum = " << sum << endl;
   return sum;
}

//_____________________________________________________________________________
Int_t RooMultiEffResModel::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
   return 0 ; // For now... problem is that RooGenConv assumes it can just add resolution & physics for conv var...
}

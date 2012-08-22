//////////////////////////////////////////////////////////////////////////////
// 
// BEGIN_HTML
// RooMultiHistEfficiency is a PDF helper class to fit efficiencies parameterized
// by a supplied functions F_i.
// 
// Given a dataset with a categories C_i that determines if a given
// event is accepted or rejected for the efficiency to be measured,
// this class evaluates the product of G_i where G_i is F_i if C is
// 'accept' and as (1-F_i) if C_i is 'reject'. Values of F_i below 0
// and above 1 are clipped. F may have an arbitrary number of
// dependents and parameters.
// END_HTML
//

#include <memory>
#include <algorithm>
#include <sstream>

#include <RooFit.h>
#include <RooRandom.h>
#include <RooStreamParser.h>
#include <RooGenContext.h>
#include <RooArgList.h>
#include <RooCategory.h>
#include <RooSuperCategory.h>
#include <RooRealVar.h>
#include <RooAbsReal.h>
#include <RooConstVar.h>

#include <RooMultiHistEfficiency.h>

ClassImp(RooMultiHistEfficiency);

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

   using std::stringstream;
   using std::make_pair;
   using std::string;
   using std::cout;
   using std::endl;
}

//_____________________________________________________________________________
RooMultiHistEfficiency::CacheElem::CacheElem(const HistEntries& entries, 
                                             const RooArgSet& iset,
                                             const char* rangeName)
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
         RooAbsReal* I = entry->efficiency()->createIntegral(observables, rangeName);
         _I.insert(make_pair(index, I));
      }
   }
}

//_____________________________________________________________________________
Double_t RooMultiHistEfficiency::CacheElem::getVal(const Int_t index) const
{
   Integrals::const_iterator it = _I.find(index);
   assert(it != _I.end());
   double val = it->second->getVal(_iset);
   // cout << "CacheElem::getVal " << it->second->GetName() << " = " << val << endl;
   return val;
}

//_____________________________________________________________________________
RooArgList RooMultiHistEfficiency::CacheElem::containedArgs(Action) 
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
RooMultiHistEfficiency::CacheElem::~CacheElem() 
{
   for (Integrals::const_iterator it = _I.begin(), end = _I.end(); it != end;
        ++it) {
      delete it->second;
   }
}

//_____________________________________________________________________________
RooMultiHistEfficiency::RooMultiHistEfficiency
(const char *name, const char *title, std::vector<HistEntry*> entries)
   : RooAbsPdf(name, title),
     _binboundaries(0),
     _prodGenCode(0),
     _super(0),
     _cacheMgr(this, 10)
{
   RooArgSet observables;
   RooArgSet categories;
   for(std::vector<HistEntry*>::const_iterator it = entries.begin(),
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
   for(std::vector<HistEntry*>::const_iterator it = entries.begin(),
          end = entries.end(); it != end; ++it) {
      RooAbsReal* eff = (*it)->efficiency()->efficiency();
      std::auto_ptr<BinBoundaries> bounds(eff->binBoundaries(*x, x->getMin(), x->getMax()));
      if (!bounds.get()) {
         continue;
      } else if (bounds->size() > most) {
         if (_binboundaries) delete _binboundaries;
         _binboundaries = bounds.release();
      }
   }

   cout << "RooMultHistEfficiency::ctor(): bins [";
   for(std::list<Double_t>::const_iterator it = _binboundaries->begin();
       it != _binboundaries->end(); ++it) {
      if (it != _binboundaries->begin()) cout << " ";
      cout << *it;
   }
   cout << "]" << endl;

   // Build entries.
   _super = makeSuper(GetName(), categories);

   typedef std::map<RooAbsCategory*, std::string> categories_t;
   categories_t signal;

   std::vector<HistEntry*> ownedEntries;
   for(std::vector<HistEntry*>::const_iterator it = entries.begin(),
          end = entries.end(); it != end; ++it) {
      HistEntry* entry = new HistEntry(**it);
      entry->setParent(this);
      ownedEntries.push_back(entry);
   }
   TString current = _super->getLabel();

   for(std::vector<HistEntry*>::const_iterator it = ownedEntries.begin(),
          end = ownedEntries.end(); it != end; ++it) {
      HistEntry* entry = *it;
      entry->select();
      Int_t index = _super->getIndex();
      std::pair<HistEntries::iterator, bool> r = _entries.insert(make_pair(index, entry));
      assert(r.second);
      _intVals.insert(make_pair(index, 0.));
   }
   _super->setLabel(current.Data());
}

//_____________________________________________________________________________
RooMultiHistEfficiency::RooMultiHistEfficiency(const RooMultiHistEfficiency& other, const char* name)
   : RooAbsPdf(other, name),
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
RooMultiHistEfficiency::~RooMultiHistEfficiency() 
{
   // Destructor
   if (_binboundaries) delete _binboundaries;
   if (_super) delete _super;
}

//_____________________________________________________________________________
std::list<Double_t>* RooMultiHistEfficiency::binBoundaries
(RooAbsRealLValue& obs, Double_t min, Double_t max) const
{
   assert (min >= _binboundaries->front());
   assert (max <= _binboundaries->back());
   std::list<Double_t>* bounds = new std::list<Double_t>;
   std::list<Double_t>::const_iterator start, end = _binboundaries->end();
   for(std::list<Double_t>::const_iterator it = _binboundaries->begin();
       it != end; ++it) {
      start = it;
      ++start;
      if (*it <= min && *start > min) {
         bounds->push_back(min);
         break;
      }
   }
   for(; start != end; ++start) {
      if (*start >= max) {
         bounds->push_back(max);
         break;
      }
      bounds->push_back(*start);
   }
   return bounds;
}

//_____________________________________________________________________________
RooAbsGenContext* RooMultiHistEfficiency::genContext(const RooArgSet &vars, const RooDataSet *prototype,
                                                     const RooArgSet* auxProto, Bool_t verbose) const
{
   // Return specialized generator context for RooEffHistProds that implements generation
   // in a more efficient way than can be done for generic correlated products
   const HistEntry* entry = _entries.begin()->second;
   const RooArgSet* observables = entry->efficiency()->observables();
   return new RooGenContext(*this, vars, prototype, auxProto, verbose, observables);
}

//_____________________________________________________________________________
Int_t RooMultiHistEfficiency::getGenerator
(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const
{
   bool all = true;
   bool none = true;

   _prodGenObs.removeAll();

   RooArgSet categories;
   for (HistEntries::const_iterator it = _entries.begin(), end = _entries.end();
        it != end; ++it) {
      categories.add(it->second->categories());
   }
   RooFIter iter = categories.fwdIterator();
   while (RooAbsArg* cat = iter.next()) {
      if (!directVars.find(*cat)) {
         all &= false;
      } else {
         none &= false;
      }
   }
   if (!(all ^ none)) {
      return 0;
   }

   RooArgSet testVars(directVars);
   testVars.remove(categories);

   Int_t prodGenCode = 0;
   RooArgSet genVars;
   for (HistEntries::const_iterator it = _entries.begin(), end = _entries.end();
        it != end; ++it) {
      const RooEffHistProd* efficiency = it->second->efficiency();
      if (genVars.getSize() == 0) {
         prodGenCode = efficiency->getGenerator(testVars, genVars, staticInitOK);
      } else {
         RooArgSet prodGenVars;
         Int_t code = efficiency->getGenerator(testVars, prodGenVars, staticInitOK);
         assert(prodGenVars.equals(genVars) && prodGenCode == code);
      }
   }

   generateVars.add(genVars);
   _prodGenCode = prodGenCode;
   _prodGenObs.add(testVars);
   if (none) {
      return 1;
   } else {
      iter = categories.fwdIterator();
      while (RooAbsArg* cat = iter.next()) {
         generateVars.add(*cat);
      }
      return 2;
   }
}

//_____________________________________________________________________________
void RooMultiHistEfficiency::initGenerator(Int_t code)
{
   // Forward one-time initialization call to component generation initialization
   // methods.
   for (HistEntries::iterator it = _entries.begin(), end = _entries.end();
        it != end; ++it) {
      it->second->efficiency()->initGenerator(_prodGenCode);
   }
   if (code == 1) {
      return;
   }
   _levels.clear();

   RooArgSet categories;
   for (HistEntries::const_iterator it = _entries.begin(), end = _entries.end();
        it != end; ++it) {
      categories.add(it->second->categories());
   }
   _super->recursiveRedirectServers(categories);

   // RooSuperCategory* super = dynamic_cast<RooSuperCategory*>(_super->absArg());
   std::auto_ptr<TIterator> superIter(_super->MakeIterator());

   TString current = _super->getLabel();
   while (TObjString* label = static_cast<TObjString*>(superIter->Next())) {
      _super->setLabel(label->String());
      Int_t index = _super->getIndex();
      HistEntries::const_iterator it = _entries.find(index);
      if (it == _entries.end()) {
         // Skip the combination for which there is no shape (all false).
         continue;
      }
      double n = it->second->relative()->getVal();
      if (!_levels.empty()) n += _levels.back().first; // cumulative
      cxcoutD(Generation) << "RooMultiHistEfficiency creating sampler for " << _prodGenObs
                          << " given " << categories
                          << " = "  << label->String() << " (level = " << n << ")" << endl;
      _levels.push_back(make_pair(n, label->String()));
   }

   // Normalise just in case, but it should be a noop.
   for (Levels::iterator i = _levels.begin(); i != _levels.end(); ++i) 
      i->first /= _levels.back().first; // normalize
   _super->setLabel(current);
}

//_____________________________________________________________________________
void RooMultiHistEfficiency::generateEvent(Int_t code)
{
   Double_t r = RooRandom::uniform();

   // RooSuperCategory* super = dynamic_cast<RooSuperCategory*>(_super->absArg());
   std::auto_ptr<TIterator> superIter(_super->MakeIterator());

   Levels::const_iterator itLevel = _levels.begin();
   // find the right generator, and generate categories at the same time...
   while (itLevel != _levels.end() && itLevel->first < r) {
      ++itLevel;
   }

   // this assigns the categories.
   _super->setLabel(itLevel->second);

   Int_t index = _super->getIndex();
   HistEntries::const_iterator itEntry = _entries.find(index);
   assert(itEntry != _entries.end());

   // now that've assigned the categories, we can use the 'real' samplers
   // which are conditional on the categories.
   itEntry->second->efficiency()->generateEvent(_prodGenCode);
}

//_____________________________________________________________________________
Bool_t	RooMultiHistEfficiency::forceAnalyticalInt(const RooAbsArg& var) const
{
   // We forward everything anyway so always return true.
   return true;
}

//_____________________________________________________________________________
Int_t RooMultiHistEfficiency::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& iset,
                                                    const char* rangeName) const 
{
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
   iset.add(allVars);
   // The underlying pdf does not depend on the categories
   // iset.remove(categories);

   bool vars = (!all && iset.getSize() != 0) || (all && iset.getSize() > categories.getSize());

   // check if we already have integrals for this combination of factors
   Int_t sterileIndex(-1);
   CacheElem* cache = (CacheElem*) _cacheMgr.getObj(&iset, &iset, &sterileIndex,
                                                    RooNameReg::ptr(rangeName));
   if (cache!=0) {
      Int_t code = _cacheMgr.lastIndex();
      code = (code + 1) << 1;
      code |= vars;
      code = code << 1;
      code |= all;
      return code;
   }

   // we don't, so we make it right here....
   RooArgSet pdfVars(iset);
   pdfVars.remove(categories);
   cache = new CacheElem(_entries, iset, rangeName);
   
   Int_t code = _cacheMgr.setObj(&iset, &iset, cache, RooNameReg::ptr(rangeName));
   // cout << "getAnalyticalIntegral cache code = " << code << endl;
   code = (code + 1) << 1;
   code |= vars;
   code = code << 1;
   code |= all;
   // cout << "getAnalyticalIntegral return code = " << code << endl;
   return code;
}

//_____________________________________________________________________________
Double_t RooMultiHistEfficiency::analyticalIntegral(Int_t code, const char* rangeName) const 
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
      std::auto_ptr<RooArgSet> params(getParameters(RooArgSet()));
      std::auto_ptr<RooArgSet> iset(_cacheMgr.nameSet2ByIndex(cacheCode - 1)->select(*params));
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

// //_____________________________________________________________________________
// Double_t RooMultiHistEfficiency::getValV(const RooArgSet* ns) const 
// {  
//    // Return p.d.f. value normalized over given set of observables
//    // cout << "RooMultiHistEfficiency::getValV " << (ns ? *ns : RooArgSet()) << endl;
//    // for (HistEntries::const_iterator it = _entries.begin(), end = _entries.end();
//    //      it != end; ++it) {
//    //    it->second->efficiency()->setNormSet(ns);
//    // }
//    return RooAbsPdf::getValV(ns);
// }

//_____________________________________________________________________________
Double_t RooMultiHistEfficiency::evaluate() const
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

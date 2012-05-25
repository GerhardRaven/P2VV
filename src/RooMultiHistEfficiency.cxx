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

#include <RooFit.h>
#include <RooRandom.h>
#include <RooStreamParser.h>
#include <RooGenContext.h>
#include <RooArgList.h>
#include <RooCategory.h>
#include <RooSuperCategory.h>
#include <RooRealVar.h>
#include <RooAbsReal.h>

#include <RooMultiHistEfficiency.h>

ClassImp(RooMultiHistEfficiency);

namespace {
   const char *makeName(const char* name, const RooArgSet& terms ) {
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
      return pname.Data();
   }

   RooSuperCategory* makeSuper(const char* name, const RooArgSet& _catVars) {
      const char *catName = makeName(name, _catVars );
      return new RooSuperCategory(catName, catName, _catVars);
   }
}

//_____________________________________________________________________________
RooMultiHistEfficiency::RooMultiHistEfficiency
(const char *name, const char *title, std::vector<MultiHistEntry> entries) :
   RooAbsPdf(name, title),
   _binboundaries(0),
   _prodGenCode(0),
   _super(0)
{  
   // Construct an N+1 dimensional efficiency p.d.f from an N-dimensional efficiency
   // function and a category cat with two states (0,1) that indicate if a given
   // event should be counted as rejected or accepted respectively
   
   RooArgSet observables;
   RooArgSet categories;
   for(std::vector<MultiHistEntry>::const_iterator it = entries.begin(),
          end = entries.end(); it != end; ++it) {
      if (observables.getSize() == 0) {
         observables.add(*(it->effProd().observables()));
         categories.add(it->categories());
      } else {
         assert(observables.equals(*(it->effProd().observables())));
         assert(categories.equals(it->categories()));
      }
   }
   
   // Observable
   RooAbsRealLValue* x = dynamic_cast<RooAbsRealLValue*>(observables.first());

   // For the binnings, we just assume that they are "matched" by the user.
   unsigned int most = 0;
   for(std::vector<MultiHistEntry>::const_iterator it = entries.begin(),
          end = entries.end(); it != end; ++it) {
      RooAbsReal* eff = it->effProd().efficiency();
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

   std::vector<MultiHistEntry*> ownedEntries;
   for(std::vector<MultiHistEntry>::const_iterator it = entries.begin(),
          end = entries.end(); it != end; ++it) {
      MultiHistEntry* entry = new MultiHistEntry(*it);
      entry->setParent(this);
      ownedEntries.push_back(entry);
   }
   TString current = _super->getLabel();
   for(std::vector<MultiHistEntry*>::const_iterator it = ownedEntries.begin(),
          end = ownedEntries.end(); it != end; ++it) {
      MultiHistEntry* entry = *it;
      entry->select();
      Int_t index = _super->getIndex();
      std::pair<HistEntries::iterator, bool> r = _entries.insert(make_pair(index, entry));
      assert(r.second);
   }
   _super->setLabel(current.Data());
}

//_____________________________________________________________________________
RooMultiHistEfficiency::RooMultiHistEfficiency(const RooMultiHistEfficiency& other, const char* name) : 
   RooAbsPdf(other, name),
   _prodGenObs(other._prodGenObs),
   _prodGenCode(other._prodGenCode),
   _levels(other._levels),
   _pdfIntCode(other._pdfIntCode),
   _pdfIntObs(other._pdfIntObs)
{
   // Copy constructor
   _binboundaries = new BinBoundaries(*other._binboundaries);
   _super = new RooSuperCategory(*other._super);

   for (HistEntries::const_iterator it = other._entries.begin(), end = other._entries.end();
        it != end; ++it) {
      MultiHistEntry* entry = new MultiHistEntry(*(it->second), this);
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
   const MultiHistEntry* entry = _entries.begin()->second;
   const RooArgSet* observables = entry->effProd().observables();
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
      const RooEffHistProd& effProd = it->second->effProd();
      if (genVars.getSize() == 0) {
         prodGenCode = effProd.getGenerator(testVars, genVars, staticInitOK);
      } else {
         RooArgSet prodGenVars;
         Int_t code = effProd.getGenerator(testVars, prodGenVars, staticInitOK);
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
      it->second->effProd().initGenerator(_prodGenCode);
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
      double n = it->second->relative();
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
   itEntry->second->effProd().generateEvent(_prodGenCode);
}

//_____________________________________________________________________________
Bool_t	RooMultiHistEfficiency::forceAnalyticalInt(const RooAbsArg& var) const
{
   assert(!_entries.empty());
   
   const RooArgSet* observables =  _entries.begin()->second->effProd().observables();
   return observables->contains(var);
}

//_____________________________________________________________________________
Int_t RooMultiHistEfficiency::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& iset,
                                                    const char* rangeName) const 
{
   cout << "RHME::getAnalyticalIntegral " << allVars << endl;

   RooArgSet intObs;

   int prodCode = 0;
   for (HistEntries::const_iterator it = _entries.begin(), end = _entries.end();
        it != end; ++it) {
      if (intObs.getSize() == 0) {
         prodCode = it->second->effProd().getAnalyticalIntegral(allVars, intObs, rangeName);
         cout << "EffProd can integrate " << intObs << " out of " << allVars << endl;
      } else {
         RooArgSet prodIntObs;
         Int_t code = it->second->effProd().getAnalyticalIntegral(allVars, prodIntObs, rangeName);
         assert(prodIntObs.equals(intObs));
         assert(code == prodCode);
      }
   }
   if (prodCode == 0) {
      return 0;
   } else {
      _pdfIntCode = prodCode;
      iset.add(intObs);
      return 1;
   }
}

//_____________________________________________________________________________
Double_t RooMultiHistEfficiency::analyticalIntegral(Int_t code, const char* rangeName) const 
{
   assert(code > 0);

   Double_t val = 0;
   for (HistEntries::const_iterator it = _entries.begin(), end = _entries.end();
        it != end; ++it) {
      const MultiHistEntry* entry = it->second;
      val += entry->relative() * entry->effProd().analyticalIntegral(_pdfIntCode, rangeName);
   }
   return val;
}

//_____________________________________________________________________________
Double_t RooMultiHistEfficiency::evaluate() const
{
   // Calculate the raw value of this p.d.f
   for (HistEntries::const_iterator it = _entries.begin(), end = _entries.end();
        it != end; ++it) {
      if (it->second->thisEntry()) {
         double val = it->second->effVal();
         // cout << "Entry for " << it->second->effProd().GetName() << " = " << val << endl;
         return val;
      }
   }
   throw std::string("The efficiency for a state is missing");
   return 0;
}

//_____________________________________________________________________________
MultiHistEntry::MultiHistEntry()
   : m_rawEff(0), m_rawRel(0), m_effProd(0), m_relative(0),
     m_index(0)
{
}

//_____________________________________________________________________________
MultiHistEntry::MultiHistEntry(const std::map<RooAbsCategory*, std::string>& categories,
                               RooEffHistProd* effProd, RooAbsReal* relative)
   : m_rawCats(categories), m_rawEff(effProd), m_rawRel(relative),
     m_effProd(0), m_relative(0), m_index(0)
{

}

//_____________________________________________________________________________
MultiHistEntry::MultiHistEntry(const MultiHistEntry& other, RooMultiHistEfficiency* parent)
   : m_rawCats(other.m_rawCats), m_rawEff(other.m_rawEff), m_rawRel(other.m_rawRel),
     m_index(other.m_index)
{
   if (!other.m_effProd) {
      m_effProd = 0;
      m_relative = 0;
      return;
   }
   m_effProd = new RooRealProxy(other.m_effProd->GetName(), parent, *other.m_effProd);
   m_relative = new RooRealProxy(other.m_relative->GetName(), parent, *other.m_relative);

   for (std::map<RooCategoryProxy*, std::string>::const_iterator it = other.m_categories.begin(),
           end = other.m_categories.end(); it != end; ++it) {
      m_categories.insert(make_pair(new RooCategoryProxy(it->first->GetName(), parent, *(it->first)),
                                    it->second));
   }

}

//_____________________________________________________________________________
MultiHistEntry::~MultiHistEntry()
{
   for (std::map<RooCategoryProxy*, std::string>::const_iterator it = m_categories.begin(),
           end = m_categories.end(); it != end; ++it) {
      if (it->first) delete it->first;
   }
   m_categories.clear();
   if (m_effProd) delete m_effProd;
   if (m_relative) delete m_relative;
}

//_____________________________________________________________________________
bool MultiHistEntry::thisEntry() const
{
   bool r = true;
   for (std::map<RooCategoryProxy*, std::string>::const_iterator it = m_categories.begin(),
           end = m_categories.end(); it != end; ++it) {
      const RooAbsCategory* cat = static_cast<const RooAbsCategory*>(it->first->absArg());
      if (strcmp(cat->getLabel(), it->second.c_str()) != 0) {
         r = false;
         break;
      }
   }
   return r;
}

//_____________________________________________________________________________
void MultiHistEntry::setParent(RooMultiHistEfficiency* parent)
{
   assert(m_effProd == 0);
   assert(m_relative == 0);

   std::string name;
   for (std::map<RooAbsCategory*, std::string>::const_iterator it = m_rawCats.begin(),
           end = m_rawCats.end(); it != end; ++it) {
      name = it->first->GetName(); name += "_proxy";
      RooCategoryProxy* proxy = new RooCategoryProxy(name.c_str(), name.c_str(), parent,
                                                     *(it->first));
      m_categories.insert(make_pair(proxy, it->second));
   }
   m_rawCats.clear();

   name = m_rawEff->GetName(); name += "_proxy";
   m_effProd = new RooRealProxy(name.c_str(), name.c_str(), parent, *m_rawEff);
   name = m_rawEff->GetName(); name += "_proxy";
   m_relative = new RooRealProxy(name.c_str(), name.c_str(), parent, *m_rawRel);

   // m_rawEff = 0;
   // m_rawRel = 0;
}

//_____________________________________________________________________________
RooArgSet MultiHistEntry::categories() const
{
   RooArgSet r;
   if (!m_rawCats.empty()) {
      for (std::map<RooAbsCategory*, std::string>::const_iterator it = m_rawCats.begin(),
              end = m_rawCats.end(); it != end; ++it) {
         if (!it->first) continue;
         r.add(*(it->first));
      }
   } else {
      for (std::map<RooCategoryProxy*, std::string>::const_iterator it = m_categories.begin(),
              end = m_categories.end(); it != end; ++it) {
         if (!it->first) continue;
         r.add(it->first->arg());
      }
   }

   return r;
}

//_____________________________________________________________________________
void MultiHistEntry::select()
{
   if (!m_rawCats.empty()) {
      for (std::map<RooAbsCategory*, std::string>::const_iterator it = m_rawCats.begin(),
              end = m_rawCats.end(); it != end; ++it) {
         if (!it->first) continue;
         RooAbsCategoryLValue* lval = dynamic_cast<RooAbsCategoryLValue*>(it->first);
         lval->setLabel(it->second.c_str());
      }
   } else {
      for (std::map<RooCategoryProxy*, std::string>::const_iterator it = m_categories.begin(),
              end = m_categories.end(); it != end; ++it) {
         RooAbsCategoryLValue* lval = dynamic_cast<RooAbsCategoryLValue*>(it->first->absArg());
         lval->setLabel(it->second.c_str());
      }
   }
}

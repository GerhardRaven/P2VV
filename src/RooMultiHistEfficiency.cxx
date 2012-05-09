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

#include <RooMultiHistEfficiency.h>

ClassImp(RooMultiHistEfficiency);

// namespace {
//    const char *makeName(const char* name, const RooArgSet& terms ) {
//       TString pname;
//       pname = name;
//       pname.Append("_");
//       RooFIter i = terms.fwdIterator();
//       RooAbsArg *arg;
//       bool first = true;
//       while((arg = i.next())) {
//          if (first) {
//             first= false;
//          } else {
//             pname.Append("_X_");
//          }
//          pname.Append(arg->GetName());
//       }
//       return pname.Data();
//    }

//    RooSuperCategory* makeSuper(const char* name, const RooArgSet& _catVars) {
//       const char *catName = makeName(name, _catVars );
//       return new RooSuperCategory(catName, catName, _catVars);
//    }
// }

//_____________________________________________________________________________
RooMultiHistEfficiency::RooMultiHistEfficiency
(const char *name, const char *title, const RooArgList& efficiencies,
 const RooArgList& categories, const TList& sigCatNames) :
   RooAbsPdf(name, title),
   _categories("categories", "categories", this),
   _efficiencies("efficiencies", "efficiencies", this)
{  
   // Construct an N+1 dimensional efficiency p.d.f from an N-dimensional efficiency
   // function and a category cat with two states (0,1) that indicate if a given
   // event should be counted as rejected or accepted respectively
   RooFIter catIter = categories.fwdIterator();
   while (RooAbsArg* cat = catIter.next()) {
      _categories.add(*cat);
   }

   std::auto_ptr<TIterator> nameIter(sigCatNames.MakeIterator());
   while (TObject* name = nameIter->Next()) {
      _sigCatNames.Add(name);
   }
   _nameIter = _sigCatNames.MakeIterator();
}

//_____________________________________________________________________________
RooMultiHistEfficiency::RooMultiHistEfficiency(const RooMultiHistEfficiency& other, const char* name) : 
   RooAbsPdf(other, name),
   _categories("categories", this, other._categories),
   _efficiencies("efficiencies", this, other._efficiencies)
{
   // Copy constructor
   _sigCatNames.Clear();
   std::auto_ptr<TIterator> nameIter(other._sigCatNames.MakeIterator());
   while (TObject* name = nameIter->Next()) {
      _sigCatNames.Add(name);
   }
   _nameIter = _sigCatNames.MakeIterator();
}

//_____________________________________________________________________________
RooMultiHistEfficiency::~RooMultiHistEfficiency() 
{
   // Destructor
   if (_nameIter) delete _nameIter;
}

// //_____________________________________________________________________________
// RooAbsGenContext* RooMultiHistEfficiency::genContext
// (const RooArgSet &vars, const RooDataSet *prototype,
//  const RooArgSet* auxProto, Bool_t verbose) const
// {
//    // use generic context explicitly allowing generation of effciency observable,
//    // since we deal with it ourselves.
//    return new RooGenContext(*this, vars, prototype, auxProto, verbose, &_observables) ;
// }

// //_____________________________________________________________________________
// void RooMultiHistEfficiency::initGenerator(Int_t code)
// {
//    // Forward one-time initialization call to component generation initialization
//    // methods.
//    pdf()->initGenerator(_pdfGenCode);
//    _levels.clear();
//    _super = makeSuper(GetName(), _categories);
//    std::auto_ptr<TIterator> superIter(_super->MakeIterator());

//    std::auto_ptr<RooAbsReal> marginal(createIntegral(_pdfGenVars));

//    TString current = _super->getLabel();
//    while (TObjString* label = static_cast<TObjString*>(superIter->Next())) {
//       _super->setLabel(label->String());
//       if (allFalse()) continue;

//       double n = marginal->getVal();
//       if (!_levels.empty()) n += _levels.back().first; // cumulative
//       cxcoutD(Generation) << "RooMultiHistEfficiency creating sampler for " << _pdfGenVars
//                           << " given " << _categories << " = "  << label->String() 
//                           << " (level = " << n << ")" << endl;
//       _levels.push_back(make_pair(n, label->String()));
//    }

//    // Given that above we properly marginalized, the next line should be a no-op.
//    for (Levels::iterator i = _levels.begin(); i != _levels.end(); ++i) 
//       i->first /= _levels.back().first; // normalize

//    _super->setLabel(current);
// }

// //_____________________________________________________________________________
// Int_t RooMultiHistEfficiency::getGenerator
// (const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const
// {
//    bool all = true;
//    bool none = true;
//    RooFIter iter = _categories.fwdIterator();
//    while (RooAbsArg* cat = iter.next()) {
//       if (!directVars.find(*cat)) {
//          all &= false;
//       } else {
//          none &= false;
//       }
//    }
//    assert(all ^ none);

//    _pdfGenVars.removeAll();
//    RooArgSet pdfVars(directVars);
//    Int_t pdfCode = pdf()->getGenerator(pdfVars, _pdfGenVars, staticInitOK) ;
//    if (pdfCode != 0) {
//       generateVars.add(_pdfGenVars);
//    }

//    iter = generateVars.fwdIterator();
//    while (RooAbsArg* cat = iter.next()) {
//       if (_categories.find(*cat)) {
//          throw EffHistProd::Exception("Cannot generate categories on which the PDF depends.");
//       }
//    }

//    // RooArgSet* pdfObs = pdf()->getObservables(_pdfGenVars);

//    if (generateVars.getSize() == 0) {
//       return 0;
//       _pdfGenCode = 0;
//    } else if (none) {
//       _pdfGenCode = pdfCode;
//       return 1;
//    } else {
//       iter = _categories.fwdIterator();
//       while (RooAbsArg* cat = iter.next()) {
//          generateVars.add(*cat);
//       }
//       _pdfGenCode = pdfCode;
//       return 2;
//    }
// }

// //_____________________________________________________________________________
// void RooMultiHistEfficiency::generateEvent(Int_t code)
// {
//    // generate values for the variables corresponding to the generation code
//    assert(code > 0);

//    // generate categories
//    Double_t r = RooRandom::uniform();
//    Levels::const_iterator it = _levels.begin(), end = _levels.end();      
//    while (it != end && it->first < r) ++it;
//    assert(it != end);
//    _super->setLabel(it->second); // this should assign _catVars...

//    while (true) {
//       // use the pdf to generate the observables.
//       pdf()->generateEvent(_pdfGenCode);
//       double val = effVal();
//       if (val > 1.) {
//          coutE(Generation) << ClassName() << "::" << GetName() 
//               << ":generateEvent: value of efficiency is larger than assumed maximum of 1."  << endl;
//          continue;
//       }
//       if (val > RooRandom::uniform()) break;
//    }
   
// }

//_____________________________________________________________________________
Double_t RooMultiHistEfficiency::evaluate() const
{
   // Calculate the raw value of this p.d.f which is the effFunc
   // value if cat==1 and it is (1-effFunc) if cat==0
   Double_t val = 0.;
   bool first = true;

   RooFIter catIter = _categories.fwdIterator();
   RooFIter effIter = _efficiencies.fwdIterator();

   _nameIter->Reset();

   while (RooAbsReal* eff = static_cast<RooAbsReal*>(effIter.next())) {
      RooAbsCategory* cat = static_cast<RooAbsCategory*>(catIter.next());
      TObjString* name = static_cast<TObjString*>(_nameIter->Next());
      double eps = 0.;
      if (0 == strcmp(name->GetString().Data(), cat->getLabel())) {
         // Accept case
         eps = eff->getVal();
      } else {
         // Reject case
         eps = 1 - eff->getVal();
      }

      cout << eff->GetName() << " = " << " val = " << eps << endl;

      // Truncate efficiency function in range 0.0-1.0
      if (eps > 1) {
         eps = 1.0;
      } else if (eps < 0) {
         eps = 0.0;
      }

      if (first) {
         val = eps;
         first = false;
      } else {
         val *= eps;
      }
   }
   return val;
}

//_____________________________________________________________________________
bool RooMultiHistEfficiency::allFalse() const
{
   RooFIter catIter = _categories.fwdIterator();
   _nameIter->Reset();
   while (RooAbsCategory* cat = static_cast<RooAbsCategory*>(catIter.next())) {
      TObjString* name = static_cast<TObjString*>(_nameIter->Next());
      if (0 == strcmp(name->GetString().Data(), cat->getLabel())) {
         // Accept
         return false;
      }
   }
   return true;
}

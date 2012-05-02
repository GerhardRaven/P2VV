//////////////////////////////////////////////////////////////////////////////
// 
// BEGIN_HTML
// RooMultiEfficiency is a PDF helper class to fit efficiencies parameterized
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

#include <RooFit.h>

#include <RooMultiEfficiency.h>
#include <RooStreamParser.h>
#include <RooArgList.h>
#include <RooCategory.h>

ClassImp(RooMultiEfficiency);

//_____________________________________________________________________________
RooMultiEfficiency::RooMultiEfficiency
(const char *name, const char *title, const RooArgList& efficiencies, const RooArgList& categories,
 const TList& sigCatNames) :
   RooAbsPdf(name,title),
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

   RooFIter effIter = efficiencies.fwdIterator();
   while (RooAbsArg* eff = effIter.next()) {
      _efficiencies.add(*eff);
   }

   TIterator *nameIter = sigCatNames.MakeIterator();
   while (TObject* name = nameIter->Next()) {
      _sigCatNames.Add(name);
   }
   delete nameIter;
}

//_____________________________________________________________________________
RooMultiEfficiency::RooMultiEfficiency(const RooMultiEfficiency& other, const char* name) : 
   RooAbsPdf(other, name),
   _categories("categories", this, other._categories),
   _efficiencies("efficiencies", this, other._efficiencies)
{
   // Copy constructor
   _sigCatNames.Clear();
   TIterator *nameIter = other._sigCatNames.MakeIterator();
   while (TObject* name = nameIter->Next()) {
      _sigCatNames.Add(name);
   }
   delete nameIter;
}

//_____________________________________________________________________________
RooMultiEfficiency::~RooMultiEfficiency() 
{
   // Destructor
}

//_____________________________________________________________________________
Double_t RooMultiEfficiency::evaluate() const
{
   // Calculate the raw value of this p.d.f which is the effFunc
   // value if cat==1 and it is (1-effFunc) if cat==0
   Double_t val = 0.;
   bool first = true;

   RooFIter catIter = _categories.fwdIterator();
   RooFIter effIter = _efficiencies.fwdIterator();
   TIterator* nameIter = _sigCatNames.MakeIterator();

   while (RooAbsReal* eff = static_cast<RooAbsReal*>(effIter.next())) {
      RooAbsCategory* cat = static_cast<RooAbsCategory*>(catIter.next());
      TObjString* name = static_cast<TObjString*>(nameIter->Next());

      double eps = 0.;
      if (0 == strcmp(name->GetString().Data(), cat->getLabel())) {
         // Accept case
         eps = eff->getVal();
      } else {
         // Reject case
         eps = 1 - eff->getVal();
      }

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
   delete nameIter;
   return val;
}

//_____________________________________________________________________________
Int_t RooMultiEfficiency::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const 
{
   if (matchArgs(allVars,analVars,_categories)) return 1 ;
   return 0 ;
}

//_____________________________________________________________________________
Double_t RooMultiEfficiency::analyticalIntegral(Int_t code, const char* /*rangeName*/) const 
{
   assert(code==1) ;
   return 1.0 ;
}

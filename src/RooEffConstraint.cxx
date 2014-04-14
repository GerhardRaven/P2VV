/***************************************************************************** 
 * Projects: RooFit, P2VV                                                    * 
 *                                                                           * 
 * Author: Roel Aaij, roel.aaij@nikhef.nl                                    * 
 *****************************************************************************/ 

// Your description goes here... 
#include <vector>
#include <iostream>
#include <cmath> 
#include <memory>

#include "TMath.h" 

#include "Riostream.h" 
#include "RooAbsArg.h"
#include "RooStringVar.h"

#include "P2VV/RooEffConstraint.h" 

namespace {
   using std::vector;
   using std::auto_ptr;
}

//_____________________________________________________________________________
RooEffConstraint::RooEffConstraint()
   : RooAbsPdf(),
     _eps_a("!eps_a", "eps a", this),
     _eps_b("!eps_b", "eps b", this)
{
}

//_____________________________________________________________________________
RooEffConstraint::RooEffConstraint(const char *name, const char *title, 
                                   const RooArgList& eps_a,
                                   const RooArgList& eps_b, const std::vector<double>& N_a,
                                   const std::vector<double>& N_b,
                                   const std::vector<double>& N_ab)
   : RooAbsPdf(name,title),
     _eps_a("!eps_a", "eps a", this),
     _eps_b("!eps_b", "eps b", this),
     _f_a(N_a.size(), 0.), _f_b(N_b.size(), 0.), _f_ab(N_ab.size(), 0.),
     _N_a(N_a), _N_b(N_b), _N_ab(N_ab)
{ 
   assert((eps_a.getSize() > 1 && eps_b.getSize() > 1 && eps_a.getSize() == eps_b.getSize())
          || (eps_a.getSize() == 1 && eps_b.getSize() > 1)
          || (eps_b.getSize() == 1 && eps_a.getSize() > 1));

   if (eps_a.getSize() == 1) {
      _eps_a.add(eps_b);
      _eps_b.add(eps_a);
   } else {
      _eps_a.add(eps_a);
      _eps_b.add(eps_b);
   }
}

//_____________________________________________________________________________
RooEffConstraint::RooEffConstraint(const RooEffConstraint& other, const char* name)
   : RooAbsPdf(other,name), 
     _eps_a("!eps_a", this, other._eps_a),
     _eps_b("!eps_b", this, other._eps_b),
     _f_a(other._f_a), _f_b(other._f_b), _f_ab(other._f_ab),
     _N_a(other._N_a), _N_b(other._N_b), _N_ab(other._N_ab)
{
} 

//_____________________________________________________________________________
RooEffConstraint::~RooEffConstraint()
{
}

//_____________________________________________________________________________
Bool_t RooEffConstraint::forceAnalyticalInt(const RooAbsArg& /*dep*/) const
{
   // Return kTRUE to force RooRealIntegral to offer all observables for internal integration
   return true;
}

//_____________________________________________________________________________
Int_t RooEffConstraint::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                                                const char* /*rangeName*/) const 
{
   analVars.add(allVars);
   return 1;
}

//_____________________________________________________________________________
Double_t RooEffConstraint::analyticalIntegral(Int_t code, const char* /*rangeName*/) const
{
   if (code != 1) {
      coutF(InputArguments) << "RooEffConstraint::analyticalIntegral("
			    << GetName() << "): integration code should be 1 (got" << code << ")"
			    << std::endl;
     assert(0);
   }

   return 1.;
}

//_____________________________________________________________________________
Double_t RooEffConstraint::getLogVal(const RooArgSet* nset) const
{
   double s = 0;
   for (int i = 0; i < _eps_a.getSize(); ++i) {
      const RooAbsReal* ea = dynamic_cast<const RooAbsReal*>(_eps_a.at(i));
      const RooAbsReal* eb = 0;
      if (_eps_b.getSize() == 1) {
         eb = dynamic_cast<const RooAbsReal*>(_eps_b.at(0));
      } else {
         eb = dynamic_cast<const RooAbsReal*>(_eps_b.at(i));
      }

      if (_N_ab.empty()) {
         // Like HLT1
         double r = eb->getVal() / ea->getVal();
         double v = _N_a[i] * log(fraction(1 / (1 + r), _f_a, i)) 
            + _N_b[i] * log(fraction(r / (1 + r), _f_b, i));
         s += v;
      } else {
         // Like HLT2
         double eav = ea->getVal();
         double ebv = eb->getVal();
         double num = eav + ebv - eav * ebv;
         
         double v = _N_a[i] * log(fraction((1 - ebv) * eav / num, _f_a, i))
            + _N_b[i] * log(fraction((1 - eav) * ebv / num, _f_b, i))
            + _N_ab[i] * log(fraction(eav * ebv / num, _f_ab, i));
         s += v;
      }
   }
   return s;
}

//_____________________________________________________________________________
Double_t RooEffConstraint::evaluate() const 
{ 
   return exp(getLogVal());
} 

//_____________________________________________________________________________
Double_t RooEffConstraint::fraction(const double f, std::vector<double>& v, const int i) const
{
   if (v[i] == 0.) {
      v[i] = f;
      return 1.;
   } else {
      return f / v[i];
   }
}

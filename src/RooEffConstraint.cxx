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
     _N_a(N_a), _N_b(N_b), _N_ab(N_ab)
{ 
   _eps_a.add(eps_a);
   _eps_b.add(eps_b);
}

//_____________________________________________________________________________
RooEffConstraint::RooEffConstraint(const RooEffConstraint& other, const char* name)
   : RooAbsPdf(other,name), 
     _eps_a("!eps_a", this, other._eps_a),
     _eps_b("!eps_b", this, other._eps_b),
     _N_a(other._N_a),
     _N_b(other._N_b),
     _N_ab(other._N_ab)
{
} 

//_____________________________________________________________________________
RooEffConstraint::~RooEffConstraint()
{
}

// //_____________________________________________________________________________
// Bool_t RooEffConstraint::forceAnalyticalInt(const RooAbsArg& /*dep*/) const
// {
//    // Return kTRUE to force RooRealIntegral to offer all observables for internal integration
//    return true;
// }

// //_____________________________________________________________________________
// Int_t RooEffConstraint::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
//                                                 const char* /*rangeName*/) const 
// {
//    analVars.add(allVars);
//    return 1;
// }

// //_____________________________________________________________________________
// Double_t RooEffConstraint::analyticalIntegral(Int_t code, const char* /*rangeName*/) const
// {
//    if (code != 1) {
//       coutF(InputArguments) << "RooEffConstraint::analyticalIntegral("
// 			    << GetName() << "): integration code should be 1 (got" << code << ")"
// 			    << std::endl;
//      assert(0);
//    }

//    return 1.;
// }

//_____________________________________________________________________________
Double_t RooEffConstraint::getLogVal(const RooArgSet* nset) const
{
   double s = 0;
   for (int i = 0; i < _eps_a.getSize(); ++i) {
      const RooAbsReal* ea = dynamic_cast<const RooAbsReal*>(_eps_a.at(i));
      const RooAbsReal* eb = dynamic_cast<const RooAbsReal*>(_eps_b.at(i));

      if (_N_ab.empty()) {
         // Like HLT1
         double r = eb->getVal() / ea->getVal();
         s += _N_a[i] * log(1 / (1 + r)) + _N_b[i] * log(r / (1 + r));
      } else {
         // Like HLT2
         double eav = ea->getVal();
         double ebv = eb->getVal();
         double num = eav + ebv - eav * ebv;
         s += _N_a[i] * log((1 - ebv) * eav / num) 
            + _N_b[i] * log((1 - eav) * ebv / num)
            + _N_ab[i] * log(eav * ebv / num);
      }
   }
   return s;
}

//_____________________________________________________________________________
Double_t RooEffConstraint::evaluate() const 
{ 
   return exp(getLogVal());
} 

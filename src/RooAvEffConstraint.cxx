/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 
#include <vector>
#include <iostream>
#include <cmath> 
#include <memory>

#include "TMath.h" 

#include "Riostream.h" 
#include "RooAbsPdf.h" 
#include "RooAbsArg.h"
#include "RooCustomizer.h"
#include "RooStringVar.h"

#include "RooAvEffConstraint.h" 

ClassImp(RooAvEffConstraint) 

using std::vector;
using std::auto_ptr;

//_____________________________________________________________________________
RooAvEffConstraint::RooAvEffConstraint(const char *name, const char *title, 
                                       RooAbsPdf& pdf, RooEffResModel& model,
                                       RooRealVar& mean, RooRealVar& sigma)
   : RooAbsPdf(name,title),
     _mean("!average_mean", "average_mean proxy", this, mean),
     _sigma("!average_sigma", "average_sigma proxy", this, sigma),
     _integral(0),
     _integrals("!integrals", 0, this),
     _efficiencies("!efficiencies", 0, this)
{ 
   auto_ptr<RooArgSet> vars(pdf.getVariables());
   RooFIter it = vars->fwdIterator();
   RooAbsArg* arg = 0;
   while ((arg = it.next())) {
      if (arg->getAttribute("Observable")) {
         RooAbsRealLValue* l = dynamic_cast<RooAbsRealLValue*>(arg);
         if (l) l->setConstant(true);
      }
   }
   // observable.setConstant(true);

   RooArgSet iset(model.convVar());
   RooAbsReal* I = pdf.createIntegral(iset);
   TString intName = TString(model.efficiency()->GetName()) + "_average_" + I->GetName();
   I->SetName(intName.Data());
   _integral = new RooRealProxy("!average_integral", "average_integral",
                                const_cast<RooAvEffConstraint*>(this), *I, false, true);

   RooRealVar& x = model.convVar(); // binboundaries not const...

   const RooArgList& ranges = model.getIntegralRanges(iset);
   it = ranges.fwdIterator();
   while (const RooStringVar* rangeName = static_cast<const RooStringVar*>(it.next())) {
      const char* range = rangeName->getVal();
      I = pdf.createIntegral(iset, range);
      _integrals.addOwned(*I);

      Double_t xmin = x.getMin(range);
      Double_t xmax = x.getMax(range);

      RooCustomizer* customizer = new RooCustomizer(*model.efficiency(),
                                                    (TString(range) + "_customizer").Data());
      RooRealVar* cv = static_cast<RooRealVar*>(x.clone(TString(x.GetName()) + "_" + range) );
      cv->setVal((xmin + xmax) / 2.);
      cv->setConstant(true);
      customizer->replaceArg(x, *cv);
      RooAbsArg *ceff = customizer->build(kFALSE);
      ceff->addOwnedComponents(*cv);
      _efficiencies.addOwned(*ceff);
      _customizers.push_back(customizer);
   }
   
} 

//_____________________________________________________________________________
RooAvEffConstraint::RooAvEffConstraint(const RooAvEffConstraint& other, const char* name)
   : RooAbsPdf(other,name), 
     _mean("!average_mean", this, other._mean),
     _sigma("!average_sigma", this, other._sigma),
     _integrals("!integrals", this, other._integrals),
     _efficiencies("!efficiencies", this, other._efficiencies)

{
   if (other._integral) {
      _integral = new RooRealProxy("!average_integral", this, *other._integral);
   } else {
      _integral = 0;
   }
} 

//_____________________________________________________________________________
RooAvEffConstraint::~RooAvEffConstraint()
{
   if (_integral) delete _integral;
   for (std::vector<RooCustomizer*>::const_iterator it = _customizers.begin(),
           end = _customizers.end(); it != end; ++it) {
      delete *it;
   }
}

//_____________________________________________________________________________
Bool_t RooAvEffConstraint::forceAnalyticalInt(const RooAbsArg& /*dep*/) const
{
   // Return kTRUE to force RooRealIntegral to offer all observables for internal integration
   return true;
}

//_____________________________________________________________________________
Int_t RooAvEffConstraint::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                                                const char* /*rangeName*/) const 
{
   analVars.add(allVars);
   return 1;
}

//_____________________________________________________________________________
Double_t RooAvEffConstraint::analyticalIntegral(Int_t code, const char* rangeName) const 
{
   assert(code == 1);
   return 1.;
}

//_____________________________________________________________________________
Double_t RooAvEffConstraint::evaluate() const 
{ 
   double av = 0;

   for (int i = 0; i < _integrals.getSize(); ++i) {
      const RooAbsReal* entry = dynamic_cast<const RooAbsReal*>(_integrals.at(i));
      assert(entry);
      const RooAbsReal* efficiency = dynamic_cast<const RooAbsReal*>(_efficiencies.at(i));
      av += efficiency->getVal() * entry->getVal();
   }

   double integral = static_cast<const RooAbsReal&>(_integral->arg()).getVal();
   av /= integral;

   Double_t arg = av - _mean;
   Double_t sig = _sigma;
   return exp(-0.5 * arg * arg / (sig * sig));
} 

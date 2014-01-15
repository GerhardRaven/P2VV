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
#include "RooAbsPdf.h" 
#include "RooAbsArg.h"
#include "RooCustomizer.h"
#include "RooStringVar.h"

#include "P2VV/RooAvEffConstraint.h" 

using std::vector;
using std::auto_ptr;

//_____________________________________________________________________________
RooAvEffConstraint::RooAvEffConstraint()
   : RooAbsPdf(), _integral(0), _pdf(0)
{
   _integrals = new RooListProxy("!integrals", "!integrals", this);
   _efficiencies = new RooListProxy("!efficiencies", "!integrals", this);
}

//_____________________________________________________________________________
RooAvEffConstraint::RooAvEffConstraint(const char *name, const char *title, 
                                       RooAbsPdf& pdf, RooEffResModel& model,
                                       RooRealVar& mean, RooRealVar& sigma)
   : RooAbsPdf(name,title),
     _mean("!average_mean", "average_mean proxy", this, mean),
     _sigma("!average_sigma", "average_sigma proxy", this, sigma),
     _integral(0),
     _model("!res_model", "!res_model", this, model),
     _pdf(0)
{ 
   _pdf = new RooRealProxy("!pdf", "!pdf", this, pdf);
   _integrals = new RooListProxy("!integrals", "!integrals", this);
   _efficiencies = new RooListProxy("!efficiencies", "!integrals", this);
} 

//_____________________________________________________________________________
RooAvEffConstraint::RooAvEffConstraint(const RooAvEffConstraint& other, const char* name)
   : RooAbsPdf(other,name), 
     _mean("!average_mean", this, other._mean),
     _sigma("!average_sigma", this, other._sigma),
     _integral(0),
     _integrals(0),
     _efficiencies(0),
     _model("!res_model", this, other._model)
{
   if (other._pdf) {
      _pdf = new RooRealProxy("!pdf", this, *other._pdf);
   } else {
      _pdf = 0;
   }

   if (other._integrals) {
      _integrals = new RooListProxy("!integrals", this, *other._integrals);
   } else {
      _integrals = 0;
   }

   if (other._efficiencies) {
      _efficiencies = new RooListProxy("!efficiencies", this, *other._efficiencies);
   } else {
      _efficiencies = 0;
   }
} 

//_____________________________________________________________________________
RooAvEffConstraint::~RooAvEffConstraint()
{
   delete _integral;
   delete _integrals;
   delete _efficiencies;
   delete _pdf;
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
Double_t RooAvEffConstraint::analyticalIntegral(Int_t code, const char* /*rangeName*/) const
{
   if (code != 1) {
      coutF(InputArguments) << "RooAvEffConstraint::analyticalIntegral("
			    << GetName() << "): integration code should be 1 (got" << code << ")"
			    << std::endl;
     assert(0);
   }

   return 1.;
}

//_____________________________________________________________________________
const RooAbsReal* RooAvEffConstraint::efficiency() const
{
   const RooEffResModel& res_model = dynamic_cast<const RooEffResModel&>(_model.arg());
   return res_model.efficiency();
}

//_____________________________________________________________________________
void RooAvEffConstraint::setPdf(RooAbsPdf& pdf)
{
   if (_pdf) {
      delete _pdf;
      _pdf = 0;
   }
   _pdf = new RooRealProxy("!pdf", "!pdf", this, pdf);
   delete _integral;
   _integral = 0;
   if(_integrals) {
      _integrals->removeAll();
      delete _integrals;
   }
   _integrals = new RooListProxy("!integrals", 0, this);

   if (_efficiencies) {
      _efficiencies->removeAll();
      delete _efficiencies;
   }
   _efficiencies = new RooListProxy("!efficiencies", 0, this);
}

//_____________________________________________________________________________
Double_t RooAvEffConstraint::evaluate() const 
{ 
   if (!_integrals || !_integrals->getSize()) {
      const_cast<RooAvEffConstraint*>(this)->initialize();
   }

   double av = 0;

   for (int i = 0; i < _integrals->getSize(); ++i) {
      const RooAbsReal* entry = dynamic_cast<const RooAbsReal*>(_integrals->at(i));
      assert(entry);
      const RooAbsReal* efficiency = dynamic_cast<const RooAbsReal*>(_efficiencies->at(i));
      double e = efficiency->getVal();
      double v = entry->getVal();
      av += e * v;
   }

   double integral = static_cast<const RooAbsReal&>(_integral->arg()).getVal();
   av /= integral;

   Double_t arg = av - _mean;
   Double_t sig = _sigma;
   return exp(-0.5 * arg * arg / (sig * sig));
} 

//_____________________________________________________________________________
void RooAvEffConstraint::initialize()
{
   const RooAbsPdf& pdf = static_cast<const RooAbsPdf&>(_pdf->arg());

   auto_ptr<RooArgSet> vars(pdf.getVariables());
   RooFIter it = vars->fwdIterator();
   RooAbsArg* arg = 0;
   while ((arg = it.next())) {
      if (arg->getAttribute("Observable")) {
         RooAbsRealLValue* l = dynamic_cast<RooAbsRealLValue*>(arg);
         if (l) l->setConstant(true);
      }
   }

   const RooEffResModel& model = static_cast<const RooEffResModel&>(_model.arg());

   RooArgSet* iset = new RooArgSet(model.convVar());
   RooAbsReal* I = pdf.createIntegral(*iset);
   TString intName = TString(model.efficiency()->GetName()) + "_average_" + I->GetName();
   I->SetName(intName.Data());

   _integral = new RooRealProxy("!average_integral", "average_integral",
                                const_cast<RooAvEffConstraint*>(this), *I, false, true);

   RooRealVar& x = model.convVar(); // binboundaries not const...

   const RooArgList& ranges = model.getIntegralRanges(*iset);
   it = ranges.fwdIterator();
   while (const RooStringVar* rangeName = static_cast<const RooStringVar*>(it.next())) {
      const char* range = rangeName->getVal();
      I = pdf.createIntegral(*iset, range);
      _integrals->addOwned(*I);

      Double_t xmin = x.getMin(range);
      Double_t xmax = x.getMax(range);

      RooCustomizer customizer(*model.efficiency(), (TString(range) + "_customizer").Data());
      RooRealVar* cv = static_cast<RooRealVar*>(x.clone(TString(x.GetName()) + "_" + range) );
      cv->setVal((xmin + xmax) / 2.);
      cv->setConstant(true);
      customizer.replaceArg(x, *cv);
      RooAbsArg *ceff = customizer.build(kFALSE);
      ceff->addOwnedComponents(*cv);
      _efficiencies->addOwned(*ceff);
   }   
}

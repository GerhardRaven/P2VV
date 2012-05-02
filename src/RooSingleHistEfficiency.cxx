/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id: RooSingleHistEfficiency.cxx 25184 2008-08-20 13:59:55Z wouter $
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
// The class RooSingleHistEfficiency implements the product of a PDF with an efficiency function.
// The normalization integral of the product is calculated numerically, but the
// event generation is handled by a specialized generator context that implements
// the event generation in a more efficient for cases where the PDF has an internal
// generator that is smarter than accept reject. 
// END_HTML
//

#include "RooFit.h"
#include "RooSingleHistEfficiency.h"
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

ClassImp(RooSingleHistEfficiency);

//_____________________________________________________________________________
RooSingleHistEfficiency::RooSingleHistEfficiency(const char *name, const char *title, 
                               RooAbsPdf& inPdf, RooAbsReal& inEff)
   : RooEffHistProd(name, title, inPdf),
     _eff("eff", "efficiency function", this, inEff),
     _observables("obs", "observables in efficiency function", this)
{  
   // Constructor of a a production of p.d.f inPdf with efficiency
   // function inEff.

   // to figure out what the observable is, we look at the overlap of
   // the variables of efficiency function and pdf.

   // FIXME: this assumes that the only overlap in variables between the PDF and the
   // efficiency are observables. If a parameter overlaps, this might fail.
  
   std::auto_ptr<RooArgSet> pdfpars(inPdf.getVariables());
   std::auto_ptr<RooArgSet> effpars(inEff.getVariables());
   std::auto_ptr<TIterator> iter(effpars->createIterator());

   RooAbsArg *effelem(0);
   while((effelem = (RooAbsArg*)iter->Next())) {
      if(pdfpars->find(effelem->GetName()))  _observables.add(*effelem);
   }

   if(_observables.getSize() == 0) {
      throw std::string("WARNING: RooSingleHistEfficiency: PDF and Efficiency function "
                        "factorise. Please use RooProd");
   } else if(_observables.getSize() > 1) {
      throw std::string("WARNING: RooSingleHistEfficiency not yet implemented for more than 1D efficiency" );
   }
  
   // an interesting hack. need to discuss with wouter. one idea: let
   // every function that is 'binned' (discrete,quantized,..) add a 
   // special binning object to its dependents.
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,32,0)
   RooRealVar* x = static_cast<RooRealVar*>(_observables.first());
   std::list<Double_t>* binboundaries = inEff.binBoundaries(*x, x->getMin(), x->getMax());
   std::copy(binboundaries->begin(), binboundaries->end(), std::back_inserter(_binboundaries));
#else
   std::list<Double_t>* binboundaries = inEff.plotSamplingHint(*x, x->getMin(), x->getMax());
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
RooSingleHistEfficiency::RooSingleHistEfficiency(const RooSingleHistEfficiency& other, const char* name):
   RooEffHistProd(other, name),
   _eff("acc", this, other._eff),
   _observables("obs", this, other._observables)
{
   // Copy constructor
}

//_____________________________________________________________________________
RooSingleHistEfficiency::~RooSingleHistEfficiency() 
{
   // Destructor
}

//_____________________________________________________________________________
RooAbsGenContext* RooSingleHistEfficiency::genContext
(const RooArgSet &vars, const RooDataSet *prototype, const RooArgSet* auxProto, Bool_t verbose) const
{
   // Return specialized generator context for RooSingleHistEfficiencys that implements generation
   // in a more efficient way than can be done for generic correlated products
   assert(pdf() != 0);
   assert(eff() != 0);
   return new RooEffGenContext(*this, *pdf(), *eff(), vars, prototype, auxProto, verbose);
}

//_____________________________________________________________________________
const RooArgSet* RooSingleHistEfficiency::effObservables() const
{
   return static_cast<const RooArgSet*>(&_observables);
}

//_____________________________________________________________________________
double RooSingleHistEfficiency::effVal() const
{
   return eff()->getVal();
}

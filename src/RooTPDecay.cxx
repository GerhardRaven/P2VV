/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Single sided decay function that can be analytically convolved
// with any RooResolutionModel implementation
// 
// Special 
// END_HTML
//
// STL
#include <algorithm>

// RooFit
#include "RooFit.h"
#include "Riostream.h"
#include "Riostream.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooGenContext.h"
#include "RooResolutionModel.h"
#include "RooParamBinning.h"
#include "RooConvGenContext.h"

// Local
#include "P2VV/RooTPConvGenContext.h"
#include "P2VV/RooTPDecay.h"

using namespace std;

//_____________________________________________________________________________
RooTPDecay::RooTPDecay(const char *name, const char *title, 
                       RooRealVar& t, RooAbsReal& tau, const RooResolutionModel& model,
                       DecayType type, const char* normRange)
   : RooDecay(name, title, t, tau, model, type)
{
   setNormRange(normRange);

   TObjArray* oa = TString(normRange).Tokenize(",");
   for (int i = 0; i < oa->GetEntries(); ++i) {
      TObjString* os = (TObjString*)(*oa)[i];
      if(!os) break;
      try {
         const RooParamBinning& binning = dynamic_cast<const RooParamBinning&>(t.getBinning(os->String(), true, false));
         const RooAbsRealLValue* bf = dynamic_cast<const RooAbsRealLValue*>(binning.lowBoundFunc());
         assert(bf && bf->getMin() == t.getMin());
         _tps.add(*bf);
         bf = dynamic_cast<const RooAbsRealLValue*>(binning.highBoundFunc());
         assert(bf && bf->getMax() == t.getMax());
         _tps.add(*bf);
      } catch (const std::bad_cast&) {
         delete oa;
         assert(false);
      }
   }
   delete oa;
}


//_____________________________________________________________________________
RooTPDecay::RooTPDecay(const RooTPDecay& other, const char* name)
   : RooDecay(other, name),
     _tps("!tps", this, other._tps)
{
}

//_____________________________________________________________________________
RooTPDecay::~RooTPDecay()
{
}

//_____________________________________________________________________________
RooAbsGenContext* RooTPDecay::genContext(const RooArgSet &vars, const RooDataSet *prototype,
                                         const RooArgSet* auxProto, Bool_t verbose) const
{
   // Create a generator context for this p.d.f. If both the p.d.f and the resolution model
   // support internal generation of the convolution observable on an infinite domain,
   // deploy a specialized convolution generator context, which generates the physics distribution
   // and the smearing separately, adding them a posteriori. If this is not possible return
   // a (slower) generic generation context that uses accept/reject sampling

   // Check if the resolution model specifies a special context to be used.
   RooResolutionModel* conv = dynamic_cast<RooResolutionModel*>(_model.absArg());
   assert(conv);

   RooArgSet* modelDep = _model.absArg()->getObservables(&vars) ;
   modelDep->remove(*convVar(),kTRUE,kTRUE) ;
   Int_t numAddDep = modelDep->getSize() ;
   delete modelDep ;

   // Check if physics PDF and resolution model can both directly generate the convolution variable
   RooArgSet dummy ;
   Bool_t pdfCanDir = (getGenerator(*convVar(),dummy) != 0) ;
   Bool_t resCanDir = conv && (conv->getGenerator(*convVar(),dummy)!=0) && conv->isDirectGenSafe(*convVar()) ;

   if (numAddDep>0 || !pdfCanDir || !resCanDir) {
      // Any resolution model with more dependents than the convolution variable
      // or pdf or resmodel do not support direct generation
      string reason ;
      if (numAddDep>0) reason += "Resolution model has more onservables that the convolution variable. " ;
      if (!pdfCanDir) reason += "PDF does not support internal generation of convolution observable. " ;
      if (!resCanDir) reason += "Resolution model does not support internal generation of convolution observable. " ;

      coutI(Generation) << "RooTPDecay::genContext(" << GetName() << ") Using regular accept/reject generator " 
                        << "for convolution p.d.f because: " << reason.c_str() << endl ;    
      return new RooGenContext(*this,vars,prototype,auxProto,verbose) ;
   } 

   RooAbsGenContext* context = conv->modelGenContext(*this, vars, prototype, auxProto, verbose);
   if (context) return context;
  
   // Any other resolution model: use specialized generator context
   // return new RooTPConvGenContext(*this, vars, _tps, prototype, auxProto, verbose);
   return new RooConvGenContext(*this, vars, prototype, auxProto, verbose);
}

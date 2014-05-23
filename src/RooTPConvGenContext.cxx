/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id: RooTPConvGenContext.cxx 28259 2009-04-16 16:21:16Z wouter $
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
// RooTPConvGenContext is an efficient implementation of the generator context
// specific for RooAbsAnaConvPdf objects. The physics model is generated
// with a truth resolution model and the requested resolution model is generated
// separately as a PDF. The convolution variable of the physics model is
// subsequently explicitly smeared with the resolution model distribution.
// END_HTML
//
#include <iostream>

#include "RooFit.h"

#include "RooMsgService.h"
#include "RooAbsAnaConvPdf.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "Riostream.h"
#include "RooRealVar.h"

#include "P2VV/RooTPConvGenContext.h"

namespace {
   using std::cout;
   using std::endl;
}

//_____________________________________________________________________________
RooTPConvGenContext::RooTPConvGenContext(const RooAbsAnaConvPdf &model, const RooArgSet &vars, const RooArgList& tps,
                                         const RooDataSet *prototype, const RooArgSet* auxProto, Bool_t verbose)
   : RooConvGenContext(model, vars, prototype, auxProto, verbose)
{
   cxcoutI(Generation) << "RooTPConvGenContext::ctor() setting up special generator context "
                       << "to apply an efficiency to the analytical convolution p.d.f. "
                       << model.GetName() << " for generation of observable(s) " << vars << endl ;

   const RooArgSet* data_vars = prototype->get();
   for (int i = 0; i < tps.getSize(); ++i) {
      RooAbsArg* tp = tps.at(i);
      RooAbsArg* tp_data = data_vars->find(tp->GetName());
      if (!tp_data) {
         cout << "tp " << tp->GetName() << " is not in proto data." << endl;
         assert(tp_data);
      }
      _tps.add(*tp_data);
   }
   _normRange = model.normRange();
}

//_____________________________________________________________________________
void RooTPConvGenContext::generateEvent(RooArgSet &theEvent, Int_t remaining)
{
   // Generate a single event

   while(1) {
      // Generate pdf and model data
      _modelGen->generateEvent(*_modelVars, remaining) ;
      _pdfGen->generateEvent(*_pdfVars, remaining) ;

      // Construct smeared convolution variable
      Double_t tval = _cvPdf->getVal() + _cvModel->getVal() ;
      // Check if it is in range.
      for (int i = 0; i < _tps.getSize(); i += 2) {
         const RooAbsReal* tp_min = static_cast<const RooAbsReal*>(_tps.at(i));
         const RooAbsReal* tp_max = static_cast<const RooAbsReal*>(_tps.at(i + 1));
         if (tval > tp_min->getVal() && tval < tp_max->getVal()) {
            theEvent = *_modelVars ;
            theEvent = *_pdfVars ;
            _cvOut->setVal(tval) ;
            theEvent.Print("v");
            return;
         }
      }
   }
}

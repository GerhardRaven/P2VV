/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   Gerhard Raven
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_CUBICBSPLINE
#define ROO_CUBICBSPLINE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooRealVar;
class RooArgList ;

class RooCubicBSpline : public RooAbsPdf {
public:

  RooCubicBSpline() ;
  RooCubicBSpline(const char *name, const char *title,
               RooRealVar& _x, const char *knotBinningName, const RooArgList& _coefList) ;

  RooCubicBSpline(const RooCubicBSpline& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooCubicBSpline(*this, newname); }
  inline virtual ~RooCubicBSpline() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName) const;

private:

  RooRealProxy _x;
  RooListProxy _coefList ;

  class _auxKnot;
  const _auxKnot *_aux; // do not persist! (but do persist the binningName used for _x!!!

  Double_t evaluate() const;

  ClassDef(RooCubicBSpline,1) // Bernstein polynomial PDF
};

#endif

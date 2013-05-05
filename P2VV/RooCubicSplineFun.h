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
#ifndef ROO_CUBICSPLINEFUN
#define ROO_CUBICSPLINEFUN

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooRealVar;
class RooArgList ;
class RooCubicSplineKnot;

class RooCubicSplineFun : public RooAbsReal {
public:

  RooCubicSplineFun() ;
  RooCubicSplineFun(const char *name, const char *title,
               RooRealVar& _x, const char *knotBinningName, const RooArgList& _coefList) ;
  ~RooCubicSplineFun() ;

  RooCubicSplineFun(const RooCubicSplineFun& other, const char* name = 0);
  TObject* clone(const char* newname) const { return new RooCubicSplineFun(*this, newname); }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName) const;

private:

  RooRealProxy _x;
  RooListProxy _coefList ;

  const RooCubicSplineKnot *_aux; // do not persist! (but do persist the binningName used for _x!!!

  Double_t evaluate() const;

  ClassDef(RooCubicSplineFun,1) // CubicSpline polynomial PDF
};

#endif

/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooCubicSplineGaussModel.h,v 1.21 2007/05/11 09:13:07 verkerke Exp $
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
#ifndef ROO_CS_GAUSS_MODEL
#define ROO_CS_GAUSS_MODEL

#include <complex>
#include <RooResolutionModel.h>
#include <RooRealProxy.h>

#include <P2VV/RooAbsEffResModel.h>
#include "P2VV/RooAbsGaussModelEfficiency.h"

class RooAbsReal;
class RooRealVar;

class RooCubicSplineGaussModel : public RooResolutionModel, public RooAbsEffResModel {
public:

  // Constructors, assignment etc
  inline RooCubicSplineGaussModel() : _flatSFInt(kFALSE)  { }
  RooCubicSplineGaussModel(const char *name, const char *title, 
            RooRealVar& x, RooAbsGaussModelEfficiency& spline,
            RooAbsReal& mean,   RooAbsReal& sigma );
  RooCubicSplineGaussModel(const char *name, const char *title, 
            RooRealVar& x, RooAbsGaussModelEfficiency& spline,
            RooAbsReal& mean,   RooAbsReal& sigma,
            RooAbsReal& meanSF, RooAbsReal& sigmaSF) ; 
  RooCubicSplineGaussModel(const RooCubicSplineGaussModel& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooCubicSplineGaussModel(*this,newname) ; }
  virtual ~RooCubicSplineGaussModel();
  
  virtual Int_t basisCode(const char* name) const ;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName) const ;

  void advertiseFlatScaleFactorIntegral(Bool_t flag) { _flatSFInt = flag ; }

  virtual const RooAbsReal* efficiency() const;
  virtual RooArgSet observables() const;

private:

  virtual Double_t evaluate() const ;

  std::complex<double> evalInt(Double_t xmin, Double_t xmax, 
                               Double_t scale, Double_t offset,
                               const std::complex<double>& z) const;
  Bool_t _flatSFInt ;
  
  RooRealProxy eff ;
  RooRealProxy mean ;
  RooRealProxy sigma ;
  RooRealProxy msf ;
  RooRealProxy ssf ;

  ClassDef(RooCubicSplineGaussModel,1)
};

#endif

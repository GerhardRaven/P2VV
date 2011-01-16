/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   JvL, Jeroen van Leerdam, Nikhef,           j.van.leerdam@nikhef.nl      *
 *   GR,  Gerhard Raven,      Nikhef & VU,      Gerhard.Raven@nikhef.nl      *
 *   PL,  Parker C Lund,      UC Irvine                                      *
 *   DK,  David Kirkby,       UC Irvine,        dkirkby@uci.edu              *
 *   WV,  Wouter Verkerke,    UC Santa Barbara, verkerke@slac.stanford.edu   *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROO_BTAGDECAY
#define ROO_BTAGDECAY

#include "RooAbsAnaConvPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"

class RooBTagDecay : public RooAbsAnaConvPdf
{

public:
  enum DecayType {SingleSided, DoubleSided, Flipped};

  inline RooBTagDecay() {}

  // constructor without flavour tags
  RooBTagDecay(const char *name, const char *title,
      RooRealVar& time, RooAbsReal& tau, RooAbsReal& dGamma, RooAbsReal& dm,
      RooAbsReal& coshCoef, RooAbsReal& sinhCoef, RooAbsReal& cosCoef,
      RooAbsReal& sinCoef, const RooResolutionModel& model, DecayType type);

  // constructor with initial state flavour tag (decay into CP eigenstate)
  RooBTagDecay(const char *name, const char* title,
    RooRealVar& time, RooAbsCategory& iTag, RooAbsReal& tau,
    RooAbsReal& dGamma, RooAbsReal& dm, RooAbsReal& dilution,
    RooAbsReal& ANorm, RooAbsReal& ATagEff, RooAbsReal& ADilMisTag,
    RooAbsReal& coshCoef, RooAbsReal& sinhCoef, RooAbsReal& cosCoef,
    RooAbsReal& sinCoef, const RooResolutionModel& model, DecayType type);

  RooBTagDecay(const RooBTagDecay& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const 
  { 
    return new RooBTagDecay(*this, newname);
  }

  virtual ~RooBTagDecay();

  virtual Double_t coefficient(Int_t basisIndex) const;
  RooArgSet*       coefVars(Int_t coefIdx) const ;

  Int_t getCoefAnalyticalIntegral(Int_t coef, RooArgSet& allVars,
      RooArgSet& analVars, const char* rangeName = 0) const ;
  Double_t coefAnalyticalIntegral(Int_t coef, Int_t code,
      const char* rangeName = 0) const ;

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars,
      Bool_t staticInitOK = kTRUE) const;
  void generateEvent(Int_t code);

protected:
  void declareBases();

  RooRealProxy     _time;
  RooCategoryProxy _iTag;
  RooCategoryProxy _fTag;
  RooRealProxy     _tau;
  RooRealProxy     _dGamma;
  RooRealProxy     _dm;
  RooRealProxy     _dilution;
  RooRealProxy     _ANorm;
  RooRealProxy     _ATagEff;
  RooRealProxy     _ADilMisTag;
  RooRealProxy     _coshCoef;
  RooRealProxy     _sinhCoef;
  RooRealProxy     _cosCoef;
  RooRealProxy     _sinCoef;
  Int_t            _coshBasis;
  Int_t            _sinhBasis;
  Int_t            _cosBasis;
  Int_t            _sinBasis;
  DecayType        _type;
  Int_t            _tags;

  ClassDef(RooBTagDecay, 1) // PDF of B decay time distribution with flavour tags
};

#endif


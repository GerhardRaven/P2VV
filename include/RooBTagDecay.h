/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl           *
 *   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl           *
 *                                                                           *
 * Copyright (c) 2011, Nikhef. All rights reserved.                          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROO_B_TAG_DECAY
#define ROO_B_TAG_DECAY

#include "RooAbsAnaConvPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooListProxy.h"

class RooCategory;
class RooRealVar;

class RooBTagDecay : public RooAbsAnaConvPdf
{

public:
  enum DecayType {SingleSided, DoubleSided, Flipped};

  inline RooBTagDecay() {}

  // constructor without flavour tags
  RooBTagDecay(const char *name, const char *title,
    RooRealVar& time, RooAbsReal& tau, RooAbsReal& dGamma, RooAbsReal& dm,
    RooAbsReal& coshCoef, RooAbsReal& sinhCoef, RooAbsReal& cosCoef,
    RooAbsReal& sinCoef, const RooResolutionModel& model, DecayType type,
    Bool_t checkTags = kTRUE);

  // constructor with both initial and final state flavour tags
  RooBTagDecay(const char *name, const char* title,
    RooRealVar& time, RooCategory& iTag, RooCategory& fTag,
    RooAbsReal& tau, RooAbsReal& dGamma, RooAbsReal& dm, RooAbsReal& dilution,
    RooAbsReal& ADilWTag, RooAbsReal& ANorm, RooAbsReal& avgCEven,
    RooAbsReal& avgCOdd, RooAbsReal& cosCoef, const RooResolutionModel& model,
    DecayType type, Bool_t checkTags = kTRUE);

  // constructor with only an initial state flavour tag
  RooBTagDecay(const char *name, const char* title,
    RooRealVar& time, RooCategory& iTag, RooAbsReal& tau,
    RooAbsReal& dGamma, RooAbsReal& dm, RooAbsReal& dilution,
    RooAbsReal& ADilWTag, RooAbsReal& avgCEven, RooAbsReal& avgCOdd,
    RooAbsReal& coshCoef, RooAbsReal& sinhCoef, RooAbsReal& cosCoef,
    RooAbsReal& sinCoef, const RooResolutionModel& model, DecayType type,
    Bool_t checkTags = kTRUE);

  // constructor with initial and final state tags and tagging categories
  RooBTagDecay(const char *name, const char* title,
    RooRealVar& time, RooCategory& iTag, RooCategory& fTag,
    RooCategory& tagCat, RooAbsReal& tau, RooAbsReal& dGamma,
    RooAbsReal& dm, RooArgList& dilutions, RooArgList& ADilWTags,
    RooAbsReal& ANorm, RooArgList& avgCEvens, RooArgList& avgCOdds,
    RooArgList& tagCatCoefs, RooAbsReal& cosCoef,
    const RooResolutionModel& model, DecayType type, Bool_t checkTags = kTRUE);

  // constructor with an initial state tag and tagging categories
  RooBTagDecay(const char *name, const char* title,
    RooRealVar& time, RooCategory& iTag, RooCategory& tagCat, RooAbsReal& tau,
    RooAbsReal& dGamma, RooAbsReal& dm, RooArgList& dilutions,
    RooArgList& ADilWTags, RooArgList& avgCEvens, RooArgList& avgCOdds,
    RooArgList& tagCatCoefs, RooAbsReal& coshCoef, RooAbsReal& sinhCoef,
    RooAbsReal& cosCoef, RooAbsReal& sinCoef, const RooResolutionModel& model,
    DecayType type, Bool_t checkTags = kTRUE);

  RooBTagDecay(const RooBTagDecay& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const 
  { 
    return new RooBTagDecay(*this, newname);
  }

  virtual ~RooBTagDecay();

  virtual Double_t coefficient(Int_t basisIndex) const;
  Double_t         tagCatCoef() const;
  Double_t         tagCatCoef(Int_t& category) const;
  RooArgSet*       coefVars(Int_t coefIdx) const;

  Int_t getCoefAnalyticalIntegral(Int_t coef, RooArgSet& allVars,
      RooArgSet& analVars, const char* rangeName = 0) const;
  Double_t coefAnalyticalIntegral(Int_t coef, Int_t code,
      const char* rangeName = 0) const;

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars,
      Bool_t staticInitOK = kTRUE) const;
  void generateEvent(Int_t code);

  Int_t getTagCatPosition(Int_t tagCatIndex) const;
  Int_t getTagCatIndex(Int_t tagCatPosition) const;

protected:
  void initTaggingCats(RooArgList& tagCatCoefs, RooArgList& dilutions,
      RooArgList& ADilWTags, RooArgList& avgCEvens, RooArgList& avgCOdds);
  void initTag(Bool_t iTag);
  void initTagCatMaps() const;
  void declareBases();

  Bool_t checkVarDep(const RooAbsArg& var, Bool_t warn = kFALSE,
      Bool_t onlyTagPars = kFALSE) const;
  Bool_t checkTag(Bool_t iTag) const;

  RooRealProxy     _time;
  RooCategoryProxy _iTag;
  RooCategoryProxy _fTag;
  RooCategoryProxy _tagCat;

  RooRealProxy _tau;
  RooRealProxy _dGamma;
  RooRealProxy _dm;

  RooListProxy _dilutions;
  RooListProxy _ADilWTags;
  RooRealProxy _ANorm;
  RooRealProxy _avgCEvenSum;
  RooRealProxy _avgCOddSum;
  RooListProxy _avgCEvens;
  RooListProxy _avgCOdds;
  RooListProxy _tagCatCoefs;

  RooRealProxy _coshCoef;
  RooRealProxy _sinhCoef;
  RooRealProxy _cosCoef;
  RooRealProxy _sinCoef;

  mutable std::map<Int_t, Int_t> _tagCatPositions; //!
  mutable std::map<Int_t, Int_t> _tagCatIndices;   //!

  Int_t _coshBasis;
  Int_t _sinhBasis;
  Int_t _cosBasis;
  Int_t _sinBasis;

  DecayType    _decayType;
  Int_t        _tagCatType;
  Int_t        _tags;
  Int_t        _iTagVal;
  Int_t        _fTagVal;
  Bool_t       _checkVars;
  RooListProxy _createdVars;

  ClassDef(RooBTagDecay, 1) // PDF of B decay time distribution with flavour tagging
};

#endif


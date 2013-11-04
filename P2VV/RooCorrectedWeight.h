/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                *
 *                                                                           *
 * Copyright (c) 2012, Nikhef. All rights reserved.                          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROO_CORRECTED_WEIGHT
#define ROO_CORRECTED_WEIGHT

#include <vector>
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooDataSet.h"

class RooAbsCategory;

class RooCorrectedWeight : public RooAbsReal
{

public:
  inline RooCorrectedWeight() {}

  RooCorrectedWeight(const char *name, const char *title,
      const RooDataSet& data, Double_t corrFactor);

  RooCorrectedWeight(const char *name, const char *title,
      RooAbsReal& origWeight, Double_t corrFactor);

  RooCorrectedWeight(const char *name, const char *title,
      const RooDataSet& data, RooAbsCategory& splitCat,
      const std::vector<Double_t> corrFactors);

  RooCorrectedWeight(const char *name, const char *title,
      RooAbsReal& origWeight, RooAbsCategory& splitCat,
      const std::vector<Double_t> corrFactors);

  RooCorrectedWeight(const RooCorrectedWeight& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const 
  { 
    return new RooCorrectedWeight(*this, newname);
  }

  ~RooCorrectedWeight();

  inline Double_t getVal(const RooArgSet* set = 0)  const {return evaluate();}
  inline Double_t getVal(const RooArgSet& set)      const {return evaluate();}
  inline Double_t getValV(const RooArgSet* set = 0) const {return evaluate();}

  Int_t position() const;
  Double_t correctionFactor() const {return _corrFactors.at(position());}
  Double_t origWeight() const {
    return _data != 0 ? _data->weight() : (Double_t)_origWeight;
  }

private:
  const RooDataSet* _data;
  RooRealProxy _origWeight;
  RooCategoryProxy _splitCat;
  std::vector<Double_t> _corrFactors;

  mutable std::map<Int_t, Int_t> _positionMap; //!

  Double_t evaluate() const;

  ClassDef(RooCorrectedWeight, 1) // corrects event weights (-ln(L)) for background dilution
};

#endif

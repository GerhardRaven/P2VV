/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                *
 *                                                                           *
 * Copyright (c) 2013, Nikhef. All rights reserved.                          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/


//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Multiplies a given event weight with a correction factor. Optionally, a
// correction factor is given for each state of the split category in a
// simultaneous fit.
// END_HTML
//

#include "Riostream.h"
#include "RooMsgService.h"
#include "RooAbsCategory.h"
#include "P2VV/RooCorrectedWeight.h"

//_____________________________________________________________________________
RooCorrectedWeight::RooCorrectedWeight(const char *name, const char* title,
    const RooDataSet& data, Double_t corrFactor) :
    RooAbsReal(name, title),
  _data(&data),
  _origWeight("origWeight", "original weight", this),
  _splitCat("splitCat", "split category", this),
  _corrFactors(std::vector<Double_t>(1, corrFactor))
{
  // constructor with data and without split category
}

//_____________________________________________________________________________
RooCorrectedWeight::RooCorrectedWeight(const char *name, const char* title,
    RooAbsReal& origWeight, Double_t corrFactor) :
    RooAbsReal(name, title),
  _data(0),
  _origWeight("origWeight", "original weight", this, origWeight),
  _splitCat("splitCat", "split category", this),
  _corrFactors(std::vector<Double_t>(1, corrFactor))
{
  // constructor with original weight and without split category
}

//_____________________________________________________________________________
RooCorrectedWeight::RooCorrectedWeight(const char *name, const char* title,
    const RooDataSet& data, RooAbsCategory& splitCat,
    const std::vector<Double_t> corrFactors) :
  RooAbsReal(name, title),
  _data(&data),
  _origWeight("origWeight", "original weight", this),
  _splitCat("splitCat", "split category", this, splitCat),
  _corrFactors(corrFactors)
{
  // constructor with data and with split category

  // check number of correction factors
  assert((Int_t)corrFactors.size() == splitCat.numTypes());
}

//_____________________________________________________________________________
RooCorrectedWeight::RooCorrectedWeight(const char *name, const char* title,
    RooAbsReal& origWeight, RooAbsCategory& splitCat,
    const std::vector<Double_t> corrFactors) :
  RooAbsReal(name, title),
  _data(0),
  _origWeight("origWeight", "original weight", this, origWeight),
  _splitCat("splitCat", "split category", this, splitCat),
  _corrFactors(corrFactors)
{
  // constructor with original weight and with split category

  // check number of correction factors
  assert((Int_t)corrFactors.size() == splitCat.numTypes());
}

//_____________________________________________________________________________
RooCorrectedWeight::RooCorrectedWeight(
    const RooCorrectedWeight& other, const char* name) :
  RooAbsReal(other, name),
  _data(other._data),
  _origWeight("origWeight", this, other._origWeight),
  _splitCat("splitCat", this, other._splitCat),
  _corrFactors(other._corrFactors),
  _positionMap(other._positionMap)
{
  // copy constructor
}

//_____________________________________________________________________________
RooCorrectedWeight::~RooCorrectedWeight()
{
  // destructor
}

Int_t RooCorrectedWeight::position() const
{
  if (_splitCat.absArg() != 0 && _positionMap.empty()) {
    // create map from split category indices to positions in factors vector
    TIterator*  catIter  = _splitCat.arg().typeIterator();
    RooCatType* catType  = 0;
    Int_t       position = 0;
    while ((catType = (RooCatType*)catIter->Next()) != 0) {
      Int_t catIndex = catType->getVal();
      _positionMap[catIndex] = position;
      ++position;
    }

    delete catIter;
  }

  // return position
  return _splitCat.absArg() != 0
      ? _positionMap.find(_splitCat.arg().getIndex())->second : 0;
}

//_____________________________________________________________________________
Double_t RooCorrectedWeight::evaluate() const
{
  // calculate and return the corrected signal sWeight
  _value = origWeight() * correctionFactor();
  return _value;
}

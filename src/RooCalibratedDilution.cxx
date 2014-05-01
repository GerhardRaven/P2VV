/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * Authors:                                                                  *
 *   RA,  Roel Aaij,          Nikhef                                         *
 *   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl           *
 *                                                                           *
 * Copyright (c) 2012, Nikhef. All rights reserved.                          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/ 

#include "RooAbsReal.h" 
#include "RooMsgService.h"
#include "P2VV/RooCalibratedDilution.h" 
#include "TMath.h"

//_____________________________________________________________________________
RooCalibratedDilution::RooCalibratedDilution(const char* name,
    const char* title, RooAbsReal& WTag) :
  RooAbsReal(name, title),
  _estWTag("WTag", "CP average wrong-tag probability", this, WTag),
  _avgEstWTag("avgEstWTag", "", this),
  _P0("P0", "", this),
  _P1("P1", "", this),
  _AP0("AP0", "", this),
  _AP1("AP1", "", this)
{}

//_____________________________________________________________________________
RooCalibratedDilution::RooCalibratedDilution(const char* name,
    const char* title, RooAbsReal& WTag, RooAbsReal& AWTag) :
  RooAbsReal(name, title),
  _estWTag("WTag", "CP average wrong-tag probability", this, WTag),
  _avgEstWTag("avgEstWTag", "", this),
  _P0("P0", "", this),
  _P1("P1", "", this),
  _AP0("AP0", "", this),
  _AP1("AWTag", "Wrong-tag probability CP asymmetry", this, AWTag)
{}

//_____________________________________________________________________________
RooCalibratedDilution::RooCalibratedDilution(const char* name,
    const char* title, RooAbsReal& estWTag, RooAbsReal& avgEstWTag,
    RooAbsReal& P0, RooAbsReal& P1) :
  RooAbsReal(name, title),
  _estWTag("estWTag", "Estimated wrong tag probability", this, estWTag),
  _avgEstWTag("avgEstWTag", "Mean estimated wrong tag", this, avgEstWTag),
  _P0("P0", "CP average calibration parameter P_0", this, P0),
  _P1("P1", "CP average calibration parameter P_1", this, P1),
  _AP0("AP0", "", this),
  _AP1("AP1", "", this)
{}

//_____________________________________________________________________________
RooCalibratedDilution::RooCalibratedDilution(const char* name,
    const char* title, RooAbsReal& estWTag, RooAbsReal& avgEstWTag,
    RooAbsReal& P0, RooAbsReal& P1, RooAbsReal& AP0, RooAbsReal& AP1) :
  RooAbsReal(name, title),
  _estWTag("estWTag", "Estimated wrong tag probability", this, estWTag),
  _avgEstWTag("avgEstWTag", "Mean estimated wrong tag", this, avgEstWTag),
  _P0("P0", "CP average calibration parameter P_0", this, P0),
  _P1("P1", "CP average calibration parameter P_1", this, P1),
  _AP0("AP0", "Calibration parameter P_0 CP asymmetry", this, AP0),
  _AP1("AP1", "Calibration parameter P_1", this, AP1)
{}

//_____________________________________________________________________________
RooCalibratedDilution::RooCalibratedDilution(
    const RooCalibratedDilution& other, const char* name) :
  RooAbsReal(other, name),
  _estWTag("estWTag", this, other._estWTag),
  _avgEstWTag("avgEstWTag", this, other._avgEstWTag),
  _P0("P0", this, other._P0),
  _P1("P1", this, other._P1),
  _AP0("AP0", this, other._AP0),
  _AP1("AP1", this, other._AP1)
{}

//_____________________________________________________________________________
Double_t RooCalibratedDilution::WTag() const
{
  if (_P0.absArg() != 0) return _P0 + _P1 * (_estWTag - _avgEstWTag);
  else return _estWTag;
}

Double_t RooCalibratedDilution::estWTag() const
{
  if (_P0.absArg() != 0) return _estWTag;
  else return -999.;
}

Double_t RooCalibratedDilution::avgEstWTag() const
{
  if (_P0.absArg() != 0) return _avgEstWTag;
  else return -999.;
}

Double_t RooCalibratedDilution::P0() const
{
  if (_P0.absArg() != 0) return _P0;
  else return -999.;
}

Double_t RooCalibratedDilution::P1() const
{
  if (_P0.absArg() != 0) return _P1;
  else return -999.;
}

Double_t RooCalibratedDilution::AP0() const
{
  if (_AP0.absArg() != 0) return _AP0;
  else return -999.;
}

Double_t RooCalibratedDilution::AP1() const
{
  if (_AP0.absArg() != 0) return _AP1;
  else return -999.;
}


//_____________________________________________________________________________
Double_t RooCalibratedDilution::evaluate() const
{
  // calculate dilution factor
  Double_t dil(1.);
  Double_t twoP1Eta(0.);
  if (_P0.absArg() == 0) {
    // without calibration parameters
    dil -= 2. * _estWTag;
  } else {
    // with calibration parameters
    twoP1Eta = 2. * _P1 * (_estWTag - _avgEstWTag);
    dil -= 2. * _P0 + twoP1Eta;
  }

  if (_AP1.absArg() != 0) {
    // check if dilution is equal to zero: asymmetry undefined
    if (TMath::Abs(dil) < 1.e-12) {
      coutF(Eval) << "RooCalibratedDilution::evaluate(" << GetName()
          << "): |dilution| < 10^-12: this should not happen..." << std::endl;
      assert(0);
    }

    // return dilution scaled wrong-tag asymmetry
    if (_P0.absArg() == 0) {
      // without calibration parameters
      return 2. * _estWTag * _AP1 / dil;
    } else {
      // with calibration parameters
      return (2. * _P0 * _AP0 + twoP1Eta * _AP1) / dil;
    }
  }

  // return dilution
  return dil;
}

/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * Authors:                                                                  *
 *   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl           *
 *   WH,  Wouter Hulsbergen,  Nikhef                                         *
 *                                                                           *
 * Copyright (c) 2011, Nikhef. All rights reserved.                          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROO_REAL_MOMENTS_H
#define ROO_REAL_MOMENTS_H

#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include <string>

//_____________________________________________________________________________
class RooAbsRealMoment
{

public:
  RooAbsRealMoment(RooAbsReal& basisFunc, Double_t norm = 1.,
      const std::string& name = std::string());
  virtual ~RooAbsRealMoment() {};

  const char* name() {return _name.c_str();}
  RooAbsReal& basisFunc() {return _basisFunc;}
  Double_t norm() {return _norm;}

  virtual RooArgSet* getObservables(const RooArgSet* set)
  {
    return basisFunc().getObservables(set);
  }
  RooArgSet* getObservables(const RooAbsData& data)
  {
    return getObservables(data.get());
  }

  virtual Double_t coefficient(Bool_t normalize = kTRUE) const;
  virtual Double_t variance(Bool_t normalize = kTRUE) const;
  virtual Double_t significance() const;
  virtual Double_t evaluate() {return _basisFunc.getVal();}
  virtual Double_t stdDev(Bool_t normalize = kTRUE) const;

  virtual void inc(Double_t weight = 1.);
  void reset() {_m0 = _m1 = _n0 = _n1 = _n2 = 0.;}

  virtual ostream& print(ostream& os, Bool_t normalize = kTRUE) const;
  void Print(Bool_t normalize = kTRUE) const {print(std::cout, normalize);}

protected:
  RooAbsReal& _basisFunc;
  std::string _name;
  Double_t _norm;
  Double_t _m0, _m1;
  Double_t _n0, _n1, _n2;
};

//_____________________________________________________________________________
class RooRealMoment : public RooAbsRealMoment
{

public:
  RooRealMoment(RooAbsReal& basisFunc, Double_t norm = 1.);

private:
};

//_____________________________________________________________________________
class RooRealEffMoment : public RooAbsRealMoment
{

public:
  RooRealEffMoment(RooAbsReal& basisFunc, Double_t norm, const RooAbsPdf& pdf,
      const RooArgSet& normSet);

  Double_t evaluate() {return _basisFunc.getVal() / _pdf.getVal(&_normSet);}
  virtual RooArgSet* getObservables(const RooArgSet* set)
  {
    return _pdf.getObservables(set);
  }

private:
  const RooAbsPdf& _pdf;
  const RooArgSet& _normSet;
};

typedef std::vector<RooAbsRealMoment*> RooRealMomentsVector;

Int_t computeRooRealMoments(RooAbsData& data, RooRealMomentsVector& moments,
    Bool_t resetFirst = kFALSE, Bool_t verbose = kTRUE);

#endif

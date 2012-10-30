/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * Authors:                                                                  *
 *   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl           *
 *   WH,  Wouter Hulsbergen,  Nikhef                                         *
 *   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl           *
 *                                                                           *
 * Copyright (c) 2011, Nikhef. All rights reserved.                          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#include <cmath>
#include "ProgressDisplay.h"
#include "RooRealMoments.h"
#include "TMath.h"
#include "RooArgSet.h"

using namespace std;

//_____________________________________________________________________________

RooAbsRealMoment::RooAbsRealMoment(RooAbsReal& basisFunc, Double_t norm,
    const std::string& name) :
  _basisFunc(basisFunc), _name(name.empty() ? _basisFunc.GetName() : name),
  _norm(norm), _m0(0.), _m1(0.), _n0(0.), _n1(0.), _n2(0.) {}

Double_t RooAbsRealMoment::coefficient(Bool_t normalize) const
{
  if (_m0 <= 0.) return -999.;

  if (normalize) return _m1 / _m0 * _norm;
  else           return _m1 / _m0;
}

Double_t RooAbsRealMoment::variance(Bool_t normalize) const
{
  // the following formulas follow either from the jackknife method
  // or from error propagation using the following error on weight_j:
  //     sigma^2( weight_j ) = weight_j^2
  // (this is also exactly how it works with s-weight).

  // jackknife: sigma2 = (N - 1)/N * sum_j ( mj1 - m )^2), where mj1 is the
  // value of m if you would leave measurement j away

  // we make one approximation: we ignore the contribution of a
  // single weight to the total in a normalization term

  // var(mu) = 1/m0^2 * sum  w_j^2 (x_j - mu)^2
  if (_m0 <= 0.) return -999.;

  Double_t mu    = coefficient(kFALSE);
  Double_t varMu = (_n2 - 2. * _n1 * mu + _n0 * mu * mu) / (_m0 * _m0);

  if (normalize) return varMu * _norm * _norm;
  else           return varMu;
}

Double_t RooAbsRealMoment::stdDev(Bool_t normalize) const
{
  Double_t var = variance(normalize);
  if (TMath::Abs(var) < 1.e-10) return 0.;
  return var < 0. ? -999. : std::sqrt(var);
}

Double_t RooAbsRealMoment::significance() const
{
  Double_t mu  = coefficient(kFALSE);
  Double_t var = variance(kFALSE);
  if (TMath::Abs(var) < 1.e-10) return TMath::Infinity();
  return var < 0. ? -999. : std::sqrt(mu * mu / var);
}

void RooAbsRealMoment::inc(Double_t weight)
{
  Double_t funcVal = evaluate();

  // TODO: make a histogram of funcVal... (two, one for accept, one for all)
  _m0 += weight;
  _m1 += weight * funcVal;

  // these we need to compute the error using the jackknife method
  _n0 += weight * weight;
  _n1 += weight * weight * funcVal;
  _n2 += weight * weight * funcVal * funcVal;
}

ostream& RooAbsRealMoment::print(ostream& os, Bool_t normalize) const
{
  return os << "moment(" << _name << ") = " << coefficient(normalize) << " +- "
      << stdDev(normalize) << " (significance: " << significance() << ")"
      << endl;
}

//_____________________________________________________________________________

RooRealMoment::RooRealMoment(RooAbsReal& basisFunc, Double_t norm) :
    RooAbsRealMoment(basisFunc, norm) {}


//_____________________________________________________________________________

RooRealEffMoment::RooRealEffMoment(RooAbsReal& basisFunc, Double_t norm,
      const RooAbsPdf& pdf, const RooArgSet& intSet, const RooArgSet& normSet) :
    RooAbsRealMoment(basisFunc, norm, std::string(basisFunc.GetName()) + "_"
        + pdf.GetName()), _pdf(pdf), _pdfInt(0), _intSet(intSet),
        _normSet(normSet)
{
  _pdfInt = _pdf.createIntegral(_intSet, _normSet);
}

//_____________________________________________________________________________

RooRealEffMoment::~RooRealEffMoment()
{
  if (_pdfInt) delete _pdfInt;
}

//_____________________________________________________________________________

Int_t computeRooRealMoments(RooAbsData& data, RooRealMomentsVector& moments,
    Bool_t resetFirst, Bool_t verbose)
{
  typedef RooRealMomentsVector::iterator RealMomIter;

  if (moments.empty()) {
    cout <<"P2VV - ERROR: computeRealMoments: moments vector is empty"<< endl;
    return -1;
  }

  if (resetFirst) {
    for (RealMomIter mom = moments.begin(); mom != moments.end(); ++mom)
      (*mom)->reset();
  }

  RooArgSet* obs = moments.front()->getObservables(data);
  
  if (verbose)
    cout << "P2VV - INFO: computeRealMoments: computing " << moments.size()
        << " moment(s) for data set '" << data.GetName() << "' with "
        << data.numEntries() << " events" << endl;

  Int_t dataIter = 0;
  ProgressDisplay *prog = verbose ? new ProgressDisplay(data.numEntries()) : 0;
  while (dataIter < data.numEntries()) {
    *obs = *data.get(dataIter++);
    for (RealMomIter mom = moments.begin(); mom != moments.end(); ++mom)
      (*mom)->inc(data.isWeighted() ? data.weight() : 1.);

    if (prog) ++*prog;
  }
  if (prog) { 
    cout << endl;
    delete prog;
  }
  return dataIter;
}

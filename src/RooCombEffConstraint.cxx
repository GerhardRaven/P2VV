/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   RA,  Roel Aaij,          Nikhef, roel.aaij@nikhef.nl                    *
 *   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                *
 *                                                                           *
 * Copyright (c) 2014, Nikhef. All rights reserved.                          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/


//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Constraint on selection-efficiency parameters from counting (weighted)
// events that pass different sets of selection criteria.
//
// Events are assumed to be selected in two stages, with two sets of selection
// criteria at each stage. The events that pass stage 1 are input to the
// selection at stage 2. It is assumed that the two selections at stage 1 (1A
// and 1B) are mutually exclusive. On the other hand, the selections at stage 2
// (2A and 2B) are assumed to be independent. That is, the probability for an
// event to be selected by 2A does not depend on whether the event was selected
// by 2B and vice versa.
//
// These assumptions give six categories of events that have passed the
// selections at both stages:
// 1) 1A; 2A and not 2B    2) 1A; 2B and not 2A    3) 1A; 2A and 2B
// 4) 1B; 2A and not 2B    5) 1B; 2B and not 2A    6) 1B; 2A and 2B
// The constraint was optimized for the case where the efficiency of selection
// 2A is close to 100%. Categories 2 and 5 will then have very few events and
// are merged.
//
// The event yields in the five remaining categories are constrained by
// assuming Poisson distributions, which are modified to have a mean equal to
// the sum of event weights in a category and a variance equal to the sum of
// squared event weights. The five yields are then expressed in terms of
// efficiency parameters:
// 1) nu_1A2A  = eps_1A * eps_2A * (1 - eps_2B) * nu
// 2) nu_2B    = (eps_1A + eps_1B) * (1 - eps_2A) * eps_2B * nu
// 3) nu_1A2AB = eps_1A * eps_2A * eps_2B * nu
// 4) nu_1B2A  = eps_1B * eps_2A * (1 - eps_2B) * nu
// 5) nu_1B2AB = eps_1B * eps_2A * eps_2B * nu
// END_HTML
//


#include <iostream>
#include <cmath> 

#include "Riostream.h" 
#include "TMath.h" 
#include "RooAbsArg.h"
#include "RooArgSet.h"
#include "RooArgList.h"

#include "P2VV/RooCombEffConstraint.h" 

namespace {
  using std::vector;
  using std::auto_ptr;
}

//_____________________________________________________________________________
RooCombEffConstraint::RooCombEffConstraint() {}

//_____________________________________________________________________________
RooCombEffConstraint::RooCombEffConstraint(const char *name, const char *title, 
    Int_t numBins, const RooArgList& nu, const RooArgList& eps1A,
    const RooArgList& eps1B, const RooArgList& eps2A, const RooArgList& eps2B,
    const vector< vector<Double_t> >& sumW,
    const vector< vector<Double_t> >& sumWSq, Int_t strategy) :
  RooAbsPdf(name, title),
  _strat(strategy),
  _nu("!nu", "nu", this),
  _eps1A("!eps1A", "epsilon 1A", this),
  _eps1B("!eps1B", "epsilon 1B", this),
  _eps2A("!eps2A", "epsilon 2A", this),
  _eps2B("!eps2B", "epsilon 2B", this),
  _nBins(numBins)
{
  assert(_strat == 0 || _strat == 1 || _strat == 2);

  assert(_nBins > 0 && (nu.getSize() == 1 || nu.getSize() == _nBins)
       && (eps1A.getSize() == 1 || eps1A.getSize() == _nBins)
       && (eps1B.getSize() == 1 || eps1B.getSize() == _nBins)
       && (eps2A.getSize() == 1 || eps2A.getSize() == _nBins)
       && (eps2B.getSize() == 1 || eps2B.getSize() == _nBins)
       && sumW.size() == 6 && sumWSq.size() == 6);

  for (Int_t bIt = 0; bIt < _nBins; ++bIt) {
    assert((nu.getSize()<=bIt || dynamic_cast<RooAbsReal*>(nu.at(bIt))!=0)
     && (eps1A.getSize()<=bIt || dynamic_cast<RooAbsReal*>(eps1A.at(bIt))!=0)
     && (eps1B.getSize()<=bIt || dynamic_cast<RooAbsReal*>(eps1B.at(bIt))!=0)
     && (eps2A.getSize()<=bIt || dynamic_cast<RooAbsReal*>(eps2A.at(bIt))!=0)
     && (eps2B.getSize()<=bIt || dynamic_cast<RooAbsReal*>(eps2B.at(bIt))!=0));
  }

  _nu.add(nu);
  _eps1A.add(eps1A);
  _eps1B.add(eps1B);
  _eps2A.add(eps2A);
  _eps2B.add(eps2B);

  for (Int_t cIt = 0; cIt < 6; ++cIt) {
    assert((Int_t)sumW[cIt].size() == _nBins
        && (Int_t)sumWSq[cIt].size() == _nBins);
    if ((_strat == 1 && cIt == 1) || (_strat == 2 && cIt == 3)) {
      assert((Int_t)sumW[_strat == 1 ? 4 : 5].size() == _nBins
          && (Int_t)sumWSq[_strat == 1 ? 4 : 5].size() == _nBins);
    } else if ((_strat == 1 && cIt == 4)
        || (_strat == 2 && (cIt == 4 || cIt == 5))) {
      continue;
    }

    vector<Double_t> sWVec;
    vector<Double_t> sWSqVec;
    Bool_t zeroWeight = kFALSE;
    for (Int_t bIt = 0; bIt < _nBins; ++bIt) {
      Double_t sW = sumW[cIt][bIt];
      Double_t sWSq = sumWSq[cIt][bIt];
      if (_strat == 1 && cIt == 1) {
        sW += sumW[4][bIt];
        sWSq += sumWSq[4][bIt];
      } else if (_strat == 2 && cIt == 3) {
        sW += sumW[5][bIt];
        sWSq += sumWSq[5][bIt];
      }
      assert((sW > 0.01 && sWSq > 0.01) || (sW == 0. && sWSq == 0.));
      if (sW == 0.) zeroWeight = kTRUE;
      sWVec.push_back(sW);
      sWSqVec.push_back(sWSq);
    }
    if (zeroWeight) {
      coutW(InputArguments) << "RooCombEffConstraint::ctor(" << GetName()
          << "): zero sums of (squared) weights found for category " << cIt + 1
          << ": constraining corresponding yield parameters to zero"
          <<  std::endl;
    }
    _sumW.push_back(sWVec);
    _sumWSq.push_back(sWSqVec);
    _normFac.push_back(vector<Double_t>(_nBins, -1.));
  }
}

//_____________________________________________________________________________
RooCombEffConstraint::RooCombEffConstraint(const RooCombEffConstraint& other,
    const char* name) :
  RooAbsPdf(other, name),
  _strat(other._strat),
  _nu("!nu", this, other._nu),
  _eps1A("!eps1A", this, other._eps1A),
  _eps1B("!eps1B", this, other._eps1B),
  _eps2A("!eps2A", this, other._eps2A),
  _eps2B("!eps2B", this, other._eps2B),
  _sumW(other._sumW),
  _sumWSq(other._sumWSq),
  _normFac(other._normFac),
  _nBins(other._nBins)
{} 

//_____________________________________________________________________________
RooCombEffConstraint::~RooCombEffConstraint() {}

//_____________________________________________________________________________
Bool_t RooCombEffConstraint::forceAnalyticalInt(const RooAbsArg& /*dep*/) const
{
  // force RooRealIntegral to offer all observables for internal integration
  return kTRUE;
}

//_____________________________________________________________________________
Int_t RooCombEffConstraint::getAnalyticalIntegral(RooArgSet& allVars,
    RooArgSet& analVars, const char* /*rangeName*/) const 
{
  coutI(Eval) << "RooCombEffConstraint::getAnalyticalIntegral(" << GetName()
      << "): called" << std::endl;
  analVars.add(allVars);
  return 1;
}

//_____________________________________________________________________________
Double_t RooCombEffConstraint::analyticalIntegral(Int_t code,
    const char* /*rangeName*/) const
{
  if (code != 1) {
    coutF(InputArguments) << "RooCombEffConstraint::analyticalIntegral("
        << GetName() << "): integration code should be 1 (got" << code << ")"
        << std::endl;
    assert(0);
  }

  return 1.;
}

//_____________________________________________________________________________
Double_t RooCombEffConstraint::getLogVal(const RooArgSet* nset) const
{
  Double_t lnL = 0.;
  for (Int_t bIt = 0; bIt < _nBins; ++bIt) {
    Double_t nuTot = _nu.getSize() == 1
        ? ((RooAbsReal*)_nu.at(0))->getVal()
        : ((RooAbsReal*)_nu.at(bIt))->getVal();
    Double_t e1A = _eps1A.getSize() == 1
        ? ((RooAbsReal*)_eps1A.at(0))->getVal()
        : ((RooAbsReal*)_eps1A.at(bIt))->getVal();
    Double_t e1B = _eps1B.getSize() == 1
        ? ((RooAbsReal*)_eps1B.at(0))->getVal()
        : ((RooAbsReal*)_eps1B.at(bIt))->getVal();
    Double_t e2A = _eps2A.getSize() == 1
        ? ((RooAbsReal*)_eps2A.at(0))->getVal()
        : ((RooAbsReal*)_eps2A.at(bIt))->getVal();
    Double_t e2B = _eps2B.getSize() == 1
        ? ((RooAbsReal*)_eps2B.at(0))->getVal()
        : ((RooAbsReal*)_eps2B.at(bIt))->getVal();

    Double_t nu2A  = e2A * (1. - e2B) * nuTot;
    Double_t nu2B  = (1. - e2A) * e2B * nuTot;
    Double_t nu2AB = e2A * e2B * nuTot;
    std::vector<Double_t> nu;
    if (_strat == 1) {
      nu.push_back(e1A * nu2A);
      nu.push_back((e1A + e1B) * nu2B);
      nu.push_back(e1A * nu2AB);
      nu.push_back(e1B * nu2A);
      nu.push_back(e1B * nu2AB);
    } else if (_strat == 2) {
      nu.push_back(e1A * nu2A);
      nu.push_back(e1A * nu2B);
      nu.push_back(e1A * nu2AB);
      nu.push_back(e1B * (nu2A + nu2AB));
    } else {
      nu.push_back(e1A * nu2A);
      nu.push_back(e1A * nu2B);
      nu.push_back(e1A * nu2AB);
      nu.push_back(e1B * nu2A);
      nu.push_back(e1B * nu2B);
      nu.push_back(e1B * nu2AB);
    }

    // construct log-likelihood sum for five nu values
    for (UInt_t cIt = 0; cIt < nu.size(); ++cIt) {
      Double_t sumWSq = _sumWSq[cIt][bIt];
      if (sumWSq > 0.) {
        // use modified poisson constraint for this category
        Double_t sumW = _sumW[cIt][bIt];
        Double_t nuN = nuNorm(nu[cIt], sumW, _normFac[cIt], bIt);
        lnL += sumW / sumWSq * (sumW * TMath::Log(nuN) - nu[cIt]);
      } else {
        // use "zero" constraint for this category
        lnL += -1.e10 * nu[cIt] * nu[cIt];
      }
    }
  }

  return lnL;
}

//_____________________________________________________________________________
Double_t RooCombEffConstraint::evaluate() const 
{ 
  return TMath::Exp(getLogVal());
} 

//_____________________________________________________________________________
Double_t RooCombEffConstraint::nuNorm(const Double_t nu, const Double_t sumW, 
    vector<Double_t>& normFac, const Int_t ind) const
{
  if (normFac[ind] < 0.) {
    if (nu > 0.) {
      normFac[ind] = TMath::Exp(nu / sumW) / nu;
    } else {
      normFac[ind] = TMath::E() / sumW;
    }
  }

  return normFac[ind] * nu;
}

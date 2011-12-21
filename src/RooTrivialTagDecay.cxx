/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooTrivialTagDecay.cxx 24286 2008-06-16 15:47:04Z wouter $
 * Authors:                                                                  *
 *   PL, Parker C Lund,   UC Irvine                                          *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/


//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Most general description of B decay time distribution with effects
// of CP violation, mixing and life time differences. This function can 
// be analytically convolved with any RooResolutionModel implementation
// END_HTML
//


#include "RooFit.h"

#include "Riostream.h"
#include "TMath.h"
#include "RooTrivialTagDecay.h"
#include "RooRealVar.h"
#include "RooProduct.h"
#include "RooArgSet.h"
#include "RooAbsCategory.h"
#include "RooRandom.h"

ClassImp(RooTrivialTagDecay);

namespace {

};

//_____________________________________________________________________________
// For now, we assume that the cos (cosh) and sin (sinh) coefficients are odd (even) wrt. tag 
// this is not quite right as soon as eg. |lambda|^2 != 1 -- then all coefficients
// get multiplied by  1 + _tag * C where C = (1-|lambda|^2)/(1+|lambda|^2)
// TODO: this can be incorporated by passing 8 coefficients, namely the
//       'even' and 'odd' for each of the 4 basis functions...
// TODO: multiply all coefficients by 
//   "(@0==0)*@1 + (@1!=0)*(1-@1)",RooArgList( _tag,_tageff )
RooTrivialTagDecay::RooTrivialTagDecay(const char *name, const char* title, 
	       RooRealVar& t, RooAbsCategory& tag, RooAbsReal& tau, RooAbsReal& dgamma, RooAbsReal& dm,  RooAbsReal& tageff,
	       RooAbsReal& fcosh /*even*/, RooAbsReal& fsinh /*even*/, RooAbsReal& fcos /*odd*/, RooAbsReal& fsin /*odd*/, 
           const RooResolutionModel& model, DecayType type) :
  RooAbsAnaConvPdf(name, title, model, t),
  _t("t", "time", this, t),
  _tag("tag", "tag", this, tag),
  _tau("tau", "Average Decay Time", this, tau),
  _dgamma("dgamma", "Delta Gamma", this, dgamma),
  _dm("dm", "Delta Mass", this, dm),
  _tageff("tageff","tageff",this,tageff),
  _fcosh("fcosh", "Cosh Coefficient", this, fcosh), 
  _fsinh("fsinh", "Sinh Coefficient", this, fsinh),
  _fcos("fcos", "q x Cos Coefficient", this, *new RooProduct("__tag_fcos__","__tag_fcos__", RooArgSet(tag,fcos))), // ,kTRUE,kFALSE,kTRUE),
  _fsin("fsin", "q x Sin Coefficient", this, *new RooProduct("__tag_fsin__","__tag_fsin__", RooArgSet(tag,fsin))), // ,kTRUE,kFALSE,kTRUE),
  _type(type)
{
  //Constructor
  switch(type)
    {
    case SingleSided:
      _basisCosh = declareBasis("exp(-@0/@1)*cosh(@0*@2/2)", RooArgList(tau,dgamma));
      _basisSinh = declareBasis("exp(-@0/@1)*sinh(@0*@2/2)", RooArgList(tau,dgamma));
      _basisCos = declareBasis("exp(-@0/@1)*cos(@0*@2)",RooArgList(tau, dm));
      _basisSin = declareBasis("exp(-@0/@1)*sin(@0*@2)",RooArgList(tau, dm));
      break;
    case Flipped:
      _basisCosh = declareBasis("exp(@0/@1)*cosh(@0*@2/2)", RooArgList(tau,dgamma));
      _basisSinh = declareBasis("exp(@0/@1)*sinh(@0*@2/2)", RooArgList(tau,dgamma));
      _basisCos = declareBasis("exp(@0/@1)*cos(@0*@2)",RooArgList(tau, dm));
      _basisSin = declareBasis("exp(@0/@1)*sin(@0*@2)",RooArgList(tau, dm));
      break;
    case DoubleSided:
      _basisCosh = declareBasis("exp(-abs(@0)/@1)*cosh(@0*@2/2)", RooArgList(tau,dgamma));
      _basisSinh = declareBasis("exp(-abs(@0)/@1)*sinh(@0*@2/2)", RooArgList(tau,dgamma));
      _basisCos = declareBasis("exp(-abs(@0)/@1)*cos(@0*@2)",RooArgList(tau, dm));
      _basisSin = declareBasis("exp(-abs(@0)/@1)*sin(@0*@2)",RooArgList(tau, dm));
      break;
    }
  cout << " RooTrivialTagDecay("<< this <<")::ctor" << endl;
}

//_____________________________________________________________________________
RooTrivialTagDecay::RooTrivialTagDecay(const RooTrivialTagDecay& other, const char* name) :
  RooAbsAnaConvPdf(other, name),
  _t("t", this, other._t),
  _tag("tag", this, other._tag),
  _tau("tau", this, other._tau),
  _dgamma("dgamma", this, other._dgamma),
  _dm("dm", this, other._dm),
  _tageff("tageff","tageff",this,other._tageff),
  _fcosh("fcosh", this, other._fcosh),
  _fsinh("fsinh", this, other._fsinh),
  _fcos("fcos", this, other._fcos),
  _fsin("fsin", this, other._fsin),
  _basisCosh(other._basisCosh),
  _basisSinh(other._basisSinh),
  _basisCos(other._basisCos),
  _basisSin(other._basisSin),
  _type(other._type)
{
  cout << " RooTrivialTagDecay("<< this <<")::copy ctor("<< &other <<"," << (name ? name : "<none>" )<< ")" << endl;
  //Copy constructor
}

TObject* RooTrivialTagDecay::clone(const char* newname) const 
{ 
    cout << "RooTrivialTagDecay("<< this <<")::clone(" << (newname? newname : "<none>" ) << ")" << endl;
    RooTrivialTagDecay *p = new RooTrivialTagDecay(*this,newname);
    cout << "RooTrivialTagDecay("<< this <<")::clone -- return "  << p << endl;
    return p;
}


//_____________________________________________________________________________
RooTrivialTagDecay::~RooTrivialTagDecay()
{
  //Destructor
}


//_____________________________________________________________________________
Double_t RooTrivialTagDecay::coefficient(Int_t basisIndex) const
{
  const RooRealProxy* p = proxy( basisIndex );
  return (p ? double(*p) : 0. ) * ( _tag ? _tageff : 1.-_tageff );
}


//_____________________________________________________________________________
RooArgSet* RooTrivialTagDecay::coefVars(Int_t basisIndex) const 
{
  const RooRealProxy* p = proxy( basisIndex );
  return  p ? p->arg().getVariables() : 0 ;
}



//_____________________________________________________________________________
Int_t RooTrivialTagDecay::getCoefAnalyticalIntegral(Int_t coef, RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
  // Note: what if we are integrating over _tag???? remove _tag from allVars!!!
  RooArgSet redVars(allVars);
  redVars.remove(_tag.arg());
  const RooRealProxy* p = proxy( coef );
  return p ? p->arg().getAnalyticalIntegral(redVars,analVars,rangeName) : 0; 
}

//_____________________________________________________________________________
Double_t RooTrivialTagDecay::coefAnalyticalIntegral(Int_t coef, Int_t code, const char* rangeName) const 
{
  // TODO: what if we are integrating over _tag???? Should not happen, as getCoefAnalyticalIntegral will not select it!
  const RooRealProxy* p = proxy( coef );
  return ( p ? p->arg().analyticalIntegral(code,rangeName) : 0 ) * ( _tag!=0 ? _tageff : 1.-_tageff );
}



//_____________________________________________________________________________
Int_t RooTrivialTagDecay::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  if (matchArgs(directVars, generateVars, _t,_tag)) return 1;
  if (matchArgs(directVars, generateVars, _t)) return 2;
  if (matchArgs(directVars, generateVars, _tag)) return 3;
  return 0;
}


//_____________________________________________________________________________
void RooTrivialTagDecay::generateEvent(Int_t code)
{
  Double_t gammamin = 1/_tau-TMath::Abs(_dgamma)/2;
  while(1) {
    // first untagged...
    if (code == 1 || code == 2 ) {
        Double_t t = -log(RooRandom::uniform())/gammamin;
        if (_type == Flipped || (_type == DoubleSided && RooRandom::uniform() <0.5) ) t *= -1;
        if ( t<_t.min() || t>_t.max() ) continue;

        Double_t dgt = _dgamma*t/2;
        Double_t ft = fabs(t);
        Double_t even = _fcosh*cosh(dgt)+_fsinh*sinh(dgt);
        Double_t f = exp(-ft/_tau)*even;
        if(f < 0) throw(string(Form( "RooTrivialTagDecay::generateEvent(%s) ERROR: PDF value less than zero" ,GetName())));
        Double_t w = 1.001*exp(-ft*gammamin)*(TMath::Abs(_fcosh)+TMath::Abs(_fsinh));
        if(w < f) throw(string(Form( "RooTrivialTagDecay::generateEvent(%s) ERROR: Envelope function less than p.d.f. " ,GetName())));
        if(w*RooRandom::uniform() > f) continue;
        _t = t;
    }
    // and now for the tagging...
    Double_t dgt = _dgamma*_t/2;
    Double_t dmt = _dm*_t;
    int tag = RooRandom::uniform() < _tageff  ? 0 : (2*RooRandom::uniform() > 1.+fabs(_fcos *cos (dmt)+_fsin *sin (dmt))/(_fcosh*cosh(dgt)+_fsinh*sinh(dgt) )  ) > 0 ? +1 : -1 ;
    if ( code == 2 && _tag != tag ) continue; // hit-miss on required answer...
    _tag = tag;
    break;
  }
}

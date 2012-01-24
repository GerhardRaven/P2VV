/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
 * Authors:                                                                  *
 *   GR, Gerhard Raven, VU Amsterdam & NIKHEF, Gerhard.Raven@nikhef.nl       *
 *                                                                           *
 * Copyright (c) 2012                                                        *
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

#include <memory>

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

//_____________________________________________________________________________
// For now, we assume that the cos (cosh) and sin (sinh) coefficients are odd (even) wrt. tag 
// this is not quite right as soon as eg. |lambda|^2 != 1 -- then all coefficients
// get multiplied by  1 + _tag * C where C = (1-|lambda|^2)/(1+|lambda|^2)
// TODO: this can be incorporated by passing 8 coefficients, namely the
//       'even' and 'odd' for each of the 4 basis functions...
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
  _fcos("fcos", "q x Cos Coefficient", this ),
  _fsin("fsin", "q x Sin Coefficient", this ),
  _type(type)
{
  string n_q_fcos = Form("%s_q_cos",GetName()) ;
  string n_q_fsin = Form("%s_q_sin",GetName()) ;

  RooProduct *q_fcos = new RooProduct(n_q_fcos.c_str(),n_q_fcos.c_str(), RooArgSet(tag,fcos)); 
  addOwnedComponents(*q_fcos);
  _fcos.setArg(*q_fcos);
  RooProduct *q_fsin = new RooProduct(n_q_fsin.c_str(),n_q_fsin.c_str(), RooArgSet(tag,fsin));
  addOwnedComponents(*q_fsin);
  _fsin.setArg(*q_fsin);
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
}

//_____________________________________________________________________________
RooTrivialTagDecay::RooTrivialTagDecay(const RooTrivialTagDecay& other, const char* name) :
  RooAbsAnaConvPdf(other, name),
  _t("t", this, other._t),
  _tag("tag", this, other._tag),
  _tau("tau", this, other._tau),
  _dgamma("dgamma", this, other._dgamma),
  _dm("dm", this, other._dm),
  _tageff("tageff",this,other._tageff),
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
  //Copy constructor
}

TObject* RooTrivialTagDecay::clone(const char* newname) const 
{ 
    return new RooTrivialTagDecay(*this,newname);
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
  assert(p!=0);
  return  (p!=0 ? double(*p) : 0. ) * ( int(_tag)!=0 ? double(_tageff)/2 : 1.-double(_tageff) );
}


//_____________________________________________________________________________
RooArgSet* RooTrivialTagDecay::coefVars(Int_t basisIndex) const 
{
  const RooRealProxy* p = proxy( basisIndex );
  assert(p!=0);
  return  p ? p->arg().getVariables() : 0 ;
}



//_____________________________________________________________________________
Int_t RooTrivialTagDecay::getCoefAnalyticalIntegral(Int_t coef, RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
  assert(rangeName==0);
  // Note: what if we are integrating over _tag???? remove _tag from allVars!!!
  RooArgSet redVars(allVars);
  redVars.remove(_tag.arg());
  const RooRealProxy* p = proxy( coef );
  assert(p!=0);
  return p ? p->arg().getAnalyticalIntegral(redVars,analVars,rangeName) : 0; 
}

//_____________________________________________________________________________
Double_t RooTrivialTagDecay::coefAnalyticalIntegral(Int_t coef, Int_t code, const char* rangeName) const 
{
  assert(rangeName==0);
  // TODO: what if we are integrating over _tag???? Should not happen, as getCoefAnalyticalIntegral will not select it!
  const RooRealProxy* p = proxy( coef );
  assert(p!=0);
  return ( p ? p->arg().analyticalIntegral(code,rangeName) : 0 ) * ( int(_tag)!=0 ? double(_tageff)/2 : 1.-double(_tageff) );
}

//_____________________________________________________________________________
Bool_t RooTrivialTagDecay::isDirectGenSafe(const RooAbsArg& arg) const 
{
  // Check if given observable can be safely generated using the
  // pdfs internal generator mechanism (if that existsP). Observables
  // on which a PDF depends via more than route are not safe
  // for use with internal generators because they introduce
  // correlations not known to the internal generator

  // Arg must be direct server of self
  if (!findServer(arg.GetName())) return kFALSE ;

  // There must be no other dependency routes
  std::auto_ptr<TIterator> sIter( serverIterator()) ;
  const RooAbsArg *server = 0;
  while((server=(const RooAbsArg*)sIter->Next())) {
    if(server == &arg) continue;
    if(server->dependsOn(arg)) {
      if ( &arg == &_tag.arg() && ( server == &_fcos.arg() || server == &_fsin.arg() ) ) continue;
      return kFALSE ;
    }
  }
  return kTRUE ;
}


//_____________________________________________________________________________
Int_t RooTrivialTagDecay::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  if (matchArgs(directVars, generateVars, _t,_tag)) return 1;
  // if (matchArgs(directVars, generateVars, _t)) return 2;
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
        Double_t ft = TMath::Abs(t);
        Double_t even = _fcosh*cosh(dgt)+_fsinh*sinh(dgt);
        Double_t f = exp(-ft/_tau)*even;
        if(f < 0) throw(string(Form( "RooTrivialTagDecay::generateEvent(%s) ERROR: PDF value less than zero" ,GetName())));
        Double_t w = 1.001*exp(-ft*gammamin)*(TMath::Abs(_fcosh)+TMath::Abs(_fsinh));
        if(w < f) throw(string(Form( "RooTrivialTagDecay::generateEvent(%s) ERROR: Envelope function less than p.d.f. " ,GetName())));
        if(w*RooRandom::uniform() > f) continue;
        _t = t;
    }
    // and now for the tagging...
    if ( RooRandom::uniform() > double(_tageff)) {
         _tag = 0;
         break;
    }
    _tag = 1;
    Double_t dgt = _dgamma*_t/2;
    Double_t dmt = _dm*_t;
    Double_t o =  _fcos *cos (dmt)+_fsin *sin (dmt);
    Double_t e =  _fcosh*cosh(dgt)+_fsinh*sinh(dgt);
    _tag =  ( 2*e*RooRandom::uniform() < (e+o) ) > 0 ? +1 : -1 ;
    break;
  }
}

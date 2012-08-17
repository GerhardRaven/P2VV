/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooEffResModel.cxx 44982 2012-07-10 08:36:13Z moneta $
 * Authors:                                                                  *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Class RooEffResModel implements a RooResolutionModel that models a EffResian
// distribution. Object of class RooEffResModel can be used
// for analytical convolutions with classes inheriting from RooAbsAnaConvPdf
// END_HTML
//

#include "RooFit.h"
#include "Riostream.h"
#include "RooEffResModel.h"
#include "RooRealConstant.h"
#include "RooCustomizer.h"
#include "RooAddition.h"

using namespace std;

ClassImp(RooEffResModel) 
;

RooEffResModel::CacheElem::~CacheElem()
{
   delete _I;
   delete _clone;
   for (std::vector<RooCustomizer*>::const_iterator it = _customizers.begin(),
           end = _customizers.end(); it != end; ++it) {
      delete *it;
   }
}

RooArgList RooEffResModel::CacheElem::containedArgs(Action) 
{
   // Return list of all RooAbsArgs in cache element
   RooArgList l(_intObs);
   l.add(*_I);
   l.add(*_clone);
   return l;
}

RooEffResModel::CacheElem::CacheElem( const RooEffResModel& parent, const RooArgSet& iset, const TNamed* rangeName )
{
   RooRealVar& x = parent.convVar(); // binboundaries not const...
   const RooAbsReal& eff = parent.eff();
   const RooAbsReal& model = parent.model();
   // the subset of iset on which the efficiency depends
   std::auto_ptr<const RooArgSet> effInt( eff.getObservables(iset) ); 

   assert(effInt->getSize() < 2); // for now, we only do 1D efficiency histograms...
   if (effInt->getSize()==0) {
      _I = parent.model().createIntegral(iset,RooNameReg::str(rangeName)); 
      return;
   }

   Double_t xmin = x.getMin(RooNameReg::str(rangeName));
   Double_t xmax = x.getMax(RooNameReg::str(rangeName));

   RooArgList effList;
   RooArgList intList;

   std::list<Double_t>* bounds = eff.binBoundaries(x, x.getMin(), x.getMax());
   std::list<Double_t>::const_iterator lo, hi = bounds->begin();
   for (unsigned int i=0; i + 1 < bounds->size();++i ) {
      lo = hi++;
      if (*hi < xmin) continue; // not there yet...
      if (*lo > xmax) break;    // past the requested interval...
      Double_t thisxmin = std::max(*lo, xmin);
      Double_t thisxmax = std::min(*hi, xmax);

      // add eff name, as it specifies the boundaries...
      TString range = TString::Format("R%d_%s_%s", i,x.GetName(),eff.GetName());

      // Add original rangeName if there is one
      if (rangeName) { 
            range.Append( "_" );
            range.Append( RooNameReg::str(rangeName) );
      }

      RooNameSet ns(iset);
      range.Append("_I_");
      range.Append(ns._nameList);

      // Create a new name for the range
      // check if already exists and matches..
      if (!x.hasRange(range)) {
        x.setRange(range, thisxmin, thisxmax);
      } else {
        assert( x.getMin(range)==thisxmin);
        assert( x.getMax(range)==thisxmax);
        cout << "re-using existing range " << range << endl;
      }
      
      intList.add(*model.createIntegral(iset, range));

      // create RooAbsReal for (average) efficiency in this range
      RooCustomizer* customizer = new RooCustomizer(eff, (range+"_customizer").Data());
      RooRealVar* cv = static_cast<RooRealVar*>(x.clone( TString(x.GetName()) + "_" + range) );
      cv->setVal((thisxmin + thisxmax) / 2.);
      cv->setConstant(true);
      customizer->replaceArg(x, *cv);
      // leak cv??
      effList.add( *customizer->build(kFALSE) );
      _customizers.push_back(customizer);
   }
   // TODO: create unique name
   TString iName( parent.GetName() );
   iName += "_integral";
   _I = new RooAddition(iName, iName, effList, intList, kTRUE);
}

//_____________________________________________________________________________
RooEffResModel::RooEffResModel(const char *name, const char *title, RooResolutionModel& model, RooAbsReal& eff) 
  : RooResolutionModel(name,title,model.convVar())
  , _model("!model","Original resolution model",this,model)
  , _eff("!efficiency","efficiency of convvar", this,eff)
  , _cacheMgr(this, 10)
{  
  // assert that efficiency is a function of convVar, and there are no overlaps...
}

//_____________________________________________________________________________
RooEffResModel::RooEffResModel(const RooEffResModel& other, const char* name) 
  : RooResolutionModel(other,name)
  , _model("!model",this,other._model)
  , _eff("!efficiency",this,other._eff)
  , _cacheMgr(other._cacheMgr,this)
{
}

//_____________________________________________________________________________
RooEffResModel::~RooEffResModel()
{
  // Destructor
}

//_____________________________________________________________________________
Int_t RooEffResModel::basisCode(const char* name) const 
{ 
   return model().basisCode(name);
} 

//_____________________________________________________________________________
Double_t RooEffResModel::evaluate() const 
{  
    Double_t mod  = model().getVal();
    Double_t eps  = eff().getVal();
    return eps*mod;
}


//_____________________________________________________________________________
Bool_t RooEffResModel::forceAnalyticalInt(const RooAbsArg& /*dep*/) const
{
  // Force RooRealIntegral to offer all observables for internal integration
   return kTRUE ;
}

//_____________________________________________________________________________
Int_t RooEffResModel::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
    if (_forceNumInt) return 0;
    analVars.add(allVars);
    getCache(&analVars,RooNameReg::ptr(rangeName));
    return 1 + _cacheMgr.lastIndex();
}

//_____________________________________________________________________________
Double_t RooEffResModel::analyticalIntegral(Int_t code, const char* rangeName) const 
{
   assert(code > 0);
   CacheElem* cache = static_cast<CacheElem*>(_cacheMgr.getObjByIndex(code - 1));
   if (!cache) {
      std::auto_ptr<RooArgSet> vars(getParameters(RooArgSet()));
      std::auto_ptr<RooArgSet> iset( _cacheMgr.nameSet1ByIndex(code - 1)->select(*vars));
      cache = getCache(iset.get(), RooNameReg::ptr(rangeName) );
      assert(cache!=0);
   }
   return cache->getVal();
}

//_____________________________________________________________________________
Int_t RooEffResModel::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  return 0 ; // For now... problem is that RooGenConv assumes it can just add resolution & physics for conv var...
}

//_____________________________________________________________________________
RooEffResModel::CacheElem* RooEffResModel::getCache(const RooArgSet *iset, const TNamed *rangeName) const 
{
   Int_t sterileIndex(-1);
   CacheElem* cache = (CacheElem*) _cacheMgr.getObj(iset, &sterileIndex, rangeName);
   if (cache) return cache;
   _cacheMgr.setObj(iset, new CacheElem( *this,  *iset,  rangeName), rangeName);
   return getCache(iset, rangeName );
}

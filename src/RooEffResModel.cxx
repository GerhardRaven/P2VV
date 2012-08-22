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

#include <memory>

#include "RooFit.h"
#include "Riostream.h"
#include "RooEffResModel.h"
#include "RooRealConstant.h"
#include "RooCustomizer.h"
#include "RooAddition.h"
using namespace std;

ClassImp(RooEffResModel) 
;

//_____________________________________________________________________________
RooEffResModel::CacheElem::~CacheElem()
{
   delete _I;
   for (std::vector<RooCustomizer*>::const_iterator it = _customizers.begin(),
           end = _customizers.end(); it != end; ++it) delete *it;
}

//_____________________________________________________________________________
RooArgList RooEffResModel::CacheElem::containedArgs(Action) 
{
   // Return list of all RooAbsArgs in cache element
   return RooArgList(*_I);
}

//_____________________________________________________________________________
RooEffResModel::CacheElem::CacheElem(const RooEffResModel& parent, const RooArgSet& iset,
                                     const TNamed* rangeName)
   : _I(0)
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

      // move range name generation (with the exception of the R%d prefix)
      // out of the loop...
      // add eff name, as it specifies the boundaries...
      TString trange = TString::Format("R%d_%s_%s", i,x.GetName(),eff.GetName());

      // Add original rangeName if there is one
      if (rangeName) { 
            trange.Append( "_" );
            trange.Append( RooNameReg::str(rangeName) );
      }

      trange.Append("_I_");
      RooNameSet ns(iset);
      trange.Append(ns._nameList);


      const char *range = RooNameReg::str( RooNameReg::ptr( trange.Data() ) );


      // Create a new name for the range
      // check if already exists and matches..
      if (!x.hasRange(range)) {
        x.setRange(range, thisxmin, thisxmax);
      }
      assert( x.getMin(range)==thisxmin);
      assert( x.getMax(range)==thisxmax);
      
      intList.add(*model.createIntegral(iset, range));

      // create RooAbsReal for (average) efficiency in this range
      RooCustomizer* customizer = new RooCustomizer(eff, (trange+"_customizer").Data());
      RooRealVar* cv = static_cast<RooRealVar*>(x.clone( TString(x.GetName()) + "_" + range) );
      cv->setVal((thisxmin + thisxmax) / 2.);
      cv->setConstant(true);
      customizer->replaceArg(x, *cv);
      RooAbsArg *ceff = customizer->build(kFALSE);
      ceff->addOwnedComponents(*cv);
      effList.add( *ceff );
      _customizers.push_back(customizer);
   }
   TString iName = TString::Format("%s_I_%s", parent.GetName(),x.GetName());
   _I = new RooAddition(iName, iName, effList, intList, kTRUE);
}

//_____________________________________________________________________________
RooEffResModel::RooEffResModel(const char *name, const char *title, RooResolutionModel& model, RooAbsReal& eff) 
   : RooResolutionModel(name,title,model.convVar())
   , _observables("observables", "observables", this)
   , _model("!model","Original resolution model",this,model)
   , _eff("!efficiency","efficiency of convvar", this,eff)
   , _cacheMgr(this, 10)
{  
   // assert that efficiency is a function of convVar, and there are no overlaps...
   _observables.add(model.convVar());
}

//_____________________________________________________________________________
RooEffResModel::RooEffResModel(const RooEffResModel& other, const char* name) 
  : RooResolutionModel(other,name)
  , _observables("observables", this, other._observables)
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
RooEffResModel* 
RooEffResModel::convolution(RooFormulaVar* inBasis, RooAbsArg* owner) const
{
  // Instantiate a clone of this resolution model representing a convolution with given
  // basis function. The owners object name is incorporated in the clones name
  // to avoid multiple convolution objects with the same name in complex PDF structures.
  // 
  // Note: The 'inBasis' formula expression must be a RooFormulaVar that encodes the formula
  // in the title of the object and this expression must be an exact match against the
  // implemented basis function strings (see derived class implementation of method basisCode()
  // for those strings

  // Check that primary variable of basis functions is our convolution variable  
  if (inBasis->getParameter(0) != x.absArg()) {
    coutE(InputArguments) << "RooEffResModel::convolution(" << GetName() << "," << this
              << ") convolution parameter of basis function and PDF don't match" << endl
              << "basis->findServer(0) = " << inBasis->findServer(0) << endl
              << "x.absArg()           = " << x.absArg() << endl ;
    return 0 ;
  }

  if (basisCode(inBasis->GetTitle())==0) {
    coutE(InputArguments) << "RooEffResModel::convolution(" << GetName() << "," << this
              << ") basis function '" << inBasis->GetTitle() << "' is not supported." << endl ;
    return 0 ;
  }

  TString newName(GetName()) ;
  newName.Append("_conv_") ;
  newName.Append(inBasis->GetName()) ;
  newName.Append("_[") ;
  newName.Append(owner->GetName()) ;
  newName.Append("]") ;

  RooResolutionModel *conv = model().convolution(inBasis,owner);

  TString newTitle(conv->GetTitle()) ;
  newTitle.Append(" convoluted with basis function ") ;
  newTitle.Append(inBasis->GetName()) ;
  conv->SetTitle(newTitle.Data()) ;

  RooEffResModel *effConv = new RooEffResModel(newName,newTitle,*conv,eff());
  effConv->addOwnedComponents(*conv);
  effConv->changeBasis(inBasis) ;
  return effConv ;
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
    // TODO: replace this by the discretized version, i.e. replace convVar by customized middle of bin...
    //       this in order to ensure evaluate & analyticalIntegral are consistent (in case eff is not discretized!!!)
    Double_t eps  = eff().getVal(); 
    return eps * mod;
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
    getCache(&analVars, RooNameReg::ptr(rangeName));
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

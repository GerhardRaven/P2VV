/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id: RooEffHistProd.cxx 25184 2008-08-20 13:59:55Z wouter $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, NIKHEF
 *   GR, Gerhard Raven, NIKHEF/VU                                            *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/


/////////////////////////////////////////////////////////////////////////////////////
// BEGIN_HTML
// The class RooEffHistProd implements the product of a PDF with an efficiency function.
// The normalization integral of the product is calculated numerically, but the
// event generation is handled by a specialized generator context that implements
// the event generation in a more efficient for cases where the PDF has an internal
// generator that is smarter than accept reject. 
// END_HTML
//

#include "RooFit.h"
#include "RooEffHistProd.h"
#include "RooEffGenContext.h"
#include "RooNameReg.h"
#include "RooRealVar.h"
#include <memory>
#include <algorithm>

ClassImp(RooEffHistProd)
  ;

//_____________________________________________________________________________
RooEffHistProd::CacheElem::CacheElem(const RooEffHistProd* parent,RooArgSet& analVars,const char *rangeName)
    : I(0),xmin(0),xmax(0)
{
  // create a function object for the corresponding integral of the underlying PDF.
  // i.e. if we have  eps(x)f(x,y)  and we get (x,y) as allVars, 
  // construct I(xmin,xmax) = int_xmin^xmax dx int dy f(x,y)
  // so that later we can do sum_i eps( (x_i+x_i+1)/2 ) * I(x_i,x_i+1)
  RooArgSet *x_ = parent->eff()->getObservables(&analVars); // the subset of analVars on which _eff depends
  const char *myRange = rangeName;
  assert(x_->getSize()<2); // for now, we only do 1D efficiency histograms...
  if (x_->getSize()==1) {
      assert(rangeName==0); // deal with ranges later -- this is a bit non-trivial.... need to clone the original...
      assert( *x_->first() == parent->x() );
      // TODO: add original rangeName in here!
      const char *name = parent->makeFPName(parent->GetName(),analVars,"_I_Range_min");
      xmin = new RooRealVar(name,name,parent->x().getMin(rangeName),parent->x().getMin(rangeName),parent->x().getMax(rangeName));
      name = parent->makeFPName(parent->GetName(),analVars,"_I_Range_max");
      xmax = new RooRealVar(name,name,parent->x().getMax(rangeName),parent->x().getMin(rangeName),parent->x().getMax(rangeName));
      myRange = parent->makeFPName(parent->GetName(),analVars,"_I_Range");
      if (parent->x().hasRange(myRange)) {
        cout << "RooEffHistProd("<<parent->GetName() << ")::CacheElem  range " << myRange << " already exists!!!" << endl;
      }
      parent->x().setRange(myRange,*xmin,*xmax);
      cout << "RooEffHistProd("<<parent->GetName() << ")::CacheElem created range " << myRange << " from " << xmin->GetName() << " to " << xmax->GetName() << endl;
  }
  I = parent->pdf()->createIntegral(analVars,analVars,myRange);
}
//_____________________________________________________________________________
RooEffHistProd::CacheElem::~CacheElem() 
{
}

//_____________________________________________________________________________
RooArgList RooEffHistProd::CacheElem::containedArgs(Action) 
{
  // Return list of all RooAbsArgs in cache element
  return RooArgList(*I,*xmin,*xmax); 
}

//_____________________________________________________________________________
RooEffHistProd::RooEffHistProd(const char *name, const char *title, 
			       RooAbsPdf& inPdf, RooAbsReal& inEff) :
  RooAbsPdf(name,title),
  _pdf("pdf","pre-efficiency pdf", this,inPdf),
  _eff("eff","efficiency function",this,inEff),
  _observables("obs","observables in efficiency function",this),
  _cacheMgr(this,10)
{  
  // Constructor of a a production of p.d.f inPdf with efficiency
  // function inEff.

  // to figure out what the observable is, we look at the overlap of
  // the variables of efficiency function and pdf.
  
  std::auto_ptr<RooArgSet> pdfpars( inPdf.getVariables() );
  std::auto_ptr<RooArgSet> effpars( inEff.getVariables() );
  std::auto_ptr<TIterator> iter(   effpars->createIterator() );

  RooAbsArg *effelem(0);
  while( (effelem = (RooAbsArg*)iter->Next()) ) {
    if( pdfpars->find( effelem->GetName() ) )  _observables.add( *effelem ) ;
  }

  if(_observables.getSize() == 0 ) {
    throw std::string("WARNING: RooEffHistProd: PDF and Efficiency function factorise. Please use RooProd") ;
  } else if(_observables.getSize()>1) {
    throw std::string("WARNING: RooEffHistProd not yet implemented for more than 1D efficiency" ) ;
  }
  
  // an interesting hack. need to discuss with wouter. one idea: let
  // every function that is 'binned' (discrete,quantized,..) add a 
  // special binning object to its dependents.
  std::list<Double_t>* binboundaries = inEff.plotSamplingHint(x(),x().getMin(),x().getMax()) ;
  if(binboundaries) {
    for( std::list<Double_t>::const_iterator it = binboundaries->begin() ;
	 it != binboundaries->end() ; ++it ) {
      double x1 = *it++, x2 = *it ;
      //std::cout << "binboundaries: " << x1 << "," << x2 << std::endl ;
      _binboundaries.push_back( 0.5*(x1+x2) ) ;
    }
    //std::copy( _binboundaries.begin(), _binboundaries.end(), ostream_iterator<Double_t>(std::cout,", ") );
    //std::cout << std::endl;
  }
  delete binboundaries ;
}

//_____________________________________________________________________________
RooEffHistProd::RooEffHistProd(const RooEffHistProd& other, const char* name) : 
  RooAbsPdf(other, name),
  _pdf("pdf",this,other._pdf),
  _eff("acc",this,other._eff),
  _observables("obs",this,other._observables), 
  _cacheMgr(other._cacheMgr,this),
  _binboundaries(other._binboundaries)
{
  // Copy constructor
}

//_____________________________________________________________________________
RooEffHistProd::~RooEffHistProd() 
{
  // Destructor
}

//_____________________________________________________________________________
Double_t RooEffHistProd::evaluate() const
{
  // Calculate and return 'raw' unnormalized value of p.d.f
  // cout << "_normSet = " << (_normSet?*_normSet:RooArgSet()) << endl ;
  return eff()->getVal() * pdf()->getVal(_normSet);
}

//_____________________________________________________________________________
RooAbsGenContext* RooEffHistProd::genContext(const RooArgSet &vars, const RooDataSet *prototype,
                                            const RooArgSet* auxProto, Bool_t verbose) const
{

  return RooAbsPdf::genContext(vars,prototype,auxProto,verbose);
  // Return specialized generator context for RooEffHistProds that implements generation
  // in a more efficient way than can be done for generic correlated products
  assert(pdf()!=0);
  assert(eff()!=0);
  return new RooEffGenContext(*this,*pdf(),*eff(),vars,prototype,auxProto,verbose) ;
}

//_____________________________________________________________________________
const char* RooEffHistProd::makeFPName(const char *prefix,const RooArgSet& terms, const char *postfix) const
{
  static TString pname;
  pname = prefix;
  if (prefix) pname.Append("_");
  std::auto_ptr<TIterator> i( terms.createIterator() );
  RooAbsArg *arg;
  Bool_t first(kTRUE);
  while((arg=(RooAbsArg*)i->Next())) {
    if (first) { first=kFALSE;}
    else pname.Append("_X_");
    pname.Append(arg->GetName());
  }
  pname.Append(postfix);
  return pname.Data();
}

//_____________________________________________________________________________
Int_t RooEffHistProd::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
  //return 0 ;

  assert(rangeName==0);
  //cout << "RooEffHistProd("<<GetName()<<")::getAnalyticalIntegral("; allVars.printValue(cout); cout << ","<< ( rangeName? rangeName : "<none>" ) <<")"<<endl;
  if (_forceNumInt) return 0;
  if (allVars.getSize()==0) return 0;
  analVars.add(allVars) ;

  Int_t sterileIndex(-1);
  CacheElem* cache = (CacheElem*) _cacheMgr.getObj(&analVars,&sterileIndex,RooNameReg::ptr(rangeName));
  if (cache!=0) {
    //cout << "RooEffHistProd("<<GetName()<<")::getAnalyticalIntegral: already have cache, returning 1+"<< _cacheMgr.lastIndex() << endl;
    return 1+_cacheMgr.lastIndex();
  }

  Int_t code = _cacheMgr.setObj(&analVars,new CacheElem(this,analVars,rangeName),RooNameReg::ptr(rangeName));
  //cout << "RooEffHistProd("<<GetName()<<")::getAnalyticalIntegral: returning 1+" << code << " for integral over " ;
  //analVars.printValue(cout); cout << " @ " << &analVars << " in range " << (rangeName?rangeName:"<none>") << endl;
  return 1+code ;
}

//_____________________________________________________________________________
Double_t RooEffHistProd::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  assert(rangeName==0);
  assert(code>0);
  CacheElem* cache = (CacheElem*) _cacheMgr.getObjByIndex(code-1);
  RooArgSet *vars = getParameters(RooArgSet());
  RooArgSet *nset = _cacheMgr.nameSet1ByIndex(code-1)->select(*vars);
  assert(cache!=0); // if this triggers, add code to repopulate the sterilized cache entry...

  // TODO: verify rangeName is consistent with the one in getAnalyticalIntegral...
  Double_t xmin = x().getMin(rangeName), xmax = x().getMax(rangeName);

  // make sure the range does is contained within the binboundaries...
  assert(_binboundaries.size()>1);
  assert(xmin<=xmax);
  assert(xmin>=_binboundaries.front());
  assert(_binboundaries.back()>=xmax);

  if (cache->xmin==0 && cache->xmax==0)  return eff()->getVal()*cache->getVal(nset); // no integral over efficiency dependant...

  double xorig = x().getVal();
  double norm1(0);
  double sum(0);
  for(BinBoundaries::const_iterator i = _binboundaries.begin(), end = _binboundaries.end(); i+1 != end; ++i ) {
    if (*(i+1)<xmin) continue; // not there yet...
    if (*i    >xmax) break;    // past the requested interval...
    Double_t thisxmin = std::max(*i,     xmin);
    Double_t thisxmax = std::min(*(i+1), xmax);
    if (thisxmin>=thisxmax) continue;
    //cout << "RooEffHistProd: setting " << x().GetName() << endl;
    x().setVal( 0.5*( thisxmin + thisxmax) ) ; // get the efficiency for this bin
    sum += eff()->getVal() * cache->getVal(thisxmin,thisxmax,nset); // *(thisxmax-thisxmin);
    //cout << "RooEffHistProd("<<GetName()<<")::analyticalIntegral: integral in bin from " << cache->xmin->getVal() << "," << cache->xmax->getVal() << " =  " << eff()->getVal(nset) << " * " << cache->getVal(thisxmin,thisxmax,nset) << endl;
  }
  //cout << "RooEffHistProd("<<GetName()<<")::analyticalIntegral: integral =  " << sum << endl;// " norm =  " << norm1 << endl;
  x().setVal(xorig);
  return  sum;
}

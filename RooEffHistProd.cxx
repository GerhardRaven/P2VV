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

ClassImp(RooEffHistProd)
  ;



//_____________________________________________________________________________
RooEffHistProd::RooEffHistProd(const char *name, const char *title, 
			       RooAbsPdf& inPdf, RooAbsReal& inEff) :
  RooAbsPdf(name,title),
  _pdf("pdf","pre-efficiency pdf", this,inPdf),
  _eff("eff","efficiency function",this,inEff),
  _observables("obs","observables in efficiency function",this),
  _nset(0),
  _fixedNset(0)
{  
  // Constructor of a a production of p.d.f inPdf with efficiency
  // function inEff.

  // to figure out what the observable is, we look at the overlap of
  // the variables of efficiency function and pdf.
  
  RooArgSet* pdfpars = inPdf.getVariables() ;
  RooArgSet* effpars = inEff.getVariables() ;
  TIterator* iter =  effpars->createIterator() ;

  RooAbsArg *effelem,*pdfelem ;
  //std::vector<RooAbsArg*> overlap ;
  while( (effelem = (RooAbsArg*)iter->Next()) ) {
    if( (pdfelem = pdfpars->find( effelem->GetName() ) ) )
      _observables.add( *effelem ) ;
    //overlap.push_back( effelem ) ;
  }
  delete pdfpars ;
  delete effpars ;
  delete iter ;

  if(_observables.getSize() == 0 ) {
    throw std::string("WARNING: RooEffHistProd: PDF and Efficiency function factorise. Please use RooProd") ;
  } else if(_observables.getSize()>1) {
    throw std::string("WARNING: RooEffHistProd not yet implemented for more than 1D efficiency" ) ;
  }
  
  //RooRealVar* rv = dynamic_cast< RooRealVar* >( overlap.front() ) ;
  //_x = RooRealProxy ("x",rv->GetName(),this,*rv ) ;
  //std::cout << "proxy: " << _x << std::endl ;
  //std::cout << "x():" << x() << std::endl ;

  // an interesting hack. need to discuss with wouter. one idea: let
  // every function that is 'binned' add a special binning object to
  // its dependents.
  std::list<Double_t>* binboundaries = inEff.plotSamplingHint(x(),x().getMin(),x().getMax()) ;
  if(binboundaries) {
    for( std::list<Double_t>::const_iterator it = binboundaries->begin() ;
	 it != binboundaries->end() ; ++it ) {
      double x1 = *it ;
      std::cout << "x1 binboundary: " << x1 << std::endl ;
      ++it ;
      double x2 = *it ;
      std::cout << "x2 binboundary: " << x2 << std::endl ;
      _binboundaries.push_back( 0.5*(x1+x2) ) ;
    }
    for( std::vector<Double_t>::const_iterator it = _binboundaries.begin() ;
	 it != _binboundaries.end() ; ++it)
      std::cout << *it << ", " << std::flush ;
    std::cout << std::endl ;
  }
  delete binboundaries ;
}




//_____________________________________________________________________________
RooEffHistProd::RooEffHistProd(const RooEffHistProd& other, const char* name) : 
  RooAbsPdf(other, name),
  _pdf("pdf",this,other._pdf),
  _eff("acc",this,other._eff),
  _observables("obs",this,other._observables), 
  _binboundaries(other._binboundaries),
  _nset(0),
  _fixedNset(0) 
{
  // Copy constructor
}




//_____________________________________________________________________________
RooEffHistProd::~RooEffHistProd() 
{
  // Destructor
}



//_____________________________________________________________________________
Double_t RooEffHistProd::getVal(const RooArgSet* set) const 
{  
  // Return p.d.f. value normalized over given set of observables

  _nset = _fixedNset ? _fixedNset : set ;
  return RooAbsPdf::getVal(set) ;
}




//_____________________________________________________________________________
Double_t RooEffHistProd::evaluate() const
{
  // Calculate and return 'raw' unnormalized value of p.d.f

  return eff()->getVal() * pdf()->getVal(_nset);
  //return eff()->evaluate() * pdf()->evaluate() ;
}


//_____________________________________________________________________________
RooAbsGenContext* RooEffHistProd::genContext(const RooArgSet &vars, const RooDataSet *prototype,
                                            const RooArgSet* auxProto, Bool_t verbose) const
{
  // Return specialized generator context for RooEffHistProds that implements generation
  // in a more efficient way than can be done for generic correlated products

  assert(pdf()!=0);
  assert(eff()!=0);
  return new RooEffGenContext(*this,*pdf(),*eff(),vars,prototype,auxProto,verbose) ;
}


//_____________________________________________________________________________
Int_t RooEffHistProd::getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars, 
					  const RooArgSet* normSet, const char* rangeName) const 
{
  // WouterV proposes to o this with getIntegral. I just don't quite
  // understand how the caching works.
  Int_t rc = pdf()->getAnalyticalIntegralWN(allVars,analVars,normSet,rangeName) ;
  if( rc != 0 ) _integralmap[rc] = analVars ;
  return rc ;
}


//_____________________________________________________________________________
Double_t RooEffHistProd::analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName) const 
{
  // Return value of integral identified by code, which should be a return value of getAnalyticalIntegralWN,
  // Code zero is always handled and signifies no integration (return value is normalized p.d.f. value)

  // Return analytical integral defined by given scenario code
  double rc = 0 ;

  // No integration scenario
  if (code==0) {
    rc =  getVal(normSet) ;
  } else {
    std::map<Int_t, RooArgList>::const_iterator mit =_integralmap.find(code) ;
    if( mit != _integralmap.end() ) {
      std::cout << "NOW COMPUTING INTEGRAL" << std::endl ;
      RooRealVar& effobs = x() ;
      //std::cout << "Value of x = " << std::flush << effobs.getVal() << std::endl;
      //RooArgList observables ;
      //observables.add(*(eff()->getVariables())) ;
      //RooRealVar& effobs = *(dynamic_cast<RooRealVar*>(&(observables[0]))) ;
      //it->second.Print() ;
      Double_t xmin = effobs.getMin(rangeName) ;
      Double_t xmax = effobs.getMax(rangeName) ;

      //std::cout << "xmin,xmax: " << xmin << " " << xmax << std::endl ;
      // We could use upper and lowerbound to find the right
      // bins. However, in most cases we anyway fit the entire
      // range. So, it's a bit a waste of time.      
      //BinBoundaries::const_iterator itbegin = std::lower_bound(_binboundaries.begin(),
      //							       _binboundaries.end(),xmin) ;
      //BinBoundaries::const_iterator itend = std::upper_bound(_binboundaries.begin(),
      //							     _binboundaries.end(),xmax) ;
      BinBoundaries::const_iterator itbegin = _binboundaries.begin() ;
      BinBoundaries::const_iterator itend = _binboundaries.end() ;

      // We could use upper and lowerbound to find the right
      // bins. However, in most cases we anyway fit the entire
      // range. So, it's a bit a waste of time.      
      double effsum(0), sum(0) ;
      BinBoundaries::const_iterator itnext = itbegin ;
      BinBoundaries::const_iterator it = itnext++ ;
      //std::cout << "binboundaries " << _binboundaries.size() << std::endl ;
      for(; itnext != itend ; ++itnext) {
	// compute the interval ranges
	Double_t thisxmin = std::max(*it,xmin) ;
	std::cout << "thisxmin " << thisxmin << std::endl ;
	Double_t thisxmax = std::min(*itnext,xmax) ;
	std::cout << "thisxmax " << thisxmax << std::endl ;
	// don't remove this: the if statement is needed, unless we
	// use upper/lowerbound to compute the correct iterartors aboce.
	if( thisxmin < thisxmax ) { 
	  // get the value of the efficiency for this bin
	  effobs.setVal( 0.5*( thisxmin + thisxmax) ) ;
	  Double_t thiseff = eff()->getVal() ;
	  std::cout << "thiseff " << thiseff << std::endl ;
	  // integrate the pdf
	  effobs.setRange(rangeName,thisxmin,thisxmax) ;
	  std::cout << "code " << code << std::endl ;
	  std::cout << "normSet " << normSet << std::endl ;
	  if (rangeName ==0){
	    std::cout << "shit " << std::endl ;
	  }
	  //std::cout << "rangeName " << rangeName << std::endl ;
	  
	  Double_t thisintegral = pdf()->analyticalIntegralWN(code,normSet,rangeName) ;
	  std::cout << "thisintegral " << thisintegral << std::endl ;
	  effsum += thiseff * thisintegral ;
	  std::cout << "effsum " << effsum << std::endl ;
	  sum += thisintegral ;
	  std::cout << "sum " << sum << std::endl ;
	  it = itnext ;
	}
      }
      // set the interval back to what it was
      effobs.setRange(rangeName,xmin,xmax) ;
      rc = sum > 0 ? effsum / sum : 0 ;
      //std::cout << "efsum, sum: " << effsum << " " << sum << std::endl ;
    } else {
      rc = getVal(normSet) ;
    }
  }
  return rc ;
}

/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOIPATIA2
#define ROOIPATIA2

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooAbsReal;
 
class RooIpatia2 : public RooAbsPdf {
public:
  RooIpatia2() {} ; 
  RooIpatia2(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _l,
	      RooAbsReal& _zeta,
	      RooAbsReal& _fb,
	      RooAbsReal& _sigma,
	      RooAbsReal& _mu,
	      RooAbsReal& _a,
	      RooAbsReal& _n,
	      RooAbsReal& _a2,
	      RooAbsReal& _n2);
  RooIpatia2(const RooIpatia2& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooIpatia2(*this,newname); }
  inline virtual ~RooIpatia2() { }

protected:

  RooRealProxy x ;
  RooRealProxy l ;
  RooRealProxy zeta ;
  RooRealProxy fb ;
  RooRealProxy sigma ;
  RooRealProxy mu ;
  RooRealProxy a ;
  RooRealProxy n ;
  RooRealProxy a2 ;
  RooRealProxy n2 ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooIpatia2,1) // Your description goes here...
};
 
#endif

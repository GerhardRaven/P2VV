/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOAVEFFCONSTRAINT
#define ROOAVEFFCONSTRAINT

#include <memory>
#include <vector>

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

#include "RooEffHistProd.h"

class RooEffHistProd;

class RooAvEffConstraint : public RooAbsPdf {
public:
   RooAvEffConstraint() {} ; 
   RooAvEffConstraint(const char *name, const char *title, 
                      RooRealVar& observable, RooEffHistProd& effProd,
                      RooRealVar& mean, RooRealVar& sigma);
   RooAvEffConstraint(const RooAvEffConstraint& other, const char* name=0) ;

   virtual TObject* clone(const char* newname) const
   {
      return new RooAvEffConstraint(*this,newname);
   }
   
   virtual Bool_t forceAnalyticalInt(const RooAbsArg& /*dep*/) const;

   virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
   virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

   virtual ~RooAvEffConstraint();
   
protected:
   
   Double_t evaluate() const;
   
private:

   RooRealProxy _observable;
   RooRealProxy _shape;
   RooRealProxy _mean;
   RooRealProxy _sigma;
   mutable RooRealProxy* _integral;
   
   ClassDef(RooAvEffConstraint, 1) // Your description goes here...
};

#endif

/*****************************************************************************
 * Projects: RooFit, P2VV                                                    *
 *                                                                           *
 * Author: Roel Aaij, roel.aaij@nikhef.nl                                    * 
 *****************************************************************************/

#ifndef ROOAVEFFCONSTRAINT
#define ROOAVEFFCONSTRAINT

#include <memory>
#include <vector>

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

#include "P2VV/RooEffResModel.h"

class RooEffResModel;
class RooCustomizer;

class RooAvEffConstraint : public RooAbsPdf {
public:

   RooAvEffConstraint();

   RooAvEffConstraint(const char *name, const char *title, 
                      RooAbsPdf& pdf, RooEffResModel& model,
                      RooRealVar& mean, RooRealVar& sigma);
   RooAvEffConstraint(const RooAvEffConstraint& other, const char* name=0) ;

   virtual TObject* clone(const char* newname) const
   {
      return new RooAvEffConstraint(*this,newname);
   }
   
   virtual Bool_t forceAnalyticalInt(const RooAbsArg& /*dep*/) const;

   virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
   virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

   virtual const RooAbsReal* efficiency() const;

   virtual void setPdf(RooAbsPdf& pdf);

   virtual ~RooAvEffConstraint();

protected:
   
   Double_t evaluate() const;
   
private:

   void initialize();

   RooRealProxy _mean;
   RooRealProxy _sigma;
   mutable RooRealProxy* _integral; //!
   RooListProxy* _integrals; //!
   RooListProxy* _efficiencies; //!
   RooRealProxy _model;
   RooRealProxy* _pdf;
   
   ClassDef(RooAvEffConstraint, 1) // Your description goes here...
};

#endif

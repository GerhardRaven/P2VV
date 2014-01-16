/*****************************************************************************
 * Projects: RooFit, P2VV                                                    *
 *                                                                           *
 * Author: Roel Aaij, roel.aaij@nikhef.nl                                    * 
 *****************************************************************************/

#ifndef ROOEFFCONSTRAINT
#define ROOEFFCONSTRAINT

#include <memory>
#include <vector>

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooEffResModel;

class RooEffConstraint : public RooAbsPdf {
public:

   RooEffConstraint();

   RooEffConstraint(const char *name, const char *title, const RooArgList& eps_a,
                    const RooArgList& eps_b, const std::vector<double>& N_a,
                    const std::vector<double>& N_b, const std::vector<double>& N_ab = std::vector<double>());

   RooEffConstraint(const RooEffConstraint& other, const char* name=0);

   virtual TObject* clone(const char* newname) const
   {
      return new RooEffConstraint(*this,newname);
   }
   
   // virtual Bool_t forceAnalyticalInt(const RooAbsArg& /*dep*/) const;
   // virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
   // virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

   virtual Double_t getLogVal(const RooArgSet* nset = 0) const;

   virtual ~RooEffConstraint();

protected:
   
   virtual Double_t evaluate() const;
   
private:

   void initialize();

   RooListProxy _eps_a;
   RooListProxy _eps_b;
   std::vector<double> _N_a;
   std::vector<double> _N_b;
   std::vector<double> _N_ab;

   ClassDef(RooEffConstraint, 1) // RooEffConstraint
};

#endif

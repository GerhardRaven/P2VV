#ifndef ROO_MULTIEFFICIENCY
#define ROO_MULTIEFFICIENCY

#include <RooAbsPdf.h>
#include <RooListProxy.h>
#include <TString.h>
#include <TList.h>

class RooArgList ;

class RooMultiEfficiency : public RooAbsPdf {
public:
  // Constructors, assignment etc
  inline RooMultiEfficiency() { 
    // Default constructor
  }
  RooMultiEfficiency(const char *name, const char *title, const RooArgList& efficiencies, const RooArgList& categories,
                const TList& sigCatNames);
  RooMultiEfficiency(const RooMultiEfficiency& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooMultiEfficiency(*this, newname); }
  virtual ~RooMultiEfficiency();

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const ;

protected:

  // Function evaluation
  virtual Double_t evaluate() const ;
  RooListProxy _categories ;   // Accept/reject categories
  RooListProxy _efficiencies ; // Efficiency modeling functions
  TList _sigCatNames ;         // Name of accept state of accept/reject categories

  ClassDef(RooMultiEfficiency,1) // Generic PDF defined by string expression and list of variables
};

#endif

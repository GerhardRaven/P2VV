/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooRelBreitWignerPdf.rdl,v 1.1 2006/04/08 21:59:57 bevan Exp $
 * Authors:                                                                  *
 * Katherine George - University of Liverpool                                * 
 * Adrian Bevan QMUL (ported from RptBaseTools to Roo - 2009)               *
 *****************************************************************************/
#ifndef AFIT_RELBREITWIGNERPDF
#define AFIT_RELBREITWIGNERPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooComplex.h"

class RooRealVar;

class RooRelBreitWigner : public RooAbsPdf {
public:
  RooRelBreitWigner() {} ;
  RooRelBreitWigner(const char *name, const char *title,
	      RooAbsReal& _x, RooAbsReal& _mean, RooAbsReal& _width, RooAbsReal& spin, RooAbsReal& radius,
              RooAbsReal& massa, RooAbsReal& massb);
  RooRelBreitWigner(const RooRelBreitWigner& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooRelBreitWigner(*this,newname); }
  inline virtual ~RooRelBreitWigner() { }

  RooComplex amplitude() const ;

protected:

  Double_t gamma() const;
  Double_t KFunction(Double_t X) const;
  Double_t FFunction(Double_t X) const;

/*
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
*/

protected:

  RooRealProxy x ;
  RooRealProxy mean ;
  RooRealProxy width ;
  RooRealProxy spin;
  RooRealProxy radius;
  RooRealProxy massa;
  RooRealProxy massb;

  Double_t evaluate() const ;

//   void initGenerator();
//   Int_t generateDependents();

private:

  ClassDef(RooRelBreitWigner,1) // Relativistic Breit Wigner PDF
};

#endif

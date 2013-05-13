/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   Gerhard Raven
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_CUBICSPLINEFUN
#define ROO_CUBICSPLINEFUN

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooCubicSplineKnot.h"
#include "RooCubicSplineGaussModel.h"

class RooRealVar;
class RooArgList ;
class TH1;

class RooCubicSplineFun : public RooAbsReal {
public:
  RooCubicSplineFun() ;
  // smooth = 0: no smoothing. As smooth becomes larger, the result will converge towards a straight line
  RooCubicSplineFun(const char* name, const char* title, RooRealVar& x,
                    const std::vector<double>& knots,
                    const std::vector<double>& values,
                    const std::vector<double>& errors = std::vector<double>(),
                    double smooth = 0, bool constCoeffs = true);
  RooCubicSplineFun(const char* name, const char* title, RooRealVar& x, const TH1* hist,
                    double smooth = 0, bool constCoeffs = true);
  RooCubicSplineFun(const char *name, const char *title, RooRealVar& x,
                    const char *knotBinningName, const RooArgList& coefList) ;
  RooCubicSplineFun(const char* name, const char* title, RooRealVar& x,
                     const std::vector<double>& knots, const RooArgList& coefList);

  ~RooCubicSplineFun() ;

  RooCubicSplineFun(const RooCubicSplineFun& other, const char* name = 0);
  TObject* clone(const char* newname) const { return new RooCubicSplineFun(*this, newname); }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName) const;

  // for use in RooCubicSplineGaussModel...
  RooComplex gaussIntegral(int i, const RooCubicSplineGaussModel::M_n& dM,
                           const RooCubicSplineGaussModel::K_n& K,
                           double offset, double* sc) const ;
  unsigned knotSize() const { return _aux->size(); }
  double u(int i) const { return _aux->u(i); }
  const std::vector<double>& knots() const { return _aux->knots(); }

  const RooArgList& coefficients() const { return _coefList; }

private:

  RooRealProxy _x;
  RooListProxy _coefList ;

  const RooCubicSplineKnot *_aux; // do not persist! (but do persist the binningName used for _x!!!

  void init(const char* name, const std::vector<double>& knots, const std::vector<double>& heights,
            const std::vector<double>& errors, double smooth, bool constCoeffs);

  Double_t evaluate() const;

  ClassDef(RooCubicSplineFun,1) // CubicSpline polynomial PDF
};

#endif

/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooCubicSplineGaussModel.h,v 1.21 2007/05/11 09:13:07 verkerke Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_CS_GAUSS_MODEL
#define ROO_CS_GAUSS_MODEL

#include <complex>
#include <RooResolutionModel.h>
#include <RooRealProxy.h>

#include <P2VV/RooAbsEffResModel.h>

class RooCubicSplineFun;
class RooAbsReal;
class RooRealVar;

class RooCubicSplineGaussModel : public RooResolutionModel, public RooAbsEffResModel {
public:

  // Constructors, assignment etc
  inline RooCubicSplineGaussModel() : _flatSFInt(kFALSE)  { }
  RooCubicSplineGaussModel(const char *name, const char *title, 
            RooRealVar& x, RooCubicSplineFun& spline,
            RooAbsReal& mean,   RooAbsReal& sigma );
  RooCubicSplineGaussModel(const char *name, const char *title, 
            RooRealVar& x, RooCubicSplineFun& spline,
            RooAbsReal& mean,   RooAbsReal& sigma,
            RooAbsReal& meanSF, RooAbsReal& sigmaSF) ; 
  RooCubicSplineGaussModel(const RooCubicSplineGaussModel& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooCubicSplineGaussModel(*this,newname) ; }
  virtual ~RooCubicSplineGaussModel();
  
  virtual Int_t basisCode(const char* name) const ;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName) const ;

  void advertiseFlatScaleFactorIntegral(Bool_t flag) { _flatSFInt = flag ; }

  virtual const RooAbsReal* efficiency() const;
  virtual std::vector<const RooAbsReal*> efficiencies() const;
  virtual const RooArgSet* observables() const;

  // TODO: push these into RooCubicSplineFUn...
  class M_n {
  public:
        M_n(double x, const std::complex<double>& z) ;
        const std::complex<double>& operator()(int i) const { assert(0<=i&&i<4); return _m[i]; }
        M_n& operator-=(const M_n& other) { for(int i=0;i<4;++i) _m[i]= _m[i]-other._m[i]; return *this; }
        M_n  operator- (const M_n& other) const { return M_n(*this)-=other; }
  private:
        std::complex<double> _m[4];
  };
  class K_n {
  public:
      K_n(const std::complex<double>& z) : _zi( std::complex<double>(1,0)/z) {}
      std::complex<double> operator()(unsigned i) const {
          assert(0<=i&&i<=3);
          switch(i) {
              case 0 : return 0.5*_zi;
              case 1 : return 0.5*_zi*_zi;
              case 2 : return _zi*(_zi*_zi+1.0);
              case 3 : std::complex<double> _zi2 = _zi*_zi; return _zi2*(3.*_zi2+3.);
          }
          assert(1==0);
          return 0;
      }
  private :
      std::complex<double> _zi;        
  };

private:

  virtual Double_t evaluate() const ;

  std::complex<double> evalInt(Double_t xmin, Double_t xmax, const std::complex<double>& z) const;

  Bool_t _flatSFInt ;
  
  RooRealProxy eff ;
  RooRealProxy mean ;
  RooRealProxy sigma ;
  RooRealProxy msf ;
  RooRealProxy ssf ;

  ClassDef(RooCubicSplineGaussModel,1)
};

#endif

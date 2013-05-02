/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooGaussModel.cxx 44982 2012-07-10 08:36:13Z moneta $
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

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Class RooGaussModel implements a RooResolutionModel that models a Gaussian
// distribution. Object of class RooGaussModel can be used
// for analytical convolutions with classes inheriting from RooAbsAnaConvPdf
// After  the convolution  is applied, one can then multiply the result
// with a cubic b-spline...
// END_HTML
//

#include "RooFit.h"

#include "TMath.h"
#include "Riostream.h"
#include "P2VV/RooCubicSplineGaussModel.h"
#include "RooMath.h"
#include "RooComplex.h"
#include "RooRealConstant.h"
#include "RooRandom.h"

using namespace std;

ClassImp(RooCubicSplineGaussModel) 
;


namespace {

    enum basisType { noBasis=0  ,      expBasis= 3
                   , sinBasis=13,  cosBasis=23
                   , sinhBasis=63, coshBasis=53 };
    static Double_t root2(sqrt(2.)) ;
    static Double_t root2pi(sqrt(2*atan2(0.,-1.))) ;
    static Double_t rootpi(sqrt(atan2(0.,-1.))) ;

    RooComplex evalCerfApprox(Double_t x, const RooComplex& z) {
      // compute exp(-x*x)erfc(-i(z-x)), 
      // use the approximation: erfc(z) = exp(-z*z)/(sqrt(pi)*z)
      // to explicitly cancel the divergent exp(y*y) behaviour of
      // CWERF for z = x + i y with large negative y
      // erfc(-i(z-x) = exp((z-x)*(z-x))/(sqrt(pi)*(-i)*(z-x))
      // 
      RooComplex mi(0,-1);
      RooComplex zp  = mi*(z-x);
      RooComplex zsq = zp*zp;
      RooComplex v = -zsq - x*x;
      RooComplex iz(z.im()+x,z.re()-x);
      return v.exp()*(zsq.exp()/(iz*rootpi) + 1)*2 ;
    }

    // Calculate Re[exp(-x^2) cwerf(i (z-x) )], taking care of numerical instabilities
    Double_t evalCerfRe(Double_t x, const RooComplex& z) {
      Double_t re =  z.re()-x;
      return (re>-4.0) ? RooMath::FastComplexErrFuncRe(RooComplex(-z.im(),re))*exp(-x*x) 
                       : evalCerfApprox(x,z).re() ;
    }

    // Calculate Im(exp(-u^2) cwerf(cwt + i(u+c))), taking care of numerical instabilities
    Double_t evalCerfIm(Double_t x, const RooComplex& z) {
      Double_t re = z.re()-x;
      return (re>-4.0) ? RooMath::FastComplexErrFuncIm(RooComplex(-z.im(),re))*exp(-x*x) 
                       : evalCerfApprox(x,z).im() ;
    }

    // Calculate exp(-x^2) cwerf(cwt + i(u+c)), taking care of numerical instabilities
    RooComplex evalCerf(Double_t x, const RooComplex& z) {
      Double_t re = z.re()-x;
      return (re>-4.0) ? RooMath::FastComplexErrFunc(RooComplex(-z.im(),re))*exp(-x*x) 
                       : evalCerfApprox(x,z) ;
    }

    RooComplex evalCerfInt(Double_t xmin, Double_t xmax, const RooComplex& z) {
      RooComplex d = ( evalCerf(xmin,z) - evalCerf(xmax,z) ) + ( RooMath::erf(xmax) - RooMath::erf(xmin) );
      return d/(z*2);
    }

}

//_____________________________________________________________________________
RooCubicSplineGaussModel::RooCubicSplineGaussModel(const char *name, const char *title, RooRealVar& xIn, 
			     RooAbsReal& _mean, RooAbsReal& _sigma) :
  RooResolutionModel(name,title,xIn), 
  _flatSFInt(kFALSE),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma),
  msf("msf","Mean Scale Factor",this,RooRealConstant::value(1)),
  ssf("ssf","Sigma Scale Factor",this,RooRealConstant::value(1))
{  
}



//_____________________________________________________________________________
RooCubicSplineGaussModel::RooCubicSplineGaussModel(const char *name, const char *title, RooRealVar& xIn, 
			     RooAbsReal& _mean, RooAbsReal& _sigma, 
			     RooAbsReal& _msSF) : 
  RooResolutionModel(name,title,xIn), 
  _flatSFInt(kFALSE),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma),
  msf("msf","Mean Scale Factor",this,_msSF),
  ssf("ssf","Sigma Scale Factor",this,_msSF)
{  
}



//_____________________________________________________________________________
RooCubicSplineGaussModel::RooCubicSplineGaussModel(const char *name, const char *title, RooRealVar& xIn, 
			     RooAbsReal& _mean, RooAbsReal& _sigma, 
			     RooAbsReal& _meanSF, RooAbsReal& _sigmaSF) : 
  RooResolutionModel(name,title,xIn), 
  _flatSFInt(kFALSE),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma),
  msf("msf","Mean Scale Factor",this,_meanSF),
  ssf("ssf","Sigma Scale Factor",this,_sigmaSF)
{  
}



//_____________________________________________________________________________
RooCubicSplineGaussModel::RooCubicSplineGaussModel(const RooCubicSplineGaussModel& other, const char* name) : 
  RooResolutionModel(other,name),
  _flatSFInt(other._flatSFInt),
  mean("mean",this,other.mean),
  sigma("sigma",this,other.sigma),
  msf("msf",this,other.msf),
  ssf("ssf",this,other.ssf)
{
}

//_____________________________________________________________________________
RooCubicSplineGaussModel::~RooCubicSplineGaussModel()
{
  // Destructor
}

//_____________________________________________________________________________
Int_t RooCubicSplineGaussModel::basisCode(const char* name) const 
{
  if (!TString("exp(-@0/@1)"              ).CompareTo(name)) return expBasis;
  if (!TString("exp(-@0/@1)*sin(@0*@2)"   ).CompareTo(name)) return sinBasis;
  if (!TString("exp(-@0/@1)*cos(@0*@2)"   ).CompareTo(name)) return cosBasis;
  if (!TString("exp(-@0/@1)*cosh(@0*@2/2)").CompareTo(name)) return coshBasis;
  if (!TString("exp(-@0/@1)*sinh(@0*@2/2)").CompareTo(name)) return sinhBasis;
  return 0 ;
} 

//_____________________________________________________________________________
Double_t RooCubicSplineGaussModel::evaluate() const 
{  
  basisType basisCode = (basisType) _basisCode ;

  Double_t tau    = (basisCode!=noBasis)                       ? ((RooAbsReal*)basis().getParameter(1))->getVal() : 0 ;
  Double_t omega  = (basisCode==sinBasis  || basisCode==cosBasis)  ? ((RooAbsReal*)basis().getParameter(2))->getVal() : 0 ;
  Double_t dGamma = (basisCode==sinhBasis || basisCode==coshBasis) ? ((RooAbsReal*)basis().getParameter(2))->getVal() : 0 ;

  if (basisCode  == coshBasis && basisCode!=noBasis && dGamma==0 ) basisCode = expBasis;

  Double_t scale = sigma*ssf*root2;
  Double_t u = (x-(mean*msf))/scale;
  // *** 1st form: Straight Gaussian, used for unconvoluted PDF or expBasis with 0 lifetime ***
  if (basisCode==noBasis || ((basisCode==expBasis || basisCode==cosBasis) && tau==0)) {
    if (verboseEval()>2) cout << "RooCubicSplineGaussModel::evaluate(" << GetName() << ") 1st form" << endl ;
    return exp(-u*u)/(scale*rootpi) ;
  }

  // *** 2nd form: 0, used for sinBasis, linBasis, and quadBasis with tau=0 ***
  if (tau==0) {
    if (verboseEval()>2) cout << "RooCubicSplineGaussModel::evaluate(" << GetName() << ") 2nd form" << endl ;
    return 0. ;
  }

  RooComplex z( double(1)/tau, -omega ); z=z*scale/2; // frac{ ( Gamma - I omega ) sigma }{ sqrt{2} }
 
  // *** 5th form: Convolution with exp(- gamma t *cos(omega*t), used for cosBasis, and expBasis (for which omega=0) ***
  if (basisCode==cosBasis || basisCode == expBasis) {
    if (verboseEval()>2) cout << "RooCubicSplineGaussModel::evaluate(" << GetName() << ") 5th form z = " << z << ", u = " << u << endl ;
    Double_t result = evalCerfRe(u,z);
    if (TMath::IsNaN(result)) cxcoutE(Tracing) << "RooCubicSplineGaussModel::evaluate(" << GetName() << ") got nan during case 4; z="<<z<<" u=" << u << " res=" << result << endl; 
    return result ;  
  }

  // *** 4th form: Convolution with exp(- gamma t )*sin(omega*t), used for sinBasis(omega<>0) ***
  if (basisCode==sinBasis) {
    if (verboseEval()>2) cout << "RooCubicSplineGaussModel::evaluate(" << GetName() << ") 4th form z = " << z << ", u = " << u << endl ;
    Double_t result = z.im()!=0 ? evalCerfIm(u,z) : 0;
    if (TMath::IsNaN(result)) cxcoutE(Tracing) << "RooCubicSplineGaussModel::evaluate(" << GetName() << ") got nan during case 3 " << endl; 
    return result ;
  }

  if (basisCode==coshBasis || basisCode ==sinhBasis) {
    // ***8th form: Convolution with exp(-|t|/tau)*cosh(dgamma*t/2), used for         coshBasisSum ***
    if (verboseEval()>2) cout << "RooCubicSplineGaussModel::evaluate(" << GetName() << ") 8th form tau = " << tau << endl ;
    int sgn = ( basisCode == coshBasis ? +1 : -1 );
    Double_t dgammah = ((RooAbsReal*)basis().getParameter(2))->getVal() / 2;
    dgammah *= scale/2;
    Double_t  result = (evalCerfRe(u,z-RooComplex(dgammah,0))+sgn*evalCerfRe(u,z+RooComplex(dgammah,0)))/2 ; 
    if (TMath::IsNaN(result)) cxcoutE(Tracing) << "RooCubicSplineGaussModel::evaluate(" << GetName() << ") got nan during case 8 " << endl; 
    return result;
  }
  assert(0) ;
  return 0 ;
}

//_____________________________________________________________________________
Int_t RooCubicSplineGaussModel::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const 
{
  switch(_basisCode) {
  // Analytical integration capability of raw PDF
  case noBasis:
    // Optionally advertise flat integral over sigma scale factor
    if (_flatSFInt && matchArgs(allVars,analVars,RooArgSet(convVar(),ssf.arg()))) return 2 ;
    if (matchArgs(allVars,analVars,convVar())) return 1 ;
    break ;

  // Analytical integration capability of convoluted PDF
  case expBasis:
  case sinBasis:
  case cosBasis:
  case coshBasis:
  case sinhBasis:
    // Optionally advertise flat integral over sigma scale factor
    if (_flatSFInt && matchArgs(allVars,analVars,RooArgSet(convVar(),ssf.arg()))) return 2 ;
    if (matchArgs(allVars,analVars,convVar())) return 1 ;
    break ;
  }
  return 0 ;
}

//_____________________________________________________________________________
Double_t RooCubicSplineGaussModel::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  // Code must be 1 or 2
  assert(code==1||code==2) ;
  Double_t ssfInt = code==2 ? (ssf.max(rangeName)-ssf.min(rangeName)) : 1.0 ;

  basisType basisCode = (basisType) _basisCode ;

  // *** 1st form: Straight Gaussian, used for unconvoluted PDF or expBasis with 0 lifetime ***
  Double_t tau = (basisCode!=noBasis)?((RooAbsReal*)basis().getParameter(1))->getVal():0 ;
  if (basisCode == coshBasis && basisCode!=noBasis ) {
     Double_t dGamma = ((RooAbsReal*)basis().getParameter(2))->getVal();
     if (dGamma==0) basisCode = expBasis;
  }

  Double_t scale = sigma*ssf*root2;
  Double_t umin = (x.min(rangeName)-(mean*msf))/scale;
  Double_t umax = (x.max(rangeName)-(mean*msf))/scale;

  if (basisCode==noBasis || ((basisCode==expBasis || basisCode==cosBasis) && tau==0)) {
    if (verboseEval()>0) cout << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") 1st form" << endl ;
    Double_t result =  0.5*(RooMath::erf( umax )-RooMath::erf( umin )) ;
    if (TMath::IsNaN(result)) { cxcoutE(Tracing) << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") got nan during case 1 " << endl; }
    return result*ssfInt ;
  }
  // *** 2nd form: unity, used for sinBasis and linBasis with tau=0 (PDF is zero) ***
  if (tau==0) {
    if (verboseEval()>0) cout << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") 2nd form" << endl ;
    return 0. ;
  }

  Double_t omega  = ((basisCode==sinBasis )||(basisCode==cosBasis )) ? ((RooAbsReal*)basis().getParameter(2))->getVal() : 0 ;
  RooComplex z( double(1)/tau, -omega ); z=z*scale/2; // frac{ ( Gamma - I omega ) sigma }{ sqrt{2} }

  // *** 3rd form: Convolution with exp(-t/tau)*cos(omega*t), used for cosBasis(omega<>0) and expBasis(omega=0) ***
  if (basisCode==cosBasis || basisCode==expBasis) {
    if (verboseEval()>0) cout << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") 3rd form z = " << z << endl ;
    Double_t result = evalCerfInt(umin,umax,z).re();
    if (TMath::IsNaN(result)) { cxcoutE(Tracing) << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") got nan during case 3 " << endl; }
    result *= scale;
    cout << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") 5th form umin=" << umin << ", umax=" << umax << " -> " << result << endl;
    return result*ssfInt;
  }
    
  // *** 4th form: Convolution with exp(-t/tau)*sin(omega*t), used for sinBasis(omega<>0,tau<>0) ***
  if (basisCode==sinBasis) {    
    if (verboseEval()>0) cout << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") 4th form omega = " << omega << ", tau = " << tau << endl ;
    Double_t result = z.im()!=0 ? evalCerfInt(umin,umax,z).im() : 0 ;
    if (TMath::IsNaN(result)) { cxcoutE(Tracing) << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") got nan during case 4 " << endl; }
    result *= scale;
    cout << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") 4th form umin=" << umin << ", umax=" << umax << " -> " << result << endl;
    return result*ssfInt;
  }

  // *** 5th form: Convolution with exp(-|t|/tau)*{cosh,sinh}(dgamma*t/2), used for {coshBasis,sinhBasis} ***
  if (basisCode==coshBasis || basisCode == sinhBasis) {
    if (verboseEval()>0) {cout << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName()                             << ") 5th form tau=" << tau << endl ; }
    Double_t dgammah =  ((RooAbsReal*)basis().getParameter(2))->getVal() / 2 ;
    dgammah *= scale;
    int sgn = ( basisCode == coshBasis ? +1 : -1 );
    Double_t result = ( evalCerfInt(umin,umax,z-RooComplex(dgammah,0)).re()+ sgn*evalCerfInt(umin,umax,z+RooComplex(dgammah,0)).re())/2;
    if (TMath::IsNaN(result)) { cxcoutE(Tracing) << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") got nan during case 5 " << endl; }
    result *= scale;
    cout << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") 8th form umin=" << umin << ", umax=" << umax << " -> " << result << endl;
    return result*ssfInt;
  }

  assert(0) ;
  return 0 ;
}

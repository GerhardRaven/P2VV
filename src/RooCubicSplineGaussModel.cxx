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

#include <memory>

#include "RooFit.h"

#include "TMath.h"
#include "Riostream.h"
#include "P2VV/RooCubicSplineGaussModel.h"
#include "P2VV/RooCubicSplineFun.h"
#include "RooMath.h"
#include "RooComplex.h"
#include "RooRealConstant.h"
#include "RooRandom.h"

using namespace std;

ClassImp(RooCubicSplineGaussModel) 
;

namespace {
    enum basisType { noBasis=0  ,  expBasis= 3
                   , sinBasis=13,  cosBasis=23
                   , sinhBasis=63, coshBasis=53 };
    static const Double_t rootpi(sqrt(TMath::Pi())) ;

    RooComplex evalApprox(Double_t x, const RooComplex& z) {
      // compute exp(-x^2)cwerf(-i(z-x)), cwerf(z) = exp(-z^2)erfc(-iz)
      // use the approximation: erfc(z) = exp(-z*z)/(sqrt(pi)*z)
      // to explicitly cancel the divergent exp(y*y) behaviour of
      // CWERF for z = x + i y with large negative y
      static const RooComplex mi(0,-1);
      RooComplex zp  = mi*(z-x);
      RooComplex zsq = zp*zp;
      RooComplex v = -zsq - x*x;
      RooComplex iz(z.im()+x,z.re()-x); // ???
      return v.exp()*(zsq.exp()/(iz*rootpi) + 1)*2 ;
    }

    // Calculate Re[exp(-x^2) cwerf(i (z-x) )], taking care of numerical instabilities
    Double_t evalRe(Double_t x, const RooComplex& z) {
      Double_t re =  z.re()-x;
      return (re>-5.0) ? RooMath::FastComplexErrFuncRe(RooComplex(-z.im(),re))*exp(-x*x) 
                       : evalApprox(x,z).re() ;
    }

    // Calculate Im[exp(-x^2) cwerf(i(z-x))], taking care of numerical instabilities
    Double_t evalIm(Double_t x, const RooComplex& z) {
      Double_t re = z.re()-x;
      return (re>-5.0) ? RooMath::FastComplexErrFuncIm(RooComplex(-z.im(),re))*exp(-x*x) 
                       : evalApprox(x,z).im() ;
    }

    // Calculate exp(-x^2) cwerf(i(z-x)), taking care of numerical instabilities
    RooComplex eval(Double_t x, const RooComplex& z) {
      Double_t re = z.re()-x;
      return (re>-5.0) ? RooMath::FastComplexErrFunc(RooComplex(-z.im(),re))*exp(-x*x) 
                       : evalApprox(x,z) ;
    }

    class L_jk {
    public:
        L_jk(double x) : _x(x) { }
        double operator()(int j, int k) const { 
            assert(0<=j&&j<4);
            assert(0<=k&&k<3);
            switch(k) {
                case 0: return j==0 ? 1 : 0 ;
                case 1: switch(j) {
                        case 0 : return  0;
                        case 1 : return -2/rootpi;
                        case 2 : return -4*_x/rootpi;
                        case 3 : return -4*(_x*_x-1)/rootpi;
                        default : assert(1==0); return 0;
                }
                case 2: switch(j) {
                        case 0 : return -1;
                        case 1 : return -2*_x;
                        case 2 : return -2*(2*_x*_x-1);
                        case 3 : return -4*_x*(2*_x*_x-3);
                        default : assert(1==0); return 0;
            }   }
            assert(1==0);
            return 0;
        }  
    private : 
        double _x;
    };

}

RooCubicSplineGaussModel::M_n::M_n(double x, const RooComplex& z) {
          RooComplex N0( RooMath::erf(x) )
                   , N1( exp(-x*x)       )
                   , N2( eval(x,z)       );
          L_jk L(x); 
          for (int i=0;i<4;++i) _m[i] = N0*L(i,0) + N1*L(i,1) + N2*L(i,2);
}


//_____________________________________________________________________________
RooCubicSplineGaussModel::RooCubicSplineGaussModel(const char *name, const char *title
                , RooRealVar& x, RooCubicSplineFun& _spline
			    , RooAbsReal& _mean, RooAbsReal& _sigma ) :
  RooResolutionModel(name,title,x), 
  _flatSFInt(kFALSE),
  spline("spline","Spline describing efficiency",this,_spline),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma),
  msf("msf","Mean Scale Factor",this,RooRealConstant::value(1.0)),
  ssf("ssf","Sigma Scale Factor",this,RooRealConstant::value(1.0))
{  
    // make sure 'x' matches the spline argument!
    std::auto_ptr<RooArgSet> svar( spline.arg().getVariables() );
    assert( svar->contains( convVar() ) );
}

//_____________________________________________________________________________
RooCubicSplineGaussModel::RooCubicSplineGaussModel(const char *name, const char *title
                , RooRealVar& x, RooCubicSplineFun& _spline
			    , RooAbsReal& _mean, RooAbsReal& _sigma
			    , RooAbsReal& _meanSF, RooAbsReal& _sigmaSF) : 
  RooResolutionModel(name,title,x), 
  _flatSFInt(kFALSE),
  spline("spline","Spline describing efficiency",this,_spline),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma),
  msf("msf","Mean Scale Factor",this,_meanSF),
  ssf("ssf","Sigma Scale Factor",this,_sigmaSF)
{  
    // make sure 'x' matches the spline argument!
    std::auto_ptr<RooArgSet> svar( spline.arg().getVariables() );
    assert( svar->contains( convVar() ) );
}

//_____________________________________________________________________________
RooCubicSplineGaussModel::RooCubicSplineGaussModel(const RooCubicSplineGaussModel& other, const char* name) : 
  RooResolutionModel(other,name),
  _flatSFInt(other._flatSFInt),
  spline("spline",this,other.spline),
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
Double_t RooCubicSplineGaussModel::efficiency() const 
{
    // TODO: provide a RooAbsReal so we can plot the efficiency -- just like createIntegral...
    return spline;
}

RooComplex RooCubicSplineGaussModel::evalInt(Double_t umin, Double_t umax, const RooComplex& z) const
{
    const RooCubicSplineFun &sp = dynamic_cast<const RooCubicSplineFun&>( spline.arg() );
    // TODO: verify we remain within [umin,umax]
    RooCubicSplineGaussModel::K_n K(z);
    Double_t scale = sigma*ssf*TMath::Sqrt2(); 
    Double_t offset = mean*msf;
    std::vector<M_n> M; M.reserve( sp.knotSize() );
    for (int i=0;i<sp.knotSize();++i) M.push_back( M_n( (sp.u(i)-offset)/scale, z ) );
    RooComplex sum(0,0);
    double sc[4]; for (int i=0;i<4;++i) sc[i] = pow(scale,i);
    for (int i=0;i<sp.knotSize()-1;++i) { //TODO: push this loop into RooCubicSplineFun... pass z,scale,coefficients
        sum = sum + sp.gaussIntegral( i,  M[i+1] - M[i], K, offset, sc );
    }
    return sum;
}

//_____________________________________________________________________________
Double_t RooCubicSplineGaussModel::evaluate() const 
{  
  basisType basisCode = (basisType) _basisCode ;
  Double_t tau    = (basisCode!=noBasis)                           ? ((RooAbsReal*)basis().getParameter(1))->getVal() : 0 ;
  Double_t omega  = (basisCode==sinBasis  || basisCode==cosBasis)  ? ((RooAbsReal*)basis().getParameter(2))->getVal() : 0 ;
  Double_t dGamma = (basisCode==sinhBasis || basisCode==coshBasis) ? ((RooAbsReal*)basis().getParameter(2))->getVal() : 0 ;
  if (basisCode  == coshBasis && basisCode!=noBasis && dGamma==0 ) basisCode = expBasis;

  Double_t scale = sigma*ssf*TMath::Sqrt2();
  Double_t u = (x-mean*msf)/scale;
  // *** 1st form: Straight Gaussian, used for unconvoluted PDF or expBasis with 0 lifetime ***
  if (basisCode==noBasis || ((basisCode==expBasis || basisCode==cosBasis) && tau==0)) {
    if (verboseEval()>2) cout << "RooCubicSplineGaussModel::evaluate(" << GetName() << ") 1st form" << endl ;
    return exp(-u*u)/(scale*rootpi) ; // ???
  }

  // *** 2nd form: 0, used for sinBasis, linBasis, and quadBasis with tau=0 ***
  if (tau==0) {
    if (verboseEval()>2) cout << "RooCubicSplineGaussModel::evaluate(" << GetName() << ") 2nd form" << endl ;
    return 0. ;
  }
  RooComplex z( double(1)/tau, -omega ); z=z*scale/2;
 
  Double_t val(0);
  if (verboseEval()>2) cout << "RooCubicSplineGaussModel::evaluate(" << GetName() << ") basisCode = " <<  basisCode << " z = " << z << ", u = " << u << endl ;

  switch (basisCode) {
    case expBasis:
    case cosBasis:
        val +=             evalRe(u,z);
        break;
    case sinBasis:
        val += z.im()!=0 ? evalIm(u,z) : 0;
        break;
    case coshBasis:
    case sinhBasis: {
        RooComplex y( scale * dGamma / 4 , 0 );
        val += (                                      evalRe(u,z-y)
               + ( basisCode == coshBasis ? +1 : -1 )*evalRe(u,z+y) )/2; 
        break;
    }
    default:
        assert(0);
  }
  if (TMath::IsNaN(val)) cxcoutE(Tracing) << "RooCubicSplineGaussModel::evaluate(" << GetName() << ") got nan during basisCode = " << basisCode << endl; 
  Double_t eff=efficiency();
  if (TMath::IsNaN(eff)) cxcoutE(Tracing) << "RooCubicSplineGaussModel::evaluate(" << GetName() << ") got nan during efficiency " << endl;
  return eff*val;
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
  Double_t ssfInt( code==2 ? (ssf.max(rangeName)-ssf.min(rangeName)) : 1.0 );

  basisType basisCode = (basisType) _basisCode ;
  Double_t tau    = (basisCode!=noBasis)                           ? ((RooAbsReal*)basis().getParameter(1))->getVal() : 0 ;
  Double_t omega  = (basisCode==sinBasis  || basisCode==cosBasis)  ? ((RooAbsReal*)basis().getParameter(2))->getVal() : 0 ;
  Double_t dGamma = (basisCode==sinhBasis || basisCode==coshBasis) ? ((RooAbsReal*)basis().getParameter(2))->getVal() : 0 ;
  if (basisCode == coshBasis && basisCode!=noBasis && dGamma==0 ) basisCode = expBasis;

  Double_t scale = sigma*ssf*TMath::Sqrt2();
  Double_t umin = (x.min(rangeName)-mean*msf)/scale;
  Double_t umax = (x.max(rangeName)-mean*msf)/scale;

  if (basisCode==noBasis || ((basisCode==expBasis || basisCode==cosBasis) && tau==0)) {
    if (verboseEval()>0) cout << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") 1st form" << endl ;
    Double_t result =  0.5*(RooMath::erf( umax )-RooMath::erf( umin )) ;
    if (TMath::IsNaN(result)) { cxcoutE(Tracing) << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") got nan during case 1 " << endl; }
    return result*ssfInt ; // ???
  }
  if (tau==0) {
    if (verboseEval()>0) cout << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") 2nd form" << endl ;
    return 0. ;
  }

  RooComplex z( double(1)/tau, -omega ); z=z*scale/2;

  if (verboseEval()>0) cout << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") basisCode = " << basisCode << " z = " << z << endl ;

  Double_t result(0);
  switch (basisCode) {
    case expBasis:
    case cosBasis:
        result +=             evalInt(umin,umax,z).re();
        break;
    case sinBasis:
        result += z.im()!=0 ? evalInt(umin,umax,z).im() : 0 ;
        break;
    case coshBasis:
    case sinhBasis: {
        RooComplex y( scale * dGamma / 4 , 0 );
        result += (                                      evalInt(umin,umax,z-y).re()
                  + ( basisCode == coshBasis ? +1 : -1 )*evalInt(umin,umax,z+y).re() )/2;
        break;
    }
    default: 
        assert(0) ;
  }
  if (TMath::IsNaN(result)) { cxcoutE(Tracing) << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") got nan for basisCode = " << basisCode << endl; }
  return scale*result*ssfInt;
}

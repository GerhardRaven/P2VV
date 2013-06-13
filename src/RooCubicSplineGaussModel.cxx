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

    std::complex<double> evalApprox(Double_t x, const std::complex<double>& z) {
      // compute exp(-x^2)cwerf(-i(z-x)), cwerf(z) = exp(-z^2)erfc(-iz)
      // use the approximation: erfc(z) = exp(-z*z)/(sqrt(pi)*z)
      // to explicitly cancel the divergent exp(y*y) behaviour of
      // CWERF for z = x + i y with large negative y
      static const std::complex<double> mi(0,-1);
      std::complex<double> zp  = mi*(z-x);
      std::complex<double> zsq = zp*zp;
      std::complex<double> v = -zsq - x*x;
      std::complex<double> iz(z.imag()+x,z.real()-x); // ???
      return exp(v)*(exp(zsq)/(iz*rootpi) + 1.)*2. ;
    }

    // Calculate Re[exp(-x^2) cwerf(i (z-x) )], taking care of numerical instabilities
    Double_t evalRe(Double_t x, const std::complex<double>& z) {
      Double_t re =  z.real()-x;
      return (re>-5.0) ? RooMath::faddeeva_fast(std::complex<double>(-z.imag(),re)).real()*exp(-x*x) 
                       : evalApprox(x,z).real() ;
    }

    // Calculate Im[exp(-x^2) cwerf(i(z-x))], taking care of numerical instabilities
    Double_t evalIm(Double_t x, const std::complex<double>& z) {
      Double_t re = z.real()-x;
      return (re>-5.0) ? RooMath::faddeeva_fast(std::complex<double>(-z.imag(),re)).imag()*exp(-x*x) 
                       : evalApprox(x,z).imag() ;
    }

    // Calculate exp(-x^2) cwerf(i(z-x)), taking care of numerical instabilities
    std::complex<double> eval(Double_t x, const std::complex<double>& z) {
      Double_t re = z.real()-x;
      return (re>-5.0) ? RooMath::faddeeva_fast(std::complex<double>(-z.imag(),re))*exp(-x*x) 
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
                        default : return 2*(*this)(j-1,2)/rootpi;
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

RooCubicSplineGaussModel::M_n::M_n(double x, const std::complex<double>& z) {
          std::complex<double> N0( RooMath::erf(x) )
                   , N1( exp(-x*x)       )
                   , N2( eval(x,z)       );
          // TODO: eliminate L_jk all together...
          L_jk L(x); 
          for (int i=0;i<4;++i) _m[i] = N0*L(i,0) + N1*L(i,1) + N2*L(i,2);
}


//_____________________________________________________________________________
RooCubicSplineGaussModel::RooCubicSplineGaussModel(const char *name, const char *title
                , RooRealVar& x, RooCubicSplineFun& _eff
			    , RooAbsReal& _mean, RooAbsReal& _sigma ) 
   : RooResolutionModel(name, title, x),
     RooAbsEffResModel(), 
     _flatSFInt(kFALSE),
     eff("eff","Spline describing efficiency",this,_eff),
     mean("mean","Mean",this,_mean),
     sigma("sigma","Width",this,_sigma),
     msf("msf","Mean Scale Factor",this,RooRealConstant::value(1.0)),
     ssf("ssf","Sigma Scale Factor",this,RooRealConstant::value(1.0))
{  
   // make sure 'x' matches the eff argument!
   std::auto_ptr<RooArgSet> svar( eff.arg().getVariables() );
   assert( svar->contains( convVar() ) );
}

//_____________________________________________________________________________
RooCubicSplineGaussModel::RooCubicSplineGaussModel(const char *name, const char *title
                                                   , RooRealVar& x, RooCubicSplineFun& _eff
                                                   , RooAbsReal& _mean, RooAbsReal& _sigma
                                                   , RooAbsReal& _meanSF, RooAbsReal& _sigmaSF)
   : RooResolutionModel(name,title,x), 
     RooAbsEffResModel(),
     _flatSFInt(kFALSE),
     eff("eff","Spline describing efficiency",this,_eff),
     mean("mean","Mean",this,_mean),
     sigma("sigma","Width",this,_sigma),
     msf("msf","Mean Scale Factor",this,_meanSF),
     ssf("ssf","Sigma Scale Factor",this,_sigmaSF)
{  
   // make sure 'x' matches the spline argument!
   std::auto_ptr<RooArgSet> svar( eff.arg().getVariables() );
   assert( svar->contains( convVar() ) );
}

//_____________________________________________________________________________
RooCubicSplineGaussModel::RooCubicSplineGaussModel(const RooCubicSplineGaussModel& other,
                                                   const char* name)
   : RooResolutionModel(other,name),
     RooAbsEffResModel(),
     _flatSFInt(other._flatSFInt),
     eff("eff",this,other.eff),
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
const RooAbsReal* RooCubicSplineGaussModel::efficiency() const 
{
    return const_cast<RooAbsReal*>(&eff.arg());
}

//_____________________________________________________________________________
std::vector<const RooAbsReal*> RooCubicSplineGaussModel::efficiencies() const { 
   return std::vector<const RooAbsReal*>(1, efficiency());
}

//_____________________________________________________________________________
const RooArgSet* RooCubicSplineGaussModel::observables() const { 
   // Return pointer to pdf in product
   return new RooArgSet(convVar());
}

//_____________________________________________________________________________
std::complex<double> RooCubicSplineGaussModel::evalInt(Double_t umin, Double_t umax, const std::complex<double>& z) const
{
    const RooCubicSplineFun &sp = dynamic_cast<const RooCubicSplineFun&>( eff.arg() );
    //TODO: verify we remain within [umin,umax]
    //TODO: push this loop into RooCubicSplineFun... pass z,scale,offset and umin,umax
    RooCubicSplineGaussModel::K_n K(z);
    Double_t scale = sigma*ssf*TMath::Sqrt2(); 
    Double_t offset = mean*msf;
    assert(sp.knotSize()>1);
    std::vector<M_n> M; M.reserve( sp.knotSize() );
    for (unsigned int i=0;i<sp.knotSize();++i) {
        double u = (sp.u(i)-offset)/scale ;
        assert( u>=umin );
        assert( u<=umax );
        M.push_back( M_n( (sp.u(i)-offset)/scale, z ) );
    }
    //cout << "knots run from " << sp.u(0) << " ( " << (sp.u(0)-offset)/scale  << " )" 
    //             <<    " to " << sp.u(sp.knotSize()-1) << " ( " << (sp.u(sp.knotSize()-1)-offset)/scale << ")" <<endl;
    //cout << "requested integral goes from " << umin << " to " << umax << endl;
    double sc[4]; for (int i=0;i<4;++i) sc[i] = pow(scale,i);
    std::complex<double> sum(0,0);
    for (unsigned i=0;i<sp.knotSize()-1;++i) {
        complex<double> I = sp.gaussIntegral( i,  M[i+1] - M[i], K, offset, sc );
        //cout << " analytical integral from " << sp.u(i) << " to " << sp.u(i+1) << " = " <<  I << endl;
        sum += I ;
    }

    double lo = scale*umin+offset;
    if (lo<sp.u(0)) {
        complex<double> I = sp.gaussIntegralE(true,  M.front() - M_n( umin,z) , K, offset, sc);
        // cout << " analytical integral from " << lo << " to " << sp.u(0) << " = " <<  I << endl;
        sum += I;
    }
    double hi = scale*umax+offset;
    if (hi>sp.u(sp.knotSize()-1)) {
        complex<double> I = sp.gaussIntegralE(false,  M_n(umax,z)-M.back() , K, offset, sc);
        // cout << " analytical integral from " << sp.u(sp.knotSize()-1) << " to " << hi << " = " <<  I << endl;
        sum += I;
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
    Double_t eff=efficiency()->getVal();
    if (TMath::IsNaN(eff)) 
       cxcoutE(Tracing) << "RooCubicSplineGaussModel::evaluate(" << GetName() 
                        << ") got nan during efficiency " << endl;
    return eff * exp(-u*u)/(scale*rootpi) ; // ???
  }

  // *** 2nd form: 0, used for sinBasis, linBasis, and quadBasis with tau=0 ***
  if (tau==0) {
    if (verboseEval()>2) cout << "RooCubicSplineGaussModel::evaluate(" << GetName() << ") 2nd form" << endl ;
    return 0. ;
  }
  std::complex<double> z( double(1)/tau, -omega ); z*=0.5*scale;
 
  Double_t val(0);
  if (verboseEval()>2) cout << "RooCubicSplineGaussModel::evaluate(" << GetName() << ") basisCode = " <<  basisCode << " z = " << z << ", u = " << u << endl ;

  switch (basisCode) {
    case expBasis:
    case cosBasis:
        val +=             evalRe(u,z);
        break;
    case sinBasis:
        val += z.imag()!=0 ? evalIm(u,z) : 0;
        break;
    case coshBasis:
    case sinhBasis: {
        std::complex<double> y( scale * dGamma / 4 , 0 );
        val += (                                      evalRe(u,z-y)
               + ( basisCode == coshBasis ? +1 : -1 )*evalRe(u,z+y) )/2; 
        break;
    }
    default:
        assert(0);
  }
  if (TMath::IsNaN(val)) 
     cxcoutE(Tracing) << "RooCubicSplineGaussModel::evaluate(" << GetName() 
                      << ") got nan during basisCode = " << basisCode << endl; 
  Double_t _eff=eff;
  if (TMath::IsNaN(_eff)) 
     cxcoutE(Tracing) << "RooCubicSplineGaussModel::evaluate(" << GetName() 
                      << ") got nan during efficiency " << endl;
  // static int i=10; if (i>0) { --i; cout << "RooCubicSplineGaussModel("<<GetName()<<")::evaluate : " << _eff << " * " << val << " = " << _eff*val << endl; }
  return _eff*val;
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

  std::complex<double> z( double(1)/tau, -omega ); z=0.5*z*scale;

  if (verboseEval()>0) cout << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") basisCode = " << basisCode << " z = " << z << endl ;

  Double_t result(0);
  switch (basisCode) {
    case expBasis:
    case cosBasis:
        result +=             evalInt(umin,umax,z).real();
        break;
    case sinBasis:
        result += z.imag()!=0 ? evalInt(umin,umax,z).imag() : 0 ;
        break;
    case coshBasis:
    case sinhBasis: {
        std::complex<double> y( scale * dGamma / 4 , 0 );
        result += 0.5 * (                                      evalInt(umin,umax,z-y).real()
                        + ( basisCode == coshBasis ? +1 : -1 )*evalInt(umin,umax,z+y).real() );
        break;
    }
    default: 
        assert(0) ;
  }
  if (TMath::IsNaN(result)) { cxcoutE(Tracing) << "RooCubicSplineGaussModel::analyticalIntegral(" << GetName() << ") got nan for basisCode = " << basisCode << endl; }
  return scale*result*ssfInt;
}

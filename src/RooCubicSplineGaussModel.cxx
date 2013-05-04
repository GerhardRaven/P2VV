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
#include "P2VV/RooCubicSplineKnot.h"
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
    static const Double_t root2(sqrt(2.)) ;
    static const Double_t rootpi(sqrt(atan2(0.,-1.))) ;

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
      return (re>-4.0) ? RooMath::FastComplexErrFuncRe(RooComplex(-z.im(),re))*exp(-x*x) 
                       : evalApprox(x,z).re() ;
    }

    // Calculate Im[exp(-x^2) cwerf(i(z-x))], taking care of numerical instabilities
    Double_t evalIm(Double_t x, const RooComplex& z) {
      Double_t re = z.re()-x;
      return (re>-4.0) ? RooMath::FastComplexErrFuncIm(RooComplex(-z.im(),re))*exp(-x*x) 
                       : evalApprox(x,z).im() ;
    }

    // Calculate exp(-x^2) cwerf(i(z-x)), taking care of numerical instabilities
    RooComplex eval(Double_t x, const RooComplex& z) {
      Double_t re = z.re()-x;
      return (re>-4.0) ? RooMath::FastComplexErrFunc(RooComplex(-z.im(),re))*exp(-x*x) 
                       : evalApprox(x,z) ;
    }


    class K_n {
    public:
        K_n(const RooComplex& z) : _zi( RooComplex(1,0)/z) {}
        RooComplex operator()(int i) const {
            assert(0<=i&&i<=3);
            switch(i) {
                case 0 : return _zi/2;
                case 1 : return _zi*_zi/2;
                case 2 : return _zi*(_zi*_zi+1);
                case 3 : RooComplex _zi2 = _zi*_zi; 
                         return _zi2*(RooComplex(3,0)+RooComplex(3,0)*_zi2);
            }
            assert(1==0);
            return 0;
        }
    private :
        RooComplex _zi;        
    };

    class N_n { 
    public:
        N_n(double x, RooComplex z) {
            _n[0] = RooMath::erf(x);
            _n[1] = exp(-x*x);
            _n[2] = eval(x,z);
        }
        const RooComplex& operator()(int i) const {
            assert(0<=i&&i<3);
            return _n[i];
        }
    private:
        RooComplex _n[3];
    };

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

    class M_n {
    public:
       M_n(double x, const RooComplex& z) {
          L_jk L(x) ;
          N_n  N(x,z) ;
          for (int i=0;i<4;++i) _m[i] = N(0)*L(i,0) 
                                      + N(1)*L(i,1) 
                                      + N(2)*L(i,2);
       }
       const RooComplex& operator()(int i) {
           assert(0<=i&&i<4); 
           return _m[i];
       }
       M_n& operator-=(const M_n& other) { for(int i=0;i<4;++i) _m[i]= _m[i]-other._m[i]; return *this; }
       M_n  operator- (const M_n& other) const { return M_n(*this)-=other; }
    private:
       RooComplex _m[4];
    };


}

//_____________________________________________________________________________
RooCubicSplineGaussModel::RooCubicSplineGaussModel(const char *name, const char *title
                , RooRealVar& xIn, const char *knotBinningName, const RooArgList& coefList
                , RooAbsReal& _mean, RooAbsReal& _sigma) :
  RooResolutionModel(name,title,xIn), 
  _flatSFInt(kFALSE),
  splineCoefficients("splineCoeffiients","List of splineCoeffiients",this),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma),
  msf("msf","Mean Scale Factor",this,RooRealConstant::value(1)),
  ssf("ssf","Sigma Scale Factor",this,RooRealConstant::value(1)),
  knots(0)
{  
  // TODO: verify coefList is consistent with knots as specified by the knotBinningName binning
  //    should be N+2 coefficients for N bins...
  const RooAbsBinning* binning = xIn.getBinningPtr(knotBinningName);
  assert( binning!=0);
  assert( coefList.getSize()==3+binning->numBins());
 
  // Constructor
  TIterator* coefIter = coefList.createIterator() ;
  RooAbsArg* coef ;
  while((coef = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      cout << "RooCubicBSpline::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName() 
	   << " is not of type RooRealVar" << endl ;
      assert(0) ;
    }
    splineCoefficients.add(*coef) ;
  }
  delete coefIter ;
  Double_t* boundaries = binning->array();
  knots = new RooCubicSplineKnot( boundaries, boundaries + binning->numBoundaries() );
}

//_____________________________________________________________________________
RooCubicSplineGaussModel::RooCubicSplineGaussModel(const char *name, const char *title
                , RooRealVar& xIn
			    , RooAbsReal& _mean, RooAbsReal& _sigma
			    , RooAbsReal& _meanSF, RooAbsReal& _sigmaSF) : 
  RooResolutionModel(name,title,xIn), 
  _flatSFInt(kFALSE),
  splineCoefficients("splineCoeffiients","List of splineCoeffiients",this),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma),
  msf("msf","Mean Scale Factor",this,_meanSF),
  ssf("ssf","Sigma Scale Factor",this,_sigmaSF),
  knots(0)
{  
}

//_____________________________________________________________________________
RooCubicSplineGaussModel::RooCubicSplineGaussModel(const RooCubicSplineGaussModel& other, const char* name) : 
  RooResolutionModel(other,name),
  _flatSFInt(other._flatSFInt),
  splineCoefficients("splineCoeffiients",this,other.splineCoefficients),
  mean("mean",this,other.mean),
  sigma("sigma",this,other.sigma),
  msf("msf",this,other.msf),
  ssf("ssf",this,other.ssf),
  knots(new RooCubicSplineKnot(*other.knots) )
{
}

//_____________________________________________________________________________
RooCubicSplineGaussModel::~RooCubicSplineGaussModel()
{
  // Destructor
  delete knots;
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
    return knots->evaluate(x,splineCoefficients);
}

RooComplex RooCubicSplineGaussModel::evalInt(Double_t umin, Double_t umax, const RooComplex& z) const
{
    K_n K(z);
    Double_t scale = sigma*ssf*root2; 
    std::vector<M_n> M; M.reserve( knots->size() );
    for (int i=0;i<knots->size();++i) M.push_back( M_n( knots->u(i)/scale, z ) );
    RooComplex sum(0,0);
    for (int i=0;i<knots->size()-1;++i) {
        RooCubicSplineKnot::S_jk S( knots->S_jk_sum( i, splineCoefficients ) );
        M_n dM( M[i+1] - M[i] );
        for (int j=0;j<4;++j) for (int k=0;k<4-j;++k) {
            sum = sum + dM(j)*S(j,k)*K(k)*pow(scale,j+k);
        }
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

  Double_t scale = sigma*ssf*root2;
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

  Double_t scale = sigma*ssf*root2;
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

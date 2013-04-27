/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooBSpline.cxx 45780 2012-08-31 15:45:27Z moneta $
 * Authors:                                                                  *
 *   Kyle Cranmer
 *                                                                           *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// BSpline basis polynomials are positive-definite in the range [0,1].
// In this implementation, we extend [0,1] to be the range of the parameter.
// There are n+1 BSpline basis polynomials of degree n.
// Thus, by providing N coefficients that are positive-definite, there 
// is a natural way to have well bahaved polynomail PDFs.
// For any n, the n+1 basis polynomials 'form a partition of unity', eg.
//  they sum to one for all values of x. See
// END_HTML
//

#include "RooFit.h"

#include "Riostream.h"
#include <math.h>
#include "TMath.h"
#include "P2VV/RooCubicBSpline.h"
#include "RooMath.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooComplex.h"
#include "RooArgList.h"

#include <iterator>
#include <algorithm>

using namespace std;

namespace _aux {
  // compute  e^{-x^2} w(I(z-x))
  //                   w(?) = FastComplexErrFunc(?) = exp(-?^2)erfc(-i?) aka. Faddeeva function
  RooComplex evalCerf(Double_t x, const RooComplex& z) {
    RooComplex _z( -z.im(), z.re()-x);
    if (_z.im()>-4.0) return RooMath::FastComplexErrFunc(_z)*exp(-x*x) ;

    // use the approximation: erf(z) = exp(-z*z)/(sqrt(pi)*z)
    // to explicitly cancel the divergent exp(y*y) behaviour of
    // CWERF for z = x + i y with large negative y
    static Double_t rootpi= sqrt(atan2(0.,-1.));
    RooComplex  iz(-z.re()-x,z.im()); // \frac{i}{\sqrt{2}}(z-x)
    RooComplex zsq=z*z;
    RooComplex v= zsq - x*x;

    return v.exp()*(-zsq.exp()/(iz*rootpi)+1)*2 ;
  }

    class L_jk {
    public:
        L_jk(double x) : _x(x) { }
        double operator()(int j, int k) const { 
            assert(0<=j&&j<4);
            assert(0<=k&&k<3);
            static const double N( double(1)/sqrt(atan2(1.0,0.0)) ); // sqrt(2/pi)
            switch(k) {
                case 0: return j==0 ? 1 : 0 ;
                case 1: switch(j) {
                        case 0 : return  0;
                        case 1 : return -N;
                        case 2 : return -N*_x;
                        case 3 : return -N*(_x*_x-1);
                        default : assert(1==0); return 0;
                }
                case 2: switch(j) {
                        case 0 : return -1;
                        case 1 : return -_x;
                        case 2 : return -(_x*_x-1);
                        case 3 : return -_x*(_x*_x-3);
                        default : assert(1==0); return 0;
            }   }
            assert(1==0);
            return 0;
        }  
    private : 
        double _x;
    };

    class K_n {
    public:
        K_n(const RooComplex& z) : _zi( RooComplex(1,0)/z) {}
        K_n(double re, double im) : _zi( RooComplex(1,0)/RooComplex(re,im) ) {}
        RooComplex operator()(int i) const {
            assert(0<=i&&i<=3);
            switch(i) {
                case 0 : return _zi;
                case 1 : return _zi*_zi;
                case 2 : return _zi*(RooComplex(1,0)+RooComplex(2,0)*_zi*_zi);
                case 3 : RooComplex _zi2 = _zi*_zi; 
                         return _zi2*(RooComplex(3,0)+RooComplex(6,0)*_zi2);
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
            static const double N( double(1)/sqrt(2.0) );
            x *= N; z = z*N;
            _n[0] = RooMath::erf(x);
            _n[1] = exp(-x*x);
            _n[2] = evalCerf(x, RooComplex(z.im(),-z.re())); // exp(-x^2) w(I(z-x) )
        }
        const RooComplex& operator()(int i) const {
            assert(0<=i&&i<3);
            return _n[i];
        }
    private:
        RooComplex _n[3];
    };

    class M_n {
    public:
       M_n(double x, const RooComplex& z) {
          _aux::L_jk L(x) ;
          _aux::N_n  N(x,z) ;
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

    class S_jk { 
    public:
        S_jk(double a, double b, double c) : t(a*b*c), d(a*b+a*c+b*c), s(a+b+c), o(1.) {}
        S_jk& operator*=(double z) { t*=z; d*=z; s*=z; o*=z; return *this; } 
        S_jk& operator/=(double z) { t/=z; d/=z; s/=z; o/=z; return *this; } 
        S_jk& operator-()          { t=-t; d=-d; s=-s; o=-o; return *this; }
        S_jk& operator+=(const S_jk& other) { t+=other.t; d+=other.d; s+=other.s; o+=other.o; return *this; } 
        S_jk& operator-=(const S_jk& other) { t-=other.t; d-=other.d; s-=other.s; o-=other.o; return *this; } 

        S_jk operator*(double z)          const { return S_jk(*this)*=z; }
        S_jk operator/(double z)          const { return S_jk(*this)/=z; }
        S_jk operator+(const S_jk& other) const { return S_jk(*this)+=other; }
        S_jk operator-(const S_jk& other) const { return S_jk(*this)-=other; }

        double operator()(int j, int k) const { 
            assert(0<=j&&j<4);
            assert(0<=k&&k<4-j); // note: for 4-j<=k<4 could return 0... but better not to invoke with those..
            if (j>k) std::swap(j,k);
            switch(3*j+k) {
                case 0: return   -t;   // (0,0) 
                case 1: return    d;   // (0,1),(1,0)
                case 2: return   -s;   // (0,2),(2,0)
                case 3: return    o;   // (0,3),(3,0)
                case 4: return -2*s;   // (1,1)
                case 5: return  3*o;   // (1,2),(2,1)
                default : assert(1==0);
            }
        }
    private:
        double t,d,s,o;
    };

}

class RooCubicBSpline::_auxKnot {
public:
    template <typename Iter> _auxKnot(Iter begin, Iter end) : _u(begin,end) {
           // P,Q,R,S only depend on the knot vector, so build at construction, and cache them...
           _PQRS.reserve(4*_u.size());
           for (int i=0;i<_u.size();++i) { 
                _PQRS.push_back( h(i+1,i-2)*h(i+1,i-1)*h(i+1,i) );
                _PQRS.push_back( h(i+1,i-1)*h(i+2,i-1)*h(i+1,i) );
                _PQRS.push_back( h(i+2,i  )*h(i+2,i-1)*h(i+1,i) );
                _PQRS.push_back( h(i+2,i  )*h(i+3,i  )*h(i+1,i) );
           }
    }

    int index(double u) const { 
        assert(u>=_u.front() && u<=_u.back());   
        std::vector<double>::const_iterator i = --std::upper_bound(_u.begin(),_u.end()-1,u);
        assert( _u.begin()<=i );
        assert( *i <= u && u<=*(i+1) );
        return std::distance(_u.begin(),i);
    };

    double evaluate(double _u, const RooArgList& b) const {
        int i = index(_u); // location in knot vector
        assert(0<=i && i+3<b.getSize());
        assert( u(i) <= _u && _u<= u(i+1) );
        return  ((RooAbsReal&)b[i  ]).getVal()*A(_u,i)
             +  ((RooAbsReal&)b[i+1]).getVal()*B(_u,i)
             +  ((RooAbsReal&)b[i+2]).getVal()*C(_u,i)
             +  ((RooAbsReal&)b[i+3]).getVal()*D(_u,i);
    }

    double analyticalIntegral(const RooArgList& b) const {
        if (_IABCD.empty()) {
           // the integrals of A,B,C,D from u(i) to u(i+1) only depend on the knot vector...
           // so we create them 'on demand' and cache the result
           _IABCD.reserve(4*_u.size());
           for (int j=0;j<_u.size();++j) { 
                _IABCD.push_back(   qua(h(j,j+1))/(4*P(j)) ) ;
                _IABCD.push_back( - cub(h(j,j+1))*(3*u(j)-4*u(j-2)+u(j+1))/(12*P(j))
                                  - sqr(h(j,j+1))*(3*sqr(u(j))-2*u(j-1)*u(j+1)+sqr(u(j+1))+u(j)*(-4*u(j-1)+2*u(j+1)-4*u(j+2)) +6*u(j-1)*u(j+2)-2*u(j+1)*u(j+2) )/(12*Q(j))
                                  + sqr(h(j,j+1))*(3*sqr(u(j+1))+sqr(u(j  ))+2*u(j)*u(j+1)-8*u(j+1)*u(j+2)-4*u(j  )*u(j+2)+6*sqr(u(j+2)))/(12*R(j)) );
                _IABCD.push_back(   sqr(h(j,j+1))*(3*sqr(u(j  ))+sqr(u(j+1))+2*u(j+1)*u(j)-8*u(j  )*u(j-1)-4*u(j-1)*u(j+1)+6*sqr(u(j-1)))/(12*Q(j))
                                  - sqr(h(j,j+1))*(3*sqr(u(j+1))+sqr(u(j))-4*u(j-1)*u(j+1)+6*u(j-1)*u(j+2)-4*u(j+1)*u(j+2)-2*u(j)*(u(j-1)-u(j+1)+u(j+2)))/(12*R(j))
                                  + cub(h(j,j+1))*(3*u(j+1)-4*u(j+3)+u(j))/(12*S(j)) );
                _IABCD.push_back(   qua(h(j,j+1))/(4*S(j)) );
            }
        }
        double norm(0);
        assert(b.getSize()-2==_u.size());
        for (int i=0; i < _u.size()-1; ++i) {
             norm += ((RooAbsReal&)b[i  ]).getVal()*_IABCD[4*i  ]
                   + ((RooAbsReal&)b[i+1]).getVal()*_IABCD[4*i+1]
                   + ((RooAbsReal&)b[i+2]).getVal()*_IABCD[4*i+2]
                   + ((RooAbsReal&)b[i+3]).getVal()*_IABCD[4*i+3];
        }
        return norm;
    }

    RooComplex analyticalIntegral(const RooComplex& z, const RooArgList& coef) const {
        std::vector<_aux::M_n> M; M.reserve(_u.size());
        for (int i=0;i<_u.size();++i) M.push_back( _aux::M_n( u(i), z ) );
        _aux::K_n K(z);
        RooComplex sum(0,0);
        for (int i=0;i<_u.size()-1;++i) {
            _aux::S_jk S( S_jk_sum( i, coef ) );
            _aux::M_n dM( M[i+1] - M[i] );
            for (int j=0;j<4;++j) for (int k=0;k<4-j;++k) sum = sum + dM(j)*S(j,k)*K(k);
        }
        return sum;
    }
private:
    double A(double _u,int i) const{ return -cub(d(_u,i+1))/P(i); }
    double B(double _u,int i) const{ return  sqr(d(_u,i+1))*d(_u,i-2)/P(i) + d(_u,i-1)*d(_u,i+2)*d(_u,i+1)/Q(i) + d(_u,i  )*sqr(d(_u,i+2))/R(i); }
    double C(double _u,int i) const{ return -sqr(d(_u,i-1))*d(_u,i+1)/Q(i) - d(_u,i  )*d(_u,i+2)*d(_u,i-1)/R(i) - d(_u,i+3)*sqr(d(_u,i  ))/S(i); }
    double D(double _u,int i) const{ return  cub(d(_u,i  ))/S(i); }

    double P(int i) const { assert(4*i  <_PQRS.size()); return  _PQRS[4*i  ]; }
    double Q(int i) const { assert(4*i+1<_PQRS.size()); return  _PQRS[4*i+1]; }
    double R(int i) const { assert(4*i+2<_PQRS.size()); return  _PQRS[4*i+2]; }
    double S(int i) const { assert(4*i+3<_PQRS.size()); return  _PQRS[4*i+3]; }



    // S matrix for i-th interval
    _aux::S_jk S_jk_sum(int i, const RooArgList& b) const { 
            if (_S_jk.empty()) {
                _S_jk.reserve(_u.size()*4);
                for(int i=0;i<_u.size();++i) {
                    _S_jk.push_back( -_aux::S_jk(u(i+1),u(i+1),u(i+1))/P(i) );
                    _S_jk.push_back(  _aux::S_jk(u(i-2),u(i+1),u(i+1))/P(i)
                                     +_aux::S_jk(u(i-1),u(i+1),u(i+2))/Q(i) 
                                     +_aux::S_jk(u(i  ),u(i+2),u(i+2))/R(i) );
                    _S_jk.push_back( -_aux::S_jk(u(i-1),u(i-1),u(i+1))/Q(i)
                                     -_aux::S_jk(u(i-1),u(i  ),u(i+2))/R(i) 
                                     -_aux::S_jk(u(i  ),u(i  ),u(i+3))/S(i) );
                    _S_jk.push_back(  _aux::S_jk(u(i  ),u(i  ),u(i  ))/S(i) );
                }
            }
            return _S_jk[4*i  ]*((RooAbsReal&)b[i  ]).getVal()
                 + _S_jk[4*i+1]*((RooAbsReal&)b[i+1]).getVal()
                 + _S_jk[4*i+2]*((RooAbsReal&)b[i+2]).getVal()
                 + _S_jk[4*i+3]*((RooAbsReal&)b[i+3]).getVal(); 
    }
    double sqr(double x) const { return x*x; }
    double cub(double x) const { return x*sqr(x); }
    double qua(double x) const { return sqr(sqr(x)); }
    double d(double _u, int j) const { return _u-u(j); }
    double h(int i, int j) const { return u(i)-u(j); }
    double u(int i) const { assert(i>-3&&i<int(_u.size()+3)); return _u[std::min(std::max(0,i),int(_u.size()-1))]; }

    const   std::vector<double> _u;
    mutable std::vector<double> _PQRS;
    mutable std::vector<double> _IABCD;
    mutable std::vector<_aux::S_jk> _S_jk;
};

ClassImp(RooCubicBSpline)
;

//_____________________________________________________________________________
RooCubicBSpline::RooCubicBSpline()
    : _aux(0)
{
}

//_____________________________________________________________________________
RooCubicBSpline::RooCubicBSpline(const char* name, const char* title, 
                           RooRealVar& x, const char* knotBinningName, const RooArgList& coefList): 
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _coefList("coefficients","List of coefficients",this),
  _aux(0)
{
  // TODO: verify coefList is consistent with knots as specified by the knotBinningName binning
  //    should be N+2 coefficients for N bins...
  const RooAbsBinning* binning = x.getBinningPtr(knotBinningName);
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
    _coefList.add(*coef) ;
  }
  delete coefIter ;
  Double_t* boundaries = binning->array();
  _aux = new _auxKnot( boundaries, boundaries + binning->numBoundaries() );
}

//_____________________________________________________________________________
RooCubicBSpline::RooCubicBSpline(const RooCubicBSpline& other, const char* name) :
  RooAbsPdf(other, name), 
  _x("x", this, other._x), 
  _coefList("coefList",this,other._coefList),
  _aux(new _auxKnot(*other._aux))
{
}

//_____________________________________________________________________________
RooCubicBSpline::~RooCubicBSpline()
{
    delete _aux;
}

//_____________________________________________________________________________
Double_t RooCubicBSpline::evaluate() const 
{
  return _aux->evaluate(_x,_coefList);
}

//_____________________________________________________________________________
Int_t RooCubicBSpline::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
  // No analytical calculation available (yet) of integrals over subranges
  if (rangeName && strlen(rangeName)) return 0 ;
  if (matchArgs(allVars, analVars, _x)) return 1;
  return 0;
}


//_____________________________________________________________________________
Double_t RooCubicBSpline::analyticalIntegral(Int_t code, const char* rangeName) const
{
  assert(code==1) ;
  return _aux->analyticalIntegral(_coefList);
}


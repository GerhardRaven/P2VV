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

  Double_t get(const RooArgList& b,int i) { return ((RooAbsReal&)b[i]).getVal() ; }
  Double_t get(const RooArgList& b,int i,int k) { return _aux::get(b,i+k); }

  template <typename T> typename T::const_reference get(const T& t, int i, int j) { return t[4*i+j]; }
  template <typename T> void push_back(T& t, const typename T::value_type& a,
                                             const typename T::value_type& b, 
                                             const typename T::value_type& c,
                                             const typename T::value_type& d) { t.push_back(a); t.push_back(b); t.push_back(c); t.push_back(d) ; }
 
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

    double evaluate(double _u, const RooArgList& b) const {
        int i = index(_u); // location in knot vector
        assert(0<=i && i+3<b.getSize());
        assert( u(i) <= _u && _u<= u(i+1) );
        return  _aux::get(b,i,0)*A(_u,i) // TODO: substitute A,B,C,D 'in situ'
             +  _aux::get(b,i,1)*B(_u,i)
             +  _aux::get(b,i,2)*C(_u,i)
             +  _aux::get(b,i,3)*D(_u,i);
    }

    double analyticalIntegral(const RooArgList& b) const {
        if (_IABCD.empty()) {
           // the integrals of A,B,C,D from u(i) to u(i+1) only depend on the knot vector...
           // so we create them 'on demand' and cache the result
           _IABCD.reserve(4*_u.size());
           for (int j=0;j<_u.size();++j) { 
               _aux::push_back(_IABCD,   qua(h(j,j+1))/(4*P(j)) 
                                     , - cub(h(j,j+1))*(3*u(j)-4*u(j-2)+u(j+1))/(12*P(j))
                                       - sqr(h(j,j+1))*(3*sqr(u(j))-2*u(j-1)*u(j+1)+sqr(u(j+1))+u(j)*(-4*u(j-1)+2*u(j+1)-4*u(j+2)) +6*u(j-1)*u(j+2)-2*u(j+1)*u(j+2) )/(12*Q(j))
                                       + sqr(h(j,j+1))*(3*sqr(u(j+1))+sqr(u(j  ))+2*u(j)*u(j+1)-8*u(j+1)*u(j+2)-4*u(j  )*u(j+2)+6*sqr(u(j+2)))/(12*R(j))
                                     ,   sqr(h(j,j+1))*(3*sqr(u(j  ))+sqr(u(j+1))+2*u(j+1)*u(j)-8*u(j  )*u(j-1)-4*u(j-1)*u(j+1)+6*sqr(u(j-1)))/(12*Q(j))
                                       - sqr(h(j,j+1))*(3*sqr(u(j+1))+sqr(u(j))-4*u(j-1)*u(j+1)+6*u(j-1)*u(j+2)-4*u(j+1)*u(j+2)-2*u(j)*(u(j-1)-u(j+1)+u(j+2)))/(12*R(j))
                                       + cub(h(j,j+1))*(3*u(j+1)-4*u(j+3)+u(j))/(12*S(j)) 
                                     ,   qua(h(j,j+1))/(4*S(j)) );
            }
        }
        assert(b.getSize()-2==_u.size());
        double norm(0);
        for (int i=0; i < _u.size()-1; ++i) for (int k=0;k<4;++k) {
            norm += _aux::get(b,i,k)*_aux::get(_IABCD,i,k) ;
        }
        return norm;
    }

private:
    int index(double u) const { 
        assert(u>=_u.front() && u<=_u.back());   
        std::vector<double>::const_iterator i = --std::upper_bound(_u.begin(),_u.end()-1,u);
        assert( _u.begin()<=i );
        assert( *i <= u && u<=*(i+1) );
        return std::distance(_u.begin(),i);
    };
    double A(double _u,int i) const{ return -cub(d(_u,i+1))/P(i); }
    double B(double _u,int i) const{ return  sqr(d(_u,i+1))*d(_u,i-2)/P(i) + d(_u,i-1)*d(_u,i+2)*d(_u,i+1)/Q(i) + d(_u,i  )*sqr(d(_u,i+2))/R(i); }
    double C(double _u,int i) const{ return -sqr(d(_u,i-1))*d(_u,i+1)/Q(i) - d(_u,i  )*d(_u,i+2)*d(_u,i-1)/R(i) - d(_u,i+3)*sqr(d(_u,i  ))/S(i); }
    double D(double _u,int i) const{ return  cub(d(_u,i  ))/S(i); }

    double P(int i) const { assert(4*i  <_PQRS.size()); return  _PQRS[4*i  ]; }
    double Q(int i) const { assert(4*i+1<_PQRS.size()); return  _PQRS[4*i+1]; }
    double R(int i) const { assert(4*i+2<_PQRS.size()); return  _PQRS[4*i+2]; }
    double S(int i) const { assert(4*i+3<_PQRS.size()); return  _PQRS[4*i+3]; }



    double sqr(double x) const { return x*x; }
    double cub(double x) const { return x*sqr(x); }
    double qua(double x) const { return sqr(sqr(x)); }
    double d(double _u, int j) const { return _u-u(j); }
    double h(int i, int j) const { return u(i)-u(j); }
    double u(int i) const { assert(i>-3&&i<int(_u.size()+3)); return _u[std::min(std::max(0,i),int(_u.size()-1))]; }

    const   std::vector<double> _u;
    mutable std::vector<double> _PQRS;
    mutable std::vector<double> _IABCD;
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


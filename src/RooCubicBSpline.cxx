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
// http://www.idav.ucdavis.edu/education/CAGDNotes/BSpline-Polynomials.pdf
// END_HTML
//

#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>
#include "TMath.h"
#include "P2VV/RooCubicBSpline.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"

#include <iterator>
#include <algorithm>

using namespace std;

class RooCubicBSpline::_auxKnot {
public:
    _auxKnot(const _auxKnot& other) : _u(other._u) {}
    _auxKnot(const RooRealVar& x, const char *binningName) {
        // create knotvector from specified binning
        const RooAbsBinning* binning = x.getBinningPtr(binningName);
        assert(binning!=0);
        int n = binning->numBoundaries();
        Double_t* bin = binning->array();
        _u.insert( _u.begin(), bin, bin+n ) ;
    }

    int index(double u) const { assert(u>=_u.front() && u<=_u.back());   
        std::vector<double>::const_iterator i = --std::upper_bound(_u.begin(),_u.end()-1,u);
        assert( _u.begin()<=i );
        assert( *i <= u && u<=*(i+1) );
        return std::distance(_u.begin(),i);
    };

    double evaluate(double _u, int i, double b0, double b1, double b2, double b3) const {
        assert( u(i) <= _u && _u<= u(i+1) );
        return b0*A(_u,i) + b1*B(_u,i) + b2*C(_u,i) + b3*D(_u,i);
    }

    double analyticalIntegral(int i, double b0, double b1, double b2, double b3) const {
        return b0*IA(i)+b1*IB(i)+b2*IC(i)+b3*ID(i); 
    }

private:

    double A(double _u,int i) const{ return -cub(d(_u,i+1))/P(i); }
    double B(double _u,int i) const{ return  sqr(d(_u,i+1))*d(_u,i-2)/P(i) + d(_u,i-1)*d(_u,i+2)*d(_u,i+1)/Q(i) + d(_u,i  )*sqr(d(_u,i+2))/R(i); }
    double C(double _u,int i) const{ return -sqr(d(_u,i-1))*d(_u,i+1)/Q(i) - d(_u,i  )*d(_u,i+2)*d(_u,i-1)/R(i) - d(_u,i+3)*sqr(d(_u,i  ))/S(i); }
    double D(double _u,int i) const{ return  cub(d(_u,i  ))/S(i); }





    // integrals from u(i) to u(i+1) -- note: these are constant one we know the knotvector, so can (should) be cached...
    double IA(int i) const{  return   qua(h(i,i+1))/(4*P(i)); }
    double IB(int i) const{  return - cub(h(i,i+1))*(3*u(i)-4*u(i-2)+u(i+1))/(12*P(i))
                                    - sqr(h(i,i+1))*(3*sqr(u(i))-2*u(i-1)*u(i+1)+sqr(u(i+1))+u(i)*(-4*u(i-1)+2*u(i+1)-4*u(i+2)) +6*u(i-1)*u(i+2)-2*u(i+1)*u(i+2) )/(12*Q(i))
                                    + sqr(h(i,i+1))*(3*sqr(u(i+1))+sqr(u(i  ))+2*u(i)*u(i+1)-8*u(i+1)*u(i+2)-4*u(i  )*u(i+2)+6*sqr(u(i+2)))/(12*R(i)); }
    double IC(int i) const{  return   sqr(h(i,i+1))*(3*sqr(u(i  ))+sqr(u(i+1))+2*u(i+1)*u(i)-8*u(i  )*u(i-1)-4*u(i-1)*u(i+1)+6*sqr(u(i-1)))/(12*Q(i))
                                    - sqr(h(i,i+1))*(3*sqr(u(i+1))+sqr(u(i))-4*u(i-1)*u(i+1)+6*u(i-1)*u(i+2)-4*u(i+1)*u(i+2)-2*u(i)*(u(i-1)-u(i+1)+u(i+2)))/(12*R(i))
                                    + cub(h(i,i+1))*(3*u(i+1)-4*u(i+3)+u(i))/(12*S(i)); }
    double ID(int i) const{  return   qua(h(i,i+1))/(4*S(i)); }

    // P,Q,R,S only depend on the knot vector, so build at construction, and cache them...
    double P(int i) const { return  h(i+1,i-2)*h(i+1,i-1)*h(i+1,i); }
    double Q(int i) const { return  h(i+1,i-1)*h(i+2,i-1)*h(i+1,i); }
    double R(int i) const { return  h(i+2,i  )*h(i+2,i-1)*h(i+1,i); }
    double S(int i) const { return  h(i+2,i  )*h(i+3,i  )*h(i+1,i); }


    class S_jk { 
    public:
        S_jk() : t(0),d(0),s(0) {}
        S_jk(double a, double b, double c) : t(a*b*c), d(a*b+a*c+b*c), s(a+b+c) {}
        S_jk& operator*=(double z) {  t*=z; d*=z; s*=z; return *this; } 
        S_jk& operator/=(double z) {  t/=z; d/=z; s/=z; return *this; } 
        S_jk& operator-()          {  t=-t; d=-d; s=-s; return *this; }
        S_jk& operator+=(const S_jk& other) {  t+=other.t; d+=other.d; s+=other.s; return *this; } 
        S_jk& operator-=(const S_jk& other) {  t-=other.t; d-=other.d; s-=other.s; return *this; } 

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
                case 3: return    1;   // (0,3),(3,0)
                case 4: return -2*s;   // (1,1)
                case 5: return    3;   // (1,2),(2,1)
                default : assert(1==0);
            }
        }
    private:
        double t,d,s ;
    };

    // S matrix for i-th interval
    S_jk S_jk_sum(int i, double b0, double b1, double b2, double b3) const {
            // TODO: make sure that the 4 S_jk are cached for each i...
            return (-S_jk(u(i+1),u(i+1),u(i+1))/P(i) )*b0  // A_i triplet * b0
                  +( S_jk(u(i-2),u(i+1),u(i+1))/P(i)+S_jk(u(i-1),u(i+1),u(i+2))/Q(i) +S_jk(u(i),u(i+2),u(i+2))/R(i))*b1   // B_i * b1
                  +(-S_jk(u(i-1),u(i-1),u(i+1))/Q(i)-S_jk(u(i-1),u(i  ),u(i+2))/R(i) -S_jk(u(i),u(i  ),u(i+3))/S(i))*b2  // C_i * b2
                  +( S_jk(u(i),u(i),u(i))/S(i) )*b3;       // D_i * b3
    }

    double sqr(double x) const { return x*x; }
    double cub(double x) const { return x*sqr(x); }
    double qua(double x) const { return sqr(sqr(x)); }
    double d(double _u, int j) const { return _u-u(j); }
    double h(int i, int j) const { return u(i)-u(j); }
    double u(int i) const { assert(i>=-3&&i<int(_u.size()+3)); i=std::min(std::max(0,i),int(_u.size()-1));  return _u[i]; }

    std::vector<double> _u;
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
  _aux = new _auxKnot( (RooRealVar&)_x.arg(), knotBinningName);
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
  double u = _x;
  int i = _aux->index(u); // location in knot vector
  assert(0<=i && i+3<_coefList.getSize());
  return _aux->evaluate(u,i, ((RooAbsReal*)_coefList.at(i  ))->getVal()
                           , ((RooAbsReal*)_coefList.at(i+1))->getVal()
                           , ((RooAbsReal*)_coefList.at(i+2))->getVal()
                           , ((RooAbsReal*)_coefList.at(i+3))->getVal() );
}


//_____________________________________________________________________________
Int_t RooCubicBSpline::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
  // No analytical calculation available (yet) of integrals over subranges
  if (rangeName && strlen(rangeName)) {
    return 0 ;
  }

  if (matchArgs(allVars, analVars, _x)) return 1;
  return 0;
}


//_____________________________________________________________________________
Double_t RooCubicBSpline::analyticalIntegral(Int_t code, const char* rangeName) const
{
  assert(code==1) ;
  Double_t  norm(0);
  for (int i=0; i < _coefList.getSize()-3; ++i) {
    norm += _aux->analyticalIntegral(i, ((RooAbsReal*)_coefList.at(i  ))->getVal()
                                      , ((RooAbsReal*)_coefList.at(i+1))->getVal()
                                      , ((RooAbsReal*)_coefList.at(i+2))->getVal()
                                      , ((RooAbsReal*)_coefList.at(i+3))->getVal() ) ;
  }
  return norm;
}


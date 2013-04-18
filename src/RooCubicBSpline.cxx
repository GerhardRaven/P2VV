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
#include "RooCubicBSpline.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"

#include <iterator>

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
        _u.resize( n+6 );
        for (int i=0;i<3;++i) { _u[i] = bin[0] ; _u[3+n+i] = bin[n-1]; }
        for (int i=0;i<n;++i) { _u[3+i] = bin[i] ; }

        //for (std::vector<double>::const_iterator i = _u.begin();i!=_u.end();++i) cout << *i << endl;
    }
    double P(int i) const{ assert(i>1&&i<_u.size()-1);  return _P(i); }
    double Q(int i) const{ assert(i>0&&i<_u.size()-2);  return _Q(i); } 
    double R(int i) const{ assert(i>0&&i<_u.size()-2 ); return _R(i); }
    double S(int i) const{ assert(i>=0&&i<_u.size()-3); return _S(i); }
    double u(int i) const{ return _u[i]; }

    int i(double u) const { assert(u>=_u.front() && u<=_u.back());   
        std::vector<double>::const_iterator _i = std::min(std::upper_bound(_u.begin(),_u.end(),u),_u.end()-4);
        assert(std::distance(_i,_u.end())>2);
        assert(std::distance(_u.begin(),_i)>2);
        assert( u<=*_i && u>=*(_i-1));
        return std::distance(_u.begin(),_i)-1;
    };

    double A(double _u,int i) const{ assert( _u>=u(i)&&_u<=u(i+1)) ; return cub(u(i+1)-_u)/P(i); }
    double B(double _u,int i) const{ assert( _u>=u(i)&&_u<=u(i+1)) ; return sqr(u(i+1)-_u)*(_u-u(i-2))/P(i) + (_u-u(i-1))*(u(i+2)-_u)*(u(i+1)-_u)/Q(i) + (_u-u(i))*pow(u(i+2)-_u,2)/R(i); }
    double C(double _u,int i) const{ assert( _u>=u(i)&&_u<=u(i+1)) ; return (u(i+1)-_u)*sqr(_u-u(i-1))/Q(i) + (_u-u(i  ))*(u(i+2)-_u)*(_u-u(i-1))/R(i) + pow(_u-u(i),2)*(u(i+3)-_u)/S(i); }
    double D(double _u,int i) const{ assert( _u>=u(i)&&_u<=u(i+1)) ; return cub(_u-u(i))/S(i); }

    // integral from u(i) to u(i+1) -- note: these are constant one we know the knotvector, so can (should) be cached...
    double IA(int i) const{ return qua(u(i)-u(i+1))/(4*P(i)); }                            //
    double IB(int i) const{ return - cub(u(i)-u(i+1))*(3*u(i)-4*u(i-2)+u(i+1))/(12*P(i))     //
                                   - sqr(u(i)-u(i+1))*(3*sqr(u(i))-2*u(i-1)*u(i+1)+sqr(u(i+1))+u(i)*(-4*u(i-1)+2*u(i+1)-4*u(i+2)) +6*u(i-1)*u(i+2)-2*u(i+1)*u(i+2) )/(12*Q(i)) //
                                   + sqr(u(i)-u(i+1))*(3*sqr(u(i+1))+sqr(u(i  ))+2*u(i)*u(i+1)-8*u(i+1)*u(i+2)-4*u(i  )*u(i+2)+6*sqr(u(i+2)))/(12*R(i)); }  //
    double IC(int i) const{ return sqr(u(i)-u(i+1))*(3*sqr(u(i  ))+sqr(u(i+1))+2*u(i+1)*u(i)-8*u(i  )*u(i-1)-4*u(i-1)*u(i+1)+6*sqr(u(i-1)))/(12*Q(i))     //
                                 - sqr(u(i)-u(i+1))*(3*sqr(u(i+1))+sqr(u(i))-4*u(i-1)*u(i+1)+6*u(i-1)*u(i+2)-4*u(i+1)*u(i+2)-2*u(i)*(u(i-1)-u(i+1)+u(i+2)))/(12*R(i)) //
                                 + cub(u(i)-u(i+1))*(3*u(i+1)-4*u(i+3)+u(i))/(12*S(i)); }  //
    double ID(int i) const{ return qua(u(i)-u(i+1))/(4*S(i)); }                            //

private:
    double sqr(double x) const { return x*x; }
    double cub(double x) const { return x*x*x; }
    double qua(double x) const { return x*x*x*x; }
    // TODO: P,Q,R,S only depend on the knot vector, so build at construction, and cache them...
    double _P(int i) const { return  (_u[i+1]-_u[i-2])*(_u[i+1]-_u[i-1]) * (_u[i+1]-_u[i]); }
    double _Q(int i) const { return  (_u[i+1]-_u[i-1])*(_u[i+2]-_u[i-1]) * (_u[i+1]-_u[i]); }
    double _R(int i) const { return  (_u[i+2]-_u[i  ])*(_u[i+2]-_u[i-1]) * (_u[i+1]-_u[i]); }
    double _S(int i) const { return  (_u[i+2]-_u[i  ])*(_u[i+3]-_u[i  ]) * (_u[i+1]-_u[i]); }

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
Double_t RooCubicBSpline::evaluate() const 
{
  double u = _x;
  int i = _aux->i(u); // location in knot vector
  assert(i>2);

  return ((RooAbsReal*)_coefList.at(i-3))->getVal()*_aux->A(u,i)
       + ((RooAbsReal*)_coefList.at(i-2))->getVal()*_aux->B(u,i) 
       + ((RooAbsReal*)_coefList.at(i-1))->getVal()*_aux->C(u,i) 
       + ((RooAbsReal*)_coefList.at(i  ))->getVal()*_aux->D(u,i) ;
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
    norm += ((RooAbsReal*)_coefList.at(i  ))->getVal()*_aux->IA(i+3);
    norm += ((RooAbsReal*)_coefList.at(i+1))->getVal()*_aux->IB(i+3) ;
    norm += ((RooAbsReal*)_coefList.at(i+2))->getVal()*_aux->IC(i+3) ;
    norm += ((RooAbsReal*)_coefList.at(i+3))->getVal()*_aux->ID(i+3) ;
  }
  return norm;
}


#include <algorithm>
#include "P2VV/RooCubicSplineKnot.h"
#include "RooAbsReal.h"

namespace RooCubicSplineKnot_aux {

  Double_t get(const RooArgList& b,int i) { return ((RooAbsReal&)b[i]).getVal() ; }
  Double_t get(const RooArgList& b,int i,int k) { return RooCubicSplineKnot_aux::get(b,i+k); }

  template <typename T> typename T::const_reference get(const T& t, int i, int j) { return t[4*i+j]; }
  template <typename T> void push_back(T& t, const typename T::value_type& a,
                                             const typename T::value_type& b, 
                                             const typename T::value_type& c,
                                             const typename T::value_type& d) { t.push_back(a); t.push_back(b); t.push_back(c); t.push_back(d) ; }
}

// on input, y contains the values at the knot locations
// on output, it contains the b-spline coefficients 
// Note: one element will be pre-pended, and one post-pended !!
void RooCubicSplineKnot::computeCoefficients(std::vector<double>& y ) const 
{
 // see http://en.wikipedia.org/wiki/Spline_interpolation
 // for the derivation of the linear system...
 // see http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
 // for the O(N) algorithm to solve the relevant linear system...
    int n = _u.size();
    assert(y.size()==_u.size());

    double bf = y.front() / A(u(0),0) ;
    double bb = y.back()  / D(u(n-1),n-2);

    y.front() = - bf * double(6) / sqr(h(1,0));
    y.back()  = - bb * double(6) / sqr(h(n-1,n-2));

    std::vector<double> c ; c.reserve(n);

    c.push_back(mc(0) / mb(0) );
    y[0] /=  mb(0);
 
    for (int i = 1; i < n; ++i) {
        double m = double(1) / (mb(i) - ma(i) * c.back() ) ;
        c.push_back( mc(i) * m );
        y[i] = (y[i] - ma(i) * y[i - 1]) * m;
    }
    for (int i = n-1 ; i-- > 0; ) y[i] -= c[i] * y[i + 1];

    y.push_back(bb);
    y.insert(y.begin(),bf); // ouch... expensive!
}

void RooCubicSplineKnot::fillPQRS() const {
           assert(_PQRS.empty());
           // P,Q,R,S only depend on the knot vector, so build at construction, and cache them...
           _PQRS.reserve(4*_u.size());
           for (int i=0;i<_u.size();++i) { 
                _PQRS.push_back( h(i+1,i-2)*h(i+1,i-1)*h(i+1,i) );
                _PQRS.push_back( h(i+1,i-1)*h(i+2,i-1)*h(i+1,i) );
                _PQRS.push_back( h(i+2,i  )*h(i+2,i-1)*h(i+1,i) );
                _PQRS.push_back( h(i+2,i  )*h(i+3,i  )*h(i+1,i) );
           }

    }

double RooCubicSplineKnot::evaluate(double _u, const RooArgList& b) const {
        using RooCubicSplineKnot_aux::get;
        int i = index(_u); // location in knot vector
        assert(0<=i && i+3<b.getSize());
        //assert( u(i) <= _u && _u<= u(i+1) );
        return  get(b,i,0)*A(_u,i) // TODO: substitute A,B,C,D 'in situ'
             +  get(b,i,1)*B(_u,i)
             +  get(b,i,2)*C(_u,i)
             +  get(b,i,3)*D(_u,i);
}

double RooCubicSplineKnot::analyticalIntegral(const RooArgList& b) const {
        using RooCubicSplineKnot_aux::push_back;
        using RooCubicSplineKnot_aux::get;
        if (_IABCD.empty()) {
           // the integrals of A,B,C,D from u(i) to u(i+1) only depend on the knot vector...
           // so we create them 'on demand' and cache the result
           _IABCD.reserve(4*_u.size());
           for (int j=0;j<_u.size();++j) { 
               push_back(_IABCD,   qua(h(j,j+1))/(4*P(j)) 
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
            norm += get(b,i,k)*RooCubicSplineKnot_aux::get(_IABCD,i,k) ;
        }
        return norm;
}

int RooCubicSplineKnot::index(double u) const 
{ 
        // assert(u>=_u.front() && u<=_u.back());   
        std::vector<double>::const_iterator i = --std::upper_bound(_u.begin(),_u.end()-1,u);
        // assert( _u.begin()<=i );
        // assert( *i <= u && u<=*(i+1) );
        return std::max(std::distance(_u.begin(),i),0L);
};

// S matrix for i-th interval
RooCubicSplineKnot::S_jk RooCubicSplineKnot::S_jk_sum(int i, const RooArgList& b) const 
{
            if (_S_jk.empty()) {
                _S_jk.reserve(_u.size()*4);
                for(int i=0;i<_u.size();++i) {
                    _S_jk.push_back( -RooCubicSplineKnot::S_jk(u(i+1),u(i+1),u(i+1))/P(i) );
                    _S_jk.push_back(  RooCubicSplineKnot::S_jk(u(i-2),u(i+1),u(i+1))/P(i)
                                     +RooCubicSplineKnot::S_jk(u(i-1),u(i+1),u(i+2))/Q(i) 
                                     +RooCubicSplineKnot::S_jk(u(i  ),u(i+2),u(i+2))/R(i) );
                    _S_jk.push_back( -RooCubicSplineKnot::S_jk(u(i-1),u(i-1),u(i+1))/Q(i)
                                     -RooCubicSplineKnot::S_jk(u(i-1),u(i  ),u(i+2))/R(i) 
                                     -RooCubicSplineKnot::S_jk(u(i  ),u(i  ),u(i+3))/S(i) );
                    _S_jk.push_back(  RooCubicSplineKnot::S_jk(u(i  ),u(i  ),u(i  ))/S(i) );
                }
            }
            using RooCubicSplineKnot_aux::get;
            return get(_S_jk,i,0)*get(b,i,0)
                 + get(_S_jk,i,1)*get(b,i,1)
                 + get(_S_jk,i,2)*get(b,i,2)
                 + get(_S_jk,i,3)*get(b,i,3);
}


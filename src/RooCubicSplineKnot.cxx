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

void RooCubicSplineKnot::smooth(std::vector<double>& y, const std::vector<double>& dy, double lambda) const {
// implementation as described in D.S.G. Pollock, Smoothing Splines" 
// see http://r.789695.n4.nabble.com/file/n905996/SPLINES.PDF
    using namespace std;
    int n = y.size();
    vector<double> uu(n-2),vv(n-3),ww(n-4), q(n-2);
    assert( dy.size()==n);
    assert( size()==n);
    // lambda = 0 : no smoothing ; lambda -> 1: straight line (the ultimate smooth curve)
    if (lambda<0|| !(lambda<1)) {
        throw std::string("RooCubicSplineKnot::smooth: smoothing parameter must be in range [0,1)");
    }
    double mu =2*lambda/(3*(1-lambda));
    for (int i=0;i<n-2;++i) {
        q[i]  = r(i)*(y[i]-y[i+1])-r(i+1)*(y[i+1]-y[i+2]);
        uu[i] = p(i+1)+mu*(sqr(r(i)*dy[i])+sqr(f(i+1)*dy[i+1])+sqr(r(i+1)*dy[i+2]) );
        if (i>n-4) continue;
        vv[i] = h(i+1)+mu*(f(i+1)*r(i+1)*sqr(dy[i+1]) + f(i+1)*r(i+2)*sqr(dy[i+2]));
        if (i>n-5) continue ;
        ww[i] =        mu*r(i+1)*r(i+2)*sqr(dy[i+2]);
    }

    // Solve A b = q  with  A(i,i)=uu(i), A(i+1,i)=A(i,i+1) = vv(i), A(i+2,i),A(i,i+2)=ww(i)
    // start with factorization: A -> L D LT
    vv[0] /= uu[0]; ww[0] /= uu[0];
    uu[1] -= uu[0]*sqr(vv[0]);
    vv[1] -= uu[0]*vv[0]*ww[0];
    vv[1] /= uu[1]; ww[1] /= uu[1];
    for (int i=2;i<n-4;++i) {
        uu[i] -= uu[i-1]*sqr(vv[i-1])+uu[i-2]*sqr(ww[i-2]);
        vv[i] -= uu[i-1]*vv[i-1]*ww[i-1];  vv[i] /= uu[i];
        ww[i] /= uu[i];
    }
    uu[n-4] -= uu[n-5]*sqr(vv[n-5])+uu[n-6]*sqr(ww[n-6]);
    vv[n-4] -= uu[n-5]*vv[n-5]*ww[n-5];  vv[n-4] /= uu[n-4];
    uu[n-3] -= uu[n-4]*sqr(vv[n-4])+uu[n-5]*sqr(ww[n-5]);

    // forward substitution (i.e. solve L ( D LT  b )  = q  for D LT b -> q
    q[1]                         -= vv[0  ]*q[0  ];
    for (int i=2;i<n-2;++i) q[i] -= vv[i-1]*q[i-1]+ww[i-2]*q[i-2];
    // rescale  (i.e. solve D (LT b ) = q  for LT b
    for (int i=0;i<n-2;++i) { q[i] /= uu[i] ; }
    // backward substitution (i.e. solve LT b = q  for b -> q )
    q[n-4]                          -= vv[n-4]*q[n-3];
    for (int i=n-5;i>=0;--i) q[i  ] -= vv[i  ]*q[i+1]+ww[i]*q[i+2];
    // solve for y...
    int i=0;
    y[i]                 -=mu*sqr(dy[i])*(r(i)*q[i]                          ); ++i;
    y[i]                 -=mu*sqr(dy[i])*(r(i)*q[i]+f(i)*q[i-1]              ); ++i;
    for (;i<n-2;++i) y[i]-=mu*sqr(dy[i])*(r(i)*q[i]+f(i)*q[i-1]+r(i-1)*q[i-2]);
    y[i]                 -=mu*sqr(dy[i])*(          f(i)*q[i-1]+r(i-1)*q[i-2]); ++i;
    y[i]                 -=mu*sqr(dy[i])*(                      r(i-1)*q[i-2]);

}

// on input, y contains the values at the knot locations
// on output, it contains the b-spline coefficients
// Note: one element will be pre-pended, and one post-pended !!
void RooCubicSplineKnot::computeCoefficients(std::vector<double>& y) const
{
 // see http://en.wikipedia.org/wiki/Spline_interpolation
 // for the derivation of the linear system...
 // see http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
 // for the O(N) algorithm to solve the relevant linear system...
    int n = size();
    assert(y.size()==size());

    double bf = y.front() / A(u(0),0) ;
    double bb = y.back()  / D(u(n-1),n-2);

    y.front() = - bf * double(6) / sqr(h(1,0));
    y.back()  = - bb * double(6) / sqr(h(n-1,n-2));

    std::vector<double> c ; c.reserve(n);
    c.push_back( mc(0) / mb(0) );
    y[0] /=  mb(0);
    for (int i = 1; i < n; ++i) {
        double m = double(1) / (mb(i) - ma(i) * c.back() ) ;
        c.push_back( mc(i) * m );
        y[i] = (y[i] - ma(i) * y[i - 1]) * m;
    }
    for (int i = n-1 ; i-- > 0; ) y[i] -= c[i] * y[i + 1];
    y.push_back(bb); y.insert(y.begin(),bf); // ouch... expensive!
}




void RooCubicSplineKnot::fillPQRS() const {
    assert(_PQRS.empty());
    // P,Q,R,S only depend on the knot vector, so build at construction, and cache them...
    _PQRS.reserve(4*size());
    for (unsigned int i=0;i<size();++i) {
        _PQRS.push_back( h(i+1,i-2)*h(i+1,i-1)*h(i+1,i) );
        _PQRS.push_back( h(i+1,i-1)*h(i+2,i-1)*h(i+1,i) );
        _PQRS.push_back( h(i+2,i  )*h(i+2,i-1)*h(i+1,i) );
        _PQRS.push_back( h(i+2,i  )*h(i+3,i  )*h(i+1,i) );
    }

}

double RooCubicSplineKnot::evaluate(double x, const RooArgList& b) const {
    using RooCubicSplineKnot_aux::get;
    int i = index(x); // location in knot vector
    assert(-1<=i && i<=size());
    // we're a 'natural' spline. Extrapolate using the derivative at the first (last) knot
    if (i==-1) {
        return evaluate(u(0),b) - d(x,0)*r(0)*(get(b,0,0)-get(b,0,1));
    }
    if (i==size()) {
        i = size()-1;
        return evaluate(u(i),b) + d(x,i)*r(i-1)*(get(b,i,2)-get(b,i,1));
    }
    assert( u(i) <= x && x<= u(i+1) );
    return get(b,i,0)*A(x,i) // TODO: substitute A,B,C,D 'in situ'
         + get(b,i,1)*B(x,i)
         + get(b,i,2)*C(x,i)
         + get(b,i,3)*D(x,i);
}

double RooCubicSplineKnot::analyticalIntegral(const RooArgList& b) const {
    using RooCubicSplineKnot_aux::push_back;
    using RooCubicSplineKnot_aux::get;
    if (_IABCD.empty()) {
        // the integrals of A,B,C,D from u(i) to u(i+1) only depend on the knot vector...
        // so we create them 'on demand' and cache the result
        _IABCD.reserve(4*size());
        for (unsigned int j=0;j<size();++j) {
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
    assert(b.getSize()-2==size());
    double norm(0);
    for (unsigned int i=0; i < size()-1; ++i) for (int k=0;k<4;++k) {
        norm += get(b,i,k)*RooCubicSplineKnot_aux::get(_IABCD,i,k) ;
    }
    return norm;
}

int RooCubicSplineKnot::index(double u) const
{
    if (u>_u.back()) return size();
    std::vector<double>::const_iterator i = --std::upper_bound(_u.begin(),_u.end()-1,u);
    return std::distance(_u.begin(),i);
};

// S matrix for i-th interval
RooCubicSplineKnot::S_jk RooCubicSplineKnot::S_jk_sum(int i, const RooArgList& b) const
{
    if (_S_jk.empty()) {
        _S_jk.reserve(size()*4);
        for(int i=0;i<size();++i) {
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


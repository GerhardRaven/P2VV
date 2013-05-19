#include "TMath.h"
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

#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixTUtils.h"
#include "TDecompLU.h"

#if 0
void RooCubicSplineKnot::smooth(std::vector<double>& y, const std::vector<double>& dy, double lambda) const {
        unsigned n=y.size();
        assert(n==size());
        assert(n==dy.size());
        assert(n >= 2);
        TMatrixD D(n-2,n);
        for (unsigned int i=0;i<n-2;++i) {
            double ih0( double(1)/h(i)) , ih1( double(1)/h(i+1) );
            D(i,i) = ih0; D(i,i+1) = -ih0-ih1; D(i,i+2) = ih1;
        }
        TMatrixDSym W(n-2);
        for (unsigned int i=0;i<n-2;++i) {
            double h0 = h(i), h1 = h(i+1);
            W(i,i) = (h0+h1)/3;
            if (i>0) { W(i,i-1) = W(i-1,i) = h0/6; }
        }
        W.Invert(); W.SimilarityT(D); W *= lambda;
        for (unsigned int i=0;i<n;++i) W(i,i) += double(1)/sqr(dy[i]);
        TVectorD vals(n);
        for (unsigned int i=0;i<n;++i) vals(i)=y[i]/dy[i];
        TDecompLU lu(W); lu.Solve( vals );
        for (unsigned int i=0;i<n;++i) y[i]=vals(i);
}
#endif

void RooCubicSplineKnot::smooth(std::vector<double>& y, const std::vector<double>& dy, double lambda) const {
// implementation as described in D.S.G. Pollock, Smoothing Splines" 
// see http://r.789695.n4.nabble.com/file/n905996/SPLINES.PDF
    using namespace std;
    int n = y.size();
    vector<double> uu(n-2),vv(n-3),ww(n-4), q(n-2);
    assert( dy.size()==n);
    assert( size()==n);
    assert(lambda>0 && lambda-1<0.0001); // lambda = 0 : no smoothing ; lambda = 1: straight line (the ultimate smooth curve)
    double mu =2*lambda/(3*(1-lambda));
    for (int i=0;i<n-2;++i) {
        q[i]  = r(i)*(y[i]-y[i+1])-r(i+1)*(y[i+1]-y[i+2]);
        uu[i] = p(i+1)+mu*(sqr(r(i))*sqr(dy[i])+sqr(f(i+1))*sqr(dy[i+1])+sqr(r(i+1))*sqr(dy[i+2]) );
        if (i>n-4) continue;
        vv[i] = h(i+1)+mu*(f(i+1)*r(i+1)*sqr(dy[i+1]) + f(i+1)*r(i+2)*sqr(dy[i+2]));
        if (i>n-5) continue ;
        ww[i] = mu*r(i+1)*r(i+2)*sqr(dy[i+2]);
    }

#if 0
    TMatrixDSym _A(n-2);
    for (int i=0;i<n-2;++i) { _A(i,i)=uu[i]; }
    for (int i=0;i<n-3;++i) { _A(i,i+1)=_A(i+1,i)=vv[i]; }
    for (int i=0;i<n-4;++i) { _A(i,i+2)=_A(i+2,i)=ww[i]; }
    _A.Print();
    TVectorD _q(n-2) ; for (int i=0;i<n-2;++i) { _q(i)=q[i]; }
#endif

    // Solve A b = q:   q: [0,n-2]
    //
    // factorisation; A -> L D LT
    //   with  A symmetric bandmatrix with bandwidth 5
    //         D diagonal, L lower triangle (LT upper triangle) with unit diagonal
    // input:   A_i,i = u[i], A_i,i+1 = A_i+1,i = v[i], A_i+2,i = A_i,i+2 = w[i]
    // output:  D_i,i = u[i], L_i,i=1, L_i,i+1 = v[i], L_i,i+2 = w[i]
    vv[0] /= uu[0]; ww[0] /= uu[0];
    uu[1] -= uu[0]*sqr(vv[0]);
    vv[1] -= uu[0]*vv[0]*ww[0];
    vv[1] /= uu[1]; ww[1] /= uu[1];
    for (int i=2;i<n-4;++i) {
        uu[i] -= uu[i-1]*sqr(vv[i-1])+uu[i-2]*sqr(ww[i-2]); // i < n-2
        vv[i] -= uu[i-1]*vv[i-1]*ww[i-1];  vv[i] /= uu[i];   // i < n-3
        ww[i] /= uu[i];                                  // i < n-4
    }
    uu[n-4] -= uu[n-5]*sqr(vv[n-5])+uu[n-6]*sqr(ww[n-6]);
    vv[n-4] -= uu[n-5]*vv[n-5]*ww[n-5];  vv[n-4] /= uu[n-4];
    uu[n-3] -= uu[n-4]*sqr(vv[n-4])+uu[n-5]*sqr(ww[n-5]);

    // forward substitution -- solve for L x = q,  q <- x
    // q[0]  = q[0];
    q[1] -= q[0]*vv[0];
    for (int i=2;i<n-2;++i) { q[i] -= vv[i-1]*q[i-1]+ww[i-2]*q[i-2]; }
    // rescale with 1/D    -- i.e. solution of L D x = q , q <- x
    for (int i=0;i<n-2;++i) { q[i] /= uu[i] ; }

    // backward substitution -- i.e. solve LT b = q
    // q[n-2]  = q[n-2]
    q[n-4] -= vv[n-4]*q[n-3];
    for (int i=n-5;i>=0;--i) { q[i]-=vv[i]*q[i+1]+ww[i]*q[i+2] ; }


#if 0
    cout << "verify solution. Residuals of linear system: " << endl;
    TVectorD _b(n-2) ; for (int i=0;i<n-2;++i) { _b(i)=q[i]; }
    TVectorD _x( _A * _b ); _x-=_q;
    _x.Print();
#endif

    // and finally, solve for y!
    int i=0;
    y[i]                 -=mu*sqr(dy[i])*(r(i)*q[i]                          ); ++i;
    y[i]                 -=mu*sqr(dy[i])*(r(i)*q[i]+f(i)*q[i-1]              );
    for (i=2;i<n-2;++i) y[i]-=mu*sqr(dy[i])*(r(i)*q[i]+f(i)*q[i-1]+r(i-1)*q[i-2]);
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
        // we're a 'natural' spline. Extrapolate using the derivative at the first knot
        if (i==-1) {
                return evaluate(u(0),b)
                      +(u(0)-x)/(u(1)-u(0))*3*(get(b,0,0)-get(b,0,1));
        }
        if (i==size()) {
                i = size()-1;
                return evaluate(u(i),b)
                      +(x-u(i))/(u(i)-u(i-1))*3*(get(b,i,2)-get(b,i,1));
        }
        assert( u(i) <= x && x<= u(i+1) );
        return  get(b,i,0)*A(x,i) // TODO: substitute A,B,C,D 'in situ'
             +  get(b,i,1)*B(x,i)
             +  get(b,i,2)*C(x,i)
             +  get(b,i,3)*D(x,i);
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


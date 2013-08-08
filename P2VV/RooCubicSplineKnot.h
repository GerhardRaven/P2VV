#ifndef ROO_CUBICSPLINEKNOT
#define ROO_CUBICSPLINEKNOT
#include <vector>
#include <algorithm>
#include "RooArgList.h"
#include "TVectorD.h"
#include "TH1.h"
#include "TGraphErrors.h"


class RooCubicSplineKnot {
public:
    RooCubicSplineKnot(const double *array, int nEntries) : _u( array, array+nEntries) { }
    template <typename Iter> RooCubicSplineKnot(Iter begin, Iter end) : _u(begin,end) { }
    RooCubicSplineKnot(const std::vector<double>& knots) : _u(knots) { } // needed to interface with python...

    double u(int i) const { 
        assert(size());
        assert(i>-3);
        assert(i<size()+3); 
        return _u[std::min(std::max(0,i),size()-1)]; }
    int size() const { return _u.size(); }
    double evaluate(double _u, const RooArgList& b) const;
    double analyticalIntegral(const RooArgList& b) const;

    void computeCoefficients(std::vector<double>& y ) const ;
    void smooth(std::vector<double>& y, const std::vector<double>& dy, double lambda) const;

    const std::vector<double>& knots() const { return _u; }

    class S_jk {
    public:
        S_jk(double a, double b, double c) : t(a*b*c), d(0.5*(a*b+a*c+b*c) ), s( 0.25*(a+b+c) ), o(0.125) { }
        S_jk(const S_jk& other, double offset=0) : t(other.t), d(other.d), s(other.s), o(other.o) {
            if (!offset) return;
            t+=offset*(-2*d+offset*(4*s-offset*o*8));
            d+=offset*(-8*s+3*offset*o*8)/2;
            s-=offset*3*o*8/4;
        }
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
            if (j>k) std::swap(j,k);
            assert(0<=j&&j<4);
            assert(0<=k&&k<4-j); // note: for 4-j<=k<4 could return 0... but better not to invoke with those..
            switch(3*j+k) {
                case 0: return   -t;   // (0,0)
                case 1: return    d;   // (0,1),(1,0)
                case 2: return   -s;   // (0,2),(2,0)
                case 3: return    o;   // (0,3),(3,0)
                case 4: return -2*s;   // (1,1)
                case 5: return  3*o;   // (1,2),(2,1)
            }
            assert(1==0);
            return 0;
        }
    private:
        double t,d,s,o;
    };
    // S matrix for i-th interval
    RooCubicSplineKnot::S_jk S_jk_sum(int i, const RooArgList& b) const ;

    class S_edge {
    public:
        S_edge(double a, double b) : alpha(0.5*a), beta(b) {}
        S_edge(const S_edge& other, double offset=0) : alpha(other.alpha), beta(other.beta) {
            assert(offset==0); // TODO: fix me!
        }
        double operator()(int j, int k) const {
            assert(j==0||j==1);
            assert(0<=(j+k) && (j+k)<2);
            return ( j+k==0 ) ? beta : alpha ;
        }
    private:
       double alpha,beta;
    };

    RooCubicSplineKnot::S_edge S_jk_edge(bool left, const RooArgList& b) const;

    // return integrals over the i-th bin of the j-th basis spline . exp(-gamma x)
    // as matrix_ij
    double expIntegral(const TH1* hist, double gamma, TVectorD& coefficients, TMatrixD& covarianceMatrix) const;

private:
    int index(double _u) const;
    double A(double _u,int i) const{ return -cub(d(_u,i+1))/P(i); }
    double B(double _u,int i) const{ return  sqr(d(_u,i+1))*d(_u,i-2)/P(i) + d(_u,i-1)*d(_u,i+2)*d(_u,i+1)/Q(i) + d(_u,i  )*sqr(d(_u,i+2))/R(i); }
    double C(double _u,int i) const{ return -sqr(d(_u,i-1))*d(_u,i+1)/Q(i) - d(_u,i  )*d(_u,i+2)*d(_u,i-1)/R(i) - d(_u,i+3)*sqr(d(_u,i  ))/S(i); }
    double D(double _u,int i) const{ return  cub(d(_u,i  ))/S(i); }

    double ma( int i) const {  // subdiagonal
        return i==size()-1 ?  double(6)/(h(i,i-2)*h(i-1) ) : A(u(i),i);
    }
    double mb( int i) const {   // diagonal
        return i==0        ?  -(double(6)/h(0)+double(6)/h(2,0))/h(0)
             : i==size()-1 ?  -(double(6)/h(i-1)+double(6)/h(i,i-2))/h(i-1)
                           : B(u(i),i) ;
    }
    double mc( int i) const {  // superdiagonal
        return i==0        ?   double(6)/(h(2,0)*h(0)) : C(u(i),i);
    }

    double P(int i) const { if (_PQRS.empty()) fillPQRS(); assert(4*i  <int(_PQRS.size())); return  _PQRS[4*i  ]; }
    double Q(int i) const { if (_PQRS.empty()) fillPQRS(); assert(4*i+1<int(_PQRS.size())); return  _PQRS[4*i+1]; }
    double R(int i) const { if (_PQRS.empty()) fillPQRS(); assert(4*i+2<int(_PQRS.size())); return  _PQRS[4*i+2]; }
    double S(int i) const { if (_PQRS.empty()) fillPQRS(); assert(4*i+3<int(_PQRS.size())); return  _PQRS[4*i+3]; }

    void fillPQRS() const;

    static double sqr(double x) { return x*x; }
    static double cub(double x) { return x*sqr(x); }
    static double qua(double x) { return sqr(sqr(x)); }
    double d(double _u, int j) const { return _u-u(j); }
    double d(double _u, int i, int j, int k) const { return d(_u,i)*d(_u,j)*d(_u,k); }
    double h(int i, int j) const { return u(i)-u(j); }
    double h(int i) const { return h(i+1,i); }
    double r(int i) const { return double(3)/h(i); }
    double f(int i) const { return -r(i-1)-r(i); }
    double p(int i) const { return 2*h(i-1)+h(i); }

    std::vector<double> _u;
    mutable std::vector<double> _PQRS;                   //!
    mutable std::vector<double> _IABCD;                  //!
    mutable std::vector<RooCubicSplineKnot::S_jk> _S_jk; //!
};

#endif

#include <vector>
#include <algorithm>
#include "RooComplex.h"
#include "RooArgList.h"


class RooKnot {
public:
    template <typename Iter> RooKnot(Iter begin, Iter end) : _u(begin,end) {
           // P,Q,R,S only depend on the knot vector, so build at construction, and cache them...
           _PQRS.reserve(4*_u.size());
           for (int i=0;i<_u.size();++i) { 
                _PQRS.push_back( h(i+1,i-2)*h(i+1,i-1)*h(i+1,i) );
                _PQRS.push_back( h(i+1,i-1)*h(i+2,i-1)*h(i+1,i) );
                _PQRS.push_back( h(i+2,i  )*h(i+2,i-1)*h(i+1,i) );
                _PQRS.push_back( h(i+2,i  )*h(i+3,i  )*h(i+1,i) );
           }
    }

    double evaluate(double _u, const RooArgList& b) const;
    double analyticalIntegral(const RooArgList& b) const;
    // RooComplex analyticalIntegral(const RooComplex& z, const RooArgList& coef) const;
private:
    int index(double u) const;
    double A(double _u,int i) const{ return -cub(d(_u,i+1))/P(i); }
    double B(double _u,int i) const{ return  sqr(d(_u,i+1))*d(_u,i-2)/P(i) + d(_u,i-1)*d(_u,i+2)*d(_u,i+1)/Q(i) + d(_u,i  )*sqr(d(_u,i+2))/R(i); }
    double C(double _u,int i) const{ return -sqr(d(_u,i-1))*d(_u,i+1)/Q(i) - d(_u,i  )*d(_u,i+2)*d(_u,i-1)/R(i) - d(_u,i+3)*sqr(d(_u,i  ))/S(i); }
    double D(double _u,int i) const{ return  cub(d(_u,i  ))/S(i); }

    double P(int i) const { assert(4*i  <_PQRS.size()); return  _PQRS[4*i  ]; }
    double Q(int i) const { assert(4*i+1<_PQRS.size()); return  _PQRS[4*i+1]; }
    double R(int i) const { assert(4*i+2<_PQRS.size()); return  _PQRS[4*i+2]; }
    double S(int i) const { assert(4*i+3<_PQRS.size()); return  _PQRS[4*i+3]; }



    class S_jk { 
    public:
        S_jk(double a, double b, double c) : t(a*b*c), d( (a*b+a*c+b*c)/2 ), s( (a+b+c)/4 ), o(double(1)/8) {}
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

    // S matrix for i-th interval
    RooKnot::S_jk S_jk_sum(int i, const RooArgList& b) const ;
    double sqr(double x) const { return x*x; }
    double cub(double x) const { return x*sqr(x); }
    double qua(double x) const { return sqr(sqr(x)); }
    double d(double _u, int j) const { return _u-u(j); }
    double h(int i, int j) const { return u(i)-u(j); }
    double u(int i) const { assert(i>-3&&i<int(_u.size()+3)); return _u[std::min(std::max(0,i),int(_u.size()-1))]; }

    const   std::vector<double> _u;
    mutable std::vector<double> _PQRS;
    mutable std::vector<double> _IABCD;
    mutable std::vector<RooKnot::S_jk> _S_jk;
};


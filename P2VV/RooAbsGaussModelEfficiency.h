#ifndef ROO_ABSGAUSSMODELEFF
#define ROO_ABSGAUSSMODELEFF
// Mixin class for Piecewise Polynomial Efficiency
//   from this class, a constant, peicewise constant and spline should be derived...
// TODO: figure out how to add inefficuence as function of some state...

// class should provide interface, and provide some of the work common to all piecewise polynomial efficiencies...

#include <complex>
#include "TMath.h"
#include "RooAbsReal.h"

class RooAbsGaussModelEfficiency : public RooAbsReal {
public:
    //C++11: using RooAbsReal::RooAbsReal;
    RooAbsGaussModelEfficiency() : RooAbsReal() {}
    ~RooAbsGaussModelEfficiency() ;
    RooAbsGaussModelEfficiency(const char *name, const char *title, const char *unit= "")   : RooAbsReal(name,title, unit) {};
    RooAbsGaussModelEfficiency(const RooAbsGaussModelEfficiency& other, const char* name=0) : RooAbsReal(other,name) {};
    virtual std::complex<double> productAnalyticalIntegral(Double_t umin, Double_t umax, 
                                                           Double_t scale, Double_t offset,
                                                           const std::complex<double>& z) const = 0;
private:
  ClassDef(RooAbsGaussModelEfficiency,1) 
};

namespace RooGaussModelAcceptance {


  class N {
    public:
        N(double x, const std::complex<double>& z) ;
        std::complex<double> operator()(unsigned i) const { return _N[i]; }
    private:
        std::complex<double> _N[3];
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
                      default : return 2*(*this)(j-1,2)/sqrt(TMath::Pi());
              }
              case 2: switch(j) {
                      case 0 : return -1;
                      case 1 : return -2*_x;
                      case 2 : return -2*(2*_x*_x-1);
                      case 3 : return -4*_x*(2*_x*_x-3);
                      // TODO: add higher orders!!!! Go up until six...
                      default : assert(1==0); return 0;
          }   }
          assert(1==0);
          return 0;
      }  
  private : 
      double _x;
  };

  template <unsigned MaxOrder=4U> class M_n {
  public:
        M_n(double x, const std::complex<double>& z) 
        {
          L_jk l(x); N n(x,z);
          for (int i=0;i<MaxOrder;++i) _m[i] = n(0)*l(i,0) + n(1)*l(i,1) + n(2)*l(i,2);
        }
        const std::complex<double>& operator()(int i) const { assert(0<=i&&i<MaxOrder); return _m[i]; }
        M_n& operator-=(const M_n& other) { for(int i=0;i<MaxOrder;++i) _m[i]= _m[i]-other._m[i]; return *this; }
        M_n  operator- (const M_n& other) const { return M_n(*this)-=other; }
  private:
        std::complex<double> _m[MaxOrder];
  };

  class K_n {
  public:
      K_n(const std::complex<double>& z) : _zi( std::complex<double>(1,0)/z) {}
      std::complex<double> operator()(unsigned i) const {
          assert(0<=i&&i<=3);
          switch(i) {
              case 0 : return 0.5*_zi;
              case 1 : return 0.5*_zi*_zi;
              case 2 : return _zi*(_zi*_zi+1.0);
              case 3 : std::complex<double> _zi2 = _zi*_zi; return _zi2*(3.*_zi2+3.);
              // TODO: add higher orders!!!! Go up until six... (product of two cubic splines)
          }
          assert(1==0);
          return 0;
      }
  private :
      std::complex<double> _zi;        
  };

}
#endif

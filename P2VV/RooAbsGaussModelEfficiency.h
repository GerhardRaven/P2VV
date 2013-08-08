#ifndef ROO_ABSGAUSSMODELEFF
#define ROO_ABSGAUSSMODELEFF
// Mixin class for Piecewise Polynomial Efficiency
// TODO: figure out how to add inefficuence as function of some state...


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

  template <unsigned MaxOrder> class M_n {
        std::complex<double> _m[MaxOrder];
  public:
        M_n(double x, const std::complex<double>& z) ;
        const std::complex<double>& operator()(unsigned i) const { assert(i<MaxOrder); return _m[i]; }
        M_n& operator-=(const M_n& other) { for(unsigned i=0;i<MaxOrder;++i) _m[i]= _m[i]-other._m[i]; return *this; }
        M_n  operator- (const M_n& other) const { return M_n(*this)-=other; }
  };

  class K_n {
        std::complex<double> _zi;
  public:
        K_n(const std::complex<double>& z) : _zi( std::complex<double>(1,0)/z) {}
        std::complex<double> operator()(unsigned i) const;
  };

}
#endif

#include "P2VV/RooAbsGaussModelEfficiency.h"
#include "RooMath.h"

namespace {
    static const Double_t rootpi(sqrt(TMath::Pi())) ;
    std::complex<double> evalApprox(Double_t x, const std::complex<double>& z) {
      // compute exp(-x^2)cwerf(-i(z-x)), cwerf(z) = exp(-z^2)erfc(-iz)
      // use the approximation: erfc(z) = exp(-z*z)/(sqrt(pi)*z)
      // to explicitly cancel the divergent exp(y*y) behaviour of
      // CWERF for z = x + i y with large negative y
      static const std::complex<double> mi(0,-1);
      std::complex<double> zp  = mi*(z-x);
      std::complex<double> zsq = zp*zp;
      std::complex<double> v = -zsq - x*x;
      std::complex<double> iz(z.imag()+x,z.real()-x); // ???
      return exp(v)*(exp(zsq)/(iz*rootpi) + 1.)*2. ;
    }

    // Calculate exp(-x^2) cwerf(i(z-x)), taking care of numerical instabilities
    std::complex<double> eval(Double_t x, const std::complex<double>& z) {
      Double_t re = z.real()-x;
      return (re>-5.0) ? RooMath::faddeeva_fast(std::complex<double>(-z.imag(),re))*exp(-x*x) 
                       : evalApprox(x,z) ;
    }
}


//RooAbsGaussModelEfficiency::RooAbsGaussModelEfficiency(const char *name, const char *title, const char *unit) 
//    : RooAbsReal( name, title, unit ) 
//{}

RooGaussModelAcceptance::N::N(double x, const std::complex<double>& z) 
{
          _N[0] =  RooMath::erf(x);
          _N[1] =  exp(-x*x);
          _N[2] =  eval(x,z);
}


// explicitly instantiate some templates...
template class RooGaussModelAcceptance::M_n<0U>;
template class RooGaussModelAcceptance::M_n<1U>;
template class RooGaussModelAcceptance::M_n<2U>;
template class RooGaussModelAcceptance::M_n<3U>;
template class RooGaussModelAcceptance::M_n<4U>;
template class RooGaussModelAcceptance::M_n<5U>;
template class RooGaussModelAcceptance::M_n<6U>;

#include <iostream>
#include <cmath>

#include "TH1.h"

double sigmaFromFT( const TH1& h1, double dMs, std::ostream& out )
{
  double sum(0), sumdt(0), sumcos(0), sumx2(0) ;
  for(int i=1; i<=h1.GetNbinsX() ; ++i) {
    double c = h1.GetBinContent( i ) ;
    double x = h1.GetBinCenter( i ) ;
    double dt = h1.GetBinWidth( i ) ;
    sum += c ;
    sumdt += c * dt ;
    sumcos += c * cos( - dMs * x ) * dt ;
    sumx2 += c * x * x ;
  }
  
  double rms = sqrt(sumx2/sum) ;
  double D = sumcos / sumdt ;
  double sigma = sqrt( -2*log(D) ) / dMs ;

  out << sum << " " << sumdt << " " << sumcos << " " << sumx2 << std::endl;
  out << "RMS of input histogram: " << rms<< std::endl
      << "If distribution were Gaussian, dilution is: " << exp(-0.5*rms*rms*dMs*dMs) << std::endl 
      << "Dilution from FT: " << D << std::endl
      << "Corresponding Gaussian resolution: " << sigma << std::endl;

  return D;
}

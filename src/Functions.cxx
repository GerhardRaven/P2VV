#include <iostream>
#include <cmath>

#include "TH1.h"

double sigmaFromFT( const TH1& h1, const double dMs, const double dMsErr, std::ostream& out)
{
   double sum(0), sumdt(0), sumdterr(0), sumcos(0), sumx2(0), error(0) ;
   for(int i=1; i<=h1.GetNbinsX() ; ++i) {
      double c = h1.GetBinContent( i ) ;
      double e = h1.GetBinError( i ) ;
      double x = h1.GetBinCenter( i ) ;
      double dt = h1.GetBinWidth( i ) ;
      sum += c ;
      sumdt += c * dt ;
      sumdterr += e * e * dt * dt;

      double cosv = cos( - dMs * x );
      sumcos += c * cosv * dt ;
      sumx2 += c * x * x ;

      double coserr = fabs(sin(dMsErr)) * dMsErr;
      error += (c * c * coserr * coserr +  cosv * cosv * e * e) * dt * dt;
   }

   error = sqrt(error / (sumdt * sumdt) + sumdterr * sumcos * sumcos / pow(sumdt, 4));

   double rms = sqrt(sumx2/sum) ;
   double D = sumcos / sumdt ;
   double sigma = sqrt( -2*log(D) ) / dMs ;

   out << sum << " " << sumdt << " " << sumcos << " " << sumx2 << std::endl;
   out << "RMS of input histogram: " << rms<< std::endl
       << "If distribution were Gaussian, dilution is: " << exp(-0.5*rms*rms*dMs*dMs) << std::endl 
       << "Dilution from FT: " << D << " +- " << error << std::endl
       << "Corresponding Gaussian resolution: " << sigma << std::endl;

   return D;
}

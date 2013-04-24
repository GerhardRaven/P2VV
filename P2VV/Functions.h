// $Id: $
#ifndef FUNCTIONS_H 
#define FUNCTIONS_H 1

#include <iostream>
class TH1l;

double sigmaFromFT( const TH1& h1, const double dMs, const double dMsErr, std::ostream& out = std::cout );

#endif // FUNCTIONS_H

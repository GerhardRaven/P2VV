// $Id: $
#ifndef FUNCTIONS_H 
#define FUNCTIONS_H 1

#include <iostream>
class TH1l;

double sigmaFromFT( const TH1& h1, double dMs, std::ostream& out = std::cout );

#endif // FUNCTIONS_H

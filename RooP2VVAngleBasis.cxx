/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   GR, Gerhard Raven,   Nikhef & VU, Gerhard.Raven@nikhef.nl
 *                                                                           *
 * Copyright (c) 2010, Nikhef & VU. All rights reserved.
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// END_HTML
//

#include "RooFit.h"
#include "Riostream.h"
#include <sstream>

#include "RooP2VVAngleBasis.h"
#include "RooLegendre.h"
#include "RooSpHarmonic.h"
#include "RooConstVar.h"

ClassImp(RooP2VVAngleBasis)
;

//_____________________________________________________________________________
RooP2VVAngleBasis::RooP2VVAngleBasis()
{
}

//_____________________________________________________________________________
// require RooRealVar observables, so that we know they are independent...
RooP2VVAngleBasis::RooP2VVAngleBasis( const char *name, const char *title
                                    , RooAbsReal& cpsi, RooAbsReal& ctheta, RooAbsReal& phi
                                    , int i, int j, int l, int m, double c )
 : RooProduct(name, title,RooArgSet())
 , _i(i), _j(j), _l(l), _m(m)
{
  if (c!=1) {
    std::stringstream C;
    C << name << ( c<0 ? "_m" : "_" ) << ( c<0?-c:c ) ;
    _compRSet.addOwned(*new RooConstVar( C.str().c_str(), C.str().c_str(), c ) );
  }
  std::stringstream P,Y;
  P << name << "_P_" << i << ( j<0 ? "_m" : "_" )  << (j<0?-j:j) ;
  _compRSet.addOwned(*new RooLegendre(   P.str().c_str(), P.str().c_str(), cpsi,i,j) );
  Y << name << "_Y_" << l << ( m<0 ? "_m" : "_" )  << (m<0?-m:m) ;
  _compRSet.addOwned(*new RooSpHarmonic( Y.str().c_str(), Y.str().c_str(), ctheta,phi,l,m) );
}

//_____________________________________________________________________________
RooP2VVAngleBasis::RooP2VVAngleBasis(const RooP2VVAngleBasis& other, const char* name) 
    : RooProduct(other, name)
    , _i(other._i), _j(other._j), _l(other._l), _m(other._m)
{
}


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
#include <math.h>

#include "RooP2VVAngleBasis.h"
#include "utils.h"

ClassImp(RooP2VVAngleBasis)
;

//_____________________________________________________________________________
RooP2VVAngleBasis::RooP2VVAngleBasis()
{
}

//_____________________________________________________________________________
// require RooRealVar observables, so that we know they are independent...
RooP2VVAngleBasis::RooP2VVAngleBasis( const char *name, const char *title
                                    , RooRealVar& cpsi, RooRealVar& ctheta, RooRealVar& phi
                                    , int i, int j, int l, int m
                                    , double c);
 : RooProduct(name, title)
 , _i(i), _j(j), _k(k), _l(l)
{
  _compRSet.addOwned( RooConstVar( ..., ..., c )
  _compRSet.addOwned( RooLegendre( ..., ..., cpsi,i,j) );
  _compRSet.addOwned( RooSpHarmonic( ..., ..., ctheta,phi,l,m) );
}

//_____________________________________________________________________________
RooP2VVAngleBasis::RooP2VVAngleBasis(const RooP2VVAngleBasis& other, const char* name) 
    : RooProduct(other, name)
    , _i(other._i), _j(other._j), _k(other._k), _l(other._l)
{
}


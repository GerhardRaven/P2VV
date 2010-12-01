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
RooP2VVAngleBasis::RooP2VVAngleBasis( const char *name, const char *title
                                    , RooAbsReal& cpsi, RooAbsReal& ctheta, RooAbsReal& phi
                                    , int i, int j, int l, int m, double c )
 : RooProduct(name, title,RooArgSet())
 , _c(c)
 , _i(i), _j(j), _l(l), _m(m)
 , _prod(false)
{
  std::stringstream P,Y;
  P << name << "_P_" << i << ( j<0 ? "_m" : "_" )  << (j<0?-j:j) ;
  _compRSet.addOwned(*new RooLegendre(   P.str().c_str(), P.str().c_str(), cpsi,i,j) );
  Y << name << "_Y_" << l << ( m<0 ? "_m" : "_" )  << (m<0?-m:m) ;
  _compRSet.addOwned(*new RooSpHarmonic( Y.str().c_str(), Y.str().c_str(), ctheta,phi,l,m) );
  if (c!=1) {
    std::stringstream C;
    C << name << ( c<0 ? "_m" : "_" ) << ( c<0?-c:c ) ;
    _compRSet.addOwned(*new RooConstVar( C.str().c_str(), C.str().c_str(), c ) );
  }
}

RooP2VVAngleBasis::RooP2VVAngleBasis( const char *name, const char *title
                                    , RooAbsReal& cpsi, RooAbsReal& ctheta, RooAbsReal& phi
                                    , int i1, int j1, int l1, int m1
                                    , int i2, int j2, int l2, int m2
                                    , double c )
 : RooProduct(name, title,RooArgSet())
 , _c(c)
 , _i(i1), _j(j1), _l(l1), _m(m1)
 , _prod(true)
{
  if (c!=1) {
    std::stringstream C;
    C << name << ( c<0 ? "_m" : "_" ) << ( c<0?-c:c ) ;
    _compRSet.addOwned(*new RooConstVar( C.str().c_str(), C.str().c_str(), c ) );
  }
  std::stringstream P,Y;
  P << name << "_P_" << i1 << ( j1<0 ? "_m" : "_" )  << (j1<0?-j1:j1) 
            <<  "_"  << i2 << ( j2<0 ? "_m" : "_" )  << (j2<0?-j2:j2) ;
  _compRSet.addOwned(*new RooLegendre(   P.str().c_str(), P.str().c_str(), cpsi,i1,j1, i2, j2) );
  Y << name << "_Y_" << l1 << ( m1<0 ? "_m" : "_" )  << (m1<0?-m1:m1) 
            <<  "_"  << l2 << ( m2<0 ? "_m" : "_" )  << (m2<0?-m2:m2) ;
  _compRSet.addOwned(*new RooSpHarmonic( Y.str().c_str(), Y.str().c_str(), ctheta,phi,l1,m1,l2,m2) );
}
//_____________________________________________________________________________
RooP2VVAngleBasis::RooP2VVAngleBasis(const RooP2VVAngleBasis& other, const char* name) 
    : RooProduct(other, name)
    , _c(other._c)
    , _i(other._i), _j(other._j), _l(other._l), _m(other._m)
    , _prod(other._prod)
{
}

RooP2VVAngleBasis* 
RooP2VVAngleBasis::createProduct(int i, int j, int l, int m, double c) const {
      std::stringstream name; name << this->GetName() << "_x_" << i << "_" << j 
                                                      << "_" << l << ( m<0 ? "_m":"_" ) << (m<0?-m:m) 
                                                      << ( c<0 ? "_m" : "_" ) << ( c<0?-c:c ) ;
      // grab first function, dynamic_cast to RooLegendre, grab its observable...
      // yes, really bad hacking...
      _compRIter->Reset();
      RooLegendre *P = dynamic_cast<RooLegendre*>(_compRIter->Next());
      RooSpHarmonic *Y = dynamic_cast<RooSpHarmonic*>(_compRIter->Next());
      assert(P!=0);
      assert(Y!=0);
      RooArgSet* Po = P->getParameters((RooAbsData*) 0);
      RooArgSet* Yo = Y->getParameters((RooAbsData*) 0);
      assert(Po->getSize()==1);  
      assert(Yo->getSize()==2);  
      TIterator* iter = Po->createIterator();
      RooAbsReal *cpsi = dynamic_cast<RooAbsReal*>(iter->Next());
      assert(cpsi!=0);
      delete iter;
      iter = Yo->createIterator();
      RooAbsReal *ctheta = dynamic_cast<RooAbsReal*>(iter->Next());
      assert(ctheta!=0);
      RooAbsReal *phi = dynamic_cast<RooAbsReal*>(iter->Next());
      assert(phi!=0);
      delete iter;


      return (!_prod) ? new  RooP2VVAngleBasis( name.str().c_str(), name.str().c_str()
                                              , *cpsi, *ctheta, *phi
                                              , _i, _j, _l, _m
                                              ,  i,  j,  l,  m
                                              , _c * c ) : 0 ;
  }

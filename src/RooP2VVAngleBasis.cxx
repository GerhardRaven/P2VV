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

#if 0
#include <fenv.h>
#include <iostream>

static int
feenableexcept (unsigned int excepts)
{
  static fenv_t fenv;
  unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
               old_excepts;  // previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = fenv.__control & FE_ALL_EXCEPT;

  // unmask
  fenv.__control &= ~new_excepts;
  // fenv.__mxcsr   &= ~(new_excepts << 7);

  std::cout  << "calling fesetenv" << std::endl;
  return ( fesetenv (&fenv) ? -1 : old_excepts );
}


int enable()  {
          std::cout << "enabling FPE" << std::endl;
          feclearexcept(FE_ALL_EXCEPT); // remove any 'stale' exceptions before switching on trapping
                               // otherwise we immediately trigger an exception...
          std::cout << "enabling FPE II" << std::endl;
          return feenableexcept(FE_ALL_EXCEPT);
}
//static int enableFPE = enable();

#endif

#include "RooFit.h"
#include "Riostream.h"
#include <sstream>
#include <memory>

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
                                    , Int_t i, Int_t j, Int_t l, Int_t m, Double_t c )
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
                                    , Int_t i1, Int_t j1, Int_t l1, Int_t m1
                                    , Int_t i2, Int_t j2, Int_t l2, Int_t m2
                                    , Double_t c )
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
RooP2VVAngleBasis::createProduct(Int_t i, Int_t j, Int_t l, Int_t m, Double_t c) const {
      std::stringstream name; name << this->GetName() << "_x_" << i << "_" << j 
                                                      << "_" << l << ( m<0 ? "_m":"_" ) << (m<0?-m:m) 
                                                      << ( c<0 ? "_m" : "_" ) << ( c<0?-c:c ) ;
      // grab first function, dynamic_cast to RooLegendre, grab its observable...
      // yes, really bad hacking...
      _compRIter->Reset();
      RooLegendre *P = dynamic_cast<RooLegendre*>(_compRIter->Next());
      assert(P!=0);
      RooArgSet* Po = P->getParameters((RooAbsData*) 0);
      assert(Po->getSize()==1);  
      std::auto_ptr<TIterator> iter( Po->createIterator() );
      RooAbsReal *cpsi = dynamic_cast<RooAbsReal*>(iter->Next());
      assert(cpsi!=0);
      RooSpHarmonic *Y = dynamic_cast<RooSpHarmonic*>(_compRIter->Next());
      assert(Y!=0);
      RooArgSet* Yo = Y->getParameters((RooAbsData*) 0);
      assert(Yo->getSize()==2);  
      iter.reset( Yo->createIterator() );
      RooAbsReal *ctheta = dynamic_cast<RooAbsReal*>(iter->Next());
      assert(ctheta!=0);
      RooAbsReal *phi = dynamic_cast<RooAbsReal*>(iter->Next());
      assert(phi!=0);

      return (!_prod) ? new  RooP2VVAngleBasis( name.str().c_str(), name.str().c_str()
                                              , *cpsi, *ctheta, *phi
                                              , _i, _j, _l, _m
                                              ,  i,  j,  l,  m
                                              , _c * c ) : 0 ;
  }

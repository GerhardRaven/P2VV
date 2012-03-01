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
#include <string>
#include <sstream>
#include <memory>
#include <algorithm>

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
 : RooProduct(name, title,RooArgList())
 , _c(c)
 , _i(i), _j(j), _l(l), _m(m)
 , _prod(false)
{
  // TODO: can we try to share these amongst all RooP2VVAngleBasis which use the same cpsi,chteta,phi???
  // actually, just make sure the names are 'shared' if two objects are the same (equivalent), then 
  // rely on the workspace conflict resolution to make them shared when adding to the workspace and
  // retrieving it...
  std::stringstream P,Y;
  P << name << "_P" << i << ( j<0 ? "m" : "" )  << (j<0?-j:j) ;
  _compRSet.addOwned(*new RooLegendre(   P.str().c_str(), P.str().c_str(), cpsi,i,j) );
  Y << name << "_Y" << l << ( m<0 ? "m" : "" )  << (m<0?-m:m) ;
  _compRSet.addOwned(*new RooSpHarmonic( Y.str().c_str(), Y.str().c_str(), ctheta,phi,l,m) );
  if (c!=1) {
    std::stringstream C;
    C << ( c<0 ? "_m" : "_" ) << ( c<0?-c:c );
    std::string strC = C.str();
    replace(strC.begin(), strC.end(), '.', 'p');
    strC.insert(0, name);
    _compRSet.addOwned(*new RooConstVar( strC.c_str(), strC.c_str(), c ) );
  }
}

RooP2VVAngleBasis::RooP2VVAngleBasis( const char *name, const char *title
                                    , RooAbsReal& cpsi, RooAbsReal& ctheta, RooAbsReal& phi
                                    , Int_t i1, Int_t j1, Int_t l1, Int_t m1
                                    , Int_t i2, Int_t j2, Int_t l2, Int_t m2
                                    , Double_t c )
 : RooProduct(name, title,RooArgList())
 , _c(c)
 , _i(i1), _j(j1), _l(l1), _m(m1)
 , _prod(true)
{
  if (c!=1) {
    std::stringstream C;
    C << ( c<0 ? "_m" : "_" ) << ( c<0?-c:c );
    std::string strC = C.str();
    replace(strC.begin(), strC.end(), '.', 'p');
    strC.insert(0, name);
    _compRSet.addOwned(*new RooConstVar( strC.c_str(), strC.c_str(), c ) );
  }
  std::stringstream P,Y;
  P << name << "_P" << i1 << ( j1<0 ? "m" : "" )  << (j1<0?-j1:j1) 
            << i2 << ( j2<0 ? "m" : "" )  << (j2<0?-j2:j2) ;
  _compRSet.addOwned(*new RooLegendre(   P.str().c_str(), P.str().c_str(), cpsi,i1,j1, i2, j2) );
  Y << name << "_Y" << l1 << ( m1<0 ? "m" : "" )  << (m1<0?-m1:m1) 
            << l2 << ( m2<0 ? "m" : "" )  << (m2<0?-m2:m2) ;
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
      std::stringstream name;
      name << "_x_" << i << j << l << ( m<0 ? "m":"" ) << (m<0?-m:m)
           << ( c<0 ? "_m" : "_" ) << ( c<0?-c:c );
      std::string strName = name.str();
      replace(strName.begin(), strName.end(), '.', 'p');
      strName.insert(0, this->GetName());

      // grab first function, dynamic_cast to RooLegendre, grab its observable...
      // yes, really bad hacking...
      _compRIter->Reset();
      RooLegendre *P = dynamic_cast<RooLegendre*>(_compRIter->Next());
      assert(P!=0);
      RooArgSet* Po = P->getVariables();
      assert(Po->getSize()==1);  
      std::auto_ptr<TIterator> iter( Po->createIterator() );
      RooAbsReal *cpsi = dynamic_cast<RooAbsReal*>(iter->Next());
      assert(cpsi!=0);
      RooSpHarmonic *Y = dynamic_cast<RooSpHarmonic*>(_compRIter->Next());
      assert(Y!=0);
      RooArgSet* Yo = Y->getVariables();
      assert(Yo->getSize()==2);  
      iter.reset( Yo->createIterator() );
      RooAbsReal *ctheta = dynamic_cast<RooAbsReal*>(iter->Next());
      assert(ctheta!=0);
      RooAbsReal *phi = dynamic_cast<RooAbsReal*>(iter->Next());
      assert(phi!=0);

      return (!_prod) ? new  RooP2VVAngleBasis( strName.c_str(), strName.c_str()
                                              , *cpsi, *ctheta, *phi
                                              , _i, _j, _l, _m
                                              ,  i,  j,  l,  m
                                              , _c * c ) : 0 ;
}

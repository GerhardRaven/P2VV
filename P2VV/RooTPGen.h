/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * Authors:                                                                  *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_TPGEN
#define ROO_TPGEN

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "Roo1DTable.h"
#include "RooConstVar.h"
#include "RooRealVar.h"

class RooGaussian;

class RooTPGen : public RooAbsPdf {
public:

   // Constructors, assignment etc
   inline RooTPGen() { }
   RooTPGen(const char *name, const char *title, RooArgList& tps, 
            const Roo1DTable& nPV, RooAbsReal& z, RooAbsPdf& zPDF);
   RooTPGen(const RooTPGen& other, const char* name=0);
   virtual TObject* clone(const char* newname) const { return new RooTPGen(*this, newname); }
   virtual ~RooTPGen();

   virtual RooAbsGenContext* genContext(const RooArgSet &vars, const RooDataSet *prototype=0, 
                                        const RooArgSet* auxProto=0, Bool_t verbose= kFALSE) const ;

   virtual Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
   virtual void initGenerator(Int_t /*code*/);
   virtual void generateEvent(Int_t code);

  
   double distance() const { return _mean.getVal(); }
   void setDistance(const double d) { _mean.setVal(d); }

   double smearing() const { return _sigma.getVal(); }
   void setSmearing(const double s) { _sigma.setVal(s); }

   void setDebug(bool d) { _debug = d; }
   bool debug() { return _debug; }

protected:

   virtual Double_t evaluate() const;

   RooListProxy _tps;
   Roo1DTable   _nPV;
   RooRealProxy _z;
   RooRealProxy _zPDF;
   mutable int _zCode;

   // default distance = 0.1 mm (0.3ps * c_light)
   RooRealVar   _d;
   RooRealVar   _mean;
   RooRealVar   _sigma;
   RooGaussian* _gauss;
   mutable int  _gCode;
   double       _min;
   double       _max;
   bool         _debug;

   ClassDef(RooTPGen,1) // General decay function p.d.f 
};

#endif

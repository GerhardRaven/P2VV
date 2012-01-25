/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOCALIBRATEDDILUTION
#define ROOCALIBRATEDDILUTION

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooCalibratedDilution : public RooAbsReal {
public:
   RooCalibratedDilution() {} ; 
   RooCalibratedDilution(const char *name, const char *title, RooAbsReal& _p0,
                         RooAbsReal& _p1, RooAbsReal& _w, RooAbsReal& _w_average);
   RooCalibratedDilution(const RooCalibratedDilution& other, const char* name=0);
   virtual TObject* clone(const char* newname) const
   { 
      return new RooCalibratedDilution(*this, newname);
   }
   inline virtual ~RooCalibratedDilution() { }

protected:

   RooRealProxy p0;
   RooRealProxy p1;
   RooRealProxy w;
   RooRealProxy w_average;
  
   Double_t evaluate() const;

private:

   ClassDef(RooCalibratedDilution, 1)
};
 
#endif

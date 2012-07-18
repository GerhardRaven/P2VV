/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooAverage.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(RooAverage) 

RooAverage::RooAverage(const char *name, const char *title, 
                       const RooArgSet& vars)
   : RooAbsReal(name,title),
   _vars("average_vars", "average_vars", this)
{ 
   RooFIter it = vars.fwdIterator();
   RooAbsArg* var = 0;
   while ((var = it.next())) {
      _vars.add(*var);
   }
} 

RooAverage::RooAverage(const RooAverage& other, const char* name)
   : RooAbsReal(other,name), 
   _vars("average_vars", this, other._vars)
{ 
} 

Double_t RooAverage::evaluate() const 
{ 
   RooFIter it = _vars.fwdIterator();
   RooAbsReal* var = 0;
   double av = 0;
   while ((var = static_cast<RooAbsReal*>(it.next()))) {
      av += var->getVal();
   }
   return av / _vars.getSize();
} 




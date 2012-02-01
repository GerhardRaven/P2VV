/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 
#include "RooRealCategory.h" 
#include "RooAbsCategory.h" 

ClassImp(RooRealCategory) 

 RooRealCategory::RooRealCategory(const char *name, const char *title, RooAbsCategory& _c) :
   RooAbsReal(name,title), 
   c("_c","_c",this,_c)
 { 
 } 

 RooRealCategory::RooRealCategory(const RooRealCategory& other, const char* name) :  
   RooAbsReal(other,name), 
   c("_c",this,other.c)
 { 
 } 

 Double_t RooRealCategory::evaluate() const 
 { 
   return c ; 
 } 

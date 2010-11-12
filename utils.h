#ifndef H_UTILS
#define H_UTILS
#ifndef _CINT_
#include <iostream>
#include "RooWorkspace.h"
#include "RooGlobalFunc.h"
#include "RooArgList.h"
#include "RooAbsReal.h"
#include "RooProduct.h"
#endif

// pitty: 'asprintf' doesn't work in CINT...
template <typename X>
char *Format(const char* fmt,const X&x) {
   char *buf = new char[1024] ; // let's hope that this is enough....
   sprintf(buf,fmt, x );
   return buf;
}   
template <typename X, typename Y>
char *Format(const char* fmt,const X&x, const Y&y) {
   char *buf = new char[1024] ; // let's hope that this is enough....
   sprintf(buf,fmt, x, y );
   return buf;
}   

template <typename X, typename Y, typename Z>
char *Format(const char* fmt,const X&x, const Y&y, const Z&z) {
   char *buf = new char[1024] ; // let's hope that this is enough....
   sprintf(buf,fmt, x, y, z );
   return buf;
}   

template <typename T>
T& get(RooWorkspace& w, T*(RooWorkspace::*fun)(const char*)const, const char* name) {
    T* x = (w.*fun)(name);
    if (x==0) {
        cout << "FAILURE: could not retrieve '" << name << "' from workspace" << endl;
        assert(1==0);
    }
    return *x;
}

template <typename T> T&   get(RooWorkspace& w, const char *name) { return dynamic_cast<T&>(*w.obj(name)); }
template<> RooAbsReal&     get<RooAbsReal>(RooWorkspace& w, const char* name) { return get<RooAbsReal>(w,&RooWorkspace::function,name); }
template<> RooAbsArg&      get<RooAbsArg>(RooWorkspace& w, const char* name)   { return get<RooAbsArg>(w,&RooWorkspace::arg,name); } // arg(w,name); }
template<> RooRealVar&     get<RooRealVar>(RooWorkspace& w, const char* name)   { return get<RooRealVar>(w,&RooWorkspace::var, name); } // (w,name); }

template <typename T>
T& import(RooWorkspace& w, const T& r, const char *n=0) {
   if (n==0) { w.import(r); return get<T>(w,r.GetName()); }
   else      { w.import(r,RooFit::RenameVariable(r.GetName(),n)); return get<T>(w,n); }
}

RooAbsReal& product(RooWorkspace& w, const RooArgList& x) {
    const char *name(0);
    for (int i=0;i<x.getSize();++i) {
        name = name ? Format("%s_%s",name,x[i].GetName()) : x[i].GetName();
    }
    RooAbsReal *p = w.function(name);
    if (p==0) p = &import(w,RooProduct(name,name,x));
    return *p;
}

RooAbsReal& product(RooWorkspace& w, RooAbsReal& x, RooAbsReal& y) {
    return product( w, RooArgList(x,y) );
}

RooAbsReal& product(RooWorkspace& w, RooAbsReal& x, RooAbsReal& y, RooAbsReal& z) {
    return product( w, RooArgList(x,y,z) );
}

RooArgList product(RooWorkspace& w, const RooArgList& x, const RooArgList& y) {
     RooArgList z;
     for (int i=0;i<x.getSize();++i) for (int j=0;j<y.getSize();++j) {
         z.add(product(w,x[i],y[j]));
     }
     return z;
}

#endif

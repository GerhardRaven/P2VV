#ifndef H_UTILS
#define H_UTILS
#ifndef _CINT_
#include <iostream>
#include "RooWorkspace.h"
#include "RooGlobalFunc.h"
#include "RooArgList.h"
#include "RooAbsReal.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
char *_Format_( const char *fmt, ... ) ;

// pitty: 'asprintf' doesn't work in CINT...
template <typename X>
inline char *Format(const char* fmt,const X&x) { return _Format_(fmt,x); }   
template <typename X, typename Y>
inline char *Format(const char* fmt,const X&x, const Y&y) { return _Format_(fmt, x, y ); }   
template <typename X, typename Y, typename Z>
inline char *Format(const char* fmt,const X&x, const Y&y, const Z&z) { return _Format_(fmt,x,y,z); }

template <typename T>
inline T& get(RooWorkspace& w, T*(RooWorkspace::*fun)(const char*)const, const char* name) {
    T* x = (w.*fun)(name);
    if (x==0) {
        cout << "FAILURE: could not retrieve '" << name << "' from workspace" << endl;
        assert(1==0);
    }
    return *x;
}

template<typename T> inline T&   get(RooWorkspace& w, const char *name) { return dynamic_cast<T&>(*w.obj(name)); } // TODO: check if T inherits from RooAbsArg -- if so, use arg, else use obj...
template<> inline RooAbsReal&     get<RooAbsReal>(RooWorkspace& w, const char* name) { return get<RooAbsReal>(w,&RooWorkspace::function,name); }
template<> inline RooAbsArg&      get<RooAbsArg>(RooWorkspace& w, const char* name)   { return get<RooAbsArg>(w,&RooWorkspace::arg,name); } // arg(w,name); }
template<> inline RooRealVar&     get<RooRealVar>(RooWorkspace& w, const char* name)   { return get<RooRealVar>(w,&RooWorkspace::var, name); } // (w,name); }

template <typename T>
inline T& import(RooWorkspace& w, const T& r, const char *n=0) {
   if (n==0) { w.import(r); return get<T>(w,r.GetName()); }
   else      { w.import(r,RooFit::RenameVariable(r.GetName(),n)); return get<T>(w,n); }
}


#endif

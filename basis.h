#ifndef H_BASIS
#define H_BASIS
#include "utils.h"
#include "RooRealVar.h"
#include "RooLegendre.h"
#include "RooSpHarmonic.h"
#include "RooP2VVAngleBasis.h"
#include "RooResolutionModel.h"

class abasis {  //TODO: make this an RooAbsReal implementation, which forwards integrals, 
                //      supports maxVal for toy generation, 
                //      and avoids making a constant to put in front of the product
                //  OR  make this class return a dedicated RooAbsReal implementation which
                //      does the above, i.e. a RooAbsReal which represents C * P_ij * Y_lm
public:
    // TODO: add version which takes RooAbsArg for cpsi, ctheta, phi (and then checks they 
    //  are in the workspace w!)
    abasis(RooWorkspace &w, const char *cpsi, const char *ctheta, const char *phi) 
        : _w(w)
        , _cpsi( get<RooRealVar>(w,cpsi) )
        , _ctheta( get<RooRealVar>(w,ctheta) )
        , _phi( get<RooRealVar>(w,phi) )
    { }

    RooAbsReal& operator()(const char* label, int i, int j, int k, int l, double c)  {
         char *name = Format("%s_%d_%d",label,i,j);
         name = Format( l<0 ? "%s_%d_m%d" : "%s_%d_%d",name,k,l<0?-l:l);
         RooAbsReal *b = _w.function(name);
         if (b==0) b = &import(_w, RooP2VVAngleBasis(name,name,_cpsi,_ctheta,_phi,i,j,k,l,c));
         return *b;
    }
private:

    RooWorkspace &_w;
    RooAbsReal &_cpsi;
    RooAbsReal &_ctheta;
    RooAbsReal &_phi;
};

class tbasis {
public:
    tbasis(RooWorkspace& w, RooResolutionModel& r, RooFormulaVar& basis, const char *label) 
        : _w(w)
        , _basis(import(w, *r.convolution(&basis,&r.convVar()), label))
    {  
    }
    
    RooAbsArg& operator()(const char* t1, const char *t2=0, const char *t3=0, const char *t4=0) {
       RooArgList l;
       const char *name(0);
       name = add( name, t1, l);
       name = add( name, t2, l);
       name = add( name, t3, l);
       name = add( name, t4, l);
       l.add( _basis);
       return product(_w,l);
    }
private:
    const char *add( const char *name, const char *t, RooArgList& l) {
        if(t==0) return name;
        l.add( get<RooAbsArg>(_w, t) );
        return name ? Format("%s_%s",name,t) : t;
    }
    RooWorkspace& _w;
    RooAbsReal&   _basis;
};
#endif

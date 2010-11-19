#ifndef H_BASIS
#define H_BASIS
#include "utils.h"
#include "RooRealVar.h"
#include "RooP2VVAngleBasis.h"

class abasis {
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

#endif

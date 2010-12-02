#include "basis.h"
#include "utils.h"
#include "RooP2VVAngleBasis.h"

ClassImp(abasis)
;

abasis::abasis(RooWorkspace &w, const char *cpsi, const char *ctheta, const char *phi) 
        : _w(w)
        , _cpsi( get<RooRealVar>(w,cpsi) )
        , _ctheta( get<RooRealVar>(w,ctheta) )
        , _phi( get<RooRealVar>(w,phi) )
    { }

abasis::~abasis() {}

RooAbsReal& abasis::operator()(const char* label, int i, int j, int k, int l, double c)  {
         char *name = Format("%s_%d_%d",label,i,j);
         name = Format( l<0 ? "%s_%d_m%d" : "%s_%d_%d",name,k,l<0?-l:l);
         RooAbsReal *b = _w.function(name);
         if (b==0) b = &import(_w, RooP2VVAngleBasis(name,name,_cpsi,_ctheta,_phi,i,j,k,l,c));
         return *b;
    }

#ifndef H_BASIS
#define H_BASIS
#include "RooAbsReal.h"
#include "RooWorkspace.h"

class abasis : public TObject {
public:
    // TODO: add version which takes RooAbsArg for cpsi, ctheta, phi (and then checks they 
    //  are in the workspace w!)
    abasis(RooWorkspace &w, const char *cpsi, const char *ctheta, const char *phi) ;
    virtual ~abasis();

    RooAbsReal& operator()(const char* label, int i, int j, int k, int l, double c);
private:

    RooWorkspace &_w;
    RooAbsReal &_cpsi;
    RooAbsReal &_ctheta;
    RooAbsReal &_phi;
    ClassDef(abasis,1) // 
};

#endif

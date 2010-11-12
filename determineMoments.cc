#include "basis.h"
#include "utils.h"
#include "TFile.h"
#include "RooAbsData.h"
#include <stdlib.h>

class moment {
public:
    moment(const RooAbsReal& x, const RooAbsPdf& pdf, const RooArgSet& nset) : _x(x),_pdf(pdf), _nset(nset), _m0(0),_m1(0),_m2(0) {}

    void inc(bool accepted=true) {
        double x = _x.getVal()/_pdf.getVal(&_nset);
        // TODO: make a histogram of x... (two, one for accept, one for all)
        ++_m0;
        if (accepted) {
            _m1 += x;
            _m2 += x*x;
        }
    }

    ostream& print(ostream& os) const {
        double mu = _m1/_m0;
        double sig2 = _m2/_m0 - mu*mu;
        return os << "moment("<< _x.GetName() <<"_" << _pdf.GetName() << ") = " << mu << " +- " << sqrt(sig2/_m0) << " significance: " << mu/sqrt(sig2/_m0) << endl;
    }
private:
    const RooAbsReal& _x;
    const RooAbsPdf& _pdf;
    const RooArgSet&  _nset;
    double _m0,_m1,_m2;
};

class eps {
public:
     eps(const RooAbsReal& cpsi, const RooAbsReal& ctheta, const RooAbsReal& phi) 
         : _cpsi(cpsi),_ctheta(ctheta),_phi(phi) 
     {}

     bool operator()() { 
         // TODO: make it 2D ;-)
         if (_ctheta.getVal()<0) return true;
         return drand48()>(0.5+0.5*(_cpsi.getVal()));
     }
private:
    const RooAbsReal& _cpsi;
    const RooAbsReal& _ctheta;
    const RooAbsReal& _phi;
};


void determineMoments() {

   TFile *f = new TFile("p2vv.root") ;
   RooWorkspace* w = (RooWorkspace*) f->Get("w") ;
   RooAbsPdf* pdf = w->pdf("pdf") ;
   RooAbsData* data = w->data("pdfData") ;

   // create moments
   RooArgSet *pdfObs = pdf->getObservables( data->get() );
   abasis ab(*w,"cpsi","ctheta","phi");
   eps efficiency( get<RooAbsReal>(*w,"cpsi")
                 , get<RooAbsReal>(*w,"ctheta")
                 , get<RooAbsReal>(*w,"phi"));
             
   std::vector<moment*> moments;
   typedef std::vector<moment*>::iterator moments_iterator; 
   for (int i=0;i<4;++i) {
     for (int l=0;l<4;++l) {
        for (int m=0;m<=l;++m) {
            moments.push_back(new moment( ab("mom",i,0,l,m,double(2*i+1)/2 ), *pdf, *pdfObs ) );
        }
     }
   }

   data->Print("V"); 
   pdfObs->Print("V");
   // loop over all data
   for (int i=0;i<data->numEntries(); ++i) {
        const RooArgSet *args = data->get(i);
        *pdfObs  = *args;
        // apply some fake efficiency, and see how it affects the moments...
        bool accept = efficiency();
        for ( moments_iterator m = moments.begin(); m!=moments.end(); ++m) (*m)->inc(accept);
   }

   // and print the results:
   for ( moments_iterator m = moments.begin(); m!=moments.end(); ++m) (*m)->print(cout);
   
}

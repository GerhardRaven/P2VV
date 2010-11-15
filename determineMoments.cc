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
        return os << "moment("<< _x.GetName() <<"_" << _pdf.GetName() << ") = " << mu << " +- " << sqrt(sig2/(_m0-1)) << " significance: " << mu/sqrt(sig2/_m0) << endl;
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
         return true;
         // TODO: make it 2D ;-)
         //if (_ctheta.getVal()<0) return true;
         //return drand48()>(0.5+0.5*(_cpsi.getVal()));
     }
private:
    const RooAbsReal& _cpsi;
    const RooAbsReal& _ctheta;
    const RooAbsReal& _phi;
};


void determineMoments(const char* fname="p2vv.root", const char* pdfName = "pdf", const char* dataName = "pdfData", const char *workspaceName = "w") {

   TFile *f = new TFile(fname);
   RooWorkspace* w = (RooWorkspace*) f->Get(workspaceName) ;
   RooAbsPdf* pdf = w->pdf(pdfName) ;
   RooAbsData* data = w->data(dataName) ;

   // create moments -- for this, we need the PDF used to generate the data
   RooArgSet *allObs = pdf->getObservables( data->get() );
   RooArgSet *marginalObs = pdf->getObservables( data->get() );
   marginalObs->remove( w->argSet("cpsi,ctheta,phi") );
   // marginalize pdf over 'the rest' so we get the normalization of the moments right...
   RooAbsPdf *pdf_marginal = pdf->createProjection(*marginalObs);

   abasis ab(*w,"cpsi","ctheta","phi");
   eps efficiency( get<RooAbsReal>(*w,"cpsi")
                 , get<RooAbsReal>(*w,"ctheta")
                 , get<RooAbsReal>(*w,"phi"));
             
   std::vector<moment*> moments;
   typedef std::vector<moment*>::iterator moments_iterator; 
   for (int i=0;i<5;++i) {
     for (int l=0;l<5;++l) {
        for (int m=0;m<=l;++m) {
            moments.push_back(new moment( ab("mom",i,0,l,m,double(2*i+1)/2 ), *pdf_marginal, *allObs ) );
        }
     }
   }

   // loop over all data
   for (int i=0;i<data->numEntries(); ++i) {
        const RooArgSet *args = data->get(i);
        *allObs  = *args;
        // apply some fake efficiency, and see how it affects the moments...
        bool accept = efficiency();
        for ( moments_iterator m = moments.begin(); m!=moments.end(); ++m) (*m)->inc(accept);
   }

   // and print the results:
   for ( moments_iterator m = moments.begin(); m!=moments.end(); ++m) (*m)->print(cout);
}

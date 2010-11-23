#include "basis.h"
#include "utils.h"
#include "TFile.h"
#include "TCanvas.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooRealSumPdf.h"
#include "RooPlot.h"
#include <stdlib.h>

class IMoment {
    public:
          IMoment(RooAbsReal &basis, const char *name=0) : _basis(basis), _m0(0),_m1(0),_m2(0), _name(name ? name : _basis.GetName() ) {}
          virtual ~IMoment() {};
          virtual void inc(bool accepted = true) {
                double x = evaluate();
                // TODO: make a histogram of x... (two, one for accept, one for all)
                ++_m0;
                if (accepted) {
                    _m1 += x;
                    _m2 += x*x;
                }
            }
          virtual ostream& print(ostream& os) const {
                double mu = _m1/_m0;
                double sig2 = _m2/_m0 - mu*mu;
                return os << "moment("<< _name << ") = " << mu << " +- " << sqrt(sig2/(_m0-1)) << " significance: " << mu/sqrt(sig2/_m0) << endl;
            }
          virtual RooAbsReal& basis() { return _basis; }
          virtual double coefficient() const = 0;
          virtual double evaluate(  ) = 0;
    protected:
          RooAbsReal &_basis;
          double _m0,_m1,_m2;
          const char *_name;
};

class Moment : public IMoment {
public:
    Moment(RooAbsReal& x, double c=1) : IMoment(x), _c(c) {}
    double evaluate() { return _basis.getVal(); }
    double coefficient() const { return _c*_m1/_m0; }
private:
    double _c;
};

class EffMoment  : public IMoment{
public:
    EffMoment(RooAbsReal& x, const RooAbsPdf& pdf, const RooArgSet& nset) : IMoment(x,Format("%s_%s",x.GetName(),pdf.GetName())),_pdf(pdf), _nset(nset)  {}

    double evaluate() { return _basis.getVal()/_pdf.getVal(&_nset); }
    double coefficient() const { return _m1/_m0; }
private:
    const RooAbsPdf& _pdf;
    const RooArgSet&  _nset;
};

class eps {
public:
     eps(const RooAbsReal& cpsi, const RooAbsReal& ctheta, const RooAbsReal& phi) 
         : _cpsi(cpsi),_ctheta(ctheta),_phi(phi) 
     {}

     bool operator()() { 
         // return true;
         bool acc = true;
         acc = acc && drand48()>(0.5+0.5*_cpsi.getVal()) ;
         acc = acc && drand48()>(0.5-0.5*_ctheta.getVal());
         double phi = _phi.getVal()/(4*acos(0.));
         phi = -1+2*phi; // [-1,1]
         //return phi<0 || drand48()>phi;
         return acc && drand48()>(0.5+0.5*fabs(phi));
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
             
   std::vector<IMoment*> moments;
   typedef std::vector<IMoment*>::iterator moments_iterator; 
   for (int i=0;i<5;++i) {
     for (int l=0;l<8;++l) {
        for (int m=-l;m<=l;++m) {
            // if we want to write it as efficiency, i.e. eps_ijk * P_i * Y_jk * PDF then we need the marginal..
            // moments.push_back(new EffMoment( ab("mom",i,0,l,m,double(2*i+1)/2 ), *pdf_marginal, *allObs ) );
            // here we effectively just want to compute the Fourier coefficients...
            moments.push_back(new Moment( ab("mom",i,0,l,m,1.), double(2*i+1)/2  ) );
        }
     }
   }
   // create some malformed input data...
   RooDataSet inEffData( "inEffData","inEffData", *allObs );
   for (int i=0;i<data->numEntries(); ++i) {
        const RooArgSet *args = data->get(i);
        *allObs  = *args;
       if (efficiency()) inEffData.add( *allObs );
   }
   //
   data = &inEffData;

   // loop over all data, determine moments
   for (int i=0;i<data->numEntries()/10; ++i) {
        const RooArgSet *args = data->get(i);
        *allObs  = *args;
        // apply some fake efficiency, and see how it affects the moments...
        bool accept = true; // efficiency();
        // if (!efficiency()) continue;
        for ( moments_iterator m = moments.begin(); m!=moments.end(); ++m) (*m)->inc(accept);
   }

   // and print the results:
   for ( moments_iterator m = moments.begin(); m!=moments.end(); ++m) (*m)->print(cout);

   // create a PDF from the moments
   RooArgList coef,fact;
   for ( moments_iterator m = moments.begin(); m!=moments.end(); ++m) {
       const char *name = Format("C_%f",(*m)->coefficient());
       w->factory(Format("%s[%f]",name,(*m)->coefficient()));
       coef.add(get<RooAbsReal>(*w,name));
       fact.add((*m)->basis());
   }
   RooAbsPdf *momsum = new RooRealSumPdf("pdf_mom","pdf_mom",fact,coef);

   const char *cvar[] = { "cpsi","ctheta","phi" };
   TCanvas *c = new TCanvas();
   c->Divide(3,1);
   for (int i = 0; i<3; ++i) {
        c->cd(1+i); RooPlot *plot = get<RooRealVar>(*w,cvar[i]).frame(); data->plotOn(plot); momsum->plotOn(plot); plot->Draw();
   }
}

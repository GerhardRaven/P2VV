#include "basis.h"
#include "utils.h"
#include "TFile.h"
#include "TCanvas.h"
#include "RooAbsArg.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooAddition_.h"
#include "RooCustomizer.h"
#include "RooPlot.h"
#include "RooGlobalFunc.h"
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
                return os << "moment("<< _name << ") = " << mu << " +- " << sqrt(sig2/(_m0-1)) << " significance: " << significance() << endl;
            }
          virtual RooAbsReal& basis() { return _basis; }
          virtual double coefficient() const = 0;
          virtual double significance() const { 
                double mu = _m1/_m0;
                double sig2 = _m2/_m0 - mu*mu;
                return fabs(mu/sqrt(sig2/(_m0-1)));
          }
          virtual double evaluate(  ) = 0;
    protected:
          RooAbsReal &_basis;
          double _m0,_m1,_m2;
          const char *_name;
};


class Moment : public IMoment {
public:
    Moment(RooAbsReal& x, double norm=1) : IMoment(x), _norm(norm) {}
    double evaluate() { return _basis.getVal(); }
    double coefficient() const { return _norm*_m1/_m0; }
private:
    double _norm;
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


void determineEfficiency(const char* fname="p2vv_4.root", const char* pdfName = "pdf", const char* dataName = "pdfData", const char *workspaceName = "w") {

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
   for (int i=0;i<4;++i) {
     for (int l=0;l<4;++l) {
        for (int m=-l;m<=l;++m) {
            // if we want to write it as efficiency, i.e. eps_ijk * P_i * Y_jk * PDF then we need the marginal..
            moments.push_back(new EffMoment( ab("mom",i,0,l,m,double(2*i+1)/2 ), *pdf_marginal, *allObs ) );
            // here we effectively just want to compute the Fourier coefficients...
            // Warning: the Y_lm are orthonormal, but the P_i are orthogonal, but the dot product is (2*i+1)/2
            // moments.push_back(new Moment( ab("mom",i,0,l,m,1.), double(2*i+1)/2  ) );
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
   // data = &inEffData;

   // loop over all data, determine moments
   for (int i=0;i<data->numEntries(); ++i) {
        const RooArgSet *args = data->get(i);
        *allObs  = *args;
        // apply some fake efficiency, and see how it affects the moments...
        bool accept = true; // efficiency();
        for ( moments_iterator m = moments.begin(); m!=moments.end(); ++m) (*m)->inc(accept);
   }

   // now we need to multiply all components of "pdf" with the efficiency...
   RooCustomizer customizer(*pdf,"_eff_");
   RooArgSet *comp = pdf->getComponents();
   TIterator *iter = comp->createIterator();
   RooAbsArg *j(0);
   while ( ( j = (RooAbsArg*)iter->Next() )!=0  ) {
        RooP2VVAngleBasis *orig = dynamic_cast<RooP2VVAngleBasis*>(j);
        if (orig==0) continue;
        TString name( orig->GetName() );
        cout << " modifying component " << name  << endl;
        name.Append("_eff");
        RooArgList s;
        for ( moments_iterator m = moments.begin(); m!=moments.end(); ++m) {
            // if (fabs((*m)->significance())<2) continue; // should _always_ use at least those moments which appear in signal pdf...
            s.add( *orig->createProduct(dynamic_cast<RooP2VVAngleBasis&>( (*m)->basis() ), (*m)->coefficient()));
        }
        
        RooAbsArg& rep = import(*w,RooAddition_( name, name, s, kTRUE)); // hand over ownership...
        cout << " new component " << rep.GetName() << endl;
        rep.Print();
        customizer.replaceArg( *orig, rep );
   }

   RooAbsPdf *pdf_eff = (RooAbsPdf*) customizer.build(kTRUE);

   const char *cvar[] = { "cpsi","ctheta","phi" };
   TCanvas *c = new TCanvas();
   c->Divide(3,1);
   for (int i = 0; i<3; ++i) {
        c->cd(1+i); 
        RooPlot *plot = get<RooRealVar>(*w,cvar[i]).frame(); 
        data->plotOn(plot); 
        pdf->plotOn(plot,RooFit::LineColor(kBlue)); 
        pdf_eff->plotOn(plot,RooFit::LineColor(kRed)); 
        plot->Draw();
   }
}

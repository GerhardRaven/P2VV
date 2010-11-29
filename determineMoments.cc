#include "basis.h"
#include "moments.h"
#include "utils.h"
#include "TFile.h"
#include "TCanvas.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooRealSumPdf.h"
#include "RooPlot.h"
#include <stdlib.h>

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


void determineMoments(const char* fname="p2vv_3.root", const char* dataName = "pdfData", const char *workspaceName = "w") {

   TFile *f = new TFile(fname);
   RooWorkspace* w = (RooWorkspace*) f->Get(workspaceName) ;
   RooAbsData* data = w->data(dataName) ;
   RooArgSet allObs = w->argSet("cpsi,ctheta,phi"); 

   abasis ab(*w,"cpsi","ctheta","phi");
   eps efficiency( get<RooAbsReal>(*w,"cpsi")
                 , get<RooAbsReal>(*w,"ctheta")
                 , get<RooAbsReal>(*w,"phi"));
             
   std::vector<IMoment*> moments;
   typedef std::vector<IMoment*>::iterator moments_iterator; 
   for (int i=0;i<4;++i) {
     for (int l=0;l<6;++l) {
        for (int m=-l;m<=l;++m) {
            // Warning: the Y_lm are orthonormal, but the P_i are orthogonal, but the dot product is (2*i+1)/2
            moments.push_back(new Moment( ab("mom",i,0,l,m,1.), double(2*i+1)/2  ) );
        }
     }
   }

   // create some malformed input data...
   RooDataSet inEffData( "inEffData","inEffData", allObs );
   for (int i=0;i<data->numEntries(); ++i) {
        allObs  = *data->get(i);
       if (efficiency()) inEffData.add( allObs );
   }
   //
   data = &inEffData;

   // loop over all data, determine moments
   for (int i=0;i<data->numEntries(); ++i) {
        allObs  = *data->get(i);
        for ( moments_iterator m = moments.begin(); m!=moments.end(); ++m) (*m)->inc();
   }

   // create a PDF from the moments
   RooArgList coef,fact;
   for ( moments_iterator m = moments.begin(); m!=moments.end(); ++m) {
       (*m)->print(cout);
       // if (fabs((*m)->significance())<2) continue; // should _always_ use at least those moments which appear in signal pdf...
       const char *name = Format("C_%f",(*m)->coefficient());
       w->factory(Format("%s[%f]",name,(*m)->coefficient()));
       coef.add(get<RooAbsReal>(*w,name));
       fact.add((*m)->basis());
   }
   cout << "using " << fact.getSize() << " moments out of " << moments.size() << endl;
   RooAbsPdf *momsum = new RooRealSumPdf("pdf_mom","pdf_mom",fact,coef);

   const char *cvar[] = { "cpsi","ctheta","phi" };
   TCanvas *c = new TCanvas();
   c->Divide(3,1);
   for (int i = 0; i<3; ++i) {
        c->cd(1+i); 
        RooPlot *plot = get<RooRealVar>(*w,cvar[i]).frame(); 
        data->plotOn(plot); 
        momsum->plotOn(plot); 
        plot->Draw();
   }
}

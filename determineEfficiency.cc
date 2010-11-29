#include "basis.h"
#include "utils.h"
#include "moments.h"
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
         acc = acc && drand48()>(0.5+0.5*fabs(phi));
         return acc;
     }
private:
    const RooAbsReal& _cpsi;
    const RooAbsReal& _ctheta;
    const RooAbsReal& _phi;
};


void determineEfficiency(const char* fname="p2vv_3.root", const char* pdfName = "pdf", const char* dataName = "pdfData", const char *workspaceName = "w") {

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
            // Warning: the Y_lm are orthonormal, but the P_i are orthogonal, but the dot product is (2*i+1)/2
            moments.push_back(new EffMoment( ab("mom",i,0,l,m,1. ),double(2*i+1)/2, *pdf_marginal, *allObs ) );
        }
     }
   }

   // create some inefficient data...
   RooDataSet inEffData( "inEffData","inEffData", *allObs );
   for (int i=0;i<data->numEntries(); ++i) {
        *allObs  = *data->get(i);
       if (efficiency()) inEffData.add( *allObs );
   }
   // replace input by inefficient data
   data = &inEffData;

   // loop over all data, determine moments
   for (int i=0;i<data->numEntries(); ++i) {
        *allObs  = *data->get(i);
        for ( moments_iterator m = moments.begin(); m!=moments.end(); ++m) (*m)->inc();
   }

   // now we need to multiply all relevant components (i.e. all RooP2VVAngleBasis ones) 
   // of "pdf" with their efficiency corrected versions...
   RooCustomizer customizer(*pdf,"_eff");
   RooArgSet *comp = pdf->getComponents();
   TIterator *iter = comp->createIterator();
   RooAbsArg *j(0);
   while ( ( j = (RooAbsArg*)iter->Next() )!=0  ) {
        RooP2VVAngleBasis *orig = dynamic_cast<RooP2VVAngleBasis*>(j);
        if (orig==0) continue;
        TString name( orig->GetName() ); name.Append("_eff");
        RooArgList s;
        for ( moments_iterator m = moments.begin(); m!=moments.end(); ++m) {
            s.add( *orig->createProduct(dynamic_cast<RooP2VVAngleBasis&>( (*m)->basis() ), (*m)->coefficient()));
        }
        RooAbsArg& rep = import(*w,RooAddition_( name, name, s, kTRUE)); // hand over ownership & put in workspace...
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

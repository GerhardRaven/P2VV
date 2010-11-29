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

// create a PDF from the moments of some reference dataset
// TODO: allow use of sWeights when computing moments
typedef std::vector<IMoment*>::iterator IMomentsIterator; 
RooAbsPdf *createMomPdf(RooWorkspace &w, const char *name, RooAbsData& data, IMomentsIterator begin, IMomentsIterator end) {
   if (begin==end) return 0; // no moments, no joy
   RooArgSet *allObs = (*begin)->basis().getObservables(data);
   for (int i=0;i<data.numEntries(); ++i) {
        *allObs  = *data.get(i);
        for ( IMomentsIterator m = begin; m!=end; ++m) (*m)->inc();
   }

   RooArgList coef,fact;
   for ( IMomentsIterator m = begin; m!=end; ++m) {
       // if (fabs((*m)->significance())<2) continue; // should _always_ use at least those moments which appear in signal pdf...
       const char *name = Format("C_%f",(*m)->coefficient());
       w.factory(Format("%s[%f]",name,(*m)->coefficient()));
       coef.add(get<RooAbsReal>(w,name));
       fact.add((*m)->basis());
   }
   return new RooRealSumPdf(name,name,fact,coef); // TODO: import into workspace...
}


void determineMoments(const char* fname="p2vv_3.root", const char* dataName = "pdfData", const char *workspaceName = "w") {

   TFile *f = new TFile(fname);
   RooWorkspace* w = (RooWorkspace*) f->Get(workspaceName) ;
   RooAbsData* data = w->data(dataName) ;
   RooArgSet allObs = w->argSet("cpsi,ctheta,phi"); 

   // create some malformed input data...
   eps efficiency( get<RooAbsReal>(*w,"cpsi")
                 , get<RooAbsReal>(*w,"ctheta")
                 , get<RooAbsReal>(*w,"phi"));
   RooDataSet inEffData( "inEffData","inEffData", allObs );
   for (int i=0;i<data->numEntries(); ++i) {
        allObs  = *data->get(i);
       if (efficiency()) inEffData.add( allObs );
   }
   // replace data with inefficient data...
   data = &inEffData;

   // decide with moments to take into account
   abasis ab(*w,"cpsi","ctheta","phi");
   std::vector<IMoment*> moments;
   for (int i=0;i<4;++i) {
     for (int l=0;l<6;++l) {
        for (int m=-l;m<=l;++m) {
            // Warning: the Y_lm are orthonormal, but the P_i are orthogonal, with dot product 2/(2*i+1)
            moments.push_back(new Moment( ab("mom",i,0,l,m,1.), double(2*i+1)/2  ) );
        }
     }
   }
   // and build a PDF using these moments...
   RooAbsPdf *pdf = createMomPdf(*w,"bkg_pdf",*data,moments.begin(),moments.end());

   // finally, make some plots to see if it worked...
   const char *cvar[] = { "cpsi","ctheta","phi" };
   TCanvas *c = new TCanvas();
   c->Divide(3,1);
   for (int i = 0; i<3; ++i) {
        c->cd(1+i); 
        RooPlot *plot = get<RooRealVar>(*w,cvar[i]).frame(); 
        data->plotOn(plot); 
        pdf->plotOn(plot); 
        plot->Draw();
   }
}

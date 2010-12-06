#ifndef _CINT_
#include "utils.h"
#include "basis.h"
#include "RooAbsArg.h"
#include "RooAddition_.h"
#include "RooBDecay.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooFormulaVar.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooMsgService.h"
#include "TCanvas.h"
#endif
using namespace std;

RooAbsPdf& _jpsikstar(RooWorkspace& w, const char* name ) 
{ 
        abasis ab(w, "cpsi", "ctheta", "phi");  // bind workspace and observables -- todo: use workspace hooks instead!
        // definition of the angular part of the PDF in terms of basis functions... 
        import(w, RooAddition_("AzAz_basis",      "AzAz_basis",       RooArgSet( ab("AzAz",       0,0,0, 0, 2.), ab("AzAz",       0,0,2,0, sqrt(1./5.)),  ab("AzAz",     0,0,2,2, -sqrt( 3./5.)) , 
                                                                                 ab("AzAz",       2,0,0, 0, 4.), ab("AzAz",       2,0,2,0, sqrt(4./5.)),  ab("AzAz",     2,0,2,2, -sqrt(12./5.)) )));
        import(w, RooAddition_("AparApar_basis",  "AparApar_basis",   RooArgSet( ab("AparApar",   2,2,0, 0, 1.), ab("AparApar",   2,2,2,0, sqrt(1./20.)),  ab("AparApar", 2,2,2,2,  sqrt( 3./20.)) ))); 
        import(w, RooAddition_("AperpAperp_basis","AperpAperp_basis", RooArgSet( ab("AperpAperp", 2,2,0, 0, 1.), ab("AperpAperp", 2,2,2,0,-sqrt(1./5.)))));
        import(w, RooAddition_("AparAperp_basis", "AparAperp_basis",  RooArgSet( ab("AparAperp",  2,2,2,-1, sqrt(3./5.)) )));
        import(w, RooAddition_("AzAperp_basis",   "AzAperp_basis",    RooArgSet( ab("AzAperp",    2,1,2, 1, sqrt(6./5.)) )));
        import(w, RooAddition_("AzApar_basis",    "AzApar_basis",     RooArgSet( ab("AzApar",     2,1,2,-2,-sqrt(6./5.)) )));
        //                                                    0    1    2       3       4      5 
        w.factory("expr::NAzAz      ('( @0 * @0 + @1 * @1 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})");
        w.factory("expr::NAparApar  ('( @4 * @4 + @5 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})");
        w.factory("expr::NAperpAperp('( @2 * @2 + @3 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})");
        w.factory("expr::ReAparAperp('( @4 * @2 + @5 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})");
        w.factory("expr::ReAzAperp  ('( @0 * @2 + @1 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})");
        w.factory("expr::ReAzApar   ('( @0 * @4 + @1 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})");
        w.factory("expr::ImAparAperp('( @4 * @3 - @5 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})");
        w.factory("expr::ImAzAperp  ('( @0 * @3 - @1 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})");

        import(w, RooFormulaVar("qmix","@0*@1",RooArgSet( get<RooCategory>(w,"qrec"),get<RooCategory>(w,"qtag") ) ) );

        w.factory("{Minus[-1],Zero[0],One[1]}");
        // In case of J/psi K*, things factorize...
        w.factory(Format("PROD::%s( RealSumPdf( { AzAz_basis , AparApar_basis, AperpAperp_basis, AparAperp_basis, AzAperp_basis, AzApar_basis} "
                                             ", { NAzAz,       NAparApar,      NAperpAperp,      ImAparAperp,     ImAzAperp,     ReAzApar } )"
                                " , BDecay(t,tau,Zero,One,Zero,qmix,Zero,dm,res,SingleSided))",name));
        return get<RooAbsPdf>(w,name);                                          
};
                                                                      
void jpsikstar() {                                                    
    RooWorkspace w("w");                                              
    // observables...
    w.factory(Format( "{ cpsi[-1,1], ctheta[-1,1], phi[0,%f], t[-1,4], qtag[bbar=+1,b=-1], qrec[jpsikstar=1,jpsikstarbar=-1} ",4*acos(0.))); // bbar=+1, so code corresponds to Bs(t=0)

    // choice: either fit for the Re&Im of the 3 amplitudes (and then
    //         constrain one phase and the sum of magnitudes)
    //         or fit in terms of angles and relative magnitudes
    // Note: initial values from arXiv:0704.0522v2 [hep-ex] BaBar PUB-07-009
    w.factory("{rz[0.556],rpar[0.211,0.,2.],rperp[0.233,0.,2.]}");
    w.factory("{deltaz[0],deltapar[-2.93,-3.15,3.15],deltaperp[2.91,-3.15,3.15]}");
    w.factory("expr::ReAz   ('rz    * cos(deltaz)',   {rz,deltaz})");
    w.factory("expr::ImAz   ('rz    * sin(deltaz)',   {rz,deltaz})");
    w.factory("expr::ReApar ('rpar  * cos(deltapar)', {rpar,deltapar})");
    w.factory("expr::ImApar ('rpar  * sin(deltapar)', {rpar,deltapar})");
    w.factory("expr::ReAperp('rperp * cos(deltaperp)',{rperp,deltaperp})");
    w.factory("expr::ImAperp('rperp * sin(deltaperp)',{rperp,deltaperp})");
    
    // TODO: write things in terms of x & y instead of dG and dm...
    w.factory("{tau[1.5,0.5,2.5],dm[0.5,0.3,0.8]]}");
    w.factory("RooGaussModel::res(t,mu[0],sigma[0.05])");
    //
    RooArgSet obs = w.argSet("qtag,qrec,cpsi,ctheta,phi,t");

    // RooCmdArg tagAsym = RooFit::Asymmetry( get<RooCategory>(w,"qtag") );
    // RooCmdArg mixAsym = RooFit::Asymmetry( get<RooCategory>(w,"qtag") );

    RooAbsPdf& pdf = _jpsikstar(w,"pdf");



    //RooMsgService::instance().addStream(RooFit::DEBUG,Topic(RooFit::Generation));
    //RooMsgService::instance().addStream(RooFit::DEBUG,Topic(RooFit::Plotting));
    RooMsgService::instance().addStream(RooFit::DEBUG,Topic(RooFit::NumIntegration));
    //RooMsgService::instance().addStream(RooFit::INFO,Topic(RooFit::Integration));

    TCanvas *c = new TCanvas;
    c->Divide(5,4);
    const char *rname[] = { "rz","rpar","rperp",0 };
    const char *dname[] = { "deltaz","deltapar","deltaperp",0 };
    for (int i=0;i<4;++i) {
        for (int k=0;k<3;++k) { 
            if (i>0) { 
                w.var(rname[k])->setVal( i==k+1 ? 1 : 0 ); 
                w.var(dname[k])->setVal( 0 );
            }
            w.var(rname[k])->Print();
            w.var(dname[k])->Print();
        }
        // if J/psi K*, we should _also_ generate qrec...
        RooAbsData *data = pdf.generate(obs,100000);
        // if (i==0) pdf.fitTo(*data,RooFit::NumCPU(7));
        RooPlot *p1 = w.var("cpsi")->frame();   data->plotOn(p1); pdf.plotOn(p1); c->cd(i*5+1); p1->Draw();
        RooPlot *p2 = w.var("ctheta")->frame(); data->plotOn(p2); pdf.plotOn(p2); c->cd(i*5+2); p2->Draw();
        RooPlot *p3 = w.var("phi")->frame();    data->plotOn(p3); pdf.plotOn(p3); c->cd(i*5+3); p3->Draw();
        RooPlot *p4 = w.var("t")->frame();      data->plotOn(p4); pdf.plotOn(p4); c->cd(i*5+4); p4->Draw();
        // RooPlot *p5 = w.var("t")->frame();      data->plotOn(p5, tagAsym ); pdf.plotOn(p5, tagAsym); c->cd(i*5+5); p5->Draw();
    }
}

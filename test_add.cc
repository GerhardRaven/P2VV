#ifndef _CINT_
#include "utils.h"
#include "basis.h"
#include "RooAbsArg.h"
#include "RooAddition_.h"
#include "RooBDecay.h"
#include "RooWorkspace.h"
#include "RooProduct.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooRealSumPdf.h"
#include "RooFormulaVar.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooMsgService.h"
#include "TCanvas.h"
#endif
using namespace std;

RooAbsReal& prod(RooWorkspace& w,RooAbsReal* a, RooAbsReal* b, RooAbsReal* c, RooAbsReal* d, RooAbsReal* e) {
    assert (a!=0) ;
    if (b==0&&c==0&&d==0) return *a;
    const char *name(0);
    name = a->GetName(); RooArgList l; l.add(*a);
    if (b) { name = Format("%s_x_%s",name,b->GetName()); l.add(*b); }
    if (c) { name = Format("%s_x_%s",name,c->GetName()); l.add(*c); }
    if (d) { name = Format("%s_x_%s",name,d->GetName()); l.add(*d); }
    if (e) { name = Format("%s_x_%s",name,e->GetName()); l.add(*e); }
    return import(w,RooProduct(name,name,l));
}

RooAbsReal& prod(RooWorkspace& w,const char* a, const char *b=0, const char *c=0, const char *d=0, const char *e=0) {
    return prod( w,      &get<RooAbsReal>(w,a)
                  ,  b ? &get<RooAbsReal>(w,b) : 0
                  ,  c ? &get<RooAbsReal>(w,c) : 0
                  ,  d ? &get<RooAbsReal>(w,d) : 0
                  ,  e ? &get<RooAbsReal>(w,e) : 0 );
}

RooAbsPdf& jpsiphi(RooWorkspace& w, const char* name 
        // , RooAbsReal& cpsi, RooAbsReal& ctheta, RooAbsReal& phi, RooAbsReal& t, RooAbsReal& qtag 
        // , RooAbsReal& ReAz, RooAbsReal& ImAz, RooAbsReal ReApar, RooAbsReal& ImApar, RooAbsReal& ReAperp, RooAbsReal& ImAperp
        // , RooAbsReal& tau, RooAbsReal& dGamma, RooAbsReal& dm
        // , RooAbsReal& C, RooAbsReal& S, RooAbsReal& D,
        // , RooAbsReal& w, RooResolutionModel& res
        ) // for now, we stick with a naming convention. In future, pass the actual RooAbsReal, RooAbsArg, ...
{ 
        // definition of the angular part of the PDF in terms of basis functions... 
        abasis ab(w, "cpsi", "ctheta", "phi");  // bind workspace and observables
        import(w, RooAddition_("AzAz_basis",      "AzAz_basis",       RooArgList( ab("AzAz",       0,0,0, 0, 2.), ab("AzAz",       0,0,2,0, sqrt(1./ 5.)),  ab("AzAz",     0,0,2,2, -sqrt( 3./ 5.)) , 
                                                                                  ab("AzAz",       2,0,0, 0, 4.), ab("AzAz",       2,0,2,0, sqrt(4./ 5.)),  ab("AzAz",     2,0,2,2, -sqrt(12./ 5.)) )));
        import(w, RooAddition_("AparApar_basis",  "AparApar_basis",   RooArgList( ab("AparApar",   2,2,0, 0, 1.), ab("AparApar",   2,2,2,0, sqrt(1./20.)),  ab("AparApar", 2,2,2,2,  sqrt( 3./20.)) ))); 
        import(w, RooAddition_("AperpAperp_basis","AperpAperp_basis", RooArgList( ab("AperpAperp", 2,2,0, 0, 1.), ab("AperpAperp", 2,2,2,0,-sqrt(1./ 5.)))));
        import(w, RooAddition_("AparAperp_basis", "AparAperp_basis",  RooArgList( ab("AparAperp",  2,2,2,-1, sqrt( 9./15.)) )));
        import(w, RooAddition_("AzAperp_basis",   "AzAperp_basis",    RooArgList( ab("AzAperp",    2,1,2, 1,-sqrt(18./15.)) )));
        import(w, RooAddition_("AzApar_basis",    "AzApar_basis",     RooArgList( ab("AzApar",     2,1,2,-2, sqrt(18./15.)) )));
        //                        0    1    2       3       4      5      6    7
        RooArgSet amp = w.argSet("ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar,qtag,C");
        import(w,RooFormulaVar("NAzAz",       "( @0 * @0 + @1 * @1 ) / ( 1+@6*@7 )", amp));
        import(w,RooFormulaVar("NAparApar",   "( @4 * @4 + @5 * @5 ) / ( 1+@6*@7 )", amp));
        import(w,RooFormulaVar("NAperpAperp", "( @2 * @2 + @3 * @3 ) / ( 1+@6*@7 )", amp));
        import(w,RooFormulaVar("ReAparAperp", "( @4 * @2 + @5 * @3 ) / ( 1+@6*@7 )", amp));
        import(w,RooFormulaVar("ReAzAperp",   "( @0 * @2 + @1 * @3 ) / ( 1+@6*@7 )", amp));
        import(w,RooFormulaVar("ReAzApar",    "( @0 * @4 + @1 * @5 ) / ( 1+@6*@7 )", amp));
        import(w,RooFormulaVar("ImAparAperp", "( @4 * @3 - @5 * @2 ) / ( 1+@6*@7 )", amp));
        import(w,RooFormulaVar("ImAzAperp",   "( @0 * @3 - @1 * @2 ) / ( 1+@6*@7 )", amp));
        import(w,RooFormulaVar("qtag_",   "@7", amp)); // make a RooAbsReal corresponding to the RooCategory 'qtag'...

        w.factory("Minus[-1]");
        import(w, RooAddition_("f_cosh","f_cosh",RooArgSet( prod(w, "NAzAz",                           "AzAz_basis"      )
                                                          , prod(w, "NAparApar",                       "AparApar_basis"  )
                                                          , prod(w, "NAperpAperp",                     "AperpAperp_basis")
                                                          , prod(w, "ImAparAperp",         "qtag_","C","AparAperp_basis")
                                                          , prod(w, "ImAzAperp",           "qtag_","C","AzAperp_basis")
                                                          , prod(w, "ReAzApar",                        "AzApar_basis") )));
        import(w, RooAddition_("f_cos" ,"f_cos", RooArgSet( prod(w, "NAzAz",               "qtag_","C","AzAz_basis")
                                                          , prod(w, "NAparApar",           "qtag_","C","AparApar_basis")
                                                          , prod(w, "NAperpAperp",         "qtag_","C","AperpAperp_basis")
                                                          , prod(w, "ImAparAperp",                     "AparAperp_basis")
                                                          , prod(w, "ImAzAperp",                       "AzAperp_basis")
                                                          , prod(w, "ReAzApar",            "qtag_","C","AzApar_basis") )));
        import(w, RooAddition_("f_sinh","f_sinh",RooArgSet( prod(w, "NAzAz",       "Minus",        "D","AzAz_basis")
                                                          , prod(w, "NAparApar",   "Minus",        "D","AparApar_basis")
                                                          , prod(w, "NAperpAperp",                 "D","AperpAperp_basis")
                                                          , prod(w, "ReAparAperp",         "qtag_","S","AparAperp_basis")
                                                          , prod(w, "ReAzAperp",           "qtag_","S","AzAperp_basis")
                                                          , prod(w, "ReAzApar",    "Minus",        "D","AzApar_basis") )));
        import(w, RooAddition_("f_sin" ,"f_sin", RooArgSet( prod(w, "NAzAz",       "Minus","qtag_","S","AzAz_basis")
                                                          , prod(w, "NAparApar",   "Minus","qtag_","S","AparApar_basis")
                                                          , prod(w, "NAperpAperp",         "qtag_","S","AperpAperp_basis")
                                                          , prod(w, "ReAparAperp", "Minus",        "D","AparAperp_basis")
                                                          , prod(w, "ReAzAperp",   "Minus",        "D","AzAperp_basis")
                                                          , prod(w, "ReAzApar",    "Minus","qtag_","S","AzApar_basis") )));
        RooResolutionModel *res = dynamic_cast<RooResolutionModel*>( w.function("res") );
        w.import( RooBDecay(name,name,   get<RooRealVar>(w,"t") , get<RooAbsReal>(w,"tau"), get<RooAbsReal>(w,"dG")
                                       , get<RooAbsReal>(w,"f_cosh"), get<RooAbsReal>(w,"f_sinh") , get<RooAbsReal>(w,"f_cos"), get<RooAbsReal>(w,"f_sin")
                                       , get<RooAbsReal>(w,"dm"), *res, RooBDecay::SingleSided) );
        return *w.pdf(name);
};

void test_add() {
    RooWorkspace w("w"); 
    // observables...
    w.factory(Format( "{ cpsi[-1,1], ctheta[-1,1], phi[0,%f], t[-1,4], qtag[bbar=+1,b=-1]} ",4*acos(0.))); // bbar=+1, so code corresponds to Bs(t=0)

    // choice: either fit for the Re&Im of the 3 amplitudes (and then
    //         constrain one phase and the sum of magnitudes)
    //         or fit in terms of angles and relative magnitudes
    // Note: initial values from arXiv:0704.0522v2 [hep-ex] BaBar PUB-07-009
    w.factory("{rz[0.556],rpar[0.211,0.,2.],rperp[0.233,0.,2.]}");
    w.factory("{deltaz[0],deltapar[-2.93,-3.15,3.15],deltaperp[2.91,-3.15,3.15]}");
    w.factory("FormulaVar::ReAz   ('rz    * cos(deltaz)',   {rz,deltaz})");
    w.factory("FormulaVar::ImAz   ('rz    * sin(deltaz)',   {rz,deltaz})");
    w.factory("FormulaVar::ReApar ('rpar  * cos(deltapar)', {rpar,deltapar})");
    w.factory("FormulaVar::ImApar ('rpar  * sin(deltapar)', {rpar,deltapar})");
    w.factory("FormulaVar::ReAperp('rperp * cos(deltaperp)',{rperp,deltaperp})");
    w.factory("FormulaVar::ImAperp('rperp * sin(deltaperp)',{rperp,deltaperp})");

    // choice: either fit for the three degrees of freedom independently
    //         i.e. make S,D,C independent parameters
    //         OR write S,D,C in terms of phi_s
    w.factory("{phis[0.8]}"); // taken as sin(0.8) \approx 0.72, close to the value of sin(2beta) ;-)
    w.factory("FormulaVar::S('sin(phis)',{phis})");
    w.factory("FormulaVar::D('cos(phis)',{phis})");
    w.factory("C[0]");
    // w.factory("{C[0],S[0],D[1]}");

    // TODO: write things in terms of x & y instead of dG and dm...
    w.factory("{tau[1.5,0.5,2.5],dm[5]}");
    w.factory("RooFormulaVar::dG('@0/@1',{dGG[0],tau})"); 
    //w.factory("RooGaussModel::res(t,mu[0],sigma[0.05])");
    w.factory("RooTruthModel::res(t)");

    RooAbsPdf& pdf = jpsiphi(w,"pdf");
    //cout << "PDF:" << endl;
    //pdf.printTree(cout);

    if (false) {
        RooAbsData *data = pdf.generate(w.argSet("qtag,cpsi,ctheta,phi,t"),10000);
        w.import(*data);
        w.writeToFile("p2vv_3.root");
        return;
    }

    if (false) { // marginalize to transversity angle pdf..
        RooAbsPdf &pdf_trans   = import(w,*pdf.createProjection(w.argSet("cpsi,phi")));
    }

    if (false) { // marginalize to untagged, time integrated 3-angle pdf..
        RooAbsPdf &pdf_untagged = import(w,*pdf.createProjection(w.argSet("qtag,t")));
        RooAbsData *data = pdf_untagged.generate(w.argSet("cpsi,ctheta,phi"),1000000);
        data->Print("V");
        w.import(*data);
        w.writeToFile("p2vv_angles.root");
        return;
    }


    RooCmdArg tagAsym = RooFit::Asymmetry( get<RooCategory>(w,"qtag") );
    TCanvas *c = new TCanvas();



    if (false) {
        c->Divide(1,3);
        RooCmdArg red = RooFit::LineColor(kRed);
        RooCmdArg blue = RooFit::LineColor(kBlue);
        RooArgSet projset =  w.argSet("cpsi,ctheta,phi,qtag"); 
        projset.add(get<RooCategory>(w,"qtag"));
        projset.Print("V");
        RooCmdArg proj = RooFit::Project(projset);

        c->cd(1); RooPlot *p6 = w.var("t")->frame();  
        pdf.plotOn(p6, proj); p6->Draw();
        // use RooAbsPdf::createProjection( RooArgSet( cpsi,ctheta,phi,qtag) ) to compare this to...

        c->cd(2); RooPlot *p5 = w.var("t")->frame();  
        pdf.plotOn(p5, proj,RooFit::Slice(get<RooCategory>(w,"qtag"),"bbar"),blue);  
        pdf.plotOn(p5, proj,RooFit::Slice(get<RooCategory>(w,"qtag"),"b"),red); 
        p5->Draw();
        c->cd(3); RooPlot *p7 = w.var("t")->frame();    pdf.plotOn(p7, tagAsym);  p7->Draw();

        return;
    }

    //RooMsgService::instance().addStream(RooFit::DEBUG,Topic(RooFit::Generation));
    //RooMsgService::instance().addStream(RooFit::DEBUG,Topic(RooFit::Plotting));
    RooMsgService::instance().addStream(RooFit::INFO,Topic(RooFit::NumIntegration));
    //RooMsgService::instance().addStream(RooFit::INFO,Topic(RooFit::Integration));

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
        RooAbsData *data = pdf.generate(w.argSet("qtag,cpsi,ctheta,phi,t"),10000);
        if (i==0) pdf.fitTo(*data,RooFit::NumCPU(7));
        RooPlot *p1 = w.var("cpsi")->frame();   data->plotOn(p1); pdf.plotOn(p1); c->cd(i*5+1); p1->Draw();
        RooPlot *p2 = w.var("ctheta")->frame(); data->plotOn(p2); pdf.plotOn(p2); c->cd(i*5+2); p2->Draw();
        RooPlot *p3 = w.var("phi")->frame();    data->plotOn(p3); pdf.plotOn(p3); c->cd(i*5+3); p3->Draw();
        RooPlot *p4 = w.var("t")->frame();      data->plotOn(p4); pdf.plotOn(p4); c->cd(i*5+4); p4->Draw();
        RooPlot *p5 = w.var("t")->frame();      data->plotOn(p5, tagAsym ); pdf.plotOn(p5, tagAsym); c->cd(i*5+5); p5->Draw();
    }
}

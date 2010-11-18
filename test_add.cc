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
        
        import(w, RooFormulaVar("qtag_","@0",RooArgSet( get<RooCategory>(w,"qtag") ) ) );
        //                                                                   0    1    2       3       4      5      6    7
        w.factory("expr::NAzAz      ('( @0 * @0 + @1 * @1 ) / ( 1+@6*@7 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar,qtag_,C})");
        w.factory("expr::NAparApar  ('( @4 * @4 + @5 * @5 ) / ( 1+@6*@7 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar,qtag_,C})");
        w.factory("expr::NAperpAperp('( @2 * @2 + @3 * @3 ) / ( 1+@6*@7 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar,qtag_,C})");
        w.factory("expr::ReAparAperp('( @4 * @2 + @5 * @3 ) / ( 1+@6*@7 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar,qtag_,C})");
        w.factory("expr::ReAzAperp  ('( @0 * @2 + @1 * @3 ) / ( 1+@6*@7 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar,qtag_,C})");
        w.factory("expr::ReAzApar   ('( @0 * @4 + @1 * @5 ) / ( 1+@6*@7 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar,qtag_,C})");
        w.factory("expr::ImAparAperp('( @4 * @3 - @5 * @2 ) / ( 1+@6*@7 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar,qtag_,C})");
        w.factory("expr::ImAzAperp  ('( @0 * @3 - @1 * @2 ) / ( 1+@6*@7 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar,qtag_,C})");

        w.factory("Minus[-1]");
        import(w, RooAddition_("f_cosh","f_cosh",RooArgSet( *w.factory("prod(NAzAz,                    AzAz_basis)")
                                                          , *w.factory("prod(NAparApar,                AparApar_basis)")
                                                          , *w.factory("prod(NAperpAperp,              AperpAperp_basis)")
                                                          , *w.factory("prod(ImAparAperp,      qtag_,C,AparAperp_basis)")
                                                          , *w.factory("prod(ImAzAperp,        qtag_,C,AzAperp_basis)")
                                                          , *w.factory("prod(ReAzApar,                 AzApar_basis)") )));
        import(w, RooAddition_("f_cos" ,"f_cos", RooArgSet( *w.factory("prod(NAzAz,            qtag_,C,AzAz_basis)")
                                                          , *w.factory("prod(NAparApar,        qtag_,C,AparApar_basis)")
                                                          , *w.factory("prod(NAperpAperp,      qtag_,C,AperpAperp_basis)")
                                                          , *w.factory("prod(ImAparAperp,              AparAperp_basis)")
                                                          , *w.factory("prod(ImAzAperp,                AzAperp_basis)")
                                                          , *w.factory("prod(ReAzApar,         qtag_,C,AzApar_basis)") )));
        import(w, RooAddition_("f_sinh","f_sinh",RooArgSet( *w.factory("prod(NAzAz,      Minus,      D,AzAz_basis)")
                                                          , *w.factory("prod(NAparApar,  Minus,      D,AparApar_basis)")
                                                          , *w.factory("prod(NAperpAperp,            D,AperpAperp_basis)")
                                                          , *w.factory("prod(ReAparAperp,      qtag_,S,AparAperp_basis)")
                                                          , *w.factory("prod(ReAzAperp,        qtag_,S,AzAperp_basis)")
                                                          , *w.factory("prod(ReAzApar,   Minus,      D,AzApar_basis)") )));
        import(w, RooAddition_("f_sin" ,"f_sin", RooArgSet( *w.factory("prod(NAzAz,      Minus,qtag_,S,AzAz_basis)")
                                                          , *w.factory("prod(NAparApar,  Minus,qtag_,S,AparApar_basis)")
                                                          , *w.factory("prod(NAperpAperp,      qtag_,S,AperpAperp_basis)")
                                                          , *w.factory("prod(ReAparAperp,Minus,      D,AparAperp_basis)")
                                                          , *w.factory("prod(ReAzAperp,  Minus,      D,AzAperp_basis)")
                                                          , *w.factory("prod(ReAzApar,   Minus,qtag_,S,AzApar_basis)") )));
        w.factory(Format("BDecay::%s(t,tau,dG,f_cosh,f_sinh,f_cos,f_sin,dm,res,SingleSided)",name));
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
    w.factory("RooGaussModel::res(t,mu[0],sigma[0.05])");

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
        RooAbsData *data = pdf.generate(w.argSet("qtag,cpsi,ctheta,phi,t"),100000);
        if (i==0) pdf.fitTo(*data,RooFit::NumCPU(7));
        RooPlot *p1 = w.var("cpsi")->frame();   data->plotOn(p1); pdf.plotOn(p1); c->cd(i*5+1); p1->Draw();
        RooPlot *p2 = w.var("ctheta")->frame(); data->plotOn(p2); pdf.plotOn(p2); c->cd(i*5+2); p2->Draw();
        RooPlot *p3 = w.var("phi")->frame();    data->plotOn(p3); pdf.plotOn(p3); c->cd(i*5+3); p3->Draw();
        RooPlot *p4 = w.var("t")->frame();      data->plotOn(p4); pdf.plotOn(p4); c->cd(i*5+4); p4->Draw();
        RooPlot *p5 = w.var("t")->frame();      data->plotOn(p5, tagAsym ); pdf.plotOn(p5, tagAsym); c->cd(i*5+5); p5->Draw();
    }
}

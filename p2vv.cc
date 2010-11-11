#ifndef _CINT_
#include <iostream>
#include "RooAbsArg.h"
#include "RooWorkspace.h"
#include "RooProduct.h"
#include "RooLegendre.h"
#include "RooSpHarmonic.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooRealSumPdf.h"
#include "RooFormulaVar.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooResolutionModel.h"
#include "RooMsgService.h"
#include "TCanvas.h"
#endif
using namespace std;
/*
gSystem->Load("libMathMore.so")
.L RooLegendre.cxx+g
.L RooSpHarmonic.cxx+g
*/

// pitty: 'asprintf' doesn't work in CINT...
template <typename X>
char *Format(const char* fmt,const X&x) {
   char *buf = new char[1024] ; // let's hope that this is enough....
   sprintf(buf,fmt, x );
   return buf;
}   
template <typename X, typename Y>
char *Format(const char* fmt,const X&x, const Y&y) {
   char *buf = new char[1024] ; // let's hope that this is enough....
   sprintf(buf,fmt, x, y );
   return buf;
}   

template <typename X, typename Y, typename Z>
char *Format(const char* fmt,const X&x, const Y&y, const Z&z) {
   char *buf = new char[1024] ; // let's hope that this is enough....
   sprintf(buf,fmt, x, y, z );
   return buf;
}   

template <typename T>
T& get(RooWorkspace& w, T*(RooWorkspace::*fun)(const char*)const, const char* name) {
    T* x = (w.*fun)(name);
    if (x==0) {
        cout << "FAILURE: could not retrieve '" << name << "' from workspace" << endl;
        assert(1==0);
    }
    return *x;
}

//TODO:  use traits to infer member function from type
// typedef  T* (RooWorkspace::*fun)(const char* name);
//    mytraits<RooAbsPdf>::fun  = RooWorkspace::pdf
//    mytraits<RooAbsReal>::fun = RooWorkspace::function

template <typename T> T&   get(RooWorkspace& w, const char *name) { return dynamic_cast<T&>(*w.obj(name)); }
template<> RooAbsReal&     get<RooAbsReal>(RooWorkspace& w, const char* name) { return get<RooAbsReal>(w,&RooWorkspace::function,name); }
template<> RooAbsArg&      get<RooAbsArg>(RooWorkspace& w, const char* name)   { return get<RooAbsArg>(w,&RooWorkspace::arg,name); } // arg(w,name); }
template<> RooRealVar&     get<RooRealVar>(RooWorkspace& w, const char* name)   { return get<RooRealVar>(w,&RooWorkspace::var, name); } // (w,name); }
template<> RooFormulaVar&  get<RooFormulaVar>(RooWorkspace& w, const char* name)   { return dynamic_cast<RooFormulaVar&>( get<RooAbsReal>(w,name)); }

template <typename T>
T& import(RooWorkspace& w, const T& r, const char *n=0) {
   if (n==0) { w.import(r); return get<T>(w,r.GetName()); }
   else      { w.import(r,RooFit::RenameVariable(r.GetName(),n)); return get<T>(w,n); }
}

RooAbsArg& product(RooWorkspace& w, const RooArgList& x) {
    const char *name(0);
    for (int i=0;i<x.getSize();++i) {
        name = name ? Format("%s_%s",name,x[i].GetName()) : x[i].GetName();
    }
    RooAbsArg *p = w.function(name);
    if (p==0) p = &import(w,RooProduct(name,name,x));
    return *p;
}

RooAbsArg& product(RooWorkspace& w, RooAbsArg& x, RooAbsArg& y) {
    return product( w, RooArgList(x,y) );
}

RooAbsArg& product(RooWorkspace& w, RooAbsArg& x, RooAbsArg& y, RooAbsArg& z) {
    return product( w, RooArgList(x,y,z) );
}

RooArgList product(RooWorkspace& w, const RooArgList& x, const RooArgList& y) {
     RooArgList z;
     for (int i=0;i<x.getSize();++i) for (int j=0;j<y.getSize();++j) {
         z.add(product(w,x[i],y[j]));
     }
     return z;
}

class abasis {  //TODO: make this an RooAbsReal implementation, which forwards integrals, 
                //      supports maxVal for toy generation, 
                //      and avoids making a constant to put in front of the product
                //  OR  make this class return a dedicated RooAbsReal implementation which
                //      does the above, i.e. a RooAbsReal which represents C * P_ij * Y_lm
public:
    // TODO: add version which takes RooAbsArg for cpsi, ctheta, phi (and then checks they 
    //  are in the workspace w!)
    abasis(RooWorkspace &w, const char *cpsi, const char *ctheta, const char *phi) 
        : _w(w)
        , _cpsi( get<RooRealVar>(w,cpsi) )
        , _ctheta( get<RooRealVar>(w,ctheta) )
        , _phi( get<RooRealVar>(w,phi) )
    { }

    RooAbsArg& operator()(const char* label, int i, int j, int k, int l, double c)  {
         char *name = Format("%s_%d_%d",label,i,j);
         name = Format( l<0 ? "%s_%d_m%d" : "%s_%d_%d",name,k,l<0?-l:l);
         _w.factory(Format("%s[%f]",name,c));
         return product( _w, *_w.var(name), Plm(i,j), Ylm(k,l));
    }
private:
    RooAbsReal& Plm(int i, int j) {
        char *name = Format( j<0 ? "P_%d_m%d" : "P_%d_%d", i, j<0?-j:j);
        RooAbsReal *P = _w.function(name);
        if (P==0) P = &import(_w, RooLegendre(name,name,_cpsi,i,j));
        return *P;
    }
    RooAbsReal& Ylm(int i, int j) {
        char *name = Format( j<0 ? "Y_%d_m%d" : "Y_%d_%d", i, j<0?-j:j);
        RooAbsReal *Y = _w.function(name);
        if (Y==0) Y = &import( _w, RooSpHarmonic(name,name,_ctheta,_phi,i,j));
        return *Y;
    }
    RooWorkspace &_w;
    RooAbsReal &_cpsi;
    RooAbsReal &_ctheta;
    RooAbsReal &_phi;
};

class tbasis {
public:
    tbasis(RooWorkspace& w, RooResolutionModel& r, RooFormulaVar& basis, const char *label) 
        : _w(w)
        , _basis(import(w, *r.convolution(&basis,&r.convVar()), label))
    {  
    }
    
    RooAbsArg& operator()(const char* t1, const char *t2=0, const char *t3=0, const char *t4=0) {
       RooArgList l;
       const char *name(0);
       name = add( name, t1, l);
       name = add( name, t2, l);
       name = add( name, t3, l);
       name = add( name, t4, l);
       l.add( _basis);
       return product(_w,l);
    }
private:
    const char *add( const char *name, const char *t, RooArgList& l) {
        if(t==0) return name;
        l.add( get<RooAbsArg>(_w, t) );
        return name ? Format("%s_%s",name,t) : t;
    }
    RooWorkspace& _w;
    RooAbsReal&   _basis;
};

RooAbsPdf& jpsiphi(RooWorkspace& w, const char* name /*, RooAbsReal& cpsi, RooAbsReal& ctheta, RooAbsReal& phi, RooAbsReal& t, RooAbsReal& qtag */) // for now, we stick with a naming convention. In future, pass the actual RooAbsReal, RooAbsArg, ...
{ 
        RooArgSet amp = w.argSet("ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar,qtag,C");
        import(w,RooFormulaVar("NAzAz",       "( @0 * @0 + @1 * @1 ) / ( 1+@6*@7 )", amp));
        import(w,RooFormulaVar("NAparApar",   "( @4 * @4 + @5 * @5 ) / ( 1+@6*@7 )", amp));
        import(w,RooFormulaVar("NAperpAperp", "( @2 * @2 + @3 * @3 ) / ( 1+@6*@7 )", amp));
        import(w,RooFormulaVar("ReAparAperp", "( @4 * @2 + @5 * @3 ) / ( 1+@6*@7 )", amp));
        import(w,RooFormulaVar("ReAzAperp",   "( @0 * @2 + @1 * @3 ) / ( 1+@6*@7 )", amp));
        import(w,RooFormulaVar("ReAzApar",    "( @0 * @4 + @1 * @5 ) / ( 1+@6*@7 )", amp));
        import(w,RooFormulaVar("ImAparAperp", "( @4 * @3 - @5 * @2 ) / ( 1+@6*@7 )", amp));
        import(w,RooFormulaVar("ImAzAperp",   "( @0 * @3 - @1 * @2 ) / ( 1+@6*@7 )", amp));

        RooResolutionModel *res = dynamic_cast<RooResolutionModel*>( w.function("res") );
        RooRealVar& t = res->convVar();

        RooArgList ttdg( t, get<RooAbsReal>(w,"tau"), get<RooAbsReal>(w,"dG"));
        RooArgList ttdm( t, get<RooAbsReal>(w,"tau"), get<RooAbsReal>(w,"dm"));
        tbasis _sinh_(w, *res, import(w,RooFormulaVar("sinhBasis","exp(-@0/@1)*sinh(@0*@2/2)",ttdg)), "sinh_t");
        tbasis _cosh_(w, *res, import(w,RooFormulaVar("coshBasis","exp(-@0/@1)*cosh(@0*@2/2)",ttdg)), "cosh_t");
        tbasis _sin_ (w, *res, import(w,RooFormulaVar("sinBasis", "exp(-@0/@1)*sin(@0*@2)",   ttdm)), "sin_t" );
        tbasis _cos_ (w, *res, import(w,RooFormulaVar("cosBasis", "exp(-@0/@1)*cos(@0*@2)",   ttdm)), "cos_t" );

        RooArgList *time[6];
        w.factory("Minus[-1]");
        time[0] = new RooArgList( _cosh_("NAzAz"                 ), _cos_( "NAzAz",      "qtag","C"), _sinh_("NAzAz",      "Minus",       "D"), _sin_("NAzAz",       "Minus","qtag","S") );
        time[1] = new RooArgList( _cosh_("NAparApar"             ), _cos_( "NAparApar",  "qtag","C"), _sinh_("NAparApar",  "Minus",       "D"), _sin_("NAparApar",   "Minus","qtag","S") );
        time[2] = new RooArgList( _cosh_("NAperpAperp"           ), _cos_( "NAperpAperp","qtag","C"), _sinh_("NAperpAperp",               "D"), _sin_("NAperpAperp",         "qtag","S") );
        time[3] = new RooArgList( _cosh_("ImAparAperp","qtag","C"), _cos_( "ImAparAperp"           ), _sinh_("ReAparAperp",        "qtag","S"), _sin_("ReAparAperp", "Minus",       "D") );
        time[4] = new RooArgList( _cosh_("ImAzAperp",  "qtag","C"), _cos_( "ImAzAperp"             ), _sinh_("ReAzAperp",          "qtag","S"), _sin_("ReAzAperp",   "Minus",       "D") );
        time[5] = new RooArgList( _cosh_("ReAzApar"              ), _cos_( "ReAzApar",   "qtag","C"), _sinh_("ReAzApar",   "Minus",       "D"), _sin_("ReAzApar",    "Minus","qtag","S") );
        //
        // for now, hardwire the max value...
        // double tmax[6] = { abs(get<RooAbsReal>(w,"NAzAz").getVal())*(1+ abs("C")+abs("D")+abs("S") )
        //                  , get<RooAbsReal>(w,"NAparApar").getVal())*(1+ abs("C")+abs("D")+abs("S") )
        //                  , sqrt(9./15.), sqrt(18./15.), sqrt(18./15.) }

        // definition of the angular part of the PDF in terms of basis functions... 
        // NOTE: the maximum value of each 'angle' function has an upper limit 
        //       given by the sum of the (absolute) values of the coefficients
        RooArgList *angle[6];
        abasis ab(w, "cpsi", "ctheta", "phi");  // bind workspace and observables
        angle[0] = new RooArgList( ab("AzAz",       0,0,0, 0, 2.), ab("AzAz",       0,0,2,0, sqrt(1./ 5.)),  ab("AzAz",     0,0,2,2, -sqrt(3. / 5.)) , 
                                   ab("AzAz",       2,0,0, 0, 4.), ab("AzAz",       2,0,2,0, sqrt(4./ 5.)),  ab("AzAz",     2,0,2,2, -sqrt(12./ 5.)) );
        angle[1] = new RooArgList( ab("AparApar",   2,2,0, 0, 1.), ab("AparApar",   2,2,2,0, sqrt(1./20.)),  ab("AparApar", 2,2,2,2,  sqrt(3. /20.)) ); 
        angle[2] = new RooArgList( ab("AperpAperp", 2,2,0, 0, 1.), ab("AperpAperp", 2,2,2,0,-sqrt(1./ 5.)));
        angle[3] = new RooArgList( ab("AparAperp",  2,2,2,-1, sqrt( 9./15.)) );
        angle[4] = new RooArgList( ab("AzAperp",    2,1,2, 1,-sqrt(18./15.)) );
        angle[5] = new RooArgList( ab("AzApar",     2,1,2,-2, sqrt(18./15.)) );

        // for now, hardwire the max value...
        // double amax[6] = { 3*(2+sqrt(1./5.)+sqrt(3./5.)), 1+sqrt(1./20.)+sqrt(12./5.),1 + sqrt(1./20.)+sqrt(3./20.)
        //                  , sqrt(9./15.), sqrt(18./15.), sqrt(18./15.) }

        RooArgList f;
        double max_(0);
        for (int i=0;i<6;++i) {
            f.add( product(w, *angle[i], *time[i] ) );
            //max_ += tmax[i]*amax[i];
        }

        w.factory("One[1]");
        RooArgList dummy;
        for (int i=0; i<f.getSize();++i) dummy.add(*w.var("One"));

        w.import( RooRealSumPdf(name,name,f,dummy) );
        return *w.pdf(name);
};

void p2vv() {
    RooWorkspace w("w",kTRUE); 
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
    cout << "PDF:" << endl;
    pdf.printTree(cout);

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

        c->cd(2); RooPlot *p5 = w.var("t")->frame();  
        pdf.plotOn(p5, proj,RooFit::Slice(get<RooCategory>(w,"qtag"),"bbar"),blue);  
        pdf.plotOn(p5, proj,RooFit::Slice(get<RooCategory>(w,"qtag"),"b"),red); 
        p5->Draw();
        c->cd(3); RooPlot *p7 = w.var("t")->frame();    pdf.plotOn(p7, tagAsym);  p7->Draw();

        return;
    }

    RooMsgService::instance().addStream(RooFit::DEBUG,Topic(RooFit::Generation));
    RooMsgService::instance().addStream(RooFit::DEBUG,Topic(RooFit::Plotting));
    RooMsgService::instance().addStream(RooFit::DEBUG,Topic(RooFit::NumIntegration));

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
        if (i==0) pdf.fitTo(*data,RooFit::NumCPU(8));
        RooPlot *p1 = w.var("cpsi")->frame();   data->plotOn(p1); pdf.plotOn(p1); c->cd(i*5+1); p1->Draw();
        RooPlot *p2 = w.var("ctheta")->frame(); data->plotOn(p2); pdf.plotOn(p2); c->cd(i*5+2); p2->Draw();
        RooPlot *p3 = w.var("phi")->frame();    data->plotOn(p3); pdf.plotOn(p3); c->cd(i*5+3); p3->Draw();
        RooPlot *p4 = w.var("t")->frame();      data->plotOn(p4); pdf.plotOn(p4); c->cd(i*5+4); p4->Draw();
        RooPlot *p5 = w.var("t")->frame();      data->plotOn(p5, tagAsym ); pdf.plotOn(p5, tagAsym); c->cd(i*5+5); p5->Draw();
        break;
    }
}

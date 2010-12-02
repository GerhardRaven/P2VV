from ROOT import *
gSystem.Load("libp2vv")

def buildAngularBasis(ws, ab) :
    #definition of the angular part of the PDF in terms of basis functions...
    ws.put(RooAddition_("AzAz_basis",            "AzAz_basis",        RooArgSet( ab("AzAz",       0,0,0, 0, 2.), ab("AzAz",       0,0,2,0, sqrt(1./ 5.)), ab("AzAz",      0,0,2,2, -sqrt( 3./ 5.))
                                                                               , ab("AzAz",       2,0,0, 0, 4.), ab("AzAz",       2,0,2,0, sqrt(4./ 5.)), ab("AzAz",      2,0,2,2, -sqrt(12./ 5.)) )))
    ws.put(RooAddition_("AparApar_basis",        "AparApar_basis",    RooArgSet( ab("AparApar",   2,2,0, 0, 1.), ab("AparApar",   2,2,2,0, sqrt(1./20.)), ab("AparApar",  2,2,2,2,  sqrt( 3./20.)) )))
    ws.put(RooAddition_("AperpAperp_basis",      "AperpAperp_basis",  RooArgSet( ab("AperpAperp", 2,2,0, 0, 1.), ab("AperpAperp", 2,2,2,0,-sqrt(1./ 5.)))))
    ws.put(RooAddition_("AparAperp_basis",       "AparAperp_basis",   RooArgSet( ab("AparAperp",  2,2,2,-1, sqrt( 9./15.)))))
    ws.put(RooAddition_("AzAperp_basis",         "AzAperp_basis",     RooArgSet( ab("AzAperp",    2,1,2, 1,-sqrt(18./15.)))))
    ws.put(RooAddition_("AzApar_basis",          "AzApar_basis",      RooArgSet( ab("AzApar",     2,1,2,-2, sqrt(18./15.)))))

def buildJpsiphi(ws, name ) :
    buildAngularBasis(ws, abasis(ws,'trcospsi','trcostheta','trphi') )

    # define the relevant combinations of strong amplitudes
    ws.factory("expr::NAzAz      ('( @0 * @0 + @1 * @1 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::NAparApar  ('( @4 * @4 + @5 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::NAperpAperp('( @2 * @2 + @3 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAparAperp('( @4 * @2 + @5 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAzAperp  ('( @0 * @2 + @1 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAzApar   ('( @0 * @4 + @1 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ImAparAperp('( @4 * @3 - @5 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ImAzAperp  ('( @0 * @3 - @1 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")


    ws.put(RooFormulaVar("qtag_","@0*(1-2*@1)",RooArgList( ws.cat('tagdecision'),ws.var('wmistag')) ) )

    #ws.factory("expr::N('1-@0*@1',{qtag_,C})") #in J/psi K*, we need to drop this factor...
    ws.factory("expr::N('1/(1+@0*@1)',{qtag_,C})") #in J/psi K*, we need to drop this factor...

    ws.factory("Minus[-1]");
    ws.factory("$Alias(Addition_,sum_)") ;

    ws.factory("sum_::fjpsiphi_cosh({ prod(N,NAzAz,                    AzAz_basis)"
                                   ", prod(N,NAparApar,                AparApar_basis)"
                                   ", prod(N,NAperpAperp,              AperpAperp_basis)"
                                   ", prod(N,ImAparAperp,      qtag_,C,AparAperp_basis)"
                                   ", prod(N,ImAzAperp,        qtag_,C,AzAperp_basis)"
                                   ", prod(N,ReAzApar,                 AzApar_basis)"
                                   "})")
    ws.factory("sum_::fjpsiphi_cos ({ prod(N,NAzAz,            qtag_,C,AzAz_basis)"
                                   ", prod(N,NAparApar,        qtag_,C,AparApar_basis)"
                                   ", prod(N,NAperpAperp,      qtag_,C,AperpAperp_basis)"
                                   ", prod(N,ImAparAperp,              AparAperp_basis)"
                                   ", prod(N,ImAzAperp,                AzAperp_basis)"
                                   ", prod(N,ReAzApar,         qtag_,C,AzApar_basis)"
                                   "})") 
    ws.factory("sum_::fjpsiphi_sinh({ prod(N,NAzAz,      Minus,      D,AzAz_basis)"
                                   ", prod(N,NAparApar,  Minus,      D,AparApar_basis)"
                                   ", prod(N,NAperpAperp,            D,AperpAperp_basis)"
                                   ", prod(N,ReAparAperp,      qtag_,S,AparAperp_basis)"
                                   ", prod(N,ReAzAperp,        qtag_,S,AzAperp_basis)"
                                   ", prod(N,ReAzApar,   Minus,      D,AzApar_basis)"
                                   "})")
    ws.factory("sum_::fjpsiphi_sin ({ prod(N,NAzAz,      Minus,qtag_,S,AzAz_basis)"
                                   ", prod(N,NAparApar,  Minus,qtag_,S,AparApar_basis)"
                                   ", prod(N,NAperpAperp,      qtag_,S,AperpAperp_basis)"
                                   ", prod(N,ReAparAperp,Minus,      D,AparAperp_basis)"
                                   ", prod(N,ReAzAperp,  Minus,      D,AzAperp_basis)"
                                   ", prod(N,ReAzApar,   Minus,qtag_,S,AzApar_basis)"
                                   "})")
    ws.factory("BDecay::%s(t,tau,dG,fjpsiphi_cosh,fjpsiphi_sinh,fjpsiphi_cos,fjpsiphi_sin,dm,res,SingleSided)" % name)
    return ws.pdf(name)







from math import pi





ws = RooWorkspace('ws')
RooWorkspace.put = getattr(RooWorkspace,'import') # import is a reserved word, so we call it 'put' from now on...

ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f], t[0.3,12], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged = 0]}"%(-pi,pi))

##choice: either fit for the Re&Im of the 3 amplitudes (and then
##        constrain one phase and the sum of magnitudes)
##        or fit in terms of angles and relative magnitudes
##         Note: initial values from arXiv:0704.0522v2 [hep-ex] BaBar PUB-07-009

ws.factory("{rz[0.556,0.,1.],rpar[0.211],rperp[0.233,0.,1.]}")
#ws.factory("{deltaz[0],deltapar[-2.93,-3.15,3.15],deltaperp[2.91,-3.15,3.15]}")
ws.factory("{deltaz[0],deltapar[-2.93],deltaperp[2.91]}")
ws.factory("expr::ReAz   ('rz    * cos(deltaz)',   {rz,deltaz})")
ws.factory("expr::ImAz   ('rz    * sin(deltaz)',   {rz,deltaz})")
ws.factory("expr::ReApar ('rpar  * cos(deltapar)', {rpar,deltapar})")
ws.factory("expr::ImApar ('rpar  * sin(deltapar)', {rpar,deltapar})")
ws.factory("expr::ReAperp('rperp * cos(deltaperp)',{rperp,deltaperp})")
ws.factory("expr::ImAperp('rperp * sin(deltaperp)',{rperp,deltaperp})")

ws.factory("{gamma[0.68,0.4,0.9],dm[17.7],dG[0.075,-0.3,0.3]}")
ws.factory("RooFormulaVar::tau('1/@0',{gamma})")

ws.factory("RooGaussModel::res(t,mu[0],sigma[0.05])")
ws.factory("{wmistag[0.5]}")

##choice: either fit for the three degrees of freedom independently
##        i.e. make S,D,C independent parameters
##ws.factory("{S[0.717,-1,1],D[0.696,-1,1],C[0,-1,1]}")
##        or write S,D,C in terms of phi_s
ws.factory("{expr::S('sin(phis)',{phis[0.]}),expr::D('cos(phis)',{phis}),C[0]}")
##        The no-CP violation case:
##ws.factory("{S[0],C[0],D[1]}")
##obs = ws.argSet("tagdecision,trcospsi,trcostheta,trphi,t")
##        For J/psi K*, C = +1/-1 depending on the final state flavour
##ws.factory("{S[0],D[0],expr::C('@0',{qrec[jpsikstar=+1,jpsikstarbar=-1]})}")
##obs.add(w.argSet("qrec"))

tagAsym = RooFit.Asymmetry(ws.cat("tagdecision"))

#######################################################################################################################################

p2vv = buildJpsiphi(ws,'jpsiphipdf')

file = TFile("WSJpsiPhiPdf.root","RECREATE")
ws.Write("workspace")
file.Close()

canvas = TCanvas('canvas','canvas')
canvas.Divide(1,4)

j=1
for i in [ 'trcospsi','trcostheta','trphi','t' ] :
   canvas.cd(j)
   j=j+1
   f = ws.var(i).frame() 
   p2vv.plotOn(f)
   f.Draw()


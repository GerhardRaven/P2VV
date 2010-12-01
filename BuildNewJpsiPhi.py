from ROOT import *
gSystem.Load("libp2vv")
from math import pi

def abasis(label,i,j,k,l,c):
    string = "%s_%d_%d"%(label,i,j)
    if l<0:
        name = "%s_%d_m%d"%(string,k,l)
    else:
        name = "%s_%d_%d"%(string,k,l)
    basisfunc = RooP2VVAngleBasis(name,name,cpsi,ctheta,phi,i,j,k,l,c)
    #ret = {basisfunc
    return basisfunc

ws = RooWorkspace('ws')

ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f], t[0.3,12], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged = 0]}"%(-pi,pi))
#ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[0,%f], t[-0.5,8], tagdecision[bbar=+1,b=-1]}"%(pi))

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

bs = ws.argSet("tagdecision,trcospsi,trcostheta,trphi,t")

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

cpsi = ws.var('trcospsi')
ctheta = ws.var('trcostheta')
phi = ws.var('trphi')

#definition of the angular part of the PDF in terms of basis functions...
AzAz1 = abasis("AzAz1",       0,0,0, 0, 2.)
AzAz2 = abasis("AzAz2",       0,0,2, 0,sqrt(1./5.))
AzAz3 = abasis("AzAz3",       0,0,2,2, -sqrt( 3./ 5.))
AzAz4 = abasis("AzAz4",       2,0,0, 0, 4.)
AzAz5 = abasis("AzAz5",       2,0,2,0, sqrt(4./ 5.))
AzAz6 = abasis("AzAz6",       2,0,2,2, -sqrt(12./ 5.))

AzAzSet = RooArgSet(AzAz1,AzAz2,AzAz3,AzAz4,AzAz5,AzAz6)

AparApar1 = abasis("AparApar1",   2,2,0, 0, 1.)
AparApar2 = abasis("AparApar2",   2,2,2,0, sqrt(1./20.))
AparApar3 = abasis("AparApar3",   2,2,2,2,  sqrt( 3./20.))

AparAparSet = RooArgSet(AparApar1,AparApar2,AparApar3)

AperpAperp1 = abasis("AperpAperp1", 2,2,0, 0, 1.)
AperpAperp2 = abasis("AperpAperp2", 2,2,2,0,-sqrt(1./ 5.))

AperpAperpSet = RooArgSet(AperpAperp1,AperpAperp2)

AparAperp = abasis("AparAperp",  2,2,2,-1, sqrt( 9./15.))

AparAperpSet = RooArgSet(AparAperp)

AzAperp = abasis("AzAperp",    2,1,2, 1,-sqrt(18./15.))

AzAperpSet = RooArgSet(AzAperp)

AzApar = abasis("AzApar",     2,1,2,-2, sqrt(18./15.))

AzAparSet = RooArgSet(AzApar)

getattr(ws,'import')(RooAddition_("AzAz_basis",            "AzAz_basis",             AzAzSet))
getattr(ws,'import')(RooAddition_("AparApar_basis",        "AparApar_basis",         AparAparSet))
getattr(ws,'import')(RooAddition_("AperpAperp_basis",      "AperpAperp_basis",       AperpAperpSet))
getattr(ws,'import')(RooAddition_("AparAperp_basis",       "AparAperp_basis",        AparAperpSet))
getattr(ws,'import')(RooAddition_("AzAperp_basis",         "AzAperp_basis",          AzAperpSet))
getattr(ws,'import')(RooAddition_("AzApar_basis",          "AzApar_basis",           AzAparSet))

ws.factory("expr::NAzAz      ('( @0 * @0 + @1 * @1 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
ws.factory("expr::NAparApar  ('( @4 * @4 + @5 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
ws.factory("expr::NAperpAperp('( @2 * @2 + @3 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
ws.factory("expr::ReAparAperp('( @4 * @2 + @5 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
ws.factory("expr::ReAzAperp  ('( @0 * @2 + @1 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
ws.factory("expr::ReAzApar   ('( @0 * @4 + @1 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
ws.factory("expr::ImAparAperp('( @4 * @3 - @5 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
ws.factory("expr::ImAzAperp  ('( @0 * @3 - @1 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")


getattr(ws,'import')(RooFormulaVar("qtag_","@0*(1-2*@1)",RooArgList( ws.cat('tagdecision'),ws.var('wmistag')) ) )

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
ws.factory("BDecay::jpsiphipdf(t,tau,dG,fjpsiphi_cosh,fjpsiphi_sinh,fjpsiphi_cos,fjpsiphi_sin,dm,res,SingleSided)")

#######################################################################################################################################

p2vv = ws.pdf('jpsiphipdf')

#data = p2vv.generate(ws.argSet("tagdecision,trcospsi,trcostheta,trphi,t"),100000)
#getattr(ws,'import')(data)

file = TFile("WSJpsiPhiPdf.root","RECREATE")
ws.Write("workspace")
file.Close()

assert False
p2vv.fitTo(data)

t = ws.var('t')
cpsi = ws.var('trcospsi')
ctheta = ws.var('trcostheta')
phi = ws.var('trphi')

tframe = t.frame()
cpsiframe = cpsi.frame()
cthetaframe = ctheta.frame()
phiframe = phi.frame()

canvas = TCanvas('canvas','canvas')
canvas.Divide(4,1)

canvas.cd(1)
data.plotOn(tframe,RooLinkedList())
p2vv.plotOn(tframe)
tframe.Draw()

canvas.cd(2)
data.plotOn(cpsiframe,RooLinkedList())
p2vv.plotOn(cpsiframe)
cpsiframe.Draw()

canvas.cd(3)
data.plotOn(cthetaframe,RooLinkedList())
p2vv.plotOn(cthetaframe)
cthetaframe.Draw()

canvas.cd(4)
data.plotOn(phiframe,RooLinkedList())
p2vv.plotOn(phiframe)
phiframe.Draw()


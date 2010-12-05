from ROOT import *
from itertools import count
gSystem.Load("libp2vv")

from ModelBuilders import buildJpsiphi,buildJpsikstar

from math import pi

ws = RooWorkspace('ws')
# create observables
ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f], t[0.3,12.], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0]}"%(-pi,pi))

##choice: either fit for the Re&Im of the 3 amplitudes (and then
##        constrain one phase and the sum of magnitudes)
##        or fit in terms of angles and relative magnitudes
##         Note: initial values from arXiv:0704.0522v2 [hep-ex] BaBar PUB-07-009

ws.factory("{rz[0.463,0.1,0.9],rpar[0.211],rperp[0.347,0.1,0.9]}")
ws.factory("{deltaz[0],deltapar[-2.93],deltaperp[2.91]}")
ws.factory("expr::ReAz   ('rz    * cos(deltaz)',   {rz,deltaz})")
ws.factory("expr::ImAz   ('rz    * sin(deltaz)',   {rz,deltaz})")
ws.factory("expr::ReApar ('rpar  * cos(deltapar)', {rpar,deltapar})")
ws.factory("expr::ImApar ('rpar  * sin(deltapar)', {rpar,deltapar})")
ws.factory("expr::ReAperp('rperp * cos(deltaperp)',{rperp,deltaperp})")
ws.factory("expr::ImAperp('rperp * sin(deltaperp)',{rperp,deltaperp})")

ws.factory("{gamma[0.68,0.4,0.9],dm[17.7],dG[0.05,-0.3,0.3]}")
ws.factory("expr::tau('1/@0',{gamma})")

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


#######################################################################################################################################

pdf = buildJpsiphi(ws,'jpsiphipdf')  ## for now we rely quite a bit on a naming convention -- 
                                     ## in future we should pass more information into the builder
                                     ## maybe a dictionary of what's what...


### let's make some nice plots to show what this PDF looks like...

obsNames =[ 'trcospsi','trcostheta','trphi','t','tagdecision' ]
obs = ws.argSet(','.join(obsNames))

## canvas = TCanvas('canvas','canvas')
## canvas.Divide(5,4)

## p = [ 'rz','rpar','rperp' ]
## for i in  range(len(p)+1) :
##     if i>0 : 
##         for j in range(len(p)) : 
##             x = ws.var(p[j])
##             x.setVal( 1 if i==j+1 else 0 )
##     for l in p + ['deltaz','deltapar','deltaperp' ] : ws.var(l).Print()

##     for (j,k) in zip(['trcospsi','trcostheta','trphi','t'],count(5*i+1)) :
##        canvas.cd(k)
##        f = ws.var(j).frame() 
##        proj = RooArgSet( obs )
##        proj.remove( ws.var(j) )
##        pdf.plotOn(f,RooFit.Project(proj))
##        f.Draw()

##     tagAsym = RooFit.Asymmetry(ws.cat("tagdecision"))
##     canvas.cd(5*i+5)
##     f = ws.var('t').frame(RooFit.Range(-.5,3.5))
##     proj = RooArgSet( obs )
##     proj.remove(ws.var('t'))
##     pdf.plotOn(f,RooFit.Project(proj),tagAsym)
##     f.Draw()

## canvas.Flush()

file = TFile("NewWSJpsiPhiPdf.root","RECREATE")
ws.Write("workspace")
file.Close()

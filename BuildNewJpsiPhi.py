from ROOT import *
gSystem.Load("libp2vv")

from ModelBuilders import buildJpsiphi,buildJpsikstar

from math import pi

ws = RooWorkspace('ws')
RooWorkspace.put = getattr(RooWorkspace,'import') # import is a reserved word, so we call it 'put' from now on...

ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f], t[-3,12], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1], qrec[jpsikstar=-1,jpsikstarbar=1]}"%(-pi,pi))

##choice: either fit for the Re&Im of the 3 amplitudes (and then
##        constrain one phase and the sum of magnitudes)
##        or fit in terms of angles and relative magnitudes
##         Note: initial values from arXiv:0704.0522v2 [hep-ex] BaBar PUB-07-009

ws.factory("{rz[0.556],rpar[0.211],rperp[0.233]}")
#ws.factory("{deltaz[0],deltapar[-2.93,-3.15,3.15],deltaperp[2.91,-3.15,3.15]}")
ws.factory("{deltaz[0],deltapar[-2.93],deltaperp[2.91]}")
ws.factory("expr::ReAz   ('rz    * cos(deltaz)',   {rz,deltaz})")
ws.factory("expr::ImAz   ('rz    * sin(deltaz)',   {rz,deltaz})")
ws.factory("expr::ReApar ('rpar  * cos(deltapar)', {rpar,deltapar})")
ws.factory("expr::ImApar ('rpar  * sin(deltapar)', {rpar,deltapar})")
ws.factory("expr::ReAperp('rperp * cos(deltaperp)',{rperp,deltaperp})")
ws.factory("expr::ImAperp('rperp * sin(deltaperp)',{rperp,deltaperp})")

ws.factory("{gamma[0.68,0.4,0.9],dm[17.7],dG[0.075,-0.3,0.3]}")
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

tagAsym = RooFit.Asymmetry(ws.cat("tagdecision"))

#######################################################################################################################################

p2vv = buildJpsikstar(ws,'jpsikstarpdf')  ## for now we rely quite a bit on a naming convention -- 
                                      ## in future we should pass more information into the builder

obsNames =[ 'trcospsi','trcostheta','trphi','t' ]
obs = ws.argSet(','.join(obsNames))

canvas = TCanvas('canvas','canvas')
canvas.Divide(5,4)


p = [ 'rz','rpar','rperp' ]
for i in  range(1) : # range(len(p)+1) :
    if i>0 : 
        for j in range(len(p)) : ws.var(p[j]).setVal( 1 if i==j+1 else 0 )

    print '*'*100
    for k in p : ws.var(k).Print()

    #data = p2vv.generate(obs,RooFit.NumEvents(100000))
    #data.Print("V")

    for (j,k) in zip(obsNames,range(5*i+1,100)) :
       pad = canvas.cd(k)
       print pad
       f = ws.var(j).frame() 
       data.plotOn(f)
       p2vv.plotOn(f)
       f.Draw()


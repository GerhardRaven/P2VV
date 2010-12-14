from ROOT import *

import math
from math import pi

from itertools import count
gSystem.Load("libp2vv")

from ModelBuilders import buildJpsiphi,buildJpsikstar

########################
### Define workspace ###
########################

ws = RooWorkspace("ws")

###################
### Observables ###
###################

ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f], t[-2,20.], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0]}"%(-pi,pi))


ws.factory("{rz2[0.460,0.4,0.7],rperp2[0.347,0.2,0.5]")
ws.factory("RooFormulaVar::rpar2('1-@0-@1',{rz2,rperp2})")
ws.factory("RooFormulaVar::rz('sqrt(@0)',{rz2})")
ws.factory("RooFormulaVar::rperp('sqrt(@0)',{rperp2})")
ws.factory("RooFormulaVar::rpar('sqrt(@0)',{rpar2})")

ws.factory("{deltaz[0.],deltapar[-2.93],deltaperp[2.91]}")

t      = ws.var('t')
ctheta = ws.var('trcostheta')
cpsi   = ws.var('trcospsi')
phi    = ws.var('trphi')


##########################
### physics parameters ###
##########################
#tau = 1.47 defines gamma:
#ws.factory("{gamma[0.68,0.4,0.9],dm[17.7],dG[0.05,-0.3,0.3]}")
ws.factory("{gamma[0.68,0.4,0.9],dm[17.7],dG[-0.05,-0.3,0.3]}")
ws.factory("expr::tau('1/@0',{gamma})")

ws.factory('{phis[0.]}')

ws.factory("{expr::S('-1*sin(phis)',{phis}),expr::D('cos(phis)',{phis}),C[0]}")
ws.factory("{expr::Sold('sin(phis)',{phis}),expr::Dold('cos(phis)',{phis}),Cold[0]}")

###############################
### Experimental parameters ###
###############################

ws.factory("RooGaussModel::res(t,mu[0],sigma[0.05])")

ws.factory("{wmistag[0.5]}")
#ws.factory("{wmistag[0.37]}")
wmistag = ws.var('wmistag')

##################################### building the NEW pdf ###################################
# create observables
ws.factory("expr::ReAz   ('rz    * cos(deltaz)',   {rz,deltaz})")
ws.factory("expr::ImAz   ('rz    * sin(deltaz)',   {rz,deltaz})")
ws.factory("expr::ReApar ('rpar  * cos(deltapar)', {rpar,deltapar})")
ws.factory("expr::ImApar ('rpar  * sin(deltapar)', {rpar,deltapar})")
ws.factory("expr::ReAperp('rperp * cos(deltaperp)',{rperp,deltaperp})")
ws.factory("expr::ImAperp('rperp * sin(deltaperp)',{rperp,deltaperp})")

#This one imports automatically in the workspace!
newpdf = buildJpsiphi(ws,'newpdf') 

############################
### Writing to workspace ###
############################

file = TFile("NewWS.root","RECREATE")
ws.Write("workspace")
file.Close()



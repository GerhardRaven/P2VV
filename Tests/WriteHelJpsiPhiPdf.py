#######################################
### Author: Daan van Eijk
### Updated on: Jun 5 11
### Description: This script generates the simplest possible signal JpsiPhi Pdf and writes it to a workspace in a rootfile
###              In HELICITY COORDINATES!!!
###              Do this to check that the choice of angles in the parameterization of the angles
###              is independent of the transversity amplitudes that we fit later.
########################################
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

#ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f], t[-2,20.], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0]}"%(-pi,pi))
#ws.defineSet("transversityangles","trcospsi,trcostheta,trphi")

ws.factory("{ helcosthetaK[-1,1], helcosthetaL[-1,1], helphi[%f,%f], t[-2,20.], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0]}"%(-pi,pi))
ws.defineSet("helicityangles","helcosthetaK,helcosthetaL,helphi")

ws.factory("{rz2[0.460,0.4,0.7],rperp2[0.347,0.2,0.5]")
ws.factory("RooFormulaVar::rpar2('1-@0-@1',{rz2,rperp2})")
ws.factory("RooFormulaVar::rz('sqrt(@0)',{rz2})")
ws.factory("RooFormulaVar::rperp('sqrt(@0)',{rperp2})")
ws.factory("RooFormulaVar::rpar('sqrt(@0)',{rpar2})")

ws.factory("{deltaz[0.],deltapar[-2.93],deltaperp[2.91]}")

##########################
### physics parameters ###
##########################
#tau = 1.47 defines gamma:
#ws.factory("{gamma[0.68,0.4,0.9],dm[17.7],dG[0.05,-0.3,0.3]}")
ws.factory("{gamma[0.68,0.4,0.9],t_sig_dm[17.7],t_sig_dG[-0.05,-0.3,0.3]}")
ws.factory("expr::t_sig_tau('1/@0',{gamma})")

ws.factory('{phis[0.]}')

ws.factory("{expr::S('-1*sin(phis)',{phis}),expr::D('cos(phis)',{phis}),C[0]}")
ws.factory("{expr::Sold('sin(phis)',{phis}),expr::Dold('cos(phis)',{phis}),Cold[0]}")

###############################
### Experimental parameters ###
###############################

ws.factory("RooGaussModel::tres_sig(t,mu[0],sigma[0.05])")

ws.factory("{wtag[0.5]}")

##################################### building the NEW pdf ###################################
# create observables
ws.factory("expr::ReAz   ('rz    * cos(deltaz)',   {rz,deltaz})")
ws.factory("expr::ImAz   ('rz    * sin(deltaz)',   {rz,deltaz})")
ws.factory("expr::ReApar ('rpar  * cos(deltapar)', {rpar,deltapar})")
ws.factory("expr::ImApar ('rpar  * sin(deltapar)', {rpar,deltapar})")
ws.factory("expr::ReAperp('rperp * cos(deltaperp)',{rperp,deltaperp})")
ws.factory("expr::ImAperp('rperp * sin(deltaperp)',{rperp,deltaperp})")

#This one imports automatically in the workspace!
newpdf = buildJpsiphi(ws,'newpdf', False) 

############################
### Writing to workspace ###
############################

file = TFile("HelWS.root","RECREATE")
ws.Write("workspace")
file.Close()



from ROOT import *

import math
from math import pi

#Old
import GaudiPython
P2VV = GaudiPython.gbl.P2VV
#to load functions (made with namespace function) like makePVVPdf:
GaudiPython.loaddict('P2VVDict')

#New
from itertools import count
gSystem.Load("libp2vv")

from ModelBuilders import buildJpsiphi,buildJpsikstar


################################
### OLDOLDOLDOLDOLDOLDOLDOLD ###
################################



########################
### Define workspace ###
########################

ws = RooWorkspace("ws")

###################
### Observables ###
###################

ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f], t[-2,20.], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0]}"%(-pi,pi))

ws.factory("{rz[0.463],rpar[0.211],rperp[0.347]}")
ws.factory("{deltaz[0],deltapar[-2.93],deltaperp[2.91]}")

t      = ws.var('t')
ctheta = ws.var('trcostheta')
cpsi   = ws.var('trcospsi')
phi    = ws.var('trphi')

myThreeAngles = P2VV.ThreeAngles(ctheta,cpsi,phi)

############ organizing the flavour tagging categories ############ 

### to generate
# bsjpsiphi (when you generate Jpsiphi, note: name = "state")
#The names of this RooCategory are inspired by the common names in the NTuple....
jpsiphiState = ws.cat('tagdecision')

#### and now map the states to the real categories in the fit pdf ####

#Flavor initial state
tagFlavCatJpsiphi = RooMappedCategory("tagFlavCatJpsiphi","tagFlavCatJpsiphi",jpsiphiState)
tagFlavCatJpsiphi.map("Bs_*",    "Bstag",   +1)
tagFlavCatJpsiphi.map("Bsbar_*", "Bsbartag",-1)

#Flavor Final state
recFlavCatJpsiphi = RooMappedCategory("recFlavCatJpsiphi","recFlavCatJpsiphi",jpsiphiState,"*_Jpsiphi",0)

#### finally also define some ranges (for practical purposes like plotting) ####
jpsiphiState.setRange("BsJpsiphi",   "Bs_Jpsiphi")
jpsiphiState.setRange("BsbarJpsiphi","Bsbar_Jpsiphi")

###########################
### Sets of Observables ###
###########################

#to generate data 
dataObsSetJpsiphi = RooArgSet(t,
                              myThreeAngles.ctheta(),
                              myThreeAngles.cpsi(),
                              myThreeAngles.phi(),
                              jpsiphiState)


####################
### Ang eff corr ###
####################

#todo: make these variables in the constructor of the list
phi1_jpsiphi = RooRealVar("phi1_jpsiphi", "phi1_jpsiphi", 32./3.*pi/3.)
phi2_jpsiphi = RooRealVar("phi2_jpsiphi", "phi2_jpsiphi", 32./3.*pi/3.)
phi3_jpsiphi = RooRealVar("phi3_jpsiphi", "phi3_jpsiphi", 32./3.*pi/3.)
phi4_jpsiphi = RooRealVar("phi4_jpsiphi", "phi4_jpsiphi",           0.)
phi5_jpsiphi = RooRealVar("phi5_jpsiphi", "phi5_jpsiphi",           0.)
phi6_jpsiphi = RooRealVar("phi6_jpsiphi", "phi6_jpsiphi",           0.)

myAngularAcceptanceListJpsiphi = P2VV.THREEANGLES.FFSS.AngularAcceptanceCorrectionList(phi1_jpsiphi,
                                                                                       phi2_jpsiphi,
                                                                                       phi3_jpsiphi,
                                                                                       phi4_jpsiphi,
                                                                                       phi5_jpsiphi,
                                                                                       phi6_jpsiphi)

#####################
### amps & phases ###
#####################

#you can (should) fix the sum to one
Azero_sq_jpsiphi = ws.var('rz')
Aperp_sq_jpsiphi = ws.var('rperp')
Apar_sq_jpsiphi  = ws.var('rpar')

phi_perp_jpsiphi = ws.var('deltaperp')
phi_par_jpsiphi  = ws.var('deltapar')
phi_zero_jpsiphi = ws.var('deltaz')

###################################
### Strong phases or parameters ###
###################################

Aperp_Azero_cos_jpsiphi = RooFormulaVar("Aperp_Azero_cos_jpsiphi","sqrt(@0)*sqrt(@1)*cos(@2-@3)",
                                        RooArgList(Azero_sq_jpsiphi, Aperp_sq_jpsiphi,phi_perp_jpsiphi,phi_zero_jpsiphi))
Aperp_Azero_sin_jpsiphi = RooFormulaVar("Aperp_Azero_sin_jpsiphi","sqrt(@0)*sqrt(@1)*sin(@2-@3)",
                                        RooArgList(Azero_sq_jpsiphi, Aperp_sq_jpsiphi,phi_perp_jpsiphi,phi_zero_jpsiphi))
Azero_Apar_cos_jpsiphi  = RooFormulaVar("Azero_Apar_cos_jpsiphi", "sqrt(@0)*sqrt(@1)*cos(@2-@3)",
                                        RooArgList(Azero_sq_jpsiphi, Apar_sq_jpsiphi, phi_par_jpsiphi ,phi_zero_jpsiphi))
Aperp_Apar_cos_jpsiphi  = RooFormulaVar("Aperp_Apar_cos_jpsiphi", "sqrt(@0)*sqrt(@1)*cos(@2-@3)",
                                        RooArgList(Apar_sq_jpsiphi, Aperp_sq_jpsiphi, phi_perp_jpsiphi,phi_par_jpsiphi))
Aperp_Apar_sin_jpsiphi  = RooFormulaVar("Aperp_Apar_sin_jpsiphi", "sqrt(@0)*sqrt(@1)*sin(@2-@3)",
                                        RooArgList(Apar_sq_jpsiphi, Aperp_sq_jpsiphi, phi_perp_jpsiphi,phi_par_jpsiphi)) 
    
###################################
### Amplitudes and phases lists ###
###################################

myAmplitudesAndPhasesListJpsiphi = P2VV.AmplitudesAndPhasesList(Azero_sq_jpsiphi,
                                                                Apar_sq_jpsiphi,
                                                                Aperp_sq_jpsiphi,
                                                                Aperp_Azero_cos_jpsiphi,
                                                                Aperp_Azero_sin_jpsiphi,
                                                                Azero_Apar_cos_jpsiphi,
                                                                Aperp_Apar_cos_jpsiphi,
                                                                Aperp_Apar_sin_jpsiphi)


##########################
### physics parameters ###
##########################
#tau = 1.47 defines gamma:
ws.factory("{gamma[0.68,0.4,0.9],dm[17.7],dG[0.05,-0.3,0.3]}")
ws.factory("expr::tau('1/@0',{gamma})")
tau = ws.function('tau')
gamma = ws.var('gamma')
dG = ws.var('dG')
dm = ws.var('dm')
############################################
### Configuring other physics parameters ###
############################################

dGG  = RooFormulaVar("dGG","#Delta#Gamma_{s}/#Gamma_{s}","@0/@1",RooArgList(dG,gamma))

ws.factory('{phis[0.8]}')
phis = ws.var('phis')

ws.factory("{expr::S('sin(phis)',{phis}),expr::D('cos(phis)',{phis}),C[0]}")
S = ws.function('S')
D = ws.function('D')

################################
### physics parameters lists ###
################################

myPhysicsParametersListJpsiphi = P2VV.PhysParamList(tau,
                                                    dGG,
                                                    dm,
                                                    S,
                                                    D)

###############################
### Experimental parameters ###
###############################

ws.factory("RooGaussModel::res(t,mu[0],sigma[0.05])")
res = ws.pdf('res')

#ws.factory("{wmistag[0.5]}")
ws.factory("{wmistag[0.37]}")
wmistag = ws.var('wmistag')

##################################### building the PDFs ###################################
oldpdf =             P2VV.Functions.makePVVPdf("oldpdf",
                                               "oldpdf",
                                               myThreeAngles,
                                               t,
                                               tagFlavCatJpsiphi,
                                               recFlavCatJpsiphi,
                                               myAmplitudesAndPhasesListJpsiphi,
                                               myPhysicsParametersListJpsiphi,
                                               res,
                                               wmistag,
                                               "Jpsiphi")

getattr(ws,'import')(oldpdf)

################################
### NEWNEWNEWNEWNEWNEWNEWNEW ###
################################
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

file = TFile("WSBoth.root","RECREATE")
ws.Write("workspace")
file.Close()

##################
### Plotjes....###
##################

obsNames =[ 'trcospsi','trcostheta','trphi','t','tagdecision' ]
obs = ws.argSet(','.join(obsNames))

canvasold = TCanvas('canvasold','canvasold')
canvasold.Divide(5,4)

p = [ 'rz','rpar','rperp' ]
for i in  range(len(p)+1) :
    if i>0 : 
        for j in range(len(p)) : 
            x = ws.var(p[j])
            x.setVal( 1 if i==j+1 else 0 )
    for l in p + ['deltaz','deltapar','deltaperp' ] : ws.var(l).Print()

    for (j,k) in zip(['trcospsi','trcostheta','trphi','t'],count(5*i+1)) :
       canvasold.cd(k)
       f = ws.var(j).frame() 
       proj = RooArgSet( obs )
       proj.remove( ws.var(j) )
       oldpdf.plotOn(f,RooFit.Project(proj),RooFit.LineColor(kRed))
       f.Draw()

    tagAsym = RooFit.Asymmetry(ws.cat("tagdecision"))
    canvasold.cd(5*i+5)
    f = ws.var('t').frame(RooFit.Range(-.5,3.5))
    proj = RooArgSet( obs )
    proj.remove(ws.var('t'))
    oldpdf.plotOn(f,RooFit.Project(proj),tagAsym,RooFit.LineColor(kRed))
    f.Draw()

canvasold.Flush()

ws.var('rz').setVal(0.463)
ws.var('rpar').setVal(0.211)
ws.var('rperp').setVal(0.347)

canvasnew = TCanvas('canvasnew','canvasnew')
canvasnew.Divide(5,4)

p = [ 'rz','rpar','rperp' ]
for i in  range(len(p)+1) :
    if i>0 : 
        for j in range(len(p)) : 
            x = ws.var(p[j])
            x.setVal( 1 if i==j+1 else 0 )
    for l in p + ['deltaz','deltapar','deltaperp' ] : ws.var(l).Print()

    for (j,k) in zip(['trcospsi','trcostheta','trphi','t'],count(5*i+1)) :
       canvasnew.cd(k)
       f = ws.var(j).frame() 
       proj = RooArgSet( obs )
       proj.remove( ws.var(j) )
       newpdf.plotOn(f,RooFit.Project(proj))
       f.Draw()

    tagAsym = RooFit.Asymmetry(ws.cat("tagdecision"))
    canvasnew.cd(5*i+5)
    f = ws.var('t').frame(RooFit.Range(-.5,3.5))
    proj = RooArgSet( obs )
    proj.remove(ws.var('t'))
    newpdf.plotOn(f,RooFit.Project(proj),tagAsym,RooFit.LineColor(kBlue))
    f.Draw()

canvasnew.Flush()

